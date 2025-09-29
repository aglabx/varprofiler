#!/usr/bin/env python3
"""
Centromere analysis utilities for VarProfiler.
"""

import numpy as np
import subprocess
from pathlib import Path
from .utils import find_binary_in_path, merge_adjacent_regions


def score_centromere_likelihood(kmer_percent, cenpb_density, array_size, is_minimum=False):
    """
    Calculate centromere likelihood score based on multiple features.
    
    Args:
        kmer_percent: K-mer diversity percentage
        cenpb_density: CENP-B boxes per kb
        array_size: Size of satellite array
        is_minimum: Whether this is a local minimum
    
    Returns:
        Tuple of (score 0-100, weights dictionary)
    """
    weights = {}
    
    # 1. K-mer diversity score (30% weight) - lower is better
    if kmer_percent < 2:
        weights['kmer'] = 30
    elif kmer_percent < 5:
        weights['kmer'] = 25
    elif kmer_percent < 10:
        weights['kmer'] = 15
    elif kmer_percent < 15:
        weights['kmer'] = 5
    else:
        weights['kmer'] = 0
    
    # 2. CENP-B box density score (40% weight) - most important
    if cenpb_density > 2.5:  # Alpha-satellite rich (>2.5/kb)
        weights['cenpb'] = 40
    elif cenpb_density > 2.0:  # High density
        weights['cenpb'] = 35
    elif cenpb_density > 1.5:  # Good density
        weights['cenpb'] = 30
    elif cenpb_density > 1.0:  # Moderate (typical mixed arrays)
        weights['cenpb'] = 25
    elif cenpb_density > 0.5:  # Low but present
        weights['cenpb'] = 15
    elif cenpb_density > 0:  # Very low (degenerate)
        weights['cenpb'] = 5
    else:
        weights['cenpb'] = 0
    
    # 3. Array size context (20% weight)
    if array_size > 5_000_000:  # >5 Mb
        weights['array_size'] = 20
    elif array_size > 2_000_000:  # >2 Mb
        weights['array_size'] = 15
    elif array_size > 1_000_000:  # >1 Mb
        weights['array_size'] = 10
    elif array_size > 500_000:  # >500 kb
        weights['array_size'] = 5
    else:
        weights['array_size'] = 2
    
    # 4. Local minimum bonus (10% weight)
    if is_minimum:
        weights['local_minimum'] = 10
    else:
        weights['local_minimum'] = 0
    
    # Calculate total score (max 100)
    total_score = sum(weights.values())
    
    return total_score, weights


def find_local_minima(windows_df, window_size=5, column='kmer_percentage'):
    """
    Find local minima in a DataFrame column.
    
    Args:
        windows_df: DataFrame with windows
        window_size: Number of adjacent windows to consider
        column: Column name to find minima in
    
    Returns:
        List of indices that are local minima
    """
    local_minima = []
    values = windows_df[column].values
    
    for i in range(len(values)):
        # Define neighborhood
        start_idx = max(0, i - window_size // 2)
        end_idx = min(len(values), i + window_size // 2 + 1)
        
        neighborhood = values[start_idx:end_idx]
        
        # Check if current is minimum in neighborhood
        if len(neighborhood) > 0 and values[i] == min(neighborhood):
            # Also check if it's significantly lower than neighborhood mean
            neighborhood_mean = np.mean(neighborhood)
            if values[i] < neighborhood_mean * 0.7:  # At least 30% lower
                local_minima.append(windows_df.index[i])
    
    return local_minima


def run_cenpb_finder(bed_file, genome_file, output_file, 
                    pattern="YTTCGTTGGAARCGGGA",
                    max_edit_distance=3, threads=1):
    """
    Run cenpb_finder tool to search for CENP-B boxes in regions.
    
    Args:
        bed_file: BED file with regions to search
        genome_file: Genome FASTA file
        output_file: Output TSV file
        pattern: CENP-B box pattern (default: canonical human pattern)
        max_edit_distance: Maximum allowed edit distance
        threads: Number of threads to use
    
    Returns:
        Dictionary with region coordinates as keys and CENP-B data as values
    """
    # Find cenpb_finder binary
    cenpb_finder = find_binary_in_path('cenpb_finder')
    
    if not cenpb_finder:
        print("Warning: cenpb_finder not found. Searched in:")
        print("  - Package directory")
        print("  - Current directory")
        print("  - System PATH")
        print("Please compile it with 'make cenpb_finder' in the VarProfiler directory.")
        print("Skipping CENP-B box analysis.")
        return {}
    
    # Build command
    cmd = [
        str(cenpb_finder),
        genome_file,
        bed_file,
        output_file,
        '-p', pattern,
        '-d', str(max_edit_distance),
        '-t', str(threads)
    ]
    
    try:
        print(f"Searching for CENP-B boxes using pattern {pattern}...")
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Warning: cenpb_finder failed: {e}")
        return {}
    
    # Parse results
    cenpb_data = {}
    try:
        with open(output_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    cenpb_count = int(parts[4])
                    cenpb_density = float(parts[5])
                    key = f"{chrom}:{start}-{end}"
                    cenpb_data[key] = {
                        'count': cenpb_count,
                        'density': cenpb_density,
                        'positions': parts[6] if len(parts) > 6 and parts[6] != '.' else None
                    }
    except Exception as e:
        print(f"Warning: Could not parse CENP-B results: {e}")
    
    return cenpb_data


def analyze_satellite_arrays(df, regions, kmer_size):
    """
    Analyze satellite arrays to find the best centromere candidates.
    
    Args:
        df: DataFrame with k-mer counts for all windows
        regions: Dictionary of low k-mer diversity regions
        kmer_size: K-mer size used in analysis
    
    Returns:
        Dictionary with satellite array analysis
    """
    satellite_arrays = {}
    
    for chrom, chrom_regions in regions.items():
        # Merge nearby regions into arrays
        merged_arrays = merge_adjacent_regions(chrom_regions, max_gap=500000)
        
        for array_idx, (array_start, array_end, array_mean) in enumerate(merged_arrays):
            array_size = array_end - array_start
            
            # Get all windows in this array
            chrom_df = df[df['chromosome'] == chrom].copy()
            if 'kmer_percentage' not in chrom_df.columns:
                chrom_df['window_size'] = chrom_df['end'] - chrom_df['start']
                chrom_df['max_kmers'] = chrom_df['window_size'] - kmer_size + 1
                chrom_df['kmer_percentage'] = (chrom_df['kmer_count'] / chrom_df['max_kmers']) * 100
            
            array_windows = chrom_df[(chrom_df['start'] >= array_start) & 
                                    (chrom_df['end'] <= array_end)]
            
            if len(array_windows) > 0:
                # Find local minima (potential centromeres)
                local_minima = find_local_minima(array_windows)
                
                array_key = f"{chrom}:{array_start}-{array_end}"
                satellite_arrays[array_key] = {
                    'chrom': chrom,
                    'start': array_start,
                    'end': array_end,
                    'size': array_size,
                    'mean_kmer_percent': array_mean,
                    'windows': array_windows,
                    'local_minima': local_minima,
                    'n_windows': len(array_windows)
                }
    
    return satellite_arrays


def select_best_centromeres(candidates, one_per_chrom=True, min_score=50):
    """
    Select the best centromere candidate(s) from all candidates.
    
    Args:
        candidates: List of centromere candidates with scores
        one_per_chrom: If True, select only one per chromosome
        min_score: Minimum score threshold
    
    Returns:
        List of selected centromeres
    """
    selected = []
    seen_chroms = set()
    
    # Sort by score
    sorted_candidates = sorted(candidates, key=lambda x: x.get('score', 0), reverse=True)
    
    for candidate in sorted_candidates:
        chrom = candidate['chrom']
        score = candidate.get('score', 0)
        
        # Check minimum score
        if score < min_score:
            continue
        
        if one_per_chrom:
            if chrom not in seen_chroms:
                selected.append(candidate)
                seen_chroms.add(chrom)
        else:
            selected.append(candidate)
    
    return selected


def classify_confidence_level(score):
    """
    Classify centromere confidence based on score.
    
    Args:
        score: Centromere likelihood score (0-100)
    
    Returns:
        Tuple of (confidence_level, color_code)
    """
    if score >= 70:
        return 'high', '255,0,0'  # Red
    elif score >= 50:
        return 'medium', '255,128,0'  # Orange
    else:
        return 'low', '255,255,0'  # Yellow


def get_centromere_statistics(centromeres):
    """
    Calculate statistics for predicted centromeres.
    
    Args:
        centromeres: List of centromere predictions
    
    Returns:
        Dictionary with statistics
    """
    if not centromeres:
        return {}
    
    scores = [c.get('score', 0) for c in centromeres]
    
    stats = {
        'total': len(centromeres),
        'avg_score': np.mean(scores),
        'median_score': np.median(scores),
        'high_confidence': sum(1 for s in scores if s >= 70),
        'medium_confidence': sum(1 for s in scores if 50 <= s < 70),
        'low_confidence': sum(1 for s in scores if s < 50),
        'chromosomes_with_centromeres': len(set(c['chrom'] for c in centromeres))
    }
    
    # CENP-B statistics if available
    cenpb_counts = [c.get('cenpb_count', 0) for c in centromeres]
    if any(count > 0 for count in cenpb_counts):
        cenpb_densities = [c.get('cenpb_density', 0) for c in centromeres if c.get('cenpb_count', 0) > 0]
        stats['with_cenpb'] = sum(1 for count in cenpb_counts if count > 0)
        stats['avg_cenpb_density'] = np.mean(cenpb_densities) if cenpb_densities else 0
        stats['max_cenpb_density'] = max(cenpb_densities) if cenpb_densities else 0
    
    return stats