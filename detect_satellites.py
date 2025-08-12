#!/usr/bin/env python3
"""
Detect regions with low k-mer variability (potential satellite DNA)
and export them as GFF annotations and FASTA sequences.
"""

import pandas as pd
import numpy as np
import argparse
import sys
import os
from collections import defaultdict


def calculate_threshold(df, percentile=5):
    """
    Calculate threshold for low variability regions.
    
    Args:
        df: DataFrame with k-mer counts
        percentile: Percentile to use as threshold (default: 5th percentile)
    
    Returns:
        dict: Thresholds per chromosome and global threshold
    """
    # Calculate k-mer percentage if not present
    if 'kmer_percentage' not in df.columns:
        df['window_size'] = df['end'] - df['start']
        df['max_kmers'] = df['window_size'] - df['kmer_size'] + 1
        df['kmer_percentage'] = (df['kmer_count'] / df['max_kmers']) * 100
    
    # Calculate global threshold
    global_threshold = np.percentile(df['kmer_percentage'], percentile)
    
    # Calculate per-chromosome thresholds
    chrom_thresholds = {}
    for chrom in df['chromosome'].unique():
        chrom_df = df[df['chromosome'] == chrom]
        chrom_thresholds[chrom] = np.percentile(chrom_df['kmer_percentage'], percentile)
    
    return global_threshold, chrom_thresholds


def merge_adjacent_windows(regions, max_gap=50000):
    """
    Merge adjacent low-variability windows into continuous regions.
    
    Args:
        regions: List of (start, end, value) tuples
        max_gap: Maximum gap between windows to merge (default: 50kb)
    
    Returns:
        List of merged (start, end, mean_value) tuples
    """
    if not regions:
        return []
    
    # Sort by start position
    regions.sort(key=lambda x: x[0])
    
    merged = []
    current_start = regions[0][0]
    current_end = regions[0][1]
    current_values = [regions[0][2]]
    
    for start, end, value in regions[1:]:
        if start - current_end <= max_gap:
            # Merge with current region
            current_end = end
            current_values.append(value)
        else:
            # Save current region and start new one
            merged.append((current_start, current_end, np.mean(current_values)))
            current_start = start
            current_end = end
            current_values = [value]
    
    # Save last region
    merged.append((current_start, current_end, np.mean(current_values)))
    
    return merged


def detect_low_variability_regions(bed_file, kmer_size, output_prefix, 
                                  percentile=5, min_length=10000,
                                  use_global_threshold=False):
    """
    Detect regions with low k-mer variability from BED file.
    
    Args:
        bed_file: Input BED file from kmer_profiler
        kmer_size: K-mer size used in analysis
        output_prefix: Prefix for output files
        percentile: Percentile threshold for low variability
        min_length: Minimum region length to report (default: 10kb)
        use_global_threshold: Use global instead of per-chromosome threshold
    """
    print(f"Reading k-mer data from {bed_file}...")
    
    # Read BED file
    col_names = ['chromosome', 'start', 'end', 'kmer_count']
    df = pd.read_csv(bed_file, sep='\t', header=None, names=col_names)
    df['kmer_size'] = kmer_size
    
    # Calculate thresholds
    global_threshold, chrom_thresholds = calculate_threshold(df, percentile)
    
    print(f"Global threshold ({percentile}th percentile): {global_threshold:.2f}%")
    
    # Detect low variability regions for each chromosome
    all_regions = defaultdict(list)
    
    for chrom in df['chromosome'].unique():
        chrom_df = df[df['chromosome'] == chrom].copy()
        
        # Calculate window size and percentage
        chrom_df['window_size'] = chrom_df['end'] - chrom_df['start']
        chrom_df['max_kmers'] = chrom_df['window_size'] - kmer_size + 1
        chrom_df['kmer_percentage'] = (chrom_df['kmer_count'] / chrom_df['max_kmers']) * 100
        
        # Choose threshold
        if use_global_threshold:
            threshold = global_threshold
        else:
            threshold = chrom_thresholds[chrom]
        
        print(f"  {chrom}: threshold = {threshold:.2f}%")
        
        # Find windows below threshold
        low_var = chrom_df[chrom_df['kmer_percentage'] < threshold]
        
        if len(low_var) > 0:
            # Extract regions
            regions = [(row['start'], row['end'], row['kmer_percentage']) 
                      for _, row in low_var.iterrows()]
            
            # Merge adjacent windows
            merged = merge_adjacent_windows(regions)
            
            # Filter by minimum length
            filtered = [(start, end, value) for start, end, value in merged 
                       if end - start >= min_length]
            
            all_regions[chrom] = filtered
            print(f"    Found {len(filtered)} regions >= {min_length}bp")
    
    return all_regions, df


def organize_regions_by_size(regions):
    """
    Organize low variability regions by their size for reporting.
    
    Note: We do NOT classify these as specific repeat types, just organize
    by size for convenient analysis. The actual repeat type depends on many
    factors including species, chromosome, and repeat unit period.
    
    Args:
        regions: Dictionary of regions per chromosome
    
    Returns:
        List of all regions with size information
    """
    all_regions = []
    
    for chrom, chrom_regions in regions.items():
        for start, end, value in chrom_regions:
            length = end - start
            region_info = {
                'chrom': chrom,
                'start': start,
                'end': end,
                'length': length,
                'mean_kmer_percentage': value,
                'size_category': get_size_category(length)
            }
            all_regions.append(region_info)
    
    # Sort by chromosome and position
    all_regions.sort(key=lambda x: (x['chrom'], x['start']))
    
    return all_regions


def get_size_category(length):
    """
    Get size category for a region (for descriptive purposes only).
    
    Args:
        length: Region length in bp
    
    Returns:
        String describing the size range
    """
    if length >= 1000000:
        return '>1Mb'
    elif length >= 100000:
        return '100kb-1Mb'
    elif length >= 10000:
        return '10-100kb'
    elif length >= 1000:
        return '1-10kb'
    else:
        return '<1kb'


def run_cenpb_finder(bed_file, genome_file, output_file, pattern="YTTCGTTGGAARCGGGA", 
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
        Dictionary with region coordinates as keys and CENP-B counts as values
    """
    import subprocess
    from pathlib import Path
    import os
    
    # Look for cenpb_finder in multiple locations
    script_dir = Path(os.path.dirname(os.path.abspath(__file__)))
    possible_paths = [
        script_dir / 'cenpb_finder',  # Same directory as this script
        Path('./cenpb_finder'),        # Current working directory
        Path('cenpb_finder'),          # In PATH
    ]
    
    cenpb_finder = None
    for path in possible_paths:
        if path.exists() and os.access(path, os.X_OK):
            cenpb_finder = path
            break
    
    # Also check if it's in PATH
    if not cenpb_finder:
        from shutil import which
        in_path = which('cenpb_finder')
        if in_path:
            cenpb_finder = Path(in_path)
    
    if not cenpb_finder:
        print(f"Warning: cenpb_finder not found. Searched in:")
        print(f"  - {script_dir}")
        print(f"  - Current directory")
        print(f"  - System PATH")
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


def score_centromere_likelihood(window, kmer_percent, cenpb_density, array_size, is_minimum=False):
    """
    Calculate centromere likelihood score based on multiple features.
    
    Returns score 0-100 and detailed weights.
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
    # Expected: ~3 boxes/kb for alpha-satellite rich regions
    # Some arrays (Y chromosome, acrocentrics) have 0-1/kb
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


def find_centromere_candidates(regions, df, kmer_size, cenpb_data=None):
    """
    Find potential centromere-containing regions using multi-signal approach:
    1. CENP-B box density (primary signal)
    2. K-mer diversity minima
    3. Large satellite array context
    4. Local minimum within array
    
    Args:
        regions: Dictionary of low variability regions per chromosome
        df: Original dataframe with all windows
        kmer_size: K-mer size
        cenpb_data: Dictionary with CENP-B box counts
    
    Returns:
        List of windows likely containing functional centromeres, scored by confidence
    """
    centromeres = []
    
    for chrom, chrom_regions in regions.items():
        # Merge nearby regions into satellite arrays
        merged_arrays = merge_adjacent_windows(chrom_regions, max_gap=500000)
        
        chrom_df = df[df['chromosome'] == chrom].copy()
        if 'kmer_percentage' not in chrom_df.columns:
            chrom_df['window_size'] = chrom_df['end'] - chrom_df['start']
            chrom_df['max_kmers'] = chrom_df['window_size'] - kmer_size + 1
            chrom_df['kmer_percentage'] = (chrom_df['kmer_count'] / chrom_df['max_kmers']) * 100
        
        # Analyze each satellite array
        for array_start, array_end, array_mean in merged_arrays:
            array_size = array_end - array_start
            
            # Skip small arrays unlikely to contain centromeres
            if array_size < 100000:
                continue
            
            # Get all windows in this array
            array_windows = chrom_df[(chrom_df['start'] >= array_start) & 
                                    (chrom_df['end'] <= array_end)]
            
            if len(array_windows) == 0:
                continue
            
            # Find local minima
            min_idx = array_windows['kmer_percentage'].idxmin()
            min_kmer_pct = array_windows.loc[min_idx, 'kmer_percentage']
            
            # Score each window in the array
            window_scores = []
            for idx, window in array_windows.iterrows():
                # Get CENP-B data
                window_key = f"{chrom}:{int(window['start'])}-{int(window['end'])}"
                cenpb_density = 0
                cenpb_count = 0
                if cenpb_data and window_key in cenpb_data:
                    cenpb_density = cenpb_data[window_key]['density']
                    cenpb_count = cenpb_data[window_key]['count']
                
                # Check if this is a local minimum
                is_minimum = (idx == min_idx) or (window['kmer_percentage'] < array_mean * 0.7)
                
                # Calculate score
                score, weights = score_centromere_likelihood(
                    window,
                    window['kmer_percentage'],
                    cenpb_density,
                    array_size,
                    is_minimum
                )
                
                window_scores.append({
                    'window': window,
                    'score': score,
                    'weights': weights,
                    'cenpb_density': cenpb_density,
                    'cenpb_count': cenpb_count,
                    'is_minimum': is_minimum
                })
            
            # Select best scoring window(s) from this array
            # If CENP-B data available, prioritize by score
            # Otherwise fall back to minimum k-mer diversity
            window_scores.sort(key=lambda x: x['score'], reverse=True)
            
            # Take the best scoring window (or multiple if score is high)
            for ws in window_scores:
                if ws['score'] >= 30:  # Minimum threshold
                    window = ws['window']
                    cent_info = {
                        'chrom': chrom,
                        'start': int(window['start']),
                        'end': int(window['end']),
                        'kmer_percentage': window['kmer_percentage'],
                        'satellite_start': array_start,
                        'satellite_end': array_end,
                        'satellite_mean': array_mean,
                        'array_size': array_size,
                        'score': ws['score'],
                        'weights': ws['weights'],
                        'cenpb_count': ws['cenpb_count'],
                        'cenpb_density': ws['cenpb_density'],
                        'is_local_minimum': ws['is_minimum']
                    }
                    centromeres.append(cent_info)
                    break  # Take only best per array (can be changed)
    
    # Sort by score and select best per chromosome
    centromeres.sort(key=lambda x: x['score'], reverse=True)
    
    # Keep only best centromere per chromosome
    best_per_chrom = {}
    for cent in centromeres:
        chrom = cent['chrom']
        if chrom not in best_per_chrom or cent['score'] > best_per_chrom[chrom]['score']:
            best_per_chrom[chrom] = cent
    
    return list(best_per_chrom.values())


def write_gff(regions, centromeres, output_file, source='VarProfiler'):
    """
    Write low k-mer variability regions to GFF3 format file.
    
    Args:
        regions: List of low variability regions
        centromeres: List of potential centromere regions
        output_file: Output GFF file path
        source: Source name for GFF
    """
    print(f"\nWriting GFF annotations to {output_file}...")
    
    with open(output_file, 'w') as f:
        # Write header
        f.write("##gff-version 3\n")
        f.write(f"##source {source}\n")
        f.write("##feature-ontology SO\n")
        f.write("# Regions with low k-mer diversity detected by VarProfiler\n")
        f.write("# These regions have reduced k-mer vocabulary (fewer unique k-mers than expected)\n")
        f.write("# May correspond to tandem repeats, satellites, or other repetitive elements\n")
        
        # ID counter
        feature_id = 1
        
        # Write low variability regions
        for region in regions:
            # Feature type - region with reduced k-mer diversity
            feature_type = 'low_kmer_diversity_region'
            
            # Color based on size (gradient from gray to red)
            if region['length'] >= 1000000:
                color = '255,0,0'  # Red for very large
            elif region['length'] >= 100000:
                color = '255,128,0'  # Orange for large
            elif region['length'] >= 10000:
                color = '255,255,0'  # Yellow for medium
            else:
                color = '128,128,128'  # Gray for small
                
            # GFF fields (1-based coordinates)
            chrom = region['chrom']
            start = region['start'] + 1  # Convert to 1-based
            end = region['end']
            score = f"{region['mean_kmer_percentage']:.2f}"
            strand = '.'
            phase = '.'
            
            # Attributes
            attributes = [
                f"ID=low_kmer_div_{feature_id}",
                f"Name=low_kmer_div_{chrom}_{start}_{end}",
                f"length={region['length']}",
                f"size_range={region['size_category']}",
                f"mean_kmer_percent={region['mean_kmer_percentage']:.2f}"
            ]
            
            # Add TRF information if available
            if 'dominant_period' in region and region['dominant_period']:
                attributes.extend([
                    f"trf_period={region['dominant_period']}",
                    f"trf_gc={region['dominant_gc']:.2f}",
                    f"trf_pmatch={region['dominant_pmatch']:.1f}",
                    f"trf_coverage={region['trf_coverage']:.2f}",
                    f"repeat_class={region['repeat_class']}"
                ])
            
            attributes.append(f"color={color}")
            
            # Write GFF line
            f.write(f"{chrom}\t{source}\t{feature_type}\t{start}\t{end}\t"
                   f"{score}\t{strand}\t{phase}\t{';'.join(attributes)}\n")
            
            feature_id += 1
        
        # Write centromere candidates
        for cent in centromeres:
            chrom = cent['chrom']
            start = cent['start'] + 1  # Convert to 1-based
            end = cent['end']
            score = f"{cent['kmer_percentage']:.2f}"
            
            # Determine confidence level based on score
            score = cent.get('score', 0)
            if score >= 70:
                confidence = 'high'
                color = '255,0,0'  # Red for high confidence
            elif score >= 50:
                confidence = 'medium'
                color = '255,128,0'  # Orange for medium
            else:
                confidence = 'low'
                color = '255,255,0'  # Yellow for low
            
            attributes = [
                f"ID=centromere_containing_region_{feature_id}",
                f"Name=centromere_region_{chrom}_{start}_{end}",
                f"score={score}",
                f"confidence={confidence}",
                f"kmer_percent={cent['kmer_percentage']:.2f}",
                f"satellite_region={cent['satellite_start']+1}..{cent['satellite_end']}",
                f"array_size={cent.get('array_size', cent['satellite_end']-cent['satellite_start'])}",
                f"satellite_mean={cent['satellite_mean']:.2f}"
            ]
            
            # Add CENP-B box information if available
            if 'cenpb_count' in cent:
                attributes.extend([
                    f"cenpb_count={cent['cenpb_count']}",
                    f"cenpb_density={cent['cenpb_density']:.3f}"
                ])
            
            # Add scoring weights if available
            if 'weights' in cent:
                weight_str = ','.join([f"{k}={v}" for k, v in cent['weights'].items()])
                attributes.append(f"score_weights={weight_str}")
            
            # Add local minimum flag
            if 'is_local_minimum' in cent:
                attributes.append(f"is_local_minimum={cent['is_local_minimum']}")
            
            attributes.append(f"color={color}")
            
            f.write(f"{chrom}\t{source}\tcentromere_containing_region\t{start}\t{end}\t"
                   f"{score}\t.\t.\t{';'.join(attributes)}\n")
            
            feature_id += 1
    
    # Print summary
    total_features = len(regions) + len(centromeres)
    print(f"  Wrote {total_features} features")
    
    # Group by size for summary
    size_counts = {}
    for region in regions:
        size_cat = region['size_category']
        size_counts[size_cat] = size_counts.get(size_cat, 0) + 1
    
    for size_cat in ['>1Mb', '100kb-1Mb', '10-100kb', '1-10kb', '<1kb']:
        if size_cat in size_counts:
            print(f"    {size_cat}: {size_counts[size_cat]} regions")
    
    if centromeres:
        print(f"    centromere_candidates: {len(centromeres)}")


def read_trf_annotations(gff_file):
    """
    Read TRF (Tandem Repeat Finder) annotations from GFF file.
    
    Args:
        gff_file: Path to GFF file with TRF annotations
    
    Returns:
        Dictionary with chromosome as key and list of (start, end, period, gc, pmatch) tuples
    """
    trf_regions = defaultdict(list)
    
    if not gff_file or not os.path.exists(gff_file):
        return trf_regions
    
    try:
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                chrom = parts[0]
                start = int(parts[3]) - 1  # Convert to 0-based
                end = int(parts[4])
                attributes = parts[8]
                
                # Parse attributes
                attr_dict = {}
                for attr in attributes.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attr_dict[key] = value
                
                # Extract period, gc, pmatch
                period = int(attr_dict.get('period', '0'))
                gc = float(attr_dict.get('gc', '0'))
                pmatch = float(attr_dict.get('pmatch', '0'))
                
                trf_regions[chrom].append({
                    'start': start,
                    'end': end,
                    'period': period,
                    'gc': gc,
                    'pmatch': pmatch,
                    'length': end - start
                })
        
        # Sort by start position
        for chrom in trf_regions:
            trf_regions[chrom].sort(key=lambda x: x['start'])
            
        print(f"  Loaded TRF annotations for {len(trf_regions)} chromosomes")
        total_regions = sum(len(v) for v in trf_regions.values())
        print(f"  Total TRF regions: {total_regions}")
        
    except Exception as e:
        print(f"  Warning: Could not read TRF annotations: {e}")
    
    return trf_regions


def find_trf_overlaps(low_div_regions, trf_regions):
    """
    Find overlaps between low k-mer diversity regions and TRF annotations.
    
    Args:
        low_div_regions: List of low k-mer diversity regions
        trf_regions: Dictionary of TRF annotations by chromosome
    
    Returns:
        Updated list of regions with TRF overlap information
    """
    for region in low_div_regions:
        chrom = region['chrom']
        region_start = region['start']
        region_end = region['end']
        
        # Find overlapping TRF regions
        overlaps = []
        if chrom in trf_regions:
            for trf in trf_regions[chrom]:
                # Check for overlap
                if trf['end'] > region_start and trf['start'] < region_end:
                    overlap_start = max(region_start, trf['start'])
                    overlap_end = min(region_end, trf['end'])
                    overlap_length = overlap_end - overlap_start
                    
                    overlaps.append({
                        'period': trf['period'],
                        'gc': trf['gc'],
                        'pmatch': trf['pmatch'],
                        'overlap_length': overlap_length,
                        'overlap_fraction': overlap_length / region['length']
                    })
        
        # Add TRF information to region
        if overlaps:
            # Sort by overlap fraction
            overlaps.sort(key=lambda x: x['overlap_fraction'], reverse=True)
            
            # Get dominant period (from largest overlap)
            region['trf_overlaps'] = len(overlaps)
            region['dominant_period'] = overlaps[0]['period']
            region['dominant_gc'] = overlaps[0]['gc']
            region['dominant_pmatch'] = overlaps[0]['pmatch']
            region['trf_coverage'] = sum(o['overlap_length'] for o in overlaps) / region['length']
            
            # Classify based on period
            period = overlaps[0]['period']
            if period <= 6:
                region['repeat_class'] = 'microsatellite'
            elif period <= 100:
                region['repeat_class'] = 'minisatellite'
            else:
                region['repeat_class'] = 'satellite'
        else:
            region['trf_overlaps'] = 0
            region['dominant_period'] = None
            region['repeat_class'] = 'unknown'
            region['trf_coverage'] = 0.0
    
    return low_div_regions


def read_fasta(filename):
    """
    Read FASTA file and return dictionary of sequences.
    
    Args:
        filename: Path to FASTA file
    
    Returns:
        Dictionary with sequence IDs as keys and sequences as values
    """
    sequences = {}
    current_id = None
    current_seq = []
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous sequence if exists
                    if current_id is not None:
                        sequences[current_id] = ''.join(current_seq)
                    # Start new sequence
                    current_id = line[1:].split()[0]  # Take first word as ID
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Save last sequence
            if current_id is not None:
                sequences[current_id] = ''.join(current_seq)
    
    except FileNotFoundError:
        return None
    
    return sequences


def extract_centromere_sequences(centromeres, genome_file, output_file):
    """
    Extract FASTA sequences specifically for predicted centromeres with detailed headers.
    
    Args:
        centromeres: List of predicted centromere regions with scores
        genome_file: Input genome FASTA file
        output_file: Output FASTA file for centromere sequences
    """
    if not centromeres:
        return
    
    print(f"\nExtracting centromere sequences from {genome_file}...")
    
    # Read genome sequences
    genome = read_fasta(genome_file)
    
    if genome is None:
        print(f"  Error: Could not read genome file {genome_file}")
        return
    
    sequences = []
    
    for i, cent in enumerate(centromeres, 1):
        chrom = cent['chrom']
        if chrom not in genome:
            print(f"  Warning: Chromosome {chrom} not found in genome")
            continue
        
        start = cent['start']
        end = cent['end']
        
        # Extract sequence
        if end <= len(genome[chrom]):
            seq = genome[chrom][start:end]
        else:
            seq = genome[chrom][start:]
        
        # Build detailed header with all information
        header_parts = [
            f"centromere_{i}",
            f"chr={chrom}",
            f"pos={start+1}-{end}",  # 1-based for display
            f"length={end-start}bp",
            f"score={cent.get('score', 0):.1f}"
        ]
        
        # Add confidence level
        score = cent.get('score', 0)
        if score >= 70:
            header_parts.append("confidence=HIGH")
        elif score >= 50:
            header_parts.append("confidence=MEDIUM")
        else:
            header_parts.append("confidence=LOW")
        
        # Add k-mer percentage
        header_parts.append(f"kmer_diversity={cent['kmer_percentage']:.2f}%")
        
        # Add CENP-B information if available
        if 'cenpb_count' in cent:
            header_parts.append(f"cenpb_count={cent['cenpb_count']}")
            header_parts.append(f"cenpb_density={cent['cenpb_density']:.3f}/kb")
        
        # Add array context
        if 'array_size' in cent:
            header_parts.append(f"array_size={cent['array_size']}bp")
        
        # Add scoring weights if available
        if 'weights' in cent:
            weight_str = ','.join([f"{k}:{v}" for k, v in cent['weights'].items()])
            header_parts.append(f"weights=({weight_str})")
        
        # Is it a local minimum?
        if cent.get('is_local_minimum'):
            header_parts.append("local_minimum=yes")
        
        header = " ".join(header_parts)
        sequences.append((header, seq))
    
    # Write FASTA file
    with open(output_file, 'w') as f:
        for header, seq in sequences:
            f.write(f">{header}\n")
            # Write sequence in 80-character lines
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")
    
    print(f"  Extracted {len(sequences)} centromere sequences to {output_file}")
    
    # Print summary
    if sequences:
        total_bp = sum(len(seq) for _, seq in sequences)
        avg_length = total_bp // len(sequences)
        print(f"  Total: {total_bp:,} bp")
        print(f"  Average length: {avg_length:,} bp")
        
        # Report by confidence
        high = sum(1 for h, _ in sequences if "confidence=HIGH" in h)
        med = sum(1 for h, _ in sequences if "confidence=MEDIUM" in h)
        low = sum(1 for h, _ in sequences if "confidence=LOW" in h)
        
        if high > 0:
            print(f"  High confidence: {high}")
        if med > 0:
            print(f"  Medium confidence: {med}")
        if low > 0:
            print(f"  Low confidence: {low}")


def extract_fasta_sequences(regions, centromeres, genome_file, output_file):
    """
    Extract FASTA sequences for detected regions.
    
    Args:
        regions: Dictionary of classified satellite regions
        centromeres: List of potential centromere regions
        genome_file: Input genome FASTA file
        output_file: Output FASTA file
    """
    print(f"\nExtracting sequences from {genome_file}...")
    
    # Read genome sequences
    genome = read_fasta(genome_file)
    
    if genome is None:
        print(f"  Error: Could not read genome file {genome_file}")
        return
    
    print(f"  Loaded {len(genome)} sequences")
    
    # Extract sequences
    sequences = []
    
    # Extract low variability region sequences
    for region in regions:
        chrom = region['chrom']
        if chrom not in genome:
            continue
        
        start = region['start']
        end = region['end']
        
        # Extract sequence substring
        if end <= len(genome[chrom]):
            seq = genome[chrom][start:end]
        else:
            # Handle case where region extends beyond sequence
            seq = genome[chrom][start:]
        
        # Create sequence ID with metadata
        seq_id = f"low_kmer_div_{chrom}:{start+1}-{end}|length={region['length']}|kmer_pct={region['mean_kmer_percentage']:.2f}|size={region['size_category']}"
        sequences.append((seq_id, seq))
    
    # Extract centromere sequences
    for cent in centromeres:
        chrom = cent['chrom']
        if chrom not in genome:
            continue
        
        start = cent['start']
        end = cent['end']
        
        # Extract sequence substring
        if end <= len(genome[chrom]):
            seq = genome[chrom][start:end]
        else:
            seq = genome[chrom][start:]
        
        seq_id = f"centromere_candidate_{chrom}:{start+1}-{end}|kmer_pct={cent['kmer_percentage']:.2f}"
        sequences.append((seq_id, seq))
    
    # Write FASTA file
    with open(output_file, 'w') as f:
        for seq_id, seq in sequences:
            f.write(f">{seq_id}\n")
            # Write sequence in 80-character lines
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")
    
    print(f"  Extracted {len(sequences)} sequences to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Detect low k-mer variability regions (satellite DNA) and export as GFF + FASTA.',
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument('bed_file', 
                       help='Input BED file from kmer_profiler')
    parser.add_argument('-k', '--kmer-size', 
                       type=int, required=True,
                       help='K-mer size used in the analysis')
    parser.add_argument('-o', '--output-prefix',
                       default='satellites',
                       help='Prefix for output files (default: satellites)')
    parser.add_argument('-p', '--percentile',
                       type=float, default=5,
                       help='Percentile threshold for low variability (default: 5)')
    parser.add_argument('-m', '--min-length',
                       type=int, default=10000,
                       help='Minimum region length in bp (default: 10000)')
    parser.add_argument('-g', '--genome',
                       help='Genome FASTA file for sequence extraction')
    parser.add_argument('--global-threshold',
                       action='store_true',
                       help='Use global threshold instead of per-chromosome')
    parser.add_argument('--find-centromeres',
                       action='store_true',
                       help='Attempt to identify centromere candidates')
    parser.add_argument('-t', '--trf-annotations',
                       help='GFF file with TRF annotations for repeat classification')
    parser.add_argument('--cenpb-pattern',
                       default='YTTCGTTGGAARCGGGA',
                       help='CENP-B box pattern (default: YTTCGTTGGAARCGGGA)')
    parser.add_argument('--cenpb-distance',
                       type=int, default=3,
                       help='Maximum edit distance for CENP-B box search (default: 3)')
    parser.add_argument('--cenpb-threads',
                       type=int, default=1,
                       help='Number of threads for CENP-B search (default: 1)')
    
    args = parser.parse_args()
    
    # Detect low variability regions
    regions, df = detect_low_variability_regions(
        args.bed_file, 
        args.kmer_size,
        args.output_prefix,
        args.percentile,
        args.min_length,
        args.global_threshold
    )
    
    if not regions:
        print("\nNo low variability regions found.")
        return
    
    # Organize regions with size information
    organized_regions = organize_regions_by_size(regions)
    
    # Load TRF annotations if provided
    if args.trf_annotations:
        print(f"\nLoading TRF annotations from {args.trf_annotations}...")
        trf_regions = read_trf_annotations(args.trf_annotations)
        
        # Find overlaps with TRF annotations
        print("Finding overlaps with TRF regions...")
        organized_regions = find_trf_overlaps(organized_regions, trf_regions)
        
        # Print overlap statistics
        with_trf = sum(1 for r in organized_regions if r['trf_overlaps'] > 0)
        print(f"  {with_trf}/{len(organized_regions)} regions overlap with TRF annotations")
        
        # Print repeat class distribution
        class_counts = defaultdict(int)
        for r in organized_regions:
            if 'repeat_class' in r:
                class_counts[r['repeat_class']] += 1
        
        print("  Repeat class distribution:")
        for cls in ['microsatellite', 'minisatellite', 'satellite', 'unknown']:
            if cls in class_counts:
                print(f"    {cls}: {class_counts[cls]}")
    
    # Search for CENP-B boxes if genome provided and finding centromeres
    cenpb_data = {}
    if args.find_centromeres and args.genome:
        # Write regions to temporary BED file for CENP-B search
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp_bed:
            for chrom, chrom_regions in regions.items():
                for start, end, value in chrom_regions:
                    tmp_bed.write(f"{chrom}\t{start}\t{end}\t{value:.2f}\n")
            tmp_bed_path = tmp_bed.name
        
        # Run CENP-B finder
        cenpb_output = f"{args.output_prefix}_cenpb.tsv"
        cenpb_data = run_cenpb_finder(
            tmp_bed_path,
            args.genome,
            cenpb_output,
            pattern=args.cenpb_pattern,
            max_edit_distance=args.cenpb_distance,
            threads=args.cenpb_threads
        )
        
        # Clean up temporary file
        import os
        os.unlink(tmp_bed_path)
        
        if cenpb_data:
            print(f"  Found CENP-B boxes in {len([v for v in cenpb_data.values() if v['count'] > 0])} regions")
    
    # Find centromere candidates if requested
    centromeres = []
    if args.find_centromeres:
        print("\nSearching for centromere candidates...")
        centromeres = find_centromere_candidates(regions, df, args.kmer_size, cenpb_data)
        print(f"  Found {len(centromeres)} potential centromeres")
        
        # Report centromeres with confidence scores
        if centromeres:
            high_conf = sum(1 for c in centromeres if c.get('score', 0) >= 70)
            med_conf = sum(1 for c in centromeres if 50 <= c.get('score', 0) < 70)
            low_conf = sum(1 for c in centromeres if c.get('score', 0) < 50)
            
            print(f"  Confidence levels:")
            if high_conf > 0:
                print(f"    High (score ≥70): {high_conf}")
            if med_conf > 0:
                print(f"    Medium (50-69): {med_conf}")
            if low_conf > 0:
                print(f"    Low (<50): {low_conf}")
            
            if cenpb_data:
                cent_with_cenpb = sum(1 for c in centromeres if 'cenpb_count' in c and c['cenpb_count'] > 0)
                if cent_with_cenpb > 0:
                    avg_cenpb = np.mean([c['cenpb_density'] for c in centromeres if c.get('cenpb_count', 0) > 0])
                    print(f"    With CENP-B boxes: {cent_with_cenpb}/{len(centromeres)} (avg density: {avg_cenpb:.3f}/kb)")
    
    # Write GFF file
    gff_file = f"{args.output_prefix}.gff3"
    write_gff(organized_regions, centromeres, gff_file)
    
    # Extract FASTA sequences if genome provided
    if args.genome:
        fasta_file = f"{args.output_prefix}.fasta"
        extract_fasta_sequences(organized_regions, centromeres, args.genome, fasta_file)
        
        # Extract centromere sequences separately
        if centromeres:
            centromere_fasta = f"{args.output_prefix}_centromeres.fasta"
            extract_centromere_sequences(centromeres, args.genome, centromere_fasta)
    
    # Write summary statistics
    summary_file = f"{args.output_prefix}_summary.txt"
    print(f"\nWriting summary to {summary_file}...")
    
    with open(summary_file, 'w') as f:
        f.write("Low K-mer Diversity Region Detection Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Input file: {args.bed_file}\n")
        f.write(f"K-mer size: {args.kmer_size}\n")
        f.write(f"Percentile threshold: {args.percentile}\n")
        f.write(f"Minimum length: {args.min_length} bp\n")
        f.write(f"Threshold type: {'Global' if args.global_threshold else 'Per-chromosome'}\n\n")
        
        f.write("Low K-mer Diversity Regions by Size:\n")
        f.write("-" * 30 + "\n")
        
        # Group by size
        size_groups = {}
        total_length = 0
        for region in organized_regions:
            size_cat = region['size_category']
            if size_cat not in size_groups:
                size_groups[size_cat] = {'count': 0, 'total_length': 0}
            size_groups[size_cat]['count'] += 1
            size_groups[size_cat]['total_length'] += region['length']
            total_length += region['length']
        
        for size_cat in ['>1Mb', '100kb-1Mb', '10-100kb', '1-10kb', '<1kb']:
            if size_cat in size_groups:
                f.write(f"{size_cat:15s}: {size_groups[size_cat]['count']:5d} regions, "
                       f"{size_groups[size_cat]['total_length']:12,d} bp total\n")
        
        if centromeres:
            cent_length = sum(c['end'] - c['start'] for c in centromeres)
            f.write(f"{'centromere_candidates':20s}: {len(centromeres):5d} regions, {cent_length:12,d} bp total\n")
            
            # Add confidence breakdown
            if any('score' in c for c in centromeres):
                f.write("\n  Centromere Confidence:\n")
                high = sum(1 for c in centromeres if c.get('score', 0) >= 70)
                med = sum(1 for c in centromeres if 50 <= c.get('score', 0) < 70)
                low = sum(1 for c in centromeres if c.get('score', 0) < 50)
                
                if high > 0:
                    f.write(f"    High confidence: {high}\n")
                if med > 0:
                    f.write(f"    Medium confidence: {med}\n")
                if low > 0:
                    f.write(f"    Low confidence: {low}\n")
                
                # Average score
                avg_score = np.mean([c.get('score', 0) for c in centromeres])
                f.write(f"    Average score: {avg_score:.1f}/100\n")
            
            total_length += cent_length
        
        f.write("-" * 30 + "\n")
        f.write(f"{'Total':15s}: {len(organized_regions):5d} regions, {total_length:12,d} bp total\n")
        
        # TRF overlap statistics if available
        if args.trf_annotations:
            f.write("\n\nTRF Annotation Overlap:\n")
            f.write("-" * 30 + "\n")
            
            with_trf = [r for r in organized_regions if r.get('trf_overlaps', 0) > 0]
            f.write(f"Regions with TRF overlap: {len(with_trf)} / {len(organized_regions)}\n")
            
            if with_trf:
                avg_coverage = np.mean([r['trf_coverage'] for r in with_trf])
                f.write(f"Average TRF coverage: {avg_coverage:.2%}\n\n")
                
                f.write("Repeat Classification (based on TRF period):\n")
                class_stats = defaultdict(lambda: {'count': 0, 'total_length': 0})
                for r in organized_regions:
                    if 'repeat_class' in r:
                        class_stats[r['repeat_class']]['count'] += 1
                        class_stats[r['repeat_class']]['total_length'] += r['length']
                
                for cls in ['microsatellite', 'minisatellite', 'satellite', 'unknown']:
                    if cls in class_stats:
                        stats = class_stats[cls]
                        f.write(f"  {cls:15s}: {stats['count']:5d} regions, {stats['total_length']:12,d} bp\n")
                
                # Period distribution
                f.write("\nPeriod Distribution (top 10):\n")
                period_counts = defaultdict(int)
                for r in organized_regions:
                    if r.get('dominant_period'):
                        period_counts[r['dominant_period']] += 1
                
                sorted_periods = sorted(period_counts.items(), key=lambda x: x[1], reverse=True)[:10]
                for period, count in sorted_periods:
                    f.write(f"  Period {period:6d}: {count:5d} regions\n")
        
        # Per-chromosome statistics
        f.write("\n\nPer-Chromosome Statistics:\n")
        f.write("-" * 30 + "\n")
        
        chrom_stats = defaultdict(lambda: {'count': 0, 'length': 0})
        for region in organized_regions:
            chrom_stats[region['chrom']]['count'] += 1
            chrom_stats[region['chrom']]['length'] += region['length']
        
        for chrom in sorted(chrom_stats.keys()):
            stats = chrom_stats[chrom]
            f.write(f"{chrom:15s}: {stats['count']:3d} regions, {stats['length']:10,d} bp\n")
    
    print(f"\nDone! Output files:")
    print(f"  - {gff_file} (GFF3 annotations)")
    if args.genome:
        print(f"  - {fasta_file} (FASTA sequences for all regions)")
        if centromeres:
            centromere_fasta = f"{args.output_prefix}_centromeres.fasta"
            print(f"  - {centromere_fasta} (FASTA sequences for centromeres)")
    print(f"  - {summary_file} (Summary statistics)")


if __name__ == '__main__':
    main()