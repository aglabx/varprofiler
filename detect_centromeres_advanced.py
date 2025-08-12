#!/usr/bin/env python3
"""
Advanced centromere detection using multiple signals:
1. Minimal k-mer diversity (most conserved regions)
2. CENP-B box density (functional centromere markers)
3. Location within large satellite arrays
4. Optional: Integration with ChIP-seq data if available
"""

import pandas as pd
import numpy as np
import argparse
import sys
from pathlib import Path
from collections import defaultdict


def score_centromere_likelihood(region, kmer_percent, cenpb_density, satellite_context):
    """
    Calculate centromere likelihood score based on multiple features.
    
    Args:
        region: Dictionary with region information
        kmer_percent: K-mer diversity percentage
        cenpb_density: CENP-B boxes per kb
        satellite_context: Information about surrounding satellite array
    
    Returns:
        Score from 0 to 100 indicating centromere likelihood
    """
    score = 0
    weights = {}
    
    # 1. K-mer diversity score (lower is better for centromeres)
    # Centromeres typically have <5% k-mer diversity
    if kmer_percent < 2:
        weights['kmer'] = 30
    elif kmer_percent < 5:
        weights['kmer'] = 20
    elif kmer_percent < 10:
        weights['kmer'] = 10
    else:
        weights['kmer'] = 0
    
    # 2. CENP-B box density score (higher is better)
    # Human centromeres typically have 0.1-0.5 CENP-B boxes per kb
    if cenpb_density > 0.3:
        weights['cenpb'] = 30
    elif cenpb_density > 0.15:
        weights['cenpb'] = 20
    elif cenpb_density > 0.05:
        weights['cenpb'] = 10
    else:
        weights['cenpb'] = 0
    
    # 3. Satellite array context (size matters)
    array_size = satellite_context.get('array_size', 0)
    if array_size > 5_000_000:  # >5 Mb array
        weights['context'] = 20
    elif array_size > 1_000_000:  # >1 Mb array
        weights['context'] = 15
    elif array_size > 500_000:  # >500 kb array
        weights['context'] = 10
    elif array_size > 100_000:  # >100 kb array
        weights['context'] = 5
    else:
        weights['context'] = 0
    
    # 4. Position within array (centromeres tend to be central)
    relative_pos = satellite_context.get('relative_position', 0.5)
    if 0.3 <= relative_pos <= 0.7:  # Central 40% of array
        weights['position'] = 10
    elif 0.2 <= relative_pos <= 0.8:  # Central 60% of array
        weights['position'] = 5
    else:
        weights['position'] = 0
    
    # 5. Conservation score (if multiple windows have similar low k-mer diversity)
    if satellite_context.get('is_local_minimum', False):
        weights['conservation'] = 10
    else:
        weights['conservation'] = 0
    
    # Calculate total score
    total_score = sum(weights.values())
    
    # Normalize to 0-100 scale
    max_possible = 100
    normalized_score = (total_score / max_possible) * 100
    
    return normalized_score, weights


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
        merged_arrays = merge_satellite_regions(chrom_regions)
        
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


def merge_satellite_regions(regions, max_gap=500000):
    """
    Merge nearby satellite regions into larger arrays.
    
    Args:
        regions: List of (start, end, value) tuples
        max_gap: Maximum gap to bridge (default: 500kb)
    
    Returns:
        List of merged arrays
    """
    if not regions:
        return []
    
    # Sort by start position
    sorted_regions = sorted(regions, key=lambda x: x[0])
    
    merged = []
    current_start = sorted_regions[0][0]
    current_end = sorted_regions[0][1]
    values = [sorted_regions[0][2]]
    
    for start, end, value in sorted_regions[1:]:
        if start - current_end <= max_gap:
            # Merge
            current_end = max(current_end, end)
            values.append(value)
        else:
            # Save current and start new
            merged.append((current_start, current_end, np.mean(values)))
            current_start = start
            current_end = end
            values = [value]
    
    # Save last
    merged.append((current_start, current_end, np.mean(values)))
    
    return merged


def find_local_minima(windows_df, window_size=5):
    """
    Find local minima in k-mer percentage within windows.
    
    Args:
        windows_df: DataFrame with windows from a satellite array
        window_size: Number of adjacent windows to consider
    
    Returns:
        List of indices that are local minima
    """
    local_minima = []
    kmer_values = windows_df['kmer_percentage'].values
    
    for i in range(len(kmer_values)):
        # Define neighborhood
        start_idx = max(0, i - window_size // 2)
        end_idx = min(len(kmer_values), i + window_size // 2 + 1)
        
        neighborhood = kmer_values[start_idx:end_idx]
        
        # Check if current is minimum in neighborhood
        if len(neighborhood) > 0 and kmer_values[i] == min(neighborhood):
            # Also check if it's significantly lower than neighborhood mean
            neighborhood_mean = np.mean(neighborhood)
            if kmer_values[i] < neighborhood_mean * 0.7:  # At least 30% lower
                local_minima.append(windows_df.index[i])
    
    return local_minima


def integrate_cenpb_data(satellite_arrays, cenpb_file):
    """
    Integrate CENP-B box data with satellite array analysis.
    
    Args:
        satellite_arrays: Dictionary of satellite arrays
        cenpb_file: Path to CENP-B finder output
    
    Returns:
        Updated satellite arrays with CENP-B data
    """
    # Read CENP-B data
    cenpb_df = pd.read_csv(cenpb_file, sep='\t', comment='#',
                          names=['chrom', 'start', 'end', 'kmer_count', 
                                'cenpb_count', 'cenpb_density', 'positions'])
    
    # Create lookup dictionary
    cenpb_lookup = {}
    for _, row in cenpb_df.iterrows():
        key = f"{row['chrom']}:{row['start']}-{row['end']}"
        cenpb_lookup[key] = {
            'count': row['cenpb_count'],
            'density': row['cenpb_density'],
            'positions': row['positions'] if pd.notna(row['positions']) else ''
        }
    
    # Add CENP-B data to windows
    for array_key, array_data in satellite_arrays.items():
        for idx, window in array_data['windows'].iterrows():
            window_key = f"{window['chromosome']}:{int(window['start'])}-{int(window['end'])}"
            if window_key in cenpb_lookup:
                array_data['windows'].loc[idx, 'cenpb_count'] = cenpb_lookup[window_key]['count']
                array_data['windows'].loc[idx, 'cenpb_density'] = cenpb_lookup[window_key]['density']
            else:
                array_data['windows'].loc[idx, 'cenpb_count'] = 0
                array_data['windows'].loc[idx, 'cenpb_density'] = 0.0
    
    return satellite_arrays


def rank_centromere_candidates(satellite_arrays):
    """
    Rank all potential centromeres across the genome.
    
    Args:
        satellite_arrays: Dictionary with satellite array analysis
    
    Returns:
        List of ranked centromere candidates
    """
    candidates = []
    
    for array_key, array_data in satellite_arrays.items():
        chrom = array_data['chrom']
        array_start = array_data['start']
        array_end = array_data['end']
        array_size = array_data['size']
        
        # Check each window in the array
        for idx, window in array_data['windows'].iterrows():
            # Calculate relative position in array
            window_center = (window['start'] + window['end']) / 2
            relative_pos = (window_center - array_start) / array_size
            
            # Check if this is a local minimum
            is_minimum = idx in array_data['local_minima']
            
            # Prepare satellite context
            satellite_context = {
                'array_size': array_size,
                'relative_position': relative_pos,
                'is_local_minimum': is_minimum,
                'array_mean_kmer': array_data['mean_kmer_percent']
            }
            
            # Get CENP-B density
            cenpb_density = window.get('cenpb_density', 0.0)
            
            # Calculate score
            score, weights = score_centromere_likelihood(
                window.to_dict(),
                window['kmer_percentage'],
                cenpb_density,
                satellite_context
            )
            
            candidates.append({
                'chrom': chrom,
                'start': int(window['start']),
                'end': int(window['end']),
                'score': score,
                'kmer_percent': window['kmer_percentage'],
                'cenpb_density': cenpb_density,
                'cenpb_count': window.get('cenpb_count', 0),
                'array_size': array_size,
                'relative_position': relative_pos,
                'is_local_minimum': is_minimum,
                'weights': weights
            })
    
    # Sort by score
    candidates.sort(key=lambda x: x['score'], reverse=True)
    
    return candidates


def select_best_centromeres(candidates, one_per_chrom=True):
    """
    Select the best centromere candidate(s) per chromosome.
    
    Args:
        candidates: List of ranked centromere candidates
        one_per_chrom: If True, select only one per chromosome
    
    Returns:
        List of selected centromeres
    """
    selected = []
    seen_chroms = set()
    
    for candidate in candidates:
        chrom = candidate['chrom']
        
        if one_per_chrom:
            if chrom not in seen_chroms:
                selected.append(candidate)
                seen_chroms.add(chrom)
        else:
            # Select all candidates with score > threshold
            if candidate['score'] > 50:  # Minimum score threshold
                selected.append(candidate)
    
    return selected


def write_results(centromeres, output_file):
    """
    Write centromere predictions to output file.
    
    Args:
        centromeres: List of predicted centromeres
        output_file: Output file path
    """
    with open(output_file, 'w') as f:
        # Write header
        f.write("# Advanced Centromere Predictions\n")
        f.write("# Score: 0-100 likelihood score based on multiple features\n")
        f.write("# Weights show contribution of each feature to the score\n")
        f.write("#chrom\tstart\tend\tscore\tkmer_percent\tcenpb_density\tcenpb_count\t")
        f.write("array_size\trelative_position\tis_local_minimum\tweights\n")
        
        for cent in centromeres:
            # Format weights
            weight_str = ','.join([f"{k}:{v}" for k, v in cent['weights'].items()])
            
            f.write(f"{cent['chrom']}\t{cent['start']}\t{cent['end']}\t")
            f.write(f"{cent['score']:.1f}\t{cent['kmer_percent']:.2f}\t")
            f.write(f"{cent['cenpb_density']:.4f}\t{cent['cenpb_count']}\t")
            f.write(f"{cent['array_size']}\t{cent['relative_position']:.3f}\t")
            f.write(f"{cent['is_local_minimum']}\t{weight_str}\n")


def write_gff(centromeres, output_file):
    """
    Write centromere predictions in GFF3 format.
    
    Args:
        centromeres: List of predicted centromeres
        output_file: Output GFF file path
    """
    with open(output_file, 'w') as f:
        f.write("##gff-version 3\n")
        f.write("##source VarProfiler_Advanced_Centromere_Detection\n")
        f.write("# Centromere predictions based on k-mer diversity and CENP-B box density\n")
        
        for i, cent in enumerate(centromeres, 1):
            attributes = [
                f"ID=centromere_{i}",
                f"Name={cent['chrom']}_centromere",
                f"score={cent['score']:.1f}",
                f"kmer_percent={cent['kmer_percent']:.2f}",
                f"cenpb_density={cent['cenpb_density']:.4f}",
                f"cenpb_count={cent['cenpb_count']}",
                f"array_size={cent['array_size']}",
                f"confidence={'high' if cent['score'] > 70 else 'medium' if cent['score'] > 50 else 'low'}",
                f"color={'255,0,0' if cent['score'] > 70 else '255,128,0' if cent['score'] > 50 else '255,255,0'}"
            ]
            
            f.write(f"{cent['chrom']}\tVarProfiler\tcentromere\t")
            f.write(f"{cent['start']+1}\t{cent['end']}\t{cent['score']:.1f}\t.\t.\t")
            f.write(";".join(attributes) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description='Advanced centromere detection using k-mer diversity and CENP-B boxes',
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument('bed_file',
                       help='BED file with k-mer counts from kmer_profiler')
    parser.add_argument('satellites_file',
                       help='BED or GFF file with satellite/low k-mer diversity regions')
    parser.add_argument('-k', '--kmer-size',
                       type=int, required=True,
                       help='K-mer size used in analysis')
    parser.add_argument('-c', '--cenpb-file',
                       help='TSV file from cenpb_finder with CENP-B box counts')
    parser.add_argument('-o', '--output',
                       default='centromeres_advanced.tsv',
                       help='Output file for centromere predictions')
    parser.add_argument('--gff',
                       help='Also output predictions in GFF3 format')
    parser.add_argument('--all',
                       action='store_true',
                       help='Report all high-scoring candidates, not just best per chromosome')
    parser.add_argument('--min-score',
                       type=float, default=50,
                       help='Minimum score threshold for centromere prediction (0-100)')
    
    args = parser.parse_args()
    
    # Read k-mer data
    print("Reading k-mer profiling data...")
    df = pd.read_csv(args.bed_file, sep='\t', header=None,
                    names=['chromosome', 'start', 'end', 'kmer_count'])
    
    # Read satellite regions
    print("Reading satellite regions...")
    regions = defaultdict(list)
    
    if args.satellites_file.endswith('.gff') or args.satellites_file.endswith('.gff3'):
        # Read GFF format
        with open(args.satellites_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    chrom = parts[0]
                    start = int(parts[3]) - 1  # Convert to 0-based
                    end = int(parts[4])
                    # Try to extract k-mer percentage from attributes
                    kmer_pct = 10.0  # Default
                    if len(parts) >= 9:
                        for attr in parts[8].split(';'):
                            if 'kmer_percent=' in attr:
                                kmer_pct = float(attr.split('=')[1])
                    regions[chrom].append((start, end, kmer_pct))
    else:
        # Read BED format
        sat_df = pd.read_csv(args.satellites_file, sep='\t', header=None)
        for _, row in sat_df.iterrows():
            chrom = row[0]
            start = row[1]
            end = row[2]
            value = row[3] if len(row) > 3 else 10.0
            regions[chrom].append((start, end, value))
    
    # Analyze satellite arrays
    print("Analyzing satellite arrays...")
    satellite_arrays = analyze_satellite_arrays(df, regions, args.kmer_size)
    print(f"  Found {len(satellite_arrays)} satellite arrays")
    
    # Integrate CENP-B data if provided
    if args.cenpb_file:
        print("Integrating CENP-B box data...")
        satellite_arrays = integrate_cenpb_data(satellite_arrays, args.cenpb_file)
        n_with_cenpb = sum(1 for array in satellite_arrays.values() 
                          if 'cenpb_density' in array['windows'].columns and 
                          array['windows']['cenpb_density'].sum() > 0)
        print(f"  {n_with_cenpb} arrays have CENP-B boxes")
    
    # Rank centromere candidates
    print("Ranking centromere candidates...")
    candidates = rank_centromere_candidates(satellite_arrays)
    print(f"  Evaluated {len(candidates)} candidate windows")
    
    # Select best centromeres
    print("Selecting best centromeres...")
    if args.all:
        # Filter by minimum score
        centromeres = [c for c in candidates if c['score'] >= args.min_score]
    else:
        centromeres = select_best_centromeres(candidates, one_per_chrom=True)
    
    print(f"  Selected {len(centromeres)} centromere(s)")
    
    # Report summary
    if centromeres:
        avg_score = np.mean([c['score'] for c in centromeres])
        high_conf = sum(1 for c in centromeres if c['score'] > 70)
        med_conf = sum(1 for c in centromeres if 50 < c['score'] <= 70)
        low_conf = sum(1 for c in centromeres if c['score'] <= 50)
        
        print(f"\nSummary:")
        print(f"  Average score: {avg_score:.1f}")
        print(f"  High confidence (>70): {high_conf}")
        print(f"  Medium confidence (50-70): {med_conf}")
        print(f"  Low confidence (<50): {low_conf}")
        
        if args.cenpb_file:
            with_cenpb = sum(1 for c in centromeres if c['cenpb_count'] > 0)
            avg_cenpb_density = np.mean([c['cenpb_density'] for c in centromeres 
                                        if c['cenpb_density'] > 0]) if with_cenpb > 0 else 0
            print(f"  With CENP-B boxes: {with_cenpb}/{len(centromeres)}")
            if avg_cenpb_density > 0:
                print(f"  Average CENP-B density: {avg_cenpb_density:.3f} per kb")
    
    # Write results
    write_results(centromeres, args.output)
    print(f"\nResults written to: {args.output}")
    
    # Write GFF if requested
    if args.gff:
        write_gff(centromeres, args.gff)
        print(f"GFF written to: {args.gff}")


if __name__ == '__main__':
    main()