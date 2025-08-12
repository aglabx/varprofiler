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


def find_centromere_candidates(regions, df, kmer_size):
    """
    Find potential centromeres as regions with minimal k-mer variability
    within satellite arrays.
    
    Args:
        regions: Dictionary of low variability regions per chromosome
        df: Original dataframe with all windows
        kmer_size: K-mer size
    
    Returns:
        List of potential centromere regions
    """
    centromeres = []
    
    for chrom, chrom_regions in regions.items():
        # Only consider large satellites (potential pericentromeric regions)
        large_regions = [r for r in chrom_regions if r[1] - r[0] >= 100000]
        
        if not large_regions:
            continue
        
        chrom_df = df[df['chromosome'] == chrom].copy()
        if 'kmer_percentage' not in chrom_df.columns:
            chrom_df['window_size'] = chrom_df['end'] - chrom_df['start']
            chrom_df['max_kmers'] = chrom_df['window_size'] - kmer_size + 1
            chrom_df['kmer_percentage'] = (chrom_df['kmer_count'] / chrom_df['max_kmers']) * 100
        
        for start, end, mean_value in large_regions:
            # Find windows within this region
            region_windows = chrom_df[(chrom_df['start'] >= start) & 
                                     (chrom_df['end'] <= end)]
            
            if len(region_windows) > 0:
                # Find minimum variability window
                min_idx = region_windows['kmer_percentage'].idxmin()
                min_window = region_windows.loc[min_idx]
                
                # Check if this is significantly lower than region mean
                if min_window['kmer_percentage'] < mean_value * 0.5:  # At least 50% lower
                    centromeres.append({
                        'chrom': chrom,
                        'start': min_window['start'],
                        'end': min_window['end'],
                        'kmer_percentage': min_window['kmer_percentage'],
                        'satellite_start': start,
                        'satellite_end': end,
                        'satellite_mean': mean_value
                    })
    
    return centromeres


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
                f"mean_kmer_percent={region['mean_kmer_percentage']:.2f}",
                f"color={color}"
            ]
            
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
            
            attributes = [
                f"ID=centromere_candidate_{feature_id}",
                f"Name=centromere_{chrom}_{start}_{end}",
                f"kmer_percent={cent['kmer_percentage']:.2f}",
                f"satellite_region={cent['satellite_start']+1}..{cent['satellite_end']}",
                f"satellite_mean={cent['satellite_mean']:.2f}",
                f"color=0,0,255"  # Blue for centromeres
            ]
            
            f.write(f"{chrom}\t{source}\tcentromere\t{start}\t{end}\t"
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
    
    # Find centromere candidates if requested
    centromeres = []
    if args.find_centromeres:
        print("\nSearching for centromere candidates...")
        centromeres = find_centromere_candidates(regions, df, args.kmer_size)
        print(f"  Found {len(centromeres)} potential centromeres")
    
    # Write GFF file
    gff_file = f"{args.output_prefix}.gff3"
    write_gff(organized_regions, centromeres, gff_file)
    
    # Extract FASTA sequences if genome provided
    if args.genome:
        fasta_file = f"{args.output_prefix}.fasta"
        extract_fasta_sequences(organized_regions, centromeres, args.genome, fasta_file)
    
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
            total_length += cent_length
        
        f.write("-" * 30 + "\n")
        f.write(f"{'Total':15s}: {len(organized_regions):5d} regions, {total_length:12,d} bp total\n")
        
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
        print(f"  - {fasta_file} (FASTA sequences)")
    print(f"  - {summary_file} (Summary statistics)")


if __name__ == '__main__':
    main()