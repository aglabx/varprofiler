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
import tempfile
from collections import defaultdict
from pathlib import Path

# Import from our new modules
from varprofiler.utils.utils import (
    read_fasta, 
    write_fasta,
    merge_adjacent_regions,
    find_binary_in_path
)
from varprofiler.utils.gff_utils import (
    write_gff_header,
    write_gff_feature,
    read_gff_annotations,
    filter_gff_by_size
)
from varprofiler.utils.centromere_utils import (
    score_centromere_likelihood,
    find_local_minima,
    run_cenpb_finder,
    select_best_centromeres,
    classify_confidence_level,
    get_centromere_statistics
)


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
    
    # Filter out gap regions (k-mer count = 0 indicates gaps filled with Ns)
    non_gap_df = df[df['kmer_count'] > 0]
    
    if len(non_gap_df) == 0:
        print("Warning: No non-gap regions found in the data")
        return 0, {}
    
    # Calculate global threshold from non-gap regions only
    global_threshold = np.percentile(non_gap_df['kmer_percentage'], percentile)
    
    # Calculate per-chromosome thresholds from non-gap regions only
    chrom_thresholds = {}
    for chrom in df['chromosome'].unique():
        chrom_non_gap = df[(df['chromosome'] == chrom) & (df['kmer_count'] > 0)]
        if len(chrom_non_gap) > 0:
            chrom_thresholds[chrom] = np.percentile(chrom_non_gap['kmer_percentage'], percentile)
        else:
            chrom_thresholds[chrom] = global_threshold  # Use global if no non-gap regions
    
    return global_threshold, chrom_thresholds


def detect_low_variability_regions(bed_file, kmer_size, output_prefix, 
                                  percentile=5, min_length=10000,
                                  use_global_threshold=False):
    """
    Detect regions with low k-mer variability from BED file.
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
        
        # Find windows below threshold, excluding gaps
        low_var = chrom_df[(chrom_df['kmer_percentage'] < threshold) & (chrom_df['kmer_count'] > 0)]
        
        if len(low_var) > 0:
            # Extract regions
            regions = [(row['start'], row['end'], row['kmer_percentage']) 
                      for _, row in low_var.iterrows()]
            
            # Merge adjacent windows
            merged = merge_adjacent_regions(regions)
            
            # Filter by minimum length
            filtered = [(start, end, value) for start, end, value in merged 
                       if end - start >= min_length]
            
            all_regions[chrom] = filtered
            print(f"    Found {len(filtered)} regions >= {min_length}bp")
    
    return all_regions, df


def organize_regions_by_size(regions):
    """
    Organize low variability regions by their size for reporting.
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
    """Get size category for a region."""
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


def find_centromere_candidates(regions, df, kmer_size, cenpb_data=None, 
                              all_candidates=False, min_score=30):
    """
    Find potential centromere-containing regions using multi-signal approach.
    """
    centromeres = []
    
    for chrom, chrom_regions in regions.items():
        # Merge nearby regions into satellite arrays
        merged_arrays = merge_adjacent_regions(chrom_regions, max_gap=500000)
        
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
            
            # Find local minima, excluding gaps
            non_gap_windows = array_windows[array_windows['kmer_count'] > 0]
            if len(non_gap_windows) == 0:
                continue
            
            # Score each window in the array
            window_scores = []
            for idx, window in array_windows.iterrows():
                # Skip gap regions
                if window['kmer_count'] == 0:
                    continue
                    
                # Get CENP-B data
                window_key = f"{chrom}:{int(window['start'])}-{int(window['end'])}"
                cenpb_density = 0
                cenpb_count = 0
                if cenpb_data and window_key in cenpb_data:
                    cenpb_density = cenpb_data[window_key]['density']
                    cenpb_count = cenpb_data[window_key]['count']
                
                # Skip windows without CENP-B boxes
                if cenpb_density == 0:
                    continue
                
                # Check if this is a local minimum
                is_minimum = window['kmer_percentage'] < array_mean * 0.7
                
                # Calculate score using utility function
                score, weights = score_centromere_likelihood(
                    window['kmer_percentage'],
                    cenpb_density,
                    array_size,
                    is_minimum
                )
                
                if score >= min_score:
                    window_scores.append({
                        'window': window,
                        'score': score,
                        'weights': weights,
                        'cenpb_density': cenpb_density,
                        'cenpb_count': cenpb_count,
                        'is_minimum': is_minimum,
                        'array_start': array_start,
                        'array_end': array_end,
                        'array_mean': array_mean,
                        'array_size': array_size
                    })
            
            # Select best scoring window(s) from this array
            window_scores.sort(key=lambda x: x['score'], reverse=True)
            
            if all_candidates:
                # Take all windows that meet the threshold
                for ws in window_scores:
                    window = ws['window']
                    cent_info = {
                        'chrom': chrom,
                        'start': int(window['start']),
                        'end': int(window['end']),
                        'kmer_percentage': window['kmer_percentage'],
                        'satellite_start': ws['array_start'],
                        'satellite_end': ws['array_end'],
                        'satellite_mean': ws['array_mean'],
                        'array_size': ws['array_size'],
                        'score': ws['score'],
                        'weights': ws['weights'],
                        'cenpb_count': ws['cenpb_count'],
                        'cenpb_density': ws['cenpb_density'],
                        'is_local_minimum': ws['is_minimum']
                    }
                    centromeres.append(cent_info)
            else:
                # Take only the best per array
                if window_scores:
                    ws = window_scores[0]
                    window = ws['window']
                    cent_info = {
                        'chrom': chrom,
                        'start': int(window['start']),
                        'end': int(window['end']),
                        'kmer_percentage': window['kmer_percentage'],
                        'satellite_start': ws['array_start'],
                        'satellite_end': ws['array_end'],
                        'satellite_mean': ws['array_mean'],
                        'array_size': ws['array_size'],
                        'score': ws['score'],
                        'weights': ws['weights'],
                        'cenpb_count': ws['cenpb_count'],
                        'cenpb_density': ws['cenpb_density'],
                        'is_local_minimum': ws['is_minimum']
                    }
                    centromeres.append(cent_info)
    
    # Sort by score
    centromeres.sort(key=lambda x: x['score'], reverse=True)
    
    if not all_candidates:
        # Keep only best per chromosome
        centromeres = select_best_centromeres(centromeres, one_per_chrom=True, min_score=min_score)
    
    return centromeres


def write_gff(regions, centromeres, output_file, source='VarProfiler'):
    """Write low k-mer variability regions to GFF3 format file."""
    print(f"\nWriting GFF annotations to {output_file}...")
    
    with open(output_file, 'w') as f:
        # Write header
        write_gff_header(f, source=source, version=3)
        f.write("# Regions with low k-mer diversity detected by VarProfiler\n")
        
        # Write low variability regions
        for i, region in enumerate(regions, 1):
            # Determine color based on size
            if region['length'] >= 1000000:
                color = '255,0,0'  # Red for very large
            elif region['length'] >= 100000:
                color = '255,128,0'  # Orange for large
            elif region['length'] >= 10000:
                color = '255,255,0'  # Yellow for medium
            else:
                color = '128,128,128'  # Gray for small
            
            attributes = {
                'ID': f"low_kmer_div_{i}",
                'Name': f"low_kmer_div_{region['chrom']}_{region['start']}_{region['end']}",
                'length': str(region['length']),
                'size_range': region['size_category'],
                'mean_kmer_percent': f"{region['mean_kmer_percentage']:.2f}",
                'color': color
            }
            
            write_gff_feature(f, region['chrom'], source, 'low_kmer_diversity_region',
                            region['start'], region['end'], 
                            score=f"{region['mean_kmer_percentage']:.2f}",
                            attributes=attributes)
        
        # Write centromere candidates
        for i, cent in enumerate(centromeres, 1):
            confidence, color = classify_confidence_level(cent.get('score', 0))
            
            attributes = {
                'ID': f"centromere_containing_region_{i}",
                'Name': f"centromere_region_{cent['chrom']}_{cent['start']}_{cent['end']}",
                'score': str(cent.get('score', 0)),
                'confidence': confidence,
                'kmer_percent': f"{cent['kmer_percentage']:.2f}",
                'satellite_region': f"{cent['satellite_start']+1}..{cent['satellite_end']}",
                'array_size': str(cent.get('array_size', 0)),
                'color': color
            }
            
            # Add CENP-B info if available
            if 'cenpb_count' in cent:
                attributes['cenpb_count'] = str(cent['cenpb_count'])
                attributes['cenpb_density'] = f"{cent['cenpb_density']:.3f}"
            
            write_gff_feature(f, cent['chrom'], source, 'centromere_containing_region',
                            cent['start'], cent['end'],
                            score=str(cent.get('score', 0)),
                            attributes=attributes)
    
    print(f"  Wrote {len(regions) + len(centromeres)} features")


def extract_sequences(regions, centromeres, genome_file, output_prefix):
    """Extract FASTA sequences for detected regions."""
    if not genome_file:
        return
        
    print(f"\nExtracting sequences from {genome_file}...")
    
    # Read genome
    genome = read_fasta(genome_file)
    if genome is None:
        print(f"  Error: Could not read genome file")
        return
    
    # Extract low variability regions
    sequences = []
    for region in regions:
        chrom = region['chrom']
        if chrom not in genome:
            continue
        
        start = region['start']
        end = region['end']
        seq = genome[chrom][start:end] if end <= len(genome[chrom]) else genome[chrom][start:]
        
        header = (f"low_kmer_div_{chrom}:{start+1}-{end}|"
                 f"length={region['length']}|"
                 f"kmer_pct={region['mean_kmer_percentage']:.2f}|"
                 f"size={region['size_category']}")
        sequences.append((header, seq))
    
    # Write all sequences
    if sequences:
        output_file = f"{output_prefix}.fasta"
        write_fasta(sequences, output_file)
        print(f"  Extracted {len(sequences)} sequences to {output_file}")
    
    # Extract centromere sequences separately
    if centromeres:
        cent_sequences = []
        for i, cent in enumerate(centromeres, 1):
            chrom = cent['chrom']
            if chrom not in genome:
                continue
            
            start = cent['start']
            end = cent['end']
            seq = genome[chrom][start:end] if end <= len(genome[chrom]) else genome[chrom][start:]
            
            score = cent.get('score', 0)
            confidence, _ = classify_confidence_level(score)
            
            header = (f"centromere_{i} chr={chrom} pos={start+1}-{end} "
                     f"length={end-start}bp score={score:.1f} confidence={confidence.upper()} "
                     f"kmer_diversity={cent['kmer_percentage']:.2f}%")
            
            if 'cenpb_count' in cent:
                header += f" cenpb_count={cent['cenpb_count']} cenpb_density={cent['cenpb_density']:.3f}/kb"
            
            cent_sequences.append((header, seq))
        
        if cent_sequences:
            cent_output = f"{output_prefix}_centromeres.fasta"
            write_fasta(cent_sequences, cent_output)
            print(f"  Extracted {len(cent_sequences)} centromere sequences to {cent_output}")


def write_summary(regions, centromeres, output_file, args):
    """Write summary statistics."""
    print(f"\nWriting summary to {output_file}...")
    
    with open(output_file, 'w') as f:
        f.write("Low K-mer Diversity Region Detection Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Input file: {args.bed_file}\n")
        f.write(f"K-mer size: {args.kmer_size}\n")
        f.write(f"Percentile threshold: {args.percentile}\n")
        f.write(f"Minimum length: {args.min_length} bp\n")
        f.write(f"Threshold type: {'Global' if args.global_threshold else 'Per-chromosome'}\n\n")
        
        # Group by size
        size_groups = defaultdict(lambda: {'count': 0, 'total_length': 0})
        total_length = 0
        
        for region in regions:
            size_cat = region['size_category']
            size_groups[size_cat]['count'] += 1
            size_groups[size_cat]['total_length'] += region['length']
            total_length += region['length']
        
        f.write("Low K-mer Diversity Regions by Size:\n")
        f.write("-" * 30 + "\n")
        
        for size_cat in ['>1Mb', '100kb-1Mb', '10-100kb', '1-10kb', '<1kb']:
            if size_cat in size_groups:
                f.write(f"{size_cat:15s}: {size_groups[size_cat]['count']:5d} regions, "
                       f"{size_groups[size_cat]['total_length']:12,d} bp total\n")
        
        if centromeres:
            # Get centromere statistics
            stats = get_centromere_statistics(centromeres)
            cent_length = sum(c['end'] - c['start'] for c in centromeres)
            
            f.write(f"{'centromere_candidates':20s}: {len(centromeres):5d} regions, "
                   f"{cent_length:12,d} bp total\n")
            
            f.write("\n  Centromere Statistics:\n")
            f.write(f"    Average score: {stats.get('avg_score', 0):.1f}/100\n")
            
            if stats.get('high_confidence', 0) > 0:
                f.write(f"    High confidence: {stats['high_confidence']}\n")
            if stats.get('medium_confidence', 0) > 0:
                f.write(f"    Medium confidence: {stats['medium_confidence']}\n")
            if stats.get('low_confidence', 0) > 0:
                f.write(f"    Low confidence: {stats['low_confidence']}\n")
            
            if 'avg_cenpb_density' in stats:
                f.write(f"    Average CENP-B density: {stats['avg_cenpb_density']:.3f}/kb\n")
            
            total_length += cent_length
        
        f.write("-" * 30 + "\n")
        f.write(f"{'Total':15s}: {len(regions)} regions, {total_length:12,d} bp total\n")
        
        # Per-chromosome statistics
        f.write("\n\nPer-Chromosome Statistics:\n")
        f.write("-" * 30 + "\n")
        
        chrom_stats = defaultdict(lambda: {'count': 0, 'length': 0})
        for region in regions:
            chrom_stats[region['chrom']]['count'] += 1
            chrom_stats[region['chrom']]['length'] += region['length']
        
        for chrom in sorted(chrom_stats.keys()):
            stats = chrom_stats[chrom]
            f.write(f"{chrom:15s}: {stats['count']:3d} regions, {stats['length']:10,d} bp\n")


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
    parser.add_argument('--all-centromeres',
                       action='store_true',
                       help='Report all centromere candidates, not just best per chromosome')
    parser.add_argument('--centromere-min-score',
                       type=int, default=30,
                       help='Minimum score for centromere candidates (0-100, default: 30)')
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
    
    # Search for CENP-B boxes and find centromeres if requested
    centromeres = []
    if args.find_centromeres and args.genome:
        print("\nSearching for centromere candidates...")
        
        # Write regions to temporary BED file for CENP-B search
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
        
        # Clean up
        os.unlink(tmp_bed_path)
        
        if cenpb_data:
            print(f"  Found CENP-B boxes in {len([v for v in cenpb_data.values() if v['count'] > 0])} regions")
            
            # Find centromere candidates
            centromeres = find_centromere_candidates(
                regions, df, args.kmer_size, cenpb_data,
                all_candidates=args.all_centromeres,
                min_score=args.centromere_min_score
            )
            
            print(f"  Found {len(centromeres)} potential centromeres")
            
            if centromeres:
                stats = get_centromere_statistics(centromeres)
                if stats:
                    print(f"  Average score: {stats['avg_score']:.1f}/100")
                    if stats.get('high_confidence', 0) > 0:
                        print(f"    High confidence: {stats['high_confidence']}")
                    if stats.get('medium_confidence', 0) > 0:
                        print(f"    Medium confidence: {stats['medium_confidence']}")
        else:
            print("  Warning: CENP-B analysis failed or no CENP-B boxes found")
    
    # Write outputs
    gff_file = f"{args.output_prefix}.gff3"
    write_gff(organized_regions, centromeres, gff_file)
    
    # Extract sequences if genome provided
    if args.genome:
        extract_sequences(organized_regions, centromeres, args.genome, args.output_prefix)
    
    # Write summary
    summary_file = f"{args.output_prefix}_summary.txt"
    write_summary(organized_regions, centromeres, summary_file, args)
    
    print(f"\nDone! Output files:")
    print(f"  - {gff_file} (GFF3 annotations)")
    if args.genome:
        print(f"  - {args.output_prefix}.fasta (FASTA sequences)")
        if centromeres:
            print(f"  - {args.output_prefix}_centromeres.fasta (Centromere sequences)")
    print(f"  - {summary_file} (Summary statistics)")


if __name__ == '__main__':
    main()