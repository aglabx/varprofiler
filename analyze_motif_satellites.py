#!/usr/bin/env python3
"""
Analyze motif enrichment in satellite DNA vs whole genome.

This script compares motif frequencies from motif_discovery output
with their occurrence in annotated satellite DNA regions.
"""

import argparse
import pandas as pd
import numpy as np
from collections import defaultdict
import sys
import re


def reverse_complement(seq):
    """Get reverse complement of DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 
                  'N': 'N', 'R': 'Y', 'Y': 'R', 'S': 'S', 
                  'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V',
                  'V': 'B', 'D': 'H', 'H': 'D'}
    return ''.join(complement.get(c.upper(), c) for c in seq[::-1])


def count_motif_in_sequence(motif, sequence, both_strands=True):
    """Count occurrences of motif in sequence."""
    seq_upper = sequence.upper()
    motif_upper = motif.upper()
    
    # Count forward strand
    count = seq_upper.count(motif_upper)
    
    # Count reverse strand if requested
    if both_strands:
        rc_motif = reverse_complement(motif_upper)
        if rc_motif != motif_upper:  # Avoid double counting palindromes
            count += seq_upper.count(rc_motif)
    
    return count


def count_motif_with_regex(motif, sequence, both_strands=True):
    """Count motif allowing for IUPAC ambiguity codes."""
    # Convert IUPAC to regex
    iupac_dict = {
        'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
        'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
        'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
    }
    
    pattern = motif.upper()
    for code, regex in iupac_dict.items():
        pattern = pattern.replace(code, regex)
    
    seq_upper = sequence.upper()
    
    # Count forward strand
    count = len(re.findall(pattern, seq_upper))
    
    # Count reverse strand if requested
    if both_strands:
        rc_motif = reverse_complement(motif.upper())
        rc_pattern = rc_motif
        for code, regex in iupac_dict.items():
            rc_pattern = rc_pattern.replace(code, regex)
        
        if rc_pattern != pattern:  # Avoid double counting palindromes
            count += len(re.findall(rc_pattern, seq_upper))
    
    return count


def load_motif_data(tsv_file):
    """Load motif discovery results from TSV file."""
    print(f"Loading motif data from {tsv_file}...")
    
    # Read TSV with proper column names
    df = pd.read_csv(tsv_file, sep='\t')
    
    # Expected columns from motif_discovery output
    expected_cols = ['motif', 'alignment', 'edit_distance', 'cigar', 
                     'count_forward', 'count_reverse', 'count_total']
    
    # Check if all expected columns are present
    missing_cols = [col for col in expected_cols if col not in df.columns]
    if missing_cols:
        print(f"Warning: Missing columns in motif file: {missing_cols}")
    
    print(f"  Loaded {len(df)} motif variants")
    print(f"  Total occurrences in genome: {df['count_total'].sum():,}")
    
    return df


def load_satellite_data(sat_file):
    """Load satellite DNA annotations."""
    print(f"Loading satellite annotations from {sat_file}...")
    
    # Read TSV file
    df = pd.read_csv(sat_file, sep='\t')
    
    # Check for required column
    if 'trf_array' not in df.columns:
        raise ValueError("Column 'trf_array' not found in satellite file")
    
    # Filter out rows without sequences
    df = df[df['trf_array'].notna()]
    df = df[df['trf_array'] != '']
    
    print(f"  Loaded {len(df)} satellite arrays")
    
    # Calculate total satellite DNA length
    total_sat_bp = df['trf_array'].str.len().sum()
    print(f"  Total satellite DNA: {total_sat_bp/1e6:.1f} Mb")
    
    return df


def analyze_motif_enrichment(motif_df, sat_df, output_file):
    """Analyze motif enrichment in satellites vs genome."""
    
    results = []
    
    print("\nAnalyzing motif enrichment in satellite DNA...")
    
    for idx, row in motif_df.iterrows():
        motif = row['motif']
        genome_count = row['count_total']
        
        # Count in all satellite arrays
        sat_count = 0
        for sat_seq in sat_df['trf_array']:
            if pd.notna(sat_seq):
                sat_count += count_motif_in_sequence(motif, sat_seq, both_strands=True)
        
        # Calculate enrichment
        # Note: This is raw count, could normalize by sequence length
        if genome_count > 0:
            sat_fraction = sat_count / genome_count
        else:
            sat_fraction = 0
        
        results.append({
            'motif': motif,
            'alignment': row.get('alignment', ''),
            'edit_distance': row.get('edit_distance', 0),
            'genome_total': genome_count,
            'satellite_total': sat_count,
            'satellite_fraction': sat_fraction,
            'enrichment': 'satellite' if sat_fraction > 0.5 else 'genome-wide'
        })
        
        # Progress indicator
        if (idx + 1) % 100 == 0:
            print(f"  Processed {idx + 1}/{len(motif_df)} motifs...")
    
    # Create results dataframe
    results_df = pd.DataFrame(results)
    
    # Sort by satellite fraction (most satellite-specific first)
    results_df = results_df.sort_values('satellite_fraction', ascending=False)
    
    # Save results
    print(f"\nSaving results to {output_file}...")
    results_df.to_csv(output_file, sep='\t', index=False)
    
    # Print summary statistics
    print("\n=== SUMMARY ===")
    print(f"Total motifs analyzed: {len(results_df)}")
    
    # Highly satellite-enriched (>80% in satellites)
    highly_enriched = results_df[results_df['satellite_fraction'] > 0.8]
    print(f"Highly satellite-enriched (>80%): {len(highly_enriched)}")
    
    # Moderately enriched (50-80%)
    moderately_enriched = results_df[(results_df['satellite_fraction'] > 0.5) & 
                                     (results_df['satellite_fraction'] <= 0.8)]
    print(f"Moderately satellite-enriched (50-80%): {len(moderately_enriched)}")
    
    # Genome-wide (<50% in satellites)
    genome_wide = results_df[results_df['satellite_fraction'] <= 0.5]
    print(f"Genome-wide distribution (<50%): {len(genome_wide)}")
    
    # Top satellite-specific motifs
    print("\n=== TOP 10 SATELLITE-SPECIFIC MOTIFS ===")
    print(f"{'Motif':<20} {'Alignment':<20} {'Genome':<12} {'Satellite':<12} {'Sat%':<8}")
    print("-" * 80)
    
    for _, row in results_df.head(10).iterrows():
        print(f"{row['motif']:<20} {row['alignment']:<20} "
              f"{row['genome_total']:<12,} {row['satellite_total']:<12,} "
              f"{row['satellite_fraction']*100:>6.1f}%")
    
    # Motifs absent from satellites
    absent_from_sat = results_df[results_df['satellite_total'] == 0]
    if len(absent_from_sat) > 0:
        print(f"\n=== MOTIFS ABSENT FROM SATELLITES ({len(absent_from_sat)} total) ===")
        print("Top 10 by genome frequency:")
        for _, row in absent_from_sat.head(10).iterrows():
            print(f"  {row['motif']:<20} (genome count: {row['genome_total']:,})")
    
    return results_df


def analyze_by_satellite_type(motif_df, sat_df, output_prefix):
    """Analyze motif distribution by satellite type/family."""
    
    # Check if satellite type columns exist
    type_columns = ['trf_type', 'trf_family', 'trf_superfamily']
    available_cols = [col for col in type_columns if col in sat_df.columns]
    
    if not available_cols:
        print("No satellite type columns found for detailed analysis")
        return
    
    for type_col in available_cols:
        print(f"\n=== Analyzing by {type_col} ===")
        
        # Group satellites by type
        sat_types = sat_df[type_col].value_counts()
        print(f"Found {len(sat_types)} different {type_col} categories")
        
        # Analyze top types
        results = []
        for sat_type in sat_types.head(10).index:
            if pd.isna(sat_type):
                continue
                
            # Get sequences for this type
            type_seqs = sat_df[sat_df[type_col] == sat_type]['trf_array']
            total_bp = type_seqs.str.len().sum()
            
            print(f"\n  {sat_type}: {len(type_seqs)} arrays, {total_bp/1e6:.1f} Mb")
            
            # Count each motif in this satellite type
            type_motif_counts = {}
            for motif in motif_df['motif'].head(100):  # Analyze top 100 motifs
                count = 0
                for seq in type_seqs:
                    if pd.notna(seq):
                        count += count_motif_in_sequence(motif, seq)
                if count > 0:
                    type_motif_counts[motif] = count
            
            # Store results
            for motif, count in type_motif_counts.items():
                results.append({
                    'satellite_type': sat_type,
                    'type_column': type_col,
                    'motif': motif,
                    'count': count,
                    'density_per_kb': (count * 1000) / total_bp if total_bp > 0 else 0
                })
        
        # Create results dataframe
        if results:
            type_df = pd.DataFrame(results)
            type_df = type_df.sort_values(['satellite_type', 'density_per_kb'], ascending=[True, False])
            
            # Save detailed results
            output_file = f"{output_prefix}_{type_col}_analysis.tsv"
            type_df.to_csv(output_file, sep='\t', index=False)
            print(f"  Saved detailed analysis to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze motif enrichment in satellite DNA vs whole genome')
    parser.add_argument('motif_tsv', 
                       help='TSV file from motif_discovery tool')
    parser.add_argument('satellite_tsv',
                       help='TSV file with satellite DNA annotations (must have trf_array column)')
    parser.add_argument('-o', '--output', default='motif_satellite_enrichment.tsv',
                       help='Output file for enrichment analysis')
    parser.add_argument('--by-type', action='store_true',
                       help='Also analyze by satellite type/family')
    parser.add_argument('--top-motifs', type=int, default=0,
                       help='Only analyze top N motifs by frequency (0=all)')
    
    args = parser.parse_args()
    
    # Load data
    motif_df = load_motif_data(args.motif_tsv)
    sat_df = load_satellite_data(args.satellite_tsv)
    
    # Filter to top motifs if requested
    if args.top_motifs > 0:
        print(f"\nAnalyzing top {args.top_motifs} motifs only...")
        motif_df = motif_df.head(args.top_motifs)
    
    # Main enrichment analysis
    results_df = analyze_motif_enrichment(motif_df, sat_df, args.output)
    
    # Additional analysis by satellite type
    if args.by_type:
        output_prefix = args.output.replace('.tsv', '')
        analyze_by_satellite_type(motif_df, sat_df, output_prefix)
    
    print("\nAnalysis complete!")


if __name__ == '__main__':
    main()