#!/usr/bin/env python3
"""
Visualize motif occurrences in monomer sequences by converting them to lowercase.

This script takes monomer decomposition data and highlights motif occurrences
by converting matched positions to lowercase in the sequences.
"""

import argparse
import pandas as pd
import sys

# Import from our new modules
from varprofiler.utils.utils import get_tqdm, write_fasta
from varprofiler.utils.motif_utils import (
    find_motif_positions,
    mark_motifs_in_sequence,
    load_motifs_from_file,
    calculate_motif_statistics
)

# Get tqdm with fallback
tqdm = get_tqdm()


def mark_motifs_with_stats(sequence, motifs):
    """
    Mark motif occurrences and return statistics.
    
    Args:
        sequence: DNA sequence
        motifs: List of motif sequences
    
    Returns:
        Tuple of (marked sequence, marked positions count, total motif occurrences)
    """
    # Use utility function to mark motifs
    marked_seq = mark_motifs_in_sequence(sequence, motifs)
    
    # Count marked positions and total motifs
    marked_count = sum(1 for c in marked_seq if c.islower())
    
    # Count total motif occurrences
    total_motif_count = 0
    for motif in motifs:
        positions = find_motif_positions(sequence, motif, both_strands=True)
        total_motif_count += len(positions)
    
    return marked_seq, marked_count, total_motif_count


def process_monomer_file(monomer_file, motifs, output_file, 
                         stats_file=None, only_monomers=False):
    """
    Process monomer decomposition file and mark motifs.
    
    Args:
        monomer_file: TSV file with monomer decomposition
        motifs: List of motif sequences to mark
        output_file: Output TSV file with marked sequences
        stats_file: Optional file to save statistics
        only_monomers: Process only MONOMER type, skip flanks
    """
    print(f"Reading monomer data from {monomer_file}...")
    
    # Read TSV file
    df = pd.read_csv(monomer_file, sep='\t')
    
    # Check for required columns
    required_cols = ['sequence_id', 'orientation', 'index', 'type', 'length', 'sequence']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    # Filter if requested
    if only_monomers:
        print("  Filtering to MONOMER type only...")
        df = df[df['type'] == 'MONOMER']
    
    print(f"  Processing {len(df)} sequences...")
    
    # Process each sequence
    marked_sequences = []
    motif_counts_list = []
    stats = []
    
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Marking motifs"):
        sequence = row['sequence']
        
        # Mark motifs and get statistics
        marked_seq, marked_count, total_motifs = mark_motifs_with_stats(sequence, motifs)
        marked_sequences.append(marked_seq)
        motif_counts_list.append(total_motifs)
        
        # Calculate coverage
        coverage = (marked_count / len(sequence)) * 100 if len(sequence) > 0 else 0
        
        # Build statistics record
        stats_record = {
            'sequence_id': row['sequence_id'],
            'orientation': row['orientation'],
            'index': row['index'],
            'type': row['type'],
            'length': row['length'],
            'marked_bp': marked_count,
            'coverage_pct': coverage,
            'total_motif_hits': total_motifs
        }
        
        # Add individual motif counts for top motifs
        for i, motif in enumerate(motifs[:10]):  # Count top 10 motifs
            positions = find_motif_positions(sequence, motif, both_strands=True)
            stats_record[f'motif_{i+1}_count'] = len(positions)
        
        stats.append(stats_record)
    
    # Add marked sequences and motif counts to dataframe
    df['marked_sequence'] = marked_sequences
    df['motif_count'] = motif_counts_list
    
    # Save output
    print(f"\nSaving marked sequences to {output_file}...")
    df.to_csv(output_file, sep='\t', index=False)
    
    # Save statistics if requested
    if stats_file:
        print(f"Saving statistics to {stats_file}...")
        stats_df = pd.DataFrame(stats)
        stats_df.to_csv(stats_file, sep='\t', index=False)
    
    # Print summary
    print_summary(df, stats)
    
    return df, pd.DataFrame(stats)


def print_summary(df, stats):
    """Print summary statistics."""
    print("\n=== SUMMARY ===")
    
    # Group by sequence_id
    stats_df = pd.DataFrame(stats)
    grouped = stats_df.groupby('sequence_id')
    print(f"Total sequence IDs: {len(grouped)}")
    
    # Count sequences with/without motifs
    with_motifs = df[df['motif_count'] > 0]
    without_motifs = df[df['motif_count'] == 0]
    print(f"\nSequences with motifs: {len(with_motifs)} ({len(with_motifs)/len(df)*100:.1f}%)")
    print(f"Sequences without motifs: {len(without_motifs)} ({len(without_motifs)/len(df)*100:.1f}%)")
    
    # Statistics by type
    type_stats = stats_df.groupby('type').agg({
        'coverage_pct': ['mean', 'std', 'min', 'max'],
        'total_motif_hits': ['sum', 'mean']
    }).round(2)
    print("\nMotif coverage by type:")
    print(type_stats)
    
    # Motif count distribution
    print("\nMotif count distribution:")
    motif_dist = df['motif_count'].value_counts().sort_index()
    for count in range(min(6, motif_dist.index.max() + 1) if len(motif_dist) > 0 else 0):
        n_seqs = motif_dist.get(count, 0)
        print(f"  {count} motifs: {n_seqs} sequences ({n_seqs/len(df)*100:.1f}%)")
    if len(motif_dist) > 0 and motif_dist.index.max() > 5:
        n_seqs = motif_dist[motif_dist.index > 5].sum()
        print(f"  >5 motifs: {n_seqs} sequences ({n_seqs/len(df)*100:.1f}%)")
    
    # Find monomers with highest motif density
    monomer_stats = stats_df[stats_df['type'] == 'MONOMER']
    if len(monomer_stats) > 0:
        top_monomers = monomer_stats.nlargest(10, 'coverage_pct')
        print("\nTop 10 monomers by motif coverage:")
        print(top_monomers[['sequence_id', 'index', 'length', 'coverage_pct', 'total_motif_hits']])


def generate_fasta(df, output_file, min_coverage=0):
    """
    Generate FASTA file from marked monomer data.
    
    Args:
        df: DataFrame with marked sequences
        output_file: Output FASTA file
        min_coverage: Minimum coverage to include in FASTA
    """
    print(f"\nGenerating FASTA file: {output_file}")
    
    fasta_sequences = []
    
    for idx, row in df.iterrows():
        # Calculate coverage if needed
        marked_seq = row.get('marked_sequence', row['sequence'])
        marked_count = sum(1 for c in marked_seq if c.islower())
        coverage = (marked_count / len(marked_seq)) * 100 if len(marked_seq) > 0 else 0
        
        if coverage < min_coverage:
            continue
        
        # Create header
        header_parts = [
            f"{row['sequence_id']}_{row['type']}_{row['index']}",
            f"ori={row['orientation']}",
            f"len={row['length']}",
            f"motif_coverage={coverage:.1f}%"
        ]
        
        header = ' '.join(header_parts)
        fasta_sequences.append((header, marked_seq))
    
    # Write using utility function
    write_fasta(fasta_sequences, output_file, wrap_length=80)
    print(f"  FASTA file written with {len(fasta_sequences)} sequences")


def main():
    parser = argparse.ArgumentParser(
        description='Visualize motif occurrences in monomer sequences')
    
    parser.add_argument('monomer_file',
                       help='TSV file with monomer decomposition')
    parser.add_argument('motif_file',
                       help='Motif file: TSV from motif_discovery OR text file with motifs')
    parser.add_argument('-o', '--output', required=True,
                       help='Output TSV file with marked sequences')
    parser.add_argument('--fasta', 
                       help='Also generate FASTA output')
    parser.add_argument('--stats',
                       help='Save statistics to file')
    parser.add_argument('--top-motifs', type=int, default=None,
                       help='Use only top N motifs')
    parser.add_argument('--min-count', type=int, default=0,
                       help='Minimum motif count to include (for TSV input)')
    parser.add_argument('--single-motif',
                       help='Mark single motif instead of reading from file')
    parser.add_argument('--only-monomers', action='store_true',
                       help='Process only MONOMER type sequences')
    parser.add_argument('--min-coverage', type=float, default=0,
                       help='Minimum coverage %% for FASTA output')
    
    args = parser.parse_args()
    
    # Load motifs
    if args.single_motif:
        # Use single motif provided on command line
        motifs = [args.single_motif]
        print(f"Using single motif: {args.single_motif}")
    else:
        # Load from file using utility function
        motifs = load_motifs_from_file(args.motif_file, 
                                       top_n=args.top_motifs, 
                                       min_count=args.min_count)
        print(f"Loaded {len(motifs)} motifs")
    
    # Process monomer file
    df, stats = process_monomer_file(
        args.monomer_file, 
        motifs, 
        args.output,
        args.stats,
        args.only_monomers
    )
    
    # Generate FASTA if requested
    if args.fasta:
        generate_fasta(df, args.fasta, args.min_coverage)
    
    print("\nVisualization complete!")
    print("Motif occurrences are shown in lowercase in the marked_sequence column.")


if __name__ == '__main__':
    main()