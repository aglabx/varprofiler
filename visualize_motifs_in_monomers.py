#!/usr/bin/env python3
"""
Visualize motif occurrences in monomer sequences by converting them to lowercase.

This script takes monomer decomposition data and highlights motif occurrences
by converting matched positions to lowercase in the sequences.
"""

import argparse
import pandas as pd
import re
import sys
from collections import defaultdict

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, total=None, desc=None):
        if desc:
            print(f"{desc}...")
        return iterable


def reverse_complement(seq):
    """Get reverse complement of DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 
                  'N': 'N', 'R': 'Y', 'Y': 'R', 'S': 'S', 
                  'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V',
                  'V': 'B', 'D': 'H', 'H': 'D'}
    return ''.join(complement.get(c.upper(), c) for c in seq[::-1])


def find_motif_positions(sequence, motif, both_strands=True):
    """
    Find all positions of a motif in a sequence.
    Returns list of (start, end) positions.
    """
    positions = []
    seq_upper = sequence.upper()
    motif_upper = motif.upper()
    motif_len = len(motif_upper)
    
    # Find forward strand matches
    start = 0
    while True:
        pos = seq_upper.find(motif_upper, start)
        if pos == -1:
            break
        positions.append((pos, pos + motif_len))
        start = pos + 1
    
    # Find reverse strand matches if requested
    if both_strands:
        rc_motif = reverse_complement(motif_upper)
        if rc_motif != motif_upper:  # Avoid double marking palindromes
            start = 0
            while True:
                pos = seq_upper.find(rc_motif, start)
                if pos == -1:
                    break
                positions.append((pos, pos + motif_len))
                start = pos + 1
    
    return positions


def mark_motifs_in_sequence(sequence, motifs):
    """
    Mark motif occurrences in sequence by converting to lowercase.
    
    Args:
        sequence: DNA sequence
        motifs: List of motif sequences
    
    Returns:
        Tuple of (modified sequence, marked positions count, total motif occurrences)
    """
    if isinstance(motifs, str):
        motifs = [motifs]
    
    # Create a position mask for lowercase conversion
    seq_len = len(sequence)
    lowercase_mask = [False] * seq_len
    
    # Find all motif positions and count occurrences
    total_motif_count = 0
    for motif in motifs:
        positions = find_motif_positions(sequence, motif, both_strands=True)
        total_motif_count += len(positions)
        for start, end in positions:
            for i in range(start, min(end, seq_len)):
                lowercase_mask[i] = True
    
    # Apply the mask to create the output sequence
    result = []
    marked_count = 0
    for i, char in enumerate(sequence):
        if lowercase_mask[i]:
            result.append(char.lower())
            marked_count += 1
        else:
            result.append(char.upper())
    
    return ''.join(result), marked_count, total_motif_count


def load_motifs(motif_file, top_n=None, min_count=0):
    """
    Load motifs from motif_discovery output file or text file.
    
    Args:
        motif_file: TSV file from motif_discovery or text file with motifs
        top_n: Load only top N motifs
        min_count: Minimum total count to include motif
    
    Returns:
        List of motif sequences
    """
    print(f"Loading motifs from {motif_file}...")
    
    try:
        # Try to load as TSV from motif_discovery
        df = pd.read_csv(motif_file, sep='\t')
        
        # Check if it has expected columns
        if 'motif' in df.columns:
            # Filter by minimum count if column exists
            if 'count_total' in df.columns and min_count > 0:
                df = df[df['count_total'] >= min_count]
            
            # Get top N if specified
            if top_n:
                df = df.head(top_n)
            
            motifs = df['motif'].tolist()
            print(f"  Loaded {len(motifs)} motifs from TSV")
            return motifs
    except:
        pass
    
    # Load as plain text file with one motif per line
    with open(motif_file, 'r') as f:
        motifs = [line.strip() for line in f if line.strip()]
    
    if top_n:
        motifs = motifs[:top_n]
    
    print(f"  Loaded {len(motifs)} motifs from text file")
    return motifs


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
        
        # Mark motifs
        marked_seq, marked_count, total_motifs = mark_motifs_in_sequence(sequence, motifs)
        marked_sequences.append(marked_seq)
        motif_counts_list.append(total_motifs)
        
        # Calculate statistics
        coverage = (marked_count / len(sequence)) * 100 if len(sequence) > 0 else 0
        
        # Count individual motif occurrences for detailed stats
        motif_counts = {}
        for motif in motifs[:10]:  # Count top 10 motifs individually
            positions = find_motif_positions(sequence, motif, both_strands=True)
            motif_counts[motif] = len(positions)
        
        stats.append({
            'sequence_id': row['sequence_id'],
            'orientation': row['orientation'],
            'index': row['index'],
            'type': row['type'],
            'length': row['length'],
            'marked_bp': marked_count,
            'coverage_pct': coverage,
            'total_motif_hits': total_motifs,
            **{f'motif_{i+1}_count': count for i, count in enumerate(motif_counts.values())}
        })
    
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
    print("\n=== SUMMARY ===")
    
    # Group by sequence_id
    grouped = pd.DataFrame(stats).groupby('sequence_id')
    print(f"Total sequence IDs: {len(grouped)}")
    
    # Count sequences with/without motifs
    with_motifs = df[df['motif_count'] > 0]
    without_motifs = df[df['motif_count'] == 0]
    print(f"\nSequences with motifs: {len(with_motifs)} ({len(with_motifs)/len(df)*100:.1f}%)")
    print(f"Sequences without motifs: {len(without_motifs)} ({len(without_motifs)/len(df)*100:.1f}%)")
    
    # Statistics by type
    type_stats = pd.DataFrame(stats).groupby('type').agg({
        'coverage_pct': ['mean', 'std', 'min', 'max'],
        'total_motif_hits': ['sum', 'mean']
    }).round(2)
    print("\nMotif coverage by type:")
    print(type_stats)
    
    # Motif count distribution
    print("\nMotif count distribution:")
    motif_dist = df['motif_count'].value_counts().sort_index()
    for count in range(min(6, motif_dist.index.max() + 1)):
        n_seqs = motif_dist.get(count, 0)
        print(f"  {count} motifs: {n_seqs} sequences ({n_seqs/len(df)*100:.1f}%)")
    if motif_dist.index.max() > 5:
        n_seqs = motif_dist[motif_dist.index > 5].sum()
        print(f"  >5 motifs: {n_seqs} sequences ({n_seqs/len(df)*100:.1f}%)")
    
    # Find monomers with highest motif density
    monomer_stats = pd.DataFrame(stats)[pd.DataFrame(stats)['type'] == 'MONOMER']
    if len(monomer_stats) > 0:
        top_monomers = monomer_stats.nlargest(10, 'coverage_pct')
        print("\nTop 10 monomers by motif coverage:")
        print(top_monomers[['sequence_id', 'index', 'length', 'coverage_pct', 'total_motif_hits']])
    
    return df, pd.DataFrame(stats)


def generate_fasta(df, output_file, min_coverage=0):
    """
    Generate FASTA file from marked monomer data.
    
    Args:
        df: DataFrame with marked sequences
        output_file: Output FASTA file
        min_coverage: Minimum coverage to include in FASTA
    """
    print(f"\nGenerating FASTA file: {output_file}")
    
    with open(output_file, 'w') as f:
        for idx, row in df.iterrows():
            # Calculate coverage if needed
            if min_coverage > 0:
                marked_count = sum(1 for c in row['marked_sequence'] if c.islower())
                coverage = (marked_count / len(row['sequence'])) * 100
                if coverage < min_coverage:
                    continue
            
            # Create header
            header = f">{row['sequence_id']}_{row['type']}_{row['index']}"
            header += f" ori={row['orientation']} len={row['length']}"
            
            # Add coverage if available
            if 'marked_sequence' in row:
                marked_count = sum(1 for c in row['marked_sequence'] if c.islower())
                coverage = (marked_count / len(row['sequence'])) * 100
                header += f" motif_coverage={coverage:.1f}%"
            
            # Write sequence
            f.write(header + '\n')
            seq = row.get('marked_sequence', row['sequence'])
            
            # Wrap at 80 characters
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')
    
    print(f"  FASTA file written")


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
        # Load from file
        motifs = load_motifs(args.motif_file, args.top_motifs, args.min_count)
    
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