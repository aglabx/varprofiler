#!/usr/bin/env python3
"""
Visualize motif occurrences in satellite sequences by converting them to lowercase.

This script takes satellite sequences and highlights motif occurrences
by converting matched positions to lowercase in the output FASTA.
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


def mark_motifs_in_sequence(sequence, motifs, max_motifs=None):
    """
    Mark motif occurrences in sequence by converting to lowercase.
    
    Args:
        sequence: DNA sequence
        motifs: List of motif sequences or single motif
        max_motifs: If specified, only process top N motifs
    
    Returns:
        Modified sequence with motifs in lowercase
    """
    if isinstance(motifs, str):
        motifs = [motifs]
    
    # Limit number of motifs if specified
    if max_motifs:
        motifs = motifs[:max_motifs]
    
    # Create a position mask for lowercase conversion
    seq_len = len(sequence)
    lowercase_mask = [False] * seq_len
    
    # Find all motif positions
    for motif in motifs:
        positions = find_motif_positions(sequence, motif, both_strands=True)
        for start, end in positions:
            for i in range(start, min(end, seq_len)):
                lowercase_mask[i] = True
    
    # Apply the mask to create the output sequence
    result = []
    for i, char in enumerate(sequence):
        if lowercase_mask[i]:
            result.append(char.lower())
        else:
            result.append(char.upper())
    
    return ''.join(result)


def load_motifs(motif_file, top_n=None, min_count=0):
    """
    Load motifs from motif_discovery output file.
    
    Args:
        motif_file: TSV file from motif_discovery
        top_n: Load only top N motifs
        min_count: Minimum total count to include motif
    
    Returns:
        List of motif sequences
    """
    print(f"Loading motifs from {motif_file}...")
    df = pd.read_csv(motif_file, sep='\t')
    
    # Filter by minimum count
    if min_count > 0:
        df = df[df['count_total'] >= min_count]
    
    # Get top N if specified
    if top_n:
        df = df.head(top_n)
    
    motifs = df['motif'].tolist()
    print(f"  Loaded {len(motifs)} motifs")
    
    return motifs


def process_satellite_file(sat_file, motifs, output_file, 
                          add_stats=False, wrap_length=80):
    """
    Process satellite annotation file and create FASTA with marked motifs.
    
    Args:
        sat_file: Satellite annotation TSV file with trf_array column
        motifs: List of motif sequences to mark
        output_file: Output FASTA file
        add_stats: Add statistics to FASTA headers
        wrap_length: Line length for sequence wrapping (0 = no wrap)
    """
    print(f"Reading satellite data from {sat_file}...")
    df = pd.read_csv(sat_file, sep='\t')
    
    # Check for required column
    if 'trf_array' not in df.columns:
        raise ValueError("Column 'trf_array' not found in satellite file")
    
    # Filter out empty sequences
    df = df[df['trf_array'].notna()]
    df = df[df['trf_array'] != '']
    
    print(f"  Processing {len(df)} satellite arrays...")
    
    # Open output file
    with open(output_file, 'w') as out:
        # Process each satellite
        for idx, row in tqdm(df.iterrows(), total=len(df), desc="Processing satellites"):
            sequence = row['trf_array']
            
            # Skip empty sequences
            if not sequence or pd.isna(sequence):
                continue
            
            # Mark motifs in sequence
            marked_seq = mark_motifs_in_sequence(sequence, motifs)
            
            # Count lowercase (marked) positions
            lowercase_count = sum(1 for c in marked_seq if c.islower())
            coverage = (lowercase_count / len(marked_seq)) * 100 if len(marked_seq) > 0 else 0
            
            # Build header
            header_parts = []
            
            # Add chromosome if available
            if 'trf_chr' in row and pd.notna(row['trf_chr']):
                header_parts.append(f"chr={row['trf_chr']}")
            
            # Add coordinates if available
            if 'trf_l_ind' in row and 'trf_r_ind' in row:
                if pd.notna(row['trf_l_ind']) and pd.notna(row['trf_r_ind']):
                    header_parts.append(f"pos={int(row['trf_l_ind'])}-{int(row['trf_r_ind'])}")
            
            # Add period if available
            if 'trf_period' in row and pd.notna(row['trf_period']):
                header_parts.append(f"period={int(row['trf_period'])}")
            
            # Add type/family if available
            if 'trf_type' in row and pd.notna(row['trf_type']) and row['trf_type'] != '':
                header_parts.append(f"type={row['trf_type']}")
            elif 'trf_family' in row and pd.notna(row['trf_family']) and row['trf_family'] != '':
                header_parts.append(f"family={row['trf_family']}")
            
            # Add statistics
            if add_stats:
                header_parts.append(f"length={len(marked_seq)}")
                header_parts.append(f"motif_bp={lowercase_count}")
                header_parts.append(f"coverage={coverage:.1f}%")
            
            # Create unique ID
            sat_id = f"sat_{idx+1}"
            if 'id' in row and pd.notna(row['id']):
                sat_id = str(row['id'])
            
            # Write header
            header = f">{sat_id} {' '.join(header_parts)}"
            out.write(header + '\n')
            
            # Write sequence (with wrapping if specified)
            if wrap_length > 0:
                for i in range(0, len(marked_seq), wrap_length):
                    out.write(marked_seq[i:i+wrap_length] + '\n')
            else:
                out.write(marked_seq + '\n')
    
    print(f"\nOutput written to {output_file}")
    
    # Print summary statistics
    print("\n=== SUMMARY ===")
    print(f"Total satellites processed: {len(df)}")
    print(f"Total motifs marked: {len(motifs)}")
    print(f"Output format: FASTA with motifs in lowercase")


def process_fasta_file(fasta_file, motifs, output_file, wrap_length=80):
    """
    Process a standard FASTA file and mark motifs.
    
    Args:
        fasta_file: Input FASTA file
        motifs: List of motif sequences to mark
        output_file: Output FASTA file
        wrap_length: Line length for sequence wrapping (0 = no wrap)
    """
    print(f"Reading FASTA from {fasta_file}...")
    
    sequences = []
    current_header = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                # Save previous sequence
                if current_header:
                    sequences.append((current_header, ''.join(current_seq)))
                
                current_header = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last sequence
        if current_header:
            sequences.append((current_header, ''.join(current_seq)))
    
    print(f"  Loaded {len(sequences)} sequences")
    
    # Process and write sequences
    with open(output_file, 'w') as out:
        for header, seq in tqdm(sequences, desc="Processing sequences"):
            # Mark motifs
            marked_seq = mark_motifs_in_sequence(seq, motifs)
            
            # Calculate statistics
            lowercase_count = sum(1 for c in marked_seq if c.islower())
            coverage = (lowercase_count / len(marked_seq)) * 100 if len(marked_seq) > 0 else 0
            
            # Write header with stats
            out.write(f">{header} motif_coverage={coverage:.1f}%\n")
            
            # Write sequence
            if wrap_length > 0:
                for i in range(0, len(marked_seq), wrap_length):
                    out.write(marked_seq[i:i+wrap_length] + '\n')
            else:
                out.write(marked_seq + '\n')
    
    print(f"\nOutput written to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Visualize motif occurrences in sequences by converting to lowercase')
    
    parser.add_argument('input_file',
                       help='Input file: satellite TSV with trf_array column OR FASTA file')
    parser.add_argument('motif_file',
                       help='Motif file: TSV from motif_discovery OR text file with one motif per line')
    parser.add_argument('-o', '--output', required=True,
                       help='Output FASTA file with motifs in lowercase')
    parser.add_argument('--top-motifs', type=int, default=None,
                       help='Use only top N motifs from file')
    parser.add_argument('--min-count', type=int, default=0,
                       help='Minimum motif count to include (for TSV input)')
    parser.add_argument('--single-motif', 
                       help='Mark single motif instead of reading from file')
    parser.add_argument('--stats', action='store_true',
                       help='Add statistics to FASTA headers')
    parser.add_argument('--wrap', type=int, default=80,
                       help='Wrap sequences at specified length (0=no wrap)')
    parser.add_argument('--fasta-input', action='store_true',
                       help='Input is FASTA file instead of satellite TSV')
    
    args = parser.parse_args()
    
    # Load motifs
    if args.single_motif:
        # Use single motif provided on command line
        motifs = [args.single_motif]
        print(f"Using single motif: {args.single_motif}")
    else:
        # Check if motif file is TSV or plain text
        try:
            # Try to load as TSV from motif_discovery
            motifs = load_motifs(args.motif_file, args.top_motifs, args.min_count)
        except:
            # Load as plain text file with one motif per line
            print(f"Loading motifs from text file {args.motif_file}...")
            with open(args.motif_file, 'r') as f:
                motifs = [line.strip() for line in f if line.strip()]
            print(f"  Loaded {len(motifs)} motifs")
    
    # Process input file
    if args.fasta_input:
        # Process as FASTA
        process_fasta_file(args.input_file, motifs, args.output, args.wrap)
    else:
        # Process as satellite TSV
        process_satellite_file(args.input_file, motifs, args.output, 
                             args.stats, args.wrap)
    
    print("\nVisualization complete!")
    print("Motif occurrences are shown in lowercase in the output FASTA.")


if __name__ == '__main__':
    main()