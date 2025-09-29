#!/usr/bin/env python3
"""
Visualize motif occurrences in satellite sequences by converting them to lowercase.

This script takes satellite sequences and highlights motif occurrences
by converting matched positions to lowercase in the output FASTA.
"""

import argparse
import pandas as pd
import sys

# Import from our new modules
from varprofiler.utils.utils import get_tqdm, write_fasta
from varprofiler.utils.motif_utils import (
    find_motif_positions,
    mark_motifs_in_sequence,
    load_motifs_from_file
)

# Get tqdm with fallback
tqdm = get_tqdm()


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
    
    # Process sequences and build FASTA data
    fasta_sequences = []
    
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
        
        # Build full header
        header = f"{sat_id} {' '.join(header_parts)}"
        fasta_sequences.append((header, marked_seq))
    
    # Write to FASTA using utility function
    write_fasta(fasta_sequences, output_file, wrap_length=wrap_length)
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
    from varprofiler.utils.utils import read_fasta
    
    print(f"Reading FASTA from {fasta_file}...")
    sequences = read_fasta(fasta_file)
    
    if sequences is None:
        print(f"Error: Could not read {fasta_file}")
        return
    
    print(f"  Loaded {len(sequences)} sequences")
    
    # Process sequences
    marked_sequences = []
    for seq_id, seq in tqdm(sequences.items(), desc="Processing sequences"):
        # Mark motifs
        marked_seq = mark_motifs_in_sequence(seq, motifs)
        
        # Calculate statistics
        lowercase_count = sum(1 for c in marked_seq if c.islower())
        coverage = (lowercase_count / len(marked_seq)) * 100 if len(marked_seq) > 0 else 0
        
        # Add coverage to header
        new_header = f"{seq_id} motif_coverage={coverage:.1f}%"
        marked_sequences.append((new_header, marked_seq))
    
    # Write to file
    write_fasta(marked_sequences, output_file, wrap_length=wrap_length)
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
        # Load from file using utility function
        motifs = load_motifs_from_file(args.motif_file, 
                                       top_n=args.top_motifs, 
                                       min_count=args.min_count)
        print(f"Loaded {len(motifs)} motifs")
    
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