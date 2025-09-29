#!/usr/bin/env python3
"""
Motif analysis utilities for VarProfiler.
"""

import re
import pandas as pd
from .utils import reverse_complement


def count_motif_in_sequence(motif, sequence, both_strands=True):
    """
    Count occurrences of motif in sequence.
    
    Args:
        motif: Motif sequence to search for
        sequence: DNA sequence to search in
        both_strands: Search both forward and reverse strands
    
    Returns:
        Total count of motif occurrences
    """
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


def find_motif_positions(sequence, motif, both_strands=True):
    """
    Find all positions of a motif in a sequence.
    
    Args:
        sequence: DNA sequence to search in
        motif: Motif sequence to search for
        both_strands: Search both forward and reverse strands
    
    Returns:
        List of (start, end) positions
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


def count_motif_with_regex(motif, sequence, both_strands=True):
    """
    Count motif allowing for IUPAC ambiguity codes.
    
    Args:
        motif: Motif with possible IUPAC codes
        sequence: DNA sequence to search in
        both_strands: Search both strands
    
    Returns:
        Total count of matches
    """
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


def load_motifs_from_tsv(motif_file, top_n=None, min_count=0):
    """
    Load motifs from motif_discovery TSV output.
    
    Args:
        motif_file: TSV file from motif_discovery
        top_n: Load only top N motifs
        min_count: Minimum total count to include motif
    
    Returns:
        DataFrame with motif data
    """
    df = pd.read_csv(motif_file, sep='\t')
    
    # Expected columns from motif_discovery output
    expected_cols = ['motif', 'alignment', 'edit_distance', 'cigar', 
                     'count_forward', 'count_reverse', 'count_total']
    
    # Check if all expected columns are present
    missing_cols = [col for col in expected_cols if col not in df.columns]
    if missing_cols:
        print(f"Warning: Missing columns in motif file: {missing_cols}")
    
    # Filter by minimum count
    if 'count_total' in df.columns and min_count > 0:
        df = df[df['count_total'] >= min_count]
    
    # Get top N if specified
    if top_n:
        df = df.head(top_n)
    
    return df


def load_motifs_from_file(filename, top_n=None, min_count=0):
    """
    Load motifs from either TSV or plain text file.
    
    Args:
        filename: Input file (TSV or text)
        top_n: Load only top N motifs
        min_count: Minimum count for TSV files
    
    Returns:
        List of motif sequences
    """
    try:
        # Try to load as TSV from motif_discovery
        df = load_motifs_from_tsv(filename, top_n, min_count)
        if 'motif' in df.columns:
            return df['motif'].tolist()
    except:
        pass
    
    # Load as plain text file with one motif per line
    with open(filename, 'r') as f:
        motifs = [line.strip() for line in f if line.strip()]
    
    if top_n:
        motifs = motifs[:top_n]
    
    return motifs


def calculate_motif_statistics(sequence, motifs):
    """
    Calculate detailed statistics for motif occurrences in a sequence.
    
    Args:
        sequence: DNA sequence
        motifs: List of motif sequences
    
    Returns:
        Dictionary with statistics
    """
    seq_len = len(sequence)
    
    # Mark motifs
    marked_seq = mark_motifs_in_sequence(sequence, motifs)
    marked_count = sum(1 for c in marked_seq if c.islower())
    coverage = (marked_count / seq_len) * 100 if seq_len > 0 else 0
    
    # Count individual motifs
    motif_counts = {}
    total_hits = 0
    for motif in motifs:
        positions = find_motif_positions(sequence, motif, both_strands=True)
        count = len(positions)
        motif_counts[motif] = count
        total_hits += count
    
    return {
        'sequence_length': seq_len,
        'marked_positions': marked_count,
        'coverage_percent': coverage,
        'total_motif_hits': total_hits,
        'unique_motifs_found': sum(1 for c in motif_counts.values() if c > 0),
        'motif_counts': motif_counts
    }


def find_enriched_motifs(motif_counts_genome, motif_counts_region):
    """
    Find motifs enriched in a specific region compared to genome-wide.
    
    Args:
        motif_counts_genome: Dictionary of genome-wide motif counts
        motif_counts_region: Dictionary of region-specific motif counts
    
    Returns:
        Dictionary with enrichment scores
    """
    enrichment = {}
    
    for motif in motif_counts_region:
        if motif in motif_counts_genome and motif_counts_genome[motif] > 0:
            # Calculate enrichment as ratio
            enrichment[motif] = {
                'region_count': motif_counts_region[motif],
                'genome_count': motif_counts_genome[motif],
                'enrichment_ratio': motif_counts_region[motif] / motif_counts_genome[motif],
                'region_fraction': motif_counts_region[motif] / sum(motif_counts_region.values())
                    if sum(motif_counts_region.values()) > 0 else 0
            }
    
    return enrichment