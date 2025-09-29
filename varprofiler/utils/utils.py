#!/usr/bin/env python3
"""
Common utilities for VarProfiler package.
"""

import sys
from pathlib import Path


def reverse_complement(seq):
    """
    Get reverse complement of DNA sequence.
    
    Args:
        seq: DNA sequence string
        
    Returns:
        Reverse complement string
    """
    complement = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'N': 'N', 'R': 'Y', 'Y': 'R', 'S': 'S',
        'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V',
        'V': 'B', 'D': 'H', 'H': 'D'
    }
    return ''.join(complement.get(c.upper(), c) for c in seq[::-1])


def get_tqdm():
    """
    Get tqdm progress bar or fallback function.
    
    Returns:
        tqdm function or fallback
    """
    try:
        from tqdm import tqdm
        return tqdm
    except ImportError:
        print("Warning: tqdm not installed. Install with: pip install tqdm")
        print("Using simple progress indicator instead.")
        
        def tqdm_fallback(iterable, total=None, desc=None, unit=None):
            if desc:
                print(f"{desc}...")
            return iterable
        
        return tqdm_fallback


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


def write_fasta(sequences, filename, wrap_length=80):
    """
    Write sequences to FASTA file.
    
    Args:
        sequences: Dictionary {id: sequence} or list of (id, sequence) tuples
        filename: Output file path
        wrap_length: Line wrap length (0 for no wrapping)
    """
    with open(filename, 'w') as f:
        # Handle both dict and list of tuples
        if isinstance(sequences, dict):
            sequences = sequences.items()
        
        for seq_id, seq in sequences:
            # Write header
            if not seq_id.startswith('>'):
                f.write(f'>{seq_id}\n')
            else:
                f.write(f'{seq_id}\n')
            
            # Write sequence with wrapping
            if wrap_length > 0:
                for i in range(0, len(seq), wrap_length):
                    f.write(seq[i:i+wrap_length] + '\n')
            else:
                f.write(seq + '\n')


def merge_adjacent_regions(regions, max_gap=50000):
    """
    Merge adjacent regions that are close together.
    
    Args:
        regions: List of (start, end, value) tuples
        max_gap: Maximum gap between regions to merge (default: 50kb)
    
    Returns:
        List of merged (start, end, mean_value) tuples
    """
    import numpy as np
    
    if not regions:
        return []
    
    # Sort by start position
    regions = sorted(regions, key=lambda x: x[0])
    
    merged = []
    current_start = regions[0][0]
    current_end = regions[0][1]
    current_values = [regions[0][2]]
    
    for start, end, value in regions[1:]:
        if start - current_end <= max_gap:
            # Merge with current region
            current_end = max(current_end, end)
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


def find_binary_in_path(binary_name, search_paths=None):
    """
    Find a binary executable in various locations.
    
    Args:
        binary_name: Name of the binary to find
        search_paths: Additional paths to search
    
    Returns:
        Path object to the binary, or None if not found
    """
    import os
    from shutil import which
    
    # Default search paths
    if search_paths is None:
        search_paths = []
    
    # Add common locations
    script_dir = Path(os.path.dirname(os.path.abspath(__file__))).parent
    default_paths = [
        script_dir / binary_name,  # Same directory as package
        Path('./') / binary_name,  # Current directory
        Path(binary_name),        # In PATH
    ]
    
    all_paths = search_paths + default_paths
    
    # Check each path
    for path in all_paths:
        if path.exists() and os.access(path, os.X_OK):
            return path
    
    # Check if it's in system PATH
    in_path = which(binary_name)
    if in_path:
        return Path(in_path)
    
    return None