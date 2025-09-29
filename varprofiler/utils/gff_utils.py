#!/usr/bin/env python3
"""
GFF file handling utilities for VarProfiler.
"""

import sys
from collections import defaultdict


def read_gff_annotations(gff_file, feature_types=None):
    """
    Read annotations from GFF/GFF3 file.
    
    Args:
        gff_file: Path to GFF file
        feature_types: List of feature types to extract, or None for all
    
    Returns:
        Dictionary with chromosome names as keys and list of annotation dicts
    """
    annotations = defaultdict(list)
    
    if not gff_file or not os.path.exists(gff_file):
        return annotations
    
    try:
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                chrom = parts[0]
                source = parts[1]
                feat_type = parts[2]
                start = int(parts[3]) - 1  # Convert to 0-based
                end = int(parts[4])
                score = parts[5]
                strand = parts[6]
                phase = parts[7]
                attributes_str = parts[8]
                
                # Filter by feature type if specified
                if feature_types and feat_type not in feature_types:
                    continue
                
                # Parse attributes
                attributes = parse_gff_attributes(attributes_str)
                
                annotation = {
                    'chrom': chrom,
                    'source': source,
                    'type': feat_type,
                    'start': start,
                    'end': end,
                    'score': score,
                    'strand': strand,
                    'phase': phase,
                    'attributes': attributes
                }
                
                annotations[chrom].append(annotation)
        
        # Sort by start position
        for chrom in annotations:
            annotations[chrom].sort(key=lambda x: x['start'])
    
    except Exception as e:
        print(f"Warning: Could not read GFF file {gff_file}: {e}", file=sys.stderr)
    
    return dict(annotations)


def parse_gff_attributes(attr_string):
    """
    Parse GFF attribute string into dictionary.
    
    Args:
        attr_string: Semicolon-separated attribute string
    
    Returns:
        Dictionary of attributes
    """
    attributes = {}
    for attr in attr_string.split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attributes[key.strip()] = value.strip()
    return attributes


def format_gff_attributes(attr_dict):
    """
    Format attributes dictionary into GFF attribute string.
    
    Args:
        attr_dict: Dictionary of attributes
    
    Returns:
        Semicolon-separated attribute string
    """
    return ';'.join(f"{k}={v}" for k, v in attr_dict.items())


def write_gff_header(file_handle, source='VarProfiler', version=3):
    """
    Write GFF header lines.
    
    Args:
        file_handle: Open file handle
        source: Source program name
        version: GFF version (2 or 3)
    """
    file_handle.write(f"##gff-version {version}\n")
    file_handle.write(f"##source {source}\n")
    if version == 3:
        file_handle.write("##feature-ontology SO\n")


def write_gff_feature(file_handle, chrom, source, feat_type, start, end,
                     score='.', strand='.', phase='.', attributes=None):
    """
    Write a single GFF feature line.
    
    Args:
        file_handle: Open file handle
        chrom: Chromosome/sequence name
        source: Source of the feature
        feat_type: Type of feature
        start: Start position (0-based, will be converted to 1-based)
        end: End position
        score: Score value or '.'
        strand: Strand (+, -, or .)
        phase: Phase for CDS features or '.'
        attributes: Dictionary of attributes or formatted string
    """
    # Convert 0-based to 1-based coordinate
    start_1based = start + 1
    
    # Format attributes if dictionary
    if isinstance(attributes, dict):
        attr_str = format_gff_attributes(attributes)
    elif attributes:
        attr_str = attributes
    else:
        attr_str = '.'
    
    # Write feature line
    file_handle.write(f"{chrom}\t{source}\t{feat_type}\t{start_1based}\t{end}\t"
                     f"{score}\t{strand}\t{phase}\t{attr_str}\n")


def filter_gff_by_type(annotations, feature_types):
    """
    Filter GFF annotations by feature type.
    
    Args:
        annotations: Dictionary from read_gff_annotations
        feature_types: List of feature types to keep
    
    Returns:
        Filtered annotations dictionary
    """
    filtered = {}
    for chrom, features in annotations.items():
        filtered_features = [f for f in features if f['type'] in feature_types]
        if filtered_features:
            filtered[chrom] = filtered_features
    return filtered


def filter_gff_by_size(annotations, min_size=None, max_size=None):
    """
    Filter GFF annotations by feature size.
    
    Args:
        annotations: Dictionary from read_gff_annotations
        min_size: Minimum feature size
        max_size: Maximum feature size
    
    Returns:
        Filtered annotations dictionary
    """
    filtered = {}
    for chrom, features in annotations.items():
        filtered_features = []
        for f in features:
            size = f['end'] - f['start']
            if (min_size is None or size >= min_size) and \
               (max_size is None or size <= max_size):
                filtered_features.append(f)
        if filtered_features:
            filtered[chrom] = filtered_features
    return filtered


def merge_overlapping_features(annotations, max_gap=0):
    """
    Merge overlapping or nearby features.
    
    Args:
        annotations: Dictionary from read_gff_annotations
        max_gap: Maximum gap to bridge when merging
    
    Returns:
        Dictionary with merged features
    """
    merged = {}
    
    for chrom, features in annotations.items():
        if not features:
            continue
        
        # Sort by start position
        sorted_features = sorted(features, key=lambda x: x['start'])
        
        merged_features = []
        current = sorted_features[0].copy()
        
        for feat in sorted_features[1:]:
            if feat['start'] - current['end'] <= max_gap:
                # Merge features
                current['end'] = max(current['end'], feat['end'])
                # Merge attributes if needed
                if 'ID' in current['attributes'] and 'ID' in feat['attributes']:
                    current['attributes']['ID'] += f",{feat['attributes']['ID']}"
            else:
                # Save current and start new
                merged_features.append(current)
                current = feat.copy()
        
        # Save last feature
        merged_features.append(current)
        merged[chrom] = merged_features
    
    return merged


def get_feature_statistics(annotations):
    """
    Calculate statistics for GFF features.
    
    Args:
        annotations: Dictionary from read_gff_annotations
    
    Returns:
        Dictionary with statistics
    """
    import numpy as np
    
    stats = {
        'total_features': 0,
        'chromosomes': set(),
        'feature_types': defaultdict(int),
        'sizes': [],
        'total_bp': 0
    }
    
    for chrom, features in annotations.items():
        stats['chromosomes'].add(chrom)
        for feat in features:
            stats['total_features'] += 1
            stats['feature_types'][feat['type']] += 1
            size = feat['end'] - feat['start']
            stats['sizes'].append(size)
            stats['total_bp'] += size
    
    if stats['sizes']:
        stats['mean_size'] = np.mean(stats['sizes'])
        stats['median_size'] = np.median(stats['sizes'])
        stats['min_size'] = min(stats['sizes'])
        stats['max_size'] = max(stats['sizes'])
    
    return stats


# Import os for path checking
import os