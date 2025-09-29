#!/usr/bin/env python3
"""
Draw sequence logo from motif_discovery output.
Creates both frequency logo and information content logo.
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib.transforms import Affine2D
from matplotlib.font_manager import FontProperties
import numpy as np
from pathlib import Path


# Classic DNA colors: A=green, C=blue, G=orange, T=red
DNA_COLORS = {
    "A": "#4daf4a",  # green
    "C": "#377eb8",  # blue
    "G": "#ff7f00",  # orange
    "T": "#e41a1c"   # red
}


def draw_frequency_logo(df, output_prefix, title="Motif Frequency Logo"):
    """Draw frequency logo where each column sums to 100%."""
    
    letters = ["A", "C", "G", "T"]
    N = len(df)
    
    # Calculate heights as percentages
    heights_pct = []
    for _, row in df.iterrows():
        # Normalize to 100% per position
        total = sum(float(row[b]) for b in letters)
        if total > 0:
            heights_pct.append({b: (float(row[b]) / total) * 100.0 for b in letters})
        else:
            heights_pct.append({b: 25.0 for b in letters})  # Equal if no data
    
    # Create figure
    fig, ax = plt.subplots(figsize=(max(8, N * 0.5), 5))
    fp = FontProperties(family="DejaVu Sans", weight="bold")
    
    # Draw letters for each position
    for i, hdict in enumerate(heights_pct, start=1):
        # Sort by height to stack properly
        items = sorted(hdict.items(), key=lambda kv: kv[1])
        y = 0.0
        
        for base, h in items:
            if h <= 1e-6:
                continue
                
            # Create text path for the letter
            tp = TextPath((0, 0), base, size=1, prop=fp)
            bb = tp.get_extents()
            tw = max(bb.width, 1e-6)
            th = max(bb.height, 1e-6)
            
            # Scale to fit width and height
            sx = 0.9 / tw  # Width scaling
            sy = h / th     # Height scaling
            
            # Transform and place
            trans = Affine2D().scale(sx, sy).translate(i - 0.5, y)
            patch = PathPatch(tp, transform=trans + ax.transData,
                            facecolor=DNA_COLORS[base], 
                            edgecolor="black", 
                            linewidth=0.5)
            ax.add_patch(patch)
            y += h
    
    # Format plot
    ax.set_xlim(0, N + 1)
    ax.set_ylim(0, 100)
    ax.set_xticks(range(1, N + 1))
    ax.set_xlabel("Position", fontsize=12)
    ax.set_ylabel("Frequency (%)", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.grid(True, axis="y", linestyle="--", linewidth=0.5, alpha=0.3)
    
    plt.tight_layout()
    
    # Save
    png_path = f"{output_prefix}_frequency.png"
    svg_path = f"{output_prefix}_frequency.svg"
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(svg_path, bbox_inches="tight")
    plt.close()
    
    return png_path, svg_path


def draw_information_logo(df, output_prefix, title="Motif Information Content"):
    """Draw information content logo where height represents bits of information."""
    
    letters = ["A", "C", "G", "T"]
    N = len(df)
    
    # Calculate information content heights
    heights_ic = []
    max_ic = 2.0  # Maximum information content for DNA (log2(4))
    
    for _, row in df.iterrows():
        # Get frequencies (already normalized in input)
        freqs = {b: float(row[b]) for b in letters}
        
        # Information content from the file if available
        if "Info_Content" in row:
            ic = float(row["Info_Content"])
        else:
            # Calculate Shannon entropy and information content
            entropy = 0
            for b in letters:
                if freqs[b] > 0:
                    entropy -= freqs[b] * np.log2(freqs[b])
            ic = max_ic - entropy
        
        # Heights proportional to frequency * information content
        heights_ic.append({b: freqs[b] * ic * 50 for b in letters})  # Scale to display
    
    # Create figure
    fig, ax = plt.subplots(figsize=(max(8, N * 0.5), 5))
    fp = FontProperties(family="DejaVu Sans", weight="bold")
    
    # Draw letters for each position
    for i, hdict in enumerate(heights_ic, start=1):
        # Sort by height to stack properly
        items = sorted(hdict.items(), key=lambda kv: kv[1])
        y = 0.0
        
        for base, h in items:
            if h <= 1e-6:
                continue
                
            # Create text path for the letter
            tp = TextPath((0, 0), base, size=1, prop=fp)
            bb = tp.get_extents()
            tw = max(bb.width, 1e-6)
            th = max(bb.height, 1e-6)
            
            # Scale to fit width and height
            sx = 0.9 / tw  # Width scaling
            sy = h / th     # Height scaling
            
            # Transform and place
            trans = Affine2D().scale(sx, sy).translate(i - 0.5, y)
            patch = PathPatch(tp, transform=trans + ax.transData,
                            facecolor=DNA_COLORS[base], 
                            edgecolor="black", 
                            linewidth=0.5)
            ax.add_patch(patch)
            y += h
    
    # Format plot
    ax.set_xlim(0, N + 1)
    ax.set_ylim(0, 100)
    ax.set_xticks(range(1, N + 1))
    ax.set_xlabel("Position", fontsize=12)
    ax.set_ylabel("Information Content (bits)", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    
    # Add y-axis labels for bits
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_yticklabels(["0", "0.5", "1.0", "1.5", "2.0"])
    
    ax.grid(True, axis="y", linestyle="--", linewidth=0.5, alpha=0.3)
    
    # Add horizontal line at 50 (1 bit)
    ax.axhline(y=50, color='gray', linestyle=':', linewidth=0.8, alpha=0.5)
    
    plt.tight_layout()
    
    # Save
    png_path = f"{output_prefix}_ic.png"
    svg_path = f"{output_prefix}_ic.svg"
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(svg_path, bbox_inches="tight")
    plt.close()
    
    return png_path, svg_path


def main():
    parser = argparse.ArgumentParser(
        description="Draw sequence logo from motif discovery results"
    )
    parser.add_argument(
        "input_file",
        help="Logo matrix TSV file (*_logo.tsv from motif_discovery)"
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output prefix for PNG/SVG files (default: input basename)"
    )
    parser.add_argument(
        "-t", "--title",
        default=None,
        help="Title for the logo plot"
    )
    parser.add_argument(
        "--type",
        choices=["frequency", "ic", "both"],
        default="both",
        help="Type of logo to generate (default: both)"
    )
    parser.add_argument(
        "--start",
        type=int,
        default=None,
        help="Start position (1-based) for logo fragment"
    )
    parser.add_argument(
        "--end",
        type=int,
        default=None,
        help="End position (1-based, inclusive) for logo fragment"
    )
    
    args = parser.parse_args()
    
    # Read the logo matrix
    try:
        df = pd.read_csv(args.input_file, sep="\t")
    except Exception as e:
        print(f"Error reading input file: {e}")
        return 1
    
    # Check required columns
    required_cols = ["Position", "A", "C", "G", "T"]
    if not all(col in df.columns for col in required_cols):
        print(f"Error: Input file must have columns: {required_cols}")
        return 1
    
    # Handle start/end positions for fragment
    if args.start or args.end:
        original_len = len(df)
        
        # Convert to 0-based indexing for DataFrame
        start_idx = (args.start - 1) if args.start else 0
        end_idx = args.end if args.end else original_len
        
        # Validate ranges
        if start_idx < 0:
            start_idx = 0
        if end_idx > original_len:
            end_idx = original_len
        if start_idx >= end_idx:
            print(f"Error: Invalid range [{args.start}, {args.end}]")
            return 1
        
        # Slice the dataframe
        df = df.iloc[start_idx:end_idx].reset_index(drop=True)
        
        # Update Position column to be sequential for display
        df["Position"] = range(1, len(df) + 1)
        
        print(f"Using positions {start_idx+1} to {end_idx} (fragment of {original_len} total positions)")
    
    # Set output prefix
    if args.output is None:
        input_path = Path(args.input_file)
        args.output = str(input_path.parent / input_path.stem)
    
    # Set title
    if args.title is None:
        if args.start or args.end:
            start_pos = args.start if args.start else 1
            end_pos = args.end if args.end else original_len
            args.title = f"Motif Logo (positions {start_pos}-{end_pos})"
        else:
            args.title = "Motif Logo"
    
    # Generate logos
    try:
        if args.type in ["frequency", "both"]:
            freq_files = draw_frequency_logo(df, args.output, f"{args.title} (Frequency)")
            print(f"Frequency logo saved: {freq_files[0]} and {freq_files[1]}")
        
        if args.type in ["ic", "both"]:
            ic_files = draw_information_logo(df, args.output, f"{args.title} (Information Content)")
            print(f"Information content logo saved: {ic_files[0]} and {ic_files[1]}")
            
    except Exception as e:
        print(f"Error generating logo: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())