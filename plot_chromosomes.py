import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse
import os
import sys

def read_gff_satellites(gff_file_path):
    """
    Read satellite DNA annotations from GFF file.
    
    Args:
        gff_file_path (str): Path to GFF file with satellite annotations.
    
    Returns:
        dict: Dictionary with chromosome names as keys and list of (start, end) tuples as values.
    """
    satellites = {}
    
    if not gff_file_path or not os.path.exists(gff_file_path):
        return satellites
    
    try:
        with open(gff_file_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) < 5:
                    continue
                    
                chrom = parts[0]
                start = int(parts[3]) - 1  # Convert to 0-based
                end = int(parts[4])
                
                # Only include satellites longer than 10kb
                if end - start > 10000:
                    if chrom not in satellites:
                        satellites[chrom] = []
                    satellites[chrom].append((start, end))
        
        # Sort satellites by start position
        for chrom in satellites:
            satellites[chrom].sort(key=lambda x: x[0])
            
        print(f"Loaded satellite annotations for {len(satellites)} chromosomes")
        
    except Exception as e:
        print(f"Warning: Could not read GFF file: {e}", file=sys.stderr)
    
    return satellites

def plot_kmer_distribution(bed_file_path, output_dir, gff_file_path=None):
    """
    Reads a BED file with k-mer counts and creates plots for each chromosome.
    
    Args:
        bed_file_path (str): Path to the input BED file.
        output_dir (str): Directory for saving plots.
        gff_file_path (str, optional): Path to GFF file with satellite DNA annotations.
    """
    try:
        # Load data with column names for clarity
        col_names = ['chromosome', 'start', 'end', 'kmer_count']
        df = pd.read_csv(bed_file_path, sep='\t', header=None, names=col_names)
    except FileNotFoundError:
        print(f"Error: File not found at path {bed_file_path}", file=sys.stderr)
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: File {bed_file_path} is empty.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while reading the file: {e}", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(output_dir):
        print(f"Creating directory for plots: {output_dir}")
        os.makedirs(output_dir)

    # Load satellite annotations if provided
    satellites = read_gff_satellites(gff_file_path) if gff_file_path else {}

    # Get list of unique chromosomes from file
    chromosomes = df['chromosome'].unique()

    for chrom in chromosomes:
        print(f"Creating plot for chromosome: {chrom}...")
        
        # Filter data for current chromosome
        chrom_df = df[df['chromosome'] == chrom].copy()
        
        # Use window midpoint for X axis and convert to megabases (Mb)
        chrom_df['position_mb'] = (chrom_df['start'] + chrom_df['end']) / 2 / 1_000_000

        # Create plot
        plt.style.use('seaborn-v0_8-whitegrid')
        fig, ax = plt.subplots(figsize=(18, 7))

        # Add satellite DNA annotations as background shading if available
        if chrom in satellites:
            for start, end in satellites[chrom]:
                start_mb = start / 1_000_000
                end_mb = end / 1_000_000
                ax.axvspan(start_mb, end_mb, alpha=0.15, color='orange', 
                          label='Satellite DNA' if start == satellites[chrom][0][0] else "")

        # Draw line and filled area beneath it
        ax.plot(chrom_df['position_mb'], chrom_df['kmer_count'], 
                label=f'Unique k-mers', color='darkcyan', linewidth=1.5)
        
        ax.fill_between(chrom_df['position_mb'], chrom_df['kmer_count'], 
                        alpha=0.2, color='darkcyan')

        # Configure plot appearance for better readability
        ax.set_title(f'Distribution of unique k-mers across chromosome {chrom}', fontsize=18, fontweight='bold', pad=20)
        ax.set_xlabel('Genomic position (Mb)', fontsize=14, labelpad=15)
        ax.set_ylabel('Number of unique k-mers', fontsize=14, labelpad=15)
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        
        # Set axis limits for clean appearance
        ax.set_xlim(0, chrom_df['position_mb'].max())
        ax.set_ylim(bottom=0) 
        
        # Format axis numbers for large values
        ax.get_yaxis().set_major_formatter(
            plt.FuncFormatter(lambda x, p: format(int(x), ',')))
        
        # Add legend if we have satellite annotations
        if chrom in satellites and len(satellites[chrom]) > 0:
            ax.legend(loc='upper right', fontsize=12, framealpha=0.9)

        # Save plot to file
        output_filename = os.path.join(output_dir, f'{chrom}_kmer_distribution.png')
        plt.savefig(output_filename, dpi=200, bbox_inches='tight')
        plt.close(fig) # Close figure to free memory

    print(f"\nDone! All plots saved in directory '{output_dir}'.")

def main():
    parser = argparse.ArgumentParser(
        description='Visualize unique k-mer counts across genome from BED file.',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('bed_file', help='Path to input BED file generated by kmer_profiler.')
    parser.add_argument(
        '-o', '--output_dir', 
        default='kmer_plots', 
        help='Directory to save plots (default: kmer_plots).'
    )
    parser.add_argument(
        '-g', '--gff',
        help='Optional GFF file with satellite DNA annotations.'
    )
    
    args = parser.parse_args()
    plot_kmer_distribution(args.bed_file, args.output_dir, args.gff)

if __name__ == '__main__':
    main()
