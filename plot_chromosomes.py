import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import plotly.graph_objects as go
import argparse
import os
import sys

def read_gff_annotations(gff_file_path, feature_type='satellite'):
    """
    Read annotations from GFF file.
    
    Args:
        gff_file_path (str): Path to GFF file with annotations.
        feature_type (str): Type of features to extract ('satellite', 'centromere', or 'both')
    
    Returns:
        dict: Dictionary with 'satellites' and/or 'centromeres' keys, each containing
              a dict with chromosome names as keys and list of (start, end, attributes) tuples.
    """
    annotations = {'satellites': {}, 'centromeres': {}}
    
    if not gff_file_path or not os.path.exists(gff_file_path):
        return annotations
    
    try:
        with open(gff_file_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                    
                chrom = parts[0]
                feat_type = parts[2]
                start = int(parts[3]) - 1  # Convert to 0-based
                end = int(parts[4])
                attributes = parts[8] if len(parts) > 8 else ''
                
                # Parse attributes
                attr_dict = {}
                for attr in attributes.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attr_dict[key] = value
                
                # Check if this is a centromere region
                if 'centromere' in feat_type.lower() or 'centromere_containing_region' in feat_type:
                    if feature_type in ['centromere', 'both']:
                        if chrom not in annotations['centromeres']:
                            annotations['centromeres'][chrom] = []
                        annotations['centromeres'][chrom].append((start, end, attr_dict))
                
                # Check if this is a satellite/low complexity region
                elif any(term in feat_type.lower() for term in ['satellite', 'tandem', 'low_kmer', 'complex']):
                    # Only include regions longer than 10kb for satellites
                    if end - start > 10000:
                        if feature_type in ['satellite', 'both']:
                            if chrom not in annotations['satellites']:
                                annotations['satellites'][chrom] = []
                            annotations['satellites'][chrom].append((start, end, attr_dict))
        
        # Sort by start position
        for ann_type in annotations:
            for chrom in annotations[ann_type]:
                annotations[ann_type][chrom].sort(key=lambda x: x[0])
        
        # Print summary
        if annotations['satellites']:
            print(f"Loaded satellite annotations for {len(annotations['satellites'])} chromosomes")
        if annotations['centromeres']:
            print(f"Loaded centromere annotations for {len(annotations['centromeres'])} chromosomes")
        
    except Exception as e:
        print(f"Warning: Could not read GFF file: {e}", file=sys.stderr)
    
    return annotations

def create_plotly_chromosome_plot(chrom_df, chrom, satellites, centromeres, output_dir, kmer_size):
    """
    Create an interactive plotly plot for a single chromosome.
    
    Args:
        chrom_df (DataFrame): Data for the chromosome
        chrom (str): Chromosome name
        satellites (dict): Satellite DNA annotations
        centromeres (dict): Centromere annotations
        output_dir (str): Directory for saving plots
        kmer_size (int): Size of k-mers used in analysis
    """
    # Create figure
    fig = go.Figure()
    
    # Add satellite DNA annotations as background shapes if available
    if chrom in satellites:
        for item in satellites[chrom]:
            # Handle both old (start, end) and new (start, end, attrs) formats
            start = item[0]
            end = item[1] if len(item) > 1 else item[0]
            start_mb = start / 1_000_000
            end_mb = end / 1_000_000
            fig.add_vrect(
                x0=start_mb, x1=end_mb,
                fillcolor="orange", opacity=0.15,
                layer="below", line_width=0,
                annotation_text="Satellite DNA" if start == satellites[chrom][0][0] else "",
                annotation_position="top left"
            )
    
    # Add centromere annotations as background shapes if available
    if chrom in centromeres:
        for item in centromeres[chrom]:
            # Handle both old (start, end) and new (start, end, attrs) formats
            start = item[0]
            end = item[1] if len(item) > 1 else item[0]
            start_mb = start / 1_000_000
            end_mb = end / 1_000_000
            fig.add_vrect(
                x0=start_mb, x1=end_mb,
                fillcolor="red", opacity=0.3,
                layer="below", line_width=0,
                annotation_text="Centromere" if start == centromeres[chrom][0][0] else "",
                annotation_position="top right"
            )
    
    # Add main k-mer distribution line and fill
    fig.add_trace(go.Scatter(
        x=chrom_df['position_mb'],
        y=chrom_df['kmer_percentage'],
        mode='lines',
        fill='tonexty' if len(satellites.get(chrom, [])) == 0 else 'tozeroy',
        fillcolor='rgba(0, 139, 139, 0.2)',  # darkcyan with transparency
        line=dict(color='darkcyan', width=2),
        name=f'Unique {kmer_size}-mers',
        hovertemplate='<b>Position:</b> %{x:.1f} Mb<br>' +
                      '<b>Unique k-mers:</b> %{y:.1f}%<br>' +
                      '<extra></extra>'
    ))
    
    # Update layout for better appearance
    fig.update_layout(
        title=dict(
            text=f'Distribution of unique {kmer_size}-mers across chromosome {chrom}',
            font=dict(size=18, family="Arial, sans-serif"),
            x=0.5,
            xanchor='center'
        ),
        xaxis=dict(
            title=dict(text='Genomic position (Mb)', font=dict(size=14)),
            tickfont=dict(size=12),
            showgrid=True,
            gridcolor='lightgray',
            gridwidth=0.5,
            range=[0, chrom_df['position_mb'].max()],
            fixedrange=False  # Allow zooming on X-axis
        ),
        yaxis=dict(
            title=dict(text=f'Unique {kmer_size}-mers (% of window)', font=dict(size=14)),
            tickfont=dict(size=12),
            showgrid=True,
            gridcolor='lightgray',
            gridwidth=0.5,
            range=[0, 110],
            ticksuffix='%',
            fixedrange=True  # Lock Y-axis, no zooming allowed
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        width=1200,
        height=500,
        margin=dict(l=80, r=50, t=80, b=80),
        showlegend=False,
        hovermode='x unified',
        dragmode='zoom'
    )
    
    # Save as HTML for interactivity
    html_filename = os.path.join(output_dir, f'{chrom}_kmer_distribution_interactive.html')
    fig.write_html(html_filename)
    
    # Skip PNG export to avoid Qt issues in headless environments
    # Users can screenshot the HTML version if needed
    
    return html_filename, None

def plot_kmer_distribution(bed_file_path, output_dir, gff_file_path, kmer_size, centromere_gff_path=None):
    """
    Reads a BED file with k-mer counts and creates plots for each chromosome.
    
    Args:
        bed_file_path (str): Path to the input BED file.
        output_dir (str): Directory for saving plots.
        gff_file_path (str or None): Path to GFF file with satellite DNA annotations.
        kmer_size (int): Size of k-mers used in analysis.
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
    satellites = {}
    centromeres = {}
    
    if gff_file_path:
        annotations = read_gff_annotations(gff_file_path, 'satellite')
        satellites = annotations['satellites']
    
    # Load centromere annotations if provided separately
    if centromere_gff_path:
        cent_annotations = read_gff_annotations(centromere_gff_path, 'both')
        if cent_annotations['centromeres']:
            centromeres = cent_annotations['centromeres']
        # Also merge satellites if present
        for chrom in cent_annotations['satellites']:
            if chrom not in satellites:
                satellites[chrom] = cent_annotations['satellites'][chrom]

    # Get list of unique chromosomes from file
    chromosomes = df['chromosome'].unique()
    
    # Calculate window sizes for normalization
    df['window_size'] = df['end'] - df['start']
    
    # Calculate maximum possible k-mers per window and normalize to percentage
    df['max_kmers'] = df['window_size'] - kmer_size + 1
    df['kmer_percentage'] = (df['kmer_count'] / df['max_kmers']) * 100

    for chrom in chromosomes:
        print(f"Creating plots for chromosome: {chrom}...")
        
        # Filter data for current chromosome
        chrom_df = df[df['chromosome'] == chrom].copy()
        
        # Use window midpoint for X axis and convert to megabases (Mb)
        chrom_df['position_mb'] = (chrom_df['start'] + chrom_df['end']) / 2 / 1_000_000
        
        # Add points at the beginning and end of chromosome for complete visualization
        if len(chrom_df) > 0:
            # Add point at position 0 with the same k-mer percentage as first window
            first_row = chrom_df.iloc[0].copy()
            first_row['position_mb'] = 0.0
            
            # Add point at chromosome end with the same k-mer percentage as last window
            last_row = chrom_df.iloc[-1].copy()
            last_row['position_mb'] = chrom_df['end'].max() / 1_000_000
            
            # Combine all points
            chrom_df = pd.concat([pd.DataFrame([first_row]), chrom_df, pd.DataFrame([last_row])], ignore_index=True)

        # Create interactive plotly plot
        html_file, png_file = create_plotly_chromosome_plot(chrom_df, chrom, satellites, centromeres, output_dir, kmer_size)
        print(f"  Interactive plot: {html_file}")
        if png_file:
            print(f"  Static plot: {png_file}")

        # Also create matplotlib plot for compatibility
        try:
            plt.style.use('seaborn-v0_8-whitegrid')
        except OSError:
            try:
                plt.style.use('seaborn-whitegrid')
            except OSError:
                plt.style.use('default')
        fig, ax = plt.subplots(figsize=(18, 7))

        # Add satellite DNA annotations as background shading if available
        if chrom in satellites:
            for start, end, attrs in satellites[chrom]:
                start_mb = start / 1_000_000
                end_mb = end / 1_000_000
                ax.axvspan(start_mb, end_mb, alpha=0.15, color='orange', 
                          label='Satellite DNA' if start == satellites[chrom][0][0] else "")
        
        # Add centromere annotations as vertical bands if available
        if chrom in centromeres:
            for start, end, attrs in centromeres[chrom]:
                start_mb = start / 1_000_000
                end_mb = end / 1_000_000
                ax.axvspan(start_mb, end_mb, alpha=0.3, color='red', 
                          label='Centromere region' if start == centromeres[chrom][0][0] else "")

        # Draw line and filled area beneath it
        ax.plot(chrom_df['position_mb'], chrom_df['kmer_percentage'], 
                label=f'Unique {kmer_size}-mers', color='darkcyan', linewidth=1.5)
        
        ax.fill_between(chrom_df['position_mb'], chrom_df['kmer_percentage'], 
                        alpha=0.2, color='darkcyan')

        # Configure plot appearance for better readability
        ax.set_title(f'Distribution of unique {kmer_size}-mers across chromosome {chrom}', fontsize=18, fontweight='bold', pad=20)
        ax.set_xlabel('Genomic position (Mb)', fontsize=14, labelpad=15)
        ax.set_ylabel(f'Unique {kmer_size}-mers (% of window)', fontsize=14, labelpad=15)
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        
        # Set axis limits for clean appearance
        ax.set_xlim(0, chrom_df['position_mb'].max())
        ax.set_ylim(0, 110)  # Fixed scale 0-110% for breathing room
        
        # Format y-axis as percentage
        ax.get_yaxis().set_major_formatter(
            plt.FuncFormatter(lambda x, p: f'{x:.0f}%'))
        
        # No legend to avoid overlapping with data

        # Save matplotlib plot to file
        output_filename = os.path.join(output_dir, f'{chrom}_kmer_distribution_matplotlib.png')
        plt.savefig(output_filename, dpi=200, bbox_inches='tight')
        plt.close(fig) # Close figure to free memory

    print(f"\nDone! All plots saved in directory '{output_dir}'.")
    print("Interactive HTML plots can be opened in a web browser for zooming and detailed exploration.")
    
    # Create karyotype-wide plot with all chromosomes
    create_karyotype_plot(df, satellites, centromeres, output_dir, kmer_size)
    
    # Create grouped karyotype plot with size-based grouping
    create_grouped_karyotype_plot(df, satellites, centromeres, output_dir, kmer_size)

def create_karyotype_plot(df, satellites, centromeres, output_dir, kmer_size):
    """
    Create a single plot showing all chromosomes with consistent scale.
    1 Mb has the same width across all chromosomes.
    Supports maternal/paternal chromosome pairs displayed side by side.
    """
    print("\nCreating karyotype-wide plot...")
    
    # Group chromosomes by base name and suffix
    def parse_chromosome(chrom):
        """Extract base chromosome name and suffix (mat/pat)"""
        is_maternal = '_mat' in chrom or '_maternal' in chrom
        is_paternal = '_pat' in chrom or '_paternal' in chrom
        
        # Remove suffixes to get base name
        base_chrom = chrom
        for suffix in ['_mat', '_maternal', '_pat', '_paternal']:
            base_chrom = base_chrom.replace(suffix, '')
        
        suffix_type = 'maternal' if is_maternal else ('paternal' if is_paternal else None)
        return base_chrom, suffix_type, chrom
    
    # Parse all chromosomes
    chrom_data = [parse_chromosome(c) for c in df['chromosome'].unique()]
    
    # Group by base chromosome
    from collections import defaultdict
    chrom_groups = defaultdict(list)
    for base, suffix, full_name in chrom_data:
        chrom_groups[base].append((suffix, full_name))
    
    # Sort base chromosomes
    def sort_base_chrom(base):
        if base.startswith('chr'):
            base = base[3:]
        if base == 'X':
            return (1, 23)
        elif base == 'Y':
            return (1, 24)
        elif base == 'Z':
            return (1, 25)
        elif base == 'W':
            return (1, 26)
        elif base == 'M' or base == 'MT':
            return (1, 27)
        else:
            try:
                return (0, int(base))
            except:
                return (2, base)
    
    sorted_bases = sorted(chrom_groups.keys(), key=sort_base_chrom)
    
    # Define sex chromosomes that should not be paired
    sex_chromosomes = {'X', 'Y', 'Z', 'W', 'chrX', 'chrY', 'chrZ', 'chrW'}
    
    # Create list of chromosome pairs (or singles)
    plot_items = []
    for base in sorted_bases:
        chroms = chrom_groups[base]
        # Sort within group: maternal first, then paternal, then no suffix
        chroms.sort(key=lambda x: 0 if x[0] == 'maternal' else (1 if x[0] == 'paternal' else 2))
        
        # Check if this is a sex chromosome
        is_sex_chrom = base in sex_chromosomes or f'chr{base}' in sex_chromosomes
        
        if not is_sex_chrom and len(chroms) == 2 and chroms[0][0] == 'maternal' and chroms[1][0] == 'paternal':
            # We have a pair - add as tuple (but not for sex chromosomes)
            plot_items.append((chroms[0][1], chroms[1][1]))
        else:
            # Add individually (sex chromosomes or unpaired autosomes)
            for _, full_name in chroms:
                plot_items.append((full_name, None))
    
    # Calculate number of rows needed (max 3 chromosome pairs per row)
    n_cols = 3  # Reduced since each pair takes 2 subplot spaces
    n_rows = (len(plot_items) + n_cols - 1) // n_cols
    
    # Create figure with subplots - double the columns for pairs
    fig, axes = plt.subplots(n_rows, n_cols * 2, figsize=(24, 4 * n_rows))
    # Ensure axes is always 2D array for consistent indexing
    if n_rows == 1 and n_cols * 2 == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols * 2 == 1:
        axes = axes.reshape(-1, 1)
    
    # Calculate position_mb for all chromosomes if not already present
    if 'position_mb' not in df.columns:
        df['position_mb'] = (df['start'] + df['end']) / 2 / 1_000_000
    
    # Calculate kmer_percentage if not already present
    if 'kmer_percentage' not in df.columns:
        df['window_size'] = df['end'] - df['start']
        df['max_kmers'] = df['window_size'] - kmer_size + 1
        df['kmer_percentage'] = (df['kmer_count'] / df['max_kmers']) * 100
    
    # Find the maximum chromosome length for consistent scaling
    max_length_mb = 0
    for item in plot_items:
        for chrom in [item[0], item[1]] if item[1] else [item[0]]:
            chrom_df = df[df['chromosome'] == chrom]
            if len(chrom_df) > 0:
                max_length_mb = max(max_length_mb, chrom_df['position_mb'].max())
    
    # Plot each chromosome or pair
    for idx, item in enumerate(plot_items):
        row = idx // n_cols
        base_col = (idx % n_cols) * 2  # Each item gets 2 columns
        
        if item[1] is not None:
            # Plot maternal/paternal pair
            for sub_idx, chrom in enumerate([item[0], item[1]]):
                ax = axes[row, base_col + sub_idx]
        
                chrom_df = df[df['chromosome'] == chrom].copy()
                
                if len(chrom_df) == 0:
                    ax.axis('off')
                    continue
                
                # Add satellite annotations if available
                if chrom in satellites:
                    for item in satellites[chrom]:
                        # Handle both old (start, end) and new (start, end, attrs) formats
                        start = item[0]
                        end = item[1]
                        start_mb = start / 1_000_000
                        end_mb = end / 1_000_000
                        ax.axvspan(start_mb, end_mb, alpha=0.15, color='orange')
                
                # Add centromere annotations if available
                if chrom in centromeres:
                    for item in centromeres[chrom]:
                        start = item[0]
                        end = item[1]
                        start_mb = start / 1_000_000
                        end_mb = end / 1_000_000
                        ax.axvspan(start_mb, end_mb, alpha=0.3, color='red')
                
                # Plot k-mer percentage
                ax.plot(chrom_df['position_mb'], chrom_df['kmer_percentage'], 
                        color='darkcyan', linewidth=0.8)
                ax.fill_between(chrom_df['position_mb'], chrom_df['kmer_percentage'], 
                                alpha=0.2, color='darkcyan')
                
                # Set consistent x-axis scale for all chromosomes
                ax.set_xlim(0, max_length_mb)
                ax.set_ylim(0, 110)  # 0-110% for breathing room
                
                # Format axes
                ax.set_title(chrom, fontsize=10, fontweight='bold')
                ax.set_xlabel('Mb', fontsize=8)
                if base_col == 0:
                    ax.set_ylabel('% unique', fontsize=8)
                ax.tick_params(axis='both', which='major', labelsize=7)
                ax.grid(True, alpha=0.3, linewidth=0.5)
                
                # Format y-axis as percentage
                ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.0f}%'))
        else:
            # Plot single chromosome (no pair)
            ax = axes[row, base_col]
            chrom = item[0]
            
            chrom_df = df[df['chromosome'] == chrom].copy()
            
            if len(chrom_df) == 0:
                ax.axis('off')
                # Also hide the paired column
                axes[row, base_col + 1].axis('off')
                continue
            
            # Add satellite annotations if available
            if chrom in satellites:
                for item in satellites[chrom]:
                    # Handle both old (start, end) and new (start, end, attrs) formats
                    start = item[0]
                    end = item[1]
                    start_mb = start / 1_000_000
                    end_mb = end / 1_000_000
                    ax.axvspan(start_mb, end_mb, alpha=0.15, color='orange')
            
            # Add centromere annotations if available
            if chrom in centromeres:
                for item in centromeres[chrom]:
                    start = item[0]
                    end = item[1]
                    start_mb = start / 1_000_000
                    end_mb = end / 1_000_000
                    ax.axvspan(start_mb, end_mb, alpha=0.3, color='red')
            
            # Plot k-mer percentage
            ax.plot(chrom_df['position_mb'], chrom_df['kmer_percentage'], 
                    color='darkcyan', linewidth=0.8)
            ax.fill_between(chrom_df['position_mb'], chrom_df['kmer_percentage'], 
                            alpha=0.2, color='darkcyan')
            
            # Set consistent x-axis scale for all chromosomes
            ax.set_xlim(0, max_length_mb)
            ax.set_ylim(0, 110)  # 0-110% for breathing room
            
            # Format axes
            ax.set_title(chrom, fontsize=10, fontweight='bold')
            ax.set_xlabel('Mb', fontsize=8)
            if base_col == 0:
                ax.set_ylabel('% unique', fontsize=8)
            ax.tick_params(axis='both', which='major', labelsize=7)
            ax.grid(True, alpha=0.3, linewidth=0.5)
            
            # Format y-axis as percentage
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.0f}%'))
            
            # Hide the paired column
            axes[row, base_col + 1].axis('off')
    
    # Hide remaining empty subplots
    for idx in range(len(plot_items), n_rows * n_cols):
        row = idx // n_cols
        base_col = (idx % n_cols) * 2
        axes[row, base_col].axis('off')
        axes[row, base_col + 1].axis('off')
    
    # Main title
    fig.suptitle(f'Karyotype-wide distribution of unique {kmer_size}-mers', 
                 fontsize=16, fontweight='bold', y=1.02)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save karyotype plot
    output_filename = os.path.join(output_dir, 'karyotype_kmer_distribution.png')
    plt.savefig(output_filename, dpi=200, bbox_inches='tight')
    plt.close(fig)
    
    print(f"Karyotype plot saved as: {output_filename}")

def create_grouped_karyotype_plot(df, satellites, centromeres, output_dir, kmer_size):
    """
    Create a multi-panel plot with chromosomes grouped by size.
    Each group has its own scale to ensure chromosomes fill at least 60% of the plot width.
    """
    print("\nCreating grouped karyotype plot...")
    
    # Calculate position_mb and kmer_percentage if not already present
    if 'position_mb' not in df.columns:
        df['position_mb'] = (df['start'] + df['end']) / 2 / 1_000_000
    
    if 'kmer_percentage' not in df.columns:
        df['window_size'] = df['end'] - df['start']
        df['max_kmers'] = df['window_size'] - kmer_size + 1
        df['kmer_percentage'] = (df['kmer_count'] / df['max_kmers']) * 100
    
    # Get chromosome sizes
    chrom_sizes = {}
    for chrom in df['chromosome'].unique():
        chrom_df = df[df['chromosome'] == chrom]
        if len(chrom_df) > 0:
            chrom_sizes[chrom] = chrom_df['position_mb'].max()
    
    # Define size groups with thresholds (in Mb)
    # We'll create groups to maximize chromosome visibility
    groups = []
    
    # Sort chromosomes by size
    sorted_chroms = sorted(chrom_sizes.items(), key=lambda x: x[1], reverse=True)
    
    # Group chromosomes by size ranges
    # Calculate optimal groupings based on size distribution
    sizes = [size for _, size in sorted_chroms]
    
    if len(sizes) > 0:
        # Define grouping thresholds - always include dot chromosomes as separate group
        if max(sizes) > 200:  # Large genomes
            thresholds = [150, 100, 50, 10, 0]
            group_names = ['Large (>150 Mb)', 'Medium (100-150 Mb)', 'Small (50-100 Mb)', 'Tiny (10-50 Mb)', 'Dot chromosomes (<10 Mb)']
        elif max(sizes) > 100:  # Medium genomes  
            thresholds = [100, 50, 25, 10, 0]
            group_names = ['Large (>100 Mb)', 'Medium (50-100 Mb)', 'Small (25-50 Mb)', 'Tiny (10-25 Mb)', 'Dot chromosomes (<10 Mb)']
        else:  # Small genomes
            thresholds = [50, 25, 10, 5, 0]
            group_names = ['Large (>50 Mb)', 'Medium (25-50 Mb)', 'Small (10-25 Mb)', 'Tiny (5-10 Mb)', 'Dot chromosomes (<5 Mb)']
        
        # Assign chromosomes to groups
        grouped_chroms = {name: [] for name in group_names}
        for chrom, size in sorted_chroms:
            for i, threshold in enumerate(thresholds[:-1]):
                if size > thresholds[i+1]:
                    grouped_chroms[group_names[i]].append((chrom, size))
                    break
        
        # Remove empty groups
        grouped_chroms = {k: v for k, v in grouped_chroms.items() if v}
        
        # Create separate figure for each group for better layout control
        for group_idx, (group_name, chroms) in enumerate(grouped_chroms.items()):
            if not chroms:
                continue
                
            # Calculate number of rows and columns for this group
            n_chroms = len(chroms)
            # More columns for dot chromosomes since they're very small
            if 'Dot' in group_name or '<10 Mb' in group_name or '<5 Mb' in group_name:
                n_cols = min(8, n_chroms)  # Max 8 dot chromosomes per row
            else:
                n_cols = min(6, n_chroms)  # Max 6 regular chromosomes per row
            n_rows = (n_chroms + n_cols - 1) // n_cols
            
            # Create figure for this group
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(24, 4 * n_rows))
            
            # Ensure axes is always 2D array
            if n_rows == 1 and n_cols == 1:
                axes = np.array([[axes]])
            elif n_rows == 1:
                axes = axes.reshape(1, -1)
            elif n_cols == 1:
                axes = axes.reshape(-1, 1)
            
            # Calculate the maximum size in this group for scaling
            group_max_size = max(size for _, size in chroms)
            
            # Plot each chromosome in this group
            for idx, (chrom, size) in enumerate(chroms):
                row = idx // n_cols
                col = idx % n_cols
                ax = axes[row, col]
                
                chrom_df = df[df['chromosome'] == chrom].copy()
                
                if len(chrom_df) == 0:
                    ax.axis('off')
                    continue
                
                # Add satellite annotations if available
                if chrom in satellites:
                    for item in satellites[chrom]:
                        # Handle both old (start, end) and new (start, end, attrs) formats
                        start = item[0]
                        end = item[1]
                        start_mb = start / 1_000_000
                        end_mb = end / 1_000_000
                        ax.axvspan(start_mb, end_mb, alpha=0.15, color='orange')
                
                # Add centromere annotations if available
                if chrom in centromeres:
                    for item in centromeres[chrom]:
                        start = item[0]
                        end = item[1]
                        start_mb = start / 1_000_000
                        end_mb = end / 1_000_000
                        ax.axvspan(start_mb, end_mb, alpha=0.3, color='red')
                
                # Plot k-mer percentage
                ax.plot(chrom_df['position_mb'], chrom_df['kmer_percentage'], 
                        color='darkcyan', linewidth=0.8)
                ax.fill_between(chrom_df['position_mb'], chrom_df['kmer_percentage'], 
                                alpha=0.2, color='darkcyan')
                
                # Set x-axis scale based on group maximum
                # We want chromosomes to fill at least 60% of the plot width
                ax.set_xlim(0, group_max_size * 1.1)  # 10% padding
                ax.set_ylim(0, 110)
                
                # Format axes
                ax.set_title(f'{chrom} ({size:.1f} Mb)', fontsize=10, fontweight='bold')
                ax.set_xlabel('Mb', fontsize=8)
                if col == 0:
                    ax.set_ylabel('% unique', fontsize=8)
                ax.tick_params(axis='both', which='major', labelsize=7)
                ax.grid(True, alpha=0.3, linewidth=0.5)
                
                # Format y-axis as percentage
                ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.0f}%'))
            
            # Hide empty subplots
            for idx in range(len(chroms), n_rows * n_cols):
                row = idx // n_cols
                col = idx % n_cols
                axes[row, col].axis('off')
            
            # Add group title
            fig.suptitle(f'{group_name} - {kmer_size}-mer distribution', 
                        fontsize=14, fontweight='bold')
            
            # Adjust layout
            plt.tight_layout()
            
            # Save grouped plot
            safe_group_name = group_name.replace('>', 'gt').replace('<', 'lt').replace(' ', '_').replace('(', '').replace(')', '')
            output_filename = os.path.join(output_dir, f'grouped_{safe_group_name}_kmer_distribution.png')
            plt.savefig(output_filename, dpi=200, bbox_inches='tight')
            plt.close(fig)
            
            print(f"  Saved {group_name}: {output_filename}")
        
        # Print grouping statistics
        print("\nChromosome grouping statistics:")
        for group_name, chroms in grouped_chroms.items():
            print(f"  {group_name}: {len(chroms)} chromosomes")
            for chrom, size in chroms[:5]:  # Show first 5
                print(f"    - {chrom}: {size:.1f} Mb")
            if len(chroms) > 5:
                print(f"    ... and {len(chroms) - 5} more")

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
    parser.add_argument(
        '-c', '--centromere-gff',
        help='Optional GFF file with centromere annotations (e.g., from detect_satellites.py).'
    )
    parser.add_argument(
        '-k', '--kmer-size',
        type=int,
        required=True,
        help='K-mer size used in the analysis (must match the value used in kmer_profiler).'
    )
    
    args = parser.parse_args()
    plot_kmer_distribution(args.bed_file, args.output_dir, args.gff, args.kmer_size, args.centromere_gff)

if __name__ == '__main__':
    main()
