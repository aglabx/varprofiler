#!/usr/bin/env python3
"""
Main entry point for VarProfiler CLI.
Provides unified interface to all tools.
"""

import sys
import os
import argparse
import subprocess
from pathlib import Path


def get_cpp_binary_path(tool_name):
    """Get path to compiled C++ binary."""
    # Check multiple locations
    package_dir = Path(__file__).parent.parent
    
    # Possible locations for binaries
    search_paths = [
        package_dir / "bin" / tool_name,  # Installed location
        package_dir / "cpp" / tool_name,  # Development location
        Path(sys.prefix) / "bin" / f"varprofiler_{tool_name}",  # System install
        Path.home() / ".local" / "bin" / f"varprofiler_{tool_name}",  # User install
    ]
    
    for path in search_paths:
        if path.exists() and os.access(path, os.X_OK):
            return str(path)
    
    return None


def run_cpp_tool(tool_name, args):
    """Run a C++ tool with given arguments."""
    binary_path = get_cpp_binary_path(tool_name)
    
    if not binary_path:
        print(f"Error: {tool_name} binary not found. Please compile C++ tools first:")
        print(f"  cd {Path(__file__).parent.parent}/cpp")
        print(f"  make")
        print(f"  make install")
        return 1
    
    # Run the tool
    cmd = [binary_path] + args
    try:
        result = subprocess.run(cmd, check=False)
        return result.returncode
    except Exception as e:
        print(f"Error running {tool_name}: {e}")
        return 1


def run_profiler(args):
    """Run kmer_profiler C++ tool."""
    # Create parser for kmer_profiler
    parser = argparse.ArgumentParser(
        prog='varprofiler profiler',
        description='Profile k-mer diversity across genomic windows'
    )
    parser.add_argument('fasta', help='Input FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output BED file')
    parser.add_argument('-k', '--kmer-size', type=int, default=21, help='K-mer size (default: 21)')
    parser.add_argument('-w', '--window-size', type=int, default=100000, help='Window size (default: 100000)')
    parser.add_argument('-s', '--step-size', type=int, default=25000, help='Step size (default: 25000)')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads (default: 1)')
    
    parsed_args = parser.parse_args(args)
    
    # Build command line for C++ tool
    cpp_args = [
        parsed_args.fasta,
        parsed_args.output,
        str(parsed_args.kmer_size),
        str(parsed_args.window_size),
        str(parsed_args.step_size),
        str(parsed_args.threads)
    ]
    
    return run_cpp_tool("kmer_profiler", cpp_args)


def run_cenpb(args):
    """Run cenpb_finder C++ tool."""
    # Create parser for cenpb_finder
    parser = argparse.ArgumentParser(
        prog='varprofiler cenpb',
        description='Search for CENP-B boxes in genomic regions'
    )
    parser.add_argument('genome', help='Genome FASTA file')
    parser.add_argument('bed', help='Input BED file with regions')
    parser.add_argument('output', help='Output TSV file')
    parser.add_argument('-p', '--pattern', default='YTTCGTTGGAARCGGGA',
                       help='CENP-B box pattern (default: YTTCGTTGGAARCGGGA)')
    parser.add_argument('-d', '--distance', type=int, default=3,
                       help='Maximum edit distance (default: 3)')
    parser.add_argument('-t', '--threads', type=int, default=1,
                       help='Number of threads (default: 1)')
    
    parsed_args = parser.parse_args(args)
    
    # Build command line for C++ tool
    cpp_args = [
        parsed_args.genome,
        parsed_args.bed,
        parsed_args.output,
        '-p', parsed_args.pattern,
        '-d', str(parsed_args.distance),
        '-t', str(parsed_args.threads)
    ]
    
    return run_cpp_tool("cenpb_finder", cpp_args)


def run_motifs(args):
    """Run motif_discovery C++ tool."""
    # Create parser for motif_discovery
    parser = argparse.ArgumentParser(
        prog='varprofiler motifs',
        description='Discover enriched motifs in sequences'
    )
    parser.add_argument('fasta', help='Input FASTA file')
    parser.add_argument('-m', '--motif', required=True, help='Motif pattern to search')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    parser.add_argument('-e', '--edit-distance', type=int, default=1,
                       help='Maximum edit distance (default: 1)')
    parser.add_argument('-t', '--threads', type=int, default=1,
                       help='Number of threads (default: 1)')
    
    parsed_args = parser.parse_args(args)
    
    # Build command line for C++ tool - positional arguments
    cpp_args = [
        parsed_args.fasta,
        parsed_args.motif,
        str(parsed_args.edit_distance),
        str(parsed_args.threads),
        parsed_args.output
    ]
    
    return run_cpp_tool("motif_discovery", cpp_args)


def run_full_pipeline(args):
    """Run the full VarProfiler pipeline."""
    from .varprofiler_pipeline import main as pipeline_main
    
    # Convert args to the format expected by pipeline
    sys.argv = ['varprofiler'] + args
    return pipeline_main()


def run_detect_satellites(args):
    """Run satellite detection."""
    from .detect_satellites import main as detect_main
    
    sys.argv = ['varprofiler'] + args
    return detect_main()


def run_plot(args):
    """Run plotting tool."""
    from .plot_chromosomes import main as plot_main
    
    sys.argv = ['varprofiler'] + args
    return plot_main()


def run_analyze_motifs(args):
    """Run motif analysis."""
    from .analyze_motif_satellites import main as analyze_main
    
    sys.argv = ['varprofiler'] + args
    return analyze_main()


def run_visualize_motifs(args):
    """Run motif visualization."""
    from .visualize_motifs_in_satellites import main as visualize_main
    
    sys.argv = ['varprofiler'] + args
    return visualize_main()


def run_visualize_monomers(args):
    """Run monomer visualization."""
    from .visualize_motifs_in_monomers import main as monomer_main
    
    sys.argv = ['varprofiler'] + args
    return monomer_main()


def run_draw_logo(args):
    """Run logo drawing tool."""
    from .draw_logo import main as logo_main
    
    sys.argv = ['varprofiler-draw-logo'] + args
    return logo_main()


def main():
    """Main entry point for VarProfiler CLI."""
    
    # Create main parser
    parser = argparse.ArgumentParser(
        prog='varprofiler',
        description='VarProfiler: Comprehensive toolkit for sequence variability analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run full pipeline
  varprofiler full genome.fasta -o results/
  
  # Run individual C++ tools
  varprofiler profiler genome.fasta -o kmers.bed -k 21 -w 100000
  varprofiler cenpb genome.fasta regions.bed output.tsv -p YTTCGTTGGAARCGGGA
  varprofiler motifs genome.fasta -m ATCGATCG -o motifs.tsv -e 2
  
  # Run Python analysis tools
  varprofiler satellites kmers.bed -k 21 -o satellites
  varprofiler plot kmers.tsv -o plots/
  varprofiler analyze-motifs motifs.tsv satellites.tsv -o enrichment.tsv
  varprofiler visualize-motifs satellites.tsv motifs.tsv -o marked.fa
  varprofiler draw-logo motifs_logo.tsv -o motif_logo --start 3 --end 18

For more help on a specific command:
  varprofiler <command> --help
        """
    )
    
    # Add subparsers for different commands
    subparsers = parser.add_subparsers(
        dest='command',
        title='Available commands',
        description='Choose a command to run',
        help='Command help',
        metavar='<command>'
    )
    
    # Full pipeline
    parser_full = subparsers.add_parser(
        'full',
        help='Run complete VarProfiler pipeline',
        description='Run the complete analysis pipeline from genome to visualizations',
        add_help=False  # Will be handled by the actual pipeline parser
    )
    
    # C++ tools
    parser_profiler = subparsers.add_parser(
        'profiler',
        help='Run k-mer profiler (C++ tool)',
        description='Profile k-mer diversity across genomic windows',
        add_help=False  # Pass through to C++ tool
    )
    
    parser_cenpb = subparsers.add_parser(
        'cenpb',
        help='Run CENP-B box finder (C++ tool)',
        description='Search for CENP-B boxes in genomic regions',
        add_help=False  # Pass through to C++ tool
    )
    
    parser_motifs = subparsers.add_parser(
        'motifs',
        help='Run motif discovery (C++ tool)',
        description='Discover enriched motifs in sequences',
        add_help=False  # Pass through to C++ tool
    )
    
    # Python analysis tools
    parser_satellites = subparsers.add_parser(
        'satellites',
        help='Detect satellite DNA regions',
        description='Identify low k-mer diversity regions (satellite DNA)',
        add_help=False
    )
    
    parser_plot = subparsers.add_parser(
        'plot',
        help='Plot k-mer profiles',
        description='Generate visualizations of k-mer diversity',
        add_help=False
    )
    
    parser_analyze = subparsers.add_parser(
        'analyze-motifs',
        help='Analyze motif enrichment',
        description='Compare motif frequencies in satellites vs genome',
        add_help=False
    )
    
    parser_vis_motifs = subparsers.add_parser(
        'visualize-motifs',
        help='Visualize motifs in satellites',
        description='Mark motif occurrences in satellite sequences',
        add_help=False
    )
    
    parser_vis_monomers = subparsers.add_parser(
        'visualize-monomers',
        help='Visualize motifs in monomers',
        description='Mark motif occurrences in monomer sequences',
        add_help=False
    )
    
    parser_logo = subparsers.add_parser(
        'draw-logo',
        help='Draw sequence logo from motif results',
        description='Generate sequence logo visualization from motif discovery output',
        add_help=False
    )
    
    # Parse args
    if len(sys.argv) == 1:
        parser.print_help()
        return 0
    
    # Get command and remaining args
    args, unknown = parser.parse_known_args()
    
    # Route to appropriate handler
    handlers = {
        'full': run_full_pipeline,
        'profiler': run_profiler,
        'cenpb': run_cenpb,
        'motifs': run_motifs,
        'satellites': run_detect_satellites,
        'plot': run_plot,
        'analyze-motifs': run_analyze_motifs,
        'visualize-motifs': run_visualize_motifs,
        'visualize-monomers': run_visualize_monomers,
        'draw-logo': run_draw_logo,
    }
    
    if args.command in handlers:
        return handlers[args.command](unknown)
    else:
        parser.print_help()
        return 1


if __name__ == '__main__':
    sys.exit(main())