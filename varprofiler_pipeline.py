#!/usr/bin/env python3
"""
VarProfiler Pipeline - Complete analysis from FASTA to visualization

This script runs the complete VarProfiler pipeline:
1. K-mer profiling with kmer_profiler
2. Satellite/centromere detection
3. Visualization of results
"""

import argparse
import os
import sys
import subprocess
import json
from pathlib import Path
import shutil
import time


class VarProfilerPipeline:
    """Main pipeline class for VarProfiler analysis."""
    
    def __init__(self, config):
        """Initialize pipeline with configuration."""
        self.config = config
        self.start_time = time.time()
        
        # Set up output directory structure
        self.setup_output_dirs()
        
        # Log configuration
        self.log_config()
    
    def setup_output_dirs(self):
        """Create output directory structure."""
        base_dir = Path(self.config['output_dir'])
        base_dir.mkdir(parents=True, exist_ok=True)
        
        # Create subdirectories
        self.dirs = {
            'base': base_dir,
            'profiles': base_dir / 'kmer_profiles',
            'satellites': base_dir / 'satellites',
            'plots': base_dir / 'plots',
            'logs': base_dir / 'logs'
        }
        
        for dir_path in self.dirs.values():
            dir_path.mkdir(exist_ok=True)
    
    def log_config(self):
        """Save configuration to JSON file."""
        config_file = self.dirs['logs'] / 'pipeline_config.json'
        with open(config_file, 'w') as f:
            json.dump(self.config, f, indent=2)
        print(f"Configuration saved to: {config_file}")
    
    def log_message(self, message, level='INFO'):
        """Log message with timestamp."""
        timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
        log_entry = f"[{timestamp}] {level}: {message}"
        print(log_entry)
        
        # Also write to log file
        log_file = self.dirs['logs'] / 'pipeline.log'
        with open(log_file, 'a') as f:
            f.write(log_entry + '\n')
    
    def run_command(self, cmd, description, check=True):
        """Run shell command and log output."""
        self.log_message(f"Running: {description}")
        self.log_message(f"Command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=check
            )
            
            # Save stdout and stderr
            if result.stdout:
                stdout_file = self.dirs['logs'] / f"{description.replace(' ', '_')}_stdout.txt"
                with open(stdout_file, 'w') as f:
                    f.write(result.stdout)
            
            if result.stderr:
                stderr_file = self.dirs['logs'] / f"{description.replace(' ', '_')}_stderr.txt"
                with open(stderr_file, 'w') as f:
                    f.write(result.stderr)
            
            self.log_message(f"Completed: {description}")
            return result
            
        except subprocess.CalledProcessError as e:
            self.log_message(f"Error in {description}: {e}", level='ERROR')
            if e.stdout:
                self.log_message(f"stdout: {e.stdout}", level='ERROR')
            if e.stderr:
                self.log_message(f"stderr: {e.stderr}", level='ERROR')
            raise
    
    def step1_kmer_profiling(self):
        """Step 1: Run k-mer profiling."""
        self.log_message("=" * 60)
        self.log_message("STEP 1: K-mer Profiling")
        self.log_message("=" * 60)
        
        # Check if kmer_profiler exists
        kmer_profiler = Path('./kmer_profiler')
        if not kmer_profiler.exists():
            self.log_message("kmer_profiler not found, attempting to compile...", level='WARNING')
            self.compile_kmer_profiler()
        
        # Output file
        bed_file = self.dirs['profiles'] / f"kmer_profile_k{self.config['kmer_size']}.bed"
        
        # Build command
        cmd = [
            str(kmer_profiler),
            self.config['input_fasta'],
            str(bed_file),
            str(self.config['kmer_size']),
            str(self.config['window_size']),
            str(self.config['step_size'])
        ]
        
        # Add threads if specified
        if self.config.get('threads'):
            cmd.append(str(self.config['threads']))
        
        # Run k-mer profiling
        self.run_command(cmd, "k-mer profiling")
        
        self.config['bed_file'] = str(bed_file)
        self.log_message(f"K-mer profile saved to: {bed_file}")
        
        return bed_file
    
    def compile_kmer_profiler(self):
        """Compile k-mer profiler if not present."""
        cpp_file = Path('kmer_profiler.cpp')
        if not cpp_file.exists():
            raise FileNotFoundError("kmer_profiler.cpp not found!")
        
        cmd = ['make']
        self.run_command(cmd, "compiling kmer_profiler")
    
    def step2_detect_satellites(self, bed_file):
        """Step 2: Detect satellite regions and potential centromeres."""
        self.log_message("=" * 60)
        self.log_message("STEP 2: Satellite and Centromere Detection")
        self.log_message("=" * 60)
        
        output_prefix = self.dirs['satellites'] / 'detected_regions'
        
        # Build command
        cmd = [
            'python', 'detect_satellites.py',
            str(bed_file),
            '-k', str(self.config['kmer_size']),
            '-o', str(output_prefix),
            '-p', str(self.config.get('percentile', 5)),
            '-m', str(self.config.get('min_satellite_length', 10000))
        ]
        
        # Add genome file if provided
        if self.config.get('input_fasta'):
            cmd.extend(['-g', self.config['input_fasta']])
        
        # Add TRF annotations if provided
        if self.config.get('trf_annotations'):
            cmd.extend(['-t', self.config['trf_annotations']])
        
        # Add centromere detection flag
        if self.config.get('find_centromeres', True):
            cmd.append('--find-centromeres')
        
        # Add global threshold flag if requested
        if self.config.get('use_global_threshold', False):
            cmd.append('--global-threshold')
        
        # Run satellite detection
        self.run_command(cmd, "satellite detection")
        
        # Save output paths
        self.config['satellite_gff'] = str(output_prefix) + '.gff3'
        self.config['satellite_fasta'] = str(output_prefix) + '.fasta'
        self.config['satellite_summary'] = str(output_prefix) + '_summary.txt'
        
        self.log_message(f"Satellite annotations saved to: {self.config['satellite_gff']}")
        
        return self.config['satellite_gff']
    
    def step3_visualization(self, bed_file, satellite_gff):
        """Step 3: Create visualization plots."""
        self.log_message("=" * 60)
        self.log_message("STEP 3: Visualization")
        self.log_message("=" * 60)
        
        output_dir = self.dirs['plots']
        
        # Build command
        cmd = [
            'python', 'plot_chromosomes.py',
            str(bed_file),
            '-k', str(self.config['kmer_size']),
            '-o', str(output_dir)
        ]
        
        # Add satellite annotations
        if satellite_gff and Path(satellite_gff).exists():
            cmd.extend(['-c', str(satellite_gff)])
        
        # Add TRF annotations if provided separately
        if self.config.get('trf_gff'):
            cmd.extend(['-g', self.config['trf_gff']])
        
        # Run visualization
        self.run_command(cmd, "visualization")
        
        self.log_message(f"Plots saved to: {output_dir}")
        
        return output_dir
    
    def generate_report(self):
        """Generate final HTML report with all results."""
        self.log_message("=" * 60)
        self.log_message("Generating Final Report")
        self.log_message("=" * 60)
        
        report_file = self.dirs['base'] / 'report.html'
        
        html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>VarProfiler Analysis Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        h1 {{ color: #333; }}
        h2 {{ color: #666; border-bottom: 2px solid #ddd; padding-bottom: 5px; }}
        .section {{ margin: 20px 0; }}
        .parameter {{ background: #f5f5f5; padding: 10px; margin: 5px 0; }}
        .file-path {{ font-family: monospace; background: #e8e8e8; padding: 2px 5px; }}
        table {{ border-collapse: collapse; width: 100%; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background: #f5f5f5; }}
        .success {{ color: green; }}
        .warning {{ color: orange; }}
        .error {{ color: red; }}
    </style>
</head>
<body>
    <h1>VarProfiler Analysis Report</h1>
    
    <div class="section">
        <h2>Analysis Parameters</h2>
        <div class="parameter"><b>Input genome:</b> <span class="file-path">{self.config['input_fasta']}</span></div>
        <div class="parameter"><b>K-mer size:</b> {self.config['kmer_size']}</div>
        <div class="parameter"><b>Window size:</b> {self.config['window_size']:,} bp</div>
        <div class="parameter"><b>Step size:</b> {self.config['step_size']:,} bp</div>
        <div class="parameter"><b>Threads:</b> {self.config.get('threads', 'auto')}</div>
        <div class="parameter"><b>Output directory:</b> <span class="file-path">{self.config['output_dir']}</span></div>
    </div>
    
    <div class="section">
        <h2>Output Files</h2>
        <table>
            <tr>
                <th>File Type</th>
                <th>Path</th>
                <th>Description</th>
            </tr>
            <tr>
                <td>K-mer Profile</td>
                <td class="file-path">{self.config.get('bed_file', 'N/A')}</td>
                <td>BED file with k-mer counts per window</td>
            </tr>
            <tr>
                <td>Satellite Annotations</td>
                <td class="file-path">{self.config.get('satellite_gff', 'N/A')}</td>
                <td>GFF3 file with detected low k-mer diversity regions</td>
            </tr>
            <tr>
                <td>Extracted Sequences</td>
                <td class="file-path">{self.config.get('satellite_fasta', 'N/A')}</td>
                <td>FASTA sequences of detected regions</td>
            </tr>
            <tr>
                <td>Summary Statistics</td>
                <td class="file-path">{self.config.get('satellite_summary', 'N/A')}</td>
                <td>Text summary of detected regions</td>
            </tr>
            <tr>
                <td>Visualization Plots</td>
                <td class="file-path">{self.dirs['plots']}</td>
                <td>PNG plots for each chromosome and karyotype views</td>
            </tr>
        </table>
    </div>
    
    <div class="section">
        <h2>Analysis Summary</h2>
"""
        
        # Read and include summary if exists
        if 'satellite_summary' in self.config:
            summary_file = Path(self.config['satellite_summary'])
            if summary_file.exists():
                with open(summary_file, 'r') as f:
                    summary_content = f.read()
                html_content += f"<pre>{summary_content}</pre>"
        
        # Add execution time
        elapsed_time = time.time() - self.start_time
        html_content += f"""
    </div>
    
    <div class="section">
        <h2>Execution Details</h2>
        <div class="parameter"><b>Total execution time:</b> {elapsed_time:.1f} seconds</div>
        <div class="parameter"><b>Pipeline version:</b> 1.0.0</div>
        <div class="parameter"><b>Report generated:</b> {time.strftime('%Y-%m-%d %H:%M:%S')}</div>
    </div>
    
    <div class="section">
        <h2>Visualization Preview</h2>
        <p>Karyotype plot: <span class="file-path">{self.dirs['plots']}/karyotype_kmer_distribution.png</span></p>
        <p>Individual chromosome plots available in: <span class="file-path">{self.dirs['plots']}/</span></p>
    </div>
</body>
</html>"""
        
        with open(report_file, 'w') as f:
            f.write(html_content)
        
        self.log_message(f"HTML report saved to: {report_file}")
        return report_file
    
    def run(self):
        """Run the complete pipeline."""
        try:
            # Step 1: K-mer profiling
            bed_file = self.step1_kmer_profiling()
            
            # Step 2: Satellite detection (optional)
            satellite_gff = None
            if self.config.get('detect_satellites', True):
                satellite_gff = self.step2_detect_satellites(bed_file)
            
            # Step 3: Visualization
            self.step3_visualization(bed_file, satellite_gff)
            
            # Generate report
            report_file = self.generate_report()
            
            # Final summary
            elapsed_time = time.time() - self.start_time
            self.log_message("=" * 60)
            self.log_message(f"Pipeline completed successfully in {elapsed_time:.1f} seconds")
            self.log_message(f"Results saved to: {self.config['output_dir']}")
            self.log_message(f"Report available at: {report_file}")
            self.log_message("=" * 60)
            
            return True
            
        except Exception as e:
            self.log_message(f"Pipeline failed: {e}", level='ERROR')
            return False


def parse_config_file(config_file):
    """Parse JSON configuration file."""
    with open(config_file, 'r') as f:
        return json.load(f)


def main():
    parser = argparse.ArgumentParser(
        description='VarProfiler Pipeline - Complete analysis from FASTA to visualization',
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # Input/output arguments
    parser.add_argument('input_fasta', 
                       help='Input genome FASTA file')
    parser.add_argument('-o', '--output-dir',
                       default='varprofiler_results',
                       help='Output directory for all results (default: varprofiler_results)')
    
    # K-mer profiling parameters
    parser.add_argument('-k', '--kmer-size',
                       type=int, default=23,
                       help='K-mer size (default: 23)')
    parser.add_argument('-w', '--window-size',
                       type=int, default=100000,
                       help='Window size in bp (default: 100000)')
    parser.add_argument('-s', '--step-size',
                       type=int, default=25000,
                       help='Step size in bp (default: 25000)')
    parser.add_argument('-t', '--threads',
                       type=int,
                       help='Number of threads (default: auto-detect)')
    
    # Satellite detection parameters
    parser.add_argument('--no-satellites',
                       action='store_true',
                       help='Skip satellite detection step')
    parser.add_argument('-p', '--percentile',
                       type=float, default=5,
                       help='Percentile threshold for low variability (default: 5)')
    parser.add_argument('-m', '--min-satellite-length',
                       type=int, default=10000,
                       help='Minimum satellite region length in bp (default: 10000)')
    parser.add_argument('--no-centromeres',
                       action='store_true',
                       help='Skip centromere detection')
    parser.add_argument('--global-threshold',
                       action='store_true',
                       help='Use global threshold for satellite detection')
    
    # Additional inputs
    parser.add_argument('--trf-annotations',
                       help='GFF file with TRF annotations')
    parser.add_argument('--trf-gff',
                       help='Additional GFF file for visualization overlay')
    
    # Configuration file
    parser.add_argument('-c', '--config',
                       help='JSON configuration file (overrides command line arguments)')
    
    # Quick profiles
    parser.add_argument('--profile',
                       choices=['quick', 'standard', 'detailed'],
                       help='Use predefined parameter profile:\n'
                            '  quick: k=15, window=200kb, step=100kb\n'
                            '  standard: k=23, window=100kb, step=25kb (default)\n'
                            '  detailed: k=31, window=50kb, step=10kb')
    
    args = parser.parse_args()
    
    # Build configuration dictionary
    config = {
        'input_fasta': args.input_fasta,
        'output_dir': args.output_dir,
        'kmer_size': args.kmer_size,
        'window_size': args.window_size,
        'step_size': args.step_size,
        'detect_satellites': not args.no_satellites,
        'find_centromeres': not args.no_centromeres,
        'percentile': args.percentile,
        'min_satellite_length': args.min_satellite_length,
        'use_global_threshold': args.global_threshold
    }
    
    # Add optional parameters
    if args.threads:
        config['threads'] = args.threads
    if args.trf_annotations:
        config['trf_annotations'] = args.trf_annotations
    if args.trf_gff:
        config['trf_gff'] = args.trf_gff
    
    # Apply profile if specified
    if args.profile:
        profiles = {
            'quick': {'kmer_size': 15, 'window_size': 200000, 'step_size': 100000},
            'standard': {'kmer_size': 23, 'window_size': 100000, 'step_size': 25000},
            'detailed': {'kmer_size': 31, 'window_size': 50000, 'step_size': 10000}
        }
        config.update(profiles[args.profile])
        print(f"Using {args.profile} profile")
    
    # Override with config file if provided
    if args.config:
        file_config = parse_config_file(args.config)
        config.update(file_config)
        print(f"Loaded configuration from {args.config}")
    
    # Check input file exists
    if not Path(config['input_fasta']).exists():
        print(f"Error: Input file {config['input_fasta']} not found!")
        sys.exit(1)
    
    # Run pipeline
    pipeline = VarProfilerPipeline(config)
    success = pipeline.run()
    
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()