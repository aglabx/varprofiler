#!/usr/bin/env python3
"""
Setup script for VarProfiler package with C++ compilation support.
"""

from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
from setuptools.command.develop import develop
from setuptools.command.build_py import build_py
from pathlib import Path
import subprocess
import sys
import os
import shutil


class CompileCppTools(build_ext):
    """Custom command to compile C++ tools."""
    
    def run(self):
        """Compile C++ tools during installation."""
        cpp_dir = Path(__file__).parent / "varprofiler" / "cpp"
        bin_dir = Path(__file__).parent / "varprofiler" / "bin"
        
        if not cpp_dir.exists():
            print("Warning: C++ source directory not found")
            return
        
        print("Compiling C++ tools...")
        
        # Create bin directory
        bin_dir.mkdir(parents=True, exist_ok=True)
        
        # Also determine system install location
        try:
            import site
            import sysconfig
            # Try to get the scripts directory
            scripts_dir = Path(sysconfig.get_path('scripts'))
            user_scripts_dir = Path(site.USER_BASE) / "bin" if site.USER_BASE else None
        except:
            scripts_dir = None
            user_scripts_dir = None
        
        # Common compilation flags (pthread for multithreading)
        base_flags = ["-std=c++17", "-O3", "-march=native", "-pthread"]
        
        # List of tools to compile
        tools = [
            ("kmer_profiler", base_flags),
            ("cenpb_finder", base_flags),
            ("motif_discovery", base_flags)
        ]
        
        # Detect compiler
        compiler = os.environ.get('CXX', 'g++')
        
        for tool, flags in tools:
            src_file = cpp_dir / f"{tool}.cpp"
            out_file = bin_dir / tool
            
            if not src_file.exists():
                print(f"  Source file {src_file} not found, skipping")
                continue
            
            cmd = [compiler] + flags + ["-o", str(out_file), str(src_file)]
            
            print(f"  Compiling {tool}...")
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                print(f"    ✓ {tool} compiled successfully")
                
                # Make executable
                out_file.chmod(0o755)
                
                # Try to install to system paths
                for install_dir in [scripts_dir, user_scripts_dir]:
                    if install_dir and install_dir.exists():
                        try:
                            dest = install_dir / f"varprofiler_{tool}"
                            shutil.copy2(out_file, dest)
                            dest.chmod(0o755)
                            print(f"    ✓ Installed to {dest}")
                            break
                        except Exception as e:
                            pass
                
            except subprocess.CalledProcessError as e:
                print(f"    ✗ Failed to compile {tool}:")
                print(f"      {e.stderr}")
                print(f"    Note: C++ tools are optional. Python tools will still work.")
            except FileNotFoundError:
                print(f"    ✗ Compiler '{compiler}' not found")
                print(f"    Please install a C++ compiler or set CXX environment variable")
                break
        
        print("C++ compilation finished")
        super().run()


class CustomInstall(install):
    """Custom install command that compiles C++ tools."""
    
    def run(self):
        """Run C++ compilation before standard install."""
        self.run_command('build_ext')
        super().run()


class CustomDevelop(develop):
    """Custom develop command that compiles C++ tools."""
    
    def run(self):
        """Run C++ compilation before development install."""
        self.run_command('build_ext')
        super().run()


# Read the contents of README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding='utf-8')

setup(
    name="varprofiler",
    version="0.1.0",
    author="aglabx",
    author_email="ad3002@gmail.com",
    description="A comprehensive toolkit for analyzing sequence variability, satellite DNA, and centromeres in genomic data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/aglabx/varprofiler",
    project_urls={
        "Bug Tracker": "https://github.com/aglabx/varprofiler/issues",
        "Documentation": "https://github.com/aglabx/varprofiler/wiki",
        "Source Code": "https://github.com/aglabx/varprofiler",
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: C++",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
    ],
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "pandas>=1.3.0",
        "matplotlib>=3.3.0",
        "plotly>=5.0.0",
        "tqdm>=4.60.0",
        "biopython>=1.79",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=3.0.0",
            "black>=22.0.0",
            "flake8>=4.0.0",
            "mypy>=0.900",
        ],
        "docs": [
            "sphinx>=4.0.0",
            "sphinx-rtd-theme>=1.0.0",
        ],
    },
    entry_points={
        "console_scripts": [
            # Main unified CLI
            "varprofiler=varprofiler.cli.main:main",
            
            # C++ tool wrappers (for direct access)
            "varprofiler_kmer_profiler=varprofiler.cli.cpp_wrapper:run_kmer_profiler",
            "varprofiler_cenpb_finder=varprofiler.cli.cpp_wrapper:run_cenpb_finder",
            "varprofiler_motif_discovery=varprofiler.cli.cpp_wrapper:run_motif_discovery",
            
            # Legacy individual tools (for backwards compatibility)
            "varprofiler-detect-satellites=varprofiler.cli.detect_satellites:main",
            "varprofiler-plot=varprofiler.cli.plot_chromosomes:main",
            "varprofiler-analyze-motifs=varprofiler.cli.analyze_motif_satellites:main",
            "varprofiler-visualize-motifs=varprofiler.cli.visualize_motifs_in_satellites:main",
            "varprofiler-visualize-monomers=varprofiler.cli.visualize_motifs_in_monomers:main",
        ],
    },
    package_data={
        "varprofiler": [
            "utils/README.md",
            "cpp/*.cpp",
            "cpp/*.h",
            "cpp/Makefile",
            "bin/*",  # Include compiled binaries if present
        ],
    },
    # Add scripts that will be tracked by pip for proper uninstall
    scripts=[],  # We'll populate this dynamically if binaries exist
    include_package_data=True,
    zip_safe=False,
    cmdclass={
        'build_ext': CompileCppTools,
        'install': CustomInstall,
        'develop': CustomDevelop,
    },
    # Even though we don't have Python extensions, we use this to trigger compilation
    ext_modules=[],
)