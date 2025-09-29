#!/usr/bin/env python3
"""
Wrapper scripts for C++ binaries to ensure proper installation/uninstallation.
"""

import sys
import subprocess
from pathlib import Path


def get_cpp_binary(name):
    """Find the C++ binary in various locations."""
    package_dir = Path(__file__).parent.parent
    
    search_paths = [
        package_dir / "bin" / name,
        package_dir / "cpp" / name,
        Path(sys.prefix) / "bin" / f"varprofiler_{name}",
        Path.home() / ".local" / "bin" / f"varprofiler_{name}",
    ]
    
    for path in search_paths:
        if path.exists() and path.is_file():
            return str(path)
    
    print(f"Error: {name} binary not found. Please compile C++ tools:")
    print(f"  cd {package_dir}/cpp")
    print(f"  make")
    sys.exit(1)


def run_kmer_profiler():
    """Run kmer_profiler binary."""
    binary = get_cpp_binary("kmer_profiler")
    sys.exit(subprocess.call([binary] + sys.argv[1:]))


def run_cenpb_finder():
    """Run cenpb_finder binary."""
    binary = get_cpp_binary("cenpb_finder")
    sys.exit(subprocess.call([binary] + sys.argv[1:]))


def run_motif_discovery():
    """Run motif_discovery binary."""
    binary = get_cpp_binary("motif_discovery")
    sys.exit(subprocess.call([binary] + sys.argv[1:]))