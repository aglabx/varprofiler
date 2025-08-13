# VarProfiler

A high-performance tool for analyzing k-mer variability across genomes, designed to identify regions of unique sequence complexity. VarProfiler helps researchers understand genomic variation patterns by computing unique k-mer counts in sliding windows across chromosomes.

## Features

### Core Tools
- **kmer_profiler**: Fast C++ tool for k-mer counting in sliding windows
  - Efficient 2-bit encoding for k-mers up to 31bp
  - Multi-threaded chromosome processing
  - Canonical k-mer handling (forward/reverse complement)
  
- **cenpb_finder**: Specialized C++ tool for CENP-B box detection
  - Edit distance search with Myers' algorithm
  - IUPAC degenerate nucleotide support
  - Both strand searching
  - Unlimited thread scaling

- **motif_discovery**: Genome-wide motif variant discovery tool
  - Finds all variants within edit distance of reference motif
  - Reports frequency statistics per strand
  - Supports IUPAC degenerate codes
  - Multi-threaded genome scanning

### Analysis Capabilities
- **Satellite DNA detection**: Automatic identification of low k-mer diversity regions
- **Advanced centromere prediction**: Multi-signal scoring system combining:
  - CENP-B box density (primary signal, 40% weight)
  - K-mer diversity minima (30% weight)
  - Satellite array size context (20% weight)
  - Local minimum detection (10% weight)
- **CENP-B box quantification**: Count and map CENP-B binding sites with biological calibration
- **Confidence scoring**: High/Medium/Low confidence levels for centromere predictions
- **Sequence extraction**: Export detected regions as FASTA sequences
- **TRF integration**: Classify repeats using Tandem Repeat Finder annotations

### Visualization & Reporting
- **Chromosome plots**: Individual plots for each chromosome
- **Karyotype view**: All chromosomes with consistent scaling
- **Grouped plots**: Size-based grouping for better visibility
- **HTML reports**: Complete analysis summary with all results
- **GFF annotations**: Standard format for genome browsers

## Installation

### Prerequisites

- C++ compiler with C++17 support (g++ or clang++)
- Python 3.7+
- Required Python packages: pandas, matplotlib

### Building from source

```bash
# Clone the repository
git clone https://github.com/aglabx/varprofiler.git
cd varprofiler

# Compile the C++ tools (kmer_profiler and cenpb_finder)
make

# Install Python dependencies
pip install -r requirements.txt

# Or install as a package
pip install .
```

## Quick Start - Complete Pipeline

The easiest way to run VarProfiler is using the unified pipeline script that handles all analysis steps automatically:

```bash
# Basic usage - runs complete analysis with default parameters
python varprofiler_pipeline.py genome.fa -o results/

# Use predefined profiles for different analysis depths
python varprofiler_pipeline.py genome.fa -o results/ --profile standard

# Use configuration file for reproducible analyses
python varprofiler_pipeline.py genome.fa -c config.json

# Using Makefile shortcuts
make standard-run GENOME=genome.fa
```

### Pipeline Profiles

VarProfiler includes three predefined analysis profiles:

| Profile | K-mer Size | Window Size | Step Size | Use Case |
|---------|------------|-------------|-----------|----------|
| quick | 15 | 200kb | 100kb | Fast overview, large genomes |
| standard | 23 | 100kb | 25kb | Default, balanced analysis |
| detailed | 31 | 50kb | 10kb | High resolution, small genomes |

### Pipeline Output

The pipeline creates an organized output directory with:
```
varprofiler_results/
├── kmer_profiles/       # BED files with k-mer counts
├── satellites/          # Detected regions (GFF, FASTA, summary)
├── plots/              # Chromosome and karyotype visualizations
├── logs/               # Pipeline logs and configuration
└── report.html         # Complete HTML report with all results
```

### Configuration File

For reproducible analyses, use a JSON configuration file (see `config_example.json`):

```json
{
  "input_fasta": "genome.fa",
  "output_dir": "varprofiler_results",
  "kmer_size": 23,
  "window_size": 100000,
  "step_size": 25000,
  "threads": 8,
  "detect_satellites": true,
  "find_centromeres": true,
  "percentile": 5,
  "min_satellite_length": 10000,
  "trf_annotations": "trf_output.gff"
}
```

Run with config: `python varprofiler_pipeline.py genome.fa -c config.json`

## Manual Usage - Step by Step

For more control over individual steps, you can run each component separately:

### Step 1: Generate k-mer profile

```bash
./kmer_profiler <input.fasta> <output.bed> <kmer_len> <window_size> <step_size> [threads]
```

**Parameters:**
- `input.fasta`: Input genome file in FASTA format
- `output.bed`: Output BED file with k-mer counts
- `kmer_len`: Length of k-mers (1-31, recommended: 23)
- `window_size`: Size of sliding window in bp (e.g., 100000)
- `step_size`: Step size for sliding window in bp (e.g., 25000)
- `threads`: (Optional) Number of threads to use (default: number of CPU cores)

**Examples:**
```bash
# Using default number of threads (auto-detect CPU cores)
./kmer_profiler genome.fa kmer_counts.bed 23 100000 25000

# Using specific number of threads (e.g., 8 threads)
./kmer_profiler genome.fa kmer_counts.bed 23 100000 25000 8
```

### Step 2: Detect satellite DNA regions (optional)

```bash
python detect_satellites.py <input.bed> -k <kmer_size> [-o output_prefix] [-g genome.fa]
```

**Parameters:**
- `input.bed`: BED file generated by kmer_profiler
- `-k, --kmer-size`: **Required** - K-mer size used in analysis
- `-o, --output-prefix`: Prefix for output files (default: satellites)
- `-p, --percentile`: Percentile threshold for low variability (default: 5)
- `-m, --min-length`: Minimum region length in bp (default: 10000)
- `-g, --genome`: Genome FASTA file for sequence extraction
- `-t, --trf-annotations`: GFF file with TRF annotations for repeat classification
- `--global-threshold`: Use global threshold instead of per-chromosome
- `--find-centromeres`: Attempt to identify centromere candidates
- `--all-centromeres`: Report all centromere candidates, not just best per chromosome
- `--centromere-min-score`: Minimum score for centromere candidates (0-100, default: 30)

**Output files:**
- `.gff3`: GFF3 annotations for all detected regions
- `.fasta`: FASTA sequences of detected regions (if genome provided)
- `_centromeres.fasta`: FASTA sequences of predicted centromeres with detailed headers
- `_summary.txt`: Summary statistics
- `_cenpb.tsv`: CENP-B box analysis results (if genome provided)

**Example:**
```bash
# Basic satellite detection
python detect_satellites.py kmer_counts.bed -k 23 -o human_satellites

# With sequence extraction and centromere detection
python detect_satellites.py kmer_counts.bed -k 23 -o human_satellites -g genome.fa --find-centromeres

# With TRF annotations for repeat classification
python detect_satellites.py kmer_counts.bed -k 23 -o human_satellites -t trf_annotations.gff

# Full analysis with all features including CENP-B box detection
python detect_satellites.py kmer_counts.bed -k 23 -o human_satellites -g genome.fa -t trf_annotations.gff --find-centromeres --cenpb-threads 4

# Find ALL potential centromeres (not just best per chromosome)
python detect_satellites.py kmer_counts.bed -k 23 -o all_centromeres -g genome.fa --find-centromeres --all-centromeres --centromere-min-score 40
```

### Step 3: Visualize results

```bash
python plot_chromosomes.py <input.bed> -k <kmer_size> [-o output_dir] [-g gff_file]
```

**Parameters:**
- `input.bed`: BED file generated by kmer_profiler
- `-k, --kmer-size`: **Required** - K-mer size used in analysis (must match kmer_profiler)
- `-o, --output_dir`: Directory to save plots (default: kmer_plots)
- `-g, --gff`: Optional GFF file with satellite DNA annotations
- `-c, --centromere-gff`: Optional GFF file with centromere annotations (e.g., from detect_satellites.py)

**Examples:**
```bash
# Basic visualization (k=23 was used in kmer_profiler)
python plot_chromosomes.py kmer_counts.bed -k 23 -o chromosome_plots

# With satellite DNA annotations overlay (k=31 was used)
python plot_chromosomes.py kmer_counts.bed -k 31 -o chromosome_plots -g satellites.gff

# With centromere regions from detect_satellites.py
python plot_chromosomes.py kmer_counts.bed -k 23 -o chromosome_plots -c satellites.gff3

# Combined: both TRF satellites and detected centromeres
python plot_chromosomes.py kmer_counts.bed -k 23 -o chromosome_plots -g trf_satellites.gff -c detected_regions.gff3

# For diploid genomes with maternal/paternal chromosomes
# Input BED file should contain chromosomes named like: chr1_mat, chr1_pat, chr2_maternal, chr2_paternal, etc.
python plot_chromosomes.py diploid_kmers.bed -k 23 -o diploid_plots
```

## Output Format

### BED file format
The tool outputs a standard BED file with the following columns:
1. Chromosome name
2. Window start position (0-based)
3. Window end position
4. Count of unique k-mers in the window

### Visualization
- Generates one PNG file per chromosome
- Additional karyotype-wide plot showing all chromosomes together
- High-resolution plots (200 DPI) suitable for publication
- X-axis: Genomic position in megabases (Mb)
- Y-axis: Percentage of unique k-mers in window (0-110%)
- Fixed Y-axis scale for easy comparison across chromosomes
- Orange shading: Satellite DNA regions from GFF annotations (when provided)
- Red shading: Centromere-containing regions (when provided via -c)
- Cyan line/area: K-mer variability profile
- Title includes k-mer size for clarity

#### Karyotype Plots

**Standard Karyotype Plot**
- Single image with all chromosomes arranged in grid (3 pairs per row for diploid genomes)
- Consistent scale: 1 Mb has same width across all chromosomes
- Natural chromosome sorting (1, 2, ..., 10, ..., 22, X, Y, Z, W, MT)
- **Diploid genome support**: 
  - For diploid assemblies, use `_mat` or `_maternal` suffix for maternal chromosomes and `_pat` or `_paternal` suffix for paternal chromosomes
  - Autosomal chromosomes will be displayed side by side as pairs
  - Sex chromosomes (X, Y, Z, W) are always displayed individually, even if they have mat/pat suffixes
- Saved as `karyotype_kmer_distribution.png`

**Grouped Karyotype Plots** 
- Chromosomes grouped by size for better visibility of smaller chromosomes
- Each size group saved as separate file with its own scale
- Chromosomes fill at least 60% of plot width in each group
- Automatic grouping based on genome size:
  - Large genomes: >150 Mb, 100-150 Mb, 50-100 Mb, 10-50 Mb, Dot chromosomes (<10 Mb)
  - Medium genomes: >100 Mb, 50-100 Mb, 25-50 Mb, 10-25 Mb, Dot chromosomes (<10 Mb)  
  - Small genomes: >50 Mb, 25-50 Mb, 10-25 Mb, 5-10 Mb, Dot chromosomes (<5 Mb)
- Dot chromosomes displayed with up to 8 per row for better space utilization
- Chromosome sizes shown in titles for easy reference
- Saved as separate files: `grouped_Large_gt150_Mb_kmer_distribution.png`, `grouped_Dot_chromosomes_lt10_Mb_kmer_distribution.png`, etc.

## Centromere Detection Strategy

VarProfiler uses a sophisticated multi-signal approach for centromere identification, moving beyond simple k-mer minima to a biologically-informed scoring system.

### Scoring System (0-100 scale)

1. **CENP-B Box Density (40% weight)** - Primary signal
   - >2.5 boxes/kb: Maximum score (alpha-satellite rich)
   - 1-2.5 boxes/kb: Moderate score (typical arrays)
   - <1 box/kb: Low score (degenerate/poor arrays)
   - Expected: ~3 boxes/kb in human alpha-satellites

2. **K-mer Diversity (30% weight)** - Conservation signal
   - <2%: Maximum score (highly conserved)
   - 2-10%: Moderate score
   - >10%: Low score

3. **Satellite Array Size (20% weight)** - Context signal
   - >5 Mb: Maximum score (large arrays)
   - 1-5 Mb: Moderate score
   - <1 Mb: Low score

4. **Local Minimum (10% weight)** - Additional confidence
   - Bonus if window is local k-mer minimum within array

### Confidence Levels
- **High (≥70)**: Red - Strong centromeric signature
- **Medium (50-69)**: Orange - Probable centromere
- **Low (<50)**: Yellow - Possible centromere

## CENP-B Box Detection

The `cenpb_finder` tool searches for CENP-B boxes - conserved ~17bp sequences that mark functional centromeres. The canonical human CENP-B box is `YTTCGTTGGAARCGGGA`.

### Step 2b: Search for CENP-B boxes (optional)

```bash
./cenpb_finder <genome.fasta> <regions.bed> <output.tsv> [options]
```

**Parameters:**
- `genome.fasta`: Input genome in FASTA format
- `regions.bed`: BED file with regions to search (e.g., from detect_satellites.py)
- `output.tsv`: Output file with CENP-B box counts and positions
- `-p PATTERN`: CENP-B box pattern with IUPAC codes (default: YTTCGTTGGAARCGGGA)
- `-d DISTANCE`: Maximum edit distance (default: 3)
- `-t THREADS`: Number of threads (default: 1, no upper limit)
- `-s`: Search forward strand only (default: both strands)
- `-v`: Verbose output

**Examples:**
```bash
# Direct usage - search in low k-mer diversity regions
./cenpb_finder genome.fa low_kmer_regions.bed cenpb_boxes.tsv -d 3 -t 192

# Custom pattern for different species
./cenpb_finder mouse.fa satellites.bed mouse_cenpb.tsv -p YTTCGTTGGAWRCGGGA -d 4

# Integrated with satellite detection (automatic)
python detect_satellites.py kmer_counts.bed -k 23 -g genome.fa --find-centromeres \
  --cenpb-pattern YTTCGTTGGAARCGGGA --cenpb-distance 3 --cenpb-threads 8
```

**Output format (TSV):**
```
#chrom  start    end      kmer_count  cenpb_count  cenpb_density  cenpb_positions
chr1    5000000  5100000  1234        15           0.15           5001234,5002345,5003456,...
chr1    5200000  5300000  1456        23           0.23           5201234,5202345,...
```

### Features
- **Edit distance search**: Uses Myers' bit-parallel algorithm for efficient approximate matching
- **IUPAC support**: Handles degenerate nucleotide codes (Y=C/T, R=A/G, W=A/T, etc.)
- **Both strand search**: Searches forward and reverse complement
- **Parallel processing**: Multi-threaded with no artificial thread limit
- **Memory efficient**: Processes regions in parallel without loading entire genome
- **Integration**: Automatically runs on low k-mer diversity regions when detecting centromeres

### Biological significance
CENP-B boxes are binding sites for CENP-B protein, which:
- Helps establish and maintain centromeric chromatin
- Is found in most human centromeres (except Y chromosome)
- Shows conservation across primates and mammals
- High density indicates functional centromeric regions

The CENP-B box count and density are included in the GFF annotations and help confirm functional centromeres.

### Step 2c: Discover CENP-B box variants (optional)

For non-human genomes or to refine the CENP-B box consensus, use the motif discovery tool:

```bash
./motif_discovery <genome.fasta> <motif> <max_distance> <threads> <output.csv>
```

**Parameters:**
- `genome.fasta`: Input genome in FASTA format
- `motif`: Reference motif with IUPAC codes (e.g., YTTCGTTGGAARCGGGA)
- `max_distance`: Maximum edit distance (0-10)
- `threads`: Number of parallel threads
- `output.csv`: Output CSV with variant statistics

**Output CSV columns:**
- `motif`: Discovered sequence variant
- `edit_distance`: Distance from reference motif
- `count_forward`: Occurrences on forward strand
- `count_reverse`: Occurrences on reverse strand
- `count_total`: Total occurrences
- `freq_forward`: Percentage on forward strand
- `freq_reverse`: Percentage on reverse strand
- `freq_total`: Overall percentage

**Examples:**
```bash
# Discover CENP-B variants in human genome
./motif_discovery human.fa YTTCGTTGGAARCGGGA 3 8 cenpb_variants.csv

# Search for variants in mouse genome with relaxed threshold
./motif_discovery mouse.fa YTTCGTTGGAWRCGGGA 4 16 mouse_cenpb.csv

# Find exact matches only
./motif_discovery genome.fa TTCGTTGGAAACGGGA 0 4 exact_matches.csv
```

**Use cases:**
1. **Species-specific optimization**: Find the most common CENP-B variant in your species
2. **Motif refinement**: Identify degenerate positions in the consensus
3. **Quality control**: Verify expected motif frequencies
4. **Comparative genomics**: Compare motif usage across species

The tool will report the top 10 most frequent variants to the console and save all variants sorted by frequency to the CSV file.

## Algorithm

VarProfiler uses an efficient sliding window approach:

1. **K-mer encoding**: Nucleotides are encoded in 2-bit format (A=00, C=01, G=10, T=11)
2. **Canonical k-mers**: For each k-mer, the lexicographically smaller of the forward and reverse complement is used
3. **Gap handling**: K-mers containing 'N' bases are skipped, preventing assembly gaps from being misidentified as low-diversity regions
4. **Sliding window**: 
   - Initial window: Count all unique k-mers
   - Subsequent windows: Add new k-mers entering the window, remove k-mers leaving
5. **Parallel processing**: Each chromosome is processed in a separate thread

## Performance Considerations

- **Memory usage**: Approximately 1GB per 100Mb of sequence data
- **Speed**: Processes human genome (~3Gb) in approximately 30-60 minutes on a modern multi-core system
- **Optimization tips**:
  - Use larger step sizes for faster processing
  - Adjust window size based on desired resolution
  - K-mer length of 23 provides good balance of specificity and speed

## Applications

- **Genome assembly validation**: Identify regions of low complexity or repetitive sequences
- **Comparative genomics**: Compare k-mer profiles between species or individuals
- **Variant detection**: Identify regions with high sequence variability
- **Primer design**: Find regions with unique k-mers for specific amplification
- **Genome annotation**: Correlate k-mer complexity with functional elements

## Citation

If you use VarProfiler in your research, please cite:
```
Marina Popova, Aleksey Komissarov (2025)
VarProfiler: A tool for k-mer variability profiling across genomes
bioRxiv (manuscript in preparation)
```

## Authors

- Marina Popova
- Aleksey Komissarov
- Claude (AI assistant) - Co-author for code implementation and documentation

## License

MIT License - See LICENSE file for details

## Contributing

Contributions are welcome! Please feel free to submit pull requests or open issues for bugs and feature requests.

## Contact

For questions and support:
- Email: Aleksey Komissarov <ad3002@gmail.com>
- GitHub Issues: https://github.com/aglabx/varprofiler/issues