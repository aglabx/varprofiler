/**
 * CENP-B Box Finder
 * 
 * Searches for CENP-B box motifs in genomic sequences using edit distance.
 * The canonical CENP-B box is YTTCGTTGGAARCGGGA (17bp) but can vary in length
 * and sequence across species. Conserved regions are typically at the ends.
 * 
 * Uses Myers' bit-parallel algorithm for efficient edit distance calculation.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cctype>
#include <cstring>
#include <iomanip>
#include <unordered_map>
#include <thread>
#include <mutex>
#include <atomic>

// IUPAC nucleotide codes for degenerate bases
std::unordered_map<char, std::string> IUPAC = {
    {'A', "A"}, {'C', "C"}, {'G', "G"}, {'T', "T"},
    {'R', "AG"}, {'Y', "CT"}, {'S', "GC"}, {'W', "AT"},
    {'K', "GT"}, {'M', "AC"}, {'B', "CGT"}, {'D', "AGT"},
    {'H', "ACT"}, {'V', "ACG"}, {'N', "ACGT"}
};

// Configuration
struct Config {
    std::string fasta_file;
    std::string bed_file;      // BED file with low k-mer diversity regions
    std::string output_file;
    std::string pattern = "YTTCGTTGGAARCGGGA";  // Default CENP-B box pattern
    int max_edit_distance = 3;  // Maximum allowed edit distance
    int threads = 1;
    bool search_both_strands = true;
    bool verbose = false;
};

// Region from BED file
struct Region {
    std::string chrom;
    size_t start;
    size_t end;
    int kmer_count;
    int cenpb_count = 0;
    std::vector<size_t> cenpb_positions;
};

// Expand IUPAC pattern to all possible sequences
std::vector<std::string> expand_iupac_pattern(const std::string& pattern) {
    std::vector<std::string> result = {""};
    
    for (char c : pattern) {
        char upper = std::toupper(c);
        std::string options = IUPAC.count(upper) ? IUPAC[upper] : std::string(1, upper);
        
        std::vector<std::string> new_result;
        for (const auto& prefix : result) {
            for (char option : options) {
                new_result.push_back(prefix + option);
            }
        }
        result = new_result;
        
        // Limit expansion to prevent combinatorial explosion
        if (result.size() > 100) {
            std::cerr << "Warning: Pattern expansion limited to first 100 variants\n";
            result.resize(100);
            break;
        }
    }
    
    return result;
}

// Reverse complement of DNA sequence
std::string reverse_complement(const std::string& seq) {
    std::string rc(seq.length(), 'N');
    for (size_t i = 0; i < seq.length(); ++i) {
        char c = seq[seq.length() - 1 - i];
        switch (c) {
            case 'A': case 'a': rc[i] = 'T'; break;
            case 'T': case 't': rc[i] = 'A'; break;
            case 'G': case 'g': rc[i] = 'C'; break;
            case 'C': case 'c': rc[i] = 'G'; break;
            default: rc[i] = 'N';
        }
    }
    return rc;
}

// Myers' bit-parallel algorithm for edit distance
// Returns minimum edit distance between pattern and any substring of text
int myers_edit_distance(const std::string& text, const std::string& pattern, int max_dist) {
    const int m = pattern.length();
    const int n = text.length();
    
    if (m > 64) {
        // Fall back to standard DP for long patterns
        return -1;  // Not implemented for simplicity
    }
    
    // Preprocessing
    uint64_t Peq[256] = {0};
    for (int i = 0; i < m; ++i) {
        Peq[(unsigned char)pattern[i]] |= (1ULL << i);
    }
    
    int min_distance = m;
    uint64_t Pv = ~0ULL;
    uint64_t Mv = 0;
    int score = m;
    
    for (int j = 0; j < n; ++j) {
        uint64_t Eq = Peq[(unsigned char)text[j]];
        uint64_t Xv = Eq | Mv;
        uint64_t Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
        
        uint64_t Ph = Mv | ~(Xh | Pv);
        uint64_t Mh = Pv & Xh;
        
        if (Ph & (1ULL << (m - 1))) ++score;
        if (Mh & (1ULL << (m - 1))) --score;
        
        Ph <<= 1;
        Mh <<= 1;
        Pv = Mh | ~(Xv | Ph);
        Mv = Ph & Xv;
        
        if (score <= max_dist) {
            min_distance = std::min(min_distance, score);
        }
    }
    
    return min_distance;
}

// Simple edit distance for when Myers' algorithm can't be used
int simple_edit_distance(const std::string& text, const std::string& pattern, int max_dist) {
    const int m = pattern.length();
    const int n = text.length();
    int min_distance = m;
    
    // Check each possible alignment
    for (int start = 0; start <= n - m + max_dist; ++start) {
        int dist = 0;
        int text_pos = start;
        
        for (int i = 0; i < m && text_pos < n; ++i) {
            if (pattern[i] != text[text_pos]) {
                dist++;
                if (dist > max_dist) break;
            }
            text_pos++;
        }
        
        // Account for remaining pattern characters if text ended
        if (text_pos < n || (text_pos == n && start + m <= n)) {
            dist += std::max(0, m - (text_pos - start));
        }
        
        if (dist <= max_dist) {
            min_distance = std::min(min_distance, dist);
        }
    }
    
    return min_distance;
}

// Find all approximate matches of pattern in sequence
std::vector<size_t> find_pattern_matches(const std::string& seq, 
                                        const std::vector<std::string>& patterns,
                                        int max_dist,
                                        bool both_strands) {
    std::vector<size_t> matches;
    
    // Convert sequence to uppercase
    std::string seq_upper(seq);
    std::transform(seq_upper.begin(), seq_upper.end(), seq_upper.begin(), ::toupper);
    
    // Search with each pattern variant
    for (const auto& pattern : patterns) {
        // Forward strand
        for (size_t i = 0; i <= seq_upper.length() - pattern.length() + max_dist; ++i) {
            if (i + pattern.length() - max_dist > seq_upper.length()) break;
            
            int window_size = std::min(pattern.length() + max_dist, seq_upper.length() - i);
            std::string window = seq_upper.substr(i, window_size);
            
            int dist;
            if (pattern.length() <= 64) {
                dist = myers_edit_distance(window, pattern, max_dist);
            } else {
                dist = simple_edit_distance(window, pattern, max_dist);
            }
            
            if (dist <= max_dist) {
                matches.push_back(i);
                i += pattern.length() / 2;  // Skip ahead to avoid overlapping matches
            }
        }
        
        // Reverse strand
        if (both_strands) {
            std::string pattern_rc = reverse_complement(pattern);
            for (size_t i = 0; i <= seq_upper.length() - pattern_rc.length() + max_dist; ++i) {
                if (i + pattern_rc.length() - max_dist > seq_upper.length()) break;
                
                int window_size = std::min(pattern_rc.length() + max_dist, seq_upper.length() - i);
                std::string window = seq_upper.substr(i, window_size);
                
                int dist;
                if (pattern_rc.length() <= 64) {
                    dist = myers_edit_distance(window, pattern_rc, max_dist);
                } else {
                    dist = simple_edit_distance(window, pattern_rc, max_dist);
                }
                
                if (dist <= max_dist) {
                    matches.push_back(i);
                    i += pattern_rc.length() / 2;  // Skip ahead
                }
            }
        }
    }
    
    // Remove duplicates and sort
    std::sort(matches.begin(), matches.end());
    matches.erase(std::unique(matches.begin(), matches.end()), matches.end());
    
    return matches;
}

// Read BED file with regions
std::vector<Region> read_bed_file(const std::string& filename) {
    std::vector<Region> regions;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open BED file: " << filename << std::endl;
        return regions;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        Region region;
        iss >> region.chrom >> region.start >> region.end;
        
        if (iss >> region.kmer_count) {
            // Optional 4th column with k-mer count
        }
        
        regions.push_back(region);
    }
    
    return regions;
}

// Read FASTA file into memory
std::unordered_map<std::string, std::string> read_fasta(const std::string& filename) {
    std::unordered_map<std::string, std::string> sequences;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open FASTA file: " << filename << std::endl;
        return sequences;
    }
    
    std::string line, current_seq_name;
    std::string current_seq;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        if (line[0] == '>') {
            // Save previous sequence
            if (!current_seq_name.empty()) {
                sequences[current_seq_name] = current_seq;
            }
            
            // Parse new sequence name
            current_seq_name = line.substr(1);
            size_t space_pos = current_seq_name.find(' ');
            if (space_pos != std::string::npos) {
                current_seq_name = current_seq_name.substr(0, space_pos);
            }
            current_seq.clear();
        } else {
            current_seq += line;
        }
    }
    
    // Save last sequence
    if (!current_seq_name.empty()) {
        sequences[current_seq_name] = current_seq;
    }
    
    return sequences;
}

// Process a batch of regions in a thread
void process_regions_thread(const std::vector<Region>& regions,
                           size_t start_idx,
                           size_t end_idx,
                           const std::unordered_map<std::string, std::string>& sequences,
                           const std::vector<std::string>& patterns,
                           const Config& config,
                           std::vector<Region>& results,
                           std::mutex& results_mutex,
                           std::atomic<size_t>& progress) {
    
    for (size_t i = start_idx; i < end_idx; ++i) {
        Region region = regions[i];
        
        // Get sequence for this region
        auto it = sequences.find(region.chrom);
        if (it == sequences.end()) {
            std::cerr << "Warning: Chromosome " << region.chrom << " not found in FASTA\n";
            continue;
        }
        
        // Extract region sequence
        if (region.end > it->second.length()) {
            std::cerr << "Warning: Region " << region.chrom << ":" << region.start 
                     << "-" << region.end << " extends beyond sequence length\n";
            continue;
        }
        
        std::string region_seq = it->second.substr(region.start, region.end - region.start);
        
        // Find CENP-B boxes
        region.cenpb_positions = find_pattern_matches(region_seq, patterns, 
                                                     config.max_edit_distance,
                                                     config.search_both_strands);
        region.cenpb_count = region.cenpb_positions.size();
        
        // Convert positions to genomic coordinates
        for (auto& pos : region.cenpb_positions) {
            pos += region.start;
        }
        
        // Store result
        {
            std::lock_guard<std::mutex> lock(results_mutex);
            results.push_back(region);
        }
        
        progress++;
        
        if (config.verbose && progress % 100 == 0) {
            std::cerr << "Processed " << progress << " regions...\r";
        }
    }
}

// Main processing function
void process_regions(const Config& config) {
    // Read input files
    std::vector<Region> regions = read_bed_file(config.bed_file);
    if (regions.empty()) {
        std::cerr << "Error: No regions found in BED file\n";
        return;
    }
    std::cout << "Loaded " << regions.size() << " regions from BED file\n";
    
    auto sequences = read_fasta(config.fasta_file);
    if (sequences.empty()) {
        std::cerr << "Error: No sequences found in FASTA file\n";
        return;
    }
    std::cout << "Loaded " << sequences.size() << " sequences from FASTA file\n";
    
    // Expand IUPAC pattern
    auto patterns = expand_iupac_pattern(config.pattern);
    std::cout << "Searching for " << patterns.size() << " pattern variant(s) with max edit distance " 
              << config.max_edit_distance << "\n";
    
    // Process regions in parallel
    std::vector<Region> results;
    std::mutex results_mutex;
    std::atomic<size_t> progress(0);
    
    size_t n_threads = config.threads;
    size_t chunk_size = (regions.size() + n_threads - 1) / n_threads;
    
    std::vector<std::thread> threads;
    for (size_t t = 0; t < n_threads; ++t) {
        size_t start_idx = t * chunk_size;
        size_t end_idx = std::min(start_idx + chunk_size, regions.size());
        
        if (start_idx < regions.size()) {
            threads.emplace_back(process_regions_thread,
                               std::ref(regions), start_idx, end_idx,
                               std::ref(sequences), std::ref(patterns),
                               std::ref(config), std::ref(results),
                               std::ref(results_mutex), std::ref(progress));
        }
    }
    
    // Wait for threads
    for (auto& t : threads) {
        t.join();
    }
    
    if (config.verbose) {
        std::cerr << "\n";
    }
    
    // Sort results by chromosome and position
    std::sort(results.begin(), results.end(), [](const Region& a, const Region& b) {
        if (a.chrom != b.chrom) return a.chrom < b.chrom;
        return a.start < b.start;
    });
    
    // Write output
    std::ofstream out(config.output_file);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot open output file: " << config.output_file << std::endl;
        return;
    }
    
    // Write header
    out << "#chrom\tstart\tend\tkmer_count\tcenpb_count\tcenpb_density\tcenpb_positions\n";
    
    // Write results
    size_t total_regions_with_cenpb = 0;
    size_t total_cenpb_boxes = 0;
    
    for (const auto& region : results) {
        double density = (double)region.cenpb_count / (region.end - region.start) * 1000;  // per kb
        
        out << region.chrom << "\t" 
            << region.start << "\t" 
            << region.end << "\t"
            << region.kmer_count << "\t"
            << region.cenpb_count << "\t"
            << std::fixed << std::setprecision(4) << density;
        
        // Output positions of CENP-B boxes
        if (!region.cenpb_positions.empty()) {
            out << "\t";
            for (size_t i = 0; i < region.cenpb_positions.size(); ++i) {
                if (i > 0) out << ",";
                out << region.cenpb_positions[i];
            }
            total_regions_with_cenpb++;
        } else {
            out << "\t.";
        }
        
        out << "\n";
        
        total_cenpb_boxes += region.cenpb_count;
    }
    
    out.close();
    
    // Print summary
    std::cout << "\nResults written to: " << config.output_file << "\n";
    std::cout << "Summary:\n";
    std::cout << "  Total regions analyzed: " << results.size() << "\n";
    std::cout << "  Regions with CENP-B boxes: " << total_regions_with_cenpb 
              << " (" << (100.0 * total_regions_with_cenpb / results.size()) << "%)\n";
    std::cout << "  Total CENP-B boxes found: " << total_cenpb_boxes << "\n";
    
    if (total_regions_with_cenpb > 0) {
        std::cout << "  Average CENP-B boxes per positive region: " 
                  << (double)total_cenpb_boxes / total_regions_with_cenpb << "\n";
    }
}

// Print usage
void print_usage(const char* program_name) {
    std::cerr << "Usage: " << program_name << " <genome.fasta> <regions.bed> <output.tsv> [options]\n\n";
    std::cerr << "Search for CENP-B box motifs in specified genomic regions using edit distance.\n\n";
    std::cerr << "Arguments:\n";
    std::cerr << "  genome.fasta   Input genome in FASTA format\n";
    std::cerr << "  regions.bed    BED file with regions to search (e.g., low k-mer diversity regions)\n";
    std::cerr << "  output.tsv     Output file with CENP-B box counts and positions\n\n";
    std::cerr << "Options:\n";
    std::cerr << "  -p PATTERN     CENP-B box pattern (default: YTTCGTTGGAARCGGGA)\n";
    std::cerr << "  -d DISTANCE    Maximum edit distance (default: 3)\n";
    std::cerr << "  -t THREADS     Number of threads (default: 1)\n";
    std::cerr << "  -s             Search forward strand only (default: both strands)\n";
    std::cerr << "  -v             Verbose output\n";
    std::cerr << "  -h             Show this help message\n\n";
    std::cerr << "Example:\n";
    std::cerr << "  " << program_name << " genome.fa low_kmer_regions.bed cenpb_boxes.tsv -d 2 -t 4\n";
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        print_usage(argv[0]);
        return 1;
    }
    
    Config config;
    config.fasta_file = argv[1];
    config.bed_file = argv[2];
    config.output_file = argv[3];
    
    // Parse optional arguments
    for (int i = 4; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "-p" && i + 1 < argc) {
            config.pattern = argv[++i];
        } else if (arg == "-d" && i + 1 < argc) {
            config.max_edit_distance = std::atoi(argv[++i]);
        } else if (arg == "-t" && i + 1 < argc) {
            config.threads = std::atoi(argv[++i]);
        } else if (arg == "-s") {
            config.search_both_strands = false;
        } else if (arg == "-v") {
            config.verbose = true;
        } else if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }
    
    // Validate configuration
    if (config.max_edit_distance < 0 || config.max_edit_distance > 10) {
        std::cerr << "Error: Edit distance must be between 0 and 10\n";
        return 1;
    }
    
    if (config.threads < 1 || config.threads > 100) {
        std::cerr << "Error: Thread count must be between 1 and 100\n";
        return 1;
    }
    
    // Run analysis
    std::cout << "CENP-B Box Finder\n";
    std::cout << "==================\n";
    std::cout << "Pattern: " << config.pattern << "\n";
    std::cout << "Max edit distance: " << config.max_edit_distance << "\n";
    std::cout << "Threads: " << config.threads << "\n";
    std::cout << "Search strands: " << (config.search_both_strands ? "both" : "forward only") << "\n\n";
    
    process_regions(config);
    
    return 0;
}