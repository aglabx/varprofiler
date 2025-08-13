/**
 * Motif Variant Discovery Tool
 * 
 * Discovers all sequence variants within a given edit distance from a reference motif
 * and reports their frequencies in the genome.
 * 
 * Usage: ./motif_discovery <genome.fasta> <motif> <max_distance> <threads> <output.csv>
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <thread>
#include <mutex>
#include <cstring>
#include <iomanip>
#include <queue>
#include <set>
#include <atomic>
#include <chrono>

// IUPAC nucleotide codes
const std::unordered_map<char, std::string> IUPAC_CODES = {
    {'A', "A"}, {'C', "C"}, {'G', "G"}, {'T', "T"},
    {'R', "AG"}, {'Y', "CT"}, {'S', "GC"}, {'W', "AT"},
    {'K', "GT"}, {'M', "AC"}, {'B', "CGT"}, {'D', "AGT"},
    {'H', "ACT"}, {'V', "ACG"}, {'N', "ACGT"}
};

// Structure to hold motif variant statistics
struct MotifStats {
    std::string sequence;
    int edit_distance;
    std::string cigar;  // CIGAR string for alignment
    size_t count_forward = 0;
    size_t count_reverse = 0;
    size_t count_total = 0;
};

// Mutex for thread-safe updates
std::mutex stats_mutex;
std::mutex progress_mutex;

// Global map to store all discovered variants
std::map<std::string, MotifStats> global_stats;

// Progress tracking
std::atomic<size_t> bases_processed(0);
size_t total_bases = 0;

// Display progress bar
void show_progress() {
    const int bar_width = 50;
    auto start_time = std::chrono::steady_clock::now();
    
    while (bases_processed < total_bases) {
        size_t current = bases_processed.load();
        double progress = (total_bases > 0) ? (double)current / total_bases : 0;
        
        // Calculate elapsed time and estimate remaining
        auto now = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start_time).count();
        int eta = (progress > 0) ? (int)(elapsed / progress - elapsed) : 0;
        
        // Create progress bar
        std::cout << "\r[";
        int pos = bar_width * progress;
        for (int i = 0; i < bar_width; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        
        // Show percentage and stats
        std::cout << "] " << std::fixed << std::setprecision(1) 
                  << (progress * 100.0) << "% "
                  << "(" << (current / 1000000) << "/" << (total_bases / 1000000) << " Mb) ";
        
        if (elapsed > 0 && progress > 0.01) {
            std::cout << "ETA: ";
            if (eta >= 3600) {
                std::cout << (eta / 3600) << "h " << ((eta % 3600) / 60) << "m";
            } else if (eta >= 60) {
                std::cout << (eta / 60) << "m " << (eta % 60) << "s";
            } else {
                std::cout << eta << "s";
            }
        }
        
        std::cout << "    " << std::flush;
        
        if (current >= total_bases) break;
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    
    // Final update
    std::cout << "\r[";
    for (int i = 0; i < bar_width; ++i) std::cout << "=";
    std::cout << "] 100.0% (" << (total_bases / 1000000) << "/" << (total_bases / 1000000) << " Mb) Done!        \n";
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

// Generate CIGAR string and calculate edit distance
std::pair<int, std::string> align_and_distance(const std::string& ref, const std::string& query) {
    int m = ref.length();
    int n = query.length();
    
    if (m != n) return {-1, ""}; // Only handle equal length for this application
    
    int distance = 0;
    std::string cigar;
    int match_count = 0;
    
    for (int i = 0; i < m; ++i) {
        if (toupper(ref[i]) == toupper(query[i])) {
            match_count++;
        } else {
            if (match_count > 0) {
                cigar += std::to_string(match_count) + "M";
                match_count = 0;
            }
            cigar += "1X";  // X for mismatch
            distance++;
        }
    }
    
    // Add final matches
    if (match_count > 0) {
        cigar += std::to_string(match_count) + "M";
    }
    
    // If perfect match, use simple CIGAR
    if (distance == 0) {
        cigar = std::to_string(m) + "M";
    }
    
    return {distance, cigar};
}

// Expand IUPAC pattern to all possible sequences
std::vector<std::string> expand_iupac(const std::string& pattern) {
    std::vector<std::string> result = {""};
    
    for (char c : pattern) {
        char uc = toupper(c);
        auto it = IUPAC_CODES.find(uc);
        std::string options = (it != IUPAC_CODES.end()) ? it->second : "N";
        
        std::vector<std::string> new_result;
        for (const auto& prefix : result) {
            for (char nucleotide : options) {
                new_result.push_back(prefix + nucleotide);
            }
        }
        result = new_result;
    }
    
    return result;
}

// Process a sequence chunk to find motif variants
void process_sequence(const std::string& seq, const std::vector<std::string>& expanded_patterns,
                      int max_distance, bool count_reverse) {
    
    std::map<std::string, MotifStats> local_stats;
    int pattern_len = expanded_patterns[0].length();
    
    // Slide through the sequence
    for (size_t i = 0; i <= seq.length() - pattern_len; ++i) {
        // Update progress every 10000 bases
        if (i % 10000 == 0) {
            bases_processed.fetch_add(10000);
        }
        std::string window = seq.substr(i, pattern_len);
        
        // Skip windows with N
        if (window.find('N') != std::string::npos || window.find('n') != std::string::npos) {
            continue;
        }
        
        // Convert to uppercase for comparison
        std::transform(window.begin(), window.end(), window.begin(), ::toupper);
        
        // Check against all expanded patterns
        int min_dist = max_distance + 1;
        std::string best_cigar;
        for (const auto& pattern : expanded_patterns) {
            auto [dist, cigar] = align_and_distance(pattern, window);
            if (dist >= 0 && dist < min_dist) {
                min_dist = dist;
                best_cigar = cigar;
            }
        }
        
        // If within threshold, record it
        if (min_dist <= max_distance) {
            if (local_stats.find(window) == local_stats.end()) {
                local_stats[window].sequence = window;
                local_stats[window].edit_distance = min_dist;
                local_stats[window].cigar = best_cigar;
            }
            
            if (count_reverse) {
                local_stats[window].count_reverse++;
            } else {
                local_stats[window].count_forward++;
            }
        }
    }
    
    // Update final progress for this sequence
    bases_processed.fetch_add(seq.length() % 10000);
    
    // Merge local stats into global
    std::lock_guard<std::mutex> lock(stats_mutex);
    for (const auto& pair : local_stats) {
        const std::string& seq = pair.first;
        const MotifStats& stats = pair.second;
        
        if (global_stats.find(seq) == global_stats.end()) {
            global_stats[seq] = stats;
        } else {
            global_stats[seq].count_forward += stats.count_forward;
            global_stats[seq].count_reverse += stats.count_reverse;
            // CIGAR and edit_distance should already be the same for same sequence
        }
    }
}

// Process a FASTA file chunk
void process_fasta_chunk(const std::vector<std::pair<std::string, std::string>>& sequences,
                         const std::vector<std::string>& expanded_patterns,
                         int max_distance) {
    
    for (const auto& seq_pair : sequences) {
        const std::string& seq = seq_pair.second;
        
        // Process forward strand
        process_sequence(seq, expanded_patterns, max_distance, false);
        
        // Process reverse strand
        std::string rc = reverse_complement(seq);
        process_sequence(rc, expanded_patterns, max_distance, true);
    }
}

// Read FASTA file
std::vector<std::pair<std::string, std::string>> read_fasta(const std::string& filename) {
    std::vector<std::pair<std::string, std::string>> sequences;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return sequences;
    }
    
    std::string line, name, seq;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        if (line[0] == '>') {
            if (!name.empty() && !seq.empty()) {
                sequences.push_back({name, seq});
            }
            name = line.substr(1);
            seq.clear();
        } else {
            seq += line;
        }
    }
    
    if (!name.empty() && !seq.empty()) {
        sequences.push_back({name, seq});
    }
    
    file.close();
    return sequences;
}

// Write results to TSV
void write_tsv(const std::string& filename, const std::map<std::string, MotifStats>& stats) {
    
    // Sort by total count (descending)
    std::vector<std::pair<std::string, MotifStats>> sorted_stats(stats.begin(), stats.end());
    std::sort(sorted_stats.begin(), sorted_stats.end(),
              [](const auto& a, const auto& b) {
                  size_t total_a = a.second.count_forward + a.second.count_reverse;
                  size_t total_b = b.second.count_forward + b.second.count_reverse;
                  return total_a > total_b;
              });
    
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot create output file " << filename << std::endl;
        return;
    }
    
    // Write header
    out << "motif\tedit_distance\tcigar\tcount_forward\tcount_reverse\tcount_total\tfreq_forward\tfreq_reverse\tfreq_total\n";
    
    // Calculate total counts for frequency calculation
    size_t total_forward = 0, total_reverse = 0;
    for (const auto& pair : sorted_stats) {
        total_forward += pair.second.count_forward;
        total_reverse += pair.second.count_reverse;
    }
    size_t total_all = total_forward + total_reverse;
    
    // Write data
    for (const auto& pair : sorted_stats) {
        const std::string& seq = pair.first;
        const MotifStats& stat = pair.second;
        size_t total = stat.count_forward + stat.count_reverse;
        
        double freq_forward = (total_forward > 0) ? (100.0 * stat.count_forward / total_forward) : 0;
        double freq_reverse = (total_reverse > 0) ? (100.0 * stat.count_reverse / total_reverse) : 0;
        double freq_total = (total_all > 0) ? (100.0 * total / total_all) : 0;
        
        out << seq << "\t"
            << stat.edit_distance << "\t"
            << stat.cigar << "\t"
            << stat.count_forward << "\t"
            << stat.count_reverse << "\t"
            << total << "\t"
            << std::fixed << std::setprecision(4)
            << freq_forward << "\t"
            << freq_reverse << "\t"
            << freq_total << "\n";
    }
    
    out.close();
    
    std::cout << "\nTop 10 most frequent variants:\n";
    std::cout << std::setw(20) << "Motif" 
              << std::setw(10) << "Distance"
              << std::setw(15) << "CIGAR"
              << std::setw(15) << "Total Count"
              << std::setw(15) << "Frequency (%)\n";
    std::cout << std::string(75, '-') << "\n";
    
    int count = 0;
    for (const auto& pair : sorted_stats) {
        if (count++ >= 10) break;
        
        size_t total = pair.second.count_forward + pair.second.count_reverse;
        double freq = (total_all > 0) ? (100.0 * total / total_all) : 0;
        
        std::cout << std::setw(20) << pair.first
                  << std::setw(10) << pair.second.edit_distance
                  << std::setw(15) << pair.second.cigar
                  << std::setw(15) << total
                  << std::setw(15) << std::fixed << std::setprecision(2) << freq << "\n";
    }
}

void print_usage(const char* program_name) {
    std::cerr << "Usage: " << program_name << " <genome.fasta> <motif> <max_distance> <threads> <output.tsv>\n\n";
    std::cerr << "Discover all motif variants within edit distance and report their frequencies.\n\n";
    std::cerr << "Arguments:\n";
    std::cerr << "  genome.fasta   Input genome in FASTA format\n";
    std::cerr << "  motif          Reference motif (supports IUPAC codes)\n";
    std::cerr << "  max_distance   Maximum edit distance from reference\n";
    std::cerr << "  threads        Number of threads to use\n";
    std::cerr << "  output.tsv     Output TSV file with variant statistics\n\n";
    std::cerr << "Output columns:\n";
    std::cerr << "  motif          Discovered sequence variant\n";
    std::cerr << "  edit_distance  Distance from reference motif\n";
    std::cerr << "  cigar          CIGAR string showing alignment\n";
    std::cerr << "  count_*        Occurrence counts\n";
    std::cerr << "  freq_*         Percentage frequencies\n\n";
    std::cerr << "Example:\n";
    std::cerr << "  " << program_name << " genome.fa YTTCGTTGGAARCGGGA 3 8 cenpb_variants.tsv\n";
}

int main(int argc, char* argv[]) {
    if (argc != 6) {
        print_usage(argv[0]);
        return 1;
    }
    
    std::string genome_file = argv[1];
    std::string motif = argv[2];
    int max_distance = std::stoi(argv[3]);
    int num_threads = std::stoi(argv[4]);
    std::string output_file = argv[5];
    
    if (max_distance < 0 || max_distance > 10) {
        std::cerr << "Error: max_distance must be between 0 and 10\n";
        return 1;
    }
    
    if (num_threads < 1) {
        std::cerr << "Error: threads must be at least 1\n";
        return 1;
    }
    
    std::cout << "Motif Variant Discovery\n";
    std::cout << "=======================\n";
    std::cout << "Reference motif: " << motif << "\n";
    std::cout << "Max edit distance: " << max_distance << "\n";
    std::cout << "Threads: " << num_threads << "\n\n";
    
    // Expand IUPAC codes in the reference motif
    std::cout << "Expanding IUPAC codes...\n";
    std::vector<std::string> expanded_patterns = expand_iupac(motif);
    std::cout << "Expanded to " << expanded_patterns.size() << " pattern(s)\n";
    
    // Read genome
    std::cout << "Reading genome from " << genome_file << "...\n";
    auto sequences = read_fasta(genome_file);
    if (sequences.empty()) {
        std::cerr << "Error: No sequences found in FASTA file\n";
        return 1;
    }
    std::cout << "Loaded " << sequences.size() << " sequence(s)\n";
    
    // Calculate total genome size
    size_t genome_size = 0;
    for (const auto& seq : sequences) {
        genome_size += seq.second.length();
    }
    std::cout << "Total genome size: " << (genome_size / 1000000) << " Mb\n\n";
    
    // Initialize progress tracking
    total_bases = genome_size * 2;  // Forward and reverse strands
    bases_processed = 0;
    
    // Process sequences in parallel
    std::cout << "Searching for motif variants...\n";
    
    // Start progress bar thread
    std::thread progress_thread(show_progress);
    
    // Divide work among threads
    size_t chunk_size = (sequences.size() + num_threads - 1) / num_threads;
    std::vector<std::thread> threads;
    
    for (int t = 0; t < num_threads; ++t) {
        size_t start = t * chunk_size;
        size_t end = std::min(start + chunk_size, sequences.size());
        
        if (start < sequences.size()) {
            threads.emplace_back([&sequences, start, end, &expanded_patterns, max_distance]() {
                std::vector<std::pair<std::string, std::string>> chunk(
                    sequences.begin() + start, sequences.begin() + end);
                process_fasta_chunk(chunk, expanded_patterns, max_distance);
            });
        }
    }
    
    // Wait for all threads
    for (auto& t : threads) {
        t.join();
    }
    
    // Mark as complete and wait for progress bar
    bases_processed = total_bases;
    progress_thread.join();
    
    // Calculate total counts for all variants
    for (auto& pair : global_stats) {
        pair.second.count_total = pair.second.count_forward + pair.second.count_reverse;
    }
    
    std::cout << "\nFound " << global_stats.size() << " unique variants\n";
    
    // Write results
    std::cout << "Writing results to " << output_file << "...\n";
    write_tsv(output_file, global_stats);
    
    std::cout << "\nDone!\n";
    
    return 0;
}