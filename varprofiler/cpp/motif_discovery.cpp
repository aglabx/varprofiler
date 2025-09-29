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
#include <unordered_set>
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
#include <cstdlib>
#include <cmath>
#ifdef __SSE2__
#include <emmintrin.h>
#endif

// IUPAC nucleotide codes
const std::unordered_map<char, std::string> IUPAC_CODES = {
    {'A', "A"}, {'C', "C"}, {'G', "G"}, {'T', "T"},
    {'R', "AG"}, {'Y', "CT"}, {'S', "GC"}, {'W', "AT"},
    {'K', "GT"}, {'M', "AC"}, {'B', "CGT"}, {'D', "AGT"},
    {'H', "ACT"}, {'V', "ACG"}, {'N', "ACGT"}
};

// Check if a nucleotide matches an IUPAC code
inline bool iupac_match(char iupac_code, char nucleotide) {
    char iupac_upper = toupper(iupac_code);
    char nuc_upper = toupper(nucleotide);
    
    // Handle N separately for any non-IUPAC characters
    if (nuc_upper == 'N') return true;  // N in sequence matches anything
    if (iupac_upper == 'N') return true;  // N in pattern matches anything
    
    // Direct match
    if (iupac_upper == nuc_upper) return true;
    
    // Check IUPAC ambiguity codes
    auto it = IUPAC_CODES.find(iupac_upper);
    if (it != IUPAC_CODES.end()) {
        return it->second.find(nuc_upper) != std::string::npos;
    }
    
    return false;
}

// Structure to hold motif variant statistics
struct MotifStats {
    std::string sequence;
    std::string alignment;  // Visual alignment (dots for matches)
    int edit_distance;
    std::string cigar;  // CIGAR string for alignment
    size_t count_forward = 0;
    size_t count_reverse = 0;
    size_t count_total = 0;
};

// Structure to hold individual motif occurrences
struct MotifOccurrence {
    std::string chr;
    size_t start;
    size_t end;
    std::string motif;
    int edit_distance;
    char strand;  // '+' or '-'
};

// Mutex for thread-safe updates
std::mutex stats_mutex;
std::mutex progress_mutex;
std::mutex occurrences_mutex;

// Global map to store all discovered variants
std::map<std::string, MotifStats> global_stats;

// Global vector to store all motif occurrences
std::vector<MotifOccurrence> global_occurrences;

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

// Structure to hold alignment result
struct AlignmentResult {
    int edit_distance;
    std::string cigar;
    std::string visual_alignment;
};

// Calculate edit distance with full alignment (including indels)
AlignmentResult align_sequences(const std::string& ref, const std::string& query, int max_dist) {
    int m = ref.length();
    int n = query.length();
    
    // Quick check: if length difference exceeds max_dist, skip
    if (abs(m - n) > max_dist) {
        return {abs(m - n), "", ""};
    }
    
    // OPTIMIZATION 1: Banded DP - only compute cells within max_dist of diagonal
    // This reduces complexity from O(m*n) to O(m*max_dist)
    
    // Use single row DP for space optimization when we only need distance
    if (m == n && max_dist <= 3) {
        // OPTIMIZATION 2: For small distances and equal lengths, use IUPAC-aware Hamming distance
        int hamming = 0;
        for (int i = 0; i < m; i++) {
            if (!iupac_match(ref[i], query[i])) hamming++;
            if (hamming > max_dist) return {hamming, "", ""};  // Early termination
        }
        // Only build alignment if within threshold
        if (hamming <= max_dist) {
            // Fast path for equal length - build CIGAR and visual
            std::string cigar, visual;
            int match_count = 0;
            for (int i = 0; i < m; i++) {
                if (iupac_match(ref[i], query[i])) {
                    match_count++;
                    visual += '.';
                } else {
                    if (match_count > 0) {
                        cigar += std::to_string(match_count) + "M";
                        match_count = 0;
                    }
                    cigar += "1X";
                    visual += toupper(query[i]);
                }
            }
            if (match_count > 0) cigar += std::to_string(match_count) + "M";
            return {hamming, cigar, visual};
        }
        return {hamming, "", ""};
    }
    
    // Full DP for unequal lengths or large distances
    // OPTIMIZATION 3: Use Ukkonen's cutoff - stop if we exceed max_dist
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, max_dist + 1));
    
    // Initialize first row and column
    for (int i = 0; i <= m && i <= max_dist; i++) dp[i][0] = i;
    for (int j = 0; j <= n && j <= max_dist; j++) dp[0][j] = j;
    
    // Fill DP table with banded optimization
    bool exceeded = false;
    for (int i = 1; i <= m; i++) {
        int start_j = std::max(1, i - max_dist);
        int end_j = std::min(n, i + max_dist);
        bool row_exceeded = true;
        
        for (int j = start_j; j <= end_j; j++) {
            if (iupac_match(ref[i-1], query[j-1])) {
                dp[i][j] = dp[i-1][j-1];  // Match
            } else {
                int del = (j > 1 && i - j + 1 <= max_dist) ? dp[i-1][j] : max_dist + 1;
                int ins = (i > 1 && j - i + 1 <= max_dist) ? dp[i][j-1] : max_dist + 1;
                int sub = dp[i-1][j-1];
                dp[i][j] = 1 + std::min({del, ins, sub});
            }
            if (dp[i][j] <= max_dist) row_exceeded = false;
        }
        
        // Early termination if entire row exceeds max_dist
        if (row_exceeded && i > max_dist) {
            exceeded = true;
            break;
        }
    }
    
    int edit_dist = dp[m][n];
    
    // Early exit if distance too large
    if (edit_dist > max_dist || exceeded) {
        return {edit_dist, "", ""};
    }
    
    // Backtrack to build CIGAR and visual alignment
    std::string cigar;
    std::string visual;
    std::vector<char> cigar_ops;
    std::vector<char> visual_chars;
    
    int i = m, j = n;
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && iupac_match(ref[i-1], query[j-1])) {
            // Match
            cigar_ops.push_back('M');
            visual_chars.push_back('.');
            i--; j--;
        } else if (i > 0 && j > 0 && dp[i][j] == dp[i-1][j-1] + 1) {
            // Substitution
            cigar_ops.push_back('X');
            visual_chars.push_back(toupper(query[j-1]));
            i--; j--;
        } else if (j > 0 && dp[i][j] == dp[i][j-1] + 1) {
            // Insertion
            cigar_ops.push_back('I');
            visual_chars.push_back(toupper(query[j-1]));
            j--;
        } else if (i > 0 && dp[i][j] == dp[i-1][j] + 1) {
            // Deletion
            cigar_ops.push_back('D');
            visual_chars.push_back('-');
            i--;
        }
    }
    
    // Reverse to get correct order
    std::reverse(cigar_ops.begin(), cigar_ops.end());
    std::reverse(visual_chars.begin(), visual_chars.end());
    
    // Build CIGAR string
    if (!cigar_ops.empty()) {
        char current_op = cigar_ops[0];
        int count = 1;
        
        for (size_t k = 1; k < cigar_ops.size(); k++) {
            if (cigar_ops[k] == current_op) {
                count++;
            } else {
                cigar += std::to_string(count) + current_op;
                current_op = cigar_ops[k];
                count = 1;
            }
        }
        cigar += std::to_string(count) + current_op;
    }
    
    // Build visual string
    for (char c : visual_chars) {
        visual += c;
    }
    
    return {edit_dist, cigar, visual};
}

// Full edit distance with indels (Levenshtein distance)
int levenshtein_distance(const std::string& s1, const std::string& s2) {
    int m = s1.length();
    int n = s2.length();
    
    // Create DP table
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1));
    
    // Initialize first row and column
    for (int i = 0; i <= m; i++) dp[i][0] = i;
    for (int j = 0; j <= n; j++) dp[0][j] = j;
    
    // Fill DP table
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            if (toupper(s1[i-1]) == toupper(s2[j-1])) {
                dp[i][j] = dp[i-1][j-1];  // Match
            } else {
                dp[i][j] = 1 + std::min({
                    dp[i-1][j],    // Deletion
                    dp[i][j-1],    // Insertion
                    dp[i-1][j-1]   // Substitution
                });
            }
        }
    }
    
    return dp[m][n];
}

// Expand IUPAC pattern to include both forward and reverse complement
// We'll use IUPAC matching in align_sequences
std::vector<std::pair<std::string, bool>> expand_iupac_with_rc(const std::string& pattern) {
    std::vector<std::pair<std::string, bool>> result;
    
    // Add forward pattern
    result.push_back({pattern, false});  // false = forward strand
    
    // Add reverse complement pattern
    std::string rc_pattern = reverse_complement(pattern);
    // Only add RC if it's different from forward (avoid palindromes counted twice)
    if (rc_pattern != pattern) {
        result.push_back({rc_pattern, true});  // true = reverse strand
    }
    
    return result;
}

// Process a sequence chunk to find motif variants
void process_sequence(const std::string& seq, const std::string& chr_name,
                      const std::vector<std::pair<std::string, bool>>& patterns_with_strand,
                      int max_distance) {
    
    std::map<std::string, MotifStats> local_stats;
    std::vector<MotifOccurrence> local_occurrences;
    
    if (patterns_with_strand.empty()) return;
    int pattern_len = patterns_with_strand[0].first.length();
    
    
    // OPTIMIZATION: Skip long stretches of Ns
    size_t skip_until = 0;
    
    // Track last match position to avoid overlaps
    size_t last_match_end = 0;
    
    // Slide through the sequence
    for (size_t i = 0; i <= seq.length() - pattern_len; ++i) {
        // Skip if we're in a skip region
        if (i < skip_until) continue;
        
        // Skip if we're still within the last match
        if (i < last_match_end) continue;
        
        // Update progress every 10000 bases
        if (i % 10000 == 0) {
            bases_processed.fetch_add(10000);
        }
        
        // Check for N-rich region and skip it
        if (seq[i] == 'N' || seq[i] == 'n') {
            size_t n_count = 0;
            size_t j = i;
            while (j < seq.length() && (seq[j] == 'N' || seq[j] == 'n')) {
                n_count++;
                j++;
            }
            if (n_count > static_cast<size_t>(pattern_len)) {
                skip_until = j;
                i = j - 1;  // Will be incremented by loop
                continue;
            }
        }
        
        // OPTIMIZATION 5: Check exact length first (most common case)
        // This avoids checking all length variations when not needed
        bool found_match = false;
        
        // Try exact length match first
        if (i + pattern_len <= seq.length()) {
            std::string window = seq.substr(i, pattern_len);
            
            // Quick N check
            bool has_n = false;
            for (char c : window) {
                if (c == 'N' || c == 'n') {
                    has_n = true;
                    break;
                }
            }
            if (!has_n) {
                // Convert to uppercase for comparison
                std::transform(window.begin(), window.end(), window.begin(), ::toupper);
                
                // Check against pattern(s) using IUPAC matching
                {
                    int min_dist = max_distance + 1;
                    std::string best_cigar;
                    std::string best_visual;
                    bool best_is_reverse = false;
                    
                    for (const auto& [pattern, is_reverse] : patterns_with_strand) {
                        auto result = align_sequences(pattern, window, max_distance);
                        if (result.edit_distance >= 0 && result.edit_distance <= max_distance && result.edit_distance < min_dist) {
                            min_dist = result.edit_distance;
                            best_cigar = result.cigar;
                            best_visual = result.visual_alignment;
                            best_is_reverse = is_reverse;
                        }
                    }
                    
                    if (min_dist <= max_distance) {
                        // For reverse strand matches, convert back to forward orientation for consensus
                        std::string consensus_seq = best_is_reverse ? reverse_complement(window) : window;
                        
                        if (local_stats.find(consensus_seq) == local_stats.end()) {
                            local_stats[consensus_seq].sequence = consensus_seq;
                            local_stats[consensus_seq].alignment = best_visual;
                            local_stats[consensus_seq].edit_distance = min_dist;
                            local_stats[consensus_seq].cigar = best_cigar;
                        }
                        if (best_is_reverse) {
                            local_stats[consensus_seq].count_reverse++;
                        } else {
                            local_stats[consensus_seq].count_forward++;
                        }
                        
                        // Record position
                        MotifOccurrence occ;
                        occ.chr = chr_name;
                        occ.start = i;
                        occ.end = i + pattern_len;
                        occ.motif = consensus_seq;  // Use forward-oriented sequence
                        occ.edit_distance = min_dist;
                        occ.strand = best_is_reverse ? '-' : '+';
                        local_occurrences.push_back(occ);
                        
                        found_match = true;
                        last_match_end = i + pattern_len;  // Skip past this match
                    }
                }
            }
        }
        
        // For handling internal indels: check windows of different lengths but verify indel is internal
        if (!found_match && max_distance > 0) {
            // Only check slightly different lengths for internal indels
            for (int len_diff = -1; len_diff <= 1; len_diff++) {
                if (len_diff == 0) continue;  // Already tried exact length
                
                int window_len = pattern_len + len_diff;
                if (window_len < pattern_len - 1 || window_len > pattern_len + 1) continue;
                if (i + window_len > seq.length()) continue;
                
                std::string window = seq.substr(i, window_len);
                
                // Quick N check
                if (window.find('N') != std::string::npos || window.find('n') != std::string::npos) {
                    continue;
                }
                
                std::transform(window.begin(), window.end(), window.begin(), ::toupper);
                
                // Check alignment and verify indel is internal (not at flanks)
                int min_dist = max_distance + 1;
                std::string best_cigar;
                std::string best_visual;
                bool best_is_reverse = false;
                
                for (const auto& [pattern, is_reverse] : patterns_with_strand) {
                    auto result = align_sequences(pattern, window, max_distance);
                    if (result.edit_distance > 0 && result.edit_distance <= max_distance && result.edit_distance < min_dist) {
                        // Check if the indel is internal by examining CIGAR
                        // Flanking indels would start/end with I or D
                        if (!result.cigar.empty()) {
                            bool starts_with_indel = (result.cigar[0] == 'I' || result.cigar[0] == 'D');
                            bool ends_with_indel = false;
                            
                            // Find last operation in CIGAR
                            size_t last_op_pos = result.cigar.find_last_of("MIDX");
                            if (last_op_pos != std::string::npos) {
                                char last_op = result.cigar[last_op_pos];
                                ends_with_indel = (last_op == 'I' || last_op == 'D');
                            }
                            
                            // Skip if indel is at flanks
                            if (starts_with_indel || ends_with_indel) {
                                continue;
                            }
                            
                            // This is an internal indel - accept it
                            min_dist = result.edit_distance;
                            best_cigar = result.cigar;
                            best_visual = result.visual_alignment;
                            best_is_reverse = is_reverse;
                        }
                    }
                }
                
                if (min_dist <= max_distance) {
                    // For reverse strand matches, convert back to forward orientation for consensus
                    std::string consensus_seq = best_is_reverse ? reverse_complement(window) : window;
                    
                    if (local_stats.find(consensus_seq) == local_stats.end()) {
                        local_stats[consensus_seq].sequence = consensus_seq;
                        local_stats[consensus_seq].alignment = best_visual;
                        local_stats[consensus_seq].edit_distance = min_dist;
                        local_stats[consensus_seq].cigar = best_cigar;
                    }
                    if (best_is_reverse) {
                        local_stats[consensus_seq].count_reverse++;
                    } else {
                        local_stats[consensus_seq].count_forward++;
                    }
                    
                    // Record position
                    MotifOccurrence occ;
                    occ.chr = chr_name;
                    occ.start = i;
                    occ.end = i + window_len;
                    occ.motif = consensus_seq;  // Use forward-oriented sequence
                    occ.edit_distance = min_dist;
                    occ.strand = best_is_reverse ? '-' : '+';
                    local_occurrences.push_back(occ);
                    
                    last_match_end = i + window_len;
                    found_match = true;
                }
                
                if (found_match) break;
            }
        }
    }
    
    // Update final progress for this sequence
    bases_processed.fetch_add(seq.length() % 10000);
    
    // Merge local stats into global
    {
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
    
    // Merge local occurrences into global
    {
        std::lock_guard<std::mutex> lock(occurrences_mutex);
        global_occurrences.insert(global_occurrences.end(), 
                                  local_occurrences.begin(), 
                                  local_occurrences.end());
    }
}

// Process a FASTA file chunk
void process_fasta_chunk(const std::vector<std::pair<std::string, std::string>>& sequences,
                         const std::vector<std::pair<std::string, bool>>& patterns_with_strand,
                         int max_distance) {
    
    for (const auto& seq_pair : sequences) {
        const std::string& chr_name = seq_pair.first;
        const std::string& seq = seq_pair.second;
        
        // Process only forward strand, searching for both forward and RC patterns
        process_sequence(seq, chr_name, patterns_with_strand, max_distance);
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

// Write BED file with all motif positions
void write_bed(const std::string& filename, const std::vector<MotifOccurrence>& occurrences) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot create BED file " << filename << std::endl;
        return;
    }
    
    // Sort by chromosome and position
    std::vector<MotifOccurrence> sorted_occ = occurrences;
    std::sort(sorted_occ.begin(), sorted_occ.end(),
              [](const MotifOccurrence& a, const MotifOccurrence& b) {
                  if (a.chr != b.chr) return a.chr < b.chr;
                  return a.start < b.start;
              });
    
    // Write BED entries
    for (const auto& occ : sorted_occ) {
        out << occ.chr << "\t" 
            << occ.start << "\t" 
            << occ.end << "\t"
            << occ.motif << "\t"
            << occ.edit_distance << "\t"
            << occ.strand << "\n";
    }
    
    out.close();
    std::cout << "Wrote " << sorted_occ.size() << " motif positions to BED file\n";
}

// Generate consensus logo matrix
void write_logo_matrix(const std::string& filename, 
                       const std::map<std::string, MotifStats>& stats,
                       int motif_length) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot create logo matrix file " << filename << std::endl;
        return;
    }
    
    // Calculate position-specific nucleotide frequencies weighted by occurrence count
    std::vector<std::map<char, double>> position_freq(motif_length);
    
    size_t total_count = 0;
    for (const auto& [seq, stat] : stats) {
        size_t count = stat.count_total;
        total_count += count;
        
        // Add this sequence's contribution to each position
        for (size_t i = 0; i < seq.length() && i < static_cast<size_t>(motif_length); ++i) {
            position_freq[i][seq[i]] += count;
        }
    }
    
    // Write header
    out << "Position\tA\tC\tG\tT\tEntropy\tInfo_Content\n";
    
    // Calculate and write frequencies for each position
    for (int pos = 0; pos < motif_length; ++pos) {
        double freq_A = position_freq[pos]['A'] / total_count;
        double freq_C = position_freq[pos]['C'] / total_count;
        double freq_G = position_freq[pos]['G'] / total_count;
        double freq_T = position_freq[pos]['T'] / total_count;
        
        // Calculate Shannon entropy
        double entropy = 0;
        for (double f : {freq_A, freq_C, freq_G, freq_T}) {
            if (f > 0) {
                entropy -= f * log2(f);
            }
        }
        
        // Information content (max entropy = 2 for 4 nucleotides)
        double info_content = 2 - entropy;
        
        out << (pos + 1) << "\t"
            << std::fixed << std::setprecision(4)
            << freq_A << "\t"
            << freq_C << "\t"
            << freq_G << "\t"
            << freq_T << "\t"
            << entropy << "\t"
            << info_content << "\n";
    }
    
    out.close();
    std::cout << "Wrote consensus logo matrix to " << filename << "\n";
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
    out << "motif\talignment\tedit_distance\tcigar\tcount_forward\tcount_reverse\tcount_total\tfreq_forward\tfreq_reverse\tfreq_total\n";
    
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
            << stat.alignment << "\t"
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
              << "  " << std::setw(20) << "Alignment"
              << "  " << std::setw(8) << "Distance"
              << "  " << std::setw(20) << "CIGAR"
              << "  " << std::setw(12) << "Total Count"
              << "  " << std::setw(12) << "Frequency (%)\n";
    std::cout << std::string(110, '-') << "\n";
    
    int count = 0;
    for (const auto& pair : sorted_stats) {
        if (count++ >= 10) break;
        
        size_t total = pair.second.count_forward + pair.second.count_reverse;
        double freq = (total_all > 0) ? (100.0 * total / total_all) : 0;
        
        std::cout << std::setw(20) << pair.first
                  << "  " << std::setw(20) << pair.second.alignment
                  << "  " << std::setw(8) << pair.second.edit_distance
                  << "  " << std::setw(20) << pair.second.cigar
                  << "  " << std::setw(12) << total
                  << "  " << std::setw(12) << std::fixed << std::setprecision(2) << freq << "\n";
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
    std::cerr << "  alignment      Visual alignment (. = match, letter = mismatch, - = gap)\n";
    std::cerr << "  edit_distance  Edit distance from reference (includes indels)\n";
    std::cerr << "  cigar          CIGAR string (M=match, X=mismatch, I=insertion, D=deletion)\n";
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
    
    // Process IUPAC codes in the reference motif
    std::cout << "Using IUPAC-aware matching for pattern: " << motif << "\n";
    std::vector<std::pair<std::string, bool>> patterns_with_strand = expand_iupac_with_rc(motif);
    std::cout << "Will search for forward pattern and its reverse complement\n";
    
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
    total_bases = genome_size;  // Only forward strand now
    bases_processed = 0;
    
    // Process sequences in parallel
    std::cout << "Searching for motif variants...\n";
    
    // Start progress bar thread
    std::thread progress_thread(show_progress);
    
    // Divide work among threads more efficiently
    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    
    // Create work queue
    std::atomic<size_t> work_index(0);
    
    for (int t = 0; t < num_threads; ++t) {
        threads.emplace_back([&sequences, &work_index, &patterns_with_strand, max_distance]() {
            while (true) {
                size_t idx = work_index.fetch_add(1);
                if (idx >= sequences.size()) break;
                
                // Process single sequence
                std::vector<std::pair<std::string, std::string>> single_seq;
                single_seq.push_back(sequences[idx]);
                process_fasta_chunk(single_seq, patterns_with_strand, max_distance);
            }
        });
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
    std::cout << "Found " << global_occurrences.size() << " total occurrences\n";
    
    // Write results
    std::cout << "Writing results...\n";
    
    // Write TSV summary
    write_tsv(output_file, global_stats);
    
    // Write BED file with positions
    std::string bed_file = output_file.substr(0, output_file.rfind('.')) + ".bed";
    write_bed(bed_file, global_occurrences);
    
    // Write logo matrix
    std::string logo_file = output_file.substr(0, output_file.rfind('.')) + "_logo.tsv";
    write_logo_matrix(logo_file, global_stats, motif.length());
    
    std::cout << "\nDone!\n";
    
    return 0;
}