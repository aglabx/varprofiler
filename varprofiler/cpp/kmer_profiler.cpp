#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <stdexcept>
#include <thread>
#include <map>
#include <mutex>
#include <future>
#include <queue>
#include <condition_variable>

// Structure for storing the result of one window
struct BedResult {
    std::string chrom_name;
    size_t start;
    size_t end;
    size_t count;
};

// --- Functions for working with k-mers in 2-bit encoding ---

// Convert nucleotide to 2-bit representation
inline uint64_t nucl_to_2bit(char c) {
    // A=00, C=01, G=10, T=11
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4; // Marker for 'N' and others
    }
}

// Get reverse complement k-mer (in 2-bit encoding)
uint64_t reverse_complement_2bit(uint64_t kmer, int k_mer_length) {
    uint64_t rc = 0;
    uint64_t temp_kmer = kmer;
    for (int i = 0; i < k_mer_length; ++i) {
        rc <<= 2;
        // (x ^ 3) inverts 2 bits: 00<>11, 01<>10
        rc |= (temp_kmer & 3) ^ 3;
        temp_kmer >>= 2;
    }
    return rc;
}

// --- Main logic ---

/**
 * @brief Processes one chromosome using an efficient sliding algorithm.
 * @param chrom_name Chromosome name.
 * @param seq Nucleotide sequence of the chromosome.
 * @param k_mer_length K-mer length.
 * @param window_size Window size.
 * @param step_size Step size.
 * @return Vector with results for each window on this chromosome.
 */
std::vector<BedResult> process_chromosome_sliding(const std::string& chrom_name, const std::string& seq, int k_mer_length, int window_size, int step_size) {
    std::vector<BedResult> results;
    if (seq.length() < static_cast<size_t>(window_size)) {
        std::cerr << "Warning: Chromosome " << chrom_name << " is too short, skipping." << std::endl;
        return results;
    }

    std::cout << "Starting processing chromosome: " << chrom_name << " in thread " << std::this_thread::get_id() << std::endl;
    
    // Mask is now calculated here, based on the passed k_mer_length
    const uint64_t K_MER_MASK = (1ULL << (k_mer_length * 2)) - 1;

    std::unordered_map<uint64_t, int> kmer_counts;
    
    // --- Initialize first WINDOW ---
    for (size_t i = 0; i < static_cast<size_t>(window_size - k_mer_length + 1); ++i) {
        std::string kmer_str = seq.substr(i, k_mer_length);
        if (kmer_str.find('N') != std::string::npos || kmer_str.find('n') != std::string::npos) continue;
        
        uint64_t current_kmer = 0;
        for(char c : kmer_str) current_kmer = (current_kmer << 2) | nucl_to_2bit(c);
        current_kmer &= K_MER_MASK;
        
        uint64_t rc_kmer = reverse_complement_2bit(current_kmer, k_mer_length);
        kmer_counts[std::min(current_kmer, rc_kmer)]++;
    }
    results.push_back({chrom_name, 0, (size_t)window_size, kmer_counts.size()});

    // --- Process remaining windows ---
    for (size_t start = static_cast<size_t>(step_size); start + static_cast<size_t>(window_size) <= seq.length(); start += static_cast<size_t>(step_size)) {
        size_t end = start + static_cast<size_t>(window_size);
        
        // Clear and recalculate k-mers for this window
        // This is simpler and avoids accumulation errors
        kmer_counts.clear();
        
        for (size_t i = start; i < end - static_cast<size_t>(k_mer_length) + 1; ++i) {
            if (i + static_cast<size_t>(k_mer_length) > seq.length()) break;
            std::string kmer_str = seq.substr(i, k_mer_length);
            if (kmer_str.find('N') != std::string::npos || kmer_str.find('n') != std::string::npos) continue;
            
            uint64_t current_kmer = 0;
            for(char c : kmer_str) current_kmer = (current_kmer << 2) | nucl_to_2bit(c);
            current_kmer &= K_MER_MASK;
            
            uint64_t rc_kmer = reverse_complement_2bit(current_kmer, k_mer_length);
            kmer_counts[std::min(current_kmer, rc_kmer)]++;
        }
        
        results.push_back({chrom_name, start, end, kmer_counts.size()});
    }
    
    std::cout << "Finished processing chromosome: " << chrom_name << std::endl;
    return results;
}


int main(int argc, char* argv[]) {
    if (argc != 6 && argc != 7) {
        std::cerr << "Usage: " << argv[0] << " <input.fasta> <output.bed> <kmer_len> <window_size> <step_size> [threads]" << std::endl;
        std::cerr << "Example: " << argv[0] << " genome.fa out.bed 23 100000 25000" << std::endl;
        std::cerr << "Example with threads: " << argv[0] << " genome.fa out.bed 23 100000 25000 8" << std::endl;
        return 1;
    }

    std::string fasta_path = argv[1];
    std::string bed_path = argv[2];
    int k_mer_length = 0;
    int window_size = 0;
    int step_size = 0;
    int num_threads = std::thread::hardware_concurrency(); // Default to number of CPU cores

    try {
        k_mer_length = std::stoi(argv[3]);
        window_size = std::stoi(argv[4]);
        step_size = std::stoi(argv[5]);
        if (argc == 7) {
            num_threads = std::stoi(argv[6]);
            if (num_threads < 1) {
                std::cerr << "Error: Number of threads must be at least 1." << std::endl;
                return 1;
            }
        }
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: Invalid format of numeric arguments." << std::endl;
        return 1;
    }

    std::cout << "Using " << num_threads << " threads for processing." << std::endl;

    if (k_mer_length > 31 || k_mer_length < 1) {
        std::cerr << "Error: K-mer length must be in range [1, 31] for 64-bit representation." << std::endl;
        return 1;
    }

    // --- Read entire FASTA file into memory ---
    std::cout << "Reading FASTA file into memory..." << std::endl;
    std::map<std::string, std::string> genome;
    std::ifstream fasta_file(fasta_path);
    if (!fasta_file) {
        std::cerr << "Critical error: Failed to open FASTA file: " << fasta_path << std::endl;
        return 1;
    }
    std::string current_chrom_name, current_chrom_seq;
    current_chrom_seq.reserve(300000000);
    std::string line;
    while (std::getline(fasta_file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!current_chrom_name.empty()) {
                genome[current_chrom_name] = std::move(current_chrom_seq);
            }
            size_t first_space = line.find(' ');
            current_chrom_name = line.substr(1, first_space != std::string::npos ? first_space - 1 : std::string::npos);
            current_chrom_seq.clear();
        } else {
            line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
            current_chrom_seq.append(line);
        }
    }
    if (!current_chrom_name.empty()) {
        genome[current_chrom_name] = std::move(current_chrom_seq);
    }
    std::cout << "Reading completed. Found " << genome.size() << " chromosomes." << std::endl;

    // --- Launch parallel processing with controlled thread pool ---
    std::vector<std::future<std::vector<BedResult>>> futures;
    std::vector<std::pair<std::string, std::string*>> chrom_queue;
    
    // Prepare chromosome queue
    for (auto& [name, seq] : genome) {
        chrom_queue.push_back({name, &seq});
    }
    
    // Process chromosomes in batches based on thread count
    size_t batch_size = std::min(static_cast<size_t>(num_threads), chrom_queue.size());
    
    for (size_t i = 0; i < chrom_queue.size(); i += batch_size) {
        // Launch a batch of threads
        size_t current_batch_size = std::min(batch_size, chrom_queue.size() - i);
        for (size_t j = 0; j < current_batch_size; ++j) {
            auto& [name, seq_ptr] = chrom_queue[i + j];
            futures.push_back(std::async(std::launch::async, process_chromosome_sliding, 
                                        name, std::ref(*seq_ptr), k_mer_length, window_size, step_size));
        }
        
        // If not the last batch, wait for current batch to complete before starting next
        if (i + batch_size < chrom_queue.size()) {
            for (size_t j = i; j < i + current_batch_size; ++j) {
                futures[j].wait();
            }
        }
    }

    // --- Collect and write results ---
    std::cout << "Collecting results from threads..." << std::endl;
    std::ofstream bed_file(bed_path);
    if (!bed_file) {
        std::cerr << "Critical error: Failed to create BED file for writing: " << bed_path << std::endl;
        return 1;
    }

    // Write results maintaining chromosome order
    for (auto& fut : futures) {
        std::vector<BedResult> chrom_results = fut.get(); // Wait for thread completion and get result
        for (const auto& res : chrom_results) {
            bed_file << res.chrom_name << "\t" << res.start << "\t" << res.end << "\t" << res.count << "\n";
        }
    }

    std::cout << "Processing fully completed. Results saved to " << bed_path << std::endl;

    return 0;
}
