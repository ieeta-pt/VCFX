#ifndef VCFX_MERGER_H
#define VCFX_MERGER_H

#include <fstream>
#include <iostream>
#include <memory>
#include <queue>
#include <string>
#include <vector>

// Forward declaration for streaming merge
struct MergeFileEntry;

// VCFX_merger: Header file for VCF file merging tool
class VCFXMerger {
  public:
    // Entry point for the tool
    int run(int argc, char *argv[]);

    // Chromosome parsing for natural comparison (public for MergeFileEntry)
    static bool parseChromNat(const std::string &chrom, std::string &prefix, long &num, std::string &suffix);

  private:
    // Displays the help message
    void displayHelp();

    // Original method: loads all variants into memory, then sorts
    void mergeVCFInMemory(const std::vector<std::string> &inputFiles, std::ostream &out);

    // Streaming method: assumes sorted inputs, uses K-way merge with O(num_files) memory
    void mergeVCFStreaming(const std::vector<std::string> &inputFiles, std::ostream &out);

    // Parse CHROM and POS from a VCF line
    static bool parseChromPos(const std::string& line, std::string& chrom, long& pos);

    // Settings
    bool assumeSorted = false;
    bool naturalChromOrder = false;
};

// Entry for the priority queue during streaming merge
struct MergeFileEntry {
    std::string chrom;
    long pos = 0;
    std::string line;
    size_t file_index;

    // Comparator for min-heap (smallest first, so use > for priority_queue)
    bool operator>(const MergeFileEntry& other) const;
};

#endif // VCFX_MERGER_H
