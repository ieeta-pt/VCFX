#ifndef VCFX_SORTER_H
#define VCFX_SORTER_H

#include <fstream>
#include <functional>
#include <memory>
#include <queue>
#include <string>
#include <vector>

// Compact record for sorting - stores only what's needed for comparison
struct SortKey {
    std::string chrom;
    int pos;
    size_t line_offset;  // Position in temp file for external sort
    size_t line_length;  // Length of the line

    // For in-memory sorting, store the full line
    std::string line;
};

// For external merge sort: a record from one of the temp files
struct MergeEntry {
    std::string chrom;
    int pos;
    std::string line;
    size_t file_index;

    // Comparator for min-heap (we want smallest first, so reverse logic)
    bool operator>(const MergeEntry& other) const;
};

class VCFXSorter {
  public:
    int run(int argc, char *argv[]);

    // Natural chromosome parsing helper - public so MergeEntry can use it
    static bool parseChromNat(const std::string &chrom, std::string &prefix, long &num, std::string &suffix);

  private:
    void displayHelp();

    // In-memory sort for small files
    void sortInMemory(std::istream &in, std::ostream &out);

    // External merge sort for large files
    void sortExternal(std::istream &in, std::ostream &out);

    // Write a sorted chunk to a temp file, returns the temp file path
    std::string writeChunk(std::vector<SortKey>& chunk, size_t chunk_num);

    // K-way merge of sorted temp files
    void mergeChunks(const std::vector<std::string>& chunk_files, std::ostream &out);

    // Parse chrom and pos from a VCF line
    static bool parseChromPos(const std::string& line, std::string& chrom, int& pos);

    // Comparison functions
    static bool lexCompare(const SortKey &a, const SortKey &b);
    static bool naturalCompare(const SortKey &a, const SortKey &b);

    // Settings
    bool naturalChromOrder = false;
    size_t chunkSizeBytes = 100 * 1024 * 1024;  // 100MB default chunk size
    std::string tempDir = "/tmp";

    // Store header lines
    std::vector<std::string> headerLines;

    // Comparison function pointer based on settings
    std::function<bool(const SortKey&, const SortKey&)> compareFunc;
};

#endif
