#ifndef VCFX_SORTER_H
#define VCFX_SORTER_H

#include <fstream>
#include <functional>
#include <memory>
#include <queue>
#include <string>
#include <unordered_map>
#include <vector>

// ============================================================================
// PERFORMANCE OPTIMIZATION: Pre-computed chromosome IDs
// ============================================================================
// Key insight: Chromosome comparison is called O(n log n) times during sort.
// By pre-computing a numeric ID once during parsing, we eliminate millions
// of string comparisons and allocations.
// ============================================================================

// Compact record for mmap-based sorting
struct CompactSortKey {
    int32_t chrom_id;     // Pre-computed chromosome ID for O(1) comparison
    int32_t pos;          // Genomic position
    size_t offset;        // Offset in mmap'd file
    uint32_t length;      // Line length (excluding newline)
    uint16_t chrom_offset; // Offset of CHROM field (always 0 for mmap)
    uint16_t chrom_len;    // Length of CHROM field (for lexicographic comparison)
};

// Original SortKey for stdin/external sort (backward compatible)
struct SortKey {
    std::string chrom;
    int32_t chrom_id;     // NEW: Pre-computed chromosome ID
    int pos;
    size_t line_offset;   // Position in temp file for external sort
    size_t line_length;   // Length of the line

    // For in-memory sorting, store the full line
    std::string line;
};

// For external merge sort: a record from one of the temp files
struct MergeEntry {
    std::string chrom;
    int32_t chrom_id;     // NEW: Pre-computed chromosome ID
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

    // OPTIMIZED: Convert chromosome string to numeric ID for fast comparison
    // This is called once per variant during parsing, not during sort comparisons
    static int32_t chromToId(const char* chrom, size_t len, bool naturalOrder);

  private:
    void displayHelp();

    // In-memory sort for small files
    void sortInMemory(std::istream &in, std::ostream &out);

    // External merge sort for large files
    void sortExternal(std::istream &in, std::ostream &out);

    // OPTIMIZED: Memory-mapped file sorting (fastest path)
    bool sortFileMmap(const char* filename, std::ostream &out);

    // Write a sorted chunk to a temp file, returns the temp file path
    std::string writeChunk(std::vector<SortKey>& chunk, size_t chunk_num);

    // K-way merge of sorted temp files
    void mergeChunks(const std::vector<std::string>& chunk_files, std::ostream &out);

    // Parse chrom and pos from a VCF line
    static bool parseChromPos(const std::string& line, std::string& chrom, int& pos);

    // OPTIMIZED: Parse chrom/pos directly from raw pointer (no string allocation)
    static bool parseChromPosFast(const char* line, size_t lineLen,
                                   const char*& chromStart, size_t& chromLen,
                                   int& pos);

    // Comparison functions
    static bool lexCompare(const SortKey &a, const SortKey &b);
    static bool naturalCompare(const SortKey &a, const SortKey &b);

    // OPTIMIZED: Comparison using pre-computed IDs (O(1) integer comparison)
    static bool compareById(const SortKey &a, const SortKey &b);

    // Settings
    bool naturalChromOrder = false;
    size_t chunkSizeBytes = 100 * 1024 * 1024;  // 100MB default chunk size
    std::string tempDir = "/tmp";

    // Store header lines
    std::vector<std::string> headerLines;

    // Comparison function pointer based on settings
    std::function<bool(const SortKey&, const SortKey&)> compareFunc;

    // Dynamic chromosome ID mapping for unknown chromosomes
    std::unordered_map<std::string, int32_t> dynamicChromIds;
    int32_t nextDynamicId = 1000;  // Start after all known chromosomes
};

#endif
