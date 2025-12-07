#ifndef VCFX_DIFF_TOOL_H
#define VCFX_DIFF_TOOL_H

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

// VCFXDiffTool: Header file for VCF Diff Tool
// Supports two modes:
// 1. Default: Load both files into hash sets (works with unsorted files)
// 2. Streaming (--assume-sorted): Two-pointer merge diff with O(1) memory
//
// Performance: Uses memory-mapped I/O with SIMD-accelerated parsing
class VCFXDiffTool {
  public:
    // Entry point for the tool
    int run(int argc, char *argv[]);

  private:
    // Displays the help message
    void displayHelp();

    // === OPTIMIZED mmap-based implementations ===
    // These are used by default for file input
    // Returns false if files cannot be opened
    bool diffInMemoryMmap(const std::string &file1Path, const std::string &file2Path);
    bool diffStreamingMmap(const std::string &file1Path, const std::string &file2Path);

    // === Fallback implementations (kept for compatibility) ===
    // Loads variants from a VCF file into a set
    bool loadVariants(const std::string &filePath, std::unordered_set<std::string> &variants);

    // Compare using hash sets (original algorithm)
    void diffInMemory(const std::string &file1Path, const std::string &file2Path);

    // Compare sorted files using O(1) memory
    void diffStreaming(const std::string &file1Path, const std::string &file2Path);

    // Generates a unique key for a variant based on chromosome, position, ref, and alt
    std::string generateVariantKey(const std::string &chrom, const std::string &pos, const std::string &ref,
                                   const std::string &alt);

    // Parse a VCF line and extract key components
    // Returns false if line is invalid or a header
    bool parseVCFLine(const std::string &line, std::string &chrom, long &pos,
                      std::string &ref, std::string &alt, std::string &key);

    // Compare two variant keys for sorting order
    // Returns <0 if a < b, 0 if a == b, >0 if a > b
    int compareKeys(const std::string &chromA, long posA, const std::string &refA, const std::string &altA,
                    const std::string &chromB, long posB, const std::string &refB, const std::string &altB);

    // Configuration
    bool assumeSorted = false;
    bool naturalChromOrder = false;
    bool quiet_ = false;
};

#endif // VCFX_DIFF_TOOL_H
