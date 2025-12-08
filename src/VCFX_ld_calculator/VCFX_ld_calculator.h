#ifndef VCFX_LD_CALCULATOR_H
#define VCFX_LD_CALCULATOR_H

#include <deque>
#include <iostream>
#include <string>
#include <vector>
#include <cstdint>

// Holds minimal data for a single variant (used in matrix mode)
struct LDVariant {
    std::string chrom;
    int pos;
    std::string id;  // ID field for output
    std::vector<int> genotype; // 0 => homRef, 1 => het, 2 => homAlt, -1 => missing
};

class VCFXLDCalculator {
  public:
    int run(int argc, char *argv[]);

  private:
    void displayHelp();
    // Parse the region string "chr1:10000-20000"
    // If none provided, regionChrom will be empty => use all
    bool parseRegion(const std::string &regionStr, std::string &regionChrom, int &regionStart, int &regionEnd);

    // The main logic: read VCF, store genotype codes for in-range variants, compute r^2, output
    // Default mode: MxM matrix (backward compatible)
    void computeLD(std::istream &in, std::ostream &out, const std::string &regionChrom, int regionStart, int regionEnd);

    // Streaming mode: sliding window, output pairs incrementally
    void computeLDStreaming(std::istream &in, std::ostream &out, const std::string &regionChrom,
                            int regionStart, int regionEnd, size_t windowSize, double threshold);

    // Optimized mmap-based versions
    void computeLDStreamingMmap(const char* data, size_t size, int outFd, const std::string &regionChrom,
                                 int regionStart, int regionEnd, size_t windowSize, double threshold,
                                 int maxDistance);
    void computeLDMatrixMmap(const char* data, size_t size, int outFd, const std::string &regionChrom,
                              int regionStart, int regionEnd, int numThreads);

    // parse a single genotype string => code
    int parseGenotype(const std::string &s);

    // compute r^2 for two variant's genotype arrays
    double computeRsq(const std::vector<int> &g1, const std::vector<int> &g2);

    // Configuration
    bool streamingMode = true;   // Default: streaming mode for performance
    bool matrixMode = false;     // Explicit matrix mode (backward compat)
    size_t windowSize = 1000;    // Default window size
    double ldThreshold = 0.0;    // Default: output all pairs
    int numThreads = 0;          // 0 = auto-detect
    int maxDistance = 0;         // 0 = no distance limit
    bool quiet = false;          // Suppress informational messages
    std::string inputFile;       // Empty = stdin
};

#endif // VCFX_LD_CALCULATOR_H
