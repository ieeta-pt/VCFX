#ifndef VCFX_HAPLOTYPE_EXTRACTOR_H
#define VCFX_HAPLOTYPE_EXTRACTOR_H

#include <iostream>
#include <string>
#include <string_view>
#include <vector>

// Forward declaration
class OutputBuffer;

// Function to display help message
void printHelp();

// Structure to represent a haplotype block (optimized with lastGenotypes cache)
struct HaplotypeBlock {
    std::string chrom;
    int start = 0;
    int end = 0;
    std::vector<std::string> haplotypes;       // One haplotype "string" per sample
    std::vector<std::string> lastGenotypes;    // Cache last GT per sample for O(1) phase check
};

// Class to handle haplotype extraction
class HaplotypeExtractor {
  public:
    HaplotypeExtractor() = default;
    ~HaplotypeExtractor() = default;

    // Runs the core logic to parse the VCF and write haplotype blocks
    // Default mode: accumulates all blocks, outputs at end (backward compatible)
    bool extractHaplotypes(std::istream &in, std::ostream &out);

    // Streaming mode: outputs blocks immediately when complete (O(1) memory per block)
    bool extractHaplotypesStreaming(std::istream &in, std::ostream &out);

    // Memory-mapped file processing (50-100x faster than stdin)
    bool extractHaplotypesMmap(const char *filename, std::ostream &out);
    bool extractHaplotypesMmapStreaming(const char *filename, std::ostream &out);

    // Set the maximum distance for grouping consecutive variants in a block
    void setBlockDistanceThreshold(int dist) { blockDistanceThreshold = dist; }

    // If true, we do a minimal consistency check across variants
    void setCheckPhaseConsistency(bool b) { checkPhaseConsistency = b; }

    // Enable or disable debug messages
    void setDebug(bool b) { debugMode = b; }

    // Enable streaming mode
    void setStreamingMode(bool b) { streamingMode = b; }

    // Enable or disable quiet mode (suppress warnings)
    void setQuiet(bool b) { quiet_ = b; }

  private:
    std::vector<std::string> sampleNames;
    size_t numSamples = 0;

    // The maximum allowed distance to remain in the same block
    int blockDistanceThreshold = 100000; // default 100 kb

    // If true, we do a simplistic cross-variant check for consistent phasing
    bool checkPhaseConsistency = false;

    // If true, print verbose debugging information
    bool debugMode = false;

    // If true, output blocks immediately when complete
    bool streamingMode = false;

    // If true, suppress warning messages
    bool quiet_ = false;

    // FORMAT field caching for performance
    std::string cachedFormat_;
    int cachedGTIndex_ = -1;

    // Reusable vector for genotype fields
    std::vector<std::string_view> genotypeFields_;

    // Parses the #CHROM line to extract sample names
    bool parseHeader(std::string_view headerLine);

    // Unified fast processing for both batch and streaming modes
    template <bool Streaming>
    bool processVariantFast(const std::vector<std::string_view> &fields,
                            std::vector<HaplotypeBlock> &haplotypeBlocks,
                            HaplotypeBlock &currentBlock,
                            bool &hasCurrentBlock,
                            OutputBuffer &out);

    // For each sample's genotype, ensures it is phased. If any unphased => return false
    bool areAllSamplesPhased(const std::vector<std::string_view> &genotypes);

    // Minimal check that new variant's genotypes are "consistent" with the existing block
    // Uses cached lastGenotypes for O(1) performance instead of O(n) rfind
    bool phaseIsConsistent(const HaplotypeBlock &block, const std::vector<std::string_view> &newGenotypes);

    // Output a single block to stream
    void outputBlock(OutputBuffer &out, const HaplotypeBlock &block);

    // Output header line
    void outputHeader(OutputBuffer &out);
};

#endif // VCFX_HAPLOTYPE_EXTRACTOR_H
