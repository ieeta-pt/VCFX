#ifndef VCFX_HAPLOTYPE_EXTRACTOR_H
#define VCFX_HAPLOTYPE_EXTRACTOR_H

#include <iostream>
#include <string>
#include <vector>

// Function to display help message
void printHelp();

// Structure to represent a haplotype block
struct HaplotypeBlock {
    std::string chrom;
    int start;
    int end;
    std::vector<std::string> haplotypes; // One haplotype "string" per sample
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

    // Set the maximum distance for grouping consecutive variants in a block
    void setBlockDistanceThreshold(int dist) { blockDistanceThreshold = dist; }

    // If true, we do a minimal consistency check across variants
    void setCheckPhaseConsistency(bool b) { checkPhaseConsistency = b; }

    // Enable or disable debug messages
    void setDebug(bool b) { debugMode = b; }

    // Enable streaming mode
    void setStreamingMode(bool b) { streamingMode = b; }

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

    // Parses the #CHROM line to extract sample names
    bool parseHeader(const std::string &headerLine);

    // Splits a string by a delimiter
    std::vector<std::string> splitString(const std::string &str, char delimiter);

    // Processes one VCF data line => update or start a haplotype block
    // Returns false if variant not fully processed
    bool processVariant(const std::vector<std::string> &fields, std::vector<HaplotypeBlock> &haplotypeBlocks);

    // Streaming version: process variant and output completed blocks immediately
    // Returns the completed block if any, or nullptr
    bool processVariantStreaming(const std::vector<std::string> &fields,
                                  HaplotypeBlock &currentBlock,
                                  bool &hasCurrentBlock,
                                  std::ostream &out,
                                  bool headerWritten);

    // For each sample's genotype, ensures it is phased. If any unphased => return false
    bool areAllSamplesPhased(const std::vector<std::string> &genotypes);

    // Minimal check that new variant's genotypes are "consistent" with the existing block
    bool phaseIsConsistent(const HaplotypeBlock &block, const std::vector<std::string> &newGenotypes);

    // Actually merges the new variant's genotypes into the last block or starts a new one
    void updateBlocks(std::vector<HaplotypeBlock> &haplotypeBlocks, const std::string &chrom, int pos,
                      const std::vector<std::string> &genotypes);

    // Output a single block to stream
    void outputBlock(std::ostream &out, const HaplotypeBlock &block);

    // Output header line
    void outputHeader(std::ostream &out);
};

#endif // VCFX_HAPLOTYPE_EXTRACTOR_H
