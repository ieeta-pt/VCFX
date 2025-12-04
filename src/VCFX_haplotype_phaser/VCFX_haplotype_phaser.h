#ifndef VCFX_HAPLOTYPE_PHASER_H
#define VCFX_HAPLOTYPE_PHASER_H

#include <deque>
#include <iostream>
#include <string>
#include <vector>

// A small struct to store a variant's key data: chromosome, position, plus the "allele sum" genotype
struct VariantData {
    std::string chrom;
    int pos;
    int index;  // Original variant index for output
    std::vector<int> genotype; // one element per sample (the sum of alleles, or -1 if missing)
};

// Stores the result of LD calculation
struct LDResult {
    double r;  // Correlation coefficient
    double r2; // Squared correlation coefficient
};

class VCFXHaplotypePhaser {
  public:
    // main runner
    int run(int argc, char *argv[]);

  private:
    // prints usage
    void displayHelp();

    // Main function that does phasing (default mode - loads all variants)
    void phaseHaplotypes(std::istream &in, std::ostream &out, double ldThreshold);

    // Streaming mode - uses sliding window, outputs blocks incrementally
    void phaseHaplotypesStreaming(std::istream &in, std::ostream &out, double ldThreshold, size_t windowSize);

    // Groups variants into haplotype blocks by naive r^2 threshold
    std::vector<std::vector<int>> groupVariants(const std::vector<VariantData> &variants, double ldThreshold);

    // calculates r^2 between two variants
    LDResult calculateLD(const VariantData &v1, const VariantData &v2);

    // Configuration
    bool streamingMode = false;
    size_t windowSize = 1000;  // Default window size for streaming mode
};

#endif // VCFX_HAPLOTYPE_PHASER_H
