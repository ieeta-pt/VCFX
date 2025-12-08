#ifndef VCFX_AF_SUBSETTER_H
#define VCFX_AF_SUBSETTER_H

#include <iostream>
#include <string>
#include <vector>

// VCFXAfSubsetter: Header file for Alternate Allele Frequency Subsetter Tool
class VCFXAfSubsetter {
  public:
    // Entry point for the tool
    int run(int argc, char *argv[]);

  private:
    // Displays the help message
    void displayHelp();

    // Subsets VCF input based on alternate allele frequency range
    void subsetByAlleleFrequency(std::istream &in, std::ostream &out, double minAF, double maxAF);

    // Process file using memory-mapped I/O (optimized for large files)
    bool processFileMmap(const char* filename, std::ostream &out, double minAF, double maxAF);

    // Parses the AF values from the INFO field (handles multi-allelic AF as comma-delimited)
    bool parseAF(const std::string &infoField, std::vector<double> &afValues);

    // Parse AF from raw pointer (mmap mode)
    bool parseAFRaw(const char* infoStart, size_t infoLen, std::vector<double> &afValues);

    // Quiet mode flag
    bool quietMode = false;
};

#endif // VCFX_AF_SUBSETTER_H
