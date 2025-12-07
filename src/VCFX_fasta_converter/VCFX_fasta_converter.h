#ifndef VCFX_FASTA_CONVERTER_H
#define VCFX_FASTA_CONVERTER_H

#include <iostream>
#include <string>
#include <vector>

// VCFXFastaConverter: Tool for converting a "variant-only" VCF into per-sample FASTA sequences
class VCFXFastaConverter {
  public:
    // Entry point for the tool
    int run(int argc, char *argv[]);

  private:
    // Displays the help message
    void displayHelp();

    // Converts VCF input to FASTA format (stdin mode - fallback)
    void convertVCFtoFasta(std::istream &in, std::ostream &out);

    // Converts VCF input to FASTA format (streaming algorithm - fastest)
    // Uses per-sample temp files for O(1) random I/O - scales to 200+ GB files
    bool convertVCFtoFastaStreaming(const char *filename, std::ostream &out);

    // Compatibility wrapper (calls convertVCFtoFastaStreaming)
    bool convertVCFtoFastaMmap(const char *filename, std::ostream &out);

    // Quiet mode - suppress warnings
    bool quiet_ = false;
};

#endif // VCFX_FASTA_CONVERTER_H
