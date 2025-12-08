#ifndef VCFX_DOSAGE_CALCULATOR_H
#define VCFX_DOSAGE_CALCULATOR_H

#include <iostream>
#include <string>
#include <vector>

// VCFX_dosage_calculator: Header file for Genotype Dosage Calculation tool
class VCFXDosageCalculator {
  public:
    // Entry point for the tool
    int run(int argc, char *argv[]);

  private:
    // Displays the help message
    void displayHelp();

    // Calculates genotype dosage from VCF input and writes output (stdin mode)
    void calculateDosage(std::istream &in, std::ostream &out);

    // Process file using memory-mapped I/O (optimized for large files)
    bool processFileMmap(const char* filename, std::ostream &out);

    // Helper function to split a string by a delimiter
    std::vector<std::string> split(const std::string &str, char delimiter);

    // Quiet mode flag (suppress warnings)
    bool quietMode = false;
};

#endif // VCFX_DOSAGE_CALCULATOR_H
