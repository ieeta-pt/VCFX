#include "vcfx_core.h" 
#include "vcfx_io.h"
// VCFX_allele_balance_calc.cpp

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// ---------------------------------------------------------------
// Header / Declarations
// ---------------------------------------------------------------
struct AlleleBalanceArguments {
    std::vector<std::string> samples;
};

// Forward declarations
void printHelp();
bool parseArguments(int argc, char *argv[], AlleleBalanceArguments &args);
bool calculateAlleleBalance(std::istream &in, std::ostream &out, const AlleleBalanceArguments &args);
double computeAlleleBalance(const std::string &genotype);

// ---------------------------------------------------------------
// Utility Functions
// ---------------------------------------------------------------
static std::vector<std::string> splitString(const std::string &str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string temp;
    while (std::getline(ss, temp, delimiter)) {
        tokens.push_back(temp);
    }
    return tokens;
}

// Print help message
void printHelp() {
    std::cout << "VCFX_allele_balance_calc\n"
              << "Usage: VCFX_allele_balance_calc [OPTIONS] < input.vcf > allele_balance.tsv\n\n"
              << "Options:\n"
              << "  --samples, -s \"Sample1 Sample2\"   Specify the sample names to calculate allele balance for.\n"
              << "  --help, -h                        Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Calculates the allele balance (ratio of reference to alternate alleles) for each sample.\n"
              << "  Allele balance is computed as (#RefAlleles / #AltAlleles), using the genotype field.\n"
              << "  This simple logic treats all non-zero alleles as 'alt' and 0 as 'ref',\n"
              << "  so multi-allelic sites are lumped into an overall alt count.\n\n"
              << "Examples:\n"
              << "  1) Calculate allele balance for SampleA and SampleB:\n"
              << "     ./VCFX_allele_balance_calc --samples \"SampleA SampleB\" < input.vcf > allele_balance.tsv\n\n"
              << "  2) Calculate allele balance for all samples:\n"
              << "     ./VCFX_allele_balance_calc < input.vcf > allele_balance_all.tsv\n\n";
}

// Parse command-line arguments
bool parseArguments(int argc, char *argv[], AlleleBalanceArguments &args) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--samples" || arg == "-s") && i + 1 < argc) {
            // Collect samples from next argument
            std::string samplesStr = argv[++i];
            args.samples = splitString(samplesStr, ' ');
            // Trim whitespace for each sample
            for (auto &sample : args.samples) {
                // left trim
                sample.erase(0, sample.find_first_not_of(" \t\n\r\f\v"));
                // right trim
                sample.erase(sample.find_last_not_of(" \t\n\r\f\v") + 1);
            }
        } else if (arg == "--help" || arg == "-h") {
            printHelp();
            std::exit(0);
        } else {
            std::cerr << "Warning: Unrecognized argument '" << arg << "'.\n";
        }
    }
    return true;
}

// Compute allele balance
// If genotype is "0/0", result= ref_count / alt_count => 2 / 0 => 0.0
// If genotype is "0/1", result= 1 / 1 => 1.0
// If genotype is "1/2", result= 0 / 2 => 0.0, etc.
double computeAlleleBalance(const std::string &genotype) {
    if (genotype.empty() || genotype == "." || genotype == "./." || genotype == ".|.") {
        return -1.0; // missing genotype
    }

    // Split by '/' or '|'. For a simple approach, we can replace '|' with '/' and split on '/':
    std::string g = genotype;
    for (auto &ch : g) {
        if (ch == '|')
            ch = '/';
    }
    auto alleles = splitString(g, '/');

    int ref_count = 0;
    int alt_count = 0;
    for (const auto &allele : alleles) {
        if (allele == "0") {
            ++ref_count;
        } else if (allele != "." && !allele.empty()) {
            ++alt_count;
        }
    }

    // If alt_count is 0 but we have some reference calls, ratio is 0.0
    if (alt_count == 0 && ref_count > 0)
        return 0.0;
    // If no meaningful data, return -1
    if (ref_count + alt_count == 0)
        return -1.0;

    return static_cast<double>(ref_count) / alt_count;
}

// Main function to read VCF, parse genotypes, and calculate allele balance
bool calculateAlleleBalance(std::istream &in, std::ostream &out, const AlleleBalanceArguments &args) {
    std::string line;
    bool headerFound = false;
    std::vector<std::string> headerFields;
    std::vector<int> sampleIndices; // columns (9-based) that we care about
    std::unordered_map<std::string, int> sampleMap;

    // Print exactly one TSV header
    out << "CHROM\tPOS\tID\tREF\tALT\tSample\tAllele_Balance\n";

    while (std::getline(in, line)) {
        // Skip empty lines
        if (line.empty()) {
            continue;
        }

        // Skip and do not print lines that start with '#'
        //   BUT parse #CHROM line to get sample columns
        if (line[0] == '#') {
            // If this is the #CHROM line, parse columns
            if (line.rfind("#CHROM", 0) == 0) {
                headerFields = splitString(line, '\t');
                // Build sample -> column index
                for (size_t i = 9; i < headerFields.size(); ++i) {
                    sampleMap[headerFields[i]] = static_cast<int>(i);
                }
                // If user specified samples, pick them out; otherwise take all
                if (!args.samples.empty()) {
                    for (const auto &sample : args.samples) {
                        auto it = sampleMap.find(sample);
                        if (it == sampleMap.end()) {
                            std::cerr << "Error: Sample '" << sample << "' not found in VCF header.\n";
                            return false;
                        }
                        sampleIndices.push_back(it->second);
                    }
                } else {
                    // all samples
                    for (size_t i = 9; i < headerFields.size(); ++i) {
                        sampleIndices.push_back(static_cast<int>(i));
                    }
                }
                headerFound = true;
            }
            continue; // do not output these lines
        }

        // We need a #CHROM header before data
        if (!headerFound) {
            std::cerr << "Error: VCF header (#CHROM) not found before records.\n";
            return false;
        }

        // Split the data line
        auto fields = splitString(line, '\t');
        if (fields.size() < 9) {
            std::cerr << "Warning: Skipping invalid VCF line with fewer than 9 fields.\n";
            continue;
        }

        // Basic VCF columns
        const std::string &chrom = fields[0];
        const std::string &pos = fields[1];
        const std::string &id = fields[2];
        const std::string &ref = fields[3];
        const std::string &alt = fields[4];
        // The genotype info starts at index 9

        // Calculate AB for each requested sample
        for (auto idx : sampleIndices) {
            if (static_cast<size_t>(idx) >= fields.size()) {
                std::cerr << "Warning: sample index " << idx << " out of range.\n";
                continue;
            }
            const std::string &genotypeStr = fields[idx];
            // genotype is typically the first colon-delimited field
            // e.g. genotypeStr = "0/1:..."
            // but we only need the actual genotype for computing ratio
            // split on ':'
            auto subFields = splitString(genotypeStr, ':');
            if (subFields.empty()) {
                continue;
            }
            std::string genotype = subFields[0]; // e.g. "0/1"

            double ab = computeAlleleBalance(genotype);
            std::string abStr = (ab < 0.0) ? "NA" : std::to_string(ab);

            // The sample name is from the VCF header
            std::string sampleName =
                (idx < static_cast<int>(headerFields.size())) ? headerFields[idx] : ("Sample_" + std::to_string(idx));

            out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\t" << sampleName << "\t"
                << abStr << "\n";
        }
    }
    return true;
}

// ---------------------------------------------------------------
// main()
// ---------------------------------------------------------------
static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_allele_balance_calc", show_help))
        return 0;
    AlleleBalanceArguments args;
    parseArguments(argc, argv, args);

    if (!args.samples.empty()) {
        std::cerr << "Info: Calculating allele balance for these samples: ";
        for (auto &s : args.samples) {
            std::cerr << s << " ";
        }
        std::cerr << "\n";
    } else {
        std::cerr << "Info: Calculating allele balance for ALL samples.\n";
    }

    bool success = calculateAlleleBalance(std::cin, std::cout, args);
    return success ? 0 : 1;
}
