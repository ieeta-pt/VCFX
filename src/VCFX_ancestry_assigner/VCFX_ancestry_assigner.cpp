#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// ---------------------------------------------------------
// A helper struct to store command-line options if needed
// ---------------------------------------------------------
struct AncestryOptions {
    bool showHelp = false;
    std::string freqFile;
};

// ---------------------------------------------------------
// Class: VCFXAncestryAssigner
// ---------------------------------------------------------
class VCFXAncestryAssigner {
  public:
    // High-level entry point
    int run(int argc, char *argv[]);

  private:
    // Command-line usage
    void displayHelp();

    // Parse freq file
    bool loadAncestralFrequencies(std::istream &in);

    // Parse one frequency line
    // Format (tab-separated): CHROM  POS  REF  ALT  POP1_FREQ  POP2_FREQ ...
    bool parseFrequencyLine(const std::string &line);

    // Assign ancestry from VCF
    void assignAncestry(std::istream &vcfIn, std::ostream &out);

  private:
    // Populations in order
    std::vector<std::string> populations;

    // Frequencies by key: "chrom:pos:ref:alt" => (pop => freq)
    std::unordered_map<std::string, std::unordered_map<std::string, double>> variantFrequencies;
};

// ---------------------------------------------------------
// VCFXAncestryAssigner::run
// ---------------------------------------------------------
int VCFXAncestryAssigner::run(int argc, char *argv[]) {
    // 1. Parse arguments
    int opt;
    bool showHelp = false;
    std::string freqFile;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'}, {"assign-ancestry", required_argument, 0, 'a'}, {0, 0, 0, 0}};

    while ((opt = getopt_long(argc, argv, "ha:", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
            break;
        case 'a':
            freqFile = std::string(optarg);
            break;
        default:
            showHelp = true;
        }
    }

    if (showHelp || freqFile.empty()) {
        displayHelp();
        return (showHelp ? 0 : 1);
    }

    // 2. Open frequency file
    std::ifstream freqStream(freqFile);
    if (!freqStream.is_open()) {
        std::cerr << "Error: Unable to open frequency file: " << freqFile << "\n";
        return 1;
    }

    // 3. Load frequencies
    if (!loadAncestralFrequencies(freqStream)) {
        std::cerr << "Error: Failed to load ancestral frequencies.\n";
        return 1;
    }

    // 4. Assign ancestry based on VCF (read from stdin, write to stdout)
    assignAncestry(std::cin, std::cout);
    return 0;
}

// ---------------------------------------------------------
// Show usage
// ---------------------------------------------------------
void VCFXAncestryAssigner::displayHelp() {
    std::cout << "VCFX_ancestry_assigner: Assign samples to ancestral populations based on variant frequencies.\n\n"
              << "Usage:\n"
              << "  VCFX_ancestry_assigner --assign-ancestry <freq_file> < input.vcf > ancestry.txt\n\n"
              << "Options:\n"
              << "  -h, --help                 Show this help message and exit\n"
              << "  -a, --assign-ancestry FILE Ancestral frequency file\n\n"
              << "Frequency File Format:\n"
              << "  The first line must be a header like:\n"
              << "    CHROM  POS  REF  ALT  POP1  POP2  ...\n"
              << "  Each subsequent line must have the same columns. For example:\n"
              << "    1   10000   A   C   0.10  0.20\n\n"
              << "Example:\n"
              << "  VCFX_ancestry_assigner --assign-ancestry ancestral_freq.tsv < input.vcf > ancestry_out.txt\n\n";
}

// ---------------------------------------------------------
// parseFrequencyLine
// Expects: CHROM  POS  REF  ALT  pop1Freq  pop2Freq ...
// ---------------------------------------------------------
bool VCFXAncestryAssigner::parseFrequencyLine(const std::string &line) {
    std::vector<std::string> fields;
    vcfx::split_tabs(line, fields);

    // Must have at least CHROM, POS, REF, ALT, plus the populations
    if (fields.size() < 4 + populations.size()) {
        return false;
    }
    const std::string &chrom = fields[0];
    int pos = 0;
    try {
        pos = std::stoi(fields[1]);
    } catch (...) {
        return false;
    }
    const std::string &ref = fields[2];
    const std::string &alt = fields[3];

    // Build freq map for this variant
    std::unordered_map<std::string, double> freqMap;
    for (size_t i = 0; i < populations.size(); ++i) {
        double freqVal = 0.0;
        try {
            freqVal = std::stod(fields[4 + i]);
        } catch (...) {
            freqVal = 0.0; // default if missing
        }
        freqMap[populations[i]] = freqVal;
    }

    // Key = CHROM:POS:REF:ALT
    const std::string key = chrom + ":" + std::to_string(pos) + ":" + ref + ":" + alt;
    variantFrequencies[key] = freqMap;
    return true;
}

// ---------------------------------------------------------
// loadAncestralFrequencies
// First line is header with columns:
//   CHROM  POS  REF  ALT  pop1  pop2 ...
// ---------------------------------------------------------
bool VCFXAncestryAssigner::loadAncestralFrequencies(std::istream &in) {
    std::string line;
    if (!std::getline(in, line)) {
        std::cerr << "Error: Frequency file is empty.\n";
        return false;
    }

    // Parse the header
    {
        std::vector<std::string> headers;
        vcfx::split_tabs(line, headers);
        // We need at least 5 columns: CHROM, POS, REF, ALT, plus 1 pop
        if (headers.size() < 5) {
            std::cerr << "Error: Frequency header must have at least 5 columns.\n";
            return false;
        }
        // The populations start at column index 4
        for (size_t i = 4; i < headers.size(); ++i) {
            populations.push_back(headers[i]);
        }
    }

    // Parse each subsequent line
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (!parseFrequencyLine(line)) {
            std::cerr << "Warning: Skipping invalid frequency line:\n" << line << "\n";
        }
    }
    return true;
}

// ---------------------------------------------------------
// assignAncestry
// Reads VCF from vcfIn, writes "Sample <tab> AssignedPopulation" to out
// ---------------------------------------------------------
void VCFXAncestryAssigner::assignAncestry(std::istream &vcfIn, std::ostream &out) {
    std::string line;
    bool haveHeader = false;
    int chrIndex = -1, posIndex = -1, refIndex = -1, altIndex = -1;
    std::vector<std::string> sampleNames;

    // For each sample: sampleScores[sample][population] => log likelihood
    std::unordered_map<std::string, std::unordered_map<std::string, double>> sampleScores;
    // For each sample: number of variants used in scoring
    std::unordered_map<std::string, int> sampleVariantCounts;

    // Reusable buffer for parsing
    std::vector<std::string> fields;
    fields.reserve(16);

    while (std::getline(vcfIn, line)) {
        if (line.empty()) {
            continue;
        }

        // VCF headers
        if (line[0] == '#') {
            // If #CHROM line, parse for sample columns
            if (line.rfind("#CHROM", 0) == 0) {
                haveHeader = true;
                std::vector<std::string> headers;
                vcfx::split_tabs(line, headers);
                // Identify columns
                for (size_t i = 0; i < headers.size(); ++i) {
                    if (headers[i] == "#CHROM" || headers[i] == "CHROM")
                        chrIndex = (int)i;
                    else if (headers[i] == "POS")
                        posIndex = (int)i;
                    else if (headers[i] == "REF")
                        refIndex = (int)i;
                    else if (headers[i] == "ALT")
                        altIndex = (int)i;
                }
                if (chrIndex < 0 || posIndex < 0 || refIndex < 0 || altIndex < 0) {
                    std::cerr << "Error: #CHROM header missing required columns CHROM POS REF ALT.\n";
                    return;
                }
                // Sample columns start at index 9
                for (size_t i = 9; i < headers.size(); ++i) {
                    sampleNames.push_back(headers[i]);
                    // Initialize scores for this sample
                    for (const auto &pop : populations) {
                        sampleScores[headers[i]][pop] = 0.0;
                    }
                    sampleVariantCounts[headers[i]] = 0;
                }
            }
            continue;
        }

        if (!haveHeader) {
            std::cerr << "Error: VCF header not found.\n";
            return;
        }

        // Parse VCF line
        fields.clear();
        vcfx::split_tabs(line, fields);

        if (fields.size() < 9) {
            std::cerr << "Warning: Skipping malformed VCF line:\n" << line << "\n";
            continue;
        }

        const std::string &chrom = fields[chrIndex];
        const std::string &pos = fields[posIndex];
        const std::string &ref = fields[refIndex];
        const std::string &alt = fields[altIndex];

        // Build key for frequency lookup
        const std::string key = chrom + ":" + pos + ":" + ref + ":" + alt;
        auto freqIt = variantFrequencies.find(key);
        if (freqIt == variantFrequencies.end()) {
            continue; // Skip variants not in frequency file
        }

        // Parse FORMAT field
        std::stringstream formatSS(fields[8]);
        std::vector<std::string> formatFields;
        std::string formatField;
        while (std::getline(formatSS, formatField, ':')) {
            formatFields.push_back(formatField);
        }

        // Find GT index in FORMAT
        int gtIndex = -1;
        for (size_t i = 0; i < formatFields.size(); ++i) {
            if (formatFields[i] == "GT") {
                gtIndex = (int)i;
                break;
            }
        }
        if (gtIndex < 0) {
            continue; // Skip if no GT field
        }

        // Process each sample
        for (size_t i = 0; i < sampleNames.size(); ++i) {
            const std::string &sample = sampleNames[i];
            if (i + 9 >= fields.size()) {
                continue;
            }

            // Parse genotype
            std::stringstream sampleSS(fields[i + 9]);
            std::vector<std::string> sampleFields;
            std::string sampleField;
            while (std::getline(sampleSS, sampleField, ':')) {
                sampleFields.push_back(sampleField);
            }

            if (gtIndex >= (int)sampleFields.size()) {
                continue; // Skip if GT field missing
            }

            const std::string &gt = sampleFields[gtIndex];
            if (gt == "./." || gt == ".|.") {
                continue; // Skip missing genotypes
            }

            // Parse alleles
            std::vector<std::string> alleles;
            std::string gtCopy = gt;
            // Replace '|' with '/' for consistency
            for (char &c : gtCopy) {
                if (c == '|')
                    c = '/';
            }
            std::stringstream gtSS(gtCopy);
            std::string allele;
            while (std::getline(gtSS, allele, '/')) {
                if (allele != ".") {
                    alleles.push_back(allele);
                }
            }
            if (alleles.size() != 2) {
                continue; // Skip invalid genotypes
            }

            // Calculate scores for each population
            for (const auto &pop : populations) {
                double altFreq = freqIt->second.at(pop);

                // Ensure frequencies are not too close to 0 or 1
                altFreq = std::max(0.001, std::min(0.999, altFreq));
                double refFreq = 1.0 - altFreq;

                // Calculate genotype probability under Hardy-Weinberg equilibrium
                double p = 0.0;
                if (gt == "0/0" || gt == "0|0") {
                    p = refFreq * refFreq;
                } else if (gt == "0/1" || gt == "1/0" || gt == "0|1" || gt == "1|0") {
                    p = 2.0 * refFreq * altFreq;
                } else if (gt == "1/1" || gt == "1|1") {
                    p = altFreq * altFreq;
                }

                // Add log-likelihood to score
                double logLikelihood = -std::log(p + 1e-12);
                sampleScores[sample][pop] += logLikelihood;
            }
            sampleVariantCounts[sample]++;
        }
    }

    // Output results
    for (const auto &sample : sampleNames) {
        if (sampleVariantCounts[sample] == 0) {
            out << sample << "\tUNKNOWN\n";
            continue;
        }

        // Find population with highest likelihood (lowest negative log-likelihood)
        std::string bestPop = populations[0];
        double bestScore = sampleScores[sample][populations[0]];

        std::cerr << "Final scores for " << sample << ":\n";
        for (const auto &pop : populations) {
            std::cerr << "  " << pop << ": " << sampleScores[sample][pop] << "\n";
            if (sampleScores[sample][pop] < bestScore) {
                bestScore = sampleScores[sample][pop];
                bestPop = pop;
            }
        }

        out << sample << "\t" << bestPop << "\n";
    }
}

// ---------------------------------------------------------
// main() - just instantiate and run
// ---------------------------------------------------------
static void show_help() {
    VCFXAncestryAssigner obj;
    char arg0[] = "VCFX_ancestry_assigner";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_ancestry_assigner", show_help))
        return 0;
    VCFXAncestryAssigner assigner;
    return assigner.run(argc, argv);
}
