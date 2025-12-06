#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <limits>
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

// ===========================================================================
// OPTIMIZED: Zero-allocation parsing helpers
// ===========================================================================

// Find GT index in FORMAT field using raw pointer
static inline int findGTIndexFast(const char* format, size_t formatLen) {
    int index = 0;
    size_t pos = 0;
    size_t start = 0;

    while (pos <= formatLen) {
        if (pos == formatLen || format[pos] == ':') {
            size_t tokenLen = pos - start;
            if (tokenLen == 2 && format[start] == 'G' && format[start + 1] == 'T') {
                return index;
            }
            index++;
            start = pos + 1;
        }
        pos++;
    }
    return -1;
}

// Extract GT field from sample data using raw pointers
// Returns pointer to start and length, -1 on failure
static inline bool extractGTField(const char* sample, size_t sampleLen,
                                  int gtIndex, const char*& gtStart, size_t& gtLen) {
    if (gtIndex < 0) return false;

    int currentIndex = 0;
    size_t pos = 0;
    size_t fieldStart = 0;

    while (pos <= sampleLen) {
        if (pos == sampleLen || sample[pos] == ':') {
            if (currentIndex == gtIndex) {
                gtStart = sample + fieldStart;
                gtLen = pos - fieldStart;
                return gtLen > 0;
            }
            currentIndex++;
            fieldStart = pos + 1;
        }
        pos++;
    }
    return false;
}

// Parse genotype directly from raw pointer
// Returns: 0 = 0/0, 1 = 0/1 or 1/0, 2 = 1/1, -1 = invalid/missing
static inline int parseGenotypeType(const char* gt, size_t gtLen) {
    if (gtLen < 3) return -1;  // Minimum: "0/0"

    // Check for missing
    if (gt[0] == '.' || (gtLen > 2 && gt[2] == '.')) return -1;

    // Parse first allele
    int a1 = 0;
    size_t pos = 0;
    while (pos < gtLen && gt[pos] >= '0' && gt[pos] <= '9') {
        a1 = a1 * 10 + (gt[pos] - '0');
        pos++;
    }

    // Skip separator (/ or |)
    if (pos >= gtLen || (gt[pos] != '/' && gt[pos] != '|')) return -1;
    pos++;

    // Check for missing second allele
    if (pos >= gtLen || gt[pos] == '.') return -1;

    // Parse second allele
    int a2 = 0;
    while (pos < gtLen && gt[pos] >= '0' && gt[pos] <= '9') {
        a2 = a2 * 10 + (gt[pos] - '0');
        pos++;
    }

    // Only handle biallelic (0 and 1)
    if (a1 > 1 || a2 > 1) return -1;

    // Return genotype type: 0=hom ref, 1=het, 2=hom alt
    return a1 + a2;
}

// ---------------------------------------------------------
// assignAncestry: OPTIMIZED - zero-allocation, dense arrays
// ---------------------------------------------------------
void VCFXAncestryAssigner::assignAncestry(std::istream &vcfIn, std::ostream &out) {
    std::string line;
    bool haveHeader = false;
    std::vector<std::string> sampleNames;

    // OPTIMIZATION: Use dense 2D vector instead of nested hash maps
    // sampleScores[sampleIdx][popIdx] => log likelihood
    std::vector<std::vector<double>> sampleScores;
    std::vector<int> sampleVariantCounts;

    size_t numPops = populations.size();

    // Reusable buffer for parsing - avoid allocations
    std::vector<std::string> fields;
    fields.reserve(16);

    // OPTIMIZATION: Pre-build key buffer to avoid repeated allocations
    std::string keyBuffer;
    keyBuffer.reserve(64);

    while (std::getline(vcfIn, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            if (line.rfind("#CHROM", 0) == 0) {
                haveHeader = true;
                vcfx::split_tabs(line, fields);

                // Sample columns start at index 9
                size_t numSamples = fields.size() > 9 ? fields.size() - 9 : 0;
                sampleNames.reserve(numSamples);
                for (size_t i = 9; i < fields.size(); ++i) {
                    sampleNames.push_back(fields[i]);
                }

                // Initialize dense arrays
                sampleScores.resize(numSamples, std::vector<double>(numPops, 0.0));
                sampleVariantCounts.resize(numSamples, 0);
            }
            continue;
        }

        if (!haveHeader) {
            std::cerr << "Error: VCF header not found.\n";
            return;
        }

        // OPTIMIZATION: Parse line with raw pointers for fixed fields
        const char* linePtr = line.c_str();
        size_t lineLen = line.size();

        // Find first 9 fields using pointer arithmetic
        const char* fieldPtrs[10];
        size_t fieldLens[10];
        int fieldCount = 0;
        size_t fieldStart = 0;

        for (size_t i = 0; i <= lineLen && fieldCount < 10; i++) {
            if (i == lineLen || linePtr[i] == '\t') {
                fieldPtrs[fieldCount] = linePtr + fieldStart;
                fieldLens[fieldCount] = i - fieldStart;
                fieldCount++;
                fieldStart = i + 1;
            }
        }

        if (fieldCount < 10) continue;

        // Build key for frequency lookup - reuse buffer
        keyBuffer.clear();
        keyBuffer.append(fieldPtrs[0], fieldLens[0]);  // CHROM
        keyBuffer.push_back(':');
        keyBuffer.append(fieldPtrs[1], fieldLens[1]);  // POS
        keyBuffer.push_back(':');
        keyBuffer.append(fieldPtrs[3], fieldLens[3]);  // REF
        keyBuffer.push_back(':');
        keyBuffer.append(fieldPtrs[4], fieldLens[4]);  // ALT

        auto freqIt = variantFrequencies.find(keyBuffer);
        if (freqIt == variantFrequencies.end()) {
            continue; // Skip variants not in frequency file
        }

        // OPTIMIZATION: Find GT index using raw pointer (no stringstream!)
        int gtIndex = findGTIndexFast(fieldPtrs[8], fieldLens[8]);
        if (gtIndex < 0) continue;

        // OPTIMIZATION: Pre-fetch frequency values for all populations (cache-friendly)
        std::vector<double> popFreqs(numPops);
        for (size_t p = 0; p < numPops; ++p) {
            double altFreq = freqIt->second.at(populations[p]);
            popFreqs[p] = std::max(0.001, std::min(0.999, altFreq));
        }

        // Process each sample - iterate through remaining tabs
        const char* samplePtr = fieldPtrs[9];
        const char* lineEnd = linePtr + lineLen;
        size_t sampleIdx = 0;

        while (samplePtr < lineEnd && sampleIdx < sampleNames.size()) {
            // Find end of current sample
            const char* sampleEnd = samplePtr;
            while (sampleEnd < lineEnd && *sampleEnd != '\t') {
                sampleEnd++;
            }
            size_t sampleLen = sampleEnd - samplePtr;

            // OPTIMIZATION: Extract GT using raw pointers (no stringstream!)
            const char* gtStart;
            size_t gtLen;
            if (extractGTField(samplePtr, sampleLen, gtIndex, gtStart, gtLen)) {
                // OPTIMIZATION: Parse genotype type directly (no string ops!)
                int genoType = parseGenotypeType(gtStart, gtLen);

                if (genoType >= 0) {
                    // Calculate log-likelihood for each population
                    for (size_t p = 0; p < numPops; ++p) {
                        double altFreq = popFreqs[p];
                        double refFreq = 1.0 - altFreq;

                        double prob = 0.0;
                        switch (genoType) {
                            case 0: prob = refFreq * refFreq; break;        // 0/0
                            case 1: prob = 2.0 * refFreq * altFreq; break;  // 0/1
                            case 2: prob = altFreq * altFreq; break;        // 1/1
                        }

                        sampleScores[sampleIdx][p] += -std::log(prob + 1e-12);
                    }
                    sampleVariantCounts[sampleIdx]++;
                }
            }

            samplePtr = (sampleEnd < lineEnd) ? sampleEnd + 1 : lineEnd;
            sampleIdx++;
        }
    }

    // Output results
    for (size_t s = 0; s < sampleNames.size(); ++s) {
        if (sampleVariantCounts[s] == 0) {
            out << sampleNames[s] << "\tUNKNOWN\n";
            continue;
        }

        // Find population with lowest negative log-likelihood (highest likelihood)
        size_t bestPopIdx = 0;
        double bestScore = sampleScores[s][0];

        std::cerr << "Final scores for " << sampleNames[s] << ":\n";
        for (size_t p = 0; p < numPops; ++p) {
            std::cerr << "  " << populations[p] << ": " << sampleScores[s][p] << "\n";
            if (sampleScores[s][p] < bestScore) {
                bestScore = sampleScores[s][p];
                bestPopIdx = p;
            }
        }

        out << sampleNames[s] << "\t" << populations[bestPopIdx] << "\n";
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
