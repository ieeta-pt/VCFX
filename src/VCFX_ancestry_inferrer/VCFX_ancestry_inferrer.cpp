#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// ----------------------------------------------------
// A helper struct for storing freq data by population
// ----------------------------------------------------
struct PopFreqKey {
    // We store CHROM, POS, REF, ALT, POP as strings
    // so we can build a single key. Alternatively, we
    // can keep them separate or use a small struct.
    // We'll build a string key "chrom:pos:ref:alt:pop"
    // for direct indexing in an unordered_map.
    std::string key;
    double frequency;
};

// We'll keep a direct structure: freqData["chr:pos:ref:alt:POP"] = frequency
typedef std::unordered_map<std::string, double> FrequencyMap;
// Map from variant key to map of pop->frequency for more efficient lookup
typedef std::unordered_map<std::string, std::unordered_map<std::string, double>> VariantPopFreqMap;

// ----------------------------------------------------
// Class: VCFXAncestryInferrer
// ----------------------------------------------------
class VCFXAncestryInferrer {
  public:
    int run(int argc, char *argv[]);

  private:
    // Show usage
    void displayHelp();

    // Load population frequencies from a file that has lines:
    // CHROM  POS  REF  ALT  POPULATION  FREQUENCY
    bool loadPopulationFrequencies(const std::string &freqFilePath);

    // Infer ancestry from VCF, returns true if successful
    bool inferAncestry(std::istream &vcfInput, std::ostream &output);

  private:
    // Frequencies keyed by "chr:pos:ref:alt:pop"
    FrequencyMap freqData;
    // More efficient structure: variant key -> (pop -> freq)
    VariantPopFreqMap variantPopFreqs;
    // Set of all known populations
    std::set<std::string> populations;
};

// ----------------------------------------------------
// main() - create the inferrer and run
// ----------------------------------------------------
static void show_help() {
    VCFXAncestryInferrer obj;
    char arg0[] = "VCFX_ancestry_inferrer";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_ancestry_inferrer", show_help))
        return 0;
    VCFXAncestryInferrer inferrer;
    return inferrer.run(argc, argv);
}

// ----------------------------------------------------
// run() - parse arguments, load freq, run inference
// ----------------------------------------------------
int VCFXAncestryInferrer::run(int argc, char *argv[]) {
    int opt;
    bool showHelp = false;
    std::string freqFilePath;

    static struct option longOpts[] = {
        {"help", no_argument, 0, 'h'}, {"frequency", required_argument, 0, 'f'}, {0, 0, 0, 0}};

    while ((opt = getopt_long(argc, argv, "hf:", longOpts, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
            break;
        case 'f':
            freqFilePath = optarg;
            break;
        default:
            showHelp = true;
        }
    }

    if (showHelp || freqFilePath.empty()) {
        displayHelp();
        // Return 0 if help was explicitly requested, otherwise 1
        return (showHelp ? 0 : 1);
    }

    // Load population frequencies
    if (!loadPopulationFrequencies(freqFilePath)) {
        std::cerr << "Error: Failed to load population frequencies from " << freqFilePath << "\n";
        return 1;
    }

    // Read VCF from stdin, write ancestry results to stdout
    if (!inferAncestry(std::cin, std::cout)) {
        return 1;
    }

    return 0;
}

// ----------------------------------------------------
// displayHelp
// ----------------------------------------------------
void VCFXAncestryInferrer::displayHelp() {
    std::cout << "VCFX_ancestry_inferrer: Infer population ancestry based on allele frequencies.\n\n"
              << "Usage:\n"
              << "  VCFX_ancestry_inferrer --frequency <freq_file> [options]\n\n"
              << "Description:\n"
              << "  Reads a VCF from standard input and outputs a 2-column table:\n"
              << "    Sample  Inferred_Population\n\n"
              << "  The frequency file must have lines of the form:\n"
              << "    CHROM  POS  REF  ALT  POPULATION  FREQUENCY\n"
              << "  (tab-separated). For multi-allelic VCF sites, an ALT allele index 1\n"
              << "  corresponds to the first item in the comma-separated ALT list,\n"
              << "  index 2 => second ALT, etc.\n\n"
              << "Example:\n"
              << "  VCFX_ancestry_inferrer --frequency pop_frequencies.txt < input.vcf > ancestry_results.txt\n";
}

// ----------------------------------------------------
// loadPopulationFrequencies
//   freq file lines: CHROM, POS, REF, ALT, POP, FREQUENCY
// ----------------------------------------------------
bool VCFXAncestryInferrer::loadPopulationFrequencies(const std::string &freqFilePath) {
    std::ifstream freqFile(freqFilePath);
    if (!freqFile.is_open()) {
        std::cerr << "Error: Cannot open frequency file: " << freqFilePath << "\n";
        return false;
    }

    int lineNum = 0;
    std::string line;
    while (std::getline(freqFile, line)) {
        lineNum++;
        if (line.empty()) {
            continue;
        }

        // Parse line
        // e.g. "chr1   12345   A   G   EUR  0.45"
        std::vector<std::string> fields;
        vcfx::split_tabs(line, fields);
        if (fields.size() < 6) {
            std::cerr << "Warning: Invalid line in frequency file (#" << lineNum << "): " << line << "\n";
            continue;
        }
        const std::string &chrom = fields[0];
        const std::string &pos = fields[1];
        const std::string &ref = fields[2];
        const std::string &alt = fields[3];
        const std::string &pop = fields[4];
        const std::string &freqStr = fields[5];

        double freq = 0.0;
        try {
            freq = std::stod(freqStr);
        } catch (...) {
            std::cerr << "Warning: Invalid frequency value in line #" << lineNum << ": " << line << "\n";
            continue;
        }

        // Build a key: "chr:pos:ref:alt:pop"
        std::string key = chrom + ":" + pos + ":" + ref + ":" + alt + ":" + pop;
        freqData[key] = freq;

        // Also store in the more efficient structure
        std::string variantKey = chrom + ":" + pos + ":" + ref + ":" + alt;
        variantPopFreqs[variantKey][pop] = freq;

        // Add to the set of known populations
        populations.insert(pop);
    }
    freqFile.close();

    if (freqData.empty()) {
        std::cerr << "Error: No valid population frequencies loaded.\n";
        return false;
    }
    return true;
}

// ----------------------------------------------------
// inferAncestry
//   1) parse the VCF header, get sample names
//   2) for each variant, parse REF + multi-ALT
//   3) for each sample's genotype, add to that sample's
//      population score = sum of freq for each ALT allele
//   4) after all lines, pick population with highest score
// ----------------------------------------------------
bool VCFXAncestryInferrer::inferAncestry(std::istream &vcfInput, std::ostream &out) {
    std::string line;
    bool foundChromHeader = false;
    std::vector<std::string> headerFields;
    std::vector<std::string> sampleNames;

    // We'll accumulate ancestry scores: sample -> (pop -> score)
    // Then we pick the highest for each sample
    std::unordered_map<std::string, std::unordered_map<std::string, double>> sampleScores;

    // Reusable buffer for parsing
    std::vector<std::string> fields;
    fields.reserve(16);

    // Check if the input stream is valid
    if (!vcfInput) {
        std::cerr << "Error: Invalid VCF input stream.\n";
        return false;
    }

    int lineNum = 0;
    while (std::getline(vcfInput, line)) {
        lineNum++;
        if (line.empty()) {
            continue;
        }

        // Check if line is a VCF header
        if (line[0] == '#') {
            // The #CHROM line is the last header line, we parse columns
            if (line.rfind("#CHROM", 0) == 0) {
                foundChromHeader = true;
                headerFields.clear();
                vcfx::split_tabs(line, headerFields);

                // Validate header has at least the required fields
                if (headerFields.size() < 9) {
                    std::cerr << "Error: Invalid VCF header format. Expected at least 9 columns.\n";
                    return false;
                }

                // sample columns start at index 9
                for (size_t c = 9; c < headerFields.size(); ++c) {
                    sampleNames.push_back(headerFields[c]);
                    // Initialize scores for each known population to 0
                    for (const auto &pop : populations) {
                        sampleScores[headerFields[c]][pop] = 0.0;
                    }
                }
            }
            continue;
        }

        // We must have #CHROM header before data
        if (!foundChromHeader) {
            std::cerr << "Error: Encountered VCF data before #CHROM header.\n";
            return false;
        }

        // Parse data line
        // Format (minimally): CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [samples...]
        fields.clear();
        vcfx::split_tabs(line, fields);

        if (fields.size() < 10) {
            // Not enough columns to have samples, skip line
            std::cerr << "Warning: Line " << lineNum << " has fewer than 10 columns, skipping.\n";
            continue;
        }

        // Indices: 0=CHROM, 1=POS, 2=ID, 3=REF, 4=ALT, 5=QUAL, 6=FILTER, 7=INFO, 8=FORMAT, 9+ = samples
        const std::string &chrom = fields[0];
        const std::string &pos = fields[1];
        // skip ID
        const std::string &ref = fields[3];
        const std::string &altStr = fields[4];
        const std::string &format = fields[8];

        // Split ALT by comma for multi-allelic
        std::vector<std::string> altAlleles;
        {
            std::stringstream altSS(altStr);
            std::string altTok;
            while (std::getline(altSS, altTok, ',')) {
                altAlleles.push_back(altTok);
            }
        }

        // Validate ALT field
        if (altAlleles.empty()) {
            std::cerr << "Warning: Line " << lineNum << " has empty ALT field, skipping.\n";
            continue;
        }

        // Find GT index in format
        // e.g. format might be "GT:DP:GQ:..."
        std::vector<std::string> formatParts;
        {
            std::stringstream fmts(format);
            std::string fmtTok;
            while (std::getline(fmts, fmtTok, ':')) {
                formatParts.push_back(fmtTok);
            }
        }
        int gtIndex = -1;
        for (size_t i = 0; i < formatParts.size(); ++i) {
            if (formatParts[i] == "GT") {
                gtIndex = (int)i;
                break;
            }
        }
        if (gtIndex < 0) {
            // no genotype in this line
            std::cerr << "Warning: Line " << lineNum << " has no GT field, skipping.\n";
            continue;
        }

        // For each sample
        for (size_t s = 0; s < sampleNames.size(); ++s) {
            // sample data is fields[9 + s]
            size_t sampleCol = 9 + s;
            if (sampleCol >= fields.size()) {
                continue;
            }
            const std::string &sampleData = fields[sampleCol];

            // Split by ':'
            std::vector<std::string> sampleParts;
            {
                std::stringstream sampSS(sampleData);
                std::string part;
                while (std::getline(sampSS, part, ':')) {
                    sampleParts.push_back(part);
                }
            }
            if (gtIndex >= (int)sampleParts.size()) {
                // no GT
                continue;
            }

            std::string genotype = sampleParts[gtIndex]; // e.g. "0/1" or "2|1"
            // unify separators
            std::replace(genotype.begin(), genotype.end(), '|', '/');
            std::vector<std::string> alleleNums;
            {
                std::stringstream gtSS(genotype);
                std::string a;
                while (std::getline(gtSS, a, '/')) {
                    alleleNums.push_back(a);
                }
            }
            // We expect diploid => 2 allele calls, but handle if 1 or 2
            if (alleleNums.empty()) {
                continue;
            }

            // For each allele: 0 => REF, 1 => altAlleles[0], 2 => altAlleles[1], etc.
            for (auto &aStr : alleleNums) {
                if (aStr.empty() || aStr == ".") {
                    continue; // missing
                }

                // Validate allele is numeric
                bool numeric = true;
                for (char c : aStr) {
                    if (!isdigit(c)) {
                        numeric = false;
                        break;
                    }
                }
                if (!numeric) {
                    continue;
                }

                int aVal = std::stoi(aStr);
                if (aVal == 0) {
                    // REF allele, skip
                    continue;
                }

                if (aVal > 0 && (size_t)aVal <= altAlleles.size()) {
                    // This is an ALT allele
                    std::string actualAlt = altAlleles[aVal - 1];

                    // Build the variant key for lookup
                    std::string variantKey = chrom + ":" + pos + ":" + ref + ":" + actualAlt;

                    // Use the more efficient lookup structure
                    auto variantIt = variantPopFreqs.find(variantKey);
                    if (variantIt != variantPopFreqs.end()) {
                        const auto &popFreqs = variantIt->second;

                        // Find the population with the highest frequency
                        double bestFreq = -1.0;
                        std::string bestPop;

                        for (const auto &popFreq : popFreqs) {
                            if (popFreq.second > bestFreq) {
                                bestFreq = popFreq.second;
                                bestPop = popFreq.first;
                            }
                        }

                        if (!bestPop.empty() && bestFreq >= 0.0) {
                            // Add the bestFreq to the sample's population score
                            sampleScores[sampleNames[s]][bestPop] += bestFreq;
                        }
                    }
                }
                // Skip if aVal > altAlleles.size()
            }
        }
    }

    // Check if we read any valid VCF data
    if (!foundChromHeader) {
        std::cerr << "Error: No valid VCF data found in input.\n";
        return false;
    }

    // Output the results
    out << "Sample\tInferred_Population\n";
    for (auto &sName : sampleNames) {
        // sampleScores[sName] => map<pop, score>
        auto it = sampleScores.find(sName);
        if (it == sampleScores.end() || it->second.empty()) {
            // no data => unknown
            out << sName << "\tUnknown\n";
            continue;
        }

        const auto &popMap = it->second;
        std::string bestPop = "Unknown";
        double bestScore = -1.0;

        for (auto &ps : popMap) {
            if (ps.second > bestScore) {
                bestScore = ps.second;
                bestPop = ps.first;
            }
        }

        out << sName << "\t" << bestPop << "\n";
    }

    return true;
}
