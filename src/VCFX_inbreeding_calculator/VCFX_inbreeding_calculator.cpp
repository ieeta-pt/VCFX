#include "VCFX_inbreeding_calculator.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <unordered_map>

// --------------------------------------------------------------------------
// A helper function to remove BOM and leading/trailing whitespace
static void sanitizeLine(std::string &line) {
    // Possibly remove a UTF-8 BOM
    static const std::string BOM = "\xEF\xBB\xBF";
    if (line.size() >= BOM.size() && line.compare(0, BOM.size(), BOM) == 0) {
        line.erase(0, BOM.size());
    }
    // Trim front
    while (!line.empty() && (line[0] == ' ' || line[0] == '\t' || line[0] == '\r' || line[0] == '\n')) {
        line.erase(line.begin());
    }
    // Trim back
    while (!line.empty() && (line.back() == ' ' || line.back() == '\t' || line.back() == '\r' || line.back() == '\n')) {
        line.pop_back();
    }
}

// --------------------------------------------------------------------------
// Split a line on '\t', removing trailing CR if present
static std::vector<std::string> splitByTab(const std::string &line) {
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string f;
    while (std::getline(ss, f, '\t')) {
        if (!f.empty() && f.back() == '\r') {
            f.pop_back();
        }
        fields.push_back(f);
    }
    return fields;
}

// --------------------------------------------------------------------------
void VCFXInbreedingCalculator::displayHelp() {
    std::cout << "VCFX_inbreeding_calculator: Compute individual inbreeding coefficients (F)\n"
              << "based on biallelic sites in a VCF.\n\n"
              << "Usage:\n"
              << "  VCFX_inbreeding_calculator [options] < input.vcf > output.txt\n\n"
              << "Description:\n"
              << "  Reads a VCF in a single pass, ignoring multi-allelic lines (ALT with commas).\n"
              << "  For each biallelic variant, we parse each sample's genotype code:\n"
              << "       0/0 => 0,   0/1 => 1,   1/1 => 2, else => -1 (ignored)\n\n"
              << "  Then, depending on --freq-mode:\n"
              << "    * excludeSample => Each sample excludes its own genotype when computing p.\n"
              << "    * global        => Compute a single global p from all samples' genotypes.\n\n"
              << "  The --skip-boundary option, if set, ignores boundary freq p=0 or p=1.\n"
              << "    BUT if you also specify --count-boundary-as-used, those boundary sites\n"
              << "    increment usedCount (forcing F=1) without contributing to sumExp.\n\n"
              << "  If sumExp=0 for a sample but usedCount>0, we output F=1.\n"
              << "  If usedCount=0, we output NA.\n\n"
              << "Options:\n"
              << "  -h, --help                Show this help.\n"
              << "  --freq-mode <mode>        'excludeSample' (default) or 'global'\n"
              << "  --skip-boundary           Skip boundary freq sites. By default, they are used.\n"
              << "  --count-boundary-as-used  If also skipping boundary, still increment usedCount.\n"
              << std::endl;
}

// --------------------------------------------------------------------------
int VCFXInbreedingCalculator::parseGenotype(const std::string &s) {
    // Check for missing or empty
    if (s.empty() || s == "." || s == "./." || s == ".|." || s == ".|" || s == "./") {
        return -1; // missing
    }
    std::string g = s;
    for (char &c : g) {
        if (c == '|')
            c = '/';
    }
    size_t slash = g.find('/');
    if (slash == std::string::npos) {
        return -1; // not diploid
    }
    std::string a1 = g.substr(0, slash);
    std::string a2 = g.substr(slash + 1);
    if (a1.empty() || a2.empty() || a1 == "." || a2 == ".") {
        return -1;
    }
    int i1, i2;
    try {
        i1 = std::stoi(a1);
        i2 = std::stoi(a2);
    } catch (...) {
        return -1;
    }
    // skip if allele≥2 => multi-allelic genotype (since we only handle 0 or 1)
    if (i1 < 0 || i1 > 1 || i2 < 0 || i2 > 1) {
        return -1;
    }
    // 0/0 =>0, 1/1 =>2, else =>1
    if (i1 == i2) {
        return (i1 == 0 ? 0 : 2);
    }
    return 1; // 0/1
}

// --------------------------------------------------------------------------
bool VCFXInbreedingCalculator::isBiallelic(const std::string &alt) { return (alt.find(',') == std::string::npos); }

// --------------------------------------------------------------------------
FrequencyMode VCFXInbreedingCalculator::parseFreqMode(const std::string &modeStr) {
    if (modeStr == "excludeSample") {
        return FrequencyMode::EXCLUDE_SAMPLE;
    } else if (modeStr == "global") {
        return FrequencyMode::GLOBAL;
    }
    // Default
    std::cerr << "Warning: unrecognized freq-mode='" << modeStr << "'. Using 'excludeSample' by default.\n";
    return FrequencyMode::EXCLUDE_SAMPLE;
}

// --------------------------------------------------------------------------
void VCFXInbreedingCalculator::calculateInbreeding(std::istream &in, std::ostream &out) {
    bool foundChrom = false;
    std::vector<std::string> sampleNames;
    int numSamples = 0;
    std::vector<InbreedingVariant> variants;

    // 1) parse lines
    while (true) {
        std::string line;
        if (!std::getline(in, line)) {
            break; // EOF
        }
        sanitizeLine(line);
        if (line.empty())
            continue;

        if (line[0] == '#') {
            // check if #CHROM
            if (!foundChrom && line.find("#CHROM") != std::string::npos) {
                foundChrom = true;
                auto tokens = splitByTab(line);
                if (tokens.size() >= 10) {
                    for (size_t i = 9; i < tokens.size(); i++) {
                        sampleNames.push_back(tokens[i]);
                    }
                    numSamples = (int)sampleNames.size();
                }
            }
            continue;
        }
        // data line
        if (!foundChrom) {
            // ignore data lines if #CHROM not found
            continue;
        }
        auto fields = splitByTab(line);
        if (fields.size() < 10) {
            continue;
        }
        if (!isBiallelic(fields[4])) {
            continue; // skip multi-allelic
        }
        InbreedingVariant var;
        var.chrom = fields[0];
        try {
            var.pos = std::stoi(fields[1]);
        } catch (...) {
            continue;
        }
        var.genotypeCodes.resize(numSamples, -1);
        for (int s = 0; s < numSamples; s++) {
            int colIndex = 9 + s;
            if ((size_t)colIndex >= fields.size())
                break;
            var.genotypeCodes[s] = parseGenotype(fields[colIndex]);
        }
        variants.push_back(var);
    }

    // If no #CHROM found or no variants, handle gracefully
    if (!foundChrom) {
        std::cerr << "Error: No #CHROM line found.\n";
        out << "Sample\tInbreedingCoefficient\n";
        return;
    }
    if (numSamples == 0) {
        std::cerr << "Error: No sample columns found.\n";
        out << "Sample\tInbreedingCoefficient\n";
        return;
    }
    if (variants.empty()) {
        // No valid (biallelic) variants
        std::cerr << "No biallelic variants found.\n";
        out << "Sample\tInbreedingCoefficient\n";
        for (int s = 0; s < numSamples; s++) {
            out << sampleNames[s] << "\tNA\n";
        }
        return;
    }

    // Prepare result arrays
    std::vector<double> sumExp(numSamples, 0.0);
    std::vector<double> obsHet(numSamples, 0.0);
    std::vector<int> usedCount(numSamples, 0);

    // 2) For each variant
    for (size_t v = 0; v < variants.size(); v++) {
        const auto &var = variants[v];
        const auto &codes = var.genotypeCodes;

        // Count altSum & nGood (i.e. how many samples have code>=0)
        int altSum = 0;
        int nGood = 0;
        for (int s = 0; s < numSamples; s++) {
            if (codes[s] >= 0) {
                altSum += codes[s];
                nGood++;
            }
        }
        // If fewer than 2 genotyped samples, skip
        if (nGood < 2) {
            continue;
        }

        // For freq-mode=GLOBAL, we compute p once
        double globalP = (double)altSum / (2.0 * nGood);

        // 3) Each sample that has a valid genotype
        for (int s = 0; s < numSamples; s++) {
            int code = codes[s];
            if (code < 0) {
                continue; // missing or invalid
            }

            double p;
            if (freqMode_ == FrequencyMode::GLOBAL) {
                p = globalP;
            } else {
                // freqMode_ == FrequencyMode::EXCLUDE_SAMPLE
                int altEx = altSum - code;
                int validEx = nGood - 1;
                if (validEx < 1) {
                    // can't define p => skip
                    continue;
                }
                p = (double)altEx / (2.0 * validEx);
            }

            // If skipBoundary_ is set and p is boundary
            if (skipBoundary_ && (p <= 0.0 || p >= 1.0)) {
                // If user wants to count boundary as used
                if (countBoundaryAsUsed_) {
                    // This increments usedCount but adds 0 to sumExp
                    usedCount[s]++;
                }
                // Then continue, so we don't add eHet
                continue;
            }

            // Mark that we used this site for sample s
            usedCount[s]++;

            // expected heterozygosity
            double eHet = 2.0 * p * (1.0 - p);
            sumExp[s] += eHet;

            // observed heterozygosity increment if 0/1 => code=1
            if (code == 1) {
                obsHet[s] += 1.0;
            }
        }
    }

    // 4) Output
    out << "Sample\tInbreedingCoefficient\n";
    for (int s = 0; s < numSamples; s++) {
        if (usedCount[s] == 0) {
            // No used sites => NA
            out << sampleNames[s] << "\tNA\n";
            continue;
        }
        double e = sumExp[s];
        if (e <= 0.0) {
            // used some sites => but sumExp=0 => fully homo => F=1
            out << sampleNames[s] << "\t1.000000\n";
            continue;
        }
        double f = 1.0 - (obsHet[s] / e);
        out << sampleNames[s] << "\t" << std::fixed << std::setprecision(6) << f << "\n";
    }
}

// --------------------------------------------------------------------------
int VCFXInbreedingCalculator::run(int argc, char *argv[]) {
    bool showHelp = false;

    static struct option long_opts[] = {{"help", no_argument, 0, 'h'},
                                        {"freq-mode", required_argument, 0, 0},
                                        {"skip-boundary", no_argument, 0, 0},
                                        {"count-boundary-as-used", no_argument, 0, 0},
                                        {0, 0, 0, 0}};

    // Parse arguments
    while (true) {
        int option_index = 0;
        int c = getopt_long(argc, argv, "h", long_opts, &option_index);
        if (c == -1)
            break;

        switch (c) {
        case 'h':
            showHelp = true;
            break;
        case 0: {
            // We got a long option
            std::string optName(long_opts[option_index].name);
            if (optName == "freq-mode") {
                freqMode_ = parseFreqMode(optarg);
            } else if (optName == "skip-boundary") {
                skipBoundary_ = true;
            } else if (optName == "count-boundary-as-used") {
                countBoundaryAsUsed_ = true;
            }
            break;
        }
        default:
            showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    // Run main logic
    calculateInbreeding(std::cin, std::cout);
    return 0;
}

// -------------------------------------------------------------------------
// main entry point
static void show_help() {
    VCFXInbreedingCalculator obj;
    char arg0[] = "VCFX_inbreeding_calculator";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_inbreeding_calculator", show_help))
        return 0;
    VCFXInbreedingCalculator calc;
    return calc.run(argc, argv);
}
