#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <sstream>

// ---------------------------------------------------------
// Data structure for command-line arguments
// ---------------------------------------------------------
struct ConcordanceArguments {
    std::string sample1;
    std::string sample2;
};

// ---------------------------------------------------------
// Print help
// ---------------------------------------------------------
static void printHelp() {
    std::cout << "VCFX_concordance_checker\n"
              << "Usage: VCFX_concordance_checker [OPTIONS] < input.vcf > concordance_report.tsv\n\n"
              << "Options:\n"
              << "  --samples, -s \"Sample1 Sample2\"  Specify exactly two sample names to compare.\n"
              << "  --help, -h                        Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Compares genotypes between two specified samples in a VCF file, including multi-allelic\n"
              << "  variants, and outputs per-variant concordance (Concordant or Discordant).\n\n"
              << "Example:\n"
              << "  ./VCFX_concordance_checker --samples \"SampleA SampleB\" < input.vcf > concordance_report.tsv\n";
}

// ---------------------------------------------------------
// Utility: split a string by delimiter
// ---------------------------------------------------------
static std::vector<std::string> splitString(const std::string &str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string temp;
    while (std::getline(ss, temp, delimiter)) {
        tokens.push_back(temp);
    }
    return tokens;
}

// ---------------------------------------------------------
// parseArguments: fill ConcordanceArguments
// ---------------------------------------------------------
static bool parseArguments(int argc, char *argv[], ConcordanceArguments &args) {
    bool showHelp = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--samples" || arg == "-s") && i + 1 < argc) {
            std::string samples_str = argv[++i];
            auto samples = splitString(samples_str, ' ');
            if (samples.size() != 2) {
                std::cerr << "Error: Please specify exactly two sample names (e.g., -s \"Sample1 Sample2\").\n";
                return false;
            }
            args.sample1 = samples[0];
            args.sample2 = samples[1];
        } else if (arg == "--help" || arg == "-h") {
            showHelp = true;
        } else {
            std::cerr << "Warning: Unrecognized argument '" << arg << "'.\n";
        }
    }

    if (showHelp) {
        printHelp();
        // Return false so main() can exit cleanly
        return false;
    }
    if (args.sample1.empty() || args.sample2.empty()) {
        std::cerr << "Error: Two sample names must be specified using --samples or -s.\n";
        std::cerr << "Use --help for usage information.\n";
        return false;
    }

    return true;
}

// ---------------------------------------------------------
// Convert a numeric allele index to a textual tag:
//   0 => "0"
//   1 => "1" (meaning altAlleles[0])
//   2 => "2" (meaning altAlleles[1]), etc.
// We'll store them as strings "0","1","2" for comparison
// ---------------------------------------------------------
static std::string normalizeDiploidGenotype(const std::string &genotypeField,
                                            const std::vector<std::string> &altAlleles) {
    // genotypeField might be "0/1:30:90" or "2|3", ".", ...
    // Extract GT (the part before the first ':')
    std::string gt = genotypeField;
    size_t colonPos = gt.find(':');
    if (colonPos != std::string::npos) {
        gt = gt.substr(0, colonPos);
    }
    // Step 1: replace '|' with '/'
    for (char &c : gt) {
        if (c == '|')
            c = '/';
    }
    // split by '/'
    auto alleles = splitString(gt, '/');
    if (alleles.size() < 2) {
        // not diploid or missing
        return "";
    }
    // parse each allele
    std::vector<int> numericAll;
    numericAll.reserve(2);
    for (auto &a : alleles) {
        if (a == "." || a.empty()) {
            // missing
            return "";
        }
        // must parse as an integer
        bool numeric = std::all_of(a.begin(), a.end(), ::isdigit);
        if (!numeric) {
            return "";
        }
        int val = std::stoi(a);
        // If val > altAlleles.size(), we can't interpret => skip
        // val==0 => reference
        // val==1 => altAlleles[0]
        // val==2 => altAlleles[1], etc.
        if (val < 0) {
            return "";
        }
        if (val > 0) {
            // if val=2 => we need altAlleles.size() >=2
            if ((size_t)val > altAlleles.size()) {
                // out of range
                return "";
            }
        }
        numericAll.push_back(val);
    }
    // We produce a string "x/y" sorted or not? Typically for genotype comparison,
    // "1/0" is same as "0/1", but let's keep the original order. It's also common
    // to store them in numeric order. We'll store them in sorted numeric order
    // so that "0/1" == "1/0". Then we join with "/".
    std::sort(numericAll.begin(), numericAll.end());
    std::ostringstream oss;
    oss << numericAll[0] << "/" << numericAll[1];
    return oss.str();
}

// ---------------------------------------------------------
// Calculate concordance
//   1) parse the #CHROM line, find sample1_index + sample2_index
//   2) for each variant, parse ALT, parse sample1 + sample2 genotype => normalized
//   3) if both genotypes are non-empty => compare => print row
//   4) track summary stats
// ---------------------------------------------------------
static bool calculateConcordance(std::istream &in, std::ostream &out, const ConcordanceArguments &args) {
    std::string line;
    bool foundChromHeader = false;
    std::vector<std::string> headerFields;
    int sample1_index = -1;
    int sample2_index = -1;

    // We'll track summary stats
    int totalVariants = 0;
    int concordant = 0;
    int discordant = 0;

    // Print a TSV header:
    // CHROM POS ID REF ALT(S) Sample1_GT Sample2_GT Concordance
    out << "CHROM\tPOS\tID\tREF\tALT\t" << args.sample1 << "_GT\t" << args.sample2 << "_GT\tConcordance\n";

    // Reusable vector for tab splitting
    std::vector<std::string> fields;
    fields.reserve(16);

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '#') {
            // If it's the #CHROM line, parse sample columns
            if (!foundChromHeader && line.rfind("#CHROM", 0) == 0) {
                vcfx::split_tabs(line, headerFields);
                // find sample1, sample2
                std::unordered_map<std::string, int> sampleMap;
                for (size_t i = 9; i < headerFields.size(); ++i) {
                    sampleMap[headerFields[i]] = (int)i;
                }
                // check existence
                auto it1 = sampleMap.find(args.sample1);
                if (it1 == sampleMap.end()) {
                    std::cerr << "Error: Sample '" << args.sample1 << "' not found in VCF header.\n";
                    return false;
                }
                auto it2 = sampleMap.find(args.sample2);
                if (it2 == sampleMap.end()) {
                    std::cerr << "Error: Sample '" << args.sample2 << "' not found in VCF header.\n";
                    return false;
                }
                sample1_index = it1->second;
                sample2_index = it2->second;
                foundChromHeader = true;
            }
            continue;
        }

        if (!foundChromHeader) {
            std::cerr << "Error: VCF data encountered before #CHROM header.\n";
            return false;
        }

        // parse data line
        vcfx::split_tabs(line, fields);
        if (fields.size() < 8) {
            // invalid
            continue;
        }
        // check we have sample columns for sample1, sample2
        if ((size_t)sample1_index >= fields.size() || (size_t)sample2_index >= fields.size()) {
            // line doesn't have enough columns
            continue;
        }

        const std::string &chrom = fields[0];
        const std::string &pos = fields[1];
        const std::string &id = fields[2];
        const std::string &ref = fields[3];
        const std::string &alt = fields[4];
        // sample columns:
        const std::string &sample1_field = fields[sample1_index];
        const std::string &sample2_field = fields[sample2_index];

        // Split the alt string (but for genotype normalization we only need the count)
        auto altAlleles = splitString(alt, ',');

        // Convert sample1 genotype to a normalized diploid string
        std::string s1_gt = normalizeDiploidGenotype(sample1_field, altAlleles);
        std::string s2_gt = normalizeDiploidGenotype(sample2_field, altAlleles);

        // If either is missing or un-parseable, skip
        if (s1_gt.empty() || s2_gt.empty()) {
            continue;
        }

        totalVariants++;
        bool same = (s1_gt == s2_gt);
        if (same) {
            concordant++;
        } else {
            discordant++;
        }
        std::string cc = same ? "Concordant" : "Discordant";

        // Print the row
        out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\t" << s1_gt << "\t" << s2_gt
            << "\t" << cc << "\n";
    }

    // Print summary to stderr so we don't break the TSV
    std::cerr << "Total Variants Compared: " << totalVariants << "\n"
              << "Concordant Genotypes: " << concordant << "\n"
              << "Discordant Genotypes: " << discordant << "\n";

    return true;
}

// ---------------------------------------------------------
// main
// ---------------------------------------------------------
static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_concordance_checker", show_help))
        return 0;
    ConcordanceArguments args;
    if (!parseArguments(argc, argv, args)) {
        // parseArguments prints error/help if needed
        return 1;
    }

    bool ok = calculateConcordance(std::cin, std::cout, args);
    return (ok ? 0 : 1);
}
