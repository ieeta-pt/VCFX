#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <getopt.h>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// --------------------------------------------------------------------------
// A helper function to print usage/help
// --------------------------------------------------------------------------
static void displayHelp() {
    std::cout << "VCFX_cross_sample_concordance: Check variant concordance across multiple samples.\n\n"
              << "Usage:\n"
              << "  VCFX_cross_sample_concordance [options] < input.vcf > concordance_results.txt\n\n"
              << "Options:\n"
              << "  -h, --help              Display this help message and exit\n"
              << "  -s, --samples LIST      Comma-separated list of samples to check\n\n"
              << "Description:\n"
              << "  Reads a multi-sample VCF from stdin, normalizes each sample's genotype\n"
              << "  (including multi-allelic variants), and determines if all samples that\n"
              << "  have a parseable genotype are in complete agreement. Prints one row per\n"
              << "  variant to stdout and a final summary to stderr.\n\n"
              << "Example:\n"
              << "  VCFX_cross_sample_concordance < input.vcf > results.tsv\n\n";
}

// --------------------------------------------------------------------------
// Utility: split a string by a delimiter (optimized - no stringstream)
// --------------------------------------------------------------------------
static void splitInto(const std::string &str, char delim, std::vector<std::string> &out) {
    out.clear();
    size_t start = 0, end;
    while ((end = str.find(delim, start)) != std::string::npos) {
        out.emplace_back(str, start, end - start);
        start = end + 1;
    }
    out.emplace_back(str, start);
}

// Convenience wrapper for compatibility
static std::vector<std::string> split(const std::string &str, char delim) {
    std::vector<std::string> tokens;
    splitInto(str, delim, tokens);
    return tokens;
}

// --------------------------------------------------------------------------
// Normalize a genotype string (e.g. "2|1", "0/1") to a canonical form
//   - Replace '|' with '/'
//   - Split by '/'
//   - Check for missing or invalid (out-of-range) => return empty
//   - Convert to int, store in vector
//   - Sort so "2/1" => "1/2"
//   - Return e.g. "1/2"
// --------------------------------------------------------------------------
static std::string normalizeGenotype(const std::string &sampleField, const std::vector<std::string> &altAlleles) {
    // sampleField might be something like "0/1:.." or just "0/1"
    // We only care about the genotype portion (first colon-delimited field).
    auto parts = split(sampleField, ':');
    if (parts.empty()) {
        return "";
    }
    std::string gt = parts[0]; // e.g. "0/1", "2|3", etc.

    // unify separators
    for (char &c : gt) {
        if (c == '|')
            c = '/';
    }
    // split by '/'
    auto alleleStrs = split(gt, '/');
    if (alleleStrs.size() < 2) {
        // might be haploid or invalid
        return "";
    }

    std::vector<int> alleleInts;
    alleleInts.reserve(alleleStrs.size());
    for (auto &a : alleleStrs) {
        if (a == "." || a.empty()) {
            return ""; // missing
        }
        // parse
        for (char c : a) {
            if (!std::isdigit(c)) {
                return "";
            }
        }
        int val = std::stoi(a);
        // If val == 0 => REF
        // If val == 1 => altAlleles[0]
        // If val == 2 => altAlleles[1], etc.
        if (val < 0) {
            return "";
        }
        if (val > 0 && (size_t)val > altAlleles.size()) {
            // out of range for the alt array
            return "";
        }
        alleleInts.push_back(val);
    }

    // For genotype comparison, sort them so "2/1" => "1/2"
    std::sort(alleleInts.begin(), alleleInts.end());
    // Build a normalized string - direct string building (no stringstream)
    std::string result;
    result.reserve(alleleInts.size() * 2);  // Pre-allocate for digits and slashes
    for (size_t i = 0; i < alleleInts.size(); ++i) {
        if (i > 0)
            result.push_back('/');
        // Fast integer to string for small values (0-99 covers all realistic allele indices)
        int val = alleleInts[i];
        if (val < 10) {
            result.push_back(static_cast<char>('0' + val));
        } else {
            result.push_back(static_cast<char>('0' + val / 10));
            result.push_back(static_cast<char>('0' + val % 10));
        }
    }
    return result;
}

// --------------------------------------------------------------------------
// The main function to process the VCF and produce cross-sample concordance
// --------------------------------------------------------------------------
static void calculateConcordance(std::istream &in, std::ostream &out, const std::vector<std::string> &subset) {
    std::string line;
    std::vector<std::string> sampleNames;
    bool gotChromHeader = false;
    std::vector<int> sampleIndices;

    // We'll print a header for our TSV
    // CHROM POS ID REF ALT Num_Samples Unique_Normalized_Genotypes Concordance_Status
    out << "CHROM\tPOS\tID\tREF\tALT\tNum_Samples\tUnique_Normalized_Genotypes\tConcordance_Status\n";

    // We track summary stats
    size_t totalVariants = 0;
    size_t allConcordantCount = 0;
    size_t discordantCount = 0;
    size_t skippedBecauseNoGenotypes = 0; // if all samples are missing/un-parseable

    // Declare fields vector outside loop for reuse
    std::vector<std::string> headers;
    headers.reserve(16);

    // 1) Read until we find the #CHROM line
    //    then parse sample names from columns 9+
    while (!gotChromHeader && std::getline(in, line)) {
        if (line.rfind("#CHROM", 0) == 0) {
            // parse columns
            vcfx::split_tabs(line, headers);
            // from index=9 onwards => sample names
            for (size_t i = 9; i < headers.size(); ++i) {
                sampleNames.push_back(headers[i]);
            }
            // Determine which sample indices to process
            if (subset.empty()) {
                for (size_t i = 0; i < sampleNames.size(); ++i) {
                    sampleIndices.push_back(static_cast<int>(i));
                }
            } else {
                std::unordered_set<std::string> want(subset.begin(), subset.end());
                for (size_t i = 0; i < sampleNames.size(); ++i) {
                    if (want.count(sampleNames[i])) {
                        sampleIndices.push_back(static_cast<int>(i));
                    }
                }
                if (sampleIndices.empty()) {
                    std::cerr << "Warning: none of the requested samples were found in the VCF header.\n";
                }
            }
            gotChromHeader = true;
            break; // proceed to reading data lines
        }
    }
    if (!gotChromHeader) {
        std::cerr << "Error: VCF header with #CHROM not found.\n";
        return;
    }

    // Declare fields vector outside loop for reuse
    std::vector<std::string> fields;
    fields.reserve(16);

    // 2) Now read the data lines
    while (std::getline(in, line)) {
        // skip header lines or empties
        if (line.empty() || line[0] == '#') {
            continue;
        }
        vcfx::split_tabs(line, fields);
        // minimal VCF => 8 columns + samples => need at least 10 if there's 1 sample
        if (fields.size() < 8) {
            // invalid line
            continue;
        }

        // parse standard columns
        const std::string &chrom = fields[0];
        const std::string &pos = fields[1];
        const std::string &id = fields[2];
        const std::string &ref = fields[3];
        const std::string &alt = fields[4];
        // we can skip qual/filter/info/format => we only need sample columns
        // sample columns start at index=9
        if (fields.size() < (9 + sampleNames.size())) {
            // line doesn't have enough columns for all samples
            // skip
            continue;
        }

        // parse alt for multi-allelic
        auto altAlleles = split(alt, ',');

        // For each sample, normalize genotype
        std::vector<std::string> normalizedGts;
        normalizedGts.reserve(sampleIndices.size());
        size_t sampleCountHere = 0;

        for (int idx : sampleIndices) {
            auto &sampleField = fields[9 + idx];
            std::string norm = normalizeGenotype(sampleField, altAlleles);
            if (!norm.empty()) {
                // we have a parseable genotype
                normalizedGts.push_back(norm);
                sampleCountHere++;
            }
        }

        // If no samples had a valid genotype, mark this variant as NO_GENOTYPES
        if (sampleCountHere == 0) {
            skippedBecauseNoGenotypes++;
            // Output with NO_GENOTYPES status instead of skipping
            out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\t" << sampleCountHere << "\t"
                << 0 // No unique genotypes
                << "\t" << "NO_GENOTYPES"
                << "\n";
            continue;
        }

        // We have at least 1 genotype
        totalVariants++;

        // gather unique genotype strings
        std::unordered_map<std::string, int> freq;
        for (auto &g : normalizedGts) {
            freq[g]++;
        }
        size_t uniqueGtCount = freq.size();
        bool allConcordant = (uniqueGtCount == 1);

        if (allConcordant) {
            allConcordantCount++;
        } else {
            discordantCount++;
        }

        std::string status = allConcordant ? "CONCORDANT" : "DISCORDANT";
        // Print row
        out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\t" << sampleCountHere << "\t"
            << uniqueGtCount << "\t" << status << "\n";
    }

    // Summaries to stderr so as not to pollute the TSV
    std::cerr << "Total Variants with >=1 parseable genotype: " << totalVariants << "\n";
    std::cerr << "   Concordant (all same genotype): " << allConcordantCount << "\n";
    std::cerr << "   Discordant (>=2 distinct genotypes): " << discordantCount << "\n";
    std::cerr << "Variants with no parseable genotypes (skipped): " << skippedBecauseNoGenotypes << "\n";
}

// --------------------------------------------------------------------------
// Command-line parsing + main
// --------------------------------------------------------------------------
static void show_help() { displayHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_cross_sample_concordance", show_help))
        return 0;
    bool showHelp = false;
    std::string samplesArg;

    static struct option longOpts[] = {
        {"help", no_argument, 0, 'h'}, {"samples", required_argument, 0, 's'}, {0, 0, 0, 0}};

    while (true) {
        int optIndex = 0;
        int c = getopt_long(argc, argv, "hs:", longOpts, &optIndex);
        if (c == -1)
            break;
        switch (c) {
        case 'h':
            showHelp = true;
            break;
        case 's':
            samplesArg = optarg;
            break;
        default:
            showHelp = true;
            break;
        }
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    std::vector<std::string> subsetSamples;
    if (!samplesArg.empty()) {
        subsetSamples = split(samplesArg, ',');
    }

    calculateConcordance(std::cin, std::cout, subsetSamples);
    return 0;
}
