#include "VCFX_duplicate_remover.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <vector>
#include <sstream>

// Utility function: Splits a string by a given delimiter (optimized).
static std::vector<std::string> splitString(const std::string &str, char delimiter) {
    std::vector<std::string> tokens;
    tokens.reserve(8);
    size_t start = 0;
    size_t end;
    while ((end = str.find(delimiter, start)) != std::string::npos) {
        tokens.emplace_back(str, start, end - start);
        start = end + 1;
    }
    tokens.emplace_back(str, start);
    return tokens;
}

// Function to display help message
void printHelp() {
    std::cout << "VCFX_duplicate_remover\n"
              << "Usage: VCFX_duplicate_remover [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h               Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Removes duplicate variants from a VCF file based on the combination of\n"
              << "  chromosome, position, REF, and ALT alleles. For multi-allelic records, the\n"
              << "  ALT field is normalized by sorting the comma-separated alleles so that the\n"
              << "  ordering does not affect duplicate detection.\n\n"
              << "Example:\n"
              << "  ./VCFX_duplicate_remover < input.vcf > unique_variants.vcf\n";
}

// Generates a unique key for a variant based on chrom, pos, ref, and alt.
// For multi-allelic ALT fields, the ALT alleles are split, sorted, and rejoined.
std::string generateNormalizedVariantKey(const std::string &chrom, const std::string &pos, const std::string &ref,
                                         const std::string &alt) {
    // Split ALT field on commas.
    std::vector<std::string> alts = splitString(alt, ',');
    // Sort alleles lexicographically.
    std::sort(alts.begin(), alts.end());
    // Rejoin sorted alleles into a single string.
    std::ostringstream oss;
    for (size_t i = 0; i < alts.size(); ++i) {
        if (i > 0) {
            oss << ",";
        }
        oss << alts[i];
    }
    std::string normalizedAlt = oss.str();
    return chrom + ":" + pos + ":" + ref + ":" + normalizedAlt;
}

// Helper function to generate a VariantKey from parsed values.
static VariantKey generateVariantKey(const std::string &chrom, const std::string &pos, const std::string &ref,
                                     const std::string &alt) {
    VariantKey key;
    key.chrom = chrom;
    try {
        key.pos = std::stoi(pos);
    } catch (...) {
        key.pos = 0;
    }
    key.ref = ref;

    // Normalize ALT: split multi-allelic values, sort them, then rejoin.  This
    // avoids parsing the generated key string, which could break for ALT
    // alleles containing ':' such as breakend notation.
    std::vector<std::string> alts = splitString(alt, ',');
    std::sort(alts.begin(), alts.end());
    std::ostringstream oss;
    for (size_t i = 0; i < alts.size(); ++i) {
        if (i > 0) {
            oss << ',';
        }
        oss << alts[i];
    }
    key.alt = oss.str();
    return key;
}

// Function to remove duplicate variants from a VCF file.
bool removeDuplicates(std::istream &in, std::ostream &out) {
    std::string line;
    // Print header lines as-is.
    // Use an unordered_set with our custom hash function to track seen variants.
    std::unordered_set<VariantKey, VariantKeyHash> seen_variants;

    // Performance: reuse vector across iterations
    std::vector<std::string> fields;
    fields.reserve(16);

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '#') {
            // Print header lines unchanged.
            out << line << "\n";
            continue;
        }

        // Performance: use find-based tab splitting instead of stringstream
        vcfx::split_tabs(line, fields);
        // Expect at least 5 columns: CHROM, POS, ID, REF, ALT
        if (fields.size() < 5) {
            std::cerr << "Warning: Skipping invalid VCF line.\n";
            continue;
        }

        // Generate a variant key that normalizes ALT ordering.
        VariantKey key = generateVariantKey(fields[0], fields[1], fields[3], fields[4]);
        if (seen_variants.find(key) == seen_variants.end()) {
            // Variant is unique; record and output the line.
            seen_variants.insert(key);
            out << line << "\n";
        } else {
            // Duplicate variant; skip.
        }
    }
    return true;
}

// ----------------------------------------------------------------------
// main: Parse command-line arguments and call removeDuplicates.
// ----------------------------------------------------------------------
static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_duplicate_remover", show_help))
        return 0;
    // Simple argument parsing: if --help or -h is provided, print help.
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
    }

    // Remove duplicates from standard input and write to standard output.
    bool success = removeDuplicates(std::cin, std::cout);
    return success ? 0 : 1;
}
