#include "VCFX_missing_detector.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

/**
 * @brief Main logic runner
 */
int VCFXMissingDetector::run(int argc, char *argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;

    static struct option long_options[] = {{"help", no_argument, 0, 'h'}, {0, 0, 0, 0}};

    while ((opt = getopt_long(argc, argv, "h", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
            break;
        default:
            showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    // Perform missing genotype detection
    detectMissingGenotypes(std::cin, std::cout);
    return 0;
}

/**
 * @brief Prints usage
 */
void VCFXMissingDetector::displayHelp() {
    std::cout << "VCFX_missing_detector: Detect variants with missing sample genotypes.\n\n"
              << "Usage:\n"
              << "  VCFX_missing_detector [options] < input.vcf > flagged.vcf\n\n"
              << "Options:\n"
              << "  -h, --help    Display this help message and exit\n\n"
              << "Description:\n"
              << "  Reads a VCF from stdin, checks each sample's genotype for missing data,\n"
              << "  and if any sample has a missing genotype, appends 'MISSING_GENOTYPES=1'\n"
              << "  in the INFO field.\n\n"
              << "Example:\n"
              << "  VCFX_missing_detector < input.vcf > flagged_output.vcf\n";
}

/**
 * @brief Helper function to check if a genotype string has missing data
 * Optimized version using string_view to avoid allocations
 *
 * We consider a genotype missing if:
 *   - The entire string is '.' or './.' or '.|.'
 *   - Either allele is '.' => e.g. '0/.', './1'
 */
static bool genotypeIsMissing(std::string_view gt_field) {
    if (gt_field.empty())
        return true;

    // Extract just the genotype part (before the first colon) - no allocation needed
    size_t colon_pos = gt_field.find(':');
    std::string_view gt = (colon_pos != std::string_view::npos)
        ? gt_field.substr(0, colon_pos)
        : gt_field;

    if (gt == "." || gt == "./." || gt == ".|.")
        return true;

    // Find separator (/ or |) - check for missing alleles directly
    for (size_t i = 0; i < gt.size(); ++i) {
        char c = gt[i];
        if (c == '/' || c == '|') {
            // Check allele before separator
            if (i == 0 || gt[i - 1] == '.') {
                return true;
            }
            // Check allele after separator
            if (i + 1 >= gt.size() || gt[i + 1] == '.') {
                return true;
            }
        }
    }

    // Check for single "." allele (haploid)
    if (gt == ".") {
        return true;
    }

    return false;
}

// ---------------------------------------------------------
// Helper: split into string_views efficiently
// ---------------------------------------------------------
static std::vector<std::string_view> split_sv(std::string_view s, char delimiter) {
    std::vector<std::string_view> tokens;
    tokens.reserve(16);
    size_t start = 0;
    size_t end;
    while ((end = s.find(delimiter, start)) != std::string_view::npos) {
        tokens.push_back(s.substr(start, end - start));
        start = end + 1;
    }
    if (start <= s.size()) {
        tokens.push_back(s.substr(start));
    }
    return tokens;
}

// ---------------------------------------------------------
// Helper: extract GT from sample field at given gtIndex
// ---------------------------------------------------------
static std::string_view extractGT(std::string_view sample, int gtIndex) {
    size_t start = 0;
    for (int i = 0; i < gtIndex && start < sample.size(); ++i) {
        size_t colon = sample.find(':', start);
        if (colon == std::string_view::npos) {
            return std::string_view();
        }
        start = colon + 1;
    }
    size_t end = sample.find(':', start);
    if (end == std::string_view::npos) {
        return sample.substr(start);
    }
    return sample.substr(start, end - start);
}

/**
 * @brief The main function to detect missing genotypes.
 *        If any sample genotype is missing, we append 'MISSING_GENOTYPES=1' to the INFO field.
 *        Optimized for minimal allocations using string_view.
 */
void VCFXMissingDetector::detectMissingGenotypes(std::istream &in, std::ostream &out) {
    // Use larger buffer for faster I/O
    constexpr size_t BUFFER_SIZE = 1 << 20; // 1MB buffer
    std::vector<char> inBuffer(BUFFER_SIZE);
    std::vector<char> outBuffer(BUFFER_SIZE);
    in.rdbuf()->pubsetbuf(inBuffer.data(), BUFFER_SIZE);
    out.rdbuf()->pubsetbuf(outBuffer.data(), BUFFER_SIZE);

    std::string line;
    line.reserve(65536); // Reserve space for large VCF lines

    // FORMAT field caching - avoid re-parsing same format
    std::string cachedFormat;
    int cachedGtIndex = -1;

    while (std::getline(in, line)) {
        if (line.empty()) {
            out << "\n";
            continue;
        }

        if (line[0] == '#') {
            out << line << "\n";
            continue;
        }

        // Use string_view for parsing (no allocation)
        std::string_view lineView(line);
        auto fields = split_sv(lineView, '\t');

        if (fields.size() < 9) {
            out << line << "\n";
            continue;
        }

        // Find GT index in FORMAT field - use cache
        std::string_view format = fields[8];
        int gtIndex;
        if (format == cachedFormat) {
            gtIndex = cachedGtIndex;
        } else {
            gtIndex = -1;
            size_t start = 0;
            int idx = 0;
            while (start < format.size()) {
                size_t end = format.find(':', start);
                std::string_view field = (end == std::string_view::npos)
                    ? format.substr(start)
                    : format.substr(start, end - start);
                if (field == "GT") {
                    gtIndex = idx;
                    break;
                }
                idx++;
                if (end == std::string_view::npos) break;
                start = end + 1;
            }
            cachedFormat = std::string(format);
            cachedGtIndex = gtIndex;
        }

        // Check sample columns for missing genotypes
        bool hasMissing = false;
        for (size_t s = 9; s < fields.size(); ++s) {
            std::string_view gt;
            if (format == "GT") {
                gt = fields[s];
            } else if (gtIndex >= 0) {
                gt = extractGT(fields[s], gtIndex);
            } else {
                continue;
            }

            if (genotypeIsMissing(gt)) {
                hasMissing = true;
                break;
            }
        }

        if (!hasMissing) {
            out << line << "\n";
            continue;
        }

        // Need to modify INFO field - find its position in the original line
        // We'll reconstruct only if there's missing data
        size_t info_start = 0;
        size_t info_end = 0;
        {
            size_t pos = 0;
            for (int i = 0; i < 7 && pos != std::string::npos; ++i) {
                pos = line.find('\t', pos);
                if (pos != std::string::npos) pos++;
            }
            info_start = pos;
            info_end = line.find('\t', info_start);
            if (info_end == std::string::npos) info_end = line.size();
        }

        std::string_view info = fields[7];

        // Output the modified line
        out.write(line.data(), info_start);
        if (info == "." || info.empty()) {
            out << "MISSING_GENOTYPES=1";
        } else {
            out << info;
            if (info.back() != ';') {
                out << ';';
            }
            out << "MISSING_GENOTYPES=1";
        }
        out.write(line.data() + info_end, line.size() - info_end);
        out << "\n";
    }
}

static void show_help() {
    VCFXMissingDetector obj;
    char arg0[] = "VCFX_missing_detector";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio

    if (vcfx::handle_common_flags(argc, argv, "VCFX_missing_detector", show_help))
        return 0;
    VCFXMissingDetector missingDetector;
    return missingDetector.run(argc, argv);
}
