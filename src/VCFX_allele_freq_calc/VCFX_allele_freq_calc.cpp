#include "vcfx_core.h"
#include <cctype>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

// ---------------------------------------------------------
// Print help message
// ---------------------------------------------------------
static void printHelp() {
    std::cout << "VCFX_allele_freq_calc\n"
              << "Usage: VCFX_allele_freq_calc [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h   Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Reads a VCF from stdin and outputs a TSV file:\n"
              << "    CHROM  POS  ID  REF  ALT  Allele_Frequency\n\n"
              << "  Allele frequency is computed as (#ALT alleles / total #alleles),\n"
              << "  counting any non-zero numeric allele (1,2,3,...) as ALT.\n\n"
              << "Example:\n"
              << "  ./VCFX_allele_freq_calc < input.vcf > allele_frequencies.tsv\n";
}

// ---------------------------------------------------------
// Utility: fast split using string_view (C++17)
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

// Legacy split for compatibility
static std::vector<std::string> split(const std::string &s, char delimiter) {
    std::vector<std::string> tokens;
    tokens.reserve(16);
    size_t start = 0;
    size_t end;
    while ((end = s.find(delimiter, start)) != std::string::npos) {
        tokens.emplace_back(s, start, end - start);
        start = end + 1;
    }
    if (start <= s.size()) {
        tokens.emplace_back(s, start);
    }
    return tokens;
}

// ---------------------------------------------------------
// parseGenotype: counts how many are alt vs total among numeric alleles
// Zero-allocation version using string_view
// ---------------------------------------------------------
static void parseGenotype(std::string_view genotype, int &alt_count, int &total_count) {
    // e.g. genotype might be "0/1" or "2|3", etc.
    size_t start = 0;
    size_t len = genotype.size();

    while (start < len) {
        // Find next separator (/ or |)
        size_t end = start;
        while (end < len && genotype[end] != '/' && genotype[end] != '|') {
            end++;
        }

        std::string_view allele = genotype.substr(start, end - start);
        start = end + 1;

        if (allele.empty() || allele == ".") {
            continue;
        }

        // Fast numeric check and value extraction
        bool isZero = true;
        bool numeric = true;
        for (char c : allele) {
            if (!std::isdigit(static_cast<unsigned char>(c))) {
                numeric = false;
                break;
            }
            if (c != '0') {
                isZero = false;
            }
        }

        if (!numeric) {
            continue;
        }

        total_count++;
        if (!isZero) {
            alt_count++;
        }
    }
}

// ---------------------------------------------------------
// Helper: find nth tab position in a line (0-indexed)
// Returns position after the tab, or npos if not found
// ---------------------------------------------------------
static size_t findNthTab(const std::string &line, int n) {
    size_t pos = 0;
    for (int i = 0; i < n && pos != std::string::npos; ++i) {
        pos = line.find('\t', pos);
        if (pos != std::string::npos) pos++;
    }
    return pos;
}

// ---------------------------------------------------------
// Helper: extract field at given index using string_view
// ---------------------------------------------------------
static std::string_view getField(const std::string &line, size_t start, size_t &next) {
    size_t end = line.find('\t', start);
    if (end == std::string::npos) {
        next = std::string::npos;
        return std::string_view(line.data() + start, line.size() - start);
    }
    next = end + 1;
    return std::string_view(line.data() + start, end - start);
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

// ---------------------------------------------------------
// Main calculation: read VCF from 'in', write results to 'out'
// Optimized for minimal allocations
// ---------------------------------------------------------
static void calculateAlleleFrequency(std::istream &in, std::ostream &out) {
    // Use larger buffer for faster I/O
    constexpr size_t BUFFER_SIZE = 1 << 20; // 1MB buffer
    std::vector<char> inBuffer(BUFFER_SIZE);
    std::vector<char> outBuffer(BUFFER_SIZE);
    in.rdbuf()->pubsetbuf(inBuffer.data(), BUFFER_SIZE);
    out.rdbuf()->pubsetbuf(outBuffer.data(), BUFFER_SIZE);

    std::string line;
    line.reserve(65536); // Reserve space for large VCF lines with many samples

    bool foundChromHeader = false;

    // FORMAT field caching - avoid re-parsing same format
    std::string cachedFormat;
    int cachedGtIndex = -1;

    // Print a single TSV header for our results
    out << "CHROM\tPOS\tID\tREF\tALT\tAllele_Frequency\n";

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            if (line.size() >= 6 && line.compare(0, 6, "#CHROM") == 0) {
                foundChromHeader = true;
            }
            continue;
        }

        if (!foundChromHeader) {
            std::cerr << "Warning: Data line encountered before #CHROM header. Skipping.\n";
            continue;
        }

        // Parse fields using string_view for efficiency
        std::string_view lineView(line);
        auto fields = split_sv(lineView, '\t');

        if (fields.size() < 9) {
            std::cerr << "Warning: Skipping invalid VCF line (fewer than 9 fields).\n";
            continue;
        }

        // Standard VCF columns via string_view (no allocation)
        std::string_view chrom = fields[0];
        std::string_view pos = fields[1];
        std::string_view id = fields[2];
        std::string_view ref = fields[3];
        std::string_view alt = fields[4];
        std::string_view fmt = fields[8];

        // Use cached GT index if FORMAT hasn't changed (common case)
        int gtIndex;
        if (fmt == cachedFormat) {
            gtIndex = cachedGtIndex;
        } else {
            gtIndex = -1;
            size_t start = 0;
            int idx = 0;
            while (start < fmt.size()) {
                size_t end = fmt.find(':', start);
                std::string_view field = (end == std::string_view::npos)
                    ? fmt.substr(start)
                    : fmt.substr(start, end - start);
                if (field == "GT") {
                    gtIndex = idx;
                    break;
                }
                idx++;
                if (end == std::string_view::npos) break;
                start = end + 1;
            }
            cachedFormat = std::string(fmt);
            cachedGtIndex = gtIndex;
        }

        if (gtIndex < 0) {
            continue;
        }

        int alt_count = 0;
        int total_count = 0;

        // Process each sample column (starting at index 9)
        for (size_t i = 9; i < fields.size(); ++i) {
            std::string_view gt = extractGT(fields[i], gtIndex);
            if (!gt.empty()) {
                parseGenotype(gt, alt_count, total_count);
            }
        }

        // Compute allele frequency
        double freq = (total_count > 0)
            ? static_cast<double>(alt_count) / static_cast<double>(total_count)
            : 0.0;

        // Print the row - need to convert string_view to output
        out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\t"
            << std::fixed << std::setprecision(4) << freq << "\n";
    }
}

// ---------------------------------------------------------
// main()
// ---------------------------------------------------------
static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    // Disable stdio synchronization for performance
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);

    if (vcfx::handle_common_flags(argc, argv, "VCFX_allele_freq_calc", show_help))
        return 0;
    // Parse arguments for help
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
    }

    // If no arguments (and no piped input), just show help
    // (This is optional; if the user pipes data, we can proceed.)
    if (argc == 1) {
        // Not strictly required, but user-friendly
        // Check if there's something on stdin, or just show help
        if (std::cin.peek() == EOF) {
            printHelp();
            return 1;
        }
    }

    // Perform allele frequency calculation
    calculateAlleleFrequency(std::cin, std::cout);
    return 0;
}
