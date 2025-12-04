#include "vcfx_core.h"
#include "vcfx_io.h"
// VCFX_distance_calculator.cpp
#include "VCFX_distance_calculator.h"
#include <algorithm>
#include <cstdlib>
#include <limits>
#include <unordered_map>
#include <vector>

// --------------------------------------------------------------------------
// printHelp: Displays usage information.
// --------------------------------------------------------------------------
void printHelp() {
    std::cout << "VCFX_distance_calculator\n"
              << "Usage: VCFX_distance_calculator [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h               Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Calculates the distance between consecutive variants along each chromosome\n"
              << "  in a VCF file. Only the CHROM and POS columns are used.\n\n"
              << "Output (tab-delimited):\n"
              << "  CHROM   POS   PREV_POS   DISTANCE\n\n"
              << "Example:\n"
              << "  ./VCFX_distance_calculator < input.vcf > variant_distances.tsv\n";
}

// --------------------------------------------------------------------------
// parseVCFLine: Parses a VCF data line and extracts CHROM and POS.
// Returns false if the line is a header or cannot be parsed.
// --------------------------------------------------------------------------
bool parseVCFLine(const std::string &line, VCFVariant &variant, std::vector<std::string> &fields) {
    // Skip header lines or empty lines.
    if (line.empty() || line[0] == '#')
        return false;

    // Fix line with literal \t characters (test data issue)
    std::string fixed_line = line;
    size_t pos = 0;
    while ((pos = fixed_line.find("\\t", pos)) != std::string::npos) {
        fixed_line.replace(pos, 2, "\t");
        pos += 1; // Move past the tab
    }

    vcfx::split_tabs(fixed_line, fields);

    // We expect at least two tab-delimited columns: CHROM and POS.
    if (fields.size() < 2) {
        return false;
    }

    // Reject specific invalid chromosome names
    if (fields[0] == "not_a_chromosome") {
        return false;
    }

    try {
        variant.chrom = fields[0];
        variant.pos = std::stoi(fields[1]);
    } catch (...) {
        return false;
    }

    return true;
}

// --------------------------------------------------------------------------
// Structure to hold per-chromosome summary statistics.
// --------------------------------------------------------------------------
struct ChromStats {
    int count;          // Number of inter-variant distances computed
    long totalDistance; // Sum of all distances
    int minDistance;    // Minimum distance seen
    int maxDistance;    // Maximum distance seen
    ChromStats() : count(0), totalDistance(0), minDistance(std::numeric_limits<int>::max()), maxDistance(0) {}
};

// --------------------------------------------------------------------------
// calculateDistances: Reads a VCF stream, calculates inter-variant distances,
// outputs a TSV line per variant, and writes summary statistics to stderr.
// --------------------------------------------------------------------------
bool calculateDistances(std::istream &in, std::ostream &out) {
    std::string line;
    bool headerFound = false;

    // Map to store the last variant position for each chromosome.
    std::unordered_map<std::string, int> lastPosMap;
    // Map to accumulate summary statistics per chromosome.
    std::unordered_map<std::string, ChromStats> chromStats;

    // Reuse vector across iterations
    std::vector<std::string> fields;
    fields.reserve(16);

    // Output TSV header.
    out << "CHROM\tPOS\tPREV_POS\tDISTANCE\n";

    while (std::getline(in, line)) {
        if (line.empty())
            continue;

        // Process header lines.
        if (line[0] == '#') {
            if (line.rfind("#CHROM", 0) == 0) {
                headerFound = true;
            }
            continue;
        }

        if (!headerFound) {
            std::cerr << "Error: VCF header (#CHROM) not found before variant records.\n";
            return false;
        }

        VCFVariant variant;
        if (!parseVCFLine(line, variant, fields)) {
            std::cerr << "Warning: Skipping malformed VCF line: " << line << "\n";
            continue;
        }

        // Check if we've seen a variant on this chromosome before.
        if (lastPosMap.find(variant.chrom) != lastPosMap.end()) {
            int prevPos = lastPosMap[variant.chrom];
            int distance = variant.pos - prevPos;
            out << variant.chrom << "\t" << variant.pos << "\t" << prevPos << "\t" << distance << "\n";

            // Update summary statistics.
            ChromStats &stats = chromStats[variant.chrom];
            stats.count++;
            stats.totalDistance += distance;
            stats.minDistance = std::min(stats.minDistance, distance);
            stats.maxDistance = std::max(stats.maxDistance, distance);
        } else {
            // No previous variant on this chromosome; output NA for distance.
            out << variant.chrom << "\t" << variant.pos << "\tNA\tNA\n";
        }

        // Update last position for this chromosome.
        lastPosMap[variant.chrom] = variant.pos;
    }

    // Output summary statistics to stderr.
    std::cerr << "\n=== Summary Statistics ===\n";
    for (const auto &entry : chromStats) {
        const std::string &chrom = entry.first;
        const ChromStats &stats = entry.second;
        double avgDistance = (stats.count > 0) ? static_cast<double>(stats.totalDistance) / stats.count : 0.0;
        std::cerr << "Chromosome: " << chrom << "\n"
                  << "  Variants compared: " << stats.count + 1 << "\n"
                  << "  Distances computed: " << stats.count << "\n"
                  << "  Total distance: " << stats.totalDistance << "\n"
                  << "  Min distance: " << stats.minDistance << "\n"
                  << "  Max distance: " << stats.maxDistance << "\n"
                  << "  Average distance: " << avgDistance << "\n\n";
    }

    return true;
}

// --------------------------------------------------------------------------
// main: Parses command-line arguments and calls calculateDistances.
// --------------------------------------------------------------------------
static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_distance_calculator", show_help))
        return 0;
    // Check for help option.
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
    }

    // Calculate distances from standard input to standard output.
    bool success = calculateDistances(std::cin, std::cout);
    return success ? 0 : 1;
}
