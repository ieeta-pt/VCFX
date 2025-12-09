#include "VCFX_distance_calculator.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <limits>
#include <unordered_map>
#include <vector>

// mmap support
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

// ============================================================================
// ChromStats implementation
// ============================================================================

ChromStats::ChromStats()
    : count(0), totalDistance(0),
      minDistance(std::numeric_limits<int>::max()), maxDistance(0) {}

// ============================================================================
// Helper functions
// ============================================================================

// Find newline in memory
static inline const char* findNewline(const char* start, const char* end) {
    return static_cast<const char*>(memchr(start, '\n', end - start));
}

// Find tab in memory
static inline const char* findTab(const char* start, const char* end) {
    return static_cast<const char*>(memchr(start, '\t', end - start));
}

// Fast integer parsing (no error checking, assumes valid input)
static inline int fastParseInt(const char* start, const char* end) {
    int result = 0;
    bool negative = false;
    const char* p = start;

    if (p < end && *p == '-') {
        negative = true;
        p++;
    }

    while (p < end && *p >= '0' && *p <= '9') {
        result = result * 10 + (*p - '0');
        p++;
    }

    return negative ? -result : result;
}

// ============================================================================
// VCFXDistanceCalculator implementation
// ============================================================================

int VCFXDistanceCalculator::run(int argc, char *argv[]) {
    int opt;
    bool showHelp = false;
    const char* inputFile = nullptr;

    optind = 1;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"input", required_argument, 0, 'i'},
        {"quiet", no_argument, 0, 'q'},
        {0, 0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "hi:q", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
            break;
        case 'i':
            inputFile = optarg;
            break;
        case 'q':
            quietMode = true;
            break;
        default:
            showHelp = true;
        }
    }

    // Check for positional argument
    if (!inputFile && optind < argc) {
        inputFile = argv[optind];
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    if (inputFile) {
        if (!processFileMmap(inputFile, std::cout)) {
            return 1;
        }
    } else {
        if (!calculateDistances(std::cin, std::cout)) {
            return 1;
        }
    }

    return 0;
}

void VCFXDistanceCalculator::displayHelp() {
    std::cout << "VCFX_distance_calculator: Calculate distances between consecutive variants.\n\n"
              << "Usage:\n"
              << "  VCFX_distance_calculator [options] [input.vcf]\n"
              << "  VCFX_distance_calculator [options] < input.vcf > output.tsv\n\n"
              << "Options:\n"
              << "  -i, --input FILE    Input VCF file (uses mmap for best performance)\n"
              << "  -q, --quiet         Suppress summary statistics to stderr\n"
              << "  -h, --help          Display this help message and exit\n\n"
              << "Description:\n"
              << "  Calculates the distance between consecutive variants along each chromosome\n"
              << "  in a VCF file. Only the CHROM and POS columns are used.\n\n"
              << "Output (tab-delimited):\n"
              << "  CHROM   POS   PREV_POS   DISTANCE\n\n"
              << "Performance:\n"
              << "  When using -i/--input, the tool uses memory-mapped I/O for\n"
              << "  ~10-15x faster processing of large files.\n\n"
              << "Example:\n"
              << "  VCFX_distance_calculator -i input.vcf > variant_distances.tsv\n"
              << "  VCFX_distance_calculator < input.vcf > variant_distances.tsv\n";
}

void VCFXDistanceCalculator::outputSummary(
    const std::unordered_map<std::string, ChromStats>& chromStats) {
    if (quietMode) return;

    std::cerr << "\n=== Summary Statistics ===\n";
    for (const auto &entry : chromStats) {
        const std::string &chrom = entry.first;
        const ChromStats &stats = entry.second;
        double avgDistance = (stats.count > 0)
            ? static_cast<double>(stats.totalDistance) / stats.count
            : 0.0;
        std::cerr << "Chromosome: " << chrom << "\n"
                  << "  Variants compared: " << stats.count + 1 << "\n"
                  << "  Distances computed: " << stats.count << "\n"
                  << "  Total distance: " << stats.totalDistance << "\n"
                  << "  Min distance: " << stats.minDistance << "\n"
                  << "  Max distance: " << stats.maxDistance << "\n"
                  << "  Average distance: " << avgDistance << "\n\n";
    }
}

bool VCFXDistanceCalculator::processFileMmap(const char* filename, std::ostream &out) {
    // Open file
    int fd = open(filename, O_RDONLY);
    if (fd < 0) {
        std::cerr << "Error: Cannot open file: " << filename << "\n";
        return false;
    }

    // Get file size
    struct stat st;
    if (fstat(fd, &st) < 0) {
        std::cerr << "Error: Cannot stat file: " << filename << "\n";
        close(fd);
        return false;
    }

    size_t fileSize = st.st_size;
    if (fileSize == 0) {
        close(fd);
        return true;
    }

    // Memory map
    const char* data = static_cast<const char*>(mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, fd, 0));
    if (data == MAP_FAILED) {
        std::cerr << "Error: Cannot mmap file: " << filename << "\n";
        close(fd);
        return false;
    }

    madvise(const_cast<char*>(data), fileSize, MADV_SEQUENTIAL);

    const char* ptr = data;
    const char* dataEnd = data + fileSize;

    // Output buffer
    std::string outputBuffer;
    outputBuffer.reserve(1024 * 1024);  // 1MB
    const size_t flushThreshold = 900 * 1024;

    // Tracking state
    std::unordered_map<std::string, int> lastPosMap;
    std::unordered_map<std::string, ChromStats> chromStats;
    lastPosMap.reserve(64);
    chromStats.reserve(64);

    bool headerFound = false;
    char outBuf[256];

    // Output TSV header
    outputBuffer.append("CHROM\tPOS\tPREV_POS\tDISTANCE\n");

    // Process lines
    while (ptr < dataEnd) {
        const char* lineEnd = findNewline(ptr, dataEnd);
        if (!lineEnd) lineEnd = dataEnd;

        size_t lineLen = lineEnd - ptr;
        if (lineLen == 0) {
            ptr = lineEnd + 1;
            continue;
        }

        // Skip header lines
        if (ptr[0] == '#') {
            if (lineLen >= 6 && strncmp(ptr, "#CHROM", 6) == 0) {
                headerFound = true;
            }
            ptr = lineEnd + 1;
            continue;
        }

        if (!headerFound) {
            if (!quietMode) {
                std::cerr << "Error: VCF header (#CHROM) not found before variant records.\n";
            }
            munmap(const_cast<char*>(data), fileSize);
            close(fd);
            return false;
        }

        // Parse CHROM (field 0)
        const char* chromStart = ptr;
        const char* tab1 = findTab(ptr, lineEnd);
        if (!tab1) {
            ptr = lineEnd + 1;
            continue;
        }
        size_t chromLen = tab1 - chromStart;

        // Skip invalid chromosome names
        if (chromLen == 16 && strncmp(chromStart, "not_a_chromosome", 16) == 0) {
            ptr = lineEnd + 1;
            continue;
        }

        // Parse POS (field 1)
        const char* posStart = tab1 + 1;
        const char* tab2 = findTab(posStart, lineEnd);
        const char* posEnd = tab2 ? tab2 : lineEnd;

        // Validate position is numeric
        bool validPos = true;
        for (const char* p = posStart; p < posEnd && validPos; p++) {
            if (*p < '0' || *p > '9') {
                if (p == posStart && *p == '-') continue;  // Allow negative
                validPos = false;
            }
        }
        if (!validPos || posStart == posEnd) {
            ptr = lineEnd + 1;
            continue;
        }

        int pos = fastParseInt(posStart, posEnd);
        std::string chrom(chromStart, chromLen);

        // Check if we've seen this chromosome before
        auto it = lastPosMap.find(chrom);
        if (it != lastPosMap.end()) {
            int prevPos = it->second;
            int distance = pos - prevPos;

            // Format output
            int len = snprintf(outBuf, sizeof(outBuf), "\t%d\t%d\t%d\n",
                               pos, prevPos, distance);
            outputBuffer.append(chrom);
            outputBuffer.append(outBuf, len);

            // Update statistics
            ChromStats &stats = chromStats[chrom];
            stats.count++;
            stats.totalDistance += distance;
            stats.minDistance = std::min(stats.minDistance, distance);
            stats.maxDistance = std::max(stats.maxDistance, distance);

            it->second = pos;
        } else {
            // First variant on this chromosome
            int len = snprintf(outBuf, sizeof(outBuf), "\t%d\tNA\tNA\n", pos);
            outputBuffer.append(chrom);
            outputBuffer.append(outBuf, len);
            lastPosMap[chrom] = pos;
        }

        ptr = lineEnd + 1;

        // Periodic flush
        if (outputBuffer.size() >= flushThreshold) {
            out.write(outputBuffer.data(), outputBuffer.size());
            outputBuffer.clear();
        }
    }

    // Final flush
    if (!outputBuffer.empty()) {
        out.write(outputBuffer.data(), outputBuffer.size());
    }

    // Cleanup
    munmap(const_cast<char*>(data), fileSize);
    close(fd);

    // Output summary
    outputSummary(chromStats);

    return true;
}

bool VCFXDistanceCalculator::calculateDistances(std::istream &in, std::ostream &out) {
    std::string line;
    bool headerFound = false;

    std::unordered_map<std::string, int> lastPosMap;
    std::unordered_map<std::string, ChromStats> chromStats;

    std::vector<std::string> fields;
    fields.reserve(16);

    char outBuf[256];
    std::string outLine;
    outLine.reserve(256);

    // Output TSV header
    out << "CHROM\tPOS\tPREV_POS\tDISTANCE\n";

    while (std::getline(in, line)) {
        if (line.empty())
            continue;

        // Process header lines
        if (line[0] == '#') {
            if (line.rfind("#CHROM", 0) == 0) {
                headerFound = true;
            }
            continue;
        }

        if (!headerFound) {
            if (!quietMode) {
                std::cerr << "Error: VCF header (#CHROM) not found before variant records.\n";
            }
            return false;
        }

        // Handle escaped tabs (rare in real data)
        const std::string *line_ptr = &line;
        std::string fixed_line;
        if (line.find("\\t") != std::string::npos) {
            fixed_line = line;
            size_t pos = 0;
            while ((pos = fixed_line.find("\\t", pos)) != std::string::npos) {
                fixed_line.replace(pos, 2, "\t");
                pos += 1;
            }
            line_ptr = &fixed_line;
        }

        vcfx::split_tabs(*line_ptr, fields);

        if (fields.size() < 2) {
            continue;
        }

        // Skip invalid chromosome names
        if (fields[0] == "not_a_chromosome") {
            continue;
        }

        int pos;
        try {
            pos = std::stoi(fields[1]);
        } catch (...) {
            continue;
        }

        const std::string& chrom = fields[0];

        // Check if we've seen this chromosome before
        if (lastPosMap.find(chrom) != lastPosMap.end()) {
            int prevPos = lastPosMap[chrom];
            int distance = pos - prevPos;

            int len = snprintf(outBuf, sizeof(outBuf), "\t%d\t%d\t%d\n",
                               pos, prevPos, distance);
            outLine.clear();
            outLine.append(chrom);
            outLine.append(outBuf, len);
            out.write(outLine.data(), outLine.size());

            ChromStats &stats = chromStats[chrom];
            stats.count++;
            stats.totalDistance += distance;
            stats.minDistance = std::min(stats.minDistance, distance);
            stats.maxDistance = std::max(stats.maxDistance, distance);
        } else {
            int len = snprintf(outBuf, sizeof(outBuf), "\t%d\tNA\tNA\n", pos);
            outLine.clear();
            outLine.append(chrom);
            outLine.append(outBuf, len);
            out.write(outLine.data(), outLine.size());
        }

        lastPosMap[chrom] = pos;
    }

    outputSummary(chromStats);
    return true;
}

// ============================================================================
// Legacy functions for backward compatibility
// ============================================================================

void printHelp() {
    std::cout << "VCFX_distance_calculator\n"
              << "Usage: VCFX_distance_calculator [OPTIONS]\n\n"
              << "Options:\n"
              << "  -i, --input FILE     Input VCF file (uses mmap for best performance)\n"
              << "  -q, --quiet          Suppress summary statistics\n"
              << "  --help, -h           Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Calculates the distance between consecutive variants along each chromosome\n"
              << "  in a VCF file. Only the CHROM and POS columns are used.\n\n"
              << "Output (tab-delimited):\n"
              << "  CHROM   POS   PREV_POS   DISTANCE\n\n"
              << "Example:\n"
              << "  ./VCFX_distance_calculator -i input.vcf > variant_distances.tsv\n"
              << "  ./VCFX_distance_calculator < input.vcf > variant_distances.tsv\n";
}

bool parseVCFLine(const std::string &line, VCFVariant &variant) {
    if (line.empty() || line[0] == '#')
        return false;

    std::vector<std::string> fields;
    vcfx::split_tabs(line, fields);

    if (fields.size() < 2)
        return false;

    if (fields[0] == "not_a_chromosome")
        return false;

    try {
        variant.chrom = fields[0];
        variant.pos = std::stoi(fields[1]);
    } catch (...) {
        return false;
    }

    return true;
}

bool calculateDistances(std::istream &in, std::ostream &out) {
    VCFXDistanceCalculator calc;
    return calc.calculateDistances(in, out);
}

// ============================================================================
// main
// ============================================================================

static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();
    if (vcfx::handle_common_flags(argc, argv, "VCFX_distance_calculator", show_help))
        return 0;

    VCFXDistanceCalculator calc;
    return calc.run(argc, argv);
}
