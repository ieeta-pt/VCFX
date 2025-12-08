#include "VCFX_af_subsetter.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <iostream>
#include <sstream>

// mmap support
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

// SIMD detection for ARM
#if defined(__aarch64__) || defined(__ARM_NEON)
#include <arm_neon.h>
#define HAS_NEON 1
#endif

int VCFXAfSubsetter::run(int argc, char *argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    double minAF = 0.0;
    double maxAF = 1.0;
    const char* inputFile = nullptr;

    // Reset getopt
    optind = 1;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"af-filter", required_argument, 0, 'a'},
        {"input", required_argument, 0, 'i'},
        {"quiet", no_argument, 0, 'q'},
        {0, 0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "ha:i:q", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
            break;
        case 'a': {
            std::string range(optarg);
            size_t dashPos = range.find('-');
            if (dashPos == std::string::npos) {
                std::cerr << "Error: Invalid AF range format. Use <minAF>-<maxAF>.\n";
                displayHelp();
                return 1;
            }
            try {
                minAF = std::stod(range.substr(0, dashPos));
                maxAF = std::stod(range.substr(dashPos + 1));
                if (minAF < 0.0 || maxAF > 1.0 || minAF > maxAF) {
                    throw std::invalid_argument("AF values out of range or minAF > maxAF.");
                }
            } catch (const std::invalid_argument &) {
                std::cerr << "Error: Invalid AF range values. Ensure they are numbers between 0.0 and 1.0 with minAF "
                             "<= maxAF.\n";
                displayHelp();
                return 1;
            }
            break;
        }
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

    // Check for positional argument (input file)
    if (!inputFile && optind < argc) {
        inputFile = argv[optind];
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    // Process input
    if (inputFile) {
        // Use mmap for file input
        if (!processFileMmap(inputFile, std::cout, minAF, maxAF)) {
            return 1;
        }
    } else {
        // Use stdin mode
        subsetByAlleleFrequency(std::cin, std::cout, minAF, maxAF);
    }

    return 0;
}

void VCFXAfSubsetter::displayHelp() {
    std::cout << "VCFX_af_subsetter: Subset variants based on alternate allele frequency (AF) ranges.\n\n"
              << "Usage:\n"
              << "  VCFX_af_subsetter [options] [input.vcf]\n"
              << "  VCFX_af_subsetter [options] < input.vcf > output.vcf\n\n"
              << "Options:\n"
              << "  -i, --input FILE               Input VCF file (uses mmap for best performance)\n"
              << "  -a, --af-filter <minAF>-<maxAF>  Specify the AF range for filtering (e.g., 0.01-0.05)\n"
              << "  -q, --quiet                    Suppress warning messages\n"
              << "  -h, --help                     Display this help message and exit\n\n"
              << "Performance:\n"
              << "  When using -i/--input, the tool uses memory-mapped I/O for\n"
              << "  ~10x faster processing of large files.\n\n"
              << "Example:\n"
              << "  VCFX_af_subsetter -i input.vcf --af-filter 0.01-0.05 > subsetted.vcf\n"
              << "  VCFX_af_subsetter --af-filter 0.01-0.05 < input.vcf > subsetted.vcf\n";
}

bool VCFXAfSubsetter::parseAF(const std::string &infoField, std::vector<double> &afValues) {
    // Find "AF=" in the INFO string
    size_t pos = infoField.find("AF=");
    if (pos == std::string::npos) {
        return false;
    }

    // Extract substring up to the next semicolon or end of string
    size_t start = pos + 3; // skip 'AF='
    size_t end = infoField.find(';', start);
    std::string afStr = (end == std::string::npos) ? infoField.substr(start) : infoField.substr(start, end - start);

    // Split by comma (to handle multi-allelic AF like "0.2,0.8")
    std::stringstream ss(afStr);
    std::string token;
    while (std::getline(ss, token, ',')) {
        try {
            double afVal = std::stod(token);
            afValues.push_back(afVal);
        } catch (...) {
            return false; // parsing error
        }
    }
    return !afValues.empty();
}

// Parse AF directly from raw pointer (for mmap mode - avoids string allocation)
bool VCFXAfSubsetter::parseAFRaw(const char* infoStart, size_t infoLen, std::vector<double> &afValues) {
    // Find "AF=" in the INFO field using memmem-like search
    const char* p = infoStart;
    const char* end = infoStart + infoLen;

    while (p + 3 <= end) {
        if (p[0] == 'A' && p[1] == 'F' && p[2] == '=') {
            p += 3;  // skip "AF="

            // Parse comma-separated values until semicolon or end
            while (p < end) {
                // Parse the number
                char* numEnd;
                double val = strtod(p, &numEnd);
                if (numEnd == p) {
                    return false;  // parsing error
                }
                afValues.push_back(val);
                p = numEnd;

                // Check what follows
                if (p >= end || *p == ';' || *p == '\t' || *p == '\n') {
                    break;  // end of AF values
                }
                if (*p == ',') {
                    p++;  // skip comma, continue parsing
                } else {
                    break;  // unexpected character
                }
            }
            return !afValues.empty();
        }
        // Move to next potential match (after semicolon or at next position)
        const char* semicolon = static_cast<const char*>(memchr(p, ';', end - p));
        if (semicolon) {
            p = semicolon + 1;
        } else {
            break;
        }
    }
    return false;  // AF= not found
}

void VCFXAfSubsetter::subsetByAlleleFrequency(std::istream &in, std::ostream &out, double minAF, double maxAF) {
    std::string line;

    // Performance: reuse vector across iterations
    std::vector<std::string> fields;
    fields.reserve(16);
    std::vector<double> afValues;
    afValues.reserve(8);

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        // Print header lines (starting with '#') as-is
        if (line[0] == '#') {
            out << line << "\n";
            continue;
        }

        // Performance: use find-based tab splitting instead of stringstream
        vcfx::split_tabs(line, fields);

        if (fields.size() < 8) {
            if (!quietMode) {
                std::cerr << "Warning: Skipping invalid VCF line (less than 8 fields): " << line << "\n";
            }
            continue;
        }

        // Parse AF values from the INFO field
        afValues.clear();
        if (!parseAF(fields[7], afValues)) {
            if (!quietMode) {
                std::cerr << "Warning: AF not found or invalid in INFO field. Skipping variant: " << line << "\n";
            }
            continue;
        }

        // If any AF is within [minAF, maxAF], keep this variant
        bool keepVariant = false;
        for (double af : afValues) {
            if (af >= minAF && af <= maxAF) {
                keepVariant = true;
                break;
            }
        }

        if (keepVariant) {
            out << line << "\n";
        }
    }
}

// Helper function to find newline in memory (uses memchr which is often SIMD-optimized)
static inline const char* findNewline(const char* start, const char* end) {
    return static_cast<const char*>(memchr(start, '\n', end - start));
}

// Find the n-th tab in the line, return pointer to character after tab
static inline const char* findNthTab(const char* start, const char* end, int n) {
    const char* p = start;
    for (int i = 0; i < n && p < end; i++) {
        p = static_cast<const char*>(memchr(p, '\t', end - p));
        if (!p) return nullptr;
        p++;  // move past the tab
    }
    return p;
}

bool VCFXAfSubsetter::processFileMmap(const char* filename, std::ostream &out, double minAF, double maxAF) {
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
        return true;  // Empty file is valid
    }

    // Memory map the file
    const char* data = static_cast<const char*>(mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, fd, 0));
    if (data == MAP_FAILED) {
        std::cerr << "Error: Cannot mmap file: " << filename << "\n";
        close(fd);
        return false;
    }

    // Hint to kernel: sequential access
    madvise(const_cast<char*>(data), fileSize, MADV_SEQUENTIAL);

    const char* ptr = data;
    const char* dataEnd = data + fileSize;

    // Output buffer for better performance
    std::string outputBuffer;
    outputBuffer.reserve(1024 * 1024);  // 1MB buffer
    const size_t flushThreshold = 900 * 1024;  // Flush at 900KB

    // Reusable vector for AF parsing
    std::vector<double> afValues;
    afValues.reserve(8);

    while (ptr < dataEnd) {
        // Find end of line
        const char* lineEnd = findNewline(ptr, dataEnd);
        if (!lineEnd) {
            lineEnd = dataEnd;
        }

        size_t lineLen = lineEnd - ptr;

        // Skip empty lines
        if (lineLen == 0) {
            ptr = lineEnd + 1;
            continue;
        }

        // Handle header lines
        if (ptr[0] == '#') {
            outputBuffer.append(ptr, lineLen);
            outputBuffer.push_back('\n');
            ptr = lineEnd + 1;

            // Periodic flush
            if (outputBuffer.size() >= flushThreshold) {
                out.write(outputBuffer.data(), outputBuffer.size());
                outputBuffer.clear();
            }
            continue;
        }

        // Data line: find INFO field (field 7, 0-indexed)
        // Need to skip 7 tabs to get to INFO
        const char* infoStart = findNthTab(ptr, lineEnd, 7);
        if (!infoStart) {
            if (!quietMode) {
                std::cerr << "Warning: Skipping invalid VCF line (less than 8 fields)\n";
            }
            ptr = lineEnd + 1;
            continue;
        }

        // Find end of INFO field (next tab or end of line)
        const char* infoEnd = static_cast<const char*>(memchr(infoStart, '\t', lineEnd - infoStart));
        if (!infoEnd) {
            infoEnd = lineEnd;
        }
        size_t infoLen = infoEnd - infoStart;

        // Parse AF values
        afValues.clear();
        if (!parseAFRaw(infoStart, infoLen, afValues)) {
            if (!quietMode) {
                std::cerr << "Warning: AF not found or invalid in INFO field. Skipping variant.\n";
            }
            ptr = lineEnd + 1;
            continue;
        }

        // Check if any AF is within range
        bool keepVariant = false;
        for (double af : afValues) {
            if (af >= minAF && af <= maxAF) {
                keepVariant = true;
                break;
            }
        }

        if (keepVariant) {
            outputBuffer.append(ptr, lineLen);
            outputBuffer.push_back('\n');
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

    return true;
}

//
// Typical main():
//
static void show_help() {
    VCFXAfSubsetter obj;
    char arg0[] = "VCFX_af_subsetter";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_af_subsetter", show_help))
        return 0;
    VCFXAfSubsetter afSubsetter;
    return afSubsetter.run(argc, argv);
}
