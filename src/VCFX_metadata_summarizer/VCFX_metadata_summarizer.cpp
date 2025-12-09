#include "VCFX_metadata_summarizer.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string_view>

// mmap support
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

// ============================================================================
// Helper functions
// ============================================================================

// Fast newline finder using memchr (often SIMD-optimized by libc)
static inline const char* findNewline(const char* start, const char* end) {
    return static_cast<const char*>(memchr(start, '\n', end - start));
}

// Fast tab finder
static inline const char* findTab(const char* start, const char* end) {
    return static_cast<const char*>(memchr(start, '\t', end - start));
}

// Extract ID from a header line like ##contig=<ID=chr1,length=...>
static std::string_view extractIDSV(std::string_view line) {
    // Find "ID="
    auto idPos = line.find("ID=");
    if (idPos == std::string_view::npos) {
        return {};
    }

    auto sub = line.substr(idPos + 3);
    // Find end: first of ',' or '>'
    size_t endPos = 0;
    for (size_t i = 0; i < sub.size(); i++) {
        if (sub[i] == ',' || sub[i] == '>') {
            endPos = i;
            break;
        }
        if (i == sub.size() - 1) {
            endPos = sub.size();
        }
    }
    return sub.substr(0, endPos);
}

// Legacy version for string
static std::string extractID(const std::string &line, const std::string &prefix) {
    auto idPos = line.find("ID=");
    if (idPos == std::string::npos) {
        return "";
    }
    auto sub = line.substr(idPos + 3);
    size_t endPos = sub.find_first_of(",>");
    if (endPos == std::string::npos) {
        return sub;
    }
    return sub.substr(0, endPos);
}

// ============================================================================
// VCFXMetadataSummarizer implementation
// ============================================================================

int VCFXMetadataSummarizer::run(int argc, char *argv[]) {
    int opt;
    bool showHelp = false;
    const char* inputFile = nullptr;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"input", required_argument, 0, 'i'},
        {0, 0, 0, 0}
    };

    optind = 1;
    while ((opt = getopt_long(argc, argv, "hi:", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
            break;
        case 'i':
            inputFile = optarg;
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
        if (!processFileMmap(inputFile)) {
            return 1;
        }
    } else {
        summarizeMetadata(std::cin);
    }

    return 0;
}

void VCFXMetadataSummarizer::displayHelp() {
    std::cout << "VCFX_metadata_summarizer: Summarize key metadata (contigs, INFO, FILTER, FORMAT, samples, variants) "
                 "from a VCF file.\n\n"
              << "Usage:\n"
              << "  VCFX_metadata_summarizer [options] [input.vcf]\n"
              << "  VCFX_metadata_summarizer [options] < input.vcf\n\n"
              << "Options:\n"
              << "  -i, --input FILE   Input VCF file (uses mmap for best performance)\n"
              << "  -h, --help         Display this help message and exit\n\n"
              << "Performance:\n"
              << "  When using -i/--input, the tool uses memory-mapped I/O for\n"
              << "  ~10-20x faster processing of large files.\n\n"
              << "Example:\n"
              << "  VCFX_metadata_summarizer -i input.vcf\n"
              << "  VCFX_metadata_summarizer < input.vcf\n";
}

// ============================================================================
// Optimized mmap-based processing
// ============================================================================

bool VCFXMetadataSummarizer::processFileMmap(const char* filename) {
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
        printSummary();
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

    // Process lines
    while (ptr < dataEnd) {
        const char* lineEnd = findNewline(ptr, dataEnd);
        if (!lineEnd) lineEnd = dataEnd;

        size_t lineLen = lineEnd - ptr;
        if (lineLen == 0) {
            ptr = lineEnd + 1;
            continue;
        }

        if (ptr[0] == '#') {
            // Header line
            if (lineLen >= 2 && ptr[1] == '#') {
                // Meta line ##...
                std::string_view line(ptr, lineLen);
                parseHeaderSV(line);
            } else if (lineLen >= 6 && strncmp(ptr, "#CHROM", 6) == 0) {
                // #CHROM line - count samples
                // Count tabs after position 8 (the 9th column starts samples)
                int tabCount = 0;
                for (size_t i = 0; i < lineLen; i++) {
                    if (ptr[i] == '\t') tabCount++;
                }
                // Fields: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = 9 fixed
                // Sample count = total fields - 9
                // Total fields = tabCount + 1
                numSamples = (tabCount + 1 > 9) ? (tabCount + 1 - 9) : 0;
            }
        } else {
            // Data line - count variant
            numVariants++;
        }

        ptr = lineEnd + 1;
    }

    // Cleanup
    munmap(const_cast<char*>(data), fileSize);
    close(fd);

    // Print summary
    printSummary();

    return true;
}

void VCFXMetadataSummarizer::parseHeaderSV(std::string_view line) {
    // Check type of line and extract ID
    if (line.find("##contig=") != std::string_view::npos) {
        auto idStr = extractIDSV(line);
        if (!idStr.empty()) {
            contigIDs.insert(std::string(idStr));
        }
    } else if (line.find("##INFO=") != std::string_view::npos) {
        auto idStr = extractIDSV(line);
        if (!idStr.empty()) {
            infoIDs.insert(std::string(idStr));
        }
    } else if (line.find("##FILTER=") != std::string_view::npos) {
        auto idStr = extractIDSV(line);
        if (!idStr.empty()) {
            filterIDs.insert(std::string(idStr));
        }
    } else if (line.find("##FORMAT=") != std::string_view::npos) {
        auto idStr = extractIDSV(line);
        if (!idStr.empty()) {
            formatIDs.insert(std::string(idStr));
        }
    }
}

// ============================================================================
// Stdin-based processing (legacy/compatibility)
// ============================================================================

void VCFXMetadataSummarizer::summarizeMetadata(std::istream &in) {
    std::string line;

    std::vector<std::string> fields;
    fields.reserve(16);

    while (true) {
        if (!std::getline(in, line))
            break;
        if (line.empty())
            continue;

        if (line[0] == '#') {
            if (line.rfind("##", 0) == 0) {
                parseHeader(line);
            } else if (line.rfind("#CHROM", 0) == 0) {
                vcfx::split_tabs(line, fields);
                if (fields.size() > 9) {
                    this->numSamples = (int)fields.size() - 9;
                } else {
                    this->numSamples = 0;
                }
            }
        } else {
            numVariants++;
        }
    }

    printSummary();
}

void VCFXMetadataSummarizer::parseHeader(const std::string &line) {
    if (line.find("##contig=") != std::string::npos) {
        auto idStr = extractID(line, "##contig=");
        if (!idStr.empty()) {
            contigIDs.insert(idStr);
        }
    } else if (line.find("##INFO=") != std::string::npos) {
        auto idStr = extractID(line, "##INFO=");
        if (!idStr.empty()) {
            infoIDs.insert(idStr);
        }
    } else if (line.find("##FILTER=") != std::string::npos) {
        auto idStr = extractID(line, "##FILTER=");
        if (!idStr.empty()) {
            filterIDs.insert(idStr);
        }
    } else if (line.find("##FORMAT=") != std::string::npos) {
        auto idStr = extractID(line, "##FORMAT=");
        if (!idStr.empty()) {
            formatIDs.insert(idStr);
        }
    }
}

void VCFXMetadataSummarizer::printSummary() const {
    std::cout << "VCF Metadata Summary:\n";
    std::cout << "---------------------\n";
    std::cout << "Number of unique contigs: " << contigIDs.size() << "\n";
    std::cout << "Number of unique INFO fields: " << infoIDs.size() << "\n";
    std::cout << "Number of unique FILTER fields: " << filterIDs.size() << "\n";
    std::cout << "Number of unique FORMAT fields: " << formatIDs.size() << "\n";
    std::cout << "Number of samples: " << numSamples << "\n";
    std::cout << "Number of variants: " << numVariants << "\n";
}

static void show_help() {
    VCFXMetadataSummarizer obj;
    char arg0[] = "VCFX_metadata_summarizer";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();
    if (vcfx::handle_common_flags(argc, argv, "VCFX_metadata_summarizer", show_help))
        return 0;
    VCFXMetadataSummarizer summarizer;
    return summarizer.run(argc, argv);
}
