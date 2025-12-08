#include "VCFX_duplicate_remover.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <sstream>
#include <vector>

// mmap support
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

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

int VCFXDuplicateRemover::run(int argc, char *argv[]) {
    int opt;
    bool showHelp = false;
    const char* inputFile = nullptr;

    // Reset getopt
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
        if (!processFileMmap(inputFile, std::cout)) {
            return 1;
        }
    } else {
        // Use stdin mode
        if (!removeDuplicates(std::cin, std::cout)) {
            return 1;
        }
    }

    return 0;
}

void VCFXDuplicateRemover::displayHelp() {
    std::cout << "VCFX_duplicate_remover: Remove duplicate variants from VCF files.\n\n"
              << "Usage:\n"
              << "  VCFX_duplicate_remover [options] [input.vcf]\n"
              << "  VCFX_duplicate_remover [options] < input.vcf > output.vcf\n\n"
              << "Options:\n"
              << "  -i, --input FILE    Input VCF file (uses mmap for best performance)\n"
              << "  -q, --quiet         Suppress warning messages\n"
              << "  -h, --help          Display this help message and exit\n\n"
              << "Description:\n"
              << "  Removes duplicate variants from a VCF file based on the combination of\n"
              << "  chromosome, position, REF, and ALT alleles. For multi-allelic records, the\n"
              << "  ALT field is normalized by sorting the comma-separated alleles so that the\n"
              << "  ordering does not affect duplicate detection.\n\n"
              << "Performance:\n"
              << "  When using -i/--input, the tool uses memory-mapped I/O for\n"
              << "  ~10-20x faster processing of large files.\n\n"
              << "Example:\n"
              << "  VCFX_duplicate_remover -i input.vcf > unique_variants.vcf\n"
              << "  VCFX_duplicate_remover < input.vcf > unique_variants.vcf\n";
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

    // Normalize ALT: split multi-allelic values, sort them, then rejoin.
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

// Generate VariantKey from raw pointers (mmap mode - avoids some allocations)
static VariantKey generateVariantKeyRaw(const char* chromStart, size_t chromLen,
                                        const char* posStart, size_t posLen,
                                        const char* refStart, size_t refLen,
                                        const char* altStart, size_t altLen) {
    VariantKey key;
    key.chrom.assign(chromStart, chromLen);

    // Parse position
    char* endPtr;
    key.pos = static_cast<int>(strtol(posStart, &endPtr, 10));

    key.ref.assign(refStart, refLen);

    // Normalize ALT: split, sort, rejoin
    std::vector<std::string> alts;
    alts.reserve(4);

    const char* p = altStart;
    const char* altEnd = altStart + altLen;
    const char* tokenStart = p;

    while (p <= altEnd) {
        if (p == altEnd || *p == ',') {
            if (p > tokenStart) {
                alts.emplace_back(tokenStart, p - tokenStart);
            }
            tokenStart = p + 1;
        }
        p++;
    }

    std::sort(alts.begin(), alts.end());

    // Rejoin
    key.alt.clear();
    for (size_t i = 0; i < alts.size(); ++i) {
        if (i > 0) {
            key.alt += ',';
        }
        key.alt += alts[i];
    }

    return key;
}

bool VCFXDuplicateRemover::removeDuplicates(std::istream &in, std::ostream &out) {
    std::string line;
    std::unordered_set<VariantKey, VariantKeyHash> seen_variants;

    // Performance: reuse vector across iterations
    std::vector<std::string> fields;
    fields.reserve(16);

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '#') {
            out << line << "\n";
            continue;
        }

        vcfx::split_tabs(line, fields);
        if (fields.size() < 5) {
            if (!quietMode) {
                std::cerr << "Warning: Skipping invalid VCF line.\n";
            }
            continue;
        }

        VariantKey key = generateVariantKey(fields[0], fields[1], fields[3], fields[4]);
        if (seen_variants.find(key) == seen_variants.end()) {
            seen_variants.insert(key);
            out << line << "\n";
        }
    }
    return true;
}

// Helper function to find newline in memory
static inline const char* findNewline(const char* start, const char* end) {
    return static_cast<const char*>(memchr(start, '\n', end - start));
}

// Find tab and return pointer to it (or nullptr)
static inline const char* findTab(const char* start, const char* end) {
    return static_cast<const char*>(memchr(start, '\t', end - start));
}

bool VCFXDuplicateRemover::processFileMmap(const char* filename, std::ostream &out) {
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

    // Hash set for deduplication
    std::unordered_set<VariantKey, VariantKeyHash> seen_variants;
    seen_variants.reserve(500000);  // Pre-allocate for typical VCF sizes

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

            if (outputBuffer.size() >= flushThreshold) {
                out.write(outputBuffer.data(), outputBuffer.size());
                outputBuffer.clear();
            }
            continue;
        }

        // Data line: extract CHROM (0), POS (1), skip ID (2), REF (3), ALT (4)
        const char* lineStart = ptr;

        // Field 0: CHROM
        const char* chromStart = ptr;
        const char* tab1 = findTab(ptr, lineEnd);
        if (!tab1) {
            if (!quietMode) {
                std::cerr << "Warning: Skipping invalid VCF line.\n";
            }
            ptr = lineEnd + 1;
            continue;
        }
        size_t chromLen = tab1 - chromStart;

        // Field 1: POS
        const char* posStart = tab1 + 1;
        const char* tab2 = findTab(posStart, lineEnd);
        if (!tab2) {
            if (!quietMode) {
                std::cerr << "Warning: Skipping invalid VCF line.\n";
            }
            ptr = lineEnd + 1;
            continue;
        }
        size_t posLen = tab2 - posStart;

        // Field 2: ID (skip)
        const char* tab3 = findTab(tab2 + 1, lineEnd);
        if (!tab3) {
            if (!quietMode) {
                std::cerr << "Warning: Skipping invalid VCF line.\n";
            }
            ptr = lineEnd + 1;
            continue;
        }

        // Field 3: REF
        const char* refStart = tab3 + 1;
        const char* tab4 = findTab(refStart, lineEnd);
        if (!tab4) {
            if (!quietMode) {
                std::cerr << "Warning: Skipping invalid VCF line.\n";
            }
            ptr = lineEnd + 1;
            continue;
        }
        size_t refLen = tab4 - refStart;

        // Field 4: ALT
        const char* altStart = tab4 + 1;
        const char* tab5 = findTab(altStart, lineEnd);
        size_t altLen = tab5 ? (tab5 - altStart) : (lineEnd - altStart);

        // Generate key and check for duplicates
        VariantKey key = generateVariantKeyRaw(chromStart, chromLen, posStart, posLen,
                                               refStart, refLen, altStart, altLen);

        if (seen_variants.find(key) == seen_variants.end()) {
            seen_variants.insert(std::move(key));
            outputBuffer.append(lineStart, lineLen);
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

// main
static void show_help() {
    VCFXDuplicateRemover obj;
    char arg0[] = "VCFX_duplicate_remover";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_duplicate_remover", show_help))
        return 0;
    VCFXDuplicateRemover remover;
    return remover.run(argc, argv);
}
