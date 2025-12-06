#include "VCFX_indexer.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

// SIMD includes for x86_64 platforms
#if defined(__x86_64__) && defined(__AVX2__)
#include <immintrin.h>
#define USE_AVX2 1
#elif defined(__x86_64__) && defined(__SSE2__)
#include <emmintrin.h>
#define USE_SSE2 1
#endif

// ============================================================================
// OPTIMIZED HELPER FUNCTIONS
// ============================================================================

#if defined(USE_AVX2)
// Find next newline using AVX2 (32 bytes at a time)
static inline const char* findNewlineSIMD(const char* p, const char* end) {
    const __m256i nl = _mm256_set1_epi8('\n');
    while (p + 32 <= end) {
        __m256i chunk = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p));
        __m256i cmp = _mm256_cmpeq_epi8(chunk, nl);
        int mask = _mm256_movemask_epi8(cmp);
        if (mask) {
            return p + __builtin_ctz(static_cast<unsigned int>(mask));
        }
        p += 32;
    }
    return static_cast<const char*>(memchr(p, '\n', static_cast<size_t>(end - p)));
}

#elif defined(USE_SSE2)
// Find next newline using SSE2 (16 bytes at a time)
static inline const char* findNewlineSIMD(const char* p, const char* end) {
    const __m128i nl = _mm_set1_epi8('\n');
    while (p + 16 <= end) {
        __m128i chunk = _mm_loadu_si128(reinterpret_cast<const __m128i*>(p));
        __m128i cmp = _mm_cmpeq_epi8(chunk, nl);
        int mask = _mm_movemask_epi8(cmp);
        if (mask) {
            return p + __builtin_ctz(static_cast<unsigned int>(mask));
        }
        p += 16;
    }
    return static_cast<const char*>(memchr(p, '\n', static_cast<size_t>(end - p)));
}

#else
// Portable fallback using memchr (still SIMD-optimized in libc)
static inline const char* findNewlineSIMD(const char* p, const char* end) {
    return static_cast<const char*>(memchr(p, '\n', static_cast<size_t>(end - p)));
}
#endif

// Extract CHROM and POS fields without any string allocation
// Returns true if valid data line, false for headers/empty/invalid
static inline bool extractChromPos(const char* line, size_t lineLen,
                                   const char*& chromStart, size_t& chromLen,
                                   int64_t& pos) {
    if (lineLen == 0) return false;

    const char* p = line;
    const char* end = line + lineLen;

    // Skip leading whitespace
    while (p < end && (*p == ' ' || *p == '\t')) p++;
    if (p >= end) return false;

    // Check for header line
    if (*p == '#') return false;

    // CHROM is first field
    chromStart = p;
    while (p < end && *p != '\t') p++;
    chromLen = static_cast<size_t>(p - chromStart);
    if (chromLen == 0 || p >= end) return false;

    p++;  // Skip tab

    // POS is second field - parse directly without stoll
    pos = 0;
    while (p < end && *p >= '0' && *p <= '9') {
        pos = pos * 10 + (*p - '0');
        p++;
    }

    return pos > 0;  // Valid if we parsed at least one digit
}

// Check if line starts with "#CHROM" (for header detection)
static inline bool isChromHeaderLine(const char* line, size_t lineLen) {
    if (lineLen < 6) return false;

    // Skip leading whitespace
    const char* p = line;
    const char* end = line + lineLen;
    while (p < end && (*p == ' ' || *p == '\t')) p++;

    size_t remaining = static_cast<size_t>(end - p);
    if (remaining < 6) return false;

    return (p[0] == '#' && p[1] == 'C' && p[2] == 'H' &&
            p[3] == 'R' && p[4] == 'O' && p[5] == 'M');
}

// ============================================================================
// LEGACY HELPERS (for stdin fallback)
// ============================================================================

// A helper to split a string by literal '\t' characters
static std::vector<std::string> splitTabs(const std::string &s) {
    std::vector<std::string> tokens;
    size_t start = 0;
    while (true) {
        size_t pos = s.find('\t', start);
        if (pos == std::string::npos) {
            tokens.push_back(s.substr(start));
            break;
        }
        tokens.push_back(s.substr(start, pos - start));
        start = pos + 1;
    }
    return tokens;
}

// Trim leading whitespace
static std::string ltrim(const std::string &str) {
    size_t startPos = 0;
    while (startPos < str.size() && std::isspace(static_cast<unsigned char>(str[startPos]))) {
        startPos++;
    }
    return str.substr(startPos);
}

// ============================================================================
// VCFX_INDEXER IMPLEMENTATION
// ============================================================================

void VCFXIndexer::displayHelp() {
    std::cout << "VCFX_indexer\n"
              << "Usage: VCFX_indexer [options] [input.vcf]\n"
              << "       VCFX_indexer [options] < input.vcf\n\n"
              << "Description:\n"
              << "  Reads a VCF from file argument or stdin and writes a 3-column index\n"
              << "  (CHROM, POS, FILE_OFFSET) to stdout. FILE_OFFSET is the byte offset\n"
              << "  from the start of the file to the beginning of each variant line.\n"
              << "  When a file is provided directly, uses memory-mapped I/O for faster processing.\n\n"
              << "Options:\n"
              << "  -h, --help    Show this help message\n\n"
              << "Example:\n"
              << "  VCFX_indexer input.vcf > index.tsv       # Fast memory-mapped mode\n"
              << "  VCFX_indexer < input.vcf > index.tsv     # Stdin mode\n";
}

int VCFXIndexer::run(int argc, char *argv[]) {
    static struct option long_opts[] = {{"help", no_argument, 0, 'h'}, {0, 0, 0, 0}};

    // Reset getopt
    optind = 1;

    while (true) {
        int c = getopt_long(argc, argv, "h", long_opts, nullptr);
        if (c == -1)
            break;
        switch (c) {
        case 'h':
            displayHelp();
            return 0;
        default:
            displayHelp();
            return 1;
        }
    }

    // Check for file argument - use fast mmap path
    if (optind < argc) {
        const char *filename = argv[optind];
        return createVCFIndexMmap(filename, std::cout);
    }

    // Fallback to stdin processing
    createVCFIndex(std::cin, std::cout);
    return 0;
}

// ============================================================================
// OPTIMIZED MMAP-BASED INDEXING
// ============================================================================

int VCFXIndexer::createVCFIndexMmap(const char *filename, std::ostream &out) {
    int fd = open(filename, O_RDONLY);
    if (fd < 0) {
        std::cerr << "Error: cannot open file: " << filename << "\n";
        return 1;
    }

    struct stat st;
    if (fstat(fd, &st) < 0) {
        close(fd);
        std::cerr << "Error: cannot stat file: " << filename << "\n";
        return 1;
    }

    size_t fileSize = static_cast<size_t>(st.st_size);
    if (fileSize == 0) {
        close(fd);
        return 0;  // Empty file
    }

    void *mapped = mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, fd, 0);
    if (mapped == MAP_FAILED) {
        close(fd);
        std::cerr << "Error: cannot mmap file: " << filename << "\n";
        return 1;
    }

    // Advise kernel we'll read sequentially
    madvise(mapped, fileSize, MADV_SEQUENTIAL);

    const char *data = static_cast<const char *>(mapped);
    const char *end = data + fileSize;

    // Output buffer (1MB with 512KB flush threshold)
    std::string outputBuffer;
    outputBuffer.reserve(1024 * 1024);

    bool foundChromHeader = false;
    bool warnedNoChromYet = false;
    bool sawAnyHeaderLine = false;

    const char *p = data;
    while (p < end) {
        // Find end of line using SIMD-optimized search
        const char *lineEnd = findNewlineSIMD(p, end);
        if (!lineEnd) {
            lineEnd = end;  // Last line without trailing newline
        }

        int64_t lineOffset = static_cast<int64_t>(p - data);
        size_t lineLen = static_cast<size_t>(lineEnd - p);

        // Handle Windows line endings
        if (lineLen > 0 && p[lineLen - 1] == '\r') {
            lineLen--;
        }

        // Skip empty lines
        if (lineLen > 0) {
            // Check for header lines
            const char *lineStart = p;
            while (lineStart < p + lineLen && (*lineStart == ' ' || *lineStart == '\t')) {
                lineStart++;
            }

            if (lineStart < p + lineLen && *lineStart == '#') {
                sawAnyHeaderLine = true;
                // Check for #CHROM header
                if (!foundChromHeader && isChromHeaderLine(p, lineLen)) {
                    foundChromHeader = true;
                    outputBuffer.append("CHROM\tPOS\tFILE_OFFSET\n");
                }
            } else if (foundChromHeader) {
                // Data line - extract CHROM and POS
                const char *chromStart;
                size_t chromLen;
                int64_t pos;

                if (extractChromPos(p, lineLen, chromStart, chromLen, pos)) {
                    // Append to output buffer
                    outputBuffer.append(chromStart, chromLen);
                    outputBuffer.push_back('\t');

                    // Fast integer to string
                    char numBuf[32];
                    int numLen = snprintf(numBuf, sizeof(numBuf), "%" PRId64, pos);
                    outputBuffer.append(numBuf, static_cast<size_t>(numLen));
                    outputBuffer.push_back('\t');

                    numLen = snprintf(numBuf, sizeof(numBuf), "%" PRId64, lineOffset);
                    outputBuffer.append(numBuf, static_cast<size_t>(numLen));
                    outputBuffer.push_back('\n');

                    // Flush buffer at 512KB threshold
                    if (outputBuffer.size() > 512 * 1024) {
                        out.write(outputBuffer.data(), static_cast<std::streamsize>(outputBuffer.size()));
                        outputBuffer.clear();
                    }
                }
            } else if (!sawAnyHeaderLine && !warnedNoChromYet) {
                std::cerr << "Error: no #CHROM header found before variant lines.\n";
                warnedNoChromYet = true;
            }
        }

        p = lineEnd + 1;
    }

    // Flush remaining buffer
    if (!outputBuffer.empty()) {
        out.write(outputBuffer.data(), static_cast<std::streamsize>(outputBuffer.size()));
    }

    munmap(mapped, fileSize);
    close(fd);
    return 0;
}

// ============================================================================
// STDIN-BASED INDEXING (Legacy fallback)
// ============================================================================

void VCFXIndexer::createVCFIndex(std::istream &in, std::ostream &out) {
    bool foundChromHeader = false;
    bool warnedNoChromYet = false;
    bool sawAnyHeaderLine = false;

    static const size_t BUF_SIZE = 64 * 1024;
    char buffer[BUF_SIZE];

    // Accumulate partial lines here
    std::string leftover;
    leftover.reserve(BUF_SIZE);

    // Tracks total raw bytes read from the start
    std::int64_t totalOffset = 0;

    // Helper lambda to process each full line
    auto processLine = [&](const std::string &line, std::int64_t lineStartOffset) {
        if (line.empty()) {
            // Skip empty lines
            return;
        }

        // Trim leading space
        std::string trimmed = ltrim(line);

        // If it starts with '#', treat as header line
        if (!trimmed.empty() && trimmed[0] == '#') {
            sawAnyHeaderLine = true;

            // Possibly #CHROM
            if (!foundChromHeader) {
                auto fields = splitTabs(trimmed);
                if (fields.size() >= 2 && fields[0] == "#CHROM") {
                    foundChromHeader = true;
                    // Once #CHROM is found, print the header
                    out << "CHROM\tPOS\tFILE_OFFSET\n";
                }
            }
            return;
        }

        // It's a data line. If we haven't seen #CHROM, that's an error
        if (!foundChromHeader) {
            if (!sawAnyHeaderLine && !warnedNoChromYet) {
                std::cerr << "Error: no #CHROM header found before variant lines.\n";
                warnedNoChromYet = true;
            }
            return;
        }

        // We do have #CHROM; parse the fields
        auto fields = splitTabs(line);
        if (fields.size() < 2) {
            return; // Not enough fields to parse
        }

        const std::string &chrom = fields[0];
        const std::string &posStr = fields[1];

        std::int64_t posVal = 0;
        try {
            posVal = std::stoll(posStr);
        } catch (...) {
            // Not a valid integer => skip
            return;
        }

        // Valid line => output the index info
        out << chrom << "\t" << posVal << "\t" << lineStartOffset << "\n";
    };

    // Main read loop
    while (true) {
        in.read(buffer, BUF_SIZE);
        std::streamsize got = in.gcount();
        if (got <= 0) {
            break;
        }

        for (std::streamsize i = 0; i < got; ++i) {
            char c = buffer[i];
            leftover.push_back(c);
            totalOffset++;

            if (c == '\n') {
                // A full line has ended
                leftover.pop_back(); // remove the '\n'
                std::int64_t lineStartOffset = totalOffset - leftover.size() - 1;

                // Handle CRLF by removing trailing '\r'
                if (!leftover.empty() && leftover.back() == '\r') {
                    leftover.pop_back();
                }

                // Process the line
                processLine(leftover, lineStartOffset);
                leftover.clear();
            }
        }

        if (!in.good()) {
            // If stream ended or had an error
            break;
        }
    }

    // Handle any leftover partial line
    if (!leftover.empty()) {
        std::int64_t lineStartOffset = totalOffset - static_cast<std::int64_t>(leftover.size());
        if (!leftover.empty() && leftover.back() == '\r') {
            leftover.pop_back();
        }
        processLine(leftover, lineStartOffset);
        leftover.clear();
    }
}

// Optional main if you build as a single executable
static void show_help() {
    VCFXIndexer obj;
    char arg0[] = "VCFX_indexer";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_indexer", show_help))
        return 0;
    VCFXIndexer idx;
    return idx.run(argc, argv);
}
