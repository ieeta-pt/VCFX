/**
 * VCFX_indel_normalizer - High-performance INDEL normalization tool
 *
 * Normalizes INDEL variants by:
 * 1. Splitting multi-allelic lines into separate lines
 * 2. Removing common leading/trailing bases
 * 3. Adjusting POS accordingly
 *
 * Optimizations:
 * - Memory-mapped I/O with MADV_SEQUENTIAL | MADV_WILLNEED
 * - SIMD-accelerated newline/tab scanning (AVX2/SSE2/NEON)
 * - Zero-copy parsing with string_view
 * - Buffered output with direct write() syscalls
 */

#include "VCFX_indel_normalizer.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <string_view>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <getopt.h>

// ============================================================================
// SIMD Detection and Intrinsics
// ============================================================================
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
    #if defined(__AVX2__)
        #include <immintrin.h>
        #define USE_AVX2 1
    #elif defined(__SSE2__)
        #include <emmintrin.h>
        #define USE_SSE2 1
    #endif
#elif defined(__aarch64__) || defined(__ARM_NEON)
    #include <arm_neon.h>
    #define USE_NEON 1
#endif

// ============================================================================
// Memory-mapped file structure
// ============================================================================
struct MappedFile {
    const char *data = nullptr;
    size_t size = 0;
    int fd = -1;

    bool open(const char *path) {
        fd = ::open(path, O_RDONLY);
        if (fd < 0) return false;
        struct stat st;
        if (fstat(fd, &st) < 0) { ::close(fd); fd = -1; return false; }
        size = static_cast<size_t>(st.st_size);
        if (size == 0) { ::close(fd); fd = -1; return true; }
        data = static_cast<const char*>(mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0));
        if (data == MAP_FAILED) { data = nullptr; ::close(fd); fd = -1; return false; }
        madvise(const_cast<char*>(data), size, MADV_SEQUENTIAL | MADV_WILLNEED);
        return true;
    }

    void close() {
        if (data && size > 0) munmap(const_cast<char*>(data), size);
        if (fd >= 0) ::close(fd);
        data = nullptr; size = 0; fd = -1;
    }
};

// ============================================================================
// Output buffer for efficient writes
// ============================================================================
class OutputBuffer {
    static constexpr size_t CAPACITY = 4 * 1024 * 1024; // 4MB buffer
    char *buf;
    size_t pos = 0;
    int fd;
public:
    explicit OutputBuffer(int fd_) : fd(fd_) { buf = new char[CAPACITY]; }
    ~OutputBuffer() { flush(); delete[] buf; }

    void flush() {
        if (pos > 0) { ::write(fd, buf, pos); pos = 0; }
    }

    void ensureSpace(size_t n) {
        if (pos + n > CAPACITY) flush();
    }

    void append(char c) {
        if (pos >= CAPACITY) flush();
        buf[pos++] = c;
    }

    void append(const char *s, size_t len) {
        while (len > 0) {
            size_t chunk = std::min(len, CAPACITY - pos);
            memcpy(buf + pos, s, chunk);
            pos += chunk; s += chunk; len -= chunk;
            if (pos >= CAPACITY) flush();
        }
    }

    void append(std::string_view sv) { append(sv.data(), sv.size()); }
    void append(const std::string &s) { append(s.data(), s.size()); }

    // Write integer efficiently
    void writeInt(int val) {
        ensureSpace(16);
        if (val == 0) { buf[pos++] = '0'; return; }
        if (val < 0) { buf[pos++] = '-'; val = -val; }
        char tmp[16]; int i = 0;
        while (val > 0) { tmp[i++] = '0' + (val % 10); val /= 10; }
        while (i > 0) buf[pos++] = tmp[--i];
    }
};

// ============================================================================
// SIMD-accelerated scanning functions
// ============================================================================
#if defined(USE_AVX2)
static inline const char* findNewlineSIMD(const char* p, const char* end) {
    const __m256i nl = _mm256_set1_epi8('\n');
    while (p + 32 <= end) {
        __m256i chunk = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p));
        __m256i cmp = _mm256_cmpeq_epi8(chunk, nl);
        int mask = _mm256_movemask_epi8(cmp);
        if (mask) return p + __builtin_ctz(mask);
        p += 32;
    }
    while (p < end && *p != '\n') ++p;
    return p;
}

static inline const char* findTabSIMD(const char* p, const char* end) {
    const __m256i tab = _mm256_set1_epi8('\t');
    while (p + 32 <= end) {
        __m256i chunk = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p));
        __m256i cmp = _mm256_cmpeq_epi8(chunk, tab);
        int mask = _mm256_movemask_epi8(cmp);
        if (mask) return p + __builtin_ctz(mask);
        p += 32;
    }
    while (p < end && *p != '\t' && *p != '\n') ++p;
    return p;
}
#elif defined(USE_SSE2)
static inline const char* findNewlineSIMD(const char* p, const char* end) {
    const __m128i nl = _mm_set1_epi8('\n');
    while (p + 16 <= end) {
        __m128i chunk = _mm_loadu_si128(reinterpret_cast<const __m128i*>(p));
        __m128i cmp = _mm_cmpeq_epi8(chunk, nl);
        int mask = _mm_movemask_epi8(cmp);
        if (mask) return p + __builtin_ctz(mask);
        p += 16;
    }
    while (p < end && *p != '\n') ++p;
    return p;
}

static inline const char* findTabSIMD(const char* p, const char* end) {
    const __m128i tab = _mm_set1_epi8('\t');
    while (p + 16 <= end) {
        __m128i chunk = _mm_loadu_si128(reinterpret_cast<const __m128i*>(p));
        __m128i cmp = _mm_cmpeq_epi8(chunk, tab);
        int mask = _mm_movemask_epi8(cmp);
        if (mask) return p + __builtin_ctz(mask);
        p += 16;
    }
    while (p < end && *p != '\t' && *p != '\n') ++p;
    return p;
}
#elif defined(USE_NEON)
static inline const char* findNewlineSIMD(const char* p, const char* end) {
    const uint8x16_t nl = vdupq_n_u8('\n');
    while (p + 16 <= end) {
        uint8x16_t chunk = vld1q_u8(reinterpret_cast<const uint8_t*>(p));
        uint8x16_t cmp = vceqq_u8(chunk, nl);
        uint64_t lo = vgetq_lane_u64(vreinterpretq_u64_u8(cmp), 0);
        uint64_t hi = vgetq_lane_u64(vreinterpretq_u64_u8(cmp), 1);
        if (lo) return p + (__builtin_ctzll(lo) >> 3);
        if (hi) return p + 8 + (__builtin_ctzll(hi) >> 3);
        p += 16;
    }
    while (p < end && *p != '\n') ++p;
    return p;
}

static inline const char* findTabSIMD(const char* p, const char* end) {
    const uint8x16_t tab = vdupq_n_u8('\t');
    while (p + 16 <= end) {
        uint8x16_t chunk = vld1q_u8(reinterpret_cast<const uint8_t*>(p));
        uint8x16_t cmp = vceqq_u8(chunk, tab);
        uint64_t lo = vgetq_lane_u64(vreinterpretq_u64_u8(cmp), 0);
        uint64_t hi = vgetq_lane_u64(vreinterpretq_u64_u8(cmp), 1);
        if (lo) return p + (__builtin_ctzll(lo) >> 3);
        if (hi) return p + 8 + (__builtin_ctzll(hi) >> 3);
        p += 16;
    }
    while (p < end && *p != '\t' && *p != '\n') ++p;
    return p;
}
#else
// Fallback using memchr
static inline const char* findNewlineSIMD(const char* p, const char* end) {
    const char* found = static_cast<const char*>(memchr(p, '\n', end - p));
    return found ? found : end;
}

static inline const char* findTabSIMD(const char* p, const char* end) {
    while (p < end && *p != '\t' && *p != '\n') ++p;
    return p;
}
#endif

// ============================================================================
// Zero-copy field extraction
// ============================================================================
static inline std::string_view getField(const char* lineStart, const char* lineEnd, int fieldIdx) {
    const char* p = lineStart;
    for (int i = 0; i < fieldIdx && p < lineEnd; ++i) {
        p = findTabSIMD(p, lineEnd);
        if (p < lineEnd) ++p;
    }
    if (p >= lineEnd) return {};
    const char* fieldEnd = findTabSIMD(p, lineEnd);
    return std::string_view(p, fieldEnd - p);
}

// ============================================================================
// Parse integer from string_view
// ============================================================================
static inline int parseInt(std::string_view sv) {
    int result = 0;
    for (char c : sv) {
        if (c < '0' || c > '9') break;
        result = result * 10 + (c - '0');
    }
    return result;
}

// ============================================================================
// Normalize a single variant (minimal left alignment)
// Returns false if after trimming, REF or ALT is empty or they are the same
// ============================================================================
static bool doNormalizeVariant(int &posInt, std::string &ref, std::string &alt) {
    if (ref == alt) return false;

    // 1) Remove leading common bases (keep at least 1)
    size_t prefixCount = 0;
    size_t minLen = std::min(ref.size(), alt.size());
    while (prefixCount < minLen && ref[prefixCount] == alt[prefixCount]) {
        prefixCount++;
    }

    if (prefixCount > 0) {
        size_t removeLeading = prefixCount - 1;
        if (removeLeading > 0) {
            ref.erase(0, removeLeading);
            alt.erase(0, removeLeading);
            posInt += static_cast<int>(removeLeading);
        }
    }

    // 2) Remove trailing common bases (keep at least 1)
    size_t rLen = ref.size();
    size_t aLen = alt.size();
    size_t suffixCount = 0;
    while (suffixCount < rLen && suffixCount < aLen) {
        if (ref[rLen - 1 - suffixCount] == alt[aLen - 1 - suffixCount]) {
            suffixCount++;
        } else {
            break;
        }
    }

    if (suffixCount > 0) {
        size_t removeS = suffixCount - 1;
        if (removeS > 0) {
            ref.erase(rLen - removeS, removeS);
            alt.erase(aLen - removeS, removeS);
        }
    }

    // Check if after trimming we have empty or identical
    if (ref.empty() || alt.empty()) return false;
    if (ref == alt) return false;

    return true;
}

// ============================================================================
// Split ALT field by commas
// ============================================================================
static void splitAlt(std::string_view alt, std::vector<std::string_view> &altList) {
    altList.clear();
    const char* p = alt.data();
    const char* end = p + alt.size();

    while (p < end) {
        const char* start = p;
        while (p < end && *p != ',') ++p;
        altList.push_back(std::string_view(start, p - start));
        if (p < end) ++p; // skip comma
    }
}

// ============================================================================
// Process mmap'd file
// ============================================================================
static void processMmap(const MappedFile &mf, OutputBuffer &out, bool quiet) {
    const char* p = mf.data;
    const char* end = mf.data + mf.size;

    bool foundChromHeader = false;
    std::vector<std::string_view> altList;
    altList.reserve(8);

    size_t lineCount = 0;
    size_t variantCount = 0;
    size_t outputCount = 0;

    while (p < end) {
        const char* lineStart = p;
        const char* lineEnd = findNewlineSIMD(p, end);

        // Handle \r\n
        const char* adjLineEnd = lineEnd;
        if (adjLineEnd > lineStart && *(adjLineEnd - 1) == '\r') --adjLineEnd;

        size_t lineLen = adjLineEnd - lineStart;

        // Empty line
        if (lineLen == 0) {
            out.append('\n');
            p = lineEnd + 1;
            continue;
        }

        // Handle comments/header
        if (*lineStart == '#') {
            out.append(lineStart, lineLen);
            out.append('\n');
            if (lineLen >= 6 &&
                lineStart[0] == '#' && lineStart[1] == 'C' && lineStart[2] == 'H' &&
                lineStart[3] == 'R' && lineStart[4] == 'O' && lineStart[5] == 'M') {
                foundChromHeader = true;
            }
            p = lineEnd + 1;
            continue;
        }

        if (!foundChromHeader) {
            if (!quiet) std::cerr << "Error: encountered data line before #CHROM header.\n";
            return;
        }

        lineCount++;

        // Extract fields
        std::string_view chrom = getField(lineStart, adjLineEnd, 0);
        std::string_view posStr = getField(lineStart, adjLineEnd, 1);
        std::string_view id = getField(lineStart, adjLineEnd, 2);
        std::string_view refSv = getField(lineStart, adjLineEnd, 3);
        std::string_view altSv = getField(lineStart, adjLineEnd, 4);

        // Check if we have enough fields
        if (refSv.empty() || altSv.empty()) {
            // Not enough columns, pass through
            out.append(lineStart, lineLen);
            out.append('\n');
            p = lineEnd + 1;
            continue;
        }

        int posInt = parseInt(posStr);

        // Get post-ALT columns (QUAL, FILTER, INFO, FORMAT, samples...)
        // Find position after ALT field
        const char* postAlt = lineStart;
        for (int i = 0; i < 5 && postAlt < adjLineEnd; ++i) {
            postAlt = findTabSIMD(postAlt, adjLineEnd);
            if (postAlt < adjLineEnd) ++postAlt;
        }
        // Back up to include the tab before QUAL
        std::string_view postCols;
        if (postAlt > lineStart && postAlt <= adjLineEnd) {
            // postAlt now points to QUAL field start
            // We want "\tQUAL\tFILTER\t..." so we need to go back one char
            postCols = std::string_view(postAlt - 1, adjLineEnd - (postAlt - 1));
        }

        // Split ALT if multi-allelic
        splitAlt(altSv, altList);
        variantCount++;

        std::string ref(refSv);

        if (altList.size() == 1) {
            // Single alt
            std::string alt(altList[0]);
            int newPos = posInt;
            std::string newRef = ref;
            std::string newAlt = alt;

            bool ok = doNormalizeVariant(newPos, newRef, newAlt);
            if (!ok) {
                // Output original line unchanged
                out.append(lineStart, lineLen);
                out.append('\n');
            } else {
                // Output normalized
                out.append(chrom);
                out.append('\t');
                out.writeInt(newPos);
                out.append('\t');
                out.append(id);
                out.append('\t');
                out.append(newRef);
                out.append('\t');
                out.append(newAlt);
                out.append(postCols);
                out.append('\n');
            }
            outputCount++;
        } else {
            // Multiple alts - produce multiple lines
            for (size_t i = 0; i < altList.size(); ++i) {
                std::string alt(altList[i]);
                int newPos = posInt;
                std::string newRef = ref;
                std::string newAlt = alt;

                bool ok = doNormalizeVariant(newPos, newRef, newAlt);
                if (!ok) {
                    // Output with original values for this alt
                    out.append(chrom);
                    out.append('\t');
                    out.writeInt(posInt);
                    out.append('\t');
                    out.append(id);
                    out.append('\t');
                    out.append(ref);
                    out.append('\t');
                    out.append(altList[i]);
                    out.append(postCols);
                    out.append('\n');
                } else {
                    out.append(chrom);
                    out.append('\t');
                    out.writeInt(newPos);
                    out.append('\t');
                    out.append(id);
                    out.append('\t');
                    out.append(newRef);
                    out.append('\t');
                    out.append(newAlt);
                    out.append(postCols);
                    out.append('\n');
                }
                outputCount++;
            }
        }

        p = lineEnd + 1;
    }

    if (!quiet) {
        std::cerr << "Processed " << variantCount << " variants, output " << outputCount << " lines\n";
    }
}

// ============================================================================
// Process stdin (fallback)
// ============================================================================
void VCFXIndelNormalizer::normalizeIndels(std::istream &in, std::ostream &out) {
    std::string line;
    bool foundChromHeader = false;
    std::vector<std::string> fields;
    fields.reserve(16);

    while (std::getline(in, line)) {
        if (line.empty()) {
            out << "\n";
            continue;
        }
        if (line[0] == '#') {
            out << line << "\n";
            if (line.rfind("#CHROM", 0) == 0) {
                foundChromHeader = true;
            }
            continue;
        }
        if (!foundChromHeader) {
            std::cerr << "Error: encountered data line before #CHROM header.\n";
            return;
        }

        // Parse the line
        fields.clear();
        size_t start = 0;
        while (start < line.size()) {
            size_t tabPos = line.find('\t', start);
            if (tabPos == std::string::npos) {
                fields.push_back(line.substr(start));
                break;
            }
            fields.push_back(line.substr(start, tabPos - start));
            start = tabPos + 1;
        }

        if (fields.size() < 10) {
            out << line << "\n";
            continue;
        }

        std::string &chrom = fields[0];
        std::string &posStr = fields[1];
        std::string &id = fields[2];
        std::string &ref = fields[3];
        std::string &alt = fields[4];

        int posInt = 0;
        try {
            posInt = std::stoi(posStr);
        } catch (...) {
            out << line << "\n";
            continue;
        }

        std::string postCols;
        {
            std::ostringstream oss;
            for (size_t c = 5; c < fields.size(); c++) {
                oss << "\t" << fields[c];
            }
            postCols = oss.str();
        }

        // Split multi-allelic
        std::vector<std::string> altList;
        {
            std::stringstream altSS(alt);
            std::string a;
            while (std::getline(altSS, a, ',')) {
                altList.push_back(a);
            }
        }

        if (altList.size() == 1) {
            std::string altOne = altList[0];
            std::string newRef = ref;
            std::string newAlt = altOne;
            int newPos = posInt;
            bool ok = doNormalizeVariant(newPos, newRef, newAlt);
            if (!ok) {
                out << line << "\n";
            } else {
                out << chrom << "\t" << newPos << "\t" << id << "\t" << newRef << "\t" << newAlt << postCols << "\n";
            }
        } else {
            for (size_t i = 0; i < altList.size(); i++) {
                std::string altOne = altList[i];
                std::string newRef = ref;
                std::string newAlt = altOne;
                int newPos = posInt;
                bool ok = doNormalizeVariant(newPos, newRef, newAlt);
                if (!ok) {
                    out << chrom << "\t" << posInt << "\t" << id << "\t" << ref << "\t" << altOne << postCols << "\n";
                } else {
                    out << chrom << "\t" << newPos << "\t" << id << "\t" << newRef << "\t" << newAlt << postCols << "\n";
                }
            }
        }
    }
}

// ============================================================================
// Help message
// ============================================================================
void VCFXIndelNormalizer::displayHelp() {
    std::cout << "VCFX_indel_normalizer v1.1 - High-performance INDEL normalization\n\n"
              << "Usage:\n"
              << "  VCFX_indel_normalizer [OPTIONS] [input.vcf]\n"
              << "  VCFX_indel_normalizer [OPTIONS] < input.vcf > output.vcf\n\n"
              << "Options:\n"
              << "  -i, --input FILE   Input VCF file (uses memory-mapping for best performance)\n"
              << "  -q, --quiet        Suppress informational messages\n"
              << "  -h, --help         Display this help message and exit\n"
              << "  -v, --version      Show program version and exit\n\n"
              << "Description:\n"
              << "  Normalizes INDEL variants by:\n"
              << "   1) Splitting multi-ALT lines into separate lines.\n"
              << "   2) Removing the longest shared prefix from REF/ALT, adjusting POS.\n"
              << "   3) Removing the largest shared suffix from REF/ALT.\n\n"
              << "  Note: true left alignment for repeated motifs requires the full reference genome.\n\n"
              << "Performance:\n"
              << "  - Memory-mapped I/O: Use -i flag for ~15-20x faster processing\n"
              << "  - SIMD acceleration for line/field scanning\n\n"
              << "Example:\n"
              << "  VCFX_indel_normalizer -i input.vcf > normalized.vcf\n"
              << "  VCFX_indel_normalizer < input.vcf > normalized.vcf\n";
}

// ============================================================================
// Run method
// ============================================================================
int VCFXIndelNormalizer::run(int argc, char *argv[]) {
    const char *inputFile = nullptr;
    bool quiet = false;

    static struct option longOpts[] = {
        {"input",   required_argument, nullptr, 'i'},
        {"quiet",   no_argument,       nullptr, 'q'},
        {"help",    no_argument,       nullptr, 'h'},
        {"version", no_argument,       nullptr, 'v'},
        {nullptr, 0, nullptr, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "i:qhv", longOpts, nullptr)) != -1) {
        switch (opt) {
            case 'i': inputFile = optarg; break;
            case 'q': quiet = true; break;
            case 'h': displayHelp(); return 0;
            case 'v': std::cout << "VCFX_indel_normalizer v1.1\n"; return 0;
            default: displayHelp(); return 1;
        }
    }

    // Check for positional argument
    if (!inputFile && optind < argc) {
        inputFile = argv[optind];
    }

    if (inputFile) {
        // Memory-mapped mode
        MappedFile mf;
        if (!mf.open(inputFile)) {
            std::cerr << "Error: Cannot open file: " << inputFile << "\n";
            return 1;
        }

        if (!quiet) {
            std::cerr << "Processing " << inputFile << " (" << (mf.size / (1024*1024)) << " MB)\n";
        }

        OutputBuffer out(STDOUT_FILENO);
        processMmap(mf, out, quiet);
        mf.close();
    } else {
        // Stdin mode
        normalizeIndels(std::cin, std::cout);
    }

    return 0;
}

// ============================================================================
// Main
// ============================================================================
static void show_help() {
    VCFXIndelNormalizer obj;
    char arg0[] = "VCFX_indel_normalizer";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);

    VCFXIndelNormalizer norm;
    return norm.run(argc, argv);
}
