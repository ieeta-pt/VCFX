/**
 * VCFX_allele_freq_calc - High-performance allele frequency calculator
 *
 * Optimizations:
 * - Memory-mapped I/O with MADV_SEQUENTIAL | MADV_WILLNEED
 * - SIMD-accelerated newline/tab scanning (AVX2/SSE2/NEON)
 * - Zero-copy parsing with string_view
 * - Buffered output with direct write() syscalls
 * - FORMAT field caching (avoid re-parsing identical FORMAT strings)
 */

#include "VCFX_allele_freq_calc.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <string_view>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
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

    // Write integer efficiently
    void writeInt(int val) {
        ensureSpace(16);
        if (val == 0) { buf[pos++] = '0'; return; }
        if (val < 0) { buf[pos++] = '-'; val = -val; }
        char tmp[16]; int i = 0;
        while (val > 0) { tmp[i++] = '0' + (val % 10); val /= 10; }
        while (i > 0) buf[pos++] = tmp[--i];
    }

    // Write double with 4 decimal places (like std::fixed << std::setprecision(4))
    void writeDouble4(double val) {
        ensureSpace(24);
        if (val < 0) { buf[pos++] = '-'; val = -val; }

        // Scale and round to 4 decimal places
        unsigned long long scaled = static_cast<unsigned long long>(val * 10000.0 + 0.5);
        unsigned long long intPart = scaled / 10000;
        unsigned long long fracPart = scaled % 10000;

        // Write integer part
        if (intPart == 0) {
            buf[pos++] = '0';
        } else {
            char tmp[20]; int i = 0;
            while (intPart > 0) { tmp[i++] = '0' + (intPart % 10); intPart /= 10; }
            while (i > 0) buf[pos++] = tmp[--i];
        }

        // Write decimal point and fractional part (always 4 digits)
        buf[pos++] = '.';
        buf[pos++] = '0' + (fracPart / 1000) % 10;
        buf[pos++] = '0' + (fracPart / 100) % 10;
        buf[pos++] = '0' + (fracPart / 10) % 10;
        buf[pos++] = '0' + fracPart % 10;
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
// Parse genotype and count alleles
// Returns (altCount, totalCount)
// ============================================================================
static inline void parseGenotypeAndCount(std::string_view gt, int &altCount, int &totalCount) {
    const char* p = gt.data();
    const char* end = p + gt.size();

    while (p < end) {
        // Skip separators
        while (p < end && (*p == '/' || *p == '|')) ++p;
        if (p >= end) break;

        // Parse allele
        const char* alleleStart = p;
        while (p < end && *p != '/' && *p != '|') ++p;

        if (alleleStart == p) continue;

        // Check for missing
        if (*alleleStart == '.') continue;

        // Check if numeric and whether it's zero
        bool isZero = true;
        bool numeric = true;
        for (const char* c = alleleStart; c < p; ++c) {
            if (*c < '0' || *c > '9') { numeric = false; break; }
            if (*c != '0') isZero = false;
        }

        if (!numeric) continue;

        totalCount++;
        if (!isZero) altCount++;
    }
}

// ============================================================================
// Find GT index in FORMAT field
// ============================================================================
static inline int findGTIndex(std::string_view format) {
    const char* p = format.data();
    const char* end = p + format.size();
    int idx = 0;

    while (p < end) {
        const char* fieldStart = p;
        while (p < end && *p != ':') ++p;

        size_t len = p - fieldStart;
        if (len == 2 && fieldStart[0] == 'G' && fieldStart[1] == 'T') {
            return idx;
        }

        idx++;
        if (p < end) ++p; // skip ':'
    }
    return -1;
}

// ============================================================================
// Extract GT from sample at given index
// ============================================================================
static inline std::string_view extractGT(std::string_view sample, int gtIndex) {
    const char* p = sample.data();
    const char* end = p + sample.size();

    // Skip to GT field
    for (int i = 0; i < gtIndex && p < end; ++i) {
        while (p < end && *p != ':') ++p;
        if (p < end) ++p;
    }

    if (p >= end) return {};

    const char* gtStart = p;
    while (p < end && *p != ':') ++p;

    return std::string_view(gtStart, p - gtStart);
}

// ============================================================================
// Process mmap'd file
// ============================================================================
static void processMmap(const MappedFile &mf, OutputBuffer &out, bool quiet) {
    const char* p = mf.data;
    const char* end = mf.data + mf.size;

    bool foundChromHeader = false;

    // FORMAT field caching
    std::string cachedFormat;
    int cachedGtIndex = -1;

    // Write header
    out.append("CHROM\tPOS\tID\tREF\tALT\tAllele_Frequency\n");

    size_t lineCount = 0;
    size_t variantCount = 0;

    while (p < end) {
        const char* lineStart = p;
        const char* lineEnd = findNewlineSIMD(p, end);

        // Handle \r\n
        const char* adjLineEnd = lineEnd;
        if (adjLineEnd > lineStart && *(adjLineEnd - 1) == '\r') --adjLineEnd;

        if (lineStart == adjLineEnd) {
            p = lineEnd + 1;
            continue;
        }

        // Handle comments/header
        if (*lineStart == '#') {
            if (adjLineEnd - lineStart >= 6 &&
                lineStart[0] == '#' && lineStart[1] == 'C' && lineStart[2] == 'H' &&
                lineStart[3] == 'R' && lineStart[4] == 'O' && lineStart[5] == 'M') {
                foundChromHeader = true;
            }
            p = lineEnd + 1;
            continue;
        }

        if (!foundChromHeader) {
            if (!quiet) std::cerr << "Warning: Data line encountered before #CHROM header. Skipping.\n";
            p = lineEnd + 1;
            continue;
        }

        lineCount++;

        // Extract fields 0-8 (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)
        std::string_view chrom = getField(lineStart, adjLineEnd, 0);
        std::string_view pos = getField(lineStart, adjLineEnd, 1);
        std::string_view id = getField(lineStart, adjLineEnd, 2);
        std::string_view ref = getField(lineStart, adjLineEnd, 3);
        std::string_view alt = getField(lineStart, adjLineEnd, 4);
        std::string_view format = getField(lineStart, adjLineEnd, 8);

        if (format.empty()) {
            p = lineEnd + 1;
            continue;
        }

        // Get GT index (with caching)
        int gtIndex;
        if (format == cachedFormat) {
            gtIndex = cachedGtIndex;
        } else {
            gtIndex = findGTIndex(format);
            cachedFormat = std::string(format);
            cachedGtIndex = gtIndex;
        }

        if (gtIndex < 0) {
            p = lineEnd + 1;
            continue;
        }

        // Count alleles across all samples
        int altCount = 0;
        int totalCount = 0;

        // Find start of sample columns (field 9)
        const char* sampleStart = lineStart;
        int fieldIdx = 0;
        while (fieldIdx < 9 && sampleStart < adjLineEnd) {
            sampleStart = findTabSIMD(sampleStart, adjLineEnd);
            if (sampleStart < adjLineEnd) {
                ++sampleStart;
                ++fieldIdx;
            }
        }

        // Process each sample
        while (sampleStart < adjLineEnd) {
            const char* sampleEnd = findTabSIMD(sampleStart, adjLineEnd);
            std::string_view sample(sampleStart, sampleEnd - sampleStart);

            std::string_view gt = extractGT(sample, gtIndex);
            if (!gt.empty()) {
                parseGenotypeAndCount(gt, altCount, totalCount);
            }

            if (sampleEnd >= adjLineEnd) break;
            sampleStart = sampleEnd + 1;
        }

        // Calculate frequency
        double freq = (totalCount > 0) ?
            static_cast<double>(altCount) / static_cast<double>(totalCount) : 0.0;

        // Write output
        out.append(chrom);
        out.append('\t');
        out.append(pos);
        out.append('\t');
        out.append(id);
        out.append('\t');
        out.append(ref);
        out.append('\t');
        out.append(alt);
        out.append('\t');
        out.writeDouble4(freq);
        out.append('\n');

        variantCount++;
        p = lineEnd + 1;
    }

    if (!quiet) {
        std::cerr << "Processed " << variantCount << " variants from " << lineCount << " data lines\n";
    }
}

// ============================================================================
// Process stdin (fallback)
// ============================================================================
static void processStdin(std::istream &in, std::ostream &outStream, bool quiet) {
    std::string line;
    line.reserve(65536);

    bool foundChromHeader = false;

    // FORMAT field caching
    std::string cachedFormat;
    int cachedGtIndex = -1;

    outStream << "CHROM\tPOS\tID\tREF\tALT\tAllele_Frequency\n";

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            if (line.size() >= 6 && line.compare(0, 6, "#CHROM") == 0) {
                foundChromHeader = true;
            }
            continue;
        }

        if (!foundChromHeader) {
            if (!quiet) std::cerr << "Warning: Data line encountered before #CHROM header. Skipping.\n";
            continue;
        }

        // Parse fields
        std::string_view lineView(line);
        std::vector<std::string_view> fields;
        fields.reserve(16);

        size_t start = 0;
        while (start < lineView.size()) {
            size_t tabPos = lineView.find('\t', start);
            if (tabPos == std::string_view::npos) {
                fields.push_back(lineView.substr(start));
                break;
            }
            fields.push_back(lineView.substr(start, tabPos - start));
            start = tabPos + 1;
        }

        if (fields.size() < 9) {
            if (!quiet) std::cerr << "Warning: Skipping invalid VCF line (fewer than 9 fields).\n";
            continue;
        }

        std::string_view format = fields[8];

        // Get GT index
        int gtIndex;
        if (format == cachedFormat) {
            gtIndex = cachedGtIndex;
        } else {
            gtIndex = findGTIndex(format);
            cachedFormat = std::string(format);
            cachedGtIndex = gtIndex;
        }

        if (gtIndex < 0) continue;

        // Count alleles
        int altCount = 0;
        int totalCount = 0;

        for (size_t i = 9; i < fields.size(); ++i) {
            std::string_view gt = extractGT(fields[i], gtIndex);
            if (!gt.empty()) {
                parseGenotypeAndCount(gt, altCount, totalCount);
            }
        }

        double freq = (totalCount > 0) ?
            static_cast<double>(altCount) / static_cast<double>(totalCount) : 0.0;

        outStream << fields[0] << "\t" << fields[1] << "\t" << fields[2] << "\t"
                  << fields[3] << "\t" << fields[4] << "\t"
                  << std::fixed << std::setprecision(4) << freq << "\n";
    }
}

// ============================================================================
// Help message
// ============================================================================
void printHelp() {
    std::cout << "VCFX_allele_freq_calc v1.1 - High-performance allele frequency calculator\n\n"
              << "Usage:\n"
              << "  VCFX_allele_freq_calc [OPTIONS] [input.vcf]\n"
              << "  VCFX_allele_freq_calc [OPTIONS] < input.vcf > output.tsv\n\n"
              << "Options:\n"
              << "  -i, --input FILE   Input VCF file (uses memory-mapping for best performance)\n"
              << "  -q, --quiet        Suppress informational messages\n"
              << "  -h, --help         Display this help message and exit\n"
              << "  -v, --version      Show program version and exit\n\n"
              << "Description:\n"
              << "  Calculates allele frequency for each variant in a VCF file.\n"
              << "  Allele frequency is computed as (#ALT alleles) / (total #alleles),\n"
              << "  counting any non-zero numeric allele (1,2,3,...) as ALT.\n\n"
              << "Output Format:\n"
              << "  CHROM  POS  ID  REF  ALT  Allele_Frequency\n\n"
              << "Performance:\n"
              << "  - Memory-mapped I/O: Use -i flag for ~15-20x faster processing\n"
              << "  - SIMD acceleration for line/field scanning\n"
              << "  - Zero-copy parsing with string_view\n\n"
              << "Examples:\n"
              << "  VCFX_allele_freq_calc -i input.vcf > frequencies.tsv\n"
              << "  VCFX_allele_freq_calc < input.vcf > frequencies.tsv\n";
}

// ============================================================================
// Main
// ============================================================================
int main(int argc, char *argv[]) {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);

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
            case 'h': printHelp(); return 0;
            case 'v': std::cout << "VCFX_allele_freq_calc v1.1\n"; return 0;
            default: printHelp(); return 1;
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
        if (std::cin.peek() == EOF) {
            printHelp();
            return 1;
        }
        processStdin(std::cin, std::cout, quiet);
    }

    return 0;
}
