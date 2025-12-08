#include "VCFX_ld_calculator.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fcntl.h>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <thread>
#include <unistd.h>
#include <vector>

// ======================================================================
// SIMD detection and includes
// ======================================================================
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#define USE_NEON 1
#elif defined(__AVX2__)
#include <immintrin.h>
#define USE_AVX2 1
#elif defined(__SSE2__)
#include <emmintrin.h>
#define USE_SSE2 1
#endif

// ======================================================================
// Memory-mapped file helper
// ======================================================================
struct MappedFile {
    const char *data = nullptr;
    size_t size = 0;
    int fd = -1;

    bool open(const char *path) {
        fd = ::open(path, O_RDONLY);
        if (fd < 0) return false;
        struct stat st;
        if (fstat(fd, &st) < 0) { ::close(fd); fd = -1; return false; }
        size = st.st_size;
        if (size == 0) { ::close(fd); fd = -1; return false; }
        data = static_cast<const char *>(mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0));
        if (data == MAP_FAILED) { ::close(fd); fd = -1; data = nullptr; return false; }
        madvise((void *)data, size, MADV_SEQUENTIAL | MADV_WILLNEED);
        return true;
    }
    void close() {
        if (data && data != MAP_FAILED) munmap((void *)data, size);
        if (fd >= 0) ::close(fd);
        data = nullptr;
        fd = -1;
    }
    ~MappedFile() { close(); }
};

// ======================================================================
// Output buffer for batched I/O
// ======================================================================
class OutputBuffer {
    std::string buf;
    int fd;
    static constexpr size_t FLUSH_THRESHOLD = 4 * 1024 * 1024; // 4MB
public:
    explicit OutputBuffer(int outFd) : fd(outFd) { buf.reserve(FLUSH_THRESHOLD + 65536); }
    ~OutputBuffer() { flush(); }
    void append(const char *s, size_t len) {
        buf.append(s, len);
        if (buf.size() >= FLUSH_THRESHOLD) flush();
    }
    void append(const std::string &s) { append(s.data(), s.size()); }
    void appendChar(char c) { buf.push_back(c); if (buf.size() >= FLUSH_THRESHOLD) flush(); }
    void flush() {
        if (!buf.empty()) {
            ::write(fd, buf.data(), buf.size());
            buf.clear();
        }
    }
};

// ======================================================================
// SIMD-accelerated newline finder
// ======================================================================
#if USE_NEON
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
#elif USE_AVX2
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
#elif USE_SSE2
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
#else
static inline const char* findNewlineSIMD(const char* p, const char* end) {
    const char* found = static_cast<const char*>(memchr(p, '\n', end - p));
    return found ? found : end;
}
#endif

// ======================================================================
// Zero-allocation helper functions for high-performance parsing
// ======================================================================

// Fast genotype parsing from raw pointer - no string allocation
// Returns: 0=0/0, 1=0/1 or 1/0, 2=1/1, -1=missing/invalid
static inline int parseGenotypeRaw(const char* s, size_t len) {
    if (len == 0) return -1;
    if (len == 1 && s[0] == '.') return -1;
    if (len == 3 && s[0] == '.' && (s[1] == '/' || s[1] == '|') && s[2] == '.') return -1;

    size_t sepPos = 0;
    for (size_t i = 0; i < len; i++) {
        if (s[i] == '/' || s[i] == '|') { sepPos = i; break; }
    }
    if (sepPos == 0 || sepPos >= len - 1) return -1;

    int a1 = 0;
    for (size_t i = 0; i < sepPos; i++) {
        char c = s[i];
        if (c == '.') return -1;
        if (c < '0' || c > '9') return -1;
        a1 = a1 * 10 + (c - '0');
    }

    int a2 = 0;
    for (size_t i = sepPos + 1; i < len; i++) {
        char c = s[i];
        if (c == '.') return -1;
        if (c < '0' || c > '9') return -1;
        a2 = a2 * 10 + (c - '0');
    }

    if (a1 > 1 || a2 > 1) return -1;
    return a1 + a2;
}

// Extract GT field from sample field (stops at first colon)
static inline bool extractGT(const char* sample, size_t sampleLen,
                              const char*& gtStart, size_t& gtLen) {
    gtStart = sample;
    gtLen = sampleLen;
    for (size_t i = 0; i < sampleLen; i++) {
        if (sample[i] == ':') { gtLen = i; break; }
    }
    return gtLen > 0;
}

// Fast integer parsing from raw pointer
static inline bool fastParseInt(const char* s, size_t len, int& result) {
    if (len == 0) return false;
    result = 0;
    for (size_t i = 0; i < len; i++) {
        char c = s[i];
        if (c < '0' || c > '9') return false;
        result = result * 10 + (c - '0');
    }
    return true;
}

// Fast double to string with 4 decimal places
static inline size_t formatR2(double r2, char* buf) {
    if (r2 <= 0.0) { memcpy(buf, "0.0000", 6); return 6; }
    if (r2 >= 1.0) { memcpy(buf, "1.0000", 6); return 6; }
    buf[0] = '0'; buf[1] = '.';
    int val = static_cast<int>(r2 * 10000.0 + 0.5);
    if (val > 9999) val = 9999;
    buf[5] = '0' + (val % 10); val /= 10;
    buf[4] = '0' + (val % 10); val /= 10;
    buf[3] = '0' + (val % 10); val /= 10;
    buf[2] = '0' + (val % 10);
    return 6;
}

// Fast integer to string
static inline size_t formatInt(int val, char* buf) {
    if (val == 0) { buf[0] = '0'; return 1; }
    char temp[12];
    int pos = 0;
    bool neg = val < 0;
    if (neg) val = -val;
    while (val > 0) { temp[pos++] = '0' + (val % 10); val /= 10; }
    size_t len = 0;
    if (neg) buf[len++] = '-';
    while (pos > 0) buf[len++] = temp[--pos];
    return len;
}

// ======================================================================
// Optimized LDVariant with pre-computed statistics for fast r² calculation
// ======================================================================
struct LDVariantOpt {
    std::string chrom;
    int pos;
    std::string id;
    std::vector<int8_t> genotype;

    // Pre-computed statistics for r² calculation
    int validCount = 0;
    long sumX = 0;
    long sumX2 = 0;
    double meanX = 0;
    double varX = 0;

    void computeStats() {
        validCount = 0; sumX = 0; sumX2 = 0;
        for (int8_t g : genotype) {
            if (g >= 0) {
                validCount++;
                sumX += g;
                sumX2 += g * g;
            }
        }
        if (validCount > 0) {
            meanX = static_cast<double>(sumX) / validCount;
            varX = static_cast<double>(sumX2) / validCount - meanX * meanX;
        } else {
            meanX = 0; varX = 0;
        }
    }
};

// ======================================================================
// SIMD-accelerated r² computation
// ======================================================================
#if USE_NEON
static inline double computeRsqSIMD(const int8_t* g1, const int8_t* g2, size_t sz) {
    int32x4_t vn = vdupq_n_s32(0);
    int32x4_t vsumX = vdupq_n_s32(0);
    int32x4_t vsumY = vdupq_n_s32(0);
    int32x4_t vsumXY = vdupq_n_s32(0);
    int32x4_t vsumX2 = vdupq_n_s32(0);
    int32x4_t vsumY2 = vdupq_n_s32(0);

    const int8x16_t negOne = vdupq_n_s8(-1);
    size_t i = 0;

    for (; i + 16 <= sz; i += 16) {
        int8x16_t x = vld1q_s8(g1 + i);
        int8x16_t y = vld1q_s8(g2 + i);

        // Valid mask: both >= 0
        uint8x16_t validX = vcgtq_s8(x, negOne);
        uint8x16_t validY = vcgtq_s8(y, negOne);
        uint8x16_t valid = vandq_u8(validX, validY);

        // Mask out invalid values
        x = vandq_s8(x, vreinterpretq_s8_u8(valid));
        y = vandq_s8(y, vreinterpretq_s8_u8(valid));

        // Process in 4 groups of 4 bytes each
        for (int j = 0; j < 4; j++) {
            int8x8_t xl = (j < 2) ? vget_low_s8(x) : vget_high_s8(x);
            int8x8_t yl = (j < 2) ? vget_low_s8(y) : vget_high_s8(y);
            uint8x8_t vl = (j < 2) ? vget_low_u8(valid) : vget_high_u8(valid);

            if (j == 1 || j == 3) {
                xl = vext_s8(xl, xl, 4);
                yl = vext_s8(yl, yl, 4);
                vl = vext_u8(vl, vl, 4);
            }

            int16x4_t x16 = vget_low_s16(vmovl_s8(xl));
            int16x4_t y16 = vget_low_s16(vmovl_s8(yl));
            uint16x4_t v16 = vget_low_u16(vmovl_u8(vl));

            int32x4_t x32 = vmovl_s16(x16);
            int32x4_t y32 = vmovl_s16(y16);
            int32x4_t v32 = vreinterpretq_s32_u32(vmovl_u16(v16));

            vn = vaddq_s32(vn, vandq_s32(v32, vdupq_n_s32(1)));
            vsumX = vaddq_s32(vsumX, x32);
            vsumY = vaddq_s32(vsumY, y32);
            vsumXY = vaddq_s32(vsumXY, vmulq_s32(x32, y32));
            vsumX2 = vaddq_s32(vsumX2, vmulq_s32(x32, x32));
            vsumY2 = vaddq_s32(vsumY2, vmulq_s32(y32, y32));
        }
    }

    // Horizontal sum
    int n = vaddvq_s32(vn);
    long sumX = vaddvq_s32(vsumX);
    long sumY = vaddvq_s32(vsumY);
    long sumXY = vaddvq_s32(vsumXY);
    long sumX2 = vaddvq_s32(vsumX2);
    long sumY2 = vaddvq_s32(vsumY2);

    // Scalar tail
    for (; i < sz; i++) {
        int8_t x = g1[i];
        int8_t y = g2[i];
        if (x >= 0 && y >= 0) {
            n++;
            sumX += x;
            sumY += y;
            sumXY += x * y;
            sumX2 += x * x;
            sumY2 += y * y;
        }
    }

    if (n < 2) return 0.0;
    double meanX = static_cast<double>(sumX) / n;
    double meanY = static_cast<double>(sumY) / n;
    double cov = static_cast<double>(sumXY) / n - meanX * meanY;
    double varX = static_cast<double>(sumX2) / n - meanX * meanX;
    double varY = static_cast<double>(sumY2) / n - meanY * meanY;
    if (varX <= 0.0 || varY <= 0.0) return 0.0;
    double r = cov / (std::sqrt(varX) * std::sqrt(varY));
    return r * r;
}
#else
// Scalar fallback with loop unrolling
static inline double computeRsqSIMD(const int8_t* g1, const int8_t* g2, size_t sz) {
    int n = 0;
    long sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, sumY2 = 0;

    size_t i = 0;
    for (; i + 4 <= sz; i += 4) {
        for (int j = 0; j < 4; j++) {
            int8_t x = g1[i + j];
            int8_t y = g2[i + j];
            if (x >= 0 && y >= 0) {
                n++;
                sumX += x;
                sumY += y;
                sumXY += x * y;
                sumX2 += x * x;
                sumY2 += y * y;
            }
        }
    }
    for (; i < sz; i++) {
        int8_t x = g1[i];
        int8_t y = g2[i];
        if (x >= 0 && y >= 0) {
            n++;
            sumX += x;
            sumY += y;
            sumXY += x * y;
            sumX2 += x * x;
            sumY2 += y * y;
        }
    }

    if (n < 2) return 0.0;
    double meanX = static_cast<double>(sumX) / n;
    double meanY = static_cast<double>(sumY) / n;
    double cov = static_cast<double>(sumXY) / n - meanX * meanY;
    double varX = static_cast<double>(sumX2) / n - meanX * meanX;
    double varY = static_cast<double>(sumY2) / n - meanY * meanY;
    if (varX <= 0.0 || varY <= 0.0) return 0.0;
    double r = cov / (std::sqrt(varX) * std::sqrt(varY));
    return r * r;
}
#endif

// Fast r² using SIMD-accelerated function
static inline double computeRsqFast(const LDVariantOpt& v1, const LDVariantOpt& v2) {
    if (v1.genotype.size() != v2.genotype.size()) return 0.0;
    if (v1.varX <= 0.0 || v2.varX <= 0.0) return 0.0;
    return computeRsqSIMD(v1.genotype.data(), v2.genotype.data(), v1.genotype.size());
}

// ----------------------------------------------------------------------
// Display Help
// ----------------------------------------------------------------------
void VCFXLDCalculator::displayHelp() {
    std::cout << "VCFX_ld_calculator: Calculate pairwise LD (r^2) for variants in a VCF region.\n"
              << "Version 2.0 - Extreme-performance with mmap, SIMD, and multi-threading.\n\n"
              << "Usage:\n"
              << "  VCFX_ld_calculator [options] < input.vcf\n"
              << "  VCFX_ld_calculator [options] -i input.vcf\n\n"
              << "Options:\n"
              << "  -i, --input FILE          Input VCF file (uses memory-mapping for best performance)\n"
              << "  -r, --region <chr:s-e>    Only compute LD for variants in [start, end] on 'chr'\n"
              << "  -w, --window <N>          Window size in variants (default: 1000)\n"
              << "  -d, --max-distance <BP>   Max base-pair distance between pairs (0=unlimited)\n"
              << "  -t, --threshold <R2>      Only output pairs with r² >= threshold (default: 0.0)\n"
              << "  -n, --threads <N>         Number of threads (default: auto)\n"
              << "  -m, --matrix              Use matrix mode (MxM output) instead of streaming\n"
              << "                            WARNING: O(M²) time - avoid for >10K variants\n"
              << "  -q, --quiet               Suppress informational messages\n"
              << "  -h, --help                Show this help message\n"
              << "  -v, --version             Show program version\n\n"
              << "Modes:\n"
              << "  Default (streaming): Outputs LD pairs incrementally using a sliding window.\n"
              << "                       Memory: O(window * samples) - constant for any file size.\n"
              << "                       Time: O(M * window) - linear in variant count.\n"
              << "  Matrix mode:         Produces an MxM matrix of all pairwise r² values.\n"
              << "                       Memory: O(M * samples) where M is number of variants.\n"
              << "                       Time: O(M²) - avoid for >10K variants!\n\n"
              << "Performance:\n"
              << "  - Memory-mapped I/O: Use -i flag for extreme speed\n"
              << "  - SIMD-accelerated r² computation (NEON/AVX2/SSE2)\n"
              << "  - Multi-threaded matrix computation\n"
              << "  - Distance-based pruning with --max-distance\n\n"
              << "Example:\n"
              << "  # Fast streaming mode with file input\n"
              << "  VCFX_ld_calculator -i input.vcf -w 500 -t 0.2 > ld_pairs.txt\n\n"
              << "  # Streaming with distance limit (biology: LD decays with distance)\n"
              << "  VCFX_ld_calculator -i input.vcf --max-distance 500000 > ld_pairs.txt\n\n"
              << "  # Matrix mode (small regions only)\n"
              << "  VCFX_ld_calculator -i input.vcf -m -r chr1:10000-20000 > ld_matrix.txt\n";
}

// ----------------------------------------------------------------------
// parseRegion
// ----------------------------------------------------------------------
bool VCFXLDCalculator::parseRegion(const std::string &regionStr, std::string &regionChrom, int &regionStart,
                                   int &regionEnd) {
    auto colonPos = regionStr.find(':');
    if (colonPos == std::string::npos) return false;
    regionChrom = regionStr.substr(0, colonPos);
    auto dashPos = regionStr.find('-', colonPos + 1);
    if (dashPos == std::string::npos) return false;
    std::string startStr = regionStr.substr(colonPos + 1, dashPos - (colonPos + 1));
    std::string endStr = regionStr.substr(dashPos + 1);
    try {
        regionStart = std::stoi(startStr);
        regionEnd = std::stoi(endStr);
    } catch (...) { return false; }
    if (regionStart > regionEnd) return false;
    return true;
}

// ----------------------------------------------------------------------
// parse a genotype string (backward compat)
// ----------------------------------------------------------------------
int VCFXLDCalculator::parseGenotype(const std::string &s) {
    if (s.empty() || s == "." || s == "./." || s == ".|.") return -1;
    std::string g(s);
    for (char &c : g) { if (c == '|') c = '/'; }
    auto slashPos = g.find('/');
    if (slashPos == std::string::npos) return -1;
    auto a1 = g.substr(0, slashPos);
    auto a2 = g.substr(slashPos + 1);
    if (a1.empty() || a2.empty() || a1 == "." || a2 == ".") return -1;
    int i1 = 0, i2 = 0;
    try { i1 = std::stoi(a1); i2 = std::stoi(a2); } catch (...) { return -1; }
    if (i1 < 0 || i2 < 0 || i1 > 1 || i2 > 1) return -1;
    if (i1 == i2) return (i1 == 0) ? 0 : 2;
    return 1;
}

// ----------------------------------------------------------------------
// compute r^2 (backward compat)
// ----------------------------------------------------------------------
double VCFXLDCalculator::computeRsq(const std::vector<int> &g1, const std::vector<int> &g2) {
    if (g1.size() != g2.size()) return 0.0;
    int n = 0;
    long sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, sumY2 = 0;
    for (size_t i = 0; i < g1.size(); i++) {
        int x = g1[i], y = g2[i];
        if (x < 0 || y < 0) continue;
        n++;
        sumX += x; sumY += y; sumXY += x * y;
        sumX2 += (x * x); sumY2 += (y * y);
    }
    if (n < 2) return 0.0;
    double meanX = (double)sumX / n, meanY = (double)sumY / n;
    double cov = ((double)sumXY / n) - (meanX * meanY);
    double varX = ((double)sumX2 / n) - (meanX * meanX);
    double varY = ((double)sumY2 / n) - (meanY * meanY);
    if (varX <= 0.0 || varY <= 0.0) return 0.0;
    double r = cov / (std::sqrt(varX) * std::sqrt(varY));
    return r * r;
}

// ----------------------------------------------------------------------
// computeLDStreamingMmap: Extreme-performance streaming mode with mmap
// ----------------------------------------------------------------------
void VCFXLDCalculator::computeLDStreamingMmap(const char* data, size_t size, int outFd,
                                               const std::string &regionChrom, int regionStart,
                                               int regionEnd, size_t windowSize, double threshold,
                                               int maxDist) {
    OutputBuffer out(outFd);
    std::deque<LDVariantOpt> window;
    bool foundChromHeader = false;
    int numSamples = 0;
    char numBuf[32];

    // Output header
    out.append("#VAR1_CHROM\tVAR1_POS\tVAR1_ID\tVAR2_CHROM\tVAR2_POS\tVAR2_ID\tR2\n");

    const char* p = data;
    const char* end = data + size;

    while (p < end) {
        const char* lineStart = p;
        const char* lineEnd = findNewlineSIMD(p, end);
        size_t lineLen = lineEnd - lineStart;
        p = (lineEnd < end) ? lineEnd + 1 : end;

        if (lineLen == 0) continue;

        if (lineStart[0] == '#') {
            if (!foundChromHeader && lineLen >= 6 &&
                lineStart[0] == '#' && lineStart[1] == 'C' && lineStart[2] == 'H' &&
                lineStart[3] == 'R' && lineStart[4] == 'O' && lineStart[5] == 'M') {
                foundChromHeader = true;
                int tabCount = 0;
                for (size_t i = 0; i < lineLen; i++) {
                    if (lineStart[i] == '\t') tabCount++;
                }
                numSamples = (tabCount >= 9) ? (tabCount - 8) : 0;
            }
            continue;
        }

        if (!foundChromHeader) {
            if (!quiet) std::cerr << "Error: data line before #CHROM\n";
            break;
        }

        // Parse fields using zero-allocation approach
        const char* fieldStarts[11];
        size_t fieldLens[11];
        int fieldCount = 0;
        size_t start = 0;

        for (size_t i = 0; i <= lineLen && fieldCount < 10; i++) {
            if (i == lineLen || lineStart[i] == '\t') {
                fieldStarts[fieldCount] = lineStart + start;
                fieldLens[fieldCount] = i - start;
                fieldCount++;
                start = i + 1;
            }
        }

        if (fieldCount < 10) continue;

        int posVal;
        if (!fastParseInt(fieldStarts[1], fieldLens[1], posVal)) continue;

        // Check region
        if (!regionChrom.empty()) {
            if (fieldLens[0] != regionChrom.size() ||
                memcmp(fieldStarts[0], regionChrom.c_str(), fieldLens[0]) != 0) continue;
            if (posVal < regionStart || posVal > regionEnd) continue;
        }

        // Create variant
        LDVariantOpt v;
        v.chrom.assign(fieldStarts[0], fieldLens[0]);
        v.pos = posVal;

        if (fieldLens[2] == 1 && fieldStarts[2][0] == '.') {
            v.id = v.chrom + ":" + std::to_string(posVal);
        } else {
            v.id.assign(fieldStarts[2], fieldLens[2]);
        }

        v.genotype.resize(numSamples, -1);

        // Parse samples
        const char* sampleStart = fieldStarts[9];
        int sampleIdx = 0;

        while (sampleStart < lineStart + lineLen && sampleIdx < numSamples) {
            const char* sampleEnd = sampleStart;
            while (sampleEnd < lineStart + lineLen && *sampleEnd != '\t') sampleEnd++;
            size_t sampleLen = sampleEnd - sampleStart;

            const char* gtStart;
            size_t gtLen;
            if (extractGT(sampleStart, sampleLen, gtStart, gtLen)) {
                v.genotype[sampleIdx] = static_cast<int8_t>(parseGenotypeRaw(gtStart, gtLen));
            }

            sampleIdx++;
            sampleStart = sampleEnd + 1;
        }

        v.computeStats();

        // Compute LD with variants in window
        for (const auto &prev : window) {
            // Distance-based pruning
            if (maxDist > 0 && v.chrom == prev.chrom) {
                if (std::abs(v.pos - prev.pos) > maxDist) continue;
            }

            double r2 = computeRsqFast(prev, v);
            if (r2 >= threshold) {
                // Build output line
                out.append(prev.chrom);
                out.appendChar('\t');
                size_t len = formatInt(prev.pos, numBuf);
                out.append(numBuf, len);
                out.appendChar('\t');
                out.append(prev.id);
                out.appendChar('\t');
                out.append(v.chrom);
                out.appendChar('\t');
                len = formatInt(v.pos, numBuf);
                out.append(numBuf, len);
                out.appendChar('\t');
                out.append(v.id);
                out.appendChar('\t');
                len = formatR2(r2, numBuf);
                out.append(numBuf, len);
                out.appendChar('\n');
            }
        }

        window.push_back(std::move(v));
        if (window.size() > windowSize) window.pop_front();
    }
}

// ----------------------------------------------------------------------
// computeLDMatrixMmap: Optimized matrix mode with multi-threading
// ----------------------------------------------------------------------
void VCFXLDCalculator::computeLDMatrixMmap(const char* data, size_t size, int outFd,
                                            const std::string &regionChrom, int regionStart,
                                            int regionEnd, int nThreads) {
    OutputBuffer out(outFd);
    std::vector<LDVariantOpt> variants;
    bool foundChromHeader = false;
    int numSamples = 0;

    const char* p = data;
    const char* end = data + size;

    // First pass: parse all variants
    while (p < end) {
        const char* lineStart = p;
        const char* lineEnd = findNewlineSIMD(p, end);
        size_t lineLen = lineEnd - lineStart;
        p = (lineEnd < end) ? lineEnd + 1 : end;

        if (lineLen == 0) continue;

        if (lineStart[0] == '#') {
            // Pass through header lines
            out.append(lineStart, lineLen);
            out.appendChar('\n');

            if (!foundChromHeader && lineLen >= 6 &&
                lineStart[1] == 'C' && lineStart[2] == 'H' && lineStart[3] == 'R' &&
                lineStart[4] == 'O' && lineStart[5] == 'M') {
                foundChromHeader = true;
                int tabCount = 0;
                for (size_t i = 0; i < lineLen; i++) {
                    if (lineStart[i] == '\t') tabCount++;
                }
                numSamples = (tabCount >= 9) ? (tabCount - 8) : 0;
            }
            continue;
        }

        if (!foundChromHeader) {
            if (!quiet) std::cerr << "Error: data line before #CHROM\n";
            break;
        }

        // Parse fields
        const char* fieldStarts[11];
        size_t fieldLens[11];
        int fieldCount = 0;
        size_t start = 0;

        for (size_t i = 0; i <= lineLen && fieldCount < 10; i++) {
            if (i == lineLen || lineStart[i] == '\t') {
                fieldStarts[fieldCount] = lineStart + start;
                fieldLens[fieldCount] = i - start;
                fieldCount++;
                start = i + 1;
            }
        }

        if (fieldCount < 10) {
            out.append(lineStart, lineLen);
            out.appendChar('\n');
            continue;
        }

        int posVal;
        if (!fastParseInt(fieldStarts[1], fieldLens[1], posVal)) {
            out.append(lineStart, lineLen);
            out.appendChar('\n');
            continue;
        }

        // Check region
        if (!regionChrom.empty()) {
            std::string chrom(fieldStarts[0], fieldLens[0]);
            if (chrom != regionChrom || posVal < regionStart || posVal > regionEnd) {
                out.append(lineStart, lineLen);
                out.appendChar('\n');
                continue;
            }
        }

        LDVariantOpt v;
        v.chrom.assign(fieldStarts[0], fieldLens[0]);
        v.pos = posVal;
        v.id.assign(fieldStarts[2], fieldLens[2]);
        v.genotype.resize(numSamples, -1);

        // Parse samples
        const char* sampleStart = fieldStarts[9];
        int sampleIdx = 0;

        while (sampleStart < lineStart + lineLen && sampleIdx < numSamples) {
            const char* sampleEnd = sampleStart;
            while (sampleEnd < lineStart + lineLen && *sampleEnd != '\t') sampleEnd++;
            size_t sampleLen = sampleEnd - sampleStart;

            const char* gtStart;
            size_t gtLen;
            if (extractGT(sampleStart, sampleLen, gtStart, gtLen)) {
                v.genotype[sampleIdx] = static_cast<int8_t>(parseGenotypeRaw(gtStart, gtLen));
            }

            sampleIdx++;
            sampleStart = sampleEnd + 1;
        }

        v.computeStats();
        variants.push_back(std::move(v));

        // Pass through data line
        out.append(lineStart, lineLen);
        out.appendChar('\n');
    }

    size_t M = variants.size();
    if (M < 2) {
        out.append("#LD_MATRIX_START\n");
        out.append("No or only one variant in the region => no pairwise LD.\n");
        out.append("#LD_MATRIX_END\n");
        return;
    }

    out.append("#LD_MATRIX_START\n");

    // Header row
    out.append("Index/Var");
    char numBuf[32];
    for (size_t j = 0; j < M; j++) {
        out.appendChar('\t');
        out.append(variants[j].chrom);
        out.appendChar(':');
        size_t len = formatInt(variants[j].pos, numBuf);
        out.append(numBuf, len);
    }
    out.appendChar('\n');

    // Compute matrix rows (potentially in parallel)
    if (nThreads <= 1 || M < 100) {
        // Single-threaded
        for (size_t i = 0; i < M; i++) {
            out.append(variants[i].chrom);
            out.appendChar(':');
            size_t len = formatInt(variants[i].pos, numBuf);
            out.append(numBuf, len);

            for (size_t j = 0; j < M; j++) {
                out.appendChar('\t');
                if (i == j) {
                    out.append("1.0000", 6);
                } else {
                    double r2 = computeRsqFast(variants[i], variants[j]);
                    len = formatR2(r2, numBuf);
                    out.append(numBuf, len);
                }
            }
            out.appendChar('\n');
        }
    } else {
        // Multi-threaded: compute rows in parallel
        std::vector<std::string> rowBuffers(M);

        auto computeRow = [&](size_t i) {
            std::string& row = rowBuffers[i];
            row.reserve(M * 8);
            char buf[32];

            row += variants[i].chrom;
            row += ':';
            size_t len = formatInt(variants[i].pos, buf);
            row.append(buf, len);

            for (size_t j = 0; j < M; j++) {
                row += '\t';
                if (i == j) {
                    row.append("1.0000", 6);
                } else {
                    double r2 = computeRsqFast(variants[i], variants[j]);
                    len = formatR2(r2, buf);
                    row.append(buf, len);
                }
            }
            row += '\n';
        };

        std::vector<std::thread> threads;
        std::atomic<size_t> nextRow(0);

        for (int t = 0; t < nThreads; t++) {
            threads.emplace_back([&]() {
                while (true) {
                    size_t row = nextRow.fetch_add(1);
                    if (row >= M) break;
                    computeRow(row);
                }
            });
        }

        for (auto& t : threads) t.join();

        // Write rows in order
        for (size_t i = 0; i < M; i++) {
            out.append(rowBuffers[i]);
        }
    }

    out.append("#LD_MATRIX_END\n");
}

// ----------------------------------------------------------------------
// computeLDStreaming: Streaming mode (stdin fallback)
// ----------------------------------------------------------------------
void VCFXLDCalculator::computeLDStreaming(std::istream &in, std::ostream &out, const std::string &regionChrom,
                                          int regionStart, int regionEnd, size_t windowSz, double threshold) {
    bool foundChromHeader = false;
    int numSamples = 0;
    std::deque<LDVariantOpt> window;
    std::string line;
    std::string outputBuffer;
    outputBuffer.reserve(1024 * 1024);
    char numBuf[32];

    outputBuffer = "#VAR1_CHROM\tVAR1_POS\tVAR1_ID\tVAR2_CHROM\tVAR2_POS\tVAR2_ID\tR2\n";

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            if (!foundChromHeader && line.rfind("#CHROM", 0) == 0) {
                foundChromHeader = true;
                const char* p = line.c_str();
                int tabCount = 0;
                while (*p) { if (*p == '\t') tabCount++; p++; }
                numSamples = (tabCount >= 9) ? (tabCount - 8) : 0;
            }
            continue;
        }

        if (!foundChromHeader) {
            std::cerr << "Error: encountered data line before #CHROM.\n";
            break;
        }

        const char* linePtr = line.c_str();
        size_t lineLen = line.size();

        const char* fieldStarts[11];
        size_t fieldLens[11];
        int fieldCount = 0;
        size_t start = 0;

        for (size_t i = 0; i <= lineLen && fieldCount < 10; i++) {
            if (i == lineLen || linePtr[i] == '\t') {
                fieldStarts[fieldCount] = linePtr + start;
                fieldLens[fieldCount] = i - start;
                fieldCount++;
                start = i + 1;
            }
        }

        if (fieldCount < 10) continue;

        int posVal;
        if (!fastParseInt(fieldStarts[1], fieldLens[1], posVal)) continue;

        if (!regionChrom.empty()) {
            if (fieldLens[0] != regionChrom.size() ||
                memcmp(fieldStarts[0], regionChrom.c_str(), fieldLens[0]) != 0) continue;
            if (posVal < regionStart || posVal > regionEnd) continue;
        }

        LDVariantOpt v;
        v.chrom.assign(fieldStarts[0], fieldLens[0]);
        v.pos = posVal;

        if (fieldLens[2] == 1 && fieldStarts[2][0] == '.') {
            v.id = v.chrom + ":" + std::to_string(posVal);
        } else {
            v.id.assign(fieldStarts[2], fieldLens[2]);
        }

        v.genotype.resize(numSamples, -1);

        const char* sampleStart = fieldStarts[9];
        int sampleIdx = 0;

        while (sampleStart < linePtr + lineLen && sampleIdx < numSamples) {
            const char* sampleEnd = sampleStart;
            while (sampleEnd < linePtr + lineLen && *sampleEnd != '\t') sampleEnd++;
            size_t sampleLen = sampleEnd - sampleStart;

            const char* gtStart;
            size_t gtLen;
            if (extractGT(sampleStart, sampleLen, gtStart, gtLen)) {
                v.genotype[sampleIdx] = static_cast<int8_t>(parseGenotypeRaw(gtStart, gtLen));
            }

            sampleIdx++;
            sampleStart = sampleEnd + 1;
        }

        v.computeStats();

        for (const auto &prev : window) {
            double r2 = computeRsqFast(prev, v);
            if (r2 >= threshold) {
                outputBuffer += prev.chrom;
                outputBuffer += '\t';
                size_t len = formatInt(prev.pos, numBuf);
                outputBuffer.append(numBuf, len);
                outputBuffer += '\t';
                outputBuffer += prev.id;
                outputBuffer += '\t';
                outputBuffer += v.chrom;
                outputBuffer += '\t';
                len = formatInt(v.pos, numBuf);
                outputBuffer.append(numBuf, len);
                outputBuffer += '\t';
                outputBuffer += v.id;
                outputBuffer += '\t';
                len = formatR2(r2, numBuf);
                outputBuffer.append(numBuf, len);
                outputBuffer += '\n';

                if (outputBuffer.size() > 900000) {
                    out << outputBuffer;
                    outputBuffer.clear();
                }
            }
        }

        window.push_back(std::move(v));
        if (window.size() > windowSz) window.pop_front();
    }

    if (!outputBuffer.empty()) out << outputBuffer;
}

// ----------------------------------------------------------------------
// computeLD: Matrix mode (stdin fallback)
// ----------------------------------------------------------------------
void VCFXLDCalculator::computeLD(std::istream &in, std::ostream &out, const std::string &regionChrom, int regionStart,
                                 int regionEnd) {
    bool foundChromHeader = false;
    std::vector<std::string> sampleNames;
    std::vector<LDVariant> variants;
    int numSamples = 0;
    std::vector<std::string> fields;
    fields.reserve(16);

    std::string line;
    while (true) {
        auto pos = in.tellg();
        if (!std::getline(in, line)) break;
        if (line.empty()) { out << line << "\n"; continue; }
        if (line[0] == '#') {
            out << line << "\n";
            if (!foundChromHeader && line.rfind("#CHROM", 0) == 0) {
                foundChromHeader = true;
                vcfx::split_tabs(line, fields);
                for (size_t c = 9; c < fields.size(); c++) sampleNames.push_back(fields[c]);
                numSamples = sampleNames.size();
            }
            continue;
        }
        if (!foundChromHeader) {
            std::cerr << "Error: encountered data line before #CHROM.\n";
            break;
        }

        vcfx::split_tabs(line, fields);
        if (fields.size() < 10) { out << line << "\n"; continue; }

        std::string &chrom = fields[0];
        int posVal = 0;
        try { posVal = std::stoi(fields[1]); } catch (...) { out << line << "\n"; continue; }

        if (!regionChrom.empty()) {
            if (chrom != regionChrom || posVal < regionStart || posVal > regionEnd) {
                out << line << "\n";
                continue;
            }
        }

        LDVariant v;
        v.chrom = chrom;
        v.pos = posVal;
        v.id = fields[2];
        v.genotype.resize(numSamples, -1);

        int sampleCol = 9;
        for (int s = 0; s < numSamples; s++) {
            if ((size_t)(sampleCol + s) >= fields.size()) break;
            v.genotype[s] = parseGenotype(fields[sampleCol + s]);
        }
        variants.push_back(std::move(v));
        out << line << "\n";
    }

    size_t M = variants.size();
    if (M < 2) {
        out << "#LD_MATRIX_START\n";
        out << "No or only one variant in the region => no pairwise LD.\n";
        out << "#LD_MATRIX_END\n";
        return;
    }

    out << "#LD_MATRIX_START\n";
    out << "Index/Var";
    for (size_t j = 0; j < M; j++) {
        out << "\t" << variants[j].chrom << ":" << variants[j].pos;
    }
    out << "\n";

    for (size_t i = 0; i < M; i++) {
        out << variants[i].chrom << ":" << variants[i].pos;
        for (size_t j = 0; j < M; j++) {
            if (i == j) {
                out << "\t1.0000";
            } else {
                double r2 = computeRsq(variants[i].genotype, variants[j].genotype);
                out << "\t" << std::fixed << std::setprecision(4) << r2;
            }
        }
        out << "\n";
    }

    out << "#LD_MATRIX_END\n";
}

// ----------------------------------------------------------------------
// main run
// ----------------------------------------------------------------------
int VCFXLDCalculator::run(int argc, char *argv[]) {
    static struct option long_opts[] = {
        {"help", no_argument, 0, 'h'},
        {"version", no_argument, 0, 'v'},
        {"input", required_argument, 0, 'i'},
        {"region", required_argument, 0, 'r'},
        {"streaming", no_argument, 0, 's'},
        {"matrix", no_argument, 0, 'm'},
        {"window", required_argument, 0, 'w'},
        {"threshold", required_argument, 0, 't'},
        {"threads", required_argument, 0, 'n'},
        {"max-distance", required_argument, 0, 'd'},
        {"quiet", no_argument, 0, 'q'},
        {0, 0, 0, 0}
    };

    bool showHelp = false;
    std::string regionStr;

    while (true) {
        int c = getopt_long(argc, argv, "hvi:r:smw:t:n:d:q", long_opts, NULL);
        if (c == -1) break;
        switch (c) {
        case 'h': showHelp = true; break;
        case 'v':
            std::cout << "VCFX_ld_calculator v2.0\n";
            return 0;
        case 'i': inputFile = optarg; break;
        case 'r': regionStr = optarg; break;
        case 's': streamingMode = true; matrixMode = false; break;
        case 'm': matrixMode = true; streamingMode = false; break;
        case 'w':
            try {
                windowSize = std::stoul(optarg);
                if (windowSize == 0) windowSize = 1;
            } catch (...) {
                std::cerr << "Error: Invalid window size '" << optarg << "'\n";
                return 1;
            }
            break;
        case 't':
            try {
                ldThreshold = std::stod(optarg);
                if (ldThreshold < 0.0) ldThreshold = 0.0;
                if (ldThreshold > 1.0) ldThreshold = 1.0;
            } catch (...) {
                std::cerr << "Error: Invalid threshold '" << optarg << "'\n";
                return 1;
            }
            break;
        case 'n':
            try {
                numThreads = std::stoi(optarg);
                if (numThreads < 0) numThreads = 0;
            } catch (...) {
                std::cerr << "Error: Invalid thread count '" << optarg << "'\n";
                return 1;
            }
            break;
        case 'd':
            try {
                maxDistance = std::stoi(optarg);
                if (maxDistance < 0) maxDistance = 0;
            } catch (...) {
                std::cerr << "Error: Invalid max-distance '" << optarg << "'\n";
                return 1;
            }
            break;
        case 'q': quiet = true; break;
        default: showHelp = true;
        }
    }

    // Check for positional argument (input file)
    if (optind < argc && inputFile.empty()) {
        inputFile = argv[optind];
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    // Parse region
    std::string regionChrom;
    int regionStart = 0, regionEnd = 0;
    if (!regionStr.empty()) {
        if (!parseRegion(regionStr, regionChrom, regionStart, regionEnd)) {
            std::cerr << "Error parsing region '" << regionStr << "'. Use e.g. chr1:10000-20000\n";
            return 1;
        }
    }

    // Auto-detect thread count
    if (numThreads == 0) {
        numThreads = std::thread::hardware_concurrency();
        if (numThreads == 0) numThreads = 4;
    }

    // Use mmap if input file specified
    if (!inputFile.empty()) {
        MappedFile mf;
        if (!mf.open(inputFile.c_str())) {
            std::cerr << "Error: cannot open file '" << inputFile << "'\n";
            return 1;
        }

        int outFd = STDOUT_FILENO;

        if (matrixMode) {
            computeLDMatrixMmap(mf.data, mf.size, outFd, regionChrom, regionStart, regionEnd, numThreads);
        } else {
            computeLDStreamingMmap(mf.data, mf.size, outFd, regionChrom, regionStart, regionEnd,
                                   windowSize, ldThreshold, maxDistance);
        }
    } else {
        // Stdin fallback
        if (matrixMode) {
            computeLD(std::cin, std::cout, regionChrom, regionStart, regionEnd);
        } else {
            computeLDStreaming(std::cin, std::cout, regionChrom, regionStart, regionEnd, windowSize, ldThreshold);
        }
    }

    return 0;
}

static void show_help() {
    VCFXLDCalculator obj;
    char arg0[] = "VCFX_ld_calculator";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();
    if (vcfx::handle_common_flags(argc, argv, "VCFX_ld_calculator", show_help))
        return 0;
    VCFXLDCalculator calc;
    return calc.run(argc, argv);
}
