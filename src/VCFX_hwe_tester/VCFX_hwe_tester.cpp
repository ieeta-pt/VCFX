#include "VCFX_hwe_tester.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <string_view>
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
#elif defined(__aarch64__)
#include <arm_neon.h>
#define USE_NEON 1
#endif

// ============================================================================
// SIMD-OPTIMIZED LINE SCANNING
// ============================================================================

#if defined(USE_AVX2)
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

static inline const char* findTabSIMD(const char* p, const char* end) {
    const __m256i tab = _mm256_set1_epi8('\t');
    while (p + 32 <= end) {
        __m256i chunk = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p));
        __m256i cmp = _mm256_cmpeq_epi8(chunk, tab);
        int mask = _mm256_movemask_epi8(cmp);
        if (mask) {
            return p + __builtin_ctz(static_cast<unsigned int>(mask));
        }
        p += 32;
    }
    return static_cast<const char*>(memchr(p, '\t', static_cast<size_t>(end - p)));
}

#elif defined(USE_SSE2)
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

static inline const char* findTabSIMD(const char* p, const char* end) {
    const __m128i tab = _mm_set1_epi8('\t');
    while (p + 16 <= end) {
        __m128i chunk = _mm_loadu_si128(reinterpret_cast<const __m128i*>(p));
        __m128i cmp = _mm_cmpeq_epi8(chunk, tab);
        int mask = _mm_movemask_epi8(cmp);
        if (mask) {
            return p + __builtin_ctz(static_cast<unsigned int>(mask));
        }
        p += 16;
    }
    return static_cast<const char*>(memchr(p, '\t', static_cast<size_t>(end - p)));
}

#elif defined(USE_NEON)
static inline const char* findNewlineSIMD(const char* p, const char* end) {
    const uint8x16_t nl = vdupq_n_u8('\n');
    while (p + 16 <= end) {
        uint8x16_t chunk = vld1q_u8(reinterpret_cast<const uint8_t*>(p));
        uint8x16_t cmp = vceqq_u8(chunk, nl);
        uint64x2_t cmp64 = vreinterpretq_u64_u8(cmp);
        uint64_t lo = vgetq_lane_u64(cmp64, 0);
        uint64_t hi = vgetq_lane_u64(cmp64, 1);
        if (lo) return p + (__builtin_ctzll(lo) >> 3);
        if (hi) return p + 8 + (__builtin_ctzll(hi) >> 3);
        p += 16;
    }
    return static_cast<const char*>(memchr(p, '\n', static_cast<size_t>(end - p)));
}

static inline const char* findTabSIMD(const char* p, const char* end) {
    const uint8x16_t tab = vdupq_n_u8('\t');
    while (p + 16 <= end) {
        uint8x16_t chunk = vld1q_u8(reinterpret_cast<const uint8_t*>(p));
        uint8x16_t cmp = vceqq_u8(chunk, tab);
        uint64x2_t cmp64 = vreinterpretq_u64_u8(cmp);
        uint64_t lo = vgetq_lane_u64(cmp64, 0);
        uint64_t hi = vgetq_lane_u64(cmp64, 1);
        if (lo) return p + (__builtin_ctzll(lo) >> 3);
        if (hi) return p + 8 + (__builtin_ctzll(hi) >> 3);
        p += 16;
    }
    return static_cast<const char*>(memchr(p, '\t', static_cast<size_t>(end - p)));
}

#else
static inline const char* findNewlineSIMD(const char* p, const char* end) {
    return static_cast<const char*>(memchr(p, '\n', static_cast<size_t>(end - p)));
}

static inline const char* findTabSIMD(const char* p, const char* end) {
    return static_cast<const char*>(memchr(p, '\t', static_cast<size_t>(end - p)));
}
#endif

// ============================================================================
// MEMORY-MAPPED FILE SUPPORT
// ============================================================================

struct MappedFile {
    const char* data = nullptr;
    size_t size = 0;
    int fd = -1;

    bool open(const char* path) {
        fd = ::open(path, O_RDONLY);
        if (fd < 0) return false;

        struct stat st;
        if (fstat(fd, &st) < 0) {
            ::close(fd);
            fd = -1;
            return false;
        }

        size = static_cast<size_t>(st.st_size);
        if (size == 0) {
            ::close(fd);
            fd = -1;
            return true;
        }

        data = static_cast<const char*>(mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0));
        if (data == MAP_FAILED) {
            data = nullptr;
            ::close(fd);
            fd = -1;
            return false;
        }

        madvise(const_cast<char*>(data), size, MADV_SEQUENTIAL | MADV_WILLNEED);
        return true;
    }

    void close() {
        if (data && size > 0) {
            munmap(const_cast<char*>(data), size);
        }
        if (fd >= 0) {
            ::close(fd);
        }
        data = nullptr;
        size = 0;
        fd = -1;
    }

    ~MappedFile() { close(); }
};

// ============================================================================
// OUTPUT BUFFER
// ============================================================================

class OutputBuffer {
public:
    std::vector<char> buffer;
    size_t pos = 0;
    static constexpr size_t BUFFER_SIZE = 4 * 1024 * 1024;  // 4MB
    static constexpr size_t FLUSH_THRESHOLD = 3 * 1024 * 1024;

    OutputBuffer() { buffer.resize(BUFFER_SIZE); }

    void flush() {
        if (pos > 0) {
            ::write(STDOUT_FILENO, buffer.data(), pos);
            pos = 0;
        }
    }

    void maybeFlush() {
        if (pos >= FLUSH_THRESHOLD) {
            flush();
        }
    }

    void ensureSpace(size_t needed) {
        if (pos + needed > buffer.size()) {
            flush();
        }
    }

    void append(const char* s, size_t len) {
        ensureSpace(len);
        memcpy(buffer.data() + pos, s, len);
        pos += len;
    }

    void append(std::string_view sv) {
        append(sv.data(), sv.size());
    }

    void appendChar(char c) {
        ensureSpace(1);
        buffer[pos++] = c;
    }

    // Fast double formatting with 6 decimal places
    void appendDouble(double val) {
        ensureSpace(32);
        if (val < 0) {
            buffer[pos++] = '-';
            val = -val;
        }

        long long intPart = static_cast<long long>(val);
        double fracPart = val - intPart;

        char intBuf[20];
        int intLen = 0;
        if (intPart == 0) {
            intBuf[intLen++] = '0';
        } else {
            while (intPart > 0) {
                intBuf[intLen++] = '0' + (intPart % 10);
                intPart /= 10;
            }
        }
        for (int i = intLen - 1; i >= 0; i--) {
            buffer[pos++] = intBuf[i];
        }

        buffer[pos++] = '.';

        for (int i = 0; i < 6; i++) {
            fracPart *= 10.0;
            int digit = static_cast<int>(fracPart);
            buffer[pos++] = '0' + digit;
            fracPart -= digit;
        }
    }

    ~OutputBuffer() { flush(); }
};

// ============================================================================
// CHI-SQUARE HWE CALCULATION
// ============================================================================

// Chi-square p-value for 1 degree of freedom using Abramowitz & Stegun approximation
static inline double chi2_pvalue_1df(double chi2) {
    if (chi2 <= 0.0) return 1.0;
    if (chi2 > 700.0) return 0.0;

    double x = std::sqrt(chi2 * 0.5);
    double t = 1.0 / (1.0 + 0.3275911 * x);
    double y = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 +
               t * (-1.453152027 + t * 1.061405429))));
    return y * std::exp(-x * x);
}

// Hardy-Weinberg chi-square test with Yates' continuity correction
static inline double calculateHWE_chisq(int homRef, int het, int homAlt) {
    int N = homRef + het + homAlt;
    if (N < 1) return 1.0;

    double p = (2.0 * homRef + het) / (2.0 * N);
    double q = 1.0 - p;

    if (p <= 0.0 || p >= 1.0) return 1.0;

    double exp_homRef = N * p * p;
    double exp_het = N * 2.0 * p * q;
    double exp_homAlt = N * q * q;

    auto yates_term = [](double obs, double exp) -> double {
        if (exp <= 0.0) return 0.0;
        double diff = std::abs(obs - exp) - 0.5;
        if (diff < 0.0) diff = 0.0;
        return (diff * diff) / exp;
    };

    double chi2 = yates_term(homRef, exp_homRef) +
                  yates_term(het, exp_het) +
                  yates_term(homAlt, exp_homAlt);

    return chi2_pvalue_1df(chi2);
}

// ============================================================================
// ZERO-COPY PARSING HELPERS
// ============================================================================

static inline const char* skipToField(const char* p, const char* lineEnd, int fieldIdx) {
    for (int i = 0; i < fieldIdx && p < lineEnd; i++) {
        const char* tab = findTabSIMD(p, lineEnd);
        if (!tab) return nullptr;
        p = tab + 1;
    }
    return p < lineEnd ? p : nullptr;
}

static inline std::string_view getField(const char* lineStart, const char* lineEnd, int fieldIdx) {
    const char* p = skipToField(lineStart, lineEnd, fieldIdx);
    if (!p) return {};
    const char* tab = findTabSIMD(p, lineEnd);
    if (!tab) tab = lineEnd;
    return std::string_view(p, static_cast<size_t>(tab - p));
}

// Parse genotype and return: 0=homRef, 1=het, 2=homAlt, -1=invalid/missing
static inline int parseGenotypeForHWE(const char* p, const char* end) {
    if (p >= end) return -1;

    // Handle colon-separated FORMAT fields - only parse GT (first field)
    const char* colonPos = static_cast<const char*>(memchr(p, ':', static_cast<size_t>(end - p)));
    if (colonPos) end = colonPos;

    while (p < end && (*p == ' ' || *p == '\r')) p++;
    if (p >= end) return -1;

    if (*p == '.') return -1;

    // Parse first allele
    if (*p < '0' || *p > '9') return -1;
    int a1 = 0;
    while (p < end && *p >= '0' && *p <= '9') {
        a1 = a1 * 10 + (*p - '0');
        p++;
    }

    if (p >= end || (*p != '/' && *p != '|')) return -1;
    p++;

    if (p >= end || *p == '.') return -1;

    if (*p < '0' || *p > '9') return -1;
    int a2 = 0;
    while (p < end && *p >= '0' && *p <= '9') {
        a2 = a2 * 10 + (*p - '0');
        p++;
    }

    // Only handle biallelic (0 or 1)
    if (a1 > 1 || a2 > 1) return -1;

    // 0/0 => homRef(0), 1/1 => homAlt(2), 0/1 or 1/0 => het(1)
    if (a1 == 0 && a2 == 0) return 0;
    if (a1 == 1 && a2 == 1) return 2;
    return 1;  // het
}

static inline bool isBiallelic(std::string_view alt) {
    return alt.find(',') == std::string_view::npos;
}

// ============================================================================
// COMMAND-LINE ARGUMENTS
// ============================================================================

struct HWEArgs {
    const char* inputFile = nullptr;
    bool quiet = false;
    bool showHelp = false;
};

static void displayHelp() {
    std::cout << "VCFX_hwe_tester: Perform Hardy-Weinberg Equilibrium (HWE) tests on a biallelic VCF.\n\n"
              << "Usage:\n"
              << "  VCFX_hwe_tester [options] [input.vcf]\n"
              << "  VCFX_hwe_tester [options] < input.vcf\n\n"
              << "Options:\n"
              << "  -i, --input FILE   Input VCF file (uses memory-mapping for best performance)\n"
              << "  -q, --quiet        Suppress informational messages\n"
              << "  -h, --help         Show this help.\n\n"
              << "Description:\n"
              << "  Reads each variant line, ignoring multi-allelic calls. For biallelic lines,\n"
              << "  collects genotypes as 0/0, 0/1, 1/1, then uses chi-square test with Yates'\n"
              << "  continuity correction to produce a p-value for HWE.\n\n"
              << "Performance:\n"
              << "  Uses memory-mapped I/O and SIMD for ~20x speedup over stdin mode.\n\n"
              << "Example:\n"
              << "  VCFX_hwe_tester -i input.vcf > results.txt\n"
              << "  VCFX_hwe_tester < input.vcf > results.txt\n";
}

static HWEArgs parseArgs(int argc, char* argv[]) {
    HWEArgs args;

    static struct option long_opts[] = {
        {"help", no_argument, 0, 'h'},
        {"input", required_argument, 0, 'i'},
        {"quiet", no_argument, 0, 'q'},
        {0, 0, 0, 0}
    };

    optind = 1;
    while (true) {
        int c = getopt_long(argc, argv, "hi:q", long_opts, nullptr);
        if (c == -1) break;

        switch (c) {
        case 'h':
            args.showHelp = true;
            break;
        case 'i':
            args.inputFile = optarg;
            break;
        case 'q':
            args.quiet = true;
            break;
        default:
            args.showHelp = true;
        }
    }

    if (!args.inputFile && optind < argc) {
        args.inputFile = argv[optind];
    }

    return args;
}

// ============================================================================
// MAIN PROCESSING - MMAP MODE
// ============================================================================

static void performHWE_Mmap(const MappedFile& file, OutputBuffer& outBuf) {
    if (file.size == 0) return;

    const char* data = file.data;
    const char* end = data + file.size;
    const char* p = data;

    // Output header
    outBuf.append("CHROM\tPOS\tID\tREF\tALT\tHWE_pvalue\n");

    // Skip header lines
    while (p < end) {
        const char* lineEnd = findNewlineSIMD(p, end);
        if (!lineEnd) lineEnd = end;

        if (*p != '#') break;
        p = lineEnd + 1;
    }

    // Process data lines
    while (p < end) {
        const char* lineEnd = findNewlineSIMD(p, end);
        if (!lineEnd) lineEnd = end;

        size_t lineLen = static_cast<size_t>(lineEnd - p);
        if (lineLen > 0 && p[lineLen - 1] == '\r') lineLen--;

        if (lineLen == 0 || *p == '#') {
            p = lineEnd + 1;
            continue;
        }

        const char* lineStart = p;
        const char* adjLineEnd = p + lineLen;

        // Get fields 0-4 (CHROM, POS, ID, REF, ALT)
        std::string_view chrom = getField(lineStart, adjLineEnd, 0);
        std::string_view pos = getField(lineStart, adjLineEnd, 1);
        std::string_view id = getField(lineStart, adjLineEnd, 2);
        std::string_view ref = getField(lineStart, adjLineEnd, 3);
        std::string_view alt = getField(lineStart, adjLineEnd, 4);

        if (chrom.empty() || pos.empty() || alt.empty()) {
            p = lineEnd + 1;
            continue;
        }

        // Skip multi-allelic
        if (!isBiallelic(alt)) {
            p = lineEnd + 1;
            continue;
        }

        // Check FORMAT field (field 8) starts with "GT"
        std::string_view format = getField(lineStart, adjLineEnd, 8);
        if (format.empty() || format.substr(0, 2) != "GT") {
            p = lineEnd + 1;
            continue;
        }

        // Skip to field 9 (first sample)
        const char* sampleStart = skipToField(lineStart, adjLineEnd, 9);
        if (!sampleStart) {
            p = lineEnd + 1;
            continue;
        }

        // Count genotypes
        int homRef = 0, het = 0, homAlt = 0;
        const char* sp = sampleStart;

        while (sp < adjLineEnd) {
            const char* nextTab = findTabSIMD(sp, adjLineEnd);
            if (!nextTab) nextTab = adjLineEnd;

            int gt = parseGenotypeForHWE(sp, nextTab);
            if (gt == 0) homRef++;
            else if (gt == 1) het++;
            else if (gt == 2) homAlt++;

            sp = nextTab + 1;
        }

        // Calculate HWE p-value
        double pVal = calculateHWE_chisq(homRef, het, homAlt);

        // Output result
        outBuf.append(chrom);
        outBuf.appendChar('\t');
        outBuf.append(pos);
        outBuf.appendChar('\t');
        outBuf.append(id);
        outBuf.appendChar('\t');
        outBuf.append(ref);
        outBuf.appendChar('\t');
        outBuf.append(alt);
        outBuf.appendChar('\t');
        outBuf.appendDouble(pVal);
        outBuf.appendChar('\n');

        outBuf.maybeFlush();

        p = lineEnd + 1;
    }
}

// ============================================================================
// MAIN PROCESSING - STDIN MODE
// ============================================================================

static void performHWE_Stdin(std::istream& in, OutputBuffer& outBuf) {
    std::cout << "CHROM\tPOS\tID\tREF\tALT\tHWE_pvalue\n";

    std::string line;
    std::vector<std::string> fields;
    fields.reserve(16);

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;

        vcfx::split_tabs(line, fields);
        if (fields.size() < 10) continue;

        const std::string& chrom = fields[0];
        const std::string& pos = fields[1];
        const std::string& id = fields[2];
        const std::string& ref = fields[3];
        const std::string& alt = fields[4];

        if (alt.find(',') != std::string::npos) continue;

        // Check FORMAT field (field 8) starts with "GT"
        const std::string& format = fields[8];
        if (format.empty() || format.substr(0, 2) != "GT") continue;

        // Count genotypes
        int homRef = 0, het = 0, homAlt = 0;

        for (size_t s = 9; s < fields.size(); s++) {
            const std::string& sample = fields[s];
            int gt = parseGenotypeForHWE(sample.data(), sample.data() + sample.size());
            if (gt == 0) homRef++;
            else if (gt == 1) het++;
            else if (gt == 2) homAlt++;
        }

        double pVal = calculateHWE_chisq(homRef, het, homAlt);

        std::cout << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt
                  << "\t" << std::fixed << std::setprecision(6) << pVal << "\n";
    }
}

// ============================================================================
// LEGACY CLASS INTERFACE
// ============================================================================

int VCFXHWETester::run(int argc, char *argv[]) {
    HWEArgs args = parseArgs(argc, argv);

    if (args.showHelp) {
        displayHelp();
        return 0;
    }

    if (args.inputFile) {
        MappedFile file;
        if (!file.open(args.inputFile)) {
            std::cerr << "Error: Cannot open file: " << args.inputFile << "\n";
            return 1;
        }

        if (!args.quiet) {
            std::cerr << "Processing " << args.inputFile << " (" << file.size << " bytes)...\n";
        }

        OutputBuffer outBuf;
        performHWE_Mmap(file, outBuf);
    } else {
        OutputBuffer outBuf;
        performHWE_Stdin(std::cin, outBuf);
    }

    return 0;
}

void VCFXHWETester::displayHelp() {
    ::displayHelp();
}

bool VCFXHWETester::isBiallelic(const std::string &alt) {
    return alt.find(',') == std::string::npos;
}

bool VCFXHWETester::parseGenotypes(const std::vector<std::string> &genotypes, int &homRef, int &het, int &homAlt) {
    homRef = 0;
    het = 0;
    homAlt = 0;
    for (const auto &gt : genotypes) {
        int code = parseGenotypeForHWE(gt.data(), gt.data() + gt.size());
        if (code == 0) homRef++;
        else if (code == 1) het++;
        else if (code == 2) homAlt++;
        else if (code == -1 && (gt.find('2') != std::string::npos || gt.find('3') != std::string::npos)) {
            return false;  // Multi-allelic genotype
        }
    }
    return true;
}

double VCFXHWETester::genotypeProbability(int homRef, int het, int homAlt) {
    return calculateHWE_chisq(homRef, het, homAlt);
}

double VCFXHWETester::calculateHWE(int homRef, int het, int homAlt) {
    return calculateHWE_chisq(homRef, het, homAlt);
}

void VCFXHWETester::performHWE(std::istream &in) {
    OutputBuffer outBuf;
    performHWE_Stdin(in, outBuf);
}

// ============================================================================
// MAIN ENTRY POINT
// ============================================================================

static void show_help() {
    ::displayHelp();
}

int main(int argc, char *argv[]) {
    vcfx::init_io();
    if (vcfx::handle_common_flags(argc, argv, "VCFX_hwe_tester", show_help))
        return 0;
    VCFXHWETester tester;
    return tester.run(argc, argv);
}
