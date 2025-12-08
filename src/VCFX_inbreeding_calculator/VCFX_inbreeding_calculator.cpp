#include "VCFX_inbreeding_calculator.h"
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
            return true;  // Empty file is valid
        }

        data = static_cast<const char*>(mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0));
        if (data == MAP_FAILED) {
            data = nullptr;
            ::close(fd);
            fd = -1;
            return false;
        }

        // Advise kernel for sequential access
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
    static constexpr size_t BUFFER_SIZE = 1024 * 1024;  // 1MB

    OutputBuffer() { buffer.resize(BUFFER_SIZE); }

    void flush() {
        if (pos > 0) {
            ::write(STDOUT_FILENO, buffer.data(), pos);
            pos = 0;
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

        // Handle integer part
        long long intPart = static_cast<long long>(val);
        double fracPart = val - intPart;

        // Format integer part
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
        // Reverse and copy
        for (int i = intLen - 1; i >= 0; i--) {
            buffer[pos++] = intBuf[i];
        }

        buffer[pos++] = '.';

        // Format 6 decimal places
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
// ZERO-COPY PARSING HELPERS
// ============================================================================

// Skip N tabs and return pointer to start of field N (0-indexed)
static inline const char* skipToField(const char* p, const char* lineEnd, int fieldIdx) {
    for (int i = 0; i < fieldIdx && p < lineEnd; i++) {
        const char* tab = findTabSIMD(p, lineEnd);
        if (!tab) return nullptr;
        p = tab + 1;
    }
    return p < lineEnd ? p : nullptr;
}

// Get field as string_view without allocation
static inline std::string_view getField(const char* lineStart, const char* lineEnd, int fieldIdx) {
    const char* p = skipToField(lineStart, lineEnd, fieldIdx);
    if (!p) return {};
    const char* tab = findTabSIMD(p, lineEnd);
    if (!tab) tab = lineEnd;
    return std::string_view(p, static_cast<size_t>(tab - p));
}

// Parse genotype code: 0/0=>0, 0/1=>1, 1/1=>2, else=>-1
// Fast inline parsing with zero allocation
static inline int parseGenotypeCode(const char* p, const char* end) {
    if (p >= end) return -1;

    // Handle colon-separated FORMAT fields - only parse GT (first field)
    const char* colonPos = static_cast<const char*>(memchr(p, ':', static_cast<size_t>(end - p)));
    if (colonPos) end = colonPos;

    // Skip leading whitespace/CR
    while (p < end && (*p == ' ' || *p == '\r')) p++;
    if (p >= end) return -1;

    // Check for missing
    if (*p == '.') return -1;

    // Parse first allele
    if (*p < '0' || *p > '9') return -1;
    int a1 = 0;
    while (p < end && *p >= '0' && *p <= '9') {
        a1 = a1 * 10 + (*p - '0');
        p++;
    }

    // Must have separator
    if (p >= end || (*p != '/' && *p != '|')) return -1;
    p++;

    // Check for missing second allele
    if (p >= end || *p == '.') return -1;

    // Parse second allele
    if (*p < '0' || *p > '9') return -1;
    int a2 = 0;
    while (p < end && *p >= '0' && *p <= '9') {
        a2 = a2 * 10 + (*p - '0');
        p++;
    }

    // Only handle biallelic (0 or 1)
    if (a1 > 1 || a2 > 1) return -1;

    // 0/0=>0, 1/1=>2, 0/1 or 1/0=>1
    if (a1 == a2) return (a1 == 0) ? 0 : 2;
    return 1;
}

// Check if ALT is biallelic (no comma)
static inline bool isBiallelic(std::string_view alt) {
    return alt.find(',') == std::string_view::npos;
}

// ============================================================================
// COMMAND-LINE ARGUMENTS
// ============================================================================

struct InbreedingArgs {
    const char* inputFile = nullptr;
    FrequencyMode freqMode = FrequencyMode::EXCLUDE_SAMPLE;
    bool skipBoundary = false;
    bool countBoundaryAsUsed = false;
    bool quiet = false;
    bool showHelp = false;
};

static void displayHelp() {
    std::cout << "VCFX_inbreeding_calculator: Compute individual inbreeding coefficients (F)\n"
              << "based on biallelic sites in a VCF.\n\n"
              << "Usage:\n"
              << "  VCFX_inbreeding_calculator [options] [input.vcf]\n"
              << "  VCFX_inbreeding_calculator [options] < input.vcf\n\n"
              << "Options:\n"
              << "  -i, --input FILE          Input VCF file (uses memory-mapping for best performance)\n"
              << "  -q, --quiet               Suppress informational messages\n"
              << "  -h, --help                Show this help.\n"
              << "  --freq-mode <mode>        'excludeSample' (default) or 'global'\n"
              << "  --skip-boundary           Skip boundary freq sites. By default, they are used.\n"
              << "  --count-boundary-as-used  If also skipping boundary, still increment usedCount.\n\n"
              << "Description:\n"
              << "  Reads a VCF in a single pass, ignoring multi-allelic lines (ALT with commas).\n"
              << "  For each biallelic variant, we parse each sample's genotype code:\n"
              << "       0/0 => 0,   0/1 => 1,   1/1 => 2, else => -1 (ignored)\n\n"
              << "  Then, depending on --freq-mode:\n"
              << "    * excludeSample => Each sample excludes its own genotype when computing p.\n"
              << "    * global        => Compute a single global p from all samples' genotypes.\n\n"
              << "  The --skip-boundary option, if set, ignores boundary freq p=0 or p=1.\n"
              << "    BUT if you also specify --count-boundary-as-used, those boundary sites\n"
              << "    increment usedCount (forcing F=1) without contributing to sumExp.\n\n"
              << "  If sumExp=0 for a sample but usedCount>0, we output F=1.\n"
              << "  If usedCount=0, we output NA.\n\n"
              << "Performance:\n"
              << "  Uses memory-mapped I/O and SIMD for ~20x speedup over stdin mode.\n"
              << "  When a file is provided directly, uses mmap for faster processing.\n\n"
              << "Example:\n"
              << "  VCFX_inbreeding_calculator -i input.vcf > inbreeding.txt\n"
              << "  VCFX_inbreeding_calculator < input.vcf > inbreeding.txt\n";
}

static InbreedingArgs parseArgs(int argc, char* argv[]) {
    InbreedingArgs args;

    static struct option long_opts[] = {
        {"help", no_argument, 0, 'h'},
        {"input", required_argument, 0, 'i'},
        {"quiet", no_argument, 0, 'q'},
        {"freq-mode", required_argument, 0, 0},
        {"skip-boundary", no_argument, 0, 0},
        {"count-boundary-as-used", no_argument, 0, 0},
        {0, 0, 0, 0}
    };

    optind = 1;
    while (true) {
        int option_index = 0;
        int c = getopt_long(argc, argv, "hi:q", long_opts, &option_index);
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
        case 0: {
            std::string optName(long_opts[option_index].name);
            if (optName == "freq-mode") {
                if (strcmp(optarg, "global") == 0) {
                    args.freqMode = FrequencyMode::GLOBAL;
                } else if (strcmp(optarg, "excludeSample") == 0) {
                    args.freqMode = FrequencyMode::EXCLUDE_SAMPLE;
                } else {
                    std::cerr << "Warning: unrecognized freq-mode='" << optarg
                              << "'. Using 'excludeSample' by default.\n";
                }
            } else if (optName == "skip-boundary") {
                args.skipBoundary = true;
            } else if (optName == "count-boundary-as-used") {
                args.countBoundaryAsUsed = true;
            }
            break;
        }
        default:
            args.showHelp = true;
        }
    }

    // Check for positional argument (input file)
    if (!args.inputFile && optind < argc) {
        args.inputFile = argv[optind];
    }

    return args;
}

// ============================================================================
// MAIN PROCESSING - MMAP MODE
// ============================================================================

static void calculateInbreedingMmap(const MappedFile& file, const InbreedingArgs& args,
                                     OutputBuffer& outBuf) {
    if (file.size == 0) {
        std::cerr << "Error: Empty file.\n";
        outBuf.append("Sample\tInbreedingCoefficient\n");
        return;
    }

    const char* data = file.data;
    const char* end = data + file.size;
    const char* p = data;

    // Parse header to get sample names
    std::vector<std::string> sampleNames;
    int numSamples = 0;
    bool foundChrom = false;

    while (p < end) {
        const char* lineEnd = findNewlineSIMD(p, end);
        if (!lineEnd) lineEnd = end;

        size_t lineLen = static_cast<size_t>(lineEnd - p);
        if (lineLen > 0 && p[lineLen - 1] == '\r') lineLen--;

        if (lineLen == 0) {
            p = lineEnd + 1;
            continue;
        }

        if (*p == '#') {
            // Check for #CHROM header line
            if (lineLen >= 6 && memcmp(p, "#CHROM", 6) == 0) {
                foundChrom = true;
                // Parse sample names from columns 9+
                const char* fieldStart = p;
                const char* fieldEnd = p + lineLen;
                int fieldIdx = 0;

                while (fieldStart < fieldEnd) {
                    const char* tab = findTabSIMD(fieldStart, fieldEnd);
                    if (!tab) tab = fieldEnd;

                    if (fieldIdx >= 9) {
                        std::string_view name(fieldStart, static_cast<size_t>(tab - fieldStart));
                        // Remove trailing CR if present
                        if (!name.empty() && name.back() == '\r') {
                            name = name.substr(0, name.size() - 1);
                        }
                        sampleNames.emplace_back(name);
                    }

                    fieldIdx++;
                    fieldStart = tab + 1;
                }
                numSamples = static_cast<int>(sampleNames.size());
            }
            p = lineEnd + 1;
            continue;
        }

        // Found first data line
        break;
    }

    if (!foundChrom || numSamples == 0) {
        std::cerr << "Error: No #CHROM line or no samples found.\n";
        outBuf.append("Sample\tInbreedingCoefficient\n");
        return;
    }

    // Per-sample accumulators
    std::vector<double> sumExp(numSamples, 0.0);
    std::vector<double> obsHet(numSamples, 0.0);
    std::vector<int> usedCount(numSamples, 0);

    // Reusable buffer for genotype codes
    std::vector<int> genotypeCodes(numSamples);
    int variantCount = 0;

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

        // Get ALT field (index 4) and check biallelic
        std::string_view alt = getField(lineStart, adjLineEnd, 4);
        if (alt.empty() || !isBiallelic(alt)) {
            p = lineEnd + 1;
            continue;
        }

        // Parse genotypes for this variant
        // Skip to field 9 (first sample)
        const char* sampleStart = skipToField(lineStart, adjLineEnd, 9);
        if (!sampleStart) {
            p = lineEnd + 1;
            continue;
        }

        int altSum = 0;
        int nGood = 0;

        const char* sp = sampleStart;
        for (int s = 0; s < numSamples && sp < adjLineEnd; s++) {
            const char* nextTab = findTabSIMD(sp, adjLineEnd);
            if (!nextTab) nextTab = adjLineEnd;

            int code = parseGenotypeCode(sp, nextTab);
            genotypeCodes[s] = code;

            if (code >= 0) {
                altSum += code;
                nGood++;
            }

            sp = nextTab + 1;
        }

        // Skip if fewer than 2 genotyped samples
        if (nGood < 2) {
            p = lineEnd + 1;
            continue;
        }
        variantCount++;

        // For freq-mode=GLOBAL, compute p once
        double globalP = static_cast<double>(altSum) / (2.0 * nGood);

        // Process each sample for this variant
        for (int s = 0; s < numSamples; s++) {
            int code = genotypeCodes[s];
            if (code < 0) continue;

            double freq;
            if (args.freqMode == FrequencyMode::GLOBAL) {
                freq = globalP;
            } else {
                // EXCLUDE_SAMPLE mode
                int altEx = altSum - code;
                int validEx = nGood - 1;
                if (validEx < 1) continue;
                freq = static_cast<double>(altEx) / (2.0 * validEx);
            }

            // Handle boundary frequencies
            if (args.skipBoundary && (freq <= 0.0 || freq >= 1.0)) {
                if (args.countBoundaryAsUsed) {
                    usedCount[s]++;
                }
                continue;
            }

            usedCount[s]++;
            double eHet = 2.0 * freq * (1.0 - freq);
            sumExp[s] += eHet;

            if (code == 1) {
                obsHet[s] += 1.0;
            }
        }

        p = lineEnd + 1;
    }

    // Output results
    outBuf.append("Sample\tInbreedingCoefficient\n");

    if (variantCount == 0) {
        if (!args.quiet) {
            std::cerr << "No biallelic variants found.\n";
        }
        for (int s = 0; s < numSamples; s++) {
            outBuf.append(sampleNames[s]);
            outBuf.append("\tNA\n");
        }
        return;
    }

    for (int s = 0; s < numSamples; s++) {
        outBuf.append(sampleNames[s]);
        outBuf.appendChar('\t');

        if (usedCount[s] == 0) {
            outBuf.append("NA\n");
            continue;
        }

        double e = sumExp[s];
        if (e <= 0.0) {
            outBuf.append("1.000000\n");
            continue;
        }

        double f = 1.0 - (obsHet[s] / e);
        outBuf.appendDouble(f);
        outBuf.appendChar('\n');
    }
}

// ============================================================================
// MAIN PROCESSING - STDIN MODE
// ============================================================================

static void calculateInbreedingStdin(std::istream& in, const InbreedingArgs& args,
                                      OutputBuffer& outBuf) {
    std::vector<std::string> sampleNames;
    int numSamples = 0;
    bool foundChrom = false;

    // Per-sample accumulators
    std::vector<double> sumExp;
    std::vector<double> obsHet;
    std::vector<int> usedCount;
    std::vector<int> genotypeCodes;

    std::vector<std::string> fields;
    fields.reserve(16);

    int variantCount = 0;
    std::string line;
    line.reserve(65536);

    while (std::getline(in, line)) {
        // Skip empty lines
        if (line.empty()) continue;

        // Remove trailing CR
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) continue;

        if (line[0] == '#') {
            // Check for #CHROM header
            if (!foundChrom && line.find("#CHROM") != std::string::npos) {
                foundChrom = true;
                vcfx::split_tabs(line, fields);
                for (size_t i = 9; i < fields.size(); i++) {
                    sampleNames.push_back(fields[i]);
                }
                numSamples = static_cast<int>(sampleNames.size());
                sumExp.resize(numSamples, 0.0);
                obsHet.resize(numSamples, 0.0);
                usedCount.resize(numSamples, 0);
                genotypeCodes.resize(numSamples);
            }
            continue;
        }

        if (!foundChrom) continue;

        // Parse variant line
        vcfx::split_tabs(line, fields);
        if (fields.size() < 10) continue;

        // Check biallelic
        if (fields[4].find(',') != std::string::npos) continue;

        // Parse genotypes
        int altSum = 0;
        int nGood = 0;

        for (int s = 0; s < numSamples; s++) {
            int colIndex = 9 + s;
            if (static_cast<size_t>(colIndex) >= fields.size()) {
                genotypeCodes[s] = -1;
            } else {
                const std::string& gt = fields[colIndex];
                genotypeCodes[s] = parseGenotypeCode(gt.data(), gt.data() + gt.size());
            }
            if (genotypeCodes[s] >= 0) {
                altSum += genotypeCodes[s];
                nGood++;
            }
        }

        if (nGood < 2) continue;
        variantCount++;

        double globalP = static_cast<double>(altSum) / (2.0 * nGood);

        for (int s = 0; s < numSamples; s++) {
            int code = genotypeCodes[s];
            if (code < 0) continue;

            double freq;
            if (args.freqMode == FrequencyMode::GLOBAL) {
                freq = globalP;
            } else {
                int altEx = altSum - code;
                int validEx = nGood - 1;
                if (validEx < 1) continue;
                freq = static_cast<double>(altEx) / (2.0 * validEx);
            }

            if (args.skipBoundary && (freq <= 0.0 || freq >= 1.0)) {
                if (args.countBoundaryAsUsed) {
                    usedCount[s]++;
                }
                continue;
            }

            usedCount[s]++;
            double eHet = 2.0 * freq * (1.0 - freq);
            sumExp[s] += eHet;

            if (code == 1) {
                obsHet[s] += 1.0;
            }
        }
    }

    // Output results
    outBuf.append("Sample\tInbreedingCoefficient\n");

    if (!foundChrom) {
        std::cerr << "Error: No #CHROM line found.\n";
        return;
    }

    if (numSamples == 0) {
        std::cerr << "Error: No sample columns found.\n";
        return;
    }

    if (variantCount == 0) {
        if (!args.quiet) {
            std::cerr << "No biallelic variants found.\n";
        }
        for (int s = 0; s < numSamples; s++) {
            outBuf.append(sampleNames[s]);
            outBuf.append("\tNA\n");
        }
        return;
    }

    for (int s = 0; s < numSamples; s++) {
        outBuf.append(sampleNames[s]);
        outBuf.appendChar('\t');

        if (usedCount[s] == 0) {
            outBuf.append("NA\n");
            continue;
        }

        double e = sumExp[s];
        if (e <= 0.0) {
            outBuf.append("1.000000\n");
            continue;
        }

        double f = 1.0 - (obsHet[s] / e);
        outBuf.appendDouble(f);
        outBuf.appendChar('\n');
    }
}

// ============================================================================
// LEGACY CLASS INTERFACE (for compatibility)
// ============================================================================

void VCFXInbreedingCalculator::displayHelp() {
    ::displayHelp();
}

int VCFXInbreedingCalculator::parseGenotype(const std::string &s) {
    return parseGenotypeCode(s.data(), s.data() + s.size());
}

bool VCFXInbreedingCalculator::isBiallelic(const std::string &alt) {
    return alt.find(',') == std::string::npos;
}

FrequencyMode VCFXInbreedingCalculator::parseFreqMode(const std::string &modeStr) {
    if (modeStr == "global") return FrequencyMode::GLOBAL;
    return FrequencyMode::EXCLUDE_SAMPLE;
}

void VCFXInbreedingCalculator::calculateInbreeding(std::istream &in, std::ostream &out) {
    InbreedingArgs args;
    args.freqMode = freqMode_;
    args.skipBoundary = skipBoundary_;
    args.countBoundaryAsUsed = countBoundaryAsUsed_;

    OutputBuffer outBuf;
    calculateInbreedingStdin(in, args, outBuf);
}

int VCFXInbreedingCalculator::run(int argc, char *argv[]) {
    InbreedingArgs args = parseArgs(argc, argv);

    if (args.showHelp) {
        displayHelp();
        return 0;
    }

    // Store settings for legacy interface
    freqMode_ = args.freqMode;
    skipBoundary_ = args.skipBoundary;
    countBoundaryAsUsed_ = args.countBoundaryAsUsed;

    OutputBuffer outBuf;

    if (args.inputFile) {
        // Memory-mapped mode
        MappedFile file;
        if (!file.open(args.inputFile)) {
            std::cerr << "Error: Cannot open file: " << args.inputFile << "\n";
            return 1;
        }

        if (!args.quiet) {
            std::cerr << "Processing " << args.inputFile << " (" << file.size << " bytes)...\n";
        }

        calculateInbreedingMmap(file, args, outBuf);
    } else {
        // Stdin mode
        calculateInbreedingStdin(std::cin, args, outBuf);
    }

    return 0;
}

// ============================================================================
// MAIN ENTRY POINT
// ============================================================================

static void show_help() {
    ::displayHelp();
}

int main(int argc, char *argv[]) {
    vcfx::init_io();
    if (vcfx::handle_common_flags(argc, argv, "VCFX_inbreeding_calculator", show_help))
        return 0;
    VCFXInbreedingCalculator calc;
    return calc.run(argc, argv);
}
