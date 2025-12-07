#include "VCFX_phase_checker.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

// mmap support
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/select.h>
#include <unistd.h>

// SIMD support detection
#if defined(__x86_64__) || defined(_M_X64)
#if defined(__AVX2__)
#include <immintrin.h>
#define USE_AVX2 1
#elif defined(__SSE2__)
#include <emmintrin.h>
#define USE_SSE2 1
#endif
#endif

// =============================================================================
// Memory-mapped file wrapper (RAII)
// =============================================================================
struct MappedFile {
    const char *data = nullptr;
    size_t size = 0;
    int fd = -1;

    bool open(const char *path) {
        fd = ::open(path, O_RDONLY);
        if (fd < 0)
            return false;

        struct stat st;
        if (fstat(fd, &st) < 0) {
            close();
            return false;
        }
        size = static_cast<size_t>(st.st_size);
        if (size == 0) {
            return true;
        }

        data = static_cast<const char *>(mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0));
        if (data == MAP_FAILED) {
            data = nullptr;
            close();
            return false;
        }

        // Advise kernel for sequential access and trigger read-ahead
        madvise(const_cast<char *>(data), size, MADV_SEQUENTIAL);
        madvise(const_cast<char *>(data), size, MADV_WILLNEED);
        return true;
    }

    void close() {
        if (data && size > 0) {
            munmap(const_cast<char *>(data), size);
            data = nullptr;
        }
        if (fd >= 0) {
            ::close(fd);
            fd = -1;
        }
        size = 0;
    }

    ~MappedFile() { close(); }
};

// =============================================================================
// Output buffer for efficient writing
// =============================================================================
class OutputBuffer {
    static constexpr size_t BUFFER_SIZE = 1024 * 1024;  // 1MB buffer
    char* buffer;
    size_t pos = 0;
    std::ostream& out;

public:
    explicit OutputBuffer(std::ostream& os) : out(os) {
        buffer = new char[BUFFER_SIZE];
    }

    ~OutputBuffer() {
        flush();
        delete[] buffer;
    }

    void write(std::string_view sv) {
        if (pos + sv.size() + 1 > BUFFER_SIZE) {
            flush();
        }
        if (sv.size() + 1 > BUFFER_SIZE) {
            // Line larger than buffer - write directly
            out.write(sv.data(), static_cast<std::streamsize>(sv.size()));
            out.put('\n');
            return;
        }
        memcpy(buffer + pos, sv.data(), sv.size());
        pos += sv.size();
        buffer[pos++] = '\n';
    }

    void writeNewline() {
        if (pos + 1 > BUFFER_SIZE) {
            flush();
        }
        buffer[pos++] = '\n';
    }

    void flush() {
        if (pos > 0) {
            out.write(buffer, static_cast<std::streamsize>(pos));
            pos = 0;
        }
    }
};

// =============================================================================
// SIMD-optimized newline detection
// =============================================================================
#if defined(USE_AVX2)
static inline const char *findNewlineSIMD(const char *p, const char *end) {
    const __m256i nl = _mm256_set1_epi8('\n');
    while (p + 32 <= end) {
        __m256i chunk = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(p));
        __m256i cmp = _mm256_cmpeq_epi8(chunk, nl);
        int mask = _mm256_movemask_epi8(cmp);
        if (mask) {
            return p + __builtin_ctz(static_cast<unsigned int>(mask));
        }
        p += 32;
    }
    return static_cast<const char *>(memchr(p, '\n', static_cast<size_t>(end - p)));
}
#elif defined(USE_SSE2)
static inline const char *findNewlineSIMD(const char *p, const char *end) {
    const __m128i nl = _mm_set1_epi8('\n');
    while (p + 16 <= end) {
        __m128i chunk = _mm_loadu_si128(reinterpret_cast<const __m128i *>(p));
        __m128i cmp = _mm_cmpeq_epi8(chunk, nl);
        int mask = _mm_movemask_epi8(cmp);
        if (mask) {
            return p + __builtin_ctz(static_cast<unsigned int>(mask));
        }
        p += 16;
    }
    return static_cast<const char *>(memchr(p, '\n', static_cast<size_t>(end - p)));
}
#else
static inline const char *findNewlineSIMD(const char *p, const char *end) {
    return static_cast<const char *>(memchr(p, '\n', static_cast<size_t>(end - p)));
}
#endif

// =============================================================================
// Find GT index in FORMAT string (zero-allocation)
// =============================================================================
static inline int findGTIndex(std::string_view format) {
    const char* p = format.data();
    const char* end = p + format.size();
    int idx = 0;
    const char* fieldStart = p;

    while (p <= end) {
        if (p == end || *p == ':') {
            size_t len = static_cast<size_t>(p - fieldStart);
            if (len == 2 && fieldStart[0] == 'G' && fieldStart[1] == 'T') {
                return idx;
            }
            idx++;
            fieldStart = p + 1;
        }
        p++;
    }
    return -1;
}

// =============================================================================
// Extract nth colon-delimited field (zero-copy)
// =============================================================================
static inline std::string_view extractNthField(std::string_view str, int n) {
    if (n < 0) return {};

    const char* p = str.data();
    const char* end = p + str.size();
    int fieldIdx = 0;
    const char* fieldStart = p;

    while (p <= end) {
        if (p == end || *p == ':') {
            if (fieldIdx == n) {
                return std::string_view(fieldStart, static_cast<size_t>(p - fieldStart));
            }
            fieldIdx++;
            fieldStart = p + 1;
        }
        p++;
    }
    return {};
}

// =============================================================================
// Optimized phasing check (zero-allocation)
// =============================================================================
static inline bool isFullyPhasedFast(std::string_view gt) {
    if (__builtin_expect(gt.empty(), 0)) return false;

    const char* data = gt.data();
    const size_t len = gt.size();

    // Fast path: 3-char genotypes like "0|1", "1|0", "0|0", "1|1" (most common case)
    if (__builtin_expect(len == 3, 1)) {
        if (data[1] == '|' && data[0] != '.' && data[2] != '.') {
            return true;
        }
        // Common missing patterns
        if (data[0] == '.') return false;  // "./.", ".|.", etc.
        if (data[1] != '|') return false;  // No pipe = not phased
        return false;
    }

    // Handle single missing "."
    if (len == 1) return false;  // Single char = haploid or missing, never phased

    // Handle common missing patterns for length > 3
    if (data[0] == '.' && len >= 3 && (data[1] == '/' || data[1] == '|') && data[2] == '.') {
        return false;
    }

    // Single-pass scan: check for '|' presence, '/' absence, and missing alleles
    // This is more efficient than multiple memchr calls
    bool hasPipe = false;
    const char* p = data;
    const char* end = data + len;
    const char* alleleStart = p;

    while (p < end) {
        const char c = *p;
        if (c == '|') {
            // Check allele before pipe is not empty or just "."
            const size_t alleleLen = static_cast<size_t>(p - alleleStart);
            if (alleleLen == 0 || (alleleLen == 1 && *alleleStart == '.')) return false;
            hasPipe = true;
            alleleStart = p + 1;
        } else if (c == '/') {
            return false;  // Unphased separator found
        }
        p++;
    }

    // Check final allele
    const size_t alleleLen = static_cast<size_t>(p - alleleStart);
    if (alleleLen == 0 || (alleleLen == 1 && *alleleStart == '.')) return false;

    return hasPipe;  // Must have at least one '|'
}

// =============================================================================
// Skip to nth tab-delimited field (zero allocation)
// =============================================================================
static inline const char* skipToField(const char* p, const char* end, int n) {
    int fieldIdx = 0;
    while (p < end && fieldIdx < n) {
        if (*p == '\t') fieldIdx++;
        p++;
    }
    return (fieldIdx == n) ? p : nullptr;
}

// =============================================================================
// Get field extent (from current position to next tab or end)
// =============================================================================
static inline std::string_view getFieldExtent(const char* p, const char* end) {
    const char* start = p;
    while (p < end && *p != '\t') p++;
    return std::string_view(start, static_cast<size_t>(p - start));
}

// =============================================================================
// Check if all samples are phased, operating directly on line memory
// Returns: 1 = all phased, 0 = not all phased, -1 = invalid line (< 10 fields)
// =============================================================================
static inline int checkAllSamplesPhasedDirect(std::string_view line,
                                               std::string& cachedFormat,
                                               int& cachedGtIndex) {
    const char* p = line.data();
    const char* end = p + line.size();

    // Skip to FORMAT field (field 8) - count tabs inline for efficiency
    int tabCount = 0;
    while (p < end && tabCount < 8) {
        if (*p == '\t') tabCount++;
        p++;
    }
    if (__builtin_expect(tabCount < 8 || p >= end, 0)) return -1;  // Invalid line

    // Get FORMAT field extent
    const char* formatStart = p;
    while (p < end && *p != '\t') p++;
    std::string_view formatStr(formatStart, static_cast<size_t>(p - formatStart));

    // Update GT index if FORMAT changed (compare length first for quick rejection)
    if (formatStr.size() != cachedFormat.size() || formatStr != cachedFormat) {
        cachedGtIndex = findGTIndex(formatStr);
        cachedFormat.assign(formatStr.data(), formatStr.size());
    }

    if (__builtin_expect(cachedGtIndex < 0, 0)) return 0;  // No GT field

    // Skip tab after FORMAT to reach first sample (field 9)
    if (p >= end || *p != '\t') return -1;  // No samples
    p++;  // Now at start of first sample

    // Check each sample
    while (p < end) {
        // Find sample extent
        const char* sampleStart = p;
        while (p < end && *p != '\t') p++;
        const size_t sampleLen = static_cast<size_t>(p - sampleStart);

        if (__builtin_expect(sampleLen == 0, 0)) {
            return 0;  // Empty sample
        }

        // Extract GT field
        std::string_view gt;
        if (__builtin_expect(cachedGtIndex == 0, 1)) {
            // Fast path: GT is first field - find first colon
            const char* colon = static_cast<const char*>(
                memchr(sampleStart, ':', sampleLen));
            gt = colon ? std::string_view(sampleStart, static_cast<size_t>(colon - sampleStart))
                       : std::string_view(sampleStart, sampleLen);
        } else {
            gt = extractNthField(std::string_view(sampleStart, sampleLen), cachedGtIndex);
        }

        if (!isFullyPhasedFast(gt)) {
            return 0;
        }

        // Move to next sample
        if (p < end) p++;  // Skip tab
    }

    return 1;
}

// =============================================================================
// VCFXPhaseChecker implementation
// =============================================================================

int VCFXPhaseChecker::run(int argc, char *argv[]) {
    bool showHelp = false;
    std::string inputFile;
    quiet_ = false;

    // Check if stdin has data available (for help display when no input)
    struct timeval tv;
    fd_set fds;
    tv.tv_sec = 0;
    tv.tv_usec = 0;
    FD_ZERO(&fds);
    FD_SET(STDIN_FILENO, &fds);
    bool hasStdinInput = select(STDIN_FILENO + 1, &fds, NULL, NULL, &tv) > 0;

    // If no arguments and no stdin input, display help
    if (argc == 1 && !hasStdinInput) {
        displayHelp();
        return 0;
    }

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
            showHelp = true;
            break;
        case 'i':
            inputFile = optarg;
            break;
        case 'q':
            quiet_ = true;
            break;
        default:
            showHelp = true;
        }
    }

    // Support positional argument for file
    if (inputFile.empty() && optind < argc) {
        inputFile = argv[optind];
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    // Use mmap fast path for file input
    if (!inputFile.empty() && inputFile != "-") {
        filterPhaseCheckedMmap(inputFile.c_str(), std::cout);
    } else {
        processVCF(std::cin, std::cout);
    }

    return 0;
}

void VCFXPhaseChecker::displayHelp() {
    std::cout <<
"VCFX_phase_checker: Output only VCF variant lines in which every sample genotype is fully phased.\n\n"
"Usage:\n"
"  VCFX_phase_checker [options] [input.vcf]\n"
"  VCFX_phase_checker [options] < input.vcf > phased_output.vcf\n\n"
"Options:\n"
"  -h, --help          Display this help message and exit\n"
"  -i, --input FILE    Input VCF file (uses fast memory-mapped I/O)\n"
"  -q, --quiet         Suppress warning messages to stderr\n\n"
"Description:\n"
"  The tool reads a VCF and checks the GT field (genotype) for each sample.\n"
"  A genotype is considered fully phased if it uses the '|' separator (e.g., \"0|1\")\n"
"  and contains no missing alleles. If every sample in a variant line is fully phased,\n"
"  the line is printed to stdout; otherwise, it is skipped with a warning to stderr.\n\n"
"Performance:\n"
"  File input (-i) uses memory-mapped I/O for 10-12x faster processing compared to stdin.\n"
"  Features include:\n"
"  - SIMD-optimized line scanning (AVX2/SSE2 on x86_64)\n"
"  - Zero-copy string parsing with string_view\n"
"  - 1MB output buffering\n"
"  - FORMAT field caching (GT index computed once per unique FORMAT)\n"
"  - Early termination on first unphased sample\n\n"
"Examples:\n"
"  VCFX_phase_checker -i input.vcf > phased.vcf       # Fast (mmap)\n"
"  VCFX_phase_checker input.vcf > phased.vcf          # Fast (mmap)\n"
"  VCFX_phase_checker < input.vcf > phased.vcf        # Slower (stdin)\n"
"  VCFX_phase_checker -q -i input.vcf > phased.vcf    # Quiet mode (no warnings)\n";
}

bool VCFXPhaseChecker::isFullyPhased(const std::string &gt) const {
    // Delegate to the fast version for compatibility
    return isFullyPhasedFast(std::string_view(gt));
}

// =============================================================================
// Memory-mapped file processing (FAST PATH)
// =============================================================================
void VCFXPhaseChecker::filterPhaseCheckedMmap(const char *filepath, std::ostream &out) {
    MappedFile file;
    if (!file.open(filepath)) {
        std::cerr << "Error: Cannot open file: " << filepath << "\n";
        return;
    }

    if (file.size == 0) return;

    const char *ptr = file.data;
    const char *end = file.data + file.size;

    // Use buffered output for efficiency
    OutputBuffer outBuf(out);

    bool headerFound = false;
    int cachedGtIndex = 0;
    std::string cachedFormat;
    cachedFormat.reserve(64);  // Pre-allocate to avoid reallocation

    while (ptr < end) {
        const char *lineEnd = findNewlineSIMD(ptr, end);
        if (!lineEnd) lineEnd = end;

        std::string_view line(ptr, static_cast<size_t>(lineEnd - ptr));

        // Handle Windows line endings
        if (!line.empty() && line.back() == '\r') {
            line.remove_suffix(1);
        }

        // Empty line
        if (line.empty()) {
            outBuf.writeNewline();
            ptr = lineEnd + 1;
            continue;
        }

        // Header lines - pass through
        if (line[0] == '#') {
            outBuf.write(line);
            if (line.size() >= 6 && line[1] == 'C' && line[2] == 'H' &&
                line[3] == 'R' && line[4] == 'O' && line[5] == 'M') {
                headerFound = true;
            }
            ptr = lineEnd + 1;
            continue;
        }

        // Data line before header
        if (__builtin_expect(!headerFound, 0)) {
            if (!quiet_) {
                std::cerr << "Warning: Data line encountered before #CHROM header; skipping line.\n";
            }
            ptr = lineEnd + 1;
            continue;
        }

        // Check all samples using zero-allocation direct scan
        // Returns: 1 = all phased, 0 = not phased, -1 = invalid line
        const int result = checkAllSamplesPhasedDirect(line, cachedFormat, cachedGtIndex);

        if (result == 1) {
            outBuf.write(line);
        } else if (result == -1) {
            if (!quiet_) {
                std::cerr << "Warning: Invalid VCF line with fewer than 10 columns; skipping line.\n";
            }
        } else if (!quiet_) {
            // Extract CHROM and POS for warning using memchr (faster than find)
            const char* tab1 = static_cast<const char*>(memchr(line.data(), '\t', line.size()));
            if (tab1) {
                const size_t remaining = line.size() - static_cast<size_t>(tab1 + 1 - line.data());
                const char* tab2 = static_cast<const char*>(memchr(tab1 + 1, '\t', remaining));
                if (tab2) {
                    std::cerr << "Unphased genotype found at CHROM="
                              << std::string_view(line.data(), static_cast<size_t>(tab1 - line.data()))
                              << ", POS="
                              << std::string_view(tab1 + 1, static_cast<size_t>(tab2 - tab1 - 1))
                              << "; line skipped.\n";
                }
            }
        }

        ptr = lineEnd + 1;
    }

    outBuf.flush();
}

// =============================================================================
// Stdin-based processing (FALLBACK PATH - still optimized)
// =============================================================================
void VCFXPhaseChecker::processVCF(std::istream &in, std::ostream &out) {
    bool headerFound = false;
    std::string line;

    // Reuse vectors across lines (preserve capacity)
    std::vector<std::string_view> fields;
    fields.reserve(16);

    // Cache GT index between lines
    std::string cachedFormat;
    cachedFormat.reserve(64);  // Pre-allocate to avoid reallocation
    int cachedGtIndex = -1;

    while (std::getline(in, line)) {
        if (__builtin_expect(line.empty(), 0)) {
            out << '\n';
            continue;
        }

        if (line[0] == '#') {
            out << line << '\n';
            if (line.size() >= 6 && line[1] == 'C' && line[2] == 'H' &&
                line[3] == 'R' && line[4] == 'O' && line[5] == 'M') {
                headerFound = true;
            }
            continue;
        }

        if (__builtin_expect(!headerFound, 0)) {
            if (!quiet_) {
                std::cerr << "Warning: Data line encountered before #CHROM header; skipping line.\n";
            }
            continue;
        }

        // Zero-copy tab splitting
        vcfx::split_tabs_view(line, fields);

        if (__builtin_expect(fields.size() < 10, 0)) {
            if (!quiet_) {
                std::cerr << "Warning: Invalid VCF line with fewer than 10 columns; skipping line.\n";
            }
            continue;
        }

        // Cache GT index if FORMAT changed (compare length first for quick rejection)
        std::string_view formatStr = fields[8];
        if (formatStr.size() != cachedFormat.size() || formatStr != cachedFormat) {
            cachedFormat.assign(formatStr.data(), formatStr.size());
            cachedGtIndex = findGTIndex(formatStr);
        }

        if (__builtin_expect(cachedGtIndex < 0, 0)) {
            if (!quiet_) {
                std::cerr << "Warning: GT field not found; skipping line.\n";
            }
            continue;
        }

        // Check all samples with optimized isFullyPhasedFast
        bool allPhased = true;
        for (size_t i = 9; i < fields.size(); ++i) {
            std::string_view gt;
            std::string_view sample = fields[i];
            if (__builtin_expect(cachedGtIndex == 0, 1)) {
                // Fast path: GT is first field - use memchr instead of find
                const char* colon = static_cast<const char*>(
                    memchr(sample.data(), ':', sample.size()));
                gt = colon ? std::string_view(sample.data(), static_cast<size_t>(colon - sample.data()))
                           : sample;
            } else {
                gt = extractNthField(sample, cachedGtIndex);
            }

            if (!isFullyPhasedFast(gt)) {
                allPhased = false;
                break;
            }
        }

        if (allPhased) {
            out << line << '\n';
        } else if (!quiet_) {
            std::cerr << "Unphased genotype found at CHROM=" << fields[0]
                      << ", POS=" << fields[1] << "; line skipped.\n";
        }
    }
}

// =============================================================================
// Entry point
// =============================================================================
static void show_help() {
    VCFXPhaseChecker obj;
    char arg0[] = "VCFX_phase_checker";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_phase_checker", show_help))
        return 0;
    VCFXPhaseChecker checker;
    return checker.run(argc, argv);
}
