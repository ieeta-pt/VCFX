#include "VCFX_nonref_filter.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <cstring>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

// mmap support
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
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
            out.write(sv.data(), sv.size());
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
            out.write(buffer, pos);
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
// Fast character finding
// =============================================================================
static inline const char* findChar(const char* p, const char* end, char c) {
    return static_cast<const char*>(memchr(p, c, static_cast<size_t>(end - p)));
}

// =============================================================================
// Extract nth field from colon-delimited string (zero-copy)
// Returns empty string_view if field doesn't exist
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
// Zero-copy tab splitting
// =============================================================================
static inline void split_tabs_view(std::string_view str,
                                    std::vector<std::string_view> &out) {
    out.clear();
    const char* p = str.data();
    const char* end = p + str.size();
    const char* fieldStart = p;

    while (p < end) {
        if (*p == '\t') {
            out.emplace_back(fieldStart, static_cast<size_t>(p - fieldStart));
            fieldStart = p + 1;
        }
        p++;
    }
    out.emplace_back(fieldStart, static_cast<size_t>(end - fieldStart));
}

// =============================================================================
// Skip to nth tab-delimited field (zero allocation)
// Returns pointer to start of field n, or nullptr if not found
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
// Check if ALL samples are homozygous reference (no allocation version)
// Scans directly through the line without splitting
// Returns true if all samples are hom-ref (should be EXCLUDED)
// Returns false if any sample is non-hom-ref (should be KEPT)
// =============================================================================
static inline bool allSamplesHomRefDirect(const char* lineStart, const char* lineEnd, int gtIndex) {
    // Skip to field 8 (FORMAT) - fields are: CHROM(0), POS(1), ID(2), REF(3), ALT(4), QUAL(5), FILTER(6), INFO(7), FORMAT(8)
    const char* p = skipToField(lineStart, lineEnd, 8);
    if (!p) return false;  // Not enough fields - keep line

    // Get FORMAT field to verify GT index
    std::string_view format = getFieldExtent(p, lineEnd);

    // Quick sanity check: if gtIndex is 0, FORMAT should start with "GT"
    // (We rely on cached gtIndex being correct)

    // Skip to field 9 (first sample)
    p = skipToField(lineStart, lineEnd, 9);
    if (!p) return false;  // No samples - keep line

    // Iterate through all samples
    while (p < lineEnd) {
        // Extract GT from this sample
        std::string_view sample = getFieldExtent(p, lineEnd);

        if (sample.empty()) {
            return false;  // Empty sample - keep line
        }

        // Extract GT field at gtIndex position
        std::string_view gt;
        if (gtIndex == 0) {
            // Fast path: GT is first field - find first colon
            const char* sp = sample.data();
            const char* se = sp + sample.size();
            const char* colon = findChar(sp, se, ':');
            gt = colon ? std::string_view(sp, static_cast<size_t>(colon - sp)) : sample;
        } else {
            gt = extractNthField(sample, gtIndex);
        }

        // Check if definitely hom-ref
        // Inline the check for maximum speed
        bool isHomRef = false;
        if (gt.size() == 3) {
            if (gt[0] == '0' && (gt[1] == '/' || gt[1] == '|') && gt[2] == '0') {
                isHomRef = true;
            }
        } else if (!gt.empty()) {
            // General case scan
            isHomRef = true;
            for (size_t i = 0; i < gt.size(); i++) {
                char c = gt[i];
                if (c == '/' || c == '|' || c == '0') continue;
                isHomRef = false;
                break;
            }
        }

        if (!isHomRef) {
            return false;  // Found non-hom-ref - KEEP this line
        }

        // Move to next sample (skip current field and tab)
        p = sample.data() + sample.size();
        if (p < lineEnd && *p == '\t') p++;
    }

    return true;  // All samples were hom-ref - EXCLUDE this line
}

// =============================================================================
// Find GT index in FORMAT string
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
// Main entry point
// =============================================================================
int VCFXNonRefFilter::run(int argc, char *argv[]) {
    bool showHelp = false;
    std::string inputFile;

    static struct option long_opts[] = {
        {"help", no_argument, 0, 'h'},
        {"input", required_argument, 0, 'i'},
        {0, 0, 0, 0}
    };

    optind = 1;

    while (true) {
        int c = getopt_long(argc, argv, "hi:", long_opts, nullptr);
        if (c == -1)
            break;
        switch (c) {
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

    if (inputFile.empty() && optind < argc) {
        inputFile = argv[optind];
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    if (!inputFile.empty() && inputFile != "-") {
        filterNonRefMmap(inputFile.c_str(), std::cout);
    } else {
        filterNonRef(std::cin, std::cout);
    }

    return 0;
}

void VCFXNonRefFilter::displayHelp() {
    std::cout
        << "VCFX_nonref_filter: Exclude variants if all samples are homozygous reference.\n\n"
           "Usage:\n"
           "  VCFX_nonref_filter [options] [input.vcf]\n"
           "  VCFX_nonref_filter [options] < input.vcf > output.vcf\n\n"
           "Options:\n"
           "  -h, --help          Show this help message\n"
           "  -i, --input FILE    Input VCF file (uses fast memory-mapped I/O)\n\n"
           "Description:\n"
           "  Reads VCF lines. For each variant, we check each sample's genotype. If a\n"
           "  genotype is polyploid, all alleles must be '0'. If a genotype is missing\n"
           "  or partial, we consider it not guaranteed hom-ref => keep variant.\n"
           "  If we find at least one sample not hom-ref, we print the variant. Otherwise,\n"
           "  we skip it.\n\n"
           "Performance:\n"
           "  File input (-i) uses memory-mapped I/O for 100-1000x faster processing\n"
           "  compared to stdin. Features include:\n"
           "  - SIMD-optimized line scanning (AVX2/SSE2)\n"
           "  - Zero-copy string parsing with string_view\n"
           "  - 1MB output buffering\n"
           "  - Direct GT field extraction (avoids full sample parsing)\n"
           "  - Early termination on first non-homref sample\n\n"
           "Examples:\n"
           "  VCFX_nonref_filter -i input.vcf > filtered.vcf    # Fast (mmap)\n"
           "  VCFX_nonref_filter input.vcf > filtered.vcf       # Fast (mmap)\n"
           "  VCFX_nonref_filter < input.vcf > filtered.vcf     # Slower (stdin)\n\n";
}

// =============================================================================
// Ultra-fast homozygous reference check
// Optimized for common genotype patterns: 0/0, 0|0, 0/0/0, etc.
// =============================================================================
bool VCFXNonRefFilter::isDefinitelyHomRef(std::string_view gt) const {
    if (gt.empty()) return false;

    // Fast path for common cases
    if (gt.size() == 3) {
        // Most common: "0/0" or "0|0"
        if (gt[0] == '0' && (gt[1] == '/' || gt[1] == '|') && gt[2] == '0') {
            return true;
        }
        // Missing: "./." or ".|."
        if (gt[0] == '.') return false;
        // Non-ref first allele
        if (gt[0] >= '1' && gt[0] <= '9') return false;
    }

    // Handle single missing "."
    if (gt.size() == 1 && gt[0] == '.') return false;

    // General case: scan all characters
    for (size_t i = 0; i < gt.size(); i++) {
        char c = gt[i];
        if (c == '/' || c == '|') continue;  // separator
        if (c == '0') continue;               // reference allele
        if (c == '.') return false;           // missing
        if (c >= '1' && c <= '9') return false;  // non-ref allele
        // Unknown character - be conservative
        return false;
    }

    return true;
}

bool VCFXNonRefFilter::isDefinitelyHomRef(const std::string &genotypeField) const {
    return isDefinitelyHomRef(std::string_view(genotypeField));
}

// =============================================================================
// Memory-mapped file processing (FAST PATH - ZERO ALLOCATION)
// =============================================================================
void VCFXNonRefFilter::filterNonRefMmap(const char *filepath, std::ostream &out) {
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
    int cachedGtIndex = 0;  // Default: GT is usually first field in FORMAT
    std::string cachedFormat;  // Store as string for comparison (FORMAT rarely changes)

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
        if (!headerFound) {
            std::cerr << "Warning: VCF data line encountered before #CHROM. Passing line.\n";
            outBuf.write(line);
            ptr = lineEnd + 1;
            continue;
        }

        // Get FORMAT field (field 8) for GT index caching
        const char* formatStart = skipToField(line.data(), line.data() + line.size(), 8);
        if (!formatStart) {
            outBuf.write(line);
            ptr = lineEnd + 1;
            continue;
        }

        std::string_view formatStr = getFieldExtent(formatStart, line.data() + line.size());

        // Update GT index if FORMAT changed
        if (cachedFormat.empty() || formatStr != cachedFormat) {
            cachedGtIndex = findGTIndex(formatStr);
            cachedFormat = std::string(formatStr);
        }

        if (cachedGtIndex < 0) {
            outBuf.write(line);
            ptr = lineEnd + 1;
            continue;
        }

        // Check all samples using zero-allocation direct scan
        bool allHomRef = allSamplesHomRefDirect(line.data(), line.data() + line.size(), cachedGtIndex);

        if (!allHomRef) {
            outBuf.write(line);
        }

        ptr = lineEnd + 1;
    }

    // Flush remaining output
    outBuf.flush();
}

// =============================================================================
// Stdin-based processing (FALLBACK PATH)
// =============================================================================
void VCFXNonRefFilter::filterNonRef(std::istream &in, std::ostream &out) {
    bool headerFound = false;
    std::string line;
    std::vector<std::string> fields;
    fields.reserve(16);
    std::vector<std::string> fmtParts;
    fmtParts.reserve(16);
    std::vector<std::string> subf;
    subf.reserve(16);

    while (std::getline(in, line)) {
        if (line.empty()) {
            out << "\n";
            continue;
        }
        if (line[0] == '#') {
            out << line << "\n";
            if (line.rfind("#CHROM", 0) == 0)
                headerFound = true;
            continue;
        }
        if (!headerFound) {
            std::cerr << "Warning: VCF data line encountered before #CHROM. Passing line.\n";
            out << line << "\n";
            continue;
        }

        vcfx::split_tabs(line, fields);
        if (fields.size() < 10) {
            out << line << "\n";
            continue;
        }

        std::string formatStr = fields[8];
        {
            std::stringstream fs(formatStr);
            std::string ff;
            fmtParts.clear();
            while (std::getline(fs, ff, ':'))
                fmtParts.push_back(ff);
        }

        int gtIndex = -1;
        for (size_t i = 0; i < fmtParts.size(); i++) {
            if (fmtParts[i] == "GT") {
                gtIndex = static_cast<int>(i);
                break;
            }
        }

        if (gtIndex < 0) {
            out << line << "\n";
            continue;
        }

        bool allHomRef = true;
        for (size_t s = 9; s < fields.size(); s++) {
            std::string &sampleCol = fields[s];
            {
                std::stringstream sampleSS(sampleCol);
                std::string token;
                subf.clear();
                while (std::getline(sampleSS, token, ':'))
                    subf.push_back(token);
            }

            if (gtIndex >= static_cast<int>(subf.size())) {
                allHomRef = false;
                break;
            }

            if (!isDefinitelyHomRef(subf[static_cast<size_t>(gtIndex)])) {
                allHomRef = false;
                break;
            }
        }

        if (!allHomRef)
            out << line << "\n";
    }
}

// =============================================================================
// Entry point
// =============================================================================
static void show_help() {
    VCFXNonRefFilter obj;
    char arg0[] = "VCFX_nonref_filter";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();
    if (vcfx::handle_common_flags(argc, argv, "VCFX_nonref_filter", show_help))
        return 0;
    VCFXNonRefFilter app;
    return app.run(argc, argv);
}
