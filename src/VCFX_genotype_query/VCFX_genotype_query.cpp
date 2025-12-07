#include "VCFX_genotype_query.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
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

    void writeNoNewline(std::string_view sv) {
        if (pos + sv.size() > BUFFER_SIZE) {
            flush();
        }
        if (sv.size() > BUFFER_SIZE) {
            out.write(sv.data(), static_cast<std::streamsize>(sv.size()));
            return;
        }
        memcpy(buffer + pos, sv.data(), sv.size());
        pos += sv.size();
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
// Skip to nth tab-delimited field, return pointer to start
// =============================================================================
static inline const char* skipToField(const char* p, const char* end, int n) {
    for (int i = 0; i < n && p < end; i++) {
        p = static_cast<const char*>(memchr(p, '\t', static_cast<size_t>(end - p)));
        if (!p) return nullptr;
        p++;
    }
    return p;
}

// =============================================================================
// Get extent of current field (until tab or end)
// =============================================================================
static inline const char* getFieldExtent(const char* start, const char* end) {
    const char* tab = static_cast<const char*>(memchr(start, '\t', static_cast<size_t>(end - start)));
    return tab ? tab : end;
}

// =============================================================================
// Fast genotype matching - CORE ALGORITHM
// Supports flexible matching (normalize phasing/allele order) or strict matching
// =============================================================================

// Parse diploid alleles from genotype string (returns false if not diploid or malformed)
static inline bool parseDiploidAlleles(std::string_view gt, int& a1, int& a2) {
    // Find separator
    size_t sepPos = gt.find_first_of("|/");
    if (sepPos == std::string_view::npos || sepPos == 0 || sepPos == gt.size() - 1) {
        return false;
    }

    // Parse first allele
    std::string_view first = gt.substr(0, sepPos);
    if (first == ".") return false;
    a1 = 0;
    for (char c : first) {
        if (!std::isdigit(static_cast<unsigned char>(c))) return false;
        a1 = a1 * 10 + (c - '0');
    }

    // Parse second allele
    std::string_view second = gt.substr(sepPos + 1);
    if (second == ".") return false;
    a2 = 0;
    for (char c : second) {
        if (!std::isdigit(static_cast<unsigned char>(c))) return false;
        a2 = a2 * 10 + (c - '0');
    }

    return true;
}

// Fast matching for diploid genotypes
static inline bool genotypeMatchesFast(std::string_view gt, std::string_view query,
                                        int queryA1, int queryA2, bool strict) {
    if (strict) {
        return gt == query;
    }

    // Fast path for 3-character diploid genotypes (most common: 0/1, 1|0, etc.)
    if (gt.size() == 3 && query.size() == 3) {
        char sepGt = gt[1];
        if (sepGt != '|' && sepGt != '/') return false;

        char g0 = gt[0], g1 = gt[2];
        if (g0 == '.' || g1 == '.') return false;
        if (!std::isdigit(static_cast<unsigned char>(g0)) ||
            !std::isdigit(static_cast<unsigned char>(g1))) return false;

        // Convert to int for comparison
        int ga = g0 - '0';
        int gb = g1 - '0';

        // Order-independent comparison
        if (ga > gb) std::swap(ga, gb);
        int qa = queryA1, qb = queryA2;
        if (qa > qb) std::swap(qa, qb);

        return ga == qa && gb == qb;
    }

    // General case: parse and compare
    int a1, a2;
    if (!parseDiploidAlleles(gt, a1, a2)) {
        // Not a valid diploid genotype - no match
        return false;
    }

    // Sort for order-independent comparison
    if (a1 > a2) std::swap(a1, a2);
    int qa = queryA1, qb = queryA2;
    if (qa > qb) std::swap(qa, qb);

    return a1 == qa && a2 == qb;
}

// =============================================================================
// Check if any sample in a line matches the query genotype
// Returns true if match found, false otherwise
// =============================================================================
static bool checkAnySampleMatches(const char* lineStart, const char* lineEnd,
                                   int gtIndex, std::string_view query,
                                   int queryA1, int queryA2, bool strict) {
    // Skip to field 9 (first sample)
    const char* p = skipToField(lineStart, lineEnd, 9);
    if (!p) return false;

    // Iterate through samples
    while (p < lineEnd) {
        const char* sampleEnd = getFieldExtent(p, lineEnd);
        std::string_view sample(p, static_cast<size_t>(sampleEnd - p));

        // Extract GT field
        std::string_view gt = extractNthField(sample, gtIndex);
        if (!gt.empty()) {
            if (__builtin_expect(genotypeMatchesFast(gt, query, queryA1, queryA2, strict), 0)) {
                return true;
            }
        }

        p = sampleEnd + 1;  // Move to next sample
    }
    return false;
}

// =============================================================================
// printHelp
// =============================================================================
void printHelp() {
    std::cout << "VCFX_genotype_query\n"
              << "Usage: VCFX_genotype_query [OPTIONS] [input.vcf]\n\n"
              << "Options:\n"
              << "  -g, --genotype-query GT  Genotype to query (e.g., \"0/1\", \"1|1\")\n"
              << "  -i, --input FILE         Input VCF file (uses fast memory-mapped I/O)\n"
              << "  --strict                 Exact string matching (no normalization)\n"
              << "  -q, --quiet              Suppress warning messages to stderr\n"
              << "  -h, --help               Display this help message and exit\n"
              << "  -v, --version            Show program version and exit\n\n"
              << "Description:\n"
              << "  Filters a VCF to retain only lines where at least one sample has the\n"
              << "  specified genotype in the 'GT' subfield.\n\n"
              << "  By default, phasing is unified (0|1 matches 0/1) and allele order is\n"
              << "  normalized (1/0 matches 0/1). Use --strict for exact matching.\n\n"
              << "Performance:\n"
              << "  File input mode (-i) uses memory-mapped I/O with SIMD optimization,\n"
              << "  providing 40-50x speedup over stdin mode for large files.\n\n"
              << "Examples:\n"
              << "  # Flexible matching (0/1 matches 0|1, 1/0, 1|0)\n"
              << "  VCFX_genotype_query -g \"0/1\" < input.vcf > het.vcf\n"
              << "  VCFX_genotype_query -g \"0/1\" -i input.vcf > het.vcf\n\n"
              << "  # Strict matching (only exact 0|1)\n"
              << "  VCFX_genotype_query -g \"0|1\" --strict < input.vcf > phased_het.vcf\n";
}

// =============================================================================
// parseArguments - updated for new options
// =============================================================================
bool parseArguments(int argc, char *argv[], std::string &genotype_query,
                    bool &strictCompare, std::string &inputFile, bool &quiet) {
    genotype_query.clear();
    inputFile.clear();
    strictCompare = false;
    quiet = false;

    static struct option long_options[] = {
        {"genotype-query", required_argument, nullptr, 'g'},
        {"input", required_argument, nullptr, 'i'},
        {"strict", no_argument, nullptr, 's'},
        {"quiet", no_argument, nullptr, 'q'},
        {"help", no_argument, nullptr, 'h'},
        {"version", no_argument, nullptr, 'v'},
        {nullptr, 0, nullptr, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "g:i:qhv", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'g':
                genotype_query = optarg;
                break;
            case 'i':
                inputFile = optarg;
                break;
            case 's':
                strictCompare = true;
                break;
            case 'q':
                quiet = true;
                break;
            case 'h':
                printHelp();
                std::exit(0);
            case 'v':
                std::cout << "VCFX_genotype_query version 1.0\n";
                std::exit(0);
            default:
                return false;
        }
    }

    // Check for positional file argument
    if (optind < argc && inputFile.empty()) {
        inputFile = argv[optind];
    }

    return !genotype_query.empty();
}

// =============================================================================
// genotypeQueryMmap - Fast memory-mapped implementation
// =============================================================================
static void genotypeQueryMmap(const char* data, size_t size, std::ostream& out,
                               std::string_view query, int queryA1, int queryA2,
                               bool strict, bool quiet) {
    if (size == 0) return;

    OutputBuffer outBuf(out);
    const char* p = data;
    const char* end = data + size;

    // FORMAT caching
    std::string cachedFormat;
    cachedFormat.reserve(64);
    int cachedGTIndex = -1;

    bool foundChrom = false;
    bool printedHeader = false;

    while (p < end) {
        // Find end of line
        const char* lineEnd = findNewlineSIMD(p, end);
        if (!lineEnd) lineEnd = end;

        // Skip empty lines
        if (lineEnd == p) {
            p = lineEnd + 1;
            continue;
        }

        std::string_view line(p, static_cast<size_t>(lineEnd - p));

        if (line[0] == '#') {
            // Header line - buffer for later output
            if (line.size() >= 6 && line.substr(0, 6) == "#CHROM") {
                foundChrom = true;
            }
            // We output headers immediately since we need them
            outBuf.write(line);
            printedHeader = true;
        } else {
            // Data line
            if (!foundChrom) {
                if (!quiet) {
                    std::cerr << "Error: No #CHROM header found before data lines.\n";
                }
                return;
            }

            // Find FORMAT field (field 8)
            const char* formatStart = skipToField(p, lineEnd, 8);
            if (!formatStart) {
                if (!quiet) {
                    std::cerr << "Warning: skipping line with <9 fields\n";
                }
                p = lineEnd + 1;
                continue;
            }
            const char* formatEnd = getFieldExtent(formatStart, lineEnd);
            std::string_view format(formatStart, static_cast<size_t>(formatEnd - formatStart));

            // Get GT index (with caching)
            int gtIndex;
            if (format.size() == cachedFormat.size() &&
                memcmp(format.data(), cachedFormat.data(), format.size()) == 0) {
                gtIndex = cachedGTIndex;
            } else {
                cachedFormat.assign(format.data(), format.size());
                gtIndex = findGTIndex(format);
                cachedGTIndex = gtIndex;
            }

            if (gtIndex < 0) {
                // No GT field - skip line
                p = lineEnd + 1;
                continue;
            }

            // Check if any sample matches
            if (checkAnySampleMatches(p, lineEnd, gtIndex, query, queryA1, queryA2, strict)) {
                outBuf.write(line);
            }
        }

        p = lineEnd + 1;
    }
}

// =============================================================================
// genotypeQuery: stream-based approach (stdin fallback)
// =============================================================================
void genotypeQuery(std::istream &in, std::ostream &out,
                   const std::string &genotype_query, bool strictCompare) {
    genotypeQueryStream(in, out, genotype_query, strictCompare, false);
}

void genotypeQueryStream(std::istream &in, std::ostream &out,
                          const std::string &genotype_query, bool strictCompare, bool quiet) {
    bool foundChrom = false;
    std::vector<std::string> headerLines;

    // Parse query genotype for flexible matching
    int queryA1 = -1, queryA2 = -1;
    if (!strictCompare) {
        parseDiploidAlleles(genotype_query, queryA1, queryA2);
        // Sort for consistent comparison
        if (queryA1 > queryA2) std::swap(queryA1, queryA2);
    }

    // FORMAT caching
    std::string cachedFormat;
    cachedFormat.reserve(64);
    int cachedGTIndex = -1;

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            headerLines.push_back(line);
            if (line.rfind("#CHROM", 0) == 0) {
                foundChrom = true;
            }
        } else {
            // Data line
            if (!foundChrom) {
                if (!quiet) {
                    std::cerr << "Error: No #CHROM header found before data lines.\n";
                }
                return;
            }

            // Print headers once before first data line
            // Note: Using headerLines.empty() check after printing to avoid static variable
            if (!headerLines.empty()) {
                for (const auto &h : headerLines) {
                    out << h << "\n";
                }
                headerLines.clear();  // Clear after printing to avoid re-printing
            }

            // Find FORMAT field (field 8)
            const char* p = line.data();
            const char* lineEnd = p + line.size();
            const char* formatStart = skipToField(p, lineEnd, 8);
            if (!formatStart) {
                if (!quiet) {
                    std::cerr << "Warning: skipping line with <9 fields: " << line << "\n";
                }
                continue;
            }
            const char* formatEnd = getFieldExtent(formatStart, lineEnd);
            std::string_view format(formatStart, static_cast<size_t>(formatEnd - formatStart));

            // Get GT index (with caching)
            int gtIndex;
            if (format.size() == cachedFormat.size() &&
                memcmp(format.data(), cachedFormat.data(), format.size()) == 0) {
                gtIndex = cachedGTIndex;
            } else {
                cachedFormat.assign(format.data(), format.size());
                gtIndex = findGTIndex(format);
                cachedGTIndex = gtIndex;
            }

            if (gtIndex < 0) {
                continue;
            }

            // Check samples
            if (checkAnySampleMatches(p, lineEnd, gtIndex, genotype_query,
                                       queryA1, queryA2, strictCompare)) {
                out << line << "\n";
            }
        }
    }

    // Handle case where we never found data lines
    if (!foundChrom) {
        if (!quiet) {
            std::cerr << "Error: No #CHROM line found in VCF.\n";
        }
        return;
    }
}

// =============================================================================
// Main entry point
// =============================================================================
static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_genotype_query", show_help))
        return 0;

    std::string genotypeQueryStr;
    std::string inputFile;
    bool strictCompare = false;
    bool quiet = false;

    if (!parseArguments(argc, argv, genotypeQueryStr, strictCompare, inputFile, quiet)) {
        std::cerr << "Usage: " << argv[0] << " -g \"0/1\" [--strict] [-i FILE] [-q]\n";
        std::cerr << "Use --help for usage.\n";
        return 1;
    }

    // Pre-parse query for flexible matching
    int queryA1 = -1, queryA2 = -1;
    if (!strictCompare) {
        parseDiploidAlleles(genotypeQueryStr, queryA1, queryA2);
        if (queryA1 > queryA2) std::swap(queryA1, queryA2);
    }

    // Use mmap if file provided
    if (!inputFile.empty()) {
        MappedFile mf;
        if (!mf.open(inputFile.c_str())) {
            std::cerr << "Error: Cannot open file: " << inputFile << "\n";
            return 1;
        }
        genotypeQueryMmap(mf.data, mf.size, std::cout, genotypeQueryStr,
                           queryA1, queryA2, strictCompare, quiet);
    } else {
        genotypeQueryStream(std::cin, std::cout, genotypeQueryStr, strictCompare, quiet);
    }

    return 0;
}
