#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <getopt.h>
#include <iostream>
#include <string>
#include <string_view>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>

// SIMD detection
#if defined(__x86_64__) || defined(_M_X64)
#if defined(__AVX2__)
#define USE_AVX2
#include <immintrin.h>
#elif defined(__SSE2__)
#define USE_SSE2
#include <emmintrin.h>
#endif
#endif

#if defined(__GNUC__) || defined(__clang__)
#define LIKELY(x) __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

// ============================================================================
// MappedFile
// ============================================================================
struct MappedFile {
    const char *data = nullptr;
    size_t size = 0;
    int fd = -1;

    bool open(const char *path) {
        fd = ::open(path, O_RDONLY);
        if (fd < 0) return false;
        struct stat st;
        if (fstat(fd, &st) < 0) { close(); return false; }
        size = static_cast<size_t>(st.st_size);
        if (size == 0) return true;
        data = static_cast<const char *>(mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0));
        if (data == MAP_FAILED) { data = nullptr; close(); return false; }
        madvise(const_cast<char *>(data), size, MADV_SEQUENTIAL | MADV_WILLNEED);
        return true;
    }

    void close() {
        if (data && size > 0) { munmap(const_cast<char *>(data), size); data = nullptr; }
        if (fd >= 0) { ::close(fd); fd = -1; }
        size = 0;
    }

    ~MappedFile() { close(); }
};

// ============================================================================
// OutputBuffer - 1MB buffered output
// ============================================================================
class OutputBuffer {
    static constexpr size_t BUFFER_SIZE = 1024 * 1024;
    char *buffer;
    size_t pos = 0;
    std::ostream &out;

  public:
    explicit OutputBuffer(std::ostream &os) : out(os) { buffer = new char[BUFFER_SIZE]; }
    ~OutputBuffer() { flush(); delete[] buffer; }

    void write(const char *data, size_t len) {
        if (pos + len > BUFFER_SIZE) flush();
        if (len > BUFFER_SIZE) { out.write(data, static_cast<std::streamsize>(len)); return; }
        memcpy(buffer + pos, data, len);
        pos += len;
    }

    void write(std::string_view sv) { write(sv.data(), sv.size()); }
    void writeChar(char c) { if (UNLIKELY(pos >= BUFFER_SIZE)) flush(); buffer[pos++] = c; }
    void flush() { if (pos > 0) { out.write(buffer, static_cast<std::streamsize>(pos)); pos = 0; } }
};

// ============================================================================
// SIMD newline scanning
// ============================================================================
#if defined(USE_AVX2)
static inline const char *findNewline(const char *p, const char *end) {
    const __m256i nl = _mm256_set1_epi8('\n');
    while (p + 32 <= end) {
        __m256i chunk = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(p));
        int mask = _mm256_movemask_epi8(_mm256_cmpeq_epi8(chunk, nl));
        if (mask) return p + __builtin_ctz(static_cast<unsigned>(mask));
        p += 32;
    }
    return static_cast<const char *>(memchr(p, '\n', static_cast<size_t>(end - p)));
}
#elif defined(USE_SSE2)
static inline const char *findNewline(const char *p, const char *end) {
    const __m128i nl = _mm_set1_epi8('\n');
    while (p + 16 <= end) {
        __m128i chunk = _mm_loadu_si128(reinterpret_cast<const __m128i *>(p));
        int mask = _mm_movemask_epi8(_mm_cmpeq_epi8(chunk, nl));
        if (mask) return p + __builtin_ctz(static_cast<unsigned>(mask));
        p += 16;
    }
    return static_cast<const char *>(memchr(p, '\n', static_cast<size_t>(end - p)));
}
#else
static inline const char *findNewline(const char *p, const char *end) {
    return static_cast<const char *>(memchr(p, '\n', static_cast<size_t>(end - p)));
}
#endif

// ============================================================================
// Zero-copy parsing helpers
// ============================================================================
static inline const char* skipToTab(const char* p, const char* end) {
    while (p < end && *p != '\t') ++p;
    return p;
}

static inline std::string_view getNthField(const char* line, const char* lineEnd, int n) {
    const char* p = line;
    for (int i = 0; i < n && p < lineEnd; ++i) {
        while (p < lineEnd && *p != '\t') ++p;
        if (p < lineEnd) ++p;
    }
    const char* start = p;
    while (p < lineEnd && *p != '\t') ++p;
    return std::string_view(start, static_cast<size_t>(p - start));
}

// ============================================================================
// Normalize diploid genotype (zero-copy)
// Returns empty string_view if invalid/missing
// ============================================================================
static std::string normalizeDiploidGenotype(std::string_view gtField, int numAltAlleles) {
    // Extract GT (before first ':')
    size_t colonPos = gtField.find(':');
    std::string_view gt = (colonPos != std::string_view::npos) ? gtField.substr(0, colonPos) : gtField;

    if (gt.empty() || gt[0] == '.') return "";

    // Parse alleles (handle both '/' and '|')
    int a1 = -1, a2 = -1;
    const char* p = gt.data();
    const char* end = p + gt.size();

    // First allele
    if (*p == '.') return "";
    while (p < end && *p >= '0' && *p <= '9') {
        if (a1 < 0) a1 = 0;
        a1 = a1 * 10 + (*p - '0');
        ++p;
    }
    if (a1 < 0 || p >= end) return "";
    if (*p != '/' && *p != '|') return "";
    ++p;

    // Second allele
    if (p >= end || *p == '.') return "";
    while (p < end && *p >= '0' && *p <= '9') {
        if (a2 < 0) a2 = 0;
        a2 = a2 * 10 + (*p - '0');
        ++p;
    }
    if (a2 < 0) return "";

    // Validate allele indices
    if (a1 > numAltAlleles || a2 > numAltAlleles) return "";

    // Sort for comparison (0/1 == 1/0)
    if (a1 > a2) std::swap(a1, a2);

    // Build normalized string
    char buf[16];
    int len = snprintf(buf, sizeof(buf), "%d/%d", a1, a2);
    return std::string(buf, static_cast<size_t>(len));
}

// Count ALT alleles (number of commas + 1, unless empty or ".")
static int countAltAlleles(std::string_view alt) {
    if (alt.empty() || alt == ".") return 0;
    int count = 1;
    for (char c : alt) {
        if (c == ',') ++count;
    }
    return count;
}

// ============================================================================
// Data structure for command-line arguments
// ============================================================================
struct ConcordanceArguments {
    std::string sample1;
    std::string sample2;
    std::string inputFile;
    bool quiet = false;
};

// ============================================================================
// Print help
// ============================================================================
static void printHelp() {
    std::cout << "VCFX_concordance_checker\n"
              << "Usage: VCFX_concordance_checker [OPTIONS] < input.vcf > concordance_report.tsv\n\n"
              << "Options:\n"
              << "  -s, --samples \"Sample1 Sample2\"  Specify exactly two sample names to compare.\n"
              << "  -i, --input FILE                  Input VCF file (uses mmap for best performance)\n"
              << "  -q, --quiet                       Suppress warning messages\n"
              << "  -h, --help                        Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Compares genotypes between two specified samples in a VCF file, including multi-allelic\n"
              << "  variants, and outputs per-variant concordance (Concordant or Discordant).\n\n"
              << "Performance:\n"
              << "  Uses memory-mapped I/O with SIMD-accelerated parsing for ~20-50x speedup.\n\n"
              << "Example:\n"
              << "  VCFX_concordance_checker -s \"SampleA SampleB\" -i input.vcf > concordance_report.tsv\n"
              << "  VCFX_concordance_checker --samples \"SampleA SampleB\" < input.vcf > report.tsv\n";
}

// ============================================================================
// Parse arguments
// ============================================================================
static bool parseArguments(int argc, char *argv[], ConcordanceArguments &args) {
    static struct option long_options[] = {
        {"samples", required_argument, nullptr, 's'},
        {"input", required_argument, nullptr, 'i'},
        {"quiet", no_argument, nullptr, 'q'},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    optind = 1;
    int opt;
    while ((opt = getopt_long(argc, argv, "s:i:qh", long_options, nullptr)) != -1) {
        switch (opt) {
        case 's': {
            std::string samples_str = optarg;
            // Split by space
            size_t space = samples_str.find(' ');
            if (space == std::string::npos) {
                std::cerr << "Error: Please specify exactly two sample names (e.g., -s \"Sample1 Sample2\").\n";
                return false;
            }
            args.sample1 = samples_str.substr(0, space);
            args.sample2 = samples_str.substr(space + 1);
            // Trim trailing spaces from sample2
            while (!args.sample2.empty() && args.sample2.back() == ' ')
                args.sample2.pop_back();
            break;
        }
        case 'i':
            args.inputFile = optarg;
            break;
        case 'q':
            args.quiet = true;
            break;
        case 'h':
            printHelp();
            return false;
        default:
            return false;
        }
    }

    // Check for positional argument
    if (args.inputFile.empty() && optind < argc) {
        args.inputFile = argv[optind];
    }

    if (args.sample1.empty() || args.sample2.empty()) {
        std::cerr << "Error: Two sample names must be specified using --samples or -s.\n";
        std::cerr << "Use --help for usage information.\n";
        return false;
    }

    return true;
}

// ============================================================================
// Calculate concordance with mmap (optimized)
// ============================================================================
static bool calculateConcordanceMmap(const char* filename, std::ostream &out, const ConcordanceArguments &args) {
    MappedFile vcf;
    if (!vcf.open(filename)) {
        std::cerr << "Error: Cannot open file " << filename << "\n";
        return false;
    }
    if (vcf.size == 0) return true;

    const char *p = vcf.data;
    const char *end = vcf.data + vcf.size;

    int sample1_index = -1, sample2_index = -1;
    bool foundHeader = false;

    // Stats
    int totalVariants = 0, concordant = 0, discordant = 0;

    OutputBuffer outBuf(out);

    // Process file
    while (p < end) {
        const char *lineEnd = findNewline(p, end);
        if (!lineEnd) lineEnd = end;
        const char *lineRealEnd = (lineEnd > p && *(lineEnd - 1) == '\r') ? lineEnd - 1 : lineEnd;

        if (*p == '#') {
            // Check for #CHROM header
            if (lineRealEnd - p >= 6 && memcmp(p, "#CHROM", 6) == 0) {
                // Parse sample columns
                int col = 0;
                const char *lp = p;
                while (lp < lineRealEnd) {
                    const char *fieldStart = lp;
                    while (lp < lineRealEnd && *lp != '\t') ++lp;
                    if (col >= 9) {
                        std::string_view sampleName(fieldStart, static_cast<size_t>(lp - fieldStart));
                        if (sampleName == args.sample1) sample1_index = col;
                        if (sampleName == args.sample2) sample2_index = col;
                    }
                    if (lp < lineRealEnd) ++lp;
                    ++col;
                }

                if (sample1_index < 0) {
                    std::cerr << "Error: Sample '" << args.sample1 << "' not found in VCF header.\n";
                    return false;
                }
                if (sample2_index < 0) {
                    std::cerr << "Error: Sample '" << args.sample2 << "' not found in VCF header.\n";
                    return false;
                }
                foundHeader = true;

                // Write output header
                outBuf.write("CHROM\tPOS\tID\tREF\tALT\t");
                outBuf.write(args.sample1);
                outBuf.write("_GT\t");
                outBuf.write(args.sample2);
                outBuf.write("_GT\tConcordance\n");
            }
            p = lineEnd + 1;
            continue;
        }

        if (!foundHeader) {
            std::cerr << "Error: VCF data encountered before #CHROM header.\n";
            return false;
        }

        // Parse data line - extract fields 0-4 and sample columns
        std::string_view chrom = getNthField(p, lineRealEnd, 0);
        std::string_view pos = getNthField(p, lineRealEnd, 1);
        std::string_view id = getNthField(p, lineRealEnd, 2);
        std::string_view ref = getNthField(p, lineRealEnd, 3);
        std::string_view alt = getNthField(p, lineRealEnd, 4);
        std::string_view sample1_field = getNthField(p, lineRealEnd, sample1_index);
        std::string_view sample2_field = getNthField(p, lineRealEnd, sample2_index);

        if (chrom.empty() || sample1_field.empty() || sample2_field.empty()) {
            p = lineEnd + 1;
            continue;
        }

        int numAlt = countAltAlleles(alt);
        std::string s1_gt = normalizeDiploidGenotype(sample1_field, numAlt);
        std::string s2_gt = normalizeDiploidGenotype(sample2_field, numAlt);

        if (s1_gt.empty() || s2_gt.empty()) {
            p = lineEnd + 1;
            continue;
        }

        totalVariants++;
        bool same = (s1_gt == s2_gt);
        if (same) concordant++;
        else discordant++;

        // Output row
        outBuf.write(chrom);
        outBuf.writeChar('\t');
        outBuf.write(pos);
        outBuf.writeChar('\t');
        outBuf.write(id);
        outBuf.writeChar('\t');
        outBuf.write(ref);
        outBuf.writeChar('\t');
        outBuf.write(alt);
        outBuf.writeChar('\t');
        outBuf.write(s1_gt);
        outBuf.writeChar('\t');
        outBuf.write(s2_gt);
        outBuf.writeChar('\t');
        outBuf.write(same ? "Concordant" : "Discordant");
        outBuf.writeChar('\n');

        p = lineEnd + 1;
    }

    outBuf.flush();

    // Print summary to stderr
    if (!args.quiet) {
        std::cerr << "Total Variants Compared: " << totalVariants << "\n"
                  << "Concordant Genotypes: " << concordant << "\n"
                  << "Discordant Genotypes: " << discordant << "\n";
    }

    return true;
}

// ============================================================================
// Calculate concordance from stdin (fallback)
// ============================================================================
static bool calculateConcordance(std::istream &in, std::ostream &out, const ConcordanceArguments &args) {
    std::string line;
    bool foundChromHeader = false;
    std::vector<std::string> headerFields;
    int sample1_index = -1;
    int sample2_index = -1;

    int totalVariants = 0;
    int concordant = 0;
    int discordant = 0;

    out << "CHROM\tPOS\tID\tREF\tALT\t" << args.sample1 << "_GT\t" << args.sample2 << "_GT\tConcordance\n";

    std::vector<std::string> fields;
    fields.reserve(16);

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            if (!foundChromHeader && line.rfind("#CHROM", 0) == 0) {
                vcfx::split_tabs(line, headerFields);
                std::unordered_map<std::string, int> sampleMap;
                for (size_t i = 9; i < headerFields.size(); ++i) {
                    sampleMap[headerFields[i]] = static_cast<int>(i);
                }
                auto it1 = sampleMap.find(args.sample1);
                if (it1 == sampleMap.end()) {
                    std::cerr << "Error: Sample '" << args.sample1 << "' not found in VCF header.\n";
                    return false;
                }
                auto it2 = sampleMap.find(args.sample2);
                if (it2 == sampleMap.end()) {
                    std::cerr << "Error: Sample '" << args.sample2 << "' not found in VCF header.\n";
                    return false;
                }
                sample1_index = it1->second;
                sample2_index = it2->second;
                foundChromHeader = true;
            }
            continue;
        }

        if (!foundChromHeader) {
            std::cerr << "Error: VCF data encountered before #CHROM header.\n";
            return false;
        }

        vcfx::split_tabs(line, fields);
        if (fields.size() < 8) continue;
        if (static_cast<size_t>(sample1_index) >= fields.size() ||
            static_cast<size_t>(sample2_index) >= fields.size()) continue;

        const std::string &alt = fields[4];
        int numAlt = countAltAlleles(alt);

        std::string s1_gt = normalizeDiploidGenotype(fields[static_cast<size_t>(sample1_index)], numAlt);
        std::string s2_gt = normalizeDiploidGenotype(fields[static_cast<size_t>(sample2_index)], numAlt);

        if (s1_gt.empty() || s2_gt.empty()) continue;

        totalVariants++;
        bool same = (s1_gt == s2_gt);
        if (same) concordant++;
        else discordant++;

        out << fields[0] << "\t" << fields[1] << "\t" << fields[2] << "\t"
            << fields[3] << "\t" << fields[4] << "\t" << s1_gt << "\t" << s2_gt
            << "\t" << (same ? "Concordant" : "Discordant") << "\n";
    }

    if (!args.quiet) {
        std::cerr << "Total Variants Compared: " << totalVariants << "\n"
                  << "Concordant Genotypes: " << concordant << "\n"
                  << "Discordant Genotypes: " << discordant << "\n";
    }

    return true;
}

// ============================================================================
// main
// ============================================================================
static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();
    if (vcfx::handle_common_flags(argc, argv, "VCFX_concordance_checker", show_help))
        return 0;

    ConcordanceArguments args;
    if (!parseArguments(argc, argv, args)) {
        return 1;
    }

    bool ok;
    if (!args.inputFile.empty()) {
        ok = calculateConcordanceMmap(args.inputFile.c_str(), std::cout, args);
    } else {
        ok = calculateConcordance(std::cin, std::cout, args);
    }

    return ok ? 0 : 1;
}
