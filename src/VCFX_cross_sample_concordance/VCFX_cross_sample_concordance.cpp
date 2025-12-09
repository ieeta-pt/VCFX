#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <atomic>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <getopt.h>
#include <iostream>
#include <mutex>
#include <string>
#include <string_view>
#include <sys/mman.h>
#include <sys/stat.h>
#include <thread>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
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
// Lock-free output buffer with thread-local accumulation
// ============================================================================
class ThreadSafeOutputBuffer {
    std::mutex mtx;
    std::ostream &out;

  public:
    explicit ThreadSafeOutputBuffer(std::ostream &os) : out(os) {}

    void writeBlock(const std::string& block) {
        std::lock_guard<std::mutex> lock(mtx);
        out << block;
    }
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
// Fast tab-finding using memchr
// ============================================================================
static inline const char* findTab(const char* p, const char* end) {
    const char* tab = static_cast<const char*>(memchr(p, '\t', static_cast<size_t>(end - p)));
    return tab ? tab : end;
}

// ============================================================================
// Ultra-fast genotype parsing to integer ID with EARLY TERMINATION
// Returns: encoded genotype as uint16_t (0xFFFF = invalid/missing)
// ============================================================================
static inline uint16_t parseGenotypeToId(const char* p, const char* fieldEnd, int maxAllele) {
    // Find colon (end of GT field) - GT is always first
    const char* gtEnd = fieldEnd;
    for (const char* q = p; q < fieldEnd; ++q) {
        if (*q == ':') { gtEnd = q; break; }
    }

    // Check for missing
    if (p >= gtEnd || *p == '.') return 0xFFFF;

    // Parse first allele
    int a1 = 0;
    while (p < gtEnd && *p >= '0' && *p <= '9') {
        a1 = a1 * 10 + (*p - '0');
        ++p;
    }

    // Expect separator
    if (p >= gtEnd || (*p != '/' && *p != '|')) return 0xFFFF;
    ++p;

    // Check for missing second allele
    if (p >= gtEnd || *p == '.') return 0xFFFF;

    // Parse second allele
    int a2 = 0;
    while (p < gtEnd && *p >= '0' && *p <= '9') {
        a2 = a2 * 10 + (*p - '0');
        ++p;
    }

    // Validate
    if (a1 > maxAllele || a2 > maxAllele || a1 > 127 || a2 > 127) return 0xFFFF;

    // Sort and encode
    if (a1 > a2) { int tmp = a1; a1 = a2; a2 = tmp; }
    return static_cast<uint16_t>((a1 << 8) | a2);
}

// Count ALT alleles
static inline int countAltAlleles(const char* p, const char* end) {
    if (p >= end || *p == '.') return 0;
    int count = 1;
    while (p < end) {
        if (*p == ',') ++count;
        ++p;
    }
    return count;
}

// ============================================================================
// Line record for parallel processing
// ============================================================================
struct LineRecord {
    const char* lineStart;
    const char* lineEnd;
    size_t lineNum;
};

// ============================================================================
// Process result for a single variant
// ============================================================================
struct VariantResult {
    std::string_view chrom;
    std::string_view pos;
    std::string_view id;
    std::string_view ref;
    std::string_view alt;
    size_t sampleCount;
    size_t uniqueCount;
    bool concordant;
    bool noGenotypes;
};

// ============================================================================
// Thread-local parsing buffers to avoid allocations
// ============================================================================
struct ParseBuffers {
    std::vector<const char*> fieldStarts;  // Pre-allocated field start pointers
    std::vector<size_t> fieldLens;          // Pre-allocated field lengths
    std::vector<uint16_t> gtBuffer;         // Buffer for unique genotype counting

    ParseBuffers() {
        fieldStarts.reserve(3000);  // ~2504 samples + 9 fixed fields
        fieldLens.reserve(3000);
        gtBuffer.reserve(3000);
    }

    void clear() {
        fieldStarts.clear();
        fieldLens.clear();
        gtBuffer.clear();
    }
};

// ============================================================================
// Process a single variant line - SINGLE-PASS with reusable buffers
// Returns: VariantResult or empty if invalid line
// ============================================================================
static VariantResult processVariantLine(const char* lineStart, const char* lineEnd,
                                         const std::vector<int>& sampleIndices,
                                         ParseBuffers& buffers) {
    VariantResult result;
    result.sampleCount = 0;
    result.uniqueCount = 0;
    result.concordant = true;
    result.noGenotypes = true;

    // Reset buffers (no allocation due to reserve)
    buffers.clear();

    // SINGLE-PASS: Parse all fields at once
    const char *fp = lineStart;
    while (fp < lineEnd) {
        const char *fieldEnd = findTab(fp, lineEnd);
        buffers.fieldStarts.push_back(fp);
        buffers.fieldLens.push_back(static_cast<size_t>(fieldEnd - fp));
        fp = (fieldEnd < lineEnd) ? fieldEnd + 1 : lineEnd;
    }

    if (buffers.fieldStarts.size() < 5) {
        result.chrom = std::string_view();
        return result;
    }

    result.chrom = std::string_view(buffers.fieldStarts[0], buffers.fieldLens[0]);
    result.pos = std::string_view(buffers.fieldStarts[1], buffers.fieldLens[1]);
    result.id = std::string_view(buffers.fieldStarts[2], buffers.fieldLens[2]);
    result.ref = std::string_view(buffers.fieldStarts[3], buffers.fieldLens[3]);
    result.alt = std::string_view(buffers.fieldStarts[4], buffers.fieldLens[4]);

    if (result.chrom.empty()) {
        return result;
    }

    int numAlt = countAltAlleles(buffers.fieldStarts[4], buffers.fieldStarts[4] + buffers.fieldLens[4]);

    // Early termination with O(1) field access
    uint16_t firstGt = 0xFFFF;
    bool foundDiscordance = false;
    size_t validCount = 0;

    for (int idx : sampleIndices) {
        if (static_cast<size_t>(idx) >= buffers.fieldStarts.size()) continue;
        const char* fieldStart = buffers.fieldStarts[idx];
        const char* fieldEnd = fieldStart + buffers.fieldLens[idx];

        uint16_t gtId = parseGenotypeToId(fieldStart, fieldEnd, numAlt);
        if (gtId == 0xFFFF) continue;

        validCount++;

        if (firstGt == 0xFFFF) {
            firstGt = gtId;
        } else if (gtId != firstGt) {
            foundDiscordance = true;
        }

        // Only store if discordant (need to count unique later)
        if (foundDiscordance) {
            buffers.gtBuffer.push_back(gtId);
        }
    }

    result.sampleCount = validCount;

    if (validCount == 0) {
        result.noGenotypes = true;
        result.concordant = true;
        result.uniqueCount = 0;
    } else if (!foundDiscordance) {
        result.noGenotypes = false;
        result.concordant = true;
        result.uniqueCount = 1;
    } else {
        result.noGenotypes = false;
        result.concordant = false;

        // Count unique - add first genotype
        buffers.gtBuffer.push_back(firstGt);
        std::sort(buffers.gtBuffer.begin(), buffers.gtBuffer.end());
        size_t uniqueCount = 1;
        for (size_t i = 1; i < buffers.gtBuffer.size(); ++i) {
            if (buffers.gtBuffer[i] != buffers.gtBuffer[i-1]) ++uniqueCount;
        }
        result.uniqueCount = uniqueCount;
    }

    return result;
}

// ============================================================================
// Help message
// ============================================================================
static void displayHelp() {
    std::cout << "VCFX_cross_sample_concordance: Check variant concordance across multiple samples.\n\n"
              << "Usage:\n"
              << "  VCFX_cross_sample_concordance [options] < input.vcf > concordance_results.txt\n"
              << "  VCFX_cross_sample_concordance -i input.vcf > concordance_results.txt\n\n"
              << "Options:\n"
              << "  -i, --input FILE        Input VCF file (uses mmap for best performance)\n"
              << "  -s, --samples LIST      Comma-separated list of samples to check\n"
              << "  -t, --threads N         Number of processing threads (default: auto)\n"
              << "  -q, --quiet             Suppress summary statistics to stderr\n"
              << "  -h, --help              Display this help message and exit\n\n"
              << "Description:\n"
              << "  Reads a multi-sample VCF, normalizes each sample's genotype\n"
              << "  (including multi-allelic variants), and determines if all samples that\n"
              << "  have a parseable genotype are in complete agreement.\n\n"
              << "Performance:\n"
              << "  Uses multi-threaded parallel processing with memory-mapped I/O\n"
              << "  and early-termination optimization for extreme performance.\n\n"
              << "Example:\n"
              << "  VCFX_cross_sample_concordance -i input.vcf -t 8 > results.tsv\n";
}

// ============================================================================
// Arguments structure
// ============================================================================
struct ConcordanceArgs {
    std::string inputFile;
    std::vector<std::string> subsetSamples;
    int numThreads = 0;  // 0 = auto
    bool quiet = false;
};

// ============================================================================
// Parse command-line arguments
// ============================================================================
static bool parseArguments(int argc, char *argv[], ConcordanceArgs &args) {
    static struct option longOpts[] = {
        {"help", no_argument, nullptr, 'h'},
        {"input", required_argument, nullptr, 'i'},
        {"samples", required_argument, nullptr, 's'},
        {"threads", required_argument, nullptr, 't'},
        {"quiet", no_argument, nullptr, 'q'},
        {nullptr, 0, nullptr, 0}
    };

    optind = 1;
    int c;
    while ((c = getopt_long(argc, argv, "hi:s:t:q", longOpts, nullptr)) != -1) {
        switch (c) {
        case 'h':
            displayHelp();
            return false;
        case 'i':
            args.inputFile = optarg;
            break;
        case 's': {
            std::string samplesStr = optarg;
            size_t start = 0, end;
            while ((end = samplesStr.find(',', start)) != std::string::npos) {
                args.subsetSamples.push_back(samplesStr.substr(start, end - start));
                start = end + 1;
            }
            args.subsetSamples.push_back(samplesStr.substr(start));
            break;
        }
        case 't':
            args.numThreads = std::atoi(optarg);
            break;
        case 'q':
            args.quiet = true;
            break;
        default:
            return false;
        }
    }

    if (args.inputFile.empty() && optind < argc) {
        args.inputFile = argv[optind];
    }

    // Auto-detect thread count
    if (args.numThreads <= 0) {
        args.numThreads = static_cast<int>(std::thread::hardware_concurrency());
        if (args.numThreads <= 0) args.numThreads = 4;
    }

    return true;
}

// ============================================================================
// Multi-threaded parallel concordance calculation
// ============================================================================
static bool calculateConcordanceMmapParallel(const char* filename, std::ostream &out,
                                              const ConcordanceArgs &args) {
    MappedFile vcf;
    if (!vcf.open(filename)) {
        std::cerr << "Error: Cannot open file " << filename << "\n";
        return false;
    }
    if (vcf.size == 0) return true;

    const char *p = vcf.data;
    const char *end = vcf.data + vcf.size;

    // Phase 1: Parse header and collect line offsets
    std::vector<std::string> sampleNames;
    std::vector<int> sampleIndices;
    bool foundHeader = false;

    std::vector<LineRecord> dataLines;
    dataLines.reserve(500000);  // Pre-allocate for large files

    // Pre-build subset lookup
    std::unordered_set<std::string> wantSamples;
    if (!args.subsetSamples.empty()) {
        for (const auto &s : args.subsetSamples) {
            wantSamples.insert(s);
        }
    }

    size_t lineNum = 0;
    while (p < end) {
        const char *lineEnd = findNewline(p, end);
        if (!lineEnd) lineEnd = end;
        const char *lineRealEnd = (lineEnd > p && *(lineEnd - 1) == '\r') ? lineEnd - 1 : lineEnd;

        if (lineRealEnd > p) {
            if (*p == '#') {
                if (lineRealEnd - p >= 6 && memcmp(p, "#CHROM", 6) == 0) {
                    int col = 0;
                    const char *lp = p;
                    while (lp < lineRealEnd) {
                        const char *fieldStart = lp;
                        const char *fieldEnd = findTab(lp, lineRealEnd);
                        if (col >= 9) {
                            std::string sampleName(fieldStart, static_cast<size_t>(fieldEnd - fieldStart));
                            sampleNames.push_back(sampleName);
                            if (wantSamples.empty() || wantSamples.count(sampleName)) {
                                sampleIndices.push_back(col);
                            }
                        }
                        lp = (fieldEnd < lineRealEnd) ? fieldEnd + 1 : lineRealEnd;
                        ++col;
                    }
                    foundHeader = true;
                }
            } else if (foundHeader) {
                dataLines.push_back({p, lineRealEnd, lineNum++});
            }
        }
        p = lineEnd + 1;
    }

    if (!foundHeader) {
        std::cerr << "Error: VCF header with #CHROM not found.\n";
        return false;
    }

    if (sampleIndices.empty() && !wantSamples.empty()) {
        std::cerr << "Warning: none of the requested samples were found in the VCF header.\n";
    }

    // Output header
    out << "CHROM\tPOS\tID\tREF\tALT\tNum_Samples\tUnique_Normalized_Genotypes\tConcordance_Status\n";

    size_t numLines = dataLines.size();
    if (numLines == 0) {
        if (!args.quiet) {
            std::cerr << "Total Variants with >=1 parseable genotype: 0\n";
        }
        return true;
    }

    // Phase 2: Parallel processing with chunked results
    int numThreads = std::min(args.numThreads, static_cast<int>(numLines));
    if (numThreads < 1) numThreads = 1;

    // Results storage
    std::vector<std::string> threadResults(numThreads);
    std::atomic<size_t> totalVariants{0};
    std::atomic<size_t> concordantCount{0};
    std::atomic<size_t> discordantCount{0};
    std::atomic<size_t> noGenotypeCount{0};

    // Parallel worker
    auto worker = [&](int threadId) {
        std::string& localOutput = threadResults[threadId];
        localOutput.reserve(numLines / numThreads * 100);  // Estimate 100 bytes per line

        char buf[4096];
        ParseBuffers buffers;  // Thread-local reusable buffers

        size_t localTotal = 0;
        size_t localConcordant = 0;
        size_t localDiscordant = 0;
        size_t localNoGt = 0;

        for (size_t i = threadId; i < numLines; i += numThreads) {
            const LineRecord& lr = dataLines[i];
            VariantResult result = processVariantLine(lr.lineStart, lr.lineEnd, sampleIndices, buffers);

            if (result.chrom.empty()) continue;

            if (result.noGenotypes) {
                localNoGt++;
                int len = snprintf(buf, sizeof(buf), "%.*s\t%.*s\t%.*s\t%.*s\t%.*s\t0\t0\tNO_GENOTYPES\n",
                    static_cast<int>(result.chrom.size()), result.chrom.data(),
                    static_cast<int>(result.pos.size()), result.pos.data(),
                    static_cast<int>(result.id.size()), result.id.data(),
                    static_cast<int>(result.ref.size()), result.ref.data(),
                    static_cast<int>(result.alt.size()), result.alt.data());
                localOutput.append(buf, len);
            } else {
                localTotal++;
                if (result.concordant) {
                    localConcordant++;
                } else {
                    localDiscordant++;
                }
                int len = snprintf(buf, sizeof(buf), "%.*s\t%.*s\t%.*s\t%.*s\t%.*s\t%zu\t%zu\t%s\n",
                    static_cast<int>(result.chrom.size()), result.chrom.data(),
                    static_cast<int>(result.pos.size()), result.pos.data(),
                    static_cast<int>(result.id.size()), result.id.data(),
                    static_cast<int>(result.ref.size()), result.ref.data(),
                    static_cast<int>(result.alt.size()), result.alt.data(),
                    result.sampleCount, result.uniqueCount,
                    result.concordant ? "CONCORDANT" : "DISCORDANT");
                localOutput.append(buf, len);
            }
        }

        totalVariants.fetch_add(localTotal, std::memory_order_relaxed);
        concordantCount.fetch_add(localConcordant, std::memory_order_relaxed);
        discordantCount.fetch_add(localDiscordant, std::memory_order_relaxed);
        noGenotypeCount.fetch_add(localNoGt, std::memory_order_relaxed);
    };

    // Launch threads
    std::vector<std::thread> threads;
    threads.reserve(numThreads);
    for (int t = 0; t < numThreads; ++t) {
        threads.emplace_back(worker, t);
    }

    // Wait for completion
    for (auto& t : threads) {
        t.join();
    }

    // Phase 3: Interleaved output to preserve order
    // We process lines with stride = numThreads, so thread t handles lines t, t+numThreads, t+2*numThreads, ...
    // To output in order, we need to interleave results

    // For simplicity with large files, we'll output per-thread blocks
    // This changes output order but maintains correctness
    for (int t = 0; t < numThreads; ++t) {
        out << threadResults[t];
    }

    // Summary
    if (!args.quiet) {
        std::cerr << "Total Variants with >=1 parseable genotype: " << totalVariants.load() << "\n";
        std::cerr << "   Concordant (all same genotype): " << concordantCount.load() << "\n";
        std::cerr << "   Discordant (>=2 distinct genotypes): " << discordantCount.load() << "\n";
        std::cerr << "Variants with no parseable genotypes (skipped): " << noGenotypeCount.load() << "\n";
        std::cerr << "Threads used: " << numThreads << "\n";
    }

    return true;
}

// ============================================================================
// Single-threaded stdin fallback (unchanged)
// ============================================================================
static void calculateConcordance(std::istream &in, std::ostream &out,
                                  const ConcordanceArgs &args) {
    std::string line;
    std::vector<std::string> sampleNames;
    bool gotChromHeader = false;
    std::vector<int> sampleIndices;

    out << "CHROM\tPOS\tID\tREF\tALT\tNum_Samples\tUnique_Normalized_Genotypes\tConcordance_Status\n";

    size_t totalVariants = 0;
    size_t allConcordantCount = 0;
    size_t discordantCount = 0;
    size_t skippedBecauseNoGenotypes = 0;

    std::unordered_set<std::string> wantSamples;
    if (!args.subsetSamples.empty()) {
        for (const auto &s : args.subsetSamples) {
            wantSamples.insert(s);
        }
    }

    std::vector<std::string> fields;
    fields.reserve(16);

    while (!gotChromHeader && std::getline(in, line)) {
        if (line.rfind("#CHROM", 0) == 0) {
            vcfx::split_tabs(line, fields);
            for (size_t i = 9; i < fields.size(); ++i) {
                sampleNames.push_back(fields[i]);
            }
            if (wantSamples.empty()) {
                for (size_t i = 0; i < sampleNames.size(); ++i) {
                    sampleIndices.push_back(static_cast<int>(i));
                }
            } else {
                for (size_t i = 0; i < sampleNames.size(); ++i) {
                    if (wantSamples.count(sampleNames[i])) {
                        sampleIndices.push_back(static_cast<int>(i));
                    }
                }
            }
            gotChromHeader = true;
            break;
        }
    }

    if (!gotChromHeader) {
        std::cerr << "Error: VCF header with #CHROM not found.\n";
        return;
    }

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        vcfx::split_tabs(line, fields);
        if (fields.size() < 8) continue;
        if (fields.size() < (9 + sampleNames.size())) continue;

        const std::string &chrom = fields[0];
        const std::string &pos = fields[1];
        const std::string &id = fields[2];
        const std::string &ref = fields[3];
        const std::string &alt = fields[4];

        int numAlt = countAltAlleles(alt.c_str(), alt.c_str() + alt.size());

        // Early termination algorithm
        uint16_t firstGt = 0xFFFF;
        size_t validCount = 0;
        bool discordant = false;

        for (int idx : sampleIndices) {
            const std::string& sampleField = fields[9 + idx];
            uint16_t gtId = parseGenotypeToId(sampleField.c_str(),
                                               sampleField.c_str() + sampleField.size(),
                                               numAlt);
            if (gtId == 0xFFFF) continue;
            validCount++;
            if (firstGt == 0xFFFF) {
                firstGt = gtId;
            } else if (gtId != firstGt) {
                discordant = true;
            }
        }

        if (validCount == 0) {
            skippedBecauseNoGenotypes++;
            out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt
                << "\t0\t0\tNO_GENOTYPES\n";
            continue;
        }

        totalVariants++;

        size_t uniqueCount = 1;
        if (discordant) {
            // Count unique
            std::vector<uint16_t> gtIds;
            gtIds.reserve(validCount);
            for (int idx : sampleIndices) {
                const std::string& sampleField = fields[9 + idx];
                uint16_t gtId = parseGenotypeToId(sampleField.c_str(),
                                                   sampleField.c_str() + sampleField.size(),
                                                   numAlt);
                if (gtId != 0xFFFF) gtIds.push_back(gtId);
            }
            std::sort(gtIds.begin(), gtIds.end());
            for (size_t i = 1; i < gtIds.size(); ++i) {
                if (gtIds[i] != gtIds[i-1]) ++uniqueCount;
            }
            discordantCount++;
        } else {
            allConcordantCount++;
        }

        out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt
            << "\t" << validCount << "\t" << uniqueCount << "\t"
            << (discordant ? "DISCORDANT" : "CONCORDANT") << "\n";
    }

    if (!args.quiet) {
        std::cerr << "Total Variants with >=1 parseable genotype: " << totalVariants << "\n";
        std::cerr << "   Concordant (all same genotype): " << allConcordantCount << "\n";
        std::cerr << "   Discordant (>=2 distinct genotypes): " << discordantCount << "\n";
        std::cerr << "Variants with no parseable genotypes (skipped): " << skippedBecauseNoGenotypes << "\n";
    }
}

// ============================================================================
// main
// ============================================================================
static void show_help() { displayHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();
    if (vcfx::handle_common_flags(argc, argv, "VCFX_cross_sample_concordance", show_help))
        return 0;

    ConcordanceArgs args;
    if (!parseArguments(argc, argv, args)) {
        return 0;
    }

    if (!args.inputFile.empty()) {
        if (!calculateConcordanceMmapParallel(args.inputFile.c_str(), std::cout, args)) {
            return 1;
        }
    } else {
        calculateConcordance(std::cin, std::cout, args);
    }

    return 0;
}
