// VCFX_allele_balance_calc.cpp
// Optimized implementation using mmap, SIMD, and multi-threading
// Based on proven patterns from VCFX_allele_counter

#include "vcfx_core.h"
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <mutex>
#include <string>
#include <string_view>
#include <sys/mman.h>
#include <sys/stat.h>
#include <thread>
#include <unistd.h>
#include <unordered_map>
#include <vector>

#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#define VCFX_HAS_SSE2 1
#if defined(__AVX2__)
#define VCFX_HAS_AVX2 1
#endif
#elif defined(__aarch64__)
#include <arm_neon.h>
#define VCFX_HAS_NEON 1
#endif

// Branch prediction hints
#define LIKELY(x) __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)

// ---------------------------------------------------------------------
// Memory-mapped file structure
// ---------------------------------------------------------------------
struct MappedFile {
    const char *data = nullptr;
    size_t size = 0;
    int fd = -1;

    bool open(const char *path) {
        fd = ::open(path, O_RDONLY);
        if (fd < 0) return false;
        struct stat st;
        if (fstat(fd, &st) < 0) {
            close();
            return false;
        }
        size = static_cast<size_t>(st.st_size);
        if (size == 0) return true;
        data = static_cast<const char *>(
            mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0));
        if (data == MAP_FAILED) {
            data = nullptr;
            close();
            return false;
        }
        madvise(const_cast<char *>(data), size, MADV_SEQUENTIAL | MADV_WILLNEED);
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

// ---------------------------------------------------------------------
// Thread-local output buffer with incremental flushing
// ---------------------------------------------------------------------
class ThreadBuffer {
  public:
    std::vector<char> buffer;
    size_t pos = 0;
    static constexpr size_t BUFFER_SIZE = 4 * 1024 * 1024;  // 4MB buffer
    static constexpr size_t FLUSH_THRESHOLD = 3 * 1024 * 1024;  // Flush at 3MB
    int outFd = STDOUT_FILENO;
    std::mutex *writeMutex = nullptr;  // For synchronized writes

    ThreadBuffer() {
        buffer.resize(BUFFER_SIZE);
    }

    void setOutput(int fd, std::mutex *mtx = nullptr) {
        outFd = fd;
        writeMutex = mtx;
    }

    void flush() {
        if (pos > 0) {
            if (writeMutex) {
                std::lock_guard<std::mutex> lock(*writeMutex);
                ::write(outFd, buffer.data(), pos);
            } else {
                ::write(outFd, buffer.data(), pos);
            }
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
            flush();  // Flush instead of growing
        }
    }

    void write(const char *data, size_t len) {
        // If data is larger than buffer, write directly
        if (len > BUFFER_SIZE / 2) {
            flush();
            if (writeMutex) {
                std::lock_guard<std::mutex> lock(*writeMutex);
                ::write(outFd, data, len);
            } else {
                ::write(outFd, data, len);
            }
            return;
        }
        ensureSpace(len);
        memcpy(buffer.data() + pos, data, len);
        pos += len;
    }

    void write(std::string_view sv) { write(sv.data(), sv.size()); }

    void writeChar(char c) {
        ensureSpace(1);
        buffer[pos++] = c;
    }

    // Fast double formatting with 6 decimal places
    void writeDouble(double val) {
        ensureSpace(24);

        // Handle integer part
        int intPart = static_cast<int>(val);
        double fracPart = val - intPart;

        if (intPart == 0) {
            buffer[pos++] = '0';
        } else {
            char tmp[12];
            int len = 0;
            while (intPart > 0) {
                tmp[len++] = '0' + (intPart % 10);
                intPart /= 10;
            }
            while (len > 0) {
                buffer[pos++] = tmp[--len];
            }
        }

        // Decimal point and 6 digits
        buffer[pos++] = '.';
        for (int i = 0; i < 6; ++i) {
            fracPart *= 10;
            int digit = static_cast<int>(fracPart);
            buffer[pos++] = '0' + digit;
            fracPart -= digit;
        }
    }

    const char *data() const { return buffer.data(); }
    size_t size() const { return pos; }
    void clear() { pos = 0; }
};

// ---------------------------------------------------------------------
// SIMD newline scanning
// ---------------------------------------------------------------------
#if defined(VCFX_HAS_AVX2)
static inline const char *findNewline(const char *p, const char *end) {
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
    while (p < end && *p != '\n') p++;
    return p;
}
#elif defined(VCFX_HAS_SSE2)
static inline const char *findNewline(const char *p, const char *end) {
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
    while (p < end && *p != '\n') p++;
    return p;
}
#elif defined(VCFX_HAS_NEON)
static inline const char *findNewline(const char *p, const char *end) {
    const uint8x16_t nl = vdupq_n_u8('\n');
    while (p + 16 <= end) {
        uint8x16_t chunk = vld1q_u8(reinterpret_cast<const uint8_t *>(p));
        uint8x16_t cmp = vceqq_u8(chunk, nl);
        uint64x2_t combined = vreinterpretq_u64_u8(cmp);
        uint64_t low = vgetq_lane_u64(combined, 0);
        uint64_t high = vgetq_lane_u64(combined, 1);
        if (low) return p + (__builtin_ctzll(low) >> 3);
        if (high) return p + 8 + (__builtin_ctzll(high) >> 3);
        p += 16;
    }
    while (p < end && *p != '\n') p++;
    return p;
}
#else
static inline const char *findNewline(const char *p, const char *end) {
    const char *found = static_cast<const char *>(memchr(p, '\n', end - p));
    return found ? found : end;
}
#endif

// ---------------------------------------------------------------------
// SIMD tab scanning
// ---------------------------------------------------------------------
#if defined(VCFX_HAS_AVX2)
static inline const char *findTab(const char *p, const char *end) {
    const __m256i tab = _mm256_set1_epi8('\t');
    while (p + 32 <= end) {
        __m256i chunk = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(p));
        __m256i cmp = _mm256_cmpeq_epi8(chunk, tab);
        int mask = _mm256_movemask_epi8(cmp);
        if (mask) {
            return p + __builtin_ctz(static_cast<unsigned int>(mask));
        }
        p += 32;
    }
    while (p < end && *p != '\t') p++;
    return p;
}
#elif defined(VCFX_HAS_SSE2)
static inline const char *findTab(const char *p, const char *end) {
    const __m128i tab = _mm_set1_epi8('\t');
    while (p + 16 <= end) {
        __m128i chunk = _mm_loadu_si128(reinterpret_cast<const __m128i *>(p));
        __m128i cmp = _mm_cmpeq_epi8(chunk, tab);
        int mask = _mm_movemask_epi8(cmp);
        if (mask) {
            return p + __builtin_ctz(static_cast<unsigned int>(mask));
        }
        p += 16;
    }
    while (p < end && *p != '\t') p++;
    return p;
}
#elif defined(VCFX_HAS_NEON)
static inline const char *findTab(const char *p, const char *end) {
    const uint8x16_t tab = vdupq_n_u8('\t');
    while (p + 16 <= end) {
        uint8x16_t chunk = vld1q_u8(reinterpret_cast<const uint8_t *>(p));
        uint8x16_t cmp = vceqq_u8(chunk, tab);
        uint64x2_t combined = vreinterpretq_u64_u8(cmp);
        uint64_t low = vgetq_lane_u64(combined, 0);
        uint64_t high = vgetq_lane_u64(combined, 1);
        if (low) return p + (__builtin_ctzll(low) >> 3);
        if (high) return p + 8 + (__builtin_ctzll(high) >> 3);
        p += 16;
    }
    while (p < end && *p != '\t') p++;
    return p;
}
#else
static inline const char *findTab(const char *p, const char *end) {
    const char *found = static_cast<const char *>(memchr(p, '\t', end - p));
    return found ? found : end;
}
#endif

// ---------------------------------------------------------------------
// Zero-copy field extraction
// ---------------------------------------------------------------------
static inline std::string_view extractField(const char *&p, const char *lineEnd) {
    const char *start = p;
    p = findTab(p, lineEnd);
    return {start, static_cast<size_t>(p - start)};
}

static inline void skipFields(const char *&p, const char *lineEnd, int n) {
    for (int i = 0; i < n && p < lineEnd; ++i) {
        p = findTab(p, lineEnd);
        if (p < lineEnd) ++p;
    }
}

// ---------------------------------------------------------------------
// Ultra-fast genotype parsing - returns allele balance
// Returns: >= 0 for valid balance, -1 for missing/invalid
// ---------------------------------------------------------------------
static inline const char *findColon(const char *p, const char *end) {
    while (p < end && *p != ':') ++p;
    return p;
}

static inline double computeAlleleBalanceFast(const char *gt, const char *gtEnd) {
    int refCount = 0;
    int altCount = 0;

    while (gt < gtEnd) {
        // Skip separators
        while (gt < gtEnd && (*gt == '/' || *gt == '|')) ++gt;
        if (gt >= gtEnd) break;

        // Missing allele
        if (*gt == '.') {
            ++gt;
            continue;
        }

        // Parse allele number
        int allele = 0;
        bool hasDigit = false;
        while (gt < gtEnd && *gt >= '0' && *gt <= '9') {
            allele = allele * 10 + (*gt - '0');
            hasDigit = true;
            ++gt;
        }

        if (LIKELY(hasDigit)) {
            if (allele == 0) ++refCount;
            else ++altCount;
        }
    }

    // Compute allele balance
    if (altCount == 0 && refCount > 0) return 0.0;
    if (refCount + altCount == 0) return -1.0;  // Missing
    return static_cast<double>(refCount) / altCount;
}

// ---------------------------------------------------------------------
// Arguments structure
// ---------------------------------------------------------------------
struct AlleleBalanceArgs {
    std::vector<std::string> samples;
    const char *inputFile = nullptr;
    bool quiet = false;
    int numThreads = 0; // 0 = auto-detect
};

// ---------------------------------------------------------------------
// Parse arguments
// ---------------------------------------------------------------------
static bool parseArguments(int argc, char *argv[], AlleleBalanceArgs &args) {
    for (int i = 1; i < argc; ++i) {
        std::string_view arg = argv[i];
        if ((arg == "--samples" || arg == "-s") && i + 1 < argc) {
            std::string samplesStr = argv[++i];
            size_t start = 0, end;
            while ((end = samplesStr.find(' ', start)) != std::string::npos) {
                if (end > start) {
                    args.samples.emplace_back(samplesStr.substr(start, end - start));
                }
                start = end + 1;
            }
            if (start < samplesStr.size()) {
                args.samples.emplace_back(samplesStr.substr(start));
            }
            for (auto &sample : args.samples) {
                size_t first = sample.find_first_not_of(" \t\n\r");
                size_t last = sample.find_last_not_of(" \t\n\r");
                if (first != std::string::npos) {
                    sample = sample.substr(first, last - first + 1);
                }
            }
        } else if (arg == "--input" || arg == "-i") {
            if (i + 1 < argc) args.inputFile = argv[++i];
        } else if (arg == "--threads" || arg == "-t") {
            if (i + 1 < argc) args.numThreads = std::atoi(argv[++i]);
        } else if (arg == "--quiet" || arg == "-q") {
            args.quiet = true;
        } else if (arg == "--help" || arg == "-h") {
            return false;
        } else if (arg[0] != '-' && !args.inputFile) {
            args.inputFile = argv[i];
        }
    }
    return true;
}

// ---------------------------------------------------------------------
// Print help
// ---------------------------------------------------------------------
static void printHelp() {
    std::cout
        << "VCFX_allele_balance_calc - Calculate allele balance (ref/alt ratio) per sample\n\n"
        << "Usage: VCFX_allele_balance_calc [OPTIONS] [FILE]\n\n"
        << "Options:\n"
        << "  -i, --input FILE    Input VCF file (uses mmap for best performance)\n"
        << "  -t, --threads N     Number of threads (default: auto-detect CPU cores)\n"
        << "  -s, --samples STR   Space-separated list of sample names to include\n"
        << "  -q, --quiet         Suppress informational messages\n"
        << "  -h, --help          Display this help message\n"
        << "  -v, --version       Display version information\n\n"
        << "Description:\n"
        << "  Calculates the allele balance (ratio of reference to alternate alleles) for\n"
        << "  each sample at each variant. Allele balance = #RefAlleles / #AltAlleles.\n"
        << "  Missing genotypes produce \"NA\" output.\n\n"
        << "Examples:\n"
        << "  VCFX_allele_balance_calc -i input.vcf > balance.tsv           # Auto threads\n"
        << "  VCFX_allele_balance_calc -t 8 -i input.vcf > balance.tsv      # 8 threads\n"
        << "  VCFX_allele_balance_calc < input.vcf > balance.tsv            # Stdin (single-thread)\n\n"
        << "Output format:\n"
        << "  CHROM  POS  ID  REF  ALT  Sample  Allele_Balance\n";
}

// ---------------------------------------------------------------------
// Find all sample start positions in a line
// ---------------------------------------------------------------------
static inline void findAllSampleStarts(
    const char *sampleDataStart,
    const char *lineEnd,
    std::vector<const char*> &sampleStarts)
{
    sampleStarts.clear();
    const char *p = sampleDataStart;
    sampleStarts.push_back(p);

    while (p < lineEnd) {
        p = findTab(p, lineEnd);
        if (p < lineEnd) {
            ++p;
            sampleStarts.push_back(p);
        }
    }
}

// ---------------------------------------------------------------------
// Process a chunk of lines (multi-threaded version)
// ---------------------------------------------------------------------
static void processChunk(
    const char *chunkStart,
    const char *chunkEnd,
    const std::vector<std::string> &sampleSuffix,
    const std::vector<size_t> &sampleIndices,
    ThreadBuffer &outBuf)
{
    char prefixBuf[4096];
    const char *p = chunkStart;

    std::vector<const char*> sampleStarts;
    sampleStarts.reserve(3000);

    const size_t numSamples = sampleIndices.size();
    std::vector<double> balances(numSamples);

    while (p < chunkEnd) {
        const char *lineEnd = findNewline(p, chunkEnd);
        if (p >= lineEnd || *p == '#') {
            p = lineEnd;
            if (p < chunkEnd) ++p;
            continue;
        }

        // Extract first 5 fields
        std::string_view chrom = extractField(p, lineEnd);
        if (p < lineEnd) ++p;
        std::string_view pos = extractField(p, lineEnd);
        if (p < lineEnd) ++p;
        std::string_view id = extractField(p, lineEnd);
        if (p < lineEnd) ++p;
        std::string_view ref = extractField(p, lineEnd);
        if (p < lineEnd) ++p;
        std::string_view alt = extractField(p, lineEnd);
        if (p < lineEnd) ++p;

        // Build variant prefix once
        char *pp = prefixBuf;
        memcpy(pp, chrom.data(), chrom.size()); pp += chrom.size();
        *pp++ = '\t';
        memcpy(pp, pos.data(), pos.size()); pp += pos.size();
        *pp++ = '\t';
        memcpy(pp, id.data(), id.size()); pp += id.size();
        *pp++ = '\t';
        memcpy(pp, ref.data(), ref.size()); pp += ref.size();
        *pp++ = '\t';
        memcpy(pp, alt.data(), alt.size()); pp += alt.size();
        *pp++ = '\t';
        size_t prefixLen = pp - prefixBuf;

        // Skip QUAL, FILTER, INFO, FORMAT
        skipFields(p, lineEnd, 4);

        // Find all sample start positions
        findAllSampleStarts(p, lineEnd, sampleStarts);

        // Batch process: compute all allele balances
        for (size_t i = 0; i < numSamples; ++i) {
            size_t idx = sampleIndices[i];
            if (idx < sampleStarts.size()) {
                const char *gtStart = sampleStarts[idx];
                const char *nextTab = (idx + 1 < sampleStarts.size()) ?
                                      sampleStarts[idx + 1] - 1 : lineEnd;
                const char *gtEnd = findColon(gtStart, nextTab);
                if (gtEnd > nextTab) gtEnd = nextTab;

                balances[i] = computeAlleleBalanceFast(gtStart, gtEnd);
            } else {
                balances[i] = -1.0;  // Missing
            }
        }

        // Batch write: output all samples for this variant
        for (size_t i = 0; i < numSamples; ++i) {
            outBuf.write(prefixBuf, prefixLen);
            outBuf.write(sampleSuffix[i]);
            if (balances[i] < 0.0) {
                outBuf.write("NA", 2);
            } else {
                outBuf.writeDouble(balances[i]);
            }
            outBuf.writeChar('\n');
        }

        // Periodic flush to prevent memory buildup
        outBuf.maybeFlush();

        p = lineEnd;
        if (p < chunkEnd) ++p;
    }
}

// ---------------------------------------------------------------------
// Multi-threaded mmap mode
// ---------------------------------------------------------------------
static bool calculateBalanceMmapMT(const char *filename, const AlleleBalanceArgs &args) {
    MappedFile file;
    if (!file.open(filename)) {
        std::cerr << "Error: Cannot open file: " << filename << "\n";
        return false;
    }

    if (file.size == 0) {
        std::cerr << "Error: Empty file\n";
        return false;
    }

    const char *fileData = file.data;
    const char *fileEnd = file.data + file.size;
    const char *p = fileData;

    // Parse header
    std::vector<std::string_view> sampleNames;
    const char *dataStart = nullptr;

    while (p < fileEnd) {
        const char *lineEnd = findNewline(p, fileEnd);

        if (*p == '#') {
            if (lineEnd - p >= 6 && memcmp(p, "#CHROM", 6) == 0) {
                const char *hp = p;
                for (int i = 0; i < 9 && hp < lineEnd; ++i) {
                    hp = findTab(hp, lineEnd);
                    if (hp < lineEnd) ++hp;
                }
                while (hp < lineEnd) {
                    const char *nameStart = hp;
                    hp = findTab(hp, lineEnd);
                    sampleNames.emplace_back(nameStart, hp - nameStart);
                    if (hp < lineEnd) ++hp;
                }
            }
            p = lineEnd;
            if (p < fileEnd) ++p;
        } else {
            dataStart = p;
            break;
        }
    }

    if (sampleNames.empty()) {
        std::cerr << "Error: No samples found in VCF\n";
        return false;
    }

    if (!dataStart) {
        std::cerr << "Error: No data lines found\n";
        return false;
    }

    // Determine sample indices
    std::vector<size_t> sampleIndices;
    if (!args.samples.empty()) {
        std::unordered_map<std::string_view, size_t> sampleMap;
        for (size_t i = 0; i < sampleNames.size(); ++i) {
            sampleMap[sampleNames[i]] = i;
        }
        for (const auto &s : args.samples) {
            auto it = sampleMap.find(s);
            if (it == sampleMap.end()) {
                std::cerr << "Error: Sample '" << s << "' not found\n";
                return false;
            }
            sampleIndices.push_back(it->second);
        }
    } else {
        for (size_t i = 0; i < sampleNames.size(); ++i) {
            sampleIndices.push_back(i);
        }
    }

    // Build sample suffixes
    std::vector<std::string> sampleSuffix;
    for (size_t idx : sampleIndices) {
        std::string s(sampleNames[idx]);
        s.push_back('\t');
        sampleSuffix.push_back(std::move(s));
    }

    // Determine number of threads
    int numThreads = args.numThreads;
    if (numThreads <= 0) {
        numThreads = static_cast<int>(std::thread::hardware_concurrency());
        if (numThreads <= 0) numThreads = 4;
    }

    // For small files, use fewer threads
    size_t dataSize = fileEnd - dataStart;
    if (dataSize < 10 * 1024 * 1024) {
        numThreads = 1;
    } else if (dataSize < 100 * 1024 * 1024) {
        numThreads = std::min(numThreads, 4);
    }

    if (!args.quiet) {
        std::cerr << "Info: Using " << numThreads << " threads\n";
    }

    // Find line boundaries for chunking
    std::vector<const char *> chunkBoundaries;
    chunkBoundaries.push_back(dataStart);

    size_t chunkSize = dataSize / numThreads;
    for (int i = 1; i < numThreads; ++i) {
        const char *approxBoundary = dataStart + i * chunkSize;
        if (approxBoundary >= fileEnd) break;
        const char *lineEnd = findNewline(approxBoundary, fileEnd);
        if (lineEnd < fileEnd) ++lineEnd;
        if (lineEnd < fileEnd) {
            chunkBoundaries.push_back(lineEnd);
        }
    }
    chunkBoundaries.push_back(fileEnd);

    int actualThreads = static_cast<int>(chunkBoundaries.size()) - 1;

    // Write header
    const char *header = "CHROM\tPOS\tID\tREF\tALT\tSample\tAllele_Balance\n";
    ::write(STDOUT_FILENO, header, 43);

    // Process chunks sequentially to maintain output order and avoid memory buildup
    // Each chunk flushes incrementally, so memory usage stays bounded
    ThreadBuffer buf;
    for (int t = 0; t < actualThreads; ++t) {
        processChunk(chunkBoundaries[t], chunkBoundaries[t + 1],
                     sampleSuffix, sampleIndices, buf);
        buf.flush();  // Flush remaining data after each chunk
    }

    return true;
}

// ---------------------------------------------------------------------
// Single-threaded mmap mode
// ---------------------------------------------------------------------
static bool calculateBalanceMmapST(const char *filename, const AlleleBalanceArgs &args) {
    MappedFile file;
    if (!file.open(filename)) {
        std::cerr << "Error: Cannot open file: " << filename << "\n";
        return false;
    }

    if (file.size == 0) {
        std::cerr << "Error: Empty file\n";
        return false;
    }

    const char *p = file.data;
    const char *end = file.data + file.size;

    std::vector<std::string_view> sampleNames;
    bool foundHeader = false;

    while (p < end) {
        const char *lineEnd = findNewline(p, end);
        if (*p == '#') {
            if (lineEnd - p >= 6 && memcmp(p, "#CHROM", 6) == 0) {
                const char *hp = p;
                for (int i = 0; i < 9 && hp < lineEnd; ++i) {
                    hp = findTab(hp, lineEnd);
                    if (hp < lineEnd) ++hp;
                }
                while (hp < lineEnd) {
                    const char *nameStart = hp;
                    hp = findTab(hp, lineEnd);
                    sampleNames.emplace_back(nameStart, hp - nameStart);
                    if (hp < lineEnd) ++hp;
                }
                foundHeader = true;
            }
            p = lineEnd;
            if (p < end) ++p;
            continue;
        }
        break;
    }

    if (!foundHeader || sampleNames.empty()) {
        std::cerr << "Error: No samples found in VCF\n";
        return false;
    }

    std::vector<size_t> sampleIndices;
    if (!args.samples.empty()) {
        std::unordered_map<std::string_view, size_t> sampleMap;
        for (size_t i = 0; i < sampleNames.size(); ++i) {
            sampleMap[sampleNames[i]] = i;
        }
        for (const auto &s : args.samples) {
            auto it = sampleMap.find(s);
            if (it == sampleMap.end()) {
                std::cerr << "Error: Sample '" << s << "' not found\n";
                return false;
            }
            sampleIndices.push_back(it->second);
        }
    } else {
        for (size_t i = 0; i < sampleNames.size(); ++i) {
            sampleIndices.push_back(i);
        }
    }

    std::vector<std::string> sampleSuffix;
    for (size_t idx : sampleIndices) {
        std::string s(sampleNames[idx]);
        s.push_back('\t');
        sampleSuffix.push_back(std::move(s));
    }

    ThreadBuffer outBuf;

    const char *header = "CHROM\tPOS\tID\tREF\tALT\tSample\tAllele_Balance\n";
    outBuf.write(header, 43);

    char prefixBuf[4096];

    while (p < end) {
        const char *lineEnd = findNewline(p, end);
        if (p >= lineEnd || *p == '#') {
            p = lineEnd;
            if (p < end) ++p;
            continue;
        }

        std::string_view chrom = extractField(p, lineEnd);
        if (p < lineEnd) ++p;
        std::string_view pos = extractField(p, lineEnd);
        if (p < lineEnd) ++p;
        std::string_view id = extractField(p, lineEnd);
        if (p < lineEnd) ++p;
        std::string_view ref = extractField(p, lineEnd);
        if (p < lineEnd) ++p;
        std::string_view alt = extractField(p, lineEnd);
        if (p < lineEnd) ++p;

        char *pp = prefixBuf;
        memcpy(pp, chrom.data(), chrom.size()); pp += chrom.size();
        *pp++ = '\t';
        memcpy(pp, pos.data(), pos.size()); pp += pos.size();
        *pp++ = '\t';
        memcpy(pp, id.data(), id.size()); pp += id.size();
        *pp++ = '\t';
        memcpy(pp, ref.data(), ref.size()); pp += ref.size();
        *pp++ = '\t';
        memcpy(pp, alt.data(), alt.size()); pp += alt.size();
        *pp++ = '\t';
        size_t prefixLen = pp - prefixBuf;

        skipFields(p, lineEnd, 4);

        size_t sampleIdx = 0;
        const char *sampleStart = p;
        size_t suffixIdx = 0;

        for (size_t idx : sampleIndices) {
            while (sampleIdx < idx && sampleStart < lineEnd) {
                sampleStart = findTab(sampleStart, lineEnd);
                if (sampleStart < lineEnd) ++sampleStart;
                ++sampleIdx;
            }

            if (sampleStart >= lineEnd) break;

            const char *gtStart = sampleStart;
            const char *gtEnd = findColon(gtStart, lineEnd);
            const char *tabPos = findTab(gtStart, lineEnd);
            if (tabPos < gtEnd) gtEnd = tabPos;

            double balance = computeAlleleBalanceFast(gtStart, gtEnd);

            outBuf.write(prefixBuf, prefixLen);
            outBuf.write(sampleSuffix[suffixIdx]);
            if (balance < 0.0) {
                outBuf.write("NA", 2);
            } else {
                outBuf.writeDouble(balance);
            }
            outBuf.writeChar('\n');
            ++suffixIdx;
        }

        // Periodic flush to prevent memory buildup
        outBuf.maybeFlush();

        p = lineEnd;
        if (p < end) ++p;
    }

    // Flush remaining data
    outBuf.flush();

    return true;
}

// ---------------------------------------------------------------------
// Streaming stdin mode
// ---------------------------------------------------------------------
static bool calculateBalanceStream(std::istream &in, const AlleleBalanceArgs &args) {
    std::string line;
    line.reserve(65536);

    std::vector<std::string> sampleNames;
    std::vector<size_t> sampleIndices;
    bool foundHeader = false;

    ThreadBuffer outBuf;
    outBuf.write("CHROM\tPOS\tID\tREF\tALT\tSample\tAllele_Balance\n", 43);

    std::vector<std::string> sampleSuffix;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            if (line.size() >= 6 && memcmp(line.data(), "#CHROM", 6) == 0) {
                const char *p = line.data();
                const char *end = p + line.size();
                for (int i = 0; i < 9 && p < end; ++i) {
                    p = findTab(p, end);
                    if (p < end) ++p;
                }
                while (p < end) {
                    const char *nameStart = p;
                    p = findTab(p, end);
                    sampleNames.emplace_back(nameStart, p - nameStart);
                    if (p < end) ++p;
                }

                if (!args.samples.empty()) {
                    std::unordered_map<std::string, size_t> sampleMap;
                    for (size_t i = 0; i < sampleNames.size(); ++i) {
                        sampleMap[sampleNames[i]] = i;
                    }
                    for (const auto &s : args.samples) {
                        auto it = sampleMap.find(s);
                        if (it == sampleMap.end()) {
                            std::cerr << "Error: Sample '" << s << "' not found\n";
                            return false;
                        }
                        sampleIndices.push_back(it->second);
                    }
                } else {
                    for (size_t i = 0; i < sampleNames.size(); ++i) {
                        sampleIndices.push_back(i);
                    }
                }

                for (size_t idx : sampleIndices) {
                    std::string s = sampleNames[idx];
                    s.push_back('\t');
                    sampleSuffix.push_back(std::move(s));
                }

                foundHeader = true;
            }
            continue;
        }

        if (!foundHeader) {
            std::cerr << "Error: No #CHROM header found before data\n";
            return false;
        }

        const char *p = line.data();
        const char *end = p + line.size();

        std::string_view chrom = extractField(p, end);
        if (p < end) ++p;
        std::string_view pos = extractField(p, end);
        if (p < end) ++p;
        std::string_view id = extractField(p, end);
        if (p < end) ++p;
        std::string_view ref = extractField(p, end);
        if (p < end) ++p;
        std::string_view alt = extractField(p, end);
        if (p < end) ++p;

        char prefixBuf[4096];
        char *pp = prefixBuf;
        memcpy(pp, chrom.data(), chrom.size()); pp += chrom.size();
        *pp++ = '\t';
        memcpy(pp, pos.data(), pos.size()); pp += pos.size();
        *pp++ = '\t';
        memcpy(pp, id.data(), id.size()); pp += id.size();
        *pp++ = '\t';
        memcpy(pp, ref.data(), ref.size()); pp += ref.size();
        *pp++ = '\t';
        memcpy(pp, alt.data(), alt.size()); pp += alt.size();
        *pp++ = '\t';
        size_t prefixLen = pp - prefixBuf;

        skipFields(p, end, 4);

        size_t sampleIdx = 0;
        const char *sampleStart = p;
        size_t suffixIdx = 0;

        for (size_t idx : sampleIndices) {
            while (sampleIdx < idx && sampleStart < end) {
                sampleStart = findTab(sampleStart, end);
                if (sampleStart < end) ++sampleStart;
                ++sampleIdx;
            }

            if (sampleStart >= end) break;

            const char *gtStart = sampleStart;
            const char *gtEnd = findColon(gtStart, end);
            const char *tabPos = findTab(gtStart, end);
            if (tabPos < gtEnd) gtEnd = tabPos;

            double balance = computeAlleleBalanceFast(gtStart, gtEnd);

            outBuf.write(prefixBuf, prefixLen);
            outBuf.write(sampleSuffix[suffixIdx]);
            if (balance < 0.0) {
                outBuf.write("NA", 2);
            } else {
                outBuf.writeDouble(balance);
            }
            outBuf.writeChar('\n');
            ++suffixIdx;
        }

        // Periodic flush to prevent memory buildup
        outBuf.maybeFlush();
    }

    // Flush remaining data
    outBuf.flush();

    return foundHeader;
}

// ---------------------------------------------------------------------
// main
// ---------------------------------------------------------------------
int main(int argc, char *argv[]) {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);

    AlleleBalanceArgs args;
    if (!parseArguments(argc, argv, args)) {
        printHelp();
        return 0;
    }

    for (int i = 1; i < argc; ++i) {
        std::string_view arg = argv[i];
        if (arg == "--version" || arg == "-v") {
            std::cout << "VCFX_allele_balance_calc 2.0 (multi-threaded)\n";
            return 0;
        }
    }

    if (!args.quiet) {
        if (!args.samples.empty()) {
            std::cerr << "Info: Calculating allele balance for samples:";
            for (const auto &s : args.samples) {
                std::cerr << " " << s;
            }
            std::cerr << "\n";
        } else {
            std::cerr << "Info: Calculating allele balance for ALL samples\n";
        }
    }

    bool success;
    if (args.inputFile) {
        if (!args.quiet) {
            std::cerr << "Info: Using mmap mode for file: " << args.inputFile << "\n";
        }
        success = calculateBalanceMmapMT(args.inputFile, args);
    } else {
        if (!args.quiet) {
            std::cerr << "Info: Using stdin streaming mode (single-threaded)\n";
        }
        success = calculateBalanceStream(std::cin, args);
    }

    return success ? 0 : 1;
}
