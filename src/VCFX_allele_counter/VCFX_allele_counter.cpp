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
// Thread-local output buffer - accumulates output for later ordered write
// ---------------------------------------------------------------------
class ThreadBuffer {
  public:
    std::vector<char> buffer;
    size_t pos = 0;

    ThreadBuffer() {
        buffer.resize(16 * 1024 * 1024); // 16MB per thread
    }

    void ensureSpace(size_t needed) {
        if (pos + needed > buffer.size()) {
            buffer.resize(buffer.size() * 2);
        }
    }

    void write(const char *data, size_t len) {
        ensureSpace(len);
        memcpy(buffer.data() + pos, data, len);
        pos += len;
    }

    void write(std::string_view sv) { write(sv.data(), sv.size()); }

    void writeChar(char c) {
        ensureSpace(1);
        buffer[pos++] = c;
    }

    void writeInt(int val) {
        ensureSpace(12);
        if (val == 0) {
            buffer[pos++] = '0';
            return;
        }
        if (val < 0) {
            buffer[pos++] = '-';
            val = -val;
        }
        char tmp[12];
        int len = 0;
        while (val > 0) {
            tmp[len++] = '0' + (val % 10);
            val /= 10;
        }
        while (len > 0) {
            buffer[pos++] = tmp[--len];
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
// Ultra-fast genotype parsing
// ---------------------------------------------------------------------
static inline void parseGenotypeRaw(const char *gt, const char *gtEnd,
                                     int &refCount, int &altCount) {
    refCount = 0;
    altCount = 0;

    while (gt < gtEnd) {
        while (gt < gtEnd && (*gt == '/' || *gt == '|')) ++gt;
        if (gt >= gtEnd) break;

        if (*gt == '.') {
            ++gt;
            continue;
        }

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
}

static inline const char *findColon(const char *p, const char *end) {
    while (p < end && *p != ':') ++p;
    return p;
}

// ---------------------------------------------------------------------
// Arguments structure
// ---------------------------------------------------------------------
struct AlleleCounterArgs {
    std::vector<std::string> samples;
    const char *inputFile = nullptr;
    bool quiet = false;
    int numThreads = 0; // 0 = auto-detect
};

// ---------------------------------------------------------------------
// Parse arguments
// ---------------------------------------------------------------------
static bool parseArguments(int argc, char *argv[], AlleleCounterArgs &args) {
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
        << "VCFX_allele_counter - Count reference and alternate alleles per sample\n\n"
        << "Usage: VCFX_allele_counter [OPTIONS] [FILE]\n\n"
        << "Options:\n"
        << "  -i, --input FILE    Input VCF file (uses mmap for best performance)\n"
        << "  -t, --threads N     Number of threads (default: auto-detect CPU cores)\n"
        << "  -s, --samples STR   Space-separated list of sample names to include\n"
        << "  -q, --quiet         Suppress informational messages\n"
        << "  -h, --help          Display this help message\n"
        << "  -v, --version       Display version information\n\n"
        << "Examples:\n"
        << "  VCFX_allele_counter -i input.vcf > counts.tsv           # Auto threads\n"
        << "  VCFX_allele_counter -t 8 -i input.vcf > counts.tsv      # 8 threads\n"
        << "  VCFX_allele_counter < input.vcf > counts.tsv            # Stdin (single-thread)\n\n"
        << "Output format:\n"
        << "  CHROM  POS  ID  REF  ALT  Sample  Ref_Count  Alt_Count\n";
}

// ---------------------------------------------------------------------
// Find all sample start positions in a line (called once per line)
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
// Process a chunk of lines (called by each thread)
// Uses batch processing for all samples per variant
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

    // Pre-allocated vectors to avoid per-line allocations
    std::vector<const char*> sampleStarts;
    sampleStarts.reserve(3000); // Most VCFs have <3000 samples

    // Pre-allocated result arrays
    const size_t numSamples = sampleIndices.size();
    std::vector<int8_t> refCounts(numSamples);
    std::vector<int8_t> altCounts(numSamples);

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

        // Find all sample start positions once
        findAllSampleStarts(p, lineEnd, sampleStarts);

        // Batch process: parse all genotypes first
        for (size_t i = 0; i < numSamples; ++i) {
            size_t idx = sampleIndices[i];
            if (idx < sampleStarts.size()) {
                const char *gtStart = sampleStarts[idx];
                const char *nextTab = (idx + 1 < sampleStarts.size()) ?
                                      sampleStarts[idx + 1] - 1 : lineEnd;
                const char *gtEnd = findColon(gtStart, nextTab);
                if (gtEnd > nextTab) gtEnd = nextTab;

                int refC = 0, altC = 0;
                parseGenotypeRaw(gtStart, gtEnd, refC, altC);
                refCounts[i] = static_cast<int8_t>(refC);
                altCounts[i] = static_cast<int8_t>(altC);
            } else {
                refCounts[i] = 0;
                altCounts[i] = 0;
            }
        }

        // Batch write: output all samples for this variant
        for (size_t i = 0; i < numSamples; ++i) {
            outBuf.write(prefixBuf, prefixLen);
            outBuf.write(sampleSuffix[i]);
            outBuf.writeInt(refCounts[i]);
            outBuf.writeChar('\t');
            outBuf.writeInt(altCounts[i]);
            outBuf.writeChar('\n');
        }

        p = lineEnd;
        if (p < chunkEnd) ++p;
    }
}

// ---------------------------------------------------------------------
// Process chunk and write to file descriptor (for streaming output)
// ---------------------------------------------------------------------
static void processChunkToFd(
    const char *chunkStart,
    const char *chunkEnd,
    const std::vector<std::string> &sampleSuffix,
    const std::vector<size_t> &sampleIndices,
    int outFd)
{
    // Use a moderate buffer that gets flushed frequently
    static constexpr size_t BUFFER_SIZE = 4 * 1024 * 1024; // 4MB
    std::vector<char> buffer(BUFFER_SIZE);
    size_t bufPos = 0;

    auto flushBuffer = [&]() {
        if (bufPos > 0) {
            size_t written = 0;
            while (written < bufPos) {
                ssize_t n = ::write(outFd, buffer.data() + written, bufPos - written);
                if (n < 0) break;
                written += n;
            }
            bufPos = 0;
        }
    };

    auto ensureSpace = [&](size_t needed) {
        if (bufPos + needed > BUFFER_SIZE) {
            flushBuffer();
        }
    };

    auto writeData = [&](const char *data, size_t len) {
        ensureSpace(len);
        memcpy(buffer.data() + bufPos, data, len);
        bufPos += len;
    };

    auto writeChar = [&](char c) {
        ensureSpace(1);
        buffer[bufPos++] = c;
    };

    auto writeInt = [&](int val) {
        ensureSpace(12);
        if (val == 0) {
            buffer[bufPos++] = '0';
            return;
        }
        char tmp[12];
        int len = 0;
        while (val > 0) {
            tmp[len++] = '0' + (val % 10);
            val /= 10;
        }
        while (len > 0) {
            buffer[bufPos++] = tmp[--len];
        }
    };

    char prefixBuf[4096];
    const char *p = chunkStart;

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

        // Build variant prefix
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

        // Process samples
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

            int refCount = 0, altCount = 0;
            parseGenotypeRaw(gtStart, gtEnd, refCount, altCount);

            writeData(prefixBuf, prefixLen);
            writeData(sampleSuffix[suffixIdx].data(), sampleSuffix[suffixIdx].size());
            writeInt(refCount);
            writeChar('\t');
            writeInt(altCount);
            writeChar('\n');
            ++suffixIdx;
        }

        p = lineEnd;
        if (p < chunkEnd) ++p;
    }

    flushBuffer();
}

// ---------------------------------------------------------------------
// Multi-threaded memory-mapped allele counting with in-memory buffers
// ---------------------------------------------------------------------
static bool countAllelesMmapMT(const char *filename, const AlleleCounterArgs &args) {
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
    if (dataSize < 10 * 1024 * 1024) { // < 10MB
        numThreads = 1;
    } else if (dataSize < 100 * 1024 * 1024) { // < 100MB
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
    const char *header = "CHROM\tPOS\tID\tREF\tALT\tSample\tRef_Count\tAlt_Count\n";
    ::write(STDOUT_FILENO, header, 48);

    if (actualThreads == 1) {
        // Single thread - use in-memory buffer
        ThreadBuffer buf;
        processChunk(chunkBoundaries[0], chunkBoundaries[1],
                     sampleSuffix, sampleIndices, buf);
        if (buf.size() > 0) {
            ::write(STDOUT_FILENO, buf.data(), buf.size());
        }
    } else {
        // Multi-threaded with in-memory buffers
        std::vector<ThreadBuffer> buffers(actualThreads);
        std::vector<std::thread> threads;

        for (int t = 0; t < actualThreads; ++t) {
            threads.emplace_back([&, t]() {
                processChunk(
                    chunkBoundaries[t],
                    chunkBoundaries[t + 1],
                    sampleSuffix,
                    sampleIndices,
                    buffers[t]
                );
            });
        }

        // Wait for all threads
        for (auto &th : threads) {
            th.join();
        }

        // Write buffers in order
        for (int t = 0; t < actualThreads; ++t) {
            if (buffers[t].size() > 0) {
                ::write(STDOUT_FILENO, buffers[t].data(), buffers[t].size());
            }
        }
    }

    return true;
}

// ---------------------------------------------------------------------
// Single-threaded mmap (for comparison/fallback)
// ---------------------------------------------------------------------
static bool countAllelesMmapST(const char *filename, const AlleleCounterArgs &args) {
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

    // Use ThreadBuffer for single-thread output
    ThreadBuffer outBuf;

    const char *header = "CHROM\tPOS\tID\tREF\tALT\tSample\tRef_Count\tAlt_Count\n";
    outBuf.write(header, 48);

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

            int refCount = 0, altCount = 0;
            parseGenotypeRaw(gtStart, gtEnd, refCount, altCount);

            outBuf.write(prefixBuf, prefixLen);
            outBuf.write(sampleSuffix[suffixIdx]);
            outBuf.writeInt(refCount);
            outBuf.writeChar('\t');
            outBuf.writeInt(altCount);
            outBuf.writeChar('\n');
            ++suffixIdx;

            // Flush periodically to avoid excessive memory
            if (outBuf.size() > 64 * 1024 * 1024) {
                ::write(STDOUT_FILENO, outBuf.data(), outBuf.size());
                outBuf.clear();
            }
        }

        p = lineEnd;
        if (p < end) ++p;
    }

    // Final flush
    if (outBuf.size() > 0) {
        ::write(STDOUT_FILENO, outBuf.data(), outBuf.size());
    }

    return true;
}

// ---------------------------------------------------------------------
// Streaming stdin mode (single-threaded)
// ---------------------------------------------------------------------
static bool countAllelesStream(std::istream &in, const AlleleCounterArgs &args) {
    std::string line;
    line.reserve(65536);

    std::vector<std::string> sampleNames;
    std::vector<size_t> sampleIndices;
    bool foundHeader = false;

    ThreadBuffer outBuf;
    outBuf.write("CHROM\tPOS\tID\tREF\tALT\tSample\tRef_Count\tAlt_Count\n", 48);

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

            int refCount = 0, altCount = 0;
            parseGenotypeRaw(gtStart, gtEnd, refCount, altCount);

            outBuf.write(prefixBuf, prefixLen);
            outBuf.write(sampleSuffix[suffixIdx]);
            outBuf.writeInt(refCount);
            outBuf.writeChar('\t');
            outBuf.writeInt(altCount);
            outBuf.writeChar('\n');
            ++suffixIdx;
        }

        // Flush periodically
        if (outBuf.size() > 64 * 1024 * 1024) {
            ::write(STDOUT_FILENO, outBuf.data(), outBuf.size());
            outBuf.clear();
        }
    }

    if (outBuf.size() > 0) {
        ::write(STDOUT_FILENO, outBuf.data(), outBuf.size());
    }

    return foundHeader;
}

// ---------------------------------------------------------------------
// main
// ---------------------------------------------------------------------
int main(int argc, char *argv[]) {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);

    AlleleCounterArgs args;
    if (!parseArguments(argc, argv, args)) {
        printHelp();
        return 0;
    }

    for (int i = 1; i < argc; ++i) {
        std::string_view arg = argv[i];
        if (arg == "--version" || arg == "-v") {
            std::cout << "VCFX_allele_counter 2.0 (multi-threaded)\n";
            return 0;
        }
    }

    if (!args.quiet) {
        if (!args.samples.empty()) {
            std::cerr << "Info: Counting alleles for samples:";
            for (const auto &s : args.samples) {
                std::cerr << " " << s;
            }
            std::cerr << "\n";
        } else {
            std::cerr << "Info: Counting alleles for ALL samples\n";
        }
    }

    bool success;
    if (args.inputFile) {
        if (!args.quiet) {
            std::cerr << "Info: Using mmap mode for file: " << args.inputFile << "\n";
        }
        // Use multi-threaded version
        success = countAllelesMmapMT(args.inputFile, args);
    } else {
        if (!args.quiet) {
            std::cerr << "Info: Using stdin streaming mode (single-threaded)\n";
        }
        success = countAllelesStream(std::cin, args);
    }

    return success ? 0 : 1;
}
