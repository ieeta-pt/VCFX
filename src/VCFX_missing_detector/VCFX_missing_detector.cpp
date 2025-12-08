/**
 * VCFX_missing_detector - Extreme-performance missing genotype detector
 *
 * ALGORITHM INSIGHT:
 * Instead of parsing each sample, we use SIMD to scan for '.' characters
 * in the sample columns region. Most lines have NO missing data, so we
 * can just copy them verbatim without any parsing!
 *
 * Optimizations:
 * - Memory-mapped I/O with MADV_SEQUENTIAL | MADV_WILLNEED
 * - SIMD '.' character search (AVX2/SSE2/NEON) - 16-32 bytes at a time
 * - Multi-threaded chunk processing
 * - Zero-copy output for lines without missing genotypes
 * - Early exit on first '.' found in sample region
 * - 4MB output buffer with direct write() syscalls
 */

#include "VCFX_missing_detector.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <string_view>
#include <iostream>
#include <sstream>
#include <vector>
#include <thread>
#include <atomic>
#include <mutex>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <getopt.h>

// ============================================================================
// SIMD Detection and Intrinsics
// ============================================================================
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
    #if defined(__AVX2__)
        #include <immintrin.h>
        #define USE_AVX2 1
    #elif defined(__SSE2__)
        #include <emmintrin.h>
        #define USE_SSE2 1
    #endif
#elif defined(__aarch64__) || defined(__ARM_NEON)
    #include <arm_neon.h>
    #define USE_NEON 1
#endif

// ============================================================================
// Memory-mapped file structure
// ============================================================================
struct MappedFile {
    const char *data = nullptr;
    size_t size = 0;
    int fd = -1;

    bool open(const char *path) {
        fd = ::open(path, O_RDONLY);
        if (fd < 0) return false;
        struct stat st;
        if (fstat(fd, &st) < 0) { ::close(fd); fd = -1; return false; }
        size = static_cast<size_t>(st.st_size);
        if (size == 0) { ::close(fd); fd = -1; return true; }
        data = static_cast<const char*>(mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0));
        if (data == MAP_FAILED) { data = nullptr; ::close(fd); fd = -1; return false; }
        madvise(const_cast<char*>(data), size, MADV_SEQUENTIAL | MADV_WILLNEED);
        return true;
    }

    void close() {
        if (data && size > 0) munmap(const_cast<char*>(data), size);
        if (fd >= 0) ::close(fd);
        data = nullptr; size = 0; fd = -1;
    }
};

// ============================================================================
// Output buffer for efficient writes
// ============================================================================
class OutputBuffer {
    static constexpr size_t CAPACITY = 4 * 1024 * 1024; // 4MB buffer
    char *buf;
    size_t pos = 0;
    int fd;
public:
    explicit OutputBuffer(int fd_) : fd(fd_) { buf = new char[CAPACITY]; }
    ~OutputBuffer() { flush(); delete[] buf; }

    void flush() {
        if (pos > 0) { ::write(fd, buf, pos); pos = 0; }
    }

    void ensureSpace(size_t n) {
        if (pos + n > CAPACITY) flush();
    }

    void append(char c) {
        if (pos >= CAPACITY) flush();
        buf[pos++] = c;
    }

    void append(const char *s, size_t len) {
        while (len > 0) {
            size_t chunk = std::min(len, CAPACITY - pos);
            memcpy(buf + pos, s, chunk);
            pos += chunk; s += chunk; len -= chunk;
            if (pos >= CAPACITY) flush();
        }
    }

    void append(std::string_view sv) { append(sv.data(), sv.size()); }
};

// ============================================================================
// SIMD-accelerated character search functions
// ============================================================================

// Find newline
#if defined(USE_AVX2)
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

static inline const char* findTabSIMD(const char* p, const char* end) {
    const __m256i tab = _mm256_set1_epi8('\t');
    while (p + 32 <= end) {
        __m256i chunk = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p));
        __m256i cmp = _mm256_cmpeq_epi8(chunk, tab);
        int mask = _mm256_movemask_epi8(cmp);
        if (mask) return p + __builtin_ctz(mask);
        p += 32;
    }
    while (p < end && *p != '\t' && *p != '\n') ++p;
    return p;
}

// CRITICAL: SIMD search for '.' in sample region - this is the key optimization!
static inline bool hasDotInRegionSIMD(const char* p, const char* end) {
    const __m256i dot = _mm256_set1_epi8('.');
    while (p + 32 <= end) {
        __m256i chunk = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p));
        __m256i cmp = _mm256_cmpeq_epi8(chunk, dot);
        int mask = _mm256_movemask_epi8(cmp);
        if (mask) return true;
        p += 32;
    }
    while (p < end) {
        if (*p == '.') return true;
        ++p;
    }
    return false;
}

#elif defined(USE_SSE2)
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

static inline const char* findTabSIMD(const char* p, const char* end) {
    const __m128i tab = _mm_set1_epi8('\t');
    while (p + 16 <= end) {
        __m128i chunk = _mm_loadu_si128(reinterpret_cast<const __m128i*>(p));
        __m128i cmp = _mm_cmpeq_epi8(chunk, tab);
        int mask = _mm_movemask_epi8(cmp);
        if (mask) return p + __builtin_ctz(mask);
        p += 16;
    }
    while (p < end && *p != '\t' && *p != '\n') ++p;
    return p;
}

static inline bool hasDotInRegionSIMD(const char* p, const char* end) {
    const __m128i dot = _mm_set1_epi8('.');
    while (p + 16 <= end) {
        __m128i chunk = _mm_loadu_si128(reinterpret_cast<const __m128i*>(p));
        __m128i cmp = _mm_cmpeq_epi8(chunk, dot);
        int mask = _mm_movemask_epi8(cmp);
        if (mask) return true;
        p += 16;
    }
    while (p < end) {
        if (*p == '.') return true;
        ++p;
    }
    return false;
}

#elif defined(USE_NEON)
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

static inline const char* findTabSIMD(const char* p, const char* end) {
    const uint8x16_t tab = vdupq_n_u8('\t');
    while (p + 16 <= end) {
        uint8x16_t chunk = vld1q_u8(reinterpret_cast<const uint8_t*>(p));
        uint8x16_t cmp = vceqq_u8(chunk, tab);
        uint64_t lo = vgetq_lane_u64(vreinterpretq_u64_u8(cmp), 0);
        uint64_t hi = vgetq_lane_u64(vreinterpretq_u64_u8(cmp), 1);
        if (lo) return p + (__builtin_ctzll(lo) >> 3);
        if (hi) return p + 8 + (__builtin_ctzll(hi) >> 3);
        p += 16;
    }
    while (p < end && *p != '\t' && *p != '\n') ++p;
    return p;
}

static inline bool hasDotInRegionSIMD(const char* p, const char* end) {
    const uint8x16_t dot = vdupq_n_u8('.');
    while (p + 16 <= end) {
        uint8x16_t chunk = vld1q_u8(reinterpret_cast<const uint8_t*>(p));
        uint8x16_t cmp = vceqq_u8(chunk, dot);
        uint64_t lo = vgetq_lane_u64(vreinterpretq_u64_u8(cmp), 0);
        uint64_t hi = vgetq_lane_u64(vreinterpretq_u64_u8(cmp), 1);
        if (lo || hi) return true;
        p += 16;
    }
    while (p < end) {
        if (*p == '.') return true;
        ++p;
    }
    return false;
}

#else
// Fallback using memchr
static inline const char* findNewlineSIMD(const char* p, const char* end) {
    const char* found = static_cast<const char*>(memchr(p, '\n', end - p));
    return found ? found : end;
}

static inline const char* findTabSIMD(const char* p, const char* end) {
    while (p < end && *p != '\t' && *p != '\n') ++p;
    return p;
}

static inline bool hasDotInRegionSIMD(const char* p, const char* end) {
    return memchr(p, '.', end - p) != nullptr;
}
#endif

// ============================================================================
// Skip to field N (0-indexed) - returns pointer to start of that field
// ============================================================================
static inline const char* skipToField(const char* p, const char* end, int fieldIdx) {
    for (int i = 0; i < fieldIdx && p < end; ++i) {
        p = findTabSIMD(p, end);
        if (p < end) ++p;
    }
    return p;
}

// ============================================================================
// Check if '.' in sample region represents a missing genotype
// We need to verify it's not in INFO field values like AF=0.5
// Strategy: Check if '.' is at position 0 of a GT field or followed by / or |
// ============================================================================
static inline bool hasMissingGenotypeInSamples(const char* sampleStart, const char* lineEnd) {
    // Quick check: if no '.' at all, definitely no missing
    if (!hasDotInRegionSIMD(sampleStart, lineEnd)) {
        return false;
    }

    // There's a '.' somewhere - need to check if it's in a genotype
    // Scan through samples looking for genotype patterns with '.'
    const char* p = sampleStart;

    while (p < lineEnd) {
        // Find start of this sample's GT (first field before :)
        const char* sampleEnd = findTabSIMD(p, lineEnd);

        // Look at the genotype portion (before first ':')
        const char* gtEnd = p;
        while (gtEnd < sampleEnd && *gtEnd != ':') ++gtEnd;

        // Check for missing genotype patterns in this GT
        // Missing patterns: ".", "./.", ".|.", "./X", "X/.", ".|X", "X|."
        const char* g = p;
        while (g < gtEnd) {
            if (*g == '.') {
                // Check context: is this a missing allele?
                // It's missing if:
                // 1. It's the entire genotype (".")
                // 2. It's followed by / or | ("./", ".|")
                // 3. It's preceded by / or | ("/.", "|.")
                // 4. It's at the start and followed by / or |

                bool prevIsSep = (g == p) || (*(g-1) == '/' || *(g-1) == '|');
                bool nextIsSep = (g+1 >= gtEnd) || (*(g+1) == '/' || *(g+1) == '|');

                if (prevIsSep || nextIsSep) {
                    return true;
                }
            }
            ++g;
        }

        // Move to next sample
        if (sampleEnd >= lineEnd) break;
        p = sampleEnd + 1;
    }

    return false;
}

// ============================================================================
// EXTREME optimization: Multi-threaded pre-scan of SAMPLE COLUMNS for '.'
// If no dots in sample columns of any line, file can be copied verbatim!
// ============================================================================
struct ScanChunkResult {
    bool hasDot = false;
    size_t lineCount = 0;
};

static void scanChunkForDots(const char* start, const char* end, ScanChunkResult &result) {
    result.hasDot = false;
    result.lineCount = 0;

    const char* p = start;
    while (p < end) {
        const char* lineEnd = findNewlineSIMD(p, end);
        if (lineEnd >= end) break;  // Partial line - skip

        result.lineCount++;

        // Skip to field 9 (samples start) - use SIMD-accelerated tab scanning
        const char* sampleStart = skipToField(p, lineEnd, 9);
        if (sampleStart < lineEnd) {
            if (hasDotInRegionSIMD(sampleStart, lineEnd)) {
                result.hasDot = true;
                return;  // Early exit!
            }
        }

        p = lineEnd + 1;
    }
}

static bool sampleColumnsHaveAnyDots(const MappedFile &mf, const char* &dataStart,
                                      size_t &lineCount, int numThreads) {
    const char* p = mf.data;
    const char* end = mf.data + mf.size;
    dataStart = p;
    lineCount = 0;

    // Skip header lines
    while (p < end && *p == '#') {
        p = findNewlineSIMD(p, end);
        if (p < end) ++p;
    }
    dataStart = p;

    // Single-threaded for small files
    size_t dataSize = end - p;
    if (dataSize < 10 * 1024 * 1024 || numThreads <= 1) {
        ScanChunkResult result;
        scanChunkForDots(p, end, result);
        lineCount = result.lineCount;
        return result.hasDot;
    }

    // Multi-threaded scan
    size_t chunkSize = dataSize / numThreads;
    if (chunkSize < 1024 * 1024) chunkSize = 1024 * 1024;

    std::vector<std::thread> threads;
    std::vector<ScanChunkResult> results(numThreads);
    std::vector<const char*> chunkStarts(numThreads);
    std::vector<const char*> chunkEnds(numThreads);
    std::atomic<bool> foundDot{false};

    // Calculate chunk boundaries (align to newlines)
    const char* chunkStart = p;
    for (int i = 0; i < numThreads; ++i) {
        chunkStarts[i] = chunkStart;

        if (i == numThreads - 1) {
            chunkEnds[i] = end;
        } else {
            const char* tentativeEnd = chunkStart + chunkSize;
            if (tentativeEnd >= end) {
                tentativeEnd = end;
            } else {
                tentativeEnd = findNewlineSIMD(tentativeEnd, end);
                if (tentativeEnd < end) ++tentativeEnd;
            }
            chunkEnds[i] = tentativeEnd;
            chunkStart = tentativeEnd;
        }
    }

    // Launch threads
    for (int i = 0; i < numThreads; ++i) {
        threads.emplace_back([i, &chunkStarts, &chunkEnds, &results, &foundDot]() {
            if (foundDot.load(std::memory_order_relaxed)) return;  // Early exit
            scanChunkForDots(chunkStarts[i], chunkEnds[i], results[i]);
            if (results[i].hasDot) {
                foundDot.store(true, std::memory_order_relaxed);
            }
        });
    }

    // Wait and collect
    for (int i = 0; i < numThreads; ++i) {
        threads[i].join();
        lineCount += results[i].lineCount;
        if (results[i].hasDot) {
            return true;
        }
    }

    return false;
}

// ============================================================================
// Ultra-fast zero-copy output: batch consecutive unchanged lines
// ============================================================================
static void processMmapZeroCopy(const MappedFile &mf, int outFd, bool quiet, int numThreads) {
    const char* dataStart = nullptr;
    size_t lineCount = 0;

    // FAST PATH: Multi-threaded pre-scan sample columns for ANY '.' characters
    // If no dots in sample columns, file can be copied verbatim!
    if (!sampleColumnsHaveAnyDots(mf, dataStart, lineCount, numThreads)) {
        // NO DOTS IN SAMPLE COLUMNS = NO MISSING GENOTYPES
        // Just write the entire file verbatim!
        if (!quiet) {
            std::cerr << "Fast path: No '.' in sample columns (scan complete)\n";
        }

        // Write entire file using a single large write for maximum throughput
        const char* p = mf.data;
        const char* end = mf.data + mf.size;
        constexpr size_t CHUNK = 16 * 1024 * 1024;  // 16MB chunks
        while (p < end) {
            size_t len = std::min(static_cast<size_t>(end - p), CHUNK);
            ::write(outFd, p, len);
            p += len;
        }

        if (!quiet) {
            std::cerr << "Processed " << lineCount << " variants, 0 with missing genotypes (0%)\n";
        }
        return;
    }

    // SLOW PATH: File has dots, need to check each line
    const char* p = mf.data;
    const char* end = mf.data + mf.size;

    size_t totalLines = 0;
    size_t missingLines = 0;

    // Track start of unchanged region for batched writes
    const char* batchStart = nullptr;

    // Small buffer for modified lines
    char modBuf[65536];

    auto flushBatch = [&]() {
        if (batchStart && batchStart < p) {
            ::write(outFd, batchStart, p - batchStart);
            batchStart = nullptr;
        }
    };

    while (p < end) {
        const char* lineStart = p;
        const char* lineEnd = findNewlineSIMD(p, end);
        const char* nextLine = (lineEnd < end) ? lineEnd + 1 : end;

        // Handle \r\n
        const char* adjLineEnd = lineEnd;
        if (adjLineEnd > lineStart && *(adjLineEnd - 1) == '\r') --adjLineEnd;

        size_t lineLen = adjLineEnd - lineStart;

        // Empty line or header - include in batch
        if (lineLen == 0 || *lineStart == '#') {
            if (!batchStart) batchStart = lineStart;
            p = nextLine;
            continue;
        }

        totalLines++;

        // Data line - check for missing genotypes
        const char* sampleStart = skipToField(lineStart, adjLineEnd, 9);

        if (sampleStart >= adjLineEnd) {
            // No sample columns - include in batch
            if (!batchStart) batchStart = lineStart;
            p = nextLine;
            continue;
        }

        // FAST PATH: Use SIMD to check for '.' in sample region
        bool hasMissing = hasMissingGenotypeInSamples(sampleStart, adjLineEnd);

        if (!hasMissing) {
            // No missing - include in batch (ZERO COPY!)
            if (!batchStart) batchStart = lineStart;
            p = nextLine;
            continue;
        }

        // Has missing - need to modify this line
        // First, flush any pending batch
        flushBatch();

        missingLines++;

        // Find INFO field (field 7)
        const char* infoStart = skipToField(lineStart, adjLineEnd, 7);
        const char* infoEnd = findTabSIMD(infoStart, adjLineEnd);
        std::string_view info(infoStart, infoEnd - infoStart);

        // Build modified line in buffer
        char* wp = modBuf;

        // Copy before INFO
        size_t beforeLen = infoStart - lineStart;
        memcpy(wp, lineStart, beforeLen);
        wp += beforeLen;

        // Write modified INFO
        if (info == "." || info.empty()) {
            memcpy(wp, "MISSING_GENOTYPES=1", 19);
            wp += 19;
        } else {
            memcpy(wp, info.data(), info.size());
            wp += info.size();
            if (info.back() != ';') *wp++ = ';';
            memcpy(wp, "MISSING_GENOTYPES=1", 19);
            wp += 19;
        }

        // Copy after INFO
        size_t afterLen = adjLineEnd - infoEnd;
        memcpy(wp, infoEnd, afterLen);
        wp += afterLen;
        *wp++ = '\n';

        ::write(outFd, modBuf, wp - modBuf);

        p = nextLine;
    }

    flushBatch();

    if (!quiet) {
        std::cerr << "Processed " << totalLines << " variants, "
                  << missingLines << " with missing genotypes ("
                  << (totalLines > 0 ? (100.0 * missingLines / totalLines) : 0.0)
                  << "%)\n";
    }
}

// ============================================================================
// Process mmap'd file - single threaded optimized version
// ============================================================================
static void processMmap(const MappedFile &mf, OutputBuffer &out, bool quiet) {
    const char* p = mf.data;
    const char* end = mf.data + mf.size;

    size_t totalLines = 0;
    size_t missingLines = 0;

    while (p < end) {
        const char* lineStart = p;
        const char* lineEnd = findNewlineSIMD(p, end);

        // Handle \r\n
        const char* adjLineEnd = lineEnd;
        if (adjLineEnd > lineStart && *(adjLineEnd - 1) == '\r') --adjLineEnd;

        size_t lineLen = adjLineEnd - lineStart;

        // Empty line
        if (lineLen == 0) {
            out.append('\n');
            p = lineEnd + 1;
            continue;
        }

        // Header lines - pass through unchanged
        if (*lineStart == '#') {
            out.append(lineStart, lineLen);
            out.append('\n');
            p = lineEnd + 1;
            continue;
        }

        totalLines++;

        // Data line - check for missing genotypes
        // Skip to field 9 (first sample column, 0-indexed)
        const char* sampleStart = skipToField(lineStart, adjLineEnd, 9);

        if (sampleStart >= adjLineEnd) {
            // No sample columns - pass through
            out.append(lineStart, lineLen);
            out.append('\n');
            p = lineEnd + 1;
            continue;
        }

        // FAST PATH: Check if any '.' in sample region
        bool hasMissing = hasMissingGenotypeInSamples(sampleStart, adjLineEnd);

        if (!hasMissing) {
            // No missing genotypes - copy line verbatim (FAST PATH)
            out.append(lineStart, lineLen);
            out.append('\n');
        } else {
            // Has missing - need to modify INFO field
            missingLines++;

            // Find INFO field (field 7)
            const char* infoStart = skipToField(lineStart, adjLineEnd, 7);
            const char* infoEnd = findTabSIMD(infoStart, adjLineEnd);

            std::string_view infoBefore(lineStart, infoStart - lineStart);
            std::string_view info(infoStart, infoEnd - infoStart);
            std::string_view infoAfter(infoEnd, adjLineEnd - infoEnd);

            // Write: before INFO + modified INFO + after INFO
            out.append(infoBefore);

            if (info == "." || info.empty()) {
                out.append("MISSING_GENOTYPES=1", 19);
            } else {
                out.append(info);
                if (info.back() != ';') {
                    out.append(';');
                }
                out.append("MISSING_GENOTYPES=1", 19);
            }

            out.append(infoAfter);
            out.append('\n');
        }

        p = lineEnd + 1;
    }

    if (!quiet) {
        std::cerr << "Processed " << totalLines << " variants, "
                  << missingLines << " with missing genotypes ("
                  << (totalLines > 0 ? (100.0 * missingLines / totalLines) : 0.0)
                  << "%)\n";
    }
}

// ============================================================================
// Multi-threaded chunk processing
// ============================================================================
struct ChunkResult {
    std::string output;
    size_t totalLines = 0;
    size_t missingLines = 0;
};

static void processChunk(const char* chunkStart, const char* chunkEnd,
                         bool isFirstChunk, ChunkResult &result) {
    result.output.reserve(chunkEnd - chunkStart + 1024);

    const char* p = chunkStart;

    while (p < chunkEnd) {
        const char* lineStart = p;
        const char* lineEnd = findNewlineSIMD(p, chunkEnd);

        // Handle partial line at end of chunk (except for last chunk)
        if (lineEnd >= chunkEnd) {
            // This line extends beyond our chunk - skip it (will be processed by next chunk)
            break;
        }

        // Handle \r\n
        const char* adjLineEnd = lineEnd;
        if (adjLineEnd > lineStart && *(adjLineEnd - 1) == '\r') --adjLineEnd;

        size_t lineLen = adjLineEnd - lineStart;

        if (lineLen == 0) {
            result.output += '\n';
            p = lineEnd + 1;
            continue;
        }

        // Header lines - only in first chunk
        if (*lineStart == '#') {
            if (isFirstChunk) {
                result.output.append(lineStart, lineLen);
                result.output += '\n';
            }
            p = lineEnd + 1;
            continue;
        }

        result.totalLines++;

        // Skip to samples (field 9)
        const char* sampleStart = skipToField(lineStart, adjLineEnd, 9);

        if (sampleStart >= adjLineEnd) {
            result.output.append(lineStart, lineLen);
            result.output += '\n';
            p = lineEnd + 1;
            continue;
        }

        bool hasMissing = hasMissingGenotypeInSamples(sampleStart, adjLineEnd);

        if (!hasMissing) {
            result.output.append(lineStart, lineLen);
            result.output += '\n';
        } else {
            result.missingLines++;

            const char* infoStart = skipToField(lineStart, adjLineEnd, 7);
            const char* infoEnd = findTabSIMD(infoStart, adjLineEnd);

            result.output.append(lineStart, infoStart - lineStart);

            std::string_view info(infoStart, infoEnd - infoStart);
            if (info == "." || info.empty()) {
                result.output.append("MISSING_GENOTYPES=1");
            } else {
                result.output.append(info);
                if (info.back() != ';') result.output += ';';
                result.output.append("MISSING_GENOTYPES=1");
            }

            result.output.append(infoEnd, adjLineEnd - infoEnd);
            result.output += '\n';
        }

        p = lineEnd + 1;
    }
}

static void processMmapMultithreaded(const MappedFile &mf, OutputBuffer &out,
                                     int numThreads, bool quiet) {
    if (mf.size < 10 * 1024 * 1024 || numThreads <= 1) {
        // Small file or single thread - use simple version
        processMmap(mf, out, quiet);
        return;
    }

    // Find where header ends
    const char* p = mf.data;
    const char* end = mf.data + mf.size;
    const char* dataStart = p;

    while (p < end) {
        if (*p != '#') {
            dataStart = p;
            break;
        }
        p = findNewlineSIMD(p, end);
        if (p < end) ++p;
    }

    // Output headers first
    if (dataStart > mf.data) {
        out.append(mf.data, dataStart - mf.data);
    }

    size_t dataSize = end - dataStart;
    size_t chunkSize = dataSize / numThreads;
    if (chunkSize < 1024 * 1024) chunkSize = 1024 * 1024; // Min 1MB chunks

    std::vector<std::thread> threads;
    std::vector<ChunkResult> results(numThreads);
    std::vector<const char*> chunkStarts(numThreads);
    std::vector<const char*> chunkEnds(numThreads);

    // Calculate chunk boundaries (align to newlines)
    const char* chunkStart = dataStart;
    for (int i = 0; i < numThreads; ++i) {
        chunkStarts[i] = chunkStart;

        if (i == numThreads - 1) {
            chunkEnds[i] = end;
        } else {
            const char* tentativeEnd = chunkStart + chunkSize;
            if (tentativeEnd >= end) {
                tentativeEnd = end;
            } else {
                // Align to next newline
                tentativeEnd = findNewlineSIMD(tentativeEnd, end);
                if (tentativeEnd < end) ++tentativeEnd;
            }
            chunkEnds[i] = tentativeEnd;
            chunkStart = tentativeEnd;
        }
    }

    // Launch threads
    for (int i = 0; i < numThreads; ++i) {
        threads.emplace_back([i, &chunkStarts, &chunkEnds, &results]() {
            processChunk(chunkStarts[i], chunkEnds[i], i == 0, results[i]);
        });
    }

    // Wait and collect
    size_t totalLines = 0, missingLines = 0;
    for (int i = 0; i < numThreads; ++i) {
        threads[i].join();
        out.append(results[i].output.data(), results[i].output.size());
        totalLines += results[i].totalLines;
        missingLines += results[i].missingLines;
    }

    if (!quiet) {
        std::cerr << "Processed " << totalLines << " variants, "
                  << missingLines << " with missing genotypes ("
                  << (totalLines > 0 ? (100.0 * missingLines / totalLines) : 0.0)
                  << "%) using " << numThreads << " threads\n";
    }
}

// ============================================================================
// Process stdin (fallback) - kept for compatibility
// ============================================================================
void VCFXMissingDetector::detectMissingGenotypes(std::istream &in, std::ostream &out) {
    std::string line;
    line.reserve(65536);

    while (std::getline(in, line)) {
        if (line.empty()) {
            out << "\n";
            continue;
        }

        if (line[0] == '#') {
            out << line << "\n";
            continue;
        }

        // Find sample start (field 9)
        const char* lineStart = line.data();
        const char* lineEnd = lineStart + line.size();
        const char* sampleStart = skipToField(lineStart, lineEnd, 9);

        if (sampleStart >= lineEnd) {
            out << line << "\n";
            continue;
        }

        bool hasMissing = hasMissingGenotypeInSamples(sampleStart, lineEnd);

        if (!hasMissing) {
            out << line << "\n";
            continue;
        }

        // Modify INFO field
        const char* infoStart = skipToField(lineStart, lineEnd, 7);
        const char* infoEnd = infoStart;
        while (infoEnd < lineEnd && *infoEnd != '\t') ++infoEnd;

        out.write(lineStart, infoStart - lineStart);

        std::string_view info(infoStart, infoEnd - infoStart);
        if (info == "." || info.empty()) {
            out << "MISSING_GENOTYPES=1";
        } else {
            out << info;
            if (info.back() != ';') out << ';';
            out << "MISSING_GENOTYPES=1";
        }

        out.write(infoEnd, lineEnd - infoEnd);
        out << "\n";
    }
}

// ============================================================================
// Help message
// ============================================================================
void VCFXMissingDetector::displayHelp() {
    std::cout << "VCFX_missing_detector v2.0 - Extreme-performance missing genotype detector\n\n"
              << "Usage:\n"
              << "  VCFX_missing_detector [OPTIONS] [input.vcf]\n"
              << "  VCFX_missing_detector [OPTIONS] < input.vcf > flagged.vcf\n\n"
              << "Options:\n"
              << "  -i, --input FILE   Input VCF file (uses memory-mapping for best performance)\n"
              << "  -t, --threads N    Number of threads (default: auto)\n"
              << "  -q, --quiet        Suppress informational messages\n"
              << "  -h, --help         Display this help message and exit\n"
              << "  -v, --version      Show program version and exit\n\n"
              << "Description:\n"
              << "  Detects variants with missing sample genotypes and flags them\n"
              << "  with 'MISSING_GENOTYPES=1' in the INFO field.\n\n"
              << "Performance:\n"
              << "  - Memory-mapped I/O: Use -i flag for extreme speed\n"
              << "  - SIMD-accelerated '.' character search (AVX2/SSE2/NEON)\n"
              << "  - Multi-threaded chunk processing\n"
              << "  - Zero-copy output for lines without missing genotypes\n\n"
              << "Example:\n"
              << "  VCFX_missing_detector -i input.vcf > flagged.vcf\n"
              << "  VCFX_missing_detector < input.vcf > flagged.vcf\n";
}

// ============================================================================
// Run method
// ============================================================================
int VCFXMissingDetector::run(int argc, char *argv[]) {
    const char *inputFile = nullptr;
    int numThreads = std::thread::hardware_concurrency();
    if (numThreads == 0) numThreads = 4;
    bool quiet = false;

    static struct option longOpts[] = {
        {"input",   required_argument, nullptr, 'i'},
        {"threads", required_argument, nullptr, 't'},
        {"quiet",   no_argument,       nullptr, 'q'},
        {"help",    no_argument,       nullptr, 'h'},
        {"version", no_argument,       nullptr, 'v'},
        {nullptr, 0, nullptr, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "i:t:qhv", longOpts, nullptr)) != -1) {
        switch (opt) {
            case 'i': inputFile = optarg; break;
            case 't': numThreads = std::max(1, atoi(optarg)); break;
            case 'q': quiet = true; break;
            case 'h': displayHelp(); return 0;
            case 'v': std::cout << "VCFX_missing_detector v2.0\n"; return 0;
            default: displayHelp(); return 1;
        }
    }

    // Check for positional argument
    if (!inputFile && optind < argc) {
        inputFile = argv[optind];
    }

    if (inputFile) {
        // Memory-mapped mode
        MappedFile mf;
        if (!mf.open(inputFile)) {
            std::cerr << "Error: Cannot open file: " << inputFile << "\n";
            return 1;
        }

        if (!quiet) {
            std::cerr << "Processing " << inputFile << " (" << (mf.size / (1024*1024)) << " MB)\n";
        }

        // Use ultra-fast zero-copy with multi-threaded pre-scan
        // This batches consecutive unchanged lines and writes directly from mmap
        processMmapZeroCopy(mf, STDOUT_FILENO, quiet, numThreads);
        mf.close();
    } else {
        // Stdin mode
        detectMissingGenotypes(std::cin, std::cout);
    }

    return 0;
}

// ============================================================================
// Main
// ============================================================================
static void show_help() {
    VCFXMissingDetector obj;
    char arg0[] = "VCFX_missing_detector";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);

    VCFXMissingDetector detector;
    return detector.run(argc, argv);
}
