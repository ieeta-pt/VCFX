#include "VCFX_haplotype_extractor.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string_view>
#include <getopt.h>

// Memory-mapped I/O includes
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

// SIMD detection
#if defined(__AVX2__)
#define USE_AVX2
#include <immintrin.h>
#elif defined(__SSE2__) || defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2)
#define USE_SSE2
#include <emmintrin.h>
#endif

// Branch prediction hints
#if defined(__GNUC__) || defined(__clang__)
#define LIKELY(x) __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

// =============================================================================
// MappedFile - RAII wrapper for memory-mapped files
// =============================================================================
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

        // Advise kernel for sequential access and read-ahead
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
// OutputBuffer - 1MB buffered output for reduced syscalls
// =============================================================================
class OutputBuffer {
    static constexpr size_t BUFFER_SIZE = 1024 * 1024;  // 1MB
    char *buffer;
    size_t pos = 0;
    std::ostream &out;

public:
    explicit OutputBuffer(std::ostream &os) : out(os) {
        buffer = new char[BUFFER_SIZE];
    }

    ~OutputBuffer() {
        flush();
        delete[] buffer;
    }

    void write(std::string_view sv) {
        if (pos + sv.size() > BUFFER_SIZE) flush();
        if (sv.size() > BUFFER_SIZE) {
            out.write(sv.data(), static_cast<std::streamsize>(sv.size()));
            return;
        }
        memcpy(buffer + pos, sv.data(), sv.size());
        pos += sv.size();
    }

    void writeChar(char c) {
        if (pos + 1 > BUFFER_SIZE) flush();
        buffer[pos++] = c;
    }

    void writeInt(int n) {
        char tmp[16];
        int len = snprintf(tmp, sizeof(tmp), "%d", n);
        if (len > 0) write(std::string_view(tmp, static_cast<size_t>(len)));
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
    // Fallback for remainder
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
// Zero-copy parsing helpers
// =============================================================================

// Find GT index in FORMAT field (zero-copy)
static inline int findGTIndex(std::string_view format) {
    const char *p = format.data();
    const char *end = p + format.size();
    int idx = 0;
    const char *fieldStart = p;

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

// Extract nth colon-delimited field (zero-copy)
static inline std::string_view extractNthField(std::string_view str, int n) {
    if (n < 0) return {};
    const char *p = str.data();
    const char *end = p + str.size();
    int fieldIdx = 0;
    const char *fieldStart = p;

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

// Split line by tabs into string_view vector (zero-copy)
static inline void splitTabsView(std::string_view line,
                                  std::vector<std::string_view> &out) {
    out.clear();
    const char *p = line.data();
    const char *end = p + line.size();
    const char *fieldStart = p;

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
// printHelp
// =============================================================================
void printHelp() {
    std::cout << "VCFX_haplotype_extractor\n"
              << "Usage: VCFX_haplotype_extractor [OPTIONS] [input.vcf]\n\n"
              << "Options:\n"
              << "  -h, --help                 Display this help message and exit.\n"
              << "  -i, --input FILE           Input VCF file (uses fast memory-mapped I/O).\n"
              << "  -b, --block-size <int>     Maximum distance for grouping consecutive variants (default 100000).\n"
              << "  -c, --check-phase-consistency  If set, try a minimal check across variants.\n"
              << "  -s, --streaming            Enable streaming mode: output blocks immediately when complete.\n"
              << "                             Uses O(block_size) memory instead of O(total_variants).\n"
              << "  -q, --quiet                Suppress warning messages.\n"
              << "  -d, --debug                Output verbose debug information.\n\n"
              << "Description:\n"
              << "  Extracts phased haplotype blocks from genotype data in a VCF file. "
              << "It reconstructs haplotypes for each sample by analyzing phased genotype fields.\n\n"
              << "Performance:\n"
              << "  File input (-i): Uses memory-mapped I/O for 50-100x faster processing.\n"
              << "  Streaming mode:  Outputs blocks immediately when complete. Enables\n"
              << "                   processing of arbitrarily large files with bounded memory.\n\n"
              << "Examples:\n"
              << "  ./VCFX_haplotype_extractor -i phased.vcf > haplotypes.tsv\n"
              << "  ./VCFX_haplotype_extractor -b 50000 < phased.vcf > haplotypes.tsv\n"
              << "  ./VCFX_haplotype_extractor -s -i large_phased.vcf > haplotypes.tsv\n"
              << "  ./VCFX_haplotype_extractor -s -b 10000 -i phased.vcf > haplotypes.tsv\n";
}

// =============================================================================
// HaplotypeExtractor implementation
// =============================================================================

// parseHeader: parse #CHROM line => sample columns (optimized with string_view)
bool HaplotypeExtractor::parseHeader(std::string_view headerLine) {
    std::vector<std::string_view> fields;
    fields.reserve(16);
    splitTabsView(headerLine, fields);

    if (fields.size() <= 9) {
        std::cerr << "Error: VCF header does not contain sample columns.\n";
        return false;
    }

    sampleNames.clear();
    sampleNames.reserve(fields.size() - 9);
    for (size_t i = 9; i < fields.size(); i++) {
        sampleNames.emplace_back(fields[i]);
    }
    numSamples = sampleNames.size();
    return true;
}

// areAllSamplesPhased: checks each genotype for '|'
bool HaplotypeExtractor::areAllSamplesPhased(const std::vector<std::string_view> &genotypes) {
    for (const auto &g : genotypes) {
        bool hasPhase = false;
        for (char c : g) {
            if (c == '|') {
                hasPhase = true;
                break;
            }
        }
        if (!hasPhase) return false;
    }
    return true;
}

// phaseIsConsistent: O(1) check using cached last genotypes
bool HaplotypeExtractor::phaseIsConsistent(const HaplotypeBlock &block,
                                            const std::vector<std::string_view> &newGenotypes) {
    if (block.lastGenotypes.size() != newGenotypes.size()) {
        return false;
    }

    if (debugMode) {
        std::cerr << "Checking phase consistency\n";
    }

    for (size_t s = 0; s < block.lastGenotypes.size(); s++) {
        const std::string &lastGT = block.lastGenotypes[s];
        const std::string_view &newGT = newGenotypes[s];

        if (debugMode) {
            std::cerr << "Sample " << s << " last GT: " << lastGT << " new GT: " << newGT << "\n";
        }

        if (lastGT.size() < 3 || newGT.size() < 3) {
            continue;
        }

        char lastAllele1 = lastGT[0];
        char lastAllele2 = lastGT[lastGT.size() - 1];
        char newAllele1 = newGT[0];
        char newAllele2 = newGT[newGT.size() - 1];

        if (debugMode) {
            std::cerr << "Comparing alleles: " << lastAllele1 << "|" << lastAllele2
                      << " vs " << newAllele1 << "|" << newAllele2 << "\n";
        }

        // Detect phase flip
        if (lastAllele1 != newAllele1 && lastAllele2 != newAllele2 &&
            lastAllele1 == newAllele2 && lastAllele2 == newAllele1) {
            if (debugMode) {
                std::cerr << "Phase flip detected in sample " << s << "\n";
            }
            return false;
        }
    }

    if (debugMode) {
        std::cerr << "All phases consistent\n";
    }
    return true;
}

// outputHeader: write header line to OutputBuffer
void HaplotypeExtractor::outputHeader(OutputBuffer &out) {
    out.write("CHROM\tSTART\tEND");
    for (const auto &s : sampleNames) {
        out.writeChar('\t');
        out.write(s);
    }
    out.writeChar('\n');
}

// outputBlock: write a single block to OutputBuffer
void HaplotypeExtractor::outputBlock(OutputBuffer &out, const HaplotypeBlock &block) {
    out.write(block.chrom);
    out.writeChar('\t');
    out.writeInt(block.start);
    out.writeChar('\t');
    out.writeInt(block.end);
    for (const auto &hap : block.haplotypes) {
        out.writeChar('\t');
        out.write(hap);
    }
    out.writeChar('\n');
}

// processVariantFast: unified fast processing for both modes
template <bool Streaming>
bool HaplotypeExtractor::processVariantFast(
    const std::vector<std::string_view> &fields,
    std::vector<HaplotypeBlock> &haplotypeBlocks,
    HaplotypeBlock &currentBlock,
    bool &hasCurrentBlock,
    OutputBuffer &out) {

    if (UNLIKELY(fields.size() < 10)) {
        if (!quiet_) {
            std::cerr << "Warning: skipping invalid VCF line (<10 fields)\n";
        }
        return false;
    }

    const std::string_view &chrom = fields[0];
    int pos = 0;

    // Fast integer parsing
    const char *p = fields[1].data();
    const char *end = p + fields[1].size();
    while (p < end) {
        if (UNLIKELY(*p < '0' || *p > '9')) {
            if (!quiet_) {
                std::cerr << "Warning: invalid POS => skip variant\n";
            }
            return false;
        }
        pos = pos * 10 + (*p - '0');
        p++;
    }

    // FORMAT caching - only recalculate GT index when FORMAT changes
    const std::string_view format = fields[8];
    int gtIndex;
    if (format.size() == cachedFormat_.size() &&
        memcmp(format.data(), cachedFormat_.data(), format.size()) == 0) {
        gtIndex = cachedGTIndex_;
    } else {
        cachedFormat_.assign(format.data(), format.size());
        gtIndex = findGTIndex(format);
        cachedGTIndex_ = gtIndex;
    }

    if (UNLIKELY(gtIndex < 0)) {
        return false;
    }

    // Extract genotypes (reuse vector)
    genotypeFields_.resize(numSamples);
    bool allPhased = true;

    for (size_t s = 0; s < numSamples; s++) {
        if (UNLIKELY(9 + s >= fields.size())) {
            genotypeFields_[s] = ".|.";
            allPhased = false;
            continue;
        }
        std::string_view gt = extractNthField(fields[9 + s], gtIndex);
        if (UNLIKELY(gt.empty())) {
            genotypeFields_[s] = ".|.";
            allPhased = false;
            continue;
        }

        // Check for phased genotype (contains '|')
        bool hasPhase = false;
        for (char c : gt) {
            if (c == '|') {
                hasPhase = true;
                break;
            }
        }
        if (!hasPhase) {
            allPhased = false;
        }
        genotypeFields_[s] = gt;
    }

    if (UNLIKELY(!allPhased)) {
        if (!quiet_) {
            std::cerr << "Warning: Not all samples phased at " << chrom << ":" << pos << ".\n";
        }
        return false;
    }

    // Process based on mode
    if constexpr (Streaming) {
        // Streaming mode: output completed blocks immediately
        if (!hasCurrentBlock) {
            // Start first block
            currentBlock.chrom.assign(chrom.data(), chrom.size());
            currentBlock.start = pos;
            currentBlock.end = pos;
            currentBlock.haplotypes.resize(numSamples);
            currentBlock.lastGenotypes.resize(numSamples);

            // Pre-reserve space to avoid reallocations
            size_t estimatedVariants = static_cast<size_t>(blockDistanceThreshold) / 50;
            for (size_t s = 0; s < numSamples; s++) {
                currentBlock.haplotypes[s].reserve(estimatedVariants * 4);
                currentBlock.haplotypes[s].assign(genotypeFields_[s].data(), genotypeFields_[s].size());
                currentBlock.lastGenotypes[s].assign(genotypeFields_[s].data(), genotypeFields_[s].size());
            }
            hasCurrentBlock = true;
        } else {
            // Check if we can extend current block
            bool canExtend = (chrom.size() == currentBlock.chrom.size() &&
                              memcmp(chrom.data(), currentBlock.chrom.data(), chrom.size()) == 0) &&
                             (pos - currentBlock.end <= blockDistanceThreshold);

            if (canExtend && checkPhaseConsistency) {
                canExtend = phaseIsConsistent(currentBlock, genotypeFields_);
            }

            if (canExtend) {
                // Extend current block
                currentBlock.end = pos;
                for (size_t s = 0; s < numSamples; s++) {
                    currentBlock.haplotypes[s] += '|';
                    currentBlock.haplotypes[s].append(genotypeFields_[s].data(), genotypeFields_[s].size());
                    currentBlock.lastGenotypes[s].assign(genotypeFields_[s].data(), genotypeFields_[s].size());
                }
            } else {
                // Output completed block
                outputBlock(out, currentBlock);

                // Start new block
                currentBlock.chrom.assign(chrom.data(), chrom.size());
                currentBlock.start = pos;
                currentBlock.end = pos;
                for (size_t s = 0; s < numSamples; s++) {
                    currentBlock.haplotypes[s].clear();
                    currentBlock.haplotypes[s].assign(genotypeFields_[s].data(), genotypeFields_[s].size());
                    currentBlock.lastGenotypes[s].assign(genotypeFields_[s].data(), genotypeFields_[s].size());
                }
            }
        }
    } else {
        // Batch mode: accumulate all blocks
        if (haplotypeBlocks.empty()) {
            HaplotypeBlock b;
            b.chrom.assign(chrom.data(), chrom.size());
            b.start = pos;
            b.end = pos;
            b.haplotypes.resize(numSamples);
            b.lastGenotypes.resize(numSamples);

            size_t estimatedVariants = static_cast<size_t>(blockDistanceThreshold) / 50;
            for (size_t s = 0; s < numSamples; s++) {
                b.haplotypes[s].reserve(estimatedVariants * 4);
                b.haplotypes[s].assign(genotypeFields_[s].data(), genotypeFields_[s].size());
                b.lastGenotypes[s].assign(genotypeFields_[s].data(), genotypeFields_[s].size());
            }
            haplotypeBlocks.push_back(std::move(b));
        } else {
            HaplotypeBlock &lastB = haplotypeBlocks.back();
            bool canExtend = (chrom.size() == lastB.chrom.size() &&
                              memcmp(chrom.data(), lastB.chrom.data(), chrom.size()) == 0) &&
                             (pos - lastB.end <= blockDistanceThreshold);

            if (canExtend && checkPhaseConsistency) {
                canExtend = phaseIsConsistent(lastB, genotypeFields_);
            }

            if (canExtend) {
                lastB.end = pos;
                for (size_t s = 0; s < numSamples; s++) {
                    lastB.haplotypes[s] += '|';
                    lastB.haplotypes[s].append(genotypeFields_[s].data(), genotypeFields_[s].size());
                    lastB.lastGenotypes[s].assign(genotypeFields_[s].data(), genotypeFields_[s].size());
                }
            } else {
                HaplotypeBlock nb;
                nb.chrom.assign(chrom.data(), chrom.size());
                nb.start = pos;
                nb.end = pos;
                nb.haplotypes.resize(numSamples);
                nb.lastGenotypes.resize(numSamples);

                size_t estimatedVariants = static_cast<size_t>(blockDistanceThreshold) / 50;
                for (size_t s = 0; s < numSamples; s++) {
                    nb.haplotypes[s].reserve(estimatedVariants * 4);
                    nb.haplotypes[s].assign(genotypeFields_[s].data(), genotypeFields_[s].size());
                    nb.lastGenotypes[s].assign(genotypeFields_[s].data(), genotypeFields_[s].size());
                }
                haplotypeBlocks.push_back(std::move(nb));
            }
        }
    }

    return true;
}

// Explicit template instantiations
template bool HaplotypeExtractor::processVariantFast<true>(
    const std::vector<std::string_view> &, std::vector<HaplotypeBlock> &,
    HaplotypeBlock &, bool &, OutputBuffer &);
template bool HaplotypeExtractor::processVariantFast<false>(
    const std::vector<std::string_view> &, std::vector<HaplotypeBlock> &,
    HaplotypeBlock &, bool &, OutputBuffer &);

// =============================================================================
// Memory-mapped file processing (fast path)
// =============================================================================
bool HaplotypeExtractor::extractHaplotypesMmap(const char *filename, std::ostream &os) {
    MappedFile mf;
    if (!mf.open(filename)) {
        std::cerr << "Error: could not open file: " << filename << "\n";
        return false;
    }

    if (mf.size == 0) {
        std::cerr << "Error: empty file\n";
        return false;
    }

    OutputBuffer out(os);
    std::vector<HaplotypeBlock> haplotypeBlocks;
    std::vector<std::string_view> fields;
    fields.reserve(16);

    HaplotypeBlock dummyBlock;
    bool dummyHasBlock = false;

    const char *p = mf.data;
    const char *end = mf.data + mf.size;
    bool foundHeader = false;

    while (p < end) {
        const char *lineEnd = findNewlineSIMD(p, end);
        if (!lineEnd) lineEnd = end;

        size_t lineLen = static_cast<size_t>(lineEnd - p);
        if (lineLen > 0 && p[lineLen - 1] == '\r') lineLen--;  // Handle \r\n

        if (lineLen == 0) {
            p = lineEnd + 1;
            continue;
        }

        std::string_view line(p, lineLen);

        if (line[0] == '#') {
            if (!foundHeader && line.size() >= 6 &&
                line[0] == '#' && line[1] == 'C' && line[2] == 'H' &&
                line[3] == 'R' && line[4] == 'O' && line[5] == 'M') {
                foundHeader = true;
                if (!parseHeader(line)) {
                    return false;
                }
            }
            p = lineEnd + 1;
            continue;
        }

        if (!foundHeader) {
            std::cerr << "Error: no #CHROM header found before data.\n";
            return false;
        }

        splitTabsView(line, fields);
        processVariantFast<false>(fields, haplotypeBlocks, dummyBlock, dummyHasBlock, out);

        p = lineEnd + 1;
    }

    // Output header and all blocks
    outputHeader(out);
    for (const auto &block : haplotypeBlocks) {
        outputBlock(out, block);
    }

    return true;
}

bool HaplotypeExtractor::extractHaplotypesMmapStreaming(const char *filename, std::ostream &os) {
    MappedFile mf;
    if (!mf.open(filename)) {
        std::cerr << "Error: could not open file: " << filename << "\n";
        return false;
    }

    if (mf.size == 0) {
        std::cerr << "Error: empty file\n";
        return false;
    }

    OutputBuffer out(os);
    std::vector<std::string_view> fields;
    fields.reserve(16);

    HaplotypeBlock currentBlock;
    bool hasCurrentBlock = false;
    std::vector<HaplotypeBlock> dummyBlocks;

    const char *p = mf.data;
    const char *end = mf.data + mf.size;
    bool foundHeader = false;

    while (p < end) {
        const char *lineEnd = findNewlineSIMD(p, end);
        if (!lineEnd) lineEnd = end;

        size_t lineLen = static_cast<size_t>(lineEnd - p);
        if (lineLen > 0 && p[lineLen - 1] == '\r') lineLen--;

        if (lineLen == 0) {
            p = lineEnd + 1;
            continue;
        }

        std::string_view line(p, lineLen);

        if (line[0] == '#') {
            if (!foundHeader && line.size() >= 6 &&
                line[0] == '#' && line[1] == 'C' && line[2] == 'H' &&
                line[3] == 'R' && line[4] == 'O' && line[5] == 'M') {
                foundHeader = true;
                if (!parseHeader(line)) {
                    return false;
                }
                outputHeader(out);
            }
            p = lineEnd + 1;
            continue;
        }

        if (!foundHeader) {
            std::cerr << "Error: no #CHROM header found before data.\n";
            return false;
        }

        splitTabsView(line, fields);
        processVariantFast<true>(fields, dummyBlocks, currentBlock, hasCurrentBlock, out);

        p = lineEnd + 1;
    }

    // Output final block
    if (hasCurrentBlock) {
        outputBlock(out, currentBlock);
    }

    return true;
}

// =============================================================================
// Stdin processing (fallback)
// =============================================================================
bool HaplotypeExtractor::extractHaplotypes(std::istream &in, std::ostream &os) {
    OutputBuffer out(os);
    std::vector<HaplotypeBlock> haplotypeBlocks;
    std::vector<std::string_view> fields;
    fields.reserve(16);

    HaplotypeBlock dummyBlock;
    bool dummyHasBlock = false;

    std::string line;
    bool foundHeader = false;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        // Strip \r if present
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        if (line[0] == '#') {
            if (!foundHeader && line.rfind("#CHROM", 0) == 0) {
                foundHeader = true;
                if (!parseHeader(line)) {
                    return false;
                }
            }
            continue;
        }

        if (!foundHeader) {
            std::cerr << "Error: no #CHROM header found before data.\n";
            return false;
        }

        splitTabsView(line, fields);
        processVariantFast<false>(fields, haplotypeBlocks, dummyBlock, dummyHasBlock, out);
    }

    // Output header and all blocks
    outputHeader(out);
    for (const auto &block : haplotypeBlocks) {
        outputBlock(out, block);
    }

    return true;
}

bool HaplotypeExtractor::extractHaplotypesStreaming(std::istream &in, std::ostream &os) {
    OutputBuffer out(os);
    std::vector<std::string_view> fields;
    fields.reserve(16);

    HaplotypeBlock currentBlock;
    bool hasCurrentBlock = false;
    std::vector<HaplotypeBlock> dummyBlocks;

    std::string line;
    bool foundHeader = false;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        if (line[0] == '#') {
            if (!foundHeader && line.rfind("#CHROM", 0) == 0) {
                foundHeader = true;
                if (!parseHeader(line)) {
                    return false;
                }
                outputHeader(out);
            }
            continue;
        }

        if (!foundHeader) {
            std::cerr << "Error: no #CHROM header found before data.\n";
            return false;
        }

        splitTabsView(line, fields);
        processVariantFast<true>(fields, dummyBlocks, currentBlock, hasCurrentBlock, out);
    }

    // Output final block
    if (hasCurrentBlock) {
        outputBlock(out, currentBlock);
    }

    return true;
}

// =============================================================================
// main
// =============================================================================
static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();
    if (vcfx::handle_common_flags(argc, argv, "VCFX_haplotype_extractor", show_help))
        return 0;

    int blockSize = 100000;
    bool doCheck = false;
    bool debug = false;
    bool streaming = false;
    bool quiet = false;
    std::string inputFile;

    static struct option long_opts[] = {
        {"help", no_argument, nullptr, 'h'},
        {"input", required_argument, nullptr, 'i'},
        {"block-size", required_argument, nullptr, 'b'},
        {"check-phase-consistency", no_argument, nullptr, 'c'},
        {"streaming", no_argument, nullptr, 's'},
        {"quiet", no_argument, nullptr, 'q'},
        {"debug", no_argument, nullptr, 'd'},
        {nullptr, 0, nullptr, 0}
    };

    optind = 1;
    while (true) {
        int c = getopt_long(argc, argv, "hi:b:csqd", long_opts, nullptr);
        if (c == -1) break;

        switch (c) {
        case 'h':
            printHelp();
            return 0;
        case 'i':
            inputFile = optarg;
            break;
        case 'b':
            blockSize = std::stoi(optarg);
            break;
        case 'c':
            doCheck = true;
            break;
        case 's':
            streaming = true;
            break;
        case 'q':
            quiet = true;
            break;
        case 'd':
            debug = true;
            break;
        default:
            printHelp();
            return 1;
        }
    }

    // Positional argument support
    if (inputFile.empty() && optind < argc) {
        inputFile = argv[optind];
    }

    HaplotypeExtractor extractor;
    extractor.setBlockDistanceThreshold(blockSize);
    extractor.setCheckPhaseConsistency(doCheck);
    extractor.setDebug(debug);
    extractor.setStreamingMode(streaming);
    extractor.setQuiet(quiet);

    bool ok;
    if (!inputFile.empty()) {
        // Fast memory-mapped path
        if (streaming) {
            ok = extractor.extractHaplotypesMmapStreaming(inputFile.c_str(), std::cout);
        } else {
            ok = extractor.extractHaplotypesMmap(inputFile.c_str(), std::cout);
        }
    } else {
        // Stdin path
        if (streaming) {
            ok = extractor.extractHaplotypesStreaming(std::cin, std::cout);
        } else {
            ok = extractor.extractHaplotypes(std::cin, std::cout);
        }
    }

    return (ok ? 0 : 1);
}
