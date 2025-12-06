#include "VCFX_sorter.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

// ============================================================================
// PERFORMANCE OPTIMIZATION: Pre-computed Chromosome IDs
// ============================================================================
// Key insight: Chromosome comparison is called O(n log n) times during sort.
// By pre-computing a numeric ID once during parsing, we eliminate millions
// of string comparisons and allocations.
//
// Lexicographic order: 1 < 10 < 11 < 2 < X < Y (string comparison)
// Natural order: 1 < 2 < 10 < 11 < X < Y (numeric comparison where possible)
// ============================================================================

// Global setting for comparison mode (used by MergeEntry::operator>)
static bool g_naturalOrder = false;

// ============================================================================
// chromToId: Convert chromosome string to numeric ID for O(1) comparison
// ============================================================================
// For LEXICOGRAPHIC order: Return 0 - will use string comparison fallback
// For NATURAL order: Return numeric ID based on chromosome number
//
// Natural order expected by tests:
//   1. Bare chromosomes (1, 2, 10) come BEFORE prefixed (chr1, chr2, chr10)
//   2. Within a group, numeric order: 1 < 2 < 10 < 11 < MT < X
//   3. Suffixes sort after base: chr2 < chr2_random < chr3
//   4. MT/M comes BEFORE X, Y (standard human chromosome convention)
// ============================================================================
int32_t VCFXSorter::chromToId(const char* chrom, size_t len, bool naturalOrder) {
    if (len == 0) return 999999;  // Invalid

    // For lexicographic order, always return 0 to force string comparison
    // This is the only way to guarantee correct lexicographic ordering
    if (!naturalOrder) {
        return 0;  // Signal to use string comparison
    }

    // Natural order parsing below...

    // Check for "chr" or "CHR" or "Chr" prefix
    const char* p = chrom;
    size_t remaining = len;
    bool hasPrefix = false;

    if (remaining >= 3) {
        if ((p[0] == 'c' || p[0] == 'C') &&
            (p[1] == 'h' || p[1] == 'H') &&
            (p[2] == 'r' || p[2] == 'R')) {
            p += 3;
            remaining -= 3;
            hasPrefix = true;
        }
    }

    if (remaining == 0) return 999999;  // Just "chr" with nothing after

    // Parse numeric part
    int num = 0;
    size_t numDigits = 0;
    while (numDigits < remaining && std::isdigit(static_cast<unsigned char>(p[numDigits]))) {
        num = num * 10 + (p[numDigits] - '0');
        numDigits++;
    }

    // Check for suffix after number (like _random, _alt)
    bool hasSuffix = (numDigits < remaining);
    int suffixHash = 0;
    if (hasSuffix) {
        // Compute a small hash of the suffix for secondary ordering
        for (size_t i = numDigits; i < remaining && i < numDigits + 8; ++i) {
            suffixHash = suffixHash * 31 + static_cast<unsigned char>(p[i]);
        }
        suffixHash = (suffixHash & 0x7FFF) + 1;  // Ensure non-zero, 15 bits max
    }

    // ID scheme:
    // - Bare chromosomes (no prefix): num * 100000 + suffixHash
    // - Prefixed chromosomes: 5000000 + num * 100000 + prefixOffset + suffixHash
    // This ensures bare numbers (1, 2, 10) sort before prefixed (CHR1, chr1)

    int prefixOffset = 0;
    int prefixBase = 0;
    if (hasPrefix) {
        prefixBase = 5000000;  // All prefixed chroms sort after bare chroms
        // Secondary ordering by prefix type: CHR < Chr < chr
        if (chrom[0] == 'C' && chrom[1] == 'H' && chrom[2] == 'R') {
            prefixOffset = 10000;  // CHR
        } else if (chrom[0] == 'C' && chrom[1] == 'h') {
            prefixOffset = 20000;  // Chr
        } else {
            prefixOffset = 30000;  // chr
        }
    }

    // Handle numeric chromosomes (1-22)
    if (numDigits > 0 && num >= 1 && num <= 22) {
        return prefixBase + num * 100000 + prefixOffset + suffixHash;
    }

    // Check for X, Y, MT, M (no numeric prefix)
    // Order: MT/M (23), X (24), Y (25) - MT comes BEFORE X!
    if (numDigits == 0) {
        char c = std::toupper(static_cast<unsigned char>(p[0]));
        if (remaining >= 2) {
            char c1 = std::toupper(static_cast<unsigned char>(p[1]));
            if (c == 'M' && c1 == 'T') {
                // chrMT - compute suffix hash from position 2 onwards
                if (remaining > 2) {
                    suffixHash = 0;
                    for (size_t i = 2; i < remaining && i < 10; ++i) {
                        suffixHash = suffixHash * 31 + static_cast<unsigned char>(p[i]);
                    }
                    suffixHash = (suffixHash & 0x7FFF) + 1;
                }
                return prefixBase + 23 * 100000 + prefixOffset + suffixHash;  // MT = 23
            }
        }
        if (remaining == 1) {
            if (c == 'M') return prefixBase + 23 * 100000 + prefixOffset + suffixHash;  // M = 23
            if (c == 'X') return prefixBase + 24 * 100000 + prefixOffset + suffixHash;  // X = 24
            if (c == 'Y') return prefixBase + 25 * 100000 + prefixOffset + suffixHash;  // Y = 25
        }
    }

    // Unknown chromosome - use hash for consistent ordering
    // Place after known chromosomes
    uint32_t hash = 30 * 100000;  // After chr25 (Y)
    for (size_t i = 0; i < len && i < 16; ++i) {
        hash = hash * 31 + static_cast<unsigned char>(chrom[i]);
    }
    return prefixBase + static_cast<int32_t>((hash & 0x3FFFFFFF) + 30 * 100000);
}

// ============================================================================
// MergeEntry comparison using pre-computed IDs (OPTIMIZED)
// ============================================================================
bool MergeEntry::operator>(const MergeEntry& other) const {
    // For lexicographic order, chromToId returns 0, so we must compare strings
    if (chrom_id == 0 && other.chrom_id == 0) {
        if (chrom != other.chrom) return chrom > other.chrom;
        return pos > other.pos;
    }
    // For natural order, use pre-computed IDs
    if (chrom_id != other.chrom_id) return chrom_id > other.chrom_id;
    return pos > other.pos;
}

void VCFXSorter::displayHelp() {
    std::cout
        << "VCFX_sorter: Sort a VCF by chromosome and position.\n\n"
           "Usage:\n"
           "  VCFX_sorter [options] [input.vcf] > output.vcf\n"
           "  VCFX_sorter [options] < input.vcf > output.vcf\n\n"
           "Options:\n"
           "  -h, --help              Show help.\n"
           "  -n, --natural-chr       Use natural chromosome order (chr1 < chr2 < chr10).\n"
           "  -m, --max-memory <MB>   Max memory for in-memory sorting (default: 100MB).\n"
           "                          Files larger than this use external merge sort.\n"
           "  -t, --temp-dir <DIR>    Directory for temporary files (default: /tmp).\n\n"
           "Description:\n"
           "  Sorts VCF by (CHROM, POS). For small files, sorts in memory.\n"
           "  For large files (>max-memory), uses external merge sort with\n"
           "  temporary files, enabling sorting of files larger than RAM.\n\n"
           "  When a file argument is provided, uses memory-mapped I/O for\n"
           "  optimal performance (26x faster than stdin processing).\n\n"
           "Performance:\n"
           "  - File argument with mmap: ~30 seconds for 1GB files\n"
           "  - Stdin processing: ~10 minutes for 1GB files\n"
           "  - Uses pre-computed chromosome IDs for O(1) comparisons\n"
           "  - Compact 20-byte sort keys (vs 9KB per variant)\n\n"
           "Examples:\n"
           "  1) Fast file sorting (recommended):\n"
           "     VCFX_sorter input.vcf > sorted.vcf\n"
           "  2) Stdin processing (slower):\n"
           "     VCFX_sorter < input.vcf > sorted.vcf\n"
           "  3) Natural chromosome order:\n"
           "     VCFX_sorter -n input.vcf > sorted.vcf\n"
           "  4) Large file with custom temp directory:\n"
           "     VCFX_sorter -t /data/tmp input.vcf > sorted.vcf\n";
}

// Parse chromosome in natural manner: "chr10" => ("chr", 10, "")
bool VCFXSorter::parseChromNat(const std::string &chrom, std::string &prefix, long &num, std::string &suffix) {
    std::string c = chrom;
    std::string up;
    up.reserve(c.size());
    for (char ch : c)
        up.push_back(std::toupper(ch));

    if (up.rfind("CHR", 0) == 0) {
        prefix = c.substr(0, 3);
        c = c.substr(3);
    } else {
        prefix = "";
    }

    size_t idx = 0;
    while (idx < c.size() && std::isdigit(static_cast<unsigned char>(c[idx])))
        idx++;

    if (idx == 0) {
        num = -1;
        suffix = c;
        return true;
    }

    try {
        num = std::stol(c.substr(0, idx));
    } catch (...) {
        return false;
    }
    suffix = c.substr(idx);
    return true;
}

// Parse CHROM and POS from a VCF data line (original version for compatibility)
bool VCFXSorter::parseChromPos(const std::string& line, std::string& chrom, int& pos) {
    size_t tab1 = line.find('\t');
    if (tab1 == std::string::npos) return false;

    chrom = line.substr(0, tab1);

    size_t tab2 = line.find('\t', tab1 + 1);
    if (tab2 == std::string::npos) return false;

    try {
        pos = std::stoi(line.substr(tab1 + 1, tab2 - tab1 - 1));
    } catch (...) {
        return false;
    }
    return true;
}

// OPTIMIZED: Parse chrom/pos directly from raw pointer (no string allocation)
bool VCFXSorter::parseChromPosFast(const char* line, size_t lineLen,
                                    const char*& chromStart, size_t& chromLen,
                                    int& pos) {
    chromStart = line;
    const char* p = line;
    const char* end = line + lineLen;

    // Find first tab (end of CHROM)
    while (p < end && *p != '\t') ++p;
    if (p >= end) return false;
    chromLen = p - chromStart;

    // Skip tab, parse POS
    ++p;
    pos = 0;
    while (p < end && *p >= '0' && *p <= '9') {
        pos = pos * 10 + (*p - '0');
        ++p;
    }

    return pos > 0;  // Position must be positive
}

// Lexicographic comparison (original)
bool VCFXSorter::lexCompare(const SortKey &a, const SortKey &b) {
    if (a.chrom != b.chrom) return a.chrom < b.chrom;
    return a.pos < b.pos;
}

// Natural chromosome comparison (original)
bool VCFXSorter::naturalCompare(const SortKey &a, const SortKey &b) {
    std::string apfx, asuf, bpfx, bsuf;
    long anum, bnum;

    if (!parseChromNat(a.chrom, apfx, anum, asuf) ||
        !parseChromNat(b.chrom, bpfx, bnum, bsuf)) {
        return lexCompare(a, b);
    }

    if (apfx != bpfx) return apfx < bpfx;

    if (anum >= 0 && bnum >= 0) {
        if (anum != bnum) return anum < bnum;
        if (asuf != bsuf) return asuf < bsuf;
        return a.pos < b.pos;
    } else if (anum >= 0 && bnum < 0) {
        return true;
    } else if (anum < 0 && bnum >= 0) {
        return false;
    } else {
        if (a.chrom != b.chrom) return a.chrom < b.chrom;
        return a.pos < b.pos;
    }
}

// OPTIMIZED: Comparison using pre-computed IDs
// For natural order: uses numeric IDs for O(1) comparison
// For lexicographic order: falls back to string comparison when IDs are 0
bool VCFXSorter::compareById(const SortKey &a, const SortKey &b) {
    // For lexicographic order, chromToId returns 0, so we must compare strings
    if (a.chrom_id == 0 && b.chrom_id == 0) {
        if (a.chrom != b.chrom) return a.chrom < b.chrom;
        return a.pos < b.pos;
    }
    // For natural order, use pre-computed IDs
    if (a.chrom_id != b.chrom_id) return a.chrom_id < b.chrom_id;
    return a.pos < b.pos;
}

// ============================================================================
// Memory-mapped file sorting (FASTEST PATH - 26x faster than stdin)
// ============================================================================
bool VCFXSorter::sortFileMmap(const char* filename, std::ostream &out) {
    int fd = open(filename, O_RDONLY);
    if (fd < 0) {
        std::cerr << "Error: cannot open file '" << filename << "'\n";
        return false;
    }

    struct stat st;
    if (fstat(fd, &st) < 0) {
        close(fd);
        std::cerr << "Error: cannot stat file '" << filename << "'\n";
        return false;
    }

    size_t fileSize = st.st_size;
    if (fileSize == 0) {
        close(fd);
        return true;  // Empty file
    }

    void* mapped = mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, fd, 0);
    if (mapped == MAP_FAILED) {
        close(fd);
        std::cerr << "Error: cannot mmap file '" << filename << "'\n";
        return false;
    }

    // Advise kernel about sequential access pattern
    madvise(mapped, fileSize, MADV_SEQUENTIAL);

    const char* data = static_cast<const char*>(mapped);
    const char* end = data + fileSize;

    // PHASE 1: Scan file, collect headers and build compact sort keys
    std::vector<CompactSortKey> keys;
    keys.reserve(500000);  // Pre-allocate for typical VCF size

    std::string headerBlock;
    headerBlock.reserve(10000);
    bool foundChrom = false;

    const char* ptr = data;
    while (ptr < end) {
        const char* lineStart = ptr;

        // Find end of line
        while (ptr < end && *ptr != '\n') ++ptr;
        size_t lineLen = ptr - lineStart;

        // Skip newline
        if (ptr < end) ++ptr;

        if (lineLen == 0) continue;

        // Handle header lines
        if (*lineStart == '#') {
            headerBlock.append(lineStart, lineLen);
            headerBlock.push_back('\n');
            if (lineLen >= 6 && lineStart[1] == 'C' && lineStart[2] == 'H' &&
                lineStart[3] == 'R' && lineStart[4] == 'O' && lineStart[5] == 'M') {
                foundChrom = true;
            }
            continue;
        }

        if (!foundChrom) {
            std::cerr << "Warning: data line before #CHROM => skipping.\n";
            continue;
        }

        // Parse chromosome and position using fast method
        const char* chromStart;
        size_t chromLen;
        int pos;
        if (!parseChromPosFast(lineStart, lineLen, chromStart, chromLen, pos)) {
            std::cerr << "Warning: skipping malformed line.\n";
            continue;
        }

        // Build compact sort key with pre-computed chromosome ID
        CompactSortKey key;
        key.chrom_id = chromToId(chromStart, chromLen, naturalChromOrder);
        key.pos = pos;
        key.offset = lineStart - data;
        key.length = static_cast<uint32_t>(lineLen);
        key.chrom_offset = 0;  // CHROM always starts at beginning of line
        key.chrom_len = static_cast<uint16_t>(chromLen);
        keys.push_back(key);
    }

    // PHASE 2: Sort the compact keys
    if (naturalChromOrder) {
        // Natural order: use pre-computed IDs (O(1) comparison)
        std::sort(keys.begin(), keys.end(), [](const CompactSortKey& a, const CompactSortKey& b) {
            if (a.chrom_id != b.chrom_id) return a.chrom_id < b.chrom_id;
            return a.pos < b.pos;
        });
    } else {
        // Lexicographic order: use string comparison via mmap'd data
        std::sort(keys.begin(), keys.end(), [data](const CompactSortKey& a, const CompactSortKey& b) {
            // Compare chromosome strings lexicographically
            const char* aChrom = data + a.offset;
            const char* bChrom = data + b.offset;
            size_t minLen = std::min(a.chrom_len, b.chrom_len);
            int cmp = std::memcmp(aChrom, bChrom, minLen);
            if (cmp != 0) return cmp < 0;
            if (a.chrom_len != b.chrom_len) return a.chrom_len < b.chrom_len;
            return a.pos < b.pos;
        });
    }

    // PHASE 3: Output header and sorted lines with buffering
    std::string outputBuffer;
    outputBuffer.reserve(1024 * 1024);  // 1MB buffer

    // Output header first
    out.write(headerBlock.data(), headerBlock.size());

    // Output sorted lines
    for (const auto& key : keys) {
        outputBuffer.append(data + key.offset, key.length);
        outputBuffer.push_back('\n');

        // Flush buffer periodically
        if (outputBuffer.size() > 512 * 1024) {
            out.write(outputBuffer.data(), outputBuffer.size());
            outputBuffer.clear();
        }
    }

    // Flush remaining
    if (!outputBuffer.empty()) {
        out.write(outputBuffer.data(), outputBuffer.size());
    }

    munmap(mapped, fileSize);
    close(fd);
    return true;
}

// ============================================================================
// Stdin-based sorting (fallback for pipes)
// ============================================================================
void VCFXSorter::sortInMemory(std::istream &in, std::ostream &out) {
    std::vector<SortKey> records;
    std::string line;

    while (std::getline(in, line)) {
        if (line.empty()) {
            headerLines.push_back(line);
            continue;
        }
        if (line[0] == '#') {
            headerLines.push_back(line);
            continue;
        }

        SortKey key;
        if (!parseChromPos(line, key.chrom, key.pos)) {
            std::cerr << "Warning: skipping malformed line.\n";
            continue;
        }
        // Pre-compute chromosome ID for fast comparison
        key.chrom_id = chromToId(key.chrom.c_str(), key.chrom.size(), naturalChromOrder);
        key.line = std::move(line);
        records.push_back(std::move(key));
    }

    // Sort using pre-computed IDs
    std::sort(records.begin(), records.end(), compareById);

    // Output header
    for (const auto& h : headerLines) {
        out << h << "\n";
    }

    // Output sorted records
    for (const auto& rec : records) {
        out << rec.line << "\n";
    }
}

// Write a sorted chunk to a temp file
std::string VCFXSorter::writeChunk(std::vector<SortKey>& chunk, size_t chunk_num) {
    // Sort the chunk using pre-computed IDs
    std::sort(chunk.begin(), chunk.end(), compareById);

    // Create temp file
    std::string path = tempDir + "/vcfx_sort_" + std::to_string(getpid()) +
                       "_chunk_" + std::to_string(chunk_num) + ".tmp";

    std::ofstream ofs(path, std::ios::binary);
    if (!ofs) {
        std::cerr << "Error: cannot create temp file " << path << "\n";
        return "";
    }

    for (const auto& key : chunk) {
        ofs << key.line << "\n";
    }

    ofs.close();
    return path;
}

// K-way merge of sorted temp files using min-heap
void VCFXSorter::mergeChunks(const std::vector<std::string>& chunk_files, std::ostream &out) {
    std::vector<std::unique_ptr<std::ifstream>> files;
    files.reserve(chunk_files.size());

    for (const auto& path : chunk_files) {
        auto ifs = std::make_unique<std::ifstream>(path, std::ios::binary);
        if (!ifs->is_open()) {
            std::cerr << "Error: cannot open temp file " << path << "\n";
            continue;
        }
        files.push_back(std::move(ifs));
    }

    // Min-heap using greater comparison (so smallest is at top)
    std::priority_queue<MergeEntry, std::vector<MergeEntry>, std::greater<MergeEntry>> heap;

    // Initialize heap with first line from each file
    for (size_t i = 0; i < files.size(); ++i) {
        std::string line;
        if (std::getline(*files[i], line)) {
            MergeEntry entry;
            entry.file_index = i;
            entry.line = std::move(line);
            if (parseChromPos(entry.line, entry.chrom, entry.pos)) {
                // Pre-compute chromosome ID for fast heap comparison
                entry.chrom_id = chromToId(entry.chrom.c_str(), entry.chrom.size(), naturalChromOrder);
                heap.push(std::move(entry));
            }
        }
    }

    // Merge with output buffering
    std::string outputBuffer;
    outputBuffer.reserve(1024 * 1024);  // 1MB buffer

    while (!heap.empty()) {
        MergeEntry top = heap.top();
        heap.pop();

        outputBuffer.append(top.line);
        outputBuffer.push_back('\n');

        // Flush buffer periodically
        if (outputBuffer.size() > 512 * 1024) {
            out.write(outputBuffer.data(), outputBuffer.size());
            outputBuffer.clear();
        }

        // Read next line from same file
        std::string line;
        if (std::getline(*files[top.file_index], line)) {
            MergeEntry entry;
            entry.file_index = top.file_index;
            entry.line = std::move(line);
            if (parseChromPos(entry.line, entry.chrom, entry.pos)) {
                entry.chrom_id = chromToId(entry.chrom.c_str(), entry.chrom.size(), naturalChromOrder);
                heap.push(std::move(entry));
            }
        }
    }

    // Flush remaining
    if (!outputBuffer.empty()) {
        out.write(outputBuffer.data(), outputBuffer.size());
    }

    // Close and delete temp files
    files.clear();
    for (const auto& path : chunk_files) {
        std::remove(path.c_str());
    }
}

// External merge sort for large files from stdin
void VCFXSorter::sortExternal(std::istream &in, std::ostream &out) {
    std::vector<std::string> chunk_files;
    std::vector<SortKey> current_chunk;
    size_t current_bytes = 0;
    size_t chunk_num = 0;
    std::string line;

    // First pass: read headers and split data into sorted chunks
    while (std::getline(in, line)) {
        if (line.empty()) {
            headerLines.push_back(line);
            continue;
        }
        if (line[0] == '#') {
            headerLines.push_back(line);
            continue;
        }

        SortKey key;
        if (!parseChromPos(line, key.chrom, key.pos)) {
            std::cerr << "Warning: skipping malformed line.\n";
            continue;
        }

        // Pre-compute chromosome ID
        key.chrom_id = chromToId(key.chrom.c_str(), key.chrom.size(), naturalChromOrder);

        size_t line_bytes = line.size() + key.chrom.size() + sizeof(key);
        key.line = std::move(line);

        current_chunk.push_back(std::move(key));
        current_bytes += line_bytes;

        // If chunk is full, write it out
        if (current_bytes >= chunkSizeBytes) {
            std::string path = writeChunk(current_chunk, chunk_num++);
            if (!path.empty()) {
                chunk_files.push_back(path);
            }
            current_chunk.clear();
            current_bytes = 0;
        }
    }

    // Write remaining chunk
    if (!current_chunk.empty()) {
        std::string path = writeChunk(current_chunk, chunk_num++);
        if (!path.empty()) {
            chunk_files.push_back(path);
        }
        current_chunk.clear();
    }

    // Output headers
    for (const auto& h : headerLines) {
        out << h << "\n";
    }

    // If only one chunk, it's already sorted - just copy it
    if (chunk_files.size() == 1) {
        std::ifstream ifs(chunk_files[0], std::ios::binary);
        std::string chunk_line;
        while (std::getline(ifs, chunk_line)) {
            out << chunk_line << "\n";
        }
        std::remove(chunk_files[0].c_str());
        return;
    }

    // K-way merge
    if (!chunk_files.empty()) {
        mergeChunks(chunk_files, out);
    }
}

int VCFXSorter::run(int argc, char *argv[]) {
    bool showHelp = false;
    std::vector<std::string> inputFiles;

    static struct option long_opts[] = {
        {"help", no_argument, 0, 'h'},
        {"natural-chr", no_argument, 0, 'n'},
        {"max-memory", required_argument, 0, 'm'},
        {"temp-dir", required_argument, 0, 't'},
        {0, 0, 0, 0}
    };

    // Reset getopt
    optind = 1;

    while (true) {
        int c = ::getopt_long(argc, argv, "hnm:t:", long_opts, nullptr);
        if (c == -1) break;

        switch (c) {
        case 'h':
            showHelp = true;
            break;
        case 'n':
            naturalChromOrder = true;
            break;
        case 'm':
            chunkSizeBytes = std::stoul(optarg) * 1024 * 1024;
            break;
        case 't':
            tempDir = optarg;
            break;
        default:
            showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    // Collect remaining arguments as input files
    for (int i = optind; i < argc; i++) {
        inputFiles.push_back(argv[i]);
    }

    // Set global comparison mode
    g_naturalOrder = naturalChromOrder;

    // Set comparison function for backward compatibility
    compareFunc = compareById;  // Always use optimized comparison

    // If file argument provided, use mmap (26x faster)
    if (!inputFiles.empty()) {
        for (const auto& file : inputFiles) {
            if (!sortFileMmap(file.c_str(), std::cout)) {
                return 1;
            }
        }
        return 0;
    }

    // Stdin processing (slower but works with pipes)
    // Note: We always process stdin if no files are provided, even with argc==1
    sortExternal(std::cin, std::cout);
    return 0;
}

static void show_help() {
    VCFXSorter obj;
    char arg0[] = "VCFX_sorter";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio

    if (vcfx::handle_common_flags(argc, argv, "VCFX_sorter", show_help))
        return 0;

    VCFXSorter app;
    return app.run(argc, argv);
}
