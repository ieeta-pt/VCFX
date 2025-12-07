#ifndef VCFX_VALIDATOR_H
#define VCFX_VALIDATOR_H

#include <iostream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <atomic>
#include <thread>
#include <mutex>
#include <functional>

// ============================================================================
// Transparent hash and equality for string_view lookups without allocation
// ============================================================================

struct StringHash {
    using is_transparent = void;
    size_t operator()(std::string_view sv) const noexcept {
        return std::hash<std::string_view>{}(sv);
    }
    size_t operator()(const std::string& s) const noexcept {
        return std::hash<std::string_view>{}(std::string_view(s));
    }
    size_t operator()(const char* s) const noexcept {
        return std::hash<std::string_view>{}(std::string_view(s));
    }
};

struct StringEqual {
    using is_transparent = void;
    bool operator()(std::string_view a, std::string_view b) const noexcept {
        return a == b;
    }
};

// Memory-mapped file wrapper for zero-copy I/O
struct MappedFile {
    const char* data = nullptr;
    size_t size = 0;
    int fd = -1;

    bool open(const char* path);
    void close();
    ~MappedFile();
};

class VCFXValidator {
  public:
    int run(int argc, char *argv[]);

  private:
    // Configuration flags
    bool strictMode = false;
    bool reportDuplicates = false;
    bool skipDuplicateCheck = false;  // --no-dup-check option
    int threadCount = 1;  // Default single-threaded
    std::string inputFile;  // Input file path (empty or "-" for stdin)
    size_t bloomSizeMB = 128;  // Default bloom filter size in MB

    // Number of columns in the #CHROM header line
    int headerColumnCount = 0;
    // Whether the header includes FORMAT/sample columns
    bool headerHasFormat = false;
    // Number of sample columns
    int sampleCount = 0;

    struct FieldDef {
        std::string number;
        std::string type;
        int numericNumber = -1;  // Cached numeric value, -1 if not a number
    };

    // Note: std::unordered_map heterogeneous lookup requires C++20 on macOS
    // Using standard maps for portability
    std::unordered_map<std::string, FieldDef> infoDefs;
    std::unordered_map<std::string, FieldDef> formatDefs;

    // Bloom filter for memory-efficient duplicate detection
    std::vector<uint64_t> bloomFilter;
    size_t bloomBitCount = 0;  // Total bits in bloom filter

    // Bloom filter operations
    void initBloomFilter(size_t sizeMB);
    inline void bloomAdd(uint64_t hash);
    inline bool bloomMayContain(uint64_t hash) const;

    // Reusable buffers to avoid allocations
    std::vector<std::string_view> fieldBuffer;
    std::vector<std::string_view> sampleBuffer;
    std::vector<std::string_view> formatPartsBuffer;
    std::vector<std::string_view> infoTokenBuffer;

    // FORMAT field caching (most VCF files use same FORMAT for all lines)
    std::string cachedFormatStr;  // Stored as string since mmap data may move
    std::vector<std::string> cachedFormatParts;
    int cachedGTIndex = -1;

    // Show usage
    void displayHelp();

    // Main validation entry points
    bool validateVCF(std::istream &in);  // For stdin
    bool validateVCFMmap(const char* filepath);  // For file with mmap

    // Validate a meta line "##" and returns true if it's correct
    bool validateMetaLine(std::string_view line, int lineNumber);

    // Check #CHROM line
    bool validateChromHeader(std::string_view line, int lineNumber);

    // Validate a data line with at least 8 columns
    bool validateDataLine(std::string_view line, int lineNumber);

    // Parse FORMAT field and cache it
    void cacheFormat(std::string_view format);

    // Fast inline helpers
    static inline void splitInto(std::string_view sv, char delim, std::vector<std::string_view> &out);
    static inline std::string_view trimView(std::string_view sv);
    static inline bool isValidDNA(std::string_view sv);
    static inline bool isValidGenotype(std::string_view sv);
    static inline uint64_t hashVariant(std::string_view chrom, std::string_view pos,
                                       std::string_view ref, std::string_view alt);
    static inline bool parsePositiveInt(std::string_view sv, int &out);
    static inline bool parseNonNegativeDouble(std::string_view sv);
    static inline size_t countChar(std::string_view sv, char c);
};

#endif
