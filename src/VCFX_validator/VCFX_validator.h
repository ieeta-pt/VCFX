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
    int threadCount = 1;  // Default single-threaded
    std::string inputFile;  // Input file path (empty or "-" for stdin)

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
    std::unordered_map<std::string, FieldDef> infoDefs;
    std::unordered_map<std::string, FieldDef> formatDefs;
    std::unordered_set<uint64_t> seenVariantHashes;  // Use hash instead of string

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
