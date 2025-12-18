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

// ============================================================================
// GATK-compatible validation types (can be individually enabled/disabled)
// ============================================================================
enum class ValidationType {
    REF,        // Check REF matches reference FASTA
    IDS,        // Check IDs exist in dbSNP
    ALLELES,    // Check ALT alleles observed in genotypes
    CHR_COUNTS, // Check AN/AC consistency
    SORTING,    // Check records are sorted
    GVCF        // GVCF-specific validation
};

class VCFXValidator {
  public:
    int run(int argc, char *argv[]);

  private:
    // Configuration flags
    bool strictMode = false;
    bool reportDuplicates = false;
    bool skipDuplicateCheck = false;
    bool allowEmpty = false;
    int threadCount = 1;
    std::string inputFile;
    size_t bloomSizeMB = 128;

    // GATK-compatible validation options
    std::string referenceFile;  // Path to reference FASTA for REF validation
    std::string dbsnpFile;      // Path to dbSNP VCF for ID validation
    bool validateRef = false;   // --validate-ref
    bool validateIds = false;   // --validate-ids (requires --dbsnp)
    bool validateChrCounts = true;  // --validate-chr-counts (default on in strict)
    bool validateSorting = true;    // --validate-sorting (default on)
    bool validateGvcf = false;      // --gvcf mode

    // Number of columns in the #CHROM header line
    int headerColumnCount = 0;
    // Whether the header includes FORMAT/sample columns
    bool headerHasFormat = false;
    // Number of sample columns
    int sampleCount = 0;

    struct FieldDef {
        std::string number;
        std::string type;
        int numericNumber = -1;
    };

    std::unordered_map<std::string, FieldDef> infoDefs;
    std::unordered_map<std::string, FieldDef> formatDefs;

    // Bloom filter for duplicate detection
    std::vector<uint64_t> bloomFilter;
    size_t bloomBitCount = 0;
    void initBloomFilter(size_t sizeMB);
    inline void bloomAdd(uint64_t hash);
    inline bool bloomMayContain(uint64_t hash) const;

    // Reference FASTA data (memory-mapped for speed)
    MappedFile refFile;
    std::unordered_map<std::string, std::pair<size_t, size_t>> contigOffsets; // contig -> (offset, length)
    bool loadReference(const char* path);
    std::string_view getRefSequence(std::string_view chrom, int pos, size_t len);

    // dbSNP ID set (bloom filter for memory efficiency)
    std::vector<uint64_t> dbsnpBloomFilter;
    size_t dbsnpBloomBitCount = 0;
    bool loadDbsnp(const char* path);
    bool isKnownId(std::string_view id) const;

    // Sorting validation state
    std::string lastChrom;
    int lastPos = 0;
    std::unordered_map<std::string, int> chromOrder; // contig -> order index

    // Reusable buffers to avoid allocations
    std::vector<std::string_view> fieldBuffer;
    std::vector<std::string_view> sampleBuffer;
    std::vector<std::string_view> formatPartsBuffer;
    std::vector<std::string_view> infoTokenBuffer;
    std::vector<bool> altAlleleObserved;
    std::vector<int> alleleIndicesBuffer;
    std::vector<int> alleleCountBuffer;  // For CHR_COUNTS validation

    // FORMAT field caching
    std::string cachedFormatStr;
    std::vector<std::string> cachedFormatParts;
    int cachedGTIndex = -1;

    // GVCF state
    int lastGvcfEnd = 0;  // Track coverage for GVCF validation

    // Show usage
    void displayHelp();

    // Main validation entry points
    bool validateVCF(std::istream &in);
    bool validateVCFMmap(const char* filepath);

    // Validation methods
    bool validateMetaLine(std::string_view line, int lineNumber);
    bool validateChromHeader(std::string_view line, int lineNumber);
    bool validateDataLine(std::string_view line, int lineNumber);

    // GATK-compatible validation methods
    bool validateRefBase(std::string_view chrom, int pos, std::string_view ref, int lineNumber);
    bool validateVariantId(std::string_view id, int lineNumber);
    bool validateChrCountsField(std::string_view info, size_t altCount, int lineNumber);
    bool validateSortOrder(std::string_view chrom, int pos, int lineNumber);
    bool validateGvcfRecord(std::string_view chrom, int pos, std::string_view alt,
                            std::string_view info, int lineNumber);

    // Fast inline helpers
    static inline void splitInto(std::string_view sv, char delim, std::vector<std::string_view> &out);
    static inline std::string_view trimView(std::string_view sv);
    static inline bool isValidDNA(std::string_view sv);
    static inline bool isValidGenotype(std::string_view sv);
    static inline uint64_t hashVariant(std::string_view chrom, std::string_view pos,
                                       std::string_view ref, std::string_view alt);
    static inline bool parsePositiveInt(std::string_view sv, int &out);
    static inline bool parseNonNegativeInt(std::string_view sv, int &out);
    static inline bool parseNonNegativeDouble(std::string_view sv);
    static inline size_t countChar(std::string_view sv, char c);
    static inline void extractAlleleIndices(std::string_view gt, std::vector<int> &indices);
    static inline uint64_t hashString(std::string_view sv);
};

#endif
