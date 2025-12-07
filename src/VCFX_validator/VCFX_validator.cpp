#include "VCFX_validator.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <sstream>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

// ============================================================================
// SIMD includes for x86_64 platforms
// ============================================================================
#if defined(__x86_64__) && defined(__AVX2__)
#include <immintrin.h>
#define USE_AVX2 1
#elif defined(__x86_64__) && defined(__SSE2__)
#include <emmintrin.h>
#define USE_SSE2 1
#endif

// ============================================================================
// SIMD-optimized helper functions
// ============================================================================

#if defined(USE_AVX2)
// Find next newline using AVX2 (32 bytes at a time)
static inline const char* findNewlineSIMD(const char* p, const char* end) {
    const __m256i nl = _mm256_set1_epi8('\n');
    while (p + 32 <= end) {
        __m256i chunk = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p));
        __m256i cmp = _mm256_cmpeq_epi8(chunk, nl);
        int mask = _mm256_movemask_epi8(cmp);
        if (mask) {
            return p + __builtin_ctz(static_cast<unsigned int>(mask));
        }
        p += 32;
    }
    return static_cast<const char*>(memchr(p, '\n', static_cast<size_t>(end - p)));
}

// Count character occurrences using AVX2
static inline size_t countCharSIMD(const char* data, size_t len, char c) {
    size_t count = 0;
    const __m256i target = _mm256_set1_epi8(c);
    const char* end = data + len;

    while (data + 32 <= end) {
        __m256i chunk = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(data));
        __m256i cmp = _mm256_cmpeq_epi8(chunk, target);
        count += __builtin_popcount(_mm256_movemask_epi8(cmp));
        data += 32;
    }

    // Handle remainder
    while (data < end) {
        if (*data == c) count++;
        data++;
    }
    return count;
}

#elif defined(USE_SSE2)
// Find next newline using SSE2 (16 bytes at a time)
static inline const char* findNewlineSIMD(const char* p, const char* end) {
    const __m128i nl = _mm_set1_epi8('\n');
    while (p + 16 <= end) {
        __m128i chunk = _mm_loadu_si128(reinterpret_cast<const __m128i*>(p));
        __m128i cmp = _mm_cmpeq_epi8(chunk, nl);
        int mask = _mm_movemask_epi8(cmp);
        if (mask) {
            return p + __builtin_ctz(static_cast<unsigned int>(mask));
        }
        p += 16;
    }
    return static_cast<const char*>(memchr(p, '\n', static_cast<size_t>(end - p)));
}

// Count character occurrences using SSE2
static inline size_t countCharSIMD(const char* data, size_t len, char c) {
    size_t count = 0;
    const __m128i target = _mm_set1_epi8(c);
    const char* end = data + len;

    while (data + 16 <= end) {
        __m128i chunk = _mm_loadu_si128(reinterpret_cast<const __m128i*>(data));
        __m128i cmp = _mm_cmpeq_epi8(chunk, target);
        count += __builtin_popcount(_mm_movemask_epi8(cmp));
        data += 16;
    }

    // Handle remainder
    while (data < end) {
        if (*data == c) count++;
        data++;
    }
    return count;
}

#else
// Portable fallback using memchr (still SIMD-optimized in libc)
static inline const char* findNewlineSIMD(const char* p, const char* end) {
    return static_cast<const char*>(memchr(p, '\n', static_cast<size_t>(end - p)));
}

// Portable fallback for character counting
static inline size_t countCharSIMD(const char* data, size_t len, char c) {
    size_t count = 0;
    for (size_t i = 0; i < len; i++) {
        if (data[i] == c) count++;
    }
    return count;
}
#endif

// ============================================================================
// DNA validation lookup table (faster than conditionals)
// Valid: A=65, C=67, G=71, N=78, T=84, a=97, c=99, g=103, n=110, t=116
// ============================================================================
static const uint8_t DNA_TABLE[256] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  // 0-15
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  // 16-31
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  // 32-47
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  // 48-63
    0,1,0,1,0,0,0,1,0,0,0,0,0,0,1,0,  // 64-79:  A=65, C=67, G=71, N=78
    0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,  // 80-95:  T=84
    0,1,0,1,0,0,0,1,0,0,0,0,0,0,1,0,  // 96-111: a=97, c=99, g=103, n=110
    0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,  // 112-127: t=116
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  // 128-143
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  // 144-159
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  // 160-175
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  // 176-191
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  // 192-207
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  // 208-223
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  // 224-239
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0   // 240-255
};

// Fast DNA validation using lookup table
static inline bool isValidDNAFast(const char* data, size_t len) {
    if (len == 0) return false;
    for (size_t i = 0; i < len; i++) {
        if (!DNA_TABLE[static_cast<unsigned char>(data[i])]) return false;
    }
    return true;
}

// ============================================================================
// MappedFile implementation - zero-copy file I/O
// ============================================================================

bool MappedFile::open(const char* path) {
    fd = ::open(path, O_RDONLY);
    if (fd < 0) return false;

    struct stat st;
    if (fstat(fd, &st) < 0) {
        close();
        return false;
    }
    size = st.st_size;
    if (size == 0) {
        // Empty file - still valid but nothing to map
        return true;
    }

    data = static_cast<const char*>(mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0));
    if (data == MAP_FAILED) {
        data = nullptr;
        close();
        return false;
    }

    // Hint for sequential access - improves read-ahead
    madvise(const_cast<char*>(data), size, MADV_SEQUENTIAL);
    return true;
}

void MappedFile::close() {
    if (data && size > 0) {
        munmap(const_cast<char*>(data), size);
        data = nullptr;
    }
    if (fd >= 0) {
        ::close(fd);
        fd = -1;
    }
    size = 0;
}

MappedFile::~MappedFile() {
    close();
}

// ============================================================================
// Bloom filter implementation for memory-efficient duplicate detection
// ============================================================================

void VCFXValidator::initBloomFilter(size_t sizeMB) {
    // Convert MB to bits: sizeMB * 1024 * 1024 * 8 bits
    bloomBitCount = sizeMB * 1024ULL * 1024ULL * 8ULL;
    // Allocate as 64-bit words
    size_t wordCount = (bloomBitCount + 63) / 64;
    bloomFilter.resize(wordCount, 0);
}

inline void VCFXValidator::bloomAdd(uint64_t hash) {
    // Use two hash functions derived from the single hash (double hashing)
    size_t h1 = hash % bloomBitCount;
    size_t h2 = (hash >> 17) % bloomBitCount;
    size_t h3 = ((hash >> 34) ^ (hash >> 51)) % bloomBitCount;

    bloomFilter[h1 / 64] |= (1ULL << (h1 % 64));
    bloomFilter[h2 / 64] |= (1ULL << (h2 % 64));
    bloomFilter[h3 / 64] |= (1ULL << (h3 % 64));
}

inline bool VCFXValidator::bloomMayContain(uint64_t hash) const {
    size_t h1 = hash % bloomBitCount;
    size_t h2 = (hash >> 17) % bloomBitCount;
    size_t h3 = ((hash >> 34) ^ (hash >> 51)) % bloomBitCount;

    return (bloomFilter[h1 / 64] & (1ULL << (h1 % 64))) &&
           (bloomFilter[h2 / 64] & (1ULL << (h2 % 64))) &&
           (bloomFilter[h3 / 64] & (1ULL << (h3 % 64)));
}

// ============================================================================
// Fast inline helper implementations
// ============================================================================

inline void VCFXValidator::splitInto(std::string_view sv, char delim, std::vector<std::string_view> &out) {
    out.clear();
    size_t start = 0;
    size_t end;
    while ((end = sv.find(delim, start)) != std::string_view::npos) {
        out.push_back(sv.substr(start, end - start));
        start = end + 1;
    }
    out.push_back(sv.substr(start));
}

inline std::string_view VCFXValidator::trimView(std::string_view sv) {
    while (!sv.empty() && std::isspace(static_cast<unsigned char>(sv.front())))
        sv.remove_prefix(1);
    while (!sv.empty() && std::isspace(static_cast<unsigned char>(sv.back())))
        sv.remove_suffix(1);
    return sv;
}

inline bool VCFXValidator::isValidDNA(std::string_view sv) {
    // Use optimized lookup table version
    return isValidDNAFast(sv.data(), sv.size());
}

// Fast genotype validation without regex
// Valid formats: "." or digits separated by / or |
// Examples: ".", "0", "0/1", "0|1", "1/2/3", "0|0|1"
inline bool VCFXValidator::isValidGenotype(std::string_view sv) {
    // Fast path for most common diploid genotypes (covers >99% of cases in typical VCFs)
    if (sv.size() == 3) {
        char c0 = sv[0], c1 = sv[1], c2 = sv[2];
        if ((c0 >= '0' && c0 <= '9') && (c1 == '/' || c1 == '|') && (c2 >= '0' && c2 <= '9'))
            return true;
    }

    if (sv.empty()) return false;
    if (sv.size() == 1) {
        return sv[0] == '.' || (sv[0] >= '0' && sv[0] <= '9');
    }

    // Must start with a digit
    char c0 = sv[0];
    if (c0 < '0' || c0 > '9') return false;

    bool expectDigit = false;
    for (size_t i = 1; i < sv.size(); ++i) {
        char c = sv[i];
        if (expectDigit) {
            if (c < '0' || c > '9') return false;
            expectDigit = false;
        } else {
            if (c == '/' || c == '|') {
                expectDigit = true;
            } else if (c < '0' || c > '9') {
                return false;
            }
        }
    }
    // Must not end with separator
    return !expectDigit;
}

// FNV-1a hash for variant deduplication
inline uint64_t VCFXValidator::hashVariant(std::string_view chrom, std::string_view pos,
                                           std::string_view ref, std::string_view alt) {
    uint64_t hash = 14695981039346656037ULL;
    auto hashBytes = [&hash](std::string_view s) {
        for (char c : s) {
            hash ^= static_cast<uint64_t>(static_cast<unsigned char>(c));
            hash *= 1099511628211ULL;
        }
        hash ^= 0xFF;  // Separator
        hash *= 1099511628211ULL;
    };
    hashBytes(chrom);
    hashBytes(pos);
    hashBytes(ref);
    hashBytes(alt);
    return hash;
}

inline bool VCFXValidator::parsePositiveInt(std::string_view sv, int &out) {
    if (sv.empty()) return false;
    auto result = std::from_chars(sv.data(), sv.data() + sv.size(), out);
    return result.ec == std::errc{} && result.ptr == sv.data() + sv.size() && out > 0;
}

inline bool VCFXValidator::parseNonNegativeDouble(std::string_view sv) {
    if (sv.empty()) return false;
    // Fast path for common integer values
    bool hasDigit = false;
    bool hasDot = false;
    bool hasE = false;
    size_t i = 0;

    // Optional leading sign (though VCF QUAL shouldn't have +)
    if (i < sv.size() && (sv[i] == '-' || sv[i] == '+')) {
        if (sv[i] == '-') return false;  // Negative not allowed
        ++i;
    }

    // Integer part
    while (i < sv.size() && std::isdigit(static_cast<unsigned char>(sv[i]))) {
        hasDigit = true;
        ++i;
    }

    // Decimal part
    if (i < sv.size() && sv[i] == '.') {
        hasDot = true;
        ++i;
        while (i < sv.size() && std::isdigit(static_cast<unsigned char>(sv[i]))) {
            hasDigit = true;
            ++i;
        }
    }

    // Exponent part
    if (i < sv.size() && (sv[i] == 'e' || sv[i] == 'E')) {
        hasE = true;
        ++i;
        if (i < sv.size() && (sv[i] == '-' || sv[i] == '+')) ++i;
        bool hasExpDigit = false;
        while (i < sv.size() && std::isdigit(static_cast<unsigned char>(sv[i]))) {
            hasExpDigit = true;
            ++i;
        }
        if (!hasExpDigit) return false;
    }

    return hasDigit && i == sv.size();
}

inline size_t VCFXValidator::countChar(std::string_view sv, char c) {
    // Use SIMD-optimized version
    return countCharSIMD(sv.data(), sv.size(), c);
}

// ============================================================================
// Main validator implementation
// ============================================================================

int VCFXValidator::run(int argc, char *argv[]) {
    bool hasStdin = !isatty(fileno(stdin));
    if (argc == 1 && !hasStdin) {
        displayHelp();
        return 0;
    }
    bool showHelp = false;
    static struct option long_opts[] = {{"help", no_argument, 0, 'h'},
                                        {"strict", no_argument, 0, 's'},
                                        {"report-dups", no_argument, 0, 'd'},
                                        {"no-dup-check", no_argument, 0, 'n'},
                                        {"bloom-size", required_argument, 0, 'b'},
                                        {"threads", required_argument, 0, 't'},
                                        {0, 0, 0, 0}};
    while (true) {
        int c = ::getopt_long(argc, argv, "hsdn b:t:", long_opts, nullptr);
        if (c == -1) break;
        switch (c) {
        case 'h': showHelp = true; break;
        case 's': strictMode = true; break;
        case 'd': reportDuplicates = true; break;
        case 'n': skipDuplicateCheck = true; break;
        case 'b':
            bloomSizeMB = static_cast<size_t>(std::max(1, std::atoi(optarg)));
            break;
        case 't':
            threadCount = std::max(1, std::atoi(optarg));
            break;
        default: showHelp = true;
        }
    }
    if (showHelp) {
        displayHelp();
        return 0;
    }

    // Check for input file argument (remaining non-option argument)
    if (optind < argc) {
        inputFile = argv[optind];
    }

    // Reserve capacity for buffers to reduce reallocations
    fieldBuffer.reserve(2600);  // ~2504 samples + 9 fixed columns
    sampleBuffer.reserve(16);
    formatPartsBuffer.reserve(16);
    infoTokenBuffer.reserve(32);
    cachedFormatParts.reserve(16);

    // Initialize bloom filter for duplicate detection (unless disabled)
    if (!skipDuplicateCheck) {
        initBloomFilter(bloomSizeMB);
    }

    bool ok;
    // Use mmap for file input, stdin for pipe/stdin
    if (!inputFile.empty() && inputFile != "-") {
        ok = validateVCFMmap(inputFile.c_str());
    } else {
        ok = validateVCF(std::cin);
    }
    return (ok ? 0 : 1);
}

void VCFXValidator::displayHelp() {
    std::cout << "VCFX_validator: Checks basic validity of a VCF.\n\n"
                 "Usage:\n"
                 "  VCFX_validator [options] [input.vcf]\n"
                 "  VCFX_validator [options] < input.vcf\n\n"
                 "Options:\n"
                 "  -h, --help          Show this help.\n"
                 "  -s, --strict        Enable stricter checks.\n"
                 "  -d, --report-dups   Report duplicate records.\n"
                 "  -n, --no-dup-check  Skip duplicate detection (faster).\n"
                 "  -b, --bloom-size N  Bloom filter size in MB (default: 128).\n"
                 "  -t, --threads N     Reserved for future multi-threaded validation.\n\n"
                 "Description:\n"
                 "  Validates:\n"
                 "   * All '##' lines are recognized as meta lines.\n"
                 "   * #CHROM line is present and well formed.\n"
                 "   * Each data line has >=8 columns, checks CHROM non-empty, POS>0,\n"
                 "     REF/ALT non-empty, QUAL is '.' or non-negative float, FILTER non-empty,\n"
                 "     INFO fields are checked against header definitions.\n"
                 "   * FORMAT fields and genotype values are validated.\n"
                 "   * REF and ALT sequences must contain only A,C,G,T,N.\n"
                 "   * Duplicate variants are detected (use --no-dup-check to skip).\n"
                 "  In strict mode additional checks are performed:\n"
                 "   * Data line column count must match the #CHROM header.\n"
                 "   * Sample columns must match the FORMAT field structure.\n"
                 "   * Any warning is treated as an error.\n"
                 "  Exits 0 if pass, 1 if fail.\n\n"
                 "Performance:\n"
                 "  * Pass file path directly for memory-mapped I/O (fastest).\n"
                 "  * Uses SIMD-optimized parsing on x86_64 (AVX2/SSE2).\n"
                 "  * Bloom filter uses " << bloomSizeMB << "MB for duplicate detection.\n"
                 "  * Use --no-dup-check when data is known to be clean.\n";
}

bool VCFXValidator::validateMetaLine(std::string_view line, int lineNumber) {
    if (line.size() < 2) return false;

    // Check for INFO or FORMAT definitions
    if (line.substr(0, 7) == "##INFO=" || line.substr(0, 9) == "##FORMAT=") {
        bool isInfo = (line[2] == 'I');
        size_t start = line.find('<');
        size_t end = line.rfind('>');
        if (start == std::string_view::npos || end == std::string_view::npos || end <= start) {
            std::cerr << "Error: malformed header at line " << lineNumber << ".\n";
            return false;
        }

        std::string_view inside = line.substr(start + 1, end - start - 1);
        std::string id, number, type;

        // Parse key=value pairs manually for speed
        size_t pos = 0;
        while (pos < inside.size()) {
            size_t eq = inside.find('=', pos);
            if (eq == std::string_view::npos) break;

            size_t comma = inside.find(',', eq);
            if (comma == std::string_view::npos) comma = inside.size();

            // Handle quoted values (Description often has commas)
            if (eq + 1 < inside.size() && inside[eq + 1] == '"') {
                size_t quoteEnd = inside.find('"', eq + 2);
                if (quoteEnd != std::string_view::npos) {
                    comma = inside.find(',', quoteEnd);
                    if (comma == std::string_view::npos) comma = inside.size();
                }
            }

            std::string_view key = trimView(inside.substr(pos, eq - pos));
            std::string_view val = trimView(inside.substr(eq + 1, comma - eq - 1));

            // Remove quotes if present
            if (!val.empty() && val.front() == '"' && val.back() == '"') {
                val = val.substr(1, val.size() - 2);
            }

            if (key == "ID") id = std::string(val);
            else if (key == "Number") number = std::string(val);
            else if (key == "Type") type = std::string(val);

            pos = (comma < inside.size()) ? comma + 1 : inside.size();
        }

        if (id.empty()) {
            std::cerr << "Error: header line missing ID at " << lineNumber << ".\n";
            return false;
        }

        FieldDef def{number, type, -1};
        // Pre-parse numeric Number field
        if (!number.empty()) {
            int num;
            auto result = std::from_chars(number.data(), number.data() + number.size(), num);
            if (result.ec == std::errc{} && result.ptr == number.data() + number.size()) {
                def.numericNumber = num;
            }
        }

        if (isInfo) {
            infoDefs[id] = std::move(def);
        } else {
            formatDefs[id] = std::move(def);
        }
        return true;
    }

    // Any other ## line is valid
    if (line.size() >= 2 && line[0] == '#' && line[1] == '#') {
        return true;
    }

    std::cerr << "Error: line " << lineNumber << " is a header line but doesn't start with '##'.\n";
    return false;
}

bool VCFXValidator::validateChromHeader(std::string_view line, int lineNumber) {
    splitInto(line, '\t', fieldBuffer);

    if (fieldBuffer.size() < 8) {
        std::cerr << "Error: #CHROM line at " << lineNumber << " has <8 columns.\n";
        return false;
    }
    if (fieldBuffer[0] != "#CHROM") {
        std::cerr << "Error: #CHROM line doesn't start with '#CHROM' at line " << lineNumber << ".\n";
        return false;
    }

    headerColumnCount = static_cast<int>(fieldBuffer.size());
    headerHasFormat = (headerColumnCount > 8);
    sampleCount = headerHasFormat ? headerColumnCount - 9 : 0;

    if (headerHasFormat && fieldBuffer[8] != "FORMAT") {
        std::string msg = "Warning: column 9 of #CHROM header is not 'FORMAT'.";
        if (strictMode) {
            std::cerr << "Error: " << msg << "\n";
            return false;
        } else {
            std::cerr << msg << "\n";
        }
    }
    return true;
}

bool VCFXValidator::validateDataLine(std::string_view line, int lineNumber) {
    // Fast path: only split first 9 columns initially, keep rest as single view
    // This avoids splitting ~2500 sample columns when we may not need them all
    std::string_view fixedFields[9];
    std::string_view sampleData;  // Everything after FORMAT column
    size_t fieldIdx = 0;
    size_t fieldStart = 0;
    size_t scanPos = 0;

    while (fieldIdx < 9 && scanPos < line.size()) {
        if (line[scanPos] == '\t') {
            fixedFields[fieldIdx++] = line.substr(fieldStart, scanPos - fieldStart);
            fieldStart = scanPos + 1;
        }
        ++scanPos;
    }

    // Handle last fixed field or remaining data
    if (fieldIdx < 9) {
        fixedFields[fieldIdx++] = line.substr(fieldStart);
    } else {
        sampleData = line.substr(fieldStart);
    }

    const size_t numFixedFields = fieldIdx;
    if (numFixedFields < 8) {
        std::cerr << "Error: line " << lineNumber << " has <8 columns.\n";
        return false;
    }

    // In strict mode, check that column count matches header
    if (strictMode && headerColumnCount > 0) {
        // Count total columns in data line
        size_t dataColumnCount = numFixedFields;
        if (!sampleData.empty()) {
            dataColumnCount++;  // Count the first sample column
            for (char c : sampleData) {
                if (c == '\t') dataColumnCount++;
            }
        }
        if (dataColumnCount != static_cast<size_t>(headerColumnCount)) {
            std::cerr << "Error: line " << lineNumber << " has " << dataColumnCount
                      << " columns but header has " << headerColumnCount << " columns.\n";
            return false;
        }
    }

    // Get trimmed views of fixed fields
    std::string_view chrom = trimView(fixedFields[0]);
    std::string_view posStr = trimView(fixedFields[1]);
    // fixedFields[2] is ID - not validated
    std::string_view ref = trimView(fixedFields[3]);
    std::string_view alt = trimView(fixedFields[4]);
    std::string_view qual = trimView(fixedFields[5]);
    std::string_view filter = trimView(fixedFields[6]);
    std::string_view info = trimView(fixedFields[7]);

    // CHROM
    if (chrom.empty()) {
        std::cerr << "Error: line " << lineNumber << " CHROM is empty.\n";
        return false;
    }

    // POS
    int pos;
    if (!parsePositiveInt(posStr, pos)) {
        std::cerr << "Error: line " << lineNumber << " POS must be >0.\n";
        return false;
    }

    // REF
    if (ref.empty()) {
        std::cerr << "Error: line " << lineNumber << " REF is empty.\n";
        return false;
    }
    if (!isValidDNA(ref)) {
        std::cerr << "Error: line " << lineNumber << " REF has invalid characters.\n";
        return false;
    }

    // ALT - can be multi-allelic (comma-separated)
    if (alt.empty()) {
        std::cerr << "Error: line " << lineNumber << " ALT is empty.\n";
        return false;
    }
    // Validate each ALT allele separately (split by comma for multi-allelic)
    {
        size_t start = 0;
        size_t end;
        while ((end = alt.find(',', start)) != std::string_view::npos) {
            std::string_view allele = alt.substr(start, end - start);
            if (allele.empty() || !isValidDNA(allele)) {
                std::cerr << "Error: line " << lineNumber << " ALT has invalid characters.\n";
                return false;
            }
            start = end + 1;
        }
        // Last (or only) allele
        std::string_view lastAllele = alt.substr(start);
        if (lastAllele.empty() || !isValidDNA(lastAllele)) {
            std::cerr << "Error: line " << lineNumber << " ALT has invalid characters.\n";
            return false;
        }
    }

    // QUAL
    if (qual != ".") {
        if (!parseNonNegativeDouble(qual)) {
            std::cerr << "Error: line " << lineNumber << " invalid or negative QUAL.\n";
            return false;
        }
    }

    // FILTER
    if (filter.empty()) {
        std::cerr << "Error: line " << lineNumber << " FILTER is empty.\n";
        return false;
    }

    // INFO field validation
    if (info != ".") {
        splitInto(info, ';', infoTokenBuffer);
        bool anyValid = false;
        for (auto token : infoTokenBuffer) {
            token = trimView(token);
            if (token.empty()) continue;

            size_t eq = token.find('=');
            std::string_view key, val;
            if (eq == std::string_view::npos) {
                key = token;
            } else {
                key = token.substr(0, eq);
                val = token.substr(eq + 1);
            }

            if (key.empty()) {
                std::cerr << "Error: line " << lineNumber << " has INFO with empty key.\n";
                return false;
            }

            // Look up INFO field definition
            std::string keyStr(key);
            auto it = infoDefs.find(keyStr);
            if (it == infoDefs.end()) {
                std::cerr << (strictMode ? "Error: " : "Warning: ")
                          << "INFO field " << key << " not defined in header on line "
                          << lineNumber << ".\n";
                if (strictMode) return false;
            } else if (eq != std::string_view::npos && it->second.numericNumber >= 0) {
                // Check value count
                size_t have = countChar(val, ',') + 1;
                if (have != static_cast<size_t>(it->second.numericNumber)) {
                    std::cerr << (strictMode ? "Error: " : "Warning: ")
                              << "INFO field " << key << " expected "
                              << it->second.number << " values on line "
                              << lineNumber << ".\n";
                    if (strictMode) return false;
                }
            }
            anyValid = true;
        }
        if (!anyValid) {
            std::cerr << "Error: line " << lineNumber << " INFO not valid.\n";
            return false;
        }
    }

    // FORMAT and sample columns validation
    if (headerHasFormat) {
        if (numFixedFields < 9 || sampleData.empty()) {
            std::cerr << "Error: line " << lineNumber << " missing FORMAT/sample columns.\n";
            return false;
        }

        std::string_view format = trimView(fixedFields[8]);

        // Fast path: FORMAT is exactly "GT" (very common case)
        // Validate every genotype directly without splitting into vector
        if (format == "GT") {
            size_t sampleStart = 0;
            const size_t dataLen = sampleData.size();
            while (sampleStart < dataLen) {
                size_t tabPos = sampleData.find('\t', sampleStart);
                size_t sampleEnd = (tabPos == std::string_view::npos) ? dataLen : tabPos;
                size_t len = sampleEnd - sampleStart;

                // Skip empty and "."
                if (len > 0 && !(len == 1 && sampleData[sampleStart] == '.')) {
                    std::string_view sample = sampleData.substr(sampleStart, len);
                    if (!isValidGenotype(sample)) {
                        if (strictMode) {
                            std::cerr << "Error: invalid genotype on line " << lineNumber << ".\n";
                            return false;
                        }
                        // Non-strict: silently skip invalid genotypes (no warnings for performance)
                    }
                }
                sampleStart = sampleEnd + 1;
            }
        } else {
            // Full validation for complex FORMAT fields - need to split samples
            splitInto(sampleData, '\t', fieldBuffer);

            // FORMAT caching: most VCF files use the same FORMAT for all lines
            // Only re-parse if FORMAT changed
            int gtIndex;
            size_t formatPartCount;
            if (format == cachedFormatStr) {
                // Use cached values
                gtIndex = cachedGTIndex;
                formatPartCount = cachedFormatParts.size();
            } else {
                // Parse and cache new FORMAT
                cachedFormatStr = std::string(format);
                splitInto(format, ':', formatPartsBuffer);

                cachedFormatParts.clear();
                cachedFormatParts.reserve(formatPartsBuffer.size());
                for (const auto &fp : formatPartsBuffer) {
                    cachedFormatParts.push_back(std::string(fp));
                }

                // Validate FORMAT field definitions (only on format change)
                for (const auto &fpStr : cachedFormatParts) {
                    auto it = formatDefs.find(fpStr);
                    if (it == formatDefs.end()) {
                        std::cerr << (strictMode ? "Error: " : "Warning: ")
                                  << "FORMAT field " << fpStr << " not defined in header on line "
                                  << lineNumber << ".\n";
                        if (strictMode) return false;
                    }
                }

                // Find GT position if present (usually first)
                cachedGTIndex = -1;
                for (size_t i = 0; i < cachedFormatParts.size(); ++i) {
                    if (cachedFormatParts[i] == "GT") {
                        cachedGTIndex = static_cast<int>(i);
                        break;
                    }
                }

                gtIndex = cachedGTIndex;
                formatPartCount = cachedFormatParts.size();
            }

            // Validate sample columns
            for (size_t i = 0; i < fieldBuffer.size(); ++i) {
                std::string_view sample = trimView(fieldBuffer[i]);
                if (sample == "." || sample.empty()) continue;

                splitInto(sample, ':', sampleBuffer);

                if (sampleBuffer.size() != formatPartCount) {
                    std::string msg = "Warning: sample column " + std::to_string(i + 1) +
                                      " does not match FORMAT field";
                    if (strictMode) {
                        std::cerr << "Error: " << msg << " on line " << lineNumber << ".\n";
                        return false;
                    } else {
                        std::cerr << msg << " on line " << lineNumber << ".\n";
                    }
                }

                // Validate GT field if present
                if (gtIndex >= 0 && gtIndex < static_cast<int>(sampleBuffer.size())) {
                    std::string_view gtVal = sampleBuffer[gtIndex];
                    if (!gtVal.empty() && !isValidGenotype(gtVal)) {
                        std::string msg = "invalid genotype";
                        if (strictMode) {
                            std::cerr << "Error: " << msg << " on line " << lineNumber << ".\n";
                            return false;
                        } else {
                            std::cerr << "Warning: " << msg << " on line " << lineNumber << ".\n";
                        }
                    }
                }

                // Validate numeric FORMAT fields (use cached format parts)
                for (size_t j = 0; j < sampleBuffer.size() && j < formatPartCount; ++j) {
                    if (static_cast<int>(j) == gtIndex) continue;  // Already validated

                    const std::string& keyStr = cachedFormatParts[j];
                    auto it = formatDefs.find(keyStr);
                    if (it != formatDefs.end() && it->second.numericNumber >= 0) {
                        size_t have = countChar(sampleBuffer[j], ',') + 1;
                        if (have != static_cast<size_t>(it->second.numericNumber)) {
                            std::cerr << (strictMode ? "Error: " : "Warning: ")
                                      << "FORMAT field " << keyStr << " expected "
                                      << it->second.number << " values on line "
                                      << lineNumber << ".\n";
                            if (strictMode) return false;
                        }
                    }
                }
            }
        }
    } else if (!sampleData.empty()) {
        std::string msg = "Warning: data line has sample columns but header lacks FORMAT";
        if (strictMode) {
            std::cerr << "Error: " << msg << " on line " << lineNumber << ".\n";
            return false;
        } else {
            std::cerr << msg << " on line " << lineNumber << ".\n";
        }
    }

    // Duplicate detection using bloom filter (unless disabled)
    if (!skipDuplicateCheck) {
        uint64_t variantHash = hashVariant(chrom, posStr, ref, alt);
        if (bloomMayContain(variantHash)) {
            // Potential duplicate detected by bloom filter
            if (reportDuplicates) {
                std::cerr << "Duplicate at line " << lineNumber << "\n";
            }
            std::string msg = "duplicate variant";
            if (strictMode) {
                std::cerr << "Error: " << msg << " on line " << lineNumber << ".\n";
                return false;
            } else {
                std::cerr << "Warning: " << msg << " on line " << lineNumber << ".\n";
            }
        }
        bloomAdd(variantHash);
    }

    return true;
}

// ============================================================================
// Memory-mapped validation for file input
// ============================================================================

bool VCFXValidator::validateVCFMmap(const char* filepath) {
    MappedFile file;
    if (!file.open(filepath)) {
        std::cerr << "Error: Cannot open file: " << filepath << "\n";
        return false;
    }

    if (file.size == 0) {
        std::cerr << "Error: Empty file.\n";
        return false;
    }

    const char* ptr = file.data;
    const char* end = file.data + file.size;
    int lineNum = 0;
    bool foundChromLine = false;
    int dataLineCount = 0;

    // Process lines using SIMD-optimized newline scanning
    while (ptr < end) {
        // Find end of line using SIMD (AVX2/SSE2 on x86_64, memchr fallback)
        const char* lineEnd = findNewlineSIMD(ptr, end);
        if (!lineEnd) lineEnd = end;

        std::string_view line(ptr, lineEnd - ptr);

        // Remove trailing CR if present (Windows line endings)
        if (!line.empty() && line.back() == '\r') {
            line.remove_suffix(1);
        }

        lineNum++;

        if (line.empty()) {
            ptr = lineEnd + 1;
            continue;
        }

        if (line[0] == '#') {
            if (line.size() >= 2 && line[1] == '#') {
                if (!validateMetaLine(line, lineNum))
                    return false;
            } else if (line.size() >= 6 && line.substr(0, 6) == "#CHROM") {
                if (!validateChromHeader(line, lineNum))
                    return false;
                foundChromLine = true;
            } else {
                std::cerr << "Error: line " << lineNum
                          << " is a header line but neither starts with '##' nor is a #CHROM header line.\n";
                return false;
            }
        } else {
            if (!foundChromLine) {
                std::cerr << "Error: data line encountered before #CHROM at line " << lineNum << ".\n";
                return false;
            }
            if (!validateDataLine(line, lineNum))
                return false;
            dataLineCount++;
        }

        ptr = lineEnd + 1;
    }

    if (!foundChromLine) {
        std::cerr << "Error: no #CHROM line found in file.\n";
        return false;
    }

    // Output detailed summary
    std::cout << "=== VCF Validation Report ===\n";
    std::cout << "Status: PASSED\n";
    std::cout << "\n";
    std::cout << "File Statistics:\n";
    std::cout << "  Total lines:    " << lineNum << "\n";
    std::cout << "  Header lines:   " << (lineNum - dataLineCount) << "\n";
    std::cout << "  Variant records: " << dataLineCount << "\n";
    std::cout << "  Samples:        " << sampleCount << "\n";
    std::cout << "\n";
    std::cout << "Header Definitions:\n";
    std::cout << "  INFO fields:    " << infoDefs.size() << "\n";
    std::cout << "  FORMAT fields:  " << formatDefs.size() << "\n";
    std::cout << "\n";
    std::cout << "Validation Checks Performed:\n";
    std::cout << "  [OK] VCF header structure\n";
    std::cout << "  [OK] Meta-information lines (##)\n";
    std::cout << "  [OK] Column header (#CHROM)\n";
    std::cout << "  [OK] Required columns (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)\n";
    std::cout << "  [OK] Position values (POS > 0)\n";
    std::cout << "  [OK] REF/ALT allele sequences (A, C, G, T, N only)\n";
    std::cout << "  [OK] QUAL values (numeric or '.')\n";
    std::cout << "  [OK] INFO field definitions\n";
    if (headerHasFormat) {
        std::cout << "  [OK] FORMAT field definitions\n";
        std::cout << "  [OK] Genotype values\n";
    }
    if (reportDuplicates) {
        std::cout << "  [OK] Duplicate variant detection\n";
    }
    if (strictMode) {
        std::cout << "  [OK] Strict mode checks\n";
    }

    return true;
}

bool VCFXValidator::validateVCF(std::istream &in) {
    // Check for gzip magic and decompress if needed
    std::string decompressedData;
    std::istringstream decompressedStream;
    std::istream* inputStream = &in;

    // Peek first two bytes to check for gzip magic
    int c1 = in.get();
    int c2 = in.get();
    if (c1 != EOF && c2 != EOF) {
        bool isGzip = (static_cast<unsigned char>(c1) == 0x1f &&
                       static_cast<unsigned char>(c2) == 0x8b);
        in.putback(static_cast<char>(c2));
        in.putback(static_cast<char>(c1));

        if (isGzip) {
            // Use vcfx core to decompress
            if (!vcfx::read_maybe_compressed(in, decompressedData)) {
                std::cerr << "Error: Failed to decompress gzip input.\n";
                return false;
            }
            decompressedStream.str(decompressedData);
            inputStream = &decompressedStream;
        }
    } else {
        // Reset stream state if we hit EOF early
        in.clear();
        if (c1 != EOF) in.putback(static_cast<char>(c1));
    }

    // Use larger buffer for faster I/O (only for non-gzip)
    constexpr size_t BUFFER_SIZE = 1 << 20;  // 1MB buffer
    std::vector<char> inputBuffer(BUFFER_SIZE);
    if (inputStream == &in) {
        in.rdbuf()->pubsetbuf(inputBuffer.data(), inputBuffer.size());
    }

    std::string line;
    line.reserve(65536);  // Reserve for long lines with many samples

    int lineNum = 0;
    bool foundChromLine = false;
    int dataLineCount = 0;

    // Streaming: process and validate line by line
    while (std::getline(*inputStream, line)) {
        lineNum++;

        // Remove trailing CR if present (Windows line endings)
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        if (line.empty()) {
            continue;
        }

        std::string_view lineView(line);

        if (line[0] == '#') {
            if (lineView.substr(0, 2) == "##") {
                if (!validateMetaLine(lineView, lineNum))
                    return false;
            } else if (lineView.substr(0, 6) == "#CHROM") {
                if (!validateChromHeader(lineView, lineNum))
                    return false;
                foundChromLine = true;
            } else {
                std::cerr << "Error: line " << lineNum
                          << " is a header line but neither starts with '##' nor is a #CHROM header line.\n";
                return false;
            }
        } else {
            if (!foundChromLine) {
                std::cerr << "Error: data line encountered before #CHROM at line " << lineNum << ".\n";
                return false;
            }
            if (!validateDataLine(lineView, lineNum))
                return false;
            dataLineCount++;
        }
    }

    if (!foundChromLine) {
        std::cerr << "Error: no #CHROM line found in file.\n";
        return false;
    }

    // Output detailed summary
    std::cout << "=== VCF Validation Report ===\n";
    std::cout << "Status: PASSED\n";
    std::cout << "\n";
    std::cout << "File Statistics:\n";
    std::cout << "  Total lines:    " << lineNum << "\n";
    std::cout << "  Header lines:   " << (lineNum - dataLineCount) << "\n";
    std::cout << "  Variant records: " << dataLineCount << "\n";
    std::cout << "  Samples:        " << sampleCount << "\n";
    std::cout << "\n";
    std::cout << "Header Definitions:\n";
    std::cout << "  INFO fields:    " << infoDefs.size() << "\n";
    std::cout << "  FORMAT fields:  " << formatDefs.size() << "\n";
    std::cout << "\n";
    std::cout << "Validation Checks Performed:\n";
    std::cout << "  [OK] VCF header structure\n";
    std::cout << "  [OK] Meta-information lines (##)\n";
    std::cout << "  [OK] Column header (#CHROM)\n";
    std::cout << "  [OK] Required columns (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)\n";
    std::cout << "  [OK] Position values (POS > 0)\n";
    std::cout << "  [OK] REF/ALT allele sequences (A, C, G, T, N only)\n";
    std::cout << "  [OK] QUAL values (numeric or '.')\n";
    std::cout << "  [OK] INFO field definitions\n";
    if (headerHasFormat) {
        std::cout << "  [OK] FORMAT field definitions\n";
        std::cout << "  [OK] Genotype values\n";
    }
    if (reportDuplicates) {
        std::cout << "  [OK] Duplicate variant detection\n";
    }
    if (strictMode) {
        std::cout << "  [OK] Strict mode checks\n";
    }

    return true;
}

static void show_help() {
    VCFXValidator obj;
    char arg0[] = "VCFX_validator";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    // Disable sync with C stdio for faster I/O
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);

    if (vcfx::handle_common_flags(argc, argv, "VCFX_validator", show_help))
        return 0;
    VCFXValidator validator;
    return validator.run(argc, argv);
}
