#include "VCFX_dosage_calculator.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <getopt.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

// SIMD detection for ARM (Apple Silicon)
#if defined(__aarch64__) || defined(__ARM_NEON)
#include <arm_neon.h>
#define USE_NEON 1
#endif

// ---------------------------------------------------------------------------
// displayHelp: Prints usage information to standard output.
// ---------------------------------------------------------------------------
void VCFXDosageCalculator::displayHelp() {
    std::cout << "VCFX_dosage_calculator: Calculate genotype dosage for each variant in a VCF file.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_dosage_calculator [options] [input.vcf]\n";
    std::cout << "  VCFX_dosage_calculator [options] < input.vcf > dosage_output.txt\n\n";
    std::cout << "Options:\n";
    std::cout << "  -i, --input FILE  Input VCF file (uses mmap for best performance)\n";
    std::cout << "  -q, --quiet       Suppress warning messages\n";
    std::cout << "  -h, --help        Display this help message and exit\n\n";
    std::cout << "Description:\n";
    std::cout << "  For each variant in the input VCF, the tool computes the dosage for each sample\n";
    std::cout << "  based on the genotype (GT) field. Dosage is defined as the number of alternate\n";
    std::cout << "  alleles (i.e. each allele > 0 counts as 1). Thus:\n";
    std::cout << "    0/0  => dosage 0\n";
    std::cout << "    0/1  => dosage 1\n";
    std::cout << "    1/1  => dosage 2\n";
    std::cout << "    1/2  => dosage 2  (each alternate, regardless of numeric value, counts as 1)\n\n";
    std::cout << "Performance:\n";
    std::cout << "  When using -i/--input, the tool uses memory-mapped I/O for\n";
    std::cout << "  ~10-15x faster processing of large files.\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_dosage_calculator -i input.vcf > dosage_output.txt\n";
    std::cout << "  VCFX_dosage_calculator < input.vcf > dosage_output.txt\n";
}

// ---------------------------------------------------------------------------
// run: Parses command-line arguments and calls calculateDosage or mmap mode.
// ---------------------------------------------------------------------------
int VCFXDosageCalculator::run(int argc, char *argv[]) {
    int opt;
    bool showHelp = false;
    const char* inputFile = nullptr;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"input", required_argument, 0, 'i'},
        {"quiet", no_argument, 0, 'q'},
        {0, 0, 0, 0}
    };

    // Reset getopt
    optind = 1;

    while ((opt = getopt_long(argc, argv, "hi:q", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
            break;
        case 'i':
            inputFile = optarg;
            break;
        case 'q':
            quietMode = true;
            break;
        default:
            showHelp = true;
        }
    }

    // Check for positional file argument
    if (!inputFile && optind < argc) {
        inputFile = argv[optind];
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    // Use mmap mode if file provided, otherwise stdin
    if (inputFile) {
        if (!processFileMmap(inputFile, std::cout)) {
            return 1;
        }
    } else {
        calculateDosage(std::cin, std::cout);
    }
    return 0;
}

// ===========================================================================
// FULLY OPTIMIZED: Zero-allocation line parsing and dosage calculation
// ===========================================================================

// Zero-allocation dosage parsing from raw pointer
// Returns dosage directly without creating any string objects
// -1 = invalid/missing
static inline int parseDosageInline(const char* gt, size_t gtLen) {
    if (gtLen == 0) return -1;

    int dosage = 0;
    size_t pos = 0;
    int alleleCount = 0;

    while (pos < gtLen) {
        // Skip separators
        while (pos < gtLen && (gt[pos] == '/' || gt[pos] == '|')) {
            pos++;
        }
        if (pos >= gtLen) break;

        // Check for missing allele
        if (gt[pos] == '.') {
            return -1;  // Missing data
        }

        // Parse numeric allele
        int allele = 0;
        bool hasDigit = false;
        while (pos < gtLen && gt[pos] >= '0' && gt[pos] <= '9') {
            allele = allele * 10 + (gt[pos] - '0');
            hasDigit = true;
            pos++;
        }

        if (!hasDigit) {
            return -1;  // Invalid non-numeric character
        }

        if (allele > 0) {
            dosage++;
        }
        alleleCount++;

        // More than 2 alleles = invalid (non-diploid)
        if (alleleCount > 2) {
            return -1;
        }
    }

    // Require exactly 2 alleles for diploid
    return (alleleCount == 2) ? dosage : -1;
}

// Zero-allocation helper: find GT index in FORMAT string using raw pointers
// Returns -1 if GT not found
static inline int findGTIndexRaw(const char* format, size_t formatLen) {
    int index = 0;
    size_t pos = 0;
    size_t start = 0;

    while (pos <= formatLen) {
        if (pos == formatLen || format[pos] == ':') {
            size_t tokenLen = pos - start;
            // Check if this token is "GT" (exactly 2 chars)
            if (tokenLen == 2 && format[start] == 'G' && format[start + 1] == 'T') {
                return index;
            }
            index++;
            start = pos + 1;
        }
        pos++;
    }
    return -1;
}

// Zero-allocation helper: extract GT field from sample data
// Returns pointer to GT start and its length
static inline bool extractGTFromSample(const char* sample, size_t sampleLen,
                                        int gtIndex, const char*& gtStart, size_t& gtLen) {
    if (gtIndex < 0) return false;

    int currentIndex = 0;
    size_t pos = 0;
    size_t fieldStart = 0;

    while (pos <= sampleLen) {
        if (pos == sampleLen || sample[pos] == ':') {
            if (currentIndex == gtIndex) {
                gtStart = sample + fieldStart;
                gtLen = pos - fieldStart;
                return gtLen > 0;
            }
            currentIndex++;
            fieldStart = pos + 1;
        }
        pos++;
    }
    return false;
}

// ---------------------------------------------------------------------------
// calculateDosage: ULTRA-OPTIMIZED - raw char buffer, no std::string in hot loop
// Uses fixed char buffer for dosages to eliminate all string operations
// ---------------------------------------------------------------------------
void VCFXDosageCalculator::calculateDosage(std::istream &in, std::ostream &out) {
    std::string line;
    bool headerParsed = false;

    // 1MB output buffer for reduced syscall overhead
    static constexpr size_t OUTPUT_BUFFER_SIZE = 1024 * 1024;
    char* outputBuffer = new char[OUTPUT_BUFFER_SIZE + 32768];  // Extra space for last line
    size_t bufPos = 0;

    // Pre-sized dosage buffer - 2 chars per sample (digit + comma) + some slack
    // For 2504 samples: 2504 * 2 = 5008 + header ~ 6000 bytes max
    static constexpr size_t DOSAGE_BUFFER_SIZE = 8192;
    char dosageBuffer[DOSAGE_BUFFER_SIZE];

    // Write header
    const char* header = "CHROM\tPOS\tID\tREF\tALT\tDosages\n";
    size_t headerLen = 29;  // strlen("CHROM\tPOS\tID\tREF\tALT\tDosages\n")
    std::memcpy(outputBuffer, header, headerLen);
    bufPos = headerLen;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            if (line.rfind("#CHROM", 0) == 0) {
                headerParsed = true;
            }
            continue;
        }

        if (!headerParsed) {
            std::cerr << "Error: VCF header (#CHROM) not found before variant records.\n";
            delete[] outputBuffer;
            return;
        }

        const char* linePtr = line.c_str();
        size_t lineLen = line.size();

        // Parse fields using zero-allocation pointer arithmetic
        const char* fields[10];
        size_t fieldLens[10];
        int fieldCount = 0;
        size_t fieldStart = 0;

        for (size_t i = 0; i <= lineLen && fieldCount < 10; i++) {
            if (i == lineLen || linePtr[i] == '\t') {
                fields[fieldCount] = linePtr + fieldStart;
                fieldLens[fieldCount] = i - fieldStart;
                fieldCount++;
                fieldStart = i + 1;
            }
        }

        if (fieldCount < 10) {
            std::cerr << "Warning: Skipping VCF line with fewer than 10 fields.\n";
            continue;
        }

        // Build header part directly into output buffer
        // CHROM
        std::memcpy(outputBuffer + bufPos, fields[0], fieldLens[0]);
        bufPos += fieldLens[0];
        outputBuffer[bufPos++] = '\t';
        // POS
        std::memcpy(outputBuffer + bufPos, fields[1], fieldLens[1]);
        bufPos += fieldLens[1];
        outputBuffer[bufPos++] = '\t';
        // ID
        std::memcpy(outputBuffer + bufPos, fields[2], fieldLens[2]);
        bufPos += fieldLens[2];
        outputBuffer[bufPos++] = '\t';
        // REF
        std::memcpy(outputBuffer + bufPos, fields[3], fieldLens[3]);
        bufPos += fieldLens[3];
        outputBuffer[bufPos++] = '\t';
        // ALT
        std::memcpy(outputBuffer + bufPos, fields[4], fieldLens[4]);
        bufPos += fieldLens[4];
        outputBuffer[bufPos++] = '\t';

        // Find GT index in FORMAT field (field 8)
        int gtIndex = findGTIndexRaw(fields[8], fieldLens[8]);

        if (gtIndex == -1) {
            outputBuffer[bufPos++] = 'N';
            outputBuffer[bufPos++] = 'A';
            outputBuffer[bufPos++] = '\n';
            if (bufPos >= OUTPUT_BUFFER_SIZE) {
                out.write(outputBuffer, bufPos);
                bufPos = 0;
            }
            continue;
        }

        // Process all samples into dosage buffer first
        char* dosagePtr = dosageBuffer;
        const char* samplePtr = fields[9];
        const char* lineEnd = linePtr + lineLen;
        bool firstSample = true;

        while (samplePtr < lineEnd) {
            // Find end of current sample
            const char* sampleEnd = samplePtr;
            while (sampleEnd < lineEnd && *sampleEnd != '\t') {
                sampleEnd++;
            }
            size_t sampleLen = sampleEnd - samplePtr;

            if (!firstSample) {
                *dosagePtr++ = ',';
            }
            firstSample = false;

            // Extract and parse GT
            const char* gtStart;
            size_t gtLen;
            if (extractGTFromSample(samplePtr, sampleLen, gtIndex, gtStart, gtLen)) {
                int dosage = parseDosageInline(gtStart, gtLen);
                if (dosage < 0) {
                    *dosagePtr++ = 'N';
                    *dosagePtr++ = 'A';
                } else {
                    *dosagePtr++ = static_cast<char>('0' + dosage);
                }
            } else {
                *dosagePtr++ = 'N';
                *dosagePtr++ = 'A';
            }

            samplePtr = (sampleEnd < lineEnd) ? sampleEnd + 1 : lineEnd;
        }
        *dosagePtr++ = '\n';

        // Copy dosages to output buffer
        size_t dosageLen = dosagePtr - dosageBuffer;
        std::memcpy(outputBuffer + bufPos, dosageBuffer, dosageLen);
        bufPos += dosageLen;

        // Flush when buffer is full
        if (bufPos >= OUTPUT_BUFFER_SIZE) {
            out.write(outputBuffer, bufPos);
            bufPos = 0;
        }
    }

    // Flush remaining buffer
    if (bufPos > 0) {
        out.write(outputBuffer, bufPos);
    }
    delete[] outputBuffer;
}

// ===========================================================================
// MEMORY-MAPPED I/O MODE: Fast processing for file inputs
// ===========================================================================

// SIMD-optimized newline finding using memchr (libc-optimized)
static inline const char* findNewline(const char* p, const char* end) {
    return static_cast<const char*>(memchr(p, '\n', static_cast<size_t>(end - p)));
}

// ---------------------------------------------------------------------------
// processFileMmap: Memory-mapped file processing for maximum performance
// Uses the same parsing logic as calculateDosage but with mmap input
// ---------------------------------------------------------------------------
bool VCFXDosageCalculator::processFileMmap(const char* filename, std::ostream &out) {
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

    size_t fileSize = static_cast<size_t>(st.st_size);
    if (fileSize == 0) {
        close(fd);
        return true;  // Empty file is valid
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
    bool headerParsed = false;

    // 1MB output buffer for reduced syscall overhead
    static constexpr size_t OUTPUT_BUFFER_SIZE = 1024 * 1024;
    char* outputBuffer = new char[OUTPUT_BUFFER_SIZE + 32768];  // Extra space for last line
    size_t bufPos = 0;

    // Pre-sized dosage buffer - 2 chars per sample (digit + comma) + some slack
    static constexpr size_t DOSAGE_BUFFER_SIZE = 8192;
    char dosageBuffer[DOSAGE_BUFFER_SIZE];

    // Write header
    const char* header = "CHROM\tPOS\tID\tREF\tALT\tDosages\n";
    size_t headerLen = 29;  // strlen("CHROM\tPOS\tID\tREF\tALT\tDosages\n")
    std::memcpy(outputBuffer, header, headerLen);
    bufPos = headerLen;

    const char* ptr = data;

    while (ptr < end) {
        // Find end of line
        const char* lineEnd = findNewline(ptr, end);
        if (!lineEnd) {
            lineEnd = end;  // Last line without trailing newline
        }

        size_t lineLen = static_cast<size_t>(lineEnd - ptr);

        // Handle Windows line endings
        if (lineLen > 0 && ptr[lineLen - 1] == '\r') {
            lineLen--;
        }

        // Skip empty lines
        if (lineLen == 0) {
            ptr = lineEnd + 1;
            continue;
        }

        // Handle header lines
        if (*ptr == '#') {
            if (lineLen >= 6 && ptr[0] == '#' && ptr[1] == 'C' && ptr[2] == 'H' &&
                ptr[3] == 'R' && ptr[4] == 'O' && ptr[5] == 'M') {
                headerParsed = true;
            }
            ptr = lineEnd + 1;
            continue;
        }

        if (!headerParsed) {
            std::cerr << "Error: VCF header (#CHROM) not found before variant records.\n";
            delete[] outputBuffer;
            munmap(mapped, fileSize);
            close(fd);
            return false;
        }

        // Parse fields using zero-allocation pointer arithmetic
        // VCF format: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLES...
        const char* linePtr = ptr;
        const char* fields[10];
        size_t fieldLens[10];
        int fieldCount = 0;
        size_t fieldStart = 0;

        for (size_t i = 0; i <= lineLen && fieldCount < 10; i++) {
            if (i == lineLen || linePtr[i] == '\t') {
                fields[fieldCount] = linePtr + fieldStart;
                fieldLens[fieldCount] = i - fieldStart;
                fieldCount++;
                fieldStart = i + 1;
            }
        }

        if (fieldCount < 10) {
            if (!quietMode) {
                std::cerr << "Warning: Skipping VCF line with fewer than 10 fields.\n";
            }
            ptr = lineEnd + 1;
            continue;
        }

        // Build output directly into buffer
        // CHROM
        std::memcpy(outputBuffer + bufPos, fields[0], fieldLens[0]);
        bufPos += fieldLens[0];
        outputBuffer[bufPos++] = '\t';
        // POS
        std::memcpy(outputBuffer + bufPos, fields[1], fieldLens[1]);
        bufPos += fieldLens[1];
        outputBuffer[bufPos++] = '\t';
        // ID
        std::memcpy(outputBuffer + bufPos, fields[2], fieldLens[2]);
        bufPos += fieldLens[2];
        outputBuffer[bufPos++] = '\t';
        // REF
        std::memcpy(outputBuffer + bufPos, fields[3], fieldLens[3]);
        bufPos += fieldLens[3];
        outputBuffer[bufPos++] = '\t';
        // ALT
        std::memcpy(outputBuffer + bufPos, fields[4], fieldLens[4]);
        bufPos += fieldLens[4];
        outputBuffer[bufPos++] = '\t';

        // Find GT index in FORMAT field (field 8)
        int gtIndex = findGTIndexRaw(fields[8], fieldLens[8]);

        if (gtIndex == -1) {
            outputBuffer[bufPos++] = 'N';
            outputBuffer[bufPos++] = 'A';
            outputBuffer[bufPos++] = '\n';
            if (bufPos >= OUTPUT_BUFFER_SIZE) {
                out.write(outputBuffer, bufPos);
                bufPos = 0;
            }
            ptr = lineEnd + 1;
            continue;
        }

        // Process all samples into dosage buffer
        char* dosagePtr = dosageBuffer;
        const char* samplePtr = fields[9];
        const char* sampleEnd = linePtr + lineLen;
        bool firstSample = true;

        while (samplePtr < sampleEnd) {
            // Find end of current sample
            const char* nextTab = samplePtr;
            while (nextTab < sampleEnd && *nextTab != '\t') {
                nextTab++;
            }
            size_t sampleLen = static_cast<size_t>(nextTab - samplePtr);

            if (!firstSample) {
                *dosagePtr++ = ',';
            }
            firstSample = false;

            // Extract and parse GT
            const char* gtStart;
            size_t gtLen;
            if (extractGTFromSample(samplePtr, sampleLen, gtIndex, gtStart, gtLen)) {
                int dosage = parseDosageInline(gtStart, gtLen);
                if (dosage < 0) {
                    *dosagePtr++ = 'N';
                    *dosagePtr++ = 'A';
                } else {
                    *dosagePtr++ = static_cast<char>('0' + dosage);
                }
            } else {
                *dosagePtr++ = 'N';
                *dosagePtr++ = 'A';
            }

            samplePtr = (nextTab < sampleEnd) ? nextTab + 1 : sampleEnd;
        }
        *dosagePtr++ = '\n';

        // Copy dosages to output buffer
        size_t dosageLen = static_cast<size_t>(dosagePtr - dosageBuffer);
        std::memcpy(outputBuffer + bufPos, dosageBuffer, dosageLen);
        bufPos += dosageLen;

        // Flush when buffer is full
        if (bufPos >= OUTPUT_BUFFER_SIZE) {
            out.write(outputBuffer, bufPos);
            bufPos = 0;
        }

        ptr = lineEnd + 1;
    }

    // Flush remaining buffer
    if (bufPos > 0) {
        out.write(outputBuffer, bufPos);
    }

    delete[] outputBuffer;
    munmap(mapped, fileSize);
    close(fd);
    return true;
}

// ---------------------------------------------------------------------------
// split: Helper function (kept for compatibility)
// ---------------------------------------------------------------------------
std::vector<std::string> VCFXDosageCalculator::split(const std::string &str, char delimiter) {
    std::vector<std::string> tokens;
    tokens.reserve(8);
    size_t start = 0;
    size_t end;
    while ((end = str.find(delimiter, start)) != std::string::npos) {
        tokens.emplace_back(str, start, end - start);
        start = end + 1;
    }
    tokens.emplace_back(str, start);
    return tokens;
}

static void show_help() {
    VCFXDosageCalculator obj;
    char arg0[] = "VCFX_dosage_calculator";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();
    if (vcfx::handle_common_flags(argc, argv, "VCFX_dosage_calculator", show_help))
        return 0;
    VCFXDosageCalculator dosageCalculator;
    return dosageCalculator.run(argc, argv);
}
