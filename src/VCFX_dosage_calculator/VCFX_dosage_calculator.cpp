#include "VCFX_dosage_calculator.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <getopt.h>

// ---------------------------------------------------------------------------
// displayHelp: Prints usage information to standard output.
// ---------------------------------------------------------------------------
void VCFXDosageCalculator::displayHelp() {
    std::cout << "VCFX_dosage_calculator: Calculate genotype dosage for each variant in a VCF file.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_dosage_calculator [options] < input.vcf > dosage_output.txt\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help    Display this help message and exit\n\n";
    std::cout << "Description:\n";
    std::cout << "  For each variant in the input VCF, the tool computes the dosage for each sample\n";
    std::cout << "  based on the genotype (GT) field. Dosage is defined as the number of alternate\n";
    std::cout << "  alleles (i.e. each allele > 0 counts as 1). Thus:\n";
    std::cout << "    0/0  => dosage 0\n";
    std::cout << "    0/1  => dosage 1\n";
    std::cout << "    1/1  => dosage 2\n";
    std::cout << "    1/2  => dosage 2  (each alternate, regardless of numeric value, counts as 1)\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_dosage_calculator < input.vcf > dosage_output.txt\n";
}

// ---------------------------------------------------------------------------
// run: Parses command-line arguments and calls calculateDosage.
// ---------------------------------------------------------------------------
int VCFXDosageCalculator::run(int argc, char *argv[]) {
    int opt;
    bool showHelp = false;

    static struct option long_options[] = {{"help", no_argument, 0, 'h'}, {0, 0, 0, 0}};

    while ((opt = getopt_long(argc, argv, "h", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
            break;
        default:
            showHelp = true;
        }
    }
    if (showHelp) {
        displayHelp();
        return 0;
    }

    calculateDosage(std::cin, std::cout);
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
