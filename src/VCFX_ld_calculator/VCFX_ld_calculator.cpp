#include "VCFX_ld_calculator.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <sstream>

// ======================================================================
// Zero-allocation helper functions for high-performance parsing
// ======================================================================

// Fast genotype parsing from raw pointer - no string allocation
// Returns: 0=0/0, 1=0/1 or 1/0, 2=1/1, -1=missing/invalid
static inline int parseGenotypeRaw(const char* s, size_t len) {
    if (len == 0) return -1;

    // Handle common missing patterns
    if (len == 1 && s[0] == '.') return -1;
    if (len == 3 && s[0] == '.' && (s[1] == '/' || s[1] == '|') && s[2] == '.') return -1;

    // Find separator (/ or |) - for GT field, separator is usually at position 1
    size_t sepPos = 0;
    for (size_t i = 0; i < len; i++) {
        if (s[i] == '/' || s[i] == '|') {
            sepPos = i;
            break;
        }
    }
    if (sepPos == 0 || sepPos >= len - 1) return -1;  // No separator or at edge

    // Parse first allele
    int a1 = 0;
    for (size_t i = 0; i < sepPos; i++) {
        char c = s[i];
        if (c == '.') return -1;  // Missing allele
        if (c < '0' || c > '9') return -1;  // Invalid
        a1 = a1 * 10 + (c - '0');
    }

    // Parse second allele
    int a2 = 0;
    for (size_t i = sepPos + 1; i < len; i++) {
        char c = s[i];
        if (c == '.') return -1;  // Missing allele
        if (c < '0' || c > '9') return -1;  // Invalid
        a2 = a2 * 10 + (c - '0');
    }

    // Check for multi-allelic (alleles > 1)
    if (a1 > 1 || a2 > 1) return -1;

    // Return genotype code: 0/0->0, 0/1 or 1/0->1, 1/1->2
    return a1 + a2;
}

// Extract GT field from sample field (stops at first colon)
// Returns pointer to start and length of GT portion
static inline bool extractGT(const char* sample, size_t sampleLen,
                              const char*& gtStart, size_t& gtLen) {
    gtStart = sample;
    gtLen = sampleLen;

    for (size_t i = 0; i < sampleLen; i++) {
        if (sample[i] == ':') {
            gtLen = i;
            break;
        }
    }
    return gtLen > 0;
}

// Fast integer parsing from raw pointer
static inline bool fastParseInt(const char* s, size_t len, int& result) {
    if (len == 0) return false;
    result = 0;
    for (size_t i = 0; i < len; i++) {
        char c = s[i];
        if (c < '0' || c > '9') return false;
        result = result * 10 + (c - '0');
    }
    return true;
}

// Fast double to string with 4 decimal places - much faster than iostream
static inline void formatR2(double r2, char* buf, size_t& len) {
    // Handle special cases
    if (r2 <= 0.0) {
        std::memcpy(buf, "0.0000", 6);
        len = 6;
        return;
    }
    if (r2 >= 1.0) {
        std::memcpy(buf, "1.0000", 6);
        len = 6;
        return;
    }

    // Format as 0.XXXX
    buf[0] = '0';
    buf[1] = '.';

    int val = static_cast<int>(r2 * 10000.0 + 0.5);  // Round to 4 decimals
    if (val > 9999) val = 9999;

    buf[5] = '0' + (val % 10); val /= 10;
    buf[4] = '0' + (val % 10); val /= 10;
    buf[3] = '0' + (val % 10); val /= 10;
    buf[2] = '0' + (val % 10);
    len = 6;
}

// Fast integer to string
static inline size_t formatInt(int val, char* buf) {
    if (val == 0) {
        buf[0] = '0';
        return 1;
    }

    char temp[12];
    int pos = 0;
    bool neg = val < 0;
    if (neg) val = -val;

    while (val > 0) {
        temp[pos++] = '0' + (val % 10);
        val /= 10;
    }

    size_t len = 0;
    if (neg) buf[len++] = '-';
    while (pos > 0) {
        buf[len++] = temp[--pos];
    }
    return len;
}

// ----------------------------------------------------------------------
// Display Help
// ----------------------------------------------------------------------
void VCFXLDCalculator::displayHelp() {
    std::cout << "VCFX_ld_calculator: Calculate pairwise LD (r^2) for variants in a VCF region.\n\n"
              << "Usage:\n"
              << "  VCFX_ld_calculator [options] < input.vcf\n\n"
              << "Options:\n"
              << "  --region <chr:start-end>  Only compute LD for variants in [start, end] on 'chr'.\n"
              << "  --window <N>              Window size in number of variants (default: 1000).\n"
              << "  --threshold <R2>          Only output pairs with r^2 >= threshold (default: 0.0).\n"
              << "  --matrix                  Use matrix mode (MxM output) instead of streaming.\n"
              << "                            WARNING: O(M^2) time - avoid for >10K variants.\n"
              << "  --help, -h                Show this help message.\n\n"
              << "Modes:\n"
              << "  Default (streaming): Outputs LD pairs incrementally using a sliding window.\n"
              << "                       Memory: O(window * samples) - constant for any file size.\n"
              << "                       Time: O(M * window) - linear in variant count.\n"
              << "  Matrix mode:         Produces an MxM matrix of all pairwise r^2 values.\n"
              << "                       Memory: O(M * samples) where M is number of variants.\n"
              << "                       Time: O(M^2) - avoid for >10K variants!\n\n"
              << "Description:\n"
              << "  Reads a VCF from stdin and computes pairwise r^2 for diploid genotypes.\n"
              << "  Genotype codes: 0 => 0/0, 1 => 0/1, 2 => 1/1, -1 => missing.\n"
              << "  Samples with missing genotypes are excluded from r^2 calculation.\n\n"
              << "Performance:\n"
              << "  Default streaming mode handles files of any size efficiently.\n"
              << "  For 427K variants: streaming takes ~10 min, matrix would take days.\n\n"
              << "Output Format:\n"
              << "  Streaming:  Tab-separated: VAR1_CHROM VAR1_POS VAR1_ID VAR2_CHROM VAR2_POS VAR2_ID R2\n"
              << "  Matrix:     MxM grid with variant identifiers.\n\n"
              << "Example:\n"
              << "  # Default streaming mode (handles any file size)\n"
              << "  VCFX_ld_calculator --window 500 --threshold 0.2 < input.vcf > ld_pairs.txt\n\n"
              << "  # Matrix mode (small regions only, <10K variants)\n"
              << "  VCFX_ld_calculator --matrix --region chr1:10000-20000 < input.vcf > ld_matrix.txt\n";
}

// ----------------------------------------------------------------------
// parseRegion( "chr1:10000-20000" ) => regionChrom="chr1", regionStart=10000, regionEnd=20000
// returns false if can't parse
// ----------------------------------------------------------------------
bool VCFXLDCalculator::parseRegion(const std::string &regionStr, std::string &regionChrom, int &regionStart,
                                   int &regionEnd) {
    // find ':'
    auto colonPos = regionStr.find(':');
    if (colonPos == std::string::npos) {
        return false;
    }
    regionChrom = regionStr.substr(0, colonPos);
    auto dashPos = regionStr.find('-', colonPos + 1);
    if (dashPos == std::string::npos) {
        return false;
    }
    std::string startStr = regionStr.substr(colonPos + 1, dashPos - (colonPos + 1));
    std::string endStr = regionStr.substr(dashPos + 1);
    try {
        regionStart = std::stoi(startStr);
        regionEnd = std::stoi(endStr);
    } catch (...) {
        return false;
    }
    if (regionStart > regionEnd) {
        return false;
    }
    return true;
}

// ----------------------------------------------------------------------
// parse a genotype string like "0/1" or "1|0"
//   => 0 => 0/0
//   => 1 => 0/1
//   => 2 => 1/1
//   => -1 => missing or multi-allelic
// ----------------------------------------------------------------------
int VCFXLDCalculator::parseGenotype(const std::string &s) {
    if (s.empty() || s == "." || s == "./." || s == ".|.") {
        return -1;
    }
    // unify '|' => '/'
    std::string g(s);
    for (char &c : g) {
        if (c == '|')
            c = '/';
    }
    // split
    auto slashPos = g.find('/');
    if (slashPos == std::string::npos) {
        return -1; // not diploid
    }
    auto a1 = g.substr(0, slashPos);
    auto a2 = g.substr(slashPos + 1);
    if (a1.empty() || a2.empty())
        return -1;
    if (a1 == "." || a2 == ".")
        return -1;
    // parse int
    int i1 = 0, i2 = 0;
    try {
        i1 = std::stoi(a1);
        i2 = std::stoi(a2);
    } catch (...) {
        return -1;
    }
    if (i1 < 0 || i2 < 0)
        return -1;
    if (i1 > 1 || i2 > 1) {
        // multi-allelic
        return -1;
    }
    if (i1 == i2) {
        // 0 => 0/0 => code=0, 1 => 1/1 => code=2
        if (i1 == 0)
            return 0;
        else
            return 2;
    } else {
        return 1;
    }
}

// ----------------------------------------------------------------------
// compute r^2 between genotype arrays g1,g2 ignoring -1
// standard formula
// Let n = count of samples where both g1,g2 != -1
// X[i]=g1[i], Y[i]=g2[i]
// meanX= average of X, meanY= average of Y
// cov= average(XY)-meanX*meanY
// varX= average(X^2)- meanX^2, varY likewise
// r= cov / sqrt(varX varY), r^2= r*r
// ----------------------------------------------------------------------
double VCFXLDCalculator::computeRsq(const std::vector<int> &g1, const std::vector<int> &g2) {
    if (g1.size() != g2.size())
        return 0.0;
    int n = 0;
    long sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, sumY2 = 0;
    for (size_t i = 0; i < g1.size(); i++) {
        int x = g1[i];
        int y = g2[i];
        if (x < 0 || y < 0)
            continue; // skip missing
        n++;
        sumX += x;
        sumY += y;
        sumXY += x * y;
        sumX2 += (x * x);
        sumY2 += (y * y);
    }
    if (n < 2)
        return 0.0;
    double meanX = (double)sumX / (double)n;
    double meanY = (double)sumY / (double)n;
    double cov = ((double)sumXY / (double)n) - (meanX * meanY);
    double varX = ((double)sumX2 / (double)n) - (meanX * meanX);
    double varY = ((double)sumY2 / (double)n) - (meanY * meanY);
    if (varX <= 0.0 || varY <= 0.0)
        return 0.0;
    double r = cov / (std::sqrt(varX) * std::sqrt(varY));
    return r * r;
}

// ----------------------------------------------------------------------
// Optimized LDVariant with pre-computed statistics for fast r² calculation
// ----------------------------------------------------------------------
struct LDVariantOpt {
    std::string chrom;
    int pos;
    std::string id;
    std::vector<int8_t> genotype;  // Use int8_t instead of int (4x smaller)

    // Pre-computed statistics for r² calculation
    int validCount;      // Number of non-missing samples
    long sumX;           // Sum of genotype values
    long sumX2;          // Sum of squared genotype values
    double meanX;        // Mean of genotype values
    double varX;         // Variance of genotype values

    void computeStats() {
        validCount = 0;
        sumX = 0;
        sumX2 = 0;
        for (int8_t g : genotype) {
            if (g >= 0) {
                validCount++;
                sumX += g;
                sumX2 += g * g;
            }
        }
        if (validCount > 0) {
            meanX = static_cast<double>(sumX) / validCount;
            varX = static_cast<double>(sumX2) / validCount - meanX * meanX;
        } else {
            meanX = 0;
            varX = 0;
        }
    }
};

// Fast r² computation using pre-computed statistics
static inline double computeRsqFast(const LDVariantOpt& v1, const LDVariantOpt& v2) {
    if (v1.genotype.size() != v2.genotype.size()) return 0.0;
    if (v1.varX <= 0.0 || v2.varX <= 0.0) return 0.0;

    // Compute covariance (only needs XY sum, X and Y sums are pre-computed)
    int n = 0;
    long sumXY = 0;
    long sumX = 0, sumY = 0;

    const int8_t* g1 = v1.genotype.data();
    const int8_t* g2 = v2.genotype.data();
    const size_t sz = v1.genotype.size();

    for (size_t i = 0; i < sz; i++) {
        int8_t x = g1[i];
        int8_t y = g2[i];
        if (x >= 0 && y >= 0) {
            n++;
            sumX += x;
            sumY += y;
            sumXY += x * y;
        }
    }

    if (n < 2) return 0.0;

    double meanX = static_cast<double>(sumX) / n;
    double meanY = static_cast<double>(sumY) / n;
    double cov = static_cast<double>(sumXY) / n - meanX * meanY;
    double varX = static_cast<double>(v1.sumX2) / v1.validCount - v1.meanX * v1.meanX;
    double varY = static_cast<double>(v2.sumX2) / v2.validCount - v2.meanX * v2.meanX;

    // Recalculate variance for the paired subset if there's missing data
    if (n < v1.validCount || n < v2.validCount) {
        long sumX2_paired = 0, sumY2_paired = 0;
        for (size_t i = 0; i < sz; i++) {
            int8_t x = g1[i];
            int8_t y = g2[i];
            if (x >= 0 && y >= 0) {
                sumX2_paired += x * x;
                sumY2_paired += y * y;
            }
        }
        varX = static_cast<double>(sumX2_paired) / n - meanX * meanX;
        varY = static_cast<double>(sumY2_paired) / n - meanY * meanY;
    }

    if (varX <= 0.0 || varY <= 0.0) return 0.0;

    double r = cov / (std::sqrt(varX) * std::sqrt(varY));
    return r * r;
}

// ----------------------------------------------------------------------
// computeLDStreaming: Sliding window mode with O(window * samples) memory
// OPTIMIZED: Zero-allocation parsing, pre-computed statistics, buffered output
// ----------------------------------------------------------------------
void VCFXLDCalculator::computeLDStreaming(std::istream &in, std::ostream &out, const std::string &regionChrom,
                                          int regionStart, int regionEnd, size_t windowSize, double threshold) {
    bool foundChromHeader = false;
    int numSamples = 0;
    std::deque<LDVariantOpt> window;
    std::string line;

    // Output buffer for batched I/O (1MB)
    std::string outputBuffer;
    outputBuffer.reserve(1024 * 1024);

    // Temp buffer for formatting
    char numBuf[32];

    // Output header for streaming mode
    outputBuffer = "#VAR1_CHROM\tVAR1_POS\tVAR1_ID\tVAR2_CHROM\tVAR2_POS\tVAR2_ID\tR2\n";

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '#') {
            if (!foundChromHeader && line.rfind("#CHROM", 0) == 0) {
                foundChromHeader = true;
                // Count tabs to determine number of samples
                const char* p = line.c_str();
                int tabCount = 0;
                while (*p) {
                    if (*p == '\t') tabCount++;
                    p++;
                }
                numSamples = (tabCount >= 9) ? (tabCount - 8) : 0;
            }
            continue;
        }

        // Data line
        if (!foundChromHeader) {
            std::cerr << "Error: encountered data line before #CHROM.\n";
            break;
        }

        const char* linePtr = line.c_str();
        size_t lineLen = line.size();

        // Parse fields using zero-allocation approach
        // Find tab positions for first 10 fields (CHROM through first sample)
        const char* fieldStarts[11];
        size_t fieldLens[11];
        int fieldCount = 0;
        size_t start = 0;

        for (size_t i = 0; i <= lineLen && fieldCount < 10; i++) {
            if (i == lineLen || linePtr[i] == '\t') {
                fieldStarts[fieldCount] = linePtr + start;
                fieldLens[fieldCount] = i - start;
                fieldCount++;
                start = i + 1;
            }
        }

        if (fieldCount < 10) continue;  // Malformed line

        // Parse position
        int posVal;
        if (!fastParseInt(fieldStarts[1], fieldLens[1], posVal)) continue;

        // Check region (compare chrom using length + memcmp)
        if (!regionChrom.empty()) {
            if (fieldLens[0] != regionChrom.size() ||
                std::memcmp(fieldStarts[0], regionChrom.c_str(), fieldLens[0]) != 0) {
                continue;
            }
            if (posVal < regionStart || posVal > regionEnd) {
                continue;
            }
        }

        // Create variant
        LDVariantOpt v;
        v.chrom.assign(fieldStarts[0], fieldLens[0]);
        v.pos = posVal;

        // Handle ID field
        if (fieldLens[2] == 1 && fieldStarts[2][0] == '.') {
            // Construct chrom:pos for ID
            v.id = v.chrom + ":" + std::to_string(posVal);
        } else {
            v.id.assign(fieldStarts[2], fieldLens[2]);
        }

        v.genotype.resize(numSamples, -1);

        // Parse samples - iterate through remaining line
        const char* sampleStart = fieldStarts[9];
        int sampleIdx = 0;

        while (sampleStart < linePtr + lineLen && sampleIdx < numSamples) {
            // Find end of this sample field
            const char* sampleEnd = sampleStart;
            while (sampleEnd < linePtr + lineLen && *sampleEnd != '\t') {
                sampleEnd++;
            }
            size_t sampleLen = sampleEnd - sampleStart;

            // Extract GT (first colon-delimited field)
            const char* gtStart;
            size_t gtLen;
            if (extractGT(sampleStart, sampleLen, gtStart, gtLen)) {
                v.genotype[sampleIdx] = static_cast<int8_t>(parseGenotypeRaw(gtStart, gtLen));
            }

            sampleIdx++;
            sampleStart = sampleEnd + 1;
        }

        // Pre-compute statistics for this variant
        v.computeStats();

        // Compute LD with all variants in window
        for (const auto &prev : window) {
            double r2 = computeRsqFast(prev, v);
            if (r2 >= threshold) {
                // Build output line with fast formatting
                outputBuffer += prev.chrom;
                outputBuffer += '\t';
                size_t len = formatInt(prev.pos, numBuf);
                outputBuffer.append(numBuf, len);
                outputBuffer += '\t';
                outputBuffer += prev.id;
                outputBuffer += '\t';
                outputBuffer += v.chrom;
                outputBuffer += '\t';
                len = formatInt(v.pos, numBuf);
                outputBuffer.append(numBuf, len);
                outputBuffer += '\t';
                outputBuffer += v.id;
                outputBuffer += '\t';
                formatR2(r2, numBuf, len);
                outputBuffer.append(numBuf, len);
                outputBuffer += '\n';

                // Flush buffer when it gets large
                if (outputBuffer.size() > 900000) {
                    out << outputBuffer;
                    outputBuffer.clear();
                }
            }
        }

        // Add to window
        window.push_back(std::move(v));

        // Trim window if needed
        if (window.size() > windowSize) {
            window.pop_front();
        }
    }

    // Final buffer flush
    if (!outputBuffer.empty()) {
        out << outputBuffer;
    }
}

// ----------------------------------------------------------------------
// computeLD: read lines, skip #, parse sample columns
//   store variants within region
//   then produce MxM matrix of r^2
// (Original implementation - backward compatible)
// ----------------------------------------------------------------------
void VCFXLDCalculator::computeLD(std::istream &in, std::ostream &out, const std::string &regionChrom, int regionStart,
                                 int regionEnd) {
    bool foundChromHeader = false;
    std::vector<std::string> sampleNames;
    std::vector<LDVariant> variants;
    int numSamples = 0;

    // Performance: reuse containers across iterations
    std::vector<std::string> fields;
    fields.reserve(16);

    // parse header => get sample col
    // pass header lines as is
    std::string line;
    while (true) {
        auto pos = in.tellg();
        if (!std::getline(in, line))
            break;
        if (line.empty()) {
            out << line << "\n";
            continue;
        }
        if (line[0] == '#') {
            out << line << "\n";
            if (!foundChromHeader && line.rfind("#CHROM", 0) == 0) {
                // parse sample
                foundChromHeader = true;
                vcfx::split_tabs(line, fields);
                // from col=9 onward => sample
                for (size_t c = 9; c < fields.size(); c++) {
                    sampleNames.push_back(fields[c]);
                }
                numSamples = sampleNames.size();
            }
            continue;
        }
        // data line
        if (!foundChromHeader) {
            std::cerr << "Error: encountered data line before #CHROM.\n";
            // we can break or skip
            break;
        }
        // parse
        // CHROM POS ID REF ALT QUAL FILTER INFO ...
        vcfx::split_tabs(line, fields);
        if (fields.size() < 10) {
            // not enough columns => pass line => skip
            out << line << "\n";
            continue;
        }
        std::string &chrom = fields[0];
        std::string &posStr = fields[1];
        int posVal = 0;
        try {
            posVal = std::stoi(posStr);
        } catch (...) {
            // pass
            out << line << "\n";
            continue;
        }
        // if regionChrom not empty => check if chrom matches
        // if regionChrom empty => we do everything
        if (!regionChrom.empty()) {
            if (chrom != regionChrom) {
                // not in region => pass line
                out << line << "\n";
                continue;
            }
            if (posVal < regionStart || posVal > regionEnd) {
                // out of range
                out << line << "\n";
                continue;
            }
        }
        // we keep this variant => parse genotype codes
        LDVariant v;
        v.chrom = chrom;
        v.pos = posVal;
        v.id = fields[2];
        v.genotype.resize(numSamples, -1);
        // sample columns => fields[9..]
        int sampleCol = 9;
        for (int s = 0; s < numSamples; s++) {
            if ((size_t)(sampleCol + s) >= fields.size())
                break;
            v.genotype[s] = parseGenotype(fields[sampleCol + s]);
        }
        variants.push_back(std::move(v));
        // pass line unchanged
        out << line << "\n";
    }

    // after reading all => we compute pairwise r^2
    // if we have M variants => produce an MxM matrix
    // we do: #LD_MATRIX_START, then we print M lines
    // first line => #VariantIndices ...
    // or we print a table with variant index as row, col + r^2
    // We'll do a simple text approach
    size_t M = variants.size();
    if (M < 2) {
        out << "#LD_MATRIX_START\n";
        out << "No or only one variant in the region => no pairwise LD.\n";
        out << "#LD_MATRIX_END\n";
        return;
    }

    out << "#LD_MATRIX_START\n";
    // Print header row => e.g. row0col0 is "", row0col i => i. We'll also print "Chrom:Pos" as col
    // first line => tab, var1, var2, ...
    out << "Index/Var";
    for (size_t j = 0; j < M; j++) {
        out << "\t" << variants[j].chrom << ":" << variants[j].pos;
    }
    out << "\n";

    // each row => i => row header => i-chrom:pos => r2 vs all j
    for (size_t i = 0; i < M; i++) {
        out << variants[i].chrom << ":" << variants[i].pos;
        for (size_t j = 0; j < M; j++) {
            if (j < i) {
                // symmetrical => we can do the same as [j][i]
                // but let's just compute anyway or store
                double r2 = computeRsq(variants[i].genotype, variants[j].genotype);
                out << "\t" << std::fixed << std::setprecision(4) << r2;
            } else if (i == j) {
                // r2 with self => 1.0
                out << "\t1.0000";
            } else {
                // i< j => compute
                double r2 = computeRsq(variants[i].genotype, variants[j].genotype);
                out << "\t" << std::fixed << std::setprecision(4) << r2;
            }
        }
        out << "\n";
    }

    out << "#LD_MATRIX_END\n";
}

// ----------------------------------------------------------------------
// main run
// ----------------------------------------------------------------------
int VCFXLDCalculator::run(int argc, char *argv[]) {
    static struct option long_opts[] = {
        {"help", no_argument, 0, 'h'},
        {"region", required_argument, 0, 'r'},
        {"streaming", no_argument, 0, 's'},
        {"matrix", no_argument, 0, 'm'},
        {"window", required_argument, 0, 'w'},
        {"threshold", required_argument, 0, 't'},
        {0, 0, 0, 0}
    };
    bool showHelp = false;
    std::string regionStr;

    while (true) {
        int c = getopt_long(argc, argv, "hr:smw:t:", long_opts, NULL);
        if (c == -1)
            break;
        switch (c) {
        case 'h':
            showHelp = true;
            break;
        case 'r':
            regionStr = optarg;
            break;
        case 's':
            streamingMode = true;
            matrixMode = false;
            break;
        case 'm':
            matrixMode = true;
            streamingMode = false;
            break;
        case 'w':
            try {
                windowSize = std::stoul(optarg);
                if (windowSize == 0) windowSize = 1;
            } catch (...) {
                std::cerr << "Error: Invalid window size '" << optarg << "'\n";
                return 1;
            }
            break;
        case 't':
            try {
                ldThreshold = std::stod(optarg);
                if (ldThreshold < 0.0) ldThreshold = 0.0;
                if (ldThreshold > 1.0) ldThreshold = 1.0;
            } catch (...) {
                std::cerr << "Error: Invalid threshold '" << optarg << "'\n";
                return 1;
            }
            break;
        default:
            showHelp = true;
        }
    }
    if (showHelp) {
        displayHelp();
        return 0;
    }

    // parse region if any
    std::string regionChrom;
    int regionStart = 0, regionEnd = 0;
    if (!regionStr.empty()) {
        if (!parseRegion(regionStr, regionChrom, regionStart, regionEnd)) {
            std::cerr << "Error parsing region '" << regionStr << "'. Use e.g. chr1:10000-20000\n";
            return 1;
        }
    }

    // Dispatch to appropriate mode
    // Default is streaming mode (efficient O(M*window) for large files)
    // Use --matrix for backward-compatible O(M²) matrix output
    if (matrixMode) {
        computeLD(std::cin, std::cout, regionChrom, regionStart, regionEnd);
    } else {
        computeLDStreaming(std::cin, std::cout, regionChrom, regionStart, regionEnd, windowSize, ldThreshold);
    }
    return 0;
}

static void show_help() {
    VCFXLDCalculator obj;
    char arg0[] = "VCFX_ld_calculator";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_ld_calculator", show_help))
        return 0;
    VCFXLDCalculator calc;
    return calc.run(argc, argv);
}
