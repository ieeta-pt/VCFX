#include "vcfx_core.h"
#include "vcfx_io.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

// ---------------------------------------------------------------------
// Structures and Declarations
// ---------------------------------------------------------------------
struct AlleleCounterArguments {
    std::vector<std::string> samples;
};

static void printHelp();
static bool parseArguments(int argc, char *argv[], AlleleCounterArguments &args);
static bool countAlleles(std::istream &in, std::ostream &out, const AlleleCounterArguments &args);

// ---------------------------------------------------------------------
// printHelp
// ---------------------------------------------------------------------
static void printHelp() {
    std::cout << "VCFX_allele_counter\n"
              << "Usage: VCFX_allele_counter [OPTIONS] < input.vcf > allele_counts.tsv\n\n"
              << "Options:\n"
              << "  --samples, -s \"Sample1 Sample2\"   Specify the sample names to include.\n"
              << "                                     If not specified, all samples are processed.\n"
              << "  --help, -h                        Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Reads a VCF from stdin and outputs a TSV file with the columns:\n"
              << "    CHROM  POS  ID  REF  ALT  Sample  Ref_Count  Alt_Count\n\n"
              << "  Each sample for each variant is listed. Alleles are determined by\n"
              << "  genotype strings (GT). This code treats any numeric allele that is\n"
              << "  not '0' as ALT, e.g. '1' or '2' or '3' => alt.\n\n"
              << "Examples:\n"
              << "  1) Count alleles for SampleA and SampleB:\n"
              << "     VCFX_allele_counter --samples \"SampleA SampleB\" < input.vcf > allele_counts.tsv\n\n"
              << "  2) Count alleles for all samples:\n"
              << "     VCFX_allele_counter < input.vcf > allele_counts_all.tsv\n";
}

// ---------------------------------------------------------------------
// parseArguments
// ---------------------------------------------------------------------
static bool parseArguments(int argc, char *argv[], AlleleCounterArguments &args) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--samples" || arg == "-s") && i + 1 < argc) {
            std::string samplesStr = argv[++i];
            // Split by space - inline
            size_t start = 0, end;
            while ((end = samplesStr.find(' ', start)) != std::string::npos) {
                if (end > start) {
                    args.samples.emplace_back(samplesStr, start, end - start);
                }
                start = end + 1;
            }
            if (start < samplesStr.size()) {
                args.samples.emplace_back(samplesStr, start);
            }
            // Trim whitespace from each sample name
            for (auto &sample : args.samples) {
                sample.erase(0, sample.find_first_not_of(" \t\n\r\f\v"));
                sample.erase(sample.find_last_not_of(" \t\n\r\f\v") + 1);
            }
        } else if (arg == "--help" || arg == "-h") {
            printHelp();
            std::exit(0);
        } else {
            std::cerr << "Warning: Unrecognized argument '" << arg << "'.\n";
        }
    }
    return true;
}

// OPTIMIZED: Zero-allocation genotype parsing
// Returns ref and alt counts directly from const char* without creating vectors
static inline void parseGenotypeInline(const char* gt, size_t gtLen, int &refCount, int &altCount) {
    refCount = 0;
    altCount = 0;
    size_t pos = 0;

    while (pos < gtLen) {
        // Skip separators and find start of allele
        while (pos < gtLen && (gt[pos] == '/' || gt[pos] == '|')) {
            pos++;
        }
        if (pos >= gtLen) break;

        // Parse allele - either a digit sequence or '.'
        if (gt[pos] == '.') {
            // Missing allele - skip
            pos++;
            continue;
        }

        // Parse numeric allele
        int allele = 0;
        bool hasDigit = false;
        while (pos < gtLen && gt[pos] >= '0' && gt[pos] <= '9') {
            allele = allele * 10 + (gt[pos] - '0');
            hasDigit = true;
            pos++;
        }

        if (hasDigit) {
            if (allele == 0) {
                refCount++;
            } else {
                altCount++;
            }
        }
    }
}

// ---------------------------------------------------------------------
// countAlleles - OPTIMIZED with buffered output and zero-allocation parsing
// ---------------------------------------------------------------------
static bool countAlleles(std::istream &in, std::ostream &out, const AlleleCounterArguments &args) {
    // Set up output buffering for high-volume output (~1.07 billion lines for 427K×2504)
    static char outBuffer[1024 * 1024];  // 1MB buffer
    out.rdbuf()->pubsetbuf(outBuffer, sizeof(outBuffer));

    std::string line;
    std::vector<std::string> headerFields;
    bool foundChromHeader = false;

    std::vector<int> sampleIndices;

    out << "CHROM\tPOS\tID\tREF\tALT\tSample\tRef_Count\tAlt_Count\n";

    // Reusable vector for tab splitting
    std::vector<std::string> fields;
    fields.reserve(2600);  // Pre-allocate for ~2504 samples + 9 fixed fields

    // Output line buffer to reduce I/O syscalls
    std::string outLine;
    outLine.reserve(256);

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            if (line.rfind("#CHROM", 0) == 0) {
                vcfx::split_tabs(line, headerFields);
                if (headerFields.size() < 9) {
                    std::cerr << "Error: #CHROM line has fewer than 9 columns.\n";
                    return false;
                }

                std::unordered_map<std::string, int> sampleMap;
                for (size_t i = 9; i < headerFields.size(); ++i) {
                    sampleMap[headerFields[i]] = static_cast<int>(i);
                }

                if (!args.samples.empty()) {
                    for (const auto &s : args.samples) {
                        auto it = sampleMap.find(s);
                        if (it == sampleMap.end()) {
                            std::cerr << "Error: Sample '" << s << "' not found in VCF header.\n";
                            return false;
                        }
                        sampleIndices.push_back(it->second);
                    }
                } else {
                    for (size_t i = 9; i < headerFields.size(); ++i) {
                        sampleIndices.push_back(static_cast<int>(i));
                    }
                }
                foundChromHeader = true;
            }
            continue;
        }

        if (!foundChromHeader) {
            std::cerr << "Error: VCF #CHROM header not found before records.\n";
            return false;
        }

        vcfx::split_tabs(line, fields);
        if (fields.size() < 9) {
            std::cerr << "Warning: Skipping invalid VCF line with fewer than 9 fields.\n";
            continue;
        }

        // Cache field references to avoid repeated indexing
        const std::string &chrom = fields[0];
        const std::string &pos = fields[1];
        const std::string &id = fields[2];
        const std::string &ref = fields[3];
        const std::string &alt = fields[4];

        // Build common prefix once per variant
        outLine.clear();
        outLine.append(chrom);
        outLine.push_back('\t');
        outLine.append(pos);
        outLine.push_back('\t');
        outLine.append(id);
        outLine.push_back('\t');
        outLine.append(ref);
        outLine.push_back('\t');
        outLine.append(alt);
        outLine.push_back('\t');
        size_t prefixLen = outLine.size();

        for (int sIndex : sampleIndices) {
            if (sIndex >= static_cast<int>(fields.size())) continue;

            const std::string &sampleField = fields[sIndex];

            // Find GT field end (first colon or end of string)
            size_t gtEnd = sampleField.find(':');
            if (gtEnd == std::string::npos) gtEnd = sampleField.size();

            // Zero-allocation genotype parsing
            int refCount = 0, altCount = 0;
            parseGenotypeInline(sampleField.c_str(), gtEnd, refCount, altCount);

            // Build output line
            outLine.resize(prefixLen);
            outLine.append(headerFields[sIndex]);
            outLine.push_back('\t');

            // Fast integer to string conversion
            char numBuf[16];
            int len = snprintf(numBuf, sizeof(numBuf), "%d\t%d\n", refCount, altCount);
            outLine.append(numBuf, len);

            out.write(outLine.data(), outLine.size());
        }
    }
    return true;
}

// ---------------------------------------------------------------------
// main()
// ---------------------------------------------------------------------
static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();
    if (vcfx::handle_common_flags(argc, argv, "VCFX_allele_counter", show_help))
        return 0;

    AlleleCounterArguments args;
    parseArguments(argc, argv, args);

    if (!args.samples.empty()) {
        std::cerr << "Info: Counting alleles for these samples:\n";
        for (auto &s : args.samples) {
            std::cerr << "  " << s << "\n";
        }
    } else {
        std::cerr << "Info: Counting alleles for ALL samples.\n";
    }

    bool success = countAlleles(std::cin, std::cout, args);
    return success ? 0 : 1;
}
