#include "VCFX_fasta_converter.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <map>
#include <sstream>
#include <unistd.h>
#include <unordered_map>

// IUPAC ambiguity lookup table - O(1) array lookup instead of map
// Index by (base1_index * 4 + base2_index) where A=0, C=1, G=2, T=3
static const char IUPAC_table[16] = {
    'A', 'M', 'R', 'W',  // A+A, A+C, A+G, A+T
    'M', 'C', 'S', 'Y',  // C+A, C+C, C+G, C+T
    'R', 'S', 'G', 'K',  // G+A, G+C, G+G, G+T
    'W', 'Y', 'K', 'T'   // T+A, T+C, T+G, T+T
};

// Map base to index (A=0, C=1, G=2, T=3, else -1)
static inline int baseToIndex(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

// Fast IUPAC combination - O(1)
static inline char combineBasesIUPAC(char b1, char b2) {
    int i1 = baseToIndex(b1);
    int i2 = baseToIndex(b2);
    if (i1 < 0 || i2 < 0) return 'N';
    return IUPAC_table[i1 * 4 + i2];
}

// Utility to convert allele index to base - returns 'N' on failure
static inline char alleleIndexToBase(int alleleIndex, const std::string &ref,
                                      const std::vector<std::string> &altAlleles) {
    if (alleleIndex == 0) {
        return (ref.size() == 1) ? static_cast<char>(std::toupper(ref[0])) : 'N';
    }
    int altPos = alleleIndex - 1;
    if (altPos < 0 || static_cast<size_t>(altPos) >= altAlleles.size()) {
        return 'N';
    }
    const std::string &a = altAlleles[altPos];
    return (a.size() == 1) ? static_cast<char>(std::toupper(a[0])) : 'N';
}

int VCFXFastaConverter::run(int argc, char *argv[]) {
    int opt;
    bool showHelp = false;

    static struct option long_options[] = {{"help", no_argument, 0, 'h'}, {0, 0, 0, 0}};

    while ((opt = getopt_long(argc, argv, "h", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'h':
        default:
            showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    convertVCFtoFasta(std::cin, std::cout);
    return 0;
}

void VCFXFastaConverter::displayHelp() {
    std::cout << "VCFX_fasta_converter: Convert a variant-only VCF into simple per-sample FASTA.\n\n"
              << "Usage:\n"
              << "  VCFX_fasta_converter [options] < input.vcf > output.fasta\n\n"
              << "Description:\n"
              << "  Reads a VCF with diploid genotypes and writes a FASTA file. Each variant\n"
              << "  line becomes one position in the FASTA alignment. For multi-allelic sites,\n"
              << "  each sample's genotype is interpreted to produce a single IUPAC base\n"
              << "  (if heterozygous with different single-base alleles) or 'N' if ambiguous.\n\n"
              << "  Indels, multi-base alleles, or complicated genotypes default to 'N'.\n\n"
              << "Example:\n"
              << "  VCFX_fasta_converter < input.vcf > output.fasta\n\n";
}

// Split string by delimiter - optimized direct parsing
static inline void splitByDelim(const std::string &s, char d, std::vector<std::string> &out) {
    out.clear();
    size_t start = 0, end;
    while ((end = s.find(d, start)) != std::string::npos) {
        out.emplace_back(s, start, end - start);
        start = end + 1;
    }
    out.emplace_back(s, start);
}

// Parse genotype from sample field and return IUPAC base
static inline char parseGenotypeToBase(const char* sampleData, size_t len, int gtIndex,
                                        const std::string &ref,
                                        const std::vector<std::string> &altAlleles) {
    // Fast path: GT is first field (most common case)
    size_t gtStart = 0, gtEnd = len;

    if (gtIndex == 0) {
        // Find first colon
        for (size_t i = 0; i < len; ++i) {
            if (sampleData[i] == ':') {
                gtEnd = i;
                break;
            }
        }
    } else {
        // Skip to gtIndex-th field
        int idx = 0;
        size_t pos = 0;
        while (idx < gtIndex && pos < len) {
            if (sampleData[pos] == ':') {
                idx++;
                if (idx == gtIndex) {
                    gtStart = pos + 1;
                }
            }
            pos++;
        }
        if (idx != gtIndex) return 'N';
        // Find end of GT field
        for (size_t i = gtStart; i < len; ++i) {
            if (sampleData[i] == ':') {
                gtEnd = i;
                break;
            }
        }
    }

    if (gtStart >= gtEnd || sampleData[gtStart] == '.') {
        return 'N';
    }

    // Parse diploid genotype inline
    int a1 = 0, a2 = 0;
    size_t pos = gtStart;

    // Parse first allele
    while (pos < gtEnd && sampleData[pos] != '/' && sampleData[pos] != '|') {
        char c = sampleData[pos];
        if (c == '.') return 'N';
        if (c < '0' || c > '9') return 'N';
        a1 = a1 * 10 + (c - '0');
        pos++;
    }

    if (pos >= gtEnd) return 'N'; // No separator found
    pos++; // Skip separator

    // Parse second allele
    while (pos < gtEnd) {
        char c = sampleData[pos];
        if (c == '.') return 'N';
        if (c < '0' || c > '9') break;
        a2 = a2 * 10 + (c - '0');
        pos++;
    }

    char b1 = alleleIndexToBase(a1, ref, altAlleles);
    char b2 = alleleIndexToBase(a2, ref, altAlleles);

    if (b1 == 'N' || b2 == 'N') return 'N';
    return (b1 == b2) ? b1 : combineBasesIUPAC(b1, b2);
}

// OPTIMIZED: Memory-efficient temp file approach
// Instead of O(V×S) memory (~1GB for 427K×2504), uses O(S) memory (~3KB)
// Writes column-major (all samples per variant) to temp file
// Reads row-major (all variants per sample) for output
void VCFXFastaConverter::convertVCFtoFasta(std::istream &in, std::ostream &out) {
    std::string line;
    std::vector<std::string> sampleNames;

    bool headerParsed = false;
    size_t numSamples = 0;
    size_t variantCount = 0;

    // Reusable buffers for parsing
    std::vector<std::string> fields;
    fields.reserve(16);
    std::vector<std::string> altAlleles;
    altAlleles.reserve(4);
    std::vector<std::string> formatFields;
    formatFields.reserve(8);

    // Create temp file for column-major storage
    char tempPath[] = "/tmp/vcfx_fasta_XXXXXX";
    int tempFd = mkstemp(tempPath);
    if (tempFd < 0) {
        std::cerr << "Error: Cannot create temp file\n";
        return;
    }

    // Buffer for variant bases (one char per sample)
    std::vector<char> variantBases;

    // PHASE 1: Single pass through VCF, write bases to temp file column-major
    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            if (line.rfind("#CHROM", 0) == 0) {
                std::vector<std::string> headers;
                vcfx::split_tabs(line, headers);
                for (size_t i = 9; i < headers.size(); ++i) {
                    sampleNames.push_back(headers[i]);
                }
                numSamples = sampleNames.size();
                variantBases.resize(numSamples);
                headerParsed = true;
            }
            continue;
        }

        if (!headerParsed) {
            std::cerr << "Error: #CHROM header not found before data lines.\n";
            close(tempFd);
            unlink(tempPath);
            return;
        }

        vcfx::split_tabs(line, fields);
        if (fields.size() < 9 + numSamples) {
            std::cerr << "Warning: Skipping malformed VCF line with insufficient columns.\n";
            continue;
        }

        const std::string &ref = fields[3];
        const std::string &altField = fields[4];
        const std::string &format = fields[8];

        // Parse alt alleles
        splitByDelim(altField, ',', altAlleles);

        // Find GT index in format
        splitByDelim(format, ':', formatFields);
        int gtIndex = -1;
        for (size_t i = 0; i < formatFields.size(); ++i) {
            if (formatFields[i] == "GT") {
                gtIndex = static_cast<int>(i);
                break;
            }
        }

        // Process each sample
        for (size_t s = 0; s < numSamples; ++s) {
            if (gtIndex < 0) {
                variantBases[s] = 'N';
            } else {
                const std::string &sampleData = fields[9 + s];
                variantBases[s] = parseGenotypeToBase(sampleData.c_str(), sampleData.size(),
                                                       gtIndex, ref, altAlleles);
            }
        }

        // Write all sample bases for this variant (column-major)
        ssize_t written = write(tempFd, variantBases.data(), numSamples);
        if (written != static_cast<ssize_t>(numSamples)) {
            std::cerr << "Error: Failed to write to temp file\n";
            close(tempFd);
            unlink(tempPath);
            return;
        }
        variantCount++;
    }

    // PHASE 2: Read row-major (transpose) and output FASTA
    // For each sample, read all their bases and output as FASTA

    // Buffer for reading
    const size_t BUFFER_SIZE = 65536;  // 64KB read buffer
    std::vector<char> readBuffer(BUFFER_SIZE);
    std::vector<char> sampleSeq(variantCount);

    for (size_t s = 0; s < numSamples; ++s) {
        // Output FASTA header
        out << ">" << sampleNames[s] << "\n";

        // Read this sample's bases from temp file (scattered reads)
        // Each variant has numSamples bytes, sample s is at offset s within each
        for (size_t v = 0; v < variantCount; ++v) {
            off_t offset = static_cast<off_t>(v * numSamples + s);
            if (lseek(tempFd, offset, SEEK_SET) < 0) {
                sampleSeq[v] = 'N';
                continue;
            }
            char c;
            if (read(tempFd, &c, 1) == 1) {
                sampleSeq[v] = c;
            } else {
                sampleSeq[v] = 'N';
            }
        }

        // Output in 60-char lines
        for (size_t i = 0; i < variantCount; i += 60) {
            size_t len = std::min(size_t(60), variantCount - i);
            out.write(sampleSeq.data() + i, static_cast<std::streamsize>(len));
            out << "\n";
        }
    }

    // Cleanup
    close(tempFd);
    unlink(tempPath);
}

static void show_help() {
    VCFXFastaConverter obj;
    char arg0[] = "VCFX_fasta_converter";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();
    if (vcfx::handle_common_flags(argc, argv, "VCFX_fasta_converter", show_help))
        return 0;
    VCFXFastaConverter app;
    return app.run(argc, argv);
}
