#include "VCFX_diff_tool.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

// ----------------------------------------------------------------------
// A helper function to split a string by a delimiter (optimized)
// ----------------------------------------------------------------------
static std::vector<std::string> split(const std::string &str, char delim) {
    std::vector<std::string> tokens;
    tokens.reserve(8);
    size_t start = 0;
    size_t end;
    while ((end = str.find(delim, start)) != std::string::npos) {
        tokens.emplace_back(str, start, end - start);
        start = end + 1;
    }
    tokens.emplace_back(str, start);
    return tokens;
}

// ----------------------------------------------------------------------
// displayHelp
// ----------------------------------------------------------------------
void VCFXDiffTool::displayHelp() {
    std::cout << "VCFX_diff_tool: Compare two VCF files and identify differences.\n\n"
              << "Usage:\n"
              << "  VCFX_diff_tool --file1 <file1.vcf> --file2 <file2.vcf> [options]\n\n"
              << "Options:\n"
              << "  -h, --help                Display this help message and exit\n"
              << "  -a, --file1 <file1.vcf>   Specify the first VCF file\n"
              << "  -b, --file2 <file2.vcf>   Specify the second VCF file\n"
              << "  -s, --assume-sorted       Assume inputs are sorted by (CHROM, POS).\n"
              << "                            Enables streaming mode with O(1) memory.\n"
              << "  -n, --natural-chr         Use natural chromosome ordering (chr1 < chr2 < chr10)\n\n"
              << "Modes:\n"
              << "  Default mode:     Loads both files into memory (works with unsorted files)\n"
              << "  Streaming mode:   Two-pointer merge diff with O(1) memory (requires sorted input)\n\n"
              << "Example:\n"
              << "  VCFX_diff_tool --file1 file1.vcf --file2 file2.vcf\n"
              << "  VCFX_diff_tool -a sorted1.vcf -b sorted2.vcf --assume-sorted\n";
}

// ----------------------------------------------------------------------
// generateVariantKey
//   - altField might be e.g. "G" or "G,T" or "T,GA"
//   - We split on commas, sort them, rejoin => e.g. "G,T"
//   - Then produce "chrom:pos:ref:thatSortedAlt"
// ----------------------------------------------------------------------
std::string VCFXDiffTool::generateVariantKey(const std::string &chrom, const std::string &pos, const std::string &ref,
                                             const std::string &altField) {
    auto alts = split(altField, ',');
    // sort them lexicographically
    std::sort(alts.begin(), alts.end());
    // rejoin
    std::ostringstream altSS;
    for (size_t i = 0; i < alts.size(); ++i) {
        if (i > 0)
            altSS << ",";
        altSS << alts[i];
    }
    std::string sortedAlt = altSS.str();

    return chrom + ":" + pos + ":" + ref + ":" + sortedAlt;
}

// ----------------------------------------------------------------------
// parseVCFLine - Extract key components from a VCF data line
// Returns false if line is a header or invalid
// ----------------------------------------------------------------------
bool VCFXDiffTool::parseVCFLine(const std::string &line, std::string &chrom, long &pos,
                                std::string &ref, std::string &alt, std::string &key) {
    if (line.empty() || line[0] == '#') {
        return false;
    }

    std::vector<std::string> fields;
    vcfx::split_tabs(line, fields);

    if (fields.size() < 5) {
        return false;
    }

    chrom = fields[0];
    try {
        pos = std::stol(fields[1]);
    } catch (...) {
        return false;
    }
    ref = fields[3];
    alt = fields[4];

    // Normalize alt field by sorting alleles
    auto alts = split(alt, ',');
    std::sort(alts.begin(), alts.end());
    std::string sortedAlt;
    for (size_t i = 0; i < alts.size(); ++i) {
        if (i > 0) sortedAlt += ",";
        sortedAlt += alts[i];
    }
    alt = sortedAlt;

    key = chrom + ":" + fields[1] + ":" + ref + ":" + alt;
    return true;
}

// ----------------------------------------------------------------------
// compareKeys - Compare two variant keys for sorting order
// Returns <0 if a < b, 0 if a == b, >0 if a > b
// ----------------------------------------------------------------------
int VCFXDiffTool::compareKeys(const std::string &chromA, long posA, const std::string &refA, const std::string &altA,
                              const std::string &chromB, long posB, const std::string &refB, const std::string &altB) {
    // Compare chromosomes
    int chromCmp;
    if (naturalChromOrder) {
        // Natural chromosome ordering (chr1 < chr2 < chr10)
        std::string prefixA, prefixB, suffixA, suffixB;
        long numA = -1, numB = -1;

        // Parse chromosome A
        std::string chrA = chromA;
        if (chrA.size() >= 3 && (chrA.substr(0, 3) == "chr" || chrA.substr(0, 3) == "Chr" || chrA.substr(0, 3) == "CHR")) {
            prefixA = chrA.substr(0, 3);
            chrA = chrA.substr(3);
        }
        size_t idxA = 0;
        while (idxA < chrA.size() && std::isdigit(chrA[idxA])) idxA++;
        if (idxA > 0) numA = std::stol(chrA.substr(0, idxA));
        suffixA = chrA.substr(idxA);

        // Parse chromosome B
        std::string chrB = chromB;
        if (chrB.size() >= 3 && (chrB.substr(0, 3) == "chr" || chrB.substr(0, 3) == "Chr" || chrB.substr(0, 3) == "CHR")) {
            prefixB = chrB.substr(0, 3);
            chrB = chrB.substr(3);
        }
        size_t idxB = 0;
        while (idxB < chrB.size() && std::isdigit(chrB[idxB])) idxB++;
        if (idxB > 0) numB = std::stol(chrB.substr(0, idxB));
        suffixB = chrB.substr(idxB);

        // Compare prefix
        chromCmp = prefixA.compare(prefixB);
        if (chromCmp != 0) return chromCmp;

        // Compare numeric part (no number = comes after numbers)
        if (numA >= 0 && numB >= 0) {
            if (numA != numB) return (numA < numB) ? -1 : 1;
        } else if (numA >= 0) {
            return -1;  // A has number, B doesn't
        } else if (numB >= 0) {
            return 1;   // B has number, A doesn't
        }

        // Compare suffix
        chromCmp = suffixA.compare(suffixB);
        if (chromCmp != 0) return chromCmp;
    } else {
        // Lexicographic ordering
        chromCmp = chromA.compare(chromB);
        if (chromCmp != 0) return chromCmp;
    }

    // Same chromosome, compare position
    if (posA != posB) return (posA < posB) ? -1 : 1;

    // Same position, compare ref
    int refCmp = refA.compare(refB);
    if (refCmp != 0) return refCmp;

    // Same ref, compare alt
    return altA.compare(altB);
}

// ----------------------------------------------------------------------
// loadVariants (for in-memory mode)
//   - Read each line, skip headers (#...)
//   - parse CHROM, POS, ID, REF, ALT
//   - generate key => add to set
// ----------------------------------------------------------------------
bool VCFXDiffTool::loadVariants(const std::string &filePath, std::unordered_set<std::string> &variants) {
    std::ifstream infile(filePath);
    if (!infile.is_open()) {
        std::cerr << "Error: Unable to open file " << filePath << "\n";
        return false;
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') {
            // skip headers and empty lines
            continue;
        }

        std::stringstream ss(line);
        std::string chrom, pos, id, ref, alt;

        // minimal parse
        if (!(ss >> chrom >> pos >> id >> ref >> alt)) {
            std::cerr << "Warning: Skipping invalid VCF line:\n" << line << "\n";
            continue;
        }

        std::string key = generateVariantKey(chrom, pos, ref, alt);
        variants.insert(key);
    }
    return true;
}

// ----------------------------------------------------------------------
// diffInMemory - Original algorithm using hash sets
// Loads both files into memory, finds set differences
// ----------------------------------------------------------------------
void VCFXDiffTool::diffInMemory(const std::string &file1Path, const std::string &file2Path) {
    std::unordered_set<std::string> variantsFile1;
    std::unordered_set<std::string> variantsFile2;

    if (!loadVariants(file1Path, variantsFile1)) {
        std::cerr << "Error: Failed to load variants from " << file1Path << "\n";
        return;
    }
    if (!loadVariants(file2Path, variantsFile2)) {
        std::cerr << "Error: Failed to load variants from " << file2Path << "\n";
        return;
    }

    // Identify differences
    std::vector<std::string> uniqueToFile1;
    std::vector<std::string> uniqueToFile2;

    for (const auto &v : variantsFile1) {
        if (variantsFile2.find(v) == variantsFile2.end()) {
            uniqueToFile1.push_back(v);
        }
    }
    for (const auto &v : variantsFile2) {
        if (variantsFile1.find(v) == variantsFile1.end()) {
            uniqueToFile2.push_back(v);
        }
    }

    // Print results
    std::cout << "Variants unique to " << file1Path << ":\n";
    for (auto &v : uniqueToFile1) {
        std::cout << v << "\n";
    }
    std::cout << "\nVariants unique to " << file2Path << ":\n";
    for (auto &v : uniqueToFile2) {
        std::cout << v << "\n";
    }
}

// ----------------------------------------------------------------------
// diffStreaming - Two-pointer merge diff with O(1) memory
// Requires sorted input files
// ----------------------------------------------------------------------
void VCFXDiffTool::diffStreaming(const std::string &file1Path, const std::string &file2Path) {
    std::ifstream file1(file1Path);
    std::ifstream file2(file2Path);

    if (!file1.is_open()) {
        std::cerr << "Error: Unable to open file " << file1Path << "\n";
        return;
    }
    if (!file2.is_open()) {
        std::cerr << "Error: Unable to open file " << file2Path << "\n";
        return;
    }

    std::vector<std::string> uniqueToFile1;
    std::vector<std::string> uniqueToFile2;

    std::string line1, line2;
    std::string chrom1, chrom2, ref1, ref2, alt1, alt2, key1, key2;
    long pos1 = 0, pos2 = 0;
    bool have1 = false, have2 = false;

    // Helper to read next valid variant from file1
    auto readNext1 = [&]() -> bool {
        while (std::getline(file1, line1)) {
            if (parseVCFLine(line1, chrom1, pos1, ref1, alt1, key1)) {
                return true;
            }
        }
        return false;
    };

    // Helper to read next valid variant from file2
    auto readNext2 = [&]() -> bool {
        while (std::getline(file2, line2)) {
            if (parseVCFLine(line2, chrom2, pos2, ref2, alt2, key2)) {
                return true;
            }
        }
        return false;
    };

    // Initialize
    have1 = readNext1();
    have2 = readNext2();

    // Two-pointer merge
    while (have1 && have2) {
        int cmp = compareKeys(chrom1, pos1, ref1, alt1, chrom2, pos2, ref2, alt2);

        if (cmp < 0) {
            // Variant in file1 but not file2
            uniqueToFile1.push_back(key1);
            have1 = readNext1();
        } else if (cmp > 0) {
            // Variant in file2 but not file1
            uniqueToFile2.push_back(key2);
            have2 = readNext2();
        } else {
            // Same variant in both files - advance both
            have1 = readNext1();
            have2 = readNext2();
        }
    }

    // Drain remaining variants from file1
    while (have1) {
        uniqueToFile1.push_back(key1);
        have1 = readNext1();
    }

    // Drain remaining variants from file2
    while (have2) {
        uniqueToFile2.push_back(key2);
        have2 = readNext2();
    }

    // Print results (same format as in-memory mode for compatibility)
    std::cout << "Variants unique to " << file1Path << ":\n";
    for (const auto &v : uniqueToFile1) {
        std::cout << v << "\n";
    }
    std::cout << "\nVariants unique to " << file2Path << ":\n";
    for (const auto &v : uniqueToFile2) {
        std::cout << v << "\n";
    }
}

// ----------------------------------------------------------------------
// run()
//   - parse args
//   - dispatch to appropriate diff algorithm
// ----------------------------------------------------------------------
int VCFXDiffTool::run(int argc, char *argv[]) {
    int opt;
    bool showHelp = false;
    std::string file1Path;
    std::string file2Path;

    static struct option long_options[] = {{"help", no_argument, 0, 'h'},
                                           {"file1", required_argument, 0, 'a'},
                                           {"file2", required_argument, 0, 'b'},
                                           {"assume-sorted", no_argument, 0, 's'},
                                           {"natural-chr", no_argument, 0, 'n'},
                                           {0, 0, 0, 0}};

    while ((opt = getopt_long(argc, argv, "ha:b:sn", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
            break;
        case 'a':
            file1Path = optarg;
            break;
        case 'b':
            file2Path = optarg;
            break;
        case 's':
            assumeSorted = true;
            break;
        case 'n':
            naturalChromOrder = true;
            break;
        default:
            showHelp = true;
        }
    }

    if (showHelp || file1Path.empty() || file2Path.empty()) {
        displayHelp();
        return (showHelp ? 0 : 1);
    }

    // Verify files exist before processing
    std::ifstream test1(file1Path);
    if (!test1.is_open()) {
        std::cerr << "Error: Unable to open file " << file1Path << "\n";
        return 1;
    }
    test1.close();

    std::ifstream test2(file2Path);
    if (!test2.is_open()) {
        std::cerr << "Error: Unable to open file " << file2Path << "\n";
        return 1;
    }
    test2.close();

    // Dispatch to appropriate algorithm
    if (assumeSorted) {
        diffStreaming(file1Path, file2Path);
    } else {
        diffInMemory(file1Path, file2Path);
    }

    return 0;
}

// ----------------------------------------------------------------------
// main
// ----------------------------------------------------------------------
static void show_help() {
    VCFXDiffTool obj;
    char arg0[] = "VCFX_diff_tool";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_diff_tool", show_help))
        return 0;
    VCFXDiffTool diffTool;
    return diffTool.run(argc, argv);
}
