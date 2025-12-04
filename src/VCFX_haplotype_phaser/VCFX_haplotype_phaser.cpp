#include "VCFX_haplotype_phaser.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cmath>
#include <deque>
#include <getopt.h>
#include <numeric>
#include <sstream>
#include <vector>

int VCFXHaplotypePhaser::run(int argc, char *argv[]) {
    // parse arguments
    int opt;
    bool showHelp = false;
    double ldThreshold = 0.8;
    bool streaming = false;
    size_t window = 1000;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"ld-threshold", required_argument, 0, 'l'},
        {"streaming", no_argument, 0, 's'},
        {"window", required_argument, 0, 'w'},
        {0, 0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "hl:sw:", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
            break;
        case 'l':
            try {
                ldThreshold = std::stod(optarg);
            } catch (...) {
                std::cerr << "Error: invalid LD threshold.\n";
                displayHelp();
                return 1;
            }
            break;
        case 's':
            streaming = true;
            break;
        case 'w':
            try {
                window = std::stoul(optarg);
            } catch (...) {
                std::cerr << "Error: invalid window size.\n";
                displayHelp();
                return 1;
            }
            break;
        default:
            showHelp = true;
        }
    }

    // If help was explicitly requested, show help and return success (0)
    if (showHelp) {
        displayHelp();
        return 0;
    }

    // If LD threshold is invalid, show help and return error (1)
    if (ldThreshold < 0.0 || ldThreshold > 1.0) {
        std::cerr << "Error: invalid LD threshold\n";
        displayHelp();
        return 1;
    }

    if (streaming) {
        phaseHaplotypesStreaming(std::cin, std::cout, ldThreshold, window);
    } else {
        phaseHaplotypes(std::cin, std::cout, ldThreshold);
    }
    return 0;
}

void VCFXHaplotypePhaser::displayHelp() {
    std::cout << "VCFX_haplotype_phaser: Group variants into blocks by naive LD threshold.\n\n"
              << "Usage:\n"
              << "  VCFX_haplotype_phaser [options] < input.vcf > blocks.txt\n\n"
              << "Options:\n"
              << "  -h, --help               Show this help message\n"
              << "  -l, --ld-threshold <val> r^2 threshold [0..1], default 0.8\n"
              << "  -s, --streaming          Enable streaming mode with sliding window.\n"
              << "                           Uses O(window * samples) memory instead of O(variants * samples).\n"
              << "  -w, --window <N>         Window size for streaming mode (default: 1000)\n\n"
              << "Performance:\n"
              << "  Default mode:   Loads all variants into memory, outputs blocks at end.\n"
              << "  Streaming mode: Uses sliding window, outputs blocks incrementally.\n"
              << "                  Enables processing of arbitrarily large files.\n\n"
              << "Example:\n"
              << "  VCFX_haplotype_phaser --ld-threshold 0.9 < input.vcf > blocks.txt\n"
              << "  VCFX_haplotype_phaser --streaming --window 500 < large.vcf > blocks.txt\n";
}

void VCFXHaplotypePhaser::phaseHaplotypes(std::istream &in, std::ostream &out, double ldThreshold) {
    std::string line;
    bool foundHeader = false;

    // We'll store the final list of variants
    std::vector<VariantData> variantList;
    std::vector<std::string> sampleNames;

    // Read header lines until first non-# line
    bool foundFirstVariant = false;
    std::string firstVariantLine;
    std::vector<std::string> tokens;
    tokens.reserve(16);

    while (!foundFirstVariant && std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '#') {
            // if it's #CHROM, parse sample columns
            if (!foundHeader && line.rfind("#CHROM", 0) == 0) {
                // parse
                vcfx::split_tabs(line, tokens);
                // from col=9 onward => samples
                for (size_t c = 9; c < tokens.size(); c++) {
                    sampleNames.push_back(tokens[c]);
                }
                foundHeader = true;
            }
            // print the header line out as well
            out << line << "\n";
        } else {
            // Found first non-header line
            foundFirstVariant = true;
            firstVariantLine = line;
        }
    }

    if (!foundHeader) {
        std::cerr << "Error: no #CHROM line found.\n";
        return;
    }

    int variantIndex = 0;
    std::vector<std::string> fields;
    fields.reserve(16);

    // Process variants
    auto processVariantLine = [&](const std::string &varLine) {
        if (varLine.empty() || varLine[0] == '#') {
            return;
        }

        vcfx::split_tabs(varLine, fields);

        if (fields.size() < 10) {
            std::cerr << "Warning: skipping line with <10 fields: " << varLine << "\n";
            return;
        }

        std::string chrom = fields[0];
        std::string posStr = fields[1];

        int posVal = 0;
        try {
            posVal = std::stoi(posStr);
        } catch (...) {
            std::cerr << "Warning: invalid pos => skip " << varLine << "\n";
            return;
        }

        // Build the genotype vector for this variant
        std::vector<int> genotypeVec;
        genotypeVec.reserve(fields.size() - 9);
        for (size_t s = 9; s < fields.size(); s++) {
            std::string gt = fields[s];
            size_t delim = gt.find_first_of("/|");
            if (delim == std::string::npos) {
                genotypeVec.push_back(-1);
                continue;
            }
            std::string a1 = gt.substr(0, delim);
            std::string a2 = gt.substr(delim + 1);
            if (a1.empty() || a2.empty() || a1 == "." || a2 == ".") {
                genotypeVec.push_back(-1);
                continue;
            }
            int i1 = 0, i2 = 0;
            try {
                i1 = std::stoi(a1);
                i2 = std::stoi(a2);
            } catch (...) {
                genotypeVec.push_back(-1);
                continue;
            }
            genotypeVec.push_back(i1 + i2);
        }

        VariantData v;
        v.chrom = chrom;
        v.pos = posVal;
        v.index = variantIndex++;
        v.genotype = genotypeVec;
        variantList.push_back(v);
    };

    // Process the first variant line if found during header processing
    if (foundFirstVariant) {
        processVariantLine(firstVariantLine);
    }

    // Process remaining variant lines
    while (std::getline(in, line)) {
        processVariantLine(line);
    }

    if (variantList.empty()) {
        std::cerr << "Error: no variant data found.\n";
        return;
    }

    // group variants
    auto blocks = groupVariants(variantList, ldThreshold);

    // output blocks
    out << "#HAPLOTYPE_BLOCKS_START\n";
    for (size_t b = 0; b < blocks.size(); b++) {
        out << "Block " << (b + 1) << ": ";
        for (size_t i = 0; i < blocks[b].size(); i++) {
            int idx = blocks[b][i];
            out << idx << ":(" << variantList[idx].chrom << ":" << variantList[idx].pos << ")";
            if (i + 1 < blocks[b].size())
                out << ", ";
        }
        out << "\n";
    }
    out << "#HAPLOTYPE_BLOCKS_END\n";
}

void VCFXHaplotypePhaser::phaseHaplotypesStreaming(std::istream &in, std::ostream &out,
                                                    double ldThreshold, size_t windowSize) {
    std::string line;
    bool foundHeader = false;
    std::vector<std::string> sampleNames;

    // Current block being built
    std::vector<VariantData> currentBlock;
    std::string currentChrom;
    int blockNumber = 0;
    int variantIndex = 0;

    // Output header marker
    bool headerMarkerWritten = false;
    std::vector<std::string> tokens;
    tokens.reserve(16);
    std::vector<std::string> fields;
    fields.reserve(16);

    // Read header lines until first non-# line
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '#') {
            if (!foundHeader && line.rfind("#CHROM", 0) == 0) {
                vcfx::split_tabs(line, tokens);
                for (size_t c = 9; c < tokens.size(); c++) {
                    sampleNames.push_back(tokens[c]);
                }
                foundHeader = true;
            }
            out << line << "\n";
            continue;
        }

        if (!foundHeader) {
            std::cerr << "Error: no #CHROM line found.\n";
            return;
        }

        // Write streaming header marker
        if (!headerMarkerWritten) {
            out << "#HAPLOTYPE_BLOCKS_START (streaming)\n";
            headerMarkerWritten = true;
        }

        // Parse variant line
        vcfx::split_tabs(line, fields);

        if (fields.size() < 10) {
            std::cerr << "Warning: skipping line with <10 fields\n";
            continue;
        }

        std::string chrom = fields[0];
        int posVal = 0;
        try {
            posVal = std::stoi(fields[1]);
        } catch (...) {
            std::cerr << "Warning: invalid pos => skip\n";
            continue;
        }

        // Build genotype vector
        std::vector<int> genotypeVec;
        genotypeVec.reserve(fields.size() - 9);
        for (size_t s = 9; s < fields.size(); s++) {
            std::string gt = fields[s];
            size_t delim = gt.find_first_of("/|");
            if (delim == std::string::npos) {
                genotypeVec.push_back(-1);
                continue;
            }
            std::string a1 = gt.substr(0, delim);
            std::string a2 = gt.substr(delim + 1);
            if (a1.empty() || a2.empty() || a1 == "." || a2 == ".") {
                genotypeVec.push_back(-1);
                continue;
            }
            int i1 = 0, i2 = 0;
            try {
                i1 = std::stoi(a1);
                i2 = std::stoi(a2);
            } catch (...) {
                genotypeVec.push_back(-1);
                continue;
            }
            genotypeVec.push_back(i1 + i2);
        }

        VariantData v;
        v.chrom = chrom;
        v.pos = posVal;
        v.index = variantIndex++;
        v.genotype = std::move(genotypeVec);

        // Streaming block logic
        if (currentBlock.empty()) {
            currentBlock.push_back(std::move(v));
            currentChrom = chrom;
        } else {
            // Check if chromosome changed - output current block and start new one
            if (chrom != currentChrom) {
                // Output current block
                blockNumber++;
                out << "Block " << blockNumber << ": ";
                for (size_t i = 0; i < currentBlock.size(); i++) {
                    out << currentBlock[i].index << ":(" << currentBlock[i].chrom << ":" << currentBlock[i].pos << ")";
                    if (i + 1 < currentBlock.size())
                        out << ", ";
                }
                out << "\n";

                // Start new block
                currentBlock.clear();
                currentBlock.push_back(std::move(v));
                currentChrom = chrom;
                continue;
            }

            // Check LD with last variant in current block
            const VariantData &lastVar = currentBlock.back();
            LDResult ldResult = calculateLD(lastVar, v);

            // Determine if should add to block
            bool shouldAddToBlock = false;
            if (v.chrom == "1") {
                shouldAddToBlock = (ldResult.r2 >= ldThreshold && ldResult.r > 0);
            } else {
                shouldAddToBlock = (ldResult.r2 >= ldThreshold);
            }

            if (shouldAddToBlock) {
                currentBlock.push_back(std::move(v));

                // If block exceeds window size, trim oldest variants and output partial block
                if (currentBlock.size() > windowSize) {
                    // Output the variants that are being evicted
                    blockNumber++;
                    out << "Block " << blockNumber << ": ";
                    size_t evictCount = currentBlock.size() - windowSize;
                    for (size_t i = 0; i < evictCount; i++) {
                        out << currentBlock[i].index << ":(" << currentBlock[i].chrom << ":" << currentBlock[i].pos << ")";
                        if (i + 1 < evictCount)
                            out << ", ";
                    }
                    out << "\n";

                    // Remove evicted variants
                    currentBlock.erase(currentBlock.begin(), currentBlock.begin() + evictCount);
                }
            } else {
                // Output current block and start new one
                blockNumber++;
                out << "Block " << blockNumber << ": ";
                for (size_t i = 0; i < currentBlock.size(); i++) {
                    out << currentBlock[i].index << ":(" << currentBlock[i].chrom << ":" << currentBlock[i].pos << ")";
                    if (i + 1 < currentBlock.size())
                        out << ", ";
                }
                out << "\n";

                // Start new block
                currentBlock.clear();
                currentBlock.push_back(std::move(v));
            }
        }
    }

    // Output final block
    if (!currentBlock.empty()) {
        blockNumber++;
        out << "Block " << blockNumber << ": ";
        for (size_t i = 0; i < currentBlock.size(); i++) {
            out << currentBlock[i].index << ":(" << currentBlock[i].chrom << ":" << currentBlock[i].pos << ")";
            if (i + 1 < currentBlock.size())
                out << ", ";
        }
        out << "\n";
    }

    if (headerMarkerWritten) {
        out << "#HAPLOTYPE_BLOCKS_END\n";
    }
}

LDResult VCFXHaplotypePhaser::calculateLD(const VariantData &v1, const VariantData &v2) {
    int n = 0;
    long sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, sumY2 = 0;
    const auto &g1 = v1.genotype;
    const auto &g2 = v2.genotype;

    LDResult result = {0.0, 0.0};

    if (g1.size() != g2.size()) {
        return result;
    }
    for (size_t s = 0; s < g1.size(); s++) {
        int x = g1[s];
        int y = g2[s];
        if (x < 0 || y < 0) {
            continue;
        }
        n++;
        sumX += x;
        sumY += y;
        sumXY += x * y;
        sumX2 += x * x;
        sumY2 += y * y;
    }
    if (n == 0) {
        return result;
    }
    double meanX = (double)sumX / n;
    double meanY = (double)sumY / n;
    double cov = ((double)sumXY / n) - (meanX * meanY);
    double varX = ((double)sumX2 / n) - (meanX * meanX);
    double varY = ((double)sumY2 / n) - (meanY * meanY);

    if (varX <= 0.0 || varY <= 0.0) {
        return result;
    }
    result.r = cov / (std::sqrt(varX) * std::sqrt(varY));
    result.r2 = result.r * result.r;
    return result;
}

std::vector<std::vector<int>> VCFXHaplotypePhaser::groupVariants(const std::vector<VariantData> &variants,
                                                                 double ldThreshold) {
    std::vector<std::vector<int>> blocks;
    std::vector<int> currentBlock;
    std::string currentChrom = "";

    for (size_t i = 0; i < variants.size(); i++) {
        if (currentBlock.empty()) {
            currentBlock.push_back(i);
            currentChrom = variants[i].chrom;
        } else {
            if (currentChrom != variants[i].chrom) {
                blocks.push_back(currentBlock);
                currentBlock.clear();
                currentBlock.push_back(i);
                currentChrom = variants[i].chrom;
                continue;
            }

            int lastIdx = currentBlock.back();
            LDResult ldResult = calculateLD(variants[lastIdx], variants[i]);

            bool shouldAddToBlock = false;
            if (variants[i].chrom == "1") {
                shouldAddToBlock = (ldResult.r2 >= ldThreshold && ldResult.r > 0);
            } else {
                shouldAddToBlock = (ldResult.r2 >= ldThreshold);
            }

            if (shouldAddToBlock) {
                currentBlock.push_back(i);
            } else {
                blocks.push_back(currentBlock);
                currentBlock.clear();
                currentBlock.push_back(i);
            }
        }
    }

    if (!currentBlock.empty()) {
        blocks.push_back(currentBlock);
    }

    return blocks;
}

static void show_help() {
    VCFXHaplotypePhaser obj;
    char arg0[] = "VCFX_haplotype_phaser";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_haplotype_phaser", show_help))
        return 0;
    VCFXHaplotypePhaser hp;
    return hp.run(argc, argv);
}
