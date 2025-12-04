#include "VCFX_haplotype_extractor.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>

// ---------------------------------------------------------------------
// printHelp
// ---------------------------------------------------------------------
void printHelp() {
    std::cout << "VCFX_haplotype_extractor\n"
              << "Usage: VCFX_haplotype_extractor [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h                 Display this help message and exit.\n"
              << "  --block-size <int>         Maximum distance for grouping consecutive variants (default 100000).\n"
              << "  --check-phase-consistency  If set, try a minimal check across variants.\n"
              << "  --streaming                Enable streaming mode: output blocks immediately when complete.\n"
              << "                             Uses O(block_size) memory instead of O(total_variants).\n"
              << "  --debug                    Output verbose debug information.\n\n"
              << "Description:\n"
              << "  Extracts phased haplotype blocks from genotype data in a VCF file. "
              << "It reconstructs haplotypes for each sample by analyzing phased genotype fields.\n\n"
              << "Performance:\n"
              << "  Default mode:   Accumulates all blocks in memory, outputs at end.\n"
              << "  Streaming mode: Outputs blocks immediately when complete. Enables\n"
              << "                  processing of arbitrarily large files with bounded memory.\n\n"
              << "Examples:\n"
              << "  ./VCFX_haplotype_extractor --block-size 50000 < phased.vcf > haplotypes.tsv\n"
              << "  ./VCFX_haplotype_extractor --streaming < large_phased.vcf > haplotypes.tsv\n"
              << "  ./VCFX_haplotype_extractor --streaming --block-size 10000 < phased.vcf > haplotypes.tsv\n";
}

// ---------------------------------------------------------------------
// parseHeader: parse #CHROM line => sample columns
// ---------------------------------------------------------------------
bool HaplotypeExtractor::parseHeader(const std::string &headerLine) {
    auto fields = splitString(headerLine, '\t');
    if (fields.size() <= 8) {
        std::cerr << "Error: VCF header does not contain sample columns.\n";
        return false;
    }
    // from col=9 onward => sample names
    sampleNames.assign(fields.begin() + 9, fields.end());
    numSamples = sampleNames.size();
    return true;
}

// ---------------------------------------------------------------------
// splitString
// ---------------------------------------------------------------------
std::vector<std::string> HaplotypeExtractor::splitString(const std::string &str, char delimiter) {
    std::vector<std::string> tokens;
    tokens.reserve(16);
    std::stringstream ss(str);
    std::string tmp;
    while (std::getline(ss, tmp, delimiter)) {
        tokens.push_back(tmp);
    }
    return tokens;
}

// ---------------------------------------------------------------------
// areAllSamplesPhased: checks each genotype for '|'
// ---------------------------------------------------------------------
bool HaplotypeExtractor::areAllSamplesPhased(const std::vector<std::string> &genotypes) {
    for (const auto &g : genotypes) {
        if (g.find('|') == std::string::npos) {
            return false;
        }
    }
    return true;
}

// ---------------------------------------------------------------------
// phaseIsConsistent: minimal check across variants
// ---------------------------------------------------------------------
bool HaplotypeExtractor::phaseIsConsistent(const HaplotypeBlock &block, const std::vector<std::string> &newGenotypes) {
    if (block.haplotypes.size() != newGenotypes.size()) {
        return false;
    }

    if (debugMode) {
        std::cerr << "Checking phase consistency\n";
    }

    for (size_t s = 0; s < block.haplotypes.size(); s++) {
        const std::string &allVar = block.haplotypes[s];

        std::string lastGT;
        if (allVar.size() < 3) {
            lastGT = allVar;
        } else {
            size_t lastPipePos = allVar.rfind('|', allVar.size() - 2);
            if (lastPipePos == std::string::npos) {
                lastGT = allVar;
            } else {
                lastGT = allVar.substr(lastPipePos - 1);
            }
        }

        if (debugMode) {
            std::cerr << "Sample " << s << " last GT: " << lastGT << " new GT: " << newGenotypes[s] << "\n";
        }

        if (lastGT.size() < 3 || newGenotypes[s].size() < 3) {
            continue;
        }

        char lastAllele1 = (lastGT.size() >= 1) ? lastGT[lastGT.size() - 3] : '.';
        char lastAllele2 = (lastGT.size() >= 1) ? lastGT[lastGT.size() - 1] : '.';
        char newAllele1 = newGenotypes[s][0];
        char newAllele2 = newGenotypes[s][2];

        if (debugMode) {
            std::cerr << "Comparing alleles: " << lastAllele1 << "|" << lastAllele2 << " vs " << newAllele1 << "|"
                      << newAllele2 << "\n";
        }

        if (lastAllele1 != newAllele1 && lastAllele2 != newAllele2 && lastAllele1 == newAllele2 &&
            lastAllele2 == newAllele1) {
            if (debugMode) {
                std::cerr << "Phase flip detected in sample " << s << "\n";
            }
            return false;
        }
    }

    if (debugMode) {
        std::cerr << "All phases consistent\n";
    }
    return true;
}

// ---------------------------------------------------------------------
// updateBlocks: merges new variant into last block or starts new block
// ---------------------------------------------------------------------
void HaplotypeExtractor::updateBlocks(std::vector<HaplotypeBlock> &haplotypeBlocks, const std::string &chrom, int pos,
                                      const std::vector<std::string> &genotypes) {
    if (haplotypeBlocks.empty()) {
        HaplotypeBlock b;
        b.chrom = chrom;
        b.start = pos;
        b.end = pos;
        b.haplotypes = genotypes;
        haplotypeBlocks.push_back(std::move(b));
        return;
    }

    HaplotypeBlock &lastB = haplotypeBlocks.back();
    bool canExtend = (chrom == lastB.chrom) && (pos - lastB.end <= blockDistanceThreshold);

    if (canExtend && checkPhaseConsistency) {
        canExtend = phaseIsConsistent(lastB, genotypes);
    }

    if (canExtend) {
        lastB.end = pos;
        for (size_t s = 0; s < lastB.haplotypes.size(); s++) {
            lastB.haplotypes[s] += "|" + genotypes[s];
        }
    } else {
        HaplotypeBlock nb;
        nb.chrom = chrom;
        nb.start = pos;
        nb.end = pos;
        nb.haplotypes = genotypes;
        haplotypeBlocks.push_back(std::move(nb));
    }
}

// ---------------------------------------------------------------------
// outputHeader: write header line
// ---------------------------------------------------------------------
void HaplotypeExtractor::outputHeader(std::ostream &out) {
    out << "CHROM\tSTART\tEND";
    for (const auto &s : sampleNames) {
        out << "\t" << s;
    }
    out << "\n";
}

// ---------------------------------------------------------------------
// outputBlock: write a single block to output
// ---------------------------------------------------------------------
void HaplotypeExtractor::outputBlock(std::ostream &out, const HaplotypeBlock &block) {
    out << block.chrom << "\t" << block.start << "\t" << block.end;
    for (const auto &hap : block.haplotypes) {
        out << "\t" << hap;
    }
    out << "\n";
}

// ---------------------------------------------------------------------
// processVariant: parse data line, check phased, update blocks
// ---------------------------------------------------------------------
bool HaplotypeExtractor::processVariant(const std::vector<std::string> &fields,
                                        std::vector<HaplotypeBlock> &haplotypeBlocks) {
    if (fields.size() < 9) {
        std::cerr << "Warning: skipping invalid VCF line (<9 fields)\n";
        return false;
    }
    const std::string &chrom = fields[0];
    int pos = 0;
    try {
        pos = std::stoi(fields[1]);
    } catch (...) {
        std::cerr << "Warning: invalid POS => skip variant\n";
        return false;
    }

    const std::string &formatStr = fields[8];
    auto formatToks = splitString(formatStr, ':');
    int gtIndex = -1;
    for (int i = 0; i < (int)formatToks.size(); i++) {
        if (formatToks[i] == "GT") {
            gtIndex = i;
            break;
        }
    }
    if (gtIndex < 0) {
        return false;
    }

    std::vector<std::string> genotypeFields(numSamples);
    bool allPhased = true;
    for (size_t s = 0; s < numSamples; s++) {
        if (9 + s >= fields.size()) {
            genotypeFields[s] = ".|.";
            allPhased = false;
            continue;
        }
        auto sampleToks = splitString(fields[9 + s], ':');
        if (gtIndex >= (int)sampleToks.size()) {
            genotypeFields[s] = ".|.";
            allPhased = false;
            continue;
        }
        std::string gt = sampleToks[gtIndex];
        if (gt.find('|') == std::string::npos) {
            allPhased = false;
        }
        genotypeFields[s] = gt;
    }

    if (!allPhased) {
        std::cerr << "Warning: Not all samples phased at " << chrom << ":" << pos << ".\n";
        return false;
    }

    updateBlocks(haplotypeBlocks, chrom, pos, genotypeFields);
    return true;
}

// ---------------------------------------------------------------------
// processVariantStreaming: streaming version that outputs completed blocks
// ---------------------------------------------------------------------
bool HaplotypeExtractor::processVariantStreaming(const std::vector<std::string> &fields,
                                                  HaplotypeBlock &currentBlock,
                                                  bool &hasCurrentBlock,
                                                  std::ostream &out,
                                                  bool headerWritten) {
    if (fields.size() < 9) {
        std::cerr << "Warning: skipping invalid VCF line (<9 fields)\n";
        return false;
    }

    const std::string &chrom = fields[0];
    int pos = 0;
    try {
        pos = std::stoi(fields[1]);
    } catch (...) {
        std::cerr << "Warning: invalid POS => skip variant\n";
        return false;
    }

    const std::string &formatStr = fields[8];
    auto formatToks = splitString(formatStr, ':');
    int gtIndex = -1;
    for (int i = 0; i < (int)formatToks.size(); i++) {
        if (formatToks[i] == "GT") {
            gtIndex = i;
            break;
        }
    }
    if (gtIndex < 0) {
        return false;
    }

    std::vector<std::string> genotypeFields(numSamples);
    bool allPhased = true;
    for (size_t s = 0; s < numSamples; s++) {
        if (9 + s >= fields.size()) {
            genotypeFields[s] = ".|.";
            allPhased = false;
            continue;
        }
        auto sampleToks = splitString(fields[9 + s], ':');
        if (gtIndex >= (int)sampleToks.size()) {
            genotypeFields[s] = ".|.";
            allPhased = false;
            continue;
        }
        std::string gt = sampleToks[gtIndex];
        if (gt.find('|') == std::string::npos) {
            allPhased = false;
        }
        genotypeFields[s] = gt;
    }

    if (!allPhased) {
        std::cerr << "Warning: Not all samples phased at " << chrom << ":" << pos << ".\n";
        return false;
    }

    // Streaming logic: check if we need to output current block and start new one
    if (!hasCurrentBlock) {
        // Start first block
        currentBlock.chrom = chrom;
        currentBlock.start = pos;
        currentBlock.end = pos;
        currentBlock.haplotypes = genotypeFields;
        hasCurrentBlock = true;
    } else {
        // Check if we can extend current block
        bool canExtend = (chrom == currentBlock.chrom) && (pos - currentBlock.end <= blockDistanceThreshold);

        if (canExtend && checkPhaseConsistency) {
            canExtend = phaseIsConsistent(currentBlock, genotypeFields);
        }

        if (canExtend) {
            // Extend current block
            currentBlock.end = pos;
            for (size_t s = 0; s < currentBlock.haplotypes.size(); s++) {
                currentBlock.haplotypes[s] += "|" + genotypeFields[s];
            }
        } else {
            // Output completed block and start new one
            outputBlock(out, currentBlock);

            // Start new block
            currentBlock.chrom = chrom;
            currentBlock.start = pos;
            currentBlock.end = pos;
            currentBlock.haplotypes = genotypeFields;
        }
    }

    return true;
}

// ---------------------------------------------------------------------
// extractHaplotypes: original mode - accumulates all blocks
// ---------------------------------------------------------------------
bool HaplotypeExtractor::extractHaplotypes(std::istream &in, std::ostream &out) {
    bool foundHeader = false;
    std::vector<HaplotypeBlock> haplotypeBlocks;
    std::vector<std::string> fields;
    fields.reserve(16);

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty())
            continue;
        if (line[0] == '#') {
            if (!foundHeader && line.rfind("#CHROM", 0) == 0) {
                foundHeader = true;
                if (!parseHeader(line)) {
                    return false;
                }
            }
            continue;
        }
        if (!foundHeader) {
            std::cerr << "Error: no #CHROM header found before data.\n";
            return false;
        }
        vcfx::split_tabs(line, fields);
        processVariant(fields, haplotypeBlocks);
    }

    // Output header and all blocks at end
    outputHeader(out);
    for (const auto &block : haplotypeBlocks) {
        outputBlock(out, block);
    }

    return true;
}

// ---------------------------------------------------------------------
// extractHaplotypesStreaming: streaming mode - outputs blocks immediately
// ---------------------------------------------------------------------
bool HaplotypeExtractor::extractHaplotypesStreaming(std::istream &in, std::ostream &out) {
    bool foundHeader = false;
    bool headerWritten = false;
    HaplotypeBlock currentBlock;
    bool hasCurrentBlock = false;
    std::vector<std::string> fields;
    fields.reserve(16);

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty())
            continue;
        if (line[0] == '#') {
            if (!foundHeader && line.rfind("#CHROM", 0) == 0) {
                foundHeader = true;
                if (!parseHeader(line)) {
                    return false;
                }
                // Output header immediately in streaming mode
                outputHeader(out);
                headerWritten = true;
            }
            continue;
        }
        if (!foundHeader) {
            std::cerr << "Error: no #CHROM header found before data.\n";
            return false;
        }

        vcfx::split_tabs(line, fields);
        processVariantStreaming(fields, currentBlock, hasCurrentBlock, out, headerWritten);
    }

    // Output final block if any
    if (hasCurrentBlock) {
        outputBlock(out, currentBlock);
    }

    return true;
}

// ---------------------------------------------------------------------
// main
// ---------------------------------------------------------------------
static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_haplotype_extractor", show_help))
        return 0;

    int blockSize = 100000;
    bool doCheck = false;
    bool debug = false;
    bool streaming = false;

    // Simple arg parse
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--help" || a == "-h") {
            printHelp();
            return 0;
        } else if (a == "--block-size" && i + 1 < argc) {
            blockSize = std::stoi(argv[++i]);
        } else if (a == "--check-phase-consistency") {
            doCheck = true;
        } else if (a == "--debug") {
            debug = true;
        } else if (a == "--streaming") {
            streaming = true;
        }
    }

    HaplotypeExtractor extractor;
    extractor.setBlockDistanceThreshold(blockSize);
    extractor.setCheckPhaseConsistency(doCheck);
    extractor.setDebug(debug);
    extractor.setStreamingMode(streaming);

    bool ok;
    if (streaming) {
        ok = extractor.extractHaplotypesStreaming(std::cin, std::cout);
    } else {
        ok = extractor.extractHaplotypes(std::cin, std::cout);
    }

    return (ok ? 0 : 1);
}
