#include "VCFX_alignment_checker.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

int VCFXAlignmentChecker::run(int argc, char *argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'}, {"alignment-discrepancy", no_argument, 0, 'a'}, {0, 0, 0, 0}};

    while ((opt = getopt_long(argc, argv, "ha", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
            break;
        case 'a':
            // Alignment-discrepancy mode (no extra actions needed here)
            break;
        default:
            showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 1;
    }

    // We expect two positional arguments after the options:
    //   1) VCF file  2) Reference FASTA file
    if ((argc - optind) < 2) {
        std::cerr << "Error: Missing required arguments: <vcf_file> <reference.fasta>\n";
        displayHelp();
        return 1;
    }

    // Get file paths
    std::string vcfFile = argv[optind];
    std::string refFile = argv[optind + 1];

    // Open VCF input
    std::ifstream vcfStream(vcfFile);
    if (!vcfStream.is_open()) {
        std::cerr << "Error: Unable to open VCF file: " << vcfFile << "\n";
        return 1;
    }

    // Load reference genome index
    if (!loadReferenceGenome(refFile)) {
        std::cerr << "Error: Failed to load reference genome.\n";
        return 1;
    }

    // Check discrepancies (results to stdout)
    checkDiscrepancies(vcfStream, std::cout);

    return 0;
}

void VCFXAlignmentChecker::displayHelp() {
    std::cout << "VCFX_alignment_checker: Identify discrepancies between VCF variants and a reference genome.\n\n"
              << "Usage:\n"
              << "  VCFX_alignment_checker --alignment-discrepancy <vcf_file> <reference.fasta>\n\n"
              << "Options:\n"
              << "  -h, --help                   Display this help message and exit\n"
              << "  -a, --alignment-discrepancy  Identify alignment discrepancies\n\n"
              << "Example:\n"
              << "  VCFX_alignment_checker --alignment-discrepancy input.vcf reference.fasta > discrepancies.txt\n";
}

bool VCFXAlignmentChecker::loadReferenceGenome(const std::string &path) {
    referencePath = path;
    referenceIndex.clear();

    referenceStream.open(path, std::ios::in);
    if (!referenceStream.is_open()) {
        std::cerr << "Error: Unable to open reference genome file: " << path << "\n";
        return false;
    }

    std::string line;
    std::string currentChrom;
    FastaIndexEntry entry;
    std::size_t seqLen = 0;

    // record file offset where we will read sequence lines
    while (std::getline(referenceStream, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '>') {
            if (!currentChrom.empty()) {
                entry.length = seqLen;
                referenceIndex[normalizeChromosome(currentChrom)] = entry;
            }

            currentChrom.clear();
            seqLen = 0;
            entry = FastaIndexEntry();

            std::size_t pos = line.find(' ');
            if (pos != std::string::npos) {
                currentChrom = line.substr(1, pos - 1);
            } else {
                currentChrom = line.substr(1);
            }

            entry.offset = referenceStream.tellg();
            entry.basesPerLine = 0;
            entry.bytesPerLine = 0;
        } else {
            if (entry.basesPerLine == 0) {
                entry.basesPerLine = line.size();
                entry.bytesPerLine = line.size() + 1; // assume single '\n'
            }
            seqLen += line.size();
        }
    }

    if (!currentChrom.empty()) {
        entry.length = seqLen;
        referenceIndex[normalizeChromosome(currentChrom)] = entry;
    }

    referenceStream.clear();
    referenceStream.seekg(0);
    return true;
}

std::string VCFXAlignmentChecker::normalizeChromosome(const std::string &chrom) {
    std::string norm = chrom;
    // convert to upper and drop leading "CHR" if present
    if (norm.size() >= 3 && (norm.rfind("chr", 0) == 0 || norm.rfind("CHR", 0) == 0)) {
        norm = norm.substr(3);
    }
    std::transform(norm.begin(), norm.end(), norm.begin(), ::toupper);
    return norm;
}

std::string VCFXAlignmentChecker::getReferenceBases(const std::string &chrom, int pos, int length) {
    auto it = referenceIndex.find(normalizeChromosome(chrom));
    if (it == referenceIndex.end()) {
        return "";
    }

    const FastaIndexEntry &entry = it->second;
    if (pos < 1 || static_cast<std::size_t>(pos - 1) >= entry.length) {
        return "";
    }

    int remaining = length;
    std::size_t currPos = static_cast<std::size_t>(pos - 1);
    std::string result;
    result.reserve(length);

    while (remaining > 0 && currPos < entry.length) {
        std::size_t lineIdx = currPos / entry.basesPerLine;
        std::size_t lineOffset = currPos % entry.basesPerLine;
        std::size_t chunk = std::min<std::size_t>(entry.basesPerLine - lineOffset, remaining);

        std::streampos filePos = entry.offset + static_cast<std::streampos>(lineIdx * entry.bytesPerLine + lineOffset);
        referenceStream.clear();
        referenceStream.seekg(filePos);
        std::string buf(chunk, '\0');
        referenceStream.read(&buf[0], chunk);
        result += buf;

        currPos += chunk;
        remaining -= static_cast<int>(chunk);
    }

    return result;
}

void VCFXAlignmentChecker::checkDiscrepancies(std::istream &vcfIn, std::ostream &out) {
    // We do not need an extra reference stream parameter; the reference is already loaded in memory.

    std::string line;
    bool headerParsed = false;
    int chrIndex = -1, posIndex = -1, refIndex = -1, altIndex = -1;

    // Reusable buffer for parsing
    std::vector<std::string> fields;
    fields.reserve(16);

    // Print a header for the discrepancies table
    out << "CHROM\tPOS\tID\tREF\tALT\tDiscrepancy_Type\tReference_Value\tVCF_Value\n";

    while (std::getline(vcfIn, line)) {
        if (line.empty()) {
            continue;
        }

        // Header lines
        if (line[0] == '#') {
            // If it's the #CHROM line, parse the column indices
            if (line.rfind("#CHROM", 0) == 0) {
                std::vector<std::string> headers;
                vcfx::split_tabs(line, headers);
                if (!headers.empty() && !headers[0].empty() && headers[0][0] == '#') {
                    headers[0].erase(0, 1); // drop leading '#'
                }
                for (size_t i = 0; i < headers.size(); ++i) {
                    if (headers[i] == "CHROM")
                        chrIndex = static_cast<int>(i);
                    else if (headers[i] == "POS")
                        posIndex = static_cast<int>(i);
                    else if (headers[i] == "REF")
                        refIndex = static_cast<int>(i);
                    else if (headers[i] == "ALT")
                        altIndex = static_cast<int>(i);
                }
                if (chrIndex == -1 || posIndex == -1 || refIndex == -1 || altIndex == -1) {
                    std::cerr << "Error: VCF header does not contain required CHROM, POS, REF, ALT fields.\n";
                    return;
                }
                headerParsed = true;
            }
            continue;
        }

        // We should have a valid header by now
        if (!headerParsed) {
            std::cerr << "Error: VCF #CHROM header line not found before data lines.\n";
            return;
        }

        // Split line into fields
        fields.clear();
        vcfx::split_tabs(line, fields);

        // Basic sanity check
        if (fields.size() < static_cast<size_t>(altIndex + 1)) {
            std::cerr << "Warning: Skipping invalid VCF line (insufficient fields): " << line << "\n";
            continue;
        }

        std::string chrom = fields[chrIndex];
        int posVal = 0;
        try {
            posVal = std::stoi(fields[posIndex]);
        } catch (...) {
            std::cerr << "Warning: Invalid POS value. Skipping line: " << line << "\n";
            continue;
        }

        std::string ref = fields[refIndex];
        std::string alt = fields[altIndex];
        std::string id = (fields.size() > 2) ? fields[2] : ".";

        // Handle multi-allelic ALTs
        std::vector<std::string> alts;
        {
            std::stringstream altSS(alt);
            std::string altAllele;
            while (std::getline(altSS, altAllele, ',')) {
                alts.push_back(altAllele);
            }
        }

        // Check each ALT allele
        for (const auto &allele : alts) {
            // If both REF and ALT are single bases, treat as SNP
            if (ref.size() == 1 && allele.size() == 1) {
                std::string ref_base = getReferenceBases(chrom, posVal, 1);
                if (ref_base.empty()) {
                    std::cerr << "Warning: Reference base not found for " << chrom << ":" << posVal << "\n";
                    continue;
                }
                // Compare REF in VCF vs reference genome
                if (ref != ref_base) {
                    out << chrom << "\t" << posVal << "\t" << id << "\t" << ref << "\t" << allele << "\t"
                        << "REF_MISMATCH"
                        << "\t" << ref_base << "\t" << ref << "\n";
                }
                // Compare ALT to the reference base at the same position.
                // Here we flag a mismatch when the ALT allele is actually the
                // same as the reference (i.e. not a true variant).
                std::string alt_base = ref_base; // The reference at that position
                if (allele == alt_base) {
                    out << chrom << "\t" << posVal << "\t" << id << "\t" << ref << "\t" << allele << "\t"
                        << "ALT_MISMATCH"
                        << "\t" << alt_base << "\t" << allele << "\n";
                }
            } else {
                // Indel or complex
                size_t ref_len = ref.size();
                size_t alt_len = allele.size();
                // Compare as many bases as the shorter string has
                size_t len = std::min(ref_len, alt_len);

                std::string ref_seq = getReferenceBases(chrom, posVal, static_cast<int>(len));
                if (ref_seq.empty()) {
                    std::cerr << "Warning: Reference sequence not found for " << chrom << ":" << posVal << "\n";
                    continue;
                }

                std::string vcf_ref = ref.substr(0, len);
                std::string vcf_alt = allele.substr(0, len);

                if (vcf_ref != ref_seq) {
                    out << chrom << "\t" << posVal << "\t" << id << "\t" << ref << "\t" << allele << "\t"
                        << "REF_DISCREPANCY"
                        << "\t" << ref_seq << "\t" << vcf_ref << "\n";
                }
                if (vcf_alt != ref_seq) {
                    out << chrom << "\t" << posVal << "\t" << id << "\t" << ref << "\t" << allele << "\t"
                        << "ALT_DISCREPANCY"
                        << "\t" << ref_seq << "\t" << vcf_alt << "\n";
                }
            }
        }
    }
}

// Typical main(), linking to run()
static void show_help() {
    VCFXAlignmentChecker obj;
    char arg0[] = "VCFX_alignment_checker";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_alignment_checker", show_help))
        return 0;
    VCFXAlignmentChecker alignmentChecker;
    return alignmentChecker.run(argc, argv);
}
