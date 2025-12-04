#include "VCFX_ref_comparator.h"
#include "vcfx_core.h" 
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <vector>

// A small helper to remove whitespace
static inline void stripSpaces(std::string &s) {
    s.erase(std::remove_if(s.begin(), s.end(), [](unsigned char c) { return std::isspace(c); }), s.end());
}

int VCFXRefComparator::run(int argc, char *argv[]) {
    if (argc == 1) {
        displayHelp();
        return 0;
    }
    bool showHelp = false;
    std::string referencePath;
    static struct option long_opts[] = {
        {"help", no_argument, 0, 'h'}, {"reference", required_argument, 0, 'r'}, {0, 0, 0, 0}};
    while (true) {
        int c = ::getopt_long(argc, argv, "hr:", long_opts, nullptr);
        if (c == -1)
            break;
        switch (c) {
        case 'h':
            showHelp = true;
            break;
        case 'r':
            referencePath = optarg;
            break;
        default:
            showHelp = true;
        }
    }
    if (showHelp) {
        displayHelp();
        return 0;
    }
    if (referencePath.empty()) {
        std::cerr << "Error: must specify --reference <FASTA>.\n";
        displayHelp();
        return 1;
    }
    if (!loadReference(referencePath)) {
        std::cerr << "Error: failed to load reference from " << referencePath << "\n";
        return 1;
    }
    compareVCF(std::cin, std::cout);
    return 0;
}

void VCFXRefComparator::displayHelp() {
    std::cout << "VCFX_ref_comparator: Compare VCF REF/ALT with a reference genome.\n\n"
                 "Usage:\n"
                 "  VCFX_ref_comparator --reference ref.fasta < input.vcf > output.vcf\n\n"
                 "Description:\n"
                 "  Reads a reference FASTA into memory. Then reads each variant line:\n"
                 "   - If chromosome or position is invalid, logs a warning and sets REF_COMPARISON=UNKNOWN_CHROM or "
                 "INVALID_POS.\n"
                 "   - Otherwise, compares the VCF's REF vs the reference substring. Then for each ALT, indicates "
                 "'REF_MATCH' if ALT= reference substring or 'NOVEL'.\n"
                 "  The result is appended to the 'INFO' field as REF_COMPARISON=...\n\n"
                 "Example:\n"
                 "  VCFX_ref_comparator --reference genome.fa < in.vcf > out.vcf\n";
}

bool VCFXRefComparator::loadReference(const std::string &referenceFastaPath) {
    std::ifstream in(referenceFastaPath);
    if (!in.is_open()) {
        std::cerr << "Error: cannot open reference " << referenceFastaPath << "\n";
        return false;
    }
    referenceGenome.clear();
    std::string line, currentChrom;
    std::ostringstream seq;
    while (true) {
        if (!std::getline(in, line))
            break;
        if (line.empty())
            continue;
        if (line[0] == '>') {
            // store old chrom
            if (!currentChrom.empty()) {
                referenceGenome[currentChrom] = seq.str();
            }
            seq.str("");
            seq.clear();
            // parse new chrom
            currentChrom = line.substr(1);
            // if there's whitespace, strip it after first token
            {
                std::istringstream iss(currentChrom);
                iss >> currentChrom;
            }
            // uppercase
            std::transform(currentChrom.begin(), currentChrom.end(), currentChrom.begin(), ::toupper);
        } else {
            stripSpaces(line);
            std::transform(line.begin(), line.end(), line.begin(), ::toupper);
            seq << line;
        }
    }
    if (!currentChrom.empty()) {
        referenceGenome[currentChrom] = seq.str();
    }
    return true;
}

void VCFXRefComparator::compareVCF(std::istream &vcfIn, std::ostream &vcfOut) {
    bool foundChromHeader = false;
    infoHeaderInserted = false;
    std::string line;
    while (true) {
        if (!std::getline(vcfIn, line))
            break;
        if (line.empty()) {
            vcfOut << line << "\n";
            continue;
        }
        if (line[0] == '#') {
            // check if #CHROM
            if (line.rfind("#CHROM", 0) == 0) {
                foundChromHeader = true;
                // Insert new INFO line BEFORE we write #CHROM
                if (!infoHeaderInserted) {
                    vcfOut << "##INFO=<ID=REF_COMPARISON,Number=1,Type=String,Description=\"Comparison of REF/ALT vs "
                              "reference genome substring\">\n";
                    infoHeaderInserted = true;
                }
                vcfOut << line << "\n";
            } else {
                vcfOut << line << "\n";
            }
            continue;
        }
        if (!foundChromHeader) {
            std::cerr << "Warning: data line encountered before #CHROM => skipping.\n";
            continue;
        }
        // parse fields
        std::stringstream ss(line);
        std::vector<std::string> fields;
        {
            std::string col;
            while (std::getline(ss, col, '\t')) {
                fields.push_back(col);
            }
        }
        if (fields.size() < 8) {
            std::cerr << "Warning: VCF line has <8 columns => skipping.\n";
            continue;
        }
        // fields: 0=CHROM,1=POS,2=ID,3=REF,4=ALT,5=QUAL,6=FILTER,7=INFO,...
        std::string &chrom = fields[0];
        std::string &posStr = fields[1];
        std::string &ref = fields[3];
        std::string &alt = fields[4];
        std::string &info = fields[7];

        // uppercase chrom
        std::transform(chrom.begin(), chrom.end(), chrom.begin(), ::toupper);

        int pos = 0;
        try {
            pos = std::stoi(posStr);
        } catch (...) {
            // out of range
            if (!info.empty() && info.back() != ';')
                info += ';';
            info += "REF_COMPARISON=INVALID_POS";
            // rewrite
            std::ostringstream newLine;
            for (int i = 0; i < 7; i++) {
                newLine << fields[i];
                if (i < 7)
                    newLine << "\t";
            }
            newLine << info;
            for (size_t i = 8; i < fields.size(); i++) {
                newLine << "\t" << fields[i];
            }
            vcfOut << newLine.str() << "\n";
            continue;
        }
        // find reference
        auto it = referenceGenome.find(chrom);
        if (it == referenceGenome.end()) {
            // unknown chrom
            if (!info.empty() && info.back() != ';')
                info += ';';
            info += "REF_COMPARISON=UNKNOWN_CHROM";
            // rewrite
            std::ostringstream newLine;
            for (int i = 0; i < 7; i++) {
                newLine << fields[i];
                if (i < 7)
                    newLine << "\t";
            }
            newLine << info;
            for (size_t i = 8; i < fields.size(); i++) {
                newLine << "\t" << fields[i];
            }
            vcfOut << newLine.str() << "\n";
            continue;
        }
        const std::string &seq = it->second;
        if (pos < 1 || pos > (int)seq.size()) {
            // invalid pos
            if (!info.empty() && info.back() != ';')
                info += ';';
            info += "REF_COMPARISON=INVALID_POS";
            // rewrite
            std::ostringstream newLine;
            for (int i = 0; i < 7; i++) {
                newLine << fields[i];
                if (i < 7)
                    newLine << "\t";
            }
            newLine << info;
            for (size_t i = 8; i < fields.size(); i++) {
                newLine << "\t" << fields[i];
            }
            vcfOut << newLine.str() << "\n";
            continue;
        }
        // ref from genome
        std::string genomeRef = seq.substr(pos - 1, ref.size()); // 1-based
        // uppercase
        std::transform(ref.begin(), ref.end(), ref.begin(), ::toupper);

        // compare ref vs genomeRef
        bool refMatch = (ref == genomeRef);

        // split alt by comma
        std::vector<std::string> altAlleles;
        {
            std::stringstream as(alt);
            std::string a;
            while (std::getline(as, a, ',')) {
                std::transform(a.begin(), a.end(), a.begin(), ::toupper);
                altAlleles.push_back(a);
            }
        }
        // for each alt, see if alt== genomeRef => "REF_MATCH" else "NOVEL"
        std::vector<std::string> comparisons;
        for (auto &a : altAlleles) {
            if (a == genomeRef)
                comparisons.push_back("REF_MATCH");
            else
                comparisons.push_back("NOVEL");
        }
        // build comparisonStr
        std::string comparisonStr;
        for (size_t i = 0; i < comparisons.size(); i++) {
            if (i > 0)
                comparisonStr += ",";
            comparisonStr += comparisons[i];
        }
        if (!info.empty() && info.back() != ';')
            info += ';';
        // e.g. "REF_COMPARISON=REF_MATCH,REF_MATCH" or "NOVEL" etc.
        if (!refMatch) {
            // optionally label mismatch => we won't specifically do that
            // the alt status is in the comparisons
        }
        info += "REF_COMPARISON=" + comparisonStr;

        // rebuild line
        std::ostringstream outLine;
        for (int i = 0; i < 7; i++) {
            outLine << fields[i];
            if (i < 6)
                outLine << "\t"; // Changed from i<7 to i<6 to avoid adding an extra tab after FILTER
        }
        outLine << "\t" << info; // Add tab after FILTER, then add INFO without spaces
        // any other columns
        for (size_t i = 8; i < fields.size(); i++) {
            outLine << "\t" << fields[i];
        }
        vcfOut << outLine.str() << "\n";
    }
}

static void show_help() {
    VCFXRefComparator obj;
    char arg0[] = "VCFX_ref_comparator";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_ref_comparator", show_help))
        return 0;
    VCFXRefComparator refComp;
    return refComp.run(argc, argv);
}
