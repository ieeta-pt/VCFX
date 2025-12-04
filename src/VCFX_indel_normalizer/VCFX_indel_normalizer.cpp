#include "VCFX_indel_normalizer.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------
// Print usage
// ---------------------------------------------------------------------
void VCFXIndelNormalizer::displayHelp() {
    std::cout << "VCFX_indel_normalizer: Normalize INDEL variants by splitting multi-allelic lines,\n"
              << "and removing common leading/trailing bases to produce a minimal left-aligned representation.\n\n"
              << "Usage:\n"
              << "  VCFX_indel_normalizer [options] < input.vcf > output.vcf\n\n"
              << "Description:\n"
              << "  This code does a simplified left alignment that:\n"
              << "   1) Splits multi-ALT lines into separate lines.\n"
              << "   2) Removes the longest shared prefix from REF/ALT, adjusting POS.\n"
              << "   3) Removes the largest shared suffix from REF/ALT.\n\n"
              << "  Note: true left alignment for repeated motifs requires the full reference genome.\n\n"
              << "Example:\n"
              << "  VCFX_indel_normalizer < input.vcf > normalized.vcf\n";
}

// ---------------------------------------------------------------------
// run
// ---------------------------------------------------------------------
int VCFXIndelNormalizer::run(int argc, char *argv[]) {
    int opt;
    bool showHelp = false;
    static struct option long_opts[] = {{"help", no_argument, 0, 'h'}, {0, 0, 0, 0}};
    while (true) {
        int c = getopt_long(argc, argv, "h", long_opts, NULL);
        if (c == -1)
            break;
        switch (c) {
        case 'h':
        default:
            showHelp = true;
            break;
        }
    }
    if (showHelp) {
        displayHelp();
        return 0;
    }

    normalizeIndels(std::cin, std::cout);
    return 0;
}

// ---------------------------------------------------------------------
// This function does the minimal left alignment approach:
//   1) remove the largest possible prefix (leading) that is identical, but keep at least 1 base
//   2) remove the largest possible suffix (trailing) that is identical, but keep at least 1 base
//   3) adjust pos by the number of removed leading bases
// returns false if after trimming, REF or ALT is empty or they are exactly the same => no variant
// ---------------------------------------------------------------------
bool VCFXIndelNormalizer::normalizeVariant(std::string &chrom, int &posInt, std::string &ref, std::string &alt) {
    if (ref == alt) {
        // no variant
        return false;
    }
    // 1) remove leading common bases
    // we find how many leading chars are identical
    // but we must keep at least 1 char in each
    int prefixCount = 0;
    {
        int minLen = std::min(ref.size(), alt.size());
        while (prefixCount < minLen) {
            if (ref[prefixCount] == alt[prefixCount]) {
                prefixCount++;
            } else {
                break;
            }
        }
    }
    // we only remove prefixCount-1 if prefixCount== length => everything the same
    // typical approach: we can remove (prefixCount-1) from the front. Because the first char of a variant must remain
    // to define the position in a VCF but let's do the approach used by bcftools: we remove (prefixCount -1) as long as
    // prefixCount>0 i.e. if there's a shared leading base, we keep exactly 1. This avoids an empty REF or ALT.
    if (prefixCount > 0) {
        int removeLeading = prefixCount - 1;
        if (removeLeading > 0) {
            ref.erase(0, removeLeading);
            alt.erase(0, removeLeading);
            // adjust pos by removeLeading
            posInt += removeLeading;
        }
    }

    // 2) remove trailing common bases
    // we do a from the end approach
    {
        int rLen = ref.size();
        int aLen = alt.size();
        int suffixCount = 0;
        while (suffixCount < rLen && suffixCount < aLen) {
            if (ref[rLen - 1 - suffixCount] == alt[aLen - 1 - suffixCount]) {
                suffixCount++;
            } else {
                break;
            }
        }
        // we keep at least 1 base, so remove (suffixCount-1) if suffixCount>0
        if (suffixCount > 0) {
            int removeS = suffixCount - 1;
            if (removeS > 0) {
                ref.erase(rLen - removeS, removeS);
                alt.erase(aLen - removeS, removeS);
            }
        }
    }

    // check if after all that, we have no difference or are empty
    if (ref.empty() || alt.empty()) {
        return false;
    }
    if (ref == alt) {
        // no variant => skip
        return false;
    }
    return true;
}

// ---------------------------------------------------------------------
// normalizeIndels: read VCF, print header lines unchanged, then for each line
//   if multi-ALT, split into separate lines. For each alt, do left-trim & right-trim
// ---------------------------------------------------------------------
void VCFXIndelNormalizer::normalizeIndels(std::istream &in, std::ostream &out) {
    std::string line;
    bool foundChromHeader = false;

    while (std::getline(in, line)) {
        if (line.empty()) {
            out << "\n";
            continue;
        }
        if (line[0] == '#') {
            out << line << "\n";
            if (line.rfind("#CHROM", 0) == 0) {
                foundChromHeader = true;
            }
            continue;
        }
        if (!foundChromHeader) {
            // error
            std::cerr << "Error: encountered data line before #CHROM header.\n";
            return;
        }

        // parse the line
        std::stringstream ss(line);
        std::vector<std::string> fields;
        {
            std::string f;
            while (std::getline(ss, f, '\t')) {
                fields.push_back(f);
            }
        }
        if (fields.size() < 10) {
            // not enough columns
            out << line << "\n";
            continue;
        }
        // CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,...samples
        std::string &chrom = fields[0];
        std::string &posStr = fields[1];
        std::string &id = fields[2];
        std::string &ref = fields[3];
        std::string &alt = fields[4];
        // The rest we keep as "rest"
        // We'll rejoin them after we handle alt
        // parse pos as int
        int posInt = 0;
        try {
            posInt = std::stoi(posStr);
        } catch (...) {
            // can't parse => skip
            out << line << "\n";
            continue;
        }
        std::string postCols;
        {
            std::ostringstream oss;
            for (size_t c = 5; c < fields.size(); c++) {
                oss << "\t" << fields[c];
            }
            postCols = oss.str();
        }

        // if alt has commas => multiple alts => split
        // produce multiple lines
        std::vector<std::string> altList;
        {
            std::stringstream altSS(alt);
            std::string a;
            while (std::getline(altSS, a, ',')) {
                altList.push_back(a);
            }
        }
        // for each alt => copy posInt, ref => norm => if success => output
        // if not => output the original or skip?
        if (altList.size() == 1) {
            // single alt => do normal
            std::string altOne = altList[0];
            std::string newRef = ref;
            std::string newAlt = altOne;
            int newPos = posInt;
            bool ok = normalizeVariant(chrom, newPos, newRef, newAlt);
            if (!ok) {
                // if not changed => just print original
                out << line << "\n";
            } else {
                // output normalized
                out << chrom << "\t" << newPos << "\t" << id << "\t" << newRef << "\t" << newAlt << postCols << "\n";
            }
        } else {
            // multiple alt => we produce multiple lines
            for (size_t i = 0; i < altList.size(); i++) {
                std::string altOne = altList[i];
                // clone ref, altOne, pos => norm
                std::string newRef = ref;
                std::string newAlt = altOne;
                int newPos = posInt;
                bool ok = normalizeVariant(chrom, newPos, newRef, newAlt);
                if (!ok) {
                    // if we fail => we output line as is but only with that alt? or skip?
                    // let's do the approach: if it fails => maybe we output the original unnorm
                    // but that would produce duplicates. Let's produce a line anyway
                    // but do not do any modifications => *maybe better is to produce the 'unmodified' altOne
                    out << chrom << "\t" << posInt << "\t" << id << "\t" << ref << "\t" << altOne << postCols << "\n";
                } else {
                    // produce a line
                    // alt is just newAlt, the rest is the same
                    out << chrom << "\t" << newPos << "\t" << id << "\t" << newRef << "\t" << newAlt << postCols
                        << "\n";
                }
            }
        }
    }
}

static void show_help() {
    VCFXIndelNormalizer obj;
    char arg0[] = "VCFX_indel_normalizer";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_indel_normalizer", show_help))
        return 0;
    VCFXIndelNormalizer norm;
    return norm.run(argc, argv);
}
