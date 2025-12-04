#include "VCFX_fasta_converter.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <getopt.h>
#include <iomanip>
#include <map>
#include <sstream>
#include <unordered_map>

// A small map from two distinct bases to an IUPAC ambiguity code
// e.g. A + G => R, C + T => Y, etc.
static const std::map<std::string, char> IUPAC_ambiguities = {{"AG", 'R'}, {"GA", 'R'}, {"CT", 'Y'}, {"TC", 'Y'},
                                                              {"AC", 'M'}, {"CA", 'M'}, {"GT", 'K'}, {"TG", 'K'},
                                                              {"AT", 'W'}, {"TA", 'W'}, {"CG", 'S'}, {"GC", 'S'}};

// Utility to convert a single numeric allele index into the corresponding base
// returns '\0' on failure
static char alleleIndexToBase(int alleleIndex, const std::string &ref, const std::vector<std::string> &altAlleles) {
    // 0 => ref, 1 => altAlleles[0], 2 => altAlleles[1], etc.
    if (alleleIndex == 0) {
        if (ref.size() == 1) {
            return std::toupper(ref[0]);
        } else {
            // multi-base or invalid for a single-locus representation
            return '\0';
        }
    } else {
        int altPos = alleleIndex - 1;
        if (altPos < 0 || (size_t)altPos >= altAlleles.size()) {
            return '\0'; // out of range
        }
        // altAlleles[altPos] must be a single base to be representable
        const std::string &a = altAlleles[altPos];
        if (a.size() == 1) {
            return std::toupper(a[0]);
        } else {
            // multi-base alt => can't represent as single base
            return '\0';
        }
    }
}

// If we have exactly two bases, see if there's a standard IUPAC code
// Otherwise returns 'N'
static char combineBasesIUPAC(char b1, char b2) {
    if (b1 == b2) {
        return b1; // e.g. A + A => A
    }
    // build 2-char string in alphabetical order
    std::string pair;
    pair.push_back(std::min(b1, b2));
    pair.push_back(std::max(b1, b2));
    auto it = IUPAC_ambiguities.find(pair);
    if (it != IUPAC_ambiguities.end()) {
        return it->second;
    }
    return 'N'; // unknown combination
}

int VCFXFastaConverter::run(int argc, char *argv[]) {
    // Parse command-line arguments
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

    // Convert VCF input from stdin to FASTA output
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

void VCFXFastaConverter::convertVCFtoFasta(std::istream &in, std::ostream &out) {
    std::string line;
    std::vector<std::string> sampleNames;
    // Each sampleName -> sequence string
    std::unordered_map<std::string, std::string> sampleSequences;

    bool headerParsed = false;

    // Reusable buffer for parsing
    std::vector<std::string> fields;
    fields.reserve(16);

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            // Parse the #CHROM header to get sample columns
            if (line.rfind("#CHROM", 0) == 0) {
                std::vector<std::string> headers;
                vcfx::split_tabs(line, headers);
                // Sample names start at column 9
                for (size_t i = 9; i < headers.size(); ++i) {
                    sampleNames.push_back(headers[i]);
                    sampleSequences[headers[i]] = "";
                }
                headerParsed = true;
            }
            continue;
        }

        if (!headerParsed) {
            std::cerr << "Error: #CHROM header not found before data lines.\n";
            return;
        }

        // Parse data line
        fields.clear();
        vcfx::split_tabs(line, fields);
        // minimal check
        if (fields.size() < (9 + sampleNames.size())) {
            // not enough columns
            std::cerr << "Warning: Skipping malformed VCF line with insufficient columns.\n";
            continue;
        }

        // VCF standard columns
        //  0:CHROM, 1:POS, 2:ID, 3:REF, 4:ALT, 5:QUAL, 6:FILTER, 7:INFO, 8:FORMAT, 9+:samples
        const std::string &chrom = fields[0];
        // const std::string &pos = fields[1]; // not strictly needed for the FASTA
        // const std::string &id  = fields[2];
        const std::string &ref = fields[3];
        const std::string &altField = fields[4];
        const std::string &format = fields[8];

        // Split alt on commas
        std::vector<std::string> altAlleles;
        {
            std::stringstream altSS(altField);
            std::string a;
            while (std::getline(altSS, a, ',')) {
                altAlleles.push_back(a);
            }
        }

        // Find GT index in format
        std::vector<std::string> formatFields;
        {
            std::stringstream fmts(format);
            std::string token;
            while (std::getline(fmts, token, ':')) {
                formatFields.push_back(token);
            }
        }
        int gtIndex = -1;
        for (size_t i = 0; i < formatFields.size(); ++i) {
            if (formatFields[i] == "GT") {
                gtIndex = static_cast<int>(i);
                break;
            }
        }
        bool hasGT = (gtIndex >= 0);

        // For each sample, figure out their genotype => one base or IUPAC or N
        for (size_t s = 0; s < sampleNames.size(); ++s) {
            const std::string &sampleName = sampleNames[s];
            // sample column is at fields[9 + s]
            if (9 + s >= fields.size()) {
                // missing sample column?
                sampleSequences[sampleName] += "N";
                continue;
            }
            const std::string &sampleData = fields[9 + s];
            if (!hasGT) {
                // no genotype => default to reference? or N?
                // We'll choose 'N' to avoid assumptions
                sampleSequences[sampleName] += "N";
                continue;
            }
            // parse sampleData by ':'
            std::vector<std::string> sampleParts;
            {
                std::stringstream sp(sampleData);
                std::string p;
                while (std::getline(sp, p, ':')) {
                    sampleParts.push_back(p);
                }
            }
            if (gtIndex >= (int)sampleParts.size()) {
                sampleSequences[sampleName] += "N";
                continue;
            }
            // genotype string
            std::string genotype = sampleParts[gtIndex];
            // unify separators
            for (char &c : genotype) {
                if (c == '|')
                    c = '/';
            }
            if (genotype.empty() || genotype == ".") {
                sampleSequences[sampleName] += "N";
                continue;
            }
            // split genotype by '/'
            std::vector<std::string> alleles;
            {
                std::stringstream gtSS(genotype);
                std::string al;
                while (std::getline(gtSS, al, '/')) {
                    alleles.push_back(al);
                }
            }
            if (alleles.size() != 2) {
                // not diploid => 'N'
                sampleSequences[sampleName] += "N";
                continue;
            }

            // Convert each allele to a single base (char), or '\0' on fail
            char b1 = '\0';
            char b2 = '\0';
            {
                // if either is '.', skip
                if (alleles[0] == "." || alleles[1] == ".") {
                    sampleSequences[sampleName] += "N";
                    continue;
                }
                // parse numeric
                bool okA1 = true, okA2 = true;
                int a1 = 0, a2 = 0;
                // parse first allele
                {
                    for (char c : alleles[0]) {
                        if (!std::isdigit(c)) {
                            okA1 = false;
                            break;
                        }
                    }
                    if (okA1)
                        a1 = std::stoi(alleles[0]);
                }
                // parse second allele
                {
                    for (char c : alleles[1]) {
                        if (!std::isdigit(c)) {
                            okA2 = false;
                            break;
                        }
                    }
                    if (okA2)
                        a2 = std::stoi(alleles[1]);
                }
                if (!okA1 || !okA2) {
                    sampleSequences[sampleName] += "N";
                    continue;
                }
                b1 = alleleIndexToBase(a1, ref, altAlleles);
                b2 = alleleIndexToBase(a2, ref, altAlleles);
            }

            if (b1 == '\0' || b2 == '\0') {
                // means we couldn't interpret at least one allele
                sampleSequences[sampleName] += "N";
                continue;
            }

            // If same base => that base
            // If different => try IUPAC
            char finalBase = '\0';
            if (b1 == b2) {
                finalBase = b1; // e.g. both 'A'
            } else {
                finalBase = combineBasesIUPAC(b1, b2); // might yield R, Y, etc., or 'N'
            }
            if (finalBase == '\0')
                finalBase = 'N';
            sampleSequences[sampleName] += finalBase;
        }
    }

    // Finally, output the sequences in FASTA format
    // e.g. >SampleName\n[sequence in 60-char lines]
    for (auto &kv : sampleSequences) {
        const std::string &sampleName = kv.first;
        const std::string &seq = kv.second;
        out << ">" << sampleName << "\n";
        // print in 60-char chunks
        for (size_t i = 0; i < seq.size(); i += 60) {
            out << seq.substr(i, 60) << "\n";
        }
    }
}

static void show_help() {
    VCFXFastaConverter obj;
    char arg0[] = "VCFX_fasta_converter";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_fasta_converter", show_help))
        return 0;
    VCFXFastaConverter app;
    return app.run(argc, argv);
}
