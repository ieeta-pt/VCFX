#include "VCFX_nonref_filter.h"
#include "vcfx_core.h" 
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

int VCFXNonRefFilter::run(int argc, char *argv[]) {
    bool showHelp = false;
    static struct option long_opts[] = {{"help", no_argument, 0, 'h'}, {0, 0, 0, 0}};
    while (true) {
        int c = getopt_long(argc, argv, "h", long_opts, nullptr);
        if (c == -1)
            break;
        switch (c) {
        case 'h':
        default:
            showHelp = true;
        }
    }
    if (showHelp) {
        displayHelp();
        return 0;
    }
    filterNonRef(std::cin, std::cout);
    return 0;
}

void VCFXNonRefFilter::displayHelp() {
    std::cout << "VCFX_nonref_filter: Exclude variants if all samples are homozygous reference.\n\n"
                 "Usage:\n"
                 "  VCFX_nonref_filter [options] < input.vcf > output.vcf\n\n"
                 "Description:\n"
                 "  Reads VCF lines. For each variant, we check each sample's genotype. If a\n"
                 "  genotype is polyploid, all alleles must be '0'. If a genotype is missing\n"
                 "  or partial, we consider it not guaranteed hom-ref => keep variant.\n"
                 "  If we find at least one sample not hom-ref, we print the variant. Otherwise,\n"
                 "  we skip it.\n\n"
                 "Example:\n"
                 "  VCFX_nonref_filter < input.vcf > filtered.vcf\n\n";
}

bool VCFXNonRefFilter::isDefinitelyHomRef(const std::string &genotypeField) const {
    if (genotypeField.empty() || genotypeField == "." || genotypeField == "./." || genotypeField == ".|.")
        return false;
    std::string g = genotypeField;
    for (char &c : g)
        if (c == '|')
            c = '/';
    // split by '/'
    std::vector<std::string> alleles;
    {
        std::stringstream ss(g);
        std::string tok;
        while (std::getline(ss, tok, '/'))
            alleles.push_back(tok);
    }
    if (alleles.empty())
        return false;
    // if any allele is not "0", => not homRef
    // if allele is "." => missing => not guaranteed homRef => false
    for (auto &al : alleles) {
        if (al != "0")
            return false;
    }
    return true;
}

void VCFXNonRefFilter::filterNonRef(std::istream &in, std::ostream &out) {
    bool headerFound = false;
    std::string line;
    std::vector<std::string> fields;
    fields.reserve(16);
    std::vector<std::string> fmtParts;
    fmtParts.reserve(16);
    std::vector<std::string> subf;
    subf.reserve(16);
    while (true) {
        if (!std::getline(in, line))
            break;
        if (line.empty()) {
            out << line << "\n";
            continue;
        }
        if (line[0] == '#') {
            out << line << "\n";
            if (line.rfind("#CHROM", 0) == 0)
                headerFound = true;
            continue;
        }
        if (!headerFound) {
            std::cerr << "Warning: VCF data line encountered before #CHROM. Passing line.\n";
            out << line << "\n";
            continue;
        }
        vcfx::split_tabs(line, fields);
        if (fields.size() < 10) {
            out << line << "\n";
            continue;
        }
        std::string formatStr = fields[8];
        {
            std::stringstream fs(formatStr);
            std::string ff;
            fmtParts.clear();
            while (std::getline(fs, ff, ':'))
                fmtParts.push_back(ff);
        }
        int gtIndex = -1;
        for (size_t i = 0; i < fmtParts.size(); i++) {
            if (fmtParts[i] == "GT") {
                gtIndex = (int)i;
                break;
            }
        }
        if (gtIndex < 0) {
            // no genotype => cannot confirm all hom-ref => we keep
            out << line << "\n";
            continue;
        }
        bool allHomRef = true;
        for (size_t s = 9; s < fields.size(); s++) {
            std::string &sampleCol = fields[s];
            {
                std::stringstream sampleSS(sampleCol);
                std::string token;
                subf.clear();
                while (std::getline(sampleSS, token, ':'))
                    subf.push_back(token);
            }
            if (gtIndex >= (int)subf.size()) {
                allHomRef = false;
                break;
            }
            if (!isDefinitelyHomRef(subf[gtIndex])) {
                allHomRef = false;
                break;
            }
        }
        if (!allHomRef)
            out << line << "\n";
    }
}

static void show_help() {
    VCFXNonRefFilter obj;
    char arg0[] = "VCFX_nonref_filter";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_nonref_filter", show_help))
        return 0;
    VCFXNonRefFilter app;
    return app.run(argc, argv);
}
