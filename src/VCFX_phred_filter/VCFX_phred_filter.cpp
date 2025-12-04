#include "VCFX_phred_filter.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>

int VCFXPhredFilter::run(int argc, char *argv[]) {
    if (argc == 1) {
        displayHelp();
        return 0;
    }
    double threshold = 30.0;
    bool showHelp = false;
    bool keepMissingAsPass = false;
    static struct option long_opts[] = {{"phred-filter", required_argument, 0, 'p'},
                                        {"keep-missing-qual", no_argument, 0, 'k'},
                                        {"help", no_argument, 0, 'h'},
                                        {0, 0, 0, 0}};
    while (true) {
        int c = ::getopt_long(argc, argv, "p:kh", long_opts, nullptr);
        if (c == -1)
            break;
        switch (c) {
        case 'p': {
            try {
                threshold = std::stod(optarg);
            } catch (...) {
                std::cerr << "Error: Invalid threshold '" << optarg << "'.\n";
                return 1;
            }
        } break;
        case 'k':
            keepMissingAsPass = true;
            break;
        case 'h':
        default:
            showHelp = true;
        }
    }
    if (showHelp) {
        displayHelp();
        return 0;
    }
    processVCF(std::cin, threshold, keepMissingAsPass);
    return 0;
}

void VCFXPhredFilter::displayHelp() {
    std::cout << "VCFX_phred_filter: Filter VCF lines by their QUAL field.\n\n"
                 "Usage:\n"
                 "  VCFX_phred_filter [options] < input.vcf > output.vcf\n\n"
                 "Options:\n"
                 "  -p, --phred-filter <VAL>      Phred QUAL threshold (default=30)\n"
                 "  -k, --keep-missing-qual       Treat '.' (missing QUAL) as pass\n"
                 "  -h, --help                    Display this help and exit\n\n"
                 "Description:\n"
                 "  Reads VCF lines from stdin. For each data line, parse the QUAL field.\n"
                 "  If QUAL >= threshold => print line. Otherwise, skip. By default, missing\n"
                 "  QUAL ('.') is treated as 0. Use --keep-missing-qual to treat '.' as pass.\n\n"
                 "Examples:\n"
                 "  1) Keep variants with QUAL>=30:\n"
                 "     VCFX_phred_filter -p 30 < in.vcf > out.vcf\n"
                 "  2) Keep missing QUAL lines:\n"
                 "     VCFX_phred_filter -p 30 --keep-missing-qual < in.vcf > out.vcf\n";
}

void VCFXPhredFilter::processVCF(std::istream &in, double threshold, bool keepMissingAsPass) {
    std::string line;
    std::vector<std::string> fields;
    fields.reserve(16);  // Performance: pre-allocate for typical VCF
    bool foundChrom = false;
    while (true) {
        if (!std::getline(in, line))
            break;
        if (line.empty()) {
            std::cout << line << "\n";
            continue;
        }
        if (line[0] == '#') {
            std::cout << line << "\n";
            if (line.rfind("#CHROM", 0) == 0)
                foundChrom = true;
            continue;
        }
        if (!foundChrom) {
            std::cerr << "Warning: data line before #CHROM => skipping line.\n";
            continue;
        }
        vcfx::split_tabs(line, fields);  // Performance: reuses vector capacity
        // we need at least CHROM,POS,ID,REF,ALT,QUAL => 6 columns
        if (fields.size() < 6) {
            std::cerr << "Warning: line has <6 columns => skipping.\n";
            continue;
        }
        double q = parseQUAL(fields[5], keepMissingAsPass);
        if (q >= threshold) {
            std::cout << line << "\n";
        }
    }
}

double VCFXPhredFilter::parseQUAL(const std::string &qualStr, bool keepMissingAsPass) {
    if (qualStr == "." || qualStr.empty()) {
        if (keepMissingAsPass)
            return 1e9;
        else
            return 0.0;
    }
    try {
        return std::stod(qualStr);
    } catch (...) {
        std::cerr << "Warning: Invalid QUAL '" << qualStr << "'. Using 0.\n";
        return 0.0;
    }
}

static void show_help() {
    VCFXPhredFilter obj;
    char arg0[] = "VCFX_phred_filter";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_phred_filter", show_help))
        return 0;
    VCFXPhredFilter pf;
    return pf.run(argc, argv);
}
