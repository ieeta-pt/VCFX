#include "VCFX_position_subsetter.h"
#include "vcfx_core.h" 
#include "vcfx_io.h"
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

int VCFXPositionSubsetter::run(int argc, char *argv[]) {
    if (argc == 1) {
        displayHelp();
        return 0;
    }
    bool showHelp = false;
    std::string regionStr;

    static struct option long_opts[] = {
        {"region", required_argument, 0, 'r'}, {"help", no_argument, 0, 'h'}, {0, 0, 0, 0}};

    while (true) {
        int c = ::getopt_long(argc, argv, "r:h", long_opts, nullptr);
        if (c == -1)
            break;
        switch (c) {
        case 'r':
            regionStr = optarg;
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
    if (regionStr.empty()) {
        std::cerr << "Error: --region <chrX:start-end> is required.\n";
        displayHelp();
        return 1;
    }
    std::string chrom;
    int start = 0, end = 0;
    if (!parseRegion(regionStr, chrom, start, end)) {
        return 1;
    }
    bool ok = subsetVCFByPosition(std::cin, std::cout, chrom, start, end);
    return ok ? 0 : 1;
}

void VCFXPositionSubsetter::displayHelp() {
    std::cout << "VCFX_position_subsetter: Subset VCF by a single genomic region.\n\n"
                 "Usage:\n"
                 "  VCFX_position_subsetter --region \"chr1:10000-20000\" < in.vcf > out.vcf\n\n"
                 "Options:\n"
                 "  -r, --region \"CHR:START-END\"   The region to keep.\n"
                 "  -h, --help                     Print this help.\n\n"
                 "Description:\n"
                 "  Reads lines from VCF input, and only prints data lines where:\n"
                 "    1) CHROM matches 'CHR' exactly, and\n"
                 "    2) POS is in [START,END].\n"
                 "  All header lines (#...) are passed unmodified.\n\n"
                 "Example:\n"
                 "  VCFX_position_subsetter --region \"chr2:500-1000\" < input.vcf > subset.vcf\n";
}

bool VCFXPositionSubsetter::parseRegion(const std::string &regionStr, std::string &chrom, int &start, int &end) {
    // find colon, dash
    auto cpos = regionStr.find(':');
    auto dpos = regionStr.find('-');
    if (cpos == std::string::npos || dpos == std::string::npos || dpos <= cpos) {
        std::cerr << "Error: invalid region: " << regionStr << ". Expected e.g. chr1:10000-20000.\n";
        return false;
    }
    chrom = regionStr.substr(0, cpos);
    std::string startStr = regionStr.substr(cpos + 1, dpos - (cpos + 1));
    std::string endStr = regionStr.substr(dpos + 1);
    try {
        start = std::stoi(startStr);
        end = std::stoi(endStr);
    } catch (...) {
        std::cerr << "Error: cannot parse region start/end.\n";
        return false;
    }
    if (start > end) {
        std::cerr << "Error: region start> end.\n";
        return false;
    }
    return true;
}

bool VCFXPositionSubsetter::subsetVCFByPosition(std::istream &in, std::ostream &out, const std::string &regionChrom,
                                                int regionStart, int regionEnd) {
    bool headerFound = false;
    std::string line;
    std::vector<std::string> fields;
    fields.reserve(16);
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
            std::cerr << "Warning: data line encountered before #CHROM => skipping.\n";
            continue;
        }
        vcfx::split_tabs(line, fields);
        if (fields.size() < 2) {
            std::cerr << "Warning: line has <2 columns => skipping.\n";
            continue;
        }
        const std::string &chrom = fields[0];
        const std::string &posStr = fields[1];
        int pos = 0;
        try {
            pos = std::stoi(posStr);
        } catch (...) {
            std::cerr << "Warning: invalid POS '" << posStr << "'. Skipping.\n";
            continue;
        }
        if (chrom == regionChrom && pos >= regionStart && pos <= regionEnd) {
            out << line << "\n";
        }
    }
    return true;
}

static void show_help() {
    VCFXPositionSubsetter obj;
    char arg0[] = "VCFX_position_subsetter";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_position_subsetter", show_help))
        return 0;
    VCFXPositionSubsetter subsetter;
    return subsetter.run(argc, argv);
}
