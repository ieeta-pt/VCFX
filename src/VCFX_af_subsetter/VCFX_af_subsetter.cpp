#include "VCFX_af_subsetter.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <sstream>

int VCFXAfSubsetter::run(int argc, char *argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    double minAF = 0.0;
    double maxAF = 1.0;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'}, {"af-filter", required_argument, 0, 'a'}, {0, 0, 0, 0}};

    while ((opt = getopt_long(argc, argv, "ha:", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
            break;
        case 'a': {
            std::string range(optarg);
            size_t dashPos = range.find('-');
            if (dashPos == std::string::npos) {
                std::cerr << "Error: Invalid AF range format. Use <minAF>-<maxAF>.\n";
                displayHelp();
                return 1;
            }
            try {
                minAF = std::stod(range.substr(0, dashPos));
                maxAF = std::stod(range.substr(dashPos + 1));
                if (minAF < 0.0 || maxAF > 1.0 || minAF > maxAF) {
                    throw std::invalid_argument("AF values out of range or minAF > maxAF.");
                }
            } catch (const std::invalid_argument &) {
                std::cerr << "Error: Invalid AF range values. Ensure they are numbers between 0.0 and 1.0 with minAF "
                             "<= maxAF.\n";
                displayHelp();
                return 1;
            }
            break;
        }
        default:
            showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 1;
    }

    // Perform AF-based subsetting from stdin to stdout
    subsetByAlleleFrequency(std::cin, std::cout, minAF, maxAF);

    return 0;
}

void VCFXAfSubsetter::displayHelp() {
    std::cout << "VCFX_af_subsetter: Subset variants based on alternate allele frequency (AF) ranges.\n\n"
              << "Usage:\n"
              << "  VCFX_af_subsetter [options] < input.vcf > output.vcf\n\n"
              << "Options:\n"
              << "  -h, --help                     Display this help message and exit\n"
              << "  -a, --af-filter <minAF>-<maxAF>  Specify the AF range for filtering (e.g., 0.01-0.05)\n\n"
              << "Example:\n"
              << "  VCFX_af_subsetter --af-filter 0.01-0.05 < input.vcf > subsetted.vcf\n";
}

bool VCFXAfSubsetter::parseAF(const std::string &infoField, std::vector<double> &afValues) {
    // Find "AF=" in the INFO string
    size_t pos = infoField.find("AF=");
    if (pos == std::string::npos) {
        return false;
    }

    // Extract substring up to the next semicolon or end of string
    size_t start = pos + 3; // skip 'AF='
    size_t end = infoField.find(';', start);
    std::string afStr = (end == std::string::npos) ? infoField.substr(start) : infoField.substr(start, end - start);

    // Split by comma (to handle multi-allelic AF like "0.2,0.8")
    std::stringstream ss(afStr);
    std::string token;
    while (std::getline(ss, token, ',')) {
        try {
            double afVal = std::stod(token);
            afValues.push_back(afVal);
        } catch (...) {
            return false; // parsing error
        }
    }
    return !afValues.empty();
}

void VCFXAfSubsetter::subsetByAlleleFrequency(std::istream &in, std::ostream &out, double minAF, double maxAF) {
    std::string line;

    // Performance: reuse vector across iterations
    std::vector<std::string> fields;
    fields.reserve(16);
    std::vector<double> afValues;
    afValues.reserve(8);

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        // Print header lines (starting with '#') as-is
        if (line[0] == '#') {
            out << line << "\n";
            continue;
        }

        // Performance: use find-based tab splitting instead of stringstream
        vcfx::split_tabs(line, fields);

        if (fields.size() < 8) {
            std::cerr << "Warning: Skipping invalid VCF line (less than 8 fields): " << line << "\n";
            continue;
        }

        // Parse AF values from the INFO field
        afValues.clear();
        if (!parseAF(fields[7], afValues)) {
            std::cerr << "Warning: AF not found or invalid in INFO field. Skipping variant: " << line << "\n";
            continue;
        }

        // If any AF is within [minAF, maxAF], keep this variant
        bool keepVariant = false;
        for (double af : afValues) {
            if (af >= minAF && af <= maxAF) {
                keepVariant = true;
                break;
            }
        }

        if (keepVariant) {
            out << line << "\n";
        }
    }
}

//
// Typical main():
//
static void show_help() {
    VCFXAfSubsetter obj;
    char arg0[] = "VCFX_af_subsetter";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_af_subsetter", show_help))
        return 0;
    VCFXAfSubsetter afSubsetter;
    return afSubsetter.run(argc, argv);
}
