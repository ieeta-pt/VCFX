#include "VCFX_population_filter.h"
#include "vcfx_core.h" 
#include "vcfx_io.h"
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

int VCFXPopulationFilter::run(int argc, char *argv[]) {
    if (argc == 1) {
        displayHelp();
        return 0;
    }
    bool showHelp = false;
    std::string populationTag;
    std::string popMapFile;

    static struct option long_opts[] = {{"help", no_argument, 0, 'h'},
                                        {"population", required_argument, 0, 'p'},
                                        {"pop-map", required_argument, 0, 'm'},
                                        {0, 0, 0, 0}};

    while (true) {
        int c = ::getopt_long(argc, argv, "hp:m:", long_opts, nullptr);
        if (c == -1)
            break;
        switch (c) {
        case 'h':
            showHelp = true;
            break;
        case 'p':
            populationTag = optarg;
            break;
        case 'm':
            popMapFile = optarg;
            break;
        default:
            showHelp = true;
        }
    }
    if (showHelp) {
        displayHelp();
        return 0;
    }
    if (populationTag.empty() || popMapFile.empty()) {
        std::cerr << "Error: Must specify --population <TAG> and --pop-map <file>.\n";
        displayHelp();
        return 1;
    }

    std::unordered_set<std::string> samplesToInclude;
    if (!loadPopulationMap(popMapFile, populationTag, samplesToInclude)) {
        std::cerr << "Error: Unable to load or parse pop map.\n";
        return 1;
    }
    if (samplesToInclude.empty()) {
        std::cerr << "Warning: No samples found for population tag: " << populationTag << "\n";
    }
    filterPopulation(std::cin, std::cout, samplesToInclude, populationTag);
    return 0;
}

void VCFXPopulationFilter::displayHelp() {
    std::cout << "VCFX_population_filter: Subset VCF to samples in specified population.\n\n"
                 "Usage:\n"
                 "  VCFX_population_filter [options] < input.vcf > output.vcf\n\n"
                 "Options:\n"
                 "  --help, -h               Print this help.\n"
                 "  --population, -p <TAG>   Population tag to keep (e.g. 'EUR','AFR', etc.)\n"
                 "  --pop-map, -m <FILE>     Tab-delimited file: 'SampleName <tab> Population'\n\n"
                 "Description:\n"
                 "  Reads the pop map, finds samples that match the chosen population.\n"
                 "  Then reads the VCF from stdin and prints lines with only those sample columns.\n"
                 "  If a sample is not in that population, it's dropped from the #CHROM header and data columns.\n\n"
                 "Example:\n"
                 "  VCFX_population_filter --population AFR --pop-map pops.txt < input.vcf > out.vcf\n";
}

bool VCFXPopulationFilter::loadPopulationMap(const std::string &popMapFile, const std::string &popTag,
                                             std::unordered_set<std::string> &samplesToInclude) {
    std::ifstream f(popMapFile);
    if (!f.is_open()) {
        std::cerr << "Error: cannot open " << popMapFile << "\n";
        return false;
    }
    std::string line;
    while (true) {
        if (!std::getline(f, line))
            break;
        if (line.empty())
            continue;
        std::stringstream ss(line);
        std::string samp, pop;
        if (!(ss >> samp >> pop)) {
            std::cerr << "Warning: popmap line invalid: " << line << "\n";
            continue;
        }
        if (pop == popTag) {
            samplesToInclude.insert(samp);
        }
    }
    return true;
}

void VCFXPopulationFilter::filterPopulation(std::istream &in, std::ostream &out,
                                            const std::unordered_set<std::string> &samplesToInclude,
                                            const std::string &popTag) {
    bool foundChromLine = false;
    std::string line;
    std::vector<std::string> finalHeader;
    std::vector<int> sampleIndices;
    std::vector<std::string> fields;
    fields.reserve(16);
    std::vector<std::string> columns;
    columns.reserve(16);

    while (true) {
        if (!std::getline(in, line))
            break;
        if (line.empty()) {
            out << line << "\n";
            continue;
        }
        if (line[0] == '#') {
            if (line.rfind("#CHROM", 0) == 0) {
                foundChromLine = true;
                vcfx::split_tabs(line, fields);
                // fields[0..8] => fixed, fields[9..] => samples
                sampleIndices.clear();
                for (size_t i = 9; i < fields.size(); i++) {
                    if (samplesToInclude.count(fields[i]) > 0) {
                        sampleIndices.push_back((int)i);
                    }
                }
                // build new #CHROM line
                std::ostringstream newChrom;
                for (int i = 0; i < 9; i++) {
                    newChrom << fields[i];
                    if (i < 8 || !sampleIndices.empty())
                        newChrom << "\t";
                }
                for (size_t i = 0; i < sampleIndices.size(); i++) {
                    newChrom << fields[sampleIndices[i]];
                    if (i + 1 < sampleIndices.size())
                        newChrom << "\t";
                }
                out << newChrom.str() << "\n";
            } else {
                out << line << "\n";
            }
            continue;
        }
        if (!foundChromLine) {
            std::cerr << "Warning: data line before #CHROM => skipping.\n";
            continue;
        }
        vcfx::split_tabs(line, columns);
        if (columns.size() < 9) {
            std::cerr << "Warning: line with fewer than 9 columns => skipping.\n";
            continue;
        }
        // build new line
        std::ostringstream newLine;
        for (int i = 0; i < 9; i++) {
            newLine << columns[i];
            if (i < 8 || !sampleIndices.empty())
                newLine << "\t";
        }
        for (size_t i = 0; i < sampleIndices.size(); i++) {
            newLine << columns[sampleIndices[i]];
            if (i + 1 < sampleIndices.size())
                newLine << "\t";
        }
        out << newLine.str() << "\n";
    }
    if (!foundChromLine) {
        std::cerr << "Error: No #CHROM header found in VCF.\n";
    }
}

static void show_help() {
    VCFXPopulationFilter obj;
    char arg0[] = "VCFX_population_filter";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_population_filter", show_help))
        return 0;
    VCFXPopulationFilter pf;
    return pf.run(argc, argv);
}
