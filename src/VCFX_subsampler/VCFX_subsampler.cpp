#include "VCFX_subsampler.h"
#include "vcfx_core.h" 
#include "vcfx_io.h"
#include <cstdlib>
#include <ctime>
#include <getopt.h>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

void VCFXSubsampler::displayHelp() {
    std::cout << "VCFX_subsampler: Randomly pick N lines from a VCF data section.\n\n"
                 "Usage:\n"
                 "  VCFX_subsampler [options] < input.vcf > output.vcf\n\n"
                 "Options:\n"
                 "  -s, --subsample <N>   Required: number of data lines (variants) to keep.\n"
                 "  --seed <INT>          Use a reproducible random seed.\n"
                 "  -h, --help            Show this help.\n\n"
                 "Description:\n"
                 "  We read all header lines (#...) first and output them as-is. Then we do\n"
                 "  reservoir sampling on subsequent lines (the data lines). If the file has\n"
                 "  fewer than N lines, we keep them all. We skip lines with <8 columns.\n\n"
                 "Example:\n"
                 "  VCFX_subsampler --subsample 1000 < big.vcf > subset.vcf\n"
                 "  VCFX_subsampler --subsample 1000 --seed 1234 < big.vcf > subset2.vcf\n";
}

int VCFXSubsampler::run(int argc, char *argv[]) {
    if (argc == 1) {
        // No arguments => show help
        displayHelp();
        return 0;
    }

    bool showHelp = false;
    int sampleSize = 0;
    unsigned int seed = static_cast<unsigned int>(std::time(nullptr));

    static struct option long_opts[] = {{"help", no_argument, 0, 'h'},
                                        {"subsample", required_argument, 0, 's'},
                                        {"seed", required_argument, 0, 1000},
                                        {0, 0, 0, 0}};

    while (true) {
        int c = ::getopt_long(argc, argv, "hs:", long_opts, nullptr);
        if (c == -1)
            break;
        switch (c) {
        case 'h':
            showHelp = true;
            break;
        case 's':
            try {
                sampleSize = std::stoi(optarg);
                if (sampleSize <= 0) {
                    throw std::invalid_argument("must be >0");
                }
            } catch (...) {
                std::cerr << "Error: invalid subsample size.\n";
                return 1;
            }
            break;
        case 1000: {
            try {
                long long val = std::stoll(optarg);
                seed = static_cast<unsigned int>(val);
            } catch (...) {
                std::cerr << "Error: invalid seed.\n";
                return 1;
            }
        } break;
        default:
            showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    if (sampleSize <= 0) {
        std::cerr << "Error: must specify --subsample <N> with N>0.\n";
        return 1;
    }

    subsampleLines(std::cin, std::cout, sampleSize, seed);
    return 0;
}

void VCFXSubsampler::subsampleLines(std::istream &in, std::ostream &out, int sampleSize, unsigned int seed) {
    std::string line;
    bool readingHeader = true; // We begin by reading header lines (#...)

    // We'll store lines in a reservoir if they are data lines
    std::vector<std::string> reservoir;
    reservoir.reserve(sampleSize);

    int count = 0; // how many data lines read in total
    std::default_random_engine rng(seed);

    // Read the file line by line
    while (true) {
        if (!std::getline(in, line))
            break; // end of file
        if (line.empty()) {
            // If an empty line, let's just output if we are still in the header,
            // or skip it otherwise. This is optional; typically, empty lines are uncommon in VCF.
            if (readingHeader) {
                out << line << "\n";
            }
            continue;
        }

        // If still reading header lines...
        if (readingHeader && line[0] == '#') {
            // This is a header/comment line => pass it through
            out << line << "\n";
            continue;
        } else {
            // We have encountered the first data line (or a line that isn't '#')
            readingHeader = false; // from now on, everything is data

            // Check column count for data lines
            {
                std::stringstream ss(line);
                std::vector<std::string> fields;
                std::string tmp;
                while (std::getline(ss, tmp, '\t')) {
                    fields.push_back(tmp);
                }
                if (fields.size() < 8) {
                    std::cerr << "Warning: skipping line with <8 columns.\n";
                    continue;
                }
            }

            // Reservoir sampling logic
            if (count < sampleSize) {
                reservoir.push_back(line);
            } else {
                // pick random int in [0..count]
                std::uniform_int_distribution<int> dist(0, count);
                int j = dist(rng);
                if (j < sampleSize) {
                    reservoir[j] = line;
                }
            }
            count++;
        }
    }

    // Output reservoir
    for (const auto &r : reservoir) {
        out << r << "\n";
    }
}

static void show_help() {
    VCFXSubsampler obj;
    char arg0[] = "VCFX_subsampler";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_subsampler", show_help))
        return 0;
    VCFXSubsampler app;
    return app.run(argc, argv);
}
