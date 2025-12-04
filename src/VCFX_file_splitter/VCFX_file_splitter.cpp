#include "VCFX_file_splitter.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <memory>
#include <sstream>
#include <unordered_map>
#include <vector>

struct ChromFile {
    std::unique_ptr<std::ofstream> ofs;
    bool headerWritten;
};

int VCFXFileSplitter::run(int argc, char *argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string outputPrefix = "split";

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'}, {"prefix", required_argument, 0, 'p'}, {0, 0, 0, 0}};

    while ((opt = getopt_long(argc, argv, "hp:", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
            break;
        case 'p':
            outputPrefix = optarg;
            break;
        default:
            showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    splitVCFByChromosome(std::cin, outputPrefix);
    return 0;
}

void VCFXFileSplitter::displayHelp() {
    std::cout << "VCFX_file_splitter: Split a VCF file into multiple files based on chromosome.\n\n"
              << "Usage:\n"
              << "  VCFX_file_splitter [options] < input.vcf\n\n"
              << "Options:\n"
              << "  -h, --help            Display this help message and exit\n"
              << "  -p, --prefix <prefix> Output file prefix (default: 'split')\n\n"
              << "Example:\n"
              << "  VCFX_file_splitter --prefix \"chr\" < input.vcf\n";
}

// Splits the VCF by chromosome, writing the full header to each file.
void VCFXFileSplitter::splitVCFByChromosome(std::istream &in, const std::string &outputPrefix) {
    std::unordered_map<std::string, ChromFile> chromFiles;

    // We'll store all lines that begin with '#' (the header lines) until we hit
    // the first data line. We also handle the possibility of additional '#' lines
    // that appear after #CHROM, replicating them to all open chromosome files.
    std::vector<std::string> initialHeaderLines;
    bool foundFirstDataLine = false;

    // We'll read the entire file line by line, but we do two main phases:
    //   1) Collect all '#' lines (the main header). Once we see a data line,
    //      we treat it as phase 2, splitting lines by chromosome.
    //   2) If we encounter extra '#' lines in phase 2, we replicate them
    //      to all already-open files.

    std::string line;
    while (true) {
        std::streampos currentPos = in.tellg();
        if (!std::getline(in, line)) {
            // end of file
            break;
        }
        if (line.empty()) {
            continue;
        }

        if (line[0] == '#') {
            if (!foundFirstDataLine) {
                // We haven't yet seen any data lines, so this is part of the "initial" header
                initialHeaderLines.push_back(line);
            } else {
                // This is an extra header line after data lines started, replicate to open files
                for (auto &kv : chromFiles) {
                    ChromFile &cf = kv.second;
                    if (cf.ofs && cf.ofs->is_open()) {
                        *(cf.ofs) << line << "\n";
                    }
                }
                // Also add to initial header lines for any future chromosome files
                initialHeaderLines.push_back(line);
            }
        } else {
            // This is a data line
            if (!foundFirstDataLine) {
                // We just encountered the first data line
                // so from now on, we are in the "splitting" phase
                foundFirstDataLine = true;
            }
            // We parse the chromosome from the data line
            std::vector<std::string> fields;
            vcfx::split_tabs(line, fields);
            if (fields.empty()) {
                std::cerr << "Warning: cannot parse CHROM from line: " << line << "\n";
                continue;
            }
            const std::string &chrom = fields[0];
            // Check or create file
            if (chromFiles.find(chrom) == chromFiles.end()) {
                // Create a new file
                std::string filename = outputPrefix + "_" + chrom + ".vcf";
                ChromFile cf;
                cf.ofs = std::make_unique<std::ofstream>(filename);
                if (!cf.ofs->is_open()) {
                    std::cerr << "Error: Unable to create file: " << filename << "\n";
                    continue;
                }
                cf.headerWritten = false;
                chromFiles[chrom] = std::move(cf);
            }
            ChromFile &cf = chromFiles[chrom];
            if (!cf.headerWritten) {
                // Write all initial header lines
                for (auto &hLine : initialHeaderLines) {
                    *(cf.ofs) << hLine << "\n";
                }
                cf.headerWritten = true;
            }
            // Write this data line
            *(cf.ofs) << line << "\n";
        }
    }

    // Close all files
    for (auto &kv : chromFiles) {
        ChromFile &cf = kv.second;
        if (cf.ofs && cf.ofs->is_open()) {
            cf.ofs->close();
        }
    }

    // If we never found any data line, it's possible the input was all header or empty
    // In that case, we do nothing or we could emit a warning
    if (!foundFirstDataLine) {
        // Possibly we can check if initialHeaderLines is non-empty => maybe we do something
        // But there's no data to split
        std::cerr << "Note: No variant data lines were found in the input.\n";
    }
}

static void show_help() {
    VCFXFileSplitter obj;
    char arg0[] = "VCFX_file_splitter";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_file_splitter", show_help))
        return 0;
    VCFXFileSplitter splitter;
    return splitter.run(argc, argv);
}
