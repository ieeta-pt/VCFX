#include "VCFX_header_parser.h"
#include "vcfx_core.h" 
#include "vcfx_io.h"
#include <iostream>
#include <sstream>

void printHelp() {
    std::cout << "VCFX_header_parser\n"
              << "Usage: VCFX_header_parser [OPTIONS]\n\n"
              << "Options:\n"
              << "  --help, -h            Display this help message and exit.\n"
              << "\n"
              << "Description:\n"
              << "  Extracts and displays the header lines from a VCF file.\n\n"
              << "Example:\n"
              << "  ./VCFX_header_parser < input.vcf > header.txt\n";
}

void processHeader(std::istream &in, std::ostream &out) {
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '#') {
            out << line << std::endl;
        } else {
            break; // Stop reading after header
        }
    }
}

static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_header_parser", show_help))
        return 0;
    // Simple argument parsing
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
    }

    // Process header
    processHeader(std::cin, std::cout);
    return 0;
}
