#include "VCFX_indexer.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm> // for std::find_if_not, if you want a neat trim
#include <cstdint>
#include <cstring>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// A helper to split a string by literal '\t' characters
static std::vector<std::string> splitTabs(const std::string &s) {
    std::vector<std::string> tokens;
    size_t start = 0;
    while (true) {
        size_t pos = s.find('\t', start);
        if (pos == std::string::npos) {
            tokens.push_back(s.substr(start));
            break;
        }
        tokens.push_back(s.substr(start, pos - start));
        start = pos + 1;
    }
    return tokens;
}

// Trim leading whitespace
static std::string ltrim(const std::string &str) {
    size_t startPos = 0;
    while (startPos < str.size() && std::isspace(static_cast<unsigned char>(str[startPos]))) {
        startPos++;
    }
    return str.substr(startPos);
}

void VCFXIndexer::displayHelp() {
    std::cout << "VCFX_indexer\n"
              << "Usage: VCFX_indexer [--help]\n\n"
              << "Description:\n"
              << "  Reads a VCF from stdin (raw bytes) and writes a 3-column index\n"
              << "  (CHROM, POS, FILE_OFFSET) to stdout. FILE_OFFSET is the byte offset\n"
              << "  from the start of the file to the beginning of each variant line.\n\n"
              << "Example:\n"
              << "  VCFX_indexer < input.vcf > index.tsv\n";
}

int VCFXIndexer::run(int argc, char *argv[]) {
    static struct option long_opts[] = {{"help", no_argument, 0, 'h'}, {0, 0, 0, 0}};

    while (true) {
        int c = getopt_long(argc, argv, "h", long_opts, nullptr);
        if (c == -1)
            break;
        switch (c) {
        case 'h':
            displayHelp();
            return 0;
        default:
            displayHelp();
            return 1;
        }
    }

    // Perform the actual indexing
    createVCFIndex(std::cin, std::cout);
    return 0;
}

void VCFXIndexer::createVCFIndex(std::istream &in, std::ostream &out) {
    bool foundChromHeader = false;
    bool warnedNoChromYet = false;
    bool sawAnyHeaderLine = false;

    static const size_t BUF_SIZE = 64 * 1024;
    char buffer[BUF_SIZE];

    // Accumulate partial lines here
    std::string leftover;
    leftover.reserve(BUF_SIZE);

    // Tracks total raw bytes read from the start
    std::int64_t totalOffset = 0;

    // Helper lambda to process each full line
    auto processLine = [&](const std::string &line, std::int64_t lineStartOffset) {
        if (line.empty()) {
            // Skip empty lines
            return;
        }

        // Trim leading space
        std::string trimmed = ltrim(line);

        // If it starts with '#', treat as header line
        if (!trimmed.empty() && trimmed[0] == '#') {
            sawAnyHeaderLine = true;

            // Possibly #CHROM
            if (!foundChromHeader) {
                auto fields = splitTabs(trimmed);
                if (fields.size() >= 2 && fields[0] == "#CHROM") {
                    foundChromHeader = true;
                    // Once #CHROM is found, print the header
                    out << "CHROM\tPOS\tFILE_OFFSET\n";
                }
            }
            return;
        }

        // It's a data line. If we haven't seen #CHROM, that's an error
        if (!foundChromHeader) {
            if (!sawAnyHeaderLine && !warnedNoChromYet) {
                std::cerr << "Error: no #CHROM header found before variant lines.\n";
                warnedNoChromYet = true;
            }
            return;
        }

        // We do have #CHROM; parse the fields
        auto fields = splitTabs(line);
        if (fields.size() < 2) {
            return; // Not enough fields to parse
        }

        const std::string &chrom = fields[0];
        const std::string &posStr = fields[1];

        std::int64_t posVal = 0;
        try {
            posVal = std::stoll(posStr);
        } catch (...) {
            // Not a valid integer => skip
            return;
        }

        // Valid line => output the index info
        out << chrom << "\t" << posVal << "\t" << lineStartOffset << "\n";
    };

    // Main read loop
    while (true) {
        in.read(buffer, BUF_SIZE);
        std::streamsize got = in.gcount();
        if (got <= 0) {
            break;
        }

        for (std::streamsize i = 0; i < got; ++i) {
            char c = buffer[i];
            leftover.push_back(c);
            totalOffset++;

            if (c == '\n') {
                // A full line has ended
                leftover.pop_back(); // remove the '\n'
                std::int64_t lineStartOffset = totalOffset - leftover.size() - 1;

                // Handle CRLF by removing trailing '\r'
                if (!leftover.empty() && leftover.back() == '\r') {
                    leftover.pop_back();
                }

                // Process the line
                processLine(leftover, lineStartOffset);
                leftover.clear();
            }
        }

        if (!in.good()) {
            // If stream ended or had an error
            break;
        }
    }

    // Handle any leftover partial line
    if (!leftover.empty()) {
        std::int64_t lineStartOffset = totalOffset - static_cast<std::int64_t>(leftover.size());
        if (!leftover.empty() && leftover.back() == '\r') {
            leftover.pop_back();
        }
        processLine(leftover, lineStartOffset);
        leftover.clear();
    }
}

// Optional main if you build as a single executable
static void show_help() {
    VCFXIndexer obj;
    char arg0[] = "VCFX_indexer";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_indexer", show_help))
        return 0;
    VCFXIndexer idx;
    return idx.run(argc, argv);
}
