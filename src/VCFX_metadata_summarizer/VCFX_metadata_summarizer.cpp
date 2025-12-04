#include "VCFX_metadata_summarizer.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <sstream>

// A helper function to extract the value of ID=... from a line like
// ##contig=<ID=chr1,length=...>
static std::string extractID(const std::string &line, const std::string &prefix) {
    // We expect something like '##contig=<ID=xxx,'
    // find prefix => e.g. "##contig=" -> line.find(prefix) might be 0
    // Then find 'ID='
    // Then read until ',' or '>' or end
    auto idPos = line.find("ID=");
    if (idPos == std::string::npos) {
        return "";
    }
    // substring from idPos+3
    auto sub = line.substr(idPos + 3);
    // if sub starts with something => we read until ',' or '>'
    // find first ',' or '>'
    size_t endPos = sub.find_first_of(",>");
    if (endPos == std::string::npos) {
        // just return entire
        return sub;
    }
    return sub.substr(0, endPos);
}

// Implementation of VCFXMetadataSummarizer
int VCFXMetadataSummarizer::run(int argc, char *argv[]) {
    int opt;
    bool showHelp = false;

    static struct option long_options[] = {{"help", no_argument, 0, 'h'}, {0, 0, 0, 0}};

    while ((opt = getopt_long(argc, argv, "h", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
            break;
        default:
            showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    // Summarize metadata from stdin
    summarizeMetadata(std::cin);

    return 0;
}

void VCFXMetadataSummarizer::displayHelp() {
    std::cout << "VCFX_metadata_summarizer: Summarize key metadata (contigs, INFO, FILTER, FORMAT, samples, variants) "
                 "from a VCF file.\n\n"
              << "Usage:\n"
              << "  VCFX_metadata_summarizer [options] < input.vcf\n\n"
              << "Options:\n"
              << "  -h, --help   Display this help message and exit\n\n"
              << "Example:\n"
              << "  VCFX_metadata_summarizer < input.vcf\n";
}

void VCFXMetadataSummarizer::summarizeMetadata(std::istream &in) {
    std::string line;

    while (true) {
        if (!std::getline(in, line))
            break;
        if (line.empty())
            continue;

        if (line[0] == '#') {
            // if it's a meta line => parse
            if (line.rfind("##", 0) == 0) {
                parseHeader(line);
            }
            // if it's #CHROM => parse sample names
            else if (line.rfind("#CHROM", 0) == 0) {
                // parse columns
                std::stringstream ss(line);
                std::vector<std::string> fields;
                {
                    std::string f;
                    while (std::getline(ss, f, '\t')) {
                        fields.push_back(f);
                    }
                }
                // from col=9 onward => samples
                if (fields.size() > 9) {
                    this->numSamples = (int)fields.size() - 9;
                } else {
                    this->numSamples = 0;
                }
            }
        } else {
            // data line => count a variant
            numVariants++;
        }
    }

    // now print summary
    printSummary();
}

void VCFXMetadataSummarizer::parseHeader(const std::string &line) {
    // check type of line => if contig, info, filter, format
    // e.g. "##contig=<ID=chr1,length=...>" => parse ID
    // "##INFO=<ID=DP,Number=1,Type=Integer,...>"
    // "##FILTER=<ID=LowQual,Description=...>"
    // "##FORMAT=<ID=GT,Number=1,Type=String,...>"

    if (line.find("##contig=") != std::string::npos) {
        auto idStr = extractID(line, "##contig=");
        if (!idStr.empty()) {
            contigIDs.insert(idStr);
        }
    } else if (line.find("##INFO=") != std::string::npos) {
        auto idStr = extractID(line, "##INFO=");
        if (!idStr.empty()) {
            infoIDs.insert(idStr);
        }
    } else if (line.find("##FILTER=") != std::string::npos) {
        auto idStr = extractID(line, "##FILTER=");
        if (!idStr.empty()) {
            filterIDs.insert(idStr);
        }
    } else if (line.find("##FORMAT=") != std::string::npos) {
        auto idStr = extractID(line, "##FORMAT=");
        if (!idStr.empty()) {
            formatIDs.insert(idStr);
        }
    }
}

void VCFXMetadataSummarizer::printSummary() const {
    std::cout << "VCF Metadata Summary:\n";
    std::cout << "---------------------\n";
    std::cout << "Number of unique contigs: " << contigIDs.size() << "\n";
    std::cout << "Number of unique INFO fields: " << infoIDs.size() << "\n";
    std::cout << "Number of unique FILTER fields: " << filterIDs.size() << "\n";
    std::cout << "Number of unique FORMAT fields: " << formatIDs.size() << "\n";
    std::cout << "Number of samples: " << numSamples << "\n";
    std::cout << "Number of variants: " << numVariants << "\n";
}

static void show_help() {
    VCFXMetadataSummarizer obj;
    char arg0[] = "VCFX_metadata_summarizer";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_metadata_summarizer", show_help))
        return 0;
    VCFXMetadataSummarizer summarizer;
    return summarizer.run(argc, argv);
}
