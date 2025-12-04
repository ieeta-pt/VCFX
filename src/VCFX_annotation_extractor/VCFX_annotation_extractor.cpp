#include "vcfx_core.h" 
#include "vcfx_io.h"
#include <algorithm>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// --------------------------------------------------------------
// A small struct to store command-line options
// --------------------------------------------------------------
struct AnnotationOptions {
    std::vector<std::string> annotations; // e.g. ["ANN", "Gene"]
};

// --------------------------------------------------------------
// Utility: split a string by a delimiter into a vector
// --------------------------------------------------------------
static std::vector<std::string> split(const std::string &str, char delim) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delim)) {
        tokens.push_back(token);
    }
    return tokens;
}

// --------------------------------------------------------------
// Utility: parse the INFO field into a map key->value
//   e.g. "ANN=xxx;Gene=YYY;DP=100" => {ANN:xxx, Gene:YYY, DP:100}
// --------------------------------------------------------------
static std::unordered_map<std::string, std::string> parseInfoToMap(const std::string &info) {
    std::unordered_map<std::string, std::string> infoMap;
    // split by ';'
    auto fields = split(info, ';');
    for (auto &f : fields) {
        if (f.empty()) {
            continue;
        }
        // e.g. f="ANN=..."
        auto eqPos = f.find('=');
        if (eqPos == std::string::npos) {
            // key without value? e.g. "SOMATIC"
            // You could store it as {SOMATIC: ""} if you want
            infoMap[f] = "";
            continue;
        }
        std::string key = f.substr(0, eqPos);
        std::string val = f.substr(eqPos + 1);
        infoMap[key] = val;
    }
    return infoMap;
}

// --------------------------------------------------------------
// Show usage/help
// --------------------------------------------------------------
static void printHelp() {
    std::cout << "VCFX_annotation_extractor: Extract variant annotations from a VCF file.\n\n"
              << "Usage:\n"
              << "  VCFX_annotation_extractor --annotation-extract \"ANN,Gene\" < input.vcf > out.tsv\n\n"
              << "Options:\n"
              << "  -a, --annotation-extract   Comma-separated list of annotations to extract (e.g., ANN,Gene)\n"
              << "  -h, --help                 Display this help message and exit\n\n"
              << "Description:\n"
              << "  Reads a VCF from stdin and prints a tab-delimited output. For multi-ALT\n"
              << "  lines, each ALT allele is printed on its own line. If an annotation field (like\n"
              << "  'ANN=') has multiple comma-separated sub-entries, we attempt to align them with\n"
              << "  the ALT alleles in order.\n\n"
              << "Example:\n"
              << "  VCFX_annotation_extractor --annotation-extract \"ANN,Gene\" < input.vcf > out.tsv\n";
}

// --------------------------------------------------------------
// parseArguments: fill in AnnotationOptions
// --------------------------------------------------------------
static bool parseArguments(int argc, char *argv[], AnnotationOptions &opts) {
    bool showHelp = false;

    static struct option long_options[] = {
        {"annotation-extract", required_argument, 0, 'a'}, {"help", no_argument, 0, 'h'}, {0, 0, 0, 0}};

    while (true) {
        int optIdx = 0;
        int c = getopt_long(argc, argv, "a:h", long_options, &optIdx);
        if (c == -1)
            break;
        switch (c) {
        case 'a': {
            // comma-separated annotation names
            auto items = split(optarg, ',');
            for (auto &it : items) {
                // trim spaces
                while (!it.empty() && (it.front() == ' ' || it.front() == '\t')) {
                    it.erase(it.begin());
                }
                while (!it.empty() && (it.back() == ' ' || it.back() == '\t')) {
                    it.pop_back();
                }
                opts.annotations.push_back(it);
            }
        } break;
        case 'h':
        default:
            showHelp = true;
            break;
        }
    }

    if (showHelp) {
        printHelp();
        // Return false to indicate we should exit.
        return false;
    }
    // If no annotations, also show help
    if (opts.annotations.empty()) {
        printHelp();
        return false;
    }
    return true;
}

// --------------------------------------------------------------
// main extraction logic
//   1) Read VCF, pass headers unchanged
//   2) For data lines, parse ALT, parse INFO
//   3) For each ALT, line up any multi-comma annotation
//   4) Print to stdout: CHROM POS ID REF ALT [annotation1, annotation2, ...]
// --------------------------------------------------------------
static void processVCF(std::istream &in, const AnnotationOptions &opts) {
    std::string line;
    bool foundChromHeader = false;

    // Print a single header row for the TSV:
    // e.g. CHROM POS ID REF ALT ANN Gene ...
    // We'll do this only after we know how many annotation fields.
    // Or we can do it right away:
    std::cout << "CHROM\tPOS\tID\tREF\tALT";
    for (auto &annName : opts.annotations) {
        std::cout << "\t" << annName;
    }
    std::cout << "\n";

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        // If header line
        if (line[0] == '#') {
            // If you want to preserve VCF headers, you can print them.
            // But here we do NOT, since we produce a new TSV.
            // If you prefer, uncomment below:
            // std::cout << line << "\n";
            if (!foundChromHeader && line.rfind("#CHROM", 0) == 0) {
                foundChromHeader = true;
            }
            continue;
        }

        // We expect #CHROM line before data
        if (!foundChromHeader) {
            std::cerr << "Warning: Data encountered before #CHROM header: skipping\n";
            continue;
        }

        // Split the line by tabs
        // Minimal VCF => 8 fields: CHROM POS ID REF ALT QUAL FILTER INFO
        // plus optional FORMAT + samples
        auto fields = split(line, '\t');
        if (fields.size() < 8) {
            std::cerr << "Warning: Invalid VCF line (fewer than 8 fields): " << line << "\n";
            continue;
        }

        const std::string &chrom = fields[0];
        const std::string &pos = fields[1];
        const std::string &id = fields[2];
        const std::string &ref = fields[3];
        const std::string &altStr = fields[4];
        // skip fields[5]=QUAL, fields[6]=FILTER
        const std::string &info = fields[7];

        // Split ALT on comma for multiple alt alleles
        std::vector<std::string> alts = split(altStr, ',');

        // Parse INFO into a map: KEY -> value
        auto infoMap = parseInfoToMap(info);

        // For each annotation requested, retrieve the value or "NA"
        // BUT if the annotation might have multiple comma-separated sub-entries (like ANN=),
        // we handle that separately:
        // We'll do a 2-phase approach:
        //   1) gather strings for each annotation
        //   2) if it's something like "ANN" that has multiple commas, we split them
        //      and try to align them with ALT alleles.
        std::unordered_map<std::string, std::string> rawAnnValues;
        for (auto &annName : opts.annotations) {
            auto it = infoMap.find(annName);
            if (it == infoMap.end()) {
                rawAnnValues[annName] = "NA";
            } else {
                // entire string for that annotation
                rawAnnValues[annName] = it->second;
            }
        }

        // For each annotation, we split on ',' ONLY if it's marked as Number=A in the VCF header
        // Currently, we know ANN is Number=A, so we'll only split that one
        // In a more complete implementation, we would parse the VCF header to know which fields are Number=A
        std::unordered_map<std::string, std::vector<std::string>> splittedAnnValues;
        for (auto &annName : opts.annotations) {
            std::string val = rawAnnValues[annName];
            if (val == "NA") {
                // no annotation => just store empty vector
                splittedAnnValues[annName] = {};
                continue;
            }
            // Only split ANN field, which we know is Number=A
            // Other fields like Gene, Impact, DP are single-value
            if (annName == "ANN") {
                auto subVals = split(val, ',');
                splittedAnnValues[annName] = subVals;
            }
        }

        // Now let's produce lines
        for (size_t altIndex = 0; altIndex < alts.size(); ++altIndex) {
            // alt allele
            const std::string &thisAlt = alts[altIndex];

            // We'll prepare a line with columns: CHROM, POS, ID, REF, ALT, then each annotation
            // For each annotation, we see if splittedAnnValues[annName].size() > altIndex
            // If yes, output that sub-value, else "NA"
            std::ostringstream outLine;
            outLine << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << thisAlt;

            // For each requested annotation
            for (auto &annName : opts.annotations) {
                outLine << "\t";
                // Check if this is a multi-value annotation (like ANN)
                auto it = splittedAnnValues.find(annName);
                if (it != splittedAnnValues.end() && !it->second.empty()) {
                    // This is a multi-value annotation
                    if (altIndex < it->second.size()) {
                        outLine << it->second[altIndex];
                    } else {
                        outLine << "NA";
                    }
                } else {
                    // This is a single-value annotation
                    auto rawIt = rawAnnValues.find(annName);
                    if (rawIt != rawAnnValues.end()) {
                        outLine << rawIt->second;
                    } else {
                        outLine << "NA";
                    }
                }
            }
            outLine << "\n";
            std::cout << outLine.str();
        }
    }
}

// --------------------------------------------------------------
// main()
// --------------------------------------------------------------
static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_annotation_extractor", show_help))
        return 0;
    AnnotationOptions opts;
    if (!parseArguments(argc, argv, opts)) {
        // parseArguments already printed help if needed
        // Return 1 to indicate error
        return 1;
    }
    // Process from stdin => produce to stdout
    processVCF(std::cin, opts);
    return 0;
}
