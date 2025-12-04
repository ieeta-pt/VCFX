#include "VCFX_impact_filter.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <getopt.h>
#include <iostream>
#include <regex>
#include <sstream>

// Helper function: convert string to uppercase
static std::string toUpper(const std::string &s) {
    std::string t(s);
    for (auto &c : t)
        c = std::toupper((unsigned char)c);
    return t;
}

// Simple classification of a (possibly extended) Impact value, e.g. "HIGH_SOMETHING"
enum class ImpactLevel { UNKNOWN, MODIFIER, LOW, MODERATE, HIGH };

static ImpactLevel classifyImpact(const std::string &rawImpact) {
    std::string u = toUpper(rawImpact);
    if (u.find("HIGH") != std::string::npos) {
        return ImpactLevel::HIGH;
    } else if (u.find("MODERATE") != std::string::npos) {
        return ImpactLevel::MODERATE;
    } else if (u.find("LOW") != std::string::npos) {
        return ImpactLevel::LOW;
    } else if (u.find("MODIFIER") != std::string::npos) {
        return ImpactLevel::MODIFIER;
    } else {
        return ImpactLevel::UNKNOWN;
    }
}

// We define a simple function: variantLevel >= targetLevel?
// Our hierarchy is: HIGH > MODERATE > LOW > MODIFIER > UNKNOWN
static bool meetsThreshold(ImpactLevel variantLevel, ImpactLevel targetLevel) {
    // Convert to int or we can do a switch-based approach
    auto rank = [](ImpactLevel lv) -> int {
        switch (lv) {
        case ImpactLevel::HIGH:
            return 4;
        case ImpactLevel::MODERATE:
            return 3;
        case ImpactLevel::LOW:
            return 2;
        case ImpactLevel::MODIFIER:
            return 1;
        default:
            return 0; // UNKNOWN
        }
    };
    return rank(variantLevel) >= rank(targetLevel);
}

// Implementation of VCFXImpactFilter
int VCFXImpactFilter::run(int argc, char *argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string targetImpact;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'}, {"filter-impact", required_argument, 0, 'i'}, {0, 0, 0, 0}};

    while ((opt = getopt_long(argc, argv, "hi:", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
            break;
        case 'i':
            targetImpact = optarg;
            break;
        default:
            showHelp = true;
        }
    }

    if (showHelp || targetImpact.empty()) {
        displayHelp();
        return 1;
    }

    // Filter from stdin to stdout
    filterByImpact(std::cin, std::cout, targetImpact);
    return 0;
}

void VCFXImpactFilter::displayHelp() {
    std::cout << "VCFX_impact_filter: Filter VCF variants based on predicted impact from annotations.\n\n"
              << "Usage:\n"
              << "  VCFX_impact_filter --filter-impact <LEVEL> < input.vcf > filtered.vcf\n\n"
              << "Options:\n"
              << "  -h, --help                   Show this help message\n"
              << "  -i, --filter-impact <LEVEL>  One of: HIGH, MODERATE, LOW, MODIFIER\n\n"
              << "Description:\n"
              << "  Looks in INFO for 'IMPACT=...' (case-insensitive), extracts that string,\n"
              << "  classifies it by whether it contains 'HIGH', 'MODERATE', 'LOW', or 'MODIFIER'.\n"
              << "  Then only outputs lines whose classification is >= the requested level.\n"
              << "  Also appends ';EXTRACTED_IMPACT=Value' to the INFO field.\n\n"
              << "Example:\n"
              << "  VCFX_impact_filter --filter-impact HIGH < input.vcf > filtered.vcf\n";
}

void VCFXImpactFilter::filterByImpact(std::istream &in, std::ostream &out, const std::string &targetImpact) {
    // interpret targetImpact
    ImpactLevel targetLevel = classifyImpact(targetImpact);
    if (targetLevel == ImpactLevel::UNKNOWN) {
        std::cerr << "Error: Unrecognized impact level \"" << targetImpact << "\".\n"
                  << "Must be one of HIGH, MODERATE, LOW, MODIFIER.\n";
        return;
    }

    // We'll store header lines, then insert our new meta-info line for EXTRACTED_IMPACT
    bool wroteHeader = false;
    bool wroteInfoMeta = false;

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            out << "\n";
            continue;
        }
        // If header
        if (line[0] == '#') {
            // If we see #CHROM and haven't inserted our meta line, do so
            if (!wroteInfoMeta && line.rfind("#CHROM", 0) == 0) {
                out << "##INFO=<ID=EXTRACTED_IMPACT,Number=1,Type=String,Description=\"Extracted from IMPACT=... in "
                       "info.\">\n";
                wroteInfoMeta = true;
            }
            out << line << "\n";
            if (line.rfind("#CHROM", 0) == 0) {
                wroteHeader = true;
            }
            continue;
        }

        if (!wroteHeader) {
            // if no #CHROM line yet => error
            std::cerr << "Error: VCF data encountered before #CHROM line.\n";
            return;
        }

        // parse columns
        std::stringstream ss(line);
        std::vector<std::string> fields;
        {
            std::string f;
            while (std::getline(ss, f, '\t')) {
                fields.push_back(f);
            }
        }
        if (fields.size() < 8) {
            // invalid line
            continue;
        }
        // We want to find "IMPACT=..." ignoring case
        // We'll do a simple case-insensitive search
        std::string info = fields[7];
        // We search for something like "IMPACT="
        // Then collect up to next ; or end
        // Alternatively, a case-insensitive regex. Let's do a simpler approach:
        // We'll do a small approach or a regex
        static const std::regex reImpact("IMPACT=([^;]+)", std::regex::icase);
        std::smatch m;
        std::string extracted = "UNKNOWN";
        if (std::regex_search(info, m, reImpact)) {
            // e.g. IMPACT=HIGH_Something
            extracted = m[1];
        }

        // classify
        ImpactLevel varLevel = classifyImpact(extracted);

        // if meets threshold => keep
        if (meetsThreshold(varLevel, targetLevel)) {
            // We append ;EXTRACTED_IMPACT=extracted to the info field
            // If info=="." => replace it with "EXTRACTED_IMPACT=extracted"
            if (info == "." || info.empty()) {
                info = "EXTRACTED_IMPACT=" + extracted;
            } else {
                info += ";EXTRACTED_IMPACT=" + extracted;
            }
            fields[7] = info;
            // rejoin line
            std::ostringstream joined;
            for (size_t i = 0; i < fields.size(); i++) {
                if (i > 0)
                    joined << "\t";
                joined << fields[i];
            }
            out << joined.str() << "\n";
        }
    }
}

static void show_help() {
    VCFXImpactFilter obj;
    char arg0[] = "VCFX_impact_filter";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_impact_filter", show_help))
        return 0;
    VCFXImpactFilter filt;
    return filt.run(argc, argv);
}
