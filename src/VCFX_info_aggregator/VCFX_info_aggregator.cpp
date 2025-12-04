#include "VCFX_info_aggregator.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cmath> // for std::isfinite
#include <cstdlib>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

// ----------------------------------------------------------------------
// Displays help
// ----------------------------------------------------------------------
void VCFXInfoAggregator::displayHelp() {
    std::cout << "VCFX_info_aggregator: Aggregate numeric INFO field values from a VCF.\n\n"
              << "Usage:\n"
              << "  VCFX_info_aggregator --aggregate-info \"DP,AF,...\" < input.vcf > output.vcf\n\n"
              << "Description:\n"
              << "  Reads a VCF from stdin, prints it unmodified, and at the end, appends a\n"
              << "  summary section of the form:\n"
              << "    #AGGREGATION_SUMMARY\n"
              << "    DP: Sum=..., Average=...\n"
              << "    AF: Sum=..., Average=...\n"
              << "  The VCF portion remains fully valid. The final lines start with '#' so most\n"
              << "  VCF parsers will ignore them.\n\n"
              << "Options:\n"
              << "  -h, --help                Print this help message.\n"
              << "  -a, --aggregate-info <fields>  Comma-separated list of INFO fields to aggregate.\n\n"
              << "Example:\n"
              << "  VCFX_info_aggregator --aggregate-info \"DP,AF\" < input.vcf > aggregated.vcf\n";
}

// ----------------------------------------------------------------------
// parse command line, call aggregator
// ----------------------------------------------------------------------
int VCFXInfoAggregator::run(int argc, char *argv[]) {
    int opt;
    bool showHelp = false;
    std::string infoFieldsStr;

    static struct option longOptions[] = {
        {"help", no_argument, 0, 'h'}, {"aggregate-info", required_argument, 0, 'a'}, {0, 0, 0, 0}};

    while ((opt = getopt_long(argc, argv, "ha:", longOptions, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
            break;
        case 'a':
            infoFieldsStr = optarg;
            break;
        default:
            showHelp = true;
            break;
        }
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }
    if (infoFieldsStr.empty()) {
        std::cerr << "Error: Must specify --aggregate-info with at least one field.\n";
        return 1;
    }

    // parse comma separated fields
    std::vector<std::string> infoFields;
    {
        std::stringstream ss(infoFieldsStr);
        std::string f;
        while (std::getline(ss, f, ',')) {
            // trim
            while (!f.empty() && isspace((unsigned char)f.back()))
                f.pop_back();
            while (!f.empty() && isspace((unsigned char)f.front()))
                f.erase(f.begin());
            if (!f.empty()) {
                infoFields.push_back(f);
            }
        }
    }
    if (infoFields.empty()) {
        std::cerr << "Error: no valid fields in --aggregate-info\n";
        return 1;
    }

    aggregateInfo(std::cin, std::cout, infoFields);
    return 0;
}

// ----------------------------------------------------------------------
// aggregator function
//   1) read lines
//   2) if line is header or data, print unmodified to out
//   3) parse data lines => parse 'INFO' col => parse each requested field => if numeric => accumulate
//   4) after reading entire file, print summary
// ----------------------------------------------------------------------
void VCFXInfoAggregator::aggregateInfo(std::istream &in, std::ostream &out,
                                       const std::vector<std::string> &infoFields) {
    // We'll store each field's collected numeric values in a vector
    std::map<std::string, std::vector<double>> collected;
    for (auto &fld : infoFields) {
        collected[fld] = {};
    }

    bool foundChromHeader = false;

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
            // pass it through unmodified
            out << line << "\n";
            if (!foundChromHeader && line.rfind("#CHROM", 0) == 0) {
                foundChromHeader = true;
            }
            continue;
        }

        if (!foundChromHeader) {
            std::cerr << "Error: encountered data line before #CHROM header.\n";
            break; // break or continue, but we'll break
        }

        // parse columns
        // minimal vcf => CHROM POS ID REF ALT QUAL FILTER INFO ...
        // We only need the 8th col => INFO
        vcfx::split_tabs(line, fields);
        // pass line unmodified
        out << line << "\n";

        if (fields.size() < 8) {
            // skip aggregator for this line
            continue;
        }
        std::string &info = fields[7];

        // parse info col => find each "key=val" or "key" flag
        // we only look for key=val with numeric val if the key is in infoFields
        // We'll do a naive parse: split by ';', then split each part at '='
        {
            std::stringstream infoSS(info);
            std::string item;
            while (std::getline(infoSS, item, ';')) {
                if (item.empty())
                    continue;
                // find '='
                size_t eqPos = item.find('=');
                if (eqPos == std::string::npos) {
                    // no '=' => skip
                    continue;
                }
                std::string key = item.substr(0, eqPos);
                std::string val = item.substr(eqPos + 1);
                if (key.empty() || val.empty())
                    continue;
                // trim spaces
                while (!key.empty() && isspace((unsigned char)key.back()))
                    key.pop_back();
                while (!key.empty() && isspace((unsigned char)key.front()))
                    key.erase(key.begin());
                while (!val.empty() && isspace((unsigned char)val.back()))
                    val.pop_back();
                while (!val.empty() && isspace((unsigned char)val.front()))
                    val.erase(val.begin());

                auto it = collected.find(key);
                if (it == collected.end()) {
                    // not one of requested fields
                    continue;
                }
                // parse numeric, skip if NaN or Inf
                try {
                    double d = std::stod(val);
                    // check for finite
                    if (std::isfinite(d)) {
                        it->second.push_back(d);
                    }
                    // else skip
                } catch (...) {
                    // skip
                }
            }
        }
    }

    // now print aggregator summary
    out << "#AGGREGATION_SUMMARY\n";
    for (auto &kv : collected) {
        const std::string &field = kv.first;
        const auto &vals = kv.second;
        double sum = 0.0;
        for (auto &v : vals) {
            sum += v;
        }
        double mean = 0.0;
        if (!vals.empty()) {
            mean = sum / static_cast<double>(vals.size());
        }
        out << field << ": Sum=" << sum << ", Average=" << mean << "\n";
    }
}

static void show_help() {
    VCFXInfoAggregator obj;
    char arg0[] = "VCFX_info_aggregator";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_info_aggregator", show_help))
        return 0;
    VCFXInfoAggregator app;
    return app.run(argc, argv);
}
