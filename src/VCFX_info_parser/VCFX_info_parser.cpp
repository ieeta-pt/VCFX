#include "VCFX_info_parser.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cstdlib>
#include <sstream>
#include <unordered_map>

// Function to display help message
void printHelp() {
    std::cout
        << "VCFX_info_parser\n"
        << "Usage: VCFX_info_parser [OPTIONS]\n\n"
        << "Options:\n"
        << "  --info, -i \"FIELD1,FIELD2\"   Specify the INFO fields to display (e.g., \"DP,AF\").\n"
        << "  --help, -h                    Display this help message and exit.\n\n"
        << "Description:\n"
        << "  Parses the INFO field of a VCF file and displays the selected INFO fields in a user-friendly format.\n\n"
        << "Examples:\n"
        << "  ./VCFX_info_parser --info \"DP,AF\" < input.vcf > output_info.tsv\n";
}

// Function to parse command-line arguments
bool parseArguments(int argc, char *argv[], std::vector<std::string> &info_fields) {
    bool foundAnyField = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if ((arg == "--info" || arg == "-i") && i + 1 < argc) {
            std::string fields_str = argv[++i];
            std::stringstream ss(fields_str);
            std::string field;
            while (std::getline(ss, field, ',')) {
                // Trim whitespace
                field.erase(0, field.find_first_not_of(" \t\n\r\f\v"));
                field.erase(field.find_last_not_of(" \t\n\r\f\v") + 1);
                if (!field.empty()) {
                    info_fields.push_back(field);
                    foundAnyField = true;
                }
            }
        } else if (arg.rfind("--info=", 0) == 0) {
            // e.g. --info=DP,AF
            std::string fields_str = arg.substr(7);
            std::stringstream ss(fields_str);
            std::string field;
            while (std::getline(ss, field, ',')) {
                // Trim
                field.erase(0, field.find_first_not_of(" \t\n\r\f\v"));
                field.erase(field.find_last_not_of(" \t\n\r\f\v") + 1);
                if (!field.empty()) {
                    info_fields.push_back(field);
                    foundAnyField = true;
                }
            }
        } else if (arg == "--help" || arg == "-h") {
            printHelp();
            std::exit(0);
        }
        // ignore other unrecognized args...
    }

    return foundAnyField;
}

// Splits a string by a delimiter (optimized with find-based approach)
std::vector<std::string> split(const std::string &s, char delimiter) {
    std::vector<std::string> tokens;
    tokens.reserve(8);
    size_t start = 0;
    size_t end;
    while ((end = s.find(delimiter, start)) != std::string::npos) {
        tokens.emplace_back(s, start, end - start);
        start = end + 1;
    }
    tokens.emplace_back(s, start);
    return tokens;
}

// Parses the INFO field and display selected fields
bool parseInfoFields(std::istream &in, std::ostream &out, const std::vector<std::string> &info_fields) {
    // If we have any fields, print header
    if (!info_fields.empty()) {
        out << "CHROM\tPOS\tID\tREF\tALT";
        for (const auto &field : info_fields) {
            out << "\t" << field;
        }
        out << "\n";
    }

    std::string line;

    // Performance: reuse containers across iterations
    std::vector<std::string> fields;
    fields.reserve(16);
    std::unordered_map<std::string, std::string> info_map;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '#') {
            // skip header lines
            continue;
        }

        // Performance: use find-based tab splitting instead of stringstream
        vcfx::split_tabs(line, fields);
        if (fields.size() < 8) {
            std::cerr << "Warning: Skipping invalid VCF line.\n";
            continue;
        }

        // parse INFO => key-value or flag
        info_map.clear();
        auto info_entries = split(fields[7], ';');
        for (auto &entry : info_entries) {
            size_t eqPos = entry.find('=');
            if (eqPos != std::string::npos) {
                std::string key = entry.substr(0, eqPos);
                std::string value = entry.substr(eqPos + 1);
                info_map[key] = value;
            } else {
                // This is a flag (no '=' => empty string)
                info_map[entry] = "";
            }
        }

        // Print row
        out << fields[0] << "\t" << fields[1] << "\t" << fields[2] << "\t" << fields[3] << "\t" << fields[4];
        for (const auto &field : info_fields) {
            auto it = info_map.find(field);
            if (it != info_map.end()) {
                // If it's a flag (empty string), print "."
                if (it->second.empty()) {
                    out << "\t.";
                } else {
                    out << "\t" << it->second;
                }
            } else {
                // Field not present
                out << "\t.";
            }
        }
        out << "\n";
    }
    return true;
}

static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_info_parser", show_help))
        return 0;
    std::vector<std::string> info_fields;

    // parse arguments
    bool hasFields = parseArguments(argc, argv, info_fields);
    if (!hasFields) {
        std::cerr << "Error: INFO fields not specified.\n"
                  << "Use --help for usage information.\n";
        return 1;
    }

    // parse and display
    bool ok = parseInfoFields(std::cin, std::cout, info_fields);
    return ok ? 0 : 1;
}
