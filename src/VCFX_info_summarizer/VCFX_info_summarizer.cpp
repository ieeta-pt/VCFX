#include "VCFX_info_summarizer.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <map>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

// Function to display help message
void printHelp() {
    std::cout << "VCFX_info_summarizer\n"
              << "Usage: VCFX_info_summarizer [OPTIONS]\n\n"
              << "Options:\n"
              << "  --info, -i \"FIELD1,FIELD2\"   Specify the INFO fields to summarize (e.g., \"DP,AF\").\n"
              << "  --help, -h                    Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Summarizes numeric fields in the INFO column of a VCF file by calculating\n"
              << "  statistics such as mean, median, and mode.\n\n"
              << "Examples:\n"
              << "  ./VCFX_info_summarizer --info \"DP,AF\" < input.vcf > summary_stats.tsv\n";
}

// Function to parse command-line arguments
bool parseArguments(int argc, char *argv[], std::vector<std::string> &info_fields) {
    bool foundAnyField = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        // Check for --info or -i with the next argument
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
        }
        // Check for --info=<FIELDS>
        else if (arg.rfind("--info=", 0) == 0) {
            std::string fields_str = arg.substr(7);
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
        } else if (arg == "--help" || arg == "-h") {
            printHelp();
            std::exit(0);
        }
        // else ignore unrecognized options
    }

    if (!foundAnyField) {
        std::cerr << "Error: INFO fields not specified.\n";
        std::cerr << "Use --help for usage information.\n";
        return false;
    }
    return true;
}

// Function to calculate mean
double calculateMean(const std::vector<double> &data) {
    if (data.empty())
        return 0.0;
    double sum = 0.0;
    for (auto val : data) {
        sum += val;
    }
    return sum / static_cast<double>(data.size());
}

// Function to calculate median
double calculateMedian(std::vector<double> data) {
    if (data.empty())
        return 0.0;
    std::sort(data.begin(), data.end());
    size_t n = data.size();
    if (n % 2 == 0) {
        return (data[n / 2 - 1] + data[n / 2]) / 2.0;
    } else {
        return data[n / 2];
    }
}

// Function to calculate mode
double calculateMode(const std::vector<double> &data) {
    if (data.empty())
        return 0.0;
    std::unordered_map<double, int> frequency;
    int maxFreq = 0;
    double modeValue = data[0];

    for (auto val : data) {
        frequency[val]++;
        if (frequency[val] > maxFreq) {
            maxFreq = frequency[val];
            modeValue = val;
        }
    }
    return modeValue;
}

// Function to parse the INFO field and collect specified fields
bool summarizeInfoFields(std::istream &in, std::ostream &out, const std::vector<std::string> &info_fields) {
    std::string line;
    bool header_found = false;

    // Map to store vectors of values for each requested INFO field
    std::map<std::string, std::vector<double>> info_data;
    for (const auto &field : info_fields) {
        info_data[field]; // ensures key is created
    }

    // Performance: reuse containers across iterations
    std::vector<std::string> fields;
    fields.reserve(16);

    while (std::getline(in, line)) {
        if (line.empty())
            continue;

        if (line[0] == '#') {
            // Check if it is #CHROM
            if (line.rfind("#CHROM", 0) == 0) {
                header_found = true;
            }
            continue;
        }

        if (!header_found) {
            std::cerr << "Error: VCF header (#CHROM) not found before records.\n";
            return false;
        }

        // parse columns
        vcfx::split_tabs(line, fields);
        if (fields.size() < 8) {
            std::cerr << "Warning: Skipping malformed VCF line: " << line << "\n";
            continue;
        }

        // Parse INFO field => store key->value
        std::unordered_map<std::string, std::string> info_map;
        {
            std::stringstream info_ss(fields[7]);
            std::string kv;
            while (std::getline(info_ss, kv, ';')) {
                if (kv.empty())
                    continue;
                size_t eq = kv.find('=');
                if (eq != std::string::npos) {
                    std::string key = kv.substr(0, eq);
                    std::string value = kv.substr(eq + 1);
                    info_map[key] = value;
                } else {
                    // A flag => store "1"
                    info_map[kv] = "1";
                }
            }
        }

        // For each requested field
        for (const auto &field : info_fields) {
            auto it = info_map.find(field);
            if (it == info_map.end()) {
                // not present
                continue;
            }
            // If present, possibly multiple comma-separated values
            std::stringstream valSS(it->second);
            std::string val;
            while (std::getline(valSS, val, ',')) {
                try {
                    double v = std::stod(val);
                    // Skip NaN / Inf
                    if (std::isnan(v) || std::isinf(v)) {
                        std::cerr << "Warning: Non-finite value for field " << field << " in line: " << line << "\n";
                        continue; // skip
                    }
                    info_data[field].push_back(v);
                } catch (...) {
                    std::cerr << "Warning: Non-numeric value for field " << field << " in line: " << line << "\n";
                }
            }
        }
    }

    // Print summary table
    out << "INFO_Field\tMean\tMedian\tMode\n";
    for (const auto &field : info_fields) {
        const auto &data = info_data.at(field);
        if (data.empty()) {
            // no numeric data => NA
            out << field << "\tNA\tNA\tNA\n";
            continue;
        }
        double mean = calculateMean(data);
        double median = calculateMedian(data);
        double mode = calculateMode(data);
        out << field << "\t" << std::fixed << std::setprecision(4) << mean << "\t" << median << "\t" << mode << "\n";
    }

    return true;
}

static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_info_summarizer", show_help))
        return 0;
    std::vector<std::string> info_fields;

    // parse arguments
    if (!parseArguments(argc, argv, info_fields)) {
        // prints error: "Error: INFO fields not specified"
        return 1;
    }

    // Summarize INFO fields
    bool success = summarizeInfoFields(std::cin, std::cout, info_fields);
    return success ? 0 : 1;
}
