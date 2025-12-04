#include "VCFX_probability_filter.h"
#include "vcfx_core.h" 
#include "vcfx_io.h"
#include <algorithm>
#include <getopt.h>
#include <regex>
#include <sstream>

// Implementation of VCFXProbabilityFilter
int VCFXProbabilityFilter::run(int argc, char *argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string condition;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'}, {"filter-probability", required_argument, 0, 'f'}, {0, 0, 0, 0}};

    while ((opt = getopt_long(argc, argv, "hf:", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
            break;
        case 'f':
            condition = optarg;
            break;
        default:
            showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    if (condition.empty()) {
        displayHelp();
        return 1;
    }

    // Perform probability filtering on stdin and output to stdout
    filterByProbability(std::cin, std::cout, condition);

    return 0;
}

void VCFXProbabilityFilter::displayHelp() {
    std::cout << "VCFX_probability_filter: Filter VCF based on genotype probability scores.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_probability_filter --filter-probability \"<CONDITION>\" [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                        Display this help message and exit\n";
    std::cout
        << "  -f, --filter-probability <cond>    Specify the genotype probability filter condition (e.g., GP>0.9)\n\n";
    std::cout << "Supported Operators: >, <, >=, <=, ==, !=\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_probability_filter --filter-probability \"GP>0.9\" < input.vcf > filtered.vcf\n";
}

void VCFXProbabilityFilter::filterByProbability(std::istream &in, std::ostream &out, const std::string &condition) {
    // Parse the filter condition using regex (e.g., "GP>0.9")
    std::regex conditionRegex(R"((\w+)\s*(>=|<=|>|<|==|!=)\s*([0-9]*\.?[0-9]+))");
    std::smatch matches;
    if (!std::regex_match(condition, matches, conditionRegex)) {
        std::cerr << "Error: Invalid filter condition format. Expected format like \"GP>0.9\".\n";
        return;
    }

    std::string field = matches[1];
    std::string op = matches[2];
    double threshold = std::stod(matches[3]);

    std::string line;
    bool headerParsed = false;
    std::vector<std::string> headerFields;
    size_t formatIndex = std::string::npos;
    size_t fieldIndex = std::string::npos;
    std::vector<std::string> fieldsVec;
    fieldsVec.reserve(16);

    while (std::getline(in, line)) {
        if (line.empty())
            continue;

        if (line[0] == '#') {
            // Handle header lines
            if (line.substr(0, 6) == "#CHROM") {
                // Parse header to identify the FORMAT column and sample columns
                out << line << "\n";
                headerParsed = true;
            } else {
                // Other header lines
                out << line << "\n";
            }
            continue;
        }

        if (!headerParsed) {
            std::cerr << "Error: VCF header line with #CHROM not found.\n";
            return;
        }

        // Parse VCF data lines
        vcfx::split_tabs(line, fieldsVec);

        if (fieldsVec.size() < 9) {
            std::cerr << "Warning: Invalid VCF line with fewer than 9 fields: " << line << "\n";
            continue;
        }

        std::string formatField = fieldsVec[8];
        std::vector<std::string> formatFields;
        std::stringstream fmt_ss(formatField);
        std::string fmt_field;
        while (std::getline(fmt_ss, fmt_field, ':')) {
            formatFields.push_back(fmt_field);
        }

        // Find the index of the specified field in the FORMAT column
        if (fieldIndex == std::string::npos) {
            for (size_t i = 0; i < formatFields.size(); ++i) {
                if (formatFields[i] == field) {
                    fieldIndex = i;
                    break;
                }
            }

            if (fieldIndex == std::string::npos) {
                std::cerr << "Error: Specified field \"" << field << "\" not found in FORMAT column.\n";
                return;
            }
        }

        bool pass = true;

        // Iterate over each sample to check the genotype probability
        for (size_t i = 9; i < fieldsVec.size(); ++i) {
            std::string sample = fieldsVec[i];
            std::stringstream samp_ss(sample);
            std::string samp_field;
            std::vector<std::string> sampleFields;

            while (std::getline(samp_ss, samp_field, ':')) {
                sampleFields.push_back(samp_field);
            }

            if (fieldIndex >= sampleFields.size()) {
                std::cerr << "Warning: Field index out of range in sample fields.\n";
                pass = false;
                break;
            }

            std::string valueStr = sampleFields[fieldIndex];
            if (valueStr.empty() || valueStr == ".") {
                pass = false;
                break;
            }

            double value;
            try {
                value = std::stod(valueStr);
            } catch (...) {
                std::cerr << "Warning: Unable to convert value \"" << valueStr << "\" to number.\n";
                pass = false;
                break;
            }

            // Apply the filter condition
            if (op == ">") {
                if (!(value > threshold)) {
                    pass = false;
                    break;
                }
            } else if (op == "<") {
                if (!(value < threshold)) {
                    pass = false;
                    break;
                }
            } else if (op == ">=") {
                if (!(value >= threshold)) {
                    pass = false;
                    break;
                }
            } else if (op == "<=") {
                if (!(value <= threshold)) {
                    pass = false;
                    break;
                }
            } else if (op == "==") {
                if (!(value == threshold)) {
                    pass = false;
                    break;
                }
            } else if (op == "!=") {
                if (!(value != threshold)) {
                    pass = false;
                    break;
                }
            }
        }

        if (pass) {
            out << line << "\n";
        }
    }
}

static void show_help() {
    VCFXProbabilityFilter obj;
    char arg0[] = "VCFX_probability_filter";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_probability_filter", show_help))
        return 0;
    VCFXProbabilityFilter probabilityFilter;
    return probabilityFilter.run(argc, argv);
}