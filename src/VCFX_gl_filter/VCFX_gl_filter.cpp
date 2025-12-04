#include "VCFX_gl_filter.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

// Implementation of VCFXGLFilter
int VCFXGLFilter::run(int argc, char *argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string filterCondition;
    bool anyMode = false; // default is 'all' mode
    bool invalidMode = false;

    static struct option long_options[] = {{"help", no_argument, 0, 'h'},
                                           {"filter", required_argument, 0, 'f'},
                                           {"mode", required_argument, 0, 'm'},
                                           {0, 0, 0, 0}};

    while ((opt = getopt_long(argc, argv, "hf:m:", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
            break;
        case 'f':
            filterCondition = optarg;
            break;
        case 'm':
            if (std::string(optarg) == "any") {
                anyMode = true;
            } else if (std::string(optarg) == "all") {
                anyMode = false;
            } else {
                std::cerr << "Error: --mode must be 'any' or 'all'.\n";
                invalidMode = true;
                showHelp = true;
            }
            break;
        default:
            showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 1;
    }

    if (filterCondition.empty()) {
        std::cerr << "Error: --filter must be specified.\n";
        displayHelp();
        return 1;
    }

    // Redirect cerr to a string stream to check for errors
    std::stringstream error_output;
    std::streambuf *cerr_buffer = std::cerr.rdbuf();
    std::cerr.rdbuf(error_output.rdbuf());

    // Filter VCF based on genotype likelihood
    filterByGL(std::cin, std::cout, filterCondition, anyMode);

    // Restore cerr
    std::cerr.rdbuf(cerr_buffer);

    // Check if there were any errors
    std::string error_msg = error_output.str();
    if (!error_msg.empty()) {
        // There was an error, output it and return error code
        std::cerr << error_msg;
        return 1;
    }

    return 0;
}

void VCFXGLFilter::displayHelp() {
    std::cout << "VCFX_gl_filter: Filter VCF based on a numeric genotype-likelihood field.\n\n"
              << "Usage:\n"
              << "  VCFX_gl_filter --filter \"<CONDITION>\" [--mode <any|all>] < input.vcf > output.vcf\n\n"
              << "Options:\n"
              << "  -h, --help                Display this help message and exit\n"
              << "  -f, --filter <CONDITION>  e.g. \"GQ>20\" or \"DP>=10.5\" or \"PL==50\"\n"
              << "  -m, --mode <any|all>      'all' => all samples must pass (default), 'any' => at least one sample "
                 "passes.\n\n"
              << "Example:\n"
              << "  VCFX_gl_filter --filter \"GQ>20.5\" --mode any < input.vcf > filtered.vcf\n\n"
              << "Description:\n"
              << "  The filter condition is a simple expression: <Field><op><value>,\n"
              << "  e.g. GQ>20 or DP!=10 or RGQ<=5.2.\n"
              << "  The 'mode' determines if all samples must satisfy the condition or\n"
              << "  if at least one sample satisfying is enough to keep the record.\n";
}

void VCFXGLFilter::filterByGL(std::istream &in, std::ostream &out, const std::string &filterCondition, bool anyMode) {
    // Regex that allows field + operator + float or int
    // e.g. "GQ>20" or "DP>=3.5" etc.
    std::regex conditionRegex(R"((\w+)\s*(>=|<=|>|<|==|!=)\s*(\d+(\.\d+)?))");
    std::smatch matches;
    if (!std::regex_match(filterCondition, matches, conditionRegex)) {
        std::cerr << "Error: Invalid filter condition format. Expected e.g. \"GQ>20\" or \"DP<=3.5\".\n";
        // Don't continue processing, just return
        return; // The error has been reported, the caller will handle the return code
    }

    std::string field = matches[1];
    std::string op = matches[2];
    double threshold = std::stod(matches[3]);

    bool headerParsed = false;
    std::vector<std::string> headerFields;
    headerFields.reserve(16);
    std::string line;

    // Reuse vectors across iterations
    std::vector<std::string> fieldsVec;
    fieldsVec.reserve(16);
    std::vector<std::string> formatTokens;
    formatTokens.reserve(8);
    std::vector<std::string> sampleTokens;
    sampleTokens.reserve(8);

    while (true) {
        if (!std::getline(in, line)) {
            break; // EOF
        }
        if (line.empty()) {
            out << "\n";
            continue;
        }

        if (line[0] == '#') {
            // Output header lines as is
            out << line << "\n";
            if (line.rfind("#CHROM", 0) == 0) {
                headerParsed = true;
            }
            continue;
        }

        // Must have #CHROM
        if (!headerParsed) {
            std::cerr << "Error: No #CHROM header found before data.\n";
            return; // Just return, the caller will handle the error status
        }

        // Parse data line
        vcfx::split_tabs(line, fieldsVec);
        if (fieldsVec.size() < 9) {
            std::cerr << "Warning: invalid VCF line (<9 fields): " << line << "\n";
            continue;
        }

        // parse the FORMAT field
        const std::string &formatStr = fieldsVec[8];
        formatTokens.clear();
        {
            size_t start = 0;
            size_t end;
            while ((end = formatStr.find(':', start)) != std::string::npos) {
                formatTokens.emplace_back(formatStr, start, end - start);
                start = end + 1;
            }
            formatTokens.emplace_back(formatStr, start);
        }

        // find field in format
        int fieldIndex = -1;
        for (int i = 0; i < (int)formatTokens.size(); i++) {
            if (formatTokens[i] == field) {
                fieldIndex = i;
                break;
            }
        }
        if (fieldIndex < 0) {
            // If not found, do we keep or skip? We skip or keep. Let's skip
            // or we might keep, but let's skip by default
            // or we can pass if user wants
            // We'll just skip
            // But let's do a note:
            // out << line << "\n"; // or skip
            continue;
        }

        // We'll check sample columns from index=9 onward
        bool recordPasses = anyMode ? false : true; // in 'any' mode we set false => flip if one sample passes
                                                    // in 'all' mode we set true => flip to false if one sample fails

        for (size_t s = 9; s < fieldsVec.size(); s++) {
            // parse sample by ':'
            const std::string &sampleStr = fieldsVec[s];
            sampleTokens.clear();
            {
                size_t start = 0;
                size_t end;
                while ((end = sampleStr.find(':', start)) != std::string::npos) {
                    sampleTokens.emplace_back(sampleStr, start, end - start);
                    start = end + 1;
                }
                sampleTokens.emplace_back(sampleStr, start);
            }

            if (fieldIndex >= (int)sampleTokens.size()) {
                // no data => fail for all mode, or do nothing for any mode
                if (!anyMode) {
                    recordPasses = false;
                }
                break;
            }
            const std::string &valStr = sampleTokens[fieldIndex];
            if (valStr.empty() || valStr == ".") {
                // treat as fail for all mode
                if (!anyMode) {
                    recordPasses = false;
                }
                break;
            }
            double val = 0.0;
            try {
                val = std::stod(valStr);
            } catch (...) {
                // not numeric => fail for all mode
                if (!anyMode) {
                    recordPasses = false;
                }
                break;
            }
            // compare
            bool samplePass = false;
            if (op == ">") {
                samplePass = (val > threshold);
            } else if (op == "<") {
                samplePass = (val < threshold);
            } else if (op == ">=") {
                samplePass = (val >= threshold);
            } else if (op == "<=") {
                samplePass = (val <= threshold);
            } else if (op == "==") {
                samplePass = (val == threshold);
            } else if (op == "!=") {
                samplePass = (val != threshold);
            }

            if (anyMode) {
                // if one sample passes => keep record => break
                if (samplePass) {
                    recordPasses = true;
                    break;
                }
            } else {
                // all mode => if one fails => skip record => break
                if (!samplePass) {
                    recordPasses = false;
                    break;
                }
            }
        }

        if (recordPasses) {
            out << line << "\n";
        }
    }
}

static void show_help() {
    VCFXGLFilter obj;
    char arg0[] = "VCFX_gl_filter";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_gl_filter", show_help))
        return 0;
    VCFXGLFilter app;
    return app.run(argc, argv);
}