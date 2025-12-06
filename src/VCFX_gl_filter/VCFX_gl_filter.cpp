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

// Zero-allocation helper: find the N-th colon-delimited field and return pointer + length
// Returns false if fieldIndex is out of bounds or field is missing/empty
static inline bool getNthField(const char* str, size_t len, int fieldIndex,
                               const char*& fieldStart, size_t& fieldLen) {
    int currentField = 0;
    size_t pos = 0;
    size_t start = 0;

    while (pos <= len) {
        if (pos == len || str[pos] == ':') {
            if (currentField == fieldIndex) {
                fieldStart = str + start;
                fieldLen = pos - start;
                return fieldLen > 0 && !(fieldLen == 1 && *fieldStart == '.');
            }
            currentField++;
            start = pos + 1;
        }
        pos++;
    }
    return false;
}

// Zero-allocation helper: find field index by name in FORMAT string
static inline int findFieldIndex(const char* format, size_t formatLen,
                                  const char* fieldName, size_t fieldNameLen) {
    int index = 0;
    size_t pos = 0;
    size_t start = 0;

    while (pos <= formatLen) {
        if (pos == formatLen || format[pos] == ':') {
            size_t tokenLen = pos - start;
            if (tokenLen == fieldNameLen &&
                std::memcmp(format + start, fieldName, fieldNameLen) == 0) {
                return index;
            }
            index++;
            start = pos + 1;
        }
        pos++;
    }
    return -1;
}

// Fast inline integer parser for common case (GQ, DP are usually integers)
static inline bool fastParseDouble(const char* str, size_t len, double& result) {
    if (len == 0) return false;

    const char* p = str;
    const char* end = str + len;

    // Handle negative
    bool negative = false;
    if (*p == '-') {
        negative = true;
        p++;
        if (p >= end) return false;
    }

    // Parse integer part
    if (*p < '0' || *p > '9') return false;

    double intPart = 0;
    while (p < end && *p >= '0' && *p <= '9') {
        intPart = intPart * 10 + (*p - '0');
        p++;
    }

    // Parse decimal part if present
    double fracPart = 0;
    if (p < end && *p == '.') {
        p++;
        double divisor = 10;
        while (p < end && *p >= '0' && *p <= '9') {
            fracPart += (*p - '0') / divisor;
            divisor *= 10;
            p++;
        }
    }

    // For multi-value fields like PL, stop at comma
    // We accept the value even if there's more data after

    result = intPart + fracPart;
    if (negative) result = -result;
    return true;
}

void VCFXGLFilter::filterByGL(std::istream &in, std::ostream &out, const std::string &filterCondition, bool anyMode) {
    // Regex that allows field + operator + float or int
    // e.g. "GQ>20" or "DP>=3.5" etc.
    std::regex conditionRegex(R"((\w+)\s*(>=|<=|>|<|==|!=)\s*(\d+(\.\d+)?))");
    std::smatch matches;
    if (!std::regex_match(filterCondition, matches, conditionRegex)) {
        std::cerr << "Error: Invalid filter condition format. Expected e.g. \"GQ>20\" or \"DP<=3.5\".\n";
        return;
    }

    std::string field = matches[1];
    std::string op = matches[2];
    double threshold = std::stod(matches[3]);

    // Pre-compute operator type for fast comparison in hot loop
    enum class OpType { GT, LT, GE, LE, EQ, NE };
    OpType opType;
    if (op == ">") opType = OpType::GT;
    else if (op == "<") opType = OpType::LT;
    else if (op == ">=") opType = OpType::GE;
    else if (op == "<=") opType = OpType::LE;
    else if (op == "==") opType = OpType::EQ;
    else opType = OpType::NE;

    const char* fieldName = field.c_str();
    size_t fieldNameLen = field.size();

    bool headerParsed = false;
    std::string line;

    // Output buffer for better I/O performance
    std::string outputBuffer;
    outputBuffer.reserve(1024 * 1024);  // 1MB buffer

    while (std::getline(in, line)) {
        if (line.empty()) {
            outputBuffer += '\n';
            continue;
        }

        if (line[0] == '#') {
            outputBuffer += line;
            outputBuffer += '\n';
            if (line.rfind("#CHROM", 0) == 0) {
                headerParsed = true;
            }
            // Flush buffer periodically
            if (outputBuffer.size() > 900000) {
                out << outputBuffer;
                outputBuffer.clear();
            }
            continue;
        }

        if (!headerParsed) {
            std::cerr << "Error: No #CHROM header found before data.\n";
            return;
        }

        const char* linePtr = line.c_str();
        size_t lineLen = line.size();

        // Find FORMAT field (column 8) and first sample start (column 9)
        // We need to find tabs 0-9 to locate FORMAT and first sample
        int tabCount = 0;
        const char* formatStart = nullptr;
        size_t formatLen = 0;
        const char* firstSampleStart = nullptr;

        size_t fieldStart = 0;
        for (size_t i = 0; i <= lineLen; i++) {
            if (i == lineLen || linePtr[i] == '\t') {
                if (tabCount == 8) {
                    // FORMAT field
                    formatStart = linePtr + fieldStart;
                    formatLen = i - fieldStart;
                } else if (tabCount == 9) {
                    // First sample starts here
                    firstSampleStart = linePtr + fieldStart;
                    break;
                }
                tabCount++;
                fieldStart = i + 1;
            }
        }

        if (tabCount < 9 || !formatStart || !firstSampleStart) {
            std::cerr << "Warning: invalid VCF line (<9 fields)\n";
            continue;
        }

        // Find field index in FORMAT (zero allocation)
        int fieldIndex = findFieldIndex(formatStart, formatLen, fieldName, fieldNameLen);
        if (fieldIndex < 0) {
            continue;  // Field not in FORMAT, skip line
        }

        // Process samples (starting at firstSampleStart)
        bool recordPasses = anyMode ? false : true;

        // Process each sample by scanning for tabs
        const char* sampleStart = firstSampleStart;
        const char* lineEnd = linePtr + lineLen;

        while (sampleStart < lineEnd) {
            // Find end of this sample (next tab or end of line)
            const char* sampleEnd = sampleStart;
            while (sampleEnd < lineEnd && *sampleEnd != '\t') {
                sampleEnd++;
            }
            size_t sampleLen = sampleEnd - sampleStart;

            // Get the N-th field from sample (zero allocation)
            const char* valueStart;
            size_t valueLen;
            if (!getNthField(sampleStart, sampleLen, fieldIndex, valueStart, valueLen)) {
                // Missing or empty value
                if (!anyMode) {
                    recordPasses = false;
                    break;
                }
                // Move to next sample
                sampleStart = (sampleEnd < lineEnd) ? sampleEnd + 1 : lineEnd;
                continue;
            }

            // Parse value (zero allocation)
            double val;
            if (!fastParseDouble(valueStart, valueLen, val)) {
                if (!anyMode) {
                    recordPasses = false;
                    break;
                }
                // Move to next sample
                sampleStart = (sampleEnd < lineEnd) ? sampleEnd + 1 : lineEnd;
                continue;
            }

            // Compare
            bool samplePass = false;
            switch (opType) {
                case OpType::GT: samplePass = (val > threshold); break;
                case OpType::LT: samplePass = (val < threshold); break;
                case OpType::GE: samplePass = (val >= threshold); break;
                case OpType::LE: samplePass = (val <= threshold); break;
                case OpType::EQ: samplePass = (val == threshold); break;
                case OpType::NE: samplePass = (val != threshold); break;
            }

            if (anyMode) {
                if (samplePass) {
                    recordPasses = true;
                    break;
                }
            } else {
                if (!samplePass) {
                    recordPasses = false;
                    break;
                }
            }

            // Move to next sample
            sampleStart = (sampleEnd < lineEnd) ? sampleEnd + 1 : lineEnd;
        }

        if (recordPasses) {
            outputBuffer += line;
            outputBuffer += '\n';

            // Flush buffer when it gets large
            if (outputBuffer.size() > 900000) {
                out << outputBuffer;
                outputBuffer.clear();
            }
        }
    }

    // Final flush
    if (!outputBuffer.empty()) {
        out << outputBuffer;
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