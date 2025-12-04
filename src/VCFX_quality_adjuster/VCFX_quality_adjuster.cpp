#include "VCFX_quality_adjuster.h"
#include "vcfx_core.h" 
#include "vcfx_io.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <sstream>

int VCFXQualityAdjuster::run(int argc, char *argv[]) {
    if (argc == 1) {
        displayHelp();
        return 0;
    }
    bool showHelp = false;
    bool clamp = true;
    std::string transformStr;
    static struct option long_opts[] = {{"help", no_argument, 0, 'h'},
                                        {"adjust-qual", required_argument, 0, 'a'},
                                        {"no-clamp", no_argument, 0, 'n'},
                                        {0, 0, 0, 0}};
    while (true) {
        int c = ::getopt_long(argc, argv, "ha:n", long_opts, nullptr);
        if (c == -1)
            break;
        switch (c) {
        case 'h':
            showHelp = true;
            break;
        case 'a':
            transformStr = optarg;
            break;
        case 'n':
            clamp = false;
            break;
        default:
            showHelp = true;
        }
    }
    if (showHelp) {
        displayHelp();
        return 0;
    }
    if (transformStr.empty()) {
        std::cerr << "Error: Must specify a transformation with --adjust-qual <FUNC>.\n";
        displayHelp();
        return 1;
    }
    initSupportedFunctions();
    std::function<double(double)> transFunc;
    if (!parseTransformationFunction(transformStr, transFunc)) {
        std::cerr << "Error: unsupported transformation '" << transformStr << "'.\n";
        return 1;
    }
    adjustQualityScores(std::cin, std::cout, transFunc, clamp);
    return 0;
}

void VCFXQualityAdjuster::displayHelp() {
    std::cout << "VCFX_quality_adjuster: Apply a transformation to the QUAL field of a VCF.\n\n"
                 "Usage:\n"
                 "  VCFX_quality_adjuster [options] < input.vcf > output.vcf\n\n"
                 "Options:\n"
                 "  -h, --help               Show this help.\n"
                 "  -a, --adjust-qual <FUNC> Required. One of: log, sqrt, square, identity.\n"
                 "  -n, --no-clamp           Do not clamp negative or large values.\n\n"
                 "Description:\n"
                 "  Reads each line from VCF. If it's a data line with >=8 columns, we parse\n"
                 "  the QUAL field (6th col). We transform it with <FUNC>, e.g.:\n"
                 "    log => log(QUAL + 1e-10)\n"
                 "    sqrt=> sqrt(QUAL)\n"
                 "    square=> (QUAL * QUAL)\n"
                 "    identity=> no change\n"
                 "  By default, negative results from e.g. log are clamped to 0, and large\n"
                 "  results are capped at 1e12. If you do not want clamping, use --no-clamp.\n\n"
                 "Examples:\n"
                 "  1) Log-transform:\n"
                 "     VCFX_quality_adjuster --adjust-qual log < in.vcf > out.vcf\n"
                 "  2) Square, keep negative or big values as is:\n"
                 "     VCFX_quality_adjuster --adjust-qual square --no-clamp < in.vcf > out.vcf\n";
}

void VCFXQualityAdjuster::initSupportedFunctions() {
    supportedFunctions.clear();
    supportedFunctions["log"] = [](double x) {
        // protect near zero with epsilon
        return std::log(x + 1e-10);
    };
    supportedFunctions["sqrt"] = [](double x) {
        return std::sqrt(std::max(0.0, x)); // can't sqrt negative
    };
    supportedFunctions["square"] = [](double x) { return x * x; };
    supportedFunctions["identity"] = [](double x) { return x; };
}

bool VCFXQualityAdjuster::parseTransformationFunction(const std::string &funcStr,
                                                      std::function<double(double)> &transFunc) {
    auto it = supportedFunctions.find(funcStr);
    if (it == supportedFunctions.end()) {
        return false;
    }
    transFunc = it->second;
    return true;
}

void VCFXQualityAdjuster::adjustQualityScores(std::istream &in, std::ostream &out,
                                              std::function<double(double)> transFunc, bool clamp) {
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
            out << line << "\n";
            continue;
        }
        // parse fields
        vcfx::split_tabs(line, fields);
        if (fields.size() < 8) {
            std::cerr << "Warning: line with <8 fields => skipping.\n";
            continue;
        }
        double oldQual = 0.0;
        bool valid = true;
        if (fields[5] == "." || fields[5].empty()) {
            // treat missing as 0?
            oldQual = 0.0;
        } else {
            try {
                oldQual = std::stod(fields[5]);
            } catch (...) {
                std::cerr << "Warning: invalid QUAL '" << fields[5] << "'. Skipping.\n";
                valid = false;
            }
        }
        if (!valid)
            continue;
        double newQual = transFunc(oldQual);
        if (clamp) {
            if (newQual < 0.0)
                newQual = 0.0;
            // clamp large values
            if (newQual > 1e12)
                newQual = 1e12;
        }
        std::string qualStr;
        if (std::isnan(newQual)) {
            // ensure consistent representation for NaN
            qualStr = "nan";
        } else {
            qualStr = std::to_string(newQual);
        }
        fields[5] = qualStr;
        std::ostringstream oss;
        for (size_t i = 0; i < fields.size(); i++) {
            if (i > 0)
                oss << "\t";
            oss << fields[i];
        }
        out << oss.str() << "\n";
    }
}

static void show_help() {
    VCFXQualityAdjuster obj;
    char arg0[] = "VCFX_quality_adjuster";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_quality_adjuster", show_help))
        return 0;
    VCFXQualityAdjuster app;
    return app.run(argc, argv);
}