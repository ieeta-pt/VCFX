#include "VCFX_sv_handler.h"
#include "vcfx_core.h" 
#include "vcfx_io.h"
#include <algorithm>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

int VCFXSvHandler::run(int argc, char *argv[]) {
    if (argc == 1) {
        displayHelp();
        return 0;
    }
    bool showHelp = false;
    bool filterOnly = false;
    bool modifySV = false;

    static struct option longOpts[] = {{"help", no_argument, 0, 'h'},
                                       {"sv-filter-only", no_argument, 0, 'f'},
                                       {"sv-modify", no_argument, 0, 'm'},
                                       {0, 0, 0, 0}};
    while (true) {
        int c = ::getopt_long(argc, argv, "hfm", longOpts, nullptr);
        if (c == -1)
            break;
        switch (c) {
        case 'h':
            showHelp = true;
            break;
        case 'f':
            filterOnly = true;
            break;
        case 'm':
            modifySV = true;
            break;
        default:
            showHelp = true;
        }
    }
    if (showHelp) {
        displayHelp();
        return 0;
    }

    handleStructuralVariants(std::cin, std::cout, filterOnly, modifySV);
    return 0;
}

void VCFXSvHandler::displayHelp() {
    std::cout << "VCFX_sv_handler: Filter or modify structural variants in a VCF.\n\n"
                 "Usage:\n"
                 "  VCFX_sv_handler [options] < input.vcf > output.vcf\n\n"
                 "Options:\n"
                 "  -h, --help           Show this help.\n"
                 "  -f, --sv-filter-only Keep only lines that have 'SVTYPE=' in their INFO.\n"
                 "  -m, --sv-modify      Modify the INFO field of structural variants.\n\n"
                 "Description:\n"
                 "  * If --sv-filter-only is set, we skip lines without structural variant.\n"
                 "  * If --sv-modify is set, we add 'SV_VALIDATED=1', 'SV_SIZE=...' for DEL/DUP.\n"
                 "    Also 'INV_TYPE=PARALLEL' for INV, 'BND_ORIENTATION=PAIR' for BND. etc.\n"
                 "  * If both are set, we do both filtering and modification.\n"
                 "  * Non-SV lines are only included if !filterOnly.\n\n"
                 "Example:\n"
                 "  1) Keep only structural variants:\n"
                 "     VCFX_sv_handler --sv-filter-only < in.vcf > out.vcf\n"
                 "  2) Modify structural variants:\n"
                 "     VCFX_sv_handler --sv-modify < in.vcf > out.vcf\n"
                 "  3) Do both:\n"
                 "     VCFX_sv_handler --sv-filter-only --sv-modify < in.vcf > out.vcf\n";
}

bool VCFXSvHandler::isStructuralVariant(const std::string &infoField) const {
    return infoField.find("SVTYPE=") != std::string::npos;
}

std::string VCFXSvHandler::parseSVType(const std::string &infoField) const {
    // find "SVTYPE="
    auto pos = infoField.find("SVTYPE=");
    if (pos == std::string::npos)
        return "";
    size_t start = pos + 7;
    // read until next ; or end
    auto semicol = infoField.find(';', start);
    if (semicol == std::string::npos) {
        return infoField.substr(start);
    } else {
        return infoField.substr(start, semicol - start);
    }
}

int VCFXSvHandler::parseEndPosition(const std::string &infoField) const {
    auto pos = infoField.find("END=");
    if (pos == std::string::npos)
        return -1;
    size_t start = pos + 4;
    auto semicol = infoField.find(';', start);
    std::string endStr =
        (semicol == std::string::npos) ? infoField.substr(start) : infoField.substr(start, semicol - start);
    try {
        return std::stoi(endStr);
    } catch (...) {
        return -1;
    }
}

int VCFXSvHandler::parsePos(const std::string &posField) const {
    try {
        return std::stoi(posField);
    } catch (...) {
        return -1;
    }
}

std::string VCFXSvHandler::manipulateSVInfo(const std::string &infoField, const std::string &svType, int pos,
                                            int endPos) const {
    std::string modified = infoField;
    if (!modified.empty() && modified.back() != ';') {
        modified += ";";
    }
    modified += "SV_VALIDATED=1"; // a sample annotation

    // if we have endPos>0 and pos>0 => maybe we compute size
    if (endPos > 0 && pos > 0 && (svType == "DEL" || svType == "DUP")) {
        int svSize = endPos - pos;
        // possibly handle negative
        if (svSize < 0)
            svSize = -svSize;
        modified += ";SV_SIZE=" + std::to_string(svSize);
    }
    if (svType == "INV") {
        modified += ";INV_TYPE=PARALLEL";
    } else if (svType == "BND") {
        modified += ";BND_ORIENTATION=PAIR";
    }
    // other types?
    return modified;
}

void VCFXSvHandler::handleStructuralVariants(std::istream &in, std::ostream &out, bool filterOnly, bool modifySV) {
    std::string line;
    // Declare fields vector outside loop for reuse
    std::vector<std::string> fields;
    fields.reserve(16);

    while (true) {
        if (!std::getline(in, line))
            break;
        if (line.empty()) {
            continue;
        }
        if (line[0] == '#') {
            out << line << "\n";
            continue;
        }
        // parse
        vcfx::split_tabs(line, fields);
        if (fields.size() < 8) {
            std::cerr << "Warning: skipping line with <8 columns.\n";
            continue;
        }
        bool isSV = isStructuralVariant(fields[7]);
        if (isSV) {
            // if filterOnly && !modify => just output
            // if modify => we parse and rewrite info
            if (filterOnly && !modifySV) {
                out << line << "\n";
                continue;
            }
            if (modifySV) {
                // parse svtype
                std::string svType = parseSVType(fields[7]);
                if (svType.empty()) {
                    std::cerr << "Warning: no SVTYPE => skipping line.\n";
                    continue;
                }
                int pos = parsePos(fields[1]);
                int endPos = parseEndPosition(fields[7]);
                if (pos < 0) {
                    std::cerr << "Warning: invalid POS => skipping.\n";
                    continue;
                }
                std::string newInfo = manipulateSVInfo(fields[7], svType, pos, endPos);
                fields[7] = newInfo;
                // rebuild
                bool first = true;
                for (size_t i = 0; i < fields.size(); i++) {
                    if (!first)
                        out << "\t";
                    else
                        first = false;
                    out << fields[i];
                }
                out << "\n";
                continue;
            }
            // if we get here => isSV, but not filterOnly or modify => just output
            out << line << "\n";
        } else {
            // not an SV
            if (!filterOnly) {
                // keep it
                out << line << "\n";
            }
            // else skip
        }
    }
}

static void show_help() {
    VCFXSvHandler obj;
    char arg0[] = "VCFX_sv_handler";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_sv_handler", show_help))
        return 0;
    VCFXSvHandler app;
    return app.run(argc, argv);
}
