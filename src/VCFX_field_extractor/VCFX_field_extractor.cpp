#include "VCFX_field_extractor.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <sstream>
#include <unordered_map>

// ------------------------------------------------------------------------
// printHelp
// ------------------------------------------------------------------------
void printHelp() {
    std::cout << "VCFX_field_extractor\n"
              << "Usage: VCFX_field_extractor --fields \"FIELD1,FIELD2,...\" [OPTIONS]\n\n"
              << "Description:\n"
              << "  Extracts specified fields from each VCF record. Fields can be:\n"
              << "    - Standard fields: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO\n"
              << "    - Subkeys in INFO (e.g. DP, AF, ANN). These are extracted from the INFO column.\n"
              << "    - Sample subfields: e.g. SampleName:GT or S2:DP, referencing the second sample's DP.\n"
              << "      You can use sample name as it appears in #CHROM line, or 'S' plus 1-based sample index.\n"
              << "If a requested field is not found or invalid, '.' is output.\n\n"
              << "Example:\n"
              << "  VCFX_field_extractor --fields \"CHROM,POS,ID,REF,ALT,DP,Sample1:GT\" < input.vcf > out.tsv\n\n"
              << "Options:\n"
              << "  --fields, -f   Comma-separated list of fields to extract\n"
              << "  --help, -h     Show this help message\n";
}

// ------------------------------------------------------------------------
// A utility function to parse "INFO" key-value pairs into a map.
// "INFO" might look like "DP=100;AF=0.5;ANN=some|stuff"
// ------------------------------------------------------------------------
static std::unordered_map<std::string, std::string> parseInfo(const std::string &infoField) {
    std::unordered_map<std::string, std::string> infoMap;
    if (infoField == "." || infoField.empty()) {
        return infoMap;
    }
    std::stringstream ss(infoField);
    std::string token;
    while (std::getline(ss, token, ';')) {
        if (token.empty()) {
            continue;
        }
        // We might have "KEY=VAL" or just "KEY" if it's a flag
        size_t eqPos = token.find('=');
        if (eqPos == std::string::npos) {
            // It's a flag
            infoMap[token] = "1";
        } else {
            std::string key = token.substr(0, eqPos);
            std::string val = token.substr(eqPos + 1);
            infoMap[key] = val;
        }
    }
    return infoMap;
}

// ------------------------------------------------------------------------
// parseLineExtract: Given a VCF data line (already split into columns),
// extracts the user-requested fields in 'fields' based on the logic:
//  - If field is one of CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO => direct from columns
//  - If it is an INFO subkey => parse info map
//  - If it is "SampleName:SUBFIELD" or "S<int>:SUBFIELD" => parse the genotype subfield
// ------------------------------------------------------------------------
static std::vector<std::string> parseLineExtract(const std::vector<std::string> &vcfCols,
                                                 const std::vector<std::string> &fields,
                                                 const std::unordered_map<std::string, int> &sampleNameToIndex) {
    // The first 8 standard columns
    // 0:CHROM,1:POS,2:ID,3:REF,4:ALT,5:QUAL,6:FILTER,7:INFO,8:FORMAT,9+:samples
    std::vector<std::string> out;
    out.reserve(fields.size());

    // parse info subkeys
    std::unordered_map<std::string, std::string> infoMap;
    if (vcfCols.size() > 7) {
        infoMap = parseInfo(vcfCols[7]); // the INFO column
    }

    // parse the FORMAT column for sample subfields
    // e.g. if format= GT:DP:GQ, then subfield "DP" is index 1
    std::vector<std::string> formatTokens;
    if (vcfCols.size() > 8) {
        std::stringstream fmts(vcfCols[8]);
        std::string fmt;
        while (std::getline(fmts, fmt, ':')) {
            formatTokens.push_back(fmt);
        }
    }

    // For each requested field
    for (auto &fld : fields) {
        std::string value = "."; // default if not found

        // Check if it's a standard field
        if (fld == "CHROM") {
            if (vcfCols.size() > 0)
                value = vcfCols[0];
        } else if (fld == "POS") {
            if (vcfCols.size() > 1)
                value = vcfCols[1];
        } else if (fld == "ID") {
            if (vcfCols.size() > 2)
                value = vcfCols[2];
        } else if (fld == "REF") {
            if (vcfCols.size() > 3)
                value = vcfCols[3];
        } else if (fld == "ALT") {
            if (vcfCols.size() > 4)
                value = vcfCols[4];
        } else if (fld == "QUAL") {
            if (vcfCols.size() > 5)
                value = vcfCols[5];
        } else if (fld == "FILTER") {
            if (vcfCols.size() > 6)
                value = vcfCols[6];
        } else if (fld == "INFO") {
            if (vcfCols.size() > 7)
                value = vcfCols[7];
        } else {
            // Possibly an INFO subkey?
            if (infoMap.find(fld) != infoMap.end()) {
                value = infoMap[fld];
            } else {
                // Possibly a sample subfield: e.g. "SampleName:GT" or "S2:DP"
                // We'll parse something like "NAME:SUBFIELD"
                // or "S<index>:SUBFIELD"
                size_t colonPos = fld.find(':');
                if (colonPos != std::string::npos) {
                    std::string sampleNameOrID = fld.substr(0, colonPos);
                    std::string subfield = fld.substr(colonPos + 1);
                    // find sample index
                    int sampleColIndex = -1;
                    if (!sampleNameOrID.empty() && sampleNameOrID[0] == 'S' &&
                        std::all_of(sampleNameOrID.begin() + 1, sampleNameOrID.end(), ::isdigit)) {
                        // format S<int>
                        int idx = std::stoi(sampleNameOrID.substr(1));
                        // sample columns start at col=9 in the VCF, but idx is 1-based
                        sampleColIndex = 9 + (idx - 1);
                    } else {
                        // sample name
                        auto itS = sampleNameToIndex.find(sampleNameOrID);
                        if (itS != sampleNameToIndex.end()) {
                            sampleColIndex = itS->second;
                        }
                    }
                    // sampleColIndex is the VCF column with that sample
                    if (sampleColIndex >= 9 && (size_t)sampleColIndex < vcfCols.size()) {
                        // parse that sample field => split by ':'
                        std::vector<std::string> sampleTokens;
                        {
                            std::stringstream sss(vcfCols[sampleColIndex]);
                            std::string tkn;
                            while (std::getline(sss, tkn, ':')) {
                                sampleTokens.push_back(tkn);
                            }
                        }
                        // find subfield in the FORMAT
                        int subIx = -1;
                        for (int i = 0; i < (int)formatTokens.size(); i++) {
                            if (formatTokens[i] == subfield) {
                                subIx = i;
                                break;
                            }
                        }
                        if (subIx >= 0 && subIx < (int)sampleTokens.size()) {
                            value = sampleTokens[subIx];
                        } else {
                            // subfield not found
                            value = ".";
                        }
                    } else {
                        // sample column not found
                        value = ".";
                    }
                } // end if colon found
            }
        } // end else standard field
        out.push_back(value);
    }

    return out;
}

// ------------------------------------------------------------------------
// extractFields:
//   1) parse #CHROM line to identify sample names => build map sample->VCFcolumnIndex
//   2) for each data line, parse columns, parse info, parse sample subfields => output
// ------------------------------------------------------------------------
void extractFields(std::istream &in, std::ostream &out, const std::vector<std::string> &fields) {
    // Print the header row (the requested field names)
    for (size_t i = 0; i < fields.size(); i++) {
        out << fields[i];
        if (i + 1 < fields.size())
            out << "\t";
    }
    out << "\n";

    std::string line;

    // We'll store sampleName -> columnIndex in this map
    std::unordered_map<std::string, int> sampleNameToIndex;
    bool foundChromHeader = false;

    // Performance: reuse vector across iterations
    std::vector<std::string> vcfCols;
    vcfCols.reserve(16);

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '#') {
            // If it's the #CHROM line, parse out sample columns
            if (!foundChromHeader && line.rfind("#CHROM", 0) == 0) {
                foundChromHeader = true;
                // parse columns using optimized split
                vcfx::split_tabs(line, vcfCols);
                // sample columns start at index=9
                for (int i = 9; i < (int)vcfCols.size(); i++) {
                    // sample name is vcfCols[i]
                    sampleNameToIndex[vcfCols[i]] = i;
                }
            }
            // skip printing header lines in the output
            continue;
        }

        // parse columns using optimized split (reuses vector capacity)
        vcfx::split_tabs(line, vcfCols);
        // parse the requested fields
        std::vector<std::string> extracted = parseLineExtract(vcfCols, fields, sampleNameToIndex);

        // print them as TSV
        for (size_t i = 0; i < extracted.size(); i++) {
            out << extracted[i];
            if (i + 1 < extracted.size())
                out << "\t";
        }
        out << "\n";
    }
}

// ------------------------------------------------------------------------
// main: parse arguments, call extractFields
// ------------------------------------------------------------------------
static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_field_extractor", show_help))
        return 0;
    std::vector<std::string> fields;
    bool showHelp = false;

    // Basic argument parsing
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            showHelp = true;
        } else if (arg.rfind("--fields", 0) == 0 || arg.rfind("-f", 0) == 0) {
            // parse next argument or after '='
            size_t eqPos = arg.find('=');
            if (eqPos != std::string::npos) {
                // e.g. --fields=CHROM,POS,ID
                std::string fieldsStr = arg.substr(eqPos + 1);
                std::stringstream sss(fieldsStr);
                std::string f;
                while (std::getline(sss, f, ',')) {
                    fields.push_back(f);
                }
            } else {
                // next argument
                if (i + 1 < argc) {
                    i++;
                    std::string fieldsStr = argv[i];
                    std::stringstream sss(fieldsStr);
                    std::string f;
                    while (std::getline(sss, f, ',')) {
                        fields.push_back(f);
                    }
                }
            }
        }
    }

    if (showHelp) {
        printHelp();
        return 0;
    }
    if (fields.empty()) {
        std::cerr << "No fields specified. Use --fields or -f to specify.\n";
        std::cerr << "Use --help for usage.\n";
        return 1;
    }

    extractFields(std::cin, std::cout, fields);
    return 0;
}
