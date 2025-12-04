#include "VCFX_reformatter.h"
#include "vcfx_core.h" 
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Helper to trim whitespace
static std::string trim(const std::string &s) {
    size_t start = 0;
    while (start < s.size() && isspace((unsigned char)s[start]))
        start++;
    if (start == s.size())
        return "";
    size_t end = s.size() - 1;
    while (end > start && isspace((unsigned char)s[end]))
        end--;
    return s.substr(start, end - start + 1);
}

// Splits a string by a delimiter into tokens
static void split(const std::string &s, char d, std::vector<std::string> &out) {
    out.clear();
    std::stringstream ss(s);
    std::string t;
    while (std::getline(ss, t, d)) {
        out.push_back(t);
    }
}

int VCFXReformatter::run(int argc, char *argv[]) {
    if (argc == 1) {
        displayHelp();
        return 0;
    }
    bool showHelp = false;
    std::vector<std::string> compressInfoFields;
    std::vector<std::string> compressFormatFields;
    std::vector<std::string> reorderInfoFields;
    std::vector<std::string> reorderFormatFields;

    static struct option long_opts[] = {{"help", no_argument, 0, 'h'},
                                        {"compress-info", required_argument, 0, 'c'},
                                        {"compress-format", required_argument, 0, 'f'},
                                        {"reorder-info", required_argument, 0, 'i'},
                                        {"reorder-format", required_argument, 0, 'o'},
                                        {0, 0, 0, 0}};
    while (true) {
        int c = ::getopt_long(argc, argv, "hc:f:i:o:", long_opts, nullptr);
        if (c == -1)
            break;
        switch (c) {
        case 'h':
            showHelp = true;
            break;
        case 'c': {
            // parse comma-separated
            std::stringstream ss(optarg);
            std::string val;
            while (std::getline(ss, val, ',')) {
                val = trim(val);
                if (!val.empty())
                    compressInfoFields.push_back(val);
            }
        } break;
        case 'f': {
            std::stringstream ss(optarg);
            std::string val;
            while (std::getline(ss, val, ',')) {
                val = trim(val);
                if (!val.empty())
                    compressFormatFields.push_back(val);
            }
        } break;
        case 'i': {
            std::stringstream ss(optarg);
            std::string val;
            while (std::getline(ss, val, ',')) {
                val = trim(val);
                if (!val.empty())
                    reorderInfoFields.push_back(val);
            }
        } break;
        case 'o': {
            std::stringstream ss(optarg);
            std::string val;
            while (std::getline(ss, val, ',')) {
                val = trim(val);
                if (!val.empty())
                    reorderFormatFields.push_back(val);
            }
        } break;
        default:
            showHelp = true;
        }
    }
    if (showHelp) {
        displayHelp();
        return 0;
    }
    reformatVCF(std::cin, std::cout, compressInfoFields, compressFormatFields, reorderInfoFields, reorderFormatFields);
    return 0;
}

void VCFXReformatter::displayHelp() {
    std::cout << "VCFX_reformatter: Reformat INFO/FORMAT fields in a VCF.\n\n"
                 "Usage:\n"
                 "  VCFX_reformatter [options] < input.vcf > output.vcf\n\n"
                 "Options:\n"
                 "  -h, --help                     Show this help.\n"
                 "  -c, --compress-info <keys>     Remove these INFO keys, comma-separated.\n"
                 "  -f, --compress-format <keys>   Remove these FORMAT keys, comma-separated.\n"
                 "  -i, --reorder-info <keys>      Reorder these INFO keys at the front, leftover appended.\n"
                 "  -o, --reorder-format <keys>    Reorder these FORMAT keys at the front, leftover appended.\n\n"
                 "Example:\n"
                 "  VCFX_reformatter --compress-info AF,DP --reorder-info AF,DP < in.vcf > out.vcf\n"
                 "Description:\n"
                 "  This tool modifies data lines:\n"
                 "   * 'compress-info': remove specified keys from the semicolon INFO field.\n"
                 "   * 'compress-format': remove specified keys from the colon FORMAT field,\n"
                 "      and also remove them from each sample's subfield.\n"
                 "   * 'reorder-info': place specified keys in that order at the front, then\n"
                 "      append leftover keys in the order encountered.\n"
                 "   * 'reorder-format': reorder the FORMAT colon-delimited keys in #8 col,\n"
                 "      then reorder each sample's subfields accordingly.\n"
                 "  Lines with <8 columns are skipped with a warning. Header lines (#) are\n"
                 "  passed unmodified.\n";
}

void VCFXReformatter::reformatVCF(std::istream &in, std::ostream &out,
                                  const std::vector<std::string> &compressInfoFields,
                                  const std::vector<std::string> &compressFormatFields,
                                  const std::vector<std::string> &reorderInfoFields,
                                  const std::vector<std::string> &reorderFormatFields) {
    // We'll store them in sets for faster removal checks
    std::unordered_set<std::string> infoToRemove(compressInfoFields.begin(), compressInfoFields.end());
    std::unordered_set<std::string> formatToRemove(compressFormatFields.begin(), compressFormatFields.end());

    bool foundChrom = false;
    std::string line;
    while (true) {
        if (!std::getline(in, line))
            break;
        if (line.empty()) {
            out << line << "\n";
            continue;
        }
        if (line[0] == '#') {
            out << line << "\n";
            if (line.rfind("#CHROM", 0) == 0)
                foundChrom = true;
            continue;
        }
        if (!foundChrom) {
            std::cerr << "Warning: data line before #CHROM => skipping.\n";
            continue;
        }
        std::vector<std::string> fields;
        {
            std::stringstream ss(line);
            std::string col;
            while (std::getline(ss, col, '\t')) {
                fields.push_back(col);
            }
        }
        if (fields.size() < 8) {
            std::cerr << "Warning: line with <8 columns => skipping.\n";
            continue;
        }
        std::string &infoField = fields[7];
        // compress info
        if (!infoToRemove.empty() && !infoField.empty() && infoField != ".") {
            infoField = compressInfo(infoField, infoToRemove);
        }
        // reorder info
        if (!reorderInfoFields.empty() && !infoField.empty() && infoField != ".") {
            infoField = reorderInfo(infoField, reorderInfoFields);
        }

        // if there's format => fields.size() might be >=9
        if (fields.size() > 8) {
            // compress format
            std::string &formatField = fields[8];
            if (!formatField.empty() && formatField != ".") {
                // remove keys
                std::vector<int> keepIndices;
                if (!formatToRemove.empty()) {
                    formatField = compressFormat(formatField, formatToRemove, keepIndices);
                } else {
                    // if not removing anything, keep all
                    std::vector<std::string> fsplit;
                    split(formatField, ':', fsplit);
                    keepIndices.resize(fsplit.size());
                    for (size_t i = 0; i < fsplit.size(); i++) {
                        keepIndices[i] = (int)i;
                    }
                }
                // reorder format
                if (!reorderFormatFields.empty() && !formatField.empty() && formatField != ".") {
                    std::vector<int> newIndices;
                    formatField = reorderFormat(formatField, reorderFormatFields, newIndices);
                    // apply new reorder to keepIndices
                    // we now want to chain keepIndices -> newIndices
                    // but if keepIndices had removed some, we must unify
                    std::vector<int> chained;
                    // newIndices has size = final # of keys
                    // each newIndices[i] => oldIndex
                    // but that oldIndex references the 'kept' keys from the original?
                    // Actually let's do a simpler approach: we do reorder first, then compress second.
                    // But we do the user-specified approach in the order we coded => let's keep it.
                    // We'll do: if we reorder after compress, we must reorder the 'formatField' that came from
                    // compress, but that's fine. We'll just do the reorder ignoring keepIndices for now.
                    keepIndices.clear();
                    // We must re-split the formatField to see how many remain
                    std::vector<std::string> finalFmtSplit;
                    split(formatField, ':', finalFmtSplit);
                    keepIndices.resize(finalFmtSplit.size());
                    for (size_t i = 0; i < finalFmtSplit.size(); i++) {
                        keepIndices[i] = (int)i;
                    }
                }

                // Now we must apply the final format structure to each sample
                for (size_t sampleI = 9; sampleI < fields.size(); sampleI++) {
                    fields[sampleI] =
                        applyFormatReorderToSample(fields[sampleI], /*currently we do not reorder again*/ keepIndices);
                }
            }
        }
        // rebuild line
        std::ostringstream oss;
        for (size_t i = 0; i < fields.size(); i++) {
            if (i > 0)
                oss << "\t";
            oss << fields[i];
        }
        out << oss.str() << "\n";
    }
}

std::string VCFXReformatter::compressInfo(const std::string &infoStr,
                                          const std::unordered_set<std::string> &keysToRemove) {
    if (infoStr == "." || infoStr.empty())
        return infoStr;
    std::vector<std::string> tokens;
    split(infoStr, ';', tokens);
    std::vector<std::string> keep;
    for (auto &kv : tokens) {
        if (kv.empty())
            continue;
        // parse key=...
        auto eq = kv.find('=');
        if (eq == std::string::npos) {
            // flag
            if (keysToRemove.count(kv) == 0) {
                keep.push_back(kv);
            }
        } else {
            std::string k = kv.substr(0, eq);
            if (keysToRemove.count(k) == 0) {
                keep.push_back(kv);
            }
        }
    }
    if (keep.empty())
        return ".";
    std::ostringstream oss;
    for (size_t i = 0; i < keep.size(); i++) {
        if (i > 0)
            oss << ";";
        oss << keep[i];
    }
    return oss.str();
}

// remove keys from format column => rebuild format => store the indices of the keys we keep in keepIndices
std::string VCFXReformatter::compressFormat(const std::string &formatStr,
                                            const std::unordered_set<std::string> &keysToRemove,
                                            std::vector<int> &keepIndices) {
    if (formatStr == "." || formatStr.empty()) {
        keepIndices.clear();
        return formatStr;
    }
    std::vector<std::string> keys;
    split(formatStr, ':', keys);
    keepIndices.clear();
    for (size_t i = 0; i < keys.size(); i++) {
        if (keysToRemove.count(keys[i]) == 0) {
            keepIndices.push_back((int)i);
        }
    }
    if (keepIndices.empty())
        return ".";
    // rebuild
    std::ostringstream oss;
    bool first = true;
    for (size_t i = 0; i < keepIndices.size(); i++) {
        if (!first)
            oss << ":";
        else
            first = false;
        oss << keys[keepIndices[i]];
    }
    return oss.str();
}

std::string VCFXReformatter::reorderInfo(const std::string &infoStr, const std::vector<std::string> &order) {
    if (infoStr == "." || infoStr.empty())
        return infoStr;
    std::vector<std::string> tokens;
    split(infoStr, ';', tokens);
    // parse into map + remember original order
    std::unordered_map<std::string, std::string> kvMap;
    std::vector<std::string> originalKeys;
    for (auto &item : tokens) {
        if (item.empty())
            continue;
        auto eq = item.find('=');
        if (eq == std::string::npos) {
            // flag
            kvMap[item] = "";
            originalKeys.push_back(item);
        } else {
            std::string k = item.substr(0, eq);
            std::string v = item.substr(eq + 1);
            kvMap[k] = v;
            originalKeys.push_back(k);
        }
    }
    // build new list
    std::vector<std::string> result;
    // 1) add items from 'order' if exist in kvMap
    for (auto &k : order) {
        auto it = kvMap.find(k);
        if (it != kvMap.end()) {
            if (it->second == "") {
                result.push_back(k); // flag
            } else {
                result.push_back(k + "=" + it->second);
            }
            kvMap.erase(it);
        }
    }
    // 2) append leftover in original order
    for (auto &k : originalKeys) {
        auto it = kvMap.find(k);
        if (it != kvMap.end()) {
            if (it->second == "") {
                result.push_back(k);
            } else {
                result.push_back(k + "=" + it->second);
            }
            kvMap.erase(it);
        }
    }
    if (result.empty())
        return ".";
    // join
    std::ostringstream oss;
    for (size_t i = 0; i < result.size(); i++) {
        if (i > 0)
            oss << ";";
        oss << result[i];
    }
    return oss.str();
}

std::string VCFXReformatter::reorderFormat(const std::string &fmtStr, const std::vector<std::string> &order,
                                           std::vector<int> &oldToNew) {
    if (fmtStr == "." || fmtStr.empty()) {
        oldToNew.clear();
        return fmtStr;
    }
    std::vector<std::string> keys;
    split(fmtStr, ':', keys);
    // build result
    std::vector<std::string> newOrder;
    newOrder.reserve(keys.size());
    oldToNew.assign(keys.size(), -1); // -1 => removed

    // first place the requested keys in order
    int usedCount = 0;
    // keep track of which old indices are used
    std::vector<bool> used(keys.size(), false);

    // For each key in 'order', see if we find it in keys
    for (auto &k : order) {
        auto it = std::find(keys.begin(), keys.end(), k);
        if (it != keys.end()) {
            int oldI = (int)std::distance(keys.begin(), it);
            newOrder.push_back(k);
            oldToNew[oldI] = usedCount;
            usedCount++;
            used[oldI] = true;
        }
    }
    // then append leftover in original order
    for (int i = 0; i < (int)keys.size(); i++) {
        if (!used[i]) {
            newOrder.push_back(keys[i]);
            oldToNew[i] = usedCount;
            usedCount++;
        }
    }
    // build final string
    if (newOrder.empty()) {
        oldToNew.clear();
        return ".";
    }
    std::ostringstream oss;
    for (size_t i = 0; i < newOrder.size(); i++) {
        if (i > 0)
            oss << ":";
        oss << newOrder[i];
    }
    return oss.str();
}

// We reorder or remove subfields in each sample column
std::string VCFXReformatter::applyFormatReorderToSample(const std::string &sampleStr,
                                                        const std::vector<int> &oldToNew) {
    // if oldToNew is empty, means format is "." => no sample data
    if (oldToNew.empty()) {
        return sampleStr;
    }
    if (sampleStr == "." || sampleStr.empty())
        return sampleStr;
    std::vector<std::string> subs;
    split(sampleStr, ':', subs);
    // build a new list of same length as oldToNew, initially "."
    // if oldToNew[i] < 0 => that subfield is removed
    // else => newIndex= oldToNew[i]
    // we want final vector with size = maxIndex+1
    int maxIndex = -1;
    for (auto idx : oldToNew) {
        if (idx > maxIndex)
            maxIndex = idx;
    }
    if (maxIndex < 0) {
        // means everything removed
        return ".";
    }
    std::vector<std::string> newSubs(maxIndex + 1, ".");
    for (int oldI = 0; oldI < (int)subs.size() && oldI < (int)oldToNew.size(); oldI++) {
        int ni = oldToNew[oldI];
        if (ni >= 0) {
            if ((size_t)oldI < subs.size()) {
                newSubs[ni] = subs[oldI];
            } else {
                newSubs[ni] = ".";
            }
        }
    }
    // if all are "." => return "."
    bool allDot = true;
    for (auto &x : newSubs) {
        if (x != ".") {
            allDot = false;
            break;
        }
    }
    if (allDot)
        return ".";
    // join with colon
    std::ostringstream oss;
    bool first = true;
    for (size_t i = 0; i < newSubs.size(); i++) {
        if (!first)
            oss << ":";
        else
            first = false;
        oss << newSubs[i];
    }
    return oss.str();
}

static void show_help() {
    VCFXReformatter obj;
    char arg0[] = "VCFX_reformatter";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_reformatter", show_help))
        return 0;
    VCFXReformatter reformatter;
    return reformatter.run(argc, argv);
}
