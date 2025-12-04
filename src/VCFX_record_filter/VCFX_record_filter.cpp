#include "VCFX_record_filter.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

// Trim leading/trailing whitespace
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

static void split(const std::string &s, char delim, std::vector<std::string> &out) {
    out.clear();
    size_t start = 0;
    size_t end;
    while ((end = s.find(delim, start)) != std::string::npos) {
        out.emplace_back(s, start, end - start);
        start = end + 1;
    }
    out.emplace_back(s, start);
}

// parse an operator token: "==", "!=", ">=", "<=", ">", "<"
static bool parseOperator(const std::string &opStr, FilterOp &op) {
    if (opStr == ">") {
        op = FilterOp::GT;
        return true;
    }
    if (opStr == ">=") {
        op = FilterOp::GE;
        return true;
    }
    if (opStr == "<") {
        op = FilterOp::LT;
        return true;
    }
    if (opStr == "<=") {
        op = FilterOp::LE;
        return true;
    }
    if (opStr == "==") {
        op = FilterOp::EQ;
        return true;
    }
    if (opStr == "!=") {
        op = FilterOp::NE;
        return true;
    }
    return false;
}

// parse single token e.g. "POS>100" or "FILTER!=PASS"
// returns false if invalid
static bool parseSingleCriterion(const std::string &token, FilterCriterion &crit) {
    // find the operator among [>=,<=,==,!=,>,<]
    // a robust approach is to search for them in descending length order
    static std::vector<std::string> ops = {">=", "<=", "==", "!=", ">", "<"};
    size_t bestPos = std::string::npos;
    std::string bestOp;
    for (auto &o : ops) {
        auto p = token.find(o);
        if (p != std::string::npos) {
            bestPos = p;
            bestOp = o;
            break; // first match
        }
    }
    if (bestPos == std::string::npos) {
        std::cerr << "Error: no operator found in '" << token << "'.\n";
        return false;
    }
    // parse fieldName= token.substr(0,bestPos)
    // parse value= token.substr(bestPos+ bestOp.size())
    crit.fieldName = trim(token.substr(0, bestPos));
    if (crit.fieldName.empty()) {
        std::cerr << "Error: empty field name in '" << token << "'.\n";
        return false;
    }
    if (!parseOperator(bestOp, crit.op)) {
        std::cerr << "Error: unknown operator '" << bestOp << "' in '" << token << "'.\n";
        return false;
    }
    std::string valPart = trim(token.substr(bestPos + bestOp.size()));
    if (valPart.empty()) {
        std::cerr << "Error: no value in '" << token << "'.\n";
        return false;
    }
    // try parse numeric
    // if parse fails => treat as string
    try {
        double d = std::stod(valPart);
        // numeric
        crit.fieldType = FieldType::NUMERIC;
        crit.numericValue = d;
        crit.stringValue.clear();
    } catch (...) {
        // treat as string
        crit.fieldType = FieldType::STRING;
        crit.stringValue = valPart;
        crit.numericValue = 0.0;
    }
    return true;
}

bool parseCriteria(const std::string &criteriaStr, std::vector<FilterCriterion> &criteria) {
    criteria.clear();
    std::vector<std::string> tokens;
    split(criteriaStr, ';', tokens);
    for (auto &t : tokens) {
        auto trimmed = trim(t);
        if (trimmed.empty())
            continue;
        FilterCriterion c;
        if (!parseSingleCriterion(trimmed, c)) {
            return false;
        }
        criteria.push_back(c);
    }
    if (criteria.empty()) {
        std::cerr << "Error: no valid criteria in '" << criteriaStr << "'.\n";
        return false;
    }
    return true;
}

// split by tab - uses vcfx::split_tabs for performance
static std::vector<std::string> tabSplit(const std::string &line) {
    return vcfx::split_tabs(line);
}

// parse info => either return key= or a flag
// if numeric, parse double
// if string => keep string
static bool getInfoValue(const std::string &infoField, const std::string &key, bool wantNumeric, double &numericOut,
                         std::string &stringOut) {
    if (infoField.empty() || infoField == ".")
        return false;
    std::vector<std::string> tokens;
    split(infoField, ';', tokens);
    for (auto &kv : tokens) {
        auto eq = kv.find('=');
        if (eq != std::string::npos) {
            std::string k = trim(kv.substr(0, eq));
            std::string v = trim(kv.substr(eq + 1));
            if (k == key) {
                if (wantNumeric) {
                    try {
                        double d = std::stod(v);
                        numericOut = d;
                        return true;
                    } catch (...) {
                        return false;
                    }
                } else {
                    stringOut = v;
                    return true;
                }
            }
        } else {
            // it's a flag
            std::string fl = trim(kv);
            if (fl == key) {
                if (wantNumeric) {
                    numericOut = 1.0;
                } else {
                    stringOut = fl;
                }
                return true;
            }
        }
    }
    return false;
}

static bool compareDouble(double x, FilterOp op, double y) {
    switch (op) {
    case FilterOp::GT:
        return (x > y);
    case FilterOp::GE:
        return (x >= y);
    case FilterOp::LT:
        return (x < y);
    case FilterOp::LE:
        return (x <= y);
    case FilterOp::EQ:
        return (x == y);
    case FilterOp::NE:
        return (x != y);
    }
    return false;
}

static bool compareString(const std::string &s, FilterOp op, const std::string &t) {
    switch (op) {
    case FilterOp::EQ:
        return (s == t);
    case FilterOp::NE:
        return (s != t);
    // if user tries < or > => invalid for strings in this code => fail
    default:
        // we do no partial compare. We fail this criterion
        return false;
    }
}

// We handle standard fields: POS(numeric), QUAL(numeric), FILTER(string?), plus info
// If fieldName== "POS", parse fields[1] => numeric
// If fieldName== "QUAL", parse fields[5] => numeric
// If fieldName== "FILTER", parse fields[6] => string
// else => treat as an INFO key
static bool evaluateCriterion(const std::vector<std::string> &fields, const FilterCriterion &c) {
    if (fields.size() < 8)
        return false;
    if (c.fieldName == "POS") {
        if (fields[1].empty())
            return false;
        try {
            double p = std::stod(fields[1]);
            return compareDouble(p, c.op, c.numericValue);
        } catch (...) {
            return false;
        }
    } else if (c.fieldName == "QUAL") {
        // index 5
        if (fields[5].empty() || fields[5] == ".") {
            // treat missing as 0? or fail?
            double q = 0.0;
            return compareDouble(q, c.op, c.numericValue);
        }
        try {
            double q = std::stod(fields[5]);
            return compareDouble(q, c.op, c.numericValue);
        } catch (...) {
            return false;
        }
    } else if (c.fieldName == "FILTER") {
        // index=6 => string compare
        if (c.fieldType == FieldType::NUMERIC) {
            // invalid usage => fail
            return false;
        }
        return compareString(fields[6], c.op, c.stringValue);
    } else {
        // treat as info key
        // c.fieldType => numeric/string
        double num = 0.0;
        std::string st;
        bool got = getInfoValue(fields[7], c.fieldName, (c.fieldType == FieldType::NUMERIC), num, st);
        if (!got)
            return false; // can't find key => fail
        if (c.fieldType == FieldType::NUMERIC) {
            return compareDouble(num, c.op, c.numericValue);
        } else {
            return compareString(st, c.op, c.stringValue);
        }
    }
}

bool recordPasses(const std::string &record, const std::vector<FilterCriterion> &criteria, bool useAndLogic) {
    // parse fields
    auto fields = tabSplit(record);
    if (fields.size() < 8)
        return false;
    bool anyPass = false;
    // check each criterion
    for (auto &crit : criteria) {
        bool pass = evaluateCriterion(fields, crit);
        if (useAndLogic) {
            if (!pass)
                return false;
        } else {
            // or logic
            if (pass)
                anyPass = true;
        }
    }
    if (useAndLogic) {
        return true; // all pass
    } else {
        return anyPass;
    }
}

void processVCF(std::istream &in, std::ostream &out, const std::vector<FilterCriterion> &criteria, bool useAndLogic) {
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
            // data line before #CHROM => skip with a warning?
            std::cerr << "Warning: data line before #CHROM => skipping.\n";
            continue;
        }
        // apply
        if (recordPasses(line, criteria, useAndLogic)) {
            out << line << "\n";
        }
    }
}

void printHelp() {
    std::cout
        << "VCFX_record_filter: Filter VCF data lines by multiple criteria.\n\n"
           "Usage:\n"
           "  VCFX_record_filter [options] --filter \"CRITERIA\"\n"
           "  < input.vcf > output.vcf\n\n"
           "Options:\n"
           "  --filter, -f \"...\"   One or more criteria separated by semicolons, e.g.\n"
           "                        \"POS>10000; QUAL>=30; AF<0.05; FILTER==PASS\"\n"
           "                        Each criterion must use an operator among >,>=,<,<=,==,!=\n\n"
           "  --logic and|or        'and' => a line must pass all criteria (default)\n"
           "                        'or'  => pass if any criterion is satisfied.\n"
           "  --help, -h            Show this help.\n\n"
           "Fields:\n"
           "  POS => numeric, QUAL => numeric, FILTER => string.\n"
           "  Others => assumed to be an INFO key. We try numeric parse if the criterion is numeric, else string.\n\n"
           "Example:\n"
           "  VCFX_record_filter --filter \"POS>=1000;FILTER==PASS;DP>10\" --logic and < in.vcf > out.vcf\n";
}

// main with typical argument parse
static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_record_filter", show_help))
        return 0;
    if (argc == 1) {
        printHelp();
        return 0;
    }
    bool showHelp = false;
    std::string criteriaStr;
    std::string logicStr = "and";
    static struct option long_opts[] = {{"help", no_argument, 0, 'h'},
                                        {"filter", required_argument, 0, 'f'},
                                        {"logic", required_argument, 0, 'l'},
                                        {0, 0, 0, 0}};
    while (true) {
        int c = ::getopt_long(argc, argv, "hf:l:", long_opts, nullptr);
        if (c == -1)
            break;
        switch (c) {
        case 'h':
            showHelp = true;
            break;
        case 'f':
            criteriaStr = optarg;
            break;
        case 'l':
            logicStr = optarg;
            break;
        default:
            showHelp = true;
        }
    }
    if (showHelp) {
        printHelp();
        return 0;
    }
    if (criteriaStr.empty()) {
        std::cerr << "Error: must provide --filter \"CRITERIA\".\n";
        printHelp();
        return 1;
    }
    bool useAndLogic = true;
    if (logicStr == "and") {
        useAndLogic = true;
    } else if (logicStr == "or") {
        useAndLogic = false;
    } else {
        std::cerr << "Error: logic must be 'and' or 'or'.\n";
        return 1;
    }
    std::vector<FilterCriterion> criteria;
    if (!parseCriteria(criteriaStr, criteria)) {
        std::cerr << "Error: failed to parse criteria.\n";
        return 1;
    }
    processVCF(std::cin, std::cout, criteria, useAndLogic);
    return 0;
}
