#include "VCFX_multiallelic_splitter.h"
#include "vcfx_core.h" 
#include "vcfx_io.h"
#include <cctype>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// OPTIMIZED: Use direct string parsing instead of stringstream
static std::vector<std::string> split(const std::string &s, char d) {
    std::vector<std::string> v;
    v.reserve(8);
    size_t start = 0, end;
    while ((end = s.find(d, start)) != std::string::npos) {
        v.emplace_back(s, start, end - start);
        start = end + 1;
    }
    v.emplace_back(s, start);
    return v;
}

static bool parseNumberEq(const std::string &line, std::string &id, std::string &num) {
    auto i = line.find("ID=");
    if (i == std::string::npos)
        return false;
    auto sub = line.substr(i + 3);
    auto e = sub.find_first_of(",>");
    if (e == std::string::npos)
        return false;
    id = sub.substr(0, e);
    auto n = line.find("Number=");
    if (n == std::string::npos)
        return false;
    auto sub2 = line.substr(n + 7);
    auto e2 = sub2.find_first_of(",>");
    if (e2 == std::string::npos)
        return false;
    num = sub2.substr(0, e2);
    return true;
}

static void parseHeaderLine(const std::string &line, VCFHeaderInfo &hdr) {
    if (line.rfind("##INFO=", 0) == 0) {
        std::string id, num;
        if (!parseNumberEq(line, id, num))
            return;
        SubfieldMeta m;
        m.isInfo = true;
        m.id = id;
        m.number = num;
        hdr.meta[id] = m;
    } else if (line.rfind("##FORMAT=", 0) == 0) {
        std::string id, num;
        if (!parseNumberEq(line, id, num))
            return;
        SubfieldMeta m;
        m.isFormat = true;
        m.id = id;
        m.number = num;
        hdr.meta[id] = m;
    }
}

void printHelp() {
    std::cout << "VCFX_multiallelic_splitter:\n"
                 "  Splits multi-allelic variants into multiple lines, rewriting GT/AD/PL.\n"
                 "Usage:\n"
                 "  VCFX_multiallelic_splitter [options] < input.vcf > output.vcf\n"
                 "  --help, -h\n";
}

static bool parseVariantLine(const std::string &line, VCFVariant &var, std::vector<std::string> &fields) {
    if (line.empty() || line[0] == '#')
        return false;
    vcfx::split_tabs(line, fields);
    if (fields.size() < 9)
        return false;
    auto &f = fields;
    var.chrom = f[0];
    try {
        var.pos = std::stoi(f[1]);
    } catch (...) {
        return false;
    }
    var.id = f[2];
    var.ref = f[3];
    {
        var.alt.clear();
        auto al = split(f[4], ',');
        for (auto &x : al)
            var.alt.push_back(x);
    }
    var.qual = f[5];
    var.filter = f[6];
    var.info = f[7];
    var.formatStr = f[8];
    var.samples.clear();
    for (size_t i = 9; i < f.size(); i++) {
        var.samples.push_back(f[i]);
    }
    return true;
}

static std::vector<std::string> splitInfo(const std::string &inf) {
    if (inf == "." || inf.empty())
        return {};
    auto seg = split(inf, ';');
    std::vector<std::string> ret;
    for (auto &x : seg) {
        if (!x.empty())
            ret.push_back(x);
    }
    return ret;
}

static std::string join(const std::vector<std::string> &v, const std::string &d) {
    if (v.empty())
        return ".";
    std::ostringstream o;
    for (size_t i = 0; i < v.size(); i++) {
        if (i > 0)
            o << d;
        o << v[i];
    }
    return o.str();
}

static bool isInteger(const std::string &s) {
    for (char c : s)
        if (!std::isdigit((unsigned char)c) && c != '-')
            return false;
    return !s.empty();
}

static std::vector<std::string> recA(const std::vector<std::string> &vals, int altIx) {
    if ((size_t)altIx >= vals.size())
        return {"."};
    return {vals[altIx]};
}

static std::vector<std::string> recR(const std::vector<std::string> &vals, int altIx) {
    if (vals.empty())
        return {"."};
    if ((size_t)altIx >= vals.size())
        return {vals[0], "."};
    return {vals[0], vals[altIx]};
}

static std::vector<std::string> recG(const std::vector<std::string> &vals, int altIx, int nAlts) {
    if (vals.size() != (size_t)((nAlts + 1) * (nAlts + 2) / 2))
        return {"."};
    auto idxOf = [&](int a, int b) { return ((2 * nAlts + 1 - a) * a) / 2 + (b - a); };
    int i00 = idxOf(0, 0), i01 = idxOf(0, altIx), i11 = idxOf(altIx, altIx);
    if (i00 < 0 || i01 < 0 || i11 < 0)
        return {"."};
    if (i00 >= (int)vals.size() || i01 >= (int)vals.size() || i11 >= (int)vals.size())
        return {"."};
    return {vals[i00], vals[i01], vals[i11]};
}

static std::string recodeINFO(const std::string &info, int altIx, int nAlts, const VCFHeaderInfo &hdr) {
    auto items = splitInfo(info);
    if (items.empty())
        return info;
    std::vector<std::string> out;
    for (auto &item : items) {
        auto eq = item.find('=');
        if (eq == std::string::npos) {
            out.push_back(item);
            continue;
        }
        auto key = item.substr(0, eq);
        auto val = item.substr(eq + 1);
        auto it = hdr.meta.find(key);
        if (it == hdr.meta.end() || !it->second.isInfo) {
            out.push_back(item);
            continue;
        }
        auto arr = split(val, ',');
        auto &num = it->second.number;
        if (num == "A") {
            auto x = recA(arr, altIx);
            if (x.size() == 1 && x[0] == ".")
                out.push_back(key + "=.");
            else
                out.push_back(key + "=" + join(x, ","));
        } else if (num == "R") {
            auto x = recR(arr, altIx + 1);
            if (x.size() == 2 && x[0] == "." && x[1] == ".")
                out.push_back(key + "=.");
            else
                out.push_back(key + "=" + join(x, ","));
        } else if (num == "G") {
            auto x = recG(arr, altIx + 1, nAlts);
            if (x.size() == 1 && x[0] == ".")
                out.push_back(key + "=.");
            else
                out.push_back(key + "=" + join(x, ","));
        } else if (num == "1") {
            out.push_back(item);
        } else if (num == "." || isInteger(num)) {
            out.push_back(item);
        } else {
            out.push_back(item);
        }
    }
    if (out.empty())
        return ".";
    std::ostringstream oss;
    for (size_t i = 0; i < out.size(); i++) {
        if (i > 0)
            oss << ";";
        oss << out[i];
    }
    return oss.str();
}

static std::vector<std::string> recField(const std::vector<std::string> &vals, int altIx, int nAlts,
                                         const std::string &key, const VCFHeaderInfo &hdr) {
    auto it = hdr.meta.find(key);
    if (it == hdr.meta.end() || !it->second.isFormat)
        return vals;
    auto &num = it->second.number;
    if (num == "A") {
        auto x = recA(vals, altIx - 1);
        if (x.size() == 1 && x[0] == ".")
            return {"."};
        return x;
    } else if (num == "R") {
        auto x = recR(vals, altIx);
        if (x.size() == 2 && x[0] == "." && x[1] == ".")
            return {"."};
        return x;
    } else if (num == "G") {
        auto x = recG(vals, altIx, nAlts);
        if (x.size() == 1 && x[0] == ".")
            return {"."};
        return x;
    }
    return vals;
}

static std::string recSample(const std::string &sample, const std::vector<std::string> &fmtKeys, int altIx, int nAlts,
                             const VCFHeaderInfo &hdr) {
    auto subs = split(sample, ':');
    if (subs.size() < fmtKeys.size())
        subs.resize(fmtKeys.size(), ".");
    for (size_t i = 0; i < fmtKeys.size(); i++) {
        if (i >= subs.size())
            break;
        auto &k = fmtKeys[i];
        if (k == "GT") {
            std::string g = subs[i];
            for (auto &c : g)
                if (c == '|')
                    c = '/';
            auto d = g.find('/');
            if (d == std::string::npos) {
                subs[i] = ".";
                continue;
            }
            auto a1 = g.substr(0, d), a2 = g.substr(d + 1);
            auto conv = [&](const std::string &x) {
                if (x == "0")
                    return "0";
                if (x == std::to_string(altIx))
                    return "1";
                return ".";
            };
            std::string na1 = conv(a1), na2 = conv(a2);
            if (na1 == "." && na2 == ".")
                subs[i] = ".";
            else
                subs[i] = na1 + "/" + na2;
        } else {
            auto arr = split(subs[i], ',');
            auto newA = recField(arr, altIx, nAlts, k, hdr);
            if (newA.size() == 1 && newA[0] == ".")
                subs[i] = ".";
            else {
                std::ostringstream o;
                for (size_t z = 0; z < newA.size(); z++) {
                    if (z > 0)
                        o << ",";
                    o << newA[z];
                }
                subs[i] = o.str();
            }
        }
    }
    std::ostringstream oss;
    for (size_t i2 = 0; i2 < subs.size(); i2++) {
        if (i2 > 0)
            oss << ":";
        oss << subs[i2];
    }
    return oss.str();
}

bool splitMultiAllelicVariants(std::istream &in, std::ostream &out) {
    VCFHeaderInfo hdr;
    bool foundChrom = false;

    // Performance: reuse containers across iterations
    std::vector<std::string> fields;
    fields.reserve(16);

    while (true) {
        auto p = in.tellg();
        if (!in.good())
            break;
        std::string line;
        if (!std::getline(in, line))
            break;
        if (line.empty()) {
            out << line << "\n";
            continue;
        }
        if (line[0] == '#') {
            out << line << "\n";
            if (line.rfind("##", 0) == 0)
                parseHeaderLine(line, hdr);
            if (line.rfind("#CHROM", 0) == 0) {
                foundChrom = true;
                break;
            }
        } else {
            in.seekg(p);
            break;
        }
    }
    if (!foundChrom) {
        while (true) {
            std::string l;
            if (!std::getline(in, l))
                break;
            out << l << "\n";
        }
        return true;
    }
    while (true) {
        std::string line;
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
        VCFVariant var;
        if (!parseVariantLine(line, var, fields)) {
            out << line << "\n";
            continue;
        }
        if (var.alt.size() <= 1) {
            out << line << "\n";
            continue;
        }
        int nAlts = var.alt.size();
        auto fmtKeys = split(var.formatStr, ':');
        auto rewriteInfo = [&](int altIx) { return recodeINFO(var.info, altIx, nAlts, hdr); };
        for (int a = 0; a < nAlts; a++) {
            std::ostringstream row;
            row << var.chrom << "\t" << var.pos << "\t" << var.id << "\t" << var.ref << "\t" << var.alt[a] << "\t"
                << var.qual << "\t" << var.filter << "\t" << rewriteInfo(a) << "\t" << var.formatStr;
            for (auto &sm : var.samples) {
                row << "\t" << recSample(sm, fmtKeys, a + 1, nAlts, hdr);
            }
            out << row.str() << "\n";
        }
    }
    return true;
}

static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_multiallelic_splitter", show_help))
        return 0;
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            printHelp();
            return 0;
        }
    }
    splitMultiAllelicVariants(std::cin, std::cout);
    return 0;
}
