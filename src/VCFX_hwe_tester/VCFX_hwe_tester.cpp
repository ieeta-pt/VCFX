#include "VCFX_hwe_tester.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

// Constructor-like run method
int VCFXHWETester::run(int argc, char *argv[]) {
    int opt;
    bool showHelp = false;
    static struct option long_options[] = {{"help", no_argument, 0, 'h'}, {0, 0, 0, 0}};
    while ((opt = getopt_long(argc, argv, "h", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'h':
        default:
            showHelp = true;
        }
    }
    if (showHelp) {
        displayHelp();
        return 0;
    }

    // Perform HWE on stdin
    performHWE(std::cin);
    return 0;
}

void VCFXHWETester::displayHelp() {
    std::cout << "VCFX_hwe_tester: Perform Hardy-Weinberg Equilibrium (HWE) tests on a biallelic VCF.\n\n"
              << "Usage:\n"
              << "  VCFX_hwe_tester [options] < input.vcf > output\n\n"
              << "Description:\n"
              << "  Reads each variant line, ignoring multi-allelic calls. For biallelic lines,\n"
              << "  collects genotypes as 0/0, 0/1, 1/1, then uses an exact test to produce\n"
              << "  a p-value for HWE.\n\n"
              << "Example:\n"
              << "  VCFX_hwe_tester < input.vcf > results.txt\n";
}

// Single definition of isBiallelic
bool VCFXHWETester::isBiallelic(const std::string &alt) {
    // If ALT has a comma => multiple alt alleles => not biallelic
    return (alt.find(',') == std::string::npos);
}

// parseGenotypes
bool VCFXHWETester::parseGenotypes(const std::vector<std::string> &genotypes, int &homRef, int &het, int &homAlt) {
    homRef = 0;
    het = 0;
    homAlt = 0;
    for (auto &gt : genotypes) {
        if (gt.empty() || gt == "." || gt == "./." || gt == ".|.") {
            continue; // missing
        }
        // unify '|' -> '/'
        std::string g = gt;
        for (char &c : g) {
            if (c == '|')
                c = '/';
        }
        size_t delim = g.find('/');
        if (delim == std::string::npos) {
            // Not diploid => skip
            continue;
        }
        std::string a1 = g.substr(0, delim);
        std::string a2 = g.substr(delim + 1);
        if (a1 == "0" && a2 == "0") {
            homRef++;
        } else if ((a1 == "0" && a2 == "1") || (a1 == "1" && a2 == "0")) {
            het++;
        } else if (a1 == "1" && a2 == "1") {
            homAlt++;
        } else {
            // e.g. a2>1 => multi-allelic => skip or fail
            return false;
        }
    }
    return true;
}

double VCFXHWETester::genotypeProbability(int homRef, int het, int homAlt) {
    int N = homRef + het + homAlt;
    int x = 2 * homAlt + het;
    int y = 2 * homRef + het;

    // We'll do enumerations with log factorial
    static const int MAXN = 20000;
    static std::vector<double> logFac;
    static int cachedN = 0;
    if (2 * N > cachedN) {
        logFac.resize(2 * N + 1, 0.0);
        double cum = 0.0;
        for (int i = 1; i <= 2 * N; i++) {
            cum += std::log((double)i);
            logFac[i] = cum;
        }
        cachedN = 2 * N;
    }

    auto logFactorial = [&](int n) -> double {
        if (n <= 0)
            return 0.0;
        if (n > cachedN)
            return 0.0;
        return logFac[n];
    };

    auto logProb = [&](int a) -> double {
        if ((x - a) < 0 || (y - a) < 0)
            return -INFINITY;
        if (((x - a) % 2) != 0 || ((y - a) % 2) != 0)
            return -INFINITY;
        int hr = (y - a) / 2;
        int ha = (x - a) / 2;
        if (hr < 0 || ha < 0)
            return -INFINITY;
        // log multinomial
        double lcoef = logFactorial(N) - (logFactorial(hr) + logFactorial(a) + logFactorial(ha));
        double p = (double)y / (double)(x + y);
        double q = (double)x / (double)(x + y);
        if (p <= 0.0 || q <= 0.0)
            return -INFINITY;
        double logTerm = hr * 2 * std::log(p) + a * std::log(2.0 * p * q) + ha * 2 * std::log(q);
        return lcoef + logTerm;
    };

    int observedHet = het;
    double logObs = logProb(observedHet);
    if (logObs == -INFINITY)
        return 1.0;

    double minLog = logObs;
    int maxA = std::min(x, y);
    std::vector<double> logVals(maxA + 1, -INFINITY);
    logVals[observedHet] = logObs;

    if (logObs < minLog)
        minLog = logObs;

    for (int a = 0; a <= maxA; a++) {
        if (a == observedHet)
            continue;
        double lp = logProb(a);
        logVals[a] = lp;
        if (lp < minLog && lp != -INFINITY)
            minLog = lp;
    }

    double sum = 0.0, obsExp = 0.0;
    for (int a = 0; a <= maxA; a++) {
        if (logVals[a] == -INFINITY)
            continue;
        double rel = logVals[a] - minLog;
        double e = std::exp(rel);
        sum += e;
        if (a == observedHet)
            obsExp = e;
    }
    double probObs = obsExp / sum;
    double pVal = 0.0;
    for (int a = 0; a <= maxA; a++) {
        if (logVals[a] == -INFINITY)
            continue;
        double rel = logVals[a] - minLog;
        double e = std::exp(rel);
        double prob = e / sum;
        if (prob <= probObs + 1e-12) {
            pVal += prob;
        }
    }
    if (pVal > 1.0)
        pVal = 1.0;
    if (pVal < 0.0)
        pVal = 0.0;
    return pVal;
}

double VCFXHWETester::calculateHWE(int homRef, int het, int homAlt) {
    int N = homRef + het + homAlt;
    if (N < 1)
        return 1.0;
    int x = 2 * homAlt + het;
    int y = 2 * homRef + het;
    if (x + y != 2 * N) {
        return 1.0;
    }
    return genotypeProbability(homRef, het, homAlt);
}

void VCFXHWETester::performHWE(std::istream &in) {
    std::string line;
    std::vector<std::string> fields;
    fields.reserve(16);
    // output header
    std::cout << "CHROM\tPOS\tID\tREF\tALT\tHWE_pvalue\n";
    while (std::getline(in, line)) {
        if (line.empty())
            continue;
        if (line[0] == '#')
            continue;
        vcfx::split_tabs(line, fields);
        if (fields.size() < 10)
            continue;
        std::string chrom = fields[0];
        std::string pos = fields[1];
        std::string id = fields[2];
        std::string ref = fields[3];
        std::string alt = fields[4];
        // skip if multi-allelic
        if (!isBiallelic(alt))
            continue;

        // find GT in format
        std::string format = fields[8];
        std::vector<std::string> fmts;
        {
            std::stringstream fs(format);
            std::string x;
            while (std::getline(fs, x, ':')) {
                fmts.push_back(x);
            }
        }
        int gt_index = -1;
        for (size_t i = 0; i < fmts.size(); i++) {
            if (fmts[i] == "GT") {
                gt_index = i;
                break;
            }
        }
        if (gt_index < 0)
            continue;

        // gather genotypes
        std::vector<std::string> gts;
        for (size_t s = 9; s < fields.size(); s++) {
            std::stringstream sampSS(fields[s]);
            std::vector<std::string> sampToks;
            {
                std::string xx;
                while (std::getline(sampSS, xx, ':')) {
                    sampToks.push_back(xx);
                }
            }
            if ((size_t)gt_index < sampToks.size()) {
                gts.push_back(sampToks[gt_index]);
            } else {
                gts.push_back(".");
            }
        }
        int hr = 0, h = 0, ha = 0;
        bool ok = parseGenotypes(gts, hr, h, ha);
        if (!ok)
            continue;
        double pVal = calculateHWE(hr, h, ha);
        std::cout << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\t" << std::fixed
                  << std::setprecision(6) << pVal << "\n";
    }
}

// actual main
static void show_help() {
    VCFXHWETester obj;
    char arg0[] = "VCFX_hwe_tester";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_hwe_tester", show_help))
        return 0;
    VCFXHWETester tester;
    return tester.run(argc, argv);
}
