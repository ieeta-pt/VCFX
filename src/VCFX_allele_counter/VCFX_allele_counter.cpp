#include "vcfx_core.h" 
#include "vcfx_io.h"
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// ---------------------------------------------------------------------
// Structures and Declarations
// ---------------------------------------------------------------------
struct AlleleCounterArguments {
    std::vector<std::string> samples;
};

static void printHelp();
static bool parseArguments(int argc, char *argv[], AlleleCounterArguments &args);
static std::vector<std::string> splitString(const std::string &str, char delimiter);
static bool countAlleles(std::istream &in, std::ostream &out, const AlleleCounterArguments &args);

// ---------------------------------------------------------------------
// printHelp
// ---------------------------------------------------------------------
static void printHelp() {
    std::cout << "VCFX_allele_counter\n"
              << "Usage: VCFX_allele_counter [OPTIONS] < input.vcf > allele_counts.tsv\n\n"
              << "Options:\n"
              << "  --samples, -s \"Sample1 Sample2\"   Specify the sample names to include.\n"
              << "                                     If not specified, all samples are processed.\n"
              << "  --help, -h                        Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Reads a VCF from stdin and outputs a TSV file with the columns:\n"
              << "    CHROM  POS  ID  REF  ALT  Sample  Ref_Count  Alt_Count\n\n"
              << "  Each sample for each variant is listed. Alleles are determined by\n"
              << "  genotype strings (GT). This code treats any numeric allele that is\n"
              << "  not '0' as ALT, e.g. '1' or '2' or '3' => alt.\n\n"
              << "Examples:\n"
              << "  1) Count alleles for SampleA and SampleB:\n"
              << "     VCFX_allele_counter --samples \"SampleA SampleB\" < input.vcf > allele_counts.tsv\n\n"
              << "  2) Count alleles for all samples:\n"
              << "     VCFX_allele_counter < input.vcf > allele_counts_all.tsv\n";
}

// ---------------------------------------------------------------------
// parseArguments
// ---------------------------------------------------------------------
static bool parseArguments(int argc, char *argv[], AlleleCounterArguments &args) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--samples" || arg == "-s") && i + 1 < argc) {
            std::string samplesStr = argv[++i];
            args.samples = splitString(samplesStr, ' ');
            // Trim whitespace from each sample name
            for (auto &sample : args.samples) {
                // left trim
                sample.erase(0, sample.find_first_not_of(" \t\n\r\f\v"));
                // right trim
                sample.erase(sample.find_last_not_of(" \t\n\r\f\v") + 1);
            }
        } else if (arg == "--help" || arg == "-h") {
            printHelp();
            std::exit(0);
        } else {
            std::cerr << "Warning: Unrecognized argument '" << arg << "'.\n";
        }
    }
    return true;
}

// ---------------------------------------------------------------------
// splitString
// ---------------------------------------------------------------------
static std::vector<std::string> splitString(const std::string &str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string temp;
    while (std::getline(ss, temp, delimiter)) {
        tokens.push_back(temp);
    }
    return tokens;
}

// ---------------------------------------------------------------------
// countAlleles
// ---------------------------------------------------------------------
static bool countAlleles(std::istream &in, std::ostream &out, const AlleleCounterArguments &args) {
    std::string line;
    std::vector<std::string> headerFields;
    bool foundChromHeader = false;

    // We will store the columns (indexes) for the samples we care about
    std::vector<int> sampleIndices;

    // Print one TSV header (not the original VCF header lines)
    out << "CHROM\tPOS\tID\tREF\tALT\tSample\tRef_Count\tAlt_Count\n";

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        // Skip lines that start with '#' (VCF headers), but parse #CHROM line
        if (line[0] == '#') {
            // If this is the #CHROM line, parse to find sample columns
            if (line.rfind("#CHROM", 0) == 0) {
                headerFields = splitString(line, '\t');
                if (headerFields.size() < 9) {
                    std::cerr << "Error: #CHROM line has fewer than 9 columns.\n";
                    return false;
                }
                // Build sample->index map
                std::unordered_map<std::string, int> sampleMap;
                for (size_t i = 9; i < headerFields.size(); ++i) {
                    sampleMap[headerFields[i]] = static_cast<int>(i);
                }

                // Decide which samples to process
                if (!args.samples.empty()) {
                    // Use only user-specified samples
                    for (const auto &s : args.samples) {
                        auto it = sampleMap.find(s);
                        if (it == sampleMap.end()) {
                            std::cerr << "Error: Sample '" << s << "' not found in VCF header.\n";
                            return false;
                        }
                        sampleIndices.push_back(it->second);
                    }
                } else {
                    // If no samples are specified, process all
                    for (size_t i = 9; i < headerFields.size(); ++i) {
                        sampleIndices.push_back((int)i);
                    }
                }
                foundChromHeader = true;
            }
            // We do NOT print the raw VCF headers to output
            continue;
        }

        // We expect a #CHROM line before data lines
        if (!foundChromHeader) {
            std::cerr << "Error: VCF #CHROM header not found before records.\n";
            return false;
        }

        // Split VCF data line
        auto fields = splitString(line, '\t');
        if (fields.size() < 9) {
            std::cerr << "Warning: Skipping invalid VCF line with fewer than 9 fields:\n" << line << "\n";
            continue;
        }

        // Basic columns
        const std::string &chrom = fields[0];
        const std::string &pos = fields[1];
        const std::string &id = fields[2];
        const std::string &ref = fields[3];
        const std::string &alt = fields[4];

        // For each chosen sample, parse genotype
        for (int sIndex : sampleIndices) {
            if (sIndex >= (int)fields.size()) {
                std::cerr << "Warning: Sample index " << sIndex << " out of range for line:\n" << line << "\n";
                continue;
            }
            const std::string &sampleField = fields[sIndex];

            // genotype is the portion before the first ':'
            size_t colonPos = sampleField.find(':');
            std::string gt = (colonPos == std::string::npos) ? sampleField : sampleField.substr(0, colonPos);

            // Replace '|' with '/' for a uniform split
            for (auto &c : gt) {
                if (c == '|')
                    c = '/';
            }
            auto alleles = splitString(gt, '/');
            if (alleles.empty()) {
                // e.g. '.' or blank
                continue;
            }

            int refCount = 0;
            int altCount = 0;
            for (auto &allele : alleles) {
                // Skip empty or '.' allele
                if (allele.empty() || allele == ".") {
                    continue;
                }
                // If numeric & "0" => ref, else alt
                bool numeric = true;
                for (char c : allele) {
                    if (!isdigit(c)) {
                        numeric = false;
                        break;
                    }
                }
                if (!numeric) {
                    // e.g. unknown character => skip
                    continue;
                }
                // If allele == "0" => ref, else alt
                if (allele == "0") {
                    refCount++;
                } else {
                    altCount++;
                }
            }

            // The sample's name
            std::string sampleName =
                (sIndex < (int)headerFields.size()) ? headerFields[sIndex] : ("Sample_" + std::to_string(sIndex));

            // Output as TSV row
            out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\t" << sampleName << "\t"
                << refCount << "\t" << altCount << "\n";
        }
    }
    return true;
}

// ---------------------------------------------------------------------
// main()
// ---------------------------------------------------------------------
static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_allele_counter", show_help))
        return 0;
    AlleleCounterArguments args;
    parseArguments(argc, argv, args);

    if (!args.samples.empty()) {
        std::cerr << "Info: Counting alleles for these samples:\n";
        for (auto &s : args.samples) {
            std::cerr << "  " << s << "\n";
        }
    } else {
        std::cerr << "Info: Counting alleles for ALL samples.\n";
    }

    bool success = countAlleles(std::cin, std::cout, args);
    return success ? 0 : 1;
}
