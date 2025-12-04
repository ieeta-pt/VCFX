#include "VCFX_missing_data_handler.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>

/**
 * @brief Displays the help message for the missing data handler tool.
 */
void printHelp() {
    std::cout
        << "VCFX_missing_data_handler\n"
        << "Usage: VCFX_missing_data_handler [OPTIONS] [files...]\n\n"
        << "Options:\n"
        << "  --fill-missing, -f            Impute missing genotypes with a default value (e.g., ./.).\n"
        << "  --default-genotype, -d GEN    Specify the default genotype for imputation (default: ./.).\n"
        << "  --help, -h                    Display this help message and exit.\n\n"
        << "Description:\n"
        << "  Flags or imputes missing genotype data in one or more VCF files. By default, missing genotypes are\n"
        << "  left as is (flagged), but can be replaced with a specified genotype using --fill-missing.\n\n"
        << "Examples:\n"
        << "  1) Flag missing data from a single file:\n"
        << "       ./VCFX_missing_data_handler < input.vcf > flagged_output.vcf\n\n"
        << "  2) Impute missing data with './.' from multiple files:\n"
        << "       ./VCFX_missing_data_handler -f --default-genotype \"./.\" file1.vcf file2.vcf > "
           "combined_output.vcf\n";
}

/**
 * @brief Splits a string by a given delimiter.
 *
 * @param str The input string to split.
 * @param delimiter The character delimiter.
 * @return A vector of split substrings.
 */
std::vector<std::string> splitString(const std::string &str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string temp;
    while (std::getline(ss, temp, delimiter)) {
        tokens.push_back(temp);
    }
    return tokens;
}

/**
 * @brief Parses command-line arguments.
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 * @param args Reference to Arguments structure to populate.
 * @return true if parsing is successful, false otherwise.
 */
bool parseArguments(int argc, char *argv[], Arguments &args) {
    static struct option long_opts[] = {{"fill-missing", no_argument, 0, 'f'},
                                        {"default-genotype", required_argument, 0, 'd'},
                                        {"help", no_argument, 0, 'h'},
                                        {0, 0, 0, 0}};

    while (true) {
        int c = getopt_long(argc, argv, "fd:h", long_opts, nullptr);
        if (c == -1)
            break;
        switch (c) {
        case 'f':
            args.fill_missing = true;
            break;
        case 'd':
            args.default_genotype = optarg;
            break;
        case 'h':
        default:
            printHelp();
            exit(0);
        }
    }

    // The remainder arguments are input files
    // if none => read from stdin
    while (optind < argc) {
        args.input_files.push_back(argv[optind++]);
    }

    return true;
}

/**
 * @brief Process a single VCF stream. Writes to out.
 *
 * If fill_missing is true, we replace missing genotypes with args.default_genotype in the GT field.
 */
static bool processVCF(std::istream &in, std::ostream &out, bool fillMissing, const std::string &defaultGT) {
    std::string line;
    bool header_found = false;

    while (std::getline(in, line)) {
        if (line.empty()) {
            out << "\n";
            continue;
        }
        if (line[0] == '#') {
            out << line << "\n";
            if (line.rfind("#CHROM", 0) == 0) {
                header_found = true;
            }
            continue;
        }
        if (!header_found) {
            std::cerr << "Error: VCF data line encountered before #CHROM header.\n";
            // could continue or skip
            // we'll just continue
        }

        // parse columns
        auto fields = splitString(line, '\t');
        if (fields.size() < 9) {
            // invalid => pass as is
            out << line << "\n";
            continue;
        }
        // the 9th col => format
        std::string &format = fields[8];
        auto fmtParts = splitString(format, ':');
        int gtIndex = -1;
        for (size_t i = 0; i < fmtParts.size(); i++) {
            if (fmtParts[i] == "GT") {
                gtIndex = i;
                break;
            }
        }
        if (gtIndex < 0) {
            // no GT => just pass line
            out << line << "\n";
            continue;
        }
        // for each sample => fields[9..]
        for (size_t s = 9; s < fields.size(); s++) {
            auto sampleParts = splitString(fields[s], ':');
            if (gtIndex >= (int)sampleParts.size()) {
                // missing data for that sample => skip
                continue;
            }
            std::string &genotype = sampleParts[gtIndex];
            bool isMissing = false;
            if (genotype.empty() || genotype == "." || genotype == "./." || genotype == ".|.") {
                isMissing = true;
            }
            if (isMissing && fillMissing) {
                sampleParts[gtIndex] = defaultGT;
            }
            // rejoin sampleParts
            std::stringstream sampSS;
            for (size_t sp = 0; sp < sampleParts.size(); sp++) {
                if (sp > 0)
                    sampSS << ":";
                sampSS << sampleParts[sp];
            }
            fields[s] = sampSS.str();
        }
        // rejoin entire line
        std::stringstream lineSS;
        for (size_t c = 0; c < fields.size(); c++) {
            if (c > 0)
                lineSS << "\t";
            lineSS << fields[c];
        }
        out << lineSS.str() << "\n";
    }

    return true;
}

/**
 * @brief Processes the VCF file(s) to handle missing genotype data,
 *        either replacing missing data with a default genotype or leaving them flagged.
 *
 * @param args Command-line arguments specifying behavior.
 * @return true if processing is successful, false otherwise.
 *
 * This function reads from each file in args.input_files, or from stdin if none specified,
 * and writes the processed lines to stdout.
 */
bool handleMissingDataAll(const Arguments &args) {
    if (args.input_files.empty()) {
        // read from stdin
        return processVCF(std::cin, std::cout, args.fill_missing, args.default_genotype);
    } else {
        // handle multiple files in sequence => we simply process each one, writing to stdout
        bool firstFile = true;
        for (size_t i = 0; i < args.input_files.size(); i++) {
            std::string &path = (std::string &)args.input_files[i];
            std::ifstream fin(path);
            if (!fin.is_open()) {
                std::cerr << "Error: cannot open file " << path << "\n";
                return false;
            }
            if (!firstFile) {
                // skip printing the #CHROM header again?
                // a naive approach: we do a small trick:
                // read lines until #CHROM => skip them, then pass the rest
                // or we can pass everything => leads to repeated headers
                bool foundChrom = false;
                std::string line;
                while (std::getline(fin, line)) {
                    if (line.rfind("#CHROM", 0) == 0) {
                        // we've found #CHROM => use processVCF for the remainder
                        // but we already read one line, so let's put it back in the stream => complicated
                        // simpler approach: we do a manual approach
                        // We'll treat that #CHROM line as data for processVCF => but that might produce error
                        // We'll do: skip header lines
                        break;
                    } else if (line.empty()) {
                        // skip
                    } else if (line[0] == '#') {
                        // skip
                    } else {
                        // we found data => put it back => complicated
                        // We'll store it in a buffer
                        std::stringstream buffer;
                        buffer << line << "\n";
                        // now process the rest
                        // reinsert?
                        // We'll do an approach: we store lines in an in-memory stream
                        std::string nextLine;
                        while (std::getline(fin, nextLine)) {
                            buffer << nextLine << "\n";
                        }
                        // now buffer holds the entire
                        // pass buffer to processVCF
                        std::istringstream iss(buffer.str());
                        processVCF(iss, std::cout, args.fill_missing, args.default_genotype);
                        fin.close();
                        goto nextFile;
                    }
                }
                // now we pass the rest
                processVCF(fin, std::cout, args.fill_missing, args.default_genotype);
            nextFile:
                fin.close();
            } else {
                // first file => pass everything
                processVCF(fin, std::cout, args.fill_missing, args.default_genotype);
                fin.close();
                firstFile = false;
            }
        }
    }
    return true;
}

/**
 * @brief Main function for the missing data handler tool.
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 * @return int Exit status.
 */
static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_missing_data_handler", show_help))
        return 0;
    Arguments args;
    parseArguments(argc, argv, args);

    if (args.fill_missing) {
        std::cerr << "Info: Missing genotypes will be imputed with genotype: " << args.default_genotype << "\n";
    } else {
        std::cerr << "Info: Missing genotypes will be left flagged.\n";
    }

    bool success = handleMissingDataAll(args);
    return success ? 0 : 1;
}
