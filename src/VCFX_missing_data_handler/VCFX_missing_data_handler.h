#ifndef VCFX_MISSING_DATA_HANDLER_H
#define VCFX_MISSING_DATA_HANDLER_H

#include <iostream>
#include <string>
#include <vector>

/**
 * @brief Structure to hold command-line arguments for the missing data handler tool.
 */
struct Arguments {
    bool fill_missing = false;            ///< Flag indicating whether to impute missing genotypes.
    std::string default_genotype = "./."; ///< Default genotype to use for imputation.
    std::vector<std::string> input_files; ///< List of input VCF files. If empty, read from stdin.
    int numThreads = 0;                   ///< Number of threads (0 = auto-detect).
};

/**
 * @brief Displays the help message for the missing data handler tool.
 */
void printHelp();

/**
 * @brief Parses command-line arguments.
 *
 * @param argc Argument count.
 * @param argv Argument vector.
 * @param args Reference to Arguments structure to populate.
 * @return true if parsing is successful, false otherwise.
 */
bool parseArguments(int argc, char *argv[], Arguments &args);

/**
 * @brief Splits a string by a given delimiter.
 *
 * @param str The input string to split.
 * @param delimiter The character delimiter.
 * @return A vector of split substrings.
 */
std::vector<std::string> splitString(const std::string &str, char delimiter);

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
bool handleMissingDataAll(const Arguments &args);

#endif // VCFX_MISSING_DATA_HANDLER_H
