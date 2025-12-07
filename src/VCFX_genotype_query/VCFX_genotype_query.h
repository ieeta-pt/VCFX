#ifndef VCFX_GENOTYPE_QUERY_H
#define VCFX_GENOTYPE_QUERY_H

#include <iostream>
#include <string>

// Parses command-line arguments to get the desired genotype query and options
// Supports -g/--genotype-query, -i/--input file, --strict, and -q/--quiet flags
bool parseArguments(int argc, char *argv[], std::string &genotype_query,
                    bool &strictCompare, std::string &inputFile, bool &quiet);

// Prints usage/help
void printHelp();

// Filters a VCF from 'in', writing only records that contain at least one sample with the requested genotype
// Original signature for backward compatibility (quiet = false)
void genotypeQuery(std::istream &in, std::ostream &out, const std::string &genotype_query, bool strictCompare);

// Stream-based implementation with quiet option
void genotypeQueryStream(std::istream &in, std::ostream &out, const std::string &genotype_query,
                          bool strictCompare, bool quiet);

#endif // VCFX_GENOTYPE_QUERY_H
