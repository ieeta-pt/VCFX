// VCFX_distance_calculator.h
#ifndef VCFX_DISTANCE_CALCULATOR_H
#define VCFX_DISTANCE_CALCULATOR_H

#include <iostream>
#include <string>
#include <unordered_map>

// Structure to hold variant information (only CHROM and POS are used)
struct VCFVariant {
    std::string chrom;
    int pos;
};

// Structure to hold per-chromosome summary statistics
struct ChromStats {
    int count;          // Number of inter-variant distances computed
    long totalDistance; // Sum of all distances
    int minDistance;    // Minimum distance seen
    int maxDistance;    // Maximum distance seen
    ChromStats();
};

// Optimized distance calculator class
class VCFXDistanceCalculator {
public:
    int run(int argc, char *argv[]);

    // Original stdin-based processing (public for legacy compatibility)
    bool calculateDistances(std::istream &in, std::ostream &out);

private:
    bool quietMode = false;

    void displayHelp();

    // Optimized mmap-based processing
    bool processFileMmap(const char* filename, std::ostream &out);

    // Output summary statistics
    void outputSummary(const std::unordered_map<std::string, ChromStats>& chromStats);
};

// Legacy function declarations for backward compatibility
void printHelp();
bool parseVCFLine(const std::string &line, VCFVariant &variant);
bool calculateDistances(std::istream &in, std::ostream &out);

#endif // VCFX_DISTANCE_CALCULATOR_H
