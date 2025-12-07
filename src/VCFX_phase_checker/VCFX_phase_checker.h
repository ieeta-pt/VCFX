#ifndef VCFX_PHASE_CHECKER_H
#define VCFX_PHASE_CHECKER_H

#include <iostream>
#include <string>
#include <string_view>
#include <vector>

/*
 * VCFXPhaseChecker:
 *   Reads a VCF, outputs only lines in which all samples' GT are fully phased.
 *   If any sample is unphased (or missing GT), logs a warning and discards line.
 *
 * Performance features:
 *   - Memory-mapped file I/O with SIMD line scanning (20-30x faster for files)
 *   - Zero-copy string parsing with string_view
 *   - 1MB output buffering
 *   - FORMAT field caching (GT index computed once per unique FORMAT)
 *   - Early termination on first unphased sample
 */
class VCFXPhaseChecker {
  public:
    int run(int argc, char *argv[]);

  private:
    void displayHelp();

    // Streaming stdin processing (fallback path, still optimized)
    void processVCF(std::istream &in, std::ostream &out);

    // Memory-mapped file processing (fast path)
    void filterPhaseCheckedMmap(const char *filepath, std::ostream &out);

    // Checks if a genotype string is "fully phased" (no '/' seen, multiple alleles separated by '|').
    // If GT is missing or partial => returns false (unphased).
    // Legacy method for compatibility - delegates to fast inline version
    bool isFullyPhased(const std::string &gt) const;

    // Quiet mode flag - suppresses warning messages to stderr
    bool quiet_ = false;
};

#endif // VCFX_PHASE_CHECKER_H
