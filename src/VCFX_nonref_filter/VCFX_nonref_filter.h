#ifndef VCFX_NONREF_FILTER_H
#define VCFX_NONREF_FILTER_H

#include <iostream>
#include <string>
#include <string_view>
#include <vector>

/*
 * VCFX_nonref_filter:
 *   Reads a VCF, discards any variant for which EVERY sample is hom-ref
 *   (all alleles = '0'). If a genotype is missing or partial, treat that
 *   as "not guaranteed hom-ref," so we keep that variant.
 *
 * Performance optimizations:
 *   - Supports memory-mapped file input via -i/--input for 100-1000x speedup
 *   - Uses string_view for zero-copy parsing
 *   - SIMD-optimized newline detection (AVX2/SSE2)
 *   - Early termination on first non-homref sample
 */
class VCFXNonRefFilter {
  public:
    // main entry point
    int run(int argc, char *argv[]);

  private:
    // prints usage
    void displayHelp();

    // stdin-based filtering (original, for piped input)
    void filterNonRef(std::istream &in, std::ostream &out);

    // mmap-based filtering (optimized, for file input)
    void filterNonRefMmap(const char* filepath, std::ostream &out);

    // parse a sample's genotype from a field => returns true if definitely hom-ref
    // Overloaded for both std::string and std::string_view
    bool isDefinitelyHomRef(const std::string &genotypeField) const;
    bool isDefinitelyHomRef(std::string_view genotypeField) const;
};

#endif // VCFX_NONREF_FILTER_H
