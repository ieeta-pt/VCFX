#ifndef VCFX_DUPLICATE_REMOVER_H
#define VCFX_DUPLICATE_REMOVER_H

#include <iostream>
#include <string>
#include <unordered_set>

// Class to encapsulate duplicate remover functionality
class VCFXDuplicateRemover {
  public:
    // Entry point for the tool
    int run(int argc, char *argv[]);

  private:
    // Displays the help message
    void displayHelp();

    // Remove duplicates from stdin
    bool removeDuplicates(std::istream &in, std::ostream &out);

    // Process file using memory-mapped I/O (optimized for large files)
    bool processFileMmap(const char* filename, std::ostream &out);

    // Quiet mode flag
    bool quietMode = false;
};

// Structure to represent a unique variant key
struct VariantKey {
    std::string chrom;
    int pos;
    std::string ref;
    std::string alt; // normalized: sorted, comma-separated alleles

    bool operator==(const VariantKey &other) const {
        return chrom == other.chrom && pos == other.pos && ref == other.ref && alt == other.alt;
    }
};

// Custom hash function for VariantKey using a hash-combine approach
struct VariantKeyHash {
    std::size_t operator()(const VariantKey &k) const {
        std::size_t h1 = std::hash<std::string>()(k.chrom);
        std::size_t h2 = std::hash<int>()(k.pos);
        std::size_t h3 = std::hash<std::string>()(k.ref);
        std::size_t h4 = std::hash<std::string>()(k.alt);
        std::size_t seed = h1;
        seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= h4 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

#endif // VCFX_DUPLICATE_REMOVER_H
