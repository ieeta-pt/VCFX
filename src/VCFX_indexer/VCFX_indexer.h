#ifndef VCFX_INDEXER_H
#define VCFX_INDEXER_H

#include <iostream>
#include <string>

class VCFXIndexer {
  public:
    int run(int argc, char *argv[]);

  private:
    void displayHelp();

    // Memory-mapped file indexing (fast path for file arguments)
    int createVCFIndexMmap(const char *filename, std::ostream &out);

    // Stdin-based indexing (fallback for pipes)
    void createVCFIndex(std::istream &in, std::ostream &out);
};

#endif // VCFX_INDEXER_H
