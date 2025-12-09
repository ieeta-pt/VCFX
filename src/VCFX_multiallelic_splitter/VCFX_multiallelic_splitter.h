#ifndef VCFX_MULTIALLELIC_SPLITTER_H
#define VCFX_MULTIALLELIC_SPLITTER_H

#include <iostream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

/* Simple structure describing whether an ID is an INFO or FORMAT field, plus its Number field */
struct SubfieldMeta {
    bool isInfo = false, isFormat = false;
    std::string id, number;
};

/* We store an ID->SubfieldMeta map for known fields from the header. */
struct VCFHeaderInfo {
    std::unordered_map<std::string, SubfieldMeta> meta;
};

/* A single variant record with possibly multiple ALT alleles. */
struct VCFVariant {
    std::string chrom;
    int pos;
    std::string id;
    std::string ref;
    std::vector<std::string> alt;
    std::string qual;
    std::string filter;
    std::string info;
    std::string formatStr;
    std::vector<std::string> samples; // one element per sample
};

/* Optimized multiallelic splitter class */
class VCFXMultiallelicSplitter {
public:
    int run(int argc, char *argv[]);

private:
    bool quietMode = false;

    void displayHelp();

    // Original stdin-based processing
    bool splitMultiAllelicVariants(std::istream &in, std::ostream &out);

    // Optimized mmap-based processing
    bool processFileMmap(const char* filename, std::ostream &out);

    // Header parsing
    void parseHeaderLine(const char* line, size_t len, VCFHeaderInfo &hdr);

    // Optimized recoding functions
    void recodeInfoField(const char* info, size_t infoLen, int altIdx, int nAlts,
                         const VCFHeaderInfo& hdr, std::string& out);

    void recodeSample(const char* sample, size_t sampleLen,
                      const std::vector<std::string_view>& fmtKeys,
                      int altIdx, int nAlts, const VCFHeaderInfo& hdr,
                      std::string& out);
};

/* Print minimal help/usage. */
void printHelp();

/* Reads a VCF from 'in', writes lines to 'out' with multi-allelic sites split,
 * rewriting subfields (GT, AD, PL) for each splitted line.
 * Header lines are passed unmodified. */
bool splitMultiAllelicVariants(std::istream &in, std::ostream &out);

#endif
