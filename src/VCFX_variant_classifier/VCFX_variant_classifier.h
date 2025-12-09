#ifndef VCFX_VARIANT_CLASSIFIER_H
#define VCFX_VARIANT_CLASSIFIER_H

#include <string>
#include <vector>
#include <string_view>

enum class VariantType { SNP, INDEL, MNV, STRUCTURAL, UNKNOWN };

class VCFXVariantClassifier {
  public:
    int run(int argc, char *argv[]);

  private:
    // Show usage
    void displayHelp();

    // Identify the user's requested output mode
    bool appendInfo = false; // if true, output a fully valid VCF with classification appended to INFO
    // otherwise produce a TSV with columns: CHROM POS ID REF ALT Classification
    bool quietMode = false;

    // The main method that reads lines from input, classifies, and writes output
    void classifyStream(std::istream &in, std::ostream &out);

    // Optimized mmap-based file processing
    bool processFileMmap(const char* filename, std::ostream &out);

    // Detect alt allele as structural if it is symbolic (<DEL>) or breakend notation ([chr or ]chr) or length
    // difference >=50
    bool isStructuralAllele(const std::string &alt) const;
    bool isStructuralAlleleSV(std::string_view alt) const;

    // Classify a single (ref, alt)
    VariantType classifyAllele(const std::string &ref, const std::string &alt) const;
    VariantType classifyAlleleSV(std::string_view ref, std::string_view alt) const;

    // From the set of alt alleles, find the final classification
    VariantType classifyVariant(const std::string &ref, const std::vector<std::string> &alts) const;

    // Stringify the variant type
    std::string typeToStr(VariantType t) const;
    const char* typeToStrFast(VariantType t) const;

    // If we are in append‐info mode, we parse line into columns, parse alt as multiple, classify, then append e.g.
    // VCF_CLASS=xxx
    std::string appendClassification(const std::string &line);

    // Splitting helper
    std::vector<std::string> split(const std::string &s, char delim) const;
};

#endif
