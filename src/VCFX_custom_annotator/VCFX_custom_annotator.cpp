#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// ---------------------------------------------------------------------------
// Class: VCFXCustomAnnotator
// ---------------------------------------------------------------------------
class VCFXCustomAnnotator {
  public:
    int run(int argc, char *argv[]);

  private:
    void displayHelp();
    bool loadAnnotations(const std::string &annotationFilePath,
                         std::unordered_map<std::string, std::string> &annotations);
    void addAnnotations(std::istream &in, std::ostream &out,
                        const std::unordered_map<std::string, std::string> &annotations);

    // Generate a unique key for a variant
    // Format: CHROM:POS:REF:ALT
    static std::string generateVariantKey(const std::string &chrom, const std::string &pos, const std::string &ref,
                                          const std::string &alt);
};

// ---------------------------------------------------------------------------
// Utility: split a string by a delimiter (optimized with find-based approach)
// ---------------------------------------------------------------------------
static std::vector<std::string> split(const std::string &str, char delimiter) {
    std::vector<std::string> tokens;
    tokens.reserve(8);
    size_t start = 0;
    size_t end;
    while ((end = str.find(delimiter, start)) != std::string::npos) {
        tokens.emplace_back(str, start, end - start);
        start = end + 1;
    }
    tokens.emplace_back(str, start);
    return tokens;
}

// ---------------------------------------------------------------------------
// displayHelp
// ---------------------------------------------------------------------------
void VCFXCustomAnnotator::displayHelp() {
    std::cout << "VCFX_custom_annotator: Add custom annotations to the INFO field in a VCF file.\n\n"
              << "Usage:\n"
              << "  VCFX_custom_annotator --add-annotation <annotations.txt> [options]\n\n"
              << "Options:\n"
              << "  -h, --help                  Display this help message and exit\n"
              << "  -a, --add-annotation <file> Specify the annotation file\n\n"
              << "Description:\n"
              << "  Reads an annotation file with lines:\n"
              << "    CHROM  POS  REF  ALT  annotation...\n"
              << "  Then for each VCF variant, if it matches CHROM:POS:REF:ALT, inserts\n"
              << "  'CustomAnnotation=...' into the INFO field.\n"
              << "  Multi-allelic ALT fields are split on commas; we attempt to annotate\n"
              << "  each ALT separately. If no annotation is found for a given ALT, 'NA'\n"
              << "  is used for that allele's slot.\n\n"
              << "Example:\n"
              << "  VCFX_custom_annotator --add-annotation annotations.txt < input.vcf > annotated.vcf\n";
}

// ---------------------------------------------------------------------------
// generateVariantKey
// ---------------------------------------------------------------------------
std::string VCFXCustomAnnotator::generateVariantKey(const std::string &chrom, const std::string &pos,
                                                    const std::string &ref, const std::string &alt) {
    return chrom + ":" + pos + ":" + ref + ":" + alt;
}

// ---------------------------------------------------------------------------
// loadAnnotations
//   Reads a file with lines: CHROM POS REF ALT annotation...
//   Stores them in a map: "chrom:pos:ref:alt" -> annotation
// ---------------------------------------------------------------------------
bool VCFXCustomAnnotator::loadAnnotations(const std::string &annotationFilePath,
                                          std::unordered_map<std::string, std::string> &annotations) {
    std::ifstream infile(annotationFilePath);
    if (!infile.is_open()) {
        std::cerr << "Error: Unable to open annotation file " << annotationFilePath << "\n";
        return false;
    }

    std::string line;
    size_t line_num = 0;
    while (std::getline(infile, line)) {
        line_num++;
        if (line.empty() || line[0] == '#') {
            // skip comments or empty lines
            continue;
        }
        std::stringstream ss(line);
        std::string chrom, pos, ref, alt;
        if (!(ss >> chrom >> pos >> ref >> alt)) {
            std::cerr << "Warning: Skipping invalid annotation line " << line_num << ": " << line << "\n";
            continue;
        }
        // The rest of the line (after the first four fields) is the annotation text
        std::string annotation;
        if (!std::getline(ss, annotation)) {
            // It's okay if there's no extra text, but we treat it as empty
            annotation = "";
        }
        // Trim leading whitespace from annotation
        if (!annotation.empty()) {
            size_t startPos = annotation.find_first_not_of(" \t");
            if (startPos == std::string::npos) {
                // it's all whitespace
                annotation.clear();
            } else {
                annotation.erase(0, startPos);
            }
        }
        // Build key
        std::string key = generateVariantKey(chrom, pos, ref, alt);
        annotations[key] = annotation;
    }
    return true;
}

// ---------------------------------------------------------------------------
// addAnnotations
//   Reads the VCF from 'in', writes to 'out' a new header line
//   and appends "CustomAnnotation=..." to each variant's INFO if found.
//   For multi-allelic ALT, we do one lookup per allele, merging them into
//   a single string "val1,val2" if there are multiple alt alleles.
//   If no annotation is found for an alt, we use "NA" for that slot.
// ---------------------------------------------------------------------------
void VCFXCustomAnnotator::addAnnotations(std::istream &in, std::ostream &out,
                                         const std::unordered_map<std::string, std::string> &annotations) {
    bool infoHeaderInserted = false;
    std::string line;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        // If it's a header line
        if (line[0] == '#') {
            // If it's #CHROM and we haven't inserted the new INFO header yet, do so
            if (!infoHeaderInserted && line.rfind("#CHROM", 0) == 0) {
                out << "##INFO=<ID=CustomAnnotation,Number=.,Type=String,Description=\"Custom annotations added by "
                       "VCFX_custom_annotator (multi-allelic)\">\n";
                infoHeaderInserted = true;
            }
            out << line << "\n";
            continue;
        }

        // parse the 8 standard fields + possibly more
        // e.g. CHROM POS ID REF ALT QUAL FILTER INFO ...
        std::stringstream ss(line);
        std::string chrom, pos, id, ref, alt, qual, filter, info;
        if (!(ss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info)) {
            std::cerr << "Warning: Skipping invalid VCF line: " << line << "\n";
            continue;
        }

        // read the remainder of the line as one string (could include FORMAT, sample columns, etc.)
        std::string rest;
        if (std::getline(ss, rest)) {
            // ' rest ' includes the leading space if any
            // We can just keep it as is
        }

        // If ALT is multi-allelic, e.g. "A,C,G", we do a separate lookup for each allele
        auto altAlleles = split(alt, ',');

        // For each alt allele, build the key, see if it's in the annotations
        // If not found, store "NA"
        std::vector<std::string> annVals;
        annVals.reserve(altAlleles.size());
        bool anyFound = false;
        for (auto &a : altAlleles) {
            std::string key = generateVariantKey(chrom, pos, ref, a);
            auto it = annotations.find(key);
            if (it == annotations.end()) {
                annVals.push_back("NA");
            } else {
                annVals.push_back(it->second.empty() ? "NA" : it->second);
                anyFound = true;
            }
        }

        // Build a single annotation field if needed:
        // e.g. "val1,val2" or "NA,val2" etc.
        // If we want to show the user which alt is which, we do them in order.
        std::string finalAnn;
        {
            // join annVals with commas
            std::ostringstream oss;
            for (size_t i = 0; i < annVals.size(); ++i) {
                if (i > 0)
                    oss << ",";
                oss << annVals[i];
            }
            finalAnn = oss.str();
        }

        // Insert into INFO if we want it even if all are "NA"
        // Some users might prefer skipping if all alt are "NA".
        // We'll keep it consistent and always add it, to see that no annotation was found.
        // If the original info is ".", we replace it; else append ";"
        if (info == ".") {
            info = "CustomAnnotation=" + finalAnn;
        } else {
            info += ";CustomAnnotation=" + finalAnn;
        }

        // Reconstruct the VCF line
        // We print CHROM POS ID REF ALT QUAL FILTER INFO then the rest
        // e.g. leftover might be " FORMAT SAMPLE1 SAMPLE2..."
        out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\t" << qual << "\t" << filter
            << "\t" << info;

        // If there's anything left in 'rest', print it
        if (!rest.empty()) {
            out << rest; // includes leading space
        }
        out << "\n";
    }
}

// ---------------------------------------------------------------------------
// run() - parse arguments, load annotation map, apply to VCF
// ---------------------------------------------------------------------------
int VCFXCustomAnnotator::run(int argc, char *argv[]) {
    // Parse command-line arguments
    int opt;
    bool showHelp = false;
    std::string annotationFilePath;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'}, {"add-annotation", required_argument, 0, 'a'}, {0, 0, 0, 0}};

    while ((opt = getopt_long(argc, argv, "ha:", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
            break;
        case 'a':
            annotationFilePath = optarg;
            break;
        default:
            showHelp = true;
        }
    }

    if (showHelp || annotationFilePath.empty()) {
        displayHelp();
        return (showHelp ? 0 : 1);
    }

    // Load annotation map
    std::unordered_map<std::string, std::string> annotations;
    if (!loadAnnotations(annotationFilePath, annotations)) {
        std::cerr << "Error: Failed to load annotations from " << annotationFilePath << "\n";
        return 1;
    }

    // Annotate from stdin to stdout
    addAnnotations(std::cin, std::cout, annotations);
    return 0;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
static void show_help() {
    VCFXCustomAnnotator obj;
    char arg0[] = "VCFX_custom_annotator";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_custom_annotator", show_help))
        return 0;
    VCFXCustomAnnotator annotator;
    return annotator.run(argc, argv);
}
