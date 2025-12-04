#include "VCFX_dosage_calculator.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <getopt.h>
#include <iomanip>
#include <sstream>

// ---------------------------------------------------------------------------
// displayHelp: Prints usage information to standard output.
// ---------------------------------------------------------------------------
void VCFXDosageCalculator::displayHelp() {
    std::cout << "VCFX_dosage_calculator: Calculate genotype dosage for each variant in a VCF file.\n\n";
    std::cout << "Usage:\n";
    std::cout << "  VCFX_dosage_calculator [options] < input.vcf > dosage_output.txt\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help    Display this help message and exit\n\n";
    std::cout << "Description:\n";
    std::cout << "  For each variant in the input VCF, the tool computes the dosage for each sample\n";
    std::cout << "  based on the genotype (GT) field. Dosage is defined as the number of alternate\n";
    std::cout << "  alleles (i.e. each allele > 0 counts as 1). Thus:\n";
    std::cout << "    0/0  => dosage 0\n";
    std::cout << "    0/1  => dosage 1\n";
    std::cout << "    1/1  => dosage 2\n";
    std::cout << "    1/2  => dosage 2  (each alternate, regardless of numeric value, counts as 1)\n\n";
    std::cout << "Example:\n";
    std::cout << "  VCFX_dosage_calculator < input.vcf > dosage_output.txt\n";
}

// ---------------------------------------------------------------------------
// run: Parses command-line arguments and calls calculateDosage.
// ---------------------------------------------------------------------------
int VCFXDosageCalculator::run(int argc, char *argv[]) {
    int opt;
    bool showHelp = false;

    static struct option long_options[] = {{"help", no_argument, 0, 'h'}, {0, 0, 0, 0}};

    while ((opt = getopt_long(argc, argv, "h", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
            break;
        default:
            showHelp = true;
        }
    }
    if (showHelp) {
        displayHelp();
        return 0;
    }

    // Process VCF input from stdin and write dosage output to stdout.
    calculateDosage(std::cin, std::cout);

    return 0;
}

// ---------------------------------------------------------------------------
// calculateDosage: Processes the VCF line by line, computes dosage for each sample,
// and outputs a tab-delimited line per variant.
// ---------------------------------------------------------------------------
void VCFXDosageCalculator::calculateDosage(std::istream &in, std::ostream &out) {
    std::string line;
    bool headerParsed = false;
    std::vector<std::string> headerFields;

    // Print output header. The output columns are:
    // CHROM, POS, ID, REF, ALT, Dosages
    // where Dosages is a comma-separated list corresponding to each sample.
    out << "CHROM\tPOS\tID\tREF\tALT\tDosages\n";

    while (std::getline(in, line)) {
        if (line.empty())
            continue;

        // Process header lines
        if (line[0] == '#') {
            // Parse header line with "#CHROM" to get sample names.
            if (line.rfind("#CHROM", 0) == 0) {
                std::stringstream ss(line);
                std::string field;
                while (std::getline(ss, field, '\t')) {
                    headerFields.push_back(field);
                }
                headerParsed = true;
            }
            // Skip printing header lines in output.
            continue;
        }

        if (!headerParsed) {
            std::cerr << "Error: VCF header (#CHROM) not found before variant records.\n";
            return;
        }

        // Split the variant line into fields (VCF standard: at least 10 fields expected)
        std::stringstream ss(line);
        std::vector<std::string> fields;
        std::string field;
        while (std::getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        // We expect at least 10 columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, and at least one sample.
        if (fields.size() < 10) {
            std::cerr << "Warning: Skipping VCF line with fewer than 10 fields: " << line << "\n";
            continue;
        }

        // Extract variant-level fields.
        std::string chrom = fields[0];
        std::string pos = fields[1];
        std::string id = fields[2];
        std::string ref = fields[3];
        std::string alt = fields[4];
        std::string format = fields[8];

        // For dosage calculation, we need to find the genotype (GT) field index.
        std::vector<std::string> formatFields = split(format, ':');
        int gtIndex = -1;
        for (size_t i = 0; i < formatFields.size(); ++i) {
            if (formatFields[i] == "GT") {
                gtIndex = static_cast<int>(i);
                break;
            }
        }
        if (gtIndex == -1) {
            std::cerr << "Warning: GT field not found in FORMAT column for variant " << chrom << ":" << pos << ".\n";
            out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\tNA\n";
            continue;
        }

        // Process each sample column (from index 9 onward)
        std::vector<std::string> dosageValues;
        for (size_t i = 9; i < fields.size(); ++i) {
            std::string sampleField = fields[i];
            // Split sample field by ':'; the first field is typically the genotype.
            std::vector<std::string> sampleParts = split(sampleField, ':');
            if (sampleParts.empty() || gtIndex >= static_cast<int>(sampleParts.size())) {
                dosageValues.push_back("NA");
                continue;
            }
            std::string genotype = sampleParts[gtIndex];
            // Normalize genotype separators: replace '|' with '/'
            std::replace(genotype.begin(), genotype.end(), '|', '/');

            if (genotype.empty() || genotype == ".") {
                dosageValues.push_back("NA");
                continue;
            }

            // Split genotype on '/'
            std::vector<std::string> alleles = split(genotype, '/');
            if (alleles.size() != 2) {
                // Non-diploid genotype or invalid format; mark as missing.
                dosageValues.push_back("NA");
                continue;
            }

            bool validGenotype = true;
            int dosage = 0;
            // For each allele, if it is "0", count 0; if it is any numeric >0, count 1.
            for (const std::string &a : alleles) {
                if (a == "." || a.empty()) {
                    validGenotype = false;
                    break;
                }
                // Ensure that the allele string is numeric.
                if (!std::all_of(a.begin(), a.end(), ::isdigit)) {
                    validGenotype = false;
                    break;
                }
                int alleleVal = std::stoi(a);
                // Count alternate alleles as 1 regardless of allele index.
                if (alleleVal > 0) {
                    dosage += 1;
                }
            }
            if (validGenotype) {
                dosageValues.push_back(std::to_string(dosage));
            } else {
                dosageValues.push_back("NA");
            }
        }

        // Join dosage values into a comma-separated string.
        std::ostringstream dosagesOSS;
        for (size_t i = 0; i < dosageValues.size(); ++i) {
            if (i > 0) {
                dosagesOSS << ",";
            }
            dosagesOSS << dosageValues[i];
        }

        // Output one line per variant with dosage information.
        out << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\t" << dosagesOSS.str() << "\n";
    }
}

// ---------------------------------------------------------------------------
// split: Helper function to split a string by a delimiter
// ---------------------------------------------------------------------------
std::vector<std::string> VCFXDosageCalculator::split(const std::string &str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

static void show_help() {
    VCFXDosageCalculator obj;
    char arg0[] = "VCFX_dosage_calculator";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_dosage_calculator", show_help))
        return 0;
    VCFXDosageCalculator dosageCalculator;
    return dosageCalculator.run(argc, argv);
}
