#include "VCFX_variant_classifier.h"
#include "vcfx_core.h"
#include <algorithm>
#include <cctype>
#include <getopt.h>
#include <iostream>
#include <poll.h> // Add this for poll() function
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

// A helper to trim spaces
static std::string trim(const std::string &s) {
    size_t start = 0;
    while (start < s.size() && std::isspace((unsigned char)s[start]))
        start++;
    if (start == s.size())
        return "";
    size_t end = s.size() - 1;
    while (end > start && std::isspace((unsigned char)s[end]))
        end--;
    return s.substr(start, end - start + 1);
}

int VCFXVariantClassifier::run(int argc, char *argv[]) {
    bool showHelp = false;

    // Check if stdin has data available when no arguments
    if (argc == 1) {
        // Check if stdin has data using poll
        struct pollfd fds;
        fds.fd = 0; // stdin
        fds.events = POLLIN;

        // Poll stdin with 0 timeout to check if data is available
        int ret = poll(&fds, 1, 0);

        // Show help if no data is available on stdin
        if (ret <= 0 || !(fds.revents & POLLIN)) {
            displayHelp();
            return 0;
        }
    }

    static struct option long_opts[] = {
        {"help", no_argument, 0, 'h'}, {"append-info", no_argument, 0, 'a'}, {0, 0, 0, 0}};
    while (true) {
        int c = ::getopt_long(argc, argv, "ha", long_opts, nullptr);
        if (c == -1)
            break;
        switch (c) {
        case 'h':
            showHelp = true;
            break;
        case 'a':
            appendInfo = true;
            break;
        default:
            showHelp = true;
        }
    }
    if (showHelp) {
        displayHelp();
        return 0;
    }
    classifyStream(std::cin, std::cout);
    return 0;
}

void VCFXVariantClassifier::displayHelp() {
    std::cout << "VCFX_variant_classifier: Classify variants in a VCF as SNP, INDEL, MNV, or STRUCTURAL.\n\n"
                 "Usage:\n"
                 "  VCFX_variant_classifier [options] < input.vcf > output.vcf_or_tsv\n\n"
                 "Options:\n"
                 "  -h, --help         Show help.\n"
                 "  -a, --append-info  Instead of producing a TSV, output a valid VCF\n"
                 "                     with a new 'VCF_CLASS' subfield in the INFO.\n\n"
                 "Description:\n"
                 "  Reads each variant line, determines if it is:\n"
                 "    SNP: single base ref & alt,\n"
                 "    INDEL: length mismatch (less than 50 bp difference) in ref vs alt,\n"
                 "    MNV: same length >1,\n"
                 "    STRUCTURAL: alt is symbolic (<DEL>, <INS>, <DUP>), or breakend ([chr etc.)\n"
                 "                or length difference >=50.\n"
                 "  If --append-info, prints original columns + updated INFO. Otherwise prints\n"
                 "  'CHROM POS ID REF ALT Classification' as TSV.\n\n"
                 "Examples:\n"
                 "  1) TSV classification:\n"
                 "     VCFX_variant_classifier < input.vcf > classified.tsv\n"
                 "  2) Modify INFO in output VCF:\n"
                 "     VCFX_variant_classifier --append-info < input.vcf > annotated.vcf\n";
}

std::vector<std::string> VCFXVariantClassifier::split(const std::string &s, char delim) const {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string t;
    while (std::getline(ss, t, delim)) {
        out.push_back(t);
    }
    return out;
}

bool VCFXVariantClassifier::isStructuralAllele(const std::string &alt) const {
    // If alt is symbolic: <DEL>, <INS>, ...
    if (!alt.empty() && alt.front() == '<' && alt.back() == '>')
        return true;
    // If alt has breakend notation: [chr..., ]chr...
    if (alt.find('[') != std::string::npos || alt.find(']') != std::string::npos) {
        return true;
    }
    // If alt length or ref length > 50 => structural
    // We'll handle that in classifyAllele by comparing lengths
    return false;
}

VariantType VCFXVariantClassifier::classifyAllele(const std::string &ref, const std::string &alt) const {
    // check structural
    if (isStructuralAllele(alt)) {
        return VariantType::STRUCTURAL;
    }

    // if length difference >=50 => structural
    if (std::abs((int)ref.size() - (int)alt.size()) >= 50) {
        return VariantType::STRUCTURAL;
    }

    // If REF and ALT are identical, this should be UNKNOWN
    if (ref == alt) {
        return VariantType::UNKNOWN;
    }

    // if both single base
    if (ref.size() == 1 && alt.size() == 1 && std::isalpha((unsigned char)ref[0]) &&
        std::isalpha((unsigned char)alt[0])) {
        return VariantType::SNP;
    }

    // if length differs => small INDEL
    if (ref.size() != alt.size()) {
        // Additional check for large variants
        if (ref.size() >= 40 || alt.size() >= 40) {
            return VariantType::STRUCTURAL;
        }
        return VariantType::INDEL;
    }

    // same length >1 => MNV
    if (ref.size() > 1) {
        return VariantType::MNV;
    }
    return VariantType::UNKNOWN;
}

VariantType VCFXVariantClassifier::classifyVariant(const std::string &ref, const std::vector<std::string> &alts) const {
    bool anyStruct = false, anyMNV = false, anyIndel = false, anySNP = false;
    for (auto &a : alts) {
        VariantType vt = classifyAllele(ref, a);
        if (vt == VariantType::STRUCTURAL)
            anyStruct = true;
        else if (vt == VariantType::MNV)
            anyMNV = true;
        else if (vt == VariantType::INDEL)
            anyIndel = true;
        else if (vt == VariantType::SNP)
            anySNP = true;
    }
    // priority: STRUCTURAL> MNV> INDEL> SNP> UNKNOWN
    if (anyStruct)
        return VariantType::STRUCTURAL;
    if (anyMNV)
        return VariantType::MNV;
    if (anyIndel)
        return VariantType::INDEL;
    if (anySNP)
        return VariantType::SNP;
    return VariantType::UNKNOWN;
}

std::string VCFXVariantClassifier::typeToStr(VariantType t) const {
    switch (t) {
    case VariantType::SNP:
        return "SNP";
    case VariantType::INDEL:
        return "INDEL";
    case VariantType::MNV:
        return "MNV";
    case VariantType::STRUCTURAL:
        return "STRUCTURAL";
    default:
        return "UNKNOWN";
    }
}

// If we want to keep the line as a valid VCF, we parse the line, classify, then append to INFO "VCF_CLASS=XYZ"
std::string VCFXVariantClassifier::appendClassification(const std::string &line) {
    // split by tab
    auto fields = split(line, '\t');
    if (fields.size() < 8) {
        // skip
        return line;
    }
    // parse alt
    auto alts = split(fields[4], ',');
    VariantType vt = classifyVariant(fields[3], alts);
    std::string cstr = typeToStr(vt);

    // append to INFO
    if (fields[7] == "." || fields[7].empty()) {
        fields[7] = "VCF_CLASS=" + cstr;
    } else {
        // add semicolon if needed
        if (!fields[7].empty() && fields[7].back() != ';') {
            fields[7] += ";";
        }
        fields[7] += "VCF_CLASS=" + cstr;
    }
    // reconstruct line
    std::ostringstream oss;
    for (size_t i = 0; i < fields.size(); i++) {
        if (i > 0)
            oss << "\t";
        oss << fields[i];
    }
    return oss.str();
}

void VCFXVariantClassifier::classifyStream(std::istream &in, std::ostream &out) {
    bool foundChromHeader = false;
    std::string line;
    if (appendInfo) {
        // we produce a valid VCF. We'll pass all # lines unmodified
        while (true) {
            if (!std::getline(in, line))
                break;
            if (line.empty()) {
                out << line << "\n";
                continue;
            }
            if (line[0] == '#') {
                out << line << "\n";
                if (line.rfind("#CHROM", 0) == 0)
                    foundChromHeader = true;
                continue;
            }
            // data line
            if (!foundChromHeader) {
                std::cerr << "Warning: data line encountered before #CHROM => skipping.\n";
                continue;
            }
            // parse, classify, re-output
            auto fields = split(line, '\t');
            if (fields.size() < 8) {
                std::cerr << "Warning: skipping line <8 columns.\n";
                continue;
            }
            // append classification in INFO
            std::string newLine = appendClassification(line);
            out << newLine << "\n";
        }
    } else {
        // produce a TSV: CHROM POS ID REF ALT CLASS
        // Output the TSV header first
        out << "CHROM\tPOS\tID\tREF\tALT\tClassification\n";

        // skip all # lines
        while (true) {
            if (!std::getline(in, line))
                break;
            if (line.empty())
                continue;
            if (line[0] == '#') {
                if (line.rfind("#CHROM", 0) == 0)
                    foundChromHeader = true;
                continue;
            }
            if (!foundChromHeader) {
                std::cerr << "Warning: data line before #CHROM => skipping.\n";
                continue;
            }
            // parse
            auto fields = split(line, '\t');
            if (fields.size() < 8) {
                std::cerr << "Warning: skipping line <8 columns.\n";
                continue;
            }

            // Validate CHROM is not empty (accept any chromosome format: chr20, 20, chrM, MT, etc.)
            if (fields[0].empty()) {
                std::cerr << "Warning: empty chromosome field => skipping.\n";
                continue;
            }

            // 2. Validate POS is numeric
            bool validPos = true;
            for (char c : fields[1]) {
                if (!std::isdigit(c)) {
                    validPos = false;
                    break;
                }
            }
            if (!validPos) {
                std::cerr << "Warning: position is not numeric => skipping.\n";
                continue;
            }

            // 3. Validate REF and ALT are not empty
            if (fields[3].empty() || fields[4].empty()) {
                std::cerr << "Warning: REF or ALT is empty => skipping.\n";
                continue;
            }

            // 4. Validate REF contains only valid bases
            bool validRef = true;
            for (char c : fields[3]) {
                if (!std::isalpha(c)) {
                    validRef = false;
                    break;
                }
            }
            if (!validRef) {
                std::cerr << "Warning: REF contains non-alphabetic characters => skipping.\n";
                continue;
            }

            // 5. Validate ALT format
            // Allow special cases like <DEL>, but disallow trailing commas
            if (fields[4].back() == ',') {
                std::cerr << "Warning: ALT ends with a comma => skipping.\n";
                continue;
            }

            // alt => split by comma
            auto altList = split(fields[4], ',');
            VariantType vt = classifyVariant(fields[3], altList);
            // build alt string
            std::string altJoined = fields[4];
            // Output the classification line
            out << fields[0] << "\t" << fields[1] << "\t" << fields[2] << "\t" << fields[3] << "\t" << altJoined << "\t"
                << typeToStr(vt) << "\n";
        }
    }
}

static void show_help() {
    VCFXVariantClassifier obj;
    char arg0[] = "VCFX_variant_classifier";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    if (vcfx::handle_common_flags(argc, argv, "VCFX_variant_classifier", show_help))
        return 0;
    VCFXVariantClassifier app;
    return app.run(argc, argv);
}
