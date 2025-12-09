#include "VCFX_variant_classifier.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <cstring>
#include <getopt.h>
#include <iostream>
#include <poll.h>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

// mmap support
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

// ============================================================================
// SIMD newline scanning
// ============================================================================

#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#define HAVE_AVX2 1
#endif

#if defined(__aarch64__) || defined(_M_ARM64)
#include <arm_neon.h>
#define HAVE_NEON 1
#endif

// Fast newline finder using memchr (often SIMD-optimized by libc)
static inline const char* findNewline(const char* start, const char* end) {
    return static_cast<const char*>(memchr(start, '\n', end - start));
}

// Fast tab finder
static inline const char* findTab(const char* start, const char* end) {
    return static_cast<const char*>(memchr(start, '\t', end - start));
}

// Fast comma finder
static inline const char* findComma(const char* start, const char* end) {
    return static_cast<const char*>(memchr(start, ',', end - start));
}

// ============================================================================
// Helper functions
// ============================================================================

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

// Check if a character is alphabetic (for validation)
static inline bool isAlphaFast(char c) {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
}

// Check if a character is a digit
static inline bool isDigitFast(char c) {
    return c >= '0' && c <= '9';
}

// ============================================================================
// VCFXVariantClassifier implementation
// ============================================================================

int VCFXVariantClassifier::run(int argc, char *argv[]) {
    bool showHelp = false;
    const char* inputFile = nullptr;

    // Check if stdin has data available when no arguments
    if (argc == 1) {
        struct pollfd fds;
        fds.fd = 0;
        fds.events = POLLIN;
        int ret = poll(&fds, 1, 0);
        if (ret <= 0 || !(fds.revents & POLLIN)) {
            displayHelp();
            return 0;
        }
    }

    static struct option long_opts[] = {
        {"help", no_argument, 0, 'h'},
        {"append-info", no_argument, 0, 'a'},
        {"input", required_argument, 0, 'i'},
        {"quiet", no_argument, 0, 'q'},
        {0, 0, 0, 0}
    };

    optind = 1;
    while (true) {
        int c = ::getopt_long(argc, argv, "hai:q", long_opts, nullptr);
        if (c == -1)
            break;
        switch (c) {
        case 'h':
            showHelp = true;
            break;
        case 'a':
            appendInfo = true;
            break;
        case 'i':
            inputFile = optarg;
            break;
        case 'q':
            quietMode = true;
            break;
        default:
            showHelp = true;
        }
    }

    // Check for positional argument
    if (!inputFile && optind < argc) {
        inputFile = argv[optind];
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    if (inputFile) {
        if (!processFileMmap(inputFile, std::cout)) {
            return 1;
        }
    } else {
        classifyStream(std::cin, std::cout);
    }
    return 0;
}

void VCFXVariantClassifier::displayHelp() {
    std::cout << "VCFX_variant_classifier: Classify variants in a VCF as SNP, INDEL, MNV, or STRUCTURAL.\n\n"
                 "Usage:\n"
                 "  VCFX_variant_classifier [options] [input.vcf]\n"
                 "  VCFX_variant_classifier [options] < input.vcf > output.vcf_or_tsv\n\n"
                 "Options:\n"
                 "  -h, --help         Show help.\n"
                 "  -i, --input FILE   Input VCF file (uses mmap for best performance).\n"
                 "  -a, --append-info  Instead of producing a TSV, output a valid VCF\n"
                 "                     with a new 'VCF_CLASS' subfield in the INFO.\n"
                 "  -q, --quiet        Suppress warnings to stderr.\n\n"
                 "Description:\n"
                 "  Reads each variant line, determines if it is:\n"
                 "    SNP: single base ref & alt,\n"
                 "    INDEL: length mismatch (less than 50 bp difference) in ref vs alt,\n"
                 "    MNV: same length >1,\n"
                 "    STRUCTURAL: alt is symbolic (<DEL>, <INS>, <DUP>), or breakend ([chr etc.)\n"
                 "                or length difference >=50.\n"
                 "  If --append-info, prints original columns + updated INFO. Otherwise prints\n"
                 "  'CHROM POS ID REF ALT Classification' as TSV.\n\n"
                 "Performance:\n"
                 "  When using -i/--input, the tool uses memory-mapped I/O for\n"
                 "  ~10-20x faster processing of large files.\n\n"
                 "Examples:\n"
                 "  1) TSV classification:\n"
                 "     VCFX_variant_classifier < input.vcf > classified.tsv\n"
                 "  2) Modify INFO in output VCF:\n"
                 "     VCFX_variant_classifier --append-info < input.vcf > annotated.vcf\n"
                 "  3) Fast file mode:\n"
                 "     VCFX_variant_classifier -i input.vcf > classified.tsv\n";
}

std::vector<std::string> VCFXVariantClassifier::split(const std::string &s, char delim) const {
    std::vector<std::string> out;
    out.reserve(16);
    size_t start = 0;
    size_t end;
    while ((end = s.find(delim, start)) != std::string::npos) {
        out.emplace_back(s, start, end - start);
        start = end + 1;
    }
    if (start <= s.size()) {
        out.emplace_back(s, start);
    }
    return out;
}

bool VCFXVariantClassifier::isStructuralAllele(const std::string &alt) const {
    if (!alt.empty() && alt.front() == '<' && alt.back() == '>')
        return true;
    if (alt.find('[') != std::string::npos || alt.find(']') != std::string::npos) {
        return true;
    }
    return false;
}

bool VCFXVariantClassifier::isStructuralAlleleSV(std::string_view alt) const {
    if (alt.empty()) return false;
    // Symbolic: <DEL>, <INS>, etc.
    if (alt.front() == '<' && alt.back() == '>')
        return true;
    // Breakend notation
    for (char c : alt) {
        if (c == '[' || c == ']') return true;
    }
    return false;
}

VariantType VCFXVariantClassifier::classifyAllele(const std::string &ref, const std::string &alt) const {
    if (isStructuralAllele(alt)) {
        return VariantType::STRUCTURAL;
    }
    if (std::abs((int)ref.size() - (int)alt.size()) >= 50) {
        return VariantType::STRUCTURAL;
    }
    if (ref == alt) {
        return VariantType::UNKNOWN;
    }
    if (ref.size() == 1 && alt.size() == 1 && std::isalpha((unsigned char)ref[0]) &&
        std::isalpha((unsigned char)alt[0])) {
        return VariantType::SNP;
    }
    if (ref.size() != alt.size()) {
        if (ref.size() >= 40 || alt.size() >= 40) {
            return VariantType::STRUCTURAL;
        }
        return VariantType::INDEL;
    }
    if (ref.size() > 1) {
        return VariantType::MNV;
    }
    return VariantType::UNKNOWN;
}

VariantType VCFXVariantClassifier::classifyAlleleSV(std::string_view ref, std::string_view alt) const {
    if (isStructuralAlleleSV(alt)) {
        return VariantType::STRUCTURAL;
    }
    int lenDiff = static_cast<int>(ref.size()) - static_cast<int>(alt.size());
    if (lenDiff < 0) lenDiff = -lenDiff;
    if (lenDiff >= 50) {
        return VariantType::STRUCTURAL;
    }
    if (ref == alt) {
        return VariantType::UNKNOWN;
    }
    if (ref.size() == 1 && alt.size() == 1 && isAlphaFast(ref[0]) && isAlphaFast(alt[0])) {
        return VariantType::SNP;
    }
    if (ref.size() != alt.size()) {
        if (ref.size() >= 40 || alt.size() >= 40) {
            return VariantType::STRUCTURAL;
        }
        return VariantType::INDEL;
    }
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

const char* VCFXVariantClassifier::typeToStrFast(VariantType t) const {
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

std::string VCFXVariantClassifier::appendClassification(const std::string &line) {
    auto fields = split(line, '\t');
    if (fields.size() < 8) {
        return line;
    }
    auto alts = split(fields[4], ',');
    VariantType vt = classifyVariant(fields[3], alts);
    std::string cstr = typeToStr(vt);

    if (fields[7] == "." || fields[7].empty()) {
        fields[7] = "VCF_CLASS=" + cstr;
    } else {
        if (!fields[7].empty() && fields[7].back() != ';') {
            fields[7] += ";";
        }
        fields[7] += "VCF_CLASS=" + cstr;
    }
    std::string result;
    result.reserve(line.size() + 32);
    for (size_t i = 0; i < fields.size(); i++) {
        if (i > 0)
            result += '\t';
        result += fields[i];
    }
    return result;
}

// ============================================================================
// Optimized mmap-based processing
// ============================================================================

bool VCFXVariantClassifier::processFileMmap(const char* filename, std::ostream &out) {
    // Open file
    int fd = open(filename, O_RDONLY);
    if (fd < 0) {
        std::cerr << "Error: Cannot open file: " << filename << "\n";
        return false;
    }

    // Get file size
    struct stat st;
    if (fstat(fd, &st) < 0) {
        std::cerr << "Error: Cannot stat file: " << filename << "\n";
        close(fd);
        return false;
    }

    size_t fileSize = st.st_size;
    if (fileSize == 0) {
        close(fd);
        return true;
    }

    // Memory map
    const char* data = static_cast<const char*>(mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, fd, 0));
    if (data == MAP_FAILED) {
        std::cerr << "Error: Cannot mmap file: " << filename << "\n";
        close(fd);
        return false;
    }

    madvise(const_cast<char*>(data), fileSize, MADV_SEQUENTIAL);

    const char* ptr = data;
    const char* dataEnd = data + fileSize;

    // Output buffer (1MB)
    std::string outputBuffer;
    outputBuffer.reserve(1024 * 1024);
    const size_t flushThreshold = 900 * 1024;

    bool headerFound = false;

    // Output header for TSV mode
    if (!appendInfo) {
        outputBuffer.append("CHROM\tPOS\tID\tREF\tALT\tClassification\n");
    }

    // Process lines
    while (ptr < dataEnd) {
        const char* lineEnd = findNewline(ptr, dataEnd);
        if (!lineEnd) lineEnd = dataEnd;

        size_t lineLen = lineEnd - ptr;
        if (lineLen == 0) {
            ptr = lineEnd + 1;
            continue;
        }

        // Handle header lines
        if (ptr[0] == '#') {
            if (lineLen >= 6 && strncmp(ptr, "#CHROM", 6) == 0) {
                headerFound = true;
            }
            if (appendInfo) {
                // Pass through header lines in append mode
                outputBuffer.append(ptr, lineLen);
                outputBuffer.push_back('\n');
            }
            ptr = lineEnd + 1;
            continue;
        }

        if (!headerFound) {
            if (!quietMode) {
                std::cerr << "Warning: data line before #CHROM => skipping.\n";
            }
            ptr = lineEnd + 1;
            continue;
        }

        // Parse fields (need CHROM, POS, ID, REF, ALT, and optionally more for append mode)
        // Fields: 0=CHROM, 1=POS, 2=ID, 3=REF, 4=ALT, 5=QUAL, 6=FILTER, 7=INFO
        const char* fieldStarts[10];
        size_t fieldLens[10];
        int numFields = 0;

        const char* fp = ptr;
        while (fp < lineEnd && numFields < 10) {
            fieldStarts[numFields] = fp;
            const char* tabPos = findTab(fp, lineEnd);
            if (!tabPos) tabPos = lineEnd;
            fieldLens[numFields] = tabPos - fp;
            numFields++;
            fp = (tabPos < lineEnd) ? tabPos + 1 : lineEnd;
        }

        if (numFields < 5) {
            ptr = lineEnd + 1;
            continue;
        }

        // Extract fields as string_view for zero-copy
        std::string_view chrom(fieldStarts[0], fieldLens[0]);
        std::string_view pos(fieldStarts[1], fieldLens[1]);
        std::string_view id(fieldStarts[2], fieldLens[2]);
        std::string_view ref(fieldStarts[3], fieldLens[3]);
        std::string_view alt(fieldStarts[4], fieldLens[4]);

        // Validate CHROM
        if (chrom.empty()) {
            if (!quietMode) {
                std::cerr << "Warning: empty chromosome field => skipping.\n";
            }
            ptr = lineEnd + 1;
            continue;
        }

        // Validate POS is numeric
        bool validPos = true;
        for (char c : pos) {
            if (!isDigitFast(c)) {
                validPos = false;
                break;
            }
        }
        if (!validPos || pos.empty()) {
            if (!quietMode) {
                std::cerr << "Warning: position is not numeric => skipping.\n";
            }
            ptr = lineEnd + 1;
            continue;
        }

        // Validate REF and ALT
        if (ref.empty() || alt.empty()) {
            if (!quietMode) {
                std::cerr << "Warning: REF or ALT is empty => skipping.\n";
            }
            ptr = lineEnd + 1;
            continue;
        }

        // Validate REF contains only valid bases
        bool validRef = true;
        for (char c : ref) {
            if (!isAlphaFast(c)) {
                validRef = false;
                break;
            }
        }
        if (!validRef) {
            if (!quietMode) {
                std::cerr << "Warning: REF contains non-alphabetic characters => skipping.\n";
            }
            ptr = lineEnd + 1;
            continue;
        }

        // Validate ALT doesn't end with comma
        if (alt.back() == ',') {
            if (!quietMode) {
                std::cerr << "Warning: ALT ends with a comma => skipping.\n";
            }
            ptr = lineEnd + 1;
            continue;
        }

        // Classify variant (handling multi-allelic)
        VariantType finalType = VariantType::UNKNOWN;
        bool anyStruct = false, anyMNV = false, anyIndel = false, anySNP = false;

        // Parse ALT alleles (comma-separated)
        const char* altStart = alt.data();
        const char* altEnd = alt.data() + alt.size();
        const char* ap = altStart;

        while (ap < altEnd) {
            const char* commaPos = findComma(ap, altEnd);
            if (!commaPos) commaPos = altEnd;

            std::string_view oneAlt(ap, commaPos - ap);
            VariantType vt = classifyAlleleSV(ref, oneAlt);

            if (vt == VariantType::STRUCTURAL) anyStruct = true;
            else if (vt == VariantType::MNV) anyMNV = true;
            else if (vt == VariantType::INDEL) anyIndel = true;
            else if (vt == VariantType::SNP) anySNP = true;

            ap = (commaPos < altEnd) ? commaPos + 1 : altEnd;
        }

        // Priority: STRUCTURAL > MNV > INDEL > SNP > UNKNOWN
        if (anyStruct) finalType = VariantType::STRUCTURAL;
        else if (anyMNV) finalType = VariantType::MNV;
        else if (anyIndel) finalType = VariantType::INDEL;
        else if (anySNP) finalType = VariantType::SNP;

        const char* typeStr = typeToStrFast(finalType);

        if (appendInfo) {
            // Append VCF_CLASS to INFO field
            if (numFields < 8) {
                // Not enough fields, just output original
                outputBuffer.append(ptr, lineLen);
                outputBuffer.push_back('\n');
            } else {
                // Output fields 0-6
                for (int f = 0; f < 7; f++) {
                    outputBuffer.append(fieldStarts[f], fieldLens[f]);
                    outputBuffer.push_back('\t');
                }
                // Modify INFO (field 7)
                std::string_view info(fieldStarts[7], fieldLens[7]);
                if (info == "." || info.empty()) {
                    outputBuffer.append("VCF_CLASS=");
                    outputBuffer.append(typeStr);
                } else {
                    outputBuffer.append(info);
                    if (info.back() != ';') {
                        outputBuffer.push_back(';');
                    }
                    outputBuffer.append("VCF_CLASS=");
                    outputBuffer.append(typeStr);
                }
                // Output remaining fields (8+)
                if (numFields > 8) {
                    outputBuffer.push_back('\t');
                    // Find position after field 7
                    const char* rest = fieldStarts[8];
                    size_t restLen = lineEnd - rest;
                    outputBuffer.append(rest, restLen);
                }
                outputBuffer.push_back('\n');
            }
        } else {
            // TSV output: CHROM POS ID REF ALT Classification
            outputBuffer.append(chrom);
            outputBuffer.push_back('\t');
            outputBuffer.append(pos);
            outputBuffer.push_back('\t');
            outputBuffer.append(id);
            outputBuffer.push_back('\t');
            outputBuffer.append(ref);
            outputBuffer.push_back('\t');
            outputBuffer.append(alt);
            outputBuffer.push_back('\t');
            outputBuffer.append(typeStr);
            outputBuffer.push_back('\n');
        }

        ptr = lineEnd + 1;

        // Periodic flush
        if (outputBuffer.size() >= flushThreshold) {
            out.write(outputBuffer.data(), outputBuffer.size());
            outputBuffer.clear();
        }
    }

    // Final flush
    if (!outputBuffer.empty()) {
        out.write(outputBuffer.data(), outputBuffer.size());
    }

    // Cleanup
    munmap(const_cast<char*>(data), fileSize);
    close(fd);

    return true;
}

void VCFXVariantClassifier::classifyStream(std::istream &in, std::ostream &out) {
    // Add I/O buffering for performance
    constexpr size_t BUFFER_SIZE = 1 << 20;  // 1MB
    std::vector<char> inBuffer(BUFFER_SIZE);
    std::vector<char> outBuffer(BUFFER_SIZE);
    in.rdbuf()->pubsetbuf(inBuffer.data(), BUFFER_SIZE);
    out.rdbuf()->pubsetbuf(outBuffer.data(), BUFFER_SIZE);

    bool foundChromHeader = false;
    std::string line;
    line.reserve(65536);
    if (appendInfo) {
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
            if (!foundChromHeader) {
                if (!quietMode) {
                    std::cerr << "Warning: data line encountered before #CHROM => skipping.\n";
                }
                continue;
            }
            auto fields = split(line, '\t');
            if (fields.size() < 8) {
                if (!quietMode) {
                    std::cerr << "Warning: skipping line <8 columns.\n";
                }
                continue;
            }
            std::string newLine = appendClassification(line);
            out << newLine << "\n";
        }
    } else {
        out << "CHROM\tPOS\tID\tREF\tALT\tClassification\n";

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
                if (!quietMode) {
                    std::cerr << "Warning: data line before #CHROM => skipping.\n";
                }
                continue;
            }
            auto fields = split(line, '\t');
            if (fields.size() < 8) {
                if (!quietMode) {
                    std::cerr << "Warning: skipping line <8 columns.\n";
                }
                continue;
            }

            if (fields[0].empty()) {
                if (!quietMode) {
                    std::cerr << "Warning: empty chromosome field => skipping.\n";
                }
                continue;
            }

            bool validPos = true;
            for (char c : fields[1]) {
                if (!std::isdigit(c)) {
                    validPos = false;
                    break;
                }
            }
            if (!validPos) {
                if (!quietMode) {
                    std::cerr << "Warning: position is not numeric => skipping.\n";
                }
                continue;
            }

            if (fields[3].empty() || fields[4].empty()) {
                if (!quietMode) {
                    std::cerr << "Warning: REF or ALT is empty => skipping.\n";
                }
                continue;
            }

            bool validRef = true;
            for (char c : fields[3]) {
                if (!std::isalpha(c)) {
                    validRef = false;
                    break;
                }
            }
            if (!validRef) {
                if (!quietMode) {
                    std::cerr << "Warning: REF contains non-alphabetic characters => skipping.\n";
                }
                continue;
            }

            if (fields[4].back() == ',') {
                if (!quietMode) {
                    std::cerr << "Warning: ALT ends with a comma => skipping.\n";
                }
                continue;
            }

            auto altList = split(fields[4], ',');
            VariantType vt = classifyVariant(fields[3], altList);
            std::string altJoined = fields[4];
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
    vcfx::init_io();
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);

    if (vcfx::handle_common_flags(argc, argv, "VCFX_variant_classifier", show_help))
        return 0;
    VCFXVariantClassifier app;
    return app.run(argc, argv);
}
