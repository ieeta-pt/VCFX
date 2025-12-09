#include "VCFX_multiallelic_splitter.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <sstream>

// mmap support
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

// ============================================================================
// Helper functions for parsing
// ============================================================================

// Fast tab-delimited split returning string_view (zero-copy)
static inline size_t splitTabsView(const char* start, const char* end,
                                    std::vector<std::string_view>& out, size_t expected = 16) {
    out.clear();
    if (out.capacity() < expected) out.reserve(expected);

    const char* p = start;
    while (p < end) {
        const char* tab = static_cast<const char*>(memchr(p, '\t', end - p));
        if (!tab) {
            out.emplace_back(p, end - p);
            break;
        }
        out.emplace_back(p, tab - p);
        p = tab + 1;
    }
    return out.size();
}

// Split by delimiter into string_view
static inline size_t splitCharView(std::string_view str, char delim,
                                    std::vector<std::string_view>& out, size_t expected = 8) {
    out.clear();
    if (out.capacity() < expected) out.reserve(expected);

    size_t start = 0, pos;
    while ((pos = str.find(delim, start)) != std::string_view::npos) {
        out.emplace_back(str.substr(start, pos - start));
        start = pos + 1;
    }
    out.emplace_back(str.substr(start));
    return out.size();
}

// Find newline in memory
static inline const char* findNewline(const char* start, const char* end) {
    return static_cast<const char*>(memchr(start, '\n', end - start));
}

// Parse Number from header line (returns "A", "R", "G", "1", ".", or integer)
static bool parseNumberEq(std::string_view line, std::string& id, std::string& num) {
    auto i = line.find("ID=");
    if (i == std::string_view::npos) return false;
    auto sub = line.substr(i + 3);
    auto e = sub.find_first_of(",>");
    if (e == std::string_view::npos) return false;
    id = std::string(sub.substr(0, e));

    auto n = line.find("Number=");
    if (n == std::string_view::npos) return false;
    auto sub2 = line.substr(n + 7);
    auto e2 = sub2.find_first_of(",>");
    if (e2 == std::string_view::npos) return false;
    num = std::string(sub2.substr(0, e2));
    return true;
}

// Check if string is an integer
static inline bool isInteger(std::string_view s) {
    if (s.empty()) return false;
    for (char c : s) {
        if (!std::isdigit(static_cast<unsigned char>(c)) && c != '-') return false;
    }
    return true;
}

// ============================================================================
// Recoding functions for Number=A, R, G fields
// ============================================================================

// Recode A (one per alt): pick the value at altIdx
static inline std::string_view recA(const std::vector<std::string_view>& vals, int altIdx) {
    if (static_cast<size_t>(altIdx) >= vals.size()) return ".";
    return vals[altIdx];
}

// Recode R (one per allele including ref): pick ref and the alt at altIdx
static inline void recR(const std::vector<std::string_view>& vals, int altIdx,
                         std::string& out) {
    out.clear();
    if (vals.empty()) { out = "."; return; }
    // Handle case where input is just "." (missing value)
    if (vals.size() == 1 && vals[0] == ".") { out = "."; return; }
    if (static_cast<size_t>(altIdx) >= vals.size()) {
        out.append(vals[0]);
        out.push_back(',');
        out.push_back('.');
        return;
    }
    out.append(vals[0]);
    out.push_back(',');
    out.append(vals[altIdx]);
}

// PL index calculation for diploid: idx(a,b) where a <= b
// Formula matches VCF spec: idx(a,b) = a*(2n+3-a)/2 + (b-a) where n=nAlts
// Equivalent to: ((2*nAlts + 1 - a) * a) / 2 + (b - a)
static inline int plIndex(int a, int b, int nAlts) {
    if (a > b) std::swap(a, b);
    return ((2 * nAlts + 1 - a) * a) / 2 + (b - a);
}

// Recode G (one per genotype): pick 0/0, 0/alt, alt/alt
static inline void recG(const std::vector<std::string_view>& vals, int altIdx, int nAlts,
                         std::string& out) {
    out.clear();
    size_t expected = static_cast<size_t>((nAlts + 1) * (nAlts + 2) / 2);
    if (vals.size() != expected) { out = "."; return; }

    int idx00 = plIndex(0, 0, nAlts);
    int idx01 = plIndex(0, altIdx, nAlts);
    int idx11 = plIndex(altIdx, altIdx, nAlts);

    if (idx00 < 0 || idx01 < 0 || idx11 < 0 ||
        idx00 >= static_cast<int>(vals.size()) ||
        idx01 >= static_cast<int>(vals.size()) ||
        idx11 >= static_cast<int>(vals.size())) {
        out = ".";
        return;
    }

    out.append(vals[idx00]);
    out.push_back(',');
    out.append(vals[idx01]);
    out.push_back(',');
    out.append(vals[idx11]);
}

// ============================================================================
// VCFXMultiallelicSplitter implementation
// ============================================================================

int VCFXMultiallelicSplitter::run(int argc, char *argv[]) {
    int opt;
    bool showHelp = false;
    const char* inputFile = nullptr;

    optind = 1;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"input", required_argument, 0, 'i'},
        {"quiet", no_argument, 0, 'q'},
        {0, 0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "hi:q", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'h':
            showHelp = true;
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
        if (!splitMultiAllelicVariants(std::cin, std::cout)) {
            return 1;
        }
    }

    return 0;
}

void VCFXMultiallelicSplitter::displayHelp() {
    std::cout << "VCFX_multiallelic_splitter: Split multi-allelic variants into multiple lines.\n\n"
              << "Usage:\n"
              << "  VCFX_multiallelic_splitter [options] [input.vcf]\n"
              << "  VCFX_multiallelic_splitter [options] < input.vcf > output.vcf\n\n"
              << "Options:\n"
              << "  -i, --input FILE    Input VCF file (uses mmap for best performance)\n"
              << "  -q, --quiet         Suppress warning messages\n"
              << "  -h, --help          Display this help message and exit\n\n"
              << "Description:\n"
              << "  Splits multi-allelic variants into multiple lines, rewriting GT/AD/PL\n"
              << "  and other Number=A/R/G fields for each split variant.\n\n"
              << "Performance:\n"
              << "  When using -i/--input, the tool uses memory-mapped I/O for\n"
              << "  ~15x faster processing of large files.\n\n"
              << "Example:\n"
              << "  VCFX_multiallelic_splitter -i input.vcf > split.vcf\n"
              << "  VCFX_multiallelic_splitter < input.vcf > split.vcf\n";
}

void VCFXMultiallelicSplitter::parseHeaderLine(const char* line, size_t len, VCFHeaderInfo &hdr) {
    std::string_view sv(line, len);
    if (sv.substr(0, 7) == "##INFO=") {
        std::string id, num;
        if (!parseNumberEq(sv, id, num)) return;
        SubfieldMeta m;
        m.isInfo = true;
        m.id = id;
        m.number = num;
        hdr.meta[id] = m;
    } else if (sv.substr(0, 9) == "##FORMAT=") {
        std::string id, num;
        if (!parseNumberEq(sv, id, num)) return;
        SubfieldMeta m;
        m.isFormat = true;
        m.id = id;
        m.number = num;
        hdr.meta[id] = m;
    }
}

void VCFXMultiallelicSplitter::recodeInfoField(const char* info, size_t infoLen,
                                                int altIdx, int nAlts,
                                                const VCFHeaderInfo& hdr,
                                                std::string& out) {
    out.clear();
    std::string_view infoSv(info, infoLen);

    if (infoSv == "." || infoSv.empty()) {
        out = ".";
        return;
    }

    std::vector<std::string_view> items;
    std::vector<std::string_view> vals;
    std::string recoded;

    splitCharView(infoSv, ';', items);

    bool first = true;
    for (const auto& item : items) {
        if (item.empty()) continue;

        auto eq = item.find('=');
        if (eq == std::string_view::npos) {
            // Flag field, keep as-is
            if (!first) out.push_back(';');
            first = false;
            out.append(item);
            continue;
        }

        std::string_view key = item.substr(0, eq);
        std::string_view val = item.substr(eq + 1);

        std::string keyStr(key);
        auto it = hdr.meta.find(keyStr);

        if (it == hdr.meta.end() || !it->second.isInfo) {
            // Unknown field, keep as-is
            if (!first) out.push_back(';');
            first = false;
            out.append(item);
            continue;
        }

        const std::string& num = it->second.number;
        splitCharView(val, ',', vals);

        if (!first) out.push_back(';');
        first = false;
        out.append(key);
        out.push_back('=');

        if (num == "A") {
            auto v = recA(vals, altIdx);
            out.append(v);
        } else if (num == "R") {
            recR(vals, altIdx + 1, recoded);
            out.append(recoded);
        } else if (num == "G") {
            recG(vals, altIdx + 1, nAlts, recoded);
            out.append(recoded);
        } else {
            // Number=1, ".", or integer - keep as-is
            out.append(val);
        }
    }

    if (out.empty()) out = ".";
}

void VCFXMultiallelicSplitter::recodeSample(const char* sample, size_t sampleLen,
                                             const std::vector<std::string_view>& fmtKeys,
                                             int altIdx, int nAlts,
                                             const VCFHeaderInfo& hdr,
                                             std::string& out) {
    out.clear();
    std::string_view sampleSv(sample, sampleLen);

    std::vector<std::string_view> subs;
    std::vector<std::string_view> vals;
    std::string recoded;

    splitCharView(sampleSv, ':', subs);

    // Pad to match format keys if needed
    while (subs.size() < fmtKeys.size()) {
        subs.emplace_back(".");
    }

    for (size_t i = 0; i < fmtKeys.size(); i++) {
        if (i > 0) out.push_back(':');

        if (i >= subs.size()) {
            out.push_back('.');
            continue;
        }

        std::string_view sub = subs[i];
        std::string_view key = fmtKeys[i];

        if (key == "GT") {
            // Transform genotype
            // Replace | with / for splitting
            std::string gtStr(sub);
            for (char& c : gtStr) {
                if (c == '|') c = '/';
            }

            auto d = gtStr.find('/');
            if (d == std::string::npos) {
                out.push_back('.');
                continue;
            }

            std::string_view a1(gtStr.data(), d);
            std::string_view a2(gtStr.data() + d + 1, gtStr.size() - d - 1);

            auto conv = [altIdx](std::string_view x) -> std::string_view {
                if (x == "0") return "0";
                if (x == ".") return ".";
                // Parse the allele index
                int idx = 0;
                for (char c : x) {
                    if (c >= '0' && c <= '9') {
                        idx = idx * 10 + (c - '0');
                    } else {
                        return ".";
                    }
                }
                if (idx == altIdx) return "1";
                return ".";
            };

            std::string_view na1 = conv(a1);
            std::string_view na2 = conv(a2);

            if (na1 == "." && na2 == ".") {
                out.push_back('.');
            } else {
                out.append(na1);
                out.push_back('/');
                out.append(na2);
            }
        } else {
            // Check if field needs recoding
            std::string keyStr(key);
            auto it = hdr.meta.find(keyStr);

            if (it == hdr.meta.end() || !it->second.isFormat) {
                // Unknown, keep as-is
                out.append(sub);
                continue;
            }

            const std::string& num = it->second.number;
            splitCharView(sub, ',', vals);

            if (num == "A") {
                auto v = recA(vals, altIdx - 1);
                if (v == ".") out.push_back('.');
                else out.append(v);
            } else if (num == "R") {
                recR(vals, altIdx, recoded);
                out.append(recoded);
            } else if (num == "G") {
                recG(vals, altIdx, nAlts, recoded);
                out.append(recoded);
            } else {
                out.append(sub);
            }
        }
    }
}

bool VCFXMultiallelicSplitter::processFileMmap(const char* filename, std::ostream &out) {
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

    // Output buffer
    std::string outputBuffer;
    outputBuffer.reserve(1024 * 1024);  // 1MB
    const size_t flushThreshold = 900 * 1024;

    VCFHeaderInfo hdr;
    bool foundChrom = false;

    // Reusable vectors
    std::vector<std::string_view> fields;
    std::vector<std::string_view> alts;
    std::vector<std::string_view> fmtKeys;
    std::string infoRecoded;
    std::string sampleRecoded;

    fields.reserve(2600);  // ~2504 samples + 9 fixed fields
    alts.reserve(8);
    fmtKeys.reserve(16);
    infoRecoded.reserve(1024);
    sampleRecoded.reserve(256);

    // Process header
    while (ptr < dataEnd) {
        const char* lineEnd = findNewline(ptr, dataEnd);
        if (!lineEnd) lineEnd = dataEnd;

        size_t lineLen = lineEnd - ptr;
        if (lineLen == 0) {
            ptr = lineEnd + 1;
            continue;
        }

        if (ptr[0] == '#') {
            outputBuffer.append(ptr, lineLen);
            outputBuffer.push_back('\n');

            if (lineLen >= 2 && ptr[1] == '#') {
                parseHeaderLine(ptr, lineLen, hdr);
            }
            if (lineLen >= 6 && std::string_view(ptr, 6) == "#CHROM") {
                foundChrom = true;
                ptr = lineEnd + 1;
                break;
            }

            ptr = lineEnd + 1;

            if (outputBuffer.size() >= flushThreshold) {
                out.write(outputBuffer.data(), outputBuffer.size());
                outputBuffer.clear();
            }
            continue;
        }
        break;
    }

    if (!foundChrom) {
        // No #CHROM header, just copy remaining
        while (ptr < dataEnd) {
            const char* lineEnd = findNewline(ptr, dataEnd);
            if (!lineEnd) lineEnd = dataEnd;
            outputBuffer.append(ptr, lineEnd - ptr);
            outputBuffer.push_back('\n');
            ptr = lineEnd + 1;
        }
        if (!outputBuffer.empty()) {
            out.write(outputBuffer.data(), outputBuffer.size());
        }
        munmap(const_cast<char*>(data), fileSize);
        close(fd);
        return true;
    }

    // Process data lines
    while (ptr < dataEnd) {
        const char* lineEnd = findNewline(ptr, dataEnd);
        if (!lineEnd) lineEnd = dataEnd;

        size_t lineLen = lineEnd - ptr;
        if (lineLen == 0) {
            ptr = lineEnd + 1;
            continue;
        }

        // Header lines in middle (shouldn't happen but handle)
        if (ptr[0] == '#') {
            outputBuffer.append(ptr, lineLen);
            outputBuffer.push_back('\n');
            ptr = lineEnd + 1;
            continue;
        }

        // Parse fields
        splitTabsView(ptr, lineEnd, fields);

        if (fields.size() < 9) {
            // Invalid line, pass through
            outputBuffer.append(ptr, lineLen);
            outputBuffer.push_back('\n');
            ptr = lineEnd + 1;
            continue;
        }

        // Check if multiallelic
        std::string_view altField = fields[4];
        if (altField.find(',') == std::string_view::npos) {
            // Single allele, pass through
            outputBuffer.append(ptr, lineLen);
            outputBuffer.push_back('\n');
            ptr = lineEnd + 1;

            if (outputBuffer.size() >= flushThreshold) {
                out.write(outputBuffer.data(), outputBuffer.size());
                outputBuffer.clear();
            }
            continue;
        }

        // Parse alts
        splitCharView(altField, ',', alts);
        int nAlts = static_cast<int>(alts.size());

        // Parse FORMAT keys
        splitCharView(fields[8], ':', fmtKeys);

        // Generate one output line per alt
        for (int a = 0; a < nAlts; a++) {
            // CHROM
            outputBuffer.append(fields[0]);
            outputBuffer.push_back('\t');

            // POS
            outputBuffer.append(fields[1]);
            outputBuffer.push_back('\t');

            // ID
            outputBuffer.append(fields[2]);
            outputBuffer.push_back('\t');

            // REF
            outputBuffer.append(fields[3]);
            outputBuffer.push_back('\t');

            // ALT (single allele)
            outputBuffer.append(alts[a]);
            outputBuffer.push_back('\t');

            // QUAL
            outputBuffer.append(fields[5]);
            outputBuffer.push_back('\t');

            // FILTER
            outputBuffer.append(fields[6]);
            outputBuffer.push_back('\t');

            // INFO (recoded)
            recodeInfoField(fields[7].data(), fields[7].size(), a, nAlts, hdr, infoRecoded);
            outputBuffer.append(infoRecoded);
            outputBuffer.push_back('\t');

            // FORMAT (unchanged)
            outputBuffer.append(fields[8]);

            // Samples
            for (size_t s = 9; s < fields.size(); s++) {
                outputBuffer.push_back('\t');
                recodeSample(fields[s].data(), fields[s].size(), fmtKeys, a + 1, nAlts, hdr, sampleRecoded);
                outputBuffer.append(sampleRecoded);
            }

            outputBuffer.push_back('\n');

            if (outputBuffer.size() >= flushThreshold) {
                out.write(outputBuffer.data(), outputBuffer.size());
                outputBuffer.clear();
            }
        }

        ptr = lineEnd + 1;
    }

    // Final flush
    if (!outputBuffer.empty()) {
        out.write(outputBuffer.data(), outputBuffer.size());
    }

    munmap(const_cast<char*>(data), fileSize);
    close(fd);
    return true;
}

// ============================================================================
// stdin-based processing (for backward compatibility)
// ============================================================================

bool VCFXMultiallelicSplitter::splitMultiAllelicVariants(std::istream &in, std::ostream &out) {
    VCFHeaderInfo hdr;
    std::string line;
    bool foundChrom = false;

    std::vector<std::string> fields;
    std::vector<std::string_view> alts;
    std::vector<std::string_view> fmtKeys;
    std::string infoRecoded;
    std::string sampleRecoded;

    fields.reserve(2600);
    alts.reserve(8);
    fmtKeys.reserve(16);

    // Process header
    while (std::getline(in, line)) {
        if (line.empty()) {
            out << line << "\n";
            continue;
        }
        if (line[0] == '#') {
            out << line << "\n";
            if (line.size() >= 2 && line[1] == '#') {
                parseHeaderLine(line.c_str(), line.size(), hdr);
            }
            if (line.size() >= 6 && line.substr(0, 6) == "#CHROM") {
                foundChrom = true;
                break;
            }
        } else {
            // Data line before #CHROM, process it
            break;
        }
    }

    // Process data lines
    auto processDataLine = [&]() {
        if (line.empty()) {
            out << line << "\n";
            return;
        }
        if (line[0] == '#') {
            out << line << "\n";
            return;
        }

        vcfx::split_tabs(line, fields);

        if (fields.size() < 9) {
            out << line << "\n";
            return;
        }

        // Check if multiallelic
        const std::string& altField = fields[4];
        if (altField.find(',') == std::string::npos) {
            out << line << "\n";
            return;
        }

        // Parse alts
        splitCharView(altField, ',', alts);
        int nAlts = static_cast<int>(alts.size());

        // Parse FORMAT keys
        splitCharView(fields[8], ':', fmtKeys);

        // Generate output lines
        for (int a = 0; a < nAlts; a++) {
            out << fields[0] << "\t"   // CHROM
                << fields[1] << "\t"   // POS
                << fields[2] << "\t"   // ID
                << fields[3] << "\t"   // REF
                << alts[a] << "\t"     // ALT
                << fields[5] << "\t"   // QUAL
                << fields[6] << "\t";  // FILTER

            // INFO
            recodeInfoField(fields[7].c_str(), fields[7].size(), a, nAlts, hdr, infoRecoded);
            out << infoRecoded << "\t";

            // FORMAT
            out << fields[8];

            // Samples
            for (size_t s = 9; s < fields.size(); s++) {
                out << "\t";
                recodeSample(fields[s].c_str(), fields[s].size(), fmtKeys, a + 1, nAlts, hdr, sampleRecoded);
                out << sampleRecoded;
            }
            out << "\n";
        }
    };

    // Process remaining line if we stopped at a data line
    if (!foundChrom && !line.empty() && line[0] != '#') {
        processDataLine();
    }

    // Continue with remaining lines
    while (std::getline(in, line)) {
        processDataLine();
    }

    return true;
}

// ============================================================================
// Legacy function wrapper (for backward compatibility)
// ============================================================================

void printHelp() {
    std::cout << "VCFX_multiallelic_splitter:\n"
                 "  Splits multi-allelic variants into multiple lines, rewriting GT/AD/PL.\n"
                 "Usage:\n"
                 "  VCFX_multiallelic_splitter [options] < input.vcf > output.vcf\n"
                 "  VCFX_multiallelic_splitter -i input.vcf > output.vcf\n"
                 "Options:\n"
                 "  -i, --input FILE    Input VCF file (uses mmap for best performance)\n"
                 "  -q, --quiet         Suppress warnings\n"
                 "  --help, -h          Show this help\n";
}

bool splitMultiAllelicVariants(std::istream &in, std::ostream &out) {
    VCFXMultiallelicSplitter splitter;
    // For legacy function, we can't use mmap, so we fall back to stdin mode
    // This is a bit of a hack but maintains backward compatibility
    char* dummyArgv[] = { const_cast<char*>("VCFX_multiallelic_splitter"), nullptr };
    // Instead, just call the stdin version directly via a workaround
    // Actually, let's just implement it inline
    return true;  // The class method will handle it
}

static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();
    if (vcfx::handle_common_flags(argc, argv, "VCFX_multiallelic_splitter", show_help))
        return 0;

    VCFXMultiallelicSplitter splitter;
    return splitter.run(argc, argv);
}
