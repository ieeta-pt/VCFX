#include "VCFX_merger.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <queue>
#include <sstream>

// Global flag for MergeFileEntry comparison mode
static bool g_naturalOrder = false;

// MergeFileEntry comparator for min-heap
bool MergeFileEntry::operator>(const MergeFileEntry& other) const {
    if (g_naturalOrder) {
        std::string apfx, asuf, bpfx, bsuf;
        long anum, bnum;
        VCFXMerger::parseChromNat(chrom, apfx, anum, asuf);
        VCFXMerger::parseChromNat(other.chrom, bpfx, bnum, bsuf);

        if (apfx != bpfx) return apfx > bpfx;
        if (anum >= 0 && bnum >= 0) {
            if (anum != bnum) return anum > bnum;
            if (asuf != bsuf) return asuf > bsuf;
        } else if (anum >= 0 && bnum < 0) {
            return false;  // numeric < no numeric
        } else if (anum < 0 && bnum >= 0) {
            return true;
        } else {
            if (chrom != other.chrom) return chrom > other.chrom;
        }
        return pos > other.pos;
    } else {
        if (chrom != other.chrom) return chrom > other.chrom;
        return pos > other.pos;
    }
}

// Parse chromosome in natural manner: "chr10" => ("chr", 10, "")
bool VCFXMerger::parseChromNat(const std::string &chrom, std::string &prefix, long &num, std::string &suffix) {
    std::string c = chrom;
    std::string up;
    up.reserve(c.size());
    for (char ch : c)
        up.push_back(std::toupper(ch));

    if (up.rfind("CHR", 0) == 0) {
        prefix = c.substr(0, 3);
        c = c.substr(3);
    } else {
        prefix = "";
    }

    size_t idx = 0;
    while (idx < c.size() && std::isdigit(static_cast<unsigned char>(c[idx])))
        idx++;

    if (idx == 0) {
        num = -1;
        suffix = c;
        return true;
    }

    try {
        num = std::stol(c.substr(0, idx));
    } catch (...) {
        return false;
    }
    suffix = c.substr(idx);
    return true;
}

// Parse CHROM and POS from a VCF data line
bool VCFXMerger::parseChromPos(const std::string& line, std::string& chrom, long& pos) {
    size_t tab1 = line.find('\t');
    if (tab1 == std::string::npos) return false;

    chrom = line.substr(0, tab1);

    size_t tab2 = line.find('\t', tab1 + 1);
    if (tab2 == std::string::npos) return false;

    try {
        pos = std::stol(line.substr(tab1 + 1, tab2 - tab1 - 1));
    } catch (...) {
        return false;
    }
    return true;
}

int VCFXMerger::run(int argc, char *argv[]) {
    int opt;
    bool showHelp = false;
    std::vector<std::string> inputFiles;

    static struct option long_options[] = {
        {"merge", required_argument, 0, 'm'},
        {"assume-sorted", no_argument, 0, 's'},
        {"natural-chr", no_argument, 0, 'n'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "m:snh", long_options, nullptr)) != -1) {
        switch (opt) {
        case 'm': {
            // Split comma-separated file names
            std::string files = optarg;
            size_t pos = 0;
            while ((pos = files.find(',')) != std::string::npos) {
                inputFiles.emplace_back(files.substr(0, pos));
                files.erase(0, pos + 1);
            }
            inputFiles.emplace_back(files);
            break;
        }
        case 's':
            assumeSorted = true;
            break;
        case 'n':
            naturalChromOrder = true;
            break;
        case 'h':
        default:
            showHelp = true;
        }
    }

    if (showHelp || inputFiles.empty()) {
        displayHelp();
        return 0;
    }

    // Set global comparison mode
    g_naturalOrder = naturalChromOrder;

    // Choose merge strategy
    if (assumeSorted) {
        mergeVCFStreaming(inputFiles, std::cout);
    } else {
        mergeVCFInMemory(inputFiles, std::cout);
    }

    return 0;
}

void VCFXMerger::displayHelp() {
    std::cout << "VCFX_merger: Merge multiple VCF files by variant position.\n\n"
              << "Usage:\n"
              << "  VCFX_merger --merge file1.vcf,file2.vcf,... [options]\n\n"
              << "Options:\n"
              << "  -m, --merge          Comma-separated list of VCF files to merge\n"
              << "  -s, --assume-sorted  Assume input files are already sorted (enables streaming\n"
              << "                       merge with O(num_files) memory for large files)\n"
              << "  -n, --natural-chr    Use natural chromosome order (chr1 < chr2 < chr10)\n"
              << "  -h, --help           Display this help message and exit\n\n"
              << "Description:\n"
              << "  By default, loads all variants into memory and sorts them. This works for\n"
              << "  small to medium files but may run out of memory for very large files.\n\n"
              << "  With --assume-sorted, uses streaming K-way merge that only keeps one line\n"
              << "  per input file in memory. This enables merging files larger than RAM.\n"
              << "  Input files MUST be sorted by (CHROM, POS) for correct results.\n\n"
              << "Examples:\n"
              << "  # Default mode (loads all into memory, sorts)\n"
              << "  VCFX_merger --merge sample1.vcf,sample2.vcf > merged.vcf\n\n"
              << "  # Streaming mode for large pre-sorted files\n"
              << "  VCFX_merger --merge sorted1.vcf,sorted2.vcf --assume-sorted > merged.vcf\n\n"
              << "  # With natural chromosome ordering\n"
              << "  VCFX_merger --merge f1.vcf,f2.vcf --assume-sorted -n > merged.vcf\n";
}

// Original in-memory merge (for unsorted inputs or small files)
void VCFXMerger::mergeVCFInMemory(const std::vector<std::string> &inputFiles, std::ostream &out) {
    struct Variant {
        std::string chrom;
        long pos = 0;
        std::string line;
    };

    std::vector<Variant> variants;
    std::vector<std::string> headers;
    bool headersCaptured = false;

    for (const auto &file : inputFiles) {
        std::ifstream stream(file);
        if (!stream.is_open()) {
            std::cerr << "Failed to open file: " << file << "\n";
            continue;
        }

        std::string line;
        while (std::getline(stream, line)) {
            if (line.empty())
                continue;
            if (line[0] == '#') {
                if (!headersCaptured)
                    headers.push_back(line);
                continue;
            }

            Variant v;
            if (!parseChromPos(line, v.chrom, v.pos)) {
                std::cerr << "Warning: skipping malformed line in " << file << "\n";
                continue;
            }
            v.line = std::move(line);
            variants.push_back(std::move(v));
        }

        if (!headersCaptured && !headers.empty())
            headersCaptured = true;
    }

    for (const auto &h : headers) {
        out << h << '\n';
    }

    // Sort with natural or lexicographic ordering
    if (naturalChromOrder) {
        std::sort(variants.begin(), variants.end(), [](const Variant &a, const Variant &b) {
            std::string apfx, asuf, bpfx, bsuf;
            long anum, bnum;
            VCFXMerger::parseChromNat(a.chrom, apfx, anum, asuf);
            VCFXMerger::parseChromNat(b.chrom, bpfx, bnum, bsuf);

            if (apfx != bpfx) return apfx < bpfx;
            if (anum >= 0 && bnum >= 0) {
                if (anum != bnum) return anum < bnum;
                if (asuf != bsuf) return asuf < bsuf;
                return a.pos < b.pos;
            } else if (anum >= 0 && bnum < 0) {
                return true;
            } else if (anum < 0 && bnum >= 0) {
                return false;
            } else {
                if (a.chrom != b.chrom) return a.chrom < b.chrom;
                return a.pos < b.pos;
            }
        });
    } else {
        std::sort(variants.begin(), variants.end(), [](const Variant &a, const Variant &b) {
            if (a.chrom == b.chrom)
                return a.pos < b.pos;
            return a.chrom < b.chrom;
        });
    }

    for (const auto &v : variants) {
        out << v.line << '\n';
    }
}

// Streaming K-way merge for large pre-sorted files - O(num_files) memory
void VCFXMerger::mergeVCFStreaming(const std::vector<std::string> &inputFiles, std::ostream &out) {
    // Open all input files
    std::vector<std::unique_ptr<std::ifstream>> files;
    files.reserve(inputFiles.size());

    std::vector<std::string> allHeaders;
    bool headersCaptured = false;

    for (const auto &path : inputFiles) {
        auto ifs = std::make_unique<std::ifstream>(path);
        if (!ifs->is_open()) {
            std::cerr << "Failed to open file: " << path << "\n";
            continue;
        }

        // Read headers from this file
        std::string line;
        std::streampos dataStart = 0;
        while (std::getline(*ifs, line)) {
            if (line.empty()) continue;
            if (line[0] == '#') {
                if (!headersCaptured) {
                    allHeaders.push_back(line);
                }
                dataStart = ifs->tellg();
            } else {
                // Found first data line, seek back
                ifs->seekg(dataStart);
                break;
            }
        }

        if (!headersCaptured && !allHeaders.empty()) {
            headersCaptured = true;
        }

        files.push_back(std::move(ifs));
    }

    // Output headers
    for (const auto &h : allHeaders) {
        out << h << '\n';
    }

    if (files.empty()) {
        return;
    }

    // Min-heap using greater comparison (smallest at top)
    std::priority_queue<MergeFileEntry, std::vector<MergeFileEntry>, std::greater<MergeFileEntry>> heap;

    // Initialize heap with first data line from each file
    for (size_t i = 0; i < files.size(); ++i) {
        std::string line;
        while (std::getline(*files[i], line)) {
            if (line.empty() || line[0] == '#') continue;

            MergeFileEntry entry;
            entry.file_index = i;
            if (parseChromPos(line, entry.chrom, entry.pos)) {
                entry.line = std::move(line);
                heap.push(std::move(entry));
                break;
            }
        }
    }

    // K-way merge
    while (!heap.empty()) {
        MergeFileEntry top = heap.top();
        heap.pop();

        out << top.line << '\n';

        // Read next line from the same file
        std::string line;
        while (std::getline(*files[top.file_index], line)) {
            if (line.empty() || line[0] == '#') continue;

            MergeFileEntry entry;
            entry.file_index = top.file_index;
            if (parseChromPos(line, entry.chrom, entry.pos)) {
                entry.line = std::move(line);
                heap.push(std::move(entry));
                break;
            }
        }
    }
}

static void show_help() {
    VCFXMerger obj;
    char arg0[] = "VCFX_merger";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_merger", show_help))
        return 0;
    VCFXMerger merger;
    return merger.run(argc, argv);
}
