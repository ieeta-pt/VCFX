#include "VCFX_sorter.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <unistd.h>

// MergeEntry comparator for min-heap
// Uses natural or lexicographic comparison based on static setting
static bool g_naturalOrder = false;

bool MergeEntry::operator>(const MergeEntry& other) const {
    if (g_naturalOrder) {
        std::string apfx, asuf, bpfx, bsuf;
        long anum, bnum;
        VCFXSorter::parseChromNat(chrom, apfx, anum, asuf);
        VCFXSorter::parseChromNat(other.chrom, bpfx, bnum, bsuf);

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

void VCFXSorter::displayHelp() {
    std::cout
        << "VCFX_sorter: Sort a VCF by chromosome and position.\n\n"
           "Usage:\n"
           "  VCFX_sorter [options] < input.vcf > output.vcf\n\n"
           "Options:\n"
           "  -h, --help              Show help.\n"
           "  -n, --natural-chr       Use natural chromosome order (chr1 < chr2 < chr10).\n"
           "  -m, --max-memory <MB>   Max memory for in-memory sorting (default: 100MB).\n"
           "                          Files larger than this use external merge sort.\n"
           "  -t, --temp-dir <DIR>    Directory for temporary files (default: /tmp).\n\n"
           "Description:\n"
           "  Sorts VCF by (CHROM, POS). For small files, sorts in memory.\n"
           "  For large files (>max-memory), uses external merge sort with\n"
           "  temporary files, enabling sorting of files larger than RAM.\n\n"
           "Examples:\n"
           "  1) Default (in-memory for small files):\n"
           "     VCFX_sorter < unsorted.vcf > sorted.vcf\n"
           "  2) Large file with natural chromosome order:\n"
           "     VCFX_sorter -n -m 500 < huge.vcf > sorted.vcf\n"
           "  3) Use custom temp directory:\n"
           "     VCFX_sorter -t /data/tmp < huge.vcf > sorted.vcf\n";
}

// Parse chromosome in natural manner: "chr10" => ("chr", 10, "")
bool VCFXSorter::parseChromNat(const std::string &chrom, std::string &prefix, long &num, std::string &suffix) {
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
bool VCFXSorter::parseChromPos(const std::string& line, std::string& chrom, int& pos) {
    // Find first tab (end of CHROM)
    size_t tab1 = line.find('\t');
    if (tab1 == std::string::npos) return false;

    chrom = line.substr(0, tab1);

    // Find second tab (end of POS)
    size_t tab2 = line.find('\t', tab1 + 1);
    if (tab2 == std::string::npos) return false;

    try {
        pos = std::stoi(line.substr(tab1 + 1, tab2 - tab1 - 1));
    } catch (...) {
        return false;
    }
    return true;
}

// Lexicographic comparison
bool VCFXSorter::lexCompare(const SortKey &a, const SortKey &b) {
    if (a.chrom != b.chrom) return a.chrom < b.chrom;
    return a.pos < b.pos;
}

// Natural chromosome comparison
bool VCFXSorter::naturalCompare(const SortKey &a, const SortKey &b) {
    std::string apfx, asuf, bpfx, bsuf;
    long anum, bnum;

    if (!parseChromNat(a.chrom, apfx, anum, asuf) ||
        !parseChromNat(b.chrom, bpfx, bnum, bsuf)) {
        return lexCompare(a, b);
    }

    if (apfx != bpfx) return apfx < bpfx;

    if (anum >= 0 && bnum >= 0) {
        if (anum != bnum) return anum < bnum;
        if (asuf != bsuf) return asuf < bsuf;
        return a.pos < b.pos;
    } else if (anum >= 0 && bnum < 0) {
        return true;  // numeric < no numeric
    } else if (anum < 0 && bnum >= 0) {
        return false;
    } else {
        if (a.chrom != b.chrom) return a.chrom < b.chrom;
        return a.pos < b.pos;
    }
}

// In-memory sort for small files
void VCFXSorter::sortInMemory(std::istream &in, std::ostream &out) {
    std::vector<SortKey> records;
    std::string line;

    while (std::getline(in, line)) {
        if (line.empty()) {
            headerLines.push_back(line);
            continue;
        }
        if (line[0] == '#') {
            headerLines.push_back(line);
            continue;
        }

        SortKey key;
        if (!parseChromPos(line, key.chrom, key.pos)) {
            std::cerr << "Warning: skipping malformed line.\n";
            continue;
        }
        key.line = std::move(line);
        records.push_back(std::move(key));
    }

    // Sort records
    std::sort(records.begin(), records.end(), compareFunc);

    // Output header
    for (const auto& h : headerLines) {
        out << h << "\n";
    }

    // Output sorted records
    for (const auto& rec : records) {
        out << rec.line << "\n";
    }
}

// Write a sorted chunk to a temp file
std::string VCFXSorter::writeChunk(std::vector<SortKey>& chunk, size_t chunk_num) {
    // Sort the chunk
    std::sort(chunk.begin(), chunk.end(), compareFunc);

    // Create temp file
    std::string path = tempDir + "/vcfx_sort_" + std::to_string(getpid()) +
                       "_chunk_" + std::to_string(chunk_num) + ".tmp";

    std::ofstream ofs(path, std::ios::binary);
    if (!ofs) {
        std::cerr << "Error: cannot create temp file " << path << "\n";
        return "";
    }

    for (const auto& key : chunk) {
        ofs << key.line << "\n";
    }

    ofs.close();
    return path;
}

// K-way merge of sorted temp files using min-heap
void VCFXSorter::mergeChunks(const std::vector<std::string>& chunk_files, std::ostream &out) {
    // Open all chunk files
    std::vector<std::unique_ptr<std::ifstream>> files;
    files.reserve(chunk_files.size());

    for (const auto& path : chunk_files) {
        auto ifs = std::make_unique<std::ifstream>(path, std::ios::binary);
        if (!ifs->is_open()) {
            std::cerr << "Error: cannot open temp file " << path << "\n";
            continue;
        }
        files.push_back(std::move(ifs));
    }

    // Min-heap using greater comparison (so smallest is at top)
    std::priority_queue<MergeEntry, std::vector<MergeEntry>, std::greater<MergeEntry>> heap;

    // Initialize heap with first line from each file
    for (size_t i = 0; i < files.size(); ++i) {
        std::string line;
        if (std::getline(*files[i], line)) {
            MergeEntry entry;
            entry.file_index = i;
            entry.line = std::move(line);
            if (parseChromPos(entry.line, entry.chrom, entry.pos)) {
                heap.push(std::move(entry));
            }
        }
    }

    // Merge
    while (!heap.empty()) {
        MergeEntry top = heap.top();
        heap.pop();

        out << top.line << "\n";

        // Read next line from same file
        std::string line;
        if (std::getline(*files[top.file_index], line)) {
            MergeEntry entry;
            entry.file_index = top.file_index;
            entry.line = std::move(line);
            if (parseChromPos(entry.line, entry.chrom, entry.pos)) {
                heap.push(std::move(entry));
            }
        }
    }

    // Close and delete temp files
    files.clear();
    for (const auto& path : chunk_files) {
        std::remove(path.c_str());
    }
}

// External merge sort for large files
void VCFXSorter::sortExternal(std::istream &in, std::ostream &out) {
    std::vector<std::string> chunk_files;
    std::vector<SortKey> current_chunk;
    size_t current_bytes = 0;
    size_t chunk_num = 0;
    std::string line;

    // First pass: read headers and split data into sorted chunks
    while (std::getline(in, line)) {
        if (line.empty()) {
            headerLines.push_back(line);
            continue;
        }
        if (line[0] == '#') {
            headerLines.push_back(line);
            continue;
        }

        SortKey key;
        if (!parseChromPos(line, key.chrom, key.pos)) {
            std::cerr << "Warning: skipping malformed line.\n";
            continue;
        }

        size_t line_bytes = line.size() + key.chrom.size() + sizeof(key);
        key.line = std::move(line);

        current_chunk.push_back(std::move(key));
        current_bytes += line_bytes;

        // If chunk is full, write it out
        if (current_bytes >= chunkSizeBytes) {
            std::string path = writeChunk(current_chunk, chunk_num++);
            if (!path.empty()) {
                chunk_files.push_back(path);
            }
            current_chunk.clear();
            current_bytes = 0;
        }
    }

    // Write remaining chunk
    if (!current_chunk.empty()) {
        std::string path = writeChunk(current_chunk, chunk_num++);
        if (!path.empty()) {
            chunk_files.push_back(path);
        }
        current_chunk.clear();
    }

    // Output headers
    for (const auto& h : headerLines) {
        out << h << "\n";
    }

    // If only one chunk, it's already sorted - just copy it
    if (chunk_files.size() == 1) {
        std::ifstream ifs(chunk_files[0], std::ios::binary);
        std::string line;
        while (std::getline(ifs, line)) {
            out << line << "\n";
        }
        std::remove(chunk_files[0].c_str());
        return;
    }

    // K-way merge
    if (!chunk_files.empty()) {
        mergeChunks(chunk_files, out);
    }
}

int VCFXSorter::run(int argc, char *argv[]) {
    bool showHelp = false;
    bool useExternal = false;  // Auto-detect by default

    static struct option long_opts[] = {
        {"help", no_argument, 0, 'h'},
        {"natural-chr", no_argument, 0, 'n'},
        {"max-memory", required_argument, 0, 'm'},
        {"temp-dir", required_argument, 0, 't'},
        {"external", no_argument, 0, 'e'},  // Force external sort
        {0, 0, 0, 0}
    };

    while (true) {
        int c = ::getopt_long(argc, argv, "hnm:t:e", long_opts, nullptr);
        if (c == -1) break;

        switch (c) {
        case 'h':
            showHelp = true;
            break;
        case 'n':
            naturalChromOrder = true;
            break;
        case 'm':
            chunkSizeBytes = std::stoul(optarg) * 1024 * 1024;  // Convert MB to bytes
            break;
        case 't':
            tempDir = optarg;
            break;
        case 'e':
            useExternal = true;
            break;
        default:
            showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    // Set global comparison mode for MergeEntry
    g_naturalOrder = naturalChromOrder;

    // Set comparison function
    if (naturalChromOrder) {
        compareFunc = naturalCompare;
    } else {
        compareFunc = lexCompare;
    }

    // Use external sort for streaming from stdin
    // (can't seek back to measure size, so we always use external sort
    //  which handles both small and large files efficiently)
    if (useExternal) {
        sortExternal(std::cin, std::cout);
    } else {
        // For stdin, we use external sort since we can't know the size upfront
        // The external sort will create only one chunk for small files anyway
        sortExternal(std::cin, std::cout);
    }

    return 0;
}

static void show_help() {
    VCFXSorter obj;
    char arg0[] = "VCFX_sorter";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio

    if (vcfx::handle_common_flags(argc, argv, "VCFX_sorter", show_help))
        return 0;

    VCFXSorter app;
    return app.run(argc, argv);
}
