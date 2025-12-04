#include "VCFX_variant_counter.h"
#include "vcfx_core.h" 
#include "vcfx_io.h"
#include <cstring>
#include <fcntl.h>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <zlib.h>

void VCFXVariantCounter::displayHelp() {
    std::cout << "VCFX_variant_counter: Counts the total number of valid variants in a VCF.\n\n"
                 "Usage:\n"
                 "  VCFX_variant_counter [options] [input.vcf]\n"
                 "  VCFX_variant_counter [options] < input.vcf\n\n"
                 "Options:\n"
                 "  -h, --help        Show this help.\n"
                 "  -s, --strict      Fail on any data line with <8 columns.\n\n"
                 "Description:\n"
                 "  Reads a VCF from file argument or stdin. For each data line,\n"
                 "  we check if it has >=8 columns; if it does, we count it; if fewer columns:\n"
                 "   * if --strict => we exit with error,\n"
                 "   * otherwise => we skip with a warning.\n"
                 "  When a file is provided directly, uses memory-mapped I/O for faster processing.\n"
                 "  Finally, we print 'Total Variants: X'.\n\n"
                 "Example:\n"
                 "  VCFX_variant_counter input.vcf          # Fast memory-mapped mode\n"
                 "  VCFX_variant_counter < input.vcf        # Stdin mode\n"
                 "  VCFX_variant_counter --strict input.vcf\n";
}

int VCFXVariantCounter::run(int argc, char *argv[]) {
    bool showHelp = false;
    static struct option long_opts[] = {{"help", no_argument, 0, 'h'}, {"strict", no_argument, 0, 's'}, {0, 0, 0, 0}};

    // Reset getopt
    optind = 1;

    // Check for options/flags
    while (true) {
        int c = ::getopt_long(argc, argv, "hs", long_opts, nullptr);
        if (c == -1)
            break;
        switch (c) {
        case 'h':
            showHelp = true;
            break;
        case 's':
            strictMode = true;
            break;
        default:
            showHelp = true;
        }
    }

    if (showHelp) {
        displayHelp();
        return 0;
    }

    int total = -1;

    // Check if a file argument was provided
    if (optind < argc) {
        const char *filename = argv[optind];
        total = countVariantsMmap(filename);
    } else {
        // Read from stdin
        auto peek1 = std::cin.peek();
        bool isEmpty = (peek1 == EOF);
        bool isGzip = false;
        if (!isEmpty) {
            int c1 = std::cin.get();
            int c2 = std::cin.get();
            if (c2 != EOF) {
                isGzip = (static_cast<unsigned char>(c1) == 0x1f && static_cast<unsigned char>(c2) == 0x8b);
                std::cin.putback(static_cast<char>(c2));
            }
            std::cin.putback(static_cast<char>(c1));
        }

        if (isEmpty) {
            total = 0;
        } else if (isGzip) {
            total = countVariantsGzip(std::cin);
        } else {
            total = countVariants(std::cin);
        }
    }

    if (total < 0) {
        return 1;
    }
    std::cout << "Total Variants: " << total << "\n";
    return 0;
}

bool VCFXVariantCounter::processLine(const std::string &line, int lineNumber, int &count) {
    if (line.empty())
        return true;
    if (line[0] == '#')
        return true;

    // Fast tab counting - no allocations needed
    int tabCount = 0;
    for (char c : line) {
        if (c == '\t') {
            tabCount++;
            if (tabCount >= 7) {
                // 8 columns = 7 tabs minimum
                count++;
                return true;
            }
        }
    }

    // Less than 8 columns
    if (strictMode) {
        std::cerr << "Error: line " << lineNumber << " has <8 columns.\n";
        return false;
    } else {
        std::cerr << "Warning: skipping line " << lineNumber << " with <8 columns.\n";
        return true;
    }
}

int VCFXVariantCounter::countVariants(std::istream &in) {
    // Use larger buffer for faster I/O
    constexpr size_t BUFFER_SIZE = 1 << 20; // 1MB buffer
    std::vector<char> buffer(BUFFER_SIZE);
    in.rdbuf()->pubsetbuf(buffer.data(), BUFFER_SIZE);

    int count = 0;
    int lineNumber = 0;
    std::string line;
    line.reserve(4096); // Pre-allocate for typical VCF line length

    while (std::getline(in, line)) {
        lineNumber++;
        if (!processLine(line, lineNumber, count))
            return -1;
    }
    return count;
}

int VCFXVariantCounter::countVariantsGzip(std::istream &in) {
    constexpr int CHUNK = 65536;  // Larger chunk for better performance
    std::vector<char> inBuf(CHUNK);
    std::vector<char> outBuf(CHUNK);
    z_stream strm;
    std::memset(&strm, 0, sizeof(strm));
    if (inflateInit2(&strm, 15 + 32) != Z_OK) {
        std::cerr << "Error: inflateInit2 failed.\n";
        return -1;
    }
    int count = 0;
    int lineNumber = 0;
    std::string buffer;
    buffer.reserve(65536);
    size_t bufferOffset = 0;  // Track consumed data offset
    int ret = Z_OK;
    do {
        in.read(inBuf.data(), CHUNK);
        strm.avail_in = static_cast<uInt>(in.gcount());
        if (strm.avail_in == 0 && in.eof())
            break;
        strm.next_in = reinterpret_cast<Bytef *>(inBuf.data());
        do {
            strm.avail_out = CHUNK;
            strm.next_out = reinterpret_cast<Bytef *>(outBuf.data());
            ret = inflate(&strm, Z_NO_FLUSH);
            if (ret == Z_STREAM_ERROR || ret == Z_NEED_DICT || ret == Z_DATA_ERROR || ret == Z_MEM_ERROR) {
                std::cerr << "Error: decompression failed.\n";
                inflateEnd(&strm);
                return -1;
            }
            size_t have = CHUNK - strm.avail_out;
            if (have > 0) {
                buffer.append(outBuf.data(), have);
                // Process complete lines using offset tracking (O(1) per line vs O(n) erase)
                size_t pos;
                while ((pos = buffer.find('\n', bufferOffset)) != std::string::npos) {
                    std::string_view line(buffer.data() + bufferOffset, pos - bufferOffset);
                    bufferOffset = pos + 1;
                    lineNumber++;
                    if (!processLineMmap(line, lineNumber, count)) {
                        inflateEnd(&strm);
                        return -1;
                    }
                }
                // Erase consumed data only when buffer gets large
                if (bufferOffset > 0 && bufferOffset > buffer.size() / 2) {
                    buffer.erase(0, bufferOffset);
                    bufferOffset = 0;
                }
            }
        } while (strm.avail_out == 0);
    } while (ret != Z_STREAM_END);

    // Process remaining data
    if (bufferOffset < buffer.size()) {
        std::string_view remaining(buffer.data() + bufferOffset, buffer.size() - bufferOffset);
        if (!remaining.empty()) {
            lineNumber++;
            if (!processLineMmap(remaining, lineNumber, count)) {
                inflateEnd(&strm);
                return -1;
            }
        }
    }
    inflateEnd(&strm);
    return count;
}

// Memory-mapped processLine - uses string_view for zero-copy
bool VCFXVariantCounter::processLineMmap(std::string_view line, int lineNumber, int &count) {
    if (line.empty())
        return true;
    if (line[0] == '#')
        return true;

    // Fast tab counting - no allocations needed
    int tabCount = 0;
    for (char c : line) {
        if (c == '\t') {
            tabCount++;
            if (tabCount >= 7) {
                count++;
                return true;
            }
        }
    }

    // Less than 8 columns
    if (strictMode) {
        std::cerr << "Error: line " << lineNumber << " has <8 columns.\n";
        return false;
    } else {
        std::cerr << "Warning: skipping line " << lineNumber << " with <8 columns.\n";
        return true;
    }
}

// Memory-mapped file counting - much faster than stdin for large files
int VCFXVariantCounter::countVariantsMmap(const char *filename) {
    int fd = open(filename, O_RDONLY);
    if (fd < 0) {
        std::cerr << "Error: cannot open file: " << filename << "\n";
        return -1;
    }

    struct stat st;
    if (fstat(fd, &st) < 0) {
        close(fd);
        std::cerr << "Error: cannot stat file: " << filename << "\n";
        return -1;
    }

    size_t fileSize = static_cast<size_t>(st.st_size);
    if (fileSize == 0) {
        close(fd);
        return 0;
    }

    void *mapped = mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, fd, 0);
    if (mapped == MAP_FAILED) {
        close(fd);
        std::cerr << "Error: cannot mmap file: " << filename << "\n";
        return -1;
    }

    // Advise kernel we'll read sequentially
    madvise(mapped, fileSize, MADV_SEQUENTIAL);

    const char *data = static_cast<const char *>(mapped);
    const char *end = data + fileSize;
    int count = 0;
    int lineNumber = 0;

    const char *lineStart = data;
    while (lineStart < end) {
        // Find end of line
        const char *lineEnd = lineStart;
        while (lineEnd < end && *lineEnd != '\n') {
            lineEnd++;
        }

        lineNumber++;
        std::string_view line(lineStart, lineEnd - lineStart);

        // Remove trailing \r if present (Windows line endings)
        if (!line.empty() && line.back() == '\r') {
            line = line.substr(0, line.size() - 1);
        }

        if (!processLineMmap(line, lineNumber, count)) {
            munmap(mapped, fileSize);
            close(fd);
            return -1;
        }

        lineStart = lineEnd + 1;
    }

    munmap(mapped, fileSize);
    close(fd);
    return count;
}

static void show_help() {
    VCFXVariantCounter obj;
    char arg0[] = "VCFX_variant_counter";
    char arg1[] = "--help";
    char *argv2[] = {arg0, arg1, nullptr};
    obj.run(2, argv2);
}

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_variant_counter", show_help))
        return 0;
    VCFXVariantCounter app;
    return app.run(argc, argv);
}
