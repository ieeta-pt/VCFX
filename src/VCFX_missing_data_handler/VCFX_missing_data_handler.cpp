#include "VCFX_missing_data_handler.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <thread>
#include <unistd.h>
#include <vector>

// ============================================================================
// RADICALLY DIFFERENT APPROACH: Memory-mapped, zero-copy, multi-threaded
// ============================================================================
// Key insight: Most lines have NO missing data. Don't parse them at all.
// Strategy:
//   1. Memory-map the input file (eliminates read syscalls)
//   2. Scan raw bytes for missing patterns (./. or .\t or .\n)
//   3. Only if pattern found, do minimal modification
//   4. Multi-thread the processing
// ============================================================================

void printHelp() {
    std::cout
        << "VCFX_missing_data_handler\n"
        << "Usage: VCFX_missing_data_handler [OPTIONS] [files...]\n\n"
        << "Options:\n"
        << "  --fill-missing, -f            Impute missing genotypes with a default value (e.g., ./.).\n"
        << "  --default-genotype, -d GEN    Specify the default genotype for imputation (default: ./.).\n"
        << "  --threads, -t NUM             Number of threads (default: auto)\n"
        << "  --help, -h                    Display this help message and exit.\n\n";
}

// ============================================================================
// SIMD-like fast pattern scan using memchr
// ============================================================================

// Check if a position contains a missing genotype pattern
// Patterns: "./." or ".|." or "." followed by \t or : or \n
static inline bool isMissingPatternAt(const char* p, const char* end) {
    if (p >= end) return false;

    if (*p != '.') return false;

    // Check "." alone (followed by delimiter)
    if (p + 1 >= end || p[1] == '\t' || p[1] == ':' || p[1] == '\n') {
        return true;
    }

    // Check "./." or ".|."
    if (p + 2 < end && (p[1] == '/' || p[1] == '|') && p[2] == '.') {
        // Make sure it's the whole genotype (followed by delimiter or end)
        if (p + 3 >= end || p[3] == '\t' || p[3] == ':' || p[3] == '\n') {
            return true;
        }
    }

    return false;
}

// Fast scan for any '.' character in sample region
static inline const char* findDotInSamples(const char* start, const char* end) {
    return (const char*)memchr(start, '.', end - start);
}

// Find the GT field position and check if it's missing
// Returns: -1 if no GT field found, -2 if GT exists but NOT missing,
//          position (>=0) if GT IS missing
static inline int findMissingGTPosition(const char* sampleStart, const char* sampleEnd, int gtIndex) {
    if (gtIndex < 0) return -1;

    const char* p = sampleStart;
    int fieldIdx = 0;
    const char* fieldStart = p;

    while (p <= sampleEnd) {
        if (p == sampleEnd || *p == ':') {
            if (fieldIdx == gtIndex) {
                // Check if this field is missing
                size_t fieldLen = p - fieldStart;
                if (fieldLen == 1 && *fieldStart == '.') {
                    return fieldStart - sampleStart; // Position of missing GT (could be 0)
                }
                if (fieldLen == 3 && fieldStart[0] == '.' &&
                    (fieldStart[1] == '/' || fieldStart[1] == '|') &&
                    fieldStart[2] == '.') {
                    return fieldStart - sampleStart; // Position of missing GT (could be 0)
                }
                return -2; // GT exists but not missing
            }
            fieldIdx++;
            fieldStart = p + 1;
        }
        p++;
    }
    return -1;
}

// ============================================================================
// Process a single line with minimal copying
// Returns true if line was modified
// ============================================================================
static bool processLineZeroCopy(const char* lineStart, size_t lineLen,
                                 int gtIndex, const std::string& replacement,
                                 std::string& output) {
    // Skip empty lines and headers
    if (lineLen == 0) {
        output.assign(lineStart, lineLen);
        output.push_back('\n');
        return false;
    }

    if (*lineStart == '#') {
        output.assign(lineStart, lineLen);
        output.push_back('\n');
        return false;
    }

    const char* lineEnd = lineStart + lineLen;

    // Quick scan: does this line even have a '.' character after the FORMAT field?
    // Find the 9th tab (start of samples)
    int tabCount = 0;
    const char* sampleRegionStart = nullptr;
    for (const char* p = lineStart; p < lineEnd; p++) {
        if (*p == '\t') {
            tabCount++;
            if (tabCount == 9) {
                sampleRegionStart = p + 1;
                break;
            }
        }
    }

    if (!sampleRegionStart || tabCount < 9) {
        // Not enough fields, pass through
        output.assign(lineStart, lineLen);
        output.push_back('\n');
        return false;
    }

    // FAST PATH: Quick scan for '.' in sample region
    const char* dotPos = findDotInSamples(sampleRegionStart, lineEnd);
    if (!dotPos) {
        // No dots at all - definitely no missing data, pass through unchanged
        output.assign(lineStart, lineLen);
        output.push_back('\n');
        return false;
    }

    // There's at least one dot - need to process samples
    // But still try to minimize work

    output.clear();
    output.reserve(lineLen + 100); // Reserve space for potential expansion

    // Copy everything up to sample region
    output.append(lineStart, sampleRegionStart - lineStart);

    bool modified = false;
    const char* sampleStart = sampleRegionStart;

    while (sampleStart < lineEnd) {
        // Find end of this sample
        const char* sampleEnd = sampleStart;
        while (sampleEnd < lineEnd && *sampleEnd != '\t') {
            sampleEnd++;
        }

        size_t sampleLen = sampleEnd - sampleStart;

        // Check if this specific sample has missing GT
        // Returns >= 0 if missing (position), -1 or -2 if not missing
        int missingPos = findMissingGTPosition(sampleStart, sampleEnd, gtIndex);

        if (missingPos >= 0) {
            // GT is missing at position missingPos within sample
            // Copy up to the missing position
            output.append(sampleStart, missingPos);

            // Insert replacement
            output.append(replacement);

            // Figure out what we skipped (./. or . or .|.)
            const char* gtStart = sampleStart + missingPos;
            size_t skipLen = 1; // At least skip the '.'
            if (gtStart + 2 < sampleEnd &&
                (gtStart[1] == '/' || gtStart[1] == '|') &&
                gtStart[2] == '.') {
                skipLen = 3;
            }

            // Copy the rest of the sample
            output.append(gtStart + skipLen, sampleEnd - (gtStart + skipLen));
            modified = true;
        } else {
            // No missing GT, copy sample as-is
            output.append(sampleStart, sampleLen);
        }

        // Add tab if not at end
        if (sampleEnd < lineEnd) {
            output.push_back('\t');
            sampleStart = sampleEnd + 1;
        } else {
            break;
        }
    }

    output.push_back('\n');
    return modified;
}

// ============================================================================
// Multi-threaded line processor
// ============================================================================
struct LineRange {
    size_t startOffset;
    size_t endOffset;
};

static void processChunk(const char* data, const std::vector<LineRange>& lines,
                         size_t startIdx, size_t endIdx,
                         int gtIndex, const std::string& replacement,
                         std::string& output) {
    std::string lineOutput;
    lineOutput.reserve(32768);
    output.clear();
    output.reserve((endIdx - startIdx) * 10000); // Estimate

    for (size_t i = startIdx; i < endIdx; i++) {
        const LineRange& lr = lines[i];
        processLineZeroCopy(data + lr.startOffset,
                           lr.endOffset - lr.startOffset,
                           gtIndex, replacement, lineOutput);
        output.append(lineOutput);
    }
}

// ============================================================================
// Main processing function with memory mapping and multi-threading
// ============================================================================
static bool processVCFMapped(const char* filename, std::ostream& out,
                              bool fillMissing, const std::string& defaultGT,
                              int numThreads) {
    // Open and memory-map the file
    int fd = open(filename, O_RDONLY);
    if (fd < 0) {
        std::cerr << "Error: cannot open file " << filename << "\n";
        return false;
    }

    struct stat sb;
    if (fstat(fd, &sb) < 0) {
        close(fd);
        std::cerr << "Error: Cannot stat file\n";
        return false;
    }

    size_t fileSize = sb.st_size;
    if (fileSize == 0) {
        close(fd);
        return true;
    }

    // Memory map the file
    char* data = (char*)mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, fd, 0);
    if (data == MAP_FAILED) {
        close(fd);
        std::cerr << "Error: Cannot mmap file\n";
        return false;
    }

    // Advise kernel for sequential access
    madvise(data, fileSize, MADV_SEQUENTIAL);

    // Phase 1: Find all line boundaries (single-threaded, very fast)
    std::vector<LineRange> lines;
    lines.reserve(500000); // Pre-allocate for typical VCF

    size_t lineStart = 0;
    int gtIndex = -1;
    bool headerDone = false;

    for (size_t i = 0; i < fileSize; i++) {
        if (data[i] == '\n') {
            LineRange lr = {lineStart, i};
            lines.push_back(lr);

            // Parse header to find GT index
            if (!headerDone && lineStart < i && data[lineStart] == '#') {
                // Check for #CHROM line
                if (i - lineStart > 6 && strncmp(data + lineStart, "#CHROM", 6) == 0) {
                    headerDone = true;
                }
            } else if (!headerDone && lineStart < i && data[lineStart] != '#') {
                // First data line - find FORMAT field to get GT index
                int tabCount = 0;
                for (size_t j = lineStart; j < i; j++) {
                    if (data[j] == '\t') {
                        tabCount++;
                        if (tabCount == 8) {
                            // Next field is FORMAT
                            size_t formatStart = j + 1;
                            size_t formatEnd = formatStart;
                            while (formatEnd < i && data[formatEnd] != '\t') formatEnd++;

                            // Find GT in format
                            int idx = 0;
                            size_t fieldStart = formatStart;
                            for (size_t k = formatStart; k <= formatEnd; k++) {
                                if (k == formatEnd || data[k] == ':') {
                                    if (k - fieldStart == 2 &&
                                        data[fieldStart] == 'G' &&
                                        data[fieldStart + 1] == 'T') {
                                        gtIndex = idx;
                                    }
                                    idx++;
                                    fieldStart = k + 1;
                                }
                            }
                            break;
                        }
                    }
                }
                headerDone = true;
            }

            lineStart = i + 1;
        }
    }

    // Handle last line without newline
    if (lineStart < fileSize) {
        lines.push_back({lineStart, fileSize});
    }

    if (!fillMissing || gtIndex < 0) {
        // Just copy the file through
        out.write(data, fileSize);
        munmap(data, fileSize);
        close(fd);
        return true;
    }

    // Phase 2: Process lines in parallel
    size_t numLines = lines.size();

    if (numThreads <= 0) {
        numThreads = std::thread::hardware_concurrency();
        if (numThreads == 0) numThreads = 4;
    }

    // For small files, don't bother with threading
    if (numLines < 10000) {
        numThreads = 1;
    }

    std::vector<std::string> outputs(numThreads);
    std::vector<std::thread> threads;

    size_t linesPerThread = (numLines + numThreads - 1) / numThreads;

    for (int t = 0; t < numThreads; t++) {
        size_t startIdx = t * linesPerThread;
        size_t endIdx = std::min(startIdx + linesPerThread, numLines);

        if (startIdx >= numLines) break;

        threads.emplace_back(processChunk, data, std::ref(lines),
                            startIdx, endIdx, gtIndex, std::ref(defaultGT),
                            std::ref(outputs[t]));
    }

    // Wait for all threads
    for (auto& t : threads) {
        t.join();
    }

    // Write outputs in order
    for (int t = 0; t < numThreads && t < (int)outputs.size(); t++) {
        if (!outputs[t].empty()) {
            out.write(outputs[t].data(), outputs[t].size());
        }
    }

    munmap(data, fileSize);
    close(fd);
    return true;
}

// ============================================================================
// Fallback for stdin processing (can't mmap)
// ============================================================================
static bool processVCFStream(std::istream& in, std::ostream& out,
                              bool fillMissing, const std::string& defaultGT) {
    std::string line;
    int gtIndex = -1;
    bool foundFormat = false;

    // Large output buffer
    static const size_t OUTPUT_BUFFER_SIZE = 4 * 1024 * 1024; // 4MB buffer
    std::string outputBuffer;
    outputBuffer.reserve(OUTPUT_BUFFER_SIZE);

    std::string processedLine;
    processedLine.reserve(32768);

    while (std::getline(in, line)) {
        if (line.empty()) {
            outputBuffer.push_back('\n');
            continue;
        }

        if (line[0] == '#') {
            outputBuffer.append(line);
            outputBuffer.push_back('\n');
            if (outputBuffer.size() >= OUTPUT_BUFFER_SIZE - 65536) {
                out.write(outputBuffer.data(), outputBuffer.size());
                outputBuffer.clear();
            }
            continue;
        }

        // Find GT index from first data line
        if (!foundFormat) {
            int tabCount = 0;
            for (size_t i = 0; i < line.size(); i++) {
                if (line[i] == '\t') {
                    tabCount++;
                    if (tabCount == 8) {
                        size_t formatStart = i + 1;
                        size_t formatEnd = line.find('\t', formatStart);
                        if (formatEnd == std::string::npos) formatEnd = line.size();

                        int idx = 0;
                        size_t fieldStart = formatStart;
                        for (size_t k = formatStart; k <= formatEnd; k++) {
                            if (k == formatEnd || line[k] == ':') {
                                if (k - fieldStart == 2 &&
                                    line[fieldStart] == 'G' &&
                                    line[fieldStart + 1] == 'T') {
                                    gtIndex = idx;
                                }
                                idx++;
                                fieldStart = k + 1;
                            }
                        }
                        break;
                    }
                }
            }
            foundFormat = true;
        }

        if (!fillMissing || gtIndex < 0) {
            outputBuffer.append(line);
            outputBuffer.push_back('\n');
        } else {
            processLineZeroCopy(line.c_str(), line.size(), gtIndex, defaultGT, processedLine);
            outputBuffer.append(processedLine);
        }

        if (outputBuffer.size() >= OUTPUT_BUFFER_SIZE - 65536) {
            out.write(outputBuffer.data(), outputBuffer.size());
            outputBuffer.clear();
        }
    }

    if (!outputBuffer.empty()) {
        out.write(outputBuffer.data(), outputBuffer.size());
    }

    return true;
}

// ============================================================================
// Argument parsing
// ============================================================================
bool parseArguments(int argc, char* argv[], Arguments& args) {
    static struct option long_opts[] = {
        {"fill-missing", no_argument, 0, 'f'},
        {"default-genotype", required_argument, 0, 'd'},
        {"threads", required_argument, 0, 't'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    args.numThreads = 0; // Auto-detect

    while (true) {
        int c = getopt_long(argc, argv, "fd:t:h", long_opts, nullptr);
        if (c == -1) break;

        switch (c) {
            case 'f':
                args.fill_missing = true;
                break;
            case 'd':
                args.default_genotype = optarg;
                break;
            case 't':
                args.numThreads = std::stoi(optarg);
                break;
            case 'h':
            default:
                printHelp();
                exit(0);
        }
    }

    while (optind < argc) {
        args.input_files.push_back(argv[optind++]);
    }

    return true;
}

// ============================================================================
// Main entry point
// ============================================================================
bool handleMissingDataAll(const Arguments& args) {
    if (args.input_files.empty()) {
        // Read from stdin - can't use mmap
        return processVCFStream(std::cin, std::cout, args.fill_missing, args.default_genotype);
    } else {
        // Process files with mmap + multi-threading
        for (const auto& path : args.input_files) {
            if (!processVCFMapped(path.c_str(), std::cout,
                                  args.fill_missing, args.default_genotype,
                                  args.numThreads)) {
                return false;
            }
        }
        return true;
    }
}

static void show_help() { printHelp(); }

int main(int argc, char* argv[]) {
    vcfx::init_io();
    if (vcfx::handle_common_flags(argc, argv, "VCFX_missing_data_handler", show_help))
        return 0;

    Arguments args;
    parseArguments(argc, argv, args);

    if (args.fill_missing) {
        std::cerr << "Info: Missing genotypes will be imputed with: " << args.default_genotype << "\n";
        std::cerr << "Info: Using " << (args.numThreads > 0 ? args.numThreads : (int)std::thread::hardware_concurrency()) << " threads\n";
    }

    return handleMissingDataAll(args) ? 0 : 1;
}
