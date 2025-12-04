#ifndef VCFX_CORE_H
#define VCFX_CORE_H

#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace vcfx {

// Trim leading and trailing whitespace from a string
std::string trim(const std::string &str);

// Split a string on the given delimiter
std::vector<std::string> split(const std::string &str, char delimiter);

// Convenience helpers for printing common messages
void print_error(const std::string &msg, std::ostream &os = std::cerr);
void print_version(const std::string &tool, const std::string &version, std::ostream &os = std::cout);

inline std::string get_version() {
#ifdef VCFX_VERSION
    return VCFX_VERSION;
#else
    return "unknown";
#endif
}

inline bool handle_version_flag(int argc, char *argv[], const std::string &tool, std::ostream &os = std::cout) {
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--version") == 0 || std::strcmp(argv[i], "-v") == 0) {
            print_version(tool, get_version(), os);
            return true;
        }
    }
    return false;
}

// Check if a specific flag (long or short form) is present
bool flag_present(int argc, char *argv[], const char *long_flag, const char *short_flag = nullptr);

// Handle the --help flag using the provided callback. Returns true if the flag
// was found and handled.
inline bool handle_help_flag(int argc, char *argv[], void (*print_help)()) {
    if (flag_present(argc, argv, "--help", "-h")) {
        if (print_help)
            print_help();
        return true;
    }
    return false;
}

// Handle both --help and --version flags. Returns true if either flag was found
// and processed (in which case the caller should exit).
inline bool handle_common_flags(int argc, char *argv[], const std::string &tool, void (*print_help)(),
                                std::ostream &os = std::cout) {
    if (handle_help_flag(argc, argv, print_help))
        return true;
    return handle_version_flag(argc, argv, tool, os);
}

// Read entire input stream, automatically decompressing if gzip/BGZF
// compressed. Returns true on success and stores the resulting text in
// 'out'.
bool read_maybe_compressed(std::istream &in, std::string &out);

// Convenience helper to read a file that may be gzip/BGZF compressed. The file
// is loaded completely into memory and stored in 'out'. Returns true on
// success.
bool read_file_maybe_compressed(const std::string &path, std::string &out);

// ------------------------------------------------------------
// StreamingGzipReader: Line-by-line reading with bounded memory
// ------------------------------------------------------------
// This class provides streaming decompression of gzip/BGZF files,
// reading one line at a time without loading the entire file into memory.
// Memory usage: O(chunk_size + line_length) ~ 64KB typical
//
// Usage:
//   std::ifstream file("data.vcf.gz", std::ios::binary);
//   vcfx::StreamingGzipReader reader(file);
//   std::string line;
//   while (reader.getline(line)) {
//       // process line
//   }
//
class StreamingGzipReader {
public:
    // Construct a streaming reader from an input stream
    // The stream should be opened in binary mode
    explicit StreamingGzipReader(std::istream &in);

    // Destructor - cleans up zlib resources
    ~StreamingGzipReader();

    // Non-copyable
    StreamingGzipReader(const StreamingGzipReader &) = delete;
    StreamingGzipReader &operator=(const StreamingGzipReader &) = delete;

    // Move constructible
    StreamingGzipReader(StreamingGzipReader &&other) noexcept;
    StreamingGzipReader &operator=(StreamingGzipReader &&other) noexcept;

    // Read the next line (without newline character)
    // Returns true if a line was read, false on EOF or error
    bool getline(std::string &line);

    // Check if the reader encountered an error
    bool error() const { return error_; }

    // Check if we've reached end of file
    bool eof() const { return eof_ && lineBuffer_.empty(); }

    // Check if the input was actually gzip compressed
    bool is_compressed() const { return isCompressed_; }

private:
    static constexpr size_t CHUNK_SIZE = 65536;  // 64KB chunks

    std::istream &in_;
    bool isCompressed_ = false;
    bool initialized_ = false;
    bool eof_ = false;
    bool error_ = false;

    // zlib stream (opaque pointer to avoid exposing zlib in header)
    void *zstrm_ = nullptr;

    // Input/output buffers
    char inBuf_[CHUNK_SIZE];
    char outBuf_[CHUNK_SIZE];
    size_t inBufPos_ = 0;
    size_t inBufLen_ = 0;

    // Line buffer for accumulating partial lines
    std::string lineBuffer_;

    // Initialize zlib stream
    bool initZlib();

    // Decompress more data into the output buffer
    // Returns number of bytes decompressed, 0 on EOF, -1 on error
    int decompressChunk();

    // Read more data from uncompressed stream
    bool readUncompressed();
};

// Helper function to create a streaming reader that auto-detects compression
// Returns a unique_ptr to StreamingGzipReader, or nullptr on error
std::unique_ptr<StreamingGzipReader> make_streaming_reader(std::istream &in);

// Helper function to create a streaming reader from a file path
// Handles both compressed and uncompressed files
std::unique_ptr<StreamingGzipReader> make_streaming_reader(const std::string &path, std::ifstream &fileStream);

} // namespace vcfx

#endif // VCFX_CORE_H
