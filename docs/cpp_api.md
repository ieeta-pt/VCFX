# C++ API

VCFX provides a core C++ library (`vcfx_core`) with utility functions and
classes for VCF file processing. The library is designed for efficient
handling of large VCF files with minimal memory usage.

## Header

Include the main header to access all functionality:

```cpp
#include "vcfx_core.h"
```

## Namespace

All functions and classes are in the `vcfx` namespace.

## String Utilities

### trim

Remove leading and trailing whitespace from a string.

```cpp
std::string trim(const std::string &str);
```

**Example:**
```cpp
std::string s = vcfx::trim("  hello world  ");
// s == "hello world"
```

### split

Split a string by a delimiter character.

```cpp
std::vector<std::string> split(const std::string &str, char delimiter);
```

**Example:**
```cpp
auto parts = vcfx::split("chr1\t100\tA\tG", '\t');
// parts == {"chr1", "100", "A", "G"}
```

## File Reading

### read_maybe_compressed

Read from an input stream, automatically decompressing if gzip/BGZF compressed.

```cpp
bool read_maybe_compressed(std::istream &in, std::string &out);
```

**Parameters:**

- `in` - Input stream (should be opened in binary mode for compressed data)
- `out` - String to store the decompressed/read content

**Returns:** `true` on success, `false` on error.

**Note:** This function loads the entire file into memory. For large files,
use `StreamingGzipReader` instead.

### read_file_maybe_compressed

Read a file by path, automatically decompressing if needed.

```cpp
bool read_file_maybe_compressed(const std::string &path, std::string &out);
```

**Parameters:**

- `path` - Path to the file (compression detected by extension or magic bytes)
- `out` - String to store the file contents

**Returns:** `true` on success, `false` on error.

**Recognized extensions:** `.gz`, `.bgz`, `.bgzf`

## Streaming Decompression

### StreamingGzipReader

A class for reading gzip/BGZF compressed files line-by-line with bounded
memory usage. This is the recommended approach for processing large VCF files.

**Memory usage:** O(chunk_size + line_length) ~ 64KB typical

```cpp
class StreamingGzipReader {
public:
    explicit StreamingGzipReader(std::istream &in);
    ~StreamingGzipReader();

    // Non-copyable, move-only
    StreamingGzipReader(const StreamingGzipReader &) = delete;
    StreamingGzipReader &operator=(const StreamingGzipReader &) = delete;
    StreamingGzipReader(StreamingGzipReader &&other) noexcept;
    StreamingGzipReader &operator=(StreamingGzipReader &&other) noexcept;

    bool getline(std::string &line);
    bool error() const;
    bool eof() const;
    bool is_compressed() const;
};
```

**Methods:**

| Method | Description |
|--------|-------------|
| `getline(line)` | Read the next line (without newline). Returns `true` if a line was read. |
| `error()` | Returns `true` if an error occurred during reading. |
| `eof()` | Returns `true` if end of file has been reached. |
| `is_compressed()` | Returns `true` if the input was gzip/BGZF compressed. |

**Example:**

```cpp
#include "vcfx_core.h"
#include <fstream>
#include <iostream>

int main() {
    std::ifstream file("data.vcf.gz", std::ios::binary);
    vcfx::StreamingGzipReader reader(file);

    std::string line;
    int variant_count = 0;

    while (reader.getline(line)) {
        if (!line.empty() && line[0] != '#') {
            variant_count++;
        }
    }

    if (reader.error()) {
        std::cerr << "Error reading file\n";
        return 1;
    }

    std::cout << "Variants: " << variant_count << "\n";
    return 0;
}
```

### make_streaming_reader

Factory functions for creating `StreamingGzipReader` instances.

```cpp
std::unique_ptr<StreamingGzipReader> make_streaming_reader(std::istream &in);

std::unique_ptr<StreamingGzipReader> make_streaming_reader(
    const std::string &path,
    std::ifstream &fileStream
);
```

**Parameters:**

- `in` - Input stream opened in binary mode
- `path` - File path to open
- `fileStream` - Reference to an ifstream that will be opened by the function

**Returns:** `unique_ptr` to the reader, or `nullptr` on error.

**Example with file path:**

```cpp
std::ifstream fileStream;
auto reader = vcfx::make_streaming_reader("large_file.vcf.gz", fileStream);

if (!reader) {
    std::cerr << "Failed to open file\n";
    return 1;
}

std::string line;
while (reader->getline(line)) {
    // process line
}
```

## Command-Line Helpers

### flag_present

Check if a command-line flag is present.

```cpp
bool flag_present(int argc, char *argv[],
                  const char *long_flag,
                  const char *short_flag = nullptr);
```

### handle_common_flags

Handle `--help` and `--version` flags automatically.

```cpp
bool handle_common_flags(int argc, char *argv[],
                         const std::string &tool,
                         void (*print_help)(),
                         std::ostream &os = std::cout);
```

**Returns:** `true` if a flag was handled (caller should exit).

### get_version

Get the VCFX version string.

```cpp
std::string get_version();
```

## Error Handling

### print_error

Print an error message to a stream.

```cpp
void print_error(const std::string &msg, std::ostream &os = std::cerr);
```

### print_version

Print version information.

```cpp
void print_version(const std::string &tool,
                   const std::string &version,
                   std::ostream &os = std::cout);
```

## Compression Support

The library supports:

- **Gzip** (.gz) - Standard gzip compression
- **BGZF** (.bgz, .bgzf) - Block gzip format used by samtools/bcftools
- **Uncompressed** - Plain text files

Compression is auto-detected using magic bytes (`0x1f 0x8b` for gzip),
so files don't need to have the correct extension.

## Thread Safety

- String utilities (`trim`, `split`) are thread-safe
- `StreamingGzipReader` instances are not thread-safe; use one per thread
- File reading functions are thread-safe if operating on different files

## Linking

Link against `vcfx_core` and zlib:

```cmake
target_link_libraries(your_target vcfx_core z)
```

Or with pkg-config:

```bash
g++ -o myapp myapp.cpp $(pkg-config --libs vcfx) -lz
```
