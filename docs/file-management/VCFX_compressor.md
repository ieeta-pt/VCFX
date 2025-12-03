# VCFX_compressor

## Overview

VCFX_compressor provides simple compression and decompression functionality for VCF files using the zlib library. It allows users to compress VCF files for storage or transfer, and decompress them for analysis.

## Usage

```bash
VCFX_compressor [OPTIONS] < input_file > output_file
```

## Options

| Option | Description |
|--------|-------------|
| `-c`, `--compress` | Compress the input VCF file (read from stdin, write to stdout) |
| `-d`, `--decompress` | Decompress the input VCF.gz file (read from stdin, write to stdout) |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_compressor is a straightforward utility that enables compression and decompression of VCF files using zlib's DEFLATE algorithm. The tool:

1. Reads input data from standard input (stdin)
2. Processes the data in memory-efficient chunks
3. Applies compression or decompression based on the specified mode
4. Writes the processed data to standard output (stdout)

The compression mode produces output compatible with gzip, while the decompression mode can handle standard gzip-compressed files. This makes the tool interoperable with widely used genomics software that expects gzip-compressed VCF files.

## Output Format

The output format depends on the chosen mode:

- **Compression mode**: Produces a gzip-compatible compressed binary file
- **Decompression mode**: Produces a plain text VCF file

## Examples

### Compressing a VCF File

```bash
# Basic compression
./VCFX_compressor --compress < input.vcf > output.vcf.gz

# Compress and view file size reduction
./VCFX_compressor --compress < input.vcf > output.vcf.gz
echo "Original size: $(wc -c < input.vcf) bytes"
echo "Compressed size: $(wc -c < output.vcf.gz) bytes"
```

### Decompressing a VCF File

```bash
# Basic decompression
./VCFX_compressor --decompress < input.vcf.gz > output.vcf

# Decompress for analysis
./VCFX_compressor --decompress < input.vcf.gz | head -n 20
```

### In a Pipeline

```bash
# Filter a VCF file, compress it, then decompress for viewing
cat input.vcf | grep -v "^#" | grep "PASS" | ./VCFX_compressor --compress > filtered.vcf.gz
./VCFX_compressor --decompress < filtered.vcf.gz | head
```

## Data Processing

VCFX_compressor processes data in chunks to maintain memory efficiency. The default chunk size is 16KB, which provides a good balance between memory usage and processing efficiency. The tool:

1. Reads input data in 16KB chunks
2. Processes each chunk using zlib's compression/decompression functions
3. Writes processed data to output immediately as each chunk is completed
4. Continues until all input data has been processed

This streaming approach allows the tool to handle files of any size without loading the entire file into memory.

## Handling Special Cases

The tool implements several strategies for handling edge cases:

1. **Empty files**: Properly handles empty input, producing valid empty output
2. **Truncated inputs**: When decompressing, detects and warns about truncated or incomplete compressed data
3. **Invalid compressed data**: Reports errors when attempting to decompress invalid or corrupted data
4. **I/O errors**: Provides error messages for issues with reading input or writing output
5. **Incorrect usage**: Enforces mutually exclusive selection of compression or decompression mode

## Performance

VCFX_compressor is designed for efficiency:

1. Processes data in chunks, maintaining a low and consistent memory footprint
2. Uses zlib's optimized compression/decompression algorithms
3. Avoids unnecessary memory copying or buffering of the entire file
4. Provides reasonable compression ratios typical of gzip compression
5. Handles large files efficiently due to its streaming architecture

## Limitations

1. **Not BGZF compatible**: Does not produce block-gzipped format required for indexed access via tabix
2. **No compression level control**: Uses zlib's default compression level with no user-configurable options
3. **Single-threaded**: Does not utilize multi-threading for potentially faster processing
4. **No integrity verification**: Does not verify the integrity of decompressed data
5. **Limited format support**: Only handles gzip compression, not other formats like bzip2 or xz
6. **No indexing support**: Does not maintain or generate indices for compressed files
7. **Standard I/O only**: Cannot directly specify input and output filenames (uses stdin/stdout) 
