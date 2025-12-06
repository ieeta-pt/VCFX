# VCFX_indexer

## Overview
`VCFX_indexer` is a high-performance utility tool for creating a byte-offset index of a VCF file. It generates a simple tab-delimited index file that maps chromosome and position to the exact byte offset in the original file, enabling efficient random access to variants without scanning the entire file. The index uses 64-bit integers for both the position and the byte offset so very large coordinates are fully supported.

## Usage

```bash
VCFX_indexer [OPTIONS] [input.vcf] > index.tsv
VCFX_indexer [OPTIONS] < input.vcf > index.tsv
```

## Options

| Option | Description |
|--------|-------------|
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

`VCFX_indexer` processes a VCF file to create a position-based index by:

1. Reading the VCF file from a file argument (fastest, uses mmap) or standard input
2. Locating the #CHROM header line to determine the start of data lines
3. For each variant line:
   - Extracting the chromosome (CHROM) and position (POS) values
   - Calculating the precise byte offset from the start of the file
4. Writing a three-column index to standard output with:
    - CHROM: The chromosome identifier from the VCF
    - POS: The position value from the VCF (stored as a 64-bit integer)
    - FILE_OFFSET: The byte offset to the start of the line in the source file (also 64-bit)

This index enables efficient random access to specific variants in large VCF files by allowing tools to seek directly to a byte offset rather than scanning the entire file. It's particularly useful for building tools that need to query specific regions of a VCF file. When a file path is provided as an argument, memory-mapped I/O is used for optimal performance on large files.

## Output Format
The index file is a tab-delimited text file with the following format:

```
CHROM   POS    FILE_OFFSET
1       100    542
1       200    621
1       300    702
```

Where:

- `CHROM` is the chromosome identifier from the VCF
- `POS` is the genomic position from the VCF (64-bit integer)
- `FILE_OFFSET` is the byte offset from the start of the VCF file (64-bit integer)

## Examples

### Basic Usage (File Argument - Fastest)

Create an index for a VCF file using memory-mapped I/O:
```bash
VCFX_indexer input.vcf > input.vcf.idx
```

### Standard Input Mode

Create an index from standard input:
```bash
VCFX_indexer < input.vcf > input.vcf.idx
```

### Using with Other Tools
Use the index to quickly extract a specific variant:

```bash
# Find the offset for position 1:12345
grep -P "^1\t12345\t" input.vcf.idx

# Use the offset (e.g., 23456) to seek directly to that variant
tail -c +23456 input.vcf | head -1
```

## Special Case Handling

### File Format Detection

- The tool automatically handles LF and CRLF line endings
- Byte offsets are calculated correctly regardless of the line ending style

### Malformed VCF Files

- Lines with unparseable position values are skipped
- If the #CHROM header is missing, an error is reported and no index entries are generated
- The tool requires variants to be tab-delimited; space-delimited files are not properly indexed

### Stream Processing

- The tool can process files from pipes as well as regular files
- For piped input, offsets are calculated relative to the beginning of the piped stream

### Empty Lines and Comments

- Empty lines and comment lines (starting with #) are properly handled and do not generate index entries

## Performance Considerations

- **File argument mode** uses memory-mapped I/O with SIMD-optimized parsing for maximum throughput
- Indexes 427,000 variants (4.3GB VCF) in ~6-9 seconds on modern hardware
- Uses zero-copy parsing with pointer arithmetic (no string allocations per line)
- 1MB output buffering reduces syscall overhead
- Memory usage is minimal: O(1) for file mode, O(line_length) for stdin mode
- The tool processes the VCF file in a single pass
- For very large files, the index itself will be much smaller than the original VCF
- The index file size scales with the number of variants, not the file size

## Limitations

- No support for indexing other fields besides CHROM and POS
- Does not validate the VCF format beyond basic column checking
- No built-in compression of the index file
- Cannot add new entries to an existing index (must regenerate the full index)
- Does not directly support query operations (must be used with other tools)
- Cannot handle compressed (gzipped) VCF files directly 
