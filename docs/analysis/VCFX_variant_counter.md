# VCFX_variant_counter

## Overview
`VCFX_variant_counter` is a high-performance utility tool that counts the total number of valid variants (data lines) in a VCF file. It supports both file arguments (using memory-mapped I/O for maximum speed) and standard input, and outputs the total count of valid variant records.

## Usage
```bash
VCFX_variant_counter [OPTIONS] [input.vcf]
VCFX_variant_counter [OPTIONS] < input.vcf
```

## Options
| Option | Description |
|--------|-------------|
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |
| `-s`, `--strict` | Fail on any data line with fewer than 8 columns |

## Description
`VCFX_variant_counter` processes a VCF file by:

1. Reading the VCF file from a file argument (fastest, uses mmap) or standard input
2. Ignoring all header lines (those starting with `#`)
3. For each data line:
   - Checking if it has at least 8 columns (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)
   - If it has 8 or more columns, counting it as a valid variant
   - If it has fewer than 8 columns:
     - In strict mode: exiting with an error
     - In normal mode: skipping the line with a warning
4. Finally, printing the total count of valid variants

This tool is useful for quickly determining the number of variants in a VCF file, which can be helpful for quality control, workflow validation, or simply getting an overview of a dataset's size. When a file path is provided as an argument, memory-mapped I/O is used for optimal performance on large files.

## VCF Format Requirements
The tool assumes a standard VCF format where:
- Header lines start with `#`
- Data lines have at least 8 tab-separated columns (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)
- Optional FORMAT and sample columns may follow

## Examples

### Basic Usage (File Argument - Fastest)
Count the variants in a VCF file using memory-mapped I/O:
```bash
VCFX_variant_counter input.vcf
```
Output:
```
Total Variants: 1234
```

### Standard Input Mode
Count variants from standard input:
```bash
VCFX_variant_counter < input.vcf
```

### Using Strict Mode
Count variants with strict validation (will fail on malformed lines):
```bash
VCFX_variant_counter --strict input.vcf
```

### Using in a Pipeline
Count variants after filtering:
```bash
cat input.vcf | grep -v "FILTER=FAIL" | VCFX_variant_counter
```

### Compressed Files
Gzip-compressed files are automatically detected when piped via stdin:
```bash
cat input.vcf.gz | VCFX_variant_counter
zcat input.vcf.gz | VCFX_variant_counter
```

## Error Handling

### Invalid Lines
- In normal mode (default):
  - Lines with fewer than 8 columns are skipped
  - A warning is printed to standard error for each skipped line
  - The count continues with valid lines

- In strict mode (`--strict`):
  - If any line has fewer than 8 columns, the program exits with an error
  - The error message includes the line number
  - The exit code is 1 to indicate failure

### Empty Lines
Empty lines are ignored and not counted.

## Performance Considerations
- **File argument mode** uses memory-mapped I/O with SIMD-optimized parsing for maximum throughput
- Processes 427,000 variants (4.3GB VCF) in ~2 seconds on modern hardware
- Uses `memchr()` for SIMD-accelerated column validation (no character-by-character parsing)
- Minimal memory requirements: O(1) for file mode, O(line_length) for stdin mode
- Performance scales linearly with the number of lines in the input file
- No external dependencies or reference files are required
- Automatic gzip detection and streaming decompression for compressed stdin input

## Limitations
- The tool only checks the number of columns, not their content
- It doesn't validate if the VCF follows the specification for field formats
- File argument mode does not support compressed files (pipe through stdin instead)
- No detailed reporting (e.g., breakdown by chromosome or variant type)
