# VCFX_missing_detector

## Overview

VCFX_missing_detector identifies and flags variants in a VCF file that contain missing genotype data in any sample. This tool helps researchers identify potentially problematic variants or samples with incomplete data that may require special handling in downstream analyses.

## Usage

```bash
VCFX_missing_detector [OPTIONS] [input.vcf]
VCFX_missing_detector [OPTIONS] < input.vcf > flagged.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `-i`, `--input FILE` | Input VCF file (uses memory-mapping for best performance) |
| `-t`, `--threads N` | Number of threads for parallel scanning (default: auto) |
| `-q`, `--quiet` | Suppress informational messages |
| `-h`, `--help` | Display help message and exit |
| `-v`, `--version` | Show program version and exit |

## Description

VCFX_missing_detector analyzes a VCF file to identify variants with missing genotype data. The tool:

1. Reads a VCF file using memory-mapped I/O for extreme performance
2. Uses SIMD-accelerated scanning to detect '.' characters in sample columns
3. Multi-threaded pre-scan determines if any missing genotypes exist
4. For each variant with missing data, adds `MISSING_GENOTYPES=1` to the INFO field
5. Uses zero-copy output for unchanged lines

This simple annotation allows researchers to easily:
- Filter variants based on missing data presence using standard VCF tools
- Identify data completeness issues that might affect analysis results
- Implement different handling strategies for variants with missing data

## Output Format

The output is a valid VCF file with the same format as the input, but with an additional INFO field annotation for variants containing missing genotypes:

```
MISSING_GENOTYPES=1
```

This annotation is appended to the existing INFO field, or replaces the `.` placeholder if the INFO field is empty.

## Examples

### File Input Mode (Recommended)

Use memory-mapped file I/O for best performance:

```bash
VCFX_missing_detector -i input.vcf > flagged.vcf
VCFX_missing_detector input.vcf > flagged.vcf
```

### Basic Usage (Stdin)

```bash
VCFX_missing_detector < input.vcf > flagged.vcf
```

### In a Pipeline with Filtering

```bash
# Flag variants with missing genotypes, then filter to keep only complete variants
VCFX_missing_detector -i input.vcf | grep -v "MISSING_GENOTYPES=1" > complete_variants.vcf
```

### Counting Missing Variants

```bash
# Count variants with missing genotypes
VCFX_missing_detector -i input.vcf | grep "MISSING_GENOTYPES=1" | wc -l
```

### Quiet Mode for Scripts

```bash
# Suppress progress messages
VCFX_missing_detector -q -i input.vcf > flagged.vcf
```

## Missing Genotype Detection

The tool uses comprehensive logic to identify various forms of missing genotype data:

1. **Completely missing genotypes**: Formats like `./.`, `.|.`, or just `.`
2. **Partially missing diploid genotypes**: When one allele is missing, like `./1` or `0/.`
3. **Multi-field format handling**: Properly extracts just the GT portion when other fields (like DP, GQ) are present
4. **Format field awareness**: Correctly identifies the GT position in the FORMAT string

The tool examines each sample column independently and flags a variant if any sample has missing data.

## Handling Special Cases

The tool implements several strategies for handling edge cases:

1. **Missing FORMAT field**: If GT is not included in the FORMAT column, the variant is passed through unchanged
2. **No sample columns**: Variants with fewer than 9 columns are passed through unchanged
3. **Empty INFO field**: If the original INFO is "." (missing), it's replaced with "MISSING_GENOTYPES=1"
4. **Non-empty INFO field**: The missing flag is appended with a semicolon separator
5. **Empty lines**: Preserved with a single newline
6. **Header lines**: Passed through unchanged
7. **Non-diploid genotypes**: The tool focuses on diploid genotypes with a single delimiter ('/' or '|')

## Performance Characteristics

The tool uses several optimizations for extreme performance:

- **Memory-mapped I/O**: Uses `mmap()` for file input to minimize syscall overhead
- **SIMD acceleration**: Uses NEON/SSE2/AVX2 instructions for fast '.' character scanning
- **Multi-threaded pre-scan**: Parallel scanning to determine if file has any missing genotypes
- **Zero-copy fast path**: Files with no missing genotypes are copied verbatim without per-line parsing
- **Batched zero-copy output**: Groups consecutive unchanged lines into single `write()` calls

### Benchmark Results (4GB VCF, chr21, 427K variants, 2504 samples)

| Mode | Time | Notes |
|------|------|-------|
| mmap (-i) | ~7s | Multi-threaded pre-scan + zero-copy output |
| stdin | ~5 min | Buffered I/O |

**Speedup: ~42x** with memory-mapped I/O

## Limitations

1. Primarily designed for diploid genotypes; may not correctly identify missing data in haploid or polyploid contexts
2. Limited to checking the GT field; does not evaluate other potential indicators of missing data
3. No built-in functionality to annotate the percentage or count of samples with missing data
4. No option to customize the INFO field tag name from the default "MISSING_GENOTYPES"
5. Cannot perform sample-specific missing data analysis
6. No threshold options (e.g., flag only if more than X% of samples have missing data)
7. Limited to binary detection (missing/not missing) without quantifying the degree of missingness
