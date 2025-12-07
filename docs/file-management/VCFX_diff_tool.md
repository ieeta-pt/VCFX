# VCFX_diff_tool

## Overview

VCFX_diff_tool compares two VCF files and identifies variants that are unique to each file, providing a simple way to detect differences between variant sets.

**New in v1.2**: Major performance update with **memory-mapped I/O** and **SIMD-accelerated parsing** delivering ~50x speedup. The tool now processes 6.8GB files in under 60 seconds.

**New in v1.1**: The tool now supports **streaming mode** for handling pre-sorted files with O(1) memory usage, enabling comparison of 50GB+ files without loading them entirely into memory.

## Usage

```bash
VCFX_diff_tool --file1 <file1.vcf> --file2 <file2.vcf> [options]
```

## Options

| Option | Description |
|--------|-------------|
| `-a`, `--file1 <FILE>` | Required. Path to the first VCF file |
| `-b`, `--file2 <FILE>` | Required. Path to the second VCF file |
| `-s`, `--assume-sorted` | Assume input files are already sorted by (CHROM, POS). Enables streaming mode with O(1) memory for large files. |
| `-n`, `--natural-chr` | Use natural chromosome ordering (chr1 < chr2 < chr10) instead of lexicographic ordering |
| `-q`, `--quiet` | Suppress warning messages (useful for batch processing) |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_diff_tool analyzes two VCF files, compares their variant content, and identifies differences between them. The tool:

1. Loads variants from both VCF files, ignoring header lines
2. Creates a normalized key for each variant based on chromosome, position, reference allele, and sorted alternate alleles
3. Identifies variants that are unique to each file by comparing these keys
4. Reports the differences in a readable format

The tool offers two comparison strategies:

### Default Mode (In-Memory)
- Loads all variants from both files into hash sets
- Uses O(1) hash lookups for efficient comparison
- Works with unsorted input files
- Memory usage scales with number of variants

### Streaming Mode (`--assume-sorted`)
- Performs two-pointer merge diff
- Only keeps current line from each file in memory
- Memory usage: O(1) regardless of file sizes
- **Requires pre-sorted input files** for correct results
- Enables comparing files larger than available RAM

This tool is particularly useful for:
- Validating VCF file transformations
- Checking tool outputs against expected results
- Comparing variant calls between different callers or pipelines
- Verifying that VCF manipulations haven't inadvertently altered variant content
- **Comparing large sorted VCF files without loading them into memory**

## Performance

### v1.2 Optimizations

The v1.2 release includes comprehensive performance optimizations:

| Optimization | Description | Impact |
|-------------|-------------|--------|
| **Memory-mapped I/O** | Uses `mmap()` with kernel read-ahead for both VCF files | 50x faster file reads |
| **SIMD line detection** | AVX2/SSE2 vectorized newline scanning (32/16 bytes per cycle) | 10x faster parsing |
| **Zero-copy parsing** | `std::string_view` throughout hot paths, no intermediate allocations | 20x less memory churn |
| **1MB output buffering** | Batched writes reduce syscall overhead | 10x faster output |

### Performance Results

| File Size | Time | Throughput |
|-----------|------|------------|
| 503 MB | 0.73s | ~690 MB/s |
| 6.8 GB × 2 | 53s | ~256 MB/s |

Previously, the tool would timeout on files larger than 1GB.

### Streaming Two-Pointer Merge Diff
When using `--assume-sorted`, the tool implements an efficient streaming comparison:

1. **Two-pointer algorithm**: Reads one line at a time from each file
2. **Constant memory**: Only current line from each file is held in memory
3. **Linear I/O**: Each line is read exactly once
4. **Immediate output**: Differences are reported as they're found

### Performance Characteristics
| Mode | Memory Usage | Requirements |
|------|--------------|--------------|
| Default | O(total_variants) | None - handles unsorted input |
| `--assume-sorted` | O(1) | Input files must be sorted |

### When to Use Streaming Mode
- Input files are already sorted (e.g., output from bcftools, VCFX_sorter)
- Working with very large files (>1GB each)
- Limited available memory
- Comparing output of sorted pipelines

## Output Format

The output consists of two sections:

```
Variants unique to file1.vcf:
chrom:pos:ref:alt
chrom:pos:ref:alt
...

Variants unique to file2.vcf:
chrom:pos:ref:alt
chrom:pos:ref:alt
...
```

Where each variant is represented as a colon-separated string with chromosome, position, reference allele, and sorted alternate alleles.

## Examples

### Basic Usage

```bash
VCFX_diff_tool --file1 original.vcf --file2 modified.vcf
```

### Comparing Variant Caller Outputs

```bash
VCFX_diff_tool --file1 caller1_output.vcf --file2 caller2_output.vcf > caller_differences.txt
```

### Validate Processing Results

```bash
# Check that filtering didn't remove variants it shouldn't have
VCFX_diff_tool --file1 expected_filtered.vcf --file2 actual_filtered.vcf
```

### Streaming Mode for Large Pre-Sorted Files

Compare large pre-sorted files with minimal memory usage:

```bash
VCFX_diff_tool --file1 sorted1.vcf --file2 sorted2.vcf --assume-sorted
```

### Streaming with Natural Chromosome Ordering

Compare pre-sorted files using natural chromosome order:

```bash
VCFX_diff_tool --file1 file1.vcf --file2 file2.vcf --assume-sorted --natural-chr
```

### Pipeline with Sorting

Pre-sort files then compare for optimal performance:

```bash
# Sort individual files first
VCFX_sorter < sample1.vcf > sample1_sorted.vcf
VCFX_sorter < sample2.vcf > sample2_sorted.vcf

# Then compare with streaming mode
VCFX_diff_tool --file1 sample1_sorted.vcf --file2 sample2_sorted.vcf --assume-sorted
```

## Handling Special Cases

- **Multi-allelic variants**: Alternate alleles are sorted alphabetically to ensure consistent comparison even if the order differs between files (e.g., "A,G" and "G,A" are treated as identical)
- **Header differences**: Header lines (starting with #) are ignored, so differences in metadata don't affect the comparison
- **Malformed VCF lines**: Invalid lines are skipped with a warning
- **Empty files**: Properly handled; will show all variants from the non-empty file as unique
- **Missing files**: Reports an error if either file cannot be opened
- **Large files**: Use `--assume-sorted` for files larger than available RAM

## Limitations

- Compares only chromosome, position, reference, and alternate alleles; ignores other fields like quality, filter, and INFO
- Cannot detect differences in sample genotypes
- No support for partial matches or fuzzy comparisons (e.g., variants that differ only in quality)
- Doesn't consider changes in INFO or FORMAT fields as differences
- Cannot compare complex structural variants represented in different ways
- **Streaming mode requires pre-sorted input** - incorrect results if inputs are not sorted

## Backward Compatibility

The tool is fully backward compatible:
- Default behavior (without `--assume-sorted`) works exactly as before
- Existing scripts and pipelines continue to function unchanged
- The `--assume-sorted` option is purely opt-in for performance optimization
