# VCFX_missing_data_handler

## Overview

VCFX_missing_data_handler identifies and processes missing genotype data in VCF files. It can either flag missing genotypes (default behavior) or impute them with a specified default value, ensuring consistent data representation for downstream analysis.

## Usage

```bash
VCFX_missing_data_handler [OPTIONS] [files...] > processed.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `--fill-missing`, `-f` | Impute missing genotypes with a default value |
| `--default-genotype`, `-d` <GEN> | Specify the default genotype for imputation (default: "./.")  |
| `--threads`, `-t` <NUM> | Number of threads for parallel processing (default: auto-detect CPU cores) |
| `--help`, `-h` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_missing_data_handler processes VCF files to identify and handle missing genotype data. The tool:

1. Reads one or more VCF files (or standard input if no files specified)
2. Identifies missing genotypes in each variant (empty, ".", "./.", or ".|.")
3. Either:
   - Leaves missing genotypes unchanged (default behavior)
   - Replaces missing genotypes with a user-specified value
4. Outputs the processed VCF data to standard output

This tool is particularly useful for:
- Preparing VCF files for tools that don't handle missing genotypes well
- Standardizing the representation of missing data
- Imputing missing genotypes with reference (0/0) or other default values
- Processing multiple VCF files with consistent handling of missing data

## Output Format

The output is a valid VCF file with the same format as the input, but with missing genotypes either left as-is or replaced with the specified default value. All header lines are preserved.

## Examples

### Basic Usage (Flag Only)

```bash
# Process a single file, leaving missing genotypes as-is
VCFX_missing_data_handler < input.vcf > flagged_output.vcf
```

### Impute Missing Data with Default Value

```bash
# Replace missing genotypes with the default value (./.):
VCFX_missing_data_handler --fill-missing < input.vcf > imputed_output.vcf
```

### Impute with Custom Genotype

```bash
# Replace missing genotypes with homozygous reference (0/0):
VCFX_missing_data_handler --fill-missing --default-genotype "0/0" < input.vcf > ref_imputed.vcf
```

### Process Multiple Files

```bash
# Process multiple files at once:
VCFX_missing_data_handler --fill-missing file1.vcf file2.vcf > combined_output.vcf
```

### Multi-threaded Processing

```bash
# Use 8 threads for processing large files:
VCFX_missing_data_handler --fill-missing --threads 8 < large_file.vcf > output.vcf

# Use single-threaded mode (disables parallelism):
VCFX_missing_data_handler --fill-missing --threads 1 < input.vcf > output.vcf
```

### In a Pipeline

```bash
# Filter a VCF file and then handle missing data:
grep -v "^#" input.vcf | grep "PASS" | \
VCFX_missing_data_handler --fill-missing --default-genotype "0/0" > filtered_imputed.vcf
```

## Missing Genotype Detection

The tool identifies the following representations of missing data:

1. Empty genotype field
2. Single dot: "."
3. Pair of dots with slash: "./."
4. Pair of dots with pipe: ".|."

## Handling Special Cases

- **No GT field in FORMAT**: If the FORMAT column does not include a GT field, the variant line is left unchanged
- **Invalid variant lines**: Lines with fewer than 9 columns are passed through unchanged
- **Multiple input files**: Processes each file in sequence, properly handling headers
- **Sample columns structure**: Carefully preserves the structure of sample columns, only modifying the GT field
- **Empty lines**: Preserved with a single newline
- **Header lines**: Passed through unchanged
- **Data before header**: Able to handle invalid VCF files where data appears before the header (with a warning)

## Performance

The tool uses advanced optimization techniques for maximum throughput:

1. **Memory-mapped I/O**: Uses `mmap` with `MADV_SEQUENTIAL` hints for optimal file reading
2. **Multi-threading**: Parallel processing across all available CPU cores (configurable with `--threads`)
3. **Zero-copy pass-through**: Lines without missing genotypes are written directly without parsing
4. **Fast pattern scanning**: Uses `memchr` for efficient detection of potential missing genotype markers
5. **Chunk-based processing**: Data is divided into chunks for parallel processing with proper line boundary handling

### Benchmark Results

On a test file with 427K variants and 2,504 samples:
- **Original implementation**: ~455 seconds
- **Optimized implementation**: ~9 seconds
- **Improvement**: ~50x faster

## Limitations

1. No option to specify which samples should have their missing data imputed
2. Cannot handle phased vs. unphased genotype distinction in imputation
3. No support for probabilistic imputation based on population frequencies
4. No ability to flag sites with a high proportion of missing data
5. Cannot process only specific regions of a VCF file
6. Imputes with the same value regardless of context or neighboring genotypes
7. No reporting of the number or percentage of imputed genotypes 
