# VCFX_af_subsetter

## Overview
`VCFX_af_subsetter` filters variants in a VCF file based on allele frequency (AF) values, allowing selection of variants within a specified frequency range. This tool is useful for focusing analysis on variants of specific population prevalence.

## Usage
```bash
# Standard mode (stdin)
VCFX_af_subsetter --af-filter "MIN-MAX" < input.vcf > filtered.vcf

# File mode (optimized for large files)
VCFX_af_subsetter -i input.vcf --af-filter "MIN-MAX" > filtered.vcf

# Positional argument mode
VCFX_af_subsetter input.vcf --af-filter "MIN-MAX" > filtered.vcf
```

## Options
| Option | Description |
|--------|-------------|
| `-a`, `--af-filter <MIN-MAX>` | Allele frequency range for filtering (e.g., `0.01-0.05`) |
| `-i`, `--input FILE` | Input VCF file (uses memory-mapped I/O for best performance) |
| `-q`, `--quiet` | Suppress warning messages |
| `-h`, `--help` | Display help message and exit |
| `-v`, `--version` | Show program version and exit |

## Description
`VCFX_af_subsetter` processes VCF files line by line and filters variants based on their allele frequency (AF) values from the INFO field. The tool:

1. Reads the VCF file from standard input
2. Parses the AF values from the INFO field of each variant record
3. Compares these values against the specified minimum and maximum thresholds
4. Retains variants with at least one allele frequency value within the range [MIN, MAX]
5. Outputs the filtered variants to standard output

This tool is particularly useful for:
- Isolating rare variants (e.g., AF < 0.01)
- Focusing on common variants (e.g., AF > 0.05)
- Selecting variants with specific population frequencies
- Removing extremely rare or fixed variants from analysis

## Output Format
The output is a standard VCF file containing:
- All original header lines from the input VCF
- Only those variant records with AF values within the specified range
- No modification to the content or format of the retained lines

## Examples

### Basic Usage
Filter for rare variants with frequency between 1% and 5%:
```bash
VCFX_af_subsetter --af-filter "0.01-0.05" < input.vcf > rare_variants.vcf
```

### Common Variants
Filter for common variants with frequency above 5%:
```bash
VCFX_af_subsetter --af-filter "0.05-1.0" < input.vcf > common_variants.vcf
```

### Extremely Rare Variants
Filter for extremely rare variants:
```bash
VCFX_af_subsetter --af-filter "0.0001-0.001" < input.vcf > very_rare_variants.vcf
```

### Specific Frequency Band
Filter for variants with a specific frequency band:
```bash
VCFX_af_subsetter --af-filter "0.4-0.6" < input.vcf > mid_frequency_variants.vcf
```

### In Pipeline
Use in a pipeline with other VCFX tools:
```bash
cat input.vcf | VCFX_af_subsetter --af-filter "0.01-0.05" | VCFX_phred_filter -p 30 > high_quality_rare_variants.vcf
```

## AF Value Parsing

### Format Requirements
The tool expects AF values in the INFO field in standard VCF format:
- As a key-value pair in the INFO column: `AF=0.123`
- For multi-allelic sites, as comma-separated values: `AF=0.01,0.05,0.1`

### Range Specification
The AF range must be specified as:
- Two numeric values between 0.0 and 1.0
- Connected by a hyphen (`-`)
- With the first value (minimum) less than or equal to the second value (maximum)

For example: `0.01-0.05`, `0.0-0.1`, `0.4-0.6`

## Handling Special Cases

### Multi-allelic Variants
For variants with multiple alternate alleles (multi-allelic):
- The INFO field may contain multiple AF values (comma-separated)
- The variant is retained if ANY of the AF values fall within the specified range
- This behavior allows for selective filtering of multi-allelic sites

### Missing AF Values
Variants without an AF annotation in the INFO field:
- Are skipped with a warning message
- Are not included in the output
- Can indicate variants where frequency information is unavailable

### Full Range
Using the range `0.0-1.0` will:
- Keep all variants with valid AF values
- Still skip variants lacking AF annotations
- Effectively function as a filter for "has valid AF information"

### Malformed Values
The tool handles several edge cases:
- Invalid range format: Reports an error if not in `MIN-MAX` format
- Out-of-range values: Ensures MIN and MAX are between 0.0 and 1.0
- Inverted ranges: Reports an error if MIN > MAX
- Non-numeric AF values: Skips variants where AF cannot be parsed as a number

## Performance

The tool offers two execution modes with different performance characteristics:

### File Mode (-i/--input)
When using the `-i` option, the tool uses memory-mapped I/O (mmap) for optimal performance on large files. This mode provides approximately **20x speedup** compared to stdin mode.

| File Size | Variants | Stdin Mode | File Mode | Speedup |
|-----------|----------|------------|-----------|---------|
| 4 GB      | 427K     | ~2m 44s    | ~8s       | ~20x    |

### Stdin Mode (default)
For smaller files or pipeline usage, the stdin mode processes input line by line with minimal memory requirements.

### Recommendations
- For files > 100MB, use `-i input.vcf` for best performance
- For pipelines or compressed input, use stdin mode with `zcat file.vcf.gz | VCFX_af_subsetter`

## Performance Characteristics
- Linear time complexity with respect to input file size
- Memory-mapped I/O for efficient large file processing
- 1MB output buffer for optimized write performance
- No preprocessing or indexing required

## Limitations
- Requires AF field to be present in the INFO column
- No way to customize the AF field name (hardcoded to "AF")
- Cannot filter based on other frequency metrics (e.g., MAF, AC, AN)
- No option to include variants with missing AF values
- No statistics provided on the number of variants filtered
- Cannot combine with other filtering criteria in a single command 
