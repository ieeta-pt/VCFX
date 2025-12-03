# VCFX_duplicate_remover

## Overview

VCFX_duplicate_remover identifies and removes duplicate variant records from a VCF file based on their essential coordinates (chromosome, position, reference allele, and alternate allele), ensuring that each unique variant is represented only once in the output.

## Usage

```bash
VCFX_duplicate_remover [OPTIONS] < input.vcf > deduplicated.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_duplicate_remover processes a VCF file to detect and remove duplicate variant records. The tool:

1. Reads a VCF file from standard input line by line
2. Passes header lines (beginning with '#') through unchanged
3. For each data line, creates a normalized variant key consisting of:
   - Chromosome name (CHROM)
   - Position (POS)
   - Reference allele (REF)
   - Normalized alternate alleles (ALT)
4. Normalizes multi-allelic variants by sorting the comma-separated ALT alleles alphabetically
5. Tracks seen variants using a hash-based data structure
6. Outputs only the first occurrence of each unique variant, discarding subsequent duplicates
7. Writes the deduplicated VCF to standard output

This tool is particularly useful for:
- Merging VCF files that may contain overlapping variants
- Cleaning up datasets that have inadvertently duplicated records
- Ensuring downstream analysis tools don't process the same variant multiple times

## Output Format

The output is a valid VCF file with the same format as the input, but with duplicate variant records removed. The first occurrence of each unique variant is preserved, maintaining all original fields (ID, QUAL, FILTER, INFO, FORMAT, samples).

## Examples

### Basic Usage

```bash
# Remove duplicate variants from a VCF file
./VCFX_duplicate_remover < input.vcf > deduplicated.vcf
```

### In a Pipeline

```bash
# Filter a VCF file and then remove duplicates
./VCFX_record_filter --quality ">20" < input.vcf | \
./VCFX_duplicate_remover > filtered_unique.vcf
```

### Checking Results

```bash
# Count variants before and after deduplication
echo "Before: $(grep -v '^#' input.vcf | wc -l) variants"
echo "After: $(grep -v '^#' deduplicated.vcf | wc -l) variants"
```

## Duplicate Detection

The tool determines uniqueness based on four key attributes:

1. **Chromosome (CHROM)**: Exact string match of chromosome name
2. **Position (POS)**: Numerical position on the chromosome
3. **Reference Allele (REF)**: The reference sequence at that position
4. **Alternate Alleles (ALT)**: The alternate alleles, normalized by sorting

For multi-allelic variants, the ALT field is normalized by:
- Splitting the comma-separated list of alternate alleles
- Sorting the alleles alphabetically
- Re-joining them with commas

This normalization ensures that variants with the same alleles but in different order are correctly identified as duplicates. For example, "A,G,T" and "T,A,G" are treated as the same variant after normalization.

## Handling Special Cases

The tool implements several strategies for handling edge cases:

1. **Multi-allelic variants**: Normalizes ALT fields by sorting to handle different orderings
2. **Empty lines**: Skipped without affecting output
3. **Malformed lines**: Lines that can't be parsed are skipped with a warning
4. **Position parsing errors**: If a POS field can't be parsed as an integer, it's set to 0
5. **Header lines**: All header lines are preserved in the output
6. **Empty files**: Properly handles empty input files, producing empty output
7. **Files with only headers**: Header lines are passed through correctly
8. **Sample columns**: Maintains all sample genotype data in the output

## Performance

VCFX_duplicate_remover is designed for efficiency:

1. Single-pass processing with O(n) time complexity where n is the number of variants
2. Uses an optimized hash-based data structure for fast variant lookups
3. Minimal memory overhead, proportional to the number of unique variants
4. Handles large VCF files with millions of variants efficiently
5. String processing optimized for performance with minimal copying

## Limitations

1. Retains the first occurrence of duplicate variants; quality or information in subsequent duplicates is discarded
2. No option to select which duplicate to keep (e.g., the one with highest quality)
3. No facility to annotate output variants with duplicate counts
4. Limited to exact matching; doesn't detect overlapping variants that might represent the same event
5. Doesn't consider INFO fields in uniqueness determination
6. Cannot handle duplicates based on sample-specific criteria
7. No option to only report duplicate variants without removing them 
