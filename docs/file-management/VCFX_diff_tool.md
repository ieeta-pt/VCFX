# VCFX_diff_tool

## Overview

VCFX_diff_tool compares two VCF files and identifies variants that are unique to each file, providing a simple way to detect differences between variant sets.

## Usage

```bash
VCFX_diff_tool --file1 <file1.vcf> --file2 <file2.vcf>
```

## Options

| Option | Description |
|--------|-------------|
| `-a`, `--file1 <FILE>` | Required. Path to the first VCF file |
| `-b`, `--file2 <FILE>` | Required. Path to the second VCF file |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_diff_tool analyzes two VCF files, compares their variant content, and identifies differences between them. The tool:

1. Loads variants from both VCF files, ignoring header lines
2. Creates a normalized key for each variant based on chromosome, position, reference allele, and sorted alternate alleles
3. Identifies variants that are unique to each file by comparing these keys
4. Reports the differences in a readable format

This tool is particularly useful for:
- Validating VCF file transformations
- Checking tool outputs against expected results
- Comparing variant calls between different callers or pipelines
- Verifying that VCF manipulations haven't inadvertently altered variant content

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
./VCFX_diff_tool --file1 original.vcf --file2 modified.vcf
```

### Comparing Variant Caller Outputs

```bash
./VCFX_diff_tool --file1 caller1_output.vcf --file2 caller2_output.vcf > caller_differences.txt
```

### Validate Processing Results

```bash
# Check that filtering didn't remove variants it shouldn't have
./VCFX_diff_tool --file1 expected_filtered.vcf --file2 actual_filtered.vcf
```

## Handling Special Cases

- **Multi-allelic variants**: Alternate alleles are sorted alphabetically to ensure consistent comparison even if the order differs between files (e.g., "A,G" and "G,A" are treated as identical)
- **Header differences**: Header lines (starting with #) are ignored, so differences in metadata don't affect the comparison
- **Malformed VCF lines**: Invalid lines are skipped with a warning
- **Empty files**: Properly handled; will show all variants from the non-empty file as unique
- **Missing files**: Reports an error if either file cannot be opened
- **Large files**: Efficiently processes files with thousands of variants using hash-based comparison

## Performance

The tool is optimized for efficiency:
- Uses hash sets for O(1) lookups when comparing variants
- Single-pass processing of each input file
- Memory usage scales with the number of unique variants in both files
- Can handle large VCF files with minimal overhead

## Limitations

- Compares only chromosome, position, reference, and alternate alleles; ignores other fields like quality, filter, and INFO
- Cannot detect differences in sample genotypes
- No support for partial matches or fuzzy comparisons (e.g., variants that differ only in quality)
- Not designed to handle VCF files with extremely large numbers of variants (hundreds of millions)
- Doesn't consider changes in INFO or FORMAT fields as differences
- Cannot compare complex structural variants represented in different ways 
