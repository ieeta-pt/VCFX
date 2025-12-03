# VCFX_allele_balance_filter

## Overview

VCFX_allele_balance_filter filters a VCF file to keep only variants where all samples have an allele balance ratio (reference alleles / total alleles) above a specified threshold, allowing for quality control and bias detection in variant calls.

## Usage

```bash
VCFX_allele_balance_filter --filter-allele-balance <THRESHOLD> < input.vcf > filtered.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `-f`, `--filter-allele-balance` <VAL> | Required. Allele balance threshold between 0.0 and 1.0 |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_allele_balance_filter examines the genotypes in each variant record and calculates an allele balance ratio for each sample. The tool:

1. Processes a VCF file line by line
2. For each variant, calculates the allele balance for each sample's genotype
3. Filters out variants where any sample has an allele balance below the specified threshold
4. Passes through all header lines unchanged
5. Outputs a filtered VCF file with only the passing variants

Allele balance is calculated as:
```
AB = (number of reference alleles) / (total number of alleles)
```

Where:
- Reference alleles are those with value "0" in the genotype field
- All other numeric alleles (1, 2, 3, etc.) are counted as alternate alleles
- Missing or non-numeric alleles are excluded from the calculation

This tool is useful for:
- Detecting potential sequencing or mapping biases
- Quality control of variant calls
- Filtering out variants with skewed allele representation
- Identifying potential sample contamination or mixed samples

## Output Format

The output is a standard VCF file with the same format as the input, but containing only the variant lines that pass the allele balance filter. All header lines are preserved.

## Examples

### Basic Usage

```bash
# Keep variants where all samples have allele balance >= 0.3
VCFX_allele_balance_filter --filter-allele-balance 0.3 < input.vcf > balanced.vcf
```

### Stringent Filtering

```bash
# Very stringent filtering (close to 50/50 balance required)
VCFX_allele_balance_filter --filter-allele-balance 0.45 < input.vcf > highly_balanced.vcf
```

### Counting Filtered Variants

```bash
# Count how many variants were filtered out
input_count=$(grep -v "^#" input.vcf | wc -l)
output_count=$(grep -v "^#" filtered.vcf | wc -l)
filtered_count=$((input_count - output_count))
echo "Filtered out $filtered_count variants based on allele balance"
```

### In a Pipeline

```bash
# Filter by quality then by allele balance
grep -v "^#" input.vcf | grep "PASS" | grep "QUAL>30" | \
VCFX_allele_balance_filter --filter-allele-balance 0.4 > high_quality_balanced.vcf
```

## Genotype Interpretation

The tool examines the GT field of each sample's genotype:

1. Extracts the GT field (before the first colon if present)
2. Treats both phased ('|') and unphased ('/') genotypes the same
3. For each allele:
   - '0' is counted as a reference allele
   - Any other number (1, 2, 3, etc.) is counted as an alternate allele
   - Non-numeric values are ignored

For example:
- "0/0" → AB = 2/2 = 1.0 (all reference)
- "0/1" → AB = 1/2 = 0.5 (half reference, half alternate)
- "1/1" → AB = 0/2 = 0.0 (all alternate)
- "0/2" → AB = 1/2 = 0.5 (half reference, half alternate)
- "./1" → AB = 0/1 = 0.0 (missing reference allele, only alternate counted)
- "./." → AB = 0/0 = 0.0 (no valid alleles)

## Handling Special Cases

- **Missing genotypes** ("./.") are treated as having AB = 0.0
- **Partial missing** ("./1") counts only the present alleles
- **Non-diploid genotypes** (e.g., "0/1/2") are handled correctly by counting alleles individually
- **Complex genotypes** (non-numeric) are skipped when calculating AB
- **Empty lines** are ignored
- **Header lines** are preserved unchanged
- **Malformed VCF lines** with insufficient columns are skipped with a warning
- **Multi-allelic variants** have all non-reference alleles (1, 2, 3, etc.) treated as alternate

## Performance

The tool is optimized for efficiency:
- Processes the VCF file in a single pass
- Minimal memory usage as it processes one variant at a time
- Constant-time computation of allele balance
- Stops calculating AB for a variant as soon as any sample fails the threshold

## Limitations

1. Uses a simple all-or-nothing approach (variant passes only if ALL samples pass)
2. No option to specify which samples to include in the filtering
3. Cannot handle sample-specific threshold values
4. No detailed reporting on which samples/variants failed and by how much
5. No option to annotate variants with their allele balance rather than filtering
6. Limited to the strict definition of allele balance (ref/total), not accounting for strand bias
7. Treats all alternate alleles equally, regardless of their identity 
