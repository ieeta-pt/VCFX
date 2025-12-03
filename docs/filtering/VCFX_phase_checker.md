# VCFX_phase_checker

## Overview

The VCFX_phase_checker tool filters a VCF file to retain only those variant lines where all sample genotypes are fully phased (using the pipe '|' phasing separator). This is particularly useful for downstream analyses that require complete phasing information.

## Usage

```bash
VCFX_phase_checker [OPTIONS] < input.vcf > phased_output.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_phase_checker reads a VCF file from standard input and examines the GT (genotype) field for every sample in each variant line. It determines whether each genotype is fully phased using the following criteria:

- A genotype is considered fully phased only if:
  - It uses the pipe character '|' as the separator between alleles (e.g., "0|1")
  - It contains no missing alleles (no ".")
  - It contains no unphased separators ("/")

The tool outputs only those variant lines where all sample genotypes meet these criteria. Lines that don't meet these criteria are skipped, and warnings are written to standard error.

## Output

The output is a valid VCF file containing:
- All header lines from the input file (unchanged)
- Only the variant lines where all samples have fully phased genotypes
- Warnings (to stderr) for each line that was skipped due to unphased genotypes

## Examples

### Basic Usage

```bash
./VCFX_phase_checker < input.vcf > phased_output.vcf
```

### Capturing Warnings

```bash
./VCFX_phase_checker < input.vcf > phased_output.vcf 2> unphased_warnings.log
```

### Counting Phased vs. Unphased Variants

```bash
# Count total variants
total=$(grep -v "^#" input.vcf | wc -l)
# Count phased variants
phased=$(./VCFX_phase_checker < input.vcf | grep -v "^#" | wc -l)
# Calculate percentage
echo "Phased variants: $phased / $total ($(echo "scale=2; 100*$phased/$total" | bc)%)"
```

## Handling Special Cases

- **Haploid genotypes** (e.g., "0" or "1"): These are not considered phased; the line will be skipped
- **Missing genotypes** (e.g., "./." or ".|."): These are not considered phased; the line will be skipped
- **Missing GT field**: Lines without a GT field in the FORMAT column are skipped with a warning
- **Multiallelic variants**: These are treated the same as biallelic variants, as long as all alleles are phased
- **Non-VCF-compliant genotype notation**: Any genotype that doesn't follow standard VCF format is not considered phased
- **Header lines**: All header lines (starting with "#") are preserved in the output
- **Samples with different ploidy levels**: Each sample is checked independently; if all are phased, the line is kept

## Performance

The tool processes files line by line with minimal memory requirements, allowing it to handle very large VCF files efficiently.

## Limitations

- No option to make a best-effort phasing assumption
- Cannot output partially phased lines or filter specific samples
- Designed to be used as part of a pipeline, not as a standalone phasing tool 
