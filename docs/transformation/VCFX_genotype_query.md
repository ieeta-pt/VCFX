# VCFX_genotype_query

## Overview

VCFX_genotype_query filters VCF files to keep only those variant lines where at least one sample has a specified genotype. It's a powerful tool for extracting variants with specific genotypic patterns.

## Usage

```bash
VCFX_genotype_query [OPTIONS] < input.vcf > filtered.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `--genotype-query`, `-g` "GENOTYPE" | Specify the genotype to query (e.g., "0/1", "1/1") |
| `--strict` | Use strict string comparison (no phasing unification or allele sorting) |
| `--help`, `-h` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_genotype_query reads a VCF file from standard input and outputs variants (plus all header lines) where at least one sample has a genotype matching the specified pattern in the 'GT' field.

By default, the tool uses flexible matching:
- Unifies phasing separators ('|' and '/' are treated the same)
- Sorts alleles numerically (e.g., "1/0" is treated the same as "0/1")

When the `--strict` option is used, the tool performs exact string matching on genotypes, maintaining phase status and allele order.

## Output Format

The output is a standard VCF file containing:
- All header lines from the input file (unchanged)
- Only variant lines where at least one sample matches the specified genotype pattern

## Examples

### Basic Usage (Flexible Matching)

```bash
# Find all variants where any sample is heterozygous (0/1)
./VCFX_genotype_query --genotype-query "0/1" < input.vcf > heterozygous.vcf
```

### Strict Matching

```bash
# Find variants where any sample has a specifically phased heterozygous genotype (0|1)
./VCFX_genotype_query --genotype-query "0|1" --strict < input.vcf > phased_heterozygous.vcf
```

### Finding Homozygous Variants

```bash
# Find variants where any sample is homozygous for the alternate allele
./VCFX_genotype_query --genotype-query "1/1" < input.vcf > homozygous_alt.vcf
```

### Finding Multi-allelic Genotypes

```bash
# Find variants where any sample has specific multi-allelic genotype
./VCFX_genotype_query --genotype-query "1/2" < input.vcf > multi_allelic.vcf
```

## Handling Special Cases

- **Phased genotypes**: By default, "0|1" and "0/1" are considered equivalent; use `--strict` to differentiate them
- **Multi-allelic variants**: Can be queried with appropriate genotypes (e.g., "1/2", "0/2")
- **Missing genotypes**: Can match against "././" or "."
- **Non-diploid genotypes**: Supported in both flexible and strict modes
- **Malformed VCF lines**: Lines with fewer than 10 columns are skipped with a warning
- **Missing GT field**: Lines without a GT field in the FORMAT column are skipped

## Performance

The tool processes VCF files line by line with minimal memory requirements, making it efficient for large files. Performance scales linearly with the number of samples in the VCF file.

## Limitations

- Only filters based on the presence of the target genotype in at least one sample
- Cannot filter based on multiple genotype patterns in a single run
- No support for complex queries (e.g., specific samples with specific genotypes)
- Skips variants where the GT field is not present in the FORMAT column 
