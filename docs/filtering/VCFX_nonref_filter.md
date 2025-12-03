# VCFX_nonref_filter

## Overview

VCFX_nonref_filter removes variants from a VCF file where all samples are homozygous reference (0/0), retaining only variants where at least one sample has an alternate allele or a missing genotype.

## Usage

```bash
VCFX_nonref_filter [OPTIONS] < input.vcf > filtered.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_nonref_filter examines each variant in a VCF file and filters out those where all samples are homozygous reference (0/0). The tool:

1. Processes a VCF file line by line
2. For each variant, examines the genotype (GT) field of all samples
3. Determines if every sample is definitively homozygous reference
4. Retains variants where at least one sample has a non-reference allele or missing data
5. Outputs a filtered VCF with only the retained variants
6. Passes through all header lines unchanged

This tool is particularly useful for:
- Removing uninformative variants where no sample has an alternate allele
- Reducing VCF file size by filtering out invariant sites
- Focusing analysis on polymorphic sites
- Preparing variant files for downstream analysis tools that expect polymorphic sites

## Output Format

The output is a standard VCF file with the same format as the input, but containing only variants where at least one sample has a non-reference allele or a missing genotype. All header lines are preserved.

## Examples

### Basic Usage

```bash
# Filter out variants where all samples are homozygous reference
VCFX_nonref_filter < input.vcf > filtered.vcf
```

### Counting Filtered Variants

```bash
# Count how many variants were removed/retained
input_count=$(grep -v "^#" input.vcf | wc -l)
output_count=$(grep -v "^#" filtered.vcf | wc -l)
removed_count=$((input_count - output_count))
echo "Removed $removed_count homozygous reference variants out of $input_count total variants"
```

### In a Pipeline

```bash
# Filter VCF file first by quality, then by non-reference status
grep -v "^#" input.vcf | grep "PASS" | \
VCFX_nonref_filter > high_quality_nonref.vcf
```

### Combining with Other Filters

```bash
# Create a pipeline of filters
cat input.vcf | \
VCFX_nonref_filter | \
VCFX_phred_filter --phred-filter 30 > filtered.vcf
```

## Homozygous Reference Detection

The tool uses comprehensive logic to identify homozygous reference genotypes:

1. If a genotype is missing (e.g., "./.", ".|.", or "."), it's considered NOT homozygous reference, and the variant is retained
2. For each specified allele in a genotype:
   - The allele must be "0" for it to be considered reference
   - Any non-"0" allele (including "1", "2", etc.) is considered alternate

For example:
- "0/0" → Homozygous reference (filtered if all samples have this)
- "0/1" → Heterozygous (variant retained)
- "1/1" → Homozygous alternate (variant retained)
- "./." → Missing genotype (variant retained)
- "0/." → Partially missing (variant retained)
- "0/0/0" → Polyploid homozygous reference (filtered if all samples have reference)

## Handling Special Cases

- **Missing genotypes**: Variants with samples having missing genotypes ("./.") are retained
- **Partial missing**: Genotypes with some missing alleles (e.g., "0/.") are considered not definitively homozygous reference, so the variant is retained
- **No GT field**: If the GT field is not present in the FORMAT column, the variant is retained
- **Empty lines**: Skipped in output
- **Header lines**: Preserved unchanged
- **Malformed VCF lines**: Lines with fewer than 10 columns (required for at least one sample) are passed through unchanged
- **Data before header**: Warning issued and lines passed through unchanged

## Performance

The tool is designed for efficiency:

1. Processes the VCF file line by line, with minimal memory requirements
2. Fast determination of homozygous reference status using string parsing
3. Early exit when a non-homozygous sample is found
4. No requirements to load the entire file into memory

## Limitations

1. No option to filter based on a subset of samples
2. Cannot retain specific homozygous reference variants based on other criteria
3. No support for filtering by percentage of non-reference samples
4. Missing data is always treated as "not definitely homozygous reference"
5. No built-in option to keep variants where the reference allele might be incorrect
6. Cannot incorporate quality values in the filtering decision
7. No reporting on the number of variants removed or statistics about filtered variants 
