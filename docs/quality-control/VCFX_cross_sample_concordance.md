# VCFX_cross_sample_concordance

## Overview
`VCFX_cross_sample_concordance` analyzes a multi-sample VCF file to determine if genotypes are consistent across all samples for each variant. It identifies variants where samples show concordance (agreement) or discordance (disagreement) in their genotype calls.

## Usage
```bash
VCFX_cross_sample_concordance [OPTIONS] < input.vcf > concordance_results.tsv
```

## Options
| Option | Description |
|--------|-------------|
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |
| `-s`, `--samples` | Comma-separated list of samples to check |

## Description
`VCFX_cross_sample_concordance` examines each variant in a multi-sample VCF file and determines if all samples with valid genotypes have the same normalized genotype. The tool:

1. Normalizes genotypes across samples by:
   - Converting phased genotypes (`|`) to unphased (`/`)
   - Handling multi-allelic variants properly
   - Sorting alleles numerically for consistent comparison (e.g., `1/0` becomes `0/1`)
2. Compares the normalized genotypes across all samples for each variant
3. Classifies each variant as:
   - `CONCORDANT`: All samples with valid genotypes have the same normalized genotype
   - `DISCORDANT`: Samples have different normalized genotypes
   - `NO_GENOTYPES`: No samples have valid genotypes for this variant
4. Outputs detailed per-variant results to standard output
5. Provides a summary of concordance statistics to standard error

This tool is particularly useful for:
- Quality control of multi-sample VCF files
- Identifying potential sample mix-ups or contamination
- Validating genotype calling consistency across technical replicates
- Assessing reliability of variant calls across different sequencing or analysis methods

## Output Format
The tool produces a tab-separated values (TSV) file with the following columns:

| Column | Description |
|--------|-------------|
| CHROM | Chromosome of the variant |
| POS | Position of the variant |
| ID | Variant identifier |
| REF | Reference allele |
| ALT | Alternate allele(s) |
| Num_Samples | Number of samples with valid genotypes |
| Unique_Normalized_Genotypes | Number of distinct normalized genotypes observed |
| Concordance_Status | Status: `CONCORDANT`, `DISCORDANT`, or `NO_GENOTYPES` |

Additionally, a summary of statistics is printed to standard error, including:
- Total number of variants processed
- Number of concordant variants
- Number of discordant variants
- Number of variants with no valid genotypes

## Examples

### Basic Usage
Check concordance in a multi-sample VCF file:
```bash
VCFX_cross_sample_concordance < input.vcf > concordance_results.tsv
```

Specify a subset of samples:

```bash
VCFX_cross_sample_concordance --samples SampleA,SampleB < input.vcf > concordance_subset.tsv
```

### Filtering for Discordant Variants
Identify only variants with discordant genotypes:
```bash
VCFX_cross_sample_concordance < input.vcf | grep "DISCORDANT" > discordant_variants.tsv
```

### Concordance Analysis in a Pipeline
Use as part of a larger analysis pipeline:
```bash
cat input.vcf | VCFX_cross_sample_concordance | awk -F'\t' '{if($8=="CONCORDANT") print $0}' > consistent_variants.tsv
```

### Saving Both Results and Summary
Capture both the detailed results and summary statistics:
```bash
VCFX_cross_sample_concordance < input.vcf > concordance_results.tsv 2> concordance_summary.txt
```

## Genotype Normalization

### Process
For each sample's genotype, the tool performs the following normalization steps:
1. Extracts the genotype field (first field in the FORMAT column)
2. Converts all phase separators (`|`) to unphased separators (`/`)
3. Splits the genotype into individual allele indices
4. Validates each allele (skips if missing or invalid)
5. Sorts allele indices in ascending order (e.g., `1/0` → `0/1`)
6. Rejoins the allele indices with `/` separators

### Example Normalizations
- `0|1` → `0/1`
- `1/0` → `0/1`
- `1|2` → `1/2`
- `2/1` → `1/2`
- `./1` → (skipped as invalid)
- `0/0` → `0/0` (unchanged)

## Handling Special Cases

### Missing Genotypes
- Genotypes with missing values (`./.`, `.`) are excluded from concordance calculation
- If no samples have valid genotypes for a variant, it is classified as `NO_GENOTYPES`

### Multi-allelic Variants
- Alt alleles are parsed from the ALT column (comma-separated)
- Genotype allele indices are validated against the number of alt alleles
- Invalid indices (exceeding the number of alt alleles) are treated as missing

### Phased Genotypes
- Phasing information is ignored for concordance calculation
- Genotypes that differ only in phasing are considered concordant

### Malformed VCF Lines
- Lines with insufficient columns are skipped
- Lines encountered before the #CHROM header are skipped with a warning
- Genotypes that cannot be parsed are treated as missing

## Performance Considerations
- Processes the VCF file line by line, requiring minimal memory
- Time complexity scales linearly with file size and number of samples
- Efficient handling of large multi-sample VCF files with normalized genotype comparisons
- No file indexing or preprocessing required

## Limitations
- Only handles diploid genotypes (e.g., `0/1`, not haploid `0` or polyploid genotypes)
- Ignores genotype phasing information in concordance assessment
- Does not consider genotype quality or other FORMAT fields in the assessment
- No option to adjust concordance thresholds (e.g., requiring 90% sample agreement)
- Cannot output detailed per-sample information for discordant variants 


