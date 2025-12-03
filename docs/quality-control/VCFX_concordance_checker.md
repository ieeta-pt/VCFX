# VCFX_concordance_checker

## Overview
`VCFX_concordance_checker` compares genotypes between two specified samples within a VCF file to determine concordance (agreement) or discordance (disagreement) for each variant. This tool is useful for comparing genotype calls between different samples, such as technical replicates or related individuals.

## Usage
```bash
VCFX_concordance_checker --samples "SAMPLE1 SAMPLE2" < input.vcf > concordance_report.tsv
```

## Options
| Option | Description |
|--------|-------------|
| `-s`, `--samples "SAMPLE1 SAMPLE2"` | Required. Names of the two samples to compare, separated by a space |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description
`VCFX_concordance_checker` analyzes a VCF file and compares the genotypes of two specified samples for each variant. The tool:

1. Normalizes genotypes by:
   - Converting phased genotypes (`|`) to unphased (`/`)
   - Sorting alleles numerically (e.g., `1/0` becomes `0/1`)
   - Validating against available alternate alleles
2. Compares the normalized genotypes between the two samples
3. Classifies each variant as:
   - `Concordant`: Both samples have identical normalized genotypes
   - `Discordant`: Samples have different normalized genotypes
4. Outputs detailed per-variant results to standard output
5. Provides a summary of concordance statistics to standard error

This tool is particularly useful for:
- Quality control of technical replicates
- Comparing genotype calls between related samples
- Validating sample identity
- Assessing reproducibility of variant calling pipelines

## Output Format
The tool produces a tab-separated values (TSV) file with the following columns:

| Column | Description |
|--------|-------------|
| CHROM | Chromosome of the variant |
| POS | Position of the variant |
| ID | Variant identifier |
| REF | Reference allele |
| ALT | Alternate allele(s) |
| SAMPLE1_GT | Normalized genotype of the first sample |
| SAMPLE2_GT | Normalized genotype of the second sample |
| Concordance | Status: `Concordant` or `Discordant` |

Additionally, a summary of statistics is printed to standard error, including:
- Total number of variants compared
- Number of concordant genotypes
- Number of discordant genotypes

## Examples

### Basic Usage
Check concordance between two samples in a VCF file:
```bash
VCFX_concordance_checker --samples "SAMPLE1 SAMPLE2" < input.vcf > concordance_report.tsv
```

### Filtering for Discordant Variants
Identify only variants with discordant genotypes:
```bash
VCFX_concordance_checker --samples "SAMPLE1 SAMPLE2" < input.vcf | grep "Discordant" > discordant_variants.tsv
```

### Calculating Concordance Rate
Count concordant variants and calculate rate:
```bash
VCFX_concordance_checker --samples "SAMPLE1 SAMPLE2" < input.vcf > report.tsv 2> stats.txt
grep -c "Concordant" report.tsv > concordant_count.txt
```

### Integration with Other Tools
Use as part of a larger analysis pipeline:
```bash
cat input.vcf | VCFX_record_filter --filter "QUAL>30" | VCFX_concordance_checker --samples "SAMPLE1 SAMPLE2" > high_quality_concordance.tsv
```

## Genotype Normalization

### Process
For each sample's genotype, the tool performs the following normalization steps:
1. Extracts the genotype field (GT)
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
- Genotypes with missing values (`./.`, `.`) are excluded from comparison
- Variants where either sample has missing genotypes are skipped

### Multi-allelic Variants
- Alt alleles are parsed from the ALT column (comma-separated)
- Genotype allele indices are validated against the number of alt alleles
- For multi-allelic variants, the tool correctly compares numerically sorted genotypes

### Phased Genotypes
- Phasing information is ignored for concordance calculation
- Genotypes that differ only in phasing are considered concordant (e.g., `0|1` and `0/1`)

### Invalid Genotypes
- Genotypes with non-numeric allele indices are skipped
- Allele indices that exceed the number of alternate alleles are treated as invalid

### Malformed VCF Lines
- Lines with insufficient columns are skipped
- Lines encountered before the #CHROM header cause an error
- VCF files without both specified samples cause an error and program termination

## Performance Considerations
- Processes the VCF file line by line, requiring minimal memory
- No preprocessing or indexing of the VCF file is required
- Linear time complexity with respect to file size
- Provides a quick summary of concordance statistics for rapid quality assessment

## Limitations
- Limited to exactly two samples for comparison
- Only processes diploid genotypes (e.g., `0/1`, not haploid `0` or polyploid genotypes)
- Ignores all FORMAT fields except for GT (genotype)
- No consideration of genotype quality or other metrics in concordance assessment
- Cannot generate detailed statistics on concordance by variant type
- Does not allow customizing concordance criteria beyond exact genotype match 
