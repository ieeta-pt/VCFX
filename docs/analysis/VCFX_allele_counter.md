# VCFX_allele_counter

## Overview
`VCFX_allele_counter` counts the number of reference and alternate alleles in each sample for each variant in a VCF file. This tool provides a simple way to quantify allele occurrences across samples.

## Usage
```bash
VCFX_allele_counter [OPTIONS] < input.vcf > allele_counts.tsv
```

## Options
| Option | Description |
|--------|-------------|
| `-s`, `--samples "Sample1 Sample2..."` | Optional. Specify sample names to calculate allele counts for (space-separated). If omitted, all samples are processed. |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description
`VCFX_allele_counter` processes a VCF file and counts reference and alternate alleles for each variant in each specified sample. The tool:

1. Reads a VCF file from standard input
2. Identifies sample columns from the VCF header
3. For each variant and each sample:
   - Extracts the genotype information
   - Counts reference alleles (0) and alternate alleles (non-0)
   - Outputs both counts in a tabular format
4. Outputs a tab-separated file with allele counts for each variant-sample combination

This tool is particularly useful for:
- Analyzing allele distribution across samples
- Quantifying the presence of specific alleles
- Preparing data for population genetics analyses
- Validating genotype calls across samples

## Output Format
The tool produces a tab-separated values (TSV) file with the following columns:

| Column | Description |
|--------|-------------|
| CHROM | Chromosome of the variant |
| POS | Position of the variant |
| ID | Variant identifier |
| REF | Reference allele |
| ALT | Alternate allele(s) |
| Sample | Sample name |
| Ref_Count | Number of reference alleles (0) in the sample's genotype |
| Alt_Count | Number of alternate alleles (non-0) in the sample's genotype |

## Examples

### Basic Usage (All Samples)
Count alleles for all samples in a VCF file:
```bash
VCFX_allele_counter < input.vcf > allele_counts_all.tsv
```

### Specific Samples
Count alleles for specific samples:
```bash
VCFX_allele_counter --samples "SAMPLE1 SAMPLE2" < input.vcf > allele_counts_subset.tsv
```

### Using with Other Tools
Process the output for further analysis:
```bash
VCFX_allele_counter < input.vcf | awk -F'\t' '$8 > 0' > samples_with_alt_alleles.tsv
```

## Allele Counting Method

### Reference Alleles
The tool counts an allele as a reference allele when it has the value "0" in the genotype field. For example:
- In genotype "0/0", there are 2 reference alleles
- In genotype "0/1", there is 1 reference allele
- In genotype "1/2", there are 0 reference alleles

### Alternate Alleles
The tool counts an allele as an alternate allele when it has any non-zero numeric value in the genotype field. For example:
- In genotype "0/0", there are 0 alternate alleles
- In genotype "0/1", there is 1 alternate allele
- In genotype "1/2", there are 2 alternate alleles
- In genotype "1/1", there are 2 alternate alleles

### Handling Special Cases
- Missing genotypes (e.g., "./.", ".|."): No counts are recorded for these samples
- Partial missing (e.g., "0/."): Only the valid allele is counted
- Non-numeric alleles: These are skipped and not counted

## Handling Special Cases

### Missing Data
- Genotypes with missing values (`./.`, `.`) are skipped
- Partial missing genotypes only count the valid alleles present

### Multi-allelic Sites
- All non-reference alleles are counted as "alternate" regardless of their specific number
- For example, in a genotype "1/2", both alleles count as alternate alleles
- The tool does not differentiate between different alternate alleles

### Phased Genotypes
- Phasing information is ignored for allele counting
- Phased genotypes (e.g., "0|1") are treated the same as unphased (e.g., "0/1")

### Invalid Genotypes
- Non-numeric allele values are skipped
- Empty genotype fields are skipped

## Performance Considerations
- Processes VCF files line by line, with minimal memory requirements
- Scales linearly with input file size and number of samples
- For very large VCF files with many samples, specifying a subset of samples can improve performance

## Limitations
- Does not distinguish between different alternate alleles (e.g., "1" vs "2")
- No options for filtering by allele count thresholds
- Cannot account for genotype quality or read depth
- Limited to processing standard VCF genotype fields
- Does not produce summary statistics or aggregate counts
- No direct integration with population genetics metrics 
