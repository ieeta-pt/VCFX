# VCFX_ancestry_inferrer

## Overview

VCFX_ancestry_inferrer infers the likely population ancestry for each sample in a VCF file by comparing sample genotypes to known population allele frequencies.

## Usage

```bash
VCFX_ancestry_inferrer --frequency <freq_file> [OPTIONS] < input.vcf > ancestry_results.txt
```

## Options

| Option | Description |
|--------|-------------|
| `--frequency <FILE>` | Required. Path to a file containing population-specific allele frequencies |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_ancestry_inferrer analyzes the genotypes of samples in a VCF file and compares them to known population-specific allele frequencies to determine the most likely ancestry for each sample. The tool:

1. Reads a frequency reference file containing population-specific allele frequencies
2. Processes the VCF file, examining each biallelic or multiallelic variant
3. For each sample, calculates ancestry scores by comparing observed genotypes to population frequency data
4. Assigns each sample to the population with the highest cumulative score
5. Outputs a simple table mapping each sample to its inferred population

The ancestry inference is based on the principle that individuals from a specific population are more likely to carry alleles at frequencies matching that population's known frequency distribution.

## Output Format

The output is a tab-delimited text file with the following columns:

```
Sample  Inferred_Population
```

Where:
- Sample is the sample name from the VCF file
- Inferred_Population is the population with the highest ancestry score

## Examples

### Basic Usage

```bash
./VCFX_ancestry_inferrer --frequency population_freqs.txt < samples.vcf > ancestry_results.txt
```

### Creating a Frequency Reference File

The frequency file should have the following tab-delimited format:
```
CHROM  POS  REF  ALT  POPULATION  FREQUENCY
1      100  A    G    EUR         0.75
1      100  A    G    AFR         0.10
1      100  A    G    EAS         0.25
```

### Using with Multi-Population Data

```bash
# Combine ancestry results with other data
./VCFX_ancestry_inferrer --frequency global_freqs.txt < diverse_cohort.vcf | \
  join -t $'\t' -1 1 -2 1 - phenotype_data.txt > annotated_results.tsv
```

## Algorithm

The ancestry inference algorithm works as follows:

1. For each variant in the VCF file:
   - For each sample with a non-reference genotype:
     - Look up the frequency of that allele in each reference population
     - Add the frequency value to that population's score for the sample
     
2. After processing all variants:
   - For each sample, find the population with the highest cumulative score
   - Assign the sample to that population

This approach assigns more weight to alleles that are common in a specific population but rare in others, making them more informative for ancestry inference.

## Handling Special Cases

- **Multi-allelic variants**: Each alternate allele is treated separately and looked up in the frequency reference
- **Phased genotypes**: Phase information is ignored; both "0|1" and "0/1" are treated identically
- **Missing genotypes**: Missing genotypes ("./.") are skipped and don't contribute to ancestry scores
- **Missing frequency data**: Variants without corresponding frequency data are skipped
- **Identical scores**: If two populations have identical scores for a sample, the first one alphabetically is assigned
- **Diploid genotypes**: Both alleles contribute independently to the ancestry score
- **Empty VCF**: Will produce no output rows (empty output file)
- **Unknown populations**: Only populations defined in the frequency file will be considered

## Performance

The tool is optimized for efficiency:
- Uses hash maps for constant-time lookups of frequency data
- Single-pass processing of the VCF file
- Memory usage scales with:
  - The number of variants in the frequency file
  - The number of reference populations
  - The number of samples in the VCF

## Limitations

- Accuracy depends on the quality and relevance of the population frequency data
- Works best with large numbers of variants (hundreds to thousands)
- Not designed for detecting admixed individuals (reports only the highest-scoring population)
- Assumes independence between variants (does not account for linkage disequilibrium)
- No confidence scores or statistical measures of assignment certainty
- Cannot handle non-biallelic complex variants (e.g., structural variants)
- Doesn't account for sample relatedness within the input VCF 
