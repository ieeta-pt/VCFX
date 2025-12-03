# VCFX_ancestry_assigner

## Overview

VCFX_ancestry_assigner assigns samples in a VCF file to ancestral populations using a likelihood-based approach based on population-specific allele frequencies.

## Usage

```bash
VCFX_ancestry_assigner --assign-ancestry <freq_file> < input.vcf > ancestry_results.txt
```

## Options

| Option | Description |
|--------|-------------|
| `-a`, `--assign-ancestry <FILE>` | Required. Path to a file containing population-specific allele frequencies |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_ancestry_assigner determines the most likely ancestral population for each sample in a VCF file by calculating genotype likelihoods across multiple populations. The tool:

1. Reads a tab-delimited file containing allele frequencies for different populations
2. Processes the genotypes for each sample in the VCF file
3. Computes likelihood scores for each possible ancestral population
4. Assigns each sample to the population with the highest likelihood score
5. Outputs a simple mapping of sample names to assigned populations

The tool uses a statistical approach that considers the probability of observing each genotype given the population-specific allele frequencies. For each genotype:
- Homozygous reference (0/0): P = (1-f)²
- Heterozygous (0/1): P = 2f(1-f) 
- Homozygous alternate (1/1): P = f²

Where f is the frequency of the alternate allele in a given population.

## Output Format

The output is a tab-delimited text file with the following columns:

```
Sample  Assigned_Population
```

Where:
- Sample is the sample name from the VCF file
- Assigned_Population is the ancestral population with the highest likelihood score

## Examples

### Basic Usage

```bash
./VCFX_ancestry_assigner --assign-ancestry population_freqs.tsv < samples.vcf > ancestry_assignments.txt
```

### Creating a Frequency Reference File

The frequency file should have the following tab-delimited format:
```
CHROM  POS  REF  ALT  EUR  ASN  AFR
chr1   10000  A    G    0.1  0.2  0.3
chr1   20000  C    T    0.2  0.3  0.4
chr2   15000  T    C    0.4  0.5  0.6
```

### Using in a Pipeline

```bash
# Process VCF and append ancestry information as a new column in a metadata file
cat input.vcf | ./VCFX_ancestry_assigner --assign-ancestry freq.tsv | \
  join -t $'\t' metadata.txt - > metadata_with_ancestry.txt
```

## Algorithm

The ancestry assignment uses a maximum likelihood approach:

1. For each sample and each variant:
   - Determine the sample's genotype (0/0, 0/1, or 1/1)
   - For each population, calculate the log-likelihood of observing that genotype:
     - Log(P(0/0|pop)) = 2 * log(1-f)
     - Log(P(0/1|pop)) = log(2) + log(f) + log(1-f)
     - Log(P(1/1|pop)) = 2 * log(f)
   - Add this log-likelihood to the population's cumulative score

2. After processing all variants:
   - For each sample, identify the population with the highest cumulative log-likelihood
   - Assign the sample to that population

This approach is statistically sound and accounts for the probability distribution of genotypes under Hardy-Weinberg equilibrium.

## Handling Special Cases

- **Missing genotypes**: Genotypes denoted as "./." are skipped and don't contribute to likelihood calculations
- **Multi-allelic variants**: Treated as biallelic by considering only the first alternate allele
- **Missing variants**: Variants present in the VCF but not in the frequency file are skipped
- **Phased genotypes**: Phase information is ignored; both "0|1" and "0/1" are treated identically
- **Equal likelihoods**: If two populations have exactly the same likelihood (rare), the first one is assigned
- **No matching variants**: If a sample has no variants that match the frequency file, it's assigned to a default population
- **Non-standard genotypes**: Any genotype other than 0/0, 0/1, or 1/1 is skipped
- **Empty VCF**: Will produce no output rows (empty output file)

## Performance

The tool is optimized for efficiency:
- Uses hash maps for fast lookup of variant frequency data
- Single-pass processing of the VCF file
- Calculates log-likelihoods to avoid numerical underflow with many variants
- Memory usage scales with:
  - The number of variants in the frequency file
  - The number of reference populations
  - The number of samples in the VCF

## Limitations

- Requires a pre-existing set of population-specific allele frequencies
- Assumes Hardy-Weinberg equilibrium for probability calculations
- Does not account for linkage disequilibrium between variants
- Cannot detect admixed individuals (assigns to a single population)
- No confidence metrics for population assignment
- Not designed for structural variants or complex multi-allelic sites
- No support for non-diploid genotypes or unusual ploidy
- Performance depends on the number and informativeness of the variants in the frequency file 
