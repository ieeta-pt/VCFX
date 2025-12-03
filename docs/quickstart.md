# Quick Start Guide

This guide provides a quick introduction to using the VCFX toolkit with practical examples.

## Overview

VCFX consists of multiple small command-line tools that are designed to be combined in pipelines. Each tool follows these principles:

- Reads from standard input and writes to standard output
- Performs a single, focused task
- Can be combined with other tools via pipes
- Has a consistent command-line interface

## Basic Usage Pattern

Most VCFX tools follow this basic pattern:

```bash
VCFX_tool_name [options] < input.vcf > output.vcf
```

For example:

```bash
VCFX_variant_classifier < input.vcf > classified.vcf
```

## Python API Usage

If you've installed VCFX via PyPI (`pip install vcfx`), you can use the Python API:

```python
import vcfx

# Count variants
n_variants = vcfx.variant_counter("input.vcf")
print(f"Total variants: {n_variants}")

# Get allele frequencies as structured data
frequencies = vcfx.allele_freq_calc("input.vcf")
for freq in frequencies:
    print(f"Chr {freq.Chromosome}, Pos {freq.Pos}: AF={freq.Allele_Frequency}")

# Check sample concordance
concordance = vcfx.concordance_checker("input.vcf", "SAMPLE1", "SAMPLE2")
discordant = [r for r in concordance if r.Concordance != "Concordant"]
print(f"Found {len(discordant)} discordant sites")
```

Note: The Python API requires the VCFX command-line tools to be installed and available in your PATH.

## Common Examples

Here are some common use cases for VCFX tools:

### Example 1: Basic Filtering

Filter for high-quality SNPs:

```bash
cat input.vcf | \
  VCFX_variant_classifier --append-info | \
  grep 'VCF_CLASS=SNP' | \
  VCFX_phred_filter --phred-filter 30 > high_quality_snps.vcf
```

### Example 2: Population Analysis

Extract European samples and calculate allele frequencies:

```bash
cat input.vcf | \
  VCFX_population_filter --population EUR --pop-map populations.txt | \
  VCFX_allele_freq_calc > eur_frequencies.tsv
```

### Example 3: Data Transformation

Normalize indels and split multiallelic variants:

```bash
cat input.vcf | \
  VCFX_indel_normalizer | \
  VCFX_multiallelic_splitter > normalized_biallelic.vcf
```

### Example 4: Quality Control

Check concordance between two samples in a single VCF file:

```bash
cat sample.vcf | VCFX_concordance_checker --samples "SAMPLE1 SAMPLE2" > concordance_report.tsv
```

## Working with Compressed Files

Most VCFX tools don't directly support compressed input/output. Use standard Unix tools:

```bash
# Reading from compressed files
zcat input.vcf.gz | VCFX_tool_name > output.vcf

# Writing to compressed files
VCFX_tool_name < input.vcf | gzip > output.vcf.gz
```

## Getting Help

All VCFX tools provide help information:

```bash
VCFX_tool_name --help
```

This will display the tool's purpose, options, and usage examples.

## Tool Categories

VCFX tools are categorized by their function:

### Data Analysis
Tools for extracting information from VCF files (e.g., `VCFX_allele_freq_calc`)

### Data Filtering
Tools for selecting variants based on criteria (e.g., `VCFX_phred_filter`)

### Data Transformation
Tools for converting or reformatting VCF data (e.g., `VCFX_indel_normalizer`)

### Quality Control
Tools for validating and checking data quality (e.g., `VCFX_validator`)

### File Management
Tools for handling VCF files (e.g., `VCFX_indexer`)

## Common Workflows

Here are some common workflows that combine multiple VCFX tools:

### Variant QC Pipeline

```bash
cat input.vcf | \
  VCFX_validator | \
  VCFX_variant_classifier --append-info | \
  VCFX_missing_detector --max-missing 0.1 | \
  VCFX_phred_filter --phred-filter 20 > qc_passed.vcf
```

### Sample Comparison

```bash
# Check concordance between two samples in a single VCF
cat input.vcf | VCFX_concordance_checker --samples "SAMPLE1 SAMPLE2" > concordance.tsv
```

### Population Structure Analysis

```bash
# Extract population-specific VCFs
cat input.vcf | VCFX_population_filter --population EUR --pop-map pop_map.txt > eur.vcf
cat input.vcf | VCFX_population_filter --population AFR --pop-map pop_map.txt > afr.vcf

# Calculate allele frequencies for each population
cat eur.vcf | VCFX_allele_freq_calc > eur_afs.tsv
cat afr.vcf | VCFX_allele_freq_calc > afr_afs.tsv
```

### Python Workflow Example

Here's the same population analysis using the Python API:

```python
import vcfx
import pandas as pd

# Filter by population and calculate frequencies
eur_vcf = vcfx.population_filter("input.vcf", population="EUR", pop_map="pop_map.txt")
afr_vcf = vcfx.population_filter("input.vcf", population="AFR", pop_map="pop_map.txt")

# Get allele frequencies as structured data
eur_freqs = vcfx.allele_freq_calc(eur_vcf)
afr_freqs = vcfx.allele_freq_calc(afr_vcf)

# Compare frequencies
eur_df = pd.DataFrame(eur_freqs)
afr_df = pd.DataFrame(afr_freqs)

# Find variants with large frequency differences
merged = eur_df.merge(afr_df, on=['Chromosome', 'Pos'], suffixes=('_EUR', '_AFR'))
merged['freq_diff'] = abs(merged['Allele_Frequency_EUR'] - merged['Allele_Frequency_AFR'])
differentiated = merged[merged['freq_diff'] > 0.2]

print(f"Found {len(differentiated)} highly differentiated variants")
```

## Next Steps

After becoming familiar with the basic usage of VCFX tools, you can:

1. Explore the [complete tool documentation](tools_overview.md) for details on each tool
2. Check the [installation guide](installation.md) if you need to install additional tools
3. Browse through the example VCF files in the repository to practice

For more complex workflows and advanced examples, refer to the individual tool documentation pages. 
