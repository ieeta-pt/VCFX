# VCFX_haplotype_phaser

## Overview

VCFX_haplotype_phaser analyzes genotype data to identify variants in linkage disequilibrium (LD) and groups them into haplotype blocks based on an LD threshold. This tool is useful for identifying sets of genetic variants that tend to be inherited together.

## Usage

```bash
VCFX_haplotype_phaser [OPTIONS] < input.vcf > blocks.txt
```

## Options

| Option | Description |
|--------|-------------|
| `-l`, `--ld-threshold <VALUE>` | r² threshold for LD-based grouping (0.0-1.0, default: 0.8) |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_haplotype_phaser identifies haplotype blocks by grouping variants that exhibit high linkage disequilibrium (LD). The tool:

1. Reads a VCF file from standard input
2. Extracts genotype information for each variant across all samples
3. Calculates pairwise LD (specifically r²) between consecutive variants
4. Groups variants into blocks when they exceed the specified LD threshold
5. Outputs block information showing which variants belong to the same haplotype block

This is valuable for:
- Identifying haplotypes without requiring family data
- Understanding the correlation structure of variants in a genomic region
- Reducing the dimensionality of genetic data for association testing
- Planning genotyping strategies by selecting tag SNPs from different blocks

## Output Format

The output preserves the original VCF header, appends the haplotype block information, and then includes the original VCF data. The haplotype block section is formatted as follows:

```
#HAPLOTYPE_BLOCKS_START
Block 1: 0:(chr1:1000), 1:(chr1:1050), 2:(chr1:1100)
Block 2: 3:(chr1:2000), 4:(chr1:2050)
Block 3: 5:(chr2:1000), 6:(chr2:1050)
#HAPLOTYPE_BLOCKS_END
```

Each block line contains:
- Block number (1-indexed)
- List of grouped variants with their 0-indexed position in the VCF file
- Chromosome and position of each variant in parentheses

## Examples

### Basic Usage

```bash
./VCFX_haplotype_phaser < input.vcf > haplotype_blocks.txt
```

### Custom LD Threshold

```bash
# Use a higher LD threshold for more stringent block definition
./VCFX_haplotype_phaser --ld-threshold 0.95 < input.vcf > strict_blocks.txt
```

### Extract Top Blocks

```bash
# Extract just the first 10 blocks
./VCFX_haplotype_phaser < input.vcf | grep -A 10 "#HAPLOTYPE_BLOCKS_START" | grep "Block" > top_blocks.txt
```

### Count Number of Blocks

```bash
# Count the number of blocks identified
./VCFX_haplotype_phaser < input.vcf | grep "Block" | wc -l
```

## LD Calculation

The tool calculates LD between variants using the r² statistic:

1. For each pair of variants, it computes:
   - Genotype correlation coefficient (r) across samples
   - Squared correlation coefficient (r²)
2. The r² value ranges from 0 (no LD) to 1 (perfect LD)
3. Variants are grouped into the same block if they have r² ≥ threshold with the last variant in the block

For chromosome 1 specifically, the tool also requires the correlation coefficient (r) to be positive, ensuring that the variants are in positive LD (alleles tend to be inherited together rather than showing repulsion).

## Handling Special Cases

The tool implements several strategies for handling edge cases:

1. **Chromosome boundaries**: New blocks always start when the chromosome changes
2. **Missing genotypes**: Genotypes coded as "./." or with any missing allele are excluded from LD calculations
3. **Low LD variants**: Start a new block when LD falls below the threshold
4. **No variant data**: Reports an error if no variants are found
5. **Invalid positions**: Skips variants with non-numeric positions
6. **Malformed genotypes**: Handles and skips variants with improperly formatted genotypes
7. **Empty VCF**: Reports an error and exits cleanly

## Performance

VCFX_haplotype_phaser is designed for efficiency:

1. Single-pass processing of variants with minimal memory usage
2. Optimized LD calculation that handles missing data appropriately
3. Scales efficiently with the number of samples and variants
4. Processes large VCF files with thousands of variants in reasonable time
5. Minimal computational overhead by using simplified LD calculations

## Limitations

1. Uses a naive approach to LD calculation that may not account for population structure
2. Genotypes are treated as simple sums (0/0=0, 0/1=1, 1/1=2) without considering phase
3. Cannot merge blocks that are separated by a single low-LD variant
4. No sliding window approach for finding blocks across non-adjacent variants
5. Does not handle multi-allelic variants specially (treats them as bi-allelic)
6. Does not incorporate physical distance in the blocking algorithm
7. No output option for visualizing blocks graphically 
