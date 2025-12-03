# VCFX_haplotype_extractor

## Overview

VCFX_haplotype_extractor reconstructs phased haplotype blocks from genotype data in a VCF file. It identifies stretches of phased variants on the same chromosome and outputs them as continuous haplotype blocks for each sample.

## Usage

```bash
VCFX_haplotype_extractor [OPTIONS] < input.vcf > haplotypes.tsv
```

## Options

| Option | Description |
|--------|-------------|
| `--block-size <SIZE>` | Maximum distance in base pairs between consecutive variants to be included in the same block (default: 100,000) |
| `--check-phase-consistency` | Enable checks for phase consistency between adjacent variants in a block |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_haplotype_extractor analyzes phased genotype data in a VCF file to reconstruct continuous haplotype blocks. The tool:

1. Reads a VCF file from standard input
2. Extracts phased genotype (GT) fields for each sample at each variant position
3. Groups consecutive phased variants into blocks based on:
   - Chromosome continuity (variants must be on the same chromosome)
   - Maximum distance threshold (default 100kb between adjacent variants)
   - Optional phase consistency checks across variants
4. Constructs haplotype strings representing the sequence of alleles on each chromosome
5. Outputs blocks of phased haplotypes in a tab-delimited format

This tool is valuable for:
- Identifying regions of continuous phasing in VCF files
- Preparing haplotype data for downstream analyses
- Reconstructing parental chromosomes from phased variant data
- Quality control of phasing algorithms

## Output Format

The output is a tab-delimited text file with columns:

```
CHROM  START  END  SAMPLE_1_HAPLOTYPES  SAMPLE_2_HAPLOTYPES  ...
```

Where:
- CHROM: Chromosome name
- START: Start position of the haplotype block
- END: End position of the haplotype block
- SAMPLE_X_HAPLOTYPES: A pipe-delimited string representing the phased genotypes for that sample

Each sample's haplotype column contains a string of pipe-separated genotypes where each genotype is itself pipe-separated (e.g., "0|1|1|0|0|1"). This represents the sequence of alleles in the phased block.

## Examples

### Basic Usage

```bash
./VCFX_haplotype_extractor < phased.vcf > haplotype_blocks.tsv
```

### Custom Block Size

```bash
# Use a smaller maximum distance (50kb) to generate more, smaller blocks
./VCFX_haplotype_extractor --block-size 50000 < phased.vcf > small_blocks.tsv
```

### With Phase Consistency Checking

```bash
# Enable checks for phase consistency between variants
./VCFX_haplotype_extractor --check-phase-consistency < phased.vcf > consistent_blocks.tsv
```

### Filtering for Large Blocks

```bash
# Extract only blocks spanning at least 10 variants
./VCFX_haplotype_extractor < phased.vcf | awk -F'|' '{if (NF >= 10) print}' > large_blocks.tsv
```

## Phase Consistency

When the `--check-phase-consistency` option is enabled, the tool performs a basic check to detect potential phase inconsistencies:

1. For each new variant, the tool examines its phased alleles for each sample
2. It compares these with the last variant added to the current block
3. If it detects a phase "flip" (e.g., changing from "0|1" to "1|0"), it may start a new block
4. This helps identify regions where phasing may be inconsistent

This basic consistency checking is useful for identifying phase switches that might indicate errors in the phasing process or real recombination events.

## Handling Special Cases

The tool implements several strategies for handling edge cases:

1. **Unphased genotypes**: Any variant with unphased genotypes ("/" delimiter) is skipped and will not be included in haplotype blocks
2. **Missing genotypes**: Variants with missing genotypes (".") are handled, but may affect block formation
3. **Multiallelic sites**: Properly processed with the actual allele codes in the haplotype strings
4. **Chromosome changes**: Automatically starts a new block when the chromosome changes
5. **Large distances**: Starts a new block when the distance between consecutive variants exceeds the threshold
6. **Empty input**: Produces no output blocks but exits cleanly
7. **Malformed VCF**: Attempts to skip malformed lines with warnings

## Performance

VCFX_haplotype_extractor is designed for efficiency:

1. Single-pass processing with O(n) time complexity where n is the number of variants
2. Memory usage scales primarily with the number of samples and the size of the largest haplotype block
3. Streaming architecture allows processing large files without loading them entirely into memory
4. Block-based approach prevents excessive memory usage for very long chromosomes

## Limitations

1. Requires phased genotypes - variants with unphased genotypes are skipped
2. Cannot join blocks across different chromosomes
3. Simple distance-based blocking may not align with biological recombination patterns
4. Basic phase consistency checking may not detect all inconsistencies
5. No ability to export or visualize the relationship between blocks
6. Does not account for potential errors in the original phasing
7. No special handling for reference gaps or known problematic regions 
