# VCFX_haplotype_extractor

## Overview

VCFX_haplotype_extractor reconstructs phased haplotype blocks from genotype data in a VCF file. It identifies stretches of phased variants on the same chromosome and outputs them as continuous haplotype blocks for each sample.

**New in v1.1**: The tool now supports **streaming mode** for bounded memory usage, enabling processing of arbitrarily large files by outputting blocks immediately when they complete.

## Usage

```bash
VCFX_haplotype_extractor [OPTIONS] < input.vcf > haplotypes.tsv
```

## Options

| Option | Description |
|--------|-------------|
| `--block-size <SIZE>` | Maximum distance in base pairs between consecutive variants to be included in the same block (default: 100,000) |
| `--check-phase-consistency` | Enable checks for phase consistency between adjacent variants in a block |
| `--streaming` | Enable streaming mode: output blocks immediately when complete. Uses O(block) memory instead of O(total_variants) |
| `--debug` | Enable verbose debug output |
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

The tool offers two processing strategies:

### Default Mode
- Accumulates all haplotype blocks in memory
- Outputs all blocks at end of file processing
- Memory usage scales with total number of blocks

### Streaming Mode (`--streaming`)
- Outputs each block immediately when complete
- Only keeps current block in memory
- Memory usage: O(samples × variants_in_block)
- **Enables processing of arbitrarily large files with bounded memory**

This tool is valuable for:
- Identifying regions of continuous phasing in VCF files
- Preparing haplotype data for downstream analyses
- Reconstructing parental chromosomes from phased variant data
- Quality control of phasing algorithms
- **Processing very large phased VCF files (50GB+)**

## Performance Improvements

### Streaming Block Output
When using `--streaming`, the tool implements immediate block output:

1. **Incremental output**: Blocks are written as soon as they complete
2. **Bounded memory**: Only current block held in memory
3. **Linear I/O**: Each line read once, output immediately when block closes
4. **No post-processing**: No need to wait for EOF

### Performance Characteristics
| Mode | Memory Usage | Output Timing |
|------|--------------|---------------|
| Default | O(total_blocks × samples) | At EOF |
| `--streaming` | O(current_block × samples) | Incremental |

### When to Use Streaming Mode
- Working with very large VCF files (>1GB)
- Limited available memory
- Need to start downstream processing before file completes
- Processing chromosome-at-a-time pipelines

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

### Streaming Mode for Large Files

```bash
./VCFX_haplotype_extractor --streaming < large_phased.vcf > haplotype_blocks.tsv
```

### Custom Block Size

```bash
# Use a smaller maximum distance (50kb) to generate more, smaller blocks
./VCFX_haplotype_extractor --block-size 50000 < phased.vcf > small_blocks.tsv
```

### Streaming with Custom Block Size

```bash
# Streaming mode with 10kb block distance
./VCFX_haplotype_extractor --streaming --block-size 10000 < large_phased.vcf > haplotypes.tsv
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
8. **Large files**: Use `--streaming` mode for files larger than available RAM

## Performance

### Default Mode
1. Single-pass processing with O(n) time complexity where n is the number of variants
2. Memory usage scales primarily with the number of samples and the total number of blocks
3. All blocks held in memory until EOF

### Streaming Mode
1. Single-pass processing with O(n) time complexity
2. Memory usage scales only with current block size (not total blocks)
3. Blocks output immediately when complete
4. Optimal for very large files or memory-constrained environments

### Memory Usage Examples

| File Size | Blocks | Samples | Default Mode | Streaming Mode |
|-----------|--------|---------|--------------|----------------|
| 100 MB | 100 | 500 | ~20 MB | ~1 MB |
| 1 GB | 1,000 | 1,000 | ~100 MB | ~5 MB |
| 10 GB | 10,000 | 2,500 | ~500 MB | ~10 MB |
| 50 GB | 50,000 | 5,000 | OOM | ~20 MB |

## Limitations

1. Requires phased genotypes - variants with unphased genotypes are skipped
2. Cannot join blocks across different chromosomes
3. Simple distance-based blocking may not align with biological recombination patterns
4. Basic phase consistency checking may not detect all inconsistencies
5. No ability to export or visualize the relationship between blocks
6. Does not account for potential errors in the original phasing
7. No special handling for reference gaps or known problematic regions

## Backward Compatibility

The tool is fully backward compatible:
- Default behavior (without `--streaming`) works exactly as before
- Existing scripts and pipelines continue to function unchanged
- The `--streaming` option is purely opt-in for large file support
- Output format is identical in both modes
