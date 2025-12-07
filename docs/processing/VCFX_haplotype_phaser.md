# VCFX_haplotype_phaser

## Overview

VCFX_haplotype_phaser analyzes genotype data to identify variants in linkage disequilibrium (LD) and groups them into haplotype blocks based on an LD threshold. This tool is useful for identifying sets of genetic variants that tend to be inherited together.

**New in v1.2**: The tool now supports **memory-mapped file I/O** for 20-50x faster processing, along with streaming mode for bounded memory usage.

## Usage

```bash
VCFX_haplotype_phaser [OPTIONS] [input.vcf]
VCFX_haplotype_phaser [OPTIONS] < input.vcf > blocks.txt
```

## Options

| Option | Description |
|--------|-------------|
| `-i`, `--input FILE` | Input VCF file (uses fast memory-mapped I/O) |
| `-l`, `--ld-threshold <VALUE>` | r² threshold for LD-based grouping (0.0-1.0, default: 0.8) |
| `-s`, `--streaming` | Enable streaming mode: uses sliding window for bounded memory usage |
| `-w`, `--window <N>` | Window size for streaming mode (default: 1000 variants) |
| `-q`, `--quiet` | Suppress warning messages (useful for batch processing) |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_haplotype_phaser identifies haplotype blocks by grouping variants that exhibit high linkage disequilibrium (LD). The tool:

1. Reads a VCF file from standard input
2. Extracts genotype information for each variant across all samples
3. Calculates pairwise LD (specifically r²) between consecutive variants
4. Groups variants into blocks when they exceed the specified LD threshold
5. Outputs block information showing which variants belong to the same haplotype block

The tool offers two processing modes:

### Default Mode
- Loads all variants into memory
- Calculates pairwise LD between consecutive variants
- Groups variants into blocks at the end of processing
- Memory usage scales with total number of variants

### Streaming Mode (`--streaming`)
- Uses a sliding window approach
- Outputs blocks incrementally as they are identified
- Memory usage: O(window × samples)
- **Enables processing of arbitrarily large files with bounded memory**

This is valuable for:
- Identifying haplotypes without requiring family data
- Understanding the correlation structure of variants in a genomic region
- Reducing the dimensionality of genetic data for association testing
- Planning genotyping strategies by selecting tag SNPs from different blocks
- **Processing very large VCF files (50GB+)**

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

### Basic Usage (Fast File Input)

```bash
./VCFX_haplotype_phaser -i input.vcf > haplotype_blocks.txt
./VCFX_haplotype_phaser input.vcf > haplotype_blocks.txt
```

### Basic Usage (Stdin)

```bash
./VCFX_haplotype_phaser < input.vcf > haplotype_blocks.txt
```

### Custom LD Threshold

```bash
# Use a higher LD threshold for more stringent block definition
./VCFX_haplotype_phaser --ld-threshold 0.95 -i input.vcf > strict_blocks.txt
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

### Streaming Mode for Large Files

```bash
# Process large VCF with bounded memory (stdin)
./VCFX_haplotype_phaser --streaming < large.vcf > blocks.txt

# Process large VCF with bounded memory (fast file input)
./VCFX_haplotype_phaser --streaming -i large.vcf > blocks.txt
```

### Streaming with Custom Window Size

```bash
# Use smaller window for lower memory usage
./VCFX_haplotype_phaser --streaming --window 500 --ld-threshold 0.8 -i large.vcf > blocks.txt
```

### Streaming with High LD Threshold

```bash
# Stricter LD threshold in streaming mode
./VCFX_haplotype_phaser --streaming --ld-threshold 0.95 -i input.vcf > strict_blocks.txt
```

### Quiet Mode for Batch Processing

```bash
# Suppress warning messages (useful in pipelines)
./VCFX_haplotype_phaser -q -i input.vcf > blocks.txt 2>/dev/null
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

VCFX_haplotype_phaser is designed for efficiency with multiple performance optimizations:

### File Input Mode (Recommended)

When using `-i/--input` or providing a file as a positional argument, the tool uses:

- **Memory-mapped file I/O** for 20-50x faster processing vs stdin
- **SIMD-optimized line scanning** (AVX2/SSE2 with memchr fallback)
- **Zero-copy string parsing** with std::string_view
- **1MB output buffering** for reduced syscalls
- **FORMAT field caching** to avoid redundant parsing
- **SIMD-optimized LD calculation** for vectorized correlation computation

```bash
# Fast mode (recommended for files)
VCFX_haplotype_phaser -i input.vcf > blocks.txt
VCFX_haplotype_phaser input.vcf > blocks.txt

# Slower mode (stdin)
VCFX_haplotype_phaser < input.vcf > blocks.txt
```

### Performance Characteristics

| Mode | Processing Speed | Best For |
|------|------------------|----------|
| File input (`-i`) | 20-50x faster | Large files, benchmarks |
| Stdin | Baseline | Pipelines, streaming sources |
| File + streaming | 25-50x faster | Very large files with memory limits |

## Performance Improvements

### Streaming Mode with Sliding Window

When using `--streaming`, the tool implements a bounded-memory approach:

1. **Sliding window**: Maintains only the last N variants in memory (configurable via `--window`)
2. **Incremental output**: Blocks are output as soon as they are complete
3. **Memory bounded**: O(window × samples) regardless of file size
4. **Same block boundaries**: Produces identical blocks as default mode for typical data

### Performance Characteristics

| Mode | Memory Usage | Output Timing |
|------|--------------|---------------|
| Default | O(total_variants × samples) | At EOF |
| `--streaming` | O(window × samples) | Incremental |

### When to Use Streaming Mode

- Working with very large VCF files (>1GB)
- Limited available memory
- Need incremental output as file processes
- Processing millions of variants

### Memory Usage Examples

| File Size | Variants | Samples | Default Mode | Streaming (window=1000) |
|-----------|----------|---------|--------------|-------------------------|
| 100 MB | 10K | 500 | ~40 MB | ~4 MB |
| 1 GB | 100K | 1,000 | ~400 MB | ~8 MB |
| 10 GB | 1M | 2,500 | ~4 GB | ~20 MB |
| 50 GB | 5M | 5,000 | OOM | ~40 MB |

## Backward Compatibility

The tool is fully backward compatible:
- Default behavior (without `--streaming`) works exactly as before
- Existing scripts and pipelines continue to function unchanged
- The `--streaming` option is purely opt-in for large file support
- Output format is identical in both modes

## Limitations

1. Uses a naive approach to LD calculation that may not account for population structure
2. Genotypes are treated as simple sums (0/0=0, 0/1=1, 1/1=2) without considering phase
3. Cannot merge blocks that are separated by a single low-LD variant
4. No sliding window approach for finding blocks across non-adjacent variants
5. Does not handle multi-allelic variants specially (treats them as bi-allelic)
6. Does not incorporate physical distance in the blocking algorithm
7. No output option for visualizing blocks graphically 
