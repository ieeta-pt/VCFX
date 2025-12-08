# VCFX_ld_calculator

## Overview

`VCFX_ld_calculator` calculates pairwise linkage disequilibrium (LD) statistics between genetic variants in a VCF file, expressed as r² values. It can analyze variants across an entire file or within a specified genomic region.

**New in v2.0**: Extreme-performance optimization with memory-mapped I/O, SIMD-accelerated r² computation, multi-threaded matrix mode, and optional distance-based pruning. Achieves 10-100x speedup over v1.2.

## Usage

```bash
VCFX_ld_calculator [OPTIONS] < input.vcf > ld_output.txt
VCFX_ld_calculator -i input.vcf [OPTIONS] > ld_output.txt
```

## Options

| Option | Description |
|--------|-------------|
| `-i`, `--input FILE` | Input VCF file (uses memory-mapping for best performance) |
| `-r`, `--region <chr:start-end>` | Only compute LD for variants in the specified region |
| `-w`, `--window <N>` | Window size in variants (default: 1000) |
| `-d`, `--max-distance <BP>` | Max base-pair distance between pairs (0=unlimited, default: 0) |
| `-t`, `--threshold <R2>` | Only output pairs with r² >= threshold (default: 0.0) |
| `-n`, `--threads <N>` | Number of threads for matrix mode (default: auto) |
| `-m`, `--matrix` | Use matrix mode (O(M²) full pairwise matrix) |
| `-q`, `--quiet` | Suppress informational messages |
| `-h`, `--help` | Display help message and exit |
| `-v`, `--version` | Show program version and exit |

## Description

`VCFX_ld_calculator` reads a VCF file and computes the pairwise linkage disequilibrium (r²) between genetic variants. Linkage disequilibrium is a measure of the non-random association between alleles at different loci, which is important for understanding genetic structure, identifying haplotype blocks, and designing association studies.

The tool operates as follows:

1. It reads the VCF file (using memory-mapping if `-i` is specified)
2. It collects diploid genotypes for each variant, encoding them as:
   - 0: Homozygous reference (0/0)
   - 1: Heterozygous (0/1 or 1/0)
   - 2: Homozygous alternate (1/1)
   - -1: Missing data or other scenarios (including multi-allelic variants)
3. For each pair of variants within the window/region, it computes pairwise r² values using SIMD acceleration
4. It outputs the LD results (format depends on mode)

## Computation Modes

### Default Mode (Streaming)
- Uses sliding window with `std::deque`
- Only keeps `window` most recent variants in memory
- Outputs pairs incrementally as they're computed
- Memory usage: O(window × samples)
- **Enables LD analysis of files larger than available RAM**

### Matrix Mode (`--matrix`)
- Use `--matrix` flag for full pairwise matrix output
- Multi-threaded computation with `--threads` option
- Loads all variants (in region) into memory
- Computes full MxM pairwise r² matrix
- Memory usage scales with number of variants
- Best for small regions where full matrix is needed

### Distance-Based Pruning (`--max-distance`)
- **New in v2.0**: Skip pairs beyond a physical distance threshold
- Biological justification: LD decays exponentially with distance
- Beyond ~1Mb, r² is typically < 0.01
- Reduces comparisons by 90-99% while capturing biologically relevant LD

## Performance Characteristics

### v2.0 Optimizations

1. **Memory-mapped I/O**: Uses `mmap()` with `MADV_SEQUENTIAL | MADV_WILLNEED` hints for extreme read performance.

2. **SIMD-accelerated r² computation**: Processes 16-32 samples simultaneously using NEON (ARM), AVX2, or SSE2 intrinsics.

3. **Multi-threaded matrix computation**: Parallel row computation with work-stealing for optimal CPU utilization.

4. **Zero-allocation parsing**: VCF lines parsed using raw pointer arithmetic, eliminating heap allocations.

5. **Compact genotype storage**: Uses `int8_t` instead of `int` for 4x smaller memory footprint.

6. **Pre-computed statistics**: Per-variant sum and variance cached to avoid redundant calculations.

7. **Buffered output**: 4MB output buffer with direct `write()` syscalls for minimal I/O overhead.

### Complexity Analysis

| Mode | Time Complexity | Memory |
|------|-----------------|--------|
| **Streaming** | O(n × w × m/32) | O(w × m) |
| **Matrix** | O(n² × m/32/T) | O(n × m) |

Where: n = variants, w = window, m = samples, T = threads

### Benchmark Results

| Mode | Dataset | Time (v1.2) | Time (v2.0) | Speedup |
|------|---------|-------------|-------------|---------|
| Streaming | 2K variants × 2504 samples | 1.13s | ~0.2s | 5x |
| Streaming | 427K variants × 2504 samples | ~10 min | ~2 min | 5x |
| Matrix | 10K variants × 2504 samples | ~32 min | ~30s | 60x |

### Memory Usage Examples

| File Size | Variants | Samples | Matrix Mode | Streaming (w=1000) |
|-----------|----------|---------|--------------|-------------------|
| 100 MB | 10,000 | 500 | ~40 MB | ~4 MB |
| 1 GB | 100,000 | 1,000 | ~800 MB | ~8 MB |
| 4 GB | 427,000 | 2,504 | ~4 GB | ~20 MB |
| 10 GB | 1,000,000 | 2,500 | OOM | ~20 MB |

## Output Format

### Streaming Mode (Pairwise List)

```
#VAR1_CHROM	VAR1_POS	VAR1_ID	VAR2_CHROM	VAR2_POS	VAR2_ID	R2
chr1	100	rs123	chr1	200	rs456	0.8500
chr1	100	rs123	chr1	300	rs789	0.4200
chr1	200	rs456	chr1	300	rs789	0.9100
```

### Matrix Mode

```
#LD_MATRIX_START
Index/Var	chr1:100	chr1:200	chr1:300
chr1:100	1.0000	0.4000	0.2000
chr1:200	0.4000	1.0000	0.6000
chr1:300	0.2000	0.6000	1.0000
#LD_MATRIX_END
```

## Examples

### File Input Mode (Recommended)

Use memory-mapped file I/O for best performance:

```bash
VCFX_ld_calculator -i input.vcf > ld_pairs.txt
VCFX_ld_calculator -i input.vcf -w 500 -t 0.2 > ld_pairs.txt
```

### Distance-Based Pruning

Skip pairs more than 500kb apart (biology: LD decays with distance):

```bash
VCFX_ld_calculator -i input.vcf --max-distance 500000 > ld_pairs.txt
```

### Region-Specific LD

Calculate LD only for variants in a specific genomic region:

```bash
VCFX_ld_calculator -i input.vcf -r chr1:10000-20000 > ld_pairs.txt
```

### Matrix Mode for Small Regions

Get full pairwise matrix with multi-threading:

```bash
VCFX_ld_calculator -i input.vcf -m -r chr1:10000-20000 -n 8 > ld_matrix.txt
```

### High LD Pairs Only

Output only pairs with r² >= 0.8:

```bash
VCFX_ld_calculator -i input.vcf -t 0.8 > high_ld_pairs.txt
```

### Combined Options

Streaming mode with region filter, window, threshold, and distance limit:

```bash
VCFX_ld_calculator -i input.vcf -r chr1:1000000-2000000 -w 1000 -t 0.5 -d 100000 > ld_pairs.txt
```

### Stdin Mode (Backward Compatible)

For piped input or integration with other tools:

```bash
cat input.vcf | VCFX_ld_calculator > ld_pairs.txt
cat input.vcf | VCFX_ld_calculator -m > ld_matrix.txt
```

## Handling Special Cases

- **Missing Genotypes**: Samples with missing genotypes (./. or .|.) are skipped when calculating LD
- **Multi-allelic Variants**: Genotypes involving alleles beyond reference and first alternate are treated as missing
- **Single Variant**: If only one variant is found, outputs message stating no pairwise LD can be calculated
- **Empty Region**: If no variants found in region, outputs appropriate message
- **Invalid Region Format**: Displays error message with correct format example
- **Chromosome Boundaries**: In streaming mode, the window continues across chromosome boundaries unless distance pruning is enabled

## Limitations

- Only supports biallelic variants; multi-allelic variants are treated as missing data
- Requires diploid genotypes; haploid genotypes will be treated as missing data
- Assumes standard VCF format with GT field in the FORMAT column
- Does not support phased vs. unphased distinction; both "/" and "|" separators are treated the same
- No built-in visualization of LD patterns
- Window is based on variant count by default; use `--max-distance` for physical distance filtering

## Backward Compatibility

**Breaking Change in v1.2**: Default mode changed from matrix to streaming.
- To get matrix output, add `--matrix` flag
- The `--streaming` flag is deprecated (streaming is now default)

**New in v2.0**: New options added:
- `-i/--input` for mmap mode
- `-n/--threads` for multi-threaded matrix
- `-d/--max-distance` for distance pruning
- `-q/--quiet` for silent operation

All existing options continue to work as before.
