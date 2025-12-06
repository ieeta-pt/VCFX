# VCFX_ld_calculator

## Overview

`VCFX_ld_calculator` calculates pairwise linkage disequilibrium (LD) statistics between genetic variants in a VCF file, expressed as r² values. It can analyze variants across an entire file or within a specified genomic region.

**New in v1.2**: **Streaming mode is now the default** for optimal performance with large files. The streaming algorithm provides 6x faster processing on 2000-variant files and up to 400x theoretical speedup on genome-wide files. Use `--matrix` for backward-compatible full matrix output.

## Usage

```bash
VCFX_ld_calculator [OPTIONS] < input.vcf > ld_output.txt
```

## Options

| Option | Description |
|--------|-------------|
| `--region <chr:start-end>` | Only compute LD for variants in the specified region |
| `--matrix` | Use matrix mode (O(M²) full pairwise matrix) - for backward compatibility |
| `--window <N>` | Window size: compute LD only between variants within N positions of each other (default: 1000) |
| `--threshold <R2>` | Only output pairs with r² >= threshold (default: 0.0, output all pairs) |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

`VCFX_ld_calculator` reads a VCF file and computes the pairwise linkage disequilibrium (r²) between genetic variants. Linkage disequilibrium is a measure of the non-random association between alleles at different loci, which is important for understanding genetic structure, identifying haplotype blocks, and designing association studies.

The tool operates as follows:

1. It reads the VCF file from standard input
2. It collects diploid genotypes for each variant, encoding them as:
   - 0: Homozygous reference (0/0)
   - 1: Heterozygous (0/1 or 1/0)
   - 2: Homozygous alternate (1/1)
   - -1: Missing data or other scenarios (including multi-allelic variants)
3. For each pair of variants within the specified region (or the entire file if no region is specified), it computes pairwise r² values, ignoring samples with missing genotypes
4. It outputs the LD results (format depends on mode)

The tool offers two computation strategies:

### Default Mode (Streaming) - NEW in v1.2
- **Now the default** for optimal performance
- Uses sliding window with `std::deque`
- Only keeps `window` most recent variants in memory
- Outputs pairs incrementally as they're computed
- Memory usage: O(window × samples)
- **Enables LD analysis of files larger than available RAM**
- **6x faster** on 2000-variant files, **400x theoretical speedup** on genome-wide files

### Matrix Mode (`--matrix`) - Backward Compatible
- Use `--matrix` flag for backward compatibility
- Loads all variants (in region) into memory
- Computes full MxM pairwise r² matrix
- Outputs symmetric matrix format
- Memory usage scales with number of variants
- Best for small regions where full matrix is needed

The r² calculation uses the standard formula:
- Let X and Y be the genotype arrays for two variants
- Calculate means of X and Y (meanX, meanY)
- Calculate covariance: cov = average(XY) - meanX × meanY
- Calculate variances: varX = average(X²) - meanX², varY similarly
- r = cov / sqrt(varX × varY)
- r² = r × r

## Performance Improvements

### Sliding Window Algorithm
When using `--streaming`, the tool implements memory-efficient LD calculation:

1. **Sliding window**: Only keeps the most recent N variants in memory
2. **Bounded memory**: O(window × samples) regardless of file size
3. **Incremental output**: LD pairs are output as they're computed
4. **Threshold filtering**: Skip pairs below r² threshold to reduce output size

### Performance Characteristics
| Mode | Memory Usage | Output Format |
|------|--------------|---------------|
| Default | O(variants × samples) | MxM matrix |
| `--streaming` | O(window × samples) | Pairwise list |

### When to Use Streaming Mode
- Working with very large VCF files (>1GB)
- Limited available memory
- Only need LD between nearby variants (not genome-wide matrix)
- Want to filter by r² threshold on-the-fly

### Recommended Settings

| Scenario | Recommended Options |
|----------|---------------------|
| Small region analysis | Default mode (no `--streaming`) |
| Large file, local LD | `--streaming --window 1000` |
| Memory-constrained | `--streaming --window 500` |
| High LD pairs only | `--streaming --threshold 0.8` |
| GWAS LD pruning | `--streaming --window 500 --threshold 0.2` |

## Output Format

### Default Mode (Matrix)

The output is a tab-delimited matrix of r² values with a header identifying the variants:

```
#LD_MATRIX_START
         chr1:100 chr1:200 chr1:300
chr1:100      1.0     0.4     0.2
chr1:200     0.4      1.0     0.6
chr1:300     0.2     0.6      1.0
```

### Streaming Mode (Pairwise List)

The output is a tab-delimited list of variant pairs with r² values:

```
#VAR1_CHROM	VAR1_POS	VAR1_ID	VAR2_CHROM	VAR2_POS	VAR2_ID	R2
chr1	100	rs123	chr1	200	rs456	0.8500
chr1	100	rs123	chr1	300	rs789	0.4200
chr1	200	rs456	chr1	300	rs789	0.9100
```

If only one or no variants are found in the region, the tool outputs a message indicating that no pairwise LD could be calculated.

## Examples

### Basic Usage (Streaming Default)

Calculate LD for all variants in a VCF file (streaming mode is default):

```bash
VCFX_ld_calculator < input.vcf > ld_pairs.txt
```

### Region-Specific LD

Calculate LD only for variants in a specific genomic region:

```bash
VCFX_ld_calculator --region chr1:10000-20000 < input.vcf > ld_pairs.txt
```

### Matrix Mode for Small Regions

Get full pairwise matrix (backward compatible):

```bash
VCFX_ld_calculator --matrix --region chr1:10000-20000 < input.vcf > ld_matrix.txt
```

### Custom Window Size

Limit LD calculation to variants within 500 positions of each other:

```bash
VCFX_ld_calculator --window 500 < input.vcf > ld_pairs.txt
```

### High LD Pairs Only

Output only pairs with r² >= 0.8:

```bash
VCFX_ld_calculator --streaming --threshold 0.8 < input.vcf > high_ld_pairs.txt
```

### Combined Options

Streaming mode with region filter, window, and threshold:

```bash
VCFX_ld_calculator --region chr1:1000000-2000000 --streaming --window 1000 --threshold 0.5 < input.vcf > ld_pairs.txt
```

### Integration with Other Tools

Filter for common variants first, then calculate LD:

```bash
cat input.vcf | VCFX_af_subsetter --af-filter '0.05-1.0' | VCFX_ld_calculator > common_variants_ld.txt
```

### LD Pruning for GWAS

Find variant pairs in high LD for pruning:

```bash
VCFX_ld_calculator --streaming --window 500 --threshold 0.2 < gwas_data.vcf > ld_prune_pairs.txt
```

## Handling Special Cases

- **Missing Genotypes**: Samples with missing genotypes (./. or .|.) are skipped when calculating LD between variant pairs
- **Multi-allelic Variants**: Genotypes involving alleles beyond the reference and first alternate (e.g., 1/2, 2/2) are treated as missing data
- **Single Variant**: If only one variant is found in the region, the tool outputs a message stating that no pairwise LD can be calculated
- **Empty Region**: If no variants are found in the specified region, the tool outputs a message stating that no pairwise LD can be calculated
- **Invalid Region Format**: If the region format is invalid, the tool will display an error message
- **Large Files**: Use `--streaming` mode for files larger than available RAM
- **Chromosome Boundaries**: In streaming mode, the window resets at chromosome boundaries

## Performance

### Optimized Implementation
The streaming mode uses several performance optimizations:

1. **Zero-allocation parsing**: VCF lines are parsed using raw pointer arithmetic instead of creating intermediate string objects, eliminating millions of heap allocations.

2. **Compact genotype storage**: Uses `int8_t` instead of `int` for genotypes (4x smaller memory footprint).

3. **Fast r² computation**: Pre-computed per-variant statistics (sum, sum of squares) reduce redundant calculations.

4. **Fast number formatting**: Custom formatters for integers and r² values replace slow iostream manipulators.

5. **Output buffering**: 1MB output buffer reduces system call overhead.

### Complexity

**Streaming Mode (Default)**
- Time complexity: O(n × w × m) where n is variants, w is window size, m is samples
- Memory usage: O(w × m) - constant with respect to file size
- Each variant is read once and evicted from window when no longer needed
- Optimal for genome-wide LD calculation with local windows

**Matrix Mode**
- Time complexity: O(n²m) where n is the number of variants and m is the number of samples
- Memory usage scales linearly with the number of variants and samples
- For large datasets, use the `--region` option to limit the analysis

### Benchmark Results

| Dataset | Variants | Samples | Pairs | Time |
|---------|----------|---------|-------|------|
| chr21 subset | 2,000 | 100 | 194,950 | 0.10s |
| chr21 subset | 2,000 | 2,504 | 194,950 | 1.13s |

Throughput: ~170,000 LD pairs/second on 2504 samples.

### Memory Usage Examples

| File Size | Variants | Samples | Matrix Mode | Streaming (w=1000) |
|-----------|----------|---------|--------------|-------------------|
| 100 MB | 10,000 | 500 | ~40 MB | ~4 MB |
| 1 GB | 100,000 | 1,000 | ~800 MB | ~8 MB |
| 10 GB | 1,000,000 | 2,500 | OOM | ~20 MB |
| 50 GB | 5,000,000 | 5,000 | OOM | ~40 MB |

## Limitations

- Only supports biallelic variants; multi-allelic variants are treated as missing data
- Requires diploid genotypes; haploid genotypes will be treated as missing data
- Assumes standard VCF format with GT field in the FORMAT column
- Does not support phased vs. unphased distinction; both "/" and "|" separators are treated the same
- No built-in visualization of LD patterns; additional tools would be needed for heatmap creation
- **Streaming mode outputs pairwise format, not matrix format**
- **Window is based on variant count, not physical distance (bp)**

## Backward Compatibility

**Breaking Change in v1.2**: The default mode has changed from matrix to streaming for better performance:

- **Old behavior** (v1.1 and earlier): `VCFX_ld_calculator < file.vcf` → MxM matrix output
- **New behavior** (v1.2+): `VCFX_ld_calculator < file.vcf` → pairwise streaming output

**Migration Guide**:
- To preserve old behavior, add `--matrix` flag: `VCFX_ld_calculator --matrix < file.vcf`
- Scripts expecting matrix format should add `--matrix`
- The `--streaming` flag is now deprecated (streaming is default); use absence of `--matrix` instead

**Why this change?**
- 6x faster on 2000-variant files
- Up to 400x faster on genome-wide files
- O(window × samples) memory instead of O(variants × samples)
- Enables LD analysis of 50GB+ files that previously caused OOM
