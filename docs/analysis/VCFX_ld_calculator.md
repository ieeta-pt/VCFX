# VCFX_ld_calculator

## Overview

`VCFX_ld_calculator` calculates pairwise linkage disequilibrium (LD) statistics between genetic variants in a VCF file, expressed as r² values. It can analyze variants across an entire file or within a specified genomic region.

**New in v1.1**: The tool now supports **streaming mode** with a sliding window for handling large files with bounded memory usage, enabling LD analysis of 50GB+ files without loading all variants into memory.

## Usage

```bash
VCFX_ld_calculator [OPTIONS] < input.vcf > ld_output.txt
```

## Options

| Option | Description |
|--------|-------------|
| `--region <chr:start-end>` | Only compute LD for variants in the specified region |
| `--streaming` | Enable streaming mode with sliding window (constant memory usage) |
| `--window <N>` | Window size for streaming mode: compute LD only between variants within N positions of each other (default: 1000) |
| `--threshold <R2>` | In streaming mode, only output pairs with r² >= threshold (default: 0.0, output all pairs) |
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

### Default Mode (MxM Matrix)
- Loads all variants (in region) into memory
- Computes full MxM pairwise r² matrix
- Outputs symmetric matrix format
- Memory usage scales with number of variants

### Streaming Mode (`--streaming`)
- Uses sliding window with `std::deque`
- Only keeps `window` most recent variants in memory
- Outputs pairs incrementally as they're computed
- Memory usage: O(window × samples)
- **Enables LD analysis of files larger than available RAM**

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

### Basic Usage

Calculate LD for all variants in a VCF file:

```bash
VCFX_ld_calculator < input.vcf > ld_matrix.txt
```

### Region-Specific LD

Calculate LD only for variants in a specific genomic region:

```bash
VCFX_ld_calculator --region chr1:10000-20000 < input.vcf > ld_matrix.txt
```

### Streaming Mode for Large Files

Calculate LD with bounded memory usage:

```bash
VCFX_ld_calculator --streaming < large_file.vcf > ld_pairs.txt
```

### Streaming with Custom Window Size

Limit LD calculation to variants within 500 positions of each other:

```bash
VCFX_ld_calculator --streaming --window 500 < input.vcf > ld_pairs.txt
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

### Default Mode
- Time complexity is O(n²m) where n is the number of variants and m is the number of samples
- Memory usage scales linearly with the number of variants and samples
- For large datasets with many variants, consider using the `--region` option to limit the analysis to specific genomic regions

### Streaming Mode
- Time complexity is O(n × w × m) where n is variants, w is window size, m is samples
- Memory usage is O(w × m) - constant with respect to file size
- Each variant is read once and evicted from window when no longer needed
- Optimal for genome-wide LD calculation with local windows

### Memory Usage Examples

| File Size | Variants | Samples | Default Mode | Streaming (w=1000) |
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

The tool is fully backward compatible:
- Default behavior (without `--streaming`) works exactly as before, producing MxM matrix
- Existing scripts and pipelines continue to function unchanged
- The `--streaming`, `--window`, and `--threshold` options are purely opt-in for large file support
