# VCFX_hwe_tester

## Overview

VCFX_hwe_tester performs Hardy-Weinberg Equilibrium (HWE) testing on biallelic variants in a VCF file, calculating and reporting p-values that measure the degree of deviation from expected genotype frequencies.

## Usage

```bash
VCFX_hwe_tester [OPTIONS] [input.vcf]
VCFX_hwe_tester [OPTIONS] < input.vcf > hwe_results.txt
```

## Options

| Option | Description |
|--------|-------------|
| `-i`, `--input FILE` | Input VCF file (uses memory-mapping for best performance) |
| `-q`, `--quiet` | Suppress informational messages |
| `-h`, `--help` | Display help message and exit |
| `-v`, `--version` | Show program version and exit |

## Description

VCFX_hwe_tester analyzes each biallelic variant in a VCF file to determine whether its genotype frequencies conform to Hardy-Weinberg Equilibrium expectations. The tool:

1. Reads the VCF file line by line
2. Filters for biallelic variants (skips sites with multiple ALT alleles)
3. Counts homozygous reference (0/0), heterozygous (0/1), and homozygous alternate (1/1) genotypes
4. Calculates a p-value for Hardy-Weinberg Equilibrium using chi-square test with Yates' correction
5. Reports results in a simple tab-delimited format

A low p-value indicates significant deviation from Hardy-Weinberg Equilibrium, which might suggest:
- Population stratification
- Selection pressure
- Non-random mating
- Genotyping errors

## Output Format

The output is a tab-delimited text file with the following columns:

```
CHROM  POS  ID  REF  ALT  HWE_pvalue
```

Where:
- CHROM, POS, ID, REF, ALT are copied from the input VCF
- HWE_pvalue is the calculated p-value for Hardy-Weinberg Equilibrium

## Examples

### File Input Mode (Recommended)

Use memory-mapped file I/O for best performance:

```bash
VCFX_hwe_tester -i input.vcf > hwe_results.txt
VCFX_hwe_tester input.vcf > hwe_results.txt
```

### Basic Usage (Stdin)

```bash
VCFX_hwe_tester < input.vcf > hwe_results.txt
```

### Filter by HWE p-value

```bash
# Extract variants with significant HWE deviation (p < 0.05)
VCFX_hwe_tester -i input.vcf | awk -F'\t' '{if(NR==1 || $6<0.05) print}' > hwe_significant.txt
```

### Check for Genotyping Errors

```bash
# Find potential genotyping errors (extremely low HWE p-values)
VCFX_hwe_tester -i input.vcf | awk -F'\t' '{if(NR==1 || $6<0.0001) print}' > potential_errors.txt
```

## Mathematical Details

The tool uses a chi-square test with Yates' continuity correction. For each variant:

1. The observed counts of genotypes (homRef, het, homAlt) are calculated
2. The expected frequencies under HWE are computed as:
   - Expected homRef = p² × N
   - Expected het = 2pq × N
   - Expected homAlt = q² × N

   where p = (2×homRef + het)/(2×N), q = 1-p, and N = total number of individuals

3. Chi-square statistic is computed with Yates' correction for continuity
4. P-value is calculated using the chi-square distribution (1 degree of freedom)

This approach provides fast and accurate p-values for large sample sizes.

## Handling Special Cases

- **Multi-allelic variants**: Skipped entirely (only biallelic variants are considered)
- **Missing genotypes**: Excluded from counts when calculating HWE
- **Phased genotypes**: Phase information is ignored; "0|1" is treated the same as "0/1"
- **Non-standard genotypes**: Any genotype other than 0/0, 0/1, or 1/1 is excluded
- **No valid genotypes**: If no valid genotypes are found, the p-value is reported as 1.0
- **Perfect equilibrium**: For variants with genotype frequencies perfectly matching HWE expectations, the p-value is 1.0
- **No GT field**: Lines without GT in FORMAT field are skipped

## Performance Characteristics

The tool uses several optimizations for high performance:

- **Memory-mapped I/O**: Uses `mmap()` for file input to minimize syscall overhead
- **SIMD acceleration**: Uses NEON/SSE2/AVX2 instructions for fast line/tab scanning
- **Zero-copy parsing**: Parses VCF fields without creating intermediate strings
- **Buffered output**: Uses 4MB output buffer with incremental flushing

### Benchmark Results (4GB VCF, chr21, 427K variants, 2504 samples)

| Mode | Time |
|------|------|
| mmap (-i) | ~17s |
| stdin | ~5 min |

**Speedup: ~18x** with memory-mapped I/O

## Limitations

- Only works with biallelic variants
- Assumes diploid genotypes
- No stratification by population or other groupings
- No correction for multiple testing
- May be less accurate for extremely rare variants with very few non-reference genotypes
- No built-in filtering for variant quality or missing data threshold
