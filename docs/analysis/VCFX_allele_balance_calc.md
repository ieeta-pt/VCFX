# VCFX_allele_balance_calc

## Overview
`VCFX_allele_balance_calc` calculates the allele balance for each sample in a VCF file, which is the ratio of reference alleles to alternate alleles in heterozygous genotypes. This metric is useful for assessing potential allelic bias in sequencing data.

## Usage
```bash
VCFX_allele_balance_calc [OPTIONS] [FILE]
```

## Options
| Option | Description |
|--------|-------------|
| `-i`, `--input FILE` | Input VCF file (uses memory-mapping for best performance) |
| `-t`, `--threads N` | Number of threads for parallel processing (default: auto-detect CPU cores) |
| `-s`, `--samples "Sample1 Sample2..."` | Optional. Specify sample names to calculate allele balance for (space-separated). If omitted, all samples are processed. |
| `-q`, `--quiet` | Suppress informational messages |
| `-h`, `--help` | Display help message and exit |
| `-v`, `--version` | Show program version and exit |

## Description
`VCFX_allele_balance_calc` processes a VCF file and calculates the allele balance for each variant in each specified sample. The tool:

1. Reads a VCF file from standard input or file
2. Identifies sample columns from the VCF header
3. For each variant and each sample:
   - Extracts the genotype information
   - Counts reference (0) and alternate (non-0) alleles
   - Calculates the allele balance as: reference allele count / alternate allele count
4. Outputs a tab-separated file with allele balance values for each variant-sample combination

This tool is particularly useful for:
- Identifying potential allele-specific biases in sequencing
- Quality control of variant calls
- Assessing imbalanced expression of alleles
- Detecting potential sample contamination

## Output Format
The tool produces a tab-separated values (TSV) file with the following columns:

| Column | Description |
|--------|-------------|
| CHROM | Chromosome of the variant |
| POS | Position of the variant |
| ID | Variant identifier |
| REF | Reference allele |
| ALT | Alternate allele(s) |
| Sample | Sample name |
| Allele_Balance | Calculated allele balance value or "NA" for missing/invalid genotypes |

## Examples

### File Input Mode (Recommended)
Use memory-mapped file I/O for best performance with multi-threaded processing:
```bash
VCFX_allele_balance_calc -i input.vcf > allele_balance.tsv
VCFX_allele_balance_calc -t 8 -i input.vcf > allele_balance.tsv  # Use 8 threads
VCFX_allele_balance_calc input.vcf > allele_balance.tsv
```

### Basic Usage (Stdin)
Calculate allele balance for all samples in a VCF file:
```bash
VCFX_allele_balance_calc < input.vcf > allele_balance_all.tsv
```

### Specific Samples
Calculate allele balance for specific samples:
```bash
VCFX_allele_balance_calc --samples "SAMPLE1 SAMPLE2" -i input.vcf > allele_balance_subset.tsv
```

### Filtering Results
Process the output to focus on imbalanced variants:
```bash
VCFX_allele_balance_calc -i input.vcf | awk -F'\t' '$7 != "NA" && $7 < 0.4' > imbalanced_variants.tsv
```

## Allele Balance Calculation

### Formula
The allele balance is calculated as:
```
Allele Balance = Number of Reference Alleles / Number of Alternate Alleles
```

Where:
- Reference alleles are those with value "0" in the genotype field
- Alternate alleles are any non-zero value (e.g., "1", "2", etc.) in the genotype field

### Interpretation
- Value of 0.0: No reference alleles (e.g., "1/1", "1/2")
- Value of 1.0: Equal number of reference and alternate alleles (e.g., "0/1")
- Value > 1.0: More reference than alternate alleles (unusual for diploid organisms)
- "NA": Missing or invalid genotype

### Special Cases
- Homozygous reference (e.g., "0/0"): Returns 0.0 (technically it would be undefined due to division by zero)
- Missing genotypes (e.g., "./.", ".|."): Returns "NA" in the output
- Partial missing (e.g., "0/."): Only valid alleles are counted
- Invalid formats: Returns "NA" in the output

## Handling Special Cases

### Missing Data
- Genotypes with missing values (`./.`, `.`) return "NA" for allele balance
- Partial missing genotypes only count the valid alleles present

### Multi-allelic Sites
- All non-reference alleles are treated as "alternate" regardless of their specific number
- For example, in a genotype "1/2", both alleles are counted as alternate alleles

### Phased Genotypes
- Phasing information is ignored for allele balance calculation
- Phased genotypes (e.g., "0|1") are treated the same as unphased (e.g., "0/1")

### Haploid Genotypes
- Not explicitly handled; the tool expects diploid or polyploid genotypes with separators

## Performance Characteristics

The tool uses several optimizations for high performance:

- **Memory-mapped I/O**: Uses `mmap()` for file input to minimize syscall overhead
- **SIMD acceleration**: Uses NEON/SSE2/AVX2 instructions for fast line/tab scanning
- **Batch processing**: Pre-computes all sample positions per line for O(1) access
- **Zero-copy parsing**: Parses VCF fields without creating intermediate strings
- **Incremental flushing**: 4MB buffered output with automatic flushing to prevent memory buildup

### Benchmark Results (4GB VCF, chr21, 427K variants)

| Samples | Time |
|---------|------|
| 1 | ~4s |
| 100 | ~8s |
| 2504 (all) | ~41s |

**Speedup vs previous version: ~50x** (35 min → 41 sec)

### Throughput

On typical hardware with multi-threading:
- **VCF parsing**: ~1.5 GB/s with mmap
- **Output generation**: ~10-15M lines/sec (with 100+ samples)
- **Scaling**: Near-linear scaling up to ~4-6 cores

## Limitations
- No option to customize the allele balance formula
- Simplified handling of multi-allelic sites (all non-reference alleles are grouped)
- No automatic filtering based on allele balance values
- Cannot account for read depth or genotype quality in calculations
- Limited to processing standard VCF genotype fields
- Does not produce summary statistics across all variants
