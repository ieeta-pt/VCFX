# VCFX_ld_calculator

## Overview

`VCFX_ld_calculator` calculates pairwise linkage disequilibrium (LD) statistics between genetic variants in a VCF file, expressed as r² values. It can analyze variants across an entire file or within a specified genomic region.

## Usage

```bash
VCFX_ld_calculator [OPTIONS] < input.vcf > ld_matrix.txt
```

## Options

| Option | Description |
|--------|-------------|
| `--region <chr:start-end>` | Only compute LD for variants in the specified region |
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
4. It outputs a matrix of r² values along with variant identifiers

The r² calculation uses the standard formula:
- Let X and Y be the genotype arrays for two variants
- Calculate means of X and Y (meanX, meanY)
- Calculate covariance: cov = average(XY) - meanX * meanY
- Calculate variances: varX = average(X²) - meanX², varY similarly
- r = cov / sqrt(varX * varY)
- r² = r * r

## Output Format

The output is a tab-delimited matrix of r² values with a header identifying the variants:

```
#LD_MATRIX_START
         chr1:100 chr1:200 chr1:300
chr1:100      1.0     0.4     0.2
chr1:200     0.4      1.0     0.6
chr1:300     0.2     0.6      1.0
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

### Integration with Other Tools

Filter for common variants first, then calculate LD:

```bash
cat input.vcf | VCFX_af_subsetter --af-filter '0.05-1.0' | VCFX_ld_calculator > common_variants_ld.txt
```

## Handling Special Cases

- **Missing Genotypes**: Samples with missing genotypes (./.  or .|.) are skipped when calculating LD between variant pairs
- **Multi-allelic Variants**: Genotypes involving alleles beyond the reference and first alternate (e.g., 1/2, 2/2) are treated as missing data
- **Single Variant**: If only one variant is found in the region, the tool outputs a message stating that no pairwise LD can be calculated
- **Empty Region**: If no variants are found in the specified region, the tool outputs a message stating that no pairwise LD can be calculated
- **Invalid Region Format**: If the region format is invalid, the tool will display an error message

## Performance

- Time complexity is O(n²m) where n is the number of variants and m is the number of samples
- Memory usage scales linearly with the number of variants and samples
- For large datasets with many variants, consider using the `--region` option to limit the analysis to specific genomic regions
- The tool processes the VCF file line by line, so it can handle large files without loading the entire file into memory

## Limitations

- Only supports biallelic variants; multi-allelic variants are treated as missing data
- Requires diploid genotypes; haploid genotypes will be treated as missing data
- Assumes standard VCF format with GT field in the FORMAT column
- Does not support phased vs. unphased distinction; both "/" and "|" separators are treated the same
- No built-in visualization of LD patterns; additional tools would be needed for heatmap creation 
