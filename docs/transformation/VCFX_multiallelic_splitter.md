# VCFX_multiallelic_splitter

## Overview

VCFX_multiallelic_splitter takes a VCF file with multi-allelic variants (variants with multiple ALT alleles) and splits them into multiple bi-allelic variant lines, while properly handling genotypes and FORMAT/INFO fields with various number specifications.

## Usage

```bash
# Standard mode (stdin)
VCFX_multiallelic_splitter [OPTIONS] < input.vcf > biallelic_output.vcf

# File mode (optimized for large files)
VCFX_multiallelic_splitter -i input.vcf > biallelic_output.vcf

# Positional argument mode
VCFX_multiallelic_splitter input.vcf > biallelic_output.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `-i`, `--input FILE` | Input VCF file (uses memory-mapped I/O for best performance) |
| `-q`, `--quiet` | Suppress warning messages |
| `-h`, `--help` | Display help message and exit |
| `-v`, `--version` | Show program version and exit |

## Description

Multi-allelic variants in VCF files (where multiple alternate alleles are specified in a comma-separated list) can complicate analysis and be incompatible with tools that require bi-allelic variants. This tool converts multi-allelic variants into the equivalent set of bi-allelic variants.

Key features:
- Maintains original VCF header information
- Correctly processes INFO fields tagged with different Number attributes (A, R, G)
- Properly adjusts genotypes and FORMAT fields for each resulting bi-allelic variant
- Preserves phasing information in genotypes
- Handles complex symbolic variants (e.g., `<DEL>`, `<INS>`)
- Correctly manages missing or malformed fields

## Output Format

The output is a standard VCF file containing:
- All header lines from the input file (unchanged)
- Bi-allelic variants only, with each multi-allelic variant split into multiple lines
- Each split variant maintains the same CHROM, POS, ID, REF, QUAL, and FILTER values
- INFO and FORMAT fields properly adjusted for each alternate allele

## Examples

### Basic Usage (stdin)

```bash
./VCFX_multiallelic_splitter < multi_allelic.vcf > biallelic.vcf
```

### File Mode (recommended for large files)

```bash
./VCFX_multiallelic_splitter -i multi_allelic.vcf > biallelic.vcf
```

### Integration with Other Tools

```bash
# Split multi-allelic variants, then run analysis requiring bi-allelic variants
cat input.vcf | \
  ./VCFX_multiallelic_splitter | \
  ./vcf_analysis_tool > results.txt
```

### Validation

```bash
# Validate that all variants in the output are indeed bi-allelic
./VCFX_multiallelic_splitter < input.vcf | \
  grep -v "^#" | awk -F'\t' '{print $5}' | grep -c ","
# Should output 0 if all variants are bi-allelic
```

## Handling Special Cases

- **INFO fields**: 
  - Number=A fields (one value per alternate allele): Each split variant gets the corresponding value
  - Number=R fields (one value per allele including reference): Values are preserved properly
  - Number=G fields (one value per genotype): Recalculated for bi-allelic case
  - Number=1 or other fixed numbers: These values are copied unchanged

- **FORMAT fields**: 
  - AD (allelic depth): Properly subset for each resulting variant
  - PL (genotype likelihoods): Recalculated for each bi-allelic output
  - GT (genotype): Adjusted to reflect the new allele indices (0/2 may become 0/1 in split variant)

- **Genotype conversion**: 
  - For each variant, genotypes are only preserved if they involve the specific alt allele
  - Genotypes not involving the current alternate allele are set to missing (./.)
  - Phased genotypes maintain their phase information

- **Edge cases**:
  - Missing data in FORMAT fields is properly handled
  - Symbolic alternate alleles are processed correctly
  - Star alleles (*) and non-ref symbolic alleles are supported

## Performance

The tool offers two execution modes with different performance characteristics:

### File Mode (-i/--input)
When using the `-i` option, the tool uses memory-mapped I/O (mmap) for optimal performance on large files. This mode provides approximately **7x speedup** compared to stdin mode.

| File Size | Variants | Samples | Stdin Mode | File Mode | Speedup |
|-----------|----------|---------|------------|-----------|---------|
| 4 GB      | 427K     | 2504    | ~2m 22s    | ~20s      | ~7x     |

### Stdin Mode (default)
For smaller files or pipeline usage, the stdin mode processes input line by line with minimal memory requirements.

### Recommendations
- For files > 100MB, use `-i input.vcf` for best performance
- For pipelines or compressed input, use stdin mode with `zcat file.vcf.gz | VCFX_multiallelic_splitter`

### Performance Characteristics

1. Single-pass processing with O(n) time complexity where n is the number of variants
2. Uses zero-copy parsing with string_view for reduced memory allocations
3. Memory-mapped I/O for efficient large file processing
4. 1MB output buffer for optimized write performance
5. Minimal memory overhead, proportional to the number of samples per line

## Limitations

- No command-line options to control the splitting behavior
- Cannot selectively split only certain multi-allelic variants
- May produce large output files when the input contains many multi-allelic variants with many samples
- Cannot reconstruct the original multi-allelic variants from the split output in all cases 
