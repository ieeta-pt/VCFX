# VCFX_multiallelic_splitter

## Overview

VCFX_multiallelic_splitter takes a VCF file with multi-allelic variants (variants with multiple ALT alleles) and splits them into multiple bi-allelic variant lines, while properly handling genotypes and FORMAT/INFO fields with various number specifications.

## Usage

```bash
VCFX_multiallelic_splitter [OPTIONS] < input.vcf > biallelic_output.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `--help`, `-h` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

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

### Basic Usage

```bash
./VCFX_multiallelic_splitter < multi_allelic.vcf > biallelic.vcf
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

The tool processes VCF files line by line with minimal memory requirements, with performance primarily dependent on:
- Number of samples in the VCF
- Number of multi-allelic sites
- Complexity of INFO and FORMAT fields

For very large VCF files with many samples, processing time scales linearly with file size.

## Limitations

- No command-line options to control the splitting behavior
- Cannot selectively split only certain multi-allelic variants
- May produce large output files when the input contains many multi-allelic variants with many samples
- Cannot reconstruct the original multi-allelic variants from the split output in all cases 
