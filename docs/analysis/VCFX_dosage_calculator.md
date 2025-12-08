# VCFX_dosage_calculator

## Overview

VCFX_dosage_calculator computes the genetic dosage (count of alternate alleles) for each sample at each variant position in a VCF file, outputting a tab-delimited summary in a convenient format for downstream analysis.

## Usage

```bash
# Standard mode (stdin)
VCFX_dosage_calculator [OPTIONS] < input.vcf > dosage_output.txt

# File mode (optimized for large files)
VCFX_dosage_calculator -i input.vcf > dosage_output.txt

# Positional argument mode
VCFX_dosage_calculator input.vcf > dosage_output.txt
```

## Options

| Option | Description |
|--------|-------------|
| `-i`, `--input FILE` | Input VCF file (uses memory-mapped I/O for best performance on large files) |
| `-q`, `--quiet` | Suppress warning messages |
| `-h`, `--help` | Display help message and exit |
| `-v`, `--version` | Show program version and exit |

## Description

VCFX_dosage_calculator reads a VCF file from standard input and calculates the genetic dosage for each sample at each variant position. Dosage is defined as the number of alternate alleles in a genotype, regardless of the specific alternate allele identifier.

The dosage calculation follows these rules:
- Reference allele (0): Counts as 0
- Any alternate allele (1, 2, 3, etc.): Counts as 1
- Missing or invalid alleles: Reported as "NA"

For diploid genotypes, this results in:
- 0/0 → dosage 0 (no alternate alleles)
- 0/1 → dosage 1 (one alternate allele)
- 1/1 → dosage 2 (two alternate alleles)
- 1/2 → dosage 2 (two alternate alleles, even though they're different alternates)
- ./. → dosage NA (missing genotype)

The tool is useful for:
- Preparing genetic data for association tests
- Quantifying genetic burden
- Converting genotypes to a simple numeric format for statistical analyses
- Creating matrix-like representations of genetic data

## Output Format

The output is a tab-delimited text file with the following columns:

```
CHROM  POS  ID  REF  ALT  Dosages
```

Where `Dosages` is a comma-separated list of dosage values for each sample, in the same order as they appear in the VCF. Values can be 0, 1, 2 (for diploid organisms), or NA for missing/invalid data.

## Examples

### Basic Usage (stdin)

```bash
./VCFX_dosage_calculator < input.vcf > dosage_results.txt
```

### File Mode (recommended for large files)

```bash
./VCFX_dosage_calculator -i input.vcf > dosage_results.txt
```

### Quiet Mode (suppress warnings)

```bash
./VCFX_dosage_calculator -q -i input.vcf > dosage_results.txt
```

### Viewing Results

```bash
# Show the first few lines of results with column headers
head -n 5 dosage_results.txt | column -t
```

### Integration with Other Tools

```bash
# Calculate dosage and use for statistical analysis
./VCFX_dosage_calculator < input.vcf | \
  awk -F'\t' '{split($6,d,","); sum=0; count=0; for(i in d) if(d[i]!="NA") {sum+=d[i]; count++} if(count>0) print $1,$2,$3,sum/count}' > avg_dosage.txt
```

## Handling Special Cases

- **Phased genotypes**: Phasing information is ignored; "0|1" and "0/1" both result in dosage 1
- **Missing genotypes**: Genotypes represented as "./.", ".|.", or "." are reported as NA
- **Multi-allelic variants**: All non-reference alleles count equally:
  - 1/2 → dosage 2 
  - 0/3 → dosage 1
  - 2/3 → dosage 2
- **Malformed genotypes**: Any genotype that doesn't follow the expected format (e.g., "0/X", "ABC") is reported as NA
- **Non-diploid organisms**: The tool assumes diploid genotypes; other ploidy levels may produce unexpected results
- **No GT field**: If the FORMAT column doesn't include GT, a warning is issued and NA is reported for all samples
- **Missing header**: If the #CHROM header line isn't found, processing stops with an error

## Performance

The tool offers two execution modes with different performance characteristics:

### File Mode (-i/--input)
When using the `-i` option, the tool uses memory-mapped I/O (mmap) with SIMD-accelerated newline detection for optimal performance on large files. This mode provides approximately **8-10x speedup** compared to stdin mode.

| File Size | Variants | Samples | Stdin Mode | File Mode | Speedup |
|-----------|----------|---------|------------|-----------|---------|
| 4 GB      | 427K     | 2504    | ~2m 41s    | ~18s      | ~9x     |

### Stdin Mode (default)
For smaller files or pipeline usage, the stdin mode processes input line by line with minimal memory requirements. Performance is primarily dependent on the number of samples in the VCF file, as each sample's genotype must be processed for every variant.

### Recommendations
- For files > 100MB, use `-i input.vcf` for best performance
- For pipelines or compressed input, use stdin mode with `zcat file.vcf.gz | VCFX_dosage_calculator`

## Limitations

- Designed primarily for diploid organisms; may not be suitable for polyploid data
- Cannot handle probabilities or fractional dosages (e.g., from imputed data)
- Treats all alternate alleles equally; cannot distinguish between different alternate alleles
- No special handling for sex chromosomes (X, Y) which may have different ploidy in males and females
- No filtering options within the tool (use other VCFX tools for pre-filtering) 
