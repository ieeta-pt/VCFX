# VCFX_allele_freq_calc

## Overview

The VCFX_allele_freq_calc tool calculates allele frequencies for variants in a VCF file. It reads a VCF file and outputs a TSV file with chromosome, position, ID, reference allele, alternate allele, and the calculated allele frequency.

## Usage

```bash
VCFX_allele_freq_calc [OPTIONS] [input.vcf]
VCFX_allele_freq_calc [OPTIONS] < input.vcf > allele_frequencies.tsv
```

## Options

| Option | Description |
|--------|-------------|
| `-i`, `--input FILE` | Input VCF file (uses memory-mapping for best performance) |
| `-q`, `--quiet` | Suppress informational messages |
| `-h`, `--help` | Display help message and exit |
| `-v`, `--version` | Show program version and exit |

## Description

VCFX_allele_freq_calc computes the allele frequency for each variant in a VCF file. The allele frequency is calculated as the number of alternate alleles divided by the total number of alleles (reference + alternate) across all samples, considering only non-missing genotypes.

The tool:
- Parses the GT (genotype) field for each sample
- Counts reference (0) and alternate (non-zero) alleles
- Calculates frequency as: `alternate_count / (reference_count + alternate_count)`
- Outputs results in a clean TSV format

## Output Format

The output is a tab-separated file with the following columns:

```
CHROM  POS  ID  REF  ALT  Allele_Frequency
```

Where `Allele_Frequency` is a value between 0.0 and 1.0, formatted with 4 decimal places.

## Examples

### File Input Mode (Recommended)

Use memory-mapped file I/O for best performance:

```bash
VCFX_allele_freq_calc -i input.vcf > allele_frequencies.tsv
VCFX_allele_freq_calc input.vcf > allele_frequencies.tsv
```

### Basic Usage (Stdin)

```bash
VCFX_allele_freq_calc < input.vcf > allele_frequencies.tsv
```

### Pipe with Other Commands

```bash
# Filter variants and calculate allele frequencies
grep -v "^#" input.vcf | grep "PASS" | VCFX_allele_freq_calc > filtered_allele_frequencies.tsv
```

## Handling Special Cases

- **Phased genotypes**: Both phased (`|`) and unphased (`/`) genotypes are handled the same way
- **Missing genotypes** (`./.`): Missing genotypes are skipped in the frequency calculation
- **Multiallelic sites**: All non-reference alleles are counted as "alternate" regardless of the specific ALT index
- **No GT field**: Variants without a GT field are skipped

## Performance Characteristics

The tool uses several optimizations for high performance:

- **Memory-mapped I/O**: Uses `mmap()` for file input to minimize syscall overhead
- **SIMD acceleration**: Uses NEON/SSE2/AVX2 instructions for fast line/tab scanning
- **Zero-copy parsing**: Parses VCF fields without creating intermediate strings
- **Buffered output**: Uses 4MB output buffer with direct `write()` syscalls

### Benchmark Results (4GB VCF, chr21, 427K variants, 2504 samples)

| Mode | Time |
|------|------|
| mmap (-i) | ~15s |
| stdin | ~5 min |

**Speedup: ~20x** with memory-mapped I/O

## Limitations

- Requires the GT field to be present in the FORMAT column
- Does not distinguish between different alternate alleles in multiallelic sites (all non-reference alleles are counted together)
- Cannot handle malformed VCF files, though it will attempt to skip invalid lines with a warning
