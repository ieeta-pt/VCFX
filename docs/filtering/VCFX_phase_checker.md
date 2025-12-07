# VCFX_phase_checker

## Overview

The VCFX_phase_checker tool filters a VCF file to retain only those variant lines where all sample genotypes are fully phased (using the pipe '|' phasing separator). This is particularly useful for downstream analyses that require complete phasing information.

## Usage

```bash
VCFX_phase_checker [OPTIONS] [input.vcf]
VCFX_phase_checker [OPTIONS] < input.vcf > phased_output.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `-h`, `--help` | Display help message and exit |
| `-v`, `--version` | Show program version and exit |
| `-i`, `--input FILE` | Input VCF file (uses fast memory-mapped I/O) |
| `-q`, `--quiet` | Suppress warning messages to stderr |

## Description

VCFX_phase_checker reads a VCF file from standard input or a file and examines the GT (genotype) field for every sample in each variant line. It determines whether each genotype is fully phased using the following criteria:

- A genotype is considered fully phased only if:
  - It uses the pipe character '|' as the separator between alleles (e.g., "0|1")
  - It contains no missing alleles (no ".")
  - It contains no unphased separators ("/")

The tool outputs only those variant lines where all sample genotypes meet these criteria. Lines that don't meet these criteria are skipped, and warnings are written to standard error (unless `-q` is specified).

## Output

The output is a valid VCF file containing:
- All header lines from the input file (unchanged)
- Only the variant lines where all samples have fully phased genotypes
- Warnings (to stderr) for each line that was skipped due to unphased genotypes (unless `-q` is used)

## Examples

### Basic Usage (stdin)

```bash
./VCFX_phase_checker < input.vcf > phased_output.vcf
```

### File Input Mode (Fast)

```bash
# Using -i flag (recommended for large files)
./VCFX_phase_checker -i input.vcf > phased_output.vcf

# Using positional argument
./VCFX_phase_checker input.vcf > phased_output.vcf
```

### Quiet Mode (Suppress Warnings)

```bash
# Useful for batch processing when you don't need warnings
./VCFX_phase_checker -q -i input.vcf > phased_output.vcf
```

### Capturing Warnings

```bash
./VCFX_phase_checker < input.vcf > phased_output.vcf 2> unphased_warnings.log
```

### Counting Phased vs. Unphased Variants

```bash
# Count total variants
total=$(grep -v "^#" input.vcf | wc -l)
# Count phased variants
phased=$(./VCFX_phase_checker -q -i input.vcf | grep -v "^#" | wc -l)
# Calculate percentage
echo "Phased variants: $phased / $total ($(echo "scale=2; 100*$phased/$total" | bc)%)"
```

## Handling Special Cases

- **Haploid genotypes** (e.g., "0" or "1"): These are not considered phased; the line will be skipped
- **Missing genotypes** (e.g., "./." or ".|."): These are not considered phased; the line will be skipped
- **Missing GT field**: Lines without a GT field in the FORMAT column are skipped with a warning
- **Multiallelic variants**: These are treated the same as biallelic variants, as long as all alleles are phased
- **Non-VCF-compliant genotype notation**: Any genotype that doesn't follow standard VCF format is not considered phased
- **Header lines**: All header lines (starting with "#") are preserved in the output
- **Samples with different ploidy levels**: Each sample is checked independently; if all are phased, the line is kept

## Performance

The tool offers two processing modes with significantly different performance characteristics:

### File Input Mode (Recommended for Large Files)

```bash
VCFX_phase_checker -i input.vcf > phased.vcf
VCFX_phase_checker input.vcf > phased.vcf
```

File input uses memory-mapped I/O with the following optimizations:
- **SIMD-optimized line scanning**: Uses AVX2 (32 bytes/cycle) or SSE2 (16 bytes/cycle) on x86_64
- **Zero-copy parsing**: Uses `string_view` to avoid string allocations
- **1MB output buffering**: Reduces write syscalls by 100-1000x
- **FORMAT caching**: GT index computed once per unique FORMAT string
- **Early termination**: Stops checking samples on first unphased genotype

Expected speedup: **10-12x** compared to stdin mode.

### Stdin Mode (Fallback)

```bash
VCFX_phase_checker < input.vcf > phased.vcf
cat input.vcf | VCFX_phase_checker > phased.vcf
```

Stdin mode is optimized but inherently slower due to streaming constraints. Still benefits from zero-copy parsing and early termination.

### Benchmarks

| Variants | Samples | stdin Mode | mmap Mode | Speedup |
|----------|---------|------------|-----------|---------|
| 1,000    | 3       | ~25ms      | ~23ms     | ~1.1x   |
| 100,000  | 3       | ~430ms     | ~36ms     | ~12x    |
| 1,000,000| 3       | ~4.3s      | ~360ms    | ~12x    |

*Benchmarks on Apple M1/Intel x86_64, actual performance may vary. Speedup increases with file size due to mmap overhead amortization.*

## Limitations

- No option to make a best-effort phasing assumption
- Cannot output partially phased lines or filter specific samples
- Designed to be used as part of a pipeline, not as a standalone phasing tool
