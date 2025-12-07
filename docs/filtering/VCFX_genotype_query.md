# VCFX_genotype_query

## Overview

VCFX_genotype_query filters a VCF file to retain only variant lines where at least one sample has a specified genotype. Supports both flexible matching (normalizes phasing and allele order) and strict exact matching.

## Usage

```bash
VCFX_genotype_query [OPTIONS] [input.vcf]
VCFX_genotype_query [OPTIONS] < input.vcf > output.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `-g`, `--genotype-query GT` | Genotype to query (e.g., "0/1", "1\|1") |
| `-i`, `--input FILE` | Input VCF file (uses fast memory-mapped I/O) |
| `--strict` | Exact string matching (no normalization) |
| `-q`, `--quiet` | Suppress warning messages to stderr |
| `-h`, `--help` | Display help message and exit |
| `-v`, `--version` | Show program version and exit |

## Description

VCFX_genotype_query reads a VCF file and examines the GT (genotype) field for every sample in each variant line. It outputs only those variant lines where at least one sample matches the specified genotype.

### Flexible Matching (Default)

By default, genotypes are normalized before comparison:
- Phasing separators are unified: `0|1` and `0/1` both match `0/1`
- Allele order is normalized: `1/0` matches `0/1`

This means querying for `0/1` will match:
- `0/1` (exact match)
- `1/0` (reversed alleles)
- `0|1` (phased equivalent)
- `1|0` (phased and reversed)

### Strict Matching

With `--strict`, genotypes must match exactly as specified:
- Querying `0|1` will NOT match `0/1`
- Querying `0/1` will NOT match `1/0`

## Output

The output is a valid VCF file containing:
- All header lines from the input file (unchanged)
- Only the variant lines where at least one sample has the matching genotype
- Warnings (to stderr) for malformed lines (unless `-q` is used)

## Examples

### Basic Usage (stdin)

```bash
# Find all variants with heterozygous samples
./VCFX_genotype_query -g "0/1" < input.vcf > het_variants.vcf
```

### File Input Mode (Fast)

```bash
# Using -i flag (recommended for large files)
./VCFX_genotype_query -g "0/1" -i input.vcf > het_variants.vcf

# Using positional argument
./VCFX_genotype_query -g "0/1" input.vcf > het_variants.vcf
```

### Strict Matching

```bash
# Only match phased heterozygous (0|1, not 0/1)
./VCFX_genotype_query -g "0|1" --strict < input.vcf > phased_het.vcf
```

### Quiet Mode (Suppress Warnings)

```bash
# Useful for batch processing
./VCFX_genotype_query -g "0/1" -q -i input.vcf > output.vcf
```

### Find Homozygous Alternate Samples

```bash
./VCFX_genotype_query -g "1/1" -i input.vcf > hom_alt.vcf
```

### Pipeline Usage

```bash
# Filter by genotype, then by quality
./VCFX_genotype_query -g "0/1" -i input.vcf | ./VCFX_quality_filter -q 30 > filtered.vcf
```

## Handling Special Cases

- **Missing genotypes** (e.g., "./.") are not matched
- **Haploid genotypes** (e.g., "0" or "1") are not matched in flexible mode
- **Multi-allelic variants** are supported (e.g., "1/2", "2/2")
- **Polyploid genotypes** are supported but must match the same ploidy as the query
- **Missing GT field**: Lines without a GT field in the FORMAT column are skipped
- **Header lines**: All header lines (starting with "#") are preserved in the output

## Performance

The tool offers two processing modes with significantly different performance characteristics:

### File Input Mode (Recommended for Large Files)

```bash
VCFX_genotype_query -g "0/1" -i input.vcf > output.vcf
```

File input uses memory-mapped I/O with the following optimizations:
- **SIMD-optimized line scanning**: Uses AVX2 (32 bytes/cycle) or SSE2 (16 bytes/cycle) on x86_64
- **Zero-copy parsing**: Uses `string_view` to avoid string allocations
- **1MB output buffering**: Reduces write syscalls by 100-1000x
- **FORMAT caching**: GT index computed once per unique FORMAT string
- **Direct genotype matching**: Fast path for common 3-char diploid genotypes
- **Early termination**: Stops checking samples on first match

Expected speedup: **40-50x** compared to stdin mode.

### Stdin Mode (Fallback)

```bash
VCFX_genotype_query -g "0/1" < input.vcf > output.vcf
cat input.vcf | VCFX_genotype_query -g "0/1" > output.vcf
```

Stdin mode is optimized but inherently slower due to streaming constraints. Still benefits from FORMAT caching and early termination.

### Benchmarks

| Variants | Samples | stdin Mode | mmap Mode | Speedup |
|----------|---------|------------|-----------|---------|
| 1,000    | 3       | ~30ms      | ~25ms     | ~1.2x   |
| 100,000  | 3       | ~500ms     | ~30ms     | ~17x    |
| 100,000  | 2,504   | ~64s       | ~1.3s     | ~50x    |
| 1,000,000| 3       | ~5s        | ~300ms    | ~17x    |

*Benchmarks on Apple M1/Intel x86_64, actual performance may vary. Speedup increases significantly with file size and sample count.*

## Limitations

- Only checks the GT (genotype) field
- Cannot filter by other FORMAT fields (DP, GQ, etc.)
- Cannot use complex expressions (e.g., "0/1 OR 0/2")
- Designed to be used as part of a pipeline, not as a standalone analysis tool
