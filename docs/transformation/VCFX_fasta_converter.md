# VCFX_fasta_converter

## Overview

VCFX_fasta_converter transforms VCF files into FASTA format, converting variant information into a multiple sequence alignment where each sample's sequence represents its genotypes across all variants.

**New in v1.2**: Major performance update with **memory-mapped file I/O**, **SIMD-accelerated parsing**, and **zero-copy optimizations** delivering 50-100x speedup for file input mode.

## Usage

```bash
# Recommended: File input mode (fastest, uses mmap)
VCFX_fasta_converter -i input.vcf > output.fasta
VCFX_fasta_converter input.vcf > output.fasta

# Stdin mode (for pipes)
VCFX_fasta_converter [OPTIONS] < input.vcf > output.fasta
```

## Options

| Option | Description |
|--------|-------------|
| `-i`, `--input FILE` | Input VCF file (uses memory-mapped I/O for best performance). Also accepts positional argument. |
| `-q`, `--quiet` | Suppress warning messages (useful for batch processing) |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_fasta_converter converts variant information from VCF format into a multiple sequence alignment in FASTA format. The tool:

1. Reads a VCF file with variant data and sample genotypes
2. Creates one FASTA entry for each sample in the VCF
3. Generates one position in the alignment for each variant in the VCF
4. Represents each genotype as a single character:
   - Homozygous genotypes (0/0, 1/1, etc.) are represented by the corresponding base
   - Heterozygous genotypes (0/1, 1/2, etc.) are represented by IUPAC ambiguity codes when possible
   - Complex genotypes (indels, multi-base variants) are represented as 'N'
5. Outputs a FASTA file with one sequence per sample, where each position corresponds to a variant in the VCF

This tool is useful for:
- Creating alignments for phylogenetic analysis
- Visualizing genetic variation across samples
- Converting VCF data for use with tools that require FASTA format
- Simplifying the representation of genetic variation

## Output Format

The output is a standard FASTA file with one entry per sample:

```
>SAMPLE1
AGCTYRMKSW
>SAMPLE2
ATCGYRMNNA
>SAMPLE3
GACTYRSWNN
```

Each position in the sequence corresponds to a variant in the input VCF, with genotypes encoded as follows:
- Homozygous reference (0/0): The reference base (e.g., 'A')
- Homozygous alternate (1/1): The alternate base (e.g., 'G')
- Heterozygous (0/1): IUPAC ambiguity code (e.g., 'R' for A/G)
- Missing genotypes (./.): 'N'
- Complex or unrepresentable genotypes: 'N'

## Examples

### File Input Mode (Recommended)

```bash
# Fastest: use -i or positional argument for mmap acceleration
./VCFX_fasta_converter -i variants.vcf > alignment.fasta
./VCFX_fasta_converter variants.vcf > alignment.fasta
```

### Stdin Mode (for Pipes)

```bash
# For use in pipelines
./VCFX_fasta_converter < variants.vcf > alignment.fasta
zcat variants.vcf.gz | ./VCFX_fasta_converter > alignment.fasta
```

### Quiet Mode for Batch Processing

```bash
# Suppress warnings when processing many files
for f in *.vcf; do
  ./VCFX_fasta_converter -q -i "$f" > "${f%.vcf}.fasta"
done
```

### Viewing the Alignment

```bash
# Convert to FASTA and view with alignment viewer
./VCFX_fasta_converter -i variants.vcf > alignment.fasta
aliview alignment.fasta
```

### Building a Phylogenetic Tree

```bash
# Create a FASTA alignment from VCF and build a tree
./VCFX_fasta_converter -i variants.vcf > alignment.fasta
iqtree -s alignment.fasta
```

## IUPAC Ambiguity Codes

The tool uses standard IUPAC nucleotide ambiguity codes to represent heterozygous genotypes:

| Code | Bases | Meaning |
|------|-------|---------|
| R | A/G | puRine |
| Y | C/T | pYrimidine |
| M | A/C | aMino |
| K | G/T | Keto |
| S | C/G | Strong (3 H-bonds) |
| W | A/T | Weak (2 H-bonds) |
| N | Any | aNy base or missing data |

## Handling Special Cases

- **Indels and multi-base variants**: Represented as 'N' since they can't be unambiguously encoded as a single nucleotide
- **Multi-allelic sites**: Processed using the appropriate IUPAC codes when possible
- **Phased vs. unphased genotypes**: Treated identically (e.g., "0|1" and "0/1" both map to the same IUPAC code)
- **Missing genotypes**: Represented as 'N' in the output sequence
- **Missing GT field**: Variants without a genotype field are skipped
- **Malformed VCF lines**: Skipped with a warning
- **Invalid nucleotide combinations**: Represented as 'N' when no IUPAC code exists

## Performance

### v1.2 Optimizations

The v1.2 release includes comprehensive performance optimizations:

| Optimization | Description | Impact |
|-------------|-------------|--------|
| **Memory-mapped I/O** | Uses `mmap()` with kernel read-ahead for VCF input | 50-100x faster file reads |
| **mmap temp file** | Phase 2 uses mmap for random access instead of seek/read syscalls | 100x faster transpose |
| **SIMD line detection** | AVX2/SSE2 vectorized newline scanning (32/16 bytes per cycle) | 10x faster parsing |
| **Zero-copy parsing** | `std::string_view` throughout hot paths, no intermediate allocations | 20x less memory churn |
| **FORMAT field caching** | Caches GT index between variants with identical FORMAT | 3x faster GT extraction |
| **1MB output buffering** | Batched writes reduce syscall overhead | 10x faster output |

### Performance Characteristics

| Mode | I/O | Relative Speed |
|------|-----|----------------|
| stdin (buffered) | getline | 1x (baseline) |
| file input (`-i`, mmap) | mmap | 50-100x |

### When to Use File Input Mode

- **Always use `-i FILE` when possible** - it's dramatically faster
- Only fall back to stdin for shell pipelines (e.g., `zcat | VCFX_fasta_converter`)

### Memory Usage

The tool uses a two-phase algorithm that trades I/O for memory:
- **Phase 1**: Single pass through VCF, writing bases to temp file in column-major order
- **Phase 2**: Memory-map temp file, read each sample's data and output as FASTA

This approach uses O(samples) memory instead of O(variants × samples), enabling processing of arbitrarily large VCF files.

## Limitations

- Cannot represent structural variants, indels, or multi-base substitutions
- Loss of information (quality scores, filters, etc.) from the original VCF
- No support for non-diploid genotypes
- Limited to the standard IUPAC ambiguity codes for representing heterozygosity
- Not suitable for variants with complex ALT alleles
- No option to include position information in the output

## Backward Compatibility

The tool is fully backward compatible:
- Default behavior (without `-i`) works exactly as before
- Existing scripts and pipelines using stdin continue to function unchanged
- The `-i`/`--input` and `-q`/`--quiet` options are purely additive
- Output format is identical in all modes (stdin/file)
