# VCFX_indel_normalizer

## Overview
`VCFX_indel_normalizer` is a tool for normalizing indel variants in VCF files. It performs left-alignment of variants by removing common prefixes and suffixes between reference and alternate alleles, and splits multi-allelic variants into separate records. This normalization is done without requiring an external reference genome.

## Usage
```bash
VCFX_indel_normalizer [OPTIONS] [input.vcf]
VCFX_indel_normalizer [OPTIONS] < input.vcf > normalized.vcf
```

## Options
| Option | Description |
|--------|-------------|
| `-i`, `--input FILE` | Input VCF file (uses memory-mapping for best performance) |
| `-q`, `--quiet` | Suppress informational messages |
| `-h`, `--help` | Display help message and exit |
| `-v`, `--version` | Show program version and exit |

## Description
`VCFX_indel_normalizer` processes a VCF file and normalizes indel variants by:

1. Reading the VCF file using memory-mapped I/O for extreme performance
2. Preserving all header lines without modification
3. For variants with multiple alternate alleles (comma-separated ALT values):
   - Splitting them into separate lines, one per alternate allele
4. For each variant:
   - Removing the longest common prefix from REF and ALT, keeping at least one base
   - Removing the longest common suffix from REF and ALT, keeping at least one base
   - Adjusting the position (POS) to account for removed leading bases
5. Writing the normalized variants using buffered output

This normalization ensures that variants are represented in a consistent, minimal left-aligned form, which is important for variant comparison, annotation, and analysis.

## Normalization Process

### Left Alignment Algorithm
The tool implements a simple but effective left-alignment approach:

1. **Prefix Removal**:
   - Identify the longest common prefix between REF and ALT
   - Remove all but one base of this common prefix
   - Adjust the variant position to account for the removed bases

2. **Suffix Removal**:
   - Identify the longest common suffix between REF and ALT
   - Remove all but one base of this common suffix

3. **Special Case Handling**:
   - If after normalization REF or ALT is empty, the variant is considered invalid
   - If after normalization REF and ALT are identical, the variant is considered invalid

### Multi-allelic Variant Handling
For variants with multiple alternate alleles:
- Each alternate allele is processed separately
- A new VCF line is generated for each alternate allele
- Each new line undergoes the same normalization process

## Examples

### File Input Mode (Recommended)

Use memory-mapped file I/O for best performance:

```bash
VCFX_indel_normalizer -i input.vcf > normalized.vcf
VCFX_indel_normalizer input.vcf > normalized.vcf
```

### Basic Usage (Stdin)

Normalize indels in a VCF file:
```bash
VCFX_indel_normalizer < input.vcf > normalized.vcf
```

### Quiet Mode for Scripts

```bash
VCFX_indel_normalizer -q -i input.vcf > normalized.vcf
```

### Example Transformations

#### Prefix Removal
```
Before: chr1 100 . ACTG AC   (deletion of TG)
After:  chr1 102 . TG   -    (adjusted position, simplified representation)
```

#### Suffix Removal
```
Before: chr1 100 . ACTG ACCC (substitution of TG with CC)
After:  chr1 100 . AC   ACC  (common suffix 'C' removed)
```

#### Multi-allelic Splitting
```
Before: chr1 100 . A    C,G,T
After:  chr1 100 . A    C
        chr1 100 . A    G
        chr1 100 . A    T
```

## Handling Special Cases

### Invalid Variants
- If after normalization a variant would have an empty REF or ALT, the original line is preserved
- If after normalization REF equals ALT (no actual variant), the original line is preserved

### Empty Lines
- Empty lines in the input are preserved as empty lines in the output

### Malformed Lines
- Lines with fewer than 10 columns (minimum for a VCF with samples) are output unchanged
- Lines with invalid position values are output unchanged

## Performance Characteristics

The tool uses several optimizations for extreme performance:

- **Memory-mapped I/O**: Uses `mmap()` for file input to minimize syscall overhead
- **SIMD acceleration**: Uses NEON/SSE2/AVX2 instructions for fast line/tab scanning
- **Zero-copy parsing**: Parses VCF fields without creating intermediate strings
- **Buffered output**: Uses 4MB output buffer with direct `write()` syscalls

### Benchmark Results (4GB VCF, chr21, 427K variants, 2504 samples)

| Mode | Time |
|------|------|
| mmap (-i) | ~4s |
| stdin | ~4.9 min |

**Speedup: ~73x** with memory-mapped I/O

## Limitations
- This tool performs simple left-alignment without checking for sequence repeats in a reference genome
- Full left-alignment of variants in repetitive regions requires a reference sequence
- Cannot handle complex structural variants beyond simple indels
- Limited to the information available in the VCF file itself
- No automatic breakup of complex variants (substitutions that could be represented as indels)
- No variant filtering capabilities
