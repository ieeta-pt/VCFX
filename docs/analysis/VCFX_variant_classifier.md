# VCFX_variant_classifier

## Overview

The VCFX_variant_classifier tool analyzes VCF files and classifies variants into various types: SNP, INDEL, MNV, or STRUCTURAL. It can either produce a TSV summary or append classifications to the original VCF file.

## Usage

```bash
VCFX_variant_classifier [OPTIONS] < input.vcf > output.vcf_or_tsv
```

## Options

| Option | Description |
|--------|-------------|
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |
| `-a`, `--append-info` | Instead of producing a TSV, output a valid VCF with a new 'VCF_CLASS' subfield in the INFO column |

## Description

VCFX_variant_classifier reads each variant line from a VCF file and determines its type based on the following criteria:

- **SNP**: Single nucleotide polymorphism, where both reference and alternate alleles are single bases
- **INDEL**: Insertions or deletions with length difference less than 50 bp
- **MNV**: Multi-nucleotide variants with the same length but multiple bases changed
- **STRUCTURAL**: Complex variants including:
  - Symbolic ALT fields (`<DEL>`, `<INS>`, etc.)
  - Breakend notation (containing `[` or `]`)
  - Variants with length difference ≥50 bp
  - Very large reference or alternate alleles (≥40 bp)
- **UNKNOWN**: Reserved for special cases like missing or identical REF/ALT

## Output Formats

### TSV Mode (Default)

By default, the tool outputs a TSV file with the following columns:
```
CHROM  POS  ID  REF  ALT  Classification
```

### VCF Mode (with --append-info)

When using the `--append-info` option, the tool:
- Preserves the original VCF format including all headers
- Adds a `VCF_CLASS=TYPE` entry to the INFO field of each variant
- Maintains all other VCF fields

## Examples

### Basic Classification to TSV

```bash
./VCFX_variant_classifier < input.vcf > classified.tsv
```

### Append Classification to VCF

```bash
./VCFX_variant_classifier --append-info < input.vcf > annotated.vcf
```

### Filtering Based on Classification

```bash
# First classify, then filter for structural variants only
./VCFX_variant_classifier < input.vcf | grep "STRUCTURAL" > structural_variants.tsv
```

## Handling Special Cases

- **Multi-allelic sites**: The most complex type among all alternates is assigned (STRUCTURAL > MNV > INDEL > SNP)
- **Malformed lines**: Lines with fewer than 8 columns are skipped with a warning
- **Missing data**: Missing ALT fields or identical REF/ALT entries are classified as UNKNOWN
- **Symbolic alleles**: Any variant with symbolic notation (e.g., `<DEL>`) is classified as STRUCTURAL

## Performance

The tool efficiently processes VCF files line by line, allowing it to handle very large files with minimal memory requirements.

## Limitations

- Classification is based on standard VCF conventions and may need adjustment for non-standard VCFs
- Cannot detect complex structural variants that aren't properly annotated in the VCF
- Edge cases like very long identical stretches may be classified in unexpected ways 
