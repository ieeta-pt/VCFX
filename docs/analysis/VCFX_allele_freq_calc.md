# VCFX_allele_freq_calc

## Overview

The VCFX_allele_freq_calc tool calculates allele frequencies for variants in a VCF file. It reads a VCF file from standard input and outputs a TSV file with chromosome, position, ID, reference allele, alternate allele, and the calculated allele frequency.

## Usage

```bash
VCFX_allele_freq_calc [OPTIONS] < input.vcf > allele_frequencies.tsv
```

## Options

| Option      | Description                                |
|-------------|--------------------------------------------|
| `--help`, `-h` | Display help message and exit (handled by `vcfx::handle_common_flags`)              |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

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

### Basic Usage

```bash
./VCFX_allele_freq_calc < input.vcf > allele_frequencies.tsv
```

### Pipe with Other Commands

```bash
# Filter variants and calculate allele frequencies
grep -v "^#" input.vcf | grep "PASS" | ./VCFX_allele_freq_calc > filtered_allele_frequencies.tsv
```

## Handling Special Cases

- **Phased genotypes**: Both phased (`|`) and unphased (`/`) genotypes are handled the same way
- **Missing genotypes** (`./.`): Missing genotypes are skipped in the frequency calculation
- **Multiallelic sites**: All non-reference alleles are counted as "alternate" regardless of the specific ALT index
- **No GT field**: Variants without a GT field are skipped

## Performance

This tool processes VCF files line by line, with minimal memory requirements. It can handle large VCF files efficiently.

## Limitations

- Requires the GT field to be present in the FORMAT column
- Does not distinguish between different alternate alleles in multiallelic sites (all non-reference alleles are counted together)
- Cannot handle malformed VCF files, though it will attempt to skip invalid lines with a warning 
