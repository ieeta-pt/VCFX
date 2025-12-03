# VCFX_ref_comparator

## Overview

VCFX_ref_comparator validates VCF variant records by comparing their REF and ALT alleles against a reference genome FASTA file, helping to identify discrepancies and annotate variants with their relation to the reference sequence.

## Usage

```bash
VCFX_ref_comparator --reference <reference.fasta> < input.vcf > annotated.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `-r`, `--reference` <FASTA> | Required. Path to reference genome in FASTA format |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_ref_comparator analyzes VCF variants by comparing them to a reference genome. The tool:

1. Loads a reference genome from a specified FASTA file
2. Processes each variant in the input VCF file
3. For each variant, compares the REF field with the corresponding sequence in the reference genome
4. Also determines if each ALT allele matches the reference sequence
5. Annotates each variant with a `REF_COMPARISON` tag in the INFO field, indicating the result of the comparison
6. Outputs an annotated VCF with all original fields preserved

This tool is particularly useful for:
- Validating the accuracy of variant calls
- Identifying potential errors in variant representation
- Distinguishing true variants from reference matching records
- Quality control of variant datasets

## Output Format

The output is a valid VCF file with the same format as the input, but with an additional `REF_COMPARISON` field added to the INFO column of each variant line. The output also includes a new header line defining the `REF_COMPARISON` INFO field.

The `REF_COMPARISON` field can have the following values:
- `REF_MISMATCH`: The REF allele does not match the reference genome
- `REF_MATCH`: The REF allele matches the reference genome
- `NOVEL`: The variant's ALT allele differs from the reference sequence
- `ALT_IS_REF`: The ALT allele matches the reference sequence (potential reference/alternate swap)
- `UNKNOWN_CHROM`: The chromosome is not found in the reference genome
- `INVALID_POS`: The position is out of bounds for the chromosome

## Examples

### Basic Usage

```bash
# Compare variants against a reference genome
VCFX_ref_comparator --reference genome.fa < input.vcf > validated.vcf
```

### Filtering for Reference Mismatches

```bash
# Find variants where the REF allele doesn't match the reference genome
VCFX_ref_comparator --reference genome.fa < input.vcf | grep "REF_MISMATCH" > mismatches.vcf
```

### Identifying ALT Alleles that Match Reference

```bash
# Find variants where the ALT allele actually matches the reference
VCFX_ref_comparator --reference genome.fa < input.vcf | grep "ALT_IS_REF" > potential_swaps.vcf
```

### Checking for Invalid Coordinates

```bash
# Identify variants with invalid chromosomes or positions
VCFX_ref_comparator --reference genome.fa < input.vcf | grep -E "UNKNOWN_CHROM|INVALID_POS" > invalid_coords.vcf
```

## Reference Comparison Process

The tool performs these steps for each variant:

1. Checks if the variant's chromosome exists in the reference genome
2. Verifies that the position is valid within the chromosome's sequence
3. Extracts the reference sequence at the specified position, matching the length of the REF allele
4. Compares the extracted sequence with the REF allele
5. For each ALT allele, determines if it matches the reference sequence

All comparisons are case-insensitive, and the reference genome is converted to uppercase during loading.

## Handling Special Cases

- **Chromosome not found**: If a chromosome in the VCF is not found in the reference genome, the variant is marked with `UNKNOWN_CHROM`
- **Position out of bounds**: If a position exceeds the length of the chromosome, the variant is marked with `INVALID_POS`
- **Multiple ALT alleles**: Each ALT allele is compared separately, and the result is included in the annotation
- **Symbolic alleles**: Not specially handled; will likely result in `REF_MATCH,NOVEL` annotations
- **Empty lines**: Preserved with a single newline
- **Header lines**: Preserved with a new INFO definition line added before the #CHROM header
- **Malformed VCF lines**: Lines with fewer than 8 columns are skipped with a warning
- **Data before header**: Skipped with a warning

## Performance

The tool is designed with the following considerations:

1. The entire reference genome is loaded into memory for fast random access
2. Chromosome names are converted to uppercase for case-insensitive matching
3. Whitespace is removed from FASTA sequences during loading
4. The VCF file is processed line by line, avoiding loading the entire file into memory
5. Only the required fields from each variant line are extracted and processed

For extremely large reference genomes, memory usage may be significant.

## Limitations

1. Requires loading the entire reference genome into memory
2. Limited to exact string comparison; no alignment is performed for complex variants
3. No special handling for symbolic alleles (like <DEL>, <INS>, etc.)
4. Does not normalize variants before comparison
5. Cannot handle reference genomes with duplicate chromosome names
6. No support for compressed reference files; FASTA must be uncompressed
7. No support for validating only a subset of variants or chromosomes 
