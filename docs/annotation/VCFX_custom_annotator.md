# VCFX_custom_annotator

## Overview

VCFX_custom_annotator adds custom annotations to VCF files by matching variants from a user-provided annotation file and inserting them into the INFO field. This tool is particularly useful for incorporating external annotations, functional predictions, or custom labels into your VCF files.

## Usage

```bash
VCFX_custom_annotator --add-annotation <annotations.txt> [OPTIONS] < input.vcf > annotated.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `-a`, `--add-annotation <file>` | Required. Path to the annotation file containing the custom annotations |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_custom_annotator is designed to enhance VCF files with custom annotations from an external source. The tool:

1. Reads an annotation file where each line contains: `CHROM POS REF ALT annotation_value`
2. Creates a lookup map using the variant coordinates and alleles as keys
3. Processes the input VCF file line by line
4. For each variant, generates a key from its coordinates and alleles
5. Looks up any matching annotations
6. Adds the annotation to the INFO field as `CustomAnnotation=value`
7. For multi-allelic variants, handles each alternate allele separately
8. Outputs the annotated VCF to standard output

This streamlined process allows for efficient annotation of VCF files with minimal computational overhead.

## Annotation File Format

The annotation file should contain tab-separated fields with the following columns:

```
CHROM POS REF ALT annotation_value
```

Example:
```
1  100  A  G  HighImpact
1  200  T  C  ModerateImpact
2  300  G  A  LowImpact
```

- **CHROM**: Chromosome name (must match the VCF)
- **POS**: Position (1-based, must match the VCF)
- **REF**: Reference allele (must match the VCF)
- **ALT**: Alternate allele
- **annotation_value**: The annotation text to add (can include spaces after the first 4 fields)

## Output Format

The output is a valid VCF file that includes:

1. All original header lines
2. A new INFO field definition: `##INFO=<ID=CustomAnnotation,Number=.,Type=String,Description="Custom annotations added by VCFX_custom_annotator (multi-allelic)">`
3. All original variant lines with the added `CustomAnnotation=value` in the INFO column
4. For multi-allelic variants, comma-separated annotation values corresponding to each ALT allele

For variants without a matching annotation, the value "NA" is used.

## Examples

### Basic Usage

```bash
# Annotate a VCF file with functional impact predictions
./VCFX_custom_annotator --add-annotation impact_predictions.txt < input.vcf > annotated.vcf
```

### Viewing Annotated Results

```bash
# Annotate and view the first few variants
./VCFX_custom_annotator --add-annotation annotations.txt < input.vcf | head -n 20
```

### Filtering Based on Annotations

```bash
# Annotate variants and filter to keep only those with "HighImpact"
./VCFX_custom_annotator --add-annotation annotations.txt < input.vcf | grep "CustomAnnotation=HighImpact" > high_impact_variants.vcf
```

## Multi-allelic Variant Handling

VCFX_custom_annotator properly handles multi-allelic variants by:

1. Parsing and splitting the ALT field in the VCF on commas
2. Looking up annotations for each REF→ALT pair separately
3. Combining the annotations into a comma-separated list in the same order as the ALT alleles
4. Using "NA" as a placeholder when no annotation is found for a specific allele

For example, if a variant has `ALT=G,C,T` and annotations exist for the G and T alleles but not C, the result will be `CustomAnnotation=annotation_G,NA,annotation_T`.

## Handling Special Cases

The tool implements several strategies for handling edge cases:

1. **Missing annotations**: Uses "NA" when no annotation is found for a variant
2. **Empty annotation file**: Results in "NA" for all variants
3. **Malformed annotation lines**: Invalid lines in the annotation file are skipped with a warning
4. **Missing annotation file**: Reports an error if the annotation file cannot be opened
5. **Empty INFO field**: If the original INFO field is "." (missing), it's replaced with the new annotation
6. **Existing INFO content**: If the INFO field has existing content, annotations are appended
7. **Multi-allelic variants**: Each allele is handled separately with proper ordering

## Performance

VCFX_custom_annotator is designed for efficiency:

1. Annotation file is loaded into memory as a hash map for O(1) lookups
2. VCF file is processed in a streaming fashion, with minimal memory overhead
3. Line-by-line processing allows handling of arbitrarily large VCF files
4. String operations are optimized to minimize unnecessary copies
5. Capable of processing thousands of variants per second on typical hardware

## Limitations

1. Requires exact matching of CHROM, POS, REF, and ALT fields
2. Cannot perform fuzzy matching or coordinate-based lookups
3. All annotations are stored in memory, which may be a limitation for extremely large annotation files
4. Limited to a single annotation per variant (although annotations can contain multiple pieces of information)
5. No built-in option to modify the name of the added INFO field from the default "CustomAnnotation"
6. No support for annotating based on overlapping regions rather than exact positions
7. Annotations with internal commas may cause parsing issues in multi-allelic contexts 
