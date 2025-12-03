# VCFX_annotation_extractor

## Overview

VCFX_annotation_extractor extracts annotation fields from a VCF file's INFO column and converts them into a tabular format. The tool is particularly useful for extracting specific annotations (such as functional impact, gene name, or any custom annotation) from VCF files into a more analysis-friendly TSV format.

## Usage

```bash
VCFX_annotation_extractor --annotation-extract "FIELD1,FIELD2,..." < input.vcf > extracted.tsv
```

## Options

| Option | Description |
|--------|-------------|
| `-a`, `--annotation-extract <FIELDS>` | Required. Comma-separated list of INFO field annotations to extract (e.g., "ANN,Gene,Impact") |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_annotation_extractor simplifies the extraction and analysis of variant annotations by:

1. Reading a VCF file from standard input
2. Parsing the INFO column to extract user-specified annotation fields
3. Handling multi-allelic variants by creating separate rows for each ALT allele
4. Aligning per-allele annotations (such as ANN) with the corresponding ALT allele
5. Producing a clean tab-delimited output with standardized columns

This tool is particularly useful for:
- Converting complex VCF annotations into a format suitable for spreadsheet applications
- Extracting specific annotation fields for focused analysis
- Preparing variant annotation data for visualization or reporting
- Working with multi-allelic variants where annotations correspond to specific alleles

## Output Format

The output is a tab-separated (TSV) file with the following columns:

```
CHROM  POS  ID  REF  ALT  <ANNOTATION1>  <ANNOTATION2>  ...
```

Where:
- The first five columns are standard VCF fields (chromosome, position, ID, reference allele, alternate allele)
- Each subsequent column contains the value of a requested annotation field
- Missing values are represented by "NA"
- Multi-allelic variants are split into multiple rows, one for each ALT allele
- Per-allele annotations (like ANN) are properly aligned with their corresponding ALT allele

## Examples

### Basic Usage - Extract Gene Annotations

```bash
./VCFX_annotation_extractor --annotation-extract "Gene" < input.vcf > genes.tsv
```

### Extract Multiple Annotation Fields

```bash
./VCFX_annotation_extractor --annotation-extract "ANN,Gene,Impact,DP" < input.vcf > annotations.tsv
```

### Process and Filter in a Pipeline

```bash
# Extract annotations from only PASS variants
grep -e "^#" -e "PASS" input.vcf | ./VCFX_annotation_extractor --annotation-extract "ANN,Gene,Impact" > pass_annotations.tsv
```

### Analyze Impact Distribution

```bash
# Extract impact annotations and count occurrences
./VCFX_annotation_extractor --annotation-extract "Impact" < input.vcf | tail -n +2 | cut -f6 | sort | uniq -c
```

## Multi-allelic Variant Handling

The tool handles multi-allelic variants specially:

1. Each ALT allele in a multi-allelic variant gets its own row in the output
2. For Number=A annotations (like ANN) that have multiple comma-separated values, each value is aligned with the corresponding ALT allele
3. For single-value annotations (like Gene, Impact), the same value is used for all ALT alleles of a variant
4. If there are more ALT alleles than annotation values, "NA" is used for the excess ALT alleles

### Example
For a variant line with `ALT=T,G,C` and `ANN=missense,stop_gained,intergenic`:
- Three rows will be generated in the output (one for each ALT)
- The annotations will be properly aligned: T→missense, G→stop_gained, C→intergenic

## Handling Special Cases

The tool implements several strategies for handling edge cases:

1. **Missing annotations**: If a requested annotation is not found, "NA" is output
2. **Malformed VCF lines**: Lines with fewer than 8 columns are skipped with a warning
3. **Empty annotations**: Empty annotation values are preserved and not replaced with "NA"
4. **Multi-value annotations**: Currently, only ANN field is treated as multi-value and split by commas
5. **Header parsing**: The tool checks for proper VCF headers before processing data
6. **Empty input**: The tool correctly handles empty input files, producing only the header line
7. **Invalid characters**: The tool preserves all characters in annotation values, including special characters

## Performance

VCFX_annotation_extractor is designed for efficiency:

1. Single-pass processing reads the VCF file line-by-line without loading the entire file into memory
2. Efficient string parsing with optimized splitting functions
3. Uses hash maps for quick annotation lookups
4. Memory usage scales with the size of individual variant lines rather than the whole file
5. Output is streamed directly without intermediate storage

## Limitations

1. Currently, only the ANN field is recognized as a per-allele (Number=A) field that needs to be split; other Number=A fields are not automatically detected
2. No VCF header parsing to automatically determine which fields are Number=A vs. Number=1
3. Cannot extract FORMAT fields or sample-specific information
4. The output does not include QUAL or FILTER columns from the input VCF
5. No wildcard or regex support for selecting annotation fields
6. Annotation fields with embedded tab or newline characters may cause issues in the output format
7. Limited error recovery for malformed INFO fields 
