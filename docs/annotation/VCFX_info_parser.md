# VCFX_info_parser

## Overview

VCFX_info_parser extracts and formats specific INFO fields from a VCF file into a tabular format for easier analysis. The tool processes a VCF file line by line, parses the INFO column, and outputs only the requested fields in a clean TSV format.

## Usage

```bash
VCFX_info_parser --info "FIELD1,FIELD2,..." < input.vcf > extracted_info.tsv
```

## Options

| Option | Description |
|--------|-------------|
| `-i`, `--info <FIELDS>` | Required. Comma-separated list of INFO fields to extract (e.g., "DP,AF,SOMATIC") |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_info_parser simplifies the process of extracting specific information from VCF files by:

1. Reading VCF data from standard input
2. Parsing the INFO column to extract user-specified fields
3. Producing a clean, tabular output with standardized headers
4. Properly handling flags, missing values, and malformed entries

This tool is particularly useful for:
- Extracting numeric values like depth (DP) or allele frequency (AF) for statistical analysis
- Converting complex VCF INFO fields into a format suitable for spreadsheet applications
- Creating simplified datasets focused on specific annotations
- Preparing data for visualization or report generation

## Output Format

The output is a tab-separated file with the following columns:

```
CHROM  POS  ID  REF  ALT  FIELD1  FIELD2  ...
```

Where:
- The first five columns are standard VCF fields (chromosome, position, ID, reference allele, alternate allele)
- Each subsequent column contains the value of a requested INFO field
- Missing values are represented by a dot (.)
- Flag fields (INFO fields with no value) are also represented by a dot (.)

## Examples

### Basic Usage - Extract Depth Information

```bash
./VCFX_info_parser --info "DP" < input.vcf > depth_data.tsv
```

### Extract Multiple Fields

```bash
./VCFX_info_parser --info "DP,AF,MQ" < input.vcf > key_metrics.tsv
```

### Working with Annotation Data

```bash
./VCFX_info_parser --info "Gene,IMPACT,Consequence" < annotated.vcf > gene_impacts.tsv
```

### Pipeline Example

```bash
# Filter a VCF file and extract specific INFO fields
cat input.vcf | grep "PASS" | ./VCFX_info_parser --info "DP,AF,SOMATIC" > filtered_annotations.tsv
```

## Handling Special Cases

The tool implements several strategies for handling edge cases:

1. **Flag fields**: INFO fields without values (flags like 'SOMATIC') are represented by a dot in the output
2. **Missing fields**: If a requested INFO field is not present in a specific variant, a dot is printed
3. **Malformed lines**: Lines that don't conform to VCF format are skipped with a warning message
4. **Empty input**: The tool correctly handles empty input files
5. **Header lines**: VCF header lines (starting with #) are skipped
6. **Line endings**: LF and CRLF line endings are supported
7. **Partial final line**: Files without a final newline character are processed correctly

## Performance

VCFX_info_parser is designed for efficiency:

1. Single-pass processing with line-by-line reading, allowing for streaming of very large files
2. Minimal memory footprint regardless of input file size
3. Efficient string parsing with no complex regular expressions
4. Fast lookup of INFO fields using hash maps

## Limitations

1. Cannot handle multi-allelic variants specially (each row is processed independently)
2. No built-in filtering capabilities (use in conjunction with other filtering tools)
3. Cannot split INFO fields with multiple values (e.g., CSQ fields from VEP)
4. Doesn't preserve VCF headers in the output
5. No option to include additional VCF columns (QUAL, FILTER) in the output
6. Cannot extract FORMAT fields or sample-specific information 
