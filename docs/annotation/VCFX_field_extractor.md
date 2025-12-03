# VCFX_field_extractor

## Overview
`VCFX_field_extractor` is a tool designed to extract and format specific fields from VCF (Variant Call Format) files. It allows users to select and output particular fields from VCF records, including standard fields, INFO subfields, and sample-specific genotype fields in a tabular format.

## Usage
```bash
VCFX_field_extractor --fields "FIELD1,FIELD2,..." [OPTIONS] < input.vcf > output.tsv
```

## Options
| Option | Description |
|--------|-------------|
| `-f`, `--fields` | Required. Comma-separated list of fields to extract (no spaces between fields) |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description
`VCFX_field_extractor` processes a VCF file and extracts only the specified fields for each variant. The tool:

1. Reads a VCF file from standard input
2. Identifies and parses the VCF header
3. For each variant line, extracts the requested fields
4. Outputs the extracted fields in a tab-separated format

The tool can extract three types of fields:
- **Standard VCF fields**: `CHROM`, `POS`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`, `INFO`
- **INFO subfields**: Any key that appears in the INFO column (e.g., `DP`, `AF`, `TYPE`)
- **Sample-specific fields**: Fields from the genotype columns, specified as either:
  - `SampleName:Subfield` (e.g., `SAMPLE1:GT` for the genotype of SAMPLE1)
  - `S<number>:Subfield` (e.g., `S1:DP` for the depth of the first sample)

## Field Types

### Standard VCF Fields
These are the eight fixed columns in the VCF format:
- `CHROM`: Chromosome
- `POS`: Position
- `ID`: Variant identifier
- `REF`: Reference allele
- `ALT`: Alternate allele(s)
- `QUAL`: Quality score
- `FILTER`: Filter status
- `INFO`: Additional information

### INFO Subfields
Any key found in the INFO column can be extracted directly by name. For example:
- `DP`: Read depth
- `AF`: Allele frequency
- `TYPE`: Variant type

### Sample Fields
Sample fields are specified using one of these formats:
- `SampleName:Subfield` where `SampleName` is the exact sample name from the VCF header
- `S<number>:Subfield` where `<number>` is the 1-based index of the sample column

Common sample subfields include:
- `GT`: Genotype
- `DP`: Read depth
- `GQ`: Genotype quality
- Other format fields defined in the VCF

## Output Format
The tool produces a tab-separated values (TSV) file with:
- A header row containing the requested field names
- One row per variant, with each requested field value
- Missing or invalid fields represented as `.`

## Examples

### Basic Standard Fields
Extract chromosome, position, ID, reference, and alternate alleles:
```bash
VCFX_field_extractor --fields "CHROM,POS,ID,REF,ALT" < input.vcf > basic_fields.tsv
```

### INFO Fields
Extract depth and allele frequency:
```bash
VCFX_field_extractor --fields "CHROM,POS,DP,AF" < input.vcf > info_fields.tsv
```

### Sample Genotype Fields
Extract genotypes for specific samples:
```bash
VCFX_field_extractor --fields "CHROM,POS,SAMPLE1:GT,SAMPLE2:GT" < input.vcf > genotypes.tsv
```

### Sample Fields by Index
Extract genotypes using sample indices:
```bash
VCFX_field_extractor --fields "CHROM,POS,S1:GT,S2:GT" < input.vcf > genotypes_by_index.tsv
```

### Mixed Field Types
Combine different field types:
```bash
VCFX_field_extractor --fields "CHROM,POS,DP,AF,SAMPLE1:GT,SAMPLE1:DP" < input.vcf > mixed_fields.tsv
```

## Handling Special Cases

### Missing Fields
- If a requested field is not found in the VCF record, a `.` is output
- This applies to missing INFO fields, invalid sample names, or non-existent format fields

### Malformed Records
- The tool attempts to handle malformed VCF lines gracefully
- For lines with too few columns, missing fields are filled with `.`
- Invalid data types in numeric fields are preserved as they appear in the input

### Header-only Files
- If a VCF file contains only headers and no variant records, the tool outputs just the header row

## Performance Considerations
- The tool processes VCF files line-by-line, with minimal memory overhead
- Extraction scales linearly with input size and number of requested fields
- For large VCF files, consider extracting only the necessary fields to improve performance

## Limitations
- Cannot filter records (only extracts fields from all records)
- Cannot perform operations or calculations on the extracted fields
- Does not support complex expressions or conditionals
- Limited to tab-separated output format
- Cannot output field descriptions or metadata from the VCF header
- No direct support for multi-allelic splitting or normalization 
