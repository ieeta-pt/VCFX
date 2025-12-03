# VCFX_info_summarizer

## Overview

VCFX_info_summarizer analyzes numeric fields in the INFO column of a VCF file and calculates summary statistics (mean, median, and mode) for each specified field. This tool enables researchers to quickly understand the distribution and central tendencies of key metrics across variants.

## Usage

```bash
VCFX_info_summarizer --info "FIELD1,FIELD2,..." < input.vcf > summary_stats.tsv
```

## Options

| Option | Description |
|--------|-------------|
| `-i`, `--info <FIELDS>` | Required. Comma-separated list of INFO fields to analyze (e.g., "DP,AF,MQ") |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_info_summarizer processes a VCF file to generate statistical summaries of specified numeric INFO fields. The tool:

1. Reads a VCF file from standard input
2. Parses the INFO column from each variant record
3. Extracts numeric values for the specified INFO fields
4. Calculates three key statistics for each field:
   - Mean (average value)
   - Median (middle value)
   - Mode (most frequently occurring value)
5. Outputs the results in a clean, tabular format

This tool is valuable for:
- Quality control assessment of sequencing data
- Understanding the distribution of metrics like depth, allele frequency, or mapping quality
- Identifying potential biases or anomalies in variant calling
- Summarizing large VCF files for reports or visualizations

## Output Format

The output is a tab-separated file with the following columns:

```
INFO_Field  Mean  Median  Mode
```

Where:
- INFO_Field: The name of the INFO field being summarized
- Mean: The arithmetic mean of all numeric values for that field
- Median: The middle value when all values are sorted
- Mode: The most frequently occurring value
- "NA" is displayed if no valid numeric values were found for a field

All numeric values are formatted with four decimal places of precision.

## Examples

### Basic Usage - Analyze Depth Statistics

```bash
./VCFX_info_summarizer --info "DP" < input.vcf > depth_stats.tsv
```

### Analyze Multiple Fields

```bash
./VCFX_info_summarizer --info "DP,AF,MQ" < input.vcf > variant_stats.tsv
```

### Analyze Complex Input with Filtering

```bash
# Get summary stats for only PASS variants
grep -e "^#" -e "PASS" input.vcf | ./VCFX_info_summarizer --info "DP,QD,FS" > pass_variant_stats.tsv
```

### Filter and Compare Multiple VCF Files

```bash
# Create a summary comparison script for multiple files
for vcf in sample1.vcf sample2.vcf sample3.vcf; do
  echo "=== $vcf ===" >> summary.txt
  ./VCFX_info_summarizer --info "DP,AF" < $vcf >> summary.txt
  echo "" >> summary.txt
done
```

## Handling Special Cases

The tool implements several strategies for handling edge cases:

1. **Non-numeric values**: Skipped with a warning to stderr, without affecting the calculations for other values
2. **Missing fields**: If a specified field is not present in a variant's INFO column, it's simply skipped
3. **Multi-value fields**: Each comma-separated value is processed individually (e.g., AF=0.1,0.2,0.3)
4. **Empty input**: Outputs "NA" for all statistics if no valid values are found
5. **Malformed VCF**: Lines that don't conform to VCF format are skipped with a warning
6. **Header validation**: Checks for the presence of a proper #CHROM header line before processing records
7. **Flag fields**: INFO flags without values are treated as having a value of "1" for statistical calculations

## Performance

VCFX_info_summarizer is designed for efficiency with large VCF files:

1. Single-pass processing with O(n) time complexity where n is the number of variants
2. O(m) memory usage where m is the number of numeric values for the specified fields
3. Efficient string parsing using streams
4. Fast statistical calculations with minimal sorting operations

## Limitations

1. Cannot process non-numeric INFO fields (strings, flags, etc.) except for converting flags to "1"
2. No ability to filter variants based on their values (must be combined with other tools for filtering)
3. Limited to basic statistics (mean, median, mode); no advanced statistics like standard deviation, quartiles, etc.
4. Does not support weighted statistics for multi-allelic variants
5. Cannot process FORMAT fields or perform sample-specific statistical summaries
6. No support for histograms or graphical representations of the data distribution 
