# VCFX_info_aggregator

## Overview

`VCFX_info_aggregator` is a tool that reads a VCF file from standard input, outputs it unmodified, and appends a summary section containing aggregated statistics (sum and average) for specified numeric INFO fields. The summary section is formatted as VCF header-like comments to maintain compatibility with VCF parsers.

## Usage

```
VCFX_info_aggregator [OPTIONS] < input.vcf > output.vcf
```

## Options

- `-a`, `--aggregate-info <fields>`: Comma-separated list of INFO fields to aggregate (required).
- `-h`, `--help`: Display help message and exit (handled by `vcfx::handle_common_flags`)
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

`VCFX_info_aggregator` processes a VCF file line by line, performing the following operations:

1. Parse command-line arguments to determine which INFO fields to aggregate.
2. Read each line from standard input.
3. Output each line unmodified to standard output.
4. For data lines (non-header lines):
   - Parse the INFO column (8th column).
   - Extract specified INFO fields if they exist and contain numeric values.
   - Accumulate these numeric values for later calculation.
5. After processing all lines, append a summary section that includes:
   - A line starting with `#AGGREGATION_SUMMARY`.
   - For each specified INFO field, a line reporting the sum and average of all valid numeric values found.

The tool is particularly useful for calculating summary statistics across an entire VCF file, such as average depth of coverage (DP) or average allele frequency (AF).

## Output Format

The output is identical to the input VCF file, with additional summary lines appended at the end:

```
... original VCF content ...
#AGGREGATION_SUMMARY
FIELD1: Sum=<total>, Average=<mean>
FIELD2: Sum=<total>, Average=<mean>
```

The summary lines start with `#` to ensure they are treated as comments by standard VCF parsers, maintaining compatibility.

## Examples

### Basic Usage

Calculate summary statistics for the depth (DP) field:

```bash
VCFX_info_aggregator --aggregate-info "DP" < input.vcf > output.vcf
```

### Multiple Fields

Calculate statistics for both depth (DP) and allele frequency (AF):

```bash
VCFX_info_aggregator --aggregate-info "DP,AF" < input.vcf > output.vcf
```

### Combining with Other Tools

Used in a pipeline to analyze filtered variants:

```bash
VCFX_record_filter --quality ">30" < input.vcf | VCFX_info_aggregator --aggregate-info "DP,AF,MQ" > filtered_with_stats.vcf
```

## Handling Special Cases

The tool implements several strategies for handling edge cases:

1. **Non-numeric values**: If a field cannot be parsed as a numeric value (e.g., "DP=abc"), it is skipped and not included in the aggregation.
2. **Missing fields**: If a specified INFO field is not present in a particular variant, it is simply skipped for that variant.
3. **Empty input**: The tool will process empty files correctly, reporting zeros for sums and averages.
4. **Malformed VCF**: If a data line is encountered before the `#CHROM` header, an error is reported.
5. **Line endings**: The tool correctly handles LF and CRLF line endings.
6. **Partial final line**: The tool properly processes files that do not end with a newline character.

## Performance

`VCFX_info_aggregator` is designed to be memory-efficient and performant:

1. It processes the VCF file in a single pass, with O(1) memory usage relative to the file size.
2. Only the specified INFO fields are parsed and accumulated, minimizing unnecessary processing.
3. The tool streams data directly from input to output, without storing the entire file in memory.

## Limitations

1. The tool only processes numeric values in INFO fields. String, flag, or complex values are not aggregated.
2. For multi-allelic variants, INFO fields like AF that may contain multiple values (comma-separated) are processed as a single value.
3. The tool does not modify or annotate the original variants; it only appends summary statistics at the end of the file.
4. The summary section, while formatted as VCF comments, may not be recognized by all downstream tools that expect a strictly-conforming VCF file.
5. There is no option to calculate additional statistics beyond sum and average (e.g., median, standard deviation). 
