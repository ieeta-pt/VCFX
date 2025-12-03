# VCFX_quality_adjuster

## Overview

`VCFX_quality_adjuster` transforms the QUAL field values in a VCF file by applying mathematical functions such as logarithm, square root, or square. This tool is useful for scaling quality scores to make them more interpretable or to prepare them for downstream analysis.

## Usage

```bash
VCFX_quality_adjuster [OPTIONS] < input.vcf > output.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `-a`, `--adjust-qual <FUNC>` | Required. The transformation function to apply. Must be one of: `log`, `sqrt`, `square`, or `identity`. |
| `-n`, `--no-clamp` | Do not clamp negative or extremely large values resulting from transformations. |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

`VCFX_quality_adjuster` processes a VCF file line by line, applying a specified mathematical transformation to the QUAL field (6th column) of each variant record. The tool:

1. Reads the VCF file from standard input.
2. Identifies header lines (beginning with #) and passes them through unchanged.
3. For data lines, extracts the QUAL value and applies the specified transformation.
4. By default, clamps negative values to 0 and extremely large values to 10^12.
5. Writes the modified VCF to standard output.

### Supported Transformations

- **log**: Applies natural logarithm (ln) to the quality score. Small constant (10^-10) is added to prevent log(0).
- **sqrt**: Applies square root to the quality score. Negative values are treated as 0.
- **square**: Multiplies the quality score by itself.
- **identity**: No transformation, passes the quality score unchanged.

## Output Format

The output is a VCF file with the same format as the input, but with transformed QUAL values. All other fields remain unchanged, maintaining full compatibility with standard VCF parsers and tools.

## Examples

### Logarithmic Transformation

```bash
# Transform quality scores using natural logarithm
VCFX_quality_adjuster --adjust-qual log < input.vcf > log_transformed.vcf
```

### Square Root Transformation

```bash
# Apply square root to quality scores
VCFX_quality_adjuster --adjust-qual sqrt < input.vcf > sqrt_transformed.vcf
```

### Square Transformation without Clamping

```bash
# Square quality scores without clamping large values
VCFX_quality_adjuster --adjust-qual square --no-clamp < input.vcf > squared_unclamped.vcf
```

### In a Pipeline

```bash
# Filter variants and then transform quality scores
VCFX_record_filter --quality ">20" < input.vcf | VCFX_quality_adjuster --adjust-qual log > filtered_log_transformed.vcf
```

## Handling Special Cases

The tool implements several strategies for handling edge cases:

1. **Missing QUAL values**: Fields marked with `.` or empty fields are treated as 0.
2. **Non-numeric QUAL values**: Lines with non-numeric QUAL values generate a warning and are skipped.
3. **Negative results**: By default, negative values resulting from transformations (e.g., log of a value < 1) are clamped to 0. This behavior can be disabled with `--no-clamp`.
4. **Very large values**: Values above 10^12 are clamped to prevent numerical issues. This behavior can be disabled with `--no-clamp`.
5. **Malformed lines**: Lines with fewer than 8 fields generate a warning and are skipped.
6. **Empty lines**: Empty lines are preserved in the output.

## Performance

`VCFX_quality_adjuster` is designed to be efficient:

1. It processes the VCF file in a single pass, requiring minimal memory footprint.
2. All transformations are simple mathematical functions with constant-time complexity.
3. The tool streams data directly from input to output without storing the entire file in memory.

## Limitations

1. Only the QUAL field is transformed; other numeric fields (like INFO fields) remain unchanged.
2. No facility to define custom transformation functions beyond the four built-in options.
3. Cannot apply different transformations to different variants or regions in a single run.
4. Lack of options for formatting the transformed values (e.g., number of decimal places).
5. No built-in option to back-transform values to their original scale. 
