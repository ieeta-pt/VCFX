# VCFX_phred_filter

## Overview

VCFX_phred_filter filters a VCF file based on the PHRED quality scores (QUAL column), removing variants that fall below a specified quality threshold to focus analysis on higher confidence variant calls.

## Usage

```bash
VCFX_phred_filter [OPTIONS] < input.vcf > filtered.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `-p`, `--phred-filter` <THRESHOLD> | Set PHRED quality score threshold (default: 30.0) |
| `-k`, `--keep-missing-qual` | Keep variants with missing quality values (represented as ".") |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_phred_filter examines the QUAL column of each variant record in a VCF file and filters out variants with quality scores below a specified threshold. The tool:

1. Processes a VCF file line by line
2. Extracts the PHRED quality score from the QUAL column (position 6) of each variant
3. Compares the quality score against the specified threshold
4. Retains variants where the quality meets or exceeds the threshold
5. Passes through all header lines unchanged

This tool is particularly useful for:
- Removing low-confidence variant calls for downstream analysis
- Quality control of variant data
- Reducing false positives in variant datasets
- Applying consistent quality standards across multiple VCF files

## Output Format

The output is a standard VCF file with the same format as the input, but containing only variants with PHRED quality scores at or above the specified threshold. All header lines are preserved.

## Examples

### Basic Usage with Default Threshold

```bash
# Filter variants using the default quality threshold (30)
VCFX_phred_filter < input.vcf > filtered.vcf
```

### Setting a Custom Threshold

```bash
# Keep only variants with quality scores of 20 or higher
VCFX_phred_filter -p 20 < input.vcf > q20_filtered.vcf
```

### Keeping Variants with Missing Quality

```bash
# Keep variants with quality ≥ 30 or missing quality values
VCFX_phred_filter -p 30 -k < input.vcf > high_quality_with_missing.vcf
```

### In a Pipeline

```bash
# Combine with other filters in a pipeline
cat input.vcf | \
VCFX_phred_filter -p 40 | \
grep "PASS" > high_quality_pass.vcf
```

### Using Long Option Format

```bash
# Using long options for clarity
VCFX_phred_filter --phred-filter 25 --keep-missing-qual < input.vcf > filtered.vcf
```

## Quality Score Handling

The tool processes quality scores as follows:

1. **Standard numeric values**: Directly compared to the threshold
2. **Missing quality values** (represented as "."): 
   - By default, treated as 0.0 (filtered out)
   - With the `-k` option, treated as extremely high (effectively infinity) to ensure they pass
3. **Invalid quality values**: Logged with a warning and treated as 0.0 (filtered out)

## Handling Special Cases

- **Missing quality values**: Can be retained with the `-k` option
- **Invalid quality formats**: Treated as 0.0 with a warning
- **Empty lines**: Preserved with a single newline
- **Header lines**: Preserved unchanged
- **Malformed VCF lines**: Lines with fewer than 6 columns are skipped with a warning
- **Data before header**: Skipped with a warning

## Performance

The tool is designed for efficiency:

1. Processes VCF files line by line, with minimal memory requirements
2. Simple numeric comparison for fast filtering decisions
3. No requirement to load the entire file into memory
4. Fast string parsing for quality values

## Limitations

1. Focuses only on the QUAL column, not on per-sample or per-genotype quality metrics
2. No support for filtering based on other quality-related fields (e.g., GQ in the genotype fields)
3. Cannot apply different thresholds to different variant types
4. No option to exclude variants based on upper quality bounds
5. Doesn't provide statistics on filtered variants
6. Cannot filter based on multiple quality metrics simultaneously
7. Lacks region-specific filtering capabilities 
