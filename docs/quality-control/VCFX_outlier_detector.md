# VCFX_outlier_detector

## Overview

VCFX_outlier_detector identifies outliers in VCF data based on numeric metrics, operating in two modes: variant mode to detect outlier variants exceeding a threshold for a specified INFO field, and sample mode to identify samples with average metrics above a threshold.

## Usage

```bash
VCFX_outlier_detector --metric <KEY> --threshold <VAL> [--variant|--sample] < input.vcf > outliers.txt
```

## Options

| Option | Description |
|--------|-------------|
| `--metric`, `-m` <KEY> | Name of the metric to use (e.g., AF, DP, GQ) |
| `--threshold`, `-t` <VAL> | Numeric threshold value for outlier detection |
| `--variant`, `-V` | Variant mode: identify variants with INFO field metrics above threshold |
| `--sample`, `-s` | Sample mode: identify samples with average genotype metrics above threshold |
| `--help`, `-h` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

**Note:** `-v` shows the version information. Use `--variant` or the short option `-V` to run in variant mode.

## Description

VCFX_outlier_detector analyzes VCF files to identify outliers based on numeric metrics. The tool operates in two distinct modes:

1. **Variant Mode (default)**: 
   - Examines each variant's specified INFO field metric
   - Reports variants where the metric exceeds the specified threshold
   - Useful for finding variants with unusual characteristics (e.g., high allele frequency, depth)

2. **Sample Mode**:
   - Calculates the average value of a specified genotype metric for each sample
   - Reports samples where the average metric exceeds the specified threshold
   - Useful for identifying samples with unusual quality characteristics

The tool processes VCF files line by line, extracting the relevant metric from either the INFO field (variant mode) or the FORMAT/genotype fields (sample mode). For sample mode, it accumulates values across all variants to calculate the per-sample averages.

## Output Format

### Variant Mode

The output is a tab-delimited file with the following columns:
```
#CHROM  POS  ID  METRIC
```
Where:
- CHROM, POS, ID: Standard VCF fields for variant identification
- METRIC: The value of the specified metric for each outlier variant

### Sample Mode

The output is a tab-delimited file with the following columns:
```
#Sample  Average_METRIC
```
Where:
- Sample: The sample name
- Average_METRIC: The average value of the metric across all variants (or "NA" if below threshold)

## Examples

### Identify Variants with High Allele Frequency

```bash
# Find variants with allele frequency > 0.05
VCFX_outlier_detector --metric AF --threshold 0.05 --variant < input.vcf > high_af_variants.txt
```

### Identify Low-Quality Samples

```bash
# Find samples with average genotype quality > 30
VCFX_outlier_detector --metric GQ --threshold 30 --sample < input.vcf > high_quality_samples.txt
```

### Detect Unusual Depth Variants

```bash
# Find variants with unusually high depth
VCFX_outlier_detector --metric DP --threshold 100 --variant < input.vcf > high_depth_variants.txt
```

### Identify Samples with High Missing Rate

```bash
# Find samples with high average missing rate
VCFX_outlier_detector --metric MISSING --threshold 0.2 --sample < input.vcf > high_missing_samples.txt
```

## Metric Extraction

The tool implements two strategies for extracting metrics:

1. **INFO field parsing (variant mode)**:
   - Extracts fields with format `KEY=VALUE` from the INFO column
   - Converts the value to a numeric type for comparison

2. **Genotype field parsing (sample mode)**:
   - First checks if the metric is directly specified with `KEY=VALUE` in the genotype field
   - Otherwise, locates the metric position in the FORMAT field and extracts the corresponding value for each sample

## Handling Special Cases

- **Missing values**: In sample mode, metrics that can't be parsed or are missing are skipped in the average calculation
- **Invalid numeric values**: Non-numeric values are ignored with appropriate warnings
- **Empty files**: Properly handled, producing an appropriate output header
- **Malformed VCF lines**: Lines with too few columns are skipped
- **Non-standard FORMAT fields**: Both standard colon-delimited formats and custom KEY=VALUE formats are supported
- **No matching metrics**: If no instances of the metric are found, a warning is issued

## Performance

VCFX_outlier_detector is designed for efficiency:

1. **Variant Mode**: 
   - Single pass through the file with O(n) time complexity where n is the number of variants
   - Minimal memory usage, regardless of file size

2. **Sample Mode**:
   - Requires a single pass but tracks running sums and counts for each sample
   - Memory usage scales with the number of samples, not the number of variants

## Limitations

1. Only supports thresholding in one direction (greater than threshold)
2. No support for statistical outlier detection (e.g., z-scores or percentile-based methods)
3. Cannot filter based on multiple metrics in a single run
4. Sample mode requires the entire file to be processed before producing output
5. No built-in options for handling multi-allelic sites differently
6. Cannot detect outliers based on metadata not present in the VCF file
7. Only numeric metrics are supported; cannot detect categorical outliers 
