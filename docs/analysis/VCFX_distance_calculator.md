# VCFX_distance_calculator

## Overview

VCFX_distance_calculator analyzes a VCF file and calculates the distance (in base pairs) between consecutive variants along each chromosome, providing insights into variant density and spacing across the genome.

## Usage

```bash
# Standard mode (stdin)
VCFX_distance_calculator [OPTIONS] < input.vcf > variant_distances.tsv

# File mode (optimized for large files)
VCFX_distance_calculator -i input.vcf > variant_distances.tsv

# Positional argument mode
VCFX_distance_calculator input.vcf > variant_distances.tsv
```

## Options

| Option | Description |
|--------|-------------|
| `-i`, `--input FILE` | Input VCF file (uses memory-mapped I/O for best performance) |
| `-q`, `--quiet` | Suppress summary statistics to stderr |
| `-h`, `--help` | Display help message and exit |
| `-v`, `--version` | Show program version and exit |

## Description

VCFX_distance_calculator processes a VCF file to measure the distance between variants on the same chromosome. The tool:

1. Reads a VCF file line-by-line
2. Extracts chromosome (CHROM) and position (POS) information from each valid variant
3. For each chromosome, tracks the position of the previous variant
4. Calculates the distance from the previous variant to the current one
5. Outputs a tab-delimited file with the results
6. Provides summary statistics to stderr, including minimum, maximum, and average distances per chromosome

This tool is useful for:
- Analyzing variant density across the genome
- Identifying regions with unusually sparse or dense variant coverage
- Quality control to detect potential issues with variant calling
- Understanding the distribution of variants in targeted sequencing

## Output Format

The output is a tab-delimited text file with the following columns:

```
CHROM  POS  PREV_POS  DISTANCE
```

Where:
- CHROM: The chromosome name
- POS: The position of the current variant
- PREV_POS: The position of the previous variant on the same chromosome (or "NA" for the first variant)
- DISTANCE: The distance in base pairs between current and previous positions (or "NA" for the first variant)

## Examples

### Basic Usage

```bash
./VCFX_distance_calculator < input.vcf > variant_distances.tsv
```

### Analyzing Specific Chromosomes

```bash
# Extract only chromosome 1 data
grep -P "^chr1\t|^CHROM" variant_distances.tsv > chr1_distances.tsv
```

### Identifying Large Gaps

```bash
# Find regions with large gaps (>100,000 bp)
./VCFX_distance_calculator < input.vcf | awk -F'\t' '$4 > 100000 {print}' > large_gaps.tsv
```

### Visualizing Distance Distribution

```bash
# Process output for visualization (e.g., with R or Python)
./VCFX_distance_calculator < input.vcf | \
  grep -v "NA" | cut -f1,4 > distances_for_plotting.tsv
```

## Summary Statistics

In addition to the main output file, VCFX_distance_calculator prints summary statistics to stderr:

```
=== Summary Statistics ===
Chromosome: chr1
  Variants compared: 501
  Distances computed: 500
  Total distance: 10000000
  Min distance: 1
  Max distance: 150000
  Average distance: 20000
```

This provides a quick overview of variant distribution patterns for each chromosome.

## Handling Special Cases

- **First variant on a chromosome**: Marked with "NA" for PREV_POS and DISTANCE
- **Unsorted VCF files**: Processes variants in the order they appear, which may result in negative distances
- **Duplicate positions**: Correctly calculates a distance of 0 between variants at the same position
- **Malformed lines**: Warns about and skips lines that don't follow VCF format
- **Missing header**: Requires a proper VCF header (#CHROM line) before processing variant records
- **Invalid chromosome names**: Skips variants with obviously invalid chromosome names
- **Non-numeric positions**: Skips variants where the position cannot be parsed as an integer

## Performance

The tool offers two execution modes with different performance characteristics:

### File Mode (-i/--input)
When using the `-i` option, the tool uses memory-mapped I/O (mmap) for optimal performance on large files. This mode provides approximately **16x speedup** compared to stdin mode.

| File Size | Variants | Stdin Mode | File Mode | Speedup |
|-----------|----------|------------|-----------|---------|
| 4 GB      | 427K     | ~2m 19s    | ~8.7s     | ~16x    |

### Stdin Mode (default)
For smaller files or pipeline usage, the stdin mode processes input line by line with minimal memory requirements.

### Recommendations
- For files > 100MB, use `-i input.vcf` for best performance
- For pipelines or compressed input, use stdin mode with `zcat file.vcf.gz | VCFX_distance_calculator`

### Performance Characteristics
- Uses hash maps for O(1) lookups of previous positions
- Memory-mapped I/O for efficient large file processing
- 1MB output buffer for optimized write performance
- Zero-copy field extraction with memchr
- Memory usage scales with the number of distinct chromosomes, not with file size

## Limitations

- Does not account for chromosome lengths (cannot detect missing regions)
- Does not distinguish between different types of variants
- Assumes variants are properly formatted according to VCF specifications
- No built-in filtering for quality or other variant attributes
- Distances are calculated based on the reference genome coordinates, not actual sequence lengths
- Does not handle structural variants in any special way (uses only the position field) 
