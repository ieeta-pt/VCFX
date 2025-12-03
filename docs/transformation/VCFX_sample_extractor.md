# VCFX_sample_extractor

## Overview

VCFX_sample_extractor is a tool that extracts a subset of samples from a VCF file, allowing you to create a smaller, focused VCF containing only the samples of interest.

## Usage

```bash
VCFX_sample_extractor [OPTIONS] < input.vcf > subset.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `-s`, `--samples` LIST | Comma or space separated list of sample names to extract |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_sample_extractor reads a VCF file from standard input, identifies the samples specified in the command line, and produces a new VCF file containing only those samples. This is useful for:

- Reducing file size by extracting only relevant samples
- Creating sample-specific VCF files for specialized analyses
- Focusing on specific cohorts or subgroups
- Compliance with data sharing permissions that allow sharing only specific samples

The tool:
1. Reads the VCF header to identify sample columns
2. Maintains all meta-information and header lines
3. Extracts only the specified samples, preserving order and data integrity
4. Warns about any requested samples that aren't found in the input VCF

## Output Format

The output is a standard VCF file containing:
- All header lines from the input file
- A modified #CHROM header line that includes only the selected samples
- All variant lines from the input with only the selected sample columns

## Examples

### Extract a Single Sample

```bash
./VCFX_sample_extractor --samples "SAMPLE1" < input.vcf > single_sample.vcf
```

### Extract Multiple Samples with Comma Delimiter

```bash
./VCFX_sample_extractor --samples "SAMPLE1,SAMPLE2,SAMPLE3" < input.vcf > subset.vcf
```

### Extract Multiple Samples with Space Delimiter

```bash
./VCFX_sample_extractor --samples "SAMPLE1 SAMPLE2 SAMPLE3" < input.vcf > subset.vcf
```

### Process Large Files

```bash
# Extract a few samples from a large compressed VCF
zcat large_file.vcf.gz | ./VCFX_sample_extractor --samples "SAMPLE1,SAMPLE2" | gzip > subset.vcf.gz
```

## Handling Special Cases

- **Missing samples**: If a requested sample isn't found in the input VCF, a warning is issued but processing continues with the samples that were found
- **No samples found**: If none of the requested samples are found in the input VCF, the output will contain only the header and variant lines with no sample columns
- **Malformed VCF**: Lines with fewer than 8 columns are skipped with a warning
- **No sample columns**: Input variant lines without sample columns (fewer than 10 columns) are skipped
- **Empty sample names**: Empty sample names in the input list are ignored

## Performance

The tool processes VCF files line by line, with minimal memory requirements even for very large VCF files. Performance scales with:
- Number of samples in the input VCF (parsing time)
- Number of samples being extracted (output size)

## Limitations

- No wildcards or regular expressions for sample name matching
- Cannot extract samples based on properties or metadata
- Cannot reorder samples in the output file (order follows the original VCF)
- No option to rename samples in the output file 
