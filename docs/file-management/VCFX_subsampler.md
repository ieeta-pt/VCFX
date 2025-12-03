# VCFX_subsampler

## Overview

`VCFX_subsampler` is a tool for randomly selecting a specified number of variants from a VCF file. It uses reservoir sampling to efficiently select a random subset of variants while preserving the VCF header information. The tool is particularly useful for creating smaller test datasets or reducing the size of large VCF files for preliminary analysis.

## Usage

```bash
VCFX_subsampler [options] < input.vcf > output.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `-s, --subsample <N>` | Required: Number of variants to keep in the output |
| `--seed <INT>` | Optional: Use a specific random seed for reproducible results |
| `-h, --help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

`VCFX_subsampler` processes a VCF file to:

1. Preserve all header lines (starting with #)
2. Randomly select N variants from the data section
3. Skip invalid lines (those with fewer than 8 columns)
4. Output the selected variants while maintaining VCF format

The tool uses reservoir sampling to ensure unbiased random selection, even when the total number of variants is unknown in advance.

## Input Requirements

- Must be a valid VCF file
- Can be piped through stdin
- Must have at least 8 columns (CHROM through INFO)
- Header lines must start with #

## Output Format

The output is a VCF file containing:
- All original VCF header lines
- N randomly selected variants (or all variants if input has fewer than N)
- Same format as input VCF
- Invalid lines (with <8 columns) are skipped with warnings

## Examples

### Basic Usage

Select 1000 random variants:

```bash
VCFX_subsampler --subsample 1000 < input.vcf > subset.vcf
```

### Reproducible Sampling

Use a fixed seed for reproducible results:

```bash
VCFX_subsampler --subsample 1000 --seed 1234 < input.vcf > subset.vcf
```

### Integration with Other Tools

Combine with other VCFX tools:

```bash
cat input.vcf | \
  VCFX_validator | \
  VCFX_subsampler --subsample 1000 | \
  VCFX_metadata_summarizer
```

## Sampling Algorithm

The tool uses reservoir sampling to:
- Process the input streamingly
- Maintain a reservoir of N variants
- Replace variants in the reservoir with probability N/count
- Ensure unbiased random selection

## Error Handling

The tool handles various error conditions:
- Missing --subsample argument
- Invalid subsample size (must be >0)
- Invalid seed value
- Invalid VCF lines (with <8 columns)
- Data lines before header

## Performance Considerations

- Memory efficient: only stores N variants in memory
- Processes input streamingly
- Preserves header information
- Skips invalid lines efficiently

## Limitations

- Only samples variants (data lines)
- Skips lines with <8 columns
- Requires at least 8 columns in VCF
- Skips data lines before #CHROM header
- No support for weighted sampling

## Common Use Cases

1. Creating test datasets
2. Reducing large VCF files for quick analysis
3. Generating random subsets for validation
4. Preparing data for development and testing
5. Creating smaller datasets for preliminary analysis

## Best Practices

1. Validate input VCF before sampling
2. Use --seed for reproducible results
3. Monitor warning messages for skipped lines
4. Consider using VCFX_validator before sampling
5. Document sampling parameters for reproducibility 
