# VCFX_region_subsampler

## Overview

`VCFX_region_subsampler` is a tool for filtering VCF variants based on genomic regions specified in a BED file. It keeps only variants whose positions fall within the specified regions, efficiently handling multiple regions and overlapping intervals.

## Usage

```bash
VCFX_region_subsampler --region-bed FILE < input.vcf > output.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `-h, --help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |
| `-b, --region-bed FILE` | BED file listing regions to keep |

## Description

`VCFX_region_subsampler` processes a VCF file and a BED file to:

1. Read and parse the BED file containing genomic regions
2. Convert 0-based BED coordinates to 1-based VCF coordinates
3. Merge overlapping or contiguous intervals for efficiency
4. Filter VCF variants to keep only those falling within specified regions
5. Preserve all VCF header information and variant details

The tool uses binary search for efficient region lookup and handles multiple regions per chromosome.

## Input Requirements

### VCF Input
- Must be a valid VCF file
- Can be piped through stdin
- Supports both VCFv4.0 and VCFv4.2 formats
- Must have at least 8 columns (CHROM through INFO)

### BED Input
- Standard BED format (chromosome, start, end)
- 0-based coordinates (automatically converted to 1-based)
- One region per line
- Supports multiple regions per chromosome
- Invalid lines are skipped with warnings

## Output Format

The output is a VCF file containing:
- All original VCF header lines
- Only variants falling within specified regions
- Original variant information preserved
- Same format as input VCF

## Examples

### Basic Usage

Filter variants using a single region:

```bash
VCFX_region_subsampler --region-bed regions.bed < input.vcf > filtered.vcf
```

### Multiple Regions

Filter using multiple regions across chromosomes:

```bash
# regions.bed:
chr1    0    100
chr2    100  200
VCFX_region_subsampler --region-bed regions.bed < input.vcf > filtered.vcf
```

### Integration with Other Tools

Combine with other VCFX tools:

```bash
cat input.vcf | \
  VCFX_validator | \
  VCFX_region_subsampler --region-bed regions.bed | \
  VCFX_metadata_summarizer
```

## Region Handling

### Coordinate System
- Input BED: 0-based coordinates
- Internal processing: 1-based coordinates
- Automatic conversion between systems

### Interval Merging
- Overlapping intervals are merged
- Contiguous intervals are combined
- Maintains efficiency for large region sets

### Region Validation
- Skips invalid BED lines
- Handles negative intervals
- Ignores zero-length intervals
- Reports warnings for invalid entries

## Error Handling

The tool handles various error conditions:
- Missing --region-bed argument
- Invalid BED file format
- Invalid VCF lines
- Missing or malformed coordinates
- Data lines before header

## Performance Considerations

- Uses binary search for region lookup
- Merges overlapping intervals for efficiency
- Processes input streamingly
- Memory efficient for large region sets
- Handles large VCF files

## Limitations

- Only filters by position (CHROM, POS)
- Does not validate VCF format (use VCFX_validator for validation)
- Requires at least 8 columns in VCF
- Skips data lines before #CHROM header
- Treats invalid BED lines as warnings, not errors

## Common Use Cases

1. Extracting variants from specific genomic regions
2. Focusing analysis on particular chromosomal segments
3. Creating region-specific VCF subsets
4. Preparing data for region-based analysis
5. Filtering variants for specific genomic features

## Best Practices

1. Validate input VCF before filtering
2. Verify BED file format and coordinates
3. Check region coverage before processing
4. Monitor warning messages for invalid regions
5. Document region selection criteria 
