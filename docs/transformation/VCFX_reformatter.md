# VCFX_reformatter

## Overview

`VCFX_reformatter` is a tool for reformatting INFO and FORMAT fields in VCF files. It provides functionality to compress (remove) specific fields and reorder fields in both INFO and FORMAT columns, making VCF files more organized and efficient.

## Usage

```bash
VCFX_reformatter [options] < input.vcf > output.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `-h, --help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |
| `-c, --compress-info <keys>` | Remove specified INFO keys (comma-separated) |
| `-f, --compress-format <keys>` | Remove specified FORMAT keys (comma-separated) |
| `-i, --reorder-info <keys>` | Reorder INFO keys (comma-separated) |
| `-o, --reorder-format <keys>` | Reorder FORMAT keys (comma-separated) |

## Description

`VCFX_reformatter` modifies VCF files in several ways:

1. **INFO Field Compression**:
   - Removes specified keys from the semicolon-separated INFO field
   - Preserves remaining fields in their original order
   - Handles both key-value pairs and flag fields

2. **FORMAT Field Compression**:
   - Removes specified keys from the colon-separated FORMAT field
   - Updates all sample columns to match the new FORMAT structure
   - Maintains data consistency across all samples

3. **INFO Field Reordering**:
   - Places specified keys at the beginning of the INFO field
   - Appends remaining keys in their original order
   - Preserves all key-value pairs and flags

4. **FORMAT Field Reordering**:
   - Reorders the FORMAT column keys
   - Updates all sample columns to match the new order
   - Maintains data alignment across all samples

## Input Requirements

- Input must be a valid VCF file
- File can be piped through stdin
- Supports both VCFv4.0 and VCFv4.2 formats
- Handles both single-sample and multi-sample VCFs

## Output Format

The output is a VCF file with:
- All header lines preserved
- Modified INFO and FORMAT fields according to specifications
- Updated sample columns to match new FORMAT structure
- Original VCF format maintained

## Examples

### Basic Usage

Remove specific INFO fields and reorder others:

```bash
VCFX_reformatter --compress-info AF,DP --reorder-info AF,DP < input.vcf > output.vcf
```

### Format Field Manipulation

Remove and reorder FORMAT fields:

```bash
VCFX_reformatter --compress-format PL,AD --reorder-format GT,DP < input.vcf > output.vcf
```

### Combined Operations

Perform multiple operations in one command:

```bash
VCFX_reformatter \
  --compress-info AF,DP \
  --compress-format PL,AD \
  --reorder-info AF,DP \
  --reorder-format GT,DP \
  < input.vcf > output.vcf
```

### Integration with Other Tools

Combine with other VCFX tools:

```bash
cat input.vcf | \
  VCFX_validator | \
  VCFX_reformatter --compress-info AF,DP | \
  VCFX_metadata_summarizer
```

## Field Handling

### INFO Field Processing
- Handles key-value pairs (e.g., "DP=10")
- Handles flag fields (e.g., "PASS")
- Preserves field separators
- Maintains field order when specified

### FORMAT Field Processing
- Updates FORMAT column structure
- Modifies all sample columns accordingly
- Preserves data alignment
- Handles missing values (".")

## Error Handling

The tool handles various error conditions:
- Malformed VCF lines
- Missing fields
- Invalid field formats
- Inconsistent sample data
- Lines with fewer than 8 columns

## Performance Considerations

- Processes input streamingly
- Efficient memory usage
- Handles large files
- Preserves original data integrity

## Limitations

- Only modifies INFO and FORMAT fields
- Does not validate VCF format (use VCFX_validator for validation)
- Does not modify other VCF columns
- Requires at least 8 columns in data lines

## Common Use Cases

1. Removing unnecessary fields to reduce file size
2. Reordering fields for better readability
3. Standardizing VCF format across different sources
4. Preparing VCF files for specific analysis tools
5. Cleaning up VCF files before processing

## Best Practices

1. Validate input VCF before reformatting
2. Back up original files before modification
3. Verify output format meets requirements
4. Use appropriate field combinations
5. Document field modifications 
