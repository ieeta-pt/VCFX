# VCFX_phase_quality_filter

## Overview

`VCFX_phase_quality_filter` is a tool for filtering variants based on their phasing quality (PQ) scores in the INFO field. It allows users to specify custom conditions for filtering variants based on their PQ values, supporting various comparison operators.

## Usage

```bash
VCFX_phase_quality_filter --filter-pq "PQ<OP><THRESHOLD>" < input.vcf > output.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `-h, --help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |
| `-f, --filter-pq` | Condition like 'PQ>30', 'PQ>=20', 'PQ!=10', etc. |

## Description

`VCFX_phase_quality_filter` reads a VCF file and filters variants based on their phasing quality scores:

1. Extracts the PQ value from the INFO field of each variant
2. Applies the specified comparison condition to the PQ value
3. Keeps variants that satisfy the condition
4. Discards variants that don't meet the criteria

The tool handles missing or invalid PQ values by treating them as 0.0, ensuring robust processing of incomplete data.

## Input Requirements

- Input must be a valid VCF file
- File can be piped through stdin
- PQ values should be present in the INFO field (optional)
- Supports both VCFv4.0 and VCFv4.2 formats

## Output Format

The output is a VCF file containing only the variants that satisfy the specified PQ condition:
- Preserves all header lines
- Maintains original VCF format
- Includes only passing variants
- Preserves all fields and annotations

## Examples

### Basic Usage

Filter variants with PQ > 30:

```bash
VCFX_phase_quality_filter --filter-pq "PQ>30" < input.vcf > high_pq.vcf
```

### Different Operators

Various comparison operators are supported:

```bash
# Keep variants with PQ >= 20
VCFX_phase_quality_filter --filter-pq "PQ>=20" < input.vcf > pq_ge_20.vcf

# Keep variants with PQ <= 15
VCFX_phase_quality_filter --filter-pq "PQ<=15" < input.vcf > pq_le_15.vcf

# Keep variants with PQ != 10
VCFX_phase_quality_filter --filter-pq "PQ!=10" < input.vcf > pq_ne_10.vcf
```

### Integration with Other Tools

Combine with other VCFX tools:

```bash
cat input.vcf | \
  VCFX_validator | \
  VCFX_phase_quality_filter --filter-pq "PQ>30" | \
  VCFX_metadata_summarizer
```

## Supported Operators

The tool supports the following comparison operators:
- `>` - Greater than
- `>=` - Greater than or equal to
- `<` - Less than
- `<=` - Less than or equal to
- `==` - Equal to
- `!=` - Not equal to

## PQ Value Handling

The tool handles PQ values in the following ways:
- Extracts PQ values from INFO field (e.g., "PQ=30")
- Treats missing PQ values as 0.0
- Handles invalid PQ values gracefully
- Supports decimal values

## Error Handling

The tool handles various error conditions:
- Missing --filter-pq argument
- Invalid condition syntax
- Malformed VCF lines
- Missing or invalid PQ values
- Invalid comparison operators

## Performance Considerations

- Processes input streamingly
- Minimal memory usage
- Efficient for both small and large VCF files
- No need to load entire file into memory

## Limitations

- Only filters based on PQ values
- Does not modify PQ values
- Does not validate VCF format (use VCFX_validator for validation)
- Treats missing PQ values as 0.0

## Common Use Cases

1. Filtering low-quality phasing results
2. Quality control of phased variants
3. Selecting high-confidence phased variants
4. Removing poorly phased variants
5. Quality-based subsetting of phased data

## Best Practices

1. Validate input VCF before filtering
2. Use appropriate PQ thresholds based on your data
3. Consider missing PQ values in your analysis
4. Combine with other quality filters
5. Document your filtering criteria 
