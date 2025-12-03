# VCFX_metadata_summarizer

## Overview

`VCFX_metadata_summarizer` is a tool that analyzes and summarizes key metadata from a VCF file. It provides a comprehensive overview of the file's structure, including counts of contigs, INFO fields, FILTER fields, FORMAT fields, samples, and variants.

## Usage

```bash
VCFX_metadata_summarizer [options] < input.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `-h, --help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

`VCFX_metadata_summarizer` reads a VCF file and generates a summary of its key structural elements:

1. Counts unique contigs defined in the header
2. Counts unique INFO fields
3. Counts unique FILTER fields
4. Counts unique FORMAT fields
5. Counts the number of samples
6. Counts the total number of variants

The tool processes the VCF file line by line, parsing both header metadata and variant records to build a complete summary of the file's structure.

## Input Requirements

- Input must be a valid VCF file
- File can be piped through stdin
- Supports both VCFv4.0 and VCFv4.2 formats
- Handles both single-sample and multi-sample VCFs

## Output Format

The output is a formatted text summary with the following structure:

```
VCF Metadata Summary:
---------------------
Number of unique contigs: <count>
Number of unique INFO fields: <count>
Number of unique FILTER fields: <count>
Number of unique FORMAT fields: <count>
Number of samples: <count>
Number of variants: <count>
```

## Examples

### Basic Usage

Summarize metadata from a VCF file:

```bash
VCFX_metadata_summarizer < input.vcf
```

### Integration with Other Tools

Combine with other VCFX tools:

```bash
cat input.vcf | VCFX_validator | VCFX_metadata_summarizer
```

### Example Output

For a minimal VCF file with one contig, one INFO field, and a single variant:

```
VCF Metadata Summary:
---------------------
Number of unique contigs: 1
Number of unique INFO fields: 1
Number of unique FILTER fields: 0
Number of unique FORMAT fields: 0
Number of samples: 0
Number of variants: 1
```

## Header Parsing

The tool parses the following types of header lines:
- `##contig=<ID=...>` - Contig definitions
- `##INFO=<ID=...>` - INFO field definitions
- `##FILTER=<ID=...>` - FILTER field definitions
- `##FORMAT=<ID=...>` - FORMAT field definitions
- `#CHROM...` - Column header line (for sample counting)

## Error Handling

The tool handles various input scenarios:
- Empty files
- Files with no header
- Files with no variants
- Files with missing metadata fields
- Files with inconsistent header structures

## Performance Considerations

- Processes input streamingly
- Memory usage scales with the number of unique metadata fields
- Efficient for both small and large VCF files
- No need to load entire file into memory

## Limitations

- Only counts presence of fields, not their values
- Does not validate field definitions
- Does not check for field consistency across variants
- Does not analyze variant content beyond counting

## Common Use Cases

1. Quick assessment of VCF file structure
2. Quality control of VCF file completeness
3. Verification of expected metadata presence
4. Sample count verification
5. Variant count verification

## Best Practices

1. Run on VCF files before processing
2. Use in combination with VCFX_validator
3. Check for expected field counts
4. Verify sample counts match expectations
5. Use as part of quality control pipelines 
