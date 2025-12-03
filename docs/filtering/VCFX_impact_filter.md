# VCFX_impact_filter

## Overview

VCFX_impact_filter filters VCF variants based on their predicted functional impact level (HIGH, MODERATE, LOW, or MODIFIER) found in the INFO field of VCF records.

## Usage

```bash
VCFX_impact_filter --filter-impact <LEVEL> < input.vcf > filtered.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `-i`, `--filter-impact <LEVEL>` | Required. Impact level threshold. Must be one of: HIGH, MODERATE, LOW, MODIFIER |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_impact_filter analyzes variant annotations in a VCF file and filters them based on a specified impact level threshold. The tool:

1. Reads a VCF file line-by-line
2. For each variant, extracts the impact level from the INFO field (looks for `IMPACT=...`)
3. Classifies the impact into one of four levels: HIGH, MODERATE, LOW, or MODIFIER
4. Keeps only variants with impact level greater than or equal to the specified threshold
5. Adds an `EXTRACTED_IMPACT` field to the INFO column of retained variants
6. Outputs the filtered VCF with the same format as the input

The impact level hierarchy used for filtering is:
HIGH > MODERATE > LOW > MODIFIER > UNKNOWN

## Output Format

The output is a valid VCF file containing only variants that meet or exceed the specified impact threshold. Each retained variant will have an additional INFO field:

```
EXTRACTED_IMPACT=<value>
```

Where `<value>` is the original impact value extracted from the INFO field.

## Examples

### Filter for HIGH impact variants only

```bash
./VCFX_impact_filter --filter-impact HIGH < input.vcf > high_impact_variants.vcf
```

### Filter for MODERATE or higher impact variants

```bash
./VCFX_impact_filter --filter-impact MODERATE < input.vcf > functional_variants.vcf
```

### Combining with other tools

```bash
# Filter by impact then convert to another format
./VCFX_impact_filter --filter-impact HIGH < input.vcf | \
  ./VCFX_format_converter --format=bed > high_impact.bed
```

## Handling Special Cases

- **Case insensitivity**: Impact values are case insensitive (e.g., "high" and "HIGH" are treated the same)
- **Extended impact values**: Values like "HIGH_MISSENSE" are recognized by looking for the presence of standard impact keywords
- **Missing IMPACT field**: Variants without an IMPACT field in the INFO column are treated as "UNKNOWN" and filtered out by default
- **Empty INFO field**: Properly handled by adding the EXTRACTED_IMPACT field as the only INFO attribute
- **Multiple impact annotations**: If multiple IMPACT fields are present, only the first one is considered
- **Invalid impact values**: Any impact value not recognized as one of the four standard levels is classified as "UNKNOWN"

## Performance

The tool is optimized for efficiency:
- Processes VCF files line-by-line with minimal memory overhead
- Uses regular expressions for reliable pattern matching
- Processes very large VCF files with linear time complexity

## Limitations

- Only extracts and analyzes the first IMPACT field found in the INFO column
- Cannot differentiate between more detailed impact subclassifications (relies on basic HIGH/MODERATE/LOW/MODIFIER keywords)
- Assumes that functional impact annotations follow standard convention with one of the four recognized impact levels
- Does not account for the specific variant type (SNP, indel, etc.) when filtering
- No built-in options to combine impact filtering with other criteria (e.g., allele frequency) 
