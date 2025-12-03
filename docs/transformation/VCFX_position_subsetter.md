# VCFX_position_subsetter

## Overview
`VCFX_position_subsetter` extracts variants from a VCF file that fall within a specified genomic region, allowing for targeted analysis of specific chromosomal segments.

## Usage
```bash
VCFX_position_subsetter --region "CHR:START-END" < input.vcf > filtered.vcf
```

## Options
| Option | Description |
|--------|-------------|
| `-r`, `--region <CHR:START-END>` | Required. Genomic region to extract in the format "chromosome:start-end" |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description
`VCFX_position_subsetter` reads a VCF file from standard input and outputs only those variants that fall within the specified genomic region. The tool:

1. Processes a VCF file line by line
2. Extracts the chromosome and position from each variant
3. Compares them to the user-specified region
4. Retains variants where:
   - The chromosome matches exactly with the specified chromosome
   - The position falls within the specified start and end coordinates (inclusive)
5. Outputs all header lines unchanged
6. Outputs only variant lines that meet the criteria

This tool is particularly useful for:
- Focusing analysis on specific genomic regions of interest
- Extracting variants in a gene or regulatory region
- Reducing file size by narrowing down to relevant regions
- Preparing data for local analysis or visualization

## Output Format
The output is a standard VCF file containing:
- All original header lines from the input VCF
- Only the variant records that fall within the specified region
- The same format and structure as the original VCF

## Examples

### Basic Usage
Extract variants on chromosome 1 between positions 1,000,000 and 2,000,000:
```bash
VCFX_position_subsetter --region "chr1:1000000-2000000" < input.vcf > chr1_region.vcf
```

### Small Region
Extract variants in a small region, such as a specific exon:
```bash
VCFX_position_subsetter --region "chr17:41245000-41245500" < input.vcf > brca1_exon.vcf
```

### Different Chromosome Format
Extract variants using a different chromosome naming format (numeric):
```bash
VCFX_position_subsetter --region "2:150000-250000" < input.vcf > chr2_region.vcf
```

### Pipeline Integration
Use as part of a longer analysis pipeline:
```bash
cat input.vcf | VCFX_position_subsetter --region "chrX:5000000-6000000" | another_tool > final_output.vcf
```

## Region Parsing

### Format Requirements
The region must be specified in the format "CHR:START-END" where:
- CHR is the chromosome name, exactly as it appears in the VCF
- START is the beginning position (inclusive)
- END is the ending position (inclusive)
- The colon and dash are required separators

For example: `chr1:10000-20000`, `X:50000-100000`, `22:30500000-31000000`

### Coordinate System
The position coordinates use the same 1-based system as VCF files, where the first base of a chromosome is position 1.

## Handling Special Cases

### Empty Results
If no variants in the input VCF fall within the specified region, the output will contain only the header lines.

### Non-Existent Chromosome
If the specified chromosome does not exist in the input VCF, the output will contain only header lines.

### Invalid Region Format
If the region is not properly formatted (missing colon, missing dash, or end smaller than start), an error is reported and no filtering is performed.

### Malformed Lines
- Lines with fewer than 2 columns are skipped with a warning
- Lines where the position cannot be parsed as an integer are skipped with a warning
- Data lines encountered before the #CHROM header are skipped with a warning

## Performance Considerations
- Processes the VCF file line by line, with minimal memory requirements
- No preprocessing or indexing is done
- Performance scales linearly with input file size
- More efficient for small regions in large files compared to manually parsing

## Limitations
- Can only extract one continuous region at a time
- No support for multiple regions in a single run
- Cannot extract by gene name or other genomic features directly
- Requires exact chromosome name matching
- No special handling of complex structural variants that may span region boundaries 
