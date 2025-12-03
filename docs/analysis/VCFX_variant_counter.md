# VCFX_variant_counter

## Overview
`VCFX_variant_counter` is a simple utility tool that counts the total number of valid variants (data lines) in a VCF file. It reads input from standard input, processes each line, and outputs the total count of valid variant records.

## Usage
```bash
VCFX_variant_counter [OPTIONS] < input.vcf
```

## Options
| Option | Description |
|--------|-------------|
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |
| `-s`, `--strict` | Fail on any data line with fewer than 8 columns |

## Description
`VCFX_variant_counter` processes a VCF file by:

1. Reading the VCF file from standard input
2. Ignoring all header lines (those starting with `#`)
3. For each data line:
   - Checking if it has at least 8 columns (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)
   - If it has 8 or more columns, counting it as a valid variant
   - If it has fewer than 8 columns:
     - In strict mode: exiting with an error
     - In normal mode: skipping the line with a warning
4. Finally, printing the total count of valid variants

This tool is useful for quickly determining the number of variants in a VCF file, which can be helpful for quality control, workflow validation, or simply getting an overview of a dataset's size.

## VCF Format Requirements
The tool assumes a standard VCF format where:
- Header lines start with `#`
- Data lines have at least 8 tab-separated columns (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)
- Optional FORMAT and sample columns may follow

## Examples

### Basic Usage
Count the variants in a VCF file:
```bash
VCFX_variant_counter < input.vcf
```
Output:
```
Total Variants: 1234
```

### Using Strict Mode
Count variants with strict validation (will fail on malformed lines):
```bash
VCFX_variant_counter --strict < input.vcf
```

### Using in a Pipeline
Count variants after filtering:
```bash
cat input.vcf | grep -v "FILTER=FAIL" | VCFX_variant_counter
```

## Error Handling

### Invalid Lines
- In normal mode (default):
  - Lines with fewer than 8 columns are skipped
  - A warning is printed to standard error for each skipped line
  - The count continues with valid lines

- In strict mode (`--strict`):
  - If any line has fewer than 8 columns, the program exits with an error
  - The error message includes the line number
  - The exit code is 1 to indicate failure

### Empty Lines
Empty lines are ignored and not counted.

## Performance Considerations
- The tool processes the VCF file line by line, with minimal memory requirements
- It performs only basic parsing and doesn't validate the content of each field
- Performance scales linearly with the number of lines in the input file
- No external dependencies or reference files are required

## Limitations
- The tool only checks the number of columns, not their content
- It doesn't validate if the VCF follows the specification for field formats
- No specific handling for compressed files (use external tools like zcat)
- No detailed reporting (e.g., breakdown by chromosome or variant type)
- Cannot handle VCF files with non-standard line endings 
