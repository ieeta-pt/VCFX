# VCFX_format_converter

## Overview
`VCFX_format_converter` is a tool for converting VCF (Variant Call Format) files to other common file formats. It currently supports conversion to BED (Browser Extensible Data) and CSV (Comma-Separated Values) formats.

## Usage
```bash
VCFX_format_converter [OPTIONS] < input.vcf > output.file
```

## Options
| Option | Description |
|--------|-------------|
| `--to-bed` | Convert the input VCF file to BED format |
| `--to-csv` | Convert the input VCF file to CSV format |
| `--help`, `-h` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description
`VCFX_format_converter` reads a VCF file from standard input and converts it to the specified output format. The tool:

1. Parses command-line arguments to determine the desired output format
2. Reads the VCF file line by line
3. Skips header lines (starting with #) in the output
4. Converts each variant record according to the selected format
5. Writes the converted data to standard output

## Conversion Details

### VCF to BED Conversion
BED format requires at least 3 columns:
1. `chrom` - The chromosome name
2. `start` - 0-based start position (VCF positions are 1-based)
3. `end` - 0-based end position (exclusive)

The converter implements the following mapping:
- `chrom`: Direct copy from VCF's CHROM column
- `start`: VCF POS - 1 (clamped to ≥ 0)
- `end`: start + length of REF allele
- `name`: Direct copy from VCF's ID column

For example, a VCF record `chr1 10000 rs123 A G ...` becomes a BED record `chr1 9999 10000 rs123`.

### VCF to CSV Conversion
The CSV conversion preserves all columns from the VCF file but changes the delimiter from tabs to commas. The tool:

1. Includes a header row from the #CHROM line (after removing the # character)
2. Properly escapes fields containing commas or quotes according to CSV standards:
   - Fields containing commas or quotes are enclosed in double quotes
   - Double quotes within fields are doubled (e.g., `"` becomes `""`)
3. Preserves all columns and their order from the original VCF

## Examples

### Converting to BED Format
```bash
VCFX_format_converter --to-bed < input.vcf > output.bed
```

This command converts each variant record in `input.vcf` to BED format, writing the results to `output.bed`.

### Converting to CSV Format
```bash
VCFX_format_converter --to-csv < input.vcf > output.csv
```

This command converts the tab-delimited `input.vcf` to comma-separated format in `output.csv`, preserving all columns.

## Handling Special Cases

### Multi-allelic Variants
- In BED format, multi-allelic variants are represented by a single interval based on the REF allele length
- In CSV format, multi-allelic variants are preserved as in the original VCF, with commas in the ALT field properly escaped

### Structural Variants
- For BED conversion, the REF allele length is used to calculate the end position
- If a structural variant includes an END tag in the INFO field, this is not currently used for BED conversion

### Special Characters in CSV
Fields containing special characters are handled according to CSV standards:
- Commas: Field is enclosed in double quotes
- Double quotes: Field is enclosed in double quotes and internal quotes are doubled
- Newlines: Preserved within quoted fields
- Fields with ID values containing commas are properly quoted

### Malformed VCF Files
The tool attempts to handle malformed input gracefully:
- Lines with too few columns are skipped
- Invalid position values are skipped
- Non-numeric QUAL values are preserved as-is in CSV output and handled appropriately for BED

### Empty Files or Headers-Only Files
- For files containing only headers (no variant records), the tool produces:
  - An empty file for BED output
  - A header row only for CSV output

## Performance Considerations
- The tool processes VCF files line by line, with minimal memory requirements
- Performance scales linearly with input file size
- No indexing is performed, allowing efficient streaming processing

## Limitations
- The BED conversion uses a simple interval representation based on REF allele length
- Structural variants with END tags are not specially handled for BED intervals
- No specific handling for insertions (which technically have zero length in reference coordinates)
- Does not create or use specialized indices
- Does not support other output formats such as GFF, TXT, or JSON
- Cannot perform filtering operations during conversion 
