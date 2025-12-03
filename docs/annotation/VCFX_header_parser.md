# VCFX_header_parser

## Overview

VCFX_header_parser is a simple utility that extracts and displays all header lines from a VCF file. This tool makes it easy to examine metadata and structural information without processing the variant data.

## Usage

```bash
VCFX_header_parser [OPTIONS] < input.vcf > header.txt
```

## Options

| Option | Description |
|--------|-------------|
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_header_parser reads a VCF file from standard input and outputs only the header lines (lines starting with "#"). The tool:

1. Reads the VCF file line by line
2. Extracts all lines beginning with "#", which include:
   - VCF version information (`##fileformat=VCFv4.2`)
   - Reference genome information (`##reference=file:///path/to/reference.fa`)
   - Contig definitions (`##contig=<ID=chr1,length=248956422>`)
   - INFO field definitions (`##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">`)
   - FILTER definitions (`##FILTER=<ID=PASS,Description="All filters passed">`)
   - FORMAT definitions (`##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`)
   - Sample column header line (`#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE1  SAMPLE2`)
3. Stops reading when it encounters the first non-header line (any line not starting with "#")
4. Outputs all collected header lines to standard output

This tool is useful for:
- Examining metadata without processing large variant datasets
- Extracting sample names from a VCF file
- Checking VCF file structure and compliance with specifications
- Creating header templates for new VCF files
- Documenting file provenance and contents

## Output Format

The output consists of all header lines from the input VCF file, in the same order they appeared in the original file:

```
##fileformat=VCFv4.2
##source=VCFX
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=248956422>
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE1  SAMPLE2
```

## Examples

### Basic Usage

```bash
./VCFX_header_parser < input.vcf > header.txt
```

### Extracting Sample Names

```bash
# Extract the sample names (all columns after FORMAT in the #CHROM line)
./VCFX_header_parser < input.vcf | grep "^#CHROM" | cut -f10- > sample_names.txt
```

### Counting Contigs

```bash
# Count the number of contigs defined in the header
./VCFX_header_parser < input.vcf | grep "##contig" | wc -l
```

### Verifying VCF Version

```bash
# Check the VCF file format version
./VCFX_header_parser < input.vcf | grep "##fileformat" | cut -d= -f2
```

## Handling Special Cases

The tool implements simple strategies for handling edge cases:

1. **Empty files**: If the input file is empty, no output is produced
2. **Files without headers**: If the file has no header lines, no output is produced
3. **Malformed headers**: All lines starting with "#" are considered header lines, even if they don't follow VCF specifications
4. **Line endings**: LF and CRLF line endings are handled correctly
5. **Partial headers**: If the file ends in the middle of the header section, all header lines up to that point are output

## Performance

VCFX_header_parser is designed for simplicity and efficiency:

1. Processes input line-by-line without loading the entire file into memory
2. Stops processing as soon as it encounters the first non-header line
3. Highly efficient for large VCF files where headers constitute a small portion of the total file size
4. Minimal memory footprint since only the current line being processed is stored in memory

## Limitations

1. No validation of header syntax or compliance with VCF specifications
2. Cannot modify or filter specific header lines
3. No ability to sort or organize header lines
4. No special handling for duplicate header entries
5. Cannot add or merge headers from multiple files 
