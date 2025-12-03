# VCFX_sorter

## Overview
`VCFX_sorter` is a utility tool for sorting VCF files by chromosome and position. It provides two sorting methods: standard lexicographic sorting and natural chromosome sorting, which handles chromosome numbering in a more intuitive way.

## Usage
```bash
VCFX_sorter [OPTIONS] < input.vcf > output.vcf
```

## Options
| Option | Description |
|--------|-------------|
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |
| `-n`, `--natural-chr` | Use natural chromosome sorting (chr1 < chr2 < chr10) instead of lexicographic sorting |

## Description
`VCFX_sorter` processes a VCF file to organize variants in a consistent order by:

1. Reading the VCF file from standard input
2. Preserving all header lines without modification
3. Loading all data lines into memory
4. Sorting the data lines by chromosome and position
5. Writing the header lines followed by the sorted data lines to standard output

The tool supports two distinct sorting methods:
- **Lexicographic sorting** (default): Sorts chromosomes alphabetically (chr1, chr10, chr2, ...)
- **Natural sorting**: Sorts chromosomes in numeric order when possible (chr1, chr2, ..., chr10, ...)

This tool is particularly useful for:
- Preparing VCF files for downstream analysis tools that expect sorted input
- Merging multiple VCF files that need consistent ordering
- Improving readability and navigation of VCF files
- Making binary searches possible on VCF data

## Sorting Details

### Lexicographic Sorting
In the default lexicographic mode:
- Chromosomes are compared as strings (e.g., 'chr2' comes after 'chr10')
- Positions are compared numerically within the same chromosome

### Natural Chromosome Sorting
When the `--natural-chr` option is used:
1. The "chr" prefix (case-insensitive) is identified and removed
2. Any leading digits are parsed as a number
3. Remaining characters are treated as a suffix
4. Sorting precedence:
   - First by chromosome prefix (if different)
   - Then by numeric part (if both have numbers)
   - Then by suffix (if both have the same number)
   - Finally by position

This results in more intuitive ordering where chr1 < chr2 < chr10, instead of chr1 < chr10 < chr2.

## Examples

### Basic Lexicographic Sorting
Sort a VCF file using standard lexicographic chromosome ordering:
```bash
VCFX_sorter < unsorted.vcf > sorted.vcf
```

### Natural Chromosome Sorting
Sort a VCF file using natural chromosome ordering:
```bash
VCFX_sorter --natural-chr < unsorted.vcf > sorted.vcf
```

## Example Transformations

### Lexicographic Sorting
```
Before:
chr2  1000  .  A  T  .  PASS  .
chr1  2000  .  G  C  .  PASS  .
chr10 1500  .  T  A  .  PASS  .

After:
chr1  2000  .  G  C  .  PASS  .
chr10 1500  .  T  A  .  PASS  .
chr2  1000  .  A  T  .  PASS  .
```

### Natural Chromosome Sorting
```
Before:
chr2  1000  .  A  T  .  PASS  .
chr1  2000  .  G  C  .  PASS  .
chr10 1500  .  T  A  .  PASS  .

After:
chr1  2000  .  G  C  .  PASS  .
chr2  1000  .  A  T  .  PASS  .
chr10 1500  .  T  A  .  PASS  .
```

## Handling Special Cases

### Malformed Lines
- Lines with fewer than 8 columns are skipped with a warning
- Lines with an invalid position value are skipped with a warning

### Empty Input
- If no input is provided, the help message is displayed

### Missing Header
- If no #CHROM header line is found in the input, a warning is issued but processing continues

### Complex Chromosome Names
- Chromosomes with non-standard naming follow sorting rules based on the selected mode
- Examples of parsing in natural mode:
  - "chr1" → prefix="chr", number=1, suffix=""
  - "chrX" → prefix="chr", number=none, suffix="X"
  - "chr10_alt" → prefix="chr", number=10, suffix="_alt"
  - "scaffold_123" → prefix="", number=none, suffix="scaffold_123"

## Performance Considerations
- The tool reads the entire VCF file into memory before sorting
- Memory usage scales with the number of variants in the input file
- Very large VCF files may require significant memory
- Processing time is dominated by the sorting operation, which is O(n log n)

## Limitations
- No support for on-disk sorting of files too large to fit in memory
- Cannot sort by other fields besides chromosome and position
- Does not validate VCF format beyond basic column counting
- No handling of compressed (gzipped) VCF files directly
- Cannot maintain the original order of variants at the same chromosome and position 
