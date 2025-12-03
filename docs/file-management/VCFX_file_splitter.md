# VCFX_file_splitter

## Overview

VCFX_file_splitter divides a VCF file into multiple smaller files based on chromosome, creating separate output files for each chromosome present in the input.

## Usage

```bash
VCFX_file_splitter [OPTIONS] < input.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `-p`, `--prefix <PREFIX>` | Output file prefix (default: "split") |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_file_splitter reads a VCF file and separates its contents into multiple files, with one file per chromosome. The tool:

1. Reads a VCF file from standard input
2. Extracts the chromosome (CHROM) information from each variant line
3. Creates a separate output file for each unique chromosome encountered
4. Writes all header lines to each output file
5. Distributes variant records to the appropriate chromosome file
6. Produces output files named using the pattern `<PREFIX>_<CHROM>.vcf`

This tool is useful for:
- Parallelizing variant processing by chromosome
- Reducing memory requirements when handling large VCF files
- Organizing variant data by chromosome for downstream analysis
- Creating chromosome-specific VCF files for targeted analysis
- Preparing data for tools that work on individual chromosomes

## Output Format

The output consists of multiple VCF files, one for each chromosome in the input. Each file contains:
- All original header lines from the input VCF
- Only the variant records for the corresponding chromosome
- The same format and structure as the original VCF

Output files are named following the pattern:
```
<PREFIX>_<CHROM>.vcf
```

For example, using the default prefix "split", the tool will generate files like:
- `split_1.vcf` (chromosome 1)
- `split_2.vcf` (chromosome 2)
- `split_X.vcf` (chromosome X)
- etc.

## Examples

### Basic Usage

```bash
./VCFX_file_splitter < input.vcf
```

This will create files like `split_1.vcf`, `split_2.vcf`, etc.

### Custom Prefix

```bash
./VCFX_file_splitter --prefix "chr" < input.vcf
```

This will create files like `chr_1.vcf`, `chr_2.vcf`, etc.

### Processing Multiple Files

```bash
# Split multiple VCF files
for file in *.vcf; do
  output_prefix="${file%.vcf}"
  ./VCFX_file_splitter --prefix "$output_prefix" < "$file"
done
```

## Handling Special Cases

- **Header Lines**: All header lines (starting with #) are included in each output file
- **Additional Headers**: If header lines appear after data lines in the input, they are replicated to all open chromosome files
- **Empty Input**: If the input file is empty or contains only headers, a warning message is displayed
- **Chromosome Naming**: Preserves chromosome names exactly as they appear in the input file, including any prefixes or special characters
- **Malformed Lines**: Lines that can't be parsed for chromosome information are skipped with a warning
- **File Creation Failures**: Reports an error if an output file cannot be created (due to permissions, disk space, etc.)
- **Large Numbers of Chromosomes**: Can handle arbitrarily many chromosomes, creating one file for each

## Performance

The splitter is optimized for efficiency:
- Single-pass processing of the input file
- Streams data directly to output files without storing records in memory
- Uses smart pointers for automatic resource management
- Efficiently handles very large VCF files with minimal memory overhead
- Output files are written incrementally as the input is processed

## Limitations

- Requires sufficient disk space to store all output files
- No built-in compression of output files
- Cannot split by other criteria (e.g., position ranges, sample names)
- Does not check for duplicate variant entries in the input
- No option to merge small chromosomes into a single output file
- Cannot control the order of variants within output files (maintains the order from the input)
- Files with many chromosomes will generate many output files 
