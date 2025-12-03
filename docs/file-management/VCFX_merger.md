# VCFX_merger

## Overview

`VCFX_merger` is a tool for merging multiple VCF files by variant position. It combines multiple VCF files while maintaining proper sorting by chromosome and position, and preserves the VCF header information from the first input file.

## Usage

```bash
VCFX_merger --merge file1.vcf,file2.vcf,... [options] > merged.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `-m, --merge` | Comma-separated list of VCF files to merge |
| `-h, --help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

`VCFX_merger` reads multiple VCF files and combines them into a single output file while:

1. Preserving the VCF header information from the first input file
2. Sorting all variants by chromosome and position
3. Maintaining the original VCF format and field structure
4. Handling multiple input files efficiently

The tool processes the input files sequentially and merges all variants while ensuring proper sorting. It is particularly useful when you need to combine multiple VCF files from different sources or samples into a single, properly sorted VCF file.

## Input Requirements

- All input files must be in valid VCF format
- Files should have consistent header structures
- The first file's header information will be used in the output
- Files can contain any number of variants

## Output Format

The output is a standard VCF file with:
- Header information from the first input file
- All variants sorted by chromosome and position
- Original VCF format preserved
- Tab-delimited fields maintained

## Examples

### Basic Usage

Merge two VCF files:

```bash
VCFX_merger --merge sample1.vcf,sample2.vcf > merged.vcf
```

### Multiple Files

Merge three or more VCF files:

```bash
VCFX_merger --merge file1.vcf,file2.vcf,file3.vcf > combined.vcf
```

### Integration with Other Tools

Merge files and then process the result:

```bash
VCFX_merger --merge sample1.vcf,sample2.vcf | VCFX_sorter | VCFX_validator > final.vcf
```

## Error Handling

The tool handles various error conditions:

- Missing input files: Reports an error if any specified input file cannot be opened
- Invalid VCF format: Preserves the original format but does not validate it
- Empty files: Handles empty input files gracefully
- Missing --merge argument: Displays help message

## Performance Considerations

- Memory usage scales with the number of variants across all input files
- Processing time depends on the total number of variants and the number of input files
- Files are processed sequentially to minimize memory usage
- Sorting is performed in memory after all variants are collected

## Limitations

- Only supports standard VCF format files
- Does not perform VCF validation (use VCFX_validator for validation)
- Preserves only the header information from the first input file
- Requires all input files to have consistent field structures

## Common Use Cases

1. Combining multiple sample VCFs into a single file
2. Merging region-specific VCF files
3. Combining results from different variant callers
4. Creating a unified VCF file from multiple analysis runs

## Best Practices

1. Validate input files before merging
2. Use consistent VCF versions across input files
3. Consider file sizes and available memory when merging many files
4. Verify the output with VCFX_validator after merging
5. Use VCFX_sorter if additional sorting is needed 
