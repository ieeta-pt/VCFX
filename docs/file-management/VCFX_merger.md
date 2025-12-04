# VCFX_merger

## Overview

`VCFX_merger` is a tool for merging multiple VCF files by variant position. It combines multiple VCF files while maintaining proper sorting by chromosome and position, and preserves the VCF header information from the first input file.

**New in v1.1**: The tool now supports **streaming K-way merge** for handling pre-sorted files with O(num_files) memory usage, enabling merging of 50GB+ files without loading them entirely into memory.

## Usage

```bash
VCFX_merger --merge file1.vcf,file2.vcf,... [options] > merged.vcf
```

## Options

| Option | Description |
|--------|-------------|
| `-m, --merge` | Comma-separated list of VCF files to merge |
| `-s, --assume-sorted` | Assume input files are already sorted by (CHROM, POS). Enables streaming K-way merge with O(num_files) memory for large files. |
| `-n, --natural-chr` | Use natural chromosome ordering (chr1 < chr2 < chr10) instead of lexicographic ordering |
| `-h, --help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

`VCFX_merger` reads multiple VCF files and combines them into a single output file while:

1. Preserving the VCF header information from the first input file
2. Sorting all variants by chromosome and position
3. Maintaining the original VCF format and field structure
4. Handling multiple input files efficiently

The tool offers two merge strategies:

### Default Mode (In-Memory)
- Loads all variants from all files into memory
- Sorts the combined variants
- Suitable for small to medium files
- Works with unsorted input files

### Streaming Mode (`--assume-sorted`)
- Performs K-way merge using a min-heap priority queue
- Only keeps one line per input file in memory
- Memory usage: O(num_files) regardless of file sizes
- **Requires pre-sorted input files** for correct results
- Enables merging files larger than available RAM

## Performance Improvements

### Streaming K-way Merge
When using `--assume-sorted`, the tool implements an efficient streaming merge:

1. **One line per file**: Only the current line from each input file is held in memory
2. **Min-heap comparison**: Uses a priority queue to always output the smallest variant
3. **Linear I/O**: Each line is read and written exactly once
4. **Constant memory**: Memory usage is O(num_files), typically ~100KB for 10 files

### Performance Characteristics
| Mode | Memory Usage | Requirements |
|------|--------------|--------------|
| Default | O(total_variants) | None - handles unsorted input |
| `--assume-sorted` | O(num_files) | Input files must be sorted |

### When to Use Streaming Mode
- Input files are already sorted (e.g., output from bcftools, VCFX_sorter)
- Working with very large files (>1GB each)
- Limited available memory
- Merging many files simultaneously

## Input Requirements

- All input files must be in valid VCF format
- Files should have consistent header structures
- The first file's header information will be used in the output
- Files can contain any number of variants
- **For streaming mode**: All input files **MUST** be sorted by (CHROM, POS)

## Output Format

The output is a standard VCF file with:
- Header information from the first input file
- All variants sorted by chromosome and position
- Original VCF format preserved
- Tab-delimited fields maintained

## Examples

### Basic Usage (Default Mode)

Merge two VCF files (unsorted input allowed):

```bash
VCFX_merger --merge sample1.vcf,sample2.vcf > merged.vcf
```

### Multiple Files

Merge three or more VCF files:

```bash
VCFX_merger --merge file1.vcf,file2.vcf,file3.vcf > combined.vcf
```

### Streaming Merge for Large Pre-Sorted Files

Merge large pre-sorted files with minimal memory usage:

```bash
VCFX_merger --merge sorted1.vcf,sorted2.vcf --assume-sorted > merged.vcf
```

### Streaming Merge with Natural Chromosome Ordering

Merge pre-sorted files using natural chromosome order:

```bash
VCFX_merger --merge file1.vcf,file2.vcf --assume-sorted --natural-chr > merged.vcf
```

### Pipeline with Sorting

Pre-sort files then merge for optimal performance:

```bash
# Sort individual files first
VCFX_sorter < sample1.vcf > sample1_sorted.vcf
VCFX_sorter < sample2.vcf > sample2_sorted.vcf

# Then merge with streaming mode
VCFX_merger --merge sample1_sorted.vcf,sample2_sorted.vcf --assume-sorted > merged.vcf
```

### Merging Many Large Files

Efficiently merge 10+ large files:

```bash
VCFX_merger --merge chr1.vcf,chr2.vcf,chr3.vcf,chr4.vcf,chr5.vcf,\
chr6.vcf,chr7.vcf,chr8.vcf,chr9.vcf,chr10.vcf --assume-sorted > combined.vcf
```

### Integration with Other Tools

Merge files and then validate the result:

```bash
VCFX_merger --merge sample1.vcf,sample2.vcf | VCFX_validator > final.vcf
```

## Error Handling

The tool handles various error conditions:

- Missing input files: Reports an error if any specified input file cannot be opened
- Invalid VCF format: Preserves the original format but does not validate it
- Empty files: Handles empty input files gracefully
- Missing --merge argument: Displays help message
- Malformed lines: Lines with unparseable CHROM/POS are skipped with a warning

## Performance Considerations

### Default Mode
- Memory usage scales with the number of variants across all input files
- Processing time depends on the total number of variants and the number of input files
- Files are processed sequentially to collect variants
- Sorting is performed in memory after all variants are collected

### Streaming Mode (`--assume-sorted`)
- Memory usage is constant: O(num_files)
- Processing time is O(n log k) where n is total variants and k is number of files
- Each file is streamed from disk without full loading
- Optimal for merging large pre-sorted files

### Recommended Usage
| Scenario | Recommended Mode |
|----------|-----------------|
| Small files, unsorted | Default mode |
| Large files, pre-sorted | `--assume-sorted` |
| Limited memory | `--assume-sorted` (pre-sort if needed) |
| Many files (>10) | `--assume-sorted` |

## Limitations

- Only supports standard VCF format files
- Does not perform VCF validation (use VCFX_validator for validation)
- Preserves only the header information from the first input file
- Requires all input files to have consistent field structures
- **Streaming mode requires pre-sorted input** - incorrect results if inputs are not sorted

## Common Use Cases

1. Combining multiple sample VCFs into a single file
2. Merging region-specific VCF files
3. Combining results from different variant callers
4. Creating a unified VCF file from multiple analysis runs
5. **Merging chromosome-split VCF files from large-scale projects**
6. **Combining VCF outputs from parallel processing pipelines**

## Best Practices

1. Validate input files before merging
2. Use consistent VCF versions across input files
3. For large files, pre-sort with `VCFX_sorter` then use `--assume-sorted`
4. Verify the output with VCFX_validator after merging
5. Use streaming mode when merging files larger than available RAM
6. Consider using natural chromosome ordering (`-n`) for human genome data

## Backward Compatibility

The tool is fully backward compatible:
- Default behavior (without `--assume-sorted`) works exactly as before
- Existing scripts and pipelines continue to function unchanged
- The `--assume-sorted` option is purely opt-in for performance optimization
