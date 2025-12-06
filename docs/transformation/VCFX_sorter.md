# VCFX_sorter

## Overview
`VCFX_sorter` is a utility tool for sorting VCF files by chromosome and position. It provides two sorting methods: standard lexicographic sorting and natural chromosome sorting, which handles chromosome numbering in a more intuitive way.

**New in v1.1**: The tool now supports **external merge sort** for handling files larger than available RAM, enabling sorting of 50GB+ VCF files with minimal memory usage.

**New in v1.2**: Added **memory-mapped I/O** for file arguments, providing ~40x faster processing for large files compared to stdin.

## Usage
```bash
VCFX_sorter [OPTIONS] [files...] > output.vcf
VCFX_sorter [OPTIONS] < input.vcf > output.vcf
```

## Options
| Option | Description |
|--------|-------------|
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |
| `-n`, `--natural-chr` | Use natural chromosome sorting (chr1 < chr2 < chr10) instead of lexicographic sorting |
| `-m`, `--max-memory <MB>` | Maximum memory for in-memory sorting in megabytes (default: 100MB). Files larger than this threshold use external merge sort with temporary files. |
| `-t`, `--temp-dir <DIR>` | Directory for temporary files when using external merge sort (default: `/tmp`) |

## Description
`VCFX_sorter` processes a VCF file to organize variants in a consistent order by:

1. Reading the VCF file from standard input
2. Preserving all header lines without modification
3. Sorting the data lines by chromosome and position using either in-memory or external merge sort
4. Writing the header lines followed by the sorted data lines to standard output

The tool supports two distinct sorting methods:
- **Lexicographic sorting** (default): Sorts chromosomes alphabetically (chr1, chr10, chr2, ...)
- **Natural sorting**: Sorts chromosomes in numeric order when possible (chr1, chr2, ..., chr10, ...)

This tool is particularly useful for:
- Preparing VCF files for downstream analysis tools that expect sorted input
- Merging multiple VCF files that need consistent ordering
- Improving readability and navigation of VCF files
- Making binary searches possible on VCF data
- **Sorting very large VCF files (50GB+) that exceed available RAM**

## Performance Improvements

### Memory-Mapped I/O (v1.2)
When file arguments are provided instead of stdin, the tool uses memory-mapped I/O (`mmap`) with several optimizations:

1. **Zero-copy Parsing**: Parses chromosome and position directly from mmap'd memory without string allocation
2. **Pre-computed Chromosome IDs**: Converts chromosome strings to numeric IDs once during parsing, eliminating O(n log n) string comparisons during sort
3. **CompactSortKey Structure**: Uses a 20-byte sort key storing only offsets instead of full VCF lines (vs ~9KB per variant with full lines)
4. **1MB Output Buffer**: Batched output reduces syscall overhead

#### Benchmark Results
On a test file with 427K variants and 2,504 samples (1.4GB):
- **stdin processing**: >30 minutes (timeout)
- **mmap file argument**: ~46 seconds
- **Improvement**: ~40x faster with file arguments

### External Merge Sort for Large Files
The tool now implements **external merge sort** for handling files larger than available RAM:

1. **Chunked Processing**: Data is split into memory-bounded chunks
2. **Sorted Chunk Files**: Each chunk is sorted in memory and written to a temporary file
3. **K-way Merge**: All sorted chunks are merged using a min-heap priority queue
4. **Automatic Cleanup**: Temporary files are automatically deleted after merging

This enables sorting of files many times larger than available memory while maintaining O(n log n) performance.

### Memory Management
- **Configurable Memory Limit**: Use `-m` to control maximum memory usage (default: 100MB)
- **Streaming Processing**: Only holds one chunk in memory at a time during external sort
- **Efficient Merging**: K-way merge uses O(num_chunks) memory regardless of file size

### Performance Characteristics
| File Size | Memory Usage | Algorithm |
|-----------|--------------|-----------|
| < 100MB (default) | O(n) | In-memory sort |
| > 100MB | O(chunk_size) | External merge sort |
| 50GB file | ~512MB peak | External merge sort |

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

### Basic Lexicographic Sorting (File Argument - Fastest)
Sort a VCF file using standard lexicographic chromosome ordering:
```bash
VCFX_sorter unsorted.vcf > sorted.vcf
```

### Basic Lexicographic Sorting (Stdin)
Sort a VCF file from standard input:
```bash
VCFX_sorter < unsorted.vcf > sorted.vcf
```

### Natural Chromosome Sorting
Sort a VCF file using natural chromosome ordering:
```bash
VCFX_sorter --natural-chr unsorted.vcf > sorted.vcf
```

### Large File Sorting (External Merge Sort)
Sort a large VCF file using external merge sort with a 500MB memory limit:
```bash
VCFX_sorter -m 500 < huge_file.vcf > sorted.vcf
```

### Custom Temporary Directory
Use a specific directory for temporary files (useful for faster storage):
```bash
VCFX_sorter -m 500 -t /fast/scratch < huge_file.vcf > sorted.vcf
```

### Large File with Natural Sorting
Sort a large file with natural chromosome ordering:
```bash
VCFX_sorter --natural-chr -m 1000 < whole_genome.vcf > sorted.vcf
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
- Lines with fewer than 2 tab-separated columns are skipped with a warning
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

### Memory Usage
- **Small files**: The tool reads the entire VCF file into memory for optimal performance
- **Large files**: External merge sort keeps memory usage bounded regardless of file size
- Use `-m` to tune the memory/performance tradeoff:
  - Higher values = fewer temporary files, faster processing
  - Lower values = less memory usage, more disk I/O

### Temporary Storage
- External merge sort creates temporary files in the specified temp directory
- Ensure sufficient disk space for approximately 1-2x the input file size
- SSD or fast storage for temp files improves performance significantly

### Processing Time
- In-memory: O(n log n) where n is the number of variants
- External merge: O(n log n) with additional disk I/O overhead
- K-way merge phase is efficient with O(n log k) comparisons where k is chunk count

### Recommended Settings
| Scenario | Recommended Options |
|----------|-------------------|
| Small files (<100MB) | Default settings |
| Medium files (100MB-1GB) | `-m 500` |
| Large files (1GB-10GB) | `-m 1000 -t /fast/ssd` |
| Very large files (>10GB) | `-m 2000 -t /fast/ssd` |

## Limitations
- Cannot sort by other fields besides chromosome and position
- Does not validate VCF format beyond basic column counting
- No handling of compressed (gzipped) VCF files directly
- Cannot maintain the original order of variants at the same chromosome and position

## Backward Compatibility
The tool is fully backward compatible. Existing commands work identically:
- Default behavior uses in-memory sorting for small files
- The external merge sort activates automatically for large files
- All existing command-line options remain unchanged
