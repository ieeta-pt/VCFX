# VCFX_gl_filter

## Overview
`VCFX_gl_filter` filters VCF records based on numeric genotype-likelihood fields in the FORMAT column, such as genotype quality (GQ), read depth (DP), or phred-scaled likelihoods (PL). This tool helps focus analysis on variants with sufficient genotype quality or other sample-level metrics.

## Usage
```bash
VCFX_gl_filter --filter "<CONDITION>" [--mode <any|all>] < input.vcf > filtered.vcf
```

## Options
| Option | Description |
|--------|-------------|
| `-f`, `--filter <CONDITION>` | Required. Filter condition (e.g., `GQ>20`, `DP>=10`, `PL<50`) |
| `-m`, `--mode <any\|all>` | Optional. Determines if all samples must pass the condition (`all`, default) or at least one sample must pass (`any`) |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description
`VCFX_gl_filter` examines numeric fields in the FORMAT column of a VCF file and filters variant records based on whether the samples satisfy the specified condition. The tool:

1. Parses the filter condition into field name, operator, and threshold value
2. Locates the specified field in the FORMAT column
3. For each variant record, evaluates sample values against the condition
4. Applies the filtering logic based on the specified mode:
   - In `all` mode (default): keeps variants where ALL samples pass the condition
   - In `any` mode: keeps variants where AT LEAST ONE sample passes the condition
5. Outputs passing records to standard output

This tool is particularly useful for:
- Removing variants with low genotype quality
- Filtering based on read depth or coverage
- Filtering on phred-scaled likelihoods or other numeric likelihood measures
- Applying consistent quality thresholds across samples

## Output Format
The output is a standard VCF file containing:
- All original header lines from the input VCF
- Only those variant records where samples satisfy the specified condition according to the mode
- No modification to the content or format of the retained lines

## Examples

### Basic Usage with Default Mode
Filter variants where all samples have genotype quality (GQ) above 20:
```bash
VCFX_gl_filter --filter "GQ>20" < input.vcf > high_quality.vcf
```

### Using 'Any' Mode
Filter variants where at least one sample has read depth (DP) of 30 or higher:
```bash
VCFX_gl_filter --filter "DP>=30" --mode any < input.vcf > high_depth.vcf
```

### Exact Value Matching
Filter variants where all samples have an exact phred-scaled likelihood value:
```bash
VCFX_gl_filter --filter "PL==50" < input.vcf > specific_pl.vcf
```

### Negative Filtering with Not Equal
Filter variants where all samples have a non-zero genotype quality:
```bash
VCFX_gl_filter --filter "GQ!=0" < input.vcf > non_zero_gq.vcf
```

### In a Pipeline
Use with other VCFX tools in a pipeline:
```bash
cat input.vcf | VCFX_gl_filter --filter "GQ>30" | VCFX_record_filter --filter "QUAL>40" > high_quality_variants.vcf
```

## Filter Condition Syntax

### Format
The filter condition must follow this syntax:
```
FIELD OPERATOR VALUE
```
Where:
- `FIELD`: Any numeric field from the FORMAT column (e.g., GQ, DP, PL)
- `OPERATOR`: One of `>`, `<`, `>=`, `<=`, `==`, `!=`
- `VALUE`: A numeric threshold (integer or decimal)

Examples of valid conditions:
- `GQ>20`: Genotype quality greater than 20
- `DP>=10.5`: Read depth greater than or equal to 10.5
- `PL<30`: Phred-scaled likelihood less than 30
- `GL!=0`: Genotype likelihood not equal to 0

### Comparison Operators
The tool supports the following comparison operators:
- `>`: Greater than
- `<`: Less than
- `>=`: Greater than or equal to
- `<=`: Less than or equal to
- `==`: Equal to
- `!=`: Not equal to

## Handling Special Cases

### Missing Fields
If the specified field is not found in the FORMAT column:
- In `all` mode: The variant is filtered out
- In `any` mode: The variant is filtered out

### Missing Values
If a sample has a missing value (`.`) for the specified field:
- In `all` mode: The variant is filtered out
- In `any` mode: The sample is treated as not passing, but the variant may be kept if other samples pass

### Empty Values
Empty values are treated similarly to missing values:
- In `all` mode: The variant is filtered out
- In `any` mode: The sample is treated as not passing

### Non-Numeric Values
If a field value cannot be converted to a number:
- In `all` mode: The variant is filtered out
- In `any` mode: The sample is treated as not passing

### Malformed VCF
- For lines with insufficient fields (less than the standard VCF format requires), the tool produces a warning and skips the line
- For data lines before the header (#CHROM) line, the tool produces a warning and skips the line

## Performance Considerations

### Optimized Implementation
The tool uses several performance optimizations for fast processing of large VCF files with many samples:

1. **Zero-allocation parsing**: Sample genotype fields are parsed using raw pointer arithmetic instead of creating intermediate string objects, eliminating millions of heap allocations on large files.

2. **Streaming sample processing**: Samples are processed one-by-one without storing positions, allowing the tool to handle files with thousands of samples efficiently.

3. **Fast numeric parsing**: A custom inline parser replaces `std::stod` to avoid exception handling overhead and unnecessary string conversions.

4. **Output buffering**: A 1MB output buffer reduces system call overhead for writing filtered records.

5. **Pre-computed operator type**: The comparison operator is converted to an enum at startup, enabling fast switch-based comparison in the hot loop.

### Complexity
- Linear time complexity with respect to file size: O(variants × samples)
- Constant memory usage regardless of sample count
- No preprocessing or indexing of the VCF file is required

### Scalability
The streaming architecture allows processing of VCF files with:
- Unlimited number of samples (tested with 5000+ samples)
- Large file sizes (multi-GB files)
- Minimal memory footprint

## Limitations
- Only works with numeric fields in the FORMAT column
- No support for filtering on string-valued FORMAT fields
- Cannot apply different conditions to different samples
- Cannot combine multiple conditions in a single filter
- No special handling for multi-allelic sites with multiple values per field
- Only evaluates the first value when a field contains multiple values 
