# VCFX_record_filter

## Overview
`VCFX_record_filter` filters VCF files based on flexible criteria applied to standard fields (POS, QUAL, FILTER) and INFO fields. It allows for complex filtering with multiple conditions using AND/OR logic.

## Usage
```bash
VCFX_record_filter --filter "CRITERIA" [OPTIONS] < input.vcf > filtered.vcf
```

## Options
| Option | Description |
|--------|-------------|
| `-f`, `--filter <CRITERIA>` | Required. One or more filtering criteria separated by semicolons (e.g., `"POS>10000;QUAL>=30;AF<0.05"`) |
| `-l`, `--logic <and\|or>` | Logic for combining multiple criteria: `and` (default) requires all criteria to pass, `or` requires any criterion to pass |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description
`VCFX_record_filter` evaluates each variant in a VCF file against specified criteria and outputs only variants that satisfy these criteria. The tool:

1. Reads a VCF file line by line from standard input
2. Passes all header lines (starting with `#`) unchanged to the output
3. For each data line, evaluates it against the specified criteria
4. If the variant satisfies the criteria, writes it to standard output
5. If the variant fails the criteria, discards it

Criteria can be specified for:
- Standard VCF fields: `POS` (numeric), `QUAL` (numeric), `FILTER` (string)
- Any key in the INFO column (automatically detected as numeric or string)

Each criterion must use one of the following operators:
- Numeric comparisons: `>`, `>=`, `<`, `<=`, `==`, `!=`
- String comparisons: `==`, `!=` (equality and inequality only)

## Output Format
The output is a standard VCF file with the same format as the input, but containing only variants that meet the specified filtering criteria. All header lines are preserved.

## Examples

### Basic Filtering
Filter variants by position:
```bash
VCFX_record_filter --filter "POS>10000" < input.vcf > filtered.vcf
```

### Quality Filtering
Filter variants with QUAL score at least 30:
```bash
VCFX_record_filter --filter "QUAL>=30" < input.vcf > filtered.vcf
```

### Multiple Criteria with AND Logic
Keep only variants that pass all criteria (default AND logic):
```bash
VCFX_record_filter --filter "POS>=1000;FILTER==PASS;DP>10" < input.vcf > filtered.vcf
```

### Multiple Criteria with OR Logic
Keep variants that pass any of the criteria:
```bash
VCFX_record_filter --filter "AF>0.1;DP>100" --logic or < input.vcf > filtered.vcf
```

### Filtering on INFO Fields
Filter based on allele frequency and depth:
```bash
VCFX_record_filter --filter "AF<0.01;DP>=50" < input.vcf > rare_variants.vcf
```

### String Comparison
Filter variants by FILTER status:
```bash
VCFX_record_filter --filter "FILTER==PASS" < input.vcf > passing_variants.vcf
```

## Criterion Parsing

### Field Types
The tool automatically determines field types:
- `POS` and `QUAL` are always treated as numeric fields
- `FILTER` is always treated as a string field
- INFO fields are parsed as numeric if possible, otherwise as strings

### Numeric Values
For numeric fields, the value specified in the criterion is converted to a double and compared using the specified operator:
```
QUAL>=30  # Passes if QUAL is at least 30
DP>10     # Passes if DP INFO field is greater than 10
```

### String Values
For string fields, only equality (`==`) and inequality (`!=`) operators are supported:
```
FILTER==PASS     # Passes if FILTER is exactly "PASS"
SVTYPE!=DEL      # Passes if SVTYPE INFO field is not "DEL"
```

## Handling Special Cases

### Missing Values
- Missing `QUAL` values (`.`) are treated as 0.0
- Missing INFO fields cause the criterion to fail
- Empty fields are handled appropriately

### Multiple Criteria
- With AND logic, a variant must pass ALL criteria to be included
- With OR logic, a variant passes if ANY criterion is satisfied

### Malformed Lines
- Lines with fewer than 8 columns are skipped
- Data lines before the `#CHROM` header are skipped with a warning

## Performance Considerations
- The tool processes VCF files line by line, requiring minimal memory
- Each line is evaluated independently, allowing for efficient processing
- For large files with many criteria, using AND logic can be more efficient as it can short-circuit on the first failing criterion

## Limitations
- String fields only support equality and inequality comparisons, not substring or pattern matching
- No built-in support for sample genotype filtering (focuses on variant-level data)
- Cannot filter based on the number of samples with a particular genotype
- No support for parentheses or complex boolean expressions beyond simple AND/OR logic
- INFO flags (without values) are treated as having a value of 1.0 for numeric comparisons 
