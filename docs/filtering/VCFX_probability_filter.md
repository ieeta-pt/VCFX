# VCFX_probability_filter

## Overview
`VCFX_probability_filter` filters a VCF file based on specified genotype probability values. It allows you to keep only variants where samples meet certain probability thresholds, using various comparison operators.

## Usage
```bash
VCFX_probability_filter --filter-probability "<CONDITION>" < input.vcf > filtered.vcf
```

## Options
| Option | Description |
|--------|-------------|
| `-f, --filter-probability <condition>` | Specify the probability filter condition (e.g., `GP>0.9`) |
| `-h, --help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description
`VCFX_probability_filter` analyzes the genotype probability fields in the FORMAT column of a VCF file and filters variants based on a user-defined condition. The tool:

1. Reads the VCF file line by line from standard input
2. Parses the specified condition (e.g., `GP>0.9`)
3. For each variant, examines the probability values in the specified field
4. Keeps variants where all samples meet the condition
5. Outputs the filtered variants to standard output

The condition must be specified in the format `FIELD OPERATOR VALUE`, where:
- `FIELD` is a valid field in the FORMAT column (e.g., `GP` for genotype probabilities)
- `OPERATOR` is one of: `>`, `<`, `>=`, `<=`, `==`, `!=`
- `VALUE` is a numeric threshold

## Output Format
The output is a standard VCF file with the same format as the input, but containing only variants that meet the specified probability condition. All header lines are preserved.

## Examples

### Basic Filtering
Filter for variants where all samples have a GP (genotype probability) value greater than 0.9:
```bash
VCFX_probability_filter --filter-probability "GP>0.9" < input.vcf > filtered.vcf
```

### Using Different Operators
Filter for variants where all samples have a GP value less than 0.1:
```bash
VCFX_probability_filter --filter-probability "GP<0.1" < input.vcf > filtered.vcf
```

### Exact Match
Filter for variants where all samples have a GP value exactly equal to 0.9:
```bash
VCFX_probability_filter --filter-probability "GP==0.9" < input.vcf > filtered.vcf
```

### Different Probability Fields
Filter using a different probability field (e.g., `PP` for posterior probability):
```bash
VCFX_probability_filter --filter-probability "PP>=0.8" < input.vcf > filtered.vcf
```

### In a Pipeline
Use as part of a processing pipeline:
```bash
cat input.vcf | VCFX_probability_filter --filter-probability "GP>0.95" | other_vcf_tool > output.vcf
```

## Probability Value Handling

### Field Parsing
The tool locates the specified probability field (e.g., `GP`) in the FORMAT column of each variant and extracts the corresponding values for each sample.

### Value Comparison
The extracted probability values are compared to the specified threshold using the given operator. A variant passes the filter only if all samples meet the condition.

### Multiple Probability Values
If a field contains multiple values (e.g., `GP` often contains three values for a biallelic variant), the filter is applied to the first value that can be successfully extracted and converted to a number.

## Handling Special Cases

### Missing Values
If a sample has a missing value (`.`) for the specified probability field, the variant is filtered out.

### Malformed Values
If a probability value cannot be converted to a number, the variant is filtered out with a warning message.

### Missing Fields
If the specified field is not found in the FORMAT column, the tool reports an error and exits.

## Performance Considerations
- The tool processes the VCF file line by line, requiring minimal memory.
- No sorting or indexing is performed, preserving the original order of variants.
- Time complexity is linear with respect to the size of the input file.

## Limitations
- The tool requires that the specified field exists in the FORMAT column of the VCF file.
- It only supports standard comparison operators and cannot handle complex conditions or multiple conditions.
- The tool does not handle cases where different samples might need different thresholds.
- Missing or malformed values cause the variant to be filtered out, which might not be desired in all cases. 
