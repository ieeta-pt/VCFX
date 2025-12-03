# VCFX_population_filter

## Overview
`VCFX_population_filter` is a utility tool for subsetting VCF files to include only samples belonging to a specified population. It filters out samples that don't belong to the chosen population group while preserving the variant data and format information.

## Usage
```bash
VCFX_population_filter [OPTIONS] < input.vcf > output.vcf
```

## Options
| Option | Description |
|--------|-------------|
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |
| `-p`, `--population <TAG>` | **Required**: Population tag to keep (e.g., 'EUR', 'AFR', 'EAS') |
| `-m`, `--pop-map <FILE>` | **Required**: Tab-delimited file mapping sample names to populations |

## Description
`VCFX_population_filter` processes a VCF file to create a population-specific subset by:

1. Reading a population map file that associates each sample with a population
2. Identifying samples that belong to the specified population
3. Reading the VCF file from standard input
4. Preserving all meta-information lines (starting with '##') without modification
5. Modifying the #CHROM header line to include only samples from the specified population
6. For each data line:
   - Keeping the first 9 mandatory columns (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)
   - Including only genotype columns for samples from the specified population
7. Writing the filtered VCF to standard output

This tool is particularly useful for:
- Creating population-specific VCF files for population genetics analysis
- Reducing file size by excluding irrelevant samples
- Focusing on specific ancestral groups for targeted studies
- Preparing data for population-stratified association studies

## Population Map Format

The population map file should be a simple tab-delimited text file with:
- Each line containing a sample name and its population designation
- The first column containing the exact sample name as it appears in the VCF header
- The second column containing the population identifier
- No header row (just data rows)

Example population map file:
```
SAMPLE1  EUR
SAMPLE2  EUR
SAMPLE3  AFR
SAMPLE4  AFR
SAMPLE5  EAS
```

## Examples

### Basic Usage
Filter a VCF file to include only European (EUR) samples:
```bash
VCFX_population_filter --population EUR --pop-map population_map.txt < input.vcf > eur_only.vcf
```

### Different Populations
Create separate files for different populations:
```bash
VCFX_population_filter --population AFR --pop-map population_map.txt < input.vcf > afr_only.vcf
VCFX_population_filter --population EAS --pop-map population_map.txt < input.vcf > eas_only.vcf
```

## Example Transformations

### Input VCF (with multiple populations)
```
##fileformat=VCFv4.2
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE1_EUR  SAMPLE2_EUR  SAMPLE3_AFR  SAMPLE4_AFR  SAMPLE5_EAS
1  100  rs123  A  T  50  PASS  AF=0.1  GT:DP  0|0:30  0|1:25  1|1:20  0|1:22  0|0:18
```

### Output VCF (filtered for EUR)
```
##fileformat=VCFv4.2
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE1_EUR  SAMPLE2_EUR
1  100  rs123  A  T  50  PASS  AF=0.1  GT:DP  0|0:30  0|1:25
```

## Special Case Handling

### Missing Samples
- If no samples match the specified population, a warning is issued
- The output VCF will contain only the first 9 mandatory columns without any sample data

### Malformed Lines
- Lines with fewer than 9 columns are skipped with a warning
- Data lines encountered before the #CHROM header are skipped with a warning

### Missing #CHROM Header
- If no #CHROM header is found in the file, an error is reported

### Invalid Population Map
- If the population map file cannot be opened, an error is reported
- Lines in the population map that don't follow the expected format are skipped with a warning

## Performance Considerations
- The tool processes the VCF file line by line, with minimal memory requirements
- Memory usage is primarily determined by the number of samples in the VCF file
- Performance scales linearly with the size of the input file
- No external dependencies or reference files are required beyond the population map

## Limitations
- Cannot filter variants based on population-specific criteria
- Does not update INFO fields like AC/AN to reflect the reduced sample set
- No support for more complex population filtering (e.g., including multiple populations)
- Cannot handle compressed (gzipped) VCF files directly
- Does not validate sample consistency between the VCF and the population map 
