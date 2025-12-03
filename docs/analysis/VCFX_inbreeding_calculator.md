# VCFX_inbreeding_calculator

## Overview

VCFX_inbreeding_calculator computes the inbreeding coefficient (F) for each sample in a VCF file, providing a measure of homozygosity relative to Hardy-Weinberg equilibrium expectations.

## Usage

```bash
VCFX_inbreeding_calculator [OPTIONS] < input.vcf > output.txt
```

## Options

| Option | Description |
|--------|-------------|
| `--freq-mode` MODE | How to calculate allele frequencies: 'excludeSample' (default) or 'global' |
| `--skip-boundary` | Skip sites with boundary frequencies (p=0 or p=1) |
| `--count-boundary-as-used` | Count boundary sites in usedCount even when skipping them |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description

VCFX_inbreeding_calculator analyzes biallelic variants in a VCF file to calculate the inbreeding coefficient (F) for each sample. The inbreeding coefficient measures the reduction in heterozygosity compared to expectations under Hardy-Weinberg equilibrium.

The tool:
1. Reads the VCF file in a single pass
2. Filters for biallelic variants (ignores sites with multiple ALT alleles)
3. Encodes genotypes as 0 (0/0), 1 (0/1), 2 (1/1), or -1 (missing/invalid)
4. Calculates allele frequencies using the specified method
5. Computes the inbreeding coefficient for each sample

### Frequency Modes

The tool offers two methods for calculating allele frequencies:

- **excludeSample**: Each sample's inbreeding coefficient is calculated using allele frequencies derived from all other samples (excluding itself)
- **global**: A single global allele frequency is calculated using all samples, and this same frequency is used for all samples

### Boundary Handling

For sites where the allele frequency (p) is 0 or 1, you have three options:

1. Use these sites normally (default)
2. Skip boundary sites completely (`--skip-boundary`)
3. Skip boundary sites for calculations but count them as used sites (`--skip-boundary --count-boundary-as-used`)

## Output Format

The output is a tab-delimited text file with the following columns:

```
Sample  InbreedingCoefficient
```

Where:
- Sample is the sample name from the VCF file
- InbreedingCoefficient is the calculated F value, or "NA" if no usable sites were found

## Examples

### Basic Usage (Default Settings)

```bash
./VCFX_inbreeding_calculator < input.vcf > inbreeding_coefficients.txt
```

### Using Global Frequency Mode

```bash
./VCFX_inbreeding_calculator --freq-mode global < input.vcf > global_inbreeding.txt
```

### Skip Boundary Frequencies

```bash
./VCFX_inbreeding_calculator --skip-boundary < input.vcf > non_boundary_inbreeding.txt
```

### Custom Boundary Handling

```bash
./VCFX_inbreeding_calculator --skip-boundary --count-boundary-as-used < input.vcf > custom_boundary.txt
```

## Formula and Calculation

The inbreeding coefficient is calculated as:

F = (O - E) / (T - E)

Where:
- O = Observed homozygosity (count of 0/0 and 1/1 genotypes)
- E = Expected homozygosity under HWE (∑(p²+q²) across sites)
- T = Total number of sites (after filtering)

In boundary cases:
- If total sites used = 0, F = NA
- If expected = observed, F = 0
- If expected = total, F = 1

## Handling Special Cases

- **Single-sample VCFs**: Always produces "NA" as inbreeding requires population context
- **Multi-allelic sites**: Skipped entirely (only biallelic variants are considered)
- **Missing genotypes**: Coded as -1 and excluded from calculations
- **Non-diploid genotypes**: Treated as missing and excluded
- **Boundary frequencies**: Special handling available via command-line options
- **Zero usable sites**: Returns "NA" for the sample
- **Small sample sizes**: May produce unreliable estimates

## Performance

The tool performs a single pass through the VCF file, giving it linear complexity with respect to file size. Memory usage scales with:
- Number of samples in the VCF
- Number of biallelic variants

## Limitations

- Only works with biallelic variants (multiallelic sites are skipped)
- Assumes diploid genotypes
- May produce unexpected results with very small sample sizes
- No built-in filtering for variant quality or other metrics
- No chromosome or region-specific analysis
- Cannot handle populations with substructure (assumes random mating) 
