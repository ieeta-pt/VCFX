# VCFX_alignment_checker

## Overview
`VCFX_alignment_checker` identifies discrepancies between VCF variant entries and a reference genome FASTA file. This tool helps validate that the reference alleles in a VCF file match the corresponding positions in the reference genome.

## Usage
```bash
VCFX_alignment_checker --alignment-discrepancy <vcf_file> <reference.fasta> > discrepancies.txt
```

## Options
| Option | Description |
|--------|-------------|
| `-a`, `--alignment-discrepancy` | Enable alignment discrepancy checking mode |
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |

## Description
`VCFX_alignment_checker` compares VCF variants against a reference genome to validate sequence consistency. The tool:

1. Loads a reference genome from a FASTA file into memory
2. Parses each variant in the input VCF file
3. For each variant:
   - Retrieves the corresponding sequence from the reference genome
   - Compares it with the REF and ALT values from the VCF
   - Reports any discrepancies found
4. Outputs a tab-separated report of all detected discrepancies

This tool is particularly useful for:
- Validating VCF files against their reference genome
- Identifying potential errors in variant calling
- Detecting misalignments in variant positions
- Quality control of VCF data before downstream analysis

## Input Requirements
- A VCF file with variant records
- A FASTA file containing the reference genome sequences
- VCF must contain standard CHROM, POS, REF, and ALT fields

## Output Format
The tool produces a tab-separated values (TSV) file with the following columns:

| Column | Description |
|--------|-------------|
| CHROM | Chromosome of the variant |
| POS | Position of the variant |
| ID | Variant identifier |
| REF | Reference allele in the VCF |
| ALT | Alternate allele in the VCF |
| Discrepancy_Type | Type of discrepancy detected (REF_DISCREPANCY or ALT_DISCREPANCY) |
| Reference_Value | The actual sequence from the reference genome |
| VCF_Value | The value from the VCF that differs from the reference |

## Examples

### Basic Usage
Check for discrepancies between variants in a VCF file and a reference genome:
```bash
VCFX_alignment_checker --alignment-discrepancy variants.vcf reference.fa > discrepancies.txt
```

### Pipeline Integration
Integrate with other tools for comprehensive validation:
```bash
VCFX_alignment_checker --alignment-discrepancy variants.vcf reference.fa | grep "REF_DISCREPANCY" > ref_errors.tsv
```

## Discrepancy Types

### REF_DISCREPANCY
Indicates that the REF field in the VCF doesn't match the reference genome at the specified position. This type of discrepancy suggests potential issues with:
- Incorrect variant calling
- Reference genome version mismatch
- Coordinate system errors

### ALT_DISCREPANCY
Indicates that the ALT field doesn't correspond to an expected variation from the reference. For SNPs, this can happen when:
- The variant quality is low
- The variant caller made an error
- There are assembly or alignment issues

## Handling Special Cases

### Chromosome Naming
The tool attempts to normalize chromosome names between the VCF and reference FASTA by:
- Adding 'chr' prefix when appropriate
- Checking for standard naming conventions (1-22, X, Y, MT)
- This normalization helps handle common mismatches between different naming conventions

### Indels and Complex Variants
For insertions, deletions, and complex variants:
- The tool compares the available bases (minimum length of REF and ALT)
- It checks if the REF allele matches the reference genome
- It also checks if the ALT allele differs from the reference as expected

### Missing Reference Sequences
If a chromosome in the VCF isn't found in the reference genome:
- A warning is issued to stderr
- The variant is skipped from discrepancy checking
- This helps identify potential naming mismatches between files

### Out-of-Range Positions
For positions beyond the reference sequence length:
- The tool issues a warning
- The variant is excluded from the discrepancy report
- This can help identify coordinate system issues

## Performance Considerations
- Loads the entire reference genome into memory for faster lookups
- Processes the VCF file sequentially, line by line
- Memory usage scales with the size of the reference genome
- Discrepancy checking is computationally efficient after loading the reference

## Limitations
- Requires loading the entire reference genome into memory, which can be memory-intensive for large genomes
- Chromosome name normalization may not handle all naming conventions
- Doesn't account for circularity in mitochondrial or bacterial genomes
- No specialized handling for structural variants
- Minimal sanity checks for FASTA format integrity
- Cannot check variants spanning multiple chromosomes 
