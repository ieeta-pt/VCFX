# VCFX Tools Overview

VCFX is a collection of C/C++ tools for processing and analyzing VCF (Variant Call Format) files, with optional WebAssembly compatibility. Each tool is an independent command-line executable that can parse input from `stdin` and write to `stdout`, enabling flexible piping and integration into bioinformatics pipelines.

The suite also includes a convenience wrapper `vcfx` so you can run commands as `vcfx <subcommand>`. For example, `vcfx variant_counter` is equivalent to running `VCFX_variant_counter`. Use `vcfx --list` or the alias `vcfx list` to see available subcommands. To view Markdown documentation for a tool, run `vcfx help <tool>`. All individual `VCFX_*` binaries remain available if you prefer calling them directly.
Every tool now supports the common flags `--help` and `--version` for quick usage or version information.

## Tool Categories

### Data Analysis

These tools help extract statistical information and insights from variant data:

- [VCFX_allele_freq_calc](analysis/VCFX_allele_freq_calc.md) - Calculate allele frequencies
- [VCFX_variant_classifier](analysis/VCFX_variant_classifier.md) - Classify variants into SNP, INDEL, MNV, or STRUCTURAL
- [VCFX_inbreeding_calculator](analysis/VCFX_inbreeding_calculator.md) - Calculate inbreeding coefficients
- [VCFX_dosage_calculator](analysis/VCFX_dosage_calculator.md) - Calculate allele dosage from genotypes
- [VCFX_hwe_tester](analysis/VCFX_hwe_tester.md) - Test for Hardy-Weinberg equilibrium
- [VCFX_distance_calculator](analysis/VCFX_distance_calculator.md) - Calculate genetic distances between samples
- [VCFX_allele_counter](analysis/VCFX_allele_counter.md) - Count alleles in VCF files
- [VCFX_allele_balance_calc](analysis/VCFX_allele_balance_calc.md) - Calculate allele balance metrics
- [VCFX_variant_counter](analysis/VCFX_variant_counter.md) - Count variants in VCF files
- [VCFX_ancestry_inferrer](analysis/VCFX_ancestry_inferrer.md) - Infer ancestry from genetic data
- [VCFX_ancestry_assigner](analysis/VCFX_ancestry_assigner.md) - Assign ancestry to samples
- [VCFX_ld_calculator](analysis/VCFX_ld_calculator.md) - Calculate pairwise linkage disequilibrium (r²) between variants

### Data Filtering

Tools for selecting variants based on specific criteria:

- [VCFX_phase_checker](filtering/VCFX_phase_checker.md) - Filter variants to keep only fully phased genotypes
- [VCFX_phred_filter](filtering/VCFX_phred_filter.md) - Filter variants based on Phred-scaled quality scores
- [VCFX_record_filter](filtering/VCFX_record_filter.md) - Filter variants based on various VCF fields
- [VCFX_gl_filter](filtering/VCFX_gl_filter.md) - Filter variants based on genotype likelihoods
- [VCFX_allele_balance_filter](filtering/VCFX_allele_balance_filter.md) - Filter variants based on allele balance
- [VCFX_population_filter](filtering/VCFX_population_filter.md) - Filter variants based on population statistics
- [VCFX_probability_filter](filtering/VCFX_probability_filter.md) - Filter variants based on probability scores
- [VCFX_nonref_filter](filtering/VCFX_nonref_filter.md) - Filter to keep only non-reference variants
- [VCFX_impact_filter](filtering/VCFX_impact_filter.md) - Filter variants based on predicted impact
- [VCFX_phase_quality_filter](filtering/VCFX_phase_quality_filter.md) - Filter variants based on phasing quality scores
- [VCFX_region_subsampler](filtering/VCFX_region_subsampler.md) - Filter variants based on genomic regions

### Data Transformation

Tools for converting or reformatting VCF data:

- [VCFX_multiallelic_splitter](transformation/VCFX_multiallelic_splitter.md) - Split multiallelic variants into biallelic records
- [VCFX_sample_extractor](transformation/VCFX_sample_extractor.md) - Extract specific samples from a VCF file
- [VCFX_position_subsetter](transformation/VCFX_position_subsetter.md) - Extract variants at specific positions
- [VCFX_format_converter](transformation/VCFX_format_converter.md) - Convert VCF files to other formats
- [VCFX_genotype_query](transformation/VCFX_genotype_query.md) - Query specific genotype patterns
- [VCFX_indel_normalizer](transformation/VCFX_indel_normalizer.md) - Normalize indel representations
- [VCFX_sv_handler](transformation/VCFX_sv_handler.md) - Handle structural variants in VCF files
- [VCFX_fasta_converter](transformation/VCFX_fasta_converter.md) - Convert VCF files to FASTA format
- [VCFX_sorter](transformation/VCFX_sorter.md) - Sort VCF files by position
- [VCFX_af_subsetter](transformation/VCFX_af_subsetter.md) - Extract variants based on allele frequency
- [VCFX_reformatter](transformation/VCFX_reformatter.md) - Reformat VCF files for better readability

### Quality Control

Tools for validating and checking data quality:

- [VCFX_concordance_checker](quality-control/VCFX_concordance_checker.md) - Check concordance between samples in a VCF file
- [VCFX_missing_detector](quality-control/VCFX_missing_detector.md) - Detect and report missing data
- [VCFX_outlier_detector](quality-control/VCFX_outlier_detector.md) - Detect outlier samples or variants
- [VCFX_alignment_checker](quality-control/VCFX_alignment_checker.md) - Check alignment of variants
- [VCFX_cross_sample_concordance](quality-control/VCFX_cross_sample_concordance.md) - Check concordance between samples
- [VCFX_validator](quality-control/VCFX_validator.md) - Validate VCF format compliance

### File Management

Tools for handling VCF files:

- [VCFX_indexer](file-management/VCFX_indexer.md) - Create an index file for random access
- [VCFX_file_splitter](file-management/VCFX_file_splitter.md) - Split VCF files into smaller chunks
- [VCFX_compressor](file-management/VCFX_compressor.md) - Compress VCF files efficiently
- [VCFX_diff_tool](file-management/VCFX_diff_tool.md) - Find differences between VCF files
- [VCFX_subsampler](file-management/VCFX_subsampler.md) - Subsample variants from a VCF file
- [VCFX_duplicate_remover](file-management/VCFX_duplicate_remover.md) - Remove duplicate variants
- [VCFX_merger](file-management/VCFX_merger.md) - Merge multiple VCF files by position

### Annotation and Reporting

Tools for annotating and extracting information from VCF files:

- [VCFX_custom_annotator](annotation/VCFX_custom_annotator.md) - Add custom annotations to VCF files
- [VCFX_info_summarizer](annotation/VCFX_info_summarizer.md) - Summarize INFO fields in VCF files
- [VCFX_header_parser](annotation/VCFX_header_parser.md) - Parse and extract information from VCF headers
- [VCFX_annotation_extractor](annotation/VCFX_annotation_extractor.md) - Extract annotations from VCF files
- [VCFX_ref_comparator](annotation/VCFX_ref_comparator.md) - Compare variants against a reference genome
- [VCFX_field_extractor](annotation/VCFX_field_extractor.md) - Extract specific fields from VCF files
- [VCFX_info_aggregator](annotation/VCFX_info_aggregator.md) - Aggregate INFO fields across variants
- [VCFX_info_parser](annotation/VCFX_info_parser.md) - Parse INFO fields in VCF files
- [VCFX_metadata_summarizer](annotation/VCFX_metadata_summarizer.md) - Summarize key metadata from VCF files

### Data Processing

Tools for processing variants and samples:

- [VCFX_missing_data_handler](processing/VCFX_missing_data_handler.md) - Handle missing data in VCF files
- [VCFX_quality_adjuster](processing/VCFX_quality_adjuster.md) - Adjust quality scores in VCF files
- [VCFX_haplotype_phaser](processing/VCFX_haplotype_phaser.md) - Phase haplotypes in VCF files
- [VCFX_haplotype_extractor](processing/VCFX_haplotype_extractor.md) - Extract haplotype information

## Common Usage Patterns

VCFX tools are designed to be combined in pipelines. Here are some common usage patterns:

### Basic Filtering and Analysis

```bash
# Extract phased variants, filter by quality, and calculate allele frequencies
cat input.vcf | \
  VCFX_phase_checker | \
  VCFX_phred_filter --phred-filter 30 | \
  VCFX_allele_freq_calc > result.tsv
```

### Variant Classification and Filtering

```bash
# Classify variants and filter for SNPs with high quality
cat input.vcf | \
  VCFX_variant_classifier --append-info | \
  grep 'VCF_CLASS=SNP' | \
  VCFX_phred_filter --phred-filter 30 > high_quality_snps.vcf
```

### Sample Comparison

```bash
# Check concordance between two samples in a single VCF
cat input.vcf | VCFX_concordance_checker --samples "SAMPLE1 SAMPLE2" > concordance.tsv
```

### Linkage Disequilibrium Analysis

```bash
# Calculate LD in a specific region after filtering for common variants
cat input.vcf | \
  VCFX_af_subsetter --af-filter '0.05-1.0' | \
  VCFX_ld_calculator --region chr1:10000-20000 > ld_matrix.txt
```

### Normalization and Splitting

```bash
# Normalize indels and split multiallelic variants
cat input.vcf | \
  VCFX_indel_normalizer | \
  VCFX_multiallelic_splitter > normalized_biallelic.vcf
```

### Population Analysis

```bash
# Extract population-specific VCFs and calculate allele frequencies
cat input.vcf | VCFX_population_filter --population EUR --pop-map pop_map.txt > eur.vcf
cat eur.vcf | VCFX_allele_freq_calc > eur_afs.tsv
```

### Quality Control Pipeline

```bash
# Validate, classify, detect missing data, and filter by quality
cat input.vcf | \
  VCFX_validator | \
  VCFX_variant_classifier --append-info | \
  VCFX_missing_detector --max-missing 0.1 | \
  VCFX_phred_filter --phred-filter 20 > qc_passed.vcf
``` 
