# VCFX: A Comprehensive VCF Manipulation Toolkit

VCFX is a collection of specialized command-line tools designed for efficient manipulation, analysis, and transformation of VCF (Variant Call Format) files used in genomic research and bioinformatics. Available via PyPI, Bioconda, and Docker.

## What is VCFX?

VCFX follows the Unix philosophy of creating small, focused tools that do one thing well and can be combined together to form powerful workflows. Each tool in the VCFX suite is optimized for a specific VCF-related task, with optional Python bindings for programmatic access.

Key features:
- 60+ specialized command-line tools
- Python API with structured data types
- Easy installation via `pip install vcfx`
- Cross-platform support (Linux, macOS)
- Composable tools for pipeline integration

VCFX enables researchers and bioinformaticians to:

- Extract specific information from VCF files
- Filter variants based on various criteria
- Transform VCF data into different formats
- Analyze genotypes and compute statistics
- Validate and check VCF file integrity
- Manipulate structural variants and complex records

## Getting Started

To begin using VCFX, first follow the [installation instructions](#installation) below, then explore the [tool categories](#tool-categories) to find the right components for your workflow.

### Installation

Choose your preferred installation method:

#### PyPI (Python users)
```bash
pip install vcfx
```

#### Bioconda (includes all tools)
```bash
conda install -c bioconda vcfx
```

#### Build from source
```bash
git clone https://github.com/ieeta-pt/VCFX.git
cd VCFX
mkdir -p build
cd build
cmake -DPYTHON_BINDINGS=ON ..
make
```

See the [full installation guide](installation.md) for more options including Docker.

### Basic Example

Here's a simple example of using VCFX to analyze variants:

```bash
# Command line usage
cat input.vcf | \
  VCFX_variant_classifier --append-info | \
  grep 'VCF_CLASS=SNP' | \
  VCFX_allele_freq_calc > snp_frequencies.tsv
```

```python
# Python API usage
import vcfx

# Count variants
count = vcfx.variant_counter("input.vcf")

# Calculate allele frequencies with structured output
freqs = vcfx.allele_freq_calc("input.vcf")
for f in freqs:
    print(f"Position {f.Pos}: AF={f.Allele_Frequency}")
```

## Tool Categories

The VCFX toolkit includes tools in the following categories:

### Data Analysis

Tools for extracting statistical information and insights from variant data:

- [VCFX_allele_freq_calc](analysis/VCFX_allele_freq_calc.md) - Calculate allele frequencies
- [VCFX_variant_classifier](analysis/VCFX_variant_classifier.md) - Classify variants into SNP, INDEL, MNV, or STRUCTURAL
- [VCFX_inbreeding_calculator](analysis/VCFX_inbreeding_calculator.md) - Calculate inbreeding coefficients
- [VCFX_dosage_calculator](analysis/VCFX_dosage_calculator.md) - Calculate allele dosage from genotypes
- [View all analysis tools...](tools_overview.md#data-analysis)

### Data Filtering

Tools for selecting variants based on specific criteria:

- [VCFX_phase_checker](filtering/VCFX_phase_checker.md) - Filter variants to keep only fully phased genotypes
- [VCFX_phred_filter](filtering/VCFX_phred_filter.md) - Filter variants based on Phred-scaled quality scores
- [VCFX_record_filter](filtering/VCFX_record_filter.md) - Filter variants based on various VCF fields
- [View all filtering tools...](tools_overview.md#data-filtering)

### Data Transformation

Tools for converting or reformatting VCF data:

- [VCFX_multiallelic_splitter](transformation/VCFX_multiallelic_splitter.md) - Split multiallelic variants into biallelic records
- [VCFX_sample_extractor](transformation/VCFX_sample_extractor.md) - Extract specific samples from a VCF file
- [VCFX_indel_normalizer](transformation/VCFX_indel_normalizer.md) - Normalize indel representations
- [View all transformation tools...](tools_overview.md#data-transformation)

### Quality Control

Tools for validating and checking data quality:

- [VCFX_concordance_checker](quality-control/VCFX_concordance_checker.md) - Check concordance between samples in a VCF file
- [VCFX_missing_detector](quality-control/VCFX_missing_detector.md) - Detect and report missing data
- [VCFX_validator](quality-control/VCFX_validator.md) - Validate VCF format compliance
- [View all quality control tools...](tools_overview.md#quality-control)

### File Management

Tools for handling VCF files:

- [VCFX_indexer](file-management/VCFX_indexer.md) - Create an index file for random access
- [VCFX_file_splitter](file-management/VCFX_file_splitter.md) - Split VCF files into smaller chunks
- [VCFX_compressor](file-management/VCFX_compressor.md) - Compress VCF files efficiently
- [View all file management tools...](tools_overview.md#file-management)

### Annotation and Reporting

Tools for annotating and extracting information from VCF files:

- [VCFX_custom_annotator](annotation/VCFX_custom_annotator.md) - Add custom annotations to VCF files
- [VCFX_info_summarizer](annotation/VCFX_info_summarizer.md) - Summarize INFO fields
- [VCFX_field_extractor](annotation/VCFX_field_extractor.md) - Extract specific fields from VCF files
- [View all annotation tools...](tools_overview.md#annotation-and-reporting)

### Data Processing

Tools for processing variants and samples:

- [VCFX_missing_data_handler](processing/VCFX_missing_data_handler.md) - Handle missing data
- [VCFX_quality_adjuster](processing/VCFX_quality_adjuster.md) - Adjust quality scores
- [VCFX_haplotype_phaser](processing/VCFX_haplotype_phaser.md) - Phase haplotypes
- [VCFX_haplotype_extractor](processing/VCFX_haplotype_extractor.md) - Extract haplotype information
- [View all processing tools...](tools_overview.md#data-processing)

For a complete list of all tools and detailed usage examples, see the [tools overview](tools_overview.md).

## Python API

VCFX provides comprehensive Python bindings that wrap all command-line tools and provide additional conveniences:

- **Structured data types**: Tool outputs are parsed into dataclasses with proper typing
- **Easy integration**: All tools accessible via `vcfx.tool_name()` functions
- **Error handling**: Clear exceptions when tools fail
- **Helper functions**: Utilities for reading compressed files, text processing, etc.

Example:
```python
import vcfx

# Tools return structured data, not just strings
results = vcfx.hwe_tester("variants.vcf")
for variant in results:
    if variant.HWE_pvalue < 0.05:
        print(f"Variant at {variant.Pos} deviates from HWE (p={variant.HWE_pvalue})")
```

See the [Python API documentation](python_api.md) for complete details.

## Who Should Use VCFX?

VCFX is designed for:

- **Bioinformaticians** working with genomic variant data
- **Researchers** analyzing VCF files from sequencing projects
- **Pipeline developers** creating reproducible genomic workflows
- **Data scientists** extracting information from genetic variants

## Key Features

- **Composability**: All tools work with standard input/output for easy pipeline integration
- **Efficiency**: Optimized for performance with large genomic datasets
- **Robustness**: Careful error handling and validation of VCF formatting
- **Flexibility**: Works with various VCF versions and extensions
- **Simplicity**: Clear, focused tools with consistent interfaces

## Common Usage Patterns

VCFX tools are designed to be used in pipelines. Here are some common usage patterns:

### Basic Filtering and Analysis

```bash
# Extract phased variants, filter by quality, and calculate allele frequencies
cat input.vcf | \
  VCFX_phase_checker | \
  VCFX_phred_filter --phred-filter 30 | \
  VCFX_allele_freq_calc > result.tsv
```

### Sample Comparison

```bash
# Check concordance between two samples in a single VCF
cat input.vcf | VCFX_concordance_checker --samples "SAMPLE1 SAMPLE2" > concordance_report.tsv
```

See the [tools overview page](tools_overview.md#common-usage-patterns) for more usage examples.

## Community and Support

- [GitHub Repository](https://github.com/ieeta-pt/VCFX)
- [Issue Tracker](https://github.com/ieeta-pt/VCFX/issues)
- [Contributing Guidelines](CONTRIBUTING.md)

## Citation

If you use VCFX in your research, please cite:

```bibtex
@inproceedings{silva2025vcfx,
  title={VCFX: A Minimalist, Modular Toolkit for Streamlined Variant Analysis},
  author={Silva, Jorge Miguel and Oliveira, Jos{\'e} Luis},
  booktitle={12th International Work-Conference on Bioinformatics and Biomedical Engineering (IWBBIO 2025)},
  year={2025},
  organization={Springer}
}
```

## License

VCFX is available under MIT License. See the [LICENSE](LICENSE.md) file for details. 
