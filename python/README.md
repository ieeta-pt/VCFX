# VCFX Python Package

Python bindings for the VCFX toolkit - a comprehensive VCF manipulation toolkit with 60+ specialized tools for genomic variant analysis.

## Installation

```bash
pip install vcfx
```

For development or building from source:
```bash
git clone https://github.com/ieeta-pt/VCFX.git
cd VCFX
pip install -e ./python
```

## Quick Start

```python
import vcfx

# Basic helper functions
text = vcfx.trim("  hello  ")  # Returns "hello"
parts = vcfx.split("A,B,C", ",")  # Returns ["A", "B", "C"]

# Read compressed files
content = vcfx.read_file_maybe_compressed("variants.vcf.gz")

# Get version
version = vcfx.get_version()
print(f"VCFX version: {version}")
```

## Tool Wrappers

VCFX provides Python wrappers for all command-line tools. These wrappers execute the tools and parse their output into structured Python objects:

```python
# Count variants
count = vcfx.variant_counter("input.vcf")
print(f"Total variants: {count}")

# Calculate allele frequencies
freqs = vcfx.allele_freq_calc("input.vcf")
for freq in freqs:
    print(f"Position {freq.Pos}: AF={freq.Allele_Frequency}")

# Check sample concordance
concordance = vcfx.concordance_checker("input.vcf", "SAMPLE1", "SAMPLE2")
for row in concordance:
    print(f"{row.Position}: {row.Concordance}")
```

## Structured Data Types

Many tool wrappers return dataclass objects with typed fields for easy access:

```python
from vcfx.results import AlleleFrequency, VariantClassification

# Allele frequency calculations return AlleleFrequency objects
freqs = vcfx.allele_freq_calc("variants.vcf")
# Access fields with type safety and IDE completion
print(freqs[0].Chromosome)  # str
print(freqs[0].Allele_Frequency)  # float

# Variant classification returns VariantClassification objects
classes = vcfx.variant_classifier("variants.vcf")
print(classes[0].Classification)  # 'SNP', 'INDEL', 'MNV', or 'STRUCTURAL'
```

## Common Workflows

### Quality Control Pipeline
```python
import vcfx

# Validate VCF format
validation_report = vcfx.validator("input.vcf")
if "ERROR" in validation_report:
    print("VCF validation failed!")
    
# Detect missing data
missing_flagged = vcfx.missing_detector("input.vcf")

# Check concordance across samples
cross_concordance = vcfx.cross_sample_concordance("input.vcf")
discordant = [r for r in cross_concordance if r.Concordance_Status != 'CONCORDANT']
```

### Filtering and Analysis
```python
# Filter by allele frequency
af_filtered = vcfx.af_subsetter("input.vcf", "0.01-0.1")

# Extract specific samples
sample_vcf = vcfx.sample_extractor("input.vcf", ["SAMPLE1", "SAMPLE2"])

# Calculate Hardy-Weinberg equilibrium
hwe_results = vcfx.hwe_tester("input.vcf")
significant = [r for r in hwe_results if r.HWE_pvalue < 0.05]
```

### Population Genetics
```python
# Infer ancestry
ancestry = vcfx.ancestry_inferrer("samples.vcf", "population_freqs.txt")
for sample in ancestry:
    print(f"{sample.Sample}: {sample.Inferred_Population}")

# Calculate inbreeding coefficients
inbreeding = vcfx.inbreeding_calculator("input.vcf", freq_mode="excludeSample")
```

## Error Handling

Tool wrappers raise `subprocess.CalledProcessError` if the underlying tool fails:

```python
try:
    result = vcfx.variant_counter("nonexistent.vcf")
except subprocess.CalledProcessError as e:
    print(f"Tool failed with exit code {e.returncode}")
    print(f"Error output: {e.stderr}")
```

## Requirements

- Python 3.10+
- For tool wrappers: VCFX command-line tools must be installed and available in PATH
  - Install via conda: `conda install -c bioconda vcfx`
  - Or build from source (see documentation)

## Available Tools

VCFX includes 60+ specialized tools organized into categories:

- **Analysis**: allele frequencies, variant classification, HWE testing, LD calculation
- **Filtering**: quality filters, population filters, missing data filters
- **Transformation**: sample extraction, multiallelic splitting, normalization
- **Quality Control**: concordance checking, validation, outlier detection
- **File Management**: indexing, compression, merging, splitting
- **Annotation**: custom annotation, INFO field processing

Use `vcfx.available_tools()` to list tools accessible in your environment.

## Documentation

- Full documentation: https://ieeta-pt.github.io/VCFX/
- Python API reference: https://ieeta-pt.github.io/VCFX/python_api/
- Tool documentation: https://ieeta-pt.github.io/VCFX/tools_overview/

## License

MIT License - see LICENSE file for details. 
