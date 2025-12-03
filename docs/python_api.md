# Python API

VCFX provides optional Python bindings exposing a subset of helper
functions from the C++ `vcfx_core` library. The bindings are built as a
native Python extension and can be enabled through CMake.

## Installation

Build the project with the `PYTHON_BINDINGS` option enabled:

```bash
mkdir build && cd build
cmake -DPYTHON_BINDINGS=ON ..
make -j
```

The compiled module will be placed in the `build/python` directory.
You can also install the package via `pip` which will invoke CMake
automatically. We recommend using a Python virtual environment because
some systems mark the system Python as "externally managed", which can
prevent direct installation. Create and activate a virtual environment
before running `pip`:

```bash
python3 -m venv venv && source venv/bin/activate
pip install --no-build-isolation ./python
```

If you are offline, ensure that the required build dependencies are already
installed because `--no-build-isolation` uses your current environment.

## Available Functions

The module exposes the following helpers:

- `trim(text)` – remove leading and trailing whitespace.
- `split(text, delimiter)` – split `text` by the given delimiter and
  return a list of strings.
- `read_file_maybe_compressed(path)` – read a plain or gzip/BGZF
  compressed file and return its contents as a string.
- `read_maybe_compressed(data)` – decompress a bytes object if it is
  gzip/BGZF compressed and return the resulting bytes.
- `get_version()` – return the VCFX version string.

## Example Usage

```python
import vcfx

print(vcfx.trim("  abc  "))
# 'abc'

print(vcfx.split("A,B,C", ","))
# ['A', 'B', 'C']

data = vcfx.read_maybe_compressed(b"hello")
print(data)

version = vcfx.get_version()
print("VCFX version:", version)
```

## Tool Wrappers

Besides the helper functions, the package provides lightweight wrappers for
all command line tools shipped with VCFX. The wrappers simply invoke the
corresponding ``VCFX_*`` executable via ``subprocess``.

For the wrappers to work, either the ``vcfx`` wrapper script or the individual
``VCFX_*`` binaries must be available on your ``PATH``. After building the
project you can source ``add_vcfx_tools_to_path.sh`` to add the build
directories to ``PATH``:

```bash
source /path/to/VCFX/add_vcfx_tools_to_path.sh
```

Use ``vcfx.available_tools()`` to see which tools are accessible on your
``PATH`` and call them either via ``vcfx.run_tool(name, *args)`` or by using
the tool name as a function:

```python
import vcfx

print(vcfx.available_tools())

# run through the generic helper
vcfx.run_tool("alignment_checker", "--help")
```

If ``VCFX_alignment_checker`` (or any other tool) is not present on
``PATH``, ``run_tool`` automatically tries to execute ``vcfx`` with the
tool name as the first argument when the ``vcfx`` wrapper is available.

For a full script demonstrating how to set up ``PATH`` with
``add_vcfx_tools_to_path.sh`` and handle errors when running tools,
see [``examples/python_usage.py``](../examples/python_usage.py).

## Convenience Wrappers

Several tools have Python helpers that run the command line program and
parse its output into structured data. Many wrappers return dataclasses
from `vcfx.results` rather than raw dictionaries. When using these
helpers numeric columns are automatically converted to ``int`` or ``float``
based on the dataclass field annotations.

```python
import vcfx

# Check alignment discrepancies and get dataclass instances
rows = vcfx.alignment_checker("tests/data/align_Y.vcf", "tests/data/align_refY.fa")
print(rows[0].Discrepancy_Type)  # 'ALT_MISMATCH'

# Count alleles for all samples
counts = vcfx.allele_counter("tests/data/allele_counter_A.vcf")
print(counts[0].Alt_Count)  # 1

# Simply count the variants in a VCF
n = vcfx.variant_counter("tests/data/variant_counter_normal.vcf")
print(n)
```

Additional wrappers are available for other tools:

```python
# Allele frequency calculation
freqs = vcfx.allele_freq_calc("tests/data/allele_freq_calc/simple.vcf")
print(freqs[0].Allele_Frequency)  # 0.5000

# Aggregate INFO fields and read the updated VCF text
annotated = vcfx.info_aggregator("tests/data/aggregator/basic.vcf", ["DP"])
print("#AGGREGATION_SUMMARY" in annotated)

# Parse INFO fields into dictionaries
info_rows = vcfx.info_parser("tests/data/info_parser/basic.vcf", ["DP"])
print(info_rows[0]["DP"])  # '10'

# Summarize INFO fields
summary = vcfx.info_summarizer("tests/data/info_summarizer/basic.vcf", ["DP"])
print(summary[0].Mean)  # 20.0000

# Convert VCF to FASTA alignment
fasta = vcfx.fasta_converter("tests/data/fasta_converter/basic.vcf")
print(fasta.splitlines()[0])  # '>SAMPLE1'
```

### Additional Wrappers

```python
# Calculate allele balance for all samples
balance = vcfx.allele_balance_calc("tests/data/allele_balance_calc_A.vcf")
print(balance[0].Allele_Balance)  # 1.000000

# Check concordance between two samples
conc = vcfx.concordance_checker(
    "tests/data/concordance_input.vcf",
    "SAMPLE1",
    "SAMPLE2",
)
print(conc[0].Concordance)  # 'Concordant'
```

```python
# Query variants with heterozygous genotypes
filtered = vcfx.genotype_query(
    "tests/data/genotype_query/sample.vcf",
    "0/1",
)
print(filtered.startswith("##"))  # True

# Remove duplicate records
dedup = vcfx.duplicate_remover("tests/data/allele_balance_calc_A.vcf")
print(dedup.splitlines()[0].startswith("#"))  # True
```

```python
# Subset variants by allele frequency range
subset = vcfx.af_subsetter("tests/data/af_subsetter_A.vcf", "0.01-0.1")
print(subset.startswith("##"))

# Detect missing genotypes
flagged = vcfx.missing_detector("tests/data/concordance_missing_data.vcf")
print("MISSING_GENOTYPES=1" in flagged)

# Hardy-Weinberg equilibrium test
hwe = vcfx.hwe_tester("tests/data/hwe_tester/basic_hwe.vcf")
print(hwe[0].HWE_pvalue)

# Inbreeding coefficient calculation
coeff = vcfx.inbreeding_calculator(
    "tests/data/inbreeding_calculator/single_sample_excludeSample_false.vcf",
    freq_mode="excludeSample",
)
print(coeff[0].InbreedingCoefficient)

# Variant classification
classes = vcfx.variant_classifier("tests/data/classifier_mixed.vcf")
print(classes[0].Classification)  # 'SNP'

# Cross-sample concordance
xconc = vcfx.cross_sample_concordance("tests/data/concordance_some_mismatch.vcf")
print(xconc[0].Concordance_Status)  # 'CONCORDANT'

# Extract specific fields
fields = vcfx.field_extractor(
    "tests/data/field_extractor_input.vcf", ["CHROM", "POS", "ID"]
)
print(fields[0]["ID"])  # 'rs123'
```

```python
# Assign ancestry to samples
assign = vcfx.ancestry_assigner(
    "tests/data/ancestry_assigner/input.vcf",
    "tests/data/ancestry_assigner/freq.tsv",
)
# ``assign`` is a list of :class:`vcfx.results.AncestryAssignment` objects
print(assign[0].Assigned_Population)  # 'EUR'
```
```python

# Calculate genotype dosages
dosage = vcfx.dosage_calculator("tests/data/dosage_calculator/basic.vcf")
print(dosage[0].Dosages)  # '0,1,2'
```

```python
# Create a position index
index_rows = vcfx.indexer("tests/data/indexer/basic.vcf")
print(index_rows[0].FILE_OFFSET)  # 255
```

### Further Wrappers

Additional helper functions mirror the rest of the ``VCFX_*`` tools. They work
like the wrappers above. A few examples:

```python
# Infer population ancestry
infer = vcfx.ancestry_inferrer(
    "tests/data/ancestry_inferrer/eur_samples.vcf",
    "tests/data/ancestry_inferrer/population_freqs.txt",
)
print(infer[0].Inferred_Population)  # 'EUR'

# Extract annotation fields
rows = vcfx.annotation_extractor(
    "tests/data/haplotype_extractor/basic.vcf",
    ["ANN", "Gene"],
)
print(rows[0]["Gene"])

# Normalize indels
normalized = vcfx.indel_normalizer("tests/data/indel_normalizer/basic.vcf")
print(normalized.startswith("##"))

# Validate a file
report = vcfx.validator("tests/data/variant_counter_normal.vcf")
print("VCF" in report)
```

## Installation from PyPI

The easiest way to install VCFX Python bindings is via PyPI:

```bash
pip install vcfx
```

This installs:
- Python bindings to the C++ helper functions
- Wrapper functions for all VCFX tools
- Type definitions and dataclasses for structured output

Note: The command-line tools themselves need to be installed separately (via Bioconda, Docker, or building from source) for the tool wrappers to work.

## Performance Considerations

The Python wrappers execute command-line tools via subprocess, which means:

1. **Overhead**: There's a small overhead for process creation
2. **Streaming**: Large files are processed in a streaming fashion by the underlying tools
3. **Memory**: Memory usage is determined by the C++ tools, not Python

For maximum performance with large datasets:
```python
# Process multiple files in parallel
from concurrent.futures import ProcessPoolExecutor
import vcfx

files = ["sample1.vcf", "sample2.vcf", "sample3.vcf"]

with ProcessPoolExecutor() as executor:
    results = list(executor.map(vcfx.variant_counter, files))
```

## Best Practices

### 1. Check tool availability
```python
import vcfx

# List available tools
tools = vcfx.available_tools()
if "variant_counter" not in tools:
    print("VCFX tools not found in PATH!")
```

### 2. Handle errors gracefully
```python
import subprocess

try:
    result = vcfx.hwe_tester("input.vcf")
except subprocess.CalledProcessError as e:
    print(f"Error: {e.stderr}")
except FileNotFoundError:
    print("VCFX tools not installed")
```

### 3. Use structured outputs
```python
# Prefer structured outputs over parsing text
frequencies = vcfx.allele_freq_calc("input.vcf")
# Good: Access typed fields
mean_af = sum(f.Allele_Frequency for f in frequencies) / len(frequencies)

# Avoid: Manual text parsing
output = vcfx.run_tool("allele_freq_calc", "input.vcf")
# Don't parse stdout manually when wrappers are available
```

### 4. Chain operations efficiently
```python
# Use pipes at the shell level for efficiency
import subprocess

# Efficient: Single pipeline
cmd = "VCFX_variant_classifier --append-info input.vcf | VCFX_phred_filter --min-qual 30"
result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

# Less efficient: Multiple Python calls
classified = vcfx.variant_classifier("input.vcf", append_info=True)
# Writing intermediate files is slower
```

## Troubleshooting

### Tools not found
If you get "command not found" errors:
1. Ensure VCFX tools are installed (via Bioconda or built from source)
2. Add tools to PATH: `export PATH=$PATH:/path/to/vcfx/build/src`
3. Or use the vcfx wrapper: `source /path/to/VCFX/add_vcfx_tools_to_path.sh`

### Import errors
If `import vcfx` fails:
1. Ensure you've installed the package: `pip install vcfx`
2. Check Python version: requires Python 3.10+
3. For development: `pip install -e /path/to/VCFX/python`

### Performance issues
For large files:
1. Ensure you have sufficient memory
2. Use streaming tools that process line-by-line
3. Consider splitting large files with `VCFX_file_splitter`
4. Process chunks in parallel when possible
