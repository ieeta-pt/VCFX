# VCFX: Comprehensive VCF Manipulation Toolkit

<p align="center">
  <img src="assets/images/VCFX.png" alt="VCFX Logo" width="180"/>
</p>

[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://ieeta-pt.github.io/VCFX/)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/vcfx.svg)](https://anaconda.org/bioconda/vcfx)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Docker](https://img.shields.io/badge/Docker-GHCR-blue)](https://ieeta-pt.github.io/VCFX/docker/)
[![PyPI version](https://img.shields.io/pypi/v/vcfx.svg)](https://pypi.org/project/vcfx/)
[![PyPI downloads](https://img.shields.io/pypi/dm/vcfx.svg)](https://pypi.org/project/vcfx/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/vcfx.svg)](https://pypi.org/project/vcfx/)

VCFX is a set of small C/C++ command line tools for manipulating and analysing Variant Call Format (VCF) files. Each tool does one job well and they can be chained together using standard streams.

## Features

- **60+ Specialized Tools** for filtering, transforming and analysing variants
- **Pipeline Ready**: tools read from stdin and write to stdout
- **Fast**: designed for large genomic datasets
- **Cross Platform**: Linux and macOS support
- **WebAssembly Builds** for browser or Node.js usage
- **Easy Installation** via PyPI, Bioconda or Docker
- **Python Bindings** for programmatic access

## Installation

### PyPI (Python Package)
```bash
pip install vcfx
```
After installing the Python bindings you can run any tool directly:
```python
import vcfx
vcfx.run_tool("alignment_checker", "--help")
```

### Bioconda
```bash
conda install -c bioconda vcfx
```

### Docker
```bash
docker pull ghcr.io/jorgemfs/vcfx:latest
docker run --rm ghcr.io/jorgemfs/vcfx:latest VCFX_tool_name --help
```

### Build from Source
```bash
git clone https://github.com/ieeta-pt/VCFX.git
cd VCFX
mkdir build && cd build
cmake .. -DPYTHON_BINDINGS=ON
make
```
Optionally run `make install` to place the tools in `~/.local/bin`.

## Quick Example
```bash
cat input.vcf | \
  VCFX_variant_classifier --append-info | \
  grep 'VCF_CLASS=SNP' | \
  VCFX_allele_freq_calc > snp_frequencies.tsv
```

## Documentation

Full documentation is available at [ieeta-pt.github.io/VCFX](https://ieeta-pt.github.io/VCFX/). The `docs` folder contains the sources and can be served locally with `mkdocs serve`.

## Repository Layout

- `src/` – C++ source code for all tools
- `python/` – optional Python bindings
- `docs/` – documentation sources
- `tests/` – shell and Python tests
- `examples/` – usage examples and workflows

## Development

Run tests with `pytest` and `ctest` from the build directory. Code style is enforced with `clang-format` and pre-commit hooks:
```bash
pre-commit install
```

## Benchmarks

Install the benchmark dependencies with conda as described in section 5 of this README. Then run:

```bash
cd benchmarks
make datasets
make tools-check
make all
```

The runtime summary will be written to `benchmarks/results.csv` and visualization notebooks are under `benchmarks/notebooks/`.
## Citation

If you use VCFX in your research please cite:
```
@inproceedings{silva2025vcfx,
  title={VCFX: A Minimalist, Modular Toolkit for Streamlined Variant Analysis},
  author={Silva, Jorge Miguel and Oliveira, José Luis},
  booktitle={12th International Work-Conference on Bioinformatics and Biomedical Engineering (IWBBIO 2025)},
  year={2025},
  organization={Springer}
}
```

## License

VCFX is distributed under the [MIT License](LICENSE).
