# Changelog

All notable changes to this project will be documented in this file.

## [1.1.4] - 2025-12-18
### Added
- GATK-compatible validation checks for VCFX_validator:
  - ALT allele observation check (ALLELES validation)
  - Empty VCF detection with `--allow-empty` override
  - Header Type/Number field validation
  - Variant sorting check (disable with `-S`)
  - AN/AC consistency check in strict mode (CHR_COUNTS)
  - REF validation against FASTA reference (`-R`)
  - dbSNP ID validation (`-D`)
  - GVCF format validation (`-g`)
- New `-i/--input` flag for validator mmap file input
- Test cases for all new validation features

### Changed
- Updated VCFX_validator documentation with comprehensive feature list
- Improved validation error messages with line numbers and context

### Performance
- Maintained ~110 MB/s throughput with all new validations enabled

## [1.0.3] - 2025-06-15
### Added
- PyPI package publishing via GitHub Actions
- Comprehensive PyPI badges (version, downloads, Python versions)
- Expanded Python README with detailed examples and workflows
- Python API best practices and troubleshooting guide
- Performance considerations documentation for Python users
- Python examples in Quick Start guide
- Structured data types documentation for all tool wrappers

### Changed
- Moved PyPI installation to primary installation method in docs
- Enhanced main documentation with Python API prominence
- Reorganized installation guide to highlight pip install
- Updated all documentation to feature Python bindings

### Fixed
- Version detection in pyproject.toml for dynamic versioning
- Linting issues (import order, trailing whitespace)
- GitHub Actions workflow for PyPI publishing
- Python type hints compatibility with older Python versions

## [1.0.2] - 2025-05-30
### Added
- Pre-commit configuration and improved linting
- Extended Python tests and wrappers
### Changed
- Updated documentation and CI scripts
### Fixed
- Compatibility issues with Python 3.10

## [1.0.1] - 2025-03-29
### Added
- Python API and wrapper scripts
- Build and packaging improvements
### Fixed
- Compilation issues on some platforms

## [1.0.0] - 2024-10-19
### Added
- Initial release with core command-line tools

