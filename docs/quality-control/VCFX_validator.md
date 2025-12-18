# VCFX_validator

## Overview

`VCFX_validator` is a high-performance, GATK-compatible VCF validation tool. It performs comprehensive checks on VCF file structure, header format, and data lines to ensure compliance with the VCF 4.x specification. The tool supports all major validation modes found in GATK's `ValidateVariants`, including reference validation, dbSNP ID checking, and GVCF format validation.

## Usage

```bash
VCFX_validator [OPTIONS] [input.vcf]
VCFX_validator [OPTIONS] -i input.vcf
VCFX_validator [OPTIONS] < input.vcf
```

## Options

### Basic Options

| Option | Description |
|--------|-------------|
| `-h`, `--help` | Display help message and exit |
| `-i`, `--input FILE` | Input VCF file (uses memory-mapped I/O) |
| `-s`, `--strict` | Enable strict mode (warnings become errors) |
| `-d`, `--report-dups` | Report duplicate records to stderr |
| `-n`, `--no-dup-check` | Skip duplicate detection (faster) |
| `-e`, `--allow-empty` | Allow VCF files with no variant records |
| `-b`, `--bloom-size N` | Bloom filter size in MB (default: 128) |
| `-t`, `--threads N` | Reserved for future multi-threaded validation |

### GATK-Compatible Validation Options

| Option | GATK Equivalent | Description |
|--------|-----------------|-------------|
| `-R`, `--reference FILE` | `--reference` | Validate REF alleles against FASTA reference |
| `-D`, `--dbsnp FILE` | `--dbsnp` | Validate variant IDs against dbSNP VCF |
| `-g`, `--gvcf` | `--gvcf` | Enable GVCF-specific validation |
| `-S`, `--no-sorting-check` | N/A | Skip variant sorting validation |
| `-C`, `--no-chr-counts` | N/A | Skip AN/AC consistency validation |

## Validation Checks

### Default Validations (Always On)

| Check | Description |
|-------|-------------|
| **VCF Structure** | Meta-info lines (##), header line (#CHROM), column count |
| **Header Definitions** | INFO/FORMAT fields have valid Type, Number, Description |
| **Position Values** | POS is a positive integer |
| **REF/ALT Sequences** | Contains only valid bases (A, C, G, T, N) |
| **QUAL Values** | Either `.` or a non-negative number |
| **FILTER Field** | Not empty |
| **Genotype Format** | Valid GT syntax (0/1, 1\|0, etc.) |
| **ALT Allele Observation** | ALT alleles are observed in sample genotypes (GATK ALLELES check) |
| **Variant Sorting** | Records are sorted by CHROM and POS (disable with `-S`) |

### Strict Mode Validations (`--strict`)

| Check | Description |
|-------|-------------|
| **Column Count** | Data lines must match header column count exactly |
| **FORMAT Field Count** | Sample fields must have correct sub-field count |
| **AN/AC Consistency** | Allele counts (AC) sum must not exceed total alleles (AN) |
| **Duplicate Detection** | Reports duplicate variants |
| **Warnings as Errors** | All warnings become errors |

### Optional Validations

| Check | Option | Description |
|-------|--------|-------------|
| **REF Validation** | `-R <fasta>` | REF allele matches reference FASTA at position |
| **ID Validation** | `-D <dbsnp.vcf>` | Variant IDs exist in dbSNP |
| **GVCF Validation** | `-g` | `<NON_REF>` allele, END tag, coverage continuity |

## Examples

### Basic Validation

```bash
# Validate from stdin
VCFX_validator < input.vcf

# Validate from file (uses mmap for better performance)
VCFX_validator input.vcf

# Validate compressed file
zcat input.vcf.gz | VCFX_validator
```

### Strict Validation

```bash
# Enable all strict checks
VCFX_validator --strict input.vcf

# Strict mode with duplicate reporting
VCFX_validator --strict --report-dups input.vcf
```

### GATK-Compatible Validation

```bash
# Validate REF against reference genome
VCFX_validator --reference hg38.fa input.vcf

# Validate IDs against dbSNP
VCFX_validator --dbsnp dbsnp156.vcf input.vcf

# Full GATK-style validation
VCFX_validator --strict --reference hg38.fa --dbsnp dbsnp156.vcf input.vcf

# GVCF validation
VCFX_validator --gvcf sample.g.vcf
```

### Pipeline Integration

```bash
# Validate before processing
VCFX_validator input.vcf && VCFX_allele_freq_calc < input.vcf > freqs.tsv

# Filter and validate
VCFX_phred_filter -p 30 < input.vcf | VCFX_validator

# Save validation errors
VCFX_validator input.vcf 2> validation_errors.txt
```

## Output

### Validation Report

When validation completes, a summary report is printed:

```
=== VCF Validation Report ===
Status: PASSED

File Statistics:
  Total lines:     1000004
  Header lines:    4
  Variant records: 1000000
  Samples:         4

Header Definitions:
  INFO fields:     1
  FORMAT fields:   1

Validation Checks Performed:
  [OK] VCF header structure
  [OK] Meta-information lines (##)
  [OK] Column header (#CHROM)
  [OK] Required columns
  [OK] Position values (POS > 0)
  [OK] REF/ALT allele sequences
  [OK] QUAL values
  [OK] INFO field definitions
  [OK] FORMAT field definitions
  [OK] Genotype values
  [OK] ALT alleles observed (GATK ALLELES)
  [OK] Variant sorting
  [OK] Duplicate detection
```

### Error Messages

```
# Missing Type in header
Error: header line missing Type at line 2.

# Invalid Number value
Error: invalid Number 'XXX' in header at line 2 (must be a non-negative integer, A, R, G, or .).

# ALT not observed (warning)
Warning: ALT allele 2 at position chr1:14804770 is not observed in any sample genotype on line 4.

# Unsorted variants (warning)
Warning: Variants not sorted: position 100 comes after 200 on chr1 at line 6.

# AN/AC mismatch (strict mode error)
Error: AC sum (99) exceeds AN (4) on line 6.

# Empty VCF
Error: VCF file contains no variant records (header-only file).
       Use --allow-empty to accept VCF files without variant data.

# REF mismatch (with -R)
Error: REF 'T' does not match reference 'A' at chr1:100 on line 5.
```

### Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Valid VCF (no errors) |
| 1 | Invalid VCF (errors found) |

## Performance

### Optimizations

- **Memory-mapped I/O**: File path arguments use mmap with `MADV_SEQUENTIAL` for optimal read-ahead
- **SIMD parsing**: AVX2/SSE2 vectorized newline detection on x86_64
- **Bloom filter**: Memory-efficient duplicate and dbSNP ID detection (default 128MB)
- **DNA lookup table**: 256-byte LUT for O(1) base validation
- **FORMAT caching**: Avoids re-parsing identical FORMAT strings
- **Zero-copy parsing**: `string_view` throughout for allocation-free parsing
- **Reference mmap**: FASTA reference is memory-mapped for fast REF validation

### Benchmarks

| Dataset | Variants | Samples | Size | Time | Throughput |
|---------|----------|---------|------|------|------------|
| 1M synthetic | 1,000,000 | 4 | 49 MB | 0.47s | 104 MB/s |
| chr21 1KG | 427,409 | 2,504 | 4.3 GB | ~40s | 110 MB/s |

### Performance Tips

```bash
# Use file path for mmap (fastest)
VCFX_validator large_file.vcf

# Skip duplicate check for known-clean data
VCFX_validator --no-dup-check large_file.vcf

# Reduce memory with smaller bloom filter
VCFX_validator --bloom-size 64 large_file.vcf

# Skip sorting check for pre-sorted files
VCFX_validator --no-sorting-check large_file.vcf
```

## Comparison with Other Tools

| Feature | VCFX_validator | GATK ValidateVariants | bcftools view |
|---------|----------------|----------------------|---------------|
| Basic VCF validation | ✓ | ✓ | ✓ |
| Header Type/Number check | ✓ | ✓ | Partial |
| ALT allele observation | ✓ | ✓ (ALLELES) | ✗ |
| REF vs FASTA | ✓ | ✓ (REF) | ✗ |
| dbSNP ID check | ✓ | ✓ (IDS) | ✗ |
| AN/AC consistency | ✓ | ✓ (CHR_COUNTS) | ✗ |
| Sorting check | ✓ | ✓ | ✗ |
| GVCF validation | ✓ | ✓ | ✗ |
| Empty file detection | ✓ | ✗ | ✗ |
| Duplicate detection | ✓ | ✗ | ✗ |
| Memory-mapped I/O | ✓ | ✗ | ✗ |
| Streaming stdin | ✓ | ✗ | ✓ |
| Throughput | ~110 MB/s | ~20 MB/s | ~80 MB/s |

## Limitations

- VCF version compatibility is not strictly enforced
- Threading (`-t`) is reserved for future implementation
- Bloom filter may produce rare false positives for duplicates (<1% at 100M variants)
- Reference validation requires indexed FASTA (simple sequential scan, not faidx)

## See Also

- [VCFX_concordance_checker](VCFX_concordance_checker.md) - Compare genotypes between VCF files
- [VCFX_missing_detector](VCFX_missing_detector.md) - Detect missing genotype data
- [VCFX_duplicate_remover](../file-management/VCFX_duplicate_remover.md) - Remove duplicate variants
