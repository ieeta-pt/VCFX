# VCFX_sv_handler

## Overview
`VCFX_sv_handler` is a utility tool for filtering and modifying structural variant (SV) records in VCF files. It can identify variants with the SVTYPE annotation, either keeping only structural variants or enhancing their annotations with additional information.

## Usage
```bash
VCFX_sv_handler [OPTIONS] < input.vcf > output.vcf
```

## Options
| Option | Description |
|--------|-------------|
| `-h`, `--help` | Display help message and exit (handled by `vcfx::handle_common_flags`) |
| `-v`, `--version` | Show program version and exit (handled by `vcfx::handle_common_flags`) |
| `-f`, `--sv-filter-only` | Keep only lines that have 'SVTYPE=' in their INFO field |
| `-m`, `--sv-modify` | Modify the INFO field of structural variants to add additional annotations |

## Description
`VCFX_sv_handler` processes a VCF file to manage structural variant records by:

1. Reading the VCF file from standard input
2. Preserving all header lines without modification
3. For each data line:
   - Checking if it contains 'SVTYPE=' in the INFO field to identify structural variants
   - If filtering is enabled, keeping only structural variant records
   - If modification is enabled, adding additional annotations to structural variant records
4. Writing the processed VCF to standard output

Structural variants (SVs) are genomic alterations that involve segments of DNA and include deletions (DEL), duplications (DUP), inversions (INV), and breakends (BND). This tool helps to specifically process these variants, which often require special handling in downstream analyses.

## Modification Details

When the `--sv-modify` option is used, the tool adds several annotations to structural variant records:

1. For all structural variants:
   - Adds `SV_VALIDATED=1` to indicate the variant has been processed

2. For deletions (DEL) and duplications (DUP):
   - Calculates and adds `SV_SIZE=<size>` based on the difference between END and POS positions

3. For inversions (INV):
   - Adds `INV_TYPE=PARALLEL`

4. For breakends (BND):
   - Adds `BND_ORIENTATION=PAIR`

These modifications can be useful for downstream analyses that rely on standardized annotations or require specific information about structural variants.

## Examples

### Filter Structural Variants
Keep only structural variant records in a VCF file:
```bash
VCFX_sv_handler --sv-filter-only < input.vcf > sv_only.vcf
```

### Modify Structural Variants
Enhance structural variant records with additional annotations:
```bash
VCFX_sv_handler --sv-modify < input.vcf > annotated.vcf
```

### Combined Operation
Filter and modify structural variants in one operation:
```bash
VCFX_sv_handler --sv-filter-only --sv-modify < input.vcf > processed_sv.vcf
```

## Example Transformations

### Filtering Structural Variants
```
Before:
chr1 100 . A T . PASS DP=30
chr1 200 . T G . PASS SVTYPE=DEL;END=300

After (with --sv-filter-only):
chr1 200 . T G . PASS SVTYPE=DEL;END=300
```

### Modifying Structural Variants
```
Before:
chr1 200 . T G . PASS SVTYPE=DEL;END=300

After (with --sv-modify):
chr1 200 . T G . PASS SVTYPE=DEL;END=300;SV_VALIDATED=1;SV_SIZE=100
```

```
Before:
chr1 400 . G C . PASS SVTYPE=INV;END=500

After (with --sv-modify):
chr1 400 . G C . PASS SVTYPE=INV;END=500;SV_VALIDATED=1;INV_TYPE=PARALLEL
```

## Special Case Handling

### Multiple Options
- When both `--sv-filter-only` and `--sv-modify` are used, records are first filtered and then modified

### Missing Fields
- If a variant has SVTYPE but no END position, it will still be processed
- For size calculations, if either POS or END is missing/invalid, SV_SIZE will not be added

### Malformed Lines
- Lines with fewer than 8 columns are skipped with a warning
- Lines with invalid POS values are skipped when attempting to modify

### Non-Structural Variants
- Non-SV records are preserved unless filtering is enabled

## Performance Considerations
- The tool processes the VCF file line by line, with minimal memory requirements
- Performance scales linearly with the size of the input file
- No external dependencies or reference files are required

## Limitations
- The tool identifies structural variants solely by the presence of "SVTYPE=" in the INFO field
- It adds fixed annotations without validating if they're appropriate for the specific variant
- No support for custom annotations or handling of specific structural variant subtypes
- Does not validate the correctness of existing structural variant annotations
- Cannot handle compressed (gzipped) VCF files directly 
