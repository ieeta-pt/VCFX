#!/usr/bin/env bash

set -e
set -o pipefail

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Get the root directory of the project
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"

echo "=== Testing VCFX_ld_calculator ==="

# Paths (work from either root or tests directory)
CALCULATOR_BIN="${ROOT_DIR}/build/src/VCFX_ld_calculator/VCFX_ld_calculator"
TEST_DATA_DIR="${SCRIPT_DIR}/data/ld_calculator"
EXPECTED_DIR="${SCRIPT_DIR}/expected/ld_calculator"
OUTPUT_DIR="${SCRIPT_DIR}/out/ld_calculator"

# Check if executable exists
if [ ! -f "$CALCULATOR_BIN" ]; then
  echo "Error: $CALCULATOR_BIN not found!"
  echo "Make sure you've built the project before running tests."
  exit 1
fi

# Create directories if they don't exist
mkdir -p "$TEST_DATA_DIR"
mkdir -p "$EXPECTED_DIR"
mkdir -p "$OUTPUT_DIR"

###############################################################################
# Test 1: Help/usage
###############################################################################
echo "Test 1: Help / usage"
if ! "$CALCULATOR_BIN" --help 2>&1 | grep -q "VCFX_ld_calculator: Calculate pairwise LD"; then
    echo "✗ Test failed: help - usage message not found."
    exit 1
fi
echo "✓ Test 1 passed"

###############################################################################
# Test 2: Single variant => no pairwise LD
###############################################################################
echo "Test 2: Single variant"

cat > "$TEST_DATA_DIR/single.vcf" <<EOF
##fileformat=VCFv4.2
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1
1       100     .       A       G       .       PASS    .       GT      0/0
EOF

# No --region => use entire file
OUT2=$(cat "$TEST_DATA_DIR/single.vcf" | "$CALCULATOR_BIN")
# We expect #LD_MATRIX_START plus message "No or only one variant..."
if ! echo "$OUT2" | grep -q "No or only one variant in the region => no pairwise LD."; then
    echo "✗ Test failed: single variant => expected 'no pairwise LD' message not found."
    echo "Output was:"
    echo "$OUT2"
    exit 1
fi
echo "✓ Test 2 passed"

###############################################################################
# Test 3: Two variants => produce 2x2 matrix
###############################################################################
echo "Test 3: Two variants"

# 1) Write a temporary file with literal \t
cat > "$TEST_DATA_DIR/two_variants_tmp.vcf" <<EOF
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1
1\t200\t.\tT\tC\t.\tPASS\t.\tGT\t0/1\t1/1
EOF
# 2) Convert \t to real tabs
sed 's/\\t/\t/g' "$TEST_DATA_DIR/two_variants_tmp.vcf" > "$TEST_DATA_DIR/two_variants.vcf"

OUT3=$(cat "$TEST_DATA_DIR/two_variants.vcf" | "$CALCULATOR_BIN")

# We expect two data lines => #LD_MATRIX_START => a 2x2 matrix
if ! echo "$OUT3" | grep -q "#LD_MATRIX_START"; then
    echo "✗ Test failed: two variants => missing #LD_MATRIX_START"
    echo "Output was:"
    echo "$OUT3"
    exit 1
fi
if ! echo "$OUT3" | grep -q "^1:100\s"; then
    echo "✗ Test failed: two variants => missing row for 1:100"
    echo "Output was:"
    echo "$OUT3"
    exit 1
fi
echo "✓ Test 3 passed"

###############################################################################
# Test 4: Region excludes everything => no pairwise LD
###############################################################################
echo "Test 4: Region excludes everything"

cat > "$TEST_DATA_DIR/out_of_range_tmp.vcf" <<EOF
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0
1\t200\t.\tT\tC\t.\tPASS\t.\tGT\t0/1
EOF
sed 's/\\t/\t/g' "$TEST_DATA_DIR/out_of_range_tmp.vcf" > "$TEST_DATA_DIR/out_of_range.vcf"

# region is chr1:300-400 => excludes POS=100 & 200 => no variants in region
OUT4=$(cat "$TEST_DATA_DIR/out_of_range.vcf" | "$CALCULATOR_BIN" --region 1:300-400)

if ! echo "$OUT4" | grep -q "No or only one variant in the region => no pairwise LD."; then
    echo "✗ Test failed: region excludes everything => expected 'no pairwise LD' message not found."
    echo "Output was:"
    echo "$OUT4"
    exit 1
fi
echo "✓ Test 4 passed"

###############################################################################
# Test 5: Region includes a subset
###############################################################################
echo "Test 5: Region includes a subset"

cat > "$TEST_DATA_DIR/partial_region_tmp.vcf" <<EOF
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
1\t150\t.\tA\tG\t.\tPASS\t.\tGT\t0/1
1\t250\t.\tT\tC\t.\tPASS\t.\tGT\t1/1
1\t350\t.\tG\tT\t.\tPASS\t.\tGT\t0/0
EOF
sed 's/\\t/\t/g' "$TEST_DATA_DIR/partial_region_tmp.vcf" > "$TEST_DATA_DIR/partial_region.vcf"

# We'll only keep variants at [200..300], so that includes POS=250 => 1 variant
OUT5=$(cat "$TEST_DATA_DIR/partial_region.vcf" | "$CALCULATOR_BIN" --region 1:200-300)

if ! echo "$OUT5" | grep -q "No or only one variant in the region => no pairwise LD."; then
    echo "✗ Test failed: region subset => expected 'no pairwise LD' message"
    echo "Output was:"
    echo "$OUT5"
    exit 1
fi
echo "✓ Test 5 passed"

###############################################################################
# Test 6: Missing / multi-allelic genotypes
###############################################################################
echo "Test 6: Missing or multi-allelic"

cat > "$TEST_DATA_DIR/missing_multi_tmp.vcf" <<EOF
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t100\tvarA\tA\tG\t.\tPASS\t.\tGT\t0/1\t./.\t1/1
1\t200\tvarB\tT\tC\t.\tPASS\t.\tGT\t1/2\t0/1\t0/0
1\t300\tvarC\tG\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0
EOF
sed 's/\\t/\t/g' "$TEST_DATA_DIR/missing_multi_tmp.vcf" > "$TEST_DATA_DIR/missing_multi.vcf"

OUT6=$(cat "$TEST_DATA_DIR/missing_multi.vcf" | "$CALCULATOR_BIN")

# We expect #LD_MATRIX_START and a 3x3 matrix
if ! echo "$OUT6" | grep -q "#LD_MATRIX_START"; then
    echo "✗ Test failed: missing/multi => no #LD_MATRIX_START found"
    echo "$OUT6"
    exit 1
fi
if ! echo "$OUT6" | grep -q "1:200"; then
    echo "✗ Test failed: missing/multi => variant 1:200 not found in LD matrix"
    echo "$OUT6"
    exit 1
fi
echo "✓ Test 6 passed"

###############################################################################
# Test 7: Data line before #CHROM => expect error message
###############################################################################
echo "Test 7: Data line before #CHROM"

cat > "$TEST_DATA_DIR/no_header_tmp.vcf" <<EOF
1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
1\t200\t.\tA\tC\t.\tPASS\t.\tGT\t0/1
EOF
sed 's/\\t/\t/g' "$TEST_DATA_DIR/no_header_tmp.vcf" > "$TEST_DATA_DIR/no_header.vcf"

ERR7=$(cat "$TEST_DATA_DIR/no_header.vcf" | "$CALCULATOR_BIN" 2>&1 || true)

# The code prints "Error: encountered data line before #CHROM."
if ! echo "$ERR7" | grep -q "Error: encountered data line before #CHROM."; then
    echo "✗ Test failed: data line before #CHROM => expected error message not found."
    echo "Stderr was:"
    echo "$ERR7"
    exit 1
fi
echo "✓ Test 7 passed"

echo "✅ All VCFX_ld_calculator tests passed!"
