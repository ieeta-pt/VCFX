#!/usr/bin/env bash

set -e
set -o pipefail

###############################################################################
# Paths & Setup
###############################################################################
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"

MERGER_BIN="${ROOT_DIR}/build/src/VCFX_merger/VCFX_merger"
TEST_DATA_DIR="${SCRIPT_DIR}/data/merger"
EXPECTED_DIR="${SCRIPT_DIR}/expected/merger"
OUTPUT_DIR="${SCRIPT_DIR}/out/merger"

mkdir -p "$TEST_DATA_DIR"
mkdir -p "$EXPECTED_DIR"
mkdir -p "$OUTPUT_DIR"

###############################################################################
# Function: makeTabs
# Replaces literal "\t" with actual ASCII tab characters
###############################################################################
function makeTabs() {
    local inFile="$1"
    local outFile="$2"
    sed 's/\\t/\t/g' "$inFile" > "$outFile"
}

# Verify the merger binary
if [ ! -f "$MERGER_BIN" ]; then
  echo "Error: $MERGER_BIN not found!"
  echo "Make sure you've built the project before running tests."
  exit 1
fi

echo "=== Testing VCFX_merger ==="

# Clean up old VCFs before we begin
rm -f "${TEST_DATA_DIR}"/*.vcf

###############################################################################
# Test 1: Help/usage
###############################################################################
echo "Test 1: Help / usage"
HELP_OUT=$("$MERGER_BIN" --help 2>&1 || true)
if ! echo "$HELP_OUT" | grep -q "VCFX_merger: Merge multiple VCF files"; then
    echo "✗ Test failed: help - usage message not found."
    echo "Output was:"
    echo "$HELP_OUT"
    exit 1
fi
echo "✓ Test 1 passed"

###############################################################################
# Test 2: No --merge => usage
###############################################################################
echo "Test 2: No --merge argument => usage"
NO_MERGE_OUT=$("$MERGER_BIN" 2>&1 || true)
if ! echo "$NO_MERGE_OUT" | grep -q "VCFX_merger: Merge multiple VCF files"; then
    echo "✗ Test failed: no --merge => expected help usage not found"
    echo "Output was:"
    echo "$NO_MERGE_OUT"
    exit 1
fi
echo "✓ Test 2 passed"

###############################################################################
# Test 3: Merge a single file
###############################################################################
echo "Test 3: Single file"

cat > "$TEST_DATA_DIR/single_tmp.vcf" <<EOF
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\t.\tA\tG\t.\tPASS\t.
1\t200\t.\tT\tC\t.\tPASS\t.
EOF

makeTabs "$TEST_DATA_DIR/single_tmp.vcf" "$TEST_DATA_DIR/single.vcf"
SINGLE_OUT=$("$MERGER_BIN" --merge "$TEST_DATA_DIR/single.vcf")

if ! echo "$SINGLE_OUT" | grep -q "##fileformat=VCFv4.2"; then
    echo "✗ Test failed: single => missing header"
    exit 1
fi
if ! echo "$SINGLE_OUT" | grep -q "^1[[:space:]]100"; then
    echo "✗ Test failed: single => missing variant pos=100"
    echo "$SINGLE_OUT"
    exit 1
fi
if ! echo "$SINGLE_OUT" | grep -q "^1[[:space:]]200"; then
    echo "✗ Test failed: single => missing variant pos=200"
    echo "$SINGLE_OUT"
    exit 1
fi
echo "✓ Test 3 passed"

###############################################################################
# Test 4: Merge two files => sorted by chromosome/pos
###############################################################################
echo "Test 4: Two files"

# Remove old fileA/fileB
rm -f "$TEST_DATA_DIR"/fileA*.vcf "$TEST_DATA_DIR"/fileB*.vcf

# fileA_tmp.vcf (positions 100,300) with \t placeholders
cat > "$TEST_DATA_DIR/fileA_tmp.vcf" <<'EOF'
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\tvarA\tA\tG\t.\tPASS\t.
1\t300\tvarC\tG\tT\t.\tPASS\t.
EOF
makeTabs "$TEST_DATA_DIR/fileA_tmp.vcf" "$TEST_DATA_DIR/fileA.vcf"

# fileB_tmp.vcf (positions 200,150) with \t placeholders
cat > "$TEST_DATA_DIR/fileB_tmp.vcf" <<'EOF'
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t200\tvarB\tT\tC\t.\tPASS\t.
1\t150\tvarM\tA\tC\t.\tPASS\t.
EOF
makeTabs "$TEST_DATA_DIR/fileB_tmp.vcf" "$TEST_DATA_DIR/fileB.vcf"

MERGE_OUT=$("$MERGER_BIN" --merge "$TEST_DATA_DIR/fileA.vcf,$TEST_DATA_DIR/fileB.vcf")
echo "$MERGE_OUT" > "$OUTPUT_DIR/merged_output.vcf"

EXPECTED_ORDER=("100" "150" "200" "300")
i=0

# We read the merged lines, skipping meta lines, expecting 4 total in sorted order
while IFS=$'\t' read -r line; do
    [[ "$line" =~ ^## ]] && continue
    if [[ "$line" =~ ^#CHROM ]]; then
        continue
    fi
    IFS=$'\t' read -r -a fields <<< "$line"
    pos="${fields[1]}"
    if [[ $i -ge 4 ]]; then
        echo "✗ Test failed: More than 4 data lines found!"
        exit 1
    fi
    if [[ "$pos" != "${EXPECTED_ORDER[$i]}" ]]; then
        echo "✗ Test failed: out of order or missing. Expected '${EXPECTED_ORDER[$i]}', got '$pos'"
        echo "Merged output was:"
        echo "$MERGE_OUT"
        exit 1
    fi
    i=$(( i + 1 ))
done < <(echo "$MERGE_OUT")

if [[ $i -ne 4 ]]; then
    echo "✗ Test failed: Expected 4 variants, found $i"
    echo "$MERGE_OUT"
    exit 1
fi
echo "✓ Test 4 passed"

###############################################################################
# Test 5: Nonexistent input => error
###############################################################################
echo "Test 5: Nonexistent input file"

BAD_OUT=$("$MERGER_BIN" --merge "does_not_exist.vcf" 2>&1 || true)
if ! echo "$BAD_OUT" | grep -q "Failed to open file: does_not_exist.vcf"; then
    echo "✗ Test failed: nonexistent => expected error"
    echo "$BAD_OUT"
    exit 1
fi
echo "✓ Test 5 passed"

###############################################################################
echo "✅ All tests for VCFX_merger passed successfully!"
