#!/bin/bash

# Exit on error
set -e

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"
DIFF_TOOL="$ROOT_DIR/build/src/VCFX_diff_tool/VCFX_diff_tool"

# Create test data directory
TEST_DIR="${SCRIPT_DIR}/data"
mkdir -p "$TEST_DIR"
cd "$TEST_DIR"

echo "=== Testing VCFX_diff_tool ==="

# Test 1: Basic functionality with identical files
echo "Test 1: Basic functionality with identical files"
cat > file1.vcf << 'EOF'
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	rs1	A	G	100	PASS	.
1	200	rs2	C	T	100	PASS	.
2	300	rs3	G	A	100	PASS	.
EOF

cp file1.vcf file2.vcf

output=$("$DIFF_TOOL" --file1 file1.vcf --file2 file2.vcf)
if [ -n "$output" ] && ! echo "$output" | grep -q "Variants unique to file1.vcf:"; then
    echo "❌ Test 1 failed: Expected no differences between identical files"
    echo "Output:"
    echo "$output"
    exit 1
fi
echo "✓ Test 1 passed"

# Test 2: Files with different variants
echo "Test 2: Files with different variants"
cat > file2.vcf << 'EOF'
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	rs1	A	G	100	PASS	.
1	200	rs2	C	T	100	PASS	.
2	300	rs3	G	A	100	PASS	.
2	400	rs4	T	C	100	PASS	.
3	500	rs5	A	G	100	PASS	.
EOF

output=$("$DIFF_TOOL" --file1 file1.vcf --file2 file2.vcf)
if ! echo "$output" | grep -q "Variants unique to file2.vcf:"; then
    echo "❌ Test 2 failed: Expected to find variants unique to file2"
    echo "Output:"
    echo "$output"
    exit 1
fi
if ! echo "$output" | grep -q "2:400:T:C"; then
    echo "❌ Test 2 failed: Expected to find variant 2:400:T:C"
    echo "Output:"
    echo "$output"
    exit 1
fi
if ! echo "$output" | grep -q "3:500:A:G"; then
    echo "❌ Test 2 failed: Expected to find variant 3:500:A:G"
    echo "Output:"
    echo "$output"
    exit 1
fi
echo "✓ Test 2 passed"

# Test 3: Multi-allelic variant handling
echo "Test 3: Multi-allelic variant handling"
cat > file1.vcf << 'EOF'
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	rs1	A	G,T	100	PASS	.
EOF

cat > file2.vcf << 'EOF'
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	rs1	A	T,G	100	PASS	.
EOF

output=$("$DIFF_TOOL" --file1 file1.vcf --file2 file2.vcf)
if [ -n "$output" ] && ! echo "$output" | grep -q "Variants unique to file1.vcf:"; then
    echo "❌ Test 3 failed: Expected no differences between files with same multi-allelic variants"
    echo "Output:"
    echo "$output"
    exit 1
fi
echo "✓ Test 3 passed"

# Test 4: Header handling
echo "Test 4: Header handling"
cat > file1.vcf << 'EOF'
##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	rs1	A	G	100	PASS	.
EOF

cat > file2.vcf << 'EOF'
##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##contig=<ID=2,length=2000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	rs1	A	G	100	PASS	.
EOF

output=$("$DIFF_TOOL" --file1 file1.vcf --file2 file2.vcf)
if [ -n "$output" ] && ! echo "$output" | grep -q "Variants unique to file1.vcf:"; then
    echo "❌ Test 4 failed: Expected no differences between files with different headers"
    echo "Output:"
    echo "$output"
    exit 1
fi
echo "✓ Test 4 passed"

# Test 5: Invalid VCF line handling
echo "Test 5: Invalid VCF line handling"
cat > file1.vcf << 'EOF'
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	rs1	A	G	100	PASS	.
invalid_line
1	200	rs2	C	T	100	PASS	.
EOF

cat > file2.vcf << 'EOF'
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	rs1	A	G	100	PASS	.
1	200	rs2	C	T	100	PASS	.
EOF

output=$("$DIFF_TOOL" --file1 file1.vcf --file2 file2.vcf)
if [ -n "$output" ] && ! echo "$output" | grep -q "Variants unique to file1.vcf:"; then
    echo "❌ Test 5 failed: Expected no differences between files with valid variants"
    echo "Output:"
    echo "$output"
    exit 1
fi
echo "✓ Test 5 passed"

# Test 6: Missing file handling
echo "Test 6: Missing file handling"
# Temporarily disable exit on error for this test
set +e
"$DIFF_TOOL" --file1 nonexistent.vcf --file2 file2.vcf > output.txt 2>&1
exit_code=$?
set -e

# Check that the exit code is non-zero
if [ $exit_code -eq 0 ]; then
    echo "❌ Test 6 failed: Expected non-zero exit code for missing file"
    exit 1
fi

# Check for the expected error message
if ! grep -q "Error: Unable to open file nonexistent.vcf" output.txt; then
    echo "❌ Test 6 failed: Expected error message for missing file"
    echo "Output:"
    cat output.txt
    exit 1
fi
echo "✓ Test 6 passed"

# Test 7: Help message
echo "Test 7: Help message"
output=$("$DIFF_TOOL" --help)
if ! echo "$output" | grep -q "VCFX_diff_tool: Compare two VCF files and identify differences"; then
    echo "❌ Test 7 failed: Expected help message"
    echo "Output:"
    echo "$output"
    exit 1
fi
echo "✓ Test 7 passed"

# Test 8: No arguments handling
echo "Test 8: No arguments handling"
# Temporarily disable exit on error for this test
set +e
"$DIFF_TOOL" > /dev/null 2>&1
exit_code=$?
set -e

if [ $exit_code -eq 0 ]; then
    echo "❌ Test 8 failed: Expected non-zero exit code for no arguments"
    exit 1
fi
echo "✓ Test 8 passed"

# Test 9: Performance with large files
echo "Test 9: Performance with large files"
# Generate large test files with 1000 variants each
> large1.vcf  # Clear or create the file
> large2.vcf  # Clear or create the file

for i in {1..1000}; do
    echo "1       $i       rs$i    A       G       100     PASS    ." >> large1.vcf
done

for i in {1..1000}; do
    if [ $i -le 500 ]; then
        echo "1       $i       rs$i    A       G       100     PASS    ." >> large2.vcf
    else
        echo "1       $i       rs$i    A       T       100     PASS    ." >> large2.vcf
    fi
done

time output=$("$DIFF_TOOL" --file1 large1.vcf --file2 large2.vcf)
if ! echo "$output" | grep -q "Variants unique to large1.vcf:"; then
    echo "❌ Test 9 failed: Expected to find variants unique to large1"
    echo "Output:"
    echo "$output"
    exit 1
fi
if ! echo "$output" | grep -q "Variants unique to large2.vcf:"; then
    echo "❌ Test 9 failed: Expected to find variants unique to large2"
    echo "Output:"
    echo "$output"
    exit 1
fi
echo "✓ Test 9 passed"

# Test 10: Streaming mode with --assume-sorted
echo "Test 10: Streaming mode with --assume-sorted"
cat > sorted1.vcf << EOF
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	rs1	A	G	100	PASS	.
1	200	rs2	C	T	100	PASS	.
1	300	rs3	G	A	100	PASS	.
EOF

cat > sorted2.vcf << EOF
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	rs1	A	G	100	PASS	.
1	200	rs2	C	T	100	PASS	.
1	400	rs4	T	C	100	PASS	.
EOF

output=$("$DIFF_TOOL" --file1 sorted1.vcf --file2 sorted2.vcf --assume-sorted)
if ! echo "$output" | grep -q "1:300:G:A"; then
    echo "❌ Test 10 failed: Expected to find 1:300:G:A unique to sorted1"
    echo "Output:"
    echo "$output"
    exit 1
fi
if ! echo "$output" | grep -q "1:400:T:C"; then
    echo "❌ Test 10 failed: Expected to find 1:400:T:C unique to sorted2"
    echo "Output:"
    echo "$output"
    exit 1
fi
echo "✓ Test 10 passed"

# Test 11: Streaming mode produces same results as default mode
echo "Test 11: Streaming mode produces same results as default mode"
"$DIFF_TOOL" --file1 sorted1.vcf --file2 sorted2.vcf > default_output.txt
"$DIFF_TOOL" --file1 sorted1.vcf --file2 sorted2.vcf --assume-sorted > streaming_output.txt
# Sort both outputs for comparison (streaming preserves order, default doesn't)
sort default_output.txt > default_sorted.txt
sort streaming_output.txt > streaming_sorted.txt
if ! diff default_sorted.txt streaming_sorted.txt > /dev/null; then
    echo "❌ Test 11 failed: Streaming and default mode produce different results"
    echo "Default:"
    cat default_output.txt
    echo "Streaming:"
    cat streaming_output.txt
    exit 1
fi
echo "✓ Test 11 passed"

# Test 12: Help shows new options
echo "Test 12: Help shows new --assume-sorted and --natural-chr options"
output=$("$DIFF_TOOL" --help)
if ! echo "$output" | grep -q "\-\-assume-sorted"; then
    echo "❌ Test 12 failed: Expected help to show --assume-sorted option"
    exit 1
fi
if ! echo "$output" | grep -q "\-\-natural-chr"; then
    echo "❌ Test 12 failed: Expected help to show --natural-chr option"
    exit 1
fi
echo "✓ Test 12 passed"

# Test 13: Empty files handling in streaming mode
echo "Test 13: Empty files handling in streaming mode"
cat > empty1.vcf << EOF
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
EOF

cat > nonempty.vcf << EOF
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	rs1	A	G	100	PASS	.
EOF

output=$("$DIFF_TOOL" --file1 empty1.vcf --file2 nonempty.vcf --assume-sorted)
if ! echo "$output" | grep -q "1:100:A:G"; then
    echo "❌ Test 13 failed: Expected to find variant unique to nonempty file"
    echo "Output:"
    echo "$output"
    exit 1
fi
echo "✓ Test 13 passed"

echo "All tests for VCFX_diff_tool passed!" 