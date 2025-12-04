#!/bin/bash

# Exit on error
set -e

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"
PHASER="$ROOT_DIR/build/src/VCFX_haplotype_phaser/VCFX_haplotype_phaser"

# Create test data directory
TEST_DIR="${SCRIPT_DIR}/data/haplotype_phaser"
mkdir -p "$TEST_DIR"
cd "$TEST_DIR"

echo "=== Testing VCFX_haplotype_phaser ==="

# Test 1: Basic functionality
echo "Test 1: Basic functionality with high LD"
cat > basic.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	100	rs1	A	G	100	PASS	.	GT	0/0	0/1	1/1
1	150	rs2	C	T	100	PASS	.	GT	0/0	0/1	1/1
1	200	rs3	G	A	100	PASS	.	GT	0/0	0/1	1/1
1	250	rs4	T	C	100	PASS	.	GT	1/1	0/1	0/0
EOF

$PHASER --ld-threshold 0.8 < basic.vcf > basic.out 2> basic.err

# Verify output
if grep -q "Block 1: 0:(1:100), 1:(1:150), 2:(1:200)" basic.out && \
   grep -q "Block 2: 3:(1:250)" basic.out; then
    echo "✓ Test 1 passed: Variants correctly grouped into blocks"
else
    echo "❌ Test 1 failed: Blocks not formed correctly"
    cat basic.out
    exit 1
fi

# Test 2: Different LD threshold
echo "Test 2: Different LD threshold"
$PHASER --ld-threshold 0.99 < basic.vcf > high_threshold.out 2> high_threshold.err

# With high threshold, should create more blocks
if grep -q "Block 1: 0:(1:100), 1:(1:150), 2:(1:200)" high_threshold.out && \
   grep -q "Block 2: 3:(1:250)" high_threshold.out; then
    echo "✓ Test 2 passed: Variants correctly grouped with high threshold"
else
    echo "❌ Test 2 failed: Blocks not formed correctly with high threshold"
    cat high_threshold.out
    exit 1
fi

# Test 3: Low LD with non-perfectly correlated variants
echo "Test 3: Low LD variants"
cat > low_ld.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3	SAMPLE4	SAMPLE5
1	100	rs1	A	G	100	PASS	.	GT	0/0	0/1	1/1	0/1	0/0
1	150	rs2	C	T	100	PASS	.	GT	0/0	0/1	1/1	0/1	0/0
1	200	rs3	G	A	100	PASS	.	GT	0/0	0/1	1/1	1/1	0/0
1	250	rs4	T	C	100	PASS	.	GT	1/1	0/1	0/0	0/0	1/1
1	300	rs5	A	C	100	PASS	.	GT	0/0	0/0	0/1	1/1	0/1
EOF

$PHASER --ld-threshold 0.5 < low_ld.vcf > low_ld.out 2> low_ld.err

# Check for more blocks with lower correlation
if grep -q "Block 1: 0:(1:100), 1:(1:150), 2:(1:200)" low_ld.out && \
   grep -q "Block 2: 3:(1:250)" low_ld.out && \
   grep -q "Block 3: 4:(1:300)" low_ld.out; then
    echo "✓ Test 3 passed: Low LD variants handled correctly"
else
    echo "❌ Test 3 failed: Low LD variants not handled correctly"
    cat low_ld.out
    exit 1
fi

# Test 4: Missing genotypes
echo "Test 4: Missing genotypes"
cat > missing.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	100	rs1	A	G	100	PASS	.	GT	0/0	0/1	1/1
1	150	rs2	C	T	100	PASS	.	GT	0/0	./.	1/1
1	200	rs3	G	A	100	PASS	.	GT	./.	0/1	./.
1	250	rs4	T	C	100	PASS	.	GT	1/1	0/1	0/0
EOF

$PHASER --ld-threshold 0.8 < missing.vcf > missing.out 2> missing.err

# Check that missing genotypes are handled correctly
if grep -q "Block 1: 0:(1:100), 1:(1:150)" missing.out && \
   grep -q "Block 2: 2:(1:200)" missing.out && \
   grep -q "Block 3: 3:(1:250)" missing.out; then
    echo "✓ Test 4 passed: Missing genotypes handled correctly"
else
    echo "❌ Test 4 failed: Missing genotypes not handled correctly"
    cat missing.out
    exit 1
fi

# Test 5: Phased genotypes
echo "Test 5: Phased genotypes"
cat > phased.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	100	rs1	A	G	100	PASS	.	GT	0|0	0|1	1|1
1	150	rs2	C	T	100	PASS	.	GT	0|0	0|1	1|1
1	200	rs3	G	A	100	PASS	.	GT	0|0	0|1	1|1
1	250	rs4	T	C	100	PASS	.	GT	1|1	0|1	0|0
EOF

$PHASER --ld-threshold 0.8 < phased.vcf > phased.out 2> phased.err

# Check that phased genotypes are handled correctly
if grep -q "Block 1: 0:(1:100), 1:(1:150), 2:(1:200)" phased.out && \
   grep -q "Block 2: 3:(1:250)" phased.out; then
    echo "✓ Test 5 passed: Phased genotypes handled correctly"
else
    echo "❌ Test 5 failed: Phased genotypes not handled correctly"
    cat phased.out
    exit 1
fi

# Test 6: Malformed VCF
echo "Test 6: Malformed VCF"
cat > malformed.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
1	100	rs1	A	G	100	PASS	.	GT	0/0	0/1
1	invalid	rs2	C	T	100	PASS	.	GT	0/0	0/1
1	200	rs3	G	A	100	PASS	.	GT	0x0	0/1
#UNEXPECTED	HEADER
1	250	rs4	T	C	100	PASS	.	GT	1/1	0/1
EOF

$PHASER --ld-threshold 0.8 < malformed.vcf > malformed.out 2> malformed.err

# Check error output and warning messages
if grep -q "Warning: invalid pos => skip" malformed.err && \
   grep -q "Block 1: 0:(1:100)" malformed.out; then
    echo "✓ Test 6 passed: Malformed VCF handled correctly"
else
    echo "❌ Test 6 failed: Malformed VCF not handled correctly"
    cat malformed.err
    cat malformed.out
    exit 1
fi

# Test 7: No CHROM header
echo "Test 7: No CHROM header"
cat > no_header.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
1	100	rs1	A	G	100	PASS	.	GT	0/0	0/1	1/1
1	150	rs2	C	T	100	PASS	.	GT	0/0	0/1	1/1
EOF

$PHASER --ld-threshold 0.8 < no_header.vcf > no_header.out 2> no_header.err

# Check error output
if grep -q "Error: no #CHROM line found" no_header.err; then
    echo "✓ Test 7 passed: No CHROM header handled correctly"
else
    echo "❌ Test 7 failed: No CHROM header not handled correctly"
    cat no_header.err
    exit 1
fi

# Test 8: Empty input
echo "Test 8: Empty input"
cat > empty.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
EOF

$PHASER --ld-threshold 0.8 < empty.vcf > empty.out 2> empty.err

# Check error output
if grep -q "Error: no variant data found" empty.err; then
    echo "✓ Test 8 passed: Empty input handled correctly"
else
    echo "❌ Test 8 failed: Empty input not handled correctly"
    cat empty.err
    exit 1
fi

# Test 9: Multiple chromosomes
echo "Test 9: Multiple chromosomes"
cat > multi_chrom.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	100	rs1	A	G	100	PASS	.	GT	0/0	0/1	1/1
1	150	rs2	C	T	100	PASS	.	GT	0/0	0/1	1/1
2	100	rs3	G	A	100	PASS	.	GT	0/0	0/1	1/1
2	150	rs4	T	C	100	PASS	.	GT	1/1	0/1	0/0
X	100	rs5	A	C	100	PASS	.	GT	0/0	0/1	1/1
EOF

$PHASER --ld-threshold 0.8 < multi_chrom.vcf > multi_chrom.out 2> multi_chrom.err

# Check that it blocks by chromosome boundaries
if grep -q "Block 1: 0:(1:100), 1:(1:150)" multi_chrom.out && \
   grep -q "Block 2: 2:(2:100), 3:(2:150)" multi_chrom.out && \
   grep -q "Block 3: 4:(X:100)" multi_chrom.out; then
    echo "✓ Test 9 passed: Multiple chromosomes handled correctly"
else
    echo "❌ Test 9 failed: Multiple chromosomes not handled correctly"
    cat multi_chrom.out
    exit 1
fi

# Test 10: Help message
echo "Test 10: Help message"
$PHASER --help > help.out 2>&1

if grep -q "VCFX_haplotype_phaser: Group variants into blocks by naive LD threshold" help.out && \
   grep -q "ld-threshold" help.out; then
    echo "✓ Test 10 passed: Help message displayed correctly"
else
    echo "❌ Test 10 failed: Help message not displayed correctly"
    cat help.out
    exit 1
fi

# Test 11: Invalid LD threshold
echo "Test 11: Invalid LD threshold"
$PHASER --ld-threshold 1.5 < basic.vcf > invalid_ld.out 2> invalid_ld.err || true

# Check for error message in stderr or usage in stdout
if grep -q "Error: invalid LD threshold" invalid_ld.err || grep -q "Usage:" invalid_ld.out; then
    echo "✓ Test 11 passed: Invalid LD threshold detected correctly"
else
    echo "❌ Test 11 failed: Invalid LD threshold not handled correctly"
    echo "stderr content:"
    cat invalid_ld.err
    echo "stdout content:"
    cat invalid_ld.out
    exit 1
fi

# Test 12: Performance with large file
echo "Test 12: Performance with large file"
> large.vcf
echo "##fileformat=VCFv4.2" >> large.vcf
echo "##source=VCFX_test" >> large.vcf
echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3	SAMPLE4	SAMPLE5" >> large.vcf

# Generate 100 variants across chromosome 1
for pos in $(seq 1000 1000 100000); do
    if [ $((RANDOM % 2)) -eq 0 ]; then
        # Highly correlated with previous variant
        echo "1	$pos	rs_$pos	A	G	100	PASS	.	GT	0/0	0/1	1/1	0/1	0/0" >> large.vcf
    else
        # Less correlated with previous variant
        echo "1	$pos	rs_$pos	A	G	100	PASS	.	GT	0/0	0/1	1/1	1/1	1/1" >> large.vcf
    fi
done

# Time the execution
time $PHASER --ld-threshold 0.8 < large.vcf > large.out 2> large.err

# Verify the tool ran without errors
if [ -s large.out ] && grep -q "Block" large.out; then
    BLOCK_COUNT=$(grep -c "Block" large.out)
    echo "✓ Test 12 passed: Large file processed successfully with $BLOCK_COUNT blocks"
else
    echo "❌ Test 12 failed: Error processing large file"
    cat large.err
    exit 1
fi

# Test 13: Streaming mode - basic functionality
echo "Test 13: Streaming mode - basic functionality"
$PHASER --streaming --ld-threshold 0.8 < basic.vcf > streaming_basic.out 2> streaming_basic.err

# Streaming output should have the same blocks as default mode
if grep -q "Block 1: 0:(1:100), 1:(1:150), 2:(1:200)" streaming_basic.out && \
   grep -q "Block 2: 3:(1:250)" streaming_basic.out; then
    echo "✓ Test 13 passed: Streaming mode produces correct blocks"
else
    echo "❌ Test 13 failed: Streaming mode block mismatch"
    cat streaming_basic.out
    exit 1
fi

# Test 14: Streaming vs default mode consistency
echo "Test 14: Streaming vs default mode output consistency"
$PHASER --ld-threshold 0.8 < multi_chrom.vcf > default_multi.out 2>/dev/null
$PHASER --streaming --ld-threshold 0.8 < multi_chrom.vcf > streaming_multi.out 2>/dev/null

# Compare block boundaries (ignoring potential whitespace differences)
default_blocks=$(grep "Block" default_multi.out | sort)
streaming_blocks=$(grep "Block" streaming_multi.out | sort)

if [ "$default_blocks" = "$streaming_blocks" ]; then
    echo "✓ Test 14 passed: Streaming output matches default mode"
else
    echo "❌ Test 14 failed: Streaming differs from default"
    echo "Default blocks:"
    echo "$default_blocks"
    echo "Streaming blocks:"
    echo "$streaming_blocks"
    exit 1
fi

# Test 15: Streaming mode with window size
echo "Test 15: Streaming mode with custom window size"
$PHASER --streaming --window 50 --ld-threshold 0.8 < low_ld.vcf > streaming_window.out 2> streaming_window.err

if grep -q "Block" streaming_window.out; then
    echo "✓ Test 15 passed: Streaming with window parameter works"
else
    echo "❌ Test 15 failed: Streaming with window did not produce output"
    cat streaming_window.err
    exit 1
fi

# Test 16: Streaming mode help shows new options
echo "Test 16: Help shows streaming options"
$PHASER --help > help_streaming.out 2>&1

if grep -q "streaming" help_streaming.out && grep -q "window" help_streaming.out; then
    echo "✓ Test 16 passed: Help shows streaming options"
else
    echo "❌ Test 16 failed: Help does not show streaming options"
    cat help_streaming.out
    exit 1
fi

echo "All tests for VCFX_haplotype_phaser passed!"
exit 0 