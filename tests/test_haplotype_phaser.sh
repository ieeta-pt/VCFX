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

# Test 17: File input mode with -i flag
echo "Test 17: File input mode with -i flag"
$PHASER -i basic.vcf > file_input.out 2> file_input.err

if grep -q "Block 1: 0:(1:100), 1:(1:150), 2:(1:200)" file_input.out && \
   grep -q "Block 2: 3:(1:250)" file_input.out; then
    echo "✓ Test 17 passed: File input mode works correctly"
else
    echo "❌ Test 17 failed: File input mode produced incorrect output"
    cat file_input.out
    exit 1
fi

# Test 18: Positional file argument
echo "Test 18: Positional file argument"
$PHASER basic.vcf > positional.out 2> positional.err

if grep -q "Block 1: 0:(1:100), 1:(1:150), 2:(1:200)" positional.out && \
   grep -q "Block 2: 3:(1:250)" positional.out; then
    echo "✓ Test 18 passed: Positional file argument works"
else
    echo "❌ Test 18 failed: Positional argument mode failed"
    cat positional.out
    exit 1
fi

# Test 19: Quiet mode suppresses warnings
echo "Test 19: Quiet mode with -q flag"
$PHASER -q < malformed.vcf > quiet.out 2> quiet.err

# In quiet mode, no warnings should appear
if ! grep -q "Warning:" quiet.err; then
    echo "✓ Test 19 passed: Quiet mode suppresses warnings"
else
    echo "❌ Test 19 failed: Quiet mode did not suppress warnings"
    cat quiet.err
    exit 1
fi

# Test 20: Output equivalence (stdin vs mmap)
echo "Test 20: Output equivalence (stdin vs mmap)"
$PHASER < multi_chrom.vcf > stdin.out 2>/dev/null
$PHASER -i multi_chrom.vcf > mmap.out 2>/dev/null

# Compare the block outputs
stdin_blocks=$(grep "Block" stdin.out | sort)
mmap_blocks=$(grep "Block" mmap.out | sort)

if [ "$stdin_blocks" = "$mmap_blocks" ]; then
    echo "✓ Test 20 passed: stdin and mmap produce identical output"
else
    echo "❌ Test 20 failed: stdin and mmap produce different output"
    echo "stdin blocks:"
    echo "$stdin_blocks"
    echo "mmap blocks:"
    echo "$mmap_blocks"
    exit 1
fi

# Test 21: File input with streaming mode
echo "Test 21: File input with streaming mode"
$PHASER --streaming -i basic.vcf > file_streaming.out 2> file_streaming.err

if grep -q "Block 1: 0:(1:100), 1:(1:150), 2:(1:200)" file_streaming.out && \
   grep -q "Block 2: 3:(1:250)" file_streaming.out; then
    echo "✓ Test 21 passed: File input with streaming mode works"
else
    echo "❌ Test 21 failed: File input streaming produced incorrect output"
    cat file_streaming.out
    exit 1
fi

# Test 22: Performance comparison (stdin vs mmap)
echo "Test 22: Performance comparison (stdin vs mmap)"
echo "  Stdin mode timing:"
time $PHASER < large.vcf > /dev/null 2>&1

echo "  Mmap mode timing:"
time $PHASER -i large.vcf > /dev/null 2>&1

echo "✓ Test 22 passed: Performance test completed"

# Test 23: Help shows new options
echo "Test 23: Help shows new -i and -q options"
$PHASER --help > help_new.out 2>&1

if grep -q "\-i, --input" help_new.out && grep -q "\-q, --quiet" help_new.out; then
    echo "✓ Test 23 passed: Help shows -i/--input and -q/--quiet options"
else
    echo "❌ Test 23 failed: Help does not show new options"
    cat help_new.out
    exit 1
fi

# Test 24: Streaming mode with very small window (circular buffer edge case)
echo "Test 24: Streaming mode with window size 2 (circular buffer edge case)"
cat > circular_test.vcf << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
1	100	.	A	G	.	.	.	GT	0/0	0/1	1/1
1	200	.	C	T	.	.	.	GT	0/0	0/1	1/1
1	300	.	G	A	.	.	.	GT	0/0	0/1	1/1
1	400	.	T	C	.	.	.	GT	0/0	0/1	1/1
1	500	.	A	T	.	.	.	GT	0/0	0/1	1/1
EOF

$PHASER --streaming --window 2 -i circular_test.vcf > circular.out 2> circular.err

# Should produce blocks despite tiny window
if grep -q "Block" circular.out; then
    echo "✓ Test 24 passed: Small window circular buffer works"
else
    echo "❌ Test 24 failed: Circular buffer with small window failed"
    cat circular.out
    cat circular.err
    exit 1
fi

# Test 25: Streaming mode output consistency (stdin vs mmap)
echo "Test 25: Streaming mode output equivalence (stdin vs mmap)"
$PHASER --streaming < multi_chrom.vcf > streaming_stdin.out 2>/dev/null
$PHASER --streaming -i multi_chrom.vcf > streaming_mmap.out 2>/dev/null

stdin_blocks=$(grep "Block" streaming_stdin.out | sort)
mmap_blocks=$(grep "Block" streaming_mmap.out | sort)

if [ "$stdin_blocks" = "$mmap_blocks" ]; then
    echo "✓ Test 25 passed: Streaming stdin and mmap produce identical output"
else
    echo "❌ Test 25 failed: Streaming stdin and mmap differ"
    echo "stdin:"
    echo "$stdin_blocks"
    echo "mmap:"
    echo "$mmap_blocks"
    exit 1
fi

# Test 26: Empty file handling with mmap
echo "Test 26: Empty file handling with mmap"
cat > empty_mmap.vcf << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
EOF

$PHASER -q -i empty_mmap.vcf > empty_mmap.out 2> empty_mmap.err

# Should handle gracefully without crashing
echo "✓ Test 26 passed: Empty file with mmap handled gracefully"

# Test 27: Large window size (boundary test)
echo "Test 27: Window size larger than variant count"
$PHASER --streaming --window 10000 -i basic.vcf > large_window.out 2> large_window.err

if grep -q "Block 1: 0:(1:100), 1:(1:150), 2:(1:200)" large_window.out && \
   grep -q "Block 2: 3:(1:250)" large_window.out; then
    echo "✓ Test 27 passed: Large window size handled correctly"
else
    echo "❌ Test 27 failed: Large window size failed"
    cat large_window.out
    exit 1
fi

echo "All tests for VCFX_haplotype_phaser passed!"
exit 0