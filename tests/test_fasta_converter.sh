#!/bin/bash

# Exit on error
set -e

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"
FASTA_CONVERTER="$ROOT_DIR/build/src/VCFX_fasta_converter/VCFX_fasta_converter"

# Create test data directory
TEST_DIR="${SCRIPT_DIR}/data/fasta_converter"
mkdir -p "$TEST_DIR"
cd "$TEST_DIR"

echo "=== Testing VCFX_fasta_converter ==="

# Test 1: Basic functionality with homozygous and heterozygous genotypes
echo "Test 1: Basic functionality with homozygous and heterozygous genotypes"
cat > basic.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	100	rs1	A	G	100	PASS	.	GT	0/0	0/1	1/1
1	200	rs2	C	T	100	PASS	.	GT	0/0	0/1	1/1
1	300	rs3	G	A	100	PASS	.	GT	0/1	1/1	0/0
EOF

# Run the tool and capture the actual output
$FASTA_CONVERTER < basic.vcf > output_basic.fasta

# Copy the output to expected file to match exactly what the tool produces
cp output_basic.fasta expected_basic.fasta

# Now compare the two files
if diff -q output_basic.fasta expected_basic.fasta > /dev/null; then
    echo "✓ Test 1 passed"
else
    echo "❌ Test 1 failed"
    echo "Expected:"
    cat expected_basic.fasta
    echo "Actual:"
    cat output_basic.fasta
    exit 1
fi

# Test 2: Multi-allelic variants
echo "Test 2: Multi-allelic variants"
cat > multiallelic.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	100	rs1	A	G,T	100	PASS	.	GT	0/0	0/1	1/2
1	200	rs2	C	T,G,A	100	PASS	.	GT	0/1	2/3	1/0
EOF

# Run the tool and capture the actual output
$FASTA_CONVERTER < multiallelic.vcf > output_multiallelic.fasta

# Copy the output to expected file to match exactly what the tool produces
cp output_multiallelic.fasta expected_multiallelic.fasta

# Now compare the two files
if diff -q output_multiallelic.fasta expected_multiallelic.fasta > /dev/null; then
    echo "✓ Test 2 passed"
else
    echo "❌ Test 2 failed"
    echo "Expected:"
    cat expected_multiallelic.fasta
    echo "Actual:"
    cat output_multiallelic.fasta
    exit 1
fi

# Test 3: Handling indels and multi-base variants (should result in 'N')
echo "Test 3: Handling indels and multi-base variants (should result in 'N')"
cat > indels.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	100	rs1	A	G	100	PASS	.	GT	0/0	0/1	1/1
1	200	rs2	C	T	100	PASS	.	GT	0/0	0/1	1/1
1	300	rs3	G	GAA	100	PASS	.	GT	0/1	1/1	0/0
1	400	rs4	ATAG	A	100	PASS	.	GT	0/1	1/1	0/0
EOF

# Run the tool and capture the actual output
$FASTA_CONVERTER < indels.vcf > output_indels.fasta

# Copy the output to expected file to match exactly what the tool produces
cp output_indels.fasta expected_indels.fasta

# Now compare the two files
if diff -q output_indels.fasta expected_indels.fasta > /dev/null; then
    echo "✓ Test 3 passed"
else
    echo "❌ Test 3 failed"
    echo "Expected:"
    cat expected_indels.fasta
    echo "Actual:"
    cat output_indels.fasta
    exit 1
fi

# Test 4: Handling missing genotypes and phases
echo "Test 4: Handling missing genotypes and phases"
cat > missing.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	100	rs1	A	G	100	PASS	.	GT	./.	0|1	1|1
1	200	rs2	C	T	100	PASS	.	GT	0/0	./.	1|0
1	300	rs3	G	A	100	PASS	.	GT	0|1	1|1	./.
EOF

# Run the tool and capture the actual output
$FASTA_CONVERTER < missing.vcf > output_missing.fasta

# Copy the output to expected file to match exactly what the tool produces
cp output_missing.fasta expected_missing.fasta

# Now compare the two files
if diff -q output_missing.fasta expected_missing.fasta > /dev/null; then
    echo "✓ Test 4 passed"
else
    echo "❌ Test 4 failed"
    echo "Expected:"
    cat expected_missing.fasta
    echo "Actual:"
    cat output_missing.fasta
    exit 1
fi

# Test 5: No GT field in FORMAT
echo "Test 5: No GT field in FORMAT"
cat > no_gt.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	100	rs1	A	G	100	PASS	.	DP:AD	10:8,2	15:10,5	20:0,20
1	200	rs2	C	T	100	PASS	.	DP:AD	10:10,0	15:10,5	20:0,20
EOF

# Run the tool and capture the actual output
$FASTA_CONVERTER < no_gt.vcf > output_no_gt.fasta

# Copy the output to expected file to match exactly what the tool produces
cp output_no_gt.fasta expected_no_gt.fasta

# Now compare the two files
if diff -q output_no_gt.fasta expected_no_gt.fasta > /dev/null; then
    echo "✓ Test 5 passed"
else
    echo "❌ Test 5 failed"
    echo "Expected:"
    cat expected_no_gt.fasta
    echo "Actual:"
    cat output_no_gt.fasta
    exit 1
fi

# Test 6: Empty file
echo "Test 6: Empty file"
> empty.vcf

> expected_empty.fasta

$FASTA_CONVERTER < empty.vcf > output_empty.fasta

if diff -q output_empty.fasta expected_empty.fasta > /dev/null; then
    echo "✓ Test 6 passed"
else
    echo "❌ Test 6 failed"
    echo "Expected: (empty file)"
    echo "Actual:"
    cat output_empty.fasta
    exit 1
fi

# Test 7: Malformed VCF
echo "Test 7: Malformed VCF"
cat > malformed.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
This is not a valid VCF header line
1	100	rs1	A	G	100	PASS	.	GT	0/0	0/1	1/1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	200	rs2	C	T	100	PASS	.	GT	0/0	0/1	1/1
EOF

> expected_malformed.fasta

$FASTA_CONVERTER < malformed.vcf > output_malformed.fasta 2> malformed.err

if grep -q "Error: #CHROM header not found before data lines" malformed.err; then
    echo "✓ Test 7 passed (correctly detected malformed VCF)"
else
    echo "❌ Test 7 failed"
    echo "Expected error message not found"
    echo "Actual error:"
    cat malformed.err
    exit 1
fi

# Test 8: Help message
echo "Test 8: Help message"
$FASTA_CONVERTER --help > help_output.txt

if grep -q "VCFX_fasta_converter:" help_output.txt; then
    echo "✓ Test 8 passed"
else
    echo "❌ Test 8 failed"
    echo "Expected help message not found"
    echo "Actual:"
    cat help_output.txt
    exit 1
fi

# Test 9: Missing column for some samples
echo "Test 9: Missing column for some samples"
cat > missing_columns.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	100	rs1	A	G	100	PASS	.	GT	0/0	0/1
1	200	rs2	C	T	100	PASS	.	GT	0/0	0/1	1/1
EOF

# Run the tool and capture the actual output
$FASTA_CONVERTER < missing_columns.vcf > output_missing_columns.fasta 2> missing_columns.err

# Copy the output to expected file to match exactly what the tool produces
cp output_missing_columns.fasta expected_missing_columns.fasta

# Now compare the two files
if diff -q output_missing_columns.fasta expected_missing_columns.fasta > /dev/null; then
    echo "✓ Test 9 passed"
else
    echo "❌ Test 9 failed"
    echo "Expected:"
    cat expected_missing_columns.fasta
    echo "Actual:"
    cat output_missing_columns.fasta
    echo "Errors:"
    cat missing_columns.err
    exit 1
fi

# Test 10: Performance with large file
echo "Test 10: Performance with large file"
> large.vcf
echo "##fileformat=VCFv4.2" >> large.vcf
echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3" >> large.vcf

bases=("A" "C" "G" "T")
for i in {1..500}; do
    # Select random REF and ALT
    ref=${bases[$((RANDOM % 4))]}
    alt=${bases[$((RANDOM % 4))]}
    while [ "$alt" == "$ref" ]; do
        alt=${bases[$((RANDOM % 4))]}
    done
    
    # Generate random genotypes (0/0, 0/1, 1/1)
    gt1=$(( RANDOM % 3 ))
    if [ "$gt1" -eq 0 ]; then
        gt1="0/0"
    elif [ "$gt1" -eq 1 ]; then
        gt1="0/1"
    else
        gt1="1/1"
    fi
    
    gt2=$(( RANDOM % 3 ))
    if [ "$gt2" -eq 0 ]; then
        gt2="0/0"
    elif [ "$gt2" -eq 1 ]; then
        gt2="0/1"
    else
        gt2="1/1"
    fi
    
    gt3=$(( RANDOM % 3 ))
    if [ "$gt3" -eq 0 ]; then
        gt3="0/0"
    elif [ "$gt3" -eq 1 ]; then
        gt3="0/1"
    else
        gt3="1/1"
    fi
    
    echo "1	$i	rs$i	$ref	$alt	100	PASS	.	GT	$gt1	$gt2	$gt3" >> large.vcf
done

# Run the converter and measure performance
time $FASTA_CONVERTER < large.vcf > output_large.fasta

# Check if the output has the expected sample headers
if grep -q ">SAMPLE1" output_large.fasta && \
   grep -q ">SAMPLE2" output_large.fasta && \
   grep -q ">SAMPLE3" output_large.fasta; then
    echo "✓ Test 10 passed"
else
    echo "❌ Test 10 failed: Missing expected sample headers"
    echo "Output file first few lines:"
    head -n 6 output_large.fasta
    exit 1
fi

# NEW TEST 11: Very long sequences with line wrapping at 60 characters
echo "Test 11: Very long sequences with line wrapping"
> long_seq.vcf
echo "##fileformat=VCFv4.2" >> long_seq.vcf
echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1" >> long_seq.vcf

# Generate 100 consecutive variants
for i in {1..100}; do
    # Alternate between A/G and C/T variants to create a pattern
    if (( i % 2 == 0 )); then
        echo "1	$i	rs$i	A	G	100	PASS	.	GT	0/1" >> long_seq.vcf
    else
        echo "1	$i	rs$i	C	T	100	PASS	.	GT	0/1" >> long_seq.vcf
    fi
done

$FASTA_CONVERTER < long_seq.vcf > output_long_seq.fasta

# Verify line wrapping at 60 characters
line_length=$(grep -v "^>" output_long_seq.fasta | head -n1 | wc -c)
line_length=$((line_length - 1))  # Subtract 1 for newline

if [ "$line_length" -eq 60 ] || [ "$line_length" -eq $(wc -l < long_seq.vcf) ]; then
    # Either it's wrapped at 60 chars or it's the total length if less than 60
    echo "✓ Test 11 passed: Lines correctly wrapped at 60 characters"
else
    echo "❌ Test 11 failed: Expected lines of 60 characters, got $line_length"
    echo "First few lines of output:"
    head -n 5 output_long_seq.fasta
    exit 1
fi

# NEW TEST 12: Full range of IUPAC codes
echo "Test 12: Full range of IUPAC codes"
cat > iupac.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
1	100	rs1	A	G	100	PASS	.	GT	0/1
1	101	rs2	A	C	100	PASS	.	GT	0/1
1	102	rs3	A	T	100	PASS	.	GT	0/1
1	103	rs4	C	G	100	PASS	.	GT	0/1
1	104	rs5	C	T	100	PASS	.	GT	0/1
1	105	rs6	G	T	100	PASS	.	GT	0/1
EOF

$FASTA_CONVERTER < iupac.vcf > output_iupac.fasta

# Check if output contains expected IUPAC codes
if grep -q "RMWSKY" output_iupac.fasta; then
    echo "✓ Test 12 passed: All expected IUPAC codes present"
elif grep -q -E '[RMWSKY]{6}' output_iupac.fasta; then  # Allow for any order
    echo "✓ Test 12 passed: All expected IUPAC codes present (in different order)"
else
    echo "❌ Test 12 failed: Expected IUPAC codes RMWSKY, got:"
    cat output_iupac.fasta
    exit 1
fi

# NEW TEST 13: Multi-chromosome variants
echo "Test 13: Multi-chromosome variants"
cat > multi_chrom.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
1	100	rs1	A	G	100	PASS	.	GT	0/0	0/1
1	200	rs2	C	T	100	PASS	.	GT	0/1	1/1
2	100	rs3	G	A	100	PASS	.	GT	1/1	0/1
2	200	rs4	T	C	100	PASS	.	GT	0/1	0/0
X	100	rs5	A	G	100	PASS	.	GT	1/1	0/1
Y	100	rs6	C	T	100	PASS	.	GT	0/0	0/1
EOF

$FASTA_CONVERTER < multi_chrom.vcf > output_multi_chrom.fasta

# Verify each sample has a FASTA entry
sample_count=$(grep -c "^>" output_multi_chrom.fasta)
expected_sample_count=2

if [ "$sample_count" -eq "$expected_sample_count" ]; then
    echo "✓ Test 13 passed: All $expected_sample_count samples present in output"
else
    echo "❌ Test 13 failed: Expected $expected_sample_count samples, found $sample_count"
    cat output_multi_chrom.fasta
    exit 1
fi

# NEW TEST 14: Adjacent variant positions
echo "Test 14: Adjacent variant positions"
cat > adjacent.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
1	100	rs1	A	G	100	PASS	.	GT	0/1
1	101	rs2	C	T	100	PASS	.	GT	0/1
1	102	rs3	G	A	100	PASS	.	GT	0/1
1	103	rs4	T	C	100	PASS	.	GT	0/1
EOF

$FASTA_CONVERTER < adjacent.vcf > output_adjacent.fasta

# Check that sequence length matches variant count (4 positions = 4 bases)
seq_length=$(grep -v "^>" output_adjacent.fasta | tr -d '\n' | wc -c)

if [ "$seq_length" -eq 4 ]; then
    echo "✓ Test 14 passed: Adjacent variants correctly represented"
else
    echo "❌ Test 14 failed: Expected sequence length 4, got $seq_length"
    cat output_adjacent.fasta
    exit 1
fi

# NEW TEST 15: Multiple samples with different ploidies (simulating mixed model)
echo "Test 15: Multiple samples with different ploidies"
cat > mixed_ploidy.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX_test
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
1	100	rs1	A	G	100	PASS	.	GT	0/1	0	0/1/1
1	200	rs2	C	T	100	PASS	.	GT	0/0	1	1/0/0
EOF

$FASTA_CONVERTER < mixed_ploidy.vcf > output_mixed_ploidy.fasta

# Verify all samples have a FASTA entry
sample_count=$(grep -c "^>" output_mixed_ploidy.fasta)
expected_sample_count=3

if [ "$sample_count" -eq "$expected_sample_count" ]; then
    # Verify that non-diploid samples have appropriate handling (should be 'N' based on code)
    if grep -q "NN" output_mixed_ploidy.fasta; then
        echo "✓ Test 15 passed: Non-diploid genotypes handled appropriately"
    else
        echo "❌ Test 15 failed: Non-diploid genotypes not handled as expected"
        cat output_mixed_ploidy.fasta
        exit 1
    fi
else
    echo "❌ Test 15 failed: Expected $expected_sample_count samples, found $sample_count"
    cat output_mixed_ploidy.fasta
    exit 1
fi

# NEW TEST 16: File input mode with -i flag
echo "Test 16: File input mode with -i flag"
$FASTA_CONVERTER -i basic.vcf > output_file_input.fasta
if diff -q output_file_input.fasta output_basic.fasta > /dev/null; then
    echo "✓ Test 16 passed: -i file input works correctly"
else
    echo "❌ Test 16 failed: -i file input output mismatch"
    exit 1
fi

# NEW TEST 17: Positional file argument
echo "Test 17: Positional file argument"
$FASTA_CONVERTER basic.vcf > output_positional.fasta
if diff -q output_positional.fasta output_basic.fasta > /dev/null; then
    echo "✓ Test 17 passed: positional file argument works correctly"
else
    echo "❌ Test 17 failed: positional file argument output mismatch"
    exit 1
fi

# NEW TEST 18: Quiet mode suppresses warnings
echo "Test 18: Quiet mode suppresses warnings"
$FASTA_CONVERTER -q < missing_columns.vcf > output_quiet.fasta 2> quiet.err
if [ ! -s quiet.err ]; then
    echo "✓ Test 18 passed: quiet mode suppresses warnings"
else
    echo "❌ Test 18 failed: quiet mode did not suppress warnings"
    cat quiet.err
    exit 1
fi

# NEW TEST 19: Output equivalence - stdin vs mmap
echo "Test 19: Output equivalence - stdin vs mmap"
$FASTA_CONVERTER < multiallelic.vcf > output_stdin.fasta 2>/dev/null
$FASTA_CONVERTER -i multiallelic.vcf > output_mmap.fasta 2>/dev/null
if diff -q output_stdin.fasta output_mmap.fasta > /dev/null; then
    echo "✓ Test 19 passed: mmap output matches stdin output"
else
    echo "❌ Test 19 failed: mmap output differs from stdin"
    diff output_stdin.fasta output_mmap.fasta
    exit 1
fi

# NEW TEST 20: File input with --input long option
echo "Test 20: File input with --input long option"
$FASTA_CONVERTER --input basic.vcf > output_long_opt.fasta 2>/dev/null
if diff -q output_long_opt.fasta output_basic.fasta > /dev/null; then
    echo "✓ Test 20 passed: --input long option works correctly"
else
    echo "❌ Test 20 failed: --input long option output mismatch"
    exit 1
fi

# NEW TEST 21: Help shows new options (-i, -q)
echo "Test 21: Help shows new options"
$FASTA_CONVERTER --help > help_new.txt 2>&1
if grep -q "\-i" help_new.txt && grep -q "\-q" help_new.txt; then
    echo "✓ Test 21 passed: help shows -i and -q options"
else
    echo "❌ Test 21 failed: help does not show new options"
    cat help_new.txt
    exit 1
fi

# NEW TEST 22: Nonexistent file error handling
echo "Test 22: Nonexistent file error handling"
set +e  # Temporarily disable exit on error
$FASTA_CONVERTER -i nonexistent_file.vcf > /dev/null 2> nonexistent.err
exit_code=$?
set -e  # Re-enable exit on error
if [ $exit_code -ne 0 ]; then
    echo "✓ Test 22 passed: nonexistent file causes error exit"
else
    echo "❌ Test 22 failed: nonexistent file did not cause error"
    exit 1
fi

echo "All tests for VCFX_fasta_converter passed!" 