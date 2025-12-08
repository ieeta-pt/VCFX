#!/usr/bin/env bash

echo "=== Testing VCFX_allele_counter ==="

# Executable paths
VCFX_EXECUTABLE="../build/src/VCFX_allele_counter/VCFX_allele_counter"

# Check if executable exists
if [ ! -f "$VCFX_EXECUTABLE" ]; then
  echo "Error: $VCFX_EXECUTABLE not found!"
  echo "Make sure you've built the project before running tests."
  exit 1
fi

# Create output directories if they don't exist
mkdir -p out
mkdir -p data
mkdir -p expected

PASSED=0
FAILED=0

run_test() {
  local name="$1"
  local result="$2"
  if [ "$result" -eq 0 ]; then
    echo "Test $name: PASSED"
    ((PASSED++))
  else
    echo "Test $name: FAILED"
    ((FAILED++))
  fi
}

# Test 1: Scenario A - stdin mode
echo "Test 1: Scenario A (stdin mode)"
$VCFX_EXECUTABLE -q < data/allele_counter_A.vcf > out/allele_counter_A_out.tsv 2>/dev/null
diff -q expected/allele_counter_A_out.tsv out/allele_counter_A_out.tsv > /dev/null 2>&1
run_test "1" $?

# Test 2: Scenario B - specific sample
echo "Test 2: Scenario B (specific sample)"
$VCFX_EXECUTABLE -q --samples "Y" < data/allele_counter_B.vcf > out/allele_counter_B_out.tsv 2>/dev/null
diff -q expected/allele_counter_B_out.tsv out/allele_counter_B_out.tsv > /dev/null 2>&1
run_test "2" $?

# Test 3: mmap mode with -i flag
echo "Test 3: mmap mode with -i flag"
$VCFX_EXECUTABLE -q -i data/allele_counter_A.vcf > out/allele_counter_A_mmap.tsv 2>/dev/null
diff -q expected/allele_counter_A_out.tsv out/allele_counter_A_mmap.tsv > /dev/null 2>&1
run_test "3" $?

# Test 4: mmap mode with positional argument
echo "Test 4: mmap mode with positional argument"
$VCFX_EXECUTABLE -q data/allele_counter_A.vcf > out/allele_counter_A_pos.tsv 2>/dev/null
diff -q expected/allele_counter_A_out.tsv out/allele_counter_A_pos.tsv > /dev/null 2>&1
run_test "4" $?

# Test 5: mmap mode with specific sample
echo "Test 5: mmap mode with specific sample"
$VCFX_EXECUTABLE -q -s "Y" -i data/allele_counter_B.vcf > out/allele_counter_B_mmap.tsv 2>/dev/null
diff -q expected/allele_counter_B_out.tsv out/allele_counter_B_mmap.tsv > /dev/null 2>&1
run_test "5" $?

# Test 6: stdin vs mmap output equivalence
echo "Test 6: stdin vs mmap output equivalence"
$VCFX_EXECUTABLE -q < data/allele_counter_A.vcf > out/allele_counter_stdin.tsv 2>/dev/null
$VCFX_EXECUTABLE -q -i data/allele_counter_A.vcf > out/allele_counter_mmap.tsv 2>/dev/null
diff -q out/allele_counter_stdin.tsv out/allele_counter_mmap.tsv > /dev/null 2>&1
run_test "6" $?

# Test 7: help flag
echo "Test 7: help flag"
$VCFX_EXECUTABLE --help 2>&1 | grep -q "Usage:" > /dev/null
run_test "7" $?

# Test 8: version flag
echo "Test 8: version flag"
$VCFX_EXECUTABLE --version 2>&1 | grep -q "2.0" > /dev/null
run_test "8" $?

# Test 9: missing file error
echo "Test 9: missing file error handling"
$VCFX_EXECUTABLE -q -i nonexistent_file.vcf > /dev/null 2>&1
if [ $? -ne 0 ]; then
  run_test "9" 0
else
  run_test "9" 1
fi

# Test 10: Create and test multi-allelic file
echo "Test 10: multi-allelic variants"
cat > out/test_multi.vcf << 'EOF'
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	.	A	T,C	.	.	.	GT	0/1	1/2
chr1	200	.	G	A	.	.	.	GT	0/0	1|1
EOF
$VCFX_EXECUTABLE -q -i out/test_multi.vcf > out/test_multi_out.tsv 2>/dev/null
# Verify output has correct counts
# S1: 0/1 -> ref=1, alt=1
# S2: 1/2 -> ref=0, alt=2 (both 1 and 2 are ALT)
grep -q "S1	1	1" out/test_multi_out.tsv && grep -q "S2	0	2" out/test_multi_out.tsv
run_test "10" $?

# Test 11: phased genotypes
echo "Test 11: phased genotypes"
cat > out/test_phased.vcf << 'EOF'
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	100	.	A	T	.	.	.	GT	0|1
chr1	200	.	G	C	.	.	.	GT	1|0
EOF
$VCFX_EXECUTABLE -q -i out/test_phased.vcf > out/test_phased_out.tsv 2>/dev/null
# Both should have ref=1, alt=1
grep -c "S1	1	1" out/test_phased_out.tsv | grep -q "2"
run_test "11" $?

# Test 12: missing genotypes
echo "Test 12: missing genotypes"
cat > out/test_missing.vcf << 'EOF'
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	100	.	A	T	.	.	.	GT	./.
chr1	200	.	G	C	.	.	.	GT	0/.
EOF
$VCFX_EXECUTABLE -q -i out/test_missing.vcf > out/test_missing_out.tsv 2>/dev/null
# First should be 0,0; second should be 1,0
grep "chr1	100" out/test_missing_out.tsv | grep -q "0	0" && \
grep "chr1	200" out/test_missing_out.tsv | grep -q "1	0"
run_test "12" $?

# Test 13: quiet mode (no info messages)
echo "Test 13: quiet mode"
$VCFX_EXECUTABLE -q -i data/allele_counter_A.vcf 2>&1 | grep -q "Info:" > /dev/null
if [ $? -eq 0 ]; then
  run_test "13" 1  # Failed - should not have Info messages
else
  run_test "13" 0  # Passed - no Info messages
fi

# Test 14: non-quiet mode (has info messages)
echo "Test 14: non-quiet mode"
$VCFX_EXECUTABLE -i data/allele_counter_A.vcf > /dev/null 2>&1
# Just verify it runs successfully
run_test "14" $?

# Summary
echo ""
echo "=== Test Summary ==="
echo "Passed: $PASSED"
echo "Failed: $FAILED"

if [ "$FAILED" -gt 0 ]; then
  echo "Some tests failed!"
  exit 1
fi

echo "All allele_counter tests passed!"
