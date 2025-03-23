#!/usr/bin/env bash

# Track failures
failures=0

echo "=== Testing VCFX_sorter ==="
mkdir -p out
mkdir -p data
mkdir -p expected

# Define paths
SORTER="../build/src/VCFX_sorter/VCFX_sorter"

# Check if executable exists
if [ ! -f "$SORTER" ]; then
  echo "Error: $SORTER not found!"
  echo "Make sure you've built the project before running tests."
  exit 1
fi

# Function to run a test and report results
run_test() {
  local test_num="$1"
  local test_desc="$2"
  local test_cmd="$3"
  local expected_file="$4"
  local output_file="$5"
  local expect_failure="${6:-false}"
  
  echo "Test $test_num: $test_desc"
  
  if [ "$expect_failure" = "true" ]; then
    if eval "$test_cmd" > "$output_file" 2>&1; then
      echo "  Error: Expected failure but got success"
      ((failures++))
    else
      echo "  Test $test_num passed (expected failure)"
    fi
  else
    eval "$test_cmd" > "$output_file" 2> /dev/null
    if ! diff -u "$expected_file" "$output_file"; then
      echo "  Test $test_num failed. Actual output:"
      cat "$output_file"
      ((failures++))
    else
      echo "  Test $test_num passed."
    fi
  fi
}

# Create test data: Unsorted VCF with various chromosome types
if [ ! -f data/unsorted.vcf ]; then
  cat > data/unsorted.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr10	1000	.	A	G	100	PASS	DP=50;AF=0.5
chr2	2000	.	C	T	90	PASS	DP=40;AF=0.4
chr1	3000	.	G	A	80	PASS	DP=30;AF=0.3
chr11	1500	.	T	C	70	PASS	DP=25;AF=0.25
chrX	500	.	A	G	60	PASS	DP=20;AF=0.2
chr1	1000	.	C	G	50	PASS	DP=15;AF=0.15
chrMT	200	.	T	A	40	PASS	DP=10;AF=0.1
chr10	500	.	G	C	30	PASS	DP=5;AF=0.05
chr10_alt	700	.	A	T	20	PASS	DP=8;AF=0.08
chr2_random	100	.	C	A	10	PASS	DP=3;AF=0.03
EOF
fi

# Create expected output for lexicographic sort
if [ ! -f expected/lex_sorted.vcf ]; then
  cat > expected/lex_sorted.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	1000	.	C	G	50	PASS	DP=15;AF=0.15
chr1	3000	.	G	A	80	PASS	DP=30;AF=0.3
chr10	500	.	G	C	30	PASS	DP=5;AF=0.05
chr10	1000	.	A	G	100	PASS	DP=50;AF=0.5
chr10_alt	700	.	A	T	20	PASS	DP=8;AF=0.08
chr11	1500	.	T	C	70	PASS	DP=25;AF=0.25
chr2	2000	.	C	T	90	PASS	DP=40;AF=0.4
chr2_random	100	.	C	A	10	PASS	DP=3;AF=0.03
chrMT	200	.	T	A	40	PASS	DP=10;AF=0.1
chrX	500	.	A	G	60	PASS	DP=20;AF=0.2
EOF
fi

# Create expected output for natural sort
if [ ! -f expected/natural_sorted.vcf ]; then
  cat > expected/natural_sorted.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	1000	.	C	G	50	PASS	DP=15;AF=0.15
chr1	3000	.	G	A	80	PASS	DP=30;AF=0.3
chr2	2000	.	C	T	90	PASS	DP=40;AF=0.4
chr2_random	100	.	C	A	10	PASS	DP=3;AF=0.03
chr10	500	.	G	C	30	PASS	DP=5;AF=0.05
chr10	1000	.	A	G	100	PASS	DP=50;AF=0.5
chr10_alt	700	.	A	T	20	PASS	DP=8;AF=0.08
chr11	1500	.	T	C	70	PASS	DP=25;AF=0.25
chrMT	200	.	T	A	40	PASS	DP=10;AF=0.1
chrX	500	.	A	G	60	PASS	DP=20;AF=0.2
EOF
fi

# Create test data with mixed chromosome prefixes
if [ ! -f data/mixed_prefixes.vcf ]; then
  cat > data/mixed_prefixes.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr10	1000	.	A	G	100	PASS	DP=50
10	1500	.	C	T	90	PASS	DP=45
Chr2	2000	.	G	A	80	PASS	DP=40
2	2500	.	T	C	70	PASS	DP=35
CHR1	3000	.	A	G	60	PASS	DP=30
1	3500	.	C	T	50	PASS	DP=25
EOF
fi

# Create expected output for mixed prefixes with natural sort
if [ ! -f expected/mixed_natural_sorted.vcf ]; then
  cat > expected/mixed_natural_sorted.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	3500	.	C	T	50	PASS	DP=25
CHR1	3000	.	A	G	60	PASS	DP=30
2	2500	.	T	C	70	PASS	DP=35
Chr2	2000	.	G	A	80	PASS	DP=40
10	1500	.	C	T	90	PASS	DP=45
chr10	1000	.	A	G	100	PASS	DP=50
EOF
fi

# Create data with invalid POS values
if [ ! -f data/invalid_pos.vcf ]; then
  cat > data/invalid_pos.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	1000	.	A	G	100	PASS	DP=50
chr2	invalid	.	C	T	90	PASS	DP=45
chr3	3000	.	G	A	80	PASS	DP=40
EOF
fi

# Run the tests
# Test 1: Lexicographic sort (default)
run_test 1 "Lexicographic chromosome sorting" \
  "$SORTER < data/unsorted.vcf" \
  "expected/lex_sorted.vcf" \
  "out/lex_sorted.vcf"

# Test 2: Natural chromosome sort
run_test 2 "Natural chromosome sorting" \
  "$SORTER --natural-chr < data/unsorted.vcf" \
  "expected/natural_sorted.vcf" \
  "out/natural_sorted.vcf"

# Test 3: Natural sort with mixed prefixes (chr/Chr/CHR/none)
run_test 3 "Natural sort with mixed chromosome prefixes" \
  "$SORTER --natural-chr < data/mixed_prefixes.vcf" \
  "expected/mixed_natural_sorted.vcf" \
  "out/mixed_natural_sorted.vcf"

# Test 4: Help message
echo "Test 4: Help message"
$SORTER --help > out/sorter_help.txt
if ! grep -q "VCFX_sorter" out/sorter_help.txt || ! grep -q "natural-chr" out/sorter_help.txt; then
  echo "  Test 4 failed: Help message doesn't contain expected content"
  ((failures++))
else
  echo "  Test 4 passed."
fi

# Test 5: Invalid inputs (expect warnings but not crashes)
echo "Test 5: Handling invalid inputs"
$SORTER < data/invalid_pos.vcf > out/invalid_pos_out.vcf 2> out/invalid_pos_warnings.txt
if [ $? -ne 0 ]; then
  echo "  Test 5 failed: Program crashed on invalid input"
  ((failures++))
else
  # The program should output some valid content and not crash
  if [ ! -s out/invalid_pos_out.vcf ]; then
    echo "  Test 5 failed: No output produced for invalid input"
    ((failures++))
  else
    # Check if the output contains the header and the expected record
    if ! grep -q "#CHROM" out/invalid_pos_out.vcf || ! grep -q "chr1" out/invalid_pos_out.vcf; then
      echo "  Test 5 failed: Output missing expected records"
      ((failures++)) 
    else
      echo "  Test 5 passed: Program handled invalid input gracefully"
    fi
  fi
fi

# Report results
if [ $failures -eq 0 ]; then
  echo "All VCFX_sorter tests passed!"
  exit 0
else
  echo "$failures VCFX_sorter tests failed."
  exit 1
fi 