#!/usr/bin/env bash
set -e

echo "=== Testing VCFX_phase_checker ==="
mkdir -p out
mkdir -p data
mkdir -p expected

# Executable paths
VCFX_EXECUTABLE="../build/src/VCFX_phase_checker/VCFX_phase_checker"

# Check if executable exists
if [ ! -f "$VCFX_EXECUTABLE" ]; then
  echo "Error: $VCFX_EXECUTABLE not found!"
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
  local expect_warnings="${6:-false}"
  
  echo "Test $test_num: $test_desc"
  
  if [ "$expect_warnings" = "true" ]; then
    # Run with stderr redirected to a separate file
    eval "$test_cmd" > "$output_file" 2>out/warnings.log
    
    # Check if we got any warnings
    if [ ! -s out/warnings.log ]; then
      echo "  Error: Expected warnings but none were generated"
      exit 1
    fi
  else
    # Run without capturing stderr
    eval "$test_cmd" > "$output_file" 2>/dev/null
  fi
  
  # Compare output with expected
  diff -u "$expected_file" "$output_file" || {
    echo "  Test $test_num failed. Actual output:"
    cat "$output_file"
    exit 1
  }
  echo "  Test $test_num passed."
}

# Create test data if it doesn't exist
if [ ! -f data/phase_all_phased.vcf ]; then
  cat > data/phase_all_phased.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
chr1	10000	rs123	A	G	100	PASS	DP=55	GT:DP:GQ	0|1:30:99	1|0:25:95	0|1:28:98
chr1	20000	rs456	C	T	80	PASS	DP=60	GT:DP:GQ	0|0:28:99	0|0:32:98	0|0:27:97
chr1	30000	rs789	G	A	50	LowQual	DP=29	GT:DP:GQ	1|1:15:45	1|1:14:52	1|1:16:48
chr2	15000	rs101	T	C	90	PASS	DP=87	GT:DP:GQ	0|1:45:99	1|0:42:99	0|1:44:99
chr2	25000	rs102	G	T	30	LowQual	DP=18	GT:DP:GQ	1|1:10:35	1|1:8:30	1|1:9:32
EOF
fi

if [ ! -f data/phase_some_unphased.vcf ]; then
  cat > data/phase_some_unphased.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
chr1	10000	rs123	A	G	100	PASS	DP=55	GT:DP:GQ	0|1:30:99	1|0:25:95	0|1:28:98
chr1	20000	rs456	C	T	80	PASS	DP=60	GT:DP:GQ	0/0:28:99	0|0:32:98	0|0:27:97
chr1	30000	rs789	G	A	50	LowQual	DP=29	GT:DP:GQ	1|1:15:45	1|1:14:52	1|1:16:48
chr2	15000	rs101	T	C	90	PASS	DP=87	GT:DP:GQ	0|1:45:99	1/0:42:99	0|1:44:99
chr2	25000	rs102	G	T	30	LowQual	DP=18	GT:DP:GQ	1|1:10:35	1|1:8:30	1|1:9:32
EOF
fi

if [ ! -f data/phase_mixed_formats.vcf ]; then
  cat > data/phase_mixed_formats.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
chr1	10000	rs123	A	G	100	PASS	DP=55	GT:DP:GQ	0|1:30:99	1|0:25:95	0|1:28:98
chr1	20000	rs456	C	T	80	PASS	DP=60	GT:DP:GQ	0|0:28:99	0|0:32:98	0|0:27:97
chr1	30000	rs789	G	A	50	LowQual	DP=29	GT:DP:GQ	.|.:15:45	1|1:14:52	1|1:16:48
chr2	15000	rs101	T	C	90	PASS	DP=87	DP	.	.	.
chr2	25000	rs102	G	T	30	LowQual	DP=18	GT:DP:GQ	1|1:10:35	./.:8:30	1|1:9:32
EOF
fi

if [ ! -f data/phase_missing_gt.vcf ]; then
  cat > data/phase_missing_gt.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
chr1	10000	rs123	A	G	100	PASS	DP=55	DP:GQ	30:99	25:95	28:98
chr1	20000	rs456	C	T	80	PASS	DP=60	DP:GQ	28:99	32:98	27:97
chr1	30000	rs789	G	A	50	LowQual	DP=29	DP:GQ	15:45	14:52	16:48
EOF
fi

if [ ! -f data/phase_malformed.vcf ]; then
  cat > data/phase_malformed.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
malformed_line
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
too_few_columns
chr1	10000	rs123	A	G	100	PASS	DP=55	GT	X/Y
chr1	20000	rs456	C	T	80	PASS	DP=60	GT	0|X
EOF
fi

if [ ! -f data/phase_multi_allelic.vcf ]; then
  cat > data/phase_multi_allelic.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
chr1	10000	rs123	A	G,T	100	PASS	DP=55	GT:DP:GQ	0|1:30:99	1|2:25:95	0|1:28:98
chr1	20000	rs456	C	G,T,A	80	PASS	DP=60	GT:DP:GQ	0|2:28:99	0|2:32:98	0|3:27:97
chr1	30000	rs789	G	A,C	50	LowQual	DP=29	GT:DP:GQ	1|2:15:45	2|1:14:52	1|2:16:48
EOF
fi

if [ ! -f data/phase_multiple_ploidy.vcf ]; then
  cat > data/phase_multiple_ploidy.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
chr1	10000	rs123	A	G	100	PASS	DP=55	GT:DP:GQ	0:30:99	0|1:25:95	0|1:28:98
chr1	20000	rs456	C	T	80	PASS	DP=60	GT:DP:GQ	0|0:28:99	0|0:32:98	0|0|0:27:97
chr1	30000	rs789	G	A	50	LowQual	DP=29	GT:DP:GQ	1|1:15:45	1|1|1|1:14:52	1:16:48
EOF
fi

# Create large test data with 1000 variants for performance testing
if [ ! -f data/phase_large.vcf ]; then
  cat > data/phase_large.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
EOF

  # Add 1000 variant lines
  for i in $(seq 1 1000); do
    if [ $((i % 5)) -eq 0 ]; then
      # Every 5th line has unphased genotypes
      echo "chr1	$((10000 + $i))	rs$i	A	G	100	PASS	DP=50	GT:DP	0/1:30	0|1:25	0|1:28" >> data/phase_large.vcf
    else
      # All others are phased
      echo "chr1	$((10000 + $i))	rs$i	A	G	100	PASS	DP=50	GT:DP	0|1:30	0|1:25	0|1:28" >> data/phase_large.vcf
    fi
  done
fi

# Expected outputs
if [ ! -f expected/phase_all_phased.vcf ]; then
  cat > expected/phase_all_phased.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
chr1	10000	rs123	A	G	100	PASS	DP=55	GT:DP:GQ	0|1:30:99	1|0:25:95	0|1:28:98
chr1	20000	rs456	C	T	80	PASS	DP=60	GT:DP:GQ	0|0:28:99	0|0:32:98	0|0:27:97
chr1	30000	rs789	G	A	50	LowQual	DP=29	GT:DP:GQ	1|1:15:45	1|1:14:52	1|1:16:48
chr2	15000	rs101	T	C	90	PASS	DP=87	GT:DP:GQ	0|1:45:99	1|0:42:99	0|1:44:99
chr2	25000	rs102	G	T	30	LowQual	DP=18	GT:DP:GQ	1|1:10:35	1|1:8:30	1|1:9:32
EOF
fi

if [ ! -f expected/phase_some_unphased.vcf ]; then
  cat > expected/phase_some_unphased.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
chr1	10000	rs123	A	G	100	PASS	DP=55	GT:DP:GQ	0|1:30:99	1|0:25:95	0|1:28:98
chr1	30000	rs789	G	A	50	LowQual	DP=29	GT:DP:GQ	1|1:15:45	1|1:14:52	1|1:16:48
chr2	25000	rs102	G	T	30	LowQual	DP=18	GT:DP:GQ	1|1:10:35	1|1:8:30	1|1:9:32
EOF
fi

if [ ! -f expected/phase_mixed_formats.vcf ]; then
  cat > expected/phase_mixed_formats.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
chr1	10000	rs123	A	G	100	PASS	DP=55	GT:DP:GQ	0|1:30:99	1|0:25:95	0|1:28:98
chr1	20000	rs456	C	T	80	PASS	DP=60	GT:DP:GQ	0|0:28:99	0|0:32:98	0|0:27:97
EOF
fi

if [ ! -f expected/phase_missing_gt.vcf ]; then
  cat > expected/phase_missing_gt.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
EOF
fi

if [ ! -f expected/phase_malformed.vcf ]; then
  cat > expected/phase_malformed.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
malformed_line
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
EOF
fi

if [ ! -f expected/phase_multi_allelic.vcf ]; then
  cat > expected/phase_multi_allelic.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
chr1	10000	rs123	A	G,T	100	PASS	DP=55	GT:DP:GQ	0|1:30:99	1|2:25:95	0|1:28:98
chr1	20000	rs456	C	G,T,A	80	PASS	DP=60	GT:DP:GQ	0|2:28:99	0|2:32:98	0|3:27:97
chr1	30000	rs789	G	A,C	50	LowQual	DP=29	GT:DP:GQ	1|2:15:45	2|1:14:52	1|2:16:48
EOF
fi

if [ ! -f expected/phase_multiple_ploidy.vcf ]; then
  cat > expected/phase_multiple_ploidy.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
chr1	20000	rs456	C	T	80	PASS	DP=60	GT:DP:GQ	0|0:28:99	0|0:32:98	0|0|0:27:97
EOF
fi

# Test 1: All genotypes are phased
run_test 1 "All genotypes phased" \
  "$VCFX_EXECUTABLE < data/phase_all_phased.vcf" \
  "expected/phase_all_phased.vcf" \
  "out/phase_all_phased.vcf"

# Test 2: Some genotypes are unphased
run_test 2 "Some genotypes unphased" \
  "$VCFX_EXECUTABLE < data/phase_some_unphased.vcf" \
  "expected/phase_some_unphased.vcf" \
  "out/phase_some_unphased.vcf" \
  true

# Test 3: Mixed formats (missing GT in some lines)
run_test 3 "Mixed formats" \
  "$VCFX_EXECUTABLE < data/phase_mixed_formats.vcf" \
  "expected/phase_mixed_formats.vcf" \
  "out/phase_mixed_formats.vcf" \
  true

# Test 4: No GT field at all
run_test 4 "No GT field" \
  "$VCFX_EXECUTABLE < data/phase_missing_gt.vcf" \
  "expected/phase_missing_gt.vcf" \
  "out/phase_missing_gt.vcf" \
  true

# Test 5: Malformed VCF
run_test 5 "Malformed VCF" \
  "$VCFX_EXECUTABLE < data/phase_malformed.vcf" \
  "expected/phase_malformed.vcf" \
  "out/phase_malformed.vcf" \
  true

# Test 6: Multi-allelic variants with phased genotypes
run_test 6 "Multi-allelic variants" \
  "$VCFX_EXECUTABLE < data/phase_multi_allelic.vcf" \
  "expected/phase_multi_allelic.vcf" \
  "out/phase_multi_allelic.vcf"

# Test 7: Mixed ploidy
run_test 7 "Mixed ploidy" \
  "$VCFX_EXECUTABLE < data/phase_multiple_ploidy.vcf" \
  "expected/phase_multiple_ploidy.vcf" \
  "out/phase_multiple_ploidy.vcf" \
  true

# Test 8: Performance with large file
echo "Test 8: Performance with large file"
time $VCFX_EXECUTABLE < data/phase_large.vcf > out/phase_large.vcf 2>out/phase_large_warnings.log

# Count how many lines were output vs input
input_count=$(grep -v "^#" data/phase_large.vcf | wc -l)
output_count=$(grep -v "^#" out/phase_large.vcf | wc -l)
warning_count=$(wc -l < out/phase_large_warnings.log)

# We expect 20% to be filtered out (every 5th line)
expected_output_count=$((input_count * 4 / 5))
expected_warning_count=$((input_count / 5))

echo "  Input variants: $input_count"
echo "  Output variants: $output_count"
echo "  Warnings generated: $warning_count"

# Allow for some tolerance in the counts
if [ "$output_count" -lt "$((expected_output_count - 10))" ] || [ "$output_count" -gt "$((expected_output_count + 10))" ]; then
  echo "  Test 8 failed. Expected around $expected_output_count variants in output, but found $output_count."
  exit 1
fi

if [ "$warning_count" -lt "$((expected_warning_count - 10))" ] || [ "$warning_count" -gt "$((expected_warning_count + 10))" ]; then
  echo "  Test 8 failed. Expected around $expected_warning_count warnings, but found $warning_count."
  exit 1
fi

echo "  Test 8 passed."

# Test 9: Help display
echo "Test 9: Verifying help display"
$VCFX_EXECUTABLE --help > out/phase_help.txt
if ! grep -q "VCFX_phase_checker" out/phase_help.txt; then
  echo "  Test 9 failed. Help message not displayed correctly."
  cat out/phase_help.txt
  exit 1
fi
echo "  Test 9 passed."

# Test 10: File input mode with -i flag
echo "Test 10: File input mode with -i flag"
$VCFX_EXECUTABLE -i data/phase_all_phased.vcf > out/phase_mmap_test.vcf 2>/dev/null
diff -u expected/phase_all_phased.vcf out/phase_mmap_test.vcf || {
    echo "  Test 10 failed. mmap output differs from expected."
    exit 1
}
echo "  Test 10 passed."

# Test 11: Positional file argument
echo "Test 11: Positional file argument"
$VCFX_EXECUTABLE data/phase_all_phased.vcf > out/phase_positional_test.vcf 2>/dev/null
diff -u expected/phase_all_phased.vcf out/phase_positional_test.vcf || {
    echo "  Test 11 failed. Positional file argument output differs from expected."
    exit 1
}
echo "  Test 11 passed."

# Test 12: Quiet mode with -q flag
echo "Test 12: Quiet mode with -q flag"
$VCFX_EXECUTABLE -q < data/phase_some_unphased.vcf > out/phase_quiet_test.vcf 2>out/phase_quiet_warnings.log
diff -u expected/phase_some_unphased.vcf out/phase_quiet_test.vcf || {
    echo "  Test 12 failed. Quiet mode output differs from expected."
    exit 1
}
# Check that no warnings were emitted
if [ -s out/phase_quiet_warnings.log ]; then
    echo "  Test 12 failed. Quiet mode should not emit warnings."
    cat out/phase_quiet_warnings.log
    exit 1
fi
echo "  Test 12 passed."

# Test 13: Output equivalence (stdin vs mmap)
echo "Test 13: Output equivalence (stdin vs mmap)"
$VCFX_EXECUTABLE < data/phase_large.vcf > out/phase_stdin_test.vcf 2>/dev/null
$VCFX_EXECUTABLE -i data/phase_large.vcf > out/phase_mmap_large_test.vcf 2>/dev/null
diff -u out/phase_stdin_test.vcf out/phase_mmap_large_test.vcf || {
    echo "  Test 13 failed. stdin and mmap outputs differ."
    exit 1
}
echo "  Test 13 passed."

# Test 14: Performance comparison with 1000 variants
echo "Test 14: Performance comparison with 1000 variants"
echo "  stdin mode:"
time $VCFX_EXECUTABLE < data/phase_large.vcf > /dev/null 2>&1
echo "  mmap mode:"
time $VCFX_EXECUTABLE -i data/phase_large.vcf > /dev/null 2>&1
echo "  Test 14 passed."

# Test 15: Large-scale performance test (100K variants) - optional benchmark
echo "Test 15: Large-scale performance test (100K variants)"
# Create large test file if it doesn't exist
if [ ! -f data/phase_100k.vcf ]; then
  echo "  Creating 100K variant test file..."
  cat > data/phase_100k.vcf << EOF
##fileformat=VCFv4.2
##source=VCFX
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2	SAMPLE3
EOF

  # Add 100000 variant lines
  for i in $(seq 1 100000); do
    if [ $((i % 5)) -eq 0 ]; then
      # Every 5th line has unphased genotypes
      echo "chr1	$((10000 + $i))	rs$i	A	G	100	PASS	DP=50	GT:DP	0/1:30	0|1:25	0|1:28" >> data/phase_100k.vcf
    else
      # All others are phased
      echo "chr1	$((10000 + $i))	rs$i	A	G	100	PASS	DP=50	GT:DP	0|1:30	0|1:25	0|1:28" >> data/phase_100k.vcf
    fi
  done
  echo "  100K variant file created."
fi

echo "  Running 100K variant performance test..."
echo "  stdin mode:"
time $VCFX_EXECUTABLE -q < data/phase_100k.vcf > /dev/null 2>&1
echo "  mmap mode:"
time $VCFX_EXECUTABLE -q -i data/phase_100k.vcf > /dev/null 2>&1

# Verify output counts match between modes
stdin_count=$($VCFX_EXECUTABLE -q < data/phase_100k.vcf 2>/dev/null | grep -v "^#" | wc -l)
mmap_count=$($VCFX_EXECUTABLE -q -i data/phase_100k.vcf 2>/dev/null | grep -v "^#" | wc -l)
if [ "$stdin_count" -ne "$mmap_count" ]; then
    echo "  Test 15 failed. stdin output count ($stdin_count) differs from mmap count ($mmap_count)."
    exit 1
fi
echo "  Output counts match: $stdin_count variants"
echo "  Test 15 passed."

echo "All VCFX_phase_checker tests passed!" 