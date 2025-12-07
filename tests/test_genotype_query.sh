#!/usr/bin/env bash
set -e

###############################################################################
# COMPLETE TEST SCRIPT FOR VCFX_genotype_query
#
# This version IGNORES trailing whitespace / blank lines in a way
# that works on both macOS (BSD sed) and GNU/Linux sed.
#
# Steps:
#   1) Creates temporary test data (VCFs) and expected outputs.
#   2) Runs the tool against each test input.
#   3) Compares the actual output to the expected file with diff,
#      but after removing trailing whitespace/blank lines from each side.
#   4) Tests help message and missing arguments usage, etc.
#
# Just run:
#   ./test_genotype_query.sh
###############################################################################

# Path to the compiled VCFX_genotype_query binary (adjust if needed):
TOOL="../build/src/VCFX_genotype_query/VCFX_genotype_query"

# Directories for test data, expected outputs, and actual output:
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
TMP_DATA_DIR="${SCRIPT_DIR}/data/genotype_query"
TMP_EXP_DIR="${SCRIPT_DIR}/expected/genotype_query"
TMP_OUT_DIR="${SCRIPT_DIR}/tmp/genotype_query"

mkdir -p "$TMP_DATA_DIR" "$TMP_EXP_DIR" "$TMP_OUT_DIR"

# Ensure the tool is built and executable
if [ ! -x "$TOOL" ]; then
  echo "Error: $TOOL not found or not executable. Please build your tool first."
  exit 1
fi

###############################################################################
# Helper function: strip_trailing
#  - Removes trailing whitespace from each line (via sed).
#  - Removes trailing blank lines at EOF (via awk).
# Usage: strip_trailing file
###############################################################################
strip_trailing() {
  local file="$1"
  # 1) Remove trailing whitespace with sed
  # 2) Remove trailing blank lines with an awk approach:
  #    - We'll store lines in memory
  #    - Then only print up through the last non-blank line
  sed 's/[[:space:]]*$//' "$file" | awk '
    { lines[NR] = $0 }
    NF { lastNonEmpty = NR }
    END {
      for (i=1; i<=lastNonEmpty; i++) {
        print lines[i]
      }
    }
  '
}

###############################################################################
# Helper: Create an input VCF file (with embedded data).
###############################################################################
create_input() {
  local file_name="$1"
  local vcf_data="$2"
  # We use 'echo -e' so that \t is interpreted as actual tabs:
  echo -e "$vcf_data" > "${TMP_DATA_DIR}/${file_name}"
}

###############################################################################
# Helper: Create an expected output VCF file.
###############################################################################
create_expected() {
  local file_name="$1"
  local vcf_data="$2"
  echo -e "$vcf_data" > "${TMP_EXP_DIR}/${file_name}"
}

###############################################################################
# run_test: runs the tool and compares against expected,
#           ignoring trailing whitespace/blank lines.
# Usage: run_test <test_name> <query> <input_file> <expected_file> [--strict]
###############################################################################
run_test() {
  local test_name="$1"
  local query="$2"
  local input_file="$3"
  local expected_file="$4"
  local strict_flag="$5"

  echo "-----------------------------------------------------"
  echo "Running test: $test_name"
  echo "  Query:    $query"
  echo "  Input:    $input_file"
  echo "  Expected: $expected_file"
  echo "  Strict:   $strict_flag"
  echo

  local output_path="${TMP_OUT_DIR}/${test_name}_output.vcf"

  # Build command
  local cmd="$TOOL --genotype-query \"$query\""
  if [ "$strict_flag" = "--strict" ]; then
    cmd="$cmd --strict"
  fi

  # Run
  eval "$cmd < \"${TMP_DATA_DIR}/${input_file}\" > \"${output_path}\""

  # Compare using strip_trailing
  if diff -u <(strip_trailing "${TMP_EXP_DIR}/${expected_file}") \
             <(strip_trailing "${output_path}"); then
    echo "✅ Test passed: $test_name"
  else
    echo "❌ Test failed: $test_name"
    echo "See diff above."
    exit 1
  fi
  echo
}

###############################################################################
# 1) CREATE INPUT VCF FILES (WITH REAL TABS)
###############################################################################

# SINGLE-SAMPLE, 3 lines
create_input "single_sample.vcf" "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tONLYSAMPLE
1\t100\trsA\tA\tG\t50\tPASS\t.\tGT\t0/1
1\t200\trsB\tA\tG\t50\tPASS\t.\tGT\t0|1
1\t300\trsC\tA\tG\t50\tPASS\t.\tGT\t1|1
"

# MULTI-SAMPLE, includes multi-allelic
create_input "multi_sample.vcf" "\
##fileformat=VCFv4.2
##contig=<ID=1>
##contig=<ID=2>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t100\trsX\tA\tG\t.\tPASS\t.\tGT\t0/0\t0|1\t1/1
1\t200\trsY\tA\tG,T\t.\tPASS\t.\tGT\t1/2\t2/2\t0/2
2\t300\trsZ\tC\tT\t.\tPASS\t.\tGT\t1|1\t1/1\t1/0
2\t400\t.\tG\tA\t.\tPASS\t.\tGT\t.\t.\t0/1
"

# MISSING/MALFORMED lines
create_input "missing_malformed.vcf" "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2
1\t150\trsM\tC\tG\t.\tPASS\t.\tGT:DP\t0/1:12\t1/1:30
1\t200\trsN\tA\tT\t.\tPASS\t.\tGT\t0/1\t.
1\t250\trsO\tA\tG\t.\tPASS\t.\t\t1/1\t1/1
chr1\t300  # <10 fields on purpose
1\t400\trsQ\tG\tA\t99\tPASS\t.\tDP\t10\t15
"

###############################################################################
# 2) CREATE EXPECTED OUTPUT FILES
###############################################################################

# 2A) SINGLE-SAMPLE, flexible query=0/1 => matches both (0/1 and 0|1)
create_expected "single_sample_flex_01.vcf" "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tONLYSAMPLE
1\t100\trsA\tA\tG\t50\tPASS\t.\tGT\t0/1
1\t200\trsB\tA\tG\t50\tPASS\t.\tGT\t0|1
"

# 2B) SINGLE-SAMPLE, strict query=0|1 => only matches the line with 0|1, not 0/1
create_expected "single_sample_strict_01.vcf" "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tONLYSAMPLE
1\t200\trsB\tA\tG\t50\tPASS\t.\tGT\t0|1
"

# 2C) MULTI-SAMPLE, query=1/1 flexible => lines that have any sample = 1/1 or 1|1
create_expected "multi_11_flexible.vcf" "\
##fileformat=VCFv4.2
##contig=<ID=1>
##contig=<ID=2>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t100\trsX\tA\tG\t.\tPASS\t.\tGT\t0/0\t0|1\t1/1
2\t300\trsZ\tC\tT\t.\tPASS\t.\tGT\t1|1\t1/1\t1/0
"

# 2D) MULTI-SAMPLE, query=1|1 strict => only lines where sample exactly has '1|1'
create_expected "multi_11_strict.vcf" "\
##fileformat=VCFv4.2
##contig=<ID=1>
##contig=<ID=2>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
2\t300\trsZ\tC\tT\t.\tPASS\t.\tGT\t1|1\t1/1\t1/0
"

# 2E) Multi-allelic, query=1/2 flexible => line 1:200 => S1=1/2 => match
create_expected "multi_12_flexible.vcf" "\
##fileformat=VCFv4.2
##contig=<ID=1>
##contig=<ID=2>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t200\trsY\tA\tG,T\t.\tPASS\t.\tGT\t1/2\t2/2\t0/2
"

# 2F) Missing/malformed, query=0/1 flexible
create_expected "missing_malformed_01.vcf" "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2
1\t150\trsM\tC\tG\t.\tPASS\t.\tGT:DP\t0/1:12\t1/1:30
1\t200\trsN\tA\tT\t.\tPASS\t.\tGT\t0/1\t.
"

# 2G) No match => only header remains
create_expected "no_match.vcf" "\
##fileformat=VCFv4.2
##contig=<ID=1>
##contig=<ID=2>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
"

###############################################################################
# 3) RUN THE TESTS
###############################################################################

# Test 1: single-sample, flexible 0/1
run_test "test_1_single_flex" \
         "0/1" \
         "single_sample.vcf" \
         "single_sample_flex_01.vcf"

# Test 2: single-sample, strict 0|1
run_test "test_2_single_strict" \
         "0|1" \
         "single_sample.vcf" \
         "single_sample_strict_01.vcf" \
         "--strict"

# Test 3: multi-sample, query=1/1 flexible
run_test "test_3_multi_11_flex" \
         "1/1" \
         "multi_sample.vcf" \
         "multi_11_flexible.vcf"

# Test 4: multi-sample, query=1|1 strict
run_test "test_4_multi_11_strict" \
         "1|1" \
         "multi_sample.vcf" \
         "multi_11_strict.vcf" \
         "--strict"

# Test 5: multi-allelic, query=1/2 flexible
run_test "test_5_multi_12_flex" \
         "1/2" \
         "multi_sample.vcf" \
         "multi_12_flexible.vcf"

# Test 6: missing/malformed, query=0/1 flexible
run_test "test_6_missing_malformed" \
         "0/1" \
         "missing_malformed.vcf" \
         "missing_malformed_01.vcf"

# Test 7: no match => only header
run_test "test_7_no_match" \
         "2/3" \
         "multi_sample.vcf" \
         "no_match.vcf"

###############################################################################
# Test 8: help message
###############################################################################
echo "-----------------------------------------------------"
echo "Running test_8_help_message"
HELP_OUT="${TMP_OUT_DIR}/help_message.txt"
"$TOOL" --help > "$HELP_OUT" 2>&1 || true
if grep -q "VCFX_genotype_query" "$HELP_OUT" && grep -q "Usage:" "$HELP_OUT"; then
    echo "✅ test_8_help_message passed"
else
    echo "❌ test_8_help_message failed"
    echo "Got:"
    cat "$HELP_OUT"
    exit 1
fi
echo

###############################################################################
# Test 9: missing arguments => usage error
###############################################################################
echo "-----------------------------------------------------"
echo "Running test_9_missing_args"
MISS_OUT="${TMP_OUT_DIR}/missing_args.txt"
if "$TOOL" 2> "$MISS_OUT"; then
    echo "❌ test_9_missing_args: the tool did not fail on missing args!"
    exit 1
else
    if grep -q "Usage:" "$MISS_OUT"; then
        echo "✅ test_9_missing_args passed"
    else
        echo "❌ test_9_missing_args failed - no usage message found"
        echo "Got:"
        cat "$MISS_OUT"
        exit 1
    fi
fi
echo

###############################################################################
# Test 10: Long option with equals
###############################################################################
echo "-----------------------------------------------------"
echo "Running test_10_long_option_equals"
EQUALS_OUT="${TMP_OUT_DIR}/long_equals_output.vcf"
"$TOOL" --genotype-query="0/1" < "${TMP_DATA_DIR}/single_sample.vcf" > "$EQUALS_OUT"

if diff -u <(strip_trailing "${TMP_EXP_DIR}/single_sample_flex_01.vcf") \
           <(strip_trailing "${EQUALS_OUT}"); then
  echo "✅ test_10_long_option_equals passed"
else
  echo "❌ test_10_long_option_equals failed"
  exit 1
fi
echo

###############################################################################
# Test 11: File input mode with -i flag
###############################################################################
echo "-----------------------------------------------------"
echo "Running test_11_file_input_mode"
FILE_OUT="${TMP_OUT_DIR}/test_11_file_input.vcf"
"$TOOL" -g "0/1" -i "${TMP_DATA_DIR}/single_sample.vcf" > "$FILE_OUT"

if diff -u <(strip_trailing "${TMP_EXP_DIR}/single_sample_flex_01.vcf") \
           <(strip_trailing "$FILE_OUT"); then
    echo "✅ test_11_file_input_mode passed"
else
    echo "❌ test_11_file_input_mode failed"
    exit 1
fi
echo

###############################################################################
# Test 12: Positional file argument
###############################################################################
echo "-----------------------------------------------------"
echo "Running test_12_positional_file"
POS_OUT="${TMP_OUT_DIR}/test_12_positional.vcf"
"$TOOL" -g "0/1" "${TMP_DATA_DIR}/single_sample.vcf" > "$POS_OUT"

if diff -u <(strip_trailing "${TMP_EXP_DIR}/single_sample_flex_01.vcf") \
           <(strip_trailing "$POS_OUT"); then
    echo "✅ test_12_positional_file passed"
else
    echo "❌ test_12_positional_file failed"
    exit 1
fi
echo

###############################################################################
# Test 13: Quiet mode
###############################################################################
echo "-----------------------------------------------------"
echo "Running test_13_quiet_mode"
QUIET_ERR="${TMP_OUT_DIR}/test_13_stderr.txt"
"$TOOL" -g "0/1" -q < "${TMP_DATA_DIR}/missing_malformed.vcf" > /dev/null 2> "$QUIET_ERR"

if [ ! -s "$QUIET_ERR" ]; then
    echo "✅ test_13_quiet_mode passed (no stderr output)"
else
    echo "❌ test_13_quiet_mode failed (stderr should be empty)"
    cat "$QUIET_ERR"
    exit 1
fi
echo

###############################################################################
# Test 14: Output equivalence (stdin vs mmap)
###############################################################################
echo "-----------------------------------------------------"
echo "Running test_14_output_equivalence"
STDIN_OUT="${TMP_OUT_DIR}/test_14_stdin.vcf"
MMAP_OUT="${TMP_OUT_DIR}/test_14_mmap.vcf"

"$TOOL" -g "1/1" < "${TMP_DATA_DIR}/multi_sample.vcf" > "$STDIN_OUT"
"$TOOL" -g "1/1" -i "${TMP_DATA_DIR}/multi_sample.vcf" > "$MMAP_OUT"

if diff -u <(strip_trailing "$STDIN_OUT") <(strip_trailing "$MMAP_OUT"); then
    echo "✅ test_14_output_equivalence passed"
else
    echo "❌ test_14_output_equivalence failed"
    exit 1
fi
echo

###############################################################################
# Test 15: Performance test (1000 variants)
###############################################################################
echo "-----------------------------------------------------"
echo "Running test_15_performance_1000"
PERF_FILE="${TMP_DATA_DIR}/perf_1000.vcf"

# Generate test file
{
    echo -e "##fileformat=VCFv4.2"
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3"
    for i in $(seq 1 1000); do
        gt=$((i % 3))
        echo -e "chr1\t${i}00\trs${i}\tA\tG\t.\tPASS\t.\tGT\t0/${gt}\t${gt}|1\t1/1"
    done
} > "$PERF_FILE"

# Time stdin mode (run 3 times and take average)
STDIN_TOTAL=0
for run in 1 2 3; do
    START=$(python3 -c "import time; print(time.time())")
    "$TOOL" -g "0/1" < "$PERF_FILE" > /dev/null 2>&1
    END=$(python3 -c "import time; print(time.time())")
    ELAPSED=$(python3 -c "print(($END - $START)*1000)")
    STDIN_TOTAL=$(python3 -c "print($STDIN_TOTAL + $ELAPSED)")
done
STDIN_TIME=$(python3 -c "print(f'{$STDIN_TOTAL/3:.1f}')")

# Time mmap mode (run 3 times and take average)
MMAP_TOTAL=0
for run in 1 2 3; do
    START=$(python3 -c "import time; print(time.time())")
    "$TOOL" -g "0/1" -i "$PERF_FILE" > /dev/null 2>&1
    END=$(python3 -c "import time; print(time.time())")
    ELAPSED=$(python3 -c "print(($END - $START)*1000)")
    MMAP_TOTAL=$(python3 -c "print($MMAP_TOTAL + $ELAPSED)")
done
MMAP_TIME=$(python3 -c "print(f'{$MMAP_TOTAL/3:.1f}')")

echo "Performance (1000 variants, 3 samples, avg of 3 runs):"
echo "  stdin mode: ${STDIN_TIME}ms"
echo "  mmap mode:  ${MMAP_TIME}ms"
echo "✅ test_15_performance_1000 completed"
echo

###############################################################################
# Test 16: Performance test (100K variants - benchmark only)
###############################################################################
echo "-----------------------------------------------------"
echo "Running test_16_performance_100k (benchmark)"
PERF_FILE_100K="${TMP_DATA_DIR}/perf_100k.vcf"

# Generate 100K variant test file
{
    echo -e "##fileformat=VCFv4.2"
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3"
    for i in $(seq 1 100000); do
        gt=$((i % 3))
        echo -e "chr1\t${i}00\t.\tA\tG\t.\tPASS\t.\tGT\t0/${gt}\t${gt}|1\t1/1"
    done
} > "$PERF_FILE_100K"

# Time stdin mode
START=$(python3 -c "import time; print(time.time())")
"$TOOL" -g "0/1" < "$PERF_FILE_100K" > /dev/null 2>&1
END=$(python3 -c "import time; print(time.time())")
STDIN_TIME_100K=$(python3 -c "print(f'{($END - $START)*1000:.1f}')")

# Time mmap mode
START=$(python3 -c "import time; print(time.time())")
"$TOOL" -g "0/1" -i "$PERF_FILE_100K" > /dev/null 2>&1
END=$(python3 -c "import time; print(time.time())")
MMAP_TIME_100K=$(python3 -c "print(f'{($END - $START)*1000:.1f}')")

# Calculate speedup
SPEEDUP=$(python3 -c "print(f'{float($STDIN_TIME_100K)/float($MMAP_TIME_100K):.1f}')")

echo "Performance (100K variants, 3 samples):"
echo "  stdin mode: ${STDIN_TIME_100K}ms"
echo "  mmap mode:  ${MMAP_TIME_100K}ms"
echo "  speedup:    ${SPEEDUP}x"
echo "✅ test_16_performance_100k completed"
echo

###############################################################################
# Test 17: Verify flexible matching (0/1 matches 1/0, 0|1, 1|0)
###############################################################################
echo "-----------------------------------------------------"
echo "Running test_17_flexible_matching"
FLEX_FILE="${TMP_DATA_DIR}/flex_test.vcf"

# Create test file with all genotype orderings
{
    echo -e "##fileformat=VCFv4.2"
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"
    echo -e "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1"
    echo -e "chr1\t200\t.\tA\tG\t.\t.\t.\tGT\t1/0"
    echo -e "chr1\t300\t.\tA\tG\t.\t.\t.\tGT\t0|1"
    echo -e "chr1\t400\t.\tA\tG\t.\t.\t.\tGT\t1|0"
} > "$FLEX_FILE"

# Query for 0/1 - should match all 4
COUNT=$("$TOOL" -g "0/1" -i "$FLEX_FILE" | grep -v "^#" | wc -l | tr -d ' ')
if [ "$COUNT" -eq 4 ]; then
    echo "✅ test_17_flexible_matching passed (matched 4/4 genotypes)"
else
    echo "❌ test_17_flexible_matching failed (matched $COUNT/4)"
    exit 1
fi
echo

###############################################################################
echo "All VCFX_genotype_query tests passed successfully!"
###############################################################################
exit 0
