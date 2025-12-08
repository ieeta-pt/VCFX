#!/usr/bin/env bash

# Test script for VCFX_haplotype_extractor
# This script tests various aspects of the haplotype extractor functionality

# Set up variables
EXEC="../build/src/VCFX_haplotype_extractor/VCFX_haplotype_extractor"
failures=0

# Check if the executable exists
if [ ! -f "$EXEC" ]; then
    echo "Error: VCFX_haplotype_extractor executable not found at $EXEC"
    echo "Please build the project before running this test."
    exit 1
fi

# Create output directories if they don't exist
mkdir -p out/haplotype_extractor

echo "=== Testing VCFX_haplotype_extractor ==="

# Function to run test and check result
run_test() {
    test_num=$1
    test_name=$2
    input_file=$3
    expected_file=$4
    additional_args=$5
    
    echo "Test $test_num: $test_name"
    $EXEC $additional_args < "$input_file" > "out/haplotype_extractor/test${test_num}.tsv" 2> "out/haplotype_extractor/test${test_num}.err"
    exit_code=$?
    
    if [ $exit_code -ne 0 ]; then
        if [[ "$test_name" == *"expected failure"* ]]; then
            echo "  Test $test_num passed: command failed as expected with exit code $exit_code"
            return 0
        else
            echo "  Test $test_num failed: command exited with code $exit_code"
            cat "out/haplotype_extractor/test${test_num}.err"
            ((failures++))
            return 1
        fi
    elif [[ "$test_name" == *"expected failure"* ]]; then
        # Check if there are warnings in stderr
        if grep -q "Warning:" "out/haplotype_extractor/test${test_num}.err"; then
            echo "  Test $test_num passed: command produced warning messages as expected"
            return 0
        else
            echo "  Test $test_num failed: command succeeded without warnings but expected issues"
            ((failures++))
            return 1
        fi
    fi
    
    if ! diff -q "$expected_file" "out/haplotype_extractor/test${test_num}.tsv" > /dev/null; then
        echo "  Test $test_num failed: output does not match expected"
        echo "  Expected:"
        cat "$expected_file"
        echo "  Actual:"
        cat "out/haplotype_extractor/test${test_num}.tsv"
        ((failures++))
        return 1
    else
        echo "  Test $test_num passed."
        return 0
    fi
}

# Function to test help display
test_help() {
    test_num=$1
    echo "Test $test_num: Help display"
    
    $EXEC --help > "out/haplotype_extractor/test${test_num}.out" 2>&1
    if ! grep -q "Usage:" "out/haplotype_extractor/test${test_num}.out"; then
        echo "  Test $test_num failed: help message does not contain 'Usage:'"
        ((failures++))
    elif ! grep -q "VCFX_haplotype_extractor" "out/haplotype_extractor/test${test_num}.out"; then
        echo "  Test $test_num failed: help message does not contain tool name"
        ((failures++))
    else
        echo "  Test $test_num passed."
    fi
}

# Test 1: Basic functionality - phased genotypes, standard block size
run_test 1 "Basic functionality" \
    "data/haplotype_extractor/basic.vcf" \
    "expected/haplotype_extractor/basic_out.tsv" ""

# Test 2: Large distance between variants - will be grouped in one block due to default threshold
run_test 2 "Large distance between variants" \
    "data/haplotype_extractor/large_distance.vcf" \
    "expected/haplotype_extractor/large_distance_out.tsv" ""

# Test 3: Small block size parameter - should create more blocks
run_test 3 "Small block size parameter" \
    "data/haplotype_extractor/large_distance.vcf" \
    "expected/haplotype_extractor/small_block_out.tsv" "--block-size 50000"

# Test 4: Phase consistency check - implementation doesn't split at phase inconsistencies
run_test 4 "Phase consistency check" \
    "data/haplotype_extractor/phase_inconsistent.vcf" \
    "expected/haplotype_extractor/consistency_check_out.tsv" "--check-phase-consistency"

# Test 5: Mixed phased/unphased genotypes - produces warnings but doesn't fail
run_test 5 "Mixed phased/unphased genotypes (expected failure)" \
    "data/haplotype_extractor/mixed_phasing.vcf" \
    "expected/haplotype_extractor/basic_out.tsv" ""

# Test 6: Missing GT fields - should produce empty haplotype blocks
run_test 6 "Missing GT fields" \
    "data/haplotype_extractor/missing_gt.vcf" \
    "expected/haplotype_extractor/empty_output.tsv" ""

# Test 7: Malformed input - should fail gracefully
run_test 7 "Malformed input (expected failure)" \
    "data/haplotype_extractor/malformed.vcf" \
    "expected/haplotype_extractor/basic_out.tsv" ""

# Test 8: Help message
test_help 8

# Test 9: Streaming mode - basic functionality
echo "Test 9: Streaming mode - basic functionality"
$EXEC --streaming < "data/haplotype_extractor/basic.vcf" > "out/haplotype_extractor/test9.tsv" 2> "out/haplotype_extractor/test9.err"
if diff -q "expected/haplotype_extractor/basic_out.tsv" "out/haplotype_extractor/test9.tsv" > /dev/null; then
    echo "  Test 9 passed."
else
    echo "  Test 9 failed: streaming mode output does not match expected"
    echo "  Expected:"
    cat "expected/haplotype_extractor/basic_out.tsv"
    echo "  Actual:"
    cat "out/haplotype_extractor/test9.tsv"
    ((failures++))
fi

# Test 10: Streaming mode - output matches default mode
echo "Test 10: Streaming vs default mode output consistency"
$EXEC < "data/haplotype_extractor/large_distance.vcf" > "out/haplotype_extractor/test10_default.tsv" 2>/dev/null
$EXEC --streaming < "data/haplotype_extractor/large_distance.vcf" > "out/haplotype_extractor/test10_streaming.tsv" 2>/dev/null
if diff -q "out/haplotype_extractor/test10_default.tsv" "out/haplotype_extractor/test10_streaming.tsv" > /dev/null; then
    echo "  Test 10 passed: streaming output matches default output"
else
    echo "  Test 10 failed: streaming output differs from default"
    echo "  Default:"
    cat "out/haplotype_extractor/test10_default.tsv"
    echo "  Streaming:"
    cat "out/haplotype_extractor/test10_streaming.tsv"
    ((failures++))
fi

# Test 11: Streaming mode with custom block size
echo "Test 11: Streaming mode with custom block size"
$EXEC --streaming --block-size 50000 < "data/haplotype_extractor/large_distance.vcf" > "out/haplotype_extractor/test11.tsv" 2>/dev/null
if diff -q "expected/haplotype_extractor/small_block_out.tsv" "out/haplotype_extractor/test11.tsv" > /dev/null; then
    echo "  Test 11 passed: streaming with block-size works correctly"
else
    echo "  Test 11 failed: streaming with block-size output mismatch"
    ((failures++))
fi

# Test 12: File input mode with -i flag
echo "Test 12: File input mode with -i flag"
$EXEC -i "data/haplotype_extractor/basic.vcf" > "out/haplotype_extractor/test12.tsv" 2> "out/haplotype_extractor/test12.err"
if diff -q "expected/haplotype_extractor/basic_out.tsv" "out/haplotype_extractor/test12.tsv" > /dev/null; then
    echo "  Test 12 passed: -i file input works correctly"
else
    echo "  Test 12 failed: -i file input output mismatch"
    echo "  Expected:"
    cat "expected/haplotype_extractor/basic_out.tsv"
    echo "  Actual:"
    cat "out/haplotype_extractor/test12.tsv"
    ((failures++))
fi

# Test 13: Positional file argument
echo "Test 13: Positional file argument"
$EXEC "data/haplotype_extractor/basic.vcf" > "out/haplotype_extractor/test13.tsv" 2> "out/haplotype_extractor/test13.err"
if diff -q "expected/haplotype_extractor/basic_out.tsv" "out/haplotype_extractor/test13.tsv" > /dev/null; then
    echo "  Test 13 passed: positional file argument works correctly"
else
    echo "  Test 13 failed: positional file argument output mismatch"
    ((failures++))
fi

# Test 14: Quiet mode suppresses warnings
echo "Test 14: Quiet mode suppresses warnings"
$EXEC -q < "data/haplotype_extractor/mixed_phasing.vcf" > "out/haplotype_extractor/test14.tsv" 2> "out/haplotype_extractor/test14.err"
if [ ! -s "out/haplotype_extractor/test14.err" ]; then
    echo "  Test 14 passed: quiet mode suppresses warnings"
else
    echo "  Test 14 failed: quiet mode did not suppress warnings"
    echo "  Stderr content:"
    cat "out/haplotype_extractor/test14.err"
    ((failures++))
fi

# Test 15: Output equivalence - stdin vs mmap
echo "Test 15: Output equivalence - stdin vs mmap"
$EXEC < "data/haplotype_extractor/large_distance.vcf" > "out/haplotype_extractor/test15_stdin.tsv" 2>/dev/null
$EXEC -i "data/haplotype_extractor/large_distance.vcf" > "out/haplotype_extractor/test15_mmap.tsv" 2>/dev/null
if diff -q "out/haplotype_extractor/test15_stdin.tsv" "out/haplotype_extractor/test15_mmap.tsv" > /dev/null; then
    echo "  Test 15 passed: mmap output matches stdin output"
else
    echo "  Test 15 failed: mmap output differs from stdin"
    diff "out/haplotype_extractor/test15_stdin.tsv" "out/haplotype_extractor/test15_mmap.tsv"
    ((failures++))
fi

# Test 16: Streaming mode with mmap file input
echo "Test 16: Streaming mode with mmap file input"
$EXEC --streaming -i "data/haplotype_extractor/basic.vcf" > "out/haplotype_extractor/test16.tsv" 2>/dev/null
if diff -q "expected/haplotype_extractor/basic_out.tsv" "out/haplotype_extractor/test16.tsv" > /dev/null; then
    echo "  Test 16 passed: streaming with mmap works correctly"
else
    echo "  Test 16 failed: streaming with mmap output mismatch"
    ((failures++))
fi

# Test 17: File input with --input long option
echo "Test 17: File input with --input long option"
$EXEC --input "data/haplotype_extractor/basic.vcf" > "out/haplotype_extractor/test17.tsv" 2>/dev/null
if diff -q "expected/haplotype_extractor/basic_out.tsv" "out/haplotype_extractor/test17.tsv" > /dev/null; then
    echo "  Test 17 passed: --input long option works correctly"
else
    echo "  Test 17 failed: --input long option output mismatch"
    ((failures++))
fi

# Test 18: Combined options - streaming + file + block-size
echo "Test 18: Combined options - streaming + file + block-size"
$EXEC --streaming -i "data/haplotype_extractor/large_distance.vcf" --block-size 50000 > "out/haplotype_extractor/test18.tsv" 2>/dev/null
if diff -q "expected/haplotype_extractor/small_block_out.tsv" "out/haplotype_extractor/test18.tsv" > /dev/null; then
    echo "  Test 18 passed: combined streaming + file + block-size works"
else
    echo "  Test 18 failed: combined options output mismatch"
    ((failures++))
fi

# Test 19: Quiet mode with file input
echo "Test 19: Quiet mode with file input"
$EXEC -q -i "data/haplotype_extractor/mixed_phasing.vcf" > "out/haplotype_extractor/test19.tsv" 2> "out/haplotype_extractor/test19.err"
if [ ! -s "out/haplotype_extractor/test19.err" ]; then
    echo "  Test 19 passed: quiet mode with file input suppresses warnings"
else
    echo "  Test 19 failed: quiet mode with file input did not suppress warnings"
    ((failures++))
fi

# Test 20: Help shows new options (-i, -q)
echo "Test 20: Help shows new options"
$EXEC --help > "out/haplotype_extractor/test20.out" 2>&1
if grep -q "\-i" "out/haplotype_extractor/test20.out" && grep -q "\-q" "out/haplotype_extractor/test20.out"; then
    echo "  Test 20 passed: help shows -i and -q options"
else
    echo "  Test 20 failed: help does not show new options"
    echo "  Help output:"
    cat "out/haplotype_extractor/test20.out"
    ((failures++))
fi

# Test 21: Phase consistency check with mmap
echo "Test 21: Phase consistency check with mmap"
$EXEC --check-phase-consistency -i "data/haplotype_extractor/phase_inconsistent.vcf" > "out/haplotype_extractor/test21.tsv" 2>/dev/null
if diff -q "expected/haplotype_extractor/consistency_check_out.tsv" "out/haplotype_extractor/test21.tsv" > /dev/null; then
    echo "  Test 21 passed: phase consistency with mmap works correctly"
else
    echo "  Test 21 failed: phase consistency with mmap output mismatch"
    ((failures++))
fi

# Test 22: Nonexistent file error handling
echo "Test 22: Nonexistent file error handling"
$EXEC -i "data/haplotype_extractor/nonexistent.vcf" > "out/haplotype_extractor/test22.tsv" 2> "out/haplotype_extractor/test22.err"
exit_code=$?
if [ $exit_code -ne 0 ]; then
    echo "  Test 22 passed: nonexistent file causes error exit"
else
    echo "  Test 22 failed: nonexistent file did not cause error"
    ((failures++))
fi

# Report results
if [ $failures -eq 0 ]; then
    echo "All VCFX_haplotype_extractor tests passed!"
    exit 0
else
    echo "$failures VCFX_haplotype_extractor tests failed."
    exit 1
fi 