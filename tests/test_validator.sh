#!/usr/bin/env bash

# Create output directories
mkdir -p out data expected

# Set the path to the executable
EXEC="../build/src/VCFX_validator/VCFX_validator"

# Check if executable exists
if [ ! -f "$EXEC" ]; then
    echo "Error: $EXEC does not exist. Please build the project first."
    exit 1
fi

# Function to run a test with expected success (stdin mode)
run_test_success() {
    local test_num=$1
    local description=$2
    local input_file=$3
    local opts="$4"

    echo -n "Test $test_num: $description... "

    # Run the command using process substitution
    local output
    local exit_code
    output=$($EXEC $opts < "$input_file" 2>&1)
    exit_code=$?

    if [ $exit_code -eq 0 ]; then
        echo "PASSED"
        return 0
    else
        echo "FAILED (Expected success, got failure)"
        echo "Exit code: $exit_code"
        echo "Command output: $output"
        return 1
    fi
}

# Function to run a test with expected success (file path mode - mmap)
run_test_success_file() {
    local test_num=$1
    local description=$2
    local input_file=$3
    local opts="$4"

    echo -n "Test $test_num: $description... "

    # Run the command with file path argument (uses mmap)
    local output
    local exit_code
    output=$($EXEC $opts "$input_file" 2>&1)
    exit_code=$?

    if [ $exit_code -eq 0 ]; then
        echo "PASSED"
        return 0
    else
        echo "FAILED (Expected success, got failure)"
        echo "Exit code: $exit_code"
        echo "Command output: $output"
        return 1
    fi
}

# Function to run a test with expected failure
run_test_failure() {
    local test_num=$1
    local description=$2
    local input_file=$3
    local expected_error=$4
    local opts="$5"

    echo -n "Test $test_num: $description... "

    # Run the command using process substitution
    local output
    local exit_code
    output=$($EXEC $opts < "$input_file" 2>&1)
    exit_code=$?

    if [ $exit_code -ne 0 ]; then
        # Check if the output contains the expected error message
        if echo "$output" | grep -q "$expected_error"; then
            echo "PASSED (Expected failure)"
            return 0
        else
            echo "FAILED (Expected error message not found)"
            echo "Expected to contain: $expected_error"
            echo "Exit code: $exit_code"
            echo "Actual output: $output"
            return 1
        fi
    else
        echo "FAILED (Expected failure, got success)"
        echo "Exit code: $exit_code"
        echo "Command output: $output"
        return 1
    fi
}

# Function to run a test with expected failure (file path mode)
run_test_failure_file() {
    local test_num=$1
    local description=$2
    local input_file=$3
    local expected_error=$4
    local opts="$5"

    echo -n "Test $test_num: $description... "

    # Run the command with file path argument
    local output
    local exit_code
    output=$($EXEC $opts "$input_file" 2>&1)
    exit_code=$?

    if [ $exit_code -ne 0 ]; then
        # Check if the output contains the expected error message
        if echo "$output" | grep -q "$expected_error"; then
            echo "PASSED (Expected failure)"
            return 0
        else
            echo "FAILED (Expected error message not found)"
            echo "Expected to contain: $expected_error"
            echo "Exit code: $exit_code"
            echo "Actual output: $output"
            return 1
        fi
    else
        echo "FAILED (Expected failure, got success)"
        echo "Exit code: $exit_code"
        echo "Command output: $output"
        return 1
    fi
}

# Function to verify output contains expected text
run_test_output_contains() {
    local test_num=$1
    local description=$2
    local input_file=$3
    local expected_text=$4
    local opts="$5"

    echo -n "Test $test_num: $description... "

    local output
    local exit_code
    output=$($EXEC $opts "$input_file" 2>&1)
    exit_code=$?

    if [ $exit_code -eq 0 ] && echo "$output" | grep -q "$expected_text"; then
        echo "PASSED"
        return 0
    else
        echo "FAILED"
        echo "Expected output to contain: $expected_text"
        echo "Exit code: $exit_code"
        echo "Actual output: $output"
        return 1
    fi
}

# Create test files (same as in test_validator.sh)
# Valid VCF file - basic
cat > data/valid.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	100	rs123	A	T	60	PASS	NS=2;DP=10	GT:GQ:DP	0/1:48:8	0/0:43:5
chr2	200	rs456	G	C	80	PASS	NS=2;DP=15	GT:GQ:DP	0/1:56:12	1/1:60:9
chr3	300	.	T	A	90	PASS	NS=2;DP=20	GT:GQ:DP	0/0:50:10	0/1:65:10
EOF

# Invalid VCF - missing CHROM header
cat > data/no_chrom_header.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
chr1	100	rs123	A	T	60	PASS	NS=2;DP=10
chr2	200	rs456	G	C	80	PASS	NS=2;DP=15
EOF

# Invalid VCF - malformed meta line
cat > data/malformed_meta.vcf << EOF
##fileformat=VCFv4.2
#INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	rs123	A	T	60	PASS	NS=2;DP=10
EOF

# Invalid VCF - data before CHROM header
cat > data/data_before_header.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
chr1	100	rs123	A	T	60	PASS	NS=2;DP=10
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr2	200	rs456	G	C	80	PASS	NS=2;DP=15
EOF

# Invalid VCF - too few columns
cat > data/too_few_columns.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER
chr1	100	rs123	A	T	60	PASS
EOF

# Invalid VCF - invalid POS (negative)
cat > data/invalid_pos.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	-10	rs123	A	T	60	PASS	NS=2;DP=10
EOF

# Invalid VCF - empty REF
cat > data/empty_ref.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	rs123		T	60	PASS	NS=2;DP=10
EOF

# Invalid VCF - empty ALT
cat > data/empty_alt.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
#CHROM	POS	ID	REF		QUAL	FILTER	INFO
chr1	100	rs123	A		60	PASS	NS=2;DP=10
EOF

# Invalid VCF - negative QUAL
cat > data/negative_qual.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	rs123	A	T	-60	PASS	NS=2;DP=10
EOF

# Invalid VCF - empty FILTER
cat > data/empty_filter.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
#CHROM	POS	ID	REF	ALT	QUAL		INFO
chr1	100	rs123	A	T	60		NS=2;DP=10
EOF

# Valid VCF - with missing values represented by dots
cat > data/valid_with_dots.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	T	.	PASS	.
chr2	200	rs456	G	C	80	PASS	NS=2;DP=15
EOF

# Header has one sample column but a data line includes two sample columns
cat > data/mismatched_columns.vcf << EOF
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
EOF
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n' >> data/mismatched_columns.vcf
printf 'chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\t0/0\n' >> data/mismatched_columns.vcf

# FORMAT expects two entries but sample has three
cat > data/format_mismatch.vcf << EOF
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
EOF
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n' >> data/format_mismatch.vcf
printf 'chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT:DP\t0/1:30:7\n' >> data/format_mismatch.vcf

# Undefined INFO field
cat > data/undefined_info.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
EOF
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' >> data/undefined_info.vcf
printf 'chr1\t100\t.\tA\tT\t60\tPASS\tDP=5;FOO=1\n' >> data/undefined_info.vcf

# Undefined FORMAT field
cat > data/undefined_format.vcf << EOF
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
EOF
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n' >> data/undefined_format.vcf
printf 'chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT:XY\t0/1:10\n' >> data/undefined_format.vcf

# Invalid genotype string
cat > data/invalid_genotype.vcf << EOF
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
EOF
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n' >> data/invalid_genotype.vcf
printf 'chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0//1\n' >> data/invalid_genotype.vcf

# Invalid REF/ALT bases
cat > data/invalid_bases.vcf << EOF
##fileformat=VCFv4.2
EOF
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' >> data/invalid_bases.vcf
printf 'chr1\t100\t.\tA\tR\t60\tPASS\t.\n' >> data/invalid_bases.vcf

# Duplicate records
cat > data/duplicate_records.vcf << EOF
##fileformat=VCFv4.2
EOF
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' >> data/duplicate_records.vcf
printf 'chr1\t100\t.\tA\tT\t60\tPASS\t.\n' >> data/duplicate_records.vcf
printf 'chr1\t100\t.\tA\tT\t60\tPASS\t.\n' >> data/duplicate_records.vcf

# gzipped input
if [ ! -f data/valid.vcf.gz ]; then
    gzip -c data/valid.vcf > data/valid.vcf.gz
fi

# ============================================================================
# Additional test files for comprehensive testing
# ============================================================================

# Multi-allelic variants
cat > data/multiallelic.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr1	100	.	A	T,C	60	PASS	DP=30	GT	0/1
chr1	200	.	G	A,T,C	80	PASS	DP=40	GT	1/2
chr1	300	.	AT	A,ATT	70	PASS	DP=25	GT	0/2
EOF

# Phased genotypes
cat > data/phased_genotypes.vcf << EOF
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	100	.	A	T	60	PASS	.	GT	0|1	1|0
chr1	200	.	G	C	80	PASS	.	GT	0|0	1|1
chr1	300	.	T	A	90	PASS	.	GT	1|1	0|1
EOF

# Large position value (edge case)
cat > data/large_pos.vcf << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	2147483647	.	A	T	60	PASS	.
EOF

# VCF with many samples
cat > data/many_samples.vcf << EOF
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
EOF
# Build header with 100 samples
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT' >> data/many_samples.vcf
for i in $(seq 1 100); do printf '\tSAMPLE%d' $i >> data/many_samples.vcf; done
printf '\n' >> data/many_samples.vcf
# Build data line with 100 genotypes
printf 'chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT' >> data/many_samples.vcf
for i in $(seq 1 100); do printf '\t0/1' >> data/many_samples.vcf; done
printf '\n' >> data/many_samples.vcf

# Empty file (should fail)
> data/empty.vcf

# File with only header (no data lines - valid)
cat > data/header_only.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
EOF

# VCF with Windows line endings (CRLF)
printf '##fileformat=VCFv4.2\r\n' > data/windows_crlf.vcf
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\r\n' >> data/windows_crlf.vcf
printf 'chr1\t100\t.\tA\tT\t60\tPASS\t.\r\n' >> data/windows_crlf.vcf

# VCF with lowercase bases (should be valid - case insensitive)
cat > data/lowercase_bases.vcf << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	a	t	60	PASS	.
chr1	200	.	ACGT	acgt	80	PASS	.
EOF

# VCF with N bases (valid)
cat > data/n_bases.vcf << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	N	A	60	PASS	.
chr1	200	.	ANA	T	80	PASS	.
chr1	300	.	G	NNN	70	PASS	.
EOF

# VCF with scientific notation QUAL
cat > data/scientific_qual.vcf << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	T	1.5e2	PASS	.
chr1	200	.	G	C	2.3E-1	PASS	.
chr1	300	.	T	A	0.0	PASS	.
EOF

# VCF with complex INFO fields
cat > data/complex_info.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total alleles">
##INFO=<ID=DB,Number=0,Type=Flag,Description="In dbSNP">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	rs123	A	T	60	PASS	AC=5;AF=0.25;AN=20;DB
chr1	200	.	G	C,T	80	PASS	AC=3,2;AF=0.15,0.10;AN=20
EOF

# Invalid: multi-allelic with empty allele
cat > data/invalid_multiallelic.vcf << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	T,	60	PASS	.
EOF

# Invalid: POS = 0
cat > data/pos_zero.vcf << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	0	.	A	T	60	PASS	.
EOF

# Polyploid genotypes (valid)
cat > data/polyploid.vcf << EOF
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr1	100	.	A	T	60	PASS	.	GT	0/0/1
chr1	200	.	G	C	80	PASS	.	GT	0/1/1/2
EOF

# Run each test separately and track failures
failures=0

# Test 1
run_test_success 1 "Basic valid VCF" "data/valid.vcf"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 2
run_test_failure 2 "Missing CHROM header" "data/no_chrom_header.vcf" "data line encountered before #CHROM"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 3
run_test_failure 3 "Malformed meta line" "data/malformed_meta.vcf" "is a header line but neither starts with '##' nor is a #CHROM"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 4
run_test_failure 4 "Data before CHROM header" "data/data_before_header.vcf" "data line encountered before #CHROM"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 5
run_test_failure 5 "Too few columns" "data/too_few_columns.vcf" "has <8 columns"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 6
run_test_failure 6 "Invalid position" "data/invalid_pos.vcf" "POS must be >0"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 7
run_test_failure 7 "Empty REF" "data/empty_ref.vcf" "REF is empty"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 8
run_test_failure 8 "Empty ALT" "data/empty_alt.vcf" "ALT is empty"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 9
run_test_failure 9 "Negative QUAL" "data/negative_qual.vcf" "negative QUAL"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 10
run_test_failure 10 "Empty FILTER" "data/empty_filter.vcf" "FILTER is empty"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 11
run_test_success 11 "Valid VCF with dots" "data/valid_with_dots.vcf"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 12 - Help message
echo -n "Test 12: Help message... "
$EXEC --help > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "PASSED"
else
    echo "FAILED"
    failures=$((failures + 1))
fi

# Test 13 - strict mode valid file
run_test_success 13 "Strict valid VCF" "data/valid.vcf" "--strict"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 14 - mismatched columns in strict mode
run_test_failure 14 "Strict mismatched columns" "data/mismatched_columns.vcf" "columns" "--strict"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 15 - FORMAT/sample mismatch in strict mode
run_test_failure 15 "Strict format mismatch" "data/format_mismatch.vcf" "FORMAT" "--strict"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 16 - undefined INFO field
run_test_failure 16 "Undefined INFO" "data/undefined_info.vcf" "INFO field" "--strict"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 17 - undefined FORMAT field
run_test_failure 17 "Undefined FORMAT" "data/undefined_format.vcf" "FORMAT field" "--strict"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 18 - invalid genotype
run_test_failure 18 "Invalid genotype" "data/invalid_genotype.vcf" "invalid genotype" "--strict"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 19 - invalid bases
run_test_failure 19 "Invalid bases" "data/invalid_bases.vcf" "invalid characters" "--strict"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 20 - duplicate detection
run_test_failure 20 "Duplicate records" "data/duplicate_records.vcf" "duplicate variant" "--strict --report-dups"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 21 - gzipped input
run_test_success 21 "Gzip input" "data/valid.vcf.gz"
[ $? -ne 0 ] && failures=$((failures + 1))

# ============================================================================
# New comprehensive tests (22+)
# ============================================================================

# Test 22 - File path mode (mmap) - valid file
run_test_success_file 22 "File path mode (mmap)" "data/valid.vcf"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 23 - File path mode with strict
run_test_success_file 23 "File path mode strict" "data/valid.vcf" "--strict"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 24 - File path mode failure detection
run_test_failure_file 24 "File path mode invalid" "data/invalid_pos.vcf" "POS must be >0"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 25 - Multi-allelic variants
run_test_success 25 "Multi-allelic variants" "data/multiallelic.vcf"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 26 - Phased genotypes (pipe separator)
run_test_success 26 "Phased genotypes" "data/phased_genotypes.vcf"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 27 - Large POS value (INT_MAX)
run_test_success 27 "Large POS value" "data/large_pos.vcf"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 28 - Many samples (100)
run_test_success 28 "Many samples (100)" "data/many_samples.vcf"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 29 - Empty file
run_test_failure 29 "Empty file" "data/empty.vcf" "no #CHROM line found"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 30 - Header only (no data lines) - should fail by default
run_test_failure 30 "Header only (no data)" "data/header_only.vcf" "no variant records"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 30b - Header only with --allow-empty flag - should pass
run_test_success "30b" "Header only with --allow-empty" "data/header_only.vcf" "--allow-empty"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 31 - Windows line endings (CRLF)
run_test_success 31 "Windows CRLF endings" "data/windows_crlf.vcf"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 32 - Lowercase bases (case insensitive)
run_test_success 32 "Lowercase bases" "data/lowercase_bases.vcf"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 33 - N bases (valid ambiguous)
run_test_success 33 "N bases (ambiguous)" "data/n_bases.vcf"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 34 - Scientific notation QUAL
run_test_success 34 "Scientific notation QUAL" "data/scientific_qual.vcf"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 35 - Complex INFO fields
run_test_success 35 "Complex INFO fields" "data/complex_info.vcf"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 36 - Invalid multi-allelic (empty allele)
run_test_failure 36 "Invalid multi-allelic" "data/invalid_multiallelic.vcf" "invalid characters"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 37 - POS = 0 (invalid)
run_test_failure 37 "POS zero" "data/pos_zero.vcf" "POS must be >0"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 38 - Polyploid genotypes
run_test_success 38 "Polyploid genotypes" "data/polyploid.vcf"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 39 - Output report contains variant count
run_test_output_contains 39 "Report shows variant count" "data/valid.vcf" "Variant records:"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 40 - Output report contains sample count
run_test_output_contains 40 "Report shows sample count" "data/valid.vcf" "Samples:"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 41 - Output report shows PASSED status
run_test_output_contains 41 "Report shows PASSED" "data/valid.vcf" "Status: PASSED"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 42 - Nonexistent file error
echo -n "Test 42: Nonexistent file error... "
output=$($EXEC "data/nonexistent_file.vcf" 2>&1)
exit_code=$?
if [ $exit_code -ne 0 ] && echo "$output" | grep -qi "cannot open\|error\|not found"; then
    echo "PASSED (Expected failure)"
else
    echo "FAILED"
    echo "Exit code: $exit_code"
    echo "Output: $output"
    failures=$((failures + 1))
fi

# ============================================================================
# New tests for GATK-like validations (43+)
# ============================================================================

# Test 43 - Malformed header with swapped Type/Number (Type is empty)
cat > data/malformed_type.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=NODE,Number=.,Integer,Type=,Description="Bad header with swapped fields">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	T	60	PASS	NODE=1
EOF
run_test_failure 43 "Malformed header Type" "data/malformed_type.vcf" "Type"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 43b - Header with unrecognized Type value
cat > data/unrecognized_type.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integerr,Description="Typo in Type">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	T	60	PASS	DP=10
EOF
run_test_failure "43b" "Unrecognized header Type" "data/unrecognized_type.vcf" "invalid Type"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 44 - Header with missing Type field
cat > data/missing_type.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Description="Missing Type field">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	T	60	PASS	DP=10
EOF
run_test_failure 44 "Missing header Type" "data/missing_type.vcf" "missing Type"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 45 - Header with invalid Number field
cat > data/invalid_number.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=invalid,Type=Integer,Description="Invalid Number">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	T	60	PASS	DP=10
EOF
run_test_failure 45 "Invalid header Number" "data/invalid_number.vcf" "invalid Number"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 46 - ALT allele not observed in genotypes (strict mode)
cat > data/alt_not_observed.vcf << EOF
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	100	.	A	T,C	60	PASS	.	GT	0/0	0/0
EOF
run_test_failure 46 "ALT not observed (strict)" "data/alt_not_observed.vcf" "ALT allele" "--strict"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 47 - ALT allele not observed produces warning (non-strict)
echo -n "Test 47: ALT not observed warning... "
output=$($EXEC < data/alt_not_observed.vcf 2>&1)
exit_code=$?
if [ $exit_code -eq 0 ] && echo "$output" | grep -q "Warning.*ALT allele"; then
    echo "PASSED"
else
    echo "FAILED"
    echo "Exit code: $exit_code"
    echo "Output: $output"
    failures=$((failures + 1))
fi

# Test 48 - ALT alleles all observed (should pass)
cat > data/alt_observed.vcf << EOF
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	100	.	A	T,C	60	PASS	.	GT	0/1	0/2
chr1	200	.	G	A	80	PASS	.	GT	1/1	0/1
EOF
run_test_success 48 "ALT alleles observed" "data/alt_observed.vcf" "--strict"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 49 - Valid header with all Type values
cat > data/valid_types.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=I1,Number=1,Type=Integer,Description="Integer type">
##INFO=<ID=F1,Number=1,Type=Float,Description="Float type">
##INFO=<ID=FL,Number=0,Type=Flag,Description="Flag type">
##INFO=<ID=C1,Number=1,Type=Character,Description="Character type">
##INFO=<ID=S1,Number=1,Type=String,Description="String type">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr1	100	.	A	T	60	PASS	I1=5;F1=0.5;FL;C1=X;S1=test	GT	0/1
EOF
run_test_success 49 "Valid header types" "data/valid_types.vcf"
[ $? -ne 0 ] && failures=$((failures + 1))

# Test 50 - Valid Number values (A, R, G, .)
cat > data/valid_numbers.vcf << EOF
##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AF,Number=R,Type=Float,Description="Allele freq">
##INFO=<ID=GP,Number=G,Type=Float,Description="Genotype prob">
##INFO=<ID=DP,Number=.,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
chr1	100	.	A	T	60	PASS	AC=5;AF=0.1,0.9;GP=0.1,0.2,0.7;DP=10	GT	0/1
EOF
run_test_success 50 "Valid Number values" "data/valid_numbers.vcf"
[ $? -ne 0 ] && failures=$((failures + 1))

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "============================================"
if [ $failures -eq 0 ]; then
    echo "All 51 tests for VCFX_validator passed!"
    exit 0
else
    echo "$failures test(s) failed out of 51."
    exit 1
fi