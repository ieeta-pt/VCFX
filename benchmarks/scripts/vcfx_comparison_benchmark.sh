#!/bin/bash
# =============================================================================
# VCFX vs bcftools/vcftools Performance Comparison Benchmark
# =============================================================================
# Outputs structured CSV with timing, memory, and CPU metrics for comparison
# between VCFX tools and industry-standard VCF processing tools.
#
# Usage: ./vcfx_comparison_benchmark.sh [--small|--full|--all]
#
# Options:
#   --small    Run on small dataset only (~5 min)
#   --full     Run on full 4.3GB dataset (~30 min)
#   --all      Run on all datasets (default)
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="$ROOT_DIR/build/src"
DATA_DIR="$ROOT_DIR/benchmarks/data"
RESULTS_DIR="$ROOT_DIR/benchmarks/results"

# Dataset configurations
TINY_FILE="$DATA_DIR/chr20_tiny.vcf"
SMALL_FILE="$DATA_DIR/chr21_small.vcf"
FULL_FILE="$DATA_DIR/chr21.1kg.phase3.v5a.vcf"

# Benchmark configuration
NUM_RUNS=3
TIMEOUT_SMALL=60      # 1 minute timeout for small files
TIMEOUT_FULL=600      # 10 minute timeout for full files

# Parse arguments
RUN_TINY=false
RUN_SMALL=false
RUN_FULL=false

case "${1:-all}" in
    --tiny)
        RUN_TINY=true
        ;;
    --small)
        RUN_SMALL=true
        ;;
    --full)
        RUN_FULL=true
        ;;
    --all|*)
        RUN_SMALL=true
        RUN_FULL=true
        ;;
esac

# Output file
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
CSV_FILE="$RESULTS_DIR/vcfx_comparison_${TIMESTAMP}.csv"

# Create results directory
mkdir -p "$RESULTS_DIR"

# =============================================================================
# Platform Detection
# =============================================================================
detect_platform() {
    if [[ "$OSTYPE" == "darwin"* ]]; then
        echo "macos"
    else
        echo "linux"
    fi
}

PLATFORM=$(detect_platform)

# =============================================================================
# Tool Availability Check
# =============================================================================
check_tools() {
    echo "=== Checking Tool Availability ==="

    # Check VCFX build
    if [ ! -d "$BUILD_DIR" ]; then
        echo "ERROR: VCFX build directory not found at $BUILD_DIR"
        echo "Please run: cmake --build build"
        exit 1
    fi

    # Check bcftools
    if command -v bcftools &> /dev/null; then
        BCFTOOLS_VERSION=$(bcftools --version | head -1)
        echo "bcftools: $BCFTOOLS_VERSION"
        HAS_BCFTOOLS=true
    else
        echo "bcftools: NOT FOUND (will skip bcftools benchmarks)"
        HAS_BCFTOOLS=false
    fi

    # Check vcftools
    if command -v vcftools &> /dev/null; then
        VCFTOOLS_VERSION=$(vcftools --version 2>&1 | head -1)
        echo "vcftools: $VCFTOOLS_VERSION"
        HAS_VCFTOOLS=true
    else
        echo "vcftools: NOT FOUND (will skip vcftools benchmarks)"
        HAS_VCFTOOLS=false
    fi

    echo ""
}

# =============================================================================
# Timing Function (Platform-Aware)
# =============================================================================
# Returns: real_time_s,user_time_s,sys_time_s,max_memory_mb,status
run_timed() {
    local cmd="$1"
    local timeout_sec="$2"
    local temp_out=$(mktemp)
    local temp_time=$(mktemp)

    # Run with timeout and capture timing
    if [[ "$PLATFORM" == "macos" ]]; then
        # macOS: use /usr/bin/time -l for memory stats
        if timeout "$timeout_sec" /usr/bin/time -l bash -c "$cmd" > /dev/null 2> "$temp_time"; then
            status="success"
        else
            status="failed"
        fi

        # Parse macOS time output
        real_time=$(grep "real" "$temp_time" 2>/dev/null | awk '{print $1}' || echo "0")
        user_time=$(grep "user" "$temp_time" 2>/dev/null | awk '{print $1}' || echo "0")
        sys_time=$(grep "sys" "$temp_time" 2>/dev/null | awk '{print $1}' || echo "0")
        # macOS reports bytes, convert to MB
        max_rss=$(grep "maximum resident set size" "$temp_time" 2>/dev/null | awk '{print $1}' || echo "0")
        max_memory_mb=$(echo "scale=2; $max_rss / 1048576" | bc 2>/dev/null || echo "0")

        # If time parsing failed, use bash time
        if [[ -z "$real_time" || "$real_time" == "0" ]]; then
            local start_time=$(python3 -c "import time; print(time.time())")
            if timeout "$timeout_sec" bash -c "$cmd" > /dev/null 2>&1; then
                status="success"
            else
                status="failed"
            fi
            local end_time=$(python3 -c "import time; print(time.time())")
            real_time=$(python3 -c "print(f'{$end_time - $start_time:.3f}')")
            user_time="$real_time"
            sys_time="0"
            max_memory_mb="0"
        fi
    else
        # Linux: use /usr/bin/time -v
        if timeout "$timeout_sec" /usr/bin/time -v bash -c "$cmd" > /dev/null 2> "$temp_time"; then
            status="success"
        else
            status="failed"
        fi

        # Parse GNU time output
        real_time=$(grep "Elapsed" "$temp_time" 2>/dev/null | sed 's/.*: //' | awk -F: '{if (NF==3) print $1*3600+$2*60+$3; else if (NF==2) print $1*60+$2; else print $1}' || echo "0")
        user_time=$(grep "User time" "$temp_time" 2>/dev/null | awk '{print $NF}' || echo "0")
        sys_time=$(grep "System time" "$temp_time" 2>/dev/null | awk '{print $NF}' || echo "0")
        # Linux reports KB, convert to MB
        max_rss=$(grep "Maximum resident" "$temp_time" 2>/dev/null | awk '{print $NF}' || echo "0")
        max_memory_mb=$(echo "scale=2; $max_rss / 1024" | bc 2>/dev/null || echo "0")
    fi

    rm -f "$temp_out" "$temp_time"

    # Calculate CPU percentage
    if [[ "$real_time" != "0" && "$real_time" != "" ]]; then
        cpu_percent=$(python3 -c "print(f'{(float($user_time) + float($sys_time)) / float($real_time) * 100:.1f}')" 2>/dev/null || echo "100")
    else
        cpu_percent="100"
    fi

    echo "$real_time,$user_time,$sys_time,$max_memory_mb,$cpu_percent,$status"
}

# =============================================================================
# Simple Timing Function (Fallback)
# =============================================================================
run_simple_timed() {
    local cmd="$1"
    local timeout_sec="$2"

    local start_time=$(python3 -c "import time; print(time.time())")

    if timeout "$timeout_sec" bash -c "$cmd" > /dev/null 2>&1; then
        local status="success"
    else
        local status="failed"
    fi

    local end_time=$(python3 -c "import time; print(time.time())")
    local real_time=$(python3 -c "print(f'{$end_time - $start_time:.3f}')")

    echo "$real_time,$real_time,0,0,100,$status"
}

# =============================================================================
# Benchmark Runner
# =============================================================================
run_benchmark() {
    local operation="$1"
    local tool="$2"
    local dataset_name="$3"
    local data_file="$4"
    local file_size_mb="$5"
    local variants="$6"
    local samples="$7"
    local cmd="$8"
    local timeout_sec="$9"

    echo -n "  $tool... "

    for run_num in $(seq 1 $NUM_RUNS); do
        result=$(run_simple_timed "$cmd" "$timeout_sec")

        IFS=',' read -r real_time user_time sys_time max_memory_mb cpu_percent status <<< "$result"

        # Write to CSV
        echo "$operation,$tool,$dataset_name,$file_size_mb,$variants,$samples,$run_num,$real_time,$user_time,$sys_time,$max_memory_mb,$cpu_percent,$status,\"$cmd\"" >> "$CSV_FILE"
    done

    # Extract last status for display
    echo "$status (${real_time}s)"
}

# =============================================================================
# Get File Stats
# =============================================================================
get_file_stats() {
    local file="$1"

    if [ ! -f "$file" ]; then
        echo "0,0,0"
        return
    fi

    local size_bytes=$(stat -f%z "$file" 2>/dev/null || stat -c%s "$file" 2>/dev/null || echo "0")
    local size_mb=$(echo "scale=2; $size_bytes / 1048576" | bc)

    # Count variants (lines not starting with #)
    local variants=$(grep -c "^[^#]" "$file" 2>/dev/null || echo "0")

    # Count samples from #CHROM header
    local samples=$(grep "^#CHROM" "$file" 2>/dev/null | awk -F'\t' '{print NF-9}' || echo "0")

    echo "$size_mb,$variants,$samples"
}

# =============================================================================
# Main Benchmark Suite
# =============================================================================
run_suite() {
    local dataset_name="$1"
    local data_file="$2"
    local timeout_sec="$3"

    echo ""
    echo "=== Benchmarking: $dataset_name ==="
    echo "File: $data_file"

    if [ ! -f "$data_file" ]; then
        echo "WARNING: File not found, skipping..."
        return
    fi

    # Get file stats
    IFS=',' read -r file_size_mb variants samples <<< "$(get_file_stats "$data_file")"
    echo "Size: ${file_size_mb}MB, Variants: $variants, Samples: $samples"
    echo ""

    # Warm up file cache
    echo "Warming up file cache..."
    cat "$data_file" > /dev/null 2>&1
    echo ""

    # -------------------------------------------------------------------------
    # 1. Variant Counting
    # -------------------------------------------------------------------------
    echo "--- Variant Counting ---"
    run_benchmark "count" "vcfx" "$dataset_name" "$data_file" "$file_size_mb" "$variants" "$samples" \
        "$BUILD_DIR/VCFX_variant_counter/VCFX_variant_counter < $data_file" "$timeout_sec"

    if $HAS_BCFTOOLS; then
        run_benchmark "count" "bcftools" "$dataset_name" "$data_file" "$file_size_mb" "$variants" "$samples" \
            "bcftools view -H $data_file | wc -l" "$timeout_sec"
    fi

    # -------------------------------------------------------------------------
    # 2. Field Extraction
    # -------------------------------------------------------------------------
    echo "--- Field Extraction ---"
    run_benchmark "field_extract" "vcfx" "$dataset_name" "$data_file" "$file_size_mb" "$variants" "$samples" \
        "$BUILD_DIR/VCFX_field_extractor/VCFX_field_extractor --fields CHROM,POS,REF,ALT < $data_file" "$timeout_sec"

    if $HAS_BCFTOOLS; then
        run_benchmark "field_extract" "bcftools" "$dataset_name" "$data_file" "$file_size_mb" "$variants" "$samples" \
            "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $data_file" "$timeout_sec"
    fi

    # -------------------------------------------------------------------------
    # 3. Quality Filtering
    # -------------------------------------------------------------------------
    echo "--- Quality Filtering (QUAL>=30) ---"
    run_benchmark "quality_filter" "vcfx" "$dataset_name" "$data_file" "$file_size_mb" "$variants" "$samples" \
        "$BUILD_DIR/VCFX_phred_filter/VCFX_phred_filter -p 30 < $data_file" "$timeout_sec"

    if $HAS_BCFTOOLS; then
        run_benchmark "quality_filter" "bcftools" "$dataset_name" "$data_file" "$file_size_mb" "$variants" "$samples" \
            "bcftools view -i 'QUAL>=30' $data_file" "$timeout_sec"
    fi

    if $HAS_VCFTOOLS; then
        run_benchmark "quality_filter" "vcftools" "$dataset_name" "$data_file" "$file_size_mb" "$variants" "$samples" \
            "vcftools --vcf $data_file --minQ 30 --recode --stdout 2>/dev/null" "$timeout_sec"
    fi

    # -------------------------------------------------------------------------
    # 4. Allele Frequency Calculation
    # -------------------------------------------------------------------------
    echo "--- Allele Frequency Calculation ---"
    run_benchmark "allele_freq" "vcfx" "$dataset_name" "$data_file" "$file_size_mb" "$variants" "$samples" \
        "$BUILD_DIR/VCFX_allele_freq_calc/VCFX_allele_freq_calc < $data_file" "$timeout_sec"

    if $HAS_VCFTOOLS; then
        run_benchmark "allele_freq" "vcftools" "$dataset_name" "$data_file" "$file_size_mb" "$variants" "$samples" \
            "vcftools --vcf $data_file --freq --stdout 2>/dev/null" "$timeout_sec"
    fi

    # -------------------------------------------------------------------------
    # 5. Validation
    # -------------------------------------------------------------------------
    echo "--- VCF Validation ---"

    # VCFX with mmap (fast path)
    run_benchmark "validate" "vcfx_mmap" "$dataset_name" "$data_file" "$file_size_mb" "$variants" "$samples" \
        "$BUILD_DIR/VCFX_validator/VCFX_validator -i $data_file --no-dup-check" "$timeout_sec"

    # VCFX with stdin (slow path)
    run_benchmark "validate" "vcfx_stdin" "$dataset_name" "$data_file" "$file_size_mb" "$variants" "$samples" \
        "$BUILD_DIR/VCFX_validator/VCFX_validator --no-dup-check < $data_file" "$timeout_sec"

    if $HAS_BCFTOOLS; then
        run_benchmark "validate" "bcftools" "$dataset_name" "$data_file" "$file_size_mb" "$variants" "$samples" \
            "bcftools view $data_file" "$timeout_sec"
    fi

    # -------------------------------------------------------------------------
    # 6. Missing Data Detection
    # -------------------------------------------------------------------------
    echo "--- Missing Data Detection ---"
    run_benchmark "missing" "vcfx" "$dataset_name" "$data_file" "$file_size_mb" "$variants" "$samples" \
        "$BUILD_DIR/VCFX_missing_detector/VCFX_missing_detector < $data_file" "$timeout_sec"

    if $HAS_VCFTOOLS; then
        run_benchmark "missing" "vcftools" "$dataset_name" "$data_file" "$file_size_mb" "$variants" "$samples" \
            "vcftools --vcf $data_file --missing-indv --stdout 2>/dev/null" "$timeout_sec"
    fi

    # -------------------------------------------------------------------------
    # 7. Region Subsetting (first 10M bases)
    # -------------------------------------------------------------------------
    echo "--- Region Subsetting (first 10Mb) ---"
    run_benchmark "region_subset" "vcfx" "$dataset_name" "$data_file" "$file_size_mb" "$variants" "$samples" \
        "$BUILD_DIR/VCFX_position_subsetter/VCFX_position_subsetter --region 21:1-10000000 < $data_file" "$timeout_sec"

    if $HAS_BCFTOOLS; then
        run_benchmark "region_subset" "bcftools" "$dataset_name" "$data_file" "$file_size_mb" "$variants" "$samples" \
            "bcftools view -r 21:1-10000000 $data_file" "$timeout_sec"
    fi

    # -------------------------------------------------------------------------
    # 8. Indexing
    # -------------------------------------------------------------------------
    echo "--- Indexing ---"
    run_benchmark "index" "vcfx" "$dataset_name" "$data_file" "$file_size_mb" "$variants" "$samples" \
        "$BUILD_DIR/VCFX_indexer/VCFX_indexer $data_file" "$timeout_sec"

    echo ""
}

# =============================================================================
# Main Execution
# =============================================================================
main() {
    echo "=============================================="
    echo "  VCFX vs bcftools/vcftools Benchmark"
    echo "=============================================="
    echo ""
    echo "Date: $(date)"
    echo "Platform: $PLATFORM"
    echo "Results: $CSV_FILE"
    echo ""

    # Check tools
    check_tools

    # Create CSV header
    echo "operation,tool,dataset,file_size_mb,variants,samples,run_num,real_time_s,user_time_s,sys_time_s,max_memory_mb,cpu_percent,status,command" > "$CSV_FILE"

    # Run benchmarks
    if $RUN_TINY && [ -f "$TINY_FILE" ]; then
        run_suite "tiny" "$TINY_FILE" "$TIMEOUT_SMALL"
    fi

    if $RUN_SMALL && [ -f "$SMALL_FILE" ]; then
        run_suite "small" "$SMALL_FILE" "$TIMEOUT_SMALL"
    fi

    if $RUN_FULL && [ -f "$FULL_FILE" ]; then
        run_suite "full" "$FULL_FILE" "$TIMEOUT_FULL"
    fi

    echo ""
    echo "=============================================="
    echo "  Benchmark Complete!"
    echo "=============================================="
    echo ""
    echo "Results saved to: $CSV_FILE"
    echo ""

    # Print summary
    echo "=== Summary ==="
    python3 << PYSCRIPT
import csv
from collections import defaultdict

results = defaultdict(lambda: defaultdict(list))

with open("$CSV_FILE") as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row['status'] == 'success':
            key = (row['operation'], row['dataset'])
            results[key][row['tool']].append(float(row['real_time_s']))

print(f"{'Operation':<20} {'Dataset':<10} {'Tool':<15} {'Avg Time (s)':<12} {'Speedup':<10}")
print("-" * 75)

for (op, dataset), tools in sorted(results.items()):
    # Find baseline (bcftools if available, else first tool)
    baseline_tool = 'bcftools' if 'bcftools' in tools else list(tools.keys())[0]
    baseline_time = sum(tools[baseline_tool]) / len(tools[baseline_tool]) if baseline_tool in tools else 1

    for tool, times in sorted(tools.items()):
        avg_time = sum(times) / len(times)
        if baseline_time > 0:
            speedup = baseline_time / avg_time
            speedup_str = f"{speedup:.2f}x" if speedup != 1.0 else "baseline"
        else:
            speedup_str = "-"
        print(f"{op:<20} {dataset:<10} {tool:<15} {avg_time:<12.3f} {speedup_str:<10}")
    print()
PYSCRIPT
}

# Run main
main "$@"
