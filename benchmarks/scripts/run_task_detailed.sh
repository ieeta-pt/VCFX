#!/usr/bin/env bash
# Enhanced benchmarking script with cross-platform performance monitoring
set -e

CMD="$1"
OUTPUT="$2"
PERF_FILE="${OUTPUT}.perf"

echo "Running detailed performance analysis for: $CMD"

# Create temporary file for time output
TEMP_TIME="/tmp/benchmark_time_$$"

# Detect platform and use appropriate time command
if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS: use /usr/bin/time -l
    /usr/bin/time -l /bin/bash -c "$CMD" > "$OUTPUT" 2> "$TEMP_TIME"

    # Parse macOS time output
    if [ -f "$TEMP_TIME" ]; then
        # Extract real time (format: "X.XX real")
        REAL_TIME=$(grep "real" "$TEMP_TIME" | awk '{print $1}')
        USER_TIME=$(grep "user" "$TEMP_TIME" | awk '{print $1}')
        SYS_TIME=$(grep "sys" "$TEMP_TIME" | awk '{print $1}')
        MAX_RSS=$(grep "maximum resident set size" "$TEMP_TIME" | awk '{print $1}')
        PAGE_FAULTS=$(grep "page faults" "$TEMP_TIME" | head -1 | awk '{print $1}')
        VOLUNTARY_SWITCHES=$(grep "voluntary context switches" "$TEMP_TIME" | awk '{print $1}')
        INVOLUNTARY_SWITCHES=$(grep "involuntary context switches" "$TEMP_TIME" | awk '{print $1}')
    fi
else
    # Linux: use /usr/bin/time -v
    /usr/bin/time -v /bin/bash -c "$CMD" > "$OUTPUT" 2> "$TEMP_TIME" || true

    # Parse GNU time output
    if [ -f "$TEMP_TIME" ]; then
        # Extract elapsed time (format: "Elapsed (wall clock) time (h:mm:ss or m:ss): X:XX.XX")
        ELAPSED_LINE=$(grep "Elapsed" "$TEMP_TIME" | head -1)
        if [[ "$ELAPSED_LINE" =~ ([0-9]+):([0-9]+)\.([0-9]+) ]]; then
            MINUTES="${BASH_REMATCH[1]}"
            SECONDS="${BASH_REMATCH[2]}"
            FRACTION="${BASH_REMATCH[3]}"
            REAL_TIME=$(echo "$MINUTES * 60 + $SECONDS.$FRACTION" | bc -l)
        elif [[ "$ELAPSED_LINE" =~ ([0-9]+):([0-9]+):([0-9]+) ]]; then
            HOURS="${BASH_REMATCH[1]}"
            MINUTES="${BASH_REMATCH[2]}"
            SECONDS="${BASH_REMATCH[3]}"
            REAL_TIME=$(echo "$HOURS * 3600 + $MINUTES * 60 + $SECONDS" | bc -l)
        else
            REAL_TIME="0"
        fi
        USER_TIME=$(grep "User time" "$TEMP_TIME" | awk '{print $NF}')
        SYS_TIME=$(grep "System time" "$TEMP_TIME" | awk '{print $NF}')
        MAX_RSS=$(grep "Maximum resident set size" "$TEMP_TIME" | awk '{print $NF}')
        PAGE_FAULTS=$(grep "Minor.*page faults" "$TEMP_TIME" | awk '{print $NF}')
        VOLUNTARY_SWITCHES=$(grep "Voluntary context switches" "$TEMP_TIME" | awk '{print $NF}')
        INVOLUNTARY_SWITCHES=$(grep "Involuntary context switches" "$TEMP_TIME" | awk '{print $NF}')
    fi
fi

# Set defaults for missing values
REAL_TIME=${REAL_TIME:-0}
USER_TIME=${USER_TIME:-0}
SYS_TIME=${SYS_TIME:-0}
MAX_RSS=${MAX_RSS:-0}
PAGE_FAULTS=${PAGE_FAULTS:-0}
VOLUNTARY_SWITCHES=${VOLUNTARY_SWITCHES:-0}
INVOLUNTARY_SWITCHES=${INVOLUNTARY_SWITCHES:-0}

# Convert real time to seconds if it's in M:S format
if [[ "$REAL_TIME" == *":"* ]]; then
    MINUTES=$(echo "$REAL_TIME" | cut -d: -f1)
    SECONDS=$(echo "$REAL_TIME" | cut -d: -f2)
    REAL_TIME_SEC=$(echo "$MINUTES * 60 + $SECONDS" | bc -l)
else
    REAL_TIME_SEC="$REAL_TIME"
fi

# Calculate additional metrics
TOTAL_CPU_TIME=$(echo "${USER_TIME:-0} + ${SYS_TIME:-0}" | bc -l 2>/dev/null || echo "0")

# Convert memory to MB (macOS reports bytes, Linux reports KB)
if [[ "$OSTYPE" == "darwin"* ]]; then
    MAX_RSS_MB=$(echo "scale=2; ${MAX_RSS:-0} / 1024 / 1024" | bc -l 2>/dev/null || echo "0")
else
    MAX_RSS_MB=$(echo "scale=2; ${MAX_RSS:-0} / 1024" | bc -l 2>/dev/null || echo "0")
fi

# Calculate CPU utilization percentage
if [ "$(echo "$REAL_TIME_SEC > 0" | bc -l 2>/dev/null || echo 0)" -eq 1 ]; then
    CPU_UTIL=$(echo "scale=2; ($TOTAL_CPU_TIME / $REAL_TIME_SEC) * 100" | bc -l 2>/dev/null || echo "0")
else
    CPU_UTIL="0"
fi

# Count results (variants or lines)
RESULT_COUNT=$(wc -l < "$OUTPUT" 2>/dev/null | tr -d ' ' || echo "0")
if [ -s "$OUTPUT" ]; then
    # Try to extract number if it's a "Total Variants: N" format
    if grep -q "Total Variants:" "$OUTPUT" 2>/dev/null; then
        RESULT_COUNT=$(grep "Total Variants:" "$OUTPUT" | grep -o '[0-9]*' | head -1)
    fi
fi

# Calculate processing speed (items/second)
if [ "$(echo "$REAL_TIME_SEC > 0" | bc -l 2>/dev/null || echo 0)" -eq 1 ] && [ "${RESULT_COUNT:-0}" -gt 0 ]; then
    PROCESSING_SPEED=$(echo "scale=2; $RESULT_COUNT / $REAL_TIME_SEC" | bc -l 2>/dev/null || echo "0")
else
    PROCESSING_SPEED="0"
fi

# Write detailed performance data
cat > "$PERF_FILE" << EOF
real_time_seconds=$REAL_TIME_SEC
user_time=$USER_TIME
system_time=$SYS_TIME
total_cpu_time=$TOTAL_CPU_TIME
cpu_utilization_percent=$CPU_UTIL
max_memory_mb=$MAX_RSS_MB
max_memory_bytes=$MAX_RSS
page_faults=$PAGE_FAULTS
voluntary_context_switches=$VOLUNTARY_SWITCHES
involuntary_context_switches=$INVOLUNTARY_SWITCHES
result_count=$RESULT_COUNT
processing_speed_per_second=$PROCESSING_SPEED
EOF

# Also write simple time for backward compatibility
echo "$REAL_TIME_SEC" > "${OUTPUT}.time"

echo "Performance data written to $PERF_FILE"

# Display summary
echo "=== PERFORMANCE SUMMARY ==="
echo "Real Time: ${REAL_TIME_SEC}s"
echo "CPU Utilization: ${CPU_UTIL}%"
echo "Peak Memory: ${MAX_RSS_MB}MB"
echo "Results: ${RESULT_COUNT}"
echo "Speed: ${PROCESSING_SPEED} items/second"

# Cleanup
rm -f "$TEMP_TIME"
