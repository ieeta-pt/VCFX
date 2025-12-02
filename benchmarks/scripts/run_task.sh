#!/usr/bin/env bash
# Simple helper to run a command and measure execution time.
set -e
CMD="$1"
OUTPUT="$2"

# Use higher precision timing with nanoseconds
START=$(date +%s.%N)
/bin/bash -c "$CMD" > "$OUTPUT"
END=$(date +%s.%N)

# Calculate duration in seconds with decimal precision
DURATION=$(echo "$END - $START" | bc -l)
echo "$DURATION" > "$OUTPUT.time"
