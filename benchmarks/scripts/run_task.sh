#!/usr/bin/env bash
# Simple helper to run a command and measure execution time.
set -e
CMD="$1"
OUTPUT="$2"
START=$(date +%s)
/bin/bash -c "$CMD" > "$OUTPUT"
END=$(date +%s)
DURATION=$((END-START))
echo "$DURATION" > "$OUTPUT.time"
