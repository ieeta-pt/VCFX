#!/usr/bin/env bash

set -e
set -o pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"
BUILD_DIR="${ROOT_DIR}/build/python_bindings"

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

cmake -DPYTHON_BINDINGS=ON ../..
make -j

# Extract the Python executable from CMake cache to ensure we use the same Python
# that the extension was built against
PYTHON_CONFIG=$(grep '_Python3_CONFIG:INTERNAL' CMakeCache.txt 2>/dev/null | cut -d= -f2)
if [ -n "$PYTHON_CONFIG" ] && [ -x "$PYTHON_CONFIG" ]; then
    # python3.13-config -> python3.13
    PYTHON_EXECUTABLE="${PYTHON_CONFIG%-config}"
fi

if [ -z "$PYTHON_EXECUTABLE" ] || [ ! -x "$PYTHON_EXECUTABLE" ]; then
    # Fallback: try to find Python from the include dir
    PYTHON_INCLUDE=$(grep '_Python3_INCLUDE_DIR:INTERNAL' CMakeCache.txt 2>/dev/null | cut -d= -f2)
    if [[ "$PYTHON_INCLUDE" =~ python([0-9]+\.[0-9]+) ]]; then
        PYTHON_VERSION="${BASH_REMATCH[1]}"
        # Try common locations
        for candidate in /usr/local/Cellar/python@${PYTHON_VERSION}/*/bin/python${PYTHON_VERSION} \
                         /usr/local/bin/python${PYTHON_VERSION} \
                         /opt/homebrew/Cellar/python@${PYTHON_VERSION}/*/bin/python${PYTHON_VERSION} \
                         /opt/homebrew/bin/python${PYTHON_VERSION}; do
            # Use eval to expand the glob
            for real_candidate in $candidate; do
                if [ -x "$real_candidate" ]; then
                    PYTHON_EXECUTABLE="$real_candidate"
                    break 2
                fi
            done
        done
    fi
fi

# Final fallback to python3 in PATH
if [ -z "$PYTHON_EXECUTABLE" ] || [ ! -x "$PYTHON_EXECUTABLE" ]; then
    PYTHON_EXECUTABLE="python3"
fi

echo "Using Python: $PYTHON_EXECUTABLE"

cd "$SCRIPT_DIR"
source "$ROOT_DIR/add_vcfx_tools_to_path.sh"

PYTHONPATH="${BUILD_DIR}/python" "$PYTHON_EXECUTABLE" - <<'PY'
import vcfx
out = vcfx.trim("  hello  ")
if out != "hello":
    raise SystemExit('trim failed')
compressed = vcfx.read_maybe_compressed(b"hello")
if compressed != b"hello":
    raise SystemExit('read_maybe_compressed failed')
print('Python bindings OK:', out)
PY
