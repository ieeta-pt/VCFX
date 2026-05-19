# Top-level make wrapper for projects/tooling that expect a Makefile in the
# repository root. The actual build system is CMake; this just shells out to
# cmake + make and flattens the resulting binaries into bin/ for consumers
# that want a single output directory.
#
# Typical usage:
#   make           # build all tools, flatten binaries into bin/
#   make clean     # remove build/ and bin/
#
# This wrapper coexists with compile_wasm.sh and the CMake/Bioconda/Docker
# flows; none of them depend on this Makefile.

BUILD_DIR := build
OUT_DIR   := bin
CMAKE     ?= cmake

all:
	mkdir -p $(BUILD_DIR) $(OUT_DIR)
	cd $(BUILD_DIR) && $(CMAKE) -DPYTHON_BINDINGS=OFF .. && $(MAKE)
	find $(BUILD_DIR)/src -maxdepth 2 -type f -perm -u+x ! -path "*/CMakeFiles/*" -exec cp {} $(OUT_DIR)/ \;

clean:
	rm -rf $(BUILD_DIR) $(OUT_DIR)

.PHONY: all clean
