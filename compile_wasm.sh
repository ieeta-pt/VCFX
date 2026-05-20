#!/usr/bin/env bash
set -e

mkdir -p build_wasm
cd build_wasm

if command -v emcmake >/dev/null 2>&1; then
    emcmake cmake -DBUILD_WASM=ON -DPYTHON_BINDINGS=OFF ..
else
    cmake -DBUILD_WASM=ON -DPYTHON_BINDINGS=OFF ..
fi

cmake --build .

echo "All VCFX tools and the vcfx wrapper built for WebAssembly in build_wasm/."
echo "Use 'ls -R build_wasm' to see output. If you want .html or .js from Emscripten, you can adjust linking flags or suffixes."
