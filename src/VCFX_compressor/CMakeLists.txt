cmake_minimum_required(VERSION 3.14)
project(VCFX_compressor LANGUAGES CXX)

add_executable(VCFX_compressor
    VCFX_compressor.cpp
    VCFX_compressor.h
)

# Link against the core library + ZLIB
# When building for WebAssembly, zlib comes from the emscripten port (enabled
# via -sUSE_ZLIB=1 in the top-level CMakeLists), not from a find_package target.
target_link_libraries(VCFX_compressor PRIVATE vcfx_core)
if(NOT BUILD_WASM)
    target_link_libraries(VCFX_compressor PRIVATE ZLIB::ZLIB)
endif()
