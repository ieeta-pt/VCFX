cmake_minimum_required(VERSION 3.14)

# Set version information
set(VCFX_VERSION_MAJOR 1)
set(VCFX_VERSION_MINOR 0)
set(VCFX_VERSION_PATCH 4)
set(VCFX_VERSION "${VCFX_VERSION_MAJOR}.${VCFX_VERSION_MINOR}.${VCFX_VERSION_PATCH}")

add_compile_definitions(VCFX_VERSION="${VCFX_VERSION}")

project(VCFX
  VERSION ${VCFX_VERSION}
  DESCRIPTION "A Comprehensive VCF Manipulation Toolkit"
  LANGUAGES CXX
)

# Set a user-friendly default install prefix when none is provided.
# This prevents installation attempts into system directories when
# running without root privileges.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local" CACHE PATH "Install prefix" FORCE)
endif()

include(GNUInstallDirs)

# Optionally allow building for WebAssembly via Emscripten
option(BUILD_WASM "Build with emscripten toolchain" OFF)
option(PYTHON_BINDINGS "Build Python bindings" ON)

if(BUILD_WASM)
    if(NOT CMAKE_TOOLCHAIN_FILE)
        if(DEFINED ENV{EMSDK} AND EXISTS "$ENV{EMSDK}/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake")
            set(CMAKE_TOOLCHAIN_FILE "$ENV{EMSDK}/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake" CACHE FILEPATH "Emscripten toolchain" FORCE)
        elseif(DEFINED ENV{EMSCRIPTEN} AND EXISTS "$ENV{EMSCRIPTEN}/cmake/Modules/Platform/Emscripten.cmake")
            set(CMAKE_TOOLCHAIN_FILE "$ENV{EMSCRIPTEN}/cmake/Modules/Platform/Emscripten.cmake" CACHE FILEPATH "Emscripten toolchain" FORCE)
        endif()
    endif()

    if(NOT EXISTS "${CMAKE_TOOLCHAIN_FILE}")
        message(FATAL_ERROR "Emscripten toolchain file not found. Please set CMAKE_TOOLCHAIN_FILE or EMSDK.")
    endif()

    message(STATUS "Building for WebAssembly (Emscripten).")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

# Enable Link-Time Optimization for Release builds (significant perf improvement)
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    include(CheckIPOSupported)
    check_ipo_supported(RESULT ipo_supported OUTPUT ipo_output)
    if(ipo_supported)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
        message(STATUS "LTO enabled for Release build")
    else()
        message(STATUS "LTO not supported: ${ipo_output}")
    endif()
endif()

# Find threading library
find_package(Threads REQUIRED)

# Locate zlib if we need it (some tools do)
# On macOS with Homebrew, you might need to specify paths. E.g.:
# list(APPEND CMAKE_PREFIX_PATH "/usr/local/opt/zlib" "/opt/homebrew/opt/zlib")
find_package(ZLIB REQUIRED)

# If building natively on ARM, e.g. Apple M1
if((NOT BUILD_WASM) AND (CMAKE_SYSTEM_PROCESSOR MATCHES "arm|aarch64"))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcpu=native -O3")
endif()

enable_testing()

# Add top-level 'src' subdirectory, which in turn references each tool subdirectory
add_subdirectory(src)

if(PYTHON_BINDINGS)
  add_subdirectory(python)
endif()

# Add the test suite
add_subdirectory(tests)

# Installation configuration

# Install header files
install(FILES 
  include/vcfx_core.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/vcfx
)

# Install documentation files if they exist
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
  install(FILES 
    README.md
    DESTINATION ${CMAKE_INSTALL_DOCDIR}
  )
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
  install(FILES 
    LICENSE
    DESTINATION ${CMAKE_INSTALL_DOCDIR}
  )
endif()

message(STATUS "VCFX version: ${VCFX_VERSION}")
message(STATUS "Installation prefix: ${CMAKE_INSTALL_PREFIX}")
