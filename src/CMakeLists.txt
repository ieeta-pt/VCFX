cmake_minimum_required(VERSION 3.14)

# Build the core library from your shared code
add_library(vcfx_core STATIC vcfx_core.cpp)
set_property(TARGET vcfx_core PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(vcfx_core PUBLIC ${CMAKE_CURRENT_LIST_DIR}/../include)
target_link_libraries(vcfx_core PUBLIC Threads::Threads)
if(NOT BUILD_WASM)
    target_link_libraries(vcfx_core PUBLIC ZLIB::ZLIB)
endif()
if(WIN32)
    target_link_libraries(vcfx_core PUBLIC getopt)
endif()

# Add all tool subdirectories
add_subdirectory(vcfx_wrapper)

# Automatically detect tool directories named "VCFX_*" and
# build/install them.
set(VCFX_TOOLS vcfx)
file(GLOB TOOL_DIRS RELATIVE ${CMAKE_CURRENT_LIST_DIR} VCFX_*)
list(REMOVE_ITEM TOOL_DIRS vcfx_wrapper) # avoid duplicate on case-insensitive filesystems
foreach(dir ${TOOL_DIRS})
    if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${dir}")
        add_subdirectory(${dir})
        list(APPEND VCFX_TOOLS ${dir})
    endif()
endforeach()

# Install the core library
install(TARGETS vcfx_core
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# Install all tool executables
install(TARGETS ${VCFX_TOOLS}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

message(STATUS "VCFX tools to be installed: ${VCFX_TOOLS}")
