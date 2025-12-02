add_executable(VCFX_validator VCFX_validator.cpp VCFX_validator.h)
target_link_libraries(VCFX_validator PRIVATE vcfx_core Threads::Threads)
