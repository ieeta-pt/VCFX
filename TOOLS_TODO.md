# VCFX Tools Optimization TODO

**Last Updated:** December 8, 2025

## Already Optimized (18 tools with mmap support)

| Tool | Optimization | Speedup | Status |
|------|-------------|---------|--------|
| **VCFX_validator** | mmap + SIMD | 1040x | ✅ Complete |
| **VCFX_variant_counter** | mmap + SIMD | 60x | ✅ Complete |
| **VCFX_fasta_converter** | mmap + SIMD + zero-copy + mmap temp | 50-100x | ✅ Complete |
| **VCFX_indexer** | mmap + SIMD | 32x | ✅ Complete |
| **VCFX_sorter** | mmap + precomputed IDs | 40x | ✅ Complete |
| **VCFX_phred_filter** | mmap | 26x | ✅ Complete |
| **VCFX_missing_data_handler** | mmap | 50x | ✅ Complete |
| **VCFX_nonref_filter** | mmap + SIMD | ~50x | ✅ Complete |
| **VCFX_genotype_query** | mmap + FORMAT caching | ~30x | ✅ Complete |
| **VCFX_phase_checker** | mmap | ~30x | ✅ Complete |
| **VCFX_haplotype_phaser** | mmap + SIMD + zero-copy | 16x | ✅ Complete |
| **VCFX_haplotype_extractor** | mmap + SIMD + zero-copy | 3.3x+ | ✅ Complete |
| **VCFX_allele_counter** | mmap + SIMD + MT + batch | 8-10x | ✅ Complete |
| **VCFX_diff_tool** | mmap + SIMD | ~20x | ✅ Complete |
| **VCFX_concordance_checker** | mmap + SIMD | ~20x | ✅ Complete |
| **VCFX_allele_balance_calc** | mmap + SIMD + incremental flush | ~50x | ✅ Complete |
| **VCFX_inbreeding_calculator** | mmap + SIMD | ~21x | ✅ Complete |
| **VCFX_hwe_tester** | mmap + SIMD | ~18x | ✅ Complete |

---

## Still Need Optimization

### Priority 1: Very Slow Tools (>5 min on 4GB file)

| Tool | Current Time | Complexity | Optimization Strategy |
|------|-------------|------------|----------------------|
| **VCFX_ld_calculator** | 32 min | O(variants²) | Already slow by design, consider window limiting |

### Priority 3: Moderate Speed Tools (3-5 min)

| Tool | Current Time | Notes |
|------|-------------|-------|
| **VCFX_allele_freq_calc** | 4.9 min | mmap would help |
| **VCFX_indel_normalizer** | 4.9 min | mmap would help |
| **VCFX_missing_detector** | 5 min | mmap would help |

---

## Fast Tools (No optimization needed, <1 second)

These tools already perform well:
- VCFX_reformatter: 0.15s
- VCFX_header_parser: 0.17s
- VCFX_merger: 0.17s
- VCFX_quality_adjuster: 0.18s
- VCFX_sv_handler: 0.19s
- VCFX_probability_filter: 0.21s
- VCFX_phase_quality_filter: 0.21s
- VCFX_outlier_detector: 0.21s
- VCFX_custom_annotator: 0.21s
- VCFX_file_splitter: 0.21s
- VCFX_subsampler: 0.22s
- VCFX_allele_balance_filter: 0.26s
- VCFX_ancestry_assigner: 0.28s

---

## Optimization Pattern (Proven)

Apply this pattern for ~30-100x speedup:

```cpp
// 1. Memory-mapped I/O
#include <sys/mman.h>
struct MappedFile {
    const char *data = nullptr;
    size_t size = 0;
    int fd = -1;
    bool open(const char *path);
    void close();
};

// 2. MADV_SEQUENTIAL hint
madvise((void*)data, size, MADV_SEQUENTIAL | MADV_WILLNEED);

// 3. SIMD line scanning (AVX2/SSE2/memchr fallback)
static inline const char* findNewlineSIMD(const char* p, const char* end);

// 4. Zero-copy parsing with string_view
static inline std::string_view extractField(const char* line, int fieldIdx);

// 5. 1MB output buffer
class OutputBuffer { /* ... */ };

// 6. CLI: -i/--input FILE, -q/--quiet
```

---

## Recommended Next Optimization

**VCFX_allele_freq_calc** - 4.9 minutes:
- Similar pattern to inbreeding_calculator
- Expected improvement: 15-20x

**VCFX_indel_normalizer** - 4.9 minutes:
- Similar pattern to other tools
- Expected improvement: 15-20x
