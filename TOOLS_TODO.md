# VCFX Tools Optimization TODO

**Last Updated:** December 8, 2025

## Summary

- **23 tools optimized** with mmap + SIMD acceleration
- **8 tools need optimization** (identified as slow on 4GB file)
- **30 tools already fast** (<1 second on 4GB file)

---

## Optimized Tools (23 total)

| Tool | Optimization | Speedup | Status |
|------|-------------|---------|--------|
| **VCFX_validator** | mmap + SIMD | 1040x | ✅ Complete |
| **VCFX_variant_counter** | mmap + SIMD | 60x | ✅ Complete |
| **VCFX_fasta_converter** | mmap + SIMD + zero-copy + mmap temp | 50-100x | ✅ Complete |
| **VCFX_indel_normalizer** | mmap + SIMD | ~73x | ✅ Complete |
| **VCFX_missing_detector** | mmap + SIMD + MT pre-scan + zero-copy | ~42x | ✅ Complete |
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
| **VCFX_allele_freq_calc** | mmap + SIMD | ~20x | ✅ Complete |
| **VCFX_ld_calculator** | mmap + SIMD r² + MT matrix + distance pruning | 5-60x | ✅ Complete |

---

## Needs Optimization (8 tools, >60s on 4GB file)

These tools were identified as slow when processing large files:

| Tool | Status | Notes |
|------|--------|-------|
| **VCFX_af_subsetter** | Needs mmap | >60s on 4GB file |
| **VCFX_cross_sample_concordance** | Needs mmap | >60s on 4GB file |
| **VCFX_distance_calculator** | Needs mmap | >60s on 4GB file |
| **VCFX_dosage_calculator** | Needs mmap | >60s on 4GB file |
| **VCFX_duplicate_remover** | Needs mmap | >60s on 4GB file |
| **VCFX_metadata_summarizer** | Needs mmap | >60s on 4GB file |
| **VCFX_multiallelic_splitter** | Needs mmap | >60s on 4GB file |
| **VCFX_variant_classifier** | Needs mmap | >60s on 4GB file |

---

## Fast Tools (Already performant, <1 second on 4GB file)

These tools already perform well without optimization:

**Original fast tools:**
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

**Newly verified fast tools (on 4GB file):**
- VCFX_alignment_checker: 0.05s
- VCFX_ancestry_inferrer: 0.02s
- VCFX_annotation_extractor: 0.03s
- VCFX_compressor: 0.02s
- VCFX_field_extractor: 0.07s
- VCFX_format_converter: 0.02s
- VCFX_gl_filter: 0.02s
- VCFX_impact_filter: 0.02s
- VCFX_info_aggregator: 0.02s
- VCFX_info_parser: 0.02s
- VCFX_info_summarizer: 0.02s
- VCFX_population_filter: 0.06s
- VCFX_position_subsetter: 0.03s
- VCFX_record_filter: 0.02s
- VCFX_ref_comparator: 0.02s
- VCFX_region_subsampler: 0.02s
- VCFX_sample_extractor: 0.02s

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

// 3. SIMD line scanning (AVX2/SSE2/NEON)
static inline const char* findNewlineSIMD(const char* p, const char* end);

// 4. Zero-copy parsing with string_view
static inline std::string_view extractField(const char* line, int fieldIdx);

// 5. 4MB output buffer
class OutputBuffer { /* ... */ };

// 6. CLI: -i/--input FILE, -q/--quiet, -t/--threads (where applicable)
```
