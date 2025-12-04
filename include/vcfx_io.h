#ifndef VCFX_IO_H
#define VCFX_IO_H

/**
 * @file vcfx_io.h
 * @brief High-performance I/O utilities for VCFX tools
 *
 * This header provides optimized functions for VCF file processing:
 * - init_io(): Disable sync_with_stdio for faster I/O
 * - split_tabs(): Fast tab-delimited splitting with vector reuse
 * - split_tabs_view(): Zero-copy splitting using string_view
 * - split_char(): Generic delimiter splitting
 *
 * Usage:
 *   #include "vcfx_io.h"
 *
 *   int main() {
 *       vcfx::init_io();  // Call first in main()
 *       std::vector<std::string> fields;
 *       fields.reserve(16);
 *       while (std::getline(std::cin, line)) {
 *           vcfx::split_tabs(line, fields);  // Reuses capacity
 *       }
 *   }
 */

#include <iostream>
#include <string>
#include <string_view>
#include <vector>

namespace vcfx {

/**
 * @brief Initialize I/O for maximum performance
 *
 * Disables synchronization with C stdio and unties cin from cout.
 * Call this at the very start of main() before any I/O operations.
 * Expected speedup: 2-5x for I/O-bound tools.
 */
inline void init_io() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
}

/**
 * @brief Split a string by tabs into a reusable vector
 *
 * This function clears the output vector and reuses its capacity,
 * avoiding repeated allocations when called in a loop.
 *
 * @param line Input string to split
 * @param out Output vector (cleared but capacity preserved)
 * @param expected Expected number of fields for initial reserve
 * @return Number of fields found
 *
 * Performance: ~5-10x faster than std::stringstream approach
 */
inline size_t split_tabs(const std::string& line,
                         std::vector<std::string>& out,
                         size_t expected = 16) {
    out.clear();
    if (out.capacity() < expected) {
        out.reserve(expected);
    }

    size_t start = 0;
    size_t end;
    while ((end = line.find('\t', start)) != std::string::npos) {
        out.emplace_back(line, start, end - start);
        start = end + 1;
    }
    // Last field (after final tab or entire string if no tabs)
    out.emplace_back(line, start);
    return out.size();
}

/**
 * @brief Split a string by tabs, returning a new vector
 *
 * Convenience overload that creates and returns a new vector.
 * Use the two-argument version in loops for better performance.
 *
 * @param line Input string to split
 * @return Vector of fields
 */
inline std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> result;
    split_tabs(line, result);
    return result;
}

/**
 * @brief Zero-copy split using string_view (fastest option)
 *
 * Returns views into the original string - no copying or allocation
 * for the field data itself. The original string must remain valid
 * while using the views.
 *
 * @param line Input string_view to split
 * @param out Output vector of string_views (cleared but capacity preserved)
 * @param expected Expected number of fields for initial reserve
 * @return Number of fields found
 *
 * Performance: Fastest option for read-only field access
 */
inline size_t split_tabs_view(std::string_view line,
                              std::vector<std::string_view>& out,
                              size_t expected = 16) {
    out.clear();
    if (out.capacity() < expected) {
        out.reserve(expected);
    }

    size_t start = 0;
    size_t end;
    while ((end = line.find('\t', start)) != std::string_view::npos) {
        out.emplace_back(line.substr(start, end - start));
        start = end + 1;
    }
    out.emplace_back(line.substr(start));
    return out.size();
}

/**
 * @brief Generic delimiter split for INFO/FORMAT fields
 *
 * Splits by any single character delimiter using string_view.
 * Useful for parsing INFO (;), FORMAT (:), and other sub-fields.
 *
 * @param str Input string_view to split
 * @param delim Delimiter character
 * @param out Output vector of string_views
 * @param expected Expected number of fields for initial reserve
 * @return Number of fields found
 */
inline size_t split_char(std::string_view str, char delim,
                         std::vector<std::string_view>& out,
                         size_t expected = 8) {
    out.clear();
    if (out.capacity() < expected) {
        out.reserve(expected);
    }

    size_t start = 0;
    size_t end;
    while ((end = str.find(delim, start)) != std::string_view::npos) {
        out.emplace_back(str.substr(start, end - start));
        start = end + 1;
    }
    out.emplace_back(str.substr(start));
    return out.size();
}

/**
 * @brief Count tabs in a line (to estimate field count)
 *
 * @param line Input string_view
 * @return Number of fields (tabs + 1)
 */
inline size_t count_fields(std::string_view line) {
    size_t count = 1;
    for (char c : line) {
        if (c == '\t') ++count;
    }
    return count;
}

/**
 * @brief VCF standard field indices
 *
 * Use these constants instead of magic numbers for clarity.
 */
namespace VCF {
    constexpr int CHROM  = 0;
    constexpr int POS    = 1;
    constexpr int ID     = 2;
    constexpr int REF    = 3;
    constexpr int ALT    = 4;
    constexpr int QUAL   = 5;
    constexpr int FILTER = 6;
    constexpr int INFO   = 7;
    constexpr int FORMAT = 8;
    constexpr int FIRST_SAMPLE = 9;
    constexpr int MIN_FIELDS = 8;  // Minimum valid VCF data line
}

} // namespace vcfx

#endif // VCFX_IO_H
