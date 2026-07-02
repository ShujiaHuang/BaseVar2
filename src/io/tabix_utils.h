/**
 * @file tabix_utils.h
 * @brief Utility functions for building tabix indexes on BGZF-compressed files.
 *
 * This header provides a thin wrapper around htslib's `tbx_index_build()`
 * with automatic format detection and error handling.  It is intentionally
 * kept separate from `iobgzf.h` (which handles raw BGZF I/O) so that
 * format-specific tabix configuration stays out of the generic I/O layer.
 *
 * @author Shujia Huang
 * @date 2025-07-02
 */
#ifndef __INCLUDE_NGSLIB_TABIX_UTILS_H__
#define __INCLUDE_NGSLIB_TABIX_UTILS_H__

#include <string>
#include <stdexcept>

#include <htslib/tbx.h>

#include "utils.h"  // ngslib::suffix_name

namespace ngslib {

    /// Standard VCF tabix configuration:
    ///   seq_col=1 (CHROM), beg_col=1 (POS), end_col=2 (END),
    ///   header_char='#', skip_lines=0, preset=TBX_VCF(1)
    inline constexpr tbx_conf_t TBX_CONF_VCF = {1, 1, 2, 0, '#', 0};

    /**
     * @brief Build a tabix index (.tbi) for a BGZF-compressed file.
     *
     * If @p filename does not end with ".gz", the function returns
     * immediately (index is only meaningful for BGZF output).
     *
     * @param filename  Path to the BGZF-compressed file.
     * @param conf      Tabix column configuration (e.g. TBX_CONF_VCF).
     * @return 0 on success.
     * @throw std::runtime_error if tbx_index_build fails.
     */
    inline int build_tabix_index(const std::string &filename,
                                 const tbx_conf_t &conf = TBX_CONF_VCF)
    {
        if (suffix_name(filename) != ".gz")
            return 0;  // Not BGZF — skip indexing silently

        if (tbx_index_build(filename.c_str(), 0, &conf))
            throw std::runtime_error(
                "[build_tabix_index] tbx_index_build failed. "
                "Is the file bgzip-compressed? File: " + filename);

        return 0;
    }

}  // namespace ngslib

#endif  // __INCLUDE_NGSLIB_TABIX_UTILS_H__
