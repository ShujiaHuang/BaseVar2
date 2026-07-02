/**
 * @file bgzf_concat.h
 * @brief BGZF block-level copy engine for fast file concatenation.
 *
 * Provides functions to copy raw BGZF blocks between files without
 * decompression/recompression.  Used by:
 *   - concat.cpp  (naive VCF concatenation)
 *   - variant_caller.cpp  (merging sub-VCFs from multi-threaded calling)
 *
 * All functions are header-only (inline) for zero-overhead integration.
 *
 * @author Shujia Huang
 * @date 2025-07-02
 */
#ifndef __INCLUDE_NGSLIB_BGZF_CONCAT_H__
#define __INCLUDE_NGSLIB_BGZF_CONCAT_H__

#include "iobgzf.h"
#include "utils.h"      // ngslib::safe_remove

#include <vector>
#include <string>
#include <stdexcept>
#include <cstdint>

namespace ngslib {

/**
 * @brief Skip the VCF text header in a BGZF-compressed file.
 *
 * Reads BGZF blocks from @p in and optionally writes header lines to @p out.
 * Stops when the first non-'#' byte is found.  Returns the byte offset within
 * the current uncompressed block where data begins, or -1 on error.
 *
 * @param in            Source BGZFile, already opened; first block must be loaded.
 * @param out           Destination BGZFile (may be used to write header).
 * @param write_header  If true, write header lines to @p out.
 * @return Offset past the header in the current block, or -1 on error.
 */
inline int skip_vcf_gz_header(BGZFile &in, BGZFile &out, bool write_header) {
    const char *buffer = in.uncompressed_data();

    if (buffer[0] != '#') {
        std::cerr << "[bgzf_concat] Error: Could not parse VCF header, expected '#', found '"
                  << buffer[0] << "'\n";
        return -1;
    }

    int nskip = 1;
    while (true) {
        if (buffer[nskip] == '\n') {
            nskip++;
            if (nskip >= in.block_length()) {
                // Need to read the next BGZF block
                if (write_header) {
                    if (out.write_raw(buffer, in.block_length()) != in.block_length()) {
                        std::cerr << "[bgzf_concat] Error: failed to write header block\n";
                        return -1;
                    }
                }
                if (in.read_block() != 0) return -1;
                if (!in.block_length()) break;
                nskip = 0;
            }
            // Check if the next character starts a data line
            if (buffer[nskip] != '#') {
                if (write_header) {
                    if (out.write_raw(buffer, nskip) != nskip) {
                        std::cerr << "[bgzf_concat] Error: failed to write header tail\n";
                        return -1;
                    }
                }
                break;
            }
        }
        nskip++;
        if (nskip >= in.block_length()) {
            if (write_header) {
                if (out.write_raw(buffer, in.block_length()) != in.block_length()) {
                    std::cerr << "[bgzf_concat] Error: failed to write header block\n";
                    return -1;
                }
            }
            if (in.read_block() != 0) return -1;
            if (!in.block_length()) break;
            nskip = 0;
        }
    }
    return nskip;
}

/**
 * @brief Copy remaining raw BGZF blocks from an opened file to an output.
 *
 * Reads raw BGZF blocks from @p fin (from the current file position) and
 * writes them to @p fout.  EOF blocks (28-byte sentinels) are silently
 * skipped so that the output receives a single EOF block from @p fout's
 * destructor.
 *
 * Typical usage:
 *   1. Open @p fin, read/skip the VCF header region.
 *   2. Flush any pending decompressed writes to @p fout.
 *   3. Call bgzf_copy_blocks_to() to stream the remaining blocks.
 *
 * Both files must remain open; @p fin is NOT closed by this function.
 *
 * @param fin   Source BGZFile, already positioned past any header region.
 * @param fout  Destination BGZFile, opened for writing.
 * @param src_name  Source filename, used only in error messages.
 */
inline void bgzf_copy_blocks_to(BGZFile &fin, BGZFile &fout,
                                const std::string &src_name)
{
    const size_t page_size = BGZF_MAX_BLOCK_SIZE;
    std::vector<uint8_t> buf(page_size);

    while (true) {
        // Read the 18-byte BGZF block header
        ssize_t nread = fin.raw_block_read(buf.data(), BGZF_HEADER_SIZE);
        if (nread == 0) break;  // clean EOF
        if (nread != BGZF_HEADER_SIZE || bgzf_block_is_valid(buf.data()) != 0)
            throw std::runtime_error("[bgzf_copy] Error: malformed BGZF block in " + src_name);

        // Extract total block size (BSIZE + 1)
        uint32_t bsize = static_cast<uint32_t>(bgzf_block_bsize(buf.data())) + 1;
        if (bsize > page_size || bsize < static_cast<uint32_t>(BGZF_HEADER_SIZE))
            throw std::runtime_error(
                "[bgzf_copy] Error: malformed BGZF block size " + std::to_string(bsize)
                + " in " + src_name);

        // Read the rest of the block
        ssize_t nrest = fin.raw_block_read(buf.data() + BGZF_HEADER_SIZE,
                                           bsize - BGZF_HEADER_SIZE);
        if (nrest < 0)
            throw std::runtime_error("[bgzf_copy] Error: I/O error reading " + src_name);

        nread += nrest;
        if (nread != static_cast<ssize_t>(bsize))
            throw std::runtime_error(
                "[bgzf_copy] Error: truncated file " + src_name
                + " (expected " + std::to_string(bsize) + " bytes, got " + std::to_string(nread) + ")");

        // Skip intermediate BGZF EOF blocks (28 bytes)
        if (is_bgzf_eof(buf.data(), nread))
            continue;

        ssize_t nwr = fout.raw_block_write(buf.data(), nread);
        if (nwr != nread)
            throw std::runtime_error("[bgzf_copy] Error: write failed for " + src_name);
    }
}

/**
 * @brief Concatenate multiple BGZF files into an already-opened output.
 *
 * For each input file, opens it, copies all raw BGZF blocks to @p fout
 * (stripping intermediate EOF blocks), then optionally deletes the input.
 *
 * Does NOT handle VCF headers — the caller should write the header to @p fout
 * before calling this function.
 *
 * Output is always BGZF-compressed (block-level copy requires BGZF format).
 *
 * @param infiles        Vector of input BGZF-compressed file paths.
 * @param fout           Output BGZFile, already opened for writing.
 * @param remove_inputs  If true, delete each input file after copying.
 */
inline void bgzf_raw_concat(const std::vector<std::string> &infiles, BGZFile &fout,
                            bool remove_inputs = false)
{
    for (const auto &infile : infiles) {
        BGZFile fin(infile, "r");
        if (!fin.is_open())
            throw std::runtime_error("[bgzf_concat] Error: cannot open input file: " + infile);

        bgzf_copy_blocks_to(fin, fout, infile);
        // fin is closed here by ~BGZFile (RAII)

        if (remove_inputs)
            safe_remove(infile);
    }
}

/**
 * @brief Naive concatenation of VCF files with header-aware block-level copy.
 *
 * Unlike bgzf_raw_concat(), this function handles VCF headers:
 *   - First file: header is written to @p fout.
 *   - Subsequent files: header is skipped (only data blocks are copied).
 *
 * All input files must have compatible headers (same samples, same order).
 * Output is always BGZF-compressed.
 *
 * @param infiles        Vector of input BGZF-compressed VCF file paths.
 * @param fout           Output BGZFile, already opened for writing.
 * @param remove_inputs  If true, delete each input file after copying.
 */
inline void bgzf_naive_concat(const std::vector<std::string> &infiles, BGZFile &fout,
                              bool remove_inputs = false)
{
    for (size_t i = 0; i < infiles.size(); ++i) {
        BGZFile fin(infiles[i], "r");
        if (!fin.is_open())
            throw std::runtime_error("[bgzf_concat] Error: cannot open input file: " + infiles[i]);

        // Read the first BGZF block (contains the header)
        if (fin.read_block() != 0 || !fin.block_length())
            throw std::runtime_error("[bgzf_concat] Error: failed to read first block of " + infiles[i]);

        // Skip the VCF header region; write it only for the first file
        int nskip = skip_vcf_gz_header(fin, fout, (i == 0));
        if (nskip < 0)
            throw std::runtime_error("[bgzf_concat] Error: failed to parse header of " + infiles[i]);

        // Write any remaining data in the current uncompressed block (after header)
        if (fin.block_length() - nskip > 0) {
            if (fout.write_raw(fin.uncompressed_data() + nskip,
                               fin.block_length() - nskip) < 0)
                throw std::runtime_error("[bgzf_concat] Error: write failed for " + infiles[i]);
        }
        fout.flush();

        // Stream the rest of the file as raw BGZF blocks
        bgzf_copy_blocks_to(fin, fout, infiles[i]);
        // fin is closed here by ~BGZFile (RAII)

        if (remove_inputs)
            safe_remove(infiles[i]);
    }
}

} // namespace ngslib

#endif // __INCLUDE_NGSLIB_BGZF_CONCAT_H__
