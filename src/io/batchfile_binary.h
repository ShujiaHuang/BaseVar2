/**
 * @file batchfile_binary.h
 *
 * @brief Binary batchfile format (.bbf / .bbi) for high-performance variant calling.
 *
 * Plan C: BGZF-compressed binary data file (.bbf) + custom binary index (.bbi).
 * Replaces the text-based batchfile to eliminate join()/split() overhead.
 *
 * .bbf format:
 *   Header: magic(4) + version(2) + sample_count(2) + sample_ids(variable)
 *   Records: per-position binary records with per-sample, per-read data
 *            (v2: each read stores its own ref_base for correct Indel normalization)
 *
 * .bbi format:
 *   Header:  magic(4) + version(2) + reserved(2) + entry_count(4)
 *   Entries: sorted array of (pos:4 + virtual_offset:8) pairs
 *   Footer:  magic(4)  -- integrity marker; if missing the file was truncated
 *
 * @author Shujia Huang
 * @date 2025
 */
#ifndef __INCLUDE_BATCHFILE_BINARY_H__
#define __INCLUDE_BATCHFILE_BINARY_H__

#include <string>
#include <vector>
#include <cstdint>

#include "iobgzf.h"
#include "caller_utils.h"  // AlignInfo, PosMap (includes robin_hood.h)
#include "basetype.h"
#include "io/utils.h"  // ngslib::GenomeRegion

// ---- Format constants ----
static const uint32_t BBF_MAGIC = 0x42424600;  // "BBF\0"
static const uint16_t BBF_VERSION = 2;  // v2: per-read ref_base (fix Indel normalization bug)

static const uint32_t BBI_MAGIC = 0x42424900;  // "BBI\0"
static const uint16_t BBI_VERSION = 1;
static const uint32_t BBI_FOOTER_MAGIC = 0x42494546;  // "BIEF" -- marks complete write

// ---- .bbf file header ----
struct BinaryBBFHeader {
    uint32_t magic;
    uint16_t version;
    uint16_t sample_count;
    // Followed by: sample_ids as comma-separated string (not null-terminated)
    // Length = file_size - 8 - 2 (header fixed part + sample_count)
};

// ---- .bbi index entry ----
struct BinaryIndexEntry {
    uint32_t pos;
    int64_t  virtual_offset;  // BGZF virtual offset (uoffset << 16 | offset)
};

// ---- .bbi file header ----
struct BinaryIndexHeader {
    uint32_t magic;
    uint16_t version;
    uint16_t reserved;       // padding for alignment
    uint32_t entry_count;
};

// =====================================================================
//  Write-side functions
// =====================================================================

/**
 * @brief Write binary header (magic + version + sample IDs) to .bbf file.
 */
void write_binary_header(ngslib::BGZFile &bf, const std::vector<std::string> &sample_ids);

/**
 * @brief Write one position record in binary format and record its virtual offset.
 *
 * @param bf              Output BGZF file handle
 * @param posinfomap_vec  Per-sample PosMap (one per sample in this batchfile)
 * @param gr              Current genomic region (chrom, start, end)
 * @param pos             Position to write (1-based)
 * @param out_index       Vector to append index entry (pos -> virtual_offset)
 */
void write_binary_record(ngslib::BGZFile &bf,
                         const PosMapVector &posinfomap_vec,
                         const ngslib::GenomeRegion &gr,
                         uint32_t pos,
                         std::vector<BinaryIndexEntry> &out_index);

/**
 * @brief Build and write .bbi index file from collected index entries.
 */
void write_binary_index(const std::string &index_path, const std::vector<BinaryIndexEntry> &entries);

// =====================================================================
//  Read-side functions
// =====================================================================

/**
 * @brief Load .bbi index into memory (sorted by position for binary search).
 *        Throws if the file is incomplete (footer missing).
 */
std::vector<BinaryIndexEntry> load_binary_index(const std::string &index_path);

/**
 * @brief Quick integrity check: verify .bbi has a valid footer without loading all entries.
 * @return true if the index file is complete and valid, false otherwise.
 */
bool validate_binary_index(const std::string &index_path);

/**
 * @brief Read sample IDs from .bbf binary header.
 */
std::vector<std::string> read_binary_sample_ids(const std::string &bbf_path);

/**
 * @brief Read one position record from .bbf at the current file position.
 *
 * @param bf          Input BGZF file handle (must be at the start of a record)
 * @param sample_count Number of samples in this batchfile
 * @param[out] ref_id Chromosome ID
 * @param[out] ref_pos Position (1-based)
 * @param[out] total_depth Total read depth across all samples
 * @param[out] samples_bi Per-sample BatchInfo (size = sample_count)
 * @return true if record was read successfully, false on EOF
 */
bool read_binary_record(ngslib::BGZFile &bf,
                        uint16_t sample_count,
                        std::string &ref_id,
                        uint32_t &ref_pos,
                        int &total_depth,
                        std::vector<BaseType::BatchInfo> &samples_bi);

#endif // __INCLUDE_BATCHFILE_BINARY_H__
