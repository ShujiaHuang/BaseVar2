/**
 * @file batchfile_binary.cpp
 *
 * @brief Implementation of binary batchfile format (.bbf / .bbi)
 *
 * Plan C: BGZF-compressed binary data + custom binary index.
 * Eliminates join()/split() overhead in the batchfile read/write pipeline.
 *
 * @author Shujia Huang
 * @date 2025
 */
#include "batchfile_binary.h"

#include <cstring>   // memcpy
#include <algorithm> // lower_bound
#include <stdexcept>

// =====================================================================
//  Write-side functions
// =====================================================================

void write_binary_header(ngslib::BGZFile &bf, const std::vector<std::string> &sample_ids) {
    uint32_t magic = BBF_MAGIC;
    uint16_t version = BBF_VERSION;
    uint16_t sample_count = static_cast<uint16_t>(sample_ids.size());

    bf.write_raw(&magic,         sizeof(magic));
    bf.write_raw(&version,       sizeof(version));
    bf.write_raw(&sample_count,  sizeof(sample_count));

    // Sample IDs as comma-separated string
    std::string ids = ngslib::join(sample_ids, ",");
    uint32_t ids_len = static_cast<uint32_t>(ids.size());
    bf.write_raw(&ids_len, sizeof(ids_len));       // store length for reader
    bf.write_raw(ids.data(), ids.size());
}

void write_binary_record(ngslib::BGZFile &bf,
                         const PosMapVector &posinfomap_vec,
                         const ngslib::GenomeRegion &gr,
                         uint32_t pos,
                         std::vector<BinaryIndexEntry> &out_index)
{
    // Record virtual offset BEFORE writing (this is what the index points to)
    int64_t voffset = bf.tell_virtual();

    size_t sn = posinfomap_vec.size();

    // ---- Position-level fields ----
    uint8_t chrom_len = static_cast<uint8_t>(gr.chrom.size());
    bf.write_raw(&chrom_len, sizeof(chrom_len));
    bf.write_raw(gr.chrom.data(), chrom_len);
    bf.write_raw(&pos, sizeof(pos));

    // Determine ref_base from the first sample that has data at this position.
    // After left-alignment, ref_base is identical across all reads at the same position.
    std::string ref_base;
    bool ref_base_found = false;
    for (size_t i = 0; i < sn && !ref_base_found; ++i) {
        auto it = posinfomap_vec[i].find(pos);
        if (it != posinfomap_vec[i].end() && !it->second.align_bases.empty()) {
            ref_base = it->second.align_bases[0].ref_base;
            ref_base_found = true;
        }
    }

    uint16_t ref_base_len = static_cast<uint16_t>(ref_base.size());
    bf.write_raw(&ref_base_len, sizeof(ref_base_len));
    if (ref_base_len > 0) {
        bf.write_raw(ref_base.data(), ref_base_len);
    }

    // ---- Per-sample data ----
    for (size_t i = 0; i < sn; ++i) {
        auto it = posinfomap_vec[i].find(pos);
        if (it != posinfomap_vec[i].end() && it->second.ref_id == gr.chrom && it->second.ref_pos == pos) {
            uint16_t depth = static_cast<uint16_t>(it->second.align_bases.size());
            bf.write_raw(&depth, sizeof(depth));

            for (const auto &ab : it->second.align_bases) {
                // read_base (variable length, includes +/- prefix for INS/DEL)
                uint16_t rb_len = static_cast<uint16_t>(ab.read_base.size());
                bf.write_raw(&rb_len, sizeof(rb_len));
                if (rb_len > 0) {
                    bf.write_raw(ab.read_base.data(), rb_len);
                }

                // qual
                bf.write_raw(&ab.base_qual, sizeof(ab.base_qual));

                // rpr (uint32, supports long reads up to ~4 billion bp)
                uint32_t rpr = static_cast<uint32_t>(ab.rpr);
                bf.write_raw(&rpr, sizeof(rpr));

                // mapq
                uint8_t mapq = static_cast<uint8_t>(ab.mapq);
                bf.write_raw(&mapq, sizeof(mapq));

                // strand: 0='+', 1='-', 2='*'
                uint8_t strand;
                if      (ab.map_strand == '+') strand = 0;
                else if (ab.map_strand == '-') strand = 1;
                else                           strand = 2;  // '*' or unknown
                bf.write_raw(&strand, sizeof(strand));
            }
        } else {
            // No data at this position for this sample
            uint16_t depth = 0;
            bf.write_raw(&depth, sizeof(depth));
        }
    }

    // Record index entry
    out_index.push_back({pos, voffset});
}

void write_binary_index(const std::string &index_path, const std::vector<BinaryIndexEntry> &entries) {
    FILE *fp = fopen(index_path.c_str(), "wb");
    if (!fp) {
        throw std::runtime_error("Cannot create binary index file: " + index_path);
    }

    // Header
    uint32_t magic   = BBI_MAGIC;
    uint16_t version = BBI_VERSION;
    uint16_t reserved = 0;
    uint32_t count   = static_cast<uint32_t>(entries.size());

    fwrite(&magic,    sizeof(magic),    1, fp);
    fwrite(&version,  sizeof(version),  1, fp);
    fwrite(&reserved, sizeof(reserved), 1, fp);
    fwrite(&count,    sizeof(count),    1, fp);

    // Entries (already sorted by position since we write positions in order)
    for (const auto &e : entries) {
        fwrite(&e.pos,            sizeof(e.pos),            1, fp);
        fwrite(&e.virtual_offset, sizeof(e.virtual_offset), 1, fp);
    }

    // Footer: integrity marker — only written after all entries are flushed
    uint32_t footer = BBI_FOOTER_MAGIC;
    fwrite(&footer, sizeof(footer), 1, fp);

    fclose(fp);
}

// =====================================================================
//  Read-side functions
// =====================================================================

std::vector<BinaryIndexEntry> load_binary_index(const std::string &index_path) {
    FILE *fp = fopen(index_path.c_str(), "rb");
    if (!fp) {
        throw std::runtime_error("Cannot open binary index file: " + index_path);
    }

    // Read and validate header
    BinaryIndexHeader hdr;
    if (fread(&hdr, sizeof(hdr), 1, fp) != 1) {
        fclose(fp);
        throw std::runtime_error("Cannot read binary index header: " + index_path);
    }
    if (hdr.magic != BBI_MAGIC) {
        fclose(fp);
        throw std::runtime_error("Invalid binary index magic: " + index_path);
    }
    if (hdr.version != BBI_VERSION) {
        fclose(fp);
        throw std::runtime_error("Unsupported binary index version " +
                                 std::to_string(hdr.version) + " in: " + index_path);
    }

    // Read all entries
    std::vector<BinaryIndexEntry> entries(hdr.entry_count);
    for (uint32_t i = 0; i < hdr.entry_count; ++i) {
        if (fread(&entries[i].pos,            sizeof(uint32_t), 1, fp) != 1 ||
            fread(&entries[i].virtual_offset, sizeof(int64_t),  1, fp) != 1) {
            fclose(fp);
            throw std::runtime_error("Truncated binary index at entry " +
                                     std::to_string(i) + " in: " + index_path);
        }
    }

    // Verify footer integrity marker
    uint32_t footer = 0;
    if (fread(&footer, sizeof(footer), 1, fp) != 1 || footer != BBI_FOOTER_MAGIC) {
        fclose(fp);
        throw std::runtime_error("Binary index file is incomplete (missing footer): " + index_path);
    }

    fclose(fp);
    return entries;  // already sorted by pos
}

bool validate_binary_index(const std::string &index_path) {
    FILE *fp = fopen(index_path.c_str(), "rb");
    if (!fp) return false;

    // Seek to the last 4 bytes (footer)
    if (fseek(fp, -4, SEEK_END) != 0) {
        fclose(fp);
        return false;
    }

    uint32_t footer = 0;
    bool ok = (fread(&footer, sizeof(footer), 1, fp) == 1) && (footer == BBI_FOOTER_MAGIC);
    fclose(fp);
    return ok;
}

std::vector<std::string> read_binary_sample_ids(const std::string &bbf_path) {
    ngslib::BGZFile f(bbf_path.c_str(), "r");
    if (!f.is_open()) {
        throw std::runtime_error("Cannot open binary batchfile: " + bbf_path);
    }

    // Read fixed header fields
    uint32_t magic;
    uint16_t version, sample_count;
    f.read_raw(&magic,        sizeof(magic));
    f.read_raw(&version,      sizeof(version));
    f.read_raw(&sample_count, sizeof(sample_count));

    if (magic != BBF_MAGIC) {
        throw std::runtime_error("Invalid binary batchfile magic: " + bbf_path);
    }
    if (version != BBF_VERSION) {
        throw std::runtime_error("Unsupported binary batchfile version " +
                                 std::to_string(version) + " in: " + bbf_path);
    }

    // Read sample IDs string length, then the string itself
    uint32_t ids_len;
    f.read_raw(&ids_len, sizeof(ids_len));

    std::string ids_str(ids_len, '\0');
    if (ids_len > 0) {
        f.read_raw(&ids_str[0], ids_len);
    }

    std::vector<std::string> sample_ids;
    ngslib::split(ids_str, sample_ids, ",", true);
    return sample_ids;
}

bool read_binary_record(ngslib::BGZFile &bf,
                        uint16_t sample_count,
                        std::string &ref_id,
                        uint32_t &ref_pos,
                        int &total_depth,
                        std::vector<BaseType::BatchInfo> &samples_bi)
{
    // ---- Position-level fields ----
    uint8_t chrom_len;
    if (bf.read_raw(&chrom_len, sizeof(chrom_len)) != 1) {
        return false;  // EOF
    }

    ref_id.resize(chrom_len);
    if (chrom_len > 0) {
        bf.read_raw(&ref_id[0], chrom_len);
    }

    bf.read_raw(&ref_pos, sizeof(ref_pos));

    uint16_t ref_base_len;
    bf.read_raw(&ref_base_len, sizeof(ref_base_len));
    std::string ref_base(ref_base_len, '\0');
    if (ref_base_len > 0) {
        bf.read_raw(&ref_base[0], ref_base_len);
    }

    // ---- Per-sample data ----
    total_depth = 0;
    samples_bi.clear();
    samples_bi.reserve(sample_count);

    for (uint16_t i = 0; i < sample_count; ++i) {
        uint16_t depth;
        bf.read_raw(&depth, sizeof(depth));
        total_depth += depth;

        BaseType::BatchInfo smp_bi(ref_id, ref_pos);

        if (depth > 0) {
            smp_bi.ref_bases.reserve(depth);
            smp_bi.align_bases.reserve(depth);
            smp_bi.align_base_quals.reserve(depth);
            smp_bi.base_pos_ranks.reserve(depth);
            smp_bi.mapqs.reserve(depth);
            smp_bi.map_strands.reserve(depth);

            for (uint16_t j = 0; j < depth; ++j) {
                // read_base (variable length)
                uint16_t rb_len;
                bf.read_raw(&rb_len, sizeof(rb_len));
                std::string read_base(rb_len, '\0');
                if (rb_len > 0) {
                    bf.read_raw(&read_base[0], rb_len);
                }

                // qual
                char qual;
                bf.read_raw(&qual, sizeof(qual));

                // rpr (uint32)
                uint32_t rpr;
                bf.read_raw(&rpr, sizeof(rpr));

                // mapq
                uint8_t mapq;
                bf.read_raw(&mapq, sizeof(mapq));

                // strand
                uint8_t strand;
                bf.read_raw(&strand, sizeof(strand));
                char map_strand;
                if      (strand == 0) map_strand = '+';
                else if (strand == 1) map_strand = '-';
                else                  map_strand = '*';

                // Populate BatchInfo vectors
                // ref_base is stored once per position; replicate for each read
                smp_bi.ref_bases.push_back(ref_base);
                smp_bi.align_bases.push_back(std::move(read_base));
                smp_bi.align_base_quals.push_back(qual);
                smp_bi.base_pos_ranks.push_back(static_cast<int>(rpr));
                smp_bi.mapqs.push_back(static_cast<int>(mapq));
                smp_bi.map_strands.push_back(map_strand);
            }
        } else {
            // No data: replicate text format behavior
            // Text format stores "N" for ref_bases and align_bases when depth=0
            smp_bi.ref_bases.push_back("N");
            smp_bi.align_bases.push_back("N");
            smp_bi.align_base_quals.push_back('!');  // BASE_Q0_CHAR = Phred+33 offset
            smp_bi.base_pos_ranks.push_back(0);
            smp_bi.mapqs.push_back(0);
            smp_bi.map_strands.push_back('*');
        }

        samples_bi.push_back(std::move(smp_bi));
    }

    return true;
}
