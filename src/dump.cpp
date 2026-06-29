/**
 * @file dump.cpp
 *
 * @brief Implementation of `basevar dump` subcommand for inspecting
 *        binary batchfile (.bbf) and binary index (.bbi) files.
 *
 * @author Shujia Huang
 * @date 2026
 */
#include "dump.h"

#include <iostream>
#include <string>
#include <vector>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <getopt.h>

#include "io/batchfile_binary.h"
#include "io/iobgzf.h"

namespace basevar {
namespace dump {

// =====================================================================
//  Usage
// =====================================================================

static void dump_usage() {
    std::cout <<
        "Usage: basevar dump <file> [options]\n"
        "\n"
        "Inspect binary batchfile (.bbf) or binary index (.bbi) files.\n"
        "\n"
        "Options:\n"
        "  --header          Show only file header (sample IDs) and exit (.bbf only)\n"
        "  -r, --region      Dump records in region chr:start-end (.bbf only)\n"
        "  -n, --limit INT   Dump at most INT records (.bbf only, default: all)\n"
        "  -v, --verbose     Show per-sample details for each record (.bbf only)\n"
        "  --entries         Show all index entries (.bbi only)\n"
        "  -h, --help        Show this help message\n"
        "\n"
        "Examples:\n"
        "  basevar dump sample.bbf                 # header + record summary\n"
        "  basevar dump sample.bbf --header        # header only\n"
        "  basevar dump sample.bbf -r chr1:1000-2000 -v\n"
        "  basevar dump sample.bbf.bbi             # index summary\n"
        "  basevar dump sample.bbf.bbi --entries   # all index entries\n"
        << std::endl;
}

// =====================================================================
//  .bbi index dump
// =====================================================================

static int dump_bbi(const std::string &path, bool show_entries) {
    FILE *fp = fopen(path.c_str(), "rb");
    if (!fp) {
        std::cerr << "[ERROR] Cannot open file: " << path << std::endl;
        return 1;
    }

    // Read header
    BinaryIndexHeader hdr;
    if (fread(&hdr, sizeof(hdr), 1, fp) != 1) {
        std::cerr << "[ERROR] Cannot read header (file too small): " << path << std::endl;
        fclose(fp);
        return 1;
    }

    // Validate magic
    bool magic_ok = (hdr.magic == BBI_MAGIC);
    bool version_ok = (hdr.version == BBI_VERSION);

    std::cout << "BBI Index: " << path << "\n";
    std::cout << "  Magic:      0x" << std::hex << hdr.magic << std::dec
              << (magic_ok ? " (valid)" : " INVALID!") << "\n";
    std::cout << "  Version:    " << hdr.version
              << (version_ok ? " (valid)" : " UNSUPPORTED!") << "\n";
    std::cout << "  Reserved:   " << hdr.reserved << "\n";
    std::cout << "  Entries:    " << hdr.entry_count << "\n";

    if (!magic_ok) {
        std::cerr << "[ERROR] Invalid BBI magic number. Not a .bbi file: " << path << std::endl;
        fclose(fp);
        return 1;
    }

    // Read entries to get position range
    uint32_t min_pos = 0, max_pos = 0;
    int64_t min_voffset = 0, max_voffset = 0;

    if (hdr.entry_count > 0) {
        // Read first entry
        BinaryIndexEntry first_entry;
        if (fread(&first_entry.pos, sizeof(uint32_t), 1, fp) != 1 ||
            fread(&first_entry.virtual_offset, sizeof(int64_t), 1, fp) != 1) {
            std::cerr << "[ERROR] Truncated index (cannot read first entry): " << path << std::endl;
            fclose(fp);
            return 1;
        }
        min_pos = first_entry.pos;
        min_voffset = first_entry.virtual_offset;

        if (hdr.entry_count > 1) {
            // Seek to last entry: position = header_size + (count-1) * entry_size
            long last_entry_offset = (long)(sizeof(BinaryIndexHeader) +
                                            (uint64_t)(hdr.entry_count - 1) * (sizeof(uint32_t) + sizeof(int64_t)));
            fseek(fp, last_entry_offset, SEEK_SET);
            BinaryIndexEntry last_entry;
            if (fread(&last_entry.pos, sizeof(uint32_t), 1, fp) != 1 ||
                fread(&last_entry.virtual_offset, sizeof(int64_t), 1, fp) != 1) {
                std::cerr << "[ERROR] Truncated index (cannot read last entry): " << path << std::endl;
                fclose(fp);
                return 1;
            }
            max_pos = last_entry.pos;
            max_voffset = last_entry.virtual_offset;
        } else {
            max_pos = min_pos;
            max_voffset = min_voffset;
        }

        std::cout << "  Position range: " << min_pos << " - " << max_pos << "\n";
        std::cout << "  Voffset range:  0x" << std::hex << min_voffset
                  << " - 0x" << max_voffset << std::dec << "\n";
    }

    // Check footer
    if (fseek(fp, -4, SEEK_END) == 0) {
        uint32_t footer = 0;
        bool footer_ok = (fread(&footer, sizeof(footer), 1, fp) == 1) &&
                         (footer == BBI_FOOTER_MAGIC);
        std::cout << "  Footer:     0x" << std::hex << footer << std::dec
                  << (footer_ok ? " (valid)" : " MISSING/INVALID - file may be truncated!") << "\n";
    }

    // Show entries if requested
    if (show_entries && hdr.entry_count > 0) {
        std::cout << "\n  Position       Virtual Offset\n";
        std::cout << "  -----------    ----------------\n";

        // Re-read from the beginning of entries
        // Header is 12 bytes, entries start right after
        fseek(fp, sizeof(BinaryIndexHeader), SEEK_SET);

        for (uint32_t i = 0; i < hdr.entry_count; ++i) {
            uint32_t pos;
            int64_t voffset;
            if (fread(&pos, sizeof(uint32_t), 1, fp) != 1 ||
                fread(&voffset, sizeof(int64_t), 1, fp) != 1) {
                std::cerr << "[ERROR] Truncated at entry " << i << std::endl;
                break;
            }
            printf("  %-13u  0x%016lx\n", pos, (unsigned long)voffset);
        }
    }

    fclose(fp);
    return 0;
}

// =====================================================================
//  .bbf data file dump
// =====================================================================

static int dump_bbf(const std::string &path,
                    bool header_only,
                    const std::string &region_str,
                    int limit,
                    bool verbose)
{
    // Open BGZF file
    ngslib::BGZFile bf(path.c_str(), "r");
    if (!bf.is_open()) {
        std::cerr << "[ERROR] Cannot open file: " << path << std::endl;
        return 1;
    }

    // Read header
    uint32_t magic;
    uint16_t version;
    uint16_t sample_count;
    bf.read_raw(&magic, sizeof(magic));
    bf.read_raw(&version, sizeof(version));
    bf.read_raw(&sample_count, sizeof(sample_count));

    bool magic_ok = (magic == BBF_MAGIC);
    bool version_ok = (version == BBF_VERSION);

    uint32_t ids_len;
    bf.read_raw(&ids_len, sizeof(ids_len));
    std::string ids_str(ids_len, '\0');
    if (ids_len > 0) {
        bf.read_raw(&ids_str[0], ids_len);
    }

    // Parse sample IDs
    std::vector<std::string> sample_ids;
    ngslib::split(ids_str, sample_ids, ",", true);

    // Print header
    std::cout << "BBF Data: " << path << "\n";
    std::cout << "  Magic:         0x" << std::hex << magic << std::dec
              << (magic_ok ? " (valid)" : " INVALID!") << "\n";
    std::cout << "  Version:       " << version
              << (version_ok ? " (valid)" : " UNSUPPORTED!") << "\n";
    std::cout << "  Sample count:  " << sample_count << "\n";
    std::cout << "  Sample IDs:    " << ids_str << "\n";

    if (!magic_ok) {
        std::cerr << "[ERROR] Invalid BBF magic number. Not a .bbf file: " << path << std::endl;
        return 1;
    }

    if (header_only) return 0;

    // Parse region if provided
    std::string query_chrom;
    uint32_t query_start = 0, query_end = UINT32_MAX;
    bool has_region = false;

    if (!region_str.empty()) {
        has_region = true;
        size_t colon_pos = region_str.find(':');
        size_t dash_pos  = region_str.find('-');
        if (colon_pos == std::string::npos) {
            // Whole chromosome
            query_chrom = region_str;
        } else {
            query_chrom = region_str.substr(0, colon_pos);
            if (dash_pos != std::string::npos && dash_pos > colon_pos) {
                query_start = std::stoul(region_str.substr(colon_pos + 1, dash_pos - colon_pos - 1));
                query_end   = std::stoul(region_str.substr(dash_pos + 1));
            } else {
                query_start = std::stoul(region_str.substr(colon_pos + 1));
            }
        }
        std::cout << "  Region filter: " << query_chrom << ":"
                  << query_start << "-" << query_end << "\n";
    }

    if (limit > 0) {
        std::cout << "  Record limit:  " << limit << "\n";
    }

    std::cout << "\n";

    // Read and display records
    uint32_t record_count = 0;
    uint32_t displayed = 0;
    std::string ref_id;
    uint32_t ref_pos;
    int total_depth;
    std::vector<BaseType::BatchInfo> samples_bi;

    // Print column header
    if (verbose) {
        std::cout << "--- Records ---\n";
    } else {
        std::cout << "Position         RefBase  TotalDepth\n";
        std::cout << "---------------  -------  ----------\n";
    }

    while (read_binary_record(bf, sample_count, ref_id, ref_pos, total_depth, samples_bi)) {
        ++record_count;

        // Apply region filter
        if (has_region) {
            if (ref_id != query_chrom) {
                // If we've passed the target chromosome, stop
                if (!query_chrom.empty() && ref_id > query_chrom) break;
                continue;
            }
            if (ref_pos < query_start) continue;
            if (ref_pos > query_end) break;
        }

        // Apply limit
        if (limit > 0 && displayed >= (uint32_t)limit) break;

        if (verbose) {
            std::cout << ref_id << ":" << ref_pos
                      << "  depth=" << total_depth << "\n";

            for (uint16_t s = 0; s < sample_count; ++s) {
                const auto &smp = samples_bi[s];
                uint16_t depth = static_cast<uint16_t>(smp.align_bases.size());

                // Check if this is a "no data" sample (single "N" entry)
                if (depth == 1 && smp.align_bases[0] == "N") {
                    std::cout << "  [" << s << "] " << sample_ids[s]
                              << ": depth=0 (no data)\n";
                } else {
                    std::cout << "  [" << s << "] " << sample_ids[s]
                              << ": depth=" << depth;
                    std::cout << "  ref_bases=[";
                    for (size_t j = 0; j < smp.ref_bases.size(); ++j) {
                        if (j > 0) std::cout << ",";
                        std::cout << smp.ref_bases[j];
                    }
                    std::cout << "]  align_bases=[";
                    for (size_t j = 0; j < smp.align_bases.size(); ++j) {
                        if (j > 0) std::cout << ",";
                        std::cout << smp.align_bases[j];
                    }
                    std::cout << "]  quals=[";
                    for (size_t j = 0; j < smp.align_base_quals.size(); ++j) {
                        if (j > 0) std::cout << ",";
                        std::cout << (int)(smp.align_base_quals[j] - 33);  // Phred quality
                    }
                    std::cout << "]  rpr=[";
                    for (size_t j = 0; j < smp.base_pos_ranks.size(); ++j) {
                        if (j > 0) std::cout << ",";
                        std::cout << smp.base_pos_ranks[j];
                    }
                    std::cout << "]  mapq=[";
                    for (size_t j = 0; j < smp.mapqs.size(); ++j) {
                        if (j > 0) std::cout << ",";
                        std::cout << smp.mapqs[j];
                    }
                    std::cout << "]  strand=[";
                    for (size_t j = 0; j < smp.map_strands.size(); ++j) {
                        if (j > 0) std::cout << ",";
                        std::cout << smp.map_strands[j];
                    }
                    std::cout << "]\n";
                }
            }
        } else {
            // Compact summary: show ref_base from first non-empty sample
            std::string ref_base_str = ".";
            for (const auto &smp : samples_bi) {
                if (!smp.ref_bases.empty() && smp.ref_bases[0] != "N") {
                    ref_base_str = smp.ref_bases[0];
                    break;
                }
            }

            printf("%-15s  %-7s  %d\n",
                   (ref_id + ":" + std::to_string(ref_pos)).c_str(),
                   ref_base_str.c_str(),
                   total_depth);
        }

        ++displayed;
    }

    // Summary
    std::cout << "\nRecords read:        " << record_count;
    if (has_region || (limit > 0 && displayed < record_count)) {
        std::cout << " (partial - filtered)";
    }
    std::cout << "\n";
    if (displayed != record_count) {
        std::cout << "Records displayed:   " << displayed << " (filtered)\n";
    }

    return 0;
}

// =====================================================================
//  Entry point
// =====================================================================

int dump_runner(int argc, char* argv[]) {
    if (argc < 2) {
        dump_usage();
        return 1;
    }

    // Parse options
    std::string input_file;
    bool header_only = false;
    std::string region_str;
    int limit = 0;
    bool verbose = false;
    bool show_entries = false;

    static struct option long_options[] = {
        {"header",  no_argument,       nullptr, 'H'},
        {"region",  required_argument, nullptr, 'r'},
        {"limit",   required_argument, nullptr, 'n'},
        {"verbose", no_argument,       nullptr, 'v'},
        {"entries", no_argument,       nullptr, 'E'},
        {"help",    no_argument,       nullptr, 'h'},
        {nullptr,   0,                 nullptr,  0 }
    };

    optind = 1;  // reset getopt
    int c;
    while ((c = getopt_long(argc, argv, "r:n:vh", long_options, nullptr)) != -1) {
        switch (c) {
            case 'H': header_only = true; break;
            case 'r': region_str = optarg; break;
            case 'n': limit = std::atoi(optarg); break;
            case 'v': verbose = true; break;
            case 'E': show_entries = true; break;
            case 'h': dump_usage(); return 0;
            default:  dump_usage(); return 1;
        }
    }

    // Remaining argument is the input file
    if (optind < argc) {
        input_file = argv[optind];
    }

    if (input_file.empty()) {
        std::cerr << "[ERROR] No input file specified.\n" << std::endl;
        dump_usage();
        return 1;
    }

    // Auto-detect file type by extension
    bool is_bbi = false;
    if (input_file.size() >= 4 && input_file.substr(input_file.size() - 4) == ".bbi") {
        is_bbi = true;
    } else if (input_file.size() >= 4 && input_file.substr(input_file.size() - 4) == ".bbf") {
        is_bbi = false;
    } else {
        // Try to detect by reading magic bytes
        FILE *fp = fopen(input_file.c_str(), "rb");
        if (!fp) {
            std::cerr << "[ERROR] Cannot open file: " << input_file << std::endl;
            return 1;
        }
        uint32_t magic = 0;
        if (fread(&magic, sizeof(magic), 1, fp) == 1) {
            if (magic == BBI_MAGIC) {
                is_bbi = true;
            } else if (magic == BBF_MAGIC) {
                is_bbi = false;
            } else {
                std::cerr << "[ERROR] Cannot determine file type. "
                          << "File should end with .bbf or .bbi, or have a valid magic number.\n"
                          << "  Got magic: 0x" << std::hex << magic << std::dec << std::endl;
                fclose(fp);
                return 1;
            }
        } else {
            std::cerr << "[ERROR] File too small to identify: " << input_file << std::endl;
            fclose(fp);
            return 1;
        }
        fclose(fp);
    }

    if (is_bbi) {
        return dump_bbi(input_file, show_entries);
    } else {
        return dump_bbf(input_file, header_only, region_str, limit, verbose);
    }
}

}  // namespace dump
}  // namespace basevar
