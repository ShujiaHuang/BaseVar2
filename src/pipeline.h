/**
 * @file pipeline.h
 *
 * @brief Pipeline generator for `basevar caller`.
 *
 *  Splits the genome into sub-regions and prints one `basevar caller`
 *  command per sub-region to stdout.
 *
 *  This module reuses:
 *    - ngslib::GenomeRegion / ngslib::split / ngslib::is_readable (io/utils.h)
 *
 *  @author Shujia Huang
 *  @date   2026-05-19
 */
#ifndef __INCLUDE_BASEVAR_PIPELINE_H__
#define __INCLUDE_BASEVAR_PIPELINE_H__

#include <string>
#include <vector>
#include <stdexcept>
#include <map>
#include <cstdint>

#include "io/utils.h"  // ngslib::GenomeRegion

namespace basevar {
namespace pipeline {

/**
 * @brief Load chromosome entries from a samtools-style .fai index file.
 *
 *  Each line of a .fai file looks like:
 *    chr1   249250621   52   60   61
 *  Only the first two columns (name, length) are used here.
 *
 * @param fai_path       Path to the .fai file.
 * @param chroms_filter  If non-empty, only chromosomes whose names appear
 *                       in this list are returned.
 * @return Vector of (chrom, 1, length) GenomeRegions in file order.
 *
 * @throws std::runtime_error if the file cannot be opened or is malformed.
 */
std::vector<ngslib::GenomeRegion> load_reference_fai(
    const std::string& fai_path,
    const std::vector<std::string>& chroms_filter = {});

/**
 * @brief Parse a single region string into a 1-based, inclusive GenomeRegion.
 *
 *  Supported formats:
 *    'chr'              -> (chr, 1, chrom_length)
 *    'chr:start'        -> (chr, start, chrom_length)
 *    'chr:start-end'    -> (chr, start, end)
 *
 * @param region_str   Single region string (no commas).
 * @param ref_fai_map  Map of {chrom -> length} from .fai for validation.
 *
 * @throws std::invalid_argument on malformed string, unknown chromosome,
 *         start < 1, start > end, or end > chrom_length.
 */
ngslib::GenomeRegion parse_region_string(
    const std::string& region_str,
    const std::map<std::string, uint32_t>& ref_fai_map);

/**
 * @brief Parse a comma-separated regions string into a list of GenomeRegion.
 *
 * @param regions_str    e.g. "chr1,chr2:1000-2000".
 * @param ref_fai_map    Chromosome length map for validation.
 * @param chroms_filter  Optional whitelist; empty means keep all.
 */
std::vector<ngslib::GenomeRegion> parse_regions(
    const std::string& regions_str,
    const std::map<std::string, uint32_t>& ref_fai_map,
    const std::vector<std::string>& chroms_filter = {});

/**
 * @brief Build the per-sub-region `basevar caller` command string.
 *
 *  The returned command string follows the format:
 *   time <exe_prog> <passthrough_str> -r <reg> -o <outdir>/<prefix>.vcf.gz \
 *        > <outdir>/<prefix>.log && echo "** <prefix> done **"
 */
std::string build_command(const std::string& exe_prog,
                          const std::string& passthrough_str,
                          const std::string& chrom,
                          uint32_t start,
                          uint32_t end,
                          const std::string& outdir);

/**
 * @brief Pipeline-specific options consumed by this subcommand.
 *        Every other argument is passed through verbatim to `basevar caller`.
 */
struct PipelineArgs {
    std::string outdir;           // -o/--outdir   (required)
    std::string ref_fai;          // --ref_fai     (required)
    uint32_t    delta = 2000000;  // -d/--delta    sub-region size in bp
    std::string chrom;            // -c/--chrom    comma-separated whitelist
};

/**
 * @brief Pipeline subcommand entry point.
 *
 *  argv layout (as received from main): argv[0] == "pipeline", then options.
 *
 * @param argc  Argument count.
 * @param argv  Argument vector.
 * @param basevar_executable  Absolute path (or argv[0]) of the running
 *                            `basevar` binary; used as the program in the
 *                            generated commands.  Falls back to "basevar"
 *                            when empty.
 *
 * @return Exit code (0 on success, non-zero on error).
 */
int pipeline_runner(int argc, char* argv[],
                    const std::string& basevar_executable = "");

}  // namespace pipeline
}  // namespace basevar

#endif
