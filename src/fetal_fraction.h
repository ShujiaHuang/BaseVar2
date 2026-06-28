/**
 * @file fetal_fraction.h
 *
 * @brief Y-chromosome based male fetal-fraction (FF) estimation for
 *        NIPT cell-free DNA (cfDNA) BAM/CRAM data.
 *
 *   This is a HIDDEN subcommand of BaseVar (`basevar fetalfrac`).  It
 *   is intentionally not advertised in the top-level `basevar --help`
 *   output: it is a private utility for the author's own NIPT studies.
 *
 *   ## Method (chrY-counting, male fetus only)
 *
 *   For ultra-low-pass NIPT WGS the canonical FF estimator is:
 *
 *      Y_ratio = N_chrY / N_autosomal
 *      FF      = max(0, Y_ratio - noise) * scaling_factor
 *
 *   where:
 *
 *     - N_autosomal counts reads on chr1..chr22 (excluding chrX/chrY/chrM
 *       and ALT/decoy contigs).  Autosomes are diploid in BOTH mother
 *       and fetus, so they cancel out and serve as a stable denominator.
 *
 *     - N_chrY counts reads on chrY OUTSIDE the pseudoautosomal regions
 *       (PAR1/PAR2).  PAR is shared with chrX and would inflate the Y
 *       count for female fetuses.
 *
 *     - scaling_factor is empirically calibrated from male reference
 *       samples.  For GRCh38 with ~2.87 Gb autosomes and ~57 Mb chrY,
 *       a 100% male fetus would yield Y_ratio ~= chrY_len / (2 *
 *       autosome_len) ~= 0.00993, so FF (frac) ~= Y_ratio * 100.7.
 *       The default of 100 is a reasonable starting point but USERS
 *       SHOULD recalibrate with their own male/female reference
 *       cohorts and pass `--scale <value>`.
 *
 *     - The female refusal threshold (`--male-threshold`) lets us
 *       detect a female fetus from background-only Y signal (index
 *       hopping, mismapping).  When Y_ratio falls below the threshold
 *       the sample is reported as FEMALE and FF is NOT computed --
 *       see the scientific rationale in the runner's `usage()` text.
 *
 *   ## What this module reuses from BaseVar
 *
 *     - ngslib::Bam / ngslib::BamRecord (src/io/bam.h)        : BAM/CRAM I/O
 *     - ngslib::BamHeader::get_sample_name (src/io/bam_header.h)
 *     - ngslib::get_firstcolumn_from_file / split (src/io/utils.h)
 *     - ThreadPool                              (src/external/thread_pool.h)
 *
 *   The original (now-removed) `ff_calculator.cpp` simulated reads with
 *   a random number generator and never touched a real BAM file; this
 *   replacement reads BAM/CRAM directly through BaseVar's I/O layer
 *   exactly like `basevar caller` and `basevar motif`.
 *
 * @author Shujia Huang
 * @date   2026-05-27
 */
#ifndef __INCLUDE_BASEVAR_FETAL_FRACTION_H__
#define __INCLUDE_BASEVAR_FETAL_FRACTION_H__

#include <functional>
#include <map>
#include <string>
#include <vector>
#include <thread>
#include <cstdint>

#include "io/bam_record.h"
#include "io/utils.h"  // ngslib::GenomeRegion (1-based, inclusive)

namespace basevar {
namespace fetalfrac {

/**
 * @brief Strip a leading "chr"/"CHR" prefix (case-insensitive).
 */
std::string strip_chr_prefix(const std::string& chrom);

/**
 * @brief True iff `chrom` is one of chrY / Y / y (case-insensitive).
 */
bool is_y_chrom(const std::string& chrom);

/**
 * @brief True iff `chrom` is one of chr1..chr22 / 1..22.
 *        Rejects chrX / chrY / chrM / chrMT and any ALT/decoy/HLA/random contig.
 */
bool is_autosome(const std::string& chrom);

/**
 * @brief Inferred fetal sex from the Y-chromosome read ratio.
 */
enum class FetalSex { MALE, FEMALE, UNDETERMINED };

/**
 * @brief Lightweight per-chromosome interval index loaded from a BED file.
 *
 *   Provides a sorted+merged set of intervals per chromosome and an
 *   O(log N) `contains(chrom, pos)` query.  Used both for:
 *     - inclusion regions (mappability BED), and
 *     - exclusion regions (PAR BED).
 *
 *   ## Coordinate convention
 *
 *   To stay consistent with the rest of BaseVar (e.g.
 *   `pipeline::parse_region_string`, `ngslib::GenomeRegion`), the
 *   index stores intervals as **1-based, inclusive** [start, end].
 *   The BED loader converts the on-disk 0-based half-open
 *   `[s, e)` representation to `[s+1, e]` automatically; the
 *   `add_interval` and `contains` API therefore both speak
 *   1-based inclusive.  Callers that hold 0-based BAM positions
 *   (e.g. `BamRecord::map_ref_start_pos`) must add 1 before
 *   querying `contains`.
 *
 *   When no BED is loaded, `contains()` always returns false and
 *   `is_loaded()` is false (the caller decides whether the absence
 *   of a BED means "accept all" or "reject all").
 */
class BedIntervalIndex {
public:
    BedIntervalIndex() : _is_loaded(false) {}

    /// Load a BED file (chrom\tstart\tend\t...).  Converts the on-disk
    /// 0-based half-open coordinates to 1-based inclusive internally.
    /// Throws on I/O error.
    void load_bed(const std::string& filepath);

    /// Add a single 1-based, inclusive interval programmatically (used
    /// for built-in PAR coordinates when no PAR BED was supplied).
    /// `start` and `end` are both inclusive; empty inputs (end < start)
    /// are silently dropped.
    void add_interval(const std::string& chrom, uint32_t start, uint32_t end);

    /// Finalise: sort + merge overlapping intervals.  Must be called
    /// after `add_interval(...)` programmatic loading.
    void finalize();

    bool is_loaded() const { return _is_loaded; }

    /// True iff `pos` (1-based) falls inside any merged interval on `chrom`.
    bool contains(const std::string& chrom, uint32_t pos) const;

    /// Read-only access to the merged-interval map (used by callers that
    /// need to compute lengths / intersections, e.g. the geometric --scale
    /// auto-computation in FetalFractionRunner).
    const std::map<std::string, std::vector<ngslib::GenomeRegion>>&
    regions() const { return _regions; }

    /// Sum of merged-interval lengths (in bp, 1-based inclusive) across
    /// every chromosome whose name satisfies `pred`.  Returns 0 when no
    /// matching intervals exist.
    uint64_t total_length(
        const std::function<bool(const std::string&)>& pred) const;

private:
    // Per-chromosome list of merged GenomeRegions.  The chrom field of
    // each stored GenomeRegion matches the map key (set redundantly so
    // that `to_string()` works on a stand-alone region during logging).
    std::map<std::string, std::vector<ngslib::GenomeRegion>> _regions;
    bool _is_loaded;

    static void _merge(std::vector<ngslib::GenomeRegion>& v);
};

/**
 * @brief Per-sample FF result.  One instance per input BAM/CRAM.
 */
struct FetalFractionResult {
    std::string sample_id;       ///< From @RG SM tag, or filename.
    std::string input_path;      ///< Source BAM/CRAM path.
    uint64_t total_reads     = 0;///< All reads scanned.
    uint64_t filtered_reads  = 0;///< Reads that failed QC filters.
    uint64_t valid_autosomal = 0;///< Autosomal reads that contributed (denominator).
    uint64_t valid_y         = 0;///< chrY-non-PAR reads that contributed (numerator).
    uint64_t y_par_excluded  = 0;///< chrY reads dropped because they fall in PAR.
    double   y_ratio         = 0.0; ///< valid_y / valid_autosomal
    FetalSex sex             = FetalSex::UNDETERMINED;
    double   fetal_fraction  = -1.0;///< -1 if not computed (FEMALE / UNDETERMINED).
};

/**
 * @brief Reference-build preset for built-in PAR coordinates.
 *
 *   When the user does not supply `--par-bed`, we apply the canonical
 *   PAR1/PAR2 coordinates for the chosen build.  Users running on a
 *   non-standard reference (T2T, custom assembly) should supply their
 *   own PAR BED instead.
 */
enum class GenomeBuild { GRCh38, GRCh37, NONE };

/**
 * @brief Command-line options for `basevar fetalfrac`.
 */
struct FetalFractionArgs {
    std::vector<std::string> input_bf;     // Positional + -L list combined.
    std::string in_filelist;               // -L
    std::string output_file;               // -o (TSV); empty => stdout-only.
    std::string regions;                   // -r (comma-separated, optional).
    std::string reference;                 // -f (required for CRAM).
    std::string mappability_bed;           // -B (optional).
    std::string par_bed;                   // --par-bed (optional).
    GenomeBuild genome_build;              // --build {hg38,hg19,none}, default hg38.
    int         min_mapq;                  // -q, default 30.
    int         thread_num;                // -t, default = hardware concurrency.
    bool        proper_pair;               // --proper-pair (PE only).
    int         max_insert_size;           // --max-insert-size; 0 = disabled.
    double      scaling_factor;            // --scale; auto if user does not set it.
    bool        scale_user_set;            // true iff --scale was explicitly given.
    double      y_noise;                   // --noise, default 0.0.
    double      male_threshold;            // --male-threshold, default 1e-4.
    bool        filename_has_samplename;   // --filename-has-samplename, default false.
    bool        calibrate;                 // --calibrate, default false.

    FetalFractionArgs()
        : genome_build(GenomeBuild::GRCh38),
          min_mapq(30),
          thread_num(static_cast<int>(std::thread::hardware_concurrency())),
          proper_pair(false),
          max_insert_size(0),
          scaling_factor(100.0),  // safe fallback if auto-compute fails
          scale_user_set(false),
          y_noise(0.0),
          male_threshold(1e-4),
          filename_has_samplename(false),
          calibrate(false)
    {
        if (thread_num < 1) thread_num = 1;
    }
};

/**
 * @brief Runner class for the `basevar fetalfrac` (hidden) subcommand.
 *
 *   Mirrors the BaseTypeRunner / MotifCounterRunner pattern: parses
 *   argv in the constructor, stores parameters, exposes a single
 *   `run()` entry point.  Each input BAM/CRAM file is processed in
 *   its own worker thread (ThreadPool) and produces an independent
 *   FetalFractionResult.
 */
class FetalFractionRunner {
public:
    FetalFractionRunner(int argc, char* argv[]);
    ~FetalFractionRunner();

    /// Help text for `basevar fetalfrac --help`.
    static const std::string usage();

    /// Execute the subcommand.  Returns 0 on success, non-zero on error.
    int run();

    // --- Read-only accessors used by unit tests ---
    const FetalFractionArgs&                args()    const { return *_args; }
    const std::vector<FetalFractionResult>& results() const { return _results; }

private:
    FetalFractionArgs* _args;
    std::string _cmdline_string;
    std::vector<FetalFractionResult> _results;
    BedIntervalIndex _mappability;  ///< inclusion regions  (optional)
    BedIntervalIndex _par_regions;  ///< exclusion regions  (built-in or BED)

    /// Resolve a sample ID for one file (filename or @RG SM tag).
    std::string _resolve_sample_id(const std::string& fn) const;

    /// Worker: read one BAM/CRAM file end-to-end and fill the result.
    /// Thread-safe: only touches local state plus the supplied result struct.
    void _process_one_file(FetalFractionResult& s) const;

    /// True if the read should contribute to the autosomal/Y counters.
    bool _passes_filters(const ngslib::BamRecord& r) const;

    /// Compute Y-ratio, sex call and FF from the per-sample counters.
    void _compute_fetal_fraction(FetalFractionResult& s) const;

    /// Populate `_par_regions` from the chosen --build preset.
    void _load_builtin_par();

    /// Compute the geometric default for `--scale` from the currently
    /// loaded mappability BED (if any) plus PAR exclusion, or from the
    /// canonical autosome / chrY-non-PAR lengths of `--build` when no
    /// `-B` was supplied.  Returns 0.0 when geometry cannot be inferred
    /// (e.g. `--build none` with no `-B`, or a `-B` BED that lacks
    /// autosomal or chrY non-PAR coverage).  Callers must fall back to
    /// the constructor-time default in that case.
    double _compute_geometric_scale() const;

    /// Write the per-sample TSV file.
    void _write_tsv(const std::string& path) const;

    /// Print a human-readable per-sample summary to `os` (English).
    void _print_summary(std::ostream& os) const;

    /// Calibration mode: process all samples, compute empirical scale
    /// from the mean y_ratio of samples passing the male threshold,
    /// and report the calibrated scale.  Returns 0 on success.
    int _run_calibrate();

    FetalFractionRunner(const FetalFractionRunner&) = delete;
    FetalFractionRunner& operator=(const FetalFractionRunner&) = delete;
};

/**
 * @brief Thin C-style entry point used by the BaseVar dispatcher in main.cpp.
 *        Mirrors `concat_runner(argc, argv)` and `motif_counter_runner`.
 */
int fetal_fraction_runner(int argc, char* argv[]);

}  // namespace fetalfrac
}  // namespace basevar

#endif  // __INCLUDE_BASEVAR_FETAL_FRACTION_H__
