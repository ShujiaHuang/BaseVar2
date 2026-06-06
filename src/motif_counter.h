/**
 * @file motif_counter.h
 *
 * @brief cfDNA end-motif (k-mer) counting subcommand for BaseVar.
 *
 *   Extracts the 5' end-motif of each cfDNA fragment from BAM/CRAM
 *   alignments and reports per-motif counts and frequencies.  Each input
 *   BAM/CRAM file is treated as one sample, and the per-sample results
 *   are emitted side-by-side in a single TSV.  Designed for
 *   non-invasive prenatal testing (NIPT) and other low-pass cfDNA studies,
 *   following the convention used by Jiang et al., PNAS 2020.
 *
 *   This module reuses BaseVar's I/O abstractions:
 *     - ngslib::Bam / ngslib::BamRecord (src/io/bam.h)
 *     - ngslib::BamHeader::get_sample_name (src/io/bam_header.h)
 *     - ngslib::get_firstcolumn_from_file / split (src/io/utils.h)
 *     - ThreadPool                              (src/external/thread_pool.h)
 *
 * @author Shujia Huang
 * @date   2026-05-26
 */
#ifndef __INCLUDE_BASEVAR_MOTIF_COUNTER_H__
#define __INCLUDE_BASEVAR_MOTIF_COUNTER_H__

#include <map>
#include <string>
#include <vector>
#include <thread>
#include <cstdint>

#include "io/bam.h"
#include "io/bam_record.h"
#include "io/bam_header.h"
#include "io/fasta.h"

namespace basevar {
namespace motif {

/**
 * @brief Reverse-complement a single DNA base.
 *        Any character that is not one of A/C/G/T is mapped to 'N'.
 */
char reverse_complement_base(char base);

/**
 * @brief Reverse-complement an entire DNA string.
 */
std::string reverse_complement(const std::string& s);

/**
 * @brief Extract the 5' end-motif (k-mer) of a cfDNA fragment from a read.
 *
 *   - Forward-mapped read : take the first  k bases of the BAM SEQ field.
 *   - Reverse-mapped read : take the last   k bases of the BAM SEQ field
 *                           and reverse-complement them.
 *
 *   The BAM SEQ field is stored in reference-forward orientation, so the
 *   above procedure recovers the first k bases of the original sequenced
 *   read, which corresponds to the 5' end of the cfDNA fragment.
 *
 * @param r  A mapped BamRecord.
 * @param k  Motif length (>= 1).
 * @return   The k-mer string, or an empty string if the read is shorter
 *           than k or unmapped.
 */
std::string extract_5p_motif(const ngslib::BamRecord& r, int k);

/**
 * @brief Reference-based variant of extract_5p_motif (Lo lab canonical path).
 *
 *   Fetch the k bases at the read's 5'-aligned position directly from the
 *   reference genome, ignoring the read's own sequenced bases.
 *
 *   - Forward-mapped read : ref[chrom, pos, pos + k)
 *   - Reverse-mapped read : reverse_complement(ref[chrom, end - k, end))
 *
 *   This is the canonical method established by the Lo lab (Jiang et al.,
 *   Cancer Discovery 2020) and faithfully reproduced by FinaleToolkit
 *   (Zheng et al., bioRxiv 2024.05.29.596414).  The original Jiang 2020
 *   source code was never released, but FinaleToolkit -- the de-facto
 *   reference open-source implementation -- always fetches end motifs
 *   from the FASTA, never from the BAM SEQ field.  The Lo group's
 *   follow-up paper (Mao et al., Cell Genomics 2026) similarly states
 *   that EM5/EM3 motifs "were deduced from the reference genome".
 *
 *   By contrast, `extract_5p_motif` reads bases directly from the BAM
 *   SEQ field; that path is BaseVar's conservative default (no FASTA
 *   required) and is itself a deviation from the canonical method.
 *
 * @param r    A mapped BamRecord.
 * @param fa   FASTA accessor for the reference genome (must be loaded).
 * @param hdr  BAM header used to translate tid to chromosome name.
 * @param k    Motif length (>= 1).
 * @return     The k-mer string (uppercase A/C/G/T/N), or an empty string
 *             on failure (out-of-range fetch, missing chromosome, etc.).
 */
std::string extract_5p_motif_from_reference(const ngslib::BamRecord& r,
                                            const ngslib::Fasta&     fa,
                                            const ngslib::BamHeader& hdr,
                                            int k);

/**
 * @brief Read-end selection policy for paired-end data.
 */
enum class ReadEnd { R1, R2, BOTH };

/**
 * @brief Per-sample motif counts and bookkeeping.  One instance is produced
 *        per input BAM/CRAM file.
 */
struct SampleResult {
    std::string sample_id;        ///< Sample ID (from @RG SM tag, or filename).
    std::string input_path;       ///< Source BAM/CRAM path.
    uint64_t total_reads     = 0; ///< All reads scanned.
    uint64_t filtered_reads  = 0; ///< Reads that failed the filters.
    uint64_t used_reads      = 0; ///< Reads that contributed to motif counts.
    uint64_t n_motifs_with_n = 0; ///< Reads whose extracted motif contained a non-ACGT base (N or IUPAC ambiguity code, e.g. M/R/Y) - excluded.
    std::map<std::string, uint64_t> motif_counts;  ///< Always populated with all 4^k keys.
    std::vector<double> fprofile_weights;          ///< F-profile weights (size 6 if computed, empty otherwise).
};

/**
 * @brief Command-line options consumed by `basevar motif`.
 */
struct MotifArgs {
    std::vector<std::string> input_bf;   // Positional + -L list combined.
    std::string in_filelist;             // -L
    std::string output_file;             // -o (TSV); empty -> only stdout summary
    std::string regions;                 // -r (comma-separated, optional)
    std::string reference;               // -f (required for CRAM)
    int         motif_length;            // -l, default 4 (range [1, 10])
    int         min_mapq;                // -q, default 30
    int         thread_num;              // -t, default = hardware concurrency
    ReadEnd     read_end;                // --reads (R1/R2/both), default R1
    bool        include_zero;            // --include-zero, default true
    bool        filename_has_samplename; // --filename-has-samplename, default false
    bool        proper_pair;             // --proper-pair: reject reads not in proper pairs (PE only)
    int         max_insert_size;         // --max-insert-size: reject |isize|>this; 0=disabled (PE only)
    bool        from_reference;          // --from-reference: extract motif from reference genome
                                         //                   instead of the read's own bases.
                                         //                   Requires --reference. Default off.
    bool        fprofile;                // --fprofile: compute F-profile decomposition (k=4 only).
    std::string fprofile_output;         // --fprofile-output FILE: write per-sample F-profile weights TSV.

    MotifArgs()
        : motif_length(4),
          min_mapq(30),
          thread_num(static_cast<int>(std::thread::hardware_concurrency())),
          read_end(ReadEnd::R1),
          include_zero(true),
          filename_has_samplename(false),
          proper_pair(false),
          max_insert_size(0),
          from_reference(false),
          fprofile(false)
    {
        if (thread_num < 1) thread_num = 1;
    }
};

/**
 * @brief Runner class for the `basevar motif` subcommand.
 *
 *   Mirrors the BaseTypeRunner pattern: parses argv in the constructor,
 *   stores parameters, and exposes a single run() entry point.  Each input
 *   BAM/CRAM file is processed in its own worker thread (ThreadPool) and
 *   produces an independent SampleResult.
 */
class MotifCounterRunner {
public:
    MotifCounterRunner(int argc, char* argv[]);
    ~MotifCounterRunner();

    /// Help text for `basevar motif --help`.
    static const std::string usage();

    /// Execute the subcommand.  Returns 0 on success, non-zero on error.
    int run();

    // --- Read-only accessors used by unit tests ---
    const MotifArgs& args() const { return *_args; }
    const std::vector<SampleResult>& results() const { return _results; }

    // Aggregate counters across all samples.
    uint64_t total_reads()      const;
    uint64_t filtered_reads()   const;
    uint64_t used_reads()       const;
    uint64_t n_motifs_with_n()  const;

private:
    MotifArgs* _args;
    std::string _cmdline_string;
    std::vector<SampleResult> _results;

    /// Resolve a sample ID for one file (filename or @RG SM tag).
    std::string _resolve_sample_id(const std::string& fn) const;

    /// Pre-populate a SampleResult's motif_counts with all 4^k zero entries.
    void _initialize_motifs_for(SampleResult& s) const;

    /// Worker: read one BAM/CRAM file end-to-end and fill the SampleResult.
    /// Thread-safe: only touches local state plus the supplied result struct.
    void _process_one_file(SampleResult& s) const;

    /// True if the read should contribute to the motif count.
    bool _passes_filters(const ngslib::BamRecord& r) const;

    /// Write the per-sample TSV file.
    void _write_tsv(const std::string& path) const;

    /// Write per-sample F-profile weights to a separate TSV file.
    void _write_fprofile_tsv(const std::string& path) const;

    /// Print a human-readable per-sample summary (English).
    void _print_summary(std::ostream& os) const;

    MotifCounterRunner(const MotifCounterRunner&) = delete;
    MotifCounterRunner& operator=(const MotifCounterRunner&) = delete;
};

/**
 * @brief Thin C-style entry point used by the BaseVar dispatcher in main.cpp.
 *        Mirrors `concat_runner(argc, argv)`.
 */
int motif_counter_runner(int argc, char* argv[]);

}  // namespace motif
}  // namespace basevar

#endif  // __INCLUDE_BASEVAR_MOTIF_COUNTER_H__
