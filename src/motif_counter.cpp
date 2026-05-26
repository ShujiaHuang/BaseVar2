/**
 * @file motif_counter.cpp
 *
 * @brief Implementation of the `basevar motif` subcommand.
 *
 * ## Canonical method (Lo lab convention)
 *
 *   The cfDNA end-motif analysis implemented here follows the method
 *   established by Prof. Y.M. Dennis Lo's group, in particular:
 *
 *     Jiang P., Sun K., Peng W., Cheng S.H., Ni M., Yeung P.C., Heung
 *     M.M.S., Xie T., Shang H., Zhou Z., Chan R.W.Y., Wong J., Wong
 *     V.W.S., Poon L.C., Leung T.Y., Lai P.B.S., Chan H.L.Y., Chiu R.W.K.,
 *     Chan K.C.A., Lo Y.M.D.  "Plasma DNA End-Motif Profiling as a
 *     Fragmentomic Marker in Cancer, Pregnancy, and Transplantation."
 *     Cancer Discovery, 2020;10(5):664-673.  PMID: 31911138.
 *
 *   Key elements of that canonical method:
 *
 *     1. The end-motif is the **first k bases (k=4 by default) at the
 *        5'-most aligned position of the cfDNA fragment, fetched from
 *        the REFERENCE genome** (NOT from the read SEQ).  The original
 *        Jiang 2020 source code was never released; the convention has
 *        been confirmed by:
 *
 *          - the Lo group's follow-up paper Mao et al., Cell Genomics
 *            2026, which states that EM5/EM3 motifs "were deduced from
 *            the reference genome";
 *          - FinaleToolkit (Zheng et al., bioRxiv 2024.05.29.596414),
 *            the open-source reference implementation, which fetches
 *            ref[chrom, start, start+k) for forward fragments and
 *            rc(ref[chrom, end-k, end)) for reverse fragments.
 *
 *        Use `--from-reference -f <fa>` to reproduce that path here.
 *     2. **Both** fragment ends per pair: in paired-end data, *each*
 *        cfDNA fragment contributes TWO end-motifs (5' of R1 AND 5'
 *        of R2).  In this tool, that is achieved with `--reads both`.
 *     3. **Properly paired** reads only: enable with `--proper-pair`.
 *     4. MAPQ >= 30 (matches our default).
 *     5. Discard chimeric / discordantly-mapped reads via an insert-size
 *        cap (the Lo lab pipeline typically uses |isize| <= ~1000 bp).
 *        Enable here with e.g. `--max-insert-size 1000`.
 *     6. Motifs containing any non-ACGT base (N or any IUPAC ambiguity
 *        code) are excluded from numerator and denominator.
 *     7. Motif Diversity Score (MDS) = Shannon entropy of the 4^k motif
 *        distribution divided by log2(4^k) = 2k.
 *
 *   ### Recommended invocation to reproduce the Lo lab convention:
 *
 *     basevar motif --from-reference -f ref.fa \
 *                   --reads both --proper-pair --max-insert-size 1000 \
 *                   -q 30 -l 4 -o out.tsv  in1.bam in2.bam ...
 *
 *   The tool's *defaults* deliberately stay conservative (read-derived
 *   bases via `extract_5p_motif`, `--reads R1`, `--proper-pair` off, no
 *   insert-size cap) so that BAM/CRAM inputs from non-cfDNA workflows
 *   (e.g. genome-wide variant calling) still produce sensible results
 *   without requiring a FASTA.  Note that the read-derived path is
 *   itself a deviation from the canonical Lo lab method; users running
 *   cfDNA / NIPT / fragmentomic analyses are strongly encouraged to use
 *   the recommended invocation above.
 *
 * @author Shujia Huang
 * @date   2026-05-26
 */
#include "motif_counter.h"

#include <getopt.h>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <functional>
#include <future>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

#include "external/thread_pool.h"
#include "io/bam_header.h"
#include "io/fasta.h"
#include "io/utils.h"
#include "version.h"

namespace basevar {
namespace motif {

// ============================================================
// Free helpers
// ============================================================

char reverse_complement_base(char base) {
    switch (base) {
        case 'A': case 'a': return 'T';
        case 'C': case 'c': return 'G';
        case 'G': case 'g': return 'C';
        case 'T': case 't': return 'A';
        default:            return 'N';
    }
}

std::string reverse_complement(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (auto it = s.rbegin(); it != s.rend(); ++it) {
        out.push_back(reverse_complement_base(*it));
    }
    return out;
}

std::string extract_5p_motif(const ngslib::BamRecord& r, int k) {
    if (!r || !r.is_mapped() || k <= 0) return "";

    // BamRecord::query_sequence() decodes the BAM SEQ field (in
    // reference-forward orientation) using ngslib::_BASES.
    std::string seq = r.query_sequence();
    if (static_cast<int>(seq.size()) < k) return "";

    if (r.is_mapped_reverse()) {
        // Reverse-mapped read: the original sequencer's 5' end of the read
        // is the LAST base in the BAM SEQ field (which has been reverse-
        // complemented when stored).  Take the last k bases and rc them
        // to recover the first k bases of the original read.
        return reverse_complement(seq.substr(seq.size() - k, k));
    }
    // Forward-mapped read: SEQ is in original sequencer orientation.
    return seq.substr(0, k);
}

std::string extract_5p_motif_from_reference(const ngslib::BamRecord& r,
                                            const ngslib::Fasta&     fa,
                                            const ngslib::BamHeader& hdr,
                                            int k) {
    if (!r || !r.is_mapped() || k <= 0) return "";

    // Translate tid -> chromosome name.
    std::string chrom = r.tid_name(const_cast<ngslib::BamHeader&>(hdr));
    if (chrom.empty()) return "";

    const uint32_t chrom_len = fa.seq_length(chrom);
    if (chrom_len == 0) return "";

    std::string ref_kmer;
    try {
        if (r.is_mapped_reverse()) {
            // Reverse-mapped read: the original 5' end of the fragment maps
            // to the *right* end of the alignment on the reference.
            // map_ref_end_pos() is half-open (0-based exclusive end).
            const hts_pos_t end = r.map_ref_end_pos();
            if (end < k) return "";
            const uint32_t s = static_cast<uint32_t>(end - k);
            const uint32_t e = static_cast<uint32_t>(end);
            if (e > chrom_len) return "";
            // Fasta::fetch is half-open [start, end).  However, htslib's
            // faidx_fetch_seq treats `end` as INCLUSIVE; ngslib::Fasta::fetch
            // wraps it directly, so pass (e - 1) to mean "last position".
            ref_kmer = fa.fetch(chrom, s, e - 1);
        } else {
            // Forward-mapped read: 5' end is the leftmost aligned position.
            const hts_pos_t pos = r.map_ref_start_pos();
            if (pos < 0) return "";
            const uint32_t s = static_cast<uint32_t>(pos);
            if (s + static_cast<uint32_t>(k) > chrom_len) return "";
            ref_kmer = fa.fetch(chrom, s, s + static_cast<uint32_t>(k) - 1);
        }
    } catch (const std::exception&) {
        return "";
    }

    if (static_cast<int>(ref_kmer.size()) < k) return "";
    // Some references are soft-masked (lowercase repeats); normalise to
    // uppercase so is_pure_acgt() and the motif-counts map see canonical keys.
    for (char& c : ref_kmer) c = std::toupper(static_cast<unsigned char>(c));

    if (r.is_mapped_reverse()) {
        return reverse_complement(ref_kmer);
    }
    return ref_kmer;
}

// ============================================================
// Argument parsing helpers
// ============================================================

static ReadEnd parse_read_end(const std::string& s) {
    std::string v;
    v.reserve(s.size());
    for (char c : s) v.push_back(std::tolower(static_cast<unsigned char>(c)));
    if (v == "r1")   return ReadEnd::R1;
    if (v == "r2")   return ReadEnd::R2;
    if (v == "both") return ReadEnd::BOTH;
    throw std::invalid_argument("[ERROR] --reads must be one of: R1, R2, both (got '" + s + "')");
}

static const char* read_end_name(ReadEnd e) {
    switch (e) {
        case ReadEnd::R1:   return "R1";
        case ReadEnd::R2:   return "R2";
        case ReadEnd::BOTH: return "both";
    }
    return "?";
}

// Pre-populate all 4^k motifs with 0 in the supplied counter map.
static void enumerate_kmers(int k, std::map<std::string, uint64_t>& dst) {
    dst.clear();
    const std::string bases = "ACGT";
    std::string current(k, 'A');
    std::function<void(int)> gen = [&](int pos) {
        if (pos == k) { dst[current] = 0; return; }
        for (char b : bases) {
            current[pos] = b;
            gen(pos + 1);
        }
    };
    gen(0);
}

// Reject any motif that is not made entirely of A/C/G/T.
//
//   `BamRecord::query_sequence()` decodes the BAM 4-bit SEQ via
//   `ngslib::_BASES[16]`.  That table maps codes 1/2/4/8/15 to
//   A/C/G/T/N respectively, but **all other codes (the IUPAC
//   ambiguity codes M/R/W/S/Y/K/V/H/D/B = 3/5/6/7/9/10/11/12/13/14)
//   map to the space character ' '**, NOT to 'N'.  A naive
//   `find('N')` check therefore lets ambiguity-coded reads slip
//   through and produces malformed rows in the output TSV.
//   This helper is the single source of truth for "clean" motifs.
static inline bool is_pure_acgt(const std::string& s) {
    for (char c : s) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') return false;
    }
    return !s.empty();
}

// ============================================================
// MotifCounterRunner
// ============================================================

static const std::string MOTIF_USAGE =
    "About: Count cfDNA end-motif (k-mer) frequencies from BAM/CRAM alignments.\n"
    "       Each input BAM/CRAM is treated as one sample; per-sample counts are\n"
    "       emitted in a single TSV (long format).\n"
    "       Method follows Jiang et al., Cancer Discovery 2020 (PMID: 31911138).\n"
    "Usage: basevar motif [options] <-o output.tsv> [-L bam.list] in1.bam [in2.bam ...]\n\n"

    "Required arguments:\n"
    "  -o, --output FILE            Output TSV file (sample, motif, count, frequency).\n\n"

    "Optional arguments:\n"
    "  -L, --align-file-list FILE   BAM/CRAM files list, one path per row.\n"
    "  -f, --reference FILE         Reference FASTA file (required for CRAM input).\n"
    "  -r, --regions REG[,...]      Restrict counting to these regions, comma-separated.\n"
    "                               Formats: chr | chr:start | chr:start-end\n"
    "  -l, --length INT             Motif length k, range [1, 10] [4]\n"
    "  -q, --mapq INT               Minimum MAPQ to keep a read [30]\n"
    "  -t, --thread INT             Number of worker threads (one file per thread)\n"
    "                               [hardware_concurrency]\n"
    "      --reads {R1|R2|both}     Which read in a pair to use for the 5' end-motif [R1].\n"
    "                               The default `R1` is conservative and works for any\n"
    "                               BAM/CRAM. To reproduce the canonical Lo lab cfDNA\n"
    "                               end-motif method (Jiang et al., Cancer Discovery 2020,\n"
    "                               PMID: 32111602), use `--reads both` so that EACH\n"
    "                               fragment contributes TWO end-motifs (5' of R1 AND\n"
    "                               5' of R2).\n"
    "      --include-zero           Emit all 4^k motifs (zeros included) in TSV (default ON)\n"
    "      --no-include-zero        Suppress motifs with zero count in TSV.\n"
    "      --filename-has-samplename\n"
    "                               Derive sample IDs from filenames instead of reading the\n"
    "                               BAM @RG SM tag.  E.g. /path/SampleA.bam -> SampleA.\n"
    "      --proper-pair            Only count reads flagged as properly paired\n"
    "                               (BAM_FPROPER_PAIR).  Required by the Lo lab convention\n"
    "                               for cfDNA analyses; silently ignored for SE data. [off]\n"
    "      --max-insert-size INT    Discard reads whose |insert size| > INT. 0 = no limit.\n"
    "                               Lo lab cfDNA pipelines typically use 1000 to drop\n"
    "                               chimeric / discordantly-mapped reads.\n"
    "                               Silently ignored for SE data. [0]\n"
    "      --from-reference         Extract motifs from the REFERENCE genome at each\n"
    "                               fragment's 5' alignment position (canonical Lo lab\n"
    "                               cfDNA method, as implemented in FinaleToolkit /\n"
    "                               Zheng et al., bioRxiv 2024.05.29.596414).  Requires\n"
    "                               -f/--reference.  Recommended for cfDNA / NIPT /\n"
    "                               fragmentomic analyses.  When OFF (default), motifs\n"
    "                               are extracted from the read's own sequenced bases\n"
    "                               (BaseVar's conservative fallback that does not\n"
    "                               require a FASTA). [off]\n"
    "  -h, --help                   Show this help message and exit.\n\n"
    
    "Recommended cfDNA / NIPT invocation (Lo lab convention):\n"
    "  basevar motif --from-reference -f ref.fa \\\n"
    "                --reads both --proper-pair --max-insert-size 1000 \\\n"
    "                -q 30 -l 4 -o out.tsv  in1.bam in2.bam ...";

const std::string MotifCounterRunner::usage() { return MOTIF_USAGE; }

// Long-only option codes (above ASCII range to avoid collisions).
enum {
    OPT_READS = 1001,
    OPT_INCLUDE_ZERO,
    OPT_NO_INCLUDE_ZERO,
    OPT_FILENAME_HAS_SAMPLENAME,
    OPT_PROPER_PAIR,      // 1005
    OPT_MAX_INSERT_SIZE,  // 1006
    OPT_FROM_REFERENCE,   // 1007
};

MotifCounterRunner::MotifCounterRunner(int argc, char* argv[])
    : _args(new MotifArgs())
{
    if (argc < 2) {
        std::cout << MOTIF_USAGE << "\n" << std::endl;
        delete _args;
        _args = nullptr;
        std::exit(EXIT_SUCCESS);
    }

    static const struct option MOTIF_LOPTS[] = {
        {"output",          required_argument, NULL, 'o'},
        {"align-file-list", required_argument, NULL, 'L'},
        {"reference",       required_argument, NULL, 'f'},
        {"regions",         required_argument, NULL, 'r'},
        {"length",          required_argument, NULL, 'l'},
        {"mapq",            required_argument, NULL, 'q'},
        {"thread",          required_argument, NULL, 't'},
        {"reads",           required_argument, NULL, OPT_READS},
        {"include-zero",          no_argument, NULL, OPT_INCLUDE_ZERO},
        {"no-include-zero",       no_argument, NULL, OPT_NO_INCLUDE_ZERO},
        {"filename-has-samplename", no_argument, NULL, OPT_FILENAME_HAS_SAMPLENAME},
        {"proper-pair",             no_argument,       NULL, OPT_PROPER_PAIR},
        {"max-insert-size",   required_argument, NULL, OPT_MAX_INSERT_SIZE},
        {"from-reference",          no_argument,       NULL, OPT_FROM_REFERENCE},
        {"help",                  no_argument, NULL, 'h'},
        {0, 0, 0, 0}
    };

    // Save the full command line for traceability (mirrors variant_caller).
    _cmdline_string = "##basevar_motif_command=";
    for (int i = 0; i < argc; ++i) {
        _cmdline_string += (i > 0 ? " " : "") + std::string(argv[i]);
    }

    int c;
    std::vector<std::string> bv;
    optind = 1;  // reset getopt state for safe re-entry within tests
    while ((c = getopt_long(argc, argv, "o:L:f:r:l:q:t:h", MOTIF_LOPTS, nullptr)) >= 0) {
        std::stringstream ss(optarg ? optarg : "");
        switch (c) {
            case 'o': _args->output_file = optarg; break;
            case 'L':
                _args->in_filelist = optarg;
                bv = ngslib::get_firstcolumn_from_file(optarg);
                _args->input_bf.insert(_args->input_bf.end(), bv.begin(), bv.end());
                break;
            case 'f': _args->reference     = optarg; break;
            case 'r': _args->regions       = optarg; break;
            case 'l': ss >> _args->motif_length;     break;
            case 'q': ss >> _args->min_mapq;         break;
            case 't': ss >> _args->thread_num;       break;

            case OPT_READS:
                _args->read_end = parse_read_end(optarg);
                break;
            case OPT_INCLUDE_ZERO:           _args->include_zero            = true;  break;
            case OPT_NO_INCLUDE_ZERO:        _args->include_zero            = false; break;
            case OPT_FILENAME_HAS_SAMPLENAME:_args->filename_has_samplename = true;  break;
            case OPT_PROPER_PAIR:            _args->proper_pair             = true;  break;
            case OPT_MAX_INSERT_SIZE:        ss >> _args->max_insert_size;           break;
            case OPT_FROM_REFERENCE:         _args->from_reference          = true;  break;

            case 'h':
                std::cout << MOTIF_USAGE << std::endl;
                std::exit(EXIT_SUCCESS);

            default:
                std::cerr << "[ERROR] Unknown argument while parsing `basevar motif`.\n";
                std::exit(EXIT_FAILURE);
        }
    }

    // Collect positional BAM/CRAM files.
    while (optind < argc) {
        _args->input_bf.push_back(argv[optind++]);
    }

    // ---- validate ----
    if (_args->input_bf.empty()) {
        throw std::invalid_argument("[ERROR] No input BAM/CRAM files provided. "
                                    "Use -L <list> or supply files positionally.");
    }
    if (_args->motif_length < 1 || _args->motif_length > 10) {
        throw std::invalid_argument("[ERROR] -l/--length must be between 1 and 10 (got "
                                    + std::to_string(_args->motif_length) + ")");
    }
    if (_args->min_mapq < 0) {
        throw std::invalid_argument("[ERROR] -q/--mapq must be >= 0 (got "
                                    + std::to_string(_args->min_mapq) + ")");
    }
    if (_args->max_insert_size < 0) {
        throw std::invalid_argument("[ERROR] --max-insert-size must be >= 0 (got "
                                    + std::to_string(_args->max_insert_size) + ")");
    }
    if (_args->from_reference && _args->reference.empty()) {
        throw std::invalid_argument("[ERROR] --from-reference requires -f/--reference "
                                    "(a FASTA file with a .fai index).");
    }
    if (_args->thread_num < 1) _args->thread_num = 1;
}

MotifCounterRunner::~MotifCounterRunner() {
    if (_args) { delete _args; _args = nullptr; }
}

// ---- aggregate accessors ----

uint64_t MotifCounterRunner::total_reads() const {
    uint64_t v = 0; for (const auto& s : _results) v += s.total_reads;     return v;
}
uint64_t MotifCounterRunner::filtered_reads() const {
    uint64_t v = 0; for (const auto& s : _results) v += s.filtered_reads;  return v;
}
uint64_t MotifCounterRunner::used_reads() const {
    uint64_t v = 0; for (const auto& s : _results) v += s.used_reads;      return v;
}
uint64_t MotifCounterRunner::n_motifs_with_n() const {
    uint64_t v = 0; for (const auto& s : _results) v += s.n_motifs_with_n; return v;
}

void MotifCounterRunner::_initialize_motifs_for(SampleResult& s) const {
    enumerate_kmers(_args->motif_length, s.motif_counts);
}

// Return the first dot-delimited token of the stripped filename, e.g.
//   /path/to/SampleA.sort.bam  ->  "SampleA"
//   /path/to/SampleB.bam       ->  "SampleB"
static std::string filename_stem(const std::string& fn) {
    std::string base = ngslib::remove_filename_extension(ngslib::basename(fn));
    size_t dot = base.find('.');
    return (dot != std::string::npos && dot > 0) ? base.substr(0, dot) : base;
}

std::string MotifCounterRunner::_resolve_sample_id(const std::string& fn) const {
    if (_args->filename_has_samplename) {
        return filename_stem(fn);
    }
    // Prefer @RG SM tag.  Fall back to the filename stem if missing.
    try {
        ngslib::BamHeader bh(fn, _args->reference);
        std::string sm = bh.get_sample_name();
        if (!sm.empty()) return sm;
    } catch (const std::exception& e) {
        std::cerr << "[WARN] Cannot read @RG SM tag from " << fn
                  << " (" << e.what() << "); falling back to filename.\n";
    }
    return filename_stem(fn);
}

bool MotifCounterRunner::_passes_filters(const ngslib::BamRecord& r) const {
    if (!r.is_mapped())       return false;
    if (r.is_secondary())     return false;
    if (r.is_supplementary()) return false;
    if (r.is_duplicate())     return false;
    if (r.is_qc_fail())       return false;
    if (r.mapq() < _args->min_mapq) return false;

    // Read-end policy.
    //   * Paired-end reads: enforce the policy via FREAD1 / FREAD2 flags.
    //   * Single-end reads: have no R1/R2 designation.  We treat them as
    //     an implicit R1 (the fragment's 5' end is well-defined), so they
    //     are accepted under `R1` and `BOTH` and rejected under `R2`
    //     (asking for R2 of single-end data is a logical no-op).
    if (r.is_paired()) {
        switch (_args->read_end) {
            case ReadEnd::R1:   if (!r.is_read1()) return false; break;
            case ReadEnd::R2:   if (!r.is_read2()) return false; break;
            case ReadEnd::BOTH: /* accept either */              break;
        }
    } else {
        if (_args->read_end == ReadEnd::R2) return false;
    }

    // Proper-pair filter (PE only).  SE data is implicitly skipped via is_paired().
    if (r.is_paired() && _args->proper_pair && !r.is_proper_pair()) return false;

    // Insert-size cap filter (PE only).  insert_size() returns signed int64_t;
    // use std::abs to handle negative values (read2 orientation).
    if (r.is_paired() && _args->max_insert_size > 0
        && std::abs(r.insert_size()) > _args->max_insert_size) return false;

    return true;
}

void MotifCounterRunner::_process_one_file(SampleResult& s) const {
    ngslib::Bam bam(s.input_path, "r", _args->reference);

    // Per-worker Fasta instance: ngslib::Fasta is documented as NOT
    // thread-safe, so each worker thread owns its own loader.  Building
    // it only loads the .fai index (cheap), no sequence data is held in
    // memory between fetches.
    std::unique_ptr<ngslib::Fasta> fa_ptr;
    if (_args->from_reference) {
        fa_ptr.reset(new ngslib::Fasta(_args->reference));
    }
    ngslib::BamHeader& hdr = bam.header();
    const int k = _args->motif_length;

    std::vector<std::string> regions;
    if (!_args->regions.empty()) {
        ngslib::split(_args->regions, regions, ",");
    }

    auto consume_record = [&](const ngslib::BamRecord& r) {
        ++s.total_reads;
        if (!_passes_filters(r)) { ++s.filtered_reads; return; }

        std::string motif = _args->from_reference
            ? extract_5p_motif_from_reference(r, *fa_ptr, hdr, k)
            : extract_5p_motif(r, k);
        if (motif.empty()) { ++s.filtered_reads; return; }

        // Any motif containing a non-ACGT base (N or IUPAC ambiguity code)
        // is excluded -- see is_pure_acgt() for why find('N') is not enough.
        if (!is_pure_acgt(motif)) {
            ++s.n_motifs_with_n;
            return;
        }
        ++s.motif_counts[motif];
        ++s.used_reads;
    };

    ngslib::BamRecord r;
    if (regions.empty()) {
        while (bam.next(r) >= 0) consume_record(r);
    } else {
        // Deduplicate reads across potentially overlapping regions.
        // Key = qname + '/1' or '/2' (paired) or '/0' (single-end).
        // Memory cost is O(unique reads in all queried regions), which is
        // acceptable for typical cfDNA targeted-region workflows.
        std::unordered_set<std::string> seen_reads;
        for (const auto& reg : regions) {
            if (reg.empty()) continue;
            bam.fetch(reg);
            while (bam.next(r) >= 0) {
                std::string key = r.qname();
                if (r.is_paired()) {
                    key += r.is_read1() ? "/1" : "/2";
                } else {
                    key += "/0";
                }
                if (!seen_reads.insert(key).second) {
                    // Already processed this read from a previous (overlapping) region.
                    continue;
                }
                consume_record(r);
            }
        }
    }
}

void MotifCounterRunner::_write_tsv(const std::string& path) const {
    std::ofstream ofs(path);
    if (!ofs) {
        throw std::runtime_error("[ERROR] Cannot open output TSV for writing: " + path);
    }
    ofs << "#sample\tmotif\tcount\tfrequency\n";
    ofs << std::fixed << std::setprecision(6);

    for (const auto& s : _results) {
        // Motifs containing N are excluded from numerator and used_reads.
        const uint64_t denom = s.used_reads;
        for (const auto& kv : s.motif_counts) {
            if (!_args->include_zero && kv.second == 0) continue;
            double freq = (denom > 0) ? static_cast<double>(kv.second) / static_cast<double>(denom) : 0.0;
            ofs << s.sample_id << "\t" << kv.first << "\t" << kv.second << "\t" << freq << "\n";
        }
    }
    ofs.close();
}

void MotifCounterRunner::_print_summary(std::ostream& os) const {
    os << "\n=== cfDNA End-Motif Statistics ===\n"
       << "Motif length (k):    " << _args->motif_length << "\n"
       << "Read end policy:     " << read_end_name(_args->read_end) << "\n"
       << "Motif source:        " << (_args->from_reference ? "reference genome" : "read sequence") << "\n"
       << "Worker threads:      " << _args->thread_num << "\n"
       << "Samples processed:   " << _results.size() << "\n";

    // Aggregate row.
    os << "\n-- Aggregate over all samples --\n"
       << "Total reads scanned:      " << total_reads()      << "\n"
       << "Filtered reads:           " << filtered_reads()   << "\n"
       << "Motifs w/ non-ACGT bases: " << n_motifs_with_n()  << " (excluded)\n"
       << "Valid reads counted:      " << used_reads()       << "\n";

    if (used_reads() == 0) {
        os << "[WARN] No reads passed the filters; nothing to report.\n";
        return;
    }

    // MDS = normalised Shannon entropy over the 4^k motif distribution.
    //   MDS = -sum(p_i * log2(p_i)) / log2(4^k)
    //       = -sum(p_i * log2(p_i)) / (2*k)
    // Range [0, 1]: 1 = perfectly uniform, 0 = all reads share one motif.
    auto compute_mds = [&](const SampleResult& s) -> double {
        if (s.used_reads == 0) return 0.0;
        const double denom = static_cast<double>(s.used_reads);
        const double norm  = 2.0 * _args->motif_length;  // log2(4^k) = 2k
        double entropy = 0.0;
        for (const auto& kv : s.motif_counts) {
            if (kv.second == 0) continue;
            double p = static_cast<double>(kv.second) / denom;
            entropy -= p * std::log2(p);
        }
        return entropy / norm;
    };

    // Per-sample one-line summary.
    // Column widths: 24 + 4*14 + 10 = 90 characters.
    os << "\n-- Per-sample summary --\n"
       << std::left
       << std::setw(24) << "Sample"
       << std::setw(14) << "Total"
       << std::setw(14) << "Filtered"
       << std::setw(14) << "Used"
       << std::setw(14) << "N-motifs"
       << std::setw(10) << "MDS" << "\n"
       << "------------------------------------------------------------------------------------------\n";
    for (const auto& s : _results) {
        os << std::left
           << std::setw(24) << s.sample_id
           << std::setw(14) << s.total_reads
           << std::setw(14) << s.filtered_reads
           << std::setw(14) << s.used_reads
           << std::setw(14) << s.n_motifs_with_n
           << std::right << std::fixed << std::setprecision(6)
           << std::setw(10) << compute_mds(s) << "\n";
        os.unsetf(std::ios::fixed);  // prevent fixed-point bleed into next iteration
    }
}

int MotifCounterRunner::run() {
    // Validate the output path up-front so a bad -o does not waste a long run.
    if (!_args->output_file.empty()) {
        std::ofstream probe(_args->output_file);
        if (!probe) {
            throw std::runtime_error("[ERROR] Cannot open output TSV for writing: "
                                     + _args->output_file);
        }
        probe.close();
    }

    const size_t n_files = _args->input_bf.size();
    int n_threads = std::min<int>(_args->thread_num, static_cast<int>(n_files));
    if (n_threads < 1) n_threads = 1;

    std::cerr << "[INFO] basevar motif v" << BASEVAR_VERSION
              << " - counting " << _args->motif_length << "-mer cfDNA end-motifs from "
              << n_files << " BAM/CRAM file(s) using " << n_threads << " worker thread(s).\n";

    // Stage 1: resolve sample IDs (sequential; cheap and gives ordered output).
    _results.resize(n_files);
    for (size_t i = 0; i < n_files; ++i) {
        _results[i].input_path = _args->input_bf[i];
        _results[i].sample_id  = _resolve_sample_id(_args->input_bf[i]);
        _initialize_motifs_for(_results[i]);
    }

    // Stage 2: process files concurrently.  Each task touches only its own
    // SampleResult, so no locking is required.
    if (n_threads <= 1 || n_files <= 1) {
        // Single-threaded path: avoids ThreadPool overhead and keeps stack
        // traces simple when the user runs with -t 1.
        for (size_t i = 0; i < n_files; ++i) {
            std::cerr << "[INFO] [" << (i + 1) << "/" << n_files << "] "
                      << _results[i].sample_id << "  <-  " << _results[i].input_path << "\n";
            _process_one_file(_results[i]);
        }
    } else {
        ThreadPool pool(n_threads);
        std::vector<std::future<void>> futures;
        futures.reserve(n_files);
        for (size_t i = 0; i < n_files; ++i) {
            futures.emplace_back(pool.submit(
                [this, i]() {
                    std::cerr << "[INFO] processing sample "
                              << _results[i].sample_id << "  <-  "
                              << _results[i].input_path << "\n";
                    _process_one_file(_results[i]);
                }));
        }
        for (auto& f : futures) f.get();
        if (pool.has_error()) {
            std::cerr << "[ERROR] one or more worker threads failed.\n";
            return 1;
        }
    }

    _print_summary(std::cout);
    if (!_args->output_file.empty()) {
        _write_tsv(_args->output_file);
        std::cerr << "[INFO] Wrote per-sample TSV to: " << _args->output_file << "\n";
    }
    return 0;
}

// ============================================================
// C-style entry point used by main.cpp dispatcher.
// ============================================================

int motif_counter_runner(int argc, char* argv[]) {
    try {
        MotifCounterRunner runner(argc, argv);
        return runner.run();
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }
}

}  // namespace motif
}  // namespace basevar
