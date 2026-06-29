/**
 * @file fetal_fraction.cpp
 *
 * @brief Implementation of the `basevar fetalfrac` (hidden) subcommand.
 *
 *   Hidden = intentionally NOT advertised in `basevar --help`.  This is
 *   a private utility for NIPT cfDNA male fetal-fraction estimation;
 *   it can still be invoked directly as `basevar fetalfrac ...`.
 *
 *   ## Scientific rationale (chrY-counting)
 *
 *   For a male fetus, the cfDNA pool in maternal plasma is a mixture:
 *
 *      mother (XX) at fraction (1 - p)   +   fetus (XY) at fraction p
 *
 *   Per molecule, the mother contributes 0 chrY copies and 2 autosomal
 *   copies; the fetus contributes 1 chrY and 2 autosomal copies.  Read
 *   density is proportional to molecule count times sequence length, so:
 *
 *      Y_reads_density   ~  p * L_chrY_uniq
 *      Auto_reads_density ~  2 * L_autosome
 *
 *   Therefore (after dropping common scale factors):
 *
 *      Y_ratio = Y / Auto = p * L_chrY_uniq / (2 * L_autosome)
 *      p       = Y_ratio * (2 * L_autosome / L_chrY_uniq)
 *              = Y_ratio * SCALE
 *
 *   For GRCh38 the unique-mappable chrY is ~28 Mb out of ~57 Mb total,
 *   and autosomes are ~2.87 Gb, so SCALE ~ 100..200 depending on which
 *   filters are applied.  The empirically-correct SCALE MUST be
 *   calibrated from the user's own male/female reference samples.
 *
 *   ## Female-fetus refusal
 *
 *   Under low-coverage WGS (lpWGS) the female fetal fraction CANNOT be
 *   measured directly via Y-counting (there is no fetal Y signal to
 *   count).  Indirect estimators (fragment-size, X-skewing) have high
 *   variance and lack molecular precision compared to SNP-based methods.
 *   To prevent "false precision", this tool deliberately REFUSES to
 *   output a female FF.  The user gets a FEMALE call and a FF of -1.
 *
 *   ## Bugs that existed in the original ff_calculator.cpp (now fixed)
 *
 *     1. The original simulated reads with std::random_device and never
 *        opened a real BAM file.  Now reads are pulled through
 *        ngslib::Bam, exactly like `basevar caller`.
 *     2. PAR1/PAR2 (pseudoautosomal regions) were taken from a per-read
 *        flag that nothing populated.  Now we apply built-in
 *        GRCh38/GRCh37 PAR coordinates (or a user-supplied BED) to
 *        actually exclude PAR positions on chrY.
 *     3. The previous default scaling factor (kappa=220) was empirically
 *        wrong for GRCh38 and not documented.  Now the default (100) is
 *        derived from genome geometry and is documented in the help text.
 *     4. There was no proper QC (secondary / supplementary / duplicate /
 *        proper-pair / insert-size).  Now the same conservative filter
 *        set used by `basevar motif` is applied.
 *     5. The autosome detector accepted contigs whose names happened
 *        not to contain X / Y / M.  Now `_is_autosome` whitelists
 *        chr1..chr22 / 1..22 and explicitly rejects ALT/decoy/HLA.
 *     6. The original program had its own main() and would have caused
 *        a duplicate-symbol link error against main.cpp.  Now the
 *        module exposes only `fetal_fraction_runner(argc, argv)` and
 *        the `FetalFractionRunner` class.
 *
 * @author Shujia Huang
 * @date   2026-05-27
 */
#include "fetal_fraction.h"

#include <getopt.h>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "external/thread_pool.h"
#include "io/bam.h"
#include "io/bam_header.h"
#include "version.h"

namespace basevar {
namespace fetalfrac {

// ============================================================
// BedIntervalIndex
// ============================================================
//
// Internal storage uses ngslib::GenomeRegion with 1-based, inclusive
// coordinates -- the same convention used everywhere else in BaseVar
// (see pipeline::parse_region_string).  load_bed() converts on-disk
// BED 0-based half-open to 1-based inclusive once, at the boundary.

void BedIntervalIndex::_merge(std::vector<ngslib::GenomeRegion>& v) {
    if (v.empty()) return;
    std::sort(v.begin(), v.end(),
              [](const ngslib::GenomeRegion& a, const ngslib::GenomeRegion& b) {
                  return a.start < b.start;
              });
    std::vector<ngslib::GenomeRegion> out;
    ngslib::GenomeRegion cur = v.front();
    for (size_t i = 1; i < v.size(); ++i) {
        // 1-based inclusive overlap: [a,b] and [c,d] overlap iff c <= b.
        if (v[i].start <= cur.end) {
            cur.end = std::max(cur.end, v[i].end);
        } else {
            out.push_back(cur);
            cur = v[i];
        }
    }
    out.push_back(cur);
    v.swap(out);
}

void BedIntervalIndex::load_bed(const std::string& filepath) {
    if (!ngslib::is_readable(filepath)) {
        throw std::runtime_error("[ERROR] Cannot open BED file: " + filepath);
    }
    std::ifstream ifs(filepath);
    if (!ifs.is_open()) {
        throw std::runtime_error("[ERROR] Cannot open BED file: " + filepath);
    }
    std::string line;
    std::vector<std::string> fields;
    while (std::getline(ifs, line)) {
        // Strip a trailing '\r' so files with CRLF line endings still parse.
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        // UCSC browser / track header lines (proper prefix match, not
        // just any line that happens to start with 't' / 'b').
        if (line.compare(0, 6, "track ")   == 0 ||
            line.compare(0, 8, "browser ") == 0) continue;

        // BED is tab-delimited; reuse BaseVar's standard tokenizer.
        ngslib::split(line, fields, "\t");
        if (fields.size() < 3) continue;          // need chrom, start, end
        if (fields[0].empty()) continue;

        // BED uses 0-based half-open coordinates [s, e).
        long long s = 0, e = 0;
        try {
            s = std::stoll(fields[1]);
            e = std::stoll(fields[2]);
        } catch (const std::exception&) {
            continue;                              // malformed numeric column
        }
        if (s < 0 || e <= s) continue;             // skip empty / inverted

        // [s, e) (0-based half-open) -> [s+1, e] (1-based inclusive).
        const uint32_t start_1b = static_cast<uint32_t>(s) + 1u;
        const uint32_t end_1b   = static_cast<uint32_t>(e);
        _regions[fields[0]].emplace_back(fields[0], start_1b, end_1b);
    }
    finalize();
}

void BedIntervalIndex::add_interval(const std::string& chrom,
                                    uint32_t start, uint32_t end) {
    // 1-based inclusive: [start, end] with start <= end.
    if (end < start) return;
    _regions[chrom].emplace_back(chrom, start, end);
}

void BedIntervalIndex::finalize() {
    for (auto& kv : _regions) _merge(kv.second);
    _is_loaded = !_regions.empty();
}

bool BedIntervalIndex::contains(const std::string& chrom, uint32_t pos) const {
    auto it = _regions.find(chrom);
    if (it == _regions.end()) return false;
    const auto& v = it->second;
    // Binary-search for the last interval with start <= pos.
    auto ub = std::upper_bound(
        v.begin(), v.end(), pos,
        [](uint32_t value, const ngslib::GenomeRegion& g) {
            return value < g.start;
        });
    if (ub == v.begin()) return false;
    --ub;
    // 1-based inclusive containment.
    return pos >= ub->start && pos <= ub->end;
}

uint64_t BedIntervalIndex::total_length(
    const std::function<bool(const std::string&)>& pred) const
{
    uint64_t total = 0;
    for (const auto& kv : _regions) {
        if (!pred(kv.first)) continue;
        for (const auto& g : kv.second) {
            // GenomeRegion::size() returns 1-based inclusive length.
            total += static_cast<uint64_t>(g.size());
        }
    }
    return total;
}

// ----------------------------------------------------------------
// Helper: total bp of intersection between two BedIntervalIndexes,
// restricted to chrY (we use this to subtract the PAR overlap from
// the chrY mappability mask when auto-computing --scale).
// Both inputs are sorted+merged inside their own classes.
// ----------------------------------------------------------------
static uint64_t intersection_length_chrY(const BedIntervalIndex& a,
                                         const BedIntervalIndex& b) {
    uint64_t total = 0;
    for (const auto& kv : a.regions()) {
        if (!is_y_chrom(kv.first)) continue;
        auto bit = b.regions().find(kv.first);
        if (bit == b.regions().end()) continue;
        const auto& av = kv.second;
        const auto& bv = bit->second;
        size_t i = 0, j = 0;
        while (i < av.size() && j < bv.size()) {
            const uint32_t lo = std::max(av[i].start, bv[j].start);
            const uint32_t hi = std::min(av[i].end, bv[j].end);
            if (lo <= hi) total += static_cast<uint64_t>(hi - lo + 1);
            if (av[i].end < bv[j].end) ++i; else ++j;
        }
    }
    return total;
}

// ============================================================
// Public chromosome-name helpers (also used by unit tests).
// ============================================================

// Strip a leading "chr" prefix (case-insensitive).
std::string strip_chr_prefix(const std::string& chrom) {
    if (chrom.size() >= 3 &&
        (chrom[0] == 'c' || chrom[0] == 'C') &&
        (chrom[1] == 'h' || chrom[1] == 'H') &&
        (chrom[2] == 'r' || chrom[2] == 'R')) {
        return chrom.substr(3);
    }
    return chrom;
}

bool is_y_chrom(const std::string& chrom) {
    std::string s = strip_chr_prefix(chrom);
    return s == "Y" || s == "y";
}

bool is_autosome(const std::string& chrom) {
    std::string s = strip_chr_prefix(chrom);
    if (s.empty()) return false;
    // Reject anything other than 1..22 (no X/Y/M/MT/EBV/ALT/decoy/HLA/random).
    for (char c : s) if (!std::isdigit(static_cast<unsigned char>(c))) return false;
    int n = std::atoi(s.c_str());
    return (n >= 1 && n <= 22);
}

// Return the first dot-delimited token of the stripped filename, e.g.
//   /path/to/SampleA.sort.bam  ->  "SampleA"
//   /path/to/SampleB.bam       ->  "SampleB"
static std::string filename_stem(const std::string& fn) {
    std::string base = ngslib::remove_filename_extension(ngslib::basename(fn));
    size_t dot = base.find('.');
    return (dot != std::string::npos && dot > 0) ? base.substr(0, dot) : base;
}

static const char* sex_name(FetalSex s) {
    switch (s) {
        case FetalSex::MALE:         return "MALE";
        case FetalSex::FEMALE:       return "FEMALE";
        case FetalSex::UNDETERMINED: return "UNDETERMINED";
    }
    return "?";
}

static GenomeBuild parse_genome_build(const std::string& s) {
    std::string v;
    v.reserve(s.size());
    for (char c : s) v.push_back(std::tolower(static_cast<unsigned char>(c)));
    if (v == "hg38" || v == "grch38") return GenomeBuild::GRCh38;
    if (v == "hg19" || v == "grch37") return GenomeBuild::GRCh37;
    if (v == "none" || v == "off")    return GenomeBuild::NONE;
    throw std::invalid_argument("[ERROR] --build must be one of: hg38, hg19, none (got '" + s + "')");
}

static const char* genome_build_name(GenomeBuild b) {
    switch (b) {
        case GenomeBuild::GRCh38: return "hg38";
        case GenomeBuild::GRCh37: return "hg19";
        case GenomeBuild::NONE:   return "none";
    }
    return "?";
}

// ============================================================
// Help text
// ============================================================

static const std::string FF_USAGE =
    "About: Estimate male fetal fraction (FF) from NIPT cfDNA BAM/CRAM data\n"
    "       using chrY read counts. HIDDEN BaseVar subcommand for private use.\n"
    "       Female fetuses are detected and reported, but FF is NOT computed for\n"
    "       them under lpWGS (use a SNP-based assay instead).\n"
    "Usage: basevar fetalfrac [options] <-o output.tsv> [-L bam.list] in1.bam [in2.bam ...]\n\n"

    "Required arguments:\n"
    "  -o, --output FILE            Output TSV file (sample, sex, FF, counters).\n\n"

    "Optional arguments:\n"
    "  -L, --align-file-list FILE   BAM/CRAM files list, one path per row.\n"
    "  -f, --reference FILE         Reference FASTA (required for CRAM).\n"
    "  -r, --regions REG[,...]      Restrict counting to these regions, comma-separated.\n"
    "                               Formats: chr | chr:start | chr:start-end\n"
    "  -B, --mappability-bed FILE   Genome-wide inclusion-mask BED.  When supplied,\n"
    "                               ONLY reads whose 5'-mapped position falls inside\n"
    "                               an interval are counted -- this filter is applied\n"
    "                               to BOTH the chrY numerator AND the chr1..22\n"
    "                               denominator BEFORE the autosome/Y split.\n"
    "                               The BED MUST therefore cover both the autosomal\n"
    "                               regions and the chrY regions you want to count;\n"
    "                               a chrY-only BED zeros out the denominator and\n"
    "                               yields UNDETERMINED.  Whenever -B is added or\n"
    "                               changed, --scale MUST be re-calibrated against\n"
    "                               male reference samples processed with the same\n"
    "                               BED (the implied 2*L_auto/L_chrY ratio changes).\n"
    "                               Strongly recommended for chrY, which is dominated\n"
    "                               by repeats / palindromes.\n"
    "      --par-bed FILE           Custom PAR exclusion BED (overrides --build).\n"
    "      --build {hg38|hg19|none} Built-in PAR coordinates [hg38].\n"
    "                               hg38: PAR1 chrY:10001-2781479,\n"
    "                                     PAR2 chrY:56887903-57217415\n"
    "                               hg19: PAR1 chrY:10001-2649520,\n"
    "                                     PAR2 chrY:59034050-59363566\n"
    "                               Use 'none' to disable PAR exclusion entirely.\n"
    "  -q, --mapq INT               Minimum MAPQ to keep a read [30].\n"
    "  -t, --thread INT             Number of worker threads (one file per thread)\n"
    "                               [hardware_concurrency].\n"
    "      --proper-pair            Only count reads flagged BAM_FPROPER_PAIR.\n"
    "                               Silently ignored for SE data. [off]\n"
    "      --max-insert-size INT    Discard reads whose |insert size| > INT.\n"
    "                               0 = no limit.  Lo-lab cfDNA pipelines typically\n"
    "                               use 1000 to drop chimeric / discordant reads.\n"
    "                               Silently ignored for SE data. [0]\n"
    "      --scale FLOAT            FF scaling factor:\n"
    "                                   FF = (Y/Auto - noise) * scale.\n"
    "                               Conceptually scale ~ 2 * L_autosome / L_chrY_uniq.\n"
    "                               If --scale is NOT supplied, it is auto-computed\n"
    "                               from genome geometry:\n"
    "                                 -B given: scale = 2 * (sum of BED autosomal bp)\n"
    "                                                  / (sum of BED chrY bp - chrY n PAR).\n"
    "                                 no -B   : scale is the canonical GRC value for\n"
    "                                           --build (~106.2 for hg38, ~102.2 for\n"
    "                                           hg19); --build none falls back to 100.\n"
    "                               The auto value is only a uniform-coverage geometric\n"
    "                               approximation -- recalibrate against male reference\n"
    "                               samples processed with the SAME -B / --build for\n"
    "                               clinical reporting.  Default: auto.\n"
    "      --noise FLOAT            Background Y/Auto ratio to subtract before\n"
    "                               scaling (e.g. index-hopping floor) [0.0].\n"
    "                               There is NO universal best value -- the noise\n"
    "                               floor is platform / library / pipeline specific\n"
    "                               and MUST be calibrated empirically against >=30\n"
    "                               confirmed-female-fetus samples processed with\n"
    "                               the SAME -B / --build / -q / library prep /\n"
    "                               sequencer.  Set --noise to the median (or 95th\n"
    "                               percentile) of the y_ratio observed in that\n"
    "                               cohort.  Reference magnitudes from the literature\n"
    "                               (post-PAR exclusion, MAPQ >= 30):\n"
    "                                 MGI DNBSEQ / BGISEQ      : 1e-5 .. 5e-5\n"
    "                                 NovaSeq + UDI            : 1e-5 .. 5e-5\n"
    "                                 HiSeq X / non-UDI dual   : 5e-5 .. 1e-4\n"
    "                                 HiSeq 2500 single-index  : 1e-4 .. 5e-4\n"
    "                                 No calibration cohort    : keep 0 (default).\n"
    "                               MGI DNBSEQ (DNBSEQ-T7 / G400, MGISEQ-2000,\n"
    "                               BGISEQ-500) uses DNB rolling-circle amplification,\n"
    "                               which is NOT susceptible to ExAmp-mediated index\n"
    "                               hopping; with proper dual indices its residual Y\n"
    "                               signal is at or below NovaSeq + UDI -- this is\n"
    "                               also the primary platform of the NIFTY cohort.\n"
    "                               Sources of residual Y signal: chrY mis-mapping\n"
    "                               (X->Y, repeats / palindromes), index hopping on\n"
    "                               patterned flowcells, and rare maternal Y\n"
    "                               mosaicism.  References: Lun 2008 PNAS; Chiu 2008\n"
    "                               PNAS / 2011 BMJ; Kim 2015 Prenat Diagn (SeqFF);\n"
    "                               Sinha 2017 bioRxiv (index hopping); Costello\n"
    "                               2018 BMC Genomics (UDI remediation); Korostin\n"
    "                               2020 PLoS One / Foox 2021 Nat Biotechnol (MGI vs\n"
    "                               Illumina); Liu 2024 Cell Genomics (NIFTY on\n"
    "                               BGISEQ/MGI). See Notes below for the coupling\n"
    "                               with --male-threshold.\n"
    "      --male-threshold FLOAT   If Y/Auto < this, the sample is called FEMALE\n"
    "                               and FF is NOT reported.  [1e-4]\n"
    "      --filename-has-samplename\n"
    "                               Derive sample IDs from filenames instead of @RG SM.\n"
    "                               E.g. /path/SampleA.bam -> SampleA.\n"
    "      --calibrate              Calibration mode: treat all input samples as known\n"
    "                               pure male (FF=1.0), compute empirical scale as\n"
    "                               1 / mean(y_ratio - noise), and report the\n"
    "                               calibrated value.  Use this scale for subsequent\n"
    "                               --scale <value> runs.\n"
    "  -h, --help                   Show this help message and exit.\n\n"

    "Output TSV columns (one row per input BAM/CRAM):\n"
    "  #sample           Sample ID (from @RG SM tag, or filename).\n"
    "  sex               MALE / FEMALE / UNDETERMINED.\n"
    "  fetal_fraction    Decimal fraction (0..1) for MALE; -1 for FEMALE/UNDETERMINED.\n"
    "  y_ratio           valid_y / valid_autosomal (raw Y/Auto ratio).\n"
    "  valid_autosomal   chr1..chr22 reads passing all filters (denominator).\n"
    "  valid_y           chrY reads outside PAR passing all filters (numerator).\n"
    "  y_par_excluded    chrY reads dropped because they fell inside PAR1/PAR2.\n"
    "  total_reads       All reads scanned across the input.\n"
    "  filtered_reads    Reads rejected by QC / region / mappability filters,\n"
    "                    or from non-informative contigs (chrX/chrM/ALT).\n"
    "  scale             Scaling factor used for this sample.\n"
    "  noise             Y/Auto noise floor used for this sample.\n"
    "  male_threshold    Y/Auto threshold used to call MALE vs FEMALE.\n\n"

    "Examples:\n"
    "  # Single BAM, GRCh38, default scaling.  Output TSV + stdout summary.\n"
    "  basevar fetalfrac -o ff.tsv  sample.bam\n\n"
    "  # Multi-sample run with mappability mask and custom calibration.\n"
    "  basevar fetalfrac -t 8 -B chrY.unique.bed --build hg38 \\\n"
    "                    --scale 105.3 --noise 1e-5 \\\n"
    "                    -o cohort.ff.tsv -L bam.list\n\n"
    "  # CRAM input.  -f is mandatory.\n"
    "  basevar fetalfrac -f hg38.fa -o ff.tsv  sample.cram\n\n"
    "  # Calibrate scale from known pure male samples.\n"
    "  basevar fetalfrac --calibrate -f hg38.fa --noise 1e-5 \\\n"
    "                    -o calibrate.tsv -L pure_male.list\n\n"

    "Notes:\n"
    "  - Autosomes are restricted to chr1..chr22 / 1..22; ALT/decoy/HLA contigs\n"
    "    are NOT counted.\n"
    "  - chrY reads inside PAR1/PAR2 are counted SEPARATELY (`y_par_excluded`)\n"
    "    and never enter the numerator -- they are shared with chrX and would\n"
    "    inflate FF for female fetuses.\n"
    "  - Apply `--proper-pair` and `--max-insert-size 1000` for cfDNA-grade\n"
    "    QC; `-q 30` is the conservative MAPQ default.\n"
    "  - `--scale` is auto-computed from genome geometry when not supplied:\n"
    "      no -B + --build hg38 -> ~106.23\n"
    "      no -B + --build hg19 -> ~102.16\n"
    "      no -B + --build none -> fallback 100.0 (with a [WARN] message).\n"
    "      -B given             -> 2 * BED_autosome_bp / (BED_chrY_bp - BED chrY n PAR).\n"
    "    The auto value is a uniform-coverage geometric approximation; for clinical\n"
    "    use, calibrate `--scale` empirically against male reference samples\n"
    "    processed with the SAME -B / --build / region settings. No reference\n"
    "    length is inferred from the BAM at runtime; only the BED and the\n"
    "    canonical GRC assembly report are used.\n"
    "  - The pair (-B, --scale) is a tightly coupled calibration unit. When you\n"
    "    let --scale auto, the program tracks -B for you.  When you pin --scale\n"
    "    explicitly, you MUST keep -B unchanged across calibration and production.\n"
    "  - --noise calibration procedure (REQUIRED for clinical use):\n"
    "      1. Run >=30 confirmed-female-fetus samples through the SAME pipeline\n"
    "         (same -B / --build / -q / --proper-pair / --max-insert-size /\n"
    "         library prep / sequencer model).\n"
    "      2. From the output TSV, take the empirical y_ratio distribution of\n"
    "         those female samples.\n"
    "      3. Set --noise to the median (or 95th percentile) of that distribution.\n"
    "      4. Set --male-threshold to clearly above the noise floor (>= 5x the\n"
    "         chosen --noise) so the MALE/FEMALE call is robust.\n"
    "      5. Recalibrate whenever ANY of the inputs in step 1 change.\n"
    "    Reference orders of magnitude: MGI DNBSEQ / BGISEQ ~ 1e-5..5e-5 (DNB\n"
    "    chemistry is immune to ExAmp index hopping and is the NIFTY-cohort\n"
    "    primary platform); NovaSeq + UDI ~ 1e-5..5e-5; HiSeq X / non-UDI ~\n"
    "    5e-5..1e-4; HiSeq 2500 single-index ~ 1e-4..5e-4.\n"
    "    Without a calibration cohort, leave --noise at 0 (the default never\n"
    "    silently biases FF; raw y_ratio remains visible in the TSV).\n"
    "  - --noise and --male-threshold are coupled: if --noise >= --male-threshold,\n"
    "    any sample with y_ratio above the threshold but below the noise floor\n"
    "    will be called MALE with FF=0 (a discontinuity). Keep --noise <\n"
    "    --male-threshold / 5 to avoid this edge case.";

const std::string FetalFractionRunner::usage() { return FF_USAGE; }

// Long-only option codes (above ASCII range to avoid collisions).
enum {
    OPT_PROPER_PAIR = 1001,
    OPT_MAX_INSERT_SIZE,
    OPT_SCALE,
    OPT_NOISE,
    OPT_MALE_THRESHOLD,
    OPT_FILENAME_HAS_SAMPLENAME,
    OPT_PAR_BED,
    OPT_BUILD,
    OPT_CALIBRATE,
};

// ============================================================
// FetalFractionRunner constructor / destructor
// ============================================================

FetalFractionRunner::FetalFractionRunner(int argc, char* argv[])
    : _args(new FetalFractionArgs())
{
    if (argc < 2) {
        // Same convention as motif_counter: print usage and exit cleanly.
        std::cout << FF_USAGE << "\n" << std::endl;
        delete _args;
        _args = nullptr;
        std::exit(EXIT_SUCCESS);
    }

    static const struct option FF_LOPTS[] = {
        {"output",                  required_argument, NULL, 'o'},
        {"align-file-list",         required_argument, NULL, 'L'},
        {"reference",               required_argument, NULL, 'f'},
        {"regions",                 required_argument, NULL, 'r'},
        {"mappability-bed",         required_argument, NULL, 'B'},
        {"mapq",                    required_argument, NULL, 'q'},
        {"thread",                  required_argument, NULL, 't'},
        {"par-bed",                 required_argument, NULL, OPT_PAR_BED},
        {"build",                   required_argument, NULL, OPT_BUILD},
        {"proper-pair",             no_argument,       NULL, OPT_PROPER_PAIR},
        {"max-insert-size",         required_argument, NULL, OPT_MAX_INSERT_SIZE},
        {"scale",                   required_argument, NULL, OPT_SCALE},
        {"noise",                   required_argument, NULL, OPT_NOISE},
        {"male-threshold",          required_argument, NULL, OPT_MALE_THRESHOLD},
        {"filename-has-samplename", no_argument,       NULL, OPT_FILENAME_HAS_SAMPLENAME},
        {"calibrate",               no_argument,       NULL, OPT_CALIBRATE},
        {"help",                    no_argument,       NULL, 'h'},
        {0, 0, 0, 0}
    };

    // Save the full command line for traceability.
    _cmdline_string = "##basevar_fetalfrac_command=";
    for (int i = 0; i < argc; ++i) {
        _cmdline_string += (i > 0 ? " " : "") + std::string(argv[i]);
    }

    int c;
    std::vector<std::string> bv;
    optind = 1;  // Reset getopt state for safe re-entry within tests.
    while ((c = getopt_long(argc, argv, "o:L:f:r:B:q:t:h", FF_LOPTS, nullptr)) >= 0) {
        std::stringstream ss(optarg ? optarg : "");
        switch (c) {
            case 'o': _args->output_file = optarg; break;
            case 'L':
                _args->in_filelist = optarg;
                bv = ngslib::get_firstcolumn_from_file(optarg);
                _args->input_bf.insert(_args->input_bf.end(), bv.begin(), bv.end());
                break;
            case 'f': _args->reference        = optarg; break;
            case 'r': _args->regions          = optarg; break;
            case 'B': _args->mappability_bed  = optarg; break;
            case 'q': ss >> _args->min_mapq;            break;
            case 't': ss >> _args->thread_num;          break;

            case OPT_PAR_BED:
                _args->par_bed = optarg;
                break;
            case OPT_BUILD:
                _args->genome_build = parse_genome_build(optarg);
                break;
            case OPT_PROPER_PAIR:
                _args->proper_pair = true;
                break;
            case OPT_MAX_INSERT_SIZE:
                ss >> _args->max_insert_size;
                break;
            case OPT_SCALE:
                ss >> _args->scaling_factor;
                _args->scale_user_set = true;
                break;
            case OPT_NOISE:
                ss >> _args->y_noise;
                break;
            case OPT_MALE_THRESHOLD:
                ss >> _args->male_threshold;
                break;
            case OPT_FILENAME_HAS_SAMPLENAME:
                _args->filename_has_samplename = true;
                break;
            case OPT_CALIBRATE:
                _args->calibrate = true;
                break;

            case 'h':
                std::cout << FF_USAGE << std::endl;
                std::exit(EXIT_SUCCESS);

            default:
                std::cerr << "[ERROR] Unknown argument while parsing `basevar fetalfrac`.\n";
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
    if (_args->min_mapq < 0) {
        throw std::invalid_argument("[ERROR] -q/--mapq must be >= 0 (got "
                                    + std::to_string(_args->min_mapq) + ")");
    }
    if (_args->max_insert_size < 0) {
        throw std::invalid_argument("[ERROR] --max-insert-size must be >= 0 (got "
                                    + std::to_string(_args->max_insert_size) + ")");
    }
    // --scale validation only applies to explicit user input.  When the
    // user does not pass --scale we auto-compute it later, after the
    // mappability / PAR indexes are loaded.
    if (_args->scale_user_set && !(_args->scaling_factor > 0.0)) {
        throw std::invalid_argument("[ERROR] --scale must be > 0 (got "
                                    + std::to_string(_args->scaling_factor) + ")");
    }
    if (_args->y_noise < 0.0) {
        throw std::invalid_argument("[ERROR] --noise must be >= 0 (got "
                                    + std::to_string(_args->y_noise) + ")");
    }
    if (_args->male_threshold < 0.0) {
        throw std::invalid_argument("[ERROR] --male-threshold must be >= 0 (got "
                                    + std::to_string(_args->male_threshold) + ")");
    }
    if (_args->thread_num < 1) _args->thread_num = 1;

    // ---- build PAR exclusion index ----
    // Priority: user --par-bed > built-in --build > none.
    if (!_args->par_bed.empty()) {
        _par_regions.load_bed(_args->par_bed);
    } else {
        _load_builtin_par();
    }

    // ---- optional mappability inclusion BED ----
    if (!_args->mappability_bed.empty()) {
        _mappability.load_bed(_args->mappability_bed);
    }

    // ---- auto-compute --scale from genome geometry ----
    // The (-B, --scale) pair is a tightly coupled calibration unit.
    // When the user does not supply --scale, derive it from genome
    // geometry so that adding / removing -B does not silently bias FF.
    // The geometric value is only an approximation (uniform-coverage
    // assumption); the user is still expected to recalibrate against
    // male reference samples for clinical reporting.
    if (!_args->scale_user_set) {
        const double geom = _compute_geometric_scale();
        if (geom > 0.0) {
            _args->scaling_factor = geom;
            std::cerr << std::fixed << std::setprecision(2)
                      << "[INFO] --scale not supplied; using geometric default "
                      << geom
                      << (_mappability.is_loaded()
                              ? " (computed from -B BED minus PAR overlap).\n"
                              : " (computed from --build canonical genome lengths).\n")
                      << "[INFO]   Recalibrate against male reference samples "
                         "processed with the same -B / --build for clinical use.\n";
        } else {
            std::cerr << "[WARN] cannot auto-compute --scale (--build none with no -B,"
                      << " or -B BED missing autosomal / chrY non-PAR coverage);"
                      << " using fallback " << _args->scaling_factor << ".\n";
        }
    }
    if (!(_args->scaling_factor > 0.0)) {
        // Defensive: should not happen with the auto-compute fallback above,
        // but guards against a future regression.
        throw std::invalid_argument("[ERROR] --scale resolved to a non-positive value ("
                                    + std::to_string(_args->scaling_factor) + ")");
    }
}

FetalFractionRunner::~FetalFractionRunner() {
    if (_args) { delete _args; _args = nullptr; }
}

// ============================================================
// Built-in PAR coordinates
// ============================================================
//
// All coordinates are 1-based inclusive [start, end], matching
// ngslib::GenomeRegion / BaseVar's `parse_region_string` convention.
// Sources:
//   - GRCh38: UCSC / Ensembl par regions (chrY:10001-2781479, chrY:56887903-57217415)
//             https://www.ncbi.nlm.nih.gov/grc/human
//   - GRCh37: UCSC / Ensembl par regions (chrY:10001-2649520, chrY:59034050-59363566)
// We register both `chrY` and `Y` keys so the index works regardless of
// the input BAM's chromosome-naming convention.
void FetalFractionRunner::_load_builtin_par() {
    if (_args->genome_build == GenomeBuild::NONE) return;

    auto add = [&](uint32_t s1, uint32_t e1, uint32_t s2, uint32_t e2) {
        // Coordinates below are already 1-based inclusive; pass through.
        for (const std::string& y : {std::string("chrY"), std::string("Y")}) {
            _par_regions.add_interval(y, s1, e1);
            _par_regions.add_interval(y, s2, e2);
        }
    };
    if (_args->genome_build == GenomeBuild::GRCh38) {
        add(10001u, 2781479u, 56887903u, 57217415u);
    } else {
        add(10001u, 2649520u, 59034050u, 59363566u);
    }
    _par_regions.finalize();
}

// ============================================================
// Geometric --scale auto-computation
// ============================================================
//
// Mathematical basis:
//
//   Y/Auto for a 100% MALE sample (uniform coverage) =
//       L_chrY_uniq / (2 * L_autosome)
//
// so the inverse:
//
//   scale = 2 * L_autosome / L_chrY_uniq
//
// gives FF_male = (Y/Auto) * scale.  The user's actual scale will
// deviate from this geometric ideal because of GC bias, repeat
// content, library prep, etc., so this is only an approximation
// to be used as a default starting point.  Empirical recalibration
// against male reference samples is REQUIRED for clinical use.
//
// Two computation paths:
//
//   (1) -B supplied:  L_autosome = sum of BED autosome bp
//                     L_chrY_uniq = sum of BED chrY bp - chrY ∩ PAR
//   (2) no -B:        L_autosome and L_chrY_uniq are taken from the
//                     canonical GRC assembly report.  Source values
//                     are the chr1..22 cumulative length and the
//                     full chrY length minus PAR1+PAR2 length, both
//                     for the chosen --build.
//
// Returns 0.0 if geometry cannot be inferred (e.g. --build none with
// no -B, or a -B BED missing autosomes or chrY-non-PAR coverage),
// in which case the caller falls back to the constructor's 100.0.
double FetalFractionRunner::_compute_geometric_scale() const {
    if (_mappability.is_loaded()) {
        const uint64_t L_auto = _mappability.total_length(
            [](const std::string& c) { return is_autosome(c); });
        const uint64_t L_y_total = _mappability.total_length(
            [](const std::string& c) { return is_y_chrom(c); });
        const uint64_t L_y_par = _par_regions.is_loaded()
            ? intersection_length_chrY(_mappability, _par_regions) : 0ULL;
        if (L_auto == 0 || L_y_total <= L_y_par) return 0.0;
        const uint64_t L_y_uniq = L_y_total - L_y_par;
        return 2.0 * static_cast<double>(L_auto) / static_cast<double>(L_y_uniq);
    }

    // No -B: use canonical genome geometry per --build.
    // Sources (NCBI / GRC official assembly reports):
    //   GRCh38: chr1..22 = 2,875,001,522 bp; chrY = 57,227,415 bp;
    //           PAR1 (chrY:10001-2781479)    = 2,771,479 bp;
    //           PAR2 (chrY:56887903-57217415) =   329,513 bp;
    //           chrY non-PAR = 57,227,415 - 2,771,479 - 329,513 = 54,126,423 bp.
    //   GRCh37: chr1..22 = 2,881,033,286 bp; chrY = 59,373,566 bp;
    //           PAR1 (chrY:10001-2649520)    = 2,639,520 bp;
    //           PAR2 (chrY:59034050-59363566) =   329,517 bp;
    //           chrY non-PAR = 59,373,566 - 2,639,520 - 329,517 = 56,404,529 bp.
    if (_args->genome_build == GenomeBuild::GRCh38) {
        return 2.0 * 2875001522.0 / 54126423.0;   // ~ 106.23
    }
    if (_args->genome_build == GenomeBuild::GRCh37) {
        return 2.0 * 2881033286.0 / 56404529.0;   // ~ 102.16
    }
    return 0.0;  // --build none and no -B: cannot infer geometry.
}

// ============================================================
// Sample-id resolution
// ============================================================

std::string FetalFractionRunner::_resolve_sample_id(const std::string& fn) const {
    if (_args->filename_has_samplename) {
        return filename_stem(fn);
    }
    // Prefer @RG SM tag; fall back to filename stem.
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

// ============================================================
// Read filter (mirrors basevar motif's policy)
// ============================================================

bool FetalFractionRunner::_passes_filters(const ngslib::BamRecord& r) const {
    if (!r.is_mapped())       return false;
    if (r.is_secondary())     return false;
    if (r.is_supplementary()) return false;
    if (r.is_duplicate())     return false;
    if (r.is_qc_fail())       return false;
    if (r.mapq() < _args->min_mapq) return false;

    if (r.is_paired() && _args->proper_pair && !r.is_proper_pair()) return false;

    if (r.is_paired() && _args->max_insert_size > 0
        && std::abs(r.insert_size()) > _args->max_insert_size) return false;

    return true;
}

// ============================================================
// Per-file worker
// ============================================================

void FetalFractionRunner::_process_one_file(FetalFractionResult& s) const {
    ngslib::Bam bam(s.input_path, "r", _args->reference);
    ngslib::BamHeader& hdr = bam.header();

    std::vector<std::string> regions;
    if (!_args->regions.empty()) {
        ngslib::split(_args->regions, regions, ",");
    }

    auto consume_record = [&](const ngslib::BamRecord& r) {
        ++s.total_reads;
        if (!_passes_filters(r)) { ++s.filtered_reads; return; }

        const std::string chrom = r.tid_name(hdr);
        if (chrom.empty()) { ++s.filtered_reads; return; }

        // 5'-mapped position on the reference (0-based).  Convert to
        // 1-based to query our BedIntervalIndex (BaseVar convention).
        const hts_pos_t pos = r.map_ref_start_pos();
        if (pos < 0) { ++s.filtered_reads; return; }
        const uint32_t pos_1b = static_cast<uint32_t>(pos) + 1u;

        // Optional inclusion mask (mappability BED).
        if (_mappability.is_loaded() && !_mappability.contains(chrom, pos_1b)) {
            ++s.filtered_reads;
            return;
        }

        if (is_y_chrom(chrom)) {
            // Exclude PAR1/PAR2: shared with chrX, would inflate Y for females.
            if (_par_regions.is_loaded() && _par_regions.contains(chrom, pos_1b)) {
                ++s.y_par_excluded;
                return;
            }
            ++s.valid_y;
        } else if (is_autosome(chrom)) {
            ++s.valid_autosomal;
        } else {
            // chrX / chrM / ALT / decoy / HLA / random: pass QC but not
            // informative for FF.  Count as filtered for bookkeeping so
            // total_reads == filtered_reads + valid_y + valid_autosomal
            // + y_par_excluded always holds.
            ++s.filtered_reads;
        }
    };

    ngslib::BamRecord r;
    if (regions.empty()) {
        while (bam.next(r) >= 0) consume_record(r);
    } else {
        for (const auto& reg : regions) {
            if (reg.empty()) continue;
            bam.fetch(reg);
            while (bam.next(r) >= 0) consume_record(r);
        }
    }

    _compute_fetal_fraction(s);
}

// ============================================================
// Sex call + FF computation
// ============================================================

void FetalFractionRunner::_compute_fetal_fraction(FetalFractionResult& s) const {
    if (s.valid_autosomal == 0) {
        s.y_ratio        = 0.0;
        s.sex            = FetalSex::UNDETERMINED;
        s.fetal_fraction = -1.0;
        return;
    }

    s.y_ratio = static_cast<double>(s.valid_y)
              / static_cast<double>(s.valid_autosomal);

    if (s.y_ratio > _args->male_threshold) {
        // MALE fetus: apply noise floor and scale.  Clamp to >= 0 so
        // sub-noise samples don't produce negative FF on the report.
        const double corrected = std::max(0.0, s.y_ratio - _args->y_noise);
        s.sex            = FetalSex::MALE;
        s.fetal_fraction = corrected * _args->scaling_factor;
    } else {
        // FEMALE fetus: under lpWGS we deliberately REFUSE to compute FF.
        s.sex            = FetalSex::FEMALE;
        s.fetal_fraction = -1.0;
    }
}

// ============================================================
// Output writers
// ============================================================

void FetalFractionRunner::_write_tsv(const std::string& path) const {
    std::ofstream ofs(path);
    if (!ofs) {
        throw std::runtime_error("[ERROR] Cannot open output TSV for writing: " + path);
    }
    ofs << "#sample\tsex\tfetal_fraction\ty_ratio\tvalid_autosomal\tvalid_y"
           "\ty_par_excluded\ttotal_reads\tfiltered_reads"
           "\tscale\tnoise\tmale_threshold\n";

    ofs << std::fixed << std::setprecision(8);
    for (const auto& s : _results) {
        ofs << s.sample_id      << "\t"
            << sex_name(s.sex)  << "\t"
            << s.fetal_fraction << "\t"
            << s.y_ratio        << "\t"
            << s.valid_autosomal<< "\t"
            << s.valid_y        << "\t"
            << s.y_par_excluded << "\t"
            << s.total_reads    << "\t"
            << s.filtered_reads << "\t"
            << _args->scaling_factor << "\t"
            << _args->y_noise        << "\t"
            << _args->male_threshold << "\n";
    }
    ofs.close();
}

void FetalFractionRunner::_print_summary(std::ostream& os) const {
    os << "\n=== NIPT cfDNA Male Fetal-Fraction Report ===\n"
       << "Genome build:        " << genome_build_name(_args->genome_build) << "\n"
       << "PAR exclusion:       "
       << (_par_regions.is_loaded() ? "ENABLED" : "disabled") << "\n"
       << "Mappability mask:    "
       << (_mappability.is_loaded() ? _args->mappability_bed : std::string("(none)")) << "\n"
       << "Min MAPQ:            " << _args->min_mapq << "\n"
       << "Proper-pair filter:  " << (_args->proper_pair ? "on" : "off") << "\n"
       << "Max insert size:     "
       << (_args->max_insert_size > 0 ? std::to_string(_args->max_insert_size) : std::string("(off)")) << "\n"
       << "Scale (Y/Auto -> FF):" << _args->scaling_factor
       << (_args->scale_user_set ? "  (user-supplied)" : "  (auto, geometric)") << "\n"
       << "Noise floor:         " << _args->y_noise        << "\n"
       << "MALE threshold:      " << _args->male_threshold << "\n"
       << "Worker threads:      " << _args->thread_num     << "\n"
       << "Samples processed:   " << _results.size()       << "\n";

    os << "\n-- Per-sample summary --\n"
       << std::left
       << std::setw(24) << "Sample"
       << std::setw(14) << "Sex"
       << std::setw(14) << "FF(%)"
       << std::setw(14) << "Y/Auto"
       << std::setw(14) << "Auto"
       << std::setw(14) << "Y(non-PAR)"
       << std::setw(12) << "Y_PAR" << "\n"
       << "----------------------------------------------------------"
          "--------------------------------------\n";

    for (const auto& s : _results) {
        std::string ff_str;
        if (s.sex == FetalSex::MALE) {
            std::ostringstream tmp;
            tmp << std::fixed << std::setprecision(2) << (s.fetal_fraction * 100.0);
            ff_str = tmp.str();
        } else {
            ff_str = "-";
        }
        std::ostringstream y_str;
        y_str << std::scientific << std::setprecision(3) << s.y_ratio;

        os << std::left
           << std::setw(24) << s.sample_id
           << std::setw(14) << sex_name(s.sex)
           << std::setw(14) << ff_str
           << std::setw(14) << y_str.str()
           << std::setw(14) << s.valid_autosomal
           << std::setw(14) << s.valid_y
           << std::setw(12) << s.y_par_excluded << "\n";
    }

    os << "\n[NOTE] Female FF is intentionally NOT computed under lpWGS;\n"
          "       use a SNP-based assay for accurate female FF measurement.\n"
          "       MALE FF reported here uses a simple linear Y/Auto model;\n"
          "       always recalibrate `--scale` against your own reference cohort.\n";
}

// ============================================================
// Calibration mode
// ============================================================

int FetalFractionRunner::_run_calibrate() {
    // Validate the output path up-front so a bad -o does not waste a long run.
    if (!_args->output_file.empty()) {
        std::ofstream probe(_args->output_file);
        if (!probe) {
            throw std::runtime_error("[ERROR] Cannot open output TSV for writing: "
                                     + _args->output_file);
        }
        probe.close();
    }

    // Process all samples first (reuse the standard pipeline).
    const size_t n_files = _args->input_bf.size();
    int n_threads = std::min<int>(_args->thread_num, static_cast<int>(n_files));
    if (n_threads < 1) n_threads = 1;

    std::cerr << "[INFO] basevar fetalfrac v" << BASEVAR_VERSION
              << " - CALIBRATION MODE: treating " << n_files
              << " sample(s) as known pure male (FF=1.0).\n";

    // Stage 1: resolve sample IDs.
    _results.resize(n_files);
    for (size_t i = 0; i < n_files; ++i) {
        _results[i].input_path = _args->input_bf[i];
        _results[i].sample_id  = _resolve_sample_id(_args->input_bf[i]);
    }

    // Stage 2: process files concurrently.
    if (n_threads <= 1 || n_files <= 1) {
        for (size_t i = 0; i < n_files; ++i) {
            std::cerr << "[INFO] [" << (i + 1) << "/" << n_files << "] "
                      << _results[i].sample_id << "  <-  "
                      << _results[i].input_path << "\n";
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

    // Stage 3: compute empirical scale from samples passing male_threshold.
    std::vector<double> valid_ratios;
    size_t n_male = 0, n_female = 0, n_undetermined = 0;
    for (const auto& s : _results) {
        if (s.sex == FetalSex::MALE) {
            ++n_male;
            valid_ratios.push_back(s.y_ratio);
        } else if (s.sex == FetalSex::FEMALE) {
            ++n_female;
        } else {
            ++n_undetermined;
        }
    }

    if (valid_ratios.empty()) {
        std::cerr << "[ERROR] No samples passed the male threshold (y_ratio > "
                  << _args->male_threshold << "). Cannot calibrate.\n";
        return 1;
    }

    // Compute mean noise-corrected y_ratio.
    // The production FF formula is: FF = max(0, y_ratio - noise) * scale.
    // For pure male (FF=1.0): scale = 1 / mean(y_ratio - noise).
    std::vector<double> corrected_ratios;
    corrected_ratios.reserve(valid_ratios.size());
    for (double r : valid_ratios) {
        corrected_ratios.push_back(std::max(0.0, r - _args->y_noise));
    }
    double sum = 0.0;
    for (double r : corrected_ratios) sum += r;
    const double mean_corrected = sum / corrected_ratios.size();

    if (mean_corrected <= 0.0) {
        std::cerr << "[ERROR] Mean noise-corrected y_ratio <= 0. "
                  << "--noise may be too high. Cannot calibrate.\n";
        return 1;
    }
    const double empirical_scale = 1.0 / mean_corrected;

    // Also report raw mean for reference.
    double raw_sum = 0.0;
    for (double r : valid_ratios) raw_sum += r;
    const double mean_raw = raw_sum / valid_ratios.size();

    // Compute standard deviation of corrected ratios for QC.
    double var_sum = 0.0;
    for (double r : corrected_ratios) {
        const double diff = r - mean_corrected;
        var_sum += diff * diff;
    }
    const double sd_ratio = std::sqrt(var_sum / corrected_ratios.size());
    const double cv = (mean_corrected > 0.0) ? (sd_ratio / mean_corrected * 100.0) : 0.0;

    // Report calibration results.
    std::cout << "\n=== Fetal Fraction Calibration Report ===\n"
              << "Mode:                CALIBRATE (known pure male samples)\n"
              << "Genome build:        " << genome_build_name(_args->genome_build) << "\n"
              << "PAR exclusion:       "
              << (_par_regions.is_loaded() ? "ENABLED" : "disabled") << "\n"
              << "Mappability mask:    "
              << (_mappability.is_loaded() ? _args->mappability_bed : std::string("(none)")) << "\n"
              << "Min MAPQ:            " << _args->min_mapq << "\n"
              << "Noise floor:         " << _args->y_noise << "\n"
              << "MALE threshold:      " << _args->male_threshold << "\n"
              << "\n"
              << "Samples processed:   " << _results.size() << "\n"
              << "  MALE (used):       " << n_male << "\n"
              << "  FEMALE (excluded): " << n_female << "\n"
              << "  UNDETERMINED:      " << n_undetermined << "\n"
              << "\n"
              << "-- Per-sample y_ratio (MALE samples) --\n";

    std::cout << std::left << std::setw(24) << "Sample"
              << std::setw(16) << "y_ratio"
              << std::setw(16) << "valid_y"
              << std::setw(16) << "valid_auto" << "\n";
    std::cout << std::string(72, '-') << "\n";

    for (const auto& s : _results) {
        if (s.sex != FetalSex::MALE) continue;
        std::ostringstream y_str;
        y_str << std::scientific << std::setprecision(4) << s.y_ratio;
        std::cout << std::left << std::setw(24) << s.sample_id
                  << std::setw(16) << y_str.str()
                  << std::setw(16) << s.valid_y
                  << std::setw(16) << s.valid_autosomal << "\n";
    }

    std::cout << "\n=== Calibration Result ===\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Noise floor:         " << _args->y_noise << "\n";
    std::cout << "Mean raw y_ratio:    " << mean_raw << "\n";
    std::cout << "Mean corrected:      " << mean_corrected
              << "  (= raw - noise)\n";
    std::cout << "SD corrected:        " << sd_ratio << "\n";
    std::cout << "CV:                  " << cv << "%\n";
    std::cout << "\n";
    std::cout << "*** EMPIRICAL SCALE: " << empirical_scale << " ***\n";
    std::cout << "\n";
    std::cout << "Usage: add --scale " << empirical_scale
              << " to your fetalfrac command for accurate FF estimation.\n";
    std::cout << "\n";
    std::cout << "Example:\n";
    std::cout << "  basevar fetalfrac --scale " << empirical_scale
              << " --noise <your_noise> -o ff.tsv -L bam.list\n";

    // Write calibration TSV if output file specified.
    if (!_args->output_file.empty()) {
        std::ofstream ofs(_args->output_file);
        if (!ofs) {
            throw std::runtime_error("[ERROR] Cannot open output TSV for writing: "
                                     + _args->output_file);
        }
        ofs << "#sample\tsex\tfetal_fraction\ty_ratio\tvalid_autosomal\tvalid_y"
               "\ty_par_excluded\ttotal_reads\tfiltered_reads"
               "\tscale\tnoise\tmale_threshold\n";
        ofs << std::fixed << std::setprecision(8);
        for (const auto& s : _results) {
            ofs << s.sample_id      << "\t"
                << sex_name(s.sex)  << "\t"
                << s.fetal_fraction << "\t"
                << s.y_ratio        << "\t"
                << s.valid_autosomal<< "\t"
                << s.valid_y        << "\t"
                << s.y_par_excluded << "\t"
                << s.total_reads    << "\t"
                << s.filtered_reads << "\t"
                << empirical_scale  << "\t"
                << _args->y_noise        << "\t"
                << _args->male_threshold << "\n";
        }
        ofs << "#calibration_mean_raw_y_ratio=" << mean_raw << "\n";
        ofs << "#calibration_mean_corrected_y_ratio=" << mean_corrected << "\n";
        ofs << "#calibration_sd_corrected_y_ratio=" << sd_ratio << "\n";
        ofs << "#calibration_cv_percent=" << cv << "\n";
        ofs << "#calibration_empirical_scale=" << empirical_scale << "\n";
        ofs.close();
        std::cerr << "[INFO] Wrote calibration TSV to: " << _args->output_file << "\n";
    }

    if (n_female > 0 || n_undetermined > 0) {
        std::cerr << "[WARN] " << (n_female + n_undetermined)
                  << " sample(s) did not pass the male threshold and were excluded "
                     "from calibration. Verify these are truly pure male samples.\n";
    }
    if (cv > 20.0) {
        std::cerr << "[WARN] High CV (" << cv << "%) in y_ratio across samples. "
                  << "Check sample quality / pipeline consistency.\n";
    }

    return 0;
}

// ============================================================
// run()
// ============================================================

int FetalFractionRunner::run() {
    // Calibration mode: different output path.
    if (_args->calibrate) {
        return _run_calibrate();
    }

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

    std::cerr << "[INFO] basevar fetalfrac v" << BASEVAR_VERSION
              << " - estimating male fetal fraction from "
              << n_files << " BAM/CRAM file(s) using "
              << n_threads << " worker thread(s).\n";

    // Stage 1: resolve sample IDs sequentially (cheap, gives ordered output).
    _results.resize(n_files);
    for (size_t i = 0; i < n_files; ++i) {
        _results[i].input_path = _args->input_bf[i];
        _results[i].sample_id  = _resolve_sample_id(_args->input_bf[i]);
    }

    // Stage 2: process files concurrently.  Each task touches only its own
    // FetalFractionResult, so no locking is required.
    if (n_threads <= 1 || n_files <= 1) {
        for (size_t i = 0; i < n_files; ++i) {
            std::cerr << "[INFO] [" << (i + 1) << "/" << n_files << "] "
                      << _results[i].sample_id << "  <-  "
                      << _results[i].input_path << "\n";
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

int fetal_fraction_runner(int argc, char* argv[]) {
    try {
        FetalFractionRunner runner(argc, argv);
        return runner.run();
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }
}

}  // namespace fetalfrac
}  // namespace basevar
