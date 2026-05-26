// Unit test for the `basevar motif` subcommand.
//
// Covers:
//   1. Pure-function tests for reverse_complement_base / reverse_complement.
//   2. extract_5p_motif behaviour for forward- and reverse-mapped reads
//      (using fixtures from tests/data/range.bam).
//   3. End-to-end runner test: invoke MotifCounterRunner with realistic
//      argv, verify the TSV file is well-formed and that motif counts sum
//      to the reported number of used reads.
//
// Run with:
//   make test_motif_counter && ./test_motif_counter
//
// Author: Shujia Huang
// Date:   2026-05-26

#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "io/bam.h"
#include "io/bam_record.h"
#include "motif_counter.h"

using basevar::motif::reverse_complement;
using basevar::motif::reverse_complement_base;
using basevar::motif::extract_5p_motif;
using basevar::motif::extract_5p_motif_from_reference;
using basevar::motif::MotifCounterRunner;

static int g_pass = 0;
static int g_fail = 0;

#define CHECK(cond, msg)                                                     \
    do {                                                                     \
        if (cond) {                                                          \
            ++g_pass;                                                        \
            std::cout << "  [ OK ] " << msg << "\n";                         \
        } else {                                                             \
            ++g_fail;                                                        \
            std::cerr << "  [FAIL] " << msg                                  \
                      << "  (at " << __FILE__ << ":" << __LINE__ << ")\n";   \
        }                                                                    \
    } while (0)

// ------------------------------------------------------------
// Test 1: reverse-complement primitives
// ------------------------------------------------------------
static void test_reverse_complement() {
    std::cout << "[Test 1] reverse_complement primitives\n";
    CHECK(reverse_complement_base('A') == 'T', "A -> T");
    CHECK(reverse_complement_base('C') == 'G', "C -> G");
    CHECK(reverse_complement_base('G') == 'C', "G -> C");
    CHECK(reverse_complement_base('T') == 'A', "T -> A");
    CHECK(reverse_complement_base('N') == 'N', "N -> N");
    CHECK(reverse_complement_base('X') == 'N', "non-ACGT -> N");
    CHECK(reverse_complement_base('a') == 'T', "lowercase a -> T");

    CHECK(reverse_complement("ACGT")  == "ACGT", "ACGT is its own RC");
    CHECK(reverse_complement("AAAA")  == "TTTT", "AAAA -> TTTT");
    CHECK(reverse_complement("CCATG") == "CATGG", "CCATG -> CATGG");
    CHECK(reverse_complement("")      == "",     "empty stays empty");
}

// ------------------------------------------------------------
// Test 2: extract_5p_motif from real BAM records
// ------------------------------------------------------------
static void test_extract_5p_motif_from_bam(const std::string& bam_path) {
    std::cout << "[Test 2] extract_5p_motif on " << bam_path << "\n";
    ngslib::Bam bam(bam_path, "r");
    ngslib::BamRecord r;

    int forward_checked = 0;
    int reverse_checked = 0;
    while (bam.next(r) >= 0) {
        if (!r.is_mapped()) continue;
        if (r.query_length() < 4) continue;

        std::string seq   = r.query_sequence();
        std::string motif = extract_5p_motif(r, 4);
        if (motif.empty()) continue;

        if (r.is_mapped_reverse()) {
            std::string expected = reverse_complement(seq.substr(seq.size() - 4, 4));
            CHECK(motif == expected,
                  "reverse-mapped read: motif equals RC of last 4 bases of SEQ");
            if (++reverse_checked >= 3) break;  // a few are enough
        } else {
            CHECK(motif == seq.substr(0, 4),
                  "forward-mapped read: motif equals first 4 bases of SEQ");
            if (++forward_checked >= 3 && reverse_checked >= 3) break;
        }
    }
    // We don't require both orientations exist, but at least one read must
    // have been examined.
    CHECK((forward_checked + reverse_checked) > 0,
          "examined at least one read for motif extraction");
}

// ------------------------------------------------------------
// Test 3: end-to-end runner
// ------------------------------------------------------------
static long count_lines(const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs) return -1;
    long n = 0;
    std::string line;
    while (std::getline(ifs, line)) ++n;
    return n;
}

static std::string first_line(const std::string& path) {
    std::ifstream ifs(path);
    std::string line;
    if (ifs) std::getline(ifs, line);
    return line;
}

static uint64_t sum_counts_in_tsv(const std::string& path) {
    std::ifstream ifs(path);
    std::string line;
    uint64_t total = 0;
    while (std::getline(ifs, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string sample, motif;
        uint64_t count;
        double freq;
        if (!(ss >> sample >> motif >> count >> freq)) continue;
        total += count;
    }
    return total;
}

// Count distinct sample IDs (column 1) in the TSV.
static size_t count_unique_samples(const std::string& path) {
    std::ifstream ifs(path);
    std::string line;
    std::vector<std::string> seen;
    while (std::getline(ifs, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string sample;
        if (!(ss >> sample)) continue;
        bool dup = false;
        for (const auto& s : seen) { if (s == sample) { dup = true; break; } }
        if (!dup) seen.push_back(sample);
    }
    return seen.size();
}

static void test_runner_end_to_end(const std::string& bam_path,
                                   const std::string& tsv_out) {
    std::cout << "[Test 3] MotifCounterRunner end-to-end on " << bam_path << "\n";

    // Build argv: argv[0] is the subcommand name (matches main.cpp dispatch).
    // -t 1 keeps the test deterministic (single-threaded path).
    std::vector<std::string> arg_storage = {
        "motif", "-l", "4", "-q", "0", "-t", "1", "-o", tsv_out, bam_path
    };
    std::vector<char*> argv;
    for (auto& s : arg_storage) argv.push_back(&s[0]);
    int argc = static_cast<int>(argv.size());

    MotifCounterRunner runner(argc, argv.data());
    int rc = runner.run();
    CHECK(rc == 0, "runner.run() returns 0");

    // Single sample: 4^4 + 1 header = 257 lines.
    long lines = count_lines(tsv_out);
    CHECK(lines == 257, "single-sample TSV has 257 lines (256 motifs + 1 header)");

    CHECK(first_line(tsv_out) == "#sample\tmotif\tcount\tfrequency",
          "TSV header is `#sample motif count frequency`");

    CHECK(count_unique_samples(tsv_out) == 1, "single input -> one sample in TSV");

    CHECK(runner.results().size() == 1, "runner stores exactly one SampleResult");
    CHECK(!runner.results().front().sample_id.empty(),
          "sample_id resolved (from @RG SM tag or filename)");

    uint64_t tsv_sum = sum_counts_in_tsv(tsv_out);
    CHECK(tsv_sum == runner.used_reads(),
          "sum of TSV counts equals runner.used_reads()");
    CHECK(runner.used_reads() > 0, "at least one read was counted");
    CHECK(runner.total_reads() >= runner.used_reads() + runner.filtered_reads()
                                  + runner.n_motifs_with_n(),
          "total_reads >= used + filtered + N-containing motifs");

    // The buggy original code mapped 4-bit code 0 ('illegal') to 'A' which
    // would inflate the AAAA count.  With the correct ngslib::_BASES table
    // the sum of counts across the 4 bases at position 0 must equal
    // used_reads (no double-counting).
    const auto& counts = runner.results().front().motif_counts;
    uint64_t prefix_a = 0, prefix_c = 0, prefix_g = 0, prefix_t = 0;
    for (const auto& kv : counts) {
        if (kv.first.empty()) continue;
        switch (kv.first[0]) {
            case 'A': prefix_a += kv.second; break;
            case 'C': prefix_c += kv.second; break;
            case 'G': prefix_g += kv.second; break;
            case 'T': prefix_t += kv.second; break;
        }
    }
    CHECK(prefix_a + prefix_c + prefix_g + prefix_t == runner.used_reads(),
          "every motif starts with one of A/C/G/T (correct base decoding)");

    // No motif key contains a non-ACGT character (sanity for the keys we emit).
    bool only_acgt = true;
    for (const auto& kv : counts) {
        for (char c : kv.first) {
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T') { only_acgt = false; break; }
        }
        if (!only_acgt) break;
    }
    CHECK(only_acgt, "motif keys contain only A/C/G/T");
}

// ------------------------------------------------------------
// Test 4: multi-sample run -> per-sample blocks + multi-thread
// ------------------------------------------------------------
static void test_runner_multi_sample(const std::string& bam_path,
                                     const std::string& tsv_out) {
    std::cout << "[Test 4] MotifCounterRunner multi-sample (-t 2) on twice the same BAM\n";

    // Two copies of the same BAM (same @RG SM tag) -> two SampleResults.
    // We force --filename-has-samplename so the two entries get distinct
    // sample IDs derived from a copied filename.  ngslib::Bam loads the
    // index eagerly, so we must copy the `.bai` alongside the `.bam`.
    const std::string bam_copy = "test_motif_counter_copy.bam";
    const std::string bai_copy = bam_copy + ".bai";
    {
        std::ifstream src(bam_path, std::ios::binary);
        std::ofstream dst(bam_copy, std::ios::binary);
        dst << src.rdbuf();
    }
    {
        std::ifstream src(bam_path + ".bai", std::ios::binary);
        std::ofstream dst(bai_copy, std::ios::binary);
        dst << src.rdbuf();
    }

    std::vector<std::string> arg_storage = {
        "motif", "-l", "4", "-q", "0", "-t", "2",
        "--filename-has-samplename",
        "-o", tsv_out,
        bam_path, bam_copy
    };
    std::vector<char*> argv;
    for (auto& s : arg_storage) argv.push_back(&s[0]);
    int argc = static_cast<int>(argv.size());

    MotifCounterRunner runner(argc, argv.data());
    int rc = runner.run();
    CHECK(rc == 0, "multi-sample runner.run() returns 0");

    // 2 samples * 256 motifs + 1 header = 513 lines.
    long lines = count_lines(tsv_out);
    CHECK(lines == 513, "two-sample TSV has 513 lines (2 x 256 motifs + 1 header)");

    CHECK(count_unique_samples(tsv_out) == 2,
          "TSV contains exactly two distinct sample IDs");

    CHECK(runner.results().size() == 2, "runner stores exactly two SampleResults");

    // Per-sample counts should match (same BAM, same filters).
    if (runner.results().size() == 2) {
        const auto& a = runner.results()[0];
        const auto& b = runner.results()[1];
        CHECK(a.used_reads == b.used_reads,
              "identical inputs produce identical per-sample used_reads");
        CHECK(a.total_reads == b.total_reads,
              "identical inputs produce identical per-sample total_reads");
        CHECK(a.sample_id != b.sample_id,
              "distinct filenames yield distinct sample IDs");
    }

    // The TSV row sum must equal the sum of per-sample used_reads.
    uint64_t tsv_sum = sum_counts_in_tsv(tsv_out);
    CHECK(tsv_sum == runner.used_reads(),
          "TSV row counts sum to aggregate used_reads (no merging across samples)");

    std::remove(bam_copy.c_str());
    std::remove(bai_copy.c_str());
}

// ------------------------------------------------------------
// Test 5 & 6: --proper-pair and --max-insert-size filter options
// ------------------------------------------------------------
// Data notes for tests/data/range.bam (confirmed via samtools view):
//   - All 112 reads are paired-end (BAM_FPAIRED set).
//   - 3 non-proper-pair reads total; only 1 is R1 (flag=97, isize=2,892,555).
//   - That same R1 read is also the sole read with |isize| > 1000.
//   - With -q 0 --reads R1 baseline: used_reads = 55.
//   - With --proper-pair:           used_reads = 54  (55 - 1).
//   - With --max-insert-size 1000:  used_reads = 54  (55 - 1, same read).
//   - With --max-insert-size 0:     used_reads = 55  (filter disabled).
//   - With --max-insert-size 100:   used_reads < 55  (many reads filtered).
static void test_filter_options(const std::string& bam_path) {
    std::cout << "[Test 5/6] --proper-pair and --max-insert-size filters on " << bam_path << "\n";

    // Helper: run MotifCounterRunner with a custom arg list and return used_reads.
    auto run_and_get_used = [&](std::vector<std::string> args) -> uint64_t {
        std::vector<char*> argv;
        for (auto& s : args) argv.push_back(&s[0]);
        MotifCounterRunner runner(static_cast<int>(argv.size()), argv.data());
        int rc = runner.run();
        CHECK(rc == 0, "runner.run() returns 0 for filter-option test");
        return runner.used_reads();
    };

    // --- Baseline (no extra filters) ---
    uint64_t used_base = run_and_get_used(
        {"motif", "-l", "4", "-q", "0", "-t", "1", bam_path});

    // --- Test 5: --proper-pair ---
    std::cout << "[Test 5] --proper-pair\n";

    uint64_t used_proper = run_and_get_used(
        {"motif", "-l", "4", "-q", "0", "-t", "1", "--proper-pair", bam_path});

    CHECK(used_proper < used_base,
          "--proper-pair reduces used_reads vs baseline");
    // Exactly 1 non-proper R1 read in range.bam (flag=97).
    CHECK(used_base - used_proper == 1,
          "--proper-pair removes exactly 1 non-proper R1 read from range.bam");

    // --- Test 6: --max-insert-size ---
    std::cout << "[Test 6] --max-insert-size\n";

    // Scenario A: 0 = disabled, result must equal baseline.
    uint64_t used_max0 = run_and_get_used(
        {"motif", "-l", "4", "-q", "0", "-t", "1", "--max-insert-size", "0", bam_path});
    CHECK(used_max0 == used_base,
          "--max-insert-size 0 disables the filter (used_reads == baseline)");

    // Scenario B: 1000 bp cap removes the one outlier R1 read (isize=2,892,555).
    uint64_t used_max1000 = run_and_get_used(
        {"motif", "-l", "4", "-q", "0", "-t", "1", "--max-insert-size", "1000", bam_path});
    CHECK(used_max1000 < used_base,
          "--max-insert-size 1000 reduces used_reads vs baseline");
    CHECK(used_base - used_max1000 == 1,
          "--max-insert-size 1000 removes exactly 1 outlier R1 read (isize=2892555)");

    // Scenario C: strict 100 bp cap — many reads are filtered.
    uint64_t used_max100 = run_and_get_used(
        {"motif", "-l", "4", "-q", "0", "-t", "1", "--max-insert-size", "100", bam_path});
    CHECK(used_max100 < used_base,
          "--max-insert-size 100 significantly reduces used_reads");

    // --proper-pair and --max-insert-size 1000 target the same read,
    // so combining them should give the same result as either alone.
    uint64_t used_both = run_and_get_used(
        {"motif", "-l", "4", "-q", "0", "-t", "1",
         "--proper-pair", "--max-insert-size", "1000", bam_path});
    CHECK(used_both == used_max1000,
          "--proper-pair + --max-insert-size 1000 yields same result as either alone "
          "(they target the same read)");
}

// ------------------------------------------------------------
// Test 7: --from-reference end-to-end
// ------------------------------------------------------------
//   Uses tests/data/ce.fa.gz, which is the FASTA matching range.bam's
//   @SQ records (CHROMOSOME_I, CHROMOSOME_II, ...).  The reference-based
//   path should:
//     1. Successfully run end-to-end with the same TSV shape as Test 3.
//     2. Reject when --from-reference is supplied without -f/--reference.
//     3. Produce motif keys made entirely of A/C/G/T (uppercase).
//     4. Have used_reads close to the read-based path (a few reads may
//        differ because of soft-clipping / errors at the fragment
//        terminus).  We only assert it is non-zero and not absurdly
//        far from the read-based count.
static void test_from_reference(const std::string& bam_path,
                                const std::string& fasta_path,
                                const std::string& tsv_out) {
    std::cout << "[Test 7] --from-reference on " << bam_path
              << " with reference " << fasta_path << "\n";

    // 7a. Validation: --from-reference without -f must throw.
    {
        std::vector<std::string> args = {
            "motif", "-l", "4", "-q", "0", "-t", "1",
            "--from-reference", "-o", tsv_out, bam_path
        };
        std::vector<char*> argv;
        for (auto& s : args) argv.push_back(&s[0]);
        bool threw = false;
        try {
            MotifCounterRunner runner(static_cast<int>(argv.size()), argv.data());
        } catch (const std::exception&) {
            threw = true;
        }
        CHECK(threw, "--from-reference without -f/--reference throws");
    }

    // 7b. Read-based baseline (for sanity comparison).
    uint64_t used_read_based = 0;
    {
        std::vector<std::string> args = {
            "motif", "-l", "4", "-q", "0", "-t", "1",
            "-o", tsv_out, bam_path
        };
        std::vector<char*> argv;
        for (auto& s : args) argv.push_back(&s[0]);
        MotifCounterRunner runner(static_cast<int>(argv.size()), argv.data());
        int rc = runner.run();
        CHECK(rc == 0, "baseline (read-based) runner returns 0");
        used_read_based = runner.used_reads();
    }

    // 7c. Reference-based run.
    std::vector<std::string> args = {
        "motif", "-l", "4", "-q", "0", "-t", "1",
        "--from-reference", "-f", fasta_path,
        "-o", tsv_out, bam_path
    };
    std::vector<char*> argv;
    for (auto& s : args) argv.push_back(&s[0]);
    MotifCounterRunner runner(static_cast<int>(argv.size()), argv.data());
    int rc = runner.run();
    CHECK(rc == 0, "--from-reference runner returns 0");

    CHECK(runner.used_reads() > 0, "--from-reference produces non-zero used_reads");

    // The two paths use the same filters and the same set of reads;
    // for an aligned read of length >= k, both extractors return a k-mer.
    // The two counts may differ only if the reference fetch fell off the
    // end of a chromosome -- which should not happen for range.bam.
    CHECK(runner.used_reads() == used_read_based,
          "--from-reference and read-based path agree on used_reads count");

    // TSV shape: 1 sample x 4^4 motifs + header = 257 lines.
    long lines = count_lines(tsv_out);
    CHECK(lines == 257, "--from-reference TSV has 257 lines (256 motifs + 1 header)");

    CHECK(first_line(tsv_out) == "#sample\tmotif\tcount\tfrequency",
          "--from-reference TSV header is correct");

    // Every motif key must be uppercase A/C/G/T (lowercase soft-masked
    // ref bases must be normalised to uppercase before counting).
    bool only_acgt = true;
    for (const auto& kv : runner.results().front().motif_counts) {
        for (char c : kv.first) {
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T') { only_acgt = false; break; }
        }
        if (!only_acgt) break;
    }
    CHECK(only_acgt, "--from-reference motif keys are uppercase A/C/G/T");

    // Sum of TSV counts equals reported used_reads (no double-counting).
    uint64_t tsv_sum = sum_counts_in_tsv(tsv_out);
    CHECK(tsv_sum == runner.used_reads(),
          "--from-reference: TSV row counts sum to used_reads");
}

int main(int argc, char* argv[]) {
    // Allow caller to override the test BAM path; default matches the
    // existing fixture used by test_bam.cpp.
    std::string bam_path = (argc > 1) ? argv[1] : "../data/range.bam";
    std::string tsv_out  = "test_motif_counter.tsv";

    std::cout << "Using BAM fixture: " << bam_path << "\n\n";

    test_reverse_complement();
    std::cout << "\n";

    try {
        test_extract_5p_motif_from_bam(bam_path);
    } catch (const std::exception& e) {
        std::cerr << "[FAIL] exception in extract_5p_motif test: " << e.what() << "\n";
        ++g_fail;
    }
    std::cout << "\n";

    try {
        test_runner_end_to_end(bam_path, tsv_out);
    } catch (const std::exception& e) {
        std::cerr << "[FAIL] exception in runner end-to-end test: " << e.what() << "\n";
        ++g_fail;
    }
    std::cout << "\n";

    try {
        test_runner_multi_sample(bam_path, "test_motif_counter_multi.tsv");
    } catch (const std::exception& e) {
        std::cerr << "[FAIL] exception in multi-sample test: " << e.what() << "\n";
        ++g_fail;
    }
    std::cout << "\n";

    try {
        test_filter_options(bam_path);
    } catch (const std::exception& e) {
        std::cerr << "[FAIL] exception in filter_options test: " << e.what() << "\n";
        ++g_fail;
    }
    std::cout << "\n";

    // Test 7 needs the FASTA companion to range.bam.  ce.fa.gz lives next
    // to range.bam in tests/data, so derive its path from the bam_path.
    {
        std::string fasta_path = "../data/ce.fa.gz";
        // If the user passed a custom bam_path, also let them override the FASTA
        // via argv[2].  Otherwise fall back to the default fixture location.
        if (argc > 2) fasta_path = argv[2];
        try {
            test_from_reference(bam_path, fasta_path, "test_motif_counter_ref.tsv");
        } catch (const std::exception& e) {
            std::cerr << "[FAIL] exception in from_reference test: " << e.what() << "\n";
            ++g_fail;
        }
    }

    // Cleanup
    std::remove(tsv_out.c_str());
    std::remove("test_motif_counter_multi.tsv");
    std::remove("test_motif_counter_ref.tsv");

    std::cout << "\n=================================\n";
    std::cout << "  passed: " << g_pass << "\n";
    std::cout << "  failed: " << g_fail << "\n";
    std::cout << "=================================\n";
    return g_fail == 0 ? 0 : 1;
}
