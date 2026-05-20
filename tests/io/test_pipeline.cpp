/**
 * @file test_pipeline.cpp
 * @brief Unit tests for the `basevar pipeline` subcommand.
 *
 *  Coverage:
 *    - load_reference_fai: full load, --chrom filter, file-not-found error.
 *    - parse_region_string: chr / chr:start / chr:start-end / error paths.
 *    - parse_regions:       multi-region, empty token, chrom filter.
 *    - build_command:       expected stdout layout matches Python script.
 *
 *  Build (from tests/io/):
 *    g++ -std=c++17 -O2 -fPIC test_pipeline.cpp \
 *        ../../src/pipeline.cpp ../../src/io/utils.cpp \
 *        -I ../../src -I ../../htslib -o test_pipeline
 *    ./test_pipeline
 */
#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <map>

#include "pipeline.h"
#include "io/utils.h"

using basevar::pipeline::load_reference_fai;
using basevar::pipeline::parse_region_string;
using basevar::pipeline::parse_regions;
using basevar::pipeline::build_command;

// ---- minimal assertion helpers --------------------------------------------
static int g_pass = 0;
static int g_fail = 0;

#define CHECK(cond) do {                                                     \
    if (cond) { ++g_pass; }                                                  \
    else {                                                                   \
        ++g_fail;                                                            \
        std::cerr << "[FAIL] " << __FILE__ << ":" << __LINE__                \
                  << "  " << #cond << "\n";                                  \
    }                                                                        \
} while (0)

#define CHECK_EQ(a, b) do {                                                  \
    auto _aa = (a); auto _bb = (b);                                          \
    if (_aa == _bb) { ++g_pass; }                                            \
    else {                                                                   \
        ++g_fail;                                                            \
        std::cerr << "[FAIL] " << __FILE__ << ":" << __LINE__                \
                  << "  expected '" << _bb << "', got '" << _aa << "'\n";    \
    }                                                                        \
} while (0)

#define CHECK_THROWS(expr) do {                                              \
    bool _threw = false;                                                     \
    try { (void)(expr); } catch (const std::exception&) { _threw = true; }   \
    if (_threw) { ++g_pass; }                                                \
    else {                                                                   \
        ++g_fail;                                                            \
        std::cerr << "[FAIL] " << __FILE__ << ":" << __LINE__                \
                  << "  expected exception from " << #expr << "\n";          \
    }                                                                        \
} while (0)

// ---- helpers --------------------------------------------------------------
// Write a temporary .fai file with the given content; return the path.
static std::string write_tmp_fai(const std::string& content) {
    static int counter = 0;
    std::string path = "/tmp/basevar_test_pipeline_" +
                       std::to_string(::getpid()) + "_" +
                       std::to_string(++counter) + ".fai";
    std::ofstream o(path);
    o << content;
    o.close();
    return path;
}

// ===========================================================================
// tests
// ===========================================================================
static void test_load_reference_fai_full() {
    // 3 chromosomes; only the first two columns are used.
    std::string fai = write_tmp_fai(
        "chr1\t1000\t6\t60\t61\n"
        "chr2\t500\t1100\t60\t61\n"
        "chrX\t250\t1700\t60\t61\n");

    auto entries = load_reference_fai(fai);
    CHECK_EQ(entries.size(), static_cast<size_t>(3));
    CHECK_EQ(entries[0].chrom, std::string("chr1"));
    CHECK_EQ(entries[0].start, 1u);
    CHECK_EQ(entries[0].end,   1000u);
    CHECK_EQ(entries[1].chrom, std::string("chr2"));
    CHECK_EQ(entries[1].end,   500u);
    CHECK_EQ(entries[2].chrom, std::string("chrX"));
    CHECK_EQ(entries[2].end,   250u);

    std::remove(fai.c_str());
}

static void test_load_reference_fai_filter() {
    std::string fai = write_tmp_fai(
        "chr1\t1000\n"
        "chr2\t500\n"
        "chrX\t250\n");

    auto entries = load_reference_fai(fai, {"chr2", "chrX"});
    CHECK_EQ(entries.size(), static_cast<size_t>(2));
    CHECK_EQ(entries[0].chrom, std::string("chr2"));
    CHECK_EQ(entries[1].chrom, std::string("chrX"));

    // Filter matching none of the chromosomes should yield an empty list.
    auto none = load_reference_fai(fai, {"chrY"});
    CHECK(none.empty());

    std::remove(fai.c_str());
}

static void test_load_reference_fai_missing() {
    CHECK_THROWS(load_reference_fai("/tmp/__does_not_exist__.fai"));
}

static void test_parse_region_string_variants() {
    std::map<std::string, uint32_t> m{{"chr1", 1000}, {"chr2", 500}};

    // 'chr' -> whole chromosome
    auto r1 = parse_region_string("chr1", m);
    CHECK_EQ(r1.chrom, std::string("chr1"));
    CHECK_EQ(r1.start, 1u);
    CHECK_EQ(r1.end,   1000u);

    // 'chr:start' -> [start, chrom_len]
    auto r2 = parse_region_string("chr1:200", m);
    CHECK_EQ(r2.start, 200u);
    CHECK_EQ(r2.end,   1000u);

    // 'chr:start-end' -> [start, end]
    auto r3 = parse_region_string("chr2:10-400", m);
    CHECK_EQ(r3.chrom, std::string("chr2"));
    CHECK_EQ(r3.start, 10u);
    CHECK_EQ(r3.end,   400u);

    // whitespace around the token should be tolerated.
    auto r4 = parse_region_string("  chr1:5-50  ", m);
    CHECK_EQ(r4.start, 5u);
    CHECK_EQ(r4.end,   50u);
}

static void test_parse_region_string_errors() {
    std::map<std::string, uint32_t> m{{"chr1", 1000}};

    // empty string
    CHECK_THROWS(parse_region_string("", m));
    // unknown chromosome
    CHECK_THROWS(parse_region_string("chrUNK", m));
    // start < 1
    CHECK_THROWS(parse_region_string("chr1:0-10", m));
    // start > end
    CHECK_THROWS(parse_region_string("chr1:500-100", m));
    // end > chrom length
    CHECK_THROWS(parse_region_string("chr1:10-2000", m));
    // non-numeric start
    CHECK_THROWS(parse_region_string("chr1:abc-10", m));
    // non-numeric end
    CHECK_THROWS(parse_region_string("chr1:10-abc", m));
}

static void test_parse_regions_multi() {
    std::map<std::string, uint32_t> m{
        {"chr1", 1000}, {"chr2", 500}, {"chr3", 200}};

    auto rs = parse_regions("chr1:100-200,chr2,chr3:50", m);
    CHECK_EQ(rs.size(), static_cast<size_t>(3));
    CHECK_EQ(rs[0].chrom, std::string("chr1"));
    CHECK_EQ(rs[0].end,   200u);
    CHECK_EQ(rs[1].chrom, std::string("chr2"));
    CHECK_EQ(rs[1].end,   500u);
    CHECK_EQ(rs[2].chrom, std::string("chr3"));
    CHECK_EQ(rs[2].start, 50u);
    CHECK_EQ(rs[2].end,   200u);

    // empty tokens (trailing comma, double comma) should be silently skipped.
    auto rs2 = parse_regions("chr1,,chr2,", m);
    CHECK_EQ(rs2.size(), static_cast<size_t>(2));
}

static void test_parse_regions_with_chrom_filter() {
    std::map<std::string, uint32_t> m{
        {"chr1", 1000}, {"chr2", 500}, {"chr3", 200}};

    // --chrom restricts AFTER parsing, so unmatched entries are dropped.
    auto rs = parse_regions("chr1,chr2:50-100,chr3", m, {"chr2"});
    CHECK_EQ(rs.size(), static_cast<size_t>(1));
    CHECK_EQ(rs[0].chrom, std::string("chr2"));
    CHECK_EQ(rs[0].start, 50u);
    CHECK_EQ(rs[0].end,   100u);
}

static void test_build_command_layout() {
    std::string cmd = build_command(
        "./bin/basevar caller",
        "-f ref.fa -L bam.list -Q 20 -B 100",
        "chr1", 1, 2000000, "/tmp/out");

    // Must contain the canonical pieces in the right order.
    CHECK(cmd.rfind("time ./bin/basevar caller", 0) == 0);
    CHECK(cmd.find(" -f ref.fa -L bam.list -Q 20 -B 100 ") != std::string::npos);
    CHECK(cmd.find(" -r chr1:1-2000000 ") != std::string::npos);
    CHECK(cmd.find("/tmp/out/chr1_1_2000000.vcf.gz") != std::string::npos);
    CHECK(cmd.find("> /tmp/out/chr1_1_2000000.log") != std::string::npos);
    CHECK(cmd.find("&& echo \"** chr1_1_2000000 done **\"") != std::string::npos);
}

static void test_build_command_no_passthrough() {
    // Empty passthrough must still produce a well-formed command (no extra spaces).
    std::string cmd = build_command("basevar caller", "",
                                    "chrX", 100, 200, "out");
    CHECK_EQ(cmd.substr(0, std::string("time basevar caller -r chrX:100-200").size()),
             std::string("time basevar caller -r chrX:100-200"));
}

// ===========================================================================
// main
// ===========================================================================
int main() {
    std::cout << "== test_pipeline ==" << std::endl;

    test_load_reference_fai_full();
    test_load_reference_fai_filter();
    test_load_reference_fai_missing();

    test_parse_region_string_variants();
    test_parse_region_string_errors();

    test_parse_regions_multi();
    test_parse_regions_with_chrom_filter();

    test_build_command_layout();
    test_build_command_no_passthrough();

    std::cout << "passed: " << g_pass << "  failed: " << g_fail << std::endl;
    return g_fail == 0 ? 0 : 1;
}
