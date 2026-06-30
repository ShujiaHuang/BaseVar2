/**
 * @file test_batchfile_binary.cpp
 * @brief Comprehensive unit tests for binary batchfile format (.bbf/.bbi)
 *
 * Tests cover:
 *   1. Header write/read roundtrip
 *   2. Single-position SNP record roundtrip
 *   3. Insertion record (read_base = "+ACGT", ref_base = "A")
 *   4. Deletion record (read_base = "-ACGT", ref_base = "ACGT")
 *   5. Empty position (depth=0 for all samples)
 *   6. Long-read rpr (>255, >65535)
 *   7. Long indel (>255bp)
 *   8. Strand '*' encoding
 *   9. Multi-position sequential read
 *  10. Index write/load/binary-search
 *  11. Multi-sample with varying depths
 *  12. Mixed variant types at same position
 */
#include <iostream>
#include <cassert>
#include <cstring>
#include <string>
#include <vector>
#include <cstdio>

#include "io/batchfile_binary.h"
#include "io/iobgzf.h"
#include "caller_utils.h"

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
    if (!(cond)) { \
        std::cerr << "[FAIL] " << msg << " (line " << __LINE__ << ")\n"; \
        ++tests_failed; \
    } else { \
        ++tests_passed; \
    } \
} while(0)

// Helper: create AlignInfo with one read
static AlignInfo make_align_info(const std::string &chrom, uint32_t pos,
                                  const std::string &ref_base, const std::string &read_base,
                                  char qual, int rpr, int mapq, char strand) {
    AlignInfo ai(chrom, pos);
    AlignBase ab;
    ab.ref_base  = ref_base;
    ab.read_base = read_base;
    ab.base_qual = qual;
    ab.rpr       = rpr;
    ab.mapq      = mapq;
    ab.map_strand = strand;
    ai.align_bases.push_back(ab);
    return ai;
}

// Helper: create empty AlignInfo (no reads)
static AlignInfo make_empty_align_info(const std::string &chrom, uint32_t pos) {
    return AlignInfo(chrom, pos);
}

// =====================================================================
//  Test 1: Header write/read roundtrip
// =====================================================================
void test_header_roundtrip() {
    std::cout << "--- Test 1: Header roundtrip ---\n";
    const std::string fn = "/tmp/test_bb_header.bbf";

    // Write
    {
        ngslib::BGZFile bf(fn.c_str(), "w");
        std::vector<std::string> sids = {"sample_A", "sample_B", "sample_C"};
        write_binary_header(bf, sids);
        bf.close();
    }

    // Read
    {
        auto sids = read_binary_sample_ids(fn);
        TEST_ASSERT(sids.size() == 3, "header: sample count");
        TEST_ASSERT(sids[0] == "sample_A", "header: sample 0");
        TEST_ASSERT(sids[1] == "sample_B", "header: sample 1");
        TEST_ASSERT(sids[2] == "sample_C", "header: sample 2");
    }

    std::remove(fn.c_str());
}

// =====================================================================
//  Test 2: SNP record roundtrip (single sample, single read)
// =====================================================================
void test_snp_roundtrip() {
    std::cout << "--- Test 2: SNP roundtrip ---\n";
    const std::string fn = "/tmp/test_bb_snp.bbf";
    const std::string chrom = "chr1";
    const uint32_t pos = 1000;

    // Build PosMapVector: 1 sample, 1 read (SNP: ref=A, read=T)
    PosMapVector pmv(1);
    pmv[0][pos] = make_align_info(chrom, pos, "A", "T", 'I', 42, 60, '+');

    ngslib::GenomeRegion gr(chrom, pos, pos);

    // Write
    std::vector<BinaryIndexEntry> idx;
    {
        ngslib::BGZFile bf(fn.c_str(), "w");
        write_binary_header(bf, {"smp1"});
        write_binary_record(bf, pmv, gr, pos, idx);
        bf.close();
    }
    TEST_ASSERT(idx.size() == 1, "SNP: index has 1 entry");
    TEST_ASSERT(idx[0].pos == pos, "SNP: index pos matches");

    // Read
    {
        ngslib::BGZFile bf(fn.c_str(), "r");
        // Skip header
        uint32_t magic; uint16_t ver, sc; uint32_t ids_len;
        bf.read_raw(&magic, 4); bf.read_raw(&ver, 2); bf.read_raw(&sc, 2);
        bf.read_raw(&ids_len, 4);
        if (ids_len > 0) { std::string dummy(ids_len, '\0'); bf.read_raw(&dummy[0], ids_len); }

        // Seek to first record
        bf.seek_virtual(idx[0].virtual_offset);

        std::string ref_id;
        uint32_t ref_pos;
        int total_depth;
        std::vector<BaseType::BatchInfo> samples;

        bool ok = read_binary_record(bf, 1, ref_id, ref_pos, total_depth, samples);
        TEST_ASSERT(ok, "SNP: read success");
        TEST_ASSERT(ref_id == chrom, "SNP: chrom matches");
        TEST_ASSERT(ref_pos == pos, "SNP: pos matches");
        TEST_ASSERT(total_depth == 1, "SNP: depth=1");
        TEST_ASSERT(samples.size() == 1, "SNP: 1 sample");
        TEST_ASSERT(samples[0].ref_bases.size() == 1, "SNP: 1 ref_base");
        TEST_ASSERT(samples[0].ref_bases[0] == "A", "SNP: ref_base=A");
        TEST_ASSERT(samples[0].align_bases[0] == "T", "SNP: align_base=T");
        TEST_ASSERT(samples[0].align_base_quals[0] == 'I', "SNP: qual=I");
        TEST_ASSERT(samples[0].base_pos_ranks[0] == 42, "SNP: rpr=42");
        TEST_ASSERT(samples[0].mapqs[0] == 60, "SNP: mapq=60");
        TEST_ASSERT(samples[0].map_strands[0] == '+', "SNP: strand=+");
        bf.close();
    }

    std::remove(fn.c_str());
}

// =====================================================================
//  Test 3: Insertion record
// =====================================================================
void test_insertion_roundtrip() {
    std::cout << "--- Test 3: Insertion roundtrip ---\n";
    const std::string fn = "/tmp/test_bb_ins.bbf";
    const std::string chrom = "chr1";
    const uint32_t pos = 2000;

    // INS: ref_base="" (empty for INS), read_base="+ACGT"
    PosMapVector pmv(1);
    pmv[0][pos] = make_align_info(chrom, pos, "", "+ACGT", 'H', 10, 50, '-');

    ngslib::GenomeRegion gr(chrom, pos, pos);
    std::vector<BinaryIndexEntry> idx;
    {
        ngslib::BGZFile bf(fn.c_str(), "w");
        write_binary_header(bf, {"smp1"});
        write_binary_record(bf, pmv, gr, pos, idx);
        bf.close();
    }

    {
        ngslib::BGZFile bf(fn.c_str(), "r");
        uint32_t magic; uint16_t ver, sc; uint32_t ids_len;
        bf.read_raw(&magic, 4); bf.read_raw(&ver, 2); bf.read_raw(&sc, 2);
        bf.read_raw(&ids_len, 4);
        if (ids_len > 0) { std::string dummy(ids_len, '\0'); bf.read_raw(&dummy[0], ids_len); }

        bf.seek_virtual(idx[0].virtual_offset);
        std::string ref_id; uint32_t ref_pos; int total_depth;
        std::vector<BaseType::BatchInfo> samples;
        bool ok = read_binary_record(bf, 1, ref_id, ref_pos, total_depth, samples);

        TEST_ASSERT(ok, "INS: read success");
        // ref_base is empty for INS (stored as ref_base_len=0)
        TEST_ASSERT(samples[0].ref_bases[0] == "", "INS: ref_base empty");
        TEST_ASSERT(samples[0].align_bases[0] == "+ACGT", "INS: read_base=+ACGT");
        TEST_ASSERT(samples[0].map_strands[0] == '-', "INS: strand=-");
        bf.close();
    }

    std::remove(fn.c_str());
}

// =====================================================================
//  Test 4: Deletion record (multi-base)
// =====================================================================
void test_deletion_roundtrip() {
    std::cout << "--- Test 4: Deletion roundtrip ---\n";
    const std::string fn = "/tmp/test_bb_del.bbf";
    const std::string chrom = "chr2";
    const uint32_t pos = 5000;

    // DEL: ref_base="ACGT", read_base="-ACGT"
    PosMapVector pmv(1);
    pmv[0][pos] = make_align_info(chrom, pos, "ACGT", "-ACGT", 'G', 100, 40, '+');

    ngslib::GenomeRegion gr(chrom, pos, pos);
    std::vector<BinaryIndexEntry> idx;
    {
        ngslib::BGZFile bf(fn.c_str(), "w");
        write_binary_header(bf, {"smp1"});
        write_binary_record(bf, pmv, gr, pos, idx);
        bf.close();
    }

    {
        ngslib::BGZFile bf(fn.c_str(), "r");
        uint32_t magic; uint16_t ver, sc; uint32_t ids_len;
        bf.read_raw(&magic, 4); bf.read_raw(&ver, 2); bf.read_raw(&sc, 2);
        bf.read_raw(&ids_len, 4);
        if (ids_len > 0) { std::string dummy(ids_len, '\0'); bf.read_raw(&dummy[0], ids_len); }

        bf.seek_virtual(idx[0].virtual_offset);
        std::string ref_id; uint32_t ref_pos; int total_depth;
        std::vector<BaseType::BatchInfo> samples;
        read_binary_record(bf, 1, ref_id, ref_pos, total_depth, samples);

        TEST_ASSERT(samples[0].ref_bases[0] == "ACGT", "DEL: ref_base=ACGT");
        TEST_ASSERT(samples[0].align_bases[0] == "-ACGT", "DEL: read_base=-ACGT");
        TEST_ASSERT(samples[0].base_pos_ranks[0] == 100, "DEL: rpr=100");
        bf.close();
    }

    std::remove(fn.c_str());
}

// =====================================================================
//  Test 5: Empty position (depth=0)
// =====================================================================
void test_empty_position() {
    std::cout << "--- Test 5: Empty position ---\n";
    const std::string fn = "/tmp/test_bb_empty.bbf";
    const std::string chrom = "chr1";
    const uint32_t pos = 3000;

    // No reads at this position
    PosMapVector pmv(1);
    // pmv[0] is empty (no entry for pos)

    ngslib::GenomeRegion gr(chrom, pos, pos);
    std::vector<BinaryIndexEntry> idx;
    {
        ngslib::BGZFile bf(fn.c_str(), "w");
        write_binary_header(bf, {"smp1"});
        write_binary_record(bf, pmv, gr, pos, idx);
        bf.close();
    }

    {
        ngslib::BGZFile bf(fn.c_str(), "r");
        uint32_t magic; uint16_t ver, sc; uint32_t ids_len;
        bf.read_raw(&magic, 4); bf.read_raw(&ver, 2); bf.read_raw(&sc, 2);
        bf.read_raw(&ids_len, 4);
        if (ids_len > 0) { std::string dummy(ids_len, '\0'); bf.read_raw(&dummy[0], ids_len); }

        bf.seek_virtual(idx[0].virtual_offset);
        std::string ref_id; uint32_t ref_pos; int total_depth;
        std::vector<BaseType::BatchInfo> samples;
        read_binary_record(bf, 1, ref_id, ref_pos, total_depth, samples);

        TEST_ASSERT(total_depth == 0, "EMPTY: depth=0");
        TEST_ASSERT(samples.size() == 1, "EMPTY: 1 sample");
        // No-data sample should have placeholder values
        TEST_ASSERT(samples[0].ref_bases[0] == "N", "EMPTY: ref_base=N");
        TEST_ASSERT(samples[0].align_bases[0] == "N", "EMPTY: align_base=N");
        TEST_ASSERT(samples[0].map_strands[0] == '*', "EMPTY: strand=*");
        bf.close();
    }

    std::remove(fn.c_str());
}

// =====================================================================
//  Test 6: Long-read rpr (>255 and >65535)
// =====================================================================
void test_long_read_rpr() {
    std::cout << "--- Test 6: Long-read rpr ---\n";
    const std::string fn = "/tmp/test_bb_lrpr.bbf";
    const std::string chrom = "chr1";
    const uint32_t pos = 4000;

    // rpr = 100000 (simulating a 100kb long read)
    PosMapVector pmv(1);
    pmv[0][pos] = make_align_info(chrom, pos, "A", "T", 'J', 100000, 30, '+');

    ngslib::GenomeRegion gr(chrom, pos, pos);
    std::vector<BinaryIndexEntry> idx;
    {
        ngslib::BGZFile bf(fn.c_str(), "w");
        write_binary_header(bf, {"smp1"});
        write_binary_record(bf, pmv, gr, pos, idx);
        bf.close();
    }

    {
        ngslib::BGZFile bf(fn.c_str(), "r");
        uint32_t magic; uint16_t ver, sc; uint32_t ids_len;
        bf.read_raw(&magic, 4); bf.read_raw(&ver, 2); bf.read_raw(&sc, 2);
        bf.read_raw(&ids_len, 4);
        if (ids_len > 0) { std::string dummy(ids_len, '\0'); bf.read_raw(&dummy[0], ids_len); }

        bf.seek_virtual(idx[0].virtual_offset);
        std::string ref_id; uint32_t ref_pos; int total_depth;
        std::vector<BaseType::BatchInfo> samples;
        read_binary_record(bf, 1, ref_id, ref_pos, total_depth, samples);

        TEST_ASSERT(samples[0].base_pos_ranks[0] == 100000, "LONG_RPR: rpr=100000 preserved");
        bf.close();
    }

    std::remove(fn.c_str());
}

// =====================================================================
//  Test 7: Long indel (>255bp)
// =====================================================================
void test_long_indel() {
    std::cout << "--- Test 7: Long indel (>255bp) ---\n";
    const std::string fn = "/tmp/test_bb_lindel.bbf";
    const std::string chrom = "chr1";
    const uint32_t pos = 6000;

    // 300bp deletion
    std::string del_bases(300, 'A');  // "AAA...A" (300 chars)
    std::string ref_base = del_bases;
    std::string read_base = "-" + del_bases;

    PosMapVector pmv(1);
    pmv[0][pos] = make_align_info(chrom, pos, ref_base, read_base, 'F', 50, 40, '-');

    ngslib::GenomeRegion gr(chrom, pos, pos);
    std::vector<BinaryIndexEntry> idx;
    {
        ngslib::BGZFile bf(fn.c_str(), "w");
        write_binary_header(bf, {"smp1"});
        write_binary_record(bf, pmv, gr, pos, idx);
        bf.close();
    }

    {
        ngslib::BGZFile bf(fn.c_str(), "r");
        uint32_t magic; uint16_t ver, sc; uint32_t ids_len;
        bf.read_raw(&magic, 4); bf.read_raw(&ver, 2); bf.read_raw(&sc, 2);
        bf.read_raw(&ids_len, 4);
        if (ids_len > 0) { std::string dummy(ids_len, '\0'); bf.read_raw(&dummy[0], ids_len); }

        bf.seek_virtual(idx[0].virtual_offset);
        std::string ref_id; uint32_t ref_pos; int total_depth;
        std::vector<BaseType::BatchInfo> samples;
        read_binary_record(bf, 1, ref_id, ref_pos, total_depth, samples);

        TEST_ASSERT(samples[0].ref_bases[0].size() == 300, "LONG_INDEL: ref_base len=300");
        TEST_ASSERT(samples[0].align_bases[0].size() == 301, "LONG_INDEL: read_base len=301 (300 + '-')");
        TEST_ASSERT(samples[0].align_bases[0][0] == '-', "LONG_INDEL: read_base starts with '-'");
        bf.close();
    }

    std::remove(fn.c_str());
}

// =====================================================================
//  Test 8: Strand '*' encoding
// =====================================================================
void test_strand_star() {
    std::cout << "--- Test 8: Strand '*' ---\n";
    const std::string fn = "/tmp/test_bb_star.bbf";
    const std::string chrom = "chr1";
    const uint32_t pos = 7000;

    PosMapVector pmv(1);
    pmv[0][pos] = make_align_info(chrom, pos, "G", "A", 'K', 5, 20, '*');

    ngslib::GenomeRegion gr(chrom, pos, pos);
    std::vector<BinaryIndexEntry> idx;
    {
        ngslib::BGZFile bf(fn.c_str(), "w");
        write_binary_header(bf, {"smp1"});
        write_binary_record(bf, pmv, gr, pos, idx);
        bf.close();
    }

    {
        ngslib::BGZFile bf(fn.c_str(), "r");
        uint32_t magic; uint16_t ver, sc; uint32_t ids_len;
        bf.read_raw(&magic, 4); bf.read_raw(&ver, 2); bf.read_raw(&sc, 2);
        bf.read_raw(&ids_len, 4);
        if (ids_len > 0) { std::string dummy(ids_len, '\0'); bf.read_raw(&dummy[0], ids_len); }

        bf.seek_virtual(idx[0].virtual_offset);
        std::string ref_id; uint32_t ref_pos; int total_depth;
        std::vector<BaseType::BatchInfo> samples;
        read_binary_record(bf, 1, ref_id, ref_pos, total_depth, samples);

        TEST_ASSERT(samples[0].map_strands[0] == '*', "STRAND_STAR: strand='*' preserved");
        bf.close();
    }

    std::remove(fn.c_str());
}

// =====================================================================
//  Test 9: Multi-position sequential read
// =====================================================================
void test_multi_position() {
    std::cout << "--- Test 9: Multi-position sequential read ---\n";
    const std::string fn = "/tmp/test_bb_multi.bbf";
    const std::string chrom = "chr1";

    // Write 5 positions: 100, 200, 300, 400, 500
    std::vector<uint32_t> positions = {100, 200, 300, 400, 500};
    ngslib::GenomeRegion gr(chrom, 100, 500);

    std::vector<BinaryIndexEntry> idx;
    {
        ngslib::BGZFile bf(fn.c_str(), "w");
        write_binary_header(bf, {"smp1"});

        for (uint32_t p : positions) {
            PosMapVector pmv(1);
            pmv[0][p] = make_align_info(chrom, p, "A", "T", 'I', p % 200, 60, '+');
            write_binary_record(bf, pmv, gr, p, idx);
        }
        bf.close();
    }

    TEST_ASSERT(idx.size() == 5, "MULTI_POS: 5 index entries");

    // Read sequentially
    {
        ngslib::BGZFile bf(fn.c_str(), "r");
        uint32_t magic; uint16_t ver, sc; uint32_t ids_len;
        bf.read_raw(&magic, 4); bf.read_raw(&ver, 2); bf.read_raw(&sc, 2);
        bf.read_raw(&ids_len, 4);
        if (ids_len > 0) { std::string dummy(ids_len, '\0'); bf.read_raw(&dummy[0], ids_len); }

        for (size_t i = 0; i < positions.size(); ++i) {
            bf.seek_virtual(idx[i].virtual_offset);
            std::string ref_id; uint32_t ref_pos; int total_depth;
            std::vector<BaseType::BatchInfo> samples;
            bool ok = read_binary_record(bf, 1, ref_id, ref_pos, total_depth, samples);

            TEST_ASSERT(ok, "MULTI_POS: read success at index " + std::to_string(i));
            TEST_ASSERT(ref_pos == positions[i], "MULTI_POS: pos=" + std::to_string(positions[i]));
            TEST_ASSERT(samples[0].base_pos_ranks[0] == (int)(positions[i] % 200),
                        "MULTI_POS: rpr matches at pos " + std::to_string(positions[i]));
        }
        bf.close();
    }

    std::remove(fn.c_str());
}

// =====================================================================
//  Test 10: Index write/load/binary-search
// =====================================================================
void test_index_binary_search() {
    std::cout << "--- Test 10: Index binary search ---\n";
    const std::string idx_fn = "/tmp/test_bb_idx.bbi";

    // Write index with positions: 100, 500, 1000, 5000, 10000
    std::vector<BinaryIndexEntry> entries = {
        {100, 0x1000}, {500, 0x2000}, {1000, 0x3000}, {5000, 0x4000}, {10000, 0x5000}
    };
    write_binary_index(idx_fn, entries);

    // Load and verify
    auto loaded = load_binary_index(idx_fn);
    TEST_ASSERT(loaded.size() == 5, "INDEX: 5 entries loaded");
    TEST_ASSERT(loaded[0].pos == 100, "INDEX: entry 0 pos=100");
    TEST_ASSERT(loaded[4].pos == 10000, "INDEX: entry 4 pos=10000");
    TEST_ASSERT(loaded[2].virtual_offset == 0x3000, "INDEX: entry 2 voffset");

    // Binary search: find first pos >= 400
    auto it = std::lower_bound(loaded.begin(), loaded.end(), (uint32_t)400,
        [](const BinaryIndexEntry &e, uint32_t pos) { return e.pos < pos; });
    TEST_ASSERT(it != loaded.end(), "INDEX_SEARCH: found entry >= 400");
    TEST_ASSERT(it->pos == 500, "INDEX_SEARCH: first pos >= 400 is 500");

    // Binary search: find first pos >= 10000
    it = std::lower_bound(loaded.begin(), loaded.end(), (uint32_t)10000,
        [](const BinaryIndexEntry &e, uint32_t pos) { return e.pos < pos; });
    TEST_ASSERT(it != loaded.end(), "INDEX_SEARCH: found entry >= 10000");
    TEST_ASSERT(it->pos == 10000, "INDEX_SEARCH: exact match at 10000");

    // Binary search: find first pos >= 10001 (beyond range)
    it = std::lower_bound(loaded.begin(), loaded.end(), (uint32_t)10001,
        [](const BinaryIndexEntry &e, uint32_t pos) { return e.pos < pos; });
    TEST_ASSERT(it == loaded.end(), "INDEX_SEARCH: no entry >= 10001");

    // Binary search: find first pos >= 1
    it = std::lower_bound(loaded.begin(), loaded.end(), (uint32_t)1,
        [](const BinaryIndexEntry &e, uint32_t pos) { return e.pos < pos; });
    TEST_ASSERT(it != loaded.end(), "INDEX_SEARCH: found entry >= 1");
    TEST_ASSERT(it->pos == 100, "INDEX_SEARCH: first pos >= 1 is 100");

    std::remove(idx_fn.c_str());
}

// =====================================================================
//  Test 11: Multi-sample with varying depths
// =====================================================================
void test_multi_sample() {
    std::cout << "--- Test 11: Multi-sample varying depths ---\n";
    const std::string fn = "/tmp/test_bb_ms.bbf";
    const std::string chrom = "chr3";
    const uint32_t pos = 8000;

    // 3 samples: smp1 has 2 reads, smp2 has 0 reads, smp3 has 1 read
    PosMapVector pmv(3);

    // smp1: 2 reads at pos
    {
        AlignInfo ai(chrom, pos);
        AlignBase ab1; ab1.ref_base = "C"; ab1.read_base = "T"; ab1.base_qual = 'I';
        ab1.rpr = 10; ab1.mapq = 60; ab1.map_strand = '+';
        ai.align_bases.push_back(ab1);
        AlignBase ab2; ab2.ref_base = "C"; ab2.read_base = "G"; ab2.base_qual = 'H';
        ab2.rpr = 20; ab2.mapq = 55; ab2.map_strand = '-';
        ai.align_bases.push_back(ab2);
        pmv[0][pos] = ai;
    }
    // smp2: no reads (empty)
    // pmv[1] left empty

    // smp3: 1 read
    pmv[2][pos] = make_align_info(chrom, pos, "C", "A", 'J', 30, 40, '+');

    ngslib::GenomeRegion gr(chrom, pos, pos);
    std::vector<BinaryIndexEntry> idx;
    {
        ngslib::BGZFile bf(fn.c_str(), "w");
        write_binary_header(bf, {"smp1", "smp2", "smp3"});
        write_binary_record(bf, pmv, gr, pos, idx);
        bf.close();
    }

    {
        ngslib::BGZFile bf(fn.c_str(), "r");
        uint32_t magic; uint16_t ver, sc; uint32_t ids_len;
        bf.read_raw(&magic, 4); bf.read_raw(&ver, 2); bf.read_raw(&sc, 2);
        bf.read_raw(&ids_len, 4);
        if (ids_len > 0) { std::string dummy(ids_len, '\0'); bf.read_raw(&dummy[0], ids_len); }

        bf.seek_virtual(idx[0].virtual_offset);
        std::string ref_id; uint32_t ref_pos; int total_depth;
        std::vector<BaseType::BatchInfo> samples;
        read_binary_record(bf, 3, ref_id, ref_pos, total_depth, samples);

        TEST_ASSERT(total_depth == 3, "MULTI_SMP: total_depth=3 (2+0+1)");
        TEST_ASSERT(samples.size() == 3, "MULTI_SMP: 3 samples");

        // smp1: 2 reads
        TEST_ASSERT(samples[0].ref_bases.size() == 2, "MULTI_SMP: smp1 has 2 reads");
        TEST_ASSERT(samples[0].align_bases[0] == "T", "MULTI_SMP: smp1 read0=T");
        TEST_ASSERT(samples[0].align_bases[1] == "G", "MULTI_SMP: smp1 read1=G");
        TEST_ASSERT(samples[0].base_pos_ranks[0] == 10, "MULTI_SMP: smp1 rpr0=10");
        TEST_ASSERT(samples[0].base_pos_ranks[1] == 20, "MULTI_SMP: smp1 rpr1=20");
        TEST_ASSERT(samples[0].map_strands[0] == '+', "MULTI_SMP: smp1 strand0=+");
        TEST_ASSERT(samples[0].map_strands[1] == '-', "MULTI_SMP: smp1 strand1=-");

        // smp2: no data
        TEST_ASSERT(samples[1].ref_bases[0] == "N", "MULTI_SMP: smp2 no-data ref=N");
        TEST_ASSERT(samples[1].map_strands[0] == '*', "MULTI_SMP: smp2 no-data strand=*");

        // smp3: 1 read
        TEST_ASSERT(samples[2].align_bases[0] == "A", "MULTI_SMP: smp3 read=A");
        TEST_ASSERT(samples[2].base_pos_ranks[0] == 30, "MULTI_SMP: smp3 rpr=30");
        bf.close();
    }

    std::remove(fn.c_str());
}

// =====================================================================
//  Test 12: Mixed variant types at same position (multi-allelic)
// =====================================================================
void test_mixed_variants() {
    std::cout << "--- Test 12: Mixed variant types ---\n";
    const std::string fn = "/tmp/test_bb_mixed.bbf";
    const std::string chrom = "chr1";
    const uint32_t pos = 9000;

    // 1 sample with 3 reads: SNP, INS, DEL at same position
    PosMapVector pmv(1);
    AlignInfo ai(chrom, pos);

    // Read 1: SNP
    AlignBase ab1; ab1.ref_base = "A"; ab1.read_base = "T"; ab1.base_qual = 'I';
    ab1.rpr = 5; ab1.mapq = 60; ab1.map_strand = '+';
    ai.align_bases.push_back(ab1);

    // Read 2: INS
    AlignBase ab2; ab2.ref_base = ""; ab2.read_base = "+CG"; ab2.base_qual = 'H';
    ab2.rpr = 15; ab2.mapq = 50; ab2.map_strand = '-';
    ai.align_bases.push_back(ab2);

    // Read 3: DEL
    AlignBase ab3; ab3.ref_base = "AC"; ab3.read_base = "-AC"; ab3.base_qual = 'G';
    ab3.rpr = 25; ab3.mapq = 40; ab3.map_strand = '+';
    ai.align_bases.push_back(ab3);

    pmv[0][pos] = ai;

    ngslib::GenomeRegion gr(chrom, pos, pos);
    std::vector<BinaryIndexEntry> idx;
    {
        ngslib::BGZFile bf(fn.c_str(), "w");
        write_binary_header(bf, {"smp1"});
        write_binary_record(bf, pmv, gr, pos, idx);
        bf.close();
    }

    {
        ngslib::BGZFile bf(fn.c_str(), "r");
        uint32_t magic; uint16_t ver, sc; uint32_t ids_len;
        bf.read_raw(&magic, 4); bf.read_raw(&ver, 2); bf.read_raw(&sc, 2);
        bf.read_raw(&ids_len, 4);
        if (ids_len > 0) { std::string dummy(ids_len, '\0'); bf.read_raw(&dummy[0], ids_len); }

        bf.seek_virtual(idx[0].virtual_offset);
        std::string ref_id; uint32_t ref_pos; int total_depth;
        std::vector<BaseType::BatchInfo> samples;
        read_binary_record(bf, 1, ref_id, ref_pos, total_depth, samples);

        TEST_ASSERT(total_depth == 3, "MIXED: depth=3");
        TEST_ASSERT(samples[0].ref_bases.size() == 3, "MIXED: 3 ref_bases");
        TEST_ASSERT(samples[0].align_bases.size() == 3, "MIXED: 3 align_bases");

        // Verify each read preserved correctly
        // v2 format: each read stores its own ref_base (per-read, not position-level)
        TEST_ASSERT(samples[0].ref_bases[0] == "A", "MIXED: ref_base0=A (SNP)");
        TEST_ASSERT(samples[0].ref_bases[1] == "", "MIXED: ref_base1='' (INS)");
        TEST_ASSERT(samples[0].ref_bases[2] == "AC", "MIXED: ref_base2=AC (DEL)");
        TEST_ASSERT(samples[0].align_bases[0] == "T", "MIXED: read0=T (SNP)");
        TEST_ASSERT(samples[0].align_bases[1] == "+CG", "MIXED: read1=+CG (INS)");
        TEST_ASSERT(samples[0].align_bases[2] == "-AC", "MIXED: read2=-AC (DEL)");
        TEST_ASSERT(samples[0].base_pos_ranks[0] == 5, "MIXED: rpr0=5");
        TEST_ASSERT(samples[0].base_pos_ranks[1] == 15, "MIXED: rpr1=15");
        TEST_ASSERT(samples[0].base_pos_ranks[2] == 25, "MIXED: rpr2=25");
        bf.close();
    }

    std::remove(fn.c_str());
}

// =====================================================================
//  Test 13: BBI footer integrity validation
// =====================================================================
void test_bbi_footer_integrity() {
    std::cout << "--- Test 13: BBI footer integrity ---\n";
    const std::string idx_fn = "/tmp/test_bb_footer.bbi";

    // Write a valid index
    std::vector<BinaryIndexEntry> entries = {
        {100, 0x1000}, {500, 0x2000}, {1000, 0x3000}
    };
    write_binary_index(idx_fn, entries);

    // Valid file: validate should return true
    TEST_ASSERT(validate_binary_index(idx_fn) == true, "FOOTER: valid index passes validation");

    // Valid file: load should succeed
    {
        auto loaded = load_binary_index(idx_fn);
        TEST_ASSERT(loaded.size() == 3, "FOOTER: load valid index succeeds");
    }

    // Truncate the file (remove last 4 bytes = footer)
    {
        FILE *fp = fopen(idx_fn.c_str(), "rb+");
        TEST_ASSERT(fp != nullptr, "FOOTER: reopen for truncation");
        fseek(fp, 0, SEEK_END);
        long full_size = ftell(fp);
        // Truncate to full_size - 4 (remove footer)
        // Use freopen trick: write partial content to a temp file
        fseek(fp, 0, SEEK_SET);
        std::vector<char> buf(full_size - 4);
        size_t n = fread(buf.data(), 1, buf.size(), fp);
        fclose(fp);

        // Rewrite truncated file
        fp = fopen(idx_fn.c_str(), "wb");
        fwrite(buf.data(), 1, n, fp);
        fclose(fp);
    }

    // Truncated file: validate should return false
    TEST_ASSERT(validate_binary_index(idx_fn) == false, "FOOTER: truncated index fails validation");

    // Truncated file: load should throw
    {
        bool threw = false;
        try {
            load_binary_index(idx_fn);
        } catch (const std::runtime_error &e) {
            threw = true;
            std::string msg(e.what());
            TEST_ASSERT(msg.find("incomplete") != std::string::npos, "FOOTER: error mentions 'incomplete'");
        }
        TEST_ASSERT(threw, "FOOTER: load throws on truncated index");
    }

    // Non-existent file: validate should return false
    TEST_ASSERT(validate_binary_index("/tmp/nonexistent_bbi_file.bbi") == false,
                "FOOTER: non-existent file returns false");

    std::remove(idx_fn.c_str());
}

// =====================================================================
//  Test 14: Per-read ref_base for Indel normalization bug fix (v2 format)
//
//  Simulates the real-world scenario after _seek_position relocation:
//  At position P-1, a sample has both SNP reads (ref_base="A") and
//  deletion reads (ref_base="AG"). The v1 format stored ref_base at
//  position level, so all reads got the same ref_base, breaking Indel
//  normalization in collect_and_normalized_allele_info.
// =====================================================================
void test_indel_refbase_fix() {
    std::cout << "--- Test 14: Per-read ref_base for Indel normalization fix ---\n";
    const std::string fn = "/tmp/test_bb_indel_fix.bbf";
    const std::string chrom = "chr14";
    const uint32_t pos = 67278837;  // leftmost_pos (anchor base position)

    PosMapVector pmv(2);  // 2 samples

    // Sample 0: 2 SNP reads at this position (ref_base = "A")
    {
        AlignInfo ai(chrom, pos);
        AlignBase ab1; ab1.ref_base = "A"; ab1.read_base = "A"; ab1.base_qual = 'I';
        ab1.rpr = 10; ab1.mapq = 60; ab1.map_strand = '+';
        ai.align_bases.push_back(ab1);
        AlignBase ab2; ab2.ref_base = "A"; ab2.read_base = "G"; ab2.base_qual = 'H';
        ab2.rpr = 20; ab2.mapq = 55; ab2.map_strand = '-';
        ai.align_bases.push_back(ab2);
        pmv[0][pos] = ai;
    }

    // Sample 1: 4 reads mixed - SNP reads (ref="A") + deletion reads (ref="AG")
    // This is exactly what _seek_position produces at the anchor position
    {
        AlignInfo ai(chrom, pos);
        // SNP reads (not relocated, ref_base stays "A")
        AlignBase ab1; ab1.ref_base = "A"; ab1.read_base = "G"; ab1.base_qual = 'I';
        ab1.rpr = 5; ab1.mapq = 60; ab1.map_strand = '+';
        ai.align_bases.push_back(ab1);
        AlignBase ab2; ab2.ref_base = "A"; ab2.read_base = "A"; ab2.base_qual = 'J';
        ab2.rpr = 15; ab2.mapq = 50; ab2.map_strand = '-';
        ai.align_bases.push_back(ab2);
        // Deletion reads (relocated from pos+1, ref_base = "A" + "G" = "AG")
        AlignBase ab3; ab3.ref_base = "AG"; ab3.read_base = "-G"; ab3.base_qual = 'H';
        ab3.rpr = 25; ab3.mapq = 40; ab3.map_strand = '+';
        ai.align_bases.push_back(ab3);
        // Insertion reads (relocated from pos+1, ref_base = "A")
        AlignBase ab4; ab4.ref_base = "A"; ab4.read_base = "+G"; ab4.base_qual = 'G';
        ab4.rpr = 35; ab4.mapq = 45; ab4.map_strand = '-';
        ai.align_bases.push_back(ab4);
        pmv[1][pos] = ai;
    }

    ngslib::GenomeRegion gr(chrom, pos, pos);
    std::vector<BinaryIndexEntry> idx;
    {
        ngslib::BGZFile bf(fn.c_str(), "w");
        write_binary_header(bf, {"smp_snp", "smp_mixed"});
        write_binary_record(bf, pmv, gr, pos, idx);
        bf.close();
    }

    {
        ngslib::BGZFile bf(fn.c_str(), "r");
        uint32_t magic; uint16_t ver, sc, ids_len;
        bf.read_raw(&magic, 4); bf.read_raw(&ver, 2); bf.read_raw(&sc, 2);
        bf.read_raw(&ids_len, 4);
        if (ids_len > 0) { std::string dummy(ids_len, '\0'); bf.read_raw(&dummy[0], ids_len); }

        bf.seek_virtual(idx[0].virtual_offset);
        std::string ref_id; uint32_t ref_pos; int total_depth;
        std::vector<BaseType::BatchInfo> samples;
        read_binary_record(bf, 2, ref_id, ref_pos, total_depth, samples);

        TEST_ASSERT(total_depth == 6, "INDEL_FIX: total_depth=6");

        // Sample 0: all SNP reads, ref_base should be "A" for all
        TEST_ASSERT(samples[0].ref_bases.size() == 2, "INDEL_FIX: smp0 has 2 reads");
        TEST_ASSERT(samples[0].ref_bases[0] == "A", "INDEL_FIX: smp0[0] ref=A");
        TEST_ASSERT(samples[0].ref_bases[1] == "A", "INDEL_FIX: smp0[1] ref=A");

        // Sample 1: mixed reads, each with its own ref_base
        TEST_ASSERT(samples[1].ref_bases.size() == 4, "INDEL_FIX: smp1 has 4 reads");
        TEST_ASSERT(samples[1].ref_bases[0] == "A", "INDEL_FIX: smp1[0] SNP ref=A");
        TEST_ASSERT(samples[1].ref_bases[1] == "A", "INDEL_FIX: smp1[1] SNP ref=A");
        TEST_ASSERT(samples[1].ref_bases[2] == "AG", "INDEL_FIX: smp1[2] DEL ref=AG (key fix!)");
        TEST_ASSERT(samples[1].ref_bases[3] == "A", "INDEL_FIX: smp1[3] INS ref=A");

        // Verify align_bases are preserved correctly
        TEST_ASSERT(samples[1].align_bases[2] == "-G", "INDEL_FIX: smp1[2] read=-G");
        TEST_ASSERT(samples[1].align_bases[3] == "+G", "INDEL_FIX: smp1[3] read=+G");

        bf.close();
    }

    std::remove(fn.c_str());
}

// =====================================================================
//  Main
// =====================================================================
int main() {
    std::cout << "=== Binary Batchfile Unit Tests ===\n\n";

    test_header_roundtrip();
    test_snp_roundtrip();
    test_insertion_roundtrip();
    test_deletion_roundtrip();
    test_empty_position();
    test_long_read_rpr();
    test_long_indel();
    test_strand_star();
    test_multi_position();
    test_index_binary_search();
    test_multi_sample();
    test_mixed_variants();
    test_indel_refbase_fix();
    test_bbi_footer_integrity();

    std::cout << "\n=== Results: " << tests_passed << " passed, "
              << tests_failed << " failed ===\n";

    return tests_failed > 0 ? 1 : 0;
}
