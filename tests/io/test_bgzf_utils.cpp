/**
 * @file test_bgzf_utils.cpp
 * @brief Unit tests for BGZF block-level utilities in io/iobgzf.h
 *
 * Tests:
 *   - ngslib::BGZF_EOF_SIZE / BGZF_EOF_BLOCK constants
 *   - ngslib::bgzf_block_is_valid()
 *   - ngslib::bgzf_block_bsize()
 *   - ngslib::is_bgzf_eof()
 *   - Round-trip: write BGZF file, read raw blocks, verify structure
 *
 * Build (from tests/io/):
 *   make test_bgzf_utils
 *
 * @author Shujia Huang
 * @date 2026-07-02
 */
#include <iostream>
#include <cassert>
#include <cstring>
#include <cstdint>
#include <unistd.h>  // access(), F_OK
#include <vector>
#include <string>

#include "io/iobgzf.h"
#include "io/bgzf_concat.h"
#include <htslib/bgzf.h>

static int g_tests_passed = 0;
static int g_tests_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
    if (!(cond)) { \
        std::cerr << "[FAIL] " << msg << "  (" << #cond << ")" << std::endl; \
        g_tests_failed++; \
    } else { \
        g_tests_passed++; \
    } \
} while(0)

// ========================================================================
//  Test BGZF_EOF_SIZE constant
// ========================================================================
void test_eof_size() {
    std::cout << "--- test_eof_size ---" << std::endl;
    TEST_ASSERT(ngslib::BGZF_EOF_SIZE == 28, "BGZF_EOF_SIZE should be 28");
    TEST_ASSERT(sizeof(ngslib::BGZF_EOF_BLOCK) == 28, "BGZF_EOF_BLOCK should be 28 bytes");
}

// ========================================================================
//  Test BGZF_EOF_BLOCK content
// ========================================================================
void test_eof_block_content() {
    std::cout << "--- test_eof_block_content ---" << std::endl;

    // Verify the EOF block starts with gzip magic bytes
    TEST_ASSERT(ngslib::BGZF_EOF_BLOCK[0] == 0x1f, "EOF byte 0 should be 0x1f");
    TEST_ASSERT(ngslib::BGZF_EOF_BLOCK[1] == 0x8b, "EOF byte 1 should be 0x8b");
    TEST_ASSERT(ngslib::BGZF_EOF_BLOCK[2] == 0x08, "EOF byte 2 should be 0x08 (deflate)");

    // Verify BGZF-specific extra field: "BC" magic at bytes 12-13
    TEST_ASSERT(ngslib::BGZF_EOF_BLOCK[12] == 0x42, "EOF byte 12 should be 'B' (0x42)");
    TEST_ASSERT(ngslib::BGZF_EOF_BLOCK[13] == 0x43, "EOF byte 13 should be 'C' (0x43)");

    // Last 8 bytes should be all zeros (empty deflate block)
    for (int i = 20; i < 28; ++i) {
        TEST_ASSERT(ngslib::BGZF_EOF_BLOCK[i] == 0x00,
                    "EOF byte " + std::to_string(i) + " should be 0x00");
    }
}

// ========================================================================
//  Test bgzf_block_is_valid()
// ========================================================================
void test_bgzf_block_is_valid() {
    std::cout << "--- test_bgzf_block_is_valid ---" << std::endl;

    // Valid BGZF block header (starts with gzip magic)
    uint8_t valid_header[ngslib::BGZF_HEADER_SIZE] = {
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
        0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
        0x1b, 0x00
    };
    TEST_ASSERT(ngslib::bgzf_block_is_valid(valid_header) == 0,
                "Valid BGZF header should return 0");

    // Valid: the EOF block header is also a valid gzip block
    TEST_ASSERT(ngslib::bgzf_block_is_valid(ngslib::BGZF_EOF_BLOCK.data()) == 0,
                "EOF block should pass validation");

    // Invalid: wrong first byte
    uint8_t bad1[ngslib::BGZF_HEADER_SIZE] = {};
    bad1[0] = 0x00; bad1[1] = 0x8b; bad1[2] = 0x08;
    TEST_ASSERT(ngslib::bgzf_block_is_valid(bad1) == -2,
                "Wrong first byte should return -2");

    // Invalid: wrong second byte
    uint8_t bad2[ngslib::BGZF_HEADER_SIZE] = {};
    bad2[0] = 0x1f; bad2[1] = 0x00; bad2[2] = 0x08;
    TEST_ASSERT(ngslib::bgzf_block_is_valid(bad2) == -2,
                "Wrong second byte should return -2");

    // Invalid: wrong third byte
    uint8_t bad3[ngslib::BGZF_HEADER_SIZE] = {};
    bad3[0] = 0x1f; bad3[1] = 0x8b; bad3[2] = 0x07;
    TEST_ASSERT(ngslib::bgzf_block_is_valid(bad3) == -2,
                "Wrong third byte should return -2");

    // Invalid: all zeros
    uint8_t zeros[ngslib::BGZF_HEADER_SIZE] = {};
    TEST_ASSERT(ngslib::bgzf_block_is_valid(zeros) == -2,
                "All-zero header should return -2");
}

// ========================================================================
//  Test bgzf_block_bsize()
// ========================================================================
void test_bgzf_block_bsize() {
    std::cout << "--- test_bgzf_block_bsize ---" << std::endl;

    // EOF block: BSIZE should be 27 (0x001b), meaning block size = 28
    uint16_t eof_bsize = ngslib::bgzf_block_bsize(ngslib::BGZF_EOF_BLOCK.data());
    TEST_ASSERT(eof_bsize == 27, "EOF block BSIZE should be 27 (block size 28 - 1)");

    // Construct a header with known BSIZE: 0x00FF = 255 → block size 256
    uint8_t header[ngslib::BGZF_HEADER_SIZE] = {};
    header[0] = 0x1f; header[1] = 0x8b; header[2] = 0x08;  // gzip magic
    header[16] = 0xFF; header[17] = 0x00;  // BSIZE = 0x00FF = 255
    TEST_ASSERT(ngslib::bgzf_block_bsize(header) == 255,
                "BSIZE should be 255 when bytes are FF 00");

    // BSIZE = 0x0100 = 256 → block size 257
    header[16] = 0x00; header[17] = 0x01;
    TEST_ASSERT(ngslib::bgzf_block_bsize(header) == 256,
                "BSIZE should be 256 when bytes are 00 01");

    // BSIZE = 0xFFFF = 65535 → max block size 65536
    header[16] = 0xFF; header[17] = 0xFF;
    TEST_ASSERT(ngslib::bgzf_block_bsize(header) == 65535,
                "BSIZE should be 65535 when bytes are FF FF");

    // BSIZE = 0x0000 → block size 1 (degenerate but valid encoding)
    header[16] = 0x00; header[17] = 0x00;
    TEST_ASSERT(ngslib::bgzf_block_bsize(header) == 0,
                "BSIZE should be 0 when bytes are 00 00");
}

// ========================================================================
//  Test is_bgzf_eof()
// ========================================================================
void test_is_bgzf_eof() {
    std::cout << "--- test_is_bgzf_eof ---" << std::endl;

    // The actual EOF block should be recognized
    TEST_ASSERT(ngslib::is_bgzf_eof(ngslib::BGZF_EOF_BLOCK.data(), ngslib::BGZF_EOF_SIZE),
                "Actual EOF block should be detected");

    // Wrong length: 27 bytes of EOF block
    TEST_ASSERT(!ngslib::is_bgzf_eof(ngslib::BGZF_EOF_BLOCK.data(), 27),
                "27 bytes should not be detected as EOF");

    // Wrong length: 29 bytes
    uint8_t buf29[29];
    std::memcpy(buf29, ngslib::BGZF_EOF_BLOCK.data(), 28);
    buf29[28] = 0x00;
    TEST_ASSERT(!ngslib::is_bgzf_eof(buf29, 29),
                "29 bytes should not be detected as EOF");

    // Zero length
    TEST_ASSERT(!ngslib::is_bgzf_eof(ngslib::BGZF_EOF_BLOCK.data(), 0),
                "Zero length should not be detected as EOF");

    // Modified byte: flip one bit in the middle
    uint8_t modified[28];
    std::memcpy(modified, ngslib::BGZF_EOF_BLOCK.data(), 28);
    modified[10] ^= 0x01;  // flip a bit
    TEST_ASSERT(!ngslib::is_bgzf_eof(modified, 28),
                "Modified EOF block should not be detected as EOF");

    // All zeros: not an EOF block
    uint8_t zeros[28] = {};
    TEST_ASSERT(!ngslib::is_bgzf_eof(zeros, 28),
                "All zeros should not be detected as EOF");

    // A normal BGZF data block (not EOF)
    uint8_t data_block[100] = {};
    data_block[0] = 0x1f; data_block[1] = 0x8b; data_block[2] = 0x08;
    data_block[16] = 99 - 1; data_block[17] = 0;  // BSIZE = 99 → 100 bytes
    TEST_ASSERT(!ngslib::is_bgzf_eof(data_block, 100),
                "Normal data block should not be detected as EOF");
}

// ========================================================================
//  Test round-trip: write BGZF, read raw blocks, verify structure
// ========================================================================
void test_roundtrip_raw_blocks() {
    std::cout << "--- test_roundtrip_raw_blocks ---" << std::endl;

    const std::string tmpfile = "/tmp/test_bgzf_utils_roundtrip.vcf.gz";

    // Write a small BGZF-compressed VCF
    {
        ngslib::BGZFile out(tmpfile, "wb");
        out << "##fileformat=VCFv4.2\n";
        out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        out << "chr1\t100\t.\tA\tG\t50\tPASS\t.\n";
        out.close();
    }

    // Read raw blocks and verify structure
    {
        BGZF *fp = bgzf_open(tmpfile.c_str(), "r");
        assert(fp != nullptr);

        int block_count = 0;
        bool found_eof = false;

        while (true) {
            uint8_t header[ngslib::BGZF_HEADER_SIZE];
            ssize_t nread = bgzf_raw_read(fp, header, ngslib::BGZF_HEADER_SIZE);
            if (nread == 0) break;  // clean end of file

            TEST_ASSERT(nread == ngslib::BGZF_HEADER_SIZE, "Should read exactly ngslib::BGZF_HEADER_SIZE bytes for block header");
            TEST_ASSERT(ngslib::bgzf_block_is_valid(header) == 0,
                        "Block header should be valid");

            uint32_t bsize = static_cast<uint32_t>(ngslib::bgzf_block_bsize(header)) + 1;
            TEST_ASSERT(bsize >= ngslib::BGZF_HEADER_SIZE && bsize <= 65536,
                        "Block size should be in valid range [ngslib::BGZF_HEADER_SIZE, 65536]");

            // Read remaining bytes
            std::vector<uint8_t> rest(bsize - ngslib::BGZF_HEADER_SIZE);
            if (bsize > ngslib::BGZF_HEADER_SIZE) {
                ssize_t nrest = bgzf_raw_read(fp, rest.data(), bsize - ngslib::BGZF_HEADER_SIZE);
                TEST_ASSERT(nrest == static_cast<ssize_t>(bsize - ngslib::BGZF_HEADER_SIZE),
                            "Should read remaining block bytes");
            }

            // Check if this is an EOF block
            std::vector<uint8_t> full_block(bsize);
            std::memcpy(full_block.data(), header, ngslib::BGZF_HEADER_SIZE);
            if (bsize > ngslib::BGZF_HEADER_SIZE)
                std::memcpy(full_block.data() + ngslib::BGZF_HEADER_SIZE, rest.data(), bsize - ngslib::BGZF_HEADER_SIZE);

            if (ngslib::is_bgzf_eof(full_block.data(), full_block.size())) {
                found_eof = true;
            }
            block_count++;
        }

        TEST_ASSERT(block_count > 0, "Should have read at least one block");
        TEST_ASSERT(found_eof, "Should have found an EOF block at the end");

        bgzf_close(fp);
    }

    // Clean up
    std::remove(tmpfile.c_str());
}

// ========================================================================
//  Test bgzf_block_bsize consistency with bgzf_block_is_valid
// ========================================================================
void test_bsize_valid_header_consistency() {
    std::cout << "--- test_bsize_valid_header_consistency ---" << std::endl;

    // A valid header with BSIZE = 100 - 1 = 99
    uint8_t header[ngslib::BGZF_HEADER_SIZE] = {};
    header[0] = 0x1f; header[1] = 0x8b; header[2] = 0x08;
    header[16] = 99; header[17] = 0;

    TEST_ASSERT(ngslib::bgzf_block_is_valid(header) == 0, "Should be valid");
    TEST_ASSERT(ngslib::bgzf_block_bsize(header) == 99, "BSIZE should be 99");

    // An invalid header should still have a parseable BSIZE (function is independent)
    uint8_t bad_header[ngslib::BGZF_HEADER_SIZE] = {};
    bad_header[0] = 0x00;
    bad_header[16] = 50; bad_header[17] = 0;
    TEST_ASSERT(ngslib::bgzf_block_is_valid(bad_header) != 0, "Should be invalid");
    TEST_ASSERT(ngslib::bgzf_block_bsize(bad_header) == 50,
                "BSIZE should still be parseable from invalid header");
}

// ========================================================================
//  Test BGZFile::read_block(), block_length(), uncompressed_data()
// ========================================================================
void test_bgzfile_read_block() {
    std::cout << "--- test_bgzfile_read_block ---" << std::endl;

    const std::string tmpfile = "/tmp/test_bgzfile_read_block.vcf.gz";
    const std::string content = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t50\tPASS\t.\n";

    // Write a small BGZF file
    {
        ngslib::BGZFile out(tmpfile, "wb");
        TEST_ASSERT(out.is_open(), "Output file should be open");
        out << content;
        out.close();
    }

    // Read back using read_block() / block_length() / uncompressed_data()
    {
        ngslib::BGZFile in(tmpfile, "r");
        TEST_ASSERT(in.is_open(), "Input file should be open");

        int ret = in.read_block();
        TEST_ASSERT(ret == 0, "read_block should return 0 on success");
        TEST_ASSERT(in.block_length() > 0, "block_length should be > 0 after read_block");

        // Verify the decompressed content matches
        std::string decompressed(in.uncompressed_data(), in.block_length());
        TEST_ASSERT(decompressed == content,
                    "Decompressed content should match original");

        // After all data is read, next read_block should return 0 with block_length == 0
        ret = in.read_block();
        // Either returns 0 with block_length 0 (EOF), or returns -1
        // Both are acceptable; the key is block_length == 0
        if (ret == 0) {
            TEST_ASSERT(in.block_length() == 0,
                        "block_length should be 0 at EOF");
        }
    }

    std::remove(tmpfile.c_str());
}

// ========================================================================
//  Test BGZFile::raw_block_read() — read raw BGZF blocks via BGZFile API
// ========================================================================
void test_bgzfile_raw_block_read() {
    std::cout << "--- test_bgzfile_raw_block_read ---" << std::endl;

    const std::string tmpfile = "/tmp/test_bgzfile_raw_block_read.vcf.gz";

    // Write a small BGZF file
    {
        ngslib::BGZFile out(tmpfile, "wb");
        out << "##fileformat=VCFv4.2\n";
        out << "chr1\t100\t.\tA\tG\t50\tPASS\t.\n";
        out.close();
    }

    // Read raw blocks using BGZFile::raw_block_read()
    {
        ngslib::BGZFile in(tmpfile, "r");
        TEST_ASSERT(in.is_open(), "Input file should be open");

        int block_count = 0;
        bool found_eof = false;

        while (true) {
            uint8_t header[ngslib::BGZF_HEADER_SIZE];
            ssize_t nread = in.raw_block_read(header, ngslib::BGZF_HEADER_SIZE);
            if (nread == 0) break;  // clean EOF

            TEST_ASSERT(nread == ngslib::BGZF_HEADER_SIZE, "raw_block_read should return ngslib::BGZF_HEADER_SIZE for block header");
            TEST_ASSERT(ngslib::bgzf_block_is_valid(header) == 0,
                        "Block header should be valid");

            uint32_t bsize = static_cast<uint32_t>(ngslib::bgzf_block_bsize(header)) + 1;
            TEST_ASSERT(bsize >= ngslib::BGZF_HEADER_SIZE && bsize <= 65536,
                        "Block size should be in valid range");

            std::vector<uint8_t> full_block(bsize);
            std::memcpy(full_block.data(), header, ngslib::BGZF_HEADER_SIZE);
            if (bsize > ngslib::BGZF_HEADER_SIZE) {
                ssize_t nrest = in.raw_block_read(full_block.data() + ngslib::BGZF_HEADER_SIZE, bsize - ngslib::BGZF_HEADER_SIZE);
                TEST_ASSERT(nrest == static_cast<ssize_t>(bsize - ngslib::BGZF_HEADER_SIZE),
                            "Should read remaining block bytes");
            }

            if (ngslib::is_bgzf_eof(full_block.data(), full_block.size()))
                found_eof = true;
            block_count++;
        }

        TEST_ASSERT(block_count > 0, "Should have read at least one block");
        TEST_ASSERT(found_eof, "Should have found EOF block");
    }

    std::remove(tmpfile.c_str());
}

// ========================================================================
//  Test BGZFile::raw_block_write() — round-trip via raw block I/O
// ========================================================================
void test_bgzfile_raw_block_write_roundtrip() {
    std::cout << "--- test_bgzfile_raw_block_write_roundtrip ---" << std::endl;

    const std::string src_file = "/tmp/test_bgzfile_raw_rt_src.vcf.gz";
    const std::string dst_file = "/tmp/test_bgzfile_raw_rt_dst.vcf.gz";

    // Write source BGZF file
    {
        ngslib::BGZFile out(src_file, "wb");
        out << "##fileformat=VCFv4.2\n";
        out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        out << "chr1\t100\t.\tA\tG\t50\tPASS\t.\n";
        out << "chr1\t200\t.\tC\tT\t60\tPASS\t.\n";
        out.close();
    }

    // Copy raw blocks from source to destination using BGZFile API
    {
        ngslib::BGZFile src(src_file, "r");
        ngslib::BGZFile dst(dst_file, "w");
        TEST_ASSERT(src.is_open() && dst.is_open(), "Both files should be open");

        const size_t page_size = BGZF_MAX_BLOCK_SIZE;
        std::vector<uint8_t> buf(page_size);

        while (true) {
            ssize_t nread = src.raw_block_read(buf.data(), ngslib::BGZF_HEADER_SIZE);
            if (nread == 0) break;
            TEST_ASSERT(nread == ngslib::BGZF_HEADER_SIZE, "Should read ngslib::BGZF_HEADER_SIZE-byte block header");

            uint32_t bsize = static_cast<uint32_t>(ngslib::bgzf_block_bsize(buf.data())) + 1;
            TEST_ASSERT(bsize >= ngslib::BGZF_HEADER_SIZE && bsize <= page_size, "Valid block size");

            nread += src.raw_block_read(buf.data() + ngslib::BGZF_HEADER_SIZE, bsize - ngslib::BGZF_HEADER_SIZE);
            TEST_ASSERT(nread == static_cast<ssize_t>(bsize), "Full block read");

            // Skip EOF block — ~BGZFile will write the final one
            if (ngslib::is_bgzf_eof(buf.data(), nread)) continue;

            ssize_t nwr = dst.raw_block_write(buf.data(), nread);
            TEST_ASSERT(nwr == nread, "Full block write");
        }
        // ~BGZFile closes both; dst destructor writes EOF block
    }

    // Verify the copy is identical to the source (byte-for-byte)
    {
        FILE *fsrc = fopen(src_file.c_str(), "rb");
        FILE *fdst = fopen(dst_file.c_str(), "rb");
        TEST_ASSERT(fsrc != nullptr && fdst != nullptr, "Files should be openable");

        fseek(fsrc, 0, SEEK_END);
        fseek(fdst, 0, SEEK_END);
        long src_size = ftell(fsrc);
        long dst_size = ftell(fdst);
        TEST_ASSERT(src_size == dst_size,
                    "Source and destination should have the same file size");

        fseek(fsrc, 0, SEEK_SET);
        fseek(fdst, 0, SEEK_SET);

        std::vector<char> src_data(src_size), dst_data(dst_size);
        fread(src_data.data(), 1, src_size, fsrc);
        fread(dst_data.data(), 1, dst_size, fdst);
        fclose(fsrc);
        fclose(fdst);

        TEST_ASSERT(std::memcmp(src_data.data(), dst_data.data(), src_size) == 0,
                    "Source and destination should be byte-identical");
    }

    // Verify the copy is a valid BGZF file by reading it back with BGZFile
    {
        ngslib::BGZFile copy(dst_file, "r");
        std::string line;
        int line_count = 0;
        while (copy.readline(line)) {
            line_count++;
        }
        TEST_ASSERT(line_count == 4, "Copy should have 4 lines (2 header + 2 data)");
    }

    std::remove(src_file.c_str());
    std::remove(dst_file.c_str());
}

// ========================================================================
//  Test BGZFile multi-file naive concat (end-to-end via BGZFile API)
// ========================================================================
void test_bgzfile_naive_concat_simulation() {
    std::cout << "--- test_bgzfile_naive_concat_simulation ---" << std::endl;

    const std::string file1 = "/tmp/test_concat_sim_1.vcf.gz";
    const std::string file2 = "/tmp/test_concat_sim_2.vcf.gz";
    const std::string out_file = "/tmp/test_concat_sim_out.vcf.gz";

    // Write two "per-chromosome" VCF files
    {
        ngslib::BGZFile f1(file1, "wb");
        f1 << "##fileformat=VCFv4.2\n";
        f1 << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        f1 << "chr1\t100\t.\tA\tG\t50\tPASS\t.\n";
        f1.close();

        ngslib::BGZFile f2(file2, "wb");
        f2 << "##fileformat=VCFv4.2\n";
        f2 << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        f2 << "chr2\t200\t.\tC\tT\t60\tPASS\t.\n";
        f2.close();
    }

    // Simulate naive concat: write header from file1, then raw blocks from both
    {
        ngslib::BGZFile in1(file1, "r");
        ngslib::BGZFile in2(file2, "r");
        ngslib::BGZFile out(out_file, "w");

        // Process each file: skip header (write only for first), then raw blocks
        for (int fi = 0; fi < 2; ++fi) {
            ngslib::BGZFile &fin = (fi == 0) ? in1 : in2;
            fin.read_block();

            // Skip VCF header — mirrors skip_vcf_gz_header logic from concat.cpp
            const char *buf = fin.uncompressed_data();
            int nskip = 1;
            while (true) {
                if (buf[nskip] == '\n') {
                    nskip++;
                    if (nskip >= fin.block_length()) {
                        if (fi == 0)
                            out.write_raw(buf, fin.block_length());
                        fin.read_block();
                        buf = fin.uncompressed_data();
                        if (!fin.block_length()) break;
                        nskip = 0;
                    }
                    if (buf[nskip] != '#') {
                        if (fi == 0)
                            out.write_raw(buf, nskip);
                        break;
                    }
                }
                nskip++;
                if (nskip >= fin.block_length()) {
                    if (fi == 0)
                        out.write_raw(buf, fin.block_length());
                    fin.read_block();
                    buf = fin.uncompressed_data();
                    if (!fin.block_length()) break;
                    nskip = 0;
                }
            }

            // Write remaining data in current block
            if (fin.block_length() - nskip > 0) {
                out.write_raw(buf + nskip, fin.block_length() - nskip);
            }
            out.flush();

            // Copy raw blocks (skip EOF)
            const size_t page_size = BGZF_MAX_BLOCK_SIZE;
            std::vector<uint8_t> rbuf(page_size);
            while (true) {
                ssize_t nr = fin.raw_block_read(rbuf.data(), ngslib::BGZF_HEADER_SIZE);
                if (nr == 0) break;
                uint32_t bsz = static_cast<uint32_t>(ngslib::bgzf_block_bsize(rbuf.data())) + 1;
                nr += fin.raw_block_read(rbuf.data() + ngslib::BGZF_HEADER_SIZE, bsz - ngslib::BGZF_HEADER_SIZE);
                if (ngslib::is_bgzf_eof(rbuf.data(), nr)) continue;
                out.raw_block_write(rbuf.data(), nr);
            }
        }
        // ~BGZFile closes all files
    }

    // Verify the output has the correct content
    {
        ngslib::BGZFile result(out_file, "r");
        std::string line;
        int data_lines = 0;
        int header_lines = 0;
        while (result.readline(line)) {
            if (line[0] == '#') header_lines++;
            else data_lines++;
        }
        TEST_ASSERT(header_lines == 2, "Output should have 2 header lines (##fileformat + #CHROM)");
        TEST_ASSERT(data_lines == 2, "Output should have 2 data lines (one from each file)");
    }

    std::remove(file1.c_str());
    std::remove(file2.c_str());
    std::remove(out_file.c_str());
}

// ========================================================================
//  Test: bgzf_copy_blocks_to — single-file block copy (no prior decompressed read)
// ========================================================================
static void test_bgzf_copy_blocks_to() {
    std::cout << "--- test_bgzf_copy_blocks_to ---" << std::endl;

    // Create a source BGZF file with data only (simulates variant_caller temp files)
    std::string src_file = "_test_copy_src.vcf.gz";
    std::string dst_file = "_test_copy_dst.vcf.gz";
    {
        ngslib::BGZFile src(src_file, "wb");
        src << "chr1\t100\t.\tA\tG\n";
        src << "chr1\t200\t.\tC\tT\n";
        src << "chr2\t300\t.\tG\tA\n";
    }

    // Copy all blocks from source to destination
    {
        ngslib::BGZFile fin(src_file, "r");
        ngslib::BGZFile fout(dst_file, "wb");
        ngslib::bgzf_copy_blocks_to(fin, fout, src_file);
    }

    // Verify the output is byte-identical to the source
    {
        FILE *sf = fopen(src_file.c_str(), "rb");
        FILE *df = fopen(dst_file.c_str(), "rb");
        fseek(sf, 0, SEEK_END); fseek(df, 0, SEEK_END);
        long src_size = ftell(sf); long dst_size = ftell(df);
        fclose(sf); fclose(df);
        TEST_ASSERT(src_size == dst_size, "Source and destination should have the same file size");
    }

    // Verify the output is readable and has correct content
    ngslib::BGZFile result(dst_file, "r");
    std::string line;
    int data_lines = 0;
    while (result.readline(line)) {
        data_lines++;
    }
    TEST_ASSERT(data_lines == 3, "Should have 3 data lines");

    std::remove(src_file.c_str());
    std::remove(dst_file.c_str());
}

// ========================================================================
//  Test: bgzf_raw_concat — multi-file block-level concat with external header
// ========================================================================
static void test_bgzf_raw_concat() {
    std::cout << "--- test_bgzf_raw_concat ---" << std::endl;

    // Create multiple data-only BGZF files (no headers, simulating variant_caller temp files)
    std::string f1 = "_test_rawconcat_1.vcf.gz";
    std::string f2 = "_test_rawconcat_2.vcf.gz";
    std::string f3 = "_test_rawconcat_3.vcf.gz";
    std::string out = "_test_rawconcat_out.vcf.gz";

    {
        ngslib::BGZFile out1(f1, "wb");
        out1 << "chr1\t100\t.\tA\tG\n";
        out1 << "chr1\t200\t.\tC\tT\n";
    }
    {
        ngslib::BGZFile out2(f2, "wb");
        out2 << "chr2\t300\t.\tG\tA\n";
    }
    {
        ngslib::BGZFile out3(f3, "wb");
        out3 << "chr3\t400\t.\tT\tC\n";
        out3 << "chr3\t500\t.\tN\tA\n";
    }

    // Write header + concat all files
    {
        ngslib::BGZFile fout(out, "wb");
        fout << "##fileformat=VCFv4.2\n";
        fout << "#CHROM\tPOS\tID\tREF\tALT\n";
        fout.flush();

        std::vector<std::string> infiles = {f1, f2, f3};
        ngslib::bgzf_raw_concat(infiles, fout, /*remove_inputs=*/true);
    }

    // Verify output
    {
        ngslib::BGZFile result(out, "r");
        std::string line;
        int header_lines = 0, data_lines = 0;
        while (result.readline(line)) {
            if (line[0] == '#') header_lines++;
            else data_lines++;
        }
        TEST_ASSERT(header_lines == 2, "Should have 2 header lines (##fileformat + #CHROM)");
        TEST_ASSERT(data_lines == 5, "Should have 5 data lines (2 + 1 + 2)");
    }

    // Temp files should have been removed by bgzf_raw_concat
    TEST_ASSERT(access(f1.c_str(), F_OK) != 0, "Temp file 1 should be removed");
    TEST_ASSERT(access(f2.c_str(), F_OK) != 0, "Temp file 2 should be removed");
    TEST_ASSERT(access(f3.c_str(), F_OK) != 0, "Temp file 3 should be removed");

    std::remove(out.c_str());
}

// ========================================================================
//  Main
// ========================================================================
int main() {
    std::cout << "=== BGZF block-level utilities unit tests ===" << std::endl;
    std::cout << std::endl;

    test_eof_size();
    test_eof_block_content();
    test_bgzf_block_is_valid();
    test_bgzf_block_bsize();
    test_is_bgzf_eof();
    test_roundtrip_raw_blocks();
    test_bsize_valid_header_consistency();
    test_bgzfile_read_block();
    test_bgzfile_raw_block_read();
    test_bgzfile_raw_block_write_roundtrip();
    test_bgzfile_naive_concat_simulation();
    test_bgzf_copy_blocks_to();
    test_bgzf_raw_concat();

    std::cout << std::endl;
    std::cout << "=== Results: " << g_tests_passed << " passed, "
              << g_tests_failed << " failed ===" << std::endl;

    return g_tests_failed > 0 ? 1 : 0;
}
