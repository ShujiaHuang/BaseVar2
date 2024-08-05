// Author: Shujia Huang
// Date: 2021-08-25
#include <iostream>
#include <string>

#include <htslib/sam.h>

#include "bam_header.h"
#include "bam_record.h"
#include "utils.h"


int main() {
    using ngslib::BamHeader;
    using ngslib::BamRecord;

    std::string fn1 = "../data/range.bam";
    std::string fn2 = "../data/range.cram";
    std::string fn3 = "../data/xx_MD.bam";
    std::string fn4 = "../data/xx_minimal.sam";

//    samFile *fp = sam_open(fn1.c_str(), "r");
    samFile *fp = sam_open(fn2.c_str(), "r");   // cram
//    samFile *fp = sam_open(fn3.c_str(), "r");   // bam
//    samFile *fp = sam_open(fn4.c_str(), "r");    // sam
    BamHeader hdr = BamHeader(fp);
    bam1_t *al = bam_init1();

    BamRecord br0;
    BamRecord br1;
    BamRecord br2 = br1;
    BamRecord br3;
    br3.init();
    BamRecord br4 = al;

    int read_count = 0;
    std::cout << hdr << "\n";
    while (br3.load_read(fp, hdr.h()) >= 0) {

        std::cout << br3 << "; bool: " << bool(br3) << "\n";
        std::cout << " * Read count: " << ++read_count

                  << "; align_length: " << br3.align_length()
                  << "; match_length('M'): " << br3.match_length()
                  << "; Read name: " << br3.qname()
                  << "; Read length: " << br3.query_length()
                  << "; query_sequence: " << br3.query_sequence()
                  << "; query_qual: " << br3.query_qual()
                  << "; mean_qqual_phred: " << br3.mean_qqual()
                  << "; query_start_pos: " << br3.query_start_pos()
                  << "; query_start_pos_reverse: " << br3.query_start_pos_reverse()
                  << "; query_end_pos: " << br3.query_end_pos()
                  << "; query_end_pos_reverse: " << br3.query_end_pos_reverse()
                  << "; read_group: " << br3.read_group()
                  << "; get_tag(NM): " << br3.get_tag("NM")
                  << "; get_tag(MD): " << br3.get_tag("MD")
                  << "; get_tag(XT): " << br3.get_tag("XT")
                  << "; get_tag(RG): " << br3.get_tag("RG")

                  << "; FLAG: " << br3.flag()
                  << "; is_paired: " << br3.is_paired()
                  << "; is_mapped: " << br3.is_mapped()
                  << "; is_mate_mapped: " << br3.is_mate_mapped()
                  << "; is_mapped_reverse: " << br3.is_mapped_reverse()
                  << "; is_mate_mapped_reverse: " << br3.is_mate_mapped_reverse()
                  << "; target_id: " << br3.tid()
                  << "; mate_target_id: " << br3.mate_tid()
                  << "; tid_name: " << br3.tid_name(hdr)
                  << "; mate_tid_name: " << br3.mate_tid_name(hdr)
                  << "; tid_length: " << br3.tid_length(hdr)
                  << "; mapped_strand: " << br3.map_strand()
                  << "; mate_mapped_strand: " << br3.mate_map_strand()

                  << "; mapped begin_position: " << br3.reference_start_pos()
                  << "; mapped end position: " << br3.reference_end_pos()
                  << "; mate mapped position: " << br3.mate_reference_start_pos()
                  << "; mapping quality: " << br3.mapq()

                  << "; insert-size: " << br3.insert_size()
                  << "; is_read1: " << br3.is_read1()
                  << "; is_read2: " << br3.is_read2()
                  << "; is_proper_pair: " << br3.is_proper_pair()
                  << "; is_secondary: " << br3.is_secondary()
                  << "; is_qc_fail: " << br3.is_qc_fail()
                  << "; is_duplicate: " << br3.is_duplicate()
                  << "; is_supplementary: " << br3.is_supplementary()

                  << "; proper_orientation: " << br3.is_proper_orientation()
                  << "; br0.is_mapped(): " << br0.is_mapped() << "\n";
    }

    br3.set_qc_fail();

    sam_close(fp);
    return 0;
}













