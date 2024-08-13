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

        std::cout << br3 << "\n" 
                  << " - Success:                 " << bool(br3) << "\n"
                  << " - Read count:              " << ++read_count << "\n"
                  << " - align_length:            " << br3.align_length() << "\n"
                  << " - match_length('M'):       " << br3.match_length() << "\n"
                  << " - mapping quality:         " << br3.mapq() << "\n"
                  << " - Read name:               " << br3.qname() << "\n"
                  << " - Read length:             " << br3.query_length() << "\n"
                  << " - query_sequence:          " << br3.query_sequence() << "\n"
                  << " - query_qual:              " << br3.query_qual() << "\n"
                  << " - mean_qqual_phred:        " << br3.mean_qqual() << "\n"
                  << " - query_start_pos:         " << br3.query_start_pos() << "\n"
                  << " - query_end_pos:           " << br3.query_end_pos() << "\n"
                  << " - query_start_pos_reverse: " << br3.query_start_pos_reverse() << "\n"
                  << " - query_end_pos_reverse:   " << br3.query_end_pos_reverse()   << "\n"
                  << " - read_group:              " << br3.read_group()  << "\n"
                  << " - get_tag(NM):             " << br3.get_tag("NM") << "\n"
                  << " - get_tag(MD):             " << br3.get_tag("MD") << "\n"
                  << " - get_tag(XT):             " << br3.get_tag("XT") << "\n"
                  << " - get_tag(RG):             " << br3.get_tag("RG") << "\n"

                  << " - FLAG:                    " << br3.flag() << "\n"
                  << " - is_paired:               " << br3.is_paired() << "\n"
                  << " - is_mapped:               " << br3.is_mapped() << "\n"
                  << " - is_mate_mapped:          " << br3.is_mate_mapped() << "\n"
                  << " - is_mapped_reverse:       " << br3.is_mapped_reverse() << "\n"
                  << " - is_mate_mapped_reverse:  " << br3.is_mate_mapped_reverse() << "\n"
                  << " - target_id:               " << br3.tid() << "\n"
                  << " - tid_name:                " << br3.tid_name(hdr) << "\n"
                  << " - tid_length:              " << br3.tid_length(hdr) << "\n"
                  << " - mate_target_id:          " << br3.mate_tid() << "\n"
                  << " - mate_tid_name:           " << br3.mate_tid_name(hdr) << "\n"
                  << " - mapped_strand:           " << br3.map_strand() << "\n"
                  << " - mapped begin_position:   " << br3.map_ref_start_pos() << "\n"
                  << " - mapped end position:     " << br3.map_ref_end_pos() << "\n"
                  << " - mate_mapped_strand:      " << br3.mate_map_strand() << "\n"
                  << " - mate mapped position:    " << br3.mate_map_ref_start_pos() << "\n"

                  << " - insert-size:             " << br3.insert_size() << "\n"
                  << " - is_read1:                " << br3.is_read1() << "\n"
                  << " - is_read2:                " << br3.is_read2() << "\n"
                  << " - is_proper_pair:          " << br3.is_proper_pair() << "\n"
                  << " - is_secondary:            " << br3.is_secondary() << "\n"
                  << " - is_qc_fail:              " << br3.is_qc_fail() << "\n"
                  << " - is_duplicate:            " << br3.is_duplicate() << "\n"
                  << " - is_supplementary:        " << br3.is_supplementary() << "\n"

                  << " - proper_orientation:      " << br3.is_proper_orientation() << "\n"
                  << " - br0.is_mapped():         " << br0.is_mapped() << "\n\n";
    }

    br3.set_qc_fail();

    sam_close(fp);
    return 0;
}













