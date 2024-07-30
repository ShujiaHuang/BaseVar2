#include <stdexcept>

#include <htslib/hts.h>
#include "bam_header.h"
#include "utils.h"

namespace ngslib {

    BamHeader::BamHeader(samFile *fp) {
        _h = sam_hdr_read(fp);
    }

    BamHeader::BamHeader(const std::string &fn) {

        if (!is_readable(fn)) {
            throw std::invalid_argument("_bam_header::BamHeader: " + fn + "not found.");
        }

        samFile *fp = hts_open(fn.c_str(), "r");
        _h = sam_hdr_read(fp);  // get a BAM header pointer on success, NULL on failure.
        sam_close(fp);
    }

    BamHeader &BamHeader::operator=(const BamHeader &bh) {

        // release _h pointer if _h is not NULL
        sam_hdr_destroy(_h);
        _h = sam_hdr_dup(bh._h);
        return *this;
    }

    BamHeader &BamHeader::operator=(const sam_hdr_t *hdr) {

        // release _h pointer if _h is not NULL
        sam_hdr_destroy(_h);
        _h = sam_hdr_dup(hdr);
        return *this;
    }

    std::ostream &operator<<(std::ostream &os, const BamHeader &hd) {

        if (hd._h)
            os << sam_hdr_str(hd._h);

        return os;
    }

    void BamHeader::destroy() {
        sam_hdr_destroy(_h);
        _h = NULL;
    }

    int BamHeader::name2id(const std::string &name) {
        int tid = sam_hdr_name2tid(_h, name.c_str());

        if (tid < 0) {
            throw std::invalid_argument(
                    "[bam_header.cpp::BamHeader:name2id] Unknown reference name or "
                    "the header not be parsed: " + name);
        }
        return tid;
    }
}  // namespace ngslib