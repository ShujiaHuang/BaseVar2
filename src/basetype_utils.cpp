#include <htslib/bgzf.h>
#include <htslib/kstring.h>

#include "basetype_utils.h"
#include "utils.h"

void merge_file_by_line(const std::vector<std::string> &infiles, const std::string &outfile, 
                        std::string header) 
{
    BGZF *OUT = bgzf_open(outfile.c_str(), "w");  // output file
    if (!OUT) throw std::runtime_error("[ERROR] " + outfile + " open failure.");

    header += "\n";
    if (bgzf_write(OUT, header.c_str(), header.length()) != header.length())
        throw std::runtime_error("[ERROR] fail to write data");

    /* Merge all files here */
    for (auto fn: infiles) {
        BGZF *f = bgzf_open(fn.c_str(), "r");
        kstring_t s; s.s = NULL; s.l = s.m = 0;
        while (bgzf_getline(f, '\n', &s) >= 0) {
            if (s.s[0] == '#') continue;  // ignore the header of subfiles.
            std::string out(s.s); out += "\n";

            if (bgzf_write(OUT, out.c_str(), out.length()) != out.length())
                throw std::runtime_error("[ERROR] fail to write data");
        }
        bgzf_close(f);
        ngslib::safe_remove(fn);
    }
    
    int is_cl = bgzf_close(OUT);
    if (is_cl < 0) throw std::runtime_error("[ERROR] " + outfile + " fail close.");
    
    return;
}

