/**
 * @file basetype_utils.h
 * 
 * @brief define the global variable and util codes specific for basetype.h
 *  
 * @author Shujia Huang
 * @date 2018-08-29
 * 
 */

#ifndef __INCLUDE_BASETYPE_UTILS_H__
#define __INCLUDE_BASETYPE_UTILS_H__

#include <string>
#include <vector>

// Getting the first column from input file, this it's used for getting 
// filename from input filelist.
std::vector<std::string> get_firstcolumn_from_file(const std::string fn);

/// Header for VCF
std::string vcf_header_define(const std::string &ref_file_path, const std::vector<std::string> &addition_info, 
                              const std::vector<std::string> &samples);
std::string cvg_header_define(const std::vector<std::string> &group_info, const std::vector<char> &BASES);
void merge_file_by_line(const std::vector<std::string> &infiles, const std::string &outfile, 
                        std::string header="#", bool is_remove_tempfile=false);

#endif