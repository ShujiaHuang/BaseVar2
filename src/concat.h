/**
 * @file concat.h
 * @brief Concat the BaseVar VCF and CVG files.
 * 
 * @author Shujia Huang
 * @date 2018-09-14
 * 
 */

#ifndef __INCLUDE_BASEVAR_CONCAT_H__
#define __INCLUDE_BASEVAR_CONCAT_H__

#include <getopt.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "caller_utils.h"  // merge_file_by_line
#include "io/iobgzf.h"
#include "io/utils.h"      // ngslib::get_firstcolumn_from_file

int _concat_basevar_outfile(const std::vector<std::string> &infiles, const std::string outfile);
int concat_runner(int argc, char *argv[]);

#endif
