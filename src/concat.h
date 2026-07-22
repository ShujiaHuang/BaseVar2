/**
 * @file concat.h
 * @brief Concatenate or combine BaseVar VCF files.
 *
 * Supports two modes:
 *   - Default: line-by-line decompression/recompression.
 *   - Naive (--naive): BGZF block-level raw concatenation (extremely fast).
 *
 * @author Shujia Huang
 * @date 2018-09-14
 *
 */

#ifndef __INCLUDE_BASEVAR_CONCAT_H__
#define __INCLUDE_BASEVAR_CONCAT_H__

#include <string>
#include <vector>

int _concat_basevar_outfile(const std::vector<std::string> &infiles, const std::string outfile);
int concat_runner(int argc, char *argv[]);

#endif
