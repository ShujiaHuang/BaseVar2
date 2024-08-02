/**
 * @file basetype.cpp
 * @brief 
 * 
 * @author Shujia Huang
 * @date 2018-08-01
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include <sstream>
#include "basetype.h"
#include "utils.h"


void BaseTypeRunner::set_arguments(int cmdline_argc, char *cmdline_argv[]) {

    if (!args) {
        throw std::invalid_argument("[basetype.cpp::BaseTypeRunner:args] 'args' must be a "
                                    "NULL pointer before it can be assigned a value.");
    }
    args = new BaseTypeArgs;  // set it to be a new BasTypeArgs stucture type.
    
    char c;
    int opt_idx;
    while((c = getopt_long(cmdline_argc, cmdline_argv, "I:L:R:m:q:B:t:r:p:G:h", BASETYPE_CMDLINE_LOPTS, &opt_idx)) >= 0) {
        std::stringstream ss(optarg ? optarg: "");  // 用于解决字符串转浮点等其他类型的转换问题
        switch (c) {
            case 'I': args->input_bf.push_back(optarg);         break;
            case 'L': args->in_bamfilelist = optarg;            break;
            case 'R': args->reference      = optarg;            break;

            case 'm': ss >> args->min_af;                       break;
            case 'q': ss >> args->mapq;                         break;
            case 'B': ss >> args->batchcount;                   break;
            case 't': ss >> args->thread_num;                   break;

            case 'r': args->regions = optarg;                   break;
            case 'p': args->in_pos_file = optarg;               break;
            case 'G': args->pop_group_file = optarg;            break;
            case '1': args->output_vcf = optarg;                break;
            case '2': args->output_cvg = optarg;                break;

            case '3': args->filename_has_samplename = true;     break;
            case '4': args->smart_rerun = true;                 break;
            case 'h': 
                std::cerr << usage() << std::endl;
                exit(1);

            default: 
                std::cerr << "Unknown argument: " << c << std::endl; 
                exit(1);
        }
    }
    // Output the commandline options
    std::cerr << 
        "[INFO] BaseVar commandline:\n"
        "basevar basetype " + (args->input_bf.empty() ? "" : "-I " + ngslib::join(args->input_bf, " -I ")) +
        (args->in_bamfilelist.empty() ? "": " -L " + args->in_bamfilelist) + " \\ \n"
        "   -R " + args->reference            + " \\ \n"
        "   -q " << args->mapq               << " \\ \n"
        "   -m " << args->min_af             << " \\ \n"
        "   -B " << args->batchcount         << " \\ \n"
        "   -t " << args->thread_num         << " \\ \n"  << (args->regions.empty() ? "": 
        "   -r " + args->regions              + " \\ \n") << (args->in_pos_file.empty() ? "": 
        "   -p " + args->in_pos_file          + " \\ \n") << (args->pop_group_file.empty() ? "": 
        "   -p " + args->pop_group_file       + " \\ \n") <<
        "   --output-vcf " + args->output_vcf + " \\ \n"
        "   --output-vcg " + args->output_cvg <<
        (args->filename_has_samplename ? "\\ \n--filename-has-samplename": "") << 
        (args->smart_rerun ? "\\ \n--smart-rerun": "") << "\n" << std::endl;
}

