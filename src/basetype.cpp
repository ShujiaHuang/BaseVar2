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
        throw std::invalid_argument("[basetype.cpp::BaseTypeRunner:args] muse be a NULL pointer.");
    }
    args = new BaseTypeArgs;
    
    char c;
    int opt_idx;
    while((c = getopt_long(cmdline_argc, cmdline_argv, "I:L:R:m:q:B:t:r:p:G", BASETYPE_CMDLINE_LOPTS, &opt_idx)) >= 0) {
        std::stringstream ss(optarg ? optarg: "");  // 用于解决非字符串类型输入的类型转换问题
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

            case '3': ss >> args->filename_has_samplename;      break;
            case '4': ss >> args->smart_rerun;                  break;
            case '?': std::cerr << BASETYPE_USAGE << std::endl; break;

            default: 
                std::cerr << "Unknown argument: " << c << std::endl; 
                exit(1);
        }
    }

    std::cerr << args->input_bf.size() << ", " << args->input_bf[0] << "\n" 
              << args->min_af << ", " << args->mapq << ", " << args->batchcount << "\n";
}

