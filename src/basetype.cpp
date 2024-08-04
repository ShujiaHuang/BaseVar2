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
#include <fstream>

#include "basetype.h"
#include "utils.h"


void BaseTypeRunner::set_arguments(int cmdline_argc, char *cmdline_argv[]) {

    if (_args) {
        throw std::invalid_argument("[basetype.cpp::BaseTypeRunner:args] 'args' must be "
                                    "a NULL pointer before it can be assigned a value.");
    }
    _args = new BaseTypeArgs;  // set it to be a new BasTypeArgs stucture type.
    
    char c;
    int opt_idx;
    /**
     * @brief Parsing the commandline options 
     * 
     */
    while((c = getopt_long(cmdline_argc, cmdline_argv, "I:L:R:m:q:B:t:r:p:G:h", BASETYPE_CMDLINE_LOPTS, &opt_idx)) >= 0) {
        std::stringstream ss(optarg ? optarg: "");  // 字符流解决命令行参数转浮点等类型的问题
        switch (c) {
            case 'I': _args->input_bf.push_back(optarg);         break;
            case 'L': _args->in_bamfilelist = optarg;            break;
            case 'R': _args->reference      = optarg;            break;

            case 'm': ss >> _args->min_af;                       break;
            case 'q': ss >> _args->mapq;                         break;
            case 'B': ss >> _args->batchcount;                   break;
            case 't': ss >> _args->thread_num;                   break;

            case 'r': _args->regions = optarg;                   break;
            case 'p': _args->in_pos_file = optarg;               break;
            case 'G': _args->pop_group_file = optarg;            break;
            case '1': _args->output_vcf = optarg;                break;
            case '2': _args->output_cvg = optarg;                break;

            case '3': _args->filename_has_samplename = true;     break;
            case '4': _args->smart_rerun = true;                 break;
            case 'h': 
                std::cerr << usage() << std::endl;
                exit(1);

            default: 
                std::cerr << "Unknown argument: " << c << std::endl; 
                exit(1);
        }
    }

    if (_args->input_bf.empty() && _args->in_bamfilelist.empty()) {
        throw std::invalid_argument("[ERROR] Missing argument '-I/--input' or '-L/--align-file-list'");
    }

    // Output the commandline options
    std::cerr << 
        "[INFO] BaseVar commandline:\n"
        "basevar basetype " + (_args->input_bf.empty() ? "" : "-I " + ngslib::join(_args->input_bf, " -I ")) +
        (_args->in_bamfilelist.empty() ? "": " -L " + _args->in_bamfilelist) + " \\ \n"
        "   -R " + _args->reference            + " \\ \n"
        "   -q " << _args->mapq               << " \\ \n"
        "   -m " << _args->min_af             << " \\ \n"
        "   -B " << _args->batchcount         << " \\ \n"
        "   -t " << _args->thread_num         << " \\ \n"  << (_args->regions.empty() ? "": 
        "   -r " + _args->regions              + " \\ \n") << (_args->in_pos_file.empty() ? "": 
        "   -p " + _args->in_pos_file          + " \\ \n") << (_args->pop_group_file.empty() ? "": 
        "   -p " + _args->pop_group_file       + " \\ \n") <<
        "   --output-vcf " + _args->output_vcf + " \\ \n"
        "   --output-vcg " + _args->output_cvg << (_args->filename_has_samplename ? " \\ \n"
        "   --filename-has-samplename": "")   << (_args->smart_rerun ? " \\ \n"
        "   --smart-rerun": "") << "\n" << std::endl;
    
    // 
    if (!_args->in_bamfilelist.empty()) _load_bamfile_list();
    _reference = _args->reference;  // load fasta
    _load_calling_interval();
}

void BaseTypeRunner::_load_bamfile_list() {
    std::ifstream i_fn(_args->in_bamfilelist.c_str());
    if (!i_fn) {
        std::cerr << "[ERROR] Cannot open file: " + _args->in_bamfilelist << std::endl;
        exit(1);
    }

    std::string tmp, fn;
    while (1) {
        i_fn >> fn;
        if (i_fn.eof()) break;
        _args->input_bf.push_back(fn);

        std::getline(i_fn, tmp, '\n');
    }
    i_fn.close();
    std::cerr << ngslib::join(_args->input_bf, ",") << std::endl;
}

void BaseTypeRunner::_load_calling_interval() {

    if (!_args->regions.empty()) {
        std::vector<std::string> region_vector;
        ngslib::split(_args->regions, region_vector, ",");

        for (size_t i(0); i < region_vector.size(); ++i) {
            _calling_intervals.push_back(_make_gregiontuple(region_vector[i]));
        }
    }

    if (!_args->in_pos_file.empty()) {

    }

    return;
}

ngslib::GenomeRegionTuple BaseTypeRunner::_make_gregiontuple(std::string gregion) {

    std::string ref_id;
    uint32_t pos_start, pos_end;  // All be 1-based

    std::vector<std::string> gr1;
    std::vector<uint32_t> gr2;

    ngslib::split(gregion, gr1, ":");     // get reference id
    ref_id = gr1[0];

    if (gr1.size() == 2) {                // 'start-end' or start
        ngslib::split(gr1[1], gr2, "-");  // get position coordinate
        pos_start = gr2[0];
        pos_end = (gr2.size() == 2) ? gr2[1] : _reference.seq_length(ref_id);
    } else {
        pos_start = 1;
        pos_end = _reference.seq_length(ref_id);  // the whole ``ref_id`` length
    }

    return std::make_tuple(ref_id, pos_start, pos_end);  // 1-based
}
