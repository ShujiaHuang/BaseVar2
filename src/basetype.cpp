/**
 * @file basetype.cpp
 * @brief 
 * 
 * @author Shujia Huang
 * @date 2018-08-01
 * 
 */
#include <sstream>
#include <fstream>

#include "basetype.h"
#include "bam_header.h"


void BaseTypeRunner::set_arguments(int cmdline_argc, char *cmdline_argv[]) {

    if (_args) {
        throw std::invalid_argument("[basetype.cpp::BaseTypeRunner:args] 'args' must be "
                                    "a NULL pointer before it can be assigned a value.");
    }
    // Inital a new BasTypeARGS and set defaut argument.
    _args = new BaseTypeARGS;
    
    // Parsing the commandline options. 
    char c;
    while((c = getopt_long(cmdline_argc, cmdline_argv, "I:L:R:m:q:B:t:r:G:h", BASETYPE_CMDLINE_LOPTS, NULL)) >= 0)
    {
        std::stringstream ss(optarg ? optarg: "");  // 字符流解决命令行参数转浮点等类型的问题
        switch (c) {
            case 'I': _args->input_bf.push_back(optarg);         break;  // 恒参，一直用
            case 'L': _args->in_bamfilelist = optarg;            break;  /* 临参 */
            case 'R': _args->reference      = optarg;            break;  /* 临参 */

            case 'm': ss >> _args->min_af;                       break;  // 恒参
            case 'q': ss >> _args->mapq;                         break;  // 恒参
            case 'B': ss >> _args->batchcount;                   break;  // 恒参
            case 't': ss >> _args->thread_num;                   break;  // 恒参

            case 'r': _args->regions = optarg;                   break;  /* 临参 */
            case 'G': _args->pop_group_file = optarg;            break;  // 
            case '1': _args->output_vcf = optarg;                break;  // 恒参
            case '2': _args->output_cvg = optarg;                break;  // 恒参

            case '3': _args->filename_has_samplename = true;     break;  // 恒参
            case '4': _args->smart_rerun = true;                 break;  // 恒参
            case 'h': 
                std::cout << usage() << std::endl;
                exit(1);

            default: 
                std::cerr << "Unknown argument: " << c << std::endl; 
                exit(1);
        }
    }
    // check the requirement argument.
    /* Make sure you have set at least one bamfile. */
    if (_args->input_bf.empty() && _args->in_bamfilelist.empty())
        throw std::invalid_argument("[ERROR] Missing argument '-I/--input' or '-L/--align-file-list'");
    if (_args->reference.empty())
        throw std::invalid_argument("[ERROR] Missing argument '-R/--reference'");

    if (_args->output_vcf.empty())
        throw std::invalid_argument("[ERROR] Missing argument '--output-vcf'");
    if (_args->output_cvg.empty())
        throw std::invalid_argument("[ERROR] Missing argument '--output-cvg'");

    // Output the commandline options
    std::cout << 
        "[INFO] BaseVar commandline:\n"
        "basevar basetype " + (_args->input_bf.empty() ? "" : "-I " + ngslib::join(_args->input_bf, " -I ")) +
        (_args->in_bamfilelist.empty() ? "": "-L " + _args->in_bamfilelist) + " \\ \n"
        "   -R " + _args->reference            + " \\ \n"
        "   -q " << _args->mapq               << " \\ \n"
        "   -m " << _args->min_af             << " \\ \n"
        "   -B " << _args->batchcount         << " \\ \n"
        "   -t " << _args->thread_num         << " \\ \n"  << (_args->regions.empty() ? "": 
        "   -r " + _args->regions              + " \\ \n") << (_args->pop_group_file.empty() ? "": 
        "   -p " + _args->pop_group_file       + " \\ \n") <<
        "   --output-vcf " + _args->output_vcf + " \\ \n"
        "   --output-vcg " + _args->output_cvg << (_args->filename_has_samplename ? " \\ \n"
        "   --filename-has-samplename": "")    << (_args->smart_rerun ? " \\ \n"
        "   --smart-rerun": "") << "\n" << std::endl;
    
    if (_args->smart_rerun) {
        std::cout << "************************************************\n"
                     "******************* WARNING ********************\n"
                     "************************************************\n"
                     ">>>>>>>> You have setted `smart rerun` <<<<<<<<<\n"
                     "Please make sure all the parameters are the same\n"
                     "with your previous commands.\n"
                     "************************************************\n"
                     "************************************************\n\n";
    }

    if (!_args->in_bamfilelist.empty()) _get_bamfile_list();
    std::cout << "[INFO] Finish loading arguments and we have " << _args->input_bf.size()
              << " BAM/CRAM files for variants calling.\n";

    // Setting the resolution of AF
    if (_args->min_af > 100.0/_args->input_bf.size())
        _args->min_af = 100.0/_args->input_bf.size();

    reference = _args->reference;  // load fasta
    _get_calling_interval();
    print_calling_interval();
    _get_sample_id_from_bam();     // get sample id from aligne_files and record into `_samples_id`

    return;
}

void BaseTypeRunner::_get_bamfile_list() {
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

    return;
}

void BaseTypeRunner::_get_sample_id_from_bam() {
    // Loading sample ID in BAM/CRMA files from RG tag.
    std::cout << "[INFO] Start loading all samples' id from alignment files.\n";
    if (_args->filename_has_samplename)
        std::cout << "[INFO] loading samples' id from filename becuase you set "
                     "--filname-has-samplename.\n";

    std::string samplename, filename;
    size_t si;
    for (size_t i(0); i < _args->input_bf.size(); ++i) {

        if ((i+1) % 1000 == 0)
            std::cout << "[INFO] loading " << i+1 << "/" << _args->input_bf.size() 
                      << " alignment files.\n";
        
        if (_args->filename_has_samplename) {
            filename = ngslib::remove_filename_extension(ngslib::basename(_args->input_bf[i]));
            si = filename.find_first_of('.');
            samplename = si > 0 && si != std::string::npos ? filename.substr(0, si) : filename;
        } else {
            // Get sampleID from BAM header, maybe time-consuming.
            ngslib::BamHeader bh(_args->input_bf[i]);
            samplename = bh.get_sample_name();
        }

        if (!samplename.empty()) {
            _samples_id.push_back(samplename);
        } else {
            throw std::invalid_argument("[BaseTypeRunner::_load_sample_id_from_bam] " + 
                                        _args->input_bf[i] + " sample ID not found.\n");
        }
    }

    return;
}

void BaseTypeRunner::_get_calling_interval() {

    // clear
    _calling_intervals.clear();

    if (!_args->regions.empty()) {
        std::vector<std::string> rg_v;
        ngslib::split(_args->regions, rg_v, ",");

        for (size_t i(0); i < rg_v.size(); ++i) {
            _calling_intervals.push_back(_make_gregiontuple(rg_v[i]));
        }
    } else {
        // calling the whole genome
        int n = reference.nseq();
        for (size_t i(0); i < n; ++i) {
            std::string ref_id = reference.iseq_name(i);
            _calling_intervals.push_back(std::make_tuple(ref_id, 1, reference.seq_length(ref_id)));
        } 
    }
    
    return;
}

void BaseTypeRunner::print_calling_interval() {

    std::cout << "---- Calling Intervals ----\n";
    for (size_t i(0); i < _calling_intervals.size(); ++i) {
        std::cout << i+1 << " - " 
                  << std::get<0>(_calling_intervals[i]) << ":" 
                  << std::get<1>(_calling_intervals[i]) << "-" 
                  << std::get<2>(_calling_intervals[i]) << "\n";
    }
    std::cout << "\n";
    return;
}

ngslib::GenomeRegionTuple BaseTypeRunner::_make_gregiontuple(std::string gregion) {

    std::string ref_id;
    uint32_t pos_start, pos_end;         // All be 1-based

    std::vector<std::string> gr;
    ngslib::split(gregion, gr, ":");     // get reference id
    ref_id = gr[0];

    if (gr.size() == 2) {                // 'start-end' or start
        std::vector<uint32_t> gs;
        ngslib::split(gr[1], gs, "-");   // get position coordinate
        pos_start = gs[0];
        pos_end = (gs.size() == 2) ? gs[1] : reference.seq_length(ref_id);
    } else {
        pos_start = 1;
        pos_end = reference.seq_length(ref_id);  // the whole ``ref_id`` length
    }

    return std::make_tuple(ref_id, pos_start, pos_end);  // 1-based
}
