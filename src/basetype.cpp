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
#include "bam.h"
#include "external/threadpool.h"


void BaseTypeRunner::set_arguments(int cmdline_argc, char *cmdline_argv[]) {

    if (cmdline_argc < 2) {
        std::cout << usage() << "\n" << std::endl;
        exit(1);
    }

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
    /* Make sure we have set at least one bamfile. */
    if (_args->input_bf.empty() && _args->in_bamfilelist.empty())
        throw std::invalid_argument("[ERROR] Missing argument '-I/--input' or '-L/--align-file-list'");
    if (_args->reference.empty())
        throw std::invalid_argument("[ERROR] Missing argument '-R/--reference'");

    if (_args->output_vcf.empty())
        throw std::invalid_argument("[ERROR] Missing argument '--output-vcf'");
    if (_args->output_cvg.empty())
        throw std::invalid_argument("[ERROR] Missing argument '--output-cvg'");
    
    if (_args->min_af <= 0)
        throw std::invalid_argument("[ERROR] '-m/--min-af' argument must be > 0.");
    if (_args->mapq <= 0)
        throw std::invalid_argument("[ERROR] '-q/--mapq' argument must be > 0.");
    if (_args->batchcount <= 0)
        throw std::invalid_argument("[ERROR] '-B/--batch-count' argument must be > 0.");
    if (_args->thread_num <= 0)
        throw std::invalid_argument("[ERROR] '-t/--thread' argument must be > 0.");

    // recovering the absolute paths for output files
    _args->output_vcf = ngslib::abspath(_args->output_vcf);
    _args->output_cvg = ngslib::abspath(_args->output_cvg);

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

    _get_sample_id_from_bam();     // get sample id from input aligne_files and record into `_samples_id`
    if (!_args->pop_group_file.empty()) _get_popgroup_info();

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
            _calling_intervals.push_back(_make_gregion_tuple(rg_v[i]));
        }
    } else {
        // Calling the whole genome
        int n = reference.nseq();
        for (size_t i(0); i < n; ++i) {
            std::string ref_id = reference.iseq_name(i);
            _calling_intervals.push_back(std::make_tuple(ref_id, 1, reference.seq_length(ref_id)));
        } 
    }
    
    return;
}

void BaseTypeRunner::print_calling_interval() {

    std::string ref_id;
    uint32_t reg_start, reg_end;
    std::cout << "---- Calling Intervals ----\n";
    for (size_t i(0); i < _calling_intervals.size(); ++i) {
        std::tie(ref_id, reg_start, reg_end) = _calling_intervals[i];
        std::cout << i+1 << " - " << ref_id << ":" << reg_start << "-" << reg_end << "\n";
    }
    std::cout << "\n";
    return;
}

void BaseTypeRunner::_get_popgroup_info() {
    // group_id => [index in _samples_id of group_id]
    std::ifstream i_fn(_args->pop_group_file.c_str());
    if (!i_fn) {
        std::cerr << "[ERROR] Cannot open file: " + _args->pop_group_file << std::endl;
        exit(1);
    }

    std::map<std::string, std::string> sample2group;
    std::string tmp, sn, gn;
    while (1) {
        // Only two columns: sample_id and group_id
        i_fn >> sn >> gn;
        if (i_fn.eof()) break;
        
        sample2group[sn] = gn;
        std::getline(i_fn, tmp, '\n');
    }
    i_fn.close();

    _groups_idx.clear();
    std::map<std::string, std::string>::iterator s2g_it;
    
    // follow the order of samples
    for (size_t i(0); i < _samples_id.size(); ++i) {
        s2g_it = sample2group.find(_samples_id[i]);
        
        // ignore all the samples which not found
        if (s2g_it != sample2group.end()) {
            // record sample index of group groups
            // group -> index of _samples_id
            _groups_idx[sample2group[_samples_id[i]]].push_back(i);
        }
    }
    // test
    // for(std::map<std::string, std::vector<size_t>>::iterator it(_groups_idx.begin()); it != _groups_idx.end(); ++it){
    //     std::cout << " - " << it->first << " " << it->second.size() << " : " << ngslib::join(it->second, ",") << std::endl;
    // }
    return;
}

ngslib::GenomeRegionTuple BaseTypeRunner::_make_gregion_tuple(std::string gregion) {

    // Genome Region
    std::string ref_id;
    uint32_t reg_start, reg_end;         // All be 1-based

    std::vector<std::string> gr;
    ngslib::split(gregion, gr, ":");     // get reference id
    ref_id = gr[0];

    if (gr.size() == 2) {                // 'start-end' or start
        std::vector<uint32_t> gs;
        ngslib::split(gr[1], gs, "-");   // get position coordinate
        reg_start = gs[0];
        reg_end = (gs.size() == 2) ? gs[1] : reference.seq_length(ref_id);
    } else {
        reg_start = 1;
        reg_end = reference.seq_length(ref_id);  // the whole ``ref_id`` length
    }

    if (reg_start > reg_end) {
        throw std::invalid_argument("[ERROR] start postion is larger than end position in "
                                    "'-r/--regions' " + gregion);
    }

    return std::make_tuple(ref_id, reg_start, reg_end);  // 1-based
}

bool __create_single_batchfile(const std::vector<std::string> batch_align_files, 
                               const std::string &fa_seq,
                               ngslib::GenomeRegionTuple genome_region,
                               std::string out_batch_file,  // output batchfile name
                               uint32_t reg_expand_size=500)
{
    std::string ref_id;
    uint32_t reg_start, reg_end;  // all be 1-based
    std::tie(ref_id, reg_start, reg_end) = genome_region;

    uint32_t exp_reg_start = reg_start > reg_expand_size ? reg_start - reg_expand_size : 1;
    uint32_t exp_reg_end   = reg_end + reg_expand_size;
    std::string exp_regstr = ref_id + ":" + ngslib::tostring(exp_reg_start) + "-" + ngslib::tostring(exp_reg_end);

    return true;
}

std::vector<std::string> BaseTypeRunner::_create_batchfiles(ngslib::GenomeRegionTuple genome_region) {
    // multiple thread 
    std::string outdir = ngslib::dirname(_args->output_vcf) ;
    std::string stem_bn = ngslib::remove_filename_extension(ngslib::basename(_args->output_vcf));
    std::string cache_outdir = outdir + "/" + stem_bn + "_cache";
    ngslib::safe_mkdir(cache_outdir);

    std::string ref_id;
    uint32_t reg_start, reg_end; 
    std::tie(ref_id, reg_start, reg_end) = genome_region;

    std::string fa_seq = reference[ref_id];  // the whole sequence of ``ref_id``
    std::string regstr = ref_id + "_" + ngslib::tostring(reg_start) + "-" + ngslib::tostring(reg_end);

    int bn = _args->input_bf.size() / _args->batchcount;
    if (_args->input_bf.size() % _args->batchcount > 0) 
        bn++;

    std::vector<std::string> batchfiles;

    ThreadPool thread_pool(_args->thread_num);
    std::vector<std::future<bool> > create_batchfile_processes;

    for (size_t i(0), j(1); i < _args->input_bf.size(); i+=_args->batchcount, ++j) {
        size_t x(i), y(i + _args->batchcount);
        std::vector<std::string> batch_align_files = ngslib::vector_slicing(_args->input_bf, x, y);

        // Set file path for batchfile
        std::string batchfile = cache_outdir + "/" + stem_bn + "." + regstr + "." + 
                                ngslib::tostring(j) + "_" + ngslib::tostring(bn) + ".bf.gz";
        batchfiles.push_back(batchfile);  // store the name of batchfile into a vector

        // Thread Pool
        create_batchfile_processes.emplace_back(
            thread_pool.enqueue(__create_single_batchfile, 
                                batch_align_files,  // 该值会变，只能拷贝，如果是引用，在多线程中会丢失变量
                                std::cref(fa_seq),  // 这个值在循环外，值不变，可以传引用
                                genome_region,
                                batchfile,
                                500)
        );
    }
    
    for (auto && p: create_batchfile_processes) {
        p.get();  // return the value of `__create_single_batchfile`
    }
    create_batchfile_processes.clear();

    return batchfiles;
}

void BaseTypeRunner::_variant_caller_process() {
    // 以区间为单位进行变异检测
    std::vector<std::string> batchfiles;
    for (size_t i(0); i < _calling_intervals.size(); ++i) {
        batchfiles = _create_batchfiles(_calling_intervals[i]);

std::cout << ngslib::join(batchfiles, "\n") << "\n\n";

    }
}

void BaseTypeRunner::run() {
    std::cout << "\n--- Running ---\n";
    _variant_caller_process();
    return;
}
