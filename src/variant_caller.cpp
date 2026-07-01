/**
 * @file basetype_utils.cpp
 * @brief 
 * 
 * @author Shujia Huang
 * @date 2024-08-01
 * 
 */
#include "variant_caller.h"
#include "io/batchfile_binary.h"
#include <cstdio>      // snprintf
#include <stdexcept>
#include <string>

const std::string BaseTypeRunner::usage() const {
    static const std::string BASETYPE_CALLER_USAGE = 
        "About: Call variants and estimate allele frequency by BaseVar.\n" 
        "Usage: basevar caller [options] <-f Fasta> <--output-vcf output_file> [-L bam.list] in1.bam [in2.bam ...] ...\n\n" 

        "Required arguments:\n" 
        "  -f, --reference FILE         Input reference fasta file.\n"
        "  -o, --output    FILE         Output VCF file.\n\n"

        "Optional options:\n"
        "  -L, --align-file-list=FILE   BAM/CRAM files list, one file per row.\n"
        "  -r, --regions=REG[,...]      Skip positions which not in these regions. This parameter could be a list\n"
        "                               of comma deleimited genome regions(e.g.: chr:start-end).\n"
        "  -G, --pop-group=FILE         Calculating the allele frequency for specific population group.\n\n"

        "  -m, --min-af=float           Setting prior precision of MAF and skip ineffective caller positions,\n"
        "                               a typical approach involves setting it to min(" + std::to_string(_args->min_af) + ", 100/x), where x \n"
        "                               represents the number of input BAM files min(" + std::to_string(_args->min_af) + ", 100/x). In most\n"
        "                               cases, users need not be overly concerned about this parameter, as it \n"
        "                               is generally handled automatically by the program.\n"
        "  -Q, --min-BQ INT             Skip bases with base quality < INT [" + std::to_string(_args->min_baseq) + "]\n"
        "  -q, --mapq=INT               Skip reads with mapping quality < INT [" + std::to_string(_args->min_mapq) + "]\n"
        "  -B, --batch-count=INT        INT simples per batchfile. [" + std::to_string(_args->batchcount) + "]\n" 
        "  -t, --thread=INT             Number of threads. [" + std::to_string(_args->thread_num) + "]\n\n"

        "  --filename-has-samplename    If the prefix name of BAMfiles/CRAMfiles start with 'Sample ID', something like 'SampleID.bam', set this\n"
        "                               argrument could save a lot of time during get the sample id from BAMfiles.\n"
        "  --gt-mode=STRING             Genotype calling mode: 'posterior' (Bayesian posterior-based GT/GQ) or\n"
        "                               'legacy' (pure likelihood argmin GT, PL-gap GQ). [posterior]\n"
        "  --ref-bias=FLOAT             Reference bias coefficient (β) for genotype likelihood calculation.\n"
        "                               For het (REF,ALT), P(R|G) = (1-β)*P(R|REF) + β*P(R|ALT).\n"
        "                               β=0.5 means no bias (default); β<0.5 corrects for alignment reference bias.\n"
        "                               Typical empirical values: 0.45-0.48. [0.5]\n"
        "  --max-alleles=INT            Maximum number of active alleles allowed at a site.\n"
        "                               Sites exceeding this threshold will be skipped. [6]\n"
        "  --smart-rerun                Rerun process by checking batchfiles.\n"
        "  -h, --help                   Show this help message and exit.\n\n"

        "Example usage:\n"
        "  [1]. basevar caller -f reference.fasta -o output.vcf.gz -Q 20 -q 30 -B 500 --filename-has-samplename -L bam.list\n"
        "  [2]. basevar caller -f reference.fasta -o output.vcf.gz -Q 20 -q 30 -B 500 --filename-has-samplename -L bam.list sample1.bam sample2.bam\n"
        "  [3]. basevar caller -f reference.fasta -o output.vcf.gz -Q 20 -q 30 -B 500 --filename-has-samplename -L bam.list -r chr1 sample1.bam sample2.bam\n"
        ; 
        
    return BASETYPE_CALLER_USAGE;
}

BaseTypeRunner::BaseTypeRunner(int cmd_argc, char *cmd_argv[]) {
    // Firstly: Initialize a new BasTypeARGS and set default argument.
    _args = new BaseTypeARGS;

    if (cmd_argc < 2) {
        std::cout << usage() << "\n" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // Parsing the commandline options. 
    static const struct option BASETYPE_CMDLINE_LOPTS[] = {
        // All options that take a value use required_argument so that
        // both "--option value" and "--option=value" syntaxes work.
        // (optional_argument would require "--option=value" exclusively.)
        {"reference",         required_argument, NULL, 'f'},
        {"output",            required_argument, NULL, 'o'},
        {"align-file-list",   required_argument, NULL, 'L'},

        {"min-af",            required_argument, NULL, 'm'},
        {"min-mapq",          required_argument, NULL, 'q'},
        {"min-BQ",            required_argument, NULL, 'Q'},
        {"batch-count",       required_argument, NULL, 'B'},
        {"thread",            required_argument, NULL, 't'},

        {"regions",           required_argument, NULL, 'r'},
        {"pop-group",         required_argument, NULL, 'G'},  // parameter for calculating allele frequency for specific population-group

        {"filename-has-samplename", no_argument, NULL, '1'},
        {"smart-rerun",             no_argument, NULL, '2'},
        {"gt-mode",           required_argument, NULL, '3'},
        {"ref-bias",          required_argument, NULL, '4'},
        {"max-alleles",       required_argument, NULL, '5'},
        {"help",                    no_argument, NULL, 'h'},

        // must set this value
        {0, 0, 0, 0}
    };

    // Save the complete command line options in VCF header.
    // This code should be run before calling `getopt_long`
    _cmdline_string = "##basevar_caller_command=";
    for (size_t i = 0; i < cmd_argc; ++i) {
        _cmdline_string += (i > 0) ? " " + std::string(cmd_argv[i]) : std::string(cmd_argv[i]);
    }

    char c;
    std::vector<std::string> bv;
    while((c = getopt_long(cmd_argc, cmd_argv, "L:f:o:m:q:Q:B:t:r:G:h", BASETYPE_CMDLINE_LOPTS, nullptr)) >= 0) {
        // 字符流解决命令行参数转浮点等类型的问题
        std::stringstream ss(optarg ? optarg: "");  
        switch (c) {
            case 'f': _args->reference  = optarg;              break;  /* 临参 */
            case 'o': _args->output_vcf = optarg;              break;  // 恒参
            case 'L':
                bv = ngslib::get_firstcolumn_from_file(optarg); 
                _args->input_bf.insert(_args->input_bf.end(), 
                                       bv.begin(), bv.end());  break;  /* 临参 */

            case 'm': ss >> _args->min_af;                     break;  // 恒参
            case 'q': ss >> _args->min_mapq;                   break;  // 恒参
            case 'Q': ss >> _args->min_baseq;                  break;  // 恒参
            case 'B': ss >> _args->batchcount;                 break;  // 恒参
            case 't': ss >> _args->thread_num;                 break;  // 恒参

            case 'r': _args->regions                 = optarg; break;  /* 临参 */
            case 'G': _args->pop_group_file          = optarg; break;  // 恒参 
            case '1': _args->filename_has_samplename = true;   break;  // 恒参
            case '2': _args->smart_rerun             = true;   break;  // 恒参
            case '3': _args->posterior_gt = (std::string(optarg) != "legacy"); break;  // 恒参
            case '4': ss >> _args->ref_bias;                   break;  // 恒参
            case '5': ss >> _args->max_alleles;                break;  // 恒参
            case 'h': 
                std::cout << usage() << std::endl; 
                exit(EXIT_SUCCESS);

            default: 
                std::cerr << "Unknown argument: " << c << std::endl; 
                abort();
        }
    }

    // Collect BAM/CRAM files
    while (optind < cmd_argc) {
        _args->input_bf.push_back(cmd_argv[optind++]);
    }

    /* Make sure we set valid arguments */
    if (_args->input_bf.empty())
        throw std::invalid_argument("[ERROR] Missing required BAM/CRAM files.");
    if (_args->reference.empty())
        throw std::invalid_argument("[ERROR] Missing argument '-f/--reference'");
    if (_args->output_vcf.empty())
        throw std::invalid_argument("[ERROR] Missing argument '-o/--output-vcf'");
    
    if (_args->min_af < 0)
        throw std::invalid_argument("[ERROR] '-m/--min-af' argument must be > 0");
    if (_args->min_baseq < 0)
        throw std::invalid_argument("[ERROR] '-Q/--min-BQ' argument must be > 0");
    if (_args->min_baseq > 60)
        throw std::invalid_argument("[ERROR] '-Q/--min-BQ' argument must be <= 60");
    if (_args->min_mapq < 0)
        throw std::invalid_argument("[ERROR] '-q/--min-mapq' argument must be > 0");
    if (_args->batchcount <= 0)
        throw std::invalid_argument("[ERROR] '-B/--batch-count' argument must be > 0");
    if (_args->thread_num <= 0)
        throw std::invalid_argument("[ERROR] '-t/--thread' argument must be > 0");
    if (_args->ref_bias < 0.0 || _args->ref_bias > 1.0)
        throw std::invalid_argument("[ERROR] '--ref-bias' argument must be between 0.0 and 1.0");
    if (_args->max_alleles < 1 || _args->max_alleles > 10)
        throw std::invalid_argument("[ERROR] '--max-alleles' argument must be between 1 and 10");

    // recovering the absolute paths of output files
    _args->output_vcf = ngslib::abspath(_args->output_vcf);

    // Output the commandline options
    std::cout << 
        "[INFO] BaseVar arguments:\n"
        "basevar caller -f "+ _args->reference + " \\ \n"
        "   -Q " << _args->min_baseq          << " \\ \n"
        "   -q " << _args->min_mapq           << " \\ \n"
        "   -m " << _args->min_af             << " \\ \n"
        "   -B " << _args->batchcount         << " \\ \n"
        "   -t " << _args->thread_num         << " \\ \n"  + (_args->regions.empty() ? "" : 
        "   -r " + _args->regions              + " \\ \n") + (_args->pop_group_file.empty() ? "" : 
        "   -G " + _args->pop_group_file       + " \\ \n") +
        "   --output-vcf " + _args->output_vcf + " \\ \n"  + (_args->filename_has_samplename ? " \\ \n"
        "   --filename-has-samplename" : "")   + (_args->smart_rerun ? " \\ \n"
        "   --smart-rerun": "")                + " " + _args->input_bf[0] + " [... there are inputting " 
        << _args->input_bf.size() << " bamfiles in total]. \n"
        << std::endl;
    
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

    std::cout << "[INFO] Finish loading arguments and we have " << _args->input_bf.size()
              << " BAM/CRAM files for variants calling."        << std::endl;

    // Setting the resolution of AF
    _args->min_af = std::min(float(100)/_args->input_bf.size(), _args->min_af);
    reference = _args->reference;  // load fasta

    // get calling interval after loading the whole `reference` genome
    _make_calling_interval();

    // keep the order of '_samples_id' as the same as input 'aligne_files'
    _get_sample_id_from_bam();  
    if (!_args->pop_group_file.empty()) 
        _get_popgroup_info();

    return;
}

void BaseTypeRunner::_get_sample_id_from_bam() {
    time_t real_start_time = time(0);

    // Loading sample ID in BAM/CRMA files from RG tag.
    if (_args->filename_has_samplename)
        std::cout << "[INFO] BaseVar will load samples id from filename directly, becuase you set "
                     "--filename-has-samplename.\n";

    std::string samplename, filename;
    size_t si;
    _samples_id.clear();
    for (size_t i(0); i < _args->input_bf.size(); ++i) {

        if (!_args->filename_has_samplename && ((i+1) % 1000 == 0 || (i+1) == _args->input_bf.size())) { // print every 1000 files
            // Time information
            time_t now = time(0);
            std::string ct(ctime(&now));
            ct.pop_back();  // rm the trailing '\n' put by `asctime`
            std::cout << "[INFO] " + ct + ". loading " << i+1 << "/" << _args->input_bf.size() 
                      << " alignment files." << std::endl;
        }
        
        if (_args->filename_has_samplename) {
            filename = ngslib::remove_filename_extension(ngslib::basename(_args->input_bf[i]));
            si = filename.find('.');
            samplename = si > 0 && si != std::string::npos ? filename.substr(0, si) : filename;
        } else {
            // Get sampleID from BAM header, a bit time-consuming.
            ngslib::BamHeader bh(_args->input_bf[i], _args->reference);
            samplename = bh.get_sample_name();
        }

        if (!samplename.empty()) {
            _samples_id.push_back(samplename);
        } else {
            throw std::invalid_argument("[BaseTypeRunner::_load_sample_id_from_bam] " + 
                                        _args->input_bf[i] + " sample ID not found.\n");
        }
    }

    // check samples duplication
    std::vector<std::string> duplicate_samples = ngslib::find_duplicates(_samples_id);
    if (!duplicate_samples.empty()) {
        std::cout << "[WARNING] Find " << duplicate_samples.size() << " duplicated samples within " 
                  << "the input bamfiles: " + ngslib::join(duplicate_samples, ",") + "\n" 
                  << std::endl;
    }

    // Time information
    time_t now = time(0);
    std::string ct(ctime(&now));
    ct.pop_back();  // rm the trailing '\n' put by `asctime`
    std::cout << "[INFO] " + ct + ". Done for loading all " + std::to_string(_samples_id.size()) + 
                 " samples' id from alignment files, " 
              << difftime(now, real_start_time) << " seconds elapsed.\n" << std::endl;

    return;
}

void BaseTypeRunner::_make_calling_interval() {

    _calling_intervals.clear();  // clear the vector of calling intervals
    if (!_args->regions.empty()) {
        std::vector<std::string> rg_v;
        ngslib::split(_args->regions, rg_v, ",");

        for (size_t i(0); i < rg_v.size(); ++i) {
            if (rg_v[i].length() == 0) continue; // ignore empty string 
            _calling_intervals.push_back(_make_gregion_region(rg_v[i]));
        }
    } else {
        // Call the whole genome
        int n = reference.nseq();
        for (size_t i(0); i < n; ++i) {
            std::string ref_id = reference.iseq_name(i);
            _calling_intervals.push_back(ngslib::GenomeRegion(ref_id, 1, reference.seq_length(ref_id)));
        } 
    }
    
    return;
}

ngslib::GenomeRegion BaseTypeRunner::_make_gregion_region(const std::string &gregion) {
    // Genome Region, 1-based
    std::string ref_id; 
    uint32_t reg_start, reg_end;

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
                                    "-r/--regions " + gregion);
    }

    return ngslib::GenomeRegion(ref_id, reg_start, reg_end);  // 1-based
}

void BaseTypeRunner::_get_popgroup_info() {
    // group_id => [index in _samples_id of group_id]
    std::ifstream i_fn(_args->pop_group_file.c_str());
    if (!i_fn) {
        std::cerr << "[ERROR] Cannot open file: " + _args->pop_group_file << std::endl;
        exit(1);
    }

    std::map<std::string, std::string> sample2group;
    std::string skip, sn, gn;
    while (1) {
        // Only two columns: sample_id and group_id
        i_fn >> sn >> gn;
        if (i_fn.eof()) break;
        
        sample2group[sn] = gn;
        std::getline(i_fn, skip, '\n');  // skip the rest information of line.
    }
    i_fn.close();

    _groups_idx.clear();
    std::map<std::string, std::string>::iterator s2g_it;
    
    // follow the order of '_samples_id'
    if (sample2group.size() > 0) {
        for (size_t i(0); i < _samples_id.size(); ++i) {
            s2g_it = sample2group.find(_samples_id[i]);
            
            // ignore all the samples which not found
            if (s2g_it != sample2group.end()) {
                // record sample index of group groups
                // group -> index of _samples_id
                _groups_idx[sample2group[_samples_id[i]]].push_back(i);
            }
        }
    }

    return;
}

// Run the processes of calling variant and output files.
int BaseTypeRunner::run() {
    // Get filepath and stem name first.
    std::string _bname = ngslib::basename(_args->output_vcf);
    size_t si = _bname.find(".vcf");
    std::string stem_bn = (si > 0 && si != std::string::npos) ? _bname.substr(0, si) : _bname;

    std::string outdir = ngslib::dirname(ngslib::abspath(_args->output_vcf));
    std::string cache_outdir = outdir + "/cache_" + stem_bn;
    ngslib::safe_mkdir(cache_outdir);  // make cache directory for batchfiles

    std::cout << "---- Start calling variants ----\n" << std::endl;

    // 以区间为单位进行变异检测, 每个区间里调用多线程
    std::vector<std::string> batchfiles, vcffiles;

    // Build group-specific INFO headers (available before region loop since _groups_idx is already populated)
    std::vector<std::string> add_group_info, group_name;
    for(std::map<std::string, std::vector<size_t>>::iterator it = _groups_idx.begin(); it != _groups_idx.end(); ++it) {
        group_name.push_back(it->first);
        add_group_info.push_back("##INFO=<ID=DP_" + it->first + ",Number=1,Type=Integer,Description="
                                 "\"Total depth of all alleles in the " + it->first + " populations\">");
    }
    for(std::map<std::string, std::vector<size_t>>::iterator it = _groups_idx.begin(); it != _groups_idx.end(); ++it) {
        // For population group AF, it is still emphasized that it is based on EM estimation, 
        // because the prior of AF may be different for different populations. 
        // Therefore, the overall prior is not combined here.
        add_group_info.push_back("##INFO=<ID=AF_" + it->first + ",Number=A,Type=Float,Description="
                                 "\"Allele frequency in the " + it->first + " populations (EM-estimated), in the range (0,1)\">");
    }
    std::string header = vcf_header_define(_args->reference, add_group_info, _samples_id, _cmdline_string);

    for (auto &gr: _calling_intervals) {
        clock_t cpu_start_time = clock();
        time_t real_start_time = time(0);

        ///////////////////////////////////////////////////////////
        ///////// Create batchfiles with multiple-thread //////////
        ///////////////////////////////////////////////////////////
        std::string rgstr  = gr.chrom + "_" + std::to_string(gr.start) + "_" + std::to_string(gr.end);
        std::string prefix = cache_outdir + "/" + stem_bn + "." + rgstr;
        batchfiles = _create_batchfiles(gr, prefix);

        // Time information
        time_t now = time(0);
        std::string ct(ctime(&now)); 
        ct.pop_back();
        std::cout << "[INFO] " + ct + ". Done for creating all " << batchfiles.size() << " batchfiles in " 
                  << gr.to_string() + " and start to call variants, " << difftime(now, real_start_time) << " (CPU time: " 
                  << (double)(clock() - cpu_start_time) / CLOCKS_PER_SEC << ") seconds elapsed in total.\n" 
                  << std::endl;

        ///////////////////////////////////////////////////////////
        // Calling variants from batchfiles with mulitple thread //
        ///////////////////////////////////////////////////////////
        real_start_time = time(0);
        cpu_start_time  = clock();

        std::string sub_vcf_fn = prefix + ".vcf.gz";
        bool is_empty = _variants_discovery(batchfiles, gr, sub_vcf_fn, header);
        vcffiles.push_back(sub_vcf_fn);

        if (is_empty) {
            std::cerr << "[WARNING] No variants found in region: " << gr.to_string() << "\n";
        }

        // Time information
        now = time(0);
        ct  = ctime(&now); 
        ct.pop_back();
        std::cout << "[INFO] " + ct + ". Done for variants detection in " + gr.to_string() + ": "
                  << sub_vcf_fn + ", " << difftime(now, real_start_time) << " (CPU time: " 
                  << (double)(clock() - cpu_start_time) / CLOCKS_PER_SEC << ") seconds elapsed in total.\n" 
                  << std::endl;
        
        if (IS_DELETE_CACHE_BATCHFILE) {
            for (auto bf: batchfiles) {
                ngslib::safe_remove(bf);
                ngslib::safe_remove(bf+".bbi");
            }
        }

        batchfiles.clear(); 
    }

    // Merge or rename VCF output
    if (vcffiles.size() == 1) {
        // Single region: rename temp file directly (avoid unnecessary read/write IO)
        std::rename(vcffiles[0].c_str(), _args->output_vcf.c_str());
    } else {
        // Multiple regions: merge all subfiles (skip existing headers from temp files)
        merge_file_by_line(vcffiles, _args->output_vcf, header, true);
    }

    const tbx_conf_t bf_tbx_conf = {1, 1, 2, 0, '#', 0};  // {preset, seq col, beg col, end col, header-char, skip-line}
    if ((ngslib::suffix_name(_args->output_vcf) == ".gz") &&          // create index 
        tbx_index_build(_args->output_vcf.c_str(), 0, &bf_tbx_conf))  // file suffix is ".tbi"
        throw std::runtime_error("tbx_index_build failed: Is the file bgzip-compressed? "
                                 "Check this file: " + _args->output_vcf + "\n");

    if (IS_DELETE_CACHE_BATCHFILE) {
        ngslib::safe_remove(cache_outdir);
    }

    return 0; // SUCCESS
}

std::vector<std::string> BaseTypeRunner::_create_batchfiles(
    const ngslib::GenomeRegion gr, 
    const std::string bf_prefix) 
{
    ThreadPool thread_pool(_args->thread_num);  // set multiple-thread
    std::vector<std::future<bool>> create_batchfile_processes;

    int bn = _args->input_bf.size() / _args->batchcount;  // number of batchfiles
    if (_args->input_bf.size() % _args->batchcount > 0)
        bn++;

    std::vector<std::string> batchfiles;
    std::string fa_seq = reference[gr.chrom];  // use the whole sequence of ``ref_id`` for simply
    for (size_t i(0), j(1); i < _args->input_bf.size(); i+=_args->batchcount, ++j) {
        // set name of batchfile and must be compressed by BGZF.
        // Binary format: .bbf (data) + .bbi (index)
        std::string batchfile = bf_prefix + "." + std::to_string(j) + "_" + std::to_string(bn) + ".bbf";
        std::string batchfile_idx = bf_prefix + "." + std::to_string(j) + "_" + std::to_string(bn) + ".bbf.bbi";
        batchfiles.push_back(batchfile);   // Store the name of batchfile into a vector

        if (_args->smart_rerun && ngslib::is_readable(batchfile) && validate_binary_index(batchfile_idx)) {
            // do not create the existed batchfile again if set `--smart-rerun`
            // validate_binary_index() checks the footer to ensure .bbi was fully written
            std::cout << "[INFO] " + batchfile + " and " + batchfile_idx + " already exist and are valid, we don't have to "
                         "create them again, when we set `--smart-rerun`.\n";
            continue;
        }

        // slicing bamfiles for a batchfile.
        size_t x(i), y(i + _args->batchcount);
        std::vector<std::string> batch_align_files = ngslib::vector_slicing(_args->input_bf, x, y);
        std::vector<std::string> batch_sample_ids  = ngslib::vector_slicing(_samples_id, x, y);

        // make Thread Pool
        // 使用 std::bind 绑定成员函数和 this 指针
        create_batchfile_processes.emplace_back(
            thread_pool.submit(std::bind(&BaseTypeRunner::_create_a_batchfile, this,
                                         batch_align_files,  // 局部变量，会变，必拷贝，不可传引用，否则线程执行时将丢失该值
                                         batch_sample_ids,   // 局部变量，会变，必拷贝，不可传引用，否则线程执行时将丢失该值
                                         std::cref(fa_seq),  // 外部变量，不变，传引用，省内存
                                         gr,
                                         batchfile))
        );
    }
    
    for (auto && p: create_batchfile_processes) {
        // Run and make sure all processes are finished
        // 一般来说，只有当 valid() 返回 true 的时候才调用 get() 去获取结果，这是 C++ 文档推荐的操作。
        if (p.valid()) {
            // get() 调用会改变其共享状态，不再可用，也就是说 get() 只能被调用一次，多次调用会触发异常。
            // 如果想要在多个线程中多次获取产出值需要使用 shared_future。
            bool x = p.get(); // retrieve the return value of `__create_a_batchfile`
        }
    }

    create_batchfile_processes.clear();  // release the thread
    return batchfiles;  // done for batchfiles created and return
}

bool BaseTypeRunner::_create_a_batchfile(const std::vector<std::string>& batch_align_files,
                                         const std::vector<std::string>& batch_sample_ids,
                                         const std::string& fa_seq,
                                         const ngslib::GenomeRegion gr,
                                         const std::string output_batch_file)
{
    clock_t cpu_start_time = clock();
    time_t real_start_time = time(0);

    // This value affected the computing memory, could be set larger than 500000, set 20 for testting
    // 这个参数是为了限制存入 `batchsamples_posinfomap_vector` 的最大读取区间，从而控制内存消耗不要太大
    static const uint32_t STEP_REGION_LEN = 500000;  

    ngslib::BGZFile obf(output_batch_file.c_str(), "w"); // output file handle of output_batch_file

    // Write binary header (magic + version + sample IDs)
    write_binary_header(obf, batch_sample_ids);

    // Index entries: pos -> BGZF virtual offset (for .bbi)
    std::vector<BinaryIndexEntry> index_entries;
    index_entries.reserve(1000000);  // pre-allocate for ~1M positions per batchfile

    PosMapVector batchsamples_posinfomap_vector;
    batchsamples_posinfomap_vector.reserve(batch_align_files.size());  //  pre-set the capacity

    bool is_empty = true, has_data = false;
    uint32_t sub_reg_beg, sub_reg_end;  // 1-based region start and end
    for (uint32_t i(gr.start), j(0); i < gr.end + 1; i += STEP_REGION_LEN, ++j) {

        // Cut smaller regions to save computing memory.
        sub_reg_beg = i;
        sub_reg_end = sub_reg_beg + STEP_REGION_LEN - 1 > gr.end ? gr.end : sub_reg_beg + STEP_REGION_LEN - 1;
        ngslib::GenomeRegion sub_gr(gr.chrom, sub_reg_beg, sub_reg_end);
        is_empty = _fetch_base_in_region(batch_align_files, fa_seq, sub_gr,
                                         batchsamples_posinfomap_vector);  // 传引用，省内存，得数据

        if (!has_data && !is_empty) {
            has_data = true;
        }

        // Write binary records for positions with data (sparse index).
        // With hundreds of samples per batchfile, almost every position has data from at least one sample,
        // so we iterate positions sequentially and use early-break when finding data.
        for (uint32_t pos = sub_gr.start; pos <= sub_gr.end; ++pos) {
            bool has_data_at_pos = false;
            for (size_t s = 0; s < batchsamples_posinfomap_vector.size(); ++s) {
                auto it = batchsamples_posinfomap_vector[s].find(pos);
                if (it != batchsamples_posinfomap_vector[s].end() &&
                    it->second.ref_id == gr.chrom && 
                    it->second.ref_pos == pos &&
                    !it->second.align_bases.empty()) 
                {
                    has_data_at_pos = true;
                    break;  // early break: only need to know if ANY sample has data
                }
            }
            if (has_data_at_pos) {
                write_binary_record(obf, batchsamples_posinfomap_vector, sub_gr, pos, index_entries);
            }
        }
        batchsamples_posinfomap_vector.clear();  // 必须清空，为下个循环做准备
    }
    
    obf.close();

    // Build binary index (.bbi) - replaces tabix index
    std::string index_file = output_batch_file + ".bbi";  // e.g., file.bbf -> file.bbf.bbi
    write_binary_index(index_file, index_entries);

    // Time information
    time_t now = time(0);
    std::string ct(ctime(&now)); 
    ct.pop_back();  // rm the trailing '\n' put by `asctime`
    std::cout << "[INFO] " + ct + ". Done for creating batchfile " 
              << output_batch_file << ", " << difftime(now, real_start_time) 
              << " (CPU time: " << (double)(clock() - cpu_start_time) / CLOCKS_PER_SEC 
              << ") seconds elapsed." << std::endl;

    return has_data;
}

bool BaseTypeRunner::_fetch_base_in_region(const std::vector<std::string> &batch_align_files,  
                                           const std::string &fa_seq,  // must be the whole chromosome sequence
                                           const ngslib::GenomeRegion gr,
                                           PosMapVector &batchsamples_posinfomap_vector)  
{
    // only using here, In case of missing the overlap reads on side position, 200bp would be enough
    static const uint32_t REG_EXPEND_SIZE = 100;

    uint32_t extend_start = gr.start > REG_EXPEND_SIZE ? gr.start - REG_EXPEND_SIZE : 1;  // 1-based
    uint32_t extend_end   = gr.end + REG_EXPEND_SIZE;
    std::string extend_rgstr = gr.chrom + ":" + std::to_string(extend_start) + "-" + std::to_string(extend_end);

    // Loop all alignment files
    bool is_empty = true;
    for(size_t i(0); i < batch_align_files.size(); ++i) {
        // 位点信息存入该变量, 且由于是按区间读取比对数据，key 值无需再包含 ref_id，因为已经不言自明。
        PosMap sample_posinfo_map;
        
        // Open BAM/CRAM file in reading mode (one sample, one bamfile)
        ngslib::Bam bf(batch_align_files[i], "r", _args->reference);
        if (bf.fetch(extend_rgstr)) { // Set 'bf' only fetch alignment reads in 'exp_rgstr'.
            hts_pos_t map_ref_start, map_ref_end;  // hts_pos_t is uint64_t
            std::vector<ngslib::BamRecord> sample_target_reads; 
            ngslib::BamRecord al;       // alignment read

            while (bf.next(al) >= 0) {  // -1 => hit the end of alignement file.
                if (al.mapq() < _args->min_mapq || al.is_duplicate() || al.is_qc_fail()) {
                    continue;
                }
                if (al.is_secondary() || al.is_supplementary()) {
                    continue;  // skip secondary and supplementary alignments
                }

                map_ref_start = al.map_ref_start_pos() + 1;  // al.map_ref_start_pos() is 0-based, convert to 1-based
                map_ref_end   = al.map_ref_end_pos();        // al.map_ref_end_pos() is 1-based

                // Only fetch reads which in [reg_start, reg_end]
                if (gr.start > map_ref_end) continue;
                if (gr.end < map_ref_start) break;

                sample_target_reads.push_back(std::move(al));  // move BamRecord to avoid deep copy (P1-2)
            }

            sample_posinfo_map.clear();  // make sure it's empty
            if (sample_target_reads.size() > 0) {
                // get alignment information of [i] sample.
                _seek_position(sample_target_reads, fa_seq, gr, sample_posinfo_map);
            }
        }

        if (is_empty && !sample_posinfo_map.empty()) { 
            // at least one sample has data in this region
            is_empty = false; 
        }

        // Push it into 'batchsamples_posinfomap_vector' even if 'sample_posinfo_map' is empty, 
        // make sure 'batchsamples_posinfomap_vector' has the same size as `batch_align_files`
        batchsamples_posinfomap_vector.push_back(sample_posinfo_map);
    }

    if (batchsamples_posinfomap_vector.size() != batch_align_files.size())
        throw std::runtime_error("[basetype.cpp::__fetch_base_in_region] 'pos_batchinfo_vector.size()' "
                                 "should be the same as 'batch_align_files.size()'");

    return is_empty;  // no cover reads in 'genome_region' if empty.
}

void BaseTypeRunner::_seek_position(const std::vector<ngslib::BamRecord> &sample_map_reads,
                     const std::string &fa_seq,   // must be the whole chromosome sequence
                     const ngslib::GenomeRegion gr,
                     PosMap &sample_posinfo_map)
{
    if (!sample_posinfo_map.empty()){
        throw std::runtime_error("[basetype.cpp::__seek_position] 'sample_posinfo_map' must be empty.");
    }

    // A vector of: (cigar_op, read position, reference position, read base, read_qual, reference base)
    std::vector<ngslib::ReadAlignedPair> aligned_pairs;
    std::set<uint32_t> indel_pos_set;
    for(auto &al: sample_map_reads) {
        AlignBase ab;
        ab.map_strand = al.map_strand();  // '*', '-' or '+'
        ab.mapq = al.mapq();

        uint32_t map_ref_pos;
        aligned_pairs = al.get_aligned_pairs(fa_seq);
        for (size_t i(0); i < aligned_pairs.size(); ++i) {
            map_ref_pos = aligned_pairs[i].ref_pos + 1;  // ref_pos is 0-based, convert to 1-based;

            if (gr.end < map_ref_pos) break;
            if (gr.start > map_ref_pos) continue;

            // 'BAM_XXX' are macros for CIGAR, which defined in 'htslib/sam.h'
            if (aligned_pairs[i].op == BAM_CMATCH ||  /* CIGAR: M */ 
                aligned_pairs[i].op == BAM_CEQUAL ||  /* CIGAR: = */
                aligned_pairs[i].op == BAM_CDIFF)     /* CIGAR: X */
            {
                // SNV. ref_base and read_base are now char (P1-1 optimization).
                ab.ref_base  = std::string(1, aligned_pairs[i].ref_base);
                ab.read_base = std::string(1, aligned_pairs[i].read_base);
                ab.base_qual = aligned_pairs[i].read_qual;
            } else if (aligned_pairs[i].op == BAM_CINS) {  /* CIGAR: I */
                // Insertion. ref_base is '\0' for insertions.
                if (aligned_pairs[i].ref_base != '\0') {
                    std::cerr << al << "\n";
                    throw std::runtime_error("[ERROR] We got reference base in insertion region.");
                }

                // do not roll back position here
                ab.ref_base  = "";  // empty reference base
                // Build insertion bases: read_base (first) + multi_base (rest)
                ab.read_base = std::string("+") + aligned_pairs[i].read_base + aligned_pairs[i].multi_base;

                // Use upstream anchor base quality (D4 fix)
                if (i > 0) {
                    ab.base_qual = aligned_pairs[i-1].read_qual;
                } else {
                    // Fallback: mean quality of insertion sequence
                    double total_qual = (aligned_pairs[i].read_qual - 33);
                    for (char q : aligned_pairs[i].multi_base) {
                        total_qual += (q - 33);
                    }
                    size_t ins_len = 1 + aligned_pairs[i].multi_base.size();
                    ab.base_qual = static_cast<char>(total_qual / ins_len + 33);
                }

            } else if (aligned_pairs[i].op == BAM_CDEL) {  /* CIGAR: D */
                // Deletion. read_base is '\0' for deletions.
                if (aligned_pairs[i].read_base != '\0') {
                    std::cerr << al << "\n";
                    throw std::runtime_error("[ERROR] We got read bases in deletion region.");
                }

                // do not roll back position here
                // Build deleted ref bases: ref_base (first) + multi_base (rest)
                ab.ref_base  = std::string(1, aligned_pairs[i].ref_base) + aligned_pairs[i].multi_base;
                ab.read_base = "-" + ab.ref_base;  // 识别用的 (后面不直接用，要替换的)，同位点可以有多类型的 DEL 

                // Use upstream anchor base quality (D4 fix)
                if (i > 0) {
                    ab.base_qual = aligned_pairs[i-1].read_qual;
                } else {
                    // Fallback: mean quality of the whole read
                    ab.base_qual = static_cast<char>(al.mean_qqual() + 33);
                }
            } else {
                // Skip in basevar for other kind of CIGAR symbols.
                continue;
            }

            // qpos is 0-based, conver to 1-based to set the rank of base on read.
            ab.rpr = aligned_pairs[i].qpos + 1;
            if (ab.base_qual < _args->min_baseq + 33) continue; // Skip low quality bases, 33 is the offset of base QUAL

            // Need to convert the ref and read_base to upper case in case of the base is a lower case.
            std::transform(ab.ref_base.begin(), ab.ref_base.end(), ab.ref_base.begin(), ::toupper);
            std::transform(ab.read_base.begin(), ab.read_base.end(), ab.read_base.begin(), ::toupper);

            // 以 map_ref_pos 为 key，将所有的 read_bases 信息存入 map 中，多个突变共享同个 ref_pos
            if (sample_posinfo_map.find(map_ref_pos) == sample_posinfo_map.end()) {
                AlignInfo pos_align_info(gr.chrom, map_ref_pos);
                sample_posinfo_map.insert({map_ref_pos, pos_align_info});
            }

            // collected all the align base infromations of the mapping read and store it into the map
            sample_posinfo_map[map_ref_pos].align_bases.push_back(ab);

            // Record the reference position for Indels
            if (ab.read_base[0] == '-' || ab.read_base[0] == '+') {
                indel_pos_set.insert(map_ref_pos);
            }
        }
    }

    // 对于 Indel 位点，按照 Indel 优先原则对位点的 align info 做修改
    for (auto pos: indel_pos_set) {
        uint32_t leftmost_pos = pos > 1 ? pos - 1 : pos; // roll back one position to the leftmost of Indels break point.
        AlignInfo raw_align_info = sample_posinfo_map[pos];
        sample_posinfo_map.erase(pos);

        AlignInfo indel_info(raw_align_info.ref_id, leftmost_pos);  // record Indels on the leftmost position
        AlignInfo non_indel_info(raw_align_info.ref_id, pos);       // raw position for non-indel variants
        for (auto ab: raw_align_info.align_bases) {
            if (ab.read_base[0] == '-' || ab.read_base[0] == '+') {
                ab.ref_base = fa_seq[leftmost_pos - 1] + ab.ref_base;  // add one leftmost ref-base
                std::transform(ab.ref_base.begin(), ab.ref_base.end(), ab.ref_base.begin(), ::toupper);
                indel_info.align_bases.push_back(ab); // 同位点可以有多类型的 DEL/INS
            } else {
                non_indel_info.align_bases.push_back(ab);
            }
        }

        if (indel_info.align_bases.empty()) {
            throw std::runtime_error("[ERROR] Must get at least one Indel here.");
        } else {
            // 如果 leftmost_pos 的位点信息将被 Indel 信息覆盖
            sample_posinfo_map[leftmost_pos] = indel_info;
        }

        // If there is still non-Indel information, we still need to keep the information
        // at the original position.
        if (non_indel_info.align_bases.size() > 0) {
            sample_posinfo_map[pos] = non_indel_info;
        }
    }

    return;
}

bool BaseTypeRunner::_variants_discovery(const std::vector<std::string> &batchfiles, 
                                         const ngslib::GenomeRegion gr,
                                         const std::string out_vcf_fn,
                                         const std::string &vcf_header) 
{
    // Split the region into multiple sub-regions for multiple-thread calling
    // Each sub-region is about 'STEP_LEN' bp length
    // Then merge all the sub-vcf files into one file.
    if (batchfiles.empty()) {
        throw std::invalid_argument("[ERROR] No batchfiles for calling variants in " + gr.to_string() + "\n");
    }

    // Pre-load all binary indexes ONCE (shared across all threads to avoid N_threads × N_batchfiles copies)
    // After Fix 1 (skip zero-depth positions), each index is small (~MB), so this is very efficient.
    std::vector<std::vector<BinaryIndexEntry>> shared_indexes;
    shared_indexes.reserve(batchfiles.size());
    for (const auto &file : batchfiles) {
        shared_indexes.push_back(load_binary_index(file + ".bbi"));
    }

    // decide how many sub-regions we should split
    int bn = _args->thread_num;
    uint32_t STEP_LEN = (gr.end - gr.start + 1) / bn;
    if ((gr.end - gr.start + 1) % bn) STEP_LEN++;
    
    // prepare multiple-thread
    ThreadPool thread_pool(_args->thread_num);  
    std::vector<std::future<bool>> call_variants_processes;
    std::vector<std::string> subvcfs;
    for (uint32_t i(gr.start), j(1); i < gr.end + 1; i += STEP_LEN, ++j) {
        std::string tmp_vcf_fn = out_vcf_fn + "." + std::to_string(j) + "_" + std::to_string(bn);
        subvcfs.push_back(tmp_vcf_fn);

        uint32_t sub_gr_start = i; // 1-based
        uint32_t sub_gr_end = (sub_gr_start+STEP_LEN-1) > gr.end ? gr.end : sub_gr_start+STEP_LEN-1;
        std::string gr_str = gr.chrom + ":" + std::to_string(sub_gr_start) + "-" + std::to_string(sub_gr_end);

        // Performance multi-thread here.
        call_variants_processes.emplace_back(
            thread_pool.submit(std::bind(&BaseTypeRunner::_variant_calling_unit, this,
                                         std::cref(batchfiles), 
                                         std::cref(_samples_id),  // 作为类变量，这个传参是多余的，之后再改吧，下同
                                         std::cref(_groups_idx),
                                         std::cref(shared_indexes),  // shared indexes (const ref, no copy)
                                         gr_str,         // 局部变量必须拷贝，会变 
                                         tmp_vcf_fn)));  // 局部变量必须拷贝，会变 
    }

    // Run and make sure all processes could be finished.
    bool is_empty = true;
    for (auto & p: call_variants_processes) {
        if (p.valid()) {
            bool x = p.get();
            if (is_empty && x) is_empty = false;
        }
    }
    call_variants_processes.clear();  // release the thread
  
    merge_file_by_line(subvcfs, out_vcf_fn, vcf_header, true);

    return is_empty;
}

// A unit for calling variants and let it run in a thread.
bool BaseTypeRunner::_variant_calling_unit(const std::vector<std::string> &batchfiles,  // total batchfiles in the region
                                           const std::vector<std::string> &sample_ids,  // total samples
                                           const std::map<std::string, std::vector<size_t>> &group_smp_idx,
                                           const std::vector<std::vector<BinaryIndexEntry>> &shared_indexes,  // pre-loaded indexes shared across threads
                                           const std::string region,    // genome region format like samtools
                                           const std::string tmp_vcf_fn) 
{
    clock_t cpu_start_time = clock();
    time_t real_start_time = time(0);

    /*************** Preparing for reading data **************/
    // Get and check the sample id from batchfiles header.
    std::vector<std::string> bf_smp_ids;
    for (const auto &file : batchfiles) {
        auto ids = read_binary_sample_ids(file);
        bf_smp_ids.insert(bf_smp_ids.end(), ids.begin(), ids.end());
    }
    if (ngslib::join(bf_smp_ids, ",") != ngslib::join(sample_ids, ","))
        throw std::runtime_error("[BUG] The order of sample ids in batchfiles must be the same as "
                                 "input bamfiles.\n" 
                                 "Sample ids in input bamfiles: " + ngslib::join(sample_ids, ",") + "\n"
                                 "Sample ids in batchfiles    : " + ngslib::join(bf_smp_ids, ",") + "\n");
    
    // Parse region string (e.g., "chr1:1000-2000")
    std::string query_chrom;
    uint32_t query_start = 0, query_end = 0;
    {
        size_t colon_pos = region.find(':');
        size_t dash_pos  = region.find('-');
        query_chrom = region.substr(0, colon_pos);
        query_start = std::stoul(region.substr(colon_pos + 1, dash_pos - colon_pos - 1));
        query_end   = std::stoul(region.substr(dash_pos + 1));
    }

    // Open .bbf files and use shared indexes (loaded once in _variants_discovery)
    struct BinaryReader {
        std::unique_ptr<ngslib::BGZFile> bf;
        const std::vector<BinaryIndexEntry> *index;  // pointer to shared index (no copy)
        uint16_t                         sample_count;
        size_t                           current_idx;
    };

    std::vector<BinaryReader> readers;
    readers.reserve(batchfiles.size());

    for (size_t fi = 0; fi < batchfiles.size(); ++fi) {
        const auto &file = batchfiles[fi];
        BinaryReader reader;
        reader.index = &shared_indexes[fi];  // point to shared index (no copy, no load)

        // Open file for reading
        reader.bf = std::make_unique<ngslib::BGZFile>(file.c_str(), "r");

        // Read header: magic(4) + version(2) + sample_count(2) + ids_len(4) + ids_string
        uint32_t magic;
        uint16_t version;
        reader.bf->read_raw(&magic, sizeof(magic));
        reader.bf->read_raw(&version, sizeof(version));
        reader.bf->read_raw(&reader.sample_count, sizeof(reader.sample_count));

        uint32_t ids_len;
        reader.bf->read_raw(&ids_len, sizeof(ids_len));
        if (ids_len > 0) {
            std::string dummy(ids_len, '\0');
            reader.bf->read_raw(&dummy[0], ids_len);
        }

        // Binary search for query_start in shared index
        auto it = std::lower_bound(reader.index->begin(), reader.index->end(), query_start,
            [](const BinaryIndexEntry &e, uint32_t pos) { return e.pos < pos; });
        reader.current_idx = (it != reader.index->end()) ? (it - reader.index->begin()) : reader.index->size();

        readers.push_back(std::move(reader));
    }

    // Seek each reader's file to its first in-range index entry.
    // With sparse indexes, different readers may start at different positions,
    // so we do NOT force all readers to the same starting position.
    // The min_pos loop below handles the different starting positions naturally.
    for (auto &r : readers) {
        if (r.current_idx < r.index->size()) {
            r.bf->seek_virtual((*r.index)[r.current_idx].virtual_offset);
        }
    }

    // Output VCF file
    ngslib::BGZFile VCF_OUT(tmp_vcf_fn.c_str(), "w");

    // Reusable buffer for position-level data (avoid per-position reallocation)
    std::vector<BaseType::BatchInfo> all_smps_bi;
    all_smps_bi.reserve(sample_ids.size());

    // Pre-constructed template for empty samples (depth=0) to avoid repeated construction
    // Note: ref_id/ref_pos are set per-position before use, so we don't set them here
    BaseType::BatchInfo empty_bi_template;
    empty_bi_template.ref_bases.push_back("N");
    empty_bi_template.align_bases.push_back("N");
    empty_bi_template.align_base_quals.push_back('!');
    empty_bi_template.base_pos_ranks.push_back(0);
    empty_bi_template.mapqs.push_back(0);
    empty_bi_template.map_strands.push_back('*');

    int n = 0, variant_count = 0;
    bool all_done = false;
    while (!all_done) {
        // Check if all readers are still in range
        all_done = true;
        for (auto &r : readers) {
            if (r.current_idx < r.index->size() && (*r.index)[r.current_idx].pos <= query_end) {
                all_done = false;
                break;
            }
        }
        if (all_done) break;

        // With sparse indexes (zero-depth positions skipped), different batchfiles
        // may have data at different positions. Find the minimum position across all
        // active readers, then only read from batchfiles at that position.
        uint32_t min_pos = UINT32_MAX;
        for (auto &r : readers) {
            if (r.current_idx < r.index->size() && (*r.index)[r.current_idx].pos <= query_end) {
                min_pos = std::min(min_pos, (*r.index)[r.current_idx].pos);
            }
        }
        if (min_pos == UINT32_MAX) break;  // no reader has data in range

        // Read one record from each batchfile at min_pos; fill empty for others
        all_smps_bi.clear();
        std::string ref_id;
        uint32_t ref_pos = min_pos;
        int total_depth = 0;

        for (auto &r : readers) {
            if (r.current_idx < r.index->size() && (*r.index)[r.current_idx].pos == min_pos) {
                // This batchfile has data at min_pos — read it
                std::vector<BaseType::BatchInfo> smp_bi;
                std::string rid;
                uint32_t rpos;
                int depth;
                if (!read_binary_record(*r.bf, r.sample_count, rid, rpos, depth, smp_bi)) {
                    throw std::runtime_error("[ERROR] Unexpected EOF reading binary batchfile");
                }

                if (ref_id.empty()) {
                    ref_id = rid;
                }
                total_depth += depth;

                // Set ref_id/ref_pos once for all samples from this reader
                for (auto &smp : smp_bi) {
                    smp.ref_id = rid;
                    smp.ref_pos = rpos;
                }

                all_smps_bi.insert(all_smps_bi.end(),
                                   std::make_move_iterator(smp_bi.begin()),
                                   std::make_move_iterator(smp_bi.end()));

                ++r.current_idx;
            } else {
                // This batchfile has no data at min_pos — fill empty BatchInfo
                for (uint16_t s = 0; s < r.sample_count; ++s) {
                    all_smps_bi.push_back(empty_bi_template);
                }
            }
        }

        // Fill ref_id/ref_pos for samples from batchfiles that had no data at this position.
        // These samples were copied from empty_bi_template with default ref_id="" and ref_pos=0,
        // but downstream _basetype_caller_unit uses smps_bi_v[0].ref_id/ref_pos for VCF output,
        // so ALL samples must have consistent values.
        for (auto &smp : all_smps_bi) {
            if (smp.ref_id.empty()) smp.ref_id = ref_id;
            if (smp.ref_pos == 0)   smp.ref_pos = ref_pos;
        }

        if (!all_smps_bi.empty() && total_depth > 0) {
            bool cc = _basevar_caller_binary(all_smps_bi, ref_id, ref_pos, total_depth,
                                              group_smp_idx, sample_ids.size(), VCF_OUT);
            if (cc) ++variant_count;
        }

        if (++n % 100000 == 0) {
            // Time information
            time_t now = time(0);
            std::string ct(ctime(&now));
            ct.pop_back();
            std::cout << "[INFO] " + ct + ". Processed " << n << " positions, " << variant_count << " variants found. "
                      << difftime(now, real_start_time) << " (CPU time: " << (double)(clock() - cpu_start_time) / CLOCKS_PER_SEC 
                      << ") seconds elapsed.\n";

        }
    }
    VCF_OUT.close();

    // Time information
    time_t now = time(0);
    std::string ct(ctime(&now));
    ct.pop_back();
    std::cout << "[INFO] " + ct + ". Found " << variant_count << " variants and has been wrote in " + tmp_vcf_fn + ", "
              << difftime(now, real_start_time) 
              << " (CPU time: " << (double)(clock() - cpu_start_time) / CLOCKS_PER_SEC 
              << ") seconds elapsed." << std::endl;

    // Files will be automatically closed when destructors are called
    return variant_count > 0;
}



/**
 * @brief Binary batchfile variant caller. Same logic as _basevar_caller() but receives
 *        pre-parsed BatchInfo data directly from binary reader (no text parsing).
 */
bool BaseTypeRunner::_basevar_caller_binary(std::vector<BaseType::BatchInfo> &all_smps_bi_vector,
                                            const std::string &ref_id,
                                            uint32_t ref_pos,
                                            int total_depth,
                                            const std::map<std::string, std::vector<size_t>> &group_smp_idx,
                                            size_t n_sample,
                                            ngslib::BGZFile &vcf_hd)
{
    if (total_depth == 0) return false;

    if (all_smps_bi_vector.size() != n_sample) {
        throw std::runtime_error(
            "[ERROR] The number of samples does not match. It should be " + std::to_string(n_sample) + 
            " but got " + std::to_string(all_smps_bi_vector.size()) + " samples in batchfiles for " + 
            ref_id + " " + std::to_string(ref_pos) + ".\n"
        );
    }
    
    // Detect variant at this position by using all samples
    auto [bt, global_vi] = _basetype_caller_unit(all_smps_bi_vector);
    
    // check if there is a variant at this position
    bool is_variant = (global_vi.ale_bases.size() > 1) || 
                      (!global_vi.ale_bases.empty() && 
                        global_vi.ale_bases[0][0] !=
                        global_vi.ref_bases[0][0]);
    if (is_variant) { 
        AlleleInfo ai = collect_and_normalized_allele_info(global_vi, all_smps_bi_vector);

        std::map<std::string, BaseType> popgroup_bt;
        if (!group_smp_idx.empty()) {
            std::vector<std::string> basecombination;
            basecombination.push_back(ai.ref);
            basecombination.insert(basecombination.end(), ai.alts.begin(), ai.alts.end());

            std::map<std::string, std::vector<size_t>>::const_iterator sub_g_it = group_smp_idx.begin();
            for (; sub_g_it != group_smp_idx.end(); ++sub_g_it) {
                auto [g_bt, g_vi] = _basetype_caller_unit(all_smps_bi_vector,
                                                          sub_g_it->second,
                                                          basecombination);
                popgroup_bt[sub_g_it->first] = g_bt;
            }
        }

        VCFRecord vcf_record = _vcfrecord_in_pos(all_smps_bi_vector, global_vi, popgroup_bt, ai);
        if (vcf_record.is_valid() && vcf_record.ref != "N") {
            vcf_hd << vcf_record.to_string() << "\n";
        } else {
            std::cerr << "[WARNING] Invalid VCF record at " << ref_id << ":" << ref_pos << " "
                      << vcf_record.ref << " " << ngslib::join(vcf_record.alt, ",")
                      << " ,skip writing this record.\n";
        }
    }

    return is_variant;
}

std::pair<BaseType, VariantInfo> BaseTypeRunner::_basetype_caller_unit(
    const std::vector<BaseType::BatchInfo> &smps_bi_v, 
    const std::vector<size_t> group_idx, 
    const std::vector<std::string> basecombination)
{
    // 如果 group_idx 为空,使用所有样本的索引
    std::vector<size_t> indices;
    if (group_idx.empty()) {
        indices.reserve(smps_bi_v.size());
        for (size_t i = 0; i < smps_bi_v.size(); ++i) {
            indices.push_back(i);
        }
    } else {
        indices = group_idx;
    }
    
    BaseType::BatchInfo all_smps_bi(smps_bi_v[0].ref_id, smps_bi_v[0].ref_pos);
    for (auto i : indices) {
        for (size_t j(0); j < smps_bi_v[i].ref_bases.size(); ++j) {
            all_smps_bi.ref_bases.push_back(smps_bi_v[i].ref_bases[j]);     // It's upper case already in batchfile.
            all_smps_bi.align_bases.push_back(smps_bi_v[i].align_bases[j]); // It's upper case already in batchfile.
            all_smps_bi.align_base_quals.push_back(smps_bi_v[i].align_base_quals[j]);
            all_smps_bi.mapqs.push_back(smps_bi_v[i].mapqs[j]);
            all_smps_bi.map_strands.push_back(smps_bi_v[i].map_strands[j]);
            all_smps_bi.base_pos_ranks.push_back(smps_bi_v[i].base_pos_ranks[j]);
        }
    }

    BaseType bt(&all_smps_bi, _args->min_af);
    bt.set_max_alleles(_args->max_alleles);  // D5 fix: configurable max alleles
    if (basecombination.empty()) {
        // If basecombination is empty, use the default reference base
        bt.lrt();
    } else {
        // If basecombination is not empty, use the provided basecombination
        bt.lrt(basecombination);
    } 

    return {bt, _get_pos_variant_info(bt, &all_smps_bi)};
}

VariantInfo BaseTypeRunner::_get_pos_variant_info(const BaseType &bt, const BaseType::BatchInfo *smp_bi) {
    VariantInfo vi(bt.get_ref_id(), bt.get_ref_pos(), bt.get_total_depth(), bt.get_var_qual()+0.5);
    int major_allele_depth = 0;
    vi.major_allele_idx    = 0;

    for (size_t i(0); i < bt.get_active_bases().size(); ++i) {
        std::string b = bt.get_active_bases()[i];
        std::string ref_base = bt.get_bases2ref().at(b);

        vi.ref_bases.push_back(ref_base); // could only be ref-base
        vi.ale_bases.push_back(b);        // could be ref or non-ref alleles
        vi.depths.push_back(bt.get_base_depth(b));
        vi.freqs.push_back(bt.get_lrt_af(b));

        if (major_allele_depth < bt.get_base_depth(b)) {
            vi.major_allele_idx = i;
            major_allele_depth  = bt.get_base_depth(b);
        }
    }

    // calculate the Strand Bias
    for (size_t i(0); i < vi.ale_bases.size(); ++i) {
        vi.strand_bias.push_back(strand_bias(vi.ale_bases[vi.major_allele_idx],
                                             vi.ale_bases[i],
                                             smp_bi->align_bases,
                                             smp_bi->map_strands));
    }

    // D1 fix: Recompute QUAL using overall LRT statistic (full model vs REF-only)
    // Uses the full (pre-step-down) model logL for correct degrees of freedom.
    // We take the MAX of global LRT QUAL and original stepwise LRT QUAL, to ensure
    // we don't lose evidence in ultra-low-coverage scenarios where the global LRT
    // may be weaker than the stepwise LRT.
    if (!vi.ale_bases.empty() && !vi.ref_bases.empty()) {
        std::string ref_base_str = vi.ref_bases[0];
        // Get the reference base character (single base for SNPs)
        char ref_char = ref_base_str.empty() ? 'N' : std::toupper(ref_base_str[0]);
        
        // Compute REF-only log-likelihood (natural log scale, matching EM's log marginal likelihood)
        double ref_only_logL = 0.0;
        for (size_t j = 0; j < smp_bi->align_bases.size(); ++j) {
            char base = std::toupper(smp_bi->align_bases[j][0]);
            if (base == 'N') continue;
            
            int Q = static_cast<int>(smp_bi->align_base_quals[j]) - 33;
            double e = (Q < 0) ? 1.0 : std::pow(10.0, -Q / 10.0);
            double p = (base == ref_char) ? (1.0 - e) : (e / 3.0);
            if (p > 0) {
                ref_only_logL += std::log(p);  // natural log to match EM scale
            }
        }
        
        // Get full model's logL (pre-step-down, D1 df fix)
        // This ensures df = N-1 matches the chi-squared statistic correctly.
        double full_logL = bt.get_full_model_logL();
        
        // LRT statistic: Lambda = 2 * (ln_L_alt - ln_L_null)
        // Both are in natural log scale, so Lambda = 2 * delta
        double delta_logL = full_logL - ref_only_logL;
        if (delta_logL > 0) {
            double lambda = 2.0 * delta_logL;
            int df = std::max(1, static_cast<int>(vi.ale_bases.size()) - 1);
            double chi_prob = chi2_test(lambda, df);
            if (std::isnan(chi_prob)) chi_prob = 1.0;
            double new_qual = (chi_prob > 0) ? -10.0 * std::log10(chi_prob) : 10000.0;
            if (new_qual < 0) new_qual = 0.0;
            int global_qual = static_cast<int>(new_qual + 0.5);
            // Take the max of global LRT QUAL and original stepwise LRT QUAL
            vi.qual = std::max(vi.qual, global_qual);
        }
        // When delta_logL <= 0, keep the original stepwise LRT QUAL (vi.qual unchanged)
    }

    return vi;
}

VCFRecord BaseTypeRunner::_vcfrecord_in_pos(
    const std::vector<BaseType::BatchInfo> &samples_batchinfo_vector, // has been normalized with ref and alt bases by 'global_variant_info'
    const VariantInfo &global_variant_info,
    const std::map<std::string, BaseType> &group_bt, 
    AlleleInfo &ai)
{
    // Return empty record if no variants
    if (samples_batchinfo_vector.empty()) return VCFRecord();

    VCFRecord vcf_record;
    vcf_record.chrom = global_variant_info.ref_id;  // all samples should have the same ref_id
    vcf_record.pos   = global_variant_info.ref_pos; // all samples should have the same ref_pos
    vcf_record.qual  = global_variant_info.qual;    // variant quality score, use the one from VariantInfo
    vcf_record.ref   = ai.ref;                      // normalized reference base (upper case) at this position
    vcf_record.alt   = ai.alts;                     // normalized alternative bases (upper case) at this position

    // Add FORMAT field
    vcf_record.format = "GT:GQ:PL:AD:DP";  // Genotype, Genotype Quality, Phred-scaled likelihoods of genotypes, Active allele depth, total coverage reads

    // =====================================================================
    // Per-sample likelihood-based genotype calling (argmin PL, no population prior)
    // Collect PL for posterior update; also collect allele counts from most-likely GT
    // =====================================================================
    std::vector<std::vector<int>> sample_pls;  // For HWE test (D7)
    sample_pls.reserve(samples_batchinfo_vector.size());
    size_t n_alt_alleles = vcf_record.alt.size();
    std::vector<int> ac_obs(n_alt_alleles, 0);  // ALT counts from argmin(PL) genotypes (distinct from VCF GT)
    int an_obs = 0;                             // Total alleles from non-missing argmin(PL) calls (distinct from VCF GT)
    for (const auto &smp_bi : samples_batchinfo_vector) {
        auto sa = process_sample_variant(vcf_record.ref, vcf_record.alt, smp_bi, -1.0, _args->ref_bias);

        // Update reads-based allele counts (used when no population prior is applied)
        for (const auto& alt : sa.sample_alts) {
            ai.allele_counts[alt]++;
            ai.total_alleles++;
        }

        // Collect allele counts from argmin(PL) GT for AC_obs/AN_obs (skip zero-coverage samples)
        if (sa.gtcode.size() == 2 && !sa.sample_alts.empty()) {
            for (int allele_code : sa.gtcode) {
                an_obs++;
                if (allele_code > 0 && static_cast<size_t>(allele_code) <= n_alt_alleles) {
                    ac_obs[allele_code - 1]++;
                }
            }
        }

        sample_pls.push_back(sa.PL);  // Collect PL for HWE
        vcf_record.samples.push_back(format_sample_string(sa));
    }

    // =====================================================================
    // Posterior genotype update: incorporate population AF as HW prior
    // AC/AN/AF = expected from posterior probabilities; AC_obs/AN_obs/AF_obs = from most-likely GT
    // =====================================================================
    if (_args->posterior_gt) {
        // Build per-allele frequency vector: [freq_REF, freq_ALT1, freq_ALT2, ...]
        size_t n_alleles = 1 + n_alt_alleles;
        std::vector<double> allele_freqs(n_alleles);
        double ref_freq = 1.0;
        for (size_t a = 0; a < n_alt_alleles; ++a) {
            double af = ai.allele_freqs[vcf_record.alt[a]];
            allele_freqs[a + 1] = af;
            ref_freq -= af;
        }
        static const double AF_MIN = 1e-6;
        allele_freqs[0] = std::max(AF_MIN, ref_freq);

        std::vector<GenotypePosterior> sample_posteriors;
        sample_posteriors.reserve(samples_batchinfo_vector.size());

        std::vector<int> ac_gt(n_alt_alleles, 0);  // ALT counts from posterior GT (matches VCF GT column)
        int an_gt = 0;                             // Total alleles from non-missing posterior GT calls

        for (size_t i = 0; i < samples_batchinfo_vector.size(); ++i) {
            // Compute genotype posterior with per-allele frequencies (multi-allelic capable)
            GenotypePosterior gp = compute_genotype_posterior(
                sample_pls[i],
                allele_freqs
            );

            // Rebuild sample annotation with posterior GT/GQ
            VCFSampleAnnotation sa;
            sa.PL = sample_pls[i];
            sa.GQ = gp.gq;
            sa.posterior = gp.posteriors;
            sa.dosage = gp.dosage;
            sa.per_allele_dosage = gp.per_allele_dosage;

            // Convert PL index to genotype codes
            auto [a1, a2] = pl_index_to_genotype(gp.best_gt_idx, n_alleles);
            sa.gtcode = {a1, a2};

            // Store allele depths (from first pass sa is gone, recompute quickly)
            // Note: allele depths don't change with posterior, re-derive from batch info
            {
                auto sa_orig = process_sample_variant(
                    vcf_record.ref, vcf_record.alt,
                    samples_batchinfo_vector[i], -1.0, _args->ref_bias);
                sa.allele_depths = sa_orig.allele_depths;
                sa.sample_alts = sa_orig.sample_alts;

                // Count alleles from posterior GT for AC_GT/AN_GT (skip zero-coverage samples)
                if (!sa_orig.sample_alts.empty()) {
                    for (int allele_code : sa.gtcode) {
                        an_gt++;
                        if (allele_code > 0 && static_cast<size_t>(allele_code) <= n_alt_alleles) {
                            ac_gt[allele_code - 1]++;
                        }
                    }
                }
            }

            sample_posteriors.push_back(gp);

            // Update sample string with posterior-based GT/GQ
            vcf_record.samples[i] = format_sample_string(sa);
        }

        // Compute dosage-based AC (per-ALT)
        auto [dosage_ac, dosage_an] = compute_dosage_ac(sample_posteriors);
        ai.dosage_counts = std::move(dosage_ac);
        ai.dosage_total_alleles = dosage_an;
        ai.gt_counts = std::move(ac_gt);
        ai.gt_total_alleles = an_gt;
    }

    // =====================================================================
    // Construct INFO field
    // =====================================================================
    std::vector<int> ac;
    std::vector<double> af;
    std::vector<double> fs, sor;
    std::vector<int> dp4({ai.strand_bias_info[ai.ref].fwd, ai.strand_bias_info[ai.ref].rev});

    bool use_dosage = (ai.dosage_total_alleles > 0);  // using posterior genotype expectations

    for (size_t a = 0; a < n_alt_alleles; ++a) {
        const auto& alt = vcf_record.alt[a];

        if (use_dosage) {
            // AC/AN/AF from posterior genotype expectations
            double ac_d = (a < ai.dosage_counts.size()) ? ai.dosage_counts[a] : 0.0;
            ac.push_back(static_cast<int>(std::round(ac_d)));
            af.push_back(ac_d / ai.dosage_total_alleles);
        } else {
            // Without population prior: AC/AN from reads-based counts, AF from EM-estimated frequency
            ac.push_back(ai.allele_counts[alt]);
            af.push_back(ai.allele_freqs[alt]);
        }

        dp4.push_back(ai.strand_bias_info[alt].fwd);  // alt_fwd
        dp4.push_back(ai.strand_bias_info[alt].rev);  // alt_rev
        fs.push_back(ai.strand_bias_info[alt].fs);    // strand bias
        sor.push_back(ai.strand_bias_info[alt].sor);  // symmetric odds ratio
    }

    int an = use_dosage ? ai.dosage_total_alleles : ai.total_alleles;

    std::vector<std::string> info = {
        "AF="  + ngslib::join(af, ","),
        "AC="  + ngslib::join(ac, ","),
        "AN="  + std::to_string(an),                  // AN, total number of alleles
        "DP="  + std::to_string(ai.total_dp),         // DP, total depth of coverage
        "DP4=" + ngslib::join(dp4, ","),
        "FS="  + ngslib::join(fs, ","),
        "SOR=" + ngslib::join(sor, ",")
    };

    // Add observed allele count fields (AC_obs, AN_obs, AF_obs) when using posterior expectations
    if (use_dosage) {
        info.push_back("AC_obs=" + ngslib::join(ac_obs, ","));
        info.push_back("AN_obs=" + std::to_string(an_obs));
        std::vector<double> af_obs(n_alt_alleles);
        for (size_t a = 0; a < n_alt_alleles; ++a) {
            af_obs[a] = (an_obs > 0) ? static_cast<double>(ac_obs[a]) / an_obs : 0.0;
        }
        info.push_back("AF_obs=" + ngslib::join(af_obs, ","));

        // Add GT-based allele count fields (AC_GT, AN_GT, AF_GT) — directly from VCF GT column
        info.push_back("AC_GT=" + ngslib::join(ai.gt_counts, ","));
        info.push_back("AN_GT=" + std::to_string(ai.gt_total_alleles));
        std::vector<double> af_gt(n_alt_alleles);
        for (size_t a = 0; a < n_alt_alleles; ++a) {
            af_gt[a] = (ai.gt_total_alleles > 0) ?
                static_cast<double>(a < ai.gt_counts.size() ? ai.gt_counts[a] : 0) / ai.gt_total_alleles : 0.0;
        }
        info.push_back("AF_GT=" + ngslib::join(af_gt, ","));
    }

    // =====================================================================
    // HWE test (D7): dosage-based, suitable for low-coverage data
    // =====================================================================
    {
        size_t n_alleles = vcf_record.alt.size() + 1;  // REF + ALTs
        // Convert PLs to likelihoods (uniform prior => posterior = likelihood)
        std::vector<std::vector<double>> sample_gt_probs;
        sample_gt_probs.reserve(sample_pls.size());
        for (const auto& pl : sample_pls) {
            sample_gt_probs.push_back(pl_to_likelihoods(pl));
        }
        double hwe_pval = hwe_dosage_test(sample_gt_probs, n_alleles);
        char buf[64];
        snprintf(buf, sizeof(buf), "%.6g", hwe_pval);
        info.push_back("HWE=" + std::string(buf));
    }

    // =====================================================================
    // RankSum annotations (D8): BaseQ, MQ, ReadPos
    // =====================================================================
    if (vcf_record.alt.size() == 1) {  // Only for bi-allelic sites
        const std::string& alt_allele = vcf_record.alt[0];
        std::vector<double> ref_bq, alt_bq, ref_mq, alt_mq, ref_rpr, alt_rpr;
        
        for (const auto& smp_bi : samples_batchinfo_vector) {
            for (size_t j = 0; j < smp_bi.align_bases.size(); ++j) {
                std::string read_base_upper = smp_bi.align_bases[j];
                std::transform(read_base_upper.begin(), read_base_upper.end(), read_base_upper.begin(), ::toupper);
                
                bool is_ref = (read_base_upper == vcf_record.ref);
                bool is_alt = (read_base_upper == alt_allele);
                
                if (is_ref || is_alt) {
                    double bq = static_cast<double>(smp_bi.align_base_quals[j]) - 33.0;
                    double mq = static_cast<double>(smp_bi.mapqs[j]);
                    double rpr = static_cast<double>(smp_bi.base_pos_ranks[j]);
                    
                    if (is_ref) {
                        ref_bq.push_back(bq); ref_mq.push_back(mq); ref_rpr.push_back(rpr);
                    } else {
                        alt_bq.push_back(bq); alt_mq.push_back(mq); alt_rpr.push_back(rpr);
                    }
                }
            }
        }
        
        // Compute Z-scores (need at least 1 read in each group)
        if (!ref_bq.empty() && !alt_bq.empty()) {
            double bq_z = wilcoxon_ranksum_zscore(ref_bq, alt_bq);
            double mq_z = wilcoxon_ranksum_zscore(ref_mq, alt_mq);
            double rpr_z = wilcoxon_ranksum_zscore(ref_rpr, alt_rpr);
            
            char buf[64];
            snprintf(buf, sizeof(buf), "%.3f", bq_z);
            info.push_back("BaseQRankSum=" + std::string(buf));
            snprintf(buf, sizeof(buf), "%.3f", mq_z);
            info.push_back("MQRankSum=" + std::string(buf));
            snprintf(buf, sizeof(buf), "%.3f", rpr_z);
            info.push_back("ReadPosRankSum=" + std::string(buf));
        }
    }

    // Add group allele frequency information if available
    if (!group_bt.empty()) {
        std::vector<std::string> group_af_info, group_dp_info;
        for (std::map<std::string, BaseType>::const_iterator it(group_bt.begin()); it != group_bt.end(); ++it) {
            group_dp_info.push_back("DP_" + it->first + "=" + std::to_string(it->second.get_total_depth())); // DP_group=xxx

            // 原本所存在群组 alt 的与全局 active_bases 不完全一样，所以会导致 ai.alts 
            // 顺序不一致的风险已经完全通过如下代码修复了。
            // D9 fix: Iterate over vcf_record.alt (VCF ALT order) instead of
            // get_active_bases() to ensure AF order matches the ALT field.
            std::vector<double> af;
            const auto& group_active = it->second.get_active_bases();
            for (const auto& alt : vcf_record.alt) {
                bool found = false;
                for (const auto& b : group_active) {
                    if (b == alt) {
                        af.push_back(it->second.get_lrt_af(b));
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    af.push_back(0.0);  // Alt not in group's active bases
                }
            }

            if (!af.empty()) {
                // 注意：这里关于不同群组的 AF，我还是保持了 EM 估计的 AF，没有进行任何的先验加权。
                // 原因是因为不同群组的 AF 的先验可能不同，因此不能进行简单的加权。
                group_af_info.push_back("AF_" + it->first + "=" + ngslib::join(af, ",")); // AF_group=xxx,xxx
            }
        }
        info.insert(info.end(), group_dp_info.begin(), group_dp_info.end());
        info.insert(info.end(), group_af_info.begin(), group_af_info.end());
    }
    vcf_record.info = ngslib::join(info, ";");

    return vcf_record;
}
