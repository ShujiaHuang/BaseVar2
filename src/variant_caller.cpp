/**
 * @file basetype_caller.cpp
 * @brief 
 * 
 * @author Shujia Huang
 * @date 2018-08-01
 * 
 */
#include "variant_caller.h"

const std::string BaseTypeRunner::usage() const {
    static const std::string BASETYPE_CALLER_USAGE = 
        "About: Call variants and estimate allele frequency by BaseVar.\n" 
        "Usage: basevar caller [options] <-f fasta> <--output-vcf out_fname> [--filename-has-samplename] [-L bam.list] in1.bam [in2.bam ...] ...\n\n" 

        "Required arguments:\n" 
        "  -f, --reference FILE         Input reference fasta file.\n"
        "  -o, --output    FILE         Output VCF file.\n\n"

        "Optional options:\n"
        "  -L, --align-file-list=FILE   BAM/CRAM files list, one file per row.\n"
        "  -r, --regions=REG[,...]      Skip positions which not in these regions. This parameter could be a list\n"
        "                               of comma deleimited genome regions(e.g.: chr:start-end).\n"
        "  -G, --pop-group=FILE         Calculating the allele frequency for specific population.\n\n"

        "  -m, --min-af=float           Setting prior precision of MAF and skip ineffective caller positions,\n"
        "                               a typical approach involves setting it to min(" + std::to_string(_args->min_af) + ", 100/x), where x \n"
        "                               represents the number of input BAM files min(" + std::to_string(_args->min_af) + ", 100/x). In most\n"
        "                               cases, users need not be overly concerned about this parameter, as it \n"
        "                               is generally handled automatically by the program.\n"
        "  -Q, --min-BQ INT             Skip bases with base quality < INT [" + std::to_string(_args->min_baseq) + "]\n"
        "  -q, --mapq=INT               Skip reads with mapping quality < INT [" + std::to_string(_args->min_mapq) + "]\n"
        "  -B, --batch-count=INT        INT base-pairs per batch. [" + std::to_string(_args->batchcount) + "]\n"  // window size
        "  -t, --thread=INT             Number of threads. [" + std::to_string(_args->thread_num) + "]\n\n"

        "  --filename-has-samplename    If the name of BAMfile/CRAMfile is something like 'SampleID.xxxx', set this\n"
        "                               argrument could save a lot of time during get the samples' id from BAMfiles.\n"
        "  -h, --help                   Show this help message and exit.\n\n"

        "Example usage:\n"
        "  [1]. basevar caller -f reference.fasta -o output.vcf.gz -Q 20 -q 30 -B 1000 --filename-has-samplename -L bam.list\n"
        "  [2]. basevar caller -f reference.fasta -o output.vcf.gz -Q 20 -q 30 -B 1000 --filename-has-samplename -L bam.list sample1.bam sample2.bam\n"
        "  [3]. basevar caller -f reference.fasta -o output.vcf.gz -Q 20 -q 30 -B 1000 --filename-has-samplename -L bam.list -r chr1\n"
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
        // Optional arguments to long style command line parameters require 'equals sign' (=). 
        // https://stackoverflow.com/questions/1052746/getopt-does-not-parse-optional-arguments-to-parameters
        {"reference",         required_argument, NULL, 'f'},
        {"output",            required_argument, NULL, 'o'},
        {"align-file-list",   optional_argument, NULL, 'L'},

        {"min-af",            optional_argument, NULL, 'm'},
        {"min-mapq",          optional_argument, NULL, 'q'},
        {"min-BQ",            optional_argument, NULL, 'Q'},
        {"batch-count",       optional_argument, NULL, 'B'},
        {"thread",            optional_argument, NULL, 't'},

        {"regions",           optional_argument, NULL, 'r'},
        {"pop-group",         optional_argument, NULL, 'G'},  // parameter for calculating allele frequency for specific population-group

        {"filename-has-samplename", no_argument, NULL, '1'},
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
        "   --filename-has-samplename" : "")   + " " + _args->input_bf[0] + " [... there are inputting " 
        << _args->input_bf.size() << " bamfile(s) in total]. \n"
        << std::endl;
    
    std::cout << "[INFO] Finish loading arguments and we have " << _args->input_bf.size()
              << " BAM/CRAM files for variants calling."        << std::endl;

    // Setting the resolution of AF
    _args->min_af = std::min(float(100)/_args->input_bf.size(), _args->min_af);
    reference = _args->reference;  // load fasta

    // get calling interval after loading the whole `reference` genome
    _make_calling_interval();

    // keep the order of '_samples_id' as the same as input 'aligne_files'
    _get_sample_id();  
    if (!_args->pop_group_file.empty()) 
        _get_popgroup_info();

    return;
}

void BaseTypeRunner::_get_sample_id() {
    time_t real_start_time = time(0);

    // Loading sample ID in BAM/CRMA files from RG tag.
    if (_args->filename_has_samplename)
        std::cout << "[INFO] Load samples' id from filename directly, "
                     "becuase you set --filename-has-samplename.\n";

    std::string samplename, filename;
    _samples_id.clear();
    for (size_t i(0); i < _args->input_bf.size(); ++i) {

        if ((i+1) % 1000 == 0)
            std::cout << "[INFO] loading "   << i+1 << "/" << _args->input_bf.size() 
                      << " alignment files." << std::endl;
        
        if (_args->filename_has_samplename) {
            filename = ngslib::remove_filename_extension(ngslib::basename(_args->input_bf[i]));
            size_t si = filename.find('.');  // Find the first index of '.'
            samplename = si > 0 && si != std::string::npos ? filename.substr(0, si) : filename;
        } else {
            // Get sampleID from BAM header, a bit time-consuming.
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

    // check samples duplication
    std::vector<std::string> duplicate_samples = ngslib::find_duplicates(_samples_id);
    if (!duplicate_samples.empty()) {
        std::cout << "[WARNING] Find " << duplicate_samples.size() << " duplicated samples within " 
                  << "the input bamfiles: " + ngslib::join(duplicate_samples, ",") + "\n\n";
    }

    // Time information
    time_t now = time(0);
    std::string ct(ctime(&now));
    ct.pop_back();  // rm the trailing '\n' put by `asctime`
    std::cout << "[INFO] " + ct + ". Done for loading all samples' id from alignment files, " 
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

    if (gr.size() == 2) {                // 'chr:start-end' or 'chr:start'
        std::vector<uint32_t> gs;
        ngslib::split(gr[1], gs, "-");   // get position coordinate
        reg_start = gs[0];
        reg_end = (gs.size() > 1) ? gs[1] : reference.seq_length(ref_id);
    } else {
        reg_start = 1;
        reg_end = reference.seq_length(ref_id);  // the whole chromosome
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
    std::vector<std::string> vcffiles;
    for (auto &gr: _calling_intervals) {
        clock_t cpu_start_time = clock();
        time_t real_start_time = time(0);

        //////////////////////////////////////////////////////////////
        // Calling variants in specific region with mulitple thread //
        //////////////////////////////////////////////////////////////
        std::string rgstr  = gr.chrom + "_" + std::to_string(gr.start) + "_" + std::to_string(gr.end);
        std::string prefix = cache_outdir + "/" + stem_bn + "." + rgstr;

        std::string sub_vcf_fn = prefix + ".vcf.gz";
        bool has_variants = _call_variants_in_region(gr, sub_vcf_fn);
        vcffiles.push_back(sub_vcf_fn);

        if (!has_variants) {
            std::cerr << "[INFO] No variants found in region: " << gr.to_string() << "\n";
        }

        // Time information
        time_t now = time(0);
        std::string ct(ctime(&now)); 
        ct.pop_back();
        std::cout << "[INFO] " + ct + ". Done for variants detection in " + gr.to_string() + ": "
                  << sub_vcf_fn + ", " << difftime(now, real_start_time) << " (CPU time: " 
                  << (double)(clock() - cpu_start_time) / CLOCKS_PER_SEC << ") seconds elapsed in total.\n" 
                  << std::endl;
    }

    // Merge all VCF subfiles
    std::vector<std::string> add_group_info, group_name;
    std::map<std::string, std::vector<size_t>>::iterator it = _groups_idx.begin();
    for(; it != _groups_idx.end(); ++it) {
        group_name.push_back(it->first);
        add_group_info.push_back("##INFO=<ID=" + it->first + "_AF,Number=A,Type=Float,Description="
                                 "\"Allele frequency in the " + it->first + " populations calculated "
                                 "base on LRT, in the range (0,1)\">");
    }

    // Merge VCF
    std::string header = vcf_header_define(_args->reference, add_group_info, _samples_id, _cmdline_string);
    merge_file_by_line(vcffiles, _args->output_vcf, header, true);

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

bool BaseTypeRunner::_call_variants_in_region(const ngslib::GenomeRegion gr, const std::string out_vcf_fn) {
    // Split region into sub-regions for multi-threading
    int bn = std::max(gr.size() / _args->batchcount, _args->thread_num);
    int STEP_LEN = (gr.size() % bn > 0) ? gr.size() / bn + 1 : gr.size() / bn;  // STEP_LEN 越长，所需内存就越大

    ThreadPool thread_pool(_args->thread_num);
    std::vector<std::future<bool>> call_variants_processes;
    std::vector<std::string> subvcfs;

    std::string fa_seq = reference[gr.chrom];  // use the whole sequence of chromosome
    for (uint32_t i(gr.start), j(1); i < gr.end + 1; i += STEP_LEN, ++j) {
        std::string tmp_vcf_fn = out_vcf_fn + "." + std::to_string(j) + "_" + std::to_string(bn);
        subvcfs.push_back(tmp_vcf_fn);

        uint32_t sub_gr_start = i;
        uint32_t sub_gr_end = std::min(sub_gr_start + STEP_LEN - 1, gr.end);
        ngslib::GenomeRegion sub_gr(gr.chrom, sub_gr_start, sub_gr_end);

        // Submit job to thread pool
        call_variants_processes.emplace_back(
            thread_pool.submit(std::bind(&BaseTypeRunner::_call_variant_unit, this,
                                         sub_gr,             // 局部变量，会变，必拷贝，不可传引用，否则线程执行时将丢失该值 
                                         std::cref(fa_seq),  // 外部变量，不变，传引用，省内存
                                         tmp_vcf_fn)));
    }

    // Wait for all threads to complete
    bool has_variants = false;
    for (auto & p: call_variants_processes) {
        // Run and make sure all processes are finished
        // 一般来说，只有当 valid() 返回 true 的时候才调用 get() 去获取结果，这也是 C++ 文档推荐的操作。
        if (p.valid()) {
            // get() 调用会改变其共享状态，不再可用，也就是说 get() 只能被调用一次，多次调用会触发异常。
            // 如果想要在多个线程中多次获取产出值需要使用 shared_future。
            bool x = p.get();
            if (!has_variants && x) has_variants = true;
        }
    }

    // Merge sub-VCFs
    std::string header = "## No need header here";
    merge_file_by_line(subvcfs, out_vcf_fn, header, true);

    return !has_variants;
}

/// Functions for calling variants outside of 'BaseTypeRunner' class 
// A unit for calling variants and let it run in a thread.
bool BaseTypeRunner::_call_variant_unit(const ngslib::GenomeRegion gr, const std::string &fa_seq, const std::string &out_vcf_fn) {    
    clock_t cpu_start_time = clock();
    time_t real_start_time = time(0);

    // Open output VCF file
    ngslib::BGZFile vcf_out(out_vcf_fn.c_str(), "w");

    // Process each sample for this position
    PosMapVector all_sample_posmaps;
    all_sample_posmaps.reserve(_args->input_bf.size());
    bool is_empty = _fetch_base_in_region(_args->input_bf, fa_seq, gr, all_sample_posmaps);  // 传引用，省内存，得数据
    if (is_empty) {
        vcf_out.close();
        return false;
    }

    std::vector<BaseType::BatchInfo> all_smps_bi_vector;  // Initialize vectors with the size of 'n_sample'
    int variant_count = 0, n_sample = all_sample_posmaps.size();
    all_smps_bi_vector.reserve(n_sample);
    for (uint32_t pos(gr.start); pos < gr.end + 1; ++pos) {

        PosMap::const_iterator smp_pos_it;  // specific sample position map
        for (size_t i = 0; i < n_sample; ++i) {
            smp_pos_it = all_sample_posmaps[i].find(pos);

            BaseType::BatchInfo smp_bi(gr.chrom, pos);  // create a BatchInfo for each sample
            if (smp_pos_it != all_sample_posmaps[i].end()) {
                if (smp_pos_it->second.ref_id != gr.chrom || smp_pos_it->second.ref_pos != pos)
                    throw std::runtime_error("[ERROR] reference id or position not match: " + gr.to_string());

                for (const auto &ab : smp_pos_it->second.align_bases) {
                    // REF may be single base for SNVs or a sub-seq for Indels
                    smp_bi.mapqs.push_back(ab.mapq);
                    smp_bi.map_strands.push_back(ab.map_strand);

                    smp_bi.ref_bases.push_back(ab.ref_base);     // It's upper case already
                    smp_bi.align_bases.push_back(ab.read_base);  // It's upper case already
                    smp_bi.align_base_quals.push_back(ab.base_qual);
                    smp_bi.base_pos_ranks.push_back(ab.rpr);
                }    
            } else {
                // Position not found in this sample
                smp_bi.mapqs.push_back(0);
                smp_bi.map_strands.push_back('*');

                smp_bi.ref_bases.push_back("N");
                smp_bi.align_bases.push_back("N");
                smp_bi.align_base_quals.push_back(BASE_Q0_CHAR);
                smp_bi.base_pos_ranks.push_back(0);
            }

            all_smps_bi_vector.push_back(smp_bi);
        }

        // check data
        if (all_smps_bi_vector.size() != n_sample) {
            throw std::runtime_error("[ERROR] The number of samples does not match. It should be " + std::to_string(n_sample) + 
                                    " but got " + std::to_string(all_smps_bi_vector.size()) + " samples in batchfiles for " + 
                                    gr.chrom + " " + std::to_string(pos) + ".\n");
        }

        // Call variants
        bool cc = _basevar_caller(all_smps_bi_vector, _groups_idx, _samples_id.size(), vcf_out);
        if (cc) variant_count++;
        all_smps_bi_vector.clear();  // clear for next position
    }
    
    vcf_out.close();

    // Time information
    time_t now = time(0);
    std::string ct(ctime(&now));
    ct.pop_back();
    std::cout << "[INFO] " + ct + ". Found " << variant_count << " variants and has been wrote in " + out_vcf_fn + ", "
              << difftime(now, real_start_time) 
              << " (CPU time: " << (double)(clock() - cpu_start_time) / CLOCKS_PER_SEC 
              << ") seconds elapsed." << std::endl;

    // Files will be automatically closed when destructors are called
    return variant_count > 0;
}

bool BaseTypeRunner::_fetch_base_in_region(const std::vector<std::string> &batch_align_files,  
                                           const std::string &fa_seq,  // must be the whole chromosome sequence
                                           const ngslib::GenomeRegion gr,
                                           PosMapVector &batchsamples_posinfomap_vector)  
{
    // only using here, In case of missing the overlap reads on side position, 100bp would be enough
    static const uint32_t REG_EXTEND_SIZE = 100;
    ngslib::GenomeRegion gr_extend(gr.chrom,                                                     // Reference id
                                   gr.start > REG_EXTEND_SIZE ? gr.start - REG_EXTEND_SIZE : 1,  // start
                                   gr.end + REG_EXTEND_SIZE);                                    // end
    std::string rg_extend_str = gr_extend.to_string(); // chr:start-end

    // Loop all alignment files
    bool is_empty = true;
    for(size_t i(0); i < batch_align_files.size(); ++i) {
        // 位点信息存入该变量, 且由于是按区间读取比对数据，key 值无需再包含 ref_id，因为已经不言自明。
        PosMap sample_posinfo_map;
        
        // Open BAM/CRAM file in reading mode (one sample, one bamfile)
        ngslib::Bam bf(batch_align_files[i], "r", _args->reference);
        if (bf.fetch(rg_extend_str)) { // Set 'bf' only fetch alignment reads in 'rg_extend_str'.
            hts_pos_t map_ref_start, map_ref_end;  // hts_pos_t is uint64_t
            std::vector<ngslib::BamRecord> sample_target_reads; 
            ngslib::BamRecord al;       // alignment read

            while (bf.next(al) >= 0) {  // -1 => hit the end of alignement file.
                if (al.mapq() < _args->min_mapq || al.is_duplicate() || al.is_qc_fail()) continue;
                if (al.is_secondary() || al.is_supplementary()) continue; // skip secondary and supplementary alignments
                
                map_ref_start = al.map_ref_start_pos() + 1;  // al.map_ref_start_pos() is 0-based, convert to 1-based
                map_ref_end   = al.map_ref_end_pos();        // al.map_ref_end_pos() is 1-based

                // Only fetch reads which in [gr.start, gr.end]
                if (map_ref_end < gr.start) continue;
                if (map_ref_start > gr.end) break;

                sample_target_reads.push_back(al);  // record the proper reads of sample
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
        throw std::runtime_error("[variant_caller.cpp::_fetch_base_in_region] 'pos_batchinfo_vector.size()' "
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

            if (map_ref_pos > gr.end) break;
            if (map_ref_pos < gr.start) continue;

            // 'BAM_XXX' are macros for CIGAR, which defined in 'htslib/sam.h'
            if (aligned_pairs[i].op == BAM_CMATCH ||  /* CIGAR: M */ 
                aligned_pairs[i].op == BAM_CEQUAL ||  /* CIGAR: = */
                aligned_pairs[i].op == BAM_CDIFF)     /* CIGAR: X */
            {
                // SNV. Only one base in 'ref_base' and 'read_base'.
                ab.ref_base  = aligned_pairs[i].ref_base[0];
                ab.read_base = aligned_pairs[i].read_base[0];
                ab.base_qual = aligned_pairs[i].read_qual[0];
            } else if (aligned_pairs[i].op == BAM_CINS) {  /* CIGAR: I */
                // Insertion. 'ref_base' is empty, 'read_base' is the inserted bases.
                if (!aligned_pairs[i].ref_base.empty()) {
                    std::cerr << al << "\n";
                    throw std::runtime_error("[ERROR] We got reference base in insertion region.");
                }

                // do not roll back position here
                ab.ref_base  = "";  // empty reference base
                ab.read_base = "+" + aligned_pairs[i].read_base;  // insertion bases

                // mean quality of the whole insertion sequence
                double total_qual = 0.0;
                for (char q : aligned_pairs[i].read_qual) {
                    total_qual += (q - 33);  // convert to int type
                }
                ab.base_qual = static_cast<char>(total_qual / aligned_pairs[i].read_qual.size() + 33); // convert to char type

            } else if (aligned_pairs[i].op == BAM_CDEL) {  /* CIGAR: D */
                // Deletion. 'read_base' is empty, 'ref_base' is the deleted bases.
                if (!aligned_pairs[i].read_base.empty()) {
                    std::cerr << al << "\n";
                    throw std::runtime_error("[ERROR] We got read bases in deletion region.");
                }

                // do not roll back position here
                ab.ref_base  = aligned_pairs[i].ref_base;
                ab.read_base = "-" + aligned_pairs[i].ref_base;  // aligned_pairs[i].ref_base 识别用的 (后面不直接用，要替换的)，同位点可以有多类型的 DEL 

                // mean quality of the whole deletion sequence
                ab.base_qual = static_cast<char>(al.mean_qqual() + 33); // convert to char type
            } else {
                // Skip in basevar for other kind of CIGAR symbals.
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

bool BaseTypeRunner::_basevar_caller(std::vector<BaseType::BatchInfo> &all_smps_bi_vector,
                                     const std::map<std::string, std::vector<size_t>> &group_smp_idx,
                                     size_t n_sample, ngslib::BGZFile &vcf_hd) 
{
    // Detect variant at this position by using all samples
    auto [bt, global_vi] = _basetype_caller_unit(all_smps_bi_vector);
    
    // check if there is a variant at this position
    bool is_variant = (global_vi.ale_bases.size() > 1) || 
                      (!global_vi.ale_bases.empty() && 
                        global_vi.ale_bases[0][0] !=
                        global_vi.ref_bases[0][0]);  // 'ale_bases' and 'ref_bases' are all upper case already.
    
    if (is_variant) { 
        // Collect and normalize allele information. Return the variant and allele information of this
        // position and replace the align_bases in samples_batchinfo_vector with the normalized bases.
        // AlleleInfo is a struct-type to collect and normalize allele information, which is used to 
        // output VCF record.
        AlleleInfo ai = collect_and_normalized_allele_info(global_vi, all_smps_bi_vector);

        std::map<std::string, BaseType> popgroup_bt;  // group_id => BaseType
        if (!group_smp_idx.empty()) { // group is not empty, call BaseType for each group
            std::vector<std::string> basecombination;
            basecombination.push_back(ai.ref); // reference base must be the first one
            basecombination.insert(basecombination.end(), ai.alts.begin(), ai.alts.end());

            std::map<std::string, std::vector<size_t>>::const_iterator sub_g_it = group_smp_idx.begin();
            for (; sub_g_it != group_smp_idx.end(); ++sub_g_it) {
                auto [g_bt, g_vi] = _basetype_caller_unit(all_smps_bi_vector,
                                                          sub_g_it->second, // group index
                                                          basecombination);
                popgroup_bt[sub_g_it->first] = g_bt;
            }
        }

        VCFRecord vcf_record = _vcfrecord_in_pos(all_smps_bi_vector, global_vi, popgroup_bt, ai);
        if (vcf_record.is_valid()) {
            // write to file
            vcf_hd << vcf_record.to_string() << "\n";
        } else {
            std::cerr << "[WARNING] Invalid VCF record at " << vcf_record.chrom << ":" << vcf_record.pos << " "
                      << vcf_record.ref << " " << ngslib::join(vcf_record.alt, ",")
                      << " ,skip writing this record.\n";
        }
    }

    // return true if there is a variant at this position
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
        if (smps_bi_v[i].ref_id != all_smps_bi.ref_id || 
            smps_bi_v[i].ref_pos != all_smps_bi.ref_pos) {
            throw std::runtime_error("Inconsistent reference position or ID");
        }
        
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
    if (basecombination.empty()) {
        // If basecombination is empty, use the default reference base
        bt.lrt();
    } else {
        // If basecombination is not empty, use the provided basecombination
        bt.lrt(basecombination);
    } 

    return {bt, get_pos_variant_info(bt, &all_smps_bi)};
}

VariantInfo BaseTypeRunner::get_pos_variant_info(const BaseType &bt, const BaseType::BatchInfo *smp_bi) {
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

    return vi;
}

VCFRecord BaseTypeRunner::_vcfrecord_in_pos(const std::vector<BaseType::BatchInfo> &samples_batchinfo_vector, // has been normalized with ref and alt bases by 'global_variant_info'
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
    for (const auto &smp_bi : samples_batchinfo_vector) { // Loop all samples in the batch file
        // For each sample, we need to collect the genotype information
        auto sa = process_sample_variant(vcf_record.ref, vcf_record.alt, smp_bi);

        // Update allele counts: include both ref and alt alleles
        for (const auto& alt : sa.sample_alts) {
            ai.allele_counts[alt]++;
            ai.total_alleles++;
        }

        // Format sample string in GT:PL:AD:DP, returned string is likely "0/1:100,50,0:30:80"
        vcf_record.samples.push_back(format_sample_string(sa));
    }

    // Construct the INFO field, like "AC=1;AN=2;AF=0.5;DP=14;DP4=10,4,0,0;"
    std::vector<int> ac;
    std::vector<double> af, caf;
    std::vector<double> fs, sor;
    std::vector<int> dp4({ai.strand_bias_info[ai.ref].fwd, ai.strand_bias_info[ai.ref].rev});
    for (const auto& alt : vcf_record.alt) { 
        // Only record the counts (AC) and frequencies (AF) of non-ref allele in INFO field
        ac.push_back(ai.allele_counts[alt]);
        af.push_back(ai.allele_freqs[alt]);  // allele frequency by LRT
        caf.push_back(ai.allele_counts[alt]/double(ai.total_alleles));

        dp4.push_back(ai.strand_bias_info[alt].fwd);  // alt_fwd
        dp4.push_back(ai.strand_bias_info[alt].rev);  // alt_rev
        fs.push_back(ai.strand_bias_info[alt].fs);    // strand bias
        sor.push_back(ai.strand_bias_info[alt].sor);  // symmetric odds ratio
    }

    std::vector<std::string> info = {
        "AF="  + ngslib::join(af, ","), 
        "CAF=" + ngslib::join(caf, ","),
        "AC="  + ngslib::join(ac, ","),
        "AN="  + std::to_string(ai.total_alleles),  // AN,total number of alleles
        "DP="  + std::to_string(ai.total_dp),       // DP, total depth of coverage
        "DP4=" + ngslib::join(dp4, ","),
        "FS="  + ngslib::join(fs, ","),
        "SOR=" + ngslib::join(sor, ",")
    };

    // Add group allele frequency information if available
    if (!group_bt.empty()) {
        std::vector<std::string> group_af_info;
        for (std::map<std::string, BaseType>::const_iterator it(group_bt.begin()); it != group_bt.end(); ++it) {

            std::vector<double> af;
            for (auto b : it->second.get_active_bases()) { // must have the same order with ai.ref + ai.alts
                if (b == ai.ref) continue;  // skip reference base

                af.push_back(it->second.get_lrt_af(b));  // 由于 active-base 不完全一样，所以存在与 ai.alts 顺序不一致的风险，要改进
            }

            if (!af.empty()) {
                group_af_info.push_back(it->first + "_AF=" + ngslib::join(af, ",")); // groupID_AF=xxx,xxx
            }
        }
        info.insert(info.end(), group_af_info.begin(), group_af_info.end());
    }

    vcf_record.info = ngslib::join(info, ";");

    return vcf_record;
}