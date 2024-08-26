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
#include <ctime>
#include <algorithm>  // std::min

#include <htslib/bgzf.h>

#include "basetype.h"
#include "external/threadpool.h"


void BaseTypeRunner::set_arguments(int cmd_argc, char *cmd_argv[]) {

    if (cmd_argc < 2) {
        std::cout << usage() << "\n" << std::endl;
        exit(1);
    }

    if (_args) {
        throw std::runtime_error("[basetype.cpp::BaseTypeRunner:args] 'args' must be "
                                 "a NULL pointer before it can be assigned a value.");
    }
    // Inital a new BasTypeARGS and set defaut argument.
    _args = new BaseTypeARGS;
    
    // Parsing the commandline options. 
    char c;
    while((c = getopt_long(cmd_argc, cmd_argv, "I:L:R:m:q:B:t:r:G:h", BASETYPE_CMDLINE_LOPTS, NULL)) >= 0) {
        // 字符流解决命令行参数转浮点等类型的问题
        std::stringstream ss(optarg ? optarg: "");  
        switch (c) {
            case 'I': _args->input_bf.push_back(optarg);         break;  // 恒参 (一直用)
            case 'L': _args->in_bamfilelist = optarg;            break;  /* 临参 */
            case 'R': _args->reference      = optarg;            break;  /* 临参 */

            case 'm': ss >> _args->min_af;                       break;  // 恒参
            case 'q': ss >> _args->mapq;                         break;  // 恒参
            case 'B': ss >> _args->batchcount;                   break;  // 恒参
            case 't': ss >> _args->thread_num;                   break;  // 恒参

            case 'r': _args->regions        = optarg;            break;  /* 临参 */
            case 'G': _args->pop_group_file = optarg;            break;  // 
            case '1': _args->output_vcf     = optarg;            break;  // 恒参
            case '2': _args->output_cvg     = optarg;            break;  // 恒参

            case '3': _args->filename_has_samplename = true;     break;  // 恒参
            case '4': _args->smart_rerun             = true;     break;  // 恒参
            case 'h': 
                std::cout << usage() << std::endl;
                exit(1);

            default: 
                std::cerr << "Unknown argument: " << c << std::endl; 
                exit(1);
        }
    }
    /* Make sure we set valid arguments */
    if (_args->input_bf.empty() && _args->in_bamfilelist.empty())
        throw std::invalid_argument("[ERROR] Missing argument '-I/--input' or '-L/--align-file-list'");
    if (_args->reference.empty())
        throw std::invalid_argument("[ERROR] Missing argument '-R/--reference'");

    if (_args->output_vcf.empty())
        throw std::invalid_argument("[ERROR] Missing argument '--output-vcf'");
    if (_args->output_cvg.empty())
        throw std::invalid_argument("[ERROR] Missing argument '--output-cvg'");
    
    if (_args->min_af <= 0)
        throw std::invalid_argument("[ERROR] '-m/--min-af' argument must be > 0");
    if (_args->mapq <= 0)
        throw std::invalid_argument("[ERROR] '-q/--mapq' argument must be > 0");
    if (_args->batchcount <= 0)
        throw std::invalid_argument("[ERROR] '-B/--batch-count' argument must be > 0");
    if (_args->thread_num <= 0)
        throw std::invalid_argument("[ERROR] '-t/--thread' argument must be > 0");

    // recovering the absolute paths of output files
    _args->output_vcf = ngslib::abspath(_args->output_vcf);
    _args->output_cvg = ngslib::abspath(_args->output_cvg);

    // Output the commandline options
    std::cout << 
        "[INFO] BaseVar commandline and arguments:\n"
        "basevar basetype -R " + _args->reference + " \\ \n" + (_args->input_bf.empty() ? "" : 
        "   -I " + ngslib::join(_args->input_bf, " -I ") + " \\ \n") + (_args->in_bamfilelist.empty() ? "" : 
        "   -L " + _args->in_bamfilelist       + " \\ \n") <<
        "   -R " + _args->reference            + " \\ \n"
        "   -q " << _args->mapq               << " \\ \n"
        "   -m " << _args->min_af             << " \\ \n"
        "   -B " << _args->batchcount         << " \\ \n"
        "   -t " << _args->thread_num         << " \\ \n"  << (_args->regions.empty() ? "" : 
        "   -r " + _args->regions              + " \\ \n") << (_args->pop_group_file.empty() ? "" : 
        "   -p " + _args->pop_group_file       + " \\ \n") <<
        "   --output-vcf " + _args->output_vcf + " \\ \n"
        "   --output-vcg " + _args->output_cvg << (_args->filename_has_samplename ? " \\ \n"
        "   --filename-has-samplename" : "")   << (_args->smart_rerun ? " \\ \n"
        "   --smart-rerun": "")                << "\n" << std::endl;
    
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
    _args->min_af = std::min(float(100)/_args->input_bf.size(), _args->min_af);

    // load fasta
    reference = _args->reference;

    _get_calling_interval();
    print_calling_interval();
    _get_sample_id_from_bam();  // get sample id from input aligne_files and record in `_samples_id`
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
    std::cout << "[INFO] Done for loading all samples' id from alignment files.\n\n";

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

ngslib::GenomeRegionTuple BaseTypeRunner::_make_gregion_tuple(std::string gregion) {

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

    return std::make_tuple(ref_id, reg_start, reg_end);  // 1-based
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

/**
 * @brief Create a batch of temp files for variant discovery (could be deleted when the jobs done).
 * 
 * @param genome_region 
 * @return std::vector<std::string> 
 * 
 */
std::vector<std::string> BaseTypeRunner::_create_batchfiles(ngslib::GenomeRegionTuple genome_region, 
                                                            std::string prefix) 
{
    std::string ref_id; uint32_t reg_start, reg_end;
    std::tie(ref_id, reg_start, reg_end) = genome_region;

    std::string regstr = ref_id + "_" + ngslib::tostring(reg_start) + "-" + ngslib::tostring(reg_end);
    std::string bf_prefix = prefix + "." + regstr;

    int bn = _args->input_bf.size() / _args->batchcount;
    if (_args->input_bf.size() % _args->batchcount > 0)
        bn++;

    ThreadPool thread_pool(_args->thread_num);  // set multiple-thread
    std::vector<std::future<bool> > create_batchfile_processes;

    std::string fa_seq = reference[ref_id];  // use the whole sequence of ``ref_id`` for simply
    std::vector<std::string> batchfiles;
    for (size_t i(0), j(1); i < _args->input_bf.size(); i+=_args->batchcount, ++j) {
        // set name for batchfile
        std::string batchfile = bf_prefix + "." + ngslib::tostring(j) + "_" + ngslib::tostring(bn) + ".bf.gz";
        batchfiles.push_back(batchfile);   // Store the name of batchfile into a vector

        if (_args->smart_rerun && ngslib::is_readable(batchfile)) {
            // do not have to create the exists batchfiles again if set `--smart-rerun`
            std::cout << "[INFO] " + batchfile + " already exists, we don't have to "
                         "create it again, when we set `--smart-rerun`.\n";
            continue;
        }

        size_t x(i), y(i + _args->batchcount);
        std::vector<std::string> batch_align_files = ngslib::vector_slicing(_args->input_bf, x, y);
        std::vector<std::string> batch_sample_ids  = ngslib::vector_slicing(_samples_id, x, y);

        // make Thread Pool
        create_batchfile_processes.emplace_back(
            thread_pool.enqueue(__create_a_batchfile, 
                                batch_align_files,  // 循环内变量，值会变，只能拷贝，不可传引用，否则多线程执行时将丢失该值
                                batch_sample_ids,   // 循环内变量，值会变，只能拷贝，不可传引用，否则多线程执行时将丢失该值
                                std::cref(fa_seq),  // 循环外变量，值不变，可传引用，省内存
                                genome_region,
                                _args->mapq,
                                batchfile));
    }
    
    for (auto && p: create_batchfile_processes) {
        // Run and make sure all processes are finished

        // 一般来说，只有当 valid() 返回 true 的时候才调用 get() 去获取结果，这也是 C++ 文档推荐的操作。
        if (p.valid()) {
            // get() 调用会改变其共享状态，不再可用，也就是说 get() 只能被调用一次，多次调用会触发异常。
            // 如果想要在多个线程中多次获取产出值需要使用 shared_future。
            bool x = p.get(); // retrieve the return value of `__create_a_batchfile`
        }
    }
    create_batchfile_processes.clear();  // release the thread

    return batchfiles;  // done for batchfiles created and return
}

void BaseTypeRunner::_variant_caller_process() {

    // Get filepath stem name
    std::string _bname = ngslib::basename(_args->output_vcf);
    size_t si = _bname.find_first_of(".vcf");
    std::string stem_bn = si > 0 && si != std::string::npos ? _bname.substr(0, si) : _bname;

    std::string outdir = ngslib::dirname(_args->output_vcf);
    std::string cache_outdir = outdir + "/cache_" + stem_bn;
    ngslib::safe_mkdir(cache_outdir);  // make cache directory for batchfiles

    if (_args->smart_rerun) {
        // Remove and rollback `thread_num` last modification files. 
        // Must do these before calling '_create_batchfiles', and DO NOT do these process in loop!
        for (size_t i(0); i < _args->thread_num; ++i)
            ngslib::safe_remove(ngslib::get_last_modification_file(cache_outdir));
    }

    // 构造 batchfile 的临时目录和输出文件名前缀
    std::string prefix = cache_outdir + "/" + stem_bn;

    // 以区间为单位进行变异检测, 每个区间里再调用多线程
    std::vector<std::string> batchfiles;
    for (size_t i(0); i < _calling_intervals.size(); ++i) {
        batchfiles = _create_batchfiles(_calling_intervals[i], prefix);
        std::cout << "[INFO] Done for creating all " << i << " - " << batchfiles.size() 
                  << " batchfiles and start to call variants.\n";
        
        // calling variants from batchfiles
        /* add codes here */

    }

    return;
}

void BaseTypeRunner::run() {
    // Run the process of calling variant and output files.
    _variant_caller_process();
    return;
}

/// Calling variant functions outside of class 'BaseTypeRunner' 
// Create temp batch file for variant discovery
void __write_record_to_batchfile(PosMapVector &batchsamples_posinfomap_vector,
                                 const std::string &fa_seq,
                                 const ngslib::GenomeRegionTuple genome_region, 
                                 BGZF *obf) 
{
    const static char BASE_Q0_ASCII = '!';  // '!' ascii code is 33

    std::string ref_id; uint32_t reg_start, reg_end;
    std::tie(ref_id, reg_start, reg_end) = genome_region;  // 1-based

    // Output columns: 
    // [CHROM, POS, REF, Depth(CoveredSample), MappingQuality, 
    //  Readbases, ReadbasesQuality, ReadPositionRank, Strand]
    size_t sn = batchsamples_posinfomap_vector.size();
    std::vector<int> mapq;                     mapq.reserve(sn);
    std::vector<std::string> map_read_bases;   map_read_bases.reserve(sn);
    std::vector<char> map_read_base_qualities; map_read_base_qualities.reserve(sn);
    std::vector<int> read_pos_ranks;           read_pos_ranks.reserve(sn);
    std::vector<char> map_strands;             map_strands.reserve(sn);

    for (uint32_t pos(reg_start); pos < reg_end+1; ++pos) {

        PosMap::const_iterator pos_it;
        uint32_t depth = 0;
        for (size_t i = 0; i < sn; i++) {
            pos_it = batchsamples_posinfomap_vector[i].find(pos);
            if (pos_it != batchsamples_posinfomap_vector[i].end()) {

                ++depth;
                if (pos_it->second.ref_id != ref_id || pos_it->second.ref_pos != pos)
                    throw std::runtime_error("[ERROR] reference id or position not match.");
                
                mapq.push_back(pos_it->second.mapq);
                map_read_bases.push_back(pos_it->second.read_base);
                map_read_base_qualities.push_back(pos_it->second.read_base_qual);
                read_pos_ranks.push_back(pos_it->second.rpr);
                map_strands.push_back(pos_it->second.map_strand);

            } else {
                mapq.push_back(0);
                map_read_bases.push_back("N");
                map_read_base_qualities.push_back(BASE_Q0_ASCII);
                read_pos_ranks.push_back(0);
                map_strands.push_back('.');
            }
        }

        std::string out = ref_id + "\t" + 
                          ngslib::tostring(pos) + "\t" + fa_seq[pos-1] + "\t" + 
                          ngslib::tostring(depth)                      + "\t" + 
                          ngslib::join(mapq, ",")                      + "\t" + 
                          ngslib::join(map_read_bases, ",")            + "\t" + 
                          ngslib::join(map_read_base_qualities, ",")   + "\t" +
                          ngslib::join(read_pos_ranks, ",")            + "\t" + 
                          ngslib::join(map_strands, ",")               + "\n";
        
        // write to file and check is successful or not.
        if (bgzf_write(obf, out.c_str(), out.length()) != out.length()) {
            throw std::runtime_error("[ERROR] fail to write data");
        }

        mapq.clear();
        map_read_bases.clear();
        map_read_base_qualities.clear();
        read_pos_ranks.clear();
        map_strands.clear();
    }

    return;
}

bool __create_a_batchfile(const std::vector<std::string> batch_align_files,  // Not a modifiable value
                          const std::vector<std::string> batch_sample_ids,   // Not a modifiable value
                          const std::string &fa_seq,                         // Not a modifiable value
                          const ngslib::GenomeRegionTuple genome_region,
                          const int mapq_thd,                   // mapping quality threshold
                          const std::string output_batch_file)  // output batchfile name
// 原为 BaseTypeRunner 的成员函数，未掌握如何将该函数指针传入 ThreadPool，遂作罢，后再改。
{   
    // This value affected the computing memory, could be set to 1000000, 20 just for test
    static const uint32_t STEP_REGION_LEN = 20;
    clock_t start_time = clock();

    std::string ref_id; uint32_t reg_start, reg_end;
    std::tie(ref_id, reg_start, reg_end) = genome_region;  // 1-based

    PosMapVector batchsamples_posinfomap_vector;
    batchsamples_posinfomap_vector.reserve(batch_align_files.size());  //  pre-set the capacity

    BGZF *obf = bgzf_open(output_batch_file.c_str(), "w"); // output file handle of output_batch_file
    if (!obf) {
        throw std::runtime_error("[ERROR] " + output_batch_file + " open failure.");
    }

    // write header inforamtion to batchfile
    std::string bf_header = "##fileformat=BaseVarBatchFile_v1.0\n" 
                            "##SampleIDs=" + ngslib::join(batch_sample_ids, ",") + "\n" + 
                            "#CHROM\tPOS\tREF\tDepth(CoveredSample)\tMappingQuality\t"
                            "Readbases\tReadbasesQuality\tReadPositionRank\tStrand\n";
    if (bgzf_write(obf, bf_header.c_str(), bf_header.length()) != bf_header.length())
        throw std::runtime_error("[ERROR] fail to write data");

    uint32_t sub_reg_start, sub_reg_end;
    bool is_not_empty = false;
    for (uint32_t i(reg_start), j(0); i < reg_end + 1; i += STEP_REGION_LEN, ++j) {
        sub_reg_start = i;
        sub_reg_end = sub_reg_start + STEP_REGION_LEN - 1 > reg_end ? reg_end : sub_reg_start + STEP_REGION_LEN - 1;

std::cout << j << " - " << ref_id << ":" << sub_reg_start << "-" << sub_reg_end << "\n";
        is_not_empty = __fetch_base_in_region(batch_align_files, fa_seq, mapq_thd, 
                                              std::make_tuple(ref_id, sub_reg_start, sub_reg_end),
                                              batchsamples_posinfomap_vector);  // 传引用，省内存，得数据

        /* Output batchfile, no matter 'batchsamples_posinfomap_vector' is empty or not. */
        __write_record_to_batchfile(batchsamples_posinfomap_vector, fa_seq,
                                    std::make_tuple(ref_id, sub_reg_start, sub_reg_end), 
                                    obf);

        batchsamples_posinfomap_vector.clear();  // 清空，为下个循环做准备
    }

    int is_cl = bgzf_close(obf);  // 关闭文件
    if (is_cl < 0) {
        std::cout << "[WARNING] " + output_batch_file + " fail close.\n";
    }

    // Time information
    time_t now = time(0);
    std::string ct(ctime(&now));
    ct.pop_back();  // rm the trailing '\n' put by `asctime`

    std::cout << "[INFO] " + ct + ". Done for creating batchfile " << output_batch_file << ", " 
              << (double)(clock() - start_time) / CLOCKS_PER_SEC   << " seconds elapsed.\n"
              << std::endl;

    return is_not_empty;
}

bool __fetch_base_in_region(const std::vector<std::string> &batch_align_files,  
                            const std::string &fa_seq,  // must be the whole chromosome sequence  
                            const int mapq_thd,         // mapping quality threshold
                            const ngslib::GenomeRegionTuple genome_region,
                            PosMapVector &batchsamples_posinfomap_vector)  
{
    static const uint32_t REG_EXPEND_SIZE = 200; // only using here, In case of missing the overlap reads

    std::string ref_id; uint32_t reg_start, reg_end;
    std::tie(ref_id, reg_start, reg_end) = genome_region;  // 1-based

    uint32_t exp_reg_start = reg_start > REG_EXPEND_SIZE ? reg_start - REG_EXPEND_SIZE : 1;
    uint32_t exp_reg_end   = reg_end + REG_EXPEND_SIZE;
    std::string exp_regstr = ref_id + ":" + ngslib::tostring(exp_reg_start) + "-" + ngslib::tostring(exp_reg_end);

    // Loop all alignment files
    bool is_not_empty = false;
    for(size_t i(0); i < batch_align_files.size(); ++i) {
        ngslib::Bam bf(batch_align_files[i], "r");  // open bamfile in reading mode (one sample, one bamfile)

        // 位点信息存入该变量, 且由于是按区间读取比对数据，key 值无需再包含 ref_id，因为已经不言自明。
        PosMap sample_posinfo_map;

        if (bf.fetch(exp_regstr)) { // Set 'bf' only fetch alignment reads in 'exp_regstr'.
uint32_t read_count = 0;
            hts_pos_t map_ref_start, map_ref_end;  // hts_pos_t is uint64_t
            std::vector<ngslib::BamRecord> sample_target_reads; 
            ngslib::BamRecord al;       // alignment read

            while (bf.next(al) >= 0) {  // -1 => hit the end of alignement file.
++read_count;
                if (al.mapq() < mapq_thd || al.is_duplicate() || al.is_qc_fail()) continue;
                map_ref_start = al.map_ref_start_pos() + 1;  // al.map_ref_start_pos() is 0-based, convert to 1-based
                map_ref_end   = al.map_ref_end_pos();        // al.map_ref_end_pos() is 1-based

                // Only fetch reads which in [reg_start, reg_end]
                if (reg_start > map_ref_end) continue;
                if (reg_end < map_ref_start) break;

                sample_target_reads.push_back(al);  // record the proper reads of sample
            }

            sample_posinfo_map.clear();  // make sure it's empty
            if (sample_target_reads.size() > 0) {
                // get alignment information of [i] sample.
                __seek_position(sample_target_reads, fa_seq, genome_region, sample_posinfo_map);
            }
std::cout << "* " + exp_regstr + " total read count: " << read_count << ". Hit read count: " << sample_target_reads.size() << "\n\n";
        }

        if (!is_not_empty && !sample_posinfo_map.empty()) { 
            // at least one sample has data in this region
            is_not_empty = true; 
        }

        // Push it into 'batchsamples_posinfomap_vector' even if 'sample_posinfo_map' is empty, 
        // make sure 'batchsamples_posinfomap_vector' has the same size as `batch_align_files`
        batchsamples_posinfomap_vector.push_back(sample_posinfo_map);
    }

    if (batchsamples_posinfomap_vector.size() != batch_align_files.size())
        throw std::runtime_error("[basetype.cpp::__fetch_base_in_region] 'pos_batchinfo_vector.size()' "
                                 "should be the same as 'batch_align_files.size()'");

    return is_not_empty;  // no cover reads in 'genome_region' if empty.
}

void __seek_position(const std::vector<ngslib::BamRecord> &sample_map_reads,
                     const std::string &fa_seq,   // must be the whole chromosome sequence
                     const ngslib::GenomeRegionTuple genome_region,
                     PosMap &sample_posinfo_map)
{
    if (!sample_posinfo_map.empty())
        throw std::runtime_error("[basetype.cpp::__seek_position] 'sample_posinfo_map' must be empty.");

    std::string ref_id; uint32_t reg_start, reg_end;
    std::tie(ref_id, reg_start, reg_end) = genome_region;  // 1-based

    AlignBaseInfo align_base_info;
    align_base_info.ref_id = ref_id;

    // A vector of: (cigar_op, read position, reference position, read base, read_qual, reference base)
    std::vector<ngslib::ReadAlignedPair> aligned_pairs;
    for(size_t i(0); i < sample_map_reads.size(); ++i) {
std::cout << sample_map_reads[i] << "\n";

        align_base_info.map_strand = sample_map_reads[i].map_strand();  // '*', '-' or '+'
        align_base_info.mapq = sample_map_reads[i].mapq();

        aligned_pairs = sample_map_reads[i].get_aligned_pairs(fa_seq);
        char mean_qqual_char = int(sample_map_reads[i].mean_qqual());
        uint32_t map_ref_pos;
        for (size_t i(0); i < aligned_pairs.size(); ++i) {
            // Todo: The data of 'align_base_info' and 'aligned_pairs[i]' is similar, 
            // could we just use 'aligned_pairs[i]' to replace 'align_base_info'?

            map_ref_pos = aligned_pairs[i].ref_pos + 1;  // ref_pos is 0-based, convert to 1-based;

            if (reg_end < map_ref_pos) break;
            if (reg_start > map_ref_pos) continue;

            // 'BAM_XXX' are macros for CIGAR, which defined in 'htslib/sam.h'
            if (aligned_pairs[i].op == BAM_CMATCH ||  /* CIGAR: M */ 
                aligned_pairs[i].op == BAM_CEQUAL ||  /* CIGAR: = */
                aligned_pairs[i].op == BAM_CDIFF)     /* CIGAR: X */
            {
                // One character
                align_base_info.ref_base       = aligned_pairs[i].ref_base[0];
                align_base_info.read_base      = aligned_pairs[i].read_base[0];
                align_base_info.read_base_qual = aligned_pairs[i].read_qual[0];
            } else if (aligned_pairs[i].op == BAM_CINS) {  /* CIGAR: I */
                if (!aligned_pairs[i].ref_base.empty()) {
                    std::cerr << sample_map_reads[i] << "\n";
                    throw std::runtime_error("[ERROR] We got reference base in insertion region.");
                }

                // roll back one position to the left side of insertion break point.
                --map_ref_pos;
                align_base_info.ref_base       = fa_seq[aligned_pairs[i].ref_pos-1]; // break point's ref base
                align_base_info.read_base      = fa_seq[aligned_pairs[i].ref_pos-1] + aligned_pairs[i].read_base;
                align_base_info.read_base_qual = mean_qqual_char;  // set to be mean quality of the whole read
            } else if (aligned_pairs[i].op == BAM_CDEL) {  /* CIGAR: D */
                if (!aligned_pairs[i].read_base.empty()) {
                    std::cerr << sample_map_reads[i] << "\n";
                    throw std::runtime_error("[ERROR] We got read bases in deletion region.");
                }

                // roll back one position to the left side of deletion break point.
                --map_ref_pos;
                align_base_info.ref_base       = fa_seq[aligned_pairs[i].ref_pos-1] + aligned_pairs[i].ref_base;
                align_base_info.read_base      = fa_seq[aligned_pairs[i].ref_pos-1]; // break point's ref base
                align_base_info.read_base_qual = mean_qqual_char;  // set to be mean quality of the whole read
            } else { 
                // Skip in basevar for other kind of CIGAR symbals.
                continue;
            }
            align_base_info.ref_pos = map_ref_pos;

            // qpos is 0-based, conver to 1-based to set the rank of base on read.
            align_base_info.rpr = aligned_pairs[i].qpos + 1;

            if (sample_posinfo_map.find(map_ref_pos) == sample_posinfo_map.end()) {
                // Just get the base from first read which aligned on this ref_pos,
                // no matter the first one it's indel or not.

                // {ref_pos (no need to add ref_id in the key) => map info}
                sample_posinfo_map.insert({map_ref_pos, align_base_info});
            }

std::cout << " - "  << aligned_pairs[i].op << " - [" << ref_id << ", " << align_base_info.ref_pos << ", " 
          << align_base_info.map_strand    << ", "   << align_base_info.ref_base << "-" << align_base_info.read_base << "] - [" 
          << aligned_pairs[i].qpos << ", " << aligned_pairs[i].read_base << ", " << aligned_pairs[i].read_qual << "]\n";
        }
std::cout << "\n";
    }

    return;
}