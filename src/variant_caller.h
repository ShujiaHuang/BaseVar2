/**
 * @file basetype_caller.h
 * 
 * @brief BaseType main caller functions
 *  
 * @author Shujia Huang
 * @date 2018-08-29
 * 
 */
#ifndef __INCLUDE_BASETYPE_CALLER_H__
#define __INCLUDE_BASETYPE_CALLER_H__

#include <getopt.h>
#include <map>
#include <fstream>
#include <sstream>
#include <ctime>      // clock
#include <algorithm>  // std::min
#include <memory>     // std::unique_ptr

#include <htslib/tbx.h>

#include "io/fasta.h"
#include "io/iobgzf.h"
#include "io/bam.h"
#include "io/utils.h"

#include "caller_utils.h"
#include "basetype.h"
#include "algorithm.h"
#include "version.h"
#include "external/thread_pool.h"

static const bool IS_DELETE_CACHE_BATCHFILE = true;
static const char BASE_Q0_CHAR = '!'; // The ascii code of '!' character is 33

/**
 * @brief BaseTypeRunner class
 */
class BaseTypeRunner {
private:
    // Commandline arguments
    struct BaseTypeARGS {
        /* Variables for all the commandline options of BaseType */
        std::vector<std::string> input_bf;  // BAM/SAM/CRAM file, a vector
        std::string reference;              // Input reference fasta file

        float min_af;                       // Setting prior precision of MAF and skip uneffective caller positions
        int min_baseq;                      // minimum base quality
        int min_mapq;                       // mapping quality
        int batchcount;                     // INT simples per batchfile
        int thread_num;                     // number of threads

        std::string regions;                // Interval regions
        std::string pop_group_file;         // Specific population
        std::string output_vcf;             // Output VCF file
        // std::string output_mode;            // e.g, "w", "wb", "wz"

        bool filename_has_samplename;       // sample name in file name

        // Set default argument
        BaseTypeARGS(): min_af(0.001), 
                        min_mapq(5), 
                        min_baseq(10), 
                        batchcount(1000), 
                        thread_num(std::thread::hardware_concurrency()), 
                        filename_has_samplename(false) {}
    };

    ngslib::Fasta reference;
    BaseTypeARGS *_args;                                     // Commandline options
    std::string _cmdline_string;                             // save the commandline options
    std::vector<std::string> _samples_id;                    // sample ID of alignment files (BAM/CRAM/SAM)
                                                             // `_samples_id` and `input_bf` have the same order 
    std::map<std::string, std::vector<size_t>> _groups_idx;  // sample group: group => samples index
    std::vector<ngslib::GenomeRegion> _calling_intervals;    // vector of calling regions

    // templary output files
    std::vector<std::string> _sub_out_vcf;

    void _make_calling_interval();  // load the calling region from input
    void _get_sample_id();
    void _get_popgroup_info();
    ngslib::GenomeRegion _make_gregion_region(const std::string &gregion);

    // Functions for calling variants in specific region
    bool _call_variants_in_region(const ngslib::GenomeRegion gr, const std::string out_vcf_fn);
    bool _call_variant_unit(const ngslib::GenomeRegion sub_gr, // 局部变量，会变，必拷贝，不可传引用，否则线程执行时将丢失该值 
                            const std::string &fa_seq,         // must be the whole chromosome sequence
                            const std::string &out_vcf_fn);

    bool _fetch_base_in_region(const std::vector<std::string> &batch_align_files,
                               const std::string &fa_seq,
                               const ngslib::GenomeRegion target_genome_region,  // 获取该区间内的 read
                               PosMapVector &batchsamples_posinfomap_vector);    // 信息存入该变量

    void _seek_position(const std::vector<ngslib::BamRecord> &sample_map_reads,  // ngslib::BamRecord include by 'bam.h'
                        const std::string &fa_seq,   // must be the whole chromosome sequence
                        const ngslib::GenomeRegion target_genome_region, // 获取该区间内所有位点的碱基比对信息，该参数和 '_fetch_base_in_region' 中一样 
                        PosMap &sample_posinfo_map);

    bool _basevar_caller(std::vector<BaseType::BatchInfo> &all_smps_bi_vector, //const std::vector<std::string> &smp_bf_line_vector, 
                         const std::map<std::string, std::vector<size_t>> &group_smp_idx,
                         size_t n_sample, 
                         ngslib::BGZFile &vcf_hd);

    std::pair<BaseType, VariantInfo> _basetype_caller_unit(
        const std::vector<BaseType::BatchInfo> &samples_batchinfo_vector, 
        const std::vector<size_t> group_idx = std::vector<size_t>(),
        const std::vector<std::string> basecombination = std::vector<std::string>());

    /**
     * @brief Get the variant information object
     * 
     * @param bt BaseType
     * @param smp_bi BaseType::BatchInfo 
     * @return VariantInfo 
     * 
     */
    VariantInfo get_pos_variant_info(const BaseType &bt, const BaseType::BatchInfo *smp_bi);
    VCFRecord _vcfrecord_in_pos(const std::vector<BaseType::BatchInfo> &samples_batchinfo_vector, 
                                const VariantInfo &global_variant_info,
                                const std::map<std::string, BaseType> &group_bt, 
                                AlleleInfo &ai);

    BaseTypeRunner(const BaseTypeRunner &) = delete;             // reject using copy constructor (C++11 style).
    BaseTypeRunner &operator=(const BaseTypeRunner &) = delete;  // reject using copy/assignment operator (C++11 style).

public:
    // default constructor
    BaseTypeRunner() : _args(NULL) {}
    BaseTypeRunner(int cmdline_argc, char *cmdline_argv[]);
    
    // Destroy the malloc'ed BasTypeArgs structure
    ~BaseTypeRunner(){ if(_args){delete _args; _args = NULL;} }

    // Common functions
    const std::string usage() const;

    // Run the whole variant calling processes
    int run();

};  // BaseTypeRunner class




#endif
