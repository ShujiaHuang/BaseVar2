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

#include <htslib/tbx.h>

#include "io/fasta.h"
#include "io/iobgzf.h"
#include "io/bam.h"
#include "io/utils.h"

#include "basetype_utils.h"
#include "basetype.h"
#include "algorithm.h"
#include "version.h"

#include "external/thread_pool.h"

static const bool IS_DELETE_CACHE_BATCHFILE = true;

/**
 * @brief BaseTypeRunner class
 */
class BaseTypeRunner {
public:
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

        bool filename_has_samplename;       // sample name in file name
        bool smart_rerun;                   // Smart rerun by checking batchfiles

        // Set default argument
        BaseTypeARGS(): min_af(0.001), 
                        min_mapq(5), 
                        min_baseq(10), 
                        batchcount(500), 
                        thread_num(std::thread::hardware_concurrency()), 
                        smart_rerun(false), 
                        filename_has_samplename(false) {}
    };
    ngslib::Fasta reference;  // public variable

    // default constructor
    BaseTypeRunner() : _args(NULL) {}
    BaseTypeRunner(int cmdline_argc, char *cmdline_argv[]);
    
    // Destroy the malloc'ed BasTypeArgs structure
    ~BaseTypeRunner(){ if(_args){delete _args; _args = NULL;} }

    // Common functions
    const std::string usage() const;
    void print_calling_interval();

    // Run the variant calling process
    void run() { _variant_caller_process(); }

private:
    std::string _cmdline_string;                             // save the commandline options
    BaseTypeARGS *_args;                                     // Commandline options
    
    std::vector<std::string> _samples_id;                    // sample ID of alignment files (BAM/CRAM/SAM)
                                                             // `_samples_id` and `input_bf` have the same order 
    std::map<std::string, std::vector<size_t>> _groups_idx;  // sample group: group => samples index
    std::vector<ngslib::GenomeRegion> _calling_intervals;    // vector of calling regions

    // templary output files
    std::vector<std::string> _sub_out_vcf;

    void _get_calling_interval();  // load the calling region from input
    void _get_sample_id_from_bam();
    void _get_popgroup_info();
    ngslib::GenomeRegion _make_gregion_region(std::string gregion);

    // For variant calling
    void _variant_caller_process();
    
    /**
     * @brief Create a batch of temp files for variant discovery (could be deleted when the jobs done).
     * 
     * @param genome_region 
     * @return std::vector<std::string> 
     * 
     */
    std::vector<std::string> _create_batchfiles(const ngslib::GenomeRegion genome_region, 
                                                const std::string bf_prefix);
    void _variants_discovery(const std::vector<std::string> &batchfiles, 
                             const ngslib::GenomeRegion genome_region,
                             const std::string sub_vcf_fn);

    BaseTypeRunner(const BaseTypeRunner &) = delete;             // reject using copy constructor (C++11 style).
    BaseTypeRunner &operator=(const BaseTypeRunner &) = delete;  // reject using copy/assignment operator (C++11 style).

};  // BaseTypeRunner class

// This function is only used by BaseTypeRunner::_create_batchfiles
bool __create_a_batchfile(const std::vector<std::string> batch_align_files, // Not a modifiable value
                          const std::vector<std::string> batch_sample_ids,  // Not a modifiable value
                          const std::string &fa_seq,                        // Not a modifiable value
                          const ngslib::GenomeRegion genome_region,         // 切割该区间
                          const int min_mapq,                               // mapping quality threshold
                          const int min_baseq,                              // base quality threshold
                          const std::string output_batch_file);             // output batchfile name

bool __fetch_base_in_region(const std::vector<std::string> &batch_align_files,
                            const std::string &fa_seq,                   
                            const int min_mapq,
                            const int min_baseq,  // minimum base quality
                            const ngslib::GenomeRegion target_genome_region,  // 获取该区间内的 read
                            PosMapVector &batchsamples_posinfomap_vector);    // 信息存入该变量

void __seek_position(const std::vector<ngslib::BamRecord> &sample_map_reads,  // ngslib::BamRecord include by 'bam.h'
                     const std::string &fa_seq,
                     const ngslib::GenomeRegion target_genome_region,         // 获取该区间内所有位点的碱基比对信息，该参数和 '__fetch_base_in_region' 中一样 
                     const int min_baseq,  // minimum base quality
                     PosMap &sample_posinfo_map);

void __write_record_to_batchfile(const PosMapVector &batchsamples_posinfomap_vector, 
                                 const std::string &fa_seq,
                                 const ngslib::GenomeRegion target_genome_region,  // 该参数和 __seek_position 中一样 
                                 ngslib::BGZFile &obf);

// Get sample id from batchfiles header.
std::vector<std::string> _get_sampleid_from_batchfiles(const std::vector<std::string> &batchfiles);

// A unit for calling variants and let it run in a thread.
bool _variant_calling_unit(const std::vector<std::string> &batchfiles, 
                           const std::vector<std::string> &sample_ids,
                           const std::map<std::string, std::vector<size_t>> & group_smp_idx,
                           const double min_af,
                           const std::string reg_str,  // genome region format like samtools
                           const std::string vcf_fn);

bool _basevar_caller(const std::vector<std::string> &smp_bf_line_vector, 
                     const std::map<std::string, std::vector<size_t>> &group_smp_idx,
                     const double min_af, 
                     size_t n_sample, 
                     ngslib::BGZFile &vcf_hd);

std::pair<BaseType, VariantInfo> _basetype_caller_unit(const std::vector<BaseType::BatchInfo> &samples_batchinfo_vector, 
                                                       const double min_af,
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
#endif
