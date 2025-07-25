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

#include <htslib/bgzf.h>
#include <htslib/tbx.h>

#include "io/fasta.h"
#include "io/bam.h"
#include "io/utils.h"
#include "external/robin_hood.h"
#include "external/thread_pool.h"

#include "basetype.h"
#include "basetype_utils.h"
#include "version.h"

static const bool IS_DELETE_CACHE_BATCHFILE = true;

struct BaseTypeARGS {
    /* Variables for all the commandline options of BaseType */
    std::vector<std::string> input_bf;  // BAM/SAM/CRAM file, a vector
    std::string in_bamfilelist;         // BAM/CRAM files list, one file per row
    std::string reference;              // Input reference fasta file

    float min_af;                       // Setting prior precision of MAF and skip uneffective caller positions
    int min_baseq;                      // minimum base quality
    int min_mapq;                       // mapping quality
    int batchcount;                     // INT simples per batchfile
    int thread_num;                     // number of threads

    std::string regions;                // Interval regions
    std::string pop_group_file;         // Specific population
    std::string output_vcf;             // Output VCF file
    std::string output_cvg;             // Output coverage file

    bool filename_has_samplename;       // sample name in file name
    bool smart_rerun;                   // Smart rerun by checking batchfiles

    // Set default argument
    BaseTypeARGS(): min_af(0.001), 
                    min_mapq(10), 
                    min_baseq(5), 
                    batchcount(200), 
                    thread_num(std::thread::hardware_concurrency()), 
                    smart_rerun(false), 
                    filename_has_samplename(false) {}
};

/**
 * @brief BaseTypeRunner class
 */
class BaseTypeRunner {

private:
    std::string _cmdline_string;                             // save the commandline options
    BaseTypeARGS *_args;                                     // Commandline options
    
    std::vector<std::string> _samples_id;                    // sample ID of alignment files (BAM/CRAM/SAM)
                                                             // `_samples_id` and `input_bf` have the same order 
    std::map<std::string, std::vector<size_t>> _groups_idx;  // sample group: group => samples index
    std::vector<ngslib::GenomeRegion> _calling_intervals;    // vector of calling regions

    // templary output files
    std::vector<std::string> _sub_out_vcf, _sub_out_cvg;

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
                             const std::string sub_vcf_fn,
                             const std::string sub_cvg_fn);

    BaseTypeRunner(const BaseTypeRunner &) = delete;             // reject using copy constructor (C++11 style).
    BaseTypeRunner &operator=(const BaseTypeRunner &) = delete;  // reject using copy/assignment operator (C++11 style).

public:
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

};  // BaseTypeRunner class

// Mainly use for create batchfile
struct AlignBaseInfo {
    std::string ref_id;
    uint32_t ref_pos;
    std::string ref_base;   // reference base

    std::string read_base;  // read base
    int mapq;               // mapping quality
    int rpr;                // read position rank
    char map_strand;        // mapping reference strand, should be one of '*', '-', '+'
    char read_base_qual;    // read base quality (get mean quality of read if Indels, 
                            // I don't care about Indels for NIPT data)
};

typedef robin_hood::unordered_map<uint32_t, AlignBaseInfo> PosMap;           // give a short name to this type
typedef std::vector<PosMap> PosMapVector;

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
                                 BGZF *obf);

// A unit for calling variants and let it run in a thread.
bool _variant_calling_unit(const std::vector<std::string> &batchfiles, 
                           const std::vector<std::string> &sample_ids,
                           const std::map<std::string, std::vector<size_t>> & group_smp_idx,
                           const double min_af,
                           const std::string region,  // genome region format like samtools
                           const std::string tmp_vcf_fn,
                           const std::string tmp_cvg_fn);

// Get sample id from batchfiles header.
std::vector<std::string> _get_sampleid_from_batchfiles(const std::vector<std::string> &batchfiles);

bool _basevar_caller(const std::vector<std::string> &smp_bf_line_vector, 
                     const std::map<std::string, std::vector<size_t>> &group_smp_idx,
                     double min_af,
                     size_t n_sample, 
                     BGZF *vcf_hd, 
                     BGZF *cvg_hd);

const BaseType __gb(const BatchInfo *smp_bi, 
              const std::vector<size_t> group_idx, 
              const std::vector<char> &basecombination, 
              double min_af);

const BatchInfo __get_group_batchinfo(const BatchInfo *smp_bi, const std::vector<size_t> group_idx);

void _out_vcf_line(const BaseType &bt, 
                   const std::map<std::string, BaseType> &group_bt, 
                   const BatchInfo *smp_bi, 
                   BGZF *vcf_hd);

void _out_cvg_line(const BatchInfo *smp_bi, 
                   const std::map<std::string, std::vector<size_t>> & group_smp_idx, 
                   BGZF *cvg_hd);

typedef std::tuple<std::string, std::map<char, int>> IndelTuple;
IndelTuple __base_depth_and_indel(const std::vector<std::string> &align_bases);

#endif
