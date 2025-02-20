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

#include "io/fasta.h"
#include "io/bam.h"
#include "io/utils.h"
#include "external/robin_hood.h"

#include "basetype.h"
#include "basetype_utils.h"

static const bool IS_DELETE_CACHE_BATCHFILE = true;

static const std::string __BASETYPE_USAGE = 
    "About: Call variants and estimate allele frequency by BaseVar.\n" 
    "Usage: basevar basetype [options] <-R Fasta> <--output-vcf> <--output-cvg> [-I input] ...\n\n" 
    "optional arguments:\n" 
    "  -I, --input=FILE             BAM/CRAM file containing reads.\n"
    "  -L, --align-file-list=FILE   BAM/CRAM files list, one file per row.\n"
    "  -R, --reference FILE         Input reference fasta file.\n\n"

    "  -m, --min-af=float           Setting prior precision of MAF and skip ineffective caller positions,\n"
    "                               a typical approach involves setting it to min(0.001, 100/x), where x \n"
    "                               represents the number of input BAM files [min(0.001, 100/x)]. In most\n"
    "                               cases, users need not be overly concerned about this parameter, as it \n"
    "                               is generally handled automatically by the program.\n"
    "  -q, --mapq=INT               Only include reads with mapping quality >= INT. [10]\n"
    "  -B, --batch-count=INT        INT simples per batchfile. [200]\n" 
    "  -t, --thread=INT             Number of threads. [4]\n\n"

    "  -G, --pop-group=FILE         Calculating the allele frequency for specific population.\n" 
    "  -r, --regions=chr:start-end  Skip positions which not in these regions. This parameter could be a list\n"
    "                               of comma deleimited genome regions(e.g.: chr:start-end).\n"
    "  --output-vcf FILE            Output VCF file.\n"
    "  --output-cvg FILE            Output position coverage file.\n\n"

    "  --filename-has-samplename    If the name of bamfile is something like 'SampleID.xxxx.bam', set this\n"
    "                               argrument could save a lot of time during get the sample id from BAMfile.\n"
    "  --smart-rerun                Rerun process by checking batchfiles.\n"
    "  -h, --help                   Show this help message and exit."; 

static const struct option BASETYPE_CMDLINE_LOPTS[] = {
    // Optional arguments to long style command line parameters require 'equals sign' (=). 
    // https://stackoverflow.com/questions/1052746/getopt-does-not-parse-optional-arguments-to-parameters
    {"input",           optional_argument, NULL, 'I'},
    {"align-file-list", optional_argument, NULL, 'L'},
    {"reference",       required_argument, NULL, 'R'},

    {"min-af",      optional_argument, NULL, 'm'},
    {"mapq",        optional_argument, NULL, 'q'},
    {"batch-count", optional_argument, NULL, 'B'},
    {"thread",      optional_argument, NULL, 't'},

    {"regions",     optional_argument, NULL, 'r'},
    {"positions",   optional_argument, NULL, 'p'},
    {"pop-group",   optional_argument, NULL, 'G'},  // Special parameter for calculating specific population allele frequence
    {"output-vcf",  required_argument, NULL, '1'},
    {"output-cvg",  required_argument, NULL, '2'},

    // {"output-batch-file", required_argument, NULL, 0},  not use?
    {"filename-has-samplename", no_argument, NULL, '3'},
    {"smart-rerun",             no_argument, NULL, '4'},
    {"help",                    no_argument, NULL, 'h'},

    // must set this value
    {0, 0, 0, 0}
};

struct BaseTypeARGS {
    /* Variables for all the commandline options of BaseType */
    std::vector<std::string> input_bf;  // BAM/SAM/CRAM file, a vector
    std::string in_bamfilelist;         // BAM/CRAM files list, one file per row
    std::string reference;              // Input reference fasta file

    float min_af;                       // Setting prior precision of MAF and skip uneffective caller positions
    int mapq;                           // mapping quality
    int batchcount;                     // INT simples per batchfile
    int thread_num;                     // number of threads

    std::string regions;                // Interval regions
    std::string pop_group_file;         // Specific population
    std::string output_vcf;             // Output VCF file
    std::string output_cvg;             // Output coverage file

    bool filename_has_samplename;       // sample name in file name
    bool smart_rerun;                   // Smart rerun by checking batchfiles

    // Set default argument
    BaseTypeARGS(): min_af(0.01), mapq(10), batchcount(200), thread_num(4), 
                    smart_rerun(false), filename_has_samplename(false) {}
};

/**
 * @brief BaseTypeRunner class
 */
class BaseTypeRunner {

private:
    BaseTypeARGS *_args;                                        // Commandline options
    std::vector<std::string> _samples_id;                       // sample ID of alignment files (BAM/CRAM/SAM)
                                                                // `_samples_id` and `input_bf` have the same order 
    std::map<std::string, std::vector<size_t>> _groups_idx;     // sample group: group => samples index
    std::vector<ngslib::GenomeRegionTuple> _calling_intervals;  // vector of calling regions

    // templary output files
    std::vector<std::string> _sub_out_vcf, _sub_out_cvg;

    void _get_calling_interval();  // load the calling region from input
    void _get_sample_id_from_bam();
    void _get_popgroup_info();
    ngslib::GenomeRegionTuple _make_gregion_tuple(std::string gregion);

    // For variant calling
    void _variant_caller_process();
    
    /**
     * @brief Create a batch of temp files for variant discovery (could be deleted when the jobs done).
     * 
     * @param genome_region 
     * @return std::vector<std::string> 
     * 
     */
    std::vector<std::string> _create_batchfiles(const ngslib::GenomeRegionTuple &genome_region, 
                                                const std::string bf_prefix);
    void _variants_discovery(const std::vector<std::string> &batchfiles, 
                             const ngslib::GenomeRegionTuple &genome_region,
                             const std::string sub_vcf_fn,
                             const std::string sub_cvg_fn);

    BaseTypeRunner(const BaseTypeRunner &) = delete;             // reject using copy constructor (C++11 style).
    BaseTypeRunner &operator=(const BaseTypeRunner &) = delete;  // reject using copy/assignment operator (C++11 style).

public:
    ngslib::Fasta reference;  // public variable

    // default constructor
    BaseTypeRunner() : _args(NULL) {}
    BaseTypeRunner(int cmdline_argc, char *cmdline_argv[]) { set_arguments(cmdline_argc, cmdline_argv); }
    
    // Destroy the malloc'ed BasTypeArgs structure
    ~BaseTypeRunner(){ if(_args){delete _args; _args = NULL;} }

    // Common functions
    std::string usage() const {return __BASETYPE_USAGE;}
    void set_arguments(int cmdline_argc, char *cmdline_argv[]);
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
bool __create_a_batchfile(const std::vector<std::string> batch_align_files,  // Not a modifiable value
                          const std::vector<std::string> batch_sample_ids,   // Not a modifiable value
                          const std::string &fa_seq,                         // Not a modifiable value
                          const ngslib::GenomeRegionTuple &genome_region,    // 切割该区间
                          const int mapq_thd,                                // mapping quality threshold
                          const std::string output_batch_file);              // output batchfile name

bool __fetch_base_in_region(const std::vector<std::string> &batch_align_files,
                            const std::string &fa_seq,                   
                            const int mapq_thd,
                            const ngslib::GenomeRegionTuple target_genome_region,  // 获取该区间内的 read
                            PosMapVector &batchsamples_posinfomap_vector);         // 信息存入该变量

void __seek_position(const std::vector<ngslib::BamRecord> &sample_map_reads,  // ngslib::BamRecord include by 'bam.h'
                     const std::string &fa_seq,
                     const ngslib::GenomeRegionTuple target_genome_region,    // 获取该区间内所有位点的碱基比对信息，该参数和 '__fetch_base_in_region' 中一样 
                     PosMap &sample_posinfo_map);

void __write_record_to_batchfile(const PosMapVector &batchsamples_posinfomap_vector, 
                                 const std::string &fa_seq,
                                 const ngslib::GenomeRegionTuple target_genome_region,  // 该参数和 __seek_position 中一样 
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
