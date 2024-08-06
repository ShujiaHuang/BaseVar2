/**
 * @file basetype.h
 * @brief 
 * 
 * @author Shujia Huang
 * @date 2018-08-01
 * 
 */
#ifndef __INCLUDE_BASETYPE_H__
#define __INCLUDE_BASETYPE_H__

#include <iostream>
#include <string>
#include <vector>
#include <getopt.h>

#include "Fasta.h"
#include "utils.h"

static const std::string __BASETYPE_USAGE = 
    "About: Calling variants by BaseVar.\n" 
    "Usage: basevar basetype [options] <-R reference.fa>|<-R reference.fa.gz> <--output-vcf> <--output-cvg> " 
    "[-I input] ...\n\n" 
    "optional arguments:\n" 
    "  -I, --input=FILE             BAM/SAM/CRAM file containing reads.\n"
    "  -L, --align-file-list=FILE   BAM/CRAM files list, one file per row.\n"
    "  -R, --reference FILE         Input reference fasta file.\n\n"

    "  -m, --min-af=float           Setting prior precision of MAF and skip uneffective caller positions.\n"
    "                               Usually you can set it to be min(0.001, 100/x), x is the number of your\n"
    "                               input BAM files.[min(0.001, 100/x, cmm.MINAF)]. Probably you don't need\n"
    "                               to take care about this parameter.\n\n" 
    "  -q, --mapq=INT               Only include reads with mapping quality >= INT. [10]\n"
    "  -B, --batch-count=INT        INT simples per batchfile. [200]\n" 
    "  -t, --thread=INT             Number of threads. [4]\n\n"

    "  -r, --regions=chr:start-end  Skip positions which not in these regions. This parameter could be a\n"
    "                               list of comma deleimited genome regions(e.g.: chr:start-end) or a\n"
    "                               file contain the list of regions.\n"
    "  -G, --pop-group=FILE         Calculating the allele frequency for specific population.\n" 
    "  --output-vcf FILE            Output VCF file. If not provide will skip variants discovery and just\n"
    "                               output position coverage file which filename is provided by --output-cvg.\n"
    "  --output-cvg FILE            Output position coverage file.\n\n"

    "  --filename-has-samplename    If the name of bamfile is something like 'SampleID.xxxx.bam',\n"
    "                               you can set this parameter to save a lot of time during get the\n"
    "                               sample id from BAM header.\n"
    "  --smart-rerun                Rerun process by checking batchfiles\n"
    "  -h, --help                   Show this help message and exit"; 

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

    // {"output-batch-file", required_argument, NULL, 0},  // not use?
    {"filename-has-samplename", no_argument, NULL, '3'},
    {"smart-rerun",             no_argument, NULL, '4'},
    {"help",                    no_argument, NULL, 'h'},
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

    bool smart_rerun;                   // Smart rerun by checking batchfiles
    bool filename_has_samplename;       // sample name in file name

    // Set default argument
    BaseTypeARGS(): min_af(0.01), mapq(10), batchcount(200), thread_num(4), 
                    smart_rerun(false), filename_has_samplename(false) {}
};

/**
 * @brief BaseTypeRunner class
 */
class BaseTypeRunner {

private:
    BaseTypeARGS *_args;                   // all the commandline options
    std::vector<std::string> _samples_id;  // sample ID of alignment files (BAM/CRAM/SAM)
    std::vector<ngslib::GenomeRegionTuple> _calling_intervals;  // vector of calling regions

    void _get_bamfile_list();
    void _get_calling_interval();  // load the calling region from input
    void _get_sample_id_from_bam();
    ngslib::GenomeRegionTuple _make_gregiontuple(std::string gregion);

    BaseTypeRunner(const BaseTypeRunner &b) = delete;             // reject using copy constructor (C++11 style).
    BaseTypeRunner &operator=(const BaseTypeRunner &b) = delete;  // reject using copy/assignment operator (C++11 style).

public:
    ngslib::Fasta reference;  // public variable

    // default constructor
    BaseTypeRunner() : _args(NULL) {}
    BaseTypeRunner(int cmdline_argc, char *cmdline_argv[]) {
        if (cmdline_argc < 2) {
            std::cerr << usage() << std::endl;
            exit(1);
        }
        set_arguments(cmdline_argc, cmdline_argv); 
    }
    // Destroy the malloc'ed BasTypeArgs structure
    ~BaseTypeRunner(){ if(_args){delete _args; _args = NULL;} }

    // Common functions
    std::string usage() const {return __BASETYPE_USAGE;}
    void set_arguments(int cmdline_argc, char *cmdline_argv[]);
    void print_calling_interval();

    // run variant calling process

};  // BaseTypeRunner class

#endif
