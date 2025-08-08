/**
 * @file basetype_utils.h
 * 
 * @brief define the global variable and util codes specific for basetype.h
 *  
 * @author Shujia Huang
 * @date 2018-08-29
 * 
 */
#ifndef __INCLUDE_BASETYPE_UTILS_H__
#define __INCLUDE_BASETYPE_UTILS_H__

#include <string>
#include <vector>
#include <iostream>

#include <htslib/bgzf.h>
#include <htslib/kstring.h>

#include "io/fasta.h"
#include "io/iobgzf.h"
#include "io/utils.h"
#include "basetype.h"

#include "external/robin_hood.h"  // robin_hood::unordered_map, robin_hood::unordered_set
#include "algorithm.h"

struct AlignBase {
    std::string ref_base;   // reference base
    std::string read_base;  // read base
    char base_qual;         // read base quality, assign mean quality of seq if it's Indels 
    int rpr;                // read position rank, 1-based, the first base in read has rpr = 1, the second base has rpr = 2, etc.

    int mapq;               // mapping quality
    char map_strand;        // mapping reference strand, should be one of '*', '-', '+'  

    std::string to_string() const {
        std::stringstream ss;
        // Required fields
        ss << ref_base  << "\t"
           << read_base << "\t"
           << base_qual << "\t"
           << rpr       << "\t"
           << mapq      << "\t"
           << map_strand;

        return ss.str();
    }
};

struct AlignInfo {
    std::string ref_id;
    uint32_t ref_pos;
    std::vector<AlignBase> align_bases;

    AlignInfo() : ref_id(""), ref_pos(0) {}
    AlignInfo(const std::string& rid, uint32_t pos) : ref_id(rid), ref_pos(pos) {}
};

typedef robin_hood::unordered_map<uint32_t, AlignInfo> PosMap;    // give a short name to this type
typedef std::vector<PosMap> PosMapVector;
typedef struct {
    int fwd, rev;
    double fs;   // Phred-scaled p-value using Fisher's exact test to detect strand bias
    double sor;  // Strand bias estimated by the Symmetric Odds Ratio test
} StrandBiasInfo;

struct VariantInfo {
    std::string ref_id;
    uint32_t ref_pos;
    int total_depth;
    int qual;  // quality score

    size_t major_allele_idx;                 // The index of major allele in `ale_bases` vector
    std::vector<std::string> ref_bases;      // REF, it's raw REF alleles  
    std::vector<std::string> ale_bases;      // ALT, it's REF and non-REF alleles
    std::vector<int> depths;                 // depth for each type of base
    std::vector<double> freqs;               // frequency for each type of base
    std::vector<StrandBiasInfo> strand_bias; // strand bias for each type of base

    VariantInfo() : ref_id(""), ref_pos(0), total_depth(0), qual(0) {};
    VariantInfo(const std::string& rid, uint32_t pos, int dp, double qs) 
        : ref_id(rid), ref_pos(pos), total_depth(dp), qual(qs) {};
    
    // Helper method to format record as VCF string
    std::string to_string() const {
        std::stringstream ss;
        ss << ref_id                       << "\t"
           << ref_pos                      << "\t"
           << total_depth                  << "\t"
           << qual                         << "\t"
           << ngslib::join(ref_bases, ",") << "\t"
           << ngslib::join(ale_bases, ",") << "\t"
           << ngslib::join(depths, ",")    << "\t"
           << ngslib::join(freqs, ",");

        return ss.str();
    };
};

struct AlleleInfo {
    std::string ref;
    std::vector<std::string> alts;  // uniq ALT alleles string
    
    int total_alleles = 0;                                  // AN, total alleles, including REF and ALT
    std::map<std::string, int> allele_counts;               // AC, allele count for each allele, key is the allele string, value is the count
    std::map<std::string, double> allele_freqs;             // AF, allele frequency, key is the allele string, value is the frequency

    std::map<std::string, int> allele_depths;               // AD, depth for each allele, key is the allele string, value is the depth
    std::map<std::string, StrandBiasInfo> strand_bias_info; // DP4, strand bias info for each allele
    int total_dp = 0;                                       // DP, total depth of all alleles, including REF and ALT
};

struct VCFSampleAnnotation {
    int GQ;                               // Genotype Quality
    std::vector<size_t> gtcode;           // Genotype
    std::vector<std::string> sample_alts; 
    std::vector<int> allele_depths;       // AD, depth
    std::vector<int> PL;                  // PL: A list of phred-scale score of genotype likelihoods
};

struct VCFRecord {
    // Required fields
    std::string chrom;     // CHROM: chromosome name
    uint32_t pos;          // POS: 1-based position
    std::string id;        // ID: variant identifier
    std::string ref;       // REF: reference allele
    std::vector<std::string> alt;  // ALT: alternate alleles
    int qual;              // QUAL: quality score
    std::string filter;    // FILTER: filter status
    std::string info;      // INFO: additional information
    
    std::string format;    // FORMAT: format string for genotype fields
    std::vector<std::string> samples;  // Sample information.

    // Constructor
    VCFRecord() : chrom(""), pos(0), id("."), ref(""), qual(0), filter("."), info(""), 
                  format("") {}; 
    
    // Helper method to validate VCF record
    bool is_valid() const {
        // Basic validation
        if (chrom.empty() || pos == 0 || ref.empty() || alt.empty()) {
            return false;
        }

        // REF must be A,C,G,T,N or * for structural variants
        for (char c : ref) {
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && 
                c != 'N' && c != '*' && 
                c != 'a' && c != 'c' && c != 'g' && c != 't' && 
                c != 'n') {
                return false;
            }
        }

        // ALT validation
        for (const auto& a : alt) {
            if (a == ref) {
                return false;  // ALT should not be same as REF
            }
        }

        return true;
    }

    // Helper method to format record as VCF string
    std::string to_string(const std::string& sep = "\t") const {
        std::stringstream ss;
        
        // Required fields
        ss << chrom << sep
           << pos   << sep
           << id    << sep
           << ref   << sep
           << (alt.empty() ? "." : ngslib::join(alt, ",")) << sep
           << (qual <= 0 ? "." : std::to_string(qual))     << sep
           << filter << sep
           << info;

        // Optional fields
        if (!samples.empty()) {
            ss << sep << format;
            for (const auto& sample : samples) {
                ss << sep << sample;
            }
        }

        return ss.str();
    }
};

// 默认保留小数点后 3 位
std::string format_double(double value, int precision = 3);

/**
 * @brief calculate the strand bias for a reference and alternative base
 * 
 * @param ref_base 
 * @param alt_bases_string 
 * @param bases 
 * @param strands 
 * @return StrandBiasInfo 
 */
StrandBiasInfo strand_bias(const std::string &ref_base, 
                           const std::string &alt_bases_string,
                           const std::vector<std::string> &bases,
                           const std::vector<char> &strands);


// Collect and normalized REF/ALT information
AlleleInfo collect_and_normalized_allele_info(VariantInfo &variant, std::vector<BaseType::BatchInfo> &smps_bi_v);

// 单个样本的信息 
VCFSampleAnnotation process_sample_variant(const std::string& upper_ref_base, const std::vector<std::string>& alts, 
                                           const BaseType::BatchInfo& smp_bi);
// Helper function: Convert PL index to genotype pair
std::pair<size_t, size_t> pl_index_to_genotype(size_t pl_idx, size_t n_alleles);


std::string format_sample_string(const VCFSampleAnnotation& sa);

/// Header for VCF
std::string vcf_header_define(const std::string &ref_file_path, const std::vector<std::string> &addition_info, 
                              const std::vector<std::string> &samples, const std::string other_comment);
std::string cvg_header_define(const std::vector<std::string> &group_info, const std::vector<std::string> &BASES);
void merge_file_by_line(const std::vector<std::string> &infiles, const std::string &outfile, 
                        std::string header="#", bool is_remove_tempfile=false);

#endif