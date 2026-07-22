/**
 * @file basetype_utils.h
 * 
 * @brief define the global variable and util codes specific for basetype.h
 *  
 * @author Shujia Huang
 * @date 2018-08-29
 * 
 */
#ifndef __INCLUDE_BASEVAR_CALLER_UTILS_H__
#define __INCLUDE_BASEVAR_CALLER_UTILS_H__

#include <string>
#include <vector>
#include <iostream>

#include "io/utils.h"
#include "basetype.h"
#include "external/robin_hood.h"  // robin_hood::unordered_map, robin_hood::unordered_set

struct AlignBase {
    std::string ref_base;   // reference base
    std::string read_base;  // read base
    char base_qual;         // read base quality, assign mean quality of seq if it's Indels 
    int rpr;                // read position rank, 1-based, the first base in read has rpr = 1, the second base has rpr = 2, etc.

    int mapq;               // mapping quality
    char map_strand;        // mapping reference strand, should be one of '*', '-', '+'  

    std::string to_string() const {
        std::string s;
        // ref_base/read_base 通常很短；32 足够容纳两个整数、两个 char 和 5 个制表符
        s.reserve(ref_base.size() + read_base.size() + 32);
        s.append(ref_base);  s.push_back('\t');
        s.append(read_base); s.push_back('\t');
        
        s.push_back(base_qual);  // char：按字符输出，与原 << 语义一致
        s.push_back('\t');

        s.append(std::to_string(rpr));  s.push_back('\t');
        s.append(std::to_string(mapq)); s.push_back('\t');
        s.push_back(map_strand);  // char：按字符输出

        return s;
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
        // freqs 是 vector<double>，走 ngslib::join 泛型模板以保持原有小数格式
        const std::string pos_s   = std::to_string(ref_pos);
        const std::string dp_s    = std::to_string(total_depth);
        const std::string qual_s  = std::to_string(qual);
        const std::string ref_s   = ngslib::join(ref_bases, ",");
        const std::string alt_s   = ngslib::join(ale_bases, ",");
        const std::string depth_s = ngslib::join(depths, ",");
        const std::string freq_s  = ngslib::join(freqs, ",");

        std::string s;
        s.reserve(ref_id.size() + pos_s.size() + dp_s.size() + qual_s.size() +
                  ref_s.size() + alt_s.size() + depth_s.size() + freq_s.size() + 7);
        s.append(ref_id);   s.push_back('\t');
        s.append(pos_s);    s.push_back('\t');
        s.append(dp_s);     s.push_back('\t');
        s.append(qual_s);   s.push_back('\t');
        s.append(ref_s);    s.push_back('\t');
        s.append(alt_s);    s.push_back('\t');
        s.append(depth_s);  s.push_back('\t');
        s.append(freq_s);

        return s;
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

    // Dosage-based counts (expected from genotype posteriors)
    std::vector<double> dosage_counts;                      // Per-ALT expected AC: [dosage_ALT1, dosage_ALT2, ...]
    int dosage_total_alleles = 0;                           // = 2 * N_samples

    // GT-based counts (from posterior GT, matches VCF GT column)
    std::vector<int> gt_counts;                             // Per-ALT count from posterior GT: [GT_ALT1, GT_ALT2, ...]
    int gt_total_alleles = 0;                               // Total alleles from non-missing posterior GT calls
};

struct VCFSampleAnnotation {
    double GQ;                             // Genotype Quality (posterior-based when population prior applied)
    std::vector<size_t> gtcode;            // Genotype
    std::vector<std::string> sample_alts;  
    std::vector<int> allele_depths;        // AD, depth
    std::vector<int> PL;                   // PL: A list of phred-scale score of genotype likelihoods (unchanged)
    // --- Bayesian posterior fields ---
    std::vector<double> posterior;         // Genotype posterior probabilities
    double dosage = 0.0;                   // Expected ALT allele dosage (total)
    std::vector<double> per_allele_dosage; // Per-ALT expected allele count: [dosage_ALT1, dosage_ALT2, ...]
};

struct VCFTextLine {
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
    VCFTextLine() : chrom(""), pos(0), id("."), ref(""), qual(0), filter("."), info(""), 
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
        const std::string pos_s  = std::to_string(pos);
        const std::string alt_s  = alt.empty() ? "." : ngslib::join(alt, ",");
        const std::string qual_s = qual <= 0 ? "." : std::to_string(qual);

        // Pre-calculate total size for a single allocation
        size_t est = chrom.size() + pos_s.size() + id.size() + ref.size() + alt_s.size() + 
                     qual_s.size() + filter.size() + info.size() + 7 * sep.size();

        if (!samples.empty()) {
            est += sep.size() + format.size();
            for (const auto& sample : samples)
                est += sep.size() + sample.size();
        }

        std::string s;
        s.reserve(est);

        // Required fields
        s.append(chrom);  s.append(sep);
        s.append(pos_s);  s.append(sep);
        s.append(id);     s.append(sep);
        s.append(ref);    s.append(sep);
        s.append(alt_s);  s.append(sep);
        s.append(qual_s); s.append(sep);
        s.append(filter); s.append(sep);
        s.append(info);

        // Optional fields
        if (!samples.empty()) {
            s.append(sep);
            s.append(format);
            for (const auto& sample : samples) {
                s.append(sep);
                s.append(sample);
            }
        }

        return s;
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
                                           const BaseType::BatchInfo& smp_bi, double af = -1.0,
                                           double ref_bias = 0.5);
// Helper function: Convert PL index to genotype pair
std::pair<size_t, size_t> pl_index_to_genotype(size_t pl_idx, size_t n_alleles);

std::string format_sample_string(const VCFSampleAnnotation& sa);

/// Header for VCF
std::string vcf_header_define(const std::string &ref_file_path, const std::vector<std::string> &addition_info, 
                              const std::vector<std::string> &samples, const std::string other_comment);
void merge_file_by_line(const std::vector<std::string> &infiles, const std::string &outfile, 
                        std::string header="#", bool is_remove_tempfile=false);

#endif