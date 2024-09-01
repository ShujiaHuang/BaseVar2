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
#include <cmath>  // use exp() function
#include <string>
#include <vector>
#include <map>

static const std::vector<char> BASES = {'A', 'C', 'G', 'T'}; // 定义这个值，限定 likelihood 数组中碱基似然值的存放顺序
static const double MLN10TO10   = -0.23025850929940458;      // ln(10)/10，把 phred-value 换成 e 为底，方便调用 exp()
static const int LRT_THRESHOLD  = 24;                        // 24 corresponding to a chi-pvalue of 10^-6
static const int QUAL_THRESHOLD = 60;                        // -10 * lg(10^-6)

// Mainly use for basevar variant callling
struct BatchInfo {

    size_t n;

    std::string ref_id, ref_base;
    uint32_t ref_pos;
    uint32_t depth; // The coverage depth on ref_pos

    std::vector<std::string> align_bases;
    std::vector<char> align_base_quals;
    // std::vector<double> align_base_quals_pvalue;  // quali phred-scale qual(int) -> base probabilty value

    std::vector<int> mapqs;
    std::vector<char> map_strands;
    std::vector<int> base_pos_ranks;

    // Set default argument
    BatchInfo(): n(0), ref_pos(0), depth(0) {}
};

/**
 * @brief Define a structure for recording allele information return by _f() 
 * in this class.
 * 
 */
typedef struct {
    std::vector<std::vector<char>> bc;
    std::vector<std::vector<double>> bp;
    std::vector<double> lr;
} AA;

// A class for calculate the base probability
class BaseType {

private:
    bool _only_call_snp;

    std::string _ref_id, _ref_base;
    uint32_t ref_pos;
    std::vector<std::string> _alt_bases;    
    double _min_af;

    
    // [A, C, G, T] likelihood vector for echo individual
    std::vector<std::vector<double>> _ind_allele_likelihood; // 2d-array, n x 4 matrix, n is sample size.
    std::map<char, size_t> _B_IDX;  // A map for recrod the BASE => index
    
    // Estimated allele frequency by EM and LRT
    std::map<std::string, double> _af_by_lrt;

    int _total_depth;              // sum depth of ACGT
    std::map<char, double> _depth; // allele depth of [A, C, G, T]

    // init the base likelihood by input bases
    std::vector<double> _set_allele_frequence(const std::vector<char> &bases);

    // Calculate population likelihood for all the combination of bases
    AA _f(const std::vector<char> &bases, int n);

    BaseType(const BaseType &b) = delete;             // reject using copy constructor (C++11 style).
    BaseType &operator=(const BaseType &b) = delete;  // reject using copy/assignment operator (C++11 style).

public:
    // Constructor
    BaseType(const BatchInfo *smp_bi, double af);
    ~BaseType() {};

    // The main function for likelihood ratio test
    /**
     * @brief 
     * 
     * @param specific_bases a 1d-array. [optional]
     * Calculating LRT for specific base combination if provided.
     */
    void lrt() { /* default */ lrt(BASES); }
    void lrt(const std::vector<char> &specific_bases);

}; // BaseType class

#endif
