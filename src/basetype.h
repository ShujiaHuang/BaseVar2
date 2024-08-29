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
#include <map>

#include <cmath>  // use exp() function

static const std::string BASES = "ACGT";  // 定义这个值，目的是为了限定 likehood 数组中碱基似然值的存放顺序
static const std::map<char, int> B_IDX = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};
static const double MLN10TO10 = -0.23025850929940458;  // 换底，把 10 换为 e 为底：log(10)/10;


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

// A class for calculate the base probability
class BaseType {

private:
    bool _only_call_snp;

    std::string _ref_id, _ref_base;
    uint32_t ref_pos;
    std::vector<std::string> _alt_bases;
    double _min_af;

    //Allele [A, C, G, T] likelihood vector for echo individual
    std::vector<std::vector<double>> ind_allele_likelihood;

    BaseType(const BaseType &b) = delete;             // reject using copy constructor (C++11 style).
    BaseType &operator=(const BaseType &b) = delete;  // reject using copy/assignment operator (C++11 style).

public:
    // Constructor
    BaseType(const BatchInfo *smp_bi, double af);

    ~BaseType() {_min_af = 0;};
    void set_minaf(double af) { _min_af = af; }

    void lrt();

}; // BaseType class

#endif
