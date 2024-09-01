/**
 * @file basetype.cpp
 * @brief 
 * 
 * @author Shujia Huang
 * @date 2018-08-01
 * 
 */

#include <numeric>  // use std::accumulate function

#include "basetype.h"
#include "algorithm.h"
#include "external/combinations.h"

//////////////////////////////////////////////////////////////
//// The codes for the member function of BaseType class /////
//////////////////////////////////////////////////////////////
BaseType::BaseType(const BatchInfo *smp_bi, double min_af) : _only_call_snp(true) {

    _min_af = min_af;
    _total_depth = 0;

    // some common value
    for (size_t i(0); i < BASES.size(); ++i) {
        _B_IDX[BASES[i]] = i;  // set [A, C, G, T] char map to array index
        _depth[BASES[i]] = 0;  // inital the base depth to 0.
    }

    // Initialized the array to {0, 0, 0, 0}, which set allele likelihood for [A, C, G, T]
    std::vector<double> allele_lh(BASES.size(), 0);
    _ind_allele_likelihood.reserve(smp_bi->n);

    double epsilon; // base error probability
    char fb;
    for (size_t i(0); i < smp_bi->n; ++i) {

        fb = smp_bi->align_bases[i][0];  // a string, get first base
        if (fb != 'N' && ((fb != '+' && fb != '-') || !_only_call_snp)) { 
            // ignore all the 'N' bases and indels if only call snp

            if (smp_bi->align_bases[i].size() != 1)
                throw std::runtime_error("[ERROR] Why dose the size of aligned base is not 1? Check: " + 
                                         smp_bi->align_bases[i]); 
            
            _depth[fb]++;  // record the depth for read base: [A, C, G, T]
            _total_depth++;

            for(size_t j(0); j < BASES.size(); ++j) {
                // convert the quality phred scale to be the base confident probabilty value
                epsilon = exp((smp_bi->align_base_quals[i] - 33) * MLN10TO10);
                allele_lh[j] = (fb == BASES[j]) ? 1.0 - epsilon : epsilon / 3;
            }
            _ind_allele_likelihood.push_back(allele_lh); // A 2d-array, n x 4 matrix, n is sample size.
        }
    }
}

std::vector<double> BaseType::_set_allele_frequence(const std::vector<char> &bases) {
    // bases 数组中 A,C,G,T 这四个碱基最多只能各出现一次
    
    double total_depth = 0;
    for (auto b: bases) total_depth += this->_depth[b];

    // Initialized the array to {0, 0, 0, 0}, which set initial observed allele likelihood for [A, C, G, T]
    std::vector<double> obs_allele_freq(BASES.size(), 0);
    if (total_depth > 0) {
        // computed by base count
        for (auto b: bases) obs_allele_freq[_B_IDX[b]] = this->_depth[b] / total_depth;
    }

    return obs_allele_freq;  // 1 x 4 vector. The allele frequence of [A, C, G, T]
}

/**
 * @brief Calculate population likelihood for all the combination of bases
 * 
 * @param bases A 1-d array. An array subset of bases from [A, C, G, T] 
 * @param n     The combination number. n must less or equal to the length of ``bases``
 * 
 * @return AA   AA.bc: An array of combinamtion bases
 *              AA.lr: Likelihood of ``bc``
 * 
 */
AA BaseType::_f(const std::vector<char> &bases, int n) {

    AA data;

    Combinations<char> c(bases, n);
    std::vector<std::vector<char>> cbs_v = c.get();  // combination bases vector
    for (size_t i = 0; i < cbs_v.size(); i++) {

        std::vector<double> obs_allele_freq = this->_set_allele_frequence(cbs_v[i]);
        double sum_freq = std::accumulate(obs_allele_freq.begin(), obs_allele_freq.end(), 0);
        if (sum_freq == 0) // Empty coverage for this type of combination, skip.
            continue;
        
        std::vector<double> log10_marginal_likelihood;
        EM(_ind_allele_likelihood, obs_allele_freq, log10_marginal_likelihood);

        double sum_log10_marginal_likelihood = std::accumulate(log10_marginal_likelihood.begin(), 
                                                               log10_marginal_likelihood.end(), 0);
        
        data.bc.push_back(cbs_v[i]);
        data.lr.push_back(sum_log10_marginal_likelihood);
        data.bp.push_back(obs_allele_freq);
    }
    
    return data;
}

// The main function for likelihood ratio test
void BaseType::lrt(const std::vector<char> &specific_bases) {
    
    std::vector<char> active_bases;
    // Get active bases which count frequence > _min_af
    for (auto b: specific_bases) {
        if (_depth[b] / _total_depth >= _min_af)
            active_bases.push_back(b);
    }

    if (active_bases.size() == 0) return;

    // init. Base combination of active_bases
    AA var = _f(active_bases, active_bases.size());

    double chi_sqrt_value = 0;
    std::vector<double> base_frq = var.bp[0];
    std::vector<double> lr_alt   = var.lr;
    std::vector<double> ltr_chivalue;

    // Find candinate altnative alleles
    for (size_t n = active_bases.size() - 1; n > 0; --n) {
        var = _f(active_bases, n);

        // ltr_chivalue = ;
    }

    return;
}