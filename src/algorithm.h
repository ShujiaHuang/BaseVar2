/**
 * @file algorthm.h
 * 
 * 
 * @brief
 *  
 * @author Shujia Huang
 * @date 2024-08-30
 * 
 */

#ifndef __INCLUDE_BASRVAR_ALIGORITHM_H__
#define __INCLUDE_BASRVAR_ALIGORITHM_H__

#include <iostream>
#include <vector>
#include <cmath>  // use 'exp' functon


/**
 * @brief Calculate the posterior probability of individual allele at each site as
 *        the four A/C/G/T bases.
 * 
 * @param obs_allele_freq        1 x 4 matrix.
 * @param ind_allele_likelihood  n x 4 matrix. n is sample size
 * @param ind_allele_post_prob   n x 4 matrix. allele posterior probabilty for each sample.
 * @param marginal_likelihood    n x 1 matrix. maginal likelihood for each sample.
 * 
 */

void e_step(const std::vector<double> &obs_allele_freq,
            const std::vector<std::vector<double>> &ind_allele_likelihood, 

            std::vector<std::vector<double>> &ind_allele_post_prob,  // return value, calculate this value in this function  
            std::vector<double> &marginal_likelihood)                // return value, calculate this value in this function 
{
    // 'likelihood' is as the same shape as `ind_allele_likelihood`
    std::vector<std::vector<double>> likelihood(ind_allele_likelihood.size(), 
                                                std::vector<double>(obs_allele_freq.size(), 0));

    size_t sample_size = ind_allele_likelihood.size(), allele_num = obs_allele_freq.size();
    for (size_t i(0); i < sample_size; ++i) {
        
        marginal_likelihood[i] = 0; // reset the raw value to be 0
        for (size_t j = 0; j < allele_num; j++) {  // for [A, C, G, T]
            likelihood[i][j] = ind_allele_likelihood[i][j] * obs_allele_freq[j];
            marginal_likelihood[i] += likelihood[i][j];
        }

        // Computed the posterior probability of A/C/G/T for each individual, change the value inplace.
        for (size_t j(0); j < allele_num; ++j) {
            // reset the posterior value
            ind_allele_post_prob[i][j] = likelihood[i][j] / marginal_likelihood[i];
        }
    }

    return;
}

/**
 * @brief Update observed allele frequency inplace.
 * 
 * @param ind_allele_post_prob 2d-array, n x 4 matrix.
 * @param obs_allele_freq      1d-array, 1 x 4. update the allele frequence inplace and return. 
 * 
 */
void m_step(const std::vector<std::vector<double>> &ind_allele_post_prob, std::vector<double> &obs_allele_freq) {

    size_t sample_size = ind_allele_post_prob.size();
    size_t allele_num  = ind_allele_post_prob[0].size();

    // Reset data
    obs_allele_freq = std::vector<double>(allele_num, 0);
    for (size_t i(0); i < sample_size; ++i) {
        for (size_t j(0); j < allele_num; ++j)
            obs_allele_freq[j] += ind_allele_post_prob[i][j];
    }

    for (size_t i(0); i < allele_num; ++i) 
        obs_allele_freq[i] /= sample_size;

    return;
}

/**
 * @brief EM algorithm
 * 
 * @param obs_allele_freq        1 x 4 vector, Observed allele frequence for [A, C, G, T]
 * @param ind_allele_likelihood  n x 4 matrix, n is sample size, 4 for [A, C, G, T]
 * @param iter_num integer, optional.  The lager EM iteration times. default: 100
 * @param epsilon  float, optional. The threshold of likelihood different for EM 
 *                 process. default 0.001.
 * 
 */
void EM(const std::vector<std::vector<double>> &ind_allele_likelihood, // n x 4 matrix, do not change the raw value.
        std::vector<double> &obs_allele_freq,   // retuen value, 1 x 4, expect allele frequence, it'll be update inplace here.
        std::vector<double> &log10_marginal_likelihood, // return value
        int iter_num=100, const float epsilon=0.001) 
{
    // n x 4 matrix, the same shape as 'ind_allele_likelihood'.
    std::vector<std::vector<double>> ind_allele_post_prob = std::vector<std::vector<double>>(
        ind_allele_likelihood.size(), std::vector<double>(obs_allele_freq.size(), 0));

    // It's a 1-d array (n x 1) one sample per value, n is sample size
    std::vector<double> marginal_likelihood = std::vector<double>(ind_allele_likelihood.size(), 0);

    // Update the value of 'ind_allele_post_prob' and compute 'marginal_likelihood' in e_step.
    e_step(obs_allele_freq, ind_allele_likelihood, ind_allele_post_prob, marginal_likelihood);

    log10_marginal_likelihood = std::vector<double>(marginal_likelihood.size(), 0);
    for (size_t i = 0; i < marginal_likelihood.size(); i++) {
        log10_marginal_likelihood[i] = log10(marginal_likelihood[i]);
    }
    
    // use 'ind_allele_post_prob' to update 'obs_allele_freq'
    m_step(ind_allele_post_prob, obs_allele_freq);

    while (iter_num--) {
        // e step: update ind_allele_post_prob
        e_step(obs_allele_freq, ind_allele_likelihood, ind_allele_post_prob, marginal_likelihood);

        // m step: update the frequence of observed alleles
        m_step(ind_allele_post_prob, obs_allele_freq); 

        double delta = 0, llh = 0;
        for (size_t i = 0; i < marginal_likelihood.size(); i++) {
            llh = log10(marginal_likelihood[i]);
            delta += abs(llh - log10_marginal_likelihood[i]);
            log10_marginal_likelihood[i] = llh;  // update
        }
        
        // Todo: too big? be careful here!!!
        if (delta < epsilon) break;
    }

    return;
}

#endif
