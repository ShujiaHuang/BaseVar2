/**
 * @file algorthm.h
 * 
 * 
 * @brief  Contain some main algorithms of BaseVar
 *  
 * @author Shujia Huang
 * @date 2018-08-30
 * 
 */
#ifndef __INCLUDE_BASRVAR_ALIGORITHM_H__
#define __INCLUDE_BASRVAR_ALIGORITHM_H__

#include <algorithm>
#include <vector>
#include <string>
#include <cmath>  // use 'log' functon
#include <numeric>
#include <utility>
#include <stdexcept>

#include <htslib/kfunc.h>

template<class ForwardIterator>
inline size_t argmin(ForwardIterator first, ForwardIterator last) {
    return std::distance(first, std::min_element(first, last));
}

template<class ForwardIterator>
inline size_t argmax(ForwardIterator first, ForwardIterator last) {
    return std::distance(first, std::max_element(first, last));
}

// sum the value for all the data which could call '+' operator
template<typename T> // T must be a numeric type
typename std::enable_if<std::is_arithmetic<T>::value, T>::type
sum(const std::vector<T> &values) {
    if (values.empty()) {
        return T(0);
    }
    return std::accumulate(values.begin(), values.end(), T(0));
}

// mean the value for all the data which could call '+' operator
template<typename T>  // T must be a numeric type
double mean(const std::vector<T> &values) {
    if (values.empty()) {
        throw std::invalid_argument("Cannot calculate mean of empty vector");
    }
    return static_cast<double>(sum(values)) / values.size();
}

// standard deviation the value for all the data which could call '+' operator
template<typename T>  // T must be a numeric type
double stddev(const std::vector<T> &values) {
    if (values.size() < 2) {
        throw std::invalid_argument("Cannot calculate standard deviation of less than 2 values");
    }
    double mean_value = mean(values);
    double sum_squared_diff = 0.0;
    for (const auto& value : values) {
        sum_squared_diff += (value - mean_value) * (value - mean_value);
    }
    return std::sqrt(sum_squared_diff / (values.size() - 1));
}

// calculate the standard error of the mean
template<typename T>  // T must be a numeric type
double standard_error(const std::vector<T> &values) {
    if (values.size() < 2) {
        throw std::invalid_argument("Cannot calculate standard error of less than 2 values");
    }
    return stddev(values) / std::sqrt(values.size());
}

// median the value for all the data which could call '+' operator
template<typename T>  // T must be a numeric type
double median(std::vector<T> &value) {
    size_t n = value.size();
    if (n == 0) return 0.0;

    std::sort(value.begin(), value.end());
    if (n % 2 == 0) {
        return (value[n/2-1] + value[n/2]) / 2.0;
    } else {
        return value[n/2];
    }
}

// Calculate the Genotype likelihoods for a given base and quality against the reference base
// ref_bias: reference bias coefficient β. For het (REF, ALT), the genotype likelihood is:
//   P(R|G) = (1-β) * P(R|REF) + β * P(R|ALT)
// Default β=0.5 means no reference bias (standard diploid model).
// Typical empirical values: 0.45–0.48 (β < 0.5 reflects alignment reference bias).
std::vector<int> calculatePL(const std::string& ref_base, 
                             const std::vector<std::string>& alt_bases, 
                             const std::vector<std::string>& align_bases, 
                             const std::vector<char>& base_quals,
                             double ref_bias = 0.5);

// Function for chi^2 test
double chi2_test(double chi_sqrt_value, double degree_of_freedom);
double norm_dist(double x);

enum class TestSide {
    LESS,
    GREATER,
    TWO_SIDED
};
/**
 *  Perform a Fisher exact test on a 2x2 contingency table. 
 *  The null hypothesis is that the true odds ratio of the 
 *  populations underlying the observations is one.
 *  
 *  Mathmatical note: https://mathworld.wolfram.com/FishersExactTest.html
 * 
 *  @param n11 Value in cell (1,1)
 *  @param n12 Value in cell (1,2)
 *  @param n21 Value in cell (2,1)
 *  @param n22 Value in cell (2,2)
 *  @param test_side Type of test to perform (left-sided, right-sided, or two-sided)
 *  @return p-value for the specified test
 * 
 *    n11  n12  | n1_
 *    n21  n22  | n2_
 *   -----------+----
 *    n_1  n_2  | n
 * 
 *  Example: https://gatk.broadinstitute.org/hc/en-us/articles/360035532152-Fisher-s-Exact-Test
 * 
 */
double fisher_exact_test(int n11, int n12, int n21, int n22, TestSide test_side=TestSide::TWO_SIDED);
struct ContingencyTable {
    int n11, n12, n21, n22;
    ContingencyTable(int a, int b, int c, int d) : n11(a), n12(b), n21(c), n22(d) {
        if (n11 < 0 || n12 < 0 || n21 < 0 || n22 < 0) {
            throw std::invalid_argument("Contingency table entries must be non-negative");
        }
    }
};
double fisher_exact_test(const ContingencyTable& table, TestSide test_side=TestSide::TWO_SIDED);

double wilcoxon_ranksum_test(const std::vector<double>& sample1, const std::vector<double>& sample2);

/**
 * @brief Wilcoxon rank sum test returning Z-score (for GATK-compatible RankSum annotations).
 * 
 * Positive Z indicates sample1 values tend to be higher; negative indicates lower.
 * 
 * @param sample1 First sample (typically REF reads)
 * @param sample2 Second sample (typically ALT reads)
 * @return Z-score (0.0 if either sample is empty)
 */
double wilcoxon_ranksum_zscore(const std::vector<double>& sample1, const std::vector<double>& sample2);

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
            std::vector<std::vector<double>> &ind_allele_post_prob,  // return value, update value inplace 
            std::vector<double> &marginal_likelihood);               // return value, calculate value inplace

/**
 * @brief Update observed allele frequency inplace.
 * 
 * @param ind_allele_post_prob 2d-array, n x 4 matrix.
 * @param obs_allele_freq      1d-array, 1 x 4. update the allele frequence inplace and return. 
 * 
 */
void m_step(const std::vector<std::vector<double>> &ind_allele_post_prob, std::vector<double> &obs_allele_freq); 

/**
 * @brief EM algorithm
 * 
 * @param ind_allele_likelihood  n x 4 matrix, n is non-N and non-Indel's sample size (same below), 4 for [A, C, G, T]
 * @param obs_allele_freq        1 x 4 vector, Observed allele frequence for [A, C, G, T]
 * @param iter_num integer, optional.  The lager EM iteration times. default: 100
 * @param epsilon  float, optional. The threshold of likelihood different for EM 
 *                 process. default 0.001.
 * 
 */
void EM(const std::vector<std::vector<double>> &ind_allele_likelihood, // n x 4 matrix, do not change the raw value.
        std::vector<double> &obs_allele_freq,          // retuen value, 1 x 4, expect allele frequence, it'll be update inplace here.
        std::vector<double> &log_marginal_likelihood,  // return value
        int iter_num=100, const float epsilon=0.001);

// ============================================================================
// Genotype posterior computation (Bayesian genotype calling)
// ============================================================================

/**
 * @brief Convert Phred-scaled PL values to normalized genotype likelihoods.
 * 
 * L(g) = 10^(-PL(g)/10), then normalize so sum(L) = 1.
 * Numerically stable: subtracts min(PL) before conversion.
 * 
 * @param PL Phred-scaled genotype likelihoods
 * @return Normalized likelihoods summing to 1
 */
std::vector<double> pl_to_likelihoods(const std::vector<int>& PL);

/**
 * @brief Compute Hardy-Weinberg genotype prior probabilities.
 * 
 * For n alleles with frequencies f_0, f_1, ..., f_{n-1}:
 *   P(i/j) = (2 - delta_ij) * f_i * f_j
 * 
 * Allele frequencies are clamped to [1e-6, 1-1e-6] to avoid degenerate priors.
 * 
 * @param n_alleles Number of alleles (including REF)
 * @param allele_freqs Allele frequency vector (length = n_alleles)
 * @return Genotype prior probabilities (length = n_alleles*(n_alleles+1)/2)
 */
std::vector<double> hw_genotype_prior(size_t n_alleles, const std::vector<double>& allele_freqs);

/**
 * @brief Dosage-based HWE chi-square test (suitable for low-coverage data).
 * 
 * Uses per-sample genotype probabilities ("soft" genotype counts) instead of
 * hard GT calls, making it robust for ultra-low-coverage sequencing.
 * 
 * @param genotype_probs Per-sample genotype probability vectors (length = n_genotypes each)
 * @param n_alleles Number of alleles (including REF), typically 2 for bi-allelic
 * @return HWE chi-square p-value (1.0 if test is not applicable)
 */
double hwe_dosage_test(const std::vector<std::vector<double>>& genotype_probs, size_t n_alleles);

/**
 * @brief Structure holding genotype posterior computation results.
 */
struct GenotypePosterior {
    size_t best_gt_idx;                  // Index of genotype with highest posterior
    std::vector<double> posteriors;      // Posterior probability for each genotype
    double dosage;                       // Total expected ALT allele count
    std::vector<double> per_allele_dosage; // Per-ALT expected allele count: [dosage_ALT1, dosage_ALT2, ...]
    double gq;                           // GQ = -10 * log10(1 - P(best_gt))
};

/**
 * @brief Compute genotype posterior probabilities using PL and per-allele population AF prior.
 * 
 * Combines data likelihood (from PL) with Hardy-Weinberg prior (from allele frequencies):
 *   P(g|D,f) propto P(D|g) * P(g|f)
 * 
 * Supports multi-allelic sites by accepting per-allele frequency vector.
 * 
 * @param PL Phred-scaled genotype likelihoods (VCF PL ordering)
 * @param allele_freqs Per-allele frequencies: [freq_REF, freq_ALT1, freq_ALT2, ...]
 * @return GenotypePosterior with best GT, posteriors, per-allele dosage, and GQ
 */
GenotypePosterior compute_genotype_posterior(
    const std::vector<int>& PL,
    const std::vector<double>& allele_freqs
);

/**
 * @brief Overload for backward compatibility: single ALT allele frequency.
 * 
 * @param PL Phred-scaled genotype likelihoods
 * @param af Population ALT allele frequency (bi-allelic: freq_ALT; REF = 1-af)
 * @return GenotypePosterior with best GT, posteriors, dosage, and GQ
 */
GenotypePosterior compute_genotype_posterior(
    const std::vector<int>& PL,
    double af
);

/**
 * @brief Compute dosage-based allele counts from genotype posteriors (multi-allelic).
 * 
 * @param sample_posts Genotype posteriors for all samples
 * @return {expected_ac_per_alt, expected_an} where
 *         expected_ac_per_alt[k] = sum of per_allele_dosage[k] for all samples,
 *         expected_an = 2 * N
 */
std::pair<std::vector<double>, int> compute_dosage_ac(
    const std::vector<GenotypePosterior>& sample_posts
);

#endif
