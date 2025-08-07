#include "algorithm.h"

// Function for chi^2 test
double chi2_test(double chi_sqrt_value, double degree_of_freedom) {
    return kf_gammaq(degree_of_freedom/2.0, chi_sqrt_value/2.0);
}

double norm_dist(double x) {
    return kf_erfc(double(x / std::sqrt(2.0))) / 2.0;
}

std::vector<int> calculatePL(const std::string& ref_base, 
                             const std::vector<std::string>& alt_bases, 
                             const std::vector<std::string>& align_bases, 
                             const std::vector<char>& base_quals) 
{
    // This function calculates the Genotype likelihoods (PL) for given bases and quality against the reference base.
    if (align_bases.size() != base_quals.size()) {
        throw std::invalid_argument("Bases and base quality vectors must be of the same length");
    }

    // Calculate number of alleles (REF + ALTs)
    size_t n_alleles = alt_bases.size() + 1;
    
    // Calculate number of possible genotypes: n(n+1)/2 where n is number of alleles
    size_t n_genotypes = (n_alleles * (n_alleles + 1)) / 2;
    
    // Create all possible genotype combinations
    std::vector<std::pair<std::string, std::string>> genotypes;
    genotypes.reserve(n_genotypes);
    
    // Generate all possible genotype pairs (including homozygous and heterozygous)
    for (size_t i = 0; i < n_alleles; i++) {
        for (size_t j = i; j < n_alleles; j++) {
            std::string allele1 = (i == 0) ? ref_base : alt_bases[i-1];
            std::string allele2 = (j == 0) ? ref_base : alt_bases[j-1];
            genotypes.push_back({allele1, allele2});
        }
    }

    // Initialize log-likelihoods for each genotype
    std::vector<double> log10_L(n_genotypes, 0.0);

    // Iterate over each read in the pileup
    for (size_t i = 0; i < align_bases.size(); ++i) {
        if (align_bases[i][0] == 'N' || align_bases[i][0] == 'n') continue; // Skip Ns

        std::string B = align_bases[i];
        int Q = static_cast<int>(base_quals[i]) - 33;
        
        double e = (Q < 0) ? 1.0 : std::pow(10.0, -Q / 10.0);
        double base_quality_prob = 1.0 - e;

        double P_correct = base_quality_prob;
        double P_error   = e / 3.0; // Probability of error (uniform across other bases)

        // For each genotype, compute P(R|G)
        for (size_t g = 0; g < n_genotypes; ++g) {
            const auto& [allele1, allele2] = genotypes[g];
            
            // Calculate P(R|allele) for both alleles in the genotype
            double P_R_given_allele1 = (B == allele1) ? P_correct : P_error;
            double P_R_given_allele2 = (B == allele2) ? P_correct : P_error;

            // P(R|G) = 0.5 * P(R|allele1) + 0.5 * P(R|allele2)
            double P_R_given_G = 0.5 * P_R_given_allele1 + 0.5 * P_R_given_allele2;

            // Add log10(P(R|G)) to the cumulative log-likelihood
            if (P_R_given_G > 0) {
                log10_L[g] += log10(P_R_given_G);
            } else {
                // Handle potential underflow (should be rare)
                log10_L[g] += -std::numeric_limits<double>::max();
            }
        }
    }

    // Find the maximum log-likelihood
    double max_log10_L = *std::max_element(log10_L.begin(), log10_L.end());

    // Calculate Phred-scaled likelihoods (PL)
    std::vector<int> pl(n_genotypes);
    for (size_t g = 0; g < n_genotypes; ++g) {
        pl[g] = static_cast<int>(std::round(-10.0 * (log10_L[g] - max_log10_L)));
    }

    return pl;
}


double fisher_exact_test(int n11, int n12, int n21, int n22, TestSide test_side) {
    // Input validation
    if (n11 < 0 || n12 < 0 || n21 < 0 || n22 < 0) {
        throw std::invalid_argument("Contingency table entries must be non-negative");
    }

    double left_pvalue, right_pvalue, twoside_pvalue;
    kt_fisher_exact(n11, n12, n21, n22, 
                    &left_pvalue, 
                    &right_pvalue, 
                    &twoside_pvalue);

    double pvalue = -1.0;
    switch(test_side) {
        case TestSide::LESS:      pvalue = left_pvalue;    break;
        case TestSide::GREATER:   pvalue = right_pvalue;   break;
        case TestSide::TWO_SIDED: pvalue = twoside_pvalue; break;
        default:
            throw std::invalid_argument("Invalid test side specification");
    }

    return pvalue; // -1 is an error
}

double fisher_exact_test(const ContingencyTable& table, TestSide test_side) {
    double left_pvalue, right_pvalue, twoside_pvalue;
    kt_fisher_exact(table.n11, table.n12, table.n21, table.n22,
                    &left_pvalue, &right_pvalue, &twoside_pvalue);
    
    double pvalue = -1.0;
    switch(test_side) {
        case TestSide::LESS:      pvalue = left_pvalue;    break;
        case TestSide::GREATER:   pvalue = right_pvalue;   break;
        case TestSide::TWO_SIDED: pvalue = twoside_pvalue; break;
        default:
            throw std::invalid_argument("Invalid test side specification");
    }

    return pvalue; // -1 is an error
}

double wilcoxon_ranksum_test(const std::vector<double>& sample1, const std::vector<double>& sample2) {

    std::vector<double> combined = sample1;
    combined.insert(combined.end(), sample2.begin(), sample2.end());

    // 排序并分配秩
    std::vector<size_t> rank_idx(combined.size());
    std::iota(rank_idx.begin(), rank_idx.end(), 0.0); // 初始化秩
    std::sort(rank_idx.begin(), rank_idx.end(), [&combined](size_t a, size_t b) { return combined[a] > combined[b]; });

    std::vector<double> rankvalues(combined.size());
    for (size_t i(0); i < rank_idx.size(); ++i) 
        rankvalues[i] = i+1; 

    // 处理重复元素
    double ranksum = 0.0, same_n = 1;
    size_t i;
    for (i = 0; i < rank_idx.size(); ++i) {
        if (i > 0 && combined[rank_idx[i]] != combined[rank_idx[i-1]]) {

            if (same_n > 1) {
                double avg_rank = ranksum / same_n; // 平均秩
                for (size_t j = i - same_n; j < i; ++j) {
                    rankvalues[j] = avg_rank; // 分配平均秩
                }
            }

            // 重置
            same_n  = 1;
            ranksum = 0;
        } else if (i > 0) {
            same_n++;
        }
        ranksum += i+1;
    }

    // 处理最后一组重复
    if (same_n > 1) {
        double avg_rank = ranksum / same_n; // 平均秩
        for (size_t j = i - same_n; j < i; ++j) {
            rankvalues[j] = avg_rank; // 分配平均秩
        }
    }

    size_t n1 = sample1.size(), n2 = sample2.size();

    // 计算样本1的秩和
    double smp1_ranksum = 0.0;
    for (size_t i = 0; i < rank_idx.size(); ++i) {
        if (rank_idx[i] < n1) {
            smp1_ranksum += rankvalues[i];
        }
    }

    double e = (double)(n1 * (n1 + n2 + 1)) / 2.0;
    double z = (smp1_ranksum - e) / std::sqrt(double(n1*n2*(n1+n2+1))/12.0);
    double p = 2 * norm_dist(std::abs(z));
    
    // 返回秩和检验 pvalue
    return p;
}

void e_step(const std::vector<double> &obs_allele_freq,
            const std::vector<std::vector<double>> &ind_allele_likelihood, 
            std::vector<std::vector<double>> &ind_allele_post_prob,  // return value, update value inplace 
            std::vector<double> &marginal_likelihood)                // return value, calculate value inplace
{
    size_t n_sample = ind_allele_likelihood.size();
    size_t n_allele = obs_allele_freq.size();

    // 'likelihood' is as the same shape as `ind_allele_likelihood`
    std::vector<std::vector<double>> likelihood(n_sample, std::vector<double>(n_allele, 0));

    // reset the raw value to be 0
    marginal_likelihood = std::vector<double>(n_sample, 0);
    for (size_t i(0); i < n_sample; ++i) {
        for (size_t j = 0; j < n_allele; j++) {  // for _UNIQ_BASES
            likelihood[i][j] = ind_allele_likelihood[i][j] * obs_allele_freq[j];
            marginal_likelihood[i] += likelihood[i][j];
        }

        // Computed the posterior probability of _UNIQ_BASES, change the value inplace.
        for (size_t j(0); j < n_allele; ++j) { 
            // reset the posterior value
            ind_allele_post_prob[i][j] = likelihood[i][j] / marginal_likelihood[i];
        }
    }

    return;
}

void m_step(const std::vector<std::vector<double>> &ind_allele_post_prob, std::vector<double> &obs_allele_freq) {
    size_t n_sample = ind_allele_post_prob.size();  // depth
    size_t n_allele = ind_allele_post_prob[0].size();

    // Reset data
    obs_allele_freq = std::vector<double>(n_allele, 0);
    for (size_t j(0); j < n_allele; ++j) {
        for (size_t i(0); i < n_sample; ++i) {
            obs_allele_freq[j] += ind_allele_post_prob[i][j];
        }
        obs_allele_freq[j] /= (double)(n_sample);  // average.
    }

    return;
}

void EM(const std::vector<std::vector<double>> &ind_allele_likelihood, // n x _UNIQ_BASES.size() matrix, do not change the raw value.
        std::vector<double> &obs_allele_freq,          // retuen value, 1 x _UNIQ_BASES.size(), expect allele frequence, it'll be update inplace here.
        std::vector<double> &log_marginal_likelihood,  // return value
        int iter_num, const float epsilon)
{
    if (iter_num <= 0) {
        throw std::invalid_argument("Iteration number must be positive");
    }
    if (epsilon <= 0) {
        throw std::invalid_argument("Epsilon must be positive");
    }

    size_t n_sample = ind_allele_likelihood.size();
    size_t n_allele = obs_allele_freq.size();

    // n x _UNIQ_BASES.size()matrix, the same shape as 'ind_allele_likelihood'.
    std::vector<std::vector<double>> ind_allele_post_prob = std::vector<std::vector<double>>(
        n_sample, std::vector<double>(n_allele, 0));

    // It's a 1-d array (n x 1) one sample per value, n is sample size
    std::vector<double> marginal_likelihood = std::vector<double>(n_sample, 0);

    // Update the value of 'ind_allele_post_prob' and compute 'marginal_likelihood' in e_step.
    e_step(obs_allele_freq, ind_allele_likelihood, ind_allele_post_prob, marginal_likelihood);

    log_marginal_likelihood = std::vector<double>(n_sample, 0);
    for (size_t i = 0; i < marginal_likelihood.size(); i++) {
        log_marginal_likelihood[i] = log(marginal_likelihood[i]);
    }
    
    // use 'ind_allele_post_prob' to update 'obs_allele_freq'
    m_step(ind_allele_post_prob, obs_allele_freq);
    while (iter_num--) {
        // e step: update ind_allele_post_prob
        e_step(obs_allele_freq, ind_allele_likelihood, ind_allele_post_prob, marginal_likelihood);

        // m step: update the frequence of observed alleles
        m_step(ind_allele_post_prob, obs_allele_freq); 

        double delta = 0, llh;
        for (size_t i = 0; i < marginal_likelihood.size(); i++) {
            llh = log(marginal_likelihood[i]);
            delta += abs(llh - log_marginal_likelihood[i]);
            log_marginal_likelihood[i] = llh;  // update
        }
        
        // Todo: be careful here!!!
        if (delta < epsilon) break;
    }
    // update the lastest expect_allele_freq
    m_step(ind_allele_post_prob, obs_allele_freq);
    
    return;
}
