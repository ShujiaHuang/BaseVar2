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

// ============================================================================
// Genotype posterior computation (Bayesian genotype calling)
// ============================================================================

std::vector<double> pl_to_likelihoods(const std::vector<int>& PL) {
    if (PL.empty()) return {};

    // Subtract min(PL) for numerical stability before conversion
    int min_pl = *std::min_element(PL.begin(), PL.end());

    std::vector<double> likelihoods(PL.size());
    double total = 0.0;
    for (size_t g = 0; g < PL.size(); ++g) {
        // L(g) = 10^(-(PL(g) - min_pl) / 10)
        // After subtracting min, the best genotype has PL=0 => L=1.0
        likelihoods[g] = std::pow(10.0, -(PL[g] - min_pl) / 10.0);
        total += likelihoods[g];
    }

    // Normalize so sum(L) = 1
    if (total > 0) {
        for (auto& l : likelihoods) l /= total;
    } else {
        // Degenerate case: all PLs are identical => uniform
        double uniform = 1.0 / PL.size();
        std::fill(likelihoods.begin(), likelihoods.end(), uniform);
    }

    return likelihoods;
}

std::vector<double> hw_genotype_prior(size_t n_alleles, const std::vector<double>& allele_freqs) {
    if (allele_freqs.size() != n_alleles) {
        throw std::invalid_argument("hw_genotype_prior: allele_freqs size must equal n_alleles");
    }
    if (n_alleles < 2) {
        throw std::invalid_argument("hw_genotype_prior: n_alleles must be >= 2");
    }

    // Clamp allele frequencies to [1e-6, 1-1e-6]
    static const double AF_MIN = 1e-6;
    static const double AF_MAX = 1.0 - 1e-6;
    std::vector<double> freqs(n_alleles);
    for (size_t i = 0; i < n_alleles; ++i) {
        freqs[i] = std::max(AF_MIN, std::min(AF_MAX, allele_freqs[i]));
    }

    // Compute HW prior in VCF PL ordering: (0,0),(0,1),(1,1),(0,2),(1,2),(2,2),...
    // For genotype (k,k): P = f_k^2
    // For genotype (k,j) where k<j: P = 2 * f_k * f_j
    size_t n_genotypes = (n_alleles * (n_alleles + 1)) / 2;
    std::vector<double> prior(n_genotypes);
    size_t idx = 0;
    for (size_t j = 0; j < n_alleles; ++j) {
        for (size_t k = 0; k <= j; ++k) {
            double multiplier = (k == j) ? 1.0 : 2.0;
            prior[idx] = multiplier * freqs[k] * freqs[j];
            ++idx;
        }
    }

    return prior;
}

// Helper: count ALT alleles in genotype at PL index
// VCF PL ordering: (0,0), (0,1), (1,1), (0,2), (1,2), (2,2), ...
// This matches hw_genotype_prior's loop: for j=0..n-1: for k=0..j:
// ALT count for genotype (i,j) = (i>0?1:0) + (j>0?1:0)
static int alt_count_at_pl_index(size_t pl_idx, size_t n_alleles) {
    // Decode PL index to allele pair (k, j) where k <= j
    // VCF ordering: j is the outer loop, k is the inner loop
    size_t idx = 0;
    for (size_t j = 0; j < n_alleles; ++j) {
        for (size_t k = 0; k <= j; ++k) {
            if (idx == pl_idx) {
                return (k > 0 ? 1 : 0) + (j > 0 ? 1 : 0);
            }
            ++idx;
        }
    }
    return 0; // should not reach here
}

GenotypePosterior compute_genotype_posterior(
    const std::vector<int>& PL,
    double af)
{
    GenotypePosterior gp;

    if (PL.empty()) {
        gp.best_gt_idx = 0;
        gp.posteriors = {};
        gp.dosage = 0.0;
        gp.gq = 0.0;
        return gp;
    }

    // For biallelic: n_alleles = 2, allele_freqs = {1-af, af}
    // For multi-allelic: af represents the ALT frequency, simplified to biallelic model
    size_t n_alleles = 2;  // Default: biallelic
    // Infer n_alleles from PL size: n_gt = n*(n+1)/2
    // n_gt=3 => n=2, n_gt=6 => n=3, n_gt=10 => n=4
    {
        size_t n_gt = PL.size();
        size_t n = 2;
        while (n * (n + 1) / 2 < n_gt) ++n;
        if (n * (n + 1) / 2 == n_gt) n_alleles = n;
    }

    // Clamp AF to valid range
    static const double AF_MIN = 1e-6;
    static const double AF_MAX = 1.0 - 1e-6;
    double clamped_af = std::max(AF_MIN, std::min(AF_MAX, af));

    // Step 1: PL -> normalized likelihoods
    std::vector<double> likelihoods = pl_to_likelihoods(PL);

    // Step 2: HW prior
    std::vector<double> prior;
    if (n_alleles == 2) {
        // Fast path for biallelic (most common case)
        prior = hw_genotype_prior(2, {1.0 - clamped_af, clamped_af});
    } else {
        // Multi-allelic: distribute AF equally among ALT alleles
        double alt_freq_each = clamped_af / (n_alleles - 1);
        std::vector<double> freqs = {1.0 - clamped_af};
        for (size_t i = 1; i < n_alleles; ++i) freqs.push_back(alt_freq_each);
        prior = hw_genotype_prior(n_alleles, freqs);
    }

    // Step 3: Posterior = likelihood * prior, normalized
    size_t n_gt = PL.size();
    gp.posteriors.resize(n_gt);
    double total = 0.0;
    for (size_t g = 0; g < n_gt; ++g) {
        gp.posteriors[g] = likelihoods[g] * prior[g];
        total += gp.posteriors[g];
    }
    if (total > 0) {
        for (auto& p : gp.posteriors) p /= total;
    } else {
        // Degenerate: use prior as posterior
        double prior_sum = std::accumulate(prior.begin(), prior.end(), 0.0);
        for (size_t g = 0; g < n_gt; ++g) {
            gp.posteriors[g] = prior[g] / prior_sum;
        }
    }

    // Step 4: Best genotype = argmax posterior
    gp.best_gt_idx = argmax(gp.posteriors.begin(), gp.posteriors.end());

    // Step 5: Dosage = sum_g P(g) * alt_count(g)
    gp.dosage = 0.0;
    for (size_t g = 0; g < n_gt; ++g) {
        gp.dosage += gp.posteriors[g] * alt_count_at_pl_index(g, n_alleles);
    }

    // Step 6: GQ = -10 * log10(1 - P(best_gt))
    double p_best = gp.posteriors[gp.best_gt_idx];
    if (p_best >= 1.0 - 1e-15) {
        gp.gq = 1000.0; // Cap at a large value
    } else if (p_best <= 1e-15) {
        gp.gq = 0.0;
    } else {
        gp.gq = -10.0 * std::log10(1.0 - p_best);
    }

    return gp;
}

std::pair<double, int> compute_dosage_ac(
    const std::vector<GenotypePosterior>& sample_posts)
{
    double expected_ac = 0.0;
    int n_samples = static_cast<int>(sample_posts.size());

    for (const auto& gp : sample_posts) {
        expected_ac += gp.dosage;
    }

    return {expected_ac, 2 * n_samples};
}
