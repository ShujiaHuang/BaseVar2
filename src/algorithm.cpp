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
                             const std::vector<char>& base_quals,
                             double ref_bias) 
{
    // This function calculates the Genotype likelihoods (PL) for given bases and quality against the reference base.
    if (align_bases.size() != base_quals.size()) {
        throw std::invalid_argument("Bases and base quality vectors must be of the same length");
    }

    // Validate ref_bias range
    if (ref_bias < 0.0 || ref_bias > 1.0) {
        throw std::invalid_argument("ref_bias must be between 0.0 and 1.0");
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

    // Reference bias correction:
    // For het genotype (REF, ALT_k), use (1-β) for REF allele and β for ALT allele.
    // For hom or non-REF het, use standard 0.5/0.5 (bias only affects REF vs ALT distinction).
    double w_ref  = 1.0 - ref_bias;  // weight for REF allele in het
    double w_alt  = ref_bias;        // weight for ALT allele in het

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

            // P(R|G): apply reference bias correction for REF/ALT heterozygotes
            double P_R_given_G;
            bool is_ref_het = (allele1 == ref_base && allele2 != ref_base);
            if (is_ref_het && ref_bias != 0.5) {
                // Het (REF, ALT): apply reference bias correction
                // P(R|G) = (1-β) * P(R|REF) + β * P(R|ALT)
                P_R_given_G = w_ref * P_R_given_allele1 + w_alt * P_R_given_allele2;
            } else {
                // Hom or non-REF het: standard diploid model (0.5, 0.5)
                P_R_given_G = 0.5 * P_R_given_allele1 + 0.5 * P_R_given_allele2;
            }

            // Add log10(P(R|G)) to the cumulative log-likelihood
            if (P_R_given_G > 0) {
                log10_L[g] += std::log10(P_R_given_G);
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

double wilcoxon_ranksum_zscore(const std::vector<double>& sample1, const std::vector<double>& sample2) {
    // Return Z-score instead of p-value (GATK convention for RankSum annotations)
    if (sample1.empty() || sample2.empty()) return 0.0;
    
    std::vector<double> combined = sample1;
    combined.insert(combined.end(), sample2.begin(), sample2.end());
    
    // Sort and assign ranks
    std::vector<size_t> rank_idx(combined.size());
    std::iota(rank_idx.begin(), rank_idx.end(), 0.0);
    std::sort(rank_idx.begin(), rank_idx.end(), [&combined](size_t a, size_t b) { return combined[a] > combined[b]; });
    
    std::vector<double> rankvalues(combined.size());
    for (size_t i = 0; i < rank_idx.size(); ++i) {
        rankvalues[i] = i + 1;
    }
    
    // Handle ties (average ranks)
    double ranksum = 0.0, same_n = 1;
    size_t i;
    for (i = 0; i < rank_idx.size(); ++i) {
        if (i > 0 && combined[rank_idx[i]] != combined[rank_idx[i-1]]) {
            if (same_n > 1) {
                double avg_rank = ranksum / same_n;
                for (size_t j = i - same_n; j < i; ++j) {
                    rankvalues[j] = avg_rank;
                }
            }
            same_n = 1;
            ranksum = 0;
        } else if (i > 0) {
            same_n++;
        }
        ranksum += i + 1;
    }
    // Handle last group of ties
    if (same_n > 1) {
        double avg_rank = ranksum / same_n;
        for (size_t j = i - same_n; j < i; ++j) {
            rankvalues[j] = avg_rank;
        }
    }
    
    // Compute rank sum for sample1
    size_t n1 = sample1.size(), n2 = sample2.size();
    double smp1_ranksum = 0.0;
    for (size_t i = 0; i < rank_idx.size(); ++i) {
        if (rank_idx[i] < n1) {
            smp1_ranksum += rankvalues[i];
        }
    }
    
    // Compute Z-score
    double e = (double)(n1 * (n1 + n2 + 1)) / 2.0;
    double var = double(n1 * n2 * (n1 + n2 + 1)) / 12.0;
    if (var < 1e-10) return 0.0;
    double z = (smp1_ranksum - e) / std::sqrt(var);
    
    return z;
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
        log_marginal_likelihood[i] = std::log(marginal_likelihood[i]);
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
            llh = std::log(marginal_likelihood[i]);
            delta += std::abs(llh - log_marginal_likelihood[i]);
            log_marginal_likelihood[i] = llh;  // update
        }
        
        // Todo: be careful here!!!
        if (delta < epsilon) break;
    }
    
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

double hwe_dosage_test(const std::vector<std::vector<double>>& genotype_probs, size_t n_alleles) {
    // Dosage-based HWE chi-square test for low-coverage data.
    // Uses "soft" genotype counts (posterior probability sums) instead of
    // hard genotype calls, following Shriner (2011) Genet Epidemiol 35:632.
    //
    // Key property: Var(O_g^soft) <= Var(O_g^hard) = N*pi_g*(1-pi_g),
    // so the standard chi-square denominator E_g is already conservative
    // (overestimates the true variance), making this test valid but
    // conservative for soft counts.
    
    if (n_alleles < 2 || genotype_probs.empty()) return 1.0;
    
    size_t n_samples = genotype_probs.size();
    size_t n_genotypes = (n_alleles * (n_alleles + 1)) / 2;
    
    // Validate input
    for (const auto& gp : genotype_probs) {
        if (gp.size() != n_genotypes) return 1.0;
    }
    
    // For bi-allelic sites (most common case)
    if (n_alleles == 2) {
        // Compute "soft" genotype counts
        // PL ordering: (0,0)=0, (0,1)=1, (1,1)=2
        double obs_hom_ref = 0.0, obs_het = 0.0, obs_hom_alt = 0.0;
        for (const auto& gp : genotype_probs) {
            obs_hom_ref += gp[0];
            obs_het     += gp[1];
            obs_hom_alt += gp[2];
        }
        
        // Estimate allele frequency from soft counts
        double total = obs_hom_ref + obs_het + obs_hom_alt;
        if (total < 1e-10) return 1.0; // No informative samples
        
        double alt_count = obs_het + 2.0 * obs_hom_alt;
        double p = alt_count / (2.0 * total); // ALT frequency
        double q = 1.0 - p;                    // REF frequency
        
        // Expected genotype counts under HWE
        double exp_hom_ref = q * q * total;
        double exp_het     = 2.0 * p * q * total;
        double exp_hom_alt = p * p * total;
        
        // Standard chi-square statistic (df=1)
        // Conservative for soft counts since Var(O_g^soft) <= E_g
        double chi2 = 0.0;
        if (exp_hom_ref > 1e-10) chi2 += (obs_hom_ref - exp_hom_ref) * (obs_hom_ref - exp_hom_ref) / exp_hom_ref;
        if (exp_het > 1e-10)       chi2 += (obs_het - exp_het) * (obs_het - exp_het) / exp_het;
        if (exp_hom_alt > 1e-10)   chi2 += (obs_hom_alt - exp_hom_alt) * (obs_hom_alt - exp_hom_alt) / exp_hom_alt;
        
        return chi2_test(chi2, 1);
    }
    
    // Multi-allelic: collapse all ALTs into bi-allelic test
    // (simplified approach; full multi-allelic HWE is more complex)
    double obs_hom_ref = 0.0, obs_het = 0.0, obs_hom_alt = 0.0;
    
    for (size_t i = 0; i < n_samples; ++i) {
        const auto& gp = genotype_probs[i];
        obs_hom_ref += gp[0]; // (0,0)
        
        // Sum all het genotypes involving REF: (0,1), (0,2), ...
        for (size_t j = 1; j < n_alleles; ++j) {
            // PL index for genotype (0,j): j*(j+1)/2
            size_t idx = j * (j + 1) / 2;
            if (idx < gp.size()) {
                obs_het += gp[idx];
            }
        }
        
        // Sum all non-REF homozygous: (1,1), (2,2), ...
        for (size_t j = 1; j < n_alleles; ++j) {
            size_t idx = j * (j + 1) / 2 + j; // (j,j)
            if (idx < gp.size()) {
                obs_hom_alt += gp[idx];
            }
        }
    }
    
    double total = obs_hom_ref + obs_het + obs_hom_alt;
    if (total < 1e-10) return 1.0;
    
    double alt_count = obs_het + 2.0 * obs_hom_alt;
    double p = alt_count / (2.0 * total);
    double q = 1.0 - p;
    
    double exp_hom_ref = q * q * total;
    double exp_het     = 2.0 * p * q * total;
    double exp_hom_alt = p * p * total;
    
    // Standard chi-square (conservative for soft counts)
    double chi2 = 0.0;
    if (exp_hom_ref > 1e-10) chi2 += (obs_hom_ref - exp_hom_ref) * (obs_hom_ref - exp_hom_ref) / exp_hom_ref;
    if (exp_het > 1e-10)       chi2 += (obs_het - exp_het) * (obs_het - exp_het) / exp_het;
    if (exp_hom_alt > 1e-10)   chi2 += (obs_hom_alt - exp_hom_alt) * (obs_hom_alt - exp_hom_alt) / exp_hom_alt;
    
    return chi2_test(chi2, 1);
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
