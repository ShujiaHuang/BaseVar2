/**
 * @file test_bayesian_caller.cpp
 * @brief Unit tests for Bayesian genotype calling functions:
 *        pl_to_likelihoods, hw_genotype_prior, compute_genotype_posterior, compute_dosage_ac
 *
 * Build:
 *   cd tests/io && make test_bayesian_caller
 *
 * Run:
 *   ./test_bayesian_caller
 */

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <cassert>
#include <iomanip>
#include <limits>
#include <stdexcept>

#include "algorithm.h"

// ============================================================================
// Simple test framework
// ============================================================================

static int g_tests_passed = 0;
static int g_tests_failed = 0;

#define CHECK(cond, msg) do { \
    if (!(cond)) { \
        std::cerr << "[FAIL] " << msg << " (line " << __LINE__ << ")" << std::endl; \
        g_tests_failed++; \
    } else { \
        g_tests_passed++; \
    } \
} while(0)

#define CHECK_CLOSE(a, b, tol, msg) CHECK(std::abs((a) - (b)) < (tol), \
    std::string(msg) + " | got=" + std::to_string(a) + " expected=" + std::to_string(b))

// ============================================================================
// 5.1 Tests: pl_to_likelihoods()
// ============================================================================

void test_pl_to_likelihoods() {
    std::cout << "=== test_pl_to_likelihoods ===" << std::endl;

    // Test 1: Standard PL {30, 10, 0}
    // After subtracting min(0): {30, 10, 0}
    // L = {10^-3, 10^-1, 10^0} = {0.001, 0.1, 1.0}
    // Normalized: total = 1.101
    // L_norm = {0.001/1.101, 0.1/1.101, 1.0/1.101}
    {
        std::vector<int> PL = {30, 10, 0};
        auto L = pl_to_likelihoods(PL);
        CHECK(L.size() == 3, "pl_to_likelihoods: size == 3");

        double total = L[0] + L[1] + L[2];
        CHECK_CLOSE(total, 1.0, 1e-10, "pl_to_likelihoods: sum == 1");

        // L[2] should be the largest (PL=0 is best)
        CHECK(L[2] > L[1], "pl_to_likelihoods: L[2] > L[1]");
        CHECK(L[1] > L[0], "pl_to_likelihoods: L[1] > L[0]");

        // Check ratios: L[1]/L[2] = 0.1/1.0 = 0.1 (before normalization, ratio preserved)
        CHECK_CLOSE(L[1] / L[2], 0.1 / 1.0, 1e-6, "pl_to_likelihoods: L[1]/L[2] ratio");
        CHECK_CLOSE(L[0] / L[2], 0.001 / 1.0, 1e-6, "pl_to_likelihoods: L[0]/L[2] ratio");
    }

    // Test 2: All zero PL {0, 0, 0} => uniform
    {
        std::vector<int> PL = {0, 0, 0};
        auto L = pl_to_likelihoods(PL);
        CHECK_CLOSE(L[0], 1.0 / 3.0, 1e-10, "pl_to_likelihoods: uniform L[0]");
        CHECK_CLOSE(L[1], 1.0 / 3.0, 1e-10, "pl_to_likelihoods: uniform L[1]");
        CHECK_CLOSE(L[2], 1.0 / 3.0, 1e-10, "pl_to_likelihoods: uniform L[2]");
    }

    // Test 3: Extreme difference {1000, 0, 500} => numerical stability
    {
        std::vector<int> PL = {1000, 0, 500};
        auto L = pl_to_likelihoods(PL);
        // After subtracting min(0): same
        // L_raw = {10^-100, 1, 10^-50}
        // L[1] should dominate
        CHECK(L[1] > 0.99, "pl_to_likelihoods: extreme L[1] dominates");
        CHECK_CLOSE(L[0] + L[1] + L[2], 1.0, 1e-10, "pl_to_likelihoods: extreme sum == 1");
    }

    // Test 4: Single genotype {0}
    {
        std::vector<int> PL = {0};
        auto L = pl_to_likelihoods(PL);
        CHECK(L.size() == 1, "pl_to_likelihoods: single size");
        CHECK_CLOSE(L[0], 1.0, 1e-10, "pl_to_likelihoods: single == 1.0");
    }

    // Test 5: Empty PL
    {
        std::vector<int> PL = {};
        auto L = pl_to_likelihoods(PL);
        CHECK(L.empty(), "pl_to_likelihoods: empty input => empty output");
    }

    std::cout << "  pl_to_likelihoods tests done." << std::endl;
}

// ============================================================================
// 5.2 Tests: hw_genotype_prior()
// ============================================================================

void test_hw_genotype_prior() {
    std::cout << "=== test_hw_genotype_prior ===" << std::endl;

    // Test 1: Biallelic f=0.5 => {0.25, 0.5, 0.25}
    {
        auto prior = hw_genotype_prior(2, {0.5, 0.5});
        CHECK(prior.size() == 3, "hw_prior: biallelic size == 3");
        CHECK_CLOSE(prior[0], 0.25, 1e-6, "hw_prior: P(0/0)=0.25");
        CHECK_CLOSE(prior[1], 0.50, 1e-6, "hw_prior: P(0/1)=0.50");
        CHECK_CLOSE(prior[2], 0.25, 1e-6, "hw_prior: P(1/1)=0.25");
    }

    // Test 2: Rare variant f=0.01
    {
        auto prior = hw_genotype_prior(2, {0.99, 0.01});
        CHECK_CLOSE(prior[0], 0.9801, 1e-6, "hw_prior: rare P(0/0)=0.9801");
        CHECK_CLOSE(prior[1], 0.0198, 1e-6, "hw_prior: rare P(0/1)=0.0198");
        CHECK_CLOSE(prior[2], 0.0001, 1e-6, "hw_prior: rare P(1/1)=0.0001");
    }

    // Test 3: Triallelic n=3, freqs={0.7, 0.2, 0.1}
    // VCF PL ordering: (0,0), (0,1), (1,1), (0,2), (1,2), (2,2)
    // Priors: 0.49, 0.28, 0.04, 0.14, 0.04, 0.01
    {
        auto prior = hw_genotype_prior(3, {0.7, 0.2, 0.1});
        CHECK(prior.size() == 6, "hw_prior: triallelic size == 6");
        CHECK_CLOSE(prior[0], 0.49, 1e-6, "hw_prior: P(0/0)=0.49");
        CHECK_CLOSE(prior[1], 0.28, 1e-6, "hw_prior: P(0/1)=0.28");
        CHECK_CLOSE(prior[2], 0.04, 1e-6, "hw_prior: P(1/1)=0.04");
        CHECK_CLOSE(prior[3], 0.14, 1e-6, "hw_prior: P(0/2)=0.14");
        CHECK_CLOSE(prior[4], 0.04, 1e-6, "hw_prior: P(1/2)=0.04");
        CHECK_CLOSE(prior[5], 0.01, 1e-6, "hw_prior: P(2/2)=0.01");
    }

    // Test 4: Normalization check for arbitrary input
    {
        auto prior = hw_genotype_prior(2, {0.3, 0.7});
        double total = prior[0] + prior[1] + prior[2];
        CHECK_CLOSE(total, 1.0, 1e-6, "hw_prior: normalization sum == 1");
    }

    // Test 5: AF=0 boundary (should be clamped)
    {
        auto prior = hw_genotype_prior(2, {0.0, 1.0});
        // After clamping: f_ref = 1e-6, f_alt = 1-1e-6
        // P(0/0) = (1e-6)^2 ≈ 0
        // P(0/1) = 2 * 1e-6 * (1-1e-6) ≈ 2e-6
        // P(1/1) = (1-1e-6)^2 ≈ 1
        CHECK(prior[2] > 0.99, "hw_prior: AF=1 clamp => P(1/1) ≈ 1");
        double total = prior[0] + prior[1] + prior[2];
        CHECK_CLOSE(total, 1.0, 1e-6, "hw_prior: AF=1 clamp normalization");
    }

    std::cout << "  hw_genotype_prior tests done." << std::endl;
}

// ============================================================================
// 5.3 Tests: compute_genotype_posterior()
// ============================================================================

void test_compute_genotype_posterior() {
    std::cout << "=== test_compute_genotype_posterior ===" << std::endl;

    // Test 1: High-depth heterozygous PL={100, 0, 100}, AF=0.5
    // Likelihood strongly favors het, prior is uninformative (0.25, 0.5, 0.25)
    // Posterior should strongly favor het
    {
        auto gp = compute_genotype_posterior({100, 0, 100}, 0.5);
        CHECK(gp.best_gt_idx == 1, "posterior: high-depth het => GT=0/1");
        CHECK(gp.gq > 20.0, "posterior: high-depth het => high GQ");
        CHECK_CLOSE(gp.dosage, 1.0, 0.01, "posterior: high-depth het => dosage≈1");
        CHECK_CLOSE(gp.posteriors[0] + gp.posteriors[1] + gp.posteriors[2], 1.0, 1e-10,
                    "posterior: sum == 1");
    }

    // Test 2: Low-depth noise + rare variant: PL={3, 0, 5}, AF=0.001
    // Likelihood slightly favors het (PL=0), but rare AF prior strongly favors hom_ref
    // Prior: P(0/0)≈0.998, P(0/1)≈0.002, P(1/1)≈0.000001
    // Likelihood: L(0/0)=10^-0.3≈0.501, L(0/1)=1, L(1/1)=10^-0.5≈0.316
    // Posterior: P(0/0)∝0.501*0.998≈0.500, P(0/1)∝1*0.002=0.002, P(1/1)∝0.316*1e-6≈3e-7
    // => GT should be 0/0
    {
        auto gp = compute_genotype_posterior({3, 0, 5}, 0.001);
        CHECK(gp.best_gt_idx == 0, "posterior: low-depth rare => GT=0/0");
        CHECK(gp.dosage < 0.05, "posterior: low-depth rare => dosage≈0");
    }

    // Test 3: Low-depth noise + common variant: PL={3, 0, 5}, AF=0.5
    // Prior: P(0/0)=0.25, P(0/1)=0.5, P(1/1)=0.25
    // Likelihood: L(0/0)≈0.501, L(0/1)=1, L(1/1)≈0.316
    // Posterior: P(0/0)∝0.125, P(0/1)∝0.5, P(1/1)∝0.079
    // => GT should still be 0/1, but GQ lower than Test 1
    {
        auto gp = compute_genotype_posterior({3, 0, 5}, 0.5);
        CHECK(gp.best_gt_idx == 1, "posterior: low-depth common => GT=0/1");
        CHECK(gp.gq < 10.0, "posterior: low-depth common => lower GQ");
        CHECK(gp.dosage > 0.5 && gp.dosage < 1.5, "posterior: low-depth common => dosage near 1");
    }

    // Test 4: AF=0 boundary (clamped)
    {
        auto gp = compute_genotype_posterior({10, 0, 5}, 0.0);
        // AF clamped to 1e-6, prior strongly favors 0/0
        // But PL strongly favors 0/1 (PL=0 vs PL=10)
        // Likelihood ratio: L(0/1)/L(0/0) = 10^(10/10) = 10
        // Prior ratio: P(0/1)/P(0/0) ≈ 2e-6
        // Posterior: 0/0 wins
        CHECK(gp.best_gt_idx == 0, "posterior: AF=0 clamp => GT=0/0");
        CHECK_CLOSE(gp.posteriors[0] + gp.posteriors[1] + gp.posteriors[2], 1.0, 1e-10,
                    "posterior: AF=0 sum == 1");
    }

    // Test 5: No information PL={0,0,0} + AF=0.3 => posterior = prior
    {
        auto gp = compute_genotype_posterior({0, 0, 0}, 0.3);
        // With uniform likelihood, posterior = prior = HW(0.7, 0.3)
        // P(0/0) = 0.49, P(0/1) = 0.42, P(1/1) = 0.09
        CHECK_CLOSE(gp.posteriors[0], 0.49, 1e-4, "posterior: no-info => P(0/0)=prior");
        CHECK_CLOSE(gp.posteriors[1], 0.42, 1e-4, "posterior: no-info => P(0/1)=prior");
        CHECK_CLOSE(gp.posteriors[2], 0.09, 1e-4, "posterior: no-info => P(1/1)=prior");
        CHECK(gp.best_gt_idx == 0, "posterior: no-info => GT=0/0 (prior max)");
        // dosage = 0*0.49 + 1*0.42 + 2*0.09 = 0.60
        CHECK_CLOSE(gp.dosage, 0.60, 1e-4, "posterior: no-info => dosage=2*AF=0.6");
    }

    // Test 6: GQ formula verification
    // PL={100, 0, 100}, AF=0.5 => very confident het
    {
        auto gp = compute_genotype_posterior({100, 0, 100}, 0.5);
        CHECK(gp.best_gt_idx == 1, "posterior: confident het => GT=0/1");
        // P(best) should be very close to 1
        double p_best = gp.posteriors[gp.best_gt_idx];
        double expected_gq = -10.0 * std::log10(1.0 - p_best);
        CHECK_CLOSE(gp.gq, expected_gq, 1e-6, "posterior: GQ formula exact match");
    }

    // Test 7: Dosage formula verification
    // PL={0, 50, 100}, AF=0.5
    {
        auto gp = compute_genotype_posterior({0, 50, 100}, 0.5);
        // Manually compute expected dosage
        double d = 0.0;
        for (size_t g = 0; g < 3; ++g) {
            int ac = 0;
            if (g == 1) ac = 1;  // (0,1)
            if (g == 2) ac = 2;  // (1,1)
            d += gp.posteriors[g] * ac;
        }
        CHECK_CLOSE(gp.dosage, d, 1e-10, "posterior: dosage formula exact");
    }

    // Test 8: Empty PL
    {
        auto gp = compute_genotype_posterior({}, 0.5);
        CHECK(gp.posteriors.empty(), "posterior: empty PL => empty posteriors");
        CHECK_CLOSE(gp.dosage, 0.0, 1e-10, "posterior: empty PL => dosage=0");
    }

    std::cout << "  compute_genotype_posterior tests done." << std::endl;
}

// ============================================================================
// 5.4 Tests: compute_dosage_ac()
// ============================================================================

void test_compute_dosage_ac() {
    std::cout << "=== test_compute_dosage_ac ===" << std::endl;

    // Test 1: 3 samples with known dosages
    {
        std::vector<GenotypePosterior> posts(3);
        posts[0].dosage = 0.01;
        posts[1].dosage = 0.95;
        posts[2].dosage = 1.80;

        auto [ac, an] = compute_dosage_ac(posts);
        CHECK_CLOSE(ac, 2.76, 1e-6, "dosage_ac: 3 samples AC=2.76");
        CHECK(an == 6, "dosage_ac: 3 samples AN=6");
    }

    // Test 2: All hom_ref (dosage ≈ 0)
    {
        std::vector<GenotypePosterior> posts(5);
        for (auto& gp : posts) gp.dosage = 0.001;

        auto [ac, an] = compute_dosage_ac(posts);
        CHECK_CLOSE(ac, 0.005, 1e-6, "dosage_ac: all hom_ref AC≈0");
        CHECK(an == 10, "dosage_ac: 5 samples AN=10");
    }

    // Test 3: All het (dosage ≈ 1)
    {
        std::vector<GenotypePosterior> posts(4);
        for (auto& gp : posts) gp.dosage = 0.99;

        auto [ac, an] = compute_dosage_ac(posts);
        CHECK_CLOSE(ac, 3.96, 1e-6, "dosage_ac: all het AC≈4");
        CHECK(an == 8, "dosage_ac: 4 samples AN=8");
    }

    // Test 4: Empty
    {
        std::vector<GenotypePosterior> posts;
        auto [ac, an] = compute_dosage_ac(posts);
        CHECK_CLOSE(ac, 0.0, 1e-10, "dosage_ac: empty AC=0");
        CHECK(an == 0, "dosage_ac: empty AN=0");
    }

    std::cout << "  compute_dosage_ac tests done." << std::endl;
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char *argv[]) {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "========================================" << std::endl;
    std::cout << "Bayesian Genotype Caller Unit Tests" << std::endl;
    std::cout << "========================================" << std::endl;

    test_pl_to_likelihoods();
    test_hw_genotype_prior();
    test_compute_genotype_posterior();
    test_compute_dosage_ac();

    std::cout << "\n========================================" << std::endl;
    std::cout << "Results: " << g_tests_passed << " passed, "
              << g_tests_failed << " failed." << std::endl;
    std::cout << "========================================" << std::endl;

    return g_tests_failed > 0 ? 1 : 0;
}
