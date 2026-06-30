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
#include "caller_utils.h"

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

    // Test 1: 3 samples with known per-allele dosages (bi-allelic)
    {
        std::vector<GenotypePosterior> posts(3);
        posts[0].per_allele_dosage = {0.01};
        posts[1].per_allele_dosage = {0.95};
        posts[2].per_allele_dosage = {1.80};

        auto [ac, an] = compute_dosage_ac(posts);
        CHECK(ac.size() == 1, "dosage_ac: bi-allelic returns 1 ALT");
        CHECK_CLOSE(ac[0], 2.76, 1e-6, "dosage_ac: 3 samples AC[0]=2.76");
        CHECK(an == 6, "dosage_ac: 3 samples AN=6");
    }

    // Test 2: All hom_ref (dosage ~ 0)
    {
        std::vector<GenotypePosterior> posts(5);
        for (auto& gp : posts) gp.per_allele_dosage = {0.001};

        auto [ac, an] = compute_dosage_ac(posts);
        CHECK_CLOSE(ac[0], 0.005, 1e-6, "dosage_ac: all hom_ref AC~0");
        CHECK(an == 10, "dosage_ac: 5 samples AN=10");
    }

    // Test 3: All het (dosage ~ 1)
    {
        std::vector<GenotypePosterior> posts(4);
        for (auto& gp : posts) gp.per_allele_dosage = {0.99};

        auto [ac, an] = compute_dosage_ac(posts);
        CHECK_CLOSE(ac[0], 3.96, 1e-6, "dosage_ac: all het AC~4");
        CHECK(an == 8, "dosage_ac: 4 samples AN=8");
    }

    // Test 4: Empty
    {
        std::vector<GenotypePosterior> posts;
        auto [ac, an] = compute_dosage_ac(posts);
        CHECK(ac.empty(), "dosage_ac: empty returns empty vector");
        CHECK(an == 0, "dosage_ac: empty AN=0");
    }

    // Test 5: Multi-allelic (2 ALTs)
    {
        std::vector<GenotypePosterior> posts(3);
        posts[0].per_allele_dosage = {0.10, 0.05};  // sample 1: low dosage for both ALTs
        posts[1].per_allele_dosage = {0.95, 0.02};  // sample 2: high dosage for ALT1
        posts[2].per_allele_dosage = {0.03, 1.85};  // sample 3: high dosage for ALT2

        auto [ac, an] = compute_dosage_ac(posts);
        CHECK(ac.size() == 2, "dosage_ac: multi-allelic returns 2 ALTs");
        CHECK_CLOSE(ac[0], 1.08, 1e-6, "dosage_ac: ALT1 AC=1.08");
        CHECK_CLOSE(ac[1], 1.92, 1e-6, "dosage_ac: ALT2 AC=1.92");
        CHECK(an == 6, "dosage_ac: 3 samples AN=6");
    }

    std::cout << "  compute_dosage_ac tests done." << std::endl;
}

// ============================================================================
// 5.5 Tests: calculatePL VCF ordering for multi-allelic (tri-allelic)
// ============================================================================

void test_calculatePL_multiallelic() {
    std::cout << "=== test_calculatePL_multiallelic ===" << std::endl;

    // Tri-allelic: REF=A, ALT1=C, ALT2=G
    // VCF PL ordering: (0,0),(0,1),(1,1),(0,2),(1,2),(2,2)
    //                   idx0  idx1  idx2  idx3  idx4  idx5

    // Test 1: Single read "C" (Q30) — verify PL values at VCF-standard positions
    // Note: ref_bias (default 0.5) makes ref-het P(R|G) = 0.5*P(R|REF) + 0.5*P(R|ALT),
    // while hom and non-ref het use P(R|G) = 0.5*P(R|a1) + 0.5*P(R|a2).
    // For a single C read: ALT1/ALT1 (PL~0) < REF/ALT1 (PL~3) because ref-het is halved.
    {
        std::vector<std::string> reads = {"C"};
        std::vector<char> quals = {static_cast<char>(30 + 33)};

        auto pl = calculatePL("A", {"C", "G"}, reads, quals, 0.5);
        CHECK(pl.size() == 6, "calcPL: triallelic has 6 PL values");

        // PL[0]=REF/REF should be high (read C doesn't match REF A at all)
        CHECK(pl[0] > 20, "calcPL: PL[0]=REF/REF is high for C read");

        // CRITICAL: verify VCF ordering — PL[2]=ALT1/ALT1 should be lower than
        // PL[3]=REF/ALT2 (C read matches ALT1 perfectly, not ALT2)
        CHECK(pl[2] < pl[3], "calcPL: PL[2]=ALT1/ALT1 < PL[3]=REF/ALT2 for C read");

        // PL[1]=REF/ALT1 should be moderate (C matches ALT1 allele)
        CHECK(pl[1] < pl[0], "calcPL: PL[1]=REF/ALT1 < PL[0]=REF/REF for C read");

        // PL[4]=ALT1/ALT2 and PL[5]=ALT2/ALT2 should be high (no G reads)
        CHECK(pl[4] > pl[2], "calcPL: PL[4]=ALT1/ALT2 > PL[2]=ALT1/ALT1 for C read");
        CHECK(pl[5] > pl[2], "calcPL: PL[5]=ALT2/ALT2 > PL[2]=ALT1/ALT1 for C read");
    }

    // Test 2: Single read "G" — verify ALT2-related genotypes have lower PL
    {
        std::vector<std::string> reads = {"G"};
        std::vector<char> quals = {static_cast<char>(30 + 33)};

        auto pl = calculatePL("A", {"C", "G"}, reads, quals, 0.5);

        // PL[3]=REF/ALT2 and PL[5]=ALT2/ALT2 should be lower than ALT1 genotypes
        CHECK(pl[3] < pl[1], "calcPL: PL[3]=REF/ALT2 < PL[1]=REF/ALT1 for G read");
        CHECK(pl[5] < pl[2], "calcPL: PL[5]=ALT2/ALT2 < PL[2]=ALT1/ALT1 for G read");
    }

    // Test 3: Reads "C" + "G" — best GT should be at VCF idx 4 (ALT1/ALT2)
    {
        std::vector<std::string> reads = {"C", "G"};
        std::vector<char> quals = {static_cast<char>(30 + 33), static_cast<char>(30 + 33)};

        auto pl = calculatePL("A", {"C", "G"}, reads, quals, 0.5);
        int min_idx = std::min_element(pl.begin(), pl.end()) - pl.begin();
        CHECK(min_idx == 4, "calcPL: C+G reads => best GT at VCF idx 4 (ALT1/ALT2)");

        // PL[4] should be 0 (best)
        CHECK(pl[4] == 0, "calcPL: PL[4]=ALT1/ALT2 is 0 for C+G reads");

        // Verify REF/ALT genotypes: PL[1] (REF/ALT1) and PL[3] (REF/ALT2) should
        // both be moderate, while PL[0] (REF/REF) should be high
        CHECK(pl[0] > pl[1], "calcPL: REF/REF > REF/ALT1 for C+G reads");
        CHECK(pl[0] > pl[3], "calcPL: REF/REF > REF/ALT2 for C+G reads");
    }

    std::cout << "  calculatePL multi-allelic tests done." << std::endl;
}

// ============================================================================
// 5.6 Tests: pl_index_to_genotype VCF ordering for multi-allelic
// ============================================================================

void test_pl_index_to_genotype() {
    std::cout << "=== test_pl_index_to_genotype ===" << std::endl;

    // Bi-allelic (n=2): should match both old and new ordering
    {
        auto [a0, b0] = pl_index_to_genotype(0, 2);
        CHECK(a0 == 0 && b0 == 0, "pl2gt: bi-allelic idx 0 = (0,0)");
        auto [a1, b1] = pl_index_to_genotype(1, 2);
        CHECK(a1 == 0 && b1 == 1, "pl2gt: bi-allelic idx 1 = (0,1)");
        auto [a2, b2] = pl_index_to_genotype(2, 2);
        CHECK(a2 == 1 && b2 == 1, "pl2gt: bi-allelic idx 2 = (1,1)");
    }

    // Tri-allelic (n=3): VCF standard ordering
    // idx 0=(0,0), idx 1=(0,1), idx 2=(1,1), idx 3=(0,2), idx 4=(1,2), idx 5=(2,2)
    {
        auto [a0, b0] = pl_index_to_genotype(0, 3);
        CHECK(a0 == 0 && b0 == 0, "pl2gt: triallelic idx 0 = (0,0) REF/REF");

        auto [a1, b1] = pl_index_to_genotype(1, 3);
        CHECK(a1 == 0 && b1 == 1, "pl2gt: triallelic idx 1 = (0,1) REF/ALT1");

        // CRITICAL: idx 2 must be (1,1) in VCF ordering, NOT (0,2) as in old ordering
        auto [a2, b2] = pl_index_to_genotype(2, 3);
        CHECK(a2 == 1 && b2 == 1, "pl2gt: triallelic idx 2 = (1,1) ALT1/ALT1 [VCF standard]");

        // CRITICAL: idx 3 must be (0,2) in VCF ordering, NOT (1,1) as in old ordering
        auto [a3, b3] = pl_index_to_genotype(3, 3);
        CHECK(a3 == 0 && b3 == 2, "pl2gt: triallelic idx 3 = (0,2) REF/ALT2 [VCF standard]");

        auto [a4, b4] = pl_index_to_genotype(4, 3);
        CHECK(a4 == 1 && b4 == 2, "pl2gt: triallelic idx 4 = (1,2) ALT1/ALT2");

        auto [a5, b5] = pl_index_to_genotype(5, 3);
        CHECK(a5 == 2 && b5 == 2, "pl2gt: triallelic idx 5 = (2,2) ALT2/ALT2");
    }

    // Tetra-allelic (n=4): verify key indices
    // VCF: (0,0)=0, (0,1)=1, (1,1)=2, (0,2)=3, (1,2)=4, (2,2)=5, (0,3)=6, (1,3)=7, (2,3)=8, (3,3)=9
    {
        auto [a2, b2] = pl_index_to_genotype(2, 4);
        CHECK(a2 == 1 && b2 == 1, "pl2gt: tetraallelic idx 2 = (1,1)");

        auto [a6, b6] = pl_index_to_genotype(6, 4);
        CHECK(a6 == 0 && b6 == 3, "pl2gt: tetraallelic idx 6 = (0,3) REF/ALT3");

        auto [a9, b9] = pl_index_to_genotype(9, 4);
        CHECK(a9 == 3 && b9 == 3, "pl2gt: tetraallelic idx 9 = (3,3) ALT3/ALT3");
    }

    std::cout << "  pl_index_to_genotype tests done." << std::endl;
}

// ============================================================================
// 5.7 Tests: compute_genotype_posterior with vector overload (multi-allelic)
// ============================================================================

void test_compute_genotype_posterior_multiallelic() {
    std::cout << "=== test_compute_genotype_posterior (multi-allelic) ===" << std::endl;

    // Test 1: Tri-allelic posterior with known PL and allele frequencies
    // PL in VCF order: (0,0)=50, (0,1)=0, (1,1)=60, (0,2)=40, (1,2)=30, (2,2)=70
    // Likelihood favors REF/ALT1 het (idx 1, PL=0)
    // Allele freqs: REF=0.6, ALT1=0.25, ALT2=0.15
    {
        std::vector<int> PL = {50, 0, 60, 40, 30, 70};
        std::vector<double> freqs = {0.6, 0.25, 0.15};
        auto gp = compute_genotype_posterior(PL, freqs);

        CHECK(gp.posteriors.size() == 6, "posterior_multi: 6 genotypes");

        // Sum of posteriors should be 1
        double sum = 0;
        for (auto p : gp.posteriors) sum += p;
        CHECK_CLOSE(sum, 1.0, 1e-10, "posterior_multi: sum == 1");

        // Best GT should be idx 1 (REF/ALT1) — strongly favored by likelihood
        CHECK(gp.best_gt_idx == 1, "posterior_multi: best GT = idx 1 (REF/ALT1)");

        // Per-allele dosage: ALT1 should be close to 1, ALT2 close to 0
        CHECK(gp.per_allele_dosage.size() == 2, "posterior_multi: 2 ALT alleles");
        CHECK(gp.per_allele_dosage[0] > 0.5, "posterior_multi: ALT1 dosage > 0.5");
        CHECK(gp.per_allele_dosage[1] < 0.5, "posterior_multi: ALT2 dosage < 0.5");

        // Total dosage should be close to 1 (one ALT allele)
        CHECK_CLOSE(gp.dosage, gp.per_allele_dosage[0] + gp.per_allele_dosage[1], 1e-10,
                    "posterior_multi: total dosage == sum of per-allele dosages");
    }

    // Test 2: Tri-allelic with PL favoring ALT1/ALT2 het (VCF idx 4)
    {
        std::vector<int> PL = {80, 30, 50, 25, 0, 40};
        std::vector<double> freqs = {0.5, 0.3, 0.2};
        auto gp = compute_genotype_posterior(PL, freqs);

        // Best GT should be idx 4 (ALT1/ALT2) — PL=0 there
        CHECK(gp.best_gt_idx == 4, "posterior_multi: ALT1/ALT2 het best GT = idx 4");

        // Both ALT alleles should have dosage > 0.5
        CHECK(gp.per_allele_dosage[0] > 0.3, "posterior_multi: ALT1 dosage > 0.3 for ALT1/ALT2 het");
        CHECK(gp.per_allele_dosage[1] > 0.3, "posterior_multi: ALT2 dosage > 0.3 for ALT1/ALT2 het");
    }

    // Test 3: PL size mismatch (bi-allelic PL=3 but freqs has 3 elements)
    // Should fall back to empty
    {
        std::vector<int> PL = {10, 0, 20};  // bi-allelic: 3 genotypes
        std::vector<double> freqs = {0.5, 0.3, 0.2};  // tri-allelic: expects 6
        auto gp = compute_genotype_posterior(PL, freqs);
        CHECK(gp.posteriors.empty(), "posterior_multi: PL size mismatch => empty");
    }

    // Test 4: Uniform prior (all freqs equal) + no-info PL
    // Posterior should equal prior = HW equilibrium
    {
        std::vector<int> PL = {0, 0, 0, 0, 0, 0};
        std::vector<double> freqs = {1.0/3, 1.0/3, 1.0/3};
        auto gp = compute_genotype_posterior(PL, freqs);

        // HW prior: P(k,j) = mult * f_k * f_j
        // P(0,0) = 1/9, P(0,1) = 2/9, P(1,1) = 1/9, P(0,2) = 2/9, P(1,2) = 2/9, P(2,2) = 1/9
        CHECK_CLOSE(gp.posteriors[0], 1.0/9, 1e-4, "posterior_multi: uniform => P(0,0)=1/9");
        CHECK_CLOSE(gp.posteriors[1], 2.0/9, 1e-4, "posterior_multi: uniform => P(0,1)=2/9");
        CHECK_CLOSE(gp.posteriors[2], 1.0/9, 1e-4, "posterior_multi: uniform => P(1,1)=1/9");
        CHECK_CLOSE(gp.posteriors[3], 2.0/9, 1e-4, "posterior_multi: uniform => P(0,2)=2/9");
        CHECK_CLOSE(gp.posteriors[4], 2.0/9, 1e-4, "posterior_multi: uniform => P(1,2)=2/9");
        CHECK_CLOSE(gp.posteriors[5], 1.0/9, 1e-4, "posterior_multi: uniform => P(2,2)=1/9");

        // Per-allele dosage: each ALT should have 2/3 (= 2 * 1/3)
        CHECK_CLOSE(gp.per_allele_dosage[0], 2.0/3, 1e-4, "posterior_multi: uniform => ALT1 dosage=2/3");
        CHECK_CLOSE(gp.per_allele_dosage[1], 2.0/3, 1e-4, "posterior_multi: uniform => ALT2 dosage=2/3");
        CHECK_CLOSE(gp.dosage, 4.0/3, 1e-4, "posterior_multi: uniform => total dosage=4/3");
    }

    std::cout << "  compute_genotype_posterior multi-allelic tests done." << std::endl;
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
    test_calculatePL_multiallelic();
    test_pl_index_to_genotype();
    test_compute_genotype_posterior_multiallelic();

    std::cout << "\n========================================" << std::endl;
    std::cout << "Results: " << g_tests_passed << " passed, "
              << g_tests_failed << " failed." << std::endl;
    std::cout << "========================================" << std::endl;

    return g_tests_failed > 0 ? 1 : 0;
}
