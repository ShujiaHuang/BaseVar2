// Test for dosage-based HWE chi-square test (Shriner 2011)
// Verifies that soft genotype counts produce correct p-values
// under various scenarios using the standard chi-square formula.
//
// Compile:
//   cd build && make -j4 && cd ..
//   g++ -std=c++17 -O2 -I src -I htslib \
//       tests/io/test_hwe_correction.cpp \
//       build/CMakeFiles/basevar.dir/src/algorithm.cpp.o \
//       build/CMakeFiles/basevar.dir/src/basetype.cpp.o \
//       build/CMakeFiles/basevar.dir/src/io/*.o \
//       htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl \
//       -o build/test_hwe_correction
//   ./build/test_hwe_correction

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <cassert>

#include "algorithm.h"

static int n_pass = 0;
static int n_fail = 0;

void check(const std::string& name, bool cond) {
    if (cond) {
        std::cout << "  [PASS] " << name << "\n";
        ++n_pass;
    } else {
        std::cout << "  [FAIL] " << name << "\n";
        ++n_fail;
    }
}

// Test 1: High-coverage hard calls (P in {0,1}) should behave like standard chi-square
void test_high_coverage_hwe() {
    std::cout << "\n=== Test 1: High-coverage HWE equilibrium ===\n";
    // 100 samples: 25 hom_ref, 50 het, 25 hom_alt → p=0.5, perfect HWE
    std::vector<std::vector<double>> gp(100, std::vector<double>(3, 0.0));
    for (size_t i = 0; i < 25; ++i) gp[i] = {1.0, 0.0, 0.0};     // hom_ref
    for (size_t i = 25; i < 75; ++i) gp[i] = {0.0, 1.0, 0.0};    // het
    for (size_t i = 75; i < 100; ++i) gp[i] = {0.0, 0.0, 1.0};   // hom_alt

    double pval = hwe_dosage_test(gp, 2);
    std::cout << "  p-value = " << pval << "\n";
    check("p-value > 0.5 (perfect HWE)", pval > 0.5);
    check("p-value <= 1.0", pval <= 1.0);
}

// Test 2: Low-coverage soft counts should not inflate chi-square
void test_low_coverage_no_inflation() {
    std::cout << "\n=== Test 2: Low-coverage soft counts (no false HWE violation) ===\n";
    // 100 samples with uniform-ish posteriors (simulating ~0.5x coverage)
    // All samples have posteriors close to (0.25, 0.5, 0.25) → AF ≈ 0.5, HWE
    std::vector<std::vector<double>> gp(100, std::vector<double>(3));
    for (size_t i = 0; i < 100; ++i) {
        gp[i] = {0.25, 0.50, 0.25};  // Uniform-like, HWE with p=0.5
    }

    double pval = hwe_dosage_test(gp, 2);
    std::cout << "  p-value = " << pval << "\n";
    // With uniform posteriors, soft counts exactly match HWE expectations
    // chi-square should be 0 → p-value = 1.0
    check("p-value > 0.5 (no false violation)", pval > 0.5);
}

// Test 3: Zero coverage (all prior) → no information, should return high p-value
void test_zero_coverage() {
    std::cout << "\n=== Test 3: Zero coverage (all prior) ===\n";
    // 50 samples, all with uniform priors p=0.3 → posteriors = (0.49, 0.42, 0.09)
    double p = 0.3;
    std::vector<std::vector<double>> gp(50, {
        (1.0-p)*(1.0-p), 2.0*p*(1.0-p), p*p
    });

    double pval = hwe_dosage_test(gp, 2);
    std::cout << "  p-value = " << pval << "\n";
    // When all posteriors equal the prior, obs = expected exactly → chi2 = 0
    check("p-value ≈ 1.0 (no information)", pval > 0.99);
}

// Test 4: Strong HWE violation with high-confidence calls
void test_strong_violation() {
    std::cout << "\n=== Test 4: Strong HWE violation (zero heterozygotes) ===\n";
    // 200 samples: 100 hom_ref + 100 hom_alt, NO heterozygotes → extreme HWE violation
    std::vector<std::vector<double>> gp(200, std::vector<double>(3, 0.0));
    for (size_t i = 0; i < 100; ++i) gp[i] = {1.0, 0.0, 0.0};     // hom_ref
    for (size_t i = 100; i < 200; ++i) gp[i] = {0.0, 0.0, 1.0};   // hom_alt

    double pval = hwe_dosage_test(gp, 2);
    std::cout << "  p-value = " << pval << "\n";
    check("p-value < 1e-10 (strong violation detected)", pval < 1e-10);
    check("p-value >= 0.0", pval >= 0.0);
}

// Test 5: Empty input → return 1.0
void test_empty_input() {
    std::cout << "\n=== Test 5: Empty input ===\n";
    std::vector<std::vector<double>> gp;
    double pval = hwe_dosage_test(gp, 2);
    std::cout << "  p-value = " << pval << "\n";
    check("p-value = 1.0 (empty input)", pval == 1.0);

    // n_alleles < 2
    pval = hwe_dosage_test({{1.0}}, 1);
    std::cout << "  p-value (n_alleles=1) = " << pval << "\n";
    check("p-value = 1.0 (n_alleles < 2)", pval == 1.0);
}

// Test 6: Uncertain genotypes with real deviation → should still detect but with corrected p-value
void test_uncertain_with_deviation() {
    std::cout << "\n=== Test 6: Uncertain genotypes with real HWE deviation ===\n";
    // 100 samples with moderate uncertainty but excess of heterozygotes
    // True AF = 0.5, but observed has 70 het instead of expected 50
    std::vector<std::vector<double>> gp(100, std::vector<double>(3));
    // 15 hom_ref with high confidence
    for (size_t i = 0; i < 15; ++i) gp[i] = {0.95, 0.04, 0.01};
    // 70 het with moderate confidence
    for (size_t i = 15; i < 85; ++i) gp[i] = {0.05, 0.90, 0.05};
    // 15 hom_alt with high confidence
    for (size_t i = 85; i < 100; ++i) gp[i] = {0.01, 0.04, 0.95};

    double pval = hwe_dosage_test(gp, 2);
    std::cout << "  p-value = " << pval << "\n";
    // With excess heterozygotes, should detect HWE deviation
    check("p-value < 0.05 (deviation detected)", pval < 0.05);
    // But corrected p-value should be larger than uncorrected would be
    check("p-value > 0.0 (not zero)", pval > 0.0);
}

// Test 7: Malformed input — gp size doesn't match expected n_genotypes
void test_malformed_input() {
    std::cout << "\n=== Test 7: Malformed input (gp size mismatch) ===\n";
    // n_alleles=2 expects gp.size()=3, but we give gp.size()=2
    std::vector<std::vector<double>> gp = {{0.5, 0.5}, {0.3, 0.7}};
    double pval = hwe_dosage_test(gp, 2);
    std::cout << "  p-value = " << pval << "\n";
    check("p-value = 1.0 (size mismatch)", pval == 1.0);
}

// Test 8: All-zero probabilities → total < 1e-10
void test_all_zero_probs() {
    std::cout << "\n=== Test 8: All-zero probabilities ===\n";
    std::vector<std::vector<double>> gp(10, {0.0, 0.0, 0.0});
    double pval = hwe_dosage_test(gp, 2);
    std::cout << "  p-value = " << pval << "\n";
    check("p-value = 1.0 (total ≈ 0)", pval == 1.0);
}

// Test 9: Multi-allelic site (n_alleles=3, 6 genotypes)
// Note: the multi-allelic branch collapses to 3 categories (hom_ref, REF-het, hom_alt)
// and DROPS non-REF het genotypes like (1,2). This is a known simplification.
void test_multiallelic() {
    std::cout << "\n=== Test 9: Multi-allelic (3 alleles, REF-only genotypes) ===\n";
    // PL ordering for 3 alleles: (0,0)=0, (0,1)=1, (1,1)=2, (0,2)=3, (1,2)=4, (2,2)=5
    // Only use REF-involving genotypes to avoid the (1,2) dropout issue
    // 100 samples: 50 (0,0), 40 (0,1), 5 (1,1), 5 (0,2), 0 (1,2), 0 (2,2)
    // Collapsed: hom_ref=50, het=45, hom_alt=5 → AF≈0.275
    std::vector<std::vector<double>> gp(100, std::vector<double>(6, 0.0));
    size_t idx = 0;
    int counts[] = {50, 40, 5, 5, 0, 0};
    for (int g = 0; g < 6; ++g) {
        for (int c = 0; c < counts[g]; ++c) {
            gp[idx][g] = 1.0;
            ++idx;
        }
    }

    double pval = hwe_dosage_test(gp, 3);
    std::cout << "  p-value = " << pval << "\n";
    check("p-value >= 0.0", pval >= 0.0);
    check("p-value <= 1.0", pval <= 1.0);
}

// Test 10: Multi-allelic with soft counts — all REF-het, no non-REF het
void test_multiallelic_soft() {
    std::cout << "\n=== Test 10: Multi-allelic soft counts (REF-only) ===\n";
    // 50 samples, posteriors concentrated on (0,0) and REF-hets
    // gp = {(0,0), (0,1), (1,1), (0,2), (1,2), (2,2)}
    // Use only REF-involving genotypes to avoid (1,2) dropout
    std::vector<std::vector<double>> gp(50, {0.60, 0.25, 0.05, 0.08, 0.01, 0.01});

    double pval = hwe_dosage_test(gp, 3);
    std::cout << "  p-value = " << pval << "\n";
    check("p-value >= 0.0", pval >= 0.0);
    check("p-value <= 1.0", pval <= 1.0);
}

// Test 11: Rare allele (AF ≈ 0.01) — exercises near-zero expected counts
void test_rare_allele() {
    std::cout << "\n=== Test 11: Rare allele (AF ≈ 0.01) ===\n";
    // 200 samples: 196 hom_ref, 4 het, 0 hom_alt → AF ≈ 0.01
    std::vector<std::vector<double>> gp(200, std::vector<double>(3, 0.0));
    for (size_t i = 0; i < 196; ++i) gp[i] = {1.0, 0.0, 0.0};
    for (size_t i = 196; i < 200; ++i) gp[i] = {0.0, 1.0, 0.0};

    double pval = hwe_dosage_test(gp, 2);
    std::cout << "  p-value = " << pval << "\n";
    // Expected hom_alt = p^2 * N ≈ 0.0001 * 200 ≈ 0.02 → skipped (E < 1e-10? no, > 1e-10)
    // Should still produce valid p-value
    check("p-value > 0.0", pval > 0.0);
    check("p-value <= 1.0", pval <= 1.0);
}

// Test 12: Verify soft count chi-square matches manual computation
void test_chi2_matches_manual() {
    std::cout << "\n=== Test 12: Chi-square matches manual computation ===\n";
    std::vector<std::vector<double>> gp(100, std::vector<double>(3));
    for (size_t i = 0; i < 40; ++i) gp[i] = {0.60, 0.30, 0.10};
    for (size_t i = 40; i < 80; ++i) gp[i] = {0.10, 0.70, 0.20};
    for (size_t i = 80; i < 100; ++i) gp[i] = {0.05, 0.25, 0.70};

    double pval_func = hwe_dosage_test(gp, 2);

    // Manually compute standard chi-square
    double obs_hr=0, obs_het=0, obs_ha=0;
    for (const auto& g : gp) { obs_hr+=g[0]; obs_het+=g[1]; obs_ha+=g[2]; }
    double total = obs_hr + obs_het + obs_ha;
    double af = (obs_het + 2*obs_ha) / (2*total);
    double q = 1.0 - af;
    double exp_hr = q*q*total, exp_het = 2*af*q*total, exp_ha = af*af*total;
    double chi2_manual = 0.0;
    if (exp_hr > 1e-10) chi2_manual += (obs_hr-exp_hr)*(obs_hr-exp_hr)/exp_hr;
    if (exp_het > 1e-10) chi2_manual += (obs_het-exp_het)*(obs_het-exp_het)/exp_het;
    if (exp_ha > 1e-10) chi2_manual += (obs_ha-exp_ha)*(obs_ha-exp_ha)/exp_ha;
    double pval_manual = chi2_test(chi2_manual, 1);

    std::cout << "  func p-value   = " << pval_func << "\n";
    std::cout << "  manual p-value = " << pval_manual << "\n";
    check("func p == manual p (within 1e-12)", std::abs(pval_func - pval_manual) < 1e-12);
}

// Test 13: Single sample edge case
void test_single_sample() {
    std::cout << "\n=== Test 13: Single sample ===\n";
    std::vector<std::vector<double>> gp = {{0.0, 1.0, 0.0}}; // 1 het sample → AF=0.5
    double pval = hwe_dosage_test(gp, 2);
    std::cout << "  p-value = " << pval << "\n";
    // 1 sample, AF=0.5, expected: 0.25 hom_ref, 0.5 het, 0.25 hom_alt
    // Obs: 0, 1, 0 → deviation but N=1, p-value still valid
    check("p-value >= 0.0", pval >= 0.0);
    check("p-value <= 1.0", pval <= 1.0);
}

// Test 14: Monomorphic site — all samples are hom_ref (AF=0)
void test_monomorphic() {
    std::cout << "\n=== Test 14: Monomorphic site (AF=0) ===\n";
    std::vector<std::vector<double>> gp(50, {1.0, 0.0, 0.0});
    double pval = hwe_dosage_test(gp, 2);
    std::cout << "  p-value = " << pval << "\n";
    // AF=0 → expected: N hom_ref, 0 het, 0 hom_alt
    // Obs matches expected → chi2 = 0, but het and hom_alt skipped (E < 1e-10)
    check("p-value ≈ 1.0 (monomorphic)", pval > 0.99);
}

// Test 15: All heterozygotes (extreme excess)
void test_all_heterozygotes() {
    std::cout << "\n=== Test 15: All heterozygotes (extreme excess) ===\n";
    std::vector<std::vector<double>> gp(100, {0.0, 1.0, 0.0});
    double pval = hwe_dosage_test(gp, 2);
    std::cout << "  p-value = " << pval << "\n";
    // AF=0.5, expected: 25 hom_ref, 50 het, 25 hom_alt
    // Obs: 0, 100, 0 → massive excess of het → strong HWE violation
    check("p-value < 1e-10 (extreme het excess)", pval < 1e-10);
}

int main() {
    std::cout << "============================================\n";
    std::cout << " HWE Dosage Test — Standard Chi-Square Tests\n";
    std::cout << "============================================\n";

    test_high_coverage_hwe();
    test_low_coverage_no_inflation();
    test_zero_coverage();
    test_strong_violation();
    test_empty_input();
    test_uncertain_with_deviation();
    test_malformed_input();
    test_all_zero_probs();
    test_multiallelic();
    test_multiallelic_soft();
    test_rare_allele();
    test_chi2_matches_manual();
    test_single_sample();
    test_monomorphic();
    test_all_heterozygotes();

    std::cout << "\n============================================\n";
    std::cout << " Results: " << n_pass << " passed, " << n_fail << " failed\n";
    std::cout << "============================================\n";

    return n_fail > 0 ? 1 : 0;
}
