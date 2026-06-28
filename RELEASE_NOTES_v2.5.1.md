# BaseVar v2.5.1 Release Notes

## Overview

BaseVar v2.5.1 is a bug-fix release that addresses all remaining known defects (D1–D9) in the `basevar caller` variant calling module, along with several additional correctness improvements discovered during thorough code audit.

## Defect Fixes in `basevar caller`

### QUAL Computation

- **D1 — Global LRT QUAL**: QUAL is now computed from the overall LRT statistic (final model vs REF-only) rather than the last stepwise LRT lambda. The maximum of global LRT QUAL and original stepwise QUAL is taken to ensure no evidence is lost in ultra-low-coverage scenarios. Degrees of freedom corrected for multi-allelic sites (`df = max(1, n_alleles - 1)`). Log-scale mismatch fixed: both REF-only and EM log-likelihoods now consistently use natural log (ln).
- **D2 — Mono-allelic QUAL**: QUAL for mono-allelic sites is now computed from the actual LRT statistic (`2 × (ln_L_alt − ln_L_ref)`) instead of being hardcoded to 10000. QUAL cap unified to 100000 across all code paths.

### Indel & Annotation

- **D4 — Indel anchor base quality**: INS and DEL quality values now use the upstream anchored base quality (`aligned_pairs[i-1].read_qual[0]`) instead of the whole-read average quality, improving EM accuracy at indel sites.
- **D7 — HWE test**: Added dosage-based Hardy-Weinberg equilibrium test (`hwe_dosage_test()`), using PL→likelihoods as soft genotype counts. Output in VCF INFO field.
- **D8 — RankSum annotations**: Added `BaseQRankSum`, `MQRankSum`, and `ReadPosRankSum` using `wilcoxon_ranksum_zscore()` (directional Z-score, GATK-compatible).

### Other Fixes

- **D5 — Configurable max alleles**: The hardcoded `BIG_N=6` limit is now a configurable `--max-alleles=INT` option (range 1–20, default 6). A warning is emitted when sites are skipped due to exceeding the threshold.
- **D6 — EM convergence**: Removed redundant `m_step()` call after EM convergence loop. The frequency estimate at convergence is now correct.
- **D9 — Group AF order**: Group allele frequency output (`AF_<group>`) now iterates over VCF ALT order instead of internal `active_bases` order, ensuring consistency with the ALT field.

## Additional Improvements

- **Copy constructor fix**: `_final_model_logL` is now properly copied in `BaseType`'s copy constructor.
- **Default constructor initialization**: All numeric members are now initialized in the default constructor.
- **Namespace compliance**: All `log10()` calls replaced with `std::log10()` for strict C++ compliance.

## New Command-Line Option

| Option | Description | Default |
|--------|-------------|---------|
| `--max-alleles=INT` | Maximum number of active alleles at a site. Sites exceeding this are skipped. | 6 |

## Files Changed

| File | Changes |
|------|---------|
| `src/variant_caller.cpp` | D1 global LRT QUAL, D4 indel quality, D7 HWE, D8 RankSum, D9 group AF order, `--max-alleles` CLI |
| `src/variant_caller.h` | `max_alleles` in `BaseTypeARGS` |
| `src/basetype.cpp` | D2 mono-allelic QUAL, D5 `_max_alleles` member, copy constructor fix |
| `src/basetype.h` | `_max_alleles` member, `set_max_alleles()` setter, default constructor init |
| `src/algorithm.cpp` | D6 EM fix, `std::log10` compliance |
| `src/caller_utils.cpp` | VCF header for HWE/RankSum, `std::log10` compliance |
| `src/caller_utils.h` | VCF output logic |
| `CMakeLists.txt` | Version bump to 2.5.1 |
| `src/version.h` | Version bump to 2.5.1 |
| `README.md` | Download links updated, `--ref-bias` and `--max-alleles` in parameter reference |
| `handoff/HANDOFF_caller_module_method.md` | D1–D9 fix status, file index, detailed descriptions |

## Defect Status Summary

All 9 defects (D1–D9) in the caller module are now fixed. The only remaining open item is **PE3** (LRT greedy algorithm does not guarantee global optimum), which is a design limitation requiring algorithmic restructuring.

| Category | Total | Fixed | Remaining |
|----------|-------|-------|-----------|
| Defects (D1–D9) | 9 | 9 | 0 |
| Principle Errors (PE1–PE4) | 4 | 2 (PE1, PE2) | 1 (PE3) + 1 mitigated (PE4) |
