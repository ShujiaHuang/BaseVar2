# BaseVar v2.5.0 Release Notes

## Highlights

**Bayesian genotype calling with two-pass architecture** — `basevar caller` now uses a Hardy-Weinberg prior (derived from LRT-estimated population AF) to compute per-sample genotype posteriors. This replaces the previous pure-likelihood `argmin(PL)` approach with a statistically principled Bayesian framework, producing better-calibrated GQ scores and more accurate genotype calls at ultra-low depth.

## New features

### Bayesian posterior genotype calling

- **Two-pass architecture**: First pass collects PL and hard-count AC/AN (unchanged, backward compatible). Second pass uses LRT AF as HW prior to compute per-sample genotype posteriors.
- **Posterior-based GT**: `GT = argmax P(g | D, f̂)` replaces `GT = argmin PL`, reducing noise-driven miscalls at ultra-low depth.
- **Posterior-based GQ**: `GQ = −10 log₁₀(1 − P(best GT))` replaces PL-gap GQ, providing a true confidence measure that reflects data uncertainty.
- **`--gt-mode` parameter**: Control genotype calling mode:
  - `posterior` (default): Bayesian posterior-based GT/GQ + dosage INFO
  - `legacy`: Revert to pre-v2.5.0 pure-likelihood behavior

### Dosage-based INFO fields

- **`INFO/AC_dosage`**: Expected ALT allele count from posterior dosages (continuous, more stable than hard counts at <1× depth)
- **`INFO/AN_dosage`**: Total allele number for dosage calculation (= 2 × N_samples)
- **`INFO/CAF_dosage`**: Dosage-based carrier frequency (= AC_dosage / AN_dosage)

These fields are emitted for bi-allelic sites when `--gt-mode posterior` (default).

## Backward compatibility

All existing VCF fields retain their original values:

| Field | Status |
|-------|--------|
| `INFO/AF`, `INFO/AC`, `INFO/CAF` | Unchanged |
| `FORMAT/PL`, `FORMAT/AD`, `FORMAT/DP` | Unchanged |
| `FORMAT/GT` | Changed (posterior-based); revert with `--gt-mode legacy` |
| `FORMAT/GQ` | Changed (posterior Phred); revert with `--gt-mode legacy` |
| `INFO/AC_dosage`, `INFO/AN_dosage`, `INFO/CAF_dosage` | **New** |

Use `--gt-mode legacy` to produce output identical to pre-v2.5.0 releases.

## Implementation details

### New modules / functions

- `algorithm.cpp`: `pl_to_likelihoods()`, `hw_genotype_prior()`, `compute_genotype_posterior()`, `compute_dosage_ac()`, `alt_count_at_pl_index()`
- `caller_utils.h`: `GenotypePosterior` struct; `VCFSampleAnnotation` gains `posterior`, `dosage` fields; `GQ` type changed from `int` to `double`
- `caller_utils.h`: `AlleleInfo` gains `dosage_counts`, `dosage_total_alleles` fields
- `variant_caller.h`: `BaseTypeARGS` gains `posterior_gt` flag

### Code organization

- Bayesian computation is isolated in `algorithm.cpp` (reusable, testable)
- Two-pass logic is in `variant_caller.cpp::_vcfrecord_in_pos()`
- Legacy code path is fully preserved (no deletion)

## Bug fixes

- Fix `alt_count_at_pl_index()` loop ordering to match VCF PL genotype ordering
- Fix `hw_genotype_prior()` to use VCF PL ordering (j outer, k inner loop)
- Add missing `#include <cstdio>` in `caller_utils.cpp`

## Testing

- New test suite: `tests/io/test_bayesian_caller.cpp` covering posterior computation, HW prior, dosage calculation, and edge cases
- All 60+ unit tests pass

## Migration guide

**No action required** for most users — the new Bayesian mode is the default and produces strictly better results at ultra-low depth.

To preserve exact pre-v2.5.0 output:

```bash
basevar caller --gt-mode legacy -f ref.fa -o out.vcf.gz -L bamfile.list
```
