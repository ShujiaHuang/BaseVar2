# BaseVar v2.6.2 Release Notes

## Overview

v2.6.2 is a correctness-focused release that fixes two latent bugs in the variant caller's PL (Phred-scaled genotype likelihood) indexing and AC/AN observation counting, and introduces key improvements to the posterior genotype calling mode. All bugs affect multi-allelic sites; bi-allelic behavior is unchanged.

## Bug Fixes

### Critical: Fix VCF 4.2 PL Ordering in `calculatePL`

**Problem**: `calculatePL()` generated genotype likelihoods in a non-standard order: `i=0..n-1` (outer), `j=i..n-1` (inner). This differs from the VCF 4.2 specification which requires lexicographic ordering: `j=0..n-1` (outer), `k=0..j` (inner).

For bi-allelic sites (n=2), both orderings produce identical results: `(0,0), (0,1), (1,1)`. However, for tri-allelic sites (n=3), indices 2 and 3 were swapped:

| PL Index | Old (non-standard) | VCF 4.2 Standard |
|:---:|:---:|:---:|
| 2 | (0,2) REF/ALT2 | **(1,1) ALT1/ALT1** |
| 3 | (1,1) ALT1/ALT1 | **(0,2) REF/ALT2** |

**Impact**: In posterior mode, PL values (non-standard order) were multiplied with HW prior values (VCF order), causing incorrect genotype posterior probabilities for multi-allelic sites. This could produce wrong GT calls and biased dosage estimates.

**Fix**: Changed `calculatePL()` genotype generation loop to use VCF 4.2 standard ordering. All downstream functions (`hw_genotype_prior`, `decode_pl_index`, `compute_genotype_posterior`) already used VCF ordering, so the fix aligns `calculatePL` with the rest of the codebase.

### Critical: Fix Multi-Allelic GT Decoding in `pl_index_to_genotype`

**Problem**: `pl_index_to_genotype()` used the same non-standard ordering as the old `calculatePL`, causing it to decode multi-allelic PL indices into wrong genotypes in the posterior second pass (where `best_gt_idx` from `compute_genotype_posterior` uses VCF ordering).

**Fix**: Changed `pl_index_to_genotype()` to VCF 4.2 standard ordering, consistent with `decode_pl_index()` and `hw_genotype_prior()`.

### Fix: AC_obs/AN_obs Counting Missing Samples as Called

**Problem**: In posterior mode, zero-coverage samples output as `./.` (missing) still had `gtcode = {0, 0}` (from `argmin` on all-zero PL). The AC_obs/AN_obs counting code treated these as called REF/REF genotypes, adding 2 REF alleles per missing sample to AN_obs.

**Example** (100-sample site, 8 called + 92 missing):
- Before fix: `AC_obs=10, AN_obs=200, AF_obs=0.05`
- After fix: `AC_obs=13, AN_obs=16, AF_obs=0.8125`

**Fix**: Added `!sa.sample_alts.empty()` guard to skip missing/uncalled samples in AC_obs/AN_obs counting.

## What's New

### Dosage-Based AC/AN/AF in Posterior Mode

In posterior mode (`--gt-mode posterior`, the default), INFO fields AC/AN/AF now report **dosage-based** values computed from genotype posterior probabilities:

- `AF` = expected allele frequency from posterior dosage (recommended by BaseVar)
- `AC` = round(expected allele count from dosage)
- `AN` = 2 × N_samples

This provides more accurate allele frequency estimates for ultra-low coverage data compared to the previous reads-based counts. Legacy mode (`--gt-mode legacy`) behavior is unchanged.

### New AC_obs/AN_obs/AF_obs INFO Fields

Three new INFO fields (posterior mode only) provide **GT-based observed allele counts** from hard genotype calls in the first pass:

```
AC_obs=13;AN_obs=16;AF_obs=0.8125
```

These represent the directly observed allele counts from called genotypes, complementing the dosage-based AC/AN/AF. They exclude zero-coverage samples.

### Removed CAF INFO Field

The `CAF` (called allele frequency) field has been removed as it was functionally redundant with `AF_obs` in posterior mode.

## Testing

### 44 New Multi-Allelic Unit Tests

Added comprehensive unit tests covering the PL ordering changes:

- **`test_calculatePL_multiallelic`** (3 test groups): Verifies tri-allelic PL values at VCF-standard positions for single-read and paired-read scenarios
- **`test_pl_index_to_genotype`** (3 test groups): Validates VCF index decoding for bi-allelic, tri-allelic, and tetra-allelic sites — critically verifying indices 2 and 3
- **`test_compute_genotype_posterior_multiallelic`** (4 test groups): Tests vector overload with tri-allelic PL, per-allele dosage computation, HW prior verification, and PL size mismatch handling

**Total: 109 unit tests passed** (65 original + 44 new).

### Real-Data Verification

Validated against 100-sample NIPT data (chr11 + chr17, 93 variant sites):
- AF = AC/AN (dosage-based) consistency: 93/93 sites ✓
- AF_obs = AC_obs/AN_obs consistency: 93/93 sites ✓
- AC_obs vs manual GT counting: 4/4 sites exact match ✓
- Legacy mode output: unchanged ✓

## Changed Files

| File | Change |
|------|--------|
| `src/algorithm.cpp` | `calculatePL`: VCF 4.2 PL ordering; `compute_genotype_posterior`: multi-allelic vector overload + `decode_pl_index` helper |
| `src/algorithm.h` | `GenotypePosterior.per_allele_dosage`; vector overload declaration; `compute_dosage_ac` return type update |
| `src/caller_utils.cpp` | `pl_index_to_genotype`: VCF 4.2 ordering; `process_sample_variant`: store `per_allele_dosage`; VCF header updates |
| `src/caller_utils.h` | `VCFSampleAnnotation.per_allele_dosage`; `AlleleInfo.dosage_counts` type change |
| `src/variant_caller.cpp` | Posterior second pass: multi-allelic support; AC_obs/AN_obs missing-sample fix; INFO field construction |
| `tests/io/Makefile` | Link `caller_utils.o` + `basetype.o` for expanded test coverage |
| `tests/io/test_bayesian_caller.cpp` | 44 new multi-allelic unit tests |
| `CMakeLists.txt` | Version 2.6.1 → 2.6.2 |

## Upgrade Recommendation

All users of posterior mode (`--gt-mode posterior`, the default) are recommended to upgrade. The VCF PL ordering fix ensures correct multi-allelic genotype calling. Legacy mode (`--gt-mode legacy`) users are unaffected.
