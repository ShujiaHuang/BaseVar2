# BaseVar v2.5.2 Release Notes

## Overview

v2.5.2 is a patch release that fixes a degrees-of-freedom (df) mismatch in the D1 global QUAL recomputation and improves code correctness.

## Bug Fixes

### 1. D1 Global QUAL Degrees-of-Freedom Mismatch (Statistical Bug)

- **Problem**: The D1 global QUAL computation in `_get_pos_variant_info` used `_final_model_logL` (the logL of the step-down reduced model) with df = `ale_bases.size() - 1`. When the step-down procedure reduced alleles from N to M (N > M > 1), the chi-squared statistic represented the N→M improvement, but df=M-1 did not match, resulting in conservative p-values.
- **Fix**: Added `_full_model_logL` member to `BaseType`, storing the full model logL **before** step-down. The D1 global QUAL now uses the full model logL with df=N-1, ensuring the chi-squared statistic and degrees of freedom are consistent.
- **Impact**: QUAL values at multi-allelic sites (>2 alleles) may change slightly. The statistical test is now mathematically correct.

### 2. `abs()` → `std::abs()` (Code Correctness)

- Changed `abs()` to `std::abs()` in the EM convergence check (`algorithm.cpp`) for proper C++ namespace compliance.

### 3. Redundant Comment Cleanup

- Removed outdated comment in `variant_caller.cpp`.

## Changed Files

| File | Change |
|------|--------|
| `src/basetype.h` | Added `_full_model_logL` member, getter, and default initialization |
| `src/basetype.cpp` | Store `_full_model_logL` before step-down loop; sync in copy constructor |
| `src/variant_caller.cpp` | D1 global QUAL now uses `bt.get_full_model_logL()` |
| `src/algorithm.cpp` | `abs()` → `std::abs()` |
| `src/version.h` | Version 2.5.1 → 2.5.2 |
| `CMakeLists.txt` | Version 2.5.1 → 2.5.2 |

## Upgrade Recommendation

All v2.5.1 users are recommended to upgrade to v2.5.2 for more accurate QUAL computation at multi-allelic sites.

## Downloads

- **Linux static binary**: `basevar-v2.5.2-linux-amd64`
- **macOS static binary**: `basevar-v2.5.2-macos-amd64`
- **Source code**: Source code (tar.gz / zip)
