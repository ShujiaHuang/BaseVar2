# BaseVar v2.6.3 Release Notes

## New Features

### AC_GT / AN_GT / AF_GT — GT-based allele count INFO fields

Added three new INFO fields that count alleles **directly from the VCF GT column** (posterior genotype calls with population prior):

| Field | Type | Description |
| ----- | ---- | ----------- |
| `AC_GT` | Integer (Number=A) | Per-ALT allele count from the VCF GT column |
| `AN_GT` | Integer (Number=1) | Total alleles in non-missing VCF GT calls |
| `AF_GT` | Float (Number=A) | Allele frequency = AC_GT / AN_GT |

These complement the existing two tiers:

- **Posterior (recommended)**: `AF`, `AC`, `AN` — expected values from posterior probabilities (with population prior, distinct from GT)
- **Observed**: `AC_obs`, `AN_obs`, `AF_obs` — discrete counts from `argmin(PL)` genotypes (without population prior, distinct from GT)

## Improvements

### VCF header description overhaul

All VCF header descriptions have been systematically audited and rewritten for clarity:

- **Removed internal terminology**: Replaced "posterior mode", "legacy mode", "hard genotype", and "Bayesian mode" with user-friendly computation-based descriptions.
- **GT-source clarification**: Added "distinct from GT" to `AC/AN/AF` and `AC_obs/AN_obs/AF_obs` descriptions, making it explicit that these are NOT derived from the VCF GT column.
- **AF recommendation**: The `AF` field description now includes "(recommended by BaseVar)" to guide users toward the preferred allele frequency metric.
- **AN description**: Updated to clarify it equals 2 × number of samples.
- **AF_group description**: Made more generic; removed LRT-specific reference since the field may use different estimation methods.
- **BaseQRankSum**: Fixed inconsistent capitalization ("Alt Vs. Ref" → "Alt vs. Ref").

### VCF 4.2 compliance fixes

- `AD`: `Type=String` → `Type=Integer` (VCF 4.2 spec compliance)
- `DP4`: `Number=A` → `Number=.` (variable-length list, not per-allele)

### Single-region VCF output optimization

When calling variants in a single region, BaseVar now directly renames the temp VCF file to the output instead of performing a read-write merge. This eliminates unnecessary I/O for single-region runs.

### Progress reporting improvements

- Loading progress messages now include timestamps.
- Variant calling progress interval increased from 10,000 to 100,000 positions, with elapsed time and variant count in the output.

### README updates

- Updated the "Bayesian genotype calling" section to reflect the current three-tier INFO field structure (posterior / observed / GT-based).
- Removed references to deprecated field names (`AC_dosage`, `AN_dosage`, `CAF_dosage`).

## Files Changed

- `src/caller_utils.cpp` — VCF header definitions, `merge_file_by_line` usage
- `src/caller_utils.h` — `AlleleInfo` struct: new `gt_counts` and `gt_total_alleles` fields
- `src/variant_caller.cpp` — AC_GT/AN_GT/AF_GT computation, INFO output, single-file rename optimization, progress reporting
- `src/variant_caller.h` — `_variants_discovery` signature update
- `CMakeLists.txt` — Version bump to 2.6.3
- `README.md` — Updated INFO field documentation
