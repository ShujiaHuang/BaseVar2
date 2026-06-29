# BaseVar v2.4.1 — `basevar motif --fprofile`: CNSLS solver replaces NNLS + post-hoc normalization

**Release date:** 2026-06-13

This is a **patch** release. It fixes two principled errors (PE1, PE2) in the
`basevar motif --fprofile` F-profile decomposition pathway introduced in
v2.4.0. No other subcommand, flag, output column, or default behavior has
been changed.

---

## Bug fixes

### PE1: NNLS + post-hoc normalization ≠ simplex-constrained optimum

**Problem.** v2.4.0 used the Lawson–Hanson NNLS solver (unconstrained
non-negative least squares) followed by a posteriori normalization
`x_i ← x_i / Σx_j` to force weights onto the probability simplex.
This two-step procedure does **not** minimize ‖Ax − b‖² on the simplex
{x ≥ 0, Σx = 1} — the normalized NNLS solution is generically
sub-optimal under the sum-to-one constraint.

**Fix.** Replaced NNLS + normalization with a **Constrained Non-Negative
Least Squares (CNSLS)** solver: projected gradient descent directly on
the probability simplex using the Michelot (1986) / Duchi et al. (2008)
sorting-based simplex projection. The constraint Σx = 1 is now enforced
at optimization time, guaranteeing the true constrained optimum.

**Impact.** Weights produced by `--fprofile` now sum to exactly 1.0 by
construction. A new `-- F-profile fit quality --` summary section reports
`Sum(x)` (always 1.000000) and `||Ax − b||²` residual as a goodness-of-fit
indicator. The old `nnls_solve()` function is retained in the source for
reference but is no longer called by Stage 3.

### PE2: F-profile reference matrix / default parameter mismatch

**Problem.** The six F-profile reference vectors (Zhou et al., PNAS 2023)
were derived using the canonical Lo lab cfDNA method:
`--from-reference --reads both --proper-pair --max-insert-size 1000`.
The tool's *defaults* use a different parameter combination
(read-derived bases, R1 only, no proper-pair filter, no insert-size cap),
creating a silent methodological mismatch that could bias decomposition
weights.

**Fix.** A runtime parameter check now detects when the user's invocation
differs from the canonical method and emits a `[WARN]` block on stderr
listing each mismatched parameter alongside the recommended invocation.
The tool does **not** silently change parameters — it alerts the user and
continues, preserving reproducibility.

---

## What's Changed (files)

| File | Change |
| ---- | ------ |
| `src/motif_counter.h` | Added `fprofile_raw_sum` and `fprofile_residual` fields to `SampleResult`. |
| `src/motif_counter.cpp` | Added `project_onto_simplex()` and `cnls_solve()` (CNSLS solver). Stage 3 now calls `cnls_solve()` instead of `nnls_solve()` + normalization. Added PE2 parameter-warning block. Added fit-quality summary section. |
| `handoff/HANDOFF_motif_counter_method.md` | Fully synchronized: all NNLS→CNSLS references, PE1/PE2 marked [已修复], new CNSLS pseudocode in §6.3, simplex projection pseudocode in §6.3.1, fit-quality fields in §3.1 data flow. Added references 6 (Michelot 1986) and 7 (Duchi et al. 2008). |
| `README.md` | Updated F-profile description: NNLS→CNSLS; FinaleToolkit comparison table; CLI help text. |
| `CMakeLists.txt` | Version bumped 2.4.0 → 2.4.1. |

---

## Tests

- **72 / 72 assertions pass** (same as v2.4.0; no regressions).
- Verified CNSLS output: `sum(x) = 1.000000` (constraint satisfied by
  construction), residual reported as fit-quality indicator.
- PE2 `[WARN]` correctly emitted when parameters deviate from canonical
  Lo lab method; suppressed when all canonical parameters are set.

---

## Backward compatibility

- All six subcommands (`caller`, `pipeline`, `concat`, `subsam`, `motif`,
  `fetalfrac`) are **unchanged** in parameters and output format.
- `--fprofile` output TSV format is **identical** to v2.4.0 (same columns,
  same header). The only difference is that weights are now computed via
  CNSLS instead of NNLS + normalization — numerical values may differ
  slightly but are more accurate.
- New stdout sections (`-- F-profile fit quality --`) are additive; no
  existing output is removed or renamed.
- PE2 warning is informational only; it does not alter program behavior
  or exit code.

---

## Downloads

| Platform | Asset | Notes |
| -------- | ----- | ----- |
| Linux (x86_64) | [`basevar-linux-static`](https://github.com/ShujiaHuang/BaseVar2/releases/download/v2.4.1/basevar-linux-static) | Requires **glibc ≥ 2.35** |
| macOS (arm64 / Intel) | [`basevar-macos-static`](https://github.com/ShujiaHuang/BaseVar2/releases/download/v2.4.1/basevar-macos-static) | Requires **macOS 12+** |

For source builds, see the [README](https://github.com/ShujiaHuang/BaseVar2#installation).

---

**Full Changelog:** https://github.com/ShujiaHuang/BaseVar2/compare/v2.4.0...v2.4.1
