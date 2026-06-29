# BaseVar v2.4.0 — `basevar motif --fprofile`: F-profile NNLS decomposition for cfDNA tissue deconvolution

**Release date:** 2026-06-06

This is a **minor** release. It adds F-profile decomposition to `basevar motif`,
enabling cfDNA tissue-of-origin deconvolution via Non-Negative Least Squares
(NNLS) against the six "founder" end-motif profiles discovered by Zhou et al.
(PNAS, 2023). No existing flag, output column, or default behavior of any
subcommand has been changed.

---

## Highlights

- **New `--fprofile` flag** for `basevar motif`: opt-in F-profile NNLS
  decomposition that resolves each sample's observed 4-mer end-motif
  frequency vector into a weighted combination of the six Zhou et al. (2023)
  founder profiles.
- **`--fprofile-output FILE`**: writes a long-format TSV (one row per sample
  per F-profile) with normalized weights summing to 1.0, interpretable as
  estimated tissue proportions.
- **Zero new dependencies**: the six 256-element reference vectors are
  embedded at compile time (`src/fprofile_data.h`); the NNLS solver
  (Lawson–Hanson with Gauss–Jordan inner loop) is implemented in-tree with
  no external linear-algebra library.
- **k ≠ 4 guard**: F-profile decomposition is only defined for 4-mers; when
  `--fprofile` is used with any other `-k` value, a `[WARN]` is emitted and
  the flag is silently disabled.
- **72 / 72 unit-test assertions pass** (was 57 in v2.3.1; +15 new).

---

## What's New

### 1. F-profile NNLS decomposition (`--fprofile`)

The six F-profiles ("F-profile I" … "F-profile VI") are the non-negative
matrix factorization (NMF) components learned by Zhou et al. from >14,000
cfDNA samples in the DELFI cohort. Each profile is a probability distribution
over the 256 possible 4-mer end motifs (lexicographic order AAAA…TTTT).

Given an observed sample frequency vector **b** (256 × 1) and the reference
matrix **A** (256 × 6), the NNLS solver finds weights **x** ≥ 0 minimizing
‖**Ax** − **b**‖². The raw weights are then normalized to sum to 1.0.

#### Reference

> Zhou, Z., Ma, R., Huang, X., …, & Sun, K. (2023).
> Cell-free DNA end motifs characterize cancer fragmentation and detect
> early cancer. *Proceedings of the National Academy of Sciences*, 120(17).
> doi:10.1073/pnas.2208138120

### 2. Output formats

**Stdout summary** (appended to the existing `basevar motif` summary):

```
-- F-profile weights --
Sample                  F-profile I   F-profile II  F-profile III F-profile IV  F-profile V   F-profile VI
ERS225193               0.000000      0.000000      0.000000      0.000000      0.000000      1.000000
```

**TSV file** (written when `--fprofile-output FILE` is used):

```
#sample	f_profile	weight
ERS225193	F-profile I	0.000000
ERS225193	F-profile II	0.000000
...
```

### 3. CLI additions

| Flag | Description |
| ---- | ----------- |
| `--fprofile` | Enable F-profile NNLS decomposition (requires k=4). |
| `--fprofile-output FILE` | Write normalized F-profile weights to FILE (TSV). Implies `--fprofile`. |

### 4. Implementation notes

- **`src/fprofile_data.h`** (new): compile-time embedding of the 256 × 6
  reference matrix, 6 profile names, and 256 kmer strings, auto-generated
  from `FinaleToolkit/src/finaletoolkit/frag/data/end_motif_f_profiles.tsv`.
- **NNLS solver** (`nnls_solve()` in `motif_counter.cpp`): Lawson–Hanson
  algorithm with a Gauss–Jordan elimination inner solver for the ≤ 6 × 6
  Gram sub-systems. Converges in at most 6 outer iterations (one per
  profile).
- **No new header or library dependencies**: the solver uses only
  `<vector>`, `<cmath>`, `<numeric>`, and `<algorithm>`.

---

## Tests

- File: `tests/io/test_motif_counter.cpp` — Test 8 (a–d), +15 new assertions.
- **72 / 72 assertions pass**, covering:
  - `--fprofile --fprofile-output` produces 6 weights summing to ~1.0;
    TSV is well-formed with correct header and row count.
  - `--fprofile-output` alone implies `--fprofile` (no explicit flag needed).
  - Without `--fprofile`, `fprofile_weights` remains empty.
  - `--fprofile` with k=5 produces a warning and empty weights (k ≠ 4 guard).

---

## Build & release

- Project version bumped **2.3.1 → 2.4.0** (`CMakeLists.txt`).
- README pre-built binary download URLs updated to `v2.4.0`.
- `src/fprofile_data.h` is a new header included via `#include` (no
  build-system change required).
- The GitHub Actions static-build workflow publishes
  `basevar-linux-static` and `basevar-macos-static` as release assets
  automatically once the `v2.4.0` tag is pushed.

---

## Backward compatibility

- All six existing subcommands (`caller`, `pipeline`, `concat`, `subsam`,
  `motif`, `fetalfrac`) are **unchanged** in semantics, parameters, and
  output format.
- `--fprofile` is an **opt-in** flag; without it, `basevar motif` behaves
  exactly as in v2.3.1.
- No new required runtime dependencies.

---

## Downloads

| Platform | Asset | Notes |
| -------- | ----- | ----- |
| Linux (x86_64) | [`basevar-linux-static`](https://github.com/ShujiaHuang/BaseVar2/releases/download/v2.4.0/basevar-linux-static) | Requires **glibc ≥ 2.35** |
| macOS (arm64 / Intel) | [`basevar-macos-static`](https://github.com/ShujiaHuang/BaseVar2/releases/download/v2.4.0/basevar-macos-static) | Requires **macOS 12+** |

For source builds, see the [README](https://github.com/ShujiaHuang/BaseVar2#installation).

---

## Citation

If you use `basevar motif --fprofile` for cfDNA end-motif tissue
deconvolution, please cite both the BaseVar paper and the F-profile paper:

> Liu, S., Liu, Y., Gu, Y., …, Huang, S. (2024). Utilizing non-invasive prenatal test sequencing data for human genetic investigation. *Cell Genomics* 4(10), 100669. [doi:10.1016/j.xgen.2024.100669](https://doi.org/10.1016/j.xgen.2024.100669)

> Zhou, Z., Ma, R., Huang, X., …, & Sun, K. (2023). Cell-free DNA end motifs characterize cancer fragmentation and detect early cancer. *PNAS* 120(17). [doi:10.1073/pnas.2208138120](https://doi.org/10.1073/pnas.2208138120)

---

**Full Changelog:** https://github.com/ShujiaHuang/BaseVar2/compare/v2.3.1...v2.4.0
