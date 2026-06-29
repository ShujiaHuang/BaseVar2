# BaseVar v2.3.1 — `basevar fetalfrac`: NIPT cfDNA male fetal-fraction estimation (hidden subcommand)

**Release date:** 2026-05-27

This is a **minor** release. It adds one new (intentionally hidden) top-level
subcommand, **`basevar fetalfrac`**, for estimating the male fetal fraction
(FF) of NIPT / lpWGS cfDNA samples directly from BAM/CRAM alignments using
the chrY/autosome read-count method. No public-facing flag, output column,
or default behavior of the five existing subcommands (`caller`, `pipeline`,
`concat`, `subsam`, `motif`) has been changed.

The subcommand is **deliberately not advertised** in `basevar --help`; it is
shipped for internal / clinical-pipeline use and is invoked by direct call
(`basevar fetalfrac ...`).

---

## Highlights

- **New `basevar fetalfrac` subcommand** for NIPT cfDNA male fetal-fraction
  estimation from chrY/autosome read counts, with proper PAR1/PAR2
  exclusion on chrY.
- **Automatic geometric `--scale`** computation when the user does not pin
  it: derived from the `-B` mappability BED if supplied, otherwise from the
  canonical GRCh38 / GRCh37 chromosome lengths.
- **Built-in PAR coordinates** for `hg38` and `hg19`, plus `--build none`
  and `--par-bed` overrides for custom assemblies.
- **File-level parallelism**, one worker thread per BAM/CRAM, no shared
  mutable state — same concurrency model as `basevar motif`.
- **Per-sample TSV output** (sex / FF / counters / parameters used) plus a
  human-readable stdout summary.
- **86 / 86 unit-test assertions pass** on macOS and Linux CI (10 new
  assertions over v2.3.0).

---

## What's New

### 1. `basevar fetalfrac` subcommand (hidden)

- **Inputs:** any number of BAM/CRAM files (each file = one sample).
- **Outputs:** a long-format TSV with one row per sample, holding sex
  (`MALE` / `FEMALE` / `UNDETERMINED`), `fetal_fraction`, the raw `y_ratio`,
  `valid_autosomal` / `valid_y` / `y_par_excluded` counters, the
  `total_reads` / `filtered_reads` audit pair, and the **exact**
  `scale` / `noise` / `male_threshold` parameters used for that sample so
  results are fully reproducible from the TSV alone.
- **Concurrency:** file-level parallelism via `-t/--thread`; worker count is
  automatically capped at `min(--thread, number_of_inputs)`.
- **Discoverability:** the subcommand is not listed in `basevar --help`
  on purpose; users have to know the name to invoke it.

### 2. Genome-wide mappability mask `-B`

When `-B/--mappability-bed` is supplied, **only** reads whose 5'-mapped
position falls inside an interval are counted. The mask is applied
**uniformly to BOTH the chrY numerator AND the chr1..22 denominator**,
*before* the autosome / Y split, so the two sides of the ratio remain
geometrically comparable.

The BED reader has been hardened to follow the rest of BaseVar:

- BED tokenization now goes through the project-standard `ngslib::split`
  tokenizer (consistent with `pipeline.cpp`, `variant_caller.cpp`,
  `motif_counter.cpp`).
- UCSC `track` / `browser` header lines are now skipped via proper prefix
  matching (a previous draft would have dropped any record whose `chrom`
  field happened to start with `t`).
- CRLF-terminated BED files are accepted (trailing `\r` is stripped).
- Numeric `start` / `end` are parsed with explicit `std::stoll` + try/catch
  so malformed lines are skipped instead of silently mis-parsed.

### 3. `--scale` is now auto-computed (geometric default)

Previously `--scale` defaulted to a fixed `100.0`, which forced every
operator to manually re-calibrate when they switched on `-B` or moved
between assemblies. v2.3.1 introduces an automatic geometric default,
selected at startup from the data the user actually provided:

| Situation | Auto-computed `scale` |
| --- | --- |
| `-B BED` supplied | `2 × Σ_BED_autosome_bp / (Σ_BED_chrY_bp − chrY ∩ PAR)` |
| no `-B`, `--build hg38` | `2 × 2,875,001,522 / 54,126,423 ≈ 106.23` |
| no `-B`, `--build hg19` | `2 × 2,881,033,286 / 56,404,529 ≈ 102.16` |
| no `-B`, `--build none` | falls back to `100.0` with a `[WARN]` line |
| `--scale FLOAT` given | the user value is honored verbatim, no auto-compute |

The chosen value is announced on stderr with `[INFO]` (or `[WARN]` on
fallback) and is also tagged in the per-run summary as either
`(auto, geometric)` or `(user-supplied)`, so downstream logs unambiguously
record which path was taken.

> The auto value is only a uniform-coverage geometric approximation. For
> clinical reporting, you should still recalibrate against male reference
> samples processed with the **same** `-B` / `--build` and pin the result
> with `--scale FLOAT`. The auto value is intended to remove the silent
> miscalibration footgun, not to replace per-cohort QC.

### 4. PAR (pseudoautosomal region) handling

chrY reads that fall inside the pseudoautosomal regions are double-counted
on chrX and would inflate FF if naively summed. v2.3.1 ships built-in PAR
coordinates for the common builds and exposes them via `--build` /
`--par-bed`:

| `--build` | PAR1 | PAR2 |
| --- | --- | --- |
| `hg38` (default) | `chrY:10001-2781479` | `chrY:56887903-57217415` |
| `hg19` | `chrY:10001-2649520` | `chrY:59034050-59363566` |
| `none` | (none — disable PAR exclusion) | |

PAR-excluded chrY reads are reported separately in the
`y_par_excluded` TSV column for full traceability.

### 5. Sex calling and counters

- A sample is called **MALE** when `y_ratio = valid_y / valid_autosomal`
  exceeds `--male-threshold` (default `1e-4`), and FF is reported.
- A sample is called **FEMALE** when `y_ratio` is below the threshold; FF
  is intentionally **not** reported (lpWGS Y-count is not a reliable FF
  estimator for female fetuses — use a SNP-based assay).
- A sample is called **UNDETERMINED** when the autosomal denominator is
  zero (e.g. a chrY-only BED was supplied, or all reads were filtered).
- All three calls are written with full counters, so the TSV is enough to
  audit any downstream filtering decision.

### 6. cfDNA-friendly read filters

Re-uses the same filter vocabulary as `basevar motif` so users do not have
to learn a second set of flags:

- `-q/--mapq` — minimum MAPQ, default `30`.
- `--proper-pair` — only count `BAM_FPROPER_PAIR` reads (silently ignored
  for SE data).
- `--max-insert-size INT` — drop reads with `|TLEN| > INT`. Lo-lab cfDNA
  pipelines typically use `1000`; default `0` (no limit).
- `--filename-has-samplename` — derive the sample ID from the filename
  instead of `@RG SM`, matching the convention used by some NIPT vendors.

---

## Recommended invocation

```bash
# Single BAM, GRCh38, default scale (auto, geometric).
basevar fetalfrac -o ff.tsv  sample.bam

# Cohort run with a mappability mask -- recommended for chrY because of
# its repeat / palindrome content.  --scale auto-computes from the BED.
basevar fetalfrac -t 8 -B chrY.unique.bed --build hg38 \
                  -o cohort.ff.tsv -L bam.list

# Pin --scale to a value calibrated on male reference samples.  When
# --scale is given, the auto-compute path is skipped.
basevar fetalfrac -t 8 -B chrY.unique.bed --build hg38 \
                  --scale 105.3 --noise 1e-5 \
                  -o cohort.ff.tsv -L bam.list

# CRAM input.  -f is mandatory.
basevar fetalfrac -f hg38.fa -o ff.tsv  sample.cram
```

---

## Tests

- New file: `tests/io/test_fetal_fraction.cpp`
- **86 / 86 assertions pass** (was 76 in v2.3.0; +10 new), covering:
  - the `BedIntervalIndex` reader (BED parsing via `ngslib::split`,
    `track` / `browser` skipping, CRLF tolerance, numeric guards).
  - the new `BedIntervalIndex::total_length(predicate)` accessor.
  - the auto-compute `--scale` path for hg38 (`~106`), hg19 (`~102`),
    `--build none` (fallback `100`), and `-B BED`-driven scenarios
    (including the chrY-entirely-inside-PAR1 fallback).
  - explicit `--scale FLOAT` is honored verbatim, no auto-compute.
  - end-to-end MALE / FEMALE / UNDETERMINED sex calls and the full TSV
    schema (`#sample`, `sex`, `fetal_fraction`, `y_ratio`,
    `valid_autosomal`, `valid_y`, `y_par_excluded`, `total_reads`,
    `filtered_reads`, `scale`, `noise`, `male_threshold`).

---

## Build & release

- Project version bumped **2.3.0 → 2.3.1** (`CMakeLists.txt`).
- README pre-built binary download URLs updated to `v2.3.1`.
- `src/main.cpp` registers the (hidden) `fetalfrac` dispatch.
- `tests/io/Makefile` now builds and links `test_fetal_fraction` alongside
  the existing `test_fasta` / `test_motif_counter` targets.
- The GitHub Actions static-build workflow publishes
  `basevar-linux-static` and `basevar-macos-static` as release assets
  automatically once the `v2.3.1` tag is pushed.

---

## Backward compatibility

- The five existing subcommands (`caller`, `pipeline`, `concat`, `subsam`,
  `motif`) are **unchanged** in semantics, parameters, and output format.
- `fetalfrac` is a pure addition and is **hidden from `basevar --help`**,
  so it cannot accidentally break any user-facing script that enumerates
  subcommands.
- No new required runtime dependencies: BED / PAR handling is in-tree;
  FASTA is only required when the input is CRAM (same rule as
  `basevar motif`).

---

## Downloads

| Platform | Asset | Notes |
| -------- | ----- | ----- |
| Linux (x86_64) | [`basevar-linux-static`](https://github.com/ShujiaHuang/BaseVar2/releases/download/v2.3.1/basevar-linux-static) | Requires **glibc ≥ 2.35** |
| macOS (arm64 / Intel) | [`basevar-macos-static`](https://github.com/ShujiaHuang/BaseVar2/releases/download/v2.3.1/basevar-macos-static) | Requires **macOS 12+** |

For source builds, see the [README](https://github.com/ShujiaHuang/BaseVar2#installation).

---

## Citation

If you use `basevar fetalfrac` for cfDNA / NIPT male-fetal-fraction
estimation, please cite the BaseVar paper:

> Liu, S., Liu, Y., Gu, Y., …, Huang, S. (2024). Utilizing non-invasive prenatal test sequencing data for human genetic investigation. *Cell Genomics* 4(10), 100669. [doi:10.1016/j.xgen.2024.100669](https://doi.org/10.1016/j.xgen.2024.100669)

---

**Full Changelog:** https://github.com/ShujiaHuang/BaseVar2/compare/v2.3.0...v2.3.1
