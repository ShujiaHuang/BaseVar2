# BaseVar v2.3.0 — `basevar motif`: cfDNA end-motif counting

**Release date:** 2026-05-26

This release introduces a new top-level subcommand, **`basevar motif`**, for counting cfDNA 5' end-motif (k-mer) frequencies from BAM/CRAM alignments. The methodology is aligned with the canonical Lo lab convention (Y.M. Dennis Lo group), as described in *Jiang P. et al., Cancer Discovery 2020* (PMID: 32111602) and as implemented in the group's open-source reference toolkit [FinaleToolkit](https://github.com/epifluidlab/FinaleToolkit) (Zheng et al., bioRxiv 2024.05.29.596414), with further confirmation in *Mao et al., Cell Genomics* 2026.

The four existing subcommands (`caller`, `pipeline`, `concat`, `subsam`) are **fully backward compatible** — no flags, output formats, or default behaviors have been changed.

---

## Highlights

- **New `basevar motif` subcommand** for cfDNA / NIPT / fragmentomic 5' end-motif analysis.
- **`--from-reference` flag** implements the canonical Lo lab path (motifs fetched from the reference FASTA at each fragment's 5' aligned position).
- **File-level parallelism** with one worker thread per BAM/CRAM, no shared mutable state.
- **Per-sample TSV output** in long format, ready for direct ingestion into ML pipelines as a fixed-shape (4^*k* rows × *N* samples) feature matrix.
- **IUPAC-ambiguity-safe** motif filtering (rejects `N` and all IUPAC ambiguity codes such as `M`/`R`/`W`/`S`/`Y`/`K`/`V`/`H`/`D`/`B`).
- **57 unit-test assertions**, all passing on macOS and Linux CI.

---

## What's New

### 1. `basevar motif` subcommand

- **Inputs:** any number of BAM/CRAM files (each file is treated as one sample).
- **Outputs:** a long-format TSV (`#sample`, `motif`, `count`, `frequency`) holding all samples side-by-side without manual merging, plus a human-readable per-sample summary on stdout (`total / filtered / used / N-motif` counts and Shannon-entropy MDS).
- **Concurrency:** file-level parallelism via `-t/--thread`; the worker count is automatically capped at `min(--thread, number_of_inputs)`.

### 2. `--from-reference` (canonical Lo lab path, recommended)

When this flag is set, the 5' k-mer is fetched from the **reference genome FASTA** at each fragment's aligned 5' position:

- forward read: `ref[chrom, pos, pos + k)`
- reverse read: `reverse_complement(ref[chrom, end - k, end))`

This matches FinaleToolkit's `region_end_motifs()` exactly and is the implementation of "motifs were deduced from the reference genome" described in *Mao et al., Cell Genomics* 2026. Each worker holds its own `ngslib::Fasta` instance to satisfy the underlying thread-safety contract.

> **Strongly recommended for cfDNA / NIPT / fragmentomic analyses.** Requires `-f/--reference`.

### 3. Default (read-based) path — conservative fallback

When `--from-reference` is **not** set, motifs are taken directly from the BAM `SEQ` field:

- forward-mapped read: first *k* bases of `SEQ` directly
- reverse-mapped read: last *k* bases of `SEQ`, reverse-complemented (recovering the original sequencer's 5' bases)

This path does not require a FASTA and is convenient for quick QC on BAM/CRAM inputs from non-cfDNA workflows. **Note:** the read-based path is itself a deviation from the canonical Lo lab convention and is not recommended for publication-quality cfDNA / NIPT studies — use `--from-reference` for those.

### 4. cfDNA-friendly read filters

- `--reads {R1|R2|both}` — which read in a pair contributes a 5' motif. Use `--reads both` so each fragment contributes both ends, matching the Lo lab convention. Single-end data is implicitly treated as R1.
- `--proper-pair` — keep only reads with `BAM_FPROPER_PAIR` set (silently ignored on SE data).
- `--max-insert-size INT` — drop reads with `|isize| > INT`. Lo lab cfDNA pipelines typically use `1000`. Silently ignored on SE data.
- `-q/--mapq` (default `30`), `-l/--length` (default `4`, range 1–10).

### 5. Sample-ID resolution

- Default: read the first `@RG` line's `SM` tag from the BAM/CRAM header.
- `--filename-has-samplename` — force using the filename stem (the part of the basename before the first `.`); also used as automatic fallback when `@RG SM` is missing.

### 6. Per-sample frequency normalization

`frequency = count / used_reads_of_that_sample` — every sample is normalized against its own usable-read denominator, **not** against a pooled total. With `--include-zero` (default ON), all 4^*k* motifs are emitted for every sample, yielding a fixed-shape, ML-ready feature matrix.

---

## Recommended invocation (canonical Lo lab cfDNA / NIPT method)

```bash
basevar motif \
    -o end_motif.lo_lab.tsv \
    --from-reference -f reference.fa \
    --reads both --proper-pair --max-insert-size 1000 \
    -q 30 -l 4 \
    -t 8 \
    -L bamfile.list
```

---

## Documentation correction

The README and source-level docstrings have been corrected to identify `--from-reference` as the **canonical Lo lab path** (earlier wording incorrectly described it as a deviation). The default read-based path is now documented as BaseVar's conservative fallback. Both `FinaleToolkit` and *Mao et al., Cell Genomics* 2026 are now cited as the authoritative sources for the reference-based extraction convention.

---

## Tests

- New file: `tests/io/test_motif_counter.cpp`
- **57 / 57 assertions pass**, covering:
  - reverse-complement and BAM-derived motif extraction primitives
  - end-to-end single-file, multi-file, and multi-sample runs
  - `--proper-pair`, `--max-insert-size`, and `--reads` filters
  - the new `--from-reference` path against `ce.fa.gz` (verified consistent with the read-based path on the bundled test data: `used_reads = 55`, `MDS = 0.689136`)

---

## Build & release

- Project version bumped **2.2.3 → 2.3.0** (`CMakeLists.txt`, `src/version.h.in`, `src/version.h`).
- Linux/macOS pre-built binary download URLs in the README updated to `v2.3.0`.
- `src/main.cpp` registers the `motif` dispatch; `basevar --help` and `basevar motif --help` show the full reference.
- The GitHub Actions static-build workflow publishes `basevar-linux-static` and `basevar-macos-static` as release assets automatically once the `v2.3.0` tag is pushed.

---

## Backward compatibility

- The four existing subcommands (`caller`, `pipeline`, `concat`, `subsam`) are unchanged in semantics, parameters, and output format.
- `motif` is a pure addition; existing scripts continue to work without any modification.
- Without `--from-reference`, the new subcommand has **no new required dependencies** (no FASTA needed).

---

## Downloads

| Platform | Asset | Notes |
| -------- | ----- | ----- |
| Linux (x86_64) | [`basevar-linux-static`](https://github.com/ShujiaHuang/BaseVar2/releases/download/v2.3.0/basevar-linux-static) | Requires **glibc ≥ 2.35** |
| macOS (arm64 / Intel) | [`basevar-macos-static`](https://github.com/ShujiaHuang/BaseVar2/releases/download/v2.3.0/basevar-macos-static) | Requires **macOS 12+** |

For source builds, see the [README](https://github.com/ShujiaHuang/BaseVar2#installation).

---

## Citation

If you use the new `basevar motif` subcommand for cfDNA end-motif analysis, please also cite the methodology references:

> Jiang P. *et al.* (2020). Plasma DNA End-Motif Profiling as a Fragmentomic Marker in Cancer, Pregnancy, and Transplantation. *Cancer Discovery* 10(5), 664–673. [doi:10.1158/2159-8290.CD-19-0622](https://doi.org/10.1158/2159-8290.CD-19-0622)

> Zheng J. *et al.* (2024). FinaleToolkit: Accelerating Cell-Free DNA Fragmentomic Analyses. *bioRxiv* 2024.05.29.596414. [doi:10.1101/2024.05.29.596414](https://doi.org/10.1101/2024.05.29.596414)

And, as always, the BaseVar paper:

> Liu, S., Liu, Y., Gu, Y., …, Huang, S. (2024). Utilizing non-invasive prenatal test sequencing data for human genetic investigation. *Cell Genomics* 4(10), 100669. [doi:10.1016/j.xgen.2024.100669](https://doi.org/10.1016/j.xgen.2024.100669)

---

**Full Changelog:** https://github.com/ShujiaHuang/BaseVar2/compare/v2.2.3...v2.3.0
