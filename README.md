<h3 align="center">
  <a href="https://github.com/ShujiaHuang/BaseVar2">
    <img height="180" src="docs/assets/images/basevar_logo.svg">
  </a>
  <br>A high-performance C++17 tool for variants calling from ultra low-depth WGS data<br>
</h3>

<p align="center">
  <a href="https://github.com/ShujiaHuang/BaseVar2/actions/workflows/build.yml"><img alt="GitHub Actions CI" src="https://github.com/ShujiaHuang/BaseVar2/actions/workflows/build.yml/badge.svg"></a> 
  <a href="https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0"><img src="https://img.shields.io/badge/License-AGPL_v3-blue.svg" alt="License: AGPL v3">
  </a>
</p>

<p align="center">
  <em>
    Download: <a href="https://github.com/ShujiaHuang/BaseVar2/releases/latest/download/basevar-linux-static">Lastest stable version</a>
  </em>
</p>

**BaseVar** is a fast, memory-efficient variant caller for ultra-low-depth (<1×) sequencing data, designed for non-invasive prenatal testing (NIPT) and population-scale genomics. It simultaneously identifies genomic variants and estimates allele frequencies in cohorts of **tens to hundreds of thousands of samples**. And the first version publication in [*Cell Genomics*](https://doi.org/10.1016/j.xgen.2024.100669).

Implemented entirely in **C++17**, BaseVar is **over 100× faster** than the [original Python implementation](https://github.com/ShujiaHuang/basevar/tree/python-version-0.6.1.1) and **5–10× faster** than BaseVar v1, while using substantially less memory. With `-B 200` and a single thread, it typically requires only **3–4 GB** of RAM, compared with **more than 20 GB** for the original Python implementation.

---

## Installation

### Option 1 — Download pre-built binary (Recommended, no compilation needed)

https://github.com/ShujiaHuang/BaseVar2/releases/latest/download/basevar-linux-static

Pre-built static binaries are available on the [GitHub Releases page](https://github.com/ShujiaHuang/BaseVar2/releases).

| Platform | Download | Notes |
| -------- | -------- | ----- |
| Linux (x86_64) | [basevar-linux-static](https://github.com/ShujiaHuang/BaseVar2/releases/latest/download/basevar-linux-static) | Requires **glibc ≥ 2.35** (see below) |
| macOS (arm64 / Intel) | [basevar-macos-static](https://github.com/ShujiaHuang/BaseVar2/releases/latest/download/basevar-macos-static) | Requires **macOS 12+** |

#### System requirements for `basevar-linux-static`

The Linux binary is a partial-static build (built on Ubuntu 22.04 / glibc 2.35). It bundles `libstdc++`, `libgcc`, `htslib`, `zlib`, `bzip2`, `xz` and `openssl` statically; only the system C library (**glibc**) is linked dynamically — and glibc symbol versions are forward-compatible only, so the binary requires the **host glibc to be ≥ 2.35**.

**Quick check on your machine:**

```bash
# If the printed glibc version is >= 2.35, basevar-linux-static will run.
ldd --version | head -1
```

A typical incompatibility error looks like:

```bash
./basevar-linux-static: /lib64/libc.so.6: version `GLIBC_2.35' not found (required by ./basevar-linux-static)
```

If you see this — or you are on CentOS / RHEL / Rocky / AlmaLinux / older Ubuntu / older Debian — please use [Option 2: compile from source](#option-2--compile-from-source). The build is straightforward and takes only a few minutes.

```bash
# Linux
wget https://github.com/ShujiaHuang/BaseVar2/releases/latest/download/basevar-linux-static
chmod +x basevar-linux-static
mv basevar-linux-static basevar
./basevar --help
```

```bash
# macOS
curl -LO https://github.com/ShujiaHuang/BaseVar2/releases/latest/download/basevar-macos-static
chmod +x basevar-macos-static
mv basevar-macos-static basevar
./basevar --help
```

---

### Option 2 — Compile from source

*Requires: C++17 compiler (GCC 7+ or Apple Clang 10+), CMake ≥ 3.12, and system libraries: zlib, bzip2, xz-utils, libcurl.*

#### Step 1 — Clone the repository (including htslib submodule)

```bash
git clone --recursive https://github.com/ShujiaHuang/basevar2.git
cd basevar2
```

> If you forgot `--recursive`, run: `git submodule update --init --recursive`

#### Step 2 — Build with CMake (standard dynamic build)

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

The executable `bin/basevar` will be produced. Verify with:

```bash
./bin/basevar --help
```

#### Step 3 (Optional) — Build a static binary locally

**macOS** (requires Homebrew):

```bash
brew install zlib bzip2 xz
cmake -B build-static -DSTATIC_BUILD=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-static
```

**Linux** (portable static via Ubuntu/glibc — same approach used in CI):

```bash
sudo apt-get install -y build-essential cmake autoconf automake \
    zlib1g-dev libbz2-dev liblzma-dev libssl-dev
cmake -B build-static -DSTATIC_BUILD=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-static
```

This bundles `libstdc++`, `libgcc`, `htslib`, and the compression libs statically; glibc remains dynamic. The resulting binary runs on the build host and on any other host with the same-or-newer glibc.

---

### Option 3 — Manual g++ compilation (fallback)

First, build htslib:

```bash
cd htslib && autoreconf -i && ./configure && make && cd ..
```

Then compile manually:

**Linux:**

```bash
cd bin/
g++ -O3 -fPIC ../src/*.cpp ../src/io/*.cpp ../htslib/libhts.a \
    -I ../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -lssl -lcrypto -o basevar
```

**macOS:**

```bash
cd bin/
g++ -O3 -fPIC ../src/*.cpp ../src/io/*.cpp ../htslib/libhts.a \
    -I ../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -o basevar
```

> **Note:** If you encounter a `test/test_khash.c` compilation error during `make` in htslib, you can safely ignore it — the required `libhts.a` archive is still produced correctly.

---

## Commands overview

```bash
Usage: basevar <command> [options]

Commands:
  caller    Call variants and estimate allele frequencies
  pipeline  Generate per-region `basevar caller` commands for whole-genome calling
  concat    Concatenate per-region VCF files into a whole-genome VCF
  subsam    Extract a subset of samples from a VCF file
  motif     Count cfDNA end-motif (k-mer) frequencies from BAM/CRAM
  dump      Inspect binary batchfile (.bbf) and binary index (.bbi) files
```

---

## `basevar caller` — Variant calling

### Full parameter reference

```bash
About: Call variants and estimate allele frequency by BaseVar.
Usage: basevar caller [options] <-f Fasta> <-o output_file> [-L bam.list/cram.list] in1.bam [in2.bam ...] ...

Required arguments:
  -f, --reference FILE         Input reference FASTA file.
  -o, --output    FILE         Output VCF file (supports .vcf.gz).

Optional options:
  -L, --align-file-list=FILE   BAM/CRAM files list, one file per row.
  -r, --regions=REG[,...]      Restrict calling to these regions (comma-separated).
                               Formats: chr  |  chr:start  |  chr:start-end
                               Example: chr1,chr2:1000000,chr3:5000000-10000000
  -G, --pop-group=FILE         Calculate allele frequency per population group.

  -m, --min-af=float           Prior MAF threshold; positions below this are skipped.
                               Default: min(0.001, 100/num_samples). Usually auto-set.
  -Q, --min-BQ INT             Minimum base quality [10]
  -q, --mapq=INT               Minimum mapping quality [5]
  -B, --batch-count=INT        Samples per batch file [500]
  -t, --thread=INT             Number of threads [14]

  --filename-has-samplename    If BAM/CRAM filenames start with the sample ID
                               (e.g. SampleID.bam), set this flag to skip reading
                               the BAM header for sample names — saves significant time.
  --gt-mode=STRING             Genotype calling mode: 'posterior' (Bayesian
                               posterior-based GT/GQ, dosage INFO) or
                               'legacy' (pure likelihood argmin GT, PL-gap GQ).
                               [posterior]
  --ref-bias=FLOAT             Reference bias coefficient (β) for genotype likelihood.
                               β=0.5 means no bias (default); β<0.5 corrects for
                               alignment reference bias. Typical: 0.45-0.48. [0.5]
  --max-alleles=INT            Maximum active alleles allowed at a site. Sites
                               exceeding this threshold will be skipped. [6]
  --smart-rerun                Skip completed batch files and resume an interrupted run.
                               Validates BBI index footer integrity to detect truncated files.
  -h, --help                   Show this help message and exit.
```

### Usage examples

**Minimal call from a list of BAM files (or CRAM files):**

```bash
basevar caller \
    -f reference.fasta \
    -o output.vcf.gz \
    -L bamfile.list
```

**Recommended call with quality filters and sample name optimization:**

```bash
basevar caller \
    -f reference.fasta \
    -o output.vcf.gz \
    -Q 20 -q 30 -B 500 -t 24 \
    --filename-has-samplename \
    -L bamfile.list
```

**Call a specific region:**

```bash
basevar caller \
    -f reference.fasta \
    -o chr11_region.vcf.gz \
    -Q 20 -q 30 -B 500 -t 24 \
    --filename-has-samplename \
    -r chr11:5246595-5248428 \
    -L bamfile.list
```

**Call multiple disjoint regions in one run:**

```bash
basevar caller \
    -f reference.fasta \
    -o multi_region.vcf.gz \
    -Q 20 -q 30 -B 500 -t 24 \
    --regions chr11:5246595-5248428,chr17:41197764-41276135 \
    -L bamfile.list
```

**Include BAM files directly on the command line:**

```bash
basevar caller \
    -f reference.fasta \
    -o output.vcf.gz \
    -Q 20 -q 30 -B 500 \
    --filename-has-samplename \
    -L bamfile.list \
    sample1.cram sample2.bam sample3.bam
```

**Per-population allele frequency calculation:**

```bash
basevar caller \
    -f reference.fasta \
    -o output.vcf.gz \
    -Q 20 -q 30 -B 500 -t 24 \
    --filename-has-samplename \
    --pop-group sample_group.info \
    -L bamfile.list
```

See the [example `sample_group.info` file](https://github.com/ShujiaHuang/BaseVar2/blob/main/tests/data/sample_group.info) for the expected format.

### Bayesian genotype calling (v2.5.0+)

Since v2.5.0, `basevar caller` uses a **two-pass Bayesian architecture** for genotype calling:

1. **First pass** (unchanged): collect per-sample PL and hard-count AC/AN (backward compatible).
2. **Second pass** (new): use the LRT-estimated population AF as a Hardy-Weinberg prior to compute per-sample genotype posteriors. GT is called as `argmax(posterior)`, GQ is `−10 log₁₀(1 − P(best GT))`, and dosage-based INFO fields (`AC_dosage`, `AN_dosage`, `CAF_dosage`) are emitted.

This is the **default** behavior. To revert to the previous pure-likelihood mode:

```bash
basevar caller --gt-mode legacy -f ref.fa -o out.vcf.gz -L bamfile.list
```

**Backward compatibility**: all existing fields (INFO/AC, INFO/CAF, INFO/AF, FORMAT/PL, FORMAT/AD, FORMAT/DP) retain their original values. The new dosage fields are additive. Use `--gt-mode legacy` to produce output identical to pre-v2.5.0.

| Field | Source | Changed? |
| ----- | ------ | -------- |
| INFO/AF | LRT EM (pooled) | Unchanged |
| FORMAT/PL | `calculatePL()` | Unchanged |
| FORMAT/GT | `argmin(PL)` → `argmax(posterior)` | Changed (revert with `--gt-mode legacy`) |
| FORMAT/GQ | PL gap → posterior Phred | Changed (revert with `--gt-mode legacy`) |
| INFO/AC_dosage | Posterior expected ALT count | **New** |
| INFO/AN_dosage | 2 × N_samples | **New** |
| INFO/CAF_dosage | AC_dosage / AN_dosage | **New** |

**Resume an interrupted run:**

```bash
basevar caller \
    -f reference.fasta \
    -o output.vcf.gz \
    -Q 20 -q 30 -B 500 -t 24 \
    --filename-has-samplename \
    --smart-rerun \
    -L bamfile.list
```

---

## `basevar pipeline` — Whole-genome pipeline generator

For whole-genome variant calling, the `pipeline` subcommand splits the genome into sub-regions and prints one `basevar caller` command per sub-region to stdout. The resulting shell script can be executed sequentially, in parallel with GNU `parallel`, or submitted to a job scheduler (SGE / SLURM / PBS).

> Since **v2.2.0**, this functionality is built directly into the `basevar` binary as a native C++ subcommand. The older `scripts/create_pipeline.py` Python script is still shipped for backward compatibility and produces identical output.

### Pipeline-specific options

| Option | Description | Default |
| ------ | ----------- | ------- |
| `-o, --outdir` | Output directory for VCF files and logs | **required** |
| `--ref_fai`    | Reference FASTA index file (`.fai`) | **required** |
| `-d, --delta`  | Size of each sub-region (bp) | `2000000` |
| `-c, --chrom`  | Restrict to comma-separated chromosome(s) | all chromosomes |

All other options (`-f`, `-L`, `-r`, `-Q`, `-q`, `-B`, `-t`, `--filename-has-samplename`, `--pop-group`, ...) are passed through verbatim to `basevar caller`. This means any new `caller` option works automatically without changes to the pipeline subcommand.

When `-r/--regions` is supplied, those regions are further split into `--delta`-sized windows; otherwise every chromosome in the `.fai` (filtered by `--chrom` if set) is processed.

### Examples

**Generate whole-genome pipeline (all chromosomes, 2 Mb windows):**

```bash
basevar pipeline \
    -o /path/to/outdir \
    --ref_fai reference.fasta.fai \
    -f reference.fasta \
    -L bamfile.list \
    -Q 20 -q 30 -B 500 -t 4 \
    --filename-has-samplename \
    --smart-rerun \
    > basevar_wgs.sh
```

**Generate pipeline for a single chromosome (5 Mb windows):**

```bash
basevar pipeline \
    -o /path/to/outdir \
    --ref_fai reference.fasta.fai \
    -c chr20 -d 5000000 \
    -f reference.fasta \
    -L bamfile.list \
    -Q 20 -q 30 -B 500 -t 4 \
    --filename-has-samplename \
    --smart-rerun \
    > basevar.chr20.sh
```

**Generate pipeline for specific regions:**

```bash
basevar pipeline \
    -o /path/to/outdir \
    --ref_fai reference.fasta.fai \
    -d 1000000 \
    -r chr11:5000000-7000000,chr17 \
    -f reference.fasta \
    -L bamfile.list \
    -Q 20 -q 30 -B 500 -t 4 \
    --filename-has-samplename \
    --smart-rerun \
    > basevar.targets.sh
```

**Run the generated pipeline:**

```bash
# Sequential (local):
bash basevar.chr20.sh

# Parallel with GNU parallel:
cat basevar_wgs.sh | parallel -j 8

# Or submit each line as a cluster job (SGE/SLURM example):
while IFS= read -r cmd; do
    echo "$cmd" | qsub -V -cwd -pe smp 4
done < basevar_wgs.sh
```

After all sub-jobs finish, concatenate the per-region VCFs by using `basevar concat` or `bcftools concat `:

```bash
ls /path/to/outdir/*.vcf.gz | sort -V > vcf.list
basevar concat -L vcf.list -o final_output.vcf.gz
```


---

## `basevar concat` — Concatenate VCF files

Concatenate per-region VCF files produced by `basevar caller` into a single VCF. The files must be provided in the correct genomic order (the tool does not sort positions).

```bash
Usage: basevar concat [options] <-o output.vcf.gz> [-L vcf.list] in1.vcf.gz [in2.vcf.gz ...]

Required:
  -o, --output=FILE      Output VCF file.

Optional:
  -L, --file-list=FILE   List of input VCF files, one per line.
```

**Example:**

```bash
# From a file list
ls outdir/*.vcf.gz | sort -V > vcf.list
basevar concat -L vcf.list -o merged.vcf.gz

# From inline file arguments
basevar concat chr1_1_2000000.vcf.gz chr1_2000001_4000000.vcf.gz -o chr1.vcf.gz
```

> You may also use `bcftools concat --naive-force` as a drop-in alternative.

---

## `basevar motif` — cfDNA end-motif counting

Extract and count the **5' end-motif** (k-mer) of every cfDNA fragment in BAM/CRAM alignments. End-motifs are a well-established feature of cell-free DNA fragmentomics — the canonical method is **Jiang P. *et al.*, *Cancer Discovery* 2020, "Plasma DNA End-Motif Profiling as a Fragmentomic Marker in Cancer, Pregnancy, and Transplantation"** ([PMID: 32111602](https://pubmed.ncbi.nlm.nih.gov/32111602/), [DOI: 10.1158/2159-8290.CD-19-0622](https://doi.org/10.1158/2159-8290.CD-19-0622)) from Prof. Y.M. Dennis Lo's group — and are particularly informative for non-invasive prenatal testing (NIPT) and other low-pass cfDNA studies.

> **Reproducing the Lo lab convention**: the canonical Jiang 2020 protocol (1) extracts the 5' k-mer **from the reference genome** at each fragment's aligned 5' position — *not* from the BAM SEQ field — as confirmed by the Lo group's open-source reference implementation [FinaleToolkit](https://github.com/epifluidlab/FinaleToolkit) (Zheng *et al.*, bioRxiv 2024.05.29.596414) and the group's follow-up Mao *et al.*, *Cell Genomics* 2026; (2) uses **both** fragment ends per pair (5' of R1 *and* 5' of R2); (3) restricts to **properly paired** reads; (4) requires **MAPQ ≥ 30**; (5) drops chimeric/discordant alignments via an insert-size cap (typically `|isize| ≤ 1000`). The recommended invocation is:
>
> ```bash
> basevar motif --from-reference -f reference.fa \
>               --reads both --proper-pair --max-insert-size 1000 \
>               -q 30 -l 4 -o out.tsv  in1.bam in2.bam ...
> ```
>
> The tool's *defaults* deliberately stay conservative (read-derived bases, `--reads R1`, `--proper-pair` off, no insert-size cap, no FASTA required) so that BAM/CRAM inputs from non-cfDNA workflows (e.g. genome-wide variant calling) still produce sensible results out of the box. Users running cfDNA / NIPT / fragmentomic analyses are strongly encouraged to use the canonical invocation above.

**Each input BAM/CRAM is treated as one sample.** Files are processed concurrently (one worker thread per file via `-t/--thread`) and results are emitted side-by-side in a single TSV without merging — every motif row carries its sample ID in the first column.

For each read passing the filters, the 5' end-motif of the underlying cfDNA fragment is extracted via one of two paths:

- **Default (read-based, no FASTA required)** — take the first *k* bases of the BAM SEQ field directly for forward-mapped reads, or the reverse-complement of the last *k* bases for reverse-mapped reads (since BAM stores SEQ in reference-forward orientation, this recovers the original sequencer's 5' bases). This path is BaseVar's conservative fallback; it is itself a deviation from the canonical Lo lab convention.
- **`--from-reference` (canonical Lo lab path, requires `-f`)** — fetch `ref[chrom, pos, pos + k)` for forward-mapped reads, or `reverse_complement(ref[chrom, end - k, end))` for reverse-mapped reads, directly from the FASTA. This matches FinaleToolkit's `region_end_motifs()` exactly and is the recommended path for cfDNA / NIPT / fragmentomic analyses.

Motifs containing **any non-ACGT base** — i.e. `N` *or* an IUPAC ambiguity code such as `M`/`R`/`W`/`S`/`Y`/`K`/`V`/`H`/`D`/`B` — are excluded from the counts (and reported separately in the summary).

### Full parameter reference of `basevar motif`

```bash
About: Count cfDNA end-motif (k-mer) frequencies from BAM/CRAM alignments.
       Each input BAM/CRAM is treated as one sample; per-sample counts are
       emitted in a single TSV (long format).
       Method follows Jiang et al., Cancer Discovery 2020 (PMID: 32111602).
Usage: basevar motif [options] <-o output.tsv> [-L bam.list] in1.bam [in2.bam ...]

Required arguments:
  -o, --output FILE            Output TSV file (sample, motif, count, frequency).

Optional arguments:
  -L, --align-file-list FILE   BAM/CRAM files list, one path per row.
  -f, --reference FILE         Reference FASTA file (required for CRAM input).
  -r, --regions REG[,...]      Restrict counting to these regions, comma-separated.
                               Formats: chr | chr:start | chr:start-end
  -l, --length INT             Motif length k, range [1, 10] [4]
  -q, --mapq INT               Minimum MAPQ to keep a read [30]
  -t, --thread INT             Number of worker threads (one file per thread)
                               [hardware_concurrency]
      --reads {R1|R2|both}     Which read in a pair to use for the 5' end-motif [R1].
                               The default `R1` is conservative and works for any
                               BAM/CRAM. To reproduce the canonical Lo lab cfDNA
                               end-motif method (Jiang et al., Cancer Discovery 2020,
                               PMID: 32111602), use `--reads both` so that EACH
                               fragment contributes TWO end-motifs (5' of R1 AND
                               5' of R2).
      --include-zero           Emit all 4^k motifs (zeros included) in TSV (default ON)
      --no-include-zero        Suppress motifs with zero count in TSV.
      --filename-has-samplename
                               Derive sample IDs from filenames instead of reading the
                               BAM @RG SM tag.  E.g. /path/SampleA.bam -> SampleA.
      --proper-pair            Only count reads flagged as properly paired
                               (BAM_FPROPER_PAIR).  Required by the Lo lab convention
                               for cfDNA analyses; silently ignored for SE data. [off]
      --max-insert-size INT    Discard reads whose |insert size| > INT. 0 = no limit.
                               Lo lab cfDNA pipelines typically use 1000 to drop
                               chimeric / discordantly-mapped reads.
                               Silently ignored for SE data. [0]
      --from-reference         Extract motifs from the REFERENCE genome at each
                               fragment's 5' alignment position (canonical Lo lab
                               cfDNA method). -f/--reference.  When OFF (default),
                               motifs are extracted from the read's own sequenced
                               bases (BaseVar's conservative fallback that does not
                               require a FASTA). [off]
      --fprofile               Compute F-profile decomposition weights using
                               the 6 reference profiles from Zhou et al. (2023)
                               PNAS (doi: 10.1073/pnas.2220982120).  Uses
                               constrained non-negative least squares (CNSLS)
                               on the probability simplex (x >= 0, sum = 1)
                               to decompose the observed 256 4-mer frequencies
                               into tissue contribution proportions.  Only
                               valid for k=4; silently disabled for other
                               values. [off]
      --fprofile-output FILE   Write per-sample F-profile weights to FILE (TSV).
                               Implies --fprofile. Output columns: sample,
                               F-profile I through F-profile VI.
  -h, --help                   Show this help message and exit.

Recommended cfDNA / NIPT invocation (Lo lab convention):
  basevar motif --from-reference -f ref.fa \
                --reads both --proper-pair --max-insert-size 1000 \
                -q 30 -l 4 -o out.tsv  in1.bam in2.bam ...
```

### Sample identification

For each input file the runner resolves a sample ID using the following order:

1. The `SM` value of the first `@RG` line in the BAM/CRAM header (default).
2. The filename stem — everything before the first `.` of the basename — when `--filename-has-samplename` is set, **or** when the `@RG` `SM` tag is missing.

If two inputs resolve to the same sample ID, that ID will appear in two distinct row blocks of the TSV; consider `--filename-has-samplename` together with disambiguating filenames (e.g. `SampleA.run1.bam`, `SampleA.run2.bam`) when this matters.

### Default filters (cfDNA-friendly)

A read contributes to the motif count only if **all** the following hold:

- The read is mapped (not `BAM_FUNMAP`).
- It is not a secondary, supplementary, duplicate, or QC-fail alignment.
- `MAPQ >= --mapq` (default `30`).
- For paired-end data, the read matches the `--reads` policy (default `R1`). Single-end reads are treated as an implicit R1 (accepted under `R1`/`both`, rejected under `R2`).
- The first *k* decoded bases contain only A/C/G/T (no `N` or IUPAC ambiguity codes).

### Output format

The TSV file uses **long format** with one row per (sample, motif) pair:

```tsv
#sample  motif   count   frequency
SampleA  AAAA    18342   0.045312
SampleA  AAAC    7521    0.018584
...
SampleB  AAAA    20131   0.048876
SampleB  AAAC    8033    0.019508
...
```

- **`sample`**: sample ID (see *Sample identification* above).
- **`motif`**: the k-mer (length controlled by `-l`).
- **`count`**: number of reads in that sample whose 5' end-motif equals this k-mer.
- **`frequency`**: `count / used_reads_of_that_sample` — frequencies are computed per-sample, not against the pooled total.

When `--include-zero` is enabled (the default), all 4<sup>*k*</sup> motifs are emitted for **every** sample — yielding a tidy, ML-ready fixed-shape feature matrix. Use `--no-include-zero` to drop zero-count rows.

A human-readable summary is also written to **stdout**, listing the per-sample totals (total / filtered / used / N-motif counts) plus an aggregate row.

### F-profile decomposition (`--fprofile`)

When `--fprofile` is enabled (k=4 only), `basevar motif` decomposes each sample's 256 4-mer frequency vector into the six **"founder" end-motif profiles (F-profiles I–VI)** discovered by **Zhou *et al.*, *PNAS* 2023** ([doi: 10.1073/pnas.2220982120](https://doi.org/10.1073/pnas.2220982120)) via non-negative matrix factorization (NMF) on large cfDNA cohorts. The decomposition uses a **constrained non-negative least squares (CNSLS)** solver — projected gradient descent on the probability simplex {x ≥ 0, Σx = 1} — so the weights are directly constrained to sum to 1.0 without post-hoc normalization.

Each F-profile represents a distinct cfDNA cleavage pattern, likely corresponding to a different nuclease or biological process. The six weights per sample can be interpreted as tissue contribution proportions — useful for fetal fraction estimation in NIPT, tumor fraction in cancer, and general tissue-of-origin deconvolution.

**F-profile output TSV** (written by `--fprofile-output`):

```tsv
sample  F-profile I  F-profile II  F-profile III  F-profile IV  F-profile V  F-profile VI
SampleA  0.234100  0.152300  0.087600  0.201200  0.183400  0.141400
SampleB  0.189200  0.210300  0.045600  0.253400  0.167800  0.133700
```

F-profile weights are also displayed in the stdout summary table alongside MDS.

> **Note:** F-profile decomposition is only defined for k=4 (256 4-mers). When `--fprofile` is used with a different motif length, a warning is printed and the feature is silently disabled.

### Concurrency

File-level parallelism is the granularity used by the motif counter — each worker thread processes a single BAM/CRAM end-to-end, and there is no shared mutable state between workers. The number of workers is automatically capped at `min(--thread, number_of_inputs)`. Pass `-t 1` to force the single-threaded path (useful for deterministic profiling).

### Motif counter examples

**Multi-sample run with 8 worker threads (4-mer end-motif, default filters):**

```bash
basevar motif \
    -o end_motif.4mer.tsv \
    -t 8 \
    -L bamfile.list
```

**Restrict to a target region with stricter MAPQ, two samples:**

```bash
basevar motif \
    -o end_motif.chr11.tsv \
    -r chr11:5246595-5248428 \
    -q 30 \
    sample1.bam sample2.bam
```

**Use both R1 and R2 (each fragment contributes two end-motifs), 5-mer:**

```bash
basevar motif \
    -o end_motif.5mer.tsv \
    -l 5 --reads both \
    -t 4 \
    -L bamfile.list
```

**Reproduce the Lo lab cfDNA end-motif method (Jiang 2020) end-to-end:**

```bash
basevar motif \
    -o end_motif.lo_lab.tsv \
    --from-reference -f reference.fa \
    --reads both --proper-pair --max-insert-size 1000 \
    -q 30 -l 4 \
    -t 8 \
    -L bamfile.list
```

**Compute F-profile decomposition weights (tissue deconvolution):**

```bash
basevar motif \
    -o end_motif.tsv \
    --from-reference -f reference.fa \
    --reads both --proper-pair --max-insert-size 1000 \
    --fprofile --fprofile-output fprofile_weights.tsv \
    -q 30 -l 4 \
    -t 8 \
    -L bamfile.list
```

The `fprofile_weights.tsv` file contains 6 columns (F-profile I–VI) per sample, with values summing to 1.0 — representing the proportional contribution of each founder cleavage profile.

**Read-based motifs (no FASTA required, BaseVar's conservative default):**

```bash
basevar motif \
    -o end_motif.read_based.tsv \
    -t 8 \
    -L bamfile.list
```

With no `--from-reference`, the 5' k-mer is taken directly from the BAM `SEQ` field (the original sequenced bases) rather than the reference. This path does **not** require a FASTA, but it is itself a deviation from the canonical Lo lab convention: the Lo group's reference open-source implementation [FinaleToolkit](https://github.com/epifluidlab/FinaleToolkit) and the Lo group's own follow-up paper (Mao et al., Cell Genomics 2026) both fetch end-motifs from the FASTA. For cfDNA / NIPT / fragmentomic analyses, prefer the recommended invocation above.

**Force filename-derived sample IDs (useful when BAM headers lack `@RG SM`):**

```bash
basevar motif \
    -o end_motif.tsv \
    --filename-has-samplename \
    -t 16 \
    -L bamfile.list
```

**Process CRAM files (reference required):**

```bash
basevar motif \
    -o end_motif.tsv \
    -f reference.fasta \
    -L cram.list
```

### Comparison with FinaleToolkit

[FinaleToolkit](https://github.com/epifluidlab/FinaleToolkit) is the open-source Python reference implementation from the Lo lab for cfDNA fragmentomics. It provides end-motif counting alongside several other features (WPS, DELFI, cleavage profile, fragment length, breakpoint motifs, and MDS). `basevar motif` focuses exclusively on end-motif counting and F-profile decomposition, but does so with significant performance and deployment advantages:

| Feature | `basevar motif` | FinaleToolkit `end-motifs` |
| ------- | :--------------: | :------------------------: |
| **Language** | C++17 (compiled native binary) | Python 3 (pysam + numpy + py2bit + tqdm) |
| **Installation** | Single static binary, zero runtime dependencies | `pip install finaletoolkit` + Python ecosystem |
| **Parallelism** | Thread-level (one thread per BAM, shared memory) | Process-level (`multiprocessing.Pool`, IPC overhead) |
| **Multi-sample batch** | Native: process N BAMs → single unified TSV | One BAM per invocation; scripting needed for batches |
| **F-profile decomposition** | Built-in (`--fprofile`): CNSLS solver (simplex-constrained), zero extra deps | Not included; requires external NMF/NNLS pipeline |
| **Read-based mode** | Yes (default; no FASTA required) | No; always requires reference (.2bit or FASTA) |
| **Reference formats** | FASTA (.fa / .fa.gz) | FASTA + 2bit (.2bit via py2bit) |
| **Output format** | Long-format TSV (sample, motif, count, freq) — ML-ready | 2-column TSV (motif, freq) per sample |
| **`--include-zero` toggle** | Yes (`--no-include-zero` to drop zeros) | No toggle; always emits all 4<sup>k</sup> motifs |
| **IUPAC filtering** | Excludes N *and* all IUPAC ambiguity codes (M/R/W/S/Y/K/V/H/D/B) | Excludes N only |
| **MDS (Shannon entropy)** | Computed and reported in stdout summary | Separate `mds` subcommand (post-hoc from file) |
| **Region restriction** | `-r chr:start-end` (inline, no BED needed) | `region_end_motifs()` API or `interval_end_motifs` (BED) |
| **Fragment length filter** | Not available | `--min-length` / `--max-length` |
| **Strand control** | `--reads {R1\|R2\|both}` | `--no-both-strands` / `--negative-strand` |
| **Proper-pair filter** | `--proper-pair` flag | Hard-coded in `AlignmentWrapper` (always enforced) |
| **Insert-size cap** | `--max-insert-size INT` | Not available |

**Key takeaways:**

- **Speed and memory**: `basevar motif` is a compiled C++17 program using thread-level parallelism with shared memory — it typically runs **5–10× faster** than FinaleToolkit's Python multiprocessing approach on the same data, with substantially lower memory overhead (no Python interpreter, no numpy arrays, no IPC serialization).

- **Deployment**: BaseVar ships as a single statically-linked binary. There is no Python environment to manage, no `pip install`, no version conflicts between pysam/numpy/py2bit. Drop the binary on any Linux or macOS machine and run.

- **Batch processing**: `basevar motif` natively processes multiple BAM/CRAM files in one invocation and emits a single unified TSV — ideal for cohort-scale analyses. With FinaleToolkit, each sample requires a separate `end-motifs` invocation and post-hoc concatenation.

- **F-profile decomposition**: BaseVar includes a built-in CNSLS solver (projected gradient descent on the probability simplex) that decomposes each sample's 256 4-mer frequency vector into the six Zhou *et al.* (2023) founder profiles, producing tissue contribution proportions directly — weights are constrained to sum to 1.0 at optimization time rather than normalized post-hoc. FinaleToolkit ships the reference data (`end_motif_f_profiles.tsv`) but does not include a decomposition solver — users must implement their own NMF/NNLS pipeline.

- **Scope**: FinaleToolkit offers a broader suite of cfDNA fragmentomics features (WPS, DELFI, cleavage profile, fragment length distributions, breakpoint motifs) that `basevar motif` does not cover. If you need WPS or DELFI alongside end-motifs, FinaleToolkit remains the appropriate tool. For high-throughput end-motif counting and F-profile decomposition at cohort scale, `basevar motif` is purpose-built.

---

## `basevar subsam` — Extract samples from VCF

Extract a subset of samples from a BaseVar VCF and output a new VCF with recalculated INFO fields (AC/AN/AF).

```bash
Usage: basevar subsam [options] -i <input.vcf[.gz]> -o <output.vcf[.gz]> [-s <samplelist>]

Options:
  -i, --input FILE      Input VCF/BCF file (required).
  -o, --output FILE     Output VCF/BCF file (required).
  -s, --sample FILE     File with sample names to keep (one per line).
  -O, --output-type     v: VCF | z: bgzipped VCF | b: BCF | u: uncompressed BCF
                        Default: inferred from output filename extension.
  --no-update-info      Do not recalculate AC/AN/AF INFO fields.
  --keep-all-site       Retain sites that become reference-only after subsetting.
```

**Examples:**

```bash
# Extract samples listed in a file
basevar subsam \
    -i full_cohort.vcf.gz \
    -o subset.vcf.gz \
    -s sample_names.txt

# Extract two specific samples, output as plain VCF
basevar subsam \
    -i full_cohort.vcf.gz \
    -o subset.vcf \
    -O v \
    SampleA SampleB

# Keep all sites (including ref-only after subsetting) and skip INFO update
basevar subsam \
    -i full_cohort.vcf.gz \
    -o subset.vcf.gz \
    -s sample_names.txt \
    --keep-all-site --no-update-info
```

---

## `basevar dump` — Inspect binary batchfile files

Inspect the intermediate binary batchfile (`.bbf`) and binary index (`.bbi`) files produced by `basevar caller`. This is useful for debugging, verifying file integrity, or examining per-sample read data at specific genomic positions.

```bash
Usage: basevar dump <file> [options]

Options:
  --header          Show only file header (sample IDs) and exit (.bbf only)
  -r, --region      Dump records in region chr:start-end (.bbf only)
  -n, --limit INT   Dump at most INT records (.bbf only, default: all)
  -v, --verbose     Show per-sample details for each record (.bbf only)
  --entries         Show all index entries (.bbi only)
  -h, --help        Show this help message
```

**Examples:**

```bash
# Inspect .bbi index summary (magic, version, entries, position range, footer integrity)
basevar dump sample.bbf.bbi

# List all index entries
basevar dump sample.bbf.bbi --entries

# Show .bbf header (sample IDs) only
basevar dump sample.bbf --header

# Show all positions in compact summary format
basevar dump sample.bbf

# Show per-sample details at a specific position
basevar dump sample.bbf -r chr11:5246595-5246595 -v

# Dump first 20 positions with full per-sample read data
basevar dump sample.bbf -n 20 -v
```

The file type is auto-detected from the extension (`.bbf` or `.bbi`). If the extension is ambiguous, the tool reads the magic bytes to determine the format.

---

## Tips and best practices

- **`-B / --batch-count`**: Controls how many samples are processed per batch. Lower values reduce per-thread memory but increase I/O. For large cohorts (>5000 samples) `-B 500` is a good starting point.
- **`--filename-has-samplename`**: If your BAM files are named `{SampleID}.bam` or `{SampleID}.cram`, always set this flag — it avoids reading every BAM header and can save hours on large cohorts.
- **`--smart-rerun`**: Safe to add on any re-run; the program checks existing batch files and validates the `.bbi` index footer integrity marker to detect truncated files before resuming.
- **Memory estimation**: `threads × batch_size / 200 × ~3–4 GB`. E.g., 24 threads, `-B 200` → ~72–96 GB total.
- **Binary batchfile format**: Since v2.5.3, BaseVar uses a BGZF-compressed binary batchfile format (`.bbf` + `.bbi` index) instead of the text-based format. This provides ~1.4× speedup over text batchfiles by eliminating join/split overhead and reducing BGZF I/O calls. Use `basevar dump` to inspect these intermediate files.
- **Sparse binary index (`.bbi`)**: The binary index skips positions where all samples have zero depth, reducing index size from ~GB to ~MB for typical cohorts. This also enables faster random access and reduces memory usage during variant calling.
- **Shared index loading**: When using multiple threads, BBI indexes are loaded once and shared across all threads via `std::cref`, eliminating redundant I/O and significantly reducing startup time for large cohorts.
- **Output compression**: Always use `.vcf.gz` as the output filename — BaseVar automatically writes bgzipped output when the extension is `.vcf.gz`.

## Citation

If you use `BaseVar` in your research work, please cite the following paper:

> Liu, S., Liu, Y., Gu, Y., Lin, X., Zhu, H., Liu, H., Xu, Z., Cheng, S., Lan, X., Li, L., Huang, M., Li, H., Nielsen, R., Davies, RW., Albrechtsen, A., Chen, GB., Qiu, X., Jin, X., **Huang, S.**, (2024). Utilizing non-invasive prenatal test sequencing data for human genetic investigation. *Cell Genomics* 4(10), 100669. [doi:10.1016/j.xgen.2024.100669](https://www.cell.com/cell-genomics/fulltext/S2666-979X(24)00288-X) (basevar version 1.3+)