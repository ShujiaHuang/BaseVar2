# Getting Started

<cite>
**Referenced Files in This Document**
- [README.md](file://README.md)
- [CMakeLists.txt](file://CMakeLists.txt)
- [src/main.cpp](file://src/main.cpp)
- [src/variant_caller.h](file://src/variant_caller.h)
- [scripts/create_pipeline.py](file://scripts/create_pipeline.py)
- [tests/data/sample_group.info](file://tests/data/sample_group.info)
- [tests/io/make.sh](file://tests/io/make.sh)
- [htslib/INSTALL](file://htslib/INSTALL)
- [src/version.h.in](file://src/version.h.in)
- [build-static/Makefile](file://build-static/Makefile)
</cite>

## Update Summary
**Changes Made**
- Removed Option 3: Manual g++ Compilation section as the manual_install.sh script has been completely removed from the codebase
- Updated installation options to focus on pre-built static binaries and CMake compilation methods only
- Revised troubleshooting guidance to reflect the removal of manual compilation option
- Updated system requirements and verification procedures to match current available installation methods

## Table of Contents
1. [Introduction](#introduction)
2. [Installation Options](#installation-options)
3. [Quick Start Examples](#quick-start-examples)
4. [Command Reference](#command-reference)
5. [Whole-Genome Pipeline](#whole-genome-pipeline)
6. [System Requirements](#system-requirements)
7. [Installation Verification](#installation-verification)
8. [Best Practices](#best-practices)
9. [Troubleshooting Guide](#troubleshooting-guide)
10. [Appendices](#appendices)

## Introduction
BaseVar2 is a C++ tool designed for variant calling from ultra-low-coverage whole-genome sequencing data (typically below 1x), with a focus on non-invasive prenatal testing (NIPT) applications. It efficiently estimates allele frequencies and produces VCF output using maximum likelihood and likelihood ratio models. This guide helps you install BaseVar2 via two methods, verify your installation, and run your first variant calling analyses from BAM files.

**Section sources**
- [README.md:1-17](file://README.md#L1-L17)

## Installation Options

BaseVar2 provides two installation methods, each suited for different use cases and environments:

### Option 1: Pre-built Static Binary (Recommended - No Compilation)
Download pre-built static binaries that require zero runtime dependencies and run on any modern Linux distribution or macOS.

**Linux (x86_64, glibc-free)**:
```bash
wget https://github.com/ShujiaHuang/BaseVar2/releases/download/v2.2.0/basevar-linux-static
chmod +x basevar-linux-static
./basevar-linux-static --help
```

**macOS (arm64 / Intel)**:
```bash
curl -LO https://github.com/ShujiaHuang/BaseVar2/releases/download/v2.2.0/basevar-macos-static
chmod +x basevar-macos-static
./basevar-macos-static --help
```

**Platform Support Matrix**:
| Platform | Binary | Runtime Dependencies | Notes |
|----------|--------|---------------------|-------|
| Linux (x86_64) | `basevar-linux-static` | Zero | Runs on CentOS 7+, Ubuntu 16.04+, Debian 9+ |
| macOS (Intel) | `basevar-macos-static` | Zero | Compatible with macOS versions shipped with BaseVar2 |
| macOS (Apple Silicon) | `basevar-macos-static` | Zero | Native arm64 support |

**Section sources**
- [README.md:19-42](file://README.md#L19-L42)

### Option 2: Compile from Source (CMake Method)
**Prerequisites**: C++17 compiler (GCC 7+ or Apple Clang 10+), CMake ≥ 3.12, and system libraries: zlib, bzip2, xz-utils, libcurl.

**Step 1**: Clone the repository (including htslib submodule)
```bash
git clone --recursive https://github.com/ShujiaHuang/basevar2.git
cd basevar2
```

**Step 2**: Build with CMake (standard dynamic build)
```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

The executable `bin/basevar` will be produced. Verify with:
```bash
./bin/basevar --help
```

**Optional**: Build a static binary locally

**macOS** (requires Homebrew):
```bash
brew install zlib bzip2 xz
cmake -B build-static -DSTATIC_BUILD=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-static
```

**Linux** (fully static via Alpine Docker):
```bash
docker run --rm -v "$PWD:/src" -w /src alpine:3.20 sh -c '
  apk add --no-cache build-base cmake autoconf automake \
    zlib-dev bzip2-dev xz-dev curl-dev openssl-dev &&
  cmake -B build-static -DSTATIC_BUILD=ON -DCMAKE_BUILD_TYPE=Release &&
  cmake --build build-static
'
```

**Section sources**
- [README.md:46-94](file://README.md#L46-L94)
- [CMakeLists.txt:22-62](file://CMakeLists.txt#L22-L62)

## Quick Start Examples

### Minimal Variant Calling
Call variants from a list of BAM files with basic parameters:
```bash
basevar caller \
    -f reference.fasta \
    -o output.vcf.gz \
    -L bamfile.list
```

### Recommended Production Settings
Recommended call with quality filters and sample name optimization:
```bash
basevar caller \
    -f reference.fasta \
    -o output.vcf.gz \
    -Q 20 -q 30 -B 500 -t 24 \
    --filename-has-samplename \
    -L bamfile.list
```

### Region-Specific Calling
Call a specific genomic region:
```bash
basevar caller \
    -f reference.fasta \
    -o chr11_region.vcf.gz \
    -Q 20 -q 30 -B 500 -t 24 \
    --filename-has-samplename \
    -r chr11:5246595-5248428 \
    -L bamfile.list
```

### Multiple Disjoint Regions
Call multiple disjoint regions in one run:
```bash
basevar caller \
    -f reference.fasta \
    -o multi_region.vcf.gz \
    -Q 20 -q 30 -B 500 -t 24 \
    --regions chr11:5246595-5248428,chr17:41197764-41276135 \
    -L bamfile.list
```

### Population-Aware Analysis
Include population groups for per-population allele frequency calculation:
```bash
basevar caller \
    -f reference.fasta \
    -o output.vcf.gz \
    -Q 20 -q 30 -B 500 -t 24 \
    --filename-has-samplename \
    --pop-group sample_group.info \
    -L bamfile.list
```

### Resume Interrupted Runs
Resume an interrupted run safely:
```bash
basevar caller \
    -f reference.fasta \
    -o output.vcf.gz \
    -Q 20 -q 30 -B 500 -t 24 \
    --filename-has-samplename \
    --smart-rerun \
    -L bamfile.list
```

**Section sources**
- [README.md:171-247](file://README.md#L171-L247)

## Command Reference

### BaseVar Commands Overview
```
Usage: basevar <command> [options]

Commands:
  caller    Call variants and estimate allele frequencies
  pipeline  Generate per-region basevar caller commands for whole-genome calling
  concat    Concatenate per-region VCF files into a whole-genome VCF
  subsam    Extract a subset of samples from a VCF file
```

**Section sources**
- [README.md:124-136](file://README.md#L124-L136)

### BaseVar Caller Parameters

**Full Parameter Reference**:
```
About: Call variants and estimate allele frequency by BaseVar.
Usage: basevar caller [options] <-f Fasta> <-o output_file> [-L bam.list] in1.bam [in2.bam ...] ...

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
  --smart-rerun                Skip completed batch files and resume an interrupted run.
  -h, --help                   Show this help message and exit.
```

**Section sources**
- [README.md:138-169](file://README.md#L138-L169)

### BaseVar Pipeline Parameters

**Pipeline Generator Options**:
```
Usage: basevar pipeline [options] <-o outdir> --ref_fai ref.fa.fai -f ref.fa -L bam.list [caller_options]

Required:
  -o, --outdir=DIR             Output directory for VCF files and logs.
  --ref_fai=FILE               Reference FASTA index file (.fai).

Optional:
  -d, --delta=INT              Sub-region size (bp) [2000000]
  -c, --chrom=STR              Comma-separated chromosomes to process
```

**Section sources**
- [README.md:249-267](file://README.md#L249-L267)

### BaseVar Concat Parameters

**Concatenate VCF Files**:
```
Usage: basevar concat [options] <-o output.vcf.gz> [-L vcf.list] in1.vcf.gz [in2.vcf.gz ...]

Required:
  -o, --output=FILE      Output VCF file.

Optional:
  -L, --file-list=FILE   List of input VCF files, one per line.
```

**Section sources**
- [README.md:334-347](file://README.md#L334-L347)

### BaseVar Subsam Parameters

**Extract Samples from VCF**:
```
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

**Section sources**
- [README.md:362-377](file://README.md#L362-L377)

## Whole-Genome Pipeline

For whole-genome variant calling, the pipeline generator creates sub-region commands that can be executed in parallel on a compute cluster or with a job scheduler.

### Pipeline-Specific Options

| Option | Description | Default |
|--------|-------------|---------|
| `-o, --outdir` | Output directory for VCF files and logs | **required** |
| `--ref_fai` | Reference FASTA index file (`.fai`) | **required** |
| `-d, --delta` | Size of each sub-region (bp) | `2000000` |
| `-c, --chrom` | Restrict to comma-separated chromosome(s) | all chromosomes |

**Section sources**
- [README.md:255-264](file://README.md#L255-L264)

### Pipeline Examples

**Generate whole-genome pipeline (all chromosomes, 2 Mb windows)**:
```bash
basevar pipeline \
    -o /path/to/outdir \
    --ref_fai reference.fasta.fai \
    -f reference.fasta \
    -L bamfile.list \
    -Q 20 -q 30 -B 500 -t 4 \
    --filename-has-samplename \
    > basevar_wgs.sh
```

**Generate pipeline for a single chromosome (5 Mb windows)**:
```bash
basevar pipeline \
    -o /path/to/outdir \
    --ref_fai reference.fasta.fai \
    -c chr20 -d 5000000 \
    -f reference.fasta \
    -L bamfile.list \
    -Q 20 -q 30 -B 500 -t 4 \
    --filename-has-samplename \
    > basevar.chr20.sh
```

**Generate pipeline for specific regions (split into 1 Mb windows)**:
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
    > basevar.targets.sh
```

**Run the pipeline**:
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

After all sub-jobs finish, concatenate the per-region VCFs:
```bash
ls /path/to/outdir/*.vcf.gz | sort -V > vcf.list
basevar concat -L vcf.list -o final_output.vcf.gz
```

**Section sources**
- [README.md:268-328](file://README.md#L268-L328)
- [scripts/create_pipeline.py:138-200](file://scripts/create_pipeline.py#L138-L200)

## System Requirements

### Compiler and Standard
- **C++17 or newer**: Required for compilation
- **CMake 3.12+**: Recommended for automated builds

### Required System Libraries
- **Linux/macOS**: zlib, bzip2, lzma, libcurl
- **Non-macOS systems**: OpenSSL and libcrypto for cloud storage support
- **Optional**: libdeflate for faster compression (highly recommended)

### Operating System Support
- **Linux**: CentOS 7+, Ubuntu 16.04+, Debian 9+
- **macOS**: All versions shipped with BaseVar2
- **Other Unix-like systems**: Supported by HTSlib

### Hardware Requirements
- **Memory**: Typically 3–4 GB per thread with `-B 200`
- **Storage**: Depends on cohort size and reference genome
- **CPU**: Multi-core systems recommended for parallel processing

**Section sources**
- [README.md:48](file://README.md#L48)
- [CMakeLists.txt:4-6](file://CMakeLists.txt#L4-L6)
- [CMakeLists.txt:75-79](file://CMakeLists.txt#L75-L79)
- [htslib/INSTALL:31-42](file://htslib/INSTALL#L31-L42)

## Installation Verification

### Verify Pre-built Binary Installation
```bash
./basevar-linux-static --help
# Should display usage information and version details
```

### Verify CMake Installation
```bash
./bin/basevar --help
# Should display the same usage information as pre-built binary
```

### Quick Smoke Test
Use a minimal BAM list and a small region to produce a small VCF output:
```bash
basevar caller \
    -f reference.fasta \
    -o test_output.vcf.gz \
    -Q 20 -q 30 -B 10 \
    --filename-has-samplename \
    -r chr1:1-100000 \
    -L small_bam.list
```

### Expected Output Verification
- VCF file should be created with proper headers
- Log output should indicate successful completion
- Memory usage should be reasonable for the test size

**Section sources**
- [README.md:171-179](file://README.md#L171-L179)
- [src/main.cpp:18-32](file://src/main.cpp#L18-L32)

## Best Practices

### Performance Optimization
- **`-B / --batch-count`**: Controls how many samples are processed per batch. Lower values reduce per-thread memory but increase I/O. For large cohorts (>5000 samples) `-B 500` is a good starting point.
- **`--filename-has-samplename`**: If your BAM files are named `{SampleID}.bam` or `{SampleID}.cram`, always set this flag — it avoids reading every BAM header and can save hours on large cohorts.
- **`--smart-rerun`**: Safe to add on any re-run; the program checks existing batch files and skips completed work.
- **Memory estimation**: `threads × batch_size / 200 × ~3–4 GB`. E.g., 24 threads, `-B 200` → ~72–96 GB total.
- **Output compression**: Always use `.vcf.gz` as the output filename — BaseVar automatically writes bgzipped output when the extension is `.vcf.gz`.

### Quality Control
- **Minimum base quality (`-Q`)**: Set to 20–30 for high-quality results
- **Minimum mapping quality (`-q`)**: Set to 30 for stringent filtering
- **Batch size (`-B`)**: Balance between memory usage and I/O efficiency
- **Threads (`-t`)**: Match to available CPU cores, typically 24 for production

### Data Organization
- **Population groups**: Use `--pop-group` for stratified analysis
- **Region specification**: Use `-r` for targeted calling or `-G` for whole-genome
- **Sample naming**: Ensure consistent naming conventions for optimal performance

**Section sources**
- [README.md:404-411](file://README.md#L404-L411)

## Troubleshooting Guide

### Common Installation Issues

**CMake fails due to missing C++17 support**
- Ensure your compiler supports C++17 and that CMake detects it
- Check compiler version: `g++ --version` or `clang --version`
- Update to GCC 7+ or Apple Clang 10+ if necessary

**HTSlib build errors**
- Install required system libraries: `zlib`, `bzip2`, `lzma`, `libcurl`
- Use autotools to configure and build HTSlib: `autoreconf -i && ./configure && make`
- On macOS, install via Homebrew: `brew install zlib bzip2 xz curl openssl`

**Linker errors on macOS**
- Adjust flags as needed; the build system adds a macOS-specific flag to accommodate platform differences
- For static builds, ensure all required libraries are available via Homebrew

**Network/cloud storage access**
- Enable libcurl and, on non-macOS systems, OpenSSL for S3/GCS support if required
- Static builds on macOS link against system frameworks for network access

### Static Build Issues

**Linux static builds**
- Use Alpine Docker container for maximum portability
- Some features (e.g., DNS via NSS) may not work in fully-static mode
- Building inside Alpine (musl) container is recommended for maximum portability

**macOS static builds**
- Apple Clang does NOT support fully static executables
- Strategy: statically link C++ runtime and third-party archives, keep system frameworks dynamic
- Only system frameworks (libSystem, etc.) remain dynamically linked

### Runtime Issues

**Zero-dependency execution**
- Pre-built static binaries have zero runtime dependencies
- Can run on any modern Linux distribution without installing libraries
- macOS binaries are universal and work on both Intel and Apple Silicon

**Memory usage concerns**
- BaseVar uses significantly less memory than the Python version
- With `-B 200` and one thread, memory usage is typically 3–4 GB per thread
- Monitor memory usage and adjust batch size accordingly

**Verification steps**
- Confirm the executable responds to the caller help option
- Run a small test with a minimal dataset to validate correctness
- Check that output VCF files are properly compressed and indexed

**Section sources**
- [README.md:48](file://README.md#L48)
- [CMakeLists.txt:46-62](file://CMakeLists.txt#L46-L62)
- [CMakeLists.txt:75-79](file://CMakeLists.txt#L75-L79)
- [htslib/INSTALL:252-278](file://htslib/INSTALL#L252-L278)
- [README.md:404-411](file://README.md#L404-L411)

## Appendices

### A. CLI Reference Overview
- **Subcommands**:
  - `caller`: Variant calling and AF estimation
  - `pipeline`: Generate per-region basevar caller commands for whole-genome calling
  - `concat`: Merge VCF files from the same sample sets
  - `subsam`: Extract variants for specified samples from a VCF

**Section sources**
- [src/main.cpp:18-32](file://src/main.cpp#L18-L32)
- [src/variant_caller.h:169-173](file://src/variant_caller.h#L169-L173)

### B. Version Metadata
- Version and author metadata are configured via CMake and templated into the version header
- Current version: 2.2.0
- Author: Shujia Huang

**Section sources**
- [src/version.h.in:1-13](file://src/version.h.in#L1-L13)
- [CMakeLists.txt:8-20](file://CMakeLists.txt#L8-L20)

### C. Test Utilities
- The tests/io/make.sh script demonstrates linking and running unit tests for IO components, which can serve as a reference for verifying your own linking
- Test data includes minimal BAM files and sample group information for validation

**Section sources**
- [tests/io/make.sh:1-24](file://tests/io/make.sh#L1-L24)
- [tests/data/sample_group.info:1-44](file://tests/data/sample_group.info#L1-L44)

### D. Development Workflow
- BaseVar is under active development
- To update to the latest version:
  ```bash
  git pull
  git submodule update --recursive
  cmake -B build -DCMAKE_BUILD_TYPE=Release
  cmake --build build
  ```

**Section sources**
- [README.md:414-424](file://README.md#L414-L424)