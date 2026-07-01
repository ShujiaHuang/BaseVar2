# Synthetic Test Data for BaseVar Caller

## 1. Overview

This directory contains a complete synthetic test data framework for validating basevar caller. It includes:

- A mini reference genome (4 chromosomes, ~5 kb total)
- 50 simulated BAM/CRAM files organized into 3 population groups
- 11 ground-truth variant sites (SNP, indel, multi-allelic)
- Per-sample genotype truth table
- Automated evaluation pipeline (`evaluate.py`)
- Systematic test suite covering all caller parameters (`run_tests.sh`)

Generated: 2026-07-01 21:29:58
Random seed: 42 (fully reproducible)

## 2. Directory Structure

```
synthetic/
├── config.py              # Simulation configuration (variants, samples, coverage)
├── reference.py           # Reference genome generator
├── reads_simulator.py     # Read simulator (variant-aware)
├── writer.py              # BAM/CRAM file writer
├── auxiliary.py           # Auxiliary file generator (truth table, sample lists)
├── generate.py            # Master script: orchestrates data generation
├── evaluate.py            # Evaluation: VCF vs ground truth comparison
├── run_tests.sh           # Systematic test suite (25 scenarios)
│
├── ref/
│   ├── mini_ref.fa        # Reference FASTA
│   └── mini_ref.fa.fai    # FASTA index
│
├── bam/
│   ├── sampleA01.bam      # BAM files (50 total)
│   └── sampleA01.bam.bai  # BAM indices
│
├── cram/
│   ├── sampleA01.cram     # CRAM files (50 total)
│   └── sampleA01.cram.crai # CRAM indices
│
├── samples.list           # Absolute paths to all BAM files (for basevar -L)
├── samples_group.info     # Sample-to-group mapping (for basevar -G)
├── ground_truth_variants.tsv  # Complete truth table
└── README.md              # This file
```

## 3. Reference Genome Design

The reference is a small synthetic genome designed for fast testing:

| Chromosome | Length | Purpose |
|:----------:|-------:|:--------|
| chr1 | 2000 bp | Main test chromosome (carries most variants) |
| chr2 | 1500 bp | Multi-chromosome region test |
| chrX | 1000 bp | Sex chromosome test |
| chrY | 500 bp | Reserved (no variants, for empty-chr testing) |

**Total:** 5000 bp, GC content ~45%

## 4. Ground Truth Variants

11 variant sites spanning multiple variant types:

| ID | Chrom | Pos | Ref | Alt | Type | Description |
|:--:|:-----:|----:|:---:|:---:|:----:|:------------|
| v1 | chr1 | 200 | C | A | SNP | Common in all groups |
| v2 | chr1 | 400 | C | A | SNP | Common in all groups |
| v3 | chr1 | 600 | G | A | SNP | GroupB-specific |
| v4 | chr1 | 800 | G | A | SNP | GroupA,GroupC enriched |
| v5 | chr1 | 1000 | C | A,G | multi-allelic | Common in all groups |
| v6 | chr1 | 1200 | TG | T | del | Common in all groups |
| v7 | chr1 | 1400 | C | CT | ins | Common in all groups |
| v8 | chr1 | 1600 | ATA | A | del | Common in all groups |
| v9 | chr2 | 300 | T | A | SNP | Common in all groups |
| v10 | chr2 | 700 | T | A | SNP | GroupA,GroupB enriched |
| v11 | chrX | 200 | C | A | SNP | Common in all groups |

### Per-group Allele Frequencies

| ID | Type | GroupA_AF | GroupB_AF | GroupC_AF |
|:--:|:----:|--------:|--------:|--------:|
| v1 | SNP | 0.50 | 0.50 | 0.50 |
| v2 | SNP | 0.10 | 0.50 | 0.80 |
| v3 | SNP | 0.00 | 0.30 | 0.00 |
| v4 | SNP | 0.30 | 0.00 | 0.60 |
| v5 | multi-allelic | 0.30,0.10 | 0.20,0.05 | 0.40,0.10 |
| v6 | del | 0.20 | 0.40 | 0.10 |
| v7 | ins | 0.20 | 0.20 | 0.50 |
| v8 | del | 0.10 | 0.30 | 0.20 |
| v9 | SNP | 0.40 | 0.40 | 0.40 |
| v10 | SNP | 0.20 | 0.60 | 0.00 |
| v11 | SNP | 0.30 | 0.50 | 0.20 |

## 5. Sample Design

50 samples organized into 3 population groups:

| Group | Samples | Count | Purpose |
|:-----:|:--------|------:|:--------|
| GroupA | sampleA01 ~ sampleA20 | 20 | Population with distinct AF profile |
| GroupB | sampleB01 ~ sampleB20 | 20 | Population with distinct AF profile |
| GroupC | sampleC01 ~ sampleC10 | 10 | Smaller population group |

### Coverage Distribution

| Coverage | Samples | Count | Purpose |
|:--------:|:--------|------:|:--------|
| 30x | sampleA01 | 1 | High-depth edge case |
| 8x | All others except below | 45 | Standard calling depth |
| 1x | sampleC01, sampleC02, sampleC03, sampleC04 | 4 | Low-depth edge case |

## 6. Auxiliary Files

**`samples.list`**
One BAM path per line (50 lines). Used with `basevar -L`.
Contains relative paths (relative to this directory); regenerate with `generate.py` if structure changes.

**`samples_group.info`**
Tab-separated: `<sample_id>\t<group_name>`. Used with `basevar -G`.

```
sampleA01    GroupA
...
```

**`ground_truth_variants.tsv`**
Tab-separated truth table with header row.

| Column group | Contents |
|:-------------|:---------|
| `variant_id, chrom, pos, ref, alt, type` | Variant definition |
| `AF_<group>...` | Per-group allele frequency |
| `<sample>_GT...` | Per-sample genotype (0/0, 0/1, 1/1, 0/2, 1/2, etc.) |

Genotypes assigned via Hardy-Weinberg equilibrium from group AF.

## 7. Usage: Testing BaseVar Caller

All commands assume you are in the `synthetic/` directory.

### 7.1 Single-sample calling

```bash
basevar caller -f ref/mini_ref.fa -o out.vcf bam/sampleA01.bam
```

### 7.2 Multi-sample calling (all samples)

```bash
basevar caller -f ref/mini_ref.fa -o out.vcf -L samples.list
```

### 7.3 With population stratification

```bash
basevar caller -f ref/mini_ref.fa -o out.vcf -L samples.list -G samples_group.info
```

### 7.4 Region-specific calling

```bash
basevar caller -f ref/mini_ref.fa -o out.vcf -L samples.list -r chr1:100-800
```

### 7.5 Posterior genotype mode

```bash
basevar caller -f ref/mini_ref.fa -o out.vcf -L samples.list --gt-mode
```

### 7.6 CRAM input

```bash
basevar caller -f ref/mini_ref.fa -o out.vcf cram/sampleA01.cram
```

### 7.7 Smart rerun (resume from cache)

```bash
basevar caller -f ref/mini_ref.fa -o out.vcf -L samples.list --smart-rerun
```

### 7.8 Custom quality filters

```bash
basevar caller -f ref/mini_ref.fa -o out.vcf -L samples.list -Q 20 -q 10 -m 0.01
```

## 8. Usage: Automated Test Suite

`run_tests.sh` runs 25 test scenarios covering all caller parameters:

```bash
bash run_tests.sh                  # Run all tests (PASS/FAIL only)
bash run_tests.sh --evaluate       # Run tests + evaluate VCF vs truth
```

| Tests | Category | Description |
|:------|:---------|:------------|
| T01–T05 | Input modes | Single BAM, multi BAM, `-L` list, CRAM, `-L`+CRAM |
| T06–T08 | Region calling | `-r` chr1, `-r` chr1:100-800, `-r` chrX |
| T09–T10 | Population stratification | `-G`, `-G` + `-r` |
| T11–T12 | GT mode | Posterior vs legacy |
| T13–T14 | Smart rerun | Fresh + resume |
| T15–T16 | Quality filters | `-Q`, `-q` |
| T17–T20 | Other parameters | `--filename-has-samplename`, `--ref-bias`, `--max-alleles`, `-B` |
| T21–T25 | Combo workflows | Full combo, concat, multi-region, etc. |

## 9. Usage: Evaluation Pipeline

`evaluate.py` compares a VCF against `ground_truth_variants.tsv`:

```bash
python3 evaluate.py --vcf out.vcf --truth ground_truth_variants.tsv
```

### Output Metrics

- **Sensitivity** (recall): fraction of truth variants found in VCF
- **Precision**: fraction of VCF variants matching truth
- **GT concordance**: per-sample genotype agreement rate
- **AF accuracy**: MAE / RMSE / Pearson *r* between called AF and truth AF
- **Per-group AF breakdown**

> Indel matching tolerates VCF left-alignment position offsets (search radius = max allele length + 1).

## 10. Regenerating Data

If you modify `config.py` (e.g., change sample count, variants, coverage):

```bash
python3 generate.py                  # Generate BAM + CRAM
python3 generate.py --skip-cram      # BAM only (faster)
```

This will regenerate all files: `ref/`, `bam/`, `cram/`, `samples.list`, `samples_group.info`, `ground_truth_variants.tsv`.

After regeneration, run the test suite:

```bash
bash run_tests.sh --evaluate
```
