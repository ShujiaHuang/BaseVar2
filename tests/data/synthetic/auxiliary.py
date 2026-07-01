"""
Auxiliary file generator.

Produces:
  - samples.list: BAM file list for -L parameter
  - samples_group.info: population group assignments for -G parameter
  - ground_truth_variants.tsv: complete truth table with per-sample genotypes
"""

import os

from config import (
    CHROMOSOMES, SAMPLE_GROUPS, LOW_COVERAGE_SAMPLES, LOW_COVERAGE,
    HIGH_COVERAGE_SAMPLES, HIGH_COVERAGE, NORMAL_COVERAGE,
    get_all_samples, get_sample_group, get_target_coverage,
)


def write_samples_list(bam_paths, output_path):
    """
    Write a samples.list file (one BAM path per line).

    Parameters
    ----------
    bam_paths : dict
        {sample_id: bam_path}
    output_path : str
    """
    all_samples = get_all_samples()
    with open(output_path, "w") as fh:
        for sample_id in all_samples:
            if sample_id in bam_paths:
                fh.write(bam_paths[sample_id] + "\n")
    print(f"  [auxiliary] Written {len(all_samples)} entries to {output_path}")


def write_samples_group(output_path):
    """
    Write a samples_group.info file (two columns: sampleID groupName).

    Parameters
    ----------
    output_path : str
    """
    with open(output_path, "w") as fh:
        for group, samples in SAMPLE_GROUPS.items():
            for sample_id in samples:
                fh.write(f"{sample_id}\t{group}\n")

    total = sum(len(s) for s in SAMPLE_GROUPS.values())
    print(f"  [auxiliary] Written {total} samples in {len(SAMPLE_GROUPS)} groups to {output_path}")


def write_ground_truth(sample_genotypes, variants, output_path):
    """
    Write ground_truth_variants.tsv with complete truth table.

    Format:
      # variant_id  chrom  pos  ref  alt  type  AF_GroupA  AF_GroupB  AF_GroupC  sampleA01_GT  ...  sampleC05_GT
      v1  chr1  200  A  G  SNP  0.50  0.50  0.50  0/1  1/1  0/0  ...  0/1

    Parameters
    ----------
    sample_genotypes : dict
        {sample_id: {variant_id: genotype_string}}
    variants : list
        Adapted variant definitions (after adapt_variants_to_ref).
    output_path : str
    """
    all_samples = get_all_samples()
    groups = list(SAMPLE_GROUPS.keys())

    with open(output_path, "w") as fh:
        # Header line
        header_parts = [
            "#variant_id", "chrom", "pos", "ref", "alt", "type"
        ]
        for group in groups:
            header_parts.append(f"AF_{group}")
        for sample_id in all_samples:
            header_parts.append(f"{sample_id}_GT")
        fh.write("\t".join(header_parts) + "\n")

        # Data lines
        for v in variants:
            vid = v["id"]
            row = [
                vid,
                v["chrom"],
                str(v["pos"]),
                v["ref"],
                v["alt"],
                v["type"],
            ]

            # AF per group
            for group in groups:
                af = v["af"][group]
                if isinstance(af, tuple):
                    # Multi-allelic: format as "af1,af2"
                    af_str = ",".join(f"{a:.2f}" for a in af)
                else:
                    af_str = f"{af:.2f}"
                row.append(af_str)

            # Per-sample genotypes
            for sample_id in all_samples:
                gt = sample_genotypes[sample_id].get(vid, "0/0")
                row.append(gt)

            fh.write("\t".join(row) + "\n")

    print(f"  [auxiliary] Written ground truth for {len(variants)} variants x {len(all_samples)} samples to {output_path}")


def write_readme(output_path, seed, bam_count, variants):
    """
    Write README.md with comprehensive Markdown documentation.

    Parameters
    ----------
    output_path : str
    seed : int
    bam_count : int
    variants : list
        Adapted variant definitions.
    """
    from datetime import datetime

    all_samples = get_all_samples()
    groups = list(SAMPLE_GROUPS.keys())
    type_labels = {"SNP": "SNP", "del": "del", "ins": "ins", "multi-allelic": "multi-allelic"}
    purposes_group = {"GroupA": "Population with distinct AF profile", "GroupB": "Population with distinct AF profile", "GroupC": "Smaller population group"}
    chrom_purposes = {"chr1": "Main test chromosome (carries most variants)", "chr2": "Multi-chromosome region test", "chrX": "Sex chromosome test", "chrY": "Reserved (no variants, for empty-chr testing)"}

    with open(output_path, "w") as fh:
        w = fh.write

        # --- Section 1: Overview ---
        w("# Synthetic Test Data for BaseVar Caller\n\n")
        w("## 1. Overview\n\n")
        w("This directory contains a complete synthetic test data framework for validating basevar caller. It includes:\n\n")
        w("- A mini reference genome (4 chromosomes, ~5 kb total)\n")
        w(f"- {bam_count} simulated BAM/CRAM files organized into {len(groups)} population groups\n")
        w(f"- {len(variants)} ground-truth variant sites (SNP, indel, multi-allelic)\n")
        w("- Per-sample genotype truth table\n")
        w("- Automated evaluation pipeline (`evaluate.py`)\n")
        w("- Systematic test suite covering all caller parameters (`run_tests.sh`)\n\n")
        w(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        w(f"Random seed: {seed} (fully reproducible)\n\n")

        # --- Section 2: Directory Structure ---
        w("## 2. Directory Structure\n\n")
        w("```\n")
        w("synthetic/\n")
        w("├── config.py              # Simulation configuration (variants, samples, coverage)\n")
        w("├── reference.py           # Reference genome generator\n")
        w("├── reads_simulator.py     # Read simulator (variant-aware)\n")
        w("├── writer.py              # BAM/CRAM file writer\n")
        w("├── auxiliary.py           # Auxiliary file generator (truth table, sample lists)\n")
        w("├── generate.py            # Master script: orchestrates data generation\n")
        w("├── evaluate.py            # Evaluation: VCF vs ground truth comparison\n")
        w("├── run_tests.sh           # Systematic test suite (25 scenarios)\n")
        w("│\n")
        w("├── ref/\n")
        w("│   ├── mini_ref.fa        # Reference FASTA\n")
        w("│   └── mini_ref.fa.fai    # FASTA index\n")
        w("│\n")
        w("├── bam/\n")
        w(f"│   ├── sampleA01.bam      # BAM files ({bam_count} total)\n")
        w("│   └── sampleA01.bam.bai  # BAM indices\n")
        w("│\n")
        w("├── cram/\n")
        w(f"│   ├── sampleA01.cram     # CRAM files ({bam_count} total)\n")
        w("│   └── sampleA01.cram.crai # CRAM indices\n")
        w("│\n")
        w("├── samples.list           # Absolute paths to all BAM files (for basevar -L)\n")
        w("├── samples_group.info     # Sample-to-group mapping (for basevar -G)\n")
        w("├── ground_truth_variants.tsv  # Complete truth table\n")
        w("└── README.md              # This file\n")
        w("```\n\n")

        # --- Section 3: Reference Genome ---
        w("## 3. Reference Genome Design\n\n")
        w("The reference is a small synthetic genome designed for fast testing:\n\n")
        w("| Chromosome | Length | Purpose |\n")
        w("|:----------:|-------:|:--------|\n")
        for chrom, length in CHROMOSOMES.items():
            w(f"| {chrom} | {length} bp | {chrom_purposes.get(chrom, '')} |\n")
        total_len = sum(CHROMOSOMES.values())
        w(f"\n**Total:** {total_len} bp, GC content ~45%\n\n")

        # --- Section 4: Ground Truth Variants ---
        w("## 4. Ground Truth Variants\n\n")
        w(f"{len(variants)} variant sites spanning multiple variant types:\n\n")
        w("| ID | Chrom | Pos | Ref | Alt | Type | Description |\n")
        w("|:--:|:-----:|----:|:---:|:---:|:----:|:------------|\n")
        for v in variants:
            group_af = []
            for g in groups:
                af = v["af"][g]
                group_af.append((g, sum(af) if isinstance(af, tuple) else af))
            active_groups = [g for g, a in group_af if a > 0]
            if len(active_groups) == len(groups):
                desc = "Common in all groups"
            elif len(active_groups) == 1:
                desc = f"{active_groups[0]}-specific"
            else:
                desc = ",".join(active_groups) + " enriched"
            w(f"| {v['id']} | {v['chrom']} | {v['pos']} | {v['ref']} | {v['alt']} | {type_labels.get(v['type'], v['type'])} | {desc} |\n")

        w("\n### Per-group Allele Frequencies\n\n")
        w("| ID | Type | " + " | ".join(f"{g}_AF" for g in groups) + " |\n")
        w("|:--:|:----:|" + "|".join("--------:" for _ in groups) + "|\n")
        for v in variants:
            af_parts = []
            for g in groups:
                af = v["af"][g]
                af_parts.append(",".join(f"{a:.2f}" for a in af) if isinstance(af, tuple) else f"{af:.2f}")
            w(f"| {v['id']} | {type_labels.get(v['type'], v['type'])} | " + " | ".join(af_parts) + " |\n")
        w("\n")

        # --- Section 5: Sample Design ---
        w("## 5. Sample Design\n\n")
        w(f"{bam_count} samples organized into {len(groups)} population groups:\n\n")
        w("| Group | Samples | Count | Purpose |\n")
        w("|:-----:|:--------|------:|:--------|\n")
        for group, samples in SAMPLE_GROUPS.items():
            w(f"| {group} | {samples[0]} ~ {samples[-1]} | {len(samples)} | {purposes_group.get(group, '')} |\n")

        w("\n### Coverage Distribution\n\n")
        w("| Coverage | Samples | Count | Purpose |\n")
        w("|:--------:|:--------|------:|:--------|\n")
        high_cov_samples = ", ".join(HIGH_COVERAGE_SAMPLES) if HIGH_COVERAGE_SAMPLES else "None"
        low_cov_samples = ", ".join(LOW_COVERAGE_SAMPLES)
        normal_count = bam_count - len(LOW_COVERAGE_SAMPLES) - len(HIGH_COVERAGE_SAMPLES)
        if HIGH_COVERAGE_SAMPLES:
            w(f"| {HIGH_COVERAGE}x | {high_cov_samples} | {len(HIGH_COVERAGE_SAMPLES)} | High-depth edge case |\n")
        w(f"| {NORMAL_COVERAGE}x | All others except below | {normal_count} | Standard calling depth |\n")
        w(f"| {LOW_COVERAGE}x | {low_cov_samples} | {len(LOW_COVERAGE_SAMPLES)} | Low-depth edge case |\n\n")

        # --- Section 6: Auxiliary Files ---
        w("## 6. Auxiliary Files\n\n")
        w("**`samples.list`**\n")
        w(f"One BAM path per line ({bam_count} lines). Used with `basevar -L`.\n")
        w("Contains absolute paths; regenerate with `generate.py` if directory changes.\n\n")
        w("**`samples_group.info`**\n")
        w("Tab-separated: `<sample_id>\\t<group_name>`. Used with `basevar -G`.\n\n")
        w("```\n")
        first_group = list(SAMPLE_GROUPS.keys())[0]
        first_sample = SAMPLE_GROUPS[first_group][0]
        w(f"{first_sample}    {first_group}\n")
        w("...\n")
        w("```\n\n")
        w("**`ground_truth_variants.tsv`**\n")
        w("Tab-separated truth table with header row.\n\n")
        w("| Column group | Contents |\n")
        w("|:-------------|:---------|\n")
        w("| `variant_id, chrom, pos, ref, alt, type` | Variant definition |\n")
        w("| `AF_<group>...` | Per-group allele frequency |\n")
        w("| `<sample>_GT...` | Per-sample genotype (0/0, 0/1, 1/1, 0/2, 1/2, etc.) |\n\n")
        w("Genotypes assigned via Hardy-Weinberg equilibrium from group AF.\n\n")

        # --- Section 7: Usage ---
        w("## 7. Usage: Testing BaseVar Caller\n\n")
        w("All commands assume you are in the `synthetic/` directory.\n\n")
        usage_examples = [
            ("7.1 Single-sample calling", "basevar caller -f ref/mini_ref.fa -o out.vcf bam/sampleA01.bam"),
            ("7.2 Multi-sample calling (all samples)", "basevar caller -f ref/mini_ref.fa -o out.vcf -L samples.list"),
            ("7.3 With population stratification", "basevar caller -f ref/mini_ref.fa -o out.vcf -L samples.list -G samples_group.info"),
            ("7.4 Region-specific calling", "basevar caller -f ref/mini_ref.fa -o out.vcf -L samples.list -r chr1:100-800"),
            ("7.5 Posterior genotype mode", "basevar caller -f ref/mini_ref.fa -o out.vcf -L samples.list --gt-mode"),
            ("7.6 CRAM input", "basevar caller -f ref/mini_ref.fa -o out.vcf cram/sampleA01.cram"),
            ("7.7 Smart rerun (resume from cache)", "basevar caller -f ref/mini_ref.fa -o out.vcf -L samples.list --smart-rerun"),
            ("7.8 Custom quality filters", "basevar caller -f ref/mini_ref.fa -o out.vcf -L samples.list -Q 20 -q 10 -m 0.01"),
        ]
        for title, cmd in usage_examples:
            w(f"### {title}\n\n")
            w("```bash\n")
            w(f"{cmd}\n")
            w("```\n\n")

        # --- Section 8: Automated Test Suite ---
        w("## 8. Usage: Automated Test Suite\n\n")
        w("`run_tests.sh` runs 25 test scenarios covering all caller parameters:\n\n")
        w("```bash\n")
        w("bash run_tests.sh                  # Run all tests (PASS/FAIL only)\n")
        w("bash run_tests.sh --evaluate       # Run tests + evaluate VCF vs truth\n")
        w("```\n\n")
        w("| Tests | Category | Description |\n")
        w("|:------|:---------|:------------|\n")
        w("| T01–T05 | Input modes | Single BAM, multi BAM, `-L` list, CRAM, `-L`+CRAM |\n")
        w("| T06–T08 | Region calling | `-r` chr1, `-r` chr1:100-800, `-r` chrX |\n")
        w("| T09–T10 | Population stratification | `-G`, `-G` + `-r` |\n")
        w("| T11–T12 | GT mode | Posterior vs legacy |\n")
        w("| T13–T14 | Smart rerun | Fresh + resume |\n")
        w("| T15–T16 | Quality filters | `-Q`, `-q` |\n")
        w("| T17–T20 | Other parameters | `--filename-has-samplename`, `--ref-bias`, `--max-alleles`, `-B` |\n")
        w("| T21–T25 | Combo workflows | Full combo, concat, multi-region, etc. |\n\n")

        # --- Section 9: Evaluation ---
        w("## 9. Usage: Evaluation Pipeline\n\n")
        w("`evaluate.py` compares a VCF against `ground_truth_variants.tsv`:\n\n")
        w("```bash\n")
        w("python3 evaluate.py --vcf out.vcf --truth ground_truth_variants.tsv\n")
        w("```\n\n")
        w("### Output Metrics\n\n")
        w("- **Sensitivity** (recall): fraction of truth variants found in VCF\n")
        w("- **Precision**: fraction of VCF variants matching truth\n")
        w("- **GT concordance**: per-sample genotype agreement rate\n")
        w("- **AF accuracy**: MAE / RMSE / Pearson *r* between called AF and truth AF\n")
        w("- **Per-group AF breakdown**\n\n")
        w("> Indel matching tolerates VCF left-alignment position offsets (search radius = max allele length + 1).\n\n")

        # --- Section 10: Regenerating Data ---
        w("## 10. Regenerating Data\n\n")
        w("If you modify `config.py` (e.g., change sample count, variants, coverage):\n\n")
        w("```bash\n")
        w("python3 generate.py                  # Generate BAM + CRAM\n")
        w("python3 generate.py --skip-cram      # BAM only (faster)\n")
        w("```\n\n")
        w("This will regenerate all files: `ref/`, `bam/`, `cram/`, `samples.list`, `samples_group.info`, `ground_truth_variants.tsv`.\n\n")
        w("After regeneration, run the test suite:\n\n")
        w("```bash\n")
        w("bash run_tests.sh --evaluate\n")
        w("```\n")

    print(f"  [auxiliary] Written README to {output_path}")
