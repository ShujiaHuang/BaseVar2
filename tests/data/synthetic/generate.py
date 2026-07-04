#!/usr/bin/env python3
"""
Synthetic test data generator for BaseVar caller.

Generates a complete set of test data including:
  - Reference genome (FASTA + index)
  - 25 BAM files (sorted, indexed) across 3 population groups
  - 25 CRAM files (same data, CRAM format)
  - Ground truth variants table with per-sample genotypes
  - Auxiliary files (samples.list, samples_group.info)

Usage:
    python3 generate.py [--seed SEED] [--output-dir DIR]

Author: BaseVar2 test data pipeline
"""

import argparse
import os
import sys
import time

# Ensure the script's directory is in the Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from config import (
    DEFAULT_SEED, VARIANTS, CHROMOSOMES, SAMPLE_GROUPS,
    get_all_samples, PAIRED_END_ENABLED, PAIRED_READ_LENGTH,
    PAIRED_INSERT_SIZE_MEAN, PAIRED_INSERT_SIZE_SD,
)
from reference import generate_reference, adapt_variants_to_ref, apply_variants_to_ref
from reads_simulator import assign_all_genotypes, generate_reads_for_sample
from writer import write_all_samples
from auxiliary import (
    write_samples_list, write_samples_group,
    write_ground_truth, write_readme,
)


def main():
    parser = argparse.ArgumentParser(
        description="Generate synthetic test data for BaseVar caller."
    )
    parser.add_argument(
        "--seed", type=int, default=DEFAULT_SEED,
        help=f"Random seed for reproducibility (default: {DEFAULT_SEED})"
    )
    parser.add_argument(
        "--output-dir", type=str, default=os.path.dirname(os.path.abspath(__file__)),
        help="Output directory (default: script directory)"
    )
    parser.add_argument(
        "--skip-cram", action="store_true",
        help="Skip CRAM generation (faster, for quick testing)"
    )
    parser.add_argument(
        "--basevar", type=str, default=None,
        help="Path to basevar executable for automatic smoke test (e.g. ../../bin/basevar)"
    )
    parser.add_argument(
        "--paired-end", action="store_true",
        help="Generate paired-end reads instead of single-end"
    )
    parser.add_argument(
        "--paired-read-length", type=int, default=None,
        help="Read length for each paired-end mate (default: config value)"
    )
    parser.add_argument(
        "--insert-size-mean", type=float, default=None,
        help="Mean template length for paired-end fragments (bp)"
    )
    parser.add_argument(
        "--insert-size-sd", type=float, default=None,
        help="Standard deviation for template length in paired-end fragments"
    )
    parser.add_argument(
        "--evaluate", action="store_true",
        help="Run evaluation after smoke test (requires --basevar)"
    )
    parser.add_argument(
        "--ancient-dna", action="store_true",
        help="Enable ancient DNA simulation (short fragments + PMD damage)"
    )
    args = parser.parse_args()

    seed = args.seed
    base_dir = args.output_dir
    ref_dir = os.path.join(base_dir, "ref")
    bam_dir = os.path.join(base_dir, "bam")
    cram_dir = os.path.join(base_dir, "cram")

    print("=" * 60)
    print("BaseVar Synthetic Test Data Generator")
    print("=" * 60)
    print(f"  Seed:       {seed}")
    print(f"  Output dir: {base_dir}")
    print(f"  Samples:    {len(get_all_samples())} ({len(SAMPLE_GROUPS)} groups)")
    print(f"  Variants:   {len(VARIANTS)}")
    print(f"  Chromosomes: {len(CHROMOSOMES)}")
    if args.ancient_dna:
        from config import (
            ANCIENT_DNA_FRAG_LENGTH_MEAN, ANCIENT_DNA_DAMAGE_RATE,
            ANCIENT_DNA_CONTAMINATION_RATE,
        )
        print(f"  Ancient DNA: ENABLED")
        print(f"    Fragment mean: {ANCIENT_DNA_FRAG_LENGTH_MEAN} bp")
        print(f"    Damage rate:   {ANCIENT_DNA_DAMAGE_RATE}")
        print(f"    Contamination: {ANCIENT_DNA_CONTAMINATION_RATE}")
    print()

    t0 = time.time()

    # ---------------------------------------------------------------
    # Step 1: Generate reference genome
    # ---------------------------------------------------------------
    print("[Step 1/6] Generating reference genome...")
    fa_path, ref_seqs = generate_reference(ref_dir, seed=seed)
    print()

    # ---------------------------------------------------------------
    # Step 2: Validate variants against reference
    # ---------------------------------------------------------------
    print("[Step 2/6] Adapting and validating variants against reference...")
    # First adapt REF/ALT to match the actual reference sequence
    import copy
    adapted_variants = copy.deepcopy(VARIANTS)
    adapt_variants_to_ref(ref_seqs, adapted_variants)
    try:
        variant_map = apply_variants_to_ref(ref_seqs, adapted_variants)
        for chrom, variants_on_chrom in variant_map.items():
            for (pos, ref, alt, vtype, vid) in variants_on_chrom:
                print(f"    {vid}: {chrom}:{pos+1} {ref}->{alt} ({vtype}) OK")
    except ValueError as e:
        print(f"  ERROR: {e}")
        print("  Please adjust variant positions in config.py to match the reference.")
        sys.exit(1)
    print()

    # ---------------------------------------------------------------
    # Step 3: Assign genotypes
    # ---------------------------------------------------------------
    print("[Step 3/6] Assigning genotypes to samples...")
    sample_genotypes = assign_all_genotypes(adapted_variants, seed=seed)

    # Print genotype summary
    for v in adapted_variants:
        vid = v["id"]
        gt_counts = {}
        for sample_id in get_all_samples():
            gt = sample_genotypes[sample_id][vid]
            gt_counts[gt] = gt_counts.get(gt, 0) + 1
        gt_str = ", ".join(f"{gt}={n}" for gt, n in sorted(gt_counts.items()))
        print(f"    {vid}: {gt_str}")
    print()

    # ---------------------------------------------------------------
    # Step 4: Generate reads
    # ---------------------------------------------------------------
    print("[Step 4/6] Generating reads for all samples...")
    all_sample_reads = {}
    for sample_id in get_all_samples():
        reads = generate_reads_for_sample(
            sample_id, ref_seqs, variant_map, sample_genotypes, seed=seed,
            paired=args.paired_end,
            paired_read_length=args.paired_read_length,
            insert_size_mean=args.insert_size_mean,
            insert_size_sd=args.insert_size_sd,
        )
        all_sample_reads[sample_id] = reads
        print(f"    {sample_id}: {len(reads)} reads")
    print()

    # ---------------------------------------------------------------
    # Step 4b (optional): Apply ancient DNA simulation
    # ---------------------------------------------------------------
    if args.ancient_dna:
        print("[Step 4b] Applying ancient DNA simulation (fragmentation + PMD)...")
        from ancient_dna import simulate_ancient_dna
        # Pass ref_seqs so the simulator can recompute NM/MD tags exactly
        all_sample_reads = simulate_ancient_dna(all_sample_reads, ref_seqs, seed=seed)
        print()

    # ---------------------------------------------------------------
    # Step 5: Write BAM/CRAM files
    # ---------------------------------------------------------------
    print("[Step 5/6] Writing BAM/CRAM files...")
    if args.skip_cram:
        # Only write BAM
        from writer import write_bam
        bam_paths = {}
        for sample_id, reads in all_sample_reads.items():
            print(f"  [writer] Writing BAM: {sample_id} ({len(reads)} reads)")
            bam_paths[sample_id] = write_bam(reads, sample_id, bam_dir, CHROMOSOMES)
        cram_paths = {}
    else:
        bam_paths, cram_paths = write_all_samples(
            all_sample_reads, bam_dir, cram_dir, fa_path, CHROMOSOMES
        )
    print()

    # ---------------------------------------------------------------
    # Step 6: Write auxiliary files
    # ---------------------------------------------------------------
    print("[Step 6/6] Writing auxiliary files...")
    write_samples_list(bam_paths, os.path.join(base_dir, "samples.list"))
    write_samples_group(os.path.join(base_dir, "samples_group.info"))
    write_ground_truth(sample_genotypes, adapted_variants, os.path.join(base_dir, "ground_truth_variants.tsv"))
    write_readme(
        os.path.join(base_dir, "README.md"),
        seed=seed,
        bam_count=len(bam_paths),
        variants=adapted_variants
    )
    print()

    elapsed = time.time() - t0
    print("=" * 60)
    print(f"Done! All files written to: {base_dir}")
    print(f"Total time: {elapsed:.1f} seconds")
    print("=" * 60)

    # ---------------------------------------------------------------
    # Optional: Smoke test with basevar caller
    # ---------------------------------------------------------------
    if args.basevar:
        import subprocess
        basevar_exe = os.path.abspath(args.basevar)
        print()
        print("=" * 60)
        print("  Smoke Test: basevar caller")
        print("=" * 60)

        test_vcf = os.path.join(base_dir, "_smoke_test.vcf")
        test_bam = bam_paths.get("sampleA01", "")
        if test_bam and os.path.exists(test_bam):
            cmd = [
                basevar_exe, "caller",
                "-f", fa_path,
                "-o", test_vcf,
                test_bam,
            ]
            print(f"  Running: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=base_dir)
            if result.returncode == 0:
                print(f"  [PASS] Smoke test succeeded. Output: {test_vcf}")

                # Optional: run evaluation
                if args.evaluate and os.path.exists(test_vcf):
                    print()
                    print("  Running evaluation...")
                    from evaluate import parse_vcf, parse_ground_truth, match_truth_to_called, compute_metrics, print_report
                    called, vcf_samples = parse_vcf(test_vcf)
                    truth, groups, truth_samples = parse_ground_truth(
                        os.path.join(base_dir, "ground_truth_variants.tsv")
                    )
                    matched, unmatched = match_truth_to_called(truth, called)
                    metrics, per_var = compute_metrics(matched, unmatched, truth, groups, truth_samples)
                    print_report(metrics, per_var, groups, truth_samples)
            else:
                print(f"  [FAIL] Smoke test failed (exit code {result.returncode})")
                print(f"  stderr: {result.stderr[:500]}")
        else:
            print("  [SKIP] No BAM file found for smoke test")


if __name__ == "__main__":
    main()
