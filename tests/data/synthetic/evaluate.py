#!/usr/bin/env python3
"""
Evaluation script for BaseVar caller output against ground truth.

Compares the called VCF against ground_truth_variants.tsv to compute:
  - Site-level sensitivity (recall): fraction of truth variants detected
  - Site-level precision: fraction of called variants that are true
  - Genotype concordance: per-sample GT agreement at truth sites
  - AF accuracy: correlation/bias of estimated AF vs truth AF
  - Per-group AF accuracy: AF_GroupA/B/C vs truth

Usage:
    python3 evaluate.py --vcf out.vcf --truth ground_truth_variants.tsv
    python3 evaluate.py --vcf out.vcf --truth ground_truth_variants.tsv --output eval_report.tsv

Author: BaseVar2 test data pipeline
"""

import argparse
import os
import sys
import gzip
from collections import defaultdict

# Ensure the script's directory is in the Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


# ---------------------------------------------------------------------------
# VCF parser (lightweight, no pysam dependency)
# ---------------------------------------------------------------------------

def _open_vcf(vcf_path):
    """
    Open a VCF file, handling plain text, gzip, and BGZF compression.

    Returns an iterable of text lines.
    """
    # Check file magic bytes
    with open(vcf_path, "rb") as f:
        magic = f.read(2)

    if magic == b'\x1f\x8b':
        # gzip/BGZF compressed
        return gzip.open(vcf_path, "rt")
    else:
        return open(vcf_path, "r")


def parse_vcf(vcf_path):
    """
    Parse a basevar VCF file and extract variant records.
    Supports plain text, gzip, and BGZF compressed VCF.

    Returns
    -------
    called_variants : dict
        {(chrom, pos): {
            "ref": str,
            "alt": [str],          # list of alt alleles
            "qual": float,
            "info": dict,          # parsed INFO fields
            "format_fields": [str],
            "samples": {sample_id: {field: value}},
        }}
    sample_order : list of str
    """
    called = {}
    sample_order = []

    with _open_vcf(vcf_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                # Header line: extract sample names
                parts = line.split("\t")
                sample_order = parts[9:]
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                continue

            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3]
            alt_list = parts[4].split(",")
            qual_str = parts[5]
            qual = float(qual_str) if qual_str != "." else 0.0
            info_str = parts[7]
            fmt_fields = parts[8].split(":")

            # Parse INFO
            info = {}
            for field in info_str.split(";"):
                if "=" in field:
                    key, val = field.split("=", 1)
                    info[key] = val

            # Parse samples
            samples = {}
            for i, sample_data in enumerate(parts[9:]):
                if i >= len(sample_order):
                    break
                sample_id = sample_order[i]
                values = sample_data.split(":")
                sample_dict = {}
                for j, fmt_key in enumerate(fmt_fields):
                    if j < len(values):
                        sample_dict[fmt_key] = values[j]
                    else:
                        sample_dict[fmt_key] = "."
                samples[sample_id] = sample_dict

            key = (chrom, pos)
            called[key] = {
                "ref": ref,
                "alt": alt_list,
                "qual": qual,
                "info": info,
                "format_fields": fmt_fields,
                "samples": samples,
            }

    return called, sample_order


# ---------------------------------------------------------------------------
# Ground truth parser
# ---------------------------------------------------------------------------

def parse_ground_truth(truth_path):
    """
    Parse ground_truth_variants.tsv.

    Returns
    -------
    truth_variants : list of dict
        Each dict has keys: id, chrom, pos, ref, alt, type,
        af_per_group: {group: af_value},
        gt_per_sample: {sample_id: genotype_string}
    groups : list of str
    samples : list of str
    """
    truth = []
    groups = []
    samples = []

    with open(truth_path, "r") as fh:
        header_line = fh.readline().strip()
        header_parts = header_line.split("\t")

        # Parse header to find groups and samples
        for part in header_parts:
            if part.startswith("AF_"):
                groups.append(part[3:])  # Remove "AF_" prefix
            elif part.endswith("_GT"):
                samples.append(part[:-3])  # Remove "_GT" suffix

        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")

            record = {
                "id": parts[0],
                "chrom": parts[1],
                "pos": int(parts[2]),
                "ref": parts[3],
                "alt": parts[4],
                "type": parts[5],
            }

            # Parse AF per group
            af_per_group = {}
            idx = 6
            for group in groups:
                af_str = parts[idx]
                if "," in af_str:
                    af_per_group[group] = tuple(float(x) for x in af_str.split(","))
                else:
                    af_per_group[group] = float(af_str)
                idx += 1
            record["af_per_group"] = af_per_group

            # Parse GT per sample
            gt_per_sample = {}
            for sample_id in samples:
                gt_per_sample[sample_id] = parts[idx]
                idx += 1
            record["gt_per_sample"] = gt_per_sample

            truth.append(record)

    return truth, groups, samples


# ---------------------------------------------------------------------------
# Matching logic
# ---------------------------------------------------------------------------

def match_truth_to_called(truth_variants, called_variants):
    """
    Match each truth variant to a called variant.

    Handles VCF indel convention where the called position may be shifted
    by 1bp relative to the truth (VCF uses anchor base + left-alignment).

    Returns
    -------
    matched : list of (truth_record, called_record or None)
    unmatched_called : list of called_record
    """
    matched = []
    matched_called_keys = set()

    for tv in truth_variants:
        chrom = tv["chrom"]
        pos = tv["pos"]
        ref = tv["ref"].upper()
        alt = tv["alt"].upper()
        truth_alts = [a.upper() for a in alt.split(",")]
        vtype = tv["type"]

        # Determine search radius based on variant type
        # Indels may shift by up to len(ref) due to VCF left-alignment
        if vtype in ("del", "ins"):
            max_offset = max(len(ref), max(len(a) for a in truth_alts)) + 1
        else:
            max_offset = 0

        # Search for a match at the exact position and nearby positions
        found = False
        for offset in range(0, max_offset + 1):
            for sign in ([0] if offset == 0 else [-1, 1]):
                test_pos = pos + sign * offset
                key = (chrom, test_pos)
                called_rec = called_variants.get(key)
                if called_rec is None:
                    continue

                called_alts = [a.upper() for a in called_rec["alt"]]

                # Check allele overlap: at least one truth alt must appear
                # in called alts, OR the called alt shares bases with truth
                # (for indels, the alt allele string may differ due to
                # different anchor bases but still represent the same event)
                if any(ta in called_alts for ta in truth_alts):
                    matched.append((tv, called_rec))
                    matched_called_keys.add(key)
                    found = True
                    break

                # For indels: check if the variant event is the same
                # by comparing the net change (inserted/deleted bases)
                if vtype in ("del", "ins"):
                    for ta in truth_alts:
                        for ca in called_alts:
                            # Same net event: ref->alt has same length change
                            truth_delta = len(ta) - len(ref)
                            called_delta = len(ca) - len(called_rec["ref"].upper())
                            if truth_delta == called_delta and truth_delta != 0:
                                matched.append((tv, called_rec))
                                matched_called_keys.add(key)
                                found = True
                                break
                        if found:
                            break
                if found:
                    break
            if found:
                break

        if not found:
            matched.append((tv, None))

    # Find unmatched called variants (potential false positives)
    unmatched_called = []
    for key, rec in called_variants.items():
        if key not in matched_called_keys:
            unmatched_called.append(rec)

    return matched, unmatched_called


# ---------------------------------------------------------------------------
# Genotype comparison
# ---------------------------------------------------------------------------

def parse_gt_string(gt_str):
    """
    Parse a VCF GT string like '0/1', '1/1', './.' into a tuple of alleles.
    Returns None for missing genotypes.
    """
    if gt_str in ("./.", ".", ""):
        return None
    sep = "/" if "/" in gt_str else "|"
    return tuple(int(x) for x in gt_str.split(sep))


def gt_concordance(gt_truth, gt_called):
    """
    Compare two GT strings. Returns True if they match, False otherwise.
    Missing called GT returns None.
    """
    t = parse_gt_string(gt_truth)
    c = parse_gt_string(gt_called)
    if t is None or c is None:
        return None
    # Compare as sorted tuples (phasing doesn't matter)
    return sorted(t) == sorted(c)


# ---------------------------------------------------------------------------
# Metrics computation
# ---------------------------------------------------------------------------

def compute_metrics(matched, unmatched_called, truth_variants, groups, samples):
    """
    Compute all evaluation metrics.

    Returns
    -------
    metrics : dict of metric_name -> value
    per_variant : list of dicts for detailed per-variant report
    """
    metrics = {}
    per_variant = []

    n_truth = len(truth_variants)
    n_matched = sum(1 for (_, c) in matched if c is not None)
    n_missed = n_truth - n_matched
    n_fp = len(unmatched_called)

    # --- Site-level metrics ---
    metrics["n_truth_variants"] = n_truth
    metrics["n_called_variants"] = n_matched + n_fp
    metrics["n_true_positives"] = n_matched
    metrics["n_false_negatives"] = n_missed
    metrics["n_false_positives"] = n_fp
    metrics["sensitivity"] = n_matched / n_truth if n_truth > 0 else float("nan")
    metrics["precision"] = n_matched / (n_matched + n_fp) if (n_matched + n_fp) > 0 else float("nan")
    if metrics["sensitivity"] + metrics["precision"] > 0:
        metrics["f1_score"] = 2 * metrics["sensitivity"] * metrics["precision"] / (metrics["sensitivity"] + metrics["precision"])
    else:
        metrics["f1_score"] = float("nan")

    # --- Genotype concordance (per-sample, across all matched sites) ---
    total_gt_comparisons = 0
    concordant = 0
    discordant = 0
    missing = 0

    # Per-sample GT concordance
    per_sample_concordance = defaultdict(lambda: {"concordant": 0, "discordant": 0, "missing": 0})

    # AF error accumulation
    af_errors = []  # (truth_af, called_af) pairs
    per_group_af_errors = {g: [] for g in groups}

    for tv, cv in matched:
        pv_report = {
            "variant_id": tv["id"],
            "chrom": tv["chrom"],
            "pos": tv["pos"],
            "truth_ref": tv["ref"],
            "truth_alt": tv["alt"],
            "truth_type": tv["type"],
            "detected": cv is not None,
        }

        if cv is None:
            pv_report["called_gt_concordance"] = "MISSED"
            per_variant.append(pv_report)
            continue

        # Per-sample GT comparison
        for sample_id in samples:
            truth_gt = tv["gt_per_sample"].get(sample_id, "0/0")
            called_sample = cv["samples"].get(sample_id, {})
            called_gt = called_sample.get("GT", "./.")

            result = gt_concordance(truth_gt, called_gt)
            if result is None:
                missing += 1
                per_sample_concordance[sample_id]["missing"] += 1
            elif result:
                concordant += 1
                total_gt_comparisons += 1
                per_sample_concordance[sample_id]["concordant"] += 1
            else:
                discordant += 1
                total_gt_comparisons += 1
                per_sample_concordance[sample_id]["discordant"] += 1

        # AF comparison (overall AF)
        truth_af_overall = _compute_overall_af(tv, groups)
        called_af_str = cv["info"].get("AF", None)
        if called_af_str is not None and truth_af_overall is not None:
            try:
                called_af = float(called_af_str)
                af_errors.append((truth_af_overall, called_af))
            except ValueError:
                pass

        # Per-group AF comparison
        for group in groups:
            truth_gaf = tv["af_per_group"].get(group)
            called_gaf_str = cv["info"].get(f"AF_{group}", None)
            if truth_gaf is not None and called_gaf_str is not None:
                try:
                    # For multi-allelic, compare first AF only
                    if isinstance(truth_gaf, tuple):
                        truth_gaf_val = truth_gaf[0]
                    else:
                        truth_gaf_val = truth_gaf
                    called_gaf = float(called_gaf_str)
                    per_group_af_errors[group].append((truth_gaf_val, called_gaf))
                except ValueError:
                    pass

        per_variant.append(pv_report)

    # Genotype concordance rate
    total_with_gt = concordant + discordant
    metrics["gt_concordance_rate"] = concordant / total_with_gt if total_with_gt > 0 else float("nan")
    metrics["gt_comparisons_total"] = total_with_gt
    metrics["gt_missing_count"] = missing

    # AF metrics
    if af_errors:
        errors = [abs(t - c) for t, c in af_errors]
        metrics["af_mae"] = sum(errors) / len(errors)  # Mean absolute error
        metrics["af_rmse"] = (sum(e**2 for e in errors) / len(errors)) ** 0.5
        metrics["af_max_error"] = max(errors)
        # Correlation
        metrics["af_correlation"] = _pearson_r([x[0] for x in af_errors], [x[1] for x in af_errors])

    # Per-group AF metrics
    for group in groups:
        gaf = per_group_af_errors[group]
        if gaf:
            errors = [abs(t - c) for t, c in gaf]
            metrics[f"af_mae_{group}"] = sum(errors) / len(errors)
            metrics[f"af_rmse_{group}"] = (sum(e**2 for e in errors) / len(errors)) ** 0.5

    # Per-sample concordance
    per_sample_metrics = {}
    for sample_id in samples:
        sc = per_sample_concordance[sample_id]
        total_s = sc["concordant"] + sc["discordant"]
        per_sample_metrics[sample_id] = {
            "concordant": sc["concordant"],
            "discordant": sc["discordant"],
            "missing": sc["missing"],
            "concordance_rate": sc["concordant"] / total_s if total_s > 0 else float("nan"),
        }
    metrics["per_sample"] = per_sample_metrics

    return metrics, per_variant


def _compute_overall_af(truth_record, groups):
    """Compute the overall AF across all groups (weighted by group size)."""
    from config import SAMPLE_GROUPS
    total_alleles = 0
    total_alt = 0.0
    for group in groups:
        n = len(SAMPLE_GROUPS.get(group, []))
        af = truth_record["af_per_group"].get(group)
        if af is not None:
            if isinstance(af, tuple):
                af = sum(af)  # Sum of alt allele frequencies for multi-allelic
            total_alleles += 2 * n
            total_alt += af * 2 * n
    return total_alt / total_alleles if total_alleles > 0 else None


def _pearson_r(x, y):
    """Compute Pearson correlation coefficient."""
    n = len(x)
    if n < 2:
        return float("nan")
    mx = sum(x) / n
    my = sum(y) / n
    sx = sum((xi - mx) ** 2 for xi in x) ** 0.5
    sy = sum((yi - my) ** 2 for yi in y) ** 0.5
    if sx == 0 or sy == 0:
        return float("nan")
    cov = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    return cov / (sx * sy)


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def print_report(metrics, per_variant, groups, samples):
    """Print a human-readable evaluation report."""
    print()
    print("=" * 70)
    print("  BaseVar Caller Evaluation Report")
    print("=" * 70)

    print()
    print("--- Site-Level Metrics ---")
    print(f"  Truth variants:        {metrics['n_truth_variants']}")
    print(f"  Called variants:       {metrics['n_called_variants']}")
    print(f"  True positives:        {metrics['n_true_positives']}")
    print(f"  False negatives:       {metrics['n_false_negatives']}")
    print(f"  False positives:       {metrics['n_false_positives']}")
    print(f"  Sensitivity (recall):  {metrics['sensitivity']:.4f}")
    print(f"  Precision:             {metrics['precision']:.4f}")
    print(f"  F1 score:              {metrics['f1_score']:.4f}")

    print()
    print("--- Genotype Concordance ---")
    print(f"  Total comparisons:     {metrics['gt_comparisons_total']}")
    print(f"  Concordant:            {metrics['gt_concordance_rate']:.4f}" if metrics['gt_concordance_rate'] == metrics['gt_concordance_rate'] else "  Concordant:            N/A")
    print(f"  Missing GTs:           {metrics['gt_missing_count']}")

    print()
    print("--- AF Accuracy ---")
    if "af_mae" in metrics:
        print(f"  Overall AF MAE:        {metrics['af_mae']:.4f}")
        print(f"  Overall AF RMSE:       {metrics['af_rmse']:.4f}")
        print(f"  Overall AF max error:  {metrics['af_max_error']:.4f}")
        print(f"  AF correlation (r):    {metrics['af_correlation']:.4f}")

    for group in groups:
        mae_key = f"af_mae_{group}"
        rmse_key = f"af_rmse_{group}"
        if mae_key in metrics:
            print(f"  {group} AF MAE:          {metrics[mae_key]:.4f}")
            print(f"  {group} AF RMSE:         {metrics[rmse_key]:.4f}")

    print()
    print("--- Per-Sample GT Concordance ---")
    print(f"  {'Sample':<14} {'Concord':>8} {'Discord':>8} {'Missing':>8} {'Rate':>8}")
    print(f"  {'-'*14} {'-'*8} {'-'*8} {'-'*8} {'-'*8}")
    per_sample = metrics.get("per_sample", {})
    for sample_id in samples:
        sm = per_sample.get(sample_id, {})
        rate_str = f"{sm['concordance_rate']:.4f}" if sm.get('concordance_rate', float('nan')) == sm.get('concordance_rate', float('nan')) else "N/A"
        print(f"  {sample_id:<14} {sm.get('concordant', 0):>8} {sm.get('discordant', 0):>8} {sm.get('missing', 0):>8} {rate_str:>8}")

    print()
    print("--- Per-Variant Detail ---")
    print(f"  {'ID':<6} {'Pos':<12} {'Type':<16} {'Detected':<10}")
    print(f"  {'-'*6} {'-'*12} {'-'*16} {'-'*10}")
    for pv in per_variant:
        detected_str = "YES" if pv["detected"] else "MISSED"
        print(f"  {pv['variant_id']:<6} {pv['chrom']}:{pv['pos']:<6} {pv['truth_type']:<16} {detected_str:<10}")

    print()
    print("=" * 70)


def write_tsv_report(metrics, per_variant, groups, samples, output_path):
    """Write a TSV evaluation report."""
    with open(output_path, "w") as fh:
        # Summary section
        fh.write("# BaseVar Caller Evaluation Report\n")
        fh.write(f"# Truth variants: {metrics['n_truth_variants']}\n")
        fh.write(f"# Called variants: {metrics['n_called_variants']}\n")
        fh.write(f"# Sensitivity: {metrics['sensitivity']:.4f}\n")
        fh.write(f"# Precision: {metrics['precision']:.4f}\n")
        fh.write(f"# F1: {metrics['f1_score']:.4f}\n")
        fh.write(f"# GT concordance: {metrics.get('gt_concordance_rate', float('nan')):.4f}\n")
        if "af_mae" in metrics:
            fh.write(f"# AF MAE: {metrics['af_mae']:.4f}\n")
            fh.write(f"# AF RMSE: {metrics['af_rmse']:.4f}\n")
        fh.write("#\n")

        # Per-variant detail
        fh.write("variant_id\tchrom\tpos\ttype\tdetected\n")
        for pv in per_variant:
            fh.write(f"{pv['variant_id']}\t{pv['chrom']}\t{pv['pos']}\t{pv['truth_type']}\t{'YES' if pv['detected'] else 'MISSED'}\n")

    print(f"  [evaluate] Report written to {output_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Evaluate BaseVar caller output against ground truth."
    )
    parser.add_argument(
        "--vcf", required=True,
        help="Path to the basevar caller VCF output."
    )
    parser.add_argument(
        "--truth", required=True,
        help="Path to ground_truth_variants.tsv."
    )
    parser.add_argument(
        "--output", default=None,
        help="Output TSV report path (optional)."
    )
    args = parser.parse_args()

    # Parse inputs
    print(f"  [evaluate] Loading VCF: {args.vcf}")
    called_variants, vcf_samples = parse_vcf(args.vcf)
    print(f"  [evaluate] Loaded {len(called_variants)} called variants, {len(vcf_samples)} samples")

    print(f"  [evaluate] Loading truth: {args.truth}")
    truth_variants, groups, truth_samples = parse_ground_truth(args.truth)
    print(f"  [evaluate] Loaded {len(truth_variants)} truth variants, {len(groups)} groups, {len(truth_samples)} samples")

    # Match truth to called
    matched, unmatched_called = match_truth_to_called(truth_variants, called_variants)

    # Compute metrics
    metrics, per_variant = compute_metrics(matched, unmatched_called, truth_variants, groups, truth_samples)

    # Print report
    print_report(metrics, per_variant, groups, truth_samples)

    # Write TSV report if requested
    if args.output:
        write_tsv_report(metrics, per_variant, groups, truth_samples, args.output)


if __name__ == "__main__":
    main()
