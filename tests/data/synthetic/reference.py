"""
Reference genome generator.

Generates a small synthetic reference FASTA with controlled GC content
and builds samtools faidx index.
"""

import os
import numpy as np
import pysam

from config import CHROMOSOMES, GC_CONTENT, DEFAULT_SEED


def _generate_sequence(length, gc_content, rng):
    """
    Generate a random DNA sequence of given length with target GC content.

    Returns a string of A/C/G/T characters.
    """
    # Probability of each base: A, C, G, T
    # GC content = P(C) + P(G), so P(C) = P(G) = gc/2
    # AT content = 1 - gc, so P(A) = P(T) = (1-gc)/2
    p_gc = gc_content / 2.0
    p_at = (1.0 - gc_content) / 2.0
    probs = [p_at, p_gc, p_gc, p_at]  # A, C, G, T
    bases = np.array(list("ACGT"))
    indices = rng.choice(4, size=length, p=probs)
    return "".join(bases[indices])


def generate_reference(output_dir, seed=DEFAULT_SEED):
    """
    Generate reference FASTA and build index.

    Parameters
    ----------
    output_dir : str
        Directory to write ref/mini_ref.fa (and .fai).
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    fa_path : str
        Path to the generated FASTA file.
    ref_seqs : dict
        {chrom_name: sequence_string} for downstream use.
    """
    rng = np.random.default_rng(seed)
    os.makedirs(output_dir, exist_ok=True)
    fa_path = os.path.join(output_dir, "mini_ref.fa")

    ref_seqs = {}
    with open(fa_path, "w") as fh:
        for chrom, length in CHROMOSOMES.items():
            seq = _generate_sequence(length, GC_CONTENT, rng)
            ref_seqs[chrom] = seq
            # Write FASTA with 60-char line wrapping
            fh.write(f">{chrom}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")

    print(f"  [reference] Written {len(ref_seqs)} chromosomes to {fa_path}")
    for chrom, seq in ref_seqs.items():
        gc = (seq.count("G") + seq.count("C")) / len(seq)
        print(f"    {chrom}: {len(seq)} bp, GC={gc:.2%}")

    # Build faidx index using pysam
    pysam.faidx(fa_path)
    print(f"  [reference] Built index: {fa_path}.fai")

    return fa_path, ref_seqs


def adapt_variants_to_ref(ref_seqs, variants):
    """
    Adapt variant definitions to match the actual reference sequence.

    For SNPs: sets REF to the actual reference base, picks a different base as ALT.
    For Indels: constructs proper REF/ALT from the reference context.

    This modifies the variants list in-place and returns it.

    Parameters
    ----------
    ref_seqs : dict
        {chrom: sequence_string} from generate_reference().
    variants : list
        From config.VARIANTS (will be modified in-place).

    Returns
    -------
    variants : list
        Updated variant definitions with correct REF/ALT matching the reference.
    """
    ALT_CHOICES = ["A", "C", "G", "T"]

    for v in variants:
        chrom = v["chrom"]
        pos_0 = v["pos"] - 1  # 1-based -> 0-based
        vtype = v["type"]

        if vtype == "SNP":
            actual_ref = ref_seqs[chrom][pos_0].upper()
            v["ref"] = actual_ref
            # Pick first alt that differs from ref
            for alt in ALT_CHOICES:
                if alt != actual_ref:
                    v["alt"] = alt
                    break

        elif vtype == "multi-allelic":
            actual_ref = ref_seqs[chrom][pos_0].upper()
            v["ref"] = actual_ref
            # Pick two alts that differ from ref and each other
            alts = []
            for alt in ALT_CHOICES:
                if alt != actual_ref:
                    alts.append(alt)
                if len(alts) == 2:
                    break
            v["alt"] = ",".join(alts)

        elif vtype == "del":
            # Deletion: REF = anchor + deleted bases, ALT = anchor
            # The anchor is the base BEFORE the deletion (VCF convention)
            # For simplicity, REF starts at pos, ALT is the first base
            actual_ref_from_pos = ref_seqs[chrom][pos_0:pos_0 + len(v["ref"])].upper()
            # Use actual reference bases
            v["ref"] = actual_ref_from_pos
            # ALT = first base of ref (anchor) only
            v["alt"] = actual_ref_from_pos[0]

        elif vtype == "ins":
            # Insertion: REF = anchor base, ALT = anchor + inserted bases
            anchor = ref_seqs[chrom][pos_0].upper()
            ins_len = len(v["alt"]) - len(v["ref"])
            # Generate insertion bases deterministically (not from reference)
            # Use a simple pattern based on position
            ins_bases = ""
            bases_list = ["A", "C", "G", "T"]
            for i in range(ins_len):
                ins_bases += bases_list[(pos_0 + i) % 4]
            v["ref"] = anchor
            v["alt"] = anchor + ins_bases

    return variants


def apply_variants_to_ref(ref_seqs, variants):
    """
    Build a variant map from adapted variant definitions, validated against
    the reference sequence.

    Parameters
    ----------
    ref_seqs : dict
        {chrom: sequence_string} from generate_reference().
    variants : list
        Adapted variant definitions (after adapt_variants_to_ref).

    Returns
    -------
    variant_map : dict
        {chrom: [(pos_0based, ref, alt, type, var_id), ...]}
        Sorted by position.
    """
    variant_map = {}
    for v in variants:
        chrom = v["chrom"]
        pos_0 = v["pos"] - 1
        ref_allele = v["ref"]
        alt_allele = v["alt"]

        # Validate ref allele matches reference
        ref_len = len(ref_allele)
        actual_ref = ref_seqs[chrom][pos_0:pos_0 + ref_len]
        if actual_ref.upper() != ref_allele.upper():
            raise ValueError(
                f"Variant {v['id']}: REF '{ref_allele}' does not match "
                f"reference at {chrom}:{v['pos']} (found '{actual_ref}')"
            )

        if chrom not in variant_map:
            variant_map[chrom] = []
        variant_map[chrom].append((pos_0, ref_allele, alt_allele, v["type"], v["id"]))

    # Sort by position
    for chrom in variant_map:
        variant_map[chrom].sort(key=lambda x: x[0])

    return variant_map
