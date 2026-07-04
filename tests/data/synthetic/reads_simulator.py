"""
Read simulator for synthetic basevar test data.

Generates reads for each sample incorporating ground-truth variants
with correct CIGAR strings for SNPs, insertions, and deletions.
"""

import numpy as np
import pysam

from config import (
    CHROMOSOMES, SAMPLE_GROUPS, READ_LENGTH,
    get_target_coverage, get_sample_group, get_variant_af,
    BASEQ_MEAN, BASEQ_SD, BASEQ_CLIP, MAPQ_RANGE, FRAGMENT_WINDOW,
    PAIRED_END_ENABLED, PAIRED_READ_LENGTH, PAIRED_INSERT_SIZE_MEAN,
    PAIRED_INSERT_SIZE_SD,
)
from reference import apply_variants_to_ref


# ---------------------------------------------------------------------------
# Genotype assignment
# ---------------------------------------------------------------------------

def _assign_genotypes_biallelic(af, n_samples, rng):
    """
    Assign diploid genotypes for a biallelic variant.

    Uses Hardy-Weinberg proportions:
      P(0/0) = (1-af)^2, P(0/1) = 2*af*(1-af), P(1/1) = af^2

    Returns
    -------
    genotypes : list of str
        Each element is "0/0", "0/1", or "1/1".
    """
    p00 = (1 - af) ** 2
    p01 = 2 * af * (1 - af)
    p11 = af ** 2
    # Normalize (handles floating point)
    total = p00 + p01 + p11
    if total < 1e-15:
        # AF is essentially 0
        return ["0/0"] * n_samples
    p00 /= total
    p01 /= total
    p11 /= total

    choices = rng.choice(3, size=n_samples, p=[p00, p01, p11])
    gt_map = {0: "0/0", 1: "0/1", 2: "1/1"}
    return [gt_map[c] for c in choices]


def _assign_genotypes_multiallelic(af_tuple, n_samples, rng):
    """
    Assign diploid genotypes for a multi-allelic variant.

    af_tuple = (af_alt1, af_alt2, ...), ref_af = 1 - sum(af_tuple)

    Genotype probabilities under HWE for 3 alleles (A0=ref, A1, A2):
      P(i/j) = 2 * p_i * p_j  for i != j
      P(i/i) = p_i^2

    Returns list of genotype strings like "0/0", "0/1", "0/2", "1/1", "1/2", "2/2".
    """
    af_list = list(af_tuple)
    ref_af = 1.0 - sum(af_list)
    if ref_af < 0:
        raise ValueError(f"Allele frequencies sum to > 1: {af_tuple}")
    allele_freqs = [ref_af] + af_list  # [p0, p1, p2, ...]
    n_alleles = len(allele_freqs)

    # Build genotype probability table
    genotypes = []
    probs = []
    for i in range(n_alleles):
        for j in range(i, n_alleles):
            gt_str = f"{i}/{j}"
            if i == j:
                p = allele_freqs[i] ** 2
            else:
                p = 2 * allele_freqs[i] * allele_freqs[j]
            genotypes.append(gt_str)
            probs.append(p)

    # Normalize
    total = sum(probs)
    if total < 1e-15:
        return ["0/0"] * n_samples
    probs = [p / total for p in probs]

    indices = rng.choice(len(genotypes), size=n_samples, p=probs)
    return [genotypes[idx] for idx in indices]


def assign_all_genotypes(variants, seed=42):
    """
    Assign genotypes for all variants across all samples.

    Returns
    -------
    sample_genotypes : dict
        {sample_id: {variant_id: genotype_string}}
    """
    rng = np.random.default_rng(seed)
    all_samples = []
    for group_samples in SAMPLE_GROUPS.values():
        all_samples.extend(group_samples)
    n_total = len(all_samples)

    # Build sample -> group mapping
    sample_to_group = {}
    for group, samples in SAMPLE_GROUPS.items():
        for s in samples:
            sample_to_group[s] = group

    # Initialize result
    sample_genotypes = {s: {} for s in all_samples}

    for v in variants:
        vid = v["id"]
        vtype = v["type"]

        # Assign genotypes per group
        group_genotypes = {}  # {group: [gt_for_each_sample_in_group]}
        for group, samples in SAMPLE_GROUPS.items():
            n = len(samples)
            af = v["af"][group]

            if vtype == "multi-allelic":
                gts = _assign_genotypes_multiallelic(af, n, rng)
            else:
                gts = _assign_genotypes_biallelic(af, n, rng)
            group_genotypes[group] = gts

        # Assign to samples
        for group, samples in SAMPLE_GROUPS.items():
            for i, s in enumerate(samples):
                sample_genotypes[s][vid] = group_genotypes[group][i]

    return sample_genotypes


# ---------------------------------------------------------------------------
# Read generation
# ---------------------------------------------------------------------------

def _generate_base_qualities(read_length, rng):
    """Generate array of base qualities (Phred scale)."""
    quals = rng.normal(BASEQ_MEAN, BASEQ_SD, size=read_length)
    quals = np.clip(quals, BASEQ_CLIP[0], BASEQ_CLIP[1]).astype(np.uint8)
    return quals


def _revcomp_seq(seq):
    """Reverse complement a DNA sequence."""
    trans = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
    return seq.translate(trans)[::-1]


def _reverse_cigar(cigar):
    """Reverse a CIGAR list for reverse-strand reads."""
    return list(reversed(cigar))


def _reverse_read(seq, quals, cigar):
    """Reverse-complement read sequence, reverse qualities and reverse CIGAR."""
    seq_rc = _revcomp_seq(seq)
    quals_rc = list(reversed(quals))
    cigar_rc = _reverse_cigar(cigar)
    return seq_rc, quals_rc, cigar_rc


def _generate_read_sequence(ref_seq, chrom_len, start_0based, read_length, rng):
    """
    Extract a reference-matching read sequence from the reference.

    Handles edge cases where read extends beyond chromosome end.
    Returns (sequence_string, actual_length).
    """
    end = min(start_0based + read_length, chrom_len)
    actual_len = end - start_0based
    if actual_len <= 0:
        return "", 0
    return ref_seq[start_0based:end].upper(), actual_len


def _inject_snp(seq, variant_offset, alt_base):
    """
    Inject a SNP into a read sequence.

    variant_offset: position of the variant within the read (0-based).
    alt_base: the alternate allele base (single character).
    """
    seq_list = list(seq)
    if 0 <= variant_offset < len(seq_list):
        seq_list[variant_offset] = alt_base
    return "".join(seq_list)


def _inject_deletion(seq, ref_seq, read_start, variant_pos_0, del_length, read_length):
    """
    Generate a read sequence and CIGAR for a deletion variant.

    The read skips `del_length` bases in the reference at the variant position.

    Returns (new_seq, new_cigar_tuples) or None if variant not in read.
    """
    # variant_pos_0 is 0-based position of deletion start in reference
    offset_in_read = variant_pos_0 - read_start

    if offset_in_read < 0 or offset_in_read >= read_length:
        return None  # Variant not in this read

    # CIGAR: matches before deletion, deletion, matches after
    before_len = offset_in_read
    after_len = read_length - before_len

    if after_len <= 0:
        return None

    # The read sequence: ref bases before deletion + ref bases after deletion
    # (skipping the deleted bases in reference)
    ref_before = ref_seq[read_start:read_start + before_len].upper()
    ref_after_start = variant_pos_0 + del_length
    ref_after_end = min(ref_after_start + after_len, len(ref_seq))
    actual_after = ref_after_end - ref_after_start
    ref_after = ref_seq[ref_after_start:ref_after_end].upper()

    new_seq = ref_before + ref_after
    new_read_len = len(new_seq)

    # CIGAR tuples: (operation, length)
    # 0=M (match), 2=D (deletion)
    cigar = []
    if before_len > 0:
        cigar.append((0, before_len))
    cigar.append((2, del_length))
    if actual_after > 0:
        cigar.append((0, actual_after))

    # Adjust if read was truncated
    if new_read_len < read_length:
        # This happens when the read extends past chromosome end after deletion
        pass

    return new_seq, cigar


def _inject_insertion(seq, ref_seq, read_start, variant_pos_0, inserted_bases, read_length):
    """
    Generate a read sequence and CIGAR for an insertion variant.

    The read has extra bases not in the reference at the variant position.

    Returns (new_seq, new_cigar_tuples) or None if variant not in read.
    """
    offset_in_read = variant_pos_0 - read_start

    if offset_in_read < 0 or offset_in_read >= read_length:
        return None

    ins_len = len(inserted_bases)
    before_len = offset_in_read
    # After insertion, we need (read_length - before_len - ins_len) more ref bases
    after_ref_len = read_length - before_len - ins_len

    if after_ref_len < 0:
        return None  # Insertion too long for this read

    # Build read sequence
    ref_before = ref_seq[read_start:read_start + before_len].upper()
    ref_after_start = variant_pos_0
    ref_after_end = min(ref_after_start + after_ref_len, len(ref_seq))
    ref_after = ref_seq[ref_after_start:ref_after_end].upper()

    new_seq = ref_before + inserted_bases.upper() + ref_after
    actual_after_ref = len(ref_after)

    # CIGAR: matches, insertion, matches
    cigar = []
    if before_len > 0:
        cigar.append((0, before_len))
    cigar.append((1, ins_len))  # 1=I (insertion)
    if actual_after_ref > 0:
        cigar.append((0, actual_after_ref))

    return new_seq, cigar


def generate_reads_for_sample(
    sample_id, ref_seqs, variant_map, sample_genotypes, seed=42,
    paired=None, paired_read_length=None,
    insert_size_mean=None, insert_size_sd=None,
):
    """
    Generate all reads for a single sample.

    Parameters
    ----------
    sample_id : str
    ref_seqs : dict
        {chrom: sequence_string}
    variant_map : dict
        {chrom: [(pos_0based, ref, alt, type, var_id), ...]}
    sample_genotypes : dict
        {sample_id: {var_id: genotype_string}}
    seed : int

    Returns
    -------
    reads : list of (chrom_name, pysam.AlignedSegment)
        All generated reads for this sample, as (chromosome, read) tuples.
    """
    # Use a deterministic seed derived from sample_id
    sample_seed = seed + hash(sample_id) % (2**31)
    rng = np.random.default_rng(sample_seed)

    target_cov = get_target_coverage(sample_id)
    genotypes = sample_genotypes[sample_id]
    all_reads = []  # list of (chrom_name, read) tuples

    if paired is None:
        paired = PAIRED_END_ENABLED
    if paired:
        read_length = paired_read_length or PAIRED_READ_LENGTH
        insert_size_mean = insert_size_mean or PAIRED_INSERT_SIZE_MEAN
        insert_size_sd = insert_size_sd or PAIRED_INSERT_SIZE_SD
    else:
        read_length = READ_LENGTH

    for chrom, chrom_len in CHROMOSOMES.items():
        chrom_seq = ref_seqs[chrom]
        if paired:
            n_reads = max(1, int(target_cov * chrom_len / (2 * read_length)))
        else:
            n_reads = max(1, int(target_cov * chrom_len / READ_LENGTH))

        # Generate uniform start positions
        if paired:
            # Precompute a conservative maximum so both mates fit.
            max_start = max(0, chrom_len - 2 * read_length)
        else:
            max_start = max(0, chrom_len - READ_LENGTH)
        if max_start == 0:
            starts = np.zeros(n_reads, dtype=int)
        else:
            starts = rng.integers(0, max_start + 1, size=n_reads)
        starts.sort()

        # Get variants on this chromosome
        chrom_variants = variant_map.get(chrom, [])

        for start_pos in starts:
            if paired:
                template_len = max(2 * read_length, int(round(rng.normal(insert_size_mean, insert_size_sd))))
                if template_len > chrom_len:
                    continue
                fragment_start = start_pos
                mate_positions = [fragment_start, fragment_start + template_len - read_length]
                query_name = f"{sample_id}_{chrom}_{fragment_start}_{rng.integers(0, 10**9)}"
            else:
                mate_positions = [start_pos]
                query_name = f"{sample_id}_{chrom}_{start_pos}_{rng.integers(0, 10**9)}"

            fragment_reads = []
            for mate_idx, read_start in enumerate(mate_positions, start=1):
                read = pysam.AlignedSegment()
                read.query_name = query_name
                read.flag = 0
                read.reference_id = -1
                read.reference_start = read_start
                read.mapping_quality = int(rng.integers(MAPQ_RANGE[0], MAPQ_RANGE[1]))

                # Generate base qualities
                quals = _generate_base_qualities(read_length, rng)

                # Start with reference sequence
                read_seq, actual_len = _generate_read_sequence(
                    chrom_seq, chrom_len, read_start, read_length, rng
                )
                if actual_len == 0:
                    continue

                read_end = read_start + actual_len
                cigar = [(0, actual_len)]  # Default: all matches
                has_indel = False

                # Check each variant on this chromosome
                for (vpos, vref, valt, vtype, vid) in chrom_variants:
                    gt = genotypes.get(vid, "0/0")

                    # Skip if read doesn't overlap variant
                    if vtype == "SNP" or vtype == "multi-allelic":
                        if not (read_start <= vpos < read_end):
                            continue
                    else:  # indel
                        if not (read_start <= vpos < read_end):
                            continue

                    # Decide if this read carries the alt allele
                    carries_alt = _read_carries_alt(gt, rng, vtype)
                    if not carries_alt:
                        continue

                    # Apply variant to read
                    if vtype == "SNP":
                        offset = vpos - read_start
                        alt_base = valt  # Single base
                        read_seq = _inject_snp(read_seq, offset, alt_base)

                    elif vtype == "multi-allelic":
                        alt_alleles = valt.split(",")
                        alt_idx = _choose_alt_allele(gt, rng)
                        if alt_idx is not None and alt_idx < len(alt_alleles):
                            offset = vpos - read_start
                            alt_base = alt_alleles[alt_idx]
                            read_seq = _inject_snp(read_seq, offset, alt_base)

                    elif vtype == "del":
                        del_len = len(vref) - len(valt)
                        result = _inject_deletion(
                            read_seq, chrom_seq, read_start, vpos, del_len, read_length
                        )
                        if result is not None:
                            read_seq, cigar = result
                            has_indel = True
                            quals = quals[:len(read_seq)]

                    elif vtype == "ins":
                        inserted = valt[len(vref):]
                        result = _inject_insertion(
                            read_seq, chrom_seq, read_start, vpos, inserted, read_length
                        )
                        if result is not None:
                            read_seq, cigar = result
                            has_indel = True
                            extra_quals = _generate_base_qualities(len(inserted), rng)
                            ins_offset = vpos - read_start
                            quals = np.concatenate([
                                quals[:ins_offset + 1],
                                extra_quals,
                                quals[ins_offset + 1:]
                            ])

                if read_seq == "" or len(read_seq) == 0:
                    continue

                if paired and mate_idx == 2:
                    # second mate is reverse-oriented in standard paired-end data
                    read_seq, quals, cigar = _reverse_read(read_seq, quals[:len(read_seq)], cigar)
                    read.is_reverse = True
                    read.is_read2 = True
                    read.is_read1 = False
                else:
                    read.is_reverse = False
                    read.is_read1 = True
                    read.is_read2 = False

                read.query_sequence = read_seq
                read.query_qualities = quals[:len(read_seq)]
                read.cigartuples = cigar
                read.is_paired = paired
                read.is_proper_pair = paired
                read.mate_is_unmapped = False
                read.next_reference_start = 0
                fragment_reads.append((read_start, read_end, read))

            if paired and len(fragment_reads) == 2:
                r1_start, r1_end, r1 = fragment_reads[0]
                r2_start, r2_end, r2 = fragment_reads[1]
                r1.next_reference_start = r2.reference_start
                r2.next_reference_start = r1.reference_start
                r1.mate_is_reverse = r2.is_reverse
                r2.mate_is_reverse = r1.is_reverse
                left = min(r1.reference_start, r2.reference_start)
                right = max(r1_end, r2_end)
                tlen = right - left
                r1.template_length = tlen if r1.reference_start == left else -tlen
                r2.template_length = tlen if r2.reference_start == left else -tlen
                all_reads.append((chrom, r1))
                all_reads.append((chrom, r2))
            elif not paired and fragment_reads:
                all_reads.append((chrom, fragment_reads[0][2]))

    return all_reads


def _read_carries_alt(genotype, rng, vtype="SNP"):
    """
    Decide if a single read carries the alt allele based on genotype.

    For 0/0: never
    For 0/1: 50% chance
    For 1/1: always
    For multi-allelic: depends on genotype
    """
    if genotype == "0/0":
        return False
    elif genotype == "1/1":
        return True
    elif genotype == "0/1":
        return rng.random() < 0.5
    elif genotype == "0/2":
        # Multi-allelic: carries alt2
        return rng.random() < 0.5
    elif genotype == "1/2":
        # Multi-allelic: carries both alt1 and alt2 - for a single read,
        # randomly pick which alt to carry
        return True
    else:
        # Generic: if any non-ref allele present, 50% chance
        return rng.random() < 0.5


def _choose_alt_allele(genotype, rng):
    """
    For multi-allelic variants, choose which alt allele a read carries.

    Returns 0-based index into alt alleles list (0 = first alt, 1 = second alt).
    Returns None if read should carry ref.
    """
    if genotype == "0/0":
        return None
    elif genotype == "0/1":
        return 0 if rng.random() < 0.5 else None
    elif genotype == "0/2":
        return 1 if rng.random() < 0.5 else None
    elif genotype == "1/1":
        return 0
    elif genotype == "2/2":
        return 1
    elif genotype == "1/2":
        # Carry either alt1 or alt2 with equal probability
        return rng.choice([0, 1])
    else:
        return None
