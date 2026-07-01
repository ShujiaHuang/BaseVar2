"""
Ancient DNA simulation module for synthetic test data.

Applies post-mortem damage (PMD) patterns to existing reads,
converting modern-style reads into ancient DNA-like reads.

Key transformations:
  1. Fragment shortening: reads are clipped to short fragments
     sampled from a geometric distribution (mean ~50bp).
  2. C→T deamination: cytosines near fragment termini are converted
     to thymines with geometrically decaying probability.
  3. Contamination: a configurable fraction of reads are left
     undamaged to simulate modern human DNA contamination.

Damage model
------------
The probability of deamination at distance *d* from a fragment end is:

    P(damage at d) = max_rate × (1 - decay) ^ d

This gives the characteristic U-shaped damage profile observed in
real ancient DNA data (Briggs et al., 2007).

Usage
-----
Called from generate.py when --ancient-dna flag is provided.
Can also be used standalone::

    from ancient_dna import simulate_ancient_dna
    ancient_reads = simulate_ancient_dna(modern_reads, seed=42)
"""

import numpy as np
import pysam

from config import (
    ANCIENT_DNA_FRAG_LENGTH_MEAN,
    ANCIENT_DNA_FRAG_LENGTH_MIN,
    ANCIENT_DNA_DAMAGE_RATE,
    ANCIENT_DNA_DAMAGE_DECAY,
    ANCIENT_DNA_CONTAMINATION_RATE,
)


# ---------------------------------------------------------------------------
# CIGAR clipping
# ---------------------------------------------------------------------------

def _clip_cigar_left(cigar, n_query_bases):
    """
    Clip *n_query_bases* query-consuming operations from the left of a CIGAR.

    Query-consuming operations: M(0), I(1), S(4), =(7), X(8).
    Reference-consuming only: D(2), N(3).

    Returns the remaining CIGAR after removing the first *n_query_bases*
    query bases.
    """
    if n_query_bases == 0:
        return list(cigar)

    QUERY_OPS = {0, 1, 4, 7, 8}  # M, I, S, =, X

    remaining = n_query_bases
    new_cigar = []
    for op, length in cigar:
        if remaining <= 0:
            new_cigar.append((op, length))
        elif op in QUERY_OPS:
            if length <= remaining:
                remaining -= length
                # Entire operation consumed; skip it
            else:
                new_cigar.append((op, length - remaining))
                remaining = 0
        else:
            # D, N: do not consume query; keep them
            new_cigar.append((op, length))

    return new_cigar


def _clip_cigar_right(cigar, n_query_bases):
    """
    Clip *n_query_bases* query-consuming operations from the right of a CIGAR.
    """
    if n_query_bases == 0:
        return list(cigar)

    QUERY_OPS = {0, 1, 4, 7, 8}

    remaining = n_query_bases
    new_cigar = []
    for op, length in reversed(cigar):
        if remaining <= 0:
            new_cigar.append((op, length))
        elif op in QUERY_OPS:
            if length <= remaining:
                remaining -= length
            else:
                new_cigar.append((op, length - remaining))
                remaining = 0
        else:
            new_cigar.append((op, length))

    new_cigar.reverse()
    return new_cigar


def _cigar_ref_bases(cigar):
    """Count total reference-consuming bases in a CIGAR list."""
    REF_OPS = {0, 2, 3, 7, 8}
    return sum(length for op, length in cigar if op in REF_OPS)


def _cigar_query_bases(cigar):
    """Count total query-consuming bases in a CIGAR list."""
    QUERY_OPS = {0, 1, 4, 7, 8}
    return sum(length for op, length in cigar if op in QUERY_OPS)


def _cigar_ref_offset_left(cigar, n_query_bases):
    """
    Count reference bases consumed by the first *n_query_bases*
    query-consuming operations in a CIGAR.

    Used to compute the reference_start offset when clipping from the left.
    """
    if n_query_bases == 0:
        return 0

    QUERY_OPS = {0, 1, 4, 7, 8}
    REF_OPS = {0, 2, 3, 7, 8}

    remaining_q = n_query_bases
    ref_offset = 0

    for op, length in cigar:
        if remaining_q <= 0:
            break
        if op in QUERY_OPS:
            if length <= remaining_q:
                remaining_q -= length
                if op in REF_OPS:
                    ref_offset += length
            else:
                consumed_q = remaining_q
                remaining_q = 0
                if op in REF_OPS:
                    ref_offset += consumed_q
        else:
            # D, N: consume reference but not query
            ref_offset += length

    return ref_offset


# ---------------------------------------------------------------------------
# Fragment length sampling
# ---------------------------------------------------------------------------

def _sample_fragment_length(rng, mean_length, min_length):
    """
    Sample a fragment length from a geometric distribution.

    The geometric distribution is parameterised so that the expected
    value equals *mean_length*.  Values below *min_length* are
    resampled to avoid unrealistically short fragments.
    """
    p = 1.0 / mean_length
    while True:
        length = rng.geometric(p)
        if length >= min_length:
            return int(length)


# ---------------------------------------------------------------------------
# PMD damage application
# ---------------------------------------------------------------------------

def _apply_pmd_damage(seq_array, frag_len, damage_rate, damage_decay, rng):
    """
    Apply C→T deamination damage to a fragment sequence (in-place).

    Damage is applied symmetrically at both termini:
      - Left end:  C → T  (mimics deamination on the 5' strand)
      - Right end: C → T  (mimics deamination on the 3' strand,
                            which appears as C→T on the forward strand)

    The probability of damage at distance *d* from the nearest terminus is::

        P(d) = damage_rate × (1 - damage_decay) ^ d

    Parameters
    ----------
    seq_array : list of str
        Single-character list of the fragment sequence (modified in-place).
    frag_len : int
        Length of the fragment.
    damage_rate : float
        Maximum damage probability at the terminal position (0-1).
    damage_decay : float
        Geometric decay parameter (0-1).  Higher values → faster decay.
    rng : numpy.random.Generator
    """
    if frag_len == 0:
        return

    half = frag_len // 2

    for d in range(min(half, frag_len)):
        prob = damage_rate * ((1.0 - damage_decay) ** d)
        if rng.random() < prob:
            # Left end: C → T
            idx = d
            if seq_array[idx] == 'C':
                seq_array[idx] = 'T'
            # Right end: C → T (symmetric)
            idx = frag_len - 1 - d
            if seq_array[idx] == 'C':
                seq_array[idx] = 'T'


# ---------------------------------------------------------------------------
# Main simulation
# ---------------------------------------------------------------------------

def simulate_ancient_dna(all_sample_reads, seed=42):
    """
    Transform modern-style reads into ancient DNA-like reads.

    For each sample, every read is:
      1. Randomly fragmented to a short length (geometric distribution).
      2. Subjected to C→T deamination at both termini.
      3. Optionally left undamaged (contamination fraction).

    Parameters
    ----------
    all_sample_reads : dict
        {sample_id: [(chrom_name, pysam.AlignedSegment), ...]}
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    ancient_reads : dict
        {sample_id: [(chrom_name, pysam.AlignedSegment), ...]}
    """
    rng = np.random.default_rng(seed)

    frag_mean = ANCIENT_DNA_FRAG_LENGTH_MEAN
    frag_min = ANCIENT_DNA_FRAG_LENGTH_MIN
    damage_rate = ANCIENT_DNA_DAMAGE_RATE
    damage_decay = ANCIENT_DNA_DAMAGE_DECAY
    contam_rate = ANCIENT_DNA_CONTAMINATION_RATE

    ancient_reads = {}

    for sample_id, reads in all_sample_reads.items():
        sample_ancient = []
        n_total = len(reads)
        n_damaged = 0
        n_contaminated = 0
        n_skipped = 0

        for chrom, read in reads:
            orig_len = read.query_length

            # Sample target fragment length
            target_len = _sample_fragment_length(rng, frag_mean, frag_min)

            # Skip reads shorter than target (shouldn't happen with
            # READ_LENGTH=100 and frag_mean=50, but be safe)
            if orig_len <= target_len:
                # Read is already short enough; still apply damage
                target_len = orig_len

            if target_len < frag_min:
                n_skipped += 1
                continue

            # Random left clip amount
            max_clip = orig_len - target_len
            left_clip = int(rng.integers(0, max_clip + 1)) if max_clip > 0 else 0
            right_clip = orig_len - target_len - left_clip

            # Extract sub-fragment sequence and qualities
            seq = read.query_sequence
            quals = list(read.query_qualities) if read.query_qualities is not None else [30] * orig_len
            orig_cigar = list(read.cigartuples) if read.cigartuples else [(0, orig_len)]

            # Compute reference offset for the left clip
            ref_offset = _cigar_ref_offset_left(orig_cigar, left_clip)

            # Clip query sequence from left and right
            frag_end = orig_len - right_clip
            seq = seq[left_clip:frag_end]
            quals = quals[left_clip:frag_end]

            # Clip CIGAR: remove left_clip query bases from left,
            # right_clip query bases from right
            inner_cigar = _clip_cigar_left(orig_cigar, left_clip)
            inner_cigar = _clip_cigar_right(inner_cigar, right_clip)
            inner_cigar = [op for op in inner_cigar if op[1] > 0]

            frag_len = len(seq)
            if frag_len < frag_min:
                n_skipped += 1
                continue

            # Decide if this read is contaminated (undamaged)
            is_contaminated = rng.random() < contam_rate

            if not is_contaminated:
                # Apply PMD damage
                seq_array = list(seq)
                _apply_pmd_damage(
                    seq_array, frag_len,
                    damage_rate, damage_decay, rng
                )
                seq = "".join(seq_array)
                n_damaged += 1
            else:
                n_contaminated += 1

            # Build new AlignedSegment
            new_read = pysam.AlignedSegment()
            new_read.query_name = read.query_name
            new_read.flag = read.flag
            new_read.reference_id = read.reference_id
            new_read.reference_start = read.reference_start + ref_offset
            new_read.mapping_quality = read.mapping_quality
            new_read.query_sequence = seq
            new_read.query_qualities = np.array(
                quals[:frag_len], dtype=np.uint8
            )

            # Build CIGAR: hard-clip the removed bases (H removes
            # from both CIGAR and query sequence, matching our
            # already-trimmed sequence).
            # Note: BAM hard clips (5) are only stored in the CIGAR;
            # the query_sequence does NOT include them.
            new_cigar = list(inner_cigar)
            new_read.cigartuples = new_cigar

            # Copy optional tags (e.g. RG)
            for tag in read.tags:
                if tag[0] not in ('RG',):  # RG will be set by writer
                    try:
                        new_read.set_tag(tag[0], tag[1])
                    except Exception:
                        pass

            sample_ancient.append((chrom, new_read))

        ancient_reads[sample_id] = sample_ancient
        print(
            f"    {sample_id}: {n_total} reads → "
            f"{len(sample_ancient)} fragments "
            f"({n_damaged} damaged, {n_contaminated} contaminated, "
            f"{n_skipped} skipped)"
        )

    return ancient_reads
