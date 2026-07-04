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

    P(damage at d) = max_rate x (1 - decay) ^ d

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
import config as synthetic_config

from config import (
    ANCIENT_DNA_FRAG_LENGTH_MEAN,
    ANCIENT_DNA_FRAG_LENGTH_MIN,
    ANCIENT_DNA_DAMAGE_RATE,
    ANCIENT_DNA_DAMAGE_DECAY,
    ANCIENT_DNA_CONTAMINATION_RATE,
    ANCIENT_DNA_STRICT_TAGS,
    ANCIENT_DNA_Q_DECAY,
    ANCIENT_DNA_SS_DNA,
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

        P(d) = damage_rate x (1 - damage_decay) ^ d

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

    # Apply damage symmetrically by default. The calling code must decide
    # whether the read was aligned in reverse orientation and pass a flag
    # (handled in the caller via applying substitutions appropriate for
    # the read orientation). This helper performs substitutions on the
    # provided seq_array in-place assuming substitutions are C->T.

    half = frag_len // 2

    for d in range(min(half, frag_len)):
        prob = damage_rate * ((1.0 - damage_decay) ** d)
        if rng.random() < prob:
            # Left end: substitution target at index d
            idx = d
            if seq_array[idx] == 'C':
                seq_array[idx] = 'T'
            # Right end: symmetric position
            idx = frag_len - 1 - d
            if seq_array[idx] == 'C':
                seq_array[idx] = 'T'


# ---------------------------------------------------------------------------
# Main simulation
# ---------------------------------------------------------------------------

def _revcomp(seq: str) -> str:
    trans = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
    return seq.translate(trans)[::-1]


def _compute_md_nm_and_set_tags(read, ref_seq_for_read):
    """
    Compute MD string and NM (edit distance) for an AlignedSegment and set tags.

    Parameters
    - read: pysam.AlignedSegment (must have reference_start/cigartuples/query_sequence)
    - ref_seq_for_read: string of reference bases covering the alignment region
    """
    cigar = read.cigartuples or [(0, len(read.query_sequence))]
    # Sequence as aligned to reference: if read.is_reverse, use reverse-complement
    seq = read.query_sequence
    is_rev = read.is_reverse
    aligned_seq = _revcomp(seq) if is_rev else seq

    ref_idx = 0
    seq_idx = 0
    md_parts = []
    matches = 0
    mismatches = 0
    indel_bases = 0

    for op, length in cigar:
        if op in (0, 7, 8):  # M, =, X
            for i in range(length):
                ref_base = ref_seq_for_read[ref_idx].upper()
                read_base = aligned_seq[seq_idx].upper()
                if read_base == ref_base:
                    matches += 1
                else:
                    if matches > 0:
                        md_parts.append(str(matches))
                        matches = 0
                    md_parts.append(ref_base)
                    mismatches += 1
                ref_idx += 1
                seq_idx += 1
        elif op == 1:  # I
            # Insertion w.r.t reference: consumes query only
            seq_idx += length
            indel_bases += length
        elif op == 2:  # D
            if matches > 0:
                md_parts.append(str(matches))
                matches = 0
            deleted = ref_seq_for_read[ref_idx:ref_idx + length].upper()
            md_parts.append('^' + deleted)
            ref_idx += length
            indel_bases += length
        elif op == 3:  # N
            # skipped region from the reference
            ref_idx += length
        elif op == 4:  # S
            seq_idx += length
        elif op == 5:  # H
            # hard clip: nothing to advance in seq
            continue
        else:
            # Unknown op: be conservative
            for i in range(length):
                ref_idx += 1
                seq_idx += 1

    if matches > 0:
        md_parts.append(str(matches))

    md_str = ''.join(md_parts) if md_parts else '0'
    nm = mismatches + indel_bases
    try:
        read.set_tag('NM', int(nm), value_type='i')
    except Exception:
        pass
    try:
        read.set_tag('MD', md_str)
    except Exception:
        pass


def simulate_ancient_dna(all_sample_reads, ref_seqs=None, seed=42):
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

    # If a contaminant BAM is provided, preload a small reservoir of contaminant reads
    contaminant_pool = None
    try:
        from config import ANCIENT_DNA_CONTAMINANT_BAM
    except Exception:
        ANCIENT_DNA_CONTAMINANT_BAM = None
    if ANCIENT_DNA_CONTAMINANT_BAM:
        try:
            cbf = pysam.AlignmentFile(ANCIENT_DNA_CONTAMINANT_BAM, 'rb')
            contaminant_pool = []
            for i, r in enumerate(cbf.fetch(until_eof=True)):
                if r.query_length is None:
                    continue
                contaminant_pool.append((r.reference_name, r))
                if i >= 10000:
                    break
            cbf.close()
        except Exception:
            contaminant_pool = None

    # If a microbial FASTA is provided, load sequences into a pool
    microbe_pool = None
    microbe_fasta = getattr(synthetic_config, 'ANCIENT_DNA_MICROBE_FASTA', None)
    microbe_rate = getattr(synthetic_config, 'ANCIENT_DNA_MICROBE_RATE', 0.0)
    if microbe_fasta:
        try:
            mfa = pysam.FastaFile(microbe_fasta)
            microbe_pool = []
            for name in mfa.references:
                seq = mfa.fetch(name).upper()
                if len(seq) > 0:
                    microbe_pool.append((name, seq))
            mfa.close()
        except Exception:
            microbe_pool = None

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
                # Apply PMD damage: handle orientation-sensitive substitutions
                seq_array = list(seq)
                if read.is_reverse:
                    # On reverse-oriented reads, C->T damage on the reference
                    # is observed as G->A on the read. Convert G->A instead.
                    # Work on the reverse-complemented sequence so that
                    # _apply_pmd_damage can operate on C->T positions.
                    rc = _revcomp(seq)
                    rc_array = list(rc)
                    _apply_pmd_damage(rc_array, frag_len, damage_rate, damage_decay, rng)
                    # Convert back to read orientation
                    seq = _revcomp(''.join(rc_array))
                else:
                    _apply_pmd_damage(seq_array, frag_len, damage_rate, damage_decay, rng)
                    seq = ''.join(seq_array)
                # Optionally reduce base quality at damaged positions
                # (simple model: subtract constant from damaged positions)
                if ANCIENT_DNA_Q_DECAY > 0:
                    # Naive approach: reduce all terminal bases by Q_DECAY
                    for i in range(min(5, frag_len)):
                        quals[i] = max(0, quals[i] - ANCIENT_DNA_Q_DECAY)
                        quals[-1 - i] = max(0, quals[-1 - i] - ANCIENT_DNA_Q_DECAY)
                n_damaged += 1
            else:
                n_contaminated += 1
                # If we have a contaminant pool, sample and use it
                if contaminant_pool:
                    chrom_c, r_c = rng.choice(contaminant_pool)
                    new_read = pysam.AlignedSegment()
                    new_read.query_name = r_c.query_name
                    new_read.flag = r_c.flag
                    new_read.reference_id = r_c.reference_id
                    new_read.reference_start = r_c.reference_start
                    new_read.mapping_quality = r_c.mapping_quality
                    new_read.query_sequence = r_c.query_sequence
                    new_read.query_qualities = list(r_c.query_qualities) if r_c.query_qualities is not None else None
                    new_read.cigartuples = list(r_c.cigartuples) if r_c.cigartuples else None
                    try:
                        new_read.set_tag('CT', 1)
                    except Exception:
                        pass
                    sample_ancient.append((chrom_c, new_read))
                    continue

            # Build new AlignedSegment
            new_read = pysam.AlignedSegment()
            new_read.query_name = read.query_name
            new_read.flag = read.flag
            new_read.reference_id = read.reference_id
            new_read.reference_start = read.reference_start + ref_offset
            new_read.mapping_quality = read.mapping_quality
            new_read.query_sequence = seq
            new_read.query_qualities = quals[:frag_len]

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

            # Recompute alignment tags if requested and reference is provided
            if ANCIENT_DNA_STRICT_TAGS and ref_seqs is not None:
                # Extract reference substring covering the alignment
                ref_seq = ref_seqs.get(chrom)
                if ref_seq is not None:
                    ref_len_needed = _cigar_ref_bases(new_cigar)
                    ref_start = new_read.reference_start
                    ref_sub = ref_seq[ref_start:ref_start + ref_len_needed]
                    _compute_md_nm_and_set_tags(new_read, ref_sub)

        # Insert microbial background reads (sampled from provided FASTA)
        if microbe_pool and ref_seqs is not None and microbe_rate > 0:
            try:
                n_microbe = int(round(len(sample_ancient) * microbe_rate))
            except Exception:
                n_microbe = 0
            ref_chroms = list(ref_seqs.keys())
            for mi in range(n_microbe):
                mname, mseq = rng.choice(microbe_pool)
                mlen = len(mseq)
                frag_m = _sample_fragment_length(rng, frag_mean, frag_min)
                if frag_m > mlen:
                    continue
                mstart = int(rng.integers(0, mlen - frag_m + 1))
                subseq = mseq[mstart:mstart + frag_m]

                host_chrom = rng.choice(ref_chroms)
                host_len = len(ref_seqs[host_chrom])
                if frag_m > host_len:
                    continue
                host_pos = int(rng.integers(0, host_len - frag_m + 1))

                mread = pysam.AlignedSegment()
                mread.query_name = f"{sample_id}_microbe_{mi}_{rng.integers(0,10**9)}"
                mread.flag = 0
                mread.reference_id = -1
                mread.reference_start = host_pos
                mread.mapping_quality = 0
                mread.query_sequence = subseq
                mread.query_qualities = [10] * frag_m
                mread.cigartuples = [(0, frag_m)]
                try:
                    mread.set_tag('MB', 1)
                except Exception:
                    pass

                # Compute strict tags if requested
                if ANCIENT_DNA_STRICT_TAGS and ref_seqs is not None:
                    ref_sub = ref_seqs[host_chrom][host_pos:host_pos + frag_m]
                    _compute_md_nm_and_set_tags(mread, ref_sub)

                sample_ancient.append((host_chrom, mread))

        ancient_reads[sample_id] = sample_ancient
        print(
            f"    {sample_id}: {n_total} reads → "
            f"{len(sample_ancient)} fragments "
            f"({n_damaged} damaged, {n_contaminated} contaminated, "
            f"{n_skipped} skipped)"
        )

    return ancient_reads
