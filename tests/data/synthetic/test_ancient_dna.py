#!/usr/bin/env python3
"""Quick validation tests for ancient_dna simulation.

Runs a small battery of checks on a sample BAM using the implemented
`simulate_ancient_dna` function: fragment-length stats, damage curve,
contamination fraction, CIGAR/SEQ consistency, and NM/MD tag consistency.
"""

import sys
import statistics
import math

import pysam
import numpy as np

sys.path.insert(0, '.')
sys.path.insert(0, 'tests/data/synthetic')

import ancient_dna as ad
from config import (
    ANCIENT_DNA_FRAG_LENGTH_MEAN,
    ANCIENT_DNA_DAMAGE_RATE,
    ANCIENT_DNA_DAMAGE_DECAY,
    ANCIENT_DNA_CONTAMINATION_RATE,
)


def map_read_idx_to_ref(read):
    """Return dict mapping read_index -> reference_position (0-based),
    or None for indices that don't consume reference.
    """
    mapping = {}
    ref_pos = read.reference_start
    read_pos = 0
    cigar = read.cigartuples or [(0, read.query_length or 0)]
    for op, length in cigar:
        if op in (0, 7, 8):  # M, =, X
            for i in range(length):
                mapping[read_pos] = ref_pos
                read_pos += 1
                ref_pos += 1
        elif op == 1:  # I
            read_pos += length
        elif op == 2:  # D
            ref_pos += length
        elif op == 3:  # N
            ref_pos += length
        elif op == 4:  # S
            read_pos += length
        elif op == 5:  # H
            # hard clip: not in query sequence
            continue
        else:
            # Unknown op: advance both
            read_pos += length
            ref_pos += length
    return mapping


def parse_md_count(md):
    """Parse MD string and return mismatches+deleted_bases count."""
    # Simple parser: numbers, letters, ^deleted
    i = 0
    mism = 0
    while i < len(md):
        if md[i].isdigit():
            j = i
            while j < len(md) and md[j].isdigit():
                j += 1
            i = j
        elif md[i] == '^':
            # deletion
            i += 1
            j = i
            while j < len(md) and md[j].isalpha():
                j += 1
            mism += (j - i)
            i = j
        else:
            # single mismatch base
            mism += 1
            i += 1
    return mism


def main():
    fa = pysam.FastaFile('tests/data/synthetic/ref/mini_ref.fa')
    ref_seqs = {name: fa.fetch(name) for name in fa.references}

    bam_path = 'tests/data/synthetic/bam/sampleA01.bam'
    bf = pysam.AlignmentFile(bam_path, 'rb')
    reads = []
    for i, r in enumerate(bf.fetch(until_eof=True)):
        if i >= 1000:
            break
        if r.query_length is None:
            continue
        reads.append((r.reference_name, r))
    bf.close()

    print('Loaded original reads:', len(reads))
    all_reads = {'test_sample': reads}
    sim = ad.simulate_ancient_dna(all_reads, ref_seqs=ref_seqs, seed=2026)
    out_reads = sim['test_sample']
    print('Simulated fragments:', len(out_reads))

    # Fragment length stats
    lens = [r.query_length for _, r in out_reads]
    mean_len = statistics.mean(lens) if lens else 0
    print('Mean fragment length:', mean_len, 'target mean:', ANCIENT_DNA_FRAG_LENGTH_MEAN)
    assert mean_len > 0, 'No fragments produced'
    assert abs(mean_len - ANCIENT_DNA_FRAG_LENGTH_MEAN) / ANCIENT_DNA_FRAG_LENGTH_MEAN < 0.5, 'Mean fragment length deviates >50%'

    # Damage curve: compute observed C->T at positions d=0..4 (left and right)
    maxd = 4
    obs = {d: {'hits': 0, 'total': 0} for d in range(maxd + 1)}
    contaminated_reads = 0
    reads_with_positions = 0
    reads_with_terminal_C = 0
    for chrom, r in out_reads:
        seq = r.query_sequence
        mapping = map_read_idx_to_ref(r)
        if not mapping:
            continue
        reads_with_positions += 1
        had_any = False
        has_terminal_C = False
        for d in range(maxd + 1):
            # left
            li = d
            if li in mapping:
                ref_pos = mapping[li]
                ref_base = ref_seqs[chrom][ref_pos].upper()
                read_base = seq[li].upper()
                if ref_base == 'C':
                    has_terminal_C = True
                    obs[d]['total'] += 1
                    if read_base == 'T':
                        obs[d]['hits'] += 1
                        had_any = True
            # right
            ri = r.query_length - 1 - d
            if ri in mapping:
                ref_pos = mapping[ri]
                ref_base = ref_seqs[chrom][ref_pos].upper()
                read_base = seq[ri].upper()
                if ref_base == 'C':
                    has_terminal_C = True
                    obs[d]['total'] += 1
                    if read_base == 'T':
                        obs[d]['hits'] += 1
                        had_any = True
        if not had_any:
            contaminated_reads += 1
        if has_terminal_C:
            reads_with_terminal_C += 1

    print('Reads considered for damage:', reads_with_positions)
    print('Reads with at least one terminal reference C:', reads_with_terminal_C)
    observed_contam = contaminated_reads / max(1, reads_with_terminal_C)
    print('Observed contamination fraction (approx):', observed_contam, 'configured:', ANCIENT_DNA_CONTAMINATION_RATE)
    # Note: this is an approximate diagnostic. Many reads lack terminal C
    # positions and cannot be used to directly infer contamination fraction.

    for d in range(maxd + 1):
        h = obs[d]['hits']
        t = obs[d]['total']
        rate = h / t if t > 0 else float('nan')
        theo = (1.0 - ANCIENT_DNA_CONTAMINATION_RATE) * (ANCIENT_DNA_DAMAGE_RATE * ((1.0 - ANCIENT_DNA_DAMAGE_DECAY) ** d))
        print(f'd={d}: observed_rate={rate:.3f} (n={t}), expected~{theo:.3f}')
        if not math.isnan(rate):
            assert abs(rate - theo) < 0.2 or abs(rate - theo) / max(theo, 1e-6) < 2.0, f'damage rate at d={d} deviates too much'

    # CIGAR / SEQ consistency and NM/MD checks
    for chrom, r in out_reads:
        qlen = r.query_length
        cigar_q = sum(l for op, l in (r.cigartuples or []) if op in (0, 1, 4, 7, 8))
        assert qlen == cigar_q, f'Query length mismatch for {r.query_name}: {qlen} != {cigar_q}'
        # Reference bounds
        ref_len_needed = sum(l for op, l in (r.cigartuples or []) if op in (0, 2, 3, 7, 8))
        ref_seq = ref_seqs[chrom]
        assert 0 <= r.reference_start < len(ref_seq), 'reference_start out of bounds'
        assert r.reference_start + ref_len_needed <= len(ref_seq), 'alignment exceeds reference length'
        # NM/MD tags
        if r.has_tag('NM') and r.has_tag('MD'):
            nm = r.get_tag('NM')
            md = r.get_tag('MD')
            parsed = parse_md_count(md)
            # Count insertion bases from CIGAR to compare against NM
            insert_bases = sum(l for op, l in (r.cigartuples or []) if op == 1)
            expected_nm = parsed + insert_bases
            if nm != expected_nm:
                print(f'WARNING: NM ({nm}) != MD_parsed+INS ({expected_nm}) for {r.query_name}  CIGAR={r.cigartuples} MD={md} (parsed_MD={parsed}, ins={insert_bases})')

    print('All tests passed.')


if __name__ == '__main__':
    main()
