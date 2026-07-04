#!/usr/bin/env python3
"""Validation script for paired-end synthetic read generation."""

import os
import sys
import tempfile

sys.path.insert(0, '.')
sys.path.insert(0, 'tests/data/synthetic')

import pysam
from collections import defaultdict

from config import CHROMOSOMES
from reads_simulator import generate_reads_for_sample
from writer import _build_header, _add_rg_to_header, _make_chrom_to_tid


def main():
    ref_seqs = {chrom: 'A' * length for chrom, length in CHROMOSOMES.items()}
    sample_genotypes = {'sampleA01': {}}

    reads = generate_reads_for_sample(
        'sampleA01', ref_seqs, {}, sample_genotypes,
        seed=42, paired=True, paired_read_length=50,
        insert_size_mean=150, insert_size_sd=10,
    )
    assert len(reads) % 2 == 0, 'Expected even number of paired reads'

    # Build header and write to a temp BAM for validation of flags and mate fields.
    header = _build_header(CHROMOSOMES)
    header = _add_rg_to_header(header, 'sampleA01')
    chrom_to_tid = _make_chrom_to_tid(header)

    with tempfile.TemporaryDirectory() as tmpdir:
        bam_path = os.path.join(tmpdir, 'paired_test.bam')
        with pysam.AlignmentFile(bam_path, 'wb', header=header) as outf:
            for chrom_name, read in reads:
                read.reference_id = chrom_to_tid[chrom_name]
                if read.is_paired and not getattr(read, 'mate_is_unmapped', False):
                    read.next_reference_id = chrom_to_tid[chrom_name]
                read.set_tag('RG', 'sampleA01')
                outf.write(read)

        pair_info = defaultdict(list)
        with pysam.AlignmentFile(bam_path, 'rb') as bf:
            for r in bf.fetch(until_eof=True):
                pair_info[r.query_name].append(r)

    assert all(len(v) == 2 for v in pair_info.values()), 'Every fragment should produce exactly two reads'

    for qname, pair in pair_info.items():
        r1, r2 = sorted(pair, key=lambda r: not r.is_read1)
        assert r1.is_read1 and r2.is_read2, f'Pair orientation wrong for {qname}'
        assert r1.is_paired and r2.is_paired
        assert r1.is_reverse != r2.is_reverse, f'Expected opposite strand orientation for {qname}'
        assert r1.next_reference_start == r2.reference_start
        assert r2.next_reference_start == r1.reference_start
        assert abs(r1.template_length) == abs(r2.template_length)
        assert r1.template_length == -r2.template_length
        assert r1.reference_start != r2.reference_start

    print('Paired-end validation passed.')


if __name__ == '__main__':
    main()
