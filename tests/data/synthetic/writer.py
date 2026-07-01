"""
BAM/CRAM writer and indexer.

Writes reads to sorted, indexed BAM and CRAM files with proper
headers including @SQ and @RG records.
"""

import os
import tempfile
import pysam

from config import CHROMOSOMES, SAMPLE_GROUPS


def _build_header(chromosomes):
    """
    Build a BAM header dict with @SQ records for all chromosomes.

    Returns a header dict suitable for pysam.AlignmentFile.
    """
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [],
    }
    for chrom, length in chromosomes.items():
        header["SQ"].append({"SN": chrom, "LN": length})
    return header


def _add_rg_to_header(header, sample_id):
    """Add a @RG (read group) record to the header."""
    header = dict(header)  # shallow copy
    header["RG"] = [{"ID": sample_id, "SM": sample_id}]
    return header


def _make_chrom_to_tid(header):
    """Build a mapping from chromosome name to reference_id."""
    chrom_to_tid = {}
    for i, sq in enumerate(header["SQ"]):
        chrom_to_tid[sq["SN"]] = i
    return chrom_to_tid


def write_bam(reads_with_chrom, sample_id, output_dir, chromosomes=CHROMOSOMES):
    """
    Write reads to a sorted, indexed BAM file.

    Parameters
    ----------
    reads_with_chrom : list of (chrom_name, pysam.AlignedSegment)
        Reads for this sample with their chromosome names.
    sample_id : str
    output_dir : str
        Directory for BAM output.
    chromosomes : dict
        {chrom_name: length} for @SQ header.

    Returns
    -------
    bam_path : str
        Path to the output BAM file.
    """
    os.makedirs(output_dir, exist_ok=True)
    bam_path = os.path.join(output_dir, f"{sample_id}.bam")

    header = _build_header(chromosomes)
    header = _add_rg_to_header(header, sample_id)
    chrom_to_tid = _make_chrom_to_tid(header)

    # Write unsorted BAM first
    tmp_path = bam_path + ".tmp.bam"
    with pysam.AlignmentFile(tmp_path, "wb", header=header) as outf:
        for chrom_name, read in reads_with_chrom:
            read.reference_id = chrom_to_tid[chrom_name]
            read.set_tag("RG", sample_id)
            outf.write(read)

    # Sort
    sorted_path = bam_path + ".sorted.bam"
    pysam.sort("-o", sorted_path, tmp_path)
    os.remove(tmp_path)

    # Move sorted to final
    os.rename(sorted_path, bam_path)

    # Index
    pysam.index(bam_path)

    return bam_path


def write_cram(reads_with_chrom, sample_id, output_dir, fa_path, chromosomes=CHROMOSOMES):
    """
    Write reads to a sorted, indexed CRAM file.

    Parameters
    ----------
    reads_with_chrom : list of (chrom_name, pysam.AlignedSegment)
    sample_id : str
    output_dir : str
    fa_path : str
        Path to reference FASTA (required for CRAM).
    chromosomes : dict

    Returns
    -------
    cram_path : str
    """
    os.makedirs(output_dir, exist_ok=True)
    cram_path = os.path.join(output_dir, f"{sample_id}.cram")

    header = _build_header(chromosomes)
    header = _add_rg_to_header(header, sample_id)
    chrom_to_tid = _make_chrom_to_tid(header)

    # Write unsorted CRAM first (use temp BAM as intermediate)
    tmp_bam = cram_path + ".tmp.bam"
    with pysam.AlignmentFile(tmp_bam, "wb", header=header) as outf:
        for chrom_name, read in reads_with_chrom:
            read.reference_id = chrom_to_tid[chrom_name]
            read.set_tag("RG", sample_id)
            outf.write(read)

    # Sort to CRAM
    sorted_cram = cram_path + ".sorted.cram"
    pysam.sort(
        "-O", "cram",
        "--reference", fa_path,
        "-o", sorted_cram,
        tmp_bam
    )
    os.remove(tmp_bam)
    os.rename(sorted_cram, cram_path)

    # Index
    pysam.index("-c", cram_path)

    return cram_path


def write_all_samples(all_sample_reads, bam_dir, cram_dir, fa_path, chromosomes=CHROMOSOMES):
    """
    Write BAM and CRAM files for all samples.

    Parameters
    ----------
    all_sample_reads : dict
        {sample_id: [(chrom_name, read), ...]}
    bam_dir : str
    cram_dir : str
    fa_path : str
    chromosomes : dict

    Returns
    -------
    bam_paths : dict
        {sample_id: bam_path}
    cram_paths : dict
        {sample_id: cram_path}
    """
    bam_paths = {}
    cram_paths = {}

    for sample_id, reads in all_sample_reads.items():
        print(f"  [writer] Writing BAM: {sample_id} ({len(reads)} reads)")
        bam_paths[sample_id] = write_bam(reads, sample_id, bam_dir, chromosomes)

        print(f"  [writer] Writing CRAM: {sample_id} ({len(reads)} reads)")
        cram_paths[sample_id] = write_cram(reads, sample_id, cram_dir, fa_path, chromosomes)

    return bam_paths, cram_paths
