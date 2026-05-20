"""A pipeline generator for BaseVar caller.
Splits the genome into regions and generates basevar caller commands for each.

Pipeline-specific options (handled by this script):
  --outdir, --ref_fai, --delta, --chrom

All other options are passed through directly to 'basevar caller'.
This design ensures that any new parameters added to basevar caller
will automatically be supported without modifying this script.

Author: Shujia Huang
Date: 2017-01-19
"""
import os
import sys
import argparse
import shutil


def load_reference_fai(in_fai, chroms=None):

    ref = []
    with open(in_fai) as fh:

        for r in fh:
            # chr2  243199373  254235646  50 51
            col = r.strip().split()
            if chroms is not None and len(chroms):
                if col[0] in chroms:
                    ref.append([col[0], 1, int(col[1])])
            else:
                ref.append([col[0], 1, int(col[1])])

    return ref


def parse_region_string(region_str, ref_fai_map):
    """Parse a single region string into (chrom, start, end) tuple.

    Supported formats:
      - 'chr'              => (chr, 1, chrom_length)
      - 'chr:start'        => (chr, start, chrom_length)
      - 'chr:start-end'    => (chr, start, end)

    Args:
        region_str: A single region string (no commas).
        ref_fai_map: dict of {chrom: length} built from the .fai file.

    Returns:
        (chrom, start, end) tuple (1-based, inclusive).

    Raises:
        ValueError: If the region string is malformed or chrom not found in fai.
    """
    region_str = region_str.strip()
    if not region_str:
        raise ValueError("Empty region string.")

    parts = region_str.split(':', 1)
    chrom = parts[0]

    if chrom not in ref_fai_map:
        raise ValueError(f"Chromosome '{chrom}' not found in reference fai.")

    chrom_len = ref_fai_map[chrom]

    if len(parts) == 1:
        # Format: 'chr'
        return (chrom, 1, chrom_len)

    coord_str = parts[1]
    coord_parts = coord_str.split('-', 1)

    try:
        start = int(coord_parts[0])
    except ValueError:
        raise ValueError(f"Invalid start position in region '{region_str}'.")

    if len(coord_parts) == 1:
        # Format: 'chr:start'
        end = chrom_len
    else:
        # Format: 'chr:start-end'
        try:
            end = int(coord_parts[1])
        except ValueError:
            raise ValueError(f"Invalid end position in region '{region_str}'.")

    if start < 1:
        raise ValueError(f"Start position must be >= 1 in region '{region_str}'.")
    if start > end:
        raise ValueError(f"Start position > end position in region '{region_str}'.")
    if end > chrom_len:
        raise ValueError(f"End position {end} exceeds chromosome length {chrom_len} "
                         f"for '{chrom}' in region '{region_str}'.")

    return (chrom, start, end)


def parse_regions(regions_str, ref_fai_map, chroms_filter=None):
    """Parse a comma-separated regions string into a list of (chrom, start, end).

    Args:
        regions_str:   Comma-separated region strings, e.g. 'chr1,chr2:1000-2000'.
        ref_fai_map:   dict of {chrom: length}.
        chroms_filter: Optional list of chromosomes to keep; None means keep all.

    Returns:
        List of (chrom, start, end) tuples.
    """
    results = []
    for token in regions_str.split(','):
        token = token.strip()
        if not token:
            continue
        chrom, start, end = parse_region_string(token, ref_fai_map)
        if chroms_filter and chrom not in chroms_filter:
            continue
        results.append((chrom, start, end))
    return results


def executable(cmd):
    """Return True if `cmd` (a command string or path) refers to an executable.
    If `cmd` contains arguments, only the program name (first token) is checked.
    """
    if not cmd:
        return False
    prog = cmd.strip().split()[0]
    # If prog contains a path separator or is a relative path, check directly
    if os.path.sep in prog or prog.startswith('.'):
        return os.path.exists(prog) and os.access(prog, os.X_OK)
        
    # Otherwise use shutil.which to find it on PATH
    return shutil.which(prog) is not None


def creat_basetype_pipe():
    # get basevar from enviroment
    basevar = os.environ.get('basevar')
    if basevar:
        exe_prog = basevar + ' caller'
    else:
        pardir = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), os.path.pardir))
        exe_prog = os.path.join(pardir, 'bin', 'basevar caller')
    
    if not executable(exe_prog):
        print(f'Error: {exe_prog} is not executable', file=sys.stderr)
        sys.exit(1)

    # Only define pipeline-specific options here.
    # All other options are passed through to 'basevar caller' as-is.
    optp = argparse.ArgumentParser(
        description="Generate basevar caller pipeline commands by splitting genome into regions.",
        epilog="Any unrecognized options are passed through directly to 'basevar caller'. "
               "See 'basevar caller --help' for the full list of supported options."
    )

    # Pipeline-specific options
    optp.add_argument('-o', '--outdir', metavar='STR', dest='outdir', required=True,
                      help='The output directory for VCF files and logs.')
    optp.add_argument('--ref_fai', metavar='FILE', dest='ref_fai', required=True,
                      help='The reference fai file (required). Used to determine chromosome '
                           'lengths and to split regions into sub-jobs.')
    optp.add_argument('-d', '--delta', metavar='INT', dest='delta', type=int,
                      help='Region size for each sub-job [2000000].', default=2000000)
    optp.add_argument('-c', '--chrom', metavar='STR', dest='chrom',
                      help='Only process comma-delimited chromosomes, e.g. chr1,chr2.',
                      default='')

    # Parse known args: pipeline-specific args are consumed here,
    # the rest are passed through to basevar caller directly.
    opt, basevar_passthrough_args = optp.parse_known_args()

    opt.outdir = os.path.abspath(opt.outdir)
    chroms_filter = [c.strip() for c in opt.chrom.split(',') if c.strip()] if opt.chrom else []

    # Build ref_fai_map: {chrom: length} for region validation and coordinate lookup
    ref_fai_entries = load_reference_fai(opt.ref_fai, chroms_filter if chroms_filter else None)
    ref_fai_map = {entry[0]: entry[2] for entry in ref_fai_entries}

    if not ref_fai_map:
        print("Error: No chromosomes loaded from --ref_fai "
              "(check --chrom filter or fai file).", file=sys.stderr)
        sys.exit(1)

    # Intercept -r/--regions from passthrough args.
    # The pipeline needs to split user-specified regions into sub-regions,
    # so -r must be extracted, parsed, then re-injected per sub-job.
    regions_str = ''
    i = 0
    while i < len(basevar_passthrough_args):
        arg = basevar_passthrough_args[i]
        if arg in ('-r', '--regions') and i + 1 < len(basevar_passthrough_args):
            regions_str = basevar_passthrough_args[i + 1]
            basevar_passthrough_args.pop(i)
            basevar_passthrough_args.pop(i)
            continue
        elif arg.startswith('--regions='):
            regions_str = arg.split('=', 1)[1]
            basevar_passthrough_args.pop(i)
            continue
        i += 1

    # Determine the working region list
    if regions_str:
        try:
            region_list = parse_regions(regions_str, ref_fai_map, chroms_filter)
        except ValueError as e:
            print(f"Error parsing -r/--regions: {e}", file=sys.stderr)
            sys.exit(1)

        if not region_list:
            print("Error: -r/--regions produced no valid regions "
                  "(check --chrom filter or region strings).", file=sys.stderr)
            sys.exit(1)
    else:
        # Use all chromosomes from fai (already filtered by --chrom if set)
        region_list = [(entry[0], entry[1], entry[2]) for entry in ref_fai_entries]

    # Build the base passthrough string (everything except -r, which we inject per sub-job)
    passthrough_str = ' '.join(basevar_passthrough_args)

    for chr_id, reg_start, reg_end in region_list:
        for i in range(reg_start - 1, reg_end, opt.delta):
            start = i + 1
            end = min(i + opt.delta, reg_end)
            reg = f'{chr_id}:{start}-{end}'

            outfile_prefix = f'{chr_id}_{start}_{end}'
            out_vcf = os.path.join(opt.outdir, f'{outfile_prefix}.vcf.gz')
            out_log = os.path.join(opt.outdir, f'{outfile_prefix}.log')

            cmd = (f'time {exe_prog} {passthrough_str} '
                   f'-r {reg} '
                   f'-o {out_vcf} '
                   f'> {out_log} && '
                   f'echo "** {outfile_prefix} done **"')
            print(cmd)


if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print("Please type: -h or --help to show the help message.\n", file=sys.stderr)
        sys.exit(1)

    creat_basetype_pipe()