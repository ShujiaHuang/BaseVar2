# BaseVar v2.6.5 Release Notes

## Performance Improvements

### Caller: BGZF block-level merge for multi-thread and multi-region output

The `basevar caller` internal VCF merging has been optimized with BGZF block-level copy, eliminating decompression/recompression overhead:

- **Thread-level merge** (previously `merge_file_by_line`): per-thread temporary `.partN.vcf.gz` files are now merged via raw BGZF block copying (`bgzf_raw_concat`). Since these temporary files contain no VCF header (headerless VCF data streams), the header is written once up-front and the remaining blocks are copied directly.

- **Multi-region merge** (`.vcf.gz` output): when the output format is BGZF-compressed, multiple region VCFs are now merged via header-aware block-level copy (`bgzf_naive_concat`), which writes the header from the first file and skips headers from subsequent files. When the output is plain-text `.vcf`, the original line-by-line merge is preserved.

### Shared BGZF block-level copy engine

The core block-level copy functions (`skip_vcf_gz_header`, `bgzf_copy_blocks_to`, `bgzf_raw_concat`, `bgzf_naive_concat`) have been extracted into `src/io/bgzf_concat.h` as a shared, header-only library. Both `concat.cpp` and `variant_caller.cpp` now use this common engine, eliminating code duplication.

## Enhancements

### `concat --naive` mode: input/output format validation

- **BCF detection**: input files starting with non-`#` bytes (BCF binary format) are rejected early with a clear error message.
- **BGZF format check**: non-BGZF-compressed inputs are detected and rejected.
- **Output suffix validation**: output file must have a `.gz` suffix since `--naive` mode always produces BGZF-compressed output.

### `basetype` subcommand: exception handling

Added try-catch wrapper around `BaseTypeRunner` to catch unhandled exceptions and print a user-friendly error message instead of an abort.

## Documentation

- `README.md`: updated `basevar concat` usage section with `--naive` mode documentation and revised examples.
- `README.md`: corrected default thread count display from `[14]` to `[auto]`.

## Tests

- `tests/test_block_merge.py`: comprehensive block-level merge validation (28 test cases covering multi-region, single-region, thread count consistency, empty regions, and mixed empty/non-empty regions).
- `tests/io/test_bgzf_utils.cpp`: extended BGZF unit tests (85 tests including `bgzf_copy_blocks_to` and `bgzf_raw_concat`).

## Files Changed

- `src/variant_caller.cpp` — Block-level merge at L913 (thread-level) and L482 (multi-region)
- `src/variant_caller.h` — Minor comment update
- `src/io/bgzf_concat.h` — New shared BGZF block-level copy engine (header-only)
- `src/concat.cpp` — Eliminated duplicated merge loop; delegates to `bgzf_naive_concat`
- `src/io/iobgzf.h` — Extended BGZFile with `raw_block_read`, `raw_block_write`, `flush`, `write_raw`, block inspection APIs
- `src/main.cpp` — try-catch wrapper for `basetype()`
- `README.md` — concat usage docs, thread default fix
- `tests/test_block_merge.py` — New block-level merge comparison test
- `tests/io/test_bgzf_utils.cpp` — Extended BGZF unit tests
- `CMakeLists.txt` — Version bump to 2.6.5
