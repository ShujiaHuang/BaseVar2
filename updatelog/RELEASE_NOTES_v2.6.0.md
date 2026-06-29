# BaseVar v2.6.0 Release Notes

## What's New

### New `basevar dump` subcommand

A new `dump` subcommand for inspecting intermediate binary batchfile (`.bbf`) and binary index (`.bbi`) files. This is useful for debugging, verifying file integrity, or examining per-sample read data at specific genomic positions.

```bash
# Inspect .bbi index summary (magic, version, entries, position range, footer integrity)
basevar dump sample.bbf.bbi

# List all index entries
basevar dump sample.bbf.bbi --entries

# Show .bbf header (sample IDs) only
basevar dump sample.bbf --header

# Show all positions in compact summary format
basevar dump sample.bbf

# Show per-sample details at a specific position
basevar dump sample.bbf -r chr11:5246595-5246595 -v

# Dump first 20 positions with full per-sample read data
basevar dump sample.bbf -n 20 -v
```

## Bug Fixes

### Critical: Fix VCF coordinate corruption in binary batchfile reader

Fixed a severe bug where `empty_bi_template` (used for samples with no data at a given position) had default `ref_id=""` and `ref_pos=0`. When the first reader in a multi-batchfile run had no data at a position, this propagated through to `_basetype_caller_unit` → `VariantInfo` → VCF output, resulting in empty chromosome names (`chrom=""`) and zero positions (`pos=0`) in the VCF.

The fix adds a backfill step after the reader loop to ensure all samples have consistent `ref_id`/`ref_pos` values before calling the variant caller.

### Fix redundant file open in binary batchfile reader

Eliminated a redundant `BGZFile::open()` call that was used solely to read `sample_count` from the header. The header is now read once from the already-open file handle, avoiding 100+ unnecessary `open/close` cycles per region.

## Performance Improvements

### Binary batchfile read path optimizations

- **Vector pre-allocation**: `all_smps_bi` is now declared outside the position loop with `reserve()`, avoiding per-position memory reallocation.
- **Pre-constructed empty template**: `empty_bi_template` is built once and copied for empty samples, replacing repeated construction of default `BatchInfo` objects.
- **Deferred ref_id assignment**: `ref_id` string is now set once per reader (after `read_binary_record()`) instead of being copied for every individual read, reducing string copies from ~1000× per position to ~100×.

### HANDOFF documentation cleanup

- Removed 82 lines of outdated bottleneck analysis (all items already implemented).
- Corrected `.bbf` and `.bbi` format descriptions to match actual implementation.
- Added `batchfile_binary.h/.cpp` to file index.

## Full Changelog

**Modified files:**
- `src/variant_caller.cpp` — binary reader optimizations + empty_bi_template ref_id/ref_pos bug fix
- `src/io/batchfile_binary.cpp` — deferred ref_id assignment optimization
- `src/main.cpp` — register `dump` subcommand
- `CMakeLists.txt` — version bump to 2.6.0
- `README.md` — document `dump` subcommand, binary format, BBI footer validation

**New files:**
- `src/dump.h` — dump subcommand header
- `src/dump.cpp` — dump subcommand implementation
