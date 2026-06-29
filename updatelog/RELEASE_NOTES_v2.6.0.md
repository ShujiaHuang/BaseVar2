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

### Binary batchfile format (Plan C)

BaseVar now uses a BGZF-compressed binary batchfile format (`.bbf` + `.bbi` index) instead of the text-based format. This provides significant performance improvements:

- **~1.4× speedup** over text batchfiles by eliminating join/split overhead
- **Reduced BGZF I/O calls** via buffer-serialized binary records
- **Compact per-sample data**: 9-byte `PackedReadFields` struct read in a single `bgzf_read()` call
- **Natural data boundaries**: `uint16_t depth` per sample prevents cross-sample contamination

### Sparse binary index (`.bbi`)

The binary index skips positions where all samples have zero depth, reducing index size from ~GB to ~MB for typical cohorts. This enables:

- Faster random access during variant calling
- Reduced memory footprint (BBI indexes are ~1000× smaller than full indexes)
- Efficient multi-region queries

### Shared index loading

When using multiple threads (`-t > 1`), BBI indexes are loaded once in `_variants_discovery` and shared across all threads via `std::cref`. This eliminates redundant I/O and significantly reduces startup time for large cohorts.

### BBI footer integrity marker

Binary index files now include a footer integrity marker (`BBI_FOOTER_MAGIC = 0x42424946`) that is validated by `--smart-rerun` to detect truncated or corrupted index files before resuming an interrupted run.

## Bug Fixes

### Critical: Fix sparse index lockstep position loss

Fixed a severe bug where the "seek all readers to the same starting position" logic would skip variant positions when different batchfiles had different first entries (due to sparse indexes). The fix removes the forced max-position seek and lets each reader independently seek to its own `lower_bound` position. The `min_pos` lockstep loop naturally handles different starting positions.

**Symptom**: With `-t 4`, specific variant positions (e.g., chr11:5247141) were missing from the output. With `-t 1` on a sub-region, the same positions were also missing.

**Root cause**: `current_pos = max(all readers' first entry positions)` advanced ALL readers to the maximum start position, skipping positions in between.

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

### Buffer-serialized binary I/O

- **Write path**: Records are serialized into a memory buffer first, then written in a single `bgzf_write()` call, eliminating thousands of small BGZF calls per position.
- **Read path**: Fixed fields are packed into a 9-byte `PackedReadFields` struct and read in one `bgzf_read()`, then unpacked in memory.

### Eliminate redundant `seek_virtual`

Removed redundant `seek_virtual()` calls in the read path. The file position is now tracked naturally after each `read_binary_record()`, avoiding unnecessary seeks.

## Code Cleanup

### Dead code removal

Removed 3 obsolete functions (~260 lines) that were superseded by the binary batchfile implementation:

- `_write_record_to_batchfile()` — old text-format batchfile writer
- `_get_sampleid_from_batchfiles()` — old text-format sample ID reader
- `_basevar_caller()` — old text-format variant caller

## Testing

### Comprehensive unit tests for binary batchfile format

Added 91 unit tests covering:

- Header roundtrip (write + read)
- SNP/insertion/deletion roundtrip
- Empty positions
- Long-read RPR (>255 bases)
- Long indels (>255 bp)
- Strand '*' handling
- Multi-position sequential read
- Index binary search
- Multi-sample varying depths
- Mixed variant types
- BBI footer integrity validation

## Full Changelog

**Modified files:**
- `src/variant_caller.cpp` — binary reader optimizations + sparse index bug fix + dead code removal
- `src/variant_caller.h` — dead code declaration removal
- `src/io/batchfile_binary.cpp` — buffer-serialized I/O + deferred ref_id assignment
- `src/io/batchfile_binary.h` — `PackedReadFields` struct + format constants
- `src/io/iobgzf.h` — `seek_virtual()` + `tell_virtual()` enhancements
- `src/main.cpp` — register `dump` subcommand
- `CMakeLists.txt` — version bump to 2.6.0
- `README.md` — document `dump` subcommand, sparse index, shared index loading

**New files:**
- `src/dump.cpp` — dump subcommand implementation
- `src/dump.h` — dump subcommand header
- `src/io/batchfile_binary.cpp` — binary batchfile read/write implementation
- `src/io/batchfile_binary.h` — binary batchfile format definitions
- `tests/io/test_batchfile_binary.cpp` — 91 unit tests for binary format
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
