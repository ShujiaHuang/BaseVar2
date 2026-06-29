# BaseVar v2.5.3 Release Notes

## Overview

v2.5.3 is a performance optimization release that significantly reduces runtime overhead in the variant calling pipeline through four targeted improvements: eliminating `ostringstream` conversion bottlenecks, replacing per-base string allocations with `char`, adding move semantics for `BamRecord`, and introducing EM hot-start for the step-down procedure.

## Performance Optimizations

### 1. P0-1: `join()` Fast-Path Specializations

- **Problem**: The generic `join<T>()` template used `ostringstream` for every element conversion (~300ns per element), even for trivially convertible types like `char` and `int`.
- **Fix**: Added three `inline` overloads for `vector<string>`, `vector<char>`, and `vector<int>` that bypass `ostringstream` entirely. The `string` overload pre-calculates total size for `reserve()`, the `char` overload uses `push_back()`, and the `int` overload uses `std::to_string()`.
- **Impact**: Eliminates the dominant overhead in batchfile text serialization, which processes billions of data elements per genome.

### 2. P1-1: `ReadAlignedPair` String-to-Char Optimization

- **Problem**: Each `ReadAlignedPair` stored `ref_base`, `read_base`, and `read_qual` as `std::string`, causing 3 heap allocations per aligned base (~150ns × 3 per read position). For a 30x genome, this amounts to ~90 billion allocations.
- **Fix**: Changed `ref_base`, `read_base`, and `read_qual` from `std::string` to `char`. Added a `multi_base` string field that is only allocated for multi-base Indels and soft-clips (rare events).
- **Impact**: Eliminates ~90 billion string allocations per genome. The `multi_base` field is empty (no heap allocation) for the vast majority of entries (M/=/X operations).

### 3. P1-2: `BamRecord` Move Semantics

- **Problem**: `push_back(al)` into `vector<BamRecord>` triggered the copy constructor, which called `bam_dup1()` to deep-copy the entire `bam1_t` structure (~250ns per read).
- **Fix**: Added move constructor and move assignment operator to `BamRecord`. Changed `push_back(al)` to `push_back(std::move(al))`. After move, `bf.next(al)` correctly reinitializes via `load_read()`.
- **Impact**: Eliminates one `bam_dup1()` deep copy per read (~250ns × 10⁹ reads ≈ 4 minutes per genome).

### 4. P1-4: EM Algorithm Hot-Start

- **Problem**: In the step-down procedure, each level's EM started from observed base frequencies, ignoring the converged frequencies already computed at the previous (higher) level.
- **Fix**: After each step-down level converges, its frequency vector is passed as the initial value for the next level. The hot-start frequency is projected onto the current combination by zeroing out alleles not in the combination and renormalizing.
- **Impact**: Reduces EM iterations at each step-down level, with cumulative savings across all variant sites genome-wide.

## Changed Files

| File | Change |
|------|--------|
| `src/io/utils.h` | Added `join()` overloads for `vector<string>`, `vector<char>`, `vector<int>` |
| `src/io/bam_record.h` | `ReadAlignedPair` fields `string` → `char`; added `multi_base`; declared move constructor/assignment |
| `src/io/bam_record.cpp` | `get_aligned_pairs()` uses char indexing; implemented move constructor/assignment |
| `src/basetype.h` | Added `hot_start_freq` parameter to `_f()` |
| `src/basetype.cpp` | Implemented hot-start frequency projection in `_f()`; pass converged frequencies through step-down levels in `lrt()` |
| `src/variant_caller.cpp` | Updated consumer code for `ReadAlignedPair` char fields; `push_back(std::move(al))` |
| `CMakeLists.txt` | Version 2.5.2 → 2.5.3 |

## Verification

All optimizations have been verified to produce identical VCF output (byte-for-byte match via `diff`) compared to v2.5.2, in both single-threaded and multi-threaded execution modes.

## Upgrade Recommendation

All v2.5.2 users are recommended to upgrade to v2.5.3 for improved runtime performance with no change in output.

## Downloads

- **Linux static binary**: `basevar-v2.5.3-linux-amd64`
- **macOS static binary**: `basevar-v2.5.3-macos-amd64`
- **Source code**: Source code (tar.gz / zip)
