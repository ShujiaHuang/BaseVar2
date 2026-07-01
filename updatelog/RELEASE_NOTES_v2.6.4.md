# BaseVar v2.6.4 Release Notes

## Bug Fixes

### Fix `getopt_long` `optional_argument` parsing issue in `caller` and `concat`

All options that accept a value in `basevar caller` and `basevar concat` have been changed from `optional_argument` to `required_argument` in the `getopt_long` definition.

**Affected options (`caller`):** `-L`, `-m`, `-q`, `-Q`, `-B`, `-t`, `-r`, `-G`
**Affected options (`concat`):** `-L`

**Root cause:** With `optional_argument`, `getopt_long` requires the value to be attached via `=` (e.g., `--thread=64`). Using space-separated syntax (e.g., `--thread 64`) caused the value to be interpreted as a positional argument, leading to errors such as:

```
terminate called after throwing an instance of 'std::runtime_error'
  what():  [bam.cpp::Bam:_open] file not found - 64
```

**After fix:** Both syntaxes now work correctly:

```bash
basevar caller --thread 64 -f ref.fa ...    # space-separated ✓
basevar caller --thread=64 -f ref.fa ...    # equals-sign ✓
```

All options remain **optional** — `required_argument` only means the value must be provided *when the option is used*, not that the option itself is mandatory.

## Files Changed

- `src/variant_caller.cpp` — Changed 8 `optional_argument` → `required_argument`; updated comment
- `src/concat.cpp` — Changed 1 `optional_argument` → `required_argument`
- `CMakeLists.txt` — Version bump to 2.6.4
