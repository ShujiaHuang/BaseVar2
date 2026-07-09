# Building BaseVar from Source

This document describes how to compile BaseVar from source code. Most users should use the [pre-built binary](https://github.com/ShujiaHuang/BaseVar2/releases) instead.

*Requires: C++17 compiler (GCC 7+ or Apple Clang 10+), CMake ≥ 3.12, and system libraries: zlib, bzip2, xz-utils, libcurl.*

---

## Option 1 — Build with CMake (standard dynamic build)

### Step 1 — Clone the repository (including htslib submodule)

```bash
git clone --recursive https://github.com/ShujiaHuang/basevar2.git
cd basevar2
```

> If you forgot `--recursive`, run: `git submodule update --init --recursive`

### Step 2 — Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

The executable `bin/basevar` will be produced. Verify with:

```bash
./bin/basevar --help
```

### Step 3 (Optional) — Build a static binary locally

**macOS** (requires Homebrew):

```bash
brew install zlib bzip2 xz
cmake -B build-static -DSTATIC_BUILD=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-static
```

**Linux** (portable static via Ubuntu/glibc):

```bash
sudo apt-get install -y build-essential cmake autoconf automake \
    zlib1g-dev libbz2-dev liblzma-dev libssl-dev
cmake -B build-static -DSTATIC_BUILD=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-static
```

This bundles `libstdc++`, `libgcc`, `htslib`, and the compression libs statically; glibc remains dynamic. The resulting binary runs on the build host and on any other host with the same-or-newer glibc.

---

## Option 2 — Manual g++ compilation (fallback)

First, build htslib:

```bash
cd htslib && autoreconf -i && ./configure && make && cd ..
```

Then compile manually:

**Linux:**

```bash
cd bin/
g++ -O3 -fPIC ../src/*.cpp ../src/io/*.cpp ../htslib/libhts.a \
    -I ../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -lssl -lcrypto -o basevar
```

**macOS:**

```bash
cd bin/
g++ -O3 -fPIC ../src/*.cpp ../src/io/*.cpp ../htslib/libhts.a \
    -I ../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -o basevar
```

> [!NOTE]
> **Note:** If you encounter a `test/test_khash.c` compilation error during `make` in htslib, you can safely ignore it — the required `libhts.a` archive is still produced correctly.
