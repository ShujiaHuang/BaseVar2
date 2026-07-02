#!/usr/bin/env python3
"""
Comprehensive test: verify block-level merge produces identical output to line-by-line merge.

Tests:
  1. Multi-region .vcf.gz (block-level merge at L482) vs multi-region .vcf (line-by-line at L489)
  2. Single-region .vcf.gz (rename path at L481)
  3. Multi-region with different thread counts (1 vs 4 vs 14) to exercise L913 thread-level merge
  4. Empty region (no variants) — header-only output
"""

import subprocess
import os
import sys
import tempfile
import hashlib
import shutil

BASEVAR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin", "basevar")
REF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tests", "data", "synthetic", "ref", "mini_ref.fa")
BAMS = [
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tests", "data", "synthetic", "bam", f"sampleA0{i}.bam")
    for i in range(1, 4)
]

def run_caller(output, regions=None, threads=None):
    """Run basevar caller and return (returncode, stdout+stderr)."""
    cmd = [BASEVAR, "caller", "-f", REF, "-o", output]
    if regions:
        cmd += ["-r", regions]
    if threads:
        cmd += ["-t", str(threads)]
    cmd += BAMS
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    return r.returncode, r.stdout + r.stderr

def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()

def decompress_vcf_gz(path):
    """Decompress .vcf.gz and return content as string."""
    r = subprocess.run(["bgzip", "-d", "-c", path], capture_output=True, text=True, timeout=30)
    return r.stdout

def strip_cmd_header(text):
    """Remove ##basevar_caller_command= line (varies by output filename/thread count)."""
    lines = text.split("\n")
    return "\n".join(l for l in lines if not l.startswith("##basevar_caller_command="))

def cleanup(path):
    """Remove file and its .tbi index if present."""
    for p in [path, path + ".tbi"]:
        if os.path.exists(p):
            os.remove(p)
    # Remove cache dir if present
    cache_dir = os.path.join(os.path.dirname(path), "cache_" + os.path.splitext(os.path.basename(path))[0])
    if os.path.exists(cache_dir):
        shutil.rmtree(cache_dir)

passed = 0
failed = 0

def check(name, condition, detail=""):
    global passed, failed
    if condition:
        passed += 1
        print(f"  [PASS] {name}")
    else:
        failed += 1
        print(f"  [FAIL] {name}  {detail}")

tmpdir = tempfile.mkdtemp(prefix="basevar_test_")
print(f"Test directory: {tmpdir}")

# ============================================================
# Test 1: Multi-region .vcf.gz (block-level) vs .vcf (line-by-line)
# ============================================================
print("\n=== Test 1: Multi-region block-level vs line-by-line ===")

vcf_gz = os.path.join(tmpdir, "multi_region.vcf.gz")
vcf_plain = os.path.join(tmpdir, "multi_region.vcf")

rc1, _ = run_caller(vcf_gz, regions="chr1:1-500,chr1:501-2000")
rc2, _ = run_caller(vcf_plain, regions="chr1:1-500,chr1:501-2000")

check("vcf.gz caller exit code", rc1 == 0, f"rc={rc1}")
check("vcf caller exit code", rc2 == 0, f"rc={rc2}")

if rc1 == 0 and rc2 == 0:
    content_gz = decompress_vcf_gz(vcf_gz)
    with open(vcf_plain) as f:
        content_plain = f.read()
    check("Content identical (block-level == line-by-line)",
          strip_cmd_header(content_gz) == strip_cmd_header(content_plain),
          f"len(gz)={len(content_gz)} vs len(plain)={len(content_plain)}")

    # Check header present exactly once
    header_lines = [l for l in content_gz.split("\n") if l.startswith("#")]
    chrom_count = sum(1 for l in header_lines if l.startswith("#CHROM"))
    check("Exactly one #CHROM line", chrom_count == 1, f"found {chrom_count}")

    # Check data lines
    data_lines = [l for l in content_gz.split("\n") if l and not l.startswith("#")]
    check("Has variant data", len(data_lines) > 0, f"found {len(data_lines)} variants")

    # BGZF validity
    r = subprocess.run(["bgzip", "-t", vcf_gz], capture_output=True, text=True)
    check("BGZF valid", r.returncode == 0, r.stderr)

    # Tabix index
    check("Tabix index exists", os.path.exists(vcf_gz + ".tbi"))
    r = subprocess.run(["tabix", vcf_gz, "chr1:1-500"], capture_output=True, text=True)
    check("Tabix query chr1:1-500", r.returncode == 0 and len(r.stdout) > 0)

cleanup(vcf_gz)
cleanup(vcf_plain)

# ============================================================
# Test 2: Single-region .vcf.gz (rename path)
# ============================================================
print("\n=== Test 2: Single-region (rename path) ===")

vcf_single = os.path.join(tmpdir, "single_region.vcf.gz")
rc, _ = run_caller(vcf_single, regions="chr1:1-2000")
check("Single-region exit code", rc == 0, f"rc={rc}")

if rc == 0:
    content = decompress_vcf_gz(vcf_single)
    data_lines = [l for l in content.split("\n") if l and not l.startswith("#")]
    check("Single-region has variants", len(data_lines) > 0)

    r = subprocess.run(["bgzip", "-t", vcf_single], capture_output=True, text=True)
    check("Single-region BGZF valid", r.returncode == 0, r.stderr)

cleanup(vcf_single)

# ============================================================
# Test 3: Different thread counts (exercises L913 thread-level merge)
# ============================================================
print("\n=== Test 3: Thread count consistency (1 vs 4 vs 14) ===")

results = {}
for t in [1, 4, 14]:
    out = os.path.join(tmpdir, f"thread{t}.vcf.gz")
    rc, _ = run_caller(out, regions="chr1:1-2000", threads=t)
    check(f"Thread-{t} exit code", rc == 0, f"rc={rc}")
    if rc == 0:
        content = decompress_vcf_gz(out)
        data = sorted([l for l in content.split("\n") if l and not l.startswith("#")])
        header = [l for l in content.split("\n") if l.startswith("#")]
        results[t] = {"data": data, "header": header, "content": content}

        r = subprocess.run(["bgzip", "-t", out], capture_output=True, text=True)
        check(f"Thread-{t} BGZF valid", r.returncode == 0, r.stderr)

if len(results) == 3:
    check("Thread-1 == Thread-4 data",
          results[1]["data"] == results[4]["data"],
          f"t1={len(results[1]['data'])} vs t4={len(results[4]['data'])}")
    check("Thread-1 == Thread-14 data",
          results[1]["data"] == results[14]["data"],
          f"t1={len(results[1]['data'])} vs t14={len(results[14]['data'])}")
    check("Thread-1 == Thread-4 header (excl. cmdline)",
          strip_cmd_header("\n".join(results[1]["header"])) == strip_cmd_header("\n".join(results[4]["header"])))

for t in [1, 4, 14]:
    cleanup(os.path.join(tmpdir, f"thread{t}.vcf.gz"))

# ============================================================
# Test 4: Empty region (header-only output)
# ============================================================
print("\n=== Test 4: Empty region (header-only) ===")

vcf_empty = os.path.join(tmpdir, "empty_region.vcf.gz")
# chr2:1-100 should have no variants in the synthetic data
rc, output = run_caller(vcf_empty, regions="chr2:1-100")
check("Empty region exit code", rc == 0, f"rc={rc}")

if rc == 0:
    content = decompress_vcf_gz(vcf_empty)
    data_lines = [l for l in content.split("\n") if l and not l.startswith("#")]
    check("Empty region has no data", len(data_lines) == 0, f"found {len(data_lines)} variants")
    header_lines = [l for l in content.split("\n") if l.startswith("#")]
    check("Empty region has header", len(header_lines) > 0)

    r = subprocess.run(["bgzip", "-t", vcf_empty], capture_output=True, text=True)
    check("Empty region BGZF valid", r.returncode == 0, r.stderr)

cleanup(vcf_empty)

# ============================================================
# Test 5: Multi-region with empty + non-empty mix
# ============================================================
print("\n=== Test 5: Mixed empty + non-empty regions ===")

vcf_mix = os.path.join(tmpdir, "mixed_region.vcf.gz")
rc, _ = run_caller(vcf_mix, regions="chr2:1-100,chr1:1-2000")
check("Mixed region exit code", rc == 0, f"rc={rc}")

if rc == 0:
    content = decompress_vcf_gz(vcf_mix)
    data_lines = [l for l in content.split("\n") if l and not l.startswith("#")]
    check("Mixed region has variants", len(data_lines) > 0)

    r = subprocess.run(["bgzip", "-t", vcf_mix], capture_output=True, text=True)
    check("Mixed region BGZF valid", r.returncode == 0, r.stderr)

    # Compare with single-region result
    if 1 in results:
        check("Mixed region data == single-region data",
              sorted(data_lines) == results[1]["data"],
              f"mix={len(data_lines)} vs single={len(results[1]['data'])}")

cleanup(vcf_mix)

# ============================================================
# Summary
# ============================================================
print(f"\n{'='*50}")
print(f"Results: {passed} passed, {failed} failed")
print(f"{'='*50}")

shutil.rmtree(tmpdir, ignore_errors=True)
sys.exit(0 if failed == 0 else 1)
