#!/usr/bin/env bash
# =============================================================================
# run_tests.sh — Systematic test suite for basevar caller
#
# Uses synthetic test data (25 samples, 3 groups, 11 variants) to validate
# all major caller parameters and workflows.
#
# Usage:
#   bash run_tests.sh              # Run all tests
#   bash run_tests.sh --quick      # Skip evaluation, just check exit codes
#   bash run_tests.sh --evaluate   # Run evaluation on key multi-sample tests
#
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASEVAR="${BASEVAR:-$SCRIPT_DIR/../../../bin/basevar}"
REF="$SCRIPT_DIR/ref/mini_ref.fa"
BAM_DIR="$SCRIPT_DIR/bam"
CRAM_DIR="$SCRIPT_DIR/cram"
TRUTH="$SCRIPT_DIR/ground_truth_variants.tsv"
SAMPLES_LIST="$SCRIPT_DIR/samples.list"
SAMPLES_GROUP="$SCRIPT_DIR/samples_group.info"
EVALUATE="$SCRIPT_DIR/evaluate.py"

TMP_DIR=$(mktemp -d "$SCRIPT_DIR/.test_output.XXXXXX")
trap 'rm -rf "$TMP_DIR"' EXIT

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
RUN_EVALUATION=false
for arg in "$@"; do
    case "$arg" in
        --quick)     RUN_EVALUATION=false ;;
        --evaluate)  RUN_EVALUATION=true ;;
        --help|-h)
            echo "Usage: bash run_tests.sh [--quick|--evaluate]"
            echo "  --quick     Skip evaluation (default)"
            echo "  --evaluate  Run evaluate.py on multi-sample tests"
            exit 0
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Counters
# ---------------------------------------------------------------------------
PASS=0
FAIL=0
SKIP=0
TOTAL=0

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
run_test() {
    local test_id="$1"
    local description="$2"
    shift 2
    local cmd=("$@")

    TOTAL=$((TOTAL + 1))
    printf "  [T%02d] %-50s " "$TOTAL" "$description"

    local log="$TMP_DIR/${test_id}.log"
    if "${cmd[@]}" > "$log" 2>&1; then
        echo "PASS"
        PASS=$((PASS + 1))
        return 0
    else
        echo "FAIL"
        FAIL=$((FAIL + 1))
        echo "       Command: ${cmd[*]}"
        echo "       Log: $log"
        # Show last few lines of stderr for debugging
        tail -3 "$log" 2>/dev/null | sed 's/^/       /'
        return 1
    fi
}

run_eval() {
    local vcf_file="$1"
    local label="$2"

    if [ ! -f "$vcf_file" ]; then
        echo "       [eval] VCF not found: $vcf_file"
        return 0  # Don't fail the script for missing VCF
    fi

    local eval_out
    eval_out=$(python3 "$EVALUATE" --vcf "$vcf_file" --truth "$TRUTH" 2>&1) || true

    local sens prec gt_conc
    sens=$(echo "$eval_out" | grep "Sensitivity" | awk '{print $NF}' || echo "N/A")
    prec=$(echo "$eval_out" | grep "Precision" | awk '{print $NF}' || echo "N/A")
    gt_conc=$(echo "$eval_out" | grep "Concordant:" | awk '{print $NF}' || echo "N/A")

    printf "       └─ %-8s  Sens=%-8s  Prec=%-8s  GT_conc=%s\n" \
        "$label" "$sens" "$prec" "$gt_conc"
}

collect_bams() {
    # Resolve relative BAM paths from samples.list against the script directory
    while IFS= read -r line; do
        [ -z "$line" ] && continue
        if [[ "$line" = /* ]]; then
            echo "$line"
        else
            echo "$SCRIPT_DIR/$line"
        fi
    done < "$SAMPLES_LIST" | tr '\n' ' '
}

# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------
echo "============================================================"
echo "  BaseVar Caller Test Suite"
echo "============================================================"
echo "  basevar:    $BASEVAR"
echo "  reference:  $REF"
echo "  output:     $TMP_DIR"
echo "  evaluate:   $RUN_EVALUATION"
echo ""

if [ ! -x "$BASEVAR" ]; then
    echo "  ERROR: basevar not found or not executable: $BASEVAR"
    echo "  Set BASEVAR env variable to the correct path."
    exit 1
fi

if [ ! -f "$REF" ]; then
    echo "  ERROR: Reference not found: $REF"
    echo "  Run 'python3 generate.py' first."
    exit 1
fi

echo "  Pre-flight: basevar version:"
"$BASEVAR" --version 2>/dev/null || "$BASEVAR" 2>&1 | head -1 | sed 's/^/    /'
echo ""

# Collect BAM paths
ALL_BAMS=$(collect_bams)
FIRST_BAM=$(echo "$ALL_BAMS" | awk '{print $1}')

# =============================================================================
echo "============================================================"
echo "  Category 1: Input Modes"
echo "================================================================"
echo ""

# T01: Single sample BAM (positional)
run_test "T01" "Single sample BAM (positional)" \
    "$BASEVAR" caller -f "$REF" -o "$TMP_DIR/T01.vcf" "$FIRST_BAM"

# T02: Multi-sample via -L samples.list
run_test "T02" "Multi-sample via -L samples.list" \
    "$BASEVAR" caller -f "$REF" -o "$TMP_DIR/T02.vcf" \
    -L "$SAMPLES_LIST" --filename-has-samplename

# T03: Multi-sample via -L + positional BAMs
run_test "T03" "Multi-sample via -L + positional BAMs" \
    "$BASEVAR" caller -f "$REF" -o "$TMP_DIR/T03.vcf" \
    -L "$SAMPLES_LIST" --filename-has-samplename $ALL_BAMS

# T08: Paired-end generation (writes paired BAMs)
run_test "T08" "Paired-end generation -- paired-end BAMs" \
    python3 "$SCRIPT_DIR/generate.py" --paired-end --skip-cram --output-dir "$TMP_DIR/paired" --seed 42

# T04: CRAM input (single sample)
FIRST_CRAM=$(ls "$CRAM_DIR"/*.cram | head -1)
run_test "T04" "CRAM input (single sample)" \
    "$BASEVAR" caller -f "$REF" -o "$TMP_DIR/T04.vcf" "$FIRST_CRAM"

# T05: CRAM input (multi-sample via -L)
# Build a CRAM list
CRAM_LIST="$TMP_DIR/cram.list"
ls "$CRAM_DIR"/*.cram > "$CRAM_LIST"
run_test "T05" "CRAM multi-sample via -L" \
    "$BASEVAR" caller -f "$REF" -o "$TMP_DIR/T05.vcf" \
    -L "$CRAM_LIST" --filename-has-samplename

# =============================================================================
echo ""
echo "================================================================"
echo "  Category 2: Region Calling (-r)"
echo "================================================================"
echo ""

# T06: Single region
run_test "T06" "Single region -r chr1:100-1100" \
    "$BASEVAR" caller -f "$REF" -o "$TMP_DIR/T06.vcf" \
    -L "$SAMPLES_LIST" --filename-has-samplename -r chr1:100-1100

# T07: Multiple regions (comma-separated)
run_test "T07" "Multiple regions -r chr1:100-600,chr2:200-800" \
    "$BASEVAR" caller -f "$REF" -o "$TMP_DIR/T07.vcf" \
    -L "$SAMPLES_LIST" --filename-has-samplename \
    -r chr1:100-600,chr2:200-800

# T08: Whole chromosome
run_test "T08" "Whole chromosome -r chr1" \
    "$BASEVAR" caller -f "$REF" -o "$TMP_DIR/T08.vcf" \
    -L "$SAMPLES_LIST" --filename-has-samplename -r chr1

# =============================================================================
echo ""
echo "================================================================"
echo "  Category 3: Population Stratification (-G)"
echo "================================================================"
echo ""

# T09: With -G group file
run_test "T09" "With -G group stratification" \
    "$BASEVAR" caller -f "$REF" -o "$TMP_DIR/T09.vcf" \
    -L "$SAMPLES_LIST" --filename-has-samplename -G "$SAMPLES_GROUP"

# T10: -G + -r combined
run_test "T10" "-G + -r combined" \
    "$BASEVAR" caller -f "$REF" -o "$TMP_DIR/T10.vcf" \
    -L "$SAMPLES_LIST" --filename-has-samplename \
    -G "$SAMPLES_GROUP" -r chr1:100-1100

# =============================================================================
echo ""
echo "================================================================"
echo "  Category 4: Genotype Mode (--gt-mode)"
echo "================================================================"
echo ""

# T11: Posterior mode (default)
run_test "T11" "--gt-mode posterior (default)" \
    "$BASEVAR" caller -f "$REF" -o "$TMP_DIR/T11.vcf" \
    -L "$SAMPLES_LIST" --filename-has-samplename --gt-mode posterior

# T12: Legacy mode
run_test "T12" "--gt-mode legacy" \
    "$BASEVAR" caller -f "$REF" -o "$TMP_DIR/T12.vcf" \
    -L "$SAMPLES_LIST" --filename-has-samplename --gt-mode legacy

# =============================================================================
echo ""
echo "================================================================"
echo "  Category 5: Smart Rerun (--smart-rerun)"
echo "================================================================"
echo ""

# T13: First run (creates batchfiles)
run_test "T13" "--smart-rerun first run" \
    "$BASEVAR" caller -f "$REF" -o "$TMP_DIR/T13.vcf" \
    -L "$SAMPLES_LIST" --filename-has-samplename --smart-rerun

# T14: Second run (should reuse batchfiles)
run_test "T14" "--smart-rerun second run (reuse)" \
    "$BASEVAR" caller -f "$REF" -o "$TMP_DIR/T14.vcf" \
    -L "$SAMPLES_LIST" --filename-has-samplename --smart-rerun

# =============================================================================
echo ""
echo "================================================================"
echo "  Category 6: Quality Filters"
echo "================================================================"
echo ""

# T15: Strict BQ/MQ
run_test "T15" "Strict filters -Q 30 -q 40" \
    "$BASEVAR" caller -f "$REF" -o "$TMP_DIR/T15.vcf" \
    -L "$SAMPLES_LIST" --filename-has-samplename -Q 30 -q 40

# T16: Relaxed BQ/MQ
run_test "T16" "Relaxed filters -Q 5 -q 10" \
    "$BASEVAR" caller -f "$REF" -o "$TMP_DIR/T16.vcf" \
    -L "$SAMPLES_LIST" --filename-has-samplename -Q 5 -q 10

# =============================================================================
echo ""
echo "================================================================"
echo "  Category 7: Other Parameters"
echo "================================================================"
echo ""

# T17: --ref-bias
run_test "T17" "--ref-bias 0.45" \
    "$BASEVAR" caller -f "$REF" -o "$TMP_DIR/T17.vcf" \
    -L "$SAMPLES_LIST" --filename-has-samplename --ref-bias 0.45

# T18: --max-alleles
run_test "T18" "--max-alleles 3" \
    "$BASEVAR" caller -f "$REF" -o "$TMP_DIR/T18.vcf" \
    -L "$SAMPLES_LIST" --filename-has-samplename --max-alleles 3

# T19: -m min-af
run_test "T19" "-m 0.01 (higher min-af)" \
    "$BASEVAR" caller -f "$REF" -o "$TMP_DIR/T19.vcf" \
    -L "$SAMPLES_LIST" --filename-has-samplename -m 0.01

# T20: Threading
run_test "T20" "-t 4 (multi-thread)" \
    "$BASEVAR" caller -f "$REF" -o "$TMP_DIR/T20.vcf" \
    -L "$SAMPLES_LIST" --filename-has-samplename -t 4

# =============================================================================
echo ""
echo "================================================================"
echo "  Category 8: Combined Workflows"
echo "================================================================"
echo ""

# T21: Full parameter combo
run_test "T21" "Full combo: -G -r -Q -q --gt-mode" \
    "$BASEVAR" caller -f "$REF" -o "$TMP_DIR/T21.vcf" \
    -L "$SAMPLES_LIST" --filename-has-samplename \
    -G "$SAMPLES_GROUP" -r chr1:100-1600 \
    -Q 15 -q 15 --gt-mode posterior

# T22: Multi-region + concat workflow
run_test "T22a" "Multi-region: chr1 only" \
    "$BASEVAR" caller -f "$REF" -o "$TMP_DIR/T22_chr1.vcf" \
    -L "$SAMPLES_LIST" --filename-has-samplename -r chr1

run_test "T22b" "Multi-region: chr2 only" \
    "$BASEVAR" caller -f "$REF" -o "$TMP_DIR/T22_chr2.vcf" \
    -L "$SAMPLES_LIST" --filename-has-samplename -r chr2

run_test "T22c" "Multi-region: chrX only" \
    "$BASEVAR" caller -f "$REF" -o "$TMP_DIR/T22_chrX.vcf" \
    -L "$SAMPLES_LIST" --filename-has-samplename -r chrX

# Concat
run_test "T22d" "Concat chr1+chr2+chrX VCFs" \
    "$BASEVAR" concat -o "$TMP_DIR/T22_concat.vcf" \
    "$TMP_DIR/T22_chr1.vcf" "$TMP_DIR/T22_chr2.vcf" "$TMP_DIR/T22_chrX.vcf"

# =============================================================================
# Evaluation
# =============================================================================
if $RUN_EVALUATION; then
    echo ""
    echo "================================================================"
    echo "  Evaluation: Compare VCF output vs Ground Truth"
    echo "================================================================"
    echo ""

    # Evaluate key multi-sample tests
    echo "  --- T02: Multi-sample via -L ---"
    run_eval "$TMP_DIR/T02.vcf" "T02"

    echo "  --- T09: With -G group stratification ---"
    run_eval "$TMP_DIR/T09.vcf" "T09"

    echo "  --- T11: Posterior mode ---"
    run_eval "$TMP_DIR/T11.vcf" "T11"

    echo "  --- T12: Legacy mode ---"
    run_eval "$TMP_DIR/T12.vcf" "T12"

    echo "  --- T21: Full combo ---"
    run_eval "$TMP_DIR/T21.vcf" "T21"

    echo "  --- T22d: Concat result ---"
    run_eval "$TMP_DIR/T22_concat.vcf" "T22"

    # Compare posterior vs legacy GT concordance
    echo ""
    echo "  --- Posterior vs Legacy comparison ---"
    echo "  (Both should produce similar sensitivity; GT concordance may differ)"
fi

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "================================================================"
echo "  Test Summary"
echo "================================================================"
echo ""
echo "  Total:  $TOTAL"
echo "  PASS:   $PASS"
echo "  FAIL:   $FAIL"
echo "  SKIP:   $SKIP"
echo ""

if [ "$FAIL" -eq 0 ]; then
    echo "  ALL TESTS PASSED"
else
    echo "  $FAIL TEST(S) FAILED"
    echo ""
    echo "  Failed test logs in: $TMP_DIR"
    # Don't clean up on failure
    trap - EXIT
fi

echo ""
echo "================================================================"

exit "$FAIL"
