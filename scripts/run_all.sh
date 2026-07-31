#!/usr/bin/env bash
#
# run_all.sh — Batch runner for cube_mesh.
#
# Usage:
#   ./run_all.sh <input_dir> [--kaggle]
#
# Runs cube_mesh on each .tif file in input_dir.
# Outputs go to <input_dir>/../output/ (sibling of input dir).
# Results logged to results.csv.
#
# --kaggle flag: wraps each run in taskset -c 0,1 to simulate
#   Kaggle 2-core environment.

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <input_dir> [--kaggle]" >&2
    exit 1
fi

INPUT_DIR="$1"
shift

KAGGLE=0
while [ $# -gt 0 ]; do
    case "$1" in
        --kaggle) KAGGLE=1 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
    shift
done

# Resolve paths
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PIPELINE="$SCRIPT_DIR/../src/cube_mesh"

if [ ! -x "$PIPELINE" ]; then
    echo "ERROR: $PIPELINE not found. Run 'make cube_mesh' first." >&2
    exit 1
fi

OUTPUT_DIR="$(dirname "$INPUT_DIR")/output"
mkdir -p "$OUTPUT_DIR"

CSV="$OUTPUT_DIR/results.csv"
echo "cube_id,total_time,s0,s1,s2,s3,s4,s5,s6,status" > "$CSV"

TOTAL_CUBES=0
TOTAL_TIME=0.0
MAX_TIME=0.0
OVER_180=0

for tif in "$INPUT_DIR"/*.tif; do
    [ -f "$tif" ] || continue

    CUBE_ID="$(basename "$tif" .tif)"
    OUT_TIF="$OUTPUT_DIR/${CUBE_ID}.tif"
    LOG="$OUTPUT_DIR/${CUBE_ID}.log"

    echo "=== Processing $CUBE_ID ==="

    # Build command
    CMD="$PIPELINE $tif $OUT_TIF"
    if [ "$KAGGLE" -eq 1 ]; then
        CMD="taskset -c 0,1 $CMD"
    fi

    # Run and capture timing + status
    START=$(date +%s%N)
    if $CMD 2>"$LOG"; then
        STATUS="OK"
    else
        STATUS="FALLBACK"
    fi
    END=$(date +%s%N)

    # Compute wall time in seconds (with decimals)
    ELAPSED=$(awk "BEGIN { printf \"%.3f\", ($END - $START) / 1000000000.0 }")

    # Parse step timings from the log
    STEP_LINE=$(grep "Step timings:" "$LOG" 2>/dev/null || echo "")
    if [ -n "$STEP_LINE" ]; then
        S0=$(echo "$STEP_LINE" | sed -n 's/.*s0=\([0-9.]*\).*/\1/p')
        S1=$(echo "$STEP_LINE" | sed -n 's/.*s1=\([0-9.]*\).*/\1/p')
        S2=$(echo "$STEP_LINE" | sed -n 's/.*s2=\([0-9.]*\).*/\1/p')
        S3=$(echo "$STEP_LINE" | sed -n 's/.*s3=\([0-9.]*\).*/\1/p')
        S4=$(echo "$STEP_LINE" | sed -n 's/.*s4=\([0-9.]*\).*/\1/p')
        S5=$(echo "$STEP_LINE" | sed -n 's/.*s5=\([0-9.]*\).*/\1/p')
        S6=$(echo "$STEP_LINE" | sed -n 's/.*s6=\([0-9.]*\).*/\1/p')
    else
        S0=""; S1=""; S2=""; S3=""; S4=""; S5=""; S6=""
    fi

    echo "$CUBE_ID,$ELAPSED,$S0,$S1,$S2,$S3,$S4,$S5,$S6,$STATUS" >> "$CSV"

    # Accumulate stats
    TOTAL_CUBES=$((TOTAL_CUBES + 1))
    TOTAL_TIME=$(awk "BEGIN { printf \"%.3f\", $TOTAL_TIME + $ELAPSED }")
    MAX_TIME=$(awk "BEGIN { m=$MAX_TIME; if ($ELAPSED > m) m=$ELAPSED; printf \"%.3f\", m }")
    if awk "BEGIN { exit ($ELAPSED > 180.0) ? 0 : 1 }"; then
        OVER_180=$((OVER_180 + 1))
    fi

    echo "  -> $STATUS  ${ELAPSED}s"
done

# Summary
echo ""
echo "==============================="
echo "Summary"
echo "==============================="
echo "  Total cubes:   $TOTAL_CUBES"
if [ "$TOTAL_CUBES" -gt 0 ]; then
    AVG_TIME=$(awk "BEGIN { printf \"%.3f\", $TOTAL_TIME / $TOTAL_CUBES }")
    echo "  Avg time:      ${AVG_TIME}s"
fi
echo "  Max time:      ${MAX_TIME}s"
echo "  Over 180s:     $OVER_180"
echo "  Results:       $CSV"
echo "==============================="
