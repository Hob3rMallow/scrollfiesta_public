#!/usr/bin/env bash
#
# run_all_with_dumps.sh — Run cube_mesh on all cubes with --dump-obj
#
# Usage: scripts/run_all_with_dumps.sh <output_dir>
#
# Runs sequentially. Resumes: skips cubes where output .tif already exists.
# Must be run from the repo root (pipeline binary runs from src/).

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <output_dir>" >&2
    exit 1
fi

PROJ_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PIPELINE="$PROJ_ROOT/src/cube_mesh"
INPUT_DIR="$PROJ_ROOT/nnunet-preds"
OUTPUT_DIR="$1"

if [ ! -x "$PIPELINE" ]; then
    echo "ERROR: $PIPELINE not found. Run 'make cube_mesh' in src/ first." >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

CSV="$OUTPUT_DIR/results.csv"
LOG_ALL="$OUTPUT_DIR/run_all.log"

# Write CSV header if file doesn't exist or is empty
if [ ! -s "$CSV" ]; then
    echo "cube_id,wall_time,s0,s1,s2,s3,s4,s5,s6,status" > "$CSV"
fi

# Count total and already-done
TOTAL=0
DONE=0
TODO=0
for tif in "$INPUT_DIR"/*.tif; do
    [ -f "$tif" ] || continue
    TOTAL=$((TOTAL + 1))
    CUBE_ID="$(basename "$tif" .tif)"
    if [ -f "$OUTPUT_DIR/$CUBE_ID/${CUBE_ID}.tif" ]; then
        DONE=$((DONE + 1))
    else
        TODO=$((TODO + 1))
    fi
done

echo "Total: $TOTAL, Already done: $DONE, Todo: $TODO" | tee -a "$LOG_ALL"
echo "Started: $(date -Iseconds)" | tee -a "$LOG_ALL"

if [ "$TODO" -eq 0 ]; then
    echo "Nothing to do — all cubes already processed."
    exit 0
fi

# Accumulators
RUN_COUNT=0
FAIL_COUNT=0
TOTAL_TIME=0

for tif in "$INPUT_DIR"/*.tif; do
    [ -f "$tif" ] || continue

    CUBE_ID="$(basename "$tif" .tif)"
    CUBE_OUT_DIR="$OUTPUT_DIR/$CUBE_ID"
    OUT_TIF="$CUBE_OUT_DIR/${CUBE_ID}.tif"
    CUBE_LOG="$OUTPUT_DIR/${CUBE_ID}.log"

    # Skip if already done
    if [ -f "$OUT_TIF" ]; then
        continue
    fi

    RUN_COUNT=$((RUN_COUNT + 1))
    mkdir -p "$CUBE_OUT_DIR"

    printf "[%d/%d] %s ... " "$((DONE + RUN_COUNT))" "$TOTAL" "$CUBE_ID"

    # Run pipeline from src/ so relative dependency and output paths stay stable
    START=$(date +%s%N)
    if (cd "$PROJ_ROOT/src" && ./cube_mesh "$tif" "$OUT_TIF" --dump-obj "$OUTPUT_DIR") 2>"$CUBE_LOG"; then
        STATUS="OK"
    else
        STATUS="FAIL"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
    END=$(date +%s%N)

    ELAPSED=$(awk "BEGIN { printf \"%.1f\", ($END - $START) / 1000000000.0 }")
    TOTAL_TIME=$(awk "BEGIN { printf \"%.1f\", $TOTAL_TIME + $ELAPSED }")

    # Parse step timings
    STEP_LINE=$(grep "Step timings:" "$CUBE_LOG" 2>/dev/null || echo "")
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
    echo "$STATUS ${ELAPSED}s"
    echo "[$((DONE + RUN_COUNT))/$TOTAL] $CUBE_ID $STATUS ${ELAPSED}s" >> "$LOG_ALL"
done

echo "" | tee -a "$LOG_ALL"
echo "===============================" | tee -a "$LOG_ALL"
echo "Summary" | tee -a "$LOG_ALL"
echo "===============================" | tee -a "$LOG_ALL"
echo "  Cubes run:     $RUN_COUNT" | tee -a "$LOG_ALL"
echo "  Failures:      $FAIL_COUNT" | tee -a "$LOG_ALL"
echo "  Total time:    ${TOTAL_TIME}s" | tee -a "$LOG_ALL"
if [ "$RUN_COUNT" -gt 0 ]; then
    AVG=$(awk "BEGIN { printf \"%.1f\", $TOTAL_TIME / $RUN_COUNT }")
    echo "  Avg time:      ${AVG}s" | tee -a "$LOG_ALL"
fi
echo "  Results CSV:   $CSV" | tee -a "$LOG_ALL"
echo "  Finished:      $(date -Iseconds)" | tee -a "$LOG_ALL"
echo "===============================" | tee -a "$LOG_ALL"
