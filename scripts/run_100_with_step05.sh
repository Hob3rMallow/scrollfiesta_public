#!/bin/bash
# Run pipeline with step0.5 (2D gap fill) preprocessing on 100 random cubes
# With --slices on gap_fill_2d and --dump-obj on pipeline for full inspection
# Usage: nohup bash scripts/run_100_with_step05.sh &

set -euo pipefail

REPO="/mnt/scroll_data/vesuvius-c"
PRED_DIR="$REPO/nnunet-preds"
GAP_FILL="$REPO/scripts/step0.5-2d-trace/gap_fill_2d"
PIPELINE="$REPO/src/cube_mesh"
OUTDIR="$REPO/output/step0.5_100cube_run"
CUBE_LIST="$OUTDIR/cube_list.txt"
RESULTS_CSV="$OUTDIR/results.csv"
LOG="$OUTDIR/run.log"
DUMP_DIR="$OUTDIR/dumps"

mkdir -p "$OUTDIR" "$DUMP_DIR"

# Select 100 random cubes (use saved list if already exists, for reproducibility)
if [ ! -f "$CUBE_LIST" ]; then
    ls "$PRED_DIR"/*.tif | shuf -n 100 | sed 's/.*\///;s/\.tif//' | sort > "$CUBE_LIST"
fi

N_CUBES=$(wc -l < "$CUBE_LIST")
echo "=== step0.5 + pipeline batch run: $N_CUBES cubes ===" | tee "$LOG"
echo "Started: $(date)" | tee -a "$LOG"
echo "Output dir: $OUTDIR" | tee -a "$LOG"
echo "Dump dir:   $DUMP_DIR" | tee -a "$LOG"
echo "" | tee -a "$LOG"
echo "cube_id,step05_rc,step05_time_s,pipeline_rc,pipeline_time_s,total_time_s,status" > "$RESULTS_CSV"

DONE=0
FAIL=0
IDX=0

while read -r CUBE; do
    IDX=$((IDX + 1))
    echo "" | tee -a "$LOG"
    echo "[$IDX/$N_CUBES] Processing $CUBE ..." | tee -a "$LOG"

    INPUT="$PRED_DIR/${CUBE}.tif"
    SLICES_DIR="$OUTDIR/${CUBE}_step0.5_slices"
    GAPFILL_OUT="$OUTDIR/${CUBE}_step0.5_gapfilled.tif"
    PIPELINE_OUT="$OUTDIR/${CUBE}_step6_output.tif"

    if [ ! -f "$INPUT" ]; then
        echo "  SKIP: input not found: $INPUT" | tee -a "$LOG"
        echo "$CUBE,-1,0,-1,0,0,SKIP_NO_INPUT" >> "$RESULTS_CSV"
        FAIL=$((FAIL + 1))
        continue
    fi

    # Step 0.5: 2D gap fill with --slices for per-slice inspection
    mkdir -p "$SLICES_DIR"
    T0=$(date +%s.%N)
    echo "  Step 0.5: gap fill (--slices $SLICES_DIR) ..." | tee -a "$LOG"
    S05_RC=0
    "$GAP_FILL" "$INPUT" "$GAPFILL_OUT" --slices "$SLICES_DIR" 2>>"$LOG" || S05_RC=$?
    T1=$(date +%s.%N)
    S05_TIME=$(echo "$T1 - $T0" | bc)
    echo "  Step 0.5: rc=$S05_RC time=${S05_TIME}s" | tee -a "$LOG"

    if [ $S05_RC -ne 0 ] || [ ! -f "$GAPFILL_OUT" ]; then
        echo "  Step 0.5 FAILED, falling back to raw input" | tee -a "$LOG"
        GAPFILL_OUT="$INPUT"
    fi

    # Pipeline (steps 0-6) with --dump-obj for mesh inspection
    T2=$(date +%s.%N)
    echo "  Pipeline: steps 0-6 (--dump-obj $DUMP_DIR) ..." | tee -a "$LOG"
    PIPE_RC=0
    (cd "$REPO/src" && "$PIPELINE" "$GAPFILL_OUT" "$PIPELINE_OUT" --dump-obj "$DUMP_DIR") 2>>"$LOG" || PIPE_RC=$?
    T3=$(date +%s.%N)
    PIPE_TIME=$(echo "$T3 - $T2" | bc)
    TOTAL_TIME=$(echo "$T3 - $T0" | bc)

    if [ $PIPE_RC -eq 0 ]; then
        STATUS="OK"
        DONE=$((DONE + 1))
    else
        STATUS="FAIL"
        FAIL=$((FAIL + 1))
    fi

    echo "  Pipeline: rc=$PIPE_RC time=${PIPE_TIME}s total=${TOTAL_TIME}s status=$STATUS" | tee -a "$LOG"
    echo "$CUBE,$S05_RC,$S05_TIME,$PIPE_RC,$PIPE_TIME,$TOTAL_TIME,$STATUS" >> "$RESULTS_CSV"

done < "$CUBE_LIST"

echo "" | tee -a "$LOG"
echo "=== SUMMARY ===" | tee -a "$LOG"
echo "Total: $N_CUBES  OK: $DONE  FAIL: $FAIL" | tee -a "$LOG"
echo "Finished: $(date)" | tee -a "$LOG"
echo "" | tee -a "$LOG"
echo "Results CSV: $RESULTS_CSV" | tee -a "$LOG"
echo "Log:         $LOG" | tee -a "$LOG"
echo "Output TIFs: $OUTDIR/*_step6_output.tif" | tee -a "$LOG"
echo "Gap-fill:    $OUTDIR/*_step0.5_gapfilled.tif" | tee -a "$LOG"
echo "Slices:      $OUTDIR/*_step0.5_slices/" | tee -a "$LOG"
echo "OBJ dumps:   $DUMP_DIR/" | tee -a "$LOG"
