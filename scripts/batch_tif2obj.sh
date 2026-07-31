#!/bin/bash
# batch_tif2obj.sh — Convert all 786 cubes to OBJ meshes
#
# For each cube, generates:
#   raw-gt-pred/{cube_id}/{cube_id}-gt.obj   (from data/train_labels/{cube_id}.tif)
#   raw-gt-pred/{cube_id}/{cube_id}-pred.obj (from nnunet-preds/{cube_id}.tif)
#   raw-gt-pred/{cube_id}/{cube_id}-raw.tif  (symlink to data/train_images/{cube_id}.tif)
#
# Resume support: skips cubes that already have both OBJs.
# Usage: scripts/batch_tif2obj.sh [--jobs N] [--force]

set -euo pipefail

PROJ_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TIF2OBJ="$PROJ_ROOT/tif2obj"
RAW_GT_PRED="$PROJ_ROOT/raw-gt-pred"
LABELS_DIR="$PROJ_ROOT/data/train_labels"
PREDS_DIR="$PROJ_ROOT/nnunet-preds"
IMAGES_DIR="$PROJ_ROOT/data/train_images"

JOBS=16
FORCE=0

while [ $# -gt 0 ]; do
    case "$1" in
        --jobs) JOBS="$2"; shift 2 ;;
        --force) FORCE=1; shift ;;
        *) echo "Usage: $0 [--jobs N] [--force]" >&2; exit 1 ;;
    esac
done

# Validate
if [ ! -x "$TIF2OBJ" ]; then
    echo "ERROR: tif2obj binary not found at $TIF2OBJ" >&2
    echo "Build it first: gcc -O2 -o tif2obj scripts/tif2obj.c -ltiff -lm" >&2
    exit 1
fi

if [ ! -d "$LABELS_DIR" ]; then
    echo "ERROR: labels dir not found: $LABELS_DIR" >&2
    exit 1
fi

if [ ! -d "$PREDS_DIR" ]; then
    echo "ERROR: predictions dir not found: $PREDS_DIR" >&2
    exit 1
fi

# Phase 1: Create directories and symlinks (serial, fast)
echo "Phase 1: Creating directories and symlinks..."
total=0
for tif in "$LABELS_DIR"/*.tif; do
    cube_id="$(basename "$tif" .tif)"
    mkdir -p "$RAW_GT_PRED/$cube_id"

    # Create raw.tif symlink if missing
    if [ ! -L "$RAW_GT_PRED/$cube_id/${cube_id}-raw.tif" ]; then
        ln -sf "$IMAGES_DIR/${cube_id}.tif" "$RAW_GT_PRED/$cube_id/${cube_id}-raw.tif"
    fi
    total=$((total + 1))
done
echo "  $total cube directories ready"

# Phase 2: Generate OBJs (parallel)
echo "Phase 2: Generating OBJ meshes with $JOBS parallel jobs..."

# Build work list
work_file=$(mktemp)
gt_todo=0
pred_todo=0

for tif in "$LABELS_DIR"/*.tif; do
    cube_id="$(basename "$tif" .tif)"

    # GT OBJ
    gt_out="$RAW_GT_PRED/$cube_id/${cube_id}-gt.obj"
    if [ "$FORCE" -eq 1 ] || [ ! -f "$gt_out" ]; then
        echo "GT $cube_id $tif $gt_out" >> "$work_file"
        gt_todo=$((gt_todo + 1))
    fi

    # Pred OBJ
    pred_tif="$PREDS_DIR/${cube_id}.tif"
    pred_out="$RAW_GT_PRED/$cube_id/${cube_id}-pred.obj"
    if [ -f "$pred_tif" ]; then
        if [ "$FORCE" -eq 1 ] || [ ! -f "$pred_out" ]; then
            echo "PRED $cube_id $pred_tif $pred_out" >> "$work_file"
            pred_todo=$((pred_todo + 1))
        fi
    fi
done

work_total=$((gt_todo + pred_todo))
echo "  GT todo: $gt_todo, Pred todo: $pred_todo, Total: $work_total"

if [ "$work_total" -eq 0 ]; then
    echo "Nothing to do — all OBJs already exist."
    rm -f "$work_file"
    exit 0
fi

# Progress tracking
done_count=0
fail_count=0
log_file="$RAW_GT_PRED/batch_tif2obj.log"
echo "# batch_tif2obj started $(date -Iseconds)" > "$log_file"
echo "# jobs=$JOBS gt_todo=$gt_todo pred_todo=$pred_todo" >> "$log_file"

process_one() {
    local kind="$1" cube_id="$2" in_tif="$3" out_obj="$4"
    local t0 t1 elapsed
    t0=$(date +%s)
    local extra_args=""
    if [ "$kind" = "GT" ]; then
        extra_args="--fg-val 1"
    fi
    if "$TIF2OBJ" "$in_tif" "$out_obj" $extra_args 2>&1; then
        t1=$(date +%s)
        elapsed=$((t1 - t0))
        echo "OK ${kind} ${cube_id} ${elapsed}s"
    else
        t1=$(date +%s)
        elapsed=$((t1 - t0))
        echo "FAIL ${kind} ${cube_id} ${elapsed}s"
    fi
}
export -f process_one
export TIF2OBJ

# Run in parallel with xargs
cat "$work_file" | xargs -P "$JOBS" -L 1 bash -c 'process_one "$@"' _ 2>&1 | \
    while IFS= read -r line; do
        echo "$line" >> "$log_file"
        done_count=$((done_count + 1))
        status=$(echo "$line" | cut -d' ' -f1)
        kind=$(echo "$line" | cut -d' ' -f2)
        cid=$(echo "$line" | cut -d' ' -f3)
        elapsed=$(echo "$line" | cut -d' ' -f4)
        if [ "$status" = "FAIL" ]; then
            fail_count=$((fail_count + 1))
        fi
        printf "\r  [%d/%d] %s %s %s %s    " "$done_count" "$work_total" "$status" "$kind" "$cid" "$elapsed"
    done

echo ""
rm -f "$work_file"

echo "# batch_tif2obj finished $(date -Iseconds)" >> "$log_file"
echo "Done. Log: $log_file"

# Summary
total_dirs=$(ls -d "$RAW_GT_PRED"/*/ 2>/dev/null | wc -l)
total_gt=$(find "$RAW_GT_PRED" -name "*-gt.obj" | wc -l)
total_pred=$(find "$RAW_GT_PRED" -name "*-pred.obj" | wc -l)
total_raw=$(find "$RAW_GT_PRED" -name "*-raw.tif" | wc -l)
echo "Summary: $total_dirs dirs, $total_gt GT OBJs, $total_pred Pred OBJs, $total_raw raw symlinks"
