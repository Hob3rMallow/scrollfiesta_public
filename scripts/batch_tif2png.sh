#!/bin/bash
# batch_tif2png.sh — Generate slice visualization PNGs for all cubes
#
# For each cube, generates:
#   raw-gt-pred/{cube_id}/{cube_id}-viz.png
#
# Inputs:
#   data/train_images/{cube_id}.tif  (raw CT)
#   data/train_labels/{cube_id}.tif  (GT labels)
#   nnunet-preds/{cube_id}.tif       (predictions)
#
# Resume support: skips cubes that already have a viz PNG.
# Usage: scripts/batch_tif2png.sh [--jobs N] [--force]

set -euo pipefail

PROJ_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TIF2PNG="$PROJ_ROOT/tif2png"
RAW_GT_PRED="$PROJ_ROOT/raw-gt-pred"
IMAGES_DIR="$PROJ_ROOT/data/train_images"
LABELS_DIR="$PROJ_ROOT/data/train_labels"
PREDS_DIR="$PROJ_ROOT/nnunet-preds"

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
if [ ! -x "$TIF2PNG" ]; then
    echo "ERROR: tif2png binary not found at $TIF2PNG" >&2
    echo "Build it first: gcc -O2 -o tif2png scripts/tif2png.c -ltiff -lpng16 -lm" >&2
    exit 1
fi

for dir in "$IMAGES_DIR" "$LABELS_DIR" "$PREDS_DIR"; do
    if [ ! -d "$dir" ]; then
        echo "ERROR: directory not found: $dir" >&2
        exit 1
    fi
done

# Build work list
work_file=$(mktemp)
todo=0
skip=0

for tif in "$LABELS_DIR"/*.tif; do
    cube_id="$(basename "$tif" .tif)"
    out_png="$RAW_GT_PRED/$cube_id/${cube_id}-viz.png"

    raw_tif="$IMAGES_DIR/${cube_id}.tif"
    pred_tif="$PREDS_DIR/${cube_id}.tif"

    # All 3 inputs must exist
    if [ ! -f "$raw_tif" ] || [ ! -f "$tif" ] || [ ! -f "$pred_tif" ]; then
        continue
    fi

    # Skip if output exists (unless --force)
    if [ "$FORCE" -eq 0 ] && [ -f "$out_png" ]; then
        skip=$((skip + 1))
        continue
    fi

    # Ensure output directory exists
    mkdir -p "$RAW_GT_PRED/$cube_id"

    echo "$cube_id $raw_tif $tif $pred_tif $out_png" >> "$work_file"
    todo=$((todo + 1))
done

echo "Todo: $todo, Skipped (already exist): $skip"

if [ "$todo" -eq 0 ]; then
    echo "Nothing to do — all viz PNGs already exist."
    rm -f "$work_file"
    exit 0
fi

# Progress tracking
log_file="$RAW_GT_PRED/batch_tif2png.log"
echo "# batch_tif2png started $(date -Iseconds)" > "$log_file"
echo "# jobs=$JOBS todo=$todo" >> "$log_file"

process_one() {
    local cube_id="$1" raw_tif="$2" gt_tif="$3" pred_tif="$4" out_png="$5"
    local t0 t1 elapsed
    t0=$(date +%s)
    if "$TIF2PNG" "$raw_tif" "$gt_tif" "$pred_tif" "$out_png" "$cube_id" 2>/dev/null; then
        t1=$(date +%s)
        elapsed=$((t1 - t0))
        echo "OK ${cube_id} ${elapsed}s"
    else
        t1=$(date +%s)
        elapsed=$((t1 - t0))
        echo "FAIL ${cube_id} ${elapsed}s"
    fi
}
export -f process_one
export TIF2PNG

echo "Generating viz PNGs with $JOBS parallel jobs..."

done_count=0
fail_count=0

cat "$work_file" | xargs -P "$JOBS" -L 1 bash -c 'process_one "$@"' _ 2>&1 | \
    while IFS= read -r line; do
        echo "$line" >> "$log_file"
        done_count=$((done_count + 1))
        status=$(echo "$line" | cut -d' ' -f1)
        cid=$(echo "$line" | cut -d' ' -f2)
        elapsed=$(echo "$line" | cut -d' ' -f3)
        if [ "$status" = "FAIL" ]; then
            fail_count=$((fail_count + 1))
        fi
        printf "\r  [%d/%d] %s %s %s    " "$done_count" "$todo" "$status" "$cid" "$elapsed"
    done

echo ""
rm -f "$work_file"

echo "# batch_tif2png finished $(date -Iseconds)" >> "$log_file"
echo "Done. Log: $log_file"

# Summary
total_png=$(find "$RAW_GT_PRED" -name "*-viz.png" | wc -l)
echo "Summary: $total_png viz PNGs in $RAW_GT_PRED"
