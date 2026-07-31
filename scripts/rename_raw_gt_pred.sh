#!/bin/bash
# rename_raw_gt_pred.sh — Rename existing raw-gt-pred files to include cube ID prefix
#
# Before: raw-gt-pred/{cube_id}/gt.obj, pred.obj, raw.tif
# After:  raw-gt-pred/{cube_id}/{cube_id}-gt.obj, {cube_id}-pred.obj, {cube_id}-raw.tif
#
# Idempotent: skips already-renamed files.
# Symlinks are re-created pointing to project-local data/train_images/{cube_id}.tif

set -euo pipefail

PROJ_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RAW_GT_PRED="$PROJ_ROOT/raw-gt-pred"

if [ ! -d "$RAW_GT_PRED" ]; then
    echo "ERROR: $RAW_GT_PRED not found" >&2
    exit 1
fi

renamed=0
skipped=0

for dir in "$RAW_GT_PRED"/*/; do
    cube_id="$(basename "$dir")"

    # Rename gt.obj -> {cube_id}-gt.obj
    if [ -f "$dir/gt.obj" ] && [ ! -f "$dir/${cube_id}-gt.obj" ]; then
        mv "$dir/gt.obj" "$dir/${cube_id}-gt.obj"
        renamed=$((renamed + 1))
    else
        skipped=$((skipped + 1))
    fi

    # Rename pred.obj -> {cube_id}-pred.obj
    if [ -f "$dir/pred.obj" ] && [ ! -f "$dir/${cube_id}-pred.obj" ]; then
        mv "$dir/pred.obj" "$dir/${cube_id}-pred.obj"
        renamed=$((renamed + 1))
    else
        skipped=$((skipped + 1))
    fi

    # Re-create raw.tif symlink -> {cube_id}-raw.tif pointing to project-local path
    if [ -L "$dir/raw.tif" ] || [ -f "$dir/raw.tif" ]; then
        rm -f "$dir/raw.tif"
    fi
    if [ ! -L "$dir/${cube_id}-raw.tif" ]; then
        ln -s "$PROJ_ROOT/data/train_images/${cube_id}.tif" "$dir/${cube_id}-raw.tif"
        renamed=$((renamed + 1))
    else
        skipped=$((skipped + 1))
    fi
done

echo "rename_raw_gt_pred.sh: renamed=$renamed skipped=$skipped"
