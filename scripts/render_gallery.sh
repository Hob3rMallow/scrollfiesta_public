#!/usr/bin/env bash
# render_gallery.sh -- render every final-stage per-cube mesh (step12_final) to a
# PNG via mesh_render, concurrently, then assemble a labelled multi-page PDF.
#
# Usage: bash scripts/render_gallery.sh <dump_dir> <gallery_dir> [size] [jobs]
set -u
DUMP="${1:?dump_dir}"
GAL="${2:?gallery_dir}"
export SIZE="${3:-1024}"
JOBS="${4:-64}"
export MR="build/Release/mesh_render.exe"
export PNGDIR="$GAL/png"
mkdir -p "$PNGDIR"

mapfile -t meshes < <(ls "$DUMP"/*/*_step12_final/*_step12_final_all.obj 2>/dev/null)
echo "rendering ${#meshes[@]} meshes at ${SIZE}px, ${JOBS}-way -> $PNGDIR  $(date)"

printf '%s\n' "${meshes[@]}" | xargs -P "$JOBS" -I {} bash -c '
  src="$1"; id=$(basename "$src" _step12_final_all.obj)
  "$MR" "$src" "$PNGDIR/$id.png" --size "$SIZE" >/dev/null 2>&1 || echo "FAIL $id"
' _ {}

npng=$(ls "$PNGDIR"/*.png 2>/dev/null | wc -l | tr -d ' ')
echo "rendered $npng PNGs  $(date)"

echo "assembling PDF..."
python scripts/make_gallery_pdf.py "$PNGDIR" "$GAL/gallery.pdf"
echo "done $(date)"
