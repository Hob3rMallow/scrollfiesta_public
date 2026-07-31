#!/usr/bin/env bash
# debug_block.sh -- isolate ONE 4x5x5 weld tile of the 4x21x21 grid into its own
# output folder and build a QC gallery, so a single block can be debugged in
# isolation (re-weldable in seconds from the cached per-cube dumps).
#
# Re-welds the chosen subgrid from the shared dump dir (default stage decim10),
# 50%-decimates it, dumps every weld stage (concat..final + per-hole OBJs), then
# renders a normal-coloured gallery: a whole-block overview (full + d50) plus one
# tile per constituent cube, assembled into a per-z-layer PDF. A clean wrap reads
# as a smooth colour gradient; a fold / bridge / crumpled sheet shows as chaotic
# speckle (see src/tools/mesh_render.c).
#
# Bash-as-runner only (CLAUDE.md): grid_weld / qslim_obj / mesh_render do the work.
#
# Usage: bash scripts/debug_block.sh <tag> <z0> <z1> <y0> <y1> <x0> <x1> [stage] [size]
#   e.g. bash scripts/debug_block.sh y2048_x1536 4352 4736 2048 2560 1536 2048
set -u
TAG="${1:?tag e.g. y2048_x1536}"
Z0="${2:?z0}"; Z1="${3:?z1}"; Y0="${4:?y0}"; Y1="${5:?y1}"; X0="${6:?x0}"; X1="${7:?x1}"
STAGE="${8:-decim10}"
SIZE="${9:-1024}"

SRC=output/grid_4x21x21_foldfix/dump
OUT=output/blk_${TAG}_debug
W=./build/Release/grid_weld.exe
Q=./build/Release/qslim_obj.exe
R=./build/Release/mesh_render.exe
PNG="$OUT/gallery/png"
mkdir -p "$OUT/weld_stages" "$PNG"

echo "=== [1/4] re-weld block $TAG  z[$Z0,$Z1] y[$Y0,$Y1] x[$X0,$X1]  stage=$STAGE  $(date) ==="
"$W" "$SRC" "$OUT/blk_${TAG}.obj" --stage "$STAGE" \
     --subgrid "$Z0" "$Z1" "$Y0" "$Y1" "$X0" "$X1" \
     --dump-stages "$OUT/weld_stages" --stage-prefix "blk_${TAG}" \
     > "$OUT/blk_${TAG}.weldlog" 2>&1
echo "  weld rc=$?  -> $OUT/blk_${TAG}.obj"
grep -E '^grid_weld: --subgrid|unpaired|non_manifold|same_dir' "$OUT/blk_${TAG}.weldlog" 2>/dev/null | head

echo "=== [2/4] qslim 50%  $(date) ==="
[ -f "$OUT/blk_${TAG}.obj" ] && \
  "$Q" "$OUT/blk_${TAG}.obj" "$OUT/blk_${TAG}_d50.obj" 0.50 > "$OUT/blk_${TAG}.qslimlog" 2>&1
echo "  d50 rc=$?"

echo "=== [3/4] render block overview + constituent cubes @ ${SIZE}px  $(date) ==="
"$R" "$OUT/blk_${TAG}.obj"      "$OUT/gallery/block_${TAG}_full.png" --size "$SIZE" >/dev/null 2>&1
[ -f "$OUT/blk_${TAG}_d50.obj" ] && \
  "$R" "$OUT/blk_${TAG}_d50.obj" "$OUT/gallery/block_${TAG}_d50.png"  --size "$SIZE" >/dev/null 2>&1

MAXJ=24; launched=0
for ((z=Z0; z<=Z1; z+=128)); do
  for ((y=Y0; y<=Y1; y+=128)); do
    for ((x=X0; x<=X1; x+=128)); do
      id=$(printf 'z%05d_y%05d_x%05d' "$z" "$y" "$x")
      obj="$SRC/$id/${id}_${STAGE}/${id}_${STAGE}_all.obj"
      [ -f "$obj" ] || continue
      nf=$(grep -c '^f ' "$obj")
      "$R" "$obj" "$PNG/${id}_${nf}f.png" --size "$SIZE" >/dev/null 2>&1 &
      launched=$((launched+1))
      [ $((launched % MAXJ)) -eq 0 ] && wait
    done
  done
done
wait
echo "  rendered $launched cube tiles + 2 block overviews"

echo "=== [4/4] assemble gallery PDF  $(date) ==="
python scripts/make_gallery_pdf.py "$PNG" "$OUT/gallery/cubes_${TAG}.pdf" --cols 5 --rows 5 --tile 380
echo "=== done -> $OUT  $(date) ==="
ls -la "$OUT"
