#!/usr/bin/env bash
# mls_radius_sweep.sh -- sweep the LOP/MLS kernel radius (MLS_RADIUS_VOX env knob
# in mls_project.c) on ONE cube and compare the result against the default R=12.
#
# For each radius it re-extracts the cube but only as far as the step1_bpa dump
# (the direct test of whether a smaller kernel keeps two close wraps apart at a
# fold), then KILLS the run -- so each radius costs ~2 min, one process at a time
# (no oversubscription, no exe lock). Renders the post-MLS cloud + step1 mesh.
#
# Bash-as-runner only: cube_mesh / mesh_render do the work.
# Usage: bash scripts/extract/mls_radius_sweep.sh <cube_id> <rad1> [rad2 ...]
set -u
c="${1:?cube_id e.g. z04352_y02048_x01920}"; shift
RADII=("$@"); [ ${#RADII[@]} -gt 0 ] || { echo "give >=1 radius"; exit 1; }
TIFF="PHerc0139-4x21x21-cubes/cubes_PRED/$c.tif"
CM=./build/Release/cube_mesh.exe
MR=./build/Release/mesh_render.exe
BASE="output/blk_y2048_x1536_debug/mls_sweep"
[ -f "$TIFF" ] || { echo "no tiff: $TIFF"; exit 1; }

for RAD in "${RADII[@]}"; do
  DD="$BASE/r$RAD/dump"; OUT="$BASE/r$RAD"; mkdir -p "$DD" "$OUT/png"
  s1="$DD/$c/${c}_step1_bpa/${c}_step1_bpa_all.obj"
  rm -f "$s1"
  echo "=== R=$RAD : extract -> step1_bpa  $(date +%H:%M:%S) ==="
  MLS_RADIUS_VOX="$RAD" "$CM" "$TIFF" "$OUT/out.tif" \
      --dump-obj "$DD" --halo 13 --no-qem --no-timeout > "$OUT/reextract.log" 2>&1 &
  pid=$!
  ok=0
  for _ in $(seq 1 360); do
    if [ -f "$s1" ]; then ok=1; sleep 3; break; fi
    kill -0 "$pid" 2>/dev/null || break
    sleep 1
  done
  kill "$pid" 2>/dev/null; sleep 1; kill -9 "$pid" 2>/dev/null; wait "$pid" 2>/dev/null
  grep -m1 "MLS.*radius override" "$OUT/reextract.log" || echo "  (no override log!)"
  nf=$(grep -c '^f ' "$s1" 2>/dev/null || echo 0)
  echo "  R=$RAD dumped=$ok step1_bpa faces=$nf"
  # render post-MLS combined cloud + step1_bpa mesh
  S0="$DD/$c/${c}_step0_mls"
  cat "$S0"/*_step0_mls_[0-9][0-9][0-9].obj 2>/dev/null > "$OUT/step0_POST_all.obj"
  "$MR" "$OUT/step0_POST_all.obj" "$OUT/png/r${RAD}_step0_POST.png" --size 1280 --ptsize 2 >/dev/null 2>&1
  "$MR" "$s1"                     "$OUT/png/r${RAD}_step1_bpa.png"  --size 1024 >/dev/null 2>&1
done
echo "=== SWEEP DONE  $(date +%H:%M:%S) ==="
# safety: no straggler cube_mesh should remain
