#!/usr/bin/env bash
# weld_decim_blocks.sh — weld the decim10 cubes of the 4x21x21 grid into
# NON-OVERLAPPING 4x5x5 tiles (y/x tile as 5+5+5+5+1, z is the full 4), then
# qslim each welded tile by 50%. Skips tiles already finished. Parallel.
#
# Output (in <grid>/blocks/): blk_y<Y0>_x<X0>.obj (welded, hole-filled) and
# blk_y<Y0>_x<X0>_d50.obj (50% decimated), + .weld_report.json / logs each.
set -u
O=output/grid_4x21x21_foldfix
W=./build/Release/grid_weld.exe
Q=./build/Release/qslim_obj.exe
B="$O/blocks"; mkdir -p "$B"
Z0=4352; Z1=4736
ys=(2048 2688 3328 3968 4608); ye=(2560 3200 3840 4480 4608)
xs=(1536 2176 2816 3456 4096); xe=(2048 2688 3328 3968 4096)
MAXJ=8; launched=0

for yi in 0 1 2 3 4; do for xi in 0 1 2 3 4; do
  y0=${ys[$yi]}; y1=${ye[$yi]}; x0=${xs[$xi]}; x1=${xe[$xi]}; tag="y${y0}_x${x0}"
  [ -f "$B/blk_${tag}_d50.obj" ] && { echo "skip done: $tag"; continue; }
  {
    "$W" "$O/dump" "$B/blk_${tag}.obj" --stage decim10 \
         --subgrid $Z0 $Z1 $y0 $y1 $x0 $x1 > "$B/blk_${tag}.weldlog" 2>&1
    [ -f "$B/blk_${tag}.obj" ] && \
      "$Q" "$B/blk_${tag}.obj" "$B/blk_${tag}_d50.obj" 0.50 > "$B/blk_${tag}.qslimlog" 2>&1
  } &
  launched=$((launched+1))
  [ $((launched % MAXJ)) -eq 0 ] && { wait; echo "  ...$launched tiles dispatched"; }
done; done
wait

echo "=== weld+decimate-50% block summary ==="
jget(){ grep -oE "\"$2\":[ ]*[0-9]+" "$1" 2>/dev/null | grep -oE '[0-9]+$'; }
for yi in 0 1 2 3 4; do for xi in 0 1 2 3 4; do
  y0=${ys[$yi]}; x0=${xs[$xi]}; tag="y${y0}_x${x0}"; rep="$B/blk_${tag}.obj.weld_report.json"
  if [ -f "$rep" ]; then
    cubes=$(jget "$rep" cubes_processed); wf=$(jget "$rep" total_unique_faces)
    nm=$(jget "$rep" non_manifold); sd=$(jget "$rep" same_dir_pairs)
    df=$(grep -oE 'f [0-9]+->[0-9]+' "$B/blk_${tag}.qslimlog" 2>/dev/null | grep -oE '[0-9]+$' | head -1)
    printf "  %-16s cubes=%-4s weld=%-8s nm=%-3s sd=%-5s -> d50=%s f\n" \
           "$tag" "${cubes:-?}" "${wf:-?}" "${nm:-?}" "${sd:-?}" "${df:--}"
  else echo "  $tag : empty / no weld (garbage slab or no cubes)"; fi
done; done