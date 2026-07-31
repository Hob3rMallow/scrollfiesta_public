#!/usr/bin/env bash
# weld_blocks.sh -- weld the 4x21x21 grid as disjoint 4x5x5 tiles via
# grid_weld --subgrid, in parallel, then audit every tile. NO whole-grid weld.
#
# Tiling: all 4 z-slabs; 4 y-bands x 4 x-bands (the last band of each axis is
# 6 wide, since 21 = 5+5+5+6). Set FORCE=1 to re-weld tiles whose .obj exists.
#
# Usage: bash scripts/weld_blocks.sh <dump_dir> <blocks_dir>
set -u
DUMP="${1:?dump_dir}"
OUTDIR="${2:?blocks_dir}"
WELD=build/Release/grid_weld.exe
MC=build/Release/manifold_check.exe
Z0=4352; Z1=4736
YB=("2048 2560" "2688 3200" "3328 3840" "3968 4608")
XB=("1536 2048" "2176 2688" "2816 3328" "3456 4096")
mkdir -p "$OUTDIR"

pids=()
launched=0
for yb in "${YB[@]}"; do
  read -r y0 y1 <<< "$yb"
  for xb in "${XB[@]}"; do
    read -r x0 x1 <<< "$xb"
    name="block_y${y0}_x${x0}"
    obj="$OUTDIR/$name.obj"
    if [ -f "$obj" ] && [ "${FORCE:-0}" != "1" ]; then echo "skip (exists): $name"; continue; fi
    ( echo "start $(date)" > "$OUTDIR/$name.weldlog"
      "$WELD" "$DUMP" "$obj" --subgrid "$Z0" "$Z1" "$y0" "$y1" "$x0" "$x1" >> "$OUTDIR/$name.weldlog" 2>&1
      echo "WELD_EXIT=$? $(date)" >> "$OUTDIR/$name.weldlog" ) &
    pids+=($!)
    launched=$((launched+1))
  done
done
echo "launched $launched block weld(s); waiting for completion..."
wait
echo "all block welds done $(date)"

echo ""
echo "================ 16-BLOCK AUDIT ================"
ls "$OUTDIR"/block_*.obj 2>/dev/null | sort > /tmp/wb_blocks.txt
"$MC" --list /tmp/wb_blocks.txt 2>/dev/null
echo ""
echo "--- per-block exit codes (relaxed gate: same_dir warns, non-manifold fails) ---"
for f in $(ls "$OUTDIR"/block_*.weldlog 2>/dev/null | sort); do
  name=$(basename "$f" .weldlog)
  ex=$(grep -oE 'WELD_EXIT=[0-9-]+' "$f" | head -1)
  ncubes=$(grep -oE 'subgrid .*-> [0-9]+ cubes' "$f" | grep -oE '[0-9]+ cubes' | head -1)
  warn=$(grep -c 'WARNING:.*same_dir' "$f")
  echo "  $name  $ex  ($ncubes)  same_dir_warn=$warn"
done
