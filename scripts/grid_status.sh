#!/usr/bin/env bash
# grid_status.sh -- race-free status snapshot for a running (or finished)
# grid_pipeline per-cube run: cumulative ETA + a single full manifold_check
# sweep over every step12_final mesh currently on disk.
#
# Usage: bash scripts/grid_status.sh <output_dir> <percube_log> [total_cubes]
set -u
OUT="${1:?output_dir}"
LOG="${2:?percube_log}"
TOTAL="${3:-1753}"
MC="build/Release/manifold_check.exe"

done=$(grep -c 'exit=' "$LOG" 2>/dev/null); done=${done:-0}
fail=$(grep 'exit=' "$LOG" 2>/dev/null | grep -vc 'exit=0'); fail=${fail:-0}
start_epoch=$(date -d "$(head -1 "$LOG" | sed 's/start //; s/ max_concurrent=.*//')" +%s 2>/dev/null)
now=$(date +%s)
rem=$((TOTAL-done))

echo "=== grid status $(date +%H:%M:%S) ==="
if [ -n "${start_epoch:-}" ]; then
  el=$((now-start_epoch))
  awk -v d=$done -v tot=$TOTAL -v f=$fail -v e=$el -v rem=$rem 'BEGIN{
    printf "done=%d/%d  fail=%d  elapsed=%.0f min\n", d, tot, f, e/60;
    if(d>0){r=d/e; printf "cumulative rate=%.1f cubes/min  ETA_remaining~%.0f min (%.2f h)\n", r*60, rem/r/60, rem/r/3600}
  }'
else
  echo "done=$done/$TOTAL  fail=$fail  (could not parse start time)"
fi

# Failed cubes + reason (distinguish expected 'no meshes' empties from real crashes)
if [ "$fail" -gt 0 ]; then
  echo "--- failed cubes (reason) ---"
  grep 'exit=' "$LOG" | grep -v 'exit=0' | grep -oE 'z[0-9]+_y[0-9]+_x[0-9]+' | while read -r c; do
    reason=$(grep -hoE 'produced no meshes|Status: FAILED|assert|Segmentation|abort' "$OUT/logs/$c.log" 2>/dev/null | head -1)
    et=$(grep -hoE 'extract=[0-9.]+' "$OUT/logs/$c.log" 2>/dev/null | head -1)
    echo "  $c  ${reason:-?}  ($et)"
  done
fi

# Race-free manifold sweep over all meshes currently on disk
ls -d "$OUT"/dump/*/*_step12_final/*_step12_final_all.obj 2>/dev/null | sort > /tmp/gs_mesh.txt
nmesh=$(wc -l < /tmp/gs_mesh.txt | tr -d ' ')
echo "--- manifold sweep ($nmesh meshes) ---"
if [ "$nmesh" -gt 0 ]; then
  "$MC" --list /tmp/gs_mesh.txt 2>/dev/null > /tmp/gs_sweep.txt
  nnm=$(grep -cE 'VERDICT=NONMANIFOLD' /tmp/gs_sweep.txt)
  nsd=$(grep -E 'same_dir=[1-9]' /tmp/gs_sweep.txt | grep -c 'V=')
  echo "  2-manifold(edge+vert): $((nmesh-nnm))/$nmesh   NON-manifold: $nnm   same_dir(foldover): $nsd"
  if [ "$nnm" -gt 0 ]; then echo "  NON-MANIFOLD cubes:"; grep 'VERDICT=NONMANIFOLD' /tmp/gs_sweep.txt | grep -oE 'z[0-9]+_y[0-9]+_x[0-9]+' | sed 's/^/    /'; fi
  if [ "$nsd" -gt 0 ]; then echo "  same_dir cubes:"; grep -E 'same_dir=[1-9]' /tmp/gs_sweep.txt | grep 'V=' | sed -E 's/_step12.*same_dir=/  same_dir=/; s/^/    /'; fi
fi
