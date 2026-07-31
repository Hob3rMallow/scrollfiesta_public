#!/usr/bin/env bash
# grid_monitor.sh -- watch a running grid_pipeline per-cube run and continuously
# audit per-cube manifoldness as cubes are produced.
#
# Usage: bash scripts/grid_monitor.sh <output_dir> <percube_log>
#   <output_dir>   e.g. output/grid_4x21x21       (holds dump/<cube>/...)
#   <percube_log>  e.g. output/grid_4x21x21_perCube.log
#
# Every INTERVAL seconds it logs: cubes done, failed (exit!=0), how many meshes
# manifold-checked, and the running count of NON-MANIFOLD/foldover cubes. New
# (not-yet-checked) meshes are diffed and only those are re-checked, so the sweep
# stays cheap. Any bad cube's manifold_check line is appended to <out>_bad.txt.
# Exits when the per-cube log records PERCUBE_EXIT=, then writes a full final
# manifold report to <out>_manifold_final.txt.
set -u
OUT="${1:?output_dir}"
LOG="${2:?percube_log}"
INTERVAL="${3:-300}"
MC="build/Release/manifold_check.exe"
GLOB_SUFFIX="_step12_final/*_step12_final_all.obj"

MON="${OUT}_monitor.log"
CHECKED="${OUT}_checked.txt"
BAD="${OUT}_bad.txt"
: > "$CHECKED"; : > "$BAD"
echo "monitor start $(date)  out=$OUT log=$LOG interval=${INTERVAL}s" > "$MON"

list_meshes() { ls -d "$OUT"/dump/*/*${GLOB_SUFFIX} 2>/dev/null | sort; }

while true; do
  done=$(grep -c 'exit=' "$LOG" 2>/dev/null); done=${done:-0}
  fail=$(grep 'exit=' "$LOG" 2>/dev/null | grep -vc 'exit=0'); fail=${fail:-0}
  list_meshes > /tmp/gm_cur.txt
  comm -13 "$CHECKED" /tmp/gm_cur.txt > /tmp/gm_new.txt
  nnew=$(wc -l < /tmp/gm_new.txt | tr -d ' ')
  if [ "$nnew" -gt 0 ]; then
    "$MC" --list /tmp/gm_new.txt 2>/dev/null \
      | grep -E 'nm_edge=[1-9][0-9]*|nm_vert=[1-9][0-9]*|same_dir=[1-9][0-9]*' >> "$BAD" || true
    cp /tmp/gm_cur.txt "$CHECKED"
  fi
  nchecked=$(wc -l < "$CHECKED" | tr -d ' ')
  nbad=$(wc -l < "$BAD" | tr -d ' ')
  echo "$(date +%H:%M:%S)  done=$done  fail=$fail  checked=$nchecked  bad=$nbad" >> "$MON"
  if grep -q 'PERCUBE_EXIT=' "$LOG" 2>/dev/null; then
    echo "monitor: per-cube run finished -> final sweep" >> "$MON"
    break
  fi
  sleep "$INTERVAL"
done

list_meshes > /tmp/gm_final.txt
"$MC" --list /tmp/gm_final.txt > "${OUT}_manifold_final.txt" 2>&1
echo "monitor done $(date)" >> "$MON"
