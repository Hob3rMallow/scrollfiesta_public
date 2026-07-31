#!/usr/bin/env bash
# grid_resume_loop.sh -- finish a partially-complete per-cube grid by re-running
# grid_pipeline --skip-existing until it CONVERGES (a pass adds no new complete
# step12_final meshes). Each pass skips already-done cubes and retries the rest,
# so transient single-cube failures get another shot; genuinely-empty edge cubes
# never complete and stop the loop. NEVER welds (--skip-weld always).
#
# Usage: bash scripts/grid_resume_loop.sh <grid_dir> <output_dir> <log> [maxpass]
set -u
GRID="${1:?grid_dir}"
OUT="${2:?output_dir}"
LOG="${3:?log}"
MAXP="${4:-5}"
EXE=build/Release/grid_pipeline.exe
TOTAL=1753

count_done() { ls -d "$OUT"/dump/*/*_step12_final/*_step12_final_all.obj 2>/dev/null | wc -l | tr -d ' '; }

echo "resume-loop start $(date)  grid=$GRID out=$OUT maxpass=$MAXP total=$TOTAL" > "$LOG"
for p in $(seq 1 "$MAXP"); do
  before=$(count_done)
  echo "=== PASS $p START $(date +%H:%M:%S)  complete_before=$before/$TOTAL ===" | tee -a "$LOG"
  "$EXE" "$GRID" "$OUT" --halo 13 --threads-per-cube 1 --max-concurrent 112 \
      --skip-existing --skip-weld --qem >> "$LOG" 2>&1
  after=$(count_done)
  echo "=== PASS $p DONE  $(date +%H:%M:%S)  complete_after=$after/$TOTAL  (+$((after-before)) new) ===" | tee -a "$LOG"
  if [ "$after" -le "$before" ]; then
    echo "CONVERGED: pass $p added 0 new complete meshes; remaining $((TOTAL-after)) cubes are un-completable (empty/edge). Stopping." | tee -a "$LOG"
    break
  fi
done
echo "resume-loop FINISHED $(date)  final_complete=$(count_done)/$TOTAL  (NO WELD by request)" | tee -a "$LOG"
