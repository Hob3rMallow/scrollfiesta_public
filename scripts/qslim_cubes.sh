#!/usr/bin/env bash
# qslim_cubes.sh <grid_output_dir> [keep_ratio=0.10] [maxjobs=32]
#
# Decimate each cube's welded-stage final (step12_final_all.obj) to keep_ratio of
# its faces via qslim_obj, written to <grid_output_dir>/decimated/<cube>_decim.obj.
# For working with lightweight individual cubes before any weld.
#
# Runner only -- the QEM lives in the C tool qslim_obj (which has --selftest).
set -u
OUT="${1:?usage: qslim_cubes.sh <grid_output_dir> [keep_ratio=0.10] [maxjobs=32]}"
RATIO="${2:-0.10}"
MAXJ="${3:-32}"
Q=./build/Release/qslim_obj.exe
DUMP="$OUT/dump"
DEC="$OUT/decimated"; mkdir -p "$DEC"
LOG="$OUT/qslim_cubes.log"
: > "$LOG"

[ -x "$Q" ] || { echo "ERROR: $Q not found"; exit 1; }
[ -d "$DUMP" ] || { echo "ERROR: no dump dir $DUMP"; exit 1; }

launched=0
for d in "$DUMP"/z*/; do
  [ -d "$d" ] || continue
  cube=$(basename "$d")
  src="$d/${cube}_step12_final/${cube}_step12_final_all.obj"
  [ -f "$src" ] || { echo "$cube SKIP no-final" >> "$LOG"; continue; }
  {
    fin=$(grep -c '^f ' "$src")
    out="$DEC/${cube}_decim.obj"
    if "$Q" "$src" "$out" "$RATIO" >/dev/null 2>&1 && [ -f "$out" ]; then
      echo "$cube $fin $(grep -c '^f ' "$out")" >> "$LOG"
    else
      echo "$cube FAIL $fin -" >> "$LOG"
    fi
  } &
  launched=$((launched+1))
  [ $((launched % MAXJ)) -eq 0 ] && { wait; echo "  ...$launched submitted"; }
done
wait

ok=$(awk '$2 ~ /^[0-9]+$/' "$LOG" | wc -l)
fail=$(grep -c ' FAIL ' "$LOG")
skip=$(grep -c ' SKIP ' "$LOG")
tin=$(awk '$2 ~ /^[0-9]+$/ {s+=$2} END{print s+0}' "$LOG")
tout=$(awk '$3 ~ /^[0-9]+$/ {s+=$3} END{print s+0}' "$LOG")
pct=$(awk "BEGIN{if($tin>0)printf \"%.1f\",100.0*$tout/$tin; else print 0}")
echo "=== qslim_cubes done: $ok decimated, $fail fail, $skip no-final | ratio=$RATIO ==="
echo "    total faces: $tin -> $tout  (${pct}% kept)"
echo "    output: $DEC  (log: $LOG)"
