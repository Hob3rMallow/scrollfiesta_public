#!/usr/bin/env bash
# split_good_bad.sh <gate_report.raw> <out_dir>
#
# Partition the developability-gate verdicts into:
#   <out>/good_welded.obj -- every gate-PASS sheet, welded across cube seams
#                            (grid_weld over a good-only step12_final dump).
#   <out>/bad.obj         -- every gate-FAIL sheet, plain-merged (world coords,
#                            face-renumbered, NOT welded).
# Bash-as-runner only: obj_merge + grid_weld do the work.
set -u
RAW="${1:?gate report .raw}"; OUT="${2:?out dir}"
M=./build/Release/obj_merge.exe
W=./build/Release/grid_weld.exe
GD="$OUT/good_dump"
mkdir -p "$OUT"; rm -rf "$GD"; mkdir -p "$GD"

# verdict + path, one per sheet
awk '/^GATE/{b=0; for(i=2;i<=NF;i++) if($i=="bad=1") b=1; print b, $NF}' "$RAW" > "$OUT/verdicts.txt"
ngood=$(awk '$1==0' "$OUT/verdicts.txt" | wc -l)
nbad=$(awk  '$1==1' "$OUT/verdicts.txt" | wc -l)
echo "sheets: good=$ngood  bad=$nbad"

# --- bad.obj: plain world-coord merge (not welded) ---
awk '$1==1{print $2}' "$OUT/verdicts.txt" | xargs "$M" "$OUT/bad.obj"

# --- good: build a good-only step12_final dump, then weld it ---
awk '$1==0{print $2}' "$OUT/verdicts.txt" | while read -r p; do
  cube=$(basename "$(dirname "$(dirname "$p")")")
  printf '%s %s\n' "$cube" "$p"
done | sort > "$OUT/good_by_cube.txt"

for cube in $(awk '{print $1}' "$OUT/good_by_cube.txt" | sort -u); do
  sd="$GD/$cube/${cube}_step12_final"; mkdir -p "$sd"
  awk -v c="$cube" '$1==c{print $2}' "$OUT/good_by_cube.txt" \
    | xargs "$M" "$sd/${cube}_step12_final_all.obj" 2>/dev/null
done
echo "good cubes in weld dump: $(ls "$GD" | wc -l) -> welding $(date +%H:%M:%S)"
"$W" "$GD" "$OUT/good_welded.obj" --stage step12_final > "$OUT/good_weld.log" 2>&1

echo "=== results ==="
printf "bad.obj          %8s v  %8s f  (51 sheets, plain merge)\n" \
       "$(grep -c '^v ' "$OUT/bad.obj")" "$(grep -c '^f ' "$OUT/bad.obj")"
printf "good_welded.obj  %8s v  %8s f  (welded)\n" \
       "$(grep -c '^v ' "$OUT/good_welded.obj")" "$(grep -c '^f ' "$OUT/good_welded.obj")"
grep -E '"cubes_processed"|"unpaired"|"non_manifold"|"same_dir' "$OUT/good_welded.obj.weld_report.json" 2>/dev/null
