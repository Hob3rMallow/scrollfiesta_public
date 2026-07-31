#!/usr/bin/env bash
# bpa_tear_census.sh — per-cube BPA tear census for one grid-run output dir.
# Emits TSV: cube_id  boundary_edges  tears(small_gaps_other)  pinch_steps  F  V
# for a chosen stage's combined OBJ of every cube under <run_dir>.
#
# Usage: bpa_tear_census.sh <run_dir> [stage=step1_bpa] [seam_gap.exe]
#   stage: step1_bpa (raw BPA) | step12_final (what ships to weld) | ...
#
# (Originally rho_ab_census.sh; recreated 2026-07-09 after the adaptive-rho
# revert deleted it. Used for the patch_repair A/B — see changelog.)
set -u
RUN="${1:?usage: bpa_tear_census.sh <run_dir> [stage] [seam_gap]}"
STAGE="${2:-step1_bpa}"
SG="${3:-build/Release/seam_gap.exe}"

jget() { echo "$1" | grep -oE "\"$2\":[0-9-]+" | head -1 | cut -d: -f2; }

# the combined OBJ is dumped twice (per-stage dir AND <cube>_all_obj/); take
# the per-stage copy only, else every cube is double-counted.
find "$RUN" -path "*_${STAGE}/*" -iname "*_${STAGE}_all.obj" | sort | while read -r f; do
  cube=$(basename "$f" | sed "s/_${STAGE}_all\.obj//")
  j=$("$SG" "$f" --band 6 2>/dev/null | grep SEAMGAP_JSON)
  [ -z "$j" ] && continue
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$cube" \
    "$(jget "$j" boundary_edges)" "$(jget "$j" small_gaps_other)" \
    "$(jget "$j" pinch_steps)" "$(jget "$j" F)" "$(jget "$j" V)"
done
