#!/usr/bin/env bash
# color_developability.sh <cube_dump_dir> <out_dir> [heatmax]
#
# Bakes DEVELOPABILITY into the per-step output OBJs as vertex colour: for each
# pipeline step's whole-mesh output (<cube>_<stage>_all.obj) it writes a
# <cube>_<stage>_dev.obj whose per-vertex RGB encodes |Gaussian curvature|
# (angle defect) -- blue = developable (K~0, what a papyrus wrap should be),
# warm/red = non-developable (cone/saddle = a bridge/tangle/fold artifact).
# Open the *_dev.obj in MeshLab/Blender (vertex colours on). Point-cloud stages
# (step0/step3 _mls -- no faces) are skipped. Prints per-step S|K| / spike count
# so the stage where non-developability is injected is legible at a glance.
#
# Bash-as-runner only: developability.exe does the measurement.
set -u
DUMP="${1:?cube_dump_dir e.g. output/.../dump/z04352_y02048_x01920}"
OUT="${2:?out_dir}"
HEATMAX="${3:-0.5}"               # |K| (rad) mapped to full red; = the spike thresh
DEV=./build/Release/developability.exe
c="$(basename "$DUMP")"
mkdir -p "$OUT"
echo "cube $c -> $OUT   (colour scale: blue 0 .. red >= ${HEATMAX} rad)"
printf "%-18s %10s %8s %8s\n" stage "S|K|" "spikes" "max|K|"
shopt -s nullglob
for sd in "$DUMP/${c}_"*/; do
  st="$(basename "$sd")"; st="${st#"${c}_"}"
  f="$sd/${c}_${st}_all.obj"
  [ -f "$f" ] || continue
  grep -q '^f ' "$f" || continue            # skip faceless (point-cloud) stages
  out="$OUT/${c}_${st}_dev.obj"
  rep="$("$DEV" "$f" --thresh "$HEATMAX" --heatmax "$HEATMAX" --heatmap "$out" 2>/dev/null)"
  sk=$(printf '%s\n' "$rep" | sed -n 's/.*curvature) = \([0-9.]*\) rad.*/\1/p')
  sp=$(printf '%s\n' "$rep" | sed -n 's/.*over thresh [0-9.]* rad: \([0-9]*\) .*/\1/p')
  mk=$(printf '%s\n' "$rep" | sed -n 's/.*max|K| = \([0-9.]*\) rad.*/\1/p')
  printf "%-18s %10s %8s %8s\n" "$st" "${sk:-?}" "${sp:-?}" "${mk:-?}"
done
echo "done -> $OUT"
