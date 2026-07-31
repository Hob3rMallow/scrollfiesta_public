#!/usr/bin/env bash
# render_stage_trace.sh <cube_dump_dir> <out_dir> [size]
#
# Renders every pipeline stage of one cube dump into a stage-ordered gallery PDF
# so an artifact can be flagged by eye at the stage it first appears.
# For each stage: one whole-cube overview (_all.obj, prefix 'A') + every
# per-component OBJ. PNG filenames are stage-ordered (s00..s12) so the gallery,
# which sorts by filename, groups them by stage.
#
# Runner only (mesh_render + make_gallery_pdf do the real work). Bash-as-runner
# is allowed by CLAUDE.md; the analysis lives in the C tool.
set -u
DUMP="${1:?usage: render_stage_trace.sh <cube_dump_dir> <out_dir> [size]}"
OUT="${2:?need out dir}"
SIZE="${3:-900}"
R=./build/Release/mesh_render.exe
PNG="$OUT/png"
mkdir -p "$PNG"

stages=(step0_mls step1_bpa step2_cc step3_mls step4_bpa step5_cc \
        step6_cc_mls step7_cc_bpa step8_holefill step9_sever \
        step10_qem step11_kibble step12_final)

cubeid="$(basename "$DUMP")"
MAXJ=24
launched=0

render() {  # <obj> <pngname>
  "$R" "$1" "$PNG/$2" --size "$SIZE" >/dev/null 2>&1 &
  launched=$((launched+1))
  if [ $((launched % MAXJ)) -eq 0 ]; then wait; fi
}

i=0
for st in "${stages[@]}"; do
  sd="$DUMP/${cubeid}_${st}"
  pfx="$(printf 's%02d' "$i")"
  i=$((i+1))
  [ -d "$sd" ] || { echo "skip (no dir): $st"; continue; }
  # whole-cube overview first ('A' sorts before component suffixes)
  allobj="$sd/${cubeid}_${st}_all.obj"
  [ -f "$allobj" ] && render "$allobj" "${pfx}_${st}_A_ALL.png"
  # every per-piece / per-component OBJ (numbered NNN, or in01_p00 / c00_p03)
  for f in "$sd/${cubeid}_${st}_"*.obj; do
    [ -e "$f" ] || continue
    base="$(basename "$f" .obj)"
    case "$base" in *_all) continue;; esac
    suf="${base#"${cubeid}_${st}_"}"   # 001 | in01_p00 | c00_p03
    nf="$(grep -c '^f ' "$f")"
    render "$f" "${pfx}_${st}_${suf}_${nf}f.png"
  done
done
wait
echo "rendered $launched PNGs -> $PNG"

python scripts/make_gallery_pdf.py "$PNG" "$OUT/stage_trace.pdf" --cols 5 --rows 4 --tile 380
echo "gallery -> $OUT/stage_trace.pdf"
