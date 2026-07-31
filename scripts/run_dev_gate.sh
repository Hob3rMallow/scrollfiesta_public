#!/usr/bin/env bash
# run_dev_gate.sh <dump_dir> <out_report> [-- extra dev_gate_tool args...]
#
# Runs the end-of-pipeline developability gate (dev_gate_tool == the same
# DevGate_classify the pipeline calls) over EVERY step12_final per-component
# sheet in a dump dir, and aggregates the verdicts. Because the gate is a pure
# function of the final mesh, this equals running it in-pipeline -- but in
# minutes, with no re-mesh. Bash-as-runner only.
set -u
DUMP="${1:?dump_dir}"; REP="${2:?out report path}"; shift 2
EXTRA=(); [ "${1:-}" = "--" ] && { shift; EXTRA=("$@"); }
T=./build/Release/dev_gate_tool.exe
mkdir -p "$(dirname "$REP")"
RAW="$REP.raw"; : > "$RAW"

mapfile -t comps < <(ls "$DUMP"/*/*_step12_final/*_step12_final_[0-9]*.obj 2>/dev/null)
echo "gating ${#comps[@]} step12_final sheets from $DUMP  $(date +%H:%M:%S)"
printf '%s\n' "${comps[@]}" | xargs -P 16 -I {} "$T" {} "${EXTRA[@]}" 2>/dev/null >> "$RAW"

awk '
/^GATE/{
  bad=0; reason="none"; nv=nf=iv=ib=orc=0; frac=0; mk=0; path=$NF
  for(i=2;i<=NF;i++){ n=index($i,"="); if(!n) continue; k=substr($i,1,n-1); v=substr($i,n+1)
    if(k=="bad")bad=v; else if(k=="reason")reason=v; else if(k=="nv")nv=v; else if(k=="nf")nf=v;
    else if(k=="int")iv=v; else if(k=="int_bad")ib=v; else if(k=="frac")frac=v; else if(k=="maxK")mk=v; else if(k=="oracle")orc=v }
  total++; rc[reason]++
  if(bad+0==1){ badc++
    s=path; sub(/.*\//,"",s); sub(/_step12_final/,"",s); sub(/\.obj$/,"",s)
    printf "%-34s nf=%-6s int=%-6s int_bad=%-5s frac=%-7s maxK=%-6s oracle=%-2s [%s]\n",
           s, nf, iv, ib, frac, mk, orc, reason >> BADF
  }
}
END{
  printf "components=%d  BAD=%d (%.1f%%)\n", total, badc, total?100.0*badc/total:0
  print "by reason:"
  for(r in rc) printf "  %-15s %d\n", r, rc[r]
}' BADF="$REP.bad" "$RAW" | tee "$REP"

echo "" | tee -a "$REP"
echo "=== BAD sheets (sorted by maxK desc) ===" | tee -a "$REP"
sort -t= -k6 -rn "$REP.bad" 2>/dev/null | tee -a "$REP"
echo "report -> $REP  ($(grep -c . "$REP.bad" 2>/dev/null) bad sheets listed)"
