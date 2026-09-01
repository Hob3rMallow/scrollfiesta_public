#!/usr/bin/env bash
# Generate per-datatype "core" wrapper TUs for an in-tree TAUCS build.
# TAUCS compiles each core source once per datatype (TAUCS_CORE_GENERAL +
# TAUCS_CORE_DOUBLE for us). The configurator normally emits these; we do it by
# hand for the double-real-only subset. See configurator/taucs_structure.h.
set -eu
WRAP="$(cd "$(dirname "$0")" && pwd)/wrap"
mkdir -p "$WRAP"
rm -f "$WRAP"/*.c

gen() { # $1 = core source basename, $2 = D|GENERAL
  if [ "$2" = "D" ]; then def="TAUCS_CORE_DOUBLE"; sfx="_D"; else def="TAUCS_CORE_GENERAL"; sfx="_GEN"; fi
  printf '/* generated: %s, %s variant */\n#define %s\n#include "../../src/%s.c"\n' \
    "$1" "$2" "$def" "$1" > "$WRAP/${1}${sfx}.c"
}

# Files compiled in BOTH generic (dispatch) and double (typed) variants.
BOTH="taucs_complex taucs_ccs_base taucs_vec_base taucs_ccs_ops taucs_sn_llt taucs_ccs_factor_llt taucs_ccs_solve_llt"
# Generic-only.
GEN_ONLY="taucs_logging taucs_timer taucs_ccs_order taucs_linsolve taucs_malloc"
# Double-only.
D_ONLY="taucs_iter"

for f in $BOTH;     do gen "$f" D; gen "$f" GENERAL; done
for f in $GEN_ONLY; do gen "$f" GENERAL; done
for f in $D_ONLY;   do gen "$f" D; done

echo "Generated wrappers in $WRAP:"
ls -1 "$WRAP"
