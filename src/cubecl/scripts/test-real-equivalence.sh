#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export HSA_ENABLE_DXG_DETECTION="${HSA_ENABLE_DXG_DETECTION:-1}"
export ROCPROFILER_REGISTER_ENABLED="${ROCPROFILER_REGISTER_ENABLED:-0}"
CLEANROOM_ROOT="${MLS_CLEANROOM_ROOT:-$HOME/mls-hip-cleanroom}"
PARITY_ROOT="$CLEANROOM_ROOT/parity-real"
RUNNER="${MLS_CLEANROOM_RUNNER:-$PARITY_ROOT/run_real_parity}"
INPUT_OBJ="${MLS_REAL_INPUT_OBJ:-$PARITY_ROOT/cpu_diag/dump/z16128_y02560_x07680/z16128_y02560_x07680_step0_mls/z16128_y02560_x07680_step0_mls_pre000.obj}"
PASSES="${MLS_REAL_PASSES:-5}"
RADIUS="${MLS_RADIUS:-12}"
ORIGIN_Z="${MLS_ORIGIN_Z:-16128}"
ORIGIN_Y="${MLS_ORIGIN_Y:-2560}"
ORIGIN_X="${MLS_ORIGIN_X:-7680}"
REPS="${MLS_BENCH_REPS:-3}"

[[ -x "$RUNNER" ]] || { echo "clean-room runner is missing: $RUNNER" >&2; exit 2; }
[[ -f "$INPUT_OBJ" ]] || { echo "real ScrollFiesta capture is missing: $INPUT_OBJ" >&2; exit 2; }
command -v cargo >/dev/null || { echo "cargo is not installed" >&2; exit 2; }
command -v rocminfo >/dev/null || { echo "rocminfo is not installed" >&2; exit 2; }

ROCMINFO_OUTPUT="$(rocminfo 2>/dev/null)"
grep -q 'Name:[[:space:]]*gfx1201' <<<"$ROCMINFO_OUTPUT" || {
    echo "gfx1201 was not found in rocminfo output" >&2
    exit 3
}

mkdir -p "$ROOT/logs"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
PREFIX="$ROOT/logs/real-equivalence-$STAMP"
REPORT="$PREFIX.log"
exec > >(tee "$REPORT") 2>&1

printf 'utc_started=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf 'input=%s\npasses=%s\nradius=%s\norigin_zyx=%s,%s,%s\nrepetitions=%s\n' \
    "$INPUT_OBJ" "$PASSES" "$RADIUS" "$ORIGIN_Z" "$ORIGIN_Y" "$ORIGIN_X" "$REPS"

printf '\n== Original ScrollFiesta CPU and clean-room HIP C++ ==\n'
set +e
MLS_BENCH_REPS="$REPS" "$RUNNER" \
    "$INPUT_OBJ" UPSTREAM_CPU "$PREFIX" "$RADIUS" real_comp001 \
    "$PASSES" "$ORIGIN_Z" "$ORIGIN_Y" "$ORIGIN_X" \
    | tee "$PREFIX.cleanroom.txt"
CLEANROOM_STATUS=${PIPESTATUS[0]}
set -e
if [[ "$CLEANROOM_STATUS" -gt 1 ]]; then
    echo "clean-room runner failed with status $CLEANROOM_STATUS" >&2
    exit "$CLEANROOM_STATUS"
fi

printf '\n== CubeCL HIP on the identical local z-y-x capture ==\n'
cd "$ROOT"
MLS_BENCH_REPS="$REPS" cargo run --release --no-default-features \
    --features cubecl-hip --bin mls-real-bench -- \
    cubecl-hip "$PREFIX.pre.positions_zyx.f32" "$PREFIX" \
    "$PASSES" "$RADIUS" "$ORIGIN_Z" "$ORIGIN_Y" "$ORIGIN_X" \
    "$PREFIX.cpu.positions_zyx.f32" "$PREFIX.cpu.normals_zyx.f32" \
    "$PREFIX.hip.positions_zyx.f32" "$PREFIX.hip.normals_zyx.f32"

printf '\nutc_finished=%s\nreport=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$REPORT"
