#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

command -v cargo >/dev/null || { echo "cargo is not installed" >&2; exit 2; }
command -v rustc >/dev/null || { echo "rustc is not installed" >&2; exit 2; }
command -v python3 >/dev/null || { echo "python3 is not installed" >&2; exit 2; }
command -v rocminfo >/dev/null || { echo "rocminfo is not installed or ROCm is not active" >&2; exit 2; }
command -v nm >/dev/null || { echo "nm is not installed" >&2; exit 2; }

export HSA_ENABLE_DXG_DETECTION="${HSA_ENABLE_DXG_DETECTION:-1}"
export ROCPROFILER_REGISTER_ENABLED="${ROCPROFILER_REGISTER_ENABLED:-0}"
ROCMINFO_OUTPUT="$(rocminfo 2>/dev/null)"
if ! grep -q 'Name:[[:space:]]*gfx1201' <<<"$ROCMINFO_OUTPUT"; then
    echo "gfx1201 was not found in rocminfo output" >&2
    exit 3
fi

mkdir -p logs
REPORT="${MLS_RUN_REPORT:-$ROOT/logs/gfx1201-test-$(date -u +%Y%m%dT%H%M%SZ).log}"
exec > >(tee "$REPORT") 2>&1

export RUST_BACKTRACE=1
export MLS_BENCH_SIDE="${MLS_BENCH_SIDE:-127}"
export MLS_RADIUS="${MLS_RADIUS:-12}"

printf 'report=%s\n' "$REPORT"
printf 'utc_started=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf 'rustc=%s\n' "$(rustc --version)"
printf 'cargo=%s\n' "$(cargo --version)"
if command -v hipcc >/dev/null; then
    printf 'hipcc=%s\n' "$(hipcc --version 2>&1 | head -n 1)"
else
    printf 'hipcc=not-on-PATH\n'
fi
printf 'gpu=%s\n' "$(awk '/Name:[[:space:]]*gfx1201/{print $2; exit}' <<<"$ROCMINFO_OUTPUT")"
printf 'bench_side=%s\n' "$MLS_BENCH_SIDE"
printf 'radius=%s\n' "$MLS_RADIUS"
if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    printf 'git_commit=%s\n' "$(git rev-parse HEAD)"
fi

printf '\n== Independent Python numerical oracle ==\n'
python3 tools/numerical_oracle.py

printf '\n== Rust CPU reference tests ==\n'
cargo test --release --features rust-cpu

printf '\n== CubeCL CPU build and comparison ==\n'
cargo run --release --no-default-features \
    --features rust-cpu,cubecl-cpu --bin mls-bench -- cubecl-cpu

printf '\n== CubeCL HIP build and gfx1201 comparison ==\n'
cargo run --release --no-default-features \
    --features rust-cpu,cubecl-hip --bin mls-bench -- cubecl-hip

printf '\n== Shared-library ABI symbol ==\n'
cargo build --release --no-default-features --features rust-cpu,cubecl-hip
nm -D target/release/libherculaneum_mls_cubecl.so | grep ' T MLS_project_verts$'

printf '\nutc_finished=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf 'all_gates=PASS\n'
printf 'report=%s\n' "$REPORT"
