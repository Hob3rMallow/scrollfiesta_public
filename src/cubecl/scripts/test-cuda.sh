#!/usr/bin/env bash
set -euo pipefail

# Self-contained CubeCL CUDA validation: builds the CUDA backend and checks it
# against the Rust CPU oracle on the synthetic corrugated-plane benchmark. This
# is the NVIDIA analogue of test-gfx1201.sh and needs only an NVIDIA GPU with a
# working CUDA toolkit; it does not require any ScrollFiesta capture.
# The real-capture cross-check against pscamillo's CUDA kernel is a separate step.

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

command -v cargo >/dev/null || { echo "cargo is not installed" >&2; exit 2; }
command -v rustc >/dev/null || { echo "rustc is not installed" >&2; exit 2; }
command -v python3 >/dev/null || { echo "python3 is not installed" >&2; exit 2; }
command -v nvidia-smi >/dev/null || { echo "nvidia-smi is not installed or no NVIDIA driver is active" >&2; exit 2; }
command -v nvcc >/dev/null || { echo "nvcc is not installed or the CUDA toolkit is not on PATH" >&2; exit 2; }
command -v nm >/dev/null || { echo "nm is not installed" >&2; exit 2; }

mkdir -p logs
REPORT="${MLS_RUN_REPORT:-$ROOT/logs/cuda-test-$(date -u +%Y%m%dT%H%M%SZ).log}"
exec > >(tee "$REPORT") 2>&1

export RUST_BACKTRACE=1
export MLS_BENCH_SIDE="${MLS_BENCH_SIDE:-127}"
export MLS_RADIUS="${MLS_RADIUS:-12}"

printf 'report=%s\n' "$REPORT"
printf 'utc_started=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf 'rustc=%s\n' "$(rustc --version)"
printf 'cargo=%s\n' "$(cargo --version)"
printf 'nvcc=%s\n' "$(nvcc --version | grep -i release || nvcc --version | tail -n 1)"
printf 'gpu=%s\n' "$(nvidia-smi --query-gpu=name --format=csv,noheader | head -n 1)"
printf 'bench_side=%s\n' "$MLS_BENCH_SIDE"
printf 'radius=%s\n' "$MLS_RADIUS"
if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    printf 'git_commit=%s\n' "$(git rev-parse HEAD)"
fi

printf '\n== Independent Python numerical oracle ==\n'
python3 tools/numerical_oracle.py

printf '\n== Rust CPU reference tests ==\n'
cargo test --release --features rust-cpu

printf '\n== CubeCL CUDA build and comparison vs Rust CPU oracle ==\n'
cargo run --release --no-default-features \
    --features rust-cpu,cubecl-cuda --bin mls-bench -- cubecl-cuda

printf '\n== Shared-library ABI symbol ==\n'
cargo build --release --no-default-features --features rust-cpu,cubecl-cuda
nm -D target/release/libherculaneum_mls_cubecl.so | grep ' T MLS_project_verts$'

printf '\nutc_finished=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf 'all_gates=PASS\n'
printf 'report=%s\n' "$REPORT"
