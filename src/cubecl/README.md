# ScrollFiesta MLS → CubeCL (portable GPU)

Clean-room MLS midpoint projection with safe Rust, CubeCL CPU, and CubeCL HIP
backends. The CubeCL HIP path is validated on an RX 9070 (`gfx1201`) against
the original ScrollFiesta CPU operator and the clean-room HIP C++ kernel.

> **Related — native HIP/ROCm variant:** [`scrollfiesta-mls-hip`](https://github.com/altommo/scrollfiesta-mls-hip)
> is the C++/HIP port of the same MLS kernel — ~6–7% faster on AMD and the parity
> reference this backend is validated against (upstream PR
> [Hob3rMallow/scrollfiesta_public#4](https://github.com/Hob3rMallow/scrollfiesta_public/pull/4)).
> This CubeCL variant trades that margin for a single Rust source that targets AMD, NVIDIA, and WGPU.

## ScrollFiesta contract

The exact exported symbol is:

```c
void MLS_project_verts(
    void *arena,
    const float *verts,
    size_t nv,
    float radius_vox,
    const float *cell_origin,
    float *out_verts,
    float *out_normals);
```

It now implements the production contract:

- exactly one projection per call;
- `verts` is both query and current support, in z-y-x order;
- the counting-sort grid is rebuilt from that support on every call;
- ScrollFiesta retains its caller-owned 5-pass extract and 20-pass resplit loops;
- ScrollFiesta's `Arena_T` is scratch state and is intentionally not dereferenced;
- input/output aliasing and a null normals output are supported.

Set `MLS_BACKEND=rust-cpu`, `cubecl-cpu`, or `cubecl-hip` before calling the
legacy symbol. `include/scrollfiesta_arena_adapter.hpp` documents the link-only
integration boundary; no private arena layout mapping is required.

The Rust-owned arena constructors remain available for standalone xyz/static-
support work. Use `MLS_cubecl_project_with_arena(...)` with those handles;
never pass a foreign arena to that explicit API.

## Numerical implementation

Each production call:

1. converts z-y-x buffers to the shared kernel's xyz grid convention;
2. builds a deterministic counting-sort grid with ScrollFiesta-compatible cell
   phase and within-cell z-y-x ordering;
3. gathers support inside the Wendland C2 radius;
4. uses compensated FP32 centroid and covariance accumulation;
5. solves the 3x3 covariance in ScrollFiesta's native basis with a scalar FP64
   Jacobi eigensolver;
6. applies the tangent-plane projection and returns one updated cloud.

The scalar eigensolver is deliberate. CubeCL's dynamic local-array form passed
small synthetic tests but selected tangent directions on the 331,013-point real
capture. Scalar state fixed that backend-lowering failure and is covered by the
real-data gate.

## Build and test

```bash
./scripts/install-rust-wsl.sh
cargo test --release --features rust-cpu
cargo test --release --no-default-features --features rust-cpu,cubecl-cpu
./scripts/test-gfx1201.sh
```

Run the identical-data production comparison when the clean-room fixtures are
available at `~/mls-hip-cleanroom`:

```bash
MLS_REAL_PASSES=5 MLS_BENCH_REPS=3 ./scripts/test-real-equivalence.sh
MLS_REAL_PASSES=20 MLS_BENCH_REPS=1 ./scripts/test-real-equivalence.sh
```

The script runs original ScrollFiesta CPU, clean-room HIP C++, and CubeCL HIP on
the same captured z-y-x array. Timings include every per-pass grid rebuild and
host/device transfer. Override `MLS_CLEANROOM_ROOT` or `MLS_REAL_INPUT_OBJ` for
another checkout/capture.

## Result

On the 331,013-vertex Component 1 capture, radius 12:

| Passes | Backend | Warm/median time | Max position vs CPU | RMS position vs CPU | Max normal vs CPU |
| ---: | --- | ---: | ---: | ---: | ---: |
| 5 | ScrollFiesta CPU | 11,014.5 ms | reference | reference | reference |
| 5 | clean HIP C++ | 869.455 ms | 0.00213777 | 0.0000102769 | 0.0248666 deg |
| 5 | CubeCL HIP | 927.511 ms | 0.00205225 | 0.000010215 | 0.0246873 deg |
| 20 | ScrollFiesta CPU | 51,757.0 ms | reference | reference | reference |
| 20 | clean HIP C++ | 4,956.23 ms | 0.0390429 | 0.000104177 | 1.69153 deg |
| 20 | CubeCL HIP | 5,177.33 ms | 0.0398853 | 0.000104474 | 1.69095 deg |

Both CubeCL runs pass the strict `<0.25` voxel weld-safety boundary and have no
zero-normal mismatches. Like clean HIP, Component 1 does not pass the strict
per-vertex `0.00124` voxel / `0.006` degree gate at a small set of ambiguous
neighborhoods. This is not relabelled as strict parity.

CubeCL HIP is about 12x faster than original CPU for five passes and 10x for 20
passes. It is currently 6.7% slower than clean HIP at five passes and 4.5%
slower at 20 passes, so the shared CubeCL implementation is credible but does
not beat the handwritten kernel yet.

See `VALIDATION.md` and the two tracked reports under `logs/` for raw evidence.

## Files

```text
src/oracle.rs                    safe Rust one-pass reference
src/cube_backend.rs              shared CubeCL CPU/HIP one-pass kernel
src/grid.rs                      deterministic counting-sort grid
src/lib.rs                       production C ABI and explicit arena API
src/bin/mls_real_bench.rs        real moving-support runner/comparator
include/mls_cubecl.h             public C header
include/scrollfiesta_arena_adapter.hpp
scripts/test-gfx1201.sh          synthetic backend hardware gate
scripts/test-real-equivalence.sh identical real-data three-way gate
```
