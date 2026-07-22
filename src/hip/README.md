# HIP (AMD ROCm) MLS backend

A drop-in HIP/ROCm implementation of `MLS_project_verts` — the MLS-midpoint (LOP, mu=0)
projection stage — so ScrollFiesta per-cube meshing runs on AMD GPUs.

Origin: a clean-room implementation, written from a functional specification rather than
translated from `common/mls_project.c`. Validated standalone on a Radeon RX 9070
(gfx1201, RDNA4): **16.76x** MLS kernel, **5.81x** multi-cube throughput (16 concurrent),
topologically-equivalent weld-safe meshes vs the CPU path. Full evidence, independent
oracle, and conformance harness: https://github.com/altommo/scrollfiesta-mls-hip

## Contents
- `impl/`    original HIP kernel + a small C API, backed by an independent FP32 oracle
- `adapter/` bridges this repo`'s `MLS_project_verts(Arena_T, ...)` ABI to the kernel.
             The `Arena_T` is scratch-only for the CPU path and is ignored here; a clean
             spatial index is built per call. z,y,x order, `cell_origin`, and the
             `< MLS_MIN_NEIGHBOURS` -> original-position + zero-normal fallback all matched.

## Build
Requires ROCm with a gfx-capable `hipcc`. From `src/`:

    make cube_mesh_hip HIP_ARCH=gfx1201

This builds `cube_mesh` with the HIP MLS backend replacing `common/mls_project.o`.

## Status (honest scope)
The HIP objects compile with `hipcc` for gfx1201 and pass the standalone conformance,
real-cube parity, and end-to-end mesh-equivalence checks linked above. The full in-tree
`cube_mesh_hip` **link** additionally needs this repo`'s vendored deps
(Triangle / Clipper2 / PoissonRecon); that end-to-end link was not exercised in the
porting environment. FP32 GPU results differ from FP32 CPU by a few ULPs — the intended
acceptance bar is topological equivalence + weld safety, not bit-exact match.

## Attribution
Ports the algorithm and ABI of ScrollFiesta (MIT, (c) 2026 Nicholas Vining). The kernel is
original work; the adapter references the MIT ABI with file:line citations
(`adapter/ABI_MAPPING.md`).
