# CubeCL WGPU backend — status: KNOWN INCORRECT (do not use for production)

The `cubecl-wgpu` feature **compiles but produces numerically incorrect MLS output**
under both of CubeCL's wgpu code-generation paths. It is retained as an opt-in,
explicitly experimental backend for debugging and upstream regression testing only.
**It must not be used to generate production MLS geometry** — it silently returns
wrong vertices, which is more dangerous than failing to compile.

It is **never** chosen by `Backend::Auto` / `MLS_BACKEND=auto`; it is reachable only via
an explicit `MLS_BACKEND=wgpu` (or `Backend::CubeWgpu`) request.

## Support matrix

| Backend | Status |
| --- | --- |
| CubeCL HIP (AMD) | Hardware validated, correct |
| CubeCL CUDA (NVIDIA) | Hardware validated, correct (incl. direct pscamillo-CUDA parity) |
| CubeCL WGPU / SPIR-V (Vulkan) | Compiles; **produces incorrect output** |
| CubeCL WGPU / WGSL (DX12, GL) | Compiles; **produces incorrect output** |
| Rust/CPU | Correctness oracle |

## Root cause (identified 2026-07-27)

The kernel's 3×3 Jacobi eigensolver runs in **f64** (double precision) for eigen stability —
correct under CUDA/HIP, which support f64 natively. Under wgpu, cubecl generates an **f64
`atan2`** in the shader, which is **invalid SPIR-V**: `spirv-val` (on the naga-29.0.4-compiled
module) rejects it —

```
error: GLSL.std.450 Atan2: expected Result Type to be a 16 or 32-bit scalar or vector float type
  %770 = OpExtInst %double %1 Atan2 %769 %767
```

`GLSL.std.450 Atan2` is only defined for 16/32-bit floats, so emitting it on a `double` yields an
invalid module — which each naga backend (SPIR-V for Vulkan, HLSL for DX12/GL) then mishandles in
its own way. This single cause explains both failure modes, and why stripping the Kahan summation
changed nothing (the defect is downstream, in the f64 eigensolve).

Filed upstream: **[tracel-ai/cubecl#1446](https://github.com/tracel-ai/cubecl/issues/1446)**.

## Evidence (2026-07-27)

The identical single CubeCL kernel is correct under the CUDA and HIP runtimes but wrong
under wgpu, on **real hardware**, not just software adapters:

| Path | Adapter(s) | `mls-bench` max vertex error | Symptom |
| --- | --- | --- | --- |
| wgpu SPIR-V | lavapipe **and real AMD Vulkan** | 1.50 vox | output ≈ input (unprojected) |
| wgpu WGSL | Mesa d3d12-GL **and real AMD DX12** | 62.1 vox | garbage |

- The two wgpu paths fail **differently** (unprojected vs garbage) — the same invalid f64-`atan2`
  shader, mishandled per naga backend (SPIR-V vs HLSL).
- Failure is **independent of the environment/adapter**: it reproduces on mature AMD
  Vulkan and DX12 drivers, not only software fallbacks.
- Failure is **independent of the Kahan accumulation**: replacing `compensated_add` with a
  plain `*sum += value` produced byte-identical wrong output — confirming the defect is the
  downstream f64 eigensolve, not the accumulation.

Confirmed **CubeCL-WGPU code-generation defect**: it emits invalid SPIR-V for an f64 transcendental
(`spirv-val`-proven) instead of supporting it or erroring clearly. Reproducer + full detail filed
upstream: [tracel-ai/cubecl#1446](https://github.com/tracel-ai/cubecl/issues/1446).

## Do not
- Do not select wgpu via `Backend::Auto`.
- Do not advertise wgpu as a production backend.
- Do not casually change the shared kernel's eigensolve precision — see Revisit.

## Revisit
Two routes:
1. an upstream cubecl fix (or a clear compile error) for f64 transcendentals — tracking
   [tracel-ai/cubecl#1446](https://github.com/tracel-ai/cubecl/issues/1446); or
2. dropping the Jacobi eigensolver to **f32** on the wgpu path — likely to make wgpu compile and
   run, but it shifts the numerics, so it would require **re-benchmarking and re-validating the
   CUDA/HIP parity** before adoption.

The failing case is retained for both.
