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

## Evidence (2026-07-27)

The identical single CubeCL kernel is correct under the CUDA and HIP runtimes but wrong
under wgpu, on **real hardware**, not just software adapters:

| Path | Adapter(s) | `mls-bench` max vertex error | Symptom |
| --- | --- | --- | --- |
| wgpu SPIR-V | lavapipe **and real AMD Vulkan** | 1.50 vox | output ≈ input (unprojected) |
| wgpu WGSL | Mesa d3d12-GL **and real AMD DX12** | 62.1 vox | garbage |

- The two wgpu paths fail **differently** (unprojected vs garbage) — so there may be more
  than one miscompiled construct, not a single small workaround.
- Failure is **independent of the environment/adapter**: it reproduces on mature AMD
  Vulkan and DX12 drivers, not only software fallbacks.
- Failure is **independent of the Kahan accumulation**: replacing `compensated_add` with a
  plain `*sum += value` produced byte-identical wrong output, so the summation is not the
  cause — the neighbour-search / output path is not executing correctly under wgpu.

This is strong evidence of a **CubeCL-WGPU code-generation defect or an unsupported-kernel
pattern**. It is not yet proven a compiler bug until a minimal reproducer rules out
undefined/unsupported CubeCL semantics and upstream confirms it — a minimal reproducer is
to be filed with the CubeCL maintainers.

## Do not
- Do not select wgpu via `Backend::Auto`.
- Do not advertise wgpu as a production backend.
- Do not rewrite the validated CUDA/HIP kernel to work around this — the expected value is
  poor and it would require full CUDA/HIP re-validation.

## Revisit
Only once the CubeCL maintainers identify the miscompiled construct or provide a fix/
workaround. The failing case is retained for that purpose.
