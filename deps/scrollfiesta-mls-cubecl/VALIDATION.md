# Validation status

## Production-equivalence result

Validated on 2026-07-23 in ROCm 7.2.4 WSL on an RX 9070 (`gfx1201`). The input
is ScrollFiesta's captured 331,013-vertex Component 1 cloud from
`z16128_y02560_x07680.tif`, radius 12, in cube-local z-y-x coordinates.

All three implementations consumed the identical starting array:

- original ScrollFiesta CPU `MLS_project_verts`;
- corrected clean-room HIP C++;
- CubeCL HIP through the one-pass current-support API.

### Five moving-support passes

Three measured repetitions after one warm-up:

| Implementation | Median | Speedup vs CPU |
| --- | ---: | ---: |
| ScrollFiesta CPU | 11,014.5 ms | 1.00x |
| clean HIP C++ | 869.455 ms | 12.67x |
| CubeCL HIP | 927.511 ms | 11.88x |

CubeCL final geometry versus CPU:

```text
max position = 0.002052250 vox
RMS position = 0.000010215 vox
max normal    = 0.024687336 degrees
zero-normal mismatches = 0
weld safety (< 0.25 vox) = PASS
strict per-vertex parity = FAIL
```

CubeCL versus clean HIP passes the strict comparison: max position
`0.000145759` vox, RMS `0.000002010` vox, max normal `0.003872668` degrees.

Raw report: `logs/real-equivalence-20260723T024844Z.log`.

### Twenty moving-support passes

One measured repetition:

| Implementation | Time | Speedup vs CPU |
| --- | ---: | ---: |
| ScrollFiesta CPU | 51,757.0 ms | 1.00x |
| clean HIP C++ | 4,956.23 ms | 10.44x |
| CubeCL HIP | 5,177.33 ms | 10.00x |

CubeCL final geometry versus CPU:

```text
max position = 0.039885289 vox
RMS position = 0.000104474 vox
max normal    = 1.690948594 degrees
zero-normal mismatches = 0
weld safety (< 0.25 vox) = PASS
strict per-vertex parity = FAIL
```

CubeCL versus clean HIP: max position `0.002622004` vox, RMS `0.000009656`
vox, max normal `0.143124914` degrees.

Raw report: `logs/real-equivalence-20260723T024942Z.log`.

## Interpretation

The production operator and caller-loop semantics are now correct. CubeCL HIP
is within 3-6% of clean HIP C++ while retaining one shared CPU/HIP kernel. It
passes the downstream weld-safety boundary after both real caller loops.

Strict Component 1 parity remains failed, consistently with the clean HIP
baseline: a small set of near-tied covariance neighborhoods is sensitive to
FP32 support-state divergence. The RMS errors are within threshold and the
maximum displacement remains below the weld limit. These results establish a
credible replacement candidate, not bitwise equivalence or a speed win over
handwritten HIP.

## Regression gates

- Safe Rust unit tests: PASS.
- C ABI explicit-arena smoke test: PASS.
- Legacy ABI accepts foreign scratch and rebuilds current support: PASS.
- CubeCL current-support CPU/Rust comparison: PASS.
- CubeCL CPU compile/execution: PASS.
- CubeCL HIP compile/execution on `gfx1201`: PASS.
- `MLS_project_verts` shared-library export: PASS.
- Five-pass real three-way comparison: weld PASS, strict parity FAIL.
- Twenty-pass real three-way comparison: weld PASS, strict parity FAIL.

The earlier fixed-support synthetic hardware report remains at
`logs/gfx1201-test-20260723T005731Z.log`; it is retained as historical evidence,
not the production-equivalence result.

The final post-redesign synthetic CPU/HIP and exported-symbol gate is recorded
at `logs/final-gfx1201-gate.log`.
