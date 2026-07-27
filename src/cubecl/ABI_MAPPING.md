# ABI mapping — CubeCL MLS backend ↔ ScrollFiesta `MLS_project_verts`

The CubeCL crate builds `libherculaneum_mls_cubecl.so`, which **exports this repo's
`MLS_project_verts` C ABI directly** — no symbol-rename adapter is needed. `cube_mesh_cubecl`
links the `.so` in place of `common/mls_project.o`.

Signature (matches `common/mls_project.h`):

```c
void MLS_project_verts(Arena_T arena, const float *verts, size_t nv,
                       float radius_vox, const float cell_origin[3],
                       float *out_verts, float *out_normals);
```

Contract, matched point-for-point (verified by direct element-wise parity against both the
CPU oracle and pscamillo's production CUDA kernel — see `CUDA_VALIDATION.md`):

| Aspect | ScrollFiesta contract | CubeCL backend |
| --- | --- | --- |
| Coordinate order | z, y, x | z, y, x |
| Scratch `arena` | scratch allocator only | ignored; the kernel builds its own spatial index per call |
| Cell hashing | `floor((coord + cell_origin) / R)`; within-cell sort by (z,y,x) | same |
| Weight kernel | Wendland C² `(1 − r/R)^4 (4r/R + 1)` | same |
| Projection | LOP (μ=0): weighted centroid + smallest-eigenvector tangent projection | same |
| Normal sign | fixed `(1,1,1)`-dot convention | same |
| `< MLS_MIN_NEIGHBOURS` | write original position + zero normal (welded downstream) | same |
| In-place aliasing | supported | supported (ping-pong at the caller, as in `extract/mesh_extract.c`) |

The AMD (HIP) and NVIDIA (CUDA) backends are hardware-validated. The wgpu backend compiles from
the same source but currently produces incorrect output — see `WGPU_STATUS.md`. Full record +
reproducible harnesses: <https://github.com/altommo/scrollfiesta-mls-cubecl>.
