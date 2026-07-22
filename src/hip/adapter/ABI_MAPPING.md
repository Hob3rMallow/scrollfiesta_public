# Phase 2 - ScrollFiesta MLS HIP Adapter ABI Mapping

## Status

OBSERVATION: Adapter files were created under `~/mls-hip-cleanroom/cleanroom/clean-side/adapter/` only:

- `mls_scrollfiesta_adapter.h`
- `mls_scrollfiesta_adapter.cpp`
- `adapter_smoke.cpp`
- `Makefile`

OBSERVATION: The adapter exports the upstream-compatible `MLS_project_verts(...)` symbol and uses `objcopy --redefine-sym MLS_project_verts=mls_clean_project_verts` on a copied clean object so the validated clean HIP implementation is not edited.

OBSERVATION: Build and smoke test passed on gfx1201:

```text
rm -rf build
mkdir -p build
g++ -O2 -std=c++17 -Wall -Wextra -Wpedantic -I../impl -c mls_scrollfiesta_adapter.cpp -o build/mls_scrollfiesta_adapter.o
g++ -O2 -std=c++17 -Wall -Wextra -Wpedantic -c adapter_smoke.cpp -o build/adapter_smoke.o
objcopy --redefine-sym MLS_project_verts=mls_clean_project_verts ../impl/build/mls_project_hip.o build/mls_project_hip_renamed.o
ar rcs build/libmls_project_clean_renamed.a build/mls_project_hip_renamed.o
HSA_ENABLE_DXG_DETECTION=1 ROCPROFILER_REGISTER_ENABLED=0 /usr/bin/hipcc --offload-arch=gfx1201 build/mls_scrollfiesta_adapter.o build/adapter_smoke.o build/libmls_project_clean_renamed.a -o build/adapter_smoke
finite=576 total=576 normal_norm_mean=1.000000 first_out_zyx=0.020196,0.000000,-0.000001
```

OBSERVATION: The corrected adapter keeps `MLS_project_verts` as one pass and
also exports `mls_scrollfiesta_lop_project(..., passes, ...)` for standalone
5/20-pass workflows. Its smoke test proves the helper is bit-identical to five
explicit one-pass calls with support rebuilt each pass (`lop5_max_delta=0`).

## Upstream ABI Mapping

| ABI element | OBSERVATION | INFERENCE / adapter action |
| --- | --- | --- |
| License and attribution | Brief identifies upstream `src/`, `include/`, `scripts/`, and `python/` as MIT, Copyright 2026 Nicholas Vining. | Adapter files carry an attribution header. No upstream code was copied into `clean-side/impl/mls_project_hip.hip`. |
| Function signature | `src/common/mls_project.h:45-50` declares `void MLS_project_verts(Arena_T arena, const float *verts, size_t nv, float radius_vox, const float cell_origin[3], float *out_verts, float *out_normals);`. | Adapter exports the same C ABI. |
| Upstream arena role | `src/common/mls_project.c:204` saves the arena, and `:444` restores it. The function allocates entries and cells from that arena between those points. `src/common/arena.h:6` defines `Arena_T` as an opaque pointer; `:17-29` exposes save/restore and allocation lifecycle. | Upstream `Arena_T` is scratch allocation, not persistent MLS grid state. Adapter ignores the incoming scratch arena and constructs a clean HIP arena from the current support cloud for each call. This is representable, but slower than a future cached arena. |
| Verts layout | `src/common/mls_project.c:223-229` reads `verts[i*3+0]` as z, `+1` as y, `+2` as x. `src/common/mls_project.h:31-35` documents output positions as `(z,y,x)`. | Adapter converts upstream z,y,x arrays to clean x,y,z arrays before calling the clean kernel, then converts positions and normals back to z,y,x. |
| Cell origin semantics | `src/common/mls_project.c:210-218` reads `cell_origin[0..2]` as z,y,x offsets. `:223-229` computes cells using `floorf((coord + origin) * inv_R)`. | Clean source uses `floor((coord - index_origin) * inv_cell)` in `~/mls-hip-cleanroom/cleanroom/clean-side/impl/mls_project_hip.hip:56-59` and `:203-208`; adapter converts origin to x,y,z and negates it so the cell index formula is equivalent. |
| Neighbor coverage | Upstream cell size is the radius because `inv_R = 1.0f / radius_vox` at `src/common/mls_project.c:196`, and the neighbor loop visits dz/dy/dx `-1..1` at `:270-302`. | Adapter creates the clean arena with `cell_size_vox = radius_vox`, making clean `span=ceil(radius/cell_size)=1` match the upstream 27-cell search. |
| Coordinate convention at call site | `src/extract/mesh_extract.c:1017-1025` passes `surf_verts`, `surf_nv`, `MLS_PROJECT_RADIUS_VOX`, and `cube_world_origin` into `MLS_project_verts`. Earlier comments at `:980-989` describe voxel-neighbor determinism across cube boundaries. Existing brief identifies the selected Phase 3 cube as a ScrollFiesta prediction cube. | The adapter treats incoming coordinates as voxel coordinates in upstream z,y,x order and preserves the same values aside from order conversion. |
| Output aliasing | Upstream header says inputs and outputs may alias at `src/common/mls_project.h:31`. Clean source rejects overlapping buffers at `~/mls-hip-cleanroom/cleanroom/clean-side/impl/mls_project_hip.hip:694-699`. | Adapter copies queries/support into private vectors and writes back after the clean call, so upstream aliasing is tolerated for position output. |
| `out_normals == NULL` | Upstream CPU has a fallback path when normals are unavailable at `src/common/mls_project.c:434-441`. Clean source requires non-null `out_normals` at `~/mls-hip-cleanroom/cleanroom/clean-side/impl/mls_project_hip.hip:679-683`. | Adapter allocates a discard normals buffer when upstream passes null. The relevant ScrollFiesta call sites pass normals for MLS projection. |

## Blockers / Unresolved Items

1. RESOLVED OBSERVATION: Upstream applies one projection per call and loops at
   the caller. `pipeline_constants.h:73-77` sets extract to 5 and resplit to 20;
   `mesh_extract.c:1019-1028` and `mesh_resplit.c:147-153` refresh support from
   the previous pass. The clean public ABI now requests one internal update.
   The adapter-level LOP helper rebuilds the clean arena from each previous
   output, so it performs 5 or 20 moving-support passes rather than 25/100
   fixed-support updates.

2. OBSERVATION: The clean kernel status is only available via `mls_clean_last_status()` and is not represented in upstream's void ABI. Upstream CPU allocation failures abort through Arena allocation rules, while clean HIP failures are status-coded.

   INFERENCE: For full pipeline integration, the adapter should expose/log clean status at the call site or fail visibly, not silently leave outputs untouched.

3. OBSERVATION: The adapter constructs the clean arena per call from the current source vertices. This mirrors the upstream CPU call contract but includes host reordering plus grid construction each time.

   INFERENCE: This is acceptable for correctness validation, but may leave performance below the Phase 4 target. A cached arena strategy would need caller cooperation because the upstream arena is scratch-only and not typed as persistent MLS state.

## Commands

```bash
cd ~/mls-hip-cleanroom/cleanroom/clean-side/adapter
make clean all
HSA_ENABLE_DXG_DETECTION=1 ROCPROFILER_REGISTER_ENABLED=0 ./build/adapter_smoke
```
