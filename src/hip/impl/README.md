# Clean-side HIP MLS projection

This directory contains an independent HIP implementation of one MLS midpoint
projection per public ABI call. It targets `gfx1201` and exports the exact
required C symbol:

```c
void MLS_project_verts(void* arena,
                       const float* verts, size_t nv,
                       float radius_vox,
                       const float* cell_origin,
                       float* out_verts,
                       float* out_normals);
```

## Compatibility boundary

The upstream `arena` type and layout are unresolved. This implementation does
not consume, recognize, or claim compatibility with an upstream arena. Its
`void* arena` argument must point to an arena returned by
`mls_clean_arena_create`. Existing callers therefore need an independently
authorized adapter or a resolved public arena contract before this object can
be called a drop-in replacement. No choice below is claimed to reproduce an
uncharacterized upstream edge case or coordinate convention.

The clean construction API takes packed FP32 `(x,y,z)` support points. Queries,
support points, and outputs use that same coordinate space and voxel units.
`index_origin` is only the origin of the uniform-grid indexing frame; it is not
added to or subtracted from geometric positions. The `cell_origin` triple on a
projection call must numerically equal the arena's construction origin. This
check makes the clean-side convention explicit and prevents a silently
misindexed call.

## Defined numerical behavior

For each valid query, one GPU thread performs the explicitly requested
fixed-support update count. The public ABI requests exactly one update. Each update:

1. re-finds samples at strict Euclidean distance `d < radius_vox`;
2. accumulates FP32 Wendland-C2 weights
   `(1-d/R)^4 * (4*d/R+1)` and the weighted centroid with compensated sums;
3. re-enumerates the same spatial cells for an FP32 covariance centered at that
   centroid;
4. solves the symmetric 3-by-3 covariance in upstream axis order using FP64 for
   the small eigensolve only, then emits an FP32 unit eigenvector; and
5. projects the current point onto the centroid plane.

The returned normal is the sign-oriented normal used by the final requested
projection. A negative dot product with `(1,1,1)`
causes a flip. When that dot is exactly zero, the first nonzero component is
made positive as a clean-side deterministic tie rule.

The spatial index is a deterministic open-addressed uniform-cell hash table.
Cells are visited in z/y/x loop order; points within a cell retain input order.
The cell width is chosen by the caller at arena construction. A width near
`radius_vox/3` through `radius_vox/2` is a useful starting point; for production
`R=12`, `4.0f` is the suggested initial choice. Calls for which
`ceil(radius/cell_width) > 64` are rejected defensively.

The wrapper copies host queries to reusable arena-owned device buffers, launches
on an arena-owned nonblocking stream, copies both output arrays back, and
synchronizes that stream before returning. Buffers grow geometrically and are
reused. An arena is bound to the active HIP device on construction. Calls on a
single arena are serialized with a host mutex; distinct arenas can run
independently. Destroying an arena concurrently with another operation is not
supported.

## Valid inputs and fail-safe cases

For `nv > 0`, all three vertex/output arrays must be aligned, valid host arrays
of `3*nv` floats and must not overlap. Radius and all coordinates must be finite,
and radius must be positive. Support points are validated during construction.
After a valid arena, radius, and matching origin are supplied, `nv == 0` is a
successful no-op and the vertex/output pointers may be null.

If one query has a non-finite coordinate, no positive-weight neighborhood, a
repeated/insufficiently separated smallest eigendirection, or another numerical
failure, that query receives a defined fail-safe result: its original finite
position (or zero for a non-finite input) and normal `(0,0,0)`. Other queries
still complete. `mls_clean_last_status` and `mls_clean_arena_get_info` expose the
aggregate status and number of affected queries without changing the required
void ABI. Argument errors are detected before launch and leave outputs
untouched. Output contents after an asynchronous HIP/runtime failure are not
specified.

These defensive policies cover behavior that remains uncharacterized upstream;
they are not upstream-parity claims.

## Build and checks

ROCm with a `gfx1201`-capable `hipcc` is required:

```sh
make
```

The default build produces `build/mls_project_hip.o`,
`build/libmls_project_clean.a`, and `build/mls_clean_runner` using
`hipcc --offload-arch=gfx1201`. The `check` target validates C-header syntax,
the ELF object, archive membership, runner presence, and the required unmangled
exported symbols. Override `HIPCC`, `ARCH`, or flags in the usual Make style if
needed.

The runner is the executable bridge for `run_hip_vs_oracle.sh`. It consumes
`MLS_FIXTURE` (informational), `MLS_SAMPLES_RAW`, `MLS_QUERIES_RAW`,
`MLS_PARAMS`, `MLS_POSITIONS_RAW`, and `MLS_NORMALS_RAW`. It validates the
declared raw sizes, little-endian FP32 xyz format, one-iteration count, radius
bit pattern, and origin; builds cells at `radius/2`; requires clean status `OK`;
and writes raw little-endian projected positions and normals. Pass
`build/mls_clean_runner` as the verifier's runner argument.

Minimal construction and use:

```c
#include "mls_project_clean.h"

float origin[3] = {0.0f, 0.0f, 0.0f};
mls_clean_arena *arena = NULL;
if (mls_clean_arena_create(points_xyz, point_count, 4.0f,
                           origin, &arena) != MLS_CLEAN_STATUS_OK) {
    /* inspect mls_clean_last_status() */
}

MLS_project_verts(arena, queries_xyz, query_count, 12.0f, origin,
                  projected_xyz, normals_xyz);
mls_clean_status status = mls_clean_last_status();
mls_clean_arena_destroy(arena);
```

Numerical and topology parity still require the authorized golden captures and
black-box characterization described in the functional specification. This
source alone cannot close those integration unknowns.
