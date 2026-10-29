#ifndef RAW_SAMPLE_INCLUDED
#define RAW_SAMPLE_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "arena.h"

/* ============================================================================
 * raw_sample.h -- sample raw CT intensity from either a grid of per-cube uint8
 * TIFFs or an uncompressed uint8 OME-Zarr v2 level 0.  Both are addressed in
 * source-space voxel coordinates and loaded lazily one chunk at a time.
 * Extracted from obj_bake_raw so RAW read-back and geometry repair share one
 * sampler. Deps: arena, tiff_io, libm.
 * ==========================================================================*/

/* Lazy cube table. slot[] states: NULL = not tried, an internal sentinel =
 * tried and missing (grid edge), else a loaded chunk^3 uint8 volume. Public
 * (not opaque) so callers can read counts and, in tests, inject slots. */
typedef struct CubeTable {
    Arena_T   arena;
    char      dir[2048];
    long      chunk;            /* cube edge length, vox */
    long      cz0, cy0, cx0;    /* min cube index covered by the table */
    long      nz, ny, nx;       /* table dims, cubes */
    uint8_t **slot;             /* [nz*ny*nx] */
    int       n_loaded, n_missing;
    size_t    n_outside;        /* requested bbox chunks outside Zarr shape */
    int       is_zarr;          /* dir/0/.zarray, uint8 + compressor:null */
    long      shape[3];         /* zarr level-0 shape; unused for TIFF cubes */
} CubeTable;

/* Size the table from the vertex bbox padded by `pad` vox (the normal-sampling
 * reach). Cubes load lazily on first fetch. Returns 0 on success, -1 on
 * degenerate/unreasonable input. All allocations from `arena`. */
int cubetable_init(CubeTable *ct, Arena_T arena, const char *dir,
                   long chunk, const float *verts, size_t nv, double pad);

/* Voxel fetch at integer (iz,iy,ix); returns 0..255, or -1 if the covering
 * cube is missing / out of the table. Loads the cube on first touch. */
int cube_fetch(CubeTable *ct, long iz, long iy, long ix);

/* Load EVERY slot in the table bbox now (present -> volume, absent -> the
 * missing sentinel). After this, cube_fetch/sample_* never mutate the table,
 * so concurrent reads from OpenMP workers are safe. Serial; returns the
 * number of cubes loaded. */
int cubetable_prewarm_all(CubeTable *ct);

/* Coverage proof after prewarm.  Expected chunks are only chunks intersecting
 * the physical Zarr shape; bbox halo outside that shape is counted separately
 * in n_outside and samples there return no datum.  A production sampler must
 * not treat an absent in-volume chunk as legitimate zero-valued CT. */
size_t cubetable_expected_chunks(const CubeTable *ct);
int cubetable_is_complete(const CubeTable *ct);

/* After prewarm: the number of distinct MISSING in-volume chunks that some
 * face can sample -- a chunk intersecting the face's (z,y,x) bbox padded by
 * `pad` vox plus the trilinear +-1, the same padding cubetable_init gives the
 * whole mesh.  A masked Zarr stores no chunk outside the scroll, so a mesh
 * bbox can be incomplete while every chunk its faces reach is present; zero
 * here is that proof.  Faces index verts[3*v..3*v+2].  reached (nullable,
 * cubetable_expected_chunks entries in slot order) is set to 1 for each
 * reached missing chunk, 0 elsewhere. */
size_t cubetable_missing_reached(const CubeTable *ct, const float *verts,
                                 const int32_t *faces, size_t nf, double pad,
                                 uint8_t *reached);

/* The faces whose padded bbox (same padding) intersects a chunk flagged in
 * `slots` (slot order, e.g. the `reached` output above): face_mask[f] = 1
 * for those, 0 otherwise.  Returns their count. */
size_t cubetable_faces_reaching(const CubeTable *ct, const float *verts,
                                const int32_t *faces, size_t nf, double pad,
                                const uint8_t *slots, uint8_t *face_mask);

/* Trilinear sample at (z,y,x); weights renormalize over available corners.
 * Returns -1.0 when no corner has data. */
double sample_trilinear(CubeTable *ct, double z, double y, double x);

/* Local 2-D structure tensor of the RAW volume on the tangent plane through p.
 * tu/tv are the intended unrolled u/v directions in source (z,y,x) space; they
 * are normalized and orthogonalized internally. A 5x5 patch spanning
 * +/-radius is sampled, with central gradients accumulated over its inner 3x3.
 * `quality` combines gradient strength, tensor coherence, 4-RoSy alignment to
 * the tu/tv axes, and a normal-profile ridge-center gate. The gate prevents a
 * tilted patch at a bright sheet shoulder from aliasing the across-sheet edge
 * into a fake tangent tensor. Quality is 0 for a flat/unsupported/off-ridge
 * patch and approaches 1 for a centered, coherent axis-aligned fiber pattern.
 * Returns 0 when enough RAW samples exist, -1 otherwise. */
typedef struct RawTangentTensor {
    double energy;          /* mean |grad I|^2 in raw^2 / vox^2 */
    double rms_gradient;    /* sqrt(energy) */
    double coherence;       /* (lambda_max-lambda_min)/(lambda_max+lambda_min) */
    double axis_alignment;  /* 1 axis-aligned, 0 diagonal (4-RoSy) */
    double normal_asymmetry;/* |I(p+0.75n)-I(p-0.75n)|, raw units */
    double ridge_center;    /* [0,1], suppresses sheet-shoulder leakage */
    double quality;         /* bounded combined score [0,1] */
} RawTangentTensor;

int sample_tangent_tensor(CubeTable *ct, const double p[3],
                          const double tu[3], const double tv[3],
                          double radius, RawTangentTensor *out);

/* Max of trilinear samples at nsteps points along +/- range*normal (the
 * normal-max multi-tap: picks up the bright papyrus core even when a boundary
 * vertex sits half in air). Single point-sample when the normal is degenerate
 * or range<=0. Returns -1.0 when nothing sampled. p and n are (z,y,x). */
double sample_vertex(CubeTable *ct, const float *p, const float *n,
                     double range, int nsteps);

#endif
