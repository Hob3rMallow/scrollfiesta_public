#ifndef RAW_SAMPLE_INCLUDED
#define RAW_SAMPLE_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "arena.h"

/* ============================================================================
 * raw_sample.h -- sample raw CT intensity from a grid of per-cube uint8 TIFFs
 * onto arbitrary 3D points. Lazy per-cube loading keyed by the cube's source-
 * space voxel origin (z#####_y#####_x#####.tif). Extracted from obj_bake_raw
 * so both the RAW read-back tool and the overlap-repair module share one
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

/* Area-weighted, normalized per-vertex normals; zero vector where undefined.
 * Returns a malloc'd [nv*3] buffer -- caller frees. */
float *vertex_normals(const float *verts, size_t nv,
                      const int32_t *faces, size_t nf);

#endif
