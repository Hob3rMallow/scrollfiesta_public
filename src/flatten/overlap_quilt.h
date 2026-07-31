#ifndef OVERLAP_QUILT_INCLUDED
#define OVERLAP_QUILT_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"

/* ============================================================================
 * overlap_quilt.h -- resolve UV overlaps in a scroll unroll by (1) splitting
 * co-located layers with a lifted multicut, (2) sampling raw CT per layer and
 * scoring a quilting (minimum-error-boundary) seam energy on fiber structure
 * against the surrounding clean sheet, (3) keeping the layer that best
 * continues that sheet, and (4) relocating the non-chosen SEPARATE-component
 * layers to a winding-consistent empty UV slot (the rest are dropped from the
 * primary texture). Steps mirror the user's proposal; energy is quilting
 * (Efros-Freeman / Kwatra / Zomet), NOT seam-carving saliency.
 *
 * Scope note: overlaps must be bounded by clean single-layer sheet (the anchor)
 * for the energy pick to be well-posed; the fully-fused core has no anchor and
 * is left as-is (flagged).
 * ==========================================================================*/

typedef struct QuiltOpts {
    float  axis_point[3];   /* (z,y,x) scroll axis point (umbilicus) */
    float  axis_dir[3];     /* (z,y,x) axis direction (normalized internally) */
    const char *raw_dir;    /* cubes_RAW dir; NULL disables texture energy
                             * (layers then picked by the radial prior alone) */
    long   chunk;           /* cube edge length, vox (default 128) */
    double radius_gate;     /* min radial gap for two faces to be distinct
                             * layers, vox. <=0 => pitch/2 (default 0) */
    double cell;            /* UV overlap-detect cell size, vox (default 2.0) */
    double tex_du, tex_dv;  /* per-layer texture raster px, vox (default 1.0) */
    double normal_range;    /* RAW normal-max reach, vox (default 2.0) */
    int    normal_samples;  /* RAW normal-max taps (default 5) */
    /* winding prior r ~= spiral_a + spiral_b * phi/2pi; pitch vox/turn.
     * From the ribbon stats JSON; used to gate/tiebreak layer choice and to
     * place relocated leftovers. spiral_b<=0 disables the radial prior. */
    double spiral_a, spiral_b, pitch;
    int    relocate;        /* 1 = Step 4 relocation on, 0 = leave/drop */
    int    cross_component_only; /* only compare faces from disconnected input
                                  * components (local diagnostic mode) */
    int    verbose;
} QuiltOpts;

/* Fill with defaults (axis = Z through umbilicus (0,3405,2878), chunk 128,
 * cell 2, tex 1, range 2 / 5 taps, relocate on). */
void QuiltOpts_default(QuiltOpts *o);

typedef struct QuiltResult {
    /* --- per original vertex --- */
    float   *uv;          /* [nv*2] repaired (relocated leftovers shifted in u) */
    /* --- per original face --- */
    uint8_t *face_keep;   /* [nf]: 1 = emit in the primary texture, else 0 */
    uint8_t *face_dec;    /* [nf]: 0 non-overlap / 1 kept / 2 relocated / 3 dropped */
    int32_t *face_region; /* [nf]: overlap region id, or -1 */
    int32_t *face_layer;  /* [nf]: layer id within its region, or -1 */
    double  *face_energy; /* [nf]: boundary energy for this layer; -1 unscored,
                           * >=1e17 means that no usable anchor seam exists */

    /* --- diagnostics --- */
    size_t n_overlap_faces;
    size_t n_regions;
    size_t n_layers_total;      /* summed over regions */
    size_t max_region_faces;
    size_t n_kept, n_relocated, n_dropped;
    size_t multi_cells_before;  /* UV cells with >1 layer, before repair */
    size_t multi_cells_after;   /* ditto, after (should drop sharply) */
    double energy_mean;         /* mean boundary quilting energy of kept layers */
    int    have_texture;        /* 1 if RAW was sampled */
} QuiltResult;

/* Run the pipeline. verts [nv*3] (z,y,x), faces [nf*3] 0-based, uv [nv*2].
 * All outputs arena-allocated. Returns 0 on success, -1 on degenerate input. */
int Quilt_run(Arena_T arena,
              const float *verts, size_t nv,
              const int32_t *faces, size_t nf,
              const float *uv,
              const QuiltOpts *opts, QuiltResult *out);

/* ---- exported subroutines (shared with the whole-grid SeamOwn pass) ------ */

struct CubeTable;   /* raw_sample.h */

/* face-adjacency pair with the shared-UV-edge length (the multicut weight) */
typedef struct QuiltAdjEdge { int32_t fa, fb; float uvlen; } QuiltAdjEdge;

/* Edge-keyed adjacency build over faces[nf*3]. Returns a malloc'd array
 * (caller frees), count in *out_n. */
QuiltAdjEdge *Quilt_build_adjacency(const int32_t *faces, size_t nf,
                                    const float *uv, size_t *out_n);

/* Axis-perpendicular radius of point p (z,y,x); e1,e2 span the plane. */
double Quilt_point_radius(const float *p, const float ap[3],
                          const double e1[3], const double e2[3]);

/* Rasterize faces flist[nfl] (global ids) into a local (gw x gh) grid at
 * origin (u0,v0), sampling RAW (normal-max) at the barycentric surface point
 * per covered pixel: img = sample, cov = 1 there. */
void Quilt_raster_faces(const int32_t *flist, size_t nfl,
                        const float *verts, const float *uv, const float *nrm,
                        const int32_t *faces,
                        struct CubeTable *ct, double range, int nsteps,
                        double u0, double v0, double du, double dv,
                        int gw, int gh, float *img, uint8_t *cov);

/* Efros min-error-boundary VALUE-domain seam energy of a layer raster against
 * an anchor raster on the same grid (offset-compensated). 1e18 = no seam. */
double Quilt_boundary_energy(const float *vL, const uint8_t *covL,
                             const float *vA, const uint8_t *covA,
                             int gw, int gh);

/* In-process unit tests (radius self-gate, 2-layer split, gradient + boundary
 * energy pick, relocation). Returns 0 if all pass, else #failures. */
int Quilt_selftest(void);

#endif
