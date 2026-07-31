#ifndef CVT_REMESH_INCLUDED
#define CVT_REMESH_INCLUDED

/*
 * cvt_remesh — uniform CVT remeshing of a triangle mesh via restricted Voronoi.
 *
 * Seed ~target_nsites points on the surface (interior area-weighted + boundary sites
 * on the loop), then iterate the CWF variational update — for each interior cell a
 * fixed-cell 3x3 solve of (lambda_NA*N_i + lambda_CVT*m_i*I) x = lambda_NA*b_i +
 * lambda_CVT*m_i*c_i (normal-anisotropy N_i,b_i + CVT centroid c_i) with a decaying
 * lambda_CVT; boundary sites relax along their loop. The output mesh is the
 * Restricted-Delaunay dual of the converged sites.
 *
 * Vertices are float[nv*3] in (z,y,x); faces int32[nf*3], 0-based. Outputs are
 * fresh arena allocations.
 */

#include "../common/arena.h"
#include "sizing_field.h"
#include <stddef.h>
#include <stdint.h>

typedef struct {
    int      n_iters;        /* iterations                       (default 50)  */
    float    lambda_na;      /* E_NA normal-anisotropy weight     (default 1.0) */
    float    lambda_cvt0;    /* initial E_CVT weight              (default 1.0) */
    float    cvt_decay;      /* lambda_CVT *= this each iter      (default 0.95)*/
    float    lambda_cvt_min; /* floor on lambda_CVT               (default 1e-3)*/
    uint32_t seed;           /* seeding RNG (0 => 12345u)                       */
    int      verbose;
} CvtOpts;

void CVT_default_opts(CvtOpts *o);

/* Remesh to a CVT of the input surface.
 *
 * If `field == NULL`: uniform density at `target_nsites` generators — the
 * original, unchanged behaviour (standalone tool / tests use this).
 *
 * If `field != NULL`: GRADED density from the sizing field (dense near the cube
 * owned-box planes for weldability, coarse interior). `target_nsites` is ignored
 * — the site count is DERIVED from the field (N = 2*integral(1/h^2 dA)). */
int CVT_remesh(Arena_T arena,
               const float *in_verts, size_t in_nv,
               const int32_t *in_faces, size_t in_nf,
               size_t target_nsites,
               const CvtOpts *opts,
               const CvtField *field,
               float **out_verts, size_t *out_nv,
               int32_t **out_faces, size_t *out_nf);

/* Self-test for the KD-tree surface projection: project_to_mesh_kd must match the
 * O(nf) brute force on every query, including meshes with wildly-varying triangle
 * sizes (validates the ball-radius superset). Returns 0 on pass, <0 on failure. */
int CVT_projection_selftest(void);

/* PINNED-BOUNDARY CVT remesh -- for remeshing a PATCH of a larger mesh so the
 * result stitches back conformally (the post-weld seam-band beautification).
 *
 * Every boundary vertex of the input patch becomes an IMMUTABLE site: not
 * resampled by arc length, never relaxed, position preserved BIT-EXACT in the
 * output. Interior sites are seeded at ~target_h spacing (count derived as
 * 2*area/h^2), kept clear of the boundary polyline by 0.85x the LOCAL boundary
 * segment length (scale-aware: a coarse junction ring gets a wide berth, a
 * fine hole rim a narrow one, so neither wedging nor depletion rings appear),
 * and CWF-relaxed + projected onto the patch surface as usual. The dual is the
 * Restricted-Delaunay of pinned + interior sites.
 *
 * Output sites [0, *out_n_pin) are the pinned boundary verts;
 * out_site_src[i] = the input vert index each pinned site reproduces (in
 * boundary-loop order). Sites >= *out_n_pin are new interior verts
 * (out_site_src[i] = -1 for those). The CALLER must verify conformity before
 * stitching (every pinned site referenced by the dual; every consecutive
 * boundary pair still an edge) and reject the patch otherwise -- the RVD dual
 * does not structurally guarantee either (fail-closed is the contract).
 *
 * Returns 0 on success; nonzero = no usable output (caller keeps the patch). */
int CVT_remesh_pinned(Arena_T arena,
                      const float *in_verts, size_t in_nv,
                      const int32_t *in_faces, size_t in_nf,
                      double target_h,
                      const CvtOpts *opts,
                      float **out_verts, size_t *out_nv,
                      int32_t **out_faces, size_t *out_nf,
                      int32_t **out_site_src, size_t *out_n_pin);

#endif /* CVT_REMESH_INCLUDED */
