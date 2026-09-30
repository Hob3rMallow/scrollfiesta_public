#ifndef ASM_FLATTEN_INCLUDED
#define ASM_FLATTEN_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"

/* ============================================================================
 * asm_flatten.h -- intrinsic flattening of ONE mesh chart.
 *
 * Tutte disk embedding (outer boundary loop on a circle, symmetric positive
 * mean-value weights, holes capped by a virtual vertex for the embedding
 * only) followed by orientation-guarded ARAP local/global iterations (Liu, Zhang, Gotsman &
 * Gortler 2008) on the real triangles.  All linear algebra is a Jacobi
 * preconditioned conjugate gradient in double on a CSR matrix built from the
 * chart's own adjacency: charts are 10^2..10^4 vertices and the flattener
 * runs inside an OpenMP loop over cubes, where TAUCS (process-global state)
 * must not be called.
 *
 * The map is orientation-preserving with respect to the face winding: a face
 * whose 2-D orientation is negative counts as flipped.  The chart frame is
 * fixed from the normalized Tutte seed (mean at the origin, principal axis
 * along +u), so two flattenings of the same chart are bitwise comparable.
 * The gauge stays fixed during optimization and encoding. Every accepted step
 * preserves positive orientation throughout its interval and decreases the
 * ARAP objective. Every candidate also guards boundary crossings across all
 * loops at the strict intrinsic-chart tolerance (1e-7 voxels), independently
 * of the assembled material-contact allowance.
 * The final float32 cache coordinates are checked separately; a folded,
 * collapsed or boundary-crossing stored map is never reported as successful.
 * ==========================================================================*/

enum {
    ASM_FLAT_OK = 0,
    ASM_FLAT_EMPTY,          /* fewer than 3 vertices or no faces */
    ASM_FLAT_BAD_INDEX,      /* face index out of range */
    ASM_FLAT_NONMANIFOLD,    /* edge with >= 3 faces or a pinch vertex */
    ASM_FLAT_NO_BOUNDARY,    /* closed surface: nothing to pin the disk to */
    ASM_FLAT_DEGENERATE,     /* zero-area geometry */
    ASM_FLAT_SOLVE,          /* the linear solve did not converge */
    ASM_FLAT_FOLDED          /* no strictly oriented, representable map */
};

typedef struct AsmFlattenStats {
    size_t nv, nf;
    size_t n_boundary_loops;
    size_t n_holes_capped;
    size_t n_flipped;        /* faces with negative 2-D orientation */
    double area3d, area_uv;  /* vox^2 (uv area is unsigned) */
    double sigma_lo, sigma_hi; /* extreme singular values over faces */
    double stress_frac;      /* faces with a singular value outside [1-band, 1+band] */
    int    arap_iters;
    int    metric_refine_iters; /* native original-metric continuation, if needed */
    int    ok;               /* 1 = usable map in out_uv */
    int    fail_reason;      /* ASM_FLAT_* */
    int    seed_mode;        /* seed of the returned map: 0 Tutte with uniform hole caps (historical), 1 capped
                              * conformal, 2 Tutte with mean-value hole caps, 3 plain conformal */
    int    fail_stage;       /* diagnostics: 1 seed, 2 ARAP, 3 orientation, 4 float32 guard, 5 self-contact,
                              * 6 valid map without the chart-metric certificate (rescue only) */
    size_t fail_contacts;    /* diagnostics: self-contact pairs at stage 5 */
} AsmFlattenStats;

/* Flatten the chart (xyz [nv*3], faces [nf*3] chart-local).  Writes uv
 * [nv*2] and the stats.  Scratch is taken from `arena` and restored before
 * returning.  Returns 0 when a usable map was produced, -1 otherwise
 * (stats->fail_reason says why; out_uv is then unspecified). */
int AsmFlatten_chart(Arena_T arena,
                     const float *xyz, size_t nv,
                     const int32_t *faces, size_t nf,
                     double band, int max_iters,
                     float *out_uv, AsmFlattenStats *stats);

/* Copy the chart's faces without its zero-area ones (the test that fails a chart as
 * ASM_FLAT_DEGENERATE) and renumber the vertices they still reference, keeping their relative
 * order.  The mesher closes zero-width slits with pairs of collinear faces; dropping them loses no
 * area and leaves a hole whose coincident copies the recovery seeds open like any pinch point.
 * out_faces [nf*3], out_vmap [nv] (new index of each vertex, -1 = unreferenced).  Returns the kept
 * face count; *out_nv is the kept vertex count. */
size_t AsmFlatten_drop_degenerate(const float *xyz, size_t nv, const int32_t *faces, size_t nf,
                                  int32_t *out_faces, int32_t *out_vmap, size_t *out_nv);

/* 1 when the face fails the flattener's zero-area test (face = 3 chart-local indices). */
int AsmFlatten_face_degenerate(const float *xyz, const int32_t *face);

/* AsmFlatten_chart for a chart that is being rescued rather than kept: the historical map is taken
 * only when it carries the audit's chart-metric certificate (every face with singular values in
 * [0.75, 1.25], 95% of the area within [0.9, 1.1]); otherwise the recovery seeds run, and a map
 * without the certificate is a failure (fail_stage 6). */
int AsmFlatten_chart_certified(Arena_T arena, const float *xyz, size_t nv, const int32_t *faces, size_t nf,
                               double band, int max_iters, float *out_uv, AsmFlattenStats *st);

/* Faces to excise from a crumpled chart so that the rest can flatten.  A cone vertex is an interior
 * vertex whose angle sum differs from 2 pi by more than cone_rad, a boundary vertex whose fan
 * exceeds 2 pi, or a pinch copy; the mesher leaves them where a fold or a fused layer ends inside a
 * sheet (PHerc0139 21x21x21: every unflattenable chart, sum |K| 20-600 rad in a few clusters).
 * Every face on a cone vertex is marked, grown by `ring` vertex rings; each marked cluster that does
 * not reach the boundary is joined to it by the faces along its shortest 3-D path (a hole around net
 * curvature cannot flatten), and a vertex left with a split fan has its fan cleared.  excise[nf]
 * (1 = excise).  seed[nf] (or NULL) marks further faces to excise like cone fans.  With cut_holes
 * every hole loop is first joined to the outer loop by a spanning tree
 * of such paths (a band cannot flatten), and clusters join that tree.  Returns the number of marked
 * faces (0: nothing to excise). */
size_t AsmFlatten_crumple_faces(Arena_T arena, const float *xyz, size_t nv, const int32_t *faces, size_t nf,
                                double cone_rad, int ring, int cut_holes, const uint8_t *seed, uint8_t *excise);

/* Last resort for a chart no seed flattens: excise its crumpled regions (ring 1; ring 1 with every
 * hole cut to the outer loop; ring 3 with the cuts; until 90% of the area is certified) and flatten
 * each remaining connected piece of at least min_piece_area with AsmFlatten_chart_certified.  The
 * excision certifying the most area wins.  piece[nf] = the face's piece (0 = the largest certified
 * piece) or -1 (excised, uncertified or small); uv[nv*2] holds each certified piece's map at its
 * vertices (pieces share no vertex); stats are the largest piece's.  Returns 0 when at least one
 * piece is certified, -1 otherwise.  A piece that fails is rescued once more from the faces its
 * conformal seed folds or puts in contact (excised like cone fans, within that piece only). */
int AsmFlatten_rescue_crumpled(Arena_T arena, const float *xyz, size_t nv, const int32_t *faces, size_t nf,
                               double band, int max_iters, double min_piece_area,
                               int32_t *piece, float *uv, size_t *n_pieces, AsmFlattenStats *st);

/* Per-face singular values of an existing (uv) map against the 3-D metric.
 * sig_lo/sig_hi are [nf]; det_sign is [nf] (+1 / -1 / 0).  Any output may be
 * NULL. */
void AsmFlatten_face_singular(const float *xyz, const int32_t *faces, size_t nf,
                              const float *uv, double *sig_lo, double *sig_hi,
                              int8_t *det_sign);

/* In-process unit tests: planar grid, half cylinder, cone, holed rectangle,
 * spherical cap (must be flagged stressed), closed tetrahedron (must fail).
 * Returns 0 on success, else the failure count. */
int AsmFlatten_selftest(void);

#endif
