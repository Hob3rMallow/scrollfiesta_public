#ifndef SEAM_BAND_CVT_INCLUDED
#define SEAM_BAND_CVT_INCLUDED

#include <stddef.h>
#include <stdint.h>
#include "../common/arena.h"
#include "seam_planes.h"

/*
 * Post-weld seam-band CVT beautification.
 *
 * The recoarsen stage collapses the transient fine weld band back toward the
 * coarse budget, but a guarded shortest-edge collapse is QEM-style decimation
 * minus even QEM's placement solve -- the band ends up collapse-scarred and
 * anisotropic next to the pristine blue-noise CVT interior. This stage
 * re-meshes the band with the SAME CVT/RVD engine the interior was built
 * with, so the welded mesh is CVT-quality everywhere.
 *
 * Mechanism: faces near a detected seam plane are grouped into connected
 * patches; each patch is remeshed by CVT_remesh_pinned -- every patch-boundary
 * vertex (the junction ring against the coarse interior, and any real open
 * boundary: holes, grid edges) is an IMMUTABLE site preserved bit-exact, so
 * the result stitches back index-conformally with no cracks and no change to
 * any open boundary loop (the hierarchical-composability requirement).
 *
 * FAIL-CLOSED per patch: the RVD dual does not structurally guarantee
 * topology, so each patch is accepted only if ALL gates pass --
 *   1. every pinned site is referenced by the dual;
 *   2. the dual's boundary-edge set EQUALS the patch's boundary-edge set
 *      (no dropped junction edge = no crack; no new boundary = no new hole);
 *   3. edge-manifold (every dual edge has <= 2 faces);
 *   3b. vertex-manifold (no bowtie/pinch vertex). A vertex-only pinch passes
 *      gate 3 (all its edges still <= 2 faces) AND gate 4 (still one
 *      vertex-connected component), but the downstream ManifoldGuard splits it
 *      by vertex duplication -- severing a same-turn connection the bridge just
 *      made. This was THE seam-split root cause (wind_audit).
 *   4. one connected component (as the patch was);
 *   5. Euler characteristic unchanged (with 2+3+3b+4: genus unchanged AND the
 *      replacement is homeomorphic to the patch -- it cannot sever a
 *      connection, the guard against RVD cells jumping between close folds).
 * A rejected patch is left exactly as the recoarsen stage produced it.
 *
 * Vertices are never moved or removed (replaced patch interiors are orphaned);
 * new interior CVT verts are appended, out_new_vert_src[i] = a nearby original
 * vert (nearest pinned site's source) for per-vertex attribute propagation.
 * No-op path (np == 0 / nothing accepted) returns the inputs by reference.
 * Returns 0 on success.
 */

typedef struct {
    double band;       /* vox from a detected plane: face eligibility        */
    double target_h;   /* band CVT target edge length (vox)                  */
    int    iters;      /* Lloyd/CWF iterations per patch                     */
    size_t min_faces;  /* patches smaller than this are left alone           */
    double tile;       /* in-plane tile size (vox): patches are split by
                        * tile so one gate failure (a core fold, a dropped
                        * junction edge) rejects ONE tile, not the whole
                        * band lattice                                       */
    int    offset_half;/* pass B: shift the tiling by tile/2 so pass A's
                        * pinned tile borders land tile-INTERIOR and get
                        * re-meshed away (no border scars)                   */
    double fine_frac;  /* process a patch only if >= this fraction of its
                        * edges is < 0.6*target_h -- skips already-CVT
                        * regions on pass B instead of churning them         */
} SeamBandCvtParams;

typedef struct {
    size_t patches;            /* candidate sub-patches                      */
    size_t accepted;
    size_t rejected;           /* total gate failures (kept as-is)           */
    size_t rej_rc;             /*   CVT returned error / empty               */
    size_t rej_pin;            /*   a pinned site unused by the dual         */
    size_t rej_bnd;            /*   boundary-edge set changed                */
    size_t rej_chi;            /*   Euler characteristic changed             */
    size_t rej_manifold;       /*   non-manifold dual edge                   */
    size_t rej_vtx;            /*   bowtie/pinch vertex (would sever @ guard) */
    size_t rej_comp;           /*   component count changed                  */
    size_t skipped_small;
    size_t skipped_clean;      /* already CVT-quality (fine_frac filter)     */
    size_t faces_band_in;      /* band faces before (accepted patches only)  */
    size_t faces_band_out;     /* band faces after  (accepted patches only)  */
    size_t verts_added;
} SeamBandCvtStats;

void SeamBandCvt_default_params(SeamBandCvtParams *p);

int SeamBandCvt_process(Arena_T arena,
                        const float *verts, size_t nv,
                        const int32_t *faces, size_t nf,
                        const SeamPlane *planes, size_t np,
                        const SeamBandCvtParams *params,
                        float **out_verts, size_t *out_nv,
                        int32_t **out_faces, size_t *out_nf,
                        int32_t **out_new_vert_src, size_t *out_n_new,
                        SeamBandCvtStats *st);

#endif
