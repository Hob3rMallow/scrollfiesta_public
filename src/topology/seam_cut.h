#ifndef SEAM_CUT_INCLUDED
#define SEAM_CUT_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"

/* ============================================================================
 * seam_cut.h -- Seamster (Sheffer & Hart, IEEE Vis 2002): cut a connected
 * triangle mesh into an open, genus-0, single-boundary-loop disk so it can be
 * planar-parameterized (LSCM/ABF).
 *
 * Faithful to the paper:
 *   Section 5   -- region distortion D_r(n) (angle deficit over a growing
 *                  region, Eqs 1-2; multi-loop avoidance; overlap filtering).
 *   Section 6.1 -- terminal selection (Algorithm 1) + all boundary vertices.
 *   Section 6.2 -- approximate minimal Steiner tree via front propagation.
 *   Section 6.3 -- genus reduction (shortest non-separating loops).
 *   cut surgery -- open the mesh along the seam tree (vertex duplication).
 *
 * Visibility V(e) is meaningless for a buried scroll, so V(e)=1 and the edge
 * weight is just length W(e)=||e|| (the paper permits restricting the view
 * domain; uniform = no occlusion preference).
 * ==========================================================================*/

typedef struct {
    double accept_frac;  /* D = accept_frac * D(M): fraction of total distortion
                          * left uncut (smaller -> more cut points). Default 0.0
                          * (cut all positive-distortion terminals). */
    int    r_max;        /* max region-radius level for D_r (paper uses 5). */
} SeamCutOpts;

typedef struct {
    float   *verts;      /* [nv*3] cut mesh (seam vertices duplicated) */
    int32_t *faces;      /* [nf*3] */
    size_t   nv, nf;
    int32_t *vmap;       /* [nv] new-vertex -> original-vertex index */
    /* diagnostics */
    long     genus_before, genus_after;
    long     loops_before, loops_after;   /* boundary-loop count */
    long     n_terminals;
    long     seam_edges;
    int      is_disk_after;
} SeamCutResult;

/* Cut a single connected component into a topological disk. All outputs are
 * arena-allocated. Returns 0 on success, -1 on failure (non-manifold beyond
 * what cutting handles, degenerate, etc.). */
int SeamCut_run(Arena_T arena, const float *verts, size_t nv,
                const int32_t *faces, size_t nf,
                const SeamCutOpts *opts, SeamCutResult *out);

/* Manifold repair: split non-manifold VERTICES (bowties: two fans meeting at a
 * point) and isolate non-manifold EDGES (>2 incident faces) by duplicating the
 * shared vertices, so each output piece is a clean 2-manifold. Idempotent on
 * already-manifold meshes (output == input). This is Seamster's 2-manifold
 * precondition step -- run it, then re-split into connected components, then cut
 * each. Outputs arena-allocated; out_vmap maps new->original vertex. Returns 0.*/
int SeamCut_repair_manifold(Arena_T arena, const float *verts, size_t nv,
                            const int32_t *faces, size_t nf,
                            float **out_verts, size_t *out_nv,
                            int32_t **out_faces, size_t *out_nf,
                            int32_t **out_vmap);

/* ---------------------------------------------------------------------------
 * Handle-loop enumeration (for bridge localization). The genus-reduction core
 * repeatedly cuts the shortest NON-SEPARATING loop until genus 0; each such loop
 * encircles one handle -- i.e. one inter-wrap bridge that fused two scroll wraps.
 * This exposes those loops WITHOUT modifying the input mesh.
 * ------------------------------------------------------------------------- */
typedef struct {
    int32_t *verts;   /* loop vertices, as indices into the INPUT mesh (edge
                       * endpoints; a shared vertex may repeat). Arena-allocated. */
    size_t   n;       /* entries in verts[] */
    double   length;  /* summed edge length of the loop (shorter = tighter neck) */
} HandleLoop;

/* Enumerate the genus-many shortest non-separating loops (one per handle) of a
 * single connected component, in INPUT-vertex index space. Does NOT alter the
 * input. Returns 0 (with *out_n == genus, 0 if genus 0); -1 on bad input. All
 * outputs arena-allocated. */
int SeamCut_handle_loops(Arena_T arena, const float *verts, size_t nv,
                         const int32_t *faces, size_t nf,
                         HandleLoop **out_loops, size_t *out_n);

/* Sever only the SHORT handles: open (cut) every non-separating loop shorter
 * than max_loop_len voxels -- the thin BPA bridges -- and leave longer scroll-
 * scale handles intact (pass max_loop_len <= 0 to cut ALL, like SeamCut_run's
 * genus reduction). Manifold-preserving (reuses the tested cut surgery). The
 * opened mesh is returned in *out_*; *out_loops_cut counts handles severed.
 * Arena-allocated. Returns 0; -1 on bad input. */
int SeamCut_sever_short_handles(Arena_T arena, const float *verts, size_t nv,
                                const int32_t *faces, size_t nf, double max_loop_len,
                                float **out_verts, size_t *out_nv,
                                int32_t **out_faces, size_t *out_nf,
                                long *out_loops_cut);

/* In-process unit tests. Returns 0 if all pass, else number of failures. */
int SeamCut_selftest(void);

#endif
