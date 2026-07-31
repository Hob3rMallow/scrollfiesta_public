#ifndef ISOTROPIC_REMESH_INCLUDED
#define ISOTROPIC_REMESH_INCLUDED

#include <stddef.h>
#include <stdint.h>
#include "../common/arena.h"

/*
 * Remesh_isotropic — incremental isotropic remeshing (Botsch-Kobbelt /
 * Surazhsky-Gotsman) to IMPROVE triangle quality at a FIXED target face count.
 * Meant to run AFTER QEM decimation: QEM is collapse-only, so a decimated
 * pinned-boundary sheet ends up with slivers, a coarse center and a fine rim,
 * and long transition triangles. This pass repairs that by redistributing to a
 * uniform target edge length L.
 *
 * Each outer iteration:  split edges > 4/3 L  ->  collapse edges < 4/5 L  ->
 *                        flip for max-min-angle  ->  tangential relax [-> reproject].
 *
 * INTERIOR-ONLY / boundary FROZEN. `pin_mask` (length in_nv, 1 = frozen) plus
 * every topological open-boundary vertex (any face_count==1 edge) is held
 * bit-identical: never moved, never collapsed away, and boundary edges are never
 * split or flipped. This preserves the weld contract — adjacent hierarchical
 * blocks share this boundary polyline and weld it at the next level.
 *
 * TOPOLOGY-SAFE by construction: split adds an interior vertex (can't merge
 * sheets); collapse is gated by the Dey-Edelsbrunner-Guha link condition, a
 * normal-fold guard, and a max-edge-length cap < the ~7 vox inter-wrap
 * clearance (so a collapse can never fuse two wraps); flip is manifold- and
 * fold-guarded; relax is tangential (on-sheet). FAIL-CLOSED: after the loop a
 * combinatorial 2-manifold audit runs, and if it is not manifold (or the face
 * count drifted out of tolerance) the function emits a verbatim manifold copy
 * of the INPUT instead. It never returns a non-manifold mesh.
 *
 * All output arrays are fresh arena allocations. Returns:
 *   0  = remeshed, manifold, within face-count tolerance
 *   1  = fell back to a verbatim copy of the input (audit failed / out of tol)
 *  <0  = hard error (bad args)
 */

typedef enum {
    REMESH_REPROJECT_NONE       = 0,  /* geometry frozen after connectivity ops   */
    REMESH_REPROJECT_TANGENTIAL = 1,  /* DEFAULT: tangential Laplacian relax only  */
    REMESH_REPROJECT_MLS        = 2   /* opt-in: MLS onto the input surface cloud  */
} RemeshReproject;

typedef struct {
    float target_edge_len;      /* <=0 => auto L = sqrt(4A / (sqrt(3)*F)), F=in_nf */
    int   n_iters;              /* outer iterations                                 */
    float split_ratio;          /* split   edges > split_ratio    * L  (4/3)        */
    float collapse_ratio;       /* collapse edges < collapse_ratio * L  (4/5)       */
    int   do_split;             /* split long edges (break up coarse regions)       */
    int   do_collapse;          /* collapse short edges                             */
    int   do_flip;              /* Surazhsky-Gotsman max-min-angle flips            */
    int   do_relax;             /* tangential Laplacian relax                       */
    RemeshReproject reproject;  /* default REMESH_REPROJECT_TANGENTIAL              */
    float relax_lambda;         /* tangential smoothing step                        */
    float max_collapse_len;     /* hard cap on collapsible edge length (vox)        */
    int   split_max_rounds;     /* inner split rounds per iteration                 */
    int   collapse_max_rounds;  /* inner collapse rounds per iteration              */
    int   flip_max_rounds;      /* inner flip rounds per iteration                  */
    float face_count_tol;       /* accept only if |nf_out - F| <= tol * F           */
} RemeshOpts;

/* Fill o with defaults sourced from pipeline_constants.h (REMESH_*). */
void Remesh_default_opts(RemeshOpts *o);

int Remesh_isotropic(Arena_T arena,
                     const float *in_verts, size_t in_nv,
                     const int32_t *in_faces, size_t in_nf,
                     const uint8_t *pin_mask,
                     const RemeshOpts *opts,
                     float **out_verts, size_t *out_nv,
                     int32_t **out_faces, size_t *out_nf);

#endif
