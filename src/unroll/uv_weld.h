#ifndef UV_WELD_INCLUDED
#define UV_WELD_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"
#include "piece_set.h"

/* ============================================================================
 * uv_weld.h -- cross-cube boundary weld on the in-RAM piece set (step4 front
 * half). Adjacent cubes' charts duplicate the same physical surface in a
 * ~4-vox boundary shell; welding those duplicates into shared vertices turns
 * the piece set into ONE connected sheet per wrap, so the banded relax can
 * couple across seams and close the quilting gutters.
 *
 * Pairing gates (conjunction; measured against the 4x21x21 audit):
 *   3D distance < r3d (3.5 vox = half the 7-vox inter-wrap clearance ==>
 *   cross-wrap welds are geometrically impossible), different cube,
 *   |du| < du_gate (6), |dv| < dv_gate (3), both verts uv-mapped (gid >= 0)
 *   and referenced by kept faces. The >=10-vox mis-registration tail FAILS
 *   the |du| gate by design: force-welding a mis-registered seam would tear
 *   the parameterization; those seams stay quilted and are counted.
 *
 * Weld = union-find clusters; the representative gets the member MEAN pos +
 * uv; faces are retargeted to representatives (nv unchanged, so sidecar
 * indices stay aligned; non-reps become degree-0). Degenerate faces from
 * transitive merges are counted -- the raster (|A2| gate) and relax (h2
 * gate) both skip them by construction.
 * ==========================================================================*/

typedef struct UVWeldOpts {
    double r3d;       /* 3D pair gate, vox (def 3.5) */
    double skin;      /* boundary-shell depth from the cube box faces (def 4.0) */
    double du_gate;   /* |u_i - u_j| max, vox (def 6.0) */
    double dv_gate;   /* |v_i - v_j| max, vox (def 3.0) */
    long   chunk;     /* cube edge, vox (def 128) */
    int    verbose;
} UVWeldOpts;

void UVWeldOpts_default(UVWeldOpts *o);

typedef struct UVWeldStats {
    size_t n_boundary;      /* shell verts considered */
    size_t n_pairs;         /* accepted weld pairs */
    size_t rej_cube, rej_gid, rej_uv;   /* rejected candidates by gate */
    size_t n_clusters;      /* weld clusters (>= 2 members) */
    size_t n_merged_verts;  /* face corner references retargeted to a rep */
    size_t n_degen_faces;   /* faces that lost a distinct corner */
    double du_absmed;       /* |du| median over accepted pairs (pre-merge) */
} UVWeldStats;

/* Weld in place (verts/uv mutated at representatives; faces retargeted).
 * Returns 0 (nothing to weld is fine), -1 on hard failure. */
int UVWeld_run(Arena_T arena, PieceSet *ps, const UVWeldOpts *opts,
               UVWeldStats *out);

int UVWeld_selftest(void);

#endif
