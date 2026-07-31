#ifndef CUBE_REGISTER_INCLUDED
#define CUBE_REGISTER_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"
#include "placed_cube.h"

/* ============================================================================
 * cube_register.h -- cross-seam registration of one cube's parameterization
 * against its already-placed neighbors, PER WINDING GROUP.
 *
 * Every cube is unwrapped independently into the SAME pinned global frame
 * (ribbon emit_global + spiral/sense pins). The radius anchor rounds ONE
 * integer offset per pair-graph winding group, so where the scroll wobbles
 * off the pinned line each group can land a whole turn off -- and different
 * groups of one cube round DIFFERENTLY (measured: seam dphi MAD ~= 2pi on
 * real 4x5x5 cubes, and the anchor drifts ~1 turn per cube height from axis
 * tilt). Registration therefore solves one (w_k, du) per GROUP from
 * boundary-skin pairs across the ~2-vox trim gap:
 *     w_k[g] = lround(median_{pairs of g}(phi_P - phi_N) / 2pi)
 *     du[g]  = median(u_P - (u_N + DeltaU(phi_N, w_k[g])))
 * Groups with too few pairs fall back to the cube-level medians. The pair
 * gate (default 3.5 vox) sits below half the 7-vox inter-wrap clearance, so
 * a pair can never join two different wraps.
 * ==========================================================================*/

typedef struct {
    PlacedReg tab;             /* per-group table (allocated from tab_arena);
                                * identity when n_groups == 0 */
    size_t  n_pairs;           /* cross-seam skin pairs used (all groups) */
    double  dphi_mad;          /* MAD of (dphi - 2pi*w_k[g]) AFTER the
                                * per-group correction, radians */
    double  du_mad;            /* MAD of the du residuals after correction */
    int     low_conf;          /* 1 = fewer than min_pairs pairs in total;
                                * table is identity (radius anchor stands) */
    int     deferred;          /* driver bookkeeping: resolved in 2nd pass */
    int32_t n_groups_direct;   /* groups solved from their own pairs */
} CubeReg;

/* u correction for a whole-turn winding shift under the pinned spiral map
 * u(phi) = a*phi + b*phi^2/(4pi):
 *     U(phi + 2pi k) - U(phi) = 2pi*a*k + b*k*phi + pi*b*k^2
 * (~= k local circumferences 2pi*r(phi), exactly what re-solving k turns
 * over would have produced under the pinned init). */
static inline double CubeReg_deltaU(double a, double b, double phi, int32_t k)
{
    const double pi = 3.14159265358979323846;
    return 2.0 * pi * a * (double)k + b * (double)k * phi
         + pi * b * (double)k * (double)k;
}

/* correction lookup: the vertex's group entry, or the cube fallback */
static inline void CubeReg_pick(const PlacedReg *r, int32_t gid,
                                int32_t *out_k, double *out_du)
{
    if (r != NULL && r->g_wk != NULL && gid >= 0 && gid < r->n_groups) {
        *out_k = r->g_wk[gid];
        *out_du = r->g_du[gid];
    } else {
        *out_k = r != NULL ? r->wk_cube : 0;
        *out_du = r != NULL ? r->du_cube : 0.0;
    }
}

/* Apply a registration to one skin vertex (u BEFORE phi: DeltaU takes the
 * raw phi). */
static inline SkinVert CubeReg_apply_vert(SkinVert s, const PlacedReg *r,
                                          double a, double b)
{
    const double pi = 3.14159265358979323846;
    int32_t k = 0;
    double du = 0.0;
    CubeReg_pick(r, s.gid, &k, &du);
    SkinVert o = s;
    o.u   = (float)((double)s.u + CubeReg_deltaU(a, b, (double)s.phi, k) + du);
    o.phi = (float)((double)s.phi + 2.0 * pi * (double)k);
    return o;
}

/* Solve the registration of cube N (raw, unregistered skin with group ids)
 * against the pooled skins of its placed neighbors (caller applies THEIR
 * registrations before pooling, via CubeReg_apply_vert). Pairs = nearest
 * pooled vertex within pair_gate of each N skin vertex (KD-tree). Groups
 * with >= min_group_pairs own pairs get their own (w_k, du); the rest fall
 * back to the cube-level medians. Fewer than min_pairs pairs in total ->
 * low_conf = 1 with an identity table (the radius anchor stands). The
 * per-group table is allocated from tab_arena (must outlive the result);
 * scratch is save/restored. Returns 0 always. */
int CubeReg_solve(Arena_T scratch, Arena_T tab_arena,
                  const SkinVert *skin_n, size_t n_n,
                  const SkinVert *pooled_p, size_t n_p,
                  double spiral_a, double spiral_b,
                  double pair_gate, size_t min_pairs, size_t min_group_pairs,
                  CubeReg *out);

/* In-process unit tests (deltaU algebra, per-group (k,du) recovery incl. a
 * split-winding cube, gate and min-pairs behavior). Returns 0 if all pass,
 * else failure count. */
int CubeReg_selftest(void);

#endif
