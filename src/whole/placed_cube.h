#ifndef PLACED_CUBE_INCLUDED
#define PLACED_CUBE_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"
#include "../flatten/ribbon.h"

/* ============================================================================
 * placed_cube.h -- one cube's unwrap into the pinned global frame, its
 * on-disk record, and the registration finalize step.
 *
 * Pass B (PlacedCube_unwrap, parallel): load the cube's leaf-stage _all.obj,
 * optionally sever fusion handles, Ribbon_run with the calibration pins
 * (emit_global), flag bad inter-wrap links, and write:
 *   <id>_mesh.obj       the EXACT geometry the sidecars index (post-sever --
 *                       severing duplicates verts along cut loops, so the
 *                       original dump OBJ no longer aligns)
 *   <id>_uvphi_raw.f32  [nv*3] float32 (u, v, phi), PRE-registration
 *   <id>_group.i32      [nv] int32 winding-group id (-1 = unmapped)
 *   <id>_facekeep.u8    [nf] 1 = keep (0 = bad fusion link OR a vertex the
 *                       transfer could not map -- uv is meaningless there)
 *   <id>_skin_raw.f32   [n] SkinVert records (28 B): uv_ok vertices within
 *                       skin_dist of the cube's nominal box faces
 *
 * Pass D (PlacedCube_finalize, parallel after registration): apply the
 * cube's (w_k, du) and write <id>_uvphi.f32, <id>_skin.f32 and
 * <id>_placed.obj (verts + registered vt, kept faces only).
 * ==========================================================================*/

/* Boundary-skin vertex: world position + parameterization + winding group.
 * All members 4 bytes -> sizeof == 28, no padding; dumped raw to .f32. */
typedef struct {
    float   p[3];      /* (z,y,x) */
    float   u, v, phi;
    int32_t gid;       /* per-cube winding-group id; -1 = unmapped */
} SkinVert;

/* A cube's registration, in the form finalize applies: per winding GROUP.
 * The radius anchor rounds one integer offset per pair-graph group, so a
 * whole-turn error is CONSTANT PER GROUP -- one correction per cube tears
 * every wrap whose group rounded differently (measured on real 4x5x5 cubes:
 * dphi MAD across a seam ~= 2pi). gid outside [0, n_groups) (or a NULL
 * table) falls back to the cube-level correction. */
typedef struct {
    const int32_t *g_wk;   /* [n_groups] whole-turn correction per group */
    const double  *g_du;   /* [n_groups] residual u shift per group */
    int32_t        n_groups;
    int32_t        wk_cube;  /* fallback correction */
    double         du_cube;
    /* Optional ArcReg per-group u-WARP: on top of the constant du, add a
     * phi-dependent delta(phi_raw), piecewise-linear over g_warp_nk[gid] knots
     * from g_warp_phi0[gid] spaced g_warp_dphi[gid] (clamped-constant outside).
     * g_warp_nk == NULL (or nk[gid] == 0) => no warp. Stored as plain arrays
     * so finalize evaluates it inline with no arc_reg.h dependency. */
    const int32_t *g_warp_nk;     /* [n_groups] knots per group (0 = none) */
    const double  *g_warp_phi0;   /* [n_groups] */
    const double  *g_warp_dphi;   /* [n_groups] */
    const double  *g_warp_delta;  /* [n_groups * g_warp_stride] heights */
    int32_t        g_warp_stride; /* g_warp_delta row stride */
} PlacedReg;

typedef struct {
    size_t nv, nf;          /* post-sever mesh size */
    size_t n_skin;
    size_t n_badface;       /* faces dropped from facekeep */
    size_t n_group_merges;  /* nearest-slice groups merged by kept topology */
    long   loops_cut;       /* severed genus handles */
    double u_min, u_max, v_min, v_max, phi_min, phi_max;  /* over uv_ok verts */
    double spiral_r2;       /* per-cube diagnostic vs the (pinned) line */
    double spiral_a, spiral_b;  /* the line used (echoes the pin, or the
                                 * per-cube fit when unpinned -- calibration
                                 * reads these off the seed cube) */
    int    w_groups;        /* wraps seen in this cube */
    size_t uv_filled;       /* neighbour-filled, excluded from facekeep */
    size_t uv_fallback;     /* verts the transfer never reached */
    double seconds;
    int    status;          /* 0 = ok; -1 = degenerate unwrap (cube skipped,
                             * no files written; never fatal to the grid) */
} PlacedStats;

/* Pass B for one cube. ropts must carry the calibration pins (emit_global=1,
 * pin_orient=1, wrap_spacing, winding_sense, spiral_a/b) -- the driver builds
 * it once from the seed calibration. origin/chunk define the nominal box for
 * the skin. sever != 0 cuts all genus handles first. cut_* are the
 * Ribbon_flag_bad_faces gates (ratio, floor_vox, len_min).
 * Returns 0 (stats->status reports per-cube failure); -1 on I/O error. */
int PlacedCube_unwrap(Arena_T arena, const char *obj_path, const char *out_dir,
                      const char *id, const RibbonOpts *ropts, int sever,
                      const int64_t origin[3], int64_t chunk, double skin_dist,
                      double cut_ratio, double cut_floor, double cut_len,
                      PlacedStats *stats);

/* Pass D for one cube: read the Pass-B record, apply the PER-GROUP
 * registration (per vertex: k = reg table at its group; phi += 2pi*k,
 * u += DeltaU(phi_raw, k) + du -- see CubeReg_deltaU) and write
 * <id>_uvphi.f32 + <id>_skin.f32 (+ <id>_placed.obj when write_obj).
 * out_ranges (nullable) receives the registered extents
 * {u_min,u_max,v_min,v_max,phi_min,phi_max} over all verts.
 * Returns 0 on success, -1 on I/O failure. */
int PlacedCube_finalize(Arena_T arena, const char *out_dir, const char *id,
                        const PlacedReg *reg,
                        double spiral_a, double spiral_b, int write_obj,
                        double out_ranges[6]);

/* Load a skin sidecar (raw != 0 -> <id>_skin_raw.f32, else <id>_skin.f32).
 * Arena-allocated. Returns 0; -1 if the file is missing/short. */
int PlacedCube_load_skin(Arena_T arena, const char *out_dir, const char *id,
                         int raw, SkinVert **out, size_t *out_n);

/* In-process unit tests (synthetic spiral cube end-to-end: files, skin
 * extraction, finalize roundtrip). Writes scratch files under
 * output/_selftest_scroll_whole/. Returns 0 if all pass, else failures. */
int PlacedCube_selftest(void);

#endif
