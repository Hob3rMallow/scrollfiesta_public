/*
 * scroll_whole.c -- whole-grid front half: unwrap every per-cube mesh into ONE
 * pinned global (u,v,phi) frame, then register each cube's integer winding and
 * residual u against its already-placed neighbors in an outgoing spiral order
 * from a mid-shell seed. One cube in memory at a time; no grid weld, no LOD.
 *
 *   Pass A  calibrate    seed cube -> pinned (spiral_a, spiral_b, sense)
 *   Pass B  unwrap       ALL cubes, parallel (independent once pinned):
 *                        <id>_mesh.obj, <id>_uvphi_raw.f32, <id>_facekeep.u8,
 *                        <id>_skin_raw.f32   (see placed_cube.h)
 *   Pass C  register     sequential spiral order over boundary skins only:
 *                        (w_k, du) per cube (see cube_register.h)
 *   Pass D  finalize     parallel: <id>_uvphi.f32, <id>_skin.f32,
 *                        <id>_placed.obj (world verts + registered vt)
 *   placed_index.json    calibration + per-cube stats + registration
 *
 * Usage:
 *   scroll_whole <dump_dir> <out_dir>
 *       [--leaf-stage step12_final] [--chunk 128]
 *       [--axis-point 0 3405 2878] [--axis-dir 1 0 0] [--wrap-spacing F]
 *       [--pair-gate 3.5] [--skin 4.0] [--min-pairs 24]
 *       [--cut-ratio 4] [--cut-floor 40] [--cut-len 0]
 *       [--seed-id z#####_y#####_x#####] [--max-concurrent 32] [--limit N]
 *       [--no-sever] [--no-obj] [--skip-existing]
 *   scroll_whole <placed_dir> --audit     cross-seam residual audit -> audit.json
 *   scroll_whole --selftest
 *
 * <dump_dir> is a grid_pipeline dump root (one z*_y*_x* dir per cube with
 * <id>_<leaf-stage>/<id>_<leaf-stage>_all.obj in world coords).
 *
 * NOTE: no TRY/EXCEPT here -- except.h's frame stack is a plain global and
 * Pass B/D run under OpenMP; an arena OOM raise simply aborts the tool.
 */
#include "../common/ves_platform.h"

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#ifdef _MSC_VER
#include <windows.h>     /* audit's *_skin.f32 scan */
#else
#include <dirent.h>
#endif

#include "../common/arena.h"
#include "../common/kdtree.h"
#include "../common/maxflow.h"
#include "../common/obj_io.h"
#include "../flatten/ribbon.h"
#include "../whole/cube_schedule.h"
#include "../whole/placed_cube.h"
#include "../whole/cube_register.h"
#include "../whole/group_graph.h"
#include "../whole/weld_track.h"
#include "../whole/weld_build.h"
#include "../whole/arc_reg.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define AUDIT_COLLIDE_GATE 5.0   /* vox: 3D-near pairs must share a wrap */
#define AUDIT_QUALITY_MIN_PAIRS 100
#define AUDIT_QUALITY_MAX_TURNOFF_FRAC 0.05
#define AUDIT_QUALITY_MIN_JOIN_LT2_FRAC 0.50

/* -1: too little evidence, 0: fail, 1: pass. */
static int audit_quality_verdict(size_t n_pairs, size_t n_turnoff,
                                 size_t n_join, size_t n_join_lt2)
{
    if (n_pairs < AUDIT_QUALITY_MIN_PAIRS ||
        n_join < AUDIT_QUALITY_MIN_PAIRS)
        return -1;

    if ((double)n_turnoff / (double)n_pairs >
            AUDIT_QUALITY_MAX_TURNOFF_FRAC ||
        (double)n_join_lt2 / (double)n_join <
            AUDIT_QUALITY_MIN_JOIN_LT2_FRAC)
        return 0;
    return 1;
}

static int audit_quality_selftest(void)
{
    int fail = 0;
    /* PHerc0139 10x10x10 failed/default and forest+uwarp audit signatures. */
    fail |= audit_quality_verdict(50216, 45794, 50216, 4128) != 0;
    fail |= audit_quality_verdict(50216, 1379, 50216, 40113) != 1;
    fail |= audit_quality_verdict(99, 90, 99, 5) != -1;
    fprintf(stderr, "audit_quality_selftest: %s\n", fail ? "FAIL" : "PASS");
    return fail;
}

/* ---- tee logging (stderr + <out>/scroll_whole.log) ------------------------ */

static FILE *g_log = NULL;

/* Escape a UTF-8/path string for use between JSON quotes. */
static int json_escape(char *dst, size_t cap, const char *value)
{
    static const char hex[] = "0123456789abcdef";
    size_t n = 0;
    if (dst == NULL || cap == 0 || value == NULL) return -1;
    for (const unsigned char *p = (const unsigned char *)value; *p; p++) {
        const char *esc = NULL;
        char unicode[7];
        switch (*p) {
        case '"':  esc = "\\\""; break;
        case '\\': esc = "\\\\"; break;
        case '\b': esc = "\\b"; break;
        case '\f': esc = "\\f"; break;
        case '\n': esc = "\\n"; break;
        case '\r': esc = "\\r"; break;
        case '\t': esc = "\\t"; break;
        default:
            if (*p < 0x20) {
                unicode[0] = '\\'; unicode[1] = 'u';
                unicode[2] = '0'; unicode[3] = '0';
                unicode[4] = hex[*p >> 4]; unicode[5] = hex[*p & 15];
                unicode[6] = '\0';
                esc = unicode;
            }
            break;
        }
        if (esc != NULL) {
            size_t k = strlen(esc);
            if (n + k >= cap) return -1;
            memcpy(dst + n, esc, k);
            n += k;
        } else {
            if (n + 1 >= cap) return -1;
            dst[n++] = (char)*p;
        }
    }
    dst[n] = '\0';
    return 0;
}

static int json_escape_selftest(void)
{
    static const char input[] = "C:\\raw\\x\"\n";
    static const char expected[] = "C:\\\\raw\\\\x\\\"\\n";
    char got[sizeof expected + 8];
    int fail = json_escape(got, sizeof got, input) != 0 ||
               strcmp(got, expected) != 0;
    fprintf(stderr, "json_escape_selftest: %s\n", fail ? "FAIL" : "PASS");
    return fail;
}
static void logf_both(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    if (g_log != NULL) {
        va_start(ap, fmt);
        vfprintf(g_log, fmt, ap);
        va_end(ap);
        fflush(g_log);
    }
}

static int file_exists(const char *p)
{
    FILE *f = fopen(p, "rb");
    if (f == NULL) return 0;
    fclose(f);
    return 1;
}

static long file_size(const char *p)
{
    FILE *f = fopen(p, "rb");
    if (f == NULL) return -1;
    if (fseek(f, 0, SEEK_END) != 0) { fclose(f); return -1; }
    long sz = ftell(f);
    fclose(f);
    return sz;
}

static void leaf_obj_path(char *buf, size_t cap, const char *dump_dir,
                          const char *id, const char *stage)
{
    snprintf(buf, cap, "%s/%s/%s_%s/%s_%s_all.obj",
             dump_dir, id, id, stage, id, stage);
}

/* ---- config ---------------------------------------------------------------- */

typedef struct {
    const char *dump_dir, *out_dir, *leaf_stage, *seed_id;
    float  axis_point[3], axis_dir[3];
    double pitch;            /* effective wrap spacing; 0 until auto-calibrated */
    int    pitch_mode;       /* 0 auto, 1 explicitly pinned, 2 legacy index */
    double pair_gate, skin_dist;
    size_t min_pairs, min_group_pairs;
    double cut_ratio, cut_floor, cut_len;
    int64_t chunk;
    int    max_concurrent;
    int    sever, write_obj, skip_existing;
    int    sweeps;           /* consistency (loop-closure) sweeps after the
                              * greedy spiral pass */
    size_t limit;            /* 0 = all cubes */
    /* --reregister (global group-graph re-solve of an existing placed dir) */
    double rr_radius_gate;   /* <=0 => 0.5*pitch */
    double rr_frac_gate;     /* <=0 => 0.5 rad */
    double rr_prior_gate;    /* <=0 => 0.6 turns */
    int    rr_min_edge_pairs;/* <=0 => GroupGraph default */
    int    rr_no_moves;
    int    rr_max_moves;     /* <=0 => GroupGraph default */
    int    rr_raw_component_gauge; /* anchor forest components to raw k chart */
    int    rr_raw_du_gauge;  /* preserve raw continuous-u chart */
    double rr_anchor_weight; /* soft physical-winding anchor in min-cut */
    int    rr_anchor_auto;   /* scale anchor by graph cycle redundancy */
    char   weld_track_path[1024];  /* --weld-track: build the Pass-1 sidecar */
    int    uwarp;            /* --uwarp: ArcReg u-warp measure pass after solve */
    int    uwarp_knots;      /* --uwarp-knots (<=0 => default 5) */
} WholeCfg;

enum {
    WHOLE_PITCH_AUTO = 0,
    WHOLE_PITCH_PINNED = 1,
    WHOLE_PITCH_LEGACY = 2
};

static const char *pitch_mode_name(int mode)
{
    if (mode == WHOLE_PITCH_PINNED) return "pinned";
    if (mode == WHOLE_PITCH_LEGACY) return "legacy";
    return "auto";
}

typedef struct {
    double spiral_a, spiral_b;
    int    sense;
    char   seed_id[48];
} WholeCal;

/* per-cube driver record */
typedef struct {
    PlacedStats st;
    CubeReg     reg;
    double      rgs[6];      /* registered u/v/phi ranges (Pass D) */
    int         active;      /* inside --limit */
    int         skipped;     /* Pass B reused existing files */
} CubeRun;

/* ---- Pass A: calibration ---------------------------------------------------- */

static int calibrate(const WholeCfg *cfg, const CubeNode *nodes, size_t n,
                     const int32_t *order, WholeCal *cal)
{
    memset(cal, 0, sizeof(*cal));
    for (size_t t = 0; t < n; t++) {
        const CubeNode *nd = &nodes[order[t]];
        char path[1024];
        leaf_obj_path(path, sizeof(path), cfg->dump_dir, nd->id, cfg->leaf_stage);

        Arena_T arena = Arena_new();
        float *verts = NULL;
        int32_t *faces = NULL;
        size_t nv = 0, nf = 0;
        int ok = 0;
        if (ObjIO_read(arena, path, &verts, &nv, &faces, &nf) == 0 && nf > 0) {
            RibbonOpts ro;
            RibbonOpts_default(&ro);
            memcpy(ro.axis_point, cfg->axis_point, sizeof(ro.axis_point));
            memcpy(ro.axis_dir, cfg->axis_dir, sizeof(ro.axis_dir));
            ro.wrap_spacing = (float)cfg->pitch;
            ro.fit_ribbon = 0;
            ro.pin_orient = 1;
            RibbonResult R;
            int rr = Ribbon_run(arena, verts, nv, faces, nf, &ro, &R);
            int pitch_ok = cfg->pitch > 0.0
                         || R.pitch_source == RIB_PITCH_ESTIMATED;
            if (rr == 0 && pitch_ok &&
                fabs(R.spiral_b) >= 1.0 && R.spiral_r2 > 0.2) {
                cal->spiral_a = R.spiral_a;
                cal->spiral_b = R.spiral_b;
                cal->sense = R.spiral_b > 0.0 ? 1 : -1;
                snprintf(cal->seed_id, sizeof(cal->seed_id), "%s", nd->id);
                logf_both("[cal] seed %s: spiral a=%.3f b=%.4f r2=%.3f "
                          "sense=%+d (pitch %s %.3f)\n", nd->id, R.spiral_a,
                          R.spiral_b, R.spiral_r2, cal->sense,
                          cfg->pitch > 0.0 ? "pinned" : "estimated",
                          fabs(R.spiral_b));
                ok = 1;
            } else {
                logf_both("[cal] %s unusable (rc=%d b=%.3f r2=%.3f "
                          "pitch_source=%d) -- trying next\n", nd->id, rr,
                          R.spiral_b, R.spiral_r2, R.pitch_source);
            }
        }
        Arena_dispose(&arena);
        if (ok) return 0;
        if (t >= 24) break;   /* two dozen bad candidates = the grid is wrong */
    }
    return -1;
}

/* ---- Pass B ----------------------------------------------------------------- */

static void pass_b(const WholeCfg *cfg, const WholeCal *cal,
                   const CubeNode *nodes, size_t n, CubeRun *runs)
{
    RibbonOpts ro;
    RibbonOpts_default(&ro);
    memcpy(ro.axis_point, cfg->axis_point, sizeof(ro.axis_point));
    memcpy(ro.axis_dir, cfg->axis_dir, sizeof(ro.axis_dir));
    ro.wrap_spacing = (float)cfg->pitch;
    ro.winding_sense = cal->sense;
    ro.spiral_a = cal->spiral_a;
    ro.spiral_b = cal->spiral_b;
    ro.emit_global = 1;
    ro.pin_orient = 1;
    ro.fit_ribbon = 0;

    int nthreads = cfg->max_concurrent;
#ifdef _OPENMP
    if (nthreads > (int)n) nthreads = (int)n > 0 ? (int)n : 1;
    ves_omp_set_threads(nthreads);
#endif
    long done = 0, skipped = 0, failed = 0;
    double t0 = ves_clock_sec();

    int i = 0;   /* MSVC OpenMP 2.0: signed loop var declared outside */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (i = 0; i < (int)n; i++) {
        CubeRun *cr = &runs[i];
        if (!cr->active) continue;
        const CubeNode *nd = &nodes[i];

        if (cfg->skip_existing) {
            char p1[1024], p2[1024], p3[1024], p4[1024];
            snprintf(p1, sizeof(p1), "%s/%s_mesh.obj", cfg->out_dir, nd->id);
            snprintf(p2, sizeof(p2), "%s/%s_uvphi_raw.f32", cfg->out_dir, nd->id);
            snprintf(p3, sizeof(p3), "%s/%s_facekeep.u8", cfg->out_dir, nd->id);
            snprintf(p4, sizeof(p4), "%s/%s_skin_raw.f32", cfg->out_dir, nd->id);
            if (file_exists(p1) && file_exists(p2) && file_exists(p3) &&
                file_exists(p4)) {
                long b2 = file_size(p2), b3 = file_size(p3), b4 = file_size(p4);
                cr->st.nv = b2 > 0 ? (size_t)b2 / (3 * sizeof(float)) : 0;
                cr->st.nf = b3 > 0 ? (size_t)b3 : 0;
                cr->st.n_skin = b4 > 0 ? (size_t)b4 / sizeof(SkinVert) : 0;
                cr->skipped = 1;
#ifdef _OPENMP
#pragma omp critical
#endif
                { skipped++; done++; }
                continue;
            }
        }

        char path[1024];
        leaf_obj_path(path, sizeof(path), cfg->dump_dir, nd->id, cfg->leaf_stage);
        Arena_T arena = Arena_new();
        int64_t origin[3] = { nd->oz, nd->oy, nd->ox };
        int rc = PlacedCube_unwrap(arena, path, cfg->out_dir, nd->id, &ro,
                                   cfg->sever, origin, cfg->chunk,
                                   cfg->skin_dist, cfg->cut_ratio,
                                   cfg->cut_floor, cfg->cut_len, &cr->st);
        Arena_dispose(&arena);
        if (rc != 0) cr->st.status = -2;

#ifdef _OPENMP
#pragma omp critical
#endif
        {
            done++;
            if (cr->st.status != 0) failed++;
            if (cr->st.status != 0 || done % 25 == 0) {
                logf_both("[unwrap %ld/%zu] %s status=%d nv=%zu skin=%zu r2=%.3f "
                          "(%.1fs elapsed)\n", done, n, nd->id, cr->st.status,
                          cr->st.nv, cr->st.n_skin, cr->st.spiral_r2,
                          ves_clock_sec() - t0);
            }
        }
    }
    logf_both("[unwrap] done: %ld cubes (%ld reused, %ld degenerate/failed) in %.1fs\n",
              done, skipped, failed, ves_clock_sec() - t0);
}

/* ---- Pass C ----------------------------------------------------------------- */

/* nonzero if any group (or the fallback) of this cube's table turns */
static int reg_has_turn(const CubeReg *r)
{
    if (r->tab.wk_cube != 0) return 1;
    for (int32_t g = 0; g < r->tab.n_groups; g++)
        if (r->tab.g_wk[g] != 0) return 1;
    return 0;
}

static void pass_c(const WholeCfg *cfg, const WholeCal *cal,
                   const CubeNode *nodes, size_t n, const int32_t *order,
                   const int32_t *comp, CubeRun *runs, Arena_T tab_arena)
{
    Arena_T arena = Arena_new();   /* skins live here for the whole pass */
    double t0 = ves_clock_sec();

    SkinVert **skins = (SkinVert **)ARENA_CALLOC(arena, n, sizeof(SkinVert *));
    size_t *nskin = (size_t *)ARENA_CALLOC(arena, n, sizeof(size_t));
    uint8_t *placed = (uint8_t *)ARENA_CALLOC(arena, n, 1);
    size_t total_skin = 0;
    for (size_t i = 0; i < n; i++) {
        if (!runs[i].active || runs[i].st.status != 0) continue;
        if (PlacedCube_load_skin(arena, cfg->out_dir, nodes[i].id, 1,
                                 &skins[i], &nskin[i]) != 0) {
            skins[i] = NULL;
            nskin[i] = 0;
        }
        total_skin += nskin[i];
    }
    logf_both("[register] %zu skins loaded (%.1f MB)\n", total_skin,
              (double)(total_skin * sizeof(SkinVert)) / 1048576.0);

    int32_t *defer = (int32_t *)ARENA_ALLOC(arena, n * sizeof(int32_t));
    size_t ndefer = 0;
    size_t nreg = 0, nlow = 0;

    for (size_t t = 0; t < n; t++) {
        int32_t i = order[t];
        CubeRun *cr = &runs[i];
        if (!cr->active || cr->st.status != 0) continue;

        /* pool the placed neighbors' skins, with their registrations applied */
        Arena_Mark mark = Arena_save(arena);
        size_t cap = 0;
        for (int e = 0; e < 6; e++) {
            int32_t nb = nodes[i].nbr[e];
            if (nb >= 0 && placed[nb]) cap += nskin[nb];
        }
        SkinVert *pool = (SkinVert *)ARENA_ALLOC(arena,
                                                 (cap ? cap : 1) * sizeof(SkinVert));
        size_t np = 0;
        for (int e = 0; e < 6; e++) {
            int32_t nb = nodes[i].nbr[e];
            if (nb < 0 || !placed[nb]) continue;
            for (size_t s = 0; s < nskin[nb]; s++)
                pool[np++] = CubeReg_apply_vert(skins[nb][s], &runs[nb].reg.tab,
                                                cal->spiral_a, cal->spiral_b);
        }
        CubeReg_solve(arena, tab_arena, skins[i], nskin[i], pool, np,
                      cal->spiral_a, cal->spiral_b,
                      cfg->pair_gate, cfg->min_pairs, cfg->min_group_pairs,
                      &cr->reg);
        Arena_restore(arena, mark);

        placed[i] = 1;
        nreg++;
        if (cr->reg.low_conf) {
            nlow++;
            if (np > 0) defer[ndefer++] = i;   /* neighbors existed, pairs didn't */
        }
        if (reg_has_turn(&cr->reg))
            logf_both("[register] %s: wk_cube=%+d du=%.2f pairs=%zu "
                      "groups=%d(direct %d) mad=%.4f/%.2f (comp %d)\n",
                      nodes[i].id, cr->reg.tab.wk_cube, cr->reg.tab.du_cube,
                      cr->reg.n_pairs, cr->reg.tab.n_groups,
                      cr->reg.n_groups_direct, cr->reg.dphi_mad,
                      cr->reg.du_mad, comp[i]);
    }

    /* deferred second pass: every neighbor is placed now -- retry the
     * low-confidence cubes against the full pooled set */
    size_t fixed = 0;
    for (size_t d = 0; d < ndefer; d++) {
        int32_t i = defer[d];
        Arena_Mark mark = Arena_save(arena);
        size_t cap = 0;
        for (int e = 0; e < 6; e++) {
            int32_t nb = nodes[i].nbr[e];
            if (nb >= 0 && placed[nb]) cap += nskin[nb];
        }
        SkinVert *pool = (SkinVert *)ARENA_ALLOC(arena,
                                                 (cap ? cap : 1) * sizeof(SkinVert));
        size_t np = 0;
        for (int e = 0; e < 6; e++) {
            int32_t nb = nodes[i].nbr[e];
            if (nb < 0 || !placed[nb]) continue;
            for (size_t s = 0; s < nskin[nb]; s++)
                pool[np++] = CubeReg_apply_vert(skins[nb][s], &runs[nb].reg.tab,
                                                cal->spiral_a, cal->spiral_b);
        }
        CubeReg r2;
        CubeReg_solve(arena, tab_arena, skins[i], nskin[i], pool, np,
                      cal->spiral_a, cal->spiral_b,
                      cfg->pair_gate, cfg->min_pairs, cfg->min_group_pairs,
                      &r2);
        Arena_restore(arena, mark);
        if (!r2.low_conf) {
            r2.deferred = 1;
            runs[i].reg = r2;
            fixed++;
            if (reg_has_turn(&r2))
                logf_both("[register/defer] %s: wk_cube=%+d du=%.2f pairs=%zu\n",
                          nodes[i].id, r2.tab.wk_cube, r2.tab.du_cube,
                          r2.n_pairs);
        }
    }

    size_t nturn = 0;
    for (size_t i = 0; i < n; i++)
        if (runs[i].active && reg_has_turn(&runs[i].reg)) nturn++;
    logf_both("[register] %zu cubes registered, %zu low-conf (%zu recovered in "
              "2nd pass), %zu whole-turn corrections, %.1fs\n",
              nreg, nlow, fixed, nturn, ves_clock_sec() - t0);

    /* --- consistency sweeps (loop closure) ---
     * The one-pass spiral chain is GREEDY: each cube sees only the neighbors
     * placed before it, so whole-turn drift accumulates along the chain and
     * the ring around the umbilicus does not close (measured: 9.6% of seam
     * pairs a turn off, worst pair 10 turns, u span inflated by a third).
     * Gauss-Seidel to a fixed point: re-solve every cube against ALL its
     * placed neighbors' CURRENT corrections, repeat until no integer (group
     * or cube) correction changes. Medians make each step robust; the
     * per-cube gauge is mutual so the ensemble cannot run away. */
    for (int sweep = 1; sweep <= cfg->sweeps; sweep++) {
        size_t changes = 0, resolved = 0;
        for (size_t t = 0; t < n; t++) {
            int32_t i = order[t];
            CubeRun *cr = &runs[i];
            if (!cr->active || cr->st.status != 0) continue;
            Arena_Mark mark = Arena_save(arena);
            size_t cap = 0;
            for (int e = 0; e < 6; e++) {
                int32_t nb = nodes[i].nbr[e];
                if (nb >= 0 && placed[nb]) cap += nskin[nb];
            }
            if (cap == 0) { Arena_restore(arena, mark); continue; }
            SkinVert *pool = (SkinVert *)ARENA_ALLOC(arena,
                                                     cap * sizeof(SkinVert));
            size_t np = 0;
            for (int e = 0; e < 6; e++) {
                int32_t nb = nodes[i].nbr[e];
                if (nb < 0 || !placed[nb]) continue;
                for (size_t s = 0; s < nskin[nb]; s++)
                    pool[np++] = CubeReg_apply_vert(skins[nb][s],
                                                    &runs[nb].reg.tab,
                                                    cal->spiral_a,
                                                    cal->spiral_b);
            }
            CubeReg rn;
            CubeReg_solve(arena, tab_arena, skins[i], nskin[i], pool, np,
                          cal->spiral_a, cal->spiral_b,
                          cfg->pair_gate, cfg->min_pairs,
                          cfg->min_group_pairs, &rn);
            Arena_restore(arena, mark);
            if (rn.low_conf) continue;
            /* integer delta between the old and new tables */
            int delta = (rn.tab.wk_cube != cr->reg.tab.wk_cube);
            int32_t gmax = rn.tab.n_groups > cr->reg.tab.n_groups
                           ? rn.tab.n_groups : cr->reg.tab.n_groups;
            for (int32_t g = 0; g < gmax && !delta; g++) {
                int32_t ko = 0, kn = 0;
                double duo = 0.0, dun = 0.0;
                CubeReg_pick(&cr->reg.tab, g, &ko, &duo);
                CubeReg_pick(&rn.tab, g, &kn, &dun);
                if (ko != kn) delta = 1;
            }
            rn.deferred = cr->reg.deferred;
            cr->reg = rn;
            resolved++;
            if (delta) changes++;
        }
        logf_both("[register] sweep %d: %zu/%zu cubes re-solved, %zu integer "
                  "changes\n", sweep, resolved, nreg, changes);
        if (changes == 0) break;
    }
    Arena_dispose(&arena);
}

/* ---- Pass D ----------------------------------------------------------------- */

static void pass_d(const WholeCfg *cfg, const WholeCal *cal,
                   const CubeNode *nodes, size_t n, CubeRun *runs)
{
    double t0 = ves_clock_sec();
    long failed = 0;
    int i = 0;   /* MSVC OpenMP 2.0: signed loop var declared outside */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (i = 0; i < (int)n; i++) {
        CubeRun *cr = &runs[i];
        if (!cr->active || cr->st.status != 0) continue;
        Arena_T arena = Arena_new();
        int rc = PlacedCube_finalize(arena, cfg->out_dir, nodes[i].id,
                                     &cr->reg.tab,
                                     cal->spiral_a, cal->spiral_b,
                                     cfg->write_obj, cr->rgs);
        Arena_dispose(&arena);
        if (rc != 0) {
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                failed++;
                logf_both("[finalize] %s FAILED\n", nodes[i].id);
            }
            cr->st.status = -2;
        }
    }
    logf_both("[finalize] done (%ld failed) in %.1fs\n", failed,
              ves_clock_sec() - t0);
}

/* ---- index ------------------------------------------------------------------ */

static void write_index(const WholeCfg *cfg, const WholeCal *cal,
                        const CubeNode *nodes, size_t n,
                        const int32_t *comp, const CubeRun *runs)
{
    char path[1024], dump_json[6144], leaf_json[768], seed_json[768];
    if (json_escape(dump_json, sizeof dump_json, cfg->dump_dir) != 0 ||
        json_escape(leaf_json, sizeof leaf_json, cfg->leaf_stage) != 0 ||
        json_escape(seed_json, sizeof seed_json, cal->seed_id) != 0) {
        logf_both("ERROR: index string is too long to JSON-escape\n");
        return;
    }
    snprintf(path, sizeof(path), "%s/placed_index.json", cfg->out_dir);
    FILE *f = fopen(path, "w");
    if (f == NULL) {
        logf_both("ERROR: cannot write %s\n", path);
        return;
    }
    size_t n_ok = 0, n_skip = 0, n_low = 0, n_turn = 0;
    double u_lo = 1e300, u_hi = -1e300, v_lo = 1e300, v_hi = -1e300;
    for (size_t i = 0; i < n; i++) {
        if (!runs[i].active) continue;
        if (runs[i].st.status != 0) { n_skip++; continue; }
        n_ok++;
        if (runs[i].reg.low_conf) n_low++;
        if (reg_has_turn(&runs[i].reg)) n_turn++;
        if (runs[i].rgs[0] < u_lo) u_lo = runs[i].rgs[0];
        if (runs[i].rgs[1] > u_hi) u_hi = runs[i].rgs[1];
        if (runs[i].rgs[2] < v_lo) v_lo = runs[i].rgs[2];
        if (runs[i].rgs[3] > v_hi) v_hi = runs[i].rgs[3];
    }
    fprintf(f, "{\n"
            "  \"tool\": \"scroll_whole\",\n"
            "  \"dump_dir\": \"%s\",\n"
            "  \"leaf_stage\": \"%s\",\n"
            "  \"chunk\": %lld,\n"
            "  \"axis_point_zyx\": [%.3f, %.3f, %.3f],\n"
            "  \"axis_dir_zyx\": [%.3f, %.3f, %.3f],\n"
            "  \"pitch\": %.4f,\n"
            "  \"pitch_mode\": \"%s\",\n"
            "  \"pair_gate\": %.2f,\n"
            "  \"skin_dist\": %.2f,\n"
            "  \"calibration\": { \"seed_id\": \"%s\", \"spiral_a\": %.6f, "
            "\"spiral_b\": %.6f, \"sense\": %d },\n"
            "  \"n_cubes\": %zu, \"n_ok\": %zu, \"n_skipped\": %zu,\n"
            "  \"n_low_conf\": %zu, \"n_turn_corrected\": %zu,\n"
            "  \"u_range\": [%.2f, %.2f], \"v_range\": [%.2f, %.2f],\n"
            "  \"cubes\": [\n",
            dump_json, leaf_json, (long long)cfg->chunk,
            (double)cfg->axis_point[0], (double)cfg->axis_point[1],
            (double)cfg->axis_point[2],
            (double)cfg->axis_dir[0], (double)cfg->axis_dir[1],
            (double)cfg->axis_dir[2],
            cfg->pitch, pitch_mode_name(cfg->pitch_mode),
            cfg->pair_gate, cfg->skin_dist,
            seed_json, cal->spiral_a, cal->spiral_b, cal->sense,
            n, n_ok, n_skip, n_low, n_turn, u_lo, u_hi, v_lo, v_hi);
    int first = 1;
    for (size_t i = 0; i < n; i++) {
        if (!runs[i].active) continue;
        const CubeRun *cr = &runs[i];
        fprintf(f, "%s    { \"id\": \"%s\", \"origin\": [%lld, %lld, %lld], "
                "\"status\": %d, \"nv\": %zu, \"nf\": %zu, \"skin\": %zu, "
                "\"bad_faces\": %zu, \"loops_cut\": %ld, \"spiral_r2\": %.4f, "
                "\"w_groups\": %d, \"uv_filled\": %zu, \"uv_fallback\": %zu, "
                "\"wk_cube\": %d, \"du_cube\": %.3f, \"n_groups\": %d, "
                "\"groups_direct\": %d, \"pairs\": %zu, "
                "\"dphi_mad\": %.4f, \"du_mad\": %.3f, "
                "\"low_conf\": %d, \"deferred\": %d, \"flood_comp\": %d, "
                "\"u\": [%.2f, %.2f], \"v\": [%.2f, %.2f], "
                "\"phi\": [%.3f, %.3f], \"sec\": %.2f }",
                first ? "" : ",\n",
                nodes[i].id, (long long)nodes[i].oz, (long long)nodes[i].oy,
                (long long)nodes[i].ox, cr->st.status, cr->st.nv, cr->st.nf,
                cr->st.n_skin, cr->st.n_badface, cr->st.loops_cut,
                cr->st.spiral_r2, cr->st.w_groups, cr->st.uv_filled,
                cr->st.uv_fallback,
                cr->reg.tab.wk_cube, cr->reg.tab.du_cube, cr->reg.tab.n_groups,
                cr->reg.n_groups_direct, cr->reg.n_pairs, cr->reg.dphi_mad,
                cr->reg.du_mad, cr->reg.low_conf, cr->reg.deferred, comp[i],
                cr->rgs[0], cr->rgs[1], cr->rgs[2], cr->rgs[3],
                cr->rgs[4], cr->rgs[5], cr->st.seconds);
        first = 0;
    }
    fprintf(f, "\n  ]\n}\n");
    fclose(f);
    logf_both("[index] %s: %zu ok / %zu skipped / %zu low-conf / %zu turn-corrected; "
              "u span %.0f vox, v span %.0f vox\n", path, n_ok, n_skip, n_low,
              n_turn, u_hi - u_lo, v_hi - v_lo);
}

/* ---- audit ------------------------------------------------------------------ */

static int cmp_dbl_sw(const void *pa, const void *pb)
{
    double a = *(const double *)pa, b = *(const double *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static double med_sw(double *buf, size_t n)
{
    if (n == 0) return 0.0;
    qsort(buf, n, sizeof(double), cmp_dbl_sw);
    return buf[n / 2];
}

/* Cross-seam residual audit over the REGISTERED skins. Every adjacent cube
 * pair contributes nearest-neighbor skin pairs within pair_gate; after a
 * correct registration dphi ~ 0 (mod nothing -- absolute), du ~ 0. */
static int run_audit(const WholeCfg *cfg)
{
    Arena_T arena = Arena_new();
    CubeNode *nodes = NULL;
    size_t n = 0;
    int32_t *order = NULL, *comp = NULL;

    /* enumerate ids from the *_skin.f32 files in the placed dir (adjacency is
     * re-derived from the ids; no dump dir or index needed) */
    typedef struct { char id[48]; } IdRec;
    IdRec *ids = NULL;
    size_t nid = 0, cap = 0;
    {
#ifdef _MSC_VER
        char pattern[1024];
        snprintf(pattern, sizeof(pattern), "%s/*_skin.f32", cfg->out_dir);
        WIN32_FIND_DATAA fd;
        HANDLE h = FindFirstFileA(pattern, &fd);
        if (h != INVALID_HANDLE_VALUE) {
            do {
                char idbuf[48];
                size_t len = strlen(fd.cFileName);
                if (len < 10 || len - 9 >= sizeof(idbuf)) continue;
                memcpy(idbuf, fd.cFileName, len - 9);   /* strip _skin.f32 */
                idbuf[len - 9] = '\0';
                if (idbuf[0] != 'z') continue;
                if (nid >= cap) {
                    size_t nc = cap == 0 ? 64 : cap * 2;
                    IdRec *ni = (IdRec *)ARENA_ALLOC(arena, nc * sizeof(IdRec));
                    if (ids) memcpy(ni, ids, nid * sizeof(IdRec));
                    ids = ni;
                    cap = nc;
                }
                snprintf(ids[nid].id, sizeof(ids[nid].id), "%s", idbuf);
                nid++;
            } while (FindNextFileA(h, &fd));
            FindClose(h);
        }
#else
        DIR *d = opendir(cfg->out_dir);
        if (d != NULL) {
            struct dirent *ent = NULL;
            while ((ent = readdir(d)) != NULL) {
                size_t len = strlen(ent->d_name);
                if (len < 10 || ent->d_name[0] != 'z') continue;
                if (strcmp(ent->d_name + len - 9, "_skin.f32") != 0) continue;
                if (len - 9 >= 48) continue;
                if (nid >= cap) {
                    size_t nc = cap == 0 ? 64 : cap * 2;
                    IdRec *ni = (IdRec *)ARENA_ALLOC(arena, nc * sizeof(IdRec));
                    if (ids) memcpy(ni, ids, nid * sizeof(IdRec));
                    ids = ni;
                    cap = nc;
                }
                memcpy(ids[nid].id, ent->d_name, len - 9);
                ids[nid].id[len - 9] = '\0';
                nid++;
            }
            closedir(d);
        }
#endif
    }
    if (nid == 0) {
        logf_both("audit: no *_skin.f32 in %s\n", cfg->out_dir);
        Arena_dispose(&arena);
        return -1;
    }

    /* nodes from ids */
    nodes = (CubeNode *)ARENA_ALLOC(arena, nid * sizeof(CubeNode));
    n = 0;
    for (size_t i = 0; i < nid; i++) {
        CubeNode nd;
        memset(&nd, 0, sizeof(nd));
        int iz = 0, iy = 0, ix = 0;
        if (sscanf(ids[i].id, "z%d_y%d_x%d", &iz, &iy, &ix) != 3) continue;
        nd.oz = iz; nd.oy = iy; nd.ox = ix;
        snprintf(nd.id, sizeof(nd.id), "%s", ids[i].id);
        nodes[n++] = nd;
    }
    if (CubeSched_link_and_order(arena, nodes, n, cfg->chunk, cfg->axis_point,
                                 cfg->pitch, 0, &order, &comp) != 0) {
        Arena_dispose(&arena);
        return -1;
    }

    /* load registered skins */
    SkinVert **skins = (SkinVert **)ARENA_CALLOC(arena, n, sizeof(SkinVert *));
    size_t *nskin = (size_t *)ARENA_CALLOC(arena, n, sizeof(size_t));
    for (size_t i = 0; i < n; i++)
        PlacedCube_load_skin(arena, cfg->out_dir, nodes[i].id, 0,
                             &skins[i], &nskin[i]);

    char path[1024];
    snprintf(path, sizeof(path), "%s/audit.json", cfg->out_dir);
    FILE *jf = fopen(path, "w");
    if (jf == NULL) { Arena_dispose(&arena); return -1; }
    fprintf(jf, "{\n  \"pairs\": [\n");

    /* residual accumulators */
    size_t all_pairs = 0, all_turnoff = 0, all_collide = 0, npairsets = 0;
    double worst_dphi = 0.0;
    char worst_id[96] = "";
    size_t hist_du[5] = { 0 };   /* <1, <2, <5, <10, >= */
    /* join completeness: |du| below 2 / 6 vox, per seam direction z/y/x --
     * "are physically-adjacent points also adjacent in u after registration" */
    size_t join_lt2[3] = { 0 }, join_lt6[3] = { 0 }, join_n[3] = { 0 };
    /* worst group-level edges by n * |k_resid| (whole-turn residual) */
    enum { WORST_CAP = 20 };
    struct {
        char a[48], b[48];
        size_t n;
        int32_t k_resid;
        double dphi_med, score;
    } worst[WORST_CAP];
    size_t nworst = 0;
    double *gphi = (double *)ARENA_ALLOC(arena, (size_t)4 * 1024 * 1024);
    double *gr = (double *)ARENA_ALLOC(arena, (size_t)4 * 1024 * 1024);
    size_t ng = 0, ngcap = 4 * 1024 * 1024 / sizeof(double);
    double *gu = (double *)ARENA_ALLOC(arena, (size_t)4 * 1024 * 1024);

    int firstpair = 1;
    for (size_t i = 0; i < n; i++) {
        for (int e = 1; e < 6; e += 2) {   /* +z, +y, +x once per pair */
            int32_t j = nodes[i].nbr[e];
            if (j < 0 || nskin[i] == 0 || nskin[(size_t)j] == 0) continue;

            Arena_Mark mark = Arena_save(arena);
            float *pts = (float *)ARENA_ALLOC(arena,
                                              nskin[(size_t)j] * 3 * sizeof(float));
            for (size_t s = 0; s < nskin[(size_t)j]; s++) {
                pts[s * 3 + 0] = skins[(size_t)j][s].p[0];
                pts[s * 3 + 1] = skins[(size_t)j][s].p[1];
                pts[s * 3 + 2] = skins[(size_t)j][s].p[2];
            }
            KDTree_T tree = KDTree_new(arena, pts, nskin[(size_t)j]);
            double *dphi = (double *)ARENA_ALLOC(arena,
                                                 (nskin[i] + 1) * sizeof(double));
            double *du = (double *)ARENA_ALLOC(arena,
                                               (nskin[i] + 1) * sizeof(double));
            double *dv = (double *)ARENA_ALLOC(arena,
                                               (nskin[i] + 1) * sizeof(double));
            size_t np = 0, turnoff = 0, collide = 0;
            float gate2 = (float)(cfg->pair_gate * cfg->pair_gate);
            float cgate2 = (float)(AUDIT_COLLIDE_GATE * AUDIT_COLLIDE_GATE);
            for (size_t s = 0; s < nskin[i]; s++) {
                float d2 = 0.0f;
                size_t bj = KDTree_nearest(tree, skins[i][s].p, &d2);
                if (d2 > cgate2) continue;
                double dp = (double)skins[(size_t)j][bj].phi
                          - (double)skins[i][s].phi;
                if (fabs(dp) > M_PI) collide++;
                if (d2 > gate2) continue;
                dphi[np] = dp;
                du[np] = (double)skins[(size_t)j][bj].u - (double)skins[i][s].u;
                dv[np] = (double)skins[(size_t)j][bj].v - (double)skins[i][s].v;
                if (fabs(dp) > M_PI) turnoff++;
                np++;
            }
            if (np > 0) {
                double *tmp = (double *)ARENA_ALLOC(arena, np * sizeof(double));
                memcpy(tmp, dphi, np * sizeof(double));
                double m_dphi = med_sw(tmp, np);
                memcpy(tmp, du, np * sizeof(double));
                double m_du = med_sw(tmp, np);
                for (size_t q = 0; q < np; q++) tmp[q] = fabs(du[q] - m_du);
                double mad_du = med_sw(tmp, np);
                memcpy(tmp, dv, np * sizeof(double));
                for (size_t q = 0; q < np; q++) tmp[q] = fabs(tmp[q]);
                double m_absdv = med_sw(tmp, np);

                for (size_t q = 0; q < np; q++) {
                    double a = fabs(du[q] - m_du);
                    hist_du[a < 1 ? 0 : a < 2 ? 1 : a < 5 ? 2 : a < 10 ? 3 : 4]++;
                }
                /* join completeness: ABSOLUTE |du| (not median-relative) --
                 * a whole-turn-off seam fails this even when internally
                 * tight around its (wrong) median */
                {
                    int dir = (e - 1) / 2;   /* 0=z 1=y 2=x */
                    for (size_t q = 0; q < np; q++) {
                        double a = fabs(du[q]);
                        join_n[dir]++;
                        if (a < 2.0) join_lt2[dir]++;
                        if (a < 6.0) join_lt6[dir]++;
                    }
                }
                /* worst-edge list by n * |whole-turn residual| */
                {
                    int32_t kres = (int32_t)lround(m_dphi / (2.0 * M_PI));
                    double score = (double)np * fabs((double)kres);
                    if (kres != 0) {
                        size_t slot = nworst;
                        if (nworst < WORST_CAP) {
                            nworst++;
                        } else {
                            slot = 0;
                            for (size_t w2 = 1; w2 < WORST_CAP; w2++)
                                if (worst[w2].score < worst[slot].score)
                                    slot = w2;
                            if (worst[slot].score >= score)
                                slot = WORST_CAP;   /* not better */
                        }
                        if (slot < WORST_CAP) {
                            snprintf(worst[slot].a, sizeof(worst[slot].a),
                                     "%s", nodes[i].id);
                            snprintf(worst[slot].b, sizeof(worst[slot].b),
                                     "%s", nodes[(size_t)j].id);
                            worst[slot].n = np;
                            worst[slot].k_resid = kres;
                            worst[slot].dphi_med = m_dphi;
                            worst[slot].score = score;
                        }
                    }
                }
                if (fabs(m_dphi) > fabs(worst_dphi)) {
                    worst_dphi = m_dphi;
                    snprintf(worst_id, sizeof(worst_id), "%s|%s",
                             nodes[i].id, nodes[(size_t)j].id);
                }
                fprintf(jf, "%s    { \"a\": \"%s\", \"b\": \"%s\", \"n\": %zu, "
                        "\"dphi_med\": %.4f, \"du_med\": %.3f, \"du_mad\": %.3f, "
                        "\"dv_absmed\": %.3f, \"turn_off\": %zu, \"collide\": %zu }",
                        firstpair ? "" : ",\n", nodes[i].id, nodes[(size_t)j].id,
                        np, m_dphi, m_du, mad_du, m_absdv, turnoff, collide);
                firstpair = 0;
                all_pairs += np;
                all_turnoff += turnoff;
                all_collide += collide;
                npairsets++;
            }
            Arena_restore(arena, mark);
        }
        /* global spiral-fit samples (registered skins; subsample) */
        for (size_t s = 0; s < nskin[i]; s += 4) {
            if (ng >= ngcap) break;
            double dy = (double)skins[i][s].p[1] - (double)cfg->axis_point[1];
            double dx = (double)skins[i][s].p[2] - (double)cfg->axis_point[2];
            gphi[ng] = (double)skins[i][s].phi;
            gr[ng] = sqrt(dy * dy + dx * dx);
            gu[ng] = (double)skins[i][s].u;
            ng++;
        }
    }

    /* global fit r ~ a + b*phi/2pi over registered skins */
    double fit_a = 0.0, fit_b = 0.0, fit_r2 = 0.0, corr_u_phi = 0.0;
    if (ng >= 16) {
        double mx = 0.0, my = 0.0;
        for (size_t q = 0; q < ng; q++) { mx += gphi[q] / (2.0 * M_PI); my += gr[q]; }
        mx /= (double)ng;
        my /= (double)ng;
        double sxx = 0.0, sxy = 0.0, syy = 0.0;
        for (size_t q = 0; q < ng; q++) {
            double x = gphi[q] / (2.0 * M_PI) - mx, y = gr[q] - my;
            sxx += x * x; sxy += x * y; syy += y * y;
        }
        fit_b = sxx > 1e-12 ? sxy / sxx : 0.0;
        fit_a = my - fit_b * mx;
        fit_r2 = (sxx > 1e-12 && syy > 1e-12) ? (sxy * sxy) / (sxx * syy) : 0.0;
        double mu = 0.0, mp = 0.0;
        for (size_t q = 0; q < ng; q++) { mu += gu[q]; mp += gphi[q]; }
        mu /= (double)ng;
        mp /= (double)ng;
        double suu = 0.0, sup = 0.0, spp = 0.0;
        for (size_t q = 0; q < ng; q++) {
            double a = gu[q] - mu, b = gphi[q] - mp;
            suu += a * a; sup += a * b; spp += b * b;
        }
        corr_u_phi = (suu > 1e-12 && spp > 1e-12) ? sup / sqrt(suu * spp) : 0.0;
    }

    size_t jl2 = join_lt2[0] + join_lt2[1] + join_lt2[2];
    size_t jl6 = join_lt6[0] + join_lt6[1] + join_lt6[2];
    size_t jn = join_n[0] + join_n[1] + join_n[2];
    int quality = audit_quality_verdict(all_pairs, all_turnoff, jn, jl2);
    const char *quality_status = quality > 0 ? "PASS" :
                                 quality == 0 ? "FAIL" : "INDETERMINATE";
    {
        fprintf(jf, "\n  ],\n"
                "  \"n_pair_sets\": %zu, \"n_pairs\": %zu,\n"
                "  \"turn_off_pairs\": %zu, \"collisions\": %zu,\n"
                "  \"du_resid_hist\": [%zu, %zu, %zu, %zu, %zu],\n"
                "  \"join_completeness\": { \"lt2\": %.4f, \"lt6\": %.4f, "
                "\"n\": %zu,\n"
                "    \"by_dir\": { \"z\": [%.4f, %.4f], \"y\": [%.4f, %.4f], "
                "\"x\": [%.4f, %.4f] } },\n"
                "  \"worst_edges\": [\n",
                npairsets, all_pairs, all_turnoff, all_collide,
                hist_du[0], hist_du[1], hist_du[2], hist_du[3], hist_du[4],
                jn ? (double)jl2 / (double)jn : 0.0,
                jn ? (double)jl6 / (double)jn : 0.0, jn,
                join_n[0] ? (double)join_lt2[0] / (double)join_n[0] : 0.0,
                join_n[0] ? (double)join_lt6[0] / (double)join_n[0] : 0.0,
                join_n[1] ? (double)join_lt2[1] / (double)join_n[1] : 0.0,
                join_n[1] ? (double)join_lt6[1] / (double)join_n[1] : 0.0,
                join_n[2] ? (double)join_lt2[2] / (double)join_n[2] : 0.0,
                join_n[2] ? (double)join_lt6[2] / (double)join_n[2] : 0.0);
        /* worst edges sorted by score desc (tiny list; selection sort) */
        for (size_t w2 = 0; w2 < nworst; w2++) {
            size_t best = w2;
            for (size_t w3 = w2 + 1; w3 < nworst; w3++)
                if (worst[w3].score > worst[best].score) best = w3;
            if (best != w2) {
                struct { char a[48], b[48]; size_t n; int32_t k_resid;
                         double dphi_med, score; } t;
                memcpy(&t, &worst[w2], sizeof(t));
                memcpy(&worst[w2], &worst[best], sizeof(t));
                memcpy(&worst[best], &t, sizeof(t));
            }
            fprintf(jf, "%s    { \"a\": \"%s\", \"b\": \"%s\", \"n\": %zu, "
                    "\"k_resid\": %d, \"dphi_med\": %.3f }",
                    w2 ? ",\n" : "", worst[w2].a, worst[w2].b, worst[w2].n,
                    worst[w2].k_resid, worst[w2].dphi_med);
        }
        fprintf(jf, "\n  ],\n"
                "  \"worst_dphi_med\": %.4f, \"worst_pair\": \"%s\",\n"
                "  \"global_fit\": { \"a\": %.3f, \"b\": %.4f, \"r2\": %.4f, "
                "\"n\": %zu },\n"
                "  \"corr_u_phi\": %.4f,\n"
                "  \"quality_gate\": { \"status\": \"%s\", \"min_pairs\": %d, "
                "\"max_turn_off_frac\": %.4f, \"min_join_lt2_frac\": %.4f }\n}\n",
                worst_dphi, worst_id, fit_a, fit_b, fit_r2, ng, corr_u_phi,
                quality_status, AUDIT_QUALITY_MIN_PAIRS,
                AUDIT_QUALITY_MAX_TURNOFF_FRAC,
                AUDIT_QUALITY_MIN_JOIN_LT2_FRAC);
        logf_both("[audit] join completeness |du|<2: %.2f%%  <6: %.2f%% "
                  "(z %.2f%%/%.2f%%  y %.2f%%/%.2f%%  x %.2f%%/%.2f%%)\n",
                  jn ? 100.0 * (double)jl2 / (double)jn : 0.0,
                  jn ? 100.0 * (double)jl6 / (double)jn : 0.0,
                  join_n[0] ? 100.0 * (double)join_lt2[0] / (double)join_n[0] : 0.0,
                  join_n[0] ? 100.0 * (double)join_lt6[0] / (double)join_n[0] : 0.0,
                  join_n[1] ? 100.0 * (double)join_lt2[1] / (double)join_n[1] : 0.0,
                  join_n[1] ? 100.0 * (double)join_lt6[1] / (double)join_n[1] : 0.0,
                  join_n[2] ? 100.0 * (double)join_lt2[2] / (double)join_n[2] : 0.0,
                  join_n[2] ? 100.0 * (double)join_lt6[2] / (double)join_n[2] : 0.0);
    }
    fclose(jf);

    logf_both("[audit] quality gate: %s (turn-off %.2f%%; limit %.2f%%, "
              "|du|<2 %.2f%%; floor %.2f%%, min pairs %d)\n",
              quality_status,
              all_pairs ? 100.0 * (double)all_turnoff / (double)all_pairs : 0.0,
              100.0 * AUDIT_QUALITY_MAX_TURNOFF_FRAC,
              jn ? 100.0 * (double)jl2 / (double)jn : 0.0,
              100.0 * AUDIT_QUALITY_MIN_JOIN_LT2_FRAC,
              AUDIT_QUALITY_MIN_PAIRS);
    logf_both("[audit] %zu adjacent pair-sets, %zu pairs: turn-off=%zu (%.3f%%) "
              "collisions(<%.0fvox)=%zu\n"
              "[audit] du residual |.|<1/2/5/10/+: %zu/%zu/%zu/%zu/%zu\n"
              "[audit] global spiral fit r ~ a + b*phi/2pi: a=%.2f b=%.3f "
              "r2=%.4f (n=%zu)  corr(u,phi)=%.4f\n"
              "[audit] worst pair dphi_med=%.3f rad (%s)\n"
              "[audit] wrote %s\n",
              npairsets, all_pairs, all_turnoff,
              all_pairs ? 100.0 * (double)all_turnoff / (double)all_pairs : 0.0,
              AUDIT_COLLIDE_GATE, all_collide,
              hist_du[0], hist_du[1], hist_du[2], hist_du[3], hist_du[4],
              fit_a, fit_b, fit_r2, ng, corr_u_phi, worst_dphi, worst_id, path);
    Arena_dispose(&arena);
    return quality == 0 ? -1 : 0;
}

/* ---- reregister -------------------------------------------------------------
 * Global group-graph re-solve of an EXISTING placed dir: reload the RAW
 * skins, build the (cube, winding-group) constraint graph with per-
 * (gidA,gidB) collision-filtered edges, solve k/du globally (max-confidence
 * spanning forest + radius gauge + exact collective-shift min-cut moves),
 * then re-finalize every cube and rewrite placed_index.json (old one kept
 * as .bak). Pass B artifacts (_mesh/_uvphi_raw/_group/_facekeep/_skin_raw)
 * are never touched -- re-runs always restart from raw truth. */

/* crude key scan of the index header (writer emits one key per pattern) */
static const char *rr_jfind(const char *s, const char *key)
{
    char pat[128];
    snprintf(pat, sizeof(pat), "\"%s\"", key);
    const char *p = strstr(s, pat);
    if (p == NULL) return NULL;
    p = strchr(p + strlen(pat), ':');
    return p != NULL ? p + 1 : NULL;
}

static int rr_hex(char c)
{
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return -1;
}

static int rr_jstr(const char *s, const char *key, char *out, size_t cap)
{
    const char *p = rr_jfind(s, key);
    if (p == NULL || out == NULL || cap == 0) return -1;
    p = strchr(p, '"');
    if (p == NULL) return -1;
    p++;
    size_t n = 0;
    while (*p != '\0' && *p != '"') {
        unsigned char c = (unsigned char)*p++;
        if (c == '\\') {
            c = (unsigned char)*p++;
            switch (c) {
            case '"': case '\\': case '/': break;
            case 'b': c = '\b'; break;
            case 'f': c = '\f'; break;
            case 'n': c = '\n'; break;
            case 'r': c = '\r'; break;
            case 't': c = '\t'; break;
            case 'u': {
                if (p[0] == '\0' || p[1] == '\0' ||
                    p[2] == '\0' || p[3] == '\0') return -1;
                int h0 = rr_hex(p[0]), h1 = rr_hex(p[1]);
                int h2 = rr_hex(p[2]), h3 = rr_hex(p[3]);
                if (h0 < 0 || h1 < 0 || h2 < 0 || h3 < 0) return -1;
                unsigned code = (unsigned)((h0 << 12) | (h1 << 8) |
                                           (h2 << 4) | h3);
                p += 4;
                if (code > 0x7f) return -1; /* writer emits UTF-8 directly */
                c = (unsigned char)code;
                break;
            }
            default: return -1;
            }
        }
        if (n + 1 >= cap) return -1;
        out[n++] = (char)c;
    }
    if (*p != '"') return -1;
    out[n] = '\0';
    return 0;
}

static int rr_jstr_selftest(void)
{
    static const char doc[] = "{\"path\":\"C:\\\\raw\\\\x\\\"\\n\"}";
    static const char expected[] = "C:\\raw\\x\"\n";
    char got[64];
    int fail = rr_jstr(doc, "path", got, sizeof got) != 0 ||
               strcmp(got, expected) != 0;
    fprintf(stderr, "rr_jstr_selftest: %s\n", fail ? "FAIL" : "PASS");
    return fail;
}

static int run_reregister(const WholeCfg *cfg_in)
{
    WholeCfg cfg = *cfg_in;
    Arena_T arena = Arena_new();
    double t_all = ves_clock_sec();

    /* ---- calibration from placed_index.json (CLI already seeded cfg; the
     * index values win for the frame parameters unless absent) ---- */
    WholeCal cal;
    memset(&cal, 0, sizeof(cal));
    char dump_dir[1024] = "";
    char leaf_stage[128] = "";
    {
        char path[1024];
        snprintf(path, sizeof(path), "%s/placed_index.json", cfg.out_dir);
        FILE *f = fopen(path, "rb");
        if (f == NULL) {
            logf_both("reregister: no placed_index.json under %s\n",
                      cfg.out_dir);
            Arena_dispose(&arena);
            return -1;
        }
        char buf[8192];
        size_t nb = fread(buf, 1, sizeof(buf) - 1, f);
        fclose(f);
        buf[nb] = '\0';
        const char *p = NULL;
        double a = 0, b = 0, c = 0;
        if ((p = rr_jfind(buf, "axis_point_zyx")) != NULL
            && sscanf(p, " [ %lf , %lf , %lf", &a, &b, &c) == 3) {
            cfg.axis_point[0] = (float)a;
            cfg.axis_point[1] = (float)b;
            cfg.axis_point[2] = (float)c;
        }
        if ((p = rr_jfind(buf, "pitch")) != NULL && sscanf(p, " %lf", &a) == 1)
            cfg.pitch = a;
        {
            char pmode[32] = "";
            if (rr_jstr(buf, "pitch_mode", pmode, sizeof(pmode)) == 0) {
                if (strcmp(pmode, "pinned") == 0)
                    cfg.pitch_mode = WHOLE_PITCH_PINNED;
                else if (strcmp(pmode, "auto") == 0)
                    cfg.pitch_mode = WHOLE_PITCH_AUTO;
                else
                    cfg.pitch_mode = WHOLE_PITCH_LEGACY;
            } else if (cfg.pitch > 0.0) {
                cfg.pitch_mode = WHOLE_PITCH_LEGACY;
            }
        }
        if ((p = rr_jfind(buf, "pair_gate")) != NULL
            && sscanf(p, " %lf", &a) == 1)
            cfg.pair_gate = a;
        if ((p = rr_jfind(buf, "skin_dist")) != NULL
            && sscanf(p, " %lf", &a) == 1)
            cfg.skin_dist = a;
        if ((p = rr_jfind(buf, "chunk")) != NULL && sscanf(p, " %lf", &a) == 1)
            cfg.chunk = (int64_t)a;
        if ((p = rr_jfind(buf, "spiral_a")) != NULL
            && sscanf(p, " %lf", &a) == 1)
            cal.spiral_a = a;
        if ((p = rr_jfind(buf, "spiral_b")) != NULL
            && sscanf(p, " %lf", &a) == 1)
            cal.spiral_b = a;
        if ((p = rr_jfind(buf, "sense")) != NULL && sscanf(p, " %lf", &a) == 1)
            cal.sense = (int)a;
        rr_jstr(buf, "seed_id", cal.seed_id, sizeof(cal.seed_id));
        rr_jstr(buf, "dump_dir", dump_dir, sizeof(dump_dir));
        rr_jstr(buf, "leaf_stage", leaf_stage, sizeof(leaf_stage));
        if (cal.spiral_b == 0.0) {
            logf_both("reregister: index has no calibration "
                      "(spiral_a/spiral_b) -- cannot re-solve\n");
            Arena_dispose(&arena);
            return -1;
        }
        if (cfg.pitch <= 0.0) {
            cfg.pitch = fabs(cal.spiral_b);
            cfg.pitch_mode = WHOLE_PITCH_LEGACY;
        }
    }

    /* ---- enumerate ids from *_skin_raw.f32 ---- */
    typedef struct { char id[48]; } RrId;
    RrId *ids = NULL;
    size_t nid = 0, cap = 0;
    {
        const char *suffix = "_skin_raw.f32";
        const size_t slen = strlen(suffix);
#ifdef _MSC_VER
        char pattern[1024];
        snprintf(pattern, sizeof(pattern), "%s/*%s", cfg.out_dir, suffix);
        WIN32_FIND_DATAA fd;
        HANDLE h = FindFirstFileA(pattern, &fd);
        if (h != INVALID_HANDLE_VALUE) {
            do {
                size_t len = strlen(fd.cFileName);
                if (len <= slen || fd.cFileName[0] != 'z') continue;
                if (len - slen >= sizeof(ids[0].id)) continue;
                if (nid >= cap) {
                    size_t nc = cap == 0 ? 64 : cap * 2;
                    RrId *ni = (RrId *)ARENA_ALLOC(arena, nc * sizeof(RrId));
                    if (ids) memcpy(ni, ids, nid * sizeof(RrId));
                    ids = ni;
                    cap = nc;
                }
                memcpy(ids[nid].id, fd.cFileName, len - slen);
                ids[nid].id[len - slen] = '\0';
                nid++;
            } while (FindNextFileA(h, &fd));
            FindClose(h);
        }
#else
        DIR *d = opendir(cfg.out_dir);
        if (d != NULL) {
            struct dirent *ent = NULL;
            while ((ent = readdir(d)) != NULL) {
                size_t len = strlen(ent->d_name);
                if (len <= slen || ent->d_name[0] != 'z') continue;
                if (strcmp(ent->d_name + len - slen, suffix) != 0) continue;
                if (len - slen >= sizeof(ids[0].id)) continue;
                if (nid >= cap) {
                    size_t nc = cap == 0 ? 64 : cap * 2;
                    RrId *ni = (RrId *)ARENA_ALLOC(arena, nc * sizeof(RrId));
                    if (ids) memcpy(ni, ids, nid * sizeof(RrId));
                    ids = ni;
                    cap = nc;
                }
                memcpy(ids[nid].id, ent->d_name, len - slen);
                ids[nid].id[len - slen] = '\0';
                nid++;
            }
            closedir(d);
        }
#endif
    }
    if (nid == 0) {
        logf_both("reregister: no *_skin_raw.f32 under %s\n", cfg.out_dir);
        Arena_dispose(&arena);
        return -1;
    }

    /* ---- nodes + adjacency ---- */
    CubeNode *nodes = (CubeNode *)ARENA_ALLOC(arena, nid * sizeof(CubeNode));
    size_t n = 0;
    for (size_t i = 0; i < nid; i++) {
        CubeNode nd;
        memset(&nd, 0, sizeof(nd));
        int iz = 0, iy = 0, ix = 0;
        if (sscanf(ids[i].id, "z%d_y%d_x%d", &iz, &iy, &ix) != 3) continue;
        nd.oz = iz;
        nd.oy = iy;
        nd.ox = ix;
        snprintf(nd.id, sizeof(nd.id), "%s", ids[i].id);
        nodes[n++] = nd;
    }
    int32_t *order = NULL, *comp = NULL;
    if (CubeSched_link_and_order(arena, nodes, n, cfg.chunk, cfg.axis_point,
                                 cfg.pitch, 0, &order, &comp) != 0) {
        Arena_dispose(&arena);
        return -1;
    }
    logf_both("[reregister] %zu cubes; pitch=%.3f (%s) pair_gate=%.2f "
              "spiral a=%.6f b=%.6f sense=%d\n", n, cfg.pitch,
              pitch_mode_name(cfg.pitch_mode), cfg.pair_gate, cal.spiral_a,
              cal.spiral_b, cal.sense);

    /* ---- raw skins ---- */
    SkinVert **skins = (SkinVert **)ARENA_CALLOC(arena, n,
                                                 sizeof(SkinVert *));
    size_t *nskin = (size_t *)ARENA_CALLOC(arena, n, sizeof(size_t));
    for (size_t i = 0; i < n; i++)
        PlacedCube_load_skin(arena, cfg.out_dir, nodes[i].id, 1,
                             &skins[i], &nskin[i]);

    /* ---- Pass-1 tracked-weld correspondences (standalone mode) ---- */
    if (cfg.weld_track_path[0] != '\0') {
        WeldBuildOpts wbo;
        WeldBuildOpts_default(&wbo);
        wbo.pair_gate = cfg.pair_gate;
        if (cfg.rr_radius_gate > 0.0) wbo.radius_gate = cfg.rr_radius_gate;
        if (cfg.rr_frac_gate > 0.0) wbo.frac_gate = cfg.rr_frac_gate;
        wbo.verbose = 1;
        WeldTrack wt;
        double tw = ves_clock_sec();
        if (WeldBuild_pass1(arena, nodes, n, skins, nskin, cfg.axis_point,
                            cal.spiral_a, cal.spiral_b, cfg.pitch, &wbo,
                            &wt) != 0) {
            logf_both("weld-track: build failed\n");
            Arena_dispose(&arena);
            return -1;
        }
        if (WeldTrack_write(&wt, cfg.weld_track_path) != 0)
            logf_both("weld-track: WARN write failed (%s)\n",
                      cfg.weld_track_path);
        logf_both("[weld-track] %zu correspondences -> %s (%.1fs)\n",
                  wt.n_corr, cfg.weld_track_path, ves_clock_sec() - tw);
        Arena_dispose(&arena);
        return 0;
    }

    /* ---- build + solve the group graph ---- */
    GroupGraphOpts go;
    GroupGraphOpts_default(&go);
    go.pair_gate = cfg.pair_gate;
    if (cfg.rr_radius_gate > 0.0) go.radius_gate = cfg.rr_radius_gate;
    if (cfg.rr_frac_gate > 0.0) go.frac_gate = cfg.rr_frac_gate;
    if (cfg.rr_prior_gate > 0.0) go.prior_gate = cfg.rr_prior_gate;
    if (cfg.rr_min_edge_pairs > 0) go.min_edge_pairs = cfg.rr_min_edge_pairs;
    if (cfg.rr_no_moves) go.do_moves = 0;
    if (cfg.rr_max_moves > 0) go.max_moves = cfg.rr_max_moves;
    go.raw_component_gauge = cfg.rr_raw_component_gauge;
    go.raw_du_gauge = cfg.rr_raw_du_gauge;
    if (cfg.rr_anchor_weight > 0.0)
        go.anchor_weight = cfg.rr_anchor_weight;
    if (cfg.rr_anchor_auto) go.anchor_redundancy_ref = 0.733;
    go.verbose = 1;

    GroupGraph gg;
    double t0 = ves_clock_sec();
    if (GroupGraph_build(arena, nodes, n, skins, nskin, cfg.axis_point,
                         cal.spiral_a, cal.spiral_b, cfg.pitch, &go,
                         &gg) != 0) {
        logf_both("reregister: graph build failed\n");
        Arena_dispose(&arena);
        return -1;
    }
    logf_both("[reregister] graph: %zu nodes, %zu edges (pairs used=%zu, "
              "rej radius/frac=%zu/%zu) (%.1fs)\n", gg.n_nodes, gg.n_edges,
              gg.pairs_used, gg.pairs_rej_radius, gg.pairs_rej_frac,
              ves_clock_sec() - t0);
    t0 = ves_clock_sec();
    GroupGraph_solve(arena, &gg, &go);
    logf_both("[reregister] solve (forest; k-gauge=%s, du-gauge=%s, "
              "anchor=%.6g): "
              "%d component(s), energy=%.1f, "
              "frustrated=%zu edges, moves=%d (%.1fs)\n",
              cfg.rr_raw_component_gauge ? "raw" : "radius",
              cfg.rr_raw_du_gauge ? "raw" : "seam",
              gg.anchor_weight_effective,
              gg.n_comp, gg.energy, gg.n_frustrated, gg.moves_applied,
              ves_clock_sec() - t0);

    /* ---- per-cube regs from the graph ---- */
    PlacedReg *regs = (PlacedReg *)ARENA_CALLOC(arena, n, sizeof(PlacedReg));
    for (size_t i = 0; i < n; i++)
        GroupGraph_cube_reg(arena, &gg, i, &regs[i]);


    /* ---- Gauss-Seidel polish from the graph seed. The graph fixes the
     * GLOBAL gauge (branch cuts, component turns); the sweeps then refine
     * per-group du with the proven per-cube median machinery and out-vote
     * any small component the radius prior mis-gauged -- from a correct
     * basin, coordinate descent defends the RIGHT majority. Integer changes
     * are counted for the convergence stop, exactly like pass_c. ---- */
    Arena_T polish_tabs = NULL;
    if (cfg.sweeps > 0 && !cfg.rr_raw_component_gauge) {
        /* tables in their OWN arena so the per-cube scratch restore can
         * never free a live table (the arena-aliasing rule: never span
         * another owner's allocations with your save/restore). Disposed at
         * the end of the run -- regs[] points into it until then. */
        Arena_T tab_arena = polish_tabs = Arena_new();
        t0 = ves_clock_sec();
        int sweep = 0;
        for (sweep = 0; sweep < cfg.sweeps; sweep++) {
            size_t changes = 0;
            for (size_t oi = 0; oi < n; oi++) {
                size_t i = (size_t)order[oi];
                if (nskin[i] == 0) continue;
                Arena_Mark mark = Arena_save(arena);
                size_t npool = 0;
                for (int e = 0; e < 6; e++) {
                    int32_t j = nodes[i].nbr[e];
                    if (j >= 0) npool += nskin[(size_t)j];
                }
                if (npool == 0) { Arena_restore(arena, mark); continue; }
                SkinVert *pool = (SkinVert *)ARENA_ALLOC(
                    arena, npool * sizeof(SkinVert));
                size_t at = 0;
                for (int e = 0; e < 6; e++) {
                    int32_t j = nodes[i].nbr[e];
                    if (j < 0) continue;
                    for (size_t s = 0; s < nskin[(size_t)j]; s++)
                        pool[at++] = CubeReg_apply_vert(
                            skins[(size_t)j][s], &regs[(size_t)j],
                            cal.spiral_a, cal.spiral_b);
                }
                CubeReg cr;
                memset(&cr, 0, sizeof(cr));
                CubeReg_solve(arena, tab_arena, skins[i], nskin[i], pool, at,
                              cal.spiral_a, cal.spiral_b, cfg.pair_gate,
                              cfg.min_pairs, cfg.min_group_pairs, &cr);
                Arena_restore(arena, mark);
                if (!cr.low_conf) {
                    int diff = cr.tab.wk_cube != regs[i].wk_cube;
                    int32_t gmax = cr.tab.n_groups < regs[i].n_groups
                                 ? cr.tab.n_groups : regs[i].n_groups;
                    for (int32_t g2 = 0; g2 < gmax && !diff; g2++)
                        if (cr.tab.g_wk[g2] != regs[i].g_wk[g2]) diff = 1;
                    if (diff) changes++;
                    regs[i] = cr.tab;   /* lives in tab_arena */
                }
            }
            logf_both("[reregister] polish sweep %d: %zu integer changes\n",
                      sweep + 1, changes);
            if (changes == 0) break;
        }
        logf_both("[reregister] polish done in %.1fs\n",
                  ves_clock_sec() - t0);
    } else if (cfg.sweeps > 0 && cfg.rr_raw_component_gauge) {
        logf_both("[reregister] polish skipped: raw component gauge is "
                  "incompatible with unconstrained per-cube integer sweeps\n");
    }

    /* ---- ArcReg u-warp (quilting fix): jointly fit a smooth phi-warp from the
     * WTRK correspondences after the discrete forest + scalar polish. The knots
     * are attached to regs[] so finalize applies u += delta(phi_raw). Corrs are
     * built from the same raw skins (cnodes indexed like skins/regs, so the WTRK
     * cube index IS the local index). ---- */
    if (cfg.uwarp) {
        double tu = ves_clock_sec();
        WeldBuildOpts wbo; WeldBuildOpts_default(&wbo);
        wbo.pair_gate = cfg.pair_gate;
        if (cfg.rr_radius_gate > 0.0) wbo.radius_gate = cfg.rr_radius_gate;
        if (cfg.rr_frac_gate > 0.0) wbo.frac_gate = cfg.rr_frac_gate;
        WeldTrack wt;
        if (WeldBuild_pass1(arena, nodes, n, skins, nskin, cfg.axis_point,
                            cal.spiral_a, cal.spiral_b, cfg.pitch, &wbo, &wt) == 0
            && wt.n_corr > 0) {
            int32_t *wt2local = (int32_t *)ARENA_ALLOC(arena,
                                    (long)(wt.n_cubes * sizeof(int32_t)));
            for (size_t i = 0; i < wt.n_cubes; i++)
                wt2local[i] = (i < n) ? (int32_t)i : -1;   /* identity order */
            ArcRegOpts ao; ArcRegOpts_default(&ao);
            if (cfg.uwarp_knots > 0) ao.n_knots = cfg.uwarp_knots;
            size_t max_w = wt.n_corr + 1, nw = max_w;
            ArcRegWarp *warps = (ArcRegWarp *)ARENA_ALLOC(arena,
                                    (long)(max_w * sizeof(ArcRegWarp)));
            ArcRegStats ast;
            ArcReg_solve(arena, &wt, wt2local, skins, nskin, regs, n,
                         cal.spiral_a, cal.spiral_b, &ao, warps, &nw, &ast);
            /* per-cube warp arrays, attached to regs (main arena outlives finalize) */
            for (size_t c = 0; c < n; c++) {
                int32_t ng = regs[c].n_groups;
                if (ng <= 0) continue;
                regs[c].g_warp_nk    = (int32_t *)ARENA_CALLOC(arena, ng, sizeof(int32_t));
                regs[c].g_warp_phi0  = (double *)ARENA_CALLOC(arena, ng, sizeof(double));
                regs[c].g_warp_dphi  = (double *)ARENA_CALLOC(arena, ng, sizeof(double));
                regs[c].g_warp_delta = (double *)ARENA_CALLOC(arena,
                                        (long)ng * ARC_REG_MAX_KNOTS, sizeof(double));
                regs[c].g_warp_stride = ARC_REG_MAX_KNOTS;
            }
            for (size_t w = 0; w < nw; w++) {
                const ArcRegWarp *aw = &warps[w];
                if (aw->cube < 0 || (size_t)aw->cube >= n) continue;
                int32_t gid = aw->gid;
                if (gid < 0 || gid >= regs[aw->cube].n_groups
                    || regs[aw->cube].g_warp_nk == NULL) continue;
                ((int32_t *)regs[aw->cube].g_warp_nk)[gid]   = aw->n_knots;
                ((double *)regs[aw->cube].g_warp_phi0)[gid]  = aw->phi0;
                ((double *)regs[aw->cube].g_warp_dphi)[gid]  = aw->dphi;
                for (int kk = 0; kk < aw->n_knots && kk < ARC_REG_MAX_KNOTS; kk++)
                    ((double *)regs[aw->cube].g_warp_delta)
                        [(size_t)gid * ARC_REG_MAX_KNOTS + kk] = aw->delta[kk];
            }
            logf_both("[uwarp] APPLIED: %zu springs, %zu nodes (%zu scalar, "
                      "%zu turn-off gated, %zu likelihood-downweighted); "
                      "mismatch median %.3f -> %.3f vox "
                      "(p95 %.3f -> %.3f); max|warp|=%.2f (%.1fs)\n",
                      ast.n_corr_used, ast.n_nodes, ast.n_scalar_nodes,
                      ast.n_gated, ast.n_downweighted,
                      ast.med_before, ast.med_after,
                      ast.p95_before, ast.p95_after, ast.max_abs_warp,
                      ves_clock_sec() - tu);
        } else {
            logf_both("[uwarp] no correspondences built; skipped\n");
        }
    }

    double (*rgs)[6] = (double (*)[6])ARENA_CALLOC(arena, n, sizeof(*rgs));
    int32_t *fin_ok = (int32_t *)ARENA_CALLOC(arena, n, sizeof(int32_t));
    t0 = ves_clock_sec();
    {
        long failed = 0;
        int i = 0;
#ifdef _OPENMP
        ves_omp_set_threads(cfg.max_concurrent > 0 ? cfg.max_concurrent : 32);
#pragma omp parallel for schedule(dynamic, 1)
#endif
        for (i = 0; i < (int)n; i++) {
            Arena_T fa = Arena_new();
            int rc = PlacedCube_finalize(fa, cfg.out_dir, nodes[i].id,
                                         &regs[i], cal.spiral_a, cal.spiral_b,
                                         cfg.write_obj, rgs[i]);
            Arena_dispose(&fa);
            fin_ok[i] = rc == 0;
            if (rc != 0) {
#ifdef _OPENMP
#pragma omp critical
#endif
                {
                    failed++;
                    logf_both("[reregister] finalize %s FAILED\n",
                              nodes[i].id);
                }
            }
        }
        logf_both("[reregister] finalize done (%ld failed) in %.1fs\n",
                  failed, ves_clock_sec() - t0);
    }

    /* ---- rewrite the index (old one -> .bak) ---- */
    {
        char oldp[1024], bakp[1024], newp[1024];
        snprintf(oldp, sizeof(oldp), "%s/placed_index.json", cfg.out_dir);
        snprintf(bakp, sizeof(bakp), "%s/placed_index.json.bak", cfg.out_dir);
        remove(bakp);
        rename(oldp, bakp);
        snprintf(newp, sizeof(newp), "%s/placed_index.json", cfg.out_dir);
        FILE *f = fopen(newp, "w");
        if (f == NULL) {
            logf_both("reregister: cannot write %s\n", newp);
            Arena_dispose(&arena);
            return -1;
        }
        char dump_json[6144], leaf_json[768], seed_json[768];
        if (json_escape(dump_json, sizeof dump_json, dump_dir) != 0 ||
            json_escape(leaf_json, sizeof leaf_json, leaf_stage) != 0 ||
            json_escape(seed_json, sizeof seed_json, cal.seed_id) != 0) {
            fclose(f);
            logf_both("reregister: index string is too long to JSON-escape\n");
            Arena_dispose(&arena);
            return -1;
        }
        size_t n_ok = 0;
        double u_lo = 1e300, u_hi = -1e300, v_lo = 1e300, v_hi = -1e300;
        for (size_t i = 0; i < n; i++) {
            if (!fin_ok[i]) continue;
            n_ok++;
            if (rgs[i][0] < u_lo) u_lo = rgs[i][0];
            if (rgs[i][1] > u_hi) u_hi = rgs[i][1];
            if (rgs[i][2] < v_lo) v_lo = rgs[i][2];
            if (rgs[i][3] > v_hi) v_hi = rgs[i][3];
        }
        fprintf(f, "{\n"
                "  \"tool\": \"scroll_whole\",\n"
                "  \"dump_dir\": \"%s\",\n"
                "  \"leaf_stage\": \"%s\",\n"
                "  \"chunk\": %lld,\n"
                "  \"axis_point_zyx\": [%.3f, %.3f, %.3f],\n"
                "  \"axis_dir_zyx\": [%.3f, %.3f, %.3f],\n"
                "  \"pitch\": %.4f,\n"
                "  \"pitch_mode\": \"%s\",\n"
                "  \"pair_gate\": %.2f,\n"
                "  \"skin_dist\": %.2f,\n"
                "  \"calibration\": { \"seed_id\": \"%s\", \"spiral_a\": %.6f, "
                "\"spiral_b\": %.6f, \"sense\": %d },\n"
                "  \"reregistered\": 1,\n"
                "  \"rereg\": { \"nodes\": %zu, \"edges\": %zu, "
                "\"components\": %d, \"energy\": %.1f, \"frustrated\": %zu, "
                "\"moves\": %d },\n"
                "  \"n_cubes\": %zu, \"n_ok\": %zu, \"n_skipped\": %zu,\n"
                "  \"n_low_conf\": 0, \"n_turn_corrected\": 0,\n"
                "  \"u_range\": [%.2f, %.2f], \"v_range\": [%.2f, %.2f],\n"
                "  \"cubes\": [\n",
                dump_json, leaf_json, (long long)cfg.chunk,
                (double)cfg.axis_point[0], (double)cfg.axis_point[1],
                (double)cfg.axis_point[2],
                (double)cfg.axis_dir[0], (double)cfg.axis_dir[1],
                (double)cfg.axis_dir[2],
                cfg.pitch, pitch_mode_name(cfg.pitch_mode),
                cfg.pair_gate, cfg.skin_dist,
                seed_json, cal.spiral_a, cal.spiral_b, cal.sense,
                gg.n_nodes, gg.n_edges, gg.n_comp, gg.energy,
                gg.n_frustrated, gg.moves_applied,
                n, n_ok, n - n_ok, u_lo, u_hi, v_lo, v_hi);
        int first = 1;
        for (size_t i = 0; i < n; i++) {
            char sp[1024];
            snprintf(sp, sizeof(sp), "%s/%s_uvphi_raw.f32", cfg.out_dir,
                     nodes[i].id);
            long usz = file_size(sp);
            snprintf(sp, sizeof(sp), "%s/%s_facekeep.u8", cfg.out_dir,
                     nodes[i].id);
            long fsz = file_size(sp);
            fprintf(f, "%s    { \"id\": \"%s\", \"origin\": [%lld, %lld, "
                    "%lld], \"status\": %d, \"nv\": %ld, \"nf\": %ld, "
                    "\"skin\": %zu, \"flood_comp\": %d,\n"
                    "      \"wk_cube\": %d, \"du_cube\": %.3f, "
                    "\"n_groups\": %d,\n"
                    "      \"g_wk\": [",
                    first ? "" : ",\n", nodes[i].id,
                    (long long)nodes[i].oz, (long long)nodes[i].oy,
                    (long long)nodes[i].ox, fin_ok[i] ? 0 : -2,
                    usz > 0 ? usz / 12 : 0, fsz > 0 ? fsz : 0,
                    nskin[i], comp[i], regs[i].wk_cube, regs[i].du_cube,
                    regs[i].n_groups);
            for (int32_t g2 = 0; g2 < regs[i].n_groups; g2++)
                fprintf(f, "%s%d", g2 ? "," : "", regs[i].g_wk[g2]);
            fprintf(f, "], \"g_du\": [");
            for (int32_t g2 = 0; g2 < regs[i].n_groups; g2++)
                fprintf(f, "%s%.3f", g2 ? "," : "", regs[i].g_du[g2]);
            fprintf(f, "],\n"
                    "      \"u\": [%.2f, %.2f], \"v\": [%.2f, %.2f], "
                    "\"phi\": [%.3f, %.3f] }",
                    rgs[i][0], rgs[i][1], rgs[i][2], rgs[i][3],
                    rgs[i][4], rgs[i][5]);
            first = 0;
        }
        fprintf(f, "\n  ]\n}\n");
        fclose(f);
        logf_both("[reregister] index rewritten (%s kept as .bak); "
                  "u span %.0f vox (%.1fs total)\n", newp, u_hi - u_lo,
                  ves_clock_sec() - t_all);
    }

    if (polish_tabs != NULL)
        Arena_dispose(&polish_tabs);
    Arena_dispose(&arena);
    return 0;
}

/* ---- main ------------------------------------------------------------------- */

static void usage(void)
{
    fprintf(stderr,
        "usage: scroll_whole <dump_dir> <out_dir> [options]\n"
        "       scroll_whole <placed_dir> --audit [--pair-gate F] [--chunk N]\n"
        "       scroll_whole <placed_dir> --reregister [--radius-gate F]\n"
        "           [--frac-gate F] [--min-edge-pairs N] [--no-moves]\n"
        "           [--raw-component-gauge|--radius-component-gauge]\n"
        "           [--raw-du-gauge]\n"
        "           [--anchor-weight F|--auto-anchor]\n"
        "           [--no-uwarp] [--uwarp-knots N]\n"
        "           [--max-moves N] [--no-obj] [--max-concurrent N] [--audit]\n"
        "       scroll_whole <placed_dir> --weld-track <out.wtrk> "
        "[--pair-gate F] [--radius-gate F]\n"
        "       scroll_whole --selftest\n"
        "options:\n"
        "  --leaf-stage S      dump stage (default step12_final)\n"
        "  --chunk N           cube size in vox (default 128)\n"
        "  --axis-point z y x  umbilicus axis point (default 0 3405 2878)\n"
        "  --axis-dir z y x    axis direction (default 1 0 0)\n"
        "  --wrap-spacing F    pin radial pitch vox/turn (default: auto-estimate;\n"
        "                      pass 0 explicitly for auto)\n"
        "  --pair-gate F       cross-seam pairing gate vox (default 3.5)\n"
        "  --skin F            boundary-skin depth vox (default 4.0)\n"
        "  --min-pairs N       registration confidence floor (default 24)\n"
        "  --min-group-pairs N per-group correction floor (default 8; groups\n"
        "                      with fewer pairs take the cube-level medians)\n"
        "  --sweeps N          loop-closure consistency sweeps (default 8)\n"
        "  --cut-ratio/--cut-floor/--cut-len   bad-link gates (4 / 40 / 0)\n"
        "  --seed-id ID        force the calibration + flood seed cube\n"
        "  --max-concurrent N  parallel unwraps (default 32)\n"
        "  --limit N           only the first N cubes of the order (debug)\n"
        "  --no-sever          skip per-cube genus severing\n"
        "  --no-obj            skip <id>_placed.obj outputs\n"
        "  --skip-existing     reuse complete Pass-B records (resume)\n");
}

int main(int argc, char **argv)
{
    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0) {
        int f = json_escape_selftest();
        f += rr_jstr_selftest();
        f += audit_quality_selftest();
        f += CubeSched_selftest();
        f += CubeReg_selftest();
        f += PlacedCube_selftest();
        f += Maxflow_selftest();
        f += GroupGraph_selftest();
        f += WeldTrack_selftest();
        f += WeldBuild_selftest();
        f += ArcReg_selftest();
        fprintf(stderr, "scroll_whole --selftest: %s (%d failures)\n",
                f == 0 ? "PASS" : "FAIL", f);
        return f == 0 ? 0 : 1;
    }
    if (argc < 3) {
        usage();
        return 2;
    }

    WholeCfg cfg;
    memset(&cfg, 0, sizeof(cfg));
    cfg.dump_dir = argv[1];
    cfg.out_dir = argv[2];
    cfg.leaf_stage = "step12_final";
    cfg.axis_point[0] = 0.0f;
    cfg.axis_point[1] = 3405.0f;
    cfg.axis_point[2] = 2878.0f;
    cfg.axis_dir[0] = 1.0f;
    cfg.pitch = 0.0;
    cfg.pitch_mode = WHOLE_PITCH_AUTO;
    cfg.pair_gate = 3.5;
    cfg.skin_dist = 4.0;
    cfg.min_pairs = 24;
    cfg.min_group_pairs = 3;   /* dphi medians quantize to whole turns, so
                                * even tiny groups vote robustly; 8 left too
                                * many minority groups on the cube fallback */
    cfg.cut_ratio = 4.0;
    cfg.cut_floor = 40.0;
    cfg.cut_len = 0.0;
    cfg.chunk = 128;
    cfg.max_concurrent = 32;
    cfg.sever = 1;
    cfg.write_obj = 1;
    cfg.sweeps = 8;
    /* Production registration: retain Ribbon's coherent raw global gauge,
     * reject false cross-sheet forest drift with a weak physical anchor, then
     * fit only the continuous seam warp.  The old radius-gauged +
     * unconstrained-polish path remains available explicitly. */
    cfg.rr_max_moves = 200;
    cfg.rr_min_edge_pairs = 3;
    cfg.rr_raw_component_gauge = 1;
    cfg.rr_anchor_weight = 0.025;
    cfg.rr_anchor_auto = 1;
    cfg.uwarp = 1;
    int audit_only = 0;
    int rereg = 0;

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--audit") == 0) {
            audit_only = 1;
            cfg.out_dir = argv[1];
        }
        else if (strcmp(argv[i], "--reregister") == 0) {
            rereg = 1;
            cfg.out_dir = argv[1];
        }
        else if (strcmp(argv[i], "--weld-track") == 0 && i + 1 < argc) {
            /* reuse the reregister loader; run_reregister early-returns after
             * building the sidecar (before the winding solve). */
            rereg = 1;
            cfg.out_dir = argv[1];
            snprintf(cfg.weld_track_path, sizeof(cfg.weld_track_path), "%s",
                     argv[++i]);
        }
        else if (strcmp(argv[i], "--uwarp") == 0)
            cfg.uwarp = 1;
        else if (strcmp(argv[i], "--no-uwarp") == 0)
            cfg.uwarp = 0;
        else if (strcmp(argv[i], "--uwarp-knots") == 0 && i + 1 < argc)
            cfg.uwarp_knots = atoi(argv[++i]);
        else if (strcmp(argv[i], "--radius-gate") == 0 && i + 1 < argc)
            cfg.rr_radius_gate = atof(argv[++i]);
        else if (strcmp(argv[i], "--frac-gate") == 0 && i + 1 < argc)
            cfg.rr_frac_gate = atof(argv[++i]);
        else if (strcmp(argv[i], "--prior-gate") == 0 && i + 1 < argc)
            cfg.rr_prior_gate = atof(argv[++i]);
        else if (strcmp(argv[i], "--min-edge-pairs") == 0 && i + 1 < argc)
            cfg.rr_min_edge_pairs = atoi(argv[++i]);
        else if (strcmp(argv[i], "--no-moves") == 0)
            cfg.rr_no_moves = 1;
        else if (strcmp(argv[i], "--raw-component-gauge") == 0)
            cfg.rr_raw_component_gauge = 1;
        else if (strcmp(argv[i], "--radius-component-gauge") == 0)
            cfg.rr_raw_component_gauge = 0;
        else if (strcmp(argv[i], "--raw-du-gauge") == 0)
            cfg.rr_raw_du_gauge = 1;
        else if (strcmp(argv[i], "--anchor-weight") == 0 && i + 1 < argc) {
            cfg.rr_anchor_weight = atof(argv[++i]);
            cfg.rr_anchor_auto = 0;
            if (cfg.rr_anchor_weight < 0.0) {
                fprintf(stderr, "ERROR: --anchor-weight must be >= 0\n");
                return 2;
            }
        }
        else if (strcmp(argv[i], "--auto-anchor") == 0) cfg.rr_anchor_auto = 1;
        else if (strcmp(argv[i], "--max-moves") == 0 && i + 1 < argc)
            cfg.rr_max_moves = atoi(argv[++i]);
        else if (strcmp(argv[i], "--leaf-stage") == 0 && i + 1 < argc)
            cfg.leaf_stage = argv[++i];
        else if (strcmp(argv[i], "--chunk") == 0 && i + 1 < argc)
            cfg.chunk = atoll(argv[++i]);
        else if (strcmp(argv[i], "--axis-point") == 0 && i + 3 < argc) {
            cfg.axis_point[0] = (float)atof(argv[++i]);
            cfg.axis_point[1] = (float)atof(argv[++i]);
            cfg.axis_point[2] = (float)atof(argv[++i]);
        }
        else if (strcmp(argv[i], "--axis-dir") == 0 && i + 3 < argc) {
            cfg.axis_dir[0] = (float)atof(argv[++i]);
            cfg.axis_dir[1] = (float)atof(argv[++i]);
            cfg.axis_dir[2] = (float)atof(argv[++i]);
        }
        else if (strcmp(argv[i], "--wrap-spacing") == 0 && i + 1 < argc) {
            cfg.pitch = atof(argv[++i]);
            if (cfg.pitch < 0.0) {
                fprintf(stderr, "ERROR: --wrap-spacing must be >= 0 "
                                "(0 means auto)\n");
                return 2;
            }
            if (cfg.pitch > 0.0 && cfg.pitch < 3.0) {
                fprintf(stderr, "ERROR: --wrap-spacing must be 0 (auto) or "
                                ">= 3 vox/turn\n");
                return 2;
            }
            cfg.pitch_mode = cfg.pitch > 0.0
                           ? WHOLE_PITCH_PINNED : WHOLE_PITCH_AUTO;
        }
        else if (strcmp(argv[i], "--pair-gate") == 0 && i + 1 < argc)
            cfg.pair_gate = atof(argv[++i]);
        else if (strcmp(argv[i], "--skin") == 0 && i + 1 < argc)
            cfg.skin_dist = atof(argv[++i]);
        else if (strcmp(argv[i], "--min-pairs") == 0 && i + 1 < argc)
            cfg.min_pairs = (size_t)atoll(argv[++i]);
        else if (strcmp(argv[i], "--min-group-pairs") == 0 && i + 1 < argc)
            cfg.min_group_pairs = (size_t)atoll(argv[++i]);
        else if (strcmp(argv[i], "--sweeps") == 0 && i + 1 < argc)
            cfg.sweeps = atoi(argv[++i]);
        else if (strcmp(argv[i], "--cut-ratio") == 0 && i + 1 < argc)
            cfg.cut_ratio = atof(argv[++i]);
        else if (strcmp(argv[i], "--cut-floor") == 0 && i + 1 < argc)
            cfg.cut_floor = atof(argv[++i]);
        else if (strcmp(argv[i], "--cut-len") == 0 && i + 1 < argc)
            cfg.cut_len = atof(argv[++i]);
        else if (strcmp(argv[i], "--seed-id") == 0 && i + 1 < argc)
            cfg.seed_id = argv[++i];
        else if (strcmp(argv[i], "--max-concurrent") == 0 && i + 1 < argc)
            cfg.max_concurrent = atoi(argv[++i]);
        else if (strcmp(argv[i], "--limit") == 0 && i + 1 < argc)
            cfg.limit = (size_t)atoll(argv[++i]);
        else if (strcmp(argv[i], "--no-sever") == 0)
            cfg.sever = 0;
        else if (strcmp(argv[i], "--no-obj") == 0)
            cfg.write_obj = 0;
        else if (strcmp(argv[i], "--skip-existing") == 0)
            cfg.skip_existing = 1;
        else if (i > 2) {
            fprintf(stderr, "ERROR: unknown arg %s\n", argv[i]);
            return 2;
        }
    }

    if (rereg) {
        char lp[1024];
        snprintf(lp, sizeof(lp), "%s/reregister.log", cfg.out_dir);
        g_log = fopen(lp, "w");
        int rc = run_reregister(&cfg);
        if (rc == 0 && audit_only) {
            /* chain the audit over the freshly registered skins */
            if (g_log) fclose(g_log);
            snprintf(lp, sizeof(lp), "%s/audit.log", cfg.out_dir);
            g_log = fopen(lp, "w");
            rc = run_audit(&cfg);
        }
        if (g_log) fclose(g_log);
        return rc == 0 ? 0 : 1;
    }
    if (audit_only) {
        char lp[1024];
        snprintf(lp, sizeof(lp), "%s/audit.log", cfg.out_dir);
        g_log = fopen(lp, "w");
        int rc = run_audit(&cfg);
        if (g_log) fclose(g_log);
        return rc == 0 ? 0 : 1;
    }

    /* out dir + log */
    {
        char probe[1024];
        snprintf(probe, sizeof(probe), "%s/x", cfg.out_dir);
        ves_ensure_parent_dir(probe);
        char lp[1024];
        snprintf(lp, sizeof(lp), "%s/scroll_whole.log", cfg.out_dir);
        g_log = fopen(lp, "w");
    }
    double t_all = ves_clock_sec();
    if (cfg.pitch > 0.0)
        logf_both("scroll_whole: dump=%s out=%s stage=%s chunk=%lld "
                  "pitch=%.3f (pinned) axis=(%.0f,%.0f,%.0f) gate=%.1f "
                  "skin=%.1f concurrent=%d\n", cfg.dump_dir, cfg.out_dir,
                  cfg.leaf_stage, (long long)cfg.chunk, cfg.pitch,
                  (double)cfg.axis_point[0], (double)cfg.axis_point[1],
                  (double)cfg.axis_point[2], cfg.pair_gate, cfg.skin_dist,
                  cfg.max_concurrent);
    else
        logf_both("scroll_whole: dump=%s out=%s stage=%s chunk=%lld "
                  "pitch=auto axis=(%.0f,%.0f,%.0f) gate=%.1f skin=%.1f "
                  "concurrent=%d\n", cfg.dump_dir, cfg.out_dir,
                  cfg.leaf_stage, (long long)cfg.chunk,
                  (double)cfg.axis_point[0], (double)cfg.axis_point[1],
                  (double)cfg.axis_point[2], cfg.pair_gate, cfg.skin_dist,
                  cfg.max_concurrent);

    Arena_T arena = Arena_new();
    CubeNode *nodes = NULL;
    size_t n = 0;
    int32_t *order = NULL, *comp = NULL;
    if (CubeSched_build(arena, cfg.dump_dir, cfg.leaf_stage, cfg.chunk,
                        cfg.axis_point, cfg.pitch, cfg.seed_id,
                        &nodes, &n, &order, &comp) != 0) {
        logf_both("ERROR: no cubes with %s under %s\n", cfg.leaf_stage,
                  cfg.dump_dir);
        return 1;
    }
    WholeCal cal;
    if (calibrate(&cfg, nodes, n, order, &cal) != 0) {
        logf_both("ERROR: calibration failed on the first candidates"
                  "%s\n", cfg.pitch > 0.0 ? "" :
                  " (auto-pitch had no supported estimate; refusing fallback)");
        return 1;
    }

    if (cfg.pitch <= 0.0) {
        cfg.pitch = fabs(cal.spiral_b);
        logf_both("[pitch] auto-estimated %.4f vox/turn from seed %s; "
                  "rebuilding spiral schedule\n", cfg.pitch, cal.seed_id);
    }

    /* Auto mode initially uses a pitch-free radial flood solely to select a
     * calibration seed. Rebuild the real winding-aware order using the measured
     * pitch (also anchors a pinned run on the seed that actually calibrated). */
    int32_t seed_idx = -1;
    for (size_t i = 0; i < n; i++)
        if (strcmp(nodes[i].id, cal.seed_id) == 0) {
            seed_idx = (int32_t)i;
            break;
        }
    if (CubeSched_link_and_order(arena, nodes, n, cfg.chunk, cfg.axis_point,
                                 cfg.pitch, seed_idx, &order, &comp) != 0) {
        logf_both("ERROR: cannot rebuild schedule with calibrated pitch\n");
        return 1;
    }

    int ncomp = 0;
    for (size_t i = 0; i < n; i++) if (comp[i] + 1 > ncomp) ncomp = comp[i] + 1;
    logf_both("[schedule] %zu cubes, %d flood component(s), seed=%s "
              "r=%.0f pitch=%.4f (%s)\n", n, ncomp, nodes[order[0]].id,
              nodes[order[0]].r, cfg.pitch, pitch_mode_name(cfg.pitch_mode));

    CubeRun *runs = (CubeRun *)ARENA_CALLOC(arena, n, sizeof(CubeRun));
    size_t nact = (cfg.limit > 0 && cfg.limit < n) ? cfg.limit : n;
    for (size_t t = 0; t < nact; t++) runs[order[t]].active = 1;
    if (nact < n)
        logf_both("[schedule] --limit %zu: %zu of %zu cubes active\n",
                  cfg.limit, nact, n);

    pass_b(&cfg, &cal, nodes, n, runs);
    pass_c(&cfg, &cal, nodes, n, order, comp, runs, arena);
    pass_d(&cfg, &cal, nodes, n, runs);
    write_index(&cfg, &cal, nodes, n, comp, runs);

    logf_both("scroll_whole: OK (total %.1fs)\n", ves_clock_sec() - t_all);
    if (g_log) fclose(g_log);
    Arena_dispose(&arena);
    return 0;
}
