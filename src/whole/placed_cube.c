/* placed_cube.c -- per-cube unwrap record + finalize. See placed_cube.h. */
#include "../common/ves_platform.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/obj_io.h"
#include "../common/union_find.h"
#include "../flatten/seam_cut.h"
#include "placed_cube.h"
#include "cube_register.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- tiny binary sidecar I/O (native little-endian x86/x64) --------------- */

static int write_blob(const char *path, const void *data, size_t nbytes)
{
    FILE *f = fopen(path, "wb");
    if (f == NULL) return -1;
    size_t w = nbytes ? fwrite(data, 1, nbytes, f) : 0;
    int bad = (w != nbytes);
    if (fclose(f) != 0) bad = 1;
    return bad ? -1 : 0;
}

static void *read_blob(Arena_T arena, const char *path, size_t *out_nbytes)
{
    *out_nbytes = 0;
    FILE *f = fopen(path, "rb");
    if (f == NULL) return NULL;
    if (fseek(f, 0, SEEK_END) != 0) { fclose(f); return NULL; }
    long sz = ftell(f);
    if (sz < 0) { fclose(f); return NULL; }
    rewind(f);
    void *buf = ARENA_ALLOC(arena, (size_t)sz + 1);
    if (sz > 0 && fread(buf, 1, (size_t)sz, f) != (size_t)sz) {
        fclose(f);
        return NULL;
    }
    fclose(f);
    *out_nbytes = (size_t)sz;
    return buf;
}

static void sidecar_path(char *buf, size_t cap, const char *out_dir,
                         const char *id, const char *suffix)
{
    snprintf(buf, cap, "%s/%s_%s", out_dir, id, suffix);
}

/* ---- Pass B ---------------------------------------------------------------- */

/* Ribbon's vertex->slice transfer labels a vertex with the pair-graph group of
 * its nearest slice sample.  That label is useful at a detached-chart boundary,
 * but it is not topology preserving: two vertices joined by an ordinary kept
 * mesh edge can select samples from different pair-graph components even when
 * their raw (u,phi) values are identical.  Treating those labels as independent
 * registration domains tears the triangle by whole turns during finalize.
 *
 * A kept face is our authoritative statement that its three vertices belong to
 * one continuous chart (bad fusion links and unmapped faces were removed just
 * above the call site).  Merge every group connected by kept topology, then
 * compact the ids so the global graph and its per-group tables see the same
 * domains as the rasterizer. */
static size_t pc_merge_kept_groups(Arena_T arena,
                                   const int32_t *faces, size_t nf,
                                   const uint8_t *keep,
                                   int32_t *group, size_t nv,
                                   int *out_n_groups)
{
    int32_t max_gid = -1;
    for (size_t i = 0; i < nv; i++)
        if (group[i] > max_gid) max_gid = group[i];
    if (max_gid < 0) {
        if (out_n_groups != NULL) *out_n_groups = 0;
        return 0;
    }

    UnionFind uf = UF_new(arena, max_gid + 1);
    uint8_t *present = (uint8_t *)ARENA_CALLOC(arena,
                                                (size_t)max_gid + 1, 1);
    for (size_t i = 0; i < nv; i++)
        if (group[i] >= 0) present[group[i]] = 1;

    for (size_t f = 0; f < nf; f++) {
        if (keep != NULL && !keep[f]) continue;
        int32_t g[3] = {
            group[faces[f * 3 + 0]],
            group[faces[f * 3 + 1]],
            group[faces[f * 3 + 2]]
        };
        if (g[0] >= 0 && g[1] >= 0) uf_union(&uf, g[0], g[1]);
        if (g[1] >= 0 && g[2] >= 0) uf_union(&uf, g[1], g[2]);
        if (g[2] >= 0 && g[0] >= 0) uf_union(&uf, g[2], g[0]);
    }

    int32_t *root_id = (int32_t *)ARENA_ALLOC(
        arena, ((size_t)max_gid + 1) * sizeof(int32_t));
    for (int32_t g = 0; g <= max_gid; g++) root_id[g] = -1;
    int32_t n_before = 0, n_after = 0;
    for (int32_t g = 0; g <= max_gid; g++) {
        if (!present[g]) continue;
        n_before++;
        int32_t r = uf_find(&uf, g);
        if (root_id[r] < 0) root_id[r] = n_after++;
    }
    for (size_t i = 0; i < nv; i++)
        if (group[i] >= 0) group[i] = root_id[uf_find(&uf, group[i])];

    if (out_n_groups != NULL) *out_n_groups = n_after;
    return (size_t)(n_before - n_after);
}

int PlacedCube_unwrap(Arena_T arena, const char *obj_path, const char *out_dir,
                      const char *id, const RibbonOpts *ropts, int sever,
                      const int64_t origin[3], int64_t chunk, double skin_dist,
                      double cut_ratio, double cut_floor, double cut_len,
                      PlacedStats *stats)
{
    assert(arena && obj_path && out_dir && id && ropts && stats);
    memset(stats, 0, sizeof(*stats));
    double t0 = ves_clock_sec();

    float *verts = NULL;
    int32_t *faces = NULL;
    size_t nv = 0, nf = 0;
    if (ObjIO_read(arena, obj_path, &verts, &nv, &faces, &nf) != 0) {
        fprintf(stderr, "placed_cube[%s]: cannot read %s\n", id, obj_path);
        return -1;
    }
    if (nv < 3 || nf < 1) {
        stats->status = -1;   /* empty leaf: skip, never fatal */
        stats->seconds = ves_clock_sec() - t0;
        return 0;
    }

    if (sever) {
        float *sv = NULL;
        int32_t *sf = NULL;
        size_t snv = 0, snf = 0;
        long cut = 0;
        if (SeamCut_sever_short_handles(arena, verts, nv, faces, nf,
                                        0.0 /* all lengths */,
                                        &sv, &snv, &sf, &snf, &cut) == 0) {
            verts = sv; nv = snv; faces = sf; nf = snf;
            stats->loops_cut = cut;
        }
    }

    RibbonResult R;
    if (Ribbon_run(arena, verts, nv, faces, nf, ropts, &R) != 0) {
        stats->status = -1;
        stats->seconds = ves_clock_sec() - t0;
        return 0;
    }
    assert(R.phi != NULL && R.group != NULL);   /* ropts must carry emit_global=1 */

    /* facekeep: drop flagged bad fusion links AND faces touching a vertex with
     * no direct slice-map transfer. Neighbor-filled u is useful for diagnostics
     * and propagation, but it produces constant-u plateaus and zero-area UV
     * triangles if treated as real parameterization geometry. */
    uint8_t *keep = (uint8_t *)ARENA_ALLOC(arena, nf ? nf : 1);
    memset(keep, 1, nf ? nf : 1);
    size_t nbad = 0;
    {
        uint8_t *bad = (uint8_t *)ARENA_CALLOC(arena, nf ? nf : 1, 1);
        size_t nflag = 0;
        Ribbon_flag_bad_faces(verts, nv, faces, nf, R.uv,
                              cut_ratio, cut_floor, cut_len, bad, &nflag);
        for (size_t t = 0; t < nf; t++) {
            int drop = bad[t];
            for (int e = 0; e < 3 && !drop; e++)
                if (!R.uv_ok[faces[t * 3 + e]]) drop = 1;
            if (drop) { keep[t] = 0; nbad++; }
        }
    }
    stats->n_badface = nbad;

    /* Registration groups must respect the topology that survives facekeep.
     * Do this before both skin extraction and sidecar emission so every later
     * consumer observes the same canonical ids. */
    {
        int n_groups = 0;
        stats->n_group_merges = pc_merge_kept_groups(
            arena, faces, nf, keep, R.group, nv, &n_groups);
        R.w_prior_groups = n_groups;
    }

    /* skin: uv_ok verts within skin_dist of the nominal box faces */
    SkinVert *skin = (SkinVert *)ARENA_ALLOC(arena,
                                             (nv ? nv : 1) * sizeof(SkinVert));
    size_t nskin = 0;
    for (size_t i = 0; i < nv; i++) {
        if (!R.uv_ok[i]) continue;
        double bd = 1e300;
        for (int ax = 0; ax < 3; ax++) {
            double lo = (double)verts[i * 3 + ax] - (double)origin[ax];
            double hi = (double)(origin[ax] + chunk) - (double)verts[i * 3 + ax];
            double d = lo < hi ? lo : hi;
            if (d < bd) bd = d;
        }
        if (bd > skin_dist) continue;
        SkinVert *s = &skin[nskin++];
        s->p[0] = verts[i * 3 + 0];
        s->p[1] = verts[i * 3 + 1];
        s->p[2] = verts[i * 3 + 2];
        s->u = R.uv[i * 2 + 0];
        s->v = R.uv[i * 2 + 1];
        s->phi = R.phi[i];
        s->gid = R.group[i];
    }
    stats->n_skin = nskin;

    /* ranges over mapped verts */
    stats->u_min = stats->v_min = stats->phi_min = 1e300;
    stats->u_max = stats->v_max = stats->phi_max = -1e300;
    for (size_t i = 0; i < nv; i++) {
        if (!R.uv_ok[i]) continue;
        double u = (double)R.uv[i * 2 + 0], v = (double)R.uv[i * 2 + 1];
        double ph = (double)R.phi[i];
        if (u < stats->u_min) stats->u_min = u;
        if (u > stats->u_max) stats->u_max = u;
        if (v < stats->v_min) stats->v_min = v;
        if (v > stats->v_max) stats->v_max = v;
        if (ph < stats->phi_min) stats->phi_min = ph;
        if (ph > stats->phi_max) stats->phi_max = ph;
    }
    stats->nv = nv;
    stats->nf = nf;
    stats->spiral_r2 = R.spiral_r2;
    stats->spiral_a = R.spiral_a;
    stats->spiral_b = R.spiral_b;
    stats->w_groups = R.w_groups;
    stats->uv_filled = R.uv_filled;
    stats->uv_fallback = R.uv_fallback;

    /* write the record */
    char path[1024];
    sidecar_path(path, sizeof(path), out_dir, id, "mesh.obj");
    if (ObjIO_write(path, verts, nv, faces, nf) != 0) return -1;

    float *uvphi = (float *)ARENA_ALLOC(arena, (nv ? nv : 1) * 3 * sizeof(float));
    for (size_t i = 0; i < nv; i++) {
        uvphi[i * 3 + 0] = R.uv[i * 2 + 0];
        uvphi[i * 3 + 1] = R.uv[i * 2 + 1];
        uvphi[i * 3 + 2] = R.phi[i];
    }
    sidecar_path(path, sizeof(path), out_dir, id, "uvphi_raw.f32");
    if (write_blob(path, uvphi, nv * 3 * sizeof(float)) != 0) return -1;

    sidecar_path(path, sizeof(path), out_dir, id, "group.i32");
    if (write_blob(path, R.group, nv * sizeof(int32_t)) != 0) return -1;

    sidecar_path(path, sizeof(path), out_dir, id, "facekeep.u8");
    if (write_blob(path, keep, nf) != 0) return -1;

    sidecar_path(path, sizeof(path), out_dir, id, "skin_raw.f32");
    if (write_blob(path, skin, nskin * sizeof(SkinVert)) != 0) return -1;

    stats->seconds = ves_clock_sec() - t0;
    return 0;
}

/* ---- Pass D ---------------------------------------------------------------- */

/* Evaluate the ArcReg per-group u-warp at raw phi (piecewise-linear over the
 * group's knots, clamped-constant outside). 0 when no warp is attached. */
static double pc_warp_eval(const PlacedReg *reg, int32_t gid, double phi)
{
    if (reg->g_warp_nk == NULL || gid < 0 || gid >= reg->n_groups) return 0.0;
    int nk = reg->g_warp_nk[gid];
    if (nk <= 0) return 0.0;
    const double *d = &reg->g_warp_delta[(size_t)gid * (size_t)reg->g_warp_stride];
    if (nk == 1) return d[0];
    double phi0 = reg->g_warp_phi0[gid], dphi = reg->g_warp_dphi[gid];
    if (dphi <= 0.0) return d[0];
    double t = (phi - phi0) / dphi;
    if (t <= 0.0) return d[0];
    if (t >= (double)(nk - 1)) return d[nk - 1];
    int i = (int)floor(t);
    double f = t - (double)i;
    return d[i] * (1.0 - f) + d[i + 1] * f;
}

int PlacedCube_finalize(Arena_T arena, const char *out_dir, const char *id,
                        const PlacedReg *reg,
                        double spiral_a, double spiral_b, int write_obj,
                        double out_ranges[6])
{
    assert(arena && out_dir && id && reg);
    Arena_Mark mark = Arena_save(arena);
    char path[1024];

    size_t nb = 0;
    sidecar_path(path, sizeof(path), out_dir, id, "uvphi_raw.f32");
    float *uvphi = (float *)read_blob(arena, path, &nb);
    if (uvphi == NULL || nb % (3 * sizeof(float)) != 0) {
        Arena_restore(arena, mark);
        return -1;
    }
    size_t nv = nb / (3 * sizeof(float));

    size_t gb = 0;
    sidecar_path(path, sizeof(path), out_dir, id, "group.i32");
    int32_t *grp = (int32_t *)read_blob(arena, path, &gb);
    if (grp == NULL || gb != nv * sizeof(int32_t)) {
        Arena_restore(arena, mark);
        return -1;
    }

    /* register per vertex via its winding group: u first (DeltaU takes the
     * raw phi), then phi */
    double rg[6] = { 1e300, -1e300, 1e300, -1e300, 1e300, -1e300 };
    for (size_t i = 0; i < nv; i++) {
        int32_t k = 0;
        double du = 0.0;
        CubeReg_pick(reg, grp[i], &k, &du);
        double phi = (double)uvphi[i * 3 + 2];
        uvphi[i * 3 + 0] = (float)((double)uvphi[i * 3 + 0]
                                   + CubeReg_deltaU(spiral_a, spiral_b, phi, k)
                                   + du + pc_warp_eval(reg, grp[i], phi));
        uvphi[i * 3 + 2] = (float)(phi + 2.0 * M_PI * (double)k);
        double uu = (double)uvphi[i * 3 + 0], vv = (double)uvphi[i * 3 + 1];
        double pp = (double)uvphi[i * 3 + 2];
        if (uu < rg[0]) rg[0] = uu;
        if (uu > rg[1]) rg[1] = uu;
        if (vv < rg[2]) rg[2] = vv;
        if (vv > rg[3]) rg[3] = vv;
        if (pp < rg[4]) rg[4] = pp;
        if (pp > rg[5]) rg[5] = pp;
    }
    if (out_ranges != NULL) memcpy(out_ranges, rg, sizeof(rg));
    sidecar_path(path, sizeof(path), out_dir, id, "uvphi.f32");
    if (write_blob(path, uvphi, nb) != 0) { Arena_restore(arena, mark); return -1; }

    /* skin */
    SkinVert *skin = NULL;
    size_t nskin = 0;
    if (PlacedCube_load_skin(arena, out_dir, id, 1, &skin, &nskin) != 0) {
        Arena_restore(arena, mark);
        return -1;
    }
    for (size_t i = 0; i < nskin; i++) {
        int32_t sgid = skin[i].gid;
        double sphi_raw = (double)skin[i].phi;
        skin[i] = CubeReg_apply_vert(skin[i], reg, spiral_a, spiral_b);
        skin[i].u = (float)((double)skin[i].u + pc_warp_eval(reg, sgid, sphi_raw));
    }
    sidecar_path(path, sizeof(path), out_dir, id, "skin.f32");
    if (write_blob(path, skin, nskin * sizeof(SkinVert)) != 0) {
        Arena_restore(arena, mark);
        return -1;
    }

    if (write_obj) {
        float *verts = NULL;
        int32_t *faces = NULL;
        size_t mnv = 0, mnf = 0;
        sidecar_path(path, sizeof(path), out_dir, id, "mesh.obj");
        if (ObjIO_read(arena, path, &verts, &mnv, &faces, &mnf) != 0 || mnv != nv) {
            Arena_restore(arena, mark);
            return -1;
        }
        size_t kb = 0;
        sidecar_path(path, sizeof(path), out_dir, id, "facekeep.u8");
        uint8_t *keep = (uint8_t *)read_blob(arena, path, &kb);
        if (keep == NULL || kb != mnf) { Arena_restore(arena, mark); return -1; }

        float *uv = (float *)ARENA_ALLOC(arena, (nv ? nv : 1) * 2 * sizeof(float));
        for (size_t i = 0; i < nv; i++) {
            uv[i * 2 + 0] = uvphi[i * 3 + 0];
            uv[i * 2 + 1] = uvphi[i * 3 + 1];
        }
        sidecar_path(path, sizeof(path), out_dir, id, "placed.obj");
        if (ObjIO_write_uv_masked(path, verts, nv, faces, mnf, uv, keep) != 0) {
            Arena_restore(arena, mark);
            return -1;
        }
    }

    Arena_restore(arena, mark);
    return 0;
}

int PlacedCube_load_skin(Arena_T arena, const char *out_dir, const char *id,
                         int raw, SkinVert **out, size_t *out_n)
{
    assert(arena && out_dir && id && out && out_n);
    *out = NULL;
    *out_n = 0;
    char path[1024];
    sidecar_path(path, sizeof(path), out_dir, id,
                 raw ? "skin_raw.f32" : "skin.f32");
    size_t nb = 0;
    SkinVert *skin = (SkinVert *)read_blob(arena, path, &nb);
    if (skin == NULL || nb % sizeof(SkinVert) != 0) return -1;
    *out = skin;
    *out_n = nb / sizeof(SkinVert);
    return 0;
}

/* ============================================================================
 * Self-test: synthetic spiral cube end-to-end through unwrap + finalize.
 * ==========================================================================*/

static int pc_check(int cond, const char *what, int *fails)
{
    if (!cond) {
        fprintf(stderr, "[placed_cube selftest]   FAIL: %s\n", what);
        (*fails)++;
    }
    return cond;
}

/* Archimedean spiral sheet around +z through (cy,cx) (st_build_spiral shape) */
static void pc_build_spiral(Arena_T arena, int nphi, int nh, double r0,
                            double g, double turns, double z0, double Hgt,
                            double cy, double cx,
                            float **out_v, size_t *out_nv,
                            int32_t **out_f, size_t *out_nf)
{
    size_t nvv = (size_t)nphi * (size_t)nh;
    float *v = (float *)ARENA_ALLOC(arena, nvv * 3 * sizeof(float));
    int32_t *f = (int32_t *)ARENA_ALLOC(arena,
                     (size_t)(nphi - 1) * (size_t)(nh - 1) * 6 * sizeof(int32_t));
    double amax = turns * 2.0 * M_PI;
    for (int j = 0; j < nh; j++) {
        double z = z0 + Hgt * (double)j / (double)(nh - 1);
        for (int i = 0; i < nphi; i++) {
            double a = amax * (double)i / (double)(nphi - 1);
            double rr = r0 + g * a;
            size_t idx = (size_t)j * (size_t)nphi + (size_t)i;
            v[idx * 3 + 0] = (float)z;
            v[idx * 3 + 1] = (float)(cy + rr * sin(a));
            v[idx * 3 + 2] = (float)(cx + rr * cos(a));
        }
    }
    size_t fi = 0;
    for (int j = 0; j < nh - 1; j++) {
        for (int i = 0; i < nphi - 1; i++) {
            int32_t a = (int32_t)((size_t)j * (size_t)nphi + (size_t)i);
            int32_t b = a + 1;
            int32_t c = (int32_t)((size_t)(j + 1) * (size_t)nphi + (size_t)i);
            int32_t d = c + 1;
            f[fi * 3 + 0] = a; f[fi * 3 + 1] = b; f[fi * 3 + 2] = c; fi++;
            f[fi * 3 + 0] = b; f[fi * 3 + 1] = d; f[fi * 3 + 2] = c; fi++;
        }
    }
    *out_v = v; *out_nv = nvv; *out_f = f; *out_nf = fi;
}

int PlacedCube_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    const char *dir = "output/_selftest_scroll_whole";
    const char *id = "z00000_y00000_x00000";
    char path[1024];
    snprintf(path, sizeof(path), "%s/x", dir);
    ves_ensure_parent_dir(path);

    /* Nearest-slice labels may disagree across a valid face.  Kept topology
     * must canonicalize them, while a dropped face must not bridge groups. */
    {
        int32_t tf[6] = { 0, 1, 2, 3, 4, 5 };
        int32_t tg[6] = { 0, 1, 2, 3, 4, 5 };
        uint8_t tk[2] = { 1, 0 };
        int ng = 0;
        size_t nm = pc_merge_kept_groups(arena, tf, 2, tk, tg, 6, &ng);
        pc_check(nm == 2 && ng == 4, "kept topology merges group domains",
                 &fails);
        pc_check(tg[0] == tg[1] && tg[1] == tg[2],
                 "kept face has one canonical group", &fails);
        pc_check(tg[3] != tg[4] && tg[4] != tg[5] && tg[3] != tg[5],
                 "dropped face does not merge groups", &fails);
    }

    /* spiral inside a 64-cube: r up to ~30 of center 32 -> skin comes from the
     * z faces (z in [2,62] -> 2 vox from both) and the outermost turn */
    const double r0 = 15.0, g = 1.2, turns = 2.0, cy = 32.0, cx = 32.0;
    const double zlo = 2.0, Hgt = 60.0;
    const int64_t origin[3] = { 0, 0, 0 };
    const int64_t chunk = 64;
    float *v = NULL;
    int32_t *f = NULL;
    size_t nvv = 0, nfc = 0;
    pc_build_spiral(arena, 500, 31, r0, g, turns, zlo, Hgt, cy, cx,
                    &v, &nvv, &f, &nfc);
    snprintf(path, sizeof(path), "%s/%s_input.obj", dir, id);
    pc_check(ObjIO_write(path, v, nvv, f, nfc) == 0, "write input obj", &fails);

    /* calibration mirrors the driver: fit the seed (here: the same cube)
     * UNPINNED first, then pin its (a,b) everywhere. Pinning the generator's
     * (r0, 2pi*g) directly is WRONG -- the radius anchor gauges phi so the
     * line's intercept lands near 0 (r ~= pitch*phi/2pi), not at r0. */
    RibbonOpts ro;
    RibbonOpts_default(&ro);
    ro.axis_point[1] = (float)cy;
    ro.axis_point[2] = (float)cx;
    ro.wrap_spacing = (float)(2.0 * M_PI * g);
    ro.emit_global = 1;
    ro.pin_orient = 1;
    ro.fit_ribbon = 0;
    PlacedStats cal;
    pc_check(PlacedCube_unwrap(arena, path, dir, id, &ro, 1, origin, chunk,
                               4.0, 4.0, 40.0, 0.0, &cal) == 0 && cal.status == 0,
             "calibration unwrap ok", &fails);
    pc_check(fabs(fabs(cal.spiral_b) - 2.0 * M_PI * g) < 0.5,
             "calibration |b| ~= 2pi*g", &fails);
    ro.winding_sense = cal.spiral_b > 0.0 ? 1 : -1;
    ro.spiral_a = cal.spiral_a;
    ro.spiral_b = cal.spiral_b;

    PlacedStats st;
    int rc = PlacedCube_unwrap(arena, path, dir, id, &ro, 1, origin, chunk,
                               4.0, 4.0, 40.0, 0.0, &st);
    fprintf(stderr, "[placed_cube selftest] unwrap rc=%d status=%d nv=%zu nf=%zu "
            "skin=%zu bad=%zu r2=%.3f u=[%.1f,%.1f] v=[%.1f,%.1f] (%.2fs)\n",
            rc, st.status, st.nv, st.nf, st.n_skin, st.n_badface, st.spiral_r2,
            st.u_min, st.u_max, st.v_min, st.v_max, st.seconds);
    pc_check(rc == 0 && st.status == 0, "unwrap ok", &fails);
    pc_check(st.nv >= nvv, "post-sever nv >= input nv", &fails);
    pc_check(st.n_skin > 0, "skin nonempty", &fails);
    pc_check(st.spiral_r2 > 0.9, "pinned line fits the spiral", &fails);
    pc_check(st.v_min > zlo - 1.0 && st.v_max < zlo + Hgt + 1.0,
             "v is absolute world z", &fails);

    /* skin sanity: every skin vert within 4 vox of the box boundary */
    SkinVert *skin = NULL;
    size_t nskin = 0;
    pc_check(PlacedCube_load_skin(arena, dir, id, 1, &skin, &nskin) == 0 &&
             nskin == st.n_skin, "skin loads back", &fails);
    int skin_ok = 1;
    for (size_t i = 0; i < nskin; i++) {
        double bd = 1e300;
        for (int ax = 0; ax < 3; ax++) {
            double lo = (double)skin[i].p[ax] - (double)origin[ax];
            double hi = (double)(origin[ax] + chunk) - (double)skin[i].p[ax];
            double d = lo < hi ? lo : hi;
            if (d < bd) bd = d;
        }
        if (bd > 4.0 + 1e-3) skin_ok = 0;
    }
    pc_check(skin_ok, "all skin verts near the box", &fails);

    /* group sidecar exists and every mapped skin vertex carries a valid gid */
    {
        char gp[1024];
        snprintf(gp, sizeof(gp), "%s/%s_group.i32", dir, id);
        FILE *gf = fopen(gp, "rb");
        long gsz = -1;
        if (gf) { fseek(gf, 0, SEEK_END); gsz = ftell(gf); fclose(gf); }
        pc_check(gsz == (long)(st.nv * sizeof(int32_t)),
                 "group sidecar sized [nv]", &fails);
        int gid_ok = 1;
        for (size_t i = 0; i < nskin; i++)
            if (skin[i].gid < 0) gid_ok = 0;
        pc_check(gid_ok, "skin gids valid", &fails);
    }

    /* finalize with a deliberate per-group (k, du) table; verify roundtrip */
    const int32_t k = 1;
    const double du = 5.0;
    int32_t t_wk[1] = { k };
    double t_du[1] = { du };
    PlacedReg reg;
    memset(&reg, 0, sizeof reg);   /* NULLs the optional g_warp_* fields */
    reg.g_wk = t_wk;
    reg.g_du = t_du;
    reg.n_groups = 1;
    reg.wk_cube = k;    /* same fallback: cube has one connected group */
    reg.du_cube = du;
    double rg[6] = { 0 };
    pc_check(PlacedCube_finalize(arena, dir, id, &reg, ro.spiral_a,
                                 ro.spiral_b, 1, rg) == 0, "finalize ok", &fails);
    pc_check(rg[1] > rg[0] && rg[3] > rg[2] && rg[5] > rg[4],
             "finalize ranges sane", &fails);
    SkinVert *sreg = NULL;
    size_t nreg = 0;
    pc_check(PlacedCube_load_skin(arena, dir, id, 0, &sreg, &nreg) == 0 &&
             nreg == nskin, "registered skin loads", &fails);
    double emax = 0.0;
    for (size_t i = 0; i < nskin && i < nreg; i++) {
        double want_phi = (double)skin[i].phi + 2.0 * M_PI * (double)k;
        double want_u = (double)skin[i].u
                      + CubeReg_deltaU(ro.spiral_a, ro.spiral_b,
                                       (double)skin[i].phi, k) + du;
        double e1 = fabs((double)sreg[i].phi - want_phi);
        double e2 = fabs((double)sreg[i].u - want_u);
        if (e1 > emax) emax = e1;
        if (e2 > emax) emax = e2;
    }
    pc_check(emax < 1e-3, "finalize applies (k,du) exactly", &fails);

    /* placed.obj exists and parses; kept faces <= nf */
    {
        float *pv = NULL;
        int32_t *pf = NULL;
        size_t pnv = 0, pnf = 0;
        snprintf(path, sizeof(path), "%s/%s_placed.obj", dir, id);
        pc_check(ObjIO_read(arena, path, &pv, &pnv, &pf, &pnf) == 0 &&
                 pnv == st.nv && pnf + st.n_badface == st.nf,
                 "placed.obj verts align, kept faces = nf - bad", &fails);
    }

    Arena_dispose(&arena);
    fprintf(stderr, "[placed_cube selftest] %s (%d failure%s)\n",
            fails == 0 ? "PASSED" : "FAILED", fails, fails == 1 ? "" : "s");
    return fails;
}
