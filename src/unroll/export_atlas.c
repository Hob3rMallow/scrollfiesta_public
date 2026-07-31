/* export_atlas.c -- see export_atlas.h. Piece enumeration + per-piece tifxyz
 * export + atlas.json manifest. (Piece planning is inline here for v1; if it
 * grows it factors out into piece_select.{c,h}.) */
#include "../common/ves_platform.h"

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "export_atlas.h"
#include "tifxyz_export.h"

void AtlasOpts_default(AtlasOpts *o)
{
    assert(o);
    memset(o, 0, sizeof(*o));
    o->piece_wraps = 1;
    o->slab_v = 4096.0;
    o->du = 1.0;
    o->dv = 1.0;
    o->write_winding = 1;
    o->max_edge3d = 0.0;
    o->stretch_ratio = 4.0;
    o->stretch_floor = 25.0;
    o->conflict_dist = 2.0;
    o->relax = 0;
    o->relax_sweeps = 30;
}

typedef struct {
    char   name[96];
    int    k0, k1;
    double v0, v1;
    size_t W, H, valid, multi, conflict;
    double origin_u, origin_v;
    double bbox_lo[3], bbox_hi[3];
} PieceRec;

static int cmp_float(const void *pa, const void *pb)
{
    float a = *(const float *)pa, b = *(const float *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

int ExportAtlas_run(Arena_T arena, const PieceSet *ps, const ScaffoldCalib *c,
                    const char *seg_root, const char *prefix,
                    const AtlasOpts *opts_in, AtlasStats *out)
{
    assert(arena && ps && c && seg_root && prefix && out);
    memset(out, 0, sizeof(*out));
    if (ps->nf == 0 || ps->nv == 0 || ps->phi == NULL) return -1;

    AtlasOpts opts;
    if (opts_in) opts = *opts_in; else AtlasOpts_default(&opts);
    if (opts.piece_wraps < 1) opts.piece_wraps = 1;
    if (opts.du <= 0.0) opts.du = 1.0;
    if (opts.dv <= 0.0) opts.dv = 1.0;
    if (opts.slab_v <= 0.0) opts.slab_v = 4096.0;

    double t_start = ves_clock_sec();

    Arena_Mark mark = Arena_save(arena);

    /* per-vertex unwound turn coordinate. */
    double *w = (double *)ARENA_ALLOC(arena, (size_t)(ps->nv * sizeof(double)));
    for (size_t v = 0; v < ps->nv; v++) w[v] = Scaffold_unwound(c, (double)ps->phi[v]);

    /* per-face: centroid k, centroid v (== world z), centroid u, wrap span. */
    int32_t *facek = (int32_t *)ARENA_ALLOC(arena, (size_t)(ps->nf * sizeof(int32_t)));
    float   *facev = (float *)ARENA_ALLOC(arena, (size_t)(ps->nf * sizeof(float)));
    float   *facecu = (float *)ARENA_ALLOC(arena, (size_t)(ps->nf * sizeof(float)));
    uint8_t *fquar = (uint8_t *)ARENA_CALLOC(arena, (size_t)ps->nf, 1L);

    int kmin = INT_MAX, kmax = INT_MIN;
    double vmin = 1e300, vmax = -1e300;
    size_t nquar = 0;
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t a = ps->faces[f*3+0], b = ps->faces[f*3+1], cc = ps->faces[f*3+2];
        double wa = w[a], wb = w[b], wc = w[cc];
        double lo = wa < wb ? (wa < wc ? wa : wc) : (wb < wc ? wb : wc);
        double hi = wa > wb ? (wa > wc ? wa : wc) : (wb > wc ? wb : wc);
        int klo = (int)floor(lo / SCAFFOLD_2PI), khi = (int)floor(hi / SCAFFOLD_2PI);
        if (khi - klo >= 2) { fquar[f] = 1; nquar++; continue; }  /* contested */
        double wc_mean = (wa + wb + wc) / 3.0;
        int kf = (int)floor(wc_mean / SCAFFOLD_2PI);
        double vf = ((double)ps->uv[a*2+1] + (double)ps->uv[b*2+1]
                     + (double)ps->uv[cc*2+1]) / 3.0;
        double uf = ((double)ps->uv[a*2+0] + (double)ps->uv[b*2+0]
                     + (double)ps->uv[cc*2+0]) / 3.0;
        facek[f] = kf; facev[f] = (float)vf; facecu[f] = (float)uf;
        if (kf < kmin) kmin = kf;
        if (kf > kmax) kmax = kf;
        if (vf < vmin) vmin = vf;
        if (vf > vmax) vmax = vf;
    }
    out->n_quarantine_faces = nquar;
    if (kmin > kmax) { Arena_restore(arena, mark); return -1; }  /* nothing kept */

    int nslab = (int)ceil((vmax - vmin + 1.0) / opts.slab_v);
    if (nslab < 1) nslab = 1;
    int n_kpieces = (kmax - kmin) / opts.piece_wraps + 1;
    size_t max_pieces = (size_t)n_kpieces * (size_t)nslab;

    PieceRec *recs = (PieceRec *)ARENA_ALLOC(arena,
                        (size_t)((max_pieces > 0 ? max_pieces : 1) * sizeof(PieceRec)));
    size_t nrec = 0;

    uint8_t *keep = (uint8_t *)ARENA_ALLOC(arena, (size_t)ps->nf);
    float   *us = (float *)ARENA_ALLOC(arena, (size_t)(ps->nf * sizeof(float)));
    size_t total_uoutliers = 0;

    for (int slab = 0; slab < nslab; slab++) {
        double v0 = vmin + (double)slab * opts.slab_v;
        double v1 = v0 + opts.slab_v;
        if (slab == nslab - 1) v1 = vmax + 1.0;

        for (int k0 = kmin; k0 <= kmax; k0 += opts.piece_wraps) {
            int k1 = k0 + opts.piece_wraps;

            /* pass 1: select by (k, v); collect centroid-u of kept faces. */
            size_t nkept = 0;
            memset(keep, 0, ps->nf);
            for (size_t f = 0; f < ps->nf; f++) {
                if (fquar[f]) continue;
                if (facek[f] < k0 || facek[f] >= k1) continue;
                if (facev[f] < v0 || facev[f] >= v1) continue;
                keep[f] = 1; us[nkept++] = facecu[f];
            }
            if (nkept == 0) { out->n_empty++; continue; }

            /* Robust u-window: a wrap's legit u-span is ~its circumference
             * (2*pi*r_k); mis-registered outlier faces (winding residual)
             * sit far outside it and would blow the canvas. Clip to
             * median_u +/- 1.6*circumference and drop faces beyond -- zero
             * drops for a clean wrap, only the outliers for a contaminated
             * one. This also keeps every piece under the TIFF cap. */
            qsort(us, nkept, sizeof(float), cmp_float);
            double med_u = (double)us[nkept / 2];
            double r_k = fabs(c->spiral_a + c->pitch * (double)k0);
            double u_half = SCAFFOLD_2PI * r_k * 1.6;
            if (u_half < 256.0) u_half = 256.0;   /* tiny inner wraps */

            /* pass 2: drop u-outliers, recompute the vertex window. */
            double wu0 = 1e300, wu1 = -1e300, wv0 = 1e300, wv1 = -1e300;
            size_t kept2 = 0, uout = 0;
            for (size_t f = 0; f < ps->nf; f++) {
                if (!keep[f]) continue;
                if (fabs((double)facecu[f] - med_u) > u_half) {
                    keep[f] = 0; uout++; continue;
                }
                kept2++;
                for (int e = 0; e < 3; e++) {
                    int32_t vi = ps->faces[f*3+e];
                    double uu = (double)ps->uv[vi*2+0], vv = (double)ps->uv[vi*2+1];
                    if (uu < wu0) wu0 = uu; if (uu > wu1) wu1 = uu;
                    if (vv < wv0) wv0 = vv; if (vv > wv1) wv1 = vv;
                }
            }
            total_uoutliers += uout;
            nkept = kept2;
            if (nkept == 0) { out->n_empty++; continue; }

            /* Pre-check the classic-TIFF band cap so an over-wide piece is
             * flagged cleanly rather than tripping the exporter's partial
             * write. (A single wrap can't be z-split to fit -- its width is
             * the circumference; the fix is a coarser du or angular sectors,
             * a follow-up. On the ownership-resolved chain the mis-registered
             * u-spread that causes this at raw placement mostly vanishes.) */
            {
                long Wp = (long)floor(wu1 / opts.du + 0.5)
                        - (long)floor(wu0 / opts.du + 0.5) + 1;
                long Hp = (long)floor(wv1 / opts.dv + 0.5)
                        - (long)floor(wv0 / opts.dv + 0.5) + 1;
                if (Wp < 1 || Hp < 1 ||
                    (double)Wp * (double)Hp > 950000000.0) {
                    out->n_over_cap++;
                    if (opts.verbose)
                        fprintf(stderr, "  [atlas] w%03d z%05ld over cap "
                                "(%ldx%ld px) -- skipped\n",
                                k0, (long)floor(v0 + 0.5), Wp, Hp);
                    continue;
                }
            }
            out->n_pieces++;

            PieceRec *pr = &recs[nrec];
            snprintf(pr->name, sizeof(pr->name), "%s_w%03d_z%05ld",
                     prefix, k0, (long)floor(v0 + 0.5));
            pr->k0 = k0; pr->k1 = k1; pr->v0 = v0; pr->v1 = v1;

            char seg_dir[1024];
            snprintf(seg_dir, sizeof(seg_dir), "%s/%s", seg_root, pr->name);

            TifxyzOpts to; TifxyzOpts_default(&to);
            to.du = opts.du; to.dv = opts.dv;
            to.max_edge3d = opts.max_edge3d;
            to.stretch_ratio = opts.stretch_ratio;
            to.stretch_floor = opts.stretch_floor;
            to.conflict_dist = opts.conflict_dist;
            to.write_winding = opts.write_winding;
            to.face_keep = keep;
            to.have_window = 1;
            to.win_u0 = wu0; to.win_u1 = wu1;
            to.win_v0 = wv0; to.win_v1 = wv1;

            TifxyzStats tst;
            int rc = TifxyzExport_run(arena, ps, seg_dir, pr->name, &to, &tst);
            if (rc != 0) {
                /* cap was pre-checked, so a failure here is zero-covered
                 * (every assigned face gated out as membrane/stretch) -- an
                 * effectively empty piece, not an error. */
                out->n_empty++;
                out->n_pieces--;
                if (opts.verbose)
                    fprintf(stderr, "  [atlas] piece %s empty after gating\n",
                            pr->name);
                continue;
            }
            out->n_written++;
            out->total_valid_px += tst.valid;
            out->total_multi_px += tst.multi;
            out->total_conflict_px += tst.conflicts;
            pr->W = tst.W; pr->H = tst.H;
            pr->valid = tst.valid; pr->multi = tst.multi; pr->conflict = tst.conflicts;
            pr->origin_u = tst.u0; pr->origin_v = tst.v0;
            memcpy(pr->bbox_lo, tst.bbox_lo, sizeof pr->bbox_lo);
            memcpy(pr->bbox_hi, tst.bbox_hi, sizeof pr->bbox_hi);
            nrec++;

            if (opts.verbose)
                fprintf(stderr, "  [atlas] %s: %zux%zu valid=%zu multi=%zu "
                        "conflict=%zu\n", pr->name, tst.W, tst.H,
                        tst.valid, tst.multi, tst.conflicts);
        }
    }

    out->n_quarantine_faces = nquar + total_uoutliers;
    out->conflict_frac = out->total_valid_px > 0
        ? (double)out->total_conflict_px / (double)out->total_valid_px : 0.0;

    /* atlas.json manifest. */
    {
        char path[1024];
        snprintf(path, sizeof(path), "%s/atlas.json", seg_root);
        if (ves_ensure_parent_dir(path) == 0) {
            FILE *mf = fopen(path, "w");
            if (mf != NULL) {
                fprintf(mf,
                    "{\n"
                    "  \"format\": \"scrollfiesta-atlas\", \"version\": 1,\n"
                    "  \"prefix\": \"%s\",\n"
                    "  \"du\": %.6f, \"dv\": %.6f,\n"
                    "  \"calibration\": { \"spiral_a\": %.6f, \"spiral_b\": %.6f, "
                    "\"pitch\": %.6f, \"sense\": %d,\n"
                    "    \"axis_point_zyx\": [%.3f, %.3f, %.3f] },\n"
                    "  \"piece_wraps\": %d, \"slab_v\": %.1f,\n"
                    "  \"k_range\": [%d, %d], \"n_pieces\": %zu, "
                    "\"n_quarantine_faces\": %zu,\n"
                    "  \"total_valid_px\": %zu, \"total_conflict_px\": %zu,\n"
                    "  \"pieces\": [\n",
                    prefix, opts.du, opts.dv,
                    c->spiral_a, c->spiral_b, c->pitch, c->sense,
                    (double)c->axis_point[0], (double)c->axis_point[1],
                    (double)c->axis_point[2],
                    opts.piece_wraps, opts.slab_v, kmin, kmax, nrec,
                    out->n_quarantine_faces,
                    out->total_valid_px, out->total_conflict_px);
                for (size_t i = 0; i < nrec; i++) {
                    PieceRec *pr = &recs[i];
                    fprintf(mf,
                        "    { \"uuid\": \"%s\", \"k\": [%d, %d], "
                        "\"v\": [%.1f, %.1f], \"W\": %zu, \"H\": %zu,\n"
                        "      \"valid_px\": %zu, \"multi_px\": %zu, "
                        "\"conflict_px\": %zu, \"origin_uv\": [%.3f, %.3f],\n"
                        "      \"bbox\": [[%.3f, %.3f, %.3f], [%.3f, %.3f, %.3f]] }%s\n",
                        pr->name, pr->k0, pr->k1, pr->v0, pr->v1, pr->W, pr->H,
                        pr->valid, pr->multi, pr->conflict, pr->origin_u, pr->origin_v,
                        pr->bbox_lo[0], pr->bbox_lo[1], pr->bbox_lo[2],
                        pr->bbox_hi[0], pr->bbox_hi[1], pr->bbox_hi[2],
                        i + 1 < nrec ? "," : "");
                }
                fprintf(mf, "  ]\n}\n");
                fclose(mf);
            }
        }
    }

    out->seconds = ves_clock_sec() - t_start;
    Arena_restore(arena, mark);
    return 0;
}

/* ------------------------------------------------------------------ selftest */

int ExportAtlas_selftest(void)
{
    int fails = 0;
    Arena_T ar = Arena_new();
    ScaffoldCalib c; Scaffold_calib_default(&c);

    /* Spiral ribbon, 4 turns, 2 z rows (same construction as scaffold's). */
    int NT = 240, WROWS = 2;
    double turns = 3.95;
    size_t nv = (size_t)NT * (size_t)WROWS;
    size_t nf = (size_t)(NT - 1) * (size_t)(WROWS - 1) * 2;
    PieceSet ps; memset(&ps, 0, sizeof ps);
    ps.verts = (float *)ARENA_ALLOC(ar, (size_t)(nv * 3 * sizeof(float)));
    ps.uv    = (float *)ARENA_ALLOC(ar, (size_t)(nv * 2 * sizeof(float)));
    ps.phi   = (float *)ARENA_ALLOC(ar, (size_t)(nv * sizeof(float)));
    ps.normals = (float *)ARENA_CALLOC(ar, (size_t)(nv * 3), (size_t)sizeof(float));
    ps.gid   = (int32_t *)ARENA_ALLOC(ar, (size_t)(nv * sizeof(int32_t)));
    ps.faces = (int32_t *)ARENA_ALLOC(ar, (size_t)(nf * 3 * sizeof(int32_t)));
    ps.face_cube = (int32_t *)ARENA_CALLOC(ar, (size_t)nf, (size_t)sizeof(int32_t));
    ps.ids = (char (*)[48])ARENA_ALLOC(ar, (size_t)(1 * 48));
    memcpy(ps.ids[0], "z00000_y00000_x00000", 21);
    ps.cube_voff = (size_t *)ARENA_ALLOC(ar, (size_t)(2 * sizeof(size_t)));
    ps.cube_voff[0] = 0; ps.cube_voff[1] = nv;
    ps.cube_org = (long (*)[3])ARENA_ALLOC(ar, (size_t)(3 * sizeof(long)));
    ps.cube_org[0][0] = 0; ps.cube_org[0][1] = 0; ps.cube_org[0][2] = 0;
    ps.n_cubes = 1; ps.nv = nv; ps.nf = nf;

    double umin = 1e9, umax = -1e9;
    for (int i = 0; i < NT; i++) {
        double t = (double)i / (double)(NT - 1) * turns;
        double theta = SCAFFOLD_2PI * t;
        double r = c.spiral_a + (-c.spiral_b) * t;
        double u = r * theta;   /* monotone arc-ish coordinate */
        if (u < umin) umin = u; if (u > umax) umax = u;
        for (int row = 0; row < WROWS; row++) {
            size_t v = (size_t)i * (size_t)WROWS + (size_t)row;
            ps.verts[v*3+0] = (float)row;
            ps.verts[v*3+1] = (float)(3405.0 + r * sin(theta));
            ps.verts[v*3+2] = (float)(2878.0 + r * cos(theta));
            ps.phi[v] = (float)(-theta);
            ps.uv[v*2+0] = (float)u; ps.uv[v*2+1] = (float)row;
            ps.gid[v] = (int32_t)floor(t);
        }
    }
    size_t fi = 0;
    for (int i = 0; i < NT - 1; i++)
        for (int row = 0; row < WROWS - 1; row++) {
            int32_t a = (int32_t)(i * WROWS + row);
            int32_t b = (int32_t)((i + 1) * WROWS + row);
            int32_t cc = (int32_t)(i * WROWS + row + 1);
            int32_t d = (int32_t)((i + 1) * WROWS + row + 1);
            ps.faces[fi*3+0]=a; ps.faces[fi*3+1]=b; ps.faces[fi*3+2]=cc; fi++;
            ps.faces[fi*3+0]=b; ps.faces[fi*3+1]=d; ps.faces[fi*3+2]=cc; fi++;
        }
    ps.u_min = umin; ps.u_max = umax; ps.v_min = 0; ps.v_max = WROWS;

    char root[1024];
    { char tmp[512]; ves_temp_dir(tmp, sizeof tmp);
      snprintf(root, sizeof root, "%s_selftest_export_atlas", tmp); }

    AtlasOpts o; AtlasOpts_default(&o);
    o.piece_wraps = 1; o.slab_v = 4096.0; o.write_winding = 1;
    AtlasStats st;
    if (ExportAtlas_run(ar, &ps, &c, root, "seg", &o, &st) != 0) {
        fprintf(stderr, "[export_atlas] run failed\n"); fails++;
    } else {
        /* k in {0,1,2,3}, 1 slab => 4 pieces written, no over-cap. */
        if (st.n_written != 4) { fprintf(stderr, "[export_atlas] wrote %zu want 4\n", st.n_written); fails++; }
        if (st.n_over_cap != 0) { fprintf(stderr, "[export_atlas] %zu over cap\n", st.n_over_cap); fails++; }
        /* a single-wrap piece must have essentially no inter-wrap conflict. */
        if (st.conflict_frac > 0.02) { fprintf(stderr, "[export_atlas] conflict_frac %.4f > 0.02\n", st.conflict_frac); fails++; }
        /* one segment's meta.json must exist. */
        char meta[1200];
        snprintf(meta, sizeof meta, "%s/seg_w000_z00000/meta.json", root);
        FILE *m = fopen(meta, "r");
        if (m == NULL) { fprintf(stderr, "[export_atlas] missing %s\n", meta); fails++; }
        else fclose(m);
    }

    Arena_dispose(&ar);
    if (fails == 0) fprintf(stderr, "[export_atlas] selftest PASSED\n");
    else fprintf(stderr, "[export_atlas] selftest FAILED (%d)\n", fails);
    return fails;
}
