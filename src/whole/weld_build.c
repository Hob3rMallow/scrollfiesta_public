/* weld_build.c -- Pass 1 tracked-weld correspondences. See weld_build.h. */
#include "weld_build.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/kdtree.h"

#define WB_PI 3.14159265358979323846

static double wb_wrap_pi(double x)
{
    return x - 2.0 * WB_PI * floor(x / (2.0 * WB_PI) + 0.5);
}

void WeldBuildOpts_default(WeldBuildOpts *o)
{
    assert(o);
    memset(o, 0, sizeof(*o));
    o->pair_gate = 3.5;
    o->radius_gate = 0.0;   /* => 0.5*pitch */
    o->frac_gate = 0.5;
    o->margin_frac = 0.5;
    o->verbose = 0;
}

int WeldBuild_pass1(Arena_T arena, const CubeNode *cnodes, size_t n_cubes,
                    SkinVert *const *skins, const size_t *nskin,
                    const float axis_point[3],
                    double spiral_a, double spiral_b, double pitch,
                    const WeldBuildOpts *opts, WeldTrack *out)
{
    assert(arena && cnodes && skins && nskin && axis_point && opts && out);
    memset(out, 0, sizeof(*out));
    if (n_cubes == 0)
        return -1;

    double rgate = opts->radius_gate > 0.0 ? opts->radius_gate : 0.5 * pitch;
    double gate2 = opts->pair_gate * opts->pair_gate;
    double margin2 = (1.0 + opts->margin_frac) * (1.0 + opts->margin_frac);

    out->version = WELD_TRACK_VERSION;
    out->spiral_a = spiral_a;
    out->spiral_b = spiral_b;
    out->pitch = pitch;
    out->axis_point[0] = axis_point[0];
    out->axis_point[1] = axis_point[1];
    out->axis_point[2] = axis_point[2];
    out->n_cubes = n_cubes;
    out->ids = (char (*)[WTRK_ID_LEN])ARENA_ALLOC(arena, n_cubes * WTRK_ID_LEN);
    for (size_t c = 0; c < n_cubes; c++)
        snprintf(out->ids[c], WTRK_ID_LEN, "%s", cnodes[c].id);

    size_t cap = 4096, nc = 0;
    WtrkCorr *cs = (WtrkCorr *)malloc(cap * sizeof(WtrkCorr));
    if (cs == NULL)
        return -1;
    size_t n_ambig = 0, per_axis[3] = { 0, 0, 0 };

    enum { NB_MAX = 16 };
    for (size_t c = 0; c < n_cubes; c++) {
        if (nskin[c] == 0)
            continue;
        for (int e6 = 1; e6 < 6; e6 += 2) {   /* +z,+y,+x: each seam once */
            int32_t j = cnodes[c].nbr[e6];
            if (j < 0 || nskin[(size_t)j] == 0)
                continue;
            uint8_t axis = e6 == 1 ? 0u : (e6 == 3 ? 1u : 2u);

            Arena_Mark mark = Arena_save(arena);
            float *jp = (float *)ARENA_ALLOC(arena, nskin[(size_t)j] * 3
                                             * sizeof(float));
            for (size_t q = 0; q < nskin[(size_t)j]; q++) {
                jp[q * 3 + 0] = skins[(size_t)j][q].p[0];
                jp[q * 3 + 1] = skins[(size_t)j][q].p[1];
                jp[q * 3 + 2] = skins[(size_t)j][q].p[2];
            }
            KDTree_T kd = KDTree_new(arena, jp, nskin[(size_t)j]);
            int32_t ball[NB_MAX];

            for (size_t q = 0; q < nskin[c]; q++) {
                const SkinVert *si = &skins[c][q];
                if (si->gid < 0)
                    continue;
                float d2 = 0.0f;
                size_t hit = KDTree_nearest(kd, si->p, &d2);
                if ((double)d2 > gate2)
                    continue;
                const SkinVert *sj = &skins[(size_t)j][hit];
                if (sj->gid < 0)
                    continue;
                double dyi = (double)si->p[1] - (double)axis_point[1];
                double dxi = (double)si->p[2] - (double)axis_point[2];
                double dyj = (double)sj->p[1] - (double)axis_point[1];
                double dxj = (double)sj->p[2] - (double)axis_point[2];
                double ri = sqrt(dyi * dyi + dxi * dxi);
                double rj = sqrt(dyj * dyj + dxj * dxj);
                double dr = fabs(ri - rj);
                if (dr > rgate)
                    continue;
                double dphi = wb_wrap_pi((double)sj->phi - (double)si->phi);
                if (fabs(dphi) > opts->frac_gate)
                    continue;

                /* ambiguity: another wrap of j within (1+margin)*d_nearest */
                uint8_t conf = 2;
                float br2 = (float)(margin2 * (double)d2);
                size_t nb = KDTree_ball_query(kd, si->p, br2, ball, NB_MAX);
                for (size_t bq = 0; bq < nb; bq++) {
                    if (ball[bq] == (int32_t)hit)
                        continue;
                    double dp = fabs((double)skins[(size_t)j][(size_t)ball[bq]]
                                         .phi - (double)sj->phi);
                    if (dp > WB_PI) { conf = 1; break; }
                }
                if (conf == 1)
                    n_ambig++;

                if (nc == cap) {
                    cap *= 2;
                    WtrkCorr *nn = (WtrkCorr *)realloc(cs,
                                       cap * sizeof(WtrkCorr));
                    if (nn == NULL) {
                        free(cs);
                        Arena_restore(arena, mark);
                        return -1;
                    }
                    cs = nn;
                }
                WtrkCorr *e = &cs[nc++];
                memset(e, 0, sizeof(*e));
                e->cube_a = (uint32_t)c;
                e->cube_b = (uint32_t)j;
                e->vert_a = (uint32_t)q;
                e->vert_b = (uint32_t)hit;
                e->gid_a = si->gid;
                e->gid_b = sj->gid;
                e->dr = (float)dr;
                e->dphi = (float)dphi;
                e->dist3d = (float)sqrt((double)d2);
                e->seam_axis = axis;
                e->conf = conf;
                per_axis[axis]++;
            }
            Arena_restore(arena, mark);
        }
    }

    out->n_corr = nc;
    out->corr = (WtrkCorr *)ARENA_ALLOC(arena, (nc + 1) * sizeof(WtrkCorr));
    memcpy(out->corr, cs, nc * sizeof(WtrkCorr));
    free(cs);
    out->self = out;

    if (opts->verbose)
        fprintf(stderr, "[weld_build] %zu correspondences (z=%zu y=%zu x=%zu) "
                "ambiguous(conf1)=%zu (%.1f%%)\n", nc, per_axis[0],
                per_axis[1], per_axis[2], n_ambig,
                nc ? 100.0 * (double)n_ambig / (double)nc : 0.0);
    return 0;
}

int WeldBuild_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    WeldBuildOpts o;
    WeldBuildOpts_default(&o);
    o.verbose = 0;

    enum { NS = 20 };
    SkinVert sa[NS + 1], sb[NS + 1];
    float ap[3] = { 0.0f, 0.0f, 0.0f };
    double pitch = 10.0, a = 1.0, b = -10.0;
    for (int i = 0; i < NS; i++) {
        sa[i].p[0] = (float)(2 * i);   /* z along the seam */
        sa[i].p[1] = 30.0f;            /* y */
        sa[i].p[2] = 40.0f;            /* x (radius ~ 50) */
        sa[i].u = 100.0f;
        sa[i].v = (float)(2 * i);
        sa[i].phi = -3.0f;
        sa[i].gid = 0;
        sb[i] = sa[i];
        sb[i].p[0] += 0.3f;            /* within the pair gate */
        sb[i].phi = (float)(-3.0 + 2.0 * WB_PI);   /* same wrap, +1 raw turn */
        sb[i].gid = 2;
    }
    /* a torn contaminant: 3D-close but a non-integer raw dphi -> frac reject */
    sa[NS] = sa[0];
    sa[NS].p[0] = 100.0f;
    sb[NS] = sa[NS];
    sb[NS].p[0] += 0.2f;
    sb[NS].phi = (float)(-3.0 + 1.0);   /* wrap_pi(1.0) = 1.0 > frac_gate */
    sb[NS].gid = 2;

    CubeNode cn[2];
    memset(cn, 0, sizeof(cn));
    snprintf(cn[0].id, sizeof(cn[0].id), "z00000_y00000_x00000");
    snprintf(cn[1].id, sizeof(cn[1].id), "z00000_y00000_x00128");
    for (int i = 0; i < 6; i++) { cn[0].nbr[i] = -1; cn[1].nbr[i] = -1; }
    cn[0].nbr[5] = 1;   /* +x neighbor */
    SkinVert *skins[2] = { sa, sb };
    size_t nskin[2] = { (size_t)NS + 1, (size_t)NS + 1 };

    WeldTrack wt;
    if (WeldBuild_pass1(arena, cn, 2, skins, nskin, ap, a, b, pitch, &o, &wt)
        != 0) {
        fprintf(stderr, "[weld_build selftest]   FAIL: pass1 rc\n");
        fails++;
    } else {
        int all_x = 1, all_gid = 1, all_clean = 1, all_conf2 = 1;
        for (size_t k = 0; k < wt.n_corr; k++) {
            if (wt.corr[k].seam_axis != 2) all_x = 0;
            if (wt.corr[k].gid_a != 0 || wt.corr[k].gid_b != 2) all_gid = 0;
            if (fabs((double)wt.corr[k].dphi) > 0.01) all_clean = 0;
            if (wt.corr[k].conf != 2) all_conf2 = 0;
        }
        if (wt.n_corr != (size_t)NS) {
            fprintf(stderr, "[weld_build selftest]   FAIL: n_corr=%zu (want %d,"
                    " contaminant should be frac-rejected)\n", wt.n_corr, NS);
            fails++;
        }
        if (!all_x) { fprintf(stderr, "[weld_build selftest]   FAIL: seam_axis"
                              " not all x\n"); fails++; }
        if (!all_gid) { fprintf(stderr, "[weld_build selftest]   FAIL: gid\n");
                        fails++; }
        if (!all_clean) { fprintf(stderr, "[weld_build selftest]   FAIL: dphi "
                          "not near 0\n"); fails++; }
        if (!all_conf2) { fprintf(stderr, "[weld_build selftest]   FAIL: clean "
                          "pairs should be conf=2\n"); fails++; }
        if (fabs(wt.spiral_b - b) > 1e-9 || wt.n_cubes != 2) {
            fprintf(stderr, "[weld_build selftest]   FAIL: header\n"); fails++;
        }
    }
    Arena_dispose(&arena);
    fprintf(stderr, "[weld_build selftest] %s (%d failure%s)\n",
            fails == 0 ? "PASSED" : "FAILED", fails, fails == 1 ? "" : "s");
    return fails;
}
