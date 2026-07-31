/* cube_register.c -- per-group cross-seam winding + residual-u registration.
 * See cube_register.h. */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/kdtree.h"
#include "cube_register.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int cmp_dbl_cr(const void *pa, const void *pb)
{
    double a = *(const double *)pa, b = *(const double *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

/* median of buf[0..n) -- SORTS buf in place */
static double med_inplace(double *buf, size_t n)
{
    assert(n > 0);
    qsort(buf, n, sizeof(double), cmp_dbl_cr);
    return buf[n / 2];
}

int CubeReg_solve(Arena_T scratch, Arena_T tab_arena,
                  const SkinVert *skin_n, size_t n_n,
                  const SkinVert *pooled_p, size_t n_p,
                  double spiral_a, double spiral_b,
                  double pair_gate, size_t min_pairs, size_t min_group_pairs,
                  CubeReg *out)
{
    assert(scratch && tab_arena && out);
    memset(out, 0, sizeof(*out));
    if (min_pairs == 0) min_pairs = 1;
    if (min_group_pairs == 0) min_group_pairs = 1;
    if (n_n == 0 || n_p == 0 || skin_n == NULL || pooled_p == NULL) {
        out->low_conf = 1;
        return 0;
    }

    /* the per-group table must be allocated BEFORE the scratch mark:
     * scratch and tab_arena may be the SAME arena, and the restore below
     * would otherwise free the table out from under the caller */
    int32_t n_groups = 0;
    for (size_t i = 0; i < n_n; i++)
        if (skin_n[i].gid + 1 > n_groups) n_groups = skin_n[i].gid + 1;
    int32_t *g_wk = NULL;
    double *g_du = NULL;
    if (n_groups > 0) {
        g_wk = (int32_t *)ARENA_ALLOC(tab_arena,
                                      (size_t)n_groups * sizeof(int32_t));
        g_du = (double *)ARENA_ALLOC(tab_arena,
                                     (size_t)n_groups * sizeof(double));
    }

    Arena_Mark mark = Arena_save(scratch);

    /* KD-tree over the pooled placed-neighbor skin positions */
    float *pts = (float *)ARENA_ALLOC(scratch, n_p * 3 * sizeof(float));
    for (size_t i = 0; i < n_p; i++) {
        pts[i * 3 + 0] = pooled_p[i].p[0];
        pts[i * 3 + 1] = pooled_p[i].p[1];
        pts[i * 3 + 2] = pooled_p[i].p[2];
    }
    KDTree_T tree = KDTree_new(scratch, pts, n_p);
    float gate_sq = (float)(pair_gate * pair_gate);

    double *dphi = (double *)ARENA_ALLOC(scratch, (n_n + 1) * sizeof(double));
    int32_t *pidx = (int32_t *)ARENA_ALLOC(scratch, (n_n + 1) * sizeof(int32_t));
    int32_t *nidx = (int32_t *)ARENA_ALLOC(scratch, (n_n + 1) * sizeof(int32_t));
    size_t np = 0;
    for (size_t i = 0; i < n_n; i++) {
        float d2 = 0.0f;
        size_t j = KDTree_nearest(tree, skin_n[i].p, &d2);
        if (d2 > gate_sq) continue;
        dphi[np] = (double)pooled_p[j].phi - (double)skin_n[i].phi;
        pidx[np] = (int32_t)j;
        nidx[np] = (int32_t)i;
        np++;
    }
    out->n_pairs = np;
    if (np < min_pairs) {
        out->low_conf = 1;
        Arena_restore(scratch, mark);
        return 0;
    }

    /* cube-level fallback correction from ALL pairs */
    double *tmp = (double *)ARENA_ALLOC(scratch, (np + 1) * sizeof(double));
    memcpy(tmp, dphi, np * sizeof(double));
    int32_t k_cube = (int32_t)lround(med_inplace(tmp, np) / (2.0 * M_PI));
    for (size_t i = 0; i < np; i++) {
        const SkinVert *sn = &skin_n[nidx[i]];
        const SkinVert *sp = &pooled_p[pidx[i]];
        tmp[i] = (double)sp->u
               - ((double)sn->u + CubeReg_deltaU(spiral_a, spiral_b,
                                                 (double)sn->phi, k_cube));
    }
    double du_cube = med_inplace(tmp, np);
    out->tab.wk_cube = k_cube;
    out->tab.du_cube = du_cube;

    /* per-group table: bucket the pairs by the N vertex's group id */
    if (n_groups > 0) {
        /* counting sort pairs by gid */
        int32_t *goff = (int32_t *)ARENA_CALLOC(scratch,
                                                (size_t)n_groups + 2,
                                                sizeof(int32_t));
        for (size_t i = 0; i < np; i++) {
            int32_t g = skin_n[nidx[i]].gid;
            if (g >= 0) goff[g + 1]++;
        }
        for (int32_t g = 0; g < n_groups; g++) goff[g + 1] += goff[g];
        int32_t *gpair = (int32_t *)ARENA_ALLOC(scratch,
                                                (np + 1) * sizeof(int32_t));
        int32_t *gcur = (int32_t *)ARENA_ALLOC(scratch,
                                               ((size_t)n_groups + 1) *
                                               sizeof(int32_t));
        memcpy(gcur, goff, (size_t)n_groups * sizeof(int32_t));
        for (size_t i = 0; i < np; i++) {
            int32_t g = skin_n[nidx[i]].gid;
            if (g >= 0) gpair[gcur[g]++] = (int32_t)i;
        }
        for (int32_t g = 0; g < n_groups; g++) {
            size_t g0 = (size_t)goff[g], g1 = (size_t)goff[g + 1];
            size_t ng = g1 - g0;
            if (ng < min_group_pairs) {
                g_wk[g] = k_cube;
                g_du[g] = du_cube;
                continue;
            }
            for (size_t q = g0; q < g1; q++) tmp[q - g0] = dphi[gpair[q]];
            int32_t kg = (int32_t)lround(med_inplace(tmp, ng) / (2.0 * M_PI));
            for (size_t q = g0; q < g1; q++) {
                const SkinVert *sn = &skin_n[nidx[gpair[q]]];
                const SkinVert *sp = &pooled_p[pidx[gpair[q]]];
                tmp[q - g0] = (double)sp->u
                            - ((double)sn->u
                               + CubeReg_deltaU(spiral_a, spiral_b,
                                                (double)sn->phi, kg));
            }
            g_wk[g] = kg;
            g_du[g] = med_inplace(tmp, ng);
            out->n_groups_direct++;
        }
        out->tab.g_wk = g_wk;
        out->tab.g_du = g_du;
        out->tab.n_groups = n_groups;
    }
    out->tab.g_warp_nk = NULL;   /* no ArcReg warp from the per-cube solver */

    /* summary MADs of the residuals AFTER the per-group correction */
    for (size_t i = 0; i < np; i++) {
        int32_t k = 0;
        double du = 0.0;
        CubeReg_pick(&out->tab, skin_n[nidx[i]].gid, &k, &du);
        tmp[i] = fabs(dphi[i] - 2.0 * M_PI * (double)k);
    }
    out->dphi_mad = med_inplace(tmp, np);
    for (size_t i = 0; i < np; i++) {
        const SkinVert *sn = &skin_n[nidx[i]];
        const SkinVert *sp = &pooled_p[pidx[i]];
        int32_t k = 0;
        double du = 0.0;
        CubeReg_pick(&out->tab, sn->gid, &k, &du);
        tmp[i] = fabs((double)sp->u
                      - ((double)sn->u + CubeReg_deltaU(spiral_a, spiral_b,
                                                        (double)sn->phi, k))
                      - du);
    }
    out->du_mad = med_inplace(tmp, np);

    Arena_restore(scratch, mark);
    return 0;
}

/* ============================================================================
 * Self-test.
 * ==========================================================================*/

static int cr_check(int cond, const char *what, int *fails)
{
    if (!cond) {
        fprintf(stderr, "[cube_register selftest]   FAIL: %s\n", what);
        (*fails)++;
    }
    return cond;
}

/* deterministic tiny LCG so the test needs no libc rand state */
static double cr_rand01(uint32_t *st)
{
    *st = *st * 1664525u + 1013904223u;
    return (double)(*st >> 8) / 16777216.0;
}

int CubeReg_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    const double a = 3.0, b = -9.5;   /* PHerc-like negative-wound pinned line */

    /* (1) deltaU algebra: U(phi + 2pi k) - U(phi) for the quadratic map */
    {
        double phi = -37.2;
        int32_t k = 3;
        double U0 = a * phi + b * phi * phi / (4.0 * M_PI);
        double p1 = phi + 2.0 * M_PI * (double)k;
        double U1 = a * p1 + b * p1 * p1 / (4.0 * M_PI);
        double d = CubeReg_deltaU(a, b, phi, k);
        cr_check(fabs((U1 - U0) - d) < 1e-9, "deltaU matches the map", &fails);
        cr_check(CubeReg_deltaU(a, b, phi, 0) == 0.0, "deltaU k=0 is 0", &fails);
    }

    /* (2) SPLIT-WINDING cube: two groups whose true corrections DIFFER
     * (k=2/du=-4.25 vs k=3/du=+2.5) -- the exact failure a cube-level
     * correction cannot fix (measured seam dphi MAD ~= 2pi on real data). */
    {
        enum { NP = 400 };
        SkinVert P[NP], N[NP];
        uint32_t st = 12345u;
        const int32_t k_true[2] = { 2, 3 };
        const double du_true[2] = { -4.25, 2.5 };
        for (int i = 0; i < NP; i++) {
            double phi = -60.0 + 0.05 * (double)i;
            double r = a + b * phi / (2.0 * M_PI);
            double u = a * phi + b * phi * phi / (4.0 * M_PI);
            /* >= 4 vox apart so a <=1-vox jitter can never swap nearest
             * partners (the u step between adjacent samples is ~r*dphi
             * ~ 4.7 vox -- a mispair would poison the du median) */
            P[i].p[0] = (float)(10.0 + 4.0 * (double)i);
            P[i].p[1] = (float)(100.0 + r);
            P[i].p[2] = (float)(100.0 + 0.03 * (double)i);
            P[i].u = (float)u;
            P[i].v = P[i].p[0];
            P[i].phi = (float)phi;
            P[i].gid = 0;
            int g = i < NP / 2 ? 0 : 1;
            N[i] = P[i];
            N[i].gid = g;
            N[i].p[0] += (float)(cr_rand01(&st) * 2.0 - 1.0);
            N[i].p[1] += (float)(cr_rand01(&st) * 2.0 - 1.0);
            N[i].p[2] += (float)(cr_rand01(&st) * 2.0 - 1.0);
            /* N's raw phi is k_true turns BELOW P; its u sits on the pinned
             * map at that lower phi, minus a residual drift du_true */
            double phin = phi - 2.0 * M_PI * (double)k_true[g];
            N[i].phi = (float)phin;
            N[i].u = (float)(u - CubeReg_deltaU(a, b, phin, k_true[g])
                             - du_true[g]);
        }
        CubeReg R;
        CubeReg_solve(arena, arena, N, NP, P, NP, a, b, 3.5, 24, 8, &R);
        fprintf(stderr, "[cube_register selftest] (2) split-winding -> "
                "groups=%d direct=%d wk={%d,%d} du={%.3f,%.3f} pairs=%zu "
                "mad(phi)=%.4f mad(u)=%.3f\n",
                R.tab.n_groups, R.n_groups_direct,
                R.tab.n_groups > 1 ? R.tab.g_wk[0] : -99,
                R.tab.n_groups > 1 ? R.tab.g_wk[1] : -99,
                R.tab.n_groups > 1 ? R.tab.g_du[0] : 0.0,
                R.tab.n_groups > 1 ? R.tab.g_du[1] : 0.0,
                R.n_pairs, R.dphi_mad, R.du_mad);
        cr_check(R.n_pairs > NP / 2, "split: enough pairs", &fails);
        cr_check(!R.low_conf, "split: confident", &fails);
        cr_check(R.tab.n_groups == 2 && R.n_groups_direct == 2,
                 "split: both groups solved directly", &fails);
        if (R.tab.n_groups == 2) {
            cr_check(R.tab.g_wk[0] == k_true[0] && R.tab.g_wk[1] == k_true[1],
                     "split: per-group k recovered", &fails);
            cr_check(fabs(R.tab.g_du[0] - du_true[0]) < 0.5 &&
                     fabs(R.tab.g_du[1] - du_true[1]) < 0.5,
                     "split: per-group du recovered", &fails);
        }
        cr_check(R.dphi_mad < 0.1, "split: post-correction phi residual ~0",
                 &fails);
        cr_check(R.du_mad < 1.0, "split: post-correction u residual tight",
                 &fails);

        /* (2b) identity: registering P against itself gives k=0, du~0 */
        CubeReg I;
        CubeReg_solve(arena, arena, P, NP, P, NP, a, b, 3.5, 24, 8, &I);
        cr_check(I.tab.wk_cube == 0 && fabs(I.tab.du_cube) < 1e-6 &&
                 I.tab.n_groups == 1 && I.tab.g_wk[0] == 0,
                 "identity is identity", &fails);

        /* (2c) gate: shove P 10 vox away -> no pairs -> low_conf identity */
        SkinVert Pf[NP];
        memcpy(Pf, P, sizeof(P));
        for (int i = 0; i < NP; i++) Pf[i].p[1] += 10.0f;
        CubeReg G;
        CubeReg_solve(arena, arena, N, NP, Pf, NP, a, b, 3.5, 24, 8, &G);
        cr_check(G.low_conf && G.n_pairs == 0 && G.tab.n_groups == 0 &&
                 G.tab.wk_cube == 0, "gate: far skins give low_conf identity",
                 &fails);

        /* (2d) min_pairs: same data, absurd threshold -> low_conf */
        CubeReg M;
        CubeReg_solve(arena, arena, N, NP, P, NP, a, b, 3.5, NP + 1, 8, &M);
        cr_check(M.low_conf && M.tab.n_groups == 0, "min_pairs gates confidence",
                 &fails);

        /* (2e) apply_vert: N corrected by the table lands on P (u and phi) */
        if (R.tab.n_groups == 2) {
            double emax = 0.0;
            for (int i = 0; i < NP; i++) {
                SkinVert c = CubeReg_apply_vert(N[i], &R.tab, a, b);
                double e1 = fabs((double)c.phi - (double)P[i].phi);
                double e2 = fabs((double)c.u - (double)P[i].u);
                if (e1 > emax) emax = e1;
                if (e2 > emax) emax = e2;
            }
            cr_check(emax < 0.51, "apply_vert lands N on P", &fails);
        }
    }

    /* (3) empty inputs never crash */
    {
        CubeReg R;
        CubeReg_solve(arena, arena, NULL, 0, NULL, 0, a, b, 3.5, 24, 8, &R);
        cr_check(R.low_conf == 1 && R.tab.wk_cube == 0 && R.n_pairs == 0,
                 "empty input -> low_conf identity", &fails);
    }

    Arena_dispose(&arena);
    fprintf(stderr, "[cube_register selftest] %s (%d failure%s)\n",
            fails == 0 ? "PASSED" : "FAILED", fails, fails == 1 ? "" : "s");
    return fails;
}
