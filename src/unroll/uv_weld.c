/* uv_weld.c -- cross-cube boundary weld. See uv_weld.h. */
#include "uv_weld.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/kdtree.h"
#include "../common/ves_platform.h"

static void *uw_xmalloc(size_t n)
{
    void *p = malloc(n > 0 ? n : 1);
    if (p == NULL) { fprintf(stderr, "uv_weld: OOM %zu\n", n); exit(1); }
    return p;
}

static void *uw_xcalloc(size_t c, size_t s)
{
    void *p = calloc(c > 0 ? c : 1, s);
    if (p == NULL) { fprintf(stderr, "uv_weld: OOM %zux%zu\n", c, s); exit(1); }
    return p;
}

void UVWeldOpts_default(UVWeldOpts *o)
{
    assert(o);
    memset(o, 0, sizeof(*o));
    o->r3d = 3.5;
    o->skin = 4.0;
    o->du_gate = 6.0;
    o->dv_gate = 3.0;
    o->chunk = 128;
    o->verbose = 0;
}

static int uw_cmp_f64(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x < y) ? -1 : (x > y);
}

int UVWeld_run(Arena_T arena, PieceSet *ps, const UVWeldOpts *opts,
               UVWeldStats *out)
{
    assert(arena && ps && opts && out);
    memset(out, 0, sizeof(*out));
    if (ps->nv == 0 || ps->nf == 0 || ps->n_cubes < 2)
        return 0;
    const UVWeldOpts *o = opts;
    double r2 = o->r3d * o->r3d;

    /* verts referenced by kept faces */
    uint8_t *has = (uint8_t *)uw_xcalloc(ps->nv, 1);
    for (size_t f = 0; f < ps->nf; f++)
        for (int k = 0; k < 3; k++)
            has[(size_t)ps->faces[f * 3 + k]] = 1;

    /* boundary shell: within `skin` of the cube's box faces */
    int32_t *bidx = (int32_t *)uw_xmalloc(ps->nv * sizeof(int32_t));
    size_t nb = 0;
    for (size_t c = 0; c < ps->n_cubes; c++) {
        long oz = ps->cube_org[c][0], oy = ps->cube_org[c][1],
             ox = ps->cube_org[c][2];
        if (oz < 0) continue;    /* non-canonical id: no box known */
        double c0[3] = { (double)oz, (double)oy, (double)ox };
        for (size_t i = ps->cube_voff[c]; i < ps->cube_voff[c + 1]; i++) {
            if (!has[i]) continue;
            if (ps->gid != NULL && ps->gid[i] < 0) continue;
            double dmin = 1e300;
            for (int m = 0; m < 3; m++) {
                double t = (double)ps->verts[i * 3 + m] - c0[m];
                double d = t < (double)o->chunk - t ? t
                                                    : (double)o->chunk - t;
                if (d < dmin) dmin = d;
            }
            if (dmin < o->skin)
                bidx[nb++] = (int32_t)i;
        }
    }
    out->n_boundary = nb;
    if (nb < 2) {
        free(has);
        free(bidx);
        return 0;
    }

    /* KD-tree over the shell */
    float *bp = (float *)uw_xmalloc(nb * 3 * sizeof(float));
    for (size_t i = 0; i < nb; i++)
        for (int m = 0; m < 3; m++)
            bp[i * 3 + m] = ps->verts[(size_t)bidx[i] * 3 + m];
    KDTree_T kd = KDTree_new(arena, bp, nb);

    /* MUTUAL-NEAREST cross-cube matching. A plain "everything within r3d"
     * rule welds LATERAL neighbours along the sheet (vertex spacing 1-2 vox
     * < r3d) and chains whole seams into one blob; the duplicate we want is
     * the counterpart at (nearly) the SAME 3D spot, i.e. each vert's nearest
     * cross-cube candidate, accepted only when the choice is mutual. */
    int32_t *best = (int32_t *)uw_xmalloc(ps->nv * sizeof(int32_t));
    float *bestd = (float *)uw_xmalloc(ps->nv * sizeof(float));
    for (size_t i = 0; i < nb; i++) {
        best[bidx[i]] = -1;
        bestd[bidx[i]] = 1e30f;
    }
    int32_t hit[64];
    for (size_t i = 0; i < nb; i++) {
        int32_t gi = bidx[i];
        float q[3] = { bp[i * 3 + 0], bp[i * 3 + 1], bp[i * 3 + 2] };
        size_t nh = KDTree_ball_query(kd, q, (float)r2, hit, 64);
        size_t ci = 0;
        {
            size_t lo = 0, hi2 = ps->n_cubes;
            while (lo + 1 < hi2) {
                size_t m = (lo + hi2) / 2;
                if (ps->cube_voff[m] <= (size_t)gi) lo = m;
                else hi2 = m;
            }
            ci = lo;
        }
        for (size_t h = 0; h < nh; h++) {
            int32_t gj = bidx[hit[h]];
            if (gj == gi) continue;
            size_t cj = 0;
            {
                size_t lo = 0, hi2 = ps->n_cubes;
                while (lo + 1 < hi2) {
                    size_t m = (lo + hi2) / 2;
                    if (ps->cube_voff[m] <= (size_t)gj) lo = m;
                    else hi2 = m;
                }
                cj = lo;
            }
            if (ci == cj) { out->rej_cube++; continue; }
            double du = fabs((double)ps->uv[(size_t)gi * 2 + 0]
                             - (double)ps->uv[(size_t)gj * 2 + 0]);
            double dv = fabs((double)ps->uv[(size_t)gi * 2 + 1]
                             - (double)ps->uv[(size_t)gj * 2 + 1]);
            if (du > o->du_gate || dv > o->dv_gate) { out->rej_uv++; continue; }
            double d2 = 0.0;
            for (int m = 0; m < 3; m++) {
                double dd = (double)ps->verts[(size_t)gi * 3 + m]
                          - (double)ps->verts[(size_t)gj * 3 + m];
                d2 += dd * dd;
            }
            if ((float)d2 < bestd[gi]) {
                bestd[gi] = (float)d2;
                best[gi] = gj;
            }
        }
    }
    int32_t *ufp = (int32_t *)uw_xmalloc(ps->nv * sizeof(int32_t));
    for (size_t i = 0; i < ps->nv; i++) ufp[i] = (int32_t)i;
    double *dus = NULL;
    size_t ndus = 0, cdus = 0;
    for (size_t i = 0; i < nb; i++) {
        int32_t gi = bidx[i];
        int32_t gj = best[gi];
        if (gj < 0 || gj < gi || best[gj] != gi)
            continue;   /* not mutual (or counted from the other side) */
        int32_t ra = gi, rb = gj;
        while (ufp[ra] != ra) { ufp[ra] = ufp[ufp[ra]]; ra = ufp[ra]; }
        while (ufp[rb] != rb) { ufp[rb] = ufp[ufp[rb]]; rb = ufp[rb]; }
        if (ra != rb) ufp[rb] = ra;
        out->n_pairs++;
        double du = fabs((double)ps->uv[(size_t)gi * 2 + 0]
                         - (double)ps->uv[(size_t)gj * 2 + 0]);
        if (ndus == cdus) {
            cdus = cdus ? cdus * 2 : 1024;
            dus = (double *)realloc(dus, cdus * sizeof(double));
            if (dus == NULL) { fprintf(stderr, "uv_weld: OOM\n"); exit(1); }
        }
        dus[ndus++] = du;
    }
    free(best);
    free(bestd);
    if (ndus > 0) {
        qsort(dus, ndus, sizeof(double), uw_cmp_f64);
        out->du_absmed = dus[ndus / 2];
    }
    free(dus);
    free(bp);

    if (out->n_pairs == 0) {
        free(has);
        free(bidx);
        free(ufp);
        return 0;
    }

    /* cluster means -> representative slot; then retarget faces */
    /* pass 1: accumulate into root (pos, u, v, phi) */
    double *acc = (double *)uw_xcalloc(ps->nv, 6 * sizeof(double));
    int32_t *cnt = (int32_t *)uw_xcalloc(ps->nv, sizeof(int32_t));
    for (size_t i = 0; i < nb; i++) {
        int32_t v = bidx[i];
        int32_t r = v;
        while (ufp[r] != r) { ufp[r] = ufp[ufp[r]]; r = ufp[r]; }
        acc[(size_t)r * 6 + 0] += ps->verts[(size_t)v * 3 + 0];
        acc[(size_t)r * 6 + 1] += ps->verts[(size_t)v * 3 + 1];
        acc[(size_t)r * 6 + 2] += ps->verts[(size_t)v * 3 + 2];
        acc[(size_t)r * 6 + 3] += ps->uv[(size_t)v * 2 + 0];
        acc[(size_t)r * 6 + 4] += ps->uv[(size_t)v * 2 + 1];
        if (ps->phi != NULL)
            acc[(size_t)r * 6 + 5] += ps->phi[(size_t)v];
        cnt[r]++;
    }
    for (size_t i = 0; i < nb; i++) {
        int32_t v = bidx[i];
        if (ufp[v] == v && cnt[v] >= 2) {
            double n = (double)cnt[v];
            ps->verts[(size_t)v * 3 + 0] = (float)(acc[(size_t)v * 6 + 0] / n);
            ps->verts[(size_t)v * 3 + 1] = (float)(acc[(size_t)v * 6 + 1] / n);
            ps->verts[(size_t)v * 3 + 2] = (float)(acc[(size_t)v * 6 + 2] / n);
            ps->uv[(size_t)v * 2 + 0] = (float)(acc[(size_t)v * 6 + 3] / n);
            ps->uv[(size_t)v * 2 + 1] = (float)(acc[(size_t)v * 6 + 4] / n);
            if (ps->phi != NULL)
                ps->phi[(size_t)v] = (float)(acc[(size_t)v * 6 + 5] / n);
            out->n_clusters++;
        }
    }
    free(acc);
    free(cnt);

    /* retarget faces; count degenerates */
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t vr[3];
        for (int k = 0; k < 3; k++) {
            int32_t r = ps->faces[f * 3 + k];
            while (ufp[r] != r) { ufp[r] = ufp[ufp[r]]; r = ufp[r]; }
            if (r != ps->faces[f * 3 + k]) out->n_merged_verts++;
            ps->faces[f * 3 + k] = r;
            vr[k] = r;
        }
        if (vr[0] == vr[1] || vr[1] == vr[2] || vr[0] == vr[2])
            out->n_degen_faces++;
    }

    if (o->verbose)
        fprintf(stderr, "[uv_weld] shell=%zu pairs=%zu (rej cube/gid/uv="
                "%zu/%zu/%zu) clusters=%zu merged-refs=%zu degen=%zu "
                "|du|med=%.2f\n", out->n_boundary, out->n_pairs,
                out->rej_cube, out->rej_gid, out->rej_uv, out->n_clusters,
                out->n_merged_verts, out->n_degen_faces, out->du_absmed);

    free(has);
    free(bidx);
    free(ufp);
    return 0;
}

/* ============================================================================
 * Self-test.
 * ==========================================================================*/

static int uw_check(int cond, const char *what, int *fails)
{
    if (!cond) {
        fprintf(stderr, "[uv_weld selftest]   FAIL: %s\n", what);
        (*fails)++;
    }
    return cond;
}

/* two abutting 4x4 patches (cubes 0/1) whose shared boundary columns sit at
 * the same 3D spot with a small uv offset; plus a mock other-wrap patch 8 vox
 * away in 3D */
int UVWeld_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();

    enum { N = 4, SIDE = N + 1 };
    size_t slab_v = SIDE * SIDE;
    float verts[3 * SIDE * SIDE * 3];
    float uv[3 * SIDE * SIDE * 2];
    int32_t gid[3 * SIDE * SIDE];
    int32_t faces[3 * N * N * 2 * 3];
    int32_t face_cube[3 * N * N * 2];
    size_t nv = 0, nf = 0;
    size_t cube_voff[4] = { 0, 0, 0, 0 };
    long org[3][3] = { { 0, 0, 0 }, { 0, 0, 128 }, { 0, 0, 0 } };

    /* cube A: x in [124,128) (right edge of box at x=128) */
    for (int iz = 0; iz < SIDE; iz++)
        for (int ix = 0; ix < SIDE; ix++) {
            size_t vi = nv + (size_t)iz * SIDE + (size_t)ix;
            verts[vi * 3 + 0] = (float)(2 * iz);
            verts[vi * 3 + 1] = 20.0f;
            verts[vi * 3 + 2] = (float)(124 + ix);
            uv[vi * 2 + 0] = (float)(124 + ix);
            uv[vi * 2 + 1] = (float)(2 * iz);
            gid[vi] = 0;
        }
    for (int iz = 0; iz < N; iz++)
        for (int ix = 0; ix < N; ix++) {
            size_t a = nv + (size_t)iz * SIDE + (size_t)ix;
            faces[nf * 3 + 0] = (int32_t)a;
            faces[nf * 3 + 1] = (int32_t)(a + 1);
            faces[nf * 3 + 2] = (int32_t)(a + SIDE);
            face_cube[nf++] = 0;
            faces[nf * 3 + 0] = (int32_t)(a + 1);
            faces[nf * 3 + 1] = (int32_t)(a + SIDE + 1);
            faces[nf * 3 + 2] = (int32_t)(a + SIDE);
            face_cube[nf++] = 0;
        }
    nv += slab_v;
    cube_voff[1] = nv;

    /* cube B: x in [128,132]; its x=128 column coincides with A's x=128
     * column in 3D, uv offset by 2 vox in u (the misregistration to weld) */
    for (int iz = 0; iz < SIDE; iz++)
        for (int ix = 0; ix < SIDE; ix++) {
            size_t vi = nv + (size_t)iz * SIDE + (size_t)ix;
            verts[vi * 3 + 0] = (float)(2 * iz);
            verts[vi * 3 + 1] = 20.0f;
            verts[vi * 3 + 2] = (float)(128 + ix);
            uv[vi * 2 + 0] = (float)(128 + ix) + 2.0f;   /* du = 2 */
            uv[vi * 2 + 1] = (float)(2 * iz);
            gid[vi] = 0;
        }
    for (int iz = 0; iz < N; iz++)
        for (int ix = 0; ix < N; ix++) {
            size_t a = nv + (size_t)iz * SIDE + (size_t)ix;
            faces[nf * 3 + 0] = (int32_t)a;
            faces[nf * 3 + 1] = (int32_t)(a + 1);
            faces[nf * 3 + 2] = (int32_t)(a + SIDE);
            face_cube[nf++] = 1;
            faces[nf * 3 + 0] = (int32_t)(a + 1);
            faces[nf * 3 + 1] = (int32_t)(a + SIDE + 1);
            faces[nf * 3 + 2] = (int32_t)(a + SIDE);
            face_cube[nf++] = 1;
        }
    nv += slab_v;
    cube_voff[2] = nv;

    /* mock other wrap in cube 0: 8 vox away in y -- must NOT weld */
    size_t wrap0 = nv;
    for (int iz = 0; iz < SIDE; iz++)
        for (int ix = 0; ix < SIDE; ix++) {
            size_t vi = nv + (size_t)iz * SIDE + (size_t)ix;
            verts[vi * 3 + 0] = (float)(2 * iz);
            verts[vi * 3 + 1] = 28.0f;
            verts[vi * 3 + 2] = (float)(124 + ix);
            uv[vi * 2 + 0] = (float)(124 + ix);
            uv[vi * 2 + 1] = (float)(2 * iz);
            gid[vi] = 1;
        }
    for (int iz = 0; iz < N; iz++)
        for (int ix = 0; ix < N; ix++) {
            size_t a = nv + (size_t)iz * SIDE + (size_t)ix;
            faces[nf * 3 + 0] = (int32_t)a;
            faces[nf * 3 + 1] = (int32_t)(a + 1);
            faces[nf * 3 + 2] = (int32_t)(a + SIDE);
            face_cube[nf++] = 0;
            faces[nf * 3 + 0] = (int32_t)(a + 1);
            faces[nf * 3 + 1] = (int32_t)(a + SIDE + 1);
            faces[nf * 3 + 2] = (int32_t)(a + SIDE);
            face_cube[nf++] = 0;
        }
    nv += slab_v;

    /* wrap patch belongs to cube A's range: extend cube A's voff instead of a
     * third cube (simplest: treat as 3 cubes with the wrap as cube 2 at the
     * same origin) */
    cube_voff[2] = wrap0;
    cube_voff[3] = nv;

    PieceSet ps;
    memset(&ps, 0, sizeof(ps));
    ps.verts = verts;
    ps.uv = uv;
    ps.gid = gid;
    ps.faces = faces;
    ps.face_cube = face_cube;
    ps.nv = nv;
    ps.nf = nf;
    ps.n_cubes = 3;
    ps.cube_voff = cube_voff;
    ps.cube_org = org;

    UVWeldOpts o;
    UVWeldOpts_default(&o);
    o.verbose = 0;

    UVWeldStats st;
    int rc = UVWeld_run(arena, &ps, &o, &st);
    uw_check(rc == 0, "rc", &fails);
    uw_check(st.n_pairs > 0 && st.n_clusters >= SIDE - 1,
             "boundary columns welded", &fails);
    /* welded rep uv = mean of the two sides: A had u=128, B had 130 -> 129 */
    int found_mean = 0;
    for (size_t i = 0; i < cube_voff[1]; i++) {
        if (fabsf(verts[i * 3 + 2] - 128.0f) < 0.6f
            && fabs((double)uv[i * 2 + 0] - 129.0) < 0.6)
            found_mean = 1;
    }
    uw_check(found_mean, "rep uv is the member mean", &fails);
    /* faces referencing B's x=128 column now share verts with A: one
     * connected component across the seam among welded columns */
    size_t cross = 0;
    for (size_t f = 0; f < nf; f++) {
        if (face_cube[f] != 1) continue;
        for (int k = 0; k < 3; k++)
            if ((size_t)faces[f * 3 + k] < cube_voff[1]) cross++;
    }
    uw_check(cross > 0, "cube B faces now reference cube A verts", &fails);
    /* the mock wrap (8 vox away) stayed unwelded */
    int wrap_ok = 1;
    for (size_t i = wrap0; i < nv; i++)
        if (fabsf(verts[i * 3 + 1] - 28.0f) > 0.1f) wrap_ok = 0;
    uw_check(wrap_ok, "other wrap untouched (3D gate)", &fails);

    /* |du| > gate stays unwelded: rerun on fresh geometry with du 8 */
    {
        for (size_t i = cube_voff[1]; i < cube_voff[2]; i++)
            uv[i * 2 + 0] += 6.0f;   /* B's du now 8 > gate 6 */
        UVWeldStats st2;
        PieceSet ps2 = ps;
        int rc2 = UVWeld_run(arena, &ps2, &o, &st2);
        uw_check(rc2 == 0 && st2.rej_uv > 0,
                 "mis-registration tail rejected by |du| gate", &fails);
    }

    Arena_dispose(&arena);
    fprintf(stderr, "[uv_weld selftest] %s (%d failure%s)\n",
            fails == 0 ? "PASSED" : "FAILED", fails, fails == 1 ? "" : "s");
    return fails;
}
