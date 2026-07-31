/* relax_grid.c -- red-black banded UV relax. See relax_grid.h. */
#include "relax_grid.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "../common/tiff_io.h"
#include "../common/ves_platform.h"
#include "../flatten/fiber_field.h"

static void *rg_xmalloc(size_t n)
{
    void *p = malloc(n > 0 ? n : 1);
    if (p == NULL) { fprintf(stderr, "relax_grid: OOM %zu\n", n); exit(1); }
    return p;
}

static void *rg_xcalloc(size_t c, size_t s)
{
    void *p = calloc(c > 0 ? c : 1, s);
    if (p == NULL) { fprintf(stderr, "relax_grid: OOM\n"); exit(1); }
    return p;
}

static int rg_cmp_i32(const void *a, const void *b)
{
    int32_t x = *(const int32_t *)a, y = *(const int32_t *)b;
    return x < y ? -1 : (x > y ? 1 : 0);
}

void RelaxGridOpts_default(RelaxGridOpts *o)
{
    assert(o);
    memset(o, 0, sizeof(*o));
    RibbonRelax_defaults(&o->base);
    o->base.sweeps = 16;
    o->band_px = 32768;
    o->halo_px = 512;
    o->rounds = 1;
    o->threads = 16;
    o->fib_grad_sigma = 1.5;
    o->fib_tensor_sigma = 6.0;
    o->fib_axis_tol = 15.0;
    o->verbose = 0;
}

typedef struct RgBandStat {
    int run, reverted, no_fiber;
    size_t moved;
    double e0, e1, s0, s1, mean_disp, max_disp;
} RgBandStat;

/* binary search local id of global vert g in sorted l2g[nlv] */
static int32_t rg_local(const int32_t *l2g, size_t nlv, int32_t g)
{
    size_t lo = 0, hi = nlv;
    while (lo < hi) {
        size_t m = (lo + hi) / 2;
        if (l2g[m] < g) lo = m + 1;
        else hi = m;
    }
    return (int32_t)lo;   /* caller guarantees presence */
}

static void rg_band_run(PieceSet *ps, const char *rawtex_tif,
                        long X0, long Y0, double du, double dv,
                        size_t W, size_t H,
                        size_t b0, size_t bw,
                        const int32_t *bfaces, size_t nbf,
                        const int32_t *owner_band, int32_t band_id,
                        const RelaxGridOpts *o, RgBandStat *st)
{
    memset(st, 0, sizeof(*st));
    if (nbf == 0) return;
    Arena_T ar = Arena_new();

    /* fiber window: band + halo + tensor support */
    long m = (long)o->halo_px + (long)ceil(3.0 * o->fib_tensor_sigma) + 2;
    long wx0 = (long)b0 - m;
    if (wx0 < 0) wx0 = 0;
    long wx1 = (long)(b0 + bw) + m;
    if (wx1 > (long)W) wx1 = (long)W;
    int ww = (int)(wx1 - wx0);

    FiberField fib;
    const FiberField *fibp = NULL;
    if (rawtex_tif != NULL && ww >= 16) {
        TiffRowReader rd;
        if (TiffIO_rows_open(ar, rawtex_tif, &rd) == 0) {
            if (rd.W == (int)W && rd.H == (int)H) {
                uint8_t *img = (uint8_t *)ARENA_ALLOC(ar,
                                    (size_t)ww * H);
                int ok = 1;
                for (int y = 0; y < (int)H && ok; y++)
                    if (TiffIO_rows_read(&rd, y, (int)wx0, ww,
                                         img + (size_t)y * (size_t)ww) != 0)
                        ok = 0;
                if (ok && FiberField_compute(ar, img, ww, (int)H,
                                             o->fib_grad_sigma,
                                             o->fib_tensor_sigma,
                                             o->fib_axis_tol, 4, &fib) == 0)
                    fibp = &fib;
            }
            TiffIO_rows_close(&rd);
        }
    }
    if (fibp == NULL) st->no_fiber = 1;

    /* submesh verts: sorted unique global ids over the band's faces */
    size_t nref = nbf * 3;
    int32_t *l2g = (int32_t *)rg_xmalloc(nref * sizeof(int32_t));
    for (size_t f = 0; f < nbf; f++)
        for (int k = 0; k < 3; k++)
            l2g[f * 3 + k] = ps->faces[(size_t)bfaces[f] * 3 + k];
    qsort(l2g, nref, sizeof(int32_t), rg_cmp_i32);
    size_t nlv = 0;
    for (size_t i = 0; i < nref; i++)
        if (nlv == 0 || l2g[i] != l2g[nlv - 1]) l2g[nlv++] = l2g[i];

    /* band-local frame: shift uv by the fiber-window origin so texel col =
     * u_loc/du lines up with the rect AND float32 stays precise */
    double u_org = (double)(X0 + wx0) * du;
    double v_org = (double)Y0 * dv;

    float *sv = (float *)rg_xmalloc(nlv * 3 * sizeof(float));
    float *suv = (float *)rg_xmalloc(nlv * 2 * sizeof(float));
    float *suv_out = (float *)rg_xmalloc(nlv * 2 * sizeof(float));
    uint8_t *pin = (uint8_t *)rg_xcalloc(nlv, 1);
    for (size_t i = 0; i < nlv; i++) {
        size_t g = (size_t)l2g[i];
        sv[i * 3 + 0] = ps->verts[g * 3 + 0];
        sv[i * 3 + 1] = ps->verts[g * 3 + 1];
        sv[i * 3 + 2] = ps->verts[g * 3 + 2];
        suv[i * 2 + 0] = (float)((double)ps->uv[g * 2 + 0] - u_org);
        suv[i * 2 + 1] = (float)((double)ps->uv[g * 2 + 1] - v_org);
        pin[i] = owner_band[g] != band_id;
    }
    int32_t *sf = (int32_t *)rg_xmalloc(nbf * 3 * sizeof(int32_t));
    for (size_t f = 0; f < nbf; f++)
        for (int k = 0; k < 3; k++)
            sf[f * 3 + k] = rg_local(l2g, nlv,
                                     ps->faces[(size_t)bfaces[f] * 3 + k]);

    RibbonRelaxOpts ro = o->base;
    ro.pin = pin;
    RibbonRelaxStats rst;
    memset(&rst, 0, sizeof(rst));
    int rc = RibbonRelax_run(ar, sv, nlv, sf, nbf, suv, NULL, fibp, du, dv,
                             &ro, suv_out, &rst);
    st->run = 1;
    st->e0 = rst.energy_before;
    st->e1 = rst.energy_after;
    st->s0 = rst.stretch_mean_before;
    st->s1 = rst.stretch_mean_after;
    st->mean_disp = rst.mean_disp;
    st->max_disp = rst.max_disp;
    if (rc != 0 || rst.reverted) {
        st->reverted = 1;
    } else {
        for (size_t i = 0; i < nlv; i++) {
            if (pin[i]) continue;
            size_t g = (size_t)l2g[i];
            double nu = (double)suv_out[i * 2 + 0] + u_org;
            double nv2 = (double)suv_out[i * 2 + 1] + v_org;
            if (fabs(nu - (double)ps->uv[g * 2 + 0]) > 1e-6
                || fabs(nv2 - (double)ps->uv[g * 2 + 1]) > 1e-6)
                st->moved++;
            ps->uv[g * 2 + 0] = (float)nu;
            ps->uv[g * 2 + 1] = (float)nv2;
        }
    }

    free(l2g);
    free(sv);
    free(suv);
    free(suv_out);
    free(pin);
    free(sf);
    Arena_dispose(&ar);
}

int RelaxGrid_run(Arena_T arena, PieceSet *ps, const char *rawtex_tif,
                  double u0, double v0, double du, double dv,
                  const RelaxGridOpts *opts, RelaxGridStats *out)
{
    assert(arena && ps && opts && out);
    memset(out, 0, sizeof(*out));
    if (ps->nv == 0 || ps->nf == 0 || du <= 0.0 || dv <= 0.0)
        return 0;
    const RelaxGridOpts *o = opts;
    double t_start = ves_clock_sec();

    /* canvas mapping (identical formula to ScrollRaster_run) */
    long X0 = (long)floor(u0 / du + 0.5);
    long Y0 = (long)floor(v0 / dv + 0.5);
    long X1 = (long)floor(ps->u_max / du + 0.5);
    long Y1 = (long)floor(ps->v_max / dv + 0.5);
    size_t W = (size_t)(X1 - X0) + 1;
    size_t H = (size_t)(Y1 - Y0) + 1;
    size_t B = (o->band_px == 0 || o->band_px >= W) ? W : o->band_px;
    size_t n_bands = (W + B - 1) / B;
    out->n_bands = n_bands;

    /* owner band per vertex (kept-face verts only; others -1) */
    int32_t *owner = (int32_t *)rg_xmalloc(ps->nv * sizeof(int32_t));
    for (size_t i = 0; i < ps->nv; i++) owner[i] = -1;
    for (size_t f = 0; f < ps->nf; f++) {
        for (int k = 0; k < 3; k++) {
            size_t g = (size_t)ps->faces[f * 3 + k];
            if (owner[g] >= 0) continue;
            long px = (long)floor((double)ps->uv[g * 2 + 0] / du + 0.5) - X0;
            if (px < 0) px = 0;
            if (px >= (long)W) px = (long)W - 1;
            owner[g] = (int32_t)((size_t)px / B);
        }
    }

    /* face -> band CSR, halo-expanded (a face enrolls in every band whose
     * halo-extended window it touches) */
    int32_t *boff = (int32_t *)rg_xcalloc(n_bands + 2, sizeof(int32_t));
    long *fb0 = (long *)rg_xmalloc(ps->nf * sizeof(long));
    long *fb1 = (long *)rg_xmalloc(ps->nf * sizeof(long));
    for (size_t f = 0; f < ps->nf; f++) {
        double ulo = 1e300, uhi = -1e300;
        for (int k = 0; k < 3; k++) {
            double u = (double)ps->uv[(size_t)ps->faces[f * 3 + k] * 2 + 0];
            if (u < ulo) ulo = u;
            if (u > uhi) uhi = u;
        }
        long px0 = (long)floor(ulo / du + 0.5) - X0 - (long)o->halo_px;
        long px1 = (long)floor(uhi / du + 0.5) - X0 + (long)o->halo_px;
        if (px0 < 0) px0 = 0;
        if (px1 >= (long)W) px1 = (long)W - 1;
        if (px1 < px0) { fb0[f] = 1; fb1[f] = 0; continue; }
        fb0[f] = px0 / (long)B;
        fb1[f] = px1 / (long)B;
        for (long b = fb0[f]; b <= fb1[f]; b++) boff[b + 1]++;
    }
    for (size_t b = 0; b < n_bands; b++) boff[b + 1] += boff[b];
    int32_t *bface = (int32_t *)rg_xmalloc(((size_t)boff[n_bands] + 1)
                                           * sizeof(int32_t));
    {
        int32_t *cur = (int32_t *)rg_xmalloc((n_bands + 1) * sizeof(int32_t));
        memcpy(cur, boff, n_bands * sizeof(int32_t));
        for (size_t f = 0; f < ps->nf; f++)
            for (long b = fb0[f]; b <= fb1[f]; b++)
                bface[cur[b]++] = (int32_t)f;
        free(cur);
    }
    free(fb0);
    free(fb1);

    RgBandStat *bst = (RgBandStat *)rg_xcalloc(n_bands, sizeof(RgBandStat));
#ifdef _OPENMP
    ves_omp_set_threads(o->threads > 0 ? o->threads : 16);
#endif
    for (int round = 0; round < (o->rounds > 0 ? o->rounds : 1); round++) {
        for (int parity = 0; parity < 2; parity++) {
            int n = (int)n_bands;
            int b = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
            for (b = parity; b < n; b += 2) {
                size_t bx0 = (size_t)b * B;
                size_t bw = (bx0 + B <= W) ? B : W - bx0;
                rg_band_run(ps, rawtex_tif, X0, Y0, du, dv, W, H, bx0, bw,
                            &bface[boff[b]],
                            (size_t)(boff[b + 1] - boff[b]),
                            owner, b, o, &bst[b]);
            }
        }
    }

    for (size_t b = 0; b < n_bands; b++) {
        const RgBandStat *s = &bst[b];
        if (!s->run) continue;
        out->bands_run++;
        if (s->reverted) out->bands_reverted++;
        if (s->no_fiber) out->bands_no_fiber++;
        out->n_moved += s->moved;
        out->energy_before += s->e0;
        out->energy_after += s->e1;
        out->stretch_before += s->s0;
        out->stretch_after += s->s1;
        out->mean_disp += s->mean_disp;
        if (s->max_disp > out->max_disp) out->max_disp = s->max_disp;
    }
    if (out->bands_run > 0) {
        out->stretch_before /= (double)out->bands_run;
        out->stretch_after /= (double)out->bands_run;
        out->mean_disp /= (double)out->bands_run;
    }
    out->seconds = ves_clock_sec() - t_start;

    if (o->verbose)
        fprintf(stderr, "[relax_grid] bands=%zu run=%zu reverted=%zu "
                "no_fiber=%zu moved=%zu E %.3g->%.3g stretch %.3f->%.3f "
                "disp mean=%.2f max=%.2f (%.1fs)\n",
                out->n_bands, out->bands_run, out->bands_reverted,
                out->bands_no_fiber, out->n_moved, out->energy_before,
                out->energy_after, out->stretch_before, out->stretch_after,
                out->mean_disp, out->max_disp, out->seconds);

    free(owner);
    free(boff);
    free(bface);
    free(bst);
    return 0;
}

/* ============================================================================
 * Self-test.
 * ==========================================================================*/

static int rg_check(int cond, const char *what, int *fails)
{
    if (!cond) {
        fprintf(stderr, "[relax_grid selftest]   FAIL: %s\n", what);
        (*fails)++;
    }
    return cond;
}

/* welded 2-patch sheet: one 24x6 grid (single cube), flat in 3D (u tracks x
 * 1:1), with a +2 vox u-step injected at x>=24 -- the quilting gutter after
 * a weld. Isometry-only relax (no rawtex) should diffuse the step. */
static void rg_fixture(float *verts, float *uv, int32_t *faces,
                       int32_t *face_cube, size_t *pnv, size_t *pnf,
                       int nx, int nz)
{
    size_t nv = 0, nf = 0;
    int sx = nx + 1, sz = nz + 1;
    for (int iz = 0; iz < sz; iz++) {
        for (int ix = 0; ix < sx; ix++) {
            size_t vi = (size_t)iz * (size_t)sx + (size_t)ix;
            verts[vi * 3 + 0] = (float)(2 * iz);
            verts[vi * 3 + 1] = 20.0f;
            verts[vi * 3 + 2] = (float)(2 * ix);
            double u = 2.0 * ix + (ix >= nx / 2 ? 2.0 : 0.0);   /* the step */
            uv[vi * 2 + 0] = (float)u;
            uv[vi * 2 + 1] = (float)(2 * iz);
            nv++;
        }
    }
    for (int iz = 0; iz < nz; iz++) {
        for (int ix = 0; ix < nx; ix++) {
            size_t a = (size_t)iz * (size_t)sx + (size_t)ix;
            faces[nf * 3 + 0] = (int32_t)a;
            faces[nf * 3 + 1] = (int32_t)(a + 1);
            faces[nf * 3 + 2] = (int32_t)(a + (size_t)sx);
            face_cube[nf++] = 0;
            faces[nf * 3 + 0] = (int32_t)(a + 1);
            faces[nf * 3 + 1] = (int32_t)(a + (size_t)sx + 1);
            faces[nf * 3 + 2] = (int32_t)(a + (size_t)sx);
            face_cube[nf++] = 0;
        }
    }
    *pnv = nv;
    *pnf = nf;
}

/* max |u(ix=hi) - u(ix=lo)| across the step columns */
static double rg_step_gap(const float *uv, int nx, int nz)
{
    int sx = nx + 1;
    double gap = 0.0;
    for (int iz = 0; iz <= nz; iz++) {
        size_t a = (size_t)iz * (size_t)sx + (size_t)(nx / 2 - 1);
        size_t b = a + 1;
        double d = fabs((double)uv[b * 2 + 0] - (double)uv[a * 2 + 0] - 2.0);
        if (d > gap) gap = d;
    }
    return gap;
}

int RelaxGrid_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();

    enum { NX = 24, NZ = 6 };
    size_t nvcap = (NX + 1) * (NZ + 1);
    float *verts = (float *)rg_xmalloc(nvcap * 3 * sizeof(float));
    float *uv = (float *)rg_xmalloc(nvcap * 2 * sizeof(float));
    float *normals = (float *)rg_xcalloc(nvcap, 3 * sizeof(float));
    int32_t *gid = (int32_t *)rg_xcalloc(nvcap, sizeof(int32_t));
    int32_t *faces = (int32_t *)rg_xmalloc((size_t)NX * NZ * 2 * 3
                                           * sizeof(int32_t));
    int32_t *face_cube = (int32_t *)rg_xmalloc((size_t)NX * NZ * 2
                                               * sizeof(int32_t));
    size_t cube_voff[2] = { 0, 0 };
    long org[1][3] = { { 0, 0, 0 } };

    RelaxGridOpts o;
    RelaxGridOpts_default(&o);
    o.threads = 2;
    o.halo_px = 4;
    o.base.sweeps = 24;
    o.verbose = 0;

    double gap_single = 0.0, gap_banded = 0.0;

    /* t1: single band closes the injected 2-vox step */
    {
        size_t nv = 0, nf = 0;
        rg_fixture(verts, uv, faces, face_cube, &nv, &nf, NX, NZ);
        cube_voff[1] = nv;
        PieceSet ps;
        memset(&ps, 0, sizeof(ps));
        ps.verts = verts; ps.uv = uv; ps.normals = normals; ps.gid = gid;
        ps.faces = faces; ps.face_cube = face_cube;
        ps.nv = nv; ps.nf = nf; ps.n_cubes = 1;
        ps.cube_voff = cube_voff; ps.cube_org = org;
        ps.u_min = 0; ps.u_max = 2.0 * NX + 2.0;
        ps.v_min = 0; ps.v_max = 2.0 * NZ;
        double g0 = rg_step_gap(uv, NX, NZ);
        RelaxGridOpts o1 = o;
        o1.band_px = 0;   /* single band */
        RelaxGridStats st;
        int rc = RelaxGrid_run(arena, &ps, NULL, ps.u_min, ps.v_min,
                               1.0, 1.0, &o1, &st);
        gap_single = rg_step_gap(uv, NX, NZ);
        rg_check(rc == 0 && st.bands_run == 1 && st.bands_reverted == 0,
                 "t1 single band runs", &fails);
        rg_check(st.n_moved > 0, "t1 verts moved", &fails);
        rg_check(gap_single < 0.6 * g0, "t1 step gap shrinks", &fails);
        rg_check(st.stretch_after <= st.stretch_before + 1e-9,
                 "t1 stretch non-increasing", &fails);
    }

    /* t2: banded red-black reaches a comparable gap (parity with t1) */
    {
        size_t nv = 0, nf = 0;
        rg_fixture(verts, uv, faces, face_cube, &nv, &nf, NX, NZ);
        cube_voff[1] = nv;
        PieceSet ps;
        memset(&ps, 0, sizeof(ps));
        ps.verts = verts; ps.uv = uv; ps.normals = normals; ps.gid = gid;
        ps.faces = faces; ps.face_cube = face_cube;
        ps.nv = nv; ps.nf = nf; ps.n_cubes = 1;
        ps.cube_voff = cube_voff; ps.cube_org = org;
        ps.u_min = 0; ps.u_max = 2.0 * NX + 2.0;
        ps.v_min = 0; ps.v_max = 2.0 * NZ;
        RelaxGridOpts o2 = o;
        o2.band_px = 16;    /* ~3 bands over ~50 px */
        o2.rounds = 2;
        RelaxGridStats st;
        int rc = RelaxGrid_run(arena, &ps, NULL, ps.u_min, ps.v_min,
                               1.0, 1.0, &o2, &st);
        gap_banded = rg_step_gap(uv, NX, NZ);
        rg_check(rc == 0 && st.n_bands >= 3 && st.bands_reverted == 0,
                 "t2 banded runs", &fails);
        rg_check(gap_banded < 1.2, "t2 banded closes the gap too", &fails);
        fprintf(stderr, "[relax_grid selftest] gap 2.0 -> single %.2f / "
                "banded %.2f\n", gap_single, gap_banded);
    }

    /* t3: RibbonRelax pin honored (pin everything -> nothing moves) */
    {
        size_t nv = 0, nf = 0;
        rg_fixture(verts, uv, faces, face_cube, &nv, &nf, NX, NZ);
        float *uv_out = (float *)rg_xmalloc(nv * 2 * sizeof(float));
        uint8_t *pin = (uint8_t *)rg_xmalloc(nv);
        memset(pin, 1, nv);
        RibbonRelaxOpts ro;
        RibbonRelax_defaults(&ro);
        ro.pin = pin;
        RibbonRelaxStats rst;
        int rc = RibbonRelax_run(arena, verts, nv, faces, nf, uv, NULL,
                                 NULL, 1.0, 1.0, &ro, uv_out, &rst);
        int same = 1;
        for (size_t i = 0; i < nv * 2; i++)
            if (fabsf(uv_out[i] - uv[i]) > 1e-6f) same = 0;
        rg_check(rc == 0 && same && rst.n_moved == 0,
                 "t3 full pin freezes uv", &fails);
        free(uv_out);
        free(pin);
    }

    /* t4: empty */
    {
        PieceSet ps;
        memset(&ps, 0, sizeof(ps));
        RelaxGridStats st;
        rg_check(RelaxGrid_run(arena, &ps, NULL, 0, 0, 1.0, 1.0, &o, &st)
                 == 0, "t4 empty clean", &fails);
    }

    free(verts); free(uv); free(normals); free(gid);
    free(faces); free(face_cube);
    Arena_dispose(&arena);
    fprintf(stderr, "[relax_grid selftest] %s (%d failure%s)\n",
            fails == 0 ? "PASSED" : "FAILED", fails, fails == 1 ? "" : "s");
    return fails;
}
