/* closeup_surface.c -- see closeup_surface.h. */
#include "closeup_surface.h"
#include "sparse_solve.h"
#include "slim_refine.h"
#include "../common/arena.h"
#include "../common/pipeline_constants.h"
#include "../common/ves_platform.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/ves_omp.h"

#define CS_INF 3.0e38f

/* ------------------------------------------------------------ observations */

void CloseupObsGate_default(CloseupObsGate *g)
{
    g->stretch_ratio = CLOSEUP_STRETCH_RATIO;
    g->stretch_floor = CLOSEUP_STRETCH_FLOOR;
    g->max_edge3d = CLOSEUP_MAX_EDGE3D_VOX;
    g->min_scale = CLOSEUP_OBS_MIN_SCALE;
    g->max_scale = CLOSEUP_OBS_MAX_SCALE;
}

void CloseupObs_free(CloseupObs *o)
{
    if (o == NULL) return;
    free(o->uv); free(o->p); free(o->w); free(o->chart);
    memset(o, 0, sizeof *o);
}

static int cs_obs_reserve(CloseupObs *o, size_t cap)
{
    float *uv = NULL, *p = NULL, *w = NULL;
    int32_t *c = NULL;
    if (cap <= o->cap) return 0;
    uv = (float *)realloc(o->uv, cap * 2 * sizeof *uv);
    if (uv == NULL) return -1;
    o->uv = uv;
    p = (float *)realloc(o->p, cap * 3 * sizeof *p);
    if (p == NULL) return -1;
    o->p = p;
    w = (float *)realloc(o->w, cap * sizeof *w);
    if (w == NULL) return -1;
    o->w = w;
    c = (int32_t *)realloc(o->chart, cap * sizeof *c);
    if (c == NULL) return -1;
    o->chart = c;
    o->cap = cap;
    return 0;
}

/* The face's pixel-centre box in the frame and its gate verdict: 1 = paint, 0 = gated, -1 = degenerate. */
static int cs_face_box(const float *verts, const double *uv, const int32_t *faces, size_t f, const CloseupObsGate *gt,
                       long *x0, long *x1, long *y0, long *y1)
{
    size_t vi[3];
    double lox = 0.0, hix = 0.0, loy = 0.0, hiy = 0.0, a2 = 0.0;
    int e = 0;
    for (int k = 0; k < 3; k++) if (faces[f*3+(size_t)k] < 0) return 0;
    vi[0] = (size_t)faces[f * 3]; vi[1] = (size_t)faces[f * 3 + 1]; vi[2] = (size_t)faces[f * 3 + 2];
    for (int k = 0; k < 3; k++) {
        for (int d = 0; d < 3; d++) if (!isfinite(verts[3*vi[k]+(size_t)d])) return 0;
        for (int d = 0; d < 2; d++) {
            double value = uv[2*vi[k]+(size_t)d];
            if (!isfinite(value) || value < LONG_MIN+2.0 || value > LONG_MAX-2.0) return 0;
        }
    }
    a2 = (uv[vi[1] * 2] - uv[vi[0] * 2]) * (uv[vi[2] * 2 + 1] - uv[vi[0] * 2 + 1]) -
         (uv[vi[1] * 2 + 1] - uv[vi[0] * 2 + 1]) * (uv[vi[2] * 2] - uv[vi[0] * 2]);
    if (!isfinite(a2) || !(fabs(a2) > 1e-12)) return -1;
    {
        double edge[3], across[3], cross[3], length2 = 0, dot = 0;
        for (int d = 0; d < 3; d++) {
            edge[d] = (double)verts[3*vi[1]+(size_t)d]-verts[3*vi[0]+(size_t)d];
            across[d] = (double)verts[3*vi[2]+(size_t)d]-verts[3*vi[0]+(size_t)d];
            length2 += edge[d]*edge[d]; dot += edge[d]*across[d];
        }
        for (int d = 0; d < 3; d++) cross[d] = edge[(d+1)%3]*across[(d+2)%3]-edge[(d+2)%3]*across[(d+1)%3];
        double length = sqrt(length2), area = hypot(hypot(cross[0],cross[1]),cross[2]);
        if (!(length > 0) || !(area > 0)) return 0;
        double along = dot/length, height = area/length;
        double a = (uv[2*vi[1]]-uv[2*vi[0]])/length;
        double c = (uv[2*vi[1]+1]-uv[2*vi[0]+1])/length;
        double b = (uv[2*vi[2]]-uv[2*vi[0]]-along*a)/height;
        double d = (uv[2*vi[2]+1]-uv[2*vi[0]+1]-along*c)/height;
        double hi = .5*(hypot(a+d,c-b)+hypot(a-d,c+b)), lo = hi > 0 ? fabs(a*d-b*c)/hi : 0;
        if (!isfinite(hi) || !isfinite(lo) || !(lo > 0) ||
            (gt->min_scale > 0 && lo < gt->min_scale) || (gt->max_scale > 0 && hi > gt->max_scale)) return 0;
    }
    for (e = 0; e < 3; e++) {
        size_t p = vi[e], q = vi[(e + 1) % 3];
        double d3 = 0.0, du = uv[q * 2] - uv[p * 2], dv = uv[q * 2 + 1] - uv[p * 2 + 1], e2 = 0.0;
        int k = 0;
        for (k = 0; k < 3; k++) {
            double dd = (double)verts[q * 3 + (size_t)k] - (double)verts[p * 3 + (size_t)k];
            d3 += dd * dd;
        }
        d3 = sqrt(d3);
        e2 = sqrt(du * du + dv * dv);
        if (gt->max_edge3d > 0.0 && d3 > gt->max_edge3d) return 0;
        if (gt->stretch_ratio > 0.0 &&
            e2 > (gt->stretch_ratio * d3 > gt->stretch_floor ? gt->stretch_ratio * d3 : gt->stretch_floor))
            return 0;
        if (!isfinite(du) || !isfinite(dv)) return -1;
    }
    lox = fmin(uv[vi[0] * 2], fmin(uv[vi[1] * 2], uv[vi[2] * 2]));
    hix = fmax(uv[vi[0] * 2], fmax(uv[vi[1] * 2], uv[vi[2] * 2]));
    loy = fmin(uv[vi[0] * 2 + 1], fmin(uv[vi[1] * 2 + 1], uv[vi[2] * 2 + 1]));
    hiy = fmax(uv[vi[0] * 2 + 1], fmax(uv[vi[1] * 2 + 1], uv[vi[2] * 2 + 1]));
    /* centres at integer + 0.5 */
    *x0 = (long)ceil(lox - 0.5 - 1e-9); *x1 = (long)floor(hix - 0.5 + 1e-9);
    *y0 = (long)ceil(loy - 0.5 - 1e-9); *y1 = (long)floor(hiy - 0.5 + 1e-9);
    return 1;
}

/* Visit (count == NULL: write) the pixel centres of face f; returns the count. */
static size_t cs_face_raster(const float *verts, const double *uv, const int32_t *faces, size_t f, long x0, long x1,
                             long y0, long y1, CloseupObs *o, size_t at, int32_t chart)
{
    size_t vi[3], n = 0;
    double ua, va, ub, vb, uc, vc, A2;
    long x = 0, y = 0;
    vi[0] = (size_t)faces[f * 3]; vi[1] = (size_t)faces[f * 3 + 1]; vi[2] = (size_t)faces[f * 3 + 2];
    ua = uv[vi[0] * 2]; va = uv[vi[0] * 2 + 1];
    ub = uv[vi[1] * 2]; vb = uv[vi[1] * 2 + 1];
    uc = uv[vi[2] * 2]; vc = uv[vi[2] * 2 + 1];
    A2 = (ub - ua) * (vc - va) - (vb - va) * (uc - ua);
    for (y = y0; y <= y1; y++) {
        double py = (double)y + 0.5;
        for (x = x0; x <= x1; x++) {
            double px = (double)x + 0.5;
            double l1 = ((uc - px) * (va - py) - (vc - py) * (ua - px)) / A2;   /* weight of b */
            double l2 = ((ua - px) * (vb - py) - (va - py) * (ub - px)) / A2;   /* weight of c */
            double l0 = 1.0 - l1 - l2;
            int k = 0;
            if (l0 < -1e-9 || l1 < -1e-9 || l2 < -1e-9) continue;
            if (o != NULL) {
                size_t w = at + n;
                o->uv[w * 2] = (float)px; o->uv[w * 2 + 1] = (float)py;
                for (k = 0; k < 3; k++)
                    o->p[w * 3 + (size_t)k] = (float)(l0 * verts[vi[0] * 3 + (size_t)k] +
                                                      l1 * verts[vi[1] * 3 + (size_t)k] +
                                                      l2 * verts[vi[2] * 3 + (size_t)k]);
                o->w[w] = 1.0f;
                o->chart[w] = chart;
            }
            n++;
        }
    }
    return n;
}

int CloseupObs_add_faces(CloseupObs *o, const float *verts, const double *uv, const int32_t *faces, size_t nf,
                         const uint8_t *face_use, const int32_t *face_chart, const CloseupObsGate *gate,
                         size_t *n_used, size_t *n_gated)
{
    size_t *cnt = NULL, total = 0, used = 0, gated = 0, f = 0, base = o->n;
    long lf = 0;
    CloseupObsGate g0;
    if (gate == NULL) { CloseupObsGate_default(&g0); gate = &g0; }
    if (nf == 0) { if (n_used) *n_used = 0; if (n_gated) *n_gated = 0; return 0; }
    cnt = (size_t *)calloc(nf + 1, sizeof *cnt);
    if (cnt == NULL) return -1;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4096)
#endif
    for (lf = 0; lf < (long)nf; lf++) {
        long x0 = 0, x1 = 0, y0 = 0, y1 = 0;
        int ok = 0;
        cnt[lf] = 0;
        if (face_use != NULL && !face_use[lf]) continue;
        ok = cs_face_box(verts, uv, faces, (size_t)lf, gate, &x0, &x1, &y0, &y1);
        if (ok == 0) { cnt[lf] = (size_t)-1; continue; }
        if (ok < 0 || x0 > x1 || y0 > y1) continue;
        cnt[lf] = cs_face_raster(verts, uv, faces, (size_t)lf, x0, x1, y0, y1, NULL, 0, -1);
    }
    for (f = 0; f < nf; f++) {
        size_t c = cnt[f];
        /* Every entry is an exclusive prefix, including rejected faces. A zero
         * here makes the preceding first face look empty to the write pass,
         * leaving its counted observations uninitialized. */
        if (c == (size_t)-1) { gated++; c = 0; }
        if (c > 0) used++;
        cnt[f] = total;
        total += c;
    }
    cnt[nf] = total;
    if (cs_obs_reserve(o, base + total) != 0) { free(cnt); return -1; }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4096)
#endif
    for (lf = 0; lf < (long)nf; lf++) {
        long x0 = 0, x1 = 0, y0 = 0, y1 = 0;
        if (cnt[lf + 1] == cnt[lf]) continue;
        if (cs_face_box(verts, uv, faces, (size_t)lf, gate, &x0, &x1, &y0, &y1) != 1) continue;
        (void)cs_face_raster(verts, uv, faces, (size_t)lf, x0, x1, y0, y1, o, base + cnt[lf],
                             face_chart != NULL ? face_chart[lf] : -1);
    }
    o->n = base + total;
    free(cnt);
    if (n_used) *n_used = used;
    if (n_gated) *n_gated = gated;
    return 0;
}

/* ------------------------------------------------------------------ grid */

/* The node arrays of a grid whose ni, nj are set: pos NaN, dist INFINITY, chart and component -1, the rest 0. */
static int cs_grid_arrays(CloseupGrid *g)
{
    size_t nn = g->ni * g->nj, k = 0;
    if (g->ni > ((size_t)1 << 22) || g->nj > ((size_t)1 << 22) || nn > ((size_t)1 << 31) || nn == 0) return -1;
    g->pos = (float *)malloc(nn * 3 * sizeof *g->pos);
    g->tu = (float *)calloc(nn * 3, sizeof *g->tu);
    g->tv = (float *)calloc(nn * 3, sizeof *g->tv);
    g->nrm = (float *)calloc(nn * 3, sizeof *g->nrm);
    g->dist = (float *)malloc(nn * sizeof *g->dist);
    g->support = (float *)calloc(nn, sizeof *g->support);
    g->chart = (int32_t *)malloc(nn * sizeof *g->chart);
    g->component = (int32_t *)malloc(nn * sizeof *g->component);
    g->state = (uint8_t *)calloc(nn, 1);
    if (g->pos == NULL || g->tu == NULL || g->tv == NULL || g->nrm == NULL || g->dist == NULL ||
        g->support == NULL || g->chart == NULL || g->component == NULL || g->state == NULL) {
        CloseupGrid_free(g);
        return -1;
    }
    for (k = 0; k < nn * 3; k++) g->pos[k] = NAN;
    for (k = 0; k < nn; k++) { g->dist[k] = INFINITY; g->chart[k] = -1; g->component[k] = -1; }
    return 0;
}

int CloseupGrid_alloc(CloseupGrid *g, double umin, double umax, double vmin, double vmax, double h)
{
    memset(g, 0, sizeof *g);
    if (!(h > 0.0) || !(umax >= umin) || !(vmax >= vmin) || !isfinite(umin) || !isfinite(vmin)) return -1;
    g->h = h;
    g->u0 = floor(umin / h) * h;
    g->v0 = floor(vmin / h) * h;
    g->ni = (size_t)ceil((umax - g->u0) / h) + 1;
    g->nj = (size_t)ceil((vmax - g->v0) / h) + 1;
    return cs_grid_arrays(g);
}

void CloseupGrid_free(CloseupGrid *g)
{
    if (g == NULL) return;
    free(g->pos); free(g->tu); free(g->tv); free(g->nrm); free(g->dist); free(g->support); free(g->chart);
    free(g->state);
    free(g->component);
    memset(g, 0, sizeof *g);
}

#define CS_GRID_MAGIC "VESCGRD1"

int CloseupGrid_save(const char *path, const CloseupGrid *g, const int32_t *meta, size_t nmeta)
{
    const size_t nn = g->ni * g->nj;
    size_t k = 0, nvalid = 0, w = 0;
    uint64_t hdr[3];
    double fr[3];
    uint32_t ver = 1, nm = (uint32_t)nmeta, *idx = NULL;
    uint8_t *st = NULL;
    int32_t *ch = NULL, *co = NULL;
    float *ps = NULL, *nr = NULL;
    FILE *fp = NULL;
    int rc = -1;
    if (nmeta > 64 || (nmeta > 0 && meta == NULL)) return -1;
    for (k = 0; k < nn; k++) nvalid += g->state[k] != CLOSEUP_NODE_NONE;
    idx = (uint32_t *)malloc((nvalid ? nvalid : 1) * sizeof *idx);
    st = (uint8_t *)malloc(nvalid ? nvalid : 1);
    ch = (int32_t *)malloc((nvalid ? nvalid : 1) * sizeof *ch);
    co = (int32_t *)malloc((nvalid ? nvalid : 1) * sizeof *co);
    ps = (float *)malloc((nvalid ? nvalid : 1) * 3 * sizeof *ps);
    nr = (float *)malloc((nvalid ? nvalid : 1) * 3 * sizeof *nr);
    if (idx == NULL || st == NULL || ch == NULL || co == NULL || ps == NULL || nr == NULL) goto done;
    for (k = 0; k < nn; k++) {
        if (g->state[k] == CLOSEUP_NODE_NONE) continue;
        idx[w] = (uint32_t)k; st[w] = g->state[k]; ch[w] = g->chart[k]; co[w] = g->component[k];
        memcpy(&ps[w * 3], &g->pos[k * 3], 3 * sizeof *ps);
        memcpy(&nr[w * 3], &g->nrm[k * 3], 3 * sizeof *nr);
        w++;
    }
    hdr[0] = (uint64_t)g->ni; hdr[1] = (uint64_t)g->nj; hdr[2] = (uint64_t)nvalid;
    fr[0] = g->u0; fr[1] = g->v0; fr[2] = g->h;
    fp = fopen(path, "wb");
    if (fp == NULL) goto done;
    if (fwrite(CS_GRID_MAGIC, 1, 8, fp) != 8 || fwrite(&ver, sizeof ver, 1, fp) != 1 || fwrite(&nm, sizeof nm, 1, fp) != 1 ||
        fwrite(hdr, sizeof *hdr, 3, fp) != 3 || fwrite(fr, sizeof *fr, 3, fp) != 3 ||
        (nmeta > 0 && fwrite(meta, sizeof *meta, nmeta, fp) != nmeta)) goto done;
    if (nvalid > 0 &&
        (fwrite(idx, sizeof *idx, nvalid, fp) != nvalid || fwrite(st, 1, nvalid, fp) != nvalid ||
         fwrite(ch, sizeof *ch, nvalid, fp) != nvalid || fwrite(co, sizeof *co, nvalid, fp) != nvalid ||
         fwrite(ps, sizeof *ps, nvalid * 3, fp) != nvalid * 3 || fwrite(nr, sizeof *nr, nvalid * 3, fp) != nvalid * 3))
        goto done;
    rc = 0;
done:
    if (fp != NULL && fclose(fp) != 0) rc = -1;
    free(idx); free(st); free(ch); free(co); free(ps); free(nr);
    return rc;
}

int CloseupGrid_load(const char *path, CloseupGrid *g, int32_t *meta, size_t cap, size_t *nmeta)
{
    char magic[8];
    uint64_t hdr[3];
    double fr[3];
    uint32_t ver = 0, nm = 0, *idx = NULL;
    int32_t mbuf[64];
    uint8_t *st = NULL;
    int32_t *ch = NULL, *co = NULL;
    float *ps = NULL, *nr = NULL;
    size_t k = 0, nvalid = 0, nn = 0;
    FILE *fp = NULL;
    int rc = -1;
    memset(g, 0, sizeof *g);
    if (nmeta != NULL) *nmeta = 0;
    fp = fopen(path, "rb");
    if (fp == NULL) return -1;
    if (fread(magic, 1, 8, fp) != 8 || memcmp(magic, CS_GRID_MAGIC, 8) != 0 || fread(&ver, sizeof ver, 1, fp) != 1 ||
        ver != 1 || fread(&nm, sizeof nm, 1, fp) != 1 || nm > 64 || fread(hdr, sizeof *hdr, 3, fp) != 3 ||
        fread(fr, sizeof *fr, 3, fp) != 3 || (nm > 0 && fread(mbuf, sizeof *mbuf, nm, fp) != nm)) goto done;
    if (hdr[0] == 0 || hdr[1] == 0 || hdr[0] > ((uint64_t)1 << 22) || hdr[1] > ((uint64_t)1 << 22) ||
        !(fr[2] > 0.0) || !isfinite(fr[0]) || !isfinite(fr[1])) goto done;
    g->ni = (size_t)hdr[0]; g->nj = (size_t)hdr[1]; g->u0 = fr[0]; g->v0 = fr[1]; g->h = fr[2];
    nn = g->ni * g->nj;
    nvalid = (size_t)hdr[2];
    if (nvalid > nn || cs_grid_arrays(g) != 0) goto done;
    idx = (uint32_t *)malloc((nvalid ? nvalid : 1) * sizeof *idx);
    st = (uint8_t *)malloc(nvalid ? nvalid : 1);
    ch = (int32_t *)malloc((nvalid ? nvalid : 1) * sizeof *ch);
    co = (int32_t *)malloc((nvalid ? nvalid : 1) * sizeof *co);
    ps = (float *)malloc((nvalid ? nvalid : 1) * 3 * sizeof *ps);
    nr = (float *)malloc((nvalid ? nvalid : 1) * 3 * sizeof *nr);
    if (idx == NULL || st == NULL || ch == NULL || co == NULL || ps == NULL || nr == NULL) goto done;
    if (nvalid > 0 &&
        (fread(idx, sizeof *idx, nvalid, fp) != nvalid || fread(st, 1, nvalid, fp) != nvalid ||
         fread(ch, sizeof *ch, nvalid, fp) != nvalid || fread(co, sizeof *co, nvalid, fp) != nvalid ||
         fread(ps, sizeof *ps, nvalid * 3, fp) != nvalid * 3 || fread(nr, sizeof *nr, nvalid * 3, fp) != nvalid * 3))
        goto done;
    for (k = 0; k < nvalid; k++) {
        const size_t n = (size_t)idx[k];
        if (n >= nn || st[k] == CLOSEUP_NODE_NONE) goto done;
        g->state[n] = st[k]; g->chart[n] = ch[k]; g->component[n] = co[k];
        memcpy(&g->pos[n * 3], &ps[k * 3], 3 * sizeof *ps);
        memcpy(&g->nrm[n * 3], &nr[k * 3], 3 * sizeof *nr);
    }
    if (nmeta != NULL) *nmeta = nm;
    for (k = 0; k < nm && k < cap && meta != NULL; k++) meta[k] = mbuf[k];
    rc = 0;
done:
    fclose(fp);
    free(idx); free(st); free(ch); free(co); free(ps); free(nr);
    if (rc != 0) CloseupGrid_free(g);
    return rc;
}

int CloseupGrid_taylor(const CloseupGrid *g, size_t i, size_t j, double u, double v, double p[3])
{
    size_t k = 0, n = 0;
    double du = 0.0, dv = 0.0;
    if (i >= g->ni || j >= g->nj) return -1;
    n = j * g->ni + i;
    if (g->state[n] == CLOSEUP_NODE_NONE) return -1;
    du = u - (g->u0 + g->h * (double)i);
    dv = v - (g->v0 + g->h * (double)j);
    for (k = 0; k < 3; k++)
        p[k] = (double)g->pos[n * 3 + k] + (double)g->tu[n * 3 + k] * du + (double)g->tv[n * 3 + k] * dv;
    return 0;
}

int CloseupGrid_eval(const CloseupGrid *g, double u, double v, double p[3], double n[3])
{
    double fx = (u - g->u0) / g->h, fy = (v - g->v0) / g->h, ax = 0.0, ay = 0.0, l = 0.0;
    long i0 = (long)floor(fx), j0 = (long)floor(fy);
    size_t c[4];
    int k = 0, q = 0;
    if (i0 < 0 || j0 < 0 || (size_t)i0 + 1 >= g->ni || (size_t)j0 + 1 >= g->nj) return -1;
    ax = fx - (double)i0; ay = fy - (double)j0;
    c[0] = (size_t)j0 * g->ni + (size_t)i0; c[1] = c[0] + 1; c[2] = c[0] + g->ni; c[3] = c[2] + 1;
    for (q = 0; q < 4; q++) if (g->state[c[q]] == CLOSEUP_NODE_NONE) return -1;
    for (k = 0; k < 3; k++) {
        double w0 = (1.0 - ax) * (1.0 - ay), w1 = ax * (1.0 - ay), w2 = (1.0 - ax) * ay, w3 = ax * ay;
        p[k] = w0 * g->pos[c[0] * 3 + (size_t)k] + w1 * g->pos[c[1] * 3 + (size_t)k] +
               w2 * g->pos[c[2] * 3 + (size_t)k] + w3 * g->pos[c[3] * 3 + (size_t)k];
        if (n != NULL)
            n[k] = w0 * g->nrm[c[0] * 3 + (size_t)k] + w1 * g->nrm[c[1] * 3 + (size_t)k] +
                   w2 * g->nrm[c[2] * 3 + (size_t)k] + w3 * g->nrm[c[3] * 3 + (size_t)k];
    }
    if (n != NULL) {
        l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        if (l > 1e-12) for (k = 0; k < 3; k++) n[k] /= l;
    }
    return 0;
}

void CloseupFitOpts_default(CloseupFitOpts *o)
{
    o->fit_px = CLOSEUP_FIT_PX;
    o->min_obs = CLOSEUP_FIT_MIN_OBS;
    o->data_reach_px = CLOSEUP_DATA_REACH_PX;
    o->fill_px = CLOSEUP_FILL_PX;
    o->fill_fit_min_px = CLOSEUP_FILL_FIT_MIN_PX;
    o->fill_fit_max_px = CLOSEUP_FILL_FIT_MAX_PX;
    o->robust_vox = CLOSEUP_ROBUST_VOX;
    o->robust_rounds = CLOSEUP_ROBUST_ROUNDS;
    o->tear_ratio = CLOSEUP_TEAR_RATIO;
    o->tear_float_vox = CLOSEUP_TEAR_FLOAT_VOX;
    o->chart_wrap = NULL;
    o->nchart_wrap = 0;
    o->bridge = 0;
    o->bridge_px = CLOSEUP_BRIDGE_PX;
    o->bridge_agree = CLOSEUP_BRIDGE_AGREE;
    o->variational_holes = 0;
}

void CloseupFitOpts_complete(CloseupFitOpts *o)
{
    CloseupFitOpts_default(o);
    o->fill_px = CLOSEUP_COMPLETE_FILL_PX;
    o->fill_fit_max_px = CLOSEUP_COMPLETE_FILL_FIT_MAX_PX;
    o->variational_holes = 1;
}

/* ------------------------------------------------------ distance transform */

/* 1-D squared distance transform (Felzenszwalb & Huttenlocher) of f over n samples with stride `st`,
 * carrying the argument of the minimum.  f == CS_INF abstains.  v/z scratch of n + 2.  The parabola
 * intersections are computed in double: q^2 passes 2^24 at q = 4,096 nodes. */
static void cs_edt1d(float *f, int32_t *arg, size_t n, size_t st, int32_t *v, double *z, float *df, int32_t *da)
{
    long k = -1;
    size_t q = 0;
    for (q = 0; q < n; q++) {
        double fq = (double)f[q * st];
        if (f[q * st] >= CS_INF) continue;
        if (k < 0) { k = 0; v[0] = (int32_t)q; z[0] = -1e300; z[1] = 1e300; continue; }
        for (;;) {
            double fv = (double)f[(size_t)v[k] * st], vk = (double)v[k], qd = (double)q;
            double s = ((fq + qd * qd) - (fv + vk * vk)) / (2.0 * qd - 2.0 * vk);
            if (s <= z[k]) { k--; if (k < 0) break; continue; }
            k++;
            v[k] = (int32_t)q; z[k] = s; z[k + 1] = 1e300;
            break;
        }
        if (k < 0) { k = 0; v[0] = (int32_t)q; z[0] = -1e300; z[1] = 1e300; }
    }
    if (k < 0) {
        for (q = 0; q < n; q++) { df[q] = CS_INF; da[q] = -1; }
    } else {
        k = 0;
        for (q = 0; q < n; q++) {
            double dq = 0.0;
            while (z[k + 1] < (double)q) k++;
            dq = (double)q - (double)v[k];
            df[q] = (float)(dq * dq + (double)f[(size_t)v[k] * st]);
            da[q] = arg[(size_t)v[k] * st];
        }
    }
    for (q = 0; q < n; q++) { f[q * st] = df[q]; arg[q * st] = da[q]; }
}

/* 2-D squared EDT over an ni x nj grid: seed[k] != 0 are the sources.  d2 in node units^2, near = source. */
static int cs_edt2d(const uint8_t *seed, size_t ni, size_t nj, float *d2, int32_t *near_node)
{
    size_t nn = ni * nj, k = 0, mx = ni > nj ? ni : nj;
    long r = 0;
    int bad = 0;
    for (k = 0; k < nn; k++) {
        d2[k] = seed[k] ? 0.0f : CS_INF;
        near_node[k] = seed[k] ? (int32_t)k : -1;
    }
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        int32_t *v = (int32_t *)malloc((mx + 2) * sizeof *v), *da = (int32_t *)malloc((mx + 1) * sizeof *da);
        double *z = (double *)malloc((mx + 2) * sizeof *z);
        float *df = (float *)malloc((mx + 1) * sizeof *df);
        if (v == NULL || da == NULL || z == NULL || df == NULL) {
#ifdef _OPENMP
#pragma omp atomic
#endif
            bad++;
        }
#ifdef _OPENMP
#pragma omp barrier
#endif
        if (!bad) {
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (r = 0; r < (long)nj; r++)
                cs_edt1d(d2 + (size_t)r * ni, near_node + (size_t)r * ni, ni, 1, v, z, df, da);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (r = 0; r < (long)ni; r++)
                cs_edt1d(d2 + (size_t)r, near_node + (size_t)r, nj, ni, v, z, df, da);
        }
        free(v); free(da); free(z); free(df);
    }
    return bad ? -1 : 0;
}

/* ------------------------------------------------------------- the fit */

/* 6x6 SPD solve for three right-hand sides (Cholesky); 0 on success. */
static int cs_chol6(double A[36], double B[18])
{
    double L[36];
    int i = 0, j = 0, k = 0, r = 0;
    memset(L, 0, sizeof L);
    for (i = 0; i < 6; i++)
        for (j = 0; j <= i; j++) {
            double s = A[i * 6 + j];
            for (k = 0; k < j; k++) s -= L[i * 6 + k] * L[j * 6 + k];
            if (i == j) {
                if (!(s > 1e-12)) return -1;
                L[i * 6 + i] = sqrt(s);
            } else {
                L[i * 6 + j] = s / L[j * 6 + j];
            }
        }
    for (r = 0; r < 3; r++) {
        double yv[6], xv[6];
        for (i = 0; i < 6; i++) {
            double s = B[r * 6 + i];
            for (k = 0; k < i; k++) s -= L[i * 6 + k] * yv[k];
            yv[i] = s / L[i * 6 + i];
        }
        for (i = 5; i >= 0; i--) {
            double s = yv[i];
            for (k = i + 1; k < 6; k++) s -= L[k * 6 + i] * xv[k];
            xv[i] = s / L[i * 6 + i];
        }
        for (i = 0; i < 6; i++) B[r * 6 + i] = xv[i];
    }
    return 0;
}

typedef struct CsQuad {
    double A[36], B[18], wsum;
    size_t n;
} CsQuad;

static void cs_quad_add(CsQuad *q, double dx, double dy, double inv, double w, const float *p)
{
    double phi[6];
    int j = 0, k = 0;
    phi[0] = 1.0; phi[1] = dx * inv; phi[2] = dy * inv;
    phi[3] = phi[1] * phi[1]; phi[4] = phi[1] * phi[2]; phi[5] = phi[2] * phi[2];
    for (j = 0; j < 6; j++) {
        double wp = w * phi[j];
        for (k = 0; k <= j; k++) q->A[j * 6 + k] += wp * phi[k];
        q->B[j] += wp * (double)p[0];
        q->B[6 + j] += wp * (double)p[1];
        q->B[12 + j] += wp * (double)p[2];
    }
    q->wsum += w;
    q->n++;
}

/* Solve the quadric; point and derivatives (per px) into the node.  0 on success. */
static int cs_quad_solve(CsQuad *q, double inv, float *pos, float *tu, float *tv)
{
    int j = 0, k = 0;
    for (j = 0; j < 6; j++)
        for (k = j + 1; k < 6; k++) q->A[j * 6 + k] = q->A[k * 6 + j];
    /* the relative ridge of the apron fit (rawtex_layers.c): one-sided border fits of thin strips stay well
     * posed without biasing the curvature */
    for (j = 3; j < 6; j++) q->A[j * 6 + j] += 1e-3 * q->A[j * 6 + j] + 1e-9 * q->wsum;
    if (cs_chol6(q->A, q->B) != 0) return -1;
    for (k = 0; k < 3; k++) {
        pos[k] = (float)q->B[k * 6];
        tu[k] = (float)(q->B[k * 6 + 1] * inv);
        tv[k] = (float)(q->B[k * 6 + 2] * inv);
    }
    return 0;
}

static void cs_cross(const float *a, const float *b, double c[3])
{
    c[0] = (double)a[1] * b[2] - (double)a[2] * b[1];
    c[1] = (double)a[2] * b[0] - (double)a[0] * b[2];
    c[2] = (double)a[0] * b[1] - (double)a[1] * b[0];
}

static void cs_set_normal(CloseupGrid *g, size_t n)
{
    double c[3], l = 0.0;
    int k = 0;
    cs_cross(&g->tu[n * 3], &g->tv[n * 3], c);
    l = sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
    for (k = 0; k < 3; k++) g->nrm[n * 3 + (size_t)k] = l > 1e-12 ? (float)(c[k] / l) : 0.0f;
}

/* Observation bins: node cells [c - h/2, c + h/2)^2, CSR. */
typedef struct CsBins {
    size_t *off;      /* [nn + 1] */
    uint32_t *idx;    /* [nobs] */
} CsBins;

static long cs_cell(const CloseupGrid *g, double u, double v)
{
    double fx = (u - g->u0) / g->h + 0.5, fy = (v - g->v0) / g->h + 0.5;
    long i = (long)floor(fx), j = (long)floor(fy);
    if (i < 0 || j < 0 || (size_t)i >= g->ni || (size_t)j >= g->nj) return -1;
    return j * (long)g->ni + i;
}

static int cs_bin(const CloseupGrid *g, const CloseupObs *o, CsBins *b)
{
    size_t nn = g->ni * g->nj, k = 0, *fill = NULL;
    memset(b, 0, sizeof *b);
    if (o->n >= (size_t)UINT32_MAX) return -1;
    b->off = (size_t *)calloc(nn + 1, sizeof *b->off);
    b->idx = (uint32_t *)malloc((o->n ? o->n : 1) * sizeof *b->idx);
    fill = (size_t *)malloc((nn + 1) * sizeof *fill);
    if (b->off == NULL || b->idx == NULL || fill == NULL) {
        free(b->off); free(b->idx); free(fill);
        memset(b, 0, sizeof *b);
        return -1;
    }
    for (k = 0; k < o->n; k++) {
        long c = cs_cell(g, o->uv[k * 2], o->uv[k * 2 + 1]);
        if (c >= 0) b->off[(size_t)c + 1]++;
    }
    for (k = 0; k < nn; k++) b->off[k + 1] += b->off[k];
    memcpy(fill, b->off, (nn + 1) * sizeof *fill);
    for (k = 0; k < o->n; k++) {
        long c = cs_cell(g, o->uv[k * 2], o->uv[k * 2 + 1]);
        if (c >= 0) b->idx[fill[c]++] = (uint32_t)k;
    }
    free(fill);
    return 0;
}

static void cs_bins_free(CsBins *b)
{
    free(b->off); free(b->idx);
    memset(b, 0, sizeof *b);
}

/* One data fit at node n from the binned observations; 0 = fitted. */
static int cs_fit_data_node(CloseupGrid *g, const CloseupObs *o, const CsBins *b, size_t i, size_t j,
                            const CloseupFitOpts *opt, int record_nearest)
{
    const double R = opt->fit_px, inv = 1.0 / R, s2 = 2.0 * (0.5 * R) * (0.5 * R);
    const long rc = (long)ceil(R / g->h) + 1;
    const double cu = g->u0 + g->h * (double)i, cv = g->v0 + g->h * (double)j;
    size_t n = j * g->ni + i;
    long ci = 0, cj = 0;
    double best = 1e300;
    int32_t best_chart = -1;
    CsQuad q;
    memset(&q, 0, sizeof q);
    for (cj = (long)j - rc; cj <= (long)j + rc; cj++) {
        if (cj < 0 || cj >= (long)g->nj) continue;
        for (ci = (long)i - rc; ci <= (long)i + rc; ci++) {
            size_t c = 0, t = 0;
            if (ci < 0 || ci >= (long)g->ni) continue;
            c = (size_t)cj * g->ni + (size_t)ci;
            for (t = b->off[c]; t < b->off[c + 1]; t++) {
                size_t k = b->idx[t];
                double dx = (double)o->uv[k * 2] - cu, dy = (double)o->uv[k * 2 + 1] - cv;
                double r2 = dx * dx + dy * dy, w = 0.0;
                if (r2 > R * R) continue;
                if (r2 < best) { best = r2; best_chart = o->chart[k]; }
                w = exp(-r2 / s2) * (double)o->w[k];
                if (w < 1e-6) continue;
                cs_quad_add(&q, dx, dy, inv, w, &o->p[k * 3]);
            }
        }
    }
    if (record_nearest) {
        g->dist[n] = best < 1e300 ? (float)sqrt(best) : INFINITY;
        g->chart[n] = best_chart;
    }
    if ((long)q.n < opt->min_obs) return -1;
    if (cs_quad_solve(&q, inv, &g->pos[n * 3], &g->tu[n * 3], &g->tv[n * 3]) != 0) return -1;
    g->support[n] = (float)q.wsum;
    return 0;
}

static int cs_cmp_float(const void *a, const void *b)
{
    float x = *(const float *)a, y = *(const float *)b;
    return x < y ? -1 : x > y ? 1 : 0;
}

/* Residuals of every observation against the surface (first-order at its nearest DATA node), Cauchy
 * weights; percentile stats over a subsample.  Observations with no DATA node keep their weight. */
static void cs_reweight(const CloseupGrid *g, CloseupObs *o, double tau, int update, CloseupFitStats *st)
{
    long k = 0;
    size_t nsamp = 0, down = 0, stride = o->n > 2000000 ? o->n / 2000000 : 1, m = 0;
    float *res = (float *)malloc((o->n / stride + 1) * sizeof *res);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+ : down)
#endif
    for (k = 0; k < (long)o->n; k++) {
        double u = o->uv[(size_t)k * 2], v = o->uv[(size_t)k * 2 + 1], p[3], r = 0.0;
        long i = (long)floor((u - g->u0) / g->h + 0.5), j = (long)floor((v - g->v0) / g->h + 0.5);
        int q = 0;
        if (i < 0 || j < 0 || (size_t)i >= g->ni || (size_t)j >= g->nj ||
            g->state[(size_t)j * g->ni + (size_t)i] != CLOSEUP_NODE_DATA) {
            if (res != NULL && (size_t)k % stride == 0) res[(size_t)k / stride] = -1.0f;
            continue;
        }
        (void)CloseupGrid_taylor(g, (size_t)i, (size_t)j, u, v, p);
        for (q = 0; q < 3; q++) { double d = p[q] - (double)o->p[(size_t)k * 3 + (size_t)q]; r += d * d; }
        r = sqrt(r);
        if (update) o->w[k] = (float)(1.0 / (1.0 + (r / tau) * (r / tau)));
        if (o->w[k] < 0.5f) down++;
        if (res != NULL && (size_t)k % stride == 0) res[(size_t)k / stride] = (float)r;
    }
    if (res != NULL) {
        for (m = 0; m <= (o->n - 1) / stride && m < o->n / stride + 1; m++)
            if (res[m] >= 0.0f) res[nsamp++] = res[m];
        if (nsamp > 0) {
            qsort(res, nsamp, sizeof *res, cs_cmp_float);
            st->resid_p50 = res[(size_t)(0.50 * (double)(nsamp - 1))];
            st->resid_p90 = res[(size_t)(0.90 * (double)(nsamp - 1))];
            st->resid_p99 = res[(size_t)(0.99 * (double)(nsamp - 1))];
        }
        free(res);
    }
    st->obs_downweighted = down;
}

/* One fill fit at node n from the DATA nodes around it; gates against the nearest data node. */
static int cs_fit_fill_node(CloseupGrid *g, size_t i, size_t j, double dist_px, size_t near_node,
                            const CloseupFitOpts *opt, const int32_t *comp)
{
    double R = 1.5 * dist_px + opt->fill_fit_min_px, inv = 0.0, s2 = 0.0, d3 = 0.0, dotn = 0.0, lu = 0.0, lv = 0.0;
    const double cu = g->u0 + g->h * (double)i, cv = g->v0 + g->h * (double)j;
    long rn = 0, st = 1, ci = 0, cj = 0;
    size_t n = j * g->ni + i;
    double c[3], l = 0.0;
    int k = 0;
    CsQuad q;
    if (R > opt->fill_fit_max_px) R = opt->fill_fit_max_px;
    if (R < opt->fill_fit_min_px) R = opt->fill_fit_min_px;
    inv = 1.0 / R; s2 = 2.0 * (0.5 * R) * (0.5 * R);
    rn = (long)ceil(R / g->h);
    st = (long)floor(R / (8.0 * g->h));
    if (st < 1) st = 1;
    memset(&q, 0, sizeof q);
    /* the stride lattice is anchored at the global origin, so neighbouring fill nodes see the same samples */
    for (cj = (long)j - rn; cj <= (long)j + rn; cj++) {
        if (cj < 0 || cj >= (long)g->nj || cj % st) continue;
        for (ci = (long)i - rn; ci <= (long)i + rn; ci++) {
            size_t m = 0;
            double dx = 0.0, dy = 0.0, r2 = 0.0;
            if (ci < 0 || ci >= (long)g->ni || ci % st) continue;
            m = (size_t)cj * g->ni + (size_t)ci;
            /* comp is the immutable DATA mask. Other workers publish FILL
             * into state[] during this loop; reading those bytes is a race. */
            if (comp != NULL ? comp[m] < 0 || comp[m] != comp[near_node] : g->state[m] != CLOSEUP_NODE_DATA) continue;
            dx = g->h * (double)(ci - (long)i); dy = g->h * (double)(cj - (long)j);
            r2 = dx * dx + dy * dy;
            if (r2 > R * R) continue;
            cs_quad_add(&q, dx, dy, inv, exp(-r2 / s2), &g->pos[m * 3]);
        }
    }
    (void)cu; (void)cv;
    if (q.n < 12) return -1;
    if (cs_quad_solve(&q, inv, &g->pos[n * 3], &g->tu[n * 3], &g->tv[n * 3]) != 0) return -1;
    /* the apron's isometry gate: tangent lengths 0.5..2 vox/px, within 1.6x the pixel distance of the nearest
     * data node, normal within 60 deg of it */
    for (k = 0; k < 3; k++) {
        double d = (double)g->pos[n * 3 + (size_t)k] - (double)g->pos[near_node * 3 + (size_t)k];
        d3 += d * d;
        lu += (double)g->tu[n * 3 + (size_t)k] * g->tu[n * 3 + (size_t)k];
        lv += (double)g->tv[n * 3 + (size_t)k] * g->tv[n * 3 + (size_t)k];
    }
    d3 = sqrt(d3); lu = sqrt(lu); lv = sqrt(lv);
    if (lu < 0.5 || lu > 2.0 || lv < 0.5 || lv > 2.0 || d3 > 1.6 * dist_px + 2.0) return -1;
    cs_cross(&g->tu[n * 3], &g->tv[n * 3], c);
    l = sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
    if (!(l > 1e-9)) return -1;
    {
        double cn[3];
        double ln = 0.0;
        cs_cross(&g->tu[near_node * 3], &g->tv[near_node * 3], cn);
        ln = sqrt(cn[0] * cn[0] + cn[1] * cn[1] + cn[2] * cn[2]);
        if (!(ln > 1e-9)) return -1;
        for (k = 0; k < 3; k++) dotn += c[k] * cn[k];
        if (dotn / (l * ln) < 0.5) return -1;
    }
    g->support[n] = (float)q.wsum;
    return 0;
}

static int32_t cs_root(int32_t *par, int32_t a)
{
    while (par[a] != a) { par[a] = par[par[a]]; a = par[a]; }
    return a;
}

/* 8-connected components of on[] over the ni x nj grid: comp[k] = the root index of k's component (-1 off). */
static void cs_components(const uint8_t *on, size_t ni, size_t nj, int32_t *comp)
{
    size_t k = 0, nn = ni * nj;
    for (k = 0; k < nn; k++) comp[k] = on[k] ? (int32_t)k : -1;
    for (k = 0; k < nn; k++) {
        size_t i = k % ni, j = k / ni;
        int q = 0;
        const long off[4][2] = { { -1, 0 }, { -1, -1 }, { 0, -1 }, { 1, -1 } };
        if (!on[k]) continue;
        for (q = 0; q < 4; q++) {
            long ii = (long)i + off[q][0], jj = (long)j + off[q][1];
            int32_t a = 0, b = 0;
            if (ii < 0 || jj < 0 || (size_t)ii >= ni) continue;
            if (!on[(size_t)jj * ni + (size_t)ii]) continue;
            a = cs_root(comp, (int32_t)k); b = cs_root(comp, (int32_t)((size_t)jj * ni + (size_t)ii));
            if (a != b) comp[a] = b;
        }
    }
    for (k = 0; k < nn; k++) if (comp[k] >= 0) comp[k] = cs_root(comp, (int32_t)k);
}

/* ---------------------------------------------------------- gated bridging */

/* unit normal of a fitted node from its tangents; -1 when degenerate */
static int cs_node_normal(const CloseupGrid *g, size_t k, double n[3])
{
    const float *tu = &g->tu[k * 3], *tv = &g->tv[k * 3];
    double l = 0.0;
    n[0] = (double)tu[1] * tv[2] - (double)tu[2] * tv[1];
    n[1] = (double)tu[2] * tv[0] - (double)tu[0] * tv[2];
    n[2] = (double)tu[0] * tv[1] - (double)tu[1] * tv[0];
    l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    if (!(l > 1e-12)) return -1;
    n[0] /= l; n[1] /= l; n[2] /= l;
    return 0;
}

/* The continuity gate of one candidate bridge between data node ka (blob A) and its nearest data node kb of blob B,
 * frame_px apart in the frame: no farther apart in 3-D than CLOSEUP_BRIDGE_STRETCH x that (+ 2 vox), and each within
 * CLOSEUP_BRIDGE_NORMAL_VOX of the other's tangent plane (two layers never bridge, whatever the index says). */
static int cs_bridge_continuous(const CloseupGrid *g, size_t ka, size_t kb, double frame_px, double *d3_out,
                                double *off_out, const char **why)
{
    const float *pa = &g->pos[ka * 3], *pb = &g->pos[kb * 3];
    double tol = CLOSEUP_BRIDGE_NORMAL_VOX > CLOSEUP_BRIDGE_NORMAL_SLOPE * frame_px ? CLOSEUP_BRIDGE_NORMAL_VOX
                                                                                 : CLOSEUP_BRIDGE_NORMAL_SLOPE * frame_px;
    if (tol > CLOSEUP_BRIDGE_NORMAL_MAX_VOX) tol = CLOSEUP_BRIDGE_NORMAL_MAX_VOX;
    double d[3], d3 = 0.0, na[3], nb[3], oa = 0.0, ob = 0.0;
    int q = 0;
    for (q = 0; q < 3; q++) { d[q] = (double)pb[q] - pa[q]; d3 += d[q] * d[q]; }
    d3 = sqrt(d3);
    *d3_out = d3; *off_out = 0.0; *why = "ok";
    if (!(d3 <= CLOSEUP_BRIDGE_STRETCH * frame_px + 2.0)) { *why = "stretch"; return 0; }
    if (cs_node_normal(g, ka, na) != 0 || cs_node_normal(g, kb, nb) != 0) { *why = "degenerate_normal"; return 0; }
    oa = fabs(na[0] * d[0] + na[1] * d[1] + na[2] * d[2]);
    ob = fabs(nb[0] * d[0] + nb[1] * d[1] + nb[2] * d[2]);
    *off_out = oa > ob ? oa : ob;
    if (oa > tol || ob > tol) { *why = "normal_offset"; return 0; }
    return 1;
}

static int cs_cmp_u64(const void *x, const void *y)
{
    const uint64_t a = *(const uint64_t *)x, b = *(const uint64_t *)y;
    return a < b ? -1 : a > b;
}

/* Gated bridging (CloseupFitOpts.bridge): relabel comp[] (the roots of cs_components over the data mask seed[]) so
 * that observed blobs on one wrap by the wrap index, within 2 x bridge_px of each other in the frame and continuing
 * each other at the gap, share one root.  The fill then fits a joint quadric across the gap and the mesh may
 * connect them; the tear test still applies.  Measured 2026-10-13 on PHerc1218's 6^3 rung: 50 of 75 consecutive
 * same-wrap pieces of the wrap-row ribbon lie within 48 px of each other. */
static void cs_bridge_groups(const CloseupGrid *g, int32_t *comp, const uint8_t *seed, const CloseupFitOpts *opt,
                             CloseupFitStats *st)
{
    const size_t nn = g->ni * g->nj;
    const double lim = 2.0 * opt->bridge_px / g->h, lim2 = lim * lim;
    size_t k = 0, nblob = 0, cap = 256, nvote = 0, nwrap = 0, nbridge = 0;
    int32_t *bidx = NULL, *broot = NULL, *bwrap = NULL, *bpar = NULL, *best_k = NULL, *nearb = NULL;
    uint64_t *votes = NULL;
    float *d2b = NULL, *best_d2 = NULL;
    uint8_t *seedb = NULL;
    if (!opt->bridge || opt->chart_wrap == NULL || opt->nchart_wrap <= 0 || !(opt->bridge_px > 0.0)) return;
    bidx = (int32_t *)malloc(nn * sizeof *bidx);
    broot = (int32_t *)malloc(cap * sizeof *broot);
    if (bidx == NULL || broot == NULL) goto done;
    for (k = 0; k < nn; k++) bidx[k] = -1;
    for (k = 0; k < nn; k++) {
        int32_t r = 0;
        if (!seed[k]) continue;
        r = comp[k];
        if (r < 0 || bidx[r] >= 0) continue;
        if (nblob == cap) {
            int32_t *t = (int32_t *)realloc(broot, 2 * cap * sizeof *t);
            if (t == NULL) goto done;
            broot = t; cap *= 2;
        }
        bidx[r] = (int32_t)nblob; broot[nblob++] = r;
    }
    if (nblob < 2) goto done;
    /* each blob's wrap: the vote of its data nodes' charts, decided at bridge_agree */
    votes = (uint64_t *)malloc(nn * sizeof *votes);
    bwrap = (int32_t *)malloc(nblob * sizeof *bwrap);
    bpar = (int32_t *)malloc(nblob * sizeof *bpar);
    best_d2 = (float *)malloc(nblob * sizeof *best_d2);
    best_k = (int32_t *)malloc(nblob * sizeof *best_k);
    if (votes == NULL || bwrap == NULL || bpar == NULL || best_d2 == NULL || best_k == NULL) goto done;
    for (k = 0; k < nn; k++) {
        int32_t c = 0, w = 0;
        if (!seed[k] || comp[k] < 0) continue;
        c = g->chart[k];
        if (c < 0 || c >= opt->nchart_wrap) continue;
        w = opt->chart_wrap[c];
        if (w == INT32_MIN) continue;
        votes[nvote++] = ((uint64_t)(uint32_t)bidx[comp[k]] << 32) | (uint32_t)(w + 1073741824);
    }
    qsort(votes, nvote, sizeof *votes, cs_cmp_u64);
    for (k = 0; k < nblob; k++) { bwrap[k] = INT32_MIN; bpar[k] = (int32_t)k; }
    k = 0;
    while (k < nvote) {
        size_t e = k, q = 0, total = 0, best = 0;
        int32_t bestw = INT32_MIN;
        const uint32_t blob = (uint32_t)(votes[k] >> 32);
        while (e < nvote && (uint32_t)(votes[e] >> 32) == blob) e++;
        total = e - k;
        for (q = k; q < e;) {
            size_t r = q;
            while (r < e && votes[r] == votes[q]) r++;
            if (r - q > best) { best = r - q; bestw = (int32_t)((uint32_t)(votes[q] & 0xffffffffu)) - 1073741824; }
            q = r;
        }
        if (total > 0 && (double)best >= opt->bridge_agree * (double)total) { bwrap[blob] = bestw; nwrap++; }
        k = e;
    }
    /* candidate pairs: the distance field of each decided blob B; the nearest data node of every other blob A of
     * the same wrap within 2 x bridge_px is the continuity gate's sample */
    d2b = (float *)malloc(nn * sizeof *d2b);
    nearb = (int32_t *)malloc(nn * sizeof *nearb);
    seedb = (uint8_t *)malloc(nn);
    if (d2b == NULL || nearb == NULL || seedb == NULL) goto done;
    for (k = 0; k < nblob; k++) {
        size_t m = 0, a = 0;
        if (bwrap[k] == INT32_MIN) continue;
        for (m = 0; m < nn; m++) seedb[m] = (uint8_t)(seed[m] && comp[m] == broot[k]);
        if (cs_edt2d(seedb, g->ni, g->nj, d2b, nearb) != 0) goto done;
        for (a = 0; a < nblob; a++) { best_d2[a] = CS_INF; best_k[a] = -1; }
        for (m = 0; m < nn; m++) {           /* the nearest data node of every other same-wrap blob (any distance) */
            int32_t ab = 0;
            if (!seed[m] || comp[m] < 0 || comp[m] == broot[k]) continue;
            ab = bidx[comp[m]];
            if (ab < 0 || bwrap[ab] != bwrap[k] || d2b[m] >= best_d2[ab]) continue;
            best_d2[ab] = d2b[m]; best_k[ab] = (int32_t)m;
        }
        for (a = 0; a < nblob; a++) {
            int32_t ra = 0, rb = 0;
            double frame_px = 0.0, d3 = 0.0, off = 0.0;
            const char *why = "";
            size_t kb = 0;
            if (a < k || best_k[a] < 0 || nearb[best_k[a]] < 0 || best_d2[a] >= CS_INF) continue;   /* each pair once */
            frame_px = sqrt((double)best_d2[a]) * g->h;
            if (frame_px > 2.0 * opt->bridge_px + 100.0) continue;
            kb = (size_t)nearb[best_k[a]];
            ra = cs_root(bpar, (int32_t)a); rb = cs_root(bpar, (int32_t)k);
            if (best_d2[a] > lim2) why = "beyond_reach";
            else if (ra == rb) why = "already_joined";
            else if (cs_bridge_continuous(g, (size_t)best_k[a], kb, frame_px, &d3, &off, &why)) { bpar[ra] = rb; nbridge++; why = "BRIDGED"; }
            fprintf(stderr, "    bridge candidate: blobs %zu + %zu (wrap %d) frame gap %.1f px, 3-D %.1f vox, normal offset %.1f vox at frame (%.0f, %.0f) xyz (%.0f, %.0f, %.0f): %s\n",
                    a, k, bwrap[k], frame_px, d3, off, g->u0 + g->h * (double)((size_t)best_k[a] % g->ni), g->v0 + g->h * (double)((size_t)best_k[a] / g->ni),
                    (double)g->pos[(size_t)best_k[a] * 3], (double)g->pos[(size_t)best_k[a] * 3 + 1], (double)g->pos[(size_t)best_k[a] * 3 + 2], why);
        }
    }
    /* every node's component becomes its blob's group root */
    for (k = 0; k < nn; k++) {
        int32_t bl = 0;
        if (comp[k] < 0) continue;
        bl = bidx[comp[k]];
        if (bl < 0) continue;
        comp[k] = broot[cs_root(bpar, bl)];
    }
done:
    if (st != NULL) { st->blobs = nblob; st->blobs_with_wrap = nwrap; st->bridges = nbridge; }
    free(bidx); free(broot); free(bwrap); free(bpar); free(best_d2); free(best_k); free(votes); free(d2b); free(nearb); free(seedb);
}

/* The tear (CLOSEUP_TEAR_*): nodes whose fit cannot be one sheet go back to NONE -- a DATA node more than
 * tear_float_vox in 3-D from every observation within data_reach_px of it, and both ends of every 4-neighbour
 * grid edge longer than tear_ratio x its frame length.  Removed outright (the caller's final pass clears their
 * geometry), never moved. */
static int cs_tear(CloseupGrid *g, const CloseupObs *o, const CsBins *b, const CloseupFitOpts *opt,
                    CloseupFitStats *st)
{
    const size_t nn = g->ni * g->nj;
    uint8_t *cut = NULL;
    size_t k = 0, nfloat = 0, nstretch = 0;
    long n = 0;
    if (!(opt->tear_ratio > 0.0) && !(opt->tear_float_vox > 0.0)) return 0;
    cut = (uint8_t *)calloc(nn ? nn : 1, 1);
    if (cut == NULL) return -1;                    /* an unmeasured surface cannot pass */
    if (opt->tear_float_vox > 0.0) {
        const double reach2 = opt->data_reach_px * opt->data_reach_px, fl2 = opt->tear_float_vox * opt->tear_float_vox;
        const long rc = (long)ceil(opt->data_reach_px / g->h) + 1;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 256)
#endif
        for (n = 0; n < (long)nn; n++) {
            const size_t i = (size_t)n % g->ni, j = (size_t)n / g->ni;
            const double cu = g->u0 + g->h * (double)i, cv = g->v0 + g->h * (double)j;
            double best = 1e300;
            long ci = 0, cj = 0;
            if (g->state[n] != CLOSEUP_NODE_DATA) continue;
            for (cj = (long)j - rc; cj <= (long)j + rc && best > fl2; cj++) {
                if (cj < 0 || cj >= (long)g->nj) continue;
                for (ci = (long)i - rc; ci <= (long)i + rc && best > fl2; ci++) {
                    size_t c = 0, t = 0;
                    if (ci < 0 || ci >= (long)g->ni) continue;
                    c = (size_t)cj * g->ni + (size_t)ci;
                    for (t = b->off[c]; t < b->off[c + 1]; t++) {
                        const size_t kk = b->idx[t];
                        const double du = (double)o->uv[kk * 2] - cu, dv = (double)o->uv[kk * 2 + 1] - cv;
                        double d2 = 0.0;
                        int q = 0;
                        if (du * du + dv * dv > reach2) continue;
                        for (q = 0; q < 3; q++) {
                            const double e = (double)o->p[kk * 3 + (size_t)q] - g->pos[(size_t)n * 3 + (size_t)q];
                            d2 += e * e;
                        }
                        if (d2 < best) best = d2;
                    }
                }
            }
            if (best > fl2) cut[n] = 1;
        }
    }
    if (opt->tear_ratio > 0.0) {
        const double lim2 = (opt->tear_ratio * g->h) * (opt->tear_ratio * g->h);
        size_t i = 0, j = 0;
        for (j = 0; j < g->nj; j++)
            for (i = 0; i < g->ni; i++) {
                const size_t a = j * g->ni + i;
                int dir = 0;
                if (g->state[a] == CLOSEUP_NODE_NONE) continue;
                for (dir = 0; dir < 2; dir++) {
                    size_t nb = 0;
                    double d2 = 0.0;
                    int q = 0;
                    if (dir == 0 ? i + 1 >= g->ni : j + 1 >= g->nj) continue;
                    nb = dir == 0 ? a + 1 : a + g->ni;
                    if (g->state[nb] == CLOSEUP_NODE_NONE) continue;
                    for (q = 0; q < 3; q++) {
                        const double e = (double)g->pos[nb * 3 + (size_t)q] - g->pos[a * 3 + (size_t)q];
                        d2 += e * e;
                    }
                    if (d2 > lim2) { cut[a] |= 2; cut[nb] |= 2; }
                }
            }
    }
    for (k = 0; k < nn; k++) {
        if (!cut[k]) continue;
        if (cut[k] & 2) nstretch++; else nfloat++;
        g->state[k] = CLOSEUP_NODE_NONE;
        g->chart[k] = -1;
    }
    free(cut);
    st->torn_stretch = nstretch;
    st->torn_float = nfloat;
    return 0;
}

#include "closeup_infill.inc"

int CloseupSurface_fit(CloseupGrid *g, CloseupObs *obs, const CloseupFitOpts *opt_in, CloseupFitStats *st)
{
    CloseupFitOpts opt;
    CloseupFitStats s0;
    CsBins b;
    size_t nn = g->ni * g->nj, k = 0;
    uint8_t *seed = NULL;
    float *d2 = NULL;
    int32_t *near_node = NULL, *comp = NULL;
    double t0 = ves_clock_sec();
    long n = 0;
    int round = 0, rc = -1;
    if (opt_in == NULL) CloseupFitOpts_default(&opt); else opt = *opt_in;
    if (st == NULL) st = &s0;
    memset(st, 0, sizeof *st);
    st->nodes = nn;
    st->obs = obs->n;
    if (obs->n == 0 || nn == 0) return -1;
    if (cs_bin(g, obs, &b) != 0) return -1;
    seed = (uint8_t *)calloc(nn, 1);
    d2 = (float *)malloc(nn * sizeof *d2);
    near_node = (int32_t *)malloc(nn * sizeof *near_node);
    if (seed == NULL || d2 == NULL || near_node == NULL) goto done;
    /* candidates: nodes with an observation in their own cell or a neighbour's (exact distance below) */
    for (k = 0; k < nn; k++) seed[k] = b.off[k + 1] > b.off[k];
    if (cs_edt2d(seed, g->ni, g->nj, d2, near_node) != 0) goto done;
    for (k = 0; k < nn; k++) {
        g->state[k] = CLOSEUP_NODE_NONE;
        g->dist[k] = d2[k] < CS_INF ? (float)(sqrt((double)d2[k]) * g->h) : INFINITY;
    }
    for (k = 0; k < obs->n; k++) obs->w[k] = 1.0f;
    for (round = 0; round <= opt.robust_rounds; round++) {
        size_t nd = 0, nfail = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 256) reduction(+ : nd, nfail)
#endif
        for (n = 0; n < (long)nn; n++) {
            size_t i = (size_t)n % g->ni, j = (size_t)n / g->ni;
            if (round == 0) {
                if (d2[n] > 2.0f) continue;    /* no observation within ~1.5 cells */
                if (cs_fit_data_node(g, obs, &b, i, j, &opt, 1) != 0 || g->dist[n] > (float)opt.data_reach_px) {
                    if (g->dist[n] <= (float)opt.data_reach_px) nfail++;
                    continue;
                }
                g->state[n] = CLOSEUP_NODE_DATA;
                nd++;
            } else {
                if (g->state[n] != CLOSEUP_NODE_DATA) continue;
                if (cs_fit_data_node(g, obs, &b, i, j, &opt, 0) != 0) {
                    g->state[n] = CLOSEUP_NODE_NONE;
                    nfail++;
                    continue;
                }
                nd++;
            }
        }
        st->data_nodes = nd;
        st->data_failed += nfail;
        cs_reweight(g, obs, opt.robust_vox, round < opt.robust_rounds, st);
    }
    /* fill: distance to the nearest DATA node; quadric of the data nodes of its component around */
    for (k = 0; k < nn; k++) seed[k] = g->state[k] == CLOSEUP_NODE_DATA;
    if (cs_edt2d(seed, g->ni, g->nj, d2, near_node) != 0) goto done;
    comp = (int32_t *)malloc(nn * sizeof *comp);
    if (comp == NULL) goto done;
    cs_components(seed, g->ni, g->nj, comp);
    cs_bridge_groups(g, comp, seed, &opt, st);
    {
        size_t nfill = 0, nrej = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 256) reduction(+ : nfill, nrej)
#endif
        for (n = 0; n < (long)nn; n++) {
            size_t i = (size_t)n % g->ni, j = (size_t)n / g->ni;
            double dpx = 0.0;
            if (g->state[n] == CLOSEUP_NODE_DATA || d2[n] >= CS_INF || near_node[n] < 0) continue;
            dpx = sqrt((double)d2[n]) * g->h;
            g->dist[n] = (float)dpx;
            if (dpx > opt.fill_px) continue;
            if (cs_fit_fill_node(g, i, j, dpx, (size_t)near_node[n], &opt, comp) != 0) {
                size_t kk = 0;
                for (kk = 0; kk < 3; kk++) g->pos[(size_t)n * 3 + kk] = NAN;
                nrej++;
                continue;
            }
            g->state[n] = CLOSEUP_NODE_FILL;
            g->chart[n] = -1;
            nfill++;
        }
        st->fill_nodes = nfill;
        st->fill_rejected = nrej;
    }
    if (cs_tear(g, obs, &b, &opt, st) != 0) goto done;
    for (k = 0; k < nn; k++) {
        if (g->state[k] == CLOSEUP_NODE_NONE) {
            size_t q = 0;
            g->component[k] = -1;
            for (q = 0; q < 3; q++) g->pos[k * 3 + q] = NAN;
            continue;
        }
        g->component[k] = g->state[k] == CLOSEUP_NODE_DATA ? comp[k] : comp[near_node[k]];
        cs_set_normal(g, k);
    }
    if (opt.variational_holes && CloseupSurface_fill_holes(g, opt.fill_px, st) != 0) goto done;
    rc = 0;
done:
    cs_bins_free(&b);
    free(seed); free(d2); free(near_node); free(comp);
    st->seconds = ves_clock_sec() - t0;
    return rc;
}

int CloseupSurface_orient(CloseupGrid *g, const float *inward, double *agree)
{
    size_t nn = g->ni * g->nj, k = 0;
    double pos = 0.0, neg = 0.0;
    int s = 0;
    for (k = 0; k < nn; k++) {
        double c[3], vote = 0.0;
        if (g->state[k] == CLOSEUP_NODE_NONE) continue;
        cs_cross(&g->tu[k * 3], &g->tv[k * 3], c);
        vote = c[0] * inward[k * 3] + c[1] * inward[k * 3 + 1] + c[2] * inward[k * 3 + 2];
        if (vote > 0.0) pos += vote; else neg -= vote;
    }
    if (pos == 0.0 && neg == 0.0) { if (agree) *agree = 0.0; return 0; }
    s = pos >= neg ? 1 : -1;
    if (agree) *agree = (s > 0 ? pos : neg) / (pos + neg);
    for (k = 0; k < nn; k++) {
        size_t q = 0;
        if (g->state[k] == CLOSEUP_NODE_NONE) continue;
        cs_set_normal(g, k);
        for (q = 0; q < 3; q++) g->nrm[k * 3 + q] *= (float)s;
    }
    return s;
}

/* ------------------------------------------------------------ flattening */

double CloseupSurface_rot_angle(double m00, double m01, double m10, double m11)
{
    return atan2(m10 - m01, m00 + m11);
}

typedef struct CsTri {
    int32_t v[3];
    double q[6];        /* rest 2-D corners (x0,y0, x1,y1, x2,y2) */
    double cot[3];      /* cot of the angle at corner k (opposite edge k+1 -> k+2) */
    double area;
} CsTri;

static int cs_dcmp(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return x < y ? -1 : x > y ? 1 : 0;
}

/* |log sigma| of the map rest -> uv over each triangle; percentiles 50/99 into p50/p99. */
static void cs_distortion(const CsTri *T, size_t nt, const double *uv, double *p50, double *p99, double *area_uv,
                          double *area_3d)
{
    double *v = (double *)malloc((nt * 2 + 1) * sizeof *v);
    size_t t = 0, m = 0;
    double au = 0.0, a3 = 0.0;
    for (t = 0; t < nt; t++) {
        const CsTri *r = &T[t];
        double x1 = r->q[2] - r->q[0], y1 = r->q[3] - r->q[1], x2 = r->q[4] - r->q[0], y2 = r->q[5] - r->q[1];
        double det = x1 * y2 - x2 * y1, u1 = 0.0, v1 = 0.0, u2 = 0.0, v2 = 0.0;
        double a = 0.0, b = 0.0, c = 0.0, d = 0.0, e = 0.0, f = 0.0, s1 = 0.0, s2 = 0.0, tr = 0.0, dt = 0.0;
        if (!(fabs(det) > 1e-12)) continue;
        u1 = uv[(size_t)r->v[1] * 2] - uv[(size_t)r->v[0] * 2]; v1 = uv[(size_t)r->v[1] * 2 + 1] - uv[(size_t)r->v[0] * 2 + 1];
        u2 = uv[(size_t)r->v[2] * 2] - uv[(size_t)r->v[0] * 2]; v2 = uv[(size_t)r->v[2] * 2 + 1] - uv[(size_t)r->v[0] * 2 + 1];
        /* J = [u1 u2; v1 v2] * inv([x1 x2; y1 y2]) */
        a = (u1 * y2 - u2 * y1) / det; b = (-u1 * x2 + u2 * x1) / det;
        c = (v1 * y2 - v2 * y1) / det; d = (-v1 * x2 + v2 * x1) / det;
        au += 0.5 * fabs(u1 * v2 - u2 * v1);
        a3 += r->area;
        /* singular values of J from J^T J */
        e = a * a + c * c; f = b * b + d * d; tr = e + f;
        dt = (a * d - b * c) * (a * d - b * c);
        s1 = sqrt(fmax(0.0, 0.5 * (tr + sqrt(fmax(0.0, tr * tr - 4.0 * dt)))));
        s2 = sqrt(fmax(0.0, 0.5 * (tr - sqrt(fmax(0.0, tr * tr - 4.0 * dt)))));
        if (s1 > 0.0) v[m++] = fabs(log(s1));
        if (s2 > 0.0) v[m++] = fabs(log(s2));
    }
    if (m > 0) {
        qsort(v, m, sizeof *v, cs_dcmp);
        *p50 = v[(size_t)(0.50 * (double)(m - 1))];
        *p99 = v[(size_t)(0.99 * (double)(m - 1))];
    } else {
        *p50 = *p99 = 0.0;
    }
    if (area_uv) *area_uv = au;
    if (area_3d) *area_3d = a3;
    free(v);
}

/* The rest shape of the coarse triangle a b c (3-D corners) as a CsTri for the distortion statistics; 0 when
 * degenerate. */
static int cs_rest_tri(const float *P0, const float *P1, const float *P2, int32_t a, int32_t b, int32_t c, CsTri *r)
{
    double e1[3], e2[3], l1 = 0.0, x = 0.0, cr[3], ar = 0.0;
    int kk = 0;
    for (kk = 0; kk < 3; kk++) { e1[kk] = (double)P1[kk] - P0[kk]; e2[kk] = (double)P2[kk] - P0[kk]; }
    l1 = sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
    if (!(l1 > 1e-9)) return 0;
    x = (e2[0] * e1[0] + e2[1] * e1[1] + e2[2] * e1[2]) / l1;
    cr[0] = e1[1] * e2[2] - e1[2] * e2[1]; cr[1] = e1[2] * e2[0] - e1[0] * e2[2]; cr[2] = e1[0] * e2[1] - e1[1] * e2[0];
    ar = sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
    if (!(ar > 1e-9)) return 0;
    r->v[0] = a; r->v[1] = b; r->v[2] = c;
    r->q[0] = 0.0; r->q[1] = 0.0; r->q[2] = l1; r->q[3] = 0.0; r->q[4] = x; r->q[5] = ar / l1;
    r->area = 0.5 * ar;
    r->cot[0] = r->cot[1] = r->cot[2] = 1.0;
    return 1;
}

/* An empty lattice node at fine (fi, fj) next to the sheet: the affine continuation of the valid fine nodes of the
 * nearest one's component within +-kstep, evaluated at the node (CLOSEUP_SLIM_FRINGE_*).  1 and out/comp_out when it
 * qualifies, else 0. */
static int cs_fringe_node(const CloseupGrid *g, long kstep, long fi, long fj, float out[3], int32_t *comp_out)
{
    const long i0 = fi - kstep < 0 ? 0 : fi - kstep, j0 = fj - kstep < 0 ? 0 : fj - kstep;
    const long i1 = fi + kstep > (long)g->ni - 1 ? (long)g->ni - 1 : fi + kstep;
    const long j1 = fj + kstep > (long)g->nj - 1 ? (long)g->nj - 1 : fj + kstep;
    double M[3][3], R[3][3], W = 0.0, su = 0.0, sv = 0.0, suu = 0.0, svv = 0.0, suv = 0.0, best = 1e300, det = 0.0;
    double cuu = 0.0, cvv = 0.0, cuv = 0.0, lmin = 0.0, p[3] = { 0.0, 0.0, 0.0 }, closest[3] = { 0.0, 0.0, 0.0 }, dd = 0.0;
    long i = 0, j = 0, cnt = 0;
    int32_t m = -1;
    int a = 0, b = 0;
    if (i0 > i1 || j0 > j1) return 0;
    for (j = j0; j <= j1; j++)                  /* the nearest valid node: its component continues here */
        for (i = i0; i <= i1; i++) {
            const size_t n = (size_t)j * g->ni + (size_t)i;
            const double d2 = (double)((i - fi) * (i - fi) + (j - fj) * (j - fj));
            if (g->state[n] == CLOSEUP_NODE_NONE || g->component[n] < 0 || !isfinite(g->pos[n * 3]) || d2 >= best) continue;
            best = d2; m = g->component[n];
            closest[0] = g->pos[n * 3]; closest[1] = g->pos[n * 3 + 1]; closest[2] = g->pos[n * 3 + 2];
        }
    if (m < 0) return 0;
    memset(M, 0, sizeof M);
    memset(R, 0, sizeof R);
    for (j = j0; j <= j1; j++)
        for (i = i0; i <= i1; i++) {
            const size_t n = (size_t)j * g->ni + (size_t)i;
            const double du = (double)(i - fi), dv = (double)(j - fj);
            const double w = 1.0 / (1.0 + (du * du + dv * dv) / (double)(kstep * kstep));
            const double x[3] = { 1.0, du, dv };
            if (g->state[n] == CLOSEUP_NODE_NONE || g->component[n] != m || !isfinite(g->pos[n * 3])) continue;
            for (a = 0; a < 3; a++) {
                for (b = 0; b < 3; b++) M[a][b] += w * x[a] * x[b];
                for (b = 0; b < 3; b++) R[a][b] += w * x[a] * (double)g->pos[n * 3 + (size_t)b];
            }
            W += w; su += w * du; sv += w * dv; suu += w * du * du; svv += w * dv * dv; suv += w * du * dv;
            cnt++;
        }
    if (cnt < CLOSEUP_SLIM_FRINGE_MIN || !(W > 0.0)) return 0;
    cuu = suu / W - (su / W) * (su / W); cvv = svv / W - (sv / W) * (sv / W); cuv = suv / W - (su / W) * (sv / W);
    lmin = 0.5 * (cuu + cvv) - sqrt(0.25 * (cuu - cvv) * (cuu - cvv) + cuv * cuv);
    if (!(lmin >= CLOSEUP_SLIM_FRINGE_SPREAD)) return 0;
    /* the intercept (the value at du = dv = 0) of M p = r, per coordinate, by Cramer's rule */
    det = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
          M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
    if (!(fabs(det) > 1e-12 * W * W * W)) return 0;
    for (b = 0; b < 3; b++)
        p[b] = (R[0][b] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (R[1][b] * M[2][2] - M[1][2] * R[2][b]) +
                M[0][2] * (R[1][b] * M[2][1] - M[1][1] * R[2][b])) / det;
    for (b = 0; b < 3; b++) dd += (p[b] - closest[b]) * (p[b] - closest[b]);
    if (!(sqrt(dd) <= CLOSEUP_SLIM_FRINGE_REACH * (double)kstep * g->h) || !isfinite(dd)) return 0;
    out[0] = (float)p[0]; out[1] = (float)p[1]; out[2] = (float)p[2];
    *comp_out = m;
    return 1;
}

/* ----------------------------------------------------------------- duplicated papyrus */

typedef struct CsKeyed { uint64_t key; int32_t node; } CsKeyed;

static int cs_keyed_cmp(const void *a, const void *b)
{
    const CsKeyed *x = (const CsKeyed *)a, *y = (const CsKeyed *)b;
    if (x->key != y->key) return x->key < y->key ? -1 : 1;
    return x->node < y->node ? -1 : (x->node > y->node ? 1 : 0);
}

static uint64_t cs_cell_key(long cz, long cy, long cx)
{
    const long off = 1L << 20;
    return ((uint64_t)(cz + off) << 42) | ((uint64_t)(cy + off) << 21) | (uint64_t)(cx + off);
}

/* The first index in sorted[0, n) whose key is >= key. */
static size_t cs_lower_bound(const CsKeyed *sorted, size_t n, uint64_t key)
{
    size_t lo = 0, hi = n;
    while (lo < hi) {
        const size_t mid = lo + (hi - lo) / 2;
        if (sorted[mid].key < key) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

int CloseupSurface_dedup_fill(CloseupGrid *g, double r_vox, double frame_px, double cos_min, CloseupDedupStats *st)
{
    CloseupDedupStats s0;
    const size_t nn = g->ni * g->nj;
    const double r2 = r_vox * r_vox, f2 = frame_px * frame_px;
    size_t k = 0, nv = 0;
    long long pairs = 0, dpairs = 0;
    CsKeyed *ks = NULL;
    uint8_t *drop = NULL;
    double t0 = ves_clock_sec();
    if (st == NULL) st = &s0;
    memset(st, 0, sizeof *st);
    if (!(r_vox > 0.0)) return -1;
    for (k = 0; k < nn; k++) nv += g->state[k] != CLOSEUP_NODE_NONE && isfinite(g->pos[k * 3]);
    ks = (CsKeyed *)malloc((nv ? nv : 1) * sizeof *ks);
    drop = (uint8_t *)calloc(nn ? nn : 1, 1);
    if (ks == NULL || drop == NULL) { free(ks); free(drop); return -1; }
    for (k = 0, nv = 0; k < nn; k++) {
        if (g->state[k] == CLOSEUP_NODE_NONE || !isfinite(g->pos[k * 3])) continue;
        ks[nv].key = cs_cell_key((long)floor(g->pos[k * 3] / r_vox), (long)floor(g->pos[k * 3 + 1] / r_vox),
                                 (long)floor(g->pos[k * 3 + 2] / r_vox));
        ks[nv].node = (int32_t)k;
        nv++;
    }
    qsort(ks, nv, sizeof *ks, cs_keyed_cmp);
    {
        long q = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4096) reduction(+:pairs, dpairs)
#endif
        for (q = 0; q < (long)nv; q++) {
            const size_t n = (size_t)ks[q].node;
            const float *p = &g->pos[n * 3], *nr = &g->nrm[n * 3];
            const long cz = (long)floor(p[0] / r_vox), cy = (long)floor(p[1] / r_vox), cx = (long)floor(p[2] / r_vox);
            const double fu = g->h * (double)(n % g->ni), fv = g->h * (double)(n / g->ni);
            long dz = 0, dy = 0, dx = 0;
            for (dz = -1; dz <= 1; dz++)
                for (dy = -1; dy <= 1; dy++)
                    for (dx = -1; dx <= 1; dx++) {
                        const uint64_t key = cs_cell_key(cz + dz, cy + dy, cx + dx);
                        size_t e = cs_lower_bound(ks, nv, key);
                        for (; e < nv && ks[e].key == key; e++) {
                            const size_t m = (size_t)ks[e].node;
                            const float *pm = &g->pos[m * 3], *nm = &g->nrm[m * 3];
                            double d2 = 0.0, du = 0.0, dv = 0.0, c = 0.0;
                            int b = 0;
                            if (m == n) continue;
                            for (b = 0; b < 3; b++) d2 += ((double)p[b] - pm[b]) * ((double)p[b] - pm[b]);
                            if (d2 > r2) continue;
                            du = g->h * (double)(m % g->ni) - fu; dv = g->h * (double)(m / g->ni) - fv;
                            if (du * du + dv * dv <= f2) continue;
                            c = fabs((double)nr[0] * nm[0] + (double)nr[1] * nm[1] + (double)nr[2] * nm[2]);
                            if (c < cos_min) continue;
                            if (n < m) {
                                pairs++;
                                dpairs += g->state[n] == CLOSEUP_NODE_DATA && g->state[m] == CLOSEUP_NODE_DATA;
                            }
                            if (g->state[n] != CLOSEUP_NODE_FILL) continue;
                            if (g->state[m] != CLOSEUP_NODE_FILL || g->dist[n] > g->dist[m] ||
                                (!(g->dist[n] < g->dist[m]) && n > m))
                                drop[n] = 1;
                        }
                    }
        }
    }
    for (k = 0; k < nn; k++) {
        if (!drop[k]) continue;
        g->state[k] = CLOSEUP_NODE_NONE;
        g->pos[k * 3] = g->pos[k * 3 + 1] = g->pos[k * 3 + 2] = NAN;
        g->chart[k] = -1;
        g->component[k] = -1;
        st->dropped_fill++;
    }
    st->pairs = (size_t)pairs;
    st->data_pairs = (size_t)dpairs;
    st->seconds = ves_clock_sec() - t0;
    free(ks); free(drop);
    return 0;
}

/* Steps from node k back along the grid line (-di, -dj) to the nearest DATA node (0 on DATA itself), for every node;
 * cap when there is none within cap - 1 steps.  The predecessor (i - di, j - dj) is always visited first: rows in the
 * order of dj, and within a row of dj == 0 the columns in the order of di. */
static void cs_line_dist(const CloseupGrid *g, int di, int dj, uint16_t cap, uint16_t *out)
{
    const long ni = (long)g->ni, nj = (long)g->nj;
    long jj = 0, ii = 0;
    for (jj = 0; jj < nj; jj++) {
        const long j = dj >= 0 ? jj : nj - 1 - jj;
        for (ii = 0; ii < ni; ii++) {
            const long i = di >= 0 ? ii : ni - 1 - ii;
            const long pi = i - di, pj = j - dj;
            const size_t k = (size_t)j * (size_t)ni + (size_t)i;
            uint16_t p = cap;
            if (g->state[k] == CLOSEUP_NODE_DATA) { out[k] = 0; continue; }
            if (pi >= 0 && pj >= 0 && pi < ni && pj < nj) p = out[(size_t)pj * (size_t)ni + (size_t)pi];
            out[k] = p >= (uint16_t)(cap - 1) ? cap : (uint16_t)(p + 1);
        }
    }
}

int CloseupSurface_trim_apron(CloseupGrid *g, double reach_px, CloseupApronStats *st)
{
    static const int axis[4][2] = { { 1, 0 }, { 0, 1 }, { 1, 1 }, { 1, -1 } };
    CloseupApronStats s0;
    const size_t nn = g->ni * g->nj;
    uint16_t *fw = NULL, *bw = NULL, cap = 0;
    uint8_t *keep = NULL;
    size_t k = 0;
    long r = 0;
    int a = 0;
    double t0 = ves_clock_sec();
    if (st == NULL) st = &s0;
    memset(st, 0, sizeof *st);
    r = lround(reach_px / g->h);
    if (r < 1 || nn == 0) return -1;
    if (r > 60000) r = 60000;
    cap = (uint16_t)(r + 1);
    fw = (uint16_t *)malloc(nn * sizeof *fw);
    bw = (uint16_t *)malloc(nn * sizeof *bw);
    keep = (uint8_t *)calloc(nn, 1);
    if (fw == NULL || bw == NULL || keep == NULL) { free(fw); free(bw); free(keep); return -1; }
    for (a = 0; a < 4; a++) {
        cs_line_dist(g, axis[a][0], axis[a][1], cap, fw);
        cs_line_dist(g, -axis[a][0], -axis[a][1], cap, bw);
        for (k = 0; k < nn; k++) keep[k] |= (uint8_t)(fw[k] <= (uint16_t)r && bw[k] <= (uint16_t)r);
    }
    for (k = 0; k < nn; k++) {
        if (g->state[k] != CLOSEUP_NODE_FILL) continue;
        st->fill_before++;
        if (keep[k]) continue;
        g->state[k] = CLOSEUP_NODE_NONE;
        g->pos[k * 3] = g->pos[k * 3 + 1] = g->pos[k * 3 + 2] = NAN;
        g->chart[k] = -1;
        g->component[k] = -1;
        st->dropped++;
    }
    st->seconds = ves_clock_sec() - t0;
    free(fw); free(bw); free(keep);
    return 0;
}

/* CloseupSurface_flatten_slim with the fringe continuation on or off (the selftest measures both). */
static int cs_flatten_slim(const CloseupGrid *g, double coarse_px, int iters, int fringe, double *uv_out, CloseupFlatStats *st)
{
    CloseupFlatStats s0;
    const size_t nn = g->ni * g->nj;
    size_t k = 0, nci = 0, ncj = 0, ncv = 0, nsrc = 0, nscf = 0, ncell = 0, cell = 0, nvalid = 0, nall = 0;
    long kstep = 0;
    long long nfringe = 0;
    uint8_t *valid = NULL, *split = NULL, *fixed = NULL, *ext = NULL, *allv = NULL;
    uint32_t *cellcnt = NULL;
    int32_t *comp = NULL, *fsrc = NULL, *fscf = NULL, *faces = NULL;
    float *xyz = NULL;
    double *uv0 = NULL, *uv = NULL;
    CsTri *T = NULL, *Tv = NULL;
    Arena_T arena = NULL;
    SlimRefineOpts o;
    SlimRefineStats ss;
    double t0 = ves_clock_sec();
    int rc = -1, solved = 0;
    if (st == NULL) st = &s0;
    memset(st, 0, sizeof *st);
    memset(&ss, 0, sizeof ss);
    for (k = 0; k < nn * 2; k++) uv_out[k] = NAN;
    if (g->ni < 2 || g->nj < 2) return -1;
    kstep = (long)floor(coarse_px / g->h + 0.5);
    if (kstep < 1) kstep = 1;
    /* the coarse lattice covers every fine node: its last row and column may lie past the grid (scaffold vertices) */
    nci = (g->ni - 1 + (size_t)kstep - 1) / (size_t)kstep + 1;
    ncj = (g->nj - 1 + (size_t)kstep - 1) / (size_t)kstep + 1;
    ncv = nci * ncj;
    ncell = (nci - 1) * (ncj - 1);
    valid = (uint8_t *)calloc(ncv, 1);
    comp = (int32_t *)malloc(ncv * sizeof *comp);
    xyz = (float *)malloc(ncv * 3 * sizeof *xyz);
    uv0 = (double *)malloc(ncv * 2 * sizeof *uv0);
    uv = (double *)malloc(ncv * 2 * sizeof *uv);
    split = (uint8_t *)calloc(ncell ? ncell : 1, 1);
    fsrc = (int32_t *)malloc((ncell ? ncell : 1) * 2 * 3 * sizeof *fsrc);
    fscf = (int32_t *)malloc((ncell ? ncell : 1) * 2 * 3 * sizeof *fscf);
    T = (CsTri *)malloc((ncell ? ncell : 1) * 2 * sizeof *T);
    allv = (uint8_t *)calloc((ncell ? ncell : 1) * 2, 1);
    ext = (uint8_t *)calloc(ncv, 1);
    cellcnt = (uint32_t *)calloc(ncell ? ncell : 1, sizeof *cellcnt);
    if (!valid || !comp || !xyz || !uv0 || !uv || !split || !fsrc || !fscf || !T || !allv || !ext || !cellcnt) goto done;
    for (k = 0; k < ncv; k++) {
        const size_t ci = k % nci, cj = k / nci, fi = ci * (size_t)kstep, fj = cj * (size_t)kstep;
        const double u = g->u0 + g->h * (double)fi, v = g->v0 + g->h * (double)fj;
        uv0[k * 2] = u; uv0[k * 2 + 1] = v;
        comp[k] = -1;
        if (fi < g->ni && fj < g->nj) {
            const size_t n = fj * g->ni + fi;
            if (g->state[n] != CLOSEUP_NODE_NONE && isfinite(g->pos[n * 3]) && isfinite(g->pos[n * 3 + 1]) &&
                isfinite(g->pos[n * 3 + 2])) {
                valid[k] = 1;
                comp[k] = g->component[n];
                memcpy(&xyz[k * 3], &g->pos[n * 3], 3 * sizeof *xyz);
                nvalid++;
                continue;
            }
        }
        xyz[k * 3] = 0.0f; xyz[k * 3 + 1] = (float)v; xyz[k * 3 + 2] = (float)u;   /* scaffold only: no 3-D role */
    }
    if (nvalid < 3) goto done;
    for (k = 0; k < nn; k++) {                   /* valid fine nodes per coarse cell (the one that carries them) */
        size_t ci = (k % g->ni) / (size_t)kstep, cj = (k / g->ni) / (size_t)kstep;
        if (g->state[k] == CLOSEUP_NODE_NONE) continue;
        if (ci > nci - 2) ci = nci - 2;
        if (cj > ncj - 2) cj = ncj - 2;
        cellcnt[cj * (nci - 1) + ci]++;
    }
    if (fringe) {                                /* the fringe: empty lattice nodes continued from the sheet */
        long q = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 512) reduction(+:nfringe)
#endif
        for (q = 0; q < (long)ncv; q++) {
            const long ci = q % (long)nci, cj = q / (long)nci;
            float p[3];
            int32_t m = -1;
            if (valid[q]) continue;
            if (cs_fringe_node(g, kstep, ci * kstep, cj * kstep, p, &m)) {
                memcpy(&xyz[(size_t)q * 3], p, sizeof p);
                comp[q] = m;
                ext[q] = 1;
                nfringe++;
            }
        }
    }
    /* every coarse cell: two triangles on the diagonal a-d (split 0) or b-c (split 1), whichever carries more surface
     * triangles (three valid or continued corners of one component, not degenerate, the cell holding data when no
     * corner is valid); the others are scaffold */
    for (cell = 0; cell < ncell; cell++) {
        const size_t ci = cell % (nci - 1), cj = cell / (nci - 1);
        const int32_t a = (int32_t)(cj * nci + ci), b = a + 1, c = a + (int32_t)nci, d = c + 1;
        int32_t tri[2][2][3];
        int good[2][2], sp = 0, q = 0, best = 0;
        tri[0][0][0] = a; tri[0][0][1] = b; tri[0][0][2] = d;       /* split 0: the diagonal a-d */
        tri[0][1][0] = a; tri[0][1][1] = d; tri[0][1][2] = c;
        tri[1][0][0] = a; tri[1][0][1] = b; tri[1][0][2] = c;       /* split 1: the diagonal b-c */
        tri[1][1][0] = b; tri[1][1][1] = d; tri[1][1][2] = c;
        for (sp = 0; sp < 2; sp++)
            for (q = 0; q < 2; q++) {
                const int32_t *t = tri[sp][q];
                CsTri tmp;
                memset(&tmp, 0, sizeof tmp);
                good[sp][q] = (valid[t[0]] || ext[t[0]]) && (valid[t[1]] || ext[t[1]]) && (valid[t[2]] || ext[t[2]]) &&
                              (valid[t[0]] || valid[t[1]] || valid[t[2]] || cellcnt[cell] > 0) &&
                              comp[t[0]] >= 0 && comp[t[0]] == comp[t[1]] &&
                              comp[t[0]] == comp[t[2]] && cs_rest_tri(&xyz[(size_t)t[0] * 3], &xyz[(size_t)t[1] * 3],
                                                                      &xyz[(size_t)t[2] * 3], t[0], t[1], t[2], &tmp);
            }
        best = (good[1][0] + good[1][1] > good[0][0] + good[0][1]) ? 1 : 0;
        split[cell] = (uint8_t)best;
        for (q = 0; q < 2; q++) {
            const int32_t *t = tri[best][q];
            if (good[best][q]) {
                (void)cs_rest_tri(&xyz[(size_t)t[0] * 3], &xyz[(size_t)t[1] * 3], &xyz[(size_t)t[2] * 3], t[0], t[1], t[2],
                                  &T[nsrc]);
                memcpy(&fsrc[nsrc * 3], t, 3 * sizeof *t);
                allv[nsrc] = (uint8_t)(valid[t[0]] && valid[t[1]] && valid[t[2]]);
                nall += allv[nsrc];
                nsrc++;
            } else {
                memcpy(&fscf[nscf * 3], t, 3 * sizeof *t);
                nscf++;
            }
        }
    }
    st->coarse_nodes = nvalid;
    st->coarse_tris = nsrc;
    st->coarse_fringe = (size_t)nfringe;
    if (nsrc == 0) goto done;
    faces = (int32_t *)malloc((nsrc + nscf) * 3 * sizeof *faces);
    fixed = (uint8_t *)calloc(ncv, 1);
    Tv = (CsTri *)malloc((nall ? nall : 1) * sizeof *Tv);   /* the statistics: triangles of valid corners only */
    if (faces == NULL || fixed == NULL || Tv == NULL) goto done;
    memcpy(faces, fsrc, nsrc * 3 * sizeof *faces);
    memcpy(faces + nsrc * 3, fscf, nscf * 3 * sizeof *faces);
    for (k = 0, cell = 0; k < nsrc; k++)
        if (allv[k]) Tv[cell++] = T[k];
    cs_distortion(Tv, nall, uv0, &st->logs_p50_before, &st->logs_p99_before, NULL, NULL);
    memcpy(uv, uv0, ncv * 2 * sizeof *uv);
    arena = Arena_new();
    if (arena != NULL) {
        SlimRefine_defaults(&o);
        o.iterations = iters;
        o.pcg_iterations = 240;
        o.strip_pack = 0;
        o.guard_boundary = 1;
        o.use_amg = 1;
        o.verbose = 0;
        o.threads = 0;
        if (nscf > 0) {               /* SCAF: the scaffold keeps its frame shape as rest; vertex 0 holds the gauge */
            o.scaffold_uv = uv0;
            o.scaffold_first_face = nsrc;
            o.local_steps = 0;
            o.pcg_tolerance = 1e-8;
            fixed[0] = 1;
            o.fixed_vertices = fixed;
        }
        if (SlimRefine_run_double(arena, xyz, ncv, faces, nsrc + nscf, uv0, &o, uv, &ss) == 0 && ss.flips_after == 0 &&
            ss.boundary_intersections_after == 0) {
            for (k = 0; k < ncv * 2 && isfinite(uv[k]); k++) {}
            solved = k == ncv * 2;
        }
    }
    if (!solved) {
        fprintf(stderr, "[closeup] SLIM flatten not accepted (flips %zu, boundary %zu): no replacement\n", ss.flips_after,
                ss.boundary_intersections_after);
        goto done;
    }
    st->components = ss.components;
    st->iterations = solved ? ss.iterations_run : 0;
    st->energy_first = ss.energy_before;
    st->energy_last = solved ? ss.energy_after : ss.energy_before;
    /* rigid alignment (rotation + translation over the surface's coarse nodes) back onto the frame */
    {
        double mu0[2] = { 0.0, 0.0 }, mu1[2] = { 0.0, 0.0 }, sxx = 0.0, sxy = 0.0, th = 0.0, cs = 0.0, sn = 0.0;
        for (k = 0; k < ncv; k++) {
            if (!valid[k]) continue;
            mu0[0] += uv0[k * 2]; mu0[1] += uv0[k * 2 + 1];
            mu1[0] += uv[k * 2]; mu1[1] += uv[k * 2 + 1];
        }
        mu0[0] /= (double)nvalid; mu0[1] /= (double)nvalid; mu1[0] /= (double)nvalid; mu1[1] /= (double)nvalid;
        for (k = 0; k < ncv; k++) {
            double ax = 0.0, ay = 0.0, bx = 0.0, by = 0.0;
            if (!valid[k]) continue;
            ax = uv[k * 2] - mu1[0]; ay = uv[k * 2 + 1] - mu1[1];
            bx = uv0[k * 2] - mu0[0]; by = uv0[k * 2 + 1] - mu0[1];
            sxx += ax * bx + ay * by;
            sxy += ax * by - ay * bx;
        }
        th = atan2(sxy, sxx);
        cs = cos(th); sn = sin(th);
        for (k = 0; k < ncv; k++) {
            const double ax = uv[k * 2] - mu1[0], ay = uv[k * 2 + 1] - mu1[1];
            uv[k * 2] = mu0[0] + cs * ax - sn * ay;
            uv[k * 2 + 1] = mu0[1] + sn * ax + cs * ay;
        }
    }
    {
        double au = 0.0, a3 = 0.0;
        cs_distortion(Tv, nall, uv, &st->logs_p50_after, &st->logs_p99_after, &au, &a3);
        st->area_ratio = a3 > 0.0 ? au / a3 : 0.0;
    }
    /* fine nodes: barycentric in the coarse triangle (surface or scaffold) that holds them in the frame */
    {
        long n = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (n = 0; n < (long)nn; n++) {
            const size_t i = (size_t)n % g->ni, j = (size_t)n / g->ni;
            size_t ci = i / (size_t)kstep, cj = j / (size_t)kstep;
            double fx = 0.0, fy = 0.0, w[3];
            int32_t a = 0, b = 0, c = 0, d = 0, t[3];
            int q = 0;
            if (g->state[n] == CLOSEUP_NODE_NONE) continue;
            if (ci > nci - 2) ci = nci - 2;
            if (cj > ncj - 2) cj = ncj - 2;
            fx = (double)((long)i - (long)(ci * (size_t)kstep)) / (double)kstep;
            fy = (double)((long)j - (long)(cj * (size_t)kstep)) / (double)kstep;
            a = (int32_t)(cj * nci + ci); b = a + 1; c = a + (int32_t)nci; d = c + 1;
            if (split[cj * (nci - 1) + ci] == 0) {
                if (fx >= fy) { t[0] = a; t[1] = b; t[2] = d; w[0] = 1.0 - fx; w[1] = fx - fy; w[2] = fy; }
                else { t[0] = a; t[1] = d; t[2] = c; w[0] = 1.0 - fy; w[1] = fx; w[2] = fy - fx; }
            } else {
                if (fx + fy <= 1.0) { t[0] = a; t[1] = b; t[2] = c; w[0] = 1.0 - fx - fy; w[1] = fx; w[2] = fy; }
                else { t[0] = b; t[1] = d; t[2] = c; w[0] = 1.0 - fy; w[1] = fx + fy - 1.0; w[2] = 1.0 - fx; }
            }
            uv_out[(size_t)n * 2] = 0.0; uv_out[(size_t)n * 2 + 1] = 0.0;
            for (q = 0; q < 3; q++) {
                uv_out[(size_t)n * 2] += w[q] * uv[(size_t)t[q] * 2];
                uv_out[(size_t)n * 2 + 1] += w[q] * uv[(size_t)t[q] * 2 + 1];
            }
        }
    }
    st->seconds = ves_clock_sec() - t0;
    rc = 0;
done:
    if (arena != NULL) Arena_dispose(&arena);
    free(valid); free(comp); free(xyz); free(uv0); free(uv); free(split); free(fsrc); free(fscf); free(T); free(faces);
    free(fixed); free(ext); free(allv); free(cellcnt); free(Tv);
    return rc;
}

int CloseupSurface_flatten_slim(const CloseupGrid *g, double coarse_px, int iters, double *uv_out, CloseupFlatStats *st)
{
    return cs_flatten_slim(g, coarse_px, iters, 1, uv_out, st);
}

int CloseupSurface_flatten(const CloseupGrid *g, double coarse_px, int iters, double tol, double *uv_out,
                           CloseupFlatStats *st)
{
    /* All callers, including the ordinary close-up seed, use the injective
     * SLIM path. Keep the old tolerance argument for source compatibility. */
    (void)tol;
    return CloseupSurface_flatten_slim(g, coarse_px, iters, uv_out, st);
}

#if 0 /* Historical ARAP control: never compiled or called by the pipeline. */
static int cs_flatten_arap_reference(const CloseupGrid *g, double coarse_px, int iters, double tol, double *uv_out,
                                      CloseupFlatStats *st)
{
    CloseupFlatStats s0;
    size_t nn = g->ni * g->nj, k = 0, nci = 0, ncj = 0, ncv = 0, nt = 0, t = 0;
    long kstep = 0;
    int32_t *cidx = NULL, *cnode = NULL;
    CsTri *T = NULL;
    double *uv0 = NULL, *uv = NULL, *rhs = NULL, *diag = NULL, *uv_best = NULL, e_best = 0.0;
    int *rows = NULL, *cols = NULL;
    double *vals = NULL;
    size_t ntrip = 0;
    SparseFactor_T F = NULL;
    double t0 = ves_clock_sec(), e_prev = 0.0, mean_diag = 0.0, eps = 0.0;
    int it = 0, rc = -1;
    if (st == NULL) st = &s0;
    memset(st, 0, sizeof *st);
    for (k = 0; k < nn * 2; k++) uv_out[k] = NAN;
    kstep = (long)floor(coarse_px / g->h + 0.5);
    if (kstep < 1) kstep = 1;
    nci = (g->ni - 1) / (size_t)kstep + 1;
    ncj = (g->nj - 1) / (size_t)kstep + 1;
    cidx = (int32_t *)malloc(nci * ncj * sizeof *cidx);
    if (cidx == NULL) return -1;
    for (k = 0; k < nci * ncj; k++) {
        size_t ci = k % nci, cj = k / nci, n = cj * (size_t)kstep * g->ni + ci * (size_t)kstep;
        cidx[k] = g->state[n] != CLOSEUP_NODE_NONE ? (int32_t)ncv++ : -1;
    }
    if (ncv < 3) { free(cidx); return -1; }
    cnode = (int32_t *)malloc(ncv * sizeof *cnode);
    uv0 = (double *)malloc(ncv * 2 * sizeof *uv0);
    uv = (double *)malloc(ncv * 2 * sizeof *uv);
    rhs = (double *)malloc(ncv * 2 * sizeof *rhs);
    uv_best = (double *)malloc(ncv * 2 * sizeof *uv_best);
    diag = (double *)calloc(ncv, sizeof *diag);
    T = (CsTri *)malloc(2 * nci * ncj * sizeof *T);
    if (cnode == NULL || uv0 == NULL || uv == NULL || rhs == NULL || uv_best == NULL || diag == NULL || T == NULL)
        goto done;
    for (k = 0; k < nci * ncj; k++) {
        size_t ci = k % nci, cj = k / nci;
        if (cidx[k] < 0) continue;
        cnode[cidx[k]] = (int32_t)(cj * (size_t)kstep * g->ni + ci * (size_t)kstep);
        uv0[(size_t)cidx[k] * 2] = g->u0 + g->h * (double)(ci * (size_t)kstep);
        uv0[(size_t)cidx[k] * 2 + 1] = g->v0 + g->h * (double)(cj * (size_t)kstep);
    }
    /* triangles: diagonal a-d when all four corners exist, else the one triangle of any three */
    for (k = 0; k + nci < nci * ncj; k++) {
        size_t ci = k % nci;
        int32_t a = 0, b = 0, c = 0, d = 0;
        int32_t tri[2][3];
        int ntri = 0, q = 0;
        if (ci + 1 >= nci) continue;
        a = cidx[k]; b = cidx[k + 1]; c = cidx[k + nci]; d = cidx[k + nci + 1];
        if (a >= 0 && b >= 0 && c >= 0 && d >= 0) {
            tri[0][0] = a; tri[0][1] = b; tri[0][2] = d;
            tri[1][0] = a; tri[1][1] = d; tri[1][2] = c;
            ntri = 2;
        } else if (a >= 0 && b >= 0 && c >= 0) { tri[0][0] = a; tri[0][1] = b; tri[0][2] = c; ntri = 1; }
        else if (b >= 0 && d >= 0 && c >= 0) { tri[0][0] = b; tri[0][1] = d; tri[0][2] = c; ntri = 1; }
        else if (a >= 0 && d >= 0 && c >= 0) { tri[0][0] = a; tri[0][1] = d; tri[0][2] = c; ntri = 1; }
        else if (a >= 0 && b >= 0 && d >= 0) { tri[0][0] = a; tri[0][1] = b; tri[0][2] = d; ntri = 1; }
        for (q = 0; q < ntri; q++) {
            CsTri *r = &T[nt];
            const float *P0 = &g->pos[(size_t)cnode[tri[q][0]] * 3];
            const float *P1 = &g->pos[(size_t)cnode[tri[q][1]] * 3];
            const float *P2 = &g->pos[(size_t)cnode[tri[q][2]] * 3];
            double e1[3], e2[3], l1 = 0.0, x = 0.0, y2 = 0.0, cr[3], ar = 0.0;
            int kk = 0, cc = 0;
            int32_t owner = g->component[cnode[tri[q][0]]];
            if (owner < 0 || owner != g->component[cnode[tri[q][1]]] || owner != g->component[cnode[tri[q][2]]]) continue;
            for (kk = 0; kk < 3; kk++) { e1[kk] = (double)P1[kk] - P0[kk]; e2[kk] = (double)P2[kk] - P0[kk]; }
            l1 = sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
            if (!(l1 > 1e-9)) continue;
            x = (e2[0] * e1[0] + e2[1] * e1[1] + e2[2] * e1[2]) / l1;
            cr[0] = e1[1] * e2[2] - e1[2] * e2[1]; cr[1] = e1[2] * e2[0] - e1[0] * e2[2]; cr[2] = e1[0] * e2[1] - e1[1] * e2[0];
            ar = sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
            y2 = ar / l1;
            if (!(ar > 1e-9)) continue;
            r->v[0] = tri[q][0]; r->v[1] = tri[q][1]; r->v[2] = tri[q][2];
            r->q[0] = 0.0; r->q[1] = 0.0; r->q[2] = l1; r->q[3] = 0.0; r->q[4] = x; r->q[5] = y2;
            r->area = 0.5 * ar;
            for (cc = 0; cc < 3; cc++) {
                /* angle at corner cc between the edges to the other two corners */
                int o1 = (cc + 1) % 3, o2 = (cc + 2) % 3;
                double ax = r->q[o1 * 2] - r->q[cc * 2], ay = r->q[o1 * 2 + 1] - r->q[cc * 2 + 1];
                double bx = r->q[o2 * 2] - r->q[cc * 2], by = r->q[o2 * 2 + 1] - r->q[cc * 2 + 1];
                double dot = ax * bx + ay * by, crs = fabs(ax * by - ay * bx), ct = crs > 1e-12 ? dot / crs : 1e3;
                if (ct < 1e-3) ct = 1e-3;
                if (ct > 1e3) ct = 1e3;
                r->cot[cc] = ct;
            }
            nt++;
        }
    }
    st->coarse_nodes = ncv;
    st->coarse_tris = nt;
    if (nt == 0) goto done;
    /* the cotangent Laplacian (edge k joins corners k+1, k+2, weighted cot(angle at k) / 2) + eps anchor */
    rows = (int *)malloc((nt * 3 + ncv) * sizeof *rows);
    cols = (int *)malloc((nt * 3 + ncv) * sizeof *cols);
    vals = (double *)malloc((nt * 3 + ncv) * sizeof *vals);
    if (rows == NULL || cols == NULL || vals == NULL) goto done;
    for (t = 0; t < nt; t++) {
        int cc = 0;
        for (cc = 0; cc < 3; cc++) {
            int32_t a = T[t].v[(cc + 1) % 3], b = T[t].v[(cc + 2) % 3];
            double w = 0.5 * T[t].cot[cc];
            diag[a] += w; diag[b] += w;
            rows[ntrip] = a > b ? a : b; cols[ntrip] = a > b ? b : a; vals[ntrip] = -w; ntrip++;
        }
    }
    for (k = 0; k < ncv; k++) mean_diag += diag[k];
    mean_diag /= (double)ncv;
    eps = 1e-7 * (mean_diag > 0.0 ? mean_diag : 1.0);
    for (k = 0; k < ncv; k++) { rows[ntrip] = (int)k; cols[ntrip] = (int)k; vals[ntrip] = diag[k] + eps; ntrip++; }
    if (Sparse_factor_spd((int)ncv, (int)ntrip, rows, cols, vals, &F) != 0) {
        fprintf(stderr, "[closeup] ARAP factorisation failed (%zu nodes)\n", ncv);
        goto done;
    }
    memcpy(uv, uv0, ncv * 2 * sizeof *uv);
    cs_distortion(T, nt, uv0, &st->logs_p50_before, &st->logs_p99_before, NULL, NULL);
    for (it = 0; it < iters; it++) {
        double E = 0.0;
        long lt = 0;
        memset(rhs, 0, ncv * 2 * sizeof *rhs);
        for (k = 0; k < ncv; k++) { rhs[k] = eps * uv0[k * 2]; rhs[ncv + k] = eps * uv0[k * 2 + 1]; }
        /* local: per triangle, the rotation closest to its weighted edge covariance S = sum w d e^T (d the current
         * UV edge, e the rest edge) -- the exact minimiser of the triangle's (clamped-cotangent) energy, so the
         * alternation never climbs; global right-hand side accumulated serially */
        for (lt = 0; lt < (long)nt; lt++) {
            const CsTri *r = &T[lt];
            double S00 = 0.0, S01 = 0.0, S10 = 0.0, S11 = 0.0, th, cs, sn;
            int cc = 0;
            for (cc = 0; cc < 3; cc++) {
                int ia = (cc + 1) % 3, ib = (cc + 2) % 3;
                double w = 0.5 * r->cot[cc];
                double ex = r->q[ia * 2] - r->q[ib * 2], ey = r->q[ia * 2 + 1] - r->q[ib * 2 + 1];
                double dx = uv[(size_t)r->v[ia] * 2] - uv[(size_t)r->v[ib] * 2];
                double dy = uv[(size_t)r->v[ia] * 2 + 1] - uv[(size_t)r->v[ib] * 2 + 1];
                S00 += w * dx * ex; S01 += w * dx * ey; S10 += w * dy * ex; S11 += w * dy * ey;
            }
            th = CloseupSurface_rot_angle(S00, S01, S10, S11);
            cs = cos(th); sn = sin(th);
            for (cc = 0; cc < 3; cc++) {
                int ia = (cc + 1) % 3, ib = (cc + 2) % 3;
                int32_t va = r->v[ia], vb = r->v[ib];
                double w = 0.5 * r->cot[cc];
                double ex = r->q[ia * 2] - r->q[ib * 2], ey = r->q[ia * 2 + 1] - r->q[ib * 2 + 1];
                double rx = cs * ex - sn * ey, ry = sn * ex + cs * ey;
                double du = uv[(size_t)va * 2] - uv[(size_t)vb * 2] - rx, dv = uv[(size_t)va * 2 + 1] - uv[(size_t)vb * 2 + 1] - ry;
                E += w * (du * du + dv * dv);
                rhs[va] += w * rx; rhs[ncv + (size_t)va] += w * ry;
                rhs[vb] -= w * rx; rhs[ncv + (size_t)vb] -= w * ry;
            }
        }
        if (it == 0) st->energy_first = E;
        st->energy_last = E;
        if (it < 64) st->energy_hist[it] = E;
        st->iterations = (size_t)it + 1;
        if (it == 0 || E < e_best) { e_best = E; memcpy(uv_best, uv, ncv * 2 * sizeof *uv_best); }
        if (it > 0 && (fabs(e_prev - E) <= tol * (e_prev > 0.0 ? e_prev : 1.0) || E > e_prev * (1.0 + 1e-9))) break;
        e_prev = E;
        {
            double *X = (double *)malloc(ncv * 2 * sizeof *X);
            if (X == NULL) goto done;
            if (Sparse_factor_solve_multi(F, rhs, X, 2) != 0) { free(X); goto done; }
            for (k = 0; k < ncv; k++) { uv[k * 2] = X[k]; uv[k * 2 + 1] = X[ncv + k]; }
            free(X);
        }
    }
    memcpy(uv, uv_best, ncv * 2 * sizeof *uv);
    st->energy_last = e_best;
    /* rigid alignment (rotation + translation, area-uniform over nodes) back onto the frame */
    {
        double mu0[2] = { 0.0, 0.0 }, mu1[2] = { 0.0, 0.0 }, sxx = 0.0, sxy = 0.0, th = 0.0, cs = 0.0, sn = 0.0;
        for (k = 0; k < ncv; k++) {
            mu0[0] += uv0[k * 2]; mu0[1] += uv0[k * 2 + 1];
            mu1[0] += uv[k * 2]; mu1[1] += uv[k * 2 + 1];
        }
        mu0[0] /= (double)ncv; mu0[1] /= (double)ncv; mu1[0] /= (double)ncv; mu1[1] /= (double)ncv;
        for (k = 0; k < ncv; k++) {
            double ax = uv[k * 2] - mu1[0], ay = uv[k * 2 + 1] - mu1[1];
            double bx = uv0[k * 2] - mu0[0], by = uv0[k * 2 + 1] - mu0[1];
            sxx += ax * bx + ay * by;
            sxy += ax * by - ay * bx;
        }
        th = atan2(sxy, sxx);
        cs = cos(th); sn = sin(th);
        for (k = 0; k < ncv; k++) {
            double ax = uv[k * 2] - mu1[0], ay = uv[k * 2 + 1] - mu1[1];
            uv[k * 2] = mu0[0] + cs * ax - sn * ay;
            uv[k * 2 + 1] = mu0[1] + sn * ax + cs * ay;
        }
    }
    {
        double au = 0.0, a3 = 0.0;
        cs_distortion(T, nt, uv, &st->logs_p50_after, &st->logs_p99_after, &au, &a3);
        st->area_ratio = a3 > 0.0 ? au / a3 : 0.0;
    }
    /* fine nodes: the coarse displacement, bilinear over the valid corners of their coarse cell (renormalised),
     * else the nearest valid coarse node within two cells */
    {
        long n = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (n = 0; n < (long)nn; n++) {
            size_t i = (size_t)n % g->ni, j = (size_t)n / g->ni, ci = i / (size_t)kstep, cj = j / (size_t)kstep;
            double fx = (double)(i % (size_t)kstep) / (double)kstep, fy = (double)(j % (size_t)kstep) / (double)kstep;
            double wsum = 0.0, dx = 0.0, dy = 0.0;
            int q = 0;
            if (g->state[n] == CLOSEUP_NODE_NONE) continue;
            for (q = 0; q < 4; q++) {
                size_t ai = ci + (size_t)(q & 1), aj = cj + (size_t)(q >> 1);
                double w = ((q & 1) ? fx : 1.0 - fx) * ((q >> 1) ? fy : 1.0 - fy);
                int32_t c = 0;
                if (ai >= nci || aj >= ncj || w <= 0.0) continue;
                c = cidx[aj * nci + ai];
                if (c < 0) continue;
                dx += w * (uv[(size_t)c * 2] - uv0[(size_t)c * 2]);
                dy += w * (uv[(size_t)c * 2 + 1] - uv0[(size_t)c * 2 + 1]);
                wsum += w;
            }
            if (!(wsum > 1e-9)) {
                long r = 0, bi = 0, bj = 0;
                double best = 1e300;
                int32_t bc = -1;
                for (r = 1; r <= 2 && bc < 0; r++)
                    for (bj = (long)cj - r; bj <= (long)cj + r + 1; bj++)
                        for (bi = (long)ci - r; bi <= (long)ci + r + 1; bi++) {
                            int32_t c = 0;
                            double ddx = 0.0, ddy = 0.0;
                            if (bi < 0 || bj < 0 || (size_t)bi >= nci || (size_t)bj >= ncj) continue;
                            c = cidx[(size_t)bj * nci + (size_t)bi];
                            if (c < 0) continue;
                            ddx = (double)bi * (double)kstep - (double)i; ddy = (double)bj * (double)kstep - (double)j;
                            if (ddx * ddx + ddy * ddy < best) { best = ddx * ddx + ddy * ddy; bc = c; }
                        }
                if (bc < 0) continue;
                dx = uv[(size_t)bc * 2] - uv0[(size_t)bc * 2];
                dy = uv[(size_t)bc * 2 + 1] - uv0[(size_t)bc * 2 + 1];
                wsum = 1.0;
            } else {
                dx /= wsum; dy /= wsum;
            }
            uv_out[(size_t)n * 2] = g->u0 + g->h * (double)i + dx;
            uv_out[(size_t)n * 2 + 1] = g->v0 + g->h * (double)j + dy;
        }
    }
    {
        /* connected components of the coarse mesh (for the report) */
        int32_t *par = (int32_t *)malloc(ncv * sizeof *par);
        if (par != NULL) {
            size_t comps = 0;
            for (k = 0; k < ncv; k++) par[k] = (int32_t)k;
            for (t = 0; t < nt; t++) {
                int cc = 0;
                for (cc = 1; cc < 3; cc++) {
                    int32_t a = T[t].v[0], b = T[t].v[cc];
                    while (par[a] != a) { par[a] = par[par[a]]; a = par[a]; }
                    while (par[b] != b) { par[b] = par[par[b]]; b = par[b]; }
                    if (a != b) par[a] = b;
                }
            }
            for (k = 0; k < ncv; k++) if (par[k] == (int32_t)k) comps++;
            st->components = comps;
            free(par);
        }
    }
    rc = 0;
done:
    if (F != NULL) Sparse_factor_free(&F);
    free(cidx); free(cnode); free(uv0); free(uv); free(rhs); free(uv_best); free(diag); free(T);
    free(rows); free(cols); free(vals);
    st->seconds = ves_clock_sec() - t0;
    return rc;
}
#endif

/* ------------------------------------------------------------------ mesh */

static void cs_mesh_triangle(const CloseupGrid *g, const int32_t *node_of, int32_t a, int32_t b,
                             int32_t c, int32_t *faces, size_t *nf)
{
    int32_t ca = g->component[node_of[a]], cb = g->component[node_of[b]], cc = g->component[node_of[c]];
    /* Nearby extrapolated pieces may meet in UV, but they must never acquire a connecting face. */
    if (ca < 0 || ca != cb || ca != cc) return;
    faces[*nf * 3] = a; faces[*nf * 3 + 1] = b; faces[*nf * 3 + 2] = c;
    (*nf)++;
}

int CloseupSurface_mesh(const CloseupGrid *g, const double *uv_nodes, float **verts, double **uv, float **nrm,
                        int32_t **faces, int32_t **node_of_vert, size_t *nv, size_t *nf)
{
    size_t nn = g->ni * g->nj, k = 0, v = 0, f = 0;
    int32_t *vid = (int32_t *)malloc(nn * sizeof *vid);
    float *V = NULL, *N = NULL;
    double *UV = NULL;
    int32_t *F = NULL, *NV = NULL;
    *verts = NULL; *uv = NULL; *nrm = NULL; *faces = NULL; *node_of_vert = NULL; *nv = 0; *nf = 0;
    if (vid == NULL) return -1;
    for (k = 0; k < nn; k++) {
        vid[k] = -1;
        if (g->state[k] != CLOSEUP_NODE_NONE && isfinite(uv_nodes[k * 2]) && isfinite(uv_nodes[k * 2 + 1]))
            vid[k] = (int32_t)v++;
    }
    V = (float *)malloc((v ? v : 1) * 3 * sizeof *V);
    N = (float *)malloc((v ? v : 1) * 3 * sizeof *N);
    UV = (double *)malloc((v ? v : 1) * 2 * sizeof *UV);
    NV = (int32_t *)malloc((v ? v : 1) * sizeof *NV);
    F = (int32_t *)malloc((2 * nn + 1) * 3 * sizeof *F);
    if (V == NULL || N == NULL || UV == NULL || NV == NULL || F == NULL) {
        free(vid); free(V); free(N); free(UV); free(NV); free(F);
        return -1;
    }
    for (k = 0; k < nn; k++) {
        int32_t w = vid[k];
        if (w < 0) continue;
        memcpy(&V[(size_t)w * 3], &g->pos[k * 3], 3 * sizeof *V);
        memcpy(&N[(size_t)w * 3], &g->nrm[k * 3], 3 * sizeof *N);
        UV[(size_t)w * 2] = uv_nodes[k * 2]; UV[(size_t)w * 2 + 1] = uv_nodes[k * 2 + 1];
        NV[w] = (int32_t)k;
    }
    for (k = 0; k + g->ni < nn; k++) {
        int32_t a = 0, b = 0, c = 0, d = 0;
        if (k % g->ni + 1 >= g->ni) continue;
        a = vid[k]; b = vid[k + 1]; c = vid[k + g->ni]; d = vid[k + g->ni + 1];
        if (a >= 0 && b >= 0 && c >= 0 && d >= 0) {
            cs_mesh_triangle(g, NV, a, b, d, F, &f); cs_mesh_triangle(g, NV, a, d, c, F, &f);
        } else if (a >= 0 && b >= 0 && c >= 0) cs_mesh_triangle(g, NV, a, b, c, F, &f);
        else if (b >= 0 && d >= 0 && c >= 0) cs_mesh_triangle(g, NV, b, d, c, F, &f);
        else if (a >= 0 && d >= 0 && c >= 0) cs_mesh_triangle(g, NV, a, d, c, F, &f);
        else if (a >= 0 && b >= 0 && d >= 0) cs_mesh_triangle(g, NV, a, b, d, F, &f);
    }
    free(vid);
    *verts = V; *uv = UV; *nrm = N; *faces = F; *node_of_vert = NV; *nv = v; *nf = f;
    return 0;
}

/* ----------------------------------------------------------------- frame */

int CloseupSurface_frame(const float *verts, const double *uv, const int32_t *faces, size_t nf,
                         const uint8_t *face_use, const double *axial, const float *inward, double cu, double cv,
                         CloseupFrame *fr)
{
    double pos = 0.0, neg = 0.0, gx = 0.0, gy = 0.0;
    size_t f = 0;
    int pass = 0;
    memset(fr, 0, sizeof *fr);
    /* pass 0: handedness; pass 1: the axial gradient in the (possibly reflected) layout */
    for (pass = 0; pass < 2; pass++) {
        for (f = 0; f < nf; f++) {
            size_t a = (size_t)faces[f * 3], b = (size_t)faces[f * 3 + 1], c = (size_t)faces[f * 3 + 2];
            double a1u = uv[b * 2] - uv[a * 2], a1v = uv[b * 2 + 1] - uv[a * 2 + 1];
            double a2u = uv[c * 2] - uv[a * 2], a2v = uv[c * 2 + 1] - uv[a * 2 + 1];
            double det = 0.0;
            if (face_use != NULL && !face_use[f]) continue;
            if (pass == 1 && fr->reflect) { a1u = -a1u; a2u = -a2u; }
            det = a1u * a2v - a2u * a1v;
            if (!(fabs(det) > 1e-9)) continue;
            if (pass == 0) {
                double e1[3], e2[3], cr[3], nin[3], vote = 0.0;
                int k = 0;
                for (k = 0; k < 3; k++) {
                    e1[k] = (double)verts[b * 3 + (size_t)k] - verts[a * 3 + (size_t)k];
                    e2[k] = (double)verts[c * 3 + (size_t)k] - verts[a * 3 + (size_t)k];
                    nin[k] = (double)inward[a * 3 + (size_t)k] + inward[b * 3 + (size_t)k] + inward[c * 3 + (size_t)k];
                }
                cr[0] = e1[1] * e2[2] - e1[2] * e2[1]; cr[1] = e1[2] * e2[0] - e1[0] * e2[2]; cr[2] = e1[0] * e2[1] - e1[1] * e2[0];
                vote = (cr[0] * nin[0] + cr[1] * nin[1] + cr[2] * nin[2]) * (det > 0.0 ? 1.0 : -1.0);
                if (vote > 0.0) pos += vote; else neg -= vote;
            } else {
                /* grad s = [s1 - s0, s2 - s0] inv([a1 a2]), area-weighted by |det| / 2 */
                double s1 = axial[b] - axial[a], s2 = axial[c] - axial[a];
                double g0 = (s1 * a2v - s2 * a1v) / det, g1 = (-s1 * a2u + s2 * a1u) / det;
                gx += 0.5 * fabs(det) * g0;
                gy += 0.5 * fabs(det) * g1;
            }
        }
        if (pass == 0) {
            if (pos == 0.0 && neg == 0.0) return -1;
            fr->hand_vote = pos - neg;
            fr->hand_agree = (pos >= neg ? pos : neg) / (pos + neg);
            fr->reflect = pos < neg;
        }
    }
    if (gx == 0.0 && gy == 0.0) return -1;
    fr->axial_gx = gx; fr->axial_gy = gy;
    fr->theta = atan2(gx, gy);          /* R(theta) (gx, gy) = (0, |g|) */
    fr->cu = fr->reflect ? -cu : cu;
    fr->cv = cv;
    return 0;
}

void CloseupFrame_apply(const CloseupFrame *fr, double u, double v, double *fu, double *fv)
{
    double x = (fr->reflect ? -u : u) - fr->cu, y = v - fr->cv, c = cos(fr->theta), s = sin(fr->theta);
    *fu = c * x - s * y;
    *fv = s * x + c * y;
}

int CloseupSurface_raster_attr(const double *uv, const int32_t *faces, size_t nf, const uint8_t *attr,
                               double umin, double vmin, size_t W, size_t H, uint8_t *out)
{
    long lf = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4096)
#endif
    for (lf = 0; lf < (long)nf; lf++) {
        size_t vi[3];
        double ua, va, ub, vb, uc, vc, A2;
        long x0, x1, y0, y1, x, y;
        vi[0] = (size_t)faces[lf * 3]; vi[1] = (size_t)faces[lf * 3 + 1]; vi[2] = (size_t)faces[lf * 3 + 2];
        ua = uv[vi[0] * 2] - umin; va = uv[vi[0] * 2 + 1] - vmin;
        ub = uv[vi[1] * 2] - umin; vb = uv[vi[1] * 2 + 1] - vmin;
        uc = uv[vi[2] * 2] - umin; vc = uv[vi[2] * 2 + 1] - vmin;
        A2 = (ub - ua) * (vc - va) - (vb - va) * (uc - ua);
        if (!(fabs(A2) > 1e-12)) continue;
        x0 = (long)ceil(fmin(ua, fmin(ub, uc)) - 0.5 - 1e-9); x1 = (long)floor(fmax(ua, fmax(ub, uc)) - 0.5 + 1e-9);
        y0 = (long)ceil(fmin(va, fmin(vb, vc)) - 0.5 - 1e-9); y1 = (long)floor(fmax(va, fmax(vb, vc)) - 0.5 + 1e-9);
        if (x0 < 0) x0 = 0;
        if (y0 < 0) y0 = 0;
        if (x1 >= (long)W) x1 = (long)W - 1;
        if (y1 >= (long)H) y1 = (long)H - 1;
        for (y = y0; y <= y1; y++)
            for (x = x0; x <= x1; x++) {
                double px = (double)x + 0.5, py = (double)y + 0.5;
                double l1 = ((uc - px) * (va - py) - (vc - py) * (ua - px)) / A2;
                double l2 = ((ua - px) * (vb - py) - (va - py) * (ub - px)) / A2;
                double l0 = 1.0 - l1 - l2;
                size_t best = vi[0];
                if (l0 < -1e-9 || l1 < -1e-9 || l2 < -1e-9) continue;
                if (l1 > l0 && l1 >= l2) best = vi[1];
                else if (l2 > l0 && l2 > l1) best = vi[2];
                out[(size_t)y * W + (size_t)x] = attr[best];
            }
    }
    return 0;
}

/* -------------------------------------------------------------- selftest */

/* the cylinder of radius R about the z axis: frame u = R theta, v = z; point (z, R sin th, R cos th) */
static void cs_cyl(double R, double u, double v, double off, float p[3])
{
    double th = u / R, rr = R - off;       /* off > 0: inward (toward the axis) */
    p[0] = (float)v; p[1] = (float)(rr * sin(th)); p[2] = (float)(rr * cos(th));
}

static int cs_selftest_fit(void)
{
    const double R = 120.0;
    CloseupObs o;
    CloseupGrid g;
    CloseupFitStats st;
    CloseupFitOpts opt;
    float *inward = NULL;
    size_t k = 0, bad = 0, nn = 0;
    double worst_a = 0.0, worst_b = 0.0, worst_hole = 0.0, worst_out = 0.0, worst_jump = 0.0;
    int x = 0, y = 0, s = 0;
    memset(&o, 0, sizeof o);
    if (cs_obs_reserve(&o, 80000) != 0) return 1;
    for (y = 0; y < 200; y++)
        for (x = 0; x < 300; x++) {
            double u = x + 0.5, v = y + 0.5, off = 0.0;
            size_t w = o.n;
            if (x == 150 || x == 151) continue;                          /* a 2 px trim gap */
            if ((u - 220.0) * (u - 220.0) + (v - 100.0) * (v - 100.0) < 100.0) continue;   /* a hole */
            off = x > 151 ? 1.2 : 0.0;                                   /* a 1.2 vox step */
            cs_cyl(R, u, v, off, &o.p[w * 3]);
            o.uv[w * 2] = (float)u; o.uv[w * 2 + 1] = (float)v; o.w[w] = 1.0f; o.chart[w] = x > 151 ? 2 : 1;
            o.n++;
        }
    for (y = 57; y < 63; y++)                   /* a stray patch of another wrap, 12 vox outward, over chart 1 */
        for (x = 57; x < 63; x++) {
            size_t w = o.n;
            cs_cyl(R, x + 0.5, y + 0.5, -12.0, &o.p[w * 3]);
            o.uv[w * 2] = (float)(x + 0.5); o.uv[w * 2 + 1] = (float)(y + 0.5); o.w[w] = 1.0f; o.chart[w] = 9;
            o.n++;
        }
    if (CloseupGrid_alloc(&g, -10.0, 310.0, -10.0, 210.0, 2.0) != 0) { CloseupObs_free(&o); return 1; }
    CloseupFitOpts_default(&opt);
    if (CloseupSurface_fit(&g, &o, &opt, &st) != 0) { CloseupGrid_free(&g); CloseupObs_free(&o); return 1; }
    nn = g.ni * g.nj;
    for (k = 0; k < nn; k++) {
        double u = g.u0 + g.h * (double)(k % g.ni), v = g.v0 + g.h * (double)(k / g.ni), d = 0.0;
        float t[3];
        int q = 0;
        if (g.state[k] == CLOSEUP_NODE_NONE) continue;
        if (u < 0.0 || u > 300.0 || v < 0.0 || v > 200.0) continue;
        cs_cyl(R, u, v, u > 151.0 ? 1.2 : 0.0, t);
        for (q = 0; q < 3; q++) d += ((double)g.pos[k * 3 + (size_t)q] - t[q]) * ((double)g.pos[k * 3 + (size_t)q] - t[q]);
        d = sqrt(d);
        if (u > 8.0 && u < 134.0 && v > 8.0 && v < 192.0 && !(u > 44.0 && u < 76.0 && v > 44.0 && v < 76.0))
            worst_a = fmax(worst_a, d);
        if (u > 168.0 && u < 292.0 && v > 8.0 && v < 192.0 && (u - 220.0) * (u - 220.0) + (v - 100.0) * (v - 100.0) > 400.0)
            worst_b = fmax(worst_b, d);
        if ((u - 220.0) * (u - 220.0) + (v - 100.0) * (v - 100.0) < 36.0) {
            worst_hole = fmax(worst_hole, d);
            if (g.state[k] != CLOSEUP_NODE_FILL) bad++;
        }
        if (u > 52.0 && u < 68.0 && v > 52.0 && v < 68.0) worst_out = fmax(worst_out, d);
    }
    for (k = 1; k + 1 < g.ni; k++) {          /* no tear across the step: second differences along v = 100 */
        size_t j = (size_t)((100.0 - g.v0) / g.h), n = j * g.ni + k;
        double d = 0.0;
        int q = 0;
        if (g.state[n - 1] == CLOSEUP_NODE_NONE || g.state[n] == CLOSEUP_NODE_NONE ||
            g.state[n + 1] == CLOSEUP_NODE_NONE)
            continue;
        for (q = 0; q < 3; q++) {
            double sd = (double)g.pos[(n + 1) * 3 + (size_t)q] - 2.0 * g.pos[n * 3 + (size_t)q] +
                       g.pos[(n - 1) * 3 + (size_t)q];
            d += sd * sd;
        }
        worst_jump = fmax(worst_jump, sqrt(d));
    }
    inward = (float *)calloc(nn * 3, sizeof *inward);
    if (inward != NULL) {
        for (k = 0; k < nn; k++) {
            double l = 0.0;
            if (g.state[k] == CLOSEUP_NODE_NONE) continue;
            l = sqrt((double)g.pos[k * 3 + 1] * g.pos[k * 3 + 1] + (double)g.pos[k * 3 + 2] * g.pos[k * 3 + 2]);
            inward[k * 3 + 1] = (float)(-g.pos[k * 3 + 1] / l); inward[k * 3 + 2] = (float)(-g.pos[k * 3 + 2] / l);
        }
        s = CloseupSurface_orient(&g, inward, NULL);
        free(inward);
    }
    fprintf(stderr, "  closeup fit: data %zu fill %zu (rejected %zu); error chart1 %.3f chart2 %.3f hole %.3f "
            "other-wrap patch %.3f vox; max 2nd difference along a row %.3f vox; resid p50 %.3f p99 %.3f; "
            "orient %d; %.2fs\n", st.data_nodes, st.fill_nodes, st.fill_rejected, worst_a, worst_b, worst_hole,
            worst_out, worst_jump, st.resid_p50, st.resid_p99, s, st.seconds);
    CloseupGrid_free(&g);
    CloseupObs_free(&o);
    if (worst_a > 0.05 || worst_b > 0.05 || worst_hole > 0.25 || worst_out > 0.25 || worst_jump > 0.5 || bad || s != 1) {
        fprintf(stderr, "FAIL closeup fit\n");
        return 1;
    }
    return 0;
}

static int cs_selftest_flatten(void)
{
    const double R = 90.0;
    CloseupObs o;
    CloseupGrid g;
    CloseupFlatStats fs;
    double *uvn = NULL, arc = 0.0, du = 0.0, dv = 0.0;
    size_t a = 0, b = 0, j = 0;
    int x = 0, y = 0, rc = 1;
    memset(&o, 0, sizeof o);
    if (cs_obs_reserve(&o, 70000) != 0) return 1;
    /* the frame stretches the right half by 10%: frame u = x for x < 150, 150 + 1.1 (x - 150) beyond */
    for (y = 0; y < 160; y++)
        for (x = 0; x < 320; x++) {
            double uf = x + 0.5, ut = uf < 150.0 ? uf : 150.0 + (uf - 150.0) / 1.1;
            size_t w = o.n;
            cs_cyl(R, ut, y + 0.5, 0.0, &o.p[w * 3]);
            o.uv[w * 2] = (float)uf; o.uv[w * 2 + 1] = (float)(y + 0.5); o.w[w] = 1.0f; o.chart[w] = 1;
            o.n++;
        }
    if (CloseupGrid_alloc(&g, 0.0, 320.0, 0.0, 160.0, 2.0) != 0) { CloseupObs_free(&o); return 1; }
    if (CloseupSurface_fit(&g, &o, NULL, NULL) != 0) goto out;
    uvn = (double *)malloc(g.ni * g.nj * 2 * sizeof *uvn);
    if (uvn == NULL || CloseupSurface_flatten(&g, 8.0, 80, 1e-9, uvn, &fs) != 0) goto out;
    j = (size_t)(80.0 / g.h);
    a = j * g.ni + (size_t)(20.0 / g.h);
    b = j * g.ni + (size_t)(300.0 / g.h);
    /* true arc length between frame u 20 and 300: 130 + 150 / 1.1 */
    arc = 130.0 + 150.0 / 1.1;
    du = uvn[b * 2] - uvn[a * 2]; dv = uvn[b * 2 + 1] - uvn[a * 2 + 1];
    fprintf(stderr, "  closeup flatten: %zu coarse nodes, %zu tris, %zu iters; |log s| p99 %.4f -> %.4f; "
            "arc %.2f (true %.2f); area ratio %.4f; %.2fs\n", fs.coarse_nodes, fs.coarse_tris, fs.iterations,
            fs.logs_p99_before, fs.logs_p99_after, sqrt(du * du + dv * dv), arc, fs.area_ratio, fs.seconds);
    {
        size_t q = 0;
        for (q = 1; q < fs.iterations && q < 64; q++)
            if (fs.energy_hist[q] > fs.energy_hist[q - 1] * (1.0 + 1e-6) + 1e-9) {
                fprintf(stderr, "FAIL closeup flatten: energy climbs at iteration %zu (%.6g -> %.6g)\n", q,
                        fs.energy_hist[q - 1], fs.energy_hist[q]);
                goto out;
            }
    }
    if (fabs(sqrt(du * du + dv * dv) / arc - 1.0) > 0.01 || fs.logs_p99_after > 0.02 || fs.logs_p99_before < 0.05) {
        fprintf(stderr, "FAIL closeup flatten\n");
        goto out;
    }
    rc = 0;
out:
    free(uvn);
    CloseupGrid_free(&g);
    CloseupObs_free(&o);
    return rc;
}

/* CloseupSurface_flatten_slim on the stretched-frame cylinder of cs_selftest_flatten with a hole and a slit (the
 * upper right part hangs on the lower part like a flap): the arc length is recovered on a row clear of both, SLIM is
 * accepted, and the fine map is injective -- no fine triangle flipped, no UV point covered twice. */
static int cs_selftest_flatten_slim(void)
{
    const double R = 90.0;
    CloseupObs o;
    CloseupGrid g;
    CloseupFlatStats fs;
    double *uvn = NULL, arc = 0.0, du = 0.0, dv = 0.0, lo[2] = { 1e300, 1e300 }, hi[2] = { -1e300, -1e300 };
    size_t a = 0, b = 0, i = 0, j = 0, flips = 0, twice = 0, nw = 0, nh = 0;
    uint16_t *cover = NULL;
    int x = 0, y = 0, rc = 1;
    memset(&o, 0, sizeof o);
    memset(&g, 0, sizeof g);
    memset(&fs, 0, sizeof fs);
    if (cs_obs_reserve(&o, 70000) != 0) return 1;
    for (y = 0; y < 160; y++)
        for (x = 0; x < 320; x++) {
            double uf = x + 0.5, ut = uf < 150.0 ? uf : 150.0 + (uf - 150.0) / 1.1;
            size_t w = o.n;
            cs_cyl(R, ut, y + 0.5, 0.0, &o.p[w * 3]);
            o.uv[w * 2] = (float)uf; o.uv[w * 2 + 1] = (float)(y + 0.5); o.w[w] = 1.0f; o.chart[w] = 1;
            o.n++;
        }
    if (CloseupGrid_alloc(&g, 0.0, 320.0, 0.0, 160.0, 2.0) != 0) { CloseupObs_free(&o); return 1; }
    if (CloseupSurface_fit(&g, &o, NULL, NULL) != 0) goto out;
    for (j = 0; j < g.nj; j++)
        for (i = 0; i < g.ni; i++) {
            const double u = g.u0 + g.h * (double)i, v = g.v0 + g.h * (double)j;
            if ((u > 60.0 && u < 100.0 && v > 40.0 && v < 70.0) || (u > 200.0 && u < 206.0 && v < 100.0))
                g.state[j * g.ni + i] = CLOSEUP_NODE_NONE;          /* a hole, and a slit down to v 100 */
        }
    uvn = (double *)malloc(g.ni * g.nj * 2 * sizeof *uvn);
    /* Exercise the ordinary seed entry point on a lacy flap: an ARAP route
     * here must not bypass the fine-triangle and overlap regression below. */
    if (uvn == NULL || CloseupSurface_flatten(&g, 8.0, 40, 1e-9, uvn, &fs) != 0) goto out;
    j = (size_t)(130.0 / g.h);
    a = j * g.ni + (size_t)(20.0 / g.h);
    b = j * g.ni + (size_t)(300.0 / g.h);
    arc = 130.0 + 150.0 / 1.1;
    du = uvn[b * 2] - uvn[a * 2]; dv = uvn[b * 2 + 1] - uvn[a * 2 + 1];
    /* injectivity: every fine cell's two triangles positive, and no 1 px sample covered twice */
    for (j = 0; j + 1 < g.nj; j++)
        for (i = 0; i + 1 < g.ni; i++) {
            const size_t c4[4] = { j * g.ni + i, j * g.ni + i + 1, (j + 1) * g.ni + i + 1, (j + 1) * g.ni + i };
            int q = 0, all = 1;
            for (q = 0; q < 4; q++) all &= g.state[c4[q]] != CLOSEUP_NODE_NONE && isfinite(uvn[c4[q] * 2]);
            if (!all) continue;
            for (q = 0; q < 4; q++) {
                lo[0] = fmin(lo[0], uvn[c4[q] * 2]); lo[1] = fmin(lo[1], uvn[c4[q] * 2 + 1]);
                hi[0] = fmax(hi[0], uvn[c4[q] * 2]); hi[1] = fmax(hi[1], uvn[c4[q] * 2 + 1]);
            }
        }
    if (!(hi[0] > lo[0])) goto out;
    nw = (size_t)(hi[0] - lo[0]) + 3; nh = (size_t)(hi[1] - lo[1]) + 3;
    cover = (uint16_t *)calloc(nw * nh, sizeof *cover);
    if (cover == NULL) goto out;
    for (j = 0; j + 1 < g.nj; j++)
        for (i = 0; i + 1 < g.ni; i++) {
            const size_t c4[4] = { j * g.ni + i, j * g.ni + i + 1, (j + 1) * g.ni + i + 1, (j + 1) * g.ni + i };
            int tris[2][3], q = 0, all = 1, tt = 0;
            tris[0][0] = 0; tris[0][1] = 1; tris[0][2] = 2; tris[1][0] = 0; tris[1][1] = 2; tris[1][2] = 3;
            for (q = 0; q < 4; q++) all &= g.state[c4[q]] != CLOSEUP_NODE_NONE && isfinite(uvn[c4[q] * 2]);
            if (!all) continue;
            for (tt = 0; tt < 2; tt++) {
                const double *P = &uvn[c4[tris[tt][0]] * 2], *Q = &uvn[c4[tris[tt][1]] * 2], *Z = &uvn[c4[tris[tt][2]] * 2];
                const double ar = (Q[0] - P[0]) * (Z[1] - P[1]) - (Q[1] - P[1]) * (Z[0] - P[0]);
                long px = 0, py = 0, y0 = 0, y1 = 0, x0 = 0, x1 = 0;
                if (!(ar > 0.0)) { flips++; continue; }
                y0 = (long)floor(fmin(P[1], fmin(Q[1], Z[1])) - lo[1]); y1 = (long)ceil(fmax(P[1], fmax(Q[1], Z[1])) - lo[1]);
                x0 = (long)floor(fmin(P[0], fmin(Q[0], Z[0])) - lo[0]); x1 = (long)ceil(fmax(P[0], fmax(Q[0], Z[0])) - lo[0]);
                for (py = y0; py <= y1; py++)
                    for (px = x0; px <= x1; px++) {
                        const double sx = lo[0] + (double)px + 0.5, sy = lo[1] + (double)py + 0.5;
                        const double w0 = ((Q[0] - sx) * (Z[1] - sy) - (Q[1] - sy) * (Z[0] - sx)) / ar;
                        const double w1 = ((Z[0] - sx) * (P[1] - sy) - (Z[1] - sy) * (P[0] - sx)) / ar;
                        if (px < 0 || py < 0 || (size_t)px >= nw || (size_t)py >= nh) continue;
                        if (w0 > 1e-9 && w1 > 1e-9 && 1.0 - w0 - w1 > 1e-9) {
                            if (cover[(size_t)py * nw + (size_t)px] == 1) twice++;
                            cover[(size_t)py * nw + (size_t)px]++;
                        }
                    }
            }
        }
    fprintf(stderr, "  closeup flatten SLIM: %zu coarse nodes, %zu tris, %zu iters; |log s| p99 %.4f -> %.4f; arc %.2f (true "
            "%.2f); %zu flipped, %zu px covered twice; %.2fs\n", fs.coarse_nodes, fs.coarse_tris, fs.iterations,
            fs.logs_p99_before, fs.logs_p99_after, sqrt(du * du + dv * dv), arc, flips, twice, fs.seconds);
    if (fs.iterations == 0 || fabs(sqrt(du * du + dv * dv) / arc - 1.0) > 0.01 || fs.logs_p99_after > 0.03 ||
        flips != 0 || twice != 0) {
        fprintf(stderr, "FAIL closeup flatten SLIM\n");
        goto out;
    }
    rc = 0;
out:
    free(uvn); free(cover);
    CloseupGrid_free(&g);
    CloseupObs_free(&o);
    return rc;
}

/* CloseupGrid_save / _load: a grid with a valid block, back node for node; a truncated file refused. */
/* The worst |log| of fine faces' UV over 3-D edge lengths among the faces with every node at frame u in [ulo, uhi]
 * (and v in (vlo, vhi)); -1 when a face there is flipped. */
static double cs_fine_stretch(const CloseupGrid *g, const double *uvn, double ulo, double uhi, double vlo, double vhi)
{
    size_t i = 0, j = 0;
    double worst = 0.0;
    for (j = 0; j + 1 < g->nj; j++)
        for (i = 0; i + 1 < g->ni; i++) {
            const size_t c4[4] = { j * g->ni + i, j * g->ni + i + 1, (j + 1) * g->ni + i + 1, (j + 1) * g->ni + i };
            const double u = g->u0 + g->h * (double)i, v = g->v0 + g->h * (double)j;
            int q = 0, all = 1;
            if (u < ulo || u + g->h > uhi || v <= vlo || v + g->h >= vhi) continue;
            for (q = 0; q < 4; q++) all &= g->state[c4[q]] != CLOSEUP_NODE_NONE && isfinite(uvn[c4[q] * 2]);
            if (!all) continue;
            {
                const double *P = &uvn[c4[0] * 2], *Q = &uvn[c4[1] * 2], *Z = &uvn[c4[2] * 2];
                if ((Q[0] - P[0]) * (Z[1] - P[1]) - (Q[1] - P[1]) * (Z[0] - P[0]) <= 0.0) return -1.0;
            }
            for (q = 0; q < 4; q++) {
                const size_t s = c4[q], t = c4[(q + 1) % 4];
                const double l3 = sqrt(((double)g->pos[s * 3] - g->pos[t * 3]) * ((double)g->pos[s * 3] - g->pos[t * 3]) +
                                       ((double)g->pos[s * 3 + 1] - g->pos[t * 3 + 1]) * ((double)g->pos[s * 3 + 1] - g->pos[t * 3 + 1]) +
                                       ((double)g->pos[s * 3 + 2] - g->pos[t * 3 + 2]) * ((double)g->pos[s * 3 + 2] - g->pos[t * 3 + 2]));
                const double lu = hypot(uvn[s * 2] - uvn[t * 2], uvn[s * 2 + 1] - uvn[t * 2 + 1]);
                if (l3 > 1e-9 && lu > 1e-12) worst = fmax(worst, fabs(log(lu / l3)));
            }
        }
    return worst;
}

/* The fringe: a plane whose frame is true for u < 290 and squeezes the sheet 2.5x beyond (1 frame px = 2.5 vox), the
 * sheet ending at u 306 -- one fine step into the coarse cell [304, 312] that its lattice node 312 leaves empty.  With
 * the fringe continuation that cell is flattened as surface and the fine faces in it come out nearly true; without it
 * they keep the frame's squeeze (scaffold).  No face flips either way. */
static int cs_selftest_flatten_fringe(void)
{
    CloseupGrid g;
    CloseupFlatStats f0s, f1s;
    double *u0n = NULL, *u1n = NULL, s0 = 0.0, s1 = 0.0;
    size_t i = 0, j = 0;
    int rc = 1;
    memset(&f0s, 0, sizeof f0s);
    memset(&f1s, 0, sizeof f1s);
    if (CloseupGrid_alloc(&g, 0.0, 320.0, 0.0, 96.0, 2.0) != 0) return 1;
    for (j = 0; j < g.nj; j++)
        for (i = 0; i < g.ni; i++) {
            const size_t n = j * g.ni + i;
            const double u = g.u0 + g.h * (double)i, v = g.v0 + g.h * (double)j;
            const double t = u < 290.0 ? u : 290.0 + 2.5 * (u - 290.0);
            if (u > 306.0) continue;
            g.state[n] = CLOSEUP_NODE_DATA; g.component[n] = 0; g.chart[n] = 0;
            g.pos[n * 3] = (float)v; g.pos[n * 3 + 1] = 0.0f; g.pos[n * 3 + 2] = (float)t;
            g.nrm[n * 3 + 1] = 1.0f;
        }
    u0n = (double *)malloc(g.ni * g.nj * 2 * sizeof *u0n);
    u1n = (double *)malloc(g.ni * g.nj * 2 * sizeof *u1n);
    if (u0n == NULL || u1n == NULL) goto out;
    if (cs_flatten_slim(&g, 8.0, 40, 0, u0n, &f0s) != 0 || cs_flatten_slim(&g, 8.0, 40, 1, u1n, &f1s) != 0) goto out;
    s0 = cs_fine_stretch(&g, u0n, 304.0, 306.0, 8.0, 88.0);
    s1 = cs_fine_stretch(&g, u1n, 304.0, 306.0, 8.0, 88.0);
    fprintf(stderr, "  closeup flatten SLIM fringe: %zu lattice nodes continued; the partly covered cell's fine faces |log s| "
            "%.3f without, %.3f with the continuation; interior |log s| p99 %.4f / %.4f\n", f1s.coarse_fringe, s0, s1,
            f0s.logs_p99_after, f1s.logs_p99_after);
    if (f1s.coarse_fringe == 0 || s0 < 0.0 || s1 < 0.0 || !(s1 < 0.25) || !(s0 > 0.5) || f1s.iterations == 0) {
        fprintf(stderr, "FAIL closeup flatten SLIM fringe\n");
        goto out;
    }
    rc = 0;
out:
    free(u0n); free(u1n);
    CloseupGrid_free(&g);
    return rc;
}

/* CloseupSurface_dedup_fill: a DATA sheet (frame u 0..100) and three FILL patches -- B far in the frame (u 200..240)
 * lying ON the sheet in 3-D (dropped), C beside the sheet in the frame continuing it (kept), D far in the frame 5 vox off
 * the sheet (kept) -- and a second DATA patch E far in the frame on the sheet (a DATA pair: counted, kept). */
static int cs_selftest_dedup(void)
{
    CloseupGrid g;
    CloseupDedupStats ds;
    size_t i = 0, j = 0, kept_b = 0, kept_c = 0, kept_d = 0, kept_e = 0;
    int rc = 1;
    memset(&ds, 0, sizeof ds);
    if (CloseupGrid_alloc(&g, 0.0, 400.0, 0.0, 40.0, 2.0) != 0) return 1;
    for (j = 0; j < g.nj; j++)
        for (i = 0; i < g.ni; i++) {
            const size_t n = j * g.ni + i;
            const double u = g.u0 + g.h * (double)i, v = g.v0 + g.h * (double)j;
            double x = 0.0, z = 0.0;
            int st = CLOSEUP_NODE_NONE;
            if (u <= 100.0) { st = CLOSEUP_NODE_DATA; x = u; }                                     /* the sheet */
            else if (u <= 120.0) { st = CLOSEUP_NODE_FILL; x = u; }                                /* C: its continuation */
            else if (u >= 200.0 && u <= 240.0) { st = CLOSEUP_NODE_FILL; x = u - 150.0; }          /* B: on the sheet */
            else if (u >= 260.0 && u <= 300.0) { st = CLOSEUP_NODE_FILL; x = u - 210.0; z = 5.0; } /* D: 5 vox off */
            else if (u >= 340.0 && u <= 360.0) { st = CLOSEUP_NODE_DATA; x = u - 330.0; }          /* E: data on data */
            if (st == CLOSEUP_NODE_NONE) continue;
            g.state[n] = (uint8_t)st; g.component[n] = 0; g.chart[n] = 0;
            g.pos[n * 3] = (float)z; g.pos[n * 3 + 1] = (float)v; g.pos[n * 3 + 2] = (float)x;
            g.nrm[n * 3] = 1.0f;
            g.dist[n] = st == CLOSEUP_NODE_DATA ? 0.0f : 10.0f;
        }
    if (CloseupSurface_dedup_fill(&g, 2.0, 30.0, 0.8, &ds) != 0) goto out;
    for (j = 0; j < g.nj; j++)
        for (i = 0; i < g.ni; i++) {
            const double u = g.u0 + g.h * (double)i;
            const int on = g.state[j * g.ni + i] != CLOSEUP_NODE_NONE;
            if (u >= 200.0 && u <= 240.0) kept_b += (size_t)on;
            else if (u > 100.0 && u <= 120.0) kept_c += (size_t)on;
            else if (u >= 260.0 && u <= 300.0) kept_d += (size_t)on;
            else if (u >= 340.0 && u <= 360.0) kept_e += (size_t)on;
        }
    fprintf(stderr, "  closeup dedup: %zu pairs (%zu DATA-DATA), %zu FILL dropped; kept B %zu, C %zu, D %zu, E %zu\n", ds.pairs,
            ds.data_pairs, ds.dropped_fill, kept_b, kept_c, kept_d, kept_e);
    if (kept_b != 0 || kept_c != 10 * g.nj || kept_d != 21 * g.nj || kept_e != 11 * g.nj || ds.dropped_fill != 21 * g.nj ||
        ds.data_pairs == 0) {
        fprintf(stderr, "FAIL closeup dedup\n");
        goto out;
    }
    rc = 0;
out:
    CloseupGrid_free(&g);
    return rc;
}

/* CloseupSurface_trim_apron: a DATA frame (a ring of DATA 20 px wide round a 40 x 40 px hole) beside a second DATA
 * block across a 30 px gap, every other node within the frame FILL: the hole (FILL surrounded by DATA) and the gap
 * (DATA on both sides along the row) stay, the margin outside the DATA (DATA on one side only) goes, and the
 * diagonal case: a FILL node off a corner whose diagonal meets DATA on one side only goes too. */
static int cs_selftest_apron(void)
{
    CloseupGrid g;
    CloseupApronStats as;
    size_t i = 0, j = 0, hole = 0, gap = 0, margin = 0, data0 = 0, data1 = 0;
    int rc = 1;
    memset(&as, 0, sizeof as);
    if (CloseupGrid_alloc(&g, 0.0, 200.0, 0.0, 120.0, 2.0) != 0) return 1;
    for (j = 0; j < g.nj; j++)
        for (i = 0; i < g.ni; i++) {
            const size_t n = j * g.ni + i;
            const double u = g.u0 + g.h * (double)i, v = g.v0 + g.h * (double)j;
            const int in_a = u >= 20.0 && u <= 100.0 && v >= 20.0 && v <= 100.0;                 /* the frame's box */
            const int in_hole = u > 40.0 && u < 80.0 && v > 40.0 && v < 80.0;
            const int in_b = u >= 130.0 && u <= 170.0 && v >= 20.0 && v <= 100.0;               /* the second block */
            int st = CLOSEUP_NODE_FILL;
            if ((in_a && !in_hole) || in_b) st = CLOSEUP_NODE_DATA;
            g.state[n] = (uint8_t)st; g.component[n] = 0; g.chart[n] = 0;
            g.pos[n * 3] = 0.0f; g.pos[n * 3 + 1] = (float)v; g.pos[n * 3 + 2] = (float)u;
            data0 += st == CLOSEUP_NODE_DATA;
        }
    if (CloseupSurface_trim_apron(&g, 30.0, &as) != 0) goto out;
    for (j = 0; j < g.nj; j++)
        for (i = 0; i < g.ni; i++) {
            const size_t n = j * g.ni + i;
            const double u = g.u0 + g.h * (double)i, v = g.v0 + g.h * (double)j;
            const int on = g.state[n] == CLOSEUP_NODE_FILL;
            data1 += g.state[n] == CLOSEUP_NODE_DATA;
            if (u > 40.0 && u < 80.0 && v > 40.0 && v < 80.0) hole += (size_t)!on;           /* must all stay */
            else if (u > 100.0 && u < 130.0 && v >= 20.0 && v <= 100.0) gap += (size_t)!on;  /* must all stay */
            else if (u < 20.0 || v < 20.0 || v > 100.0 || u > 170.0) margin += (size_t)on;   /* must all go */
        }
    fprintf(stderr, "  closeup apron: %zu of %zu FILL dropped; hole lost %zu, gap lost %zu, margin kept %zu, DATA %zu -> %zu\n",
            as.dropped, as.fill_before, hole, gap, margin, data0, data1);
    if (hole != 0 || gap != 0 || margin != 0 || data0 != data1 || as.dropped == 0) {
        fprintf(stderr, "FAIL closeup apron\n");
        goto out;
    }
    rc = 0;
out:
    CloseupGrid_free(&g);
    return rc;
}

static int cs_selftest_grid_file(void)
{
    char dir[512], path[700];
    CloseupGrid g, h;
    int32_t meta[4] = { 7, -1, 1, 0 }, back[4] = { 0, 0, 0, 0 };
    size_t k = 0, nm = 0, nn = 0;
    int fails = 0;
    FILE *fp = NULL;
    memset(&g, 0, sizeof g);
    memset(&h, 0, sizeof h);
    if (ves_selftest_dir(dir, sizeof dir, "closeup_grid_file") != 0) return 1;
    snprintf(path, sizeof path, "%s/g.cgrid", dir);
    if (CloseupGrid_alloc(&g, -3.0, 11.0, 2.0, 9.0, 1.5) != 0) return 1;
    nn = g.ni * g.nj;
    for (k = 0; k < nn; k++) {
        if (k % 3 == 1) continue;
        g.state[k] = (uint8_t)(k % 2 ? CLOSEUP_NODE_DATA : CLOSEUP_NODE_FILL);
        g.chart[k] = (int32_t)k - 5; g.component[k] = (int32_t)(k % 4);
        g.pos[k * 3] = (float)k; g.pos[k * 3 + 1] = 0.5f * (float)k; g.pos[k * 3 + 2] = -(float)k;
        g.nrm[k * 3] = 0.0f; g.nrm[k * 3 + 1] = 0.6f; g.nrm[k * 3 + 2] = 0.8f;
    }
    fails += CloseupGrid_save(path, &g, meta, 3) != 0;
    fails += CloseupGrid_load(path, &h, back, 4, &nm) != 0;
    if (!fails) {
        fails += h.ni != g.ni || h.nj != g.nj || fabs(h.u0 - g.u0) > 0.0 || fabs(h.v0 - g.v0) > 0.0 || fabs(h.h - g.h) > 0.0;
        fails += nm != 3 || back[0] != 7 || back[1] != -1 || back[2] != 1 || back[3] != 0;
        for (k = 0; k < nn && !fails; k++) {
            if (h.state[k] != g.state[k]) { fails++; break; }
            if (g.state[k] == CLOSEUP_NODE_NONE) { fails += isfinite(h.pos[k * 3]) || h.component[k] != -1; continue; }
            fails += h.chart[k] != g.chart[k] || h.component[k] != g.component[k] ||
                     memcmp(&h.pos[k * 3], &g.pos[k * 3], 3 * sizeof *g.pos) != 0 ||
                     memcmp(&h.nrm[k * 3], &g.nrm[k * 3], 3 * sizeof *g.nrm) != 0;
        }
    }
    CloseupGrid_free(&h);
    fp = fopen(path, "r+b");                       /* truncated: the node block cut short */
    if (fp == NULL) fails++;
    else {
        fclose(fp);
        fp = fopen(path, "rb");
        if (fp != NULL) {
            char buf[256];
            size_t got = fread(buf, 1, sizeof buf, fp);
            fclose(fp);
            fp = fopen(path, "wb");
            if (fp != NULL) { fwrite(buf, 1, got / 2, fp); fclose(fp); }
        }
        fails += CloseupGrid_load(path, &h, back, 4, &nm) == 0;
        CloseupGrid_free(&h);
    }
    CloseupGrid_free(&g);
    fprintf(stderr, "  closeup grid file: %s\n", fails ? "FAIL" : "ok");
    if (fails == 0) (void)ves_selftest_clean(dir);
    return fails ? 1 : 0;
}

static int cs_selftest_raster(void)
{
    /* a 4 x 3 px rectangle split in two triangles: 12 centres, each exactly once */
    const float V[12] = { 0, 0, 0,  0, 0, 4,  0, 3, 4,  0, 3, 0 };
    const double UV[8] = { 0.0, 0.0,  4.0, 0.0,  4.0, 3.0,  0.0, 3.0 };
    const int32_t F[6] = { 0, 1, 2,  0, 2, 3 };
    CloseupObs o;
    CloseupObsGate gt;
    size_t used = 0, gated = 0;
    int rc = 0;
    memset(&o, 0, sizeof o);
    CloseupObsGate_default(&gt);
    if (CloseupObs_add_faces(&o, V, UV, F, 2, NULL, NULL, &gt, &used, &gated) != 0) return 1;
    /* centres on the diagonal belong to both triangles (edge tolerance): 12 + the diagonal hits */
    if (o.n < 12 || o.n > 14 || used != 2 || gated != 0) rc = 1;
    fprintf(stderr, "  closeup raster: %zu observations from 2 faces (12 centres + diagonal ties)\n", o.n);
    CloseupObs_free(&o);
    {
        /* Accepted, rejected, accepted: the rejected face must retain the
         * current prefix. Poison spare capacity so a counted but unwritten
         * observation fails deterministically, including on an append. */
        float broken[12]; memcpy(broken, V, sizeof broken); broken[11] = 100.0f;
        const int32_t faces[9] = {0,1,2, 0,2,3, 0,1,2};
        if (cs_obs_reserve(&o, 128) != 0) return 1;
        for (size_t k = 0; k < o.cap * 3; k++) o.p[k] = -12345.0f;
        for (int pass = 0; pass < 2; pass++) {
            size_t base = o.n;
            if (CloseupObs_add_faces(&o, broken, UV, faces, 3, NULL, NULL, &gt, &used, &gated) ||
                used != 2 || gated != 1 || o.n <= base) rc = 1;
            for (size_t k = base; k < o.n; k++)
                if (o.p[3*k] != 0 || !(o.p[3*k+1] >= 0 && o.p[3*k+1] <= 3) ||
                    !(o.p[3*k+2] >= 0 && o.p[3*k+2] <= 4)) rc = 1;
        }
        CloseupObs_free(&o);
    }
    {
        /* Edge ratios alone miss the 100x stretch across a thin triangle. */
        float thin[9] = {0,0,0, 0,0,10, 0,.01f,5};
        double shear[6] = {0,0, 10,0, 5,1}; int32_t face[3] = {0,1,2};
        long x0=0,x1=0,y0=0,y1=0;
        if (cs_face_box(thin,shear,face,0,&gt,&x0,&x1,&y0,&y1) != 0) rc = 1;
        shear[5] = .0001;
        if (cs_face_box(thin,shear,face,0,&gt,&x0,&x1,&y0,&y1) != 0) rc = 1;
        shear[5] = (double)thin[7];
        if (cs_face_box(thin,shear,face,0,&gt,&x0,&x1,&y0,&y1) != 1) rc = 1;
        thin[0] = NAN;
        if (cs_face_box(thin,shear,face,0,&gt,&x0,&x1,&y0,&y1) != 0) rc = 1;
    }
    if (fabs(CloseupSurface_rot_angle(cos(0.3) * 1.2, -sin(0.3) * 0.9, sin(0.3) * 1.2, cos(0.3) * 0.9) - 0.3) > 0.05)
        rc = 1;
    if (rc) fprintf(stderr, "FAIL closeup raster/rotation\n");
    return rc;
}

static int cs_selftest_bridge(int complete)
{
    /* piece A: plane z = 0 on u in [0, 100); piece B: plane z = 30 on u in [132, 232) (another wrap laid out 32 px
     * to the right).  No fill node may take geometry from both: every fill node stays within 1 vox of one plane. */
    CloseupObs o;
    CloseupGrid g;
    CloseupFitOpts opt;
    size_t k = 0, bad = 0, nfill = 0, nv = 0, nf = 0;
    float *verts = NULL, *nrm = NULL;
    double *uvn = NULL, *uv = NULL;
    int32_t *faces = NULL, *nodes = NULL;
    int x = 0, y = 0;
    memset(&o, 0, sizeof o);
    memset(&g, 0, sizeof g);
    if (complete) CloseupFitOpts_complete(&opt); else CloseupFitOpts_default(&opt);
    if (cs_obs_reserve(&o, 2 * 100 * 60) != 0) return 1;
    for (y = 0; y < 60; y++)
        for (x = 0; x < 232; x++) {
            size_t w = o.n;
            if (x >= 100 && x < 132) continue;
            o.uv[w * 2] = (float)(x + 0.5); o.uv[w * 2 + 1] = (float)(y + 0.5);
            o.p[w * 3] = x < 100 ? 0.0f : 30.0f; o.p[w * 3 + 1] = (float)(y + 0.5); o.p[w * 3 + 2] = (float)(x + 0.5);
            o.w[w] = 1.0f; o.chart[w] = x < 100 ? 1 : 2;
            o.n++;
        }
    if (CloseupGrid_alloc(&g, -10.0, 242.0, -10.0, 70.0, 2.0) != 0 || CloseupSurface_fit(&g, &o, &opt, NULL) != 0) {
        CloseupObs_free(&o); CloseupGrid_free(&g); return 1;
    }
    for (k = 0; k < g.ni * g.nj; k++) {
        double z = 0.0;
        if (g.state[k] != CLOSEUP_NODE_FILL) continue;
        nfill++;
        z = g.pos[k * 3];
        if (fabs(z) > 1.0 && fabs(z - 30.0) > 1.0) bad++;
    }
    fprintf(stderr, "  closeup bridge (complete %d): %zu fill nodes, %zu mix the two pieces\n", complete, nfill, bad);
    uvn = (double *)malloc(g.ni * g.nj * 2 * sizeof *uvn);
    if (uvn == NULL) bad++;
    else {
        for (k = 0; k < g.ni * g.nj; k++) {
            uvn[k * 2] = g.u0 + (double)(k % g.ni) * g.h;
            uvn[k * 2 + 1] = g.v0 + (double)(k / g.ni) * g.h;
        }
        if (CloseupSurface_mesh(&g, uvn, &verts, &uv, &nrm, &faces, &nodes, &nv, &nf) != 0 || nf == 0) bad++;
        else for (k = 0; k < nf; k++) {
            double za = verts[(size_t)faces[k * 3] * 3], zb = verts[(size_t)faces[k * 3 + 1] * 3];
            double zc = verts[(size_t)faces[k * 3 + 2] * 3];
            if (fabs(za - zb) > 1.0 || fabs(za - zc) > 1.0) bad++;
        }
        fprintf(stderr, "  closeup bridge mesh (complete %d): %zu triangles, %zu mixing/errors\n", complete, nf, bad);
    }
    free(uvn); free(verts); free(uv); free(nrm); free(faces); free(nodes);
    CloseupObs_free(&o);
    CloseupGrid_free(&g);
    if (bad) { fprintf(stderr, "FAIL closeup bridge\n"); return 1; }
    return 0;
}

static int cs_selftest_bridge_gated(void)
{
    /* The gated bridge on the bridge fixture (pieces A and B 32 px apart in the frame): (0) one plane, one wrap by
     * the index -> one component, one bridge; (1) one plane, different wraps -> apart; (2) one wrap by the index but
     * B 30 vox off in z (two layers mislabelled) -> apart by the continuity gate. */
    static const int32_t same[3] = { INT32_MIN, 5, 5 }, diff[3] = { INT32_MIN, 5, 6 };
    int bad = 0, cas = 0;
    for (cas = 0; cas < 3; cas++) {
        CloseupObs o;
        CloseupGrid g;
        CloseupFitOpts opt;
        CloseupFitStats st;
        size_t k = 0, comps = 0;
        int32_t c0 = -1, c1 = -1;
        int x = 0, y = 0;
        const double zb = cas == 2 ? 30.0 : 0.0;
        memset(&o, 0, sizeof o);
        memset(&g, 0, sizeof g);
        CloseupFitOpts_default(&opt);
        opt.bridge = 1; opt.chart_wrap = cas == 1 ? diff : same; opt.nchart_wrap = 3;
        if (cs_obs_reserve(&o, 2 * 100 * 60) != 0) return 1;
        for (y = 0; y < 60; y++)
            for (x = 0; x < 232; x++) {
                size_t w = o.n;
                if (x >= 100 && x < 132) continue;
                o.uv[w * 2] = (float)(x + 0.5); o.uv[w * 2 + 1] = (float)(y + 0.5);
                o.p[w * 3] = x < 100 ? 0.0f : (float)zb; o.p[w * 3 + 1] = (float)(y + 0.5); o.p[w * 3 + 2] = (float)(x + 0.5);
                o.w[w] = 1.0f; o.chart[w] = x < 100 ? 1 : 2;
                o.n++;
            }
        if (CloseupGrid_alloc(&g, -10.0, 242.0, -10.0, 70.0, 2.0) != 0 || CloseupSurface_fit(&g, &o, &opt, &st) != 0) {
            CloseupObs_free(&o); CloseupGrid_free(&g); return 1;
        }
        for (k = 0; k < g.ni * g.nj; k++) {
            double u = g.u0 + (double)(k % g.ni) * g.h;
            if (g.state[k] != CLOSEUP_NODE_DATA) continue;
            if (u < 90.0) c0 = g.component[k]; else if (u > 142.0) c1 = g.component[k];
        }
        comps = (c0 >= 0 && c1 >= 0 && c0 == c1) ? 1 : 2;
        if (cas == 0 && (comps != 1 || st.bridges != 1)) bad++;
        if (cas == 1 && (comps != 2 || st.bridges != 0)) bad++;
        if (cas == 2 && (comps != 2 || st.bridges != 0)) bad++;
        fprintf(stderr, "  closeup gated bridge (case %d): blobs %zu, with a wrap %zu, bridges %zu, components of the two pieces %zu\n",
                cas, st.blobs, st.blobs_with_wrap, st.bridges, comps);
        CloseupObs_free(&o);
        CloseupGrid_free(&g);
    }
    if (bad) { fprintf(stderr, "FAIL closeup gated bridge\n"); return 1; }
    return 0;
}

static int cs_selftest_tear(void)
{
    /* Two planes laid out side by side in the frame with no gap (u in [0, 100) and [100, 200)): (a) one sheet
     * (z = 0 on both) -- nothing tears; (b) the right one 40 vox off in z (another sheet placed next to the
     * first, as the PHerc1203 holes close-up's bridge) -- the junction tears, no remaining edge is longer than
     * tear_ratio x h, and both planes keep their interiors. */
    int pass = 0, fails = 0;
    for (pass = 0; pass < 2; pass++) {
        CloseupObs o;
        CloseupGrid g;
        CloseupFitStats st;
        size_t k = 0, too_long = 0, interior_lost = 0;
        double worst = 0.0;
        int x = 0, y = 0;
        memset(&o, 0, sizeof o);
        memset(&g, 0, sizeof g);
        if (cs_obs_reserve(&o, 200 * 60) != 0) return 1;
        for (y = 0; y < 60; y++)
            for (x = 0; x < 200; x++) {
                size_t w = o.n;
                o.uv[w * 2] = (float)(x + 0.5); o.uv[w * 2 + 1] = (float)(y + 0.5);
                o.p[w * 3] = (pass == 1 && x >= 100) ? 40.0f : 0.0f;
                o.p[w * 3 + 1] = (float)(y + 0.5); o.p[w * 3 + 2] = (float)(x + 0.5);
                o.w[w] = 1.0f; o.chart[w] = x < 100 ? 1 : 2;
                o.n++;
            }
        if (CloseupGrid_alloc(&g, -10.0, 210.0, -10.0, 70.0, 2.0) != 0 || CloseupSurface_fit(&g, &o, NULL, &st) != 0) {
            CloseupObs_free(&o); CloseupGrid_free(&g); return 1;
        }
        for (k = 0; k < g.ni * g.nj; k++) {
            const size_t i = k % g.ni, j = k / g.ni;
            const double u = g.u0 + g.h * (double)i, v = g.v0 + g.h * (double)j;
            int dir = 0;
            if (g.state[k] == CLOSEUP_NODE_NONE) {
                if (v > 6.0 && v < 54.0 && ((u > 6.0 && u < 84.0) || (u > 116.0 && u < 194.0))) interior_lost++;
                continue;
            }
            for (dir = 0; dir < 2; dir++) {
                size_t nb = 0;
                double d2 = 0.0;
                int q = 0;
                if (dir == 0 ? i + 1 >= g.ni : j + 1 >= g.nj) continue;
                nb = dir == 0 ? k + 1 : k + g.ni;
                if (g.state[nb] == CLOSEUP_NODE_NONE) continue;
                for (q = 0; q < 3; q++) {
                    const double e = (double)g.pos[nb * 3 + (size_t)q] - g.pos[k * 3 + (size_t)q];
                    d2 += e * e;
                }
                worst = fmax(worst, sqrt(d2) / g.h);
                if (sqrt(d2) > CLOSEUP_TEAR_RATIO * g.h) too_long++;
            }
        }
        fprintf(stderr, "  closeup tear (%s): torn stretch %zu, float %zu; longest edge %.2f x h; interior nodes lost %zu\n",
                pass ? "two sheets side by side" : "one sheet", st.torn_stretch, st.torn_float, worst, interior_lost);
        if (pass == 0 && (st.torn_stretch + st.torn_float != 0 || interior_lost != 0)) fails++;
        if (pass == 1 && (st.torn_stretch + st.torn_float == 0 || too_long != 0 || interior_lost != 0)) fails++;
        CloseupObs_free(&o);
        CloseupGrid_free(&g);
    }
    if (fails) fprintf(stderr, "FAIL closeup tear\n");
    return fails ? 1 : 0;
}

static int cs_selftest_complete(void)
{
    /* A known smooth sheet with a 120-px square hole and an open border notch. Verify actual
     * withheld geometry, including the unsupported centre, rather than just the count of fills. */
    CloseupObs o;
    CloseupGrid g;
    CloseupFitOpts opt;
    size_t k = 0, nfill = 0, centre = 0;
    double maxerr = 0.0, centre_err = 0.0;
    int x = 0, y = 0, rc = 1;
    memset(&o, 0, sizeof o);
    memset(&g, 0, sizeof g);
    if (cs_obs_reserve(&o, 400 * 320) != 0) goto done;
    for (y = 0; y < 320; y++)
        for (x = 0; x < 400; x++) {
            double u = x + 0.5, v = y + 0.5;
            size_t w = o.n;
            if ((x >= 140 && x < 260 && y >= 100 && y < 220) ||
                (x >= 300 && y >= 110 && y < 210)) continue;
            o.uv[w * 2] = (float)u; o.uv[w * 2 + 1] = (float)v;
            o.p[w * 3] = (float)(0.0001 * (u - 200.0) * (u - 200.0) + 0.00015 * (v - 160.0) * (v - 160.0));
            o.p[w * 3 + 1] = (float)v; o.p[w * 3 + 2] = (float)u;
            o.w[w] = 1.0f; o.chart[w] = 1;
            o.n++;
        }
    if (CloseupGrid_alloc(&g, 0.0, 400.0, 0.0, 320.0, 2.0) != 0 ||
        CloseupSurface_fit(&g, &o, NULL, NULL) != 0) goto done;
    centre = 80 * g.ni + 100;
    if (g.state[centre] != CLOSEUP_NODE_NONE) goto done;
    CloseupFitOpts_complete(&opt);
    if (CloseupSurface_fit(&g, &o, &opt, NULL) != 0 || g.state[centre] != CLOSEUP_NODE_FILL) goto done;
    for (k = 0; k < g.ni * g.nj; k++) {
        double u = g.u0 + (double)(k % g.ni) * g.h, v = g.v0 + (double)(k / g.ni) * g.h;
        double z = 0.0001 * (u - 200.0) * (u - 200.0) + 0.00015 * (v - 160.0) * (v - 160.0);
        double dz = 0.0, dy = 0.0, dx = 0.0, err = 0.0;
        if (g.state[k] != CLOSEUP_NODE_FILL) continue;
        nfill++;
        dz = g.pos[k * 3] - z; dy = g.pos[k * 3 + 1] - v; dx = g.pos[k * 3 + 2] - u;
        err = sqrt(dz * dz + dy * dy + dx * dx);
        if (k == centre) centre_err = err;
        if (err > maxerr) maxerr = err;
    }
    fprintf(stderr, "  closeup complete: %zu filled nodes, withheld geometry max error %.6f vox, centre %.6f vox\n",
            nfill, maxerr, centre_err);
    /* The existing curvature ridge biases one-sided notch fits: require sub-half-voxel accuracy
     * there, and the tighter symmetric-hole accuracy at the centre. */
    if (nfill < 5000 || maxerr > 0.5 || centre_err > 0.05) goto done;
    rc = 0;
done:
    CloseupObs_free(&o);
    CloseupGrid_free(&g);
    if (rc) fprintf(stderr, "FAIL closeup complete\n");
    return rc;
}

/* The frame of a cylinder patch: handedness from the inward normals, rotation from the axial gradient.  Negating the
 * axial field (a scan whose z runs from the scroll's bottom to its top: ink_closeup "z_top_to_bottom": false,
 * PHerc0343) must turn the frame by exactly 180 degrees and keep its reflection (a rotation preserves handedness);
 * a mirror-imaged layout must flip the reflection under either sense. */
static int cs_selftest_frame(void)
{
    enum { NU = 11, NV = 6 };
    const double R = 100.0, cu = 20.0, cv = 10.0, pi = 3.14159265358979;
    float verts[NU * NV * 3], inward[NU * NV * 3];
    double uv[NU * NV * 2], uvm[NU * NV * 2], ax[NU * NV], axn[NU * NV], fu = 0.0, fv = 0.0, gu = 0.0, gv = 0.0, d = 0.0;
    int32_t faces[(NU - 1) * (NV - 1) * 6];
    size_t nf = 0;
    CloseupFrame a, b, ma, mb;
    int i = 0, j = 0, rc = 0;
    for (j = 0; j < NV; j++)
        for (i = 0; i < NU; i++) {
            size_t k = (size_t)j * NU + (size_t)i;
            double u = 4.0 * i, v = 4.0 * j, th = u / R;
            cs_cyl(R, u, v, 0.0, &verts[k * 3]);
            inward[k * 3] = 0.0f; inward[k * 3 + 1] = (float)-sin(th); inward[k * 3 + 2] = (float)-cos(th);
            uv[k * 2] = u; uv[k * 2 + 1] = v;
            uvm[k * 2] = -u; uvm[k * 2 + 1] = v;
            ax[k] = v; axn[k] = -v;
        }
    for (j = 0; j + 1 < NV; j++)
        for (i = 0; i + 1 < NU; i++) {
            int32_t p = (int32_t)(j * NU + i);
            faces[nf * 3] = p; faces[nf * 3 + 1] = p + 1; faces[nf * 3 + 2] = p + NU + 1; nf++;
            faces[nf * 3] = p; faces[nf * 3 + 1] = p + NU + 1; faces[nf * 3 + 2] = p + NU; nf++;
        }
    if (CloseupSurface_frame(verts, uv, faces, nf, NULL, ax, inward, cu, cv, &a) != 0 ||
        CloseupSurface_frame(verts, uv, faces, nf, NULL, axn, inward, cu, cv, &b) != 0 ||
        CloseupSurface_frame(verts, uvm, faces, nf, NULL, ax, inward, -cu, cv, &ma) != 0 ||
        CloseupSurface_frame(verts, uvm, faces, nf, NULL, axn, inward, -cu, cv, &mb) != 0) {
        fprintf(stderr, "FAIL closeup frame: no frame\n");
        return 1;
    }
    d = fmod(fabs(b.theta - a.theta), 2.0 * pi);
    CloseupFrame_apply(&a, cu, cv + 1.0, &fu, &fv);          /* one voxel up the axis */
    CloseupFrame_apply(&b, cu, cv + 1.0, &gu, &gv);
    fprintf(stderr, "  closeup frame: reflect %d/%d (mirrored layout %d/%d), rotation %.4f -> %.4f rad (diff %.6f), "
            "axis step -> (%.3f, %.3f) / (%.3f, %.3f)\n", a.reflect, b.reflect, ma.reflect, mb.reflect, a.theta, b.theta,
            d, fu, fv, gu, gv);
    if (a.reflect != b.reflect || ma.reflect != mb.reflect || ma.reflect == a.reflect || fabs(d - pi) > 1e-9 ||
        fabs(fu) > 1e-9 || fabs(fv - 1.0) > 1e-9 || fabs(gu) > 1e-9 || fabs(gv + 1.0) > 1e-9) {
        fprintf(stderr, "FAIL closeup frame\n");
        rc = 1;
    }
    return rc;
}

int CloseupSurface_selftest(void)
{
    int f = 0;
    f += cs_selftest_raster();
    f += cs_selftest_infill();
    f += cs_selftest_bridge(0);
    f += cs_selftest_bridge_gated();
    f += cs_selftest_bridge(1);
    f += cs_selftest_complete();
    f += cs_selftest_frame();
    f += cs_selftest_fit();
    f += cs_selftest_tear();
    f += cs_selftest_flatten();
    f += cs_selftest_flatten_slim();
    f += cs_selftest_grid_file();
    f += cs_selftest_flatten_fringe();
    f += cs_selftest_dedup();
    f += cs_selftest_apron();
    fprintf(stderr, "[closeup_surface] selftest %s\n", f ? "FAILED" : "passed");
    return f ? 1 : 0;
}
