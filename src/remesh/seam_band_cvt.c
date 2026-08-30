/*
 * seam_band_cvt.c -- post-weld seam-band CVT beautification (see header).
 *
 * Scratch is malloc/free (local to one call, seam_weld convention); grown
 * vert/face arrays and the per-patch CVT results are arena allocations. No
 * Arena_save/restore across CVT calls (outputs must not sit above a mark).
 */
#include "seam_band_cvt.h"
#include "cvt_remesh.h"
#include "../common/mesh_manifold.h"
#include "../common/pipeline_constants.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void SeamBandCvt_default_params(SeamBandCvtParams *p)
{
    assert(p);
    p->band        = SEAM_BANDCVT_BAND_VOX;
    p->target_h    = SEAM_BANDCVT_H_VOX;
    p->iters       = SEAM_BANDCVT_ITERS;
    p->min_faces   = SEAM_BANDCVT_MIN_FACES;
    p->tile        = SEAM_BANDCVT_TILE_VOX;
    p->offset_half = 0;
    p->fine_frac   = SEAM_BANDCVT_FINE_FRAC;
}

/* union-find with path halving */
static int32_t uf_find(int32_t *par, int32_t x)
{
    while (par[x] != x) { par[x] = par[par[x]]; x = par[x]; }
    return x;
}

typedef struct { int64_t k; int32_t count; } EdgeRun;

static int i64p_cmp(const void *pa, const void *pb)
{
    int64_t a = *(const int64_t *)pa, b = *(const int64_t *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static int64_t ekey(int32_t a, int32_t b)
{
    int32_t u = a < b ? a : b, v = a < b ? b : a;
    return ((int64_t)u << 32) | (int64_t)(uint32_t)v;
}

/* Sorted unique-edge table with multiplicities for a face list. Caller frees
 * *keys. Returns unique count; *out_run1 = number of run-1 (boundary) edges. */
static size_t edge_table(const int32_t *F, size_t nf,
                         int64_t **keys_out, size_t *out_run1)
{
    size_t ne = nf * 3, uniq = 0, run1 = 0, i;
    int64_t *keys = (int64_t *)malloc((ne ? ne : 1) * sizeof(int64_t));
    for (size_t t = 0; t < nf; t++)
        for (int e = 0; e < 3; e++)
            keys[t*3+(size_t)e] = ekey(F[t*3+e], F[t*3+(e+1)%3]);
    qsort(keys, ne, sizeof(int64_t), i64p_cmp);
    for (i = 0; i < ne; ) {
        size_t j = i + 1;
        while (j < ne && keys[j] == keys[i]) j++;
        keys[uniq++] = keys[i];          /* compact uniques in place */
        if (j - i == 1) run1++;
        i = j;
    }
    *keys_out = keys;
    if (out_run1) *out_run1 = run1;
    return uniq;
}

/* boundary (run-1) edge keys only, sorted. Caller frees. */
static size_t boundary_edges_sorted(const int32_t *F, size_t nf, int64_t **out)
{
    size_t ne = nf * 3, nb = 0, i;
    int64_t *keys = (int64_t *)malloc((ne ? ne : 1) * sizeof(int64_t));
    int64_t *bnd  = (int64_t *)malloc((ne ? ne : 1) * sizeof(int64_t));
    for (size_t t = 0; t < nf; t++)
        for (int e = 0; e < 3; e++)
            keys[t*3+(size_t)e] = ekey(F[t*3+e], F[t*3+(e+1)%3]);
    qsort(keys, ne, sizeof(int64_t), i64p_cmp);
    for (i = 0; i < ne; ) {
        size_t j = i + 1;
        while (j < ne && keys[j] == keys[i]) j++;
        if (j - i == 1) bnd[nb++] = keys[i];
        i = j;
    }
    free(keys);
    *out = bnd;
    return nb;
}

/* max edge multiplicity (manifold check) */
static int max_edge_run(const int32_t *F, size_t nf)
{
    size_t ne = nf * 3, i;
    int maxr = 0;
    int64_t *keys = (int64_t *)malloc((ne ? ne : 1) * sizeof(int64_t));
    for (size_t t = 0; t < nf; t++)
        for (int e = 0; e < 3; e++)
            keys[t*3+(size_t)e] = ekey(F[t*3+e], F[t*3+(e+1)%3]);
    qsort(keys, ne, sizeof(int64_t), i64p_cmp);
    for (i = 0; i < ne; ) {
        size_t j = i + 1;
        while (j < ne && keys[j] == keys[i]) j++;
        if ((int)(j - i) > maxr) maxr = (int)(j - i);
        i = j;
    }
    free(keys);
    return maxr;
}

/* (patch root, face id) pair for grouping band faces into patches */
typedef struct { int32_t root; int32_t face; } FR;
static int fr_cmp(const void *pa, const void *pb)
{
    const FR *a = (const FR *)pa, *b = (const FR *)pb;
    if (a->root != b->root) return a->root < b->root ? -1 : 1;
    return a->face < b->face ? -1 : (a->face > b->face ? 1 : 0);
}

/* (tile key, face id) pair for the tile sort */
typedef struct { int64_t k; int32_t face; } KFc;
static int kf64_cmp(const void *pa, const void *pb)
{
    const KFc *a = (const KFc *)pa, *b = (const KFc *)pb;
    if (a->k != b->k) return a->k < b->k ? -1 : 1;
    return a->face < b->face ? -1 : (a->face > b->face ? 1 : 0);
}

static int i32_cmp(const void *pa, const void *pb)
{
    int32_t a = *(const int32_t *)pa, b = *(const int32_t *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

/* index of v in the sorted array a[n] (must be present) */
static int32_t i32_pos(const int32_t *a, size_t n, int32_t v)
{
    size_t lo = 0, hi = n;
    while (lo < hi) { size_t mid = lo + (hi - lo)/2; if (a[mid] < v) lo = mid+1; else hi = mid; }
    return (int32_t)lo;
}

/* connected components over faces via shared verts */
static size_t face_components(const int32_t *F, size_t nf, size_t nv)
{
    int32_t *par = (int32_t *)malloc((nv ? nv : 1) * sizeof(int32_t));
    uint8_t *used = (uint8_t *)calloc(nv ? nv : 1, 1);
    size_t i, n = 0;
    for (i = 0; i < nv; i++) par[i] = (int32_t)i;
    for (size_t t = 0; t < nf; t++) {
        int32_t a = F[t*3+0], b = F[t*3+1], c = F[t*3+2];
        used[a] = used[b] = used[c] = 1;
        int32_t ra = uf_find(par, a), rb = uf_find(par, b), rc = uf_find(par, c);
        if (rb != ra) par[rb] = ra;
        rc = uf_find(par, c);
        if (rc != uf_find(par, a)) par[rc] = uf_find(par, a);
    }
    for (i = 0; i < nv; i++)
        if (used[i] && uf_find(par, (int32_t)i) == (int32_t)i) n++;
    free(par); free(used);
    return n;
}

int SeamBandCvt_process(Arena_T arena,
                        const float *verts, size_t nv,
                        const int32_t *faces, size_t nf,
                        const SeamPlane *planes, size_t np,
                        const SeamBandCvtParams *params,
                        float **out_verts, size_t *out_nv,
                        int32_t **out_faces, size_t *out_nf,
                        int32_t **out_new_vert_src, size_t *out_n_new,
                        SeamBandCvtStats *st)
{
    SeamBandCvtParams p;
    assert(arena);
    assert(out_verts && out_nv && out_faces && out_nf);
    assert(out_new_vert_src && out_n_new);

    if (params) p = *params; else SeamBandCvt_default_params(&p);
    if (st) memset(st, 0, sizeof *st);

    /* rollback/debug knob: disable the vertex-manifold (bowtie) gate. With it
     * OFF, band-CVT reverts to the pre-fix behavior that severed same-turn seams
     * (an A/B handle, not for production). */
    int vtx_gate = (getenv("SEAM_BANDCVT_NO_VTXGATE") == NULL);

    *out_verts = (float *)verts;   *out_nv = nv;
    *out_faces = (int32_t *)faces; *out_nf = nf;
    *out_new_vert_src = NULL;      *out_n_new = 0;
    if (nf == 0 || nv == 0 || np == 0 || planes == NULL) return 0;
    if (p.target_h <= 0.0 || p.iters <= 0) return 0;

    /* 1) band faces -> (nearest plane, in-plane tile) groups -> connected
     * sub-patches. Tiling bounds the fail-closed blast radius: one gate
     * failure (a core fold the RVD would jump, a dropped junction edge)
     * rejects ONE tile's sub-patch, not the whole band lattice -- measured on
     * the 4x5x5, the untiled version fused the band into one grid-spanning
     * patch that always tripped a gate somewhere and rejected wholesale.
     * Pass B (offset_half) shifts the tiling by tile/2 so pass A's pinned
     * tile borders land tile-INTERIOR and are re-meshed away. */
    uint8_t *in_band = (uint8_t *)calloc(nf, 1);
    int64_t *fkey = (int64_t *)malloc((nf ? nf : 1) * sizeof(int64_t));
    if (!in_band || !fkey) { free(in_band); free(fkey); return -1; }
    size_t n_band = 0;
    double off = p.offset_half ? 0.5 * p.tile : 0.0;
    for (size_t f = 0; f < nf; f++) {
        int near_plane = 0;
        for (int e = 0; e < 3; e++) {
            if (SeamPlanes_vert_dist(verts, faces[f*3+e], planes, np) <= p.band) {
                near_plane = 1; break;
            }
        }
        if (!near_plane) continue;
        in_band[f] = 1; n_band++;
        {
            double c[3];
            for (int k = 0; k < 3; k++)
                c[k] = ((double)verts[(size_t)faces[f*3+0]*3+(size_t)k]
                      + (double)verts[(size_t)faces[f*3+1]*3+(size_t)k]
                      + (double)verts[(size_t)faces[f*3+2]*3+(size_t)k]) / 3.0;
            size_t bp = 0; double bd = 1e300;
            for (size_t q = 0; q < np; q++) {
                double d = fabs(c[planes[q].axis] - planes[q].coord);
                if (d < bd) { bd = d; bp = q; }
            }
            int ax = planes[bp].axis;
            int a1 = (ax + 1) % 3, a2 = (ax + 2) % 3;
            long tu = (long)floor((c[a1] + off) / p.tile);
            long tv = (long)floor((c[a2] + off) / p.tile);
            fkey[f] = ((int64_t)bp << 42)
                    | (((int64_t)(tu + (1L << 20)) & 0x1FFFFF) << 21)
                    |  ((int64_t)(tv + (1L << 20)) & 0x1FFFFF);
        }
    }
    if (n_band == 0) { free(in_band); free(fkey); return 0; }

    /* sort band faces by tile key */
    KFc *kf = (KFc *)malloc(n_band * sizeof(KFc));
    if (!kf) { free(in_band); free(fkey); return -1; }
    {
        size_t w = 0;
        for (size_t f = 0; f < nf; f++)
            if (in_band[f]) { kf[w].k = fkey[f]; kf[w].face = (int32_t)f; w++; }
    }
    free(fkey);
    qsort(kf, n_band, sizeof(KFc), kf64_cmp);

    /* per tile: local union-find over shared verts -> global sub-patch ids */
    FR *fr = (FR *)malloc(n_band * sizeof(FR));
    if (!fr) { free(in_band); free(kf); return -1; }
    {
        int32_t next_patch = 0;
        size_t ti = 0;
        while (ti < n_band) {
            size_t tj = ti + 1;
            while (tj < n_band && kf[tj].k == kf[ti].k) tj++;
            size_t tn = tj - ti;
            /* gather + sort unique vert ids of this tile's faces */
            int32_t *vl = (int32_t *)malloc(tn * 3 * sizeof(int32_t));
            for (size_t x = 0; x < tn; x++)
                for (int e = 0; e < 3; e++)
                    vl[x*3+(size_t)e] = faces[(size_t)kf[ti+x].face*3+(size_t)e];
            qsort(vl, tn * 3, sizeof(int32_t), i32_cmp);
            size_t nvl = 0;
            for (size_t x = 0; x < tn * 3; x++)
                if (x == 0 || vl[x] != vl[x-1]) vl[nvl++] = vl[x];
            int32_t *par = (int32_t *)malloc(nvl * sizeof(int32_t));
            int32_t *rmap = (int32_t *)malloc(nvl * sizeof(int32_t));
            for (size_t x = 0; x < nvl; x++) { par[x] = (int32_t)x; rmap[x] = -1; }
            for (size_t x = 0; x < tn; x++) {
                int32_t l0 = i32_pos(vl, nvl, faces[(size_t)kf[ti+x].face*3+0]);
                int32_t l1 = i32_pos(vl, nvl, faces[(size_t)kf[ti+x].face*3+1]);
                int32_t l2 = i32_pos(vl, nvl, faces[(size_t)kf[ti+x].face*3+2]);
                int32_t r0 = uf_find(par, l0), r1 = uf_find(par, l1);
                if (r1 != r0) par[r1] = r0;
                int32_t r2 = uf_find(par, l2);
                r0 = uf_find(par, l0);
                if (r2 != r0) par[r2] = r0;
            }
            for (size_t x = 0; x < tn; x++) {
                int32_t l0 = i32_pos(vl, nvl, faces[(size_t)kf[ti+x].face*3+0]);
                int32_t r = uf_find(par, l0);
                if (rmap[r] < 0) rmap[r] = next_patch++;
                fr[ti+x].root = rmap[r];
                fr[ti+x].face = kf[ti+x].face;
            }
            free(vl); free(par); free(rmap);
            ti = tj;
        }
    }
    free(kf);
    qsort(fr, n_band, sizeof(FR), fr_cmp);

    /* accepted-patch replacement bookkeeping */
    uint8_t *face_replaced = (uint8_t *)calloc(nf, 1);
    int32_t *g2l = (int32_t *)malloc(nv * sizeof(int32_t));   /* global->local */
    if (!face_replaced || !g2l) {
        free(in_band); free(fr); free(face_replaced); free(g2l);
        return -1;
    }
    for (size_t v = 0; v < nv; v++) g2l[v] = -1;

    /* growing output for accepted dual faces (global indices) + new verts */
    size_t rep_cap = 1024, rep_nf = 0;
    int32_t *rep_f = (int32_t *)malloc(rep_cap * 3 * sizeof(int32_t));
    size_t newv_cap = 1024, n_newv = 0;
    float  *newv   = (float *)malloc(newv_cap * 3 * sizeof(float));
    int32_t *newsrc = (int32_t *)malloc(newv_cap * sizeof(int32_t));
    int rc_all = 0;

    size_t pi = 0;
    while (pi < n_band && rc_all == 0) {
        size_t pj = pi + 1;
        while (pj < n_band && fr[pj].root == fr[pi].root) pj++;
        size_t pf_n = pj - pi;
        if (st) st->patches++;
        if (pf_n < p.min_faces) {
            if (st) st->skipped_small++;
            pi = pj; continue;
        }

        /* DESPUR: iteratively shed faces with >= 2 patch-boundary edges (ear
         * faces / spur tips). The any-vert band rule + centroid tiling leave
         * sawtooth borders with one-face spurs, and a spur-tip pin's cell
         * touches only its 2 loop neighbours -> no dual corner -> the pin
         * gate rejects the whole patch. Shed faces just stay un-beautified. */
        int32_t *plist = (int32_t *)malloc(pf_n * sizeof(int32_t));
        size_t pn = pf_n;
        for (size_t k = 0; k < pf_n; k++) plist[k] = fr[pi + k].face;
        for (int dspass = 0; dspass < 16 && pn >= 3; dspass++) {
            size_t ne = pn * 3, w = 0;
            int64_t *ek = (int64_t *)malloc(ne * sizeof(int64_t));
            for (size_t k = 0; k < pn; k++)
                for (int e = 0; e < 3; e++)
                    ek[k*3+(size_t)e] = ekey(faces[(size_t)plist[k]*3+(size_t)e],
                                             faces[(size_t)plist[k]*3+(size_t)((e+1)%3)]);
            int64_t *sk = (int64_t *)malloc(ne * sizeof(int64_t));
            memcpy(sk, ek, ne * sizeof(int64_t));
            qsort(sk, ne, sizeof(int64_t), i64p_cmp);
            for (size_t k = 0; k < pn; k++) {
                int nb = 0;
                for (int e = 0; e < 3; e++) {
                    int64_t key = ek[k*3+(size_t)e];
                    size_t lo = 0, hi = ne;
                    while (lo < hi) { size_t mid = lo + (hi-lo)/2;
                                      if (sk[mid] < key) lo = mid+1; else hi = mid; }
                    size_t run = 0;
                    while (lo+run < ne && sk[lo+run] == key) run++;
                    if (run == 1) nb++;
                }
                if (nb < 2) plist[w++] = plist[k];
            }
            free(ek); free(sk);
            if (w == pn) break;
            pn = w;
        }
        if (pn < p.min_faces) {
            if (st) st->skipped_small++;
            free(plist); pi = pj; continue;
        }

        /* already CVT-quality? (pass B skips pass A's accepted interiors) */
        if (p.fine_frac > 0.0) {
            size_t fine = 0, tote = 0;
            double thr2 = 0.36 * p.target_h * p.target_h;   /* (0.6h)^2 */
            for (size_t k = 0; k < pn; k++) {
                const int32_t *t = &faces[(size_t)plist[k]*3];
                for (int e = 0; e < 3; e++) {
                    int32_t u = t[e], v = t[(e+1)%3];
                    double dz=(double)verts[(size_t)u*3+0]-verts[(size_t)v*3+0];
                    double dy=(double)verts[(size_t)u*3+1]-verts[(size_t)v*3+1];
                    double dx=(double)verts[(size_t)u*3+2]-verts[(size_t)v*3+2];
                    if (dz*dz+dy*dy+dx*dx < thr2) fine++;
                    tote++;
                }
            }
            if ((double)fine < p.fine_frac * (double)tote) {
                if (st) st->skipped_clean++;
                free(plist); pi = pj; continue;
            }
        }

        /* local submesh */
        int32_t *lf = (int32_t *)malloc(pn * 3 * sizeof(int32_t));
        int32_t *l2g = (int32_t *)malloc(pn * 3 * sizeof(int32_t));
        size_t lnv = 0, pf_used = pn;
        for (size_t k = 0; k < pn; k++) {
            int32_t gf = plist[k];
            for (int e = 0; e < 3; e++) {
                int32_t gvid = faces[(size_t)gf*3+(size_t)e];
                if (g2l[gvid] < 0) { g2l[gvid] = (int32_t)lnv; l2g[lnv] = gvid; lnv++; }
                lf[k*3+(size_t)e] = g2l[gvid];
            }
        }
        pf_n = pf_used;
        float *lv = (float *)malloc(lnv * 3 * sizeof(float));
        for (size_t v = 0; v < lnv; v++) {
            lv[v*3+0] = verts[(size_t)l2g[v]*3+0];
            lv[v*3+1] = verts[(size_t)l2g[v]*3+1];
            lv[v*3+2] = verts[(size_t)l2g[v]*3+2];
        }

        /* input patch topology (for the gates) */
        int64_t *in_bnd = NULL;
        size_t in_nb = boundary_edges_sorted(lf, pf_n, &in_bnd);
        int64_t *in_keys = NULL; size_t in_r1 = 0;
        size_t in_E = edge_table(lf, pf_n, &in_keys, &in_r1);
        long chi_in = (long)lnv - (long)in_E + (long)pf_n;
        size_t comp_in = face_components(lf, pf_n, lnv);
        free(in_keys);

        /* CVT the patch. Everything the CVT arena-allocates (double copy,
         * KD tree, sites, the dual) is scratch beyond this patch: the accept
         * path below copies verts/faces/src into malloc buffers, so the mark
         * is restored at the bottom of the loop and the weld arena stays
         * bounded at ONE patch's working set (~4k patches on the 4x5x5 would
         * otherwise pile up gigabytes of dead CVT scratch). */
        Arena_Mark cvt_mark = Arena_save(arena);
        CvtOpts co; CVT_default_opts(&co);
        co.n_iters = p.iters;
        float *cv = NULL; int32_t *cf = NULL; int32_t *csrc = NULL;
        size_t cnv = 0, cnf = 0, cpin = 0;
        int rc = CVT_remesh_pinned(arena, lv, lnv, lf, pf_n, p.target_h, &co,
                                   &cv, &cnv, &cf, &cnf, &csrc, &cpin);

        int ok = (rc == 0 && cnf > 0 && cpin >= 3);
        int why = ok ? 0 : 1;             /* 1=rc 2=pin 3=bnd 4=chi 5=nm 6=comp 7=vtx */

        /* gate 1: every pinned site referenced */
        if (ok) {
            uint8_t *pin_used = (uint8_t *)calloc(cnv, 1);
            for (size_t t = 0; t < cnf*3; t++) pin_used[cf[t]] = 1;
            for (size_t i = 0; i < cpin && ok; i++)
                if (!pin_used[i]) {
                    ok = 0; why = 2;
                    fprintf(stderr, "    [bandcvt] reject: pin %zu/%zu unused "
                            "at (%.1f, %.1f, %.1f)\n", i, cpin,
                            (double)cv[i*3+0], (double)cv[i*3+1],
                            (double)cv[i*3+2]);
                }
            free(pin_used);
        }
        /* gates 2+5: boundary set equal (mapped to input-local ids) + chi */
        if (ok) {
            /* map dual faces to input-local ids where pinned; interior sites
             * get temp ids lnv + (site - cpin) */
            int32_t *mf = (int32_t *)malloc(cnf * 3 * sizeof(int32_t));
            for (size_t t = 0; t < cnf*3; t++) {
                int32_t s = cf[t];
                mf[t] = (s < (int32_t)cpin) ? csrc[s]
                                            : (int32_t)lnv + (s - (int32_t)cpin);
            }
            int64_t *out_bnd = NULL;
            size_t out_nb = boundary_edges_sorted(mf, cnf, &out_bnd);
            if (out_nb != in_nb ||
                (in_nb && memcmp(in_bnd, out_bnd, in_nb*sizeof(int64_t)) != 0)) {
                ok = 0; why = 3;
                /* report the first asymmetric edge, both directions */
                for (size_t x = 0, y = 0; x < in_nb || y < out_nb; ) {
                    int64_t a = x < in_nb ? in_bnd[x] : INT64_MAX;
                    int64_t b = y < out_nb ? out_bnd[y] : INT64_MAX;
                    if (a == b) { x++; y++; continue; }
                    int in_only = a < b;
                    int64_t k = in_only ? a : b;
                    int32_t u = (int32_t)(k >> 32), v = (int32_t)(k & 0xffffffff);
                    const float *pu = u < (int32_t)lnv ? &lv[(size_t)u*3] : NULL;
                    fprintf(stderr, "    [bandcvt] reject: bnd edge %s dual: "
                            "(%d,%d)%s%.1f,%.1f,%.1f)\n",
                            in_only ? "missing from" : "invented by", u, v,
                            pu ? " at (" : " (interior-site edge; (",
                            pu ? (double)pu[0] : 0.0,
                            pu ? (double)pu[1] : 0.0,
                            pu ? (double)pu[2] : 0.0);
                    break;
                }
            }
            if (ok) {
                int64_t *out_keys = NULL; size_t out_r1 = 0;
                size_t out_E = edge_table(mf, cnf, &out_keys, &out_r1);
                /* used verts of the dual: cnv minus empty-cell orphans */
                uint8_t *vu = (uint8_t *)calloc(cnv, 1);
                size_t used_v = 0;
                for (size_t t = 0; t < cnf*3; t++)
                    if (!vu[cf[t]]) { vu[cf[t]] = 1; used_v++; }
                free(vu);
                long chi_out = (long)used_v - (long)out_E + (long)cnf;
                if (chi_out != chi_in) { ok = 0; why = 4; }
                free(out_keys);
            }
            /* gate 3: edge-manifold */
            if (ok && max_edge_run(mf, cnf) > 2) { ok = 0; why = 5; }
            /* gate 3b: vertex-manifold -- reject a dual with a bowtie/pinch
             * vertex. A vertex-only pinch keeps every EDGE at <=2 faces so the
             * per-edge gate above passes it, and it keeps ONE vertex-connected
             * component so gate 4 passes it too -- but the downstream
             * ManifoldGuard splits the pinch by vertex duplication, severing a
             * same-turn connection the bridge had just made. That is THE
             * seam-split root cause (wind_audit: band-CVT tripled splits-at-seam
             * 18.5%->62.6%). nm_verts is exactly the bowtie count; audited on
             * the dual in isolation (cf/cnv), which is index-isomorphic to the
             * assembled patch, and gate 2 (boundary-edge set identical) means a
             * pinned boundary vertex cannot gain a NEW fan from the rest of the
             * mesh -- so an interior bowtie is the only way to sever. */
            if (ok && vtx_gate) {
                MeshManifoldStats vms = MeshManifold_audit(arena, cnv, cf, cnf);
                if (vms.nm_verts > 0) { ok = 0; why = 7; }
            }
            /* gate 4: component count unchanged (despur can split a patch) */
            if (ok && face_components(mf, cnf, lnv + (cnv - cpin)) != comp_in) {
                ok = 0; why = 6;
            }
            free(out_bnd);
            free(mf);
        }

        if (ok) {
            /* accept: mark faces replaced, append new verts + mapped faces */
            for (size_t k = 0; k < pf_n; k++)
                face_replaced[plist[k]] = 1;
            /* interior sites -> new global verts; provenance = nearest pinned
             * site's source vert */
            size_t base_new = n_newv;
            for (size_t s = cpin; s < cnv; s++) {
                if (n_newv == newv_cap) {
                    newv_cap *= 2;
                    newv = (float *)realloc(newv, newv_cap * 3 * sizeof(float));
                    newsrc = (int32_t *)realloc(newsrc, newv_cap * sizeof(int32_t));
                    if (!newv || !newsrc) { rc_all = -1; break; }
                }
                newv[n_newv*3+0] = cv[s*3+0];
                newv[n_newv*3+1] = cv[s*3+1];
                newv[n_newv*3+2] = cv[s*3+2];
                {
                    double best = 1e300; int32_t bsrc = l2g[csrc[0]];
                    for (size_t i = 0; i < cpin; i++) {
                        double dz=(double)cv[s*3+0]-cv[i*3+0];
                        double dy=(double)cv[s*3+1]-cv[i*3+1];
                        double dx=(double)cv[s*3+2]-cv[i*3+2];
                        double d2=dz*dz+dy*dy+dx*dx;
                        if (d2 < best) { best = d2; bsrc = l2g[csrc[i]]; }
                    }
                    newsrc[n_newv] = bsrc;
                }
                n_newv++;
            }
            for (size_t t = 0; t < cnf && rc_all == 0; t++) {
                if (rep_nf == rep_cap) {
                    rep_cap *= 2;
                    rep_f = (int32_t *)realloc(rep_f, rep_cap * 3 * sizeof(int32_t));
                    if (!rep_f) { rc_all = -1; break; }
                }
                for (int e = 0; e < 3; e++) {
                    int32_t s = cf[t*3+e];
                    rep_f[rep_nf*3+(size_t)e] =
                        (s < (int32_t)cpin)
                            ? l2g[csrc[s]]
                            : (int32_t)(nv + base_new + ((size_t)s - cpin));
                }
                rep_nf++;
            }
            if (st) {
                st->accepted++;
                st->faces_band_in += pf_n;
                st->faces_band_out += cnf;
            }
        } else {
            if (st) {
                st->rejected++;
                switch (why) {
                    case 1: st->rej_rc++; break;
                    case 2: st->rej_pin++; break;
                    case 3: st->rej_bnd++; break;
                    case 4: st->rej_chi++; break;
                    case 5: st->rej_manifold++; break;
                    case 6: st->rej_comp++; break;
                    case 7: st->rej_vtx++; break;
                    default: break;
                }
            }
            /* debug: dump the rejected patch (input + dual) as OBJs */
            { const char *dd = getenv("SEAM_BANDCVT_DEBUG_DIR");
              if (dd) {
                static int dump_id = 0;
                char pth[1024];
                snprintf(pth, sizeof pth, "%s/reject%03d_why%d_in.obj",
                         dd, dump_id, why);
                FILE *fp = fopen(pth, "w");
                if (fp) {
                    for (size_t v = 0; v < lnv; v++)
                        fprintf(fp, "v %f %f %f\n", (double)lv[v*3+2],
                                (double)lv[v*3+1], (double)lv[v*3+0]);
                    for (size_t t = 0; t < pf_n; t++)
                        fprintf(fp, "f %d %d %d\n", lf[t*3+0]+1, lf[t*3+1]+1,
                                lf[t*3+2]+1);
                    fclose(fp);
                }
                if (rc == 0 && cnf > 0) {
                    snprintf(pth, sizeof pth, "%s/reject%03d_why%d_dual.obj",
                             dd, dump_id, why);
                    fp = fopen(pth, "w");
                    if (fp) {
                        for (size_t v = 0; v < cnv; v++)
                            fprintf(fp, "v %f %f %f\n", (double)cv[v*3+2],
                                    (double)cv[v*3+1], (double)cv[v*3+0]);
                        for (size_t t = 0; t < cnf; t++)
                            fprintf(fp, "f %d %d %d\n", cf[t*3+0]+1,
                                    cf[t*3+1]+1, cf[t*3+2]+1);
                        fclose(fp);
                    }
                }
                dump_id++;
              } }
        }

        /* reset g2l for this patch's verts; release this patch's CVT scratch
         * (accepted results were copied out to malloc buffers above) */
        Arena_restore(arena, cvt_mark);
        for (size_t v = 0; v < lnv; v++) g2l[l2g[v]] = -1;
        free(lf); free(l2g); free(lv); free(in_bnd); free(plist);
        pi = pj;
    }

    free(in_band); free(fr); free(g2l);

    if (rc_all != 0 || rep_nf == 0) {
        free(face_replaced); free(rep_f); free(newv); free(newsrc);
        return rc_all;   /* 0 with nothing accepted: no-op outputs stand */
    }

    /* 2) rebuild: surviving faces + replacements; grown vert array */
    size_t keep = 0;
    for (size_t f = 0; f < nf; f++) if (!face_replaced[f]) keep++;
    size_t fin_nf = keep + rep_nf;
    size_t fin_nv = nv + n_newv;
    float *fv = (float *)ARENA_ALLOC(arena, (long)(fin_nv*3*sizeof(float)));
    int32_t *ff = (int32_t *)ARENA_ALLOC(arena, (long)(fin_nf*3*sizeof(int32_t)));
    int32_t *fsrc = (int32_t *)ARENA_ALLOC(arena,
                        (long)((n_newv ? n_newv : 1)*sizeof(int32_t)));
    memcpy(fv, verts, nv*3*sizeof(float));
    memcpy(fv + nv*3, newv, n_newv*3*sizeof(float));
    memcpy(fsrc, newsrc, n_newv*sizeof(int32_t));
    {
        size_t w = 0;
        for (size_t f = 0; f < nf; f++) {
            if (face_replaced[f]) continue;
            ff[w*3+0]=faces[f*3+0]; ff[w*3+1]=faces[f*3+1]; ff[w*3+2]=faces[f*3+2];
            w++;
        }
        memcpy(ff + w*3, rep_f, rep_nf*3*sizeof(int32_t));
    }
    free(face_replaced); free(rep_f); free(newv); free(newsrc);

    *out_verts = fv; *out_nv = fin_nv;
    *out_faces = ff; *out_nf = fin_nf;
    *out_new_vert_src = fsrc; *out_n_new = n_newv;
    if (st) st->verts_added = n_newv;
    return 0;
}
