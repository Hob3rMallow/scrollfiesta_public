/*
 * orient_weld.c -- post-weld cross-component orientation. See orient_weld.h.
 *
 * Internal scratch uses malloc/free (local to one call, same convention as
 * ball_pivot.c / seam_weld.c); the arena parameter is unused.
 */
#include "orient_weld.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Robust radial reference (matches obj_orient_audit's criterion): per-face
 * radial votes are AREA-weighted and gated by this normal-vs-radial clarity,
 * and a component counts as "radially decisive" only when this fraction of its
 * gated area agrees on one side. The old fallback summed unweighted per-vertex
 * normals against the largest component's per-vertex sum; on a folded,
 * radially-mixed core that reference went noisy and stranded clean detached
 * components with the wrong global winding. */
#define ORIENT_WELD_MIN_COS       0.3
#define ORIENT_WELD_RADIAL_AGREE  0.8

static int uf_find(int *uf, int x) { while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; } return x; }

typedef struct { long nf; int id; } CompOrder;
static int cmp_comporder(const void *pa, const void *pb)
{
    const CompOrder *a = (const CompOrder *)pa, *b = (const CompOrder *)pb;
    if (a->nf != b->nf) return a->nf < b->nf ? 1 : -1;   /* face count descending */
    return a->id < b->id ? -1 : (a->id > b->id);          /* tie-break: id ascending */
}

/* In-plane radial unit vector at p relative to the axis; returns 0 if the
 * point is on the axis (radial direction undefined), else 1. ad must be a
 * unit vector. */
static int radial_unit(const float *p, const float *axis_point,
                       const double ad[3], double rhat[3])
{
    double d[3], dax = 0.0, len = 0.0;
    int k = 0;
    for (k = 0; k < 3; k++) d[k] = (double)p[k] - (double)axis_point[k];
    dax = d[0]*ad[0] + d[1]*ad[1] + d[2]*ad[2];
    for (k = 0; k < 3; k++) d[k] -= dax * ad[k];
    len = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
    if (len < 1.0) return 0;
    for (k = 0; k < 3; k++) rhat[k] = d[k] / len;
    return 1;
}

int OrientWeld_components(Arena_T arena,
                          const float *verts, size_t nv,
                          int32_t *faces, size_t nf,
                          float radius,
                          size_t *out_flipped)
{
    return OrientWeld_components_axis(arena, verts, nv, faces, nf, radius,
                                      NULL, NULL, out_flipped, NULL);
}

int OrientWeld_components_axis(Arena_T arena,
                               const float *verts, size_t nv,
                               int32_t *faces, size_t nf,
                               float radius,
                               const float *axis_point,
                               const float *axis_dir,
                               size_t *out_flipped,
                               size_t *out_radial)
{
    (void)arena;
    double ad[3] = { 0.0, 0.0, 0.0 };
    int have_axis = 0;
    if (out_flipped) *out_flipped = 0;
    if (out_radial) *out_radial = 0;
    if (nf == 0 || nv == 0) return 0;
    if (axis_point != NULL && axis_dir != NULL) {
        double L = sqrt((double)axis_dir[0]*axis_dir[0]
                        + (double)axis_dir[1]*axis_dir[1]
                        + (double)axis_dir[2]*axis_dir[2]);
        if (L > 1e-9) {
            ad[0] = axis_dir[0]/L; ad[1] = axis_dir[1]/L; ad[2] = axis_dir[2]/L;
            have_axis = 1;
        }
    }
    int n = (int)nv;

    /* 1) Vertex-connected components (union the 3 verts of each face). */
    int *uf = (int *)malloc((size_t)n * sizeof(int));
    for (int i = 0; i < n; i++) uf[i] = i;
    for (size_t f = 0; f < nf; f++) {
        int a = faces[f*3+0], b = faces[f*3+1], c = faces[f*3+2];
        int ra = uf_find(uf,a), rb = uf_find(uf,b);
        if (ra != rb) uf[rb] = ra;
        ra = uf_find(uf,a); int rc = uf_find(uf,c);
        if (ra != rc) uf[rc] = ra;
    }
    int *cid = (int *)malloc((size_t)n * sizeof(int));
    for (int i = 0; i < n; i++) cid[i] = -1;
    int ncomp = 0;
    for (int i = 0; i < n; i++) { int r = uf_find(uf,i); if (cid[r] < 0) cid[r] = ncomp++; }
    for (int i = 0; i < n; i++) cid[i] = cid[uf_find(uf,i)];   /* vertex -> component */
    free(uf);
    if (ncomp < 2) { free(cid); return 0; }                   /* nothing to anchor against */

    long *cnf = (long *)calloc((size_t)ncomp, sizeof(long));
    for (size_t f = 0; f < nf; f++) cnf[cid[faces[f*3+0]]]++;

    /* 1a) CSR vertex + face lists per component, so the vote/flip loops walk
     *     ONLY the component's own elements. The old code swept ALL n verts
     *     once per component (and all nf faces per flip); edge collapses
     *     orphan tens of thousands of single-vertex components, so that was
     *     an ncomp x n quadratic — measured 148 s of a 351 s 4x5x5 weld.
     *     Lists are filled in ascending vertex/face order (stable counting
     *     sort), so the vote accumulation and flip order — and therefore the
     *     output — are bit-identical to the full sweeps. */
    long *cvoff = (long *)calloc((size_t)ncomp + 1, sizeof(long));
    long *cfoff = (long *)calloc((size_t)ncomp + 1, sizeof(long));
    for (int i = 0; i < n; i++) cvoff[cid[i] + 1]++;
    for (size_t f = 0; f < nf; f++) cfoff[cid[faces[f*3+0]] + 1]++;
    for (int c = 0; c < ncomp; c++) { cvoff[c+1] += cvoff[c]; cfoff[c+1] += cfoff[c]; }
    int *cvlist = (int *)malloc((size_t)n * sizeof(int));
    int32_t *cflist = (int32_t *)malloc((nf > 0 ? nf : 1) * sizeof(int32_t));
    {
        long *vcur = (long *)malloc((size_t)ncomp * sizeof(long));
        long *fcur = (long *)malloc((size_t)ncomp * sizeof(long));
        memcpy(vcur, cvoff, (size_t)ncomp * sizeof(long));
        memcpy(fcur, cfoff, (size_t)ncomp * sizeof(long));
        for (int i = 0; i < n; i++) cvlist[vcur[cid[i]]++] = i;
        for (size_t f = 0; f < nf; f++) cflist[fcur[cid[faces[f*3+0]]]++] = (int32_t)f;
        free(vcur); free(fcur);
    }

    /* 1b) Per-component AREA-weighted, min-cos-gated radial votes (identical
     *     measure to obj_orient_audit): area_out/area_in over faces whose normal
     *     points away from / toward the scroll axis. This is the robust radial
     *     reference for the fallback below; only needed when an axis is given. */
    double *caout = (double *)calloc((size_t)ncomp, sizeof(double));
    double *cain  = (double *)calloc((size_t)ncomp, sizeof(double));
    if (have_axis) {
        for (size_t f = 0; f < nf; f++) {
            int i0 = faces[f*3+0], i1 = faces[f*3+1], i2 = faces[f*3+2];
            int c = cid[i0];
            const float *p0 = &verts[(size_t)i0*3], *p1 = &verts[(size_t)i1*3], *p2 = &verts[(size_t)i2*3];
            double e1[3], e2[3], gn[3];
            float cenf[3];
            double rhat[3], area = 0.0, cosr = 0.0;
            int k = 0;
            for (k = 0; k < 3; k++) {
                e1[k] = (double)p1[k] - p0[k];
                e2[k] = (double)p2[k] - p0[k];
                cenf[k] = (float)(((double)p0[k] + p1[k] + p2[k]) / 3.0);
            }
            gn[0] = e1[1]*e2[2] - e1[2]*e2[1];
            gn[1] = e1[2]*e2[0] - e1[0]*e2[2];
            gn[2] = e1[0]*e2[1] - e1[1]*e2[0];
            area = 0.5 * sqrt(gn[0]*gn[0] + gn[1]*gn[1] + gn[2]*gn[2]);
            if (area < 1e-12) continue;
            if (!radial_unit(cenf, axis_point, ad, rhat)) continue;  /* on-axis: skip */
            cosr = (gn[0]*rhat[0] + gn[1]*rhat[1] + gn[2]*rhat[2]) / (2.0 * area);
            if (cosr >= ORIENT_WELD_MIN_COS)      caout[c] += area;
            else if (cosr <= -ORIENT_WELD_MIN_COS) cain[c]  += area;
        }
    }

    /* 2) Area-weighted vertex normals from the current winding (so flipping a
     *    face's winding negates its contribution). Internally consistent
     *    components therefore get a coherent per-vertex normal field. */
    double *vn = (double *)calloc((size_t)n * 3, sizeof(double));
    for (size_t f = 0; f < nf; f++) {
        int i0 = faces[f*3+0], i1 = faces[f*3+1], i2 = faces[f*3+2];
        const float *p0 = &verts[(size_t)i0*3], *p1 = &verts[(size_t)i1*3], *p2 = &verts[(size_t)i2*3];
        double a0 = p1[0]-p0[0], a1 = p1[1]-p0[1], a2 = p1[2]-p0[2];
        double b0 = p2[0]-p0[0], b1 = p2[1]-p0[1], b2 = p2[2]-p0[2];
        double c0 = a1*b2 - a2*b1, c1 = a2*b0 - a0*b2, c2 = a0*b1 - a1*b0;  /* area-weighted */
        int idx[3] = { i0, i1, i2 };
        for (int k = 0; k < 3; k++) { vn[idx[k]*3+0]+=c0; vn[idx[k]*3+1]+=c1; vn[idx[k]*3+2]+=c2; }
    }
    for (int i = 0; i < n; i++) {
        double *p = &vn[i*3], L = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
        if (L > 1e-20) { p[0]/=L; p[1]/=L; p[2]/=L; }
    }

    /* 3) Uniform spatial grid over all verts (cell grows if the volume is huge). */
    double gmin[3] = {1e30,1e30,1e30}, gmax[3] = {-1e30,-1e30,-1e30};
    for (int i = 0; i < n; i++) for (int k = 0; k < 3; k++) {
        double x = verts[i*3+k]; if (x < gmin[k]) gmin[k] = x; if (x > gmax[k]) gmax[k] = x;
    }
    double cs = radius > 1e-6 ? radius : 1.0;
    int dim[3]; size_t ncell;
    for (;;) {
        for (int k = 0; k < 3; k++) { dim[k] = (int)((gmax[k]-gmin[k])/cs) + 1; if (dim[k] < 1) dim[k] = 1; }
        ncell = (size_t)dim[0]*(size_t)dim[1]*(size_t)dim[2];
        if (ncell <= 50000000u) break;
        cs *= 2.0;
    }
    int reach = (int)ceil((double)radius / cs); if (reach < 1) reach = 1;
    int *cellof = (int *)malloc((size_t)n * sizeof(int));
    long *coff = (long *)calloc(ncell + 1, sizeof(long));
    int planeYZ = dim[1]*dim[2];
    for (int i = 0; i < n; i++) {
        int g[3];
        for (int k = 0; k < 3; k++) { int q = (int)((verts[i*3+k]-gmin[k])/cs); if (q<0)q=0; if (q>=dim[k])q=dim[k]-1; g[k]=q; }
        int ic = (g[0]*dim[1] + g[1])*dim[2] + g[2];
        cellof[i] = ic; coff[ic+1]++;
    }
    for (size_t i = 0; i < ncell; i++) coff[i+1] += coff[i];
    int *cellv = (int *)malloc((size_t)coff[ncell] * sizeof(int));
    long *cur = (long *)malloc(ncell * sizeof(long));
    for (size_t i = 0; i < ncell; i++) cur[i] = coff[i];
    for (int i = 0; i < n; i++) cellv[cur[cellof[i]]++] = i;
    free(cur);

    /* 4) Components largest-first; the biggest is the trusted anchor. */
    CompOrder *ord = (CompOrder *)malloc((size_t)ncomp * sizeof(CompOrder));
    for (int c = 0; c < ncomp; c++) { ord[c].nf = cnf[c]; ord[c].id = c; }
    qsort(ord, (size_t)ncomp, sizeof(CompOrder), cmp_comporder);

    uint8_t *oriented = (uint8_t *)calloc((size_t)ncomp, 1);
    double r2 = (double)radius * radius;
    size_t flipped = 0, radial_decided = 0;

    /* Robust radial reference = the anchor (largest) component's area-weighted,
     * min-cos-gated radial sign. Everything is oriented to the anchor (spatial
     * chains to it; radial compares to this sign), so all decisions stay
     * mutually consistent. */
    int anchor_c = ord[0].id, radial_ref = 0, have_ref = 0;
    if (have_axis && (caout[anchor_c] + cain[anchor_c]) > 0.0) {
        radial_ref = (caout[anchor_c] >= cain[anchor_c]) ? 1 : -1;
        have_ref = 1;
    }

    for (int oi = 0; oi < ncomp; oi++) {
        int c = ord[oi].id;
        if (oi == 0) { oriented[c] = 1; continue; }   /* anchor: largest comp */

        /* Face-less components (collapse-orphaned verts) sort to the tail and
         * can be skipped outright: every face-bearing component was processed
         * before any of them (so they never served as a vote target), their
         * own vertex normals are zero (net == 0 -> never flip), and they have
         * no faces to flip. Output-identical to processing them. */
        if (ord[oi].nf == 0) break;

        /* Vote: sum dot(n_v, n_w) over near vertex pairs to already-oriented comps. */
        double net = 0.0;
        size_t pairs = 0, nvc = (size_t)(cvoff[c+1] - cvoff[c]);
        for (long vli = cvoff[c]; vli < cvoff[c+1]; vli++) {
            int v = cvlist[vli];
            int g0 = cellof[v]/planeYZ, rem = cellof[v]%planeYZ, g1 = rem/dim[2], g2 = rem%dim[2];
            const double *nvv = &vn[v*3];
            double vz = verts[v*3+0], vy = verts[v*3+1], vx = verts[v*3+2];
            for (int dz=-reach; dz<=reach; dz++) for (int dy=-reach; dy<=reach; dy++) for (int dx=-reach; dx<=reach; dx++) {
                int n0=g0+dz, n1=g1+dy, n2=g2+dx;
                if (n0<0||n1<0||n2<0||n0>=dim[0]||n1>=dim[1]||n2>=dim[2]) continue;
                size_t nc = ((size_t)n0*(size_t)dim[1]+(size_t)n1)*(size_t)dim[2]+(size_t)n2;
                for (long t = coff[nc]; t < coff[nc+1]; t++) {
                    int w = cellv[t];
                    if (cid[w] == c || !oriented[cid[w]]) continue;
                    double ddz=verts[w*3+0]-vz, ddy=verts[w*3+1]-vy, ddx=verts[w*3+2]-vx;
                    if (ddz*ddz+ddy*ddy+ddx*ddx > r2) continue;
                    const double *nw = &vn[w*3];
                    net += nvv[0]*nw[0] + nvv[1]*nw[1] + nvv[2]*nw[2];
                    pairs++;
                }
            }
        }
        int do_flip = 0;
        if (getenv("ORIENT_WELD_DEBUG") != NULL) {
            double out = caout[c], in = cain[c], tot = out + in;
            double agr = (tot > 0.0) ? (out > in ? out : in) / tot : 0.0;
            fprintf(stderr, "[orient_weld] rank=%d nf=%ld pairs=%zu net=%.2f "
                    "area_out=%.1f area_in=%.1f agr=%.2f ref=%d\n",
                    oi, ord[oi].nf, pairs, net, out, in, agr, radial_ref);
        }
        /* Spatial evidence is DECISIVE only with real mass: >=100 near pairs,
         * or average >=1 pair per vertex (tiny comps in genuine full contact).
         * A glancing touch must not pre-empt the radial vote -- on the real
         * scroll a 39k-vert detached shell collected 14 incidental pairs
         * (net +12) against a radial vote of -37k and stayed backward. */
        int spatial_decisive = (pairs >= 100)
                               || (nvc > 0 && pairs >= nvc);
        if (spatial_decisive || (pairs > 0 && !(have_axis && have_ref))) {
            do_flip = (net < 0.0);                     /* spatial vote, as before */
        } else if (have_axis && have_ref) {
            /* Weak or no spatial evidence: RADIAL fallback. Compare this
             * component's AREA-weighted, min-cos-gated radial facing (robust on
             * folded cores, unlike the old per-vertex sum) to the anchor's. */
            double out = caout[c], in = cain[c], tot = out + in;
            if (tot > 0.0) {
                double agr = (out > in ? out : in) / tot;
                int csign = (out >= in) ? 1 : -1;
                if (agr >= ORIENT_WELD_RADIAL_AGREE) {
                    do_flip = (csign != radial_ref);   /* decisive disagreement */
                    radial_decided++;
                } else if (pairs > 0) {
                    do_flip = (net < 0.0);   /* radially ambiguous: spatial */
                }
                /* else radially ambiguous + no contact: leave untouched */
            } else if (pairs > 0) {
                do_flip = (net < 0.0);       /* radial mute (debris): spatial */
            }
        }
        if (do_flip) {                                 /* opposes its neighbours -> flip */
            for (long fli = cfoff[c]; fli < cfoff[c+1]; fli++) {
                int32_t f = cflist[fli];
                int32_t tmp = faces[f*3+1]; faces[f*3+1] = faces[f*3+2]; faces[f*3+2] = tmp;
            }
            for (long vli = cvoff[c]; vli < cvoff[c+1]; vli++) {
                int v = cvlist[vli];
                vn[v*3+0]=-vn[v*3+0]; vn[v*3+1]=-vn[v*3+1]; vn[v*3+2]=-vn[v*3+2];
            }
            flipped++;
        }
        oriented[c] = 1;
    }

    if (out_flipped) *out_flipped = flipped;
    if (out_radial) *out_radial = radial_decided;
    free(cvoff); free(cfoff); free(cvlist); free(cflist);
    free(cid); free(cnf); free(caout); free(cain); free(vn); free(cellof); free(coff); free(cellv); free(ord); free(oriented);
    return 0;
}
