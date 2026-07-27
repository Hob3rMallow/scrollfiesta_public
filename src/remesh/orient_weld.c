/*
 * orient_weld.c -- post-weld cross-component orientation. See orient_weld.h.
 *
 * Internal scratch uses malloc/free (local to one call, same convention as
 * ball_pivot.c / seam_weld.c); the arena parameter is unused.
 */
#include "orient_weld.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int uf_find(int *uf, int x) { while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; } return x; }

typedef struct { long nf; int id; } CompOrder;
static int cmp_comporder(const void *pa, const void *pb)
{
    const CompOrder *a = (const CompOrder *)pa, *b = (const CompOrder *)pb;
    if (a->nf != b->nf) return a->nf < b->nf ? 1 : -1;   /* face count descending */
    return a->id < b->id ? -1 : (a->id > b->id);          /* tie-break: id ascending */
}

int OrientWeld_components(Arena_T arena,
                          const float *verts, size_t nv,
                          int32_t *faces, size_t nf,
                          float radius,
                          size_t *out_flipped)
{
    (void)arena;
    if (out_flipped) *out_flipped = 0;
    if (nf == 0 || nv == 0) return 0;
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
    size_t flipped = 0;
    for (int oi = 0; oi < ncomp; oi++) {
        int c = ord[oi].id;
        if (oi == 0) { oriented[c] = 1; continue; }   /* anchor: largest component */

        /* Vote: sum dot(n_v, n_w) over near vertex pairs to already-oriented comps. */
        double net = 0.0;
        for (int v = 0; v < n; v++) {
            if (cid[v] != c) continue;
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
                }
            }
        }
        if (net < 0.0) {                               /* opposes its neighbours -> flip */
            for (size_t f = 0; f < nf; f++) {
                if (cid[faces[f*3+0]] != c) continue;
                int32_t tmp = faces[f*3+1]; faces[f*3+1] = faces[f*3+2]; faces[f*3+2] = tmp;
            }
            for (int v = 0; v < n; v++) if (cid[v] == c) { vn[v*3+0]=-vn[v*3+0]; vn[v*3+1]=-vn[v*3+1]; vn[v*3+2]=-vn[v*3+2]; }
            flipped++;
        }
        oriented[c] = 1;
    }

    if (out_flipped) *out_flipped = flipped;
    free(cid); free(cnf); free(vn); free(cellof); free(coff); free(cellv); free(ord); free(oriented);
    return 0;
}
