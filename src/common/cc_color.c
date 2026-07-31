/* cc_color.c -- connected-component vertex coloring (see cc_color.h).
 *
 * Extracted from the obj_cc_color tool so the pipeline stage dumps and grid_weld
 * can bake the coloring straight into their OBJs. Self-contained: libc + math
 * only. NOTE: the union-find helpers are file-static and PREFIXED (ccuf_*) on
 * purpose -- cube_mesh and grid_weld also link src/common/union_find.c, so a
 * plain uf_find/uf_union here would be a duplicate-symbol link error.
 */
#include "cc_color.h"

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

/* ---------- HSV -> RGB (h in [0,360), s,v in [0,1]) ---------- */
static void cc_hsv2rgb(double h, double s, double v, float *r, float *g, float *b)
{
    double c = v * s;
    double hp = h / 60.0;
    double x = c * (1.0 - fabs(fmod(hp, 2.0) - 1.0));
    double r1 = 0, g1 = 0, b1 = 0;
    if      (hp < 1) { r1 = c; g1 = x; }
    else if (hp < 2) { r1 = x; g1 = c; }
    else if (hp < 3) { g1 = c; b1 = x; }
    else if (hp < 4) { g1 = x; b1 = c; }
    else if (hp < 5) { r1 = x; b1 = c; }
    else             { r1 = c; b1 = x; }
    double m = v - c;
    *r = (float)(r1 + m); *g = (float)(g1 + m); *b = (float)(b1 + m);
}

/* ---------- union-find (path halving + union by size), prefixed ---------- */
static int32_t ccuf_find(int32_t *p, int32_t x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; }
static void ccuf_union(int32_t *p, int32_t *sz, int32_t a, int32_t b)
{
    a = ccuf_find(p, a); b = ccuf_find(p, b);
    if (a == b) return;
    if (sz[a] < sz[b]) { int32_t t = a; a = b; b = t; }
    p[b] = a; sz[a] += sz[b];
}

typedef struct { int32_t root; size_t faces; float rgb[3]; } CCComp;
static int cc_comp_cmp(const void *pa, const void *pb)
{
    const CCComp *a = (const CCComp *)pa, *b = (const CCComp *)pb;
    if (a->faces != b->faces) return a->faces < b->faces ? 1 : -1;  /* descending */
    return a->root < b->root ? -1 : 1;
}

void CCColor_default_opts(CCColorOpts *o)
{
    if (!o) return;
    o->sat = 0.62;
    o->min_faces = 0;
}

static void fill_grey(float *out_rgb, size_t nv)
{
    for (size_t i = 0; i < nv * 3; i++) out_rgb[i] = 0.5f;
}

size_t CCColor_compute(size_t nv, const int32_t *faces, size_t nf,
                       const CCColorOpts *opts, float *out_rgb, CCColorStats *stats)
{
    CCColorOpts d; if (!opts) { CCColor_default_opts(&d); opts = &d; }
    if (stats) { stats->ncomp = stats->nv = stats->nf = 0; stats->largest_faces = 0;
                 stats->cover50 = stats->cover90 = stats->cover99 = 0;
                 stats->dust_comps = stats->dust_faces = 0; }
    if (!out_rgb || nv == 0) return 0;
    if (!faces || nf == 0) { fill_grey(out_rgb, nv); return 0; }

    int32_t *par = (int32_t *)malloc(nv * sizeof(int32_t));
    int32_t *usz = (int32_t *)malloc(nv * sizeof(int32_t));
    if (!par || !usz) { free(par); free(usz); fill_grey(out_rgb, nv); return 0; }
    for (size_t i = 0; i < nv; i++) { par[i] = (int32_t)i; usz[i] = 1; }
    for (size_t f = 0; f < nf; f++) {
        int32_t a = faces[f*3+0], b = faces[f*3+1], c = faces[f*3+2];
        if (a < 0 || b < 0 || c < 0 || (size_t)a >= nv || (size_t)b >= nv || (size_t)c >= nv) continue;
        ccuf_union(par, usz, a, b); ccuf_union(par, usz, b, c);
    }

    size_t *root_faces = (size_t *)calloc(nv, sizeof(size_t));
    if (!root_faces) { free(par); free(usz); fill_grey(out_rgb, nv); return 0; }
    for (size_t f = 0; f < nf; f++) {
        int32_t a = faces[f*3+0];
        if (a < 0 || (size_t)a >= nv) continue;
        root_faces[ccuf_find(par, a)]++;
    }
    size_t ncomp = 0;
    for (size_t i = 0; i < nv; i++) if (root_faces[i] > 0) ncomp++;

    CCComp *comp = (CCComp *)malloc((ncomp ? ncomp : 1) * sizeof(CCComp));
    int32_t *root2rank = (int32_t *)malloc(nv * sizeof(int32_t));
    if (!comp || !root2rank) { free(comp); free(root2rank); free(root_faces); free(par); free(usz);
                               fill_grey(out_rgb, nv); return 0; }
    size_t ci = 0;
    for (size_t i = 0; i < nv; i++) if (root_faces[i] > 0) { comp[ci].root = (int32_t)i; comp[ci].faces = root_faces[i]; ci++; }
    qsort(comp, ncomp, sizeof(CCComp), cc_comp_cmp);

    for (size_t i = 0; i < nv; i++) root2rank[i] = -1;
    for (size_t r = 0; r < ncomp; r++) {
        root2rank[comp[r].root] = (int32_t)r;
        if (comp[r].faces < opts->min_faces) {
            comp[r].rgb[0] = comp[r].rgb[1] = comp[r].rgb[2] = 0.28f;   /* dust -> dim grey */
        } else {
            double hue = fmod((double)r * 0.61803398875, 1.0) * 360.0;  /* golden angle */
            double val = (r % 2) ? 0.98 : 0.80;                         /* separate collisions */
            cc_hsv2rgb(hue, opts->sat, val, &comp[r].rgb[0], &comp[r].rgb[1], &comp[r].rgb[2]);
        }
    }

    for (size_t i = 0; i < nv; i++) {
        int32_t rk = root2rank[ccuf_find(par, (int32_t)i)];
        if (rk < 0) { out_rgb[i*3+0] = out_rgb[i*3+1] = out_rgb[i*3+2] = 0.5f; }
        else { out_rgb[i*3+0] = comp[rk].rgb[0]; out_rgb[i*3+1] = comp[rk].rgb[1]; out_rgb[i*3+2] = comp[rk].rgb[2]; }
    }

    if (stats) {
        stats->ncomp = ncomp; stats->nv = nv; stats->nf = nf;
        stats->largest_faces = ncomp ? comp[0].faces : 0;
        size_t cum = 0;
        for (size_t r = 0; r < ncomp; r++) {
            cum += comp[r].faces;
            if (!stats->cover50 && cum * 100 >= nf * 50) stats->cover50 = r + 1;
            if (!stats->cover90 && cum * 100 >= nf * 90) stats->cover90 = r + 1;
            if (!stats->cover99 && cum * 100 >= nf * 99) stats->cover99 = r + 1;
            if (comp[r].faces < opts->min_faces) { stats->dust_comps++; stats->dust_faces += comp[r].faces; }
        }
    }

    free(comp); free(root2rank); free(root_faces); free(par); free(usz);
    return ncomp;
}

/* ---------- selftest ---------- */
int CCColor_selftest(void)
{
    int fails = 0;
    /* verts 0..7; tri A=(0,1,2), tri B=(3,4,5), tri C=(2,6,7) shares vertex 2 with
     * A (a bowtie), so {A,C}=one component (verts 0,1,2,6,7) and B=another. -> 2. */
    int32_t faces[] = { 0,1,2,  3,4,5,  2,6,7 };
    float rgb[8*3];
    CCColorStats st;
    size_t nc = CCColor_compute(8, faces, 3, NULL, rgb, &st);
    if (nc != 2)          { fprintf(stderr, "  [cc_color] selftest: %zu comps (want 2)\n", nc); fails++; }
    if (st.ncomp != 2)    { fprintf(stderr, "  [cc_color] selftest: stats ncomp=%zu (want 2)\n", st.ncomp); fails++; }
    /* the bowtie-joined verts share one color; the disjoint triangle differs */
    int same_A = (rgb[0*3+0]==rgb[2*3+0] && rgb[0*3+1]==rgb[2*3+1] && rgb[0*3+2]==rgb[2*3+2]) &&
                 (rgb[2*3+0]==rgb[6*3+0] && rgb[2*3+1]==rgb[6*3+1] && rgb[2*3+2]==rgb[6*3+2]);
    int diff_B = !(rgb[0*3+0]==rgb[3*3+0] && rgb[0*3+1]==rgb[3*3+1] && rgb[0*3+2]==rgb[3*3+2]);
    if (!same_A) { fprintf(stderr, "  [cc_color] selftest: bowtie verts not same color\n"); fails++; }
    if (!diff_B) { fprintf(stderr, "  [cc_color] selftest: disjoint comp shares color\n"); fails++; }
    /* empty / degenerate inputs must not crash and must grey-fill */
    { float g[3*3]; size_t z = CCColor_compute(3, NULL, 0, NULL, g, NULL);
      if (z != 0 || g[0] != 0.5f) { fprintf(stderr, "  [cc_color] selftest: empty-faces not greyed\n"); fails++; } }
    fprintf(stderr, "cc_color selftest: %s\n", fails ? "FAIL" : "PASS");
    return fails;
}
