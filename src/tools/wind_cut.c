/* wind_cut.c -- winding-driven fusion cutter for scroll meshes.
 *
 * The CVT remesher is winding-BLIND: with CVT_INTERIOR_EDGE=13 vox > pitch 9.5
 * it triangulates ACROSS the inter-wrap gap, stitching wrap N to wrap N+1 into
 * one giant component (the "red monster": 4x5x5 comp #0 = 54% of faces, spans
 * 36.6 wrap-pitches radially, 2897 full-turn short-circuits). BPA and the weld
 * both gate cross-wrap edges out; CVT puts them back.
 *
 * The winding coordinate w = r/pitch - theta/2pi is CONSTANT along a wrap and
 * jumps by ~1 across the gap, so an edge with |dw| ~ 1 IS a radial short-circuit.
 * This tool removes every face carrying such an edge -- branch-cut-free (local
 * pairwise dw), so it works even at the r=12 core where no edge-length threshold
 * can. The fused component falls apart into its per-wrap sheets along the cut.
 *
 *   wind_cut <in.obj> <out.obj> --umb-y Y --umb-x X [--pitch P=9.5] [--tol T=0.7]
 *            [--min-comp-faces M=0]
 *   wind_cut --selftest
 * Exit: 0 ok / selftest pass, 1 IO, 2 usage, 3 selftest fail.
 */
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define TWO_PI (2.0 * M_PI)

/* ---------- growing arrays ---------- */
typedef struct { float *v; size_t n, cap; } FVec;
typedef struct { int32_t *f; size_t n, cap; } IVec;

static int fv_push(FVec *a, float x, float y, float z)
{
    if (a->n == a->cap) { size_t nc = a->cap ? a->cap * 2 : 1u << 16;
        float *np = (float *)realloc(a->v, nc * 3 * sizeof(float)); if (!np) return -1; a->v = np; a->cap = nc; }
    a->v[a->n*3+0] = x; a->v[a->n*3+1] = y; a->v[a->n*3+2] = z; a->n++; return 0;
}
static int iv_push(IVec *a, int32_t x, int32_t y, int32_t z)
{
    if (a->n == a->cap) { size_t nc = a->cap ? a->cap * 2 : 1u << 16;
        int32_t *np = (int32_t *)realloc(a->f, nc * 3 * sizeof(int32_t)); if (!np) return -1; a->f = np; a->cap = nc; }
    a->f[a->n*3+0] = x; a->f[a->n*3+1] = y; a->f[a->n*3+2] = z; a->n++; return 0;
}

static int read_obj(const char *path, FVec *V, IVec *F)
{
    FILE *fp = fopen(path, "rb");
    if (!fp) { fprintf(stderr, "wind_cut: cannot open %s\n", path); return -1; }
    size_t lcap = 1u << 16; char *line = (char *)malloc(lcap);
    int rc = 0;
    if (!line) { fclose(fp); return -1; }
    while (fgets(line, (int)lcap, fp)) {
        while (!strchr(line, '\n') && !feof(fp)) {
            size_t len = strlen(line);
            char *nl = (char *)realloc(line, lcap * 2); if (!nl) { rc = -1; goto done; }
            line = nl; lcap *= 2;
            if (!fgets(line + len, (int)(lcap - len), fp)) break;
        }
        if (line[0] == 'v' && (line[1] == ' ' || line[1] == '\t')) {
            double x = 0, y = 0, z = 0;
            if (sscanf(line + 2, "%lf %lf %lf", &x, &y, &z) == 3)
                if (fv_push(V, (float)x, (float)y, (float)z) != 0) { rc = -1; goto done; }
        } else if (line[0] == 'f' && (line[1] == ' ' || line[1] == '\t')) {
            char *tok = strtok(line + 2, " \t\r\n");
            int32_t idx[3]; int n = 0;
            while (tok && n < 3) { idx[n++] = (int32_t)atol(tok); tok = strtok(NULL, " \t\r\n"); }
            if (n == 3) {
                int32_t a = idx[0], b = idx[1], c = idx[2];
                if (a < 0) a = (int32_t)V->n + a + 1;
                if (b < 0) b = (int32_t)V->n + b + 1;
                if (c < 0) c = (int32_t)V->n + c + 1;
                if (iv_push(F, a - 1, b - 1, c - 1) != 0) { rc = -1; goto done; }
            }
        }
    }
done:
    free(line); fclose(fp);
    return rc;
}

static int write_obj(const char *path, const FVec *V, const int32_t *F, size_t nf)
{
    FILE *fp = fopen(path, "wb");
    if (!fp) { fprintf(stderr, "wind_cut: cannot write %s\n", path); return -1; }
    for (size_t i = 0; i < V->n; i++)
        fprintf(fp, "v %.6f %.6f %.6f\n", V->v[i*3+0], V->v[i*3+1], V->v[i*3+2]);
    for (size_t i = 0; i < nf; i++)
        fprintf(fp, "f %d %d %d\n", F[i*3+0]+1, F[i*3+1]+1, F[i*3+2]+1);
    fclose(fp);
    return 0;
}

/* ---------- union-find ---------- */
static int32_t uf_find(int32_t *p, int32_t x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; }
static void uf_union(int32_t *p, int32_t *sz, int32_t a, int32_t b)
{
    a = uf_find(p, a); b = uf_find(p, b);
    if (a == b) return;
    if (sz[a] < sz[b]) { int32_t t = a; a = b; b = t; }
    p[b] = a; sz[a] += sz[b];
}

/* ---------- winding (coord order z,y,x; umbilicus in y,x) ---------- */
static double wrap_pmpi(double a) { while (a > M_PI) a -= TWO_PI; while (a <= -M_PI) a += TWO_PI; return a; }
static double vradius(const float *V, size_t i, double uy, double ux)
{ double dy = (double)V[i*3+1] - uy, dx = (double)V[i*3+2] - ux; return hypot(dy, dx); }
static double dwind(const float *V, size_t u, size_t v, double uy, double ux, double pitch)
{
    double dyu = (double)V[u*3+1] - uy, dxu = (double)V[u*3+2] - ux;
    double dyv = (double)V[v*3+1] - uy, dxv = (double)V[v*3+2] - ux;
    double dr  = hypot(dyv, dxv) - hypot(dyu, dxu);
    double dth = wrap_pmpi(atan2(dyv, dxv) - atan2(dyu, dxu));
    return dr / pitch - dth / TWO_PI;
}

/* max |dw| over a face's three edges */
static double face_max_dw(const float *V, const int32_t *f, double uy, double ux, double pitch)
{
    double d0 = fabs(dwind(V, (size_t)f[0], (size_t)f[1], uy, ux, pitch));
    double d1 = fabs(dwind(V, (size_t)f[1], (size_t)f[2], uy, ux, pitch));
    double d2 = fabs(dwind(V, (size_t)f[2], (size_t)f[0], uy, ux, pitch));
    double m = d0 > d1 ? d0 : d1; return m > d2 ? m : d2;
}

/* Component report over an active-face mask: n comps (>= min faces), biggest
 * comp face count + radial turn-span (r-span/pitch). */
typedef struct { size_t n_comp; size_t big_faces; double big_turns; } CompStat;
static CompStat comp_stats(const FVec *V, const int32_t *F, const uint8_t *active, size_t nf,
                           double uy, double ux, double pitch, size_t min_faces)
{
    CompStat cs = {0,0,0.0};
    size_t nv = V->n;
    int32_t *p = (int32_t *)malloc(nv * sizeof(int32_t));
    int32_t *sz = (int32_t *)malloc(nv * sizeof(int32_t));
    if (!p || !sz) { free(p); free(sz); return cs; }
    for (size_t i = 0; i < nv; i++) { p[i] = (int32_t)i; sz[i] = 1; }
    for (size_t i = 0; i < nf; i++) {
        if (!active[i]) continue;
        uf_union(p, sz, F[i*3+0], F[i*3+1]);
        uf_union(p, sz, F[i*3+1], F[i*3+2]);
    }
    /* per-root: face count + rmin/rmax */
    size_t *fc = (size_t *)calloc(nv, sizeof(size_t));
    double *rmin = (double *)malloc(nv * sizeof(double));
    double *rmax = (double *)malloc(nv * sizeof(double));
    if (!fc || !rmin || !rmax) { free(p); free(sz); free(fc); free(rmin); free(rmax); return cs; }
    for (size_t i = 0; i < nv; i++) { rmin[i] = 1e30; rmax[i] = -1e30; }
    for (size_t i = 0; i < nf; i++) {
        if (!active[i]) continue;
        int32_t r = uf_find(p, F[i*3+0]);
        fc[r]++;
        for (int k = 0; k < 3; k++) {
            double rr = vradius(V->v, (size_t)F[i*3+k], uy, ux);
            if (rr < rmin[r]) rmin[r] = rr;
            if (rr > rmax[r]) rmax[r] = rr;
        }
    }
    for (size_t i = 0; i < nv; i++) {
        if (fc[i] < min_faces || fc[i] == 0) continue;
        cs.n_comp++;
        if (fc[i] > cs.big_faces) { cs.big_faces = fc[i]; cs.big_turns = (rmax[i] - rmin[i]) / pitch; }
    }
    free(p); free(sz); free(fc); free(rmin); free(rmax);
    return cs;
}

/* ---------- cut ---------- */
static int run(const char *in, const char *out, double uy, double ux, double pitch,
               double tol, size_t min_comp_faces)
{
    FVec V = {0}; IVec F = {0};
    if (read_obj(in, &V, &F) != 0) { free(V.v); free(F.f); return 1; }
    if (V.n == 0 || F.n == 0) { fprintf(stderr, "wind_cut: empty mesh %s\n", in); free(V.v); free(F.f); return 1; }

    uint8_t *keep = (uint8_t *)malloc(F.n);
    uint8_t *all  = (uint8_t *)malloc(F.n);
    if (!keep || !all) { free(V.v); free(F.f); free(keep); free(all); return 1; }
    size_t n_cut = 0;
    for (size_t i = 0; i < F.n; i++) {
        all[i] = 1;
        int cross = face_max_dw(V.v, &F.f[i*3], uy, ux, pitch) > tol;
        keep[i] = cross ? 0 : 1;
        if (cross) n_cut++;
    }

    CompStat before = comp_stats(&V, F.f, all, F.n, uy, ux, pitch, 1);

    /* compact kept faces */
    int32_t *KF = (int32_t *)malloc(F.n * 3 * sizeof(int32_t));
    if (!KF) { free(V.v); free(F.f); free(keep); free(all); return 1; }
    size_t nkf = 0;
    for (size_t i = 0; i < F.n; i++)
        if (keep[i]) { KF[nkf*3+0]=F.f[i*3+0]; KF[nkf*3+1]=F.f[i*3+1]; KF[nkf*3+2]=F.f[i*3+2]; nkf++; }

    /* optional: drop small components (dust from shredding) */
    size_t n_drop_faces = 0, n_drop_comp = 0;
    if (min_comp_faces > 0 && nkf > 0) {
        size_t nv = V.n;
        int32_t *p = (int32_t *)malloc(nv * sizeof(int32_t));
        int32_t *sz = (int32_t *)malloc(nv * sizeof(int32_t));
        size_t *fc = (size_t *)calloc(nv, sizeof(size_t));
        if (p && sz && fc) {
            for (size_t i = 0; i < nv; i++) { p[i]=(int32_t)i; sz[i]=1; }
            for (size_t i = 0; i < nkf; i++) { uf_union(p,sz,KF[i*3+0],KF[i*3+1]); uf_union(p,sz,KF[i*3+1],KF[i*3+2]); }
            for (size_t i = 0; i < nkf; i++) fc[uf_find(p,KF[i*3+0])]++;
            size_t w = 0;
            for (size_t i = 0; i < nkf; i++) {
                if (fc[uf_find(p,KF[i*3+0])] >= min_comp_faces) {
                    KF[w*3+0]=KF[i*3+0]; KF[w*3+1]=KF[i*3+1]; KF[w*3+2]=KF[i*3+2]; w++;
                } else n_drop_faces++;
            }
            for (size_t i = 0; i < nv; i++) if (fc[i] > 0 && fc[i] < min_comp_faces) n_drop_comp++;
            nkf = w;
        }
        free(p); free(sz); free(fc);
    }

    /* after-stats over the surviving faces */
    uint8_t *kact = (uint8_t *)malloc(nkf ? nkf : 1);
    for (size_t i = 0; i < nkf; i++) kact[i] = 1;
    CompStat after = comp_stats(&V, KF, kact, nkf, uy, ux, pitch, 1);
    free(kact);

    int rc = write_obj(out, &V, KF, nkf) == 0 ? 0 : 1;

    printf("wind_cut: %s -> %s\n", in, out);
    printf("  umbilicus=(y %.1f, x %.1f) pitch=%.2f  cut tol |dw|>%.2f\n", uy, ux, pitch, tol);
    printf("  faces: in=%zu  cut(cross-wrap)=%zu (%.1f%%)  dropped(small comp)=%zu  out=%zu\n",
           F.n, n_cut, 100.0*(double)n_cut/(double)F.n, n_drop_faces, nkf);
    printf("  components: %zu -> %zu (+%zd)\n", before.n_comp, after.n_comp,
           (ptrdiff_t)after.n_comp - (ptrdiff_t)before.n_comp);
    printf("  biggest component: %zu f spanning %.1f turns  ->  %zu f spanning %.1f turns\n",
           before.big_faces, before.big_turns, after.big_faces, after.big_turns);
    if (min_comp_faces > 0) printf("  (dropped %zu small comps < %zu faces)\n", n_drop_comp, min_comp_faces);

    free(V.v); free(F.f); free(keep); free(all); free(KF);
    return rc;
}

/* ---------- selftest: two concentric rings fused by one bridge -> cut -> 2 comps ---------- */
static int selftest(void)
{
    /* Build two concentric octagon rings (r=10 and r=20 about origin in y,x),
     * each a closed strip of quads (triangulated), pitch=10 so they differ by
     * dw~1. Add ONE bridge quad linking a vertex of the inner ring to the outer
     * ring (the fusion). Cutting |dw|>0.5 must remove the bridge -> 2 comps. */
    int fails = 0;
    FVec V = {0}; IVec F = {0};
    const int NSEG = 8; const double zlo = 0.0, zhi = 4.0;
    /* inner ring verts: 0..2*NSEG-1 (bottom,top per seg); outer: next block */
    for (int ring = 0; ring < 2; ring++) {
        double r = ring == 0 ? 10.0 : 20.0;
        for (int s = 0; s < NSEG; s++) {
            double a = TWO_PI * s / NSEG;
            float y = (float)(r*cos(a)), x = (float)(r*sin(a));
            fv_push(&V, (float)zlo, y, x);
            fv_push(&V, (float)zhi, y, x);
        }
    }
    /* strips */
    for (int ring = 0; ring < 2; ring++) {
        int base = ring * 2 * NSEG;
        for (int s = 0; s < NSEG; s++) {
            int s2 = (s + 1) % NSEG;
            int a = base + 2*s, b = base + 2*s + 1, c = base + 2*s2, d = base + 2*s2 + 1;
            iv_push(&F, a, b, d); iv_push(&F, a, d, c);
        }
    }
    size_t nf_clean = F.n;
    /* one bridge quad: inner seg0 (verts 0,1) to outer seg0 (verts 2*NSEG, 2*NSEG+1) */
    int io = 2*NSEG;
    iv_push(&F, 0, 1, io+1); iv_push(&F, 0, io+1, io);

    double uy = 0.0, ux = 0.0, pitch = 10.0, tol = 0.5;
    /* count comps before (all faces) */
    uint8_t *all = (uint8_t *)malloc(F.n); for (size_t i=0;i<F.n;i++) all[i]=1;
    CompStat before = comp_stats(&V, F.f, all, F.n, uy, ux, pitch, 1);
    /* cut */
    size_t n_cut = 0; uint8_t *keep = (uint8_t *)malloc(F.n);
    for (size_t i = 0; i < F.n; i++) { int cr = face_max_dw(V.v,&F.f[i*3],uy,ux,pitch)>tol; keep[i]=!cr; if(cr)n_cut++; }
    int32_t *KF = (int32_t*)malloc(F.n*3*sizeof(int32_t)); size_t nkf=0;
    for (size_t i=0;i<F.n;i++) if(keep[i]){KF[nkf*3+0]=F.f[i*3+0];KF[nkf*3+1]=F.f[i*3+1];KF[nkf*3+2]=F.f[i*3+2];nkf++;}
    uint8_t *ka=(uint8_t*)malloc(nkf); for(size_t i=0;i<nkf;i++)ka[i]=1;
    CompStat after = comp_stats(&V, KF, ka, nkf, uy, ux, pitch, 1);

    printf("  [selftest] before: comps=%zu (expect 1, bridge fuses the rings)\n", before.n_comp);
    printf("  [selftest] cut %zu bridge face(s) (expect 2)\n", n_cut);
    printf("  [selftest] after:  comps=%zu (expect 2, rings separated)\n", after.n_comp);
    printf("  [selftest] kept clean strip faces: %zu / %zu\n", nkf, nf_clean);
    if (before.n_comp != 1) { printf("  FAIL: rings not fused before cut\n"); fails++; }
    if (n_cut != 2)         { printf("  FAIL: expected 2 bridge faces cut, got %zu\n", n_cut); fails++; }
    if (after.n_comp != 2)  { printf("  FAIL: rings not separated after cut\n"); fails++; }
    if (nkf != nf_clean)    { printf("  FAIL: clean strip faces changed (%zu != %zu)\n", nkf, nf_clean); fails++; }

    free(all); free(keep); free(KF); free(ka); free(V.v); free(F.f);
    printf("wind_cut selftest: %s\n", fails ? "FAIL" : "PASS");
    return fails ? 3 : 0;
}

int main(int argc, char **argv)
{
    if (argc == 2 && !strcmp(argv[1], "--selftest")) return selftest();
    if (argc < 3) {
        fprintf(stderr,
          "usage: %s <in.obj> <out.obj> --umb-y Y --umb-x X [--pitch P=9.5] [--tol T=0.7]\n"
          "         [--min-comp-faces M=0]\n"
          "       %s --selftest\n", argv[0], argv[0]);
        return 2;
    }
    const char *in = argv[1], *out = argv[2];
    double uy = 0, ux = 0, pitch = 9.5, tol = 0.7; size_t minc = 0;
    int have_y = 0, have_x = 0;
    for (int i = 3; i < argc; i++) {
        if      (!strcmp(argv[i], "--umb-y") && i+1 < argc) { uy = atof(argv[++i]); have_y = 1; }
        else if (!strcmp(argv[i], "--umb-x") && i+1 < argc) { ux = atof(argv[++i]); have_x = 1; }
        else if (!strcmp(argv[i], "--pitch") && i+1 < argc) pitch = atof(argv[++i]);
        else if (!strcmp(argv[i], "--tol") && i+1 < argc) tol = atof(argv[++i]);
        else if (!strcmp(argv[i], "--min-comp-faces") && i+1 < argc) minc = (size_t)atoll(argv[++i]);
        else { fprintf(stderr, "wind_cut: unknown arg %s\n", argv[i]); return 2; }
    }
    if (!have_y || !have_x) { fprintf(stderr, "wind_cut: --umb-y and --umb-x are required\n"); return 2; }
    if (pitch <= 0 || tol <= 0) { fprintf(stderr, "wind_cut: pitch/tol must be > 0\n"); return 2; }
    return run(in, out, uy, ux, pitch, tol, minc);
}
