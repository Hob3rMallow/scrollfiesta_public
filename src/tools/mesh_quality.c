/* mesh_quality — triangle-quality + valence report for an OBJ mesh.
 *
 *   mesh_quality <mesh.obj>
 *   mesh_quality --selftest
 *
 * Reports the metrics that separate a CVT/RVD remesh from a QEM decimation:
 *   - min-angle histogram (the sliver signature): worst angle, mean/median of
 *     per-triangle minimum angle, and the fraction of triangles below 10/20/30 deg.
 *   - radius-ratio q = 2*r_in/r_circ in [0,1] (1 = equilateral, ->0 = sliver): mean.
 *   - interior valence-6 fraction (CVT drives interior vertices toward valence 6).
 *
 * Geometry only; coordinate order (z,y,x vs x,y,z) is irrelevant to angles/lengths.
 * Companion to manifold_check (which covers topology). */
#include "../common/obj_io.h"
#include "../common/arena.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Per-triangle minimum interior angle (deg) and radius ratio q=2 r_in/r_circ.
 * Degenerate (zero-area / zero-length edge) triangles report min_ang=0, q=0. */
static void tri_metrics(const float *V, int32_t a, int32_t b, int32_t c,
                        double *min_ang_deg, double *q) {
    const float *A = &V[(size_t)a*3], *B = &V[(size_t)b*3], *C = &V[(size_t)c*3];
    double ab[3] = { B[0]-A[0], B[1]-A[1], B[2]-A[2] };
    double bc[3] = { C[0]-B[0], C[1]-B[1], C[2]-B[2] };
    double ca[3] = { A[0]-C[0], A[1]-C[1], A[2]-C[2] };
    double la = sqrt(bc[0]*bc[0]+bc[1]*bc[1]+bc[2]*bc[2]);  /* opposite A */
    double lb = sqrt(ca[0]*ca[0]+ca[1]*ca[1]+ca[2]*ca[2]);  /* opposite B */
    double lc = sqrt(ab[0]*ab[0]+ab[1]*ab[1]+ab[2]*ab[2]);  /* opposite C */
    if (la <= 0.0 || lb <= 0.0 || lc <= 0.0) { *min_ang_deg = 0.0; *q = 0.0; return; }
    /* angle at each vertex via law of cosines, clamped for float safety */
    double cA = (lb*lb + lc*lc - la*la) / (2.0*lb*lc);
    double cB = (la*la + lc*lc - lb*lb) / (2.0*la*lc);
    double cC = (la*la + lb*lb - lc*lc) / (2.0*la*lb);
    if (cA < -1.0) cA = -1.0; if (cA > 1.0) cA = 1.0;
    if (cB < -1.0) cB = -1.0; if (cB > 1.0) cB = 1.0;
    if (cC < -1.0) cC = -1.0; if (cC > 1.0) cC = 1.0;
    double angA = acos(cA), angB = acos(cB), angC = acos(cC);
    double mn = angA < angB ? (angA < angC ? angA : angC) : (angB < angC ? angB : angC);
    *min_ang_deg = mn * 180.0 / M_PI;
    /* radius ratio: r_in = area/s, r_circ = (la*lb*lc)/(4 area) => q = 8 area^2 / (s la lb lc) */
    double s = 0.5*(la+lb+lc);
    double area2 = s*(s-la)*(s-lb)*(s-lc);
    if (area2 <= 0.0) { *q = 0.0; return; }
    double area = sqrt(area2);
    double denom = s * la * lb * lc;
    *q = denom > 0.0 ? (8.0 * area * area) / denom : 0.0;
}

/* qsort comparator for doubles (ascending) */
static int cmp_dbl(const void *pa, const void *pb) {
    double a = *(const double *)pa, b = *(const double *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

/* Report the quality of a mesh already in memory. Returns 0. */
static int report(const char *name, const float *V, size_t nv,
                  const int32_t *F, size_t nf, Arena_T arena) {
    if (nf == 0) { printf("=== mesh_quality: %s ===\n  (no faces)\n", name); return 0; }

    double *mins = (double *)ARENA_ALLOC(arena, (long)(nf*sizeof(double)));
    int hist[6] = {0,0,0,0,0,0};   /* [0,10)[10,20)[20,30)[30,40)[40,50)[50,60] */
    int n_degen = 0;
    double worst = 180.0, sum_min = 0.0, sum_q = 0.0;
    size_t lt10 = 0, lt20 = 0, lt30 = 0;
    for (size_t t = 0; t < nf; t++) {
        double ma, q;
        tri_metrics(V, F[t*3+0], F[t*3+1], F[t*3+2], &ma, &q);
        mins[t] = ma; sum_min += ma; sum_q += q;
        if (ma <= 0.0) n_degen++;
        if (ma < worst) worst = ma;
        if (ma < 10.0) lt10++;
        if (ma < 20.0) lt20++;
        if (ma < 30.0) lt30++;
        int b = (int)(ma / 10.0); if (b > 5) b = 5; if (b < 0) b = 0;
        hist[b]++;
    }
    qsort(mins, nf, sizeof(double), cmp_dbl);
    double median = mins[nf/2];

    /* Vertex valence + boundary detection via a directed-edge count.
     * Build undirected edge multiplicity: an edge used by exactly 1 triangle is
     * a boundary edge; its endpoints are boundary vertices. Valence = # distinct
     * neighbours; we approximate with incident undirected-edge count via a hash. */
    size_t *val = (size_t *)ARENA_CALLOC(arena, (long)nv, (long)sizeof(size_t));
    unsigned char *is_bnd = (unsigned char *)ARENA_CALLOC(arena, (long)nv, (long)sizeof(unsigned char));
    /* edge hash: open-addressing on (min<<32|max) -> use count table sized 2*3*nf */
    size_t ecap = 1; while (ecap < nf*6 + 16) ecap <<= 1;
    int64_t *ekey = (int64_t *)ARENA_ALLOC(arena, (long)(ecap*sizeof(int64_t)));
    int32_t *ecnt = (int32_t *)ARENA_CALLOC(arena, (long)ecap, (long)sizeof(int32_t));
    for (size_t i = 0; i < ecap; i++) ekey[i] = -1;
    for (size_t t = 0; t < nf; t++) {
        int32_t v[3] = { F[t*3+0], F[t*3+1], F[t*3+2] };
        for (int e = 0; e < 3; e++) {
            int32_t p = v[e], qy = v[(e+1)%3];
            int32_t lo = p < qy ? p : qy, hi = p < qy ? qy : p;
            int64_t k = ((int64_t)lo << 32) | (uint32_t)hi;
            size_t h = (size_t)((uint64_t)k * 1099511628211ULL) & (ecap-1);
            while (ekey[h] != -1 && ekey[h] != k) h = (h+1) & (ecap-1);
            if (ekey[h] == -1) { ekey[h] = k; val[lo]++; val[hi]++; }
            ecnt[h]++;
        }
    }
    /* mark boundary vertices from multiplicity-1 edges */
    for (size_t h = 0; h < ecap; h++) {
        if (ekey[h] == -1 || ecnt[h] != 1) continue;
        int32_t lo = (int32_t)(ekey[h] >> 32), hi = (int32_t)(uint32_t)ekey[h];
        is_bnd[lo] = 1; is_bnd[hi] = 1;
    }
    size_t n_int = 0, n_int_v6 = 0;
    for (size_t i = 0; i < nv; i++) {
        if (is_bnd[i] || val[i] == 0) continue;
        n_int++;
        if (val[i] == 6) n_int_v6++;
    }

    printf("=== mesh_quality: %s ===\n", name);
    printf("  V=%zu  F=%zu  degenerate=%d\n", nv, nf, n_degen);
    printf("  min-angle: worst=%.2f  mean=%.2f  median=%.2f (deg)\n",
           worst, sum_min/(double)nf, median);
    printf("  slivers:   <10deg=%.2f%%  <20deg=%.2f%%  <30deg=%.2f%%\n",
           100.0*(double)lt10/(double)nf, 100.0*(double)lt20/(double)nf,
           100.0*(double)lt30/(double)nf);
    printf("  radius-ratio q (1=equilateral): mean=%.3f\n", sum_q/(double)nf);
    printf("  min-angle histogram (deg):\n");
    const char *lbl[6] = {" 0-10","10-20","20-30","30-40","40-50","50-60"};
    for (int b = 0; b < 6; b++)
        printf("    %s: %7d  (%.1f%%)\n", lbl[b], hist[b],
               100.0*(double)hist[b]/(double)nf);
    if (n_int > 0)
        printf("  interior valence-6: %.1f%% (%zu/%zu interior verts)\n",
               100.0*(double)n_int_v6/(double)n_int, n_int_v6, n_int);
    return 0;
}

static int selftest(void) {
    Arena_T a = Arena_new();
    int fails = 0;
    /* equilateral triangle -> min angle 60, q=1 */
    {
        float V[9] = { 0,0,0,  0,1,0,  0,0.5f,0.8660254f };
        int32_t F[3] = {0,1,2};
        double ma, q; tri_metrics(V, F[0], F[1], F[2], &ma, &q);
        if (fabs(ma - 60.0) > 1e-3) { fprintf(stderr, "selftest: equilateral min-angle %.4f != 60\n", ma); fails++; }
        if (fabs(q - 1.0) > 1e-3)   { fprintf(stderr, "selftest: equilateral q %.4f != 1\n", q); fails++; }
    }
    /* right isoceles -> min angle 45 */
    {
        float V[9] = { 0,0,0,  0,1,0,  0,0,1 };
        int32_t F[3] = {0,1,2};
        double ma, q; tri_metrics(V, F[0], F[1], F[2], &ma, &q);
        if (fabs(ma - 45.0) > 1e-3) { fprintf(stderr, "selftest: right-iso min-angle %.4f != 45\n", ma); fails++; }
        (void)q;
    }
    /* needle sliver -> tiny min angle, small q */
    {
        float V[9] = { 0,0,0,  0,1,0,  0,0.5f,0.001f };
        int32_t F[3] = {0,1,2};
        double ma, q; tri_metrics(V, F[0], F[1], F[2], &ma, &q);
        if (ma >= 1.0) { fprintf(stderr, "selftest: needle min-angle %.4f not < 1\n", ma); fails++; }
        if (q >= 0.05) { fprintf(stderr, "selftest: needle q %.4f not < 0.05\n", q); fails++; }
    }
    /* degenerate (zero-length edge) -> 0,0 */
    {
        float V[9] = { 0,0,0,  0,0,0,  0,0,1 };
        int32_t F[3] = {0,1,2};
        double ma, q; tri_metrics(V, F[0], F[1], F[2], &ma, &q);
        if (ma != 0.0 || q != 0.0) { fprintf(stderr, "selftest: degenerate not (0,0): %.4f %.4f\n", ma, q); fails++; }
    }
    /* two-triangle quad: shared edge interior, all four corners on boundary,
     * report() must run without touching out-of-range memory */
    {
        float V[12] = { 0,0,0,  0,1,0,  0,1,1,  0,0,1 };
        int32_t F[6] = { 0,1,2,  0,2,3 };
        report("selftest-quad", V, 4, F, 2, a);
    }
    Arena_dispose(&a);
    fprintf(stderr, "mesh_quality selftest: %s\n", fails==0 ? "OK" : "FAIL");
    return fails == 0 ? 0 : 1;
}

int main(int argc, char **argv) {
    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0) return selftest();
    if (argc < 2) {
        fprintf(stderr, "usage: %s <mesh.obj>\n       %s --selftest\n", argv[0], argv[0]);
        return 2;
    }
    Arena_T a = Arena_new();
    float *V; size_t nv; int32_t *F; size_t nf;
    if (ObjIO_read(a, argv[1], &V, &nv, &F, &nf) != 0) {
        fprintf(stderr, "error: cannot read %s\n", argv[1]);
        Arena_dispose(&a);
        return 1;
    }
    report(argv[1], V, nv, F, nf, a);
    Arena_dispose(&a);
    return 0;
}
