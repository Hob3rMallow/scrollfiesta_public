/* rvd_clip_test — validates the restricted-Voronoi clipper via the partition
 * property: the RVD partitions the whole surface among the sites, so the sum of
 * restricted-cell areas must equal the mesh area. Also a hand-checkable two-site
 * split. Returns 0 on pass. */
#include "rvd_clip.h"
#include "../common/arena.h"

#include <stdio.h>
#include <math.h>

/* M x M grid of quads on the unit square in the z=0 plane; verts (z,y,x). */
static void build_grid(Arena_T a, int M,
                       double **PV, size_t *PNV, int32_t **PF, size_t *PNF) {
    int side = M + 1;
    size_t nv = (size_t)side * (size_t)side;
    size_t nf = (size_t)M * (size_t)M * 2;
    double  *V = (double *)ARENA_ALLOC(a, (long)(nv * 3 * sizeof(double)));
    int32_t *F = (int32_t *)ARENA_ALLOC(a, (long)(nf * 3 * sizeof(int32_t)));
    for (int i = 0; i < side; i++)
        for (int j = 0; j < side; j++) {
            size_t idx = (size_t)i * side + j;
            V[idx*3+0] = 0.0;                 /* z */
            V[idx*3+1] = (double)i / M;       /* y */
            V[idx*3+2] = (double)j / M;       /* x */
        }
    size_t fi = 0;
    for (int i = 0; i < M; i++)
        for (int j = 0; j < M; j++) {
            int32_t v00 = (int32_t)((size_t)i*side + j);
            int32_t v01 = (int32_t)((size_t)i*side + j + 1);
            int32_t v10 = (int32_t)((size_t)(i+1)*side + j);
            int32_t v11 = (int32_t)((size_t)(i+1)*side + j + 1);
            F[fi*3+0]=v00; F[fi*3+1]=v10; F[fi*3+2]=v11; fi++;
            F[fi*3+0]=v00; F[fi*3+1]=v11; F[fi*3+2]=v01; fi++;
        }
    *PV = V; *PNV = nv; *PF = F; *PNF = nf;
}

int main(void) {
    Arena_T a = Arena_new();
    double *V; size_t nv; int32_t *F; size_t nf;
    build_grid(a, 20, &V, &nv, &F, &nf);

    Rvd_T r = Rvd_new(a, V, nv, F, nf);
    double area = Rvd_mesh_area(r);
    printf("mesh area = %.9f (expect 1.0)\n", area);

    int ok = 1;

    /* Test 1: two sites at x=0.25 and x=0.75 -> bisector x=0.5 -> two equal halves. */
    {
        double s[6] = { 0,0.5,0.25,  0,0.5,0.75 };
        double ar[2], ce[6];
        if (Rvd_accumulate(r, s, 2, ar, ce, NULL, NULL, NULL, NULL, NULL) != 0) { printf("FAIL accumulate(2)\n"); ok = 0; }
        printf("2-site areas = %.9f %.9f (expect 0.5 0.5)\n", ar[0], ar[1]);
        printf("2-site x     = %.6f %.6f (expect 0.25 0.75)\n", ce[2], ce[5]);
        if (fabs(ar[0]-0.5) > 1e-6 || fabs(ar[1]-0.5) > 1e-6) { printf("FAIL 2-site area\n"); ok = 0; }
        if (fabs((ar[0]+ar[1]) - area) > 1e-9)                { printf("FAIL 2-site partition\n"); ok = 0; }
        if (fabs(ce[2]-0.25) > 1e-3 || fabs(ce[5]-0.75) > 1e-3){ printf("FAIL 2-site centroid\n"); ok = 0; }
    }

    /* Test 2: 9 interior sites -> partition of the whole square. */
    {
        double s[27]; int k = 0;
        for (int i = 1; i <= 3; i++)
            for (int j = 1; j <= 3; j++) { s[k*3+0]=0; s[k*3+1]=i/4.0; s[k*3+2]=j/4.0; k++; }
        double ar[9], ce[27];
        int32_t *tris = NULL; size_t nt = 0;
        if (Rvd_accumulate(r, s, 9, ar, ce, NULL, NULL, &tris, &nt, NULL) != 0) { printf("FAIL accumulate(9)\n"); ok = 0; }
        double sum = 0; for (int i = 0; i < 9; i++) sum += ar[i];
        printf("9-site area sum = %.9f (expect %.9f)\n", sum, area);
        if (fabs(sum - area) > 1e-8) { printf("FAIL 9-site partition\n"); ok = 0; }
        printf("9-site dual triangles = %zu\n", nt);
        int idx_ok = 1;
        for (size_t t = 0; t < nt*3; t++) if (tris[t] < 0 || tris[t] >= 9) idx_ok = 0;
        if (!idx_ok) { printf("FAIL dual index range\n"); ok = 0; }
        if (nt == 0)  { printf("FAIL dual produced no triangles\n"); ok = 0; }
    }

    Arena_dispose(&a);
    printf("%s\n", ok ? "rvd_clip_test: OK" : "rvd_clip_test: FAIL");
    return ok ? 0 : 1;
}
