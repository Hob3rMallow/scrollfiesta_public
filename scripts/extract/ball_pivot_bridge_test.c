/*
 * ball_pivot_bridge_test.c — unit test for BallPivot_bridge (init-front
 * seam weld). Builds two independently triangulated half-sheets separated
 * by a gap, drops dense band "stepping-stone" point columns in the gap,
 * and asks BallPivot_bridge to weld the two open boundaries across it.
 *
 * Success criteria match exactly what grid_weld audits on a real weld:
 *   - manifold: no edge in >2 faces
 *   - seam closed: every facing-boundary edge of both halves becomes paired
 *
 * Two scenarios:
 *   narrow gap (2 vox): rho=1.5 can span directly, band points optional.
 *   wide gap   (4 vox): direct span impossible, so a closed seam PROVES the
 *                       band stepping-stones were used — the property the
 *                       grid_weld seam weld relies on.
 *
 * The sheet is gently curved along x with correct per-vertex normals,
 * mirroring real LOP data (a flat sparse patch makes BPA pivots degenerate).
 *
 * Compile (from repo root, in a VS dev shell):
 *   cl /nologo /W4 /Zi /Fe:output\tools\ball_pivot_bridge_test.exe ^
 *      /Fo:output\tools\ ^
 *      scripts\step0-mesh-extract\ball_pivot_bridge_test.c ^
 *      src\remesh\ball_pivot.c src\common\arena.c src\common\except.c
 */
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../src/common/arena.h"
#include "../../src/remesh/ball_pivot.h"

#define NROW 5
#define MAXCOL 8
#define VIDX(col, row) ((col) * NROW + (row))

static int edge_run(const int32_t *faces, size_t nf, int32_t u, int32_t v)
{
    int run = 0;
    for (size_t f = 0; f < nf; f++) {
        int32_t a = faces[f*3+0], b = faces[f*3+1], c = faces[f*3+2];
        if ((a==u&&b==v)||(b==u&&c==v)||(c==u&&a==v) ||
            (a==v&&b==u)||(b==v&&c==u)||(c==v&&a==u)) run++;
    }
    return run;
}

static int max_edge_run(const int32_t *faces, size_t nf)
{
    int worst = 0;
    for (size_t f = 0; f < nf; f++)
        for (int e = 0; e < 3; e++) {
            int r = edge_run(faces, nf, faces[f*3+e], faces[f*3+((e+1)%3)]);
            if (r > worst) worst = r;
        }
    return worst;
}

static int vert_face_count(const int32_t *faces, size_t nf, int32_t v)
{
    int c = 0;
    for (size_t f = 0; f < nf; f++)
        if (faces[f*3+0]==v || faces[f*3+1]==v || faces[f*3+2]==v) c++;
    return c;
}

static size_t collect_boundary_edges(const int32_t *faces, size_t nf,
                                     BpaInitEdge *out, size_t out_cap)
{
    size_t n = 0;
    for (size_t f = 0; f < nf; f++) {
        int32_t v[3] = { faces[f*3+0], faces[f*3+1], faces[f*3+2] };
        for (int e = 0; e < 3; e++) {
            int32_t u = v[e], w = v[(e+1)%3], opp = v[(e+2)%3];
            if (edge_run(faces, nf, u, w) == 1 && n < out_cap) {
                out[n].va = u; out[n].vb = w; out[n].v_opp = opp; n++;
            }
        }
    }
    return n;
}

/* Regular-grid triangulation of columns [c0, c0+1] over all rows. */
static size_t tri_strip(int c0, int32_t *out)
{
    size_t n = 0;
    for (int r = 0; r + 1 < NROW; r++) {
        int a = VIDX(c0, r),   b = VIDX(c0, r+1);
        int d = VIDX(c0+1, r), e = VIDX(c0+1, r+1);
        out[n*3+0]=a; out[n*3+1]=d; out[n*3+2]=e; n++;
        out[n*3+0]=a; out[n*3+1]=e; out[n*3+2]=b; n++;
    }
    return n;
}

/* One bridging scenario. left mesh = cols {0,1}; right mesh = cols
 * {rc, rc+1}; band = the columns strictly between 1 and rc.
 * `band_required` => assert no band column is wholly orphaned. */
static int run_scenario(const char *name, int rc, int band_required)
{
    Arena_T arena = Arena_new();
    int ncol = rc + 2;
    int nv = ncol * NROW;
    printf("== %s: left{0,1} right{%d,%d} band{2..%d} (%d cols) ==\n",
           name, rc, rc+1, rc-1, ncol);

    float *verts   = (float *)malloc((size_t)nv*3*sizeof(float));
    float *normals = (float *)malloc((size_t)nv*3*sizeof(float));
    for (int col = 0; col < ncol; col++) {
        double fx  = 0.4 * sin(1.1 * col);
        double dfx = 0.4 * 1.1 * cos(1.1 * col);
        double nz = 1.0, ny = 0.0, nx = -dfx;
        double nl = sqrt(nz*nz + ny*ny + nx*nx);
        for (int row = 0; row < NROW; row++) {
            int i = VIDX(col, row);
            verts[i*3+0] = (float)fx;
            verts[i*3+1] = (float)row;
            verts[i*3+2] = (float)col;
            normals[i*3+0] = (float)(nz/nl);
            normals[i*3+1] = (float)(ny/nl);
            normals[i*3+2] = (float)(nx/nl);
        }
    }

    int32_t left_faces[2*(NROW-1)*3];
    int32_t right_faces[2*(NROW-1)*3];
    size_t left_nf  = tri_strip(0,  left_faces);
    size_t right_nf = tri_strip(rc, right_faces);

    BpaInitEdge init[256];
    size_t n_init = 0;
    n_init += collect_boundary_edges(left_faces,  left_nf,  init+n_init, 256-n_init);
    n_init += collect_boundary_edges(right_faces, right_nf, init+n_init, 256-n_init);

    int32_t *bridge_faces = NULL;
    size_t bridge_nf = 0;
    int rcc = BallPivot_bridge(arena, verts, normals, (size_t)nv, 1.5f, 1.5f,
                               init, n_init, NULL /* gate off */,
                               &bridge_faces, &bridge_nf);
    printf("  bridge rc=%d, %zu faces\n", rcc, bridge_nf);

    size_t total_nf = left_nf + right_nf + bridge_nf;
    int32_t *all = (int32_t *)malloc(total_nf*3*sizeof(int32_t));
    size_t off = 0;
    memcpy(all+off, left_faces,   left_nf*3*sizeof(int32_t));  off += left_nf*3;
    memcpy(all+off, right_faces,  right_nf*3*sizeof(int32_t)); off += right_nf*3;
    memcpy(all+off, bridge_faces, bridge_nf*3*sizeof(int32_t));

    int failures = 0;
    if (rcc != 0)       { printf("  FAIL: rc=%d\n", rcc); failures++; }
    if (bridge_nf == 0) { printf("  FAIL: no bridge faces\n"); failures++; }

    int worst = max_edge_run(all, total_nf);
    printf("  max edge multiplicity: %d\n", worst);
    if (worst > 2) { printf("  FAIL: non-manifold (run=%d)\n", worst); failures++; }

    int unclosed = 0;
    for (int r = 0; r + 1 < NROW; r++) {
        if (edge_run(all, total_nf, VIDX(1,r),  VIDX(1,r+1))  != 2) unclosed++;
        if (edge_run(all, total_nf, VIDX(rc,r), VIDX(rc,r+1)) != 2) unclosed++;
    }
    printf("  unclosed facing edges: %d/%d\n", unclosed, 2*(NROW-1));
    if (unclosed != 0) { printf("  FAIL: seam not closed\n"); failures++; }

    if (band_required) {
        /* A wholly orphaned band column would mean the ball tunneled past
         * it — only possible if the gap were spannable, which it isn't. */
        int dead_cols = 0;
        for (int col = 2; col < rc; col++) {
            int used = 0;
            for (int r = 0; r < NROW; r++)
                if (vert_face_count(all, total_nf, VIDX(col,r)) > 0) used = 1;
            if (!used) dead_cols++;
        }
        printf("  fully-orphaned band columns: %d\n", dead_cols);
        if (dead_cols != 0) { printf("  FAIL: band column unused\n"); failures++; }
    }

    free(all); free(verts); free(normals);
    Arena_dispose(&arena);
    printf("  -> %s\n", failures==0 ? "PASS" : "FAIL");
    return failures;
}

int main(void)
{
    int failures = 0;
    failures += run_scenario("narrow gap (2 vox, band optional)", 3, 0);
    failures += run_scenario("wide gap (4 vox, band required)",   5, 1);

    if (failures == 0) { printf("\nALL BRIDGE TESTS PASSED\n"); return 0; }
    printf("\n%d BRIDGE TEST(S) FAILED\n", failures);
    return 1;
}
