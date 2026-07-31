/* patch_repair_test.c — unit tests for the tangent-plane Delaunay tear repair
 * (src/remesh/patch_repair.c). Whitebox: #includes the .c to reach statics.
 *
 * Covers the 2026-07-09 lightning-bolt fix:
 *   - a jagged interior tear in a flat sheet heals to a single-perimeter
 *     2-manifold (no interior loops, no bad-degree boundary verts);
 *   - a PINCHED tear (two holes sharing a vertex, boundary degree 4) heals;
 *   - CROSS-SHEET IMPOSSIBILITY: with two stacked sheets 0.8 vox apart,
 *     repairing a tear in one adds ZERO edges between the sheets and leaves
 *     them separate components (the user-mandated invariant);
 *   - a tear on a sharply folded (non-planar) patch is SKIPPED;
 *   - a bite out of the sheet PERIMETER is not treated as a tear;
 *   - determinism: two runs on the same input are byte-identical;
 *   - trivial inputs (nf==0) return a clean copy.
 *
 * Exit 0 = all pass, 1 = a failure.
 */
#include "remesh/patch_repair.c"   /* whitebox include */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int g_fails = 0;
#define CHECK(cond, msg) do { \
    if (!(cond)) { printf("  FAIL: %s\n", (msg)); g_fails++; } \
    else         { printf("  ok:   %s\n", (msg)); } \
} while (0)

/* ---------- mesh analysis (mirrors ball_pivot_test.c) -------------------- */
typedef struct {
    int max_edge_mult;
    int boundary_edges;
    int boundary_loops;    /* valid iff boundary_bad_deg == 0 */
    int boundary_bad_deg;
    int components;
    int cross_set_edges;   /* edges straddling the index `split` boundary */
} MeshStats;

static int tuf_find(int *p, int x){ while(p[x]!=x){p[x]=p[p[x]];x=p[x];} return x; }
static void tuf_union(int *p,int a,int b){ a=tuf_find(p,a); b=tuf_find(p,b); if(a!=b)p[a]=b; }

static MeshStats analyze(int nv, const int32_t *F, size_t nf, int split)
{
    MeshStats s; memset(&s, 0, sizeof s);
    int *mult = (int *)calloc((size_t)nv*(size_t)nv, sizeof(int));
    int *uf   = (int *)malloc((size_t)nv*sizeof(int));
    int *bdeg = (int *)calloc((size_t)nv, sizeof(int));
    int *bn1  = (int *)malloc((size_t)nv*sizeof(int));
    int *bn2  = (int *)malloc((size_t)nv*sizeof(int));
    for (int i=0;i<nv;i++){ uf[i]=i; bn1[i]=-1; bn2[i]=-1; }
    for (size_t f=0; f<nf; f++){
        int v[3] = { F[f*3+0], F[f*3+1], F[f*3+2] };
        int e[3][2] = {{v[0],v[1]},{v[1],v[2]},{v[0],v[2]}};
        for (int k=0;k<3;k++){
            int u=e[k][0], w=e[k][1]; if(u>w){int t=u;u=w;w=t;}
            mult[(size_t)u*nv+w]++;
            tuf_union(uf,u,w);
        }
    }
    for (int u=0;u<nv;u++) for (int w=u+1;w<nv;w++){
        int m = mult[(size_t)u*nv+w];
        if (m==0) continue;
        if (m > s.max_edge_mult) s.max_edge_mult = m;
        if ((u<split) != (w<split)) s.cross_set_edges++;
        if (m==1){
            s.boundary_edges++;
            if (bdeg[u]==0) bn1[u]=w; else if (bdeg[u]==1) bn2[u]=w; bdeg[u]++;
            if (bdeg[w]==0) bn1[w]=u; else if (bdeg[w]==1) bn2[w]=u; bdeg[w]++;
        }
    }
    {
        int *seen=(int*)calloc((size_t)nv,sizeof(int));
        for (size_t f=0; f<nf; f++) for(int k=0;k<3;k++){
            int r=tuf_find(uf,F[f*3+k]); if(!seen[r]){seen[r]=1;s.components++;}
        }
        free(seen);
    }
    for (int i=0;i<nv;i++) if (bdeg[i]!=0 && bdeg[i]!=2) s.boundary_bad_deg++;
    if (s.boundary_bad_deg==0){
        int *vis=(int*)calloc((size_t)nv,sizeof(int));
        for (int i=0;i<nv;i++) if (bdeg[i]==2 && !vis[i]){
            int cur=i, prev=-1, guard=0;
            while (cur!=-1 && !vis[cur] && guard++ < nv+2){
                vis[cur]=1;
                int nx = (bn1[cur]!=prev)?bn1[cur]:bn2[cur];
                prev=cur; cur=nx;
            }
            s.boundary_loops++;
        }
        free(vis);
    }
    free(mult); free(uf); free(bdeg); free(bn1); free(bn2);
    return s;
}

/* ---------- synthetic sheet builders ------------------------------------- */
/* R x C vertex grid at z=zc (coords z,y,x), spacing 1.0, split-diagonal
 * triangulation, consistent winding. Verts written at [off..off+R*C),
 * faces appended to F (each cell -> 2 tris). Returns vert count added. */
static int build_grid(float *V, int off, int R, int C, float zc,
                      int32_t *F, int *nf)
{
    for (int r=0;r<R;r++) for (int c=0;c<C;c++){
        int i = off + r*C + c;
        V[i*3+0]=zc; V[i*3+1]=(float)r; V[i*3+2]=(float)c;
    }
    for (int r=0;r<R-1;r++) for (int c=0;c<C-1;c++){
        int a=off+r*C+c, b=off+r*C+c+1, d=off+(r+1)*C+c, e=off+(r+1)*C+c+1;
        F[(*nf)*3+0]=a; F[(*nf)*3+1]=d; F[(*nf)*3+2]=b; (*nf)++;
        F[(*nf)*3+0]=b; F[(*nf)*3+1]=d; F[(*nf)*3+2]=e; (*nf)++;
    }
    return R*C;
}

/* delete faces whose all-3 verts are in cellset: helper via predicate */
typedef int (*face_pred)(const int32_t *fv, const float *V);
static size_t drop_faces(int32_t *F, size_t nf, const float *V, face_pred kill)
{
    size_t w = 0;
    for (size_t f=0; f<nf; f++){
        if (kill(&F[f*3], V)) continue;
        F[w*3+0]=F[f*3+0]; F[w*3+1]=F[f*3+1]; F[w*3+2]=F[f*3+2];
        w++;
    }
    return w;
}

/* jagged interior hole: kill faces whose centroid lies in a zigzag strip
 * around y~5..6, x in 4..8 (of sheet at z==0 only) */
static int kill_zigzag(const int32_t *fv, const float *V)
{
    double cz=0, cy=0, cx=0;
    for (int k=0;k<3;k++){ cz+=V[fv[k]*3+0]; cy+=V[fv[k]*3+1]; cx+=V[fv[k]*3+2]; }
    cz/=3; cy/=3; cx/=3;
    if (fabs(cz) > 0.1) return 0;                 /* sheet A only            */
    if (cx < 4.0 || cx > 8.0) return 0;
    double yc = 5.0 + ((int)floor(cx) % 2 ? 0.6 : -0.6);   /* zigzag        */
    return fabs(cy - yc) < 1.2;
}

/* pinched pair: two small holes sharing vertex (6,6) */
static int kill_pinch(const int32_t *fv, const float *V)
{
    double cy=0, cx=0, cz=0;
    for (int k=0;k<3;k++){ cz+=V[fv[k]*3+0]; cy+=V[fv[k]*3+1]; cx+=V[fv[k]*3+2]; }
    cz/=3; cy/=3; cx/=3;
    if (fabs(cz) > 0.1) return 0;
    int hole1 = (cy > 5.0 && cy < 6.0 && cx > 5.0 && cx < 6.0);
    int hole2 = (cy > 6.0 && cy < 7.0 && cx > 6.0 && cx < 7.0);
    return hole1 || hole2;
}

/* bite out of the perimeter: kill faces near a border corner */
static int kill_perimeter_bite(const int32_t *fv, const float *V)
{
    double cy=0, cx=0;
    for (int k=0;k<3;k++){ cy+=V[fv[k]*3+1]; cx+=V[fv[k]*3+2]; }
    cy/=3; cx/=3;
    return (cy < 1.5 && cx < 2.5);
}

/* ---------- tests --------------------------------------------------------- */

static void test_trivial(void)
{
    printf("[trivial inputs]\n");
    Arena_T a = Arena_new();
    int32_t *of = NULL; size_t onf = 777;
    PatchRepairStats st;
    int rc = PatchRepair_repair(a, NULL, 0, NULL, 0, &of, &onf, &st);
    CHECK(rc == 0 && onf == 0 && of != NULL, "nf==0 -> clean empty copy");
    Arena_dispose(&a);
}

static void test_tear_heal(void)
{
    printf("[jagged interior tear heals]\n");
    enum { R = 12, C = 13 };
    int nv = R*C;
    float *V = (float *)malloc((size_t)nv*3*sizeof(float));
    int32_t *F = (int32_t *)malloc((size_t)(R-1)*(C-1)*2*3*sizeof(int32_t));
    int nf = 0;
    build_grid(V, 0, R, C, 0.0f, F, &nf);
    size_t nf2 = drop_faces(F, (size_t)nf, V, kill_zigzag);
    MeshStats before = analyze(nv, F, nf2, nv);
    /* perimeter of a full R x C grid is 2(R-1)+2(C-1) edges; more boundary
     * than that == an interior tear exists (loops uncountable when the tear
     * pinches, so count edges instead) */
    CHECK(before.boundary_edges > 2*(R-1) + 2*(C-1), "input really has an interior tear");

    Arena_T a = Arena_new();
    int32_t *of = NULL; size_t onf = 0;
    PatchRepairStats st;
    int rc = PatchRepair_repair(a, V, (size_t)nv, F, nf2, &of, &onf, &st);
    printf("    tears=%zu patches=%zu repaired=%zu (+%zu -%zu) skips g%zu r%zu p%zu t%zu\n",
           st.n_tear_clusters, st.n_patches, st.n_repaired,
           st.faces_added, st.faces_deleted,
           st.n_skip_gate, st.n_skip_ring, st.n_skip_project, st.n_skip_tri);
    CHECK(rc == 0, "repair returns 0");
    CHECK(st.n_repaired >= 1, "at least one patch repaired");
    MeshStats after = analyze(nv, of, onf, nv);
    printf("    after: bE=%d loops=%d baddeg=%d maxmult=%d comps=%d\n",
           after.boundary_edges, after.boundary_loops,
           after.boundary_bad_deg, after.max_edge_mult, after.components);
    CHECK(after.max_edge_mult <= 2, "edge-manifold after repair");
    CHECK(after.boundary_bad_deg == 0, "no bad-degree boundary verts");
    CHECK(after.boundary_loops == 1, "single boundary loop (tear gone)");
    CHECK(after.components == 1, "still one component");
    Arena_dispose(&a);
    free(V); free(F);
}

static void test_pinched_tear(void)
{
    printf("[pinched (figure-8) tear heals]\n");
    enum { R = 13, C = 13 };
    int nv = R*C;
    float *V = (float *)malloc((size_t)nv*3*sizeof(float));
    int32_t *F = (int32_t *)malloc((size_t)(R-1)*(C-1)*2*3*sizeof(int32_t));
    int nf = 0;
    build_grid(V, 0, R, C, 0.0f, F, &nf);
    size_t nf2 = drop_faces(F, (size_t)nf, V, kill_pinch);
    MeshStats before = analyze(nv, F, nf2, nv);
    CHECK(before.boundary_edges > (R-1+C-1)*2, "input has interior boundary");

    Arena_T a = Arena_new();
    int32_t *of = NULL; size_t onf = 0;
    PatchRepairStats st;
    int rc = PatchRepair_repair(a, V, (size_t)nv, F, nf2, &of, &onf, &st);
    printf("    tears=%zu repaired=%zu skips g%zu r%zu p%zu t%zu\n",
           st.n_tear_clusters, st.n_repaired,
           st.n_skip_gate, st.n_skip_ring, st.n_skip_project, st.n_skip_tri);
    CHECK(rc == 0 && st.n_repaired >= 1, "pinched tear repaired");
    MeshStats after = analyze(nv, of, onf, nv);
    CHECK(after.max_edge_mult <= 2 && after.boundary_bad_deg == 0 &&
          after.boundary_loops == 1, "healed to single-perimeter manifold");
    Arena_dispose(&a);
    free(V); free(F);
}

static void test_cross_sheet_impossible(void)
{
    printf("[cross-sheet impossibility: stacked sheets 0.8 vox apart]\n");
    enum { R = 12, C = 13 };
    int half = R*C, nv = 2*half;
    float *V = (float *)malloc((size_t)nv*3*sizeof(float));
    int32_t *F = (int32_t *)malloc((size_t)(R-1)*(C-1)*2*2*3*sizeof(int32_t));
    int nf = 0;
    build_grid(V, 0,    R, C, 0.0f, F, &nf);
    build_grid(V, half, R, C, 0.8f, F, &nf);
    size_t nf2 = drop_faces(F, (size_t)nf, V, kill_zigzag);   /* tear in A only */

    Arena_T a = Arena_new();
    int32_t *of = NULL; size_t onf = 0;
    PatchRepairStats st;
    int rc = PatchRepair_repair(a, V, (size_t)nv, F, nf2, &of, &onf, &st);
    CHECK(rc == 0 && st.n_repaired >= 1, "tear in sheet A repaired");
    MeshStats after = analyze(nv, of, onf, half);
    printf("    after: cross=%d comps=%d loops=%d\n",
           after.cross_set_edges, after.components, after.boundary_loops);
    CHECK(after.cross_set_edges == 0, "ZERO edges between the sheets");
    CHECK(after.components == 2, "sheets remain separate components");
    CHECK(after.boundary_loops == 2, "two perimeters, no interior loops");
    /* sheet B faces byte-identical: every output face fully in B equals an
     * input face fully in B (kept order) */
    size_t inB_in = 0, inB_out = 0;
    for (size_t f=0; f<nf2;  f++) if (F[f*3]>=half && F[f*3+1]>=half && F[f*3+2]>=half) inB_in++;
    for (size_t f=0; f<onf; f++) if (of[f*3]>=half && of[f*3+1]>=half && of[f*3+2]>=half) inB_out++;
    CHECK(inB_in == inB_out, "sheet B untouched (face count identical)");
    Arena_dispose(&a);
    free(V); free(F);
}

static void test_planarity_skip(void)
{
    printf("[non-planar patch skipped]\n");
    enum { R = 12, C = 13 };
    int nv = R*C;
    float *V = (float *)malloc((size_t)nv*3*sizeof(float));
    int32_t *F = (int32_t *)malloc((size_t)(R-1)*(C-1)*2*3*sizeof(int32_t));
    int nf = 0;
    build_grid(V, 0, R, C, 0.0f, F, &nf);
    /* tear FIRST on the flat sheet, then fold sharply along x=6:
     * z = 2.5 * |x-6| (patch PCA thickness >> gate). The tear straddles the
     * fold, so its patch is decisively non-planar. */
    size_t nf2 = drop_faces(F, (size_t)nf, V, kill_zigzag);
    for (int i=0;i<nv;i++) V[i*3+0] = 2.5f * (float)fabs((double)V[i*3+2] - 6.0);
    MeshStats before = analyze(nv, F, nf2, nv);

    Arena_T a = Arena_new();
    int32_t *of = NULL; size_t onf = 0;
    PatchRepairStats st;
    int rc = PatchRepair_repair(a, V, (size_t)nv, F, nf2, &of, &onf, &st);
    printf("    tears=%zu repaired=%zu skip_gate=%zu\n",
           st.n_tear_clusters, st.n_repaired, st.n_skip_gate);
    CHECK(rc == 0, "repair returns 0");
    CHECK(st.n_tear_clusters >= 1, "fold tear is detected (not vacuous)");
    CHECK(st.n_repaired == 0 && onf == nf2, "fold patch skipped, mesh unchanged");
    MeshStats after = analyze(nv, of, onf, nv);
    CHECK(after.boundary_edges == before.boundary_edges, "boundary identical");
    Arena_dispose(&a);
    free(V); free(F);
}

static void test_perimeter_bite(void)
{
    printf("[perimeter bite is not a tear]\n");
    enum { R = 12, C = 13 };
    int nv = R*C;
    float *V = (float *)malloc((size_t)nv*3*sizeof(float));
    int32_t *F = (int32_t *)malloc((size_t)(R-1)*(C-1)*2*3*sizeof(int32_t));
    int nf = 0;
    build_grid(V, 0, R, C, 0.0f, F, &nf);
    size_t nf2 = drop_faces(F, (size_t)nf, V, kill_perimeter_bite);

    Arena_T a = Arena_new();
    int32_t *of = NULL; size_t onf = 0;
    PatchRepairStats st;
    int rc = PatchRepair_repair(a, V, (size_t)nv, F, nf2, &of, &onf, &st);
    printf("    tears=%zu repaired=%zu\n", st.n_tear_clusters, st.n_repaired);
    CHECK(rc == 0, "repair returns 0");
    CHECK(st.n_tear_clusters == 0, "bite merges with perimeter cluster -> no tear");
    CHECK(onf == nf2, "mesh unchanged");
    Arena_dispose(&a);
    free(V); free(F);
}

static void test_determinism(void)
{
    printf("[determinism]\n");
    enum { R = 12, C = 13 };
    int nv = R*C;
    float *V = (float *)malloc((size_t)nv*3*sizeof(float));
    int32_t *F = (int32_t *)malloc((size_t)(R-1)*(C-1)*2*3*sizeof(int32_t));
    int nf = 0;
    build_grid(V, 0, R, C, 0.0f, F, &nf);
    size_t nf2 = drop_faces(F, (size_t)nf, V, kill_zigzag);

    Arena_T a1 = Arena_new(), a2 = Arena_new();
    int32_t *o1 = NULL, *o2 = NULL; size_t n1 = 0, n2 = 0;
    PatchRepair_repair(a1, V, (size_t)nv, F, nf2, &o1, &n1, NULL);
    PatchRepair_repair(a2, V, (size_t)nv, F, nf2, &o2, &n2, NULL);
    CHECK(n1 == n2 && n1 > 0, "same face count");
    CHECK(memcmp(o1, o2, n1*3*sizeof(int32_t)) == 0, "byte-identical output");
    Arena_dispose(&a1); Arena_dispose(&a2);
    free(V); free(F);
}

int main(void)
{
    setvbuf(stdout, NULL, _IONBF, 0);
    printf("=== patch_repair_test ===\n");
    test_trivial();
    test_tear_heal();
    test_pinched_tear();
    test_cross_sheet_impossible();
    test_planarity_skip();
    test_perimeter_bite();
    test_determinism();
    printf("=== %s (%d failure%s) ===\n",
           g_fails==0 ? "PASS" : "FAIL", g_fails, g_fails==1?"":"s");
    return g_fails==0 ? 0 : 1;
}
