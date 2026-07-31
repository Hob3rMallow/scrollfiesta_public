/* mesh_resplit_test.c -- unit tests for MeshResplit_run (src/extract/mesh_resplit.c).
 *
 * Builds an input ComponentMesh holding two DISCONNECTED flat grid sheets in one
 * cloud and checks that MeshResplit separates them into two clean, edge-manifold
 * meshes; that a single-sheet input passes through (one component); and that a
 * NULL cloud passes everything through unchanged.
 *
 * Exit 0 = all pass, 1 = a failure.
 */
#include "extract/mesh_resplit.h"
#include "common/arena.h"
#include "common/mesh_types.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int g_fails = 0;
#define CHECK(cond, msg) do { \
    if (!(cond)) { printf("  FAIL: %s\n", (msg)); g_fails++; } \
    else         { printf("  ok:   %s\n", (msg)); } \
} while (0)

/* Fill an R x C flat grid sheet at z = zc into V/F at vert/face offsets. */
static void build_grid(float *V, int32_t *F, int voff, int foff,
                       int R, int C, float zc, size_t *pnv, size_t *pnf)
{
    int nv = 0;
    for (int r = 0; r < R; r++) for (int c = 0; c < C; c++) {
        int i = voff + nv;
        V[i*3+0] = zc; V[i*3+1] = (float)r; V[i*3+2] = (float)c;
        nv++;
    }
    int nf = 0;
    for (int r = 0; r < R-1; r++) for (int c = 0; c < C-1; c++) {
        int a = voff + r*C + c,     b = voff + r*C + (c+1);
        int d = voff + (r+1)*C + c, e = voff + (r+1)*C + (c+1);
        F[(foff+nf)*3+0]=a; F[(foff+nf)*3+1]=b; F[(foff+nf)*3+2]=d; nf++;
        F[(foff+nf)*3+0]=b; F[(foff+nf)*3+1]=e; F[(foff+nf)*3+2]=d; nf++;
    }
    *pnv = (size_t)nv; *pnf = (size_t)nf;
}

/* max edge multiplicity (>2 => non-manifold edge) for a small mesh. */
static int max_edge_mult(const ComponentMesh *m)
{
    size_t nv = m->nv;
    int *mult = (int *)calloc(nv*nv, sizeof(int));
    int mx = 0;
    for (size_t f = 0; f < m->nf; f++) {
        int v[3] = { m->faces[f*3+0], m->faces[f*3+1], m->faces[f*3+2] };
        int e[3][2] = {{v[0],v[1]},{v[1],v[2]},{v[0],v[2]}};
        for (int k = 0; k < 3; k++) {
            int a=e[k][0], b=e[k][1]; if(a>b){int t=a;a=b;b=t;}
            int c = ++mult[(size_t)a*nv+b];
            if (c > mx) mx = c;
        }
    }
    free(mult);
    return mx;
}

int main(void)
{
    printf("=== mesh_resplit_test ===\n");
    const int R = 16, C = 16;                 /* 256 verts/sheet (> RESPLIT_MIN_COMP_VERTS) */
    size_t per_nv = (size_t)R*C, per_nf = (size_t)2*(R-1)*(C-1);

    /* ---- input: two disconnected flat sheets (z=0 and z=20) in one mesh ---- */
    size_t nv = 2*per_nv, nf = 2*per_nf;
    float   *V = (float *)malloc(nv*3*sizeof(float));
    int32_t *F = (int32_t *)malloc(nf*3*sizeof(int32_t));
    size_t a_nv,a_nf,b_nv,b_nf;
    build_grid(V, F, 0,            0,             R, C, 0.0f,  &a_nv, &a_nf);
    build_grid(V, F, (int)per_nv,  (int)per_nf,   R, C, 20.0f, &b_nv, &b_nf);

    ComponentMesh in; memset(&in, 0, sizeof in);
    in.verts = V; in.faces = F; in.nv = nv; in.nf = nf; in.comp_id = 1; in.self = &in;

    /* cloud: originals == LOP images == the verts (identity for the test) */
    MeshResplitCloud cloud; memset(&cloud, 0, sizeof cloud);
    cloud.orig_pts = V; cloud.lop_pts = V; cloud.n = nv;
    cloud.cell_origin[0]=0; cloud.cell_origin[1]=0; cloud.cell_origin[2]=0;

    Arena_T arena = Arena_new();

    /* ---- test 1: two sheets -> two clean components ---- */
    printf("[two disconnected sheets]\n");
    {
        ComponentMesh *outm = NULL; size_t n_out = 0;
        int rc = MeshResplit_run(arena, &in, 1, &cloud, 1, NULL, &outm, &n_out, NULL);
        CHECK(rc == 0, "MeshResplit_run returns 0");
        printf("    n_out=%zu\n", n_out);
        CHECK(n_out == 2, "two disconnected sheets -> 2 components");
        int allmf = 1, allvalid = 1; float zmin=1e9f, zmax=-1e9f;
        for (size_t i = 0; i < n_out; i++) {
            if (!ComponentMesh_valid(&outm[i])) allvalid = 0;
            if (max_edge_mult(&outm[i]) > 2) allmf = 0;
            float zc = outm[i].centroid[0];
            if (zc < zmin) zmin = zc; if (zc > zmax) zmax = zc;
        }
        CHECK(allvalid, "both outputs are valid ComponentMeshes");
        CHECK(allmf, "both outputs are edge-manifold");
        CHECK(zmax - zmin > 15.0f, "the two outputs sit on the two separated sheets");
    }

    /* ---- test 2: single sheet -> pass through (one component) ---- */
    printf("[single sheet]\n");
    {
        ComponentMesh one; memset(&one, 0, sizeof one);
        one.verts = V; one.faces = F; one.nv = per_nv; one.nf = per_nf;
        one.comp_id = 1; one.self = &one;
        MeshResplitCloud c1; memset(&c1, 0, sizeof c1);
        c1.orig_pts = V; c1.lop_pts = V; c1.n = per_nv;
        ComponentMesh *outm = NULL; size_t n_out = 0;
        int rc = MeshResplit_run(arena, &one, 1, &c1, 1, NULL, &outm, &n_out, NULL);
        CHECK(rc == 0 && n_out == 1, "single sheet passes through as 1 component");
        CHECK(n_out == 1 && outm[0].verts == V, "pass-through reuses the input geometry");
    }

    /* ---- test 3: NULL cloud -> everything passes through ---- */
    printf("[NULL cloud]\n");
    {
        ComponentMesh *outm = NULL; size_t n_out = 0;
        int rc = MeshResplit_run(arena, &in, 1, NULL, 1, NULL, &outm, &n_out, NULL);
        CHECK(rc == 0 && n_out == 1, "NULL cloud -> input passes through unchanged");
    }

    /* ---- test 4: remesh_pieces (step 4 core) -> re-LOP given sub-meshes ---- */
    printf("[remesh_pieces: 2 pieces + parent cloud]\n");
    {
        /* the two sheets as standalone sub-meshes (piece 1's faces rebased) */
        ComponentMesh p[2]; memset(p, 0, sizeof p);
        int32_t *F1 = (int32_t *)malloc(per_nf*3*sizeof(int32_t));
        for (size_t t = 0; t < per_nf*3; t++) F1[t] = F[per_nf*3 + t] - (int32_t)per_nv;
        p[0].verts = V;            p[0].faces = F;  p[0].nv = per_nv; p[0].nf = per_nf; p[0].self = &p[0];
        p[1].verts = V + per_nv*3; p[1].faces = F1; p[1].nv = per_nv; p[1].nf = per_nf; p[1].self = &p[1];

        ComponentMesh *outm = NULL; size_t n_out = 0;
        int rc = MeshResplit_remesh_pieces(arena, &cloud, p, 2, NULL, NULL,
                                           &outm, &n_out, NULL);
        CHECK(rc == 0 && n_out == 2, "remesh_pieces: 2 pieces -> 2 meshes");
        int allmf = 1; float zmin = 1e9f, zmax = -1e9f;
        for (size_t i = 0; i < n_out; i++) {
            float vz0 = 1e9f, vz1 = -1e9f;
            for (size_t v = 0; v < outm[i].nv; v++) {
                float z = outm[i].verts[v*3+0];
                if (z < vz0) vz0 = z; if (z > vz1) vz1 = z;
            }
            float zc = 0.5f*(vz0+vz1);     /* geometric mid-z (centroid unset on fallback) */
            printf("    out[%zu]: nv=%zu nf=%zu vert_z[%.2f,%.2f] centroid_z=%.2f\n",
                   i, outm[i].nv, outm[i].nf, (double)vz0, (double)vz1,
                   (double)outm[i].centroid[0]);
            if (max_edge_mult(&outm[i]) > 2) allmf = 0;
            if (zc < zmin) zmin = zc; if (zc > zmax) zmax = zc;
        }
        CHECK(allmf, "remesh_pieces outputs are edge-manifold");
        CHECK(zmax - zmin > 15.0f, "remesh_pieces keeps the two sheets separated");
        free(F1);
    }

    /* ---- test 5: margin-gate drops seam-ambiguous points (no vacuuming) ----
     * Two CLOSE sheets (z=0 and z=2.5) as pieces; the cloud adds a full block of
     * ambiguous midpoint points (z=1.25) that are ~equidistant to both pieces.
     * The assignment margin must DROP them, not vacuum them into a piece (which
     * would inflate + fold its re-BPA'd sheet -- the step7_cc_bpa regression). */
    printf("[margin-gate: ambiguous midpoints dropped, not vacuumed]\n");
    {
        const float zB = 2.5f, zMid = 1.25f;
        float   *BV = (float *)malloc(per_nv*3*sizeof(float));
        int32_t *FB = (int32_t *)malloc(per_nf*3*sizeof(int32_t));
        for (size_t i=0;i<per_nv;i++){ BV[i*3+0]=zB; BV[i*3+1]=V[i*3+1]; BV[i*3+2]=V[i*3+2]; }
        for (size_t t=0;t<per_nf*3;t++) FB[t] = F[per_nf*3+t] - (int32_t)per_nv;

        ComponentMesh p[2]; memset(p,0,sizeof p);
        p[0].verts=V;  p[0].faces=F;  p[0].nv=per_nv; p[0].nf=per_nf; p[0].self=&p[0]; /* A @ z=0   */
        p[1].verts=BV; p[1].faces=FB; p[1].nv=per_nv; p[1].nf=per_nf; p[1].self=&p[1]; /* B @ z=2.5 */

        size_t cn = 3*per_nv;
        float *CV = (float *)malloc(cn*3*sizeof(float));
        for (size_t i=0;i<per_nv;i++){ CV[i*3+0]=0.0f; CV[i*3+1]=V[i*3+1]; CV[i*3+2]=V[i*3+2]; }
        for (size_t i=0;i<per_nv;i++){ size_t o=per_nv+i;   CV[o*3+0]=zB;   CV[o*3+1]=V[i*3+1]; CV[o*3+2]=V[i*3+2]; }
        for (size_t i=0;i<per_nv;i++){ size_t o=2*per_nv+i; CV[o*3+0]=zMid; CV[o*3+1]=V[i*3+1]; CV[o*3+2]=V[i*3+2]; }
        MeshResplitCloud cc; memset(&cc,0,sizeof cc);
        cc.orig_pts=CV; cc.lop_pts=CV; cc.n=cn;

        ComponentMesh *outm=NULL; size_t n_out=0;
        int rc = MeshResplit_remesh_pieces(arena, &cc, p, 2, NULL, NULL, &outm, &n_out, NULL);
        CHECK(rc==0 && n_out==2, "margin: 2 close pieces -> 2 meshes");
        size_t tot=0; int allmf=1;
        for (size_t i=0;i<n_out;i++){ tot+=outm[i].nv; if(max_edge_mult(&outm[i])>2) allmf=0; }
        printf("    total out nv=%zu (2 sheets=%zu; vacuumed midblock would push toward %zu)\n",
               tot, (size_t)(2*per_nv), (size_t)(3*per_nv));
        CHECK(allmf, "margin: outputs edge-manifold");
        CHECK(tot <= (size_t)((double)(2*per_nv)*1.25),
              "margin: midpoint block DROPPED not vacuumed (output not inflated)");
        free(BV); free(FB); free(CV);
    }

    /* ---- test 6: resurface_own re-LOP+re-BPAs from the piece's OWN verts ----
     * A sinusoidally-bumped grid sheet: re-surfacing must flatten the z
     * (through-thickness) variance, come out edge-manifold, preserve comp_id, and
     * NOT inflate the vertex count (it uses only the piece's own verts, so there
     * is no cloud to vacuum -- the step7_cc_bpa fold cannot occur). */
    printf("[resurface_own: re-LOP+re-BPA own verts -> clean, flatter, not inflated]\n");
    {
        float   *LV = (float *)malloc(per_nv*3*sizeof(float));
        int32_t *LF = (int32_t *)malloc(per_nf*3*sizeof(int32_t));
        memcpy(LF, F, per_nf*3*sizeof(int32_t));     /* sheet-A faces (0-based) */
        for (int r=0;r<R;r++) for (int c=0;c<C;c++) {
            size_t i=(size_t)(r*C+c);
            float bump = 1.0f*sinf((float)r*0.9f)*cosf((float)c*0.9f);
            LV[i*3+0]=bump; LV[i*3+1]=(float)r; LV[i*3+2]=(float)c;
        }
        ComponentMesh lm; memset(&lm,0,sizeof lm);
        lm.verts=LV; lm.faces=LF; lm.nv=per_nv; lm.nf=per_nf; lm.comp_id=7; lm.self=&lm;

        double mz0=0; for(size_t i=0;i<per_nv;i++) mz0+=LV[i*3]; mz0/=(double)per_nv;
        double zvar0=0; for(size_t i=0;i<per_nv;i++){ double d=LV[i*3]-mz0; zvar0+=d*d; } zvar0/=(double)per_nv;

        int rc = MeshResplit_resurface_own(arena, &lm, NULL, 10);
        CHECK(rc==0, "resurface_own returns 0");
        CHECK(ComponentMesh_valid(&lm) && lm.comp_id==7, "valid mesh, comp_id preserved");

        double mz1=0; for(size_t i=0;i<lm.nv;i++) mz1+=lm.verts[i*3]; mz1/=(double)lm.nv;
        double zvar1=0; for(size_t i=0;i<lm.nv;i++){ double d=lm.verts[i*3]-mz1; zvar1+=d*d; } zvar1/=(double)lm.nv;
        printf("    nv %zu -> %zu, z-var %.4f -> %.4f\n", (size_t)per_nv, lm.nv, zvar0, zvar1);

        int nan=0; for(size_t i=0;i<lm.nv*3;i++) if(lm.verts[i]!=lm.verts[i]) nan=1;
        CHECK(!nan, "no NaN after resurface");
        CHECK(max_edge_mult(&lm) <= 2, "re-surfaced mesh is edge-manifold");
        CHECK(zvar1 < zvar0*0.7, "through-thickness lumps flattened");
        CHECK(lm.nv <= (size_t)((double)per_nv*1.3), "not inflated (own verts only, no vacuuming)");
        free(LV); free(LF);
    }

    Arena_dispose(&arena);
    free(V); free(F);
    printf("=== %s (%d failure%s) ===\n",
           g_fails==0 ? "PASS" : "FAIL", g_fails, g_fails==1?"":"s");
    return g_fails==0 ? 0 : 1;
}
