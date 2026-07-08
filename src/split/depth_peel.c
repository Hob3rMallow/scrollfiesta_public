/*
 * depth_peel.c -- Stacked-wrap separator. See depth_peel.h. Pure splitter:
 * project verts onto the PCA normal to get a depth w, mark the seam where a
 * triangle edge JUMPS in depth by more than min_gap (the inter-wrap bridge edges),
 * remove the seam band (+ geometric dilation), return the connected survivors.
 * Mirrors dev_cut.c's cut-zone -> connected-components construction exactly; the
 * only difference is the seam test (per-edge depth-jump vs per-vertex Crane energy).
 */
#define _USE_MATH_DEFINES
#include "depth_peel.h"

#include "../common/kdtree.h"
#include "../common/union_find.h"
#include "../common/pca.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- sub-mesh construction (mirrors dev_cut.c) ---------------------------- */

static void finalize_submesh(ComponentMesh *m, float *v, size_t nv,
                             int32_t *f, size_t nf)
{
    memset(m, 0, sizeof(*m));
    m->verts = v; m->faces = f; m->nv = nv; m->nf = nf;
    PCA_normal(v, nv, m->pca_normal, m->centroid);
    m->pin_mask = NULL; m->vert_normals = NULL; m->nv_pre_fill = 0;
    m->self = m;
}

/* Pass `mesh` through unchanged as a single piece (shallow struct copy -- the
 * verts/faces buffers live in the same arena and outlive the call). */
static int emit_single(Arena_T arena, const ComponentMesh *mesh,
                       ComponentMesh **out, size_t *n)
{
    ComponentMesh *o = (ComponentMesh *)ARENA_ALLOC(arena, (long)sizeof(ComponentMesh));
    *o = *mesh; o->self = o;
    *out = o; *n = 1;
    return 0;
}

/* Build the sub-mesh of faces all three of whose vertices have vert_comp == c
 * (original positions, compacted). Returns 0 on success, -1 if degenerate. */
static int build_component(Arena_T arena, const ComponentMesh *mesh,
                           const int32_t *vert_comp, int32_t c, ComponentMesh *out)
{
    size_t N = mesh->nv, nf = mesh->nf;
    size_t cnf = 0;
    for (size_t f = 0; f < nf; f++) {
        int32_t a = mesh->faces[f*3+0], b = mesh->faces[f*3+1], d = mesh->faces[f*3+2];
        if (vert_comp[a] == c && vert_comp[b] == c && vert_comp[d] == c) cnf++;
    }
    if (cnf == 0) return -1;

    int32_t *o2n = (int32_t *)ARENA_ALLOC(arena, (long)(N*sizeof(int32_t)));
    for (size_t i = 0; i < N; i++) o2n[i] = -1;
    for (size_t f = 0; f < nf; f++) {
        int32_t a = mesh->faces[f*3+0], b = mesh->faces[f*3+1], d = mesh->faces[f*3+2];
        if (vert_comp[a] == c && vert_comp[b] == c && vert_comp[d] == c) {
            o2n[a] = 0; o2n[b] = 0; o2n[d] = 0;
        }
    }
    int32_t cnv = 0;
    for (size_t i = 0; i < N; i++) if (o2n[i] != -1) o2n[i] = cnv++;
    if (cnv < 3) return -1;

    float *vv = (float *)ARENA_ALLOC(arena, (long)((size_t)cnv*3*sizeof(float)));
    for (size_t i = 0; i < N; i++) if (o2n[i] >= 0) {
        size_t ni = (size_t)o2n[i];
        vv[ni*3+0] = mesh->verts[i*3+0];
        vv[ni*3+1] = mesh->verts[i*3+1];
        vv[ni*3+2] = mesh->verts[i*3+2];
    }
    int32_t *ff = (int32_t *)ARENA_ALLOC(arena, (long)(cnf*3*sizeof(int32_t)));
    size_t k = 0;
    for (size_t f = 0; f < nf; f++) {
        int32_t a = mesh->faces[f*3+0], b = mesh->faces[f*3+1], d = mesh->faces[f*3+2];
        if (vert_comp[a] == c && vert_comp[b] == c && vert_comp[d] == c) {
            ff[k*3+0] = o2n[a]; ff[k*3+1] = o2n[b]; ff[k*3+2] = o2n[d]; k++;
        }
    }
    finalize_submesh(out, vv, (size_t)cnv, ff, cnf);
    return 0;
}

/* ---- DepthPeel_process ---------------------------------------------------- */

int DepthPeel_process(Arena_T arena, const ComponentMesh *mesh,
                      double min_gap, double gap_depth, size_t min_comp_verts,
                      ComponentMesh **out_meshes, size_t *out_count)
{
    assert(arena && out_meshes && out_count);
    *out_meshes = NULL; *out_count = 0;
    if (!mesh || !mesh->verts || !mesh->faces) return -1;
    if (mesh->nv < 3 || mesh->nf < 1)
        return emit_single(arena, mesh, out_meshes, out_count);

    size_t nv = mesh->nv, nf = mesh->nf;
    Arena_Mark scratch = Arena_save(arena);

    /* 1. depth along the component thickness axis (PCA normal). */
    float *w = (float *)ARENA_ALLOC(arena, (long)(nv*sizeof(float)));
    float wmin = 0.0f, wmax = 0.0f;
    PCA_project(mesh->verts, nv, mesh->pca_normal, w, &wmin, &wmax);

    /* 2. seam vertices: an endpoint of any triangle edge whose depth jumps by
     * more than min_gap. Only an inter-wrap bridge edge (a long, near-normal
     * step across the empty gap) makes such a jump; intra-sheet edges (~0.6 vox)
     * never do, so a single sheet -- however folded -- marks nothing here. */
    unsigned char *seam = (unsigned char *)ARENA_CALLOC(arena, (long)nv, 1L);
    size_t n_seam = 0;
    float fmin_gap = (float)min_gap;
    for (size_t f = 0; f < nf; f++) {
        int32_t idx[3] = { mesh->faces[f*3+0], mesh->faces[f*3+1], mesh->faces[f*3+2] };
        for (int k = 0; k < 3; k++) {
            int32_t a = idx[k], b = idx[(k+1)%3];
            if (fabsf(w[a] - w[b]) > fmin_gap) {
                if (!seam[a]) { seam[a] = 1; n_seam++; }
                if (!seam[b]) { seam[b] = 1; n_seam++; }
            }
        }
    }
    if (n_seam == 0) {                              /* single sheet -> no cut */
        Arena_restore(arena, scratch);
        return emit_single(arena, mesh, out_meshes, out_count);
    }

    /* 3. exclusion band: seam verts dilated by gap_depth (KD-tree, like
     * bridge_cut's CUT_GAP_DEPTH). gap_depth <= 0 -> seam verts only. */
    unsigned char *excluded = (unsigned char *)ARENA_CALLOC(arena, (long)nv, 1L);
    if (gap_depth > 0.0) {
        float *spts = (float *)ARENA_ALLOC(arena, (long)(n_seam*3*sizeof(float)));
        size_t si = 0;
        for (size_t v = 0; v < nv; v++) if (seam[v]) {
            spts[si*3+0] = mesh->verts[v*3+0];
            spts[si*3+1] = mesh->verts[v*3+1];
            spts[si*3+2] = mesh->verts[v*3+2];
            si++;
        }
        KDTree_T tree = KDTree_new(arena, spts, n_seam);
        float gap_sq = (float)(gap_depth*gap_depth);
        for (size_t v = 0; v < nv; v++) {
            float d2 = 0.0f;
            (void)KDTree_nearest(tree, &mesh->verts[v*3], &d2);
            if (d2 < gap_sq) excluded[v] = 1;
        }
    } else {
        for (size_t v = 0; v < nv; v++) excluded[v] = seam[v];
    }

    /* 4. connected components of the SURVIVING faces (none of whose verts are
     * excluded). The removed band severs the mesh where the bridge spanned it. */
    UnionFind uf = UF_new(arena, (int32_t)nv);
    unsigned char *refd = (unsigned char *)ARENA_CALLOC(arena, (long)nv, 1L);
    for (size_t f = 0; f < nf; f++) {
        int32_t a = mesh->faces[f*3+0], b = mesh->faces[f*3+1], c = mesh->faces[f*3+2];
        if (excluded[a] || excluded[b] || excluded[c]) continue;
        uf_union(&uf, a, b); uf_union(&uf, b, c);
        refd[a] = 1; refd[b] = 1; refd[c] = 1;
    }

    /* 5. component sizes over referenced verts; promote big ones to real pieces. */
    int32_t *vcount = (int32_t *)ARENA_CALLOC(arena, (long)nv, (long)sizeof(int32_t));
    for (size_t v = 0; v < nv; v++) if (refd[v]) vcount[uf_find(&uf, (int32_t)v)]++;
    int32_t *comp_of_root = (int32_t *)ARENA_ALLOC(arena, (long)(nv*sizeof(int32_t)));
    for (size_t v = 0; v < nv; v++) comp_of_root[v] = -1;
    int32_t n_real = 0;
    for (size_t v = 0; v < nv; v++) {
        int32_t r = uf_find(&uf, (int32_t)v);
        if ((size_t)r == v && vcount[v] >= (int32_t)min_comp_verts)
            comp_of_root[v] = n_real++;
    }

    /* 6. per-vertex final component label (survives the arena restore via malloc). */
    int32_t *vert_comp = (int32_t *)malloc(nv*sizeof(int32_t));
    if (!vert_comp) { Arena_restore(arena, scratch); return -1; }
    for (size_t v = 0; v < nv; v++) {
        if (!refd[v]) { vert_comp[v] = -1; continue; }
        vert_comp[v] = comp_of_root[uf_find(&uf, (int32_t)v)];
    }

    Arena_restore(arena, scratch);

    if (n_real < 2) {                               /* nothing separated */
        free(vert_comp);
        return emit_single(arena, mesh, out_meshes, out_count);
    }

    /* 7. materialize the pieces. */
    ComponentMesh *o = (ComponentMesh *)ARENA_ALLOC(arena,
                          (long)((size_t)n_real*sizeof(ComponentMesh)));
    size_t cnt = 0;
    for (int32_t c = 0; c < n_real; c++) {
        if (build_component(arena, mesh, vert_comp, c, &o[cnt]) == 0) {
            o[cnt].comp_id = mesh->comp_id;
            cnt++;
        }
    }
    free(vert_comp);

    if (cnt < 2)                                    /* a piece failed to build */
        return emit_single(arena, mesh, out_meshes, out_count);

    *out_meshes = o; *out_count = cnt;
    return 0;
}

/* ============================================================================
 * Self-test: a flat sheet and a gently-curved sheet self-gate (one piece, no
 * large depth-jump), two parallel sheets welded by a wall bridge split into two,
 * and trivial input is a clean no-op.
 * ==========================================================================*/

static void make_grid(Arena_T arena, int ncol, int nrow, const float *z,
                      float **ov, size_t *onv, int32_t **of, size_t *onf)
{
    size_t nv = (size_t)ncol*(size_t)nrow;
    size_t nf = (size_t)(ncol-1)*(size_t)(nrow-1)*2;
    float   *v = (float *)ARENA_ALLOC(arena, (long)(nv*3*sizeof(float)));
    int32_t *f = (int32_t *)ARENA_ALLOC(arena, (long)(nf*3*sizeof(int32_t)));
    for (int r = 0; r < nrow; r++) for (int c = 0; c < ncol; c++) {
        size_t idx = (size_t)r*(size_t)ncol + (size_t)c;
        v[idx*3+0] = (float)c;
        v[idx*3+1] = (float)r;
        v[idx*3+2] = z[idx];
    }
    size_t fi = 0;
    for (int r = 0; r < nrow-1; r++) for (int c = 0; c < ncol-1; c++) {
        int32_t a = (int32_t)((size_t)r*(size_t)ncol + (size_t)c);
        int32_t b = (int32_t)((size_t)r*(size_t)ncol + (size_t)c + 1);
        int32_t cc = (int32_t)((size_t)(r+1)*(size_t)ncol + (size_t)c);
        int32_t d = (int32_t)((size_t)(r+1)*(size_t)ncol + (size_t)c + 1);
        f[fi*3+0]=a; f[fi*3+1]=b; f[fi*3+2]=cc; fi++;
        f[fi*3+0]=b; f[fi*3+1]=d; f[fi*3+2]=cc; fi++;
    }
    *ov = v; *onv = nv; *of = f; *onf = fi;
}

static void wrap_mesh(ComponentMesh *cm, float *v, size_t nv, int32_t *f, size_t nf)
{
    memset(cm, 0, sizeof(*cm));
    cm->verts = v; cm->faces = f; cm->nv = nv; cm->nf = nf;
    PCA_normal(v, nv, cm->pca_normal, cm->centroid);
    cm->self = cm;
}

int DepthPeel_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    fprintf(stderr, "=== DepthPeel selftest ===\n");

    const double min_gap = 1.5;   /* jump > 1.5 vox = bridge (intra-sheet ~1.0) */
    const double gap     = 0.5;   /* grid spacing 1.0 -> removes only seam verts */
    const size_t minv    = 10;

    /* (1) Flat sheet: no depth-jump -> one piece. */
    {
        int N = 11; size_t nz = (size_t)N*(size_t)N;
        float *z = (float *)ARENA_CALLOC(arena, (long)nz, (long)sizeof(float));
        float *v; int32_t *f; size_t nv, nf;
        make_grid(arena, N, N, z, &v, &nv, &f, &nf);
        ComponentMesh cm; wrap_mesh(&cm, v, nv, f, nf);
        ComponentMesh *out = NULL; size_t n = 0;
        int rc = DepthPeel_process(arena, &cm, min_gap, gap, minv, &out, &n);
        int ok = (rc == 0) && (n == 1);
        fprintf(stderr, "  [flat] rc=%d pieces=%zu -> %s\n", rc, n, ok?"ok":"FAIL");
        if (!ok) fails++;
    }

    /* (2) Gently curved single sheet (shallow parabola, max adjacent jump ~0.9 <
     * min_gap): curved but one wrap -> must NOT split. */
    {
        int N = 11; size_t nz = (size_t)N*(size_t)N;
        float *z = (float *)ARENA_ALLOC(arena, (long)(nz*sizeof(float)));
        for (int r = 0; r < N; r++) for (int c = 0; c < N; c++)
            z[(size_t)r*(size_t)N+(size_t)c] = 0.1f*(float)((c-5)*(c-5));
        float *v; int32_t *f; size_t nv, nf;
        make_grid(arena, N, N, z, &v, &nv, &f, &nf);
        ComponentMesh cm; wrap_mesh(&cm, v, nv, f, nf);
        ComponentMesh *out = NULL; size_t n = 0;
        int rc = DepthPeel_process(arena, &cm, min_gap, gap, minv, &out, &n);
        int ok = (rc == 0) && (n == 1);
        fprintf(stderr, "  [curve] rc=%d pieces=%zu -> %s\n", rc, n, ok?"ok":"FAIL");
        if (!ok) fails++;
    }

    /* (3) Two parallel sheets 3 vox apart, welded by a wall bridge: the bridge
     * edges jump 3 vox in depth -> cut -> two pieces, each ~one sheet. */
    {
        int N = 11; size_t per = (size_t)N*(size_t)N;
        size_t nv = 2*per;
        size_t cap_nf = 2*(size_t)(N-1)*(size_t)(N-1)*2 + 4;
        float   *v = (float *)ARENA_ALLOC(arena, (long)(nv*3*sizeof(float)));
        int32_t *f = (int32_t *)ARENA_ALLOC(arena, (long)(cap_nf*3*sizeof(int32_t)));
        /* sheet A at z=0 (verts 0..per-1), sheet B at z=3 (verts per..2per-1) */
        for (int s = 0; s < 2; s++) for (int r = 0; r < N; r++) for (int c = 0; c < N; c++) {
            size_t idx = (size_t)s*per + (size_t)r*(size_t)N + (size_t)c;
            v[idx*3+0] = (float)c;
            v[idx*3+1] = (float)r;
            v[idx*3+2] = (s==0) ? 0.0f : 3.0f;
        }
        size_t fi = 0;
        for (int s = 0; s < 2; s++) {
            int32_t base = (int32_t)((size_t)s*per);
            for (int r = 0; r < N-1; r++) for (int c = 0; c < N-1; c++) {
                int32_t a = base + (int32_t)((size_t)r*(size_t)N + (size_t)c);
                int32_t b = base + (int32_t)((size_t)r*(size_t)N + (size_t)c + 1);
                int32_t cc= base + (int32_t)((size_t)(r+1)*(size_t)N + (size_t)c);
                int32_t d = base + (int32_t)((size_t)(r+1)*(size_t)N + (size_t)c + 1);
                f[fi*3+0]=a; f[fi*3+1]=b; f[fi*3+2]=cc; fi++;
                f[fi*3+0]=b; f[fi*3+1]=d; f[fi*3+2]=cc; fi++;
            }
        }
        /* wall bridge: two triangles welding corner (0,0) of A to corner of B. */
        int32_t a0 = 0, a1 = 1, aN = (int32_t)N;
        int32_t b0 = (int32_t)per, b1 = (int32_t)per + 1;
        f[fi*3+0]=a0; f[fi*3+1]=b0; f[fi*3+2]=a1; fi++;
        f[fi*3+0]=b0; f[fi*3+1]=b1; f[fi*3+2]=a1; fi++;
        (void)aN;
        ComponentMesh cm; wrap_mesh(&cm, v, nv, f, fi);
        ComponentMesh *out = NULL; size_t n = 0;
        int rc = DepthPeel_process(arena, &cm, min_gap, gap, minv, &out, &n);
        int ok = (rc == 0) && (n == 2);
        if (ok) for (size_t i = 0; i < n; i++) if (out[i].nv < minv) ok = 0;
        fprintf(stderr, "  [stack] rc=%d pieces=%zu -> %s\n", rc, n, ok?"ok":"FAIL");
        if (!ok) fails++;
    }

    /* (4) Trivial input: empty mesh -> graceful, no crash. */
    {
        ComponentMesh cm; memset(&cm, 0, sizeof(cm)); cm.self = &cm;
        ComponentMesh *out = NULL; size_t n = 0;
        int rc = DepthPeel_process(arena, &cm, min_gap, gap, minv, &out, &n);
        int ok = (rc != 0) || (n <= 1);
        fprintf(stderr, "  [empty] rc=%d pieces=%zu -> %s\n", rc, n, ok?"ok":"FAIL");
        if (!ok) fails++;
    }

    if (fails == 0) fprintf(stderr, "[DepthPeel selftest] ok\n");
    else            fprintf(stderr, "[DepthPeel selftest] %d FAILURE(s)\n", fails);
    Arena_dispose(&arena);
    return fails;
}
