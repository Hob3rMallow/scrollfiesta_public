/*
 * dev_cut.c -- Developability cut (Stein/Grinspun/Crane 2018, sec 4.4). See
 * dev_cut.h. Pure splitter: threshold the per-vertex Crane energy on the raw
 * mesh, remove the seam band (+ geometric dilation), return the connected
 * survivors. Mirrors the cut-zone construction of split_mesh in bridge_cut.c,
 * but the barrier is the seam-vertex SET (lambda_v > eps), not a max-flow cut.
 */
#define _USE_MATH_DEFINES
#include "dev_cut.h"

#include "../common/kdtree.h"
#include "../common/union_find.h"
#include "../common/pca.h"
#include "../topology/developability.h"   /* Develop_vertex_energy */

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- boundary / used detection (undirected edge multiplicity) ------------- */

static int cmp_u64(const void *a, const void *b)
{ uint64_t x=*(const uint64_t*)a, y=*(const uint64_t*)b; return (x<y)?-1:(x>y)?1:0; }

/* is_boundary[v] = 1 if v touches an edge used by exactly one triangle.
 * used[v] = 1 if v appears in any face (may be NULL). Both caller-allocated[nv];
 * is_boundary is zeroed here, used must be pre-zeroed by the caller. */
static void mark_boundary(Arena_T arena, const int32_t *faces, size_t nf, size_t nv,
                          unsigned char *is_boundary, unsigned char *used)
{
    Arena_Mark mark = Arena_save(arena);
    size_t ne = nf*3;
    uint64_t *keys = (uint64_t *)ARENA_ALLOC(arena, (long)(ne*sizeof(uint64_t)));
    size_t m = 0;
    for (size_t t = 0; t < nf; t++) {
        int32_t v[3] = { faces[t*3+0], faces[t*3+1], faces[t*3+2] };
        if (used) { used[v[0]] = 1; used[v[1]] = 1; used[v[2]] = 1; }
        for (int k = 0; k < 3; k++) {
            int32_t a = v[k], b = v[(k+1)%3];
            int32_t lo = a<b?a:b, hi = a<b?b:a;
            keys[m++] = (uint64_t)lo * (uint64_t)nv + (uint64_t)hi;
        }
    }
    qsort(keys, m, sizeof(uint64_t), cmp_u64);
    for (size_t i = 0; i < nv; i++) is_boundary[i] = 0;
    size_t i = 0;
    while (i < m) {
        size_t j = i+1;
        while (j < m && keys[j] == keys[i]) j++;
        if (j - i == 1) {                          /* boundary edge */
            uint64_t key = keys[i];
            int32_t lo = (int32_t)(key / nv), hi = (int32_t)(key % nv);
            is_boundary[lo] = 1; is_boundary[hi] = 1;
        }
        i = j;
    }
    Arena_restore(arena, mark);
}

/* ---- sub-mesh construction ------------------------------------------------ */

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

/* ---- DevCut_process ------------------------------------------------------- */

int DevCut_process(Arena_T arena, const ComponentMesh *mesh,
                   double eps, double gap_depth, size_t min_comp_verts,
                   ComponentMesh **out_meshes, size_t *out_count)
{
    assert(arena && out_meshes && out_count);
    *out_meshes = NULL; *out_count = 0;
    if (!mesh || !mesh->verts || !mesh->faces) return -1;
    if (mesh->nv < 3 || mesh->nf < 1)
        return emit_single(arena, mesh, out_meshes, out_count);

    size_t nv = mesh->nv, nf = mesh->nf;
    Arena_Mark scratch = Arena_save(arena);

    /* 1. per-vertex Crane energy (boundary verts come back 0). */
    double *lam = (double *)ARENA_ALLOC(arena, (long)(nv*sizeof(double)));
    if (Develop_vertex_energy(arena, mesh->verts, nv, mesh->faces, nf, lam) != 0) {
        Arena_restore(arena, scratch);
        return emit_single(arena, mesh, out_meshes, out_count);
    }

    /* 2. seam vertices. */
    unsigned char *seam = (unsigned char *)ARENA_CALLOC(arena, (long)nv, 1L);
    size_t n_seam = 0;
    for (size_t v = 0; v < nv; v++) if (lam[v] > eps) { seam[v] = 1; n_seam++; }
    if (n_seam == 0) {                              /* developable -> no cut */
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
     * excluded). The removed band severs the mesh where the seam spanned it. */
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

/* ---- DevCut_nondev_fraction (gate metric) --------------------------------- */

double DevCut_nondev_fraction(Arena_T arena, const ComponentMesh *mesh,
                              double eps, double *out_mean_lambda)
{
    if (out_mean_lambda) *out_mean_lambda = 0.0;
    if (!mesh || !mesh->verts || !mesh->faces || mesh->nv < 3 || mesh->nf < 1)
        return 0.0;

    size_t nv = mesh->nv, nf = mesh->nf;
    Arena_Mark mark = Arena_save(arena);

    double *lam = (double *)ARENA_ALLOC(arena, (long)(nv*sizeof(double)));
    unsigned char *is_b = (unsigned char *)ARENA_ALLOC(arena, (long)nv);
    unsigned char *used = (unsigned char *)ARENA_CALLOC(arena, (long)nv, 1L);
    mark_boundary(arena, mesh->faces, nf, nv, is_b, used);
    if (Develop_vertex_energy(arena, mesh->verts, nv, mesh->faces, nf, lam) != 0) {
        Arena_restore(arena, mark);
        return 0.0;
    }

    size_t n_int = 0, n_bad = 0;
    double sum = 0.0;
    for (size_t v = 0; v < nv; v++) {
        if (!used[v] || is_b[v]) continue;          /* interior only */
        n_int++;
        sum += lam[v];
        if (lam[v] > eps) n_bad++;
    }
    double frac = n_int ? (double)n_bad / (double)n_int : 0.0;
    if (out_mean_lambda && n_int) *out_mean_lambda = sum / (double)n_int;

    Arena_restore(arena, mark);
    return frac;
}

/* ============================================================================
 * Self-test: a flat sheet self-gates (no cut), a dome is detected but does not
 * separate, a vertex-spanning non-developable seam splits into two more-
 * developable halves, and trivial input is a clean no-op.
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

int DevCut_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    fprintf(stderr, "=== DevCut selftest ===\n");

    const double eps  = 0.01;   /* test threshold (well below seam, above flat) */
    const double gap  = 0.5;    /* grid spacing 1.0 -> removes only seam verts   */
    const size_t minv = 10;

    /* (1) Flat sheet: developable -> no cut, fraction 0. */
    {
        int N = 11; size_t nz = (size_t)N*(size_t)N;
        float *z = (float *)ARENA_CALLOC(arena, (long)nz, (long)sizeof(float));
        float *v; int32_t *f; size_t nv, nf;
        make_grid(arena, N, N, z, &v, &nv, &f, &nf);
        ComponentMesh cm; wrap_mesh(&cm, v, nv, f, nf);

        double frac = DevCut_nondev_fraction(arena, &cm, eps, NULL);
        ComponentMesh *out = NULL; size_t n = 0;
        int rc = DevCut_process(arena, &cm, eps, gap, minv, &out, &n);
        int ok = (rc == 0) && (n == 1) && (frac < 1e-6);
        fprintf(stderr, "  [flat] rc=%d pieces=%zu frac=%.4f -> %s\n",
                rc, n, frac, ok?"ok":"FAIL");
        if (!ok) fails++;
    }

    /* (2) Cylinder: curved but DEVELOPABLE (a rolled sheet, like papyrus). Its
     * facet normals are all perpendicular to the axis, so they are coplanar ->
     * lambda ~ 0 at every vertex -> it must NOT be cut, however curved. This is
     * the guard that the splitter never severs a genuine scroll wrap. */
    {
        int nu = 13, nvr = 11; double R = 5.0;        /* half-cylinder grid */
        size_t nv = (size_t)nu*(size_t)nvr;
        size_t nf = (size_t)(nu-1)*(size_t)(nvr-1)*2;
        float   *v = (float *)ARENA_ALLOC(arena, (long)(nv*3*sizeof(float)));
        int32_t *f = (int32_t *)ARENA_ALLOC(arena, (long)(nf*3*sizeof(int32_t)));
        for (int r = 0; r < nvr; r++) for (int c = 0; c < nu; c++) {
            double u = M_PI * (double)c/(nu-1);
            size_t idx = (size_t)r*(size_t)nu + (size_t)c;
            v[idx*3+0] = (float)(R*cos(u));
            v[idx*3+1] = (float)(R*sin(u));
            v[idx*3+2] = (float)r;
        }
        size_t fi = 0;
        for (int r = 0; r < nvr-1; r++) for (int c = 0; c < nu-1; c++) {
            int32_t a  = (int32_t)((size_t)r*(size_t)nu + (size_t)c);
            int32_t b  = (int32_t)((size_t)r*(size_t)nu + (size_t)c + 1);
            int32_t cc = (int32_t)((size_t)(r+1)*(size_t)nu + (size_t)c);
            int32_t d  = (int32_t)((size_t)(r+1)*(size_t)nu + (size_t)c + 1);
            f[fi*3+0]=a; f[fi*3+1]=b; f[fi*3+2]=cc; fi++;
            f[fi*3+0]=b; f[fi*3+1]=d; f[fi*3+2]=cc; fi++;
        }
        ComponentMesh cm; wrap_mesh(&cm, v, nv, f, fi);

        double *lam = (double *)ARENA_ALLOC(arena, (long)(nv*sizeof(double)));
        Develop_vertex_energy(arena, v, nv, f, fi, lam);
        double lmax = 0.0; for (size_t i = 0; i < nv; i++) if (lam[i] > lmax) lmax = lam[i];
        double frac = DevCut_nondev_fraction(arena, &cm, eps, NULL);

        ComponentMesh *out = NULL; size_t n = 0;
        int rc = DevCut_process(arena, &cm, eps, gap, minv, &out, &n);
        int ok = (rc == 0) && (n == 1) && (lmax < eps) && (frac < 1e-6);
        fprintf(stderr, "  [cyl] rc=%d pieces=%zu maxlambda=%.6f frac=%.4f -> %s\n",
                rc, n, lmax, frac, ok?"ok":"FAIL");
        if (!ok) fails++;
    }

    /* (3) Spanning non-developable seam: a zig-zag middle row (alternating +/-amp)
     * makes every middle vertex's normals non-coplanar (high lambda) while the two
     * flat halves stay developable. The cut must split into 2 more-developable
     * pieces. */
    {
        int ncol = 11, nrow = 15, mid = 7; double amp = 1.0;
        size_t nz = (size_t)ncol*(size_t)nrow;
        float *z = (float *)ARENA_CALLOC(arena, (long)nz, (long)sizeof(float));
        for (int c = 0; c < ncol; c++)
            z[(size_t)mid*(size_t)ncol+(size_t)c] = (float)((c & 1) ? amp : -amp);
        float *v; int32_t *f; size_t nv, nf;
        make_grid(arena, ncol, nrow, z, &v, &nv, &f, &nf);
        ComponentMesh cm; wrap_mesh(&cm, v, nv, f, nf);

        double pfrac = DevCut_nondev_fraction(arena, &cm, eps, NULL);
        ComponentMesh *out = NULL; size_t n = 0;
        int rc = DevCut_process(arena, &cm, eps, gap, minv, &out, &n);
        int ok = (rc == 0) && (n == 2);
        if (ok) for (size_t i = 0; i < n; i++) {
            double cf = DevCut_nondev_fraction(arena, &out[i], eps, NULL);
            if (out[i].nv < minv) ok = 0;
            if (!(cf < pfrac)) ok = 0;
        }
        fprintf(stderr, "  [seam] rc=%d pieces=%zu parent_frac=%.4f -> %s\n",
                rc, n, pfrac, ok?"ok":"FAIL");
        if (!ok) fails++;
    }

    /* (4) Trivial input: empty mesh -> graceful, no crash. */
    {
        ComponentMesh cm; memset(&cm, 0, sizeof(cm)); cm.self = &cm;
        ComponentMesh *out = NULL; size_t n = 0;
        int rc = DevCut_process(arena, &cm, eps, gap, minv, &out, &n);
        int ok = (rc != 0) || (n <= 1);
        fprintf(stderr, "  [empty] rc=%d pieces=%zu -> %s\n", rc, n, ok?"ok":"FAIL");
        if (!ok) fails++;
    }

    if (fails == 0) fprintf(stderr, "[DevCut selftest] ok\n");
    else            fprintf(stderr, "[DevCut selftest] %d FAILURE(s)\n", fails);
    Arena_dispose(&arena);
    return fails;
}
