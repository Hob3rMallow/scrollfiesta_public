/* manifold_guard.c -- final per-cube manifold repair (see manifold_guard.h). */
#include "manifold_guard.h"

#include "pinhole_fill.h"      /* PinholeFill_split_pinches */
#include "orient_mesh.h"       /* OrientMesh_consistent */
#include "../common/mesh_manifold.h"

#include <stdlib.h>            /* qsort */
#include <string.h>           /* memset, memcpy */
#include <math.h>             /* sqrt */

/* ---- triangle area (verts are float (z,y,x); area is order-independent) ---- */
static double tri_area(const float *V, int32_t a, int32_t b, int32_t c) {
    double ux = (double)V[(size_t)b*3+0] - V[(size_t)a*3+0];
    double uy = (double)V[(size_t)b*3+1] - V[(size_t)a*3+1];
    double uz = (double)V[(size_t)b*3+2] - V[(size_t)a*3+2];
    double vx = (double)V[(size_t)c*3+0] - V[(size_t)a*3+0];
    double vy = (double)V[(size_t)c*3+1] - V[(size_t)a*3+1];
    double vz = (double)V[(size_t)c*3+2] - V[(size_t)a*3+2];
    double cx = uy*vz - uz*vy, cy = uz*vx - ux*vz, cz = ux*vy - uy*vx;
    return 0.5 * sqrt(cx*cx + cy*cy + cz*cz);
}

/* ---- directed half-edge carrying its face id and forward flag ---- */
typedef struct { int32_t lo, hi, face; int fwd; } HEF;
static int cmp_hef(const void *pa, const void *pb) {
    const HEF *a = (const HEF *)pa, *b = (const HEF *)pb;
    if (a->lo != b->lo) return (a->lo < b->lo) ? -1 : 1;
    if (a->hi != b->hi) return (a->hi < b->hi) ? -1 : 1;
    return 0;
}

/* Resolve >2-face edges: keep the largest-area forward face and the largest-area
 * backward face (a clean opposite-winding pair, or just one if all same
 * direction); delete the rest. Compacts cm->faces in place. */
static void resolve_nm_edges(Arena_T arena, ComponentMesh *cm,
                             size_t *io_nm, size_t *io_del) {
    size_t nf = cm->nf;
    if (nf == 0) return;
    int32_t *F = cm->faces;
    const float *V = cm->verts;

    double *area = (double *)ARENA_ALLOC(arena, (long)(nf * sizeof(double)));
    for (size_t f = 0; f < nf; f++)
        area[f] = tri_area(V, F[f*3+0], F[f*3+1], F[f*3+2]);

    size_t hn = nf * 3;
    HEF *he = (HEF *)ARENA_ALLOC(arena, (long)(hn * sizeof(HEF)));
    size_t k = 0;
    for (size_t f = 0; f < nf; f++) {
        int32_t v[3] = { F[f*3+0], F[f*3+1], F[f*3+2] };
        for (int e = 0; e < 3; e++) {
            int32_t s = v[e], d = v[(e+1)%3];
            HEF *h = &he[k++];
            if (s < d) { h->lo = s; h->hi = d; h->fwd = 1; }
            else       { h->lo = d; h->hi = s; h->fwd = 0; }
            h->face = (int32_t)f;
        }
    }
    qsort(he, hn, sizeof(HEF), cmp_hef);

    uint8_t *mark = (uint8_t *)ARENA_CALLOC(arena, (long)nf, (long)sizeof(uint8_t));
    size_t i = 0, nm = 0;
    while (i < hn) {
        size_t j = i + 1;
        while (j < hn && he[j].lo == he[i].lo && he[j].hi == he[i].hi) j++;
        size_t run = j - i;
        if (run > 2) {
            nm++;
            int32_t bestF = -1, bestB = -1;
            double aF = -1.0, aB = -1.0;
            for (size_t t = i; t < j; t++) {
                int32_t f = he[t].face;
                if (he[t].fwd) { if (area[f] > aF) { aF = area[f]; bestF = f; } }
                else           { if (area[f] > aB) { aB = area[f]; bestB = f; } }
            }
            for (size_t t = i; t < j; t++) {
                int32_t f = he[t].face;
                if (f != bestF && f != bestB) mark[f] = 1;
            }
        }
        i = j;
    }

    size_t w = 0, del = 0;
    for (size_t f = 0; f < nf; f++) {
        if (mark[f]) { del++; continue; }
        if (w != f) {
            F[w*3+0] = F[f*3+0]; F[w*3+1] = F[f*3+1]; F[w*3+2] = F[f*3+2];
        }
        w++;
    }
    cm->nf = w;
    if (io_nm)  *io_nm  += nm;
    if (io_del) *io_del += del;
}

int ManifoldGuard_process(Arena_T arena, ComponentMesh *meshes, size_t n_meshes,
                          int reorient, ManifoldGuardStats *out) {
    ManifoldGuardStats st;
    memset(&st, 0, sizeof st);

    /* 1. edge-manifold every component */
    for (size_t i = 0; i < n_meshes; i++) {
        if (meshes[i].nf == 0) continue;
        resolve_nm_edges(arena, &meshes[i], &st.nm_edges_resolved, &st.faces_deleted);
    }
    /* 2. vertex-manifold (split bowties), now that input is edge-manifold */
    PinholeFill_split_pinches(arena, meshes, n_meshes, &st.pinch_splits);
    /* 3. re-assert consistent winding (caller's choice -- resolve+split above
     * preserve winding, so a caller that already oriented can skip this). */
    if (reorient) {
        for (size_t i = 0; i < n_meshes; i++) {
            if (meshes[i].nf == 0) continue;
            size_t fl = 0, co = 0, rs = 0;
            OrientMesh_consistent(arena, meshes[i].verts, meshes[i].nv, NULL,
                                  meshes[i].faces, meshes[i].nf, &fl, &co, &rs);
            st.orient_flips += fl;
        }
    }

    if (out) *out = st;
    return 0;
}

/* ---- selftest ---- */
static ComponentMesh mk_mesh(Arena_T a, const float *V, size_t nv,
                             const int32_t *Fin, size_t nf) {
    ComponentMesh cm;
    memset(&cm, 0, sizeof cm);
    cm.verts = (float *)ARENA_ALLOC(a, (long)(nv * 3 * sizeof(float)));
    memcpy(cm.verts, V, nv * 3 * sizeof(float));
    cm.faces = (int32_t *)ARENA_ALLOC(a, (long)(nf * 3 * sizeof(int32_t)));
    memcpy(cm.faces, Fin, nf * 3 * sizeof(int32_t));
    cm.nv = nv; cm.nf = nf; cm.comp_id = 1;
    return cm;
}

int ManifoldGuard_selftest(Arena_T arena) {
    int fails = 0;

    /* case A: non-manifold edge {0,1} shared by 3 faces (2 fwd, 1 bwd) ->
     * resolve keeps a clean opposite-winding pair, drops the third. */
    {
        float V[] = { 0,0,0,  0,0,2,  0,1,0,  0,-1,0,  0,2,0 };
        int32_t F[] = { 0,1,2,  1,0,3,  0,1,4 };
        ComponentMesh cm = mk_mesh(arena, V, 5, F, 3); cm.self = &cm;
        ManifoldGuardStats st;
        ManifoldGuard_process(arena, &cm, 1, 1, &st);
        MeshManifoldStats ms = MeshManifold_audit(arena, cm.nv, cm.faces, cm.nf);
        if (!(MeshManifold_ok(&ms) && cm.nf < 3 && st.nm_edges_resolved == 1)) fails++;
    }
    /* case B: bowtie vertex 0 shared by two edge-disjoint fans -> split. */
    {
        float V[] = { 0,0,0,  0,0,1,  0,1,0,  1,0,0,  1,0,1 };
        int32_t F[] = { 0,1,2,  0,3,4 };
        ComponentMesh cm = mk_mesh(arena, V, 5, F, 2); cm.self = &cm;
        MeshManifoldStats pre = MeshManifold_audit(arena, cm.nv, cm.faces, cm.nf);
        ManifoldGuardStats st;
        ManifoldGuard_process(arena, &cm, 1, 1, &st);
        MeshManifoldStats ms = MeshManifold_audit(arena, cm.nv, cm.faces, cm.nf);
        if (!(pre.nm_verts == 1 && MeshManifold_ok(&ms) &&
              st.pinch_splits >= 1 && cm.nv > 5)) fails++;
    }
    /* case C: clean two-triangle disk -> untouched, still manifold. */
    {
        float V[] = { 0,0,0,  0,0,1,  0,1,0,  0,1,1 };
        int32_t F[] = { 0,1,2,  1,3,2 };
        ComponentMesh cm = mk_mesh(arena, V, 4, F, 2); cm.self = &cm;
        ManifoldGuardStats st;
        ManifoldGuard_process(arena, &cm, 1, 1, &st);
        MeshManifoldStats ms = MeshManifold_audit(arena, cm.nv, cm.faces, cm.nf);
        if (!(MeshManifold_ok(&ms) && cm.nf == 2 &&
              st.nm_edges_resolved == 0 && st.pinch_splits == 0)) fails++;
    }

    return fails;
}
