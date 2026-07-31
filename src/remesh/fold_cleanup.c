#include "fold_cleanup.h"

#include <math.h>
#include <string.h>
#include <stdint.h>

/* See fold_cleanup.h for the algorithm. */

/* For a same-dir edge with ALIGNED normals (dot>=0, a winding defect rather
 * than a geometric fold), only delete a face if the worse of the two is a
 * clear low-coherence outlier flap below this threshold. At/above it the patch
 * is a balanced non-orientable knot and is left untouched. Folds (dot<0) are
 * always cut regardless of this value. */
#define FOLD_OUTLIER_COHER 0.5

static void face_normal(const float *V, const int32_t *F, size_t f,
                        double n[3], double *area) {
    int32_t ia = F[f * 3 + 0], ib = F[f * 3 + 1], ic = F[f * 3 + 2];
    double e1[3] = { (double)V[ib*3+0] - V[ia*3+0],
                     (double)V[ib*3+1] - V[ia*3+1],
                     (double)V[ib*3+2] - V[ia*3+2] };
    double e2[3] = { (double)V[ic*3+0] - V[ia*3+0],
                     (double)V[ic*3+1] - V[ia*3+1],
                     (double)V[ic*3+2] - V[ia*3+2] };
    n[0] = e1[1]*e2[2] - e1[2]*e2[1];
    n[1] = e1[2]*e2[0] - e1[0]*e2[2];
    n[2] = e1[0]*e2[1] - e1[1]*e2[0];
    double len = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    *area = 0.5 * len;
    if (len > 1e-20) { n[0] /= len; n[1] /= len; n[2] /= len; }
}

/* Build vertex -> incident-face CSR (off[nv+1], inc[3*nf]). */
static void build_vert_faces(Arena_T arena, const int32_t *faces, size_t nf,
                             size_t nv, size_t **out_off, int **out_inc) {
    size_t *off = (size_t *)ARENA_CALLOC(arena, (size_t)(nv + 1),
                                         (long)sizeof(size_t));
    for (size_t f = 0; f < nf; f++) {
        for (int k = 0; k < 3; k++) { off[(size_t)faces[f * 3 + k] + 1]++; }
    }
    for (size_t v = 0; v < nv; v++) { off[v + 1] += off[v]; }
    int *inc = (int *)ARENA_ALLOC(arena, (size_t)(off[nv] * sizeof(int)));
    size_t *cur = (size_t *)ARENA_ALLOC(arena, (size_t)(nv * sizeof(size_t)));
    memcpy(cur, off, nv * sizeof(size_t));
    for (size_t f = 0; f < nf; f++) {
        for (int k = 0; k < 3; k++) {
            int32_t v = faces[f * 3 + k];
            inc[cur[v]++] = (int)f;
        }
    }
    *out_off = off; *out_inc = inc;
}

/* Mean dot of face f's normal with each edge-neighbour's normal, EXCLUDING the
 * partner across edge {pu,pv}. Higher = f sits more coherently in its fan. */
static double coherence(const float *V, const int32_t *F, size_t f,
                        int32_t pu, int32_t pv,
                        const size_t *off, const int *inc) {
    double nf_n[3], a; face_normal(V, F, f, nf_n, &a);
    int32_t plo = pu < pv ? pu : pv, phi = pu < pv ? pv : pu;
    double sum = 0.0; int cnt = 0;
    for (int e = 0; e < 3; e++) {
        int32_t u = F[f * 3 + e], v = F[f * 3 + (e + 1) % 3];
        int32_t lo = u < v ? u : v, hi = u < v ? v : u;
        if (lo == plo && hi == phi) { continue; }    /* skip the fold partner */
        for (size_t k = off[lo]; k < off[lo + 1]; k++) {
            size_t g = (size_t)inc[k];
            if (g == f) { continue; }
            int32_t ga = F[g*3], gb = F[g*3+1], gc = F[g*3+2];
            int has_hi = (ga == hi || gb == hi || gc == hi);
            if (!has_hi) { continue; }                /* g must share {lo,hi} */
            double gn[3], ga2; face_normal(V, F, g, gn, &ga2);
            sum += nf_n[0]*gn[0] + nf_n[1]*gn[1] + nf_n[2]*gn[2];
            cnt++;
        }
    }
    return cnt ? sum / (double)cnt : 0.0;
}

/* One pass over a component: mark the outlier face of every fold edge for
 * deletion, then compact. Returns faces removed this pass. */
static size_t fold_pass(Arena_T arena, ComponentMesh *cm) {
    size_t nv = cm->nv, nf = cm->nf;
    int32_t *faces = cm->faces;
    const float *V = cm->verts;

    size_t *off = NULL; int *inc = NULL;
    build_vert_faces(arena, faces, nf, nv, &off, &inc);

    uint8_t *del = (uint8_t *)ARENA_CALLOC(arena, (size_t)nf, 1);
    size_t marked = 0;

    for (size_t f = 0; f < nf; f++) {
        for (int e = 0; e < 3; e++) {
            int32_t u = faces[f * 3 + e], v = faces[f * 3 + (e + 1) % 3];
            /* find a DISTINCT face g (g>f) that traverses u->v the same way */
            for (size_t k = off[u]; k < off[u + 1]; k++) {
                size_t g = (size_t)inc[k];
                if (g <= f) { continue; }
                int same = 0;
                for (int e2 = 0; e2 < 3; e2++) {
                    if (faces[g*3+e2] == u && faces[g*3+(e2+1)%3] == v) {
                        same = 1; break;
                    }
                }
                if (!same) { continue; }
                double nf1[3], nf2[3], a1, a2;
                face_normal(V, faces, f, nf1, &a1);
                face_normal(V, faces, g, nf2, &a2);
                double dot = nf1[0]*nf2[0] + nf1[1]*nf2[1] + nf1[2]*nf2[2];
                if (del[f] || del[g]) { continue; }   /* already handled */
                double cf = coherence(V, faces, f, u, v, off, inc);
                double cg = coherence(V, faces, g, u, v, off, inc);
                double vcoher = (cf <= cg) ? cf : cg;
                /* A same-dir edge with opposing normals (dot<0) is always a
                 * fold -> drop the worse face. With aligned normals (dot>=0)
                 * it is a winding defect BFS could not flip away; only drop a
                 * face if it is a clear low-coherence outlier flap, otherwise
                 * the patch is a balanced non-orientable knot (leave it). */
                if (dot >= 0.0 && vcoher >= FOLD_OUTLIER_COHER) { continue; }
                size_t victim = (cf <= cg) ? f : g;
                del[victim] = 1; marked++;
            }
        }
    }

    if (marked == 0) { return 0; }

    /* compact faces in place (new arena buffer; old is past a mark we keep) */
    size_t nf_new = nf - marked;
    int32_t *nf_faces = (int32_t *)ARENA_ALLOC(arena,
                              (long)(nf_new * 3 * sizeof(int32_t)));
    size_t w = 0;
    for (size_t f = 0; f < nf; f++) {
        if (del[f]) { continue; }
        nf_faces[w*3+0] = faces[f*3+0];
        nf_faces[w*3+1] = faces[f*3+1];
        nf_faces[w*3+2] = faces[f*3+2];
        w++;
    }
    cm->faces = nf_faces;
    cm->nf = nf_new;
    return marked;
}

int FoldCleanup_process(Arena_T arena,
                        ComponentMesh *meshes, size_t n_meshes,
                        int max_passes,
                        size_t *out_faces_removed) {
    size_t tot = 0;
    if (max_passes < 1) { max_passes = 1; }
    for (size_t i = 0; i < n_meshes; i++) {
        ComponentMesh *cm = &meshes[i];
        if (cm->nf < 1 || cm->nv < 3) { continue; }
        for (int p = 0; p < max_passes; p++) {
            size_t r = fold_pass(arena, cm);
            tot += r;
            if (r == 0) { break; }
        }
    }
    if (out_faces_removed) { *out_faces_removed = tot; }
    return 0;
}
