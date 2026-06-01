#include "orient_mesh.h"

#include <string.h>

/* Recompute shared-edge consistency from LIVE face arrays at visit time so
 * flips are never stale; anchor each component's global sign to the MLS
 * normal field. See orient_mesh.h for the full rationale. */

static size_t next_pow2(size_t n) {
    size_t p = 1;
    while (p < n) { p <<= 1; }
    return p;
}

static uint64_t edge_hash(int32_t lo, int32_t hi) {
    uint64_t h = (uint64_t)(uint32_t)lo * 0x9E3779B97F4A7C15ULL
               + (uint64_t)(uint32_t)hi;
    h ^= h >> 33; h *= 0xff51afd7ed558ccdULL; h ^= h >> 33;
    return h;
}

int OrientMesh_consistent(Arena_T arena,
                          const float *verts, size_t nv,
                          const float *normals,
                          int32_t *faces, size_t nf,
                          size_t *out_flips,
                          size_t *out_components,
                          size_t *out_residual_same_dir) {
    (void)nv;
    if (out_flips) { *out_flips = 0; }
    if (out_components) { *out_components = 0; }
    if (out_residual_same_dir) { *out_residual_same_dir = 0; }
    if (nf == 0 || faces == NULL) { return 0; }

    Arena_Mark mark = Arena_save(arena);

    /* ---- Phase A: build manifold edge -> face adjacency ---- */
    size_t ne_max = nf * 3;
    size_t hsz = next_pow2(ne_max * 2);
    if (hsz < 1024) { hsz = 1024; }
    uint64_t hmask = (uint64_t)(hsz - 1);

    int *bk = (int *)ARENA_ALLOC(arena, (long)(hsz * sizeof(int)));
    memset(bk, 0xFF, hsz * sizeof(int)); /* -1 in every slot */

    typedef struct { int32_t lo, hi, f0, f1, next; } Edge;
    Edge *ed = (Edge *)ARENA_ALLOC(arena, (long)(ne_max * sizeof(Edge)));
    size_t en = 0;

    for (size_t f = 0; f < nf; f++) {
        int32_t v[3] = { faces[f * 3 + 0], faces[f * 3 + 1], faces[f * 3 + 2] };
        for (int e = 0; e < 3; e++) {
            int32_t a = v[e], b = v[(e + 1) % 3];
            int32_t lo = (a < b) ? a : b;
            int32_t hi = (a < b) ? b : a;
            size_t bb = (size_t)(edge_hash(lo, hi) & hmask);
            int idx = bk[bb];
            while (idx >= 0) {
                if (ed[idx].lo == lo && ed[idx].hi == hi) { break; }
                idx = ed[idx].next;
            }
            if (idx < 0) {
                idx = (int)en++;
                ed[idx].lo = lo; ed[idx].hi = hi;
                ed[idx].f0 = (int32_t)f; ed[idx].f1 = -1;
                ed[idx].next = bk[bb]; bk[bb] = idx;
            } else if (ed[idx].f1 == -1) {
                ed[idx].f1 = (int32_t)f;
            } else {
                ed[idx].f1 = -2; /* >2 faces: non-manifold edge, not crossed */
            }
        }
    }

    /* CSR face adjacency over exactly-2-face (manifold) edges */
    int *deg = (int *)ARENA_CALLOC(arena, (long)nf, (long)sizeof(int));
    for (size_t i = 0; i < en; i++) {
        if (ed[i].f1 >= 0) { deg[ed[i].f0]++; deg[ed[i].f1]++; }
    }
    size_t *off = (size_t *)ARENA_ALLOC(arena,
                                        (long)((nf + 1) * sizeof(size_t)));
    off[0] = 0;
    for (size_t i = 0; i < nf; i++) { off[i + 1] = off[i] + (size_t)deg[i]; }
    size_t adj_n = off[nf];
    int *adjf = (int *)ARENA_ALLOC(arena, (long)(adj_n * sizeof(int)));
    size_t *cur = (size_t *)ARENA_ALLOC(arena, (long)(nf * sizeof(size_t)));
    memcpy(cur, off, nf * sizeof(size_t));
    for (size_t i = 0; i < en; i++) {
        if (ed[i].f1 >= 0) {
            adjf[cur[ed[i].f0]++] = ed[i].f1;
            adjf[cur[ed[i].f1]++] = ed[i].f0;
        }
    }

    /* ---- Phase B: BFS orient (live windings); Phase C: anchor sign ---- */
    uint8_t *vis = (uint8_t *)ARENA_CALLOC(arena, (long)nf, 1);
    int *q = (int *)ARENA_ALLOC(arena, (long)(nf * sizeof(int)));
    size_t flips = 0, comps = 0;

    for (size_t s = 0; s < nf; s++) {
        if (vis[s]) { continue; }
        comps++;
        size_t qh = 0, qt = 0, comp_flips = 0;
        q[qt++] = (int)s; vis[s] = 1;
        while (qh < qt) {
            int fc = q[qh++];
            int32_t a = faces[fc * 3 + 0];
            int32_t b = faces[fc * 3 + 1];
            int32_t c = faces[fc * 3 + 2];
            int32_t pe[3][2] = { { a, b }, { b, c }, { c, a } };
            for (size_t k = off[fc]; k < off[(size_t)fc + 1]; k++) {
                int fn = adjf[k];
                if (vis[fn]) { continue; }
                int32_t na = faces[fn * 3 + 0];
                int32_t nb = faces[fn * 3 + 1];
                int32_t nc = faces[fn * 3 + 2];
                int32_t nedg[3][2] = { { na, nb }, { nb, nc }, { nc, na } };
                int same = 0;
                for (int pi = 0; pi < 3 && !same; pi++) {
                    for (int ni = 0; ni < 3; ni++) {
                        if (pe[pi][0] == nedg[ni][0] &&
                            pe[pi][1] == nedg[ni][1]) { same = 1; break; }
                    }
                }
                if (same) {
                    int32_t t = faces[fn * 3 + 1];
                    faces[fn * 3 + 1] = faces[fn * 3 + 2];
                    faces[fn * 3 + 2] = t;
                    flips++; comp_flips++;
                }
                vis[fn] = 1;
                q[qt++] = fn;
            }
        }

        /* Phase C: pick the component's global sign.
         *
         * normals != NULL (e.g. post-BPA): anchor to the MLS normal field —
         * a world-coordinate property, so neighbouring cubes that share the
         * geometry choose the same sign (needed for the seam weld).
         *
         * normals == NULL (e.g. post-QEM, no per-vertex normals left): keep
         * the orientation of the MAJORITY of the component's original faces.
         * The BFS made everything agree with face `s`; if that meant flipping
         * more than half the faces, then `s` itself was the minority (the few
         * faces a collapse mis-wound) and the whole component is now inverted
         * — flip it back. */
        if (normals == NULL) {
            if (comp_flips * 2 > qt) {
                for (size_t qi = 0; qi < qt; qi++) {
                    int fc = q[qi];
                    int32_t t = faces[fc * 3 + 1];
                    faces[fc * 3 + 1] = faces[fc * 3 + 2];
                    faces[fc * 3 + 2] = t;
                    flips++;
                }
            }
        } else {
            double align = 0.0;
            for (size_t qi = 0; qi < qt; qi++) {
                int fc = q[qi];
                int32_t fa = faces[fc * 3 + 0];
                int32_t fb = faces[fc * 3 + 1];
                int32_t fcc = faces[fc * 3 + 2];
                float e1[3] = { verts[fb * 3 + 0] - verts[fa * 3 + 0],
                                verts[fb * 3 + 1] - verts[fa * 3 + 1],
                                verts[fb * 3 + 2] - verts[fa * 3 + 2] };
                float e2[3] = { verts[fcc * 3 + 0] - verts[fa * 3 + 0],
                                verts[fcc * 3 + 1] - verts[fa * 3 + 1],
                                verts[fcc * 3 + 2] - verts[fa * 3 + 2] };
                float gn[3] = { e1[1] * e2[2] - e1[2] * e2[1],
                                e1[2] * e2[0] - e1[0] * e2[2],
                                e1[0] * e2[1] - e1[1] * e2[0] };
                float vn[3] = {
                    normals[fa * 3 + 0] + normals[fb * 3 + 0] + normals[fcc * 3 + 0],
                    normals[fa * 3 + 1] + normals[fb * 3 + 1] + normals[fcc * 3 + 1],
                    normals[fa * 3 + 2] + normals[fb * 3 + 2] + normals[fcc * 3 + 2] };
                align += (double)(gn[0] * vn[0] + gn[1] * vn[1] + gn[2] * vn[2]);
            }
            if (align < 0.0) {
                for (size_t qi = 0; qi < qt; qi++) {
                    int fc = q[qi];
                    int32_t t = faces[fc * 3 + 1];
                    faces[fc * 3 + 1] = faces[fc * 3 + 2];
                    faces[fc * 3 + 2] = t;
                    flips++;
                }
            }
        }
    }

    /* ---- residual same-direction interior-edge count (diagnostic) ---- */
    if (out_residual_same_dir) {
        memset(bk, 0xFF, hsz * sizeof(int));
        size_t en2 = 0, same = 0;
        for (size_t f = 0; f < nf; f++) {
            int32_t v[3] = { faces[f * 3 + 0], faces[f * 3 + 1],
                             faces[f * 3 + 2] };
            for (int e = 0; e < 3; e++) {
                int32_t a = v[e], b = v[(e + 1) % 3];
                int32_t lo = (a < b) ? a : b;
                int32_t hi = (a < b) ? b : a;
                int32_t dir = (a < b) ? 1 : 0;
                size_t bb = (size_t)(edge_hash(lo, hi) & hmask);
                int idx = bk[bb];
                while (idx >= 0) {
                    if (ed[idx].lo == lo && ed[idx].hi == hi) { break; }
                    idx = ed[idx].next;
                }
                if (idx < 0) {
                    idx = (int)en2++;
                    ed[idx].lo = lo; ed[idx].hi = hi;
                    ed[idx].f0 = dir;  /* dir0 */
                    ed[idx].f1 = 1;    /* count */
                    ed[idx].next = bk[bb]; bk[bb] = idx;
                } else {
                    ed[idx].f1++;
                    if (ed[idx].f1 == 2 && ed[idx].f0 == dir) { same++; }
                }
            }
        }
        *out_residual_same_dir = same;
    }

    if (out_flips) { *out_flips = flips; }
    if (out_components) { *out_components = comps; }

    Arena_restore(arena, mark);
    return 0;
}
