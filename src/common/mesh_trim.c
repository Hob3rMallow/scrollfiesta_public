#include "mesh_trim.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int Mesh_trim_to_owned_box(Arena_T arena,
                           const float *verts_in, size_t nv_in,
                           const int32_t *faces_in, size_t nf_in,
                           const uint8_t *pin_mask_in,
                           float owned_lo, float owned_hi,
                           float **out_verts, size_t *out_nv,
                           int32_t **out_faces, size_t *out_nf,
                           uint8_t **out_pin_mask)
{
    assert(arena);
    assert(verts_in && faces_in);
    assert(out_verts && out_nv && out_faces && out_nf);

    *out_verts = NULL;
    *out_nv = 0;
    *out_faces = NULL;
    *out_nf = 0;
    if (out_pin_mask) *out_pin_mask = NULL;
    if (nf_in == 0 || nv_in == 0) return 0;

    /* Pass 1: keep faces with ALL THREE vertices inside [owned_lo, owned_hi).
     * Strict (all-vertex), not centroid: a centroid test let a vertex poke
     * ~edge/2 past the boundary, which shrank the inter-cube gap at an insetted
     * seam back into z-fight range (~0.5 vox). All-vertex containment makes the
     * inset exact, so adjacent cubes end up a clean 2*inset apart. */
    int32_t *kept_faces = (int32_t *)ARENA_ALLOC(arena,
                              (long)nf_in * 3L * (long)sizeof(int32_t));
    size_t nf_kept = 0;
    for (size_t f = 0; f < nf_in; f++) {
        int32_t a = faces_in[f * 3 + 0];
        int32_t b = faces_in[f * 3 + 1];
        int32_t c = faces_in[f * 3 + 2];
        if ((size_t)a >= nv_in || (size_t)b >= nv_in || (size_t)c >= nv_in) {
            continue;  /* bogus index, skip defensively */
        }
        int32_t tri[3] = { a, b, c };
        int inside = 1;
        for (int t = 0; t < 3; t++) {
            float vz = verts_in[tri[t]*3+0];
            float vy = verts_in[tri[t]*3+1];
            float vx = verts_in[tri[t]*3+2];
            if (vz < owned_lo || vz >= owned_hi ||
                vy < owned_lo || vy >= owned_hi ||
                vx < owned_lo || vx >= owned_hi) { inside = 0; break; }
        }
        if (!inside) continue;
        kept_faces[nf_kept * 3 + 0] = a;
        kept_faces[nf_kept * 3 + 1] = b;
        kept_faces[nf_kept * 3 + 2] = c;
        nf_kept++;
    }
    if (nf_kept == 0) return 0;

    /* Pass 2: mark referenced vertices. */
    uint8_t *used = (uint8_t *)ARENA_CALLOC(arena, (long)nv_in, 1L);
    for (size_t f = 0; f < nf_kept; f++) {
        used[kept_faces[f * 3 + 0]] = 1;
        used[kept_faces[f * 3 + 1]] = 1;
        used[kept_faces[f * 3 + 2]] = 1;
    }

    /* Pass 3: build remap and compact verts. */
    int32_t *remap = (int32_t *)ARENA_ALLOC(arena,
                         (long)nv_in * (long)sizeof(int32_t));
    size_t nv_kept = 0;
    for (size_t v = 0; v < nv_in; v++) {
        if (used[v]) {
            remap[v] = (int32_t)nv_kept;
            nv_kept++;
        } else {
            remap[v] = -1;
        }
    }

    float *new_verts = (float *)ARENA_ALLOC(arena,
                          (long)nv_kept * 3L * (long)sizeof(float));
    size_t wi = 0;
    for (size_t v = 0; v < nv_in; v++) {
        if (used[v]) {
            new_verts[wi * 3 + 0] = verts_in[v * 3 + 0];
            new_verts[wi * 3 + 1] = verts_in[v * 3 + 1];
            new_verts[wi * 3 + 2] = verts_in[v * 3 + 2];
            wi++;
        }
    }

    int32_t *new_faces = (int32_t *)ARENA_ALLOC(arena,
                            (long)nf_kept * 3L * (long)sizeof(int32_t));
    for (size_t f = 0; f < nf_kept; f++) {
        new_faces[f * 3 + 0] = remap[kept_faces[f * 3 + 0]];
        new_faces[f * 3 + 1] = remap[kept_faces[f * 3 + 1]];
        new_faces[f * 3 + 2] = remap[kept_faces[f * 3 + 2]];
    }

    *out_verts = new_verts;
    *out_nv = nv_kept;
    *out_faces = new_faces;
    *out_nf = nf_kept;

    /* Pin mask: pure index remap (no vert position changes). */
    if (pin_mask_in && out_pin_mask) {
        uint8_t *new_pins = (uint8_t *)ARENA_CALLOC(arena,
                                (long)nv_kept, 1L);
        for (size_t v = 0; v < nv_in; v++) {
            if (used[v] && pin_mask_in[v]) {
                new_pins[remap[v]] = 1;
            }
        }
        *out_pin_mask = new_pins;
    }

    return 0;
}

/* ================= cut-at-plane trim ================= */
/* Scratch is malloc/free (local to one call, seam_weld convention); only the
 * returned arrays are arena-allocated. */

/* One corner of the polygon being clipped. `edge_u/edge_v` name the ORIGINAL
 * mesh edge this corner lies on (-1 = none): an original vertex lies on
 * "itself" (edge_u == edge_v == its index); a cut point on an original edge
 * carries that edge (u < v); a cut point on a clip-generated segment carries
 * (-1,-1) and is face-private. This support tracking is what lets a second
 * plane cutting the SAME original edge still share its cut vertex with the
 * neighboring face (watertight cut boundary). */
typedef struct {
    double  p[3];
    int32_t vid;               /* output vertex index (>=0 once emitted)   */
    int32_t edge_u, edge_v;    /* original-edge support, or (-1,-1)        */
} CutCorner;

/* Open-addressing hash: (original edge u<v, plane 0..5) -> output vert id. */
typedef struct { int32_t u, v; int8_t plane; int32_t vid; } CutHashEnt;

typedef struct {
    CutHashEnt *ent;
    size_t      cap;           /* power of two */
} CutHash;

static size_t cuthash_slot(const CutHash *h, int32_t u, int32_t v, int plane)
{
    uint64_t k = (uint64_t)(uint32_t)u * 73856093u
               ^ (uint64_t)(uint32_t)v * 19349663u
               ^ (uint64_t)(unsigned)plane * 83492791u;
    k ^= k >> 33; k *= 0xff51afd7ed558ccdULL; k ^= k >> 33;
    return (size_t)(k & (h->cap - 1));
}

/* Find-or-reserve the cut vertex for (edge u<v, plane). Returns the slot;
 * caller checks ent[slot].vid < 0 to know it must create the vertex. */
static size_t cuthash_find(CutHash *h, int32_t u, int32_t v, int plane)
{
    size_t s = cuthash_slot(h, u, v, plane);
    for (;;) {
        CutHashEnt *e = &h->ent[s];
        if (e->u == -1) {                       /* empty: reserve */
            e->u = u; e->v = v; e->plane = (int8_t)plane; e->vid = -1;
            return s;
        }
        if (e->u == u && e->v == v && e->plane == (int8_t)plane) return s;
        s = (s + 1) & (h->cap - 1);
    }
}

/* Growable output vertex array (doubles for construction accuracy). */
typedef struct { double *p; size_t n, cap; } VGrow;

static int vgrow_push(VGrow *g, const double p[3])
{
    if (g->n == g->cap) {
        size_t nc = g->cap ? g->cap * 2 : 1024;
        double *np = (double *)realloc(g->p, nc * 3 * sizeof(double));
        if (!np) return -1;
        g->p = np; g->cap = nc;
    }
    g->p[g->n*3+0] = p[0]; g->p[g->n*3+1] = p[1]; g->p[g->n*3+2] = p[2];
    g->n++;
    return 0;
}

/* The original-edge support of the polygon segment A->B (u<v), or 0 if none.
 * Cases: both corners orig verts of one edge; one orig vert + a cut point on
 * an edge containing it; two cut points on the same edge. */
static int segment_support(const CutCorner *a, const CutCorner *b,
                           int32_t *out_u, int32_t *out_v)
{
    int32_t u = -1, v = -1;
    int a_orig = (a->edge_u >= 0 && a->edge_u == a->edge_v);
    int b_orig = (b->edge_u >= 0 && b->edge_u == b->edge_v);
    if (a_orig && b_orig) {
        u = a->edge_u; v = b->edge_u;
    } else if (a_orig && b->edge_u >= 0) {
        if (a->edge_u == b->edge_u || a->edge_u == b->edge_v) {
            u = b->edge_u; v = b->edge_v;
        }
    } else if (b_orig && a->edge_u >= 0) {
        if (b->edge_u == a->edge_u || b->edge_u == a->edge_v) {
            u = a->edge_u; v = a->edge_v;
        }
    } else if (a->edge_u >= 0 && b->edge_u >= 0 &&
               a->edge_u == b->edge_u && a->edge_v == b->edge_v) {
        u = a->edge_u; v = a->edge_v;
    }
    if (u < 0 || v < 0 || u == v) return 0;
    if (u > v) { int32_t t = u; u = v; v = t; }
    *out_u = u; *out_v = v;
    return 1;
}

int Mesh_trim_cut_to_owned_box(Arena_T arena,
                               const float *verts_in, size_t nv_in,
                               const int32_t *faces_in, size_t nf_in,
                               float owned_lo, float owned_hi,
                               float snap_eps,
                               float **out_verts, size_t *out_nv,
                               int32_t **out_faces, size_t *out_nf,
                               size_t *out_faces_cut)
{
    assert(arena);
    assert(out_verts && out_nv && out_faces && out_nf);
    assert(owned_hi > owned_lo);

    *out_verts = NULL; *out_nv = 0;
    *out_faces = NULL; *out_nf = 0;
    if (out_faces_cut) *out_faces_cut = 0;
    if (nf_in == 0 || nv_in == 0) return 0;
    if (verts_in == NULL || faces_in == NULL) return -1;

    const double lo = (double)owned_lo, hi = (double)owned_hi;
    const double eps = (double)snap_eps;

    /* Pass 0: snapped double copy of the input verts. Any coordinate within
     * snap_eps of a box plane is clamped exactly onto it, so borderline
     * jitter (CVT's double->float rounding) can no longer flip a vertex to
     * the outside, and grazing cuts pre-degenerate instead of minting
     * sub-eps strips. All faces see the same snapped coordinates. */
    VGrow vg = { NULL, 0, 0 };
    vg.p = (double *)malloc((nv_in ? nv_in : 1) * 3 * sizeof(double));
    if (!vg.p) return -1;
    vg.n = nv_in; vg.cap = nv_in ? nv_in : 1;
    for (size_t v = 0; v < nv_in; v++) {
        for (int a = 0; a < 3; a++) {
            double c = (double)verts_in[v*3+(size_t)a];
            if (fabs(c - lo) <= eps) c = lo;
            else if (fabs(c - hi) <= eps) c = hi;
            vg.p[v*3+(size_t)a] = c;
        }
    }

    /* Cut-vertex hash: worst case ~2 shared cuts per input edge; size to the
     * next power of two >= 4 * 3 * nf (load factor <= ~0.5). */
    CutHash hash;
    hash.cap = 64;
    while (hash.cap < nf_in * 12) hash.cap <<= 1;
    hash.ent = (CutHashEnt *)malloc(hash.cap * sizeof(CutHashEnt));
    if (!hash.ent) { free(vg.p); return -1; }
    for (size_t i = 0; i < hash.cap; i++) hash.ent[i].u = -1;

    /* Output face accumulator (malloc, moved to arena at the end). A clipped
     * triangle fans into at most 7 tris (9-gon). */
    size_t   fcap = nf_in * 2 + 64;
    size_t   nf_out = 0;
    int32_t *fout = (int32_t *)malloc(fcap * 3 * sizeof(int32_t));
    if (!fout) { free(vg.p); free(hash.ent); return -1; }

    size_t faces_cut = 0;
    CutCorner poly[16], next[16];
    int rc = 0;

    for (size_t f = 0; f < nf_in && rc == 0; f++) {
        int32_t t[3] = { faces_in[f*3+0], faces_in[f*3+1], faces_in[f*3+2] };
        if ((size_t)t[0] >= nv_in || (size_t)t[1] >= nv_in ||
            (size_t)t[2] >= nv_in) continue;   /* bogus index, skip */

        /* Trivial classification against all 6 planes at once. */
        int all_in = 1;
        int out_mask = 0x3F;   /* per-plane "all verts outside" bits */
        for (int k = 0; k < 3; k++) {
            const double *p = &vg.p[(size_t)t[k]*3];
            int m = 0;
            for (int a = 0; a < 3; a++) {
                if (p[a] < lo) m |= 1 << (a*2);
                if (p[a] > hi) m |= 1 << (a*2+1);
            }
            if (m) all_in = 0;
            out_mask &= m;
        }
        if (out_mask) continue;               /* fully outside one plane */
        if (all_in) {                          /* untouched: emit as-is */
            if (nf_out == fcap) {
                fcap = fcap * 2;
                int32_t *nf2 = (int32_t *)realloc(fout, fcap * 3 * sizeof(int32_t));
                if (!nf2) { rc = -1; break; }
                fout = nf2;
            }
            fout[nf_out*3+0] = t[0]; fout[nf_out*3+1] = t[1]; fout[nf_out*3+2] = t[2];
            nf_out++;
            continue;
        }

        /* Sutherland-Hodgman against the 6 half-spaces. */
        size_t pn = 3;
        for (int k = 0; k < 3; k++) {
            poly[k].p[0] = vg.p[(size_t)t[k]*3+0];
            poly[k].p[1] = vg.p[(size_t)t[k]*3+1];
            poly[k].p[2] = vg.p[(size_t)t[k]*3+2];
            poly[k].vid = t[k];
            poly[k].edge_u = poly[k].edge_v = t[k];   /* orig-vert marker */
        }
        for (int plane = 0; plane < 6 && pn >= 3; plane++) {
            int    ax   = plane / 2;
            int    is_hi = plane & 1;
            double pc   = is_hi ? hi : lo;
            size_t nn = 0;
            for (size_t i = 0; i < pn; i++) {
                const CutCorner *A = &poly[i];
                const CutCorner *B = &poly[(i + 1) % pn];
                double da = is_hi ? (pc - A->p[ax]) : (A->p[ax] - pc);
                double db = is_hi ? (pc - B->p[ax]) : (B->p[ax] - pc);
                int ina = (da >= 0.0), inb = (db >= 0.0);
                if (ina) next[nn++] = *A;
                if (ina == inb) continue;
                /* Crossing. Compute on the ORIGINAL edge when the segment has
                 * one (canonical direction, bit-identical across the two
                 * incident faces); else from the local segment. */
                int32_t su = -1, sv = -1;
                int shared = segment_support(A, B, &su, &sv);
                CutCorner c;
                c.edge_u = shared ? su : -1;
                c.edge_v = shared ? sv : -1;
                c.vid = -1;
                if (shared) {
                    const double *pu = &vg.p[(size_t)su*3];
                    const double *pv = &vg.p[(size_t)sv*3];
                    double denom = pv[ax] - pu[ax];
                    if (fabs(denom) < 1e-300) continue;   /* parallel: no cross */
                    double tt = (pc - pu[ax]) / denom;
                    if (tt < 0.0) tt = 0.0;
                    if (tt > 1.0) tt = 1.0;
                    for (int a = 0; a < 3; a++)
                        c.p[a] = pu[a] + tt * (pv[a] - pu[a]);
                } else {
                    double denom = B->p[ax] - A->p[ax];
                    if (fabs(denom) < 1e-300) continue;
                    double tt = (pc - A->p[ax]) / denom;
                    for (int a = 0; a < 3; a++)
                        c.p[a] = A->p[a] + tt * (B->p[a] - A->p[a]);
                }
                c.p[ax] = pc;                  /* exactly on the plane */
                /* Endpoint reuse: a crossing within snap_eps of a corner IS
                 * that corner (no sub-eps sliver strips). */
                {
                    double dA = 0.0, dB = 0.0;
                    for (int a = 0; a < 3; a++) {
                        double ea = c.p[a] - A->p[a]; dA += ea * ea;
                        double eb = c.p[a] - B->p[a]; dB += eb * eb;
                    }
                    if (dA <= eps * eps) { if (!ina) next[nn++] = *A; continue; }
                    if (dB <= eps * eps) continue;   /* B emitted next iter if inside */
                }
                if (shared) {
                    size_t slot = cuthash_find(&hash, su, sv, plane);
                    if (hash.ent[slot].vid >= 0) {
                        c.vid = hash.ent[slot].vid;
                        c.p[0] = vg.p[(size_t)c.vid*3+0];   /* exact shared coords */
                        c.p[1] = vg.p[(size_t)c.vid*3+1];
                        c.p[2] = vg.p[(size_t)c.vid*3+2];
                    } else {
                        if (vgrow_push(&vg, c.p) != 0) { rc = -1; break; }
                        c.vid = (int32_t)(vg.n - 1);
                        hash.ent[slot].vid = c.vid;
                    }
                }
                next[nn++] = c;
                if (nn >= 15) break;           /* paranoia: cannot happen (<=9-gon) */
            }
            if (rc != 0) break;
            memcpy(poly, next, nn * sizeof(CutCorner));
            pn = nn;
        }
        if (rc != 0) break;

        /* Drop consecutive duplicates (same emitted vid or coincident pos). */
        if (pn >= 3) {
            CutCorner ded[16];
            size_t dn = 0;
            for (size_t i = 0; i < pn; i++) {
                const CutCorner *A = &poly[i];
                const CutCorner *B = &ded[dn ? dn - 1 : 0];
                if (dn) {
                    if (A->vid >= 0 && B->vid >= 0 && A->vid == B->vid) continue;
                    double d2 = 0.0;
                    for (int a = 0; a < 3; a++) {
                        double e = A->p[a] - B->p[a]; d2 += e * e;
                    }
                    if (d2 <= 1e-24) continue;
                }
                ded[dn++] = *A;
            }
            /* wrap-around duplicate */
            if (dn >= 2) {
                const CutCorner *A = &ded[0], *B = &ded[dn-1];
                double d2 = 0.0;
                for (int a = 0; a < 3; a++) { double e = A->p[a] - B->p[a]; d2 += e * e; }
                if ((A->vid >= 0 && B->vid >= 0 && A->vid == B->vid) || d2 <= 1e-24)
                    dn--;
            }
            memcpy(poly, ded, dn * sizeof(CutCorner));
            pn = dn;
        }
        if (pn < 3) { faces_cut++; continue; }  /* clipped away (or degenerate) */

        /* Materialize private cut corners as verts, then fan-triangulate. */
        for (size_t i = 0; i < pn && rc == 0; i++) {
            if (poly[i].vid < 0) {
                if (vgrow_push(&vg, poly[i].p) != 0) { rc = -1; break; }
                poly[i].vid = (int32_t)(vg.n - 1);
            }
        }
        if (rc != 0) break;
        faces_cut++;
        for (size_t i = 1; i + 1 < pn; i++) {
            int32_t a = poly[0].vid, b = poly[i].vid, c = poly[i+1].vid;
            if (a == b || b == c || a == c) continue;
            /* degenerate-area guard (colinear fans from grazing cuts) */
            {
                double e1[3], e2[3], cr[3];
                for (int k = 0; k < 3; k++) {
                    e1[k] = vg.p[(size_t)b*3+k] - vg.p[(size_t)a*3+k];
                    e2[k] = vg.p[(size_t)c*3+k] - vg.p[(size_t)a*3+k];
                }
                cr[0] = e1[1]*e2[2] - e1[2]*e2[1];
                cr[1] = e1[2]*e2[0] - e1[0]*e2[2];
                cr[2] = e1[0]*e2[1] - e1[1]*e2[0];
                double a2 = cr[0]*cr[0] + cr[1]*cr[1] + cr[2]*cr[2];
                if (a2 <= 1e-24) continue;
            }
            if (nf_out == fcap) {
                fcap = fcap * 2;
                int32_t *nf2 = (int32_t *)realloc(fout, fcap * 3 * sizeof(int32_t));
                if (!nf2) { rc = -1; break; }
                fout = nf2;
            }
            fout[nf_out*3+0] = a; fout[nf_out*3+1] = b; fout[nf_out*3+2] = c;
            nf_out++;
        }
    }

    if (rc != 0 || nf_out == 0) {
        free(vg.p); free(hash.ent); free(fout);
        return rc;   /* rc==0 with empty output: everything clipped away */
    }

    /* Compact used verts into arena floats, remap faces (mirrors the
     * drop-only trim's passes 2-3). */
    uint8_t *used = (uint8_t *)calloc(vg.n, 1);
    if (!used) { free(vg.p); free(hash.ent); free(fout); return -1; }
    for (size_t f = 0; f < nf_out; f++) {
        used[fout[f*3+0]] = 1; used[fout[f*3+1]] = 1; used[fout[f*3+2]] = 1;
    }
    int32_t *remap = (int32_t *)malloc(vg.n * sizeof(int32_t));
    if (!remap) { free(vg.p); free(hash.ent); free(fout); free(used); return -1; }
    size_t nv_kept = 0;
    for (size_t v = 0; v < vg.n; v++)
        remap[v] = used[v] ? (int32_t)nv_kept++ : -1;

    float *new_verts = (float *)ARENA_ALLOC(arena,
                          (long)nv_kept * 3L * (long)sizeof(float));
    size_t wi = 0;
    for (size_t v = 0; v < vg.n; v++) {
        if (!used[v]) continue;
        new_verts[wi*3+0] = (float)vg.p[v*3+0];
        new_verts[wi*3+1] = (float)vg.p[v*3+1];
        new_verts[wi*3+2] = (float)vg.p[v*3+2];
        wi++;
    }
    int32_t *new_faces = (int32_t *)ARENA_ALLOC(arena,
                            (long)nf_out * 3L * (long)sizeof(int32_t));
    for (size_t f = 0; f < nf_out; f++) {
        new_faces[f*3+0] = remap[fout[f*3+0]];
        new_faces[f*3+1] = remap[fout[f*3+1]];
        new_faces[f*3+2] = remap[fout[f*3+2]];
    }

    free(vg.p); free(hash.ent); free(fout); free(used); free(remap);

    *out_verts = new_verts;
    *out_nv = nv_kept;
    *out_faces = new_faces;
    *out_nf = nf_out;
    if (out_faces_cut) *out_faces_cut = faces_cut;
    return 0;
}
