/*
 * seam_weld.c -- bridge cube-boundary seams of a concatenated multi-cube BPA
 * mesh by re-running Ball-Pivoting across the seam. See seam_weld.h.
 *
 * DIRECT BRIDGE (no peelback). Each cube is BPA'd on its own region (step0 trims
 * the cloud to the cube face), so adjacent cubes meet at the seam plane with a
 * ~1-vox gap. We collect the open boundary half-edges that lie IN a detected
 * seam plane (SeamPlanes_edge_in -- grid-edge perimeter boundaries are excluded)
 * and feed them as the BPA init front; the directed-front glue rolls a ball
 * across the gap at an escalating, capped radius, emitting bridge faces wound
 * consistently with the existing surface. The radius cap (2*rho_max below the
 * inter-wrap clearance) guarantees NO inter-wrap mergers. (An earlier version
 * peeled a ring of seam triangles first to widen the gap; that receded the
 * boundary and starved the bridge -- on real dense seams the direct bridge
 * closes far better, so the peel was removed.)
 *
 * Scratch memory is malloc/free (local to one call); only the returned
 * combined face array is arena-allocated.
 */
#include "seam_weld.h"
#include "../common/run_ctx.h"
#include "ball_pivot.h"
#include "seam_planes.h"
#include "../common/pipeline_constants.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846  /* MSVC math.h omits it without _USE_MATH_DEFINES */
#endif

/* SeamPlane + the plane detection/classification queries (detect, edge-in,
 * face-straddle, vert-dist) live in seam_planes.{h,c}, shared with the
 * weld-time band refine/recoarsen stages. Bodies moved there verbatim. */

static void compute_vertex_normals(const float *verts, size_t nv,
                                   const int32_t *faces, size_t nf,
                                   float *normals)
{
    memset(normals, 0, nv * 3 * sizeof(float));
    for (size_t f = 0; f < nf; f++) {
        int32_t i0 = faces[f*3+0], i1 = faces[f*3+1], i2 = faces[f*3+2];
        const float *p0 = &verts[(size_t)i0*3];
        const float *p1 = &verts[(size_t)i1*3];
        const float *p2 = &verts[(size_t)i2*3];
        double e1[3] = { (double)p1[0]-p0[0], (double)p1[1]-p0[1], (double)p1[2]-p0[2] };
        double e2[3] = { (double)p2[0]-p0[0], (double)p2[1]-p0[1], (double)p2[2]-p0[2] };
        /* area-weighted face normal (cross magnitude = 2*area) */
        double n0 = e1[1]*e2[2] - e1[2]*e2[1];
        double n1 = e1[2]*e2[0] - e1[0]*e2[2];
        double n2 = e1[0]*e2[1] - e1[1]*e2[0];
        normals[(size_t)i0*3+0] += (float)n0; normals[(size_t)i0*3+1] += (float)n1; normals[(size_t)i0*3+2] += (float)n2;
        normals[(size_t)i1*3+0] += (float)n0; normals[(size_t)i1*3+1] += (float)n1; normals[(size_t)i1*3+2] += (float)n2;
        normals[(size_t)i2*3+0] += (float)n0; normals[(size_t)i2*3+1] += (float)n1; normals[(size_t)i2*3+2] += (float)n2;
    }
    for (size_t v = 0; v < nv; v++) {
        double z = normals[v*3+0], y = normals[v*3+1], x = normals[v*3+2];
        double len = sqrt(z*z + y*y + x*x);
        if (len > 1e-20) {
            normals[v*3+0] = (float)(z/len);
            normals[v*3+1] = (float)(y/len);
            normals[v*3+2] = (float)(x/len);
        }
    }
}

typedef struct {
    int32_t key_u, key_v;      /* sorted endpoints (dedup key) */
    int32_t va, vb, opp;       /* directed half-edge + opposite vert */
} HalfEdge;

static int cmp_he(const void *pa, const void *pb)
{
    const HalfEdge *a = (const HalfEdge *)pa;
    const HalfEdge *b = (const HalfEdge *)pb;
    if (a->key_u != b->key_u) return (a->key_u < b->key_u) ? -1 : 1;
    if (a->key_v != b->key_v) return (a->key_v < b->key_v) ? -1 : 1;
    return 0;
}

static int cmp_double_asc(const void *pa, const void *pb)
{
    double a = *(const double *)pa, b = *(const double *)pb;
    return (a < b) ? -1 : (a > b ? 1 : 0);
}

static double edge_len2(const float *verts, int32_t a, int32_t b)
{
    double dz = (double)verts[(size_t)a*3+0] - (double)verts[(size_t)b*3+0];
    double dy = (double)verts[(size_t)a*3+1] - (double)verts[(size_t)b*3+1];
    double dx = (double)verts[(size_t)a*3+2] - (double)verts[(size_t)b*3+2];
    return dz*dz + dy*dy + dx*dx;
}

/* Open-addressing edge->face-count hash for the manifold-guarded merge of
 * bridge faces onto the existing mesh. The bridge cloud necessarily holds
 * existing-surface verts (init-edge endpoints + near-seam verts), so the
 * rolling ball can recreate an existing triangle or put a third face on an
 * interior edge; this hash lets us reject any bridge face that would push an
 * edge past two faces. */
typedef struct { int64_t key; int cnt; } ECell;

static size_t ec_pow2(size_t n)
{
    size_t p = 1;
    while (p < n) p <<= 1;
    return p;
}

static int64_t ec_key(int32_t a, int32_t b)
{
    int32_t lo = (a < b) ? a : b, hi = (a < b) ? b : a;
    return ((int64_t)lo << 32) | (int64_t)(uint32_t)hi;
}

static size_t ec_slot(const ECell *t, size_t mask, int64_t key, int *found)
{
    size_t i = (size_t)((uint64_t)key * 0x9E3779B97F4A7C15ULL) & mask;
    for (;;) {
        if (t[i].key == key) { *found = 1; return i; }
        if (t[i].key == -1)  { *found = 0; return i; }
        i = (i + 1) & mask;
    }
}

static int ec_count(const ECell *t, size_t mask, int64_t key)
{
    int f = 0; size_t i = ec_slot(t, mask, key, &f);
    return f ? t[i].cnt : 0;
}

static void ec_inc(ECell *t, size_t mask, int64_t key)
{
    int f = 0; size_t i = ec_slot(t, mask, key, &f);
    if (!f) { t[i].key = key; t[i].cnt = 0; }
    t[i].cnt++;
}

/* Debug: dump the BPA init front so a missed/under-detected seam is visible
 * BEFORE the bridge runs. <prefix>.front_tris.obj = the (va,vb,v_opp) incident
 * triangles BPA hinges off of (the literal start front). <prefix>.front_edges.obj
 * = EVERY run-1 boundary edge as a colored l-segment: green = selected (passed
 * SeamPlanes_edge_in, fed to BPA), red = excluded (run-1 boundary but NOT in a
 * seam plane). A seam edge that should bridge but shows red is a detection miss
 * (SeamPlanes_detect / _edge_in), not a BPA failure. Inline per-segment
 * verts (not shared) so a vertex on both a selected and an excluded edge can't
 * blend colors. Gated by SEAM_DUMP_FRONT; mirrors SEAM_DUMP_BRIDGE below. */
static void dump_seam_front(const char *prefix,
                            const float *verts, size_t nv,
                            const BpaInitEdge *init, size_t n_init,
                            const BpaInitEdge *excl, size_t n_excl)
{
    char path[1024];

    snprintf(path, sizeof(path), "%s.front_tris.obj", prefix);
    FILE *tf = fopen(path, "w");
    if (tf) {
        for (size_t v = 0; v < nv; v++)
            fprintf(tf, "v %.6f %.6f %.6f\n",
                    (double)verts[v*3+0], (double)verts[v*3+1],
                    (double)verts[v*3+2]);
        for (size_t i = 0; i < n_init; i++)
            fprintf(tf, "f %d %d %d\n",
                    init[i].va + 1, init[i].vb + 1, init[i].v_opp + 1);
        fclose(tf);
        fprintf(stderr, "  [SEAM_DUMP_FRONT] wrote %s (%zu front triangles)\n",
                path, n_init);
    }

    snprintf(path, sizeof(path), "%s.front_edges.obj", prefix);
    FILE *ef = fopen(path, "w");
    if (ef) {
        fprintf(ef, "# seam BPA init front: green = selected (bridged), "
                    "red = excluded (run-1 boundary edge not in a seam plane)\n");
        size_t vidx = 1;
        for (size_t i = 0; i < n_init; i++) {
            const float *a = &verts[(size_t)init[i].va * 3];
            const float *b = &verts[(size_t)init[i].vb * 3];
            fprintf(ef, "v %.6f %.6f %.6f 0.0 1.0 0.0\n",
                    (double)a[0], (double)a[1], (double)a[2]);
            fprintf(ef, "v %.6f %.6f %.6f 0.0 1.0 0.0\n",
                    (double)b[0], (double)b[1], (double)b[2]);
            fprintf(ef, "l %zu %zu\n", vidx, vidx + 1); vidx += 2;
        }
        for (size_t i = 0; i < n_excl; i++) {
            const float *a = &verts[(size_t)excl[i].va * 3];
            const float *b = &verts[(size_t)excl[i].vb * 3];
            fprintf(ef, "v %.6f %.6f %.6f 1.0 0.0 0.0\n",
                    (double)a[0], (double)a[1], (double)a[2]);
            fprintf(ef, "v %.6f %.6f %.6f 1.0 0.0 0.0\n",
                    (double)b[0], (double)b[1], (double)b[2]);
            fprintf(ef, "l %zu %zu\n", vidx, vidx + 1); vidx += 2;
        }
        fclose(ef);
        fprintf(stderr, "  [SEAM_DUMP_FRONT] wrote %s "
                "(%zu selected green, %zu excluded red)\n",
                path, n_init, n_excl);
    }
}

/* --- pre-bridge sliver cull ------------------------------------------------ */

/* Min altitude (2*area / longest edge) -- the "thinnest dimension" of a
 * triangle. A needle/sliver has a tiny value even at non-zero area. */
static double tri_min_alt3(const float *p0, const float *p1, const float *p2)
{
    double ax = p1[0]-p0[0], ay = p1[1]-p0[1], az = p1[2]-p0[2];
    double bx = p2[0]-p0[0], by = p2[1]-p0[1], bz = p2[2]-p0[2];
    double cz = ax*by - ay*bx, cy = az*bx - ax*bz, cx = ay*bz - az*by;
    double area2 = sqrt(cx*cx + cy*cy + cz*cz);            /* = 2*area */
    double l0 = sqrt(ax*ax+ay*ay+az*az);
    double l1 = sqrt(bx*bx+by*by+bz*bz);
    double dx = p2[0]-p1[0], dy = p2[1]-p1[1], dz = p2[2]-p1[2];
    double l2 = sqrt(dx*dx+dy*dy+dz*dz);
    double lmax = l0 > l1 ? (l0 > l2 ? l0 : l2) : (l1 > l2 ? l1 : l2);
    return lmax < 1e-12 ? 0.0 : area2 / lmax;
}

/* undirected edge with owning face + local edge index, for boundary detection */
typedef struct { int32_t u, v; int32_t face; int32_t le; } FaceEdge;
static int cmp_faceedge(const void *pa, const void *pb)
{
    const FaceEdge *a = (const FaceEdge *)pa, *b = (const FaceEdge *)pb;
    if (a->u != b->u) return a->u < b->u ? -1 : 1;
    if (a->v != b->v) return a->v < b->v ? -1 : 1;
    return 0;
}

/* Mark (in del[], length nf) two kinds of pre-bridge garbage near the seam:
 *
 *   (A) SLIVER: a triangle with a seam-plane boundary edge whose min altitude is
 *       below min_alt -- a needle the rolling ball would prime off of.
 *
 *   (B) SEAM-WARD DANGLING TIP: a triangle attached to the rest of the mesh by
 *       only ONE edge (its other two edges are open boundary), whose free vertex
 *       sits in the seam band AND is the closest of the three to the seam. This is
 *       exactly the leftover the bridge would run an edge straight across -- a
 *       T-junction. (A) tends to CREATE these: removing a thin boundary triangle
 *       orphans its shared vertex into such a tip, which is why the caller
 *       iterates -- a tip only appears once its sliver neighbours are gone.
 *
 * Returns the number of faces marked this pass. */
static size_t mark_sliver_and_tip_tris(const float *verts,
        const int32_t *faces, size_t nf,
        const SeamPlane *planes, size_t np, double band,
        double min_alt, uint8_t *del)
{
    size_t hn = nf * 3;
    FaceEdge *fe = (FaceEdge *)malloc((hn ? hn : 1) * sizeof(FaceEdge));
    for (size_t f = 0; f < nf; f++) {
        int32_t t[3] = { faces[f*3+0], faces[f*3+1], faces[f*3+2] };
        for (int e = 0; e < 3; e++) {
            int32_t a = t[e], b = t[(e+1)%3];
            FaceEdge *h = &fe[f*3+(size_t)e];
            h->u = (a < b) ? a : b; h->v = (a < b) ? b : a;
            h->face = (int32_t)f; h->le = (int32_t)e;
        }
    }
    qsort(fe, hn, sizeof(FaceEdge), cmp_faceedge);
    /* per-face 3-bit mask of which local edges are open boundary (run==1) */
    uint8_t *bmask = (uint8_t *)calloc(nf ? nf : 1, 1);
    for (size_t i = 0; i < hn; ) {
        size_t j = i + 1;
        while (j < hn && fe[j].u == fe[i].u && fe[j].v == fe[i].v) j++;
        if (j - i == 1) bmask[fe[i].face] |= (uint8_t)(1 << fe[i].le);
        i = j;
    }
    free(fe);

    size_t ndel = 0;
    for (size_t f = 0; f < nf; f++) {
        if (del[f] || bmask[f] == 0) continue;
        int32_t t[3] = { faces[f*3+0], faces[f*3+1], faces[f*3+2] };
        const float *p0 = &verts[(size_t)t[0]*3];
        const float *p1 = &verts[(size_t)t[1]*3];
        const float *p2 = &verts[(size_t)t[2]*3];
        /* (A) sliver with a seam-plane boundary edge */
        int has_seam_bnd = 0;
        for (int e = 0; e < 3; e++)
            if ((bmask[f] & (1<<e)) &&
                SeamPlanes_edge_in(verts, t[e], t[(e+1)%3], planes, np, band)) { has_seam_bnd = 1; break; }
        if (has_seam_bnd && tri_min_alt3(p0, p1, p2) < min_alt) { del[f] = 1; ndel++; continue; }
        /* (B) seam-ward dangling tip: >=2 boundary edges; the vertex they share
         *     pokes toward the seam (in band + closest of the three). */
        int nb = (bmask[f]&1) + ((bmask[f]>>1)&1) + ((bmask[f]>>2)&1);
        if (nb >= 2) {
            double d0 = SeamPlanes_vert_dist(verts, t[0], planes, np);
            double d1 = SeamPlanes_vert_dist(verts, t[1], planes, np);
            double d2 = SeamPlanes_vert_dist(verts, t[2], planes, np);
            int tip;
            if (nb == 3)                                  tip = (d0<=d1&&d0<=d2)?0:(d1<=d2?1:2);
            else if ((bmask[f]&1) && (bmask[f]&4))        tip = 0;   /* edges 0=(0,1),2=(2,0) -> v0 */
            else if ((bmask[f]&1) && (bmask[f]&2))        tip = 1;   /* edges 0,1 -> v1 */
            else                                          tip = 2;   /* edges 1,2 -> v2 */
            double dt = tip==0?d0:(tip==1?d1:d2);
            double dmin_other = 1e30;
            for (int k=0;k<3;k++) if (k!=tip){ double dk=k==0?d0:(k==1?d1:d2); if (dk<dmin_other) dmin_other=dk; }
            if (dt < band && dt < dmin_other) { del[f] = 1; ndel++; continue; }
        }
    }
    free(bmask);
    return ndel;
}

int SeamWeld_bridge(Arena_T arena,
                    const float *verts, size_t nv,
                    const int32_t *faces, size_t nf,
                    float cube_size, float rho, float rho_max_in, float band,
                    const uint8_t *want_mask,
                    const BpaBridgeGate *gate,
                    int32_t **out_faces, size_t *out_nf,
                    size_t *out_n_bridge)
{
    *out_faces = NULL;
    *out_nf = 0;
    if (out_n_bridge) *out_n_bridge = 0;
    if (nf == 0) return 0;

    /* 1) Per-vertex normals from faces (QEM mesh carries none). */
    float *normals = (float *)malloc(nv * 3 * sizeof(float));
    compute_vertex_normals(verts, nv, faces, nf, normals);

    /* vert-used (by any face) bitmap, for plane detection + cloud build. */
    uint8_t *used_any = (uint8_t *)calloc(nv, 1);
    for (size_t f = 0; f < nf; f++) {
        used_any[faces[f*3+0]] = 1;
        used_any[faces[f*3+1]] = 1;
        used_any[faces[f*3+2]] = 1;
    }

    /* 2) Detect seam planes. */
    SeamPlane planes[64];
    size_t np = SeamPlanes_detect(verts, nv, used_any, (double)cube_size,
                                  (double)band, planes, 64);
    if (np == 0) {
        /* Nothing to bridge -- return faces unchanged. */
        int32_t *out = (int32_t *)ARENA_ALLOC(arena, (size_t)(nf*3*sizeof(int32_t)));
        memcpy(out, faces, nf*3*sizeof(int32_t));
        *out_faces = out; *out_nf = nf;
        free(normals); free(used_any);
        return 0;
    }

    /* rho_max bounds the escalating bridge radius AND the over-long bridge-face
     * cutoff (2*rho_max < ~CUT_GAP_DEPTH => NO inter-wrap mergers). Caller arg
     * (a phase-2 restricted re-weld passes a WIDER cap to span divots); <=0
     * falls back to BRIDGE_RHO_MAX; SEAM_RHO_MAX env still overrides. */
    double rho_max = (rho_max_in > 0.0f) ? (double)rho_max_in : (double)BRIDGE_RHO_MAX;
    { const char *e = sf_env("SEAM_RHO_MAX"); if (e) rho_max = atof(e); }

    /* 2.5) Pre-bridge sliver + tip cull: drop sliver boundary triangles AND the
     *      seam-ward dangling tips their removal exposes, so the bridge front
     *      never primes off a needle NOR runs an edge across a leftover tip (a
     *      T-junction / doubled strip). ITERATED to a fixpoint: removing a sliver
     *      orphans its shared vertex into a tip, which the next pass removes.
     *      Rebind faces/nf to the cleaned list so the front build, bridge, and
     *      output all operate on (and emit) the cleaned mesh. SEAM_SLIVER_MIN_ALT=0
     *      disables. */
    int32_t *sliver_owned = NULL;
    {
        double sliver_alt = (double)SEAM_SLIVER_MIN_ALT;
        const char *se = sf_env("SEAM_SLIVER_MIN_ALT"); if (se) sliver_alt = atof(se);
        if (sliver_alt > 0.0) {
            size_t nf0 = nf, total_del = 0;
            for (int pass = 0; pass < 8; pass++) {
                uint8_t *del = (uint8_t *)calloc(nf ? nf : 1, 1);
                size_t ndel = mark_sliver_and_tip_tris(verts, faces, nf, planes, np,
                                                       (double)band, sliver_alt, del);
                if (ndel == 0) { free(del); break; }
                size_t nfc = nf - ndel;
                int32_t *fc = (int32_t *)malloc((nfc ? nfc : 1) * 3 * sizeof(int32_t));
                size_t k = 0;
                for (size_t f = 0; f < nf; f++) {
                    if (del[f]) continue;
                    fc[k*3+0] = faces[f*3+0]; fc[k*3+1] = faces[f*3+1]; fc[k*3+2] = faces[f*3+2];
                    k++;
                }
                free(del);
                free(sliver_owned);          /* free the previous pass's array (NULL ok) */
                faces = fc; nf = nfc; sliver_owned = fc;
                total_del += ndel;
            }
            if (total_del > 0)
                compute_vertex_normals(verts, nv, faces, nf, normals);  /* refresh */
            fprintf(stderr, "  [seam] sliver+tip cull: dropped %zu of %zu faces "
                    "(min_alt < %.2f vox, iterated)\n", total_del, nf0, sliver_alt);
        }
    }

    /* 3) Collect boundary half-edges (undirected run==1) over the faces whose
     *    edge sits near a seam plane -> the BPA init front. (va,vb) is in the
     *    owning face's winding order, so the directed-front glue emits bridge
     *    faces wound consistently with the existing surface. */
    size_t hn = nf * 3;
    HalfEdge *he = (HalfEdge *)malloc((hn ? hn : 1) * sizeof(HalfEdge));
    for (size_t f = 0; f < nf; f++) {
        int32_t t[3] = { faces[f*3+0], faces[f*3+1], faces[f*3+2] };
        for (int e = 0; e < 3; e++) {
            int32_t a = t[e], b = t[(e+1)%3], o = t[(e+2)%3];
            HalfEdge *h = &he[f*3+(size_t)e];
            h->va = a; h->vb = b; h->opp = o;
            h->key_u = (a < b) ? a : b;
            h->key_v = (a < b) ? b : a;
        }
    }
    qsort(he, hn, sizeof(HalfEdge), cmp_he);

    /* Boundary edges: undirected run length 1, both endpoints near a plane. */
    BpaInitEdge *init = (BpaInitEdge *)malloc((hn ? hn : 1) * sizeof(BpaInitEdge));
    size_t n_init = 0;
    /* Run-1 boundary edges that FAIL the seam-plane test -- kept only for the
     * SEAM_DUMP_FRONT diagnostic, so a missed seam edge shows up red. */
    BpaInitEdge *excl = (BpaInitEdge *)malloc((hn ? hn : 1) * sizeof(BpaInitEdge));
    size_t n_excl = 0;
    for (size_t i = 0; i < hn; ) {
        size_t j = i + 1;
        while (j < hn && he[j].key_u == he[i].key_u && he[j].key_v == he[i].key_v) j++;
        if (j - i == 1) {
            int32_t va = he[i].va, vb = he[i].vb, opp = he[i].opp;
            if (SeamPlanes_edge_in(verts, va, vb, planes, np, (double)band)) {
                init[n_init].va = va; init[n_init].vb = vb; init[n_init].v_opp = opp;
                n_init++;
            } else {
                excl[n_excl].va = va; excl[n_excl].vb = vb; excl[n_excl].v_opp = opp;
                n_excl++;
            }
        }
        i = j;
    }

    /* 3b) Grazing-seam promotion. SeamPlanes_edge_in's parallel test rejects
     * boundary edges that RUN ACROSS the plane -- right for grid-perimeter
     * edges, wrong at a GRAZING seam (the wrap running parallel to the cube
     * face, e.g. the umbilicus-aligned band of every seam plane): there the
     * per-cube centroid trim interleaves the two charts' coverage, the open
     * boundaries wander obliquely off-plane, and since the bridge cloud is
     * init-edge verts ONLY, an excluded edge's verts are invisible to the
     * ball at ANY radius -- the gap can never close. Promote an excluded
     * edge when some boundary edge on the OPPOSITE side of the plane lies
     * within bridge reach (2*rho_max) of its midpoint: a real grazing seam
     * has an opposing open boundary; a perimeter edge has nothing across
     * from it.
     *
     * SAME-WRAP guard: at the compacted core, DIFFERENT wraps' boundaries can
     * also face each other within reach across a plane -- promoting those
     * re-admits cross-wrap mergers (4x5x5 A/B: handles 7 -> 14). So a pair
     * only counts as opposing if it is same-wrap by the winding phase
     * w = r/pitch - theta/2pi about the umbilicus (|dw| <= tol), taken from the
     * same BpaBridgeGate the bridge uses. Without that scroll knowledge the
     * promotion CANNOT tell wraps apart and stays OFF -- which is also the case
     * for a phase-2 pair re-weld (gate == NULL): there the cloud is already
     * restricted to two confirmed sheets, so promotion is neither needed nor
     * safe. SEAM_NO_GRAZING_PROMOTE=1 also disables. */
    uint8_t *promoted_vert = NULL;  /* marks verts of promoted (grazing) edges;
                                     * their bridge faces bypass the straddle
                                     * test in the merge filter below. */
    double grz_umb_y = 0.0, grz_umb_x = 0.0, grz_pitch = 0.0, grz_tol = 0.0;
    {
        if (gate && gate->pitch > 0.0 && (gate->umb_y != 0.0 || gate->umb_x != 0.0)) {
            grz_umb_y = gate->umb_y; grz_umb_x = gate->umb_x; grz_pitch = gate->pitch;
            grz_tol = (gate->tol > 0.0) ? gate->tol : 0.25;
        }
    }
    if (n_excl > 0 && grz_tol > 0.0 && !sf_env("SEAM_NO_GRAZING_PROMOTE")) {
        double reach = 2.0 * rho_max;
        size_t n_all = n_init + n_excl;
        /* midpoints for every run-1 edge (init first, then excl). */
        float *mid = (float *)malloc(n_all * 3 * sizeof(float));
        for (size_t i = 0; i < n_all; i++) {
            const BpaInitEdge *e = (i < n_init) ? &init[i] : &excl[i - n_init];
            mid[i*3+0] = 0.5f*(verts[(size_t)e->va*3+0] + verts[(size_t)e->vb*3+0]);
            mid[i*3+1] = 0.5f*(verts[(size_t)e->va*3+1] + verts[(size_t)e->vb*3+1]);
            mid[i*3+2] = 0.5f*(verts[(size_t)e->va*3+2] + verts[(size_t)e->vb*3+2]);
        }
        /* uniform grid hash over midpoints, cell = reach, chained buckets. */
        size_t hcap = ec_pow2(n_all * 2 + 16);
        int64_t *hkey = (int64_t *)malloc(hcap * sizeof(int64_t));
        int32_t *hidx = (int32_t *)malloc(hcap * sizeof(int32_t));
        int32_t *hchain = (int32_t *)malloc(n_all * sizeof(int32_t));
        for (size_t i = 0; i < hcap; i++) { hkey[i] = -1; hidx[i] = -1; }
        size_t hmask = hcap - 1;
        double inv_cell = 1.0 / reach;
        #define GRZ_CELL_KEY(cz, cy, cx) \
            ( (((int64_t)(cz) & 0x1FFFFF) << 42) | (((int64_t)(cy) & 0x1FFFFF) << 21) \
              | ((int64_t)(cx) & 0x1FFFFF) )
        for (size_t i = 0; i < n_all; i++) {
            int64_t cz = (int64_t)floor((double)mid[i*3+0] * inv_cell);
            int64_t cy = (int64_t)floor((double)mid[i*3+1] * inv_cell);
            int64_t cx = (int64_t)floor((double)mid[i*3+2] * inv_cell);
            int64_t key = GRZ_CELL_KEY(cz, cy, cx);
            size_t s = (size_t)key & hmask;
            while (hkey[s] != -1 && hkey[s] != key) s = (s + 1) & hmask;
            hkey[s] = key;
            hchain[i] = hidx[s];
            hidx[s] = (int32_t)i;
        }
        size_t n_promoted = 0;
        for (size_t p = 0; p < np; p++) {
            int ax = planes[p].axis;
            double co = planes[p].coord;
            for (size_t x = 0; x < n_excl; x++) {
                if (excl[x].va < 0) continue;                 /* already promoted */
                size_t i = n_init + x;   /* mid[] index (excl block is stable) */
                double mc = (double)mid[i*3+(size_t)ax];
                if (fabs(mc - co) >= (double)band) continue;  /* not this seam */
                int side = (mc > co);
                int found = 0;
                int64_t bz = (int64_t)floor((double)mid[i*3+0] * inv_cell);
                int64_t by = (int64_t)floor((double)mid[i*3+1] * inv_cell);
                int64_t bx = (int64_t)floor((double)mid[i*3+2] * inv_cell);
                for (int64_t dz = -1; dz <= 1 && !found; dz++)
                for (int64_t dy = -1; dy <= 1 && !found; dy++)
                for (int64_t dx = -1; dx <= 1 && !found; dx++) {
                    int64_t key = GRZ_CELL_KEY(bz+dz, by+dy, bx+dx);
                    size_t s = (size_t)key & hmask;
                    while (hkey[s] != -1 && hkey[s] != key) s = (s + 1) & hmask;
                    if (hkey[s] != key) continue;
                    for (int32_t k = hidx[s]; k != -1; k = hchain[k]) {
                        if ((size_t)k == i) continue;
                        double kc = (double)mid[(size_t)k*3+(size_t)ax];
                        if (fabs(kc - co) >= (double)band) continue;
                        if ((kc > co) == side) continue;      /* same side */
                        double ddz = (double)mid[(size_t)k*3+0] - (double)mid[i*3+0];
                        double ddy = (double)mid[(size_t)k*3+1] - (double)mid[i*3+1];
                        double ddx = (double)mid[(size_t)k*3+2] - (double)mid[i*3+2];
                        if (ddz*ddz + ddy*ddy + ddx*ddx > reach*reach) continue;
                        /* same-wrap check: winding phase about the umbilicus */
                        {
                            double iy = (double)mid[i*3+1] - grz_umb_y;
                            double ix = (double)mid[i*3+2] - grz_umb_x;
                            double ky = (double)mid[(size_t)k*3+1] - grz_umb_y;
                            double kx = (double)mid[(size_t)k*3+2] - grz_umb_x;
                            double dr  = hypot(ky, kx) - hypot(iy, ix);
                            double dth = atan2(ky, kx) - atan2(iy, ix);
                            while (dth >  M_PI) dth -= 2.0*M_PI;
                            while (dth < -M_PI) dth += 2.0*M_PI;
                            double dw = dr/grz_pitch - dth/(2.0*M_PI);
                            if (fabs(dw) > grz_tol) continue; /* different wrap */
                        }
                        found = 1; break;
                    }
                }
                if (found) {
                    /* append to init (hn-sized: n_init+n_excl <= hn). Do NOT
                     * touch mid[]'s excl block -- indices stay stable. */
                    if (!promoted_vert)
                        promoted_vert = (uint8_t *)calloc(nv, 1);
                    promoted_vert[excl[x].va] = 1;
                    promoted_vert[excl[x].vb] = 1;
                    init[n_init + n_promoted] = excl[x];
                    n_promoted++;
                    excl[x].va = -1;                          /* mark promoted */
                }
            }
        }
        #undef GRZ_CELL_KEY
        if (n_promoted > 0) {
            n_init += n_promoted;
            size_t w = 0;                     /* compact excl[] for the dump */
            for (size_t x = 0; x < n_excl; x++)
                if (excl[x].va >= 0) excl[w++] = excl[x];
            n_excl = w;
            fprintf(stderr, "  [seam] grazing promote: %zu excluded edge(s) "
                    "with an opposing boundary within %.1f vox joined the front\n",
                    n_promoted, reach);
        }
        free(hchain); free(hidx); free(hkey); free(mid);
    }

    /* Phase-2 restriction: keep only front edges wholly inside want_mask (the
     * two confirmed sheets of a pair re-weld). The cloud (built from init verts
     * below) is thereby restricted to those two sheets, so a gate-off permissive
     * weld can only join them. NULL want_mask = the primary weld (all edges). */
    if (want_mask) {
        size_t w = 0;
        for (size_t i = 0; i < n_init; i++) {
            int32_t a = init[i].va, b = init[i].vb, o = init[i].v_opp;
            if (want_mask[a] && want_mask[b] && want_mask[o]) init[w++] = init[i];
        }
        n_init = w;
    }

    /* Emit the init-front diagnostic BEFORE any early-out, so an empty or
     * under-detected front (n_init == 0) is itself visible in the dump. */
    {
        const char *fpfx = sf_env("SEAM_DUMP_FRONT");
        if (fpfx) dump_seam_front(fpfx, verts, nv, init, n_init, excl, n_excl);
    }

    if (n_init == 0) {
        int32_t *out = (int32_t *)ARENA_ALLOC(arena, (size_t)((nf ? nf : 1)*3*sizeof(int32_t)));
        memcpy(out, faces, nf*3*sizeof(int32_t));
        *out_faces = out; *out_nf = nf;
        free(normals); free(used_any); free(he); free(init); free(excl);
        free(sliver_owned); free(promoted_vert);
        return 0;
    }

    /* 4) Bridge cloud = the boundary-edge verts (va, vb, opp) only. It
     *    deliberately excludes the deeper interior: re-admitting already-
     *    triangulated interior verts would only let the rolling ball re-
     *    triangulate existing surface (rejected by the merge guard -- wasted
     *    work) or over-reach. opp supplies the hinge + same-side context so the
     *    ball rolls OUTWARD into the gap, not back over the surface. Vertex-
     *    index order is preserved for determinism. */
    int32_t *g2l = (int32_t *)malloc(nv * sizeof(int32_t));
    for (size_t v = 0; v < nv; v++) g2l[v] = -1;
    uint8_t *want = (uint8_t *)calloc(nv, 1);
    for (size_t i = 0; i < n_init; i++) {
        want[init[i].va] = 1; want[init[i].vb] = 1; want[init[i].v_opp] = 1;
    }
    size_t ln = 0;
    for (size_t v = 0; v < nv; v++) if (want[v]) ln++;

    int32_t *l2g = (int32_t *)malloc((ln ? ln : 1) * sizeof(int32_t));
    float *lv = (float *)malloc((ln ? ln : 1) * 3 * sizeof(float));
    float *lnrm = (float *)malloc((ln ? ln : 1) * 3 * sizeof(float));
    size_t li = 0;
    for (size_t v = 0; v < nv; v++) {
        if (!want[v]) continue;
        g2l[v] = (int32_t)li;
        l2g[li] = (int32_t)v;
        lv[li*3+0] = verts[v*3+0]; lv[li*3+1] = verts[v*3+1]; lv[li*3+2] = verts[v*3+2];
        lnrm[li*3+0] = normals[v*3+0]; lnrm[li*3+1] = normals[v*3+1]; lnrm[li*3+2] = normals[v*3+2];
        li++;
    }
    /* Remap init edges to local indices. */
    BpaInitEdge *init_local = (BpaInitEdge *)malloc(n_init * sizeof(BpaInitEdge));
    for (size_t i = 0; i < n_init; i++) {
        init_local[i].va = g2l[init[i].va];
        init_local[i].vb = g2l[init[i].vb];
        init_local[i].v_opp = g2l[init[i].v_opp];
    }

    /* 5) Adaptive + escalating + capped ball radius. base = K * median near-
     *    seam boundary-edge length, clamped to [MIN, MAX], with the caller's
     *    rho (SEAM_RHO) as a floor; the front escalates up to rho_max (= MAX,
     *    env SEAM_RHO_MAX) so post-QEM sparse boundaries still close while
     *    2*rho_max stays below the inter-wrap clearance (NO inter-wrap merge). */
    double *elen = (double *)malloc(n_init * sizeof(double));
    for (size_t i = 0; i < n_init; i++)
        elen[i] = sqrt(edge_len2(verts, init[i].va, init[i].vb));
    qsort(elen, n_init, sizeof(double), cmp_double_asc);
    double s_med = elen[n_init/2];
    free(elen);

    double base = (double)BRIDGE_RHO_K * s_med;
    if (base < (double)BRIDGE_RHO_MIN) base = (double)BRIDGE_RHO_MIN;
    if (base > (double)BRIDGE_RHO_MAX) base = (double)BRIDGE_RHO_MAX;
    if ((double)rho > base) base = (double)rho;          /* caller/env floor */
    if (rho_max < base) rho_max = base;   /* rho_max read above for the peel cap */

    fprintf(stderr, "  [seam] planes=%zu init=%zu cloud=%zu verts; boundary-edge "
            "median=%.2f -> rho=%.2f..%.2f\n",
            np, n_init, ln, s_med, base, rho_max);

    /* Coarse-rim detector: the one-line tripwire for THIS bridge's silent
     * failure class. A front triangle primes only when its circumradius fits
     * the ball (<= rho_max); a rim spaced much beyond ~1.2*rho_max cannot
     * prime anywhere and the seam stays open with a clean-looking audit
     * (measured: 11.7-vox rims -> 2 bridge faces, 23.5k unpaired). If this
     * fires, an upstream densifier is missing (grid_weld's seam refine off?)
     * or a producer shipped a coarse rim. */
    if (s_med > 1.5 * rho_max) {
        fprintf(stderr,
            "  [seam] WARNING: boundary-edge median %.2f vox >> bridge reach "
            "(rho_max %.2f): front triangles cannot prime, expect a near-zero "
            "bridge. Refine the seam band before bridging (SEAM_NO_REFINE "
            "unset?) or ship finer rims.\n", s_med, rho_max);
    }

    /* 6) Bridge (single call; radius escalation happens inside BallPivot_bridge,
     *    where all radii share one edge store so the glue stays manifold). */
    int32_t *bridge_local = NULL;
    size_t n_bridge = 0;
    BallPivot_bridge(arena, lv, lnrm, ln, (float)base, (float)rho_max,
                     init_local, n_init, gate, &bridge_local, &n_bridge);
    fprintf(stderr, "  [seam] BPA produced %zu bridge faces\n", n_bridge);

    /* 7) Combine all original faces + bridge faces, manifold-guarded. Drop a
     *    bridge face that is degenerate, over-long (longest edge > 2*rho_max =>
     *    possible inter-wrap span), or would push any edge past two faces. */
    double maxedge2 = (2.0 * rho_max) * (2.0 * rho_max);
    size_t ecap = ec_pow2((nf + n_bridge) * 3 * 2 + 16);
    ECell *ec = (ECell *)malloc(ecap * sizeof(ECell));
    for (size_t i = 0; i < ecap; i++) { ec[i].key = -1; ec[i].cnt = 0; }
    size_t emask = ecap - 1;
    for (size_t f = 0; f < nf; f++) {   /* seed the edge-count hash from the input
                                         * faces so the 3-fan guard rejects any
                                         * bridge face that over-fills an edge. */
        ec_inc(ec, emask, ec_key(faces[f*3+0], faces[f*3+1]));
        ec_inc(ec, emask, ec_key(faces[f*3+1], faces[f*3+2]));
        ec_inc(ec, emask, ec_key(faces[f*3+2], faces[f*3+0]));
    }

    size_t total = nf + n_bridge;
    int32_t *out = (int32_t *)ARENA_ALLOC(arena, (size_t)((total ? total : 1) * 3 * sizeof(int32_t)));
    memcpy(out, faces, nf * 3 * sizeof(int32_t));
    size_t outn = nf, accepted = 0, rej_long = 0, rej_fold = 0;
    for (size_t f = 0; f < n_bridge; f++) {
        int32_t g0 = l2g[bridge_local[f*3+0]];
        int32_t g1 = l2g[bridge_local[f*3+1]];
        int32_t g2 = l2g[bridge_local[f*3+2]];
        if (g0 == g1 || g1 == g2 || g0 == g2) continue;        /* degenerate */
        if (SeamPlanes_face_straddle(verts, g0, g1, g2, planes, np, (double)band) < 0) {
            /* The straddle test assumes the sheet CROSSES the plane; at a
             * grazing seam legit closure faces can sit entirely on one side.
             * Bypassing it for faces touching a promoted vert was tried and
             * REJECTED as a default: on the 4x5x5 it admitted fold-backs
             * (same_dir 0 -> 34, handles 7 -> 14) for only ~50 fewer open
             * edges. SEAM_GRAZING_BYPASS=1 re-enables for experiments. */
            static int grz_bypass = -1;
            if (grz_bypass < 0)
                grz_bypass = sf_env("SEAM_GRAZING_BYPASS") ? 1 : 0;
            int grazing = grz_bypass && promoted_vert &&
                          (promoted_vert[g0] || promoted_vert[g1] ||
                           promoted_vert[g2]);
            if (!grazing) {
                rej_fold++; continue;  /* all one side -> fold-back, not a bridge */
            }
        }
        double em = edge_len2(verts, g0, g1);
        double e1 = edge_len2(verts, g1, g2); if (e1 > em) em = e1;
        double e2 = edge_len2(verts, g2, g0); if (e2 > em) em = e2;
        if (em > maxedge2) { rej_long++; continue; }           /* inter-wrap guard */
        int64_t k0 = ec_key(g0, g1), k1 = ec_key(g1, g2), k2 = ec_key(g2, g0);
        if (ec_count(ec, emask, k0) >= 2 ||
            ec_count(ec, emask, k1) >= 2 ||
            ec_count(ec, emask, k2) >= 2) continue;            /* would 3-fan */
        ec_inc(ec, emask, k0); ec_inc(ec, emask, k1); ec_inc(ec, emask, k2);
        out[outn*3+0] = g0; out[outn*3+1] = g1; out[outn*3+2] = g2;
        outn++; accepted++;
    }
    free(ec);

    *out_faces = out;
    *out_nf = outn;
    if (out_n_bridge) *out_n_bridge = accepted;
    fprintf(stderr, "  [seam] merged: %zu bridge faces kept (%zu over-long, "
            "%zu fold-back rejected), %zu -> %zu faces\n",
            accepted, rej_long, rej_fold, nf, outn);

    /* Debug: dump the accepted bridge faces (global verts) in isolation. */
    {
        const char *dp = sf_env("SEAM_DUMP_BRIDGE");
        if (dp) {
            FILE *bf = fopen(dp, "w");
            if (bf) {
                for (size_t v = 0; v < nv; v++)
                    fprintf(bf, "v %.6f %.6f %.6f\n",
                            (double)verts[v*3+0], (double)verts[v*3+1],
                            (double)verts[v*3+2]);
                for (size_t f = nf; f < outn; f++)
                    fprintf(bf, "f %d %d %d\n",
                            out[f*3+0]+1, out[f*3+1]+1, out[f*3+2]+1);
                fclose(bf);
                fprintf(stderr, "  [SEAM_DUMP_BRIDGE] wrote %s (%zu bridge faces)\n",
                        dp, outn - nf);
            }
        }
    }

    free(normals); free(used_any); free(he); free(init); free(excl); free(promoted_vert);
    free(g2l); free(want); free(l2g); free(lv); free(lnrm); free(init_local);
    free(sliver_owned);
    return 0;
}
