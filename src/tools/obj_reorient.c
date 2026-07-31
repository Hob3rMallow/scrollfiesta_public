/* ============================================================================
 * obj_reorient -- fix cross-component winding flips and non-orientable knots on
 * a welded OBJ surface mesh (e.g. the whole-grid LOD monolith).
 *
 * Two defects show up on a welded scroll mesh even when it is 2-manifold:
 *   1. Whole connected COMPONENTS wound backward relative to the bulk -- each
 *      internally consistent but globally flipped (they render the wrong colour
 *      in a two-sided shader). The post-weld OrientWeld pass misses spatially
 *      isolated ones when its radial reference is noisy.
 *   2. A residual `same_dir` count > 0: a genuine non-orientable knot (an odd
 *      cycle the winding BFS cannot satisfy). Flipping cannot fix it; the knot
 *      must be CUT.
 *
 * Pipeline (matches the weld's orient tail, then a knot cut):
 *   OrientMesh_consistent  -- make each component internally consistent
 *   OrientWeld_components_axis (FIXED radial reference) -- flip backward comps
 *   find_and_cut_knots      -- delete the frustrated face of each same_dir edge
 *   convention guard        -- keep the mesh's global in/out sense unchanged
 *
 * NOTE: OrientMesh is run BEFORE OrientWeld only. Running it AFTER would re-pick
 * an arbitrary per-component global sign and undo OrientWeld's radial fix.
 *
 * Axis convention: verts and axis_point are (z,y,x); default axis is the scroll
 * umbilicus (0, 3405, 2878) along +Z -- the same frame obj_orient_audit uses.
 * ==========================================================================*/

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/arena.h"
#include "../common/obj_io.h"
#include "../common/mesh_manifold.h"
#include "../common/mesh_types.h"
#include "../remesh/orient_mesh.h"
#include "../remesh/orient_weld.h"
#include "../remesh/pinhole_fill.h"     /* PinholeFill_split_pinches (manifold cut) */

/* Area-weighted radial-clarity gate, identical to obj_orient_audit's default. */
#define OBJ_REORIENT_MIN_COS 0.3

typedef struct {
    size_t nf_in, nf_out;
    size_t om_flips, om_components;         /* OrientMesh */
    size_t ow_flipped, ow_radial;           /* OrientWeld (spatial-aware) */
    size_t audit_flipped;                    /* audit-criterion component flips */
    size_t region_rounds;                    /* cut-and-flip repair rounds run */
    size_t region_cuts;                      /* interface edges seam-cut (0 faces deleted) */
    size_t region_flipped;                   /* components flipped by region fix */
    size_t region_dups;                      /* verts duplicated (seam + pinch) */
    size_t knots_cut;                        /* faces deleted to break knots */
    size_t pinch_splits;                     /* pinch verts split to stay manifold */
    size_t same_dir_before, same_dir_after;
    size_t nm_before, nm_after;              /* non-manifold EDGES */
    size_t nmv_before, nmv_after;            /* non-manifold (pinch) VERTS */
    int    global_flip;                      /* whole mesh reversed to keep sense */
} ReorientStats;

/* ------------------------------------------------------------- geometry */

static double face_area(const float *verts, const int32_t *faces, size_t f)
{
    const float *p0 = &verts[(size_t)faces[f * 3 + 0] * 3];
    const float *p1 = &verts[(size_t)faces[f * 3 + 1] * 3];
    const float *p2 = &verts[(size_t)faces[f * 3 + 2] * 3];
    double e1[3], e2[3], gn[3];
    int k = 0;
    for (k = 0; k < 3; k++) { e1[k] = (double)p1[k] - p0[k]; e2[k] = (double)p2[k] - p0[k]; }
    gn[0] = e1[1] * e2[2] - e1[2] * e2[1];
    gn[1] = e1[2] * e2[0] - e1[0] * e2[2];
    gn[2] = e1[0] * e2[1] - e1[1] * e2[0];
    return 0.5 * sqrt(gn[0] * gn[0] + gn[1] * gn[1] + gn[2] * gn[2]);
}

/* Per-face area-weighted radial vote: +1 = normal points AWAY from the axis
 * (out), -1 = toward (in), 0 = neutral (degenerate, on-axis, or |cos| below
 * the clarity gate). *out_area is always the face area. Identical measure to
 * obj_orient_audit's. */
static int face_radial(const float *verts, const int32_t *faces, size_t f,
                       const float axp[3], const double ad[3], double *out_area)
{
    size_t a = (size_t)faces[f * 3 + 0], b = (size_t)faces[f * 3 + 1],
           c = (size_t)faces[f * 3 + 2];
    double e1[3], e2[3], gn[3], cen[3], rv[3], area, dax, rlen, cosr;
    int k = 0;
    for (k = 0; k < 3; k++) {
        e1[k] = (double)verts[b * 3 + k] - verts[a * 3 + k];
        e2[k] = (double)verts[c * 3 + k] - verts[a * 3 + k];
        cen[k] = ((double)verts[a * 3 + k] + verts[b * 3 + k] + verts[c * 3 + k]) / 3.0
                 - (double)axp[k];
    }
    gn[0] = e1[1] * e2[2] - e1[2] * e2[1];
    gn[1] = e1[2] * e2[0] - e1[0] * e2[2];
    gn[2] = e1[0] * e2[1] - e1[1] * e2[0];
    area = 0.5 * sqrt(gn[0] * gn[0] + gn[1] * gn[1] + gn[2] * gn[2]);
    if (out_area) *out_area = area;
    if (area < 1e-12) return 0;
    dax = cen[0] * ad[0] + cen[1] * ad[1] + cen[2] * ad[2];
    for (k = 0; k < 3; k++) rv[k] = cen[k] - dax * ad[k];
    rlen = sqrt(rv[0] * rv[0] + rv[1] * rv[1] + rv[2] * rv[2]);
    if (rlen < 1.0) return 0;   /* on-axis: radial dir undefined */
    cosr = (gn[0] * rv[0] + gn[1] * rv[1] + gn[2] * rv[2]) / (2.0 * area * rlen);
    if (cosr >= OBJ_REORIENT_MIN_COS)  return 1;
    if (cosr <= -OBJ_REORIENT_MIN_COS) return -1;
    return 0;
}

/* Whole-mesh area-weighted, min-cos-gated radial sign: +1 = normals mostly
 * point AWAY from the axis (out), -1 = toward (in), 0 = undecidable. */
static int mesh_radial_sign(const float *verts, const int32_t *faces, size_t nf,
                            const float axp[3], const double ad[3])
{
    double out = 0.0, in = 0.0, area = 0.0;
    size_t f = 0;
    for (f = 0; f < nf; f++) {
        int s = face_radial(verts, faces, f, axp, ad, &area);
        if (s > 0) out += area;
        else if (s < 0) in += area;
    }
    if (out == 0.0 && in == 0.0) return 0;
    return (out >= in) ? 1 : -1;
}

static void flip_all(int32_t *faces, size_t nf)
{
    size_t f = 0;
    for (f = 0; f < nf; f++) {
        int32_t t = faces[f * 3 + 1];
        faces[f * 3 + 1] = faces[f * 3 + 2];
        faces[f * 3 + 2] = t;
    }
}

static int32_t uf_find_local(int32_t *uf, int32_t x)
{
    while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; }
    return x;
}

/* Flip whole components whose area-weighted, min-cos-gated radial sign disagrees
 * with the largest (by area) component's -- obj_orient_audit's EXACT
 * WINDING-FLIPPED criterion (flip iff agr >= agree_min && sign != R0).
 *
 * Flipping a whole connected component reverses only its INTERNAL edges (both
 * faces of every shared edge reverse together) and it shares no edge with any
 * other component, so this can NEVER change same_dir or manifoldness -- it is a
 * pure, safe orientation normalization. This is what drives WINDING-FLIPPED to 0
 * (OrientWeld's spatial vote is more conservative and leaves some behind).
 * Returns the number of components flipped. */
static size_t flip_audit_winding(const float *verts, size_t nv, int32_t *faces,
                                 size_t nf, const float axp[3], const double ad[3],
                                 double agree_min)
{
    int32_t *uf = NULL, *root2dense = NULL, *comp_of_vert = NULL;
    double *aout = NULL, *ain = NULL, *atot = NULL;
    uint8_t *flip = NULL;
    size_t i = 0, f = 0, nflip = 0;
    int32_t ncomp = 0, best = 0, R0 = 0;
    double bestA = -1.0;

    if (nv == 0 || nf == 0 || nv > (size_t)INT32_MAX) return 0;
    uf = (int32_t *)malloc(nv * sizeof(int32_t));
    root2dense = (int32_t *)malloc(nv * sizeof(int32_t));
    comp_of_vert = (int32_t *)malloc(nv * sizeof(int32_t));
    if (!uf || !root2dense || !comp_of_vert) { free(uf); free(root2dense); free(comp_of_vert); return 0; }

    for (i = 0; i < nv; i++) { uf[i] = (int32_t)i; root2dense[i] = -1; }
    for (f = 0; f < nf; f++) {
        int32_t a = faces[f * 3 + 0], b = faces[f * 3 + 1], c = faces[f * 3 + 2];
        int32_t ra = uf_find_local(uf, a), rb = uf_find_local(uf, b), rc;
        if (ra != rb) uf[rb] = ra;
        ra = uf_find_local(uf, a); rc = uf_find_local(uf, c);
        if (ra != rc) uf[rc] = ra;
    }
    for (f = 0; f < nf; f++) {
        int32_t r = uf_find_local(uf, faces[f * 3 + 0]);
        if (root2dense[r] < 0) root2dense[r] = ncomp++;
    }
    if (ncomp < 1) { free(uf); free(root2dense); free(comp_of_vert); return 0; }
    for (i = 0; i < nv; i++) comp_of_vert[i] = root2dense[uf_find_local(uf, (int32_t)i)];
    free(uf); free(root2dense);

    aout = (double *)calloc((size_t)ncomp, sizeof(double));
    ain  = (double *)calloc((size_t)ncomp, sizeof(double));
    atot = (double *)calloc((size_t)ncomp, sizeof(double));
    flip = (uint8_t *)calloc((size_t)ncomp, 1);
    if (!aout || !ain || !atot || !flip) { free(comp_of_vert); free(aout); free(ain); free(atot); free(flip); return 0; }

    for (f = 0; f < nf; f++) {
        int32_t ci = comp_of_vert[faces[f * 3 + 0]];
        double area = 0.0;
        int s = 0;
        if (ci < 0) continue;
        s = face_radial(verts, faces, f, axp, ad, &area);
        atot[ci] += area;
        if (s > 0) aout[ci] += area;
        else if (s < 0) ain[ci] += area;
    }

    for (i = 0; i < (size_t)ncomp; i++)
        if (atot[i] > bestA) { bestA = atot[i]; best = (int32_t)i; }
    R0 = (aout[best] >= ain[best]) ? 1 : -1;

    for (i = 0; i < (size_t)ncomp; i++) {
        double tot = aout[i] + ain[i];
        double agr = (tot > 0.0) ? (aout[i] > ain[i] ? aout[i] : ain[i]) / tot : 0.0;
        int sign = (aout[i] >= ain[i]) ? 1 : -1;
        if (tot > 0.0 && agr >= agree_min && sign != R0) { flip[i] = 1; nflip++; }
    }
    if (nflip > 0) {
        for (f = 0; f < nf; f++) {
            int32_t ci = comp_of_vert[faces[f * 3 + 0]];
            if (ci >= 0 && flip[ci]) {
                int32_t t = faces[f * 3 + 1];
                faces[f * 3 + 1] = faces[f * 3 + 2];
                faces[f * 3 + 2] = t;
            }
        }
    }
    free(comp_of_vert); free(aout); free(ain); free(atot); free(flip);
    return nflip;
}

/* --------------------------------------------------- non-orientable knot cut */

typedef struct { uint64_t key; int32_t face; uint8_t dir; } HalfEdge;

static int cmp_he(const void *pa, const void *pb)
{
    const HalfEdge *a = (const HalfEdge *)pa, *b = (const HalfEdge *)pb;
    if (a->key < b->key) return -1;
    if (a->key > b->key) return 1;
    return 0;
}

/* Locate `same_dir` interior edges (2 faces traversing the shared edge the SAME
 * way -- a winding frustration) and delete one face of each until same_dir == 0.
 * Victim selection is a greedy vertex cover: prefer the face touching the MOST
 * frustrated edges (so one deletion clears several knots), breaking ties by
 * smaller area then lower index. Deleting a face can only REDUCE edge
 * multiplicities, so this is monotone and terminates. Faces are compacted in
 * place; *pnf shrinks. Leaves a tiny triangular notch per cut -- negligible, and
 * re-capping the same 3 verts would only re-introduce the knot. Returns the
 * number of faces deleted. */
static size_t find_and_cut_knots(const float *verts, int32_t *faces,
                                 size_t *pnf, int max_iters)
{
    size_t total = 0;
    int iter = 0;
    for (iter = 0; iter < max_iters; iter++) {
        size_t nf = *pnf, i = 0, hn = 0, nvict = 0, w = 0, nse = 0;
        HalfEdge *he = NULL;
        uint8_t *victim = NULL;
        int32_t *ea = NULL, *eb = NULL, *deg = NULL;

        if (nf == 0) break;
        he = (HalfEdge *)malloc(nf * 3 * sizeof(HalfEdge));
        if (he == NULL) break;
        for (i = 0; i < nf; i++) {
            int e = 0;
            for (e = 0; e < 3; e++) {
                uint32_t a = (uint32_t)faces[i * 3 + e];
                uint32_t b = (uint32_t)faces[i * 3 + (e + 1) % 3];
                uint32_t lo = a < b ? a : b, hi = a < b ? b : a;
                he[hn].key = ((uint64_t)lo << 32) | (uint64_t)hi;
                he[hn].dir = (uint8_t)(a < b ? 0 : 1);
                he[hn].face = (int32_t)i;
                hn++;
            }
        }
        qsort(he, hn, sizeof(HalfEdge), cmp_he);

        /* Pass 1: collect same_dir edge face-pairs + per-face frustration degree. */
        ea  = (int32_t *)malloc((nf * 3 / 2 + 1) * sizeof(int32_t));
        eb  = (int32_t *)malloc((nf * 3 / 2 + 1) * sizeof(int32_t));
        deg = (int32_t *)calloc(nf, sizeof(int32_t));
        if (ea == NULL || eb == NULL || deg == NULL) { free(he); free(ea); free(eb); free(deg); break; }
        for (i = 0; i + 1 < hn; ) {
            size_t j = i + 1;
            while (j < hn && he[j].key == he[i].key) j++;
            if (j - i == 2 && he[i].dir == he[i + 1].dir) {   /* same_dir edge */
                ea[nse] = he[i].face; eb[nse] = he[i + 1].face; nse++;
                deg[he[i].face]++; deg[he[i + 1].face]++;
            }
            i = j;
        }
        free(he);
        if (nse == 0) { free(ea); free(eb); free(deg); break; }

        /* Pass 2: greedy cover -- one victim per still-uncovered frustrated edge. */
        victim = (uint8_t *)calloc(nf, 1);
        if (victim == NULL) { free(ea); free(eb); free(deg); break; }
        for (i = 0; i < nse; i++) {
            int32_t fa = ea[i], fb = eb[i], vic = fa;
            if (victim[fa] || victim[fb]) continue;   /* edge already covered */
            if (deg[fa] != deg[fb]) {
                vic = (deg[fa] > deg[fb]) ? fa : fb;
            } else {
                double aa = face_area(verts, faces, (size_t)fa);
                double ab = face_area(verts, faces, (size_t)fb);
                vic = (aa < ab) ? fa : (ab < aa) ? fb : (fa < fb ? fa : fb);
            }
            victim[vic] = 1; nvict++;
        }
        free(ea); free(eb); free(deg);

        for (i = 0; i < nf; i++) {
            if (victim[i]) continue;
            if (w != i) {
                faces[w * 3 + 0] = faces[i * 3 + 0];
                faces[w * 3 + 1] = faces[i * 3 + 1];
                faces[w * 3 + 2] = faces[i * 3 + 2];
            }
            w++;
        }
        *pnf = w;
        total += nvict;
        free(victim);
    }
    return total;
}

/* ---------------- intra-component radial-region diagnostics ----------------
 * Component-level flips cannot fix a backward area INSIDE a connected
 * component: consistent winding ties it to the bulk through whatever fold or
 * fused neck it hangs from, so it inherits the opposite radial facing while
 * same_dir stays 0. Segment decisive-sign faces into edge-connected same-sign
 * REGIONS and report the ones facing against the mesh majority -- the residual
 * "wrong-colour sheets" a two-sided shader still shows after all component-
 * level repair. crease = manifold edges straight into the opposite sign (a
 * hard 180-degree interface); neutral = edges into gate-failed faces (a fold
 * band or axis-parallel wall the region is connected through). */
typedef struct {
    int32_t id;
    int32_t sign;          /* +1 radially out, -1 in */
    size_t  nfaces;
    double  area;
    double  bb[6];         /* zmin,ymin,xmin, zmax,ymax,xmax */
    size_t  crease;
    size_t  neutral;
} RadRegion;

static int cmp_rr_area(const void *pa, const void *pb)
{
    const RadRegion *a = (const RadRegion *)pa, *b = (const RadRegion *)pb;
    if (a->area != b->area) return a->area < b->area ? 1 : -1;
    return a->id < b->id ? -1 : (a->id > b->id);
}

/* Shared segmentation state for the report + the cut-and-flip repair. */
typedef struct {
    int8_t   *sg;      /* per-face radial sign (+1/-1/0) */
    double   *farea;   /* per-face area */
    int32_t  *rg;      /* per-face region id (-1 = neutral face) */
    int32_t  *ea, *eb; /* manifold-edge face pairs */
    size_t    npair;
    RadRegion *rr;     /* UNSORTED: rr[rg[f]] is f's region */
    size_t    nrr;
    double    aout, ain;
} RadSeg;

static void radseg_free(RadSeg *S)
{
    free(S->sg); free(S->farea); free(S->rg); free(S->ea); free(S->eb); free(S->rr);
    memset(S, 0, sizeof *S);
}

static int radseg_build(const float *verts, size_t nv, const int32_t *faces,
                        size_t nf, const float axp[3], const double ad[3],
                        RadSeg *S)
{
    HalfEdge *he = NULL;
    int32_t *stk = NULL, *deg = NULL, *adj = NULL;
    size_t *off = NULL;
    size_t i = 0, f = 0, hn = 0, caprr = 0;
    int ok = -1;

    memset(S, 0, sizeof *S);
    if (nf == 0 || nv == 0) return -1;
    S->sg = (int8_t *)malloc(nf);
    S->farea = (double *)malloc(nf * sizeof(double));
    S->rg = (int32_t *)malloc(nf * sizeof(int32_t));
    S->ea = (int32_t *)malloc((nf * 3 / 2 + 1) * sizeof(int32_t));
    S->eb = (int32_t *)malloc((nf * 3 / 2 + 1) * sizeof(int32_t));
    he = (HalfEdge *)malloc(nf * 3 * sizeof(HalfEdge));
    if (!S->sg || !S->farea || !S->rg || !S->ea || !S->eb || !he) goto done;

    for (f = 0; f < nf; f++) {
        S->sg[f] = (int8_t)face_radial(verts, faces, f, axp, ad, &S->farea[f]);
        if (S->sg[f] > 0) S->aout += S->farea[f];
        else if (S->sg[f] < 0) S->ain += S->farea[f];
        S->rg[f] = -1;
    }

    /* manifold edge pairs (mesh is edge-manifold here; runs of 2) */
    for (f = 0; f < nf; f++) {
        int e = 0;
        for (e = 0; e < 3; e++) {
            uint32_t a = (uint32_t)faces[f * 3 + e];
            uint32_t b = (uint32_t)faces[f * 3 + (e + 1) % 3];
            uint32_t lo = a < b ? a : b, hi = a < b ? b : a;
            he[hn].key = ((uint64_t)lo << 32) | (uint64_t)hi;
            he[hn].dir = 0;
            he[hn].face = (int32_t)f;
            hn++;
        }
    }
    qsort(he, hn, sizeof(HalfEdge), cmp_he);
    for (i = 0; i + 1 < hn; ) {
        size_t j = i + 1;
        while (j < hn && he[j].key == he[i].key) j++;
        if (j - i == 2) {
            S->ea[S->npair] = he[i].face;
            S->eb[S->npair] = he[i + 1].face;
            S->npair++;
        }
        i = j;
    }
    free(he); he = NULL;

    /* CSR face adjacency */
    deg = (int32_t *)calloc(nf, sizeof(int32_t));
    off = (size_t *)malloc((nf + 1) * sizeof(size_t));
    adj = (int32_t *)malloc((S->npair * 2 + 1) * sizeof(int32_t));
    stk = (int32_t *)malloc(nf * sizeof(int32_t));
    if (!deg || !off || !adj || !stk) goto done;
    for (i = 0; i < S->npair; i++) { deg[S->ea[i]]++; deg[S->eb[i]]++; }
    off[0] = 0;
    for (f = 0; f < nf; f++) off[f + 1] = off[f] + (size_t)deg[f];
    memset(deg, 0, nf * sizeof(int32_t));
    for (i = 0; i < S->npair; i++) {
        adj[off[S->ea[i]] + (size_t)deg[S->ea[i]]++] = S->eb[i];
        adj[off[S->eb[i]] + (size_t)deg[S->eb[i]]++] = S->ea[i];
    }

    /* grow same-sign regions */
    for (f = 0; f < nf; f++) {
        size_t top = 0;
        RadRegion R;
        if (S->sg[f] == 0 || S->rg[f] >= 0) continue;
        memset(&R, 0, sizeof R);
        R.id = (int32_t)S->nrr;
        R.sign = S->sg[f];
        R.bb[0] = R.bb[1] = R.bb[2] = 1e30;
        R.bb[3] = R.bb[4] = R.bb[5] = -1e30;
        stk[top++] = (int32_t)f;
        S->rg[f] = R.id;
        while (top > 0) {
            int32_t cf = stk[--top];
            size_t t = 0;
            int k = 0;
            R.nfaces++;
            R.area += S->farea[cf];
            for (k = 0; k < 3; k++) {
                const float *p = &verts[(size_t)faces[(size_t)cf * 3 + k] * 3];
                int d = 0;
                for (d = 0; d < 3; d++) {
                    if ((double)p[d] < R.bb[d]) R.bb[d] = p[d];
                    if ((double)p[d] > R.bb[3 + d]) R.bb[3 + d] = p[d];
                }
            }
            for (t = off[cf]; t < off[cf] + (size_t)deg[cf]; t++) {
                int32_t nb = adj[t];
                if (S->rg[nb] >= 0 || S->sg[nb] != R.sign) continue;
                S->rg[nb] = R.id;
                stk[top++] = nb;
            }
        }
        if (S->nrr == caprr) {
            RadRegion *nrr_p = NULL;
            caprr = caprr ? caprr * 2 : 256;
            nrr_p = (RadRegion *)realloc(S->rr, caprr * sizeof(RadRegion));
            if (nrr_p == NULL) goto done;
            S->rr = nrr_p;
        }
        S->rr[S->nrr++] = R;
    }

    /* attachment stats */
    for (i = 0; i < S->npair; i++) {
        int32_t fa2 = S->ea[i], fb2 = S->eb[i];
        if (S->sg[fa2] != 0 && S->sg[fb2] != 0 && S->sg[fa2] != S->sg[fb2]) {
            S->rr[S->rg[fa2]].crease++;
            S->rr[S->rg[fb2]].crease++;
        } else if (S->sg[fa2] != 0 && S->sg[fb2] == 0) S->rr[S->rg[fa2]].neutral++;
        else if (S->sg[fb2] != 0 && S->sg[fa2] == 0) S->rr[S->rg[fb2]].neutral++;
    }
    ok = 0;

done:
    free(he); free(deg); free(off); free(adj); free(stk);
    if (ok != 0) radseg_free(S);
    return ok;
}

/* Returns the fraction of decisive area facing against R0 (0 when clean).
 * Prints the top backward regions; optionally writes a per-vertex-coloured
 * OBJ (backward=red, majority=grey, neutral=dark). */
static double radial_regions(const float *verts, size_t nv, const int32_t *faces,
                             size_t nf, const float axp[3], const double ad[3],
                             int R0, size_t max_rows, const char *dump_path)
{
    RadSeg S;
    size_t i = 0, f = 0, shown = 0, nback = 0, nback_big = 0;
    double dec = 0.0, opp = 0.0, ret = 0.0;

    if (radseg_build(verts, nv, faces, nf, axp, ad, &S) != 0) return 0.0;
    dec = S.aout + S.ain;
    opp = (R0 > 0) ? S.ain : S.aout;
    for (i = 0; i < S.nrr; i++) {
        if ((int)S.rr[i].sign != R0) {
            nback++;
            if (S.rr[i].nfaces >= 1000) nback_big++;
        }
    }
    fprintf(stderr,
        "[regions] decisive area out=%.3g in=%.3g (majority=%s) | BACKWARD %.2f%% "
        "in %zu region(s), %zu with >=1k faces\n",
        S.aout, S.ain, R0 > 0 ? "out" : "in",
        dec > 0.0 ? 100.0 * opp / dec : 0.0, nback, nback_big);
    if (S.nrr > 0 && max_rows > 0) {
        qsort(S.rr, S.nrr, sizeof(RadRegion), cmp_rr_area);   /* report-only */
        for (i = 0; i < S.nrr && shown < max_rows; i++) {
            if ((int)S.rr[i].sign == R0) continue;   /* backward rows only */
            fprintf(stderr,
                "  [bwd %2zu] %s f=%7zu area=%9.0f (%5.2f%% dec) z[%.0f,%.0f] "
                "y[%.0f,%.0f] x[%.0f,%.0f] crease=%zu neutral=%zu\n",
                shown, S.rr[i].sign > 0 ? "out" : "in ",
                S.rr[i].nfaces, S.rr[i].area,
                dec > 0.0 ? 100.0 * S.rr[i].area / dec : 0.0,
                S.rr[i].bb[0], S.rr[i].bb[3], S.rr[i].bb[1], S.rr[i].bb[4],
                S.rr[i].bb[2], S.rr[i].bb[5], S.rr[i].crease, S.rr[i].neutral);
            shown++;
        }
    }

    if (dump_path != NULL) {
        float *cols = (float *)malloc(nv * 3 * sizeof(float));
        if (cols != NULL) {
            for (i = 0; i < nv; i++) {
                cols[i * 3 + 0] = 0.38f; cols[i * 3 + 1] = 0.38f; cols[i * 3 + 2] = 0.42f;
            }
            for (f = 0; f < nf; f++) {   /* majority first, backward overrides */
                int k = 0;
                if (S.sg[f] == 0 || (int)S.sg[f] != R0) continue;
                for (k = 0; k < 3; k++) {
                    size_t vi = (size_t)faces[f * 3 + k];
                    cols[vi * 3 + 0] = 0.70f; cols[vi * 3 + 1] = 0.70f; cols[vi * 3 + 2] = 0.66f;
                }
            }
            for (f = 0; f < nf; f++) {
                int k = 0;
                if (S.sg[f] == 0 || (int)S.sg[f] == R0) continue;
                for (k = 0; k < 3; k++) {
                    size_t vi = (size_t)faces[f * 3 + k];
                    cols[vi * 3 + 0] = 0.92f; cols[vi * 3 + 1] = 0.13f; cols[vi * 3 + 2] = 0.13f;
                }
            }
            if (ObjIO_write_per_vertex_color(dump_path, verts, nv, faces, nf, cols) == 0)
                fprintf(stderr, "[regions] wrote %s (backward=red, majority=grey, "
                        "neutral=dark)\n", dump_path);
            free(cols);
        }
    }
    ret = dec > 0.0 ? opp / dec : 0.0;
    radseg_free(&S);
    return ret;
}

static int cmp_u64(const void *pa, const void *pb)
{
    uint64_t a = *(const uint64_t *)pa, b = *(const uint64_t *)pb;
    return a < b ? -1 : (a > b);
}

/* Undirected-edge key (lo<<32|hi) of the edge shared by manifold-pair faces
 * fa,fb; 0 if they share fewer than 2 vertices (cannot happen for a pair). */
static uint64_t pair_edge_key(const int32_t *faces, int32_t fa, int32_t fb)
{
    int32_t s[2] = { 0, 0 };
    int ns = 0, i = 0, j = 0;
    for (i = 0; i < 3 && ns < 2; i++) {
        int32_t va = faces[(size_t)fa * 3 + i];
        for (j = 0; j < 3; j++) {
            if (faces[(size_t)fb * 3 + j] == va) { s[ns++] = va; break; }
        }
    }
    if (ns < 2) return 0;
    {
        uint32_t lo = (uint32_t)(s[0] < s[1] ? s[0] : s[1]);
        uint32_t hi = (uint32_t)(s[0] < s[1] ? s[1] : s[0]);
        return ((uint64_t)lo << 32) | (uint64_t)hi;
    }
}

static int key_in(const uint64_t *keys, size_t n, uint64_t k)
{
    size_t lo = 0, hi = n;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (keys[mid] < k) lo = mid + 1; else hi = mid;
    }
    return lo < n && keys[lo] == k;
}

/* Fans of v's incident faces connected only through NON-cut edges at v.
 * fuf is a face-indexed union-find scratch (size nf); returns the fan count
 * with roots left in fuf for the incident slots. */
static int fans_at_vertex(const int32_t *faces, int32_t v,
                          const int32_t *inc, size_t s, size_t e,
                          const uint64_t *cutk, size_t ncut, int32_t *fuf)
{
    size_t j = 0, q = 0;
    int fans = 0;
    for (j = s; j < e; j++) fuf[inc[j]] = inc[j];
    for (j = s; j < e; j++) {
        int32_t fi = inc[j];
        for (q = j + 1; q < e; q++) {
            int32_t fj = inc[q];
            int k = 0, joined = 0;
            for (k = 0; k < 3 && !joined; k++) {
                int32_t w = faces[(size_t)fi * 3 + k];
                int t = 0;
                if (w == v) continue;
                for (t = 0; t < 3; t++) {
                    if (faces[(size_t)fj * 3 + t] == w) {
                        uint32_t lo = (uint32_t)(v < w ? v : w);
                        uint32_t hi = (uint32_t)(v < w ? w : v);
                        if (!key_in(cutk, ncut,
                                    ((uint64_t)lo << 32) | (uint64_t)hi))
                            joined = 1;
                        break;
                    }
                }
            }
            if (joined) {
                int32_t ra = uf_find_local(fuf, fi), rb = uf_find_local(fuf, fj);
                if (ra != rb) fuf[rb] = ra;
            }
        }
    }
    for (j = s; j < e; j++)
        if (uf_find_local(fuf, inc[j]) == inc[j]) fans++;
    return fans;
}

static int cmp_i32(const void *pa, const void *pb)
{
    int32_t a = *(const int32_t *)pa, b = *(const int32_t *)pb;
    return a < b ? -1 : (a > b);
}

/* Vertex-split CUT along a set of edges: disconnect face adjacency across
 * every cut edge WITHOUT deleting a face. For each vertex incident to a cut
 * edge, group its incident faces into fans connected only through non-cut
 * edges at that vertex; the first fan keeps the vertex, every other fan gets a
 * coincident duplicate (verts grown from the arena; geometry untouched).
 * Standard mesh seam cut -- preserves edge-manifoldness. Returns the number of
 * duplicated vertices. */
static size_t cut_edges_split(Arena_T arena, float **pverts, size_t *pnv,
                              int32_t *faces, size_t nf,
                              const uint64_t *cutk, size_t ncut)
{
    size_t nv = *pnv, na = 0, m = 0, extra = 0, next_new = 0;
    size_t i = 0, f = 0, a = 0, j = 0, ndup = 0;
    int32_t *A = NULL, *ainc = NULL, *fuf = NULL;
    size_t *aoff = NULL, *acur = NULL;
    float *nvert = NULL;

    if (ncut == 0 || nf == 0) return 0;

    /* affected vertices = endpoints of cut edges, sorted unique */
    A = (int32_t *)malloc(ncut * 2 * sizeof(int32_t));
    if (A == NULL) return 0;
    for (i = 0; i < ncut; i++) {
        A[na++] = (int32_t)(cutk[i] >> 32);
        A[na++] = (int32_t)(cutk[i] & 0xffffffffu);
    }
    qsort(A, na, sizeof(int32_t), cmp_i32);
    for (i = 0; i < na; i++)
        if (m == 0 || A[i] != A[m - 1]) A[m++] = A[i];
    na = m;

    /* incidence lists for affected vertices only */
    aoff = (size_t *)calloc(na + 1, sizeof(size_t));
    acur = (size_t *)malloc((na + 1) * sizeof(size_t));
    fuf  = (int32_t *)malloc(nf * sizeof(int32_t));
    if (aoff == NULL || acur == NULL || fuf == NULL) goto done;
    for (f = 0; f < nf; f++) {
        int k = 0;
        for (k = 0; k < 3; k++) {
            int32_t v = faces[f * 3 + k];
            size_t lo = 0, hi = na;
            while (lo < hi) {
                size_t mid = lo + (hi - lo) / 2;
                if (A[mid] < v) lo = mid + 1; else hi = mid;
            }
            if (lo < na && A[lo] == v) aoff[lo + 1]++;
        }
    }
    for (i = 0; i < na; i++) aoff[i + 1] += aoff[i];
    ainc = (int32_t *)malloc((aoff[na] + 1) * sizeof(int32_t));
    if (ainc == NULL) goto done;
    memcpy(acur, aoff, (na + 1) * sizeof(size_t));
    for (f = 0; f < nf; f++) {
        int k = 0;
        for (k = 0; k < 3; k++) {
            int32_t v = faces[f * 3 + k];
            size_t lo = 0, hi = na;
            while (lo < hi) {
                size_t mid = lo + (hi - lo) / 2;
                if (A[mid] < v) lo = mid + 1; else hi = mid;
            }
            if (lo < na && A[lo] == v) ainc[acur[lo]++] = (int32_t)f;
        }
    }

    /* pass 1: count extra fans */
    for (a = 0; a < na; a++) {
        int fans = fans_at_vertex(faces, A[a], ainc, aoff[a], aoff[a + 1],
                                  cutk, ncut, fuf);
        if (fans >= 2) extra += (size_t)(fans - 1);
    }
    if (extra == 0) goto done;

    /* pass 2: grow verts, recompute fans, assign ids, repoint corners */
    nvert = (float *)ARENA_ALLOC(arena, (long)((nv + extra) * 3 * sizeof(float)));
    memcpy(nvert, *pverts, nv * 3 * sizeof(float));
    next_new = nv;
    for (a = 0; a < na; a++) {
        int32_t v = A[a];
        size_t s = aoff[a], e = aoff[a + 1];
        int32_t seen_root[64], seen_id[64];
        int nseen = 0;
        (void)fans_at_vertex(faces, v, ainc, s, e, cutk, ncut, fuf);
        for (j = s; j < e; j++) {
            int32_t fi = ainc[j];
            int32_t r = uf_find_local(fuf, fi);
            int32_t oid = -1;
            int t = 0;
            for (t = 0; t < nseen; t++)
                if (seen_root[t] == r) { oid = seen_id[t]; break; }
            if (oid < 0) {
                if (nseen == 0 || nseen >= 64) {
                    oid = v;                       /* first fan keeps v */
                } else {
                    oid = (int32_t)next_new++;     /* extra fan -> duplicate */
                    nvert[(size_t)oid * 3 + 0] = nvert[(size_t)v * 3 + 0];
                    nvert[(size_t)oid * 3 + 1] = nvert[(size_t)v * 3 + 1];
                    nvert[(size_t)oid * 3 + 2] = nvert[(size_t)v * 3 + 2];
                }
                if (nseen < 64) {
                    seen_root[nseen] = r;
                    seen_id[nseen] = oid;
                    nseen++;
                }
            }
            if (oid != v) {
                for (t = 0; t < 3; t++) {
                    if (faces[(size_t)fi * 3 + t] == v) {
                        faces[(size_t)fi * 3 + t] = oid;
                        break;
                    }
                }
            }
        }
    }
    *pverts = nvert;
    *pnv = next_new;
    ndup = next_new - nv;

done:
    free(A); free(aoff); free(acur); free(ainc); free(fuf);
    return ndup;
}

/* Cut-and-flip repair for LARGE intra-component backward regions. A wrap-scale
 * connected backward area is an orientation-propagation error wired to the
 * bulk through a fold line or fused neck (a genuine fold flap shows up as many
 * SMALL regions, not one grid-spanning web -- hence the min_faces gate).
 * Repair: flood-fill from majority-sign faces (never entering a picked region)
 * to classify BULK vs FREED -- neutral fold-band appendages reachable only
 * through the region travel WITH it instead of shedding as debris -- then
 * vertex-split CUT the freed side along its interface (NO face deleted,
 * geometry untouched), split remaining vertex-only contacts, and flip the
 * freed component(s) via the audit criterion. Each structural backward web
 * becomes exactly ONE new component. The coincident seam cannot be re-welded:
 * gluing a flipped side back would recreate the same_dir defect the flip
 * removed (a fold interface is an orientability obstruction -- one component
 * through it and consistent facing across it are mutually exclusive).
 * Iterates until no qualifying region remains (max 4 rounds). */
static int fix_backward_regions(Arena_T arena, float **pverts, size_t *pnv,
                                int32_t *faces, size_t *pnf,
                                const float axp[3], const double ad[3],
                                size_t min_faces, double agree,
                                size_t *out_cuts, size_t *out_flipped,
                                size_t *out_dups)
{
    int round = 0;
    for (round = 0; round < 4; round++) {
        RadSeg S;
        uint8_t *picked = NULL, *freed = NULL, *reach = NULL;
        int32_t *deg = NULL, *adj = NULL, *stk = NULL;
        size_t *off = NULL;
        uint64_t *cutk = NULL;
        size_t i = 0, f = 0, npick = 0, nck = 0, nflip = 0, sp = 0, ndup = 0;
        size_t nf = *pnf;
        int R0m = mesh_radial_sign(*pverts, faces, nf, axp, ad);
        if (R0m == 0) break;
        if (radseg_build(*pverts, *pnv, faces, nf, axp, ad, &S) != 0) break;

        picked = (uint8_t *)calloc(S.nrr ? S.nrr : 1, 1);
        freed = (uint8_t *)calloc(nf, 1);
        reach = (uint8_t *)calloc(nf, 1);
        stk = (int32_t *)malloc(nf * sizeof(int32_t));
        deg = (int32_t *)calloc(nf, sizeof(int32_t));
        off = (size_t *)malloc((nf + 1) * sizeof(size_t));
        adj = (int32_t *)malloc((S.npair * 2 + 1) * sizeof(int32_t));
        cutk = (uint64_t *)malloc((S.npair + 1) * sizeof(uint64_t));
        if (!picked || !freed || !reach || !stk || !deg || !off || !adj || !cutk) {
            free(picked); free(freed); free(reach); free(stk);
            free(deg); free(off); free(adj); free(cutk);
            radseg_free(&S);
            break;
        }
        for (i = 0; i < S.nrr; i++) {
            if ((int)S.rr[i].sign != R0m && S.rr[i].nfaces >= min_faces) {
                picked[i] = 1;
                npick++;
            }
        }
        if (npick == 0) {
            free(picked); free(freed); free(reach); free(stk);
            free(deg); free(off); free(adj); free(cutk);
            radseg_free(&S);
            break;
        }
        for (f = 0; f < nf; f++)
            if (S.rg[f] >= 0 && picked[S.rg[f]]) freed[f] = 1;

        /* CSR adjacency from the seg's manifold pairs */
        for (i = 0; i < S.npair; i++) { deg[S.ea[i]]++; deg[S.eb[i]]++; }
        off[0] = 0;
        for (f = 0; f < nf; f++) off[f + 1] = off[f] + (size_t)deg[f];
        memset(deg, 0, nf * sizeof(int32_t));
        for (i = 0; i < S.npair; i++) {
            adj[off[S.ea[i]] + (size_t)deg[S.ea[i]]++] = S.eb[i];
            adj[off[S.eb[i]] + (size_t)deg[S.eb[i]]++] = S.ea[i];
        }

        /* flood: BULK = reachable from majority-sign faces without entering
         * a picked region */
        for (f = 0; f < nf; f++) {
            size_t top = 0;
            if (freed[f] || reach[f] || S.sg[f] != (int8_t)R0m) continue;
            reach[f] = 1;
            stk[top++] = (int32_t)f;
            while (top > 0) {
                int32_t cf = stk[--top];
                size_t t = 0;
                for (t = off[cf]; t < off[cf] + (size_t)deg[cf]; t++) {
                    int32_t nb = adj[t];
                    if (reach[nb] || freed[nb]) continue;
                    reach[nb] = 1;
                    stk[top++] = nb;
                }
            }
        }
        /* captives (appendages unreachable from the bulk) travel with the region */
        for (f = 0; f < nf; f++)
            if (!reach[f] && !freed[f]) freed[f] = 1;

        /* cut edges = freed/bulk interface */
        for (i = 0; i < S.npair; i++) {
            if (freed[S.ea[i]] != freed[S.eb[i]]) {
                uint64_t k = pair_edge_key(faces, S.ea[i], S.eb[i]);
                if (k != 0) cutk[nck++] = k;
            }
        }
        qsort(cutk, nck, sizeof(uint64_t), cmp_u64);
        *out_cuts += nck;

        free(picked); free(freed); free(reach); free(stk);
        free(deg); free(off); free(adj);
        radseg_free(&S);

        ndup = cut_edges_split(arena, pverts, pnv, faces, nf, cutk, nck);
        *out_dups += ndup;
        free(cutk);

        {   /* vertex-only contacts (seam endpoints, isolated touches) */
            ComponentMesh cm;
            memset(&cm, 0, sizeof cm);
            cm.verts = *pverts; cm.faces = faces; cm.nv = *pnv; cm.nf = nf;
            cm.self = &cm;
            PinholeFill_split_pinches(arena, &cm, 1, &sp);
            *pverts = cm.verts;
            *pnv = cm.nv;
            *out_dups += sp;
        }
        nflip = flip_audit_winding(*pverts, *pnv, faces, nf, axp, ad, agree);
        *out_flipped += nflip;
        if (nck == 0 && nflip == 0) break;
    }
    return round;
}

/* ------------------------------------------------------------------ driver */

static int reorient_run(Arena_T arena, float **pverts, size_t *pnv, int32_t *faces,
                        size_t *pnf, const float axp[3], const float axd[3],
                        float radius, int orient_first, int do_cut,
                        size_t fix_min, ReorientStats *st)
{
    float *verts = *pverts;
    size_t nv = *pnv;
    double ad[3] = { 0.0, 0.0, 0.0 };
    double L = 0.0;
    int have_axis = 0, input_sign = 0, out_sign = 0;
    MeshManifoldStats s0, s1;

    memset(st, 0, sizeof *st);
    st->nf_in = *pnf;
    if (nv == 0 || *pnf == 0) return 0;

    L = sqrt((double)axd[0] * axd[0] + (double)axd[1] * axd[1] + (double)axd[2] * axd[2]);
    if (L > 1e-9) { ad[0] = axd[0] / L; ad[1] = axd[1] / L; ad[2] = axd[2] / L; have_axis = 1; }

    s0 = MeshManifold_audit(arena, nv, faces, *pnf);
    st->same_dir_before = s0.same_dir_edges;
    st->nm_before = s0.nm_edges;
    st->nmv_before = s0.nm_verts;
    if (have_axis) input_sign = mesh_radial_sign(verts, faces, *pnf, axp, ad);

    if (orient_first)
        OrientMesh_consistent(arena, verts, nv, NULL, faces, *pnf,
                              &st->om_flips, &st->om_components, NULL);

    OrientWeld_components_axis(arena, verts, nv, faces, *pnf, radius,
                              have_axis ? axp : NULL, have_axis ? axd : NULL,
                              &st->ow_flipped, &st->ow_radial);

    /* Guarantee the audit criterion (WINDING-FLIPPED -> 0): flip every component
     * whose radial sign decisively disagrees with the bulk. Safe after OrientWeld
     * -- whole-component flips never touch same_dir or manifoldness. */
    if (have_axis)
        st->audit_flipped = flip_audit_winding(verts, nv, faces, *pnf, axp, ad, 0.8);

    /* Intra-component backward regions (wrap-scale webs wired through folds /
     * fused necks): seam-cut them free (no faces deleted) and flip them. */
    if (fix_min > 0 && have_axis) {
        st->region_rounds = (size_t)fix_backward_regions(arena, &verts, &nv,
                                faces, pnf, axp, ad, fix_min, 0.8,
                                &st->region_cuts, &st->region_flipped,
                                &st->region_dups);
    }

    if (do_cut) {
        st->knots_cut = find_and_cut_knots(verts, faces, pnf, 4);
        /* Deleting a knot face can orphan a vertex into a bowtie (pinch). Split
         * those so the result stays 2-manifold (nm_vert=0) rather than trading a
         * cosmetic same_dir for a hard-fail pinch. split_pinch_verts grows verts
         * from the arena and repoints faces in place (face count unchanged). */
        if (st->knots_cut > 0) {
            ComponentMesh cm;
            memset(&cm, 0, sizeof cm);
            cm.verts = verts; cm.faces = faces; cm.nv = nv; cm.nf = *pnf;
            cm.self = &cm;
            PinholeFill_split_pinches(arena, &cm, 1, &st->pinch_splits);
            verts = cm.verts;   /* may have grown */
            nv = cm.nv;
        }
    }

    /* Keep the mesh's global in/out convention: a whole-mesh flip preserves
     * winding consistency (every edge's two faces reverse together), so this is
     * cosmetic-only and never touches same_dir. */
    if (have_axis) {
        out_sign = mesh_radial_sign(verts, faces, *pnf, axp, ad);
        if (input_sign != 0 && out_sign != 0 && out_sign != input_sign) {
            flip_all(faces, *pnf);
            st->global_flip = 1;
        }
    }

    s1 = MeshManifold_audit(arena, nv, faces, *pnf);
    st->same_dir_after = s1.same_dir_edges;
    st->nm_after = s1.nm_edges;
    st->nmv_after = s1.nm_verts;
    st->nf_out = *pnf;
    *pverts = verts;
    *pnv = nv;
    return 0;
}

static void print_stats(const ReorientStats *st)
{
    fprintf(stderr,
        "[obj_reorient] faces %zu -> %zu | OrientMesh flips=%zu comps=%zu | "
        "OrientWeld flipped=%zu (radial %zu) | audit flipped=%zu | "
        "region-fix rounds=%zu cut_edges=%zu flips=%zu dup_verts=%zu | "
        "knots cut=%zu (pinch splits=%zu) | global_flip=%d\n"
        "               same_dir %zu -> %zu | nm_edge %zu -> %zu | nm_vert %zu -> %zu\n",
        st->nf_in, st->nf_out, st->om_flips, st->om_components,
        st->ow_flipped, st->ow_radial, st->audit_flipped,
        st->region_rounds, st->region_cuts, st->region_flipped,
        st->region_dups, st->knots_cut,
        st->pinch_splits, st->global_flip,
        st->same_dir_before, st->same_dir_after, st->nm_before, st->nm_after,
        st->nmv_before, st->nmv_after);
}

/* -------------------------------------------------------------- selftest */

#define CHECK(cond, msg) \
    do { if (!(cond)) { fprintf(stderr, "  FAIL: %s\n", (msg)); nfail++; } \
         else { fprintf(stderr, "  ok: %s\n", (msg)); } } while (0)

static void *xgrow(void *p, size_t need, size_t *cap, size_t elem)
{
    size_t ncap = *cap;
    if (need <= ncap) return p;
    ncap = (ncap == 0) ? 1024 : ncap;
    while (ncap < need) ncap *= 2;
    p = realloc(p, ncap * elem);
    if (p == NULL) { fprintf(stderr, "OOM\n"); exit(1); }
    *cap = ncap;
    return p;
}

/* Cylinder-arc sheet (z,y,x)=(zofs+j, rad cos th, rad sin th); outward winding
 * unless flip. Appends; returns first new vertex index in *out_base. */
static void build_arc(float rad, double th0, double th1, int nth, int nz,
                      float zofs, int flip,
                      float **verts, int32_t **faces, size_t *nv, size_t *nf,
                      size_t *cap_v, size_t *cap_f, size_t *out_base)
{
    int i = 0, j = 0;
    size_t base = *nv;
    *out_base = base;
    for (j = 0; j <= nz; j++) {
        for (i = 0; i <= nth; i++) {
            double th = th0 + (th1 - th0) * (double)i / (double)nth;
            *verts = (float *)xgrow(*verts, (*nv + 1) * 3, cap_v, sizeof(float));
            (*verts)[*nv * 3 + 0] = zofs + (float)j;
            (*verts)[*nv * 3 + 1] = rad * (float)cos(th);
            (*verts)[*nv * 3 + 2] = rad * (float)sin(th);
            (*nv)++;
        }
    }
    for (j = 0; j < nz; j++) {
        for (i = 0; i < nth; i++) {
            int32_t a = (int32_t)(base + (size_t)j * (size_t)(nth + 1) + (size_t)i);
            int32_t b = a + 1, c = a + nth + 1, d = c + 1;
            int32_t t1[3], t2[3];
            t1[0] = a; t1[1] = b; t1[2] = c;
            t2[0] = b; t2[1] = d; t2[2] = c;
            if (flip) {
                int32_t t = t1[1]; t1[1] = t1[2]; t1[2] = t;
                t = t2[1]; t2[1] = t2[2]; t2[2] = t;
            }
            *faces = (int32_t *)xgrow(*faces, (*nf + 2) * 3, cap_f, sizeof(int32_t));
            memcpy(&(*faces)[*nf * 3], t1, sizeof t1); (*nf)++;
            memcpy(&(*faces)[*nf * 3], t2, sizeof t2); (*nf)++;
        }
    }
}

/* Ribbon: ONE parametric sheet over a (theta, r) column path -- used to build a
 * hairpin (outbound leg, radial turn, return leg) whose return leg is a
 * winding-CONSISTENT backward region: same_dir = 0 by construction because the
 * whole strip is a single parametric grid. */
static void build_ribbon(const double *ths, const double *rs, int ncols, int nz,
                         float **verts, int32_t **faces, size_t *nv, size_t *nf,
                         size_t *cap_v, size_t *cap_f)
{
    int i = 0, j = 0;
    size_t base = *nv;
    for (j = 0; j <= nz; j++) {
        for (i = 0; i < ncols; i++) {
            *verts = (float *)xgrow(*verts, (*nv + 1) * 3, cap_v, sizeof(float));
            (*verts)[*nv * 3 + 0] = (float)j;
            (*verts)[*nv * 3 + 1] = (float)(rs[i] * cos(ths[i]));
            (*verts)[*nv * 3 + 2] = (float)(rs[i] * sin(ths[i]));
            (*nv)++;
        }
    }
    for (j = 0; j < nz; j++) {
        for (i = 0; i + 1 < ncols; i++) {
            int32_t a = (int32_t)(base + (size_t)j * (size_t)ncols + (size_t)i);
            int32_t b = a + 1, c = a + ncols, d = c + 1;
            int32_t t1[3], t2[3];
            t1[0] = a; t1[1] = b; t1[2] = c;
            t2[0] = b; t2[1] = d; t2[2] = c;
            *faces = (int32_t *)xgrow(*faces, (*nf + 2) * 3, cap_f, sizeof(int32_t));
            memcpy(&(*faces)[*nf * 3], t1, sizeof t1); (*nf)++;
            memcpy(&(*faces)[*nf * 3], t2, sizeof t2); (*nf)++;
        }
    }
}

static double radial_sign_of(const float *verts, const int32_t *faces,
                             size_t nf, size_t v0, size_t v1)
{
    double s = 0.0;
    size_t f = 0;
    for (f = 0; f < nf; f++) {
        size_t a = (size_t)faces[f * 3 + 0], b = (size_t)faces[f * 3 + 1],
               c = (size_t)faces[f * 3 + 2];
        double e1[3], e2[3], gn[3], cy, cx, rl;
        int k = 0;
        if (a < v0 || a >= v1) continue;
        for (k = 0; k < 3; k++) {
            e1[k] = (double)verts[b * 3 + k] - verts[a * 3 + k];
            e2[k] = (double)verts[c * 3 + k] - verts[a * 3 + k];
        }
        gn[0] = e1[1] * e2[2] - e1[2] * e2[1];
        gn[1] = e1[2] * e2[0] - e1[0] * e2[2];
        gn[2] = e1[0] * e2[1] - e1[1] * e2[0];
        cy = ((double)verts[a * 3 + 1] + verts[b * 3 + 1] + verts[c * 3 + 1]) / 3.0;
        cx = ((double)verts[a * 3 + 2] + verts[b * 3 + 2] + verts[c * 3 + 2]) / 3.0;
        rl = sqrt(cy * cy + cx * cx);
        if (rl < 1.0) continue;
        s += (gn[1] * cy + gn[2] * cx) / rl;
    }
    return s;
}

static int selftest(void)
{
    int nfail = 0;
    const float AXP[3] = { 0.0f, 0.0f, 0.0f };
    const float AXD[3] = { 1.0f, 0.0f, 0.0f };
    Arena_T arena = Arena_new();

    fprintf(stderr, "[selftest] obj_reorient\n");

    /* t1: knot cut in isolation. A 2x3 strip (4 tris) is orientable; flipping
     * one corner triangle (2 interior edges) creates same_dir=2, nm_edge=0.
     * find_and_cut_knots must delete that 1 face and reach same_dir=0. */
    {
        float v[18] = { 0,0,0,  0,0,1,  0,0,2,   0,1,0,  0,1,1,  0,1,2 };
        int32_t f[12] = { 0,1,4,  0,4,3,  1,2,5,  1,5,4 };
        size_t nf = 4;
        MeshManifoldStats a0, a1;
        /* flip tri 0 -> (0,4,1) */
        f[1] = 4; f[2] = 1;
        a0 = MeshManifold_audit(arena, 6, f, nf);
        CHECK(a0.same_dir_edges == 2 && a0.nm_edges == 0,
              "t1 setup: flipped corner tri -> same_dir=2, nm_edge=0");
        {
            size_t cut = find_and_cut_knots(v, f, &nf, 4);
            CHECK(cut == 1 && nf == 3, "t1 cut deletes 1 face");
        }
        a1 = MeshManifold_audit(arena, 6, f, nf);
        CHECK(a1.same_dir_edges == 0, "t1 same_dir 2 -> 0 after cut");
    }

    /* t2: cut is a no-op on an already-clean mesh. */
    {
        float *v = NULL; int32_t *fc = NULL;
        size_t nv = 0, nf = 0, cv = 0, cf = 0, base = 0, cut = 0;
        build_arc(50.0f, 0.2, 2.2, 12, 6, 0.0f, 0, &v, &fc, &nv, &nf, &cv, &cf, &base);
        cut = find_and_cut_knots(v, fc, &nf, 4);
        CHECK(cut == 0, "t2 no cut on clean mesh");
        free(v); free(fc);
    }

    /* t3: end-to-end. Two coaxial arcs, the far one detached + flipped. The full
     * pipeline must correct it (consistent winding, same_dir=0) and keep the
     * anchor's global sense. */
    {
        float *v = NULL; int32_t *fc = NULL;
        size_t nv = 0, nf = 0, cv = 0, cf = 0, baseA = 0, baseB = 0;
        ReorientStats st;
        build_arc(50.0f, 0.2, 2.2, 24, 8, 0.0f, 0, &v, &fc, &nv, &nf, &cv, &cf, &baseA);
        build_arc(50.0f, 3.4, 4.4, 10, 4, 200.0f, 1, &v, &fc, &nv, &nf, &cv, &cf, &baseB);
        CHECK(radial_sign_of(v, fc, nf, baseB, nv) < 0.0, "t3 setup: far comp inward");
        reorient_run(arena, &v, &nv, fc, &nf, AXP, AXD, 3.0f, 0, 1, 0, &st);
        CHECK(st.ow_flipped + st.audit_flipped >= 1, "t3 backward comp flipped");
        CHECK(st.same_dir_after == 0 && st.nm_after == 0, "t3 clean manifold out");
        CHECK(radial_sign_of(v, fc, nf, baseA, baseB) > 0.0
              && radial_sign_of(v, fc, nf, baseB, nv) > 0.0,
              "t3 both comps now outward (consistent, sense preserved)");
        {
            double adn[3] = { 1.0, 0.0, 0.0 };
            int R0m = mesh_radial_sign(v, fc, nf, AXP, adn);
            CHECK(R0m == 1
                  && radial_regions(v, nv, fc, nf, AXP, adn, R0m, 4, NULL) < 1e-12,
                  "t3b regions: no backward area after fix");
        }
        free(v); free(fc);
    }

    /* t4: idempotence -- a second run changes nothing. */
    {
        float *v = NULL; int32_t *fc = NULL;
        size_t nv = 0, nf = 0, cv = 0, cf = 0, baseA = 0, baseB = 0;
        ReorientStats st1, st2;
        build_arc(50.0f, 0.2, 2.2, 24, 8, 0.0f, 0, &v, &fc, &nv, &nf, &cv, &cf, &baseA);
        build_arc(50.0f, 3.4, 4.4, 10, 4, 200.0f, 1, &v, &fc, &nv, &nf, &cv, &cf, &baseB);
        reorient_run(arena, &v, &nv, fc, &nf, AXP, AXD, 3.0f, 0, 1, 0, &st1);
        reorient_run(arena, &v, &nv, fc, &nf, AXP, AXD, 3.0f, 0, 1, 0, &st2);
        CHECK(st2.ow_flipped == 0 && st2.audit_flipped == 0 && st2.knots_cut == 0,
              "t4 second run is a no-op");
        free(v); free(fc);
    }

    /* t6: hairpin ribbon -- a winding-CONSISTENT backward region inside ONE
     * component (the return leg of the hairpin), i.e. the intra-component case
     * that no component-level pass can reach. --fix-regions must rim-cut it
     * free, separate it, and flip it, keeping the mesh manifold. */
    {
        double ths[31], rs[31], adn[3] = { 1.0, 0.0, 0.0 };
        float *v = NULL, *v0 = NULL;
        int32_t *fc = NULL;
        size_t nv = 0, nf = 0, cv = 0, cf = 0;
        int i = 0, R0m = 0;
        double bad0 = 0.0, bad1 = 0.0;
        ReorientStats st;
        for (i = 0; i < 20; i++) { ths[i] = 0.2 + 2.0 * i / 19.0; rs[i] = 50.0; }
        for (i = 20; i < 23; i++) { ths[i] = 2.2 + 0.02 * (i - 19); rs[i] = 50.0 + 0.9 * (i - 19); }
        for (i = 23; i < 31; i++) { ths[i] = 2.2 - 1.3 * (i - 22) / 8.0; rs[i] = 53.0; }
        build_ribbon(ths, rs, 31, 6, &v, &fc, &nv, &nf, &cv, &cf);
        v0 = v;
        {
            MeshManifoldStats a = MeshManifold_audit(arena, nv, fc, nf);
            CHECK(a.same_dir_edges == 0 && a.nm_edges == 0,
                  "t6 setup: hairpin consistent + manifold");
        }
        R0m = mesh_radial_sign(v, fc, nf, AXP, adn);
        bad0 = radial_regions(v, nv, fc, nf, AXP, adn, R0m, 0, NULL);
        CHECK(bad0 > 0.15, "t6 setup: sizeable backward region");
        reorient_run(arena, &v, &nv, fc, &nf, AXP, AXD, 3.0f, 0, 1, 10, &st);
        CHECK(st.region_flipped >= 1 && st.region_cuts > 0,
              "t6 region cut+flip fired");
        CHECK(st.nf_in == st.nf_out, "t6 NO faces deleted (seam cut only)");
        CHECK(st.same_dir_after == 0 && st.nm_after == 0 && st.nmv_after == 0,
              "t6 output manifold");
        R0m = mesh_radial_sign(v, fc, nf, AXP, adn);
        bad1 = radial_regions(v, nv, fc, nf, AXP, adn, R0m, 0, NULL);
        CHECK(bad1 < 0.02, "t6 backward area repaired");
        free(v0); free(fc);   /* v may have moved into the arena via pinch-split */
    }

    Arena_dispose(&arena);
    fprintf(stderr, "[selftest] %s (%d failures)\n",
            nfail == 0 ? "ALL PASS" : "FAILURES", nfail);
    return nfail;
}

/* ------------------------------------------------------------------ main */

static void usage(void)
{
    fprintf(stderr,
        "usage: obj_reorient <in.obj> <out.obj> [options]\n"
        "       obj_reorient --selftest\n"
        "  Fix backward-wound components + non-orientable knots on a welded mesh.\n"
        "options:\n"
        "  --axis-point Z Y X   scroll axis point (default 0 3405 2878)\n"
        "  --axis-dir Z Y X     scroll axis direction (default 1 0 0)\n"
        "  --radius F           cross-component vote ball, vox (default 3.0)\n"
        "  --keep-knots         do NOT cut residual same_dir knots\n"
        "  --no-orient-first    skip the pre-pass OrientMesh consistency step\n"
        "  --dry-run            analyse + report, do not write out.obj\n"
        "  --regions            report intra-component backward radial regions\n"
        "  --dump-colored P     write per-vertex-coloured OBJ (backward=red)\n"
        "  --fix-regions N      cut+flip backward regions with >= N faces\n"
        "                       (vertex-split seam cut; NO faces deleted)\n");
}

int main(int argc, char **argv)
{
    const char *in_path = NULL, *out_path = NULL;
    float axp[3] = { 0.0f, 3405.0f, 2878.0f };
    float axd[3] = { 1.0f, 0.0f, 0.0f };
    float radius = 3.0f;
    int do_cut = 1, orient_first = 1, dry_run = 0, show_regions = 0, i = 0;
    size_t fix_min = 0;
    const char *dump_colored = NULL;
    Arena_T arena = NULL;
    float *verts = NULL; int32_t *faces = NULL;
    size_t nv = 0, nf = 0;
    ReorientStats st;

    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0)
        return selftest() == 0 ? 0 : 1;
    if (argc < 3) { usage(); return 1; }
    in_path = argv[1];
    out_path = argv[2];
    for (i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--axis-point") == 0 && i + 3 < argc) {
            axp[0] = (float)strtod(argv[++i], NULL);
            axp[1] = (float)strtod(argv[++i], NULL);
            axp[2] = (float)strtod(argv[++i], NULL);
        } else if (strcmp(argv[i], "--axis-dir") == 0 && i + 3 < argc) {
            axd[0] = (float)strtod(argv[++i], NULL);
            axd[1] = (float)strtod(argv[++i], NULL);
            axd[2] = (float)strtod(argv[++i], NULL);
        } else if (strcmp(argv[i], "--radius") == 0 && i + 1 < argc) {
            radius = (float)strtod(argv[++i], NULL);
        } else if (strcmp(argv[i], "--keep-knots") == 0) {
            do_cut = 0;
        } else if (strcmp(argv[i], "--no-orient-first") == 0) {
            orient_first = 0;
        } else if (strcmp(argv[i], "--dry-run") == 0) {
            dry_run = 1;
        } else if (strcmp(argv[i], "--regions") == 0) {
            show_regions = 1;
        } else if (strcmp(argv[i], "--dump-colored") == 0 && i + 1 < argc) {
            dump_colored = argv[++i];
        } else if (strcmp(argv[i], "--fix-regions") == 0 && i + 1 < argc) {
            fix_min = (size_t)strtoul(argv[++i], NULL, 10);
        } else {
            fprintf(stderr, "unknown/incomplete option: %s\n", argv[i]);
            usage();
            return 1;
        }
    }

    arena = Arena_new();
    if (ObjIO_read(arena, in_path, &verts, &nv, &faces, &nf) != 0) {
        fprintf(stderr, "obj_reorient: failed to read %s\n", in_path);
        Arena_dispose(&arena);
        return 1;
    }
    fprintf(stderr, "[obj_reorient] read %s: %zu verts, %zu faces\n", in_path, nv, nf);
    fprintf(stderr, "[obj_reorient] axis=(%.0f,%.0f,%.0f)+t(%.0f,%.0f,%.0f) radius=%.1f "
            "cut=%d orient_first=%d\n",
            (double)axp[0], (double)axp[1], (double)axp[2],
            (double)axd[0], (double)axd[1], (double)axd[2],
            (double)radius, do_cut, orient_first);

    reorient_run(arena, &verts, &nv, faces, &nf, axp, axd, radius,
                 orient_first, do_cut, fix_min, &st);
    print_stats(&st);

    if (show_regions || dump_colored != NULL) {
        double L2 = sqrt((double)axd[0] * axd[0] + (double)axd[1] * axd[1]
                         + (double)axd[2] * axd[2]);
        if (L2 > 1e-9) {
            double adn[3] = { axd[0] / L2, axd[1] / L2, axd[2] / L2 };
            int R0m = mesh_radial_sign(verts, faces, nf, axp, adn);
            radial_regions(verts, nv, faces, nf, axp, adn, R0m, 24, dump_colored);
        }
    }

    if (st.nm_after > 0 || st.nmv_after > st.nmv_before)
        fprintf(stderr, "[obj_reorient] WARNING: NON-MANIFOLD (nm_edge=%zu nm_vert=%zu)\n",
                st.nm_after, st.nmv_after);
    if (st.same_dir_after > 0)
        fprintf(stderr, "[obj_reorient] WARNING: %zu same_dir edge(s) remain "
                "(raise --radius or inspect)\n", st.same_dir_after);

    if (dry_run) {
        fprintf(stderr, "[obj_reorient] dry-run: not writing output\n");
    } else if (ObjIO_write(out_path, verts, nv, faces, nf) != 0) {
        fprintf(stderr, "obj_reorient: failed to write %s\n", out_path);
        Arena_dispose(&arena);
        return 1;
    } else {
        fprintf(stderr, "[obj_reorient] wrote %s: %zu verts, %zu faces\n",
                out_path, nv, nf);
    }

    Arena_dispose(&arena);
    return 0;
}
