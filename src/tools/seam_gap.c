/* seam_gap.c -- Census of boundary-loop "gaps" in a welded mesh, localized to
 * cube-boundary seam planes.
 *
 * A weld that "misses obvious triangles" leaves the mesh OPEN at the seam: a
 * ring of boundary edges (edges owned by exactly one face) that should have been
 * closed by the BPA bridge / pinhole / hole-fill. The smallest such defect is a
 * single-triangle gap -- a 3-edge boundary loop, literally one missing triangle.
 *
 * grid_weld's own report only prints a single `unpaired` (= total boundary edge)
 * count, which lumps genuine seam holes together with the entire OUTER perimeter
 * of the block. This tool instead TRACES the boundary edges into loops, measures
 * each loop's size + location, and reports how many small gaps sit ON a seam
 * plane vs. elsewhere -- the number the weld audit actually needs.
 *
 * Definitions:
 *   boundary edge   -- undirected edge incident to exactly 1 face.
 *   non-manifold    -- undirected edge incident to >=3 faces (reported, not traced).
 *   loop            -- a closed chain of boundary edges (traced via the single
 *                      incident face's winding, so the hole interior is consistent).
 *   single-tri gap  -- a loop of exactly 3 edges (== 3 verts): one missing triangle.
 *   small gap       -- a loop of 4..SMALL_MAX edges (default 8).
 *   near-seam       -- every vertex of the loop is within `band` vox of ONE detected
 *                      seam plane (a multiple of `cube`). Isolates seam holes from
 *                      the block's outer trim perimeter (which runs far from any
 *                      interior seam plane).
 *
 * Usage:
 *   seam_gap <mesh.obj> [--cube <size=128>] [--band <vox=6>] [--small-max <8>]
 *            [--plane <axis> <coord>]... [--out <gaps.obj>] [--verbose]
 *   seam_gap --selftest
 *
 * axis: 0=Z, 1=Y, 2=X (OBJ vertex order in this pipeline is Z Y X).
 * If no --plane is given, seam planes are auto-detected (mirrors seam_weld.c).
 *
 * stdout ends with a machine-parseable line:
 *   SEAMGAP_JSON: { ... }
 * so a runner can log the census without re-parsing the human report.
 *
 * Exit: 0 = no near-seam gaps (clean) OR selftest pass,
 *       1 = near-seam gaps found OR selftest fail,
 *       2 = usage / IO error.
 *
 * Build: standalone (seam_gap.vcxproj compiles only this file).
 */
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SMALL_MAX_DEFAULT 8
#define MAX_PLANES 64

/* ===================================================================
 * Minimal OBJ reader (self-contained; identical convention to seam_audit.c).
 * First three floats of each "v" line; 0-based indices of each triangular "f"
 * line (tolerates v/vt/vn slashes).
 * =================================================================== */
typedef struct { float *v; size_t nv; int32_t *f; size_t nf; } Mesh;

static int obj_read(const char *path, Mesh *m)
{
    FILE *fp = fopen(path, "r");
    if (!fp) return -1;
    size_t vcap = 1024, fcap = 1024;
    m->v = (float *)malloc(vcap * 3 * sizeof(float));
    m->f = (int32_t *)malloc(fcap * 3 * sizeof(int32_t));
    m->nv = 0; m->nf = 0;
    char line[512];
    while (fgets(line, sizeof line, fp)) {
        if (line[0] == 'v' && line[1] == ' ') {
            if (m->nv >= vcap) { vcap *= 2; m->v = (float *)realloc(m->v, vcap * 3 * sizeof(float)); }
            double a = 0, b = 0, c = 0;
            sscanf(line + 2, "%lf %lf %lf", &a, &b, &c);
            m->v[m->nv * 3 + 0] = (float)a; m->v[m->nv * 3 + 1] = (float)b; m->v[m->nv * 3 + 2] = (float)c;
            m->nv++;
        } else if (line[0] == 'f' && line[1] == ' ') {
            int ia = 0, ib = 0, ic = 0;
            if (sscanf(line + 2, "%d/%*d/%*d %d/%*d/%*d %d/%*d/%*d", &ia, &ib, &ic) == 3 ||
                sscanf(line + 2, "%d//%*d %d//%*d %d//%*d", &ia, &ib, &ic) == 3 ||
                sscanf(line + 2, "%d/%*d %d/%*d %d/%*d", &ia, &ib, &ic) == 3 ||
                sscanf(line + 2, "%d %d %d", &ia, &ib, &ic) == 3) {
                if (m->nf >= fcap) { fcap *= 2; m->f = (int32_t *)realloc(m->f, fcap * 3 * sizeof(int32_t)); }
                m->f[m->nf * 3 + 0] = ia - 1; m->f[m->nf * 3 + 1] = ib - 1; m->f[m->nf * 3 + 2] = ic - 1;
                m->nf++;
            }
        }
    }
    fclose(fp);
    return 0;
}
static void mesh_free(Mesh *m) { free(m->v); free(m->f); m->v = NULL; m->f = NULL; }

/* ===================================================================
 * Seam-plane detection (mirrors seam_weld.c::detect_planes / seam_audit.c).
 * A plane at coord c on axis a is real only if the mesh has material on BOTH
 * sides (so the outer bounding faces of the block are not mistaken for seams).
 * =================================================================== */
#define SEAM_MIN_SIDE 2.0
typedef struct { int axis; double coord; } Plane;

static size_t detect_planes(const Mesh *m, double cube, double band, Plane *planes, size_t cap)
{
    double mn[3] = { 1e30, 1e30, 1e30 }, mx[3] = { -1e30, -1e30, -1e30 };
    for (size_t v = 0; v < m->nv; v++)
        for (int a = 0; a < 3; a++) { double c = m->v[v * 3 + (size_t)a]; if (c < mn[a]) mn[a] = c; if (c > mx[a]) mx[a] = c; }
    size_t np = 0;
    for (int a = 0; a < 3 && np < cap; a++) {
        long k0 = (long)ceil(mn[a] / cube), k1 = (long)floor(mx[a] / cube);
        for (long k = k0; k <= k1 && np < cap; k++) {
            double c = (double)k * cube; int lo = 0, hi = 0;
            for (size_t v = 0; v < m->nv; v++) {
                double cc = m->v[v * 3 + (size_t)a];
                if (cc > c - band && cc < c - SEAM_MIN_SIDE) lo = 1;
                if (cc > c + SEAM_MIN_SIDE && cc < c + band) hi = 1;
                if (lo && hi) break;
            }
            if (lo && hi) { planes[np].axis = a; planes[np].coord = c; np++; }
        }
    }
    return np;
}

/* Distance of coord to the nearest detected plane on the same axis; also reports
 * whether the point is within `band` of that plane. Returns 1e30 if no plane on
 * that axis. */
static double nearest_plane_dist(const Plane *planes, size_t np, int axis, double coord)
{
    double best = 1e30;
    for (size_t i = 0; i < np; i++) {
        if (planes[i].axis != axis) continue;
        double d = fabs(coord - planes[i].coord);
        if (d < best) best = d;
    }
    return best;
}

/* ===================================================================
 * Directed edge list -> undirected multiplicity -> boundary edges.
 * A boundary edge is an undirected edge that appears in exactly one face; its
 * single directed instance (x->y in that face's winding) is the trace direction.
 * =================================================================== */
typedef struct { int32_t s, d; } DEdge;   /* directed src->dst */

static int cmp_undirected(const void *pa, const void *pb)
{
    const DEdge *a = (const DEdge *)pa, *b = (const DEdge *)pb;
    int32_t a0 = a->s < a->d ? a->s : a->d, a1 = a->s < a->d ? a->d : a->s;
    int32_t b0 = b->s < b->d ? b->s : b->d, b1 = b->s < b->d ? b->d : b->s;
    if (a0 != b0) return a0 < b0 ? -1 : 1;
    if (a1 != b1) return a1 < b1 ? -1 : 1;
    return 0;
}

/* Sort directed boundary edges by src ascending, so a vertex's outgoing edges
 * are contiguous and binary-searchable. */
static int cmp_dedge_src(const void *pa, const void *pb)
{
    const DEdge *a = (const DEdge *)pa, *b = (const DEdge *)pb;
    return a->s < b->s ? -1 : (a->s > b->s ? 1 : 0);
}

/* ===================================================================
 * A traced boundary loop.
 * =================================================================== */
#define LOOP_VID_CAP 64
typedef struct {
    size_t nedges;          /* == number of verts in a simple cycle */
    double cen[3];          /* centroid of the loop's verts */
    double mn[3], mx[3];    /* loop bbox */
    double extent;          /* bbox diagonal (proxy for hole diameter) */
    int    near_seam;       /* every vert within band of ONE plane */
    int    seam_axis;       /* which axis' plane it hugs (-1 if not near-seam) */
    double seam_coord;      /* that plane's coord */
    int    straddles;       /* has verts on both sides of that plane */
    int    pinched;         /* walk passed a vertex with >1 outgoing boundary edge */
    int32_t vids[LOOP_VID_CAP]; /* loop vertex indices in walk order (small loops only) */
    int    nvids;           /* min(nedges, LOOP_VID_CAP) */
} Loop;

typedef struct {
    size_t boundary_edges;
    size_t nm_edges;        /* undirected edges with >=3 faces */
    Loop  *loops;
    size_t nloops;
    size_t pinch_starts;    /* boundary verts with >1 outgoing boundary edge */
} GapReport;

/* Trace boundary loops. next-pointer walk: from the current directed boundary
 * edge x->y, continue with an unvisited boundary edge leaving y. A clean hole
 * has exactly one such edge per vertex; a pinch (figure-8) has several -- we
 * take any unvisited one and count the pinch, which splits the figure-8 into
 * sub-loops (correct for a size census). */
static void trace_loops(const Mesh *m, const Plane *planes, size_t np, double band,
                        GapReport *rep)
{
    rep->boundary_edges = 0; rep->nm_edges = 0; rep->pinch_starts = 0;
    rep->loops = NULL; rep->nloops = 0;
    if (m->nf == 0) return;

    size_t hn = m->nf * 3;
    DEdge *he = (DEdge *)malloc(hn * sizeof(DEdge));
    for (size_t f = 0; f < m->nf; f++) {
        int32_t a = m->f[f * 3 + 0], b = m->f[f * 3 + 1], c = m->f[f * 3 + 2];
        he[f * 3 + 0].s = a; he[f * 3 + 0].d = b;
        he[f * 3 + 1].s = b; he[f * 3 + 1].d = c;
        he[f * 3 + 2].s = c; he[f * 3 + 2].d = a;
    }
    qsort(he, hn, sizeof(DEdge), cmp_undirected);

    /* Collect the boundary directed edges (undirected multiplicity == 1). */
    DEdge *bnd = (DEdge *)malloc(hn * sizeof(DEdge));
    size_t nb = 0;
    for (size_t i = 0; i < hn;) {
        size_t j = i + 1;
        int32_t a0 = he[i].s < he[i].d ? he[i].s : he[i].d;
        int32_t a1 = he[i].s < he[i].d ? he[i].d : he[i].s;
        while (j < hn) {
            int32_t b0 = he[j].s < he[j].d ? he[j].s : he[j].d;
            int32_t b1 = he[j].s < he[j].d ? he[j].d : he[j].s;
            if (b0 != a0 || b1 != a1) break;
            j++;
        }
        size_t run = j - i;
        if (run == 1) { bnd[nb++] = he[i]; }
        else if (run > 2) { rep->nm_edges++; }
        i = j;
    }
    rep->boundary_edges = nb;
    free(he);
    if (nb == 0) { free(bnd); return; }

    /* Sort the boundary edges by src so a vertex's outgoing edges are contiguous
     * and binary-searchable. `cur` below indexes this sorted array; used[] tracks
     * consumed edges by that same position. */
    qsort(bnd, nb, sizeof(DEdge), cmp_dedge_src);
    int32_t *src_sorted = (int32_t *)malloc(nb * sizeof(int32_t));
    for (size_t i = 0; i < nb; i++) src_sorted[i] = bnd[i].s;

    /* Per-vertex boundary degree = #boundary edges touching it. A simple loop
     * vertex has degree 2; degree != 2 means a pinch (the boundary passes
     * through the vertex more than once) -- exactly the test seam_hole_fill /
     * pinhole use to REJECT a loop as non-simple, so it flags loops the fillers
     * cannot close even though the tracer walks them as a cycle. */
    int32_t *bdeg = (int32_t *)calloc(m->nv, sizeof(int32_t));
    for (size_t i = 0; i < nb; i++) { bdeg[bnd[i].s]++; bdeg[bnd[i].d]++; }

    uint8_t *used = (uint8_t *)calloc(nb, 1);

    /* find range [lo,hi) in sorted bnd whose src == key */
    #define FIND_SRC(key, LO, HI) do {                                   \
        size_t _l = 0, _r = nb;                                          \
        while (_l < _r) { size_t _mch = (_l + _r) / 2;                   \
            if (src_sorted[_mch] < (key)) _l = _mch + 1; else _r = _mch; } \
        (LO) = _l;                                                       \
        size_t _l2 = 0, _r2 = nb;                                        \
        while (_l2 < _r2) { size_t _mch = (_l2 + _r2) / 2;               \
            if (src_sorted[_mch] <= (key)) _l2 = _mch + 1; else _r2 = _mch; } \
        (HI) = _l2;                                                      \
    } while (0)

    size_t loops_cap = 64;
    rep->loops = (Loop *)malloc(loops_cap * sizeof(Loop));
    rep->nloops = 0;

    for (size_t start = 0; start < nb; start++) {
        if (used[start]) continue;
        /* Walk from bnd[start]. */
        size_t cur = start;
        Loop lp; memset(&lp, 0, sizeof lp);
        lp.mn[0] = lp.mn[1] = lp.mn[2] = 1e30;
        lp.mx[0] = lp.mx[1] = lp.mx[2] = -1e30;
        lp.nvids = 0;
        double sum[3] = { 0, 0, 0 };
        size_t nedges = 0;
        int32_t loop_start_vert = bnd[start].s;
        int guard = 0;
        while (1) {
            if (used[cur]) break;
            used[cur] = 1;
            int32_t sv = bnd[cur].s, dv = bnd[cur].d;
            if (bdeg[sv] != 2) lp.pinched = 1;   /* pinch vertex (fillers reject) */
            /* accumulate the src vertex of this edge (each vert counted once
             * around a simple loop). */
            for (int k = 0; k < 3; k++) {
                double p = m->v[(size_t)sv * 3 + (size_t)k];
                sum[k] += p;
                if (p < lp.mn[k]) lp.mn[k] = p;
                if (p > lp.mx[k]) lp.mx[k] = p;
            }
            if (lp.nvids < LOOP_VID_CAP) lp.vids[lp.nvids++] = sv;
            nedges++;
            if (dv == loop_start_vert) break;   /* closed */
            /* find an unvisited outgoing boundary edge from dv */
            size_t lo, hi;
            FIND_SRC(dv, lo, hi);
            size_t next = (size_t)-1, avail = 0;
            for (size_t t = lo; t < hi; t++) {
                if (!used[t]) { if (next == (size_t)-1) next = t; avail++; }
            }
            if (avail > 1) { rep->pinch_starts++; lp.pinched = 1; }
            if (next == (size_t)-1) break;      /* dead end (open chain / consumed) */
            cur = next;
            if (++guard > (int)nb + 4) break;   /* safety */
        }
        if (nedges == 0) continue;
        lp.nedges = nedges;
        for (int k = 0; k < 3; k++) lp.cen[k] = sum[k] / (double)nedges;
        {
            double dz = lp.mx[0] - lp.mn[0], dy = lp.mx[1] - lp.mn[1], dx = lp.mx[2] - lp.mn[2];
            lp.extent = sqrt(dz * dz + dy * dy + dx * dx);
        }

        /* near-seam: every vert within band of ONE plane. We test with the
         * loop bbox: for a plane on axis a at coord c, all verts within band
         * <=> the bbox on axis a lies within [c-band, c+band]. */
        lp.near_seam = 0; lp.seam_axis = -1; lp.seam_coord = 0; lp.straddles = 0;
        for (size_t pi = 0; pi < np; pi++) {
            int a = planes[pi].axis; double c = planes[pi].coord;
            if (lp.mn[a] >= c - band && lp.mx[a] <= c + band) {
                lp.near_seam = 1; lp.seam_axis = a; lp.seam_coord = c;
                lp.straddles = (lp.mn[a] < c && lp.mx[a] > c);
                break;
            }
        }
        if (rep->nloops >= loops_cap) {
            loops_cap *= 2;
            rep->loops = (Loop *)realloc(rep->loops, loops_cap * sizeof(Loop));
        }
        rep->loops[rep->nloops++] = lp;
    }
    #undef FIND_SRC
    free(src_sorted); free(bnd); free(used); free(bdeg);
}

/* Optional: write a colored OBJ of the gap loops as line segments.
 *   red    = single-triangle gap on a seam
 *   orange = small gap (4..small_max) on a seam
 *   yellow = larger boundary loop on a seam
 *   grey   = boundary loop NOT near a seam (outer perimeter etc.)
 * Only near-seam loops (+ small non-seam gaps) are emitted; the full outer
 * perimeter is skipped unless --verbose scale is wanted. */
static void write_gap_obj(const char *path, const Mesh *m, const GapReport *rep,
                          size_t small_max)
{
    FILE *fp = fopen(path, "w");
    if (!fp) { fprintf(stderr, "seam_gap: cannot write %s\n", path); return; }
    fprintf(fp, "# seam_gap loops: red=1-tri seam, orange=small seam, yellow=big seam, grey=non-seam small\n");
    /* We need the loop verts again; re-trace is costly, so instead re-walk the
     * boundary here would duplicate logic. Simplify: emit each loop's bbox
     * diagonal as a marker plus its centroid point. That is enough to locate a
     * gap in a viewer; the welded OBJ itself shows the actual open ring. */
    size_t vidx = 1;
    for (size_t i = 0; i < rep->nloops; i++) {
        const Loop *lp = &rep->loops[i];
        int emit = lp->near_seam || lp->nedges <= small_max;
        if (!emit) continue;
        double r, g, b;
        if (lp->near_seam && lp->nedges == 3)        { r = 1.0; g = 0.0; b = 0.0; }
        else if (lp->near_seam && lp->nedges <= small_max) { r = 1.0; g = 0.55; b = 0.0; }
        else if (lp->near_seam)                      { r = 1.0; g = 1.0; b = 0.0; }
        else                                         { r = 0.6; g = 0.6; b = 0.6; }
        /* a small cross at the centroid */
        double cx = lp->cen[0], cy = lp->cen[1], cz = lp->cen[2];
        double s = 1.0;
        fprintf(fp, "v %.4f %.4f %.4f %.3f %.3f %.3f\n", cx - s, cy, cz, r, g, b);
        fprintf(fp, "v %.4f %.4f %.4f %.3f %.3f %.3f\n", cx + s, cy, cz, r, g, b);
        fprintf(fp, "v %.4f %.4f %.4f %.3f %.3f %.3f\n", cx, cy - s, cz, r, g, b);
        fprintf(fp, "v %.4f %.4f %.4f %.3f %.3f %.3f\n", cx, cy + s, cz, r, g, b);
        fprintf(fp, "v %.4f %.4f %.4f %.3f %.3f %.3f\n", cx, cy, cz - s, r, g, b);
        fprintf(fp, "v %.4f %.4f %.4f %.3f %.3f %.3f\n", cx, cy, cz + s, r, g, b);
        fprintf(fp, "l %zu %zu\n", vidx + 0, vidx + 1);
        fprintf(fp, "l %zu %zu\n", vidx + 2, vidx + 3);
        fprintf(fp, "l %zu %zu\n", vidx + 4, vidx + 5);
        vidx += 6;
    }
    fclose(fp);
    (void)m;
}

/* ===================================================================
 * Report
 * =================================================================== */
static int report(const Mesh *m, const GapReport *rep, size_t small_max,
                  const Plane *planes, size_t np, int verbose, int dump_loops)
{
    size_t tri_seam = 0, tri_other = 0;
    size_t small_seam = 0, small_other = 0;
    size_t big_seam = 0;
    size_t seam_loops = 0;
    size_t max_loop = 0;
    size_t straddle_tri = 0;
    /* Merger-safe fill candidates: a near-seam loop that straddles the plane,
     * is small (<= small_max edges), and whose extent (bbox diagonal) is under
     * the inter-wrap clearance (7 vox) cannot join two different wraps -- both
     * its sides are the SAME sheet crossing the ~1-vox seam. These are exactly
     * the holes a relaxed seam-hole fill could safely close. */
    const double SAFE_EXTENT = 6.0;
    size_t safe_fill = 0;
    double ext_max_seam_small = 0.0;

    for (size_t i = 0; i < rep->nloops; i++) {
        const Loop *lp = &rep->loops[i];
        if (lp->nedges > max_loop) max_loop = lp->nedges;
        if (lp->near_seam) seam_loops++;
        if (lp->nedges == 3) {
            if (lp->near_seam) { tri_seam++; if (lp->straddles) straddle_tri++; }
            else tri_other++;
        } else if (lp->nedges <= small_max) {
            if (lp->near_seam) small_seam++; else small_other++;
        } else if (lp->near_seam) {
            big_seam++;
        }
        if (lp->near_seam && lp->nedges <= small_max) {
            if (lp->extent > ext_max_seam_small) ext_max_seam_small = lp->extent;
            if (lp->straddles && lp->extent <= SAFE_EXTENT) safe_fill++;
        }
    }

    printf("=== seam_gap ===\n");
    printf("  V=%zu F=%zu   seam planes=%zu   band-classified\n", m->nv, m->nf, np);
    for (size_t i = 0; i < np; i++)
        printf("    plane[%zu] axis=%d (%s) coord=%.1f\n", i, planes[i].axis,
               planes[i].axis == 0 ? "Z" : (planes[i].axis == 1 ? "Y" : "X"),
               planes[i].coord);
    printf("  boundary edges=%zu  non-manifold edges=%zu  loops=%zu  pinch-steps=%zu\n",
           rep->boundary_edges, rep->nm_edges, rep->nloops, rep->pinch_starts);
    printf("  SEAM gaps:  single-triangle=%zu (straddling=%zu)  small(4..%zu)=%zu  larger=%zu\n",
           tri_seam, straddle_tri, small_max, small_seam, big_seam);
    printf("  non-seam:   single-triangle=%zu  small=%zu  (larger loops incl. outer perimeter omitted)\n",
           tri_other, small_other);
    printf("  merger-safe fill candidates (straddling small, extent<=%.1f vox): %zu"
           "   (max seam-small extent seen: %.2f vox)\n",
           SAFE_EXTENT, safe_fill, ext_max_seam_small);

    if (verbose) {
        printf("  --- near-seam gap loops (<= %zu edges) ---\n", small_max);
        for (size_t i = 0; i < rep->nloops; i++) {
            const Loop *lp = &rep->loops[i];
            if (!lp->near_seam || lp->nedges > small_max) continue;
            printf("    len=%-2zu  cen(z,y,x)=(%.1f,%.1f,%.1f)  axis=%s coord=%.1f  straddle=%d  extent=%.2f\n",
                   lp->nedges, lp->cen[0], lp->cen[1], lp->cen[2],
                   lp->seam_axis == 0 ? "Z" : (lp->seam_axis == 1 ? "Y" : "X"),
                   lp->seam_coord, lp->straddles, lp->extent);
        }
    }

    if (dump_loops) {
        printf("  --- near-seam small-loop vertices (idx: z y x | side) ---\n");
        for (size_t i = 0; i < rep->nloops; i++) {
            const Loop *lp = &rep->loops[i];
            if (!lp->near_seam || lp->nedges > small_max) continue;
            printf("  loop %zu: len=%zu axis=%s coord=%.1f extent=%.2f straddle=%d %s\n",
                   i, lp->nedges,
                   lp->seam_axis == 0 ? "Z" : (lp->seam_axis == 1 ? "Y" : "X"),
                   lp->seam_coord, lp->extent, lp->straddles,
                   lp->pinched ? "PINCHED" : "simple");
            for (int k = 0; k < lp->nvids; k++) {
                int32_t v = lp->vids[k];
                double c = m->v[(size_t)v * 3 + (size_t)lp->seam_axis];
                double sd = c - lp->seam_coord;
                printf("      v%-7d  %.3f %.3f %.3f   %s%.2f\n",
                       v, (double)m->v[(size_t)v*3+0], (double)m->v[(size_t)v*3+1],
                       (double)m->v[(size_t)v*3+2], sd >= 0 ? "+" : "", sd);
            }
        }
    }

    /* machine-parseable summary */
    printf("SEAMGAP_JSON: {\"V\":%zu,\"F\":%zu,\"boundary_edges\":%zu,\"nm_edges\":%zu,"
           "\"loops\":%zu,\"seam_loops\":%zu,\"tri_gaps_seam\":%zu,\"tri_straddle\":%zu,"
           "\"tri_gaps_other\":%zu,\"small_gaps_seam\":%zu,\"small_gaps_other\":%zu,"
           "\"big_seam_loops\":%zu,\"max_loop\":%zu,\"pinch_steps\":%zu,\"planes\":%zu,"
           "\"safe_fill\":%zu,\"ext_max_seam_small\":%.2f}\n",
           m->nv, m->nf, rep->boundary_edges, rep->nm_edges, rep->nloops, seam_loops,
           tri_seam, straddle_tri, tri_other, small_seam, small_other, big_seam,
           max_loop, rep->pinch_starts, np, safe_fill, ext_max_seam_small);

    return (tri_seam + small_seam + big_seam) > 0 ? 1 : 0;
}

/* ===================================================================
 * Selftest -- synthetic meshes with known boundary topology.
 * =================================================================== */
static void tetra(Mesh *m, int drop_face)   /* closed tetra, or open if drop_face>=0 */
{
    static float V[12] = { 0,0,0,  10,0,0,  0,10,0,  0,0,10 };
    static int32_t Fall[12] = { 0,1,2,  0,3,1,  0,2,3,  1,3,2 };
    m->v = (float *)malloc(sizeof V); memcpy(m->v, V, sizeof V); m->nv = 4;
    m->f = (int32_t *)malloc(sizeof Fall);
    m->nf = 0;
    for (int fi = 0; fi < 4; fi++) {
        if (fi == drop_face) continue;
        m->f[m->nf * 3 + 0] = Fall[fi * 3 + 0];
        m->f[m->nf * 3 + 1] = Fall[fi * 3 + 1];
        m->f[m->nf * 3 + 2] = Fall[fi * 3 + 2];
        m->nf++;
    }
}

static int selftest(void)
{
    int fails = 0;
    Plane planes[MAX_PLANES]; size_t np;
    GapReport rep;

    /* 1. Closed tetra: 0 boundary edges, 0 loops. */
    {
        Mesh m; tetra(&m, -1);
        np = 0;   /* no planes */
        trace_loops(&m, planes, np, 6.0, &rep);
        if (rep.boundary_edges != 0 || rep.nloops != 0) {
            fprintf(stderr, "[selftest] closed tetra: be=%zu loops=%zu -> FAIL\n",
                    rep.boundary_edges, rep.nloops); fails++;
        } else fprintf(stderr, "[selftest] closed tetra -> ok (watertight)\n");
        free(rep.loops); mesh_free(&m);
    }

    /* 2. Open tetra (drop face 3): one triangular hole = loop of 3 edges. */
    {
        Mesh m; tetra(&m, 3);
        np = 0;
        trace_loops(&m, planes, np, 6.0, &rep);
        int ok = (rep.nloops == 1 && rep.loops && rep.loops[0].nedges == 3 &&
                  rep.boundary_edges == 3);
        if (!ok) {
            fprintf(stderr, "[selftest] open tetra: loops=%zu be=%zu len0=%zu -> FAIL\n",
                    rep.nloops, rep.boundary_edges, rep.nloops ? rep.loops[0].nedges : 0);
            fails++;
        } else fprintf(stderr, "[selftest] open tetra -> ok (1 tri hole)\n");
        free(rep.loops); mesh_free(&m);
    }

    /* 3. Single triangle: its 3 edges are all boundary -> 1 loop len 3. */
    {
        Mesh m;
        static float V[9] = { 0,0,0, 10,0,0, 0,10,0 };
        static int32_t F[3] = { 0,1,2 };
        m.v = (float *)malloc(sizeof V); memcpy(m.v, V, sizeof V); m.nv = 3;
        m.f = (int32_t *)malloc(sizeof F); memcpy(m.f, F, sizeof F); m.nf = 1;
        np = 0;
        trace_loops(&m, planes, np, 6.0, &rep);
        if (!(rep.nloops == 1 && rep.loops[0].nedges == 3)) {
            fprintf(stderr, "[selftest] single tri: loops=%zu -> FAIL\n", rep.nloops); fails++;
        } else fprintf(stderr, "[selftest] single tri -> ok (len-3 loop)\n");
        free(rep.loops); mesh_free(&m);
    }

    /* 4. Two tris sharing an edge (a quad): shared edge interior, 4 boundary
     *    edges -> 1 loop of len 4 (a small gap). */
    {
        Mesh m;
        static float V[12] = { 0,0,0, 10,0,0, 10,10,0, 0,10,0 };
        static int32_t F[6] = { 0,1,2,  0,2,3 };
        m.v = (float *)malloc(sizeof V); memcpy(m.v, V, sizeof V); m.nv = 4;
        m.f = (int32_t *)malloc(sizeof F); memcpy(m.f, F, sizeof F); m.nf = 2;
        np = 0;
        trace_loops(&m, planes, np, 6.0, &rep);
        if (!(rep.nloops == 1 && rep.loops[0].nedges == 4 && rep.boundary_edges == 4)) {
            fprintf(stderr, "[selftest] quad: loops=%zu be=%zu len=%zu -> FAIL\n",
                    rep.nloops, rep.boundary_edges, rep.nloops ? rep.loops[0].nedges : 0);
            fails++;
        } else fprintf(stderr, "[selftest] quad -> ok (len-4 small gap)\n");
        free(rep.loops); mesh_free(&m);
    }

    /* 5. Non-manifold edge: 3 triangles fanning a shared edge (0-1). The edge
     *    0-1 has 3 incident faces -> 1 nm edge. */
    {
        Mesh m;
        static float V[15] = { 0,0,0, 10,0,0, 0,10,0, 0,-10,0, 5,0,10 };
        static int32_t F[9] = { 0,1,2,  0,1,3,  0,1,4 };
        m.v = (float *)malloc(sizeof V); memcpy(m.v, V, sizeof V); m.nv = 5;
        m.f = (int32_t *)malloc(sizeof F); memcpy(m.f, F, sizeof F); m.nf = 3;
        np = 0;
        trace_loops(&m, planes, np, 6.0, &rep);
        if (rep.nm_edges != 1) {
            fprintf(stderr, "[selftest] nm fan: nm_edges=%zu -> FAIL\n", rep.nm_edges); fails++;
        } else fprintf(stderr, "[selftest] nm fan -> ok (1 non-manifold edge)\n");
        free(rep.loops); mesh_free(&m);
    }

    /* 6. Seam localization: open tetra hole placed straddling plane X=0 within
     *    band. Its verts span x in [-5,5] -> near-seam + straddles. */
    {
        Mesh m; tetra(&m, 3);
        /* shift so the dropped-face hole (verts 1,3,2 = (10,0,0),(0,0,10),(0,10,0))
         * straddles X (axis 2). Vertex order is Z Y X, so axis 2 is the 3rd coord.
         * Move x-coords to straddle 0: subtract 5 from the 3rd coordinate. */
        for (size_t v = 0; v < m.nv; v++) m.v[v * 3 + 2] -= 5.0f;
        planes[0].axis = 2; planes[0].coord = 0.0; np = 1;
        trace_loops(&m, planes, np, 20.0, &rep);
        int ok = (rep.nloops == 1 && rep.loops[0].nedges == 3 &&
                  rep.loops[0].near_seam == 1);
        if (!ok) {
            fprintf(stderr, "[selftest] seam tri: loops=%zu near=%d -> FAIL\n",
                    rep.nloops, rep.nloops ? rep.loops[0].near_seam : -1); fails++;
        } else fprintf(stderr, "[selftest] seam tri -> ok (near-seam len-3)\n");
        free(rep.loops); mesh_free(&m);
    }

    /* 7. detect_planes: two stacked unit sheets straddling X=128 find that plane. */
    {
        Mesh m;
        /* two 2-tri sheets, one at x in [125,127], one at [129,131]; the gap in
         * between makes X=128 a real seam. Vertex order Z,Y,X. */
        static float V[24] = {
            0,0,125,  0,10,125,  10,0,127,   10,10,127,   /* low side */
            0,0,129,  0,10,129,  10,0,131,   10,10,131    /* high side */
        };
        static int32_t F[12] = { 0,1,2, 1,3,2,  4,5,6, 5,7,6 };
        m.v = (float *)malloc(sizeof V); memcpy(m.v, V, sizeof V); m.nv = 8;
        m.f = (int32_t *)malloc(sizeof F); memcpy(m.f, F, sizeof F); m.nf = 4;
        np = detect_planes(&m, 128.0, 6.0, planes, MAX_PLANES);
        int found = 0;
        for (size_t i = 0; i < np; i++) if (planes[i].axis == 2 && fabs(planes[i].coord - 128.0) < 1e-6) found = 1;
        if (!found) {
            fprintf(stderr, "[selftest] detect_planes: np=%zu, X=128 not found -> FAIL\n", np); fails++;
        } else fprintf(stderr, "[selftest] detect_planes -> ok (X=128 seam)\n");
        mesh_free(&m);
    }

    fprintf(stderr, "=== seam_gap selftest %s (%d failure%s) ===\n",
            fails ? "FAILED" : "PASSED", fails, fails == 1 ? "" : "s");
    return fails ? 1 : 0;
}

/* ===================================================================
 * Fold / overlap detector.
 * An interior (2-face) edge whose two faces fold back onto each other -- normals
 * nearly opposite -- is a flap: the two triangles overlap in space (z-fighting),
 * a defect the manifold/winding audits are blind to (it is edge-manifold and can
 * be consistently wound). Reports near-seam folds with location + dihedral.
 * =================================================================== */
typedef struct { int32_t lo, hi, f; } EF;
static int cmp_ef(const void *pa, const void *pb)
{
    const EF *a = (const EF *)pa, *b = (const EF *)pb;
    if (a->lo != b->lo) return a->lo < b->lo ? -1 : 1;
    if (a->hi != b->hi) return a->hi < b->hi ? -1 : 1;
    return 0;
}
static void face_norm(const Mesh *m, int32_t f, double n[3])
{
    int32_t a = m->f[f*3+0], b = m->f[f*3+1], c = m->f[f*3+2];
    double e1[3] = { m->v[b*3+0]-m->v[a*3+0], m->v[b*3+1]-m->v[a*3+1], m->v[b*3+2]-m->v[a*3+2] };
    double e2[3] = { m->v[c*3+0]-m->v[a*3+0], m->v[c*3+1]-m->v[a*3+1], m->v[c*3+2]-m->v[a*3+2] };
    n[0] = e1[1]*e2[2]-e1[2]*e2[1]; n[1] = e1[2]*e2[0]-e1[0]*e2[2]; n[2] = e1[0]*e2[1]-e1[1]*e2[0];
    double L = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
    if (L > 1e-20) { n[0]/=L; n[1]/=L; n[2]/=L; }
}
static int detect_folds(const Mesh *m, const Plane *planes, size_t np, double band,
                        double dot_thresh, int near_seam_only)
{
    size_t hn = m->nf * 3;
    EF *ef = (EF *)malloc(hn * sizeof(EF));
    for (size_t f = 0; f < m->nf; f++) for (int e = 0; e < 3; e++) {
        int32_t a = m->f[f*3+e], b = m->f[f*3+(e+1)%3];
        ef[f*3+e].lo = a<b?a:b; ef[f*3+e].hi = a<b?b:a; ef[f*3+e].f = (int32_t)f;
    }
    qsort(ef, hn, sizeof(EF), cmp_ef);
    int nfold = 0, printed = 0;
    for (size_t i = 0; i < hn;) {
        size_t j = i + 1;
        while (j < hn && ef[j].lo==ef[i].lo && ef[j].hi==ef[i].hi) j++;
        if (j - i == 2) {
            double n1[3], n2[3]; face_norm(m, ef[i].f, n1); face_norm(m, ef[i+1].f, n2);
            double d = n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];
            if (d < dot_thresh) {
                double mid[3];
                for (int k = 0; k < 3; k++) mid[k] = 0.5*(m->v[(size_t)ef[i].lo*3+k]+m->v[(size_t)ef[i].hi*3+k]);
                int near = 0;
                for (size_t p = 0; p < np; p++)
                    if (fabs(mid[planes[p].axis] - planes[p].coord) <= band) { near = 1; break; }
                if (near_seam_only && !near) { i = j; continue; }
                nfold++;
                if (printed < 20) {
                    printf("  FOLD dot=%+.2f  edge (z,y,x)=(%.1f,%.1f,%.1f)  verts %d-%d%s\n",
                           d, mid[0], mid[1], mid[2], ef[i].lo, ef[i].hi, near ? "  [near-seam]" : "");
                    printed++;
                }
            }
        }
        i = j;
    }
    printf("  fold/overlap edges (normal dot < %.2f%s): %d\n",
           dot_thresh, near_seam_only ? ", near-seam" : "", nfold);
    free(ef);
    return nfold;
}

/* ===================================================================
 * Main
 * =================================================================== */
int main(int argc, char **argv)
{
    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0) return selftest();
    if (argc < 2) {
        fprintf(stderr,
            "Usage: %s <mesh.obj> [--cube <128>] [--band <6>] [--small-max <8>]\n"
            "       [--plane <axis> <coord>]... [--out <gaps.obj>] [--verbose] [--dump-loops]\n"
            "       %s --selftest\n"
            "  axis: 0=Z 1=Y 2=X. If no --plane given, seam planes are auto-detected.\n"
            "  Prints a census of boundary-loop gaps; ends with SEAMGAP_JSON: {...}.\n",
            argv[0], argv[0]);
        return 2;
    }
    const char *path = argv[1];
    double cube = 128.0, band = 6.0;
    size_t small_max = SMALL_MAX_DEFAULT;
    const char *out_path = NULL;
    int verbose = 0;
    int dump_loops = 0;
    int folds = 0; double folds_thresh = -0.5; int folds_near_only = 0;
    int crop_set = 0; double crop_c[3] = {0,0,0}, crop_r = 0;  /* --crop z y x r */
    Plane user_planes[MAX_PLANES]; size_t n_user = 0;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--cube") && i + 1 < argc) cube = atof(argv[++i]);
        else if (!strcmp(argv[i], "--band") && i + 1 < argc) band = atof(argv[++i]);
        else if (!strcmp(argv[i], "--small-max") && i + 1 < argc) small_max = (size_t)atol(argv[++i]);
        else if (!strcmp(argv[i], "--out") && i + 1 < argc) out_path = argv[++i];
        else if (!strcmp(argv[i], "--verbose")) verbose = 1;
        else if (!strcmp(argv[i], "--dump-loops")) dump_loops = 1;
        else if (!strcmp(argv[i], "--folds") || !strcmp(argv[i], "--folds-near")) {
            folds = 1; folds_near_only = (strcmp(argv[i], "--folds-near") == 0);
            /* optional numeric threshold (may be negative, e.g. -0.8) */
            if (i+1 < argc) { char c0 = argv[i+1][0], c1 = argv[i+1][1];
                if ((c0>='0'&&c0<='9') || c0=='.' ||
                    (c0=='-' && (c1=='.' || (c1>='0'&&c1<='9'))))
                    folds_thresh = atof(argv[++i]); }
        }
        else if (!strcmp(argv[i], "--crop") && i + 4 < argc) {
            crop_set = 1; crop_c[0]=atof(argv[++i]); crop_c[1]=atof(argv[++i]);
            crop_c[2]=atof(argv[++i]); crop_r=atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "--plane") && i + 2 < argc) {
            if (n_user < MAX_PLANES) {
                user_planes[n_user].axis = atoi(argv[++i]);
                user_planes[n_user].coord = atof(argv[++i]);
                n_user++;
            }
        } else {
            fprintf(stderr, "seam_gap: unknown arg %s\n", argv[i]);
            return 2;
        }
    }

    Mesh m;
    if (obj_read(path, &m) != 0) {
        fprintf(stderr, "seam_gap: cannot read %s\n", path);
        return 2;
    }

    /* --crop: write ONLY the faces with a vertex within crop_r of crop_c, remapped,
     * to out_path -- a small local neighbourhood for close-up rendering. Boundary
     * edges of a gap show as the open ring; the crop keeps it visible. */
    if (crop_set) {
        if (!out_path) { fprintf(stderr, "seam_gap: --crop needs --out\n"); mesh_free(&m); return 2; }
        int32_t *remap = (int32_t *)malloc(m.nv * sizeof(int32_t));
        for (size_t v = 0; v < m.nv; v++) remap[v] = -1;
        double r2 = crop_r * crop_r;
        FILE *fp = fopen(out_path, "w");
        if (!fp) { fprintf(stderr, "seam_gap: cannot write %s\n", out_path); free(remap); mesh_free(&m); return 2; }
        size_t out_nv = 0, out_nf = 0;
        /* first pass: mark verts of kept faces */
        uint8_t *keep = (uint8_t *)calloc(m.nf, 1);
        for (size_t f = 0; f < m.nf; f++) {
            int in = 0;
            for (int k = 0; k < 3; k++) {
                int32_t vi = m.f[f*3+k];
                double dz = m.v[(size_t)vi*3+0]-crop_c[0];
                double dy = m.v[(size_t)vi*3+1]-crop_c[1];
                double dx = m.v[(size_t)vi*3+2]-crop_c[2];
                if (dz*dz+dy*dy+dx*dx <= r2) { in = 1; break; }
            }
            if (in) { keep[f]=1; for (int k=0;k<3;k++){ int32_t vi=m.f[f*3+k]; if(remap[vi]<0) remap[vi]=(int32_t)out_nv++; } }
        }
        for (size_t v = 0; v < m.nv; v++) if (remap[v] >= 0)
            fprintf(fp, "v %.4f %.4f %.4f\n", (double)m.v[v*3+0], (double)m.v[v*3+1], (double)m.v[v*3+2]);
        for (size_t f = 0; f < m.nf; f++) if (keep[f]) {
            fprintf(fp, "f %d %d %d\n", remap[m.f[f*3+0]]+1, remap[m.f[f*3+1]]+1, remap[m.f[f*3+2]]+1);
            out_nf++;
        }
        fclose(fp);
        fprintf(stderr, "seam_gap: crop -> %s (%zu verts, %zu faces within %.1f vox)\n",
                out_path, out_nv, out_nf, crop_r);
        free(keep); free(remap); mesh_free(&m);
        return 0;
    }

    Plane planes[MAX_PLANES]; size_t np;
    if (n_user > 0) {
        for (size_t i = 0; i < n_user; i++) planes[i] = user_planes[i];
        np = n_user;
    } else {
        np = detect_planes(&m, cube, band > 6.0 ? band : 6.0, planes, MAX_PLANES);
    }

    if (folds) {
        printf("=== fold/overlap edges ===\n");
        detect_folds(&m, planes, np, band, folds_thresh, folds_near_only);
    }

    GapReport rep;
    trace_loops(&m, planes, np, band, &rep);
    int has_gap = report(&m, &rep, small_max, planes, np, verbose, dump_loops);
    if (out_path) write_gap_obj(out_path, &m, &rep, small_max);

    free(rep.loops);
    mesh_free(&m);
    return has_gap;
}
