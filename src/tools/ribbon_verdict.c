/* ribbon_verdict.c -- mechanical acceptance gate for a fitted quad ribbon.
 *
 * WHY THIS EXISTS (2026-09-01).  The 4x5x5 "state of the art"
 * (output/pherc0139_4x5x5_qfit_nopack_20260901) shipped with a changelog
 * claim that its cross-sections "show one monotone blue->red spiral with
 * every wrap tracked separately".  They do not: the ribbon draws 37,820
 * faces whose edges jump between physical wraps, and the claim was made from
 * summary statistics that could not see them.  Every number this tool
 * reports is one a human would otherwise have to take on trust, and each is
 * a property of the DELIVERABLE, never of an intermediate fill percentage.
 *
 * The gates, in the order a reader notices the defect:
 *
 *   TANGLE        no emitted edge exceeds the physical wrap gate, and no
 *                 emitted edge steps a wrap radially.  These are the radial
 *                 "spokes" and the right-angle staircase jogs.
 *   WINDING       no emitted face joins two different reconstruction lanes,
 *                 and no lattice cell is covered twice.  A duplicate cover is
 *                 one wrap drawn on top of another.
 *   DISTORTION    the 3-D length a quad edge represents per unit of UV, which
 *                 must not vary wildly inside one triangle.  Reported as
 *                 anisotropy quantiles plus the fraction beyond 4x, matching
 *                 the Sander threshold the bakes are already judged on.
 *   COVERAGE      occupied lattice cells against the atlas rectangle, and the
 *                 measured (not interpolated) share of them.
 *   FRAGMENTATION component count and the largest component's share of the
 *                 occupied cells.  A sheet nobody can follow is not a sheet.
 *
 * Connectivity is the MESH's own face adjacency, never a lattice-neighbour
 * rule -- the same doctrine as ribbon_largest_strip, so the two agree on what
 * a component is.
 *
 *   ribbon_verdict <ribbon.vmesh> [--labels recon.i32] [--support support.u8]
 *                  [--umb-y F --umb-x F] [--pitch F=9.5]
 *                  [--du F=2] [--dv F=1] [--report out.json]
 *   ribbon_verdict --selftest
 *
 * Exit: 0 all gates pass / selftest pass, 1 IO, 2 usage, 3 selftest fail,
 *       4 a gate FAILED (the deliverable is not shippable).
 */
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/arena.h"
#include "../common/mesh_bin.h"
#include "../common/union_find.h"

/* ---- gate thresholds ---------------------------------------------------
 * Every value is a measured property of the shipped 4x5x5 lineage, not a
 * preference.  They are REGIME-BOUND (CLAUDE.md rule 26): revalidate when the
 * rung's radial band or pitch changes. */

/* Max 3-D step between adjacent ribbon lattice cells.  Mirrors RIB_WRAP_GATE
 * in src/flatten/ribbon.c; below the 7-vox inter-wrap clearance, so an edge
 * under it cannot span two wraps. */
#define RV_WRAP_GATE          6.0

/* An emitted edge crosses a wrap only if its RADIAL step reaches the
 * inter-wrap clearance.  With RV_WRAP_GATE at 6.0 vox and a 9.5-vox pitch that
 * is impossible by construction, which is the point: this gate is the
 * REGIME check on the wrap gate (CLAUDE.md rule 26).  On a rung whose pitch is
 * small enough that a 6-vox step does span a wrap, it fires.  Steeply radial
 * SHORT edges are ordinary sheet geometry -- a fold, a turn -- and half a
 * pitch wrongly flagged 139 of them on the 4x5x5. */
#define RV_WRAP_DR_FRACTION   0.7

/* Sander-equivalent, in the project's own language: the bakes report ">= 4x
 * stretch" and the shipped lineage runs 0.00%-0.04% there.
 *
 * The quantity must be SCALE relative to the sheet's own median, not the ratio
 * between a triangle's own edges: u is angle x a constant here, not arclength,
 * so the u- and v-direction scales differ everywhere by construction.  That is
 * anisotropy the metric campaign owns, not warping a reader would see.  What
 * damages readability is a patch whose scale departs from the rest of the
 * sheet, which is what area-per-UV-area against the median measures. */
#define RV_ANISO_LIMIT        4.0
#define RV_ANISO_FRAC_MAX     0.005

/* Component COUNT is not a defect measure.  Measured 2026-09-01 on the 4x5x5:
 * the ribbon's 21 occupied column runs cover steadily larger radii (79-303,
 * 176-341, 202-367, ... 404-457) and get steadily shorter, which is the box
 * CROP -- an outer wrap leaves the 512x640x640 grid part-way round, so only
 * part of each turn exists.  Joining two such bands would fuse material ~16
 * wraps apart, exactly the inter-wrap merger the changelog closed.  The
 * baseline's 58.3%-in-one-piece was manufactured by cross-wrap spokes.
 *
 * What IS a defect is an AVOIDABLE break: two components that sit within the
 * wrap gate in 3-D *and* within a couple of lattice columns in u, i.e. the
 * emitter had every reason to join them and did not. */
#define RV_CONTACT_U_COLS     2.0   /* |du| in lattice columns for a contact */

/* The measure that IS actionable: inside one contiguous run of occupied
 * lattice columns the data has no crop gap, so the sheet there should be one
 * piece.  Comparing the largest component against the largest such RUN, rather
 * than against the whole atlas, asks "is the ribbon coherent where the data is
 * continuous?" and is immune to how much of a turn the grid box happens to
 * contain. */
#define RV_RUN_COHERENCE_MIN  0.75

/* Degenerate UV edges carry no metric information; below this they leave the
 * anisotropy population rather than dividing by ~0. */
#define RV_UV_EPS             1e-6

#define RV_PITCH_DEFAULT      9.5

/* Umbilicus core exclusion.  The shipped view already drops r < 64 (the
 * wrap-separation campaign's verdict: turns cram as the circumference goes to
 * zero, the prediction there is a fused blob at ~8 extractable verts per 1k
 * foreground voxels, and in-place recovery was falsified for StrokeStrip,
 * SLIM and quilting alike).  Measuring the deliverable means measuring what
 * ships, so the gates honour the same mask -- and the excluded fraction is
 * always reported so the exclusion can never hide a regression. */
#define RV_CORE_RADIUS_DEFAULT 0.0

/* Coverage against the SOURCE mesh.  Every other gate here asks whether what
 * the ribbon drew is correct; none asks whether it drew everything it should.
 * Looking at the 4x5x5 cross-sections showed the ribbon dropping the innermost
 * one or two wraps and a band at the bottom -- 14% of the source curve length
 * -- while all six gates passed.  A source vertex counts as covered when a
 * ribbon vertex sits within half the wrap gate of it. */
/* Depth-seam (plateau) detection.  The 2026-09-01 red-box review: blocks of
 * claims sit 1-2 vox deeper or shallower than their surroundings with sharp
 * lattice-aligned edges -- the claim DP hopping between the FRONT and BACK
 * face of the ~1.2-vox prediction shell (occasionally a full wrap).  Every
 * hop is under the wrap gate and invisible to the other gates, but the bake
 * samples a different fibre layer of the papyrus and paints rectangles of
 * misregistered texture.  Detection: high-pass the per-cell radius with a
 * local, face-connected-component median (removes legitimate slope -- the raw
 * shell band flags 21%% of edges, the high-passed field 2.3%%), mark emitted
 * axis-aligned mesh-edge steps over RV_DEPTH_STEP, and count runs of
 * RV_DEPTH_RUN_MIN or more.  Merely adjacent atlas charts are not a surface
 * seam: a provenance cut must therefore remove the edge from this population.
 * Champion baseline at introduction: 5,991 seam runs, 0.874%% seam-edge
 * share -- this gate is BORN FAILING by design (metric-first development);
 * the claim-depth fix is judged by driving it under RV_DEPTH_SEAM_MAX. */
#define RV_DEPTH_STEP        0.75   /* vox: high-passed step per lattice edge */
#define RV_DEPTH_RUN_MIN     4      /* cells: minimum axis-aligned seam run */
#define RV_DEPTH_MEDIAN_WIN  9      /* cells: local median window (odd) */
#define RV_DEPTH_SEAM_MAX    0.001  /* max seam-edge share of lattice edges */
/* The depth sweep holds six arrays over the whole (u,v) lattice, 14 bytes a cell.  Past this
 * many cells it is skipped -- and until 2026-09-17 it was skipped SILENTLY, so the 21x5x5
 * (2,767 x 95,328 = 2.64e8 cells, just over the cap) reported its eighth gate as absent for
 * want of an umbilicus it actually had.  The reason is now carried to the gate line. */
#define RV_DEPTH_MAX_CELLS   2.5e8

#define RV_COVER_RADIUS       3.0
#define RV_COVER_MIN          0.85

typedef struct {
    /* tangle */
    size_t edges_total;
    size_t edges_long;          /* > RV_WRAP_GATE in 3-D */
    size_t edges_cross_wrap;    /* radial step >= RV_WRAP_DR_FRACTION*pitch */
    double edge_len_max;
    double edge_dr_max;
    /* winding */
    size_t faces_cross_lane;
    size_t faces_cross_material;
    size_t cells_duplicated;
    int    have_labels;
    int    have_materials;
    int    have_axis;
    int    have_support;
    /* distortion */
    size_t core_faces_excluded;
    double core_radius;
    double aniso_p50, aniso_p90, aniso_p99, aniso_max;
    double aniso_frac_over;
    size_t aniso_samples;
    /* coverage */
    size_t cells_occupied;
    size_t cells_rect;
    size_t verts_supported;
    double u_lo, u_hi, v_lo, v_hi;
    /* fragmentation */
    size_t n_components;
    size_t largest_cells;
    size_t avoidable_break_pairs;   /* component pairs in 3-D + u contact */
    size_t avoidable_break_contacts;/* vertex pairs supporting them */
    size_t source_verts;            /* source mesh vertices, if supplied */
    size_t source_covered;          /* ... within RV_COVER_RADIUS of a ribbon vertex */
    int    have_source;
    size_t column_runs;             /* maximal runs of occupied columns */
    size_t best_run_cells;          /* cells in the largest run x row band */
    size_t best_run_largest_comp;   /* the largest component inside it */
    size_t row_bands;               /* unjoinable row bands in the largest run */
    /* depth seams */
    size_t depth_seam_runs;
    size_t depth_seam_edges;
    size_t depth_edges_total;
    int    depth_measured;
    char   depth_skip[160];         /* why it was not measured, when it was not */
} RvReport;

/* ---- small helpers ----------------------------------------------------- */

static double rv_d3(const float *a, const float *b)
{
    double d0 = (double)a[0] - (double)b[0];
    double d1 = (double)a[1] - (double)b[1];
    double d2 = (double)a[2] - (double)b[2];
    return sqrt(d0*d0 + d1*d1 + d2*d2);
}

static int rv_cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    if (x < y) return -1;
    if (x > y) return 1;
    return 0;
}

static double rv_quantile(const double *sorted, size_t n, double q)
{
    double pos = 0.0;
    size_t idx = 0;
    if (n == 0) return 0.0;
    pos = q * (double)(n - 1);
    idx = (size_t)(pos + 0.5);
    if (idx >= n) idx = n - 1;
    return sorted[idx];
}

/* Radius about the scroll axis.  Ribbon vertices are (z,y,x); the axis is the
 * z direction through (umb_y, umb_x), which is how every other tool in this
 * lane addresses PHerc0139. */
static double rv_radius(const float *p, double umb_y, double umb_x)
{
    double dy = (double)p[1] - umb_y;
    double dx = (double)p[2] - umb_x;
    return sqrt(dy*dy + dx*dx);
}

/* Emit one row per edge in a retained depth-seam run.  This is deliberately
 * optional: the acceptance report stays compact, while a failed Stage 3 can
 * be localized without reimplementing the topology-aware metric in Python.
 * `u_edge` means the edge joins neighbouring U columns and the run extends in
 * V; otherwise it joins neighbouring V rows and extends in U. */
static void rv_write_depth_run(FILE *fp, size_t run_id, int u_edge,
                               size_t fixed, size_t start, size_t len,
                               size_t wcols, long lo_col, long lo_row,
                               double du, double dv, const float *V,
                               const float *UV, const float *cr,
                               const float *hp, const int32_t *cell_vertex,
                               const int32_t *cell_component,
                               const int32_t *labels,
                               const int32_t *materials,
                               const unsigned char *support)
{
    size_t k = 0;
    if (fp == NULL || cell_vertex == NULL) return;
    for (k = 0; k < len; k++) {
        size_t rr = u_edge ? start + k : fixed;
        size_t cc = u_edge ? fixed : start + k;
        size_t a = rr * wcols + cc;
        size_t b = u_edge ? a + 1 : a + wcols;
        int32_t va = cell_vertex[a], vb = cell_vertex[b];
        int32_t la = va >= 0 && labels != NULL ? labels[va] : -1;
        int32_t lb = vb >= 0 && labels != NULL ? labels[vb] : -1;
        int32_t ma = va >= 0 && materials != NULL ? materials[va] : -1;
        int32_t mb = vb >= 0 && materials != NULL ? materials[vb] : -1;
        int sa = va >= 0 && support != NULL ? (int)support[va] : -1;
        int sb = vb >= 0 && support != NULL ? (int)support[vb] : -1;
        double ua = va >= 0 ? (double)UV[(size_t)va*2+0]
                            : (double)(lo_col + (long)cc) * du;
        double va_uv = va >= 0 ? (double)UV[(size_t)va*2+1]
                               : (double)(lo_row + (long)rr) * dv;
        double ub = vb >= 0 ? (double)UV[(size_t)vb*2+0]
                            : (double)(lo_col + (long)cc + (u_edge ? 1 : 0))*du;
        double vb_uv = vb >= 0 ? (double)UV[(size_t)vb*2+1]
                               : (double)(lo_row + (long)rr + (u_edge ? 0 : 1))*dv;
        const float *pa = va >= 0 ? &V[(size_t)va*3] : NULL;
        const float *pb = vb >= 0 ? &V[(size_t)vb*3] : NULL;
        fprintf(fp,
                "%zu,%s,%zu,%zu,%d,%d,%d,%d,%d,%d,%d,%d,%d,"
                "%.6f,%.6f,%.6f,%.6f,"
                "%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,"
                "%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                run_id, u_edge ? "U" : "V", len, k,
                va, vb, cell_component[a], la, lb, ma, mb, sa, sb,
                ua, va_uv, ub, vb_uv,
                pa ? (double)pa[0] : 0.0, pa ? (double)pa[1] : 0.0,
                pa ? (double)pa[2] : 0.0, pb ? (double)pb[0] : 0.0,
                pb ? (double)pb[1] : 0.0, pb ? (double)pb[2] : 0.0,
                (double)cr[a], (double)cr[b], (double)hp[a], (double)hp[b],
                fabs((double)hp[a] - (double)hp[b]),
                fabs((double)cr[a] - (double)cr[b]));
    }
}

/* One undirected edge of the emitted triangulation, canonicalized. */
typedef struct { int32_t a, b; } RvEdge;

static int rv_cmp_edge(const void *pa, const void *pb)
{
    const RvEdge *x = (const RvEdge *)pa;
    const RvEdge *y = (const RvEdge *)pb;
    if (x->a != y->a) return x->a < y->a ? -1 : 1;
    if (x->b != y->b) return x->b < y->b ? -1 : 1;
    return 0;
}

/* A lattice cell, optionally tagged with the component that owns it. */
typedef struct { int32_t comp, col, row; } RvCell;

static int rv_cmp_cell(const void *pa, const void *pb)
{
    const RvCell *x = (const RvCell *)pa;
    const RvCell *y = (const RvCell *)pb;
    if (x->comp != y->comp) return x->comp < y->comp ? -1 : 1;
    if (x->col  != y->col)  return x->col  < y->col  ? -1 : 1;
    if (x->row  != y->row)  return x->row  < y->row  ? -1 : 1;
    return 0;
}

/* ---- the measurement --------------------------------------------------- */

/* `labels` and `support` may be NULL; the corresponding gates then report as
 * ABSENT rather than silently passing -- an unmeasured gate must never look
 * like a satisfied one. */
static int rv_measure(Arena_T arena, const MeshBinData *mesh,
                      const int32_t *labels, const int32_t *materials,
                      const unsigned char *support, const MeshBinData *source,
                      int have_axis, double umb_y, double umb_x, double pitch,
                      double core_radius, double du, double dv,
                      FILE *depth_csv, RvReport *out)
{
    const float *V = mesh->verts;
    const float *UV = mesh->uv;
    size_t nv = mesh->nv, nf = mesh->nf;
    RvEdge *edges = NULL;
    RvCell *cells = NULL;
    double *aniso = NULL;
    UnionFind uf;
    int32_t *dense = NULL, *label = NULL;
    size_t *comp_cells = NULL;
    size_t ne = 0, na = 0, i = 0, f = 0, ncomp = 0;
    double dr_gate = RV_WRAP_DR_FRACTION * pitch;

    memset(out, 0, sizeof(*out));
    out->have_labels = labels != NULL;
    out->have_materials = materials != NULL;
    out->have_support = support != NULL;
    out->have_axis = have_axis;
    out->core_radius = core_radius;
    if (nv == 0 || nv > (size_t)INT32_MAX || nf == 0) return -1;
    if (UV == NULL || du <= 0.0 || dv <= 0.0) return -1;
    /* the arena request is a 32-bit long; a pathological ribbon (a fit that
     * emitted its whole lattice) must be refused with a reason, not an OOM
     * from a wrapped size */
    if (nf > ((size_t)LONG_MAX) / (3u * sizeof(RvEdge))) {
        fprintf(stderr, "ERROR: %zu faces exceed the verdict allocation "
                        "budget (max %zu); the ribbon itself is pathological",
                nf, ((size_t)LONG_MAX) / (3u * sizeof(RvEdge)));
        fputc(10, stderr);
        return -1;
    }

    /* ---- unique emitted edges, and per-triangle anisotropy ---- */
    edges = (RvEdge *)ARENA_ALLOC(arena, (size_t)(nf * 3 * sizeof(RvEdge)));
    aniso = (double *)ARENA_ALLOC(arena, (size_t)(nf * sizeof(double)));
    for (f = 0; f < nf; f++) {
        int32_t t[3];
        double ratio[3];
        int e = 0, nratio = 0;
        for (e = 0; e < 3; e++) t[e] = mesh->faces[f*3 + (size_t)e];
        for (e = 0; e < 3; e++) {
            int32_t a = t[e], b = t[(e + 1) % 3];
            RvEdge *dst = &edges[ne++];
            double l3 = rv_d3(&V[(size_t)a*3], &V[(size_t)b*3]);
            double dU = (double)UV[(size_t)a*2+0] - (double)UV[(size_t)b*2+0];
            double dVv = (double)UV[(size_t)a*2+1] - (double)UV[(size_t)b*2+1];
            double luv = sqrt(dU*dU + dVv*dVv);
            dst->a = a < b ? a : b;
            dst->b = a < b ? b : a;
            if (luv > RV_UV_EPS) ratio[nratio++] = l3 / luv;
        }
        /* Scale of this triangle: 3-D area per unit of UV area.  Compared
         * against the sheet's median below, so a uniform global stretch costs
         * nothing and only a locally warped patch is flagged. */
        (void)nratio;
        {
            const float *p0 = &V[(size_t)t[0]*3];
            const float *p1 = &V[(size_t)t[1]*3];
            const float *p2 = &V[(size_t)t[2]*3];
            double e1[3], e2[3], cx, cy, cz, a3, auv;
            int d = 0;
            for (d = 0; d < 3; d++) {
                e1[d] = (double)p1[d] - (double)p0[d];
                e2[d] = (double)p2[d] - (double)p0[d];
            }
            cx = e1[1]*e2[2] - e1[2]*e2[1];
            cy = e1[2]*e2[0] - e1[0]*e2[2];
            cz = e1[0]*e2[1] - e1[1]*e2[0];
            a3 = 0.5 * sqrt(cx*cx + cy*cy + cz*cz);
            {
                double u0 = (double)UV[(size_t)t[0]*2+0];
                double v0 = (double)UV[(size_t)t[0]*2+1];
                double u1 = (double)UV[(size_t)t[1]*2+0];
                double v1 = (double)UV[(size_t)t[1]*2+1];
                double u2 = (double)UV[(size_t)t[2]*2+0];
                double v2 = (double)UV[(size_t)t[2]*2+1];
                auv = 0.5 * fabs((u1-u0)*(v2-v0) - (u2-u0)*(v1-v0));
            }
            if (auv > RV_UV_EPS && a3 > 0.0) {
                /* Skip the umbilicus core when a mask is in force. */
                int in_core = 0;
                if (have_axis && core_radius > 0.0) {
                    double rr = (rv_radius(p0, umb_y, umb_x) +
                                 rv_radius(p1, umb_y, umb_x) +
                                 rv_radius(p2, umb_y, umb_x)) / 3.0;
                    in_core = rr < core_radius;
                }
                if (in_core) out->core_faces_excluded++;
                else aniso[na++] = a3 / auv;
            }
        }
    }
    qsort(edges, ne, sizeof(*edges), rv_cmp_edge);
    {
        size_t k = 0;
        for (i = 0; i < ne; i++) {
            if (i > 0 && edges[i].a == edges[i-1].a &&
                edges[i].b == edges[i-1].b) continue;
            edges[k++] = edges[i];
        }
        ne = k;
    }
    out->edges_total = ne;
    for (i = 0; i < ne; i++) {
        const float *pa = &V[(size_t)edges[i].a * 3];
        const float *pb = &V[(size_t)edges[i].b * 3];
        double len = rv_d3(pa, pb);
        if (len > out->edge_len_max) out->edge_len_max = len;
        if (len > RV_WRAP_GATE) out->edges_long++;
        if (have_axis) {
            double dr = fabs(rv_radius(pa, umb_y, umb_x) -
                             rv_radius(pb, umb_y, umb_x));
            if (dr > out->edge_dr_max) out->edge_dr_max = dr;
            if (dr >= dr_gate) out->edges_cross_wrap++;
        }
    }

    if (na > 0) {
        size_t over = 0;
        double med = 0.0;
        qsort(aniso, na, sizeof(*aniso), rv_cmp_double);
        med = rv_quantile(aniso, na, 0.50);
        if (!(med > RV_UV_EPS)) med = 1.0;
        /* Normalize by the sheet's own scale: a globally stretched but
         * uniform parameterization is perfectly readable. */
        for (i = 0; i < na; i++) aniso[i] /= med;
        out->aniso_samples = na;
        out->aniso_p50 = rv_quantile(aniso, na, 0.50);
        out->aniso_p90 = rv_quantile(aniso, na, 0.90);
        out->aniso_p99 = rv_quantile(aniso, na, 0.99);
        out->aniso_max = aniso[na - 1];
        for (i = 0; i < na; i++)
            if (aniso[i] >= RV_ANISO_LIMIT || aniso[i] <= 1.0/RV_ANISO_LIMIT)
                over++;
        out->aniso_frac_over = (double)over / (double)na;
    }

    /* ---- faces that join two reconstruction lanes ---- */
    if (labels != NULL) {
        for (f = 0; f < nf; f++) {
            int32_t l0 = labels[mesh->faces[f*3+0]];
            int32_t l1 = labels[mesh->faces[f*3+1]];
            int32_t l2 = labels[mesh->faces[f*3+2]];
            /* -1 is "unresolved derived cell", not a distinct lane: an
             * interpolated vertex inherits no identity and must not be
             * reported as a lane crossing. */
            if (l0 < 0 || l1 < 0 || l2 < 0) continue;
            if (l0 != l1 || l1 != l2) out->faces_cross_lane++;
        }
    }
    if (materials != NULL) {
        for (f = 0; f < nf; f++) {
            int32_t m0 = materials[mesh->faces[f*3+0]];
            int32_t m1 = materials[mesh->faces[f*3+1]];
            int32_t m2 = materials[mesh->faces[f*3+2]];
            if (m0 < 0 || m1 < 0 || m2 < 0) continue;
            if (m0 != m1 || m1 != m2) out->faces_cross_material++;
        }
    }
    if (support != NULL)
        for (i = 0; i < nv; i++) if (support[i] != 0) out->verts_supported++;

    /* ---- components by the mesh's own face adjacency ---- */
    uf = UF_new(arena, (int32_t)nv);
    for (f = 0; f < nf; f++) {
        uf_union(&uf, mesh->faces[f*3+0], mesh->faces[f*3+1]);
        uf_union(&uf, mesh->faces[f*3+1], mesh->faces[f*3+2]);
    }
    dense = (int32_t *)ARENA_ALLOC(arena, (size_t)(nv * sizeof(int32_t)));
    label = (int32_t *)ARENA_ALLOC(arena, (size_t)(nv * sizeof(int32_t)));
    for (i = 0; i < nv; i++) dense[i] = -1;
    for (i = 0; i < nv; i++) {
        int32_t r = uf_find(&uf, (int32_t)i);
        if (dense[r] < 0) dense[r] = (int32_t)ncomp++;
        label[i] = dense[r];
    }
    out->n_components = ncomp;

    /* ---- lattice cells: occupancy, duplicates, per-component area ----
     * One keyed array serves both questions.  Ranking components by DISTINCT
     * cells (not vertices) is the ribbon_largest_strip rule: a finely
     * tessellated sliver must not outrank a coarse long strip. */
    cells = (RvCell *)ARENA_ALLOC(arena, (size_t)(nv * sizeof(RvCell)));
    comp_cells = (size_t *)ARENA_ALLOC(arena, (size_t)(ncomp * sizeof(size_t)));
    for (i = 0; i < ncomp; i++) comp_cells[i] = 0;
    out->u_lo = out->v_lo = 1e300;
    out->u_hi = out->v_hi = -1e300;
    for (i = 0; i < nv; i++) {
        double u = (double)UV[i*2+0], v = (double)UV[i*2+1];
        if (u < out->u_lo) out->u_lo = u;
        if (u > out->u_hi) out->u_hi = u;
        if (v < out->v_lo) out->v_lo = v;
        if (v > out->v_hi) out->v_hi = v;
        cells[i].comp = label[i];
        cells[i].col = (int32_t)lround(u / du);
        cells[i].row = (int32_t)lround(v / dv);
    }
    qsort(cells, nv, sizeof(*cells), rv_cmp_cell);
    for (i = 0; i < nv; i++) {
        if (i > 0 && cells[i].comp == cells[i-1].comp &&
            cells[i].col == cells[i-1].col && cells[i].row == cells[i-1].row)
            continue;
        comp_cells[cells[i].comp]++;
    }
    for (i = 0; i < ncomp; i++) {
        out->cells_occupied += comp_cells[i];
        if (comp_cells[i] > out->largest_cells)
            out->largest_cells = comp_cells[i];
    }
    /* A duplicated cell is one covered by two vertices that the mesh does not
     * join -- i.e. two sheets stacked at the same UV.  Re-key without the
     * component so a cross-component stack is what gets counted. */
    for (i = 0; i < nv; i++) cells[i].comp = 0;
    qsort(cells, nv, sizeof(*cells), rv_cmp_cell);
    {
        size_t distinct = 0;
        for (i = 0; i < nv; i++) {
            if (i > 0 && cells[i].col == cells[i-1].col &&
                cells[i].row == cells[i-1].row) continue;
            distinct++;
        }
        out->cells_duplicated = nv - distinct;
    }
    /* ---- avoidable breaks -------------------------------------------------
     * A uniform spatial hash at the wrap gate; a bucket plus its 26 neighbours
     * hold every candidate.  A contact is two vertices of DIFFERENT components
     * within the gate in 3-D and within RV_CONTACT_U_COLS lattice columns in
     * u: adjacent in space and in parameter, so nothing but a cut kept them
     * apart.  Pairs are counted once via a sorted (min,max) pass. */
    {
        const size_t nb = 1u << 20;
        int32_t *bhead = (int32_t *)ARENA_ALLOC(arena, (size_t)(nb*sizeof(int32_t)));
        int32_t *bnext = (int32_t *)ARENA_ALLOC(arena, (size_t)(nv*sizeof(int32_t)));
        RvEdge *pair = (RvEdge *)ARENA_ALLOC(arena, (size_t)(nv*sizeof(RvEdge)));
        size_t npair = 0;
        double cellsz = RV_WRAP_GATE;
        double ulim = RV_CONTACT_U_COLS * du;
        double blo[3] = { 1e300, 1e300, 1e300 };
        for (i = 0; i < nb; i++) bhead[i] = -1;
        for (i = 0; i < nv; i++)
            for (f = 0; f < 3; f++)
                if ((double)V[i*3+f] < blo[f]) blo[f] = (double)V[i*3+f];
        for (i = 0; i < nv; i++) {
            long long g0 = (long long)floor(((double)V[i*3+0]-blo[0])/cellsz);
            long long g1 = (long long)floor(((double)V[i*3+1]-blo[1])/cellsz);
            long long g2 = (long long)floor(((double)V[i*3+2]-blo[2])/cellsz);
            unsigned long long h = (unsigned long long)g0 * 0x9E3779B97F4A7C15ull
                                 ^ (unsigned long long)g1 * 0xBF58476D1CE4E5B9ull
                                 ^ (unsigned long long)g2 * 0x94D049BB133111EBull;
            size_t b;
            h ^= h >> 31;
            b = (size_t)h & (nb - 1);
            bnext[i] = bhead[b];
            bhead[b] = (int32_t)i;
        }
        for (i = 0; i < nv && npair < nv; i++) {
            long long g0 = (long long)floor(((double)V[i*3+0]-blo[0])/cellsz);
            long long g1 = (long long)floor(((double)V[i*3+1]-blo[1])/cellsz);
            long long g2 = (long long)floor(((double)V[i*3+2]-blo[2])/cellsz);
            int dz, dy, dx;
            for (dz = -1; dz <= 1; dz++)
            for (dy = -1; dy <= 1; dy++)
            for (dx = -1; dx <= 1; dx++) {
                unsigned long long h =
                      (unsigned long long)(g0+dz) * 0x9E3779B97F4A7C15ull
                    ^ (unsigned long long)(g1+dy) * 0xBF58476D1CE4E5B9ull
                    ^ (unsigned long long)(g2+dx) * 0x94D049BB133111EBull;
                int32_t it;
                int scanned = 0;
                size_t b;
                h ^= h >> 31;
                b = (size_t)h & (nb - 1);
                for (it = bhead[b]; it >= 0 && scanned < 64;
                     it = bnext[it], scanned++) {
                    size_t j = (size_t)it;
                    if (j <= i) continue;
                    if (label[i] == label[j]) continue;
                    if (rv_d3(&V[i*3], &V[j*3]) > RV_WRAP_GATE) continue;
                    if (fabs((double)UV[i*2+0] - (double)UV[j*2+0]) > ulim)
                        continue;
                    /* A break is only AVOIDABLE if a triangle could span
                     * it.  A mesh cannot join two vertices without a third,
                     * so look for one inside the gate of both; without it the
                     * two are honestly unjoinable and the break belongs to
                     * the data, not to the emitter. */
                    {
                        int has_third = 0;
                        int dz2, dy2, dx2;
                        for (dz2 = -1; dz2 <= 1 && !has_third; dz2++)
                        for (dy2 = -1; dy2 <= 1 && !has_third; dy2++)
                        for (dx2 = -1; dx2 <= 1 && !has_third; dx2++) {
                            unsigned long long h2 =
                                  (unsigned long long)(g0+dz2) * 0x9E3779B97F4A7C15ull
                                ^ (unsigned long long)(g1+dy2) * 0xBF58476D1CE4E5B9ull
                                ^ (unsigned long long)(g2+dx2) * 0x94D049BB133111EBull;
                            int32_t it2;
                            int sc2 = 0;
                            size_t b2;
                            h2 ^= h2 >> 31;
                            b2 = (size_t)h2 & (nb - 1);
                            for (it2 = bhead[b2]; it2 >= 0 && sc2 < 64;
                                 it2 = bnext[it2], sc2++) {
                                size_t m = (size_t)it2;
                                if (m == i || m == j) continue;
                                if (rv_d3(&V[m*3], &V[i*3]) > RV_WRAP_GATE)
                                    continue;
                                if (rv_d3(&V[m*3], &V[j*3]) > RV_WRAP_GATE)
                                    continue;
                                has_third = 1;
                                break;
                            }
                        }
                        if (!has_third) continue;
                    }
                    if (npair >= nv) break;
                    pair[npair].a = label[i] < label[j] ? label[i] : label[j];
                    pair[npair].b = label[i] < label[j] ? label[j] : label[i];
                    npair++;
                }
            }
        }
        /* Coverage: how much of the SOURCE the ribbon actually drew.  The
         * bucket grid is already loaded with every ribbon vertex, so this is
         * one extra pass over the source. */
        if (source != NULL && source->nv > 0) {
            size_t sv = 0;
            out->have_source = 1;
            out->source_verts = source->nv;
            for (sv = 0; sv < source->nv; sv++) {
                const float *sp = &source->verts[sv*3];
                /* The shipped view drops the umbilicus core, so coverage is
                 * reported against the same domain the deliverable uses. */
                if (have_axis && core_radius > 0.0 &&
                    rv_radius(sp, umb_y, umb_x) < core_radius) {
                    out->source_verts--;
                    continue;
                }
                long long g0 = (long long)floor(((double)sp[0]-blo[0])/cellsz);
                long long g1 = (long long)floor(((double)sp[1]-blo[1])/cellsz);
                long long g2 = (long long)floor(((double)sp[2]-blo[2])/cellsz);
                int dz, dy, dx, hit = 0;
                for (dz = -1; dz <= 1 && !hit; dz++)
                for (dy = -1; dy <= 1 && !hit; dy++)
                for (dx = -1; dx <= 1 && !hit; dx++) {
                    unsigned long long h =
                          (unsigned long long)(g0+dz) * 0x9E3779B97F4A7C15ull
                        ^ (unsigned long long)(g1+dy) * 0xBF58476D1CE4E5B9ull
                        ^ (unsigned long long)(g2+dx) * 0x94D049BB133111EBull;
                    int32_t it;
                    int sc = 0;
                    size_t b;
                    h ^= h >> 31;
                    b = (size_t)h & (nb - 1);
                    for (it = bhead[b]; it >= 0 && sc < 96;
                         it = bnext[it], sc++) {
                        if (rv_d3(&V[(size_t)it*3], sp) <= RV_COVER_RADIUS) {
                            hit = 1;
                            break;
                        }
                    }
                }
                if (hit) out->source_covered++;
            }
        }
        out->avoidable_break_contacts = npair;
        if (npair > 0) {
            qsort(pair, npair, sizeof(*pair), rv_cmp_edge);
            for (i = 0; i < npair; i++)
                if (i == 0 || pair[i].a != pair[i-1].a ||
                    pair[i].b != pair[i-1].b)
                    out->avoidable_break_pairs++;
        }
    }

    /* ---- coherence inside one contiguous column run ---------------------- */
    {
        long lo_col = (long)lround(out->u_lo / du);
        long hi_col = (long)lround(out->u_hi / du);
        size_t ncol = (size_t)(hi_col - lo_col + 1);
        unsigned char *occ = NULL;
        size_t *run_of_col = NULL;
        size_t *run_cells = NULL;
        size_t nrun = 0, best = 0, c = 0;
        if (ncol > 0 && ncol < (size_t)1 << 28) {
            occ = (unsigned char *)ARENA_ALLOC(arena, (size_t)ncol);
            run_of_col = (size_t *)ARENA_ALLOC(arena, (size_t)(ncol*sizeof(size_t)));
            memset(occ, 0, ncol);
            for (i = 0; i < nv; i++) {
                long cc = (long)lround((double)UV[i*2+0] / du) - lo_col;
                if (cc >= 0 && (size_t)cc < ncol) occ[cc] = 1;
            }
            for (c = 0; c < ncol; c++) {
                if (!occ[c]) { run_of_col[c] = (size_t)-1; continue; }
                if (c == 0 || !occ[c-1]) nrun++;
                run_of_col[c] = nrun - 1;
            }
            out->column_runs = nrun;
            if (nrun > 0) {
                run_cells = (size_t *)ARENA_ALLOC(arena, (size_t)(nrun*sizeof(size_t)));
                for (i = 0; i < nrun; i++) run_cells[i] = 0;
                for (i = 0; i < nv; i++) {
                    long cc = (long)lround((double)UV[i*2+0] / du) - lo_col;
                    if (cc >= 0 && (size_t)cc < ncol &&
                        run_of_col[cc] != (size_t)-1)
                        run_cells[run_of_col[cc]]++;
                }
                for (i = 0; i < nrun; i++)
                    if (run_cells[i] > run_cells[best]) best = i;
                out->best_run_cells = run_cells[best];
                {
                    size_t *inrun = (size_t *)ARENA_ALLOC(
                        arena, (size_t)(ncomp * sizeof(size_t)));
                    /* Restrict the denominator to the largest run's dominant
                     * ROW BAND.  Claim peeling stacks one band of lattice
                     * rows per layer, and rows of adjacent bands are
                     * lattice-adjacent but geometrically unrelated (band 1
                     * row 0 sits at the same z as band 0 row 0), so no legal
                     * face can join them, and a denominator spanning every
                     * band caps coherence at ~1/layers by construction.
                     * Row-run contiguity CANNOT find the boundary -- band
                     * occupancy is contiguous straight across it -- so the
                     * boundary is detected from the deliverable itself: an
                     * adjacent occupied row pair splits the denominator only
                     * when it is UNJOINABLE, i.e. the same-column vertex
                     * pairs it shares sit mostly farther apart in 3-D than
                     * the wrap gate (or it shares no column at all).  A real
                     * tear -- adjacent rows within the gate but in different
                     * components -- stays inside one band and keeps hurting
                     * coherence, which is the point of the gate. */
                    long lo_row = (long)lround(out->v_lo / dv);
                    long hi_row = (long)lround(out->v_hi / dv);
                    size_t nrow = (size_t)(hi_row - lo_row + 1);
                    int banded = 0;
                    for (i = 0; i < ncomp; i++) inrun[i] = 0;
                    if (nrow > 1 && nrow < (size_t)1 << 22 &&
                        (double)nrow * (double)ncol < 2.5e8) {
                        int32_t *cellv = (int32_t *)ARENA_ALLOC(
                            arena, (size_t)(nrow * ncol * sizeof(int32_t)));
                        size_t *bandof = (size_t *)ARENA_ALLOC(
                            arena, (size_t)(nrow * sizeof(size_t)));
                        size_t rr = 0, nband = 0;
                        for (rr = 0; rr < nrow * ncol; rr++) cellv[rr] = -1;
                        for (i = 0; i < nv; i++) {
                            long cc = (long)lround((double)UV[i*2+0]/du) - lo_col;
                            long rw2 = (long)lround((double)UV[i*2+1]/dv) - lo_row;
                            if (cc >= 0 && (size_t)cc < ncol &&
                                run_of_col[cc] == best &&
                                rw2 >= 0 && (size_t)rw2 < nrow)
                                cellv[(size_t)rw2 * ncol + (size_t)cc] =
                                    (int32_t)i;
                        }
                        for (rr = 0; rr < nrow; rr++) {
                            int rowocc = 0;
                            size_t c2 = 0;
                            for (c2 = 0; c2 < ncol; c2++)
                                if (cellv[rr * ncol + c2] >= 0) {
                                    rowocc = 1;
                                    break;
                                }
                            if (!rowocc) { bandof[rr] = (size_t)-1; continue; }
                            if (nband > 0 && rr > 0 &&
                                bandof[rr-1] != (size_t)-1) {
                                size_t common = 0, farcnt = 0;
                                for (c2 = 0; c2 < ncol; c2++) {
                                    int32_t a = cellv[(rr-1) * ncol + c2];
                                    int32_t b = cellv[rr * ncol + c2];
                                    if (a < 0 || b < 0) continue;
                                    common++;
                                    if (rv_d3(&V[(size_t)a*3],
                                              &V[(size_t)b*3]) > RV_WRAP_GATE)
                                        farcnt++;
                                }
                                if (common > 0 && farcnt * 2 <= common) {
                                    bandof[rr] = bandof[rr-1];
                                    continue;
                                }
                            }
                            bandof[rr] = nband++;
                        }
                        if (nband > 0) {
                            size_t *bcells = (size_t *)ARENA_ALLOC(
                                arena, (size_t)(nband * sizeof(size_t)));
                            size_t bbest = 0;
                            for (rr = 0; rr < nband; rr++) bcells[rr] = 0;
                            for (rr = 0; rr < nrow; rr++) {
                                size_t c2 = 0;
                                if (bandof[rr] == (size_t)-1) continue;
                                for (c2 = 0; c2 < ncol; c2++)
                                    if (cellv[rr * ncol + c2] >= 0)
                                        bcells[bandof[rr]]++;
                            }
                            for (rr = 0; rr < nband; rr++)
                                if (bcells[rr] > bcells[bbest]) bbest = rr;
                            out->best_run_cells = bcells[bbest];
                            out->row_bands = nband;
                            for (rr = 0; rr < nrow; rr++) {
                                size_t c2 = 0;
                                if (bandof[rr] != bbest) continue;
                                for (c2 = 0; c2 < ncol; c2++) {
                                    int32_t a = cellv[rr * ncol + c2];
                                    if (a >= 0) inrun[label[a]]++;
                                }
                            }
                            banded = 1;
                        }
                    }
                    if (!banded) {
                        out->row_bands = 1;
                        for (i = 0; i < nv; i++) {
                            long cc = (long)lround((double)UV[i*2+0]/du) - lo_col;
                            if (cc >= 0 && (size_t)cc < ncol &&
                                run_of_col[cc] == best)
                                inrun[label[i]]++;
                        }
                    }
                    for (i = 0; i < ncomp; i++)
                        if (inrun[i] > out->best_run_largest_comp)
                            out->best_run_largest_comp = inrun[i];
                }
            }
        }
    }

    /* ---- depth seams (see RV_DEPTH_* above) ------------------------- */
    if (!have_axis)
        snprintf(out->depth_skip, sizeof out->depth_skip, "no axis: needs --umb-y/--umb-x");
    if (have_axis) {
        long lo_col2 = (long)lround(out->u_lo / du);
        long hi_col2 = (long)lround(out->u_hi / du);
        long lo_row2 = (long)lround(out->v_lo / dv);
        long hi_row2 = (long)lround(out->v_hi / dv);
        size_t wcols = (size_t)(hi_col2 - lo_col2 + 1);
        size_t wrows = (size_t)(hi_row2 - lo_row2 + 1);
        if (!(wcols >= 2 && wrows >= 2))
            snprintf(out->depth_skip, sizeof out->depth_skip,
                     "the lattice is %zu x %zu cells: degenerate", wrows, wcols);
        else if (!((double)wcols * (double)wrows < RV_DEPTH_MAX_CELLS))
            snprintf(out->depth_skip, sizeof out->depth_skip,
                     "the lattice is %zu x %zu = %.3g cells, over the %.3g cap (not an axis problem)",
                     wrows, wcols, (double)wrows * (double)wcols, RV_DEPTH_MAX_CELLS);
        if (wcols >= 2 && wrows >= 2 &&
            (double)wcols * (double)wrows < RV_DEPTH_MAX_CELLS) {
            size_t ncell2 = wcols * wrows;
            float *cr = (float *)ARENA_ALLOC(
                arena, (size_t)(ncell2 * sizeof(float)));
            float *hp = (float *)ARENA_ALLOC(
                arena, (size_t)(ncell2 * sizeof(float)));
            unsigned char *mk = (unsigned char *)ARENA_ALLOC(
                arena, (size_t)ncell2);
            unsigned char *topo = (unsigned char *)ARENA_CALLOC(
                arena, (size_t)ncell2, 1);
            int32_t *cell_component = (int32_t *)ARENA_ALLOC(
                arena, (size_t)(ncell2 * sizeof(int32_t)));
            int32_t *cell_vertex = depth_csv != NULL
                ? (int32_t *)ARENA_ALLOC(
                    arena, (size_t)(ncell2 * sizeof(int32_t)))
                : NULL;
            size_t ci2 = 0;
            for (ci2 = 0; ci2 < ncell2; ci2++) {
                cr[ci2] = -1.0f;
                cell_component[ci2] = -1;
                if (cell_vertex != NULL) cell_vertex[ci2] = -1;
            }
            for (i = 0; i < nv; i++) {
                long cc = (long)lround((double)UV[i*2+0]/du) - lo_col2;
                long rw2 = (long)lround((double)UV[i*2+1]/dv) - lo_row2;
                double dy2 = (double)V[i*3+1] - umb_y;
                double dx2 = (double)V[i*3+2] - umb_x;
                if (cc >= 0 && (size_t)cc < wcols &&
                    rw2 >= 0 && (size_t)rw2 < wrows)
                    cr[(size_t)rw2 * wcols + (size_t)cc] =
                        (float)sqrt(dy2*dy2 + dx2*dx2);
                if (cc >= 0 && (size_t)cc < wcols &&
                    rw2 >= 0 && (size_t)rw2 < wrows)
                    cell_component[(size_t)rw2 * wcols + (size_t)cc] =
                        label[i];
                if (cell_vertex != NULL && cc >= 0 && (size_t)cc < wcols &&
                    rw2 >= 0 && (size_t)rw2 < wrows)
                    cell_vertex[(size_t)rw2 * wcols + (size_t)cc] = (int32_t)i;
            }
            /* Encode only axis-aligned edges that the emitted triangulation
             * actually owns: bit 0 joins this cell to +U, bit 1 to +V. */
            for (i = 0; i < ne; i++) {
                size_t a = (size_t)edges[i].a, b = (size_t)edges[i].b;
                long ca = (long)lround((double)UV[a*2+0]/du) - lo_col2;
                long ra = (long)lround((double)UV[a*2+1]/dv) - lo_row2;
                long cb = (long)lround((double)UV[b*2+0]/du) - lo_col2;
                long rb = (long)lround((double)UV[b*2+1]/dv) - lo_row2;
                if (ca < 0 || cb < 0 || ra < 0 || rb < 0 ||
                    (size_t)ca >= wcols || (size_t)cb >= wcols ||
                    (size_t)ra >= wrows || (size_t)rb >= wrows)
                    continue;
                if (ra == rb && labs(ca - cb) == 1) {
                    long left = ca < cb ? ca : cb;
                    topo[(size_t)ra * wcols + (size_t)left] |= 1u;
                } else if (ca == cb && labs(ra - rb) == 1) {
                    long top = ra < rb ? ra : rb;
                    topo[(size_t)top * wcols + (size_t)ca] |= 2u;
                }
            }
            {
                const int hw = RV_DEPTH_MEDIAN_WIN / 2;
                size_t rr2 = 0;
                for (rr2 = 0; rr2 < wrows; rr2++) {
                    size_t cc2 = 0;
                    for (cc2 = 0; cc2 < wcols; cc2++) {
                        float win[RV_DEPTH_MEDIAN_WIN * RV_DEPTH_MEDIAN_WIN];
                        int nw = 0, aa = 0;
                        size_t at2 = rr2 * wcols + cc2;
                        int32_t owner = cell_component[at2];
                        if (cr[at2] < 0.0f) { hp[at2] = 0.0f; continue; }
                        for (aa = -hw; aa <= hw; aa++) {
                            int bb = 0;
                            long rw3 = (long)rr2 + aa;
                            if (rw3 < 0 || rw3 >= (long)wrows) continue;
                            for (bb = -hw; bb <= hw; bb++) {
                                long cc3 = (long)cc2 + bb;
                                float v2;
                                if (cc3 < 0 || cc3 >= (long)wcols) continue;
                                v2 = cr[(size_t)rw3 * wcols + (size_t)cc3];
                                if (v2 < 0.0f) continue;
                                if (cell_component[(size_t)rw3 * wcols +
                                                   (size_t)cc3] != owner)
                                    continue;
                                /* insertion into sorted prefix */
                                {
                                    int qq = nw++;
                                    while (qq > 0 && win[qq-1] > v2) {
                                        win[qq] = win[qq-1];
                                        qq--;
                                    }
                                    win[qq] = v2;
                                }
                            }
                        }
                        hp[at2] = nw > 0 ? cr[at2] - win[nw/2] : 0.0f;
                    }
                }
            }
            {
                size_t rr2 = 0, cc2 = 0, runlen = 0, run_start = 0;
                /* vertical seams: steps between column pairs, runs down rows */
                for (ci2 = 0; ci2 < ncell2; ci2++) mk[ci2] = 0;
                for (rr2 = 0; rr2 < wrows; rr2++)
                    for (cc2 = 0; cc2 + 1 < wcols; cc2++) {
                        size_t at2 = rr2 * wcols + cc2;
                        if (cr[at2] < 0.0f || cr[at2+1] < 0.0f) continue;
                        if (!(topo[at2] & 1u)) continue;
                        out->depth_edges_total++;
                        if (fabs((double)hp[at2] - (double)hp[at2+1]) >=
                            RV_DEPTH_STEP)
                            mk[at2] = 1;
                    }
                for (cc2 = 0; cc2 + 1 < wcols; cc2++) {
                    runlen = 0;
                    for (rr2 = 0; rr2 < wrows; rr2++) {
                        if (mk[rr2 * wcols + cc2]) {
                            if (runlen == 0) run_start = rr2;
                            runlen++;
                        }
                        else {
                            if (runlen >= RV_DEPTH_RUN_MIN) {
                                out->depth_seam_runs++;
                                out->depth_seam_edges += runlen;
                                rv_write_depth_run(
                                    depth_csv, out->depth_seam_runs, 1, cc2,
                                    run_start, runlen, wcols, lo_col2, lo_row2,
                                    du, dv, V, UV, cr, hp, cell_vertex,
                                    cell_component, labels, materials, support);
                            }
                            runlen = 0;
                        }
                    }
                    if (runlen >= RV_DEPTH_RUN_MIN) {
                        out->depth_seam_runs++;
                        out->depth_seam_edges += runlen;
                        rv_write_depth_run(
                            depth_csv, out->depth_seam_runs, 1, cc2,
                            run_start, runlen, wcols, lo_col2, lo_row2,
                            du, dv, V, UV, cr, hp, cell_vertex,
                            cell_component, labels, materials, support);
                    }
                }
                /* horizontal seams: steps between row pairs, runs along cols */
                for (ci2 = 0; ci2 < ncell2; ci2++) mk[ci2] = 0;
                for (rr2 = 0; rr2 + 1 < wrows; rr2++)
                    for (cc2 = 0; cc2 < wcols; cc2++) {
                        size_t at2 = rr2 * wcols + cc2;
                        if (cr[at2] < 0.0f || cr[at2 + wcols] < 0.0f)
                            continue;
                        if (!(topo[at2] & 2u)) continue;
                        out->depth_edges_total++;
                        if (fabs((double)hp[at2] -
                                 (double)hp[at2 + wcols]) >= RV_DEPTH_STEP)
                            mk[at2] = 1;
                    }
                for (rr2 = 0; rr2 + 1 < wrows; rr2++) {
                    runlen = 0;
                    for (cc2 = 0; cc2 < wcols; cc2++) {
                        if (mk[rr2 * wcols + cc2]) {
                            if (runlen == 0) run_start = cc2;
                            runlen++;
                        }
                        else {
                            if (runlen >= RV_DEPTH_RUN_MIN) {
                                out->depth_seam_runs++;
                                out->depth_seam_edges += runlen;
                                rv_write_depth_run(
                                    depth_csv, out->depth_seam_runs, 0, rr2,
                                    run_start, runlen, wcols, lo_col2, lo_row2,
                                    du, dv, V, UV, cr, hp, cell_vertex,
                                    cell_component, labels, materials, support);
                            }
                            runlen = 0;
                        }
                    }
                    if (runlen >= RV_DEPTH_RUN_MIN) {
                        out->depth_seam_runs++;
                        out->depth_seam_edges += runlen;
                        rv_write_depth_run(
                            depth_csv, out->depth_seam_runs, 0, rr2,
                            run_start, runlen, wcols, lo_col2, lo_row2,
                            du, dv, V, UV, cr, hp, cell_vertex,
                            cell_component, labels, materials, support);
                    }
                }
            }
            out->depth_measured = 1;
        }
    }

    {
        double cols = (out->u_hi - out->u_lo) / du + 1.0;
        double rows = (out->v_hi - out->v_lo) / dv + 1.0;
        if (cols < 1.0) cols = 1.0;
        if (rows < 1.0) rows = 1.0;
        out->cells_rect = (size_t)(cols * rows);
    }
    return 0;
}

/* ---- gates ------------------------------------------------------------- */

typedef struct { const char *name; int pass; int measured; char detail[192]; } RvGate;

/* Evaluate every gate.  `gates` must hold at least 6 entries.  Returns the
 * number of FAILED gates; an unmeasured gate is reported and counts as a
 * failure, because shipping on an unmeasured criterion is how the 4x5x5
 * cross-section claim happened in the first place. */
static int rv_gates(const RvReport *r, RvGate *g, size_t *ngate)
{
    size_t n = 0;
    int fails = 0;
    double share = r->cells_occupied
                 ? (double)r->largest_cells / (double)r->cells_occupied : 0.0;
    double fill = r->cells_rect
                ? (double)r->cells_occupied / (double)r->cells_rect : 0.0;

    g[n].name = "TANGLE/long-edges"; g[n].measured = 1;
    g[n].pass = r->edges_long == 0;
    snprintf(g[n].detail, sizeof(g[n].detail),
             "%zu of %zu emitted edges over %.1f vox (max %.2f)",
             r->edges_long, r->edges_total, RV_WRAP_GATE, r->edge_len_max);
    n++;

    g[n].name = "TANGLE/cross-wrap"; g[n].measured = r->have_axis;
    g[n].pass = r->have_axis && r->edges_cross_wrap == 0;
    if (r->have_axis)
        snprintf(g[n].detail, sizeof(g[n].detail),
                 "%zu emitted edges step >= half a pitch radially (max dr %.2f)",
                 r->edges_cross_wrap, r->edge_dr_max);
    else
        snprintf(g[n].detail, sizeof(g[n].detail),
                 "NOT MEASURED: pass --umb-y/--umb-x");
    n++;

    /* A cross-LANE face is not a defect: every emitted edge is already under
     * the wrap gate, comfortably inside the 7-vox inter-wrap clearance, so
     * such a face is the same physical sheet measured by two reconstruction
     * lanes -- a chart join, and the changelog records that discarding those
     * shatters the sheet 39 -> 1,945 components.  What must never happen is
     * joining two MATERIAL lineages: that is the global stitcher's decision,
     * never a raster coincidence.  Lane crossings are reported below. */
    g[n].name = "WINDING/material-purity"; g[n].measured = r->have_materials;
    g[n].pass = r->have_materials && r->faces_cross_material == 0;
    if (r->have_materials)
        snprintf(g[n].detail, sizeof(g[n].detail),
                 "%zu faces join two material lineages (%zu join two lanes, "
                 "which is a chart join, not a defect)",
                 r->faces_cross_material, r->faces_cross_lane);
    else
        snprintf(g[n].detail, sizeof(g[n].detail),
                 "NOT MEASURED: pass --materials <material_identity.i32>");
    n++;

    g[n].name = "WINDING/single-cover"; g[n].measured = 1;
    g[n].pass = r->cells_duplicated == 0;
    snprintf(g[n].detail, sizeof(g[n].detail),
             "%zu vertices share a lattice cell with another",
             r->cells_duplicated);
    n++;

    g[n].name = "DISTORTION/anisotropy"; g[n].measured = r->aniso_samples > 0;
    g[n].pass = r->aniso_samples > 0 && r->aniso_frac_over <= RV_ANISO_FRAC_MAX;
    snprintf(g[n].detail, sizeof(g[n].detail),
             "%.4f%% of triangles beyond %.1fx OR under 1/%.1fx the sheet "
             "median scale (p50 %.2f p90 %.2f p99 %.2f max %.1f)",
             100.0 * r->aniso_frac_over, RV_ANISO_LIMIT, RV_ANISO_LIMIT,
             r->aniso_p50, r->aniso_p90, r->aniso_p99, r->aniso_max);
    n++;

    {
        double cov = r->source_verts
                   ? (double)r->source_covered / (double)r->source_verts : 0.0;
        g[n].name = "COVERAGE/source"; g[n].measured = r->have_source;
        g[n].pass = r->have_source && cov >= RV_COVER_MIN;
        if (r->have_source)
            snprintf(g[n].detail, sizeof(g[n].detail),
                     "%.1f%% of %zu source vertices lie within %.1f vox of the "
                     "ribbon", 100.0 * cov, r->source_verts, RV_COVER_RADIUS);
        else
            snprintf(g[n].detail, sizeof(g[n].detail),
                     "NOT MEASURED: pass --source <welded.vmesh>");
        n++;
    }
    {
        double coh = r->best_run_cells
                   ? (double)r->best_run_largest_comp / (double)r->best_run_cells
                   : 0.0;
        g[n].name = "FRAGMENTATION/coherence"; g[n].measured = r->best_run_cells > 0;
        g[n].pass = r->best_run_cells > 0 && coh >= RV_RUN_COHERENCE_MIN;
        snprintf(g[n].detail, sizeof(g[n].detail),
                 "largest run x band holds %zu cells, %.1f%% in one component; "
                 "%zu runs, %zu bands, %zu components, largest %.1f%% overall, "
                 "fill %.1f%%",
                 r->best_run_cells, 100.0 * coh, r->column_runs, r->row_bands,
                 r->n_components, 100.0 * share, 100.0 * fill);
        n++;
    }

    {
        double share = r->depth_edges_total
            ? (double)r->depth_seam_edges / (double)r->depth_edges_total
            : 0.0;
        g[n].name = "TEXTURE/depth-seams";
        g[n].measured = r->depth_measured;
        g[n].pass = r->depth_measured && share <= RV_DEPTH_SEAM_MAX;
        if (r->depth_measured)
            snprintf(g[n].detail, sizeof(g[n].detail),
                     "%zu plateau seam runs; %.3f%% of %zu emitted lattice edges "
                     "step >= %.2f vox after slope removal (max %.3f%%)",
                     r->depth_seam_runs, 100.0 * share,
                     r->depth_edges_total, RV_DEPTH_STEP,
                     100.0 * RV_DEPTH_SEAM_MAX);
        else
            snprintf(g[n].detail, sizeof(g[n].detail), "NOT MEASURED: %s",
                     r->depth_skip[0] ? r->depth_skip : "no axis: needs --umb-y/--umb-x");
        n++;
    }

    for (*ngate = 0, *ngate = n; n-- > 0; )
        if (!g[n].pass) fails++;
    return fails;
}

static void rv_print(const RvReport *r, const RvGate *g, size_t ngate,
                     int fails)
{
    size_t i = 0;
    fprintf(stderr, "[ribbon_verdict]\n");
    for (i = 0; i < ngate; i++)
        fprintf(stderr, "  %-24s %s  %s\n", g[i].name,
                g[i].pass ? "PASS" : (g[i].measured ? "FAIL" : "ABSENT"),
                g[i].detail);
    if (r->have_support)
        fprintf(stderr, "  %-24s      %zu measured vertices\n",
                "support", r->verts_supported);
    fprintf(stderr, "  VERDICT: %s (%d gate(s) failed)\n",
            fails == 0 ? "SHIPPABLE" : "NOT SHIPPABLE", fails);
}

static int rv_write_report(const char *path, const char *input,
                           const RvReport *r, const RvGate *g, size_t ngate,
                           int fails)
{
    FILE *fp = fopen(path, "wb");
    size_t i = 0;
    double share = r->cells_occupied
                 ? (double)r->largest_cells / (double)r->cells_occupied : 0.0;
    double fill = r->cells_rect
                ? (double)r->cells_occupied / (double)r->cells_rect : 0.0;
    if (fp == NULL) return -1;
    fprintf(fp, "{\n  \"schema\": \"vesuvius-ribbon-verdict-v1\",\n");
        fputs("  \"input\": \"", fp);
    /* JSON-escape the path: Windows backslashes made every report on this
     * platform unparseable (quadribbon then ran blind: "verdict unavailable"
     * at every stage, 2026-09-02) */
    for (const char *c = input; *c; c++) {
        if (*c == 92 || *c == 34) fputc(92, fp);
        fputc(*c, fp);
    }
    fputs("\",\n", fp);
    fprintf(fp, "  \"verdict\": \"%s\",\n  \"gates_failed\": %d,\n",
            fails == 0 ? "SHIPPABLE" : "NOT_SHIPPABLE", fails);
    fprintf(fp, "  \"gates\": [");
    for (i = 0; i < ngate; i++)
        fprintf(fp, "%s\n    { \"name\": \"%s\", \"pass\": %s, "
                    "\"measured\": %s, \"detail\": \"%s\" }",
                i ? "," : "", g[i].name, g[i].pass ? "true" : "false",
                g[i].measured ? "true" : "false", g[i].detail);
    fprintf(fp, "\n  ],\n");
    fprintf(fp, "  \"tangle\": { \"edges_total\": %zu, \"edges_long\": %zu, "
                "\"edges_cross_wrap\": %zu, \"edge_len_max\": %.4f, "
                "\"edge_dr_max\": %.4f },\n",
            r->edges_total, r->edges_long, r->edges_cross_wrap,
            r->edge_len_max, r->edge_dr_max);
    fprintf(fp, "  \"winding\": { \"faces_cross_lane\": %zu, "
                "\"faces_cross_material\": %zu, "
                "\"cells_duplicated\": %zu },\n",
            r->faces_cross_lane, r->faces_cross_material,
            r->cells_duplicated);
    fprintf(fp, "  \"distortion\": { \"samples\": %zu, \"p50\": %.4f, "
                "\"p90\": %.4f, \"p99\": %.4f, \"max\": %.4f, "
                "\"frac_over_limit\": %.6f, \"core_radius\": %.1f, "
                "\"core_faces_excluded\": %zu },\n",
            r->aniso_samples, r->aniso_p50, r->aniso_p90, r->aniso_p99,
            r->aniso_max, r->aniso_frac_over,
            r->core_radius, r->core_faces_excluded);
    fprintf(fp, "  \"coverage\": { \"cells_occupied\": %zu, "
                "\"cells_rect\": %zu, \"fill\": %.6f, "
                "\"verts_supported\": %zu,\n"
                "                 \"source_verts\": %zu, \"source_covered\": %zu, "
                "\"source_frac\": %.6f,\n"
                "                 \"u_lo\": %.3f, \"u_hi\": %.3f, "
                "\"v_lo\": %.3f, \"v_hi\": %.3f },\n",
            r->cells_occupied, r->cells_rect, fill, r->verts_supported,
            r->source_verts, r->source_covered,
            r->source_verts ? (double)r->source_covered / (double)r->source_verts : 0.0,
            r->u_lo, r->u_hi, r->v_lo, r->v_hi);
    fprintf(fp, "  \"fragmentation\": { \"n_components\": %zu, "
                "\"largest_cells\": %zu, \"largest_share\": %.6f, "
                "\"column_runs\": %zu, \"row_bands\": %zu, "
                "\"best_run_cells\": %zu, \"best_run_largest_comp\": %zu },\n",
            r->n_components, r->largest_cells, share,
            r->column_runs, r->row_bands,
            r->best_run_cells, r->best_run_largest_comp);
    fprintf(fp, "  \"depth_seams\": { \"runs\": %zu, \"seam_edges\": %zu, "
                "\"edges_total\": %zu, \"seam_share\": %.6f, \"measured\": %s }\n}\n",
            r->depth_seam_runs, r->depth_seam_edges, r->depth_edges_total,   /* skip reason: see the gate line */
            r->depth_edges_total ? (double)r->depth_seam_edges / (double)r->depth_edges_total : 0.0,
            r->depth_measured ? "true" : "false");
    fclose(fp);
    return 0;
}

/* ---- sidecar loading --------------------------------------------------- */

/* Sidecars are raw payloads with exactly one entry per VMESH vertex, so the
 * file length IS the contract.  A wrong-length sidecar is a stale artifact
 * from another run and must be refused, never truncated to fit. */
static void *rv_load_sidecar(const char *path, size_t nv, size_t esize)
{
    FILE *fp = fopen(path, "rb");
    void *buf = NULL;
    long len = 0;
    if (fp == NULL) {
        fprintf(stderr, "ERROR: cannot open %s\n", path);
        return NULL;
    }
    if (fseek(fp, 0, SEEK_END) != 0) { fclose(fp); return NULL; }
    len = ftell(fp);
    if (len < 0 || (size_t)len != nv * esize) {
        fprintf(stderr, "ERROR: %s is %ld bytes, expected %zu "
                        "(%zu verts x %zu)\n",
                path, len, nv * esize, nv, esize);
        fclose(fp);
        return NULL;
    }
    rewind(fp);
    buf = malloc(nv * esize);
    if (buf == NULL || fread(buf, esize, nv, fp) != nv) {
        fprintf(stderr, "ERROR: short read on %s\n", path);
        free(buf);
        fclose(fp);
        return NULL;
    }
    fclose(fp);
    return buf;
}

/* Lattice pitch straight from the deliverable: the median gap between
 * consecutive DISTINCT values of one UV component.  The emitted ribbon is
 * lattice-quantized by construction, so this recovers du/dv exactly; taking
 * it from the data removes the one flag that was silently wrong for every
 * verdict before 2026-09-01 (--dv 1 against a slice_h=2 lattice made every
 * occupied row its own band). */
static double rv_detect_pitch(const float *uv, size_t nv, int comp)
{
    double *vals = NULL, *gaps = NULL, med = 0.0;
    size_t i = 0, nun = 0, ng = 0;
    if (nv < 2) return 1.0;
    vals = (double *)malloc(nv * sizeof *vals);
    if (vals == NULL) return 1.0;
    for (i = 0; i < nv; i++) vals[i] = (double)uv[i * 2 + comp];
    qsort(vals, nv, sizeof *vals, rv_cmp_double);
    for (i = 1, nun = 1; i < nv; i++)
        if (vals[i] > vals[nun - 1] + 1e-6) vals[nun++] = vals[i];
    if (nun < 2) { free(vals); return 1.0; }
    gaps = (double *)malloc((nun - 1) * sizeof *gaps);
    if (gaps == NULL) { free(vals); return 1.0; }
    for (i = 1; i < nun; i++) gaps[ng++] = vals[i] - vals[i - 1];
    qsort(gaps, ng, sizeof *gaps, rv_cmp_double);
    med = gaps[ng / 2];
    free(gaps);
    free(vals);
    return med > 1e-6 ? med : 1.0;
}

/* ---- selftest ---------------------------------------------------------- */

static int rv_check(int cond, const char *what, int *fails)
{
    fprintf(stderr, "[ribbon_verdict selftest] %s -> %s\n",
            what, cond ? "ok" : "FAIL");
    if (!cond) (*fails)++;
    return cond;
}

/* Two lattice rows of `ncol` columns at radius `r0`, plus (when `spoke` is
 * set) one column displaced outward by `dr` -- the exact shape of the defect
 * this tool exists to catch. */
static void rv_build_strip(float *verts, float *uv, int32_t *faces,
                           int ncol, double r0, int spoke, double dr,
                           double umb_y, double umb_x)
{
    int k = 0, j = 0, nf = 0;
    for (k = 0; k < 2; k++) {
        for (j = 0; j < ncol; j++) {
            size_t vi = (size_t)(k * ncol + j);
            double ang = 0.02 * j;
            double rr = r0 + ((spoke && j == ncol / 2) ? dr : 0.0);
            verts[vi*3+0] = (float)(2.0 * k);
            verts[vi*3+1] = (float)(umb_y + rr * cos(ang));
            verts[vi*3+2] = (float)(umb_x + rr * sin(ang));
            uv[vi*2+0] = (float)(2.0 * j);
            uv[vi*2+1] = (float)(1.0 * k);
        }
    }
    for (j = 0; j + 1 < ncol; j++) {
        int32_t a = (int32_t)j, b = (int32_t)(j + 1);
        int32_t c = (int32_t)(ncol + j), d = (int32_t)(ncol + j + 1);
        faces[nf*3+0] = a; faces[nf*3+1] = b; faces[nf*3+2] = c; nf++;
        faces[nf*3+0] = b; faces[nf*3+1] = d; faces[nf*3+2] = c; nf++;
    }
}

/* An interior three-column depth plateau.  `cut` removes the two chart-boundary
 * edge strips while leaving their UV columns adjacent, exercising the rule
 * that atlas adjacency without emitted topology is not a surface seam. */
static size_t rv_build_depth_plateau(float *verts, float *uv, int32_t *faces,
                                     int nrow, int ncol, int cut,
                                     double umb_y, double umb_x)
{
    size_t nf = 0;
    int k = 0, j = 0;
    const int p0 = 4, p1 = 6;
    for (k = 0; k < nrow; k++) {
        for (j = 0; j < ncol; j++) {
            size_t vi = (size_t)(k * ncol + j);
            double ang = 0.01 * (double)j;
            double radius = j >= p0 && j <= p1 ? 52.0 : 50.0;
            verts[vi*3+0] = (float)(2.0 * k);
            verts[vi*3+1] = (float)(umb_y + radius * cos(ang));
            verts[vi*3+2] = (float)(umb_x + radius * sin(ang));
            uv[vi*2+0] = (float)(2.0 * j);
            uv[vi*2+1] = (float)(2.0 * k);
        }
    }
    for (k = 0; k + 1 < nrow; k++) {
        for (j = 0; j + 1 < ncol; j++) {
            int32_t a, b, c, d;
            if (cut && (j == p0 - 1 || j == p1)) continue;
            a = (int32_t)(k * ncol + j);
            b = a + 1;
            c = (int32_t)((k + 1) * ncol + j);
            d = c + 1;
            faces[nf*3+0] = a; faces[nf*3+1] = b;
            faces[nf*3+2] = c; nf++;
            faces[nf*3+0] = b; faces[nf*3+1] = d;
            faces[nf*3+2] = c; nf++;
        }
    }
    return nf;
}

static int rv_selftest(void)
{
    int fails = 0;
    const double umb_y = 100.0, umb_x = 100.0, pitch = 9.5;
    const int ncol = 8;
    Arena_T arena = Arena_new();
    float verts[2*8*3], uv[2*8*2];
    int32_t faces[2*7*3];
    MeshBinData mesh;
    RvReport rep;
    RvGate g[10];
    size_t ngate = 0;

    memset(&mesh, 0, sizeof(mesh));
    mesh.verts = verts; mesh.uv = uv; mesh.faces = faces;
    mesh.nv = (size_t)(2 * ncol); mesh.nf = (size_t)(2 * (ncol - 1));

    /* (1) a clean strip passes the tangle gates */
    rv_build_strip(verts, uv, faces, ncol, 50.0, 0, 0.0, umb_y, umb_x);
    rv_check(rv_measure(arena, &mesh, NULL, NULL, NULL, NULL, 1, umb_y, umb_x, pitch, 0.0,
                        2.0, 1.0, NULL, &rep) == 0,
             "clean strip measures", &fails);
    rv_check(rep.edges_long == 0, "clean strip has no long edge", &fails);
    rv_check(rep.edges_cross_wrap == 0, "clean strip has no cross-wrap edge",
             &fails);
    rv_check(rep.cells_duplicated == 0, "clean strip covers each cell once",
             &fails);
    rv_check(rep.n_components == 1, "clean strip is one component", &fails);

    /* (2) a radial spoke is caught by BOTH tangle gates.  One displaced
     *     column touches 2 columns x 2 rows worth of edges. */
    rv_build_strip(verts, uv, faces, ncol, 50.0, 1, 16.0, umb_y, umb_x);
    rv_check(rv_measure(arena, &mesh, NULL, NULL, NULL, NULL, 1, umb_y, umb_x, pitch, 0.0,
                        2.0, 1.0, NULL, &rep) == 0,
             "spoke strip measures", &fails);
    rv_check(rep.edges_long > 0, "spoke raises the long-edge gate", &fails);
    rv_check(rep.edges_cross_wrap > 0, "spoke raises the cross-wrap gate",
             &fails);
    rv_check(rep.edge_dr_max > 15.0, "spoke dr is reported (~16 vox)", &fails);
    {
        int f2 = rv_gates(&rep, g, &ngate) > 0;
        rv_check(f2, "spoke strip is NOT SHIPPABLE", &fails);
        rv_check(ngate == 8, "eight gates are evaluated", &fails);
    }

    /* (3) an unmeasured gate is ABSENT, never a silent pass */
    rv_build_strip(verts, uv, faces, ncol, 50.0, 0, 0.0, umb_y, umb_x);
    rv_check(rv_measure(arena, &mesh, NULL, NULL, NULL, NULL, 0, 0.0, 0.0, pitch, 0.0,
                        2.0, 1.0, NULL, &rep) == 0,
             "no-axis measures", &fails);
    (void)rv_gates(&rep, g, &ngate);
    rv_check(!g[1].measured && !g[1].pass,
             "cross-wrap gate without an axis is ABSENT and fails", &fails);
    rv_check(!g[2].measured && !g[2].pass,
             "material gate without materials is ABSENT and fails", &fails);
    /* an absent gate must SAY WHY: "needs --umb" and "over the cell cap" are different
     * failures, and a pile with both a table and an umbilicus was being told it had neither */
    rv_check(strstr(rep.depth_skip, "no axis") != NULL,
             "an axis-less depth gate names the axis as its reason", &fails);
    {
        size_t k = 0;
        for (k = 0; k < ngate; k++)
            if (strcmp(g[k].name, "TEXTURE/depth-seams") == 0) break;
        rv_check(k < ngate && !g[k].measured && strstr(g[k].detail, "no axis") != NULL,
                 "the depth gate carries its skip reason to the gate line", &fails);
    }

    /* (4) lane purity: a face spanning two lanes is caught */
    {
        int32_t labels[2*8];
        size_t i = 0;
        for (i = 0; i < mesh.nv; i++) labels[i] = (int32_t)(i % 2);
        rv_check(rv_measure(arena, &mesh, labels, NULL, NULL, NULL, 1, umb_y, umb_x,
                            pitch, 0.0, 2.0, 1.0, NULL, &rep) == 0,
                 "labelled mesh measures", &fails);
        rv_check(rep.faces_cross_lane == mesh.nf,
                 "alternating lanes are counted as crossings", &fails);
        for (i = 0; i < mesh.nv; i++) labels[i] = 7;
        (void)rv_measure(arena, &mesh, labels, NULL, NULL, NULL, 1, umb_y, umb_x, pitch,
                          0.0, 2.0, 1.0, NULL, &rep);
        rv_check(rep.faces_cross_lane == 0, "one lane yields no crossing",
                 &fails);
        /* -1 is "derived", not a lane: it must not count as a crossing */
        for (i = 0; i < mesh.nv; i++) labels[i] = (i == 0) ? -1 : 7;
        (void)rv_measure(arena, &mesh, labels, NULL, NULL, NULL, 1, umb_y, umb_x, pitch,
                          0.0, 2.0, 1.0, NULL, &rep);
        rv_check(rep.faces_cross_lane == 0,
                 "derived (-1) vertices are not lane crossings", &fails);
    }

    /* (5) duplicate cover is caught: collapse two columns onto one UV cell */
    rv_build_strip(verts, uv, faces, ncol, 50.0, 0, 0.0, umb_y, umb_x);
    uv[1*2+0] = uv[0*2+0];
    uv[1*2+1] = uv[0*2+1];
    (void)rv_measure(arena, &mesh, NULL, NULL, NULL, NULL, 1, umb_y, umb_x, pitch, 0.0,
                     2.0, 1.0, NULL, &rep);
    rv_check(rep.cells_duplicated == 1, "a stacked cell is counted once",
             &fails);

    /* (6) peel bands: lattice-adjacent rows the wrap gate cannot join must
     *     not dilute coherence (claim peeling stacks such bands), while a
     *     real tear INSIDE one band must still hurt it */
    {
        float v6[26*3], u6[26*2];
        int32_t f6[14*3];
        MeshBinData m6;
        int tear = 0;
        memset(&m6, 0, sizeof(m6));
        m6.verts = v6; m6.uv = u6; m6.faces = f6;
        for (tear = 0; tear < 2; tear++) {
            size_t vi = 0;
            int k = 0, j = 0, nf6 = 0;
            /* band 0: 2 rows x 8 cols at r=50 -- one sheet (or torn in two) */
            for (k = 0; k < 2; k++)
                for (j = 0; j < 8; j++, vi++) {
                    double ang = 0.02 * j;
                    v6[vi*3+0] = (float)(2.0 * k);
                    v6[vi*3+1] = (float)(umb_y + 50.0 * cos(ang));
                    v6[vi*3+2] = (float)(umb_x + 50.0 * sin(ang));
                    u6[vi*2+0] = (float)(2.0 * j);
                    u6[vi*2+1] = (float)k;
                }
            /* band 1: 2 rows x 5 cols at r=150 in the SAME lattice columns,
             * lattice-adjacent rows, ~100 vox away in 3-D: unjoinable
             * face-less confetti, exactly what a peel layer contributes */
            for (k = 0; k < 2; k++)
                for (j = 0; j < 5; j++, vi++) {
                    double ang = 0.02 * j;
                    v6[vi*3+0] = (float)(4.0 + 2.0 * k);
                    v6[vi*3+1] = (float)(umb_y + 150.0 * cos(ang));
                    v6[vi*3+2] = (float)(umb_x + 150.0 * sin(ang));
                    u6[vi*2+0] = (float)(2.0 * j);
                    u6[vi*2+1] = (float)(2 + k);
                }
            m6.nv = vi;
            for (j = 0; j + 1 < 8; j++) {
                int32_t a = (int32_t)j, b = (int32_t)(j + 1);
                int32_t c2 = (int32_t)(8 + j), d = (int32_t)(8 + j + 1);
                if (tear && j == 3) continue;
                f6[nf6*3+0] = a; f6[nf6*3+1] = b; f6[nf6*3+2] = c2; nf6++;
                f6[nf6*3+0] = b; f6[nf6*3+1] = d; f6[nf6*3+2] = c2; nf6++;
            }
            m6.nf = (size_t)nf6;
            rv_check(rv_measure(arena, &m6, NULL, NULL, NULL, NULL, 1, umb_y,
                                umb_x, pitch, 0.0, 2.0, 1.0, NULL, &rep) == 0,
                     tear ? "torn banded mesh measures"
                          : "banded mesh measures", &fails);
            rv_check(rep.row_bands == 2,
                     "unjoinable adjacent rows split into 2 bands", &fails);
            rv_check(rep.best_run_cells == 16,
                     "denominator is the dominant band only", &fails);
            if (!tear)
                rv_check(rep.best_run_largest_comp == 16,
                         "confetti band does not dilute coherence", &fails);
            else
                rv_check(rep.best_run_largest_comp == 8,
                         "a tear inside one band still hurts coherence",
                         &fails);
        }
    }

    /* (7) the depth gate follows emitted topology.  A real, connected plateau
     * is detected; cutting its two provenance boundaries makes three honest
     * adjacent atlas charts and removes those boundaries from the seam gate. */
    {
        enum { DR = 8, DC = 12 };
        float vd[DR*DC*3], ud[DR*DC*2];
        int32_t fd[2*(DR-1)*(DC-1)*3];
        MeshBinData md;
        memset(&md, 0, sizeof(md));
        md.verts = vd; md.uv = ud; md.faces = fd;
        md.nv = DR * DC;
        md.nf = rv_build_depth_plateau(
            vd, ud, fd, DR, DC, 0, umb_y, umb_x);
        rv_check(rv_measure(arena, &md, NULL, NULL, NULL, NULL, 1, umb_y,
                            umb_x, pitch, 0.0, 2.0, 2.0, NULL, &rep) == 0,
                 "connected depth plateau measures", &fails);
        rv_check(rep.depth_seam_runs >= 2 &&
                 rep.depth_seam_edges >= 2 * DR,
                 "connected depth plateau raises the seam gate", &fails);
        md.nf = rv_build_depth_plateau(
            vd, ud, fd, DR, DC, 1, umb_y, umb_x);
        rv_check(rv_measure(arena, &md, NULL, NULL, NULL, NULL, 1, umb_y,
                            umb_x, pitch, 0.0, 2.0, 2.0, NULL, &rep) == 0,
                 "cut depth plateau measures", &fails);
        rv_check(rep.depth_seam_runs == 0 && rep.depth_seam_edges == 0,
                 "cut chart boundaries are not depth seams", &fails);
    }

    Arena_dispose(&arena);
    fprintf(stderr, "ribbon_verdict --selftest: %s (%d failure(s))\n",
            fails ? "FAIL" : "PASS", fails);
    return fails ? 3 : 0;
}

/* ---- main -------------------------------------------------------------- */

static void rv_usage(void)
{
    fprintf(stderr,
        "usage: ribbon_verdict <ribbon.vmesh> [--labels recon.i32]\n"
        "                      [--source welded.vmesh]\n"
        "                      [--materials material_identity.i32]\n"
        "                      [--support support.u8] [--umb-y F --umb-x F]\n"
        "                      [--pitch F=9.5] [--du F=2] [--dv F=1]\n"
        "                      [--report out.json] [--depth-runs out.csv]\n"
        "       ribbon_verdict --selftest\n");
}

int main(int argc, char **argv)
{
    const char *in_path = NULL, *report = NULL, *depth_runs_path = NULL;
    const char *labels_path = NULL, *support_path = NULL;
    const char *materials_path = NULL, *source_path = NULL;
    double umb_y = 0.0, umb_x = 0.0, pitch = RV_PITCH_DEFAULT;
    double core_radius = RV_CORE_RADIUS_DEFAULT;
    double du = 0.0, dv = 0.0;   /* 0 = detect from the lattice itself */
    int have_y = 0, have_x = 0, i = 0, fails = 0, rc = 1;
    Arena_T arena = NULL;
    MeshBinData mesh;
    int32_t *labels = NULL, *materials = NULL;
    unsigned char *support = NULL;
    MeshBinData src;
    int have_source = 0;
    FILE *depth_csv = NULL;
    RvReport rep;
    RvGate gates[10];
    size_t ngate = 0;

    if (argc == 2 && strcmp(argv[1], "--selftest") == 0) return rv_selftest();

    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--labels") == 0 && i + 1 < argc)
            labels_path = argv[++i];
        else if (strcmp(argv[i], "--materials") == 0 && i + 1 < argc)
            materials_path = argv[++i];
        else if (strcmp(argv[i], "--source") == 0 && i + 1 < argc)
            source_path = argv[++i];
        else if (strcmp(argv[i], "--support") == 0 && i + 1 < argc)
            support_path = argv[++i];
        else if (strcmp(argv[i], "--report") == 0 && i + 1 < argc)
            report = argv[++i];
        else if (strcmp(argv[i], "--depth-runs") == 0 && i + 1 < argc)
            depth_runs_path = argv[++i];
        else if (strcmp(argv[i], "--umb-y") == 0 && i + 1 < argc)
            { umb_y = atof(argv[++i]); have_y = 1; }
        else if (strcmp(argv[i], "--umb-x") == 0 && i + 1 < argc)
            { umb_x = atof(argv[++i]); have_x = 1; }
        else if (strcmp(argv[i], "--pitch") == 0 && i + 1 < argc)
            pitch = atof(argv[++i]);
        else if (strcmp(argv[i], "--core-radius") == 0 && i + 1 < argc)
            core_radius = atof(argv[++i]);
        else if (strcmp(argv[i], "--du") == 0 && i + 1 < argc)
            du = atof(argv[++i]);
        else if (strcmp(argv[i], "--dv") == 0 && i + 1 < argc)
            dv = atof(argv[++i]);
        else if (argv[i][0] == '-') { rv_usage(); return 2; }
        else if (in_path == NULL) in_path = argv[i];
        else { rv_usage(); return 2; }
    }
    if (in_path == NULL || pitch <= 0.0 || du < 0.0 || dv < 0.0) {
        rv_usage();
        return 2;
    }

    memset(&mesh, 0, sizeof(mesh));
    if (MeshBin_read_malloc(in_path, &mesh) != 0) {
        fprintf(stderr, "ERROR: cannot read %s\n", in_path);
        return 1;
    }
    fprintf(stderr, "[ribbon_verdict] %s: %zu verts, %zu faces\n",
            in_path, mesh.nv, mesh.nf);
    if (mesh.uv == NULL) {
        fprintf(stderr, "ERROR: %s carries no UV; a ribbon verdict needs "
                        "the parameterization\n", in_path);
        goto done;
    }
    if (du <= 0.0) du = rv_detect_pitch(mesh.uv, mesh.nv, 0);
    if (dv <= 0.0) dv = rv_detect_pitch(mesh.uv, mesh.nv, 1);
    fprintf(stderr, "[ribbon_verdict] lattice pitch du=%.3f dv=%.3f\n",
            du, dv);
    if (labels_path != NULL) {
        labels = (int32_t *)rv_load_sidecar(labels_path, mesh.nv,
                                            sizeof(int32_t));
        if (labels == NULL) goto done;
    }
    if (materials_path != NULL) {
        materials = (int32_t *)rv_load_sidecar(materials_path, mesh.nv,
                                               sizeof(int32_t));
        if (materials == NULL) goto done;
    }
    if (support_path != NULL) {
        support = (unsigned char *)rv_load_sidecar(support_path, mesh.nv, 1);
        if (support == NULL) goto done;
    }

    memset(&src, 0, sizeof(src));
    if (source_path != NULL) {
        if (MeshBin_read_malloc(source_path, &src) != 0) {
            fprintf(stderr, "ERROR: cannot read source %s\n", source_path);
            goto done;
        }
        have_source = 1;
        fprintf(stderr, "  source: %zu verts\n", src.nv);
    }
    if (depth_runs_path != NULL) {
        depth_csv = fopen(depth_runs_path, "wb");
        if (depth_csv == NULL) {
            fprintf(stderr, "ERROR: cannot write %s\n", depth_runs_path);
            goto done;
        }
        fputs("run_id,orientation,run_length,edge_index,vertex_a,vertex_b,"
              "mesh_component,reconstruction_a,reconstruction_b,material_a,"
              "material_b,support_a,support_b,u_a,v_a,u_b,v_b,z_a,y_a,x_a,"
              "z_b,y_b,x_b,radius_a,radius_b,hp_a,hp_b,abs_hp_step,"
              "abs_radius_step\n", depth_csv);
    }
    arena = Arena_new();
    if (rv_measure(arena, &mesh, labels, materials, support,
                   have_source ? &src : NULL,
                   have_y && have_x,
                    umb_y, umb_x, pitch, core_radius, du, dv,
                    depth_csv, &rep) != 0) {
        fprintf(stderr, "ERROR: %s is empty or malformed\n", in_path);
        goto done;
    }
    fails = rv_gates(&rep, gates, &ngate);
    rv_print(&rep, gates, ngate, fails);
    if (report != NULL) {
        if (rv_write_report(report, in_path, &rep, gates, ngate, fails) != 0) {
            fprintf(stderr, "ERROR: cannot write %s\n", report);
            goto done;
        }
        fprintf(stderr, "  wrote %s\n", report);
    }
    if (depth_csv != NULL) {
        fclose(depth_csv);
        depth_csv = NULL;
        fprintf(stderr, "  wrote %s\n", depth_runs_path);
    }
    rc = fails == 0 ? 0 : 4;
done:
    if (depth_csv != NULL) fclose(depth_csv);
    if (arena != NULL) Arena_dispose(&arena);
    if (have_source) MeshBin_dispose(&src);
    free(labels);
    free(materials);
    free(support);
    MeshBin_dispose(&mesh);
    return rc;
}
