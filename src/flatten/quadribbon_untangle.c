/* ============================================================================
 * quadribbon_untangle.c -- see quadribbon_untangle.h.
 *
 * Extraction of solid_quad_ribbon.c's collision machinery (turn-order repair,
 * elastic-shell contact solve, settle-toward-rest), generalized from the
 * H x W structured lattice to any parameterized quadribbon:
 *
 *   lattice quad cell        -> occupied UV bucket (QrGrid), CSR vertex sets
 *   v = r*W + c neighbours   -> mesh vertex adjacency CSR
 *   face/2 -> cell           -> face centroid UV -> bucket
 *   red/black Gauss-Seidel   -> two-buffer Jacobi, sweeps doubled
 *   phase-based pitch        -> row-geometry inference (the original file's
 *                               own fallback, promoted to primary) with an
 *                               optional |pitch| prior from the config
 *
 * ACCEPTANCE GATES PERMANENTLY DISARMED (user directive 2026-08-30): every
 * trial is exact-audited, its gated verdict computed and logged verbatim
 * ("accept (gate-would-reject)" marks a trial the old machinery would have
 * refused), and then committed unconditionally.  The final state ships --
 * never a lex-best snapshot, never a rollback.  Rejections and failures are
 * research data judged by inspection (bakes + cross sections), not stop
 * conditions.  The enforcement code (preflight rejection, rollback ladders,
 * lexicographic partial publication, return-path adoption, transaction
 * rollback) is preserved under #if 0 -- deliberately NOT behind a runtime
 * flag, so it cannot re-arm by accident.
 *
 * Allocation is malloc/free (matching the source region); motion is radial
 * about an axis parallel to +Z in zyx coordinates.
 * ==========================================================================*/
#include "quadribbon_untangle.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../remesh/intersection_cleanup.h"

/* ---- exact-audit wall-clock accounting ------------------------------------ */

static size_t qr_audit_calls = 0;
static double qr_audit_ms = 0.0;

static void qr_audit_reset(void) { qr_audit_calls = 0; qr_audit_ms = 0.0; }

static void qr_audit_report(const char *tag)
{
    fprintf(stderr,
            "[audit] %s totals: calls=%zu wall=%.0f ms (mean %.0f ms)\n",
            tag, qr_audit_calls, qr_audit_ms,
            qr_audit_calls ? qr_audit_ms / (double)qr_audit_calls : 0.0);
}

/* ---- conflict scan (verbatim shape from the lattice version) -------------- */

typedef struct QrConflict {
    size_t face_lo, face_hi;
    double u_lo, u_hi;
    int hit_kind;
} QrConflict;

typedef struct QrScan {
    const float *uv;
    const int32_t *faces;
    double minimum_u_separation;
    size_t *local_face_degree;
    QrConflict *pair;
    size_t count, capacity;
    int allocation_failed;
} QrScan;

static double qr_face_mean_u(const float *uv, const int32_t *face)
{
    return ((double)uv[(size_t)face[0] * 2] + (double)uv[(size_t)face[1] * 2] +
            (double)uv[(size_t)face[2] * 2]) / 3.0;
}

static double qr_face_mean_radius(const float *verts, const int32_t *face,
                                  double axis_y, double axis_x)
{
    double sum = 0.0;
    for (int k = 0; k < 3; k++) {
        size_t v = (size_t)face[k];
        sum += hypot((double)verts[v * 3 + 1] - axis_y,
                     (double)verts[v * 3 + 2] - axis_x);
    }
    return sum / 3.0;
}

static int qr_conflict_visitor(size_t face_a, size_t face_b, int hit_kind,
                               void *context)
{
    QrScan *s = (QrScan *)context;
    const int32_t *fa = &s->faces[face_a * 3], *fb = &s->faces[face_b * 3];
    double ua = qr_face_mean_u(s->uv, fa), ub = qr_face_mean_u(s->uv, fb);
    if (fabs(ub - ua) < s->minimum_u_separation) {
        if (s->local_face_degree) {
            s->local_face_degree[face_a]++;
            s->local_face_degree[face_b]++;
        }
        return 0;
    }
    if (s->count == s->capacity) {
        size_t cap = s->capacity ? s->capacity * 2u : 1024u;
        QrConflict *p = NULL;
        if (cap < s->capacity || cap > SIZE_MAX / sizeof *s->pair) {
            s->allocation_failed = 1;
            return -1;
        }
        p = (QrConflict *)realloc(s->pair, cap * sizeof *s->pair);
        if (!p) { s->allocation_failed = 1; return -1; }
        s->pair = p;
        s->capacity = cap;
    }
    {
        QrConflict *p = &s->pair[s->count++];
        if (ua <= ub) {
            p->face_lo = face_a; p->face_hi = face_b;
            p->u_lo = ua; p->u_hi = ub;
        } else {
            p->face_lo = face_b; p->face_hi = face_a;
            p->u_lo = ub; p->u_hi = ua;
        }
        p->hit_kind = hit_kind;
    }
    return 0;
}

static int qr_audit(const float *verts, size_t nv, const int32_t *faces,
                    size_t nf, const float *uv, double minimum_u_separation,
                    size_t *face_degree, size_t *local_face_degree,
                    QrScan *scan, IntersectionCleanupStats *stats)
{
    IntersectionCleanupParams p;
    clock_t t0 = clock();
    int rc = 0;
    memset(scan, 0, sizeof *scan);
    scan->uv = uv;
    scan->faces = faces;
    scan->minimum_u_separation = minimum_u_separation;
    scan->local_face_degree = local_face_degree;
    if (local_face_degree)
        memset(local_face_degree, 0, nf * sizeof *local_face_degree);
    IntersectionCleanup_default_params(&p);
    p.gap_max = 1.0;
    p.parallel_angle_deg = 20.0;
    p.max_conflicts = nf > 1000000u ? nf : 1000000u;
    p.include_hinges = 0;
    rc = IntersectionCleanup_audit_visit_parallel(
        verts, nv, faces, nf, NULL, &p, face_degree, qr_conflict_visitor,
        scan, stats);
    {
        double ms = 1000.0 * (double)(clock() - t0) / (double)CLOCKS_PER_SEC;
        qr_audit_calls++;
        qr_audit_ms += ms;
        fprintf(stderr,
                "[audit] faces=%zu conflicts=%zu candidates=%zu ms=%.0f\n",
                nf, stats ? stats->conflicts : 0,
                stats ? stats->candidate_pairs : 0, ms);
    }
    if (rc != 0 || scan->allocation_failed) {
        fprintf(stderr,
                "[qr-untangle] exact audit incomplete: rc=%d "
                "retained-long=%zu allocation-failed=%d\n",
                rc, scan->count, scan->allocation_failed);
        free(scan->pair);
        memset(scan, 0, sizeof *scan);
        return -1;
    }
    return 0;
}

/* ---- QrGrid: rows, UV buckets, vertex adjacency --------------------------- */

typedef struct QrGrid {
    size_t nv, nf;
    const int32_t *faces;
    const float *uv;
    /* vertex adjacency CSR (undirected mesh edges) */
    int32_t *adj_off, *adj;
    /* UV buckets ("cells"): sorted occupied keys, CSR cell -> vertices */
    double bu, bv, u0, v0;
    int64_t *cell_key;
    size_t ncell;
    int32_t *cell_voff, *cell_verts;
    int32_t *vertex_cell;
    int32_t *face_cell;
    double grid_du;   /* median in-row u step (the lattice column width) */
} QrGrid;

static int64_t qr_cell_key(long vbin, long ubin)
{
    return ((int64_t)(vbin + (1L << 20)) << 32) |
           (int64_t)(ubin + (1L << 30));
}

static void qr_cell_bins(int64_t key, long *vbin, long *ubin)
{
    *vbin = (long)(key >> 32) - (1L << 20);
    *ubin = (long)(key & 0xffffffffLL) - (1L << 30);
}

static int qr_cmp_i64(const void *a, const void *b)
{
    int64_t x = *(const int64_t *)a, y = *(const int64_t *)b;
    return x < y ? -1 : x > y ? 1 : 0;
}

static int qr_cmp_dbl(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return x < y ? -1 : x > y ? 1 : 0;
}

static size_t qr_cell_find(const QrGrid *g, long vbin, long ubin)
{
    int64_t want = qr_cell_key(vbin, ubin);
    size_t lo = 0, hi = g->ncell;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (g->cell_key[mid] < want) lo = mid + 1; else hi = mid;
    }
    return lo < g->ncell && g->cell_key[lo] == want ? lo : SIZE_MAX;
}

static void qr_grid_free(QrGrid *g)
{
    free(g->adj_off); free(g->adj);
    free(g->cell_key); free(g->cell_voff); free(g->cell_verts);
    free(g->vertex_cell); free(g->face_cell);
    memset(g, 0, sizeof *g);
}

static int qr_grid_build(QrGrid *g, size_t nv, const int32_t *faces,
                         size_t nf, const float *uv)
{
    memset(g, 0, sizeof *g);
    g->nv = nv; g->nf = nf; g->faces = faces; g->uv = uv;

    /* --- vertex adjacency from face edges (deduped via sorted i64 keys) --- */
    {
        int64_t *ek = (int64_t *)malloc((nf * 3 + 1) * sizeof *ek);
        size_t nek = 0, nuniq = 0;
        if (!ek) return -1;
        for (size_t f = 0; f < nf; f++)
            for (int k = 0; k < 3; k++) {
                int32_t a = faces[f * 3 + k], b = faces[f * 3 + (k + 1) % 3];
                int64_t lo = a < b ? a : b, hi = a < b ? b : a;
                if (a < 0 || b < 0 || (size_t)a >= nv || (size_t)b >= nv ||
                    a == b)
                    continue;
                ek[nek++] = (lo << 31) | hi;
            }
        qsort(ek, nek, sizeof *ek, qr_cmp_i64);
        for (size_t k = 0; k < nek; k++)
            if (k == 0 || ek[k] != ek[nuniq - 1]) ek[nuniq++] = ek[k];
        g->adj_off = (int32_t *)calloc(nv + 1, sizeof *g->adj_off);
        g->adj = (int32_t *)malloc((2 * nuniq + 1) * sizeof *g->adj);
        if (!g->adj_off || !g->adj) { free(ek); qr_grid_free(g); return -1; }
        for (size_t k = 0; k < nuniq; k++) {
            g->adj_off[(ek[k] >> 31) + 1]++;
            g->adj_off[(ek[k] & 0x7fffffff) + 1]++;
        }
        for (size_t v = 0; v < nv; v++) g->adj_off[v + 1] += g->adj_off[v];
        {
            int32_t *cur = (int32_t *)malloc((nv + 1) * sizeof *cur);
            if (!cur) { free(ek); qr_grid_free(g); return -1; }
            memcpy(cur, g->adj_off, nv * sizeof *cur);
            for (size_t k = 0; k < nuniq; k++) {
                int32_t a = (int32_t)(ek[k] >> 31);
                int32_t b = (int32_t)(ek[k] & 0x7fffffff);
                g->adj[cur[a]++] = b;
                g->adj[cur[b]++] = a;
            }
            free(cur);
        }
        free(ek);
    }

    /* --- bucket sizes: median in-row u step and median row pitch ---------- */
    {
        size_t ns = nv < 200000 ? nv : 200000;
        double *du = (double *)malloc((ns + 1) * sizeof *du);
        double *dv = (double *)malloc((ns + 1) * sizeof *dv);
        size_t ndu = 0, ndv = 0;
        if (!du || !dv) { free(du); free(dv); qr_grid_free(g); return -1; }
        for (size_t v = 0; v < nv && (ndu < ns || ndv < ns); v++) {
            for (int32_t e = g->adj_off[v]; e < g->adj_off[v + 1]; e++) {
                size_t w = (size_t)g->adj[e];
                double ddu = fabs((double)uv[w * 2] - (double)uv[v * 2]);
                double ddv = fabs((double)uv[w * 2 + 1] -
                                  (double)uv[v * 2 + 1]);
                if (ddv < 1e-6 && ddu > 1e-6 && ndu < ns) du[ndu++] = ddu;
                if (ddv > 1e-6 && ndv < ns) dv[ndv++] = ddv;
            }
        }
        if (ndu < 8 || ndv < 8) { free(du); free(dv); qr_grid_free(g); return -1; }
        qsort(du, ndu, sizeof *du, qr_cmp_dbl);
        qsort(dv, ndv, sizeof *dv, qr_cmp_dbl);
        g->grid_du = du[ndu / 2];
        g->bu = g->grid_du > 0.25 ? g->grid_du : 0.25;
        g->bv = dv[ndv / 2] > 0.25 ? dv[ndv / 2] : 0.25;
        free(du); free(dv);
    }
    {
        double ulo = DBL_MAX, vlo = DBL_MAX;
        for (size_t v = 0; v < nv; v++) {
            if ((double)uv[v * 2] < ulo) ulo = (double)uv[v * 2];
            if ((double)uv[v * 2 + 1] < vlo) vlo = (double)uv[v * 2 + 1];
        }
        g->u0 = ulo; g->v0 = vlo;
    }

    /* --- occupied buckets, CSR cell->verts, vertex/face->cell ------------- */
    {
        int64_t *keys = (int64_t *)malloc((nv + 1) * sizeof *keys);
        g->vertex_cell = (int32_t *)malloc((nv + 1) * sizeof *g->vertex_cell);
        if (!keys || !g->vertex_cell) { free(keys); qr_grid_free(g); return -1; }
        for (size_t v = 0; v < nv; v++) {
            long ub = (long)floor(((double)uv[v * 2] - g->u0) / g->bu);
            long vb = (long)floor(((double)uv[v * 2 + 1] - g->v0) / g->bv +
                                  0.5);
            keys[v] = qr_cell_key(vb, ub);
        }
        g->cell_key = (int64_t *)malloc((nv + 1) * sizeof *g->cell_key);
        if (!g->cell_key) { free(keys); qr_grid_free(g); return -1; }
        memcpy(g->cell_key, keys, nv * sizeof *keys);
        qsort(g->cell_key, nv, sizeof *g->cell_key, qr_cmp_i64);
        {
            size_t n = 0;
            for (size_t k = 0; k < nv; k++)
                if (k == 0 || g->cell_key[k] != g->cell_key[n - 1])
                    g->cell_key[n++] = g->cell_key[k];
            g->ncell = n;
        }
        g->cell_voff = (int32_t *)calloc(g->ncell + 1, sizeof *g->cell_voff);
        g->cell_verts = (int32_t *)malloc((nv + 1) * sizeof *g->cell_verts);
        if (!g->cell_voff || !g->cell_verts) {
            free(keys); qr_grid_free(g); return -1;
        }
        for (size_t v = 0; v < nv; v++) {
            long vb, ub;
            size_t c;
            qr_cell_bins(keys[v], &vb, &ub);
            c = qr_cell_find(g, vb, ub);
            g->vertex_cell[v] = (int32_t)c;
            g->cell_voff[c + 1]++;
        }
        for (size_t c = 0; c < g->ncell; c++)
            g->cell_voff[c + 1] += g->cell_voff[c];
        {
            int32_t *cur = (int32_t *)malloc((g->ncell + 1) * sizeof *cur);
            if (!cur) { free(keys); qr_grid_free(g); return -1; }
            memcpy(cur, g->cell_voff, g->ncell * sizeof *cur);
            for (size_t v = 0; v < nv; v++)
                g->cell_verts[cur[g->vertex_cell[v]]++] = (int32_t)v;
            free(cur);
        }
        free(keys);
    }
    g->face_cell = (int32_t *)malloc((nf + 1) * sizeof *g->face_cell);
    if (!g->face_cell) { qr_grid_free(g); return -1; }
    for (size_t f = 0; f < nf; f++) {
        double mu = 0.0, mv = 0.0;
        long ub = 0, vb = 0;
        size_t c;
        for (int k = 0; k < 3; k++) {
            size_t v = (size_t)faces[f * 3 + k];
            mu += (double)uv[v * 2];
            mv += (double)uv[v * 2 + 1];
        }
        mu /= 3.0; mv /= 3.0;
        ub = (long)floor((mu - g->u0) / g->bu);
        vb = (long)floor((mv - g->v0) / g->bv + 0.5);
        c = qr_cell_find(g, vb, ub);
        if (c == SIZE_MAX) {
            /* centroid bucket unoccupied by any vertex: use a corner's */
            c = (size_t)g->vertex_cell[(size_t)faces[f * 3]];
        }
        g->face_cell[f] = (int32_t)c;
    }
    return 0;
}

/* ---- orientation preflight (transcribed; lattice feather -> adjacency) ---- */

static int qr_face_orientation_bad(const float *p, const float *anchor,
                                   size_t nv, const int32_t *face)
{
    int32_t ii[3] = { face[0], face[1], face[2] };
    size_t a, b, c;
    double oe0[3], oe1[3], ce0[3], ce1[3], on[3], cn[3], ref2, dot;
    if (ii[0] < 0 || ii[1] < 0 || ii[2] < 0 || (size_t)ii[0] >= nv ||
        (size_t)ii[1] >= nv || (size_t)ii[2] >= nv)
        return 1;
    a = (size_t)ii[0]; b = (size_t)ii[1]; c = (size_t)ii[2];
    for (int k = 0; k < 3; k++) {
        oe0[k] = (double)anchor[b * 3 + (size_t)k] - anchor[a * 3 + (size_t)k];
        oe1[k] = (double)anchor[c * 3 + (size_t)k] - anchor[a * 3 + (size_t)k];
        ce0[k] = (double)p[b * 3 + (size_t)k] - p[a * 3 + (size_t)k];
        ce1[k] = (double)p[c * 3 + (size_t)k] - p[a * 3 + (size_t)k];
    }
    on[0] = oe0[1] * oe1[2] - oe0[2] * oe1[1];
    on[1] = oe0[2] * oe1[0] - oe0[0] * oe1[2];
    on[2] = oe0[0] * oe1[1] - oe0[1] * oe1[0];
    cn[0] = ce0[1] * ce1[2] - ce0[2] * ce1[1];
    cn[1] = ce0[2] * ce1[0] - ce0[0] * ce1[2];
    cn[2] = ce0[0] * ce1[1] - ce0[1] * ce1[0];
    ref2 = on[0] * on[0] + on[1] * on[1] + on[2] * on[2];
    dot = on[0] * cn[0] + on[1] * cn[1] + on[2] * cn[2];
    return !(ref2 > 1e-18) || !(dot > 0.05 * ref2);
}

/* Faces already degenerate or inverted in the transaction/round input can
 * never be repaired by zeroing candidate motion: their identity trial fails
 * the same preflight, every alpha dies on them, and the whole transaction
 * aborts with "no exact-audit-safe response remained" (full-sheet v120).
 * Flag them once per round against the round input; candidate gating then
 * reads "no NEW orientation violations". */
static size_t qr_orientation_baseline(const float *anchor, size_t nv,
                                      const int32_t *faces, size_t nf,
                                      uint8_t *input_bad)
{
    size_t bad = 0;
    ptrdiff_t f;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:bad) schedule(static)
#endif
    for (f = 0; f < (ptrdiff_t)nf; f++) {
        uint8_t b = (uint8_t)qr_face_orientation_bad(
            anchor, anchor, nv, &faces[(size_t)f * 3]);
        input_bad[(size_t)f] = b;
        bad += b;
    }
    return bad;
}

static size_t qr_orientation_violations_list(
    const float *p, const float *anchor, size_t nv, const int32_t *faces,
    const int32_t *face_list, size_t nfl, const uint8_t *input_bad,
    float *motion_scale, size_t *newly_zeroed)
{
    size_t bad = 0, nnew = 0;
    if (!motion_scale) {
        ptrdiff_t i;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:bad) schedule(static)
#endif
        for (i = 0; i < (ptrdiff_t)nfl; i++) {
            size_t f = (size_t)face_list[i];
            if (input_bad && input_bad[f]) continue;
            bad += (size_t)qr_face_orientation_bad(p, anchor, nv,
                                                   &faces[f * 3]);
        }
        if (newly_zeroed) *newly_zeroed = 0;
        return bad;
    }
    for (size_t i = 0; i < nfl; i++) {
        size_t f = (size_t)face_list[i];
        if (input_bad && input_bad[f]) continue;
        if (qr_face_orientation_bad(p, anchor, nv, &faces[f * 3])) {
            const int32_t *ii = &faces[f * 3];
            bad++;
            for (int k = 0; k < 3; k++) {
                size_t v = (size_t)ii[k];
                if (motion_scale[v] > 0.0f) { motion_scale[v] = 0.0f; nnew++; }
            }
        }
    }
    if (newly_zeroed) *newly_zeroed = nnew;
    return bad;
}

static size_t qr_orientation_violations(
    const float *p, const float *anchor, size_t nv, const int32_t *faces,
    size_t nf, const uint8_t *input_bad, float *motion_scale,
    size_t *newly_zeroed)
{
    size_t bad = 0, nnew = 0;
    if (!motion_scale) {
        ptrdiff_t f;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:bad) schedule(static)
#endif
        for (f = 0; f < (ptrdiff_t)nf; f++) {
            if (input_bad && input_bad[(size_t)f]) continue;
            bad += (size_t)qr_face_orientation_bad(p, anchor, nv,
                                                   &faces[(size_t)f * 3]);
        }
        if (newly_zeroed) *newly_zeroed = 0;
        return bad;
    }
    for (size_t f = 0; f < nf; f++) {
        if (input_bad && input_bad[f]) continue;
        if (qr_face_orientation_bad(p, anchor, nv, &faces[f * 3])) {
            const int32_t *ii = &faces[f * 3];
            bad++;
            for (int k = 0; k < 3; k++) {
                size_t v = (size_t)ii[k];
                if (motion_scale[v] > 0.0f) { motion_scale[v] = 0.0f; nnew++; }
            }
        }
    }
    if (newly_zeroed) *newly_zeroed = nnew;
    return bad;
}

/* min-plus dilation of the freeze mask over the vertex graph (the lattice
 * version dilated over the 4-neighbour grid); +0.125 per ring, 8 rings */
static void qr_feather_motion_scale(const QrGrid *g, float *scale,
                                    float *temporary, int rings)
{
    for (int ring = 0; ring < rings; ring++) {
        memcpy(temporary, scale, g->nv * sizeof *temporary);
        ptrdiff_t v;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (v = 0; v < (ptrdiff_t)g->nv; v++) {
            float a = scale[(size_t)v];
            for (int32_t e = g->adj_off[v]; e < g->adj_off[v + 1]; e++) {
                float b = scale[(size_t)g->adj[e]] + 0.125f;
                if (b < a) a = b;
            }
            temporary[(size_t)v] = a < 1.0f ? a : 1.0f;
        }
        memcpy(scale, temporary, g->nv * sizeof *scale);
    }
}

/* ---- winding pitch inference (row geometry; original fallback promoted) --- */

static int qr_infer_pitch(const QrGrid *g, const float *verts,
                          const uint8_t *movable, double axis_y, double axis_x,
                          const QrScan *scan, double pitch_hint,
                          double *out_pitch)
{
    size_t cap = scan->count > g->nv ? scan->count : g->nv;
    double *span = NULL;
    double *sample = NULL;
    size_t ns = 0;
    double du_period = 0.0;
    /* The SIGN of dr/du is a global property of the parameterization (the
     * spiral's handedness).  The row-hop sampler below measures it by
     * comparing single vertices one conflict-period apart in u, which
     * crosses island-placement gaps -- on the fragmented 4x5x5 soup that
     * gave 57% sign agreement (a coin flip) and silently disarmed the whole
     * untangler.  Measure the sign instead as the least-squares slope of r
     * against u over ALL mesh edges: within-wrap edges never cross
     * placement gaps, and millions of them average the +-2 vox surface
     * roughness away.  The t-statistic gates confidence; the magnitude
     * comes from the config prior.  The row-hop path remains the hint-less
     * fallback. */
    if (pitch_hint > 0.0) {
        double suu = 0.0, sur = 0.0, srr = 0.0;
        size_t nedge = 0;
        for (size_t v = 0; v < g->nv; v++) {
            double uvv = (double)g->uv[v * 2];
            double rv = hypot((double)verts[v * 3 + 1] - axis_y,
                              (double)verts[v * 3 + 2] - axis_x);
            if (!isfinite(rv)) continue;
            for (int32_t e = g->adj_off[v]; e < g->adj_off[v + 1]; e++) {
                size_t w = (size_t)g->adj[e];
                double du = 0.0, dr = 0.0, rw = 0.0;
                if (w <= v) continue;          /* each undirected edge once */
                rw = hypot((double)verts[w * 3 + 1] - axis_y,
                           (double)verts[w * 3 + 2] - axis_x);
                if (!isfinite(rw)) continue;
                du = (double)g->uv[w * 2] - uvv;
                if (fabs(du) > 8.0 * g->grid_du) continue; /* not lattice-local */
                dr = rw - rv;
                suu += du * du;
                sur += du * dr;
                srr += dr * dr;
                nedge++;
            }
        }
        if (suu > 0.0 && srr > 0.0 && nedge >= 1024) {
            double corr = sur / sqrt(suu * srr);
            double tstat = corr * sqrt((double)nedge);
            fprintf(stderr,
                    "[qr-untangle] edge-slope pitch sign: slope=%.5f "
                    "corr=%.4f t=%.1f edges=%zu -> pitch=%.3f\n",
                    sur / suu, corr, tstat, nedge,
                    (sur < 0.0 ? -pitch_hint : pitch_hint));
            if (fabs(tstat) >= 6.0) {
                *out_pitch = sur < 0.0 ? -pitch_hint : pitch_hint;
                return 0;
            }
        }
        fprintf(stderr,
                "[qr-untangle] edge-slope sign inconclusive; trying "
                "row-hop sampler\n");
    }
    span = (double *)malloc((cap + 1) * sizeof *span);
    sample = (double *)malloc((cap + 1) * sizeof *sample);
    if (!span || !sample) { free(span); free(sample); return -1; }
    /* median conflict u-period */
    for (size_t q = 0; q < scan->count; q++)
        span[q] = scan->pair[q].u_hi - scan->pair[q].u_lo;
    qsort(span, scan->count, sizeof *span, qr_cmp_dbl);
    du_period = span[scan->count / 2];
    if (!(du_period > 2.0 * g->grid_du)) du_period = 2.0 * g->grid_du;
    /* same-row geometry at that period: rows are the v buckets, vertices
     * within a cell row scanned via the cell CSR (cells are sorted by
     * (vbin, ubin), so consecutive cells of one row are adjacent) */
    for (size_t c = 0; c < g->ncell && ns + 4 < cap; c++) {
        long vb, ub, vb2, ub2;
        size_t c2;
        qr_cell_bins(g->cell_key[c], &vb, &ub);
        c2 = qr_cell_find(g, vb, ub + (long)lround(du_period / g->bu));
        if (c2 == SIZE_MAX) continue;
        qr_cell_bins(g->cell_key[c2], &vb2, &ub2);
        (void)vb2; (void)ub2;
        {
            size_t a = (size_t)g->cell_verts[g->cell_voff[c]];
            size_t b = (size_t)g->cell_verts[g->cell_voff[c2]];
            double ra, rb;
            if (movable && (movable[a] || movable[b])) continue;
            ra = hypot((double)verts[a * 3 + 1] - axis_y,
                       (double)verts[a * 3 + 2] - axis_x);
            rb = hypot((double)verts[b * 3 + 1] - axis_y,
                       (double)verts[b * 3 + 2] - axis_x);
            if (isfinite(ra) && isfinite(rb)) sample[ns++] = rb - ra;
        }
    }
    if (ns < 16 && movable != NULL) {
        /* movable-gated pass found too few anchored pairs; retry ungated
         * (the whole initializer, logged as such below) */
        ns = 0;
        for (size_t c = 0; c < g->ncell && ns + 4 < cap; c++) {
            long vb, ub;
            size_t c2;
            qr_cell_bins(g->cell_key[c], &vb, &ub);
            c2 = qr_cell_find(g, vb, ub + (long)lround(du_period / g->bu));
            if (c2 == SIZE_MAX) continue;
            {
                size_t a = (size_t)g->cell_verts[g->cell_voff[c]];
                size_t b = (size_t)g->cell_verts[g->cell_voff[c2]];
                double ra = hypot((double)verts[a * 3 + 1] - axis_y,
                                  (double)verts[a * 3 + 2] - axis_x);
                double rb = hypot((double)verts[b * 3 + 1] - axis_y,
                                  (double)verts[b * 3 + 2] - axis_x);
                if (isfinite(ra) && isfinite(rb)) sample[ns++] = rb - ra;
            }
        }
    }
    if (ns < 16) { free(span); free(sample); return -1; }
    qsort(sample, ns, sizeof *sample, qr_cmp_dbl);
    {
        double pitch = sample[ns / 2];
        size_t agrees = 0;
        double agreement;
        for (size_t k = 0; k < ns; k++) agrees += (sample[k] * pitch) > 0.0;
        agreement = (double)agrees / (double)ns;
        fprintf(stderr,
                "[qr-untangle] inferred u-period=%.1f vox; pitch "
                "p10/p50/p90=%.3f/%.3f/%.3f vox, sign agreement %.1f%%, "
                "samples=%zu\n",
                du_period, sample[ns / 10], pitch, sample[(ns * 9) / 10],
                100.0 * agreement, ns);
        if (!isfinite(pitch) || agreement < 0.60) {
            free(span); free(sample);
            return -1;
        }
        if (fabs(pitch) < 2.0) {
            if (pitch_hint > 0.0) {
                /* magnitude from the config prior, sign from the data */
                pitch = pitch >= 0.0 ? pitch_hint : -pitch_hint;
                fprintf(stderr,
                        "[qr-untangle] measured |pitch| ambiguous; using "
                        "config magnitude with measured sign: %.3f\n", pitch);
            } else {
                free(span); free(sample);
                return -1;
            }
        }
        *out_pitch = pitch;
    }
    free(span);
    free(sample);
    return 0;
}

/* ---- turn-order repair (pass 1) ------------------------------------------- */

static size_t qr_mark_face_seed(const int32_t *face, const uint8_t *movable,
                                double value, double weight,
                                double *seed_sum, double *seed_weight)
{
    size_t n = 0;
    for (int k = 0; k < 3; k++) n += (size_t)(movable[(size_t)face[k]] != 0);
    if (!n) return 0;
    for (int k = 0; k < 3; k++) {
        size_t v = (size_t)face[k];
        if (!movable[v]) continue;
        seed_sum[v] += weight * value;
        seed_weight[v] += weight;
    }
    return n;
}

static void qr_dilate_active(const QrGrid *g, uint8_t *active,
                             const uint8_t *movable, int rings)
{
    int32_t *frontier = (int32_t *)malloc((g->nv + 1) * sizeof *frontier);
    int32_t *next = (int32_t *)malloc((g->nv + 1) * sizeof *next);
    size_t nfrontier = 0;
    if (!frontier || !next) { free(frontier); free(next); return; }
    for (size_t v = 0; v < g->nv; v++)
        if (active[v]) frontier[nfrontier++] = (int32_t)v;
    for (int ring = 0; ring < rings && nfrontier; ring++) {
        size_t nnext = 0;
        for (size_t i = 0; i < nfrontier; i++) {
            size_t v = (size_t)frontier[i];
            for (int32_t e = g->adj_off[v]; e < g->adj_off[v + 1]; e++) {
                size_t w = (size_t)g->adj[e];
                if (movable[w] && !active[w]) {
                    active[w] = 1;
                    next[nnext++] = (int32_t)w;
                }
            }
        }
        {
            int32_t *swap = frontier;
            frontier = next;
            next = swap;
        }
        nfrontier = nnext;
    }
    free(frontier);
    free(next);
}

static int qr_build_correction(const QrGrid *g, const float *verts,
                               const uint8_t *movable, double axis_y,
                               double axis_x, const QrScan *scan, double pitch,
                               double desired_gap, float *correction,
                               size_t *out_seed, size_t *out_active,
                               int *out_iterations, double *out_max)
{
    size_t nv = g->nv;
    double *seed_sum = (double *)calloc(nv, sizeof *seed_sum);
    double *seed_weight = (double *)calloc(nv, sizeof *seed_weight);
    uint8_t *active = (uint8_t *)calloc(nv, 1);
    float *next = (float *)calloc(nv, sizeof *next);
    const double order = pitch >= 0.0 ? 1.0 : -1.0;
    const double deficit_cap = fmax(desired_gap, 1.5 * fabs(pitch));
    size_t seeds = 0, nactive = 0;
    int iterations = 0;
    if (!seed_sum || !seed_weight || !active || !next) {
        free(seed_sum); free(seed_weight); free(active); free(next);
        return -1;
    }
    for (size_t q = 0; q < scan->count; q++) {
        const int32_t *lo = &g->faces[scan->pair[q].face_lo * 3];
        const int32_t *hi = &g->faces[scan->pair[q].face_hi * 3];
        double rlo = qr_face_mean_radius(verts, lo, axis_y, axis_x);
        double rhi = qr_face_mean_radius(verts, hi, axis_y, axis_x);
        double deficit = desired_gap - order * (rhi - rlo);
        size_t nlo = 0, nhi = 0;
        double lo_share, hi_share, severity;
        if (!(deficit > 0.0)) continue;
        if (deficit > deficit_cap) deficit = deficit_cap;
        for (int k = 0; k < 3; k++) {
            nlo += (size_t)(movable[(size_t)lo[k]] != 0);
            nhi += (size_t)(movable[(size_t)hi[k]] != 0);
        }
        if (!nlo && !nhi) continue;
        lo_share = nlo ? (nhi ? 0.5 : 1.0) : 0.0;
        hi_share = nhi ? (nlo ? 0.5 : 1.0) : 0.0;
        severity = 1.0 + deficit / fmax(desired_gap, 1.0);
        if (lo_share)
            qr_mark_face_seed(lo, movable, -order * deficit * lo_share,
                              severity, seed_sum, seed_weight);
        if (hi_share)
            qr_mark_face_seed(hi, movable, order * deficit * hi_share,
                              severity, seed_sum, seed_weight);
    }
    for (size_t v = 0; v < nv; v++)
        if (seed_weight[v] > 0.0) { active[v] = 1; seeds++; }
    if (!seeds) {
        free(seed_sum); free(seed_weight); free(active); free(next);
        return 1;
    }
    qr_dilate_active(g, active, movable, 32);
    for (size_t v = 0; v < nv; v++) nactive += active[v] != 0;
    memset(correction, 0, nv * sizeof *correction);
    {
        const double screen = 0.002;
        /* Jacobi over the vertex graph replaces the lattice red/black GS;
         * iteration budget doubled to compensate its slower mixing */
        for (iterations = 0; iterations < 768; iterations++) {
            double maximum_change = 0.0;
            ptrdiff_t v;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (v = 0; v < (ptrdiff_t)nv; v++) {
                double sum = 0.0, degree = 0.0, sw, target, value;
                if (!active[(size_t)v] || !movable[(size_t)v]) {
                    next[(size_t)v] = correction[(size_t)v];
                    continue;
                }
                for (int32_t e = g->adj_off[v]; e < g->adj_off[v + 1]; e++) {
                    sum += (double)correction[(size_t)g->adj[e]];
                    degree += 1.0;
                }
                sw = seed_weight[(size_t)v] > 0.0
                   ? 20.0 * fmin(seed_weight[(size_t)v], 5.0) : 0.0;
                target = seed_weight[(size_t)v] > 0.0
                       ? seed_sum[(size_t)v] / seed_weight[(size_t)v] : 0.0;
                value = (sum + sw * target) / (degree + screen + sw);
                next[(size_t)v] = (float)value;
            }
            for (size_t w = 0; w < nv; w++) {
                double change = fabs((double)next[w] - (double)correction[w]);
                if (change > maximum_change) maximum_change = change;
                correction[w] = next[w];
            }
            if (maximum_change <= 1e-5) { iterations++; break; }
        }
    }
    {
        double mx = 0.0;
        for (size_t v = 0; v < nv; v++)
            if (fabs((double)correction[v]) > mx)
                mx = fabs((double)correction[v]);
        *out_seed = seeds;
        *out_active = nactive;
        *out_iterations = iterations;
        *out_max = mx;
    }
    free(seed_sum); free(seed_weight); free(active); free(next);
    return 0;
}

static void qr_apply_radial(float *out, const float *base, size_t nv,
                            const uint8_t *movable, const float *correction,
                            const float *keep, double scale, double axis_y,
                            double axis_x, double *out_rms, double *out_max)
{
    double ss = 0.0, mx = 0.0;
    size_t moved = 0;
    memcpy(out, base, nv * 3 * sizeof *out);
    for (size_t v = 0; v < nv; v++) {
        double dy, dx, retained, radius, d, factor;
        if (!movable[v] || correction[v] == 0.0f) continue;
        dy = (double)base[v * 3 + 1] - axis_y;
        dx = (double)base[v * 3 + 2] - axis_x;
        retained = keep ? (double)keep[v] : 1.0;
        radius = hypot(dy, dx);
        d = scale * retained * (double)correction[v];
        if (!(radius > 1e-9) || radius + d <= 1e-6) continue;
        factor = (radius + d) / radius;
        out[v * 3 + 1] = (float)(axis_y + factor * dy);
        out[v * 3 + 2] = (float)(axis_x + factor * dx);
        ss += d * d;
        if (fabs(d) > mx) mx = fabs(d);
        moved++;
    }
    *out_rms = moved ? sqrt(ss / (double)moved) : 0.0;
    *out_max = mx;
}

static int qr_repair_turn_order(const QrGrid *g, float *verts,
                                const uint8_t *movable, double axis_y,
                                double axis_x, double minimum_u_separation,
                                double pitch_hint,
                                QuadribbonUntangleStats *s)
{
    size_t nv = g->nv, nf = g->nf;
    const int32_t *faces = g->faces;
    const float *uv = g->uv;
    float *base = (float *)malloc(nv * 3 * sizeof *base);
    float *candidate = (float *)malloc(nv * 3 * sizeof *candidate);
    float *correction = (float *)malloc(nv * sizeof *correction);
    float *keep = (float *)malloc(nv * sizeof *keep);
    float *keep_tmp = (float *)malloc(nv * sizeof *keep_tmp);
    size_t *local_before = (size_t *)calloc(nf, sizeof *local_before);
    size_t *local_after = (size_t *)calloc(nf, sizeof *local_after);
    uint8_t *orientation_input_bad = (uint8_t *)malloc(nf);
    int rc = -1;
    if (!base || !candidate || !correction || !keep || !keep_tmp ||
        !local_before || !local_after || !orientation_input_bad)
        goto done;
    qr_audit_reset();
    for (int round = 0; round < 10; round++) {
        QrScan before_scan;
        IntersectionCleanupStats before;
        double pitch = 0.0, desired_gap;
        size_t seeds = 0, active = 0;
        int solve_iterations = 0, crc, accepted = 0;
        double correction_max = 0.0;
        memset(&before, 0, sizeof before);
        if (qr_audit(verts, nv, faces, nf, uv, minimum_u_separation, NULL,
                     local_before, &before_scan, &before) != 0)
            goto done;
        if (round == 0) {
            s->turn_order_input_long = before_scan.count;
        }
        fprintf(stderr,
                "[turn-order] round %d exact conflicts=%zu (long-range=%zu), "
                "folds=%zu\n",
                round + 1, before.conflicts, before_scan.count,
                before.fold_pairs);
        s->turn_order_output_long = before_scan.count;
        if (before_scan.count == 0) { free(before_scan.pair); rc = 0; goto done; }
        if (qr_infer_pitch(g, verts, movable, axis_y, axis_x, &before_scan,
                           pitch_hint, &pitch) != 0) {
            fprintf(stderr,
                    "[turn-order] refusing ambiguous winding direction/pitch\n");
            free(before_scan.pair);
            rc = 0;
            goto done;
        }
        desired_gap = fmax(3.0, fmin(16.0, 0.75 * fabs(pitch)));
        crc = qr_build_correction(g, verts, movable, axis_y, axis_x,
                                  &before_scan, pitch, desired_gap, correction,
                                  &seeds, &active, &solve_iterations,
                                  &correction_max);
        if (crc != 0) {
            if (crc > 0)
                fprintf(stderr,
                        "[turn-order] no movable vertices in long-range pairs\n");
            free(before_scan.pair);
            rc = crc < 0 ? -1 : 0;
            goto done;
        }
        fprintf(stderr,
                "[turn-order] sparse radial solve: desired gap %.3f, "
                "seeds=%zu, active=%zu/%zu, iterations=%d, max correction "
                "%.3f vox\n",
                desired_gap, seeds, active, nv, solve_iterations,
                correction_max);
        memcpy(base, verts, nv * 3 * sizeof *base);
        {
            size_t input_bad_faces = qr_orientation_baseline(
                base, nv, faces, nf, orientation_input_bad);
            if (input_bad_faces)
                fprintf(stderr,
                        "[turn-order] input orientation-bad faces=%zu "
                        "excluded from candidate preflight\n",
                        input_bad_faces);
        }
        /* Dense intersection fronts often need a sub-quarter transaction once
         * the large easy crossings have left the sheet.  Try those scales
         * before escalating; a rejected trial is read-only and
         * exact-audited. */
        {
            static const double trial_scale_dense[7] =
                { 0.25, 0.125, 0.0625, 0.375, 0.5, 0.75, 1.0 };
            static const double trial_scale_sparse[7] =
                { 0.375, 0.5, 0.25, 0.125, 0.0625, 0.75, 1.0 };
            const double *trial_scale = before_scan.count > 3000
                                      ? trial_scale_dense : trial_scale_sparse;
            const int maximum_local_rollbacks = 12;
            for (int trial = 0; trial < 7; trial++) {
                for (size_t v = 0; v < nv; v++) keep[v] = 1.0f;
                for (int rollback = 0; rollback <= maximum_local_rollbacks;
                     rollback++) {
                    double move_rms = 0.0, move_max = 0.0;
                    size_t orientation_bad;
                    QrScan after_scan;
                    IntersectionCleanupStats after;
                    size_t before_local, after_local;
                    int safe;
                    qr_apply_radial(candidate, base, nv, movable, correction,
                                    keep, trial_scale[trial], axis_y, axis_x,
                                    &move_rms, &move_max);
                    orientation_bad = qr_orientation_violations(
                        candidate, base, nv, faces, nf, orientation_input_bad,
                        NULL, NULL);
                    /* GATES DISARMED (user directive 2026-08-30): every
                     * acceptance predicate in this module is now
                     * measurement-only -- computed and logged in full, never
                     * enforced.  The algorithm's raw effect is judged by
                     * inspection (bakes + cross sections); rejections and
                     * failures are research data, not stop conditions.  The
                     * enforcement code is preserved under #if 0 -- NOT behind
                     * a runtime flag -- so it cannot re-arm by accident. */
                    if (orientation_bad)
                        fprintf(stderr,
                                "[turn-order] preflight measurement alpha "
                                "%.3f: orientation violations=%zu (gate "
                                "disarmed, proceeding to audit)\n",
                                trial_scale[trial], orientation_bad);
#if 0               /* disarmed: preflight rejection + feathered rollback */
                    if (orientation_bad) {
                        size_t orientation_marked = 0;
                        fprintf(stderr,
                                "[turn-order] reject trial alpha %.3f "
                                "rollback %d before exact audit: orientation "
                                "violations=%zu\n",
                                trial_scale[trial], rollback, orientation_bad);
                        if (rollback == maximum_local_rollbacks) break;
                        (void)qr_orientation_violations(
                            candidate, base, nv, faces, nf,
                            orientation_input_bad, keep, &orientation_marked);
                        if (orientation_marked == 0) break;
                        qr_feather_motion_scale(g, keep, keep_tmp, 8);
                        fprintf(stderr,
                                "[turn-order] preflight rolled back %zu "
                                "orientation vertices and feathered 8 rings\n",
                                orientation_marked);
                        continue;
                    }
#endif
                    memset(&after, 0, sizeof after);
                    if (qr_audit(candidate, nv, faces, nf, uv,
                                 minimum_u_separation, NULL, local_after,
                                 &after_scan, &after) != 0) {
                        free(before_scan.pair);
                        goto done;
                    }
                    before_local = before.conflicts - before_scan.count;
                    after_local = after.conflicts - after_scan.count;
                    safe = after.conflicts < before.conflicts &&
                           after_scan.count < before_scan.count &&
                           after_local <= before_local &&
                           after.fold_pairs <= before.fold_pairs &&
                           orientation_bad == 0;
                    fprintf(stderr,
                            "[turn-order] %s trial alpha %.3f rollback %d: "
                            "move rms/max %.3f/%.3f; conflicts %zu->%zu, long "
                            "%zu->%zu, local %zu->%zu, folds %zu->%zu\n",
                            safe ? "accept" : "accept (gate-would-reject)",
                            trial_scale[trial],
                            rollback, move_rms, move_max, before.conflicts,
                            after.conflicts, before_scan.count,
                            after_scan.count, before_local, after_local,
                            before.fold_pairs, after.fold_pairs);
                    /* GATE DISARMED: commit every exact-audited trial. */
                    memcpy(verts, candidate, nv * 3 * sizeof *verts);
                    accepted = 1;
                    s->turn_order_rounds++;
                    {
                        size_t *tmp = local_before;
                        local_before = local_after;
                        local_after = tmp;
                    }
                    free(after_scan.pair);
                    break;
#if 0               /* disarmed: lexicographic acceptance + rollback ladder */
                    if (safe) {
                        memcpy(verts, candidate, nv * 3 * sizeof *verts);
                        accepted = 1;
                        s->turn_order_rounds++;
                        {
                            size_t *tmp = local_before;
                            local_before = local_after;
                            local_after = tmp;
                        }
                        free(after_scan.pair);
                        break;
                    }
                    if (rollback == maximum_local_rollbacks ||
                        after.conflicts >= before.conflicts ||
                        after_scan.count >= before_scan.count) {
                        free(after_scan.pair);
                        break;
                    }
                    {
                        size_t marked = 0, orientation_marked = 0;
                        for (size_t f = 0; f < nf; f++)
                            if (local_after[f] > local_before[f])
                                for (int k = 0; k < 3; k++) {
                                    size_t v = (size_t)faces[f * 3 + (size_t)k];
                                    if (keep[v] > 0.0f) {
                                        keep[v] = 0.0f;
                                        marked++;
                                    }
                                }
                        (void)qr_orientation_violations(
                            candidate, base, nv, faces, nf,
                            orientation_input_bad, keep, &orientation_marked);
                        marked += orientation_marked;
                        free(after_scan.pair);
                        if (marked == 0) break;
                        qr_feather_motion_scale(g, keep, keep_tmp, 8);
                        fprintf(stderr,
                                "[turn-order] locally rolled back %zu "
                                "seed/boundary vertices and feathered 8 "
                                "rings\n",
                                marked);
                    }
#endif
                }
                if (accepted) break;
            }
        }
        free(before_scan.pair);
        if (!accepted) { rc = 0; goto done; }
    }
    rc = 0;
done:
    if (rc == 0) audit_totals: qr_audit_report("turn-order");
    (void)0;
    free(base); free(candidate); free(correction); free(keep); free(keep_tmp);
    free(local_before); free(local_after); free(orientation_input_bad);
    return rc;
}

/* ---- elastic-shell contact solve (pass 2) --------------------------------- */

typedef struct QrContact {
    size_t cell_lo, cell_hi;
} QrContact;

static int qr_contact_cmp(const void *aa, const void *bb)
{
    const QrContact *a = (const QrContact *)aa;
    const QrContact *b = (const QrContact *)bb;
    if (a->cell_lo < b->cell_lo) return -1;
    if (a->cell_lo > b->cell_lo) return 1;
    if (a->cell_hi < b->cell_hi) return -1;
    if (a->cell_hi > b->cell_hi) return 1;
    return 0;
}

static int qr_cells_within_collar(const QrGrid *g, const float *verts,
                                  size_t ca, size_t cb, double collar)
{
    double distance2 = 0.0;
    for (int k = 0; k < 3; k++) {
        double alo = DBL_MAX, ahi = -DBL_MAX, blo = DBL_MAX, bhi = -DBL_MAX;
        for (int32_t i = g->cell_voff[ca]; i < g->cell_voff[ca + 1]; i++) {
            double a = (double)verts[(size_t)g->cell_verts[i] * 3 + (size_t)k];
            if (a < alo) alo = a;
            if (a > ahi) ahi = a;
        }
        for (int32_t i = g->cell_voff[cb]; i < g->cell_voff[cb + 1]; i++) {
            double b = (double)verts[(size_t)g->cell_verts[i] * 3 + (size_t)k];
            if (b < blo) blo = b;
            if (b > bhi) bhi = b;
        }
        {
            double gap = ahi < blo ? blo - ahi : bhi < alo ? alo - bhi : 0.0;
            distance2 += gap * gap;
        }
    }
    return distance2 <= collar * collar;
}

/* Inflate every exact contact by a small UV neighbourhood, then retain only
 * cell pairs whose 3-D AABBs lie inside the speculative collar.  Collision
 * engines do this before response: without the collar an intersection can
 * simply migrate from the constrained pair to its unconstrained neighbour.
 * The two sides use the same row offset because param V is axial world Z. */
static QrContact *qr_build_contacts(const QrGrid *g, const float *verts,
                                    const QrScan *scan, double collar,
                                    int inflate, size_t *out_count)
{
    size_t cap, count = 0;
    QrContact *pair = NULL;
    if (out_count) *out_count = 0;
    if (!scan || scan->count == 0 || inflate < 1) return NULL;
    {
        size_t reserve_multiplier = scan->count < 1024 ? 64u : 8u;
        if (scan->count > SIZE_MAX / reserve_multiplier) return NULL;
        cap = scan->count * reserve_multiplier;
        if (cap < 1024) cap = 1024;
    }
    pair = (QrContact *)malloc(cap * sizeof *pair);
    if (!pair) return NULL;
    for (size_t q = 0; q < scan->count; q++) {
        size_t qlo = (size_t)g->face_cell[scan->pair[q].face_lo];
        size_t qhi = (size_t)g->face_cell[scan->pair[q].face_hi];
        long rlo, clo, rhi, chi;
        qr_cell_bins(g->cell_key[qlo], &rlo, &clo);
        qr_cell_bins(g->cell_key[qhi], &rhi, &chi);
        for (int dr = -inflate; dr <= inflate; dr++) {
            for (int dl = -inflate; dl <= inflate; dl++)
                for (int dh = -inflate; dh <= inflate; dh++) {
                    size_t cell_lo = qr_cell_find(g, rlo + dr, clo + dl);
                    size_t cell_hi = qr_cell_find(g, rhi + dr, chi + dh);
                    if (cell_lo == SIZE_MAX || cell_hi == SIZE_MAX) continue;
                    if (!qr_cells_within_collar(g, verts, cell_lo, cell_hi,
                                                collar))
                        continue;
                    if (count == cap) {
                        size_t nxt = cap * 2u;
                        QrContact *grown;
                        if (nxt < cap || nxt > SIZE_MAX / sizeof *pair) {
                            free(pair);
                            return NULL;
                        }
                        grown = (QrContact *)realloc(pair, nxt * sizeof *pair);
                        if (!grown) { free(pair); return NULL; }
                        pair = grown;
                        cap = nxt;
                    }
                    pair[count].cell_lo = cell_lo;
                    pair[count].cell_hi = cell_hi;
                    count++;
                }
        }
    }
    if (!count) { free(pair); return NULL; }
    qsort(pair, count, sizeof *pair, qr_contact_cmp);
    {
        size_t unique = 0;
        for (size_t q = 0; q < count; q++) {
            if (unique && pair[q].cell_lo == pair[unique - 1].cell_lo &&
                pair[q].cell_hi == pair[unique - 1].cell_hi)
                continue;
            pair[unique++] = pair[q];
        }
        if (out_count) *out_count = unique;
    }
    return pair;
}

static int qr_merge_ledger(const QrContact *old_contact, size_t nold,
                           const QrContact *new_contact, size_t nnew,
                           QrContact **out_contact, size_t *out_count)
{
    QrContact *merged = NULL;
    size_t i = 0, j = 0, n = 0;
    if (!out_contact || !out_count || (nold && !old_contact) ||
        (nnew && !new_contact) || nold > SIZE_MAX - nnew ||
        (nold + nnew) > SIZE_MAX / sizeof(QrContact))
        return -1;
    merged = (QrContact *)malloc((nold + nnew ? nold + nnew : 1u) *
                                 sizeof *merged);
    if (!merged) return -1;
    while (i < nold || j < nnew) {
        const QrContact *take = NULL;
        if (j >= nnew) take = &old_contact[i++];
        else if (i >= nold) take = &new_contact[j++];
        else {
            const QrContact *a = &old_contact[i], *b = &new_contact[j];
            if (a->cell_lo < b->cell_lo ||
                (a->cell_lo == b->cell_lo && a->cell_hi < b->cell_hi))
                take = &old_contact[i++];
            else if (b->cell_lo < a->cell_lo ||
                     (b->cell_lo == a->cell_lo && b->cell_hi < a->cell_hi))
                take = &new_contact[j++];
            else { take = &old_contact[i++]; j++; }
        }
        if (n && merged[n - 1].cell_lo == take->cell_lo &&
            merged[n - 1].cell_hi == take->cell_hi)
            continue;
        merged[n++] = *take;
    }
    *out_contact = merged;
    *out_count = n;
    return 0;
}

static uint8_t *qr_active_contact_band(const QrGrid *g,
                                       const QrContact *contact,
                                       size_t ncontact,
                                       const uint8_t *movable, int rings,
                                       size_t *out_active)
{
    uint8_t *active = (uint8_t *)calloc(g->nv, 1);
    if (!active) return NULL;
    for (size_t q = 0; q < ncontact; q++)
        for (int side = 0; side < 2; side++) {
            size_t cell = side ? contact[q].cell_hi : contact[q].cell_lo;
            for (int32_t i = g->cell_voff[cell]; i < g->cell_voff[cell + 1];
                 i++)
                active[(size_t)g->cell_verts[i]] = 1;
        }
    qr_dilate_active(g, active, movable, rings);
    if (out_active) {
        size_t n = 0;
        for (size_t v = 0; v < g->nv; v++) n += active[v] != 0;
        *out_active = n;
    }
    return active;
}

static int32_t *qr_active_list(const uint8_t *movable, const uint8_t *active,
                               size_t nv, size_t *out_n)
{
    size_t count = 0, n = 0;
    int32_t *list = NULL;
    for (size_t v = 0; v < nv; v++) count += (size_t)(movable[v] && active[v]);
    list = (int32_t *)malloc((count ? count : 1u) * sizeof *list);
    if (!list) return NULL;
    for (size_t v = 0; v < nv; v++)
        if (movable[v] && active[v]) list[n++] = (int32_t)v;
    *out_n = n;
    return list;
}

static int32_t *qr_band_face_list(const QrGrid *g, const uint8_t *movable,
                                  const uint8_t *active, size_t *out_count)
{
    size_t count = 0, n = 0;
    int32_t *list = NULL;
    for (size_t f = 0; f < g->nf; f++)
        for (int k = 0; k < 3; k++) {
            size_t v = (size_t)g->faces[f * 3 + (size_t)k];
            if (v < g->nv && movable[v] && active[v]) { count++; break; }
        }
    list = (int32_t *)malloc((count ? count : 1u) * sizeof *list);
    if (!list) return NULL;
    for (size_t f = 0; f < g->nf; f++)
        for (int k = 0; k < 3; k++) {
            size_t v = (size_t)g->faces[f * 3 + (size_t)k];
            if (v < g->nv && movable[v] && active[v]) {
                list[n++] = (int32_t)f;
                break;
            }
        }
    if (out_count) *out_count = n;
    return list;
}

static double qr_vertex_radius(const float *verts, size_t v, double axis_y,
                               double axis_x)
{
    return hypot((double)verts[v * 3 + 1] - axis_y,
                 (double)verts[v * 3 + 2] - axis_x);
}

static void qr_move_vertex_radial(float *verts, size_t v, double amount,
                                  double axis_y, double axis_x)
{
    double dy = (double)verts[v * 3 + 1] - axis_y;
    double dx = (double)verts[v * 3 + 2] - axis_x;
    double radius = hypot(dy, dx), scale;
    if (!(radius > 1e-9) || radius + amount <= 1e-6) return;
    scale = (radius + amount) / radius;
    verts[v * 3 + 1] = (float)(axis_y + scale * dy);
    verts[v * 3 + 2] = (float)(axis_x + scale * dx);
}

static size_t qr_project_radial_constraints(
    const QrGrid *g, float *displacement, const double *rest_radius,
    const QrContact *contact, size_t ncontact, const uint8_t *movable,
    const uint8_t *active, double pitch, double clearance, double relaxation,
    double maximum_impulse, double maximum_displacement, int reverse,
    size_t *hard_contacts, double *maximum_correction)
{
    const double order = pitch >= 0.0 ? 1.0 : -1.0;
    size_t projected = 0, hard = 0;
    double mx = 0.0;
    for (size_t at = 0; at < ncontact; at++) {
        size_t q = reverse ? ncontact - 1u - at : at;
        /* All quad-vertex inequalities are equivalent to the single support
         * inequality
         *   min(order*r_hi) >= max(order*r_lo) + clearance.
         * Project its currently active extrema.  Re-selecting them every
         * sweep is the usual active-support contact solve and avoids testing
         * duplicate triangle pairs billions of times on a whole scroll. */
        size_t vl = SIZE_MAX, vh = SIZE_MAX;
        double support_lo = -DBL_MAX, support_hi = DBL_MAX;
        double C, wl, wh, impulse, dl, dh;
        for (int32_t i = g->cell_voff[contact[q].cell_lo];
             i < g->cell_voff[contact[q].cell_lo + 1]; i++) {
            size_t v = (size_t)g->cell_verts[i];
            double sup = order * (rest_radius[v] + (double)displacement[v]);
            if (sup > support_lo) { support_lo = sup; vl = v; }
        }
        for (int32_t i = g->cell_voff[contact[q].cell_hi];
             i < g->cell_voff[contact[q].cell_hi + 1]; i++) {
            size_t v = (size_t)g->cell_verts[i];
            double sup = order * (rest_radius[v] + (double)displacement[v]);
            if (sup < support_hi) { support_hi = sup; vh = v; }
        }
        if (vl == SIZE_MAX || vh == SIZE_MAX) continue;
        C = support_hi - support_lo - clearance;
        if (C >= 0.0) continue;
        wl = (movable[vl] && active[vl]) ? 1.0 : 0.0;
        wh = (movable[vh] && active[vh]) ? 1.0 : 0.0;
        if (wl + wh <= 0.0) { hard++; continue; }
        impulse = relaxation * (-C) / (wl + wh + 1e-6);
        if (impulse > maximum_impulse) impulse = maximum_impulse;
        dl = -order * wl * impulse;
        dh = order * wh * impulse;
        if (wl) {
            double value = (double)displacement[vl] + dl;
            if (value > maximum_displacement) value = maximum_displacement;
            if (value < -maximum_displacement) value = -maximum_displacement;
            displacement[vl] = (float)value;
            if (fabs(dl) > mx) mx = fabs(dl);
        }
        if (wh) {
            double value = (double)displacement[vh] + dh;
            if (value > maximum_displacement) value = maximum_displacement;
            if (value < -maximum_displacement) value = -maximum_displacement;
            displacement[vh] = (float)value;
            if (fabs(dh) > mx) mx = fabs(dh);
        }
        projected++;
    }
    if (hard_contacts) *hard_contacts += hard;
    if (maximum_correction && mx > *maximum_correction)
        *maximum_correction = mx;
    return projected;
}

/* Jacobi over the vertex graph replaces the lattice red/black GS; called
 * twice per original sweep to compensate its slower mixing. */
static void qr_smooth_radial_displacement(const QrGrid *g, float *displacement,
                                          float *scratch, const int32_t *list,
                                          size_t nlist, double screen,
                                          double relaxation)
{
    for (int rep = 0; rep < 2; rep++) {
        ptrdiff_t i;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (i = 0; i < (ptrdiff_t)nlist; i++) {
            size_t v = (size_t)list[i];
            double sum = 0.0, degree = 0.0, target;
            for (int32_t e = g->adj_off[v]; e < g->adj_off[v + 1]; e++) {
                sum += (double)displacement[(size_t)g->adj[e]];
                degree += 1.0;
            }
            target = sum / (degree + screen);
            scratch[i] = (float)((1.0 - relaxation) * (double)displacement[v] +
                                 relaxation * target);
        }
        for (size_t k = 0; k < nlist; k++)
            displacement[(size_t)list[k]] = scratch[k];
    }
}

static void qr_apply_displacement(float *out, const float *rest, size_t nv,
                                  const float *displacement, double axis_y,
                                  double axis_x)
{
    memcpy(out, rest, nv * 3 * sizeof *out);
    for (size_t v = 0; v < nv; v++)
        if (displacement[v] != 0.0f)
            qr_move_vertex_radial(out, v, (double)displacement[v], axis_y,
                                  axis_x);
}

static void qr_constraint_violation(const QrGrid *g, const float *verts,
                                    const QrContact *contact, size_t ncontact,
                                    double pitch, double clearance,
                                    double axis_y, double axis_x,
                                    size_t *out_count, double *out_rms,
                                    double *out_maximum)
{
    const double order = pitch >= 0.0 ? 1.0 : -1.0;
    size_t count = 0;
    double ss = 0.0, mx = 0.0;
    for (size_t q = 0; q < ncontact; q++) {
        double support_lo = -DBL_MAX, support_hi = DBL_MAX, deficit;
        for (int32_t i = g->cell_voff[contact[q].cell_lo];
             i < g->cell_voff[contact[q].cell_lo + 1]; i++) {
            double sup = order * qr_vertex_radius(
                verts, (size_t)g->cell_verts[i], axis_y, axis_x);
            if (sup > support_lo) support_lo = sup;
        }
        for (int32_t i = g->cell_voff[contact[q].cell_hi];
             i < g->cell_voff[contact[q].cell_hi + 1]; i++) {
            double sup = order * qr_vertex_radius(
                verts, (size_t)g->cell_verts[i], axis_y, axis_x);
            if (sup < support_hi) support_hi = sup;
        }
        deficit = clearance - (support_hi - support_lo);
        if (!(deficit > 0.0)) continue;
        ss += deficit * deficit;
        if (deficit > mx) mx = deficit;
        count++;
    }
    if (out_count) *out_count = count;
    if (out_rms) *out_rms = count ? sqrt(ss / (double)count) : 0.0;
    if (out_maximum) *out_maximum = mx;
}

static void qr_displacement_stats(const float *p, const float *anchor,
                                  size_t nv, double *rms, double *maximum)
{
    double ss = 0.0, mx = 0.0;
    for (size_t v = 0; v < nv; v++) {
        double d2 = 0.0, d;
        for (int k = 0; k < 3; k++) {
            double dd = (double)p[v * 3 + (size_t)k] -
                        anchor[v * 3 + (size_t)k];
            d2 += dd * dd;
        }
        d = sqrt(d2);
        ss += d2;
        if (d > mx) mx = d;
    }
    *rms = nv ? sqrt(ss / (double)nv) : 0.0;
    *maximum = mx;
}

static void qr_radial_field_stats(const float *displacement, size_t nv,
                                  double *rms, double *maximum)
{
    double ss = 0.0, mx = 0.0;
    for (size_t v = 0; v < nv; v++) {
        double d = fabs((double)displacement[v]);
        ss += d * d;
        if (d > mx) mx = d;
    }
    if (rms) *rms = nv ? sqrt(ss / (double)nv) : 0.0;
    if (maximum) *maximum = mx;
}

/* Front-stall escalation: when accepted rounds stop shrinking the long-range
 * set, the intersection curve is travelling slower than the speculative
 * constraints extend.  Widen the collar/inflation/band so constraints stay
 * ahead of the front; drop back to base after a decisive improvement. */
static int qr_stall_escalation_step(size_t prev_long, size_t new_long,
                                    double collar_base, double collar_max,
                                    int *stall_rounds, double *collar,
                                    int *inflate, int *rings)
{
    int changed = 0;
    double improvement = prev_long
                       ? 1.0 - (double)new_long / (double)prev_long : 1.0;
    if (improvement >= 0.20) {
        if (*collar != collar_base || *inflate != 3 || *rings != 24)
            changed = 1;
        *collar = collar_base;
        *inflate = 3;
        *rings = 24;
        *stall_rounds = 0;
    } else if (improvement < 0.02) {
        (*stall_rounds)++;
        if (*stall_rounds >= 3) {
            *stall_rounds = 0;
            if (*collar < collar_max) {
                *collar = fmin(*collar + 4.0, collar_max);
                changed = 1;
            }
            if (*inflate < 5) { *inflate = 5; changed = 1; }
            if (*rings < 48) { *rings = 48; changed = 1; }
        }
    } else {
        *stall_rounds = 0;
    }
    return changed;
}

/* Settle toward the measured coil.  Untangling deliberately over-opens the
 * sheet (crop evidence: ~10 vox radial RMS from the pre-escape ribbon while
 * the ledger clearance only requires ~1.25), and the downstream CT snap then
 * starts far off the data -- the 6.49%->17.10% dark regression.  Once the
 * shell is conflict-free, repeatedly pull the radial displacement toward
 * zero (beta per round), re-project the persistent contact ledger to keep
 * the inter-turn clearance, and accept only candidates that a full exact
 * audit certifies as still conflict-free with strictly decreasing radial
 * RMS.  A rejected pull weakens beta ((1+beta)/2) up to three times before
 * the settle stops; the published state is always the last accepted
 * conflict-free geometry. */
static void qr_settle_toward_rest(
    const QrGrid *g, float *verts, const uint8_t *movable,
    const float *rest_target, const double *rest_radius,
    QrContact **contact_ledger, size_t *ncontact_ledger, double pitch,
    double clearance, double axis_y, double axis_x,
    double minimum_u_separation, size_t transaction_input_fold,
    int settle_rounds, double settle_beta, double displacement_bound,
    float *base, float *response, float *trial, float *displacement,
    float *scratch, uint8_t *orientation_input_bad, QrScan *cached_scan,
    IntersectionCleanupStats *cached_stats, QuadribbonUntangleStats *s)
{
    size_t nv = g->nv, nf = g->nf;
    const int32_t *faces = g->faces;
    const float *uv = g->uv;
    double beta = settle_beta;
    int betas_weakened = 0, accepted_rounds = 0, stale_merges = 0;
    size_t previous_reject_conflicts = SIZE_MAX;
    double previous_rms = -1.0;
    if (settle_rounds < 1 || !contact_ledger || !ncontact_ledger ||
        *ncontact_ledger == 0 || !(settle_beta > 0.0 && settle_beta < 1.0))
        return;
    for (int round = 0; round < settle_rounds; round++) {
        const QrContact *contact = *contact_ledger;
        size_t ncontact = *ncontact_ledger;
        size_t nactive = 0, nlist = 0, nband_face = 0;
        uint8_t *active = NULL;
        int32_t *list = NULL, *band_face = NULL;
        double rms_before = 0.0, max_before = 0.0;
        memcpy(base, verts, nv * 3 * sizeof *base);
        active = qr_active_contact_band(g, contact, ncontact, movable, 24,
                                        &nactive);
        if (!active) return;
        for (size_t v = 0; v < nv; v++) {
            double current_r = qr_vertex_radius(base, v, axis_y, axis_x);
            displacement[v] = (float)(current_r - rest_radius[v]);
            if (movable[v] && fabs((double)displacement[v]) > 1e-5)
                active[v] = 1;
        }
        qr_radial_field_stats(displacement, nv, &rms_before, &max_before);
        if (round == 0) {
            s->settle_rms_before = rms_before;
            s->settle_rms_after = rms_before;
            previous_rms = rms_before;
        }
        if (!(rms_before > 1e-6)) { free(active); break; }
        list = qr_active_list(movable, active, nv, &nlist);
        band_face = qr_band_face_list(g, movable, active, &nband_face);
        if (!list || !band_face) {
            free(list); free(band_face); free(active);
            return;
        }
        (void)qr_orientation_baseline(base, nv, faces, nf,
                                      orientation_input_bad);
        {
            ptrdiff_t bi;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (bi = 0; bi < (ptrdiff_t)nlist; bi++) {
                size_t v = (size_t)list[bi];
                displacement[v] = (float)((double)displacement[v] * beta);
            }
        }
        {
            size_t hard = 0, projected = 0;
            double max_impulse = 0.0;
            for (int sweep = 0; sweep < 96; sweep++) {
                qr_smooth_radial_displacement(g, displacement, scratch, list,
                                              nlist, 0.05, 0.70);
                for (int cs = 0; cs < 2; cs++)
                    projected += qr_project_radial_constraints(
                        g, displacement, rest_radius, contact, ncontact,
                        movable, active, pitch, clearance, 0.65, 0.35,
                        displacement_bound, (sweep + cs) & 1, &hard,
                        &max_impulse);
            }
            for (int cs = 0; cs < 32; cs++)
                projected += qr_project_radial_constraints(
                    g, displacement, rest_radius, contact, ncontact, movable,
                    active, pitch, clearance, 0.90, 0.5, displacement_bound,
                    cs & 1, &hard, &max_impulse);
            (void)projected;
        }
        qr_apply_displacement(response, rest_target, nv, displacement, axis_y,
                              axis_x);
        memcpy(trial, base, nv * 3 * sizeof *trial);
        {
            ptrdiff_t bi;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (bi = 0; bi < (ptrdiff_t)nlist; bi++) {
                size_t v = (size_t)list[bi];
                for (int k = 0; k < 3; k++)
                    trial[v * 3 + (size_t)k] = response[v * 3 + (size_t)k];
            }
        }
        {
            double rms_after = 0.0, max_after = 0.0;
            size_t orientation_bad;
            QrScan after_scan;
            IntersectionCleanupStats after;
            memset(&after_scan, 0, sizeof after_scan);
            memset(&after, 0, sizeof after);
            qr_radial_field_stats(displacement, nv, &rms_after, &max_after);
            orientation_bad = qr_orientation_violations_list(
                trial, base, nv, faces, band_face, nband_face,
                orientation_input_bad, NULL, NULL);
            /* GATE DISARMED (2026-08-30): every settle pull is audited and
             * committed; the gated verdict (the original predicate,
             * verbatim) is computed for the log only.  NOTE: the return-path
             * contact adoption in the #if 0 region below was reachable only
             * through rejection -- with gates disarmed it is dead code.  If
             * inspection shows settle re-colliding, that loss is the first
             * suspect. */
            if (qr_audit(trial, nv, faces, nf, uv, minimum_u_separation,
                         NULL, NULL, &after_scan, &after) != 0) {
                free(list); free(band_face); free(active);
                return;
            }
            {
                int gated = orientation_bad == 0 &&
                            rms_after + 1e-9 < 0.995 * previous_rms &&
                            after_scan.count == 0 && after.stab_pairs == 0 &&
                            after.fold_pairs <= transaction_input_fold;
                fprintf(stderr,
                        "[elastic-shell] settle round %d: beta %.4f, "
                        "displacement rms %.4f -> %.4f, exact "
                        "long/stab/fold=%zu/%zu/%zu, orientation=%zu (%s)\n",
                        round + 1, beta, rms_before, rms_after,
                        after_scan.count, after.stab_pairs, after.fold_pairs,
                        orientation_bad,
                        gated ? "accept" : "accept (gate-would-reject)");
            }
            s->settle_rounds_run++;
            (void)betas_weakened; (void)previous_reject_conflicts;
            (void)stale_merges;
            memcpy(verts, trial, nv * 3 * sizeof *verts);
            free(cached_scan->pair);
            *cached_scan = after_scan;
            *cached_stats = after;
            s->settle_rms_after = rms_after;
            accepted_rounds++;
            previous_rms = rms_after;
            previous_reject_conflicts = SIZE_MAX;
            stale_merges = 0;
            free(list); free(band_face); free(active);
            continue;
#if 0       /* disarmed: settle rejection handling (return-path contact
             * adoption + beta weakening), unreachable behind the
             * unconditional commit above */
            if (audited && (after_scan.count > 0 || after.stab_pairs > 0)) {
                /* The return path crossed pairs the escape ledger never saw --
                 * the pull re-collides where the coil originally
                 * interpenetrated.  Adopt the failed trial's exact conflicts
                 * as contacts and retry: the projection then holds clearance
                 * exactly there, so the next pull relaxes only through
                 * genuinely free space.  Only when merging stops shrinking
                 * the rejection do we weaken the pull. */
                int made_progress = after.conflicts < previous_reject_conflicts;
                size_t nnew = 0;
                QrContact *new_contact;
                previous_reject_conflicts = after.conflicts;
                new_contact = qr_build_contacts(g, trial, &after_scan, 4.0, 3,
                                                &nnew);
                free(after_scan.pair);
                if (new_contact) {
                    QrContact *merged = NULL;
                    size_t nmerged = 0;
                    if (qr_merge_ledger(*contact_ledger, *ncontact_ledger,
                                        new_contact, nnew, &merged,
                                        &nmerged) == 0) {
                        free(*contact_ledger);
                        *contact_ledger = merged;
                        *ncontact_ledger = nmerged;
                        fprintf(stderr,
                                "[elastic-shell] settle adopted %zu "
                                "return-path contacts (ledger %zu)\n",
                                nnew, nmerged);
                    }
                    free(new_contact);
                }
                free(list); free(band_face); free(active);
                if (made_progress) { stale_merges = 0; continue; }
                if (++stale_merges < 2) continue;
                stale_merges = 0;
                betas_weakened++;
                if (betas_weakened > 3) break;
                beta = (1.0 + beta) * 0.5;
                continue;
            }
            if (audited) free(after_scan.pair);
            free(list); free(band_face); free(active);
            /* A rejected pull is either blocked by the ledger clearance floor
             * or pulled too hard through it; weaken toward the identity and
             * retry, then stop -- the accepted state already satisfies every
             * gate. */
            betas_weakened++;
            if (betas_weakened > 3) break;
            beta = (1.0 + beta) * 0.5;
#endif
        }
    }
    s->settle_accepted = accepted_rounds;
    fprintf(stderr,
            "[elastic-shell] settle result: rounds=%d accepted=%d, "
            "displacement rms %.4f->%.4f\n",
            s->settle_rounds_run, accepted_rounds, s->settle_rms_before,
            s->settle_rms_after);
}

static int qr_resolve_self_collisions(
    const QrGrid *g, float *verts, const uint8_t *movable,
    const float *rest_target, double axis_y, double axis_x,
    double minimum_u_separation, const QuadribbonUntangleOpts *o,
    QuadribbonUntangleStats *s)
{
    size_t nv = g->nv, nf = g->nf;
    const int32_t *faces = g->faces;
    const float *uv = g->uv;
    float *base = (float *)malloc(nv * 3 * sizeof *base);
    float *response = (float *)malloc(nv * 3 * sizeof *response);
    float *trial = (float *)malloc(nv * 3 * sizeof *trial);
    float *transaction_start = (float *)malloc(nv * 3 *
                                               sizeof *transaction_start);
    float *displacement = (float *)calloc(nv, sizeof *displacement);
    float *scratch = (float *)malloc(nv * sizeof *scratch);
    double *rest_radius = (double *)malloc(nv * sizeof *rest_radius);
    float *keep = (float *)malloc(nv * sizeof *keep);
    float *keep_tmp = (float *)malloc(nv * sizeof *keep_tmp);
    size_t *local_before = (size_t *)calloc(nf, sizeof *local_before);
    size_t *local_after = (size_t *)calloc(nf, sizeof *local_after);
    uint8_t *orientation_input_bad = (uint8_t *)malloc(nf);
    double pitch = 0.0, clearance = 1.25;
    size_t transaction_input_fold = 0;
    size_t best_long = SIZE_MAX, best_violation_count = SIZE_MAX;
    IntersectionCleanupStats best_stats;
    int best_initialized = 0, rounds_since_best = 0;
    const double collar_base = o->collision_collar > 0.0
                             ? o->collision_collar : 4.0;
    const double collar_max = collar_base + 8.0;
    double collar = collar_base;
    int contact_inflate = 3, band_rings = 24, stall_rounds = 0;
    size_t prev_accepted_long = 0;
    int have_prev_accepted_long = 0;
    QrContact *contact_ledger = NULL;
    size_t ncontact_ledger = 0;
    QrScan cached_scan;
    IntersectionCleanupStats cached_stats;
    int cached_audit = 0;
    int rc = -1;
    memset(&best_stats, 0, sizeof best_stats);
    memset(&cached_scan, 0, sizeof cached_scan);
    memset(&cached_stats, 0, sizeof cached_stats);
    if (!base || !response || !trial || !transaction_start || !displacement ||
        !scratch || !rest_radius || !keep || !keep_tmp || !local_before ||
        !local_after || !orientation_input_bad)
        goto fail;
    qr_audit_reset();
    /* Strict mode retains the entry state here for atomic rollback.  Default
     * inspectable-output mode updates the same buffer only when a fully
     * audited candidate improves the exact collision objective. */
    memcpy(transaction_start, verts, nv * 3 * sizeof *transaction_start);
    for (size_t v = 0; v < nv; v++)
        rest_radius[v] = qr_vertex_radius(rest_target, v, axis_y, axis_x);
    for (int round = 0; round < o->collision_rounds; round++) {
        QrScan before_scan;
        IntersectionCleanupStats before;
        size_t nnew_contact = 0, nmerged_contact = 0, nactive = 0;
        QrContact *new_contact = NULL, *merged_contact = NULL;
        const QrContact *contact = NULL;
        size_t ncontact = 0, nlist = 0, nband_face = 0;
        uint8_t *active = NULL;
        int32_t *list = NULL, *band_face = NULL;
        int accepted = 0;
        memset(&before, 0, sizeof before);
        if (cached_audit) {
            before_scan = cached_scan;
            before = cached_stats;
            memset(&cached_scan, 0, sizeof cached_scan);
            memset(&cached_stats, 0, sizeof cached_stats);
            cached_audit = 0;
            fprintf(stderr,
                    "[cloth-collision] reusing accepted exact audit for next "
                    "active set\n");
        } else if (qr_audit(verts, nv, faces, nf, uv, minimum_u_separation,
                            NULL, local_before, &before_scan, &before) != 0) {
            goto fail;
        }
        if (round == 0) {
            s->input_conflicts = before.conflicts;
            s->input_long_conflicts = before_scan.count;
            transaction_input_fold = before.fold_pairs;
            best_long = before_scan.count;
            best_stats = before;
            best_initialized = 1;
            if (before_scan.count &&
                qr_infer_pitch(g, verts, movable, axis_y, axis_x, &before_scan,
                               o->wrap_pitch_hint, &pitch) != 0) {
                fprintf(stderr,
                        "[cloth-collision] refusing ambiguous phase ordering\n");
                free(before_scan.pair);
                break;
            }
            s->pitch = pitch;
            /* Numerical mid-surface clearance, deliberately much smaller than
             * the measured inter-turn pitch.  The rest energy, not clearance,
             * determines where the metal settles after it opens. */
            clearance = 1.25;
            s->clearance = clearance;
        }
        fprintf(stderr,
                "[cloth-collision] active-set round %d: exact=%zu, long=%zu "
                "[overlap=%zu stab=%zu fold=%zu], clearance=%.3f\n",
                round + 1, before.conflicts, before_scan.count,
                before.overlap_pairs, before.stab_pairs, before.fold_pairs,
                clearance);
        if (before_scan.count == 0) {
            cached_scan = before_scan;
            cached_stats = before;
            cached_audit = 1;
            break;
        }
        new_contact = qr_build_contacts(g, verts, &before_scan, collar,
                                        contact_inflate, &nnew_contact);
        if (!new_contact) { free(before_scan.pair); goto fail; }
        if (qr_merge_ledger(contact_ledger, ncontact_ledger, new_contact,
                            nnew_contact, &merged_contact,
                            &nmerged_contact) != 0) {
            free(new_contact);
            free(before_scan.pair);
            goto fail;
        }
        free(new_contact);
        free(contact_ledger);
        contact_ledger = merged_contact;
        ncontact_ledger = nmerged_contact;
        /* Ledger pruning (the v144 plateau lever): the merge above only ever
         * GROWS the ledger, so support inequalities from pairs that have long
         * since separated keep constraining every later solve.  Drop entries
         * whose cell AABBs now sit beyond a retention collar; a pruned pair
         * that re-approaches is re-adopted from that round's scan by
         * qr_build_contacts, and the exact audit still gates every
         * acceptance. */
        {
            size_t kept = 0;
            double retain = collar * 3.0;
            for (size_t q = 0; q < ncontact_ledger; q++) {
                if (qr_cells_within_collar(g, verts,
                                           contact_ledger[q].cell_lo,
                                           contact_ledger[q].cell_hi,
                                           retain))
                    contact_ledger[kept++] = contact_ledger[q];
            }
            if (kept < ncontact_ledger)
                fprintf(stderr,
                        "[cloth-collision] ledger prune: %zu -> %zu "
                        "(retain %.1f vox)\n",
                        ncontact_ledger, kept, retain);
            ncontact_ledger = kept;
        }
        contact = contact_ledger;
        ncontact = ncontact_ledger;
        active = qr_active_contact_band(g, contact, ncontact, movable,
                                        band_rings, &nactive);
        if (!active) { free(before_scan.pair); goto fail; }
        memcpy(base, verts, nv * 3 * sizeof *base);
        {
            size_t input_bad_faces = qr_orientation_baseline(
                base, nv, faces, nf, orientation_input_bad);
            if (input_bad_faces)
                fprintf(stderr,
                        "[cloth-collision] input orientation-bad faces=%zu "
                        "excluded from candidate preflight\n",
                        input_bad_faces);
        }
        nactive = 0;
        for (size_t v = 0; v < nv; v++) {
            double current_r = qr_vertex_radius(base, v, axis_y, axis_x);
            displacement[v] = (float)(current_r - rest_radius[v]);
            /* The coarse escape is only an initializer.  Every vertex it
             * moved belongs to the elastic return solve even if it is no
             * longer on today's exact contact curve. */
            if (movable[v] && fabs(current_r - rest_radius[v]) > 1e-5)
                active[v] = 1;
            nactive += active[v] != 0;
        }
        list = qr_active_list(movable, active, nv, &nlist);
        band_face = qr_band_face_list(g, movable, active, &nband_face);
        if (!list || !band_face) {
            free(list); free(band_face);
            free(active);
            free(before_scan.pair);
            goto fail;
        }
        {
            size_t hard = 0, projected = 0;
            double max_impulse = 0.0;
            const double displacement_bound = o->displacement_bound > 0.0
                                            ? o->displacement_bound : 128.0;
            for (int sweep = 0; sweep < 512; sweep++) {
                qr_smooth_radial_displacement(g, displacement, scratch, list,
                                              nlist, 0.006, 0.70);
                for (int cs = 0; cs < 2; cs++)
                    projected += qr_project_radial_constraints(
                        g, displacement, rest_radius, contact, ncontact,
                        movable, active, pitch, clearance, 0.65, 0.35,
                        displacement_bound, (sweep + cs) & 1, &hard,
                        &max_impulse);
            }
            /* Finish on the inequalities, not on a smoothing sweep. */
            for (int cs = 0; cs < 128; cs++)
                projected += qr_project_radial_constraints(
                    g, displacement, rest_radius, contact, ncontact, movable,
                    active, pitch, clearance, 0.90, 0.5, displacement_bound,
                    cs & 1, &hard, &max_impulse);
            qr_apply_displacement(response, rest_target, nv, displacement,
                                  axis_y, axis_x);
            s->contacts_built += ncontact;
            s->hard_contacts += hard;
            {
                double full_rms = 0.0, full_max = 0.0;
                size_t violation_count_before = 0;
                double violation_rms_before = 0.0, violation_max_before = 0.0;
                qr_displacement_stats(response, base, nv, &full_rms,
                                      &full_max);
                qr_constraint_violation(g, base, contact, ncontact, pitch,
                                        clearance, axis_y, axis_x,
                                        &violation_count_before,
                                        &violation_rms_before,
                                        &violation_max_before);
                fprintf(stderr,
                        "[cloth-collision] contacts=%zu, active=%zu/%zu, "
                        "projections=%zu, hard=%zu, response rms/max="
                        "%.6f/%.6f; constraint violations=%zu rms/max "
                        "%.6f/%.6f\n",
                        ncontact, nactive, nv, projected, hard, full_rms,
                        full_max, violation_count_before,
                        violation_rms_before, violation_max_before);
                /* alpha ladder with feathered local rollbacks */
                {
                    static const double alpha[7] =
                        { 1.0, 0.75, 0.5, 0.375, 0.25, 0.125, 0.0625 };
                    /* Whole-sheet evidence shows that a large, useful response
                     * can be down to only a handful of unsafe front vertices
                     * at rollback six.  Cutting it off there discards the
                     * response and replaces one local repair with several more
                     * global contact rounds. */
                    const int maximum_local_rollbacks = 12;
                    memcpy(trial, base, nv * 3 * sizeof *trial);
                    for (int attempt = 0; attempt < 7; attempt++) {
                        for (size_t v = 0; v < nv; v++) keep[v] = 1.0f;
                        for (int rollback = 0;
                             rollback <= maximum_local_rollbacks; rollback++) {
                            size_t orientation_bad;
                            QrScan after_scan;
                            IntersectionCleanupStats after;
                            size_t violation_count_after = 0;
                            double violation_rms_after = 0.0;
                            double violation_max_after = 0.0;
                            size_t churn_cap;
                            int phase_progress, bounded_churn, safe;
                            {
                                ptrdiff_t bi;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                                for (bi = 0; bi < (ptrdiff_t)nlist; bi++) {
                                    size_t v = (size_t)list[bi];
                                    for (int k = 0; k < 3; k++)
                                        trial[v * 3 + (size_t)k] =
                                            (float)((double)base[v * 3 + (size_t)k] +
                                                alpha[attempt] *
                                                (double)keep[v] *
                                                ((double)response[v * 3 + (size_t)k] -
                                                 base[v * 3 + (size_t)k]));
                                }
                            }
                            orientation_bad = qr_orientation_violations_list(
                                trial, base, nv, faces, band_face, nband_face,
                                orientation_input_bad, NULL, NULL);
                            /* GATE DISARMED (2026-08-30): preflight is
                             * measurement-only; see the turn-order note. */
                            if (orientation_bad)
                                fprintf(stderr,
                                        "[cloth-collision] preflight "
                                        "measurement alpha %.4f: orientation "
                                        "violations=%zu (gate disarmed, "
                                        "proceeding to audit)\n",
                                        alpha[attempt], orientation_bad);
#if 0                       /* disarmed: preflight rejection + rollback */
                            if (orientation_bad) {
                                size_t orientation_marked = 0;
                                fprintf(stderr,
                                        "[cloth-collision] reject alpha %.4f "
                                        "rollback %d before exact audit: "
                                        "orientation violations=%zu\n",
                                        alpha[attempt], rollback,
                                        orientation_bad);
                                if (rollback == maximum_local_rollbacks)
                                    break;
                                (void)qr_orientation_violations_list(
                                    trial, base, nv, faces, band_face,
                                    nband_face, orientation_input_bad, keep,
                                    &orientation_marked);
                                if (orientation_marked == 0) break;
                                qr_feather_motion_scale(g, keep, keep_tmp, 8);
                                continue;
                            }
#endif
                            memset(&after, 0, sizeof after);
                            if (qr_audit(trial, nv, faces, nf, uv,
                                         minimum_u_separation, NULL,
                                         local_after, &after_scan,
                                         &after) != 0) {
                                free(list); free(band_face);
                                free(active);
                                free(before_scan.pair);
                                goto fail;
                            }
                            qr_constraint_violation(
                                g, trial, contact, ncontact, pitch, clearance,
                                axis_y, axis_x, &violation_count_after,
                                &violation_rms_after, &violation_max_after);
                            /* Existing penetration has no collision-free
                             * continuous endpoint step: its intersection
                             * curve must sweep across neighbouring
                             * tessellation faces before leaving the sheet.
                             * Permit that pair identity/count churn inside
                             * this atomic transaction; phase violation and
                             * orientation, rather than pair monotonicity, are
                             * the progress/safety measures. */
                            churn_cap = before.conflicts <= SIZE_MAX / 8u
                                      ? before.conflicts * 8u
                                      : SIZE_MAX - before.conflicts;
                            if (churn_cap < 4096u) churn_cap = 4096u;
                            phase_progress =
                                violation_count_after < violation_count_before ||
                                violation_rms_after + 1e-6 <
                                    0.995 * violation_rms_before ||
                                after_scan.count < before_scan.count;
                            bounded_churn = after.conflicts <=
                                            before.conflicts + churn_cap;
                            safe = phase_progress && bounded_churn &&
                                   after.fold_pairs <= before.fold_pairs &&
                                   orientation_bad == 0;
                            fprintf(stderr,
                                    "[cloth-collision] %s alpha %.4f rollback "
                                    "%d: exact %zu->%zu, long %zu->%zu, "
                                    "folds %zu->%zu; phase violations "
                                    "%zu->%zu rms %.6f->%.6f\n",
                                    safe ? "accept"
                                         : "accept (gate-would-reject)",
                                    alpha[attempt],
                                    rollback, before.conflicts,
                                    after.conflicts, before_scan.count,
                                    after_scan.count, before.fold_pairs,
                                    after.fold_pairs, violation_count_before,
                                    violation_count_after,
                                    violation_rms_before,
                                    violation_rms_after);
                            /* GATE DISARMED: commit every audited trial. */
                            {
                                memcpy(verts, trial, nv * 3 * sizeof *verts);
                                accepted = 1;
                                s->accepted_rounds++;
                                {
                                    double rms = 0.0, mx = 0.0;
                                    qr_displacement_stats(verts, base, nv,
                                                          &rms, &mx);
                                    s->movement_rms = rms;
                                    s->movement_max = mx;
                                }
                                if (!o->require_collision_free) {
                                    int lex_better =
                                        after_scan.count < best_long ||
                                        (after_scan.count == best_long &&
                                         (after.stab_pairs <
                                              best_stats.stab_pairs ||
                                          (after.stab_pairs ==
                                               best_stats.stab_pairs &&
                                           after.conflicts <
                                               best_stats.conflicts)));
                                    /* The escape is deliberately non-monotone
                                     * in pair counts: while the intersection
                                     * front travels, phase violations are the
                                     * trustworthy progress signal (full-sheet
                                     * v139 collapsed 661k->89k violations in
                                     * one round yet was stopped by
                                     * lexicographic-only patience).  Either
                                     * kind of new best resets patience; only
                                     * the lexicographic best is published. */
                                    int phase_better = violation_count_after <
                                                       best_violation_count;
                                    if (lex_better) {
                                        memcpy(transaction_start, trial,
                                               nv * 3 * sizeof *trial);
                                        best_long = after_scan.count;
                                        best_stats = after;
                                    }
                                    if (lex_better || phase_better)
                                        rounds_since_best = 0;
                                    else
                                        rounds_since_best++;
                                    if (phase_better)
                                        best_violation_count =
                                            violation_count_after;
                                    fprintf(stderr,
                                            "[cloth-collision] patience "
                                            "%d/%d; best long/stab/total="
                                            "%zu/%zu/%zu, best "
                                            "phase-violations=%zu\n",
                                            rounds_since_best,
                                            o->collision_patience, best_long,
                                            best_stats.stab_pairs,
                                            best_stats.conflicts,
                                            best_violation_count);
                                }
                                if (have_prev_accepted_long &&
                                    qr_stall_escalation_step(
                                        prev_accepted_long, after_scan.count,
                                        collar_base, collar_max,
                                        &stall_rounds, &collar,
                                        &contact_inflate, &band_rings))
                                    fprintf(stderr,
                                            "[cloth-collision] front stall "
                                            "response: collar=%.1f inflate=%d "
                                            "rings=%d\n",
                                            collar, contact_inflate,
                                            band_rings);
                                prev_accepted_long = after_scan.count;
                                have_prev_accepted_long = 1;
                                /* Carry the accepted audit into the next
                                 * active set instead of rebuilding the same
                                 * BVH and rediscovering the same contacts. */
                                cached_scan = after_scan;
                                cached_stats = after;
                                cached_audit = 1;
                                {
                                    size_t *tmp = local_before;
                                    local_before = local_after;
                                    local_after = tmp;
                                }
                                break;
                            }
#if 0                       /* disarmed: rejection + local rollback ladder */
                            if (rollback == maximum_local_rollbacks ||
                                !phase_progress || !bounded_churn) {
                                free(after_scan.pair);
                                break;
                            }
                            /* A full-sheet contact move can safely untangle
                             * one winding while folding a few cells at the
                             * moving intersection front.  Freeze only
                             * vertices on newly-created local conflicts or
                             * orientation failures, feather that freeze into
                             * the response, and re-audit.  Long-range pair
                             * churn is deliberately not marked: it is the
                             * intersection curve travelling to the edge. */
                            {
                                size_t marked = 0, orientation_marked = 0;
                                for (size_t f = 0; f < nf; f++)
                                    if (local_after[f] > local_before[f])
                                        for (int k = 0; k < 3; k++) {
                                            size_t v = (size_t)
                                                faces[f * 3 + (size_t)k];
                                            if (keep[v] > 0.0f) {
                                                keep[v] = 0.0f;
                                                marked++;
                                            }
                                        }
                                (void)qr_orientation_violations_list(
                                    trial, base, nv, faces, band_face,
                                    nband_face, orientation_input_bad, keep,
                                    &orientation_marked);
                                marked += orientation_marked;
                                free(after_scan.pair);
                                if (marked == 0) break;
                                qr_feather_motion_scale(g, keep, keep_tmp, 8);
                                fprintf(stderr,
                                        "[cloth-collision] locally rolled "
                                        "back %zu fold/orientation vertices "
                                        "and feathered 8 rings\n",
                                        marked);
                            }
#endif
                        }
                        if (accepted) break;
                    }
                }
            }
        }
        free(list);
        free(band_face);
        free(active);
        if (accepted) free(before_scan.pair);
        if (accepted && !o->require_collision_free &&
            o->collision_patience > 0 &&
            rounds_since_best >= o->collision_patience) {
            fprintf(stderr,
                    "[cloth-collision] stopping after %d accepted rounds "
                    "without a better exact-audited topology (patience=%d); "
                    "the FINAL state ships (gates disarmed 2026-08-30)\n",
                    rounds_since_best, o->collision_patience);
            break;
        }
        if (!accepted) {
            fprintf(stderr,
                    "[cloth-collision] no exact-audit-safe response remained\n");
            cached_scan = before_scan;
            cached_stats = before;
            cached_audit = 1;
            break;
        }
    }
    /* GATE DISARMED (2026-08-30): settle runs whenever a contact ledger
     * exists.  The original entry additionally required a fully clean
     * audited state:
     *   cached_scan.count == 0 && cached_stats.stab_pairs == 0        */
    if (cached_audit && ncontact_ledger > 0 && pitch != 0.0 &&
        (cached_scan.count != 0 || cached_stats.stab_pairs != 0))
        fprintf(stderr,
                "[elastic-shell] settle entry with unresolved long=%zu "
                "stab=%zu (clean-entry gate disarmed)\n",
                cached_scan.count, cached_stats.stab_pairs);
    if (cached_audit &&
        ncontact_ledger > 0 && pitch != 0.0)
        qr_settle_toward_rest(g, verts, movable, rest_target, rest_radius,
                              &contact_ledger, &ncontact_ledger, pitch,
                              clearance, axis_y, axis_x, minimum_u_separation,
                              transaction_input_fold, o->settle_rounds,
                              o->settle_beta,
                              o->displacement_bound > 0.0
                                  ? o->displacement_bound : 128.0,
                              base, response, trial, displacement, scratch,
                              orientation_input_bad, &cached_scan,
                              &cached_stats, s);
    {
        QrScan final_scan;
        IntersectionCleanupStats final_stats;
        int complete;
        memset(&final_scan, 0, sizeof final_scan);
        memset(&final_stats, 0, sizeof final_stats);
        if (cached_audit) {
            final_scan = cached_scan;
            final_stats = cached_stats;
            memset(&cached_scan, 0, sizeof cached_scan);
            cached_audit = 0;
        } else if (qr_audit(verts, nv, faces, nf, uv, minimum_u_separation,
                            NULL, NULL, &final_scan, &final_stats) != 0) {
            goto fail;
        }
        complete = final_scan.count == 0 && final_stats.stab_pairs == 0 &&
                   final_stats.fold_pairs <= transaction_input_fold;
        /* GATE DISARMED (2026-08-30): the FINAL state ships; the
         * lexicographic-best snapshot is logged as measurement only. */
        if (!complete && best_initialized)
            fprintf(stderr,
                    "[cloth-collision] gate disarmed: retaining FINAL state "
                    "(long=%zu stab=%zu folds=%zu); lex-best snapshot would "
                    "have been long=%zu stab=%zu folds=%zu (input folds=%zu)\n",
                    final_scan.count, final_stats.stab_pairs,
                    final_stats.fold_pairs, best_long, best_stats.stab_pairs,
                    best_stats.fold_pairs, transaction_input_fold);
#if 0   /* disarmed: lexicographic-best partial restore */
        if (!complete && !o->require_collision_free && best_initialized) {
            free(final_scan.pair);
            memset(&final_scan, 0, sizeof final_scan);
            memcpy(verts, transaction_start, nv * 3 * sizeof *verts);
            final_stats = best_stats;
            final_scan.count = best_long;
            s->retained_partial = 1;
            fprintf(stderr,
                    "[cloth-collision] PARTIAL: retained best exact-audited "
                    "state; remaining long=%zu, stab=%zu, folds=%zu (input "
                    "folds=%zu)\n",
                    final_scan.count, final_stats.stab_pairs,
                    final_stats.fold_pairs, transaction_input_fold);
        }
#endif
        s->complete = complete;
        s->output_conflicts = final_stats.conflicts;
        s->output_long_conflicts = final_scan.count;
        fprintf(stderr,
                "[cloth-collision] result: accepted rounds=%d, exact "
                "%zu->%zu, long %zu->%zu [overlap=%zu stab=%zu fold=%zu]\n",
                s->accepted_rounds, s->input_conflicts, s->output_conflicts,
                s->input_long_conflicts, s->output_long_conflicts,
                final_stats.overlap_pairs, final_stats.stab_pairs,
                final_stats.fold_pairs);
        qr_audit_report("elastic-shell");
        /* GATE DISARMED (2026-08-30): the full-transaction rollback is
         * suppressed; an armed require_collision_free is logged only. */
        if (!complete && o->require_collision_free)
            fprintf(stderr,
                    "[cloth-collision] gate disarmed: require_collision_free "
                    "rollback suppressed; shipping the final state\n");
#if 0   /* disarmed: full-transaction rollback */
        if (!complete && o->require_collision_free) {
            fprintf(stderr,
                    "[cloth-collision] ROLLBACK: the elastic-shell "
                    "transaction did not reach a non-interpenetrating state; "
                    "restoring its input\n");
            memcpy(verts, transaction_start, nv * 3 * sizeof *verts);
            free(final_scan.pair);
            rc = 1;
            goto cleanup;
        }
#endif
        free(final_scan.pair);
    }
    rc = 0;
cleanup:
    free(cached_scan.pair);
    free(contact_ledger);
    free(base); free(response); free(trial); free(transaction_start);
    free(displacement); free(scratch); free(rest_radius);
    free(keep); free(keep_tmp);
    free(local_before); free(local_after);
    free(orientation_input_bad);
    return rc;
fail:
    rc = -1;
    goto cleanup;
}

/* ---- public API ----------------------------------------------------------- */

void QuadribbonUntangle_defaults(QuadribbonUntangleOpts *o)
{
    memset(o, 0, sizeof *o);
    o->require_collision_free = 0;
    o->collision_patience = 12;
    o->collision_rounds = 64;
    o->collision_collar = 4.0;
    o->displacement_bound = 128.0;
    o->settle_rounds = 16;
    o->settle_beta = 0.70;
    o->minimum_u_separation = 0.0;
    o->wrap_pitch_hint = 0.0;
}

/* Characterize the residual long-range conflicts (measurement only, touches
 * nothing): a pair at shell-thickness distance with near-parallel faces is a
 * prediction-shell DUPLICATE (front/back surface of one thick predicted
 * wrap -- resolvable by dedup labeling, never by displacement); a crossing
 * pair is a true tangle.  Prints the split so the plateau's composition is
 * a number, not a guess. */
static void qr_characterize_residual(const QrGrid *g, float *verts,
                                     double axis_y, double axis_x,
                                     double min_u_sep)
{
    size_t nf = g->nf, nv = g->nv;
    const int32_t *faces = g->faces;
    QrScan scan;
    IntersectionCleanupStats st;
    size_t *facedeg = (size_t *)calloc(nf, sizeof *facedeg);
    double *dist = NULL;
    size_t ndup = 0, npar = 0, ncross = 0;
    memset(&scan, 0, sizeof scan);
    memset(&st, 0, sizeof st);
    if (facedeg == NULL) return;
    if (qr_audit(verts, nv, faces, nf, g->uv, min_u_sep, NULL, facedeg,
                 &scan, &st) != 0) {
        free(facedeg);
        return;
    }
    if (scan.count == 0) {
        free(scan.pair);
        free(facedeg);
        return;
    }
    dist = (double *)malloc((scan.count + 1) * sizeof *dist);
    if (dist == NULL) {
        free(scan.pair);
        free(facedeg);
        return;
    }
    for (size_t q = 0; q < scan.count; q++) {
        const int32_t *fa = &faces[(size_t)scan.pair[q].face_lo * 3];
        const int32_t *fb = &faces[(size_t)scan.pair[q].face_hi * 3];
        double ca[3] = { 0.0, 0.0, 0.0 }, cb[3] = { 0.0, 0.0, 0.0 };
        double na[3] = { 0.0, 0.0, 0.0 }, nb[3] = { 0.0, 0.0, 0.0 };
        double d2 = 0.0, ndot = 0.0, la = 0.0, lb = 0.0;
        for (int k = 0; k < 3; k++)
            for (int j = 0; j < 3; j++) {
                ca[j] += (double)verts[(size_t)fa[k] * 3 + j] / 3.0;
                cb[j] += (double)verts[(size_t)fb[k] * 3 + j] / 3.0;
            }
        {
            double ea[3], eb[3], fa2[3], fb2[3];
            for (int j = 0; j < 3; j++) {
                ea[j] = (double)verts[(size_t)fa[1] * 3 + j] -
                        (double)verts[(size_t)fa[0] * 3 + j];
                fa2[j] = (double)verts[(size_t)fa[2] * 3 + j] -
                         (double)verts[(size_t)fa[0] * 3 + j];
                eb[j] = (double)verts[(size_t)fb[1] * 3 + j] -
                        (double)verts[(size_t)fb[0] * 3 + j];
                fb2[j] = (double)verts[(size_t)fb[2] * 3 + j] -
                         (double)verts[(size_t)fb[0] * 3 + j];
            }
            na[0] = ea[1] * fa2[2] - ea[2] * fa2[1];
            na[1] = ea[2] * fa2[0] - ea[0] * fa2[2];
            na[2] = ea[0] * fa2[1] - ea[1] * fa2[0];
            nb[0] = eb[1] * fb2[2] - eb[2] * fb2[1];
            nb[1] = eb[2] * fb2[0] - eb[0] * fb2[2];
            nb[2] = eb[0] * fb2[1] - eb[1] * fb2[0];
        }
        for (int j = 0; j < 3; j++) {
            double dd = ca[j] - cb[j];
            d2 += dd * dd;
            la += na[j] * na[j];
            lb += nb[j] * nb[j];
            ndot += na[j] * nb[j];
        }
        dist[q] = sqrt(d2);
        ndot = (la > 1e-20 && lb > 1e-20)
             ? fabs(ndot) / sqrt(la * lb) : 0.0;
        if (dist[q] <= 4.5 && ndot >= 0.80) ndup++;
        else if (ndot >= 0.80) npar++;
        else ncross++;
    }
    qsort(dist, scan.count, sizeof *dist, qr_cmp_dbl);
    fprintf(stderr,
            "[qr-untangle] residual characterization: %zu long pair(s) -> "
            "%zu shell-duplicate-class (d<=4.5, |n.n|>=0.8), %zu "
            "parallel-far, %zu crossing; dist p10/p50/p90 = "
            "%.2f/%.2f/%.2f vox\n",
            scan.count, ndup, npar, ncross,
            dist[scan.count / 10], dist[scan.count / 2],
            dist[(scan.count * 9) / 10]);
    {
        /* attribution evidence: where do the residuals LIVE?  |du| separates
         * duplicate placements from genuine cross-wrap junk; radius r pins
         * whether they concentrate in the fused umbilicus core. */
        double *duv = dist;   /* reuse the buffer */
        size_t nd2 = 0;
        for (size_t q = 0; q < scan.count; q++)
            duv[nd2++] = fabs(scan.pair[q].u_hi - scan.pair[q].u_lo);
        qsort(duv, nd2, sizeof *duv, qr_cmp_dbl);
        fprintf(stderr,
                "[qr-untangle] residual |du| p10/p50/p90 = "
                "%.0f/%.0f/%.0f vox",
                duv[nd2 / 10], duv[nd2 / 2], duv[(nd2 * 9) / 10]);
        nd2 = 0;
        for (size_t q = 0; q < scan.count; q++) {
            const int32_t *fa = &faces[(size_t)scan.pair[q].face_lo * 3];
            double cy = 0.0, cx2 = 0.0;
            for (int k = 0; k < 3; k++) {
                cy += (double)verts[(size_t)fa[k] * 3 + 1] / 3.0;
                cx2 += (double)verts[(size_t)fa[k] * 3 + 2] / 3.0;
            }
            duv[nd2++] = hypot(cy - axis_y, cx2 - axis_x);
        }
        qsort(duv, nd2, sizeof *duv, qr_cmp_dbl);
        fprintf(stderr,
                "; r p10/p50/p90 = %.0f/%.0f/%.0f vox\n",
                duv[nd2 / 10], duv[nd2 / 2], duv[(nd2 * 9) / 10]);
    }
    free(dist);
    free(scan.pair);
    free(facedeg);
}

int QuadribbonUntangle_run(float *verts, size_t nv, const int32_t *faces,
                           size_t nf, const float *uv, const uint8_t *movable,
                           double axis_y, double axis_x,
                           const QuadribbonUntangleOpts *opts,
                           QuadribbonUntangleStats *out)
{
    QuadribbonUntangleOpts o;
    QuadribbonUntangleStats s;
    QrGrid g;
    uint8_t *movable_all = NULL;
    float *rest_target = NULL;
    double min_u_sep;
    int rc = -1;
    memset(&s, 0, sizeof s);
    s.attempted = 1;
    if (out) *out = s;
    if (!verts || !faces || !uv || nv == 0 || nf == 0) return -1;
    if (opts) o = *opts; else QuadribbonUntangle_defaults(&o);
    if (qr_grid_build(&g, nv, faces, nf, uv) != 0) {
        fprintf(stderr, "[qr-untangle] cannot build UV grid (degenerate "
                "parameterization?)\n");
        return -1;
    }
    if (!movable) {
        movable_all = (uint8_t *)malloc(nv);
        if (!movable_all) { qr_grid_free(&g); return -1; }
        memset(movable_all, 1, nv);
        movable = movable_all;
    }
    rest_target = (float *)malloc(nv * 3 * sizeof *rest_target);
    if (!rest_target) { free(movable_all); qr_grid_free(&g); return -1; }
    memcpy(rest_target, verts, nv * 3 * sizeof *rest_target);
    min_u_sep = o.minimum_u_separation > 0.0
              ? o.minimum_u_separation
              : fmax(100.0, 64.0 * g.grid_du);
    fprintf(stderr,
            "[qr-untangle] nv=%zu nf=%zu cells=%zu grid_du=%.3f bv=%.3f "
            "min_u_sep=%.1f axis=(%.1f,%.1f)\n",
            nv, nf, g.ncell, g.grid_du, g.bv, min_u_sep, axis_y, axis_x);
    if (qr_repair_turn_order(&g, verts, movable, axis_y, axis_x, min_u_sep,
                             o.wrap_pitch_hint, &s) != 0)
        goto done;
    rc = qr_resolve_self_collisions(&g, verts, movable, rest_target, axis_y,
                                    axis_x, min_u_sep, &o, &s);
    if (s.output_long_conflicts > 0)
        qr_characterize_residual(&g, verts, axis_y, axis_x, min_u_sep);
done:
    free(rest_target);
    free(movable_all);
    qr_grid_free(&g);
    if (out) *out = s;
    return rc;
}

/* ---- selftests ------------------------------------------------------------ */

/* Build a synthetic ribbon: `turns` concentric arc layers, each rows x cols,
 * axis at origin along +Z (zyx frame), row pitch 2 in z, u continuous along
 * the arc with a big per-layer offset (so cross-layer conflicts classify as
 * long-range).  Layer L sits at radius base_r + L * dr, except where `bulge`
 * pushes layer 0 outward past layer 1 to create a genuine interpenetration. */
static void qr_selftest_ribbon(float *verts, float *uv, int32_t *faces,
                               int rows, int cols, int turns, double base_r,
                               double dr, double bulge)
{
    const double arc0 = 0.3, arc1 = 2.7;
    size_t nvl = (size_t)rows * (size_t)cols;
    for (int L = 0; L < turns; L++)
        for (int r = 0; r < rows; r++)
            for (int c = 0; c < cols; c++) {
                size_t v = (size_t)L * nvl + (size_t)r * (size_t)cols +
                           (size_t)c;
                double t = arc0 + (arc1 - arc0) * (double)c /
                           (double)(cols - 1);
                double rad = base_r + (double)L * dr;
                if (L == 0 && bulge > 0.0) {
                    double mid = 0.5 * (arc0 + arc1);
                    double w = exp(-8.0 * (t - mid) * (t - mid));
                    rad += bulge * w;
                }
                verts[v * 3 + 0] = (float)(2 * r);
                verts[v * 3 + 1] = (float)(rad * sin(t));
                verts[v * 3 + 2] = (float)(rad * cos(t));
                uv[v * 2 + 0] = (float)(t * base_r + (double)L * 1000.0);
                uv[v * 2 + 1] = (float)(2 * r);
            }
    {
        size_t nfc = 0;
        for (int L = 0; L < turns; L++)
            for (int r = 0; r + 1 < rows; r++)
                for (int c = 0; c + 1 < cols; c++) {
                    int32_t a = (int32_t)((size_t)L * nvl +
                                          (size_t)r * (size_t)cols +
                                          (size_t)c);
                    int32_t b = a + 1, d = a + cols, e = d + 1;
                    faces[nfc * 3 + 0] = a;
                    faces[nfc * 3 + 1] = b;
                    faces[nfc * 3 + 2] = d;
                    nfc++;
                    faces[nfc * 3 + 0] = b;
                    faces[nfc * 3 + 1] = e;
                    faces[nfc * 3 + 2] = d;
                    nfc++;
                }
    }
}

int QuadribbonUntangle_selftest(void)
{
    int fails = 0;
    enum { ROWS = 10, COLS = 60, TURNS = 2 };
    size_t nv = (size_t)ROWS * COLS * TURNS;
    size_t nf = (size_t)(ROWS - 1) * (COLS - 1) * 2 * TURNS;
    float *verts = (float *)malloc(nv * 3 * sizeof *verts);
    float *uv = (float *)malloc(nv * 2 * sizeof *uv);
    int32_t *faces = (int32_t *)malloc(nf * 3 * sizeof *faces);
    if (!verts || !uv || !faces) {
        free(verts); free(uv); free(faces);
        fprintf(stderr, "[qr-untangle selftest] OOM\n");
        return -1;
    }

    /* (u1) clean two-layer ribbon: no conflicts, zero displacement */
    qr_selftest_ribbon(verts, uv, faces, ROWS, COLS, TURNS, 30.0, 9.5, 0.0);
    {
        QuadribbonUntangleOpts o;
        QuadribbonUntangleStats st;
        float *before = (float *)malloc(nv * 3 * sizeof *before);
        double drift = 0.0;
        QuadribbonUntangle_defaults(&o);
        o.minimum_u_separation = 300.0;
        o.wrap_pitch_hint = 9.5;
        memcpy(before, verts, nv * 3 * sizeof *before);
        if (QuadribbonUntangle_run(verts, nv, faces, nf, uv, NULL, 0.0, 0.0,
                                   &o, &st) != 0) {
            fprintf(stderr, "  FAIL: (u1) clean run rc\n");
            fails++;
        }
        for (size_t i = 0; i < nv * 3; i++) {
            double d = fabs((double)verts[i] - before[i]);
            if (d > drift) drift = d;
        }
        if (drift > 1e-6 || st.input_long_conflicts != 0) {
            fprintf(stderr,
                    "  FAIL: (u1) clean ribbon moved (drift=%.3g long=%zu)\n",
                    drift, st.input_long_conflicts);
            fails++;
        }
        free(before);
    }

    /* (u2) interpenetrating bulge: layer 0 pushes 12 vox through layer 1
     * (pitch 9.5); the shell must remove every long-range conflict and keep
     * orientation clean */
    qr_selftest_ribbon(verts, uv, faces, ROWS, COLS, TURNS, 30.0, 9.5, 12.0);
    {
        QuadribbonUntangleOpts o;
        QuadribbonUntangleStats st;
        QuadribbonUntangle_defaults(&o);
        o.minimum_u_separation = 300.0;
        o.wrap_pitch_hint = 9.5;
        o.collision_rounds = 24;
        o.settle_rounds = 6;
        if (QuadribbonUntangle_run(verts, nv, faces, nf, uv, NULL, 0.0, 0.0,
                                   &o, &st) != 0) {
            fprintf(stderr, "  FAIL: (u2) untangle rc\n");
            fails++;
        }
        fprintf(stderr,
                "[qr-untangle selftest] (u2) turn-order long %zu->%zu, shell "
                "long %zu->%zu complete=%d accepted=%d settle=%d/%d\n",
                st.turn_order_input_long, st.turn_order_output_long,
                st.input_long_conflicts, st.output_long_conflicts,
                st.complete, st.accepted_rounds, st.settle_accepted,
                st.settle_rounds_run);
        if (st.turn_order_input_long == 0 && st.input_long_conflicts == 0) {
            fprintf(stderr,
                    "  FAIL: (u2) fixture created no interpenetration\n");
            fails++;
        }
        if (st.output_long_conflicts != 0) {
            fprintf(stderr,
                    "  FAIL: (u2) untangle left %zu long-range conflicts\n",
                    st.output_long_conflicts);
            fails++;
        }
    }

    free(verts);
    free(uv);
    free(faces);
    fprintf(stderr, "[qr-untangle selftest] %s (%d failure(s))\n",
            fails == 0 ? "PASS" : "FAIL", fails);
    return fails == 0 ? 0 : -1;
}
