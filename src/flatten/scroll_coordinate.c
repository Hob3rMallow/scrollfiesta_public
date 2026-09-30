#include "scroll_coordinate.h"
#include "sparse_solve.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* A node's offset is its integer potential relative to its parent. Always
 * attaching the larger root to the smaller root gives a source-stable gauge.
 * Path compression is iterative and never depends on thread scheduling. */
typedef struct {
    int32_t *parent;
    int64_t *offset;
} ScForest;

static int sc_add(int64_t a, int64_t b, int64_t *out)
{
    if ((b > 0 && a > INT64_MAX - b) ||
        (b < 0 && a < INT64_MIN - b)) return -1;
    *out = a + b;
    return 0;
}

static int sc_find(ScForest *f, int32_t v, int32_t *root, int64_t *offset)
{
    int32_t r = v;
    int64_t total = 0;
    while (f->parent[r] != r) {
        if (sc_add(total, f->offset[r], &total) != 0) return -1;
        r = f->parent[r];
    }
    *root = r;
    *offset = total;
    while (f->parent[v] != v) {
        int32_t next = f->parent[v];
        int64_t step = f->offset[v];
        f->parent[v] = r;
        f->offset[v] = total;
        if (step == INT64_MIN || sc_add(total, -step, &total) != 0) return -1;
        v = next;
    }
    return 0;
}

/* Return 1 for a contradictory cycle, 0 for a compatible equality. */
static int sc_join(ScForest *f, int32_t a, int32_t b, int64_t delta)
{
    int32_t ra = 0, rb = 0;
    int64_t da = 0, db = 0, relative = 0;
    if (sc_find(f, a, &ra, &da) != 0 ||
        sc_find(f, b, &rb, &db) != 0 || db == INT64_MIN ||
        sc_add(delta, da, &relative) != 0 ||
        sc_add(relative, -db, &relative) != 0) return -1;
    if (ra == rb) return relative != 0;
    if (ra < rb) {
        f->parent[rb] = ra;
        f->offset[rb] = relative;
    } else {
        if (relative == INT64_MIN) return -1;
        f->parent[ra] = rb;
        f->offset[ra] = -relative;
    }
    return 0;
}

static int sc_compare(const void *pa, const void *pb)
{
    const ScrollCoordinateEdge *a = (const ScrollCoordinateEdge *)pa;
    const ScrollCoordinateEdge *b = (const ScrollCoordinateEdge *)pb;
    if (a->exact != b->exact) return a->exact > b->exact ? -1 : 1;
    if (a->a != b->a) return a->a < b->a ? -1 : 1;
    if (a->b != b->b) return a->b < b->b ? -1 : 1;
    if (a->delta != b->delta) return a->delta < b->delta ? -1 : 1;
    if (a->weight != b->weight) return a->weight < b->weight ? -1 : 1;
    return 0;
}

int ScrollCoordinate_solve(Arena_T arena, size_t nnodes,
                           const ScrollCoordinateEdge *edges, size_t nedges,
                           int64_t **out_winding,
                           ScrollCoordinateReport *report)
{
    ScrollCoordinateReport local;
    ScrollCoordinateEdge *sorted = NULL;
    int64_t *answer = NULL, *relative = NULL;
    int32_t *group = NULL, *group_of_root = NULL, *unknown = NULL;
    ScForest exact = {0}, graph = {0};
    size_t ngroups = 0, nunknown = 0;
    int rc = -1;
    if (out_winding == NULL) return -1;
    *out_winding = NULL;
    if (report == NULL) report = &local;
    memset(report, 0, sizeof *report);
    report->nodes = nnodes;
    if (arena == NULL || nnodes > (size_t)INT32_MAX ||
        nedges > (size_t)INT_MAX / 3 || (nedges > 0 && edges == NULL)) return -1;
    if (nnodes == 0) return nedges == 0 ? 0 : -1;

    /* Allocate the result before scratch so only the compact answer survives. */
    answer = (int64_t *)ARENA_ALLOC(arena, nnodes * sizeof *answer);
    Arena_Mark mark = Arena_save(arena);
    sorted = (ScrollCoordinateEdge *)ARENA_ALLOC(
        arena, (nedges > 0 ? nedges : 1) * sizeof *sorted);
    exact.parent = (int32_t *)ARENA_ALLOC(arena, nnodes * sizeof(int32_t));
    exact.offset = (int64_t *)ARENA_CALLOC(arena, nnodes, sizeof(int64_t));
    relative = (int64_t *)ARENA_ALLOC(arena, nnodes * sizeof *relative);
    group = (int32_t *)ARENA_ALLOC(arena, nnodes * sizeof *group);
    group_of_root = (int32_t *)ARENA_ALLOC(arena, nnodes * sizeof *group_of_root);
    for (size_t i = 0; i < nnodes; i++) {
        exact.parent[i] = (int32_t)i;
        group_of_root[i] = -1;
    }
    for (size_t i = 0; i < nedges; i++) {
        ScrollCoordinateEdge e = edges[i];
        if (e.a < 0 || e.b < 0 || (size_t)e.a >= nnodes ||
            (size_t)e.b >= nnodes || !isfinite(e.weight) ||
            !(e.weight > 0.0) || e.delta == INT32_MIN) goto done;
        e.exact = e.exact != 0;
        if (e.a > e.b) {
            int32_t swap = e.a;
            e.a = e.b; e.b = swap; e.delta = -e.delta;
        }
        sorted[i] = e;
    }
    qsort(sorted, nedges, sizeof *sorted, sc_compare);
    for (size_t i = 0; i < nedges; i++) {
        const ScrollCoordinateEdge *e = &sorted[i];
        if (!e->exact) { report->soft_edges++; continue; }
        report->exact_edges++;
        int joined = sc_join(&exact, e->a, e->b, (int64_t)e->delta);
        if (joined < 0) goto done;
        if (joined > 0) report->contradictory_exact_edges++;
    }
    if (report->contradictory_exact_edges > 0) goto done;
    for (size_t i = 0; i < nnodes; i++) {
        int32_t root = 0;
        if (sc_find(&exact, (int32_t)i, &root, &relative[i]) != 0) goto done;
        if (group_of_root[root] < 0) group_of_root[root] = (int32_t)ngroups++;
        group[i] = group_of_root[root];
    }
    report->exact_groups = ngroups;
    graph.parent = (int32_t *)ARENA_ALLOC(arena, ngroups * sizeof(int32_t));
    graph.offset = (int64_t *)ARENA_CALLOC(arena, ngroups, sizeof(int64_t));
    unknown = (int32_t *)ARENA_ALLOC(arena, ngroups * sizeof *unknown);
    for (size_t i = 0; i < ngroups; i++) graph.parent[i] = (int32_t)i;
    for (size_t i = 0; i < nedges; i++) {
        const ScrollCoordinateEdge *e = &sorted[i];
        if (sc_join(&graph, group[e->a], group[e->b], 0) != 0) goto done;
    }
    for (size_t i = 0; i < ngroups; i++) {
        int32_t root = 0;
        int64_t unused = 0;
        if (sc_find(&graph, (int32_t)i, &root, &unused) != 0) goto done;
        if ((int32_t)i == root) {
            unknown[i] = -1;
            report->gauge_components++;
        } else unknown[i] = (int32_t)nunknown++;
    }
    double *x = (double *)ARENA_CALLOC(
        arena, nunknown > 0 ? nunknown : 1, sizeof *x);
    if (nunknown > 0) {
        size_t capacity = report->soft_edges * 3, nt = 0;
        int *rows = (int *)ARENA_ALLOC(arena, capacity * sizeof *rows);
        int *cols = (int *)ARENA_ALLOC(arena, capacity * sizeof *cols);
        double *values = (double *)ARENA_ALLOC(arena, capacity * sizeof *values);
        double *rhs = (double *)ARENA_CALLOC(arena, nunknown, sizeof *rhs);
        for (size_t i = 0; i < nedges; i++) {
            const ScrollCoordinateEdge *e = &sorted[i];
            if (e->exact || group[e->a] == group[e->b]) continue;
            int a = unknown[group[e->a]], b = unknown[group[e->b]];
            double delta = (double)e->delta + (double)relative[e->a]
                         - (double)relative[e->b];
            if (a >= 0) {
                rows[nt] = a; cols[nt] = a; values[nt++] = e->weight;
                rhs[a] -= e->weight * delta;
            }
            if (b >= 0) {
                rows[nt] = b; cols[nt] = b; values[nt++] = e->weight;
                rhs[b] += e->weight * delta;
            }
            if (a >= 0 && b >= 0) {
                rows[nt] = a > b ? a : b;
                cols[nt] = a > b ? b : a;
                values[nt++] = -e->weight;
            }
        }
        if (nt > (size_t)INT_MAX ||
            Sparse_solve_sym((int)nunknown, (int)nt, rows, cols, values,
                             rhs, x, SPARSE_SPD) != 0) goto done;
    }
    for (size_t i = 0; i < nnodes; i++) {
        int index = unknown[group[i]];
        double correction = index < 0 ? 0.0 : x[index];
        /* Keep conversion well inside int64 range; neither casts at 2^63 nor
         * cancellation of huge integer gauges is useful physical evidence. */
        if (!isfinite(correction) || fabs(correction) > 2147483647.0) goto done;
        int64_t rounded = (int64_t)round(correction);
        if (sc_add(relative[i], rounded, &answer[i]) != 0) goto done;
    }
    for (size_t i = 0; i < nedges; i++) {
        const ScrollCoordinateEdge *e = &sorted[i];
        int64_t residual = 0;
        if (answer[e->a] == INT64_MIN ||
            sc_add(answer[e->b], -answer[e->a], &residual) != 0 ||
            sc_add(residual, -(int64_t)e->delta, &residual) != 0 ||
            residual == INT64_MIN) goto done;
        if (e->exact) {
            if (residual != 0) { report->contradictory_exact_edges++; goto done; }
        } else {
            int64_t magnitude = residual < 0 ? -residual : residual;
            if (magnitude > report->max_soft_residual)
                report->max_soft_residual = magnitude;
            if (residual != 0) report->unsatisfied_soft_edges++;
            report->soft_squared_error += e->weight * (double)residual * (double)residual;
        }
    }
    *out_winding = answer;
    rc = 0;
done:
    Arena_restore(arena, mark);
    return rc;
}

static int sc_test(int condition, const char *name)
{
    fprintf(stderr, "[scroll coordinate] %s: %s\n", name, condition ? "PASS" : "FAIL");
    return condition ? 0 : 1;
}

int ScrollCoordinate_selftest(void)
{
    Arena_T arena = Arena_new();
    ScrollCoordinateReport report;
    int64_t *answer = NULL;
    int fails = 0;
    ScrollCoordinateEdge loop[] = {
        {0, 1, 3, 1.0, 1}, {0, 2, 7, 1.0, 1},
        {1, 3, 6, 1.0, 1}, {2, 3, 2, 1.0, 1}
    };
    int rc = ScrollCoordinate_solve(arena, 4, loop, 4, &answer, &report);
    fails += sc_test(rc == 0 && answer && answer[0] == 0 && answer[1] == 3 &&
                     answer[2] == 7 && answer[3] == 9 && report.exact_groups == 1,
                     "all four sides of an exact chart loop");
    loop[3].delta = 3;
    rc = ScrollCoordinate_solve(arena, 4, loop, 4, &answer, &report);
    fails += sc_test(rc != 0 && answer == NULL && report.contradictory_exact_edges > 0,
                     "contradictory non-tree edge cannot be hidden");
    loop[3].delta = 2;
    for (size_t i = 0; i < 4; i++) loop[i].exact = 0;
    rc = ScrollCoordinate_solve(arena, 5, loop, 4, &answer, &report);
    fails += sc_test(rc == 0 && answer && answer[3] == 9 && answer[4] == 0 &&
                     report.gauge_components == 2 && report.unsatisfied_soft_edges == 0,
                     "soft graph solve and disconnected gauge");
    int64_t expected[5] = {0};
    if (rc == 0) memcpy(expected, answer, sizeof expected);
    ScrollCoordinateEdge reversed[4];
    for (size_t i = 0; i < 4; i++) {
        reversed[i] = loop[3 - i];
        int32_t swap = reversed[i].a;
        reversed[i].a = reversed[i].b; reversed[i].b = swap;
        reversed[i].delta = -reversed[i].delta;
    }
    rc = ScrollCoordinate_solve(arena, 5, reversed, 4, &answer, &report);
    fails += sc_test(rc == 0 && answer && memcmp(expected, answer, sizeof expected) == 0,
                     "edge enumeration and orientation independence");
    ScrollCoordinateEdge mixed[] = {
        {0, 1, 5, 1.0, 1}, {1, 2, 3, 100.0, 0},
        {0, 2, 8, 100.0, 0}, {0, 2, 10, 0.1, 0}
    };
    rc = ScrollCoordinate_solve(arena, 3, mixed, 4, &answer, &report);
    fails += sc_test(rc == 0 && answer && answer[1] == 5 && answer[2] == 8 &&
                     report.unsatisfied_soft_edges == 1 && report.max_soft_residual == 2,
                     "hard contraction and visible uncertain evidence residual");
    rc = ScrollCoordinate_solve(arena, 0, NULL, 0, &answer, &report);
    fails += sc_test(rc == 0 && answer == NULL, "empty graph");
    mixed[0].weight = NAN;
    rc = ScrollCoordinate_solve(arena, 3, mixed, 4, &answer, &report);
    fails += sc_test(rc != 0 && answer == NULL, "nonfinite evidence rejected");
    Arena_dispose(&arena);
    return fails;
}
