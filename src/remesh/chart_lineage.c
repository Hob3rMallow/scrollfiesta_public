/* chart_lineage.c -- recover pre-split chart provenance by exact position. */
#include "chart_lineage.h"

#include "sheet_reweld.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    uint32_t bits[3];
    int32_t component;            /* -1 means coordinate was ambiguous */
} ChartPoint;

static uint32_t chart_float_bits(float value)
{
    uint32_t bits;
    memcpy(&bits, &value, sizeof bits);
    if ((bits << 1) == 0) bits = 0;   /* canonicalize -0 */
    return bits;
}

static int chart_point_cmp(const void *pa, const void *pb)
{
    const ChartPoint *a = (const ChartPoint *)pa;
    const ChartPoint *b = (const ChartPoint *)pb;
    for (int k = 0; k < 3; k++) {
        if (a->bits[k] != b->bits[k])
            return a->bits[k] < b->bits[k] ? -1 : 1;
    }
    if (a->component != b->component)
        return a->component < b->component ? -1 : 1;
    return 0;
}

static int chart_same_position(const ChartPoint *a, const ChartPoint *b)
{
    return a->bits[0] == b->bits[0] &&
           a->bits[1] == b->bits[1] &&
           a->bits[2] == b->bits[2];
}

static int32_t chart_lookup(const ChartLineage *lineage, const float *p)
{
    const ChartPoint *points = (const ChartPoint *)lineage->points;
    uint32_t key[3] = {
        chart_float_bits(p[0]), chart_float_bits(p[1]),
        chart_float_bits(p[2])
    };
    size_t lo = 0, hi = lineage->n_points;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        int cmp = 0;
        for (int k = 0; k < 3; k++) {
            if (points[mid].bits[k] != key[k]) {
                cmp = points[mid].bits[k] < key[k] ? -1 : 1;
                break;
            }
        }
        if (cmp < 0) lo = mid + 1;
        else hi = mid;
    }
    if (lo >= lineage->n_points) return -1;
    for (int k = 0; k < 3; k++)
        if (points[lo].bits[k] != key[k]) return -1;
    return points[lo].component;
}

int ChartLineage_capture(Arena_T arena,
                         const float *verts, size_t nv,
                         const int32_t *faces, size_t nf,
                         ChartLineage *out)
{
    int32_t *component = NULL;
    size_t ncomponents = 0, npoints = 0;
    ChartPoint *points;
    if (!arena || !verts || !faces || !out || nv == 0 || nf == 0)
        return -1;
    memset(out, 0, sizeof *out);
    if (SheetReweld_label(arena, faces, nf, nv,
                          &component, &ncomponents) != 0)
        return -1;
    points = (ChartPoint *)ARENA_ALLOC(
        arena, ((nv ? nv : 1) * sizeof(ChartPoint)));
    for (size_t v = 0; v < nv; v++) {
        if (component[v] < 0) continue;
        for (int k = 0; k < 3; k++)
            points[npoints].bits[k] = chart_float_bits(verts[v*3+(size_t)k]);
        points[npoints].component = component[v];
        npoints++;
    }
    qsort(points, npoints, sizeof *points, chart_point_cmp);
    {
        size_t write = 0;
        for (size_t i = 0; i < npoints; ) {
            size_t j = i + 1;
            int32_t group = points[i].component;
            while (j < npoints && chart_same_position(&points[i], &points[j])) {
                if (points[j].component != group) group = -1;
                j++;
            }
            points[write] = points[i];
            points[write].component = group;
            write++;
            i = j;
        }
        npoints = write;
    }
    out->points = points;
    out->n_points = npoints;
    out->n_components = ncomponents;
    return 0;
}

int ChartLineage_assign(Arena_T arena, const ChartLineage *lineage,
                        const float *verts, size_t nv,
                        const int32_t *faces, size_t nf,
                        int32_t **out_merge_group,
                        ChartLineageStats *stats)
{
    int32_t *current = NULL, *component_lineage, *group;
    uint8_t *conflict;
    size_t ncurrent = 0;
    ChartLineageStats st;
    if (!arena || !lineage || !lineage->points || !verts || !faces ||
        !out_merge_group || nv == 0 || nf == 0)
        return -1;
    memset(&st, 0, sizeof st);
    if (SheetReweld_label(arena, faces, nf, nv,
                          &current, &ncurrent) != 0)
        return -1;
    component_lineage = (int32_t *)ARENA_ALLOC(
        arena, ((ncurrent ? ncurrent : 1) * sizeof(int32_t)));
    conflict = (uint8_t *)ARENA_CALLOC(
        arena, (ncurrent ? ncurrent : 1), 1L);
    for (size_t c = 0; c < ncurrent; c++) component_lineage[c] = -1;
    for (size_t v = 0; v < nv; v++) {
        int32_t c = current[v];
        int32_t source;
        if (c < 0) continue;
        source = chart_lookup(lineage, &verts[v*3]);
        if (source < 0) continue;
        st.matched_vertices++;
        if (component_lineage[c] < 0) component_lineage[c] = source;
        else if (component_lineage[c] != source) conflict[c] = 1;
    }
    group = (int32_t *)ARENA_ALLOC(
        arena, (nv * sizeof(int32_t)));
    for (size_t c = 0; c < ncurrent; c++) {
        if (conflict[c]) {
            component_lineage[c] = -1;
            st.ambiguous_components++;
        } else if (component_lineage[c] >= 0) {
            st.components_with_lineage++;
        }
    }
    for (size_t v = 0; v < nv; v++) {
        int32_t c = current[v];
        if (c < 0) group[v] = -1;
        else if (component_lineage[c] >= 0) group[v] = component_lineage[c];
        else group[v] = (int32_t)(lineage->n_components + (size_t)c);
    }
    st.current_components = ncurrent;
    *out_merge_group = group;
    if (stats) *stats = st;
    return 0;
}

int ChartLineage_selftest(void)
{
    /* Two fans meet at vertex 0 before repair.  Splitting that vertex produces
     * two disconnected charts, but both must retain lineage 0. */
    static const float before_v[] = {
        0,0,0,  0,1,0,  0,0,1,  0,-1,0,  0,0,-1
    };
    static const int32_t before_f[] = { 0,1,2,  0,3,4 };
    static const float after_v[] = {
        0,0,0,  0,1,0,  0,0,1,  0,-1,0,  0,0,-1,  0,0,0
    };
    static const int32_t after_f[] = { 0,1,2,  5,3,4 };
    Arena_T arena = Arena_new();
    ChartLineage lineage;
    ChartLineageStats st;
    int32_t *group = NULL;
    int fail = 0;
    if (ChartLineage_capture(arena, before_v, 5, before_f, 2,
                             &lineage) != 0 ||
        ChartLineage_assign(arena, &lineage, after_v, 6, after_f, 2,
                            &group, &st) != 0 ||
        group[0] != group[5] || st.current_components != 2 ||
        st.components_with_lineage != 2) {
        fprintf(stderr, "[selftest] chart lineage split recovery -> FAIL\n");
        fail = 1;
    } else {
        fprintf(stderr, "[selftest] chart lineage split recovery -> ok\n");
    }
    Arena_dispose(&arena);
    return fail;
}
