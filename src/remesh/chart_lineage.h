#ifndef CHART_LINEAGE_INCLUDED
#define CHART_LINEAGE_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"

/* Persistent chart identity across topology-preserving vertex-fan splits.
 *
 * Pinch repair duplicates vertices and therefore turns one pre-split chart into
 * several ordinary connected components.  Recomputing component labels after
 * that operation loses the fact that those fragments came from one chart and
 * prevents legitimate boundary zippers from rejoining them.  ChartLineage
 * captures the component partition immediately before the split and recovers it
 * later from bit-identical source coordinates.  New fill vertices inherit the
 * unique lineage of their current component. */
typedef struct {
    void *points;                 /* private sorted coordinate table */
    size_t n_points;
    size_t n_components;
} ChartLineage;

typedef struct {
    size_t current_components;
    size_t components_with_lineage;
    size_t ambiguous_components;
    size_t matched_vertices;
} ChartLineageStats;

int ChartLineage_capture(Arena_T arena,
                         const float *verts, size_t nv,
                         const int32_t *faces, size_t nf,
                         ChartLineage *out);

int ChartLineage_assign(Arena_T arena, const ChartLineage *lineage,
                        const float *verts, size_t nv,
                        const int32_t *faces, size_t nf,
                        int32_t **out_merge_group,
                        ChartLineageStats *stats);

int ChartLineage_selftest(void);

#endif
