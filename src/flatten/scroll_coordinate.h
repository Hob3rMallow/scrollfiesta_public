#ifndef SCROLL_COORDINATE_INCLUDED
#define SCROLL_COORDINATE_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"

/* Integer chart synchronization for a FIXED scroll model. Node indices are
 * assigned by stable source-chart order before this function is called.
 * Every edge says winding[b] - winding[a] = delta. Exact edges are contracted
 * with integer potentials; a contradictory exact cycle fails the solve.
 * Soft edges enter one graph least-squares problem, including non-tree edges.
 * No query rectangle, parent chart, or crop-local origin enters this API. */
typedef struct {
    int32_t a, b;
    int32_t delta;
    double weight;               /* finite and positive for soft evidence */
    int exact;                   /* exact correspondence, not high confidence */
} ScrollCoordinateEdge;

typedef struct {
    size_t nodes, exact_edges, soft_edges;
    size_t exact_groups, gauge_components;
    size_t contradictory_exact_edges, unsatisfied_soft_edges;
    double soft_squared_error;
    int64_t max_soft_residual;
} ScrollCoordinateReport;

/* Output is arena-owned, with one integer per node. The least stable-index
 * node in each disconnected graph component fixes that component's zero.
 * On failure *out_winding is NULL; report remains available for diagnosis.
 * All hard equalities are verified again after the numeric solve. Soft
 * residuals are exposed even when the solve succeeds: callers must apply
 * their physical-evidence acceptance policy before publishing a model.
 * Empty input succeeds with a NULL output. No input arrays are modified. */
int ScrollCoordinate_solve(Arena_T arena, size_t nnodes,
                           const ScrollCoordinateEdge *edges, size_t nedges,
                           int64_t **out_winding,
                           ScrollCoordinateReport *report);

int ScrollCoordinate_selftest(void);

#endif
