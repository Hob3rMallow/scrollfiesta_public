#ifndef DEV_GATE_INCLUDED
#define DEV_GATE_INCLUDED

#include "../common/arena.h"
#include "../common/mesh_types.h"

/* Developability gate -- end-of-pipeline "bad sheet" detector.
 *
 * A finished component (sheet) is flagged BAD if EITHER:
 *   (1) it has a significant number of non-developable INTERIOR vertices
 *       (|angle-defect K| > k_thresh on >= min_bad interior verts AND a fraction
 *       > frac_thresh). Boundary vertices are excluded -- their K is ill-defined.
 *   (2) the UV-raycast oracle still reports more than one sheet after the
 *       splitters have run (i.e. a tangle the bridge/overlap cut could not
 *       separate). Set oracle_grid <= 0 to skip this criterion.
 *
 * A clean papyrus wrap is developable and single-sheet, so it passes both. */

enum { DEVGATE_NONDEV = 1, DEVGATE_ORACLE = 2 };

typedef struct {
    double k_thresh;     /* |K| (rad) above which an interior vertex is non-developable */
    double frac_thresh;  /* flag if (bad interior verts / interior verts) exceeds this   */
    size_t min_bad;      /* ...and at least this many bad interior verts (absolute floor) */
    int    oracle_grid;  /* UV grid size for criterion 2; <=0 disables it               */
} DevGateParams;

typedef struct {
    int    bad;                 /* 1 if flagged */
    int    reason;              /* bitmask of DEVGATE_*                  */
    size_t n_interior;
    size_t n_interior_bad;
    double frac_interior_bad;
    double max_abs_k;
    int    oracle_sheets;       /* -1 if not run */
} DevGateVerdict;

void DevGate_default_params(DevGateParams *p);

/* Classify one component. arena state is unchanged on return (oracle scratch is
 * saved/restored internally). Returns 0. */
int DevGate_classify(Arena_T arena, const ComponentMesh *cm,
                     const DevGateParams *p, DevGateVerdict *v);

#endif
