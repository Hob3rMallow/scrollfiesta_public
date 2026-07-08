/* dev_gate.c -- end-of-pipeline developability "bad sheet" gate. See dev_gate.h. */
#include "dev_gate.h"

#include "../common/dev_metric.h"
#include "../common/pipeline_constants.h"
#include "../split/oracle.h"

#include <string.h>

void DevGate_default_params(DevGateParams *p)
{
    p->k_thresh    = 0.5;                 /* rad (~29 deg) per-vertex spike      */
    p->frac_thresh = 0.005;               /* 0.5% of interior verts             */
    p->min_bad     = 10;                  /* absolute floor so tiny sheets pass  */
    p->oracle_grid = ORACLE_UV_GRID_SIZE; /* same grid the splitter oracle uses  */
}

int DevGate_classify(Arena_T arena, const ComponentMesh *cm,
                     const DevGateParams *p, DevGateVerdict *v)
{
    memset(v, 0, sizeof(*v));
    v->oracle_sheets = -1;
    if (!cm || cm->nv == 0 || cm->nf == 0) return 0;

    /* Criterion 1: non-developable interior. */
    DevMetricStats s;
    DevMetric_interior_stats(cm->verts, cm->nv, cm->faces, cm->nf, p->k_thresh, &s);
    v->n_interior        = s.n_interior;
    v->n_interior_bad    = s.n_interior_bad;
    v->frac_interior_bad = s.frac_interior_bad;
    v->max_abs_k         = s.max_abs_k;
    if (s.n_interior_bad >= p->min_bad && s.frac_interior_bad > p->frac_thresh)
        v->reason |= DEVGATE_NONDEV;

    /* Criterion 2: oracle still sees multiple sheets after the split round. */
    if (p->oracle_grid > 0) {
        Arena_Mark mark = Arena_save(arena);
        v->oracle_sheets = Oracle_count_sheets(arena, cm, p->oracle_grid);
        Arena_restore(arena, mark);
        if (v->oracle_sheets > 1) v->reason |= DEVGATE_ORACLE;
    }

    v->bad = v->reason ? 1 : 0;
    return 0;
}
