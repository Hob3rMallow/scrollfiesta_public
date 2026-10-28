#ifndef CLOSEUP_LOCAL_GLOBAL_INCLUDED
#define CLOSEUP_LOCAL_GLOBAL_INCLUDED

#include "closeup_angular_metric_angular.h"

/* Positive-handed, centered complete material domain; unchanged XYZ/faces.
 * reference_uv defines the immutable 25% area/count budgets. seed must already
 * satisfy the per-group 5% bad10-area budget. Caller owns restoration to the
 * packed frame and complete-sheet material checks through material_guard.
 * Returns the last fully checked iterate on a later numerical refusal; status
 * and termination explicitly distinguish that from convergence. */
int CloseupLocalGlobal_run(Arena_T arena,const float *xyz,size_t nv,
    const int32_t *faces,size_t nf,const double *seed,const PMAOptions *options,
    double *output,PMAStats *stats);
int CloseupLocalGlobal_selftest(void);

#endif
