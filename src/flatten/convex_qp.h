#ifndef VES_CONVEX_QP_INCLUDED
#define VES_CONVEX_QP_INCLUDED
#include <stddef.h>
/* Same equilibrated OSQP proposal as c3_seam_restoration.solve_restoration_qp.
 * H is lower COO; A is general COO with upper bounds. Duplicate entries sum.
 * Return 1 for a finite solved/inaccurate/iteration-limited proposal. The
 * caller must check its actual constraints and original geometry. */
int ConvexQp_solve(int n,int hn,const int *hr,const int *hc,const double *hv,
    const double *target,int m,int an,const int *ar,const int *ac,const double *av,
    const double *upper,double *result);
/* Whole-chart initialization matches translated_hull_seed's unscaled OSQP
 * units, tolerances and 20,000-iteration budget; a capped proposal is refused. */
int ConvexQp_translation_seed(int n,int hn,const int *hr,const int *hc,const double *hv,
    const double *target,int m,int an,const int *ar,const int *ac,const double *av,
    const double *upper,double *result);
/* A wall-clock limit for every following solve (seconds; 0 = none).  A
 * solve that reaches it is reported as not solved. */
void ConvexQp_time_limit(double seconds);
/* Absolute monotonic deadline shared across nested/fallback solves; 0 disables.
 * Process-local, like time_limit: callers must serialize QP use. */
void ConvexQp_deadline(double deadline);
int ConvexQp_selftest(void);
#endif
