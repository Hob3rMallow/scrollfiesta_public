#ifndef DEV_METRIC_INCLUDED
#define DEV_METRIC_INCLUDED

#include <stddef.h>
#include <stdint.h>

/* Discrete developability metric: per-vertex angle-defect Gaussian curvature.
 *
 *   K_v = (boundary ? PI : 2*PI) - sum(incident triangle angles at v)
 *
 * A developable surface (e.g. a papyrus wrap) has K_v ~ 0 at every interior
 * vertex; cone/saddle points (bridges, tangles, folds) concentrate |K|.
 *
 * This is the single definition of "developability" shared by the QC tool
 * (developability.exe), the end-of-pipeline gate (dev_gate), and any future
 * caller. Uses malloc for scratch -- COLD paths only (QC / end-of-pipeline),
 * never a hot inner loop. */

/* Fill per-vertex arrays. Caller allocates defect[nv] (double),
 * boundary[nv] and used[nv] (unsigned char each). A vertex is `boundary` if it
 * touches an edge used by exactly one triangle; `used` if any face references
 * it. Unused vertices get defect 0. */
void DevMetric_compute(const float *verts, size_t nv,
                       const int32_t *faces, size_t nf,
                       double *defect, unsigned char *boundary,
                       unsigned char *used);

typedef struct {
    size_t n_interior;          /* used, non-boundary vertices               */
    size_t n_interior_bad;      /* of those, with |K| > k_thresh             */
    size_t n_boundary;          /* used boundary vertices (K excluded)       */
    double max_abs_k;           /* max |K| over interior vertices (rad)      */
    double sum_abs_k_interior;  /* total absolute Gaussian curvature (rad)   */
    double frac_interior_bad;   /* n_interior_bad / n_interior (0 if none)   */
} DevMetricStats;

/* Convenience: compute interior-vertex developability stats in one call.
 * Boundary vertices are excluded (their angle defect is ill-defined). */
void DevMetric_interior_stats(const float *verts, size_t nv,
                              const int32_t *faces, size_t nf,
                              double k_thresh, DevMetricStats *out);

#endif
