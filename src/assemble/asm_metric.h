#ifndef ASM_METRIC_INCLUDED
#define ASM_METRIC_INCLUDED
#include "asm_types.h"
#include "../common/pipeline_constants.h"

/* Shared by admission, deformation, encoding and independent auditing.
 * Keep raw singular values in reports; only comparisons absorb roundoff. */
static inline int AsmMetric_in_band(double lo,double hi,double lower,double upper)
{ return lo<=hi+ASM_METRIC_ROUNDOFF && lo>=lower-ASM_METRIC_ROUNDOFF && hi<=upper+ASM_METRIC_ROUNDOFF; }
static inline int AsmMetric_within10(double lo,double hi)
{ return AsmMetric_in_band(lo,hi,.90,1.10); }
static inline int AsmMetric_within25(double lo,double hi)
{ return AsmMetric_in_band(lo,hi,.75,1.25); }

/* Immutable original-triangle metric. The second coordinate uses the norm of
 * a 3-D cross product, avoiding cancellation in a thin triangle's Gram matrix. */
typedef struct AsmMetricFace { double length, along, height, area; } AsmMetricFace;
typedef struct AsmMetricStats {
    double area, within10, within25, minimum, maximum;
    size_t invalid;
} AsmMetricStats;

/* Caller owns ref[nf]. Source geometry and UVs are never modified. */
int AsmMetric_prepare(const AsmChart *chart, AsmMetricFace *ref);
void AsmMetric_measure(const AsmChart *chart, const AsmMetricFace *ref,
                       const float *uv, AsmMetricStats *out);
/* Measure the actual encoded global coordinates. The layout's explicit
 * mirror parity is allowed; additional flips or collapsed faces are not. */
void AsmMetric_measure_encoded(const AsmChart *chart, const AsmMetricFace *ref,
                               const float *uv, AsmMetricStats *out);
void AsmMetric_measure_encoded_uv64(const AsmChart *chart, const AsmMetricFace *ref,
                                    const double *uv, AsmMetricStats *out);
/* Preserve passing original-chart limits. Existing failures may improve or
 * remain unchanged; a large neighbour cannot supply another chart's budget. */
int AsmMetric_preserved(const AsmMetricStats *before, const AsmMetricStats *after);
/* Guard a local-UV face throughout a deformation, including the individual
 * strict metric certificate that chart area fractions cannot replace. */
int AsmMetric_local_face_preserved(const AsmChart *chart, size_t face,
                                  const AsmMetricFace *ref, const float *trial);
/* Every previously strict-valid original face must remain strict-valid in the
 * actual encoded field; a chart's area fractions cannot hide a lost face. */
int AsmMetric_encoded_faces_preserved(const AsmChart *chart, const AsmMetricFace *ref, const double *encoded);
/* The same guard when the encoded coordinates carry a different mirror parity
 * than the chart's current one (a whole-piece presentation reflection that
 * also toggles ASM_CHART_MIRROR): encoded_mirror is the parity they must show. */
int AsmMetric_encoded_faces_preserved_parity(const AsmChart *chart, const AsmMetricFace *ref, const double *encoded,
                                             int encoded_mirror);
int AsmMetric_selftest(void);
#endif
