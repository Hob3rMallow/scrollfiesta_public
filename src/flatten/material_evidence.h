#ifndef MATERIAL_EVIDENCE_INCLUDED
#define MATERIAL_EVIDENCE_INCLUDED

#include "winding_mrf.h"

/* Physical interval evidence, not vertex/triangle vote counts. unit identifies
 * a canonical support interval in the SOURCE hierarchy (not an output face).
 * begin/end are positions along that interval in voxels. Refinement must retain
 * the unit and split these intervals; duplicate observations must not add mass.
 * Units must be disjoint physical support; assigning every sample a fresh unit
 * would violate this API's contract. The caller records the unit provenance.
 *
 * A unit contributes at most one effective observation, shared by competing
 * targets. length_scale is the physical support needed for one full unit.
 * quality is a bounded observation score, NOT a calibrated probability. */
typedef struct {
    int32_t a, b, target;
    int kind;                   /* WINDING_MRF_EQUAL / AT_LEAST / AT_MOST */
    uint64_t unit;
    double begin, end, quality;
    uint64_t source_face_a, source_face_b; /* ORIGINAL face provenance; retained
                                           * through subdivision/coarsening */
    uint64_t source_edge_a, source_edge_b;
    uint64_t region_a, region_b; /* physical correlation blocks on BOTH sources;
                                  * retained through subdivision, never per sample */
    double begin_b, end_b;      /* paired source-edge coordinates; may decrease */
} MaterialEvidenceInterval;

typedef struct {
    WindingMRFEdge edge;
    int kind;
    size_t units;
    size_t regions_a, regions_b; /* units = min(distinct blocks on each side) */
    double covered_length, effective_support;
    double pair_margin;        /* least local pair-energy increase for +/-1;
                                  * zero for one-sided evidence, not a posterior */
    int significant;            /* publication evidence, not a hard constraint */
} MaterialEvidenceFactor;

typedef struct {
    size_t intervals, units, factors, significant_factors;
    double covered_length, effective_support;
} MaterialEvidenceReport;

/* Exact upper envelope of interval quality within each (pair,kind,unit,target).
 * Sum targets with a SHARED per-unit budget, then cap total influence to one
 * per original face AND physical correlation block on EACH side before
 * accumulating by hypothesis. A triangle
 * cannot gain influence by contacting many triangles, units or hypotheses.
 * The per-source-face cap bounds influence for a single corrupted original
 * triangle. It is a robustness device, not a claim of statistical independence
 * between neighboring faces. Blocks bound correlated influence; source-scale
 * calibration remains required. units counts distinct physical blocks, not
 * source triangles or generated interval records. The less-supported side
 * limits significance. Every mass contributes to both sides' shared budgets.
 * For equality evidence, significant also requires a local pair-energy margin
 * greater than twice one bounded unit, against BOTH neighboring labels. For
 * these reduced L1 scores, an adversarial replacement of one unit cannot then
 * change this pair's preferred target. This is NOT a global graph guarantee,
 * a calibrated probability, or permission to accept a material interface.
 * All nonzero hypotheses survive as soft MRF factors, even when insufficient
 * to publish a continuation. A local contradiction never vetoes a whole graph.
 * Inputs unmodified; outputs arena-owned. Units cannot change pair/family.
 * Return -1 on malformed/overflowing input, with *out=NULL and *count=0. */
int MaterialEvidence_reduce(Arena_T arena,
                            const MaterialEvidenceInterval *intervals, size_t n,
                            double length_scale, MaterialEvidenceFactor **out,
                            size_t *count, MaterialEvidenceReport *report);

int MaterialEvidence_selftest(void);

#endif
