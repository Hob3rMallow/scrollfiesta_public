#ifndef MATERIAL_FRONT_INCLUDED
#define MATERIAL_FRONT_INCLUDED

#include "material_evidence.h"
#include "scroll_source.h"

/* Original boundary segments are an evidence representation, not stitched
 * geometry. face is a stable original-source face identity; descendants keep
 * it. source_edge identifies the oriented original halfedge for diagnostics. */
typedef struct {
    float xyz[6], normal[3], outward[3];
    double phase[2];
    uint64_t face, source_edge;
    int32_t chart;
    int32_t source_vertices[2];
    uint64_t curve, region_base;
    double arc;                /* source-boundary arclength at xyz[0..2] */
} MaterialFront;

typedef struct {
    size_t fronts, candidates, supported_intervals, contacts, ambiguous_phase;
    size_t degenerate_intervals; /* projected support no larger than roundoff */
    double supported_length;
} MaterialFrontReport;

/* Extract one original cube. source_prefix occupies the high 32 bits of face
 * and halfedge identities. chart_base is a run-local graph index, never UV.
 * Nonmanifold/orientation-invalid front vertices cannot propose continuation.
 * Core/cyclic classification is retained by the observation, not deleted here. */
int MaterialFront_extract(Arena_T arena, const MeshBinData *mesh,
                          const ScrollSource *source, uint32_t source_prefix,
                          int32_t chart_base, MaterialFront **out, size_t *count);

/* Trace original oriented boundary topology, without geometric joining.
 * Assign stable per-cube curve anchors, arclength and physical block ranges.
 * Open chains start at their source endpoint; loops at their least halfedge.
 * Input must be ONE source cube with unique oriented edges and nonbranching
 * boundary vertices. Permuting records does not change the assigned frame.
 * Modifies only curve/region_base/arc. Extract calls this automatically. */
int MaterialFront_trace(Arena_T arena, MaterialFront *fronts, size_t n);

/* Measure the continuous projected overlap of mutually facing segments, with
 * endpoints and midpoint tested. BVH broadphase has no neighbor-prefix cap.
 * Does NOT choose a winner, weld edges, or establish material identity.
 * Bounds are physical, independent of vertex density. Intervals keep original
 * face provenance and pass through MaterialEvidence_reduce before the MRF.
 * Inputs must have been traced. Overlaps are split at BOTH sources' physical
 * block boundaries, retaining the monotone source-edge correspondence. */
int MaterialFront_collect(Arena_T arena, const MaterialFront *fronts, size_t n,
                          MaterialEvidenceInterval **out, size_t *count,
                          MaterialFrontReport *report);

int MaterialFront_selftest(void);

#endif
