#ifndef ATLAS_CANDIDATES_INCLUDED
#define ATLAS_CANDIDATES_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "atlas_strip.h"

/*
 * Geometry-first adapter from a disjoint union of retained triangle patches
 * to the coarse ordered-stroke problem.
 *
 * The mesh is intersected with one global family of planes perpendicular to
 * axis_dir.  Plane/mesh intersections become ordered strokes.  Candidate
 * cross-sections use only 3D position and tangent agreement; raw u is used
 * solely to orient each stroke and initialize its independent gauge.
 */

typedef struct {
    double slice_spacing;       /* global axial plane spacing */
    double min_stroke_length;   /* discard shorter intersection fragments */
    double sample_spacing;      /* coarse arclength spacing along strokes */

    double match_radius;        /* in-plane cross-slice candidate radius */
    double match_angle_deg;     /* oriented tangent gate */
    int    max_slice_stride;    /* bridge this many empty axial planes */
    int    max_candidates;      /* alternatives retained per source sample */

    double continuation_radius; /* same-slice endpoint gap gate */
    double continuation_angle_deg;

    /*
     * Pinned spiral calibration for winding-aware admission.  pitch <= 0
     * disables every geometric mechanism (legacy raw-u behavior, and the
     * winding audit is skipped).  The relation gate and audits use only
     * PAIRWISE winding deltas -- local differentials that are common-mode
     * immune to axis wander and second-order in pitch error.  The gauge and
     * prior use the ABSOLUTE radius-anchored winding, which is fragile on
     * eccentric or mis-centred geometry, so it is trust-gated per stroke:
     * where the local spiral model does not fit, the prior stays silent.
     */
    double pitch;        /* vox per turn */
    double spiral_a;     /* u(phi) = a*phi + b*phi^2/(4pi) */
    double spiral_b;
    int    sense;        /* winding sense, +1 or -1 */
    int    wind_gate;    /* reject relations with |dwind| > wind_tol */
    double wind_tol;     /* turns; also the per-stroke prior trust fraction */
    int    geo_gauge;    /* orient/gauge strokes from the geometric winding */
    int    wind_cut;     /* segment slice chains at wrap crossings: on a
                          * coarse mesh whose interior edges exceed the
                          * pitch, one chain can snake across wraps through
                          * in-mesh fusions and carry the collapse inside a
                          * single stroke */
} AtlasCandidateOptions;

void AtlasCandidateOptions_default(AtlasCandidateOptions *opts);

/*
 * Exact provenance of one coarse sample.  The point lies on the input mesh
 * edge mesh_vertex[0] -> mesh_vertex[1] at mesh_t, so the solved coarse value
 * can be imposed on a later triangle-FEM extension without a nearest-vertex
 * approximation.
 */
typedef struct {
    int32_t mesh_vertex[2];
    double  mesh_t;
    int32_t cube;
    int32_t plane;
    int32_t mesh_component; /* retained-triangle connected component */
    double  raw_u;
    double  axial;
    /* Geometry about the calibrated axis; zero when pitch <= 0. */
    double  radius;         /* in-plane radius */
    double  theta;          /* in-plane angle, driver reference frame */
    double  wind;           /* sense*r/pitch - theta/2pi (continuous) */
    double  u_geo;          /* spiral_u at the radius-anchored winding */
} AtlasCandidateSampleRef;

typedef enum {
    ATLAS_CANDIDATE_INACTIVE = 0,
    ATLAS_CANDIDATE_SOURCE = 1,
    ATLAS_CANDIDATE_BASELINE = 2,
    ATLAS_CANDIDATE_TOPOLOGY = 3,
    ATLAS_CANDIDATE_REFINED_DELAMINATION = 4,
    ATLAS_CANDIDATE_REFINED_SHEET_SEAM = 5
} AtlasCandidateState;

typedef enum {
    ATLAS_BUNDLE_UNCLASSIFIED = 0,
    ATLAS_BUNDLE_SEPARATE_SHEET = 1,
    ATLAS_BUNDLE_DELAMINATION = 2,
    ATLAS_BUNDLE_SAME_SHEET_SEAM = 3
} AtlasCandidateBundleDecision;

/* One refinement hypothesis joins two retained-triangle components. */
typedef struct {
    int32_t component0;
    int32_t component1;
    size_t links;
    size_t cross_links;
    size_t continuation_links;
    size_t planes;
    int32_t plane_min;
    int32_t plane_max;

    int32_t component0_lifespan;
    int32_t component1_lifespan;
    int32_t component0_extension;
    int32_t component1_extension;
    double contact_overlap_fraction;

    double initial_residual_median;
    double initial_residual_mad;
    double frozen_residual_median;
    double frozen_residual_mad;
    double prior_median;
    double evidence;
    double top_prior_fraction;
    double topology_alternative_fraction;
    double arc_correlation;
    double first_plane_gap;
    double last_plane_gap;
    double component0_runner_up_ratio;
    double component1_runner_up_ratio;

    int mutual_best;
    AtlasCandidateBundleDecision decision;
} AtlasCandidateBundle;

typedef struct {
    double support_floor;
    double gross_initial_residual;
    double initial_coherence_mad_limit;
    double overlap_residual_limit;
    double overlap_coherence_mad_limit;
    double persistent_contact_fraction;
    double minimum_topology_alternative_fraction;
    double ambiguity_ratio_limit;
    double minimum_prior;
    double minimum_top_prior_fraction;
    double minimum_arc_correlation;
    int minimum_links;
    int minimum_planes;

    int front_minimum_links;
    int front_minimum_planes;
    int front_maximum_extension;
    double front_close_gap;
    double front_minimum_contact_overlap;
} AtlasCandidateRefineOptions;

void AtlasCandidateRefineOptions_default(AtlasCandidateRefineOptions *opts);

typedef struct {
    size_t mesh_components;
    size_t candidate_links;
    size_t topology_candidate_links;
    size_t topology_selected_cross_sections;
    size_t topology_ambiguous_cross_sections;
    size_t topology_selected_continuations;
    size_t inactive_cross_component_links;
    size_t inactive_cross_component_continuations;
    size_t unmatched_cross_sections;

    size_t refinement_bundles;
    size_t refinement_mutual_bundles;
    size_t refinement_separate_sheet_bundles;
    size_t refinement_delamination_bundles;
    size_t refinement_sheet_seam_bundles;
    size_t refinement_inconclusive_bundles;
    size_t refinement_pruned_candidate_links;
    size_t refinement_pruned_continuations;
    size_t refinement_remaining_active_cross_sections;
} AtlasCandidateSelectionStats;

typedef struct {
    size_t input_faces;
    size_t intersection_nodes;
    size_t intersection_segments;
    size_t branch_nodes;
    size_t short_strokes_dropped;

    size_t strokes;
    size_t samples;
    size_t continuations;
    size_t cross_sections;
    size_t members;
    size_t matched_source_samples;
    size_t ambiguous_cross_sections;
    size_t truncated_ball_queries;

    size_t support_components;
    size_t largest_component_samples;
    double source_match_fraction;

    /*
     * Winding audit over legacy-admissible relations (pitch > 0 only).
     * Histogram bins on |dwind|: [0,.1) [.1,.35) [.35,.65) [.65,1.5) [1.5,+).
     * "fusing" = |dwind - round(dwind)| <= 0.25 with round != 0, i.e. a
     * relation that would identify two different wraps.  "rejected" counts
     * only what the wind gate actually dropped.
     */
    size_t member_wind_hist[5];
    size_t continuation_wind_hist[5];
    size_t member_wind_fusing;
    size_t continuation_wind_fusing;
    size_t member_wind_rejected;
    size_t continuation_wind_rejected;
    size_t strokes_wind_bridge;     /* EMITTED strokes spanning >= 0.5 turns
                                     * (with wind_cut on this must be ~0) */
    size_t strokes_prior_untrusted; /* u_geo drifts off arclength; no prior */
    double stroke_wind_span_max;    /* turns, over emitted strokes */
    size_t chains_wind_split;       /* slice chains cut at wrap crossings */
    size_t chain_wind_cuts;         /* total cuts applied */
} AtlasCandidateStats;

typedef struct {
    AtlasStripProblem problem;

    AtlasStripSample *samples;
    AtlasStripStroke *strokes;
    AtlasStripMember *members;
    AtlasStripCrossSection *cross_sections;
    AtlasStripContinuation *continuations;
    MonotoneQpAnchor *anchors;

    AtlasCandidateSampleRef *sample_ref;
    int32_t *vertex_mesh_component; /* [input mesh vertex], -1 if unused */
    double *initial_u;          /* [problem.nsamples], strictly feasible */
    double *geo_u;              /* [problem.nsamples] prior targets, NaN on
                                 * untrusted strokes; NULL when pitch <= 0 */
    double *member_prior;       /* immutable geometry admission scores */
    double *continuation_prior;
    uint8_t *member_state;      /* AtlasCandidateState, [nmembers] */
    uint8_t *continuation_state;/* AtlasCandidateState, [ncontinuations] */
    int32_t *member_bundle;     /* refinement bundle, or -1 */
    int32_t *continuation_bundle;

    AtlasCandidateBundle *bundles;
    size_t nbundles;
    AtlasCandidateSelectionStats selection;

    int32_t plane_min;
    int32_t plane_max;
    size_t mesh_components;
    AtlasCandidateStats stats;
} AtlasCandidateSet;

/*
 * face_cube may be NULL (all faces are source cube 0).  vertex_u supplies the
 * local/raw chart only for stroke orientation and initial gauge.  No u/phi
 * mismatch participates in candidate admission.  The returned active problem
 * exactly matches the original HEAD baseline: every admitted assignment and
 * continuation is active.  Refinement happens only after that global solve.
 */
int AtlasCandidates_build(Arena_T arena,
                          const float *vertices,
                          size_t nvertices,
                          const int32_t *triangles,
                          const int32_t *face_cube,
                          size_t ntriangles,
                          const float *vertex_u,
                          const float axis_point[3],
                          const float axis_dir[3],
                          const AtlasCandidateOptions *opts,
                          AtlasCandidateSet *out);

/* Restore the exact original all-admitted HEAD assignment graph. */
int AtlasCandidates_select_baseline(Arena_T arena,
                                    double support_floor,
                                    AtlasCandidateSet *set);

/*
 * Conservative pass: retain at most one candidate per source, and only when
 * source and target already belong to the same retained-triangle component.
 * Cross-component alternatives remain admitted but inactive. "Unmatched" is
 * therefore an explicit and common outcome.
 */
int AtlasCandidates_select_conservative(Arena_T arena,
                                        double support_floor,
                                        AtlasCandidateSet *set);

/*
 * Post-baseline overlap classification. frozen_u is the globally positioned
 * HEAD solution. Coherent cross-component assignments which started on grossly
 * different gauges but were fitted into the same atlas location are bundled.
 * Persistent bundles supported by mutual endpoint continuations are retained
 * as same-sheet cube seams. Bounded fronts are retained as delaminations. A
 * persistent bundle is cut as a missed sheet only when its source sections
 * also contain a same-component alternative. Inconclusive bundles remain
 * untouched.
 */
int AtlasCandidates_select_refined(Arena_T arena,
                                   const double *frozen_u,
                                   const AtlasCandidateRefineOptions *opts,
                                   AtlasCandidateSet *set);

/*
 * Recompute support components and their exact gauge anchors from the CURRENT
 * active link set (member_state / continuation_state).
 *
 * Anything that deactivates links has to call this.  A support component is a
 * connected piece of the candidate graph and carries exactly one exact anchor;
 * cutting links splits components, and a piece that separates from the one
 * holding its anchor is then left with no gauge pinned at all, so the strip
 * system is singular and the solve fails with a reduced-SPD error rather than
 * anything that names the cause.
 *
 * Newly separated components are anchored at their raw per-stroke gauge, which
 * is the same footing the HEAD baseline starts every component on; where that
 * gauge is wrong it is exactly what the radial-order evidence exists to fix.
 */
int AtlasCandidates_reanchor(Arena_T arena, AtlasCandidateSet *set,
                             double support_floor);

/* Small disjoint-mesh seam fixture, including the continuation solve. */
int AtlasCandidates_selftest(void);

#endif
