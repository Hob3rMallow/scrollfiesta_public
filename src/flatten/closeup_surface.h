/* ============================================================================
 * closeup_surface.h -- one continuous sheet surface over an ink region, fitted
 * to every chart that observes it, and its isometric flattening.
 *
 * The assembler's layout is a patchwork: per-cube charts abut with 1-3 vox
 * steps and trim gaps, the ribbon's rotation drifts from chart to chart, and
 * holes cut letters in half.  For a close-up of one ink region (ink_closeup)
 * the charts become OBSERVATIONS of a single surface instead:
 *
 *   observations  every chart pixel of the layout (1 px = 1 vox), rasterised
 *                 in a FRAME (the layout UV rotated so the scroll axis runs
 *                 down the rows, reflected to the team's reading handedness),
 *                 carrying its 3-D point and chart id;
 *   fused grid    nodes every h px; a DATA node (an observation within
 *                 data_reach_px) takes a Gaussian-weighted local quadric of the
 *                 observations within fit_px, reweighted by a Cauchy weight of
 *                 the residuals (robust_rounds), so a chart step becomes one
 *                 smooth surface and a stray other-wrap observation weighs
 *                 nothing; a FILL node (holes and borders, up to fill_px from
 *                 the data) takes a quadric of the data nodes around it, which
 *                 continues the sheet's curvature; then the tear: a node
 *                 with an edge longer than tear_ratio x its frame length, or
 *                 a DATA node floating more than tear_float_vox off every
 *                 observation near it (a fit pulled between two sheets laid
 *                 out side by side), goes back to no surface;
 *   flattening    as-rigid-as-possible on a coarse copy of the grid (TAUCS,
 *                 factor once), carried to the fine nodes bilinearly: the
 *                 flattening of the fused 3-D surface itself, whatever the
 *                 layout's seams did to the frame.
 *
 * Coordinates: world points (z, y, x) float; frame UV in px (1 px = 1 vox of
 * papyrus).  Node (i, j) sits at frame UV (u0 + h i, v0 + h j); arrays are
 * node-major with i fastest: index = j * ni + i.
 * ==========================================================================*/
#ifndef CLOSEUP_SURFACE_H
#define CLOSEUP_SURFACE_H

#include <stddef.h>
#include <stdint.h>

/* ---- observations -------------------------------------------------------- */

typedef struct CloseupObs {
    size_t   n, cap;
    float   *uv;       /* [n*2] frame UV, px */
    float   *p;        /* [n*3] world (z, y, x) */
    float   *w;        /* [n] robust weight in (0, 1] */
    int32_t *chart;    /* [n] source chart id (-1 unknown) */
} CloseupObs;

typedef struct CloseupObsGate {
    double stretch_ratio, stretch_floor;   /* the bake's UV-stretch gate */
    double max_edge3d;                     /* longest 3-D edge of an observed face, vox */
    double min_scale, max_scale;           /* singular values of the complete 3-D -> UV Jacobian; <=0: off */
} CloseupObsGate;

void CloseupObsGate_default(CloseupObsGate *g);
void CloseupObs_free(CloseupObs *o);

/* Rasterise faces into observations at the frame's pixel centres (x + 0.5, y + 0.5): every centre inside
 * a face (edge tolerance 1e-9) takes the barycentric 3-D point.  uv: frame UV per vertex [nv*2] (double).
 * face_use (nullable) selects faces; face_chart (nullable) gives each face's chart id.  Faces failing the
 * gate are counted in *n_gated and skipped.  Appends to o.  0 on success, -1 on no memory. */
int CloseupObs_add_faces(CloseupObs *o, const float *verts, const double *uv, const int32_t *faces, size_t nf,
                         const uint8_t *face_use, const int32_t *face_chart, const CloseupObsGate *gate,
                         size_t *n_used, size_t *n_gated);

/* ---- the fused grid --------------------------------------------------------- */

enum { CLOSEUP_NODE_NONE = 0, CLOSEUP_NODE_DATA = 1, CLOSEUP_NODE_FILL = 2 };

typedef struct CloseupGrid {
    size_t   ni, nj;       /* nodes along u (columns) and v (rows) */
    double   u0, v0, h;    /* node (i, j) at frame UV (u0 + h i, v0 + h j) */
    float   *pos;          /* [nn*3] surface point; NaN where state is NONE */
    float   *tu, *tv;      /* [nn*3] dX/du, dX/dv per px */
    float   *nrm;          /* [nn*3] unit normal: normalize(tu x tv), sign set by CloseupSurface_orient */
    float   *dist;         /* [nn] px to the nearest observation (INFINITY: none in the grid) */
    float   *support;      /* [nn] summed fit weight (DATA: observations; FILL: data nodes) */
    int32_t *chart;        /* [nn] chart of the observation nearest the node, -1 for FILL/NONE */
    int32_t *component;    /* [nn] observed component supporting DATA/FILL, -1 for NONE */
    uint8_t *state;        /* [nn] CLOSEUP_NODE_* */
} CloseupGrid;

/* Nodes covering [umin, umax] x [vmin, vmax] (frame px) every h px.  0, -1 on bad input / no memory. */
int  CloseupGrid_alloc(CloseupGrid *g, double umin, double umax, double vmin, double vmax, double h);
void CloseupGrid_free(CloseupGrid *g);

/* The grid as a file (2026-10-15: a growth's grid before its final flattening, so a later flattening can be tried on a
 * finished growth without growing it again): the frame, and per valid node (state != NONE) its index, state, chart,
 * component, pos and nrm, plus nmeta int32 of the caller's.  tu, tv, dist and support are not kept (a loaded grid has
 * them 0 / INFINITY).  0 on success. */
int  CloseupGrid_save(const char *path, const CloseupGrid *g, const int32_t *meta, size_t nmeta);
/* A CloseupGrid_save file into a fresh grid (CloseupGrid_free it); up to cap meta values into meta and their count into
 * *nmeta.  0 on success, -1 on a bad or truncated file or no memory (g is then empty). */
int  CloseupGrid_load(const char *path, CloseupGrid *g, int32_t *meta, size_t cap, size_t *nmeta);

/* First-order evaluation near node (i, j): p = pos + tu du + tv dv at frame (u, v); 0 if the node is valid. */
int  CloseupGrid_taylor(const CloseupGrid *g, size_t i, size_t j, double u, double v, double p[3]);
/* Bilinear point and normal at frame (u, v) over the four surrounding nodes; 0 when all four are valid. */
int  CloseupGrid_eval(const CloseupGrid *g, double u, double v, double p[3], double n[3]);

typedef struct CloseupFitOpts {
    double fit_px;           /* data fit radius (Gaussian sigma = fit_px / 2) */
    int    min_obs;          /* observations a data fit needs */
    double data_reach_px;    /* nearest observation within this -> DATA node */
    double fill_px;          /* FILL nodes up to this distance from the data */
    double fill_fit_min_px, fill_fit_max_px;
    double robust_vox;       /* Cauchy scale */
    int    robust_rounds;
    double tear_ratio;       /* a node with an edge longer than this x its frame length goes to NONE (<= 0: off) */
    double tear_float_vox;   /* ... and a DATA node this far from every observation within data_reach_px (<= 0: off) */
    /* Gated bridging (2026-10-13): two observed blobs whose charts the wrap index puts on ONE wrap (chart_wrap[chart],
     * INT32_MIN = unknown; a blob's wrap is the vote of its data nodes, decided at bridge_agree) and whose data lie
     * within 2 x bridge_px of each other in the frame share one component before the fill, so the fill fits a joint
     * quadric across the gap and the mesh may connect them; the tear test still cuts a bridge whose 3-D length
     * exceeds tear_ratio x its frame length.  bridge = 0 (the default) keeps every blob its own component. */
    const int32_t *chart_wrap;
    int32_t nchart_wrap;
    int    bridge;
    double bridge_px, bridge_agree;
    int    variational_holes; /* fixed DATA collar, Cartesian thin-plate fill of enclosed holes */
} CloseupFitOpts;

void CloseupFitOpts_default(CloseupFitOpts *o);
/* Wider segmentation-gap completion; keeps the ordinary data fit and all component/geometry gates. */
void CloseupFitOpts_complete(CloseupFitOpts *o);

typedef struct CloseupFitStats {
    size_t nodes, data_nodes, fill_nodes, fill_rejected, data_failed;
    size_t obs, obs_downweighted;          /* weight < 0.5 after the last round */
    double resid_p50, resid_p90, resid_p99; /* |observation - surface|, vox (last round) */
    double seconds;
    size_t torn_stretch, torn_float;       /* nodes the tear returned to NONE (stretched edge / floating data);
                                            * data_nodes / fill_nodes count the fit BEFORE the tear */
    size_t blobs, blobs_with_wrap, bridges; /* observed blobs; those with a decided wrap; same-wrap pairs within
                                            * 2 x bridge_px joined into one component (gated bridging) */
    size_t variational_holes, variational_filled, variational_refused, variational_nodes;
} CloseupFitStats;

/* Fill enclosed holes by minimizing squared second derivatives of all three
 * Cartesian coordinates. DATA and its boundary derivatives are fixed. Exterior
 * gaps and boundaries spanning different components are never joined. Each
 * hole is accepted as a whole after metric/normal checks; provenance is FILL.
 * max_distance_px limits unsupported reach. Stats accumulate in st. */
int CloseupSurface_fill_holes(CloseupGrid *g, double max_distance_px, CloseupFitStats *st);

/* Complete a retained regular grid without refitting or overwriting any valid
 * DATA/FILL node. Importers supply actual face-connectivity component labels,
 * finite positions, state (adopted DATA is 1), chart and saved normals. Reach is
 * recomputed to DATA. A separate geometric orientation reference supports saved
 * inward normals of either sign. Exterior and mixed-component holes are refused.
 * Existing valid arrays stay byte-exact; new nodes alone become FILL. Stats
 * accumulate in st. The caller must preserve its original mesh faces on export. */
int CloseupSurface_complete_fixed(CloseupGrid *g, double max_distance_px, CloseupFitStats *st);

/* As above, with an existing SLIM parameterization per node [nn*2]. Retained
 * UVs are immutable. The same Cartesian thin-plate operator continues UV into
 * holes as a feasibility seed, not as a metric flattening substitute. The
 * unchanged metric guard uses the full actual 3-D/UV Jacobian and requires
 * positive UV area. Accepted new nodes receive UVs; refused nodes retain their
 * input UVs. NULL selects the regular frame. For a cell with one new corner,
 * export only the complement of its retained three-corner face, respecting the
 * original diagonal. Other completed cells use diagonal (0,0)--(1,1).
 * Run final SLIM on the enlarged mesh.
 * A node cannot represent multiple UV copies at a cut; import before cutting. */
int CloseupSurface_complete_fixed_uv(CloseupGrid *g, double *uv_nodes, double max_distance_px, CloseupFitStats *st);

typedef struct CloseupProposalStats {
    size_t holes, proposed, refused, contact_refused, nodes;
    int metric_qualification_pending; /* always 1: a proposal is not an accepted completion */
} CloseupProposalStats;

/* Detached working grid only: stage finite, oriented Cartesian patches while
 * postponing their metric decision to joint SLIM with geometry-preserving cuts.
 * Retained nodes are immutable. Reach/components/world orientation and the
 * existing nonlocal node-contact guard remain enforced, rejecting whole new
 * patches. Use the original regular frame for the enlarged mesh's injective
 * SLIM seed. *patch_ids is malloc'd [nn], -1 retained/absent, otherwise a dense
 * proposal ID. Free it with free(). No proposal is qualified by this API. */
int CloseupSurface_propose_fixed(CloseupGrid *g, double max_distance_px, int32_t **patch_ids,
                                 CloseupProposalStats *st);

typedef struct CloseupProposalMetricStats {
    size_t triangles, patches, passed_patches, refused_patches;
    size_t metric_failures, orientation_failures, invalid_faces, empty_patches;
    double min_squared_scale, max_squared_scale;
    int metric_qualified; /* metric only; source/contact/topology/fiber audits are separate */
} CloseupProposalMetricStats;

/* Final actual-UV metric gate, after joint SLIM and optional cuts. face_patch
 * [nf] is -1 for retained faces, otherwise the staged patch ID; cut duplicates
 * retain their patch identity. Every appended triangle must satisfy the same
 * [0.25,4] squared singular-value bounds and the original frame's UV area sign. A failure marks
 * the entire patch in patch_pass [npatches], without modifying any geometry.
 * Empty patch IDs cannot qualify. 0 means audit ran, not that it passed. */
int CloseupSurface_audit_proposal_metric(const float *verts, size_t nv, const double *uv, const double *frame_uv,
                                        const int32_t *faces, size_t nf, const int32_t *face_patch,
                                        size_t npatches, uint8_t *patch_pass, CloseupProposalMetricStats *st);

typedef struct CloseupRebuildResult {
    size_t nodes, faces, patches;
    int32_t *node_patch;       /* [nn], -1 for observed/absent, otherwise inferred patch */
    int32_t *component_group;  /* [nn], supported grouping; original DATA component is unchanged */
    int32_t *face_nodes;       /* [faces*3], inferred triangles in grid-node indices */
    int32_t *face_patch;       /* [faces], carries through subsequent vertex-only cuts */
} CloseupRebuildResult;

typedef struct CloseupRebuildWorldStats {
    size_t pairs, observed_observed, observed_inferred, inferred_inferred;
    size_t overlaps, stabs, folds;
    int complete;
} CloseupRebuildWorldStats;

typedef struct CloseupRebuildStats {
    size_t data_nodes, old_fill_nodes, domain_nodes, enclosed_gap_nodes;
    size_t source_components, supported_groups, certified_joins, refused_joins;
    size_t patches, proposed_patches, refused_patches, proposed_nodes, proposed_faces;
    size_t reach_refused, orientation_refused, contact_refused, unsupported_nodes;
    size_t separated_cells;
    size_t alternate_diagonals;
    size_t screened_nodes;
    double tangent_screen_energy;
    size_t old_fill_start_nodes, target_start_patches, motion_sweeps, motion_vertices;
    size_t coupled_passes, coupled_blocks, coupled_vertices, coupled_solves;
    size_t orientation_before_motion, orientation_after_motion;
    size_t transport_components, transport_anchored_components, transport_unanchored_components;
    size_t transport_unanchored_faces, transport_degenerate_faces, transport_bad_edges;
    size_t transport_winding_edges, transport_nonmanifold_edges;
    CloseupRebuildWorldStats world_exact, world_proximity;
    double target_error_before_motion, target_error_after_motion;
    int metric_qualification_pending; /* always 1; source/contact/topology/fiber proofs also required */
} CloseupRebuildStats;

typedef struct CloseupRebuildDomainStats {
    size_t required_nodes, required_enclosed_gap_nodes, added_nodes, nodes, required_beyond_reach, data_donors;
    double reach_px, source_frame_bounds[4];
} CloseupRebuildDomainStats;

/* Explicit supported-regrowth-domain-v2: in-place union of the prior required
 * mask and incident-DATA frame reach, bounded by the import's one-cell outer
 * padding. This is only a proposal mask. It cannot merge source groups or
 * bypass physical DATA witnesses, local orientation, contact or metric gates. */
int CloseupSurface_rebuild_regrowth_domain(const CloseupGrid *g, uint8_t *domain,
                                           double max_distance_px, CloseupRebuildDomainStats *st);

/* Separate DATA-preserving repair of prior inferred material. Input is a
 * detached original regular grid with adopted DATA mapped to state 1 and
 * actual source face-connectivity labels. domain[nn] is exactly the previous
 * valid footprint plus its enclosed gaps, excluding the exterior apron.
 * observed_faces[3*observed_nf] are ALL fully observed source faces, in grid
 * node indices. Their diagonals and winding constrain the returned topology.
 * DATA arrays remain byte-exact, including saved normals and component IDs.
 * Only FILL is discarded/rebuilt. Natural-boundary Cartesian thin plates are
 * anchored exclusively to DATA; a continuous DATA-tangent screen fixes otherwise
 * free affine modes with lambda=(h*h/(h*h+frame_distance_squared))^2. The same
 * supported component supplies every screening reference. No unknown-count
 * heuristic rejects a solve.
 * Different source components connect only after same-chart boundary witnesses
 * all pass physical tangent, normal-offset, orientation and length tests.
 * Unsupported groups stay separated. New triangles, including all-DATA caps,
 * are returned explicitly so exporter and geometric gates share one topology.
 * Original observed vertices/faces must remain mapped immutable prefixes in
 * the exported mesh; preserve the original full artifact and replacement
 * counts. Final joint SLIM/cuts and unchanged proposal metric audit are
 * mandatory. Successful return stages a proposal, never qualifies it.
 * On runtime error g is unchanged. On success old FILL may be absent where a
 * whole inferred patch was refused. Result arrays use malloc; free below. */
int CloseupSurface_rebuild_generated(CloseupGrid *g, const uint8_t *domain,
                                     const int32_t *observed_faces, size_t observed_nf,
                                     double max_distance_px, CloseupRebuildResult *result,
                                     CloseupRebuildStats *st);
void CloseupRebuildResult_free(CloseupRebuildResult *result);

/* Fit the grid to the observations (see the file comment).  obs->w is overwritten by the robust rounds.
 * Normals are normalize(tu x tv) (unoriented).  0 on success, -1 on no memory / no observations. */
int CloseupSurface_fit(CloseupGrid *g, CloseupObs *obs, const CloseupFitOpts *o, CloseupFitStats *st);

/* CloseupSurface_fit, recording its robust data rounds in cache, so that the next fit of the same grid with the same
 * options -- after observations were only APPENDED to obs (a completion's adopt rounds) -- refits only the nodes the
 * new observations can reach through the rounds and takes every other node's rounds from the record.  The result is
 * CloseupSurface_fit's, bit for bit; a cache that does not match the grid, options or observations is rebuilt. */
typedef struct CloseupFitCache CloseupFitCache;
CloseupFitCache *CloseupFitCache_new(void);
void CloseupFitCache_free(CloseupFitCache **cache);
int CloseupSurface_fit_cached(CloseupGrid *g, CloseupObs *obs, const CloseupFitOpts *o, CloseupFitStats *st,
                              CloseupFitCache *cache);

/* Orientation: inward[k*3..] is the unit direction toward the scroll axis at node k (any array of nn
 * entries; zero vectors abstain).  Every normal is set to s * normalize(tu x tv) with s the sign of the
 * summed vote sum (tu x tv) . inward over valid nodes, so the normals point inward.  Returns s (+1: the
 * frame shows the recto in the team's reading handedness -- (dP/dcol x dP/drow) . n_inward > 0 -- and
 * -1: mirror-imaged), 0 when nothing votes.  *agree (nullable): fraction of |vote| agreeing with s. */
int CloseupSurface_orient(CloseupGrid *g, const float *inward, double *agree);

typedef struct CloseupFlatStats {
    size_t coarse_nodes, coarse_tris, components, iterations;
    size_t coarse_fringe;                      /* SLIM: empty lattice nodes continued from the sheet (CLOSEUP_SLIM_FRINGE_*) */
    double energy_first, energy_last;
    double energy_hist[64];                    /* the energy of the first 64 iterations (report diagnostics) */
    double logs_p50_before, logs_p99_before;  /* |log sigma| of the frame over coarse triangles */
    double logs_p50_after, logs_p99_after;    /* ... of the flattening */
    double area_ratio;                         /* UV area / 3-D area after */
    double seconds;
} CloseupFlatStats;

/* SLIM isometric flattening (see the file comment): uv_out [nn*2] receives the flattened UV of every valid node
 * (NaN elsewhere, and for nodes no coarse triangle reaches), rigidly aligned (area-weighted Procrustes) to
 * the frame so the orientation is kept. coarse_px is rounded to a multiple of h. tol is retained for source
 * compatibility and ignored. No ARAP or frame fallback. 0 on success. */
int CloseupSurface_flatten(const CloseupGrid *g, double coarse_px, int iters, double tol, double *uv_out,
                           CloseupFlatStats *st);

/* Injective flattening (2026-10-15; the user, on the grown candidates' filaments: "we're doing SLIM for the flatten,
 * right?"): SLIM (src/flatten/slim_refine.c: symmetric Dirichlet, flip-free line search) on the coarse grid, started
 * from the frame itself (node (i, j) at frame UV (u0 + h i, v0 + h j): injective by construction), with SCAF scaffold
 * triangles over every coarse cell the surface does not cover (holes, the gaps between components, the rest of the
 * grid's rectangle), so no part can slide over another; every fine node is then carried by barycentric interpolation
 * in the coarse triangle -- surface or scaffold -- that holds it in the frame: a piecewise-linear map of flip-free
 * triangles, injective. When SLIM fails or leaves a flip this returns an error without a replacement. uv_out as
 * CloseupSurface_flatten's, rigidly aligned to the frame; 0 on success. */
int CloseupSurface_flatten_slim(const CloseupGrid *g, double coarse_px, int iters, double *uv_out, CloseupFlatStats *st);

typedef struct CloseupDedupStats {
    size_t pairs;          /* node pairs within r_vox in 3-D, normals aligned, more than frame_px apart in the frame */
    size_t data_pairs;     /* ... of them DATA with DATA (counted, never dropped) */
    size_t dropped_fill;   /* FILL nodes dropped */
    double seconds;
} CloseupDedupStats;

/* Duplicated papyrus (2026-10-15): a FILL node within r_vox in 3-D of another valid node that lies more than frame_px
 * away in the frame (normals within acos(cos_min)) extrapolates over papyrus the surface already holds elsewhere -- the
 * same papyrus twice in the flattened sheet.  It goes (state NONE, pos NaN, chart and component -1); of two such FILL
 * nodes the one farther from the observations (dist; the higher index on a tie).  DATA is never dropped.  0 on
 * success, -1 on no memory. */
int CloseupSurface_dedup_fill(CloseupGrid *g, double r_vox, double frame_px, double cos_min, CloseupDedupStats *st);

typedef struct CloseupApronStats {
    size_t fill_before;    /* FILL nodes before the trim */
    size_t dropped;        /* FILL nodes with observed papyrus on one side only: the outward apron, dropped */
    double seconds;
} CloseupApronStats;

/* The outward apron (2026-10-15): a FILL node stays only where it lies BETWEEN observed papyrus -- DATA within reach_px
 * on both sides along at least one of the four grid lines through it (row, column, both diagonals): a gap or a hole the
 * sheet surrounds.  FILL with DATA on one side only extends the sheet past its observed edge, where the extrapolation
 * cuts across the neighbouring wraps; it goes (state NONE, pos NaN, chart and component -1).  DATA is never dropped.
 * 0 on success, -1 on no memory or a reach under one grid step. */
int CloseupSurface_trim_apron(CloseupGrid *g, double reach_px, CloseupApronStats *st);

/* Triangulate the valid cells whose four (or three) corners have finite uv_nodes: every such node becomes a
 * vertex (point, uv, normal, node index).  Arrays are malloc'd; free with free().  0 on success. */
int CloseupSurface_mesh(const CloseupGrid *g, const double *uv_nodes, float **verts, double **uv, float **nrm,
                        int32_t **faces, int32_t **node_of_vert, size_t *nv, size_t *nf);

/* ---- the frame -------------------------------------------------------------- */

/* A layout UV (u, v) maps to the frame by an optional reflection (u -> -u), then a rotation by theta about
 * (cu, cv) (the reflected coordinates): the scroll axis runs down the rows (+v = increasing axial
 * coordinate) and the recto shows in the team's reading handedness. */
typedef struct CloseupFrame {
    int    reflect;
    double theta, cu, cv;
    double hand_vote;        /* area-weighted handedness vote in the layout (+ team, - mirror-imaged) */
    double hand_agree;       /* fraction of |vote| agreeing with its sign */
    double axial_gx, axial_gy;   /* mean gradient of the axial coordinate (after reflection, before rotation) */
} CloseupFrame;

/* The frame of the faces selected by face_use (nullable = all) of a UV mesh: axial[v] is the axis arc length
 * at vertex v, inward[v*3..] the unit direction toward the axis.  (cu, cv) is the rotation centre in the
 * layout's own coordinates (it is reflected with u when reflect is set).  0 on success, -1 when no face votes. */
int  CloseupSurface_frame(const float *verts, const double *uv, const int32_t *faces, size_t nf,
                          const uint8_t *face_use, const double *axial, const float *inward, double cu, double cv,
                          CloseupFrame *fr);
void CloseupFrame_apply(const CloseupFrame *fr, double u, double v, double *fu, double *fv);

/* Rasterise a per-vertex uint8 attribute of a UV mesh on a W x H raster whose pixel (x, y) is centred at
 * (umin + x + 0.5, vmin + y + 0.5): every covered pixel takes the attribute of the corner with the largest
 * barycentric weight; uncovered pixels are left untouched.  0 on success. */
int  CloseupSurface_raster_attr(const double *uv, const int32_t *faces, size_t nf, const uint8_t *attr,
                                double umin, double vmin, size_t W, size_t H, uint8_t *out);

/* 2x2 closest rotation angle of M = [m00 m01; m10 m11] (polar decomposition). */
double CloseupSurface_rot_angle(double m00, double m01, double m10, double m11);

int CloseupSurface_selftest(void);
int CloseupSurface_rebuild_selftest(void);

#endif /* CLOSEUP_SURFACE_H */
