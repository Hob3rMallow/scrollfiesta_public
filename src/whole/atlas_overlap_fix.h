#ifndef ATLAS_OVERLAP_FIX_INCLUDED
#define ATLAS_OVERLAP_FIX_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "../unroll/piece_set.h"
#include "../unroll/scaffold.h"

/* ============================================================================
 * atlas_overlap_fix.h -- move welded sheet fragments by whole wraps, and park
 * the ones that cannot be placed.
 *
 * TWO LEVELS OF STRUCTURE.
 *
 *   CHART   one connected component of the retained mesh.  PieceSet never
 *           shares vertices across cubes, so a chart is exactly one cube's
 *           sheet.
 *
 *   GROUP   one connected component of the AUTHENTIC-NEIGHBOUR graph over
 *           charts.  Two charts are authentic neighbours when they sit within
 *           a weld distance of each other in XYZ with agreeing normals -- that
 *           is, when the XYZ seam weld would join them into one surface.  A
 *           group is therefore one physically welded sheet fragment, spanning
 *           however many cubes it spans.
 *
 * THE GROUP IS THE ONLY THING THAT MOVES.  Charts inside a group are held
 * together by a real weld; shifting one of them alone would tear it.  So the
 * unit of rigid motion is the whole group, and the output is a small atlas of
 * welded fragments rather than a confetti of independently-placed charts.
 *
 * The premise is deliberately narrow.  The incoming registered field already
 * places most groups correctly; where it fails, the failure is a WHOLE TURN of
 * arc length, because that is the only quantity an independently-gauged cube
 * can get wrong.  So the only move available is
 *
 *     u[group] += k * (2*pi*r_group),   k integer,
 *
 * chosen so the group stops intersecting anything else in (u, v).  No shape is
 * changed, no group is deformed, and no two charts are ever merged.
 *
 * Four rules, in order:
 *
 *   1. Kill triangles with hopeless UV/XYZ edge stretch.  Comparing each UV
 *      edge to its source edge catches parameterization damage without
 *      deleting a triangle merely because it was already skinny in XYZ.
 *   2. Group the charts by weldability.  Weldable neighbours belong together
 *      in u because that is how scroll geometry works.
 *   3. A group with no overlap against any OTHER group is SETTLED -- it is
 *      already right, so it never moves.
 *   4. Every other group is tried at +/-1, +/-2, ... wraps, furthest-from-axis
 *      first, and the first placement that clears its overlaps wins.  A group
 *      that no shift can place is PARKED whole into its own empty v band below
 *      the atlas, which is overlap-free by construction and makes the residue
 *      visible instead of silently interleaved.
 *
 * Overlap BETWEEN two charts of the same group is reported separately as
 * intra_overlap: no rigid move can fix it, because they move together.  A high
 * intra_overlap count means the group itself is mis-registered internally and
 * is the signal to cut it, not to shift it.
 * ==========================================================================*/

typedef struct {
    double aspect_ratio_kill;    /* UV/XYZ edge stretch above which face dies (8) */
    size_t floater_max_faces;     /* kill whole charts at or below this size
                                   * before fitting (0 disables) */
    size_t floater_group_max_faces; /* after authentic welds, kill whole groups
                                     * at or below this total face count before
                                     * fitting (0 disables) */
    double chart_qc_kill;         /* kill a whole chart when any retained face's
                                   * Sander sigma1/sigma2 exceeds this
                                   * (0 disables) */
    /* 6 vox: what the XYZ seam weld itself bridges (2*rho_max).  Must stay
     * below the 7-vox inter-wrap clearance or a "neighbour" may be the next
     * wrap in, which would weld two wraps into one immovable group. */
    double neighbour_distance;
    double neighbour_normal_dot; /* |dot| of the two CONTACT normals (0.7) */
    size_t neighbour_min_shared; /* close vertex pairs required to accept (3) */
    int    max_shift_wraps;      /* search +/-1 .. +/-this many wraps (4) */
    double overlap_cell_size;    /* broad-phase cell floor, u/v units (4) */
    double park_margin;          /* blank v gap around a parked group (16) */
    int    park_unplaceable;     /* 1 = park what no shift can place (1) */

    /* Per-chart wind reconciliation inside each welded group. */
    int    wind_correct;         /* 1 = integrate weld wind errors (1) */
    int    wind_only;            /* stop after one global reconciliation pass;
                                  * used to seed bounded tiled optimization */
    double wind_residual_limit;  /* turns; above this a weld is not gauge
                                  * evidence and is excluded from the tree
                                  * (0.25) */
    int    wind_cap;             /* reject |wraps| beyond this as noise (8) */
    /* Experimental ablation: derive each weld's required relative depth from
     * the registered phase sampled at the fixed XYZ contact, rather than from
     * the current u placement.  The phase relation is source-stable when u is
     * changed by the optimizer. */
    int    phase_targets;        /* 1 = source-phase weld targets (0 while
                                  * the controlled ablation is evaluated) */
    /* Controlled ablation: for a geometrically graded group, make each weld's
     * relative-depth target agree with the same weld-chain winding potential
     * used by the graded move.  The unwind then remains an intact topology
     * rather than being immediately classified as a collection of tears. */
    int    graded_targets;
    int    phase_graded;         /* derive the graded correction from lifted
                                  * source phase minus the winding already in u,
                                  * instead of absolute-radius differences */
    int    phase_all_groups;     /* allow that correction in every connected,
                                  * cycle-consistent group, not only one already
                                  * diagnosed with wrong-wind overlap */
    int    phase_flat_lock;      /* if a fully connected, cycle-consistent
                                  * source-phase lift assigns one relative bin
                                  * to a group, forbid chart/neighbourhood moves
                                  * that would split it; rigid moves remain */
    int    phase_lateral_targets;/* set trusted lateral dk targets from the
                                  * difference of lifted-source and incoming-u
                                  * winds, rather than assuming every pair is
                                  * already in the same gauge */
    int    phase_regions;        /* turn strong, consistent cross-group phase
                                   * links into one compound gauge region;
                                   * relative group gauges stay fixed while the
                                   * region searches a common winding */
    int    lock_graded;          /* after a graded move, preserve that group's
                                  * relative winding; rigid translations and
                                  * another graded assignment remain legal */
    int    postgraded_refine;    /* protect a graded group until its first
                                  * coherent assignment, then permit exact
                                  * chart/neighbourhood cleanup (experimental) */

    /* Tabu re-gauge.  When set, the rigid group placement is replaced by a
     * tabu search that assigns every CHART its own integer winding depth
     * k in [-tabu_span, +tabu_span] (0 = where the field put it).  The energy
     * is  lambda_ov * (shared occupancy cells between charts)
     *   + lambda_nb * sum over weld edges of contacts-weighted |dk - target|
     *   + lambda_rad * huber(wind implied by u  -  wind implied by radius)
     *   + lambda_gauge * 0.5 * (k - phase_bin)^2.
     * Welds may TEAR where the accumulated drift demands it -- that is the
     * point: a rigid group move cannot fix a spiral whose u under-advances,
     * a per-chart re-gauge can.  Nothing is parked in tabu mode. */
    int    tabu;                 /* 1 = per-chart tabu search (0) */
    int    tabu_span;            /* depth range +/- this many wraps (4) */
    int    tabu_tenure;          /* iterations a departed depth stays banned (40) */
    int    tabu_max_iters;       /* hard iteration cap (2000) */
    int    tabu_stall;           /* stop after this many non-improving (300) */
    double tabu_cell;            /* collision-proxy cell size, vox (4) */
    int    tabu_ladder_override; /* reuse one global field-ladder calibration
                                  * across independently solved atlas tiles */
    double tabu_ladder_step;     /* signed du / (2*pi*r) per outward wind */
    double tabu_ladder_u0;       /* global u intercept for ladder inversion */
    int    tabu_active_core;     /* block-coordinate solve: only charts whose
                                  * cube origin lies in this half-open box move */
    long   tabu_core_z_lo, tabu_core_z_hi;
    long   tabu_core_y_lo, tabu_core_y_hi;
    long   tabu_core_x_lo, tabu_core_x_hi;
    double lambda_ov;            /* energy per shared occupancy cell (1) */
    double lambda_nb;            /* energy per contact-weighted depth tear (4) */
    double lambda_rad;           /* energy per turn of radial deviation (2) */
    double lambda_gauge;         /* keep a graded group's free common gauge
                                  * near the incoming atlas without penalizing
                                  * its relative phase bins (0 = disabled) */
    double nb_unhappy_w;         /* steering-only downweight of weld edges to a
                                  * currently-colliding neighbour (0.25);
                                  * 1.0 disables the happiness steering */
    /* Lateral coupling: same-wind neighbours the weld graph cannot see.  Two
     * charts at the SAME radius (|dr| < 0.6 pitch), in adjacent cubes, and
     * side by side in u are one sheet interrupted by a hole or a trimmed
     * seam -- no weld exists, so without this term the search tears them
     * apart for free (measured: 60% of 620 such pairs torn).  Same radius
     * means same wind regardless of any shared gauge error, so coupling
     * their depths (target dk = 0) is safe evidence. */
    double lateral_w;            /* contact-equivalent weight (8; 0 = off;
                                  * 16 trades overlap for more coherence) */
    double lateral_u_gap;        /* max |du| to call charts side-by-side (200) */
    /* A HAPPY FAMILY is a welded group with zero internal overlap: it is
     * already assembled correctly and only its placement is wrong.  Edges
     * between two members of a happy family are boosted by this factor (and
     * exempted from the unhappy-neighbour steering discount -- a layered
     * family is collectively unhappy exactly when it must move TOGETHER), so
     * per-chart defection becomes expensive while the whole-group move,
     * which pays no internal cost, stays free.  (<= 1 disables.) */
    double family_w;             /* boost on intra-family edges (8) */

    /* Stage C of the alternation: after tabu fixes the WINDINGS, re-register
     * the continuous per-chart u shifts with those windings held fixed --
     * intact welds become closure springs, same-depth laterals keep-together
     * springs, and radial pairs enforce the measured per-wind separation so
     * the continuous stage cannot re-collapse what tabu just separated.
     * Tabu owns the integers, the register owns the metric. */
    int    relayout;             /* 1 = re-register after tabu (0) */
    double relayout_seam_cap;    /* max deviation from the metric solution
                                  * across one intact topology edge (25;
                                  * 0 disables the topology-first budget) */
} AtlasOverlapFixOptions;

void AtlasOverlapFixOptions_default(AtlasOverlapFixOptions *opts);

typedef enum {
    ATLAS_OVERLAP_FIX_SETTLED = 0,  /* no cross-group overlap; never moved */
    ATLAS_OVERLAP_FIX_SHIFTED = 1,  /* moved a whole number of wraps in u */
    ATLAS_OVERLAP_FIX_PARKED  = 2,  /* no shift worked; moved out of the atlas */
    ATLAS_OVERLAP_FIX_STUCK   = 3   /* still overlapping; parking disabled */
} AtlasOverlapFixAction;

typedef struct {
    int32_t group;
    size_t  ncharts;
    size_t  nfaces;
    size_t  nvertices;
    double  mean_radius;      /* distance from the scroll axis */
    double  radius_spread;    /* max-min: how rigid a constant shift really is */
    double  circumference;    /* 2*pi*mean_radius: one wrap of u here */
    double  shift_u;
    double  shift_v;
    int32_t wraps;            /* shift_u / circumference, signed */
    size_t  overlap_raw;      /* cross-group pairs as the field arrived */
    size_t  overlap_before;   /* cross-group face pairs, post-wind */
    size_t  overlap_after;
    size_t  intra_overlap;    /* chart-vs-chart inside this group: unfixable */
    size_t  intra_wrongwind;  /* of those, pairs >= 0.75 pitch apart RADIALLY:
                               * genuinely mis-wound internals.  The rest are
                               * folds -- doubled geometry, correctly
                               * parameterized -- and must not cost a group
                               * its happy-family protection. */
    int     action;           /* AtlasOverlapFixAction */
} AtlasOverlapFixGroup;

/* Per-chart record of the tabu decision (tabu mode only). */
typedef struct {
    int32_t chart;
    int32_t group;
    int32_t neighbourhood;    /* lateral same-sheet component */
    int32_t phase_flat;       /* group's source-phase lift has one bin */
    int32_t phase_region;     /* compound relative-gauge region, -1 if none */
    int32_t region_offset;    /* fixed group base inside that region */
    int32_t cube;            /* index into ps->ids */
    size_t  nfaces;
    size_t  nvertices;
    double  qc_max;          /* worst initial Sander sigma1/sigma2 */
    double  radius_med;      /* median vertex distance from the axis */
    double  u_med;           /* median post-wind u */
    double  w_phi0;          /* wind implied by u_med, continuous turns */
    double  w_rad;           /* wind implied by radius, recentred */
    double  w_phase;         /* source angle after weld-graph integer gauge */
    int32_t phase_gauge;     /* integer added to sense*phi/(2*pi) */
    int32_t phase_bin;       /* diagnostic correction: w_phase - w_phi0,
                              * median-centred within the welded group */
    double  wind;            /* absolute continuous wind, turns: the source
                              * phase lift anchored per gauge island to the
                              * final field (w_phi0 + k) by rounded median.
                              * THIS is the winding a consumer should trust;
                              * k is only the round-local correction. */
    int32_t k;               /* assigned winding depth (0 = unmoved) */
    int32_t gbin;            /* weld-chain wind bin within the group */
    double  shift_u;         /* the u shift that depth produced */
    size_t  relayout_rank;   /* fixed left-to-right order after tabu */
    double  relayout_metric_shift_u; /* unconstrained metric preference */
    double  relayout_shift_u;/* accepted continuous post-tabu U slide */
    int64_t cells_initial;   /* shared occupancy cells at k = 0 */
    int64_t cells_final;     /* shared occupancy cells at the final depth */
} AtlasOverlapFixChart;

/* One accepted tabu move.  Deltas are canonical energy in shared-cell units,
 * so the trace shows which term paid for each move. */
typedef struct {
    int32_t iter;
    int32_t kind;            /* 0 chart, 1 welded group, 2 lateral
                              * neighbourhood, 3 graded group, 4 compound
                              * relative-gauge region */
    int32_t id;              /* id in the move class above */
    int32_t from;            /* chart: departed depth; compound: always 0 */
    int32_t to;              /* chart: new depth; compound: applied delta */
    double  d_overlap;
    double  d_neighbour;
    double  d_radial;
    double  energy_after;
} AtlasOverlapFixTabuMove;

/* One archived assignment selected for exact SAT verification. */
typedef struct {
    size_t state;              /* accepted-search state; 0 is the input */
    size_t exact_cross;
    size_t exact_intra;
    size_t layout_abs;         /* sum |k|, diagnostic tie-break */
    size_t weld_jumps_25;      /* usable welds with |final contact du| > 25 */
    size_t weld_jumps_100;     /* usable welds with |final contact du| > 100 */
    double search_overlap;     /* exact face-pair count tracked in search */
    double neighbour;
    double radial;
    double energy;
    double exact_score;
    double weld_mean_abs;      /* mean |du + shift[a] - shift[b]|, vox */
    double weld_rms;           /* RMS of that physical weld mismatch, vox */
    double weld_max;           /* worst physical weld mismatch, vox */
    int    winner;
} AtlasOverlapFixTabuCandidate;

/* One usable topology edge in the final tabu assignment. */
typedef struct {
    int32_t edge;
    int32_t chart0;
    int32_t chart1;
    int32_t target;
    int32_t final_delta;
    int32_t lateral;
    int32_t family;
    double  weight;
    double  du;               /* contact mismatch before tabu */
    double  final_du;         /* contact mismatch after the winning shifts */
    double  phase_turn;       /* (w_phase-w_phi0)[0] - same[1] */
    double  phase_residual;   /* distance from phase_turn to an integer */
} AtlasOverlapFixTabuEdge;

/* One chart pair with at least one exact SAT face-pair intersection in the
 * winning discrete assignment.  This is deliberately a diagnostic rather
 * than an energy term: the continuous pass decides whether the metric can
 * remove it without changing the fixed winding topology. */
typedef struct {
    int32_t chart0;
    int32_t chart1;
    uint32_t face_pairs;
    int32_t edge;             /* direct topology edge, or -1 */
    int32_t target;           /* direct edge target; 0 when edge == -1 */
    int32_t final_delta;
    int32_t lateral;          /* 0 weld, 1 lateral, -1 no direct edge */
    double radius_delta;
    double min_xyz;
} AtlasOverlapFixTabuOverlap;

typedef struct {
    int32_t chart_lo;
    int32_t chart_hi;
    int32_t face_lo;
    int32_t face_hi;
    double lower;
    double final_delta;
    double slack;
    uint32_t updates;
    int topology;
    int collision;
    int initial;
    int active;
} AtlasOverlapFixRelayoutBound;

typedef struct {
    int32_t chart0;
    int32_t chart1;
    int32_t face0;
    int32_t face1;
    int32_t allowed;
    int32_t reason; /* 0 hard, 1 direct, 2 XYZ-near, 3 self, 4 topology */
    double min_xyz; /* closest source-space vertex pair on these triangles */
    double radius_delta;
} AtlasOverlapFixRelayoutResidual;


/* Filled by the solve in tabu mode; arrays are arena-allocated.  Optional. */
typedef struct {
    AtlasOverlapFixChart *chart;
    size_t nchart;
    AtlasOverlapFixTabuMove *move;
    size_t nmove;
    const int8_t *state_k;          /* [nstate * nchart], accepted search
                                     * states; state 0 is the input */
    size_t nstate;
    AtlasOverlapFixTabuCandidate *candidate;
    size_t ncandidate;
    AtlasOverlapFixTabuEdge *edge;
    size_t nedge;
    AtlasOverlapFixTabuOverlap *overlap;
    size_t noverlap;
    AtlasOverlapFixRelayoutBound *relayout_bound;
    size_t nrelayout_bound;
    AtlasOverlapFixRelayoutResidual *relayout_residual;
    size_t nrelayout_residual;
    const int32_t *vertex_chart;   /* [ps->nv] chart per vertex, -1 unused */
    const double *u_post_tabu;     /* [ps->nv] the winner placement BEFORE
                                    * the relayout touched it: the pure tabu
                                    * decision, for per-layer inspection */
} AtlasOverlapFixTabuReport;

typedef struct {
    size_t charts_before_cull;
    size_t total_charts;
    size_t total_groups;
    size_t total_faces;
    size_t faces_killed;
    size_t groups_before_cull;
    size_t groups_killed_small;
    size_t faces_killed_groups;
    size_t group_size_hist[12];
    size_t faces_killed_stretch;
    size_t faces_killed_charts;
    size_t charts_killed_small;
    size_t charts_killed_qc;
    /* Initial chart populations after the face stretch cull, before whole-chart
     * culling.  Size bins: 1, 2-3, 4-7, 8-15, 16-31, 32-63, 64-127,
     * 128-255, 256-511, 512+.  QC bins: <=2, <=4, <=8, <=16, <=32, <=64,
     * <=128, <=256, >256. */
    size_t chart_size_hist[10];
    size_t chart_qc_hist[9];
    size_t authentic_neighbour_edges;
    size_t largest_group_charts;

    /* Weld gauge evidence.  weld_residual_hist bins |du/circumference - round|
     * at 0.05 / 0.15 / 0.25 / 0.35 / above: a pile-up in bin 0 is what makes
     * "the error is a whole number of turns" a measurement rather than an
     * assumption. */
    size_t weld_residual_hist[5];
    size_t weld_edges_nonzero_wrap;
    /* Same diagnostic for the source phase target
     * -sense*(phi[chart0]-phi[chart1])/(2*pi). */
    size_t weld_phase_residual_hist[5];
    size_t weld_phase_edges_nonzero;
    size_t weld_edges_trusted;
    size_t weld_edges_untrusted;
    size_t weld_edges_oversized;
    size_t gauge_tree_edges;
    size_t gauge_cycle_edges;
    size_t gauge_cycle_consistent;
    size_t gauge_cycle_inconsistent;   /* holonomy: unfixable by any gauge */
    size_t gauge_islands;
    size_t charts_wind_corrected;
    int32_t wind_wraps_min;
    int32_t wind_wraps_max;

    /* Overlap at each of the three stages: as it arrived, after the per-chart
     * wind reconciliation, and after the rigid group placement. */
    size_t overlap_pairs_raw;
    size_t intra_overlap_pairs_raw;
    size_t overlap_pairs_before;      /* cross-group, post-wind */
    size_t overlap_pairs_after;
    size_t intra_overlap_pairs;       /* inside one group: no shift can fix */
    /* |dr|/pitch histogram of the post-wind intra-group pairs, bins at
     * <0.25 / <0.75 / <1.25 / <1.75 / <2.5 / more.  Bins 0-1 are one surface
     * folded onto itself, which no re-gauge can fix; bin 2 is the NEXT WRAP
     * sitting on this one's u -- exactly what a per-chart wind reassignment
     * can fix.  This is the ceiling measurement for any gauge-based repair. */
    size_t intra_radial_hist[6];
    size_t groups_overlapping_before;
    size_t groups_overlapping_after;

    size_t settled_groups;
    size_t shifted_groups;
    size_t parked_groups;
    size_t stuck_groups;

    double park_v_base;       /* atlas v_min; parked groups live below this */

    /* Tabu re-gauge (tabu mode only). */
    size_t tabu_iterations;
    size_t tabu_moves_chart;
    size_t tabu_moves_group;
    size_t tabu_moves_nbhd;          /* whole weld+lateral components moved */
    size_t tabu_moves_graded;        /* graded rewinds applied */
    size_t tabu_moves_region;        /* compound relative-gauge regions */
    size_t tabu_graded_groups;       /* groups eligible for graded rewinds */
    size_t tabu_charts_moved;        /* final depth != 0 */
    size_t tabu_active_charts;       /* movable charts in block-coordinate mode */
    size_t tabu_welds_torn;          /* usable weld edges left off-target */
    size_t tabu_edges_excluded;      /* untrusted edges carrying no term */
    size_t tabu_best_iter;           /* iteration of canonical proxy minimum */
    size_t tabu_hash_rebuilds;
    size_t tabu_exact_queries;       /* chart-placement pair queries */
    size_t tabu_exact_cache_entries; /* memoized AABB-surviving pairs */
    size_t tabu_exact_cache_hits;
    size_t tabu_exact_cache_misses;
    size_t tabu_exact_aabb_rejects;
    size_t tabu_exact_sat_tests;     /* triangle pairs reaching narrow phase */
    size_t tabu_archive_states;       /* every accepted state, never evicted */
    size_t tabu_leaderboard;         /* candidate assignments SAT-verified */
    size_t tabu_winner;              /* index within SAT-verified subset */
    size_t tabu_winner_iter;         /* archived search iteration selected */
    double tabu_exact_score;         /* exact-overlap canonical objective */
    double tabu_energy_initial;      /* canonical energy, shared-cell units */
    double tabu_energy_final;
    /* The field's measured wind ladder: signed du per (2*pi*r) of one wind
     * outward.  0 means too few radial pairs and the spiral-map fallback. */
    double tabu_ladder_step;
    double tabu_ladder_u0;
    int tabu_ladder_fixed;
    size_t tabu_ladder_pairs;
    size_t tabu_phase_gauge_islands;
    size_t tabu_phase_flat_groups;    /* trusted groups with one relative bin */
    size_t tabu_phase_gauge_inconsistent;
    size_t tabu_lateral_edges;       /* same-wind pairs coupled beyond welds */
    size_t tabu_phase_lateral_trusted;
    size_t tabu_phase_lateral_nonzero;
    size_t tabu_phase_regions;       /* active compound gauge regions */
    size_t tabu_phase_region_links;  /* accepted group-pair constraints */
    size_t tabu_lateral_torn;        /* of those, left at dk != 0 */
    size_t tabu_family_edges;        /* edges boosted inside happy families */
    size_t tabu_family_torn;         /* of those, left at dk != target */

    /* Relayout (windings-fixed continuous re-registration). */
    size_t relayout_edges_weld;
    size_t relayout_edges_lateral;
    size_t relayout_edges_radial;
    double relayout_weld_rms_before;
    double relayout_weld_rms_after;
    double relayout_radial_rms_before;
    double relayout_radial_rms_after;
    double relayout_metric_shift_rms; /* unconstrained spring preference */
    double relayout_metric_shift_max;
    double relayout_shift_rms;
    double relayout_shift_max;
    double relayout_atlas_width_before;
    double relayout_atlas_width_after;
    size_t relayout_moved_charts;
    size_t relayout_exact_pairs_before;
    size_t relayout_exact_pairs_after;
    size_t relayout_cross_pairs_before;
    size_t relayout_cross_pairs_after;
    size_t relayout_allowed_pairs_before;
    size_t relayout_allowed_pairs_after;
    size_t relayout_hard_pairs_before;
    size_t relayout_hard_pairs_after;
    size_t relayout_topology_limited_before;
    size_t relayout_topology_limited_after;
    size_t relayout_total_bounds;
    size_t relayout_topology_bounds;
    size_t relayout_collision_bounds_rejected;
    size_t relayout_collision_bounds;
    size_t relayout_collision_bounds_added;
    size_t relayout_outer_iterations;
    int    relayout_qp_iterations;
    int    relayout_qp_active_bounds;
    double relayout_qp_objective;
    double relayout_qp_stationarity;
    double relayout_qp_max_violation;
    int    relayout_topology_infeasible;
    int    relayout_collision_failed;
    int    relayout_failed;          /* solve failed; u left as tabu placed */
    int    relayout_reverted;        /* independent exact check rejected it */
    /* Final exact SAT measurement of what tabu could not fix, and its radial
     * signature (same bins as intra_radial_hist). */
    size_t intra_overlap_pairs_after;
    size_t intra_radial_hist_after[6];
} AtlasOverlapFixStats;

/*
 * out_u / out_v receive the corrected per-vertex parameterization.  v starts
 * as the axial coordinate and changes only for parked groups.  face_keep[f] is
 * set to 0 for killed faces; the caller compacts with PieceSet_apply_facekeep.
 *
 * out_vertex_group  [ps->nv]  group id per vertex, -1 where unused.  Optional.
 * out_group/out_ngroups       the per-group decision records.  Optional.
 * tabu_report                 per-chart depths + accepted-move trace, filled
 *                             only in tabu mode.  Optional.
 */
int AtlasOverlapFix_solve(Arena_T arena,
                          const PieceSet *ps,
                          const ScaffoldCalib *cal,
                          const float *registered_u,
                          const AtlasOverlapFixOptions *opts,
                          double *out_u,
                          double *out_v,
                          uint8_t *face_keep,
                          int32_t *out_vertex_group,
                          AtlasOverlapFixGroup **out_group,
                          size_t *out_ngroups,
                          AtlasOverlapFixTabuReport *tabu_report,
                          AtlasOverlapFixStats *stats);

int AtlasOverlapFix_selftest(void);

#endif
