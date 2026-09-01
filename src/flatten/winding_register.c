#define _USE_MATH_DEFINES
#include "winding_register.h"
#include "winding_mrf.h"

#include "../common/union_find.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define WR_TWO_PI             (2.0 * M_PI)
/* Same-V radial neighbours are the strongest observable winding relation on
 * a scroll crop.  Keep their bins close to the mesh sampling scale: the old
 * 8-voxel x 5.6-degree buckets merged unrelated rows/rays and, paradoxically,
 * left small true neighbours without the four consistent votes required to
 * enter the graph.  At 2 vox x 1.4 degrees a relation is local enough to carry
 * a large binary penalty without confusing nearby plies/materials. */
#define WR_PHASE_BINS         256
#define WR_AXIAL_BIN          2.0
#define WR_ORDER_BINARY_WEIGHT 32.0
/* REGIME-BOUND (changelog rule): 2M was tuned when the biggest input was the
 * welded 4x5x5.  The 10x pre-weld soup is 12.4M vertices; at stride 7 a
 * continuation PAIR only survives sampling ~2% of the time (1/49), the
 * relation graph starves (obs/comp 666 -> 37), almost nothing merges
 * (materials 8331 -> 5211 instead of ~200:1), and the atlas shreds into
 * confetti.  16M covered the 10x pre-weld soup, but the genuinely welded
 * 10x control is 18.19M vertices -- just over the cap, so it silently fell
 * back to stride 2 and lost ~3/4 of every continuation pair.  32M keeps the
 * welded 10x at stride 1.  The cap only bounds the sample array; capacity is
 * still nvertices/stride, so a smaller rung allocates no more than before. */
#define WR_MAX_SAMPLES        32000000u
#define WR_CELL_SCAN_LIMIT    96
#define WR_CONT_WTOL          0.30
#define WR_CONT_PHASE_TOL     0.25
#define WR_ORDER_PHASE_TOL    0.25
/* MEASURED BAD 2026-08-31, kept compiled out: multi-pitch order relations
 * (same ray, k>=2 wraps apart, local-pitch snapped) flooded the graph with
 * cross-material junk -- at k>=2 a ray's radial neighbour is frequently a
 * DIFFERENT sheet, so "k wraps out = +k turns" fails without material
 * identity (pile fit: order sat 0.998 -> 0.239, conflicts 14 -> 977, lift
 * span 15.1 -> 21.1 turns; welded fit: sat 0.142, span -> 32.2t).  Re-arm
 * only with a same-material constraint (e.g. the solidified lattice's
 * extended identity). */
#define WR_ORDER_MULTI        0
#define WR_ORDER_MAX_WRAPS    12     /* multi-pitch order reach, in wraps */
#define WR_ORDER_SNAP         0.25   /* x local pitch; abstain outside */
#define WR_PACK_GUTTER_TURNS  0.10
/* Integer loop closure over the relation graph (port of the proven
 * track_assembly._repair_subtree_shifts): a mis-snapped forest edge displaces
 * its whole subtree by a constant, and the non-forest "closer" relations that
 * straddle the cut all disagree by exactly that constant.  Rank candidate
 * forest cuts by the mode-observation mass of currently-conflicted crossing
 * relations, try small integer shifts on the subtree, accept only a STRICT
 * net reduction of conflicted evidence mass (forest edges priced in), best
 * gain wins per round.  MAX_ROUNDS 0 compiles the pass out (A/B control). */
#define WR_REPAIR_MAX_ROUNDS  64
#define WR_REPAIR_MAX_SHIFT   2
#define WR_REPAIR_RANK_CAP    256
#define WR_REPAIR_WORK_CAP    (UINT64_C(1) << 30)
#define WR_MRF_ROUNDS          4
#define WR_MRF_DELTA_RADIUS    2
#define WR_BOUNDARY_AXIAL_SLACK (0.25 * WR_AXIAL_BIN)
#define WR_CUBE_EDGE          128.0
/* Conflict exclusion is an OUTER loop around the ordinary relation MRF.
 * A UV claim is one fine same-V/same-phase winding bin.  This first guarded
 * pass deliberately targets the catastrophic failure we can identify without
 * ambiguity: four or more claimants, including spatially remote pieces,
 * occupying the same UV.  Pitch-scale
 * separation can also be an incorrect neighbouring wrap, but on the raw BPA
 * mesh it is confounded with sheet thickness, cracks, and radial roughness and
 * must not become a license to walk every component through the integer frame.
 *
 * The best-supported remote claimant keeps its current integer gauge; every
 * repeatedly observed loser receives a near-hard sparse unary forbidding that
 * exact gauge on the next graph cut.  Ordinary two/three-layer occupancy stays
 * quadribbon's responsibility.  Crucially, each outer solve freezes every
 * component which is not sitting on a forbidden assignment.  The losing batch
 * can therefore move together and use its relations to the already-good UVs,
 * but it cannot drag those good UVs into a global cascade.  Moves stay within
 * the same +/-8-turn trust region as the original recentered MRF. */
#define WR_CONFLICT_MAX_ROUNDS              16
#define WR_CONFLICT_MIN_BINS                 1
#define WR_CONFLICT_MIN_CLAIMANTS            4
#define WR_CONFLICT_RADIUS_FRACTION          0.55
#define WR_CONFLICT_REMOTE_PITCHES           8.0
#define WR_CONFLICT_REMOTE_MIN_VOX          48.0
#define WR_CONFLICT_TRUST_RADIUS              8
#define WR_CONFLICT_PRIOR_WEIGHT              2.0
#define WR_CONFLICT_UNARY_WEIGHT          65536.0
/* Absolute Gaussian prior used only by the first relation-MRF stage.  It is
 * intentionally stronger than the sparse conflict-retry prior: seam evidence
 * may split a locally contradicted continuation cycle, but must not translate
 * a macroscopic fraction of the scroll through four recentered rounds. */
#define WR_MRF_BASE_PRIOR_WEIGHT              32.0
/* Continuation is the binary same-winding certificate.  In particular, a
 * run of adjacent observations at the same V should move as a unit.  The
 * Kruskal score carries a +1000 ordering sentinel; after removing it, retain
 * enough evidence mass that a Gaussian trust unary cannot cheaply staircase
 * a coherent run one component at a time.  The Artifact C exact fixture puts
 * the lower bound at 10.5; 16 leaves a useful margin while still allowing a
 * weak false continuation to be cut by two independent seam-order factors. */
#define WR_MRF_CONTINUATION_SCALE              16.0

typedef struct {
    int32_t axial_bin, base_axial_bin, phase_bin, component, raw_turn;
    uint32_t count;
    double radius, q;
} WrStrand;

enum {
    WR_OBS_CONTINUATION = 0,
    WR_OBS_ORDER = 1,
    WR_OBS_SEAM_ORDER = 2
};

static const char *wr_observation_kind_name(int kind)
{
    return kind == WR_OBS_CONTINUATION ? "continuation" :
           kind == WR_OBS_SEAM_ORDER ? "order_seam" : "order";
}

/* target means correction[b] - correction[a] == target. */
typedef struct {
    int32_t a, b, target;
    float residual;
    uint8_t kind;
} WrObservation;

typedef struct {
    int32_t a, b, target;
    size_t observations, mode_observations;
    double agreement, residual, weight;
    uint8_t kind, eligible, selected;
} WrRelation;

typedef struct { int32_t other, target, next; } WrAdjacency;

typedef struct {
    int32_t x, y, z, head;
    uint8_t used;
} WrCell;

typedef struct {
    int32_t root;
    int64_t vertices;
    double qmin, qmax;
} WrPackGroup;

typedef struct {
    int valid;
    int32_t target;
    size_t observations, mode_observations;
    double agreement, residual, weight;
    uint8_t kind, eligible;
} WrCandidate;

typedef struct {
    int32_t axial_bin, phase_bin, component, root;
    int64_t absolute_turn;
    uint32_t count;
    double radius;
} WrConflictClaim;

typedef struct {
    int32_t component;
    int32_t correction;
} WrForbiddenWinding;

/* Forensic trace for real-data conflict-resolution audits.  This is kept
 * behind an environment switch because a fragmented full-scroll solve can
 * legitimately contain thousands of components; the ordinary production log
 * should remain compact.  The trace is observational only. */
static int wr_trace_conflicts(void)
{
    const char *value = getenv("VES_WINDING_TRACE_CONFLICTS");
    return value != NULL && value[0] != '\0' && strcmp(value, "0") != 0;
}

static double wr_huber(double value, double delta);

/* Opt-in forensic dump of the actual MRF presented to graph cut.  This is
 * intentionally an environment switch rather than production output: a full
 * scroll can have many thousands of sites, while the ordinary pipeline only
 * needs the compact aggregate statistics above.  The caller creates the
 * directory named by VES_WINDING_DUMP_TERMS. */
static const char *wr_dump_terms_dir(void)
{
    const char *value = getenv("VES_WINDING_DUMP_TERMS");
    return value != NULL && value[0] != '\0' && strcmp(value, "0") != 0
         ? value : NULL;
}

static FILE *wr_open_term_dump(
    const char *stem, int solve_index, int round_index)
{
    const char *dir = wr_dump_terms_dir();
    char path[4096];
    int count;
    if (dir == NULL || stem == NULL) return NULL;
    if (solve_index >= 0 && round_index >= 0)
        count = snprintf(path, sizeof path, "%s/%s_s%02d_r%02d.csv",
                         dir, stem, solve_index, round_index);
    else
        count = snprintf(path, sizeof path, "%s/%s", dir, stem);
    if (count < 0 || (size_t)count >= sizeof path) return NULL;
    return fopen(path, "wb");
}

static double wr_debug_site_unary(
    const WindingMRFSite *site, int32_t label,
    const WindingMRFOptions *options)
{
    /* Must mirror winding_mrf.c::wm_unary_energy exactly. */
    if (site->fixed)
        return label == site->initial_label ? 0.0
             : (double)(INT32_MAX / 32);
    if (!(site->weight > 0.0) || !(site->sigma > 0.0) ||
        !isfinite(site->center) || !isfinite(site->sigma) ||
        !isfinite(site->weight))
        return 0.0;
    return site->weight * wr_huber(
        ((double)label - site->center) / site->sigma,
        options->unary_huber_delta);
}

static double wr_debug_sparse_unary(
    const WindingMRFUnaryPenalty *penalty, size_t npenalty,
    int32_t component, int32_t label)
{
    double total = 0.0;
    for (size_t i = 0; i < npenalty; i++)
        if (penalty[i].site == component && penalty[i].label == label)
            total += penalty[i].weight;
    return total;
}

static int32_t wr_debug_quantize_cost(double value, double scale)
{
    const int32_t cap = INT32_MAX / 32;
    double cost = value * scale;
    if (!(cost > 0.0)) return 0;
    if (!isfinite(cost) || cost >= (double)cap) return cap;
    return (int32_t)floor(cost + 0.5);
}

static void wr_dump_mrf_round(
    int solve_index, int round_index,
    const WindingMRFSite *site, int32_t ncomponents,
    const WrRelation *relation, size_t nrelation,
    const WindingMRFEdge *edge, size_t nedge,
    const WindingMRFUnaryPenalty *penalty, size_t npenalty,
    const WindingMRFOptions *options,
    const int32_t *correction, const int32_t *delta,
    const float *confidence, const WindingMRFStats *stats)
{
    FILE *sites = wr_open_term_dump("mrf_sites", solve_index, round_index);
    FILE *edges = wr_open_term_dump("mrf_edges", solve_index, round_index);
    FILE *summary = wr_open_term_dump("mrf_summary", solve_index, round_index);
    if (sites != NULL) {
        fprintf(sites,
            "component,label,correction_before,candidate_correction,"
            "center,sigma,weight,fixed,initial_label,base_unary,"
            "sparse_unary,total_unary,quantized_unary,solved_delta,"
            "selected,confidence\n");
        for (int32_t c = 0; c < ncomponents; c++) {
            for (int32_t label = options->label_min;
                 label <= options->label_max; label++) {
                double base = wr_debug_site_unary(&site[c], label, options);
                double sparse = wr_debug_sparse_unary(
                    penalty, npenalty, c, label);
                double total = base + sparse;
                fprintf(sites,
                    "%d,%d,%d,%lld,%.17g,%.17g,%.17g,%d,%d,"
                    "%.17g,%.17g,%.17g,%d,%d,%d,%.9g\n",
                    c, label, correction[c],
                    (long long)((int64_t)correction[c] + label),
                    site[c].center, site[c].sigma, site[c].weight,
                    site[c].fixed, site[c].initial_label,
                    base, sparse, total,
                    wr_debug_quantize_cost(total, options->cost_scale),
                    delta[c], label == delta[c],
                    confidence != NULL ? confidence[c] : NAN);
            }
        }
        fclose(sites);
    }
    if (edges != NULL) {
        size_t at = 0;
        fprintf(edges,
            "relation,a,b,kind,absolute_target,correction_a,correction_b,"
            "solved_before,residual_target,weight,quantized_weight,"
            "mode_observations,observations,agreement,measurement_residual,"
            "forest,delta_a,delta_b,solved_after,residual_before,"
            "residual_after,energy_before,energy_after\n");
        for (size_t i = 0; i < nrelation; i++) {
            const WrRelation *current = &relation[i];
            const WindingMRFEdge *factor;
            int64_t solved_before, solved_after, residual_before, residual_after;
            if (!current->eligible) continue;
            if (at >= nedge) break;
            factor = &edge[at++];
            solved_before = (int64_t)correction[current->b] -
                            correction[current->a];
            solved_after = solved_before +
                (int64_t)delta[current->b] - delta[current->a];
            residual_before = solved_before - current->target;
            residual_after = solved_after - current->target;
            if (residual_before < 0) residual_before = -residual_before;
            if (residual_after < 0) residual_after = -residual_after;
            fprintf(edges,
                "%zu,%d,%d,%s,%d,%d,%d,%lld,%d,%.17g,%d,%zu,%zu,"
                "%.17g,%.17g,%d,%d,%d,%lld,%lld,%lld,%.17g,%.17g\n",
                i, current->a, current->b,
                wr_observation_kind_name(current->kind),
                current->target, correction[current->a], correction[current->b],
                (long long)solved_before, factor->target, factor->weight,
                wr_debug_quantize_cost(factor->weight, options->cost_scale),
                current->mode_observations, current->observations,
                current->agreement, current->residual, current->selected,
                delta[current->a], delta[current->b], (long long)solved_after,
                (long long)residual_before, (long long)residual_after,
                factor->weight * (double)residual_before,
                factor->weight * (double)residual_after);
        }
        fclose(edges);
    }
    if (summary != NULL) {
        fprintf(summary, "key,value\n");
        fprintf(summary, "solve,%d\nround,%d\n", solve_index, round_index);
        fprintf(summary, "sites,%d\nedges,%zu\n", ncomponents, nedge);
        fprintf(summary, "labels,%lld\n",
                (long long)options->label_max - options->label_min + 1);
        fprintf(summary, "cost_scale,%.17g\n", options->cost_scale);
        fprintf(summary, "energy_before,%.17g\n", stats->energy_before);
        fprintf(summary, "energy_after,%.17g\n", stats->energy_after);
        fprintf(summary, "quantized_energy_before,%lld\n",
                (long long)stats->quantized_energy_before);
        fprintf(summary, "quantized_energy_after,%lld\n",
                (long long)stats->quantized_energy_after);
        fprintf(summary, "changed_labels,%zu\n", stats->changed_labels);
        fprintf(summary, "mean_confidence,%.17g\n", stats->mean_confidence);
        fclose(summary);
    }
}

static int wr_mrf_dump_serial = 0;

static void wr_dump_final_terms(
    const float *vertices, size_t nvertices, const double *q,
    const int32_t *component, int32_t ncomponents,
    const int32_t *component_size, int32_t anchor_component,
    const int32_t *correction,
    const WindingRegisterFieldUnary *field_unary,
    const WrRelation *relation, size_t nrelation, UnionFind *graph)
{
    FILE *components = NULL, *relations = NULL, *raw = NULL;
    double *range;
    if (wr_dump_terms_dir() == NULL || ncomponents <= 0) return;
    if ((size_t)ncomponents > SIZE_MAX / (8 * sizeof(double))) return;
    range = (double *)malloc((size_t)ncomponents * 8 * sizeof(double));
    if (range == NULL) return;
    for (int32_t c = 0; c < ncomponents; c++) {
        range[(size_t)c*8+0] = DBL_MAX;
        range[(size_t)c*8+1] = -DBL_MAX;
        range[(size_t)c*8+2] = DBL_MAX;
        range[(size_t)c*8+3] = -DBL_MAX;
        range[(size_t)c*8+4] = DBL_MAX;
        range[(size_t)c*8+5] = -DBL_MAX;
        range[(size_t)c*8+6] = DBL_MAX;
        range[(size_t)c*8+7] = -DBL_MAX;
    }
    for (size_t i = 0; i < nvertices; i++) {
        int32_t c = component[i];
        double *r;
        if (c < 0 || c >= ncomponents) continue;
        r = &range[(size_t)c*8];
        if (vertices[i*3+0] < r[0]) r[0] = vertices[i*3+0];
        if (vertices[i*3+0] > r[1]) r[1] = vertices[i*3+0];
        if (vertices[i*3+1] < r[2]) r[2] = vertices[i*3+1];
        if (vertices[i*3+1] > r[3]) r[3] = vertices[i*3+1];
        if (vertices[i*3+2] < r[4]) r[4] = vertices[i*3+2];
        if (vertices[i*3+2] > r[5]) r[5] = vertices[i*3+2];
        if (q[i] < r[6]) r[6] = q[i];
        if (q[i] > r[7]) r[7] = q[i];
    }
    components = wr_open_term_dump("components.csv", -1, -1);
    if (components != NULL) {
        fprintf(components,
            "component,vertices,anchor,graph_root,z_min,z_max,y_min,y_max,"
            "x_min,x_max,q_min,q_max,correction,corrected_q_min,"
            "corrected_q_max,field_center,field_sigma,field_weight\n");
        for (int32_t c = 0; c < ncomponents; c++) {
            const double *r = &range[(size_t)c*8];
            double field_center = field_unary != NULL
                                ? field_unary->center[c] : NAN;
            double field_sigma = field_unary != NULL
                               ? field_unary->sigma[c] : NAN;
            double field_weight = field_unary != NULL
                                ? field_unary->weight[c] : NAN;
            fprintf(components,
                "%d,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
                "%.17g,%.17g,%d,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                c, component_size[c], c == anchor_component,
                uf_find(graph, c),
                r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7],
                correction[c], r[6] + correction[c],
                r[7] + correction[c], field_center, field_sigma, field_weight);
        }
        fclose(components);
    }
    relations = wr_open_term_dump("relations_final.csv", -1, -1);
    if (relations != NULL) {
        fprintf(relations,
            "relation,a,b,kind,target,correction_a,correction_b,solved,"
            "integer_residual,weight,mode_observations,observations,agreement,"
            "measurement_residual,eligible,forest,graph_root_a,graph_root_b\n");
        for (size_t i = 0; i < nrelation; i++) {
            const WrRelation *current = &relation[i];
            int64_t solved = (int64_t)correction[current->b] -
                             correction[current->a];
            int64_t residual = solved - current->target;
            if (residual < 0) residual = -residual;
            fprintf(relations,
                "%zu,%d,%d,%s,%d,%d,%d,%lld,%lld,%.17g,%zu,%zu,"
                "%.17g,%.17g,%d,%d,%d,%d\n",
                i, current->a, current->b,
                wr_observation_kind_name(current->kind),
                current->target, correction[current->a], correction[current->b],
                (long long)solved, (long long)residual, current->weight,
                current->mode_observations, current->observations,
                current->agreement, current->residual, current->eligible,
                current->selected, uf_find(graph, current->a),
                uf_find(graph, current->b));
        }
        fclose(relations);
    }
    raw = wr_open_term_dump("mesh_component.i32", -1, -1);
    if (raw != NULL) {
        fwrite(component, sizeof(*component), nvertices, raw);
        fclose(raw);
    }
    raw = wr_open_term_dump("turn_correction.i32", -1, -1);
    if (raw != NULL) {
        fwrite(correction, sizeof(*correction), (size_t)ncomponents, raw);
        fclose(raw);
    }
    free(range);
}

static double wr_wrap_to_pi(double angle)
{
    double value = fmod(angle + M_PI, WR_TWO_PI);
    if (value < 0.0) value += WR_TWO_PI;
    return value - M_PI;
}

static double wr_huber(double value, double delta)
{
    value = fabs(value);
    if (value <= delta) return 0.5 * value * value;
    return delta * (value - 0.5 * delta);
}

static int wr_compare_strand_key(const void *pa, const void *pb)
{
    const WrStrand *a = (const WrStrand *)pa;
    const WrStrand *b = (const WrStrand *)pb;
    if (a->axial_bin != b->axial_bin)
        return a->axial_bin < b->axial_bin ? -1 : 1;
    if (a->phase_bin != b->phase_bin)
        return a->phase_bin < b->phase_bin ? -1 : 1;
    if (a->component != b->component)
        return a->component < b->component ? -1 : 1;
    if (a->raw_turn != b->raw_turn)
        return a->raw_turn < b->raw_turn ? -1 : 1;
    if (a->base_axial_bin != b->base_axial_bin)
        return a->base_axial_bin < b->base_axial_bin ? -1 : 1;
    if (a->radius != b->radius) return a->radius < b->radius ? -1 : 1;
    return 0;
}

static int wr_compare_strand_radius(const void *pa, const void *pb)
{
    const WrStrand *a = (const WrStrand *)pa;
    const WrStrand *b = (const WrStrand *)pb;
    if (a->axial_bin != b->axial_bin)
        return a->axial_bin < b->axial_bin ? -1 : 1;
    if (a->phase_bin != b->phase_bin)
        return a->phase_bin < b->phase_bin ? -1 : 1;
    if (a->radius != b->radius) return a->radius < b->radius ? -1 : 1;
    if (a->component != b->component)
        return a->component < b->component ? -1 : 1;
    if (a->raw_turn != b->raw_turn)
        return a->raw_turn < b->raw_turn ? -1 : 1;
    return a->base_axial_bin < b->base_axial_bin ? -1 :
           (a->base_axial_bin > b->base_axial_bin ? 1 : 0);
}

static int wr_compare_observation(const void *pa, const void *pb)
{
    const WrObservation *a = (const WrObservation *)pa;
    const WrObservation *b = (const WrObservation *)pb;
    if (a->a != b->a) return a->a < b->a ? -1 : 1;
    if (a->b != b->b) return a->b < b->b ? -1 : 1;
    if (a->kind != b->kind) return a->kind < b->kind ? -1 : 1;
    if (a->target != b->target) return a->target < b->target ? -1 : 1;
    return a->residual < b->residual ? -1 :
           (a->residual > b->residual ? 1 : 0);
}

static int wr_compare_relation_weight(const void *pa, const void *pb)
{
    const WrRelation *a = (const WrRelation *)pa;
    const WrRelation *b = (const WrRelation *)pb;
    if (a->eligible != b->eligible) return a->eligible ? -1 : 1;
    if (a->weight != b->weight) return a->weight > b->weight ? -1 : 1;
    if (a->kind != b->kind) return a->kind < b->kind ? -1 : 1;
    if (a->a != b->a) return a->a < b->a ? -1 : 1;
    return a->b < b->b ? -1 : (a->b > b->b ? 1 : 0);
}

static int wr_compare_pack_group(const void *pa, const void *pb)
{
    const WrPackGroup *a = (const WrPackGroup *)pa;
    const WrPackGroup *b = (const WrPackGroup *)pb;
    if (a->vertices != b->vertices)
        return a->vertices > b->vertices ? -1 : 1;
    return a->root < b->root ? -1 : (a->root > b->root ? 1 : 0);
}

static int wr_compare_conflict_claim(const void *pa, const void *pb)
{
    const WrConflictClaim *a = (const WrConflictClaim *)pa;
    const WrConflictClaim *b = (const WrConflictClaim *)pb;
    if (a->axial_bin != b->axial_bin)
        return a->axial_bin < b->axial_bin ? -1 : 1;
    if (a->absolute_turn != b->absolute_turn)
        return a->absolute_turn < b->absolute_turn ? -1 : 1;
    if (a->phase_bin != b->phase_bin)
        return a->phase_bin < b->phase_bin ? -1 : 1;
    if (a->root != b->root) return a->root < b->root ? -1 : 1;
    if (a->radius != b->radius) return a->radius < b->radius ? -1 : 1;
    return a->component < b->component ? -1 :
           (a->component > b->component ? 1 : 0);
}

typedef struct {
    int64_t weight;      /* crossing conflicted evidence mass */
    int32_t component;   /* subtree seed (cut below its parent edge) */
} WrCut;

static int wr_compare_cut(const void *pa, const void *pb)
{
    const WrCut *a = (const WrCut *)pa;
    const WrCut *b = (const WrCut *)pb;
    if (a->weight != b->weight) return a->weight > b->weight ? -1 : 1;
    return a->component < b->component ? -1 :
           (a->component > b->component ? 1 : 0);
}

static uint64_t wr_hash3(int32_t x, int32_t y, int32_t z)
{
    uint64_t hash = UINT64_C(1469598103934665603);
    uint32_t value[3] = { (uint32_t)x, (uint32_t)y, (uint32_t)z };
    for (int k = 0; k < 3; k++) {
        hash ^= (uint64_t)value[k];
        hash *= UINT64_C(1099511628211);
        hash ^= hash >> 32;
    }
    return hash;
}

static WrCell *wr_cell_find(WrCell *cell, size_t mask,
                            int32_t x, int32_t y, int32_t z, int create)
{
    size_t slot = (size_t)wr_hash3(x, y, z) & mask;
    for (size_t guard = 0; guard <= mask; guard++) {
        WrCell *current = &cell[slot];
        if (!current->used) {
            if (!create) return NULL;
            current->used = 1;
            current->x = x; current->y = y; current->z = z;
            current->head = -1;
            return current;
        }
        if (current->x == x && current->y == y && current->z == z)
            return current;
        slot = (slot + 1) & mask;
    }
    return NULL;
}

static int wr_cell_coordinate(float value, double inverse, int32_t *out)
{
    double coordinate = floor((double)value * inverse);
    if (!isfinite(coordinate) || coordinate < (double)INT32_MIN ||
        coordinate > (double)INT32_MAX)
        return -1;
    *out = (int32_t)coordinate;
    return 0;
}

static void wr_append_observation(
    WrObservation *observation, size_t capacity,
    size_t *count, size_t *dropped,
    int32_t a, int32_t b, int32_t target, double residual, int kind)
{
    if (a == b || !isfinite(residual)) return;
    if (a > b) {
        int32_t swap = a; a = b; b = swap; target = -target;
    }
    if (*count >= capacity) {
        (*dropped)++;
        return;
    }
    observation[*count].a = a;
    observation[*count].b = b;
    observation[*count].target = target;
    observation[*count].residual = (float)residual;
    observation[*count].kind = (uint8_t)kind;
    (*count)++;
}

static WrCandidate wr_observation_mode(
    const WrObservation *observation, size_t first, size_t last, int kind)
{
    WrCandidate candidate;
    memset(&candidate, 0, sizeof(candidate));
    candidate.valid = first < last;
    candidate.kind = (uint8_t)kind;
    candidate.observations = last - first;
    size_t best_first = first, best_count = 0;
    int32_t best_target = 0;
    for (size_t mode_first = first; mode_first < last;) {
        size_t mode_last = mode_first + 1;
        while (mode_last < last &&
               observation[mode_last].target ==
               observation[mode_first].target)
            mode_last++;
        size_t count = mode_last - mode_first;
        int take = count > best_count;
        if (count == best_count) {
            int64_t current_abs = observation[mode_first].target < 0
                ? -(int64_t)observation[mode_first].target
                :  (int64_t)observation[mode_first].target;
            int64_t best_abs = best_target < 0
                ? -(int64_t)best_target : (int64_t)best_target;
            take = current_abs < best_abs;
        }
        if (take) {
            best_first = mode_first;
            best_count = count;
            best_target = observation[mode_first].target;
        }
        mode_first = mode_last;
    }
    candidate.target = best_target;
    candidate.mode_observations = best_count;
    candidate.agreement = candidate.observations
        ? (double)best_count / (double)candidate.observations : 0.0;
    candidate.residual = best_count
        ? (double)observation[best_first + best_count / 2].residual : DBL_MAX;
    if (kind == WR_OBS_CONTINUATION) {
        candidate.eligible = best_count >= 3 &&
            candidate.agreement >= 0.75 && candidate.residual <= 0.20;
        candidate.weight = 1000.0 +
            log1p((double)best_count) * candidate.agreement /
            (1.0 + 4.0 * candidate.residual);
    } else {
        candidate.eligible = best_count >= 4 &&
            candidate.agreement >= 0.60 && candidate.residual <= 0.20;
        /* This is a shifted equality factor, not a weak radial heuristic:
         * within one fine V/ray bucket, adjacent layers must satisfy the
         * integer delta implied by their lifted local phases.  Price it well
         * above field unaries and ordinary loop closers.  Continuation keeps
         * its +1000 hierarchy sentinel and therefore still wins any genuine
         * same-sheet contradiction. */
        candidate.weight = WR_ORDER_BINARY_WEIGHT *
            log1p((double)best_count) * candidate.agreement /
            (1.0 + 4.0 * candidate.residual);
    }
    return candidate;
}

static int wr_fill_strands(
    WrStrand *strand, size_t nsample,
    const int32_t *sample_vertex, size_t nvertices,
    const double *axial, const double *radius, const double *q,
    const int32_t *component, int32_t ncomponents,
    double axial_min, double axial_offset)
{
    if (!isfinite(axial_min) || !isfinite(axial_offset)) return -1;
    for (size_t si = 0; si < nsample; si++) {
        int32_t vertex = sample_vertex[si];
        if (vertex < 0 || (size_t)vertex >= nvertices) return -1;
        size_t i = (size_t)vertex;
        int32_t c = component[i];
        double qi = q[i];
        if (c < 0 || c >= ncomponents || !isfinite(qi) ||
            !isfinite(radius[i]) || !isfinite(axial[i]))
            return -1;
        double phase = qi - floor(qi);
        double base_axial_value = floor(
            (axial[i] - axial_min) / WR_AXIAL_BIN);
        double axial_value = floor(
            (axial[i] - axial_min - axial_offset) / WR_AXIAL_BIN);
        double raw_turn_value = floor(qi + 1e-8);
        if (!isfinite(base_axial_value) || !isfinite(axial_value) ||
            base_axial_value < (double)INT32_MIN ||
            base_axial_value > (double)INT32_MAX ||
            axial_value < (double)INT32_MIN ||
            axial_value > (double)INT32_MAX ||
            raw_turn_value < (double)INT32_MIN ||
            raw_turn_value > (double)INT32_MAX)
            return -1;
        int phase_bin = (int)floor(phase * (double)WR_PHASE_BINS);
        if (phase_bin < 0) phase_bin = 0;
        if (phase_bin >= WR_PHASE_BINS) phase_bin = WR_PHASE_BINS - 1;
        strand[si].axial_bin = (int32_t)axial_value;
        strand[si].base_axial_bin = (int32_t)base_axial_value;
        strand[si].phase_bin = phase_bin;
        strand[si].component = c;
        strand[si].raw_turn = (int32_t)raw_turn_value;
        strand[si].count = 1;
        strand[si].radius = radius[i];
        strand[si].q = qi;
    }
    return 0;
}

static int wr_build_samples(
    Arena_T arena, size_t nvertices, const double *axial,
    const double *radius, const double *q, const int32_t *component,
    int32_t ncomponents, double axial_min,
    WrStrand **out_strand, int32_t **out_vertex, size_t *out_count)
{
    size_t stride = nvertices > WR_MAX_SAMPLES
                  ? (nvertices + WR_MAX_SAMPLES - 1) / WR_MAX_SAMPLES : 1;
    size_t capacity = (nvertices + stride - 1) / stride;
    if ((size_t)ncomponents > SIZE_MAX - capacity) return -1;
    capacity += (size_t)ncomponents; /* at least one sample per component */
    if (capacity > (size_t)LONG_MAX / sizeof(WrStrand) ||
        capacity > (size_t)LONG_MAX / sizeof(int32_t))
        return -1;
    WrStrand *strand = (WrStrand *)ARENA_ALLOC(
        arena, (long)(capacity * sizeof(WrStrand)));
    int32_t *sample_vertex = (int32_t *)ARENA_ALLOC(
        arena, (long)(capacity * sizeof(int32_t)));
    uint8_t *sampled = (uint8_t *)ARENA_CALLOC(
        arena, (long)ncomponents, (long)sizeof(uint8_t));
    size_t count = 0;
    for (size_t i = 0; i < nvertices; i++) {
        int32_t c = component[i];
        if (c < 0 || c >= ncomponents) return -1;
        if (i % stride != 0 && sampled[c]) continue;
        if (i > (size_t)INT32_MAX) return -1;
        sampled[c] = 1;
        sample_vertex[count] = (int32_t)i;
        count++;
    }
    if (count == 0 || wr_fill_strands(
            strand, count, sample_vertex, nvertices,
            axial, radius, q, component, ncomponents,
            axial_min, 0.0) != 0)
        return -1;
    *out_strand = strand;
    *out_vertex = sample_vertex;
    *out_count = count;
    return 0;
}

static int wr_component_axial_extents(
    Arena_T arena, const double *axial, const int32_t *component,
    size_t nvertices, int32_t ncomponents,
    double **out_minimum, double **out_maximum)
{
    double *minimum = (double *)ARENA_ALLOC(
        arena, (long)((size_t)ncomponents * sizeof(double)));
    double *maximum = (double *)ARENA_ALLOC(
        arena, (long)((size_t)ncomponents * sizeof(double)));
    for (int32_t c = 0; c < ncomponents; c++) {
        minimum[c] = DBL_MAX;
        maximum[c] = -DBL_MAX;
    }
    for (size_t i = 0; i < nvertices; i++) {
        int32_t c = component[i];
        if (c < 0 || c >= ncomponents || !isfinite(axial[i])) return -1;
        if (axial[i] < minimum[c]) minimum[c] = axial[i];
        if (axial[i] > maximum[c]) maximum[c] = axial[i];
    }
    for (int32_t c = 0; c < ncomponents; c++)
        if (minimum[c] == DBL_MAX || maximum[c] == -DBL_MAX) return -1;
    *out_minimum = minimum;
    *out_maximum = maximum;
    return 0;
}

static void wr_collect_continuations(
    const float *vertices, const double *radius, const double *theta,
    const double *q, const int32_t *component,
    const int32_t *sample_vertex, size_t nsample,
    double pitch, int winding_sense,
    WrObservation *observation, size_t capacity,
    size_t continuation_limit, size_t *count, size_t *dropped)
{
    if (!(pitch > 1e-6) || !isfinite(pitch) ||
        (winding_sense != -1 && winding_sense != 1))
        return;
    double gap = 0.75 * pitch;
    if (gap < 3.0) gap = 3.0;
    if (gap > 8.0) gap = 8.0;
    double gap2 = gap * gap, inverse = 1.0 / gap;
    if (nsample > SIZE_MAX / 2) return;
    size_t required = 2 * nsample, slots = 1;
    while (slots < required) {
        if (slots > SIZE_MAX / 2) return;
        slots <<= 1;
    }
    if (slots > SIZE_MAX / sizeof(WrCell) ||
        nsample > SIZE_MAX / sizeof(int32_t))
        return;
    WrCell *cell = (WrCell *)calloc(slots, sizeof(WrCell));
    int32_t *next = (int32_t *)malloc(nsample * sizeof(int32_t));
    if (cell == NULL || next == NULL) {
        free(next); free(cell);
        return; /* fail-soft: ray-order evidence remains available */
    }
    for (size_t si = 0; si < nsample; si++) {
        int32_t vertex = sample_vertex[si], ix, iy, iz;
        next[si] = -1;
        if (wr_cell_coordinate(vertices[(size_t)vertex*3+0],
                               inverse, &iz) != 0 ||
            wr_cell_coordinate(vertices[(size_t)vertex*3+1],
                               inverse, &iy) != 0 ||
            wr_cell_coordinate(vertices[(size_t)vertex*3+2],
                               inverse, &ix) != 0)
            continue;
        WrCell *slot = wr_cell_find(cell, slots - 1, ix, iy, iz, 1);
        if (slot == NULL) continue;
        next[si] = slot->head;
        slot->head = (int32_t)si;
    }
    for (size_t si = 0; si < nsample && *count < continuation_limit; si++) {
        int32_t vi = sample_vertex[si], ix, iy, iz;
        if (wr_cell_coordinate(vertices[(size_t)vi*3+0],
                               inverse, &iz) != 0 ||
            wr_cell_coordinate(vertices[(size_t)vi*3+1],
                               inverse, &iy) != 0 ||
            wr_cell_coordinate(vertices[(size_t)vi*3+2],
                               inverse, &ix) != 0)
            continue;
        for (int dz = -1; dz <= 1; dz++)
        for (int dy = -1; dy <= 1; dy++)
        for (int dx = -1; dx <= 1; dx++) {
            WrCell *slot = wr_cell_find(
                cell, slots - 1, ix + dx, iy + dy, iz + dz, 0);
            if (slot == NULL) continue;
            int scanned = 0;
            for (int32_t sj = slot->head;
                 sj >= 0 && scanned < WR_CELL_SCAN_LIMIT;
                 sj = next[sj], scanned++) {
                if ((size_t)sj <= si) continue;
                int32_t vj = sample_vertex[(size_t)sj];
                int32_t ci = component[vi], cj = component[vj];
                if (ci == cj) continue;
                double d0 = (double)vertices[(size_t)vi*3+0] -
                            (double)vertices[(size_t)vj*3+0];
                double d1 = (double)vertices[(size_t)vi*3+1] -
                            (double)vertices[(size_t)vj*3+1];
                double d2 = (double)vertices[(size_t)vi*3+2] -
                            (double)vertices[(size_t)vj*3+2];
                if (d0*d0 + d1*d1 + d2*d2 > gap2) continue;
                double dtheta = wr_wrap_to_pi(theta[vj] - theta[vi]);
                double helix = (radius[vj] - radius[vi]) / pitch -
                    (double)winding_sense * dtheta / WR_TWO_PI;
                if (fabs(helix) > WR_CONT_WTOL) continue;
                long target_long = lround(q[vi] - q[vj]);
                if (target_long < INT32_MIN || target_long > INT32_MAX)
                    continue;
                double residual =
                    fabs(q[vj] + (double)target_long - q[vi]);
                if (residual > WR_CONT_PHASE_TOL) continue;
                wr_append_observation(
                    observation, capacity, count, dropped, ci, cj,
                    (int32_t)target_long, residual, WR_OBS_CONTINUATION);
            }
        }
    }
    free(next);
    free(cell);
}

static size_t wr_collapse_strands(WrStrand *strand, size_t nsample)
{
    qsort(strand, nsample, sizeof(WrStrand), wr_compare_strand_key);
    size_t nstrand = 0;
    for (size_t first = 0; first < nsample;) {
        size_t last = first + 1;
        while (last < nsample &&
               strand[last].axial_bin == strand[first].axial_bin &&
               strand[last].phase_bin == strand[first].phase_bin &&
               strand[last].component == strand[first].component &&
               strand[last].raw_turn == strand[first].raw_turn &&
               strand[last].base_axial_bin ==
                   strand[first].base_axial_bin)
            last++;
        WrStrand merged = strand[first];
        uint64_t total = 0;
        merged.radius = 0.0;
        merged.q = 0.0;
        for (size_t i = first; i < last; i++) {
            total += strand[i].count;
            merged.radius += strand[i].radius * (double)strand[i].count;
            merged.q += strand[i].q * (double)strand[i].count;
        }
        if (total == 0) return 0;
        merged.radius /= (double)total;
        merged.q /= (double)total;
        merged.count = total > UINT32_MAX ? UINT32_MAX : (uint32_t)total;
        strand[nstrand++] = merged;
        first = last;
    }
    qsort(strand, nstrand, sizeof(WrStrand), wr_compare_strand_radius);
    return nstrand;
}

static int wr_compare_gap(const void *pa, const void *pb)
{
    double a = *(const double *)pa, b = *(const double *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

/* A shifted-lattice order vote may override continuation closure only when
 * the mesh itself certifies an owned 128-voxel cube-face break: one component
 * terminates on the lower side and the other begins on the upper side.
 * Merely crossing one of the ubiquitous 2-voxel bin boundaries while both
 * components continue through it is not seam evidence. */
static int wr_is_axial_boundary_pair(
    const WrStrand *a, const WrStrand *b, double axial_min,
    const double *component_axial_min,
    const double *component_axial_max)
{
    const WrStrand *lower = a, *upper = b;
    if (a->base_axial_bin == b->base_axial_bin) return 0;
    if (lower->base_axial_bin > upper->base_axial_bin) {
        lower = b;
        upper = a;
    }
    if ((int64_t)upper->base_axial_bin - lower->base_axial_bin != 1)
        return 0;
    double boundary = axial_min +
        ((double)lower->base_axial_bin + 1.0) * WR_AXIAL_BIN;
    double cube_face = nearbyint(boundary / WR_CUBE_EDGE) * WR_CUBE_EDGE;
    if (fabs(boundary - cube_face) > WR_BOUNDARY_AXIAL_SLACK) return 0;
    return component_axial_max[lower->component] <=
               boundary + WR_BOUNDARY_AXIAL_SLACK &&
           component_axial_min[upper->component] >=
               boundary - WR_BOUNDARY_AXIAL_SLACK;
}

static size_t wr_collect_order_lattice(
    const WrStrand *strand, size_t nstrand, double pitch,
    int seam_lattice, double axial_min,
    const double *component_axial_min,
    const double *component_axial_max,
    WrObservation *observation, size_t capacity,
    size_t *count, size_t *dropped)
{
    double gap_min = pitch > 1e-6 ? 0.35 * pitch : 0.0;
    double gap_max = pitch > 1e-6 ? 2.25 * pitch : DBL_MAX;
    size_t bins = 0;
    for (size_t first = 0; first < nstrand;) {
        size_t last = first + 1;
        double adjacent_gap[64];
        size_t nadjacent = 0;
        double local_pitch = pitch;
        while (last < nstrand &&
               strand[last].axial_bin == strand[first].axial_bin &&
               strand[last].phase_bin == strand[first].phase_bin)
            last++;
        bins++;
        for (size_t i = first; i + 1 < last; i++) {
            const WrStrand *inner = &strand[i];
            const WrStrand *outer = &strand[i + 1];
            if (inner->component == outer->component) continue;
            double gap = outer->radius - inner->radius;
            if (gap >= gap_min && gap <= 1.6 * pitch && nadjacent < 64)
                adjacent_gap[nadjacent++] = gap;
            if (gap < gap_min || gap > gap_max) continue;
            long target_long = lround(inner->q + 1.0 - outer->q);
            if (target_long < INT32_MIN || target_long > INT32_MAX) continue;
            double residual = fabs(
                outer->q + (double)target_long - inner->q - 1.0);
            if (residual > WR_ORDER_PHASE_TOL) continue;
            wr_append_observation(
                observation, capacity, count, dropped,
                inner->component, outer->component, (int32_t)target_long,
                residual,
                seam_lattice && wr_is_axial_boundary_pair(
                    inner, outer, axial_min,
                    component_axial_min, component_axial_max)
                    ? WR_OBS_SEAM_ORDER : WR_OBS_ORDER);
        }
        /* Multi-pitch order: strands k>=2 wraps apart on the SAME ray carry
         * the only evidence binding islands whose junction offsets exceed
         * the 2.25-pitch adjacency window -- exactly where the lift drift
         * lives (measured 2026-08-31: junction dr of 2-8 pitches with NO
         * relation at all).  The wrap-separation campaign proved a scalar
         * 9.5 gauge manufactures counts at multi-pitch range, so k snaps
         * against the LOCAL pitch (median of this ray's known-adjacent
         * gaps) and abstains outside a 0.25-pitch residual.  These flow
         * through the ordinary relation build + forest + loop closure, so
         * conflicting cycles are repaired or reported, never chained
         * blindly. */
        if (WR_ORDER_MULTI && pitch > 1e-6) {
            if (nadjacent >= 2) {
                qsort(adjacent_gap, nadjacent, sizeof(double),
                      wr_compare_gap);
                local_pitch = adjacent_gap[nadjacent / 2];
            }
            for (size_t i = first; i < last; i++) {
                for (size_t j = i + 1; j < last; j++) {
                    const WrStrand *inner = &strand[i];
                    const WrStrand *outer = &strand[j];
                    double gap = outer->radius - inner->radius;
                    long wraps = 0;
                    double snap_residual = 0.0, residual = 0.0;
                    long target_long = 0;
                    if (gap > WR_ORDER_MAX_WRAPS * local_pitch) break;
                    if (inner->component == outer->component) continue;
                    if (gap <= gap_max) continue;   /* k=1 handled above */
                    wraps = lround(gap / local_pitch);
                    if (wraps < 2) continue;
                    snap_residual = fabs(gap - (double)wraps * local_pitch);
                    if (snap_residual > WR_ORDER_SNAP * local_pitch)
                        continue;
                    target_long = lround(inner->q + (double)wraps - outer->q);
                    if (target_long < INT32_MIN || target_long > INT32_MAX)
                        continue;
                    residual = fabs(outer->q + (double)target_long -
                                    inner->q - (double)wraps);
                    if (residual > WR_ORDER_PHASE_TOL) continue;
                    wr_append_observation(
                        observation, capacity, count, dropped,
                        inner->component, outer->component,
                        (int32_t)target_long, residual,
                        seam_lattice && wr_is_axial_boundary_pair(
                            inner, outer, axial_min,
                            component_axial_min, component_axial_max)
                            ? WR_OBS_SEAM_ORDER : WR_OBS_ORDER);
                }
            }
        }
        first = last;
    }
    return bins;
}

/* Boundary-complete cube meshes meet geometrically but are independently
 * triangulated, so opposite sides of a face need not contribute a vertex to
 * the exact same 1/256 phase bucket.  On the shifted axial lattice, match
 * adjacent phase buckets across the two original axial bins.  Each lower-half
 * strand votes for at most one nearest-ray candidate on either radial side;
 * this preserves adjacent-wrap semantics without admitting all O(n^2) pairs.
 * Exact-phase seam pairs remain handled by wr_collect_order_lattice. */
static void wr_collect_seam_order(
    const WrStrand *strand, size_t nstrand, double pitch, double axial_min,
    const double *component_axial_min,
    const double *component_axial_max,
    WrObservation *observation, size_t capacity,
    size_t *count, size_t *dropped)
{
    if (!(pitch > 1e-6) || !isfinite(pitch)) return;
    const double gap_min = 0.35 * pitch;
    const double gap_max = 2.25 * pitch;
    for (size_t axial_first = 0; axial_first < nstrand;) {
        size_t axial_last = axial_first + 1;
        size_t phase_first[WR_PHASE_BINS], phase_last[WR_PHASE_BINS];
        while (axial_last < nstrand &&
               strand[axial_last].axial_bin ==
                   strand[axial_first].axial_bin)
            axial_last++;
        for (int p = 0; p < WR_PHASE_BINS; p++) {
            phase_first[p] = SIZE_MAX;
            phase_last[p] = SIZE_MAX;
        }
        for (size_t i = axial_first; i < axial_last; i++) {
            int p = strand[i].phase_bin;
            if (phase_first[p] == SIZE_MAX) phase_first[p] = i;
            phase_last[p] = i + 1;
        }
        for (size_t i = axial_first; i < axial_last; i++) {
            const WrStrand *source = &strand[i];
            size_t best[2] = { SIZE_MAX, SIZE_MAX };
            double best_residual[2] = { DBL_MAX, DBL_MAX };
            double best_pitch_error[2] = { DBL_MAX, DBL_MAX };
            if (source->base_axial_bin != source->axial_bin) continue;
            int neighbor_phase[2] = {
                source->phase_bin == 0
                    ? WR_PHASE_BINS - 1 : source->phase_bin - 1,
                source->phase_bin + 1 == WR_PHASE_BINS
                    ? 0 : source->phase_bin + 1
            };
            for (int np = 0; np < 2; np++) {
                int p = neighbor_phase[np];
                if (phase_first[p] == SIZE_MAX) continue;
                for (size_t j = phase_first[p]; j < phase_last[p]; j++) {
                    const WrStrand *candidate = &strand[j];
                    if (candidate->base_axial_bin !=
                             source->base_axial_bin + 1 ||
                        candidate->component == source->component ||
                        !wr_is_axial_boundary_pair(
                            source, candidate, axial_min,
                            component_axial_min, component_axial_max))
                        continue;
                    const WrStrand *inner = source;
                    const WrStrand *outer = candidate;
                    if (inner->radius > outer->radius) {
                        inner = candidate;
                        outer = source;
                    }
                    double gap = outer->radius - inner->radius;
                    if (gap < gap_min || gap > gap_max) continue;
                    long target_long = lround(
                        inner->q + 1.0 - outer->q);
                    if (target_long < INT32_MIN ||
                        target_long > INT32_MAX)
                        continue;
                    double residual = fabs(
                        outer->q + (double)target_long - inner->q - 1.0);
                    if (residual > WR_ORDER_PHASE_TOL) continue;
                    int side = candidate->radius < source->radius ? 0 : 1;
                    double pitch_error = fabs(gap - pitch);
                    int take = residual < best_residual[side];
                    if (residual == best_residual[side]) {
                        take = pitch_error < best_pitch_error[side];
                        if (pitch_error == best_pitch_error[side] &&
                            best[side] != SIZE_MAX) {
                            const WrStrand *old = &strand[best[side]];
                            take = candidate->component < old->component ||
                                (candidate->component == old->component &&
                                 candidate->raw_turn < old->raw_turn);
                        }
                    }
                    if (take) {
                        best[side] = j;
                        best_residual[side] = residual;
                        best_pitch_error[side] = pitch_error;
                    }
                }
            }
            for (int side = 0; side < 2; side++) {
                if (best[side] == SIZE_MAX) continue;
                const WrStrand *inner = source;
                const WrStrand *outer = &strand[best[side]];
                if (inner->radius > outer->radius) {
                    const WrStrand *swap = inner; inner = outer; outer = swap;
                }
                long target_long = lround(inner->q + 1.0 - outer->q);
                double residual = fabs(
                    outer->q + (double)target_long - inner->q - 1.0);
                wr_append_observation(
                    observation, capacity, count, dropped,
                    inner->component, outer->component,
                    (int32_t)target_long, residual, WR_OBS_SEAM_ORDER);
            }
        }
        axial_first = axial_last;
    }
}

/* Collect order on two interleaved axial lattices.  A cube face which is a
 * boundary in the original 2-voxel lattice lies at the centre of the
 * half-bin-shifted lattice, so strands on opposite sides receive the same
 * ordinary radial-order test.  The two passes deliberately share every gap,
 * phase, eligibility, and mode-aggregation rule; only the partition changes.
 *
 * Reuse the raw strand buffer rather than retaining a second nsample-sized
 * copy: collect the shifted lattice, reconstruct the original lattice, then
 * leave its collapsed strands in place for the existing conflict detector. */
static int wr_collect_order(
    WrStrand *strand, size_t nsample,
    const int32_t *sample_vertex, size_t nvertices,
    const double *axial, const double *radius, const double *q,
    const int32_t *component, int32_t ncomponents,
    const double *component_axial_min,
    const double *component_axial_max,
    double axial_min, double pitch,
    WrObservation *observation, size_t capacity,
    size_t *count, size_t *dropped,
    size_t *out_nstrand, size_t *out_bins)
{
    if (wr_fill_strands(
            strand, nsample, sample_vertex, nvertices,
            axial, radius, q, component, ncomponents,
            axial_min, 0.5 * WR_AXIAL_BIN) != 0)
        return -1;
    size_t shifted_nstrand = wr_collapse_strands(strand, nsample);
    if (shifted_nstrand == 0) return -1;
    size_t bins = wr_collect_order_lattice(
        strand, shifted_nstrand, pitch, 1, axial_min,
        component_axial_min, component_axial_max, observation, capacity,
        count, dropped);
    wr_collect_seam_order(
        strand, shifted_nstrand, pitch, axial_min,
        component_axial_min, component_axial_max,
        observation, capacity, count, dropped);

    if (wr_fill_strands(
            strand, nsample, sample_vertex, nvertices,
            axial, radius, q, component, ncomponents,
            axial_min, 0.0) != 0)
        return -1;
    size_t nstrand = wr_collapse_strands(strand, nsample);
    if (nstrand == 0) return -1;
    bins += wr_collect_order_lattice(
        strand, nstrand, pitch, 0, axial_min,
        component_axial_min, component_axial_max,
        observation, capacity, count, dropped);
    *out_nstrand = nstrand;
    *out_bins = bins;
    return 0;
}

static int wr_build_relations(
    Arena_T arena, WrObservation *observation, size_t nobservation,
    WrRelation **out_relation, size_t *out_count)
{
    qsort(observation, nobservation, sizeof(WrObservation),
          wr_compare_observation);
    if (nobservation > ((size_t)LONG_MAX / sizeof(WrRelation)) - 1)
        return -1;
    WrRelation *relation = (WrRelation *)ARENA_ALLOC(
        arena, (long)((nobservation + 1) * sizeof(WrRelation)));
    size_t nrelation = 0;
    for (size_t first = 0; first < nobservation;) {
        size_t last = first + 1;
        while (last < nobservation &&
               observation[last].a == observation[first].a &&
               observation[last].b == observation[first].b)
            last++;
        WrCandidate continuation, order, seam_order;
        memset(&continuation, 0, sizeof(continuation));
        memset(&order, 0, sizeof(order));
        memset(&seam_order, 0, sizeof(seam_order));
        for (size_t kind_first = first; kind_first < last;) {
            size_t kind_last = kind_first + 1;
            while (kind_last < last &&
                   observation[kind_last].kind ==
                   observation[kind_first].kind)
                kind_last++;
            WrCandidate candidate = wr_observation_mode(
                observation, kind_first, kind_last,
                observation[kind_first].kind);
            if (candidate.kind == WR_OBS_CONTINUATION)
                continuation = candidate;
            else if (candidate.kind == WR_OBS_SEAM_ORDER)
                seam_order = candidate;
            else
                order = candidate;
            kind_first = kind_last;
        }
        WrCandidate chosen_order;
        if (seam_order.valid && seam_order.eligible)
            chosen_order = seam_order;
        else if (order.valid && order.eligible)
            chosen_order = order;
        else if (seam_order.valid &&
                 (!order.valid || seam_order.mode_observations >=
                                  order.mode_observations))
            chosen_order = seam_order;
        else
            chosen_order = order;
        WrCandidate chosen;
        if (continuation.valid && continuation.eligible)
            chosen = continuation;
        else if (chosen_order.valid && chosen_order.eligible)
            chosen = chosen_order;
        else if (continuation.valid &&
                 (!chosen_order.valid ||
                  continuation.mode_observations >=
                      chosen_order.mode_observations))
            chosen = continuation;
        else
            chosen = chosen_order;
        if (chosen.valid) {
            WrRelation *current = &relation[nrelation++];
            memset(current, 0, sizeof(*current));
            current->a = observation[first].a;
            current->b = observation[first].b;
            current->target = chosen.target;
            current->kind = chosen.kind;
            current->observations = chosen.observations;
            current->mode_observations = chosen.mode_observations;
            current->agreement = chosen.agreement;
            current->residual = chosen.residual;
            current->weight = chosen.weight;
            current->eligible = (uint8_t)chosen.eligible;
        }
        first = last;
    }
    qsort(relation, nrelation, sizeof(WrRelation),
          wr_compare_relation_weight);
    *out_relation = relation;
    *out_count = nrelation;
    return 0;
}

/* Integer loop closure (see the WR_REPAIR_* constants above).  Operates on
 * the BFS-chained corrections BEFORE anchoring/packing, per forest tree.
 * Hierarchy (relation weight) chose the forest; measured support
 * (mode_observations) votes on repairs -- with the continuation +1000 weight
 * offset, weight-ranked cuts would price order-kind closers at ~zero exactly
 * where the 10x drift lives (fit-stage conflicts are 100% order-kind). */
static int wr_repair_subtree_shifts(
    Arena_T arena, const WrRelation *relation, size_t nrelation,
    int32_t ncomponents,
    const int32_t *head, const WrAdjacency *adjacency,
    int32_t *correction, WindingRegisterStats *stats)
{
    int32_t *order = (int32_t *)ARENA_ALLOC(
        arena, (long)((size_t)ncomponents * sizeof(int32_t)));
    int32_t *pre = (int32_t *)ARENA_ALLOC(
        arena, (long)((size_t)ncomponents * sizeof(int32_t)));
    int32_t *sub = (int32_t *)ARENA_ALLOC(
        arena, (long)((size_t)ncomponents * sizeof(int32_t)));
    int32_t *parent = (int32_t *)ARENA_ALLOC(
        arena, (long)((size_t)ncomponents * sizeof(int32_t)));
    int32_t *tree_of = (int32_t *)ARENA_ALLOC(
        arena, (long)((size_t)ncomponents * sizeof(int32_t)));
    int32_t *stack = (int32_t *)ARENA_ALLOC(
        arena, (long)((size_t)ncomponents * sizeof(int32_t)));
    int32_t ntree = 0, npre = 0;
    uint64_t work = 0;
    size_t neligible = 0;
    for (size_t i = 0; i < nrelation; i++)
        if (relation[i].eligible &&
            (relation[i].kind != WR_OBS_SEAM_ORDER || relation[i].selected))
            neligible++;
    if (neligible == 0) return 0;

    /* preorder tour over the forest adjacency; trees are contiguous spans */
    for (int32_t c = 0; c < ncomponents; c++) { pre[c] = -1; parent[c] = -1; }
    int32_t *tree_first_pre = (int32_t *)ARENA_ALLOC(
        arena, (long)(((size_t)ncomponents + 1) * sizeof(int32_t)));
    for (int32_t seed = 0; seed < ncomponents; seed++) {
        size_t top = 0;
        if (pre[seed] >= 0) continue;
        tree_first_pre[ntree] = npre;
        stack[top++] = seed;
        pre[seed] = npre;
        order[npre++] = seed;
        tree_of[seed] = ntree;
        while (top > 0) {
            int32_t a = stack[--top];
            for (int32_t edge = head[a]; edge >= 0;
                 edge = adjacency[edge].next) {
                int32_t b = adjacency[edge].other;
                if (pre[b] >= 0) continue;
                parent[b] = a;
                pre[b] = npre;
                order[npre++] = b;
                tree_of[b] = ntree;
                stack[top++] = b;
            }
        }
        ntree++;
    }
    tree_first_pre[ntree] = npre;
    /* subtree sizes: reverse preorder accumulation */
    for (int32_t c = 0; c < ncomponents; c++) sub[c] = 1;
    for (int32_t i = npre - 1; i >= 0; i--) {
        int32_t c = order[i];
        if (parent[c] >= 0) sub[parent[c]] += sub[c];
    }

    /* eligible edges, grouped per tree (counting sort keeps the
     * weight-sorted relation order inside each group -- deterministic) */
    int32_t *e_a = (int32_t *)ARENA_ALLOC(
        arena, (long)(neligible * sizeof(int32_t)));
    int32_t *e_b = (int32_t *)ARENA_ALLOC(
        arena, (long)(neligible * sizeof(int32_t)));
    int32_t *e_target = (int32_t *)ARENA_ALLOC(
        arena, (long)(neligible * sizeof(int32_t)));
    int64_t *e_w = (int64_t *)ARENA_ALLOC(
        arena, (long)(neligible * sizeof(int64_t)));
    int64_t *e_implied = (int64_t *)ARENA_ALLOC(
        arena, (long)(neligible * sizeof(int64_t)));
    uint8_t *e_conflict = (uint8_t *)ARENA_ALLOC(
        arena, (long)(neligible * sizeof(uint8_t)));
    size_t *tree_edge_first = (size_t *)ARENA_CALLOC(
        arena, (long)ntree + 1, (long)sizeof(size_t));
    size_t *edge_of_tree = (size_t *)ARENA_ALLOC(
        arena, (long)(neligible * sizeof(size_t)));
    {
        size_t e = 0;
        for (size_t i = 0; i < nrelation; i++) {
            const WrRelation *r = &relation[i];
            /* A contradictory half-bin seam observation is a probabilistic
             * cycle closer, not a deterministic subtree-repair vote.  A seam
             * edge selected solely to connect two otherwise disconnected
             * conservative islands remains part of the forest and is priced
             * here like every other tree edge. */
            if (!r->eligible ||
                (r->kind == WR_OBS_SEAM_ORDER && !r->selected))
                continue;
            e_a[e] = r->a; e_b[e] = r->b; e_target[e] = r->target;
            e_w[e] = r->mode_observations > 0
                   ? (int64_t)r->mode_observations : 1;
            if (!r->selected) stats->repair_closers++;
            tree_edge_first[tree_of[r->a] + 1]++;
            e++;
        }
        for (int32_t t = 0; t < ntree; t++)
            tree_edge_first[t + 1] += tree_edge_first[t];
        size_t *cursor = (size_t *)ARENA_ALLOC(
            arena, (long)((size_t)ntree * sizeof(size_t)));
        for (int32_t t = 0; t < ntree; t++) cursor[t] = tree_edge_first[t];
        for (size_t k = 0; k < e; k++)
            edge_of_tree[cursor[tree_of[e_a[k]]]++] = k;
    }
    if (stats->repair_closers == 0) return 0;

    WrCut *cut = (WrCut *)ARENA_ALLOC(
        arena, (long)((size_t)ncomponents * sizeof(WrCut)));
    size_t *conflicted = (size_t *)ARENA_ALLOC(
        arena, (long)(neligible * sizeof(size_t)));
    size_t *crossing = (size_t *)ARENA_ALLOC(
        arena, (long)(neligible * sizeof(size_t)));
    uint8_t *crossing_inb = (uint8_t *)ARENA_ALLOC(
        arena, (long)(neligible * sizeof(uint8_t)));

    for (int32_t t = 0; t < ntree; t++) {
        size_t efirst = tree_edge_first[t], elast = tree_edge_first[t + 1];
        size_t tree_edges = elast - efirst;
        int32_t span_first = tree_first_pre[t];
        int32_t span_count = tree_first_pre[t + 1] - span_first;
        if (tree_edges == 0) continue;
        /* a tree's forest edges number span_count-1 exactly; anything above
         * that is a cycle-closing relation worth checking */
        if (tree_edges <= (size_t)(span_count - 1)) continue;
        for (int round = 0; round < WR_REPAIR_MAX_ROUNDS; round++) {
            size_t nconf = 0, ncand = 0;
            int64_t best_gain = 0;
            int32_t best_component = -1, best_shift = 0;
            for (size_t k = efirst; k < elast; k++) {
                size_t e = edge_of_tree[k];
                e_implied[e] = (int64_t)correction[e_b[e]] -
                               (int64_t)correction[e_a[e]];
                e_conflict[e] = e_implied[e] != (int64_t)e_target[e];
                if (e_conflict[e]) conflicted[nconf++] = e;
            }
            if (round == 0) stats->repair_conflicts_pre += nconf;
            if (nconf == 0) break;
            work += (uint64_t)span_count * (uint64_t)nconf +
                    (uint64_t)WR_REPAIR_RANK_CAP * (uint64_t)tree_edges;
            if (work > WR_REPAIR_WORK_CAP) {
                stats->repair_capped_roots++;
                break;
            }
            /* rank candidate cuts by crossing conflicted evidence mass */
            for (int32_t i = 0; i < span_count; i++) {
                int32_t c = order[span_first + i];
                int64_t weight = 0;
                if (parent[c] < 0) continue;
                for (size_t j = 0; j < nconf; j++) {
                    size_t e = conflicted[j];
                    int ina = pre[e_a[e]] >= pre[c] &&
                              pre[e_a[e]] < pre[c] + sub[c];
                    int inb = pre[e_b[e]] >= pre[c] &&
                              pre[e_b[e]] < pre[c] + sub[c];
                    if (ina != inb) weight += e_w[e];
                }
                if (weight > 0) {
                    cut[ncand].weight = weight;
                    cut[ncand].component = c;
                    ncand++;
                }
            }
            if (ncand == 0) break;
            qsort(cut, ncand, sizeof(WrCut), wr_compare_cut);
            if (ncand > (size_t)WR_REPAIR_RANK_CAP)
                ncand = (size_t)WR_REPAIR_RANK_CAP;
            for (size_t ci = 0; ci < ncand; ci++) {
                int32_t c = cut[ci].component;
                size_t ncross = 0;
                int64_t local_now = 0;
                static const int shifts[4] = { -2, -1, 1, 2 };
                if (cut[ci].weight <= best_gain) break;
                for (size_t k = efirst; k < elast; k++) {
                    size_t e = edge_of_tree[k];
                    int ina = pre[e_a[e]] >= pre[c] &&
                              pre[e_a[e]] < pre[c] + sub[c];
                    int inb = pre[e_b[e]] >= pre[c] &&
                              pre[e_b[e]] < pre[c] + sub[c];
                    if (ina == inb) continue;
                    crossing[ncross] = e;
                    crossing_inb[ncross] = (uint8_t)inb;
                    ncross++;
                    if (e_conflict[e]) local_now += e_w[e];
                }
                for (int si = 0; si < 4; si++) {
                    int s = shifts[si];
                    int64_t after = 0, gain = 0;
                    if (s > WR_REPAIR_MAX_SHIFT || s < -WR_REPAIR_MAX_SHIFT)
                        continue;
                    for (size_t j = 0; j < ncross; j++) {
                        size_t e = crossing[j];
                        int64_t implied = e_implied[e] +
                                          (crossing_inb[j] ? s : -s);
                        if (implied != (int64_t)e_target[e]) after += e_w[e];
                    }
                    gain = local_now - after;
                    if (gain > best_gain) {
                        best_gain = gain;
                        best_component = c;
                        best_shift = s;
                    }
                }
            }
            if (best_component < 0) break;
            for (int32_t i = pre[best_component];
                 i < pre[best_component] + sub[best_component]; i++) {
                int32_t c = order[i];
                int64_t value = (int64_t)correction[c] + best_shift;
                if (value < INT32_MIN || value > INT32_MAX) return -1;
                correction[c] = (int32_t)value;
            }
            stats->repair_shifts++;
        }
    }
    return 0;
}

/* The forest gives a globally coherent initialization but necessarily throws
 * away every cycle closer.  Refine it with all eligible relations as an MRF.
 * Labels are small residual shifts around the current solution, not absolute
 * turn numbers: five labels keep both memory and alpha-beta swap work bounded
 * on the 100k-component full-scroll regime.  Recentring after each accepted
 * round is the hierarchy; four rounds can repair an eight-turn accumulated
 * drift without ever constructing an enormous absolute-label table. */
static int wr_mrf_refine(
    Arena_T arena, const WrRelation *relation, size_t nrelation,
    int32_t ncomponents, UnionFind *graph, int32_t *correction,
    const WindingRegisterFieldUnary *field_unary,
    const WrForbiddenWinding *forbidden, size_t nforbidden,
    const WrForbiddenWinding *locked, size_t nlocked,
    const int32_t *trust_center, int32_t trust_radius, double trust_weight,
    float *out_confidence, WindingRegisterStats *stats)
{
    size_t nedge = 0;
    for (size_t i = 0; i < nrelation; i++)
        if (relation[i].eligible) nedge++;
    if (nedge == 0 || ncomponents < 2) return 0;
    if (trust_center != NULL &&
        (!(trust_weight > 0.0) || !isfinite(trust_weight)))
        return -1;
    Arena_Mark mark = Arena_save(arena);
    WindingMRFSite *site = (WindingMRFSite *)ARENA_ALLOC(
        arena, (size_t)ncomponents * sizeof *site);
    WindingMRFEdge *edge = (WindingMRFEdge *)ARENA_ALLOC(
        arena, nedge * sizeof *edge);
    size_t penalty_capacity = nforbidden;
    const size_t residual_labels = 2u * WR_MRF_DELTA_RADIUS + 1u;
    if (trust_center != NULL) {
        if (trust_radius < 0 ||
            (size_t)ncomponents >
                (SIZE_MAX - penalty_capacity) / residual_labels) {
            Arena_restore(arena, mark);
            return -1;
        }
        penalty_capacity += (size_t)ncomponents * residual_labels;
    }
    WindingMRFUnaryPenalty *penalty = penalty_capacity > 0
        ? (WindingMRFUnaryPenalty *)ARENA_ALLOC(
            arena, penalty_capacity * sizeof *penalty) : NULL;
    int32_t *lock_correction = (int32_t *)ARENA_ALLOC(
        arena, (size_t)ncomponents * sizeof *lock_correction);
    uint8_t *is_locked = (uint8_t *)ARENA_CALLOC(
        arena, (size_t)ncomponents, sizeof *is_locked);
    for (size_t i = 0; i < nlocked; i++) {
        int32_t c = locked[i].component;
        if (c < 0 || c >= ncomponents ||
            (is_locked[c] && lock_correction[c] != locked[i].correction)) {
            Arena_restore(arena, mark);
            return -1;
        }
        is_locked[c] = 1;
        lock_correction[c] = locked[i].correction;
    }
    WindingRegisterFieldUnary calibrated;
    memset(&calibrated, 0, sizeof calibrated);
    if (field_unary != NULL) {
        double *center = (double *)ARENA_ALLOC(
            arena, (size_t)ncomponents * sizeof *center);
        double *sigma = (double *)ARENA_ALLOC(
            arena, (size_t)ncomponents * sizeof *sigma);
        double *weight = (double *)ARENA_CALLOC(
            arena, (size_t)ncomponents, sizeof *weight);
        double *sw = (double *)ARENA_CALLOC(
            arena, (size_t)ncomponents, sizeof *sw);
        double *sx = (double *)ARENA_CALLOC(
            arena, (size_t)ncomponents, sizeof *sx);
        double *sy = (double *)ARENA_CALLOC(
            arena, (size_t)ncomponents, sizeof *sy);
        double *sxx = (double *)ARENA_CALLOC(
            arena, (size_t)ncomponents, sizeof *sxx);
        double *sxy = (double *)ARENA_CALLOC(
            arena, (size_t)ncomponents, sizeof *sxy);
        double *syy = (double *)ARENA_CALLOC(
            arena, (size_t)ncomponents, sizeof *syy);
        int32_t *count = (int32_t *)ARENA_CALLOC(
            arena, (size_t)ncomponents, sizeof *count);
        double *slope = (double *)ARENA_CALLOC(
            arena, (size_t)ncomponents, sizeof *slope);
        double *intercept = (double *)ARENA_CALLOC(
            arena, (size_t)ncomponents, sizeof *intercept);
        double *root_r2 = (double *)ARENA_CALLOC(
            arena, (size_t)ncomponents, sizeof *root_r2);
        double *root_rmse = (double *)ARENA_CALLOC(
            arena, (size_t)ncomponents, sizeof *root_rmse);
        for (int32_t c = 0; c < ncomponents; c++) {
            center[c] = NAN;
            sigma[c] = 1.0;
            if (!isfinite(field_unary->center[c]) ||
                !isfinite(field_unary->sigma[c]) ||
                !isfinite(field_unary->weight[c]) ||
                !(field_unary->sigma[c] > 0.0) ||
                !(field_unary->weight[c] > 0.0))
                continue;
            int32_t root = uf_find(graph, c);
            double w = field_unary->weight[c];
            double x = field_unary->center[c];
            double y = correction[c];
            sw[root] += w;
            sx[root] += w*x;
            sy[root] += w*y;
            sxx[root] += w*x*x;
            sxy[root] += w*x*y;
            syy[root] += w*y*y;
            count[root]++;
        }
        double r2_sum = 0.0;
        for (int32_t root = 0; root < ncomponents; root++) {
            if (count[root] < 3 || !(sw[root] > 0.0)) continue;
            double vx = sw[root]*sxx[root] - sx[root]*sx[root];
            double vy = sw[root]*syy[root] - sy[root]*sy[root];
            double cov = sw[root]*sxy[root] - sx[root]*sy[root];
            if (!(vx > 1e-8) || !(vy > 1e-8)) continue;
            slope[root] = cov / vx;
            root_r2[root] = cov*cov / (vx*vy);
            if (!(slope[root] > 0.125) || slope[root] > 8.0 ||
                root_r2[root] < 0.25 || !isfinite(root_r2[root])) {
                slope[root] = 0.0;
                continue;
            }
            intercept[root] = (sy[root] - slope[root]*sx[root]) / sw[root];
            stats->mrf_field_calibrated_roots++;
            r2_sum += root_r2[root];
        }
        for (int32_t c = 0; c < ncomponents; c++) {
            int32_t root = uf_find(graph, c);
            if (!(slope[root] > 0.0) ||
                !isfinite(field_unary->center[c]) ||
                !(field_unary->weight[c] > 0.0))
                continue;
            double predicted = intercept[root] +
                               slope[root]*field_unary->center[c];
            double residual = predicted - correction[c];
            root_rmse[root] += field_unary->weight[c] * residual * residual;
        }
        for (int32_t root = 0; root < ncomponents; root++)
            if (slope[root] > 0.0)
                root_rmse[root] = sqrt(root_rmse[root] / sw[root]);
        for (int32_t c = 0; c < ncomponents; c++) {
            int32_t root = uf_find(graph, c);
            if (!(slope[root] > 0.0) ||
                !isfinite(field_unary->center[c]) ||
                !isfinite(field_unary->sigma[c]) ||
                !(field_unary->weight[c] > 0.0))
                continue;
            center[c] = intercept[root] + slope[root]*field_unary->center[c];
            double measurement = slope[root]*field_unary->sigma[c];
            sigma[c] = sqrt(measurement*measurement +
                            root_rmse[root]*root_rmse[root]);
            if (sigma[c] < 1.0) sigma[c] = 1.0;
            /* The field is a gauge-scale/tie-break observation.  Eligible
             * continuation/order evidence remains the primary authority. */
            weight[c] = fmin(0.25, field_unary->weight[c] * root_r2[root]);
            stats->mrf_field_calibrated_sites++;
            if (wr_trace_conflicts())
                fprintf(stderr,
                        "  winding conflict trace: field c%d correction=%+d "
                        "center=%.9g residual=%+.9g sigma=%.9g weight=%.9g\n",
                        c, correction[c], center[c], center[c] - correction[c],
                        sigma[c], weight[c]);
        }
        if (stats->mrf_field_calibrated_roots > 0)
            stats->mrf_field_calibration_r2 = r2_sum /
                (double)stats->mrf_field_calibrated_roots;
        calibrated.center = center;
        calibrated.sigma = sigma;
        calibrated.weight = weight;
        field_unary = &calibrated;
    }
    WindingMRFOptions options;
    WindingMRF_default_options(&options);
    options.label_min = -WR_MRF_DELTA_RADIUS;
    options.label_max = WR_MRF_DELTA_RADIUS;
    options.max_iterations = 6;
    options.abstain_threshold = 0.35;
    int dump_solve = wr_dump_terms_dir() != NULL
                   ? wr_mrf_dump_serial++ : -1;

    for (int round = 0; round < WR_MRF_ROUNDS; round++) {
        for (int32_t c = 0; c < ncomponents; c++) {
            int have_field = field_unary != NULL &&
                isfinite(field_unary->center[c]) &&
                isfinite(field_unary->sigma[c]) &&
                isfinite(field_unary->weight[c]) &&
                field_unary->sigma[c] > 0.0 &&
                field_unary->weight[c] > 0.0;
            /* Keep the calibrated absolute observation during conflict
             * resolution.  The old trust-center branch replaced it with a
             * symmetric prior, so forbidding the current label left +/-1
             * exactly tied on weakly attached batches.  The trust term is an
             * additional sparse unary below; it must not erase the only
             * evidence which can choose the correct side of the old label. */
            site[c].center = have_field
                ? field_unary->center[c] - (double)correction[c] : 0.0;
            site[c].sigma = have_field ? field_unary->sigma[c] : 1.0;
            /* Break exact ties toward the existing hierarchy when the field
             * abstains; a valid field unary gets its measurement. */
            site[c].weight = have_field ? field_unary->weight[c] : 0.02;
            site[c].initial_label = 0;
            site[c].fixed = uf_find(graph, c) == c;
            if (is_locked[c]) {
                int64_t lock_delta = (int64_t)lock_correction[c] -
                                     (int64_t)correction[c];
                if (lock_delta < options.label_min ||
                    lock_delta > options.label_max) {
                    Arena_restore(arena, mark);
                    return -1;
                }
                site[c].initial_label = (int32_t)lock_delta;
                site[c].fixed = 1;
            }
        }
        size_t at = 0;
        for (size_t i = 0; i < nrelation; i++) {
            const WrRelation *current = &relation[i];
            if (!current->eligible) continue;
            int64_t solved = (int64_t)correction[current->b] -
                             correction[current->a];
            int64_t residual_target = (int64_t)current->target - solved;
            if (residual_target < INT32_MIN || residual_target > INT32_MAX) {
                Arena_restore(arena, mark);
                return -1;
            }
            edge[at].a = current->a;
            edge[at].b = current->b;
            edge[at].target = (int32_t)residual_target;
            /* The historical +1000 is a Kruskal hierarchy sentinel, not
             * evidence mass.  Remove it; loop votes otherwise use measured
             * support quality. */
            edge[at].weight = current->kind == WR_OBS_CONTINUATION
                ? WR_MRF_CONTINUATION_SCALE *
                    fmax(current->weight - 1000.0, 0.25)
                : fmax(current->weight, 0.25);
            at++;
        }
        size_t np = 0;
        for (size_t i = 0; i < nforbidden; i++) {
            int32_t c = forbidden[i].component;
            if (c < 0 || c >= ncomponents) {
                Arena_restore(arena, mark);
                return -1;
            }
            int64_t delta64 = (int64_t)forbidden[i].correction -
                              (int64_t)correction[c];
            if (delta64 < options.label_min || delta64 > options.label_max)
                continue;
            penalty[np].site = c;
            penalty[np].label = (int32_t)delta64;
            penalty[np].weight = WR_CONFLICT_UNARY_WEIGHT;
            np++;
        }
        if (trust_center != NULL) {
            for (int32_t c = 0; c < ncomponents; c++) {
                int64_t low = (int64_t)trust_center[c] - trust_radius;
                int64_t high = (int64_t)trust_center[c] + trust_radius;
                for (int32_t label = options.label_min;
                     label <= options.label_max; label++) {
                    int64_t candidate = (int64_t)correction[c] + label;
                    double value = trust_weight * wr_huber(
                        (double)(candidate - trust_center[c]),
                        options.unary_huber_delta);
                    if (candidate < low || candidate > high)
                        value += WR_CONFLICT_UNARY_WEIGHT;
                    if (!(value > 0.0)) continue;
                    penalty[np].site = c;
                    penalty[np].label = label;
                    penalty[np].weight = value;
                    np++;
                }
            }
        }
        int32_t *delta = NULL;
        float *confidence = NULL;
        WindingMRFStats mrf;
        if (WindingMRF_solve_with_penalties(
                arena, site, (size_t)ncomponents, edge, at,
                penalty, np, &options, &delta, &confidence, &mrf) != 0) {
            Arena_restore(arena, mark);
            return -1;
        }
        if (dump_solve >= 0)
            wr_dump_mrf_round(
                dump_solve, round, site, ncomponents,
                relation, nrelation, edge, at, penalty, np, &options,
                correction, delta, confidence, &mrf);
        if (stats->mrf_rounds == 0)
            stats->mrf_energy_before = mrf.energy_before;
        stats->mrf_energy_after = mrf.energy_after;
        stats->mrf_rounds++;
        stats->mrf_mean_confidence = mrf.mean_confidence;
        stats->mrf_abstained_sites = mrf.abstained_sites;
        if (out_confidence != NULL)
            memcpy(out_confidence, confidence,
                   (size_t)ncomponents * sizeof *out_confidence);
        size_t changed = 0;
        for (int32_t c = 0; c < ncomponents; c++) {
            if (delta[c] == 0) continue;
            int64_t value = (int64_t)correction[c] + delta[c];
            if (value < INT32_MIN || value > INT32_MAX) {
                Arena_restore(arena, mark);
                return -1;
            }
            correction[c] = (int32_t)value;
            changed++;
        }
        stats->mrf_label_changes += changed;
        if (changed == 0) break;
    }
    Arena_restore(arena, mark);
    return 0;
}

/* Detect mutually exclusive claims in the current winding cover.  The key is
 * exactly the discretization used to construct same-V order evidence:
 * (axial bin, absolute integer turn, wrapped phase bin).  Within one key,
 * close-radius fragments are one physical claimant.  Only clusters separated
 * by at least eight pitches (and at least 48 voxels) are called mutually
 * exclusive here: that is the unambiguous, spatially remote duplicate-UV
 * signature.  Claimants in different measured relation islands do not compete
 * here--unobserved islands are separated by the later deterministic atlas pack.
 *
 * Winner selection is deliberately global rather than greedy-by-radius.  A
 * component scores by the mass of every currently satisfied relation to the
 * rest of the solution, with a small size tie-break.  Thus the claimant that
 * agrees best with already coherent UVs keeps the assignment.  Losing
 * components must meet the configured UV-bin evidence threshold before their
 * current integer correction is added to the sparse forbidden list. */
static int wr_collect_conflict_exclusions(
    Arena_T arena, const WrStrand *strand, size_t nstrand,
    const WrRelation *relation, size_t nrelation,
    int32_t ncomponents, const int32_t *component_size,
    int32_t anchor_component, double pitch, UnionFind *graph,
    const int32_t *correction, const uint8_t *movable,
    WrForbiddenWinding *forbidden, size_t forbidden_capacity,
    size_t *nforbidden,
    WrForbiddenWinding *locked, size_t locked_capacity,
    size_t *nlocked, size_t *out_conflict_bins,
    size_t *out_losing_claims, size_t *out_new_exclusions,
    size_t *out_new_locks)
{
    if (arena == NULL || graph == NULL || correction == NULL ||
        nforbidden == NULL || out_conflict_bins == NULL ||
        nlocked == NULL || out_losing_claims == NULL ||
        out_new_exclusions == NULL || out_new_locks == NULL ||
        ncomponents <= 0 || component_size == NULL ||
        anchor_component < 0 || anchor_component >= ncomponents)
        return -1;
    *out_conflict_bins = 0;
    *out_losing_claims = 0;
    *out_new_exclusions = 0;
    *out_new_locks = 0;
    if (strand == NULL || nstrand < 2 || !(pitch > 0.0) || !isfinite(pitch))
        return 0;

    Arena_Mark mark = Arena_save(arena);
    WrConflictClaim *claim = (WrConflictClaim *)ARENA_ALLOC(
        arena, nstrand * sizeof *claim);
    double *score = (double *)ARENA_ALLOC(
        arena, (size_t)ncomponents * sizeof *score);
    uint32_t *losing_bins = (uint32_t *)ARENA_CALLOC(
        arena, (size_t)ncomponents, sizeof *losing_bins);
    uint32_t *winning_bins = (uint32_t *)ARENA_CALLOC(
        arena, (size_t)ncomponents, sizeof *winning_bins);
    for (int32_t c = 0; c < ncomponents; c++)
        score[c] = log1p((double)(component_size[c] > 0
                               ? component_size[c] : 0));
    for (size_t i = 0; i < nrelation; i++) {
        const WrRelation *current = &relation[i];
        if (!current->eligible) continue;
        int64_t solved = (int64_t)correction[current->b] -
                         (int64_t)correction[current->a];
        if (solved != current->target) continue;
        double mass = (double)(current->mode_observations > 0
                             ? current->mode_observations : 1) *
                      current->agreement /
                      (1.0 + 4.0 * current->residual);
        if (current->kind == WR_OBS_CONTINUATION) mass *= 4.0;
        score[current->a] += mass;
        score[current->b] += mass;
    }
    /* Every measured island needs one immutable gauge.  Making its union-find
     * root the strongest claimant prevents an exclusion from fighting the
     * fixed-site constraint inside the residual graph cut. */
    for (int32_t c = 0; c < ncomponents; c++)
        if (uf_find(graph, c) == c) score[c] += 1.0e12;
    score[anchor_component] += 1.0e13;

    for (size_t i = 0; i < nstrand; i++) {
        int32_t c = strand[i].component;
        if (c < 0 || c >= ncomponents) {
            Arena_restore(arena, mark);
            return -1;
        }
        int64_t absolute_turn = (int64_t)strand[i].raw_turn + correction[c];
        claim[i].axial_bin = strand[i].axial_bin;
        claim[i].phase_bin = strand[i].phase_bin;
        claim[i].component = c;
        claim[i].root = uf_find(graph, c);
        claim[i].absolute_turn = absolute_turn;
        claim[i].count = strand[i].count;
        claim[i].radius = strand[i].radius;
    }
    qsort(claim, nstrand, sizeof *claim, wr_compare_conflict_claim);

    double same_sheet_radius = WR_CONFLICT_RADIUS_FRACTION * pitch;
    if (same_sheet_radius < 3.0) same_sheet_radius = 3.0;
    double remote_radius = WR_CONFLICT_REMOTE_PITCHES * pitch;
    if (remote_radius < WR_CONFLICT_REMOTE_MIN_VOX)
        remote_radius = WR_CONFLICT_REMOTE_MIN_VOX;
    for (size_t first = 0; first < nstrand;) {
        size_t last = first + 1;
        while (last < nstrand &&
               claim[last].axial_bin == claim[first].axial_bin &&
               claim[last].absolute_turn == claim[first].absolute_turn &&
               claim[last].phase_bin == claim[first].phase_bin &&
               claim[last].root == claim[first].root)
            last++;
        if (last - first < 2) {
            first = last;
            continue;
        }

        size_t cluster_count = 0, winner_first = first, winner_last = first;
        double winner_score = -DBL_MAX;
        int winner_trusted = -1;
        for (size_t cluster_first = first; cluster_first < last;) {
            size_t cluster_last = cluster_first + 1;
            while (cluster_last < last &&
                   claim[cluster_last].radius -
                   claim[cluster_last - 1].radius <= same_sheet_radius)
                cluster_last++;
            double cluster_score = 0.0;
            int cluster_trusted = 0;
            for (size_t i = cluster_first; i < cluster_last; i++) {
                cluster_score += score[claim[i].component];
                if (movable != NULL && !movable[claim[i].component])
                    cluster_trusted = 1;
            }
            /* Once an outer episode has identified its losing batch, that
             * provenance is immutable.  A displaced loser may never become
             * the authority which evicts previously good context. */
            if (cluster_trusted > winner_trusted ||
                (cluster_trusted == winner_trusted &&
                 cluster_score > winner_score)) {
                winner_score = cluster_score;
                winner_trusted = cluster_trusted;
                winner_first = cluster_first;
                winner_last = cluster_last;
            }
            cluster_count++;
            cluster_first = cluster_last;
        }
        if (cluster_count >= WR_CONFLICT_MIN_CLAIMANTS) {
            double winner_min = claim[winner_first].radius;
            double winner_max = claim[winner_last - 1].radius;
            int bin_conflict = 0;
            for (size_t cluster_first = first;
                 cluster_first < last;) {
                size_t cluster_last = cluster_first + 1;
                while (cluster_last < last &&
                       claim[cluster_last].radius -
                       claim[cluster_last - 1].radius <= same_sheet_radius)
                    cluster_last++;
                if (cluster_first != winner_first) {
                    double cluster_min = claim[cluster_first].radius;
                    double cluster_max = claim[cluster_last - 1].radius;
                    double separation = cluster_max < winner_min
                        ? winner_min - cluster_max
                        : cluster_min > winner_max
                        ? cluster_min - winner_max : 0.0;
                    if (separation >= remote_radius) {
                        bin_conflict = 1;
                        for (size_t i = cluster_first;
                             i < cluster_last; i++) {
                            int32_t c = claim[i].component;
                            if (losing_bins[c] < UINT32_MAX)
                                losing_bins[c]++;
                            (*out_losing_claims)++;
                        }
                    }
                }
                cluster_first = cluster_last;
            }
            if (bin_conflict) {
                (*out_conflict_bins)++;
                for (size_t i = winner_first; i < winner_last; i++) {
                    int32_t c = claim[i].component;
                    if (winning_bins[c] < UINT32_MAX) winning_bins[c]++;
                }
            }
        }
        first = last;
    }

    if (forbidden != NULL) {
        for (int32_t c = 0; c < ncomponents; c++) {
            if (movable != NULL && !movable[c]) continue;
            if (losing_bins[c] < WR_CONFLICT_MIN_BINS ||
                winning_bins[c] != 0)
                continue;
            int already = 0;
            for (size_t i = 0; i < *nlocked; i++) {
                if (locked[i].component == c &&
                    locked[i].correction == correction[c]) {
                    already = 1;
                    break;
                }
            }
            for (size_t i = 0; i < *nforbidden; i++) {
                if (forbidden[i].component == c &&
                    forbidden[i].correction == correction[c]) {
                    already = 1;
                    break;
                }
            }
            if (already) continue;
            if (*nforbidden >= forbidden_capacity) {
                Arena_restore(arena, mark);
                return -1;
            }
            forbidden[*nforbidden].component = c;
            forbidden[*nforbidden].correction = correction[c];
            (*nforbidden)++;
            (*out_new_exclusions)++;
            if (wr_trace_conflicts())
                fprintf(stderr,
                        "  winding conflict trace: forbid c%d correction=%+d "
                        "losing_bins=%u winning_bins=%u score=%.9g root=%d "
                        "vertices=%d\n",
                        c, correction[c], (unsigned)losing_bins[c],
                        (unsigned)winning_bins[c], score[c],
                        uf_find(graph, c), component_size[c]);
        }
    }
    if (locked != NULL) {
        for (int32_t c = 0; c < ncomponents; c++) {
            if (winning_bins[c] < WR_CONFLICT_MIN_BINS ||
                losing_bins[c] != 0)
                continue;
            int already = 0;
            for (size_t i = 0; i < *nforbidden; i++) {
                if (forbidden[i].component == c &&
                    forbidden[i].correction == correction[c]) {
                    already = 1;
                    break;
                }
            }
            for (size_t i = 0; i < *nlocked; i++) {
                if (locked[i].component == c) {
                    already = 1;
                    break;
                }
            }
            if (already) continue;
            if (*nlocked >= locked_capacity) {
                Arena_restore(arena, mark);
                return -1;
            }
            locked[*nlocked].component = c;
            locked[*nlocked].correction = correction[c];
            (*nlocked)++;
            (*out_new_locks)++;
            if (wr_trace_conflicts())
                fprintf(stderr,
                        "  winding conflict trace: winner c%d correction=%+d "
                        "winning_bins=%u losing_bins=%u score=%.9g root=%d "
                        "vertices=%d\n",
                        c, correction[c], (unsigned)winning_bins[c],
                        (unsigned)losing_bins[c], score[c],
                        uf_find(graph, c), component_size[c]);
        }
    }
    Arena_restore(arena, mark);
    return 0;
}

static int wr_forest_corrections(
    Arena_T arena, WrRelation *relation, size_t nrelation,
    int32_t ncomponents, int32_t anchor_component,
    const int32_t *component_size, const int32_t *component,
    const double *q, size_t nvertices,
    const WrStrand *strand, size_t nstrand, double pitch,
    const WindingRegisterFieldUnary *field_unary,
    int enable_conflict_exclusion,
    int32_t **out_correction, int32_t **out_relation_island,
    float **out_component_confidence,
    UnionFind *out_graph,
    WindingRegisterStats *stats)
{
    UnionFind graph = UF_new(arena, ncomponents);
    size_t nforest = 0;
    /* Stage 1 is deliberately conservative.  Build the maximum-weight forest
     * from continuation and ordinary within-bin order first.  The shifted
     * lattice is then allowed to connect components which remain genuinely
     * disconnected, but a contradictory seam edge that closes an existing
     * path cannot choose an absolute gauge before the MRF and its unaries are
     * present. */
    for (int seam_pass = 0; seam_pass < 2; seam_pass++) {
        for (size_t i = 0; i < nrelation; i++) {
            WrRelation *current = &relation[i];
            if (!current->eligible) continue;
            if ((current->kind == WR_OBS_SEAM_ORDER) != seam_pass) continue;
            if (seam_pass == 0) stats->eligible_relations++;
            int32_t a = uf_find(&graph, current->a);
            int32_t b = uf_find(&graph, current->b);
            if (a == b) continue;
            uf_union(&graph, a, b);
            current->selected = 1;
            nforest++;
        }
    }
    /* Seam-pass relations are eligible even when they close a cycle. */
    for (size_t i = 0; i < nrelation; i++)
        if (relation[i].eligible && relation[i].kind == WR_OBS_SEAM_ORDER)
            stats->eligible_relations++;
    stats->forest_relations = nforest;

    int32_t *head = (int32_t *)ARENA_ALLOC(
        arena, (long)((size_t)ncomponents * sizeof(int32_t)));
    for (int32_t c = 0; c < ncomponents; c++) head[c] = -1;
    WrAdjacency *adjacency = (WrAdjacency *)ARENA_ALLOC(
        arena, (long)((2 * nforest + 1) * sizeof(WrAdjacency)));
    size_t nadjacency = 0;
    for (size_t i = 0; i < nrelation; i++) {
        WrRelation *current = &relation[i];
        if (!current->selected) continue;
        adjacency[nadjacency].other = current->b;
        adjacency[nadjacency].target = current->target;
        adjacency[nadjacency].next = head[current->a];
        head[current->a] = (int32_t)nadjacency++;
        adjacency[nadjacency].other = current->a;
        adjacency[nadjacency].target = -current->target;
        adjacency[nadjacency].next = head[current->b];
        head[current->b] = (int32_t)nadjacency++;
    }

    int32_t *correction = (int32_t *)ARENA_ALLOC(
        arena, (long)((size_t)ncomponents * sizeof(int32_t)));
    float *component_confidence = out_component_confidence != NULL
        ? (float *)ARENA_CALLOC(arena, (size_t)ncomponents, sizeof(float))
        : NULL;
    int32_t *queue = (int32_t *)ARENA_ALLOC(
        arena, (long)((size_t)ncomponents * sizeof(int32_t)));
    uint8_t *known = (uint8_t *)ARENA_CALLOC(
        arena, (long)ncomponents, (long)sizeof(uint8_t));
    for (int32_t root = 0; root < ncomponents; root++) {
        if (known[root]) continue;
        size_t head_index = 0, tail_index = 0;
        known[root] = 1;
        correction[root] = 0;
        queue[tail_index++] = root;
        while (head_index < tail_index) {
            int32_t a = queue[head_index++];
            for (int32_t edge = head[a]; edge >= 0;
                 edge = adjacency[edge].next) {
                int32_t b = adjacency[edge].other;
                int64_t value = (int64_t)correction[a] +
                                adjacency[edge].target;
                if (value < INT32_MIN || value > INT32_MAX) return -1;
                if (!known[b]) {
                    known[b] = 1;
                    correction[b] = (int32_t)value;
                    queue[tail_index++] = b;
                }
            }
        }
    }

    int32_t anchor_root = uf_find(&graph, anchor_component);
    {
        /* drift metric BEFORE closure: lifted-q span of the anchored island
         * (correction_min/max stay pack-inclusive and are NOT a drift
         * measure) */
        double span_min = DBL_MAX, span_max = -DBL_MAX;
        for (size_t i = 0; i < nvertices; i++) {
            int32_t c = component[i];
            double value = 0.0;
            if (uf_find(&graph, c) != anchor_root) continue;
            value = q[i] + (double)correction[c];
            if (value < span_min) span_min = value;
            if (value > span_max) span_max = value;
        }
        stats->anchor_span_pre_turns =
            span_max >= span_min ? span_max - span_min : 0.0;
    }
    if (WR_REPAIR_MAX_ROUNDS > 0 && nrelation > 0 && nforest > 0) {
        if (wr_repair_subtree_shifts(arena, relation, nrelation, ncomponents,
                                     head, adjacency, correction,
                                     stats) != 0)
            return -1;
    }
    /* Keep the conservative forest/loop-closure certificate as an absolute
     * Gaussian trust centre during every recentered MRF round.  Without an
     * absolute centre, the nominal fallback unary is centred at zero anew on
     * each round: a large connected batch can therefore take four cheap
     * +/-2 steps and drift eight turns while satisfying a collection of newly
     * exposed seam factors.  The trust centre still permits a locally
     * supported one-turn split (the intended Artifact C move), but charges a
     * growing cost for moving an ever larger part of the scroll or repeatedly
     * walking the same sites away from the first-stage certificate. */
    int32_t *mrf_trust_center = (int32_t *)ARENA_ALLOC(
        arena, (long)((size_t)ncomponents * sizeof *mrf_trust_center));
    memcpy(mrf_trust_center, correction,
           (size_t)ncomponents * sizeof *mrf_trust_center);
    if (wr_mrf_refine(arena, relation, nrelation, ncomponents, &graph,
                      correction, field_unary, NULL, 0,
                      NULL, 0,
                      mrf_trust_center, WR_CONFLICT_TRUST_RADIUS,
                      WR_MRF_BASE_PRIOR_WEIGHT,
                      component_confidence, stats) != 0)
        return -1;
    if (enable_conflict_exclusion &&
        strand != NULL && nstrand > 1 && pitch > 0.0) {
        if ((size_t)ncomponents > SIZE_MAX / WR_CONFLICT_MAX_ROUNDS)
            return -1;
        size_t forbidden_capacity =
            (size_t)ncomponents * WR_CONFLICT_MAX_ROUNDS;
        WrForbiddenWinding *forbidden = (WrForbiddenWinding *)ARENA_ALLOC(
            arena, forbidden_capacity * sizeof *forbidden);
        WrForbiddenWinding *round_winner = (WrForbiddenWinding *)ARENA_ALLOC(
            arena, (size_t)ncomponents * sizeof *round_winner);
        WrForbiddenWinding *round_locked = (WrForbiddenWinding *)ARENA_ALLOC(
            arena, (size_t)ncomponents * sizeof *round_locked);
        uint8_t *movable = (uint8_t *)ARENA_CALLOC(
            arena, (size_t)ncomponents, sizeof *movable);
        uint8_t *active = (uint8_t *)ARENA_ALLOC(
            arena, (size_t)ncomponents * sizeof *active);
        int32_t *before = (int32_t *)ARENA_ALLOC(
            arena, (size_t)ncomponents * sizeof *before);
        int32_t *best = (int32_t *)ARENA_ALLOC(
            arena, (size_t)ncomponents * sizeof *best);
        int32_t *trust_center = (int32_t *)ARENA_ALLOC(
            arena, (size_t)ncomponents * sizeof *trust_center);
        float *trial_confidence = component_confidence != NULL
            ? (float *)ARENA_ALLOC(
                arena, (size_t)ncomponents * sizeof *trial_confidence)
            : NULL;
        float *best_confidence = component_confidence != NULL
            ? (float *)ARENA_ALLOC(
                arena, (size_t)ncomponents * sizeof *best_confidence)
            : NULL;
        memcpy(trust_center, correction,
               (size_t)ncomponents * sizeof *trust_center);
        size_t nforbidden = 0;
        size_t initial_bins = 0;
        size_t best_bins = SIZE_MAX, best_losing = SIZE_MAX;
        size_t best_forbidden = 0, best_locked = 0;
        size_t attempted_rounds = 0, attempted_changes = 0;
        WindingRegisterStats best_stats = *stats;
        int best_valid = 0;
        int stopped = 0;
        int retry_pending = 0;
        for (int round = 0; round < WR_CONFLICT_MAX_ROUNDS; round++) {
            size_t bins = 0, losing = 0, added = 0, added_locks = 0;
            size_t nwinner = 0;
            int have_retry = retry_pending;
            retry_pending = 0;
            if (wr_collect_conflict_exclusions(
                    arena, strand, nstrand, relation, nrelation,
                    ncomponents, component_size, anchor_component, pitch,
                    &graph, correction, round == 0 ? NULL : movable,
                    forbidden, forbidden_capacity,
                    &nforbidden, round_winner, (size_t)ncomponents, &nwinner,
                    &bins, &losing, &added, &added_locks) != 0)
                return -1;
            /* The first detector pass establishes immutable provenance for
             * this correction episode.  A loser may be forbidden at several
             * successive labels, but it can never become a winner which
             * displaces an initially good component. */
            if (round == 0) {
                for (size_t i = 0; i < nforbidden; i++) {
                    int32_t c = forbidden[i].component;
                    if (c >= 0 && c < ncomponents) movable[c] = 1;
                }
            }
            memset(active, 0, (size_t)ncomponents * sizeof *active);
            for (size_t i = 0; i < nforbidden; i++) {
                int32_t c = forbidden[i].component;
                if (c >= 0 && c < ncomponents && movable[c] &&
                    forbidden[i].correction == correction[c])
                    active[c] = 1;
            }
            /* An explicit unambiguous winner remains context even if a stale
             * exclusion happens to mention its current gauge. */
            for (size_t i = 0; i < nwinner; i++)
                active[round_winner[i].component] = 0;
            size_t nround_locked = 0;
            for (int32_t c = 0; c < ncomponents; c++) {
                if (active[c]) continue;
                round_locked[nround_locked].component = c;
                round_locked[nround_locked].correction = correction[c];
                nround_locked++;
            }
            if (round == 0) {
                stats->mrf_conflict_bins_before = bins;
                initial_bins = bins;
                memcpy(best, correction,
                       (size_t)ncomponents * sizeof *best);
                if (component_confidence != NULL)
                    memcpy(best_confidence, component_confidence,
                           (size_t)ncomponents * sizeof *best_confidence);
                best_bins = bins;
                best_losing = losing;
                best_forbidden = 0;
                best_locked = 0;
                best_stats = *stats;
                best_valid = 1;
            }
            stats->mrf_conflict_bins_after = bins;
            stats->mrf_conflict_losing_claims = losing;
            stats->mrf_conflict_exclusions = nforbidden;
            stats->mrf_conflict_winner_locks = nround_locked;
            fprintf(stderr,
                    "  winding conflict exclusion: round=%d bins=%zu "
                    "losing_claims=%zu new_unaries=%zu total=%zu "
                    "explicit_winners=%zu context_locks=%zu\n",
                    round, bins, losing, added, nforbidden,
                    added_locks, nround_locked);
            if (bins == 0) {
                stats->mrf_conflict_converged = 1;
                stopped = 1;
                break;
            }
            if (added == 0 && !have_retry) {
                stats->mrf_conflict_converged = bins == 0;
                stopped = 1;
                break;
            }
            memcpy(before, correction,
                   (size_t)ncomponents * sizeof *before);
            WindingRegisterStats trial_stats = *stats;
            if (wr_mrf_refine(
                    arena, relation, nrelation, ncomponents, &graph,
                    correction, field_unary, forbidden, nforbidden,
                    round_locked, nround_locked,
                    trust_center, WR_CONFLICT_TRUST_RADIUS,
                    WR_CONFLICT_PRIOR_WEIGHT,
                    trial_confidence, &trial_stats) != 0)
                return -1;
            size_t changed = 0;
            for (int32_t c = 0; c < ncomponents; c++) {
                if (before[c] == correction[c]) continue;
                changed++;
                if (wr_trace_conflicts()) {
                    double before_energy = 0.0, after_energy = 0.0;
                    size_t before_ok = 0, after_ok = 0, degree = 0;
                    for (size_t i = 0; i < nrelation; i++) {
                        const WrRelation *current = &relation[i];
                        if (!current->eligible ||
                            (current->a != c && current->b != c))
                            continue;
                        int64_t old_solved = (int64_t)before[current->b] -
                                             before[current->a];
                        int64_t new_solved =
                            (int64_t)correction[current->b] -
                            correction[current->a];
                        double weight = current->kind == WR_OBS_CONTINUATION
                            ? 4.0 * fmax(current->weight - 1000.0, 0.25)
                            : fmax(current->weight, 0.25);
                        before_energy += weight * fabs(
                            (double)old_solved - current->target);
                        after_energy += weight * fabs(
                            (double)new_solved - current->target);
                        before_ok += old_solved == current->target;
                        after_ok += new_solved == current->target;
                        degree++;
                    }
                    fprintf(stderr,
                            "  winding conflict trace: move round=%d c%d "
                            "%+d->%+d active=%d relations=%zu sat=%zu->%zu "
                            "pair_energy=%.9g->%.9g\n",
                            round, c, before[c], correction[c],
                            active[c] != 0, degree, before_ok, after_ok,
                            before_energy, after_energy);
                    for (size_t i = 0; i < nrelation; i++) {
                        const WrRelation *current = &relation[i];
                        if (!current->eligible ||
                            (current->a != c && current->b != c))
                            continue;
                        int64_t old_solved = (int64_t)before[current->b] -
                                             before[current->a];
                        int64_t new_solved =
                            (int64_t)correction[current->b] -
                            correction[current->a];
                        fprintf(stderr,
                                "    relation c%d->c%d kind=%s target=%+d "
                                "solved=%+lld->%+lld support=%zu/%zu "
                                "agreement=%.6g residual=%.6g weight=%.9g%s\n",
                                current->a, current->b,
                                wr_observation_kind_name(current->kind),
                                current->target,
                                (long long)old_solved,
                                (long long)new_solved,
                                current->mode_observations,
                                current->observations, current->agreement,
                                current->residual, current->weight,
                                current->selected ? " forest" : "");
                    }
                }
            }
            attempted_rounds++;
            attempted_changes += changed;
            if (changed == 0) {
                stats->mrf_conflict_converged = 0;
                stopped = 1;
                break;
            }
            size_t after_bins = 0, after_losing = 0;
            size_t audit_added = 0, audit_locks = 0;
            size_t audit_count = nforbidden;
            size_t audit_locked = 0;
            if (wr_collect_conflict_exclusions(
                    arena, strand, nstrand, relation, nrelation,
                    ncomponents, component_size, anchor_component, pitch,
                    &graph, correction, movable, NULL, 0, &audit_count,
                    NULL, 0, &audit_locked,
                    &after_bins, &after_losing, &audit_added,
                    &audit_locks) != 0)
                return -1;
            trial_stats.mrf_conflict_bins_after = after_bins;
            trial_stats.mrf_conflict_losing_claims = after_losing;
            trial_stats.mrf_conflict_exclusions = nforbidden;
            trial_stats.mrf_conflict_winner_locks = nround_locked;
            fprintf(stderr,
                    "  winding conflict exclusion: solve round=%d "
                    "bins=%zu->%zu losing=%zu->%zu changed=%zu\n",
                    round, bins, after_bins, losing, after_losing, changed);
            /* A correction pass is a descent method over the detector which
             * requested it.  Accepting 2 -> 123 bins was the exact start of
             * the observed 4x5x5 cascade.  If graph cut cannot find a strict
             * lexicographic improvement, abstain and retain the last sound
             * certificate instead of manufacturing new conflicts. */
            if (after_bins > bins ||
                (after_bins == bins && after_losing >= losing)) {
                size_t rejected_added = 0;
                for (int32_t c = 0; c < ncomponents; c++) {
                    if (!active[c] || before[c] == correction[c]) continue;
                    int already = 0;
                    for (size_t i = 0; i < nforbidden; i++) {
                        if (forbidden[i].component == c &&
                            forbidden[i].correction == correction[c]) {
                            already = 1;
                            break;
                        }
                    }
                    if (already) continue;
                    if (nforbidden >= forbidden_capacity) return -1;
                    forbidden[nforbidden].component = c;
                    forbidden[nforbidden].correction = correction[c];
                    nforbidden++;
                    rejected_added++;
                }
                memcpy(correction, before,
                       (size_t)ncomponents * sizeof *correction);
                fprintf(stderr,
                        "  winding conflict exclusion: reject non-improving "
                        "iterate bins=%zu->%zu losing=%zu->%zu; "
                        "learned_unaries=%zu\n",
                        bins, after_bins, losing, after_losing,
                        rejected_added);
                /* Conflict-directed constraint generation: the graph cut
                 * supplied a coherent alternative, but the independent UV
                 * certificate proved it worse.  Forbid exactly that proposal
                 * for the movable batch and ask graph cut for its next-best
                 * solution.  This searches both sides without a hand-written
                 * winding-direction heuristic. */
                if (rejected_added > 0) {
                    retry_pending = 1;
                    continue;
                }
                stats->mrf_conflict_converged = 0;
                stopped = 1;
                break;
            }
            *stats = trial_stats;
            if (component_confidence != NULL)
                memcpy(component_confidence, trial_confidence,
                       (size_t)ncomponents * sizeof *component_confidence);
            if (after_bins < best_bins ||
                (after_bins == best_bins && after_losing < best_losing)) {
                memcpy(best, correction,
                       (size_t)ncomponents * sizeof *best);
                if (component_confidence != NULL)
                    memcpy(best_confidence, component_confidence,
                           (size_t)ncomponents * sizeof *best_confidence);
                best_bins = after_bins;
                best_losing = after_losing;
                best_forbidden = nforbidden;
                best_locked = nround_locked;
                best_stats = *stats;
            }
            if (after_bins == 0) {
                stats->mrf_conflict_converged = 1;
                stopped = 1;
                break;
            }
        }
        if (!stopped) {
            size_t bins = 0, losing = 0, added = 0, added_locks = 0;
            size_t audit_count = nforbidden;
            size_t audit_locked = 0;
            if (wr_collect_conflict_exclusions(
                    arena, strand, nstrand, relation, nrelation,
                    ncomponents, component_size, anchor_component, pitch,
                    &graph, correction, movable, NULL, 0, &audit_count,
                    NULL, 0, &audit_locked, &bins, &losing, &added,
                    &added_locks) != 0)
                return -1;
            stats->mrf_conflict_bins_after = bins;
            stats->mrf_conflict_losing_claims = losing;
            stats->mrf_conflict_converged = bins == 0;
        }
        if (best_valid) {
            size_t terminal_bins = stats->mrf_conflict_bins_after;
            memcpy(correction, best,
                   (size_t)ncomponents * sizeof *correction);
            *stats = best_stats;
            if (component_confidence != NULL)
                memcpy(component_confidence, best_confidence,
                       (size_t)ncomponents * sizeof *component_confidence);
            stats->mrf_conflict_rounds = attempted_rounds;
            stats->mrf_conflict_label_changes = attempted_changes;
            stats->mrf_conflict_bins_before = initial_bins;
            stats->mrf_conflict_bins_after = best_bins;
            stats->mrf_conflict_losing_claims = best_losing;
            stats->mrf_conflict_exclusions = best_forbidden;
            stats->mrf_conflict_winner_locks = best_locked;
            stats->mrf_conflict_converged = best_bins == 0;
            if (terminal_bins != best_bins)
                fprintf(stderr,
                        "  winding conflict exclusion: restore best "
                        "bounded iterate bins=%zu (terminal=%zu, "
                        "unaries=%zu locks=%zu)\n",
                        best_bins, terminal_bins,
                        best_forbidden, best_locked);
        }
    } else {
        stats->mrf_conflict_converged = 1;
    }
    int64_t anchor_shift = -(int64_t)correction[anchor_component];
    for (int32_t c = 0; c < ncomponents; c++) {
        if (uf_find(&graph, c) != anchor_root) continue;
        int64_t value = (int64_t)correction[c] + anchor_shift;
        if (value < INT32_MIN || value > INT32_MAX) return -1;
        correction[c] = (int32_t)value;
    }

    double *root_min = (double *)ARENA_ALLOC(
        arena, (long)((size_t)ncomponents * sizeof(double)));
    double *root_max = (double *)ARENA_ALLOC(
        arena, (long)((size_t)ncomponents * sizeof(double)));
    int64_t *root_vertices = (int64_t *)ARENA_CALLOC(
        arena, (long)ncomponents, (long)sizeof(int64_t));
    int32_t *root_head = (int32_t *)ARENA_ALLOC(
        arena, (long)((size_t)ncomponents * sizeof(int32_t)));
    int32_t *root_next = (int32_t *)ARENA_ALLOC(
        arena, (long)((size_t)ncomponents * sizeof(int32_t)));
    for (int32_t c = 0; c < ncomponents; c++) {
        root_min[c] = DBL_MAX;
        root_max[c] = -DBL_MAX;
        root_head[c] = -1;
    }
    for (size_t i = 0; i < nvertices; i++) {
        int32_t c = component[i], root = uf_find(&graph, c);
        double value = q[i] + (double)correction[c];
        if (value < root_min[root]) root_min[root] = value;
        if (value > root_max[root]) root_max[root] = value;
    }
    stats->anchor_span_turns =
        root_max[anchor_root] >= root_min[anchor_root]
            ? root_max[anchor_root] - root_min[anchor_root] : 0.0;
    for (int32_t c = 0; c < ncomponents; c++) {
        int32_t root = uf_find(&graph, c);
        root_vertices[root] += component_size[c];
        root_next[c] = root_head[root];
        root_head[root] = c;
    }

    WrPackGroup *pack = (WrPackGroup *)ARENA_ALLOC(
        arena, (long)((size_t)ncomponents * sizeof(WrPackGroup)));
    size_t npack = 0;
    for (int32_t c = 0; c < ncomponents; c++) {
        if (uf_find(&graph, c) != c || root_vertices[c] <= 0) continue;
        stats->relation_components++;
        if (c == anchor_root) continue;
        pack[npack].root = c;
        pack[npack].vertices = root_vertices[c];
        pack[npack].qmin = root_min[c];
        pack[npack].qmax = root_max[c];
        npack++;
    }
    stats->packed_relation_components = npack;
    qsort(pack, npack, sizeof(WrPackGroup), wr_compare_pack_group);

    /* Preserve which integer-gauge offsets are measured relative to one
     * another.  q still needs a deterministic serial cover layout so every
     * fragment has a unique identity, but downstream metric parameterization
     * must not mistake those artificial cover gaps for material arclength. */
    int32_t *relation_island = NULL;
    if (out_relation_island != NULL) {
        relation_island = (int32_t *)ARENA_ALLOC(
            arena, (long)((size_t)ncomponents * sizeof(int32_t)));
        for (int32_t c = 0; c < ncomponents; c++)
            relation_island[c] = uf_find(&graph, c) == anchor_root ? 0 : -1;
    }
    double cursor = root_max[anchor_root];
    for (size_t p = 0; p < npack; p++) {
        double needed = cursor + WR_PACK_GUTTER_TURNS - pack[p].qmin;
        int64_t shift = needed > 0.0
                      ? (int64_t)ceil(needed - 1e-12) : 0;
        if (shift < INT32_MIN || shift > INT32_MAX) return -1;
        /* Components were indexed by relation root once above.  Scanning all
         * ncomponents for every packed island made this O(C^2), which is
         * catastrophic for fragmented full-scroll inputs (106k components).
         * Every component now appears in exactly one root list, so all packed
         * updates together are O(C). */
        for (int32_t c = root_head[pack[p].root]; c >= 0;
             c = root_next[c]) {
            if (relation_island != NULL)
                relation_island[c] = (int32_t)p + 1;
            int64_t value = (int64_t)correction[c] + shift;
            if (value < INT32_MIN || value > INT32_MAX) return -1;
            correction[c] = (int32_t)value;
            stats->packed_mesh_components++;
        }
        double high = pack[p].qmax + (double)shift;
        if (high > cursor) cursor = high;
    }

    *out_correction = correction;
    if (out_relation_island != NULL)
        *out_relation_island = relation_island;
    if (out_component_confidence != NULL)
        *out_component_confidence = component_confidence;
    *out_graph = graph;
    return 0;
}

int WindingRegister_run(
    Arena_T arena,
    const float *vertices, size_t nvertices,
    const double *axial, const double *radius, const double *theta,
    const double *q,
    const int32_t *component, int32_t ncomponents,
    const int32_t *component_size, int32_t anchor_component,
    double axial_min, double pitch, int winding_sense,
    int32_t **out_correction, int32_t **out_relation_island,
    int32_t **out_continuation_island,
    WindingRegisterStats *stats)
{
    return WindingRegister_run_with_field(
        arena, vertices, nvertices, axial, radius, theta, q,
        component, ncomponents, component_size, anchor_component,
        axial_min, pitch, winding_sense, NULL, 1,
        out_correction, out_relation_island, out_continuation_island,
        NULL, stats);
}

int WindingRegister_run_with_field(
    Arena_T arena,
    const float *vertices, size_t nvertices,
    const double *axial, const double *radius, const double *theta,
    const double *q,
    const int32_t *component, int32_t ncomponents,
    const int32_t *component_size, int32_t anchor_component,
    double axial_min, double pitch, int winding_sense,
    const WindingRegisterFieldUnary *field_unary,
    int enable_conflict_exclusion,
    int32_t **out_correction, int32_t **out_relation_island,
    int32_t **out_continuation_island,
    float **out_component_confidence,
    WindingRegisterStats *stats)
{
    if (arena == NULL || vertices == NULL || nvertices == 0 ||
        axial == NULL || radius == NULL || theta == NULL || q == NULL ||
        component == NULL || ncomponents <= 0 || component_size == NULL ||
        anchor_component < 0 || anchor_component >= ncomponents ||
        (winding_sense != -1 && winding_sense != 1) ||
        out_correction == NULL || stats == NULL ||
        (field_unary != NULL &&
         (field_unary->center == NULL || field_unary->sigma == NULL ||
          field_unary->weight == NULL)))
        return -1;
    memset(stats, 0, sizeof(*stats));

    if ((size_t)ncomponents > (size_t)LONG_MAX / sizeof(int32_t))
        return -1;
    if (ncomponents == 1) {
        int32_t *correction = (int32_t *)ARENA_CALLOC(
            arena, 1, (long)sizeof(int32_t));
        *out_correction = correction;
        if (out_relation_island != NULL) {
            int32_t *relation_island = (int32_t *)ARENA_CALLOC(
                arena, 1, (long)sizeof(int32_t));
            *out_relation_island = relation_island;
        }
        if (out_continuation_island != NULL) {
            int32_t *continuation_island = (int32_t *)ARENA_CALLOC(
                arena, 1, (long)sizeof(int32_t));
            *out_continuation_island = continuation_island;
        }
        if (out_component_confidence != NULL) {
            float *confidence = (float *)ARENA_ALLOC(arena, sizeof(float));
            confidence[0] = 1.0f;
            *out_component_confidence = confidence;
        }
        stats->continuation_components = 1;
        stats->relation_components = 1;
        stats->continuation_satisfaction = 1.0;
        stats->order_satisfaction = 1.0;
        return 0;
    }

    WrStrand *strand = NULL;
    int32_t *sample_vertex = NULL;
    size_t nsample = 0;
    double *component_axial_min = NULL, *component_axial_max = NULL;
    if (wr_component_axial_extents(
            arena, axial, component, nvertices, ncomponents,
            &component_axial_min, &component_axial_max) != 0)
        return -1;
    if (wr_build_samples(
            arena, nvertices, axial, radius, q, component, ncomponents,
            axial_min, &strand, &sample_vertex, &nsample) != 0)
        return -1;
    if (nsample > (SIZE_MAX - 1024) / 4) return -1;
    size_t observation_capacity = 4 * nsample + 1024;
    if (observation_capacity >
        (size_t)LONG_MAX / sizeof(WrObservation))
        return -1;
    WrObservation *observation = (WrObservation *)ARENA_ALLOC(
        arena, (long)(observation_capacity * sizeof(WrObservation)));
    size_t nobservation = 0, dropped = 0;
    size_t order_reserve = 2 * nsample;
    size_t continuation_limit = observation_capacity > order_reserve
                              ? observation_capacity - order_reserve : 0;
    wr_collect_continuations(
        vertices, radius, theta, q, component, sample_vertex, nsample,
        pitch, winding_sense, observation, observation_capacity,
        continuation_limit, &nobservation, &dropped);

    size_t nstrand = 0;
    if (wr_collect_order(
            strand, nsample, sample_vertex, nvertices,
            axial, radius, q, component, ncomponents,
            component_axial_min, component_axial_max,
            axial_min, pitch, observation, observation_capacity,
            &nobservation, &dropped, &nstrand, &stats->bins) != 0)
        return -1;
    stats->strands = nstrand;
    stats->observations_dropped = dropped;
    for (size_t i = 0; i < nobservation; i++) {
        if (observation[i].kind == WR_OBS_CONTINUATION)
            stats->continuation_observations++;
        else
            stats->order_observations++;
    }

    WrRelation *relation = NULL;
    size_t nrelation = 0;
    if (wr_build_relations(
            arena, observation, nobservation, &relation, &nrelation) != 0)
        return -1;
    stats->relations = nrelation;

    /* Same-sheet continuation is an equivalence relation.  Once its closure
     * says two components share q, their radial adjacency cannot also mean
     * "one turn apart"; that is a coarse-bin duplicate, not a contradiction. */
    UnionFind continuation_graph = UF_new(arena, ncomponents);
    for (size_t i = 0; i < nrelation; i++) {
        WrRelation *current = &relation[i];
        if (current->eligible &&
            current->kind == WR_OBS_CONTINUATION)
            uf_union(&continuation_graph, current->a, current->b);
    }
    {
        int32_t *root_label = (int32_t *)ARENA_ALLOC(
            arena, (long)((size_t)ncomponents * sizeof(int32_t)));
        int32_t *labels = out_continuation_island != NULL
                        ? (int32_t *)ARENA_ALLOC(
                            arena, (long)((size_t)ncomponents * sizeof(int32_t)))
                        : NULL;
        for (int32_t c = 0; c < ncomponents; c++) root_label[c] = -1;
        int32_t anchor_root = uf_find(&continuation_graph, anchor_component);
        root_label[anchor_root] = 0;
        int32_t next_label = 1;
        for (int32_t c = 0; c < ncomponents; c++) {
            int32_t root = uf_find(&continuation_graph, c);
            if (root_label[root] < 0) root_label[root] = next_label++;
            if (labels != NULL) labels[c] = root_label[root];
        }
        stats->continuation_components = (size_t)next_label;
        if (out_continuation_island != NULL)
            *out_continuation_island = labels;
    }
    for (size_t i = 0; i < nrelation; i++) {
        WrRelation *current = &relation[i];
        if (!current->eligible || current->kind != WR_OBS_ORDER) continue;
        if (uf_find(&continuation_graph, current->a) !=
            uf_find(&continuation_graph, current->b))
            continue;
        current->eligible = 0;
        stats->order_relations_suppressed++;
    }

    int32_t *correction = NULL;
    UnionFind graph;
    if (wr_forest_corrections(
            arena, relation, nrelation, ncomponents, anchor_component,
             component_size, component, q, nvertices,
             strand, nstrand, pitch,
             field_unary, enable_conflict_exclusion,
             &correction, out_relation_island,
             out_component_confidence, &graph, stats) != 0)
        return -1;

    stats->correction_min = INT_MAX;
    stats->correction_max = INT_MIN;
    for (int32_t c = 0; c < ncomponents; c++) {
        if (correction[c] < stats->correction_min)
            stats->correction_min = correction[c];
        if (correction[c] > stats->correction_max)
            stats->correction_max = correction[c];
    }
    size_t continuation_total = 0, continuation_ok = 0;
    size_t order_total = 0, order_ok = 0;
    for (size_t i = 0; i < nobservation; i++) {
        const WrObservation *current = &observation[i];
        int64_t solved = (int64_t)correction[current->b] -
                         correction[current->a];
        int satisfied = solved == current->target;
        if (current->kind == WR_OBS_CONTINUATION) {
            continuation_total++;
            if (satisfied) continuation_ok++;
        } else {
            if (current->kind == WR_OBS_ORDER &&
                uf_find(&continuation_graph, current->a) ==
                uf_find(&continuation_graph, current->b)) {
                stats->order_observations_suppressed++;
                continue;
            }
            order_total++;
            if (satisfied) order_ok++;
        }
    }
    stats->continuation_satisfaction = continuation_total
        ? (double)continuation_ok / (double)continuation_total : 1.0;
    stats->order_satisfaction = order_total
        ? (double)order_ok / (double)order_total : 1.0;
    for (size_t i = 0; i < nrelation; i++) {
        const WrRelation *current = &relation[i];
        if (!current->eligible) continue;
        int64_t solved = (int64_t)correction[current->b] -
                         correction[current->a];
        if (solved != current->target) {
            stats->relation_conflicts++;
            fprintf(stderr,
                    "  winding conflict: c%d<->c%d kind=%s target=%d "
                    "solved=%lld support=%zu/%zu residual=%.3f%s\n",
                    current->a, current->b,
                    wr_observation_kind_name(current->kind),
                    current->target, (long long)solved,
                    current->mode_observations, current->observations,
                    current->residual,
            current->selected ? " (forest)" : "");
        }
    }
    wr_dump_final_terms(
        vertices, nvertices, q, component, ncomponents, component_size,
        anchor_component, correction, field_unary, relation, nrelation,
        &graph);

    if (ncomponents <= 64) {
        double *component_min = (double *)ARENA_ALLOC(
            arena, (long)((size_t)ncomponents * sizeof(double)));
        double *component_max = (double *)ARENA_ALLOC(
            arena, (long)((size_t)ncomponents * sizeof(double)));
        uint16_t *continuation_degree = (uint16_t *)ARENA_CALLOC(
            arena, (long)ncomponents, (long)sizeof(uint16_t));
        uint16_t *order_degree = (uint16_t *)ARENA_CALLOC(
            arena, (long)ncomponents, (long)sizeof(uint16_t));
        for (int32_t c = 0; c < ncomponents; c++) {
            component_min[c] = DBL_MAX;
            component_max[c] = -DBL_MAX;
        }
        for (size_t i = 0; i < nvertices; i++) {
            int32_t c = component[i];
            if (q[i] < component_min[c]) component_min[c] = q[i];
            if (q[i] > component_max[c]) component_max[c] = q[i];
        }
        for (size_t i = 0; i < nrelation; i++) {
            const WrRelation *current = &relation[i];
            if (!current->eligible) continue;
            uint16_t *degree = current->kind == WR_OBS_CONTINUATION
                             ? continuation_degree : order_degree;
            if (degree[current->a] < UINT16_MAX) degree[current->a]++;
            if (degree[current->b] < UINT16_MAX) degree[current->b]++;
        }
        int32_t anchor_root = uf_find(&graph, anchor_component);
        fprintf(stderr,
                "  winding components: id vertices raw_q corrected_q "
                "correction cont/order provenance\n");
        for (int32_t c = 0; c < ncomponents; c++) {
            const char *provenance =
                uf_find(&graph, c) != anchor_root ? "packed" :
                continuation_degree[c] ? "continuation" :
                order_degree[c] ? "ordered" :
                c == anchor_component ? "anchor" : "forest";
            fprintf(stderr,
                    "    c%-2d %8d [%7.3f,%7.3f] [%7.3f,%7.3f] "
                    "%+4d %u/%u %s\n",
                    c, component_size[c],
                    component_min[c], component_max[c],
                    component_min[c] + (double)correction[c],
                    component_max[c] + (double)correction[c],
                    correction[c], (unsigned)continuation_degree[c],
                    (unsigned)order_degree[c], provenance);
        }
    }

    *out_correction = correction;
    fprintf(stderr,
            "  winding gauge: sense=%+d samples=%zu bins=%zu strands=%zu "
            "obs(cont=%zu order=%zu suppressed=%zu) "
            "relations=%zu/%zu suppressed=%zu forest=%zu\n"
            "                 graph_islands=%zu packed=%zu comps=%zu "
            "corr=[%d,%d] sat(cont=%.3f order=%.3f) conflicts=%zu drop=%zu "
            "repair(closers=%zu pre=%zu shifts=%zu capped=%zu post=%zu) "
            "mrf(rounds=%zu changes=%zu abstain=%zu E=%.3g->%.3g conf=%.3f "
            "field=%zu/%zu r2=%.3f) "
            "exclude(rounds=%zu bins=%zu->%zu claims=%zu unaries=%zu locks=%zu "
            "changes=%zu converged=%s) "
            "span=%.1f->%.1ft\n",
            winding_sense, nsample, stats->bins, stats->strands,
            stats->continuation_observations, stats->order_observations,
            stats->order_observations_suppressed,
            stats->eligible_relations, stats->relations,
            stats->order_relations_suppressed,
            stats->forest_relations, stats->relation_components,
            stats->packed_relation_components, stats->packed_mesh_components,
            stats->correction_min, stats->correction_max,
            stats->continuation_satisfaction, stats->order_satisfaction,
            stats->relation_conflicts, stats->observations_dropped,
            stats->repair_closers, stats->repair_conflicts_pre,
            stats->repair_shifts, stats->repair_capped_roots,
            stats->relation_conflicts,
            stats->mrf_rounds, stats->mrf_label_changes,
            stats->mrf_abstained_sites, stats->mrf_energy_before,
            stats->mrf_energy_after, stats->mrf_mean_confidence,
            stats->mrf_field_calibrated_sites,
            stats->mrf_field_calibrated_roots,
            stats->mrf_field_calibration_r2,
            stats->mrf_conflict_rounds,
            stats->mrf_conflict_bins_before,
            stats->mrf_conflict_bins_after,
            stats->mrf_conflict_losing_claims,
            stats->mrf_conflict_exclusions,
            stats->mrf_conflict_winner_locks,
            stats->mrf_conflict_label_changes,
            stats->mrf_conflict_converged ? "yes" : "no",
            stats->anchor_span_pre_turns, stats->anchor_span_turns);
    return 0;
}

static void wr_selftest_check(int condition, const char *message, int *fails)
{
    if (condition) return;
    fprintf(stderr, "[winding register selftest] FAIL: %s\n", message);
    (*fails)++;
}

int WindingRegister_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();

    /* Same-sheet duplicates have different local integer gauges.  Proximity
     * continuation must align them without consulting absolute radius. */
    {
        enum { N = 24, NV = 2 * N };
        float vertex[NV * 3];
        double axial[NV], radius[NV], theta[NV], q[NV];
        int32_t component[NV], size[2] = {N, N};
        const double pitch = 9.5;
        for (int copy = 0; copy < 2; copy++) {
            for (int i = 0; i < N; i++) {
                size_t v = (size_t)copy * N + (size_t)i;
                double phase = 0.20 * (double)i / (double)(N - 1);
                double global_q = 3.0 + phase;
                double angle = WR_TWO_PI * phase;
                double rr = 50.0 + pitch * global_q +
                            4.0 * sin(3.0 * angle);
                vertex[v*3+0] = (float)(0.1 * (double)i);
                vertex[v*3+1] = (float)(rr * sin(angle));
                vertex[v*3+2] = (float)(rr * cos(angle));
                axial[v] = vertex[v*3+0];
                radius[v] = rr;
                theta[v] = angle;
                q[v] = global_q - (copy ? 3.0 : 0.0);
                component[v] = copy;
            }
        }
        int32_t *correction = NULL, *continuation_island = NULL;
        WindingRegisterStats stats;
        int rc = WindingRegister_run(
            arena, vertex, NV, axial, radius, theta, q, component, 2,
            size, 0, 0.0, pitch, 1, &correction, NULL,
            &continuation_island, &stats);
        wr_selftest_check(rc == 0, "continuation solve returns success", &fails);
        if (rc == 0) {
            wr_selftest_check(correction[0] == 0 && correction[1] == 3,
                              "continuation recovers integer gauge +3", &fails);
            wr_selftest_check(stats.relation_components == 1 &&
                              stats.continuation_components == 1 &&
                              continuation_island[0] == continuation_island[1] &&
                              stats.continuation_satisfaction > 0.99,
                              "continuation graph is connected/consistent",
                              &fails);
        }
    }

    /* Three non-touching layers share rays.  Radial order, not r/pitch as an
     * absolute coordinate, must produce the compact gauges 0,1,2. */
    {
        enum { NA = 8, NP = 8, NC = 3, NV = NA * NP * NC };
        float vertex[NV * 3];
        double axial[NV], radius[NV], theta[NV], q[NV];
        int32_t component[NV], size[NC] = {NA*NP, NA*NP, NA*NP};
        const double pitch = 9.5;
        size_t v = 0;
        for (int c = 0; c < NC; c++) {
            for (int a = 0; a < NA; a++) {
                for (int p = 0; p < NP; p++, v++) {
                    double phase = ((double)p + 0.25) / (double)NP;
                    double angle = WR_TWO_PI * phase;
                    double rr = 50.0 + pitch * ((double)c + phase);
                    vertex[v*3+0] = (float)(8.0 * (double)a + 0.25);
                    vertex[v*3+1] = (float)(rr * sin(angle));
                    vertex[v*3+2] = (float)(rr * cos(angle));
                    axial[v] = vertex[v*3+0];
                    radius[v] = rr;
                    theta[v] = angle;
                    q[v] = phase;
                    component[v] = c;
                }
            }
        }
        int32_t *correction = NULL, *continuation_island = NULL;
        WindingRegisterStats stats;
        int rc = WindingRegister_run(
            arena, vertex, NV, axial, radius, theta, q, component, NC,
            size, 0, 0.0, pitch, 1, &correction, NULL,
            &continuation_island, &stats);
        wr_selftest_check(rc == 0, "order solve returns success", &fails);
        if (rc == 0) {
            wr_selftest_check(correction[0] == 0 &&
                              correction[1] == 1 &&
                              correction[2] == 2,
                              "same-ray order recovers compact layer ranks",
                              &fails);
            wr_selftest_check(stats.order_satisfaction > 0.99 &&
                              stats.relation_components == 1 &&
                              stats.continuation_components == 3 &&
                              continuation_island[0] != continuation_island[1] &&
                              continuation_island[1] != continuation_island[2] &&
                              continuation_island[0] != continuation_island[2],
                              "order graph is connected/consistent", &fails);
        }
    }

    /* A cube face exactly on an original axial-bin boundary used to suppress
     * all order evidence across that face.  Four matching phase bins on the
     * two sides must remain separated in the original lattice, co-bin in the
     * half-shifted lattice, and aggregate into an eligible +1 relation. */
    {
        enum { NPHASE = 4, NV = 2 * NPHASE };
        const double axial_min = 4352.018066;
        const double face = 4608.018066;
        const double pitch = 9.5;
        WrStrand strand[NV];
        WrObservation observation[2 * NV];
        double axial[NV], radius[NV], q[NV];
        double component_axial_min[2] = {
            face - 0.243164, face + 0.091309
        };
        double component_axial_max[2] = {
            face - 0.243164, face + 0.091309
        };
        int32_t component[NV], sample_vertex[NV];
        size_t nobservation = 0, dropped = 0, nstrand = 0, bins = 0;
        for (int p = 0; p < NPHASE; p++) {
            double phase = 0.10 + 0.10 * (double)p;
            size_t inner = (size_t)p;
            size_t outer = (size_t)NPHASE + (size_t)p;
            axial[inner] = face - 0.243164;
            axial[outer] = face + 0.091309;
            radius[inner] = 100.0;
            radius[outer] = 106.1;
            q[inner] = q[outer] = phase;
            component[inner] = 0;
            component[outer] = 1;
            sample_vertex[inner] = (int32_t)inner;
            sample_vertex[outer] = (int32_t)outer;
        }
        int rc = wr_collect_order(
            strand, NV, sample_vertex, NV,
            axial, radius, q, component, 2,
            component_axial_min, component_axial_max,
            axial_min, pitch,
            observation, 2 * NV, &nobservation, &dropped,
            &nstrand, &bins);
        wr_selftest_check(
            rc == 0 && nobservation == NPHASE && dropped == 0,
            "half-shifted axial bins recover four cross-face votes", &fails);
        if (rc == 0) {
            Arena_Mark mark = Arena_save(arena);
            WrRelation *relation = NULL;
            size_t nrelation = 0;
            rc = wr_build_relations(
                arena, observation, nobservation, &relation, &nrelation);
            wr_selftest_check(
                rc == 0 && nrelation == 1 && relation[0].a == 0 &&
                relation[0].b == 1 && relation[0].target == 1 &&
                relation[0].kind == WR_OBS_SEAM_ORDER &&
                relation[0].eligible &&
                relation[0].mode_observations == NPHASE,
                "cross-face votes aggregate into eligible order +1", &fails);
            Arena_restore(arena, mark);
        }
        {
            WrStrand lower, upper;
            double overlap_min[2] = { face - 2.0, face - 2.0 };
            double overlap_max[2] = { face + 2.0, face + 2.0 };
            memset(&lower, 0, sizeof lower);
            memset(&upper, 0, sizeof upper);
            lower.base_axial_bin = 127;
            lower.component = 0;
            upper.base_axial_bin = 128;
            upper.component = 1;
            wr_selftest_check(
                !wr_is_axial_boundary_pair(
                    &lower, &upper, axial_min,
                    overlap_min, overlap_max),
                "interior-crossing components are not promoted as a seam",
                &fails);
            lower.base_axial_bin = 126;
            upper.base_axial_bin = 127;
            wr_selftest_check(
                !wr_is_axial_boundary_pair(
                    &lower, &upper, axial_min,
                    component_axial_min, component_axial_max),
                "non-cube axial boundaries are not promoted as a seam",
                &fails);
        }
    }

    /* Absolute radius is intentionally adversarial: sampling bias makes the
     * OUTER component's median radius smaller by ~50 vox, so an r/pitch anchor
     * would put it several turns inward.  Common-ray ordering is unchanged and
     * must still recover one layer. */
    {
        enum {
            NA = 8, NP = 8, NSHARED = NA * NP,
            NBIAS = 128, NPER = NSHARED + NBIAS, NV = 2 * NPER
        };
        float vertex[NV * 3];
        double axial[NV], radius[NV], theta[NV], q[NV];
        int32_t component[NV], size[2] = {NPER, NPER};
        const double pitch = 9.5;
        size_t v = 0;
        for (int c = 0; c < 2; c++) {
            for (int a = 0; a < NA; a++) {
                for (int p = 0; p < NP; p++, v++) {
                    double phase = ((double)p + 0.25) / (double)NP;
                    double angle = WR_TWO_PI * phase;
                    double z = 8.0 * (double)a + 0.25;
                    double deformation =
                        30.0 * sin(angle) + 12.0 * cos(0.31 * z);
                    double rr = 100.0 + pitch * (double)c + deformation;
                    vertex[v*3+0] = (float)z;
                    vertex[v*3+1] = (float)(rr * sin(angle));
                    vertex[v*3+2] = (float)(rr * cos(angle));
                    axial[v] = z; radius[v] = rr; theta[v] = angle;
                    q[v] = phase; component[v] = c;
                }
            }
            for (int i = 0; i < NBIAS; i++, v++) {
                double phase = c == 0 ? 0.25 : 0.75;
                double angle = WR_TWO_PI * phase;
                double z = 100.0 + 0.01 * (double)i;
                double rr = 100.0 + pitch * (double)c +
                            30.0 * sin(angle);
                vertex[v*3+0] = (float)z;
                vertex[v*3+1] = (float)(rr * sin(angle));
                vertex[v*3+2] = (float)(rr * cos(angle));
                axial[v] = z; radius[v] = rr; theta[v] = angle;
                q[v] = phase; component[v] = c;
            }
        }
        int32_t *correction = NULL;
        WindingRegisterStats stats;
        int rc = WindingRegister_run(
            arena, vertex, NV, axial, radius, theta, q, component, 2,
            size, 0, 0.0, pitch, 1, &correction, NULL, NULL, &stats);
        wr_selftest_check(radius[NSHARED] >
                          radius[NPER + NSHARED],
                          "adversarial medians reverse radial component order",
                          &fails);
        wr_selftest_check(rc == 0, "deformed order solve returns success",
                          &fails);
        if (rc == 0) {
            wr_selftest_check(correction[0] == 0 && correction[1] == 1,
                              "ray order survives non-Archimedean bias",
                              &fails);
            wr_selftest_check(stats.order_satisfaction > 0.99,
                              "deformed order observations stay consistent",
                              &fails);
        }
    }

    /* With no constraint path the second gauge is unknowable; it must be
     * reported and packed, not inferred from radius. */
    {
        float vertex[6] = {0,0,10, 1000,0,10};
        double axial[2] = {0,1000}, radius[2] = {10,10};
        double theta[2] = {0,0}, q[2] = {0,0};
        int32_t component[2] = {0,1}, size[2] = {1,1};
        int32_t *correction = NULL, *island = NULL;
        WindingRegisterStats stats;
        int rc = WindingRegister_run(
            arena, vertex, 2, axial, radius, theta, q, component, 2,
            size, 0, 0.0, 9.5, 1, &correction, &island, NULL, &stats);
        wr_selftest_check(rc == 0, "unobservable solve returns success", &fails);
        if (rc == 0) {
            wr_selftest_check(stats.relation_components == 2 &&
                              stats.packed_relation_components == 1 &&
                              correction[1] >= 1 && island != NULL &&
                              island[0] == 0 && island[1] == 1,
                              "unobservable island is declared and packed",
                              &fails);
        }
    }

    /* A fragmented scroll can contain tens of thousands of unobservable
     * relation islands.  Their deterministic packing must remain linear in
     * the component count and preserve the root ordering. */
    {
        enum { NC = 8192 };
        Arena_Mark mark = Arena_save(arena);
        int32_t *component = (int32_t *)ARENA_ALLOC(
            arena, (long)(NC * sizeof(int32_t)));
        int32_t *size = (int32_t *)ARENA_ALLOC(
            arena, (long)(NC * sizeof(int32_t)));
        double *q = (double *)ARENA_CALLOC(
            arena, NC, (long)sizeof(double));
        for (int32_t c = 0; c < NC; c++) {
            component[c] = c;
            size[c] = 1;
        }
        int32_t *correction = NULL, *island = NULL;
        UnionFind graph;
        WindingRegisterStats stats;
        memset(&stats, 0, sizeof(stats));
        int rc = wr_forest_corrections(
            arena, NULL, 0, NC, 0, size, component, q, NC,
            NULL, 0, 9.5,
            NULL, 1, &correction, &island, NULL, &graph, &stats);
        wr_selftest_check(rc == 0,
                          "many-island packing returns success", &fails);
        if (rc == 0) {
            wr_selftest_check(
                stats.relation_components == NC &&
                stats.packed_relation_components == NC - 1 &&
                stats.packed_mesh_components == NC - 1 &&
                correction[0] == 0 && correction[1] == 1 &&
                correction[NC - 1] == NC - 1 &&
                island[0] == 0 && island[1] == 1 &&
                island[NC - 1] == NC - 1,
                "many-island packing is compact and deterministic", &fails);
        }
        Arena_restore(arena, mark);
    }

    /* (R1) loop closure repairs a wrong heavy forest edge: 5-chain whose
     * r2 target is +2 (true +1) with two straddling closers.  The cut at
     * component 3 (crossing conflicted mass 16) shifts subtree {3,4} by -1
     * (gain 13: fixes both closers at mass 16, sacrifices r2 at mass 3). */
    {
        Arena_Mark mark = Arena_save(arena);
        int32_t component[5] = { 0, 1, 2, 3, 4 };
        int32_t size[5] = { 1, 1, 1, 1, 1 };
        double q[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
        WrRelation relation[6];
        int32_t *correction = NULL, *island = NULL;
        UnionFind graph;
        WindingRegisterStats stats;
        int rc = 0;
        memset(relation, 0, sizeof relation);
        relation[0].a = 0; relation[0].b = 1; relation[0].target = 1;
        relation[0].weight = 100.9; relation[0].mode_observations = 10;
        relation[1].a = 1; relation[1].b = 2; relation[1].target = 1;
        relation[1].weight = 100.8; relation[1].mode_observations = 10;
        relation[2].a = 2; relation[2].b = 3; relation[2].target = 2;
        relation[2].weight = 100.7; relation[2].mode_observations = 3;
        relation[3].a = 3; relation[3].b = 4; relation[3].target = 1;
        relation[3].weight = 100.6; relation[3].mode_observations = 10;
        relation[4].a = 1; relation[4].b = 3; relation[4].target = 2;
        relation[4].weight = 1.5; relation[4].mode_observations = 8;
        relation[5].a = 2; relation[5].b = 4; relation[5].target = 2;
        relation[5].weight = 1.4; relation[5].mode_observations = 8;
        for (int i = 0; i < 6; i++) relation[i].eligible = 1;
        memset(&stats, 0, sizeof stats);
        rc = wr_forest_corrections(
            arena, relation, 6, 5, 0, size, component, q, 5,
            NULL, 0, 9.5,
            NULL, 1, &correction, &island, NULL, &graph, &stats);
        wr_selftest_check(rc == 0, "(R1) closure returns success", &fails);
        if (rc == 0) {
            wr_selftest_check(
                correction[0] == 0 && correction[1] == 1 &&
                correction[2] == 2 && correction[3] == 3 &&
                correction[4] == 4,
                "(R1) subtree shift restores the true chain", &fails);
            wr_selftest_check(
                stats.forest_relations == 4 &&
                stats.repair_closers == 2 &&
                stats.repair_conflicts_pre == 2 &&
                stats.repair_shifts == 1 &&
                stats.repair_capped_roots == 0,
                "(R1) repair stats", &fails);
            wr_selftest_check(
                correction[3] - correction[1] == relation[4].target &&
                correction[4] - correction[2] == relation[5].target &&
                correction[3] - correction[2] != relation[2].target,
                "(R1) closers satisfied, wrong forest edge is the residual",
                &fails);
            wr_selftest_check(
                fabs(stats.anchor_span_pre_turns - 5.0) < 1e-9 &&
                fabs(stats.anchor_span_turns - 4.0) < 1e-9,
                "(R1) anchored-island span shrinks 5 -> 4 turns", &fails);
        }
        Arena_restore(arena, mark);
    }

    /* (R2) balanced evidence must NOT repair: the only cut with any effect
     * fixes mass 3 and breaks mass 3 (gain exactly 0); strict > 0 refuses. */
    {
        Arena_Mark mark = Arena_save(arena);
        int32_t component[4] = { 0, 1, 2, 3 };
        int32_t size[4] = { 1, 1, 1, 1 };
        double q[4] = { 0.0, 0.0, 0.0, 0.0 };
        WrRelation relation[4];
        int32_t *correction = NULL, *island = NULL;
        UnionFind graph;
        WindingRegisterStats stats;
        int rc = 0;
        memset(relation, 0, sizeof relation);
        relation[0].a = 0; relation[0].b = 1; relation[0].target = 1;
        relation[0].weight = 100.9; relation[0].mode_observations = 10;
        relation[1].a = 1; relation[1].b = 2; relation[1].target = 2;
        relation[1].weight = 100.8; relation[1].mode_observations = 3;
        relation[2].a = 2; relation[2].b = 3; relation[2].target = 1;
        relation[2].weight = 100.7; relation[2].mode_observations = 10;
        relation[3].a = 0; relation[3].b = 3; relation[3].target = 3;
        relation[3].weight = 1.5; relation[3].mode_observations = 3;
        for (int i = 0; i < 4; i++) relation[i].eligible = 1;
        memset(&stats, 0, sizeof stats);
        rc = wr_forest_corrections(
            arena, relation, 4, 4, 0, size, component, q, 4,
            NULL, 0, 9.5,
            NULL, 1, &correction, &island, NULL, &graph, &stats);
        wr_selftest_check(rc == 0, "(R2) closure returns success", &fails);
        if (rc == 0) {
            wr_selftest_check(
                correction[0] == 0 && correction[1] == 1 &&
                correction[2] == 3 && correction[3] == 4,
                "(R2) balanced evidence leaves corrections unchanged",
                &fails);
            wr_selftest_check(
                stats.repair_closers == 1 &&
                stats.repair_conflicts_pre == 1 &&
                stats.repair_shifts == 0,
                "(R2) strict-gain gate refuses the zero-gain shift", &fails);
        }
        Arena_restore(arena, mark);
    }

    /* (R3) determinism: two fresh runs of (R1) agree bit-for-bit. */
    {
        Arena_Mark mark = Arena_save(arena);
        int32_t component[5] = { 0, 1, 2, 3, 4 };
        int32_t size[5] = { 1, 1, 1, 1, 1 };
        double q[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
        int32_t run_correction[2][5];
        size_t run_shifts[2] = { 0, 0 };
        int ok = 1;
        for (int run = 0; run < 2; run++) {
            WrRelation relation[6];
            int32_t *correction = NULL, *island = NULL;
            UnionFind graph;
            WindingRegisterStats stats;
            memset(relation, 0, sizeof relation);
            relation[0].a = 0; relation[0].b = 1; relation[0].target = 1;
            relation[0].weight = 100.9; relation[0].mode_observations = 10;
            relation[1].a = 1; relation[1].b = 2; relation[1].target = 1;
            relation[1].weight = 100.8; relation[1].mode_observations = 10;
            relation[2].a = 2; relation[2].b = 3; relation[2].target = 2;
            relation[2].weight = 100.7; relation[2].mode_observations = 3;
            relation[3].a = 3; relation[3].b = 4; relation[3].target = 1;
            relation[3].weight = 100.6; relation[3].mode_observations = 10;
            relation[4].a = 1; relation[4].b = 3; relation[4].target = 2;
            relation[4].weight = 1.5; relation[4].mode_observations = 8;
            relation[5].a = 2; relation[5].b = 4; relation[5].target = 2;
            relation[5].weight = 1.4; relation[5].mode_observations = 8;
            for (int i = 0; i < 6; i++) relation[i].eligible = 1;
            memset(&stats, 0, sizeof stats);
            if (wr_forest_corrections(
                    arena, relation, 6, 5, 0, size, component, q, 5,
                    NULL, 0, 9.5,
                    NULL, 1, &correction, &island, NULL, &graph, &stats) != 0) {
                ok = 0;
                break;
            }
            memcpy(run_correction[run], correction, 5 * sizeof(int32_t));
            run_shifts[run] = stats.repair_shifts;
        }
        wr_selftest_check(
            ok &&
            memcmp(run_correction[0], run_correction[1],
                   5 * sizeof(int32_t)) == 0 &&
            run_shifts[0] == run_shifts[1],
            "(R3) closure is deterministic across runs", &fails);
        Arena_restore(arena, mark);
    }

    /* (R4) conflict exclusion: false continuations put four remotely distant
     * components in the same four UV bins.  The anchored claimant keeps the
     * assignment; sparse near-hard unaries force the other components to a
     * neighbouring integer label, after which catastrophic occupancy is gone. */
    {
        Arena_Mark mark = Arena_save(arena);
        int32_t component[4] = { 0, 1, 2, 3 };
        int32_t size[4] = { 100, 20, 20, 20 };
        double q[4] = { 0.10, 1.10, 2.10, 3.10 };
        WrRelation relation[3];
        WrStrand strand[16];
        int32_t *correction = NULL, *island = NULL;
        UnionFind graph;
        WindingRegisterStats stats;
        memset(relation, 0, sizeof relation);
        for (int c = 1; c < 4; c++) {
            relation[c - 1].a = 0;
            relation[c - 1].b = c;
            relation[c - 1].target = -c;
            relation[c - 1].weight = 1002.0 - c;
            relation[c - 1].observations = 8;
            relation[c - 1].mode_observations = 8;
            relation[c - 1].agreement = 1.0;
            relation[c - 1].eligible = 1;
            relation[c - 1].kind = WR_OBS_CONTINUATION;
        }
        for (int z = 0; z < 4; z++) {
            for (int c = 0; c < 4; c++) {
                WrStrand *s = &strand[4*z + c];
                memset(s, 0, sizeof *s);
                s->axial_bin = z;
                s->phase_bin = 25;
                s->component = c;
                s->raw_turn = c;
                s->count = 1;
                s->radius = 10.0 + 90.0*c;
                s->q = 0.10 + c;
            }
        }
        memset(&stats, 0, sizeof stats);
        {
            int32_t *initial = NULL, *initial_island = NULL;
            UnionFind initial_graph;
            WindingRegisterStats initial_stats;
            memset(&initial_stats, 0, sizeof initial_stats);
            int initial_rc = wr_forest_corrections(
                arena, relation, 3, 4, 0, size, component, q, 4,
                strand, 16, 9.5, NULL, 0, &initial, &initial_island, NULL,
                &initial_graph, &initial_stats);
            wr_selftest_check(
                initial_rc == 0 && initial[0] == 0 && initial[1] == -1 &&
                initial[2] == -2 && initial[3] == -3 &&
                initial_stats.mrf_conflict_rounds == 0 &&
                initial_stats.mrf_conflict_exclusions == 0,
                "(R4) initial-certificate mode stops before exclusion",
                &fails);
            for (int i = 0; i < 3; i++) relation[i].selected = 0;
        }
        int rc = wr_forest_corrections(
            arena, relation, 3, 4, 0, size, component, q, 4,
            strand, 16, 9.5, NULL, 1, &correction, &island, NULL,
            &graph, &stats);
        wr_selftest_check(rc == 0,
                          "(R4) conflict exclusion returns success", &fails);
        if (rc == 0) {
            wr_selftest_check(
                correction[0] == 0 && correction[1] != -1 &&
                correction[2] != -2 && correction[3] != -3,
                "(R4) losers leave their forbidden winding", &fails);
            wr_selftest_check(
                stats.mrf_conflict_bins_before == 4 &&
                stats.mrf_conflict_bins_after == 0 &&
                stats.mrf_conflict_exclusions == 3 &&
                stats.mrf_conflict_label_changes >= 3 &&
                stats.mrf_conflict_converged,
                "(R4) exclusion converges and reports occupancy closure",
                &fails);
        }
        Arena_restore(arena, mark);
    }

    /* (R5) monotone constraint generation: the first equally plausible move
     * of the three losing claims lands on a different, previously good set of
     * three claims.  That proposal must be rejected, the exact trial labels
     * must be learned as additional forbidden unaries, and a later proposal
     * must resolve the original collision without moving the good context. */
    {
        Arena_Mark mark = Arena_save(arena);
        int32_t component[7] = { 0, 1, 2, 3, 4, 5, 6 };
        int32_t size[7] = { 100, 20, 20, 20, 30, 30, 30 };
        double q[7] = { 0.10, 1.10, 2.10, 3.10, 0.20, 1.20, 2.20 };
        WrRelation relation[6];
        WrStrand strand[28];
        int32_t *correction = NULL, *island = NULL;
        UnionFind graph;
        WindingRegisterStats stats;
        memset(relation, 0, sizeof relation);
        for (int c = 1; c < 7; c++) {
            relation[c - 1].a = 0;
            relation[c - 1].b = c;
            relation[c - 1].target = -(c <= 3 ? c : c - 3);
            relation[c - 1].weight = 1002.0 - c;
            relation[c - 1].observations = 8;
            relation[c - 1].mode_observations = 8;
            relation[c - 1].agreement = 1.0;
            relation[c - 1].eligible = 1;
            relation[c - 1].kind = WR_OBS_CONTINUATION;
        }
        for (int z = 0; z < 4; z++) {
            for (int c = 0; c < 7; c++) {
                WrStrand *s = &strand[7*z + c];
                memset(s, 0, sizeof *s);
                s->axial_bin = z;
                s->phase_bin = 25;
                s->component = c;
                s->raw_turn = c <= 3 ? c : c - 4;
                s->count = 1;
                s->radius = 10.0 + 90.0*c;
                s->q = q[c];
            }
        }
        memset(&stats, 0, sizeof stats);
        int rc = wr_forest_corrections(
            arena, relation, 6, 7, 0, size, component, q, 7,
            strand, 28, 9.5, NULL, 1, &correction, &island, NULL,
            &graph, &stats);
        wr_selftest_check(rc == 0,
                          "(R5) guarded retry returns success", &fails);
        if (rc == 0) {
            wr_selftest_check(
                correction[0] == 0 &&
                correction[4] == -1 && correction[5] == -2 &&
                correction[6] == -3,
                "(R5) conflicting retry never displaces good context",
                &fails);
            wr_selftest_check(
                correction[1] != -1 && correction[2] != -2 &&
                correction[3] != -3 &&
                stats.mrf_conflict_bins_before == 4 &&
                stats.mrf_conflict_bins_after == 0 &&
                stats.mrf_conflict_rounds >= 2 &&
                stats.mrf_conflict_exclusions >= 6 &&
                stats.mrf_conflict_converged,
                "(R5) rejected cascade is learned and alternate converges",
                &fails);
        }
        Arena_restore(arena, mark);
    }

    Arena_dispose(&arena);
    if (fails == 0)
        fprintf(stderr, "[winding register selftest] all tests OK\n");
    return fails;
}
