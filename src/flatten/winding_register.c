#define _USE_MATH_DEFINES
#include "winding_register.h"

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
#define WR_PHASE_BINS         64
#define WR_AXIAL_BIN          8.0
/* REGIME-BOUND (changelog rule): 2M was tuned when the biggest input was the
 * welded 4x5x5.  The 10x pre-weld soup is 12.4M vertices; at stride 7 a
 * continuation PAIR only survives sampling ~2% of the time (1/49), the
 * relation graph starves (obs/comp 666 -> 37), almost nothing merges
 * (materials 8331 -> 5211 instead of ~200:1), and the atlas shreds into
 * confetti.  16M keeps every rung to date at stride 1. */
#define WR_MAX_SAMPLES        16000000u
#define WR_CELL_SCAN_LIMIT    96
#define WR_CONT_WTOL          0.30
#define WR_CONT_PHASE_TOL     0.25
#define WR_ORDER_PHASE_TOL    0.25
#define WR_PACK_GUTTER_TURNS  0.10

typedef struct {
    int32_t axial_bin, phase_bin, component, raw_turn;
    uint32_t count;
    double radius, q;
} WrStrand;

enum { WR_OBS_CONTINUATION = 0, WR_OBS_ORDER = 1 };

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

static double wr_wrap_to_pi(double angle)
{
    double value = fmod(angle + M_PI, WR_TWO_PI);
    if (value < 0.0) value += WR_TWO_PI;
    return value - M_PI;
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
    return a->raw_turn < b->raw_turn ? -1 :
           (a->raw_turn > b->raw_turn ? 1 : 0);
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
        candidate.weight =
            log1p((double)best_count) * candidate.agreement /
            (1.0 + 4.0 * candidate.residual);
    }
    return candidate;
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
        sampled[c] = 1;
        double qi = q[i];
        double phase = qi - floor(qi);
        double axial_value = floor((axial[i] - axial_min) / WR_AXIAL_BIN);
        double raw_turn_value = floor(qi + 1e-8);
        if (!isfinite(qi) || !isfinite(radius[i]) ||
            !isfinite(axial_value) ||
            axial_value < (double)INT32_MIN ||
            axial_value > (double)INT32_MAX ||
            raw_turn_value < (double)INT32_MIN ||
            raw_turn_value > (double)INT32_MAX)
            return -1;
        int phase_bin = (int)floor(phase * (double)WR_PHASE_BINS);
        if (phase_bin < 0) phase_bin = 0;
        if (phase_bin >= WR_PHASE_BINS) phase_bin = WR_PHASE_BINS - 1;
        strand[count].axial_bin = (int32_t)axial_value;
        strand[count].phase_bin = phase_bin;
        strand[count].component = c;
        strand[count].raw_turn = (int32_t)raw_turn_value;
        strand[count].count = 1;
        strand[count].radius = radius[i];
        strand[count].q = qi;
        sample_vertex[count] = (int32_t)i;
        count++;
    }
    *out_strand = strand;
    *out_vertex = sample_vertex;
    *out_count = count;
    return count ? 0 : -1;
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
               strand[last].raw_turn == strand[first].raw_turn)
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

static size_t wr_collect_order(
    const WrStrand *strand, size_t nstrand, double pitch,
    WrObservation *observation, size_t capacity,
    size_t *count, size_t *dropped)
{
    double gap_min = pitch > 1e-6 ? 0.35 * pitch : 0.0;
    double gap_max = pitch > 1e-6 ? 2.25 * pitch : DBL_MAX;
    size_t bins = 0;
    for (size_t first = 0; first < nstrand;) {
        size_t last = first + 1;
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
            if (gap < gap_min || gap > gap_max) continue;
            long target_long = lround(inner->q + 1.0 - outer->q);
            if (target_long < INT32_MIN || target_long > INT32_MAX) continue;
            double residual = fabs(
                outer->q + (double)target_long - inner->q - 1.0);
            if (residual > WR_ORDER_PHASE_TOL) continue;
            wr_append_observation(
                observation, capacity, count, dropped,
                inner->component, outer->component, (int32_t)target_long,
                residual, WR_OBS_ORDER);
        }
        first = last;
    }
    return bins;
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
        WrCandidate continuation, order;
        memset(&continuation, 0, sizeof(continuation));
        memset(&order, 0, sizeof(order));
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
            else
                order = candidate;
            kind_first = kind_last;
        }
        WrCandidate chosen;
        if (continuation.valid && continuation.eligible)
            chosen = continuation;
        else if (order.valid && order.eligible)
            chosen = order;
        else if (continuation.valid &&
                 (!order.valid ||
                  continuation.mode_observations >= order.mode_observations))
            chosen = continuation;
        else
            chosen = order;
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

static int wr_forest_corrections(
    Arena_T arena, WrRelation *relation, size_t nrelation,
    int32_t ncomponents, int32_t anchor_component,
    const int32_t *component_size, const int32_t *component,
    const double *q, size_t nvertices,
    int32_t **out_correction, int32_t **out_relation_island,
    UnionFind *out_graph,
    WindingRegisterStats *stats)
{
    UnionFind graph = UF_new(arena, ncomponents);
    size_t nforest = 0;
    for (size_t i = 0; i < nrelation; i++) {
        WrRelation *current = &relation[i];
        if (!current->eligible) continue;
        stats->eligible_relations++;
        int32_t a = uf_find(&graph, current->a);
        int32_t b = uf_find(&graph, current->b);
        if (a == b) continue;
        uf_union(&graph, a, b);
        current->selected = 1;
        nforest++;
    }
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
    if (arena == NULL || vertices == NULL || nvertices == 0 ||
        axial == NULL || radius == NULL || theta == NULL || q == NULL ||
        component == NULL || ncomponents <= 0 || component_size == NULL ||
        anchor_component < 0 || anchor_component >= ncomponents ||
        (winding_sense != -1 && winding_sense != 1) ||
        out_correction == NULL || stats == NULL)
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
        stats->continuation_components = 1;
        stats->relation_components = 1;
        stats->continuation_satisfaction = 1.0;
        stats->order_satisfaction = 1.0;
        return 0;
    }

    WrStrand *strand = NULL;
    int32_t *sample_vertex = NULL;
    size_t nsample = 0;
    if (wr_build_samples(
            arena, nvertices, axial, radius, q, component, ncomponents,
            axial_min, &strand, &sample_vertex, &nsample) != 0)
        return -1;
    if (nsample > (SIZE_MAX - 1024) / 2) return -1;
    size_t observation_capacity = 2 * nsample + 1024;
    if (observation_capacity >
        (size_t)LONG_MAX / sizeof(WrObservation))
        return -1;
    WrObservation *observation = (WrObservation *)ARENA_ALLOC(
        arena, (long)(observation_capacity * sizeof(WrObservation)));
    size_t nobservation = 0, dropped = 0;
    size_t continuation_limit = observation_capacity > nsample
                              ? observation_capacity - nsample : 0;
    wr_collect_continuations(
        vertices, radius, theta, q, component, sample_vertex, nsample,
        pitch, winding_sense, observation, observation_capacity,
        continuation_limit, &nobservation, &dropped);

    size_t nstrand = wr_collapse_strands(strand, nsample);
    if (nstrand == 0) return -1;
    stats->strands = nstrand;
    stats->bins = wr_collect_order(
        strand, nstrand, pitch, observation, observation_capacity,
        &nobservation, &dropped);
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
            &correction, out_relation_island, &graph, stats) != 0)
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
            if (uf_find(&continuation_graph, current->a) ==
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
                    "solved=%lld support=%zu/%zu residual=%.3f\n",
                    current->a, current->b,
                    current->kind == WR_OBS_CONTINUATION
                        ? "continuation" : "order",
                    current->target, (long long)solved,
                    current->mode_observations, current->observations,
                    current->residual);
        }
    }

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
            "corr=[%d,%d] sat(cont=%.3f order=%.3f) conflicts=%zu drop=%zu\n",
            winding_sense, nsample, stats->bins, stats->strands,
            stats->continuation_observations, stats->order_observations,
            stats->order_observations_suppressed,
            stats->eligible_relations, stats->relations,
            stats->order_relations_suppressed,
            stats->forest_relations, stats->relation_components,
            stats->packed_relation_components, stats->packed_mesh_components,
            stats->correction_min, stats->correction_max,
            stats->continuation_satisfaction, stats->order_satisfaction,
            stats->relation_conflicts, stats->observations_dropped);
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
            &correction, &island, &graph, &stats);
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

    Arena_dispose(&arena);
    if (fails == 0)
        fprintf(stderr, "[winding register selftest] all tests OK\n");
    return fails;
}
