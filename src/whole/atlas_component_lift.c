#include "atlas_component_lift.h"

#include "../common/union_find.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int32_t to;
    int32_t next;
    double target;
} AclAdj;

typedef struct {
    double value;
    double weight;
} AclGaugePrior;

static int acl_compare_gauge_prior(const void *pa, const void *pb)
{
    const AclGaugePrior *a = (const AclGaugePrior *)pa;
    const AclGaugePrior *b = (const AclGaugePrior *)pb;
    return a->value < b->value ? -1 : (a->value > b->value ? 1 : 0);
}

static int acl_valid_component(const AtlasComponentLiftProblem *p, int32_t c)
{
    return c >= 0 && (size_t)c < p->ncomponents;
}

int AtlasComponentLift_solve(Arena_T arena,
                             const AtlasComponentLiftProblem *p,
                             double *out_shift,
                             int32_t *out_group,
                             AtlasComponentLiftStats *stats)
{
    if (arena == NULL || p == NULL || p->ncomponents == 0 ||
        p->ncomponents > (size_t)INT32_MAX ||
        p->component_weight == NULL || out_shift == NULL ||
        stats == NULL || !isfinite(p->cycle_tolerance) ||
        p->cycle_tolerance < 0.0 ||
        (p->nprotected_edges > 0 && p->protected_edges == NULL) ||
        (p->nseparation_edges > 0 && p->separation_edges == NULL))
        return -1;
    memset(stats, 0, sizeof *stats);
    for (size_t i = 0; i < p->ncomponents; i++) {
        if (!isfinite(p->component_weight[i]) ||
            p->component_weight[i] < 0.0)
            return -1;
        if (p->component_prior_shift != NULL &&
            !isfinite(p->component_prior_shift[i]))
            return -1;
        out_shift[i] = 0.0;
    }

    UnionFind protect = UF_new(arena, (int32_t)p->ncomponents);
    for (size_t i = 0; i < p->nprotected_edges; i++) {
        const AtlasComponentLiftEdge *e = &p->protected_edges[i];
        if (!acl_valid_component(p, e->component0) ||
            !acl_valid_component(p, e->component1) ||
            !isfinite(e->target) || !isfinite(e->confidence) ||
            e->confidence < 0.0)
            return -1;
        /* Protection is an equivalence relation. */
        uf_union(&protect, e->component0, e->component1);
    }

    int32_t *root_group = (int32_t *)ARENA_ALLOC(
        arena, p->ncomponents * sizeof(int32_t));
    int32_t *component_group = (int32_t *)ARENA_ALLOC(
        arena, p->ncomponents * sizeof(int32_t));
    for (size_t i = 0; i < p->ncomponents; i++) root_group[i] = -1;
    size_t ngroup = 0;
    for (size_t i = 0; i < p->ncomponents; i++) {
        int32_t root = uf_find(&protect, (int32_t)i);
        if (root_group[root] < 0)
            root_group[root] = (int32_t)ngroup++;
        component_group[i] = root_group[root];
        if (out_group != NULL) out_group[i] = component_group[i];
    }
    stats->protected_groups = ngroup;

    double *group_weight = (double *)ARENA_CALLOC(
        arena, ngroup, sizeof(double));
    for (size_t i = 0; i < p->ncomponents; i++)
        group_weight[component_group[i]] += p->component_weight[i];

    int32_t *head = (int32_t *)ARENA_ALLOC(
        arena, ngroup * sizeof(int32_t));
    uint8_t *internal_conflict = (uint8_t *)ARENA_CALLOC(
        arena, ngroup, sizeof(uint8_t));
    for (size_t i = 0; i < ngroup; i++) head[i] = -1;
    AclAdj *adj = (AclAdj *)ARENA_ALLOC(
        arena, (2 * p->nseparation_edges + 1) * sizeof(AclAdj));
    size_t nadj = 0;
    for (size_t i = 0; i < p->nseparation_edges; i++) {
        const AtlasComponentLiftEdge *e = &p->separation_edges[i];
        if (!acl_valid_component(p, e->component0) ||
            !acl_valid_component(p, e->component1) ||
            !isfinite(e->target) || !isfinite(e->confidence) ||
            e->confidence < 0.0)
            return -1;
        int32_t a = component_group[e->component0];
        int32_t b = component_group[e->component1];
        if (a == b) {
            double residual = fabs(e->target);
            if (residual > stats->max_cycle_residual)
                stats->max_cycle_residual = residual;
            if (residual > p->cycle_tolerance) {
                internal_conflict[a] = 1;
                stats->internal_conflicts++;
            }
            continue;
        }
        adj[nadj].to = b;
        adj[nadj].target = e->target;
        adj[nadj].next = head[a];
        head[a] = (int32_t)nadj++;
        adj[nadj].to = a;
        adj[nadj].target = -e->target;
        adj[nadj].next = head[b];
        head[b] = (int32_t)nadj++;
    }

    uint8_t *visited = (uint8_t *)ARENA_CALLOC(
        arena, ngroup, sizeof(uint8_t));
    double *relative = (double *)ARENA_CALLOC(
        arena, ngroup, sizeof(double));
    double *group_shift = (double *)ARENA_CALLOC(
        arena, ngroup, sizeof(double));
    int32_t *queue = (int32_t *)ARENA_ALLOC(
        arena, ngroup * sizeof(int32_t));
    int32_t *graph_mark = (int32_t *)ARENA_ALLOC(
        arena, ngroup * sizeof(int32_t));
    for (size_t i = 0; i < ngroup; i++) graph_mark[i] = -1;
    AclGaugePrior *gauge_prior = (AclGaugePrior *)ARENA_ALLOC(
        arena, p->ncomponents * sizeof(AclGaugePrior));

    for (size_t seed = 0; seed < ngroup; seed++) {
        if (visited[seed] || (head[seed] < 0 && !internal_conflict[seed]))
            continue;
        size_t first = 0, last = 0;
        queue[last++] = (int32_t)seed;
        visited[seed] = 1;
        relative[seed] = 0.0;
        int conflict = internal_conflict[seed] != 0;
        double graph_max_residual = 0.0;

        while (first < last) {
            int32_t a = queue[first++];
            if (internal_conflict[a]) conflict = 1;
            for (int32_t ei = head[a]; ei >= 0; ei = adj[ei].next) {
                int32_t b = adj[ei].to;
                double predicted = relative[a] + adj[ei].target;
                if (!visited[b]) {
                    visited[b] = 1;
                    relative[b] = predicted;
                    queue[last++] = b;
                } else {
                    double residual = fabs(predicted - relative[b]);
                    if (residual > graph_max_residual)
                        graph_max_residual = residual;
                    if (residual > p->cycle_tolerance) conflict = 1;
                }
            }
        }
        stats->separation_graphs++;
        if (graph_max_residual > stats->max_cycle_residual)
            stats->max_cycle_residual = graph_max_residual;
        if (conflict) {
            stats->conflicting_graphs++;
            continue;
        }

        /* Keep the largest protected chart exactly at HEAD. */
        int32_t anchor = queue[0];
        double anchor_weight = group_weight[anchor];
        for (size_t i = 1; i < last; i++) {
            int32_t g = queue[i];
            if (group_weight[g] > anchor_weight ||
                (group_weight[g] == anchor_weight && g < anchor)) {
                anchor = g;
                anchor_weight = group_weight[g];
            }
        }
        double gauge = -relative[anchor];
        if (p->component_prior_shift != NULL) {
            for (size_t i = 0; i < last; i++)
                graph_mark[queue[i]] = (int32_t)stats->separation_graphs;
            size_t nprior = 0;
            double total_weight = 0.0;
            for (size_t i = 0; i < p->ncomponents; i++) {
                int32_t g = component_group[i];
                if (graph_mark[g] != (int32_t)stats->separation_graphs ||
                    p->component_weight[i] <= 0.0)
                    continue;
                gauge_prior[nprior].value =
                    p->component_prior_shift[i] - relative[g];
                gauge_prior[nprior].weight = p->component_weight[i];
                total_weight += gauge_prior[nprior].weight;
                nprior++;
            }
            if (nprior > 0 && total_weight > 0.0) {
                qsort(gauge_prior, nprior, sizeof(AclGaugePrior),
                      acl_compare_gauge_prior);
                double cumulative = 0.0;
                for (size_t i = 0; i < nprior; i++) {
                    cumulative += gauge_prior[i].weight;
                    if (2.0 * cumulative >= total_weight) {
                        gauge = gauge_prior[i].value;
                        break;
                    }
                }
            }
        }
        for (size_t i = 0; i < last; i++) {
            int32_t g = queue[i];
            group_shift[g] = relative[g] + gauge;
        }
        stats->solved_graphs++;
    }

    for (size_t i = 0; i < p->ncomponents; i++) {
        double shift = group_shift[component_group[i]];
        out_shift[i] = shift;
        if (fabs(shift) > 1e-10) {
            stats->moved_components++;
            stats->moved_weight += p->component_weight[i];
            if (fabs(shift) > stats->max_abs_shift)
                stats->max_abs_shift = fabs(shift);
        }
    }
    for (size_t g = 0; g < ngroup; g++)
        if (fabs(group_shift[g]) > 1e-10) stats->moved_groups++;
    return 0;
}

static void acl_check(int condition, const char *what, int *fails)
{
    if (condition) return;
    fprintf(stderr, "[atlas_component_lift selftest] FAIL: %s\n", what);
    (*fails)++;
}

int AtlasComponentLift_selftest(void)
{
    int fails = 0;
    {
        Arena_T arena = Arena_new();
        double weight[3] = {100.0, 10.0, 5.0};
        AtlasComponentLiftEdge protect = {1, 2, 0.0, 1.0, 0};
        AtlasComponentLiftEdge separate = {0, 1, 20.0, 1.0, 1};
        AtlasComponentLiftProblem p = {
            3, weight, NULL, &protect, 1, &separate, 1, 1e-8
        };
        double shift[3];
        int32_t group[3];
        AtlasComponentLiftStats stats;
        int rc = AtlasComponentLift_solve(
            arena, &p, shift, group, &stats);
        acl_check(rc == 0, "coherent solve", &fails);
        acl_check(fabs(shift[0]) < 1e-12 &&
                  fabs(shift[1] - 20.0) < 1e-12 &&
                  fabs(shift[2] - 20.0) < 1e-12,
                  "largest chart anchored and protected tail moved", &fails);
        acl_check(group[1] == group[2] && group[0] != group[1],
                  "protected quotient", &fails);
        acl_check(stats.solved_graphs == 1 &&
                  stats.moved_components == 2,
                  "coherent statistics", &fails);
        Arena_dispose(&arena);
    }
    {
        Arena_T arena = Arena_new();
        double weight[3] = {1.0, 1.0, 1.0};
        AtlasComponentLiftEdge edge[3] = {
            {0, 1, 10.0, 1.0, 0},
            {1, 2, 10.0, 1.0, 1},
            {0, 2, 25.0, 1.0, 2}
        };
        AtlasComponentLiftProblem p = {
            3, weight, NULL, NULL, 0, edge, 3, 1e-8
        };
        double shift[3];
        AtlasComponentLiftStats stats;
        int rc = AtlasComponentLift_solve(
            arena, &p, shift, NULL, &stats);
        acl_check(rc == 0, "conflict solve returns diagnostics", &fails);
        acl_check(fabs(shift[0]) < 1e-12 &&
                  fabs(shift[1]) < 1e-12 &&
                  fabs(shift[2]) < 1e-12,
                  "conflicting cycle stays at HEAD", &fails);
        acl_check(stats.conflicting_graphs == 1 &&
                  stats.max_cycle_residual >= 5.0 - 1e-12,
                  "conflicting cycle reported", &fails);
        Arena_dispose(&arena);
    }
    {
        Arena_T arena = Arena_new();
        double weight[2] = {100.0, 10.0};
        double prior[2] = {-20.0, 0.0};
        AtlasComponentLiftEdge edge = {0, 1, 20.0, 1.0, 0};
        AtlasComponentLiftProblem p = {
            2, weight, prior, NULL, 0, &edge, 1, 1e-8
        };
        double shift[2];
        AtlasComponentLiftStats stats;
        int rc = AtlasComponentLift_solve(
            arena, &p, shift, NULL, &stats);
        acl_check(rc == 0 && fabs(shift[0] + 20.0) < 1e-12 &&
                  fabs(shift[1]) < 1e-12,
                  "prior chooses the already-correct absolute side", &fails);
        Arena_dispose(&arena);
    }
    fprintf(stderr, "[atlas_component_lift selftest] %s\n",
            fails == 0 ? "PASSED" : "FAILED");
    return fails;
}
