#include "atlas_boxcut.h"

#include "../split/multicut_wrap.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define ABC_REPULSIVE_WEIGHT 1.0e6

typedef struct {
    uint64_t key;
    int32_t face;
    double uv_length;
} AbcEdgeRecord;

static int abc_compare_edge_record(const void *pa, const void *pb)
{
    const AbcEdgeRecord *a = (const AbcEdgeRecord *)pa;
    const AbcEdgeRecord *b = (const AbcEdgeRecord *)pb;
    if (a->key != b->key) return a->key < b->key ? -1 : 1;
    return a->face < b->face ? -1 : (a->face > b->face ? 1 : 0);
}

static int abc_compare_u64(const void *pa, const void *pb)
{
    uint64_t a = *(const uint64_t *)pa;
    uint64_t b = *(const uint64_t *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static uint64_t abc_vertex_edge_key(int32_t a, int32_t b)
{
    uint32_t lo = (uint32_t)(a < b ? a : b);
    uint32_t hi = (uint32_t)(a < b ? b : a);
    return ((uint64_t)lo << 32) | (uint64_t)hi;
}

static uint64_t abc_face_pair_key(int32_t a, int32_t b)
{
    return abc_vertex_edge_key(a, b);
}

static size_t abc_lower_edge(const AbcEdgeRecord *edge, size_t count,
                             uint64_t key)
{
    size_t lo = 0, hi = count;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (edge[mid].key < key) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

static int abc_append_key(uint64_t **key, size_t *count, size_t *capacity,
                          uint64_t value)
{
    if (*count == *capacity) {
        size_t next = *capacity ? 2 * *capacity : 256;
        if (next < *capacity || next > SIZE_MAX / sizeof(uint64_t)) return -1;
        uint64_t *grown = (uint64_t *)realloc(
            *key, next * sizeof(uint64_t));
        if (grown == NULL) return -1;
        *key = grown;
        *capacity = next;
    }
    (*key)[(*count)++] = value;
    return 0;
}

static double abc_face_radius(const AtlasBoxcutProblem *p, size_t face,
                              const double axis[3], double axis_norm2)
{
    double radius = 0.0;
    for (int k = 0; k < 3; k++) {
        int32_t vertex = p->faces[face * 3 + (size_t)k];
        const float *x = &p->vertices[(size_t)vertex * 3];
        double d[3] = {(double)x[0] - (double)p->axis_point[0],
                       (double)x[1] - (double)p->axis_point[1],
                       (double)x[2] - (double)p->axis_point[2]};
        double along = (d[0] * axis[0] + d[1] * axis[1] +
                        d[2] * axis[2]) / axis_norm2;
        double q0 = d[0] - along * axis[0];
        double q1 = d[1] - along * axis[1];
        double q2 = d[2] - along * axis[2];
        radius += sqrt(q0 * q0 + q1 * q1 + q2 * q2);
    }
    return radius / 3.0;
}

int AtlasBoxcut_partition(Arena_T arena,
                          const AtlasBoxcutProblem *p,
                          int32_t **out_face_label,
                          AtlasBoxcutStats *stats)
{
    if (arena == NULL || p == NULL || out_face_label == NULL ||
        stats == NULL || p->vertices == NULL || p->nvertices < 3 ||
        p->nvertices > (size_t)INT32_MAX || p->faces == NULL ||
        p->nfaces == 0 || p->nfaces > (size_t)INT32_MAX || p->u == NULL ||
        p->v == NULL || p->axis_point == NULL || p->axis_dir == NULL ||
        p->overlap == NULL || p->radius_gate < 0.0 ||
        !isfinite(p->radius_gate))
        return -1;
    memset(stats, 0, sizeof *stats);
    stats->faces = p->nfaces;
    stats->exact_overlap_pairs = p->overlap->npairs;
    stats->radius_gate = p->radius_gate;

    size_t nrecord = p->nfaces * 3;
    AbcEdgeRecord *record = (AbcEdgeRecord *)ARENA_ALLOC(
        arena, nrecord * sizeof(AbcEdgeRecord));
    for (size_t f = 0; f < p->nfaces; f++) {
        for (int e = 0; e < 3; e++) {
            int32_t a = p->faces[f * 3 + (size_t)e];
            int32_t b = p->faces[f * 3 + (size_t)((e + 1) % 3)];
            if (a < 0 || b < 0 || (size_t)a >= p->nvertices ||
                (size_t)b >= p->nvertices || !isfinite(p->u[a]) ||
                !isfinite(p->u[b]) || !isfinite(p->v[a]) ||
                !isfinite(p->v[b]))
                return -1;
            double du = p->u[b] - p->u[a];
            double dv = p->v[b] - p->v[a];
            AbcEdgeRecord *r = &record[f * 3 + (size_t)e];
            r->key = abc_vertex_edge_key(a, b);
            r->face = (int32_t)f;
            r->uv_length = sqrt(du * du + dv * dv);
        }
    }
    qsort(record, nrecord, sizeof(AbcEdgeRecord), abc_compare_edge_record);

    size_t nintrinsic = 0;
    for (size_t first = 0; first < nrecord;) {
        size_t last = first + 1;
        while (last < nrecord && record[last].key == record[first].key)
            last++;
        size_t n = last - first;
        if (n > 1) {
            if (n > SIZE_MAX / (n - 1) ||
                nintrinsic > SIZE_MAX - n * (n - 1) / 2)
                return -1;
            nintrinsic += n * (n - 1) / 2;
        }
        first = last;
    }

    uint64_t *seam_key = NULL;
    size_t nseam_key = 0, seam_capacity = 0;
    if (p->seam != NULL && p->seam->npairs > 1 &&
        p->seam->nbundles > 0) {
        size_t *offset = (size_t *)ARENA_CALLOC(
            arena, p->seam->nbundles + 1, sizeof(size_t));
        for (size_t i = 0; i < p->seam->npairs; i++) {
            int32_t b = p->seam->pairs[i].bundle;
            if (b < 0 || (size_t)b >= p->seam->nbundles) return -1;
            if (p->eligible_seam_bundle == NULL ||
                p->eligible_seam_bundle[b])
                offset[(size_t)b + 1]++;
        }
        for (size_t i = 0; i < p->seam->nbundles; i++)
            offset[i + 1] += offset[i];
        size_t *cursor = (size_t *)ARENA_ALLOC(
            arena, p->seam->nbundles * sizeof(size_t));
        memcpy(cursor, offset, p->seam->nbundles * sizeof(size_t));
        int32_t *index = (int32_t *)ARENA_ALLOC(
            arena, p->seam->npairs * sizeof(int32_t));
        for (size_t i = 0; i < p->seam->npairs; i++) {
            int32_t b = p->seam->pairs[i].bundle;
            if (p->eligible_seam_bundle != NULL &&
                !p->eligible_seam_bundle[b])
                continue;
            index[cursor[b]++] = (int32_t)i;
        }
        for (size_t b = 0; b < p->seam->nbundles; b++) {
            for (size_t ia = offset[b]; ia < offset[b + 1]; ia++) {
                const AtlasSeamPair *a =
                    &p->seam->pairs[index[ia]];
                for (size_t ib = ia + 1; ib < offset[b + 1]; ib++) {
                    const AtlasSeamPair *c =
                        &p->seam->pairs[index[ib]];
                    if (a->vertex0 == c->vertex0 ||
                        a->vertex1 == c->vertex1)
                        continue;
                    uint64_t key0 = abc_vertex_edge_key(
                        a->vertex0, c->vertex0);
                    uint64_t key1 = abc_vertex_edge_key(
                        a->vertex1, c->vertex1);
                    size_t f0 = abc_lower_edge(record, nrecord, key0);
                    size_t f1 = abc_lower_edge(record, nrecord, key1);
                    if (f0 == nrecord || record[f0].key != key0 ||
                        f1 == nrecord || record[f1].key != key1)
                        continue;
                    size_t l0 = f0 + 1, l1 = f1 + 1;
                    while (l0 < nrecord && record[l0].key == key0) l0++;
                    while (l1 < nrecord && record[l1].key == key1) l1++;
                    for (size_t i = f0; i < l0; i++) {
                        for (size_t j = f1; j < l1; j++) {
                            if (record[i].face == record[j].face) continue;
                            if (abc_append_key(
                                    &seam_key, &nseam_key, &seam_capacity,
                                    abc_face_pair_key(record[i].face,
                                                      record[j].face)) != 0) {
                                free(seam_key);
                                return -1;
                            }
                        }
                    }
                }
            }
        }
    }
    qsort(seam_key, nseam_key, sizeof(uint64_t), abc_compare_u64);
    size_t nseam = 0;
    for (size_t i = 0; i < nseam_key; i++)
        if (i == 0 || seam_key[i] != seam_key[i - 1])
            seam_key[nseam++] = seam_key[i];

    size_t nadjacency = nintrinsic + nseam;
    if (nadjacency > (size_t)INT32_MAX) {
        free(seam_key);
        return -1;
    }
    int32_t *adj_from = (int32_t *)ARENA_ALLOC(
        arena, (nadjacency ? nadjacency : 1) * sizeof(int32_t));
    int32_t *adj_to = (int32_t *)ARENA_ALLOC(
        arena, (nadjacency ? nadjacency : 1) * sizeof(int32_t));
    double *adj_length = (double *)ARENA_ALLOC(
        arena, (nadjacency ? nadjacency : 1) * sizeof(double));
    size_t na = 0;
    double length_sum = 0.0;
    for (size_t first = 0; first < nrecord;) {
        size_t last = first + 1;
        while (last < nrecord && record[last].key == record[first].key)
            last++;
        for (size_t i = first; i < last; i++) {
            for (size_t j = i + 1; j < last; j++) {
                adj_from[na] = record[i].face;
                adj_to[na] = record[j].face;
                adj_length[na] = 0.5 *
                    (record[i].uv_length + record[j].uv_length);
                length_sum += adj_length[na];
                na++;
            }
        }
        first = last;
    }
    double average_length = na ? length_sum / (double)na : 1.0;
    if (!isfinite(average_length) || average_length <= 1e-12)
        average_length = 1.0;
    for (size_t i = 0; i < nseam; i++) {
        adj_from[na] = (int32_t)(seam_key[i] >> 32);
        adj_to[na] = (int32_t)(seam_key[i] & UINT32_MAX);
        adj_length[na] = average_length;
        na++;
    }
    free(seam_key);
    if (na != nadjacency) return -1;

    double axis[3] = {(double)p->axis_dir[0],
                      (double)p->axis_dir[1],
                      (double)p->axis_dir[2]};
    double axis_norm2 = axis[0] * axis[0] + axis[1] * axis[1] +
                        axis[2] * axis[2];
    if (!isfinite(axis_norm2) || axis_norm2 <= 1e-15) return -1;
    double *face_radius = (double *)ARENA_ALLOC(
        arena, p->nfaces * sizeof(double));
    for (size_t f = 0; f < p->nfaces; f++)
        face_radius[f] = abc_face_radius(p, f, axis, axis_norm2);

    size_t nrepulsive = 0;
    for (size_t i = 0; i < p->overlap->npairs; i++) {
        const AtlasOverlapPair *pair = &p->overlap->pairs[i];
        if (pair->face0 < 0 || pair->face1 < 0 ||
            (size_t)pair->face0 >= p->nfaces ||
            (size_t)pair->face1 >= p->nfaces)
            return -1;
        if (p->allowed_overlap_pair != NULL &&
            p->allowed_overlap_pair[i]) {
            stats->allowed_delamination_pairs++;
            continue;
        }
        /*
         * A zero gate is the strict, fail-closed mode: every positive-area
         * UV overlap supplied by the exact audit is repulsive, including
         * pairs whose face-mean radii happen to be equal.  A positive gate
         * retains the older exploratory near-radius exemption.
         */
        if (p->radius_gate == 0.0 ||
            fabs(face_radius[pair->face1] -
                 face_radius[pair->face0]) > p->radius_gate)
            nrepulsive++;
        else
            stats->near_radius_pairs_ignored++;
    }
    if (nadjacency > (size_t)INT32_MAX - nrepulsive) return -1;
    size_t nlifted = nadjacency + nrepulsive;
    int32_t *lifted_from = (int32_t *)ARENA_ALLOC(
        arena, (nlifted ? nlifted : 1) * sizeof(int32_t));
    int32_t *lifted_to = (int32_t *)ARENA_ALLOC(
        arena, (nlifted ? nlifted : 1) * sizeof(int32_t));
    double *lifted_weight = (double *)ARENA_ALLOC(
        arena, (nlifted ? nlifted : 1) * sizeof(double));
    for (size_t i = 0; i < nadjacency; i++) {
        lifted_from[i] = adj_from[i];
        lifted_to[i] = adj_to[i];
        adj_length[i] /= average_length;
        lifted_weight[i] = adj_length[i];
    }
    size_t nr = 0;
    for (size_t i = 0; i < p->overlap->npairs; i++) {
        const AtlasOverlapPair *pair = &p->overlap->pairs[i];
        if (p->allowed_overlap_pair != NULL &&
            p->allowed_overlap_pair[i])
            continue;
        if (p->radius_gate != 0.0 &&
            fabs(face_radius[pair->face1] -
                 face_radius[pair->face0]) <= p->radius_gate)
            continue;
        lifted_from[nadjacency + nr] = pair->face0;
        lifted_to[nadjacency + nr] = pair->face1;
        lifted_weight[nadjacency + nr] = -ABC_REPULSIVE_WEIGHT;
        nr++;
    }
    if (nr != nrepulsive) return -1;

    int32_t *label = (int32_t *)ARENA_ALLOC(
        arena, p->nfaces * sizeof(int32_t));
    int32_t nlabel = 0;
    int rc = LiftedMulticut_kernighan_lin(
        (int32_t)p->nfaces, (int32_t)nadjacency, adj_from, adj_to,
        (int32_t)nlifted, lifted_from, lifted_to, lifted_weight,
        0, label, &nlabel);
    if (rc != 0 || nlabel <= 0) return -1;
    for (size_t i = 0; i < nrepulsive; i++) {
        int32_t a = lifted_from[nadjacency + i];
        int32_t b = lifted_to[nadjacency + i];
        if (label[a] == label[b]) stats->repulsive_edges_joined++;
        else stats->repulsive_edges_cut++;
    }

    stats->intrinsic_adjacency_edges = nintrinsic;
    stats->seam_adjacency_edges = nseam;
    stats->original_graph_edges = nadjacency;
    stats->radial_repulsive_edges = nrepulsive;
    stats->charts = (size_t)nlabel;
    stats->average_uv_edge_length = average_length;
    stats->adjacency_face0 = adj_from;
    stats->adjacency_face1 = adj_to;
    stats->adjacency_weight = adj_length;
    *out_face_label = label;
    return 0;
}
