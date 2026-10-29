/*
 * topology_invariants.c -- scalable component topology and H_1 generators.
 *
 * The implementation deliberately stays combinatorial.  Geometry is used only
 * to weight generator representatives; all certification decisions depend on
 * the triangle complex.
 */
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "topology_invariants.h"
#include "../common/arena.h"
#include "../common/mesh_manifold.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int32_t a, b;
    int32_t f0, f1;
    uint32_t count;
    uint8_t dir0, dir1;
    uint8_t tree, cotree;
} AuditEdge;

typedef struct {
    uint64_t *key;
    int32_t *edge;
    size_t cap, used;
} EdgeHash;

typedef struct {
    int32_t root;
    size_t faces;
} ComponentSeed;

typedef struct {
    double key;
    int32_t edge;
} WeightedEdge;

typedef struct {
    int32_t vertex;
    int32_t edge;
} VertexAdj;

typedef struct {
    double distance;
    int32_t vertex;
} HeapNode;

typedef struct {
    HeapNode *data;
    size_t size, cap;
} MinHeap;

typedef struct {
    int32_t edge;
    int32_t component;
    size_t nvertices;
    double length;
    double rooted_length;
} GeneratorCandidate;

typedef struct {
    int32_t *entry;
    size_t count;
} SparseRow;

static int cmp_boundary_invariant(const void *pa, const void *pb)
{
    const TopologyBoundaryInvariant *a =
        (const TopologyBoundaryInvariant *)pa;
    const TopologyBoundaryInvariant *b =
        (const TopologyBoundaryInvariant *)pb;
    if (a->component != b->component)
        return a->component < b->component ? -1 : 1;
    if (a->component_perimeter != b->component_perimeter)
        return a->component_perimeter ? -1 : 1;
    if (a->length != b->length) return a->length > b->length ? -1 : 1;
    if (a->edges != b->edges) return a->edges > b->edges ? -1 : 1;
    return a->root_vertex < b->root_vertex ? -1 :
           (a->root_vertex > b->root_vertex);
}

static int size_mul(size_t a, size_t b, size_t *out)
{
    if (a != 0 && b > SIZE_MAX / a) return -1;
    *out = a * b;
    return 0;
}

static void *checked_malloc(size_t count, size_t size)
{
    size_t bytes = 0;
    if (size_mul(count, size, &bytes) != 0) return NULL;
    if (bytes == 0) bytes = 1;
    return malloc(bytes);
}

static void *checked_calloc(size_t count, size_t size)
{
    if (count != 0 && size > SIZE_MAX / count) return NULL;
    if (count == 0) count = 1;
    return calloc(count, size);
}

static uint64_t edge_key(int32_t a, int32_t b)
{
    uint32_t lo = (uint32_t)(a < b ? a : b);
    uint32_t hi = (uint32_t)(a < b ? b : a);
    return ((uint64_t)lo << 32) | (uint64_t)hi;
}

static size_t hash_u64(uint64_t key, size_t mask)
{
    uint64_t x = key + UINT64_C(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    x ^= x >> 31;
    return (size_t)x & mask;
}

static int edge_hash_init(EdgeHash *hash, size_t hint)
{
    size_t cap = 16;
    while (cap < hint) {
        if (cap > SIZE_MAX / 2) return -1;
        cap *= 2;
    }
    hash->key = (uint64_t *)checked_malloc(cap, sizeof(uint64_t));
    hash->edge = (int32_t *)checked_malloc(cap, sizeof(int32_t));
    if (!hash->key || !hash->edge) {
        free(hash->key); free(hash->edge);
        memset(hash, 0, sizeof(*hash));
        return -1;
    }
    memset(hash->key, 0xff, cap * sizeof(uint64_t));
    hash->cap = cap;
    hash->used = 0;
    return 0;
}

static void edge_hash_dispose(EdgeHash *hash)
{
    free(hash->key);
    free(hash->edge);
    memset(hash, 0, sizeof(*hash));
}

static int edge_hash_rehash(EdgeHash *hash)
{
    EdgeHash next;
    size_t newcap;
    memset(&next, 0, sizeof(next));
    if (hash->cap > SIZE_MAX / 2) return -1;
    newcap = hash->cap * 2;
    if (edge_hash_init(&next, newcap) != 0) return -1;
    for (size_t i = 0; i < hash->cap; i++) {
        uint64_t key = hash->key[i];
        size_t slot;
        if (key == UINT64_MAX) continue;
        slot = hash_u64(key, next.cap - 1);
        while (next.key[slot] != UINT64_MAX)
            slot = (slot + 1) & (next.cap - 1);
        next.key[slot] = key;
        next.edge[slot] = hash->edge[i];
        next.used++;
    }
    edge_hash_dispose(hash);
    *hash = next;
    return 0;
}

static int32_t edge_hash_find(const EdgeHash *hash, uint64_t key)
{
    size_t slot = hash_u64(key, hash->cap - 1);
    for (;;) {
        if (hash->key[slot] == UINT64_MAX) return -1;
        if (hash->key[slot] == key) return hash->edge[slot];
        slot = (slot + 1) & (hash->cap - 1);
    }
}

/* Returns 1 for a new slot, 0 for an existing slot, -1 on allocation failure. */
static int edge_hash_insert(EdgeHash *hash, uint64_t key, int32_t value,
                            int32_t *out_value)
{
    size_t slot;
    if ((hash->used + 1) * 10 > hash->cap * 7) {
        if (edge_hash_rehash(hash) != 0) return -1;
    }
    slot = hash_u64(key, hash->cap - 1);
    for (;;) {
        if (hash->key[slot] == UINT64_MAX) {
            hash->key[slot] = key;
            hash->edge[slot] = value;
            hash->used++;
            *out_value = value;
            return 1;
        }
        if (hash->key[slot] == key) {
            *out_value = hash->edge[slot];
            return 0;
        }
        slot = (slot + 1) & (hash->cap - 1);
    }
}

static int32_t local_find(int32_t *parent, int32_t x)
{
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

static int local_union(int32_t *parent, int32_t *size, int32_t a, int32_t b)
{
    a = local_find(parent, a);
    b = local_find(parent, b);
    if (a == b) return 0;
    if (size[a] < size[b]) {
        int32_t t = a; a = b; b = t;
    }
    parent[b] = a;
    size[a] += size[b];
    return 1;
}

static int dual_union(int32_t *parent, uint8_t *rank, int32_t a, int32_t b)
{
    int32_t ra = a, rb = b;
    while (parent[ra] != ra) {
        parent[ra] = parent[parent[ra]];
        ra = parent[ra];
    }
    while (parent[rb] != rb) {
        parent[rb] = parent[parent[rb]];
        rb = parent[rb];
    }
    if (ra == rb) return 0;
    if (rank[ra] < rank[rb]) {
        int32_t t = ra; ra = rb; rb = t;
    }
    parent[rb] = ra;
    if (rank[ra] == rank[rb]) rank[ra]++;
    return 1;
}

static int cmp_component_seed(const void *pa, const void *pb)
{
    const ComponentSeed *a = (const ComponentSeed *)pa;
    const ComponentSeed *b = (const ComponentSeed *)pb;
    if (a->faces != b->faces) return a->faces < b->faces ? 1 : -1;
    return a->root < b->root ? -1 : (a->root > b->root);
}

static int cmp_weight_desc(const void *pa, const void *pb)
{
    const WeightedEdge *a = (const WeightedEdge *)pa;
    const WeightedEdge *b = (const WeightedEdge *)pb;
    if (a->key < b->key) return 1;
    if (a->key > b->key) return -1;
    return a->edge < b->edge ? -1 : (a->edge > b->edge);
}

static int cmp_generator_candidate(const void *pa, const void *pb)
{
    const GeneratorCandidate *a = (const GeneratorCandidate *)pa;
    const GeneratorCandidate *b = (const GeneratorCandidate *)pb;
    if (a->component != b->component)
        return a->component < b->component ? -1 : 1;
    if (a->length < b->length) return -1;
    if (a->length > b->length) return 1;
    return a->edge < b->edge ? -1 : (a->edge > b->edge);
}

static int cmp_i32(const void *pa, const void *pb)
{
    int32_t a = *(const int32_t *)pa;
    int32_t b = *(const int32_t *)pb;
    return a < b ? -1 : (a > b);
}

static int sparse_xor(int32_t **left, size_t *nleft,
                      const int32_t *right, size_t nright)
{
    int32_t *result;
    size_t i = 0, j = 0, n = 0;
    if (*nleft > SIZE_MAX - nright) return -1;
    result = (int32_t *)checked_malloc(*nleft + nright, sizeof(int32_t));
    if (!result) return -1;
    while (i < *nleft || j < nright) {
        if (j == nright || (i < *nleft && (*left)[i] < right[j])) {
            result[n++] = (*left)[i++];
        } else if (i == *nleft || right[j] < (*left)[i]) {
            result[n++] = right[j++];
        } else {
            i++; j++; /* equal entries cancel over F_2 */
        }
    }
    free(*left);
    *left = result;
    *nleft = n;
    return 0;
}

/* Exact beta_2 is the nullity of the face-to-edge boundary map.  Rows with
 * one incident face force that face coefficient to zero; rows with two faces
 * identify the two coefficients.  Only edges with 3+ incident faces remain as
 * sparse parity equations, so Gaussian elimination is confined to the actual
 * non-manifold core instead of the multi-million-face mesh. */
static int compute_exact_betti(size_t nf, const int32_t *faces,
                               const int32_t *face_edges,
                               const AuditEdge *edges, size_t ne,
                               const int32_t *vertex_component,
                               TopologyAuditReport *report)
{
    int32_t *face_parent = NULL, *face_size = NULL;
    uint8_t *forced = NULL;
    int32_t *root_to_variable = NULL;
    size_t *variable_count = NULL, *equation_rank = NULL;
    size_t *nm_offset = NULL, *nm_cursor = NULL;
    int32_t *nm_face = NULL;
    SparseRow *pivot = NULL;
    size_t total_variables = 0, total_nm_incidence = 0;
    int rc = -1;

    face_parent = (int32_t *)checked_malloc(nf, sizeof(int32_t));
    face_size = (int32_t *)checked_malloc(nf, sizeof(int32_t));
    forced = (uint8_t *)checked_calloc(nf, sizeof(uint8_t));
    root_to_variable = (int32_t *)checked_malloc(nf, sizeof(int32_t));
    variable_count = (size_t *)checked_calloc(report->face_components,
                                              sizeof(size_t));
    equation_rank = (size_t *)checked_calloc(report->face_components,
                                             sizeof(size_t));
    if (!face_parent || !face_size || !forced || !root_to_variable ||
        !variable_count || !equation_rank) goto cleanup;
    for (size_t f = 0; f < nf; f++) {
        face_parent[f] = (int32_t)f;
        face_size[f] = 1;
        root_to_variable[f] = -1;
    }
    for (size_t e = 0; e < ne; e++) {
        if (edges[e].count == 2)
            local_union(face_parent, face_size, edges[e].f0, edges[e].f1);
    }
    for (size_t e = 0; e < ne; e++) {
        if (edges[e].count == 1)
            forced[local_find(face_parent, edges[e].f0)] = 1;
    }
    for (size_t f = 0; f < nf; f++) {
        int32_t root = local_find(face_parent, (int32_t)f);
        int32_t component;
        if (root != (int32_t)f || forced[root]) continue;
        if (total_variables > (size_t)INT32_MAX) goto cleanup;
        component = vertex_component[faces[f*3]];
        if (component < 0) goto cleanup;
        root_to_variable[root] = (int32_t)total_variables++;
        variable_count[component]++;
    }

    nm_offset = (size_t *)checked_malloc(ne + 1, sizeof(size_t));
    nm_cursor = (size_t *)checked_malloc(ne, sizeof(size_t));
    if (!nm_offset || !nm_cursor) goto cleanup;
    nm_offset[0] = 0;
    for (size_t e = 0; e < ne; e++) {
        size_t add = edges[e].count > 2 ? (size_t)edges[e].count : 0;
        if (nm_offset[e] > SIZE_MAX - add) goto cleanup;
        nm_offset[e+1] = nm_offset[e] + add;
        nm_cursor[e] = nm_offset[e];
    }
    total_nm_incidence = nm_offset[ne];
    nm_face = (int32_t *)checked_malloc(total_nm_incidence, sizeof(int32_t));
    if (!nm_face) goto cleanup;
    for (size_t f = 0; f < nf; f++) {
        for (int k = 0; k < 3; k++) {
            int32_t e = face_edges[f*3+(size_t)k];
            if (edges[e].count > 2) {
                if (nm_cursor[e] >= nm_offset[(size_t)e+1]) goto cleanup;
                nm_face[nm_cursor[e]++] = (int32_t)f;
            }
        }
    }
    for (size_t e = 0; e < ne; e++)
        if (nm_cursor[e] != nm_offset[e+1]) goto cleanup;

    pivot = (SparseRow *)checked_calloc(total_variables, sizeof(SparseRow));
    if (!pivot) goto cleanup;
    for (size_t e = 0; e < ne; e++) {
        int32_t *row;
        size_t count = 0, reduced = 0;
        int32_t component;
        if (edges[e].count <= 2) continue;
        row = (int32_t *)checked_malloc((size_t)edges[e].count,
                                        sizeof(int32_t));
        if (!row) goto cleanup;
        for (size_t p = nm_offset[e]; p < nm_offset[e+1]; p++) {
            int32_t root = local_find(face_parent, nm_face[p]);
            int32_t variable;
            if (forced[root]) continue;
            variable = root_to_variable[root];
            if (variable < 0) { free(row); goto cleanup; }
            row[count++] = variable;
        }
        qsort(row, count, sizeof(int32_t), cmp_i32);
        for (size_t i = 0; i < count; ) {
            size_t j = i + 1;
            while (j < count && row[j] == row[i]) j++;
            if (((j - i) & 1u) != 0) row[reduced++] = row[i];
            i = j;
        }
        count = reduced;
        component = vertex_component[edges[e].a];
        while (count > 0) {
            int32_t lead = row[count - 1];
            if (pivot[lead].count == 0) {
                pivot[lead].entry = row;
                pivot[lead].count = count;
                row = NULL;
                equation_rank[component]++;
                break;
            }
            if (sparse_xor(&row, &count, pivot[lead].entry,
                           pivot[lead].count) != 0) {
                free(row);
                goto cleanup;
            }
        }
        free(row);
    }

    report->beta_1 = 0;
    report->beta_2 = 0;
    report->betti_complete = 1;
    report->minimal_generator_rank = 0;
    for (size_t c = 0; c < report->face_components; c++) {
        TopologyComponentInvariant *component = &report->component[c];
        if (equation_rank[c] > variable_count[c]) goto cleanup;
        component->beta_2 =
            (int64_t)(variable_count[c] - equation_rank[c]);
        component->beta_1 =
            1 - component->euler_characteristic + component->beta_2;
        if (component->beta_1 < 0) goto cleanup;
        component->minimal_generator_rank = (size_t)component->beta_1;
        report->beta_1 += component->beta_1;
        report->beta_2 += component->beta_2;
        report->minimal_generator_rank += component->minimal_generator_rank;
    }
    rc = 0;

cleanup:
    if (pivot) {
        for (size_t i = 0; i < total_variables; i++)
            free(pivot[i].entry);
    }
    free(pivot);
    free(face_parent); free(face_size); free(forced);
    free(root_to_variable); free(variable_count); free(equation_rank);
    free(nm_offset); free(nm_cursor); free(nm_face);
    return rc;
}

static double audit_edge_length(const float *verts, const AuditEdge *edge)
{
    if (!verts) return 1.0;
    {
        double dx = (double)verts[(size_t)edge->a*3+0] -
                    (double)verts[(size_t)edge->b*3+0];
        double dy = (double)verts[(size_t)edge->a*3+1] -
                    (double)verts[(size_t)edge->b*3+1];
        double dz = (double)verts[(size_t)edge->a*3+2] -
                    (double)verts[(size_t)edge->b*3+2];
        return sqrt(dx*dx + dy*dy + dz*dz);
    }
}

static int heap_push(MinHeap *heap, double distance, int32_t vertex)
{
    size_t i;
    if (heap->size == heap->cap) {
        size_t cap = heap->cap ? heap->cap * 2 : 1024;
        HeapNode *next;
        if (cap < heap->cap) return -1;
        next = (HeapNode *)realloc(heap->data, cap * sizeof(HeapNode));
        if (!next) return -1;
        heap->data = next;
        heap->cap = cap;
    }
    i = heap->size++;
    while (i > 0) {
        size_t p = (i - 1) / 2;
        if (heap->data[p].distance < distance ||
            (heap->data[p].distance == distance &&
             heap->data[p].vertex <= vertex)) break;
        heap->data[i] = heap->data[p];
        i = p;
    }
    heap->data[i].distance = distance;
    heap->data[i].vertex = vertex;
    return 0;
}

static HeapNode heap_pop(MinHeap *heap)
{
    HeapNode result = heap->data[0];
    HeapNode tail = heap->data[--heap->size];
    size_t i = 0;
    if (heap->size == 0) return result;
    while (i * 2 + 1 < heap->size) {
        size_t child = i * 2 + 1;
        if (child + 1 < heap->size &&
            (heap->data[child+1].distance < heap->data[child].distance ||
             (heap->data[child+1].distance == heap->data[child].distance &&
              heap->data[child+1].vertex < heap->data[child].vertex)))
            child++;
        if (heap->data[child].distance > tail.distance ||
            (heap->data[child].distance == tail.distance &&
             heap->data[child].vertex >= tail.vertex)) break;
        heap->data[i] = heap->data[child];
        i = child;
    }
    heap->data[i] = tail;
    return result;
}

static void set_error(TopologyAuditReport *report, const char *message)
{
    if (!report) return;
    strncpy(report->error, message, sizeof(report->error) - 1);
    report->error[sizeof(report->error) - 1] = '\0';
}

void TopologyAudit_options_default(TopologyAuditOptions *options)
{
    if (!options) return;
    options->emit_generators = 0;
    options->max_generators_per_component = 0;
}

void TopologyAudit_dispose(TopologyAuditReport *report)
{
    if (!report) return;
    if (report->generator) {
        for (size_t i = 0; i < report->emitted_generators; i++)
            free(report->generator[i].vertices);
    }
    free(report->generator);
    free(report->component);
    free(report->boundary);
    free(report->vertex_component);
    memset(report, 0, sizeof(*report));
}

static int fail_report(TopologyAuditReport *report, const char *message)
{
    char saved[256];
    strncpy(saved, message, sizeof(saved) - 1);
    saved[sizeof(saved) - 1] = '\0';
    TopologyAudit_dispose(report);
    set_error(report, saved);
    return -1;
}

/* Count and, optionally, materialize the simple tree path closed by edge. */
static int generator_path(const float *verts, const AuditEdge *edge,
                          const int32_t *parent,
                          const int32_t *depth, const double *distance,
                          size_t *out_count, double *out_length,
                          int32_t *out_vertices)
{
    int32_t u = edge->a, v = edge->b;
    int32_t lca_u = u, lca_v = v;
    size_t left = 1, right = 0;
    while (depth[lca_u] > depth[lca_v]) {
        lca_u = parent[lca_u]; left++;
    }
    while (depth[lca_v] > depth[lca_u]) {
        lca_v = parent[lca_v]; right++;
    }
    while (lca_u != lca_v) {
        lca_u = parent[lca_u];
        lca_v = parent[lca_v];
        left++; right++;
    }
    *out_count = left + right;
    *out_length = distance[u] + distance[v] -
                  2.0 * distance[lca_u] + audit_edge_length(verts, edge);

    if (out_vertices) {
        int32_t x = u, y = v;
        size_t li = 0, ri = left + right;
        while (x != lca_u) {
            out_vertices[li++] = x;
            x = parent[x];
        }
        out_vertices[li++] = lca_u;
        while (y != lca_u) {
            out_vertices[--ri] = y;
            y = parent[y];
        }
    }
    return 0;
}

static int build_generators(const float *verts, size_t nv, size_t nf,
                            AuditEdge *edges, size_t ne,
                            const int32_t *vertex_component,
                            const TopologyAuditOptions *options,
                            TopologyAuditReport *report)
{
    uint8_t *eligible = NULL;
    size_t *degree = NULL, *offset = NULL, *cursor = NULL;
    VertexAdj *adj = NULL;
    double *distance = NULL;
    int32_t *parent = NULL, *parent_edge = NULL, *depth = NULL;
    MinHeap heap;
    WeightedEdge *weighted = NULL;
    size_t nweighted = 0;
    int32_t *dual_parent = NULL;
    uint8_t *dual_rank = NULL;
    size_t *cotree_count = NULL, *tree_count = NULL, *reached_count = NULL;
    GeneratorCandidate *candidate = NULL;
    size_t ncandidate = 0;
    size_t ncomponents = report->face_components;
    size_t nadj = 0, ndual;
    int rc = -1;

    memset(&heap, 0, sizeof(heap));
    eligible = (uint8_t *)checked_calloc(ncomponents, sizeof(uint8_t));
    tree_count = (size_t *)checked_calloc(ncomponents, sizeof(size_t));
    reached_count = (size_t *)checked_calloc(ncomponents, sizeof(size_t));
    cotree_count = (size_t *)checked_calloc(ncomponents, sizeof(size_t));
    if (!eligible || !tree_count || !reached_count || !cotree_count) goto cleanup;

    for (size_t c = 0; c < ncomponents; c++) {
        TopologyComponentInvariant *component = &report->component[c];
        component->generator_basis_complete = component->surface_valid ? 1 : 0;
        if (component->surface_valid && component->beta_1 > 0) eligible[c] = 1;
    }

    /* A multi-source Dijkstra is one independent shortest-path tree per
     * connected component because no adjacency crosses component boundaries. */
    degree = (size_t *)checked_calloc(nv, sizeof(size_t));
    if (!degree) goto cleanup;
    for (size_t e = 0; e < ne; e++) {
        int32_t c = vertex_component[edges[e].a];
        if (c >= 0 && eligible[c]) {
            degree[edges[e].a]++;
            degree[edges[e].b]++;
            nadj += 2;
        }
    }
    offset = (size_t *)checked_malloc(nv + 1, sizeof(size_t));
    cursor = (size_t *)checked_malloc(nv, sizeof(size_t));
    adj = (VertexAdj *)checked_malloc(nadj, sizeof(VertexAdj));
    distance = (double *)checked_malloc(nv, sizeof(double));
    parent = (int32_t *)checked_malloc(nv, sizeof(int32_t));
    parent_edge = (int32_t *)checked_malloc(nv, sizeof(int32_t));
    depth = (int32_t *)checked_calloc(nv, sizeof(int32_t));
    if (!offset || !cursor || !adj || !distance || !parent ||
        !parent_edge || !depth) goto cleanup;
    offset[0] = 0;
    for (size_t v = 0; v < nv; v++) {
        offset[v+1] = offset[v] + degree[v];
        cursor[v] = offset[v];
        distance[v] = DBL_MAX;
        parent[v] = -1;
        parent_edge[v] = -1;
    }
    for (size_t e = 0; e < ne; e++) {
        int32_t a = edges[e].a, b = edges[e].b;
        int32_t c = vertex_component[a];
        if (c < 0 || !eligible[c]) continue;
        adj[cursor[a]].vertex = b;
        adj[cursor[a]].edge = (int32_t)e;
        cursor[a]++;
        adj[cursor[b]].vertex = a;
        adj[cursor[b]].edge = (int32_t)e;
        cursor[b]++;
    }
    for (size_t c = 0; c < ncomponents; c++) {
        int32_t root;
        if (!eligible[c]) continue;
        root = report->component[c].root_vertex;
        if (root < 0) goto cleanup;
        distance[root] = 0.0;
        if (heap_push(&heap, 0.0, root) != 0) goto cleanup;
    }
    while (heap.size) {
        HeapNode node = heap_pop(&heap);
        int32_t v = node.vertex;
        if (node.distance != distance[v]) continue;
        for (size_t p = offset[v]; p < offset[v+1]; p++) {
            int32_t w = adj[p].vertex;
            int32_t e = adj[p].edge;
            double next = node.distance + audit_edge_length(verts, &edges[e]);
            if (next < distance[w]) {
                distance[w] = next;
                parent[w] = v;
                parent_edge[w] = e;
                depth[w] = depth[v] + 1;
                if (heap_push(&heap, next, w) != 0) goto cleanup;
            }
        }
    }
    for (size_t v = 0; v < nv; v++) {
        int32_t c = vertex_component[v];
        if (c < 0 || !eligible[c]) continue;
        if (distance[v] < DBL_MAX) reached_count[c]++;
        if (parent_edge[v] >= 0) {
            edges[parent_edge[v]].tree = 1;
            tree_count[c]++;
        }
    }

    for (size_t c = 0; c < ncomponents; c++) {
        TopologyComponentInvariant *component = &report->component[c];
        if (!eligible[c]) continue;
        if (reached_count[c] != component->vertices ||
            tree_count[c] + 1 != component->vertices) {
            component->generator_basis_complete = 0;
            component->defect_mask |= TOPOLOGY_DEFECT_GENERATOR_FAILURE;
        }
    }

    /* The extended dual has one outside node per component.  Boundary primal
     * edges connect their incident face to that outside node. */
    if (nf > (size_t)INT32_MAX ||
        ncomponents > (size_t)INT32_MAX - nf) goto cleanup;
    ndual = nf + ncomponents;
    dual_parent = (int32_t *)checked_malloc(ndual, sizeof(int32_t));
    dual_rank = (uint8_t *)checked_calloc(ndual, sizeof(uint8_t));
    weighted = (WeightedEdge *)checked_malloc(ne, sizeof(WeightedEdge));
    if (!dual_parent || !dual_rank || !weighted) goto cleanup;
    for (size_t i = 0; i < ndual; i++) dual_parent[i] = (int32_t)i;

    for (size_t e = 0; e < ne; e++) {
        AuditEdge *edge = &edges[e];
        int32_t c = vertex_component[edge->a];
        if (c < 0 || !eligible[c] || edge->tree) continue;
        if (edge->count != 1 && edge->count != 2) continue;
        weighted[nweighted].edge = (int32_t)e;
        weighted[nweighted].key = distance[edge->a] +
                                  audit_edge_length(verts, edge) +
                                  distance[edge->b];
        nweighted++;
    }
    qsort(weighted, nweighted, sizeof(WeightedEdge), cmp_weight_desc);
    for (size_t i = 0; i < nweighted; i++) {
        AuditEdge *edge = &edges[weighted[i].edge];
        int32_t c = vertex_component[edge->a];
        int32_t d0 = edge->f0;
        int32_t d1 = edge->count == 2 ? edge->f1 :
                     (int32_t)(nf + (size_t)c);
        if (dual_union(dual_parent, dual_rank, d0, d1)) {
            edge->cotree = 1;
            cotree_count[c]++;
        }
    }
    for (size_t c = 0; c < ncomponents; c++) {
        TopologyComponentInvariant *component = &report->component[c];
        size_t expected;
        if (!eligible[c]) continue;
        expected = component->faces - (component->boundary_loops == 0 ? 1u : 0u);
        if (cotree_count[c] != expected) {
            component->generator_basis_complete = 0;
            component->defect_mask |= TOPOLOGY_DEFECT_GENERATOR_FAILURE;
        }
    }

    candidate = (GeneratorCandidate *)checked_malloc(ne, sizeof(GeneratorCandidate));
    if (!candidate) goto cleanup;
    for (size_t e = 0; e < ne; e++) {
        AuditEdge *edge = &edges[e];
        int32_t c = vertex_component[edge->a];
        size_t nvertices;
        double length;
        if (c < 0 || !eligible[c] || edge->tree || edge->cotree) continue;
        generator_path(verts, edge, parent, depth, distance,
                       &nvertices, &length, NULL);
        candidate[ncandidate].edge = (int32_t)e;
        candidate[ncandidate].component = c;
        candidate[ncandidate].nvertices = nvertices;
        candidate[ncandidate].length = length;
        candidate[ncandidate].rooted_length =
            distance[edge->a] + audit_edge_length(verts, edge) +
            distance[edge->b];
        ncandidate++;
    }

    {
        size_t *basis_count = (size_t *)checked_calloc(ncomponents, sizeof(size_t));
        if (!basis_count) goto cleanup;
        for (size_t i = 0; i < ncandidate; i++)
            basis_count[candidate[i].component]++;
        for (size_t c = 0; c < ncomponents; c++) {
            TopologyComponentInvariant *component = &report->component[c];
            if (!eligible[c]) continue;
            if ((int64_t)basis_count[c] != component->beta_1) {
                component->generator_basis_complete = 0;
                component->defect_mask |= TOPOLOGY_DEFECT_GENERATOR_FAILURE;
            }
        }
        free(basis_count);
    }

    qsort(candidate, ncandidate, sizeof(GeneratorCandidate),
          cmp_generator_candidate);

    if (options->emit_generators && ncandidate > 0) {
        size_t emit_count = 0;
        size_t last_component = SIZE_MAX, ordinal = 0;
        for (size_t i = 0; i < ncandidate; i++) {
            size_t c = (size_t)candidate[i].component;
            if (c != last_component) {
                last_component = c;
                ordinal = 0;
            }
            if (options->max_generators_per_component == 0 ||
                ordinal < options->max_generators_per_component)
                emit_count++;
            ordinal++;
        }
        report->generator = (TopologyGenerator *)
            checked_calloc(emit_count, sizeof(TopologyGenerator));
        if (!report->generator) goto cleanup;
        /* Set the owned slot count before filling it so a mid-loop allocation
         * failure can dispose every already-materialized vertex sequence. */
        report->emitted_generators = emit_count;

        last_component = SIZE_MAX;
        ordinal = 0;
        emit_count = 0;
        for (size_t i = 0; i < ncandidate; i++) {
            GeneratorCandidate *source = &candidate[i];
            size_t c = (size_t)source->component;
            TopologyGenerator *target;
            AuditEdge *edge;
            if (c != last_component) {
                last_component = c;
                ordinal = 0;
            }
            if (options->max_generators_per_component != 0 &&
                ordinal >= options->max_generators_per_component) {
                ordinal++;
                continue;
            }
            target = &report->generator[emit_count++];
            edge = &edges[source->edge];
            target->component = c;
            target->ordinal = ordinal;
            target->closing_edge_a = edge->a;
            target->closing_edge_b = edge->b;
            target->length = source->length;
            target->rooted_length = source->rooted_length;
            target->nvertices = source->nvertices;
            target->vertices = (int32_t *)
                checked_malloc(source->nvertices, sizeof(int32_t));
            if (!target->vertices) goto cleanup;
            generator_path(verts, edge, parent, depth, distance,
                           &target->nvertices, &target->length,
                           target->vertices);
            report->component[c].emitted_generators++;
            ordinal++;
        }
        report->emitted_generators = emit_count;
    }

    report->generator_basis_complete = report->betti_complete;
    for (size_t c = 0; c < ncomponents; c++) {
        TopologyComponentInvariant *component = &report->component[c];
        if (!component->generator_basis_complete)
            report->generator_basis_complete = 0;
        if (component->emitted_generators < component->minimal_generator_rank &&
            options->emit_generators)
            report->generators_truncated = 1;
    }
    rc = 0;

cleanup:
    free(eligible);
    free(degree); free(offset); free(cursor); free(adj);
    free(distance); free(parent); free(parent_edge); free(depth);
    free(heap.data);
    free(weighted);
    free(dual_parent); free(dual_rank);
    free(cotree_count); free(tree_count); free(reached_count);
    free(candidate);
    return rc;
}

int TopologyAudit_analyze(const float *verts, size_t nv,
                          const int32_t *faces, size_t nf,
                          const TopologyAuditOptions *options_in,
                          TopologyAuditReport *report)
{
    TopologyAuditOptions defaults;
    const TopologyAuditOptions *options = options_in;
    int32_t *parent = NULL, *usize = NULL;
    uint8_t *used = NULL;
    size_t *root_faces = NULL;
    ComponentSeed *seed = NULL;
    int32_t *root_to_component = NULL, *vertex_component = NULL;
    EdgeHash hash;
    AuditEdge *edges = NULL;
    size_t ne = 0, edge_cap = 0;
    int32_t *face_edges = NULL;
    int32_t *boundary_parent = NULL, *boundary_size = NULL;
    uint32_t *boundary_degree = NULL;
    int32_t *boundary_root_to_index = NULL;
    int32_t *component_perimeter_index = NULL;
    size_t *boundary_vertex_offset = NULL, *boundary_vertex_cursor = NULL;
    int32_t *boundary_vertices = NULL;
    uint8_t *nonmanifold_vertex = NULL;
    int8_t *face_flip = NULL;
    int32_t *face_queue = NULL;
    Arena_T manifold_arena = NULL;
    size_t ncomponents = 0;
    int rc = -1;
    const char *failure = "out of memory";

    memset(&hash, 0, sizeof(hash));
    if (!report) return -1;
    memset(report, 0, sizeof(*report));
    if (!options) {
        TopologyAudit_options_default(&defaults);
        options = &defaults;
    }
    if (!faces || nv == 0 || nf == 0) {
        set_error(report, "empty triangle mesh");
        return -1;
    }
    if (nv > (size_t)INT32_MAX || nf > (size_t)INT32_MAX) {
        set_error(report, "mesh exceeds int32 topology limits");
        return -1;
    }

    parent = (int32_t *)checked_malloc(nv, sizeof(int32_t));
    usize = (int32_t *)checked_malloc(nv, sizeof(int32_t));
    used = (uint8_t *)checked_calloc(nv, sizeof(uint8_t));
    root_faces = (size_t *)checked_calloc(nv, sizeof(size_t));
    if (!parent || !usize || !used || !root_faces) goto cleanup;
    for (size_t v = 0; v < nv; v++) {
        parent[v] = (int32_t)v;
        usize[v] = 1;
        if (verts && (!isfinite((double)verts[v*3+0]) ||
                      !isfinite((double)verts[v*3+1]) ||
                      !isfinite((double)verts[v*3+2]))) {
            failure = "non-finite vertex coordinate";
            goto cleanup;
        }
    }
    for (size_t f = 0; f < nf; f++) {
        int32_t a = faces[f*3+0], b = faces[f*3+1], c = faces[f*3+2];
        if (a < 0 || b < 0 || c < 0 ||
            (size_t)a >= nv || (size_t)b >= nv || (size_t)c >= nv) {
            failure = "face index outside the vertex array";
            goto cleanup;
        }
        if (a == b || b == c || c == a) {
            failure = "degenerate triangle with repeated vertex index";
            goto cleanup;
        }
        used[a] = used[b] = used[c] = 1;
        local_union(parent, usize, a, b);
        local_union(parent, usize, b, c);
    }
    for (size_t f = 0; f < nf; f++) {
        int32_t root = local_find(parent, faces[f*3]);
        root_faces[root]++;
    }
    for (size_t v = 0; v < nv; v++)
        if (root_faces[v]) ncomponents++;
    if (ncomponents == 0) {
        failure = "mesh has no referenced triangles";
        goto cleanup;
    }

    seed = (ComponentSeed *)checked_malloc(ncomponents, sizeof(ComponentSeed));
    report->component = (TopologyComponentInvariant *)
        checked_calloc(ncomponents, sizeof(TopologyComponentInvariant));
    root_to_component = (int32_t *)checked_malloc(nv, sizeof(int32_t));
    vertex_component = (int32_t *)checked_malloc(nv, sizeof(int32_t));
    if (!seed || !report->component || !root_to_component ||
        !vertex_component) goto cleanup;
    {
        size_t k = 0;
        for (size_t v = 0; v < nv; v++) {
            root_to_component[v] = -1;
            vertex_component[v] = -1;
            if (root_faces[v]) {
                seed[k].root = (int32_t)v;
                seed[k].faces = root_faces[v];
                k++;
            }
        }
    }
    qsort(seed, ncomponents, sizeof(ComponentSeed), cmp_component_seed);
    for (size_t c = 0; c < ncomponents; c++) {
        TopologyComponentInvariant *component = &report->component[c];
        root_to_component[seed[c].root] = (int32_t)c;
        component->rank = c;
        component->faces = seed[c].faces;
        component->has_geometry = verts != NULL;
        for (int axis = 0; axis < 3; axis++) {
            component->bbox_min[axis] = DBL_MAX;
            component->bbox_max[axis] = -DBL_MAX;
            component->centroid[axis] = 0.0;
        }
        component->root_vertex = -1;
        component->orientable = 1;
        component->input_winding_consistent = 1;
        component->beta_0 = 1;
        component->beta_1 = -1;
        component->beta_2 = -1;
        component->orientable_genus = -1;
        component->crosscap_number = -1;
    }
    for (size_t v = 0; v < nv; v++) {
        if (used[v]) {
            int32_t root = local_find(parent, (int32_t)v);
            int32_t c = root_to_component[root];
            vertex_component[v] = c;
            report->component[c].vertices++;
            if (verts) {
                TopologyComponentInvariant *component = &report->component[c];
                for (int axis = 0; axis < 3; axis++) {
                    double value = (double)verts[v*3+(size_t)axis];
                    if (value < component->bbox_min[axis])
                        component->bbox_min[axis] = value;
                    if (value > component->bbox_max[axis])
                        component->bbox_max[axis] = value;
                    component->centroid[axis] += value;
                }
            }
        } else {
            report->isolated_vertices++;
        }
    }
    if (verts) {
        for (size_t c = 0; c < ncomponents; c++) {
            TopologyComponentInvariant *component = &report->component[c];
            double denom = component->vertices ? (double)component->vertices : 1.0;
            for (int axis = 0; axis < 3; axis++)
                component->centroid[axis] /= denom;
        }
    }

    /* Build a dense unique-edge table while retaining the first two incident
     * faces and their directed half-edge orientation. */
    {
        size_t hint;
        if (nf > (SIZE_MAX - 16) / 2) goto cleanup;
        hint = nf * 2 + 16;
        if (edge_hash_init(&hash, hint) != 0) goto cleanup;
        edge_cap = nf * 2 + 16;
        edges = (AuditEdge *)checked_malloc(edge_cap, sizeof(AuditEdge));
        if (!edges) goto cleanup;
    }
    for (size_t f = 0; f < nf; f++) {
        int32_t tri[3] = {
            faces[f*3+0], faces[f*3+1], faces[f*3+2]
        };
        for (int k = 0; k < 3; k++) {
            int32_t from = tri[k], to = tri[(k+1)%3];
            int32_t lo = from < to ? from : to;
            int32_t hi = from < to ? to : from;
            int32_t index;
            int inserted = edge_hash_insert(&hash, edge_key(lo, hi),
                                             (int32_t)ne, &index);
            if (inserted < 0) goto cleanup;
            if (inserted) {
                if (ne == edge_cap) {
                    size_t next_cap = edge_cap + edge_cap / 2 + 16;
                    AuditEdge *next;
                    if (next_cap <= edge_cap ||
                        next_cap > (size_t)INT32_MAX) goto cleanup;
                    next = (AuditEdge *)realloc(edges,
                                                next_cap * sizeof(AuditEdge));
                    if (!next) goto cleanup;
                    edges = next;
                    edge_cap = next_cap;
                }
                if (ne > (size_t)INT32_MAX) goto cleanup;
                memset(&edges[ne], 0, sizeof(AuditEdge));
                edges[ne].a = lo;
                edges[ne].b = hi;
                edges[ne].f0 = (int32_t)f;
                edges[ne].f1 = -1;
                edges[ne].count = 1;
                edges[ne].dir0 = (uint8_t)(from == lo);
                ne++;
            } else {
                AuditEdge *edge = &edges[index];
                if (edge->count == 1) {
                    edge->f1 = (int32_t)f;
                    edge->dir1 = (uint8_t)(from == lo);
                }
                if (edge->count == UINT32_MAX) {
                    failure = "edge incidence counter overflow";
                    goto cleanup;
                }
                edge->count++;
            }
        }
    }
    face_edges = (int32_t *)checked_malloc(nf * 3, sizeof(int32_t));
    if (!face_edges) goto cleanup;
    for (size_t f = 0; f < nf; f++) {
        for (int k = 0; k < 3; k++) {
            int32_t a = faces[f*3+(size_t)k];
            int32_t b = faces[f*3+(size_t)((k+1)%3)];
            int32_t edge = edge_hash_find(&hash, edge_key(a, b));
            if (edge < 0) {
                failure = "internal edge-table lookup failure";
                goto cleanup;
            }
            face_edges[f*3+(size_t)k] = edge;
        }
    }
    edge_hash_dispose(&hash);

    boundary_parent = (int32_t *)checked_malloc(nv, sizeof(int32_t));
    boundary_size = (int32_t *)checked_malloc(nv, sizeof(int32_t));
    boundary_degree = (uint32_t *)checked_calloc(nv, sizeof(uint32_t));
    if (!boundary_parent || !boundary_size || !boundary_degree) goto cleanup;
    for (size_t v = 0; v < nv; v++) {
        boundary_parent[v] = (int32_t)v;
        boundary_size[v] = 1;
    }
    for (size_t e = 0; e < ne; e++) {
        AuditEdge *edge = &edges[e];
        int32_t c = vertex_component[edge->a];
        TopologyComponentInvariant *component;
        if (c < 0 || c != vertex_component[edge->b]) {
            failure = "edge crosses connected-component assignment";
            goto cleanup;
        }
        component = &report->component[c];
        component->edges++;
        if (edge->count == 1) {
            component->boundary_edges++;
            boundary_degree[edge->a]++;
            boundary_degree[edge->b]++;
            local_union(boundary_parent, boundary_size, edge->a, edge->b);
            if (component->root_vertex < 0 ||
                edge->a < component->root_vertex)
                component->root_vertex = edge->a;
            if (edge->b < component->root_vertex)
                component->root_vertex = edge->b;
        } else if (edge->count > 2) {
            component->nonmanifold_edges++;
        } else if (edge->dir0 == edge->dir1) {
            component->same_direction_edges++;
            component->input_winding_consistent = 0;
        }
    }
    for (size_t v = 0; v < nv; v++) {
        int32_t c = vertex_component[v];
        if (c < 0) continue;
        if (report->component[c].root_vertex < 0 ||
            (report->component[c].boundary_edges == 0 &&
             (int32_t)v < report->component[c].root_vertex))
            report->component[c].root_vertex = (int32_t)v;
        if (boundary_degree[v] > 0) {
            if (boundary_degree[v] != 2)
                report->component[c].boundary_irregular_vertices++;
            if (local_find(boundary_parent, (int32_t)v) == (int32_t)v)
                report->component[c].boundary_loops++;
        }
    }

    /* Materialize the boundary graph components and their physical scale.
     * This turns a bare Betti failure into an actionable distinction between
     * tiny punctures, large missing-support bays, and branched boundaries. */
    for (size_t c = 0; c < ncomponents; c++)
        report->boundary_components += report->component[c].boundary_loops;
    if (report->boundary_components) {
        const size_t exact_diameter_limit = 256;
        size_t nb = report->boundary_components;
        report->boundary = (TopologyBoundaryInvariant *)
            checked_calloc(nb, sizeof(TopologyBoundaryInvariant));
        boundary_root_to_index = (int32_t *)checked_malloc(nv, sizeof(int32_t));
        component_perimeter_index = (int32_t *)
            checked_malloc(ncomponents, sizeof(int32_t));
        if (!report->boundary || !boundary_root_to_index ||
            !component_perimeter_index) goto cleanup;
        for (size_t v = 0; v < nv; v++) boundary_root_to_index[v] = -1;
        for (size_t c = 0; c < ncomponents; c++)
            component_perimeter_index[c] = -1;

        size_t next_boundary = 0;
        for (size_t v = 0; v < nv; v++) {
            int32_t root, bi, c;
            TopologyBoundaryInvariant *boundary;
            if (boundary_degree[v] == 0) continue;
            root = local_find(boundary_parent, (int32_t)v);
            bi = boundary_root_to_index[root];
            if (bi < 0) {
                if (next_boundary >= nb) {
                    failure = "boundary-component count mismatch";
                    goto cleanup;
                }
                bi = (int32_t)next_boundary++;
                boundary_root_to_index[root] = bi;
                boundary = &report->boundary[bi];
                boundary->component = (size_t)vertex_component[v];
                boundary->root_vertex = (int32_t)v;
                boundary->has_geometry = verts != NULL;
                for (int axis = 0; axis < 3; axis++) {
                    boundary->bbox_min[axis] = DBL_MAX;
                    boundary->bbox_max[axis] = -DBL_MAX;
                }
            }
            boundary = &report->boundary[bi];
            c = vertex_component[v];
            if (c < 0 || boundary->component != (size_t)c) {
                failure = "boundary crosses connected-component assignment";
                goto cleanup;
            }
            boundary->vertices++;
            if (boundary_degree[v] != 2) boundary->irregular_vertices++;
            if ((int32_t)v < boundary->root_vertex)
                boundary->root_vertex = (int32_t)v;
            if (verts) {
                for (int axis = 0; axis < 3; axis++) {
                    double value = (double)verts[v*3+(size_t)axis];
                    if (value < boundary->bbox_min[axis])
                        boundary->bbox_min[axis] = value;
                    if (value > boundary->bbox_max[axis])
                        boundary->bbox_max[axis] = value;
                }
            }
        }
        if (next_boundary != nb) {
            failure = "boundary-component materialization mismatch";
            goto cleanup;
        }

        for (size_t e = 0; e < ne; e++) {
            AuditEdge *edge = &edges[e];
            TopologyBoundaryInvariant *boundary;
            int32_t root, bi;
            double edge_length = 1.0;
            if (edge->count != 1) continue;
            root = local_find(boundary_parent, edge->a);
            bi = boundary_root_to_index[root];
            if (bi < 0 || (size_t)bi >= nb) {
                failure = "boundary-edge component lookup failure";
                goto cleanup;
            }
            boundary = &report->boundary[bi];
            boundary->edges++;
            if (verts) {
                double dz = (double)verts[(size_t)edge->a*3+0] -
                            (double)verts[(size_t)edge->b*3+0];
                double dy = (double)verts[(size_t)edge->a*3+1] -
                            (double)verts[(size_t)edge->b*3+1];
                double dx = (double)verts[(size_t)edge->a*3+2] -
                            (double)verts[(size_t)edge->b*3+2];
                edge_length = sqrt(dz*dz + dy*dy + dx*dx);
            }
            boundary->length += edge_length;
            if (edge_length > boundary->max_edge_length)
                boundary->max_edge_length = edge_length;
        }

        boundary_vertex_offset = (size_t *)checked_calloc(nb + 1,
                                                           sizeof(size_t));
        boundary_vertex_cursor = (size_t *)checked_malloc(nb, sizeof(size_t));
        if (!boundary_vertex_offset || !boundary_vertex_cursor) goto cleanup;
        for (size_t b = 0; b < nb; b++)
            boundary_vertex_offset[b+1] =
                boundary_vertex_offset[b] + report->boundary[b].vertices;
        boundary_vertices = (int32_t *)checked_malloc(
            boundary_vertex_offset[nb], sizeof(int32_t));
        if (!boundary_vertices) goto cleanup;
        memcpy(boundary_vertex_cursor, boundary_vertex_offset,
               nb * sizeof(size_t));
        for (size_t v = 0; v < nv; v++) {
            int32_t root, bi;
            if (boundary_degree[v] == 0) continue;
            root = local_find(boundary_parent, (int32_t)v);
            bi = boundary_root_to_index[root];
            boundary_vertices[boundary_vertex_cursor[bi]++] = (int32_t)v;
        }

        for (size_t b = 0; b < nb; b++) {
            TopologyBoundaryInvariant *boundary = &report->boundary[b];
            double bbox2 = 0.0;
            boundary->simple_cycle = boundary->irregular_vertices == 0 &&
                                     boundary->edges == boundary->vertices;
            if (!verts) {
                for (int axis = 0; axis < 3; axis++)
                    boundary->bbox_min[axis] = boundary->bbox_max[axis] = 0.0;
                boundary->diameter = 0.0;
                boundary->diameter_exact = 0;
            } else {
                for (int axis = 0; axis < 3; axis++) {
                    double extent = boundary->bbox_max[axis] -
                                    boundary->bbox_min[axis];
                    bbox2 += extent * extent;
                }
                boundary->diameter = sqrt(bbox2);
                boundary->diameter_exact = 0;
                if (boundary->vertices <= exact_diameter_limit) {
                    double diameter2 = 0.0;
                    size_t begin = boundary_vertex_offset[b];
                    size_t end = boundary_vertex_offset[b+1];
                    for (size_t i = begin; i < end; i++) {
                        int32_t a = boundary_vertices[i];
                        for (size_t j = i + 1; j < end; j++) {
                            int32_t d = boundary_vertices[j];
                            double dz = (double)verts[(size_t)a*3+0] -
                                        (double)verts[(size_t)d*3+0];
                            double dy = (double)verts[(size_t)a*3+1] -
                                        (double)verts[(size_t)d*3+1];
                            double dx = (double)verts[(size_t)a*3+2] -
                                        (double)verts[(size_t)d*3+2];
                            double dist2 = dz*dz + dy*dy + dx*dx;
                            if (dist2 > diameter2) diameter2 = dist2;
                        }
                    }
                    boundary->diameter = sqrt(diameter2);
                    boundary->diameter_exact = 1;
                }
            }
            {
                size_t c = boundary->component;
                int32_t prior = component_perimeter_index[c];
                if (prior < 0 || boundary->length >
                    report->boundary[prior].length)
                    component_perimeter_index[c] = (int32_t)b;
            }
        }
        for (size_t c = 0; c < ncomponents; c++) {
            int32_t perimeter = component_perimeter_index[c];
            if (perimeter >= 0)
                report->boundary[perimeter].component_perimeter = 1;
        }
        qsort(report->boundary, nb, sizeof(*report->boundary),
              cmp_boundary_invariant);
        {
            size_t prior_component = SIZE_MAX, ordinal = 0;
            for (size_t b = 0; b < nb; b++) {
                if (report->boundary[b].component != prior_component) {
                    prior_component = report->boundary[b].component;
                    ordinal = 0;
                }
                report->boundary[b].rank = b;
                report->boundary[b].ordinal = ordinal++;
            }
        }
    }

    /* Edge checks cannot see a bowtie/pinch vertex.  Reuse the pipeline's
     * validated fan audit so the disk certificate covers both manifold axioms. */
    nonmanifold_vertex = (uint8_t *)checked_calloc(nv, sizeof(uint8_t));
    if (!nonmanifold_vertex) goto cleanup;
    manifold_arena = Arena_new();
    MeshManifold_mark_nonmanifold_vertices(manifold_arena, nv, faces, nf,
                                           nonmanifold_vertex);
    Arena_dispose(&manifold_arena);
    for (size_t v = 0; v < nv; v++) {
        if (nonmanifold_vertex[v] && vertex_component[v] >= 0)
            report->component[vertex_component[v]].nonmanifold_vertices++;
    }

    /* Intrinsic orientability: solve the face-flip parity constraints over the
     * dual graph.  Same-direction input edges merely request one local flip;
     * only a contradictory parity cycle proves non-orientability. */
    face_flip = (int8_t *)checked_malloc(nf, sizeof(int8_t));
    face_queue = (int32_t *)checked_malloc(nf, sizeof(int32_t));
    if (!face_flip || !face_queue) goto cleanup;
    memset(face_flip, 0xff, nf * sizeof(int8_t));
    for (size_t start = 0; start < nf; start++) {
        size_t qh = 0, qt = 0;
        if (face_flip[start] >= 0) continue;
        face_flip[start] = 0;
        face_queue[qt++] = (int32_t)start;
        while (qh < qt) {
            int32_t f = face_queue[qh++];
            int32_t c = vertex_component[faces[(size_t)f*3]];
            for (int k = 0; k < 3; k++) {
                AuditEdge *edge = &edges[face_edges[(size_t)f*3+(size_t)k]];
                int32_t other;
                int relation, expected;
                if (edge->count != 2) continue;
                other = edge->f0 == f ? edge->f1 : edge->f0;
                relation = edge->dir0 == edge->dir1 ? 1 : 0;
                expected = face_flip[f] ^ relation;
                if (face_flip[other] < 0) {
                    face_flip[other] = (int8_t)expected;
                    face_queue[qt++] = other;
                } else if (face_flip[other] != expected) {
                    report->component[c].orientable = 0;
                }
            }
        }
    }

    report->vertices = nv;
    report->edges = ne;
    report->faces = nf;
    report->face_components = ncomponents;
    report->topological_components = ncomponents + report->isolated_vertices;
    report->euler_characteristic = (int64_t)nv - (int64_t)ne + (int64_t)nf;
    report->beta_0 = (int64_t)report->topological_components;
    report->betti_complete = 1;
    report->all_components_are_disks =
        report->isolated_vertices == 0 ? 1 : 0;
    report->generator_basis_complete = 1;

    for (size_t c = 0; c < ncomponents; c++) {
        TopologyComponentInvariant *component = &report->component[c];
        component->euler_characteristic =
            (int64_t)component->vertices - (int64_t)component->edges +
            (int64_t)component->faces;
        component->surface_valid =
            component->nonmanifold_edges == 0 &&
            component->nonmanifold_vertices == 0 &&
            component->boundary_irregular_vertices == 0;
    }
    if (compute_exact_betti(nf, faces, face_edges, edges, ne,
                            vertex_component, report) != 0) {
        failure = "exact F2 boundary-rank calculation failed";
        goto cleanup;
    }

    for (size_t c = 0; c < ncomponents; c++) {
            TopologyComponentInvariant *component = &report->component[c];
            if (component->surface_valid) {
                int64_t expected_beta2 =
                    component->boundary_loops == 0 ? 1 : 0;
                int64_t expected_beta1 =
                    1 + expected_beta2 -
                    component->euler_characteristic;
                int64_t genus_numerator =
                    2 - (int64_t)component->boundary_loops -
                    component->euler_characteristic;
                int algebra_ok =
                    component->beta_1 == expected_beta1 &&
                    component->beta_2 == expected_beta2;
                if (component->orientable) {
                    if (genus_numerator < 0 || (genus_numerator & 1))
                        algebra_ok = 0;
                    else
                        component->orientable_genus = genus_numerator / 2;
                } else {
                    if (genus_numerator < 1) algebra_ok = 0;
                    else component->crosscap_number = genus_numerator;
                }
                if (!algebra_ok) {
                    component->surface_valid = 0;
                    component->orientable_genus =
                        component->crosscap_number = -1;
                    component->defect_mask |=
                        TOPOLOGY_DEFECT_ALGEBRA_MISMATCH;
                }
            }

            if (component->nonmanifold_edges)
                component->defect_mask |=
                    TOPOLOGY_DEFECT_NONMANIFOLD_EDGE;
            if (component->nonmanifold_vertices)
                component->defect_mask |=
                    TOPOLOGY_DEFECT_NONMANIFOLD_VERTEX;
            if (component->boundary_irregular_vertices)
                component->defect_mask |=
                    TOPOLOGY_DEFECT_IRREGULAR_BOUNDARY;
            if (!component->orientable)
                component->defect_mask |= TOPOLOGY_DEFECT_NONORIENTABLE;
            if (component->boundary_loops != 1)
                component->defect_mask |= TOPOLOGY_DEFECT_BOUNDARY_COUNT;
            if (component->beta_1 > 0)
                component->defect_mask |= TOPOLOGY_DEFECT_NONTRIVIAL_H1;
            if (component->beta_2 > 0)
                component->defect_mask |= TOPOLOGY_DEFECT_NONTRIVIAL_H2;
            if (!component->input_winding_consistent)
                component->defect_mask |=
                    TOPOLOGY_DEFECT_INCONSISTENT_WINDING;

            component->homeomorphic_to_disk =
                component->surface_valid && component->orientable &&
                component->boundary_loops == 1 &&
                component->euler_characteristic == 1 &&
                component->beta_0 == 1 &&
                component->beta_1 == 0 &&
                component->beta_2 == 0;
            if (component->homeomorphic_to_disk)
                report->disk_components++;
            else {
                report->nondisk_components++;
                report->all_components_are_disks = 0;
            }
            if (!component->surface_valid) {
                report->invalid_surface_components++;
                report->generator_basis_complete = 0;
            }
            if (!component->orientable)
                report->nonorientable_components++;
            if (!component->input_winding_consistent)
                report->inconsistent_winding_components++;
    }

    if (build_generators(verts, nv, nf, edges, ne, vertex_component,
                         options, report) != 0) {
        failure = "generator-basis construction failed";
        goto cleanup;
    }

    report->vertex_component = vertex_component;
    vertex_component = NULL;
    rc = 0;

cleanup:
    if (manifold_arena) Arena_dispose(&manifold_arena);
    edge_hash_dispose(&hash);
    free(parent); free(usize); free(used); free(root_faces);
    free(seed); free(root_to_component); free(vertex_component);
    free(edges); free(face_edges);
    free(boundary_parent); free(boundary_size); free(boundary_degree);
    free(boundary_root_to_index); free(component_perimeter_index);
    free(boundary_vertex_offset); free(boundary_vertex_cursor);
    free(boundary_vertices);
    free(nonmanifold_vertex); free(face_flip); free(face_queue);
    if (rc != 0) return fail_report(report, failure);
    return 0;
}

static int selftest_check(const char *name, int condition)
{
    fprintf(stderr, "  [topology %s] %s\n", name,
            condition ? "ok" : "FAIL");
    return condition ? 0 : 1;
}

static int selftest_disk(void)
{
    static const float v[] = {
        0,0,0,  1,0,0,  0,1,0,  1,1,0
    };
    static const int32_t f[] = { 0,1,2,  1,3,2 };
    TopologyAuditOptions options;
    TopologyAuditReport report;
    int ok;
    TopologyAudit_options_default(&options);
    options.emit_generators = 1;
    if (TopologyAudit_analyze(v, 4, f, 2, &options, &report) != 0)
        return selftest_check("disk/analyze", 0);
    ok = report.face_components == 1 && report.beta_0 == 1 &&
         report.beta_1 == 0 && report.beta_2 == 0 &&
         report.component[0].boundary_loops == 1 &&
         report.component[0].homeomorphic_to_disk &&
         report.minimal_generator_rank == 0 &&
         report.generator_basis_complete;
    TopologyAudit_dispose(&report);
    return selftest_check("disk", ok);
}

static int selftest_annulus(void)
{
    enum { NU = 10, NV = NU * 2, NF = NU * 2 };
    float v[NV * 3];
    int32_t f[NF * 3];
    TopologyAuditOptions options;
    TopologyAuditReport report;
    size_t fi = 0;
    int ok;
    for (int row = 0; row < 2; row++) {
        for (int i = 0; i < NU; i++) {
            double a = 6.2831853071795864769 * (double)i / (double)NU;
            size_t p = (size_t)(row * NU + i) * 3;
            v[p+0] = (float)cos(a);
            v[p+1] = (float)sin(a);
            v[p+2] = (float)row;
        }
    }
    for (int i = 0; i < NU; i++) {
        int n = (i + 1) % NU;
        int32_t a = i, b = n, c = NU + i, d = NU + n;
        f[fi*3+0]=a; f[fi*3+1]=b; f[fi*3+2]=c; fi++;
        f[fi*3+0]=b; f[fi*3+1]=d; f[fi*3+2]=c; fi++;
    }
    TopologyAudit_options_default(&options);
    options.emit_generators = 1;
    if (TopologyAudit_analyze(v, NV, f, NF, &options, &report) != 0)
        return selftest_check("annulus/analyze", 0);
    ok = report.beta_1 == 1 && report.beta_2 == 0 &&
         report.component[0].boundary_loops == 2 &&
         report.component[0].orientable_genus == 0 &&
         !report.component[0].homeomorphic_to_disk &&
         report.minimal_generator_rank == 1 &&
         report.emitted_generators == 1 &&
         report.generator[0].nvertices >= 3 &&
         report.generator_basis_complete;
    TopologyAudit_dispose(&report);
    return selftest_check("annulus", ok);
}

static int selftest_sphere(void)
{
    static const float v[] = {
        1,1,1,  -1,-1,1,  -1,1,-1,  1,-1,-1
    };
    static const int32_t f[] = {
        0,2,1,  0,1,3,  0,3,2,  1,2,3
    };
    TopologyAuditReport report;
    int ok;
    if (TopologyAudit_analyze(v, 4, f, 4, NULL, &report) != 0)
        return selftest_check("sphere/analyze", 0);
    ok = report.beta_0 == 1 && report.beta_1 == 0 &&
         report.beta_2 == 1 &&
         report.component[0].boundary_loops == 0 &&
         report.component[0].orientable_genus == 0 &&
         !report.component[0].homeomorphic_to_disk;
    TopologyAudit_dispose(&report);
    return selftest_check("sphere", ok);
}

static int selftest_punctured_torus(void)
{
    enum { NU = 6, NJ = 5, NV = NU * NJ, NF_FULL = NU * NJ * 2 };
    float v[NV * 3];
    int32_t f[(NF_FULL - 1) * 3];
    TopologyAuditOptions options;
    TopologyAuditReport report;
    size_t fi = 0, full_index = 0;
    int ok;
    for (int j = 0; j < NJ; j++) {
        double vv = 6.2831853071795864769 * (double)j / (double)NJ;
        for (int i = 0; i < NU; i++) {
            double uu = 6.2831853071795864769 * (double)i / (double)NU;
            double radius = 2.0 + 0.65 * cos(vv);
            size_t p = (size_t)(j * NU + i) * 3;
            v[p+0] = (float)(radius * cos(uu));
            v[p+1] = (float)(radius * sin(uu));
            v[p+2] = (float)(0.65 * sin(vv));
        }
    }
    for (int j = 0; j < NJ; j++) {
        int jn = (j + 1) % NJ;
        for (int i = 0; i < NU; i++) {
            int in = (i + 1) % NU;
            int32_t tri[2][3] = {
                { j*NU+i, j*NU+in, jn*NU+i },
                { j*NU+in, jn*NU+in, jn*NU+i }
            };
            for (int t = 0; t < 2; t++, full_index++) {
                if (full_index == 0) continue; /* one triangular boundary */
                f[fi*3+0] = tri[t][0];
                f[fi*3+1] = tri[t][1];
                f[fi*3+2] = tri[t][2];
                fi++;
            }
        }
    }
    TopologyAudit_options_default(&options);
    options.emit_generators = 1;
    if (TopologyAudit_analyze(v, NV, f, NF_FULL - 1,
                              &options, &report) != 0)
        return selftest_check("punctured-torus/analyze", 0);
    ok = report.component[0].surface_valid &&
         report.component[0].orientable &&
         report.component[0].boundary_loops == 1 &&
         report.component[0].euler_characteristic == -1 &&
         report.component[0].orientable_genus == 1 &&
         report.component[0].beta_1 == 2 &&
         report.minimal_generator_rank == 2 &&
         report.emitted_generators == 2 &&
         report.generator_basis_complete;
    TopologyAudit_dispose(&report);
    return selftest_check("punctured-torus", ok);
}

static int selftest_disconnected_disks(void)
{
    static const float v[] = {
        0,0,0, 1,0,0, 0,1,0,
        3,0,0, 4,0,0, 3,1,0
    };
    static const int32_t f[] = { 0,1,2,  3,4,5 };
    TopologyAuditReport report;
    int ok;
    if (TopologyAudit_analyze(v, 6, f, 2, NULL, &report) != 0)
        return selftest_check("two-disks/analyze", 0);
    ok = report.face_components == 2 && report.beta_0 == 2 &&
         report.beta_1 == 0 && report.beta_2 == 0 &&
         report.disk_components == 2 &&
         report.all_components_are_disks;
    TopologyAudit_dispose(&report);
    return selftest_check("two-disks", ok);
}

static int selftest_bowtie(void)
{
    static const float v[] = {
        0,0,0, 1,0,0, 0,1,0, -1,0,0, 0,-1,0
    };
    static const int32_t f[] = { 0,1,2,  0,3,4 };
    TopologyAuditReport report;
    int ok;
    if (TopologyAudit_analyze(v, 5, f, 2, NULL, &report) != 0)
        return selftest_check("bowtie/analyze", 0);
    ok = report.face_components == 1 &&
         report.invalid_surface_components == 1 &&
         report.component[0].nonmanifold_edges == 0 &&
         report.component[0].nonmanifold_vertices == 1 &&
         !report.component[0].homeomorphic_to_disk &&
         report.betti_complete &&
         report.component[0].beta_1 == 0 &&
         report.component[0].beta_2 == 0;
    TopologyAudit_dispose(&report);
    return selftest_check("bowtie", ok);
}

static int selftest_nonmanifold_rank(void)
{
    static const float v[] = { 0,0,0,  1,0,0,  0,1,0 };
    /* Three distinct 2-cells with the same boundary.  The boundary matrix has
     * three equal columns: rank(d2)=1, so beta_2=3-1=2. */
    static const int32_t f[] = {
        0,1,2,  0,1,2,  0,1,2
    };
    TopologyAuditReport report;
    int ok;
    if (TopologyAudit_analyze(v, 3, f, 3, NULL, &report) != 0)
        return selftest_check("nonmanifold-rank/analyze", 0);
    ok = report.invalid_surface_components == 1 &&
         report.component[0].nonmanifold_edges == 3 &&
         report.component[0].beta_0 == 1 &&
         report.component[0].beta_1 == 0 &&
         report.component[0].beta_2 == 2 &&
         report.beta_1 == 0 && report.beta_2 == 2 &&
         report.betti_complete;
    TopologyAudit_dispose(&report);
    return selftest_check("nonmanifold-rank", ok);
}

static int selftest_mobius(void)
{
    enum { NU = 9, NW = 3, NV = NU * NW, NF = NU * (NW - 1) * 2 };
    float v[NV * 3];
    int32_t f[NF * 3];
    TopologyAuditOptions options;
    TopologyAuditReport report;
    size_t fi = 0;
    int ok;
    for (int i = 0; i < NU; i++) {
        double u = 6.2831853071795864769 * (double)i / (double)NU;
        for (int j = 0; j < NW; j++) {
            double t = 0.7 * ((double)j / (double)(NW-1) - 0.5);
            double r = 2.0 + t * cos(0.5 * u);
            size_t p = (size_t)(i * NW + j) * 3;
            v[p+0] = (float)(r * cos(u));
            v[p+1] = (float)(r * sin(u));
            v[p+2] = (float)(t * sin(0.5 * u));
        }
    }
    for (int i = 0; i < NU; i++) {
        int seam = (i + 1 == NU);
        for (int j = 0; j < NW - 1; j++) {
            int32_t a = i*NW+j, c = i*NW+j+1;
            int32_t b = seam ? (NW-1-j) : (i+1)*NW+j;
            int32_t d = seam ? (NW-2-j) : (i+1)*NW+j+1;
            f[fi*3+0]=a; f[fi*3+1]=b; f[fi*3+2]=c; fi++;
            f[fi*3+0]=b; f[fi*3+1]=d; f[fi*3+2]=c; fi++;
        }
    }
    TopologyAudit_options_default(&options);
    options.emit_generators = 1;
    if (TopologyAudit_analyze(v, NV, f, NF, &options, &report) != 0)
        return selftest_check("mobius/analyze", 0);
    ok = report.component[0].surface_valid &&
         !report.component[0].orientable &&
         report.component[0].boundary_loops == 1 &&
         report.component[0].beta_1 == 1 &&
         report.component[0].beta_2 == 0 &&
         report.component[0].crosscap_number == 1 &&
         report.minimal_generator_rank == 1 &&
         report.emitted_generators == 1 &&
         report.generator_basis_complete;
    TopologyAudit_dispose(&report);
    return selftest_check("mobius", ok);
}

int TopologyAudit_selftest(void)
{
    int failures = 0;
    failures += selftest_disk();
    failures += selftest_annulus();
    failures += selftest_sphere();
    failures += selftest_punctured_torus();
    failures += selftest_disconnected_disks();
    failures += selftest_bowtie();
    failures += selftest_nonmanifold_rank();
    failures += selftest_mobius();
    fprintf(stderr, "  topology invariant suite: %s (%d failure%s)\n",
            failures ? "FAIL" : "PASS", failures,
            failures == 1 ? "" : "s");
    return failures;
}
