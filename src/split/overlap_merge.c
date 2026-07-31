/*
 * overlap_merge.c -- bounded delamination-heal hypothesis for OverlapSep.
 *
 * This is deliberately fail-closed.  The lifted multicut remains the source of
 * the two labels and the ordinary split remains the fallback.  A merge is only
 * constructed for one small two-label overlap event, and only survives if a
 * locally planar union re-triangulation is a single orientable manifold with
 * the same simple boundary-loop count as the input and the caller's production
 * overlap audit reports zero remaining overlaps.
 */

#include "overlap_merge.h"

#include "../common/kdtree.h"
#include "../common/csr.h"
#include "../common/mesh_manifold.h"
#include "../common/obj_io.h"
#include "../common/pca.h"
#include "../remesh/orient_mesh.h"

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Triangle library (unconstrained Delaunay; the projected source-triangle
 * union is applied as a strict post-filter). */
#define ANSI_DECLARATORS
#define REAL double
#define VOID int
#include "triangle.h"

extern jmp_buf triangle_jmpbuf;
extern int triangle_jmpbuf_set;

#define OM_MIN_PATCH_QUALITY 0.03
#define OM_MIN_FACE_AREA     1.0e-8

typedef struct {
    uint64_t key;
    int32_t  face;
    int8_t   direction;
    uint8_t  patch;
} OMEdge;

typedef struct {
    uint64_t key;
    int32_t  face;
    int32_t  label;
    int32_t  a;
    int32_t  b;
} OMSourceEdge;

typedef struct {
    double   min_x;
    double   min_y;
    double   cell;
    int32_t  nx;
    int32_t  ny;
    int32_t *offset;
    int32_t *faces;
    const double  *uv;
    const int32_t *mesh_faces;
} OMDomain;

typedef struct {
    size_t face_components;
    size_t boundary_loops;
    int    simple_boundary;
} OMShape;

static OMShape mesh_shape(Arena_T arena,
                          const int32_t *faces, size_t nf, size_t nv);

static void debug_write_mask_obj(Arena_T arena, const char *name,
                                 const float *verts, size_t nv,
                                 const int32_t *faces, size_t nf,
                                 const uint8_t *mask)
{
    const char *dir = getenv("VES_OVERLAP_DEBUG_DIR");
    if (!dir || !*dir) return;
    Arena_Mark mark = Arena_save(arena);
    uint8_t *used = (uint8_t *)ARENA_CALLOC(arena, nv, sizeof(*used));
    size_t out_nf = 0;
    for (size_t fi = 0; fi < nf; fi++) {
        if (!mask[fi]) continue;
        out_nf++;
        for (int k = 0; k < 3; k++) used[faces[fi * 3 + (size_t)k]] = 1;
    }
    if (out_nf == 0) {
        Arena_restore(arena, mark);
        return;
    }
    int32_t *remap = (int32_t *)ARENA_ALLOC(arena, nv * sizeof(*remap));
    size_t out_nv = 0;
    for (size_t vi = 0; vi < nv; vi++)
        remap[vi] = used[vi] ? (int32_t)out_nv++ : -1;
    float *out_verts = (float *)ARENA_ALLOC(
        arena, out_nv * 3 * sizeof(*out_verts));
    for (size_t vi = 0; vi < nv; vi++) {
        if (remap[vi] >= 0)
            memcpy(&out_verts[(size_t)remap[vi] * 3], &verts[vi * 3],
                   3 * sizeof(*out_verts));
    }
    int32_t *out_faces = (int32_t *)ARENA_ALLOC(
        arena, out_nf * 3 * sizeof(*out_faces));
    size_t cursor = 0;
    for (size_t fi = 0; fi < nf; fi++) {
        if (!mask[fi]) continue;
        for (int k = 0; k < 3; k++)
            out_faces[cursor * 3 + (size_t)k] =
                remap[faces[fi * 3 + (size_t)k]];
        cursor++;
    }
    char path[1024];
    snprintf(path, sizeof(path), "%s/%s", dir, name);
    ObjIO_write(path, out_verts, out_nv, out_faces, out_nf);
    fprintf(stderr, "      debug: wrote %s (nv=%zu nf=%zu)\n",
            path, out_nv, out_nf);
    Arena_restore(arena, mark);
}

static int compare_edges(const void *a, const void *b)
{
    const OMEdge *ea = (const OMEdge *)a;
    const OMEdge *eb = (const OMEdge *)b;
    if (ea->key < eb->key) return -1;
    if (ea->key > eb->key) return 1;
    if (ea->face < eb->face) return -1;
    if (ea->face > eb->face) return 1;
    return 0;
}

static int compare_source_edges(const void *a, const void *b)
{
    const OMSourceEdge *ea = (const OMSourceEdge *)a;
    const OMSourceEdge *eb = (const OMSourceEdge *)b;
    if (ea->key < eb->key) return -1;
    if (ea->key > eb->key) return 1;
    if (ea->label < eb->label) return -1;
    if (ea->label > eb->label) return 1;
    if (ea->face < eb->face) return -1;
    if (ea->face > eb->face) return 1;
    return 0;
}

static int compare_double(const void *a, const void *b)
{
    double da = *(const double *)a;
    double db = *(const double *)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

static int compare_uint64(const void *a, const void *b)
{
    uint64_t ua = *(const uint64_t *)a;
    uint64_t ub = *(const uint64_t *)b;
    if (ua < ub) return -1;
    if (ua > ub) return 1;
    return 0;
}

static void fill_edge_records(const int32_t *faces, size_t nf,
                              size_t patch_start, OMEdge *edges);

/* Lift the planar triangulation as a scalar height field over its common
 * plane.  A free 3-D Laplacian is not appropriate here: the target samples
 * come from two delaminated plies, so fixing the derived union boundary in
 * XYZ lets the solve shear the Delaunay laterally into long fans.  Keeping
 * (u,v) fixed preserves the validated planar embedding; only signed height
 * along the common normal is smoothed toward the source samples. */
static float *laplacian_heightfield_lift_patch(
                                   Arena_T arena,
                                   const float *planar,
                                   const float *targets,
                                   const float normal[3],
                                   size_t nv,
                                   const int32_t *faces, size_t nf,
                                   int iterations, float target_alpha)
{
    float *target_height = (float *)ARENA_CALLOC(
        arena, nv, sizeof(*target_height));
    float *current = (float *)ARENA_CALLOC(
        arena, nv, sizeof(*current));
    float *next = (float *)ARENA_CALLOC(arena, nv, sizeof(*next));
    float *lifted = (float *)ARENA_ALLOC(
        arena, nv * 3 * sizeof(*lifted));
    uint8_t *used = (uint8_t *)ARENA_CALLOC(arena, nv, sizeof(*used));
    for (size_t fi = 0; fi < nf; fi++)
        for (int k = 0; k < 3; k++) used[faces[fi * 3 + (size_t)k]] = 1;
    for (size_t vi = 0; vi < nv; vi++) {
        if (!used[vi]) continue;
        double dz = (double)targets[vi * 3] - planar[vi * 3];
        double dy = (double)targets[vi * 3 + 1] - planar[vi * 3 + 1];
        double dx = (double)targets[vi * 3 + 2] - planar[vi * 3 + 2];
        target_height[vi] = (float)(dz * normal[0] +
                                    dy * normal[1] +
                                    dx * normal[2]);
        current[vi] = target_height[vi];
    }

    CSR_T adjacency = CSR_from_faces(arena, faces, nf, nv);
    const int32_t *offset = CSR_offset(adjacency);
    const int32_t *neighbor = CSR_target(adjacency);
    if (iterations < 1) iterations = 1;
    if (target_alpha < 0.0f) target_alpha = 0.0f;
    if (target_alpha > 1.0f) target_alpha = 1.0f;
    for (int iteration = 0; iteration < iterations; iteration++) {
        memcpy(next, current, nv * sizeof(*next));
        for (size_t vi = 0; vi < nv; vi++) {
            if (!used[vi]) continue;
            int32_t begin = offset[vi], end = offset[vi + 1];
            int32_t degree = end - begin;
            if (degree <= 0) continue;
            double mean = 0.0;
            for (int32_t pos = begin; pos < end; pos++)
                mean += current[(size_t)neighbor[pos]];
            mean /= (double)degree;
            next[vi] = target_alpha * target_height[vi] +
                       (1.0f - target_alpha) * (float)mean;
        }
        float *swap = current;
        current = next;
        next = swap;
    }
    memcpy(lifted, planar, nv * 3 * sizeof(*lifted));
    for (size_t vi = 0; vi < nv; vi++) {
        if (!used[vi]) continue;
        for (int d = 0; d < 3; d++)
            lifted[vi * 3 + (size_t)d] += current[vi] * normal[d];
    }
    return lifted;
}

static int32_t om_find(int32_t *parent, int32_t x)
{
    int32_t root = x;
    while (parent[root] != root) root = parent[root];
    while (parent[x] != x) {
        int32_t next = parent[x];
        parent[x] = root;
        x = next;
    }
    return root;
}

static void om_union(int32_t *parent, uint8_t *rank, int32_t a, int32_t b)
{
    int32_t ra = om_find(parent, a);
    int32_t rb = om_find(parent, b);
    if (ra == rb) return;
    if (rank[ra] < rank[rb]) {
        parent[ra] = rb;
    } else {
        parent[rb] = ra;
        if (rank[ra] == rank[rb]) rank[ra]++;
    }
}

static uint64_t edge_key(int32_t a, int32_t b)
{
    uint32_t lo = (uint32_t)((a < b) ? a : b);
    uint32_t hi = (uint32_t)((a < b) ? b : a);
    return ((uint64_t)lo << 32) | (uint64_t)hi;
}

static void fill_edge_records(const int32_t *faces, size_t nf,
                              size_t patch_start, OMEdge *edges)
{
    for (size_t fi = 0; fi < nf; fi++) {
        const int32_t *tri = &faces[fi * 3];
        for (int e = 0; e < 3; e++) {
            int32_t a = tri[e];
            int32_t b = tri[(e + 1) % 3];
            OMEdge *rec = &edges[fi * 3 + (size_t)e];
            rec->key = edge_key(a, b);
            rec->face = (int32_t)fi;
            rec->direction = (int8_t)((a < b) ? 1 : -1);
            rec->patch = (uint8_t)(fi >= patch_start);
        }
    }
    qsort(edges, nf * 3, sizeof(*edges), compare_edges);
}

/*
 * The first constrained triangulation is an arrangement of both projected
 * source disks.  `keep` initially marks triangles covered by at least one
 * source disk.  Flood uncovered triangles from the convex-hull boundary;
 * every uncovered component not reached by that flood is a bounded gap in the
 * local Boolean union and belongs to the healed sheet.
 */
static size_t fill_bounded_arrangement_gaps(Arena_T arena,
                                            const int32_t *faces, size_t nf,
                                            uint8_t *keep,
                                            size_t *out_components)
{
    *out_components = 0;
    if (nf == 0) return 0;

    OMEdge *edges = (OMEdge *)ARENA_ALLOC(
        arena, nf * 3 * sizeof(*edges));
    fill_edge_records(faces, nf, nf, edges);
    int32_t *neighbor = (int32_t *)ARENA_ALLOC(
        arena, nf * 3 * sizeof(*neighbor));
    uint8_t *hull_face = (uint8_t *)ARENA_CALLOC(
        arena, nf, sizeof(*hull_face));
    for (size_t i = 0; i < nf * 3; i++) neighbor[i] = -1;

    size_t cursor = 0;
    while (cursor < nf * 3) {
        size_t end = cursor + 1;
        while (end < nf * 3 && edges[end].key == edges[cursor].key) end++;
        if (end - cursor == 1) {
            hull_face[(size_t)edges[cursor].face] = 1;
        } else {
            for (size_t a = cursor; a < end; a++) {
                int32_t fa = edges[a].face;
                for (size_t b = cursor; b < end; b++) {
                    int32_t fb = edges[b].face;
                    if (fa == fb) continue;
                    for (int slot = 0; slot < 3; slot++) {
                        if (neighbor[(size_t)fa * 3 + (size_t)slot] < 0) {
                            neighbor[(size_t)fa * 3 + (size_t)slot] = fb;
                            break;
                        }
                    }
                }
            }
        }
        cursor = end;
    }

    uint8_t *outside = (uint8_t *)ARENA_CALLOC(
        arena, nf, sizeof(*outside));
    int32_t *queue = (int32_t *)ARENA_ALLOC(
        arena, nf * sizeof(*queue));
    size_t head = 0, tail = 0;
    for (size_t fi = 0; fi < nf; fi++) {
        if (!keep[fi] && hull_face[fi]) {
            outside[fi] = 1;
            queue[tail++] = (int32_t)fi;
        }
    }
    while (head < tail) {
        int32_t fi = queue[head++];
        for (int slot = 0; slot < 3; slot++) {
            int32_t fj = neighbor[(size_t)fi * 3 + (size_t)slot];
            if (fj < 0 || keep[fj] || outside[fj]) continue;
            outside[fj] = 1;
            queue[tail++] = fj;
        }
    }

    uint8_t *seen = (uint8_t *)ARENA_CALLOC(
        arena, nf, sizeof(*seen));
    size_t filled = 0;
    for (size_t seed = 0; seed < nf; seed++) {
        if (keep[seed] || outside[seed] || seen[seed]) continue;
        (*out_components)++;
        head = 0;
        tail = 0;
        seen[seed] = 1;
        queue[tail++] = (int32_t)seed;
        while (head < tail) {
            int32_t fi = queue[head++];
            keep[fi] = 1;
            filled++;
            for (int slot = 0; slot < 3; slot++) {
                int32_t fj = neighbor[(size_t)fi * 3 + (size_t)slot];
                if (fj < 0 || keep[fj] || outside[fj] || seen[fj]) continue;
                seen[fj] = 1;
                queue[tail++] = fj;
            }
        }
    }
    return filled;
}

/* Extract the one true exterior contour of a kept face set as elementary
 * arrangement edges.  The caller has already filled bounded gaps, so a valid
 * quotient domain must be a single disk. */
static int extract_kept_boundary(Arena_T arena,
                                 const int32_t *faces, size_t nf,
                                 const uint8_t *keep, size_t nv,
                                 int32_t **out_segments,
                                 size_t *out_n_segments)
{
    size_t kept_nf = 0;
    for (size_t fi = 0; fi < nf; fi++) if (keep[fi]) kept_nf++;
    if (kept_nf == 0) return -1;

    int32_t *kept_faces = (int32_t *)ARENA_ALLOC(
        arena, kept_nf * 3 * sizeof(*kept_faces));
    size_t kf = 0;
    for (size_t fi = 0; fi < nf; fi++) {
        if (!keep[fi]) continue;
        memcpy(&kept_faces[kf * 3], &faces[fi * 3],
               3 * sizeof(*kept_faces));
        kf++;
    }
    OMShape shape = mesh_shape(arena, kept_faces, kept_nf, nv);
    if (shape.face_components != 1 || shape.boundary_loops != 1 ||
        !shape.simple_boundary) return -1;

    OMEdge *edges = (OMEdge *)ARENA_ALLOC(
        arena, kept_nf * 3 * sizeof(*edges));
    fill_edge_records(kept_faces, kept_nf, kept_nf, edges);
    size_t boundary_edges = 0;
    size_t cursor = 0;
    while (cursor < kept_nf * 3) {
        size_t end = cursor + 1;
        while (end < kept_nf * 3 && edges[end].key == edges[cursor].key) end++;
        if (end - cursor == 1) boundary_edges++;
        cursor = end;
    }
    if (boundary_edges < 3) return -1;

    int32_t *segments = (int32_t *)ARENA_ALLOC(
        arena, boundary_edges * 2 * sizeof(*segments));
    size_t si = 0;
    cursor = 0;
    while (cursor < kept_nf * 3) {
        size_t end = cursor + 1;
        while (end < kept_nf * 3 && edges[end].key == edges[cursor].key) end++;
        if (end - cursor == 1) {
            segments[si * 2] = (int32_t)(edges[cursor].key >> 32);
            segments[si * 2 + 1] =
                (int32_t)(edges[cursor].key & UINT32_C(0xffffffff));
            si++;
        }
        cursor = end;
    }
    *out_segments = segments;
    *out_n_segments = boundary_edges;
    return 0;
}

/* A planar union can switch from one source frontier to the other at a segment
 * intersection.  Until the surrounding collar is included, that switch leaves
 * a small unmatched bounded loop in the stitched 3-D mesh.  Keep the largest
 * `expected_loops` source boundary components and mark retained faces on every
 * additional loop for absorption into the next union pass. */
static int mark_extra_boundary_loop_collar(Arena_T arena,
                                           const int32_t *faces, size_t nf,
                                           size_t nv, size_t kept_nf,
                                           const int32_t *kept_source_faces,
                                           size_t expected_loops,
                                           uint8_t *out_absorb,
                                           size_t *out_loops,
                                           size_t *out_absorb_faces)
{
    OMEdge *edges = (OMEdge *)ARENA_ALLOC(
        arena, nf * 3 * sizeof(*edges));
    fill_edge_records(faces, nf, kept_nf, edges);
    int32_t *parent = (int32_t *)ARENA_ALLOC(
        arena, nv * sizeof(*parent));
    uint8_t *rank = (uint8_t *)ARENA_CALLOC(
        arena, nv, sizeof(*rank));
    for (size_t vi = 0; vi < nv; vi++) parent[vi] = -1;

    size_t cursor = 0;
    while (cursor < nf * 3) {
        size_t end = cursor + 1;
        while (end < nf * 3 && edges[end].key == edges[cursor].key) end++;
        if (end - cursor == 1) {
            int32_t a = (int32_t)(edges[cursor].key >> 32);
            int32_t b = (int32_t)(edges[cursor].key & UINT32_C(0xffffffff));
            if (parent[a] < 0) parent[a] = a;
            if (parent[b] < 0) parent[b] = b;
            om_union(parent, rank, a, b);
        }
        cursor = end;
    }

    size_t *loop_edges = (size_t *)ARENA_CALLOC(
        arena, nv, sizeof(*loop_edges));
    cursor = 0;
    while (cursor < nf * 3) {
        size_t end = cursor + 1;
        while (end < nf * 3 && edges[end].key == edges[cursor].key) end++;
        if (end - cursor == 1) {
            int32_t a = (int32_t)(edges[cursor].key >> 32);
            int32_t root = om_find(parent, a);
            loop_edges[root]++;
        }
        cursor = end;
    }
    size_t nloops = 0;
    for (size_t vi = 0; vi < nv; vi++)
        if (loop_edges[vi] > 0) nloops++;
    *out_loops = nloops;
    *out_absorb_faces = 0;
    if (nloops <= expected_loops) return 0;

    uint8_t *keep_loop = (uint8_t *)ARENA_CALLOC(
        arena, nv, sizeof(*keep_loop));
    for (size_t k = 0; k < expected_loops; k++) {
        size_t best = SIZE_MAX, best_edges = 0;
        for (size_t vi = 0; vi < nv; vi++) {
            if (!keep_loop[vi] && loop_edges[vi] > best_edges) {
                best = vi;
                best_edges = loop_edges[vi];
            }
        }
        if (best != SIZE_MAX) keep_loop[best] = 1;
    }

    cursor = 0;
    while (cursor < nf * 3) {
        size_t end = cursor + 1;
        while (end < nf * 3 && edges[end].key == edges[cursor].key) end++;
        if (end - cursor == 1) {
            int32_t a = (int32_t)(edges[cursor].key >> 32);
            int32_t root = om_find(parent, a);
            int32_t raw_face = edges[cursor].face;
            if (!keep_loop[root] && raw_face >= 0 &&
                (size_t)raw_face < kept_nf) {
                int32_t source_face = kept_source_faces[raw_face];
                if (!out_absorb[source_face]) {
                    out_absorb[source_face] = 1;
                    (*out_absorb_faces)++;
                }
            }
        }
        cursor = end;
    }
    return *out_absorb_faces > 0 ? 1 : -1;
}

static int interface_orientation(Arena_T arena,
                                 const int32_t *faces, size_t nf,
                                 size_t patch_start,
                                 size_t *out_same,
                                 size_t *out_opposite)
{
    Arena_Mark mark = Arena_save(arena);
    OMEdge *edges = (OMEdge *)ARENA_ALLOC(arena, nf * 3 * sizeof(*edges));
    fill_edge_records(faces, nf, patch_start, edges);

    size_t same = 0;
    size_t opposite = 0;
    size_t i = 0;
    while (i < nf * 3) {
        size_t j = i + 1;
        while (j < nf * 3 && edges[j].key == edges[i].key) j++;
        if (j - i == 2 && edges[i].patch != edges[i + 1].patch) {
            if (edges[i].direction == edges[i + 1].direction) same++;
            else opposite++;
        }
        i = j;
    }

    Arena_restore(arena, mark);
    *out_same = same;
    *out_opposite = opposite;
    return (same + opposite > 0) ? 0 : -1;
}

static OMShape mesh_shape(Arena_T arena,
                          const int32_t *faces, size_t nf, size_t nv)
{
    OMShape shape;
    memset(&shape, 0, sizeof(shape));
    shape.simple_boundary = 1;
    if (nf == 0 || nv == 0 || nf > (size_t)INT32_MAX ||
        nv > (size_t)INT32_MAX) {
        shape.simple_boundary = 0;
        return shape;
    }

    Arena_Mark mark = Arena_save(arena);
    OMEdge *edges = (OMEdge *)ARENA_ALLOC(arena, nf * 3 * sizeof(*edges));
    fill_edge_records(faces, nf, nf, edges);

    int32_t *face_parent = (int32_t *)ARENA_ALLOC(
        arena, nf * sizeof(*face_parent));
    uint8_t *face_rank = (uint8_t *)ARENA_CALLOC(
        arena, nf, sizeof(*face_rank));
    for (size_t fi = 0; fi < nf; fi++) face_parent[fi] = (int32_t)fi;

    int32_t *boundary_parent = (int32_t *)ARENA_ALLOC(
        arena, nv * sizeof(*boundary_parent));
    uint8_t *boundary_rank = (uint8_t *)ARENA_CALLOC(
        arena, nv, sizeof(*boundary_rank));
    int32_t *boundary_degree = (int32_t *)ARENA_CALLOC(
        arena, nv, sizeof(*boundary_degree));
    for (size_t vi = 0; vi < nv; vi++) boundary_parent[vi] = -1;

    size_t i = 0;
    while (i < nf * 3) {
        size_t j = i + 1;
        while (j < nf * 3 && edges[j].key == edges[i].key) j++;
        for (size_t k = i + 1; k < j; k++) {
            om_union(face_parent, face_rank, edges[i].face, edges[k].face);
        }
        if (j - i == 1) {
            int32_t a = (int32_t)(edges[i].key >> 32);
            int32_t b = (int32_t)(edges[i].key & 0xffffffffu);
            if (boundary_parent[a] < 0) boundary_parent[a] = a;
            if (boundary_parent[b] < 0) boundary_parent[b] = b;
            boundary_degree[a]++;
            boundary_degree[b]++;
            om_union(boundary_parent, boundary_rank, a, b);
        }
        i = j;
    }

    for (size_t fi = 0; fi < nf; fi++) {
        if (om_find(face_parent, (int32_t)fi) == (int32_t)fi)
            shape.face_components++;
    }
    for (size_t vi = 0; vi < nv; vi++) {
        if (boundary_parent[vi] < 0) continue;
        if (boundary_degree[vi] != 2) shape.simple_boundary = 0;
        if (om_find(boundary_parent, (int32_t)vi) == (int32_t)vi)
            shape.boundary_loops++;
    }

    Arena_restore(arena, mark);
    return shape;
}

static double mask_diameter(const ComponentMesh *mesh,
                            const uint8_t *face_mask)
{
    double lo[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
    double hi[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
    int any = 0;
    for (size_t fi = 0; fi < mesh->nf; fi++) {
        if (!face_mask[fi]) continue;
        for (int k = 0; k < 3; k++) {
            int32_t vi = mesh->faces[fi * 3 + (size_t)k];
            for (int d = 0; d < 3; d++) {
                double value = (double)mesh->verts[(size_t)vi * 3 + (size_t)d];
                if (value < lo[d]) lo[d] = value;
                if (value > hi[d]) hi[d] = value;
            }
            any = 1;
        }
    }
    if (!any) return 0.0;
    double d0 = hi[0] - lo[0];
    double d1 = hi[1] - lo[1];
    double d2 = hi[2] - lo[2];
    return sqrt(d0 * d0 + d1 * d1 + d2 * d2);
}

static void grow_mask_same_label(Arena_T arena,
                                 size_t nf,
                                 const int32_t *face_labels,
                                 const int32_t *adj_fa,
                                 const int32_t *adj_fb,
                                 size_t n_adj,
                                 int rings,
                                 uint8_t *mask)
{
    if (rings < 0) rings = 0;
    if (rings > 64) rings = 64;
    uint8_t *grown = (uint8_t *)ARENA_ALLOC(
        arena, nf * sizeof(*grown));
    for (int ring = 0; ring < rings; ring++) {
        memcpy(grown, mask, nf * sizeof(*grown));
        for (size_t i = 0; i < n_adj; i++) {
            int32_t fa = adj_fa[i], fb = adj_fb[i];
            if (face_labels[fa] != face_labels[fb]) continue;
            if (mask[fa]) grown[fb] = 1;
            if (mask[fb]) grown[fa] = 1;
        }
        memcpy(mask, grown, nf * sizeof(*mask));
    }
}

static OMShape mask_shape(Arena_T arena, const ComponentMesh *mesh,
                          const uint8_t *mask, size_t *out_count)
{
    size_t count = 0;
    for (size_t fi = 0; fi < mesh->nf; fi++)
        if (mask[fi]) count++;
    int32_t *faces = (int32_t *)ARENA_ALLOC(
        arena, (count > 0 ? count : 1) * 3 * sizeof(*faces));
    size_t cursor = 0;
    for (size_t fi = 0; fi < mesh->nf; fi++) {
        if (!mask[fi]) continue;
        memcpy(&faces[cursor * 3], &mesh->faces[fi * 3],
               3 * sizeof(*faces));
        cursor++;
    }
    if (out_count) *out_count = count;
    return mesh_shape(arena, faces, count, mesh->nv);
}

/* The exact overlap seed has one patch on each local ply.  Same-label ring
 * growth therefore produces two source-topological disks even though the
 * repair hypothesis is ONE healed sheet.  Join those disks through the
 * shortest bounded path in the ORIGINAL face adjacency, then (only if needed)
 * thicken the connector until the selected repair domain is one simple disk.
 *
 * More than two selected components is not a delamination merge: it is the
 * stacked/bridged multi-chart case and must remain a multicut. */
static int connect_repair_disks(Arena_T arena,
                                const ComponentMesh *mesh,
                                const int32_t *adj_fa,
                                const int32_t *adj_fb,
                                size_t n_adj,
                                int max_path_rings,
                                uint8_t *selected,
                                size_t *out_added,
                                const char **out_reason)
{
    const size_t nf = mesh->nf;
    int32_t *uf = (int32_t *)ARENA_ALLOC(arena, nf * sizeof(*uf));
    uint8_t *rank = (uint8_t *)ARENA_CALLOC(
        arena, nf, sizeof(*rank));
    for (size_t fi = 0; fi < nf; fi++)
        uf[fi] = selected[fi] ? (int32_t)fi : -1;
    for (size_t i = 0; i < n_adj; i++) {
        int32_t a = adj_fa[i], b = adj_fb[i];
        if (selected[a] && selected[b]) om_union(uf, rank, a, b);
    }

    int32_t roots[3] = {-1, -1, -1};
    size_t n_components = 0;
    for (size_t fi = 0; fi < nf; fi++) {
        if (!selected[fi]) continue;
        int32_t root = om_find(uf, (int32_t)fi);
        int seen = 0;
        for (size_t c = 0; c < n_components && c < 3; c++)
            if (roots[c] == root) seen = 1;
        if (!seen) {
            if (n_components < 3) roots[n_components] = root;
            n_components++;
        }
    }

    if (n_components > 2) {
        *out_reason = "repair-has-more-than-two-plies";
        return -1;
    }

    size_t original_count = 0;
    OMShape shape = mask_shape(arena, mesh, selected, &original_count);
    if (n_components == 1 && shape.boundary_loops == 1 &&
        shape.simple_boundary) {
        *out_added = 0;
        return 0;
    }
    if (n_components != 2) {
        *out_reason = "repair-components-invalid";
        return -1;
    }

    int32_t *degree = (int32_t *)ARENA_CALLOC(
        arena, nf + 1, sizeof(*degree));
    for (size_t i = 0; i < n_adj; i++) {
        degree[(size_t)adj_fa[i] + 1]++;
        degree[(size_t)adj_fb[i] + 1]++;
    }
    for (size_t fi = 0; fi < nf; fi++) degree[fi + 1] += degree[fi];
    int32_t *cursor = (int32_t *)ARENA_ALLOC(
        arena, nf * sizeof(*cursor));
    memcpy(cursor, degree, nf * sizeof(*cursor));
    int32_t *neighbors = (int32_t *)ARENA_ALLOC(
        arena, n_adj * 2 * sizeof(*neighbors));
    for (size_t i = 0; i < n_adj; i++) {
        int32_t a = adj_fa[i], b = adj_fb[i];
        neighbors[cursor[a]++] = b;
        neighbors[cursor[b]++] = a;
    }

    int32_t *parent = (int32_t *)ARENA_ALLOC(
        arena, nf * sizeof(*parent));
    int32_t *distance = (int32_t *)ARENA_ALLOC(
        arena, nf * sizeof(*distance));
    int32_t *queue = (int32_t *)ARENA_ALLOC(
        arena, nf * sizeof(*queue));
    for (size_t fi = 0; fi < nf; fi++) {
        parent[fi] = -1;
        distance[fi] = -1;
    }
    size_t qh = 0, qt = 0;
    for (size_t fi = 0; fi < nf; fi++) {
        if (selected[fi] && om_find(uf, (int32_t)fi) == roots[0]) {
            distance[fi] = 0;
            parent[fi] = -2;
            queue[qt++] = (int32_t)fi;
        }
    }

    int32_t last = -1;
    int found_depth = -1;
    int path_limit = max_path_rings > 0 ? max_path_rings : 12;
    while (qh < qt && last < 0) {
        int32_t face = queue[qh++];
        if (distance[face] >= path_limit) continue;
        for (int32_t e = degree[face]; e < degree[(size_t)face + 1]; e++) {
            int32_t next = neighbors[e];
            if (selected[next]) {
                if (om_find(uf, next) == roots[1]) {
                    last = face;
                    found_depth = distance[face];
                    break;
                }
                continue;
            }
            if (distance[next] >= 0) continue;
            distance[next] = distance[face] + 1;
            parent[next] = face;
            queue[qt++] = next;
        }
    }
    if (last < 0) {
        *out_reason = "repair-disks-not-locally-connected";
        return -1;
    }

    uint8_t *connector = (uint8_t *)ARENA_CALLOC(
        arena, nf, sizeof(*connector));
    size_t path_faces = 0;
    for (int32_t f = last; f >= 0 && parent[f] != -2; f = parent[f]) {
        if (!selected[f]) {
            selected[f] = 1;
            connector[f] = 1;
            path_faces++;
        }
    }

    /* A one-face-wide path can touch an existing disk at a vertex as well as
     * its intended edge, making a pinched boundary.  Thicken locally, never
     * more than three rings, until the repair domain is a true disk. */
    uint8_t *grown = (uint8_t *)ARENA_ALLOC(
        arena, nf * sizeof(*grown));
    int thickened = 0;
    for (;;) {
        shape = mask_shape(arena, mesh, selected, NULL);
        if (shape.face_components == 1 && shape.boundary_loops == 1 &&
            shape.simple_boundary) break;
        if (thickened >= 3) {
            *out_reason = "connected-repair-not-a-disk";
            return -1;
        }
        memcpy(grown, selected, nf * sizeof(*grown));
        for (size_t i = 0; i < n_adj; i++) {
            int32_t a = adj_fa[i], b = adj_fb[i];
            if (connector[a]) grown[b] = 1;
            if (connector[b]) grown[a] = 1;
        }
        for (size_t fi = 0; fi < nf; fi++)
            if (grown[fi] && !selected[fi]) connector[fi] = 1;
        memcpy(selected, grown, nf * sizeof(*selected));
        thickened++;
    }

    size_t final_count = 0;
    shape = mask_shape(arena, mesh, selected, &final_count);
    *out_added = final_count - original_count;
    fprintf(stderr,
            "      repair connector: components=2 path_depth=%d "
            "path_faces=%zu thicken=%d added=%zu -> components=%zu "
            "loops=%zu simple=%s\n",
            found_depth, path_faces, thickened, *out_added,
            shape.face_components, shape.boundary_loops,
            shape.simple_boundary ? "yes" : "no");
    if (getenv("VES_OVERLAP_DEBUG_DIR"))
        debug_write_mask_obj(arena, "overlap_merge_repair_connector.obj",
                             mesh->verts, mesh->nv,
                             mesh->faces, mesh->nf, connector);
    return 0;
}

static double mask_plane_rms(Arena_T arena,
                             const ComponentMesh *mesh,
                             const uint8_t *face_mask,
                             const float normal[3],
                             const float center[3])
{
    Arena_Mark mark = Arena_save(arena);
    uint8_t *seen = (uint8_t *)ARENA_CALLOC(
        arena, mesh->nv, sizeof(*seen));
    double sum = 0.0;
    size_t count = 0;
    for (size_t fi = 0; fi < mesh->nf; fi++) {
        if (!face_mask[fi]) continue;
        for (int k = 0; k < 3; k++) {
            int32_t vi = mesh->faces[fi * 3 + (size_t)k];
            if (seen[vi]) continue;
            seen[vi] = 1;
            const float *p = &mesh->verts[(size_t)vi * 3];
            double depth = ((double)p[0] - center[0]) * normal[0] +
                           ((double)p[1] - center[1]) * normal[1] +
                           ((double)p[2] - center[2]) * normal[2];
            sum += depth * depth;
            count++;
        }
    }
    double result = count > 0 ? sqrt(sum / (double)count) : DBL_MAX;
    Arena_restore(arena, mark);
    return result;
}

static int fit_mask_plane(Arena_T arena,
                          const ComponentMesh *mesh,
                          const uint8_t *face_mask,
                          float normal[3], float center[3],
                          double *out_rms)
{
    uint8_t *vertex_mask = (uint8_t *)ARENA_CALLOC(
        arena, mesh->nv, sizeof(*vertex_mask));
    size_t n_points = 0;
    for (size_t fi = 0; fi < mesh->nf; fi++) {
        if (!face_mask[fi]) continue;
        for (int k = 0; k < 3; k++) {
            int32_t vi = mesh->faces[fi * 3 + (size_t)k];
            if (!vertex_mask[vi]) {
                vertex_mask[vi] = 1;
                n_points++;
            }
        }
    }
    if (n_points < 3) return -1;

    float *points = (float *)ARENA_ALLOC(
        arena, n_points * 3 * sizeof(*points));
    size_t cursor = 0;
    for (size_t vi = 0; vi < mesh->nv; vi++) {
        if (!vertex_mask[vi]) continue;
        memcpy(&points[cursor * 3], &mesh->verts[vi * 3],
               3 * sizeof(float));
        cursor++;
    }
    if (PCA_normal(points, n_points, normal, center) != 0) return -1;
    double alignment = (double)normal[0] * mesh->pca_normal[0] +
                       (double)normal[1] * mesh->pca_normal[1] +
                       (double)normal[2] * mesh->pca_normal[2];
    if (alignment < 0.0) {
        normal[0] = -normal[0];
        normal[1] = -normal[1];
        normal[2] = -normal[2];
    }

    double sum = 0.0;
    for (size_t i = 0; i < n_points; i++) {
        double dz = (double)points[i * 3] - center[0];
        double dy = (double)points[i * 3 + 1] - center[1];
        double dx = (double)points[i * 3 + 2] - center[2];
        double depth = dz * normal[0] + dy * normal[1] + dx * normal[2];
        sum += depth * depth;
    }
    *out_rms = sqrt(sum / (double)n_points);
    return 0;
}

/* A local sheet can be long and narrow enough that point-PCA chooses the
 * wrong axis, especially when the two plies have some separation.  For the
 * merge construction we want the common tangent plane: average the selected
 * faces' oriented area normals (aligning the second ply by sign), while using
 * the unique selected vertices only for the plane origin and RMS diagnostic. */
static int fit_mask_tangent_plane(Arena_T arena,
                                  const ComponentMesh *mesh,
                                  const uint8_t *face_mask,
                                  float normal[3], float center[3],
                                  double *out_rms)
{
    uint8_t *seen = (uint8_t *)ARENA_CALLOC(
        arena, mesh->nv, sizeof(*seen));
    double csum[3] = {0.0, 0.0, 0.0};
    size_t n_points = 0;
    for (size_t fi = 0; fi < mesh->nf; fi++) {
        if (!face_mask[fi]) continue;
        for (int k = 0; k < 3; k++) {
            int32_t vi = mesh->faces[fi * 3 + (size_t)k];
            if (seen[vi]) continue;
            seen[vi] = 1;
            csum[0] += mesh->verts[(size_t)vi * 3];
            csum[1] += mesh->verts[(size_t)vi * 3 + 1];
            csum[2] += mesh->verts[(size_t)vi * 3 + 2];
            n_points++;
        }
    }
    if (n_points < 3) return -1;
    center[0] = (float)(csum[0] / (double)n_points);
    center[1] = (float)(csum[1] / (double)n_points);
    center[2] = (float)(csum[2] / (double)n_points);

    double reference[3] = {0.0, 0.0, 0.0};
    int have_reference = 0;
    double sum[3] = {0.0, 0.0, 0.0};
    for (size_t fi = 0; fi < mesh->nf; fi++) {
        if (!face_mask[fi]) continue;
        const int32_t *tri = &mesh->faces[fi * 3];
        const float *a = &mesh->verts[(size_t)tri[0] * 3];
        const float *b = &mesh->verts[(size_t)tri[1] * 3];
        const float *c = &mesh->verts[(size_t)tri[2] * 3];
        double ab[3] = {(double)b[0] - a[0], (double)b[1] - a[1],
                        (double)b[2] - a[2]};
        double ac[3] = {(double)c[0] - a[0], (double)c[1] - a[1],
                        (double)c[2] - a[2]};
        double cross[3] = {
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0]
        };
        double length = sqrt(cross[0] * cross[0] + cross[1] * cross[1] +
                             cross[2] * cross[2]);
        if (length <= 1.0e-12) continue;
        if (!have_reference) {
            reference[0] = cross[0];
            reference[1] = cross[1];
            reference[2] = cross[2];
            have_reference = 1;
        }
        double alignment = cross[0] * reference[0] +
                           cross[1] * reference[1] +
                           cross[2] * reference[2];
        if (alignment < 0.0) {
            cross[0] = -cross[0];
            cross[1] = -cross[1];
            cross[2] = -cross[2];
        }
        sum[0] += cross[0];
        sum[1] += cross[1];
        sum[2] += cross[2];
    }
    double length = sqrt(sum[0] * sum[0] + sum[1] * sum[1] +
                         sum[2] * sum[2]);
    if (!have_reference || length <= 1.0e-12) return -1;
    normal[0] = (float)(sum[0] / length);
    normal[1] = (float)(sum[1] / length);
    normal[2] = (float)(sum[2] / length);
    *out_rms = mask_plane_rms(arena, mesh, face_mask, normal, center);
    return 0;
}

static int om_faces_share_vertex(const int32_t *faces,
                                 int32_t fa, int32_t fb)
{
    const int32_t *a = &faces[(size_t)fa * 3];
    const int32_t *b = &faces[(size_t)fb * 3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (a[i] == b[j]) return 1;
        }
    }
    return 0;
}

static int triangles_overlap_strict(const double *uv,
                                    const int32_t *faces,
                                    int32_t fa, int32_t fb)
{
    int32_t ids[2][3] = {
        {faces[(size_t)fa * 3], faces[(size_t)fa * 3 + 1],
         faces[(size_t)fa * 3 + 2]},
        {faces[(size_t)fb * 3], faces[(size_t)fb * 3 + 1],
         faces[(size_t)fb * 3 + 2]}
    };
    for (int tri = 0; tri < 2; tri++) {
        for (int edge = 0; edge < 3; edge++) {
            const double *p = &uv[(size_t)ids[tri][edge] * 2];
            const double *q = &uv[(size_t)ids[tri][(edge + 1) % 3] * 2];
            double axis_x = -(q[1] - p[1]);
            double axis_y = q[0] - p[0];
            double length = sqrt(axis_x * axis_x + axis_y * axis_y);
            if (length <= 1.0e-12) continue;
            axis_x /= length;
            axis_y /= length;

            double min_a = DBL_MAX, max_a = -DBL_MAX;
            double min_b = DBL_MAX, max_b = -DBL_MAX;
            for (int k = 0; k < 3; k++) {
                const double *a = &uv[(size_t)ids[0][k] * 2];
                const double *b = &uv[(size_t)ids[1][k] * 2];
                double pa = a[0] * axis_x + a[1] * axis_y;
                double pb = b[0] * axis_x + b[1] * axis_y;
                if (pa < min_a) min_a = pa;
                if (pa > max_a) max_a = pa;
                if (pb < min_b) min_b = pb;
                if (pb > max_b) max_b = pb;
            }
            const double eps = 1.0e-8;
            if (max_a <= min_b + eps || max_b <= min_a + eps) return 0;
        }
    }
    return 1;
}

static int face_near_vertex_count(const ComponentMesh *mesh,
                                  int32_t face,
                                  KDTree_T other_tree,
                                  float max_distance_sq)
{
    int near_count = 0;
    for (int corner = 0; corner < 3; corner++) {
        int32_t vi = mesh->faces[(size_t)face * 3 + (size_t)corner];
        const float *p = &mesh->verts[(size_t)vi * 3];
        float distance_sq = FLT_MAX;
        (void)KDTree_nearest(other_tree, p, &distance_sq);
        if (distance_sq <= max_distance_sq) near_count++;
    }
    return near_count;
}

/* The production multicut seed is defined in the component's global PCA
 * projection.  That is an excellent detector, but a poor plane estimator for
 * a narrow delamination: its repulsive faces need not be co-located with the
 * physical seam between the two resulting charts.  Once multicut has supplied
 * exactly two labels and one bounded event, recover their global physical
 * contact as faces whose three vertices lie close to the other label in 3-D.
 * The face-count and diameter gates below make a broad stack fail closed. */
static int find_physical_contact_seed(Arena_T arena,
                                      const ComponentMesh *mesh,
                                      const int32_t *face_labels,
                                      int32_t label_a, int32_t label_b,
                                      double contact_distance,
                                      uint8_t *seed,
                                      size_t *out_faces_a,
                                      size_t *out_faces_b)
{
    uint8_t *vertex_a = (uint8_t *)ARENA_CALLOC(
        arena, mesh->nv, sizeof(*vertex_a));
    uint8_t *vertex_b = (uint8_t *)ARENA_CALLOC(
        arena, mesh->nv, sizeof(*vertex_b));
    for (size_t fi = 0; fi < mesh->nf; fi++) {
        uint8_t *mask = face_labels[fi] == label_a ? vertex_a :
                        face_labels[fi] == label_b ? vertex_b : NULL;
        if (!mask) continue;
        for (int corner = 0; corner < 3; corner++)
            mask[mesh->faces[fi * 3 + (size_t)corner]] = 1;
    }

    size_t n_vertices_a = 0, n_vertices_b = 0;
    for (size_t vi = 0; vi < mesh->nv; vi++) {
        if (vertex_a[vi]) n_vertices_a++;
        if (vertex_b[vi]) n_vertices_b++;
    }
    if (n_vertices_a == 0 || n_vertices_b == 0) return -1;
    float *points_a = (float *)ARENA_ALLOC(
        arena, n_vertices_a * 3 * sizeof(*points_a));
    float *points_b = (float *)ARENA_ALLOC(
        arena, n_vertices_b * 3 * sizeof(*points_b));
    size_t cursor_a = 0, cursor_b = 0;
    for (size_t vi = 0; vi < mesh->nv; vi++) {
        if (vertex_a[vi]) {
            memcpy(&points_a[cursor_a * 3], &mesh->verts[vi * 3],
                   3 * sizeof(*points_a));
            cursor_a++;
        }
        if (vertex_b[vi]) {
            memcpy(&points_b[cursor_b * 3], &mesh->verts[vi * 3],
                   3 * sizeof(*points_b));
            cursor_b++;
        }
    }
    KDTree_T tree_a = KDTree_new(arena, points_a, n_vertices_a);
    KDTree_T tree_b = KDTree_new(arena, points_b, n_vertices_b);

    memset(seed, 0, mesh->nf * sizeof(*seed));
    float max_distance_sq = (float)(contact_distance * contact_distance);
    size_t faces_a = 0, faces_b = 0;
    size_t search_faces_a = 0, search_faces_b = 0;
    size_t near_hist_a[4] = {0, 0, 0, 0};
    size_t near_hist_b[4] = {0, 0, 0, 0};
    for (size_t fi = 0; fi < mesh->nf; fi++) {
        if (face_labels[fi] == label_a) {
            int near = face_near_vertex_count(
                mesh, (int32_t)fi, tree_b, max_distance_sq);
            search_faces_a++;
            near_hist_a[near]++;
            if (near == 3) {
                seed[fi] = 1;
                faces_a++;
            }
        } else if (face_labels[fi] == label_b) {
            int near = face_near_vertex_count(
                mesh, (int32_t)fi, tree_a, max_distance_sq);
            search_faces_b++;
            near_hist_b[near]++;
            if (near == 3) {
                seed[fi] = 1;
                faces_b++;
            }
        }
    }
    fprintf(stderr,
            "      contact search: global faces=%zu+%zu verts=%zu+%zu "
            "near[0..3]=%zu/%zu/%zu/%zu + %zu/%zu/%zu/%zu\n",
            search_faces_a, search_faces_b, n_vertices_a, n_vertices_b,
            near_hist_a[0], near_hist_a[1], near_hist_a[2], near_hist_a[3],
            near_hist_b[0], near_hist_b[1], near_hist_b[2], near_hist_b[3]);
    *out_faces_a = faces_a;
    *out_faces_b = faces_b;
    return faces_a >= 8 && faces_b >= 8 ? 0 : -1;
}

static double orient2d(const double *a, const double *b, const double *c)
{
    return (b[0] - a[0]) * (c[1] - a[1]) -
           (b[1] - a[1]) * (c[0] - a[0]);
}

static int point_in_triangle(const double *p,
                             const double *a,
                             const double *b,
                             const double *c)
{
    double scale = 1.0;
    double values[6] = {a[0], a[1], b[0], b[1], c[0], c[1]};
    for (int i = 0; i < 6; i++) {
        double av = fabs(values[i]);
        if (av > scale) scale = av;
    }
    double eps = 1.0e-9 * scale * scale;
    double o0 = orient2d(a, b, p);
    double o1 = orient2d(b, c, p);
    double o2 = orient2d(c, a, p);
    int has_neg = (o0 < -eps) || (o1 < -eps) || (o2 < -eps);
    int has_pos = (o0 > eps) || (o1 > eps) || (o2 > eps);
    return !(has_neg && has_pos);
}

static int domain_build(Arena_T arena,
                        const double *uv,
                        const int32_t *faces, size_t nf,
                        const uint8_t *selected,
                        OMDomain *domain)
{
    memset(domain, 0, sizeof(*domain));
    double min_x = DBL_MAX, min_y = DBL_MAX;
    double max_x = -DBL_MAX, max_y = -DBL_MAX;
    double edge_sum = 0.0;
    size_t edge_count = 0;
    size_t selected_count = 0;

    for (size_t fi = 0; fi < nf; fi++) {
        if (!selected[fi]) continue;
        selected_count++;
        const int32_t *tri = &faces[fi * 3];
        for (int k = 0; k < 3; k++) {
            int32_t a = tri[k];
            int32_t b = tri[(k + 1) % 3];
            double x = uv[(size_t)a * 2];
            double y = uv[(size_t)a * 2 + 1];
            double dx = x - uv[(size_t)b * 2];
            double dy = y - uv[(size_t)b * 2 + 1];
            edge_sum += sqrt(dx * dx + dy * dy);
            edge_count++;
            if (x < min_x) min_x = x;
            if (x > max_x) max_x = x;
            if (y < min_y) min_y = y;
            if (y > max_y) max_y = y;
        }
    }
    if (selected_count == 0 || edge_count == 0) return -1;

    double cell = edge_sum / (double)edge_count;
    if (cell < 0.25) cell = 0.25;
    int32_t nx = (int32_t)ceil((max_x - min_x) / cell) + 1;
    int32_t ny = (int32_t)ceil((max_y - min_y) / cell) + 1;
    if (nx < 1) nx = 1;
    if (ny < 1) ny = 1;
    if (nx > 2048 || ny > 2048) {
        double span_x = max_x - min_x;
        double span_y = max_y - min_y;
        double need_x = span_x / 2047.0;
        double need_y = span_y / 2047.0;
        if (need_x > cell) cell = need_x;
        if (need_y > cell) cell = need_y;
        nx = (int32_t)ceil(span_x / cell) + 1;
        ny = (int32_t)ceil(span_y / cell) + 1;
    }
    if (nx < 1 || ny < 1 || nx > 2048 || ny > 2048) return -1;

    size_t n_cells = (size_t)nx * (size_t)ny;
    int32_t *count = (int32_t *)ARENA_CALLOC(
        arena, n_cells, sizeof(*count));

    for (size_t fi = 0; fi < nf; fi++) {
        if (!selected[fi]) continue;
        double fmin_x = DBL_MAX, fmin_y = DBL_MAX;
        double fmax_x = -DBL_MAX, fmax_y = -DBL_MAX;
        for (int k = 0; k < 3; k++) {
            int32_t vi = faces[fi * 3 + (size_t)k];
            double x = uv[(size_t)vi * 2];
            double y = uv[(size_t)vi * 2 + 1];
            if (x < fmin_x) fmin_x = x;
            if (x > fmax_x) fmax_x = x;
            if (y < fmin_y) fmin_y = y;
            if (y > fmax_y) fmax_y = y;
        }
        int32_t ix0 = (int32_t)floor((fmin_x - min_x) / cell);
        int32_t ix1 = (int32_t)floor((fmax_x - min_x) / cell);
        int32_t iy0 = (int32_t)floor((fmin_y - min_y) / cell);
        int32_t iy1 = (int32_t)floor((fmax_y - min_y) / cell);
        if (ix0 < 0) ix0 = 0;
        if (iy0 < 0) iy0 = 0;
        if (ix1 >= nx) ix1 = nx - 1;
        if (iy1 >= ny) iy1 = ny - 1;
        for (int32_t iy = iy0; iy <= iy1; iy++) {
            for (int32_t ix = ix0; ix <= ix1; ix++) {
                size_t ci = (size_t)iy * (size_t)nx + (size_t)ix;
                if (count[ci] == INT32_MAX) return -1;
                count[ci]++;
            }
        }
    }

    int32_t *offset = (int32_t *)ARENA_ALLOC(
        arena, (n_cells + 1) * sizeof(*offset));
    offset[0] = 0;
    for (size_t ci = 0; ci < n_cells; ci++) {
        if (offset[ci] > INT32_MAX - count[ci]) return -1;
        offset[ci + 1] = offset[ci] + count[ci];
    }
    int32_t total = offset[n_cells];
    int32_t *cell_faces = (int32_t *)ARENA_ALLOC(
        arena, (size_t)total * sizeof(*cell_faces));
    int32_t *cursor = (int32_t *)ARENA_CALLOC(
        arena, n_cells, sizeof(*cursor));

    for (size_t fi = 0; fi < nf; fi++) {
        if (!selected[fi]) continue;
        double fmin_x = DBL_MAX, fmin_y = DBL_MAX;
        double fmax_x = -DBL_MAX, fmax_y = -DBL_MAX;
        for (int k = 0; k < 3; k++) {
            int32_t vi = faces[fi * 3 + (size_t)k];
            double x = uv[(size_t)vi * 2];
            double y = uv[(size_t)vi * 2 + 1];
            if (x < fmin_x) fmin_x = x;
            if (x > fmax_x) fmax_x = x;
            if (y < fmin_y) fmin_y = y;
            if (y > fmax_y) fmax_y = y;
        }
        int32_t ix0 = (int32_t)floor((fmin_x - min_x) / cell);
        int32_t ix1 = (int32_t)floor((fmax_x - min_x) / cell);
        int32_t iy0 = (int32_t)floor((fmin_y - min_y) / cell);
        int32_t iy1 = (int32_t)floor((fmax_y - min_y) / cell);
        if (ix0 < 0) ix0 = 0;
        if (iy0 < 0) iy0 = 0;
        if (ix1 >= nx) ix1 = nx - 1;
        if (iy1 >= ny) iy1 = ny - 1;
        for (int32_t iy = iy0; iy <= iy1; iy++) {
            for (int32_t ix = ix0; ix <= ix1; ix++) {
                size_t ci = (size_t)iy * (size_t)nx + (size_t)ix;
                cell_faces[offset[ci] + cursor[ci]++] = (int32_t)fi;
            }
        }
    }

    domain->min_x = min_x;
    domain->min_y = min_y;
    domain->cell = cell;
    domain->nx = nx;
    domain->ny = ny;
    domain->offset = offset;
    domain->faces = cell_faces;
    domain->uv = uv;
    domain->mesh_faces = faces;
    return 0;
}

static int domain_contains(const OMDomain *domain, const double p[2])
{
    int32_t ix = (int32_t)floor((p[0] - domain->min_x) / domain->cell);
    int32_t iy = (int32_t)floor((p[1] - domain->min_y) / domain->cell);
    if (ix < 0 || iy < 0 || ix >= domain->nx || iy >= domain->ny)
        return 0;
    size_t ci = (size_t)iy * (size_t)domain->nx + (size_t)ix;
    int32_t begin = domain->offset[ci];
    int32_t end = domain->offset[ci + 1];
    for (int32_t k = begin; k < end; k++) {
        int32_t fi = domain->faces[k];
        const int32_t *tri = &domain->mesh_faces[(size_t)fi * 3];
        const double *a = &domain->uv[(size_t)tri[0] * 2];
        const double *b = &domain->uv[(size_t)tri[1] * 2];
        const double *c = &domain->uv[(size_t)tri[2] * 2];
        if (point_in_triangle(p, a, b, c)) return 1;
    }
    return 0;
}

/* Mod-2 coverage of the projected selected source patch.  Internal source
 * edges cancel in pairs; its constrained topological frontier is therefore
 * exactly the boundary across which this parity changes.  A delaminated
 * double-covered region becomes a hole rather than a second coincident ply. */
static int domain_coverage(const OMDomain *domain, const double p[2])
{
    int32_t ix = (int32_t)floor((p[0] - domain->min_x) / domain->cell);
    int32_t iy = (int32_t)floor((p[1] - domain->min_y) / domain->cell);
    if (ix < 0 || iy < 0 || ix >= domain->nx || iy >= domain->ny) return 0;
    int32_t cell = iy * domain->nx + ix;
    int count = 0;
    for (int32_t pos = domain->offset[cell];
         pos < domain->offset[cell + 1]; pos++) {
        int32_t fi = domain->faces[pos];
        const int32_t *tri = &domain->mesh_faces[(size_t)fi * 3];
        const double *a = &domain->uv[(size_t)tri[0] * 2];
        const double *b = &domain->uv[(size_t)tri[1] * 2];
        const double *c = &domain->uv[(size_t)tri[2] * 2];
        if (point_in_triangle(p, a, b, c)) count++;
    }
    return count;
}

static void free_triangle_output(struct triangulateio *out)
{
    if (out->pointlist) trifree((VOID *)out->pointlist);
    if (out->pointattributelist) trifree((VOID *)out->pointattributelist);
    if (out->pointmarkerlist) trifree((VOID *)out->pointmarkerlist);
    if (out->trianglelist) trifree((VOID *)out->trianglelist);
    if (out->triangleattributelist)
        trifree((VOID *)out->triangleattributelist);
    if (out->neighborlist) trifree((VOID *)out->neighborlist);
    if (out->segmentlist) trifree((VOID *)out->segmentlist);
    if (out->segmentmarkerlist) trifree((VOID *)out->segmentmarkerlist);
    if (out->edgelist) trifree((VOID *)out->edgelist);
    if (out->edgemarkerlist) trifree((VOID *)out->edgemarkerlist);
    if (out->normlist) trifree((VOID *)out->normlist);
}

static int delaunay_points(Arena_T arena,
                           const double *points, size_t n_points,
                           const int32_t *segments, size_t n_segments,
                           double max_area,
                           double **out_points, size_t *out_n_points,
                           int32_t **out_faces, size_t *out_nf)
{
    *out_points = NULL;
    *out_n_points = 0;
    *out_faces = NULL;
    *out_nf = 0;
    if (n_points < 3 || n_points > (size_t)INT_MAX) return -1;

    struct triangulateio input;
    struct triangulateio output;
    memset(&input, 0, sizeof(input));
    memset(&output, 0, sizeof(output));
    input.numberofpoints = (int)n_points;
    input.pointlist = (REAL *)malloc(n_points * 2 * sizeof(REAL));
    if (!input.pointlist) return -1;
    memcpy(input.pointlist, points, n_points * 2 * sizeof(REAL));
    if (n_segments > 0) {
        if (n_segments > (size_t)INT_MAX) {
            free(input.pointlist);
            return -1;
        }
        input.numberofsegments = (int)n_segments;
        input.segmentlist = (int *)malloc(n_segments * 2 * sizeof(int));
        if (!input.segmentlist) {
            free(input.pointlist);
            return -1;
        }
        for (size_t i = 0; i < n_segments * 2; i++)
            input.segmentlist[i] = (int)segments[i];
    }

    triangle_jmpbuf_set = 1;
    if (setjmp(triangle_jmpbuf) != 0) {
        triangle_jmpbuf_set = 0;
        free(input.pointlist);
        free(input.segmentlist);
        return -1;
    }
    /* p=PSLG, Y=no discretionary boundary splitting, z=zero based.  The first
     * arrangement pass is unrefined; the final healed disk gets an interior
     * area cap so sparse union samples cannot form multi-voxel collar fans. */
    char switches[64];
    if (max_area > 0.0)
        snprintf(switches, sizeof(switches), "pYq20a%.9gzQ", max_area);
    else
        snprintf(switches, sizeof(switches), "pYzQ");
    triangulate(switches, &input, &output, NULL);
    triangle_jmpbuf_set = 0;
    free(input.pointlist);
    free(input.segmentlist);

    fprintf(stderr,
            "      constrained Delaunay: input=%zu points/%zu segments "
            "output=%d points/%d triangles/%d segments\n",
            n_points, n_segments, output.numberofpoints,
            output.numberoftriangles, output.numberofsegments);

    if (output.numberoftriangles <= 0 || !output.trianglelist ||
        output.numberofpoints < (int)n_points || !output.pointlist) {
        free_triangle_output(&output);
        return -1;
    }

    size_t result_nv = (size_t)output.numberofpoints;
    double *result_points = (double *)ARENA_ALLOC(
        arena, result_nv * 2 * sizeof(*result_points));
    memcpy(result_points, output.pointlist,
           result_nv * 2 * sizeof(*result_points));
    size_t nf = (size_t)output.numberoftriangles;
    int32_t *faces = (int32_t *)ARENA_ALLOC(
        arena, nf * 3 * sizeof(*faces));
    for (size_t fi = 0; fi < nf; fi++) {
        for (int k = 0; k < 3; k++) {
            int index = output.trianglelist[fi * 3 + (size_t)k];
            if (index < 0 || (size_t)index >= result_nv) {
                free_triangle_output(&output);
                return -1;
            }
            faces[fi * 3 + (size_t)k] = (int32_t)index;
        }
    }
    free_triangle_output(&output);
    *out_points = result_points;
    *out_n_points = result_nv;
    *out_faces = faces;
    *out_nf = nf;
    return 0;
}

static double triangle_quality(const float *verts, const int32_t tri[3],
                               double *out_area)
{
    const float *a = &verts[(size_t)tri[0] * 3];
    const float *b = &verts[(size_t)tri[1] * 3];
    const float *c = &verts[(size_t)tri[2] * 3];
    double ab[3], ac[3], bc[3];
    for (int d = 0; d < 3; d++) {
        ab[d] = (double)b[d] - (double)a[d];
        ac[d] = (double)c[d] - (double)a[d];
        bc[d] = (double)c[d] - (double)b[d];
    }
    double cross[3] = {
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0]
    };
    double area = 0.5 * sqrt(cross[0] * cross[0] +
                             cross[1] * cross[1] +
                             cross[2] * cross[2]);
    double sum_l2 = 0.0;
    for (int d = 0; d < 3; d++) {
        sum_l2 += ab[d] * ab[d] + ac[d] * ac[d] + bc[d] * bc[d];
    }
    *out_area = area;
    if (sum_l2 <= 0.0) return 0.0;
    return 4.0 * sqrt(3.0) * area / sum_l2;
}

/* Construct an injective chart for the source-topological repair disk.
 * Orthogonal projection is deliberately not used as the domain boundary: a
 * delamination can make a perfectly simple source loop cross itself in that
 * projection.  The loop is mapped to a convex ellipse with source edge-length
 * spacing, and the interior receives a uniform Tutte (harmonic) embedding.
 * These coordinates define topology only; the loop's real 3-D positions stay
 * pinned when the replacement is stitched and fitted to the volume. */
static int parameterize_repair_disk(Arena_T arena,
                                    const ComponentMesh *mesh,
                                    const uint8_t *selected,
                                    size_t selected_count,
                                    const double *projected_uv,
                                    double *chart_uv)
{
    const size_t nv = mesh->nv;
    int32_t *disk_faces = (int32_t *)ARENA_ALLOC(
        arena, selected_count * 3 * sizeof(*disk_faces));
    size_t disk_nf = 0;
    uint8_t *disk_vertex = (uint8_t *)ARENA_CALLOC(
        arena, nv, sizeof(*disk_vertex));
    for (size_t fi = 0; fi < mesh->nf; fi++) {
        if (!selected[fi]) continue;
        memcpy(&disk_faces[disk_nf * 3], &mesh->faces[fi * 3],
               3 * sizeof(*disk_faces));
        for (int k = 0; k < 3; k++) disk_vertex[mesh->faces[fi * 3 + k]] = 1;
        disk_nf++;
    }
    if (disk_nf != selected_count) return -1;

    OMShape disk_shape = mesh_shape(arena, disk_faces, disk_nf, nv);
    fprintf(stderr,
            "      repair source topology: components=%zu loops=%zu "
            "simple=%s\n",
            disk_shape.face_components, disk_shape.boundary_loops,
            disk_shape.simple_boundary ? "yes" : "no");
    if (disk_shape.face_components != 1 ||
        disk_shape.boundary_loops != 1 || !disk_shape.simple_boundary)
        return -1;

    OMEdge *edges = (OMEdge *)ARENA_ALLOC(
        arena, disk_nf * 3 * sizeof(*edges));
    fill_edge_records(disk_faces, disk_nf, disk_nf, edges);
    int32_t *neighbor0 = (int32_t *)ARENA_ALLOC(
        arena, nv * sizeof(*neighbor0));
    int32_t *neighbor1 = (int32_t *)ARENA_ALLOC(
        arena, nv * sizeof(*neighbor1));
    uint8_t *boundary_degree = (uint8_t *)ARENA_CALLOC(
        arena, nv, sizeof(*boundary_degree));
    for (size_t vi = 0; vi < nv; vi++) {
        neighbor0[vi] = -1;
        neighbor1[vi] = -1;
    }

    size_t boundary_edges = 0;
    size_t cursor = 0;
    int32_t start = -1;
    while (cursor < disk_nf * 3) {
        size_t end = cursor + 1;
        while (end < disk_nf * 3 && edges[end].key == edges[cursor].key)
            end++;
        if (end - cursor == 1) {
            int32_t a = (int32_t)(edges[cursor].key >> 32);
            int32_t b = (int32_t)(edges[cursor].key & UINT32_C(0xffffffff));
            if (boundary_degree[a] >= 2 || boundary_degree[b] >= 2)
                return -1;
            if (boundary_degree[a]++ == 0) neighbor0[a] = b;
            else neighbor1[a] = b;
            if (boundary_degree[b]++ == 0) neighbor0[b] = a;
            else neighbor1[b] = a;
            if (start < 0) start = a;
            boundary_edges++;
        }
        cursor = end;
    }
    if (boundary_edges < 3 || start < 0) return -1;
    for (size_t vi = 0; vi < nv; vi++) {
        if (boundary_degree[vi] != 0 && boundary_degree[vi] != 2)
            return -1;
    }

    int32_t *loop = (int32_t *)ARENA_ALLOC(
        arena, boundary_edges * sizeof(*loop));
    int32_t previous = -1;
    int32_t current = start;
    for (size_t i = 0; i < boundary_edges; i++) {
        if (i > 0 && current == start) return -1;
        loop[i] = current;
        int32_t next = neighbor0[current] != previous ?
            neighbor0[current] : neighbor1[current];
        if (next < 0) return -1;
        previous = current;
        current = next;
    }
    if (current != start) return -1;

    double *arc = (double *)ARENA_ALLOC(
        arena, (boundary_edges + 1) * sizeof(*arc));
    arc[0] = 0.0;
    double center_uv[2] = {0.0, 0.0};
    for (size_t i = 0; i < boundary_edges; i++) {
        int32_t a = loop[i];
        int32_t b = loop[(i + 1) % boundary_edges];
        const float *pa = &mesh->verts[(size_t)a * 3];
        const float *pb = &mesh->verts[(size_t)b * 3];
        double d0 = (double)pb[0] - pa[0];
        double d1 = (double)pb[1] - pa[1];
        double d2 = (double)pb[2] - pa[2];
        arc[i + 1] = arc[i] + sqrt(d0 * d0 + d1 * d1 + d2 * d2);
        center_uv[0] += projected_uv[(size_t)a * 2];
        center_uv[1] += projected_uv[(size_t)a * 2 + 1];
    }
    double perimeter = arc[boundary_edges];
    if (perimeter <= 1.0e-6) return -1;
    center_uv[0] /= (double)boundary_edges;
    center_uv[1] /= (double)boundary_edges;

    double cov00 = 0.0, cov01 = 0.0, cov11 = 0.0;
    for (size_t i = 0; i < boundary_edges; i++) {
        int32_t vi = loop[i];
        double x = projected_uv[(size_t)vi * 2] - center_uv[0];
        double y = projected_uv[(size_t)vi * 2 + 1] - center_uv[1];
        cov00 += x * x;
        cov01 += x * y;
        cov11 += y * y;
    }
    cov00 /= (double)boundary_edges;
    cov01 /= (double)boundary_edges;
    cov11 /= (double)boundary_edges;
    double trace = cov00 + cov11;
    double disc = sqrt(fmax(0.0,
        (cov00 - cov11) * (cov00 - cov11) + 4.0 * cov01 * cov01));
    double lambda0 = 0.5 * (trace + disc);
    double lambda1 = 0.5 * (trace - disc);
    double axis0[2];
    if (fabs(cov01) > 1.0e-12) {
        axis0[0] = cov01;
        axis0[1] = lambda0 - cov00;
    } else if (cov00 >= cov11) {
        axis0[0] = 1.0; axis0[1] = 0.0;
    } else {
        axis0[0] = 0.0; axis0[1] = 1.0;
    }
    double axis_length = sqrt(axis0[0] * axis0[0] + axis0[1] * axis0[1]);
    if (axis_length <= 1.0e-12) return -1;
    axis0[0] /= axis_length;
    axis0[1] /= axis_length;
    double axis1[2] = {-axis0[1], axis0[0]};

    const double pi = 3.14159265358979323846;
    double base_radius = perimeter / (2.0 * pi);
    double radius0 = sqrt(fmax(0.0, 2.0 * lambda0));
    double radius1 = sqrt(fmax(0.0, 2.0 * lambda1));
    if (radius0 < 0.75 * base_radius) radius0 = 0.75 * base_radius;
    if (radius0 > 1.75 * base_radius) radius0 = 1.75 * base_radius;
    if (radius1 < 0.55 * base_radius) radius1 = 0.55 * base_radius;
    if (radius1 > radius0) radius1 = radius0;

    int best_sign = 1;
    double best_phase = 0.0;
    double best_error = DBL_MAX;
    for (int sign_index = 0; sign_index < 2; sign_index++) {
        int sign = sign_index == 0 ? 1 : -1;
        for (int sample = 0; sample < 256; sample++) {
            double phase = 2.0 * pi * (double)sample / 256.0;
            double error = 0.0;
            for (size_t i = 0; i < boundary_edges; i++) {
                double theta = phase + sign * 2.0 * pi * arc[i] / perimeter;
                double c = cos(theta), s = sin(theta);
                double q0 = center_uv[0] + radius0 * c * axis0[0] +
                            radius1 * s * axis1[0];
                double q1 = center_uv[1] + radius0 * c * axis0[1] +
                            radius1 * s * axis1[1];
                int32_t vi = loop[i];
                double e0 = q0 - projected_uv[(size_t)vi * 2];
                double e1 = q1 - projected_uv[(size_t)vi * 2 + 1];
                error += e0 * e0 + e1 * e1;
            }
            if (error < best_error) {
                best_error = error;
                best_phase = phase;
                best_sign = sign;
            }
        }
    }

    double *current_uv = (double *)ARENA_ALLOC(
        arena, nv * 2 * sizeof(*current_uv));
    double *next_uv = (double *)ARENA_ALLOC(
        arena, nv * 2 * sizeof(*next_uv));
    uint8_t *boundary_vertex = (uint8_t *)ARENA_CALLOC(
        arena, nv, sizeof(*boundary_vertex));
    memcpy(current_uv, projected_uv, nv * 2 * sizeof(*current_uv));
    for (size_t vi = 0; vi < nv; vi++) {
        if (!disk_vertex[vi]) continue;
        current_uv[vi * 2] = center_uv[0];
        current_uv[vi * 2 + 1] = center_uv[1];
    }
    for (size_t i = 0; i < boundary_edges; i++) {
        int32_t vi = loop[i];
        double theta = best_phase + best_sign * 2.0 * pi * arc[i] / perimeter;
        double c = cos(theta), s = sin(theta);
        current_uv[(size_t)vi * 2] =
            center_uv[0] + radius0 * c * axis0[0] + radius1 * s * axis1[0];
        current_uv[(size_t)vi * 2 + 1] =
            center_uv[1] + radius0 * c * axis0[1] + radius1 * s * axis1[1];
        boundary_vertex[vi] = 1;
    }

    CSR_T adjacency = CSR_from_faces(arena, disk_faces, disk_nf, nv);
    const int32_t *offset = CSR_offset(adjacency);
    const int32_t *target = CSR_target(adjacency);
    int iterations = 0;
    for (; iterations < 4000; iterations++) {
        memcpy(next_uv, current_uv, nv * 2 * sizeof(*next_uv));
        double max_delta = 0.0;
        for (size_t vi = 0; vi < nv; vi++) {
            if (!disk_vertex[vi] || boundary_vertex[vi]) continue;
            int32_t begin = offset[vi], end = offset[vi + 1];
            if (end <= begin) return -1;
            double q0 = 0.0, q1 = 0.0;
            for (int32_t pos = begin; pos < end; pos++) {
                int32_t vj = target[pos];
                q0 += current_uv[(size_t)vj * 2];
                q1 += current_uv[(size_t)vj * 2 + 1];
            }
            q0 /= (double)(end - begin);
            q1 /= (double)(end - begin);
            double d0 = q0 - current_uv[vi * 2];
            double d1 = q1 - current_uv[vi * 2 + 1];
            double delta = sqrt(d0 * d0 + d1 * d1);
            if (delta > max_delta) max_delta = delta;
            next_uv[vi * 2] = q0;
            next_uv[vi * 2 + 1] = q1;
        }
        double *swap = current_uv;
        current_uv = next_uv;
        next_uv = swap;
        if (max_delta < 1.0e-7 * fmax(1.0, base_radius)) break;
    }
    memcpy(chart_uv, current_uv, nv * 2 * sizeof(*chart_uv));

    size_t positive = 0, negative = 0, degenerate = 0;
    for (size_t fi = 0; fi < disk_nf; fi++) {
        int32_t a = disk_faces[fi * 3];
        int32_t b = disk_faces[fi * 3 + 1];
        int32_t c = disk_faces[fi * 3 + 2];
        double area2 = orient2d(&chart_uv[(size_t)a * 2],
                                &chart_uv[(size_t)b * 2],
                                &chart_uv[(size_t)c * 2]);
        if (area2 > 1.0e-10) positive++;
        else if (area2 < -1.0e-10) negative++;
        else degenerate++;
    }
    fprintf(stderr,
            "      repair disk chart: boundary=%zu iterations=%d "
            "faces(+/−/0)=%zu/%zu/%zu radii=%.2f/%.2f\n",
            boundary_edges, iterations, positive, negative, degenerate,
            radius0, radius1);
    if (degenerate != 0 || (positive != 0 && negative != 0)) return -1;
    return 0;
}

static int build_candidate(Arena_T arena,
                           const ComponentMesh *mesh,
                           const uint8_t *selected,
                           size_t selected_count,
                           const int32_t *face_labels,
                           int32_t label_a,
                           int32_t label_b,
                           const float common_normal[3],
                           const float common_center[3],
                           double common_plane_rms,
                           const OverlapMergeConfig *config,
                           const OMShape *input_shape,
                           uint8_t *out_collar_absorb,
                           ComponentMesh *out,
                           OverlapMergeStats *stats,
                           const char **out_reason)
{
    const size_t nv = mesh->nv;
    const size_t nf = mesh->nf;
    if (out_collar_absorb)
        memset(out_collar_absorb, 0, nf * sizeof(*out_collar_absorb));
    /* Labels identify the two local plies, but a merge candidate is not a
     * multicut candidate: preserve the original, single exterior mesh.  Only
     * wholly interior selected vertices retain one sample per local ply. */
    uint8_t *selected_label_mask = (uint8_t *)ARENA_CALLOC(
        arena, nv, sizeof(*selected_label_mask));
    uint8_t *kept_vertex = (uint8_t *)ARENA_CALLOC(
        arena, nv, sizeof(*kept_vertex));
    size_t kept_nf = 0;

    for (size_t fi = 0; fi < nf; fi++) {
        const int32_t *tri = &mesh->faces[fi * 3];
        uint8_t bit = (uint8_t)(face_labels[fi] == label_a ? 1u :
                                face_labels[fi] == label_b ? 2u : 0u);
        assert(bit != 0);
        if (selected[fi]) {
            for (int k = 0; k < 3; k++)
                selected_label_mask[tri[k]] |= bit;
        } else {
            kept_nf++;
            for (int k = 0; k < 3; k++) kept_vertex[tri[k]] = 1;
        }
    }
    if (kept_nf == 0 || selected_count == 0) {
        *out_reason = "empty-exterior";
        return -1;
    }
    {
        int32_t *kept_faces_diag = (int32_t *)ARENA_ALLOC(
            arena, kept_nf * 3 * sizeof(*kept_faces_diag));
        size_t kept_cursor_diag = 0;
        for (size_t fi = 0; fi < nf; fi++) {
            if (selected[fi]) continue;
            memcpy(&kept_faces_diag[kept_cursor_diag * 3],
                   &mesh->faces[fi * 3], 3 * sizeof(*kept_faces_diag));
            kept_cursor_diag++;
        }
        OMShape kept_shape_diag = mesh_shape(
            arena, kept_faces_diag, kept_nf, nv);
        fprintf(stderr,
                "      retained complement topology: components=%zu "
                "loops=%zu simple=%s\n",
                kept_shape_diag.face_components,
                kept_shape_diag.boundary_loops,
                kept_shape_diag.simple_boundary ? "yes" : "no");
    }

    size_t patch_nv = 0;
    for (size_t vi = 0; vi < nv; vi++) {
        if (selected_label_mask[vi]) patch_nv++;
    }
    if (patch_nv < 3) {
        *out_reason = "too-few-patch-vertices";
        return -1;
    }

    float normal[3] = {
        common_normal[0], common_normal[1], common_normal[2]
    };
    float center[3] = {
        common_center[0], common_center[1], common_center[2]
    };
    float axis_u[3], axis_v[3];
    PCA_orthonormal_basis(normal, axis_u, axis_v);

    double *uv = (double *)ARENA_ALLOC(arena, nv * 2 * sizeof(*uv));
    for (size_t vi = 0; vi < nv; vi++) {
        double dz = (double)mesh->verts[vi * 3] - center[0];
        double dy = (double)mesh->verts[vi * 3 + 1] - center[1];
        double dx = (double)mesh->verts[vi * 3 + 2] - center[2];
        uv[vi * 2] = dz * axis_u[0] + dy * axis_u[1] + dx * axis_u[2];
        uv[vi * 2 + 1] =
            dz * axis_v[0] + dy * axis_v[1] + dx * axis_v[2];
    }
    double *chart_uv = (double *)ARENA_ALLOC(
        arena, nv * 2 * sizeof(*chart_uv));
    /* The two source neighborhoods are intentionally separate disks: their
     * overlap in this common tangent plane is the delamination hypothesis we
     * are testing.  Their planar UNION, not a source-topology connector, is
     * the healed domain. */
    memcpy(chart_uv, uv, nv * 2 * sizeof(*chart_uv));

    if (getenv("VES_OVERLAP_DEBUG_DIR")) {
        uint8_t *piece_a = (uint8_t *)ARENA_CALLOC(
            arena, nf, sizeof(*piece_a));
        uint8_t *piece_b = (uint8_t *)ARENA_CALLOC(
            arena, nf, sizeof(*piece_b));
        float *plane_verts = (float *)ARENA_ALLOC(
            arena, nv * 3 * sizeof(*plane_verts));
        float *chart_verts = (float *)ARENA_ALLOC(
            arena, nv * 3 * sizeof(*chart_verts));
        for (size_t fi = 0; fi < nf; fi++) {
            if (!selected[fi]) continue;
            if (face_labels[fi] == label_a) piece_a[fi] = 1;
            else if (face_labels[fi] == label_b) piece_b[fi] = 1;
        }
        for (size_t vi = 0; vi < nv; vi++) {
            double u = uv[vi * 2], v = uv[vi * 2 + 1];
            for (int d = 0; d < 3; d++)
                plane_verts[vi * 3 + (size_t)d] = center[d] +
                    (float)(u * axis_u[d] + v * axis_v[d]);
            u = chart_uv[vi * 2];
            v = chart_uv[vi * 2 + 1];
            for (int d = 0; d < 3; d++)
                chart_verts[vi * 3 + (size_t)d] = center[d] +
                    (float)(u * axis_u[d] + v * axis_v[d]);
        }
        debug_write_mask_obj(arena, "overlap_merge_ring5_piece_0.obj",
                             mesh->verts, nv, mesh->faces, nf, piece_a);
        debug_write_mask_obj(arena, "overlap_merge_ring5_piece_1.obj",
                             mesh->verts, nv, mesh->faces, nf, piece_b);
        debug_write_mask_obj(
            arena, "overlap_merge_ring5_piece_0_projected.obj",
            plane_verts, nv, mesh->faces, nf, piece_a);
        debug_write_mask_obj(
            arena, "overlap_merge_ring5_piece_1_projected.obj",
            plane_verts, nv, mesh->faces, nf, piece_b);
        debug_write_mask_obj(
            arena, "overlap_merge_repair_disk_chart.obj",
            chart_verts, nv, mesh->faces, nf, selected);
    }

    size_t n_pinned = 0;
    size_t n_interior = 0;
    for (size_t vi = 0; vi < nv; vi++) {
        uint8_t selected_mask = selected_label_mask[vi];
        if (!selected_mask) continue;
        if (kept_vertex[vi]) {
            n_pinned++;
        } else {
            if (selected_mask & 1u) n_interior++;
            if (selected_mask & 2u) n_interior++;
        }
    }
    if (n_pinned < 3) {
        *out_reason = "too-few-frontier-vertices";
        return -1;
    }

    int32_t *kept_map = (int32_t *)ARENA_ALLOC(
        arena, nv * sizeof(*kept_map));
    for (size_t vi = 0; vi < nv; vi++) {
        kept_map[vi] = -1;
    }
    size_t n_kept_copies = 0;
    for (size_t vi = 0; vi < nv; vi++) {
        if (kept_vertex[vi]) kept_map[vi] = (int32_t)n_kept_copies++;
    }

    int32_t *pinned_vertex = (int32_t *)ARENA_ALLOC(
        arena, n_pinned * sizeof(*pinned_vertex));
    int32_t *pinned_global = (int32_t *)ARENA_ALLOC(
        arena, n_pinned * sizeof(*pinned_global));
    int32_t *interior_vertex = (int32_t *)ARENA_ALLOC(
        arena, n_interior * sizeof(*interior_vertex));
    uint8_t *interior_label = (uint8_t *)ARENA_ALLOC(
        arena, n_interior * sizeof(*interior_label));
    size_t pi = 0, ii = 0;
    for (size_t vi = 0; vi < nv; vi++) {
        uint8_t selected_mask = selected_label_mask[vi];
        if (!selected_mask) continue;
        if (kept_vertex[vi]) {
            pinned_vertex[pi] = (int32_t)vi;
            pinned_global[pi++] = kept_map[vi];
        } else {
            if (selected_mask & 1u) {
                interior_vertex[ii] = (int32_t)vi;
                interior_label[ii++] = 1u;
            }
            if (selected_mask & 2u) {
                interior_vertex[ii] = (int32_t)vi;
                interior_label[ii++] = 2u;
            }
        }
    }
    assert(pi == n_pinned && ii == n_interior);

    /*
     * A still-double frontier means the selected patch has not grown beyond
     * the overlap yet.  Do not weld fixed exterior samples; let the adaptive
     * ring loop grow until each frontier is single-cover.
     */
    /* Only true numerical duplicates are unusable fixed Delaunay points.
     * Near neighbours are legitimate frontier samples; whether the frontier
     * is still double-covered is decided by the zero-overlap audit after the
     * candidate is built, not by an arbitrary proximity surrogate. */
    double frontier_eps2 = 1.0e-12;
    for (size_t a = 0; a < n_pinned; a++) {
        int32_t va = pinned_vertex[a];
        for (size_t b = a + 1; b < n_pinned; b++) {
            int32_t vb = pinned_vertex[b];
            double du = chart_uv[(size_t)va * 2] -
                        chart_uv[(size_t)vb * 2];
            double dv = chart_uv[(size_t)va * 2 + 1] -
                        chart_uv[(size_t)vb * 2 + 1];
            if (du * du + dv * dv <= frontier_eps2) {
                *out_reason = "double-covered-frontier";
                return -1;
            }
        }
    }

    int32_t *cluster_parent = (int32_t *)ARENA_ALLOC(
        arena, n_interior * sizeof(*cluster_parent));
    uint8_t *cluster_rank = (uint8_t *)ARENA_CALLOC(
        arena, n_interior, sizeof(*cluster_rank));
    for (size_t i = 0; i < n_interior; i++) cluster_parent[i] = (int32_t)i;

    /* Preserve both source samplings.  Only mutually nearest samples from
     * opposite pieces may collapse; a transitive radius cluster could merge
     * several vertices from one ply through a single vertex on the other. */
    int32_t *nearest = (int32_t *)ARENA_ALLOC(
        arena, n_interior * sizeof(*nearest));
    double *nearest_d2 = (double *)ARENA_ALLOC(
        arena, n_interior * sizeof(*nearest_d2));
    for (size_t i = 0; i < n_interior; i++) {
        nearest[i] = -1;
        nearest_d2[i] = DBL_MAX;
    }
    double merge_d2 = config->merge_distance * config->merge_distance;
    for (size_t a = 0; a < n_interior; a++) {
        int32_t va = interior_vertex[a];
        for (size_t b = a + 1; b < n_interior; b++) {
            if (interior_label[a] == interior_label[b]) continue;
            int32_t vb = interior_vertex[b];
            double du = uv[(size_t)va * 2] - uv[(size_t)vb * 2];
            double dv = uv[(size_t)va * 2 + 1] - uv[(size_t)vb * 2 + 1];
            double d2 = du * du + dv * dv;
            if (d2 < nearest_d2[a]) {
                nearest_d2[a] = d2;
                nearest[a] = (int32_t)b;
            }
            if (d2 < nearest_d2[b]) {
                nearest_d2[b] = d2;
                nearest[b] = (int32_t)a;
            }
        }
    }
    for (size_t a = 0; a < n_interior; a++) {
        int32_t b = nearest[a];
        if (b >= 0 && (size_t)b > a && nearest[b] == (int32_t)a &&
            nearest_d2[a] <= merge_d2) {
            om_union(cluster_parent, cluster_rank, (int32_t)a, b);
        }
    }

    int32_t *root_group = (int32_t *)ARENA_ALLOC(
        arena, n_interior * sizeof(*root_group));
    for (size_t i = 0; i < n_interior; i++) root_group[i] = -1;
    size_t n_groups = 0;
    for (size_t i = 0; i < n_interior; i++) {
        int32_t root = om_find(cluster_parent, (int32_t)i);
        if (root_group[root] < 0) root_group[root] = (int32_t)n_groups++;
    }

    size_t n_points = n_pinned + n_groups;
    if (n_points < 3 || n_points > (size_t)INT32_MAX) {
        *out_reason = "invalid-delaunay-point-count";
        return -1;
    }
    double *points = (double *)ARENA_CALLOC(
        arena, n_points * 2, sizeof(*points));
    float *point_targets = (float *)ARENA_CALLOC(
        arena, n_points * 3, sizeof(*point_targets));
    int32_t *point_global = (int32_t *)ARENA_ALLOC(
        arena, n_points * sizeof(*point_global));
    int32_t *point_for_a = (int32_t *)ARENA_ALLOC(
        arena, nv * sizeof(*point_for_a));
    int32_t *point_for_b = (int32_t *)ARENA_ALLOC(
        arena, nv * sizeof(*point_for_b));
    for (size_t vi = 0; vi < nv; vi++) {
        point_for_a[vi] = -1;
        point_for_b[vi] = -1;
    }
    for (size_t i = 0; i < n_pinned; i++) {
        int32_t vi = pinned_vertex[i];
        points[i * 2] = chart_uv[(size_t)vi * 2];
        points[i * 2 + 1] = chart_uv[(size_t)vi * 2 + 1];
        memcpy(&point_targets[i * 3], &mesh->verts[(size_t)vi * 3],
               3 * sizeof(float));
        point_global[i] = pinned_global[i];
        if (selected_label_mask[vi] & 1u) point_for_a[vi] = (int32_t)i;
        if (selected_label_mask[vi] & 2u) point_for_b[vi] = (int32_t)i;
    }

    int32_t *group_count = (int32_t *)ARENA_CALLOC(
        arena, n_groups, sizeof(*group_count));
    for (size_t i = 0; i < n_interior; i++) {
        int32_t root = om_find(cluster_parent, (int32_t)i);
        int32_t group = root_group[root];
        int32_t vi = interior_vertex[i];
        size_t p = n_pinned + (size_t)group;
        points[p * 2] += chart_uv[(size_t)vi * 2];
        points[p * 2 + 1] += chart_uv[(size_t)vi * 2 + 1];
        for (int d = 0; d < 3; d++)
            point_targets[p * 3 + (size_t)d] +=
                mesh->verts[(size_t)vi * 3 + (size_t)d];
        group_count[group]++;
        if (interior_label[i] == 1u) point_for_a[vi] = (int32_t)p;
        else point_for_b[vi] = (int32_t)p;
    }
    for (size_t g = 0; g < n_groups; g++) {
        size_t p = n_pinned + g;
        points[p * 2] /= (double)group_count[g];
        points[p * 2 + 1] /= (double)group_count[g];
        for (int d = 0; d < 3; d++)
            point_targets[p * 3 + (size_t)d] /= (float)group_count[g];
        point_global[p] = (int32_t)(n_kept_copies + g);
    }

    /* Constrain the boundary of the selected patch in the original topology,
     * not the boundaries of two multicut pieces.  Keeping those two artificial
     * contours is equivalent to doing the split first and creates a three-face
     * seam wherever their planar domains overlap. */
    OMSourceEdge *source_edges = (OMSourceEdge *)ARENA_ALLOC(
        arena, nf * 3 * sizeof(*source_edges));
    for (size_t fi = 0; fi < nf; fi++) {
        const int32_t *tri = &mesh->faces[fi * 3];
        for (int e = 0; e < 3; e++) {
            OMSourceEdge *rec = &source_edges[fi * 3 + (size_t)e];
            rec->a = tri[e];
            rec->b = tri[(e + 1) % 3];
            rec->key = edge_key(rec->a, rec->b);
            rec->face = (int32_t)fi;
            rec->label = face_labels[fi];
        }
    }
    qsort(source_edges, nf * 3, sizeof(*source_edges), compare_source_edges);

    uint64_t *constraint_keys = (uint64_t *)ARENA_ALLOC(
        arena, selected_count * 3 * sizeof(*constraint_keys));
    size_t n_constraint_keys = 0;
    size_t edge_cursor = 0;
    while (edge_cursor < nf * 3) {
        size_t edge_end = edge_cursor + 1;
        while (edge_end < nf * 3 &&
               source_edges[edge_end].key == source_edges[edge_cursor].key)
            edge_end++;
        int source_boundary = edge_end - edge_cursor == 1;
        size_t selected_incident = 0;
        size_t kept_incident = 0;
        const OMSourceEdge *sample = NULL;
        for (size_t r = edge_cursor; r < edge_end; r++) {
            if (selected[source_edges[r].face]) {
                selected_incident++;
                sample = &source_edges[r];
            } else {
                kept_incident++;
            }
        }
        if (selected_incident > 0 &&
            (kept_incident > 0 || source_boundary)) {
            int32_t pa = -1, pb = -1;
            if (kept_vertex[sample->a]) {
                pa = point_for_a[sample->a] >= 0 ?
                    point_for_a[sample->a] : point_for_b[sample->a];
            } else {
                int32_t *point_for = sample->label == label_a ?
                    point_for_a : point_for_b;
                pa = point_for[sample->a];
            }
            if (kept_vertex[sample->b]) {
                pb = point_for_a[sample->b] >= 0 ?
                    point_for_a[sample->b] : point_for_b[sample->b];
            } else {
                int32_t *point_for = sample->label == label_a ?
                    point_for_a : point_for_b;
                pb = point_for[sample->b];
            }
            if (pa < 0 || pb < 0 || pa == pb) {
                *out_reason = "collapsed-frontier-segment";
                return -1;
            }
            constraint_keys[n_constraint_keys++] = edge_key(pa, pb);
        }
        edge_cursor = edge_end;
    }
    if (n_constraint_keys < 3) {
        *out_reason = "too-few-frontier-segments";
        return -1;
    }
    qsort(constraint_keys, n_constraint_keys, sizeof(*constraint_keys),
          compare_uint64);
    size_t n_segments = 0;
    for (size_t i = 0; i < n_constraint_keys; i++) {
        if (i == 0 || constraint_keys[i] != constraint_keys[i - 1])
            constraint_keys[n_segments++] = constraint_keys[i];
    }
    int32_t *segments = (int32_t *)ARENA_ALLOC(
        arena, n_segments * 2 * sizeof(*segments));
    for (size_t i = 0; i < n_segments; i++) {
        segments[i * 2] = (int32_t)(constraint_keys[i] >> 32);
        segments[i * 2 + 1] =
            (int32_t)(constraint_keys[i] & UINT32_C(0xffffffff));
    }
    fprintf(stderr, "      constrained five-ring frontier: %zu segments\n",
            n_segments);

    OMDomain domain;
    if (domain_build(arena, chart_uv, mesh->faces, nf, selected, &domain) != 0) {
        *out_reason = "domain-index-failed";
        return -1;
    }

    double *tri_points = NULL;
    size_t n_tri_points = 0;
    int32_t *delaunay_faces = NULL;
    size_t n_delaunay = 0;
    if (delaunay_points(arena, points, n_points, segments, n_segments, 0.0,
                        &tri_points, &n_tri_points,
                        &delaunay_faces, &n_delaunay) != 0) {
        *out_reason = "delaunay-failed";
        return -1;
    }
    for (size_t i = 0; i < n_points; i++) {
        if (fabs(tri_points[i * 2] - points[i * 2]) > 1.0e-9 ||
            fabs(tri_points[i * 2 + 1] - points[i * 2 + 1]) > 1.0e-9) {
            *out_reason = "delaunay-reordered-input";
            return -1;
        }
    }

    float *tri_targets = (float *)ARENA_ALLOC(
        arena, n_tri_points * 3 * sizeof(*tri_targets));
    memcpy(tri_targets, point_targets, n_points * 3 * sizeof(*tri_targets));
    for (size_t i = n_points; i < n_tri_points; i++) {
        const double *q = &tri_points[i * 2];
        double target_sum[3] = {0.0, 0.0, 0.0};
        size_t target_count = 0;
        for (size_t si = 0; si < n_segments; si++) {
            int32_t ia = segments[si * 2];
            int32_t ib = segments[si * 2 + 1];
            const double *a = &points[(size_t)ia * 2];
            const double *b = &points[(size_t)ib * 2];
            double dx = b[0] - a[0], dy = b[1] - a[1];
            double length2 = dx * dx + dy * dy;
            if (length2 <= 1.0e-18) continue;
            double t = ((q[0] - a[0]) * dx + (q[1] - a[1]) * dy) /
                       length2;
            if (t < -1.0e-9 || t > 1.0 + 1.0e-9) continue;
            double ex = q[0] - (a[0] + t * dx);
            double ey = q[1] - (a[1] + t * dy);
            if (ex * ex + ey * ey > 1.0e-14 * fmax(1.0, length2))
                continue;
            if (t < 0.0) t = 0.0;
            if (t > 1.0) t = 1.0;
            for (int d = 0; d < 3; d++)
                target_sum[d] +=
                    (1.0 - t) * point_targets[(size_t)ia * 3 + (size_t)d] +
                    t * point_targets[(size_t)ib * 3 + (size_t)d];
            target_count++;
        }
        if (target_count == 0) {
            size_t nearest_point = 0;
            double closest_point_d2 = DBL_MAX;
            for (size_t p = 0; p < n_points; p++) {
                double dx = q[0] - points[p * 2];
                double dy = q[1] - points[p * 2 + 1];
                double d2 = dx * dx + dy * dy;
                if (d2 < closest_point_d2) {
                    closest_point_d2 = d2;
                    nearest_point = p;
                }
            }
            memcpy(&tri_targets[i * 3], &point_targets[nearest_point * 3],
                   3 * sizeof(float));
        } else {
            for (int d = 0; d < 3; d++)
                tri_targets[i * 3 + (size_t)d] =
                    (float)(target_sum[d] / (double)target_count);
        }
    }

    /* Pass one is only the source-frontier arrangement.  Convert its
     * coverage to a true Boolean union, close all bounded gaps, and extract
     * exactly one exterior contour. */
    uint8_t *arrangement_keep = (uint8_t *)ARENA_CALLOC(
        arena, n_delaunay, sizeof(*arrangement_keep));
    for (size_t fi = 0; fi < n_delaunay; fi++) {
        int32_t ia = delaunay_faces[fi * 3];
        int32_t ib = delaunay_faces[fi * 3 + 1];
        int32_t ic = delaunay_faces[fi * 3 + 2];
        const double *a = &tri_points[(size_t)ia * 2];
        const double *b = &tri_points[(size_t)ib * 2];
        const double *c = &tri_points[(size_t)ic * 2];
        if (fabs(orient2d(a, b, c)) <= 1.0e-12) continue;
        double centroid[2] = {
            (a[0] + b[0] + c[0]) / 3.0,
            (a[1] + b[1] + c[1]) / 3.0
        };
        if (domain_coverage(&domain, centroid) > 0)
            arrangement_keep[fi] = 1;
    }
    size_t arrangement_gap_components = 0;
    size_t arrangement_gap_faces = fill_bounded_arrangement_gaps(
        arena, delaunay_faces, n_delaunay, arrangement_keep,
        &arrangement_gap_components);
    int32_t *outer_segments = NULL;
    size_t n_outer_segments = 0;
    if (extract_kept_boundary(
            arena, delaunay_faces, n_delaunay, arrangement_keep,
            n_tri_points, &outer_segments, &n_outer_segments) != 0) {
        *out_reason = "boolean-union-not-one-disk";
        return -1;
    }
    fprintf(stderr,
            "      Boolean union: exterior=%zu segments, filled=%zu faces "
            "in %zu bounded gap(s)\n",
            n_outer_segments, arrangement_gap_faces,
            arrangement_gap_components);

    /* Pass two deliberately forgets the two source contours.  The only PSLG
     * constraint is the healed disk's exterior contour; all source samples
     * become ordinary interior observations of the same sheet. */
    double *healed_points = NULL;
    size_t n_healed_points = 0;
    int32_t *healed_faces = NULL;
    size_t n_healed_faces = 0;
    if (delaunay_points(arena, tri_points, n_tri_points,
                        outer_segments, n_outer_segments, 0.35,
                        &healed_points, &n_healed_points,
                        &healed_faces, &n_healed_faces) != 0) {
        *out_reason = "healed-delaunay-failed";
        return -1;
    }
    for (size_t i = 0; i < n_tri_points; i++) {
        if (fabs(healed_points[i * 2] - tri_points[i * 2]) > 1.0e-9 ||
            fabs(healed_points[i * 2 + 1] - tri_points[i * 2 + 1]) > 1.0e-9) {
            *out_reason = "healed-delaunay-reordered-input";
            return -1;
        }
    }
    float *healed_targets = (float *)ARENA_ALLOC(
        arena, n_healed_points * 3 * sizeof(*healed_targets));
    memcpy(healed_targets, tri_targets,
           n_tri_points * 3 * sizeof(*healed_targets));
    for (size_t i = n_tri_points; i < n_healed_points; i++) {
        size_t nearest_point = 0;
        double best_d2 = DBL_MAX;
        for (size_t p = 0; p < n_tri_points; p++) {
            double du = healed_points[i * 2] - tri_points[p * 2];
            double dv = healed_points[i * 2 + 1] - tri_points[p * 2 + 1];
            double d2 = du * du + dv * dv;
            if (d2 < best_d2) {
                best_d2 = d2;
                nearest_point = p;
            }
        }
        memcpy(&healed_targets[i * 3], &tri_targets[nearest_point * 3],
               3 * sizeof(*healed_targets));
    }
    tri_points = healed_points;
    n_tri_points = n_healed_points;
    delaunay_faces = healed_faces;
    n_delaunay = n_healed_faces;
    tri_targets = healed_targets;

    float *tri_planar = (float *)ARENA_ALLOC(
        arena, n_tri_points * 3 * sizeof(*tri_planar));
    for (size_t p = 0; p < n_tri_points; p++) {
        double u = tri_points[p * 2], v = tri_points[p * 2 + 1];
        for (int d = 0; d < 3; d++)
            tri_planar[p * 3 + (size_t)d] = center[d] +
                (float)(u * axis_u[d] + v * axis_v[d]);
    }

    size_t n_crossing_points = n_tri_points - n_points;
    int32_t *tri_point_global = (int32_t *)ARENA_ALLOC(
        arena, n_tri_points * sizeof(*tri_point_global));
    memcpy(tri_point_global, point_global,
           n_points * sizeof(*tri_point_global));
    for (size_t i = n_points; i < n_tri_points; i++)
        tri_point_global[i] =
            (int32_t)(n_kept_copies + n_groups + (i - n_points));

    size_t raw_nv = n_kept_copies + n_groups + n_crossing_points;
    float *raw_verts = (float *)ARENA_ALLOC(
        arena, raw_nv * 3 * sizeof(*raw_verts));
    for (size_t vi = 0; vi < nv; vi++) {
        if (kept_map[vi] >= 0) {
            memcpy(&raw_verts[(size_t)kept_map[vi] * 3], &mesh->verts[vi * 3],
                   3 * sizeof(float));
        }
    }
    for (size_t p = n_pinned; p < n_tri_points; p++) {
        int32_t global = tri_point_global[p];
        double u = tri_points[p * 2];
        double v = tri_points[p * 2 + 1];
        raw_verts[(size_t)global * 3] =
            center[0] + (float)(u * axis_u[0] + v * axis_v[0]);
        raw_verts[(size_t)global * 3 + 1] =
            center[1] + (float)(u * axis_u[1] + v * axis_v[1]);
        raw_verts[(size_t)global * 3 + 2] =
            center[2] + (float)(u * axis_u[2] + v * axis_v[2]);
    }

    int32_t *raw_faces = (int32_t *)ARENA_ALLOC(
        arena, (kept_nf + n_delaunay) * 3 * sizeof(*raw_faces));
    int32_t *kept_source_faces = (int32_t *)ARENA_ALLOC(
        arena, kept_nf * sizeof(*kept_source_faces));
    size_t raw_nf = 0;
    for (size_t fi = 0; fi < nf; fi++) {
        if (selected[fi]) continue;
        const int32_t *src = &mesh->faces[fi * 3];
        raw_faces[raw_nf * 3] = kept_map[src[0]];
        raw_faces[raw_nf * 3 + 1] = kept_map[src[1]];
        raw_faces[raw_nf * 3 + 2] = kept_map[src[2]];
        kept_source_faces[raw_nf] = (int32_t)fi;
        assert(raw_faces[raw_nf * 3] >= 0 &&
               raw_faces[raw_nf * 3 + 1] >= 0 &&
               raw_faces[raw_nf * 3 + 2] >= 0);
        raw_nf++;
    }
    assert(raw_nf == kept_nf);

    uint8_t *union_keep = (uint8_t *)ARENA_CALLOC(
        arena, n_delaunay, sizeof(*union_keep));
    for (size_t fi = 0; fi < n_delaunay; fi++) {
        int32_t ia = delaunay_faces[fi * 3];
        int32_t ib = delaunay_faces[fi * 3 + 1];
        int32_t ic = delaunay_faces[fi * 3 + 2];
        const double *a = &tri_points[(size_t)ia * 2];
        const double *b = &tri_points[(size_t)ib * 2];
        const double *c = &tri_points[(size_t)ic * 2];
        if (fabs(orient2d(a, b, c)) <= 1.0e-12) continue;
        union_keep[fi] = 1;
    }
    size_t healed_gap_components = 0;
    size_t healed_gap_faces = fill_bounded_arrangement_gaps(
        arena, delaunay_faces, n_delaunay, union_keep,
        &healed_gap_components);
    int32_t *healed_boundary = NULL;
    size_t n_healed_boundary = 0;
    if (extract_kept_boundary(
            arena, delaunay_faces, n_delaunay, union_keep,
            n_tri_points, &healed_boundary, &n_healed_boundary) != 0) {
        *out_reason = "healed-domain-not-one-disk";
        return -1;
    }
    fprintf(stderr,
            "      healed CDT domain: boundary=%zu filled=%zu faces/%zu gaps\n",
            n_healed_boundary, healed_gap_faces, healed_gap_components);

    int32_t *union_faces = (int32_t *)ARENA_ALLOC(
        arena, n_delaunay * 3 * sizeof(*union_faces));
    int32_t *evenodd_faces = (int32_t *)ARENA_ALLOC(
        arena, n_delaunay * 3 * sizeof(*evenodd_faces));
    size_t n_union_faces = 0;
    size_t n_evenodd_faces = 0;
    for (size_t fi = 0; fi < n_delaunay; fi++) {
        int32_t ia = delaunay_faces[fi * 3];
        int32_t ib = delaunay_faces[fi * 3 + 1];
        int32_t ic = delaunay_faces[fi * 3 + 2];
        const double *a = &tri_points[(size_t)ia * 2];
        const double *b = &tri_points[(size_t)ib * 2];
        const double *c = &tri_points[(size_t)ic * 2];
        double area2 = orient2d(a, b, c);
        if (fabs(area2) <= 1.0e-12) continue;

        double centroid[2] = {
            (a[0] + b[0] + c[0]) / 3.0,
            (a[1] + b[1] + c[1]) / 3.0
        };
        int coverage = domain_coverage(&domain, centroid);
        if (!union_keep[fi]) continue;

        int32_t local[3] = {ia, ib, ic};
        if (area2 < 0.0) {
            local[1] = ic;
            local[2] = ib;
        }
        memcpy(&union_faces[n_union_faces * 3], local,
               3 * sizeof(*union_faces));
        n_union_faces++;
        if (coverage & 1) {
            memcpy(&evenodd_faces[n_evenodd_faces * 3], local,
                   3 * sizeof(*evenodd_faces));
            n_evenodd_faces++;
        }

        int32_t ga = tri_point_global[ia];
        int32_t gb = tri_point_global[ib];
        int32_t gc = tri_point_global[ic];
        if (ga == gb || gb == gc || gc == ga) continue;
        if (area2 > 0.0) {
            raw_faces[raw_nf * 3] = ga;
            raw_faces[raw_nf * 3 + 1] = gb;
            raw_faces[raw_nf * 3 + 2] = gc;
        } else {
            raw_faces[raw_nf * 3] = ga;
            raw_faces[raw_nf * 3 + 1] = gc;
            raw_faces[raw_nf * 3 + 2] = gb;
        }
        raw_nf++;
    }
    if (getenv("VES_OVERLAP_DEBUG_DIR")) {
        const char *dir = getenv("VES_OVERLAP_DEBUG_DIR");
        char path[1024];
        snprintf(path, sizeof(path),
                 "%s/overlap_merge_planar_union_patch.obj", dir);
        ObjIO_write(path, tri_planar, n_tri_points,
                    union_faces, n_union_faces);
        snprintf(path, sizeof(path),
                 "%s/overlap_merge_planar_evenodd_patch.obj", dir);
        ObjIO_write(path, tri_planar, n_tri_points,
                    evenodd_faces, n_evenodd_faces);

        float *union_lifted = laplacian_heightfield_lift_patch(
            arena, tri_planar, tri_targets, normal, n_tri_points,
            union_faces, n_union_faces, 600, 0.005f);
        snprintf(path, sizeof(path),
                 "%s/overlap_merge_heightfield_union_patch.obj", dir);
        ObjIO_write(path, union_lifted, n_tri_points,
                    union_faces, n_union_faces);
        float *evenodd_lifted = laplacian_heightfield_lift_patch(
            arena, tri_planar, tri_targets, normal, n_tri_points,
            evenodd_faces, n_evenodd_faces, 600, 0.005f);
        snprintf(path, sizeof(path),
                 "%s/overlap_merge_heightfield_evenodd_patch.obj", dir);
        ObjIO_write(path, evenodd_lifted, n_tri_points,
                    evenodd_faces, n_evenodd_faces);
        fprintf(stderr,
                "      debug: planar/heightfield union=%zu faces "
                "evenodd=%zu faces\n",
                n_union_faces, n_evenodd_faces);
    }
    if (raw_nf == kept_nf) {
        *out_reason = "empty-retriangulation";
        return -1;
    }

    /* Every repair frontier edge must have exactly one core triangle and one
     * retained collar triangle.  Never delete exterior/collar faces to make a
     * projected union fit: that turns a stitching error into real holes. */
    OMEdge *union_edges = (OMEdge *)ARENA_ALLOC(
        arena, raw_nf * 3 * sizeof(*union_edges));
    fill_edge_records(raw_faces, raw_nf, kept_nf, union_edges);
    size_t union_edge_cursor = 0;
    size_t conflict_edges = 0;
    size_t conflict_faces = 0;
    while (union_edge_cursor < raw_nf * 3) {
        size_t union_edge_end = union_edge_cursor + 1;
        while (union_edge_end < raw_nf * 3 &&
               union_edges[union_edge_end].key ==
                   union_edges[union_edge_cursor].key)
            union_edge_end++;
        if (union_edge_end - union_edge_cursor > 2) {
            size_t patch_incident = 0, kept_incident = 0;
            int32_t kept_face = -1;
            for (size_t r = union_edge_cursor; r < union_edge_end; r++) {
                if (union_edges[r].patch) patch_incident++;
                else {
                    kept_incident++;
                    kept_face = union_edges[r].face;
                }
            }
            if (patch_incident == 2 && kept_incident == 1 && kept_face >= 0) {
                conflict_edges++;
                int32_t source_face = kept_source_faces[kept_face];
                if (out_collar_absorb && !out_collar_absorb[source_face]) {
                    out_collar_absorb[source_face] = 1;
                    conflict_faces++;
                }
            } else {
                *out_reason = "unresolvable-union-edge";
                return -1;
            }
        }
        union_edge_cursor = union_edge_end;
    }
    if (conflict_edges > 0) {
        fprintf(stderr,
                "      collar constraint failure: %zu frontier edge(s) "
                "covered on both core sides; absorb %zu retained face(s)\n",
                conflict_edges, conflict_faces);
        *out_reason = "collar-frontier-double-covered";
        return -1;
    }

    MeshManifoldStats raw_manifold = MeshManifold_audit(
        arena, raw_nv, raw_faces, raw_nf);
    fprintf(stderr,
            "      raw constrained candidate: nv=%zu nf=%zu patch=%zu "
            "nm_edges=%zu nm_verts=%zu boundary=%zu same_dir=%zu\n",
            raw_nv, raw_nf, raw_nf - kept_nf,
            raw_manifold.nm_edges, raw_manifold.nm_verts,
            raw_manifold.boundary_edges, raw_manifold.same_dir_edges);

    if (raw_manifold.nm_edges != 0) {
        *out_reason = "non-manifold-union-edge";
        return -1;
    }
    if (raw_manifold.nm_verts > 0) {
        uint8_t *pinch = (uint8_t *)ARENA_CALLOC(
            arena, raw_nv, sizeof(*pinch));
        MeshManifold_mark_nonmanifold_vertices(
            arena, raw_nv, raw_faces, raw_nf, pinch);
        size_t absorb_faces = 0, mapped_pinches = 0;
        for (size_t raw_vi = 0; raw_vi < n_kept_copies; raw_vi++) {
            if (!pinch[raw_vi]) continue;
            int32_t source_vi = -1;
            for (size_t vi = 0; vi < nv; vi++) {
                if (kept_map[vi] == (int32_t)raw_vi) {
                    source_vi = (int32_t)vi;
                    break;
                }
            }
            if (source_vi < 0) continue;
            mapped_pinches++;
            for (size_t fi = 0; fi < nf; fi++) {
                if (selected[fi]) continue;
                const int32_t *tri = &mesh->faces[fi * 3];
                if (tri[0] != source_vi && tri[1] != source_vi &&
                    tri[2] != source_vi) continue;
                if (out_collar_absorb && !out_collar_absorb[fi]) {
                    out_collar_absorb[fi] = 1;
                    absorb_faces++;
                }
            }
        }
        fprintf(stderr,
                "      collar vertex constraint: pinches=%zu mapped=%zu "
                "absorb %zu retained one-ring face(s)\n",
                raw_manifold.nm_verts, mapped_pinches, absorb_faces);
        if (absorb_faces > 0) {
            *out_reason = "collar-vertex-pinch";
            return -1;
        }
        *out_reason = "unresolved-patch-vertex-pinch";
        return -1;
    }

    /* Always preserve the complete topologically manifold hypothesis before
     * the boundary-loop acceptance gate.  Rejected candidates are often the
     * most useful visual diagnostic and must not disappear just because they
     * were correctly prevented from entering the pipeline. */
    if (getenv("VES_OVERLAP_DEBUG_DIR")) {
        const char *dir = getenv("VES_OVERLAP_DEBUG_DIR");
        char path[1024];
        snprintf(path, sizeof(path),
                 "%s/overlap_merge_candidate_pre_fit.obj", dir);
        ObjIO_write(path, raw_verts, raw_nv, raw_faces, raw_nf);
        snprintf(path, sizeof(path),
                 "%s/overlap_merge_candidate_patch_pre_fit.obj", dir);
        ObjIO_write(path, raw_verts, raw_nv, &raw_faces[kept_nf * 3],
                    raw_nf - kept_nf);
    }

    size_t raw_boundary_loops = 0, loop_absorb_faces = 0;
    int loop_collar = mark_extra_boundary_loop_collar(
        arena, raw_faces, raw_nf, raw_nv, kept_nf, kept_source_faces,
        input_shape->boundary_loops, out_collar_absorb,
        &raw_boundary_loops, &loop_absorb_faces);
    fprintf(stderr,
            "      collar boundary constraint: loops=%zu expected=%zu "
            "absorb %zu retained face(s)\n",
            raw_boundary_loops, input_shape->boundary_loops,
            loop_absorb_faces);
    if (loop_collar > 0) {
        *out_reason = "collar-boundary-hole";
        return -1;
    }
    if (loop_collar < 0) {
        *out_reason = "unresolved-collar-boundary-hole";
        return -1;
    }

    if (getenv("VES_OVERLAP_DEBUG_DIR")) {
        const char *dir = getenv("VES_OVERLAP_DEBUG_DIR");
        char path[1024];
        snprintf(path, sizeof(path), "%s/overlap_merge_evenodd_raw.obj", dir);
        ObjIO_write(path, raw_verts, raw_nv, raw_faces, raw_nf);
        snprintf(path, sizeof(path), "%s/overlap_merge_evenodd_patch.obj", dir);
        ObjIO_write(path, raw_verts, raw_nv, &raw_faces[kept_nf * 3],
                    raw_nf - kept_nf);
    }

    size_t orient_flips = 0, orient_components = 0, orient_residual = 0;
    OrientMesh_consistent(arena, raw_verts, raw_nv, NULL,
                          raw_faces, raw_nf,
                          &orient_flips, &orient_components,
                          &orient_residual);
    fprintf(stderr,
            "      orientation: flips=%zu components=%zu residual=%zu\n",
            orient_flips, orient_components, orient_residual);
    if (orient_residual > 0) {
        OMEdge *diag_edges = (OMEdge *)ARENA_ALLOC(
            arena, raw_nf * 3 * sizeof(*diag_edges));
        fill_edge_records(raw_faces, raw_nf, kept_nf, diag_edges);
        size_t di = 0, shown = 0;
        while (di < raw_nf * 3) {
            size_t dj = di + 1;
            while (dj < raw_nf * 3 &&
                   diag_edges[dj].key == diag_edges[di].key) dj++;
            if (dj - di == 2 &&
                diag_edges[di].direction == diag_edges[di + 1].direction) {
                fprintf(stderr,
                        "        same-dir edge %u-%u faces=%d(%s),%d(%s)\n",
                        (unsigned)(diag_edges[di].key >> 32),
                        (unsigned)(diag_edges[di].key & UINT32_C(0xffffffff)),
                        diag_edges[di].face,
                        diag_edges[di].patch ? "patch" : "kept",
                        diag_edges[di + 1].face,
                        diag_edges[di + 1].patch ? "patch" : "kept");
                shown++;
            }
            di = dj;
        }
        assert(shown == orient_residual);
    }
    if (orient_residual != 0) {
        *out_reason = "candidate-not-orientable";
        return -1;
    }

    size_t same = 0, opposite = 0;
    if (interface_orientation(arena, raw_faces, raw_nf, kept_nf,
                              &same, &opposite) != 0 || same != 0) {
        *out_reason = "interface-orientation-failed";
        return -1;
    }

    uint8_t *used = (uint8_t *)ARENA_CALLOC(
        arena, raw_nv, sizeof(*used));
    for (size_t fi = 0; fi < raw_nf; fi++) {
        used[raw_faces[fi * 3]] = 1;
        used[raw_faces[fi * 3 + 1]] = 1;
        used[raw_faces[fi * 3 + 2]] = 1;
    }
    int32_t *remap = (int32_t *)ARENA_ALLOC(
        arena, raw_nv * sizeof(*remap));
    size_t compact_nv = 0;
    for (size_t vi = 0; vi < raw_nv; vi++) {
        if (used[vi]) remap[vi] = (int32_t)compact_nv++;
        else remap[vi] = -1;
    }

    float *compact_verts = (float *)ARENA_ALLOC(
        arena, compact_nv * 3 * sizeof(*compact_verts));
    uint8_t *compact_pin = (uint8_t *)ARENA_CALLOC(
        arena, compact_nv, sizeof(*compact_pin));
    for (size_t vi = 0; vi < raw_nv; vi++) {
        if (remap[vi] < 0) continue;
        size_t dst = (size_t)remap[vi];
        memcpy(&compact_verts[dst * 3], &raw_verts[vi * 3],
               3 * sizeof(float));
        if (vi < n_kept_copies) compact_pin[dst] = 1;
    }
    int32_t *compact_faces = (int32_t *)ARENA_ALLOC(
        arena, raw_nf * 3 * sizeof(*compact_faces));
    for (size_t fi = 0; fi < raw_nf; fi++) {
        compact_faces[fi * 3] = remap[raw_faces[fi * 3]];
        compact_faces[fi * 3 + 1] = remap[raw_faces[fi * 3 + 1]];
        compact_faces[fi * 3 + 2] = remap[raw_faces[fi * 3 + 2]];
    }
    uint8_t *compact_role = (uint8_t *)ARENA_CALLOC(
        arena, compact_nv, sizeof(*compact_role));
    for (size_t fi = kept_nf; fi < raw_nf; fi++) {
        for (int k = 0; k < 3; k++) {
            int32_t vi = compact_faces[fi * 3 + (size_t)k];
            compact_role[vi] = compact_pin[vi]
                ? OVERLAP_MERGE_VERTEX_ANCHOR
                : OVERLAP_MERGE_VERTEX_REPAIR;
        }
    }
    float *compact_chart_uv = (float *)ARENA_ALLOC(
        arena, compact_nv * 2 * sizeof(*compact_chart_uv));
    for (size_t vi = 0; vi < compact_nv; vi++) {
        double dz = (double)compact_verts[vi * 3] - center[0];
        double dy = (double)compact_verts[vi * 3 + 1] - center[1];
        double dx = (double)compact_verts[vi * 3 + 2] - center[2];
        compact_chart_uv[vi * 2] = (float)(
            dz * axis_u[0] + dy * axis_u[1] + dx * axis_u[2]);
        compact_chart_uv[vi * 2 + 1] = (float)(
            dz * axis_v[0] + dy * axis_v[1] + dx * axis_v[2]);
    }

    MeshManifoldStats manifold = MeshManifold_audit(
        arena, compact_nv, compact_faces, raw_nf);
    if (!MeshManifold_ok(&manifold)) {
        fprintf(stderr,
                "      merge candidate topology: nm_edges=%zu nm_verts=%zu "
                "boundary=%zu same_dir=%zu\n",
                manifold.nm_edges, manifold.nm_verts,
                manifold.boundary_edges, manifold.same_dir_edges);
        *out_reason = "non-manifold-candidate";
        return -1;
    }
    if (manifold.same_dir_edges != 0) {
        *out_reason = "inconsistent-candidate-winding";
        return -1;
    }

    OMShape shape = mesh_shape(arena, compact_faces, raw_nf, compact_nv);
    if (shape.face_components != 1) {
        *out_reason = "candidate-not-one-component";
        return -1;
    }
    if (!shape.simple_boundary) {
        *out_reason = "candidate-boundary-not-simple";
        return -1;
    }
    fprintf(stderr,
            "      candidate boundary loops: before=%zu after=%zu\n",
            input_shape->boundary_loops, shape.boundary_loops);
    if (shape.boundary_loops != input_shape->boundary_loops) {
        *out_reason = "candidate-boundary-loop-change";
        return -1;
    }

    size_t patch_added = raw_nf - kept_nf;
    double *qualities = (double *)ARENA_ALLOC(
        arena, patch_added * sizeof(*qualities));
    for (size_t fi = kept_nf; fi < raw_nf; fi++) {
        double area = 0.0;
        double quality = triangle_quality(
            compact_verts, &compact_faces[fi * 3], &area);
        if (area < OM_MIN_FACE_AREA) {
            *out_reason = "degenerate-patch-face";
            return -1;
        }
        qualities[fi - kept_nf] = quality;
    }
    qsort(qualities, patch_added, sizeof(*qualities), compare_double);
    size_t q01_index = patch_added / 100;
    if (qualities[q01_index] < OM_MIN_PATCH_QUALITY) {
        *out_reason = "patch-quality-floor";
        return -1;
    }

    memset(out, 0, sizeof(*out));
    out->verts = compact_verts;
    out->faces = compact_faces;
    out->pin_mask = NULL;
    out->vert_normals = NULL;
    out->nv = compact_nv;
    out->nf = raw_nf;
    out->comp_id = mesh->comp_id;
    out->nv_pre_fill = 0;
    PCA_normal(out->verts, out->nv, out->pca_normal, out->centroid);
    out->self = out;

    stats->patch_faces_removed = selected_count;
    stats->patch_faces_added = patch_added;
    stats->patch_face_start = kept_nf;
    stats->patch_vertices_coalesced = n_interior - n_groups;
    stats->interface_edges = opposite;
    stats->boundary_loops_before = input_shape->boundary_loops;
    stats->boundary_loops_after = shape.boundary_loops;
    stats->plane_rms = common_plane_rms;
    stats->plane_normal[0] = normal[0];
    stats->plane_normal[1] = normal[1];
    stats->plane_normal[2] = normal[2];
    stats->plane_center[0] = center[0];
    stats->plane_center[1] = center[1];
    stats->plane_center[2] = center[2];
    stats->exterior_source_faces = kept_source_faces;
    stats->vertex_role = compact_role;
    stats->chart_uv = compact_chart_uv;
    *out_reason = "accepted";
    return 0;
}

void OverlapMerge_default_config(OverlapMergeConfig *config)
{
    assert(config);
    config->min_rings = 5;
    config->max_rings = 5;
    config->local_search_rings = 12;
    config->max_overlap_pairs = 1024;
    config->max_seed_faces = 768;
    config->max_patch_faces = 4096;
    config->max_patch_fraction = 0.055;
    config->max_seed_diameter = 48.0;
    config->max_patch_diameter = 64.0;
    config->max_plane_rms = 2.5;
    config->contact_distance = 0.90;
    config->merge_distance = 0.10;
}

int OverlapMerge_try(Arena_T arena,
                     const ComponentMesh *mesh,
                     const int32_t *face_labels,
                     int32_t num_labels,
                     const int32_t *adj_fa,
                     const int32_t *adj_fb,
                     size_t n_adj,
                     const int32_t *ovl_fa,
                     const int32_t *ovl_fb,
                     size_t n_ovl,
                     const OverlapMergeConfig *config,
                     ComponentMesh *out,
                     OverlapMergeStats *stats)
{
    assert(arena && mesh && face_labels && config && out && stats);
    memset(out, 0, sizeof(*out));
    memset(stats, 0, sizeof(*stats));
    stats->reason = "not-attempted";

    if (!adj_fa || !adj_fb || !ovl_fa || !ovl_fb ||
        num_labels <= 0 || mesh->nf > (size_t)INT32_MAX ||
        mesh->nv > (size_t)INT32_MAX) {
        stats->reason = "invalid-input";
        return -1;
    }

    int32_t *label_count = (int32_t *)ARENA_CALLOC(
        arena, (size_t)num_labels, sizeof(*label_count));
    for (size_t fi = 0; fi < mesh->nf; fi++) {
        int32_t label = face_labels[fi];
        if (label < 0 || label >= num_labels) {
            stats->reason = "invalid-face-label";
            return -1;
        }
        label_count[label]++;
    }

    int32_t labels[2] = {-1, -1};
    int n_nonempty = 0;
    for (int32_t label = 0; label < num_labels; label++) {
        if (label_count[label] == 0) continue;
        if (n_nonempty < 2) labels[n_nonempty] = label;
        n_nonempty++;
    }
    if (n_nonempty != 2) {
        stats->reason = "requires-exactly-two-labels";
        return -1;
    }

    if (n_ovl == 0 || n_ovl > config->max_overlap_pairs) {
        stats->reason = "overlap-pair-budget";
        return -1;
    }
    stats->overlap_pairs = n_ovl;

    uint8_t *seed = (uint8_t *)ARENA_CALLOC(
        arena, mesh->nf, sizeof(*seed));
    for (size_t i = 0; i < n_ovl; i++) {
        int32_t fa = ovl_fa[i];
        int32_t fb = ovl_fb[i];
        if (fa < 0 || fb < 0 ||
            (size_t)fa >= mesh->nf || (size_t)fb >= mesh->nf) {
            stats->reason = "invalid-overlap-pair";
            return -1;
        }
        int32_t la = face_labels[fa];
        int32_t lb = face_labels[fb];
        if (la == lb) {
            stats->reason = "joined-repulsive-edge";
            return -1;
        }
        if (!((la == labels[0] && lb == labels[1]) ||
              (la == labels[1] && lb == labels[0]))) {
            stats->reason = "overlap-outside-label-pair";
            return -1;
        }
        seed[fa] = 1;
        seed[fb] = 1;
    }

    size_t seed_count = 0;
    for (size_t fi = 0; fi < mesh->nf; fi++) {
        if (seed[fi]) seed_count++;
    }
    stats->seed_faces = seed_count;
    if (seed_count == 0 || seed_count > config->max_seed_faces) {
        stats->reason = "seed-face-budget";
        return -1;
    }
    stats->seed_diameter = mask_diameter(mesh, seed);
    if (stats->seed_diameter > config->max_seed_diameter) {
        stats->reason = "seed-diameter-budget";
        return -1;
    }

    /*
     * Count connected PRODUCTION overlap events before replacing the detector's
     * seed by the physical-contact seed.  Distant contacts on the same two
     * stacked charts are not one small delamination, while tiny disconnected
     * near-contact specks inside this bounded event must not create false
     * extra events.
     */
    int32_t *event_parent = (int32_t *)ARENA_ALLOC(
        arena, mesh->nf * sizeof(*event_parent));
    uint8_t *event_rank = (uint8_t *)ARENA_CALLOC(
        arena, mesh->nf, sizeof(*event_rank));
    for (size_t fi = 0; fi < mesh->nf; fi++)
        event_parent[fi] = seed[fi] ? (int32_t)fi : -1;
    for (size_t i = 0; i < n_adj; i++) {
        int32_t fa = adj_fa[i], fb = adj_fb[i];
        if (seed[fa] && seed[fb])
            om_union(event_parent, event_rank, fa, fb);
    }
    for (size_t i = 0; i < n_ovl; i++)
        om_union(event_parent, event_rank, ovl_fa[i], ovl_fb[i]);

    size_t event_count = 0;
    for (size_t fi = 0; fi < mesh->nf; fi++) {
        if (seed[fi] &&
            om_find(event_parent, (int32_t)fi) == (int32_t)fi)
            event_count++;
    }
    if (event_count != 1) {
        stats->reason = "multiple-overlap-events";
        return -1;
    }

    OMShape input_shape = mesh_shape(
        arena, mesh->faces, mesh->nf, mesh->nv);
    if (input_shape.face_components != 1 || !input_shape.simple_boundary) {
        stats->reason = "input-topology-not-simple";
        return -1;
    }

    size_t fraction_cap = (size_t)floor(
        config->max_patch_fraction * (double)mesh->nf);
    size_t patch_cap = config->max_patch_faces;
    if (fraction_cap < patch_cap) patch_cap = fraction_cap;
    if (patch_cap < seed_count) {
        stats->reason = "patch-face-budget";
        return -1;
    }

    uint8_t *selected = (uint8_t *)ARENA_ALLOC(
        arena, mesh->nf * sizeof(*selected));
    uint8_t *grown = (uint8_t *)ARENA_ALLOC(
        arena, mesh->nf * sizeof(*grown));
    memcpy(selected, seed, mesh->nf * sizeof(*selected));
    int min_rings = config->min_rings;
    int max_rings = config->max_rings;
    if (min_rings < 0) min_rings = 0;
    if (max_rings < min_rings) max_rings = min_rings;
    grow_mask_same_label(arena, mesh->nf, face_labels,
                         adj_fa, adj_fb, n_adj, min_rings, selected);
    size_t selected_count = 0;
    for (size_t fi = 0; fi < mesh->nf; fi++)
        if (selected[fi]) selected_count++;
    uint8_t *collar_absorb = (uint8_t *)ARENA_ALLOC(
        arena, mesh->nf * sizeof(*collar_absorb));
    const char *last_reason = "no-ring-cleared";

    for (int rings = min_rings; rings <= max_rings; rings++) {
        for (int collar_pass = 0; collar_pass <= 24; collar_pass++) {
            if (selected_count > patch_cap) {
                last_reason = "patch-face-budget";
                break;
            }
            double diameter = mask_diameter(mesh, selected);
            stats->patch_diameter = diameter;
            if (diameter > config->max_patch_diameter) {
                last_reason = "patch-diameter-budget";
                break;
            }

            Arena_Mark attempt = Arena_save(arena);
            ComponentMesh candidate;
            memset(&candidate, 0, sizeof(candidate));
            OverlapMergeStats candidate_stats = *stats;
            const char *reason = "construction-failed";
            float common_normal[3] = {0.0f, 0.0f, 0.0f};
            float common_center[3] = {0.0f, 0.0f, 0.0f};
            double common_plane_rms = 0.0;
            int plane_ok = 1;
            if (fit_mask_tangent_plane(arena, mesh, selected,
                                       common_normal, common_center,
                                       &common_plane_rms) != 0) {
                reason = "local-tangent-plane-fit-failed";
                plane_ok = 0;
            } else if (common_plane_rms > config->max_plane_rms) {
                reason = "local-overlap-event-not-planar";
                plane_ok = 0;
            }
            fprintf(stderr,
                    "      overlap patch ring %d collar %d: plane_rms=%.3f "
                    "n=(%.3f,%.3f,%.3f)\n",
                    rings, collar_pass, common_plane_rms,
                    common_normal[0], common_normal[1], common_normal[2]);
            candidate_stats.plane_rms = common_plane_rms;
            candidate_stats.plane_normal[0] = common_normal[0];
            candidate_stats.plane_normal[1] = common_normal[1];
            candidate_stats.plane_normal[2] = common_normal[2];
            int built = -1;
            memset(collar_absorb, 0,
                   mesh->nf * sizeof(*collar_absorb));
            if (plane_ok) {
                built = build_candidate(
                    arena, mesh, selected, selected_count, face_labels,
                    labels[0], labels[1], common_normal, common_center,
                    common_plane_rms, config, &input_shape, collar_absorb,
                    &candidate, &candidate_stats, &reason);
            }
            if (built == 0) {
                candidate_stats.rings = rings;
                candidate_stats.patch_diameter = diameter;
                candidate_stats.reason = "topology-accepted";
                *out = candidate;
                out->self = out;
                *stats = candidate_stats;
                return 0;
            }

            size_t absorb_count = 0;
            if (strcmp(reason, "collar-frontier-double-covered") == 0 ||
                strcmp(reason, "collar-vertex-pinch") == 0 ||
                strcmp(reason, "collar-boundary-hole") == 0) {
                for (size_t fi = 0; fi < mesh->nf; fi++)
                    if (collar_absorb[fi] && !selected[fi]) absorb_count++;
            }
            fprintf(stderr,
                    "      delamination ring %d collar %d: selected=%zu "
                    "diameter=%.2f -> %s; absorb=%zu\n",
                    rings, collar_pass, selected_count, diameter, reason,
                    absorb_count);
            last_reason = reason;
            Arena_restore(arena, attempt);

            if (absorb_count > 0 && collar_pass < 24) {
                for (size_t fi = 0; fi < mesh->nf; fi++) {
                    if (collar_absorb[fi] && !selected[fi]) {
                        selected[fi] = 1;
                        selected_count++;
                    }
                }
                continue;
            }
            break;
        }

        if (strcmp(last_reason, "patch-face-budget") == 0 ||
            strcmp(last_reason, "patch-diameter-budget") == 0)
            break;

        if (rings == max_rings) break;
        memcpy(grown, selected, mesh->nf * sizeof(*grown));
        for (size_t i = 0; i < n_adj; i++) {
            int32_t fa = adj_fa[i], fb = adj_fb[i];
            if (face_labels[fa] != face_labels[fb]) continue;
            if (selected[fa]) grown[fb] = 1;
            if (selected[fb]) grown[fa] = 1;
        }
        size_t next_count = 0;
        for (size_t fi = 0; fi < mesh->nf; fi++) {
            if (grown[fi]) next_count++;
        }
        if (next_count == selected_count) {
            last_reason = "ring-growth-stalled";
            break;
        }
        memcpy(selected, grown, mesh->nf * sizeof(*selected));
        selected_count = next_count;
    }

    stats->reason = last_reason;
    return -1;
}
