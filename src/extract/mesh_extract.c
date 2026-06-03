#include "../common/ves_platform.h"

#include "mesh_extract.h"
#include "../common/arena.h"
#include "../common/tiff_io.h"
#include "../common/halo_loader.h"
#include "../common/union_find.h"
#include "../common/pca.h"
#include "../common/mls_project.h"
#include "../common/vert_weld.h"
#include "../common/obj_io.h"
#include "../common/dump_obj.h"
#include "../common/pipeline_constants.h"
#include "../common/mesh_types.h"
#include "../common/qem.h"
#include "../common/obj_colors.h"
#include "../common/csr.h"
#include "../remesh/ball_pivot.h"
#include "../remesh/orient_mesh.h"
#include "../remesh/normal_orient.h"

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>

/* 6-connectivity neighbor offsets (dz, dy, dx) */
static const int NBR6[6][3] = {
    {-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}
};

/* ================================================================
 * Phase 2: 3D Connected Components via BFS
 * ================================================================ */

static int comp_info_cmp_desc(const void *a, const void *b)
{
    const CompInfo *ca = (const CompInfo *)a;
    const CompInfo *cb = (const CompInfo *)b;
    if (cb->count > ca->count) return 1;
    if (cb->count < ca->count) return -1;
    return 0;
}

int cc_label_3d(Arena_T arena, const uint8_t *vol,
                int D, int H, int W,
                int32_t *labels,
                CompInfo **out_comps, int32_t *out_n_comps)
{
    size_t vol_size = (size_t)D * (size_t)H * (size_t)W;
    int32_t HW = H * W;

    int32_t *queue = (int32_t *)ARENA_ALLOC(arena,
                                            (long)(vol_size * sizeof(int32_t)));
    int32_t max_comps = 4096;
    CompInfo *comps = (CompInfo *)ARENA_ALLOC(arena,
                                              (long)((size_t)max_comps * sizeof(CompInfo)));
    int32_t n_comps = 0;

    for (int32_t idx = 0; (size_t)idx < vol_size; idx++) {
        if (vol[idx] == 0 || labels[idx] != 0) continue;

        int32_t comp_label = n_comps + 1;
        int32_t qh = 0, qt = 0;
        queue[qt++] = idx;
        labels[idx] = comp_label;
        int32_t count = 0;

        while (qh < qt) {
            int32_t cur = queue[qh++];
            count++;

            int z = cur / HW;
            int rem = cur - z * HW;
            int y = rem / W;
            int x = rem - y * W;

            for (int n = 0; n < 6; n++) {
                int nz = z + NBR6[n][0];
                int ny = y + NBR6[n][1];
                int nx = x + NBR6[n][2];
                if (nz < 0 || nz >= D || ny < 0 || ny >= H || nx < 0 || nx >= W)
                    continue;
                int32_t ni = (int32_t)((size_t)nz * (size_t)HW +
                                       (size_t)ny * (size_t)W + (size_t)nx);
                if (vol[ni] == 0 || labels[ni] != 0) continue;
                labels[ni] = comp_label;
                queue[qt++] = ni;
            }
        }

        if (n_comps >= max_comps) {
            int32_t new_max = max_comps * 2;
            CompInfo *nc = (CompInfo *)ARENA_ALLOC(arena,
                                                    (long)((size_t)new_max * sizeof(CompInfo)));
            memcpy(nc, comps, (size_t)n_comps * sizeof(CompInfo));
            comps = nc;
            max_comps = new_max;
        }
        comps[n_comps].label = comp_label;
        comps[n_comps].count = count;
        n_comps++;
    }

    *out_comps = comps;
    *out_n_comps = n_comps;
    return 0;
}

/* Propagate MLS-normal sign consistency through mesh-edge BFS so adjacent
 * vertices share the same sign. Intended to eliminate the diagonal tears
 * the (1,1,1)-hemisphere bias in mls_project.c:382-400 leaves on curved
 * sheets where the local sheet normal crosses the (1,1,1)-perpendicular
 * plane. Algorithm: BFS from the smallest-(z,y,x) seed; whenever we cross
 * an edge (u, v) for the first time, flip n_v if dot(n_u, n_v) < 0.
 *
 * STATUS (2026-05-27): EXPERIMENTAL, DISABLED BY DEFAULT in the call
 * site. The first deployment caught a regression: for connected components
 * thicker than the LOP kernel radius (R=12 vox), the kernel cannot span
 * the through-thickness, so front- and back-side vertices have genuinely
 * opposite local-outward directions. BFS forces them to agree, which
 * makes backface_cull_per_vert keep both sides of the envelope. On the
 * z04480_y03328 + x02816/x02944 2-cube test this doubled one component's
 * face count (50,851 -> 100,820) and added 403 open-boundary edges to the
 * welded grid (4443 -> 4846) with no measurable tear reduction. Enable
 * via the SIGN_PROPAGATE env var for further experimentation; a different
 * cull strategy (e.g. spatial pairing of co-located opposite-winding
 * faces) is probably the right fix.
 *
 * Returns the number of vertices whose sign was flipped. */
static size_t propagate_normal_signs(Arena_T arena,
                                     const float *verts, size_t nv,
                                     const int32_t *faces, size_t nf,
                                     float *normals)
{
    if (nv == 0 || nf == 0) return 0;

    Arena_Mark mark = Arena_save(arena);

    CSR_T adj = CSR_from_faces(arena, faces, nf, nv);
    const int32_t *off = CSR_offset(adj);
    const int32_t *tgt = CSR_target(adj);

    uint8_t *visited = (uint8_t *)ARENA_CALLOC(arena, (long)nv, 1L);
    int32_t *queue = (int32_t *)ARENA_ALLOC(arena,
                        (long)nv * (long)sizeof(int32_t));

    size_t n_flips = 0;

    for (;;) {
        /* Seed: smallest (z, y, x) among unvisited verts with a defined
         * normal. Halo-deterministic — bit-identical position triples in
         * cubes A and B compare the same way. Verts with undefined MLS
         * normal are marked visited as we scan past them. */
        int32_t seed = -1;
        float best_z = 0.0f, best_y = 0.0f, best_x = 0.0f;
        for (size_t v = 0; v < nv; v++) {
            if (visited[v]) continue;
            float nz = normals[v * 3 + 0];
            float ny = normals[v * 3 + 1];
            float nx = normals[v * 3 + 2];
            float nmag2 = nz*nz + ny*ny + nx*nx;
            if (nmag2 < 1e-18f) {
                visited[v] = 1;
                continue;
            }
            float z = verts[v * 3 + 0];
            float y = verts[v * 3 + 1];
            float x = verts[v * 3 + 2];
            if (seed < 0 ||
                z <  best_z ||
                (z == best_z && y <  best_y) ||
                (z == best_z && y == best_y && x < best_x)) {
                seed = (int32_t)v;
                best_z = z;
                best_y = y;
                best_x = x;
            }
        }
        if (seed < 0) break;

        visited[seed] = 1;
        int32_t q_head = 0;
        int32_t q_tail = 0;
        queue[q_tail++] = seed;

        while (q_head < q_tail) {
            int32_t u = queue[q_head++];
            float nuz = normals[u * 3 + 0];
            float nuy = normals[u * 3 + 1];
            float nux = normals[u * 3 + 2];

            for (int32_t j = off[u]; j < off[u + 1]; j++) {
                int32_t w = tgt[j];
                if (visited[w]) continue;

                float nwz = normals[w * 3 + 0];
                float nwy = normals[w * 3 + 1];
                float nwx = normals[w * 3 + 2];
                float wmag2 = nwz*nwz + nwy*nwy + nwx*nwx;
                if (wmag2 < 1e-18f) {
                    visited[w] = 1;
                    continue;
                }

                float d = nuz*nwz + nuy*nwy + nux*nwx;
                if (d < 0.0f) {
                    normals[w * 3 + 0] = -nwz;
                    normals[w * 3 + 1] = -nwy;
                    normals[w * 3 + 2] = -nwx;
                    n_flips++;
                }
                visited[w] = 1;
                queue[q_tail++] = w;
            }
        }
    }

    Arena_restore(arena, mark);
    return n_flips;
}

/* ---- Spatial-pairing cull helpers ---- */

typedef struct {
    uint64_t key;
    int32_t  face_idx;
    float    cz, cy, cx;   /* centroid (kept for tie-break + radius check) */
} PairHashEntry;

static int cmp_pair_hash(const void *a, const void *b)
{
    const PairHashEntry *ea = (const PairHashEntry *)a;
    const PairHashEntry *eb = (const PairHashEntry *)b;
    if (ea->key < eb->key) return -1;
    if (ea->key > eb->key) return 1;
    /* Tiebreak by position then face_idx — keeps the order deterministic
     * across cubes (same physical centroid sorts identically regardless of
     * each cube's face indexing). */
    if (ea->cz < eb->cz) return -1;
    if (ea->cz > eb->cz) return 1;
    if (ea->cy < eb->cy) return -1;
    if (ea->cy > eb->cy) return 1;
    if (ea->cx < eb->cx) return -1;
    if (ea->cx > eb->cx) return 1;
    if (ea->face_idx < eb->face_idx) return -1;
    if (ea->face_idx > eb->face_idx) return 1;
    return 0;
}

static inline uint64_t pair_cell_key(int32_t cz, int32_t cy, int32_t cx)
{
    uint64_t uz = (uint64_t)(cz + (1 << 20)) & 0x1FFFFF;
    uint64_t uy = (uint64_t)(cy + (1 << 20)) & 0x1FFFFF;
    uint64_t ux = (uint64_t)(cx + (1 << 20)) & 0x1FFFFF;
    return (uz << 42) | (uy << 21) | ux;
}

/* Spatial-pairing cull: thin the MC double envelope by finding co-located
 * opposite-winding face pairs and dropping one of each pair. Replaces
 * backface_cull_per_vert's directional dot test.
 *
 * For each face f, compute centroid and unit normal. Build a uniform
 * spatial hash on centroids (cell size = PAIRING_RADIUS_VOX). For each
 * face, search the 27-cell ball around its centroid for any face within
 * PAIRING_RADIUS_VOX whose normal dots more negative than
 * PAIRING_DOT_THRESHOLD — that is the spatial twin. Tie-break: keep the
 * face whose centroid is lex-larger by (z, y, x); drop the other. Faces
 * with no twin (single-sided sheets, closed manifolds thicker than
 * PAIRING_RADIUS_VOX, slab side-walls perpendicular to closing pairs) are
 * kept.
 *
 * Why no global sign convention: the predicate is symmetric in the two
 * faces of a pair — neither needs to point "outward" by any external
 * reference — so curved sheets cannot tear at any direction-bias
 * boundary. The dropped face is whichever of the pair has the smaller
 * centroid, which is a property of the geometry, not of any chosen
 * reference frame.
 *
 * Halo-determinism: cells are world-aligned via `cell_origin` (so
 * adjacent cubes hash the same physical centroid to the same cell);
 * sort tiebreak is by position (so bit-identical centroids in cubes A
 * and B sort identically regardless of face indexing); the keep/drop
 * decision is by position. Same physical pair → same drop choice in
 * both cubes.
 *
 * Compacts faces in-place; returns the kept-face count. */
static size_t spatial_pairing_cull(Arena_T arena,
                                   const float *verts, size_t nv,
                                   int32_t *faces, size_t nf,
                                   const float cell_origin[3])
{
    if (nv == 0 || nf == 0) return nf;

    Arena_Mark mark = Arena_save(arena);

    float R = PAIRING_RADIUS_VOX;
    float R2 = R * R;
    float inv_R = 1.0f / R;

    float ocz = cell_origin ? cell_origin[0] : 0.0f;
    float ocy = cell_origin ? cell_origin[1] : 0.0f;
    float ocx = cell_origin ? cell_origin[2] : 0.0f;

    /* Per-face centroid + unit normal. Store inline in the hash entries
     * for tiebreak + radius check; also keep a flat normals[] for fast
     * dot-product lookup by face index. */
    float *centroids = (float *)ARENA_ALLOC(arena,
                          (long)nf * 3L * (long)sizeof(float));
    float *normals = (float *)ARENA_ALLOC(arena,
                        (long)nf * 3L * (long)sizeof(float));
    PairHashEntry *entries = (PairHashEntry *)ARENA_ALLOC(arena,
                              (long)nf * (long)sizeof(PairHashEntry));
    uint8_t *keep = (uint8_t *)ARENA_ALLOC(arena, (long)nf);

    size_t n_degenerate = 0;
    for (size_t f = 0; f < nf; f++) {
        int32_t a = faces[f*3+0];
        int32_t b = faces[f*3+1];
        int32_t c = faces[f*3+2];
        float az = verts[a*3+0], ay = verts[a*3+1], ax = verts[a*3+2];
        float bz = verts[b*3+0], by = verts[b*3+1], bx = verts[b*3+2];
        float cz = verts[c*3+0], cy = verts[c*3+1], cx = verts[c*3+2];

        float cenZ = (az + bz + cz) * (1.0f / 3.0f);
        float cenY = (ay + by + cy) * (1.0f / 3.0f);
        float cenX = (ax + bx + cx) * (1.0f / 3.0f);
        centroids[f*3+0] = cenZ;
        centroids[f*3+1] = cenY;
        centroids[f*3+2] = cenX;

        float abz_ = bz - az, aby_ = by - ay, abx_ = bx - ax;
        float acz_ = cz - az, acy_ = cy - ay, acx_ = cx - ax;
        float nz = aby_ * acx_ - abx_ * acy_;
        float ny = abx_ * acz_ - abz_ * acx_;
        float nx = abz_ * acy_ - aby_ * acz_;
        float nmag = sqrtf(nz*nz + ny*ny + nx*nx);
        if (nmag > 1e-9f) {
            float inv = 1.0f / nmag;
            normals[f*3+0] = nz * inv;
            normals[f*3+1] = ny * inv;
            normals[f*3+2] = nx * inv;
            keep[f] = 1;
        } else {
            normals[f*3+0] = 0.0f;
            normals[f*3+1] = 0.0f;
            normals[f*3+2] = 0.0f;
            keep[f] = 0;    /* drop degenerate face */
            n_degenerate++;
        }

        entries[f].cz = cenZ;
        entries[f].cy = cenY;
        entries[f].cx = cenX;
        entries[f].face_idx = (int32_t)f;
        entries[f].key = pair_cell_key(
            (int32_t)floorf((cenZ + ocz) * inv_R),
            (int32_t)floorf((cenY + ocy) * inv_R),
            (int32_t)floorf((cenX + ocx) * inv_R));
    }

    qsort(entries, nf, sizeof(PairHashEntry), cmp_pair_hash);

    /* Build parallel key array for binary-search lookup. */
    uint64_t *keys = (uint64_t *)ARENA_ALLOC(arena,
                        (long)nf * (long)sizeof(uint64_t));
    for (size_t i = 0; i < nf; i++) keys[i] = entries[i].key;

    /* For each face, search 27 neighbour cells for a co-located opposite-
     * winding face whose centroid is lex-greater than this one's. If any
     * exists, drop this face (the higher-centroid one of each pair wins).
     * If no such "winning twin" exists, keep this face. */
    for (size_t f = 0; f < nf; f++) {
        if (!keep[f]) continue;    /* already-dropped degenerate */

        float fz = centroids[f*3+0];
        float fy = centroids[f*3+1];
        float fx = centroids[f*3+2];
        float fnz = normals[f*3+0];
        float fny = normals[f*3+1];
        float fnx = normals[f*3+2];

        int32_t ciz = (int32_t)floorf((fz + ocz) * inv_R);
        int32_t ciy = (int32_t)floorf((fy + ocy) * inv_R);
        int32_t cix = (int32_t)floorf((fx + ocx) * inv_R);

        int dropped = 0;
        for (int dz = -1; dz <= 1 && !dropped; dz++) {
            for (int dy = -1; dy <= 1 && !dropped; dy++) {
                for (int dx = -1; dx <= 1 && !dropped; dx++) {
                    uint64_t target = pair_cell_key(ciz+dz, ciy+dy, cix+dx);
                    size_t lo = 0, hi = nf;
                    while (lo < hi) {
                        size_t mid = lo + (hi - lo) / 2;
                        if (keys[mid] < target) lo = mid + 1;
                        else hi = mid;
                    }
                    for (size_t k = lo; k < nf && keys[k] == target; k++) {
                        int32_t f2 = entries[k].face_idx;
                        if ((size_t)f2 == f) continue;
                        if (!keep[f2]) continue;    /* degenerate twin */

                        float t2z = entries[k].cz;
                        float t2y = entries[k].cy;
                        float t2x = entries[k].cx;
                        float dz2 = t2z - fz;
                        float dy2 = t2y - fy;
                        float dx2 = t2x - fx;
                        if (dz2*dz2 + dy2*dy2 + dx2*dx2 >= R2) continue;

                        float d = fnz * normals[f2*3+0]
                                + fny * normals[f2*3+1]
                                + fnx * normals[f2*3+2];
                        if (d > PAIRING_DOT_THRESHOLD) continue;

                        /* Opposite-winding twin in range. Tie-break:
                         * keep the face with the larger (z, y, x) centroid;
                         * fall back to face_idx if positions coincide bit-
                         * exact (degenerate, shouldn't happen for distinct
                         * MC triangles). */
                        int twin_wins = 0;
                        if (t2z > fz)                            twin_wins = 1;
                        else if (t2z == fz && t2y > fy)          twin_wins = 1;
                        else if (t2z == fz && t2y == fy && t2x > fx) twin_wins = 1;
                        else if (t2z == fz && t2y == fy && t2x == fx
                                 && f2 > (int32_t)f)             twin_wins = 1;

                        if (twin_wins) {
                            keep[f] = 0;
                            dropped = 1;
                            break;
                        }
                        /* else: f wins this pair; keep looking — f might
                         * have another twin that beats it. */
                    }
                }
            }
        }
    }

    /* Compact faces in place. */
    size_t kept = 0;
    for (size_t f = 0; f < nf; f++) {
        if (!keep[f]) continue;
        if (kept != f) {
            faces[kept*3+0] = faces[f*3+0];
            faces[kept*3+1] = faces[f*3+1];
            faces[kept*3+2] = faces[f*3+2];
        }
        kept++;
    }

    (void)n_degenerate;

    Arena_restore(arena, mark);
    return kept;
}

/* ================================================================
 * Phase 5: Cleanup — vertex compaction, UF split
 * ================================================================ */

/* Remove unreferenced vertices, remap face indices.
 * Allocates new vertex array from arena; modifies faces in-place. */
static void clean_unreferenced(Arena_T arena,
                               const float *old_verts, size_t old_nv,
                               int32_t *faces, size_t nf,
                               float **new_verts, size_t *new_nv)
{
    int32_t *remap = (int32_t *)ARENA_CALLOC(arena, (long)old_nv,
                                              (long)sizeof(int32_t));
    /* Mark used vertices */
    for (size_t i = 0; i < nf * 3; i++) {
        assert(faces[i] >= 0 && (size_t)faces[i] < old_nv);
        remap[faces[i]] = 1;
    }

    /* Build remap: remap[old] = new index */
    int32_t count = 0;
    for (size_t i = 0; i < old_nv; i++) {
        if (remap[i]) {
            remap[i] = count++;
        } else {
            remap[i] = -1;
        }
    }

    /* Allocate compact vertex array */
    float *verts = (float *)ARENA_ALLOC(arena,
                                         (long)((size_t)count * 3 * sizeof(float)));
    for (size_t i = 0; i < old_nv; i++) {
        if (remap[i] >= 0) {
            verts[remap[i] * 3 + 0] = old_verts[i * 3 + 0];
            verts[remap[i] * 3 + 1] = old_verts[i * 3 + 1];
            verts[remap[i] * 3 + 2] = old_verts[i * 3 + 2];
        }
    }

    /* Remap face indices */
    for (size_t i = 0; i < nf * 3; i++) {
        faces[i] = remap[faces[i]];
    }

    *new_verts = verts;
    *new_nv = (size_t)count;
}

/* Union-Find mesh split: remove sub-components < MIN_FRAGMENT_FACES.
 * Returns new face count; compacts faces in-place. */
static size_t mesh_split_remove_small(Arena_T arena,
                                      size_t nv,
                                      int32_t *faces, size_t nf)
{
    if (nf == 0) return 0;

    UnionFind uf = UF_new(arena, (int32_t)nv);

    for (size_t i = 0; i < nf; i++) {
        int32_t a = faces[i * 3 + 0];
        int32_t b = faces[i * 3 + 1];
        int32_t c = faces[i * 3 + 2];
        uf_union(&uf, a, b);
        uf_union(&uf, b, c);
        uf_union(&uf, a, c);
    }

    /* Count faces per UF root */
    int32_t *face_count = (int32_t *)ARENA_CALLOC(arena, (long)nv,
                                                    (long)sizeof(int32_t));
    for (size_t i = 0; i < nf; i++) {
        int32_t root = uf_find(&uf, faces[i * 3]);
        face_count[root]++;
    }

    /* Keep faces whose component has >= MIN_FRAGMENT_FACES */
    size_t kept = 0;
    for (size_t i = 0; i < nf; i++) {
        int32_t root = uf_find(&uf, faces[i * 3]);
        if (face_count[root] >= MIN_FRAGMENT_FACES) {
            faces[kept * 3 + 0] = faces[i * 3 + 0];
            faces[kept * 3 + 1] = faces[i * 3 + 1];
            faces[kept * 3 + 2] = faces[i * 3 + 2];
            kept++;
        }
    }
    return kept;
}

/* ================================================================
 * MeshExtract_run — Main entry point
 * ================================================================ */

int MeshExtract_run(Arena_T          arena,
                    const char      *tiff_path,
                    const char      *pred_dir,
                    int              halo_voxels,
                    size_t           cube_D,
                    size_t           cube_H,
                    size_t           cube_W,
                    int              n_threads,
                    const char      *dump_cube_dir,
                    const char      *cube_id,
                    int              skip_qem,
                    ComponentMesh  **out_meshes,
                    size_t          *out_n_meshes,
                    MeshResplitCloud **out_clouds)
{
    assert(arena);
    assert(out_meshes && out_n_meshes);
    (void)cube_H; (void)cube_W;
    (void)n_threads;

    *out_meshes = NULL;
    *out_n_meshes = 0;
    if (out_clouds) *out_clouds = NULL;

    if (halo_voxels < 0) {
        fprintf(stderr, "MeshExtract: halo_voxels=%d must be >= 0\n",
                halo_voxels);
        return -1;
    }

#ifdef VESUVIUS_DEBUG
    double t_start, t_phase;
    t_start = ves_clock_sec();
#endif

    /* Phase 1: Load volume (with halo if requested). */
    uint8_t *vol = NULL;
    int D = 0, H = 0, W = 0;
    /* Cube origin in world voxel coords (cube-local (0,0,0) → world). Used
     * later by MLS_project_verts as cell_origin so spatial-hash cells align
     * across adjacent cubes (same physical point gets the same cell key in
     * cube A and cube B). Zero for non-halo mode. */
    float cube_world_origin[3] = {0.0f, 0.0f, 0.0f};
    if (halo_voxels > 0) {
        if (!pred_dir || !cube_id) {
            fprintf(stderr, "MeshExtract: halo requires pred_dir and cube_id\n");
            return -1;
        }
        int p_size = 0;
        int64_t origin[3] = {0, 0, 0};
        if (HaloLoader_load(arena, pred_dir, cube_id,
                            (int)cube_D, halo_voxels,
                            &vol, &p_size, origin) != 0) {
            fprintf(stderr, "MeshExtract: HaloLoader_load failed for %s\n",
                    cube_id);
            return -1;
        }
        D = H = W = p_size;
        /* HaloLoader returns the padded-buffer origin (cube_origin - halo).
         * cube-local (0,0,0) is at world = origin + halo_voxels. */
        cube_world_origin[0] = (float)(origin[0] + halo_voxels);
        cube_world_origin[1] = (float)(origin[1] + halo_voxels);
        cube_world_origin[2] = (float)(origin[2] + halo_voxels);
    } else {
        assert(tiff_path);
        if (TiffIO_load(arena, tiff_path, &vol, &D, &H, &W) != 0) {
            fprintf(stderr, "MeshExtract: failed to load TIFF: %s\n",
                    tiff_path);
            return -1;
        }
    }

    size_t vol_size = (size_t)D * (size_t)H * (size_t)W;
    if (vol_size == 0) return 0;

    /* Threshold to binary */
    for (size_t i = 0; i < vol_size; i++) {
        vol[i] = (vol[i] > 0) ? 1 : 0;
    }

#ifdef VESUVIUS_DEBUG
    t_phase = ves_clock_sec();
    fprintf(stderr, "  Extract: loaded TIFF %dx%dx%d (%.3f s)\n", D, H, W,
            t_phase - t_start);
#endif

    /* Phase 2: 3D connected components */
    int32_t *labels = (int32_t *)ARENA_CALLOC(arena, (long)vol_size,
                                               (long)sizeof(int32_t));
    CompInfo *comps = NULL;
    int32_t n_raw_comps = 0;
    cc_label_3d(arena, vol, D, H, W, labels, &comps, &n_raw_comps);

    /* Filter < MIN_CC_SIZE, sort descending, keep top MAX_COMPONENTS */
    int32_t n_valid = 0;
    for (int32_t i = 0; i < n_raw_comps; i++) {
        if (comps[i].count >= MIN_CC_SIZE) {
            comps[n_valid++] = comps[i];
        }
    }
    qsort(comps, (size_t)n_valid, sizeof(CompInfo), comp_info_cmp_desc);
    if (n_valid > MAX_COMPONENTS) n_valid = MAX_COMPONENTS;

#ifdef VESUVIUS_DEBUG
    {
        double t_cc = ves_clock_sec();
        fprintf(stderr, "  Extract: %d raw CCs, %d after filter (min=%d) (%.3f s)\n",
                n_raw_comps, n_valid, MIN_CC_SIZE,
                t_cc - t_start);
    }
#endif

    if (n_valid == 0) return 0;

    /* Allocate output array */
    ComponentMesh *meshes = (ComponentMesh *)ARENA_ALLOC(arena,
                                                          (long)((size_t)n_valid * sizeof(ComponentMesh)));
    memset(meshes, 0, (size_t)n_valid * sizeof(ComponentMesh));

    /* Optional per-component raw+LOP cloud snapshot for the connectivity
     * re-split (MeshResplit). Index-aligned with `meshes`; filled at finalize. */
    MeshResplitCloud *clouds = NULL;
    if (out_clouds) {
        clouds = (MeshResplitCloud *)ARENA_ALLOC(arena,
                     (long)((size_t)n_valid * sizeof(MeshResplitCloud)));
        memset(clouds, 0, (size_t)n_valid * sizeof(MeshResplitCloud));
    }

    /* Scratch arena for per-component temporaries */
    Arena_T scratch = Arena_new();
    size_t out_count = 0;
    int HW = H * W;

    /* Halo bounds in vol coords: owned region = [halo_voxels, halo_voxels + cube_D).
     * When halo == 0, owned region = whole volume. */
    int owned_lo = halo_voxels;
    int owned_hi = halo_voxels + (int)cube_D;

    for (int32_t ci = 0; ci < n_valid; ci++) {
        Arena_free(scratch);
        int32_t comp_label = comps[ci].label;
        /* Connectivity-resplit snapshot (main arena: scratch is reset next iter,
         * but the resplit runs after all components return). Filled after LOP. */
        float *snap_orig = NULL, *snap_lop = NULL; size_t snap_n = 0;

#ifdef VESUVIUS_DEBUG
        double t_c0, t_c1;
        t_c0 = ves_clock_sec();
#endif

        /* Find bounding box of this component */
        int zmin = D, zmax = -1, ymin = H, ymax = -1, xmin = W, xmax = -1;
        for (int z = 0; z < D; z++) {
            for (int y = 0; y < H; y++) {
                for (int x = 0; x < W; x++) {
                    if (labels[(size_t)z * (size_t)HW +
                               (size_t)y * (size_t)W + (size_t)x] == comp_label) {
                        if (z < zmin) zmin = z;
                        if (z > zmax) zmax = z;
                        if (y < ymin) ymin = y;
                        if (y > ymax) ymax = y;
                        if (x < xmin) xmin = x;
                        if (x > xmax) xmax = x;
                    }
                }
            }
        }
        if (zmax < 0) continue;

        /* When halo_voxels > 0, drop components whose bounding box lies
         * entirely outside the owned region. These belong to neighbor cubes
         * and would be emitted by them — running the full pipeline on them
         * here wastes compute and the centroid trim would drop them anyway. */
        if (halo_voxels > 0) {
            int overlap_z = (zmax >= owned_lo) && (zmin < owned_hi);
            int overlap_y = (ymax >= owned_lo) && (ymin < owned_hi);
            int overlap_x = (xmax >= owned_lo) && (xmin < owned_hi);
            if (!overlap_z || !overlap_y || !overlap_x) {
#ifdef VESUVIUS_DEBUG
                fprintf(stderr,
                    "    comp %d: bbox [%d-%d,%d-%d,%d-%d] outside owned"
                    " [%d,%d) -> halo-only, drop\n",
                    comp_label, zmin, zmax, ymin, ymax, xmin, xmax,
                    owned_lo, owned_hi);
#endif
                continue;
            }
        }

#ifdef VESUVIUS_DEBUG
        t_c1 = ves_clock_sec();
        fprintf(stderr, "    comp %d: bbox [%d-%d, %d-%d, %d-%d] (%.3f s)\n",
                comp_label, zmin, zmax, ymin, ymax, xmin, xmax,
                t_c1 - t_c0);
#endif

        /* Extract component mask into tight bounding box, then pad with
         * 1 voxel of zeros on all 6 faces (matching Python np.pad approach).
         * The 1-voxel border gives the LOP/BPA neighbourhood a defined
         * value just outside the bbox so the projection and pivot behave
         * the same at the bbox edge as in the interior. */
        int bD = zmax - zmin + 1;
        int bH = ymax - ymin + 1;
        int bW = xmax - xmin + 1;

        /* Padded dimensions: tight bbox + 1 voxel border on each side */
        int pD = bD + 2;
        int pH = bH + 2;
        int pW = bW + 2;

        /* Allocate zero-filled padded subvolume (scratch) */
        size_t pad_size = (size_t)pD * (size_t)pH * (size_t)pW;
        uint8_t *padded = (uint8_t *)ARENA_CALLOC(scratch, (long)pad_size, 1L);

        /* Fill the padded volume covering [zmin-1, zmax+1] x ... so that
         * MC at the bbox border sees the right value:
         *
         *  - If the border voxel is INSIDE the working volume, sample
         *    labels[] there. This matters when the same component's
         *    voxels actually extend one step past the tight bbox (which
         *    is impossible by definition of the tight bbox, so this case
         *    contributes nothing new but is correct).
         *  - If the border voxel is OUTSIDE the working volume, REPLICATE
         *    the inner-edge value (i.e. set it equal to the labeled
         *    value at the nearest in-bounds voxel). This is the
         *    cross-cube-determinism fix: zero-padding at the working-
         *    volume edge generated phantom MC closing triangles that
         *    the neighbor cube's MC didn't produce (the neighbor saw
         *    real data on the other side), so meshes diverged at every
         *    seam. Replicating instead leaves the surface OPEN at the
         *    working-volume edge; the neighboring cube's owned-region
         *    mesh supplies that surface and grid_weld merges seams.
         */
        for (int z = zmin - 1; z <= zmax + 1; z++) {
            for (int y = ymin - 1; y <= ymax + 1; y++) {
                for (int x = xmin - 1; x <= xmax + 1; x++) {
                    int iz = z, iy = y, ix = x;
                    if (iz < 0) iz = 0; else if (iz >= D) iz = D - 1;
                    if (iy < 0) iy = 0; else if (iy >= H) iy = H - 1;
                    if (ix < 0) ix = 0; else if (ix >= W) ix = W - 1;
                    if (labels[(size_t)iz * (size_t)HW +
                               (size_t)iy * (size_t)W +
                               (size_t)ix] == comp_label) {
                        size_t pi = (size_t)(z - zmin + 1) * (size_t)pH * (size_t)pW +
                                    (size_t)(y - ymin + 1) * (size_t)pW +
                                    (size_t)(x - xmin + 1);
                        padded[pi] = 1;
                    }
                }
            }
        }

#ifdef VESUVIUS_DEBUG
        t_c1 = ves_clock_sec();
        fprintf(stderr, "    comp %d: subvol %dx%dx%d (%.3f s)\n",
                comp_label, pD, pH, pW,
                t_c1 - t_c0);
#endif

        /* Source the point cloud by voxel-center sampling: one point at
         * the center of every foreground voxel. This gives a single-layer
         * cloud directly — no double envelope to undo via a cull. Faces
         * are built later by BPA. */
        float *surf_verts = NULL;
        int32_t *surf_faces = NULL;
        size_t surf_nv = 0, surf_nf = 0;
        {
            /* Count FG voxels in padded, then allocate and fill. */
            size_t fg = 0;
            for (size_t pi = 0; pi < pad_size; pi++) if (padded[pi]) fg++;
            surf_nv = fg;
            surf_nf = 0;
            if (fg > 0) {
                surf_verts = (float *)ARENA_ALLOC(scratch,
                                (long)(fg * 3 * sizeof(float)));
                size_t out_idx = 0;
                for (int pz = 0; pz < pD; pz++) {
                    for (int py = 0; py < pH; py++) {
                        for (int px = 0; px < pW; px++) {
                            size_t pi = (size_t)pz * (size_t)pH * (size_t)pW +
                                        (size_t)py * (size_t)pW + (size_t)px;
                            if (!padded[pi]) continue;
                            /* Voxel center at +0.5; the shift block below
                             * applies the cube-local offset so positions
                             * land at (vol_z + 0.5 - halo). */
                            surf_verts[out_idx * 3 + 0] = (float)pz + 0.5f;
                            surf_verts[out_idx * 3 + 1] = (float)py + 0.5f;
                            surf_verts[out_idx * 3 + 2] = (float)px + 0.5f;
                            out_idx++;
                        }
                    }
                }
            }
        }

#ifdef VESUVIUS_DEBUG
        t_c1 = ves_clock_sec();
        fprintf(stderr, "    comp %d: %s -> %zu verts, %zu faces (%.3f s)\n",
                comp_label, "voxel-centers",
                surf_nv, surf_nf, t_c1 - t_c0);
#endif

        if (surf_nv < (size_t)MLS_MIN_NEIGHBOURS) continue;

        /* Convert MC vertices from padded-local to cube-local coordinates.
         * Padding added +1 to all coords (subtract 1), then add bounding box
         * origin in vol space (zmin/ymin/xmin), then subtract halo_voxels to
         * shift the cube's owned region to [0, cube_D). When halo == 0 this
         * is identical to the legacy single-cube behavior. */
        for (size_t v = 0; v < surf_nv; v++) {
            surf_verts[v * 3 + 0] += (float)(zmin - 1 - halo_voxels);
            surf_verts[v * 3 + 1] += (float)(ymin - 1 - halo_voxels);
            surf_verts[v * 3 + 2] += (float)(xmin - 1 - halo_voxels);
        }

        /* Cross-cube MC determinism: in halo mode, drop verts whose
         * cube-local coord is at the WORKING-VOLUME EDGE. These verts
         * come from MC processing on the 1-voxel padding around the
         * bbox; with the replicate-edge padding fix the values look
         * consistent with the inner-edge but they're still based on
         * SAMPLED-NOT-REAL data at the working-volume boundary. The
         * neighbour cube has real data at those source positions and
         * is the authoritative source. Dropping these verts here
         * gives bit-exact cross-cube MC overlap (synthetic test:
         * 100% A-match, 99.23% B-match where the 0.77% B-only verts
         * are correct B-unique contributions A can't see).
         *
         * Working-volume range in cube-local coords is
         *   [-halo_voxels, cube_D + halo_voxels) for z, similarly y,x.
         * Verts at z <= -halo_voxels OR z >= cube_D + halo_voxels are
         * dropped (and their incident faces along with them). */
        if (halo_voxels > 0 && surf_nv > 0) {
            float lo = -(float)halo_voxels;
            float hi_z = (float)cube_D + (float)halo_voxels;
            float hi_y = (float)cube_H + (float)halo_voxels;
            float hi_x = (float)cube_W + (float)halo_voxels;
            int32_t *remap = (int32_t *)ARENA_ALLOC(scratch,
                                (long)surf_nv * (long)sizeof(int32_t));
            size_t nv_kept = 0, nv_dropped = 0;
            for (size_t v = 0; v < surf_nv; v++) {
                float zc = surf_verts[v * 3 + 0];
                float yc = surf_verts[v * 3 + 1];
                float xc = surf_verts[v * 3 + 2];
                if (zc <= lo || zc >= hi_z ||
                    yc <= lo || yc >= hi_y ||
                    xc <= lo || xc >= hi_x) {
                    remap[v] = -1;
                    nv_dropped++;
                } else {
                    remap[v] = (int32_t)nv_kept;
                    if (nv_kept != v) {
                        surf_verts[nv_kept * 3 + 0] = zc;
                        surf_verts[nv_kept * 3 + 1] = yc;
                        surf_verts[nv_kept * 3 + 2] = xc;
                    }
                    nv_kept++;
                }
            }
            if (nv_dropped > 0) {
                size_t nf_kept = 0;
                for (size_t f = 0; f < surf_nf; f++) {
                    int32_t v0 = surf_faces[f * 3 + 0];
                    int32_t v1 = surf_faces[f * 3 + 1];
                    int32_t v2 = surf_faces[f * 3 + 2];
                    if (remap[v0] < 0 || remap[v1] < 0 || remap[v2] < 0) continue;
                    surf_faces[nf_kept * 3 + 0] = remap[v0];
                    surf_faces[nf_kept * 3 + 1] = remap[v1];
                    surf_faces[nf_kept * 3 + 2] = remap[v2];
                    nf_kept++;
                }
#ifdef VESUVIUS_DEBUG
                fprintf(stderr,
                        "    comp %d: working-edge filter dropped %zu/%zu v, "
                        "%zu/%zu f\n",
                        comp_label, nv_dropped, surf_nv,
                        surf_nf - nf_kept, surf_nf);
#endif
                surf_nv = nv_kept;
                surf_nf = nf_kept;
                /* For voxel-center mode we have no faces; nf==0 is fine. */
                if (surf_nv == 0) continue;
            }
        }

        /* Stable palette color for this component, keyed off the original
         * component label so the same sheet keeps the same color through
         * every diagnostic stage and the final dump. */
        const float *diag_color =
            OBJ_COLORS[(comp_label - 1 + OBJ_NUM_COLORS) % OBJ_NUM_COLORS];

        /* MLS-midpoint collapse of the MC double envelope to its centerline
         * (LOP with μ=0, Lipman et al. SIGGRAPH 2007), then weld coincident
         * verts and drop degenerate triangles. Replaces the legacy
         * PCA-normal + backface_cull pair: that approach picked one face
         * of the envelope based on a single locally-estimated eigenvector,
         * which flipped between cubes whose vertex clouds differed even
         * slightly in the halo region. MLS-midpoint depends only on a
         * finite kernel radius (MLS_PROJECT_RADIUS_VOX) of neighbours, so
         * it is bit-deterministic across cube boundaries as long as the
         * halo width exceeds the kernel radius.
         *
         * See plan: the-main-problem-with-composed-badger.md */
        float *mls_verts = (float *)ARENA_ALLOC(scratch,
                              (long)(surf_nv * 3 * sizeof(float)));
        float *mls_normals = (float *)ARENA_ALLOC(scratch,
                                (long)(surf_nv * 3 * sizeof(float)));
        /* Ping-pong buffer for multi-iter LOP. MLS_project_verts reads from
         * `src` and writes to `dst`; if they were the same buffer, the
         * Gauss-Seidel-style read-after-write would make the centroid for
         * vertex i depend on whether vertices 0..i-1 have already been
         * updated in this iter — destroying halo determinism. Two buffers
         * + explicit swap gives clean Jacobi-style iteration.
         *
         * Normals are computed EVERY iter because MLS_project_verts
         * uses them to do tangent-plane projection (V' = V − ((V−C)·N)N)
         * — the projection direction is the local sheet normal, which
         * the function won't compute when out_normals is NULL. Passing
         * NULL would fall back to centroid replacement, which shrinks
         * the sheet in-plane with each iter (catastrophic at high iter
         * counts). The covariance pass roughly doubles per-iter cost
         * vs centroid-only but is essential for non-shrinking iteration. */
        float *mls_scratch_normals = (float *)ARENA_ALLOC(scratch,
                                        (long)(surf_nv * 3 * sizeof(float)));
        float *mls_scratch = (MLS_PROJECT_ITERS >= 2)
            ? (float *)ARENA_ALLOC(scratch,
                                   (long)(surf_nv * 3 * sizeof(float)))
            : NULL;
        const float *src = surf_verts;
        float *dst = mls_verts;
        for (int it = 0; it < MLS_PROJECT_ITERS; it++) {
            int is_last = (it == MLS_PROJECT_ITERS - 1);
            float *normals_dst = is_last ? mls_normals : mls_scratch_normals;
            MLS_project_verts(scratch, src, surf_nv,
                              MLS_PROJECT_RADIUS_VOX,
                              cube_world_origin,
                              dst, normals_dst);
            src = dst;
            dst = (dst == mls_verts) ? mls_scratch : mls_verts;
        }
        /* Final result lives in `src` (the last write target). If that
         * isn't mls_verts (even ITERS >= 2), copy it back. */
        if (src != mls_verts) {
            memcpy(mls_verts, src, surf_nv * 3 * sizeof(float));
        }

        /* Snapshot the raw voxel-center cloud (surf_verts, preserved as the LOP's
         * read-only src) and its 1:1 LOP image (mls_verts), BEFORE the weld
         * renumbers them, for the connectivity re-split. Main arena so it
         * outlives this component's scratch. */
        if (clouds) {
            snap_n = surf_nv;
            snap_orig = (float *)ARENA_ALLOC(arena, (long)(surf_nv * 3 * sizeof(float)));
            snap_lop  = (float *)ARENA_ALLOC(arena, (long)(surf_nv * 3 * sizeof(float)));
            memcpy(snap_orig, surf_verts, surf_nv * 3 * sizeof(float));
            memcpy(snap_lop,  mls_verts,  surf_nv * 3 * sizeof(float));
        }

#ifdef VESUVIUS_DEBUG
        t_c1 = ves_clock_sec();
        fprintf(stderr, "    comp %d: MLS-midpoint (%d iter, R=%.1f vox) "
                "-> %zu verts (%.3f s)\n",
                comp_label, MLS_PROJECT_ITERS,
                (double)MLS_PROJECT_RADIUS_VOX, surf_nv, t_c1 - t_c0);
#endif

        /* Diagnostic dump: LOP-only (post MLS, BEFORE
         * the backface cull). Useful for tracking where holes are
         * introduced — comparing this against the post-cull dump shows
         * exactly which triangles the cull removed. */
        if (dump_cube_dir && cube_id && getenv("EXTRACT_DIAG")) {
            char sub[1024];
            snprintf(sub, sizeof(sub), "%s/%s_step0_post_lop",
                     dump_cube_dir, cube_id);
            DumpObj_ensure_dir(sub);
            char path[1024];
            snprintf(path, sizeof(path),
                     "%s/%s_post_lop_comp%03zu.obj",
                     sub, cube_id, out_count + 1);
            DumpObj_write_one_world(scratch, path, cube_id,
                                    mls_verts, surf_nv, surf_faces, surf_nf,
                                    diag_color);

            /* Also dump the pure oriented point cloud (verts + MLS
             * normals, no faces). Useful for assessing whether the
             * LOP-projected points + normals are clean enough on their
             * own, without committing to any particular triangulation. */
            char psub[1024];
            snprintf(psub, sizeof(psub),
                     "%s/%s_step0_post_lop_points",
                     dump_cube_dir, cube_id);
            DumpObj_ensure_dir(psub);
            char ppath[1024];
            snprintf(ppath, sizeof(ppath),
                     "%s/%s_post_lop_points_comp%03zu.obj",
                     psub, cube_id, out_count + 1);
            DumpObj_write_points_world(scratch, ppath, cube_id,
                                       mls_verts, mls_normals, surf_nv);
        }

        /* Triangulate the LOP'd cloud with BPA to produce surf_faces. The
         * input is a single-layer cloud by construction, so there is no
         * double envelope to cull. The output mesh is manifold by the
         * BPA pre-check; downstream hole_fill closes the small holes
         * BPA leaves where the pivot got stuck on tight curvature.
         *
         * Pre-BPA weld: LOP can collapse multiple voxel-center inputs
         * to the same output position (e.g. a 2-vox-thick FG region's
         * top and bottom layers both projecting to the midline).
         * Without welding, BPA picks both as distinct inputs and emits
         * sliver triangles connecting them — sub-millivoxel edges,
         * aspect ratios > 100. MLS_WELD_EPS_VOX = 0.25 vox merges
         * those duplicates before BPA. */
        {
            size_t welded_nv = 0;
            float *welded_verts = NULL, *welded_normals = NULL;
            /* nf=0 because we have no faces yet; the weld pass becomes
             * a pure vertex deduplication. */
            int32_t dummy_face = 0;
            size_t dummy_out_nf = 0;
            float weld_eps = MLS_WELD_EPS_VOX;
            /* Iterate weld to fixed point — Weld_verts isn't idempotent
             * for our data (cell-boundary effects can leave near-dup
             * pairs after the first pass). Most cubes converge in 2-3
             * passes; cap at 5 as safety. */
            const float *src_v = mls_verts;
            const float *src_n = mls_normals;
            size_t src_nv = surf_nv;
            welded_verts = NULL; welded_normals = NULL; welded_nv = src_nv;
            for (int iter = 0; iter < 5; iter++) {
                size_t prev_nv = welded_nv;
                Weld_verts(scratch,
                           src_v, src_nv, src_n,
                           &dummy_face, 0, &dummy_out_nf,
                           weld_eps,
                           &welded_verts, &welded_nv, &welded_normals);
                if (welded_nv == prev_nv && iter > 0) break;
                src_v = welded_verts;
                src_n = welded_normals;
                src_nv = welded_nv;
            }
#ifdef VESUVIUS_DEBUG
            fprintf(stderr,
                "    comp %d: pre-BPA weld %zu -> %zu verts (eps=%.2f, "
                "iterated)\n",
                comp_label, surf_nv, welded_nv, (double)weld_eps);
#endif
            mls_verts   = welded_verts;
            mls_normals = welded_normals;
            surf_nv       = welded_nv;

            /* Dump the welded cloud BEFORE the owned trim (points only), world
             * coords. Diff against step0_post_lop_trimmed (AFTER trim) to isolate
             * EXACTLY what the owned/halo trim removes, with the weld held constant
             * (so coords match 1:1; no weld-merge confound). */
            if (dump_cube_dir && cube_id && getenv("EXTRACT_DIAG")) {
                char wsub[1024];
                snprintf(wsub, sizeof(wsub), "%s/%s_step0_post_weld",
                         dump_cube_dir, cube_id);
                DumpObj_ensure_dir(wsub);
                char wpath[1024];
                snprintf(wpath, sizeof(wpath),
                         "%s/%s_post_weld_comp%03zu.obj",
                         wsub, cube_id, out_count + 1);
                DumpObj_write_points_world(scratch, wpath, cube_id,
                                           mls_verts, mls_normals, surf_nv);
            }

            /* Pre-BPA owned-region trim (halo mode): drop POINTS within
             * BPA_OWNED_TRIM_INSET vox of every cube face so BPA's open boundary
             * lands strictly INSIDE the face. LOP has already run with full
             * two-sided halo support, so kept owned verts keep their denoised,
             * cross-cube-stable positions; we only remove triangulation input at
             * and beyond the face. Insetting (not splitting exactly at the shared
             * plane) keeps adjacent cubes' charts from TOUCHING/doubling where a
             * wrap grazes the seam; the 2*INSET gap is spanned by the seam bridge.
             * No faces yet -> compact verts + normals in place. */
            if (halo_voxels > 0 && surf_nv > 0) {
                float ins  = (float)BPA_OWNED_TRIM_INSET;
                float lo   = ins;
                float hi_z = (float)cube_D - ins;
                float hi_y = (float)cube_H - ins;
                float hi_x = (float)cube_W - ins;
                size_t kept = 0;
                for (size_t v = 0; v < surf_nv; v++) {
                    float zc = mls_verts[v*3+0];
                    float yc = mls_verts[v*3+1];
                    float xc = mls_verts[v*3+2];
                    if (zc < lo || zc > hi_z ||
                        yc < lo || yc > hi_y ||
                        xc < lo || xc > hi_x) continue;
                    if (kept != v) {
                        mls_verts[kept*3+0] = zc;
                        mls_verts[kept*3+1] = yc;
                        mls_verts[kept*3+2] = xc;
                        mls_normals[kept*3+0] = mls_normals[v*3+0];
                        mls_normals[kept*3+1] = mls_normals[v*3+1];
                        mls_normals[kept*3+2] = mls_normals[v*3+2];
                    }
                    kept++;
                }
#ifdef VESUVIUS_DEBUG
                fprintf(stderr,
                    "    comp %d: pre-BPA owned trim (inset=%.2f) %zu -> %zu verts\n",
                    comp_label, (double)BPA_OWNED_TRIM_INSET, surf_nv, kept);
#endif
                surf_nv = kept;
                /* Dump the LOP cloud AFTER the owned trim (points only, no BPA),
                 * world coords. Compare against step0_post_lop_points (PRE-trim)
                 * to see exactly which points the owned/halo trim removed. */
                if (dump_cube_dir && cube_id && getenv("EXTRACT_DIAG")) {
                    char tsub[1024];
                    snprintf(tsub, sizeof(tsub), "%s/%s_step0_post_lop_trimmed",
                             dump_cube_dir, cube_id);
                    DumpObj_ensure_dir(tsub);
                    char tpath[1024];
                    snprintf(tpath, sizeof(tpath),
                             "%s/%s_post_lop_trimmed_comp%03zu.obj",
                             tsub, cube_id, out_count + 1);
                    DumpObj_write_points_world(scratch, tpath, cube_id,
                                               mls_verts, mls_normals, surf_nv);
                }
                if (surf_nv < 3) continue;   /* nothing left to triangulate */
            }

            /* Hoppe 1992 §3.3 globally-consistent orientation before BPA: replace
             * mls_project's (1,1,1) per-vertex normal sign with a MST-propagated,
             * majority-re-anchored consistent sign, so BPA's normal gate stops
             * stranding pivots at sign-flip bands (the front-split into anti-wound
             * components). DEFAULT ON (validated: fixes the splits, no within-cube
             * regression); set HOPPE_NO_ORIENT to disable. HOPPE_RADIUS overrides
             * the neighbour radius. */
            if (!getenv("HOPPE_NO_ORIENT")) {
                float hr = 1.2f;
                const char *he = getenv("HOPPE_RADIUS");
                if (he) { double t = atof(he); if (t > 0.0) hr = (float)t; }
                size_t hf = 0;
                NormalOrient_consistent(scratch, mls_verts, mls_normals,
                                        surf_nv, hr, &hf);
                fprintf(stderr,
                    "    comp %d: Hoppe orient flipped %zu/%zu normals (r=%.2f)\n",
                    comp_label, hf, surf_nv, (double)hr);
            }

            int32_t *bpa_faces = NULL;
            size_t   bpa_nf = 0;
            int rc = BallPivot_reconstruct(scratch,
                                            mls_verts, mls_normals, surf_nv,
                                            BPA_RHO_VOX,
                                            &bpa_faces, &bpa_nf);
            if (rc != 0 || bpa_nf == 0) {
                fprintf(stderr,
                    "    comp %d: BPA failed (rc=%d, nf=%zu) — skip\n",
                    comp_label, rc, bpa_nf);
                continue;
            }
            surf_faces = bpa_faces;
            surf_nf    = bpa_nf;
            surf_verts = mls_verts;
#ifdef VESUVIUS_DEBUG
            t_c1 = ves_clock_sec();
            fprintf(stderr,
                "    comp %d: BPA rho=%.1f -> %zu faces (%zu verts, %.2f s)\n",
                comp_label, (double)BPA_RHO_VOX, surf_nf, surf_nv, t_c1 - t_c0);
#endif
        }

        /* Sign-propagation BFS — DISABLED by default after the 2026-05-27
         * test run showed it backfires on thick (>R_kernel) connected
         * components: the kernel can't span the thickness so front- and
         * back-side vertices have genuinely opposite local outward
         * directions, and the BFS forces them to agree, making the cull
         * keep both sides. Net effect on the 2-cube test was +50k extra
         * faces and slightly MORE open-boundary edges (4846 vs 4443),
         * with no measurable tear reduction. Set SIGN_PROPAGATE=1 to
         * re-enable for further experimentation. See changelog. */
        size_t n_sign_flips = 0;
        if (getenv("SIGN_PROPAGATE") != NULL) {
            n_sign_flips = propagate_normal_signs(scratch,
                                                  mls_verts, surf_nv,
                                                  surf_faces, surf_nf,
                                                  mls_normals);
            fprintf(stderr,
                    "    comp %d: sign-propagate flipped %zu/%zu vert "
                    "normals\n",
                    comp_label, n_sign_flips, surf_nv);
        }
        (void)n_sign_flips;

        /* BPA assigns each triangle's winding from its seed vertex's MLS
         * normal sign, which carries the (1,1,1)-hemisphere bias of
         * mls_project.c and therefore flips across the scroll wrap — ~6%
         * of interior edges come out wound the SAME direction by their two
         * faces (inconsistent winding). Run a deterministic per-component
         * BFS that recomputes each face's shared-edge direction from the
         * LIVE (already-flipped) face arrays at visit time and flips any
         * neighbour winding a shared edge the same way as its parent, then
         * anchors each component's global sign against the MLS normal field.
         * This drives same-direction interior edges to ~0 (residuals are
         * odd-cycle pinch verts, cleared by the pinhole pass). Must run
         * per-cube here, before export: grid_weld's own winding repair has
         * a stale-edge bug and cannot converge. */
        {
            size_t or_flips = 0, or_comps = 0, or_resid = 0;
            OrientMesh_consistent(scratch,
                                  mls_verts, surf_nv, mls_normals,
                                  surf_faces, surf_nf,
                                  &or_flips, &or_comps, &or_resid);
#ifdef VESUVIUS_DEBUG
            fprintf(stderr,
                "    comp %d: orient BFS %zu comps, %zu flips, "
                "%zu residual same-dir edges\n",
                comp_label, or_comps, or_flips, or_resid);
#endif
        }

        surf_verts = mls_verts;

        (void)Weld_verts;             /* kept available for future weld pass */
        (void)spatial_pairing_cull;   /* attempted 2026-05-27, see changelog */

        if (surf_nf == 0 || surf_nv == 0) continue;

        /* Mesh cleanup: UF split + remove small fragments */
        surf_nf = mesh_split_remove_small(scratch, surf_nv,
                                        surf_faces, surf_nf);
        if (surf_nf == 0) continue;

        /* Clean unreferenced vertices */
        float *clean_verts = NULL;
        size_t clean_nv = 0;
        clean_unreferenced(scratch, surf_verts, surf_nv, surf_faces, surf_nf,
                           &clean_verts, &clean_nv);
        surf_verts = clean_verts;
        surf_nv = clean_nv;

#ifdef VESUVIUS_DEBUG
        t_c1 = ves_clock_sec();
        fprintf(stderr, "    comp %d: cleanup -> %zu verts, %zu faces (%.3f s)\n",
                comp_label, surf_nv, surf_nf,
                t_c1 - t_c0);
#endif

        /* No micro-hole fill here: BPA already handles small holes via
         * re-seeding, and the centroid-insertion fill_micro_holes uses for
         * 5-6 vert loops lands ~0.001 vox from existing verts on near-flat
         * sheets, producing slivers (aspect ratio 100+). The boundary loops
         * BPA leaves are legitimate sheet edges, not micro-holes to patch. */

#ifdef VESUVIUS_DEBUG
        t_c1 = ves_clock_sec();
        fprintf(stderr, "    comp %d: post-cleanup -> %zu verts, %zu faces (%.3f s)\n",
                comp_label, surf_nv, surf_nf,
                t_c1 - t_c0);
#endif

        if (surf_nf == 0 || surf_nv == 0) continue;

        /* Phase B: tighten the EXPORT halo to EXPORT_HALO_VOX (typically 1).
         * The wide read halo (halo_voxels) was needed for MC determinism
         * and for LOP's R=12 kernel to have enough neighbour data; but the
         * exported per-cube OBJ only needs ~1 vox of overlap with each
         * neighbour so grid_weld's 1e-4 hash can find duplicates. The
         * intermediate halo verts are noise that makes downstream OBJs
         * huge and visually hard to read.
         *
         * Same drop pattern as Phase A above, just at a tighter threshold. */
        if (halo_voxels > EXPORT_HALO_VOX && surf_nv > 0) {
            float export_lo = -(float)EXPORT_HALO_VOX;
            float export_hi_z = (float)cube_D + (float)EXPORT_HALO_VOX;
            float export_hi_y = (float)cube_H + (float)EXPORT_HALO_VOX;
            float export_hi_x = (float)cube_W + (float)EXPORT_HALO_VOX;
            int32_t *remap_e = (int32_t *)ARENA_ALLOC(scratch,
                                (long)surf_nv * (long)sizeof(int32_t));
            size_t nv_kept_e = 0;
            for (size_t v = 0; v < surf_nv; v++) {
                float zc = surf_verts[v * 3 + 0];
                float yc = surf_verts[v * 3 + 1];
                float xc = surf_verts[v * 3 + 2];
                if (zc < export_lo || zc > export_hi_z ||
                    yc < export_lo || yc > export_hi_y ||
                    xc < export_lo || xc > export_hi_x) {
                    remap_e[v] = -1;
                } else {
                    remap_e[v] = (int32_t)nv_kept_e;
                    if (nv_kept_e != v) {
                        surf_verts[nv_kept_e * 3 + 0] = zc;
                        surf_verts[nv_kept_e * 3 + 1] = yc;
                        surf_verts[nv_kept_e * 3 + 2] = xc;
                    }
                    nv_kept_e++;
                }
            }
            size_t nf_kept_e = 0;
            for (size_t f = 0; f < surf_nf; f++) {
                int32_t v0 = surf_faces[f * 3 + 0];
                int32_t v1 = surf_faces[f * 3 + 1];
                int32_t v2 = surf_faces[f * 3 + 2];
                if (remap_e[v0] < 0 || remap_e[v1] < 0 || remap_e[v2] < 0) continue;
                surf_faces[nf_kept_e * 3 + 0] = remap_e[v0];
                surf_faces[nf_kept_e * 3 + 1] = remap_e[v1];
                surf_faces[nf_kept_e * 3 + 2] = remap_e[v2];
                nf_kept_e++;
            }
#ifdef VESUVIUS_DEBUG
            fprintf(stderr,
                "    comp %d: export-halo trim (read=%d, export=%d) "
                "%zu->%zu v, %zu->%zu f\n",
                comp_label, halo_voxels, EXPORT_HALO_VOX,
                surf_nv, nv_kept_e, surf_nf, nf_kept_e);
#endif
            surf_nv = nv_kept_e;
            surf_nf = nf_kept_e;
            if (surf_nf == 0 || surf_nv == 0) continue;
        }

        /* Dump BPA surface output (stage name still "step0_pre_simplify";
         * QEM has been moved
         * to post-Step-6 in main.c so there is no in-Step-0 QEM call
         * anymore). The `skip_qem` parameter is now ignored here. */
        (void)skip_qem;
        if (dump_cube_dir && cube_id && getenv("EXTRACT_DIAG")) {
            char sub[1024];
            snprintf(sub, sizeof(sub), "%s/%s_step0_pre_simplify",
                     dump_cube_dir, cube_id);
            DumpObj_ensure_dir(sub);
            char path[1024];
            snprintf(path, sizeof(path), "%s/%s_pre_simplify_comp%03zu.obj",
                     sub, cube_id, out_count + 1);
            DumpObj_write_one_world(scratch, path, cube_id,
                                    surf_verts, surf_nv, surf_faces, surf_nf,
                                    diag_color);
        }

        /* Copy final results to main arena */
        float *final_verts = (float *)ARENA_ALLOC(arena,
                                                    (long)(surf_nv * 3 * sizeof(float)));
        memcpy(final_verts, surf_verts, surf_nv * 3 * sizeof(float));

        int32_t *final_faces = (int32_t *)ARENA_ALLOC(arena,
                                                        (long)(surf_nf * 3 * sizeof(int32_t)));
        memcpy(final_faces, surf_faces, surf_nf * 3 * sizeof(int32_t));

        /* TEMP (matches the weld-skipping diagnostic above): we skip the
         * second MLS pass for normals — those are only used by welding
         * and downstream, neither of which runs in MC+LOP-only mode.
         * Saves ~50s per cube in Debug. */
        float *final_normals = NULL;

        /* Per-mesh PCA normal/centroid. The mesh is now single-sided
         * (MLS + weld collapsed the double envelope), so this PCA is
         * unambiguous — both eigenvectors-of-this-mesh and the welded
         * vert_normals agree in direction, and the (1,1,1)-orientation
         * convention in PCA_normal() gives the same sign as MLS does. */
        float pca_normal[3] = {0, 0, 1};
        float centroid[3] = {0, 0, 0};
        PCA_normal(final_verts, surf_nv, pca_normal, centroid);

        /* Fill ComponentMesh */
        ComponentMesh *cm = &meshes[out_count];
        cm->verts = final_verts;
        cm->faces = final_faces;
        cm->vert_normals = final_normals;
        cm->nv = surf_nv;
        cm->nf = surf_nf;
        cm->comp_id = (int)(out_count + 1);
        cm->pca_normal[0] = pca_normal[0];
        cm->pca_normal[1] = pca_normal[1];
        cm->pca_normal[2] = pca_normal[2];
        cm->centroid[0] = centroid[0];
        cm->centroid[1] = centroid[1];
        cm->centroid[2] = centroid[2];
        cm->self = cm;

        /* No pin classification: the cross-cube weld is a BPA pass over the
         * merged sheets (SeamWeld_bridge), not a bit-exact seam join, so there
         * is no pin band to compute. The raw-MC pin-snapback above still keeps
         * seam verts deterministic for the bridge; downstream stages no longer
         * read a pin_mask, so leave it NULL. */
        cm->pin_mask = NULL;

        if (clouds) {
            clouds[out_count].orig_pts = snap_orig;
            clouds[out_count].lop_pts  = snap_lop;
            clouds[out_count].cell_origin[0] = cube_world_origin[0];
            clouds[out_count].cell_origin[1] = cube_world_origin[1];
            clouds[out_count].cell_origin[2] = cube_world_origin[2];
            clouds[out_count].n = snap_n;
        }

        assert(ComponentMesh_valid(cm));
        out_count++;

#ifdef VESUVIUS_DEBUG
        fprintf(stderr, "  Extract: comp %d -> %zu verts, %zu faces, "
                "normal=(%.3f,%.3f,%.3f)\n",
                cm->comp_id, cm->nv, cm->nf,
                (double)cm->pca_normal[0], (double)cm->pca_normal[1],
                (double)cm->pca_normal[2]);
#endif
    }

    Arena_dispose(&scratch);

#ifdef VESUVIUS_DEBUG
    {
        double t_end = ves_clock_sec();
        fprintf(stderr, "  Extract: %zu components, total %.3f s\n", out_count,
                t_end - t_start);
    }
#endif

    *out_meshes = meshes;
    *out_n_meshes = out_count;
    if (out_clouds) *out_clouds = clouds;
    return 0;
}
