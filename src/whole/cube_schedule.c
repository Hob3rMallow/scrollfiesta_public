/* cube_schedule.c -- dump-grid discovery + spiral registration order.
 * See cube_schedule.h for the design. */
#include "../common/ves_platform.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cube_schedule.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- id parsing ----------------------------------------------------------- */

static int parse_cube_origin(const char *cube_id,
                             int64_t *vz, int64_t *vy, int64_t *vx)
{
    int iz = 0, iy = 0, ix = 0;
    if (sscanf(cube_id, "z%d_y%d_x%d", &iz, &iy, &ix) != 3) return -1;
    *vz = iz;
    *vy = iy;
    *vx = ix;
    return 0;
}

/* ---- directory scan (grid_weld.c pattern) --------------------------------- */

typedef struct {
    CubeNode *nodes;
    size_t    n, cap;
} NodeList;

static void nodelist_push(Arena_T arena, NodeList *nl, const CubeNode *nd)
{
    if (nl->n >= nl->cap) {
        size_t new_cap = nl->cap == 0 ? 64 : nl->cap * 2;
        CubeNode *nn = (CubeNode *)ARENA_ALLOC(arena, new_cap * sizeof(CubeNode));
        if (nl->nodes != NULL) memcpy(nn, nl->nodes, nl->n * sizeof(CubeNode));
        nl->nodes = nn;
        nl->cap = new_cap;
    }
    nl->nodes[nl->n++] = *nd;
}

static int leaf_obj_exists(const char *dump_dir, const char *id,
                           const char *leaf_stage)
{
    char path[1024];
    snprintf(path, sizeof(path), "%s/%s/%s_%s/%s_%s_all.obj",
             dump_dir, id, id, leaf_stage, id, leaf_stage);
    FILE *f = fopen(path, "rb");
    if (f == NULL) return 0;
    fclose(f);
    return 1;
}

#ifdef _MSC_VER
#include <windows.h>

static int scan_cubes(Arena_T arena, const char *dump_dir,
                      const char *leaf_stage, NodeList *out)
{
    char pattern[1024];
    snprintf(pattern, sizeof(pattern), "%s/*", dump_dir);
    WIN32_FIND_DATAA fd;
    HANDLE h = FindFirstFileA(pattern, &fd);
    if (h == INVALID_HANDLE_VALUE) return -1;
    do {
        if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) continue;
        if (fd.cFileName[0] != 'z') continue;
        CubeNode nd;
        memset(&nd, 0, sizeof(nd));
        if (parse_cube_origin(fd.cFileName, &nd.oz, &nd.oy, &nd.ox) != 0) continue;
        if (!leaf_obj_exists(dump_dir, fd.cFileName, leaf_stage)) continue;
        snprintf(nd.id, sizeof(nd.id), "%s", fd.cFileName);
        nodelist_push(arena, out, &nd);
    } while (FindNextFileA(h, &fd));
    FindClose(h);
    return 0;
}
#else
#include <dirent.h>

static int scan_cubes(Arena_T arena, const char *dump_dir,
                      const char *leaf_stage, NodeList *out)
{
    DIR *d = opendir(dump_dir);
    if (d == NULL) return -1;
    struct dirent *ent = NULL;
    while ((ent = readdir(d)) != NULL) {
        if (ent->d_name[0] != 'z') continue;
        CubeNode nd;
        memset(&nd, 0, sizeof(nd));
        if (parse_cube_origin(ent->d_name, &nd.oz, &nd.oy, &nd.ox) != 0) continue;
        if (!leaf_obj_exists(dump_dir, ent->d_name, leaf_stage)) continue;
        snprintf(nd.id, sizeof(nd.id), "%s", ent->d_name);
        nodelist_push(arena, out, &nd);
    }
    closedir(d);
    return 0;
}
#endif

/* ---- neighbor linking (sorted origins + binary search) -------------------- */

typedef struct { int64_t oz, oy, ox; int32_t idx; } OriginKey;

static int cmp_origin(const void *pa, const void *pb)
{
    const OriginKey *a = (const OriginKey *)pa, *b = (const OriginKey *)pb;
    if (a->oz != b->oz) return a->oz < b->oz ? -1 : 1;
    if (a->oy != b->oy) return a->oy < b->oy ? -1 : 1;
    if (a->ox != b->ox) return a->ox < b->ox ? -1 : 1;
    return 0;
}

static int32_t origin_find(const OriginKey *keys, size_t n,
                           int64_t oz, int64_t oy, int64_t ox)
{
    size_t lo = 0, hi = n;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        const OriginKey *k = &keys[mid];
        int c = 0;
        if (k->oz != oz)      c = k->oz < oz ? -1 : 1;
        else if (k->oy != oy) c = k->oy < oy ? -1 : 1;
        else if (k->ox != ox) c = k->ox < ox ? -1 : 1;
        if (c == 0) return k->idx;
        if (c < 0) lo = mid + 1; else hi = mid;
    }
    return -1;
}

static int cmp_dbl_cs(const void *pa, const void *pb)
{
    double a = *(const double *)pa, b = *(const double *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

int CubeSched_link_and_order(Arena_T arena, CubeNode *nodes, size_t n,
                             int64_t chunk, const float axis_point[3],
                             double pitch, int32_t seed_idx,
                             int32_t **out_order, int32_t **out_comp)
{
    assert(nodes && out_order && out_comp);
    if (n == 0) return -1;
    const int radial_order = pitch <= 0.0;

    /* polar geometry at cube centers (axis = +z through axis_point in (z,y,x)) */
    for (size_t i = 0; i < n; i++) {
        double cy = (double)nodes[i].oy + (double)chunk * 0.5;
        double cx = (double)nodes[i].ox + (double)chunk * 0.5;
        double dy = cy - (double)axis_point[1];
        double dx = cx - (double)axis_point[2];
        nodes[i].r       = sqrt(dy * dy + dx * dx);
        nodes[i].theta   = atan2(dx, dy);
        nodes[i].w_phase = radial_order
                         ? nodes[i].r
                         : nodes[i].r / pitch - nodes[i].theta / (2.0 * M_PI);
    }

    /* results outlive the scratch mark below */
    int32_t *order = (int32_t *)ARENA_ALLOC(arena, n * sizeof(int32_t));
    int32_t *comp  = (int32_t *)ARENA_ALLOC(arena, n * sizeof(int32_t));

    /* 6-neighbor links */
    Arena_Mark mark = Arena_save(arena);
    OriginKey *keys = (OriginKey *)ARENA_ALLOC(arena, n * sizeof(OriginKey));
    for (size_t i = 0; i < n; i++) {
        keys[i].oz = nodes[i].oz; keys[i].oy = nodes[i].oy;
        keys[i].ox = nodes[i].ox; keys[i].idx = (int32_t)i;
    }
    qsort(keys, n, sizeof(OriginKey), cmp_origin);
    for (size_t i = 0; i < n; i++) {
        const int64_t oz = nodes[i].oz, oy = nodes[i].oy, ox = nodes[i].ox;
        nodes[i].nbr[0] = origin_find(keys, n, oz - chunk, oy, ox);
        nodes[i].nbr[1] = origin_find(keys, n, oz + chunk, oy, ox);
        nodes[i].nbr[2] = origin_find(keys, n, oz, oy - chunk, ox);
        nodes[i].nbr[3] = origin_find(keys, n, oz, oy + chunk, ox);
        nodes[i].nbr[4] = origin_find(keys, n, oz, oy, ox - chunk);
        nodes[i].nbr[5] = origin_find(keys, n, oz, oy, ox + chunk);
    }

    /* seed: median center radius (a clean mid-shell band; the core is the
     * LEAST reliable place to anchor a registration chain) */
    if (seed_idx < 0) {
        double *rs = (double *)ARENA_ALLOC(arena, n * sizeof(double));
        for (size_t i = 0; i < n; i++) rs[i] = nodes[i].r;
        qsort(rs, n, sizeof(double), cmp_dbl_cs);
        double rmed = rs[n / 2];
        double best = 1e300;
        for (size_t i = 0; i < n; i++) {
            double d = fabs(nodes[i].r - rmed);
            if (d < best) { best = d; seed_idx = (int32_t)i; }
        }
    }
    assert(seed_idx >= 0 && (size_t)seed_idx < n);
    double wp0 = nodes[seed_idx].w_phase;

    /* best-first flood: min |w_phase - wp0| among unplaced cubes adjacent to a
     * placed one; re-seed (new flood component) when the frontier empties.
     * O(n^2) selection -- n is a few thousand cubes, microseconds. */
    uint8_t *placed = (uint8_t *)ARENA_CALLOC(arena, n, 1);
    uint8_t *elig   = (uint8_t *)ARENA_CALLOC(arena, n, 1);
    int32_t ncomp = 0;
    elig[seed_idx] = 1;
    for (size_t t = 0; t < n; t++) {
        int32_t pick = -1;
        double bestk = 1e300;
        for (size_t i = 0; i < n; i++) {
            if (placed[i] || !elig[i]) continue;
            double k = fabs(nodes[i].w_phase - wp0);
            if (k < bestk) { bestk = k; pick = (int32_t)i; }
        }
        if (pick < 0) {
            /* disconnected component: re-seed at the min-key unplaced cube */
            for (size_t i = 0; i < n; i++) {
                if (placed[i]) continue;
                double k = fabs(nodes[i].w_phase - wp0);
                if (k < bestk) { bestk = k; pick = (int32_t)i; }
            }
            ncomp++;
        }
        assert(pick >= 0);
        placed[pick] = 1;
        order[t] = pick;
        comp[pick] = ncomp;
        for (int e = 0; e < 6; e++) {
            int32_t nb = nodes[pick].nbr[e];
            if (nb >= 0 && !placed[nb]) elig[nb] = 1;
        }
    }

    Arena_restore(arena, mark);
    *out_order = order;
    *out_comp = comp;
    return 0;
}

int CubeSched_build(Arena_T arena, const char *dump_dir, const char *leaf_stage,
                    int64_t chunk, const float axis_point[3], double pitch,
                    const char *seed_id,
                    CubeNode **out_nodes, size_t *out_n,
                    int32_t **out_order, int32_t **out_comp)
{
    assert(dump_dir && leaf_stage && out_nodes && out_n && out_order && out_comp);
    NodeList nl;
    memset(&nl, 0, sizeof(nl));
    if (scan_cubes(arena, dump_dir, leaf_stage, &nl) != 0) return -1;
    if (nl.n == 0) return -1;

    int32_t seed_idx = -1;
    if (seed_id != NULL && seed_id[0] != '\0') {
        for (size_t i = 0; i < nl.n; i++)
            if (strcmp(nl.nodes[i].id, seed_id) == 0) { seed_idx = (int32_t)i; break; }
        if (seed_idx < 0) {
            fprintf(stderr, "cube_schedule: WARN seed id %s not in dump; "
                    "using median-radius seed\n", seed_id);
        }
    }
    int rc = CubeSched_link_and_order(arena, nl.nodes, nl.n, chunk, axis_point,
                                      pitch, seed_idx, out_order, out_comp);
    if (rc != 0) return rc;
    *out_nodes = nl.nodes;
    *out_n = nl.n;
    return 0;
}

/* ============================================================================
 * Self-test.
 * ==========================================================================*/

static int cs_check(int cond, const char *what, int *fails)
{
    if (!cond) {
        fprintf(stderr, "[cube_schedule selftest]   FAIL: %s\n", what);
        (*fails)++;
    }
    return cond;
}

/* build an ny*nx single-layer grid of nodes; skip (sy,sx) when skip_one */
static size_t cs_make_grid(CubeNode *nodes, int ny, int nx, int64_t chunk,
                           int skip_one, int sy, int sx)
{
    size_t n = 0;
    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
            if (skip_one && y == sy && x == sx) continue;
            CubeNode *nd = &nodes[n++];
            memset(nd, 0, sizeof(*nd));
            nd->oz = 0;
            nd->oy = (int64_t)y * chunk;
            nd->ox = (int64_t)x * chunk;
            snprintf(nd->id, sizeof(nd->id), "z%05d_y%05d_x%05d",
                     0, (int)nd->oy, (int)nd->ox);
        }
    }
    return n;
}

/* invariant: every cube in the same flood component as an earlier cube must be
 * 6-adjacent to some already-placed cube at its turn */
static int cs_order_valid(const CubeNode *nodes, size_t n,
                          const int32_t *order, const int32_t *comp)
{
    uint8_t placed[4096];
    if (n > sizeof(placed)) return 0;
    memset(placed, 0, n);
    for (size_t t = 0; t < n; t++) {
        int32_t i = order[t];
        if (i < 0 || (size_t)i >= n || placed[i]) return 0;
        int has_nbr = 0, first_of_comp = 1;
        for (size_t s = 0; s < t; s++)
            if (comp[order[s]] == comp[i]) { first_of_comp = 0; break; }
        for (int e = 0; e < 6; e++) {
            int32_t nb = nodes[i].nbr[e];
            if (nb >= 0 && placed[nb]) has_nbr = 1;
        }
        if (!first_of_comp && !has_nbr) return 0;
        placed[i] = 1;
    }
    return 1;
}

int CubeSched_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    const int64_t chunk = 128;
    const float axis[3] = { 0.0f, 320.0f, 320.0f };   /* center of a 5x5 grid */

    /* (0) pitch-free scheduling is genuinely pitch-free: radial distance is
     * the provisional phase, not a hidden dataset-specific pitch fallback. */
    {
        CubeNode nodes[4];
        size_t n = cs_make_grid(nodes, 2, 2, chunk, 0, 0, 0);
        int32_t *order = NULL, *comp = NULL;
        int rc = CubeSched_link_and_order(arena, nodes, n, chunk, axis, 0.0,
                                          -1, &order, &comp);
        int radial = rc == 0;
        for (size_t i = 0; i < n && radial; i++)
            if (fabs(nodes[i].w_phase - nodes[i].r) > 1e-12) radial = 0;
        cs_check(rc == 0 && radial, "auto order has no hidden pitch", &fails);
    }

    /* (1) full 5x5 grid: all placed, one component, adjacency invariant,
     * seed is the median-radius cube, core (min-r) is placed LATE */
    {
        CubeNode nodes[25];
        size_t n = cs_make_grid(nodes, 5, 5, chunk, 0, 0, 0);
        int32_t *order = NULL, *comp = NULL;
        int rc = CubeSched_link_and_order(arena, nodes, n, chunk, axis, 10.0,
                                          -1, &order, &comp);
        cs_check(rc == 0, "grid rc", &fails);
        cs_check(cs_order_valid(nodes, n, order, comp), "grid order invariant", &fails);
        int allc0 = 1;
        for (size_t i = 0; i < n; i++) if (comp[i] != 0) allc0 = 0;
        cs_check(allc0, "grid single flood component", &fails);
        /* the min-radius (core) cube must not be among the first placed */
        size_t core = 0;
        for (size_t i = 1; i < n; i++) if (nodes[i].r < nodes[core].r) core = i;
        size_t core_rank = 0;
        for (size_t t = 0; t < n; t++) if ((size_t)order[t] == core) core_rank = t;
        fprintf(stderr, "[cube_schedule selftest] (1) 5x5: core rank %zu/%zu\n",
                core_rank, n);
        cs_check(core_rank > n / 2, "core registered late", &fails);
        /* neighbor sanity: interior cube has all 4 in-plane links, no z links */
        size_t mid = 0;
        for (size_t i = 0; i < n; i++)
            if (nodes[i].oy == 2 * chunk && nodes[i].ox == 2 * chunk) mid = i;
        cs_check(nodes[mid].nbr[0] == -1 && nodes[mid].nbr[1] == -1,
                 "no z neighbors in single layer", &fails);
        cs_check(nodes[mid].nbr[2] >= 0 && nodes[mid].nbr[3] >= 0 &&
                 nodes[mid].nbr[4] >= 0 && nodes[mid].nbr[5] >= 0,
                 "interior 4-linked", &fails);
    }

    /* (2) grid with a hole still orders everything */
    {
        CubeNode nodes[25];
        size_t n = cs_make_grid(nodes, 5, 5, chunk, 1, 2, 1);
        int32_t *order = NULL, *comp = NULL;
        int rc = CubeSched_link_and_order(arena, nodes, n, chunk, axis, 10.0,
                                          -1, &order, &comp);
        cs_check(rc == 0 && cs_order_valid(nodes, n, order, comp),
                 "hole grid ordered", &fails);
    }

    /* (3) detached island gets its own flood component */
    {
        CubeNode nodes[26];
        size_t n = cs_make_grid(nodes, 5, 5, chunk, 0, 0, 0);
        CubeNode *isl = &nodes[n++];
        memset(isl, 0, sizeof(*isl));
        isl->oz = 0; isl->oy = 20 * chunk; isl->ox = 20 * chunk;
        snprintf(isl->id, sizeof(isl->id), "z%05d_y%05d_x%05d",
                 0, (int)isl->oy, (int)isl->ox);
        int32_t *order = NULL, *comp = NULL;
        int rc = CubeSched_link_and_order(arena, nodes, n, chunk, axis, 10.0,
                                          -1, &order, &comp);
        cs_check(rc == 0 && cs_order_valid(nodes, n, order, comp),
                 "island grid ordered", &fails);
        cs_check(comp[n - 1] == 1, "island in flood component 1", &fails);
        cs_check((size_t)order[n - 1] == n - 1, "island placed last", &fails);
    }

    Arena_dispose(&arena);
    fprintf(stderr, "[cube_schedule selftest] %s (%d failure%s)\n",
            fails == 0 ? "PASSED" : "FAILED", fails, fails == 1 ? "" : "s");
    return fails;
}
