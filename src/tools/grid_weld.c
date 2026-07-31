/*
 * grid_weld.c -- Weld per-cube OBJ meshes into a single grid-spanning OBJ.
 *
 * The cubes are PLAIN-CONCATENATED (verts are already in source-space world
 * coords; no cross-cube vertex fusing) and then welded by re-running BPA
 * across each detected cube-boundary plane (SeamWeld_bridge). There is no
 * bitwise hash-join: the BPA mesh path does not pin seam verts, and step0
 * trims each cube to its owned box, so adjacent cubes never share coincident
 * seam verts -- the BPA bridge is what fuses them. After the bridge the tool
 * runs winding repair, tiny-component cull, fold cleanup, a pinhole reclose,
 * and a manifold audit, then writes the welded OBJ + JSON report.
 *
 * Usage:
 *   grid_weld <grid_obj_dir> <output.obj>
 *
 * grid_obj_dir layout (matches DumpObj_write_meshes):
 *   <grid_obj_dir>/<cube_id>/<cube_id>_step12_final/<cube_id>_step12_final_all.obj
 *
 * Output:
 *   <output.obj>                  -- welded single OBJ
 *   <output.obj>.weld_report.json -- vert/face/edge counts + manifold audit
 */
#include "../common/ves_platform.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846  /* MSVC math.h omits it without _USE_MATH_DEFINES */
#endif

#include "../common/arena.h"
#include "../common/except.h"
#include "../common/obj_io.h"
#include "../common/cc_color.h"
#include "../common/mesh_types.h"
#include "../common/mesh_manifold.h"
#include "../common/pipeline_constants.h"
#include "../remesh/seam_weld.h"
#include "../remesh/seam_planes.h"
#include "../remesh/seam_refine.h"
#include "../remesh/seam_band_cvt.h"
#include "../remesh/sheet_reweld.h"
#include "../remesh/seam_hole_fill.h"
#include "../remesh/fold_cleanup.h"
#include "../remesh/pinhole_fill.h"
#include "../remesh/manifold_guard.h"
#include "../remesh/orient_mesh.h"
#include "../remesh/orient_weld.h"
#include "../remesh/weld_cleanup.h"
#include "../holefill/hole_fill.h"
#include "../common/vert_weld.h"
#include "../flatten/developability.h"

#define CUBE_SIZE_VOX 128.0f    /* Cube boundary planes are integer multiples */
#define SEAM_AXIS_EPS 0.6f      /* vert is "on" a cube-boundary plane if any
                                 * coord is within this distance of an
                                 * integer multiple of CUBE_SIZE_VOX */

/* 1 if any coordinate is within `zone` vox of a cube-boundary plane (a multiple
 * of CUBE_SIZE_VOX) -- i.e. the vertex sits in the cross-cube seam/weld zone, as
 * opposed to a per-cube interior. Used by the lambda gate to target welds. */
static int near_cube_boundary(const float *v, float zone)
{
    for (int k = 0; k < 3; k++) {
        float m = v[k] - CUBE_SIZE_VOX * floorf(v[k] / CUBE_SIZE_VOX);
        if (m < zone || m > CUBE_SIZE_VOX - zone) return 1;
    }
    return 0;
}

/* "Nice" palette: 16 saturated but distinguishable hues (avoids mud and
 * pure green which is reserved for seam lines). Hand-picked rather than
 * generated so cubes have stable, visually distinct colors across runs. */
static const float CUBE_PALETTE[16][3] = {
    { 0.894f, 0.102f, 0.110f },  /* red       */
    { 0.215f, 0.494f, 0.722f },  /* blue      */
    { 1.000f, 0.498f, 0.000f },  /* orange    */
    { 0.596f, 0.306f, 0.639f },  /* purple    */
    { 0.651f, 0.337f, 0.157f },  /* brown     */
    { 0.969f, 0.506f, 0.749f },  /* pink      */
    { 0.301f, 0.686f, 0.290f },  /* leaf green (distinguishable from seam) */
    { 0.498f, 0.498f, 0.498f },  /* gray      */
    { 0.984f, 0.890f, 0.435f },  /* yellow    */
    { 0.400f, 0.761f, 0.647f },  /* teal      */
    { 0.988f, 0.553f, 0.384f },  /* coral     */
    { 0.553f, 0.627f, 0.796f },  /* lavender  */
    { 0.906f, 0.541f, 0.765f },  /* magenta   */
    { 0.706f, 0.871f, 0.412f },  /* lime-ish  */
    { 0.553f, 0.310f, 0.310f },  /* maroon    */
    { 0.165f, 0.631f, 0.596f }   /* sea green */
};
static const float SEAM_COLOR[3] = { 0.000f, 1.000f, 0.000f };  /* #00FF00 bridge weld verts */
static const float PINHOLE_COLOR[3] = { 0.100f, 0.450f, 1.000f }; /* bright blue: pinhole/holefill-added verts (red carries a "bad" connotation) */

/* ===================================================================
 * Cube ID parsing -- "z{vz:05d}_y{vy:05d}_x{vx:05d}" -> origin offsets.
 * =================================================================== */

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

/* Existence check for the node-OBJ resolver's nested-vs-flat fallback. */
static int gw_file_exists(const char *p)
{
    FILE *f = fopen(p, "rb");
    if (f) { fclose(f); return 1; }
    return 0;
}

/* ===================================================================
 * Directory enumeration -- list <cube_id> subdirectories of grid_obj_dir.
 * =================================================================== */

typedef struct {
    char **ids;
    size_t n;
    size_t cap;
} CubeList;

static void cubelist_push(Arena_T arena, CubeList *cl, const char *id)
{
    if (cl->n >= cl->cap) {
        size_t new_cap = cl->cap == 0 ? 16 : cl->cap * 2;
        char **new_ids = (char **)ARENA_ALLOC(arena,
                            (long)(new_cap * sizeof(char *)));
        if (cl->ids) {
            memcpy(new_ids, cl->ids, cl->n * sizeof(char *));
        }
        cl->ids = new_ids;
        cl->cap = new_cap;
    }
    size_t len = strlen(id);
    char *copy = (char *)ARENA_ALLOC(arena, (long)(len + 1));
    memcpy(copy, id, len + 1);
    cl->ids[cl->n++] = copy;
}

#ifdef _MSC_VER
#include <windows.h>

static int enumerate_cube_dirs(Arena_T arena, const char *base_dir,
                               CubeList *out)
{
    char pattern[1024];
    snprintf(pattern, sizeof(pattern), "%s/*", base_dir);
    WIN32_FIND_DATAA fd;
    HANDLE h = FindFirstFileA(pattern, &fd);
    if (h == INVALID_HANDLE_VALUE) return -1;
    do {
        if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) continue;
        if (strcmp(fd.cFileName, ".") == 0) continue;
        if (strcmp(fd.cFileName, "..") == 0) continue;
        /* Heuristic: cube IDs match "z*_y*_x*" pattern (start with z). */
        if (fd.cFileName[0] != 'z') continue;
        cubelist_push(arena, out, fd.cFileName);
    } while (FindNextFileA(h, &fd));
    FindClose(h);
    return 0;
}
#else
#include <dirent.h>

static int enumerate_cube_dirs(Arena_T arena, const char *base_dir,
                               CubeList *out)
{
    DIR *d = opendir(base_dir);
    if (!d) return -1;
    struct dirent *ent;
    while ((ent = readdir(d)) != NULL) {
        if (ent->d_name[0] != 'z') continue;
        if (strchr(ent->d_name, '_') == NULL) continue;
        cubelist_push(arena, out, ent->d_name);
    }
    closedir(d);
    return 0;
}
#endif

/* Keep only cubes whose origin falls inside an inclusive bbox (--subgrid).
 * Lets grid_weld stitch one rectangular sub-block of a larger dump dir without
 * symlink/junction forests (whose recursive deletion could clobber the real
 * per-cube dumps). Compacts the id array in place. */
static void cubelist_filter_bbox(CubeList *cl,
                                 int64_t z0, int64_t z1, int64_t y0, int64_t y1,
                                 int64_t x0, int64_t x1)
{
    size_t w = 0;
    for (size_t i = 0; i < cl->n; i++) {
        int64_t vz = 0, vy = 0, vx = 0;
        if (parse_cube_origin(cl->ids[i], &vz, &vy, &vx) != 0) continue;
        if (vz >= z0 && vz <= z1 && vy >= y0 && vy <= y1 &&
            vx >= x0 && vx <= x1) {
            cl->ids[w++] = cl->ids[i];
        }
    }
    cl->n = w;
}

/* Unit test for the --subgrid bbox filter (run via --selftest). */
static int grid_weld_selftest(void)
{
    Arena_T a = Arena_new();
    CubeList cl = {0, 0, 0};
    cubelist_push(a, &cl, "z04352_y03328_x02816"); /* in  */
    cubelist_push(a, &cl, "z04352_y02048_x01536"); /* out: y below  */
    cubelist_push(a, &cl, "z04736_y03840_x03328"); /* in  (upper corner) */
    cubelist_push(a, &cl, "z04352_y03328_x03456"); /* out: x above  */
    cubelist_push(a, &cl, "notacube");             /* out: unparseable */
    /* umbilicus block bbox: z[4352,4736] y[3328,3840] x[2816,3328] */
    cubelist_filter_bbox(&cl, 4352, 4736, 3328, 3840, 2816, 3328);
    int fails = 0;
    if (cl.n != 2) {
        fprintf(stderr, "[selftest] survivors: got %zu expect 2 -> FAIL\n", cl.n);
        fails++;
    } else {
        if (strcmp(cl.ids[0], "z04352_y03328_x02816") != 0) {
            fprintf(stderr, "[selftest] survivor[0]=%s -> FAIL\n", cl.ids[0]); fails++;
        }
        if (strcmp(cl.ids[1], "z04736_y03840_x03328") != 0) {
            fprintf(stderr, "[selftest] survivor[1]=%s -> FAIL\n", cl.ids[1]); fails++;
        }
    }
    if (!fails) fprintf(stderr, "[selftest] subgrid-bbox-filter -> ok (2 survivors)\n");
    Arena_dispose(&a);
    fprintf(stderr, "=== grid_weld selftest %s (%d failure%s) ===\n",
            fails ? "FAILED" : "PASSED", fails, fails == 1 ? "" : "s");
    return fails ? 1 : 0;
}

/* ===================================================================
 * Face dedup -- canonical edge representation = (min, max) sorted.
 * For manifold audit, count edges and per-edge face counts.
 *
 * Faces are deduped by sorted vert triple (a < b < c). After all cubes
 * loaded, identical triples (same a, b, c in any order) are merged into
 * one face -- this catches the case where two cubes' halos both emit
 * the same triangle at the seam.
 * =================================================================== */

typedef struct {
    int32_t a, b, c;       /* sorted ascending */
    int32_t orig0, orig1, orig2;  /* original winding (for emit) */
} SortedFace;

static void sort3i(int32_t v[3])
{
    if (v[0] > v[1]) { int32_t t = v[0]; v[0] = v[1]; v[1] = t; }
    if (v[1] > v[2]) { int32_t t = v[1]; v[1] = v[2]; v[2] = t; }
    if (v[0] > v[1]) { int32_t t = v[0]; v[0] = v[1]; v[1] = t; }
}

static int cmp_sorted_face(const void *pa, const void *pb)
{
    const SortedFace *fa = (const SortedFace *)pa;
    const SortedFace *fb = (const SortedFace *)pb;
    if (fa->a != fb->a) return (fa->a < fb->a) ? -1 : 1;
    if (fa->b != fb->b) return (fa->b < fb->b) ? -1 : 1;
    if (fa->c != fb->c) return (fa->c < fb->c) ? -1 : 1;
    return 0;
}

/* ===================================================================
 * Manifold audit -- build directed half-edges, sort by undirected key,
 * count run sizes: 1=unpaired, 2=manifold, >2=non_manifold; for runs of
 * 2, check same-direction = winding inversion.
 * =================================================================== */

typedef struct {
    int32_t src, dst;
} DHE;

static int cmp_dhe_undirected(const void *pa, const void *pb)
{
    const DHE *a = (const DHE *)pa;
    const DHE *b = (const DHE *)pb;
    int32_t a0 = (a->src < a->dst) ? a->src : a->dst;
    int32_t a1 = (a->src < a->dst) ? a->dst : a->src;
    int32_t b0 = (b->src < b->dst) ? b->src : b->dst;
    int32_t b1 = (b->src < b->dst) ? b->dst : b->src;
    if (a0 != b0) return (a0 < b0) ? -1 : 1;
    if (a1 != b1) return (a1 < b1) ? -1 : 1;
    return 0;
}

static int cmp_int64(const void *pa, const void *pb)
{
    int64_t a = *(const int64_t *)pa;
    int64_t b = *(const int64_t *)pb;
    if (a < b) return -1;
    if (a > b) return 1;
    return 0;
}

typedef struct {
    size_t unpaired;
    size_t non_manifold;
    size_t same_dir_pairs;
    size_t manifold_pairs;
} ManifoldStats;

static ManifoldStats manifold_audit(Arena_T arena,
                                    const int32_t *faces, size_t nf,
                                    /* If non-NULL, emit problem-edge
                                     * locations as line segments. */
                                    const float *verts_for_diag,
                                    FILE *diag_fp)
{
    ManifoldStats s = {0, 0, 0, 0};
    if (nf == 0) return s;
    size_t hn = nf * 3;
    DHE *he = (DHE *)ARENA_ALLOC(arena, (long)hn * (long)sizeof(DHE));
    for (size_t f = 0; f < nf; f++) {
        int32_t v0 = faces[f * 3 + 0];
        int32_t v1 = faces[f * 3 + 1];
        int32_t v2 = faces[f * 3 + 2];
        he[f * 3 + 0].src = v0; he[f * 3 + 0].dst = v1;
        he[f * 3 + 1].src = v1; he[f * 3 + 1].dst = v2;
        he[f * 3 + 2].src = v2; he[f * 3 + 2].dst = v0;
    }
    qsort(he, hn, sizeof(DHE), cmp_dhe_undirected);
    size_t diag_vidx = 1;
    if (diag_fp) {
        fprintf(diag_fp, "# manifold-audit problem edges\n");
        fprintf(diag_fp, "# unpaired = red, non_manifold = magenta, same_dir = yellow\n");
    }
    size_t i = 0;
    while (i < hn) {
        size_t j = i + 1;
        int32_t a0 = (he[i].src < he[i].dst) ? he[i].src : he[i].dst;
        int32_t a1 = (he[i].src < he[i].dst) ? he[i].dst : he[i].src;
        while (j < hn) {
            int32_t b0 = (he[j].src < he[j].dst) ? he[j].src : he[j].dst;
            int32_t b1 = (he[j].src < he[j].dst) ? he[j].dst : he[j].src;
            if (b0 != a0 || b1 != a1) break;
            j++;
        }
        size_t run = j - i;
        int kind = -1;  /* 0 unpaired, 1 non_manifold, 2 same_dir */
        if (run == 1) { s.unpaired++; kind = 0; }
        else if (run > 2) { s.non_manifold++; kind = 1; }
        else {
            int dir0 = (he[i].src < he[i].dst) ? 0 : 1;
            int dir1 = (he[i + 1].src < he[i + 1].dst) ? 0 : 1;
            if (dir0 == dir1) { s.same_dir_pairs++; kind = 2; }
            else s.manifold_pairs++;
        }
        if (diag_fp && verts_for_diag && (kind == 1 || kind == 2)) {
            /* Emit only the bug categories (skip outer-boundary unpaired
             * which would flood the diag). */
            float r = (kind == 1) ? 1.0f : 1.0f;
            float g = (kind == 1) ? 0.0f : 1.0f;
            float b = (kind == 1) ? 1.0f : 0.0f;
            int32_t p0 = a0, p1 = a1;
            fprintf(diag_fp, "v %.4f %.4f %.4f %.3f %.3f %.3f\n",
                (double)verts_for_diag[p0 * 3 + 0],
                (double)verts_for_diag[p0 * 3 + 1],
                (double)verts_for_diag[p0 * 3 + 2],
                (double)r, (double)g, (double)b);
            fprintf(diag_fp, "v %.4f %.4f %.4f %.3f %.3f %.3f\n",
                (double)verts_for_diag[p1 * 3 + 0],
                (double)verts_for_diag[p1 * 3 + 1],
                (double)verts_for_diag[p1 * 3 + 2],
                (double)r, (double)g, (double)b);
            fprintf(diag_fp, "l %zu %zu\n", diag_vidx, diag_vidx + 1);
            diag_vidx += 2;
        }
        i = j;
    }
    return s;
}

/* ===================================================================
 * Non-manifold neighbourhood OBJ -- triangles touching any non-manifold
 * edge, plus their 1-ring neighbours, written as a regular mesh. The
 * non-manifold faces are coloured magenta; the 1-ring is grey for
 * context. Lets you actually see what's going on at each bug, vs the
 * line-segment-only bad_edges.obj.
 * =================================================================== */
static void emit_nonmanifold_neighborhood(Arena_T arena,
                                          const int32_t *faces, size_t nf,
                                          const float *verts, size_t nv,
                                          const char *out_path)
{
    if (nf == 0 || nv == 0) return;

    /* 1) Rebuild the same sorted half-edge index manifold_audit uses, but
     *    keep face indices so we can map runs back to faces. */
    typedef struct { int32_t a0, a1; int32_t face; } DHE2;
    size_t hn = nf * 3;
    DHE2 *he = (DHE2 *)ARENA_ALLOC(arena, (long)hn * (long)sizeof(DHE2));
    for (size_t f = 0; f < nf; f++) {
        int32_t v[3] = { faces[f*3+0], faces[f*3+1], faces[f*3+2] };
        for (int e = 0; e < 3; e++) {
            int32_t s = v[e], d = v[(e+1)%3];
            he[f*3 + e].a0 = (s < d) ? s : d;
            he[f*3 + e].a1 = (s < d) ? d : s;
            he[f*3 + e].face = (int32_t)f;
        }
    }
    /* Sort by undirected edge. */
    qsort(he, hn, sizeof(DHE2), cmp_dhe_undirected);

    /* 2) Mark every face that owns at least one non-manifold edge. */
    uint8_t *nm_face = (uint8_t *)ARENA_CALLOC(arena, (long)nf, 1L);
    size_t i = 0;
    while (i < hn) {
        size_t j = i + 1;
        while (j < hn && he[j].a0 == he[i].a0 && he[j].a1 == he[i].a1) j++;
        size_t run = j - i;
        if (run > 2) {
            for (size_t k = i; k < j; k++) nm_face[he[k].face] = 1;
        }
        i = j;
    }

    /* 3) Vert use map: 1 = belongs to a non-manifold face (priority),
     *    2 = belongs to a 1-ring neighbour face only.
     *    A vert can be reached by both — keep the lower (1 wins). */
    uint8_t *vert_use = (uint8_t *)ARENA_CALLOC(arena, (long)nv, 1L);
    for (size_t f = 0; f < nf; f++) {
        if (!nm_face[f]) continue;
        for (int k = 0; k < 3; k++) vert_use[faces[f*3+k]] = 1;
    }
    /* 1-ring neighbours: faces that share at least one vert with a
     *  non-manifold face. */
    uint8_t *ring_face = (uint8_t *)ARENA_CALLOC(arena, (long)nf, 1L);
    for (size_t f = 0; f < nf; f++) {
        if (nm_face[f]) continue;
        for (int k = 0; k < 3; k++) {
            if (vert_use[faces[f*3+k]] == 1) { ring_face[f] = 1; break; }
        }
    }
    for (size_t f = 0; f < nf; f++) {
        if (!ring_face[f]) continue;
        for (int k = 0; k < 3; k++) {
            int32_t v = faces[f*3+k];
            if (vert_use[v] == 0) vert_use[v] = 2;
        }
    }

    /* 4) Compact verts that we'll emit, build old->new remap. */
    int32_t *remap = (int32_t *)ARENA_ALLOC(arena, (long)nv * (long)sizeof(int32_t));
    size_t out_nv = 0;
    for (size_t v = 0; v < nv; v++) {
        if (vert_use[v]) { remap[v] = (int32_t)(out_nv++); }
        else             { remap[v] = -1; }
    }
    if (out_nv == 0) return;

    FILE *fp = fopen(out_path, "w");
    if (!fp) return;
    fprintf(fp, "# non-manifold neighbourhood\n");
    fprintf(fp, "# magenta = faces on a non-manifold edge\n");
    fprintf(fp, "# grey    = 1-ring context faces\n");

    /* 5) Write verts. magenta for non-manifold-owned, grey for ring. */
    for (size_t v = 0; v < nv; v++) {
        if (!vert_use[v]) continue;
        float r, g, b;
        if (vert_use[v] == 1) { r = 1.0f; g = 0.0f; b = 1.0f; }
        else                  { r = 0.5f; g = 0.5f; b = 0.5f; }
        fprintf(fp, "v %.4f %.4f %.4f %.3f %.3f %.3f\n",
                (double)verts[v*3+0], (double)verts[v*3+1], (double)verts[v*3+2],
                (double)r, (double)g, (double)b);
    }
    /* 6) Write faces (1-indexed OBJ). */
    for (size_t f = 0; f < nf; f++) {
        if (!nm_face[f] && !ring_face[f]) continue;
        int32_t a = remap[faces[f*3+0]];
        int32_t b2 = remap[faces[f*3+1]];
        int32_t c = remap[faces[f*3+2]];
        fprintf(fp, "f %d %d %d\n", a + 1, b2 + 1, c + 1);
    }
    fclose(fp);
}

/* ===================================================================
 * Winding repair -- BFS-propagate orientation per connected component.
 *
 * For each connected component of the face-edge graph: pick the first
 * face as the orientation anchor and visit neighbors in BFS order. For
 * each neighbor sharing edge (u, v) with a parent face: the parent had
 * the edge in some direction (u -> v or v -> u); the neighbor must have
 * the opposite direction for them to be a manifold pair. If the
 * neighbor's direction is the same as the parent's, flip the neighbor
 * (reverse v1 <-> v2).
 *
 * Modifies flat_faces in place. Non-manifold edges (3+ faces) are
 * skipped during propagation (BFS only follows one neighbor per edge).
 * Disconnected components each get independently oriented; their
 * absolute orientation may not match a global "outward" convention but
 * within each component, all manifold edges become true pairs.
 *
 * Returns the number of faces flipped.
 * =================================================================== */

typedef struct {
    int32_t v0, v1;
    int32_t face;
} EdgeEntry;

static int cmp_edge_undirected(const void *pa, const void *pb)
{
    const EdgeEntry *a = (const EdgeEntry *)pa;
    const EdgeEntry *b = (const EdgeEntry *)pb;
    int32_t au = (a->v0 < a->v1) ? a->v0 : a->v1;
    int32_t av = (a->v0 < a->v1) ? a->v1 : a->v0;
    int32_t bu = (b->v0 < b->v1) ? b->v0 : b->v1;
    int32_t bv = (b->v0 < b->v1) ? b->v1 : b->v0;
    if (au != bu) return (au < bu) ? -1 : 1;
    if (av != bv) return (av < bv) ? -1 : 1;
    return 0;
}

static size_t repair_winding(Arena_T arena,
                             int32_t *faces, size_t nf,
                             size_t *out_components)
{
    if (nf == 0) {
        *out_components = 0;
        return 0;
    }
    /* Build per-face edge list with face-id back-pointer. */
    size_t hn = nf * 3;
    EdgeEntry *edges = (EdgeEntry *)ARENA_ALLOC(arena,
                          (long)hn * (long)sizeof(EdgeEntry));
    for (size_t f = 0; f < nf; f++) {
        int32_t v0 = faces[f * 3 + 0];
        int32_t v1 = faces[f * 3 + 1];
        int32_t v2 = faces[f * 3 + 2];
        edges[f * 3 + 0].v0 = v0; edges[f * 3 + 0].v1 = v1; edges[f * 3 + 0].face = (int32_t)f;
        edges[f * 3 + 1].v0 = v1; edges[f * 3 + 1].v1 = v2; edges[f * 3 + 1].face = (int32_t)f;
        edges[f * 3 + 2].v0 = v2; edges[f * 3 + 2].v1 = v0; edges[f * 3 + 2].face = (int32_t)f;
    }
    qsort(edges, hn, sizeof(EdgeEntry), cmp_edge_undirected);

    /* For each face, find its 3 adjacent (via shared-edge) faces via
     * the sorted edge list. Adjacency map: face -> up to 3 neighbors. */
    int32_t *adj = (int32_t *)ARENA_ALLOC(arena,
                       (long)nf * 3L * (long)sizeof(int32_t));
    /* Edge keeper: which edge does each adjacency entry use, expressed
     * as the (v0, v1) of the parent face's edge. */
    int32_t *adj_edge_u = (int32_t *)ARENA_ALLOC(arena,
                              (long)nf * 3L * (long)sizeof(int32_t));
    int32_t *adj_edge_v = (int32_t *)ARENA_ALLOC(arena,
                              (long)nf * 3L * (long)sizeof(int32_t));
    int32_t *adj_count = (int32_t *)ARENA_CALLOC(arena,
                             (long)nf, (long)sizeof(int32_t));
    for (size_t i = 0; i < hn; ) {
        size_t j = i + 1;
        int32_t au = (edges[i].v0 < edges[i].v1) ? edges[i].v0 : edges[i].v1;
        int32_t av = (edges[i].v0 < edges[i].v1) ? edges[i].v1 : edges[i].v0;
        while (j < hn) {
            int32_t bu = (edges[j].v0 < edges[j].v1) ? edges[j].v0 : edges[j].v1;
            int32_t bv = (edges[j].v0 < edges[j].v1) ? edges[j].v1 : edges[j].v0;
            if (bu != au || bv != av) break;
            j++;
        }
        size_t run = j - i;
        if (run == 2) {
            int32_t fA = edges[i].face;
            int32_t fB = edges[i + 1].face;
            if (adj_count[fA] < 3 && adj_count[fB] < 3) {
                adj[fA * 3 + adj_count[fA]] = fB;
                adj_edge_u[fA * 3 + adj_count[fA]] = edges[i].v0;
                adj_edge_v[fA * 3 + adj_count[fA]] = edges[i].v1;
                adj_count[fA]++;
                adj[fB * 3 + adj_count[fB]] = fA;
                adj_edge_u[fB * 3 + adj_count[fB]] = edges[i + 1].v0;
                adj_edge_v[fB * 3 + adj_count[fB]] = edges[i + 1].v1;
                adj_count[fB]++;
            }
        }
        /* run==1: boundary edge, no adjacency.
         * run>2: non-manifold, skip (propagation through these is undefined). */
        i = j;
    }

    /* BFS per connected component. */
    uint8_t *visited = (uint8_t *)ARENA_CALLOC(arena, (long)nf, 1L);
    int32_t *queue = (int32_t *)ARENA_ALLOC(arena,
                         (long)nf * (long)sizeof(int32_t));
    size_t n_flipped = 0;
    size_t n_components = 0;

    for (size_t seed = 0; seed < nf; seed++) {
        if (visited[seed]) continue;
        n_components++;
        size_t qh = 0, qt = 0;
        queue[qt++] = (int32_t)seed;
        visited[seed] = 1;
        while (qh < qt) {
            int32_t fcur = queue[qh++];
            int32_t cv0 = faces[fcur * 3 + 0];
            int32_t cv1 = faces[fcur * 3 + 1];
            int32_t cv2 = faces[fcur * 3 + 2];
            for (int k = 0; k < adj_count[fcur]; k++) {
                int32_t fnext = adj[fcur * 3 + k];
                if (visited[fnext]) continue;
                /* Shared edge, stored undirected (ea, eb). Recompute BOTH
                 * faces' current traversal direction from the live faces[]
                 * array. Caching the parent's directed edge at build time was
                 * a bug: once a parent is flipped its cached direction goes
                 * stale, so the neighbor test used the wrong sign and same_dir
                 * never dropped despite hundreds of flips. */
                int32_t ea = adj_edge_u[fcur * 3 + k];
                int32_t eb = adj_edge_v[fcur * 3 + k];
                int parent_ab = ((cv0 == ea && cv1 == eb) ||
                                 (cv1 == ea && cv2 == eb) ||
                                 (cv2 == ea && cv0 == eb));
                int32_t nv0 = faces[fnext * 3 + 0];
                int32_t nv1 = faces[fnext * 3 + 1];
                int32_t nv2 = faces[fnext * 3 + 2];
                int neigh_ab = ((nv0 == ea && nv1 == eb) ||
                                (nv1 == ea && nv2 == eb) ||
                                (nv2 == ea && nv0 == eb));
                /* Consistent winding => the two faces traverse the shared edge
                 * in OPPOSITE directions. Same direction => flip neighbor. */
                if (parent_ab == neigh_ab) {
                    int32_t tmp = faces[fnext * 3 + 1];
                    faces[fnext * 3 + 1] = faces[fnext * 3 + 2];
                    faces[fnext * 3 + 2] = tmp;
                    n_flipped++;
                }
                visited[fnext] = 1;
                queue[qt++] = fnext;
            }
        }
    }

    *out_components = n_components;
    return n_flipped;
}

/* Remove faces belonging to connected components smaller than min_verts unique
 * vertices, compacting faces in place. The seam re-BPA + eat-back sheds tiny
 * floating slivers (1-2 triangles) that are not part of any sheet; deleting
 * them clears their spurious boundary loops without touching the real surface
 * (sheets here are >=100 verts, garbage is <=4 -- any threshold between is
 * safe). Returns faces removed; sets *out_removed_comps to the CC count culled.
 * Face-graph adjacency built exactly as repair_winding (manifold edges only;
 * isolated/non-manifold-edge faces fall into their own components). */
static size_t cull_tiny_components(Arena_T arena, int32_t *faces, size_t nf,
                                   size_t min_verts, size_t *out_nf,
                                   size_t *out_removed_comps) {
    *out_removed_comps = 0;
    if (nf == 0) { *out_nf = 0; return 0; }

    size_t hn = nf * 3;
    EdgeEntry *edges = (EdgeEntry *)ARENA_ALLOC(arena,
                          (long)hn * (long)sizeof(EdgeEntry));
    for (size_t f = 0; f < nf; f++) {
        int32_t v0 = faces[f*3+0], v1 = faces[f*3+1], v2 = faces[f*3+2];
        edges[f*3+0].v0=v0; edges[f*3+0].v1=v1; edges[f*3+0].face=(int32_t)f;
        edges[f*3+1].v0=v1; edges[f*3+1].v1=v2; edges[f*3+1].face=(int32_t)f;
        edges[f*3+2].v0=v2; edges[f*3+2].v1=v0; edges[f*3+2].face=(int32_t)f;
    }
    qsort(edges, hn, sizeof(EdgeEntry), cmp_edge_undirected);

    int32_t *adj = (int32_t *)ARENA_ALLOC(arena,
                       (long)nf * 3L * (long)sizeof(int32_t));
    int32_t *adj_count = (int32_t *)ARENA_CALLOC(arena,
                             (long)nf, (long)sizeof(int32_t));
    for (size_t i = 0; i < hn; ) {
        size_t j = i + 1;
        int32_t au = (edges[i].v0 < edges[i].v1) ? edges[i].v0 : edges[i].v1;
        int32_t av = (edges[i].v0 < edges[i].v1) ? edges[i].v1 : edges[i].v0;
        while (j < hn) {
            int32_t bu = (edges[j].v0 < edges[j].v1) ? edges[j].v0 : edges[j].v1;
            int32_t bv = (edges[j].v0 < edges[j].v1) ? edges[j].v1 : edges[j].v0;
            if (bu != au || bv != av) break;
            j++;
        }
        if (j - i == 2) {
            int32_t fA = edges[i].face, fB = edges[i + 1].face;
            if (adj_count[fA] < 3 && adj_count[fB] < 3) {
                adj[fA * 3 + adj_count[fA]++] = fB;
                adj[fB * 3 + adj_count[fB]++] = fA;
            }
        }
        i = j;
    }

    /* BFS components; record each face's component and the comp's vert count. */
    int32_t *comp = (int32_t *)ARENA_ALLOC(arena,
                        (long)nf * (long)sizeof(int32_t));
    for (size_t f = 0; f < nf; f++) { comp[f] = -1; }
    int32_t *queue = (int32_t *)ARENA_ALLOC(arena,
                         (long)nf * (long)sizeof(int32_t));
    /* vert -> last comp that counted it, so each vert is tallied once per comp */
    size_t max_v = 0;
    for (size_t k = 0; k < nf * 3; k++) {
        if ((size_t)faces[k] + 1 > max_v) { max_v = (size_t)faces[k] + 1; }
    }
    int32_t *vseen = (int32_t *)ARENA_ALLOC(arena,
                         (long)max_v * (long)sizeof(int32_t));
    for (size_t v = 0; v < max_v; v++) { vseen[v] = -1; }
    size_t *comp_verts = (size_t *)ARENA_ALLOC(arena,
                             (long)nf * (long)sizeof(size_t)); /* <=nf comps */

    size_t n_comp = 0;
    for (size_t seed = 0; seed < nf; seed++) {
        if (comp[seed] != -1) { continue; }
        int32_t cid = (int32_t)n_comp;
        size_t vcount = 0, qh = 0, qt = 0;
        queue[qt++] = (int32_t)seed; comp[seed] = cid;
        while (qh < qt) {
            int32_t fcur = queue[qh++];
            for (int e = 0; e < 3; e++) {
                int32_t v = faces[fcur*3+e];
                if (vseen[v] != cid) { vseen[v] = cid; vcount++; }
            }
            for (int k = 0; k < adj_count[fcur]; k++) {
                int32_t fnext = adj[fcur * 3 + k];
                if (comp[fnext] != -1) { continue; }
                comp[fnext] = cid; queue[qt++] = fnext;
            }
        }
        comp_verts[n_comp++] = vcount;
    }

    /* Compact: drop faces whose component is too small. */
    size_t w = 0, removed_comps = 0;
    for (size_t c = 0; c < n_comp; c++) {
        if (comp_verts[c] < min_verts) { removed_comps++; }
    }
    for (size_t f = 0; f < nf; f++) {
        if (comp_verts[comp[f]] < min_verts) { continue; }
        faces[w*3+0] = faces[f*3+0];
        faces[w*3+1] = faces[f*3+1];
        faces[w*3+2] = faces[f*3+2];
        w++;
    }
    *out_nf = w;
    *out_removed_comps = removed_comps;
    return nf - w;
}

/* Set an environment variable only if it is not already set (an explicit
 * caller-set value always wins). Used so --pair auto-enables the SEAM_DUMP_FRONT
 * init-front diagnostic without the user having to set it by hand. */
static void set_env_if_unset(const char *name, const char *val)
{
    if (getenv(name)) return;
#ifdef _MSC_VER
    _putenv_s(name, val);
#else
    setenv(name, val, 1);
#endif
}

/* ===================================================================
 * Per-stage OBJ dump (--dump-stages). Writes the welded mesh as it is
 * after one weld stage to <dir>/<prefix>_<NN>_<name>.obj, colored by the
 * originating cube; verts past color_nv (added by pinhole/holefill) are
 * light grey. No-op when dir is NULL. Each stage is written as it
 * completes, so a hang/crash in a later stage still leaves every earlier
 * stage on disk. Uses a scratch arena mark for the color array.
 * =================================================================== */
static void dump_stage(Arena_T arena, const char *dir, const char *prefix,
                       int idx, const char *name,
                       const float *verts, size_t nv,
                       const int32_t *faces, size_t nf,
                       const int16_t *vert_cube_idx,
                       const int8_t *cube_palette,
                       size_t n_cubes, size_t color_nv)
{
    if (!dir) return;
    static const float GREY[3] = { 0.85f, 0.85f, 0.85f };
    Arena_Mark m = Arena_save(arena);
    float *colors = (float *)ARENA_ALLOC(arena,
                       (long)(nv * 3L * (long)sizeof(float)));
    /* Default: colour by CONNECTED COMPONENT (wrap shatter / fusion shows on
     * load). Set GW_COLOR_BY=cube for the old per-cube provenance palette
     * (useful for seam-weld debugging). */
    const char *gw_color_by = getenv("GW_COLOR_BY");
    if (gw_color_by && !strcmp(gw_color_by, "cube")) {
        for (size_t v = 0; v < nv; v++) {
            const float *c = GREY;
            if (v < color_nv) {
                int cube = vert_cube_idx[v];
                int pal = (cube >= 0 && cube < (int)n_cubes) ? cube_palette[cube] : 0;
                c = CUBE_PALETTE[pal];
            }
            colors[v * 3 + 0] = c[0];
            colors[v * 3 + 1] = c[1];
            colors[v * 3 + 2] = c[2];
        }
    } else {
        CCColorOpts opts; CCColor_default_opts(&opts); opts.sat = 0.75;
        CCColor_compute(nv, faces, nf, &opts, colors, NULL);
    }
    char path[1024];
    snprintf(path, sizeof(path), "%s/%s_%02d_%s.obj", dir, prefix, idx, name);
    ves_ensure_parent_dir(path);
    int rc = ObjIO_write_per_vertex_color(path, verts, nv, faces, nf, colors);
    fprintf(stderr, "  [dump-stage %02d] %-11s %zu v %zu f -> %s%s\n",
            idx, name, nv, nf, path, rc == 0 ? "" : "  [WRITE FAILED]");
    fflush(stderr);
    Arena_restore(arena, m);
}

/* ===================================================================
 * Main
 * =================================================================== */

/* Phase wall-clock: prints the time since the previous checkpoint. The weld
 * is single-process and serial, so a running phase profile in the log is the
 * whole story ("full logging" house rule; grid_weld.log keeps it). */
static double gw_phase_prev = -1.0;
static void gw_phase(const char *name)
{
    double now = ves_clock_sec();
    if (gw_phase_prev >= 0.0) {
        fprintf(stderr, "  [time] %-12s %7.2fs\n", name, now - gw_phase_prev);
    }
    gw_phase_prev = now;
}

int main(int argc, char **argv)
{
    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0) return grid_weld_selftest();
    if (argc < 3) {
        fprintf(stderr,
            "Usage: %s <grid_obj_dir> <output.obj> [--stage <name>] [--emit-weld-verts <path>]\n"
            "       %s <grid_obj_dir> <output.obj> --pair <cubeA_id> <compA> <cubeB_id> <compB> [--stage <name>]\n"
            "\n"
            "Reads <grid_obj_dir>/<cube_id>/<cube_id>_<stage>/<cube_id>_<stage>_all.obj\n"
            "for every <cube_id> subdirectory. Plain-concatenates the cubes and\n"
            "welds them by re-running BPA across each cube-boundary seam plane.\n"
            "Writes welded mesh to <output.obj> and a JSON audit report\n"
            "to <output.obj>.weld_report.json.\n"
            "\n"
            "--pair <cubeA_id> <compA> <cubeB_id> <compB>: debug mode -- weld just\n"
            "                two individual components (the _<comp>.obj files, comp\n"
            "                1-based) instead of the whole grid. Auto-dumps the\n"
            "                BPA init front to <output.obj>.front_{tris,edges}.obj.\n"
            "\n"
            "--stage <name>: which dump stage to read (default: step12_final).\n"
            "                e.g. step1_bpa, step7_cc_bpa, step12_final.\n"
            "\n"
            "--subgrid z0 z1 y0 y1 x0 x1: weld only cubes whose origin lies in this\n"
            "                inclusive source-voxel bbox -- stitches one rectangular\n"
            "                sub-block of a larger dump dir (e.g. a 4x5x5 tile).\n"
            "\n"
            "--dump-stages <dir>: write each weld stage to its own colored OBJ\n"
            "                <dir>/<prefix>_NN_<stage>.obj (concat..final), plus a\n"
            "                per-hole input/result OBJ to <dir>/holes/ for every\n"
            "                interior hole filled. Shows which stage -- or which\n"
            "                hole -- breaks. (Verbose: use on a small region.)\n"
            "--stage-prefix <id>: filename prefix for --dump-stages (default weld).\n"
            "\n"
            "Optional --emit-weld-verts <path> writes a packed binary sidecar:\n"
            "  uint32 nv; uint8 is_weld[nv];\n"
            "Used by the winding diagnostic to compute its headline weld-vert\n"
            "disagreement rate.\n",
            argv[0], argv[0]);
        return 1;
    }
    const char *grid_dir = argv[1];
    const char *out_path = argv[2];
    const char *weld_verts_out = NULL;
    const char *stage = "step12_final";
    const char *dump_stages_dir = NULL; /* --dump-stages <dir>: one OBJ per weld stage */
    const char *stage_prefix = "weld";  /* --stage-prefix <id>: dump filename prefix */
    float seam_cube = CUBE_SIZE_VOX;    /* seam-plane spacing (voxels) */
    int pair_mode = 0;                  /* --pair: weld two components only */
    const char *pair_cube[2] = { NULL, NULL };
    int pair_comp[2] = { 0, 0 };        /* 1-based component ids */
    int no_bridge = 0;                  /* --no-bridge: concat only (stage 1) */
    int no_pinhole = (getenv("SEAM_NO_PINHOLE") != NULL); /* --no-pinhole: stage 2 */
    int no_cleanup = (getenv("SEAM_NO_CLEANUP") != NULL); /* --no-cleanup: skip post-weld flip+collapse */
    int no_holefill = (getenv("SEAM_NO_HOLEFILL") != NULL); /* --no-holefill: skip interior-hole fill */
    int subgrid = 0;            /* --subgrid z0 z1 y0 y1 x0 x1: weld one origin-bbox block only */
    int64_t sg_z0 = 0, sg_z1 = 0, sg_y0 = 0, sg_y1 = 0, sg_x0 = 0, sg_x1 = 0;
    /* --vert-cap/--face-cap: override the node-count-derived capacity. The
     * default estimate (1.5M vert / 2M face per input node) is calibrated for
     * leaf-sized cubes; an UNDECIMATED upper-level weld (hierarchical_weld
     * --no-decimate) folds the whole scroll into a handful of huge nodes, so a
     * 4-node terminal weld needs far more than 4*2M faces. The caller passes
     * the real budget here; 0 = use the estimate. Still clamped to the 32-bit
     * alloc ceilings (150M vert / 80M face). */
    size_t vert_cap_override = 0, face_cap_override = 0;
    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--emit-weld-verts") == 0 && i + 1 < argc) {
            weld_verts_out = argv[++i];
        } else if (strcmp(argv[i], "--no-bridge") == 0) {
            no_bridge = 1;
        } else if (strcmp(argv[i], "--no-pinhole") == 0) {
            no_pinhole = 1;
        } else if (strcmp(argv[i], "--no-cleanup") == 0) {
            no_cleanup = 1;
        } else if (strcmp(argv[i], "--no-holefill") == 0) {
            no_holefill = 1;
        } else if (strcmp(argv[i], "--stage") == 0 && i + 1 < argc) {
            stage = argv[++i];
        } else if (strcmp(argv[i], "--dump-stages") == 0 && i + 1 < argc) {
            dump_stages_dir = argv[++i];
        } else if (strcmp(argv[i], "--stage-prefix") == 0 && i + 1 < argc) {
            stage_prefix = argv[++i];
        } else if (strcmp(argv[i], "--pair") == 0 && i + 4 < argc) {
            pair_mode = 1;
            pair_cube[0] = argv[++i];
            pair_comp[0] = atoi(argv[++i]);
            pair_cube[1] = argv[++i];
            pair_comp[1] = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--seam-weld") == 0) {
            /* Deprecated no-op: the BPA seam-weld is unconditional now (the
             * bitwise hash-join weld it replaced is gone). Accepted so
             * existing callers/scripts don't error. */
        } else if (strcmp(argv[i], "--cube-size") == 0 && i + 1 < argc) {
            seam_cube = (float)atof(argv[++i]);
        } else if (strcmp(argv[i], "--subgrid") == 0 && i + 6 < argc) {
            subgrid = 1;
            sg_z0 = atoll(argv[++i]); sg_z1 = atoll(argv[++i]);
            sg_y0 = atoll(argv[++i]); sg_y1 = atoll(argv[++i]);
            sg_x0 = atoll(argv[++i]); sg_x1 = atoll(argv[++i]);
        } else if (strcmp(argv[i], "--vert-cap") == 0 && i + 1 < argc) {
            vert_cap_override = (size_t)atoll(argv[++i]);
        } else if (strcmp(argv[i], "--face-cap") == 0 && i + 1 < argc) {
            face_cap_override = (size_t)atoll(argv[++i]);
        } else {
            fprintf(stderr, "grid_weld: unknown arg %s\n", argv[i]);
            return 1;
        }
    }

    /* --dump-stages also turns on per-hole OBJ dumps inside the interior
     * hole-fill: each hole's input boundary is written BEFORE the fill (so a
     * hole that hangs Triangle leaves its polygon as the last *_in.obj on disk),
     * and the filled patch is written on success. hole_fill.c reads these env
     * vars; an explicit caller-set value always wins. */
    char hole_dump_dir[1024];
    if (dump_stages_dir) {
        snprintf(hole_dump_dir, sizeof(hole_dump_dir), "%s/holes", dump_stages_dir);
        set_env_if_unset("SEAM_HOLE_DUMP_DIR", hole_dump_dir);
        set_env_if_unset("SEAM_HOLE_DUMP_PREFIX", stage_prefix);
    }

    Arena_T arena = Arena_new();
    int ok = 0;
    /* Final accumulators (sized after cube enumeration). */
    float *out_verts = NULL;
    size_t out_nv = 0;
    SortedFace *all_faces = NULL;
    size_t all_nf = 0;
    size_t face_cap = 0;
    ManifoldStats ms = {0, 0, 0, 0};
    size_t n_unique_faces = 0;
    CubeList cubes = {0, 0, 0};

    TRY
        gw_phase("start");
        /* Cube set. Pair mode welds exactly two named components; grid mode
         * enumerates every cube subdirectory. Pushing the two cube IDs into the
         * same CubeList lets every sizing/palette/coloring/JSON path below run
         * unchanged with cubes.n == 2. */
        if (pair_mode) {
            if (pair_comp[0] <= 0 || pair_comp[1] <= 0) {
                fprintf(stderr, "grid_weld: --pair comp ids must be >= 1\n");
                RAISE(IO_Failed);
            }
            cubelist_push(arena, &cubes, pair_cube[0]);
            cubelist_push(arena, &cubes, pair_cube[1]);
            fprintf(stderr,
                "grid_weld: pair mode -- %s comp %03d + %s comp %03d (stage %s)\n",
                pair_cube[0], pair_comp[0], pair_cube[1], pair_comp[1], stage);
        } else {
            if (enumerate_cube_dirs(arena, grid_dir, &cubes) != 0) {
                fprintf(stderr, "grid_weld: cannot enumerate %s\n", grid_dir);
                RAISE(IO_Failed);
            }
            if (cubes.n == 0) {
                fprintf(stderr, "grid_weld: no cube directories under %s\n",
                        grid_dir);
                RAISE(IO_Failed);
            }
            fprintf(stderr, "grid_weld: %zu cubes under %s\n", cubes.n, grid_dir);
            if (subgrid) {
                cubelist_filter_bbox(&cubes, sg_z0, sg_z1, sg_y0, sg_y1, sg_x0, sg_x1);
                fprintf(stderr,
                    "grid_weld: --subgrid z[%lld,%lld] y[%lld,%lld] x[%lld,%lld] -> %zu cubes\n",
                    (long long)sg_z0, (long long)sg_z1, (long long)sg_y0,
                    (long long)sg_y1, (long long)sg_x0, (long long)sg_x1, cubes.n);
                if (cubes.n == 0) {
                    fprintf(stderr, "grid_weld: --subgrid matched no cubes\n");
                    RAISE(IO_Failed);
                }
            }
        }

        /* Output vert array. Plain concatenation of every cube's verts (the
         * bitwise hash-join weld is gone), so size for the no-reuse worst
         * case. Per-cube dumps are ~30K-180K verts, but grid_weld is also used
         * to stitch already-welded *block* meshes: for hierarchical_weld LOD
         * tiers the boundary-pinned decimated blocks stay large (300K-700K+
         * verts each at upper levels), so a weld of just 2 such tiles blew the
         * old 2*500K=1M cap. 1.5M/input covers 2-input upper-level welds with
         * margin, plus hole-fill Steiner verts. The 150M ceiling still bounds a
         * full grid. */
        size_t vert_per_cube_est = 1500000;
        size_t vert_cap = vert_cap_override ? vert_cap_override
                                            : cubes.n * vert_per_cube_est;
        if (vert_cap < (1 << 16)) vert_cap = (1 << 16);
        /* Cap at 32-bit alloc limit: vert_cap * 12 bytes < 2 GB --> 178M. */
        if (vert_cap > 150000000) vert_cap = 150000000;
        out_verts = (float *)ARENA_ALLOC(arena,
                       (long)(vert_cap * 3L * (long)sizeof(float)));
        /* Track which cube first inserted each vert (for per-cube
         * color assignment in the output OBJ). int16 supports up to 32K
         * cubes -- well beyond any realistic grid we'd put in memory. */
        int16_t *vert_cube_idx = (int16_t *)ARENA_ALLOC(arena,
                                     (long)vert_cap * (long)sizeof(int16_t));
        /* Weld bitmap: 1 if this vert was inserted by one cube and later
         * reused by a different cube (i.e., an actual cross-cube weld).
         * Distinct from "vert is near a cube-boundary plane" -- the latter
         * also includes outer-boundary verts where no neighbor exists and
         * verts on standalone components that don't extend across the
         * boundary. We only want the former. */
        uint8_t *vert_is_weld = (uint8_t *)ARENA_CALLOC(arena,
                                    (long)vert_cap, 1L);
        /* Cube-index -> palette-index map. Use the integer cube-grid
         * coordinates (cz, cy, cx) with primes per axis that don't share
         * factors with the palette size (16). Adjacent cubes (differ by
         * one in any axis) always get a different palette index, so
         * neighbors visually contrast. Deterministic across runs and
         * grids. */
        int8_t *cube_palette = (int8_t *)ARENA_ALLOC(arena,
                                   (long)cubes.n * (long)sizeof(int8_t));
        for (size_t i = 0; i < cubes.n; i++) {
            int64_t vz = 0, vy = 0, vx = 0;
            (void)parse_cube_origin(cubes.ids[i], &vz, &vy, &vx);
            /* Convert source-voxel offsets to integer cube indices. */
            int32_t cz = (int32_t)(vz / 128);
            int32_t cy = (int32_t)(vy / 128);
            int32_t cx = (int32_t)(vx / 128);
            /* 7, 3, 1: each coprime to 16 and to each other for axis
             * independence. Adjacent cubes differ by exactly one of
             * {7, 3, 1} mod 16, all nonzero. */
            uint32_t idx = (uint32_t)(cz * 7 + cy * 3 + cx);
            cube_palette[i] = (int8_t)(idx & 0xFu);
        }

        /* Face accumulator. ~400K/cube for per-cube dumps; raised to 2M to also
         * cover the large boundary-pinned LOD tiles hierarchical_weld welds at
         * upper levels (~300-700K each, matching the vert bump above) plus the
         * bridge + hole-fill faces the weld adds. The 80M ceiling still bounds a
         * grid. */
        face_cap = face_cap_override ? face_cap_override : cubes.n * 2000000;
        if (face_cap < (1 << 16)) face_cap = (1 << 16);
        /* Cap at 32-bit alloc limit: face_cap * 24 bytes < 2 GB --> 87M. */
        if (face_cap > 80000000) face_cap = 80000000;
        all_faces = (SortedFace *)ARENA_ALLOC(arena,
                       (long)(face_cap * (long)sizeof(SortedFace)));

        /* Process each cube. */
        Arena_T scratch = Arena_new();
        size_t total_in_verts = 0;
        size_t total_in_faces = 0;
        size_t total_unique_verts = 0;

        for (size_t ci = 0; ci < cubes.n; ci++) {
            const char *cube_id = cubes.ids[ci];
            int64_t vz = 0, vy = 0, vx = 0;
            if (parse_cube_origin(cube_id, &vz, &vy, &vx) != 0) {
                fprintf(stderr, "grid_weld: bad cube_id %s (skip)\n", cube_id);
                continue;
            }
            char obj_path[1024];
            if (pair_mode) {
                /* One specific component (_NNN.obj), not the merged _all.obj. */
                snprintf(obj_path, sizeof(obj_path),
                    "%s/%s/%s_%s/%s_%s_%03d.obj",
                    grid_dir, cube_id, cube_id, stage, cube_id, stage,
                    pair_comp[ci]);
            } else {
                snprintf(obj_path, sizeof(obj_path),
                    "%s/%s/%s_%s/%s_%s_all.obj",
                    grid_dir, cube_id, cube_id, stage, cube_id, stage);
                if (!gw_file_exists(obj_path)) {
                    /* Flat hierarchical-LOD layout (hierarchical_weld):
                     * <dir>/<id>/<id>_<stage>_all.obj */
                    snprintf(obj_path, sizeof(obj_path),
                        "%s/%s/%s_%s_all.obj",
                        grid_dir, cube_id, cube_id, stage);
                }
            }

            Arena_free(scratch);
            float *vin = NULL;
            int32_t *fin = NULL;
            size_t nv_in = 0, nf_in = 0;
            if (ObjIO_read(scratch, obj_path, &vin, &nv_in, &fin, &nf_in) != 0) {
                fprintf(stderr, "  cube %s: OBJ read failed (%s) -- %s\n",
                        cube_id, obj_path, pair_mode ? "FATAL" : "skip");
                if (pair_mode) RAISE(IO_Failed);  /* both comps required */
                continue;
            }
            total_in_verts += nv_in;
            total_in_faces += nf_in;

            /* Plain concatenation: append every vert with an offset remap.
             * The bitwise hash-join weld is gone -- step0 trims each cube to
             * its owned [0,cube] box, so adjacent cubes' seam verts were never
             * coincident anyway; SeamWeld_bridge (BPA) fuses the seam below.
             * Per-cube OBJs are already in source-voxel WORLD coords. */
            (void)vz; (void)vy; (void)vx;
            int32_t *remap = (int32_t *)ARENA_ALLOC(scratch,
                                (long)nv_in * (long)sizeof(int32_t));
            if (out_nv + nv_in > vert_cap) {
                fprintf(stderr,
                    "grid_weld: vert capacity %zu exceeded\n", vert_cap);
                RAISE(IO_Failed);
            }
            for (size_t v = 0; v < nv_in; v++) {
                int32_t new_idx = (int32_t)out_nv;
                out_verts[out_nv * 3 + 0] = vin[v * 3 + 0];
                out_verts[out_nv * 3 + 1] = vin[v * 3 + 1];
                out_verts[out_nv * 3 + 2] = vin[v * 3 + 2];
                vert_cube_idx[out_nv] = (int16_t)ci;
                out_nv++;
                remap[v] = new_idx;
            }
            total_unique_verts += nv_in;

            /* Append remapped faces (sorted ascending for later dedup). */
            for (size_t f = 0; f < nf_in; f++) {
                int32_t a = remap[fin[f * 3 + 0]];
                int32_t b = remap[fin[f * 3 + 1]];
                int32_t c = remap[fin[f * 3 + 2]];
                if (a == b || b == c || a == c) continue;  /* degenerate */
                int32_t sorted[3] = { a, b, c };
                sort3i(sorted);
                if (all_nf >= face_cap) {
                    fprintf(stderr,
                        "grid_weld: face capacity %zu exceeded\n", face_cap);
                    RAISE(IO_Failed);
                }
                all_faces[all_nf].a = sorted[0];
                all_faces[all_nf].b = sorted[1];
                all_faces[all_nf].c = sorted[2];
                all_faces[all_nf].orig0 = a;
                all_faces[all_nf].orig1 = b;
                all_faces[all_nf].orig2 = c;
                all_nf++;
            }

            fprintf(stderr, "  cube %s: %zu v, %zu f (concatenated)\n",
                    cube_id, nv_in, nf_in);
        }
        Arena_dispose(&scratch);

        fprintf(stderr,
            "grid_weld: aggregate: %zu v_in -> %zu unique, %zu f_in pre-dedup\n",
            total_in_verts, out_nv, all_nf);
        gw_phase("load");

        /* Pair mode: print each component's world bbox so it is obvious whether
         * the two straddle a cube-boundary plane (a multiple of seam_cube on one
         * axis). If they do not, detect_planes finds nothing and the weld is a
         * no-op -- the warning after the bridge call confirms that case. */
        if (pair_mode) {
            for (size_t ci2 = 0; ci2 < cubes.n; ci2++) {
                float lo[3] = { 1e30f, 1e30f, 1e30f };
                float hi[3] = { -1e30f, -1e30f, -1e30f };
                size_t cnt = 0;
                for (size_t v = 0; v < out_nv; v++) {
                    if (vert_cube_idx[v] != (int16_t)ci2) continue;
                    for (int a = 0; a < 3; a++) {
                        float c = out_verts[v * 3 + a];
                        if (c < lo[a]) lo[a] = c;
                        if (c > hi[a]) hi[a] = c;
                    }
                    cnt++;
                }
                fprintf(stderr,
                    "  pair[%zu] %s comp %03d: %zu v, bbox "
                    "z[%.1f,%.1f] y[%.1f,%.1f] x[%.1f,%.1f]\n",
                    ci2, cubes.ids[ci2], pair_comp[ci2], cnt,
                    (double)lo[0], (double)hi[0], (double)lo[1], (double)hi[1],
                    (double)lo[2], (double)hi[2]);
            }
        }

        /* Sort + dedup faces by (a, b, c). */
        qsort(all_faces, all_nf, sizeof(SortedFace), cmp_sorted_face);
        size_t write_i = 0;
        for (size_t i = 0; i < all_nf; ) {
            size_t j = i + 1;
            while (j < all_nf &&
                   all_faces[j].a == all_faces[i].a &&
                   all_faces[j].b == all_faces[i].b &&
                   all_faces[j].c == all_faces[i].c) {
                j++;
            }
            if (write_i != i) all_faces[write_i] = all_faces[i];
            write_i++;
            i = j;
        }
        n_unique_faces = write_i;
        fprintf(stderr, "grid_weld: %zu unique faces (deduped from %zu)\n",
                n_unique_faces, all_nf);
        gw_phase("facesort");

        /* Run manifold audit on welded mesh. */
        int32_t *flat_faces = (int32_t *)ARENA_ALLOC(arena,
                                  (long)n_unique_faces * 3L * (long)sizeof(int32_t));
        for (size_t f = 0; f < n_unique_faces; f++) {
            flat_faces[f * 3 + 0] = all_faces[f].orig0;
            flat_faces[f * 3 + 1] = all_faces[f].orig1;
            flat_faces[f * 3 + 2] = all_faces[f].orig2;
        }

        /* Seam-weld stage (unconditional -- this IS the weld now). The BPA
         * mesh path does not pin seam verts and the cubes are plain-
         * concatenated, so the overlap is two independent, overlapping
         * sheets. Re-run BPA across each detected cube-boundary plane to fuse
         * them into one manifold (user's "treat the halo as a new BPA weld"),
         * then close the micro-holes the greedy roll leaves with the
         * dedicated pinhole pass (NOT CDT). This runs BEFORE winding repair +
         * audit, which then operate on the fused mesh. `color_nv` freezes the
         * pre-seam count so the per-cube vertex palette below stays in range
         * if pinhole adds split verts. */
        size_t color_nv = out_nv;
        /* Stage 00: plain concatenation + face dedup, before any weld. */
        dump_stage(arena, dump_stages_dir, stage_prefix, 0, "concat",
                   out_verts, out_nv, flat_faces, n_unique_faces,
                   vert_cube_idx, cube_palette, cubes.n, color_nv);
        size_t recoarsen_collapses = 0, recoarsen_faces_in = 0,
               recoarsen_faces_out = 0;
        size_t bandcvt_accepted = 0, bandcvt_rejected = 0,
               bandcvt_faces_in = 0, bandcvt_faces_out = 0;
        if (!no_bridge) {              /* --no-bridge: stage-1 concat, no weld */
            /* band 6 (was 4): the bridge front considers boundary edges within
             * this many vox of the seam plane. At 4, near-seam material whose
             * boundary sits 4-6 vox out (a slightly wider gap or off-centre seam)
             * was never fed to the bridge, leaving an unfilled strip. 6 catches it
             * while staying < the 7-vox inter-wrap clearance, and the over-long
             * guard (span <= 2*rho_max = 6) + fold guard make a merger impossible.
             * Measured (2x1x1): bridge faces 1059->1075, seam-band open edges
             * 259->239, overlap/fold/non-manifold all still 0. SEAM_BAND overrides. */
            float rho = 1.5f, band = 6.0f;
            const char *e;
            if ((e = getenv("SEAM_RHO")))  rho  = (float)atof(e);
            if ((e = getenv("SEAM_BAND"))) band = (float)atof(e);
            /* Phase-1 winding gate, from env (the gate is now a param, not
             * env-read inside the bridge). Armed only when the umbilicus +
             * pitch are supplied; grid runners pass them. */
            BpaBridgeGate seam_gate;
            memset(&seam_gate, 0, sizeof seam_gate);
            if ((e = getenv("SEAM_UMBILICUS_Y"))) seam_gate.umb_y = atof(e);
            if ((e = getenv("SEAM_UMBILICUS_X"))) seam_gate.umb_x = atof(e);
            if ((e = getenv("SEAM_WRAP_PITCH"))) { double v = atof(e); if (v > 0) seam_gate.pitch = v; }
            if (seam_gate.pitch > 0.0 && (seam_gate.umb_y != 0.0 || seam_gate.umb_x != 0.0)) {
                seam_gate.tol  = SEAM_WIND_TOL_DEFAULT_TURNS;
                seam_gate.hard = SEAM_WIND_HARD_TOL_DEFAULT_TURNS;
                if ((e = getenv("SEAM_WIND_TOL")))      { double v = atof(e); if (v > 0) seam_gate.tol  = v; }
                if ((e = getenv("SEAM_WIND_HARD_TOL"))) { double v = atof(e); if (v > 0) seam_gate.hard = v; }
            }
            /* Pair-debug: always emit the BPA init front next to the output so
             * missed seam edges (red) are visible without setting env by hand.
             * An explicitly-set SEAM_DUMP_FRONT wins. */
            if (pair_mode) set_env_if_unset("SEAM_DUMP_FRONT", out_path);

            /* Weld-time seam-band refinement (default ON): subdivide + flip
             * the band next to each detected seam plane into well-shaped
             * ~3.5-vox triangles so the bridge can prime on uniform-coarse
             * CVT cubes (a ~13-vox triangle's circumradius ~7 > rho_max=3
             * fails sphere_center -> zero bridges; see seam_refine.h). New
             * verts are chord midpoints carrying their source vert's cube
             * provenance; the band is collapsed back out by the recoarsen
             * stage after all closers. On dense or graded-rim inputs the
             * band edges are already <= target, so this no-ops. SEAM_NO_REFINE
             * disables; SEAM_REFINE_TARGET tunes. */
            if (!getenv("SEAM_NO_REFINE")) {
                uint8_t *rf_used = (uint8_t *)ARENA_CALLOC(arena, (long)out_nv, 1L);
                for (size_t f = 0; f < n_unique_faces; f++) {
                    rf_used[flat_faces[f*3+0]] = 1;
                    rf_used[flat_faces[f*3+1]] = 1;
                    rf_used[flat_faces[f*3+2]] = 1;
                }
                SeamPlane fplanes[64];
                size_t fnp = SeamPlanes_detect(out_verts, out_nv, rf_used,
                                               (double)seam_cube, (double)band,
                                               fplanes, 64);
                if (fnp > 0) {
                    SeamRefineParams sp; SeamRefine_default_params(&sp);
                    sp.band = band;
                    if ((e = getenv("SEAM_REFINE_TARGET"))) {
                        double v = atof(e); if (v > 0) sp.target_len = (float)v;
                    }
                    SeamRefineStats sst;
                    float *rf_v = NULL; int32_t *rf_f = NULL;
                    int32_t *rf_src = NULL;
                    size_t rf_nv = 0, rf_nf = 0, rf_nnew = 0;
                    if (SeamRefine_process(arena, out_verts, out_nv,
                                           flat_faces, n_unique_faces,
                                           fplanes, fnp, &sp,
                                           &rf_v, &rf_nv, &rf_f, &rf_nf,
                                           &rf_src, &rf_nnew, &sst) == 0
                        && rf_nnew > 0) {
                        if (rf_nv > vert_cap) {
                            fprintf(stderr, "grid_weld: vert capacity %zu "
                                    "exceeded by seam refine (%zu)\n",
                                    vert_cap, rf_nv);
                            RAISE(IO_Failed);
                        }
                        /* src[i] always indexes an already-provenanced vert
                         * (< out_nv + i), so one ordered pass suffices. */
                        for (size_t v2 = 0; v2 < rf_nnew; v2++)
                            vert_cube_idx[out_nv + v2] =
                                vert_cube_idx[rf_src[v2]];
                        out_verts = rf_v; out_nv = rf_nv;
                        flat_faces = rf_f; n_unique_faces = rf_nf;
                        color_nv = out_nv;  /* refine verts carry cube colors */
                        fprintf(stderr,
                            "Seam refine: %zu bnd + %zu int split(s), %zu "
                            "flip(s), %zu round(s) -> +%zu verts +%zu faces "
                            "(target %.1f vox, %zu plane(s))\n",
                            sst.bnd_splits, sst.int_splits, sst.flips,
                            sst.rounds, sst.verts_added, sst.faces_added,
                            (double)sp.target_len, fnp);
                    }
                }
                gw_phase("refine");
                dump_stage(arena, dump_stages_dir, stage_prefix, 1, "refine",
                           out_verts, out_nv, flat_faces, n_unique_faces,
                           vert_cube_idx, cube_palette, cubes.n, color_nv);
            }

            fprintf(stderr,
                "Seam-weld: cube=%.0f rho=%.2f band=%.2f over %zu faces\n",
                (double)seam_cube, (double)rho, (double)band,
                n_unique_faces);

            /* Phase-2 sheet labeling: connected components of the PRE-weld
             * mesh (cubes concatenated, no cross-cube faces yet), done before
             * the bridge so vert indices stay stable through it. Skipped when
             * the phase-2 reweld is disabled. */
            /* DEFAULT OFF (2026-07-09): the 4x5x5 audit found the phase-2 reweld
             * net-catastrophic -- it uses a rho=6 ball (2*rho=12 vox) with the
             * winding gate OFF, breaking the phase-1 invariant 2*BRIDGE_RHO_MAX=6
             * < 7-vox clearance. In the crumpled core a matched sheet folds
             * within 12 vox of itself, so the over-wide gate-off ball welds the
             * sheet to ITSELF across a fold -> overlapping flaps (and boundary
             * tangles -> MORE gaps). The two-sheet cloud restriction only guards
             * against a THIRD sheet, not intra-sheet fold-welds. Grid-wide,
             * disabling it drops gaps 515->115, flaps 1436->84, same_dir 39->4.
             * Opt back in with SEAM_REWELD=1 (e.g. to test a gated rewrite). */
            int reweld_off = (getenv("SEAM_REWELD") == NULL);
            int32_t *vert_sheet = NULL; size_t n_sheets = 0;
            if (!reweld_off)
                SheetReweld_label(arena, flat_faces, n_unique_faces, out_nv,
                                  &vert_sheet, &n_sheets);

            int32_t *sw_faces = NULL; size_t sw_nf = 0, n_bridge = 0;
            SeamWeld_bridge(arena, out_verts, out_nv, flat_faces,
                            n_unique_faces, seam_cube, rho, 0.0f, band,
                            NULL /* want_mask: primary weld = all */,
                            &seam_gate,
                            &sw_faces, &sw_nf, &n_bridge);
            fprintf(stderr, "  bridge: %zu -> %zu faces (+%zu bridge)\n",
                    n_unique_faces, sw_nf, n_bridge);
            gw_phase("bridge");

            /* Phase 2: sheet-correspondence permissive re-weld. Recover legit
             * same-sheet closures the conservative phase-1 winding gate left as
             * holes: confirm a 1:1 sheet correspondence across each seam
             * (phase-1 bridge votes + near-seam geometric overlap, mutual-best
             * + margin) and re-weld each confirmed pair on a cloud restricted
             * to those two sheets (gate OFF, wider ball). The safety is the
             * cloud restriction -- the ball cannot reach a third sheet.
             * SEAM_NO_REWELD -> phase-1-only baseline. */
            size_t n_bridge_total = n_bridge;
            if (!reweld_off && vert_sheet && n_sheets >= 2) {
                SheetReweldParams rwp;
                SheetReweld_default_params(&rwp);
                rwp.band = band;   /* match the phase-1 seam band */
                if ((e = getenv("SHEET_REWELD_RHOMAX"))) {
                    float v = (float)atof(e); if (v > 0.0f) rwp.rho_max = v;
                }
                int32_t *rw_faces = NULL; size_t rw_nf = 0;
                SheetReweldStats rwst;
                SheetReweld_process(arena, out_verts, out_nv, sw_faces, sw_nf,
                                    n_bridge, vert_sheet, n_sheets,
                                    seam_cube, &seam_gate, &rwp,
                                    &rw_faces, &rw_nf, &rwst);
                fprintf(stderr, "  reweld(phase2): %zu sheets, %zu candidate "
                        "pair(s), %zu confirmed, +%zu faces -> %zu\n",
                        rwst.n_sheets, rwst.n_candidates, rwst.n_pairs,
                        rwst.n_faces_added, rw_nf);
                sw_faces = rw_faces; sw_nf = rw_nf;
                n_bridge_total = n_bridge + rwst.n_faces_added;
            }
            if (pair_mode && n_bridge == 0)
                fprintf(stderr,
                    "WARNING: 0 bridge faces -- detect_planes found no shared "
                    "seam plane between the two components. Check the bboxes "
                    "above: they must straddle a multiple of %.0f on one axis "
                    "(see the front_edges.obj dump for the boundary edges).\n",
                    (double)seam_cube);

            /* Mark the verts touched by the appended bridge faces as weld
             * verts (the BPA model's notion of a weld is the bridge, not a
             * fused vertex). These indices are < color_nv and survive the
             * fold/pinhole/cull below, so the green seam highlight + the
             * weld-vert report stat below stay meaningful. */
            for (size_t f = sw_nf - n_bridge_total; f < sw_nf; f++) {
                for (int k = 0; k < 3; k++) {
                    int32_t vv = sw_faces[f * 3 + k];
                    if (vv >= 0 && (size_t)vv < out_nv) vert_is_weld[vv] = 1;
                }
            }

            /* Stage 01: raw BPA bridge output (pre-orientation). */
            dump_stage(arena, dump_stages_dir, stage_prefix, 2, "bridge",
                       out_verts, out_nv, sw_faces, sw_nf,
                       vert_cube_idx, cube_palette, cubes.n, color_nv);

            /* Orient the fused mesh BEFORE fold cleanup with the live-winding
             * BFS (OrientMesh_consistent), not the stale-edge repair_winding:
             * the latter's cached edge directions go stale once a parent face
             * flips, leaving residual same_dir that FoldCleanup then misreads as
             * folds and deletes valid faces. The directed-glue bridge is already
             * consistently wound, so this pass only reconciles the two cubes'
             * independent per-cube orientations across the seam. The outer
             * OrientMesh below re-runs (idempotent) after pinhole. */
            size_t sw_comp = 0, sw_flipped = 0, sw_resid = 0;
            OrientMesh_consistent(arena, out_verts, out_nv, NULL,
                                  sw_faces, sw_nf,
                                  &sw_flipped, &sw_comp, &sw_resid);
            fprintf(stderr, "  pre-orient: %zu components, %zu flipped, "
                    "%zu residual same_dir\n", sw_comp, sw_flipped, sw_resid);
            gw_phase("preorient");
            /* Stage 02: after pre-fold orientation reconciliation. */
            dump_stage(arena, dump_stages_dir, stage_prefix, 3, "preorient",
                       out_verts, out_nv, sw_faces, sw_nf,
                       vert_cube_idx, cube_palette, cubes.n, color_nv);

            ComponentMesh cm;
            memset(&cm, 0, sizeof cm);
            cm.verts = out_verts; cm.faces = sw_faces;
            cm.nv = out_nv; cm.nf = sw_nf; cm.comp_id = 1; cm.self = &cm;

            /* Cut BPA fold flaps (same_dir interior edges BFS cannot fix:
             * dot<0 creases, plus clear dot>=0 low-coherence outliers). Leaves
             * a boundary notch the pinhole pass below recloses. Must precede
             * PinholeFill per fold_cleanup.h. */
            size_t fold_removed = 0;
            FoldCleanup_process(arena, &cm, 1, 8, &fold_removed);
            fprintf(stderr, "  foldcleanup: removed %zu fold faces -> %zu faces\n",
                    fold_removed, cm.nf);
            gw_phase("foldclean");
            /* Stage 03: after fold-flap removal. */
            dump_stage(arena, dump_stages_dir, stage_prefix, 4, "foldcleanup",
                       cm.verts, cm.nv, cm.faces, cm.nf,
                       vert_cube_idx, cube_palette, cubes.n, color_nv);

            size_t splits = 0, filled = 0, added = 0, skipped = 0;
            /* SPLIT-ONLY here. The bridge leaves bowtie VERTICES but no
             * non-manifold edges; split_pinch_verts resolves them cleanly
             * (duplicate per fan, monotone -- never makes a >2-face edge). The
             * old call also FILLED 3-loops, and that fill re-closed the boundary
             * the split had just opened, minting the seam NM edges (889 bowties
             * -> 38 NM edges). Tiny-hole filling is left to the outer pinhole
             * pass (line ~1300), which now runs on a vertex-manifold mesh.
             * --no-pinhole / SEAM_NO_PINHOLE skips it (bridge-only diagnostic). */
            if (!no_pinhole)
            PinholeFill_split_pinches(arena, &cm, 1, &splits);
            fprintf(stderr,
                "  pinhole: splits=%zu filled=%zu tris+=%zu skipped=%zu"
                " -> %zu faces, %zu verts\n",
                splits, filled, added, skipped, cm.nf, cm.nv);
            gw_phase("pinchsplit");
            /* Stage 04: after pinhole micro-hole reclose. */
            dump_stage(arena, dump_stages_dir, stage_prefix, 5, "pinhole",
                       cm.verts, cm.nv, cm.faces, cm.nf,
                       vert_cube_idx, cube_palette, cubes.n, color_nv);

            /* Cull the tiny floating slivers the re-bridge sheds at the seam
             * (1-2 triangle fragments are not part of any sheet). */
            size_t culled_nf = cm.nf, culled_comps = 0;
            size_t culled = cull_tiny_components(arena, cm.faces, cm.nf, 16,
                                                 &culled_nf, &culled_comps);
            cm.nf = culled_nf;
            fprintf(stderr, "  cull: removed %zu faces in %zu tiny comps -> %zu faces\n",
                    culled, culled_comps, cm.nf);
            gw_phase("cull");
            /* Stage 05: after tiny-component cull. */
            dump_stage(arena, dump_stages_dir, stage_prefix, 6, "cull",
                       cm.verts, cm.nv, cm.faces, cm.nf,
                       vert_cube_idx, cube_palette, cubes.n, color_nv);

            /* Post-weld sliver / T-junction cleanup: Surazhsky-Gotsman flips
             * first, then a guarded short-edge collapse for the residue. This
             * is the only edge-collapse the bridge faces ever see (per-cube QEM
             * ran before the weld). Vertices are not moved; collapsed verts are
             * orphaned so the color arrays below stay valid by index. */
            if (!no_cleanup) {
                WeldCleanupParams wcp;
                WeldCleanupStats wcs;
                WeldCleanup_default_params(&wcp);
                WeldCleanup_process(arena, &cm, &wcp, &wcs);
                fprintf(stderr,
                    "  cleanup: %zu flips, %zu collapses -> %zu faces"
                    " (sliver/degenerate targets %zu -> %zu)\n",
                    wcs.n_flips, wcs.n_collapses, cm.nf,
                    wcs.targets_in, wcs.targets_out);
            }
            gw_phase("cleanup");

            /* Stage 06: after post-weld flip + short-edge collapse cleanup. */
            dump_stage(arena, dump_stages_dir, stage_prefix, 7, "cleanup",
                       cm.verts, cm.nv, cm.faces, cm.nf,
                       vert_cube_idx, cube_palette, cubes.n, color_nv);

            /* Capstone of the weld cycle: detect + fill INTERIOR holes the join
             * created. A boundary bay the weld could close has, by being
             * bridged on its open side, BECOME an interior hole; the geometric
             * (signed-area) interior test fills exactly those and leaves the
             * outer perimeter + still-open bays alone. CDT/Liepa fill. */
            if (!no_holefill) {
                size_t hl_loops = 0, hl_interior = 0, hl_filled = 0;
                HoleFill_process_ex(arena, &cm.verts, &cm.faces, &cm.nv, &cm.nf,
                                    NULL, 1 /* interior_only */,
                                    &hl_loops, &hl_interior, &hl_filled);
                fprintf(stderr,
                    "  holefill: %zu loops, %zu interior -> %zu filled"
                    " -> %zu faces, %zu verts\n",
                    hl_loops, hl_interior, hl_filled, cm.nf, cm.nv);
            }
            gw_phase("holefill");

            /* Stage 07: after interior-hole fill (the stage that hangs on the
             * full grid -- per-hole dumps in <dir>/holes/ catch the culprit). */
            dump_stage(arena, dump_stages_dir, stage_prefix, 8, "holefill",
                       cm.verts, cm.nv, cm.faces, cm.nf,
                       vert_cube_idx, cube_palette, cubes.n, color_nv);

            flat_faces = cm.faces;
            n_unique_faces = cm.nf;
            out_verts = cm.verts;   /* pinhole may have grown the vert array */
            out_nv = cm.nv;
        }

        /* Pre-repair audit. */
        ms = manifold_audit(arena, flat_faces, n_unique_faces, NULL, NULL);
        {
            MeshManifoldStats vms0 = MeshManifold_audit(arena, out_nv,
                                                        flat_faces, n_unique_faces);
            fprintf(stderr,
                "Pre-repair audit: unpaired=%zu non_manifold=%zu same_dir=%zu "
                "manifold=%zu pinch_verts=%zu\n",
                ms.unpaired, ms.non_manifold, ms.same_dir_pairs,
                ms.manifold_pairs, vms0.nm_verts);
        }
        gw_phase("audit_pre");

        /* Final winding pass: OrientMesh_consistent's live-winding BFS drives
         * same_dir to 0 even across the seam (repair_winding's stale-edge cache
         * could not). residual_same_dir > 0 here flags a genuine non-orientable
         * knot (cleared upstream by pinhole splitting), not a welder defect. */
        size_t n_components = 0, n_flipped = 0, n_resid = 0;
        OrientMesh_consistent(arena, out_verts, out_nv, NULL,
                              flat_faces, n_unique_faces,
                              &n_flipped, &n_components, &n_resid);
        fprintf(stderr,
            "Winding repair (OrientMesh): %zu components, %zu flipped, "
            "%zu residual same_dir\n", n_components, n_flipped, n_resid);
        gw_phase("orient");

        /* Post-weld cross-component orientation: OrientMesh makes each component
         * internally consistent but anchors its global sign against winding-
         * derived normals, so a backward-wound block (consistent with its OWN
         * flipped normals) survives. Flip whole components whose normals oppose
         * their SPATIAL neighbours instead. Radius 3 vox < the >=7-vox inter-wrap
         * clearance, so a component only ever votes against the same sheet across
         * the seam gap, never an adjacent wrap.
         *
         * DETACHED components (no other geometry within 3 vox anywhere -- shell
         * fragments whose connection was never meshed) get no spatial vote and
         * would keep the per-cube (1,1,1) anchor's azimuth-dependent sign. When
         * the umbilicus env is set (same SEAM_UMBILICUS_Y/X as the seam gate),
         * those fall back to a RADIAL vote against the largest component. */
        size_t ow_flipped = 0, ow_radial = 0;
        float ow_axp[3] = { 0.0f, 0.0f, 0.0f };
        float ow_axd[3] = { 1.0f, 0.0f, 0.0f };   /* scroll axis = Z in (z,y,x) */
        const float *ow_paxp = NULL, *ow_paxd = NULL;
        {
            const char *ey = getenv("SEAM_UMBILICUS_Y");
            const char *ex = getenv("SEAM_UMBILICUS_X");
            if (ey != NULL && ex != NULL) {
                ow_axp[1] = (float)atof(ey);
                ow_axp[2] = (float)atof(ex);
                ow_paxp = ow_axp;
                ow_paxd = ow_axd;
            }
        }
        OrientWeld_components_axis(arena, out_verts, out_nv, flat_faces,
                                   n_unique_faces, 3.0f, ow_paxp, ow_paxd,
                                   &ow_flipped, &ow_radial);
        fprintf(stderr, "Component orientation: %zu backward component(s) flipped"
                " (%zu decided by radial fallback%s)\n", ow_flipped, ow_radial,
                ow_paxp != NULL ? "" : " -- fallback OFF, no umbilicus env");
        gw_phase("orientweld");

        /* Final pinhole pass over the WHOLE welded mesh. The earlier pass (inside
         * the seam block) ran before cull / cleanup-collapse / holefill / the
         * orientation pass -- each of which can leave or expose a single-triangle
         * hole ANYWHERE on the mesh, not just at the seam. PinholeFill has no
         * seam-zone restriction (its only gates are the no-merger diameter, the
         * no-bubble component, and the degenerate-sliver guards), so re-running it
         * on the whole mesh fills the pinholes those later stages left, wherever
         * they are. --no-pinhole / SEAM_NO_PINHOLE skips it. */
        if (!no_pinhole) {
            ComponentMesh pcm; memset(&pcm, 0, sizeof pcm);
            pcm.verts = out_verts; pcm.faces = flat_faces;
            pcm.nv = out_nv; pcm.nf = n_unique_faces; pcm.comp_id = 1; pcm.self = &pcm;
            size_t f_sp = 0, f_fl = 0, f_ad = 0, f_sk = 0;
            PinholeFill_process(arena, &pcm, 1, 0, &f_sp, &f_fl, &f_ad, &f_sk);
            fprintf(stderr,
                "Final pinhole pass (whole mesh): filled=%zu tris+=%zu skipped=%zu"
                " -> %zu faces, %zu verts\n", f_fl, f_ad, f_sk, pcm.nf, pcm.nv);
            flat_faces = pcm.faces; n_unique_faces = pcm.nf;
            out_verts = pcm.verts; out_nv = pcm.nv;
        }
        gw_phase("pinhole2");

        /* Micro-weld: post-concat stages (seam bridge, CDT fill) can CREATE a
         * vertex exactly on an existing vertex's position (observed: 4 such
         * pairs on a 2-cube weld, each pinching a small boundary loop into a
         * figure-8 slit that nothing can close -- the CDT path drops its < 4-
         * vert pinch sub-loops, and PinholeFill sees no 3-cycle because the
         * coincident pair has two distinct indices). Welding the pair zips the
         * slit shut with no new faces. eps is far below any legitimate vertex
         * spacing; the concat-time dedup only catches input duplicates.
         * SEAM_NO_MICROWELD=1 disables. */
        if (!getenv("SEAM_NO_MICROWELD")) {
            size_t w_nf = 0, w_nv = 0;
            float *w_verts = NULL;
            Weld_verts(arena, out_verts, out_nv, NULL,
                       flat_faces, n_unique_faces, &w_nf,
                       1e-3f, /*guard_orient=*/true, &w_verts, &w_nv, NULL);
            if (w_nv != out_nv || w_nf != n_unique_faces) {
                fprintf(stderr, "Micro-weld: %zu -> %zu verts, %zu -> %zu faces "
                        "(coincident stage-created duplicates)\n",
                        out_nv, w_nv, n_unique_faces, w_nf);
            }
            out_verts = w_verts; out_nv = w_nv; n_unique_faces = w_nf;
        }
        gw_phase("microweld");

        /* Fill fixpoint. Each fill pass reshapes the boundary structure (a
         * pinhole fill splits or shrinks a larger loop; an interior fill
         * exposes fresh 3-loops), and a single pass leaves those orphans open
         * -- observed as 4-6-edge slots at a grazing seam surviving both the
         * interior fill AND the final pinhole pass. Re-run interior-fill +
         * pinhole until neither makes progress (bounded). */
        if (!no_holefill && !no_pinhole) {
            for (int round = 1; round <= 3; round++) {
                size_t r_loops = 0, r_int = 0, r_filled = 0;
                HoleFill_process_ex(arena, &out_verts, &flat_faces,
                                    &out_nv, &n_unique_faces,
                                    NULL, 1 /* interior_only */,
                                    &r_loops, &r_int, &r_filled);
                ComponentMesh rcm; memset(&rcm, 0, sizeof rcm);
                rcm.verts = out_verts; rcm.faces = flat_faces;
                rcm.nv = out_nv; rcm.nf = n_unique_faces;
                rcm.comp_id = 1; rcm.self = &rcm;
                size_t r_sp = 0, r_fl = 0, r_ad = 0, r_sk = 0;
                PinholeFill_process(arena, &rcm, 1, 0, &r_sp, &r_fl, &r_ad, &r_sk);
                out_verts = rcm.verts; flat_faces = rcm.faces;
                out_nv = rcm.nv; n_unique_faces = rcm.nf;
                if (r_filled + r_fl == 0) break;
                fprintf(stderr, "Fill fixpoint round %d: interior=%zu pinhole=%zu"
                        " -> %zu faces, %zu verts\n",
                        round, r_filled, r_fl, n_unique_faces, out_nv);
            }
        }
        gw_phase("fixpoint");

        /* Lambda gate -- "don't weld if it produces high lambda". The seam bridge
         * can fuse two DIFFERENT wraps where they pass within ~rho at the core; the
         * bridge faces then carry high Crane energy lambda (a crease a single
         * developable wrap never makes). Sever them: drop faces touching a high-
         * lambda vert that sits in the cross-cube SEAM zone, reopening the wrong
         * merger. A correct same-wrap weld is developable (lambda ~ 0) and is kept.
         * OFF by default (2026-07-08 A/B: weak on real mergers -- cut 2 of 9 handles
         * as an add-on, left 7 of 8 alone -- and a same-wrap FOLD crossing a seam
         * also carries high lambda, so it risks intra-sheet splits, the worse
         * failure). Enable with SEAM_LAMBDA_GATE=1 for diagnostics/experiments;
         * SEAM_LAMBDA_MAX / SEAM_LAMBDA_ZONE tune. */
        if (getenv("SEAM_LAMBDA_GATE")) {
            double lmax = 0.05f; float zone = 4.0f;
            { const char *e = getenv("SEAM_LAMBDA_MAX");  if (e) { double v=atof(e); if (v>0) lmax=v; } }
            { const char *e = getenv("SEAM_LAMBDA_ZONE"); if (e) { double v=atof(e); if (v>0) zone=(float)v; } }
            double *lam = (double *)ARENA_ALLOC(arena, (long)(out_nv*sizeof(double)));
            if (Develop_vertex_energy(arena, out_verts, out_nv, flat_faces,
                                      n_unique_faces, lam) == 0) {
                unsigned char *cut = (unsigned char *)ARENA_CALLOC(arena, (long)out_nv, 1L);
                size_t ncut = 0;
                for (size_t v = 0; v < out_nv; v++)
                    if (lam[v] > lmax && near_cube_boundary(&out_verts[v*3], zone)) {
                        cut[v] = 1; ncut++;
                    }
                size_t w = 0, removed = 0;
                for (size_t f = 0; f < n_unique_faces; f++) {
                    int32_t a = flat_faces[f*3+0], b = flat_faces[f*3+1], c = flat_faces[f*3+2];
                    if (cut[a] || cut[b] || cut[c]) { removed++; continue; }
                    flat_faces[w*3+0]=a; flat_faces[w*3+1]=b; flat_faces[w*3+2]=c; w++;
                }
                n_unique_faces = w;
                fprintf(stderr,
                    "Lambda gate: %zu high-lambda seam vert(s) -> %zu weld face(s) "
                    "severed (lambda > %.3f within %.1f vox of a seam)\n",
                    ncut, removed, lmax, (double)zone);
            }
        }

        /* Phase-jump sever (fusion-line cutter EXPERIMENT -- default off,
         * SEAM_PHASE_SEVER=1 arms it; needs the gate's umbilicus/pitch env).
         * Target: junctions the m7 PREDICTION itself contains (delamination
         * pairs ~half a pitch apart, the crushed-core web) -- no bridge gate
         * can touch them; the mesh must be CUT along the contact line.
         *
         * STATUS 2026-07-09, measured on the red/pink corner puddle: local
         * thresholds DO NOT WORK. v1 (edge |dw| > tol AND lambda) and v2
         * (below: clusters of faces whose own phase span > tol, small +
         * lambda-hot) both left the junction connected -- the membrane
         * crosses the 0.5-turn gap through dozens of ~0.025-turn micro-steps
         * (max face span there: 0.211; median crossing-face span 0.025), so
         * any per-face/per-edge threshold either misses it or shreds
         * innocent geometry (659 faces > 0.10 span in that one box alone).
         * What DOES work (proven offline on the exemplar): PLATEAU
         * MEMBERSHIP -- the two surfaces are local phase plateaus (w = 28.05
         * and 28.55 there); cutting the 679 faces whose mid-phase lies
         * BETWEEN the plateaus disconnects the layers surgically. v3 =
         * per-vertex plateau assignment from neighborhood phase modes, cut
         * between-plateau fabric; folds (one shared plateau), grazing sheets
         * (single mode), and the exempt core are safe by construction.
         * The v2 code below is retained as the experiment scaffold. Knobs:
         * SEAM_PHASE_SEVER_TOL (0.25 turn), SEAM_PHASE_SEVER_LAMBDA (0.05),
         * SEAM_PHASE_SEVER_RMIN_PITCHES (2.0), SEAM_PHASE_SEVER_MAXC (2500). */
        if (getenv("SEAM_PHASE_SEVER")) {
            double ps_umb_y = 0.0, ps_umb_x = 0.0, ps_pitch = 0.0;
            { const char *e = getenv("SEAM_UMBILICUS_Y"); if (e) ps_umb_y = atof(e); }
            { const char *e = getenv("SEAM_UMBILICUS_X"); if (e) ps_umb_x = atof(e); }
            { const char *e = getenv("SEAM_WRAP_PITCH"); if (e) { double v = atof(e); if (v > 0) ps_pitch = v; } }
            double ps_tol = 0.25, ps_lam = 0.05, ps_rmin_p = 2.0;
            { const char *e = getenv("SEAM_PHASE_SEVER_TOL"); if (e) { double v = atof(e); if (v > 0) ps_tol = v; } }
            { const char *e = getenv("SEAM_PHASE_SEVER_LAMBDA"); if (e) { double v = atof(e); if (v > 0) ps_lam = v; } }
            { const char *e = getenv("SEAM_PHASE_SEVER_RMIN_PITCHES"); if (e) { double v = atof(e); if (v > 0) ps_rmin_p = v; } }
            if (ps_pitch <= 0.0 || (ps_umb_y == 0.0 && ps_umb_x == 0.0)) {
                fprintf(stderr, "Phase sever: SKIPPED (needs SEAM_UMBILICUS_Y/X "
                        "+ SEAM_WRAP_PITCH)\n");
            } else {
                double *lam = (double *)ARENA_ALLOC(arena,
                                  (long)(out_nv * sizeof(double)));
                double *ww = (double *)ARENA_ALLOC(arena,
                                  (long)(out_nv * sizeof(double)));
                double *rr = (double *)ARENA_ALLOC(arena,
                                  (long)(out_nv * sizeof(double)));
                if (Develop_vertex_energy(arena, out_verts, out_nv, flat_faces,
                                          n_unique_faces, lam) == 0) {
                    double rmin = ps_rmin_p * ps_pitch;
                    size_t ps_maxc = 2500;   /* cluster-size cap: membranes are
                                              * compact strips; a huge phase-
                                              * mixing region is real geometry */
                    { const char *e = getenv("SEAM_PHASE_SEVER_MAXC");
                      if (e) { long v = atol(e); if (v > 0) ps_maxc = (size_t)v; } }
                    for (size_t v = 0; v < out_nv; v++) {
                        double dy = (double)out_verts[v*3+1] - ps_umb_y;
                        double dx = (double)out_verts[v*3+2] - ps_umb_x;
                        rr[v] = hypot(dy, dx);
                        ww[v] = rr[v] / ps_pitch - atan2(dy, dx) / (2.0 * M_PI);
                    }
                    /* Candidate faces: own phase span > tol (the crossing strip
                     * between two surfaces), fully outside the core exemption.
                     * The span uses the same wrapped-dtheta phase as the gate:
                     * evaluate all 3 edges, take the max |dw|. */
                    uint8_t *cand = (uint8_t *)ARENA_CALLOC(arena,
                                        (long)n_unique_faces, 1L);
                    for (size_t f = 0; f < n_unique_faces; f++) {
                        int32_t t[3] = { flat_faces[f*3+0], flat_faces[f*3+1],
                                         flat_faces[f*3+2] };
                        if (rr[t[0]] < rmin || rr[t[1]] < rmin || rr[t[2]] < rmin)
                            continue;
                        double span = 0.0;
                        for (int e = 0; e < 3; e++) {
                            int32_t a = t[e], b = t[(e+1)%3];
                            double dth = atan2((double)out_verts[a*3+1] - ps_umb_y,
                                               (double)out_verts[a*3+2] - ps_umb_x)
                                       - atan2((double)out_verts[b*3+1] - ps_umb_y,
                                               (double)out_verts[b*3+2] - ps_umb_x);
                            while (dth >  M_PI) dth -= 2.0*M_PI;
                            while (dth < -M_PI) dth += 2.0*M_PI;
                            double dw = fabs((rr[a] - rr[b]) / ps_pitch
                                             - dth / (2.0 * M_PI));
                            if (dw > span) span = dw;
                        }
                        if (span > ps_tol) cand[f] = 1;
                    }
                    /* Cluster candidates via shared verts (union-find over
                     * faces); cut a cluster iff it is SMALL (<= maxc faces)
                     * and its vertex set touches high lambda -- the attach
                     * lines of a junction membrane are non-developable. */
                    int32_t *fpar = (int32_t *)ARENA_ALLOC(arena,
                                        (long)(n_unique_faces * sizeof(int32_t)));
                    for (size_t f = 0; f < n_unique_faces; f++)
                        fpar[f] = (int32_t)f;
                    /* map vert -> one candidate face, to union share-a-vert faces */
                    int32_t *vf = (int32_t *)ARENA_ALLOC(arena,
                                      (long)(out_nv * sizeof(int32_t)));
                    for (size_t v = 0; v < out_nv; v++) vf[v] = -1;
                    for (size_t f = 0; f < n_unique_faces; f++) {
                        if (!cand[f]) continue;
                        for (int k = 0; k < 3; k++) {
                            int32_t v = flat_faces[f*3+k];
                            if (vf[v] < 0) { vf[v] = (int32_t)f; continue; }
                            /* union f with vf[v] */
                            int32_t x = (int32_t)f, y = vf[v];
                            while (fpar[x] != x) x = fpar[x] = fpar[fpar[x]];
                            while (fpar[y] != y) y = fpar[y] = fpar[fpar[y]];
                            if (x != y) fpar[x] = y;
                        }
                    }
                    /* cluster stats */
                    size_t removed = 0, n_clusters = 0, n_cut_clusters = 0;
                    /* count sizes + lambda touch per root (two passes) */
                    int32_t *croot = (int32_t *)ARENA_ALLOC(arena,
                                         (long)(n_unique_faces * sizeof(int32_t)));
                    size_t *csize = (size_t *)ARENA_CALLOC(arena,
                                        (long)n_unique_faces,
                                        (long)sizeof(size_t));
                    uint8_t *chot = (uint8_t *)ARENA_CALLOC(arena,
                                        (long)n_unique_faces, 1L);
                    for (size_t f = 0; f < n_unique_faces; f++) {
                        if (!cand[f]) { croot[f] = -1; continue; }
                        int32_t x = (int32_t)f;
                        while (fpar[x] != x) x = fpar[x] = fpar[fpar[x]];
                        croot[f] = x;
                        if (csize[x]++ == 0) n_clusters++;
                        for (int k = 0; k < 3; k++)
                            if (lam[flat_faces[f*3+k]] > ps_lam) chot[x] = 1;
                    }
                    size_t w = 0;
                    for (size_t f = 0; f < n_unique_faces; f++) {
                        int cutf = 0;
                        if (croot[f] >= 0) {
                            size_t sz = csize[croot[f]];
                            if (sz <= ps_maxc && chot[croot[f]]) cutf = 1;
                        }
                        if (cutf) { removed++; continue; }
                        flat_faces[w*3+0] = flat_faces[f*3+0];
                        flat_faces[w*3+1] = flat_faces[f*3+1];
                        flat_faces[w*3+2] = flat_faces[f*3+2];
                        w++;
                    }
                    for (size_t f = 0; f < n_unique_faces; f++)
                        if (croot[f] >= 0 && csize[croot[f]] <= ps_maxc
                            && chot[croot[f]] && croot[f] == (int32_t)f)
                            n_cut_clusters++;
                    n_unique_faces = w;
                    fprintf(stderr,
                        "Phase sever: %zu face(s) in %zu cluster(s) cut "
                        "(of %zu candidate clusters; span > %.2f turn, "
                        "cluster <= %zu faces, lambda > %.3f, r > %.0f)\n",
                        removed, n_cut_clusters, n_clusters, ps_tol,
                        ps_maxc, ps_lam, rmin);
                }
            }
        }

        /* Seam-hole fill (audit output/weld_audit_4x5x5/AUDIT.md): the BPA bridge
         * zips most of each grazing seam but leaves a dotted line of tiny
         * straddling punctures; pinhole's 4.5-vox diameter gate and interior-only
         * hole-fill leave many open. Close them here, AFTER every other closer,
         * with three merger-safe gates: near-seam + straddle, small extent, and
         * winding-phase coherence (armed by the same umbilicus/pitch as the seam
         * gate -- so it can never fill across the sub-clearance core wraps). The
         * final manifold guard below re-verifies the result. SEAM_NO_SEAMFILL
         * disables. */
        if (!getenv("SEAM_NO_SEAMFILL")) {
            ComponentMesh scm; memset(&scm, 0, sizeof scm);
            scm.verts = out_verts; scm.faces = flat_faces;
            scm.nv = out_nv; scm.nf = n_unique_faces; scm.comp_id = 1; scm.self = &scm;
            SeamHoleFillParams shp; SeamHoleFill_default_params(&shp);
            shp.cube = seam_cube;
            const char *e;
            if ((e = getenv("SEAM_UMBILICUS_Y"))) shp.umb_y = atof(e);
            if ((e = getenv("SEAM_UMBILICUS_X"))) shp.umb_x = atof(e);
            if ((e = getenv("SEAM_WRAP_PITCH"))) { double v = atof(e); if (v > 0) shp.pitch = v; }
            if ((e = getenv("SEAM_WIND_TOL"))) { double v = atof(e); if (v > 0) shp.wind_tol_turns = v; }
            if ((e = getenv("SEAM_SEAMFILL_MAXLOOP"))) { int v = atoi(e); if (v >= 3) shp.max_loop = v; }
            if ((e = getenv("SEAM_SEAMFILL_EXTENT"))) { double v = atof(e); if (v > 0) shp.max_extent = v; }
            SeamHoleFillStats shs;
            SeamHoleFill_process(arena, &scm, &shp, &shs);
            fprintf(stderr,
                "Seam-hole fill: %zu candidate(s), %zu filled (+%zu tris); "
                "skip[phase=%zu extent=%zu geom=%zu manifold=%zu bubble=%zu]%s\n",
                shs.seam_candidates, shs.filled, shs.tris_added,
                shs.skip_phase, shs.skip_extent, shs.skip_geom,
                shs.skip_manifold, shs.skip_bubble,
                shp.pitch > 0.0 ? "" : "  (phase gate OFF -- no umbilicus/pitch)");
            flat_faces = scm.faces; n_unique_faces = scm.nf;
            out_verts = scm.verts; out_nv = scm.nv;
        }
        gw_phase("seamfill");

        /* Final manifold guard. The strict seam bridge leaves only a handful of
         * residual non-manifold edges (vs many under the old relaxed zip).
         * Resolve >2-face edges + split any residual pinch so the welded mesh is
         * a 2-manifold by construction -- the same guard the per-cube path runs
         * after trim. reorient=0: the weld's own OrientMesh/OrientWeld already
         * fixed winding and anchored component signs to spatial neighbours;
         * resolve+split preserve winding, so re-orienting here would only risk
         * re-flipping component signs with no manifold benefit. */
        if (!no_cleanup) {
            ComponentMesh wcm; memset(&wcm, 0, sizeof wcm);
            wcm.verts = out_verts; wcm.faces = flat_faces;
            wcm.nv = out_nv; wcm.nf = n_unique_faces; wcm.comp_id = 1; wcm.self = &wcm;
            ManifoldGuardStats mg;
            ManifoldGuard_process(arena, &wcm, 1, 0 /*reorient*/, &mg);
            fprintf(stderr,
                "Manifold guard: %zu NM-edge(s) (-%zu faces), %zu pinch split(s)\n",
                mg.nm_edges_resolved, mg.faces_deleted, mg.pinch_splits);
            flat_faces = wcm.faces; n_unique_faces = wcm.nf;
            out_verts = wcm.verts; out_nv = wcm.nv;
        }
        gw_phase("guard");

        /* Post-guard convergence. The manifold guard drops faces to resolve NM
         * edges, which OPENS small boundary notches that every earlier closer
         * (all of which ran before it) never sees -- observed as residual len-3/4
         * seam gaps in the 4x5x5 audit. Re-run the merger-safe closers (pinhole
         * 3-loop + bowtie, seam-hole ear-clip) and the guard until none makes
         * progress (bounded 3 rounds). Every closer is already merger-safe
         * (phase-gated / a 3-loop adds no new edge); the guard runs LAST each
         * round so the mesh stays a 2-manifold by construction. SEAM_NO_POSTGUARD
         * disables. */
        if (!no_cleanup && !getenv("SEAM_NO_POSTGUARD")) {
            SeamHoleFillParams pgp; SeamHoleFill_default_params(&pgp);
            pgp.cube = seam_cube;
            { const char *e;
              if ((e = getenv("SEAM_UMBILICUS_Y"))) pgp.umb_y = atof(e);
              if ((e = getenv("SEAM_UMBILICUS_X"))) pgp.umb_x = atof(e);
              if ((e = getenv("SEAM_WRAP_PITCH")))  { double v = atof(e); if (v > 0) pgp.pitch = v; }
              if ((e = getenv("SEAM_WIND_TOL")))    { double v = atof(e); if (v > 0) pgp.wind_tol_turns = v; }
              if ((e = getenv("SEAM_SEAMFILL_EXTENT"))) { double v = atof(e); if (v > 0) pgp.max_extent = v; } }
            /* 3 rounds is the empirical fixed point. Rounds 4-6 were tried
             * (2026-07-18) and enter a LIMIT CYCLE at the residual fold
             * slits: seamfill caps the slit, the guard's pinch audit splits
             * it back open, repeat -- net faces oscillate +-2 with no
             * convergence. The residue (11-12 sub-6-vox loops on the 4x5x5)
             * is loops containing bit-coincident vertex PAIRS at folds,
             * where the orient-guarded micro-weld refuses the merge
             * (recto-verso protection) and any triangulation would need a
             * zero-area triangle. Those stay open by design. */
            for (int gr = 1; gr <= 3; gr++) {
                /* Zip coincident fill-created vertex pairs FIRST each round:
                 * the late closers (CDT interior fill, pinhole, seam fill) can
                 * mint a vertex exactly on an existing one, leaving a
                 * zero-width slit no triangulator can cap (observed: a 5-loop
                 * with two bit-coincident verts surviving every closer). The
                 * main micro-weld stage ran BEFORE these verts existed. Same
                 * call + eps; SEAM_NO_MICROWELD disables both. */
                if (!getenv("SEAM_NO_MICROWELD")) {
                    size_t pw_nf = 0, pw_nv = 0;
                    float *pw_v = NULL;
                    Weld_verts(arena, out_verts, out_nv, NULL,
                               flat_faces, n_unique_faces, &pw_nf,
                               1e-3f, /*guard_orient=*/true, &pw_v, &pw_nv, NULL);
                    if (pw_nv != out_nv) {
                        fprintf(stderr, "  postguard micro-weld: %zu -> %zu verts\n",
                                out_nv, pw_nv);
                    }
                    out_verts = pw_v; out_nv = pw_nv; n_unique_faces = pw_nf;
                }
                ComponentMesh pg; memset(&pg, 0, sizeof pg);
                pg.verts = out_verts; pg.faces = flat_faces;
                pg.nv = out_nv; pg.nf = n_unique_faces; pg.comp_id = 1; pg.self = &pg;
                size_t g_sp=0, g_fl=0, g_ad=0, g_sk=0;
                PinholeFill_process(arena, &pg, 1, 0, &g_sp, &g_fl, &g_ad, &g_sk);
                SeamHoleFillStats pgs; SeamHoleFill_process(arena, &pg, &pgp, &pgs);
                ManifoldGuardStats pmg; ManifoldGuard_process(arena, &pg, 1, 0, &pmg);
                out_verts = pg.verts; flat_faces = pg.faces;
                out_nv = pg.nv; n_unique_faces = pg.nf;
                fprintf(stderr,
                    "Post-guard round %d: pinhole=%zu seamfill=%zu guard(-%zu f, "
                    "%zu split) -> %zu faces, %zu verts\n",
                    gr, g_fl, pgs.filled, pmg.faces_deleted, pmg.pinch_splits,
                    n_unique_faces, out_nv);
                if (g_fl == 0 && pgs.filled == 0 &&
                    pmg.faces_deleted == 0 && pmg.pinch_splits == 0) break;
            }
        }
        gw_phase("postguard");

        /* Seam-band recoarsen: collapse the temporarily-fine seam band (thin
         * graded CVT rim / weld-time refinement) back toward the coarse budget.
         * Placed AFTER every closer -- interior fill, fixpoint, seam-hole fill,
         * manifold guard, post-guard convergence -- because those close holes
         * at DENSE band scale, where every loop is small enough for the
         * merger-safe fill gates (pinhole 6.5 / seamfill extent caps). Running
         * recoarsen earlier coarsened residual slit rims past those gates and
         * left them open (measured: seam loops 9 -> 188 on the 4x5x5 with an
         * aggressive band). Here the seam is as closed as it will get, and
         * recoarsening cannot reopen it: collapses contract existing edges
         * only, boundary loops exit bit-identical (weld_cleanup.h) -- which
         * also leaves unbridged holes and a hierarchical level's outer faces
         * untouched, letting hierarchical_weld stack this per level. The final
         * winding repair + manifold audits below re-verify the result.
         * SEAM_NO_RECOARSEN=1 skips; SEAM_RECOARSEN_BELOW/_BAND tune. */
        if (!no_bridge && !no_cleanup && !getenv("SEAM_NO_RECOARSEN")) {
            /* Same seam band default the bridge used (its `band` local is
             * scoped to the bridge block above; re-read env, same default). */
            double r_band = 6.0;
            { const char *rbe = getenv("SEAM_BAND");
              if (rbe) { double v = atof(rbe); if (v > 0) r_band = v; } }
            uint8_t *r_used = (uint8_t *)ARENA_CALLOC(arena, (long)out_nv, 1L);
            for (size_t f = 0; f < n_unique_faces; f++) {
                r_used[flat_faces[f*3+0]] = 1;
                r_used[flat_faces[f*3+1]] = 1;
                r_used[flat_faces[f*3+2]] = 1;
            }
            SeamPlane rplanes[64];
            size_t rnp = SeamPlanes_detect(out_verts, out_nv, r_used,
                                           (double)seam_cube, r_band,
                                           rplanes, 64);
            if (rnp > 0) {
                ComponentMesh rcm; memset(&rcm, 0, sizeof rcm);
                rcm.verts = out_verts; rcm.faces = flat_faces;
                rcm.nv = out_nv; rcm.nf = n_unique_faces;
                rcm.comp_id = 1; rcm.self = &rcm;
                WeldRecoarsenParams rp;
                WeldCleanup_default_recoarsen_params(&rp);
                { const char *re;
                  if ((re = getenv("SEAM_RECOARSEN_BELOW"))) {
                      double v = atof(re); if (v > 0) rp.collapse_below = v; }
                  if ((re = getenv("SEAM_RECOARSEN_BAND"))) {
                      double v = atof(re); if (v > 0) rp.band = v; } }
                WeldRecoarsenStats rs;
                WeldCleanup_recoarsen_seam(arena, &rcm, rplanes, rnp, &rp, &rs);
                flat_faces = rcm.faces; n_unique_faces = rcm.nf;
                recoarsen_collapses = rs.n_collapses;
                recoarsen_faces_in = rs.faces_in;
                recoarsen_faces_out = rs.faces_out;
                fprintf(stderr,
                    "Recoarsen: %zu collapse(s), %zu flip(s), %zu -> %zu faces "
                    "(band %.1f, below %.1f vox, %zu plane(s))\n",
                    rs.n_collapses, rs.n_flips, rs.faces_in, rs.faces_out,
                    rp.band, rp.collapse_below, rnp);
            }
            gw_phase("recoarsen");
            dump_stage(arena, dump_stages_dir, stage_prefix, 9, "recoarsen",
                       out_verts, out_nv, flat_faces, n_unique_faces,
                       vert_cube_idx, cube_palette, cubes.n, color_nv);
        }

        /* Seam-band CVT beautification: re-mesh the recoarsened band with the
         * SAME CVT/RVD engine the per-cube interiors were built with, so the
         * weld is blue-noise CVT quality everywhere instead of collapse-
         * scarred. Pinned-boundary (junction rings + hole rims bit-exact ->
         * conforming stitch, open boundaries unchanged) and FAIL-CLOSED per
         * patch (conformity + Euler/manifold/connectivity gates; a rejected
         * patch -- e.g. tight core folds where RVD cells could jump -- keeps
         * its recoarsen geometry). Runs before the final winding repair so
         * dual-face winding is reconciled, and before the final audits which
         * re-verify everything. SEAM_NO_BANDCVT=1 skips; SEAM_BANDCVT_H
         * tunes the band target edge length. */
        if (!no_bridge && !no_cleanup && !getenv("SEAM_NO_BANDCVT")) {
            double bc_band = 6.0;
            { const char *bce = getenv("SEAM_BAND");
              if (bce) { double v = atof(bce); if (v > 0) bc_band = v; } }
            uint8_t *bc_used = (uint8_t *)ARENA_CALLOC(arena, (long)out_nv, 1L);
            for (size_t f = 0; f < n_unique_faces; f++) {
                bc_used[flat_faces[f*3+0]] = 1;
                bc_used[flat_faces[f*3+1]] = 1;
                bc_used[flat_faces[f*3+2]] = 1;
            }
            SeamPlane bplanes[64];
            size_t bnp = SeamPlanes_detect(out_verts, out_nv, bc_used,
                                           (double)seam_cube, bc_band,
                                           bplanes, 64);
            if (bnp > 0) {
                SeamBandCvtParams bp; SeamBandCvt_default_params(&bp);
                { const char *bce;
                  if ((bce = getenv("SEAM_BANDCVT_H"))) {
                      double v = atof(bce); if (v > 0) bp.target_h = v; }
                  if ((bce = getenv("SEAM_BANDCVT_BAND"))) {
                      double v = atof(bce); if (v > 0) bp.band = v; } }
                /* Two passes: A on the aligned tiling, B on a half-tile
                 * OFFSET tiling so A's pinned tile borders land tile-interior
                 * and are re-meshed away (the fine_frac filter keeps B off
                 * A's accepted interiors). */
                for (int bc_pass = 0; bc_pass < 2; bc_pass++) {
                    bp.offset_half = bc_pass;
                    SeamBandCvtStats bst;
                    float *bc_v = NULL; int32_t *bc_f = NULL; int32_t *bc_src = NULL;
                    size_t bc_nv = 0, bc_nf = 0, bc_nnew = 0;
                    if (SeamBandCvt_process(arena, out_verts, out_nv,
                                            flat_faces, n_unique_faces,
                                            bplanes, bnp, &bp,
                                            &bc_v, &bc_nv, &bc_f, &bc_nf,
                                            &bc_src, &bc_nnew, &bst) == 0
                        && bst.accepted > 0) {
                        if (bc_nv > vert_cap) {
                            fprintf(stderr, "grid_weld: vert capacity %zu "
                                    "exceeded by band CVT (%zu)\n",
                                    vert_cap, bc_nv);
                            RAISE(IO_Failed);
                        }
                        for (size_t v2 = 0; v2 < bc_nnew; v2++)
                            vert_cube_idx[out_nv + v2] = vert_cube_idx[bc_src[v2]];
                        out_verts = bc_v; out_nv = bc_nv;
                        flat_faces = bc_f; n_unique_faces = bc_nf;
                        color_nv = out_nv;
                        bandcvt_accepted += bst.accepted;
                        bandcvt_rejected += bst.rejected;
                        bandcvt_faces_in += bst.faces_band_in;
                        bandcvt_faces_out += bst.faces_band_out;
                    }
                    fprintf(stderr,
                        "Band CVT pass %c: %zu/%zu patch(es) accepted "
                        "(rej rc=%zu pin=%zu bnd=%zu chi=%zu nm=%zu vtx=%zu comp=%zu; "
                        "%zu small, %zu clean), band %zu -> %zu faces, "
                        "+%zu verts (h %.1f, tile %.0f)\n",
                        bc_pass ? 'B' : 'A', bst.accepted, bst.patches,
                        bst.rej_rc, bst.rej_pin, bst.rej_bnd, bst.rej_chi,
                        bst.rej_manifold, bst.rej_vtx, bst.rej_comp,
                        bst.skipped_small, bst.skipped_clean,
                        bst.faces_band_in, bst.faces_band_out,
                        bst.verts_added, bp.target_h, bp.tile);
                }
                /* pinch insurance: the per-patch gates prove edge-
                 * manifoldness, not vertex-links; let the guard mop any
                 * bowtie before the final audits (idle when clean). */
                if (bandcvt_accepted > 0) {
                    ComponentMesh gm; memset(&gm, 0, sizeof gm);
                    gm.verts = out_verts; gm.faces = flat_faces;
                    gm.nv = out_nv; gm.nf = n_unique_faces;
                    gm.comp_id = 1; gm.self = &gm;
                    ManifoldGuardStats mg2;
                    ManifoldGuard_process(arena, &gm, 1, 0, &mg2);
                    out_verts = gm.verts; flat_faces = gm.faces;
                    out_nv = gm.nv; n_unique_faces = gm.nf;
                    if (mg2.nm_edges_resolved || mg2.pinch_splits)
                        fprintf(stderr,
                            "  band-cvt guard: %zu NM edge(s), %zu pinch "
                            "split(s)\n", mg2.nm_edges_resolved,
                            mg2.pinch_splits);
                }
            }
            gw_phase("bandcvt");
            dump_stage(arena, dump_stages_dir, stage_prefix, 10, "bandcvt",
                       out_verts, out_nv, flat_faces, n_unique_faces,
                       vert_cube_idx, cube_palette, cubes.n, color_nv);
        }

        /* Re-run the intra-component winding BFS AFTER the late fills. Every
         * closer above (final pinhole, bowtie, seam-hole ear-clip, post-guard)
         * appends faces AFTER the main OrientMesh pass at line ~1426, and a fill
         * wound to the loop's own Newell normal can disagree with the sheet it
         * patches -> same_dir edges (0 in v2 while the holes stayed open, 468 once
         * they were filled). OrientMesh_consistent's live-winding BFS drives
         * same_dir back to 0 without changing topology -- it only flips face
         * winding, preserving component structure (so it cannot mask a merger:
         * the component COUNT is unchanged, verified by the post-repair audit). */
        if (!no_cleanup) {
            size_t ro_comp = 0, ro_flip = 0, ro_resid = 0;
            OrientMesh_consistent(arena, out_verts, out_nv, NULL,
                                  flat_faces, n_unique_faces,
                                  &ro_flip, &ro_comp, &ro_resid);
            fprintf(stderr, "Post-fill winding repair: %zu components, %zu "
                    "flipped, %zu residual same_dir\n", ro_comp, ro_flip, ro_resid);
        }
        gw_phase("orient2");

        /* Stage 11: final mesh after the last winding pass (== welded.obj). */
        dump_stage(arena, dump_stages_dir, stage_prefix, 11, "final",
                   out_verts, out_nv, flat_faces, n_unique_faces,
                   vert_cube_idx, cube_palette, cubes.n, color_nv);

        /* Post-repair audit + emit diagnostic OBJ at <out>.bad_edges.obj. */
        char diag_path[1024];
        snprintf(diag_path, sizeof(diag_path), "%s.bad_edges.obj", out_path);
        FILE *diag_fp = fopen(diag_path, "w");
        ms = manifold_audit(arena, flat_faces, n_unique_faces,
                            out_verts, diag_fp);
        if (diag_fp) {
            fclose(diag_fp);
            fprintf(stderr, "Wrote %s\n", diag_path);
        }

        /* Also emit a focused non-manifold-only OBJ with the actual
         * triangles touching each non-manifold edge plus their 1-ring
         * neighbours. This is the "let me see what is going on at this
         * specific bug" view that bad_edges.obj (line segments only)
         * cannot give. */
        if (ms.non_manifold > 0) {
            char nm_path[1024];
            snprintf(nm_path, sizeof(nm_path), "%s.nonmanifold.obj", out_path);
            emit_nonmanifold_neighborhood(arena, flat_faces, n_unique_faces,
                                          out_verts, out_nv, nm_path);
            fprintf(stderr, "Wrote %s (non-manifold neighbourhood)\n",
                    nm_path);
        }
        /* Vertex-manifold (pinch/bowtie) audit. The edge-multiplicity audit
         * above is BLIND to a vertex where two triangle fans meet at a single
         * point with no shared edge -- exactly what the seam-bridge BPA leaves.
         * Check it explicitly via the shared vertex-fan auditor. */
        MeshManifoldStats vms = MeshManifold_audit(arena, out_nv,
                                                   flat_faces, n_unique_faces);
        size_t pinch_verts = vms.nm_verts;
        fprintf(stderr,
            "Post-repair audit: unpaired=%zu non_manifold=%zu same_dir=%zu "
            "manifold=%zu pinch_verts=%zu\n",
            ms.unpaired, ms.non_manifold, ms.same_dir_pairs, ms.manifold_pairs,
            pinch_verts);
        gw_phase("audit_post");

        /* Count weld verts (verts touched by a BPA seam-bridge face).
         * Populated right after SeamWeld_bridge, so only meaningful for the
         * first color_nv verts -- any verts pinhole split in beyond that are
         * never bridge verts. */
        size_t n_weld_verts = 0;
        for (size_t v = 0; v < color_nv; v++) {
            if (vert_is_weld[v]) n_weld_verts++;
        }
        fprintf(stderr, "Welded verts (shared across >=2 cubes): %zu\n",
                n_weld_verts);

        /* Write OBJ: per-vert color (cube palette for normal verts,
         * bright green for weld verts). Faces interpolate vert colors,
         * so welded edges between cubes will appear as green seams. */
        FILE *ofp = fopen(out_path, "w");
        if (!ofp) {
            fprintf(stderr, "grid_weld: cannot open %s\n", out_path);
            RAISE(IO_Failed);
        }
        /* Default: colour by CONNECTED COMPONENT so wrap shatter / cross-wrap
         * fusion is visible on load (mesh_render --vcolor, no obj_cc_color step).
         * GW_COLOR_BY=cube keeps the old per-cube / green-weld / blue-pinhole
         * provenance view for seam-weld debugging. */
        const char *gw_color_by = getenv("GW_COLOR_BY");
        int by_cube = (gw_color_by && !strcmp(gw_color_by, "cube"));
        float *cc = NULL;
        if (!by_cube) {
            cc = (float *)malloc(out_nv * 3 * sizeof(float));
            if (cc) { CCColorOpts opts; CCColor_default_opts(&opts); opts.sat = 0.75;
                      CCColor_compute(out_nv, flat_faces, n_unique_faces, &opts, cc, NULL); }
        }
        fprintf(ofp, "# grid_weld output: %s vertex colors\n",
                cc ? "connected-component" : "per-cube provenance");
        fprintf(ofp, "# %zu verts (%zu weld), %zu faces\n",
                out_nv, n_weld_verts, n_unique_faces);
        for (size_t v = 0; v < out_nv; v++) {
            float cr, cg, cb;
            if (cc) {
                cr = cc[v*3+0]; cg = cc[v*3+1]; cb = cc[v*3+2];
            } else {
                const float *c;
                int is_weld = (v < color_nv) ? vert_is_weld[v] : 0;
                if (is_weld)            c = SEAM_COLOR;      /* green: bridge weld vert */
                else if (v >= color_nv) c = PINHOLE_COLOR;   /* blue: pinhole/holefill vert */
                else { int cube = vert_cube_idx[v];
                       int pal = (cube >= 0 && cube < (int)cubes.n) ? cube_palette[cube] : 0;
                       c = CUBE_PALETTE[pal]; }
                cr = c[0]; cg = c[1]; cb = c[2];
            }
            fprintf(ofp, "v %.6f %.6f %.6f %.4f %.4f %.4f\n",
                (double)out_verts[v * 3 + 0],
                (double)out_verts[v * 3 + 1],
                (double)out_verts[v * 3 + 2],
                (double)cr, (double)cg, (double)cb);
        }
        free(cc);
        for (size_t f = 0; f < n_unique_faces; f++) {
            fprintf(ofp, "f %d %d %d\n",
                flat_faces[f * 3 + 0] + 1,
                flat_faces[f * 3 + 1] + 1,
                flat_faces[f * 3 + 2] + 1);
        }
        fclose(ofp);
        fprintf(stderr, "Wrote %s (%zu verts, %zu green weld verts, %zu faces)\n",
                out_path, out_nv, n_weld_verts, n_unique_faces);
        gw_phase("write_obj");

        /* Optional weld-vert sidecar (consumed by the winding
         * diagnostic to compute its headline weld-vert disagreement
         * rate). Format: uint32 nv; uint8 is_weld[nv]. */
        if (weld_verts_out) {
            FILE *wfp = fopen(weld_verts_out, "wb");
            if (!wfp) {
                fprintf(stderr,
                    "grid_weld: cannot open weld-verts file %s\n",
                    weld_verts_out);
            } else {
                /* Only the first color_nv verts can be weld verts (bridge uses
                 * pre-existing verts; any pinhole-split verts beyond color_nv
                 * are never marked). Clamp so we never read past the
                 * color_nv-sized vert_is_weld bitmap. */
                uint32_t nv32 = (uint32_t)color_nv;
                fwrite(&nv32, sizeof(uint32_t), 1, wfp);
                fwrite(vert_is_weld, sizeof(uint8_t), color_nv, wfp);
                fclose(wfp);
                fprintf(stderr,
                    "Wrote %s (%zu weld-vert flags)\n",
                    weld_verts_out, color_nv);
            }
        }

        /* Write weld report. */
        char report_path[1024];
        snprintf(report_path, sizeof(report_path), "%s.weld_report.json",
                 out_path);
        FILE *rp = fopen(report_path, "w");
        if (rp) {
            fprintf(rp,
                "{\n"
                "  \"cubes_processed\": %zu,\n"
                "  \"total_input_verts\": %zu,\n"
                "  \"total_unique_verts\": %zu,\n"
                "  \"total_input_faces\": %zu,\n"
                "  \"total_unique_faces\": %zu,\n"
                "  \"recoarsen\": {\n"
                "    \"collapses\": %zu,\n"
                "    \"faces_in\":  %zu,\n"
                "    \"faces_out\": %zu\n"
                "  },\n"
                "  \"band_cvt\": {\n"
                "    \"patches_accepted\": %zu,\n"
                "    \"patches_rejected\": %zu,\n"
                "    \"band_faces_in\":  %zu,\n"
                "    \"band_faces_out\": %zu\n"
                "  },\n"
                "  \"manifold_audit\": {\n"
                "    \"unpaired\":     %zu,\n"
                "    \"non_manifold\": %zu,\n"
                "    \"same_dir_pairs\": %zu,\n"
                "    \"manifold_pairs\": %zu,\n"
                "    \"pinch_verts\":   %zu\n"
                "  }\n"
                "}\n",
                cubes.n, total_in_verts, out_nv, total_in_faces,
                n_unique_faces,
                recoarsen_collapses, recoarsen_faces_in, recoarsen_faces_out,
                bandcvt_accepted, bandcvt_rejected,
                bandcvt_faces_in, bandcvt_faces_out,
                ms.unpaired, ms.non_manifold,
                ms.same_dir_pairs, ms.manifold_pairs, pinch_verts);
            fclose(rp);
            fprintf(stderr, "Wrote %s\n", report_path);
        }

        /* Manifold-by-construction guarantee: every edge must have exactly
         * 2 incident faces, except OUTER-boundary edges (the "unpaired" count).
         * NON-MANIFOLD edges (>2 faces) or pinch/bowtie vertices are genuine
         * topology breakage -> HARD failure (exit 1). A residual same_dir edge
         * (a doubled or oppositely-wound face that survives the winding repair)
         * is still edge+vertex manifold -- a cosmetic blemish, so it WARNS but
         * does NOT fail the exit code. This lets a batch of welds distinguish
         * real failures (crash/OOM/non-manifold) from these minor residuals.
         * bad_edges.obj localizes whatever the audit flagged. */
        if (ms.non_manifold > 0 || pinch_verts > 0) {
            fprintf(stderr,
                "ERROR: manifold audit failed -- NON-MANIFOLD (non_manifold=%zu "
                "pinch_verts=%zu same_dir=%zu)\n",
                ms.non_manifold, pinch_verts, ms.same_dir_pairs);
            ok = 0;
        } else if (ms.same_dir_pairs > 0) {
            fprintf(stderr,
                "WARNING: %zu same_dir edge(s) (doubled/flipped face; edge+vertex "
                "manifold) -- cosmetic, NOT failing exit. See bad_edges.obj\n",
                ms.same_dir_pairs);
            ok = 1;
        } else {
            ok = 1;
        }

    EXCEPT(Arena_Failed)
        fprintf(stderr, "grid_weld: arena OOM\n");
    EXCEPT(IO_Failed)
        fprintf(stderr, "grid_weld: I/O failure\n");
    END_TRY;

    Arena_dispose(&arena);
    return ok ? 0 : 1;
}
