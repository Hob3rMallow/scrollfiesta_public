/* reeb_probe.c -- per-component bridge diagnosis for welded scroll meshes.
 *
 * Confirms (or refutes) that a welded component still contains inter-wrap
 * BRIDGES that the UV-difference oracle is blind to, and LOCALIZES them, by
 * cross-checking three independent signals per connected component:
 *
 *   genus       (mesh_topo.h)  -- a bridge that forms a handle gives genus >= 1
 *                                 (handle-only alarm; a clean wrap annulus is g0).
 *   reeb necks  (split/reeb.h) -- a thin topological neck = interior minimum of
 *                                 the geodesic cut profile (catches the thin-neck
 *                                 MERGE that does NOT raise genus, too).
 *   side-facing (split/neck_probe.h) -- a triangle whose normal faces side-to-
 *                                 side relative to its neighbourhood ("one thin
 *                                 triangle facing side-to-side").
 *
 * A high-confidence bridge face is one that is BOTH on a reeb neck band AND
 * side-facing. The colored dump paints: gray = body, orange = reeb neck,
 * blue = side-facing, RED = both (the bridge).
 *
 * Usage:
 *   reeb_probe <mesh.obj> [--largest | --comp N] [--min-faces M]
 *                         [--dump-bridges out.obj] [--csv stats.csv]
 *   reeb_probe --selftest
 *
 *   --largest      analyze (and dump) only the largest component.
 *   --comp N       analyze only component N (1-based, descending face count).
 *   --min-faces M  ignore components with < M faces (default 64).
 *   --dump-bridges write the analyzed components recolored by signal.
 *   --csv          append a one-row-per-component CSV.
 *
 * Exit: 0 = analysis ran (a detected bridge is a FINDING, not a tool error),
 *       1 = IO error, 2 = usage error, 3 = selftest fail.
 *
 * Project-linked (calls the real Reeb_analyze / NeckProbe_scan / MeshTopo_analyze
 * the pipeline will use), modeled on obj_components.c / manifold_check.c.
 */
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#include "../common/arena.h"
#include "../common/obj_io.h"
#include "../common/union_find.h"
#include "../flatten/mesh_topo.h"
#include "../split/reeb.h"
#include "../split/neck_probe.h"
#include "../flatten/seam_cut.h"

/* category priorities for the colored dump */
enum { CAT_NONE = 0, CAT_REEB = 1, CAT_SIDE = 2, CAT_BOTH = 3 };
static const float CAT_COL[4][3] = {
    { 0.45f, 0.45f, 0.45f },   /* body  */
    { 0.95f, 0.55f, 0.10f },   /* reeb neck   */
    { 0.20f, 0.45f, 0.95f },   /* side-facing */
    { 0.95f, 0.12f, 0.12f }    /* both = bridge */
};

typedef struct {
    int     id;
    size_t  nv, nf;
    double  genus;
    long    boundary_loops, nm_edges;
    size_t  reeb_necks;
    int32_t bottleneck;
    size_t  sidefacing;
    size_t  n_sf_clusters;  /* connected groups of side-facing triangles */
    size_t  n_bridges;      /* clusters whose removal cuts a handle (genus drop) */
    double  genus_cut;      /* component genus after removing all bridge clusters */
    size_t  bridge_faces;   /* faces in confirmed bridge clusters */
    size_t  n_handle_loops; /* topological handle loops (one per genus) */
} CompReport;

/* distinct bright colors for the handle loops (cycled) */
static const float LOOP_PAL[6][3] = {
    { 0.95f, 0.10f, 0.10f },   /* red    */
    { 0.10f, 0.90f, 0.30f },   /* green  */
    { 0.95f, 0.85f, 0.10f },   /* yellow */
    { 0.20f, 0.55f, 1.00f },   /* blue   */
    { 0.85f, 0.30f, 0.95f },   /* magenta*/
    { 0.10f, 0.90f, 0.90f }    /* cyan   */
};

/* Build a compacted sub-mesh from a list of global face indices, recording the
 * map back to global face / vertex indices for the dump. All arena-allocated. */
typedef struct {
    size_t   nv, nf;
    float   *verts;     /* [nv*3] */
    int32_t *faces;     /* [nf*3] local */
    int32_t *gface;     /* [nf] local face -> global face */
    int32_t *gvert;     /* [nv] local vert -> global vert */
} SubMesh;

static void build_submesh(Arena_T arena, const float *V, size_t gnv,
                          const int32_t *F, const int32_t *gfaces, size_t nfc,
                          SubMesh *out)
{
    int32_t *remap = (int32_t *)ARENA_ALLOC(arena, (long)(gnv * sizeof(int32_t)));
    for (size_t i = 0; i < gnv; i++) remap[i] = -1;

    size_t lnv = 0;
    for (size_t i = 0; i < nfc; i++) {
        int32_t gf = gfaces[i];
        for (int k = 0; k < 3; k++) {
            int32_t gv = F[(size_t)gf * 3 + (size_t)k];
            if (gv >= 0 && (size_t)gv < gnv && remap[gv] < 0)
                remap[gv] = (int32_t)lnv++;
        }
    }
    float   *verts = (float *)  ARENA_ALLOC(arena, (long)(lnv * 3 * sizeof(float)));
    int32_t *gvert = (int32_t *)ARENA_ALLOC(arena, (long)(lnv * sizeof(int32_t)));
    for (size_t gv = 0; gv < gnv; gv++)
        if (remap[gv] >= 0) {
            size_t lv = (size_t)remap[gv];
            verts[lv * 3 + 0] = V[gv * 3 + 0];
            verts[lv * 3 + 1] = V[gv * 3 + 1];
            verts[lv * 3 + 2] = V[gv * 3 + 2];
            gvert[lv] = (int32_t)gv;
        }
    int32_t *faces = (int32_t *)ARENA_ALLOC(arena, (long)(nfc * 3 * sizeof(int32_t)));
    int32_t *gface = (int32_t *)ARENA_ALLOC(arena, (long)(nfc * sizeof(int32_t)));
    for (size_t i = 0; i < nfc; i++) {
        int32_t gf = gfaces[i];
        for (int k = 0; k < 3; k++)
            faces[i * 3 + (size_t)k] = remap[F[(size_t)gf * 3 + (size_t)k]];
        gface[i] = gf;
    }
    out->nv = lnv; out->nf = nfc;
    out->verts = verts; out->faces = faces; out->gface = gface; out->gvert = gvert;
}

/* True total genus of the mesh with the faces in `drop` removed, AND its
 * connected-component count (via *out_C). The component count matters: the
 * (2 - loops - euler)/2 genus formula assumes ONE component, so removing a
 * cluster that DISCONNECTS a piece would spuriously look like a genus drop.
 * For C components the correct total genus is (2C - loops - euler_ref)/2, where
 * euler_ref discounts vertices left unreferenced by the removal. A handle cut
 * keeps C unchanged and drops genus; a separating cut raises C and preserves it. */
static double genus_without(Arena_T arena, const float *V, size_t nv,
                            const int32_t *F, size_t nf, const uint8_t *drop,
                            int *out_C)
{
    Arena_Mark m = Arena_save(arena);
    int32_t *FF = (int32_t *)ARENA_ALLOC(arena, (long)(nf * 3 * sizeof(int32_t)));
    size_t k = 0;
    for (size_t fi = 0; fi < nf; fi++) {
        if (drop[fi]) continue;
        FF[k * 3 + 0] = F[fi * 3 + 0];
        FF[k * 3 + 1] = F[fi * 3 + 1];
        FF[k * 3 + 2] = F[fi * 3 + 2];
        k++;
    }
    double g = 0.0;
    int C = 0;
    if (k > 0) {
        MeshTopoInfo t;
        if (MeshTopo_analyze(arena, V, nv, FF, k, &t) == 0) {
            /* connected components over the remaining faces' vertices */
            UnionFind uf = UF_new(arena, (int32_t)nv);
            uint8_t *ref = (uint8_t *)ARENA_CALLOC(arena, (long)nv, 1L);
            for (size_t i = 0; i < k; i++) {
                int32_t a = FF[i * 3 + 0], b = FF[i * 3 + 1], c = FF[i * 3 + 2];
                uf_union(&uf, a, b); uf_union(&uf, a, c);
                ref[a] = 1; ref[b] = 1; ref[c] = 1;
            }
            uint8_t *seen = (uint8_t *)ARENA_CALLOC(arena, (long)nv, 1L);
            for (size_t v = 0; v < nv; v++)
                if (ref[v]) { int32_t r = uf_find(&uf, (int32_t)v); if (!seen[r]) { seen[r] = 1; C++; } }
            double euler_ref = (double)t.euler - (double)t.n_unref_verts;
            g = ((double)(2 * C) - (double)t.n_boundary_loops - euler_ref) / 2.0;
        }
    }
    if (out_C) *out_C = C;
    Arena_restore(arena, m);
    return g;
}

/* Analyze one sub-mesh: genus + reeb necks + side-facing, then DEFINITIVELY
 * localize handle bridges -- cluster the side-facing triangles and remove each
 * cluster; a cluster whose removal drops the component genus is cutting a handle,
 * i.e. an inter-wrap bridge. If `prio`/`col` non-NULL, paint global vertices:
 * RED = confirmed bridge, blue = side-facing, orange = reeb neck. */
static void analyze_sub(Arena_T arena, const SubMesh *s, const int32_t *Fglobal,
                        CompReport *rep, uint8_t *prio, float *col, int verbose,
                        int genus_only)
{
    MeshTopoInfo topo;
    int trc = MeshTopo_analyze(arena, s->verts, s->nv, s->faces, s->nf, &topo);

    /* Fast census: genus is the bridge COUNT; skip the (slow) neck/cluster/side-
     * facing work. When dumping, still localize + color the handle LOOPS (cheap:
     * the tree-cotree absorbs the boundary, so only 2*genus generators are
     * tested) -- a clean highlight with no geometric noise. */
    if (genus_only) {
        rep->nv = s->nv; rep->nf = s->nf;
        rep->genus          = (trc == 0) ? topo.genus : -1.0;
        rep->boundary_loops = (trc == 0) ? topo.n_boundary_loops : -1;
        rep->nm_edges       = (trc == 0) ? topo.n_nonmanifold_edges : -1;
        size_t nloops = (trc == 0 && topo.genus >= 1.0) ? (size_t)(topo.genus + 0.5) : 0;

        if (prio && col && trc == 0 && topo.genus >= 1.0) {
            HandleLoop *loops = NULL; size_t nl = 0;
            if (SeamCut_handle_loops(arena, s->verts, s->nv, s->faces, s->nf,
                                     &loops, &nl) == 0) {
                nloops = nl;
                /* Mark each loop's verts with its palette index (1-based), then
                 * dilate a few rings over the faces so a thin loop becomes a
                 * visible patch on the big mesh. */
                uint8_t *lcidx = (uint8_t *)ARENA_CALLOC(arena, (long)s->nv, 1L);
                for (size_t li = 0; li < nl; li++) {
                    if (verbose)
                        printf("           handle loop %zu: length=%.1f vox, %zu verts\n",
                               li, loops[li].length, loops[li].n);
                    uint8_t tag = (uint8_t)((li % 6) + 1);
                    for (size_t t = 0; t < loops[li].n; t++) {
                        int32_t lv = loops[li].verts[t];
                        if (lv >= 0 && (size_t)lv < s->nv) lcidx[lv] = tag;
                    }
                }
                for (int ring = 0; ring < 3; ring++)        /* ~3-ring dilation */
                    for (size_t fi = 0; fi < s->nf; fi++) {
                        int32_t a = s->faces[fi*3], b = s->faces[fi*3+1], c = s->faces[fi*3+2];
                        uint8_t t = 0;
                        if (lcidx[a] > t) t = lcidx[a];
                        if (lcidx[b] > t) t = lcidx[b];
                        if (lcidx[c] > t) t = lcidx[c];
                        if (!t) continue;
                        if (!lcidx[a]) lcidx[a] = t;
                        if (!lcidx[b]) lcidx[b] = t;
                        if (!lcidx[c]) lcidx[c] = t;
                    }
                for (size_t lv = 0; lv < s->nv; lv++) {
                    if (!lcidx[lv]) continue;
                    int32_t gv = s->gvert[lv];
                    if (gv < 0) continue;
                    const float *lc = LOOP_PAL[(lcidx[lv] - 1) % 6];
                    prio[gv] = 255;
                    col[(size_t)gv * 3 + 0] = lc[0];
                    col[(size_t)gv * 3 + 1] = lc[1];
                    col[(size_t)gv * 3 + 2] = lc[2];
                }
            }
        }
        rep->n_handle_loops = nloops;
        if (verbose)
            printf("  comp %2d: V=%-7zu F=%-7zu genus=%5.1f loops=%-4ld nmE=%-4ld | "
                   "bridges(genus)=%zu\n",
                   rep->id, rep->nv, rep->nf, rep->genus, rep->boundary_loops,
                   rep->nm_edges, rep->n_handle_loops);
        return;
    }

    ReebResult rr;
    int rrc = Reeb_analyze(arena, s->verts, s->nv, s->faces, s->nf, NULL, &rr);

    NeckProbeResult np;
    int nprc = NeckProbe_scan(arena, s->verts, s->nv, s->faces, s->nf, NULL, &np);

    rep->nv = s->nv; rep->nf = s->nf;
    rep->genus          = (trc == 0) ? topo.genus : -1.0;
    rep->boundary_loops = (trc == 0) ? topo.n_boundary_loops : -1;
    rep->nm_edges       = (trc == 0) ? topo.n_nonmanifold_edges : -1;
    rep->reeb_necks     = (rrc == 0) ? rr.n_necks : 0;
    rep->bottleneck     = (rrc == 0) ? rr.bottleneck : -1;
    rep->sidefacing     = (nprc == 0) ? np.n_sidefacing : 0;

    const uint8_t *sd = (nprc == 0) ? np.face_flag : NULL;
    const uint8_t *rk = (rrc == 0) ? rr.face_flag : NULL;

    /* Cluster the side-facing triangles (union faces sharing a vertex). */
    uint8_t *bridgeface = (uint8_t *)ARENA_CALLOC(arena, (long)s->nf, 1L);
    size_t n_clusters = 0, n_bridges = 0, n_bridge_faces = 0;
    double genus_cut = rep->genus;

    if (sd && rep->genus >= 1.0) {
        UnionFind cuf = UF_new(arena, (int32_t)s->nf);
        int32_t *firstf = (int32_t *)ARENA_ALLOC(arena, (long)(s->nv * sizeof(int32_t)));
        for (size_t v = 0; v < s->nv; v++) firstf[v] = -1;
        for (size_t lf = 0; lf < s->nf; lf++) {
            if (!sd[lf]) continue;
            for (int k = 0; k < 3; k++) {
                int32_t v = s->faces[lf * 3 + (size_t)k];
                if (v < 0) continue;
                if (firstf[v] < 0) firstf[v] = (int32_t)lf;
                else uf_union(&cuf, (int32_t)lf, firstf[v]);
            }
        }
        /* Index clusters by root; bucket their faces (CSR). */
        int32_t *clidx = (int32_t *)ARENA_ALLOC(arena, (long)(s->nf * sizeof(int32_t)));
        for (size_t lf = 0; lf < s->nf; lf++) clidx[lf] = -1;
        for (size_t lf = 0; lf < s->nf; lf++) {
            if (!sd[lf]) continue;
            int32_t r = uf_find(&cuf, (int32_t)lf);
            if (clidx[r] < 0) clidx[r] = (int32_t)n_clusters++;
        }
        int32_t *cloff = (int32_t *)ARENA_CALLOC(arena, (long)(n_clusters + 1), (long)sizeof(int32_t));
        for (size_t lf = 0; lf < s->nf; lf++)
            if (sd[lf]) cloff[clidx[uf_find(&cuf, (int32_t)lf)] + 1]++;
        for (size_t c = 1; c <= n_clusters; c++) cloff[c] += cloff[c - 1];
        int32_t *clf = (int32_t *)ARENA_ALLOC(arena, (long)(((size_t)np.n_sidefacing ? np.n_sidefacing : 1) * sizeof(int32_t)));
        int32_t *clcur = (int32_t *)ARENA_ALLOC(arena, (long)(n_clusters * sizeof(int32_t)));
        for (size_t c = 0; c < n_clusters; c++) clcur[c] = cloff[c];
        for (size_t lf = 0; lf < s->nf; lf++)
            if (sd[lf]) { int32_t c = clidx[uf_find(&cuf, (int32_t)lf)]; clf[clcur[c]++] = (int32_t)lf; }

        /* Test each cluster: removing it must DROP the genus while keeping the
         * component connected (C==1) -- that is a handle (inter-wrap bridge) cut,
         * not a separating cut that merely lops off a fin. */
        uint8_t *drop = (uint8_t *)ARENA_CALLOC(arena, (long)s->nf, 1L);
        uint8_t *alldrop = (uint8_t *)ARENA_CALLOC(arena, (long)s->nf, 1L);
        for (size_t c = 0; c < n_clusters; c++) {
            size_t a = (size_t)cloff[c], b = (size_t)cloff[c + 1];
            for (size_t t = a; t < b; t++) drop[clf[t]] = 1;
            int cc = 0;
            double gc = genus_without(arena, s->verts, s->nv, s->faces, s->nf, drop, &cc);
            int is_bridge = (cc == 1 && rep->genus - gc >= 0.5);
            if (is_bridge) {
                n_bridges++;
                for (size_t t = a; t < b; t++) {
                    bridgeface[clf[t]] = 1;
                    alldrop[clf[t]] = 1;
                    n_bridge_faces++;
                }
                if (verbose)
                    printf("           bridge cluster %zu: %zu side-facing tris, "
                           "genus %.1f -> %.1f (components=%d)\n",
                           c, b - a, rep->genus, gc, cc);
            }
            for (size_t t = a; t < b; t++) drop[clf[t]] = 0;
        }
        if (n_bridges > 0) {
            int cc = 0;
            genus_cut = genus_without(arena, s->verts, s->nv, s->faces, s->nf, alldrop, &cc);
        }
    }

    rep->n_sf_clusters = n_clusters;
    rep->n_bridges     = n_bridges;
    rep->genus_cut     = genus_cut;
    rep->bridge_faces  = n_bridge_faces;

    /* Color global vertices by category (RED bridge > blue side > orange reeb). */
    if (prio && col) {
        for (size_t lf = 0; lf < s->nf; lf++) {
            int cat = CAT_NONE;
            if (bridgeface[lf])             cat = CAT_BOTH;
            else if (sd && sd[lf])          cat = CAT_SIDE;
            else if (rk && rk[lf])          cat = CAT_REEB;
            if (cat == CAT_NONE) continue;
            int32_t gf = s->gface[lf];
            for (int k = 0; k < 3; k++) {
                int32_t gv = Fglobal[(size_t)gf * 3 + (size_t)k];
                if (gv < 0) continue;
                if ((uint8_t)cat > prio[gv]) {
                    prio[gv] = (uint8_t)cat;
                    col[(size_t)gv * 3 + 0] = CAT_COL[cat][0];
                    col[(size_t)gv * 3 + 1] = CAT_COL[cat][1];
                    col[(size_t)gv * 3 + 2] = CAT_COL[cat][2];
                }
            }
        }
    }

    /* Topological handle loops -- ALL of them, one per genus (this is the
     * localizer that finds the bridges the side-facing heuristic misses). Painted
     * LAST with a per-loop palette so each handle is visually distinct (overrides). */
    HandleLoop *loops = NULL; size_t nloops = 0;
    if (rep->genus >= 1.0 &&
        SeamCut_handle_loops(arena, s->verts, s->nv, s->faces, s->nf, &loops, &nloops) == 0) {
        for (size_t li = 0; li < nloops; li++) {
            const float *lc = LOOP_PAL[li % 6];
            if (verbose)
                printf("           handle loop %zu: length=%.1f vox, %zu loop-verts\n",
                       li, loops[li].length, loops[li].n);
            if (prio && col)
                for (size_t t = 0; t < loops[li].n; t++) {
                    int32_t lv = loops[li].verts[t];
                    if (lv < 0 || (size_t)lv >= s->nv) continue;
                    int32_t gv = s->gvert[lv];
                    if (gv < 0) continue;
                    prio[gv] = 255;
                    col[(size_t)gv * 3 + 0] = lc[0];
                    col[(size_t)gv * 3 + 1] = lc[1];
                    col[(size_t)gv * 3 + 2] = lc[2];
                }
        }
    }
    rep->n_handle_loops = nloops;

    if (verbose) {
        printf("  comp %2d: V=%-7zu F=%-7zu genus=%5.1f loops=%-3ld nmE=%-3ld | "
               "HANDLE LOOPS=%zu (the bridges) | side-facing=%zu in %zu clusters "
               "(%zu cut a handle) | reeb necks=%zu\n",
               rep->id, rep->nv, rep->nf, rep->genus, rep->boundary_loops,
               rep->nm_edges, rep->n_handle_loops, rep->sidefacing,
               rep->n_sf_clusters, rep->n_bridges, rep->reeb_necks);
    }
}

/* ---- selftest ---- */
static int selftest(void)
{
    int f = 0;
    printf("=== reeb_probe selftest ===\n");
    f += MeshTopo_selftest();
    f += NeckProbe_selftest();
    f += Reeb_selftest();
    f += SeamCut_selftest();
    printf("=== reeb_probe selftest %s (%d failure%s) ===\n",
           f ? "FAILED" : "PASSED", f, f == 1 ? "" : "s");
    return f ? 3 : 0;
}

int main(int argc, char **argv)
{
    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0) return selftest();
    if (argc < 2) {
        fprintf(stderr,
            "Usage: %s <mesh.obj> [--largest|--comp N] [--min-faces M]\n"
            "                 [--dump-bridges out.obj] [--csv stats.csv]\n"
            "                 [--sever out.obj] [--max-len L=60]\n"
            "       %s --selftest\n", argv[0], argv[0]);
        return 2;
    }
    const char *path = argv[1];
    const char *dump = NULL, *csv = NULL, *sever_out = NULL;
    int only_largest = 0, only_comp = -1, genus_only = 0;
    size_t min_faces = 64;
    double sever_max_len = 60.0;
    for (int i = 2; i < argc; i++) {
        if      (!strcmp(argv[i], "--largest")) only_largest = 1;
        else if (!strcmp(argv[i], "--comp") && i + 1 < argc) only_comp = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--min-faces") && i + 1 < argc) min_faces = (size_t)strtoul(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--dump-bridges") && i + 1 < argc) dump = argv[++i];
        else if (!strcmp(argv[i], "--csv") && i + 1 < argc) csv = argv[++i];
        else if (!strcmp(argv[i], "--sever") && i + 1 < argc) sever_out = argv[++i];
        else if (!strcmp(argv[i], "--max-len") && i + 1 < argc) sever_max_len = atof(argv[++i]);
        else if (!strcmp(argv[i], "--genus-only")) genus_only = 1;
        else { fprintf(stderr, "unknown arg: %s\n", argv[i]); return 2; }
    }

    Arena_T arena = Arena_new();
    float *V = NULL; int32_t *F = NULL; size_t nv = 0, nf = 0;
    if (ObjIO_read(arena, path, &V, &nv, &F, &nf) != 0) {
        fprintf(stderr, "ERROR: cannot read %s\n", path);
        Arena_dispose(&arena);
        return 1;
    }
    printf("=== reeb_probe: %s ===\n", path);
    printf("  V=%zu  F=%zu\n", nv, nf);
    if (nv == 0 || nf == 0) { Arena_dispose(&arena); return 0; }

    clock_t t0 = clock();

    /* connected components via union-find over face edges */
    UnionFind uf = UF_new(arena, (int32_t)nv);
    for (size_t fi = 0; fi < nf; fi++) {
        int32_t a = F[fi * 3 + 0], b = F[fi * 3 + 1], c = F[fi * 3 + 2];
        if (a >= 0 && b >= 0 && (size_t)a < nv && (size_t)b < nv) uf_union(&uf, a, b);
        if (a >= 0 && c >= 0 && (size_t)a < nv && (size_t)c < nv) uf_union(&uf, a, c);
    }
    /* faces per root */
    int32_t *facecnt = (int32_t *)ARENA_CALLOC(arena, (long)nv, (long)sizeof(int32_t));
    for (size_t fi = 0; fi < nf; fi++) {
        int32_t a = F[fi * 3 + 0];
        if (a >= 0 && (size_t)a < nv) facecnt[uf_find(&uf, a)]++;
    }
    /* collect roots with faces, sort by descending count */
    size_t ncomp = 0;
    for (size_t v = 0; v < nv; v++) if (facecnt[v] > 0) ncomp++;
    int32_t *roots = (int32_t *)ARENA_ALLOC(arena, (long)((ncomp ? ncomp : 1) * sizeof(int32_t)));
    size_t ci = 0;
    for (size_t v = 0; v < nv; v++) if (facecnt[v] > 0) roots[ci++] = (int32_t)v;
    /* simple insertion sort by facecnt desc (ncomp is small in practice) */
    for (size_t i = 1; i < ncomp; i++) {
        int32_t key = roots[i]; size_t j = i;
        while (j > 0 && facecnt[roots[j - 1]] < facecnt[key]) { roots[j] = roots[j - 1]; j--; }
        roots[j] = key;
    }
    /* root -> component index */
    int32_t *root2idx = (int32_t *)ARENA_ALLOC(arena, (long)(nv * sizeof(int32_t)));
    for (size_t v = 0; v < nv; v++) root2idx[v] = -1;
    for (size_t k = 0; k < ncomp; k++) root2idx[roots[k]] = (int32_t)k;
    /* bucket global faces by component (CSR) */
    int32_t *coff = (int32_t *)ARENA_CALLOC(arena, (long)(ncomp + 1), (long)sizeof(int32_t));
    for (size_t fi = 0; fi < nf; fi++) {
        int32_t a = F[fi * 3 + 0];
        if (a < 0 || (size_t)a >= nv) continue;
        coff[root2idx[uf_find(&uf, a)] + 1]++;
    }
    for (size_t k = 1; k <= ncomp; k++) coff[k] += coff[k - 1];
    int32_t *cfaces = (int32_t *)ARENA_ALLOC(arena, (long)(((size_t)coff[ncomp] ? (size_t)coff[ncomp] : 1) * sizeof(int32_t)));
    int32_t *cur = (int32_t *)ARENA_ALLOC(arena, (long)(ncomp * sizeof(int32_t)));
    for (size_t k = 0; k < ncomp; k++) cur[k] = coff[k];
    for (size_t fi = 0; fi < nf; fi++) {
        int32_t a = F[fi * 3 + 0];
        if (a < 0 || (size_t)a >= nv) continue;
        int32_t cidx = root2idx[uf_find(&uf, a)];
        cfaces[cur[cidx]++] = (int32_t)fi;
    }

    printf("  components: %zu (>= %zu faces analyzed)\n", ncomp, min_faces);

    /* dump buffers (whole mesh, gray base) */
    uint8_t *prio = NULL; float *col = NULL;
    if (dump) {
        prio = (uint8_t *)ARENA_CALLOC(arena, (long)nv, 1L);
        col  = (float *)  ARENA_ALLOC(arena, (long)(nv * 3 * sizeof(float)));
        for (size_t i = 0; i < nv; i++) {
            col[i * 3 + 0] = CAT_COL[CAT_NONE][0];
            col[i * 3 + 1] = CAT_COL[CAT_NONE][1];
            col[i * 3 + 2] = CAT_COL[CAT_NONE][2];
        }
    }

    FILE *cf = NULL;
    if (csv) {
        cf = fopen(csv, "w");
        if (cf) fprintf(cf, "comp,nv,nf,genus,boundary_loops,nm_edges,handle_loops,sidefacing,sf_clusters,sf_bridges,reeb_necks\n");
        else fprintf(stderr, "  WARN: cannot write csv %s\n", csv);
    }

    size_t analyzed = 0, total_bridge_faces = 0, comps_with_bridge = 0;
    for (size_t k = 0; k < ncomp; k++) {
        size_t cnf = (size_t)(coff[k + 1] - coff[k]);
        int is_largest = (k == 0);
        if (only_largest && !is_largest) continue;
        if (only_comp >= 0 && (int)k != only_comp - 1) continue;
        if (cnf < min_faces) continue;

        SubMesh s;
        Arena_Mark m = Arena_save(arena);
        build_submesh(arena, V, nv, F, &cfaces[coff[k]], cnf, &s);

        CompReport rep; memset(&rep, 0, sizeof(rep));
        rep.id = (int)k + 1;
        analyze_sub(arena, &s, F, &rep, prio, col, 1, genus_only);
        analyzed++;
        total_bridge_faces += rep.n_handle_loops;
        if (rep.n_handle_loops > 0) comps_with_bridge++;

        if (cf)
            fprintf(cf, "%d,%zu,%zu,%.1f,%ld,%ld,%zu,%zu,%zu,%zu,%zu\n",
                    rep.id, rep.nv, rep.nf, rep.genus, rep.boundary_loops,
                    rep.nm_edges, rep.n_handle_loops, rep.sidefacing,
                    rep.n_sf_clusters, rep.n_bridges, rep.reeb_necks);

        Arena_restore(arena, m);
    }
    if (cf) fclose(cf);

    double secs = (double)(clock() - t0) / (double)CLOCKS_PER_SEC;
    printf("\n  SUMMARY: analyzed %zu component(s) in %.2fs; %zu component(s) contain "
           "bridges; %zu handle loop(s) (= bridges) total.\n",
           analyzed, secs, comps_with_bridge, total_bridge_faces);

    if (dump) {
        if (ObjIO_write_per_vertex_color(dump, V, nv, F, nf, col) == 0)
            printf("  wrote colored dump -> %s\n"
                   "   (gray=body, blue=side-facing; each HANDLE LOOP painted a "
                   "distinct bright color = the bridges)\n", dump);
        else
            fprintf(stderr, "  ERROR: cannot write dump %s\n", dump);
    }

    /* --sever: cut the SHORT handles (thin BPA bridges, loop < --max-len) on the
     * largest component and write the opened mesh. Verifies genus + manifold. */
    if (sever_out && ncomp > 0) {
        size_t cnf0 = (size_t)(coff[1] - coff[0]);
        SubMesh s; build_submesh(arena, V, nv, F, &cfaces[coff[0]], cnf0, &s);
        MeshTopoInfo tb; MeshTopo_analyze(arena, s.verts, s.nv, s.faces, s.nf, &tb);

        float *ov = NULL; int32_t *of = NULL; size_t onv = 0, onf = 0; long cut = 0;
        int rc = SeamCut_sever_short_handles(arena, s.verts, s.nv, s.faces, s.nf,
                                             sever_max_len, &ov, &onv, &of, &onf, &cut);
        if (rc == 0) {
            MeshTopoInfo ta; MeshTopo_analyze(arena, ov, onv, of, onf, &ta);
            printf("\n  SEVER (largest comp, short handles < %.0f vox):\n"
                   "    genus %.0f -> %.0f, boundary loops %ld -> %ld, nm_edges %ld -> %ld, "
                   "%ld handle(s) severed, V %zu -> %zu\n",
                   sever_max_len, tb.genus, ta.genus, tb.n_boundary_loops,
                   ta.n_boundary_loops, tb.n_nonmanifold_edges, ta.n_nonmanifold_edges,
                   cut, s.nv, onv);
            if (ObjIO_write(sever_out, ov, onv, of, onf) == 0)
                printf("    wrote severed mesh -> %s\n", sever_out);
            else
                fprintf(stderr, "    ERROR: cannot write %s\n", sever_out);
        } else {
            fprintf(stderr, "  SEVER failed (rc=%d)\n", rc);
        }
    }

    Arena_dispose(&arena);
    return 0;
}
