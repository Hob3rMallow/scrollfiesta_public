/*
 * api_dlopen.c — the canonical consumer test.
 *
 * Links NOTHING from ScrollFiesta. Loads the shared library at runtime,
 * resolves the single "sf_get_api" symbol, and exercises the function-
 * pointer table — exactly what an embedding host (VC3D) does. Usage:
 *
 *     api_dlopen <path-to-scrollfiesta-shared-library>
 *
 * Exit code 0 = all checks pass.
 */
#include "scrollfiesta.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
  #ifndef WIN32_LEAN_AND_MEAN
    #define WIN32_LEAN_AND_MEAN
  #endif
  #include <windows.h>
  typedef HMODULE sf_lib_handle;
  static sf_lib_handle lib_open(const char *path) { return LoadLibraryA(path); }
  static void *lib_sym(sf_lib_handle h, const char *name)
  { return (void *)GetProcAddress(h, name); }
#else
  #include <dlfcn.h>
  typedef void *sf_lib_handle;
  static sf_lib_handle lib_open(const char *path) { return dlopen(path, RTLD_NOW); }
  static void *lib_sym(sf_lib_handle h, const char *name) { return dlsym(h, name); }
#endif

static int g_fails = 0;

static int cancel_immediately(void *user, const char *stage, double fraction)
{
    (void)user; (void)stage; (void)fraction;
    return 1;   /* request cancellation */
}

#define CHECK(cond, msg)                                                     \
    do {                                                                     \
        if (!(cond)) {                                                       \
            fprintf(stderr, "FAIL: %s (%s:%d)\n", msg, __FILE__, __LINE__);  \
            g_fails++;                                                       \
        }                                                                    \
    } while (0)

/* Build an (n x n)-vertex flat grid sheet in the z=0 plane: vertices at
 * (i, j, 0), two CCW triangles per cell. A clean manifold disk. */
static void make_grid_sheet(int n, sf_mesh *m)
{
    size_t nv = (size_t)n * (size_t)n;
    size_t nf = 2u * (size_t)(n - 1) * (size_t)(n - 1);
    memset(m, 0, sizeof *m);
    m->vertices = malloc(nv * 3 * sizeof(float));
    m->faces    = malloc(nf * 3 * sizeof(int32_t));
    m->n_vertices = nv;
    m->n_faces    = nf;
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            size_t v = (size_t)j * (size_t)n + (size_t)i;
            m->vertices[v * 3 + 0] = (float)i;
            m->vertices[v * 3 + 1] = (float)j;
            m->vertices[v * 3 + 2] = 0.0f;
        }
    }
    size_t f = 0;
    for (int j = 0; j + 1 < n; j++) {
        for (int i = 0; i + 1 < n; i++) {
            int32_t v00 = j * n + i;
            int32_t v10 = j * n + i + 1;
            int32_t v01 = (j + 1) * n + i;
            int32_t v11 = (j + 1) * n + i + 1;
            m->faces[f * 3 + 0] = v00;
            m->faces[f * 3 + 1] = v10;
            m->faces[f * 3 + 2] = v01;
            f++;
            m->faces[f * 3 + 0] = v10;
            m->faces[f * 3 + 1] = v11;
            m->faces[f * 3 + 2] = v01;
            f++;
        }
    }
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "usage: %s <scrollfiesta shared library>\n", argv[0]);
        return 2;
    }

    sf_lib_handle h = lib_open(argv[1]);
    if (!h) {
        fprintf(stderr, "FAIL: could not load '%s'\n", argv[1]);
        return 1;
    }

    sf_get_api_fn get_api = (sf_get_api_fn)lib_sym(h, "sf_get_api");
    CHECK(get_api != NULL, "resolve sf_get_api");
    if (!get_api)
        return 1;

    /* Wrong ABI must be refused. */
    CHECK(get_api(0xdeadbeefu) == NULL, "wrong ABI returns NULL");

    const sf_api *sf = get_api(SCROLLFIESTA_ABI_VERSION);
    CHECK(sf != NULL, "matching ABI returns the table");
    if (!sf)
        return 1;
    CHECK(sf->abi_version == SCROLLFIESTA_ABI_VERSION, "table abi_version");
    CHECK(sf->struct_size == (uint32_t)sizeof(sf_api),
          "table struct_size matches this header");

    int maj = -1, min = -1, pat = -1;
    sf->version(&maj, &min, &pat);
    printf("loaded scrollfiesta %s (%d.%d.%d, abi %u)\n",
           sf->version_string(), maj, min, pat, sf->abi_version);

    /* ── topology audit on a clean sheet ─────────────────────────────── */
    sf_mesh sheet;
    make_grid_sheet(8, &sheet);   /* caller-owned buffers */

    sf_topology_report topo;
    sf_status rc = sf->topology_audit(&sheet, NULL, &topo);
    CHECK(rc == SF_OK, "topology_audit rc");
    CHECK(topo.n_verts == 64, "audit vert count");
    CHECK(topo.n_faces == 98, "audit face count");
    CHECK(topo.n_components == 1, "audit one component");
    CHECK(topo.n_nonmanifold_edges == 0, "audit no non-manifold edges");
    CHECK(topo.n_boundary_loops == 1, "audit one boundary loop");
    CHECK(topo.is_disk == 1, "audit is a disk");

    /* ── cleanup round-trip: a clean sheet must come back intact with
     *    identity provenance ─────────────────────────────────────────── */
    sf_cleanup_config ccfg = sf->cleanup_config_default();
    sf_mesh           cleaned;
    sf_cleanup_report crep;
    rc = sf->cleanup(&sheet, &ccfg, &cleaned, &crep);
    CHECK(rc == SF_OK, "cleanup rc");
    if (rc == SF_OK) {
        CHECK(cleaned.n_vertices == sheet.n_vertices, "cleanup keeps verts");
        CHECK(cleaned.n_faces == sheet.n_faces, "cleanup keeps faces");
        CHECK(cleaned.vmap != NULL, "cleanup emits provenance");
        int identity = 1;
        int positions_intact = 1;
        for (size_t i = 0; cleaned.vmap && i < cleaned.n_vertices; i++) {
            if (cleaned.vmap[i] != (int32_t)i)
                identity = 0;
            if (cleaned.vmap[i] >= 0) {
                const float *a = &cleaned.vertices[i * 3];
                const float *b = &sheet.vertices[(size_t)cleaned.vmap[i] * 3];
                if (memcmp(a, b, 3 * sizeof(float)) != 0)
                    positions_intact = 0;
            }
        }
        CHECK(identity, "cleanup provenance is identity on a clean sheet");
        CHECK(positions_intact, "provenanced positions bit-identical");
        sf->mesh_free(&cleaned);
    }

    /* ── depth-peel self-gate: one sheet stays one piece ─────────────── */
    sf_split_peel_config pcfg = sf->split_peel_config_default();
    pcfg.min_comp_verts = 4;   /* tiny test mesh */
    sf_mesh_list    pieces;
    sf_split_report srep;
    rc = sf->split_depth_peel(&sheet, &pcfg, &pieces, &srep);
    CHECK(rc == SF_OK, "split_depth_peel rc");
    if (rc == SF_OK) {
        CHECK(pieces.count == 1, "flat sheet self-gates to one piece");
        CHECK(srep.self_gated == 1, "self_gated flag");
        sf->mesh_list_free(&pieces);
    }

    /* ── cancellation: a callback that cancels at the first report ───── */
    {
        sf_cleanup_config cc = sf->cleanup_config_default();
        cc.common.progress      = cancel_immediately;
        cc.common.progress_user = NULL;
        sf_mesh           cm2;
        sf_cleanup_report cr2;
        sf_status rc2 = sf->cleanup(&sheet, &cc, &cm2, &cr2);
        CHECK(rc2 == SF_CANCELLED, "cancel via progress callback");
        CHECK(cm2.vertices == NULL, "cancelled op leaves output empty");
    }

    /* ── weld is experimental: table entry may be NULL — probe pattern ─ */
    if (sf->weld == NULL)
        printf("weld: not available in this build (expected for 0.9)\n");

    /* ── error-string sanity ──────────────────────────────────────────── */
    CHECK(strcmp(sf->status_str(SF_OK), "ok") == 0, "status_str(SF_OK)");
    CHECK(sf->status_str(SF_CANCELLED) != NULL, "status_str(SF_CANCELLED)");

    free(sheet.vertices);
    free(sheet.faces);

    if (g_fails == 0) {
        printf("api_dlopen: all checks passed\n");
        return 0;
    }
    fprintf(stderr, "api_dlopen: %d check(s) FAILED\n", g_fails);
    return 1;
}
