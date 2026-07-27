/*
 * api_smoke.c — direct-link exercise of the public C API.
 *
 * Complements tests/api_dlopen.c (the runtime-loading consumer test) with
 * broader op coverage on synthetic fixtures: audit, cleanup, hole fill,
 * decimation, fairing, orientation, splits, handle loops, MLS, BPA, the
 * volume pipeline, timeouts, and the xyz round trip. Builds fixtures
 * in-code; no data files. Exit code = number of failed checks.
 */
#include "scrollfiesta.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_fails = 0;

#define CHECK(cond, msg)                                                     \
    do {                                                                     \
        if (!(cond)) {                                                       \
            fprintf(stderr, "FAIL: %s (%s:%d)\n", msg, __FILE__, __LINE__);  \
            g_fails++;                                                       \
        }                                                                    \
    } while (0)

/* (n x n) flat grid sheet at z = z0, spacing `step`, two CCW tris per cell. */
static void make_sheet(int n, float step, float z0, sf_mesh *m)
{
    size_t nv = (size_t)n * (size_t)n;
    size_t nf = 2u * (size_t)(n - 1) * (size_t)(n - 1);
    memset(m, 0, sizeof *m);
    m->vertices = malloc(nv * 3 * sizeof(float));
    m->faces = malloc(nf * 3 * sizeof(int32_t));
    m->n_vertices = nv;
    m->n_faces = nf;
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            size_t v = (size_t)j * n + i;
            m->vertices[v * 3 + 0] = i * step;
            m->vertices[v * 3 + 1] = j * step;
            m->vertices[v * 3 + 2] = z0;
        }
    }
    size_t f = 0;
    for (int j = 0; j + 1 < n; j++) {
        for (int i = 0; i + 1 < n; i++) {
            int32_t v00 = j * n + i, v10 = j * n + i + 1;
            int32_t v01 = (j + 1) * n + i, v11 = (j + 1) * n + i + 1;
            m->faces[f * 3 + 0] = v00; m->faces[f * 3 + 1] = v10;
            m->faces[f * 3 + 2] = v01; f++;
            m->faces[f * 3 + 0] = v10; m->faces[f * 3 + 1] = v11;
            m->faces[f * 3 + 2] = v01; f++;
        }
    }
}

static void free_mesh_caller(sf_mesh *m)
{
    free(m->vertices);
    free(m->faces);
    memset(m, 0, sizeof *m);
}

static int cancel_immediately(void *user, const char *stage, double fraction)
{
    (void)user; (void)stage; (void)fraction;
    return 1;
}

int main(void)
{
    printf("scrollfiesta %s api smoke\n", sf_version_string());
    CHECK(sf_get_api(SCROLLFIESTA_ABI_VERSION) != NULL, "sf_get_api matches");
    CHECK(sf_get_api(SCROLLFIESTA_ABI_VERSION + 100) == NULL, "wrong ABI NULL");

    sf_mesh sheet;
    make_sheet(10, 1.f, 0.f, &sheet);

    /* ── audit ─────────────────────────────────────────────────────────── */
    sf_topology_report topo;
    CHECK(sf_topology_audit(&sheet, NULL, &topo) == SF_OK, "audit rc");
    CHECK(topo.n_components == 1 && topo.is_disk == 1, "sheet is a disk");
    CHECK(topo.n_nonmanifold_edges == 0, "no NM edges");

    /* ── hole fill: knock one interior vertex's faces out ─────────────── */
    {
        /* Rebuild the sheet without any face touching vertex (5,5). */
        sf_mesh holed;
        make_sheet(10, 1.f, 0.f, &holed);
        size_t keep = 0;
        const int32_t center = 5 * 10 + 5;
        for (size_t f = 0; f < holed.n_faces; f++) {
            if (holed.faces[f * 3] == center || holed.faces[f * 3 + 1] == center ||
                holed.faces[f * 3 + 2] == center)
                continue;
            memcpy(&holed.faces[keep * 3], &holed.faces[f * 3],
                   3 * sizeof(int32_t));
            keep++;
        }
        holed.n_faces = keep;

        sf_topology_report before;
        CHECK(sf_topology_audit(&holed, NULL, &before) == SF_OK, "holed audit");
        CHECK(before.n_boundary_loops == 2, "hole adds a boundary loop");

        sf_holefill_config cfg = sf_holefill_config_default();
        sf_mesh filled;
        sf_holefill_report rep;
        CHECK(sf_fill_holes(&holed, &cfg, &filled, &rep) == SF_OK, "fill rc");
        CHECK(rep.n_filled >= 1, "one hole filled");
        sf_topology_report after;
        CHECK(sf_topology_audit(&filled, NULL, &after) == SF_OK, "post audit");
        CHECK(after.n_boundary_loops == 1, "hole closed");
        /* originals keep identity provenance */
        CHECK(filled.vmap != NULL, "fill emits vmap");
        for (size_t i = 0; filled.vmap && i < holed.n_vertices; i++)
            if (filled.vmap[i] != (int32_t)i) {
                CHECK(0, "fill vmap identity prefix");
                break;
            }
        sf_mesh_free(&filled);
        free_mesh_caller(&holed);
    }

    fprintf(stderr, "== stage decimate ==\n"); fflush(stderr);
    /* ── decimate ──────────────────────────────────────────────────────── */
    {
        sf_mesh big;
        make_sheet(40, 1.f, 0.f, &big);
        sf_decimate_config cfg = sf_decimate_config_default();
        cfg.target_ratio = 0.25f;
        sf_mesh out;
        sf_decimate_report rep;
        CHECK(sf_decimate(&big, &cfg, &out, &rep) == SF_OK, "decimate rc");
        CHECK(rep.faces_out < rep.faces_in, "decimate reduced");
        CHECK(out.vmap == NULL, "decimate has no provenance by contract");
        sf_mesh_free(&out);
        free_mesh_caller(&big);
    }

    fprintf(stderr, "== stage fair ==\n"); fflush(stderr);
    /* ── fairing (few iters) keeps counts, moves nothing on a flat sheet ─ */
    {
        sf_fair_config cfg = sf_fair_config_default();
        cfg.max_iters = 5;
        sf_mesh out;
        sf_fair_report rep;
        CHECK(sf_fair_developable(&sheet, &cfg, &out, &rep) == SF_OK, "fair rc");
        CHECK(out.n_vertices == sheet.n_vertices, "fair keeps verts");
        CHECK(out.vmap && out.vmap[3] == 3, "fair identity vmap");
        sf_mesh_free(&out);
    }

    fprintf(stderr, "== stage orient ==\n"); fflush(stderr);
    /* ── orient: flip half the faces, expect consistency restored ─────── */
    {
        sf_mesh flipped;
        make_sheet(10, 1.f, 0.f, &flipped);
        for (size_t f = 0; f < flipped.n_faces; f += 2) {
            int32_t t = flipped.faces[f * 3 + 1];
            flipped.faces[f * 3 + 1] = flipped.faces[f * 3 + 2];
            flipped.faces[f * 3 + 2] = t;
        }
        sf_orient_config cfg = sf_orient_config_default();
        sf_mesh out;
        sf_orient_report rep;
        CHECK(sf_orient(&flipped, &cfg, &out, &rep) == SF_OK, "orient rc");
        CHECK(rep.residual_same_dir == 0, "orient leaves no same-dir edges");
        sf_mesh_free(&out);
        free_mesh_caller(&flipped);
    }

    fprintf(stderr, "== stage peel ==\n"); fflush(stderr);
    /* ── depth peel: two stacked sheets joined by one wall face split ──── */
    {
        sf_mesh a, b;
        /* Sheets must be wide relative to their separation so the component
         * PCA normal is the stacking axis (in-plane variance ~21 vs z
         * variance 2.25 here); the wall edge then jumps 3 > min_gap 1.5. */
        make_sheet(16, 1.f, 0.f, &a);
        make_sheet(16, 1.f, 3.f, &b);
        sf_mesh joined;
        memset(&joined, 0, sizeof joined);
        joined.n_vertices = a.n_vertices + b.n_vertices;
        joined.n_faces = a.n_faces + b.n_faces + 1;
        joined.vertices = malloc(joined.n_vertices * 3 * sizeof(float));
        joined.faces = malloc(joined.n_faces * 3 * sizeof(int32_t));
        memcpy(joined.vertices, a.vertices, a.n_vertices * 3 * sizeof(float));
        memcpy(joined.vertices + a.n_vertices * 3, b.vertices,
               b.n_vertices * 3 * sizeof(float));
        memcpy(joined.faces, a.faces, a.n_faces * 3 * sizeof(int32_t));
        for (size_t f = 0; f < b.n_faces * 3; f++)
            joined.faces[a.n_faces * 3 + f] =
                b.faces[f] + (int32_t)a.n_vertices;
        /* wall face bridging the two layers */
        joined.faces[(joined.n_faces - 1) * 3 + 0] = 0;
        joined.faces[(joined.n_faces - 1) * 3 + 1] = 1;
        joined.faces[(joined.n_faces - 1) * 3 + 2] = (int32_t)a.n_vertices;

        sf_split_peel_config cfg = sf_split_peel_config_default();
        cfg.min_comp_verts = 8;
        sf_mesh_list pieces;
        sf_split_report rep;
        CHECK(sf_split_depth_peel(&joined, &cfg, &pieces, &rep) == SF_OK,
              "peel rc");
        CHECK(pieces.count == 2, "peel separates the two layers");
        for (size_t i = 0; i < pieces.count; i++)
            CHECK(pieces.items[i].vmap != NULL, "peel pieces carry vmap");
        sf_mesh_list_free(&pieces);
        free_mesh_caller(&joined);
        free_mesh_caller(&a);
        free_mesh_caller(&b);
    }

    fprintf(stderr, "== stage mls ==\n"); fflush(stderr);
    /* ── MLS smooths a noisy cloud toward the plane ────────────────────── */
    {
        const int n = 24;
        size_t np = (size_t)n * n;
        float *pts = malloc(np * 3 * sizeof(float));
        unsigned seed = 12345u;
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) {
                seed = seed * 1664525u + 1013904223u;
                float noise = ((int)((seed >> 16) % 200) - 100) / 400.f; /* +-0.25 */
                size_t v = (size_t)j * n + i;
                pts[v * 3 + 0] = (float)i;
                pts[v * 3 + 1] = (float)j;
                pts[v * 3 + 2] = noise;
            }
        float *proj = NULL, *nrm = NULL;
        CHECK(sf_mls_project(pts, np, 4.f, NULL, &proj, &nrm) == SF_OK,
              "mls rc");
        double before = 0, after = 0;
        for (size_t i = 0; i < np; i++) {
            before += fabs(pts[i * 3 + 2]);
            after += fabs(proj[i * 3 + 2]);
        }
        CHECK(after < before, "mls reduced through-thickness noise");
        sf_free(proj);
        sf_free(nrm);
        free(pts);
    }

    fprintf(stderr, "== stage bpa ==\n"); fflush(stderr);
    /* ── BPA reconstructs a plane from points + normals ────────────────── */
    {
        /* Pipeline-realistic cloud: ~0.77 vox spacing (the voxel-shell
         * density BPA's default rho 1.2 is tuned for), slight jitter, and a
         * gentle z curve. NB an EXACTLY coplanar cloud (all z identical)
         * crashes BallPivot_reconstruct -- a pre-existing robustness bug on
         * degenerate input the pipeline never produces; tracked separately,
         * deliberately not exercised by this smoke test. */
        const int n = 16;
        const float step = 0.77f;
        size_t np = (size_t)n * n;
        float *pts = malloc(np * 3 * sizeof(float));
        float *nrm = malloc(np * 3 * sizeof(float));
        unsigned rng = 7u;
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) {
                size_t v = (size_t)j * n + i;
                rng = rng * 1664525u + 1013904223u;
                float jx = ((int)((rng >> 16) % 100) - 50) / 2000.f; /* +-0.025 */
                rng = rng * 1664525u + 1013904223u;
                float jy = ((int)((rng >> 16) % 100) - 50) / 2000.f;
                pts[v * 3 + 0] = i * step + jx;
                pts[v * 3 + 1] = j * step + jy;
                pts[v * 3 + 2] = 0.15f * sinf(i * 0.4f);
                nrm[v * 3 + 0] = 0.f;
                nrm[v * 3 + 1] = 0.f;
                nrm[v * 3 + 2] = 1.f;
            }
        sf_bpa_config cfg = sf_bpa_config_default();
        sf_mesh out;
        sf_bpa_report rep;
        CHECK(sf_reconstruct_bpa(pts, nrm, np, &cfg, &out, &rep) == SF_OK,
              "bpa rc");
        CHECK(out.n_faces > np, "bpa triangulated the grid");
        sf_topology_report t2;
        CHECK(sf_topology_audit(&out, NULL, &t2) == SF_OK, "bpa audit");
        CHECK(t2.n_nonmanifold_edges == 0, "bpa manifold");
        sf_mesh_free(&out);
        free(pts);
        free(nrm);
    }

    fprintf(stderr, "== stage pipeline ==\n"); fflush(stderr);
    /* ── volume pipeline on a synthetic curved sheet ───────────────────── */
    {
        const int cube = 64, halo = 0, p = cube + 2 * halo;
        uint8_t *vol = calloc((size_t)p * p * p, 1);
        /* thick-ish sheet: z = 24 + 6*sin(x/10), 3 voxels thick */
        for (int y = 4; y < p - 4; y++)
            for (int x = 4; x < p - 4; x++) {
                int z0 = 24 + (int)(6.0 * sin(x / 10.0));
                for (int dz = 0; dz < 3; dz++) {
                    int z = z0 + dz;
                    if (z >= 0 && z < p)
                        vol[(size_t)z * p * p + (size_t)y * p + x] = 1;
                }
            }
        sf_volume v = { vol, p, p, p };
        sf_pipeline_config cfg = sf_pipeline_config_default();
        cfg.cube_size = cube;
        cfg.halo_voxels = halo;
        cfg.skip_qem = 1;   /* keep the smoke test fast */
        sf_mesh_list full, trimmed;
        sf_pipeline_report rep;
        sf_status rc = sf_pipeline_run(&v, &cfg, &full, &trimmed, &rep);
        CHECK(rc == SF_OK, "pipeline rc");
        if (rc == SF_OK) {
            CHECK(full.count >= 1, "pipeline produced a mesh");
            size_t total_faces = 0;
            for (size_t i = 0; i < full.count; i++)
                total_faces += full.items[i].n_faces;
            CHECK(total_faces > 100, "pipeline mesh is non-trivial");
            sf_mesh_list_free(&full);
            sf_mesh_list_free(&trimmed);
        }
        free(vol);
    }

    fprintf(stderr, "== stage weld ==\n"); fflush(stderr);
    /* ── weld: two sheets astride the x=128 cube boundary bridge into one
     *    component ─────────────────────────────────────────────────────── */
    {
        /* Pipeline-like sheets: 0.77 vox spacing, gentle z curve, ending
         * ~0.5 vox short of the seam plane on each side (gap ~1 vox, well
         * inside the bridge ball's reach). */
        const int n = 20;
        const float step = 0.77f;
        sf_mesh sides[2];
        for (int sidx = 0; sidx < 2; sidx++) {
            const float x0 = sidx == 0 ? 128.f - 0.5f - (n - 1) * step
                                       : 128.f + 0.5f;
            sf_mesh *m = &sides[sidx];
            make_sheet(n, step, 0.f, m);
            for (size_t v = 0; v < m->n_vertices; v++) {
                m->vertices[v * 3 + 0] += x0;
                m->vertices[v * 3 + 1] += 100.f;
                m->vertices[v * 3 + 2] =
                    40.f + 0.2f * sinf(m->vertices[v * 3 + 0] * 0.3f);
            }
        }

        sf_topology_report t_pre;
        {
            /* Sanity: concatenated but unwelded = 2 components. */
            sf_mesh cat;
            memset(&cat, 0, sizeof cat);
            cat.n_vertices = sides[0].n_vertices + sides[1].n_vertices;
            cat.n_faces = sides[0].n_faces + sides[1].n_faces;
            cat.vertices = malloc(cat.n_vertices * 3 * sizeof(float));
            cat.faces = malloc(cat.n_faces * 3 * sizeof(int32_t));
            memcpy(cat.vertices, sides[0].vertices,
                   sides[0].n_vertices * 3 * sizeof(float));
            memcpy(cat.vertices + sides[0].n_vertices * 3, sides[1].vertices,
                   sides[1].n_vertices * 3 * sizeof(float));
            memcpy(cat.faces, sides[0].faces,
                   sides[0].n_faces * 3 * sizeof(int32_t));
            for (size_t f = 0; f < sides[1].n_faces * 3; f++)
                cat.faces[sides[0].n_faces * 3 + f] =
                    sides[1].faces[f] + (int32_t)sides[0].n_vertices;
            CHECK(sf_topology_audit(&cat, NULL, &t_pre) == SF_OK,
                  "pre-weld audit");
            CHECK(t_pre.n_components == 2, "two components before weld");
            free_mesh_caller(&cat);
        }

        sf_weld_config wcfg = sf_weld_config_default();
        sf_mesh welded;
        sf_weld_report wrep;
        sf_status wrc = sf_weld(sides, 2, &wcfg, &welded, &wrep);
        CHECK(wrc == SF_OK, "weld rc");
        if (wrc == SF_OK) {
            CHECK(wrep.bridge_faces > 0, "weld bridged the seam");
            sf_topology_report t_post;
            CHECK(sf_topology_audit(&welded, NULL, &t_post) == SF_OK,
                  "post-weld audit");
            CHECK(t_post.n_components == 1, "one component after weld");
            CHECK(t_post.n_nonmanifold_edges == 0, "weld stays manifold");
            sf_mesh_free(&welded);
        }
        free_mesh_caller(&sides[0]);
        free_mesh_caller(&sides[1]);
    }

    fprintf(stderr, "== stage cancel ==\n"); fflush(stderr);
    /* ── cancellation + timeout ────────────────────────────────────────── */
    {
        sf_cleanup_config cfg = sf_cleanup_config_default();
        cfg.common.progress = cancel_immediately;
        sf_mesh out;
        CHECK(sf_cleanup(&sheet, &cfg, &out, NULL) == SF_CANCELLED,
              "cancel via callback");
        CHECK(out.vertices == NULL, "cancelled output empty");

        sf_fair_config fcfg = sf_fair_config_default();
        fcfg.max_iters = 100000;
        fcfg.tol_grad = 0.0;
        fcfg.tol_e = 0.0;
        fcfg.common.timeout_sec = 0.02;
        sf_mesh out2;
        sf_status rc = sf_fair_developable(&sheet, &fcfg, &out2, NULL);
        /* tiny fixture may converge before the deadline fires; both are fine,
         * but a timeout must never yield a partially-filled output */
        CHECK(rc == SF_ERROR_TIMEOUT || rc == SF_OK, "timeout path sane");
        if (rc == SF_ERROR_TIMEOUT)
            CHECK(out2.vertices == NULL, "timed-out output empty");
        else
            sf_mesh_free(&out2);
    }

    fprintf(stderr, "== stage obj ==\n"); fflush(stderr);
    /* ── OBJ round trip (true-xyz convention) ──────────────────────────── */
    {
        const char *path = "api_smoke_roundtrip.obj";
        CHECK(sf_mesh_save_obj(path, &sheet) == SF_OK, "obj save");
        sf_mesh back;
        CHECK(sf_mesh_load_obj(path, &back) == SF_OK, "obj load");
        CHECK(back.n_vertices == sheet.n_vertices, "obj verts");
        CHECK(back.n_faces == sheet.n_faces, "obj faces");
        sf_mesh_free(&back);
        remove(path);
    }

    fprintf(stderr, "== stage selftest ==\n"); fflush(stderr);
    /* ── module self-tests ─────────────────────────────────────────────── */
    CHECK(sf_selftest() == 0, "module selftests");

    free_mesh_caller(&sheet);

    if (g_fails == 0) {
        printf("api_smoke: all checks passed\n");
        return 0;
    }
    fprintf(stderr, "api_smoke: %d check(s) FAILED\n", g_fails);
    return g_fails;
}
