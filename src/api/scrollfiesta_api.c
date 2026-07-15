/*
 * scrollfiesta_api.c — implementation of the public C API (include/
 * scrollfiesta.h): argument validation, the per-op guard (mutex + arena +
 * run context + exception containment), the op wrappers over the internal
 * modules, config defaults, and the sf_get_api function-pointer table.
 *
 * Wrapper discipline (see sf_internal.h): everything that can RAISE (arena
 * allocations, module calls, cancellation polls) happens first; the final
 * malloc-backed copy-out never raises, so outputs are either untouched or
 * complete. No longjmp ever crosses the API boundary.
 */
#include "sf_internal.h"

#include "../common/mls_project.h"
#include "../common/obj_io.h"
#include "../common/pipeline_constants.h"
#include "../common/qem.h"
#include "../holefill/hole_fill.h"
#include "../remesh/ball_pivot.h"
#include "../remesh/component_cull.h"
#include "../remesh/manifold_guard.h"
#include "../remesh/orient_mesh.h"
#include "../remesh/orient_weld.h"
#include "../remesh/pinhole_fill.h"
#include "../remesh/weld_cleanup.h"
#include "../split/bridge_cut.h"
#include "../split/depth_peel.h"
#include "../split/dev_cut.h"
#include "../split/oracle.h"
#include "../split/overlap_sep.h"
#include "../topology/developability.h"
#include "../topology/mesh_topo.h"
#include "../topology/seam_cut.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

/* ══════════════════════════════════════════════════════════════════════════
 * Guard
 * ════════════════════════════════════════════════════════════════════════ */

ves_mutex_t g_sf_op_mutex = VES_MUTEX_INITIALIZER;

sf_status sf_op_begin(SfOpGuard *g, const sf_common_opts *opts)
{
    volatile sf_status rc = SF_OK;

    memset(g, 0, sizeof *g);
    ves_mutex_lock(&g_sf_op_mutex);
    g->locked = 1;

    if (opts) {
        g->ctx.progress      = opts->progress;
        g->ctx.progress_user = opts->progress_user;
        g->ctx.log           = opts->log;
        g->ctx.log_user      = opts->log_user;
        /* sf_kv and SfRunKV are layout-identical {const char*, const char*}
         * pairs; the public header cannot include run_ctx.h. */
        g->ctx.tuning        = (const SfRunKV *)opts->tuning;
        g->ctx.n_tuning      = opts->n_tuning;
        if (opts->timeout_sec > 0.0)
            g->ctx.deadline_sec = ves_clock_sec() + opts->timeout_sec;
        ves_omp_set_threads(opts->n_threads);
    }

    TRY
        g->arena = Arena_new();
    EXCEPT(Arena_Failed)
        rc = SF_ERROR_OOM;
    END_TRY;

    if (rc != SF_OK) {
        ves_mutex_unlock(&g_sf_op_mutex);
        g->locked = 0;
        return rc;
    }
    RunCtx_push(&g->ctx);
    return SF_OK;
}

void sf_op_end(SfOpGuard *g)
{
    RunCtx_pop();
    if (g->arena)
        Arena_dispose(&g->arena);
    if (g->locked)
        ves_mutex_unlock(&g_sf_op_mutex);
    g->locked = 0;
}

static sf_status validate_in(const sf_mesh *in)
{
    return sf_mesh_validate(in, NULL, 0);
}

/* Raw point arrays (BPA/MLS inputs) get the same finiteness guarantee as
 * meshes: a single garbage coordinate would otherwise blow up the spatial
 * grids downstream instead of failing cleanly. */
static int points_finite(const float *xyz, size_t n)
{
    for (size_t i = 0; i < n * 3; i++) {
        if (!isfinite(xyz[i]))
            return 0;
    }
    return 1;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Version / status / defaults
 * ════════════════════════════════════════════════════════════════════════ */

SF_API void sf_version(int *major, int *minor, int *patch)
{
    if (major) *major = SCROLLFIESTA_VERSION_MAJOR;
    if (minor) *minor = SCROLLFIESTA_VERSION_MINOR;
    if (patch) *patch = SCROLLFIESTA_VERSION_PATCH;
}

SF_API const char *sf_version_string(void)
{
    return "0.9.0";
}

SF_API const char *sf_status_str(sf_status s)
{
    switch (s) {
    case SF_OK:                return "ok";
    case SF_ERROR:             return "internal failure";
    case SF_ERROR_BAD_ARG:     return "invalid argument";
    case SF_ERROR_OOM:         return "out of memory";
    case SF_ERROR_IO:          return "I/O failure";
    case SF_ERROR_TIMEOUT:     return "timed out";
    case SF_CANCELLED:         return "cancelled";
    case SF_ERROR_NO_TIFF:     return "built without TIFF support";
    case SF_ERROR_ABI:         return "ABI mismatch";
    case SF_ERROR_UNSUPPORTED: return "not supported in this build";
    }
    return "unknown status";
}

SF_API sf_common_opts sf_common_opts_default(void)
{
    sf_common_opts o;
    memset(&o, 0, sizeof o);
    return o;
}

SF_API sf_cleanup_config sf_cleanup_config_default(void)
{
    sf_cleanup_config c;
    memset(&c, 0, sizeof c);
    c.manifold_repair    = 1;
    c.reorient           = 1;
    c.fill_pinholes      = 1;
    c.sliver_cleanup     = 0;
    c.sliver_min_alt     = WELD_CLEANUP_SLIVER_MIN_ALT;
    c.degen_area         = WELD_CLEANUP_DEGEN_AREA;
    c.max_collapse_len   = WELD_CLEANUP_MAX_COLLAPSE_LEN;
    c.cull_min_area_frac = 0.0f;   /* pipeline uses KIBBLE_AREA_FRAC = 0.02 */
    c.respect_pins       = 1;
    return c;
}

SF_API sf_holefill_config sf_holefill_config_default(void)
{
    sf_holefill_config c;
    memset(&c, 0, sizeof c);
    c.interior_only = 1;
    return c;
}

SF_API sf_decimate_config sf_decimate_config_default(void)
{
    sf_decimate_config c;
    memset(&c, 0, sizeof c);
    c.target_ratio = QEM_TARGET_RATIO;
    c.respect_pins = 1;
    return c;
}

SF_API sf_fair_config sf_fair_config_default(void)
{
    sf_fair_config c;
    memset(&c, 0, sizeof c);
    c.max_iters    = 100;
    c.tol_grad     = 1e-6;
    c.tol_e        = 1e-6;
    c.max_disp_vox = 1.0;
    c.pin_boundary = 1;
    return c;
}

SF_API sf_orient_config sf_orient_config_default(void)
{
    sf_orient_config c;
    memset(&c, 0, sizeof c);
    c.align_components = 1;
    c.align_radius     = 3.0f;
    return c;
}

SF_API sf_split_peel_config sf_split_peel_config_default(void)
{
    sf_split_peel_config c;
    memset(&c, 0, sizeof c);
    c.min_gap        = DEPTH_PEEL_MIN_GAP_VOX;
    c.gap_depth      = DEPTH_PEEL_GAP_DEPTH;
    c.min_comp_verts = DEPTH_PEEL_MIN_COMP_VERTS;
    return c;
}

SF_API sf_split_dev_config sf_split_dev_config_default(void)
{
    sf_split_dev_config c;
    memset(&c, 0, sizeof c);
    c.eps            = DEV_CUT_EPS;
    c.gap_depth      = DEV_CUT_GAP_DEPTH;
    c.min_comp_verts = DEV_CUT_MIN_COMP_VERTS;
    c.apply_gate     = 1;
    c.gate_margin    = DEV_CUT_GATE_MARGIN;
    c.gate_abs_ceil  = DEV_CUT_GATE_ABS_CEIL;
    return c;
}

SF_API sf_split_bridge_config sf_split_bridge_config_default(void)
{
    sf_split_bridge_config c;
    memset(&c, 0, sizeof c);
    c.oracle_grid_size = ORACLE_UV_GRID_SIZE;
    c.component_timeout_sec = 30.0;
    return c;
}

SF_API sf_split_overlap_config sf_split_overlap_config_default(void)
{
    sf_split_overlap_config c;
    memset(&c, 0, sizeof c);
    c.sheet_count = 0;
    c.component_timeout_sec = 30.0;
    return c;
}

SF_API sf_detangle_config sf_detangle_config_default(void)
{
    sf_detangle_config c;
    memset(&c, 0, sizeof c);
    c.enable_peel    = 1;
    c.enable_dev     = 1;
    c.enable_bridge  = 1;
    c.enable_overlap = 1;
    c.min_comp_verts = RESPLIT_MIN_COMP_VERTS;
    c.component_timeout_sec = 30.0;
    return c;
}

SF_API sf_seamcut_config sf_seamcut_config_default(void)
{
    sf_seamcut_config c;
    memset(&c, 0, sizeof c);
    c.accept_frac = 0.0;
    c.r_max       = 5;
    return c;
}

SF_API sf_bpa_config sf_bpa_config_default(void)
{
    sf_bpa_config c;
    memset(&c, 0, sizeof c);
    c.rho = BPA_RHO_VOX;
    return c;
}

SF_API sf_pipeline_config sf_pipeline_config_default(void)
{
    sf_pipeline_config c;
    memset(&c, 0, sizeof c);
    c.cube_size         = 128;
    c.halo_voxels       = 13;
    c.enable_depth_peel = 1;
    c.enable_dev_cut    = 1;
    c.enable_sever      = 1;
    return c;
}

SF_API sf_weld_config sf_weld_config_default(void)
{
    sf_weld_config c;
    memset(&c, 0, sizeof c);
    c.cube_size  = 128.0f;
    c.rho        = BRIDGE_RHO_MIN;
    c.rho_max    = BRIDGE_RHO_MAX;
    c.weld_eps   = 1e-4f;
    c.fill_holes = 1;
    c.cleanup    = 1;
    return c;
}

/* ══════════════════════════════════════════════════════════════════════════
 * OBJ convenience I/O (true xyz — not the internal zyx tool format)
 * ════════════════════════════════════════════════════════════════════════ */

SF_API sf_status sf_mesh_load_obj(const char *path, sf_mesh *out)
{
    SfOpGuard          g;
    volatile sf_status rc;

    if (!path || !out)
        return SF_ERROR_BAD_ARG;
    memset(out, 0, sizeof *out);

    rc = sf_op_begin(&g, NULL);
    if (rc != SF_OK)
        return rc;

    TRY {
        float   *v  = NULL;
        int32_t *f  = NULL;
        size_t   nv = 0, nf = 0;
        if (ObjIO_read(g.arena, path, &v, &nv, &f, &nf) != 0 || nv == 0 || nf == 0)
            rc = SF_ERROR_IO;
        else
            rc = sf_mesh_out_raw(out, v, nv, f, nf, NULL, NULL, NULL);
    }
    SF_EXCEPT_TAIL(g, rc)
    END_TRY;

    sf_op_end(&g);
    return rc;
}

SF_API sf_status sf_mesh_save_obj(const char *path, const sf_mesh *m)
{
    if (!path || validate_in(m) != SF_OK)
        return SF_ERROR_BAD_ARG;
    if (ObjIO_write(path, m->vertices, m->n_vertices, m->faces, m->n_faces) != 0)
        return SF_ERROR_IO;
    return SF_OK;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Diagnostics
 * ════════════════════════════════════════════════════════════════════════ */

SF_API sf_status sf_topology_audit(const sf_mesh *in, const sf_common_opts *opts,
                                   sf_topology_report *out)
{
    SfOpGuard          g;
    volatile sf_status rc;

    if (!out || validate_in(in) != SF_OK)
        return SF_ERROR_BAD_ARG;
    memset(out, 0, sizeof *out);

    rc = sf_op_begin(&g, opts);
    if (rc != SF_OK)
        return rc;

    TRY {
        int32_t     *faces = sf_faces_swap(g.arena, in->faces, in->n_faces);
        MeshTopoInfo info;
        if (MeshTopo_analyze(g.arena, NULL, in->n_vertices, faces,
                             in->n_faces, &info) != 0) {
            rc = SF_ERROR;
        } else {
            out->n_verts            = info.n_verts;
            out->n_faces            = info.n_faces;
            out->n_edges            = info.n_edges;
            out->n_boundary_edges   = info.n_boundary_edges;
            out->n_manifold_edges   = info.n_manifold_edges;
            out->n_nonmanifold_edges = info.n_nonmanifold_edges;
            out->n_boundary_loops   = info.n_boundary_loops;
            out->euler              = info.euler;
            out->genus              = info.genus;
            out->n_unref_verts      = info.n_unref_verts;
            out->is_disk            = info.is_disk;
            out->n_components = sf_count_components(g.arena, in->n_vertices,
                                                    faces, in->n_faces);
        }
    }
    SF_EXCEPT_TAIL(g, rc)
    END_TRY;

    sf_op_end(&g);
    return rc;
}

SF_API sf_status sf_developability_energy(const sf_mesh *in, double eps,
                                          const sf_common_opts *opts,
                                          double **out_lambda,
                                          sf_developability_report *rep)
{
    SfOpGuard          g;
    volatile sf_status rc;

    if (validate_in(in) != SF_OK)
        return SF_ERROR_BAD_ARG;
    if (out_lambda) *out_lambda = NULL;
    if (rep) memset(rep, 0, sizeof *rep);
    if (eps <= 0.0) eps = DEV_CUT_EPS;

    rc = sf_op_begin(&g, opts);
    if (rc != SF_OK)
        return rc;

    TRY {
        ComponentMesh cm;
        sf_cm_from_mesh(g.arena, in, &cm);

        double *lam = ARENA_ALLOC(g.arena, in->n_vertices * sizeof *lam);
        if (Develop_vertex_energy(g.arena, cm.verts, cm.nv, cm.faces, cm.nf,
                                  lam) != 0) {
            rc = SF_ERROR;
        } else {
            if (rep) {
                double mean = 0.0;
                rep->nondev_fraction =
                    DevCut_nondev_fraction(g.arena, &cm, eps, &mean);
                rep->mean_lambda = mean;
            }
            if (out_lambda) {
                double *copy = sf_malloc(in->n_vertices * sizeof *copy);
                if (!copy) {
                    rc = SF_ERROR_OOM;
                } else {
                    memcpy(copy, lam, in->n_vertices * sizeof *copy);
                    *out_lambda = copy;
                }
            }
        }
    }
    SF_EXCEPT_TAIL(g, rc)
    END_TRY;

    sf_op_end(&g);
    return rc;
}

SF_API sf_status sf_count_sheets(const sf_mesh *in, int grid_size,
                                 const sf_common_opts *opts, int *out_sheets)
{
    SfOpGuard          g;
    volatile sf_status rc;

    if (!out_sheets || validate_in(in) != SF_OK)
        return SF_ERROR_BAD_ARG;
    *out_sheets = 0;
    if (grid_size <= 0)
        grid_size = ORACLE_UV_GRID_SIZE;

    rc = sf_op_begin(&g, opts);
    if (rc != SF_OK)
        return rc;

    TRY {
        ComponentMesh cm;
        sf_cm_from_mesh(g.arena, in, &cm);
        *out_sheets = Oracle_count_sheets(g.arena, &cm, grid_size);
    }
    SF_EXCEPT_TAIL(g, rc)
    END_TRY;

    sf_op_end(&g);
    return rc;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Cleanup family
 * ════════════════════════════════════════════════════════════════════════ */

SF_API sf_status sf_cleanup(const sf_mesh *in, const sf_cleanup_config *cfg,
                            sf_mesh *out, sf_cleanup_report *rep)
{
    SfOpGuard          g;
    volatile sf_status rc;
    sf_cleanup_config  def;

    if (!out || validate_in(in) != SF_OK)
        return SF_ERROR_BAD_ARG;
    memset(out, 0, sizeof *out);
    if (rep) memset(rep, 0, sizeof *rep);
    if (!cfg) { def = sf_cleanup_config_default(); cfg = &def; }

    rc = sf_op_begin(&g, &cfg->common);
    if (rc != SF_OK)
        return rc;

    TRY {
        ComponentMesh cm;
        sf_cm_from_mesh(g.arena, in, &cm);

        const float *orig_zyx = cm.verts;   /* arena copies persist */
        const size_t orig_nv  = cm.nv;

        if (cfg->manifold_repair) {
            ManifoldGuardStats mg;
            memset(&mg, 0, sizeof mg);
            if (ManifoldGuard_process(g.arena, &cm, 1, cfg->reorient, &mg) != 0)
                rc = SF_ERROR;
            else if (rep) {
                rep->nm_edges_resolved = mg.nm_edges_resolved;
                rep->faces_deleted     = mg.faces_deleted;
                rep->pinch_splits      = mg.pinch_splits;
                rep->orient_flips      = mg.orient_flips;
            }
            RunCtx_progress_check("cleanup", 0.25);
        }

        if (rc == SF_OK && cfg->fill_pinholes) {
            size_t splits = 0, filled = 0, tris = 0, skipped = 0;
            if (PinholeFill_process(g.arena, &cm, 1, cfg->respect_pins,
                                    &splits, &filled, &tris, &skipped) != 0)
                rc = SF_ERROR;
            else if (rep) {
                rep->pinch_splits   += splits;
                rep->pinholes_filled = filled;
                rep->tris_added      = tris;
                rep->loops_skipped   = skipped;
            }
            RunCtx_progress_check("cleanup", 0.5);
        }

        if (rc == SF_OK && cfg->sliver_cleanup) {
            WeldCleanupParams p;
            WeldCleanupStats  ws;
            WeldCleanup_default_params(&p);
            p.sliver_min_alt   = cfg->sliver_min_alt;
            p.degen_area       = cfg->degen_area;
            p.max_collapse_len = cfg->max_collapse_len;
            memset(&ws, 0, sizeof ws);
            if (WeldCleanup_process(g.arena, &cm, &p, &ws) != 0)
                rc = SF_ERROR;
            else if (rep) {
                rep->n_flips     = ws.n_flips;
                rep->n_collapses = ws.n_collapses;
            }
            RunCtx_progress_check("cleanup", 0.75);
        }

        ComponentMesh *final_cm = &cm;
        ComponentMesh  merged;
        if (rc == SF_OK && cfg->cull_min_area_frac > 0.0f) {
            ComponentMesh *culled  = NULL;
            size_t         n_out   = 0;
            size_t         n_culled = 0;
            if (ComponentCull_by_area(g.arena, &cm, 1, cfg->cull_min_area_frac,
                                      &culled, &n_out, &n_culled) != 0) {
                rc = SF_ERROR;
            } else if (n_out == 0) {
                rc = SF_ERROR;   /* everything culled — refuse to emit empty */
            } else {
                if (rep) rep->components_culled = n_culled;
                if (n_out == 1) {
                    final_cm = &culled[0];
                } else {
                    if (sf_concat_cms(g.arena, culled, n_out, &merged) != 0)
                        rc = SF_ERROR;
                    else
                        final_cm = &merged;
                }
            }
        }

        if (rc == SF_OK) {
            const int32_t *vmap = sf_vmap_by_position(
                g.arena, orig_zyx, orig_nv, final_cm->verts, final_cm->nv);
            RunCtx_progress_check("cleanup", 1.0);
            rc = sf_mesh_out_cm(out, final_cm, vmap);
        }
    }
    SF_EXCEPT_TAIL(g, rc)
    END_TRY;

    if (rc != SF_OK)
        sf_mesh_free(out);
    sf_op_end(&g);
    return rc;
}

SF_API sf_status sf_fill_holes(const sf_mesh *in, const sf_holefill_config *cfg,
                               sf_mesh *out, sf_holefill_report *rep)
{
    SfOpGuard          g;
    volatile sf_status rc;
    sf_holefill_config def;

    if (!out || validate_in(in) != SF_OK)
        return SF_ERROR_BAD_ARG;
    memset(out, 0, sizeof *out);
    if (rep) memset(rep, 0, sizeof *rep);
    if (!cfg) { def = sf_holefill_config_default(); cfg = &def; }

    rc = sf_op_begin(&g, &cfg->common);
    if (rc != SF_OK)
        return rc;

    TRY {
        float   *verts = sf_xyz_to_zyx(g.arena, in->vertices, in->n_vertices);
        int32_t *faces = sf_faces_swap(g.arena, in->faces, in->n_faces);
        size_t   nv    = in->n_vertices;
        size_t   nf    = in->n_faces;
        size_t   n_loops = 0, n_interior = 0, n_filled = 0;

        if (HoleFill_process_ex(g.arena, &verts, &faces, &nv, &nf, NULL,
                                cfg->interior_only,
                                &n_loops, &n_interior, &n_filled) != 0) {
            rc = SF_ERROR;
        } else {
            if (rep) {
                rep->n_loops    = n_loops;
                rep->n_interior = n_interior;
                rep->n_filled   = n_filled;
            }
            /* Originals keep their indices; appended fill vertices are new. */
            int32_t *vmap = ARENA_ALLOC(g.arena, nv * sizeof *vmap);
            for (size_t i = 0; i < nv; i++)
                vmap[i] = (i < in->n_vertices) ? (int32_t)i : -1;
            rc = sf_mesh_out(out, verts, nv, faces, nf, NULL, NULL, vmap);
        }
    }
    SF_EXCEPT_TAIL(g, rc)
    END_TRY;

    if (rc != SF_OK)
        sf_mesh_free(out);
    sf_op_end(&g);
    return rc;
}

SF_API sf_status sf_decimate(const sf_mesh *in, const sf_decimate_config *cfg,
                             sf_mesh *out, sf_decimate_report *rep)
{
    SfOpGuard          g;
    volatile sf_status rc;
    sf_decimate_config def;

    if (!out || validate_in(in) != SF_OK)
        return SF_ERROR_BAD_ARG;
    memset(out, 0, sizeof *out);
    if (rep) memset(rep, 0, sizeof *rep);
    if (!cfg) { def = sf_decimate_config_default(); cfg = &def; }

    rc = sf_op_begin(&g, &cfg->common);
    if (rc != SF_OK)
        return rc;

    TRY {
        float   *verts = sf_xyz_to_zyx(g.arena, in->vertices, in->n_vertices);
        int32_t *faces = sf_faces_swap(g.arena, in->faces, in->n_faces);

        size_t target_nf = cfg->target_faces;
        if (target_nf == 0) {
            float ratio = (cfg->target_ratio > 0.0f && cfg->target_ratio <= 1.0f)
                              ? cfg->target_ratio : QEM_TARGET_RATIO;
            target_nf = (size_t)((double)in->n_faces * (double)ratio);
            if (target_nf < 1)
                target_nf = 1;
        }

        const uint8_t *pins = (cfg->respect_pins && in->pin_mask)
                                  ? in->pin_mask : NULL;
        float   *ov = NULL;
        int32_t *of = NULL;
        uint8_t *op = NULL;
        size_t   onv = 0, onf = 0;

        if (QEM_simplify_pinned(g.arena, verts, in->n_vertices,
                                faces, in->n_faces,
                                pins, NULL, target_nf,
                                &ov, &onv, &of, &onf,
                                pins ? &op : NULL) != 0) {
            rc = SF_ERROR;
        } else {
            if (rep) {
                rep->faces_in  = in->n_faces;
                rep->faces_out = onf;
                rep->verts_out = onv;
            }
            /* Decimation moves surviving vertices: no provenance (vmap NULL). */
            rc = sf_mesh_out(out, ov, onv, of, onf, NULL, op, NULL);
        }
    }
    SF_EXCEPT_TAIL(g, rc)
    END_TRY;

    if (rc != SF_OK)
        sf_mesh_free(out);
    sf_op_end(&g);
    return rc;
}

SF_API sf_status sf_fair_developable(const sf_mesh *in, const sf_fair_config *cfg,
                                     sf_mesh *out, sf_fair_report *rep)
{
    SfOpGuard          g;
    volatile sf_status rc;
    sf_fair_config     def;

    if (!out || validate_in(in) != SF_OK)
        return SF_ERROR_BAD_ARG;
    memset(out, 0, sizeof *out);
    if (rep) memset(rep, 0, sizeof *rep);
    if (!cfg) { def = sf_fair_config_default(); cfg = &def; }

    rc = sf_op_begin(&g, &cfg->common);
    if (rc != SF_OK)
        return rc;

    TRY {
        float   *verts   = sf_xyz_to_zyx(g.arena, in->vertices, in->n_vertices);
        int32_t *faces   = sf_faces_swap(g.arena, in->faces, in->n_faces);
        float   *normals = in->vertex_normals
                               ? sf_xyz_to_zyx(g.arena, in->vertex_normals,
                                               in->n_vertices)
                               : NULL;

        DevelopOpts  opts;
        DevelopStats stats;
        memset(&opts, 0, sizeof opts);
        memset(&stats, 0, sizeof stats);
        opts.max_iters    = cfg->max_iters;
        opts.tol_grad     = cfg->tol_grad;
        opts.tol_e        = cfg->tol_e;
        opts.max_disp_vox = cfg->max_disp_vox;
        opts.pin_boundary = cfg->pin_boundary;

        if (Develop_optimize(g.arena, verts, in->n_vertices,
                             faces, in->n_faces,
                             in->pin_mask, normals, &opts, &stats) != 0) {
            rc = SF_ERROR;
        } else {
            if (rep) {
                rep->e_initial = stats.e_initial;
                rep->e_final   = stats.e_final;
                rep->max_disp  = stats.max_disp;
                rep->iters     = stats.iters;
                rep->converged = stats.converged;
            }
            int32_t *vmap = sf_vmap_identity(g.arena, in->n_vertices);
            rc = sf_mesh_out(out, verts, in->n_vertices, faces, in->n_faces,
                             normals, in->pin_mask, vmap);
        }
    }
    SF_EXCEPT_TAIL(g, rc)
    END_TRY;

    if (rc != SF_OK)
        sf_mesh_free(out);
    sf_op_end(&g);
    return rc;
}

SF_API sf_status sf_orient(const sf_mesh *in, const sf_orient_config *cfg,
                           sf_mesh *out, sf_orient_report *rep)
{
    SfOpGuard          g;
    volatile sf_status rc;
    sf_orient_config   def;

    if (!out || validate_in(in) != SF_OK)
        return SF_ERROR_BAD_ARG;
    memset(out, 0, sizeof *out);
    if (rep) memset(rep, 0, sizeof *rep);
    if (!cfg) { def = sf_orient_config_default(); cfg = &def; }

    rc = sf_op_begin(&g, &cfg->common);
    if (rc != SF_OK)
        return rc;

    TRY {
        float   *verts   = sf_xyz_to_zyx(g.arena, in->vertices, in->n_vertices);
        int32_t *faces   = sf_faces_swap(g.arena, in->faces, in->n_faces);
        float   *normals = in->vertex_normals
                               ? sf_xyz_to_zyx(g.arena, in->vertex_normals,
                                               in->n_vertices)
                               : NULL;
        size_t flips = 0, comps = 0, resid = 0, flipped = 0;

        if (OrientMesh_consistent(g.arena, verts, in->n_vertices, normals,
                                  faces, in->n_faces,
                                  &flips, &comps, &resid) != 0) {
            rc = SF_ERROR;
        } else {
            if (cfg->align_components) {
                if (OrientWeld_components(g.arena, verts, in->n_vertices,
                                          faces, in->n_faces,
                                          cfg->align_radius, &flipped) != 0)
                    rc = SF_ERROR;
            }
        }
        if (rc == SF_OK) {
            if (rep) {
                rep->flips              = flips;
                rep->components         = comps;
                rep->residual_same_dir  = resid;
                rep->components_flipped = flipped;
            }
            int32_t *vmap = sf_vmap_identity(g.arena, in->n_vertices);
            rc = sf_mesh_out(out, verts, in->n_vertices, faces, in->n_faces,
                             normals, in->pin_mask, vmap);
        }
    }
    SF_EXCEPT_TAIL(g, rc)
    END_TRY;

    if (rc != SF_OK)
        sf_mesh_free(out);
    sf_op_end(&g);
    return rc;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Detangle / split family
 * ════════════════════════════════════════════════════════════════════════ */

/* Copy-out a piece array as an sf_mesh_list with per-piece provenance vmaps.
 * vmaps are computed first (arena; may raise); the malloc emit never raises. */
static sf_status emit_pieces(Arena_T a,
                             const float *in_zyx, size_t in_nv,
                             const ComponentMesh *pieces, size_t n,
                             sf_mesh_list *out)
{
    int32_t **vmaps = ARENA_ALLOC(a, n * sizeof *vmaps);
    for (size_t i = 0; i < n; i++)
        vmaps[i] = sf_vmap_by_position(a, in_zyx, in_nv,
                                       pieces[i].verts, pieces[i].nv);

    sf_mesh_list lst;
    memset(&lst, 0, sizeof lst);
    lst.items = sf_malloc(n * sizeof *lst.items);
    if (!lst.items)
        return SF_ERROR_OOM;
    memset(lst.items, 0, n * sizeof *lst.items);
    lst.count = n;

    for (size_t i = 0; i < n; i++) {
        sf_status rc = sf_mesh_out_cm(&lst.items[i], &pieces[i], vmaps[i]);
        if (rc != SF_OK) {
            sf_mesh_list_free(&lst);
            return rc;
        }
    }
    *out = lst;
    return SF_OK;
}

SF_API sf_status sf_split_depth_peel(const sf_mesh *in,
                                     const sf_split_peel_config *cfg,
                                     sf_mesh_list *out, sf_split_report *rep)
{
    SfOpGuard            g;
    volatile sf_status   rc;
    sf_split_peel_config def;

    if (!out || validate_in(in) != SF_OK)
        return SF_ERROR_BAD_ARG;
    memset(out, 0, sizeof *out);
    if (rep) memset(rep, 0, sizeof *rep);
    if (!cfg) { def = sf_split_peel_config_default(); cfg = &def; }

    rc = sf_op_begin(&g, &cfg->common);
    if (rc != SF_OK)
        return rc;

    TRY {
        ComponentMesh cm;
        sf_cm_from_mesh(g.arena, in, &cm);

        ComponentMesh *pieces = NULL;
        size_t         n      = 0;
        if (DepthPeel_process(g.arena, &cm, cfg->min_gap, cfg->gap_depth,
                              cfg->min_comp_verts, &pieces, &n) != 0 || n == 0) {
            rc = SF_ERROR;
        } else {
            if (rep) {
                rep->n_pieces   = n;
                rep->self_gated = (n == 1);
            }
            rc = emit_pieces(g.arena, cm.verts, cm.nv, pieces, n, out);
        }
    }
    SF_EXCEPT_TAIL(g, rc)
    END_TRY;

    if (rc != SF_OK)
        sf_mesh_list_free(out);
    sf_op_end(&g);
    return rc;
}

SF_API sf_status sf_split_developability(const sf_mesh *in,
                                         const sf_split_dev_config *cfg,
                                         sf_mesh_list *out, sf_split_report *rep)
{
    SfOpGuard           g;
    volatile sf_status  rc;
    sf_split_dev_config def;

    if (!out || validate_in(in) != SF_OK)
        return SF_ERROR_BAD_ARG;
    memset(out, 0, sizeof *out);
    if (rep) memset(rep, 0, sizeof *rep);
    if (!cfg) { def = sf_split_dev_config_default(); cfg = &def; }

    rc = sf_op_begin(&g, &cfg->common);
    if (rc != SF_OK)
        return rc;

    TRY {
        ComponentMesh cm;
        sf_cm_from_mesh(g.arena, in, &cm);

        double parent_frac = 0.0;
        if (cfg->apply_gate || rep)
            parent_frac = DevCut_nondev_fraction(g.arena, &cm, cfg->eps, NULL);

        ComponentMesh *pieces = NULL;
        size_t         n      = 0;
        if (DevCut_process(g.arena, &cm, cfg->eps, cfg->gap_depth,
                           cfg->min_comp_verts, &pieces, &n) != 0 || n == 0) {
            rc = SF_ERROR;
        } else {
            int    accept        = 1;
            double worst_piece   = 0.0;
            if (n > 1 && cfg->apply_gate) {
                for (size_t i = 0; i < n; i++) {
                    double f = DevCut_nondev_fraction(g.arena, &pieces[i],
                                                      cfg->eps, NULL);
                    if (f > worst_piece)
                        worst_piece = f;
                    if (f > parent_frac - cfg->gate_margin ||
                        f > cfg->gate_abs_ceil)
                        accept = 0;
                }
            }
            if (n > 1 && !accept) {
                /* Gate rejected the split: pass the input through unchanged. */
                pieces = &cm;
                n      = 1;
            }
            if (rep) {
                rep->n_pieces         = n;
                rep->self_gated       = (n == 1);
                rep->nondev_before    = parent_frac;
                rep->nondev_after_max = worst_piece;
            }
            rc = emit_pieces(g.arena, cm.verts, cm.nv, pieces, n, out);
        }
    }
    SF_EXCEPT_TAIL(g, rc)
    END_TRY;

    if (rc != SF_OK)
        sf_mesh_list_free(out);
    sf_op_end(&g);
    return rc;
}

SF_API sf_status sf_split_bridges(const sf_mesh *in,
                                  const sf_split_bridge_config *cfg,
                                  sf_mesh_list *out, sf_split_report *rep)
{
    SfOpGuard              g;
    volatile sf_status     rc;
    sf_split_bridge_config def;

    if (!out || validate_in(in) != SF_OK)
        return SF_ERROR_BAD_ARG;
    memset(out, 0, sizeof *out);
    if (rep) memset(rep, 0, sizeof *rep);
    if (!cfg) { def = sf_split_bridge_config_default(); cfg = &def; }

    rc = sf_op_begin(&g, &cfg->common);
    if (rc != SF_OK)
        return rc;

    TRY {
        ComponentMesh cm;
        sf_cm_from_mesh(g.arena, in, &cm);

        ComponentMesh *pieces = NULL;
        size_t         n      = 0;
        if (BridgeCut_process(g.arena, &cm, cfg->oracle_grid_size,
                              cfg->component_timeout_sec, &pieces, &n) != 0 ||
            n == 0) {
            rc = SF_ERROR;
        } else {
            if (rep) {
                rep->n_pieces   = n;
                rep->self_gated = (n == 1);
            }
            rc = emit_pieces(g.arena, cm.verts, cm.nv, pieces, n, out);
        }
    }
    SF_EXCEPT_TAIL(g, rc)
    END_TRY;

    if (rc != SF_OK)
        sf_mesh_list_free(out);
    sf_op_end(&g);
    return rc;
}

SF_API sf_status sf_split_overlap(const sf_mesh *in,
                                  const sf_split_overlap_config *cfg,
                                  sf_mesh_list *out, sf_split_report *rep)
{
    SfOpGuard               g;
    volatile sf_status      rc;
    sf_split_overlap_config def;

    if (!out || validate_in(in) != SF_OK)
        return SF_ERROR_BAD_ARG;
    memset(out, 0, sizeof *out);
    if (rep) memset(rep, 0, sizeof *rep);
    if (!cfg) { def = sf_split_overlap_config_default(); cfg = &def; }

    rc = sf_op_begin(&g, &cfg->common);
    if (rc != SF_OK)
        return rc;

    TRY {
        ComponentMesh cm;
        sf_cm_from_mesh(g.arena, in, &cm);

        ComponentMesh *pieces = NULL;
        size_t         n      = 0;
        int mrc = OverlapSep_process(g.arena, &cm, cfg->sheet_count,
                                     cfg->common.n_threads,
                                     cfg->component_timeout_sec, &pieces, &n);
        /* -1 means "no separation needed/possible" with a pass-through copy
         * in pieces[0]; that is a successful single-piece result here. */
        if ((mrc != 0 && mrc != -1) || n == 0 || !pieces) {
            rc = SF_ERROR;
        } else {
            if (rep) {
                rep->n_pieces   = n;
                rep->self_gated = (n == 1);
            }
            rc = emit_pieces(g.arena, cm.verts, cm.nv, pieces, n, out);
        }
    }
    SF_EXCEPT_TAIL(g, rc)
    END_TRY;

    if (rc != SF_OK)
        sf_mesh_list_free(out);
    sf_op_end(&g);
    return rc;
}

SF_API sf_status sf_detangle(const sf_mesh *in, const sf_detangle_config *cfg,
                             sf_mesh_list *out, sf_detangle_report *rep)
{
    SfOpGuard           g;
    volatile sf_status  rc;
    sf_detangle_config  def;

    if (!out || validate_in(in) != SF_OK)
        return SF_ERROR_BAD_ARG;
    memset(out, 0, sizeof *out);
    if (rep) memset(rep, 0, sizeof *rep);
    if (!cfg) { def = sf_detangle_config_default(); cfg = &def; }

    rc = sf_op_begin(&g, &cfg->common);
    if (rc != SF_OK)
        return rc;

    TRY {
        ComponentMesh cm;
        sf_cm_from_mesh(g.arena, in, &cm);

        /* Current piece set, replaced stage by stage. */
        ComponentMesh *cur   = &cm;
        size_t         n_cur = 1;

        struct {
            const char *name;
            int         enabled;
            int         which;   /* 0=peel 1=dev 2=bridge 3=overlap */
        } stages[4] = {
            { "peel",    cfg->enable_peel,    0 },
            { "devcut",  cfg->enable_dev,     1 },
            { "bridge",  cfg->enable_bridge,  2 },
            { "overlap", cfg->enable_overlap, 3 },
        };
        size_t stage_splits[4] = { 0, 0, 0, 0 };

        for (int s = 0; s < 4 && rc == SF_OK; s++) {
            if (!stages[s].enabled)
                continue;
            RunCtx_progress_check("detangle", (double)s / 4.0);

            ComponentMesh **piece_out  = ARENA_ALLOC(g.arena,
                                            n_cur * sizeof *piece_out);
            size_t         *piece_n    = ARENA_ALLOC(g.arena,
                                            n_cur * sizeof *piece_n);
            size_t          total      = 0;

            for (size_t p = 0; p < n_cur && rc == SF_OK; p++) {
                ComponentMesh *pieces = NULL;
                size_t         n      = 0;
                int            mrc    = 0;

                switch (stages[s].which) {
                case 0:
                    mrc = DepthPeel_process(g.arena, &cur[p],
                                            DEPTH_PEEL_MIN_GAP_VOX,
                                            DEPTH_PEEL_GAP_DEPTH,
                                            cfg->min_comp_verts, &pieces, &n);
                    break;
                case 1: {
                    double parent = DevCut_nondev_fraction(g.arena, &cur[p],
                                                           DEV_CUT_EPS, NULL);
                    mrc = DevCut_process(g.arena, &cur[p], DEV_CUT_EPS,
                                         DEV_CUT_GAP_DEPTH,
                                         cfg->min_comp_verts, &pieces, &n);
                    if (mrc == 0 && n > 1) {
                        int accept = 1;
                        for (size_t i = 0; i < n; i++) {
                            double f = DevCut_nondev_fraction(g.arena,
                                                              &pieces[i],
                                                              DEV_CUT_EPS,
                                                              NULL);
                            if (f > parent - DEV_CUT_GATE_MARGIN ||
                                f > DEV_CUT_GATE_ABS_CEIL)
                                accept = 0;
                        }
                        if (!accept) {
                            pieces = &cur[p];
                            n      = 1;
                        }
                    }
                    break;
                }
                case 2:
                    mrc = BridgeCut_process(g.arena, &cur[p],
                                            ORACLE_UV_GRID_SIZE,
                                            cfg->component_timeout_sec,
                                            &pieces, &n);
                    break;
                case 3:
                    mrc = OverlapSep_process(g.arena, &cur[p], 0,
                                             cfg->common.n_threads,
                                             cfg->component_timeout_sec,
                                             &pieces, &n);
                    if (mrc == -1 && n >= 1 && pieces)
                        mrc = 0;   /* pass-through copy */
                    break;
                }

                if (mrc != 0 || n == 0 || !pieces) {
                    /* A stage failing on one piece keeps that piece as-is. */
                    piece_out[p] = &cur[p];
                    piece_n[p]   = 1;
                    total       += 1;
                } else {
                    piece_out[p] = pieces;
                    piece_n[p]   = n;
                    total       += n;
                    if (n > 1)
                        stage_splits[s] += n - 1;
                }
            }

            if (rc == SF_OK) {
                ComponentMesh *next = ARENA_ALLOC(g.arena,
                                                  total * sizeof *next);
                size_t k = 0;
                for (size_t p = 0; p < n_cur; p++) {
                    for (size_t i = 0; i < piece_n[p]; i++) {
                        next[k] = piece_out[p][i];
                        next[k].self = &next[k];
                        k++;
                    }
                }
                cur   = next;
                n_cur = total;
            }
        }

        if (rc == SF_OK) {
            if (rep) {
                rep->n_pieces       = n_cur;
                rep->peel_splits    = stage_splits[0];
                rep->dev_splits     = stage_splits[1];
                rep->bridge_splits  = stage_splits[2];
                rep->overlap_splits = stage_splits[3];
            }
            RunCtx_progress_check("detangle", 1.0);
            rc = emit_pieces(g.arena, cm.verts, cm.nv, cur, n_cur, out);
        }
    }
    SF_EXCEPT_TAIL(g, rc)
    END_TRY;

    if (rc != SF_OK)
        sf_mesh_list_free(out);
    sf_op_end(&g);
    return rc;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Topology surgery
 * ════════════════════════════════════════════════════════════════════════ */

/* Split a raw (verts,faces) mesh into face-connectivity components with
 * local->input vertex maps. Arena-backed; may raise. */
typedef struct {
    float   *verts;
    size_t   nv;
    int32_t *faces;
    size_t   nf;
    int32_t *orig_of_local;   /* [nv] local vertex -> input vertex index */
} SfRawComp;

static int32_t rawcomp_find(int32_t *parent, int32_t x)
{
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

static size_t split_raw_components(Arena_T a, const float *verts, size_t nv,
                                   const int32_t *faces, size_t nf,
                                   SfRawComp **out_comps)
{
    int32_t *parent = ARENA_ALLOC(a, nv * sizeof *parent);
    for (size_t i = 0; i < nv; i++)
        parent[i] = (int32_t)i;
    for (size_t f = 0; f < nf; f++) {
        int32_t r0 = rawcomp_find(parent, faces[f * 3 + 0]);
        int32_t r1 = rawcomp_find(parent, faces[f * 3 + 1]);
        if (r1 != r0) parent[r1] = r0;
        r0 = rawcomp_find(parent, faces[f * 3 + 0]);
        int32_t r2 = rawcomp_find(parent, faces[f * 3 + 2]);
        if (r2 != r0) parent[r2] = r0;
    }

    /* Assign dense component ids to roots of referenced vertices. */
    int32_t *comp_of_root = ARENA_ALLOC(a, nv * sizeof *comp_of_root);
    for (size_t i = 0; i < nv; i++)
        comp_of_root[i] = -1;
    size_t n_comps = 0;
    for (size_t f = 0; f < nf; f++) {
        int32_t r = rawcomp_find(parent, faces[f * 3]);
        if (comp_of_root[r] < 0)
            comp_of_root[r] = (int32_t)n_comps++;
    }
    if (n_comps == 0) {
        *out_comps = NULL;
        return 0;
    }

    SfRawComp *comps = ARENA_CALLOC(a, n_comps, sizeof *comps);

    /* Count faces + verts per component. */
    size_t *fcount = ARENA_CALLOC(a, n_comps, sizeof *fcount);
    for (size_t f = 0; f < nf; f++)
        fcount[comp_of_root[rawcomp_find(parent, faces[f * 3])]]++;

    int32_t *local_of_orig = ARENA_ALLOC(a, nv * sizeof *local_of_orig);
    for (size_t i = 0; i < nv; i++)
        local_of_orig[i] = -1;

    size_t *vcount = ARENA_CALLOC(a, n_comps, sizeof *vcount);
    for (size_t f = 0; f < nf; f++) {
        size_t c = (size_t)comp_of_root[rawcomp_find(parent, faces[f * 3])];
        for (int k = 0; k < 3; k++) {
            int32_t v = faces[f * 3 + k];
            if (local_of_orig[v] < 0) {
                local_of_orig[v] = (int32_t)vcount[c]++;
            }
        }
    }
    /* local_of_orig is per-component; reset and redo per comp allocation. */
    for (size_t i = 0; i < nv; i++)
        local_of_orig[i] = -1;

    for (size_t c = 0; c < n_comps; c++) {
        comps[c].verts = ARENA_ALLOC(a, vcount[c] * 3 * sizeof(float));
        comps[c].faces = ARENA_ALLOC(a, fcount[c] * 3 * sizeof(int32_t));
        comps[c].orig_of_local = ARENA_ALLOC(a, vcount[c] * sizeof(int32_t));
        comps[c].nv = 0;
        comps[c].nf = 0;
    }
    for (size_t f = 0; f < nf; f++) {
        size_t c = (size_t)comp_of_root[rawcomp_find(parent, faces[f * 3])];
        SfRawComp *rcmp = &comps[c];
        for (int k = 0; k < 3; k++) {
            int32_t v = faces[f * 3 + k];
            if (local_of_orig[v] < 0) {
                int32_t lv = (int32_t)rcmp->nv++;
                local_of_orig[v] = lv;
                memcpy(&rcmp->verts[lv * 3], &verts[(size_t)v * 3],
                       3 * sizeof(float));
                rcmp->orig_of_local[lv] = v;
            }
            rcmp->faces[rcmp->nf * 3 + k] = local_of_orig[v];
        }
        rcmp->nf++;
    }

    *out_comps = comps;
    return n_comps;
}

SF_API sf_status sf_cut_to_disk(const sf_mesh *in, const sf_seamcut_config *cfg,
                                sf_mesh *out, sf_seamcut_report *rep)
{
    SfOpGuard          g;
    volatile sf_status rc;
    sf_seamcut_config  def;

    if (!out || validate_in(in) != SF_OK)
        return SF_ERROR_BAD_ARG;
    memset(out, 0, sizeof *out);
    if (rep) memset(rep, 0, sizeof *rep);
    if (!cfg) { def = sf_seamcut_config_default(); cfg = &def; }

    rc = sf_op_begin(&g, &cfg->common);
    if (rc != SF_OK)
        return rc;

    TRY {
        float   *verts0 = sf_xyz_to_zyx(g.arena, in->vertices, in->n_vertices);
        int32_t *faces0 = sf_faces_swap(g.arena, in->faces, in->n_faces);

        /* 1. Manifold repair (bowtie/NM-edge duplication). */
        float   *rv = NULL; size_t rnv = 0;
        int32_t *rf = NULL; size_t rnf = 0;
        int32_t *vmap_r = NULL;
        if (SeamCut_repair_manifold(g.arena, verts0, in->n_vertices,
                                    faces0, in->n_faces,
                                    &rv, &rnv, &rf, &rnf, &vmap_r) != 0) {
            rc = SF_ERROR;
        } else {
            RunCtx_progress_check("seamcut", 0.25);

            /* 2. Per-component Seamster cut. */
            SfRawComp *comps  = NULL;
            size_t     ncomps = split_raw_components(g.arena, rv, rnv, rf, rnf,
                                                     &comps);
            if (ncomps == 0) {
                rc = SF_ERROR;
            } else {
                SeamCutOpts opts;
                opts.accept_frac = cfg->accept_frac;
                opts.r_max       = cfg->r_max;

                /* Cut each component; on failure keep it uncut. */
                float   **cv   = ARENA_ALLOC(g.arena, ncomps * sizeof *cv);
                int32_t **cf   = ARENA_ALLOC(g.arena, ncomps * sizeof *cf);
                int32_t **cvm  = ARENA_ALLOC(g.arena, ncomps * sizeof *cvm);
                size_t   *cnv  = ARENA_ALLOC(g.arena, ncomps * sizeof *cnv);
                size_t   *cnf  = ARENA_ALLOC(g.arena, ncomps * sizeof *cnf);
                int       all_disks = 1;

                for (size_t c = 0; c < ncomps; c++) {
                    SeamCutResult res;
                    memset(&res, 0, sizeof res);
                    if (SeamCut_run(g.arena, comps[c].verts, comps[c].nv,
                                    comps[c].faces, comps[c].nf,
                                    &opts, &res) == 0) {
                        cv[c]  = res.verts;
                        cf[c]  = res.faces;
                        cnv[c] = res.nv;
                        cnf[c] = res.nf;
                        /* compose: cut-local -> comp-local -> repaired -> input */
                        int32_t *m = ARENA_ALLOC(g.arena,
                                                 res.nv * sizeof *m);
                        for (size_t i = 0; i < res.nv; i++) {
                            int32_t comp_local = res.vmap ? res.vmap[i]
                                                          : (int32_t)i;
                            int32_t repaired =
                                comps[c].orig_of_local[comp_local];
                            m[i] = vmap_r ? vmap_r[repaired] : repaired;
                        }
                        cvm[c] = m;
                        if (rep) {
                            rep->genus_before += res.genus_before;
                            rep->genus_after  += res.genus_after;
                            rep->loops_before += res.loops_before;
                            rep->loops_after  += res.loops_after;
                            rep->n_terminals  += res.n_terminals;
                            rep->seam_edges   += res.seam_edges;
                        }
                        if (!res.is_disk_after)
                            all_disks = 0;
                    } else {
                        /* Keep the component uncut. */
                        cv[c]  = comps[c].verts;
                        cf[c]  = comps[c].faces;
                        cnv[c] = comps[c].nv;
                        cnf[c] = comps[c].nf;
                        int32_t *m = ARENA_ALLOC(g.arena,
                                                 comps[c].nv * sizeof *m);
                        for (size_t i = 0; i < comps[c].nv; i++) {
                            int32_t repaired = comps[c].orig_of_local[i];
                            m[i] = vmap_r ? vmap_r[repaired] : repaired;
                        }
                        cvm[c] = m;
                        all_disks = 0;
                    }
                    RunCtx_progress_check("seamcut",
                                          0.25 + 0.7 * (double)(c + 1) /
                                                     (double)ncomps);
                }

                /* 3. Concatenate. */
                size_t tot_nv = 0, tot_nf = 0;
                for (size_t c = 0; c < ncomps; c++) {
                    tot_nv += cnv[c];
                    tot_nf += cnf[c];
                }
                float   *fv = ARENA_ALLOC(g.arena, tot_nv * 3 * sizeof *fv);
                int32_t *ff = ARENA_ALLOC(g.arena, tot_nf * 3 * sizeof *ff);
                int32_t *fm = ARENA_ALLOC(g.arena, tot_nv * sizeof *fm);
                size_t voff = 0, foff = 0;
                for (size_t c = 0; c < ncomps; c++) {
                    memcpy(fv + voff * 3, cv[c], cnv[c] * 3 * sizeof(float));
                    memcpy(fm + voff, cvm[c], cnv[c] * sizeof(int32_t));
                    for (size_t i = 0; i < cnf[c] * 3; i++)
                        ff[foff * 3 + i] = cf[c][i] + (int32_t)voff;
                    voff += cnv[c];
                    foff += cnf[c];
                }
                if (rep) {
                    rep->n_components  = (long)ncomps;
                    rep->is_disk_after = all_disks;
                }
                rc = sf_mesh_out(out, fv, tot_nv, ff, tot_nf, NULL, NULL, fm);
            }
        }
    }
    SF_EXCEPT_TAIL(g, rc)
    END_TRY;

    if (rc != SF_OK)
        sf_mesh_free(out);
    sf_op_end(&g);
    return rc;
}

SF_API sf_status sf_sever_short_handles(const sf_mesh *in, double max_loop_len,
                                        const sf_common_opts *opts,
                                        sf_mesh *out, long *out_loops_cut)
{
    SfOpGuard          g;
    volatile sf_status rc;

    if (!out || validate_in(in) != SF_OK)
        return SF_ERROR_BAD_ARG;
    memset(out, 0, sizeof *out);
    if (out_loops_cut) *out_loops_cut = 0;

    rc = sf_op_begin(&g, opts);
    if (rc != SF_OK)
        return rc;

    TRY {
        float   *verts0 = sf_xyz_to_zyx(g.arena, in->vertices, in->n_vertices);
        int32_t *faces0 = sf_faces_swap(g.arena, in->faces, in->n_faces);

        SfRawComp *comps  = NULL;
        size_t     ncomps = split_raw_components(g.arena, verts0,
                                                 in->n_vertices,
                                                 faces0, in->n_faces, &comps);
        if (ncomps == 0) {
            rc = SF_ERROR;
        } else {
            float   **cv  = ARENA_ALLOC(g.arena, ncomps * sizeof *cv);
            int32_t **cf  = ARENA_ALLOC(g.arena, ncomps * sizeof *cf);
            size_t   *cnv = ARENA_ALLOC(g.arena, ncomps * sizeof *cnv);
            size_t   *cnf = ARENA_ALLOC(g.arena, ncomps * sizeof *cnf);
            long      total_cut = 0;

            for (size_t c = 0; c < ncomps; c++) {
                long cut = 0;
                if (SeamCut_sever_short_handles(g.arena,
                                                comps[c].verts, comps[c].nv,
                                                comps[c].faces, comps[c].nf,
                                                max_loop_len,
                                                &cv[c], &cnv[c],
                                                &cf[c], &cnf[c], &cut) != 0) {
                    cv[c]  = comps[c].verts;
                    cf[c]  = comps[c].faces;
                    cnv[c] = comps[c].nv;
                    cnf[c] = comps[c].nf;
                    cut    = 0;
                }
                total_cut += cut;
            }

            size_t tot_nv = 0, tot_nf = 0;
            for (size_t c = 0; c < ncomps; c++) {
                tot_nv += cnv[c];
                tot_nf += cnf[c];
            }
            float   *fv = ARENA_ALLOC(g.arena, tot_nv * 3 * sizeof *fv);
            int32_t *ff = ARENA_ALLOC(g.arena, tot_nf * 3 * sizeof *ff);
            size_t voff = 0, foff = 0;
            for (size_t c = 0; c < ncomps; c++) {
                memcpy(fv + voff * 3, cv[c], cnv[c] * 3 * sizeof(float));
                for (size_t i = 0; i < cnf[c] * 3; i++)
                    ff[foff * 3 + i] = cf[c][i] + (int32_t)voff;
                voff += cnv[c];
                foff += cnf[c];
            }

            /* Severed vertices are duplicated bit-verbatim: recover
             * provenance by exact position. */
            const int32_t *vmap = sf_vmap_by_position(g.arena, verts0,
                                                      in->n_vertices,
                                                      fv, tot_nv);
            if (out_loops_cut)
                *out_loops_cut = total_cut;
            rc = sf_mesh_out(out, fv, tot_nv, ff, tot_nf, NULL, NULL, vmap);
        }
    }
    SF_EXCEPT_TAIL(g, rc)
    END_TRY;

    if (rc != SF_OK)
        sf_mesh_free(out);
    sf_op_end(&g);
    return rc;
}

SF_API sf_status sf_find_handle_loops(const sf_mesh *in,
                                      const sf_common_opts *opts,
                                      sf_handle_loop **out_loops, size_t *out_n)
{
    SfOpGuard          g;
    volatile sf_status rc;

    if (!out_loops || !out_n || validate_in(in) != SF_OK)
        return SF_ERROR_BAD_ARG;
    *out_loops = NULL;
    *out_n     = 0;

    rc = sf_op_begin(&g, opts);
    if (rc != SF_OK)
        return rc;

    TRY {
        float   *verts0 = sf_xyz_to_zyx(g.arena, in->vertices, in->n_vertices);
        int32_t *faces0 = sf_faces_swap(g.arena, in->faces, in->n_faces);

        SfRawComp *comps  = NULL;
        size_t     ncomps = split_raw_components(g.arena, verts0,
                                                 in->n_vertices,
                                                 faces0, in->n_faces, &comps);

        /* Gather loops across components (input-index space). */
        HandleLoop *all   = NULL;
        size_t      n_all = 0;

        for (size_t c = 0; c < ncomps && rc == SF_OK; c++) {
            HandleLoop *loops = NULL;
            size_t      n     = 0;
            if (SeamCut_handle_loops(g.arena, comps[c].verts, comps[c].nv,
                                     comps[c].faces, comps[c].nf,
                                     &loops, &n) != 0)
                continue;   /* skip a failing component */
            if (n == 0)
                continue;
            HandleLoop *grown = ARENA_ALLOC(g.arena,
                                            (n_all + n) * sizeof *grown);
            if (n_all)
                memcpy(grown, all, n_all * sizeof *grown);
            for (size_t i = 0; i < n; i++) {
                HandleLoop hl = loops[i];
                int32_t *remapped = ARENA_ALLOC(g.arena,
                                                hl.n * sizeof *remapped);
                for (size_t k = 0; k < hl.n; k++)
                    remapped[k] = comps[c].orig_of_local[hl.verts[k]];
                hl.verts = remapped;
                grown[n_all + i] = hl;
            }
            all   = grown;
            n_all += n;
        }

        if (rc == SF_OK && n_all > 0) {
            sf_handle_loop *result = sf_malloc(n_all * sizeof *result);
            if (!result) {
                rc = SF_ERROR_OOM;
            } else {
                memset(result, 0, n_all * sizeof *result);
                size_t done = 0;
                for (; done < n_all; done++) {
                    result[done].n      = all[done].n;
                    result[done].length = all[done].length;
                    result[done].vertex_indices =
                        sf_malloc(all[done].n * sizeof(int32_t));
                    if (!result[done].vertex_indices)
                        break;
                    memcpy(result[done].vertex_indices, all[done].verts,
                           all[done].n * sizeof(int32_t));
                }
                if (done < n_all) {
                    sf_handle_loops_free(result, done);
                    rc = SF_ERROR_OOM;
                } else {
                    *out_loops = result;
                    *out_n     = n_all;
                }
            }
        }
    }
    SF_EXCEPT_TAIL(g, rc)
    END_TRY;

    sf_op_end(&g);
    return rc;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Reconstruction
 * ════════════════════════════════════════════════════════════════════════ */

SF_API sf_status sf_reconstruct_bpa(const float *points_xyz,
                                    const float *normals_xyz,
                                    size_t n_points, const sf_bpa_config *cfg,
                                    sf_mesh *out, sf_bpa_report *rep)
{
    SfOpGuard          g;
    volatile sf_status rc;
    sf_bpa_config      def;

    if (!points_xyz || !normals_xyz || n_points < 3 || !out)
        return SF_ERROR_BAD_ARG;
    if (!points_finite(points_xyz, n_points) ||
        !points_finite(normals_xyz, n_points))
        return SF_ERROR_BAD_ARG;
    memset(out, 0, sizeof *out);
    if (rep) memset(rep, 0, sizeof *rep);
    if (!cfg) { def = sf_bpa_config_default(); cfg = &def; }

    rc = sf_op_begin(&g, &cfg->common);
    if (rc != SF_OK)
        return rc;

    TRY {
        float *verts   = sf_xyz_to_zyx(g.arena, points_xyz, n_points);
        float *normals = sf_xyz_to_zyx(g.arena, normals_xyz, n_points);
        float  rho     = cfg->rho > 0.0f ? cfg->rho : BPA_RHO_VOX;

        int32_t *faces = NULL;
        size_t   nf    = 0;
        if (BallPivot_reconstruct(g.arena, verts, normals, n_points, rho,
                                  &faces, &nf) != 0 || nf == 0) {
            rc = SF_ERROR;
        } else {
            if (rep)
                rep->n_faces = nf;
            int32_t *vmap = sf_vmap_identity(g.arena, n_points);
            rc = sf_mesh_out(out, verts, n_points, faces, nf,
                             normals, NULL, vmap);
        }
    }
    SF_EXCEPT_TAIL(g, rc)
    END_TRY;

    if (rc != SF_OK)
        sf_mesh_free(out);
    sf_op_end(&g);
    return rc;
}

SF_API sf_status sf_mls_project(const float *points_xyz, size_t n_points,
                                float radius_vox, const sf_common_opts *opts,
                                float **out_points_xyz, float **out_normals_xyz)
{
    SfOpGuard          g;
    volatile sf_status rc;

    if (!points_xyz || n_points == 0 || !out_points_xyz)
        return SF_ERROR_BAD_ARG;
    if (!points_finite(points_xyz, n_points))
        return SF_ERROR_BAD_ARG;
    *out_points_xyz = NULL;
    if (out_normals_xyz) *out_normals_xyz = NULL;
    if (radius_vox <= 0.0f)
        radius_vox = MLS_PROJECT_RADIUS_VOX;

    rc = sf_op_begin(&g, opts);
    if (rc != SF_OK)
        return rc;

    TRY {
        float *verts = sf_xyz_to_zyx(g.arena, points_xyz, n_points);
        float *pv    = ARENA_ALLOC(g.arena, n_points * 3 * sizeof *pv);
        float *pn    = out_normals_xyz
                           ? ARENA_ALLOC(g.arena, n_points * 3 * sizeof *pn)
                           : NULL;
        const float origin[3] = { 0.0f, 0.0f, 0.0f };

        MLS_project_verts(g.arena, verts, n_points, radius_vox, origin, pv, pn);
        /* The OMP loop inside flag-skips on cancel; raise here (serial). */
        RunCtx_check();

        float *op = sf_malloc(n_points * 3 * sizeof *op);
        float *on = out_normals_xyz ? sf_malloc(n_points * 3 * sizeof *on)
                                    : NULL;
        if (!op || (out_normals_xyz && !on)) {
            free(op);
            free(on);
            rc = SF_ERROR_OOM;
        } else {
            for (size_t i = 0; i < n_points; i++) {
                op[i * 3 + 0] = pv[i * 3 + 2];
                op[i * 3 + 1] = pv[i * 3 + 1];
                op[i * 3 + 2] = pv[i * 3 + 0];
                if (on) {
                    on[i * 3 + 0] = pn[i * 3 + 2];
                    on[i * 3 + 1] = pn[i * 3 + 1];
                    on[i * 3 + 2] = pn[i * 3 + 0];
                }
            }
            *out_points_xyz = op;
            if (out_normals_xyz)
                *out_normals_xyz = on;
        }
    }
    SF_EXCEPT_TAIL(g, rc)
    END_TRY;

    sf_op_end(&g);
    return rc;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Self-test
 * ════════════════════════════════════════════════════════════════════════ */

SF_API int sf_selftest(void)
{
    SfOpGuard    g;
    volatile int fails = 0;

    if (sf_op_begin(&g, NULL) != SF_OK)
        return -1;

    TRY {
        fails += MeshTopo_selftest();
        fails += Develop_selftest();
        fails += DepthPeel_selftest();
        fails += DevCut_selftest();
        fails += SeamCut_selftest();
        fails += ManifoldGuard_selftest(g.arena);
    }
    EXCEPT(Arena_Failed) fails = -1;
    EXCEPT(IO_Failed)    fails = -1;
    EXCEPT(Timeout)      fails = -1;
    EXCEPT(Sf_Cancelled) fails = -1;
    END_TRY;

    sf_op_end(&g);
    return fails;
}

/* ══════════════════════════════════════════════════════════════════════════
 * The API table
 * ════════════════════════════════════════════════════════════════════════ */

static const sf_api g_sf_api = {
    /* abi_version */ SCROLLFIESTA_ABI_VERSION,
    /* struct_size */ (uint32_t)sizeof(sf_api),

    sf_version,
    sf_version_string,
    sf_status_str,

    sf_malloc,
    sf_free,
    sf_mesh_free,
    sf_mesh_list_free,
    sf_handle_loops_free,
    sf_mesh_validate,

    sf_mesh_load_obj,
    sf_mesh_save_obj,

    sf_common_opts_default,
    sf_cleanup_config_default,
    sf_holefill_config_default,
    sf_decimate_config_default,
    sf_fair_config_default,
    sf_orient_config_default,
    sf_split_peel_config_default,
    sf_split_dev_config_default,
    sf_split_bridge_config_default,
    sf_split_overlap_config_default,
    sf_detangle_config_default,
    sf_seamcut_config_default,
    sf_bpa_config_default,
    sf_pipeline_config_default,
    sf_weld_config_default,

    sf_topology_audit,
    sf_developability_energy,
    sf_count_sheets,

    sf_cleanup,
    sf_fill_holes,
    sf_decimate,
    sf_fair_developable,
    sf_orient,

    sf_split_depth_peel,
    sf_split_developability,
    sf_split_bridges,
    sf_split_overlap,
    sf_detangle,

    sf_cut_to_disk,
    sf_sever_short_handles,
    sf_find_handle_loops,

    sf_reconstruct_bpa,
    sf_mls_project,

    sf_pipeline_run,

    /* weld: experimental, unimplemented in 0.9 — NULL by contract. */
    NULL,

    sf_selftest,
};

SF_API const sf_api *sf_get_api(uint32_t requested_abi)
{
    if (requested_abi != SCROLLFIESTA_ABI_VERSION)
        return NULL;
    return &g_sf_api;
}
