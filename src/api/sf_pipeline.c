/*
 * sf_pipeline.c — the volume->mesh pipeline entry of the public C API
 * (sf_pipeline_run over pipeline_process_cube's in-memory path), plus the
 * sf_weld placeholder (experimental; unimplemented in 0.9).
 */
#include "sf_internal.h"

#include "../pipeline/pipeline_cube.h"

#include <stdio.h>
#include <string.h>

/* Emit a ComponentMesh array as a malloc'd sf_mesh_list without provenance
 * (reconstruction outputs have none). Never raises. */
static sf_status emit_cm_list(const ComponentMesh *ms, size_t n,
                              sf_mesh_list *out)
{
    sf_mesh_list lst;
    memset(&lst, 0, sizeof lst);
    memset(out, 0, sizeof *out);
    if (n == 0)
        return SF_OK;                 /* empty list is a valid result */
    lst.items = sf_malloc(n * sizeof *lst.items);
    if (!lst.items)
        return SF_ERROR_OOM;
    memset(lst.items, 0, n * sizeof *lst.items);
    lst.count = n;
    for (size_t i = 0; i < n; i++) {
        sf_status rc = sf_mesh_out_cm(&lst.items[i], &ms[i], NULL);
        if (rc != SF_OK) {
            sf_mesh_list_free(&lst);
            return rc;
        }
    }
    *out = lst;
    return SF_OK;
}

SF_API sf_status sf_pipeline_run(const sf_volume *vol,
                                 const sf_pipeline_config *cfg,
                                 sf_mesh_list *out_full,
                                 sf_mesh_list *out_trimmed,
                                 sf_pipeline_report *rep)
{
    SfOpGuard          g;
    volatile sf_status rc;
    sf_pipeline_config def;

    if (!vol || !vol->data || !out_full)
        return SF_ERROR_BAD_ARG;
    if (!cfg) { def = sf_pipeline_config_default(); cfg = &def; }
    if (cfg->cube_size <= 0 || cfg->halo_voxels < 0)
        return SF_ERROR_BAD_ARG;
    {
        int expect = cfg->cube_size + 2 * cfg->halo_voxels;
        if (vol->nx != expect || vol->ny != expect || vol->nz != expect)
            return SF_ERROR_BAD_ARG;
    }
    memset(out_full, 0, sizeof *out_full);
    if (out_trimmed) memset(out_trimmed, 0, sizeof *out_trimmed);
    if (rep) memset(rep, 0, sizeof *rep);

    rc = sf_op_begin(&g, &cfg->common);
    if (rc != SF_OK)
        return rc;

    /* Extend the per-call tuning table with the stage toggles (the modules
     * read them through sf_env). User-provided entries come first so they
     * win the linear lookup. */
    {
        size_t   n_user = g.ctx.n_tuning;
        SfRunKV *merged = NULL;

        TRY {
            merged = ARENA_ALLOC(g.arena, (n_user + 3) * sizeof *merged);
        }
        EXCEPT(Arena_Failed)
            rc = SF_ERROR_OOM;
        END_TRY;

        if (rc != SF_OK) {
            sf_op_end(&g);
            return rc;
        }
        if (n_user)
            memcpy(merged, g.ctx.tuning, n_user * sizeof *merged);
        size_t k = n_user;
        if (!cfg->enable_depth_peel) {
            merged[k].key = "VES_DEPTHPEEL_OFF"; merged[k].value = "1"; k++;
        }
        if (!cfg->enable_dev_cut) {
            merged[k].key = "VES_DEVCUT_OFF"; merged[k].value = "1"; k++;
        }
        if (!cfg->enable_sever) {
            merged[k].key = "VES_SEVER_OFF"; merged[k].value = "1"; k++;
        }
        g.ctx.tuning   = merged;
        g.ctx.n_tuning = k;
    }

    TRY {
        char cube_id[64];
        snprintf(cube_id, sizeof cube_id, "z%05lld_y%05lld_x%05lld",
                 (long long)cfg->origin_xyz[2],
                 (long long)cfg->origin_xyz[1],
                 (long long)cfg->origin_xyz[0]);

        PipelineInput  pin;
        PipelineOutput pout;
        memset(&pin, 0, sizeof pin);
        memset(&pout, 0, sizeof pout);

        pin.tiff_path        = NULL;
        pin.pred_dir         = NULL;
        pin.cube_id          = cube_id;
        pin.halo_voxels      = cfg->halo_voxels;
        pin.cube_D           = cfg->cube_size;
        pin.cube_H           = cfg->cube_size;
        pin.cube_W           = cfg->cube_size;
        pin.n_threads        = cfg->common.n_threads > 0 ? cfg->common.n_threads : 1;
        pin.qem_target_ratio = cfg->qem_target_ratio;
        pin.dump_dir         = cfg->dump_dir;
        pin.skip_qem         = cfg->skip_qem;
        pin.vol_in           = vol->data;
        pin.p_size_in        = vol->nx;
        pin.cube_origin_zyx[0] = cfg->origin_xyz[2];
        pin.cube_origin_zyx[1] = cfg->origin_xyz[1];
        pin.cube_origin_zyx[2] = cfg->origin_xyz[0];

        RunCtx_progress_check("pipeline", 0.0);

        if (pipeline_process_cube(g.arena, &pin, &pout) != 0) {
            rc = SF_ERROR;
        } else {
            if (rep) {
                rep->t_extract    = pout.t_extract;
                rep->t_qem        = pout.t_qem;
                rep->t_trim       = pout.t_trim;
                rep->t_dump       = pout.t_dump;
                rep->n_bad_sheets = pout.n_bad_sheets;
            }
            RunCtx_progress_check("pipeline", 0.95);

            rc = emit_cm_list(pout.meshes, pout.n_meshes, out_full);
            if (rc == SF_OK && out_trimmed)
                rc = emit_cm_list(pout.trimmed, pout.n_trimmed, out_trimmed);
        }
    }
    SF_EXCEPT_TAIL(g, rc)
    END_TRY;

    if (rc != SF_OK) {
        sf_mesh_list_free(out_full);
        if (out_trimmed)
            sf_mesh_list_free(out_trimmed);
    }
    sf_op_end(&g);
    return rc;
}

/* Experimental — composition of the weld-side modules is planned but not yet
 * wired; the sf_api table carries NULL for this entry in 0.9 and the direct
 * symbol reports SF_ERROR_UNSUPPORTED. */
SF_API sf_status sf_weld(const sf_mesh *meshes, size_t n_meshes,
                         const sf_weld_config *cfg,
                         sf_mesh *out, sf_weld_report *rep)
{
    (void)meshes; (void)n_meshes; (void)cfg;
    if (out) memset(out, 0, sizeof *out);
    if (rep) memset(rep, 0, sizeof *rep);
    return SF_ERROR_UNSUPPORTED;
}
