/*
 * sf_pipeline.c — the volume->mesh pipeline entry of the public C API
 * (sf_pipeline_run over pipeline_process_cube's in-memory path), plus the
 * sf_weld placeholder (experimental; unimplemented in 0.9).
 */
#include "sf_internal.h"

#include "../common/pca.h"
#include "../common/vert_weld.h"
#include "../holefill/hole_fill.h"
#include "../pipeline/pipeline_cube.h"
#include "../remesh/manifold_guard.h"
#include "../remesh/orient_weld.h"
#include "../remesh/pinhole_fill.h"
#include "../remesh/seam_weld.h"
#include "../remesh/weld_cleanup.h"

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

/* ── sf_weld ────────────────────────────────────────────────────────────────
 * Composition of the library weld modules in grid_weld's stage order:
 * concat -> Weld_verts (duplicate-vertex merge across the halo overlap) ->
 * face dedup -> SeamWeld_bridge (BPA across cube-boundary planes) ->
 * PinholeFill -> interior HoleFill -> OrientWeld_components ->
 * WeldCleanup -> ManifoldGuard(reorient=0). The grid_weld TOOL additionally
 * has env-gated experimental stages (micro-weld, lambda gate, seam-fill,
 * phase-sever) that this call deliberately does not run.
 * ── face dedup ── after the vertex merge, faces duplicated in the halo
 * overlap become identical index triples (up to rotation/winding); keep the
 * first occurrence of each unordered triple. */

typedef struct {
    int32_t a, b, c;   /* sorted ascending */
    size_t  idx;       /* original face index */
} SfFaceKey;

static int face_key_cmp(const void *pa, const void *pb)
{
    const SfFaceKey *x = (const SfFaceKey *)pa;
    const SfFaceKey *y = (const SfFaceKey *)pb;
    if (x->a != y->a) return x->a < y->a ? -1 : 1;
    if (x->b != y->b) return x->b < y->b ? -1 : 1;
    if (x->c != y->c) return x->c < y->c ? -1 : 1;
    if (x->idx != y->idx) return x->idx < y->idx ? -1 : 1;
    return 0;
}

/* Compacts faces[] in place (first occurrence kept, input order preserved);
 * returns the new face count. */
static size_t sf_dedup_faces(Arena_T arena, int32_t *faces, size_t nf)
{
    if (nf < 2)
        return nf;
    Arena_Mark mark = Arena_save(arena);
    SfFaceKey *keys = ARENA_ALLOC(arena, nf * sizeof *keys);
    uint8_t   *drop = ARENA_CALLOC(arena, nf, 1);
    for (size_t i = 0; i < nf; i++) {
        int32_t v0 = faces[i * 3 + 0];
        int32_t v1 = faces[i * 3 + 1];
        int32_t v2 = faces[i * 3 + 2];
        int32_t t;
        if (v0 > v1) { t = v0; v0 = v1; v1 = t; }
        if (v1 > v2) { t = v1; v1 = v2; v2 = t; }
        if (v0 > v1) { t = v0; v0 = v1; v1 = t; }
        keys[i].a = v0; keys[i].b = v1; keys[i].c = v2; keys[i].idx = i;
    }
    qsort(keys, nf, sizeof *keys, face_key_cmp);
    for (size_t i = 1; i < nf; i++) {
        if (keys[i].a == keys[i - 1].a && keys[i].b == keys[i - 1].b &&
            keys[i].c == keys[i - 1].c)
            drop[keys[i].idx] = 1;   /* sorted ties break by idx: keep first */
    }
    size_t keep = 0;
    for (size_t i = 0; i < nf; i++) {
        if (drop[i])
            continue;
        if (keep != i) {
            faces[keep * 3 + 0] = faces[i * 3 + 0];
            faces[keep * 3 + 1] = faces[i * 3 + 1];
            faces[keep * 3 + 2] = faces[i * 3 + 2];
        }
        keep++;
    }
    Arena_restore(arena, mark);
    return keep;
}

SF_API sf_status sf_weld(const sf_mesh *meshes, size_t n_meshes,
                         const sf_weld_config *cfg,
                         sf_mesh *out, sf_weld_report *rep)
{
    SfOpGuard          g;
    volatile sf_status rc;
    sf_weld_config     def;

    if (!meshes || n_meshes == 0 || !out)
        return SF_ERROR_BAD_ARG;
    for (size_t i = 0; i < n_meshes; i++) {
        if (sf_mesh_validate(&meshes[i], NULL, 0) != SF_OK)
            return SF_ERROR_BAD_ARG;
    }
    memset(out, 0, sizeof *out);
    if (rep) memset(rep, 0, sizeof *rep);
    if (!cfg) { def = sf_weld_config_default(); cfg = &def; }

    rc = sf_op_begin(&g, &cfg->common);
    if (rc != SF_OK)
        return rc;

    /* rho_max reaches the bridge module through its documented env knob. */
    {
        char rho_max_buf[32];
        snprintf(rho_max_buf, sizeof rho_max_buf, "%g",
                 cfg->rho_max > 0.f ? (double)cfg->rho_max : 3.0);
        size_t   n_user = g.ctx.n_tuning;
        SfRunKV *merged = NULL;
        TRY
            merged = ARENA_ALLOC(g.arena, (n_user + 1) * sizeof *merged);
        EXCEPT(Arena_Failed)
            rc = SF_ERROR_OOM;
        END_TRY;
        if (rc != SF_OK) {
            sf_op_end(&g);
            return rc;
        }
        if (n_user)
            memcpy(merged, g.ctx.tuning, n_user * sizeof *merged);
        /* the buffer must outlive the call: arena-copy the value */
        char *val = NULL;
        TRY {
            val = ARENA_ALLOC(g.arena, strlen(rho_max_buf) + 1);
            strcpy(val, rho_max_buf);
        }
        EXCEPT(Arena_Failed)
            rc = SF_ERROR_OOM;
        END_TRY;
        if (rc != SF_OK) {
            sf_op_end(&g);
            return rc;
        }
        merged[n_user].key   = "SEAM_RHO_MAX";
        merged[n_user].value = val;
        g.ctx.tuning   = merged;
        g.ctx.n_tuning = n_user + 1;
    }

    TRY {
        /* 1. Concatenate all inputs into one (z,y,x) soup. */
        size_t cat_nv = 0, cat_nf = 0;
        for (size_t i = 0; i < n_meshes; i++) {
            cat_nv += meshes[i].n_vertices;
            cat_nf += meshes[i].n_faces;
        }
        float   *cat_verts = ARENA_ALLOC(g.arena, cat_nv * 3 * sizeof *cat_verts);
        int32_t *cat_faces = ARENA_ALLOC(g.arena, cat_nf * 3 * sizeof *cat_faces);
        {
            size_t voff = 0, foff = 0;
            for (size_t i = 0; i < n_meshes; i++) {
                const sf_mesh *m = &meshes[i];
                for (size_t v = 0; v < m->n_vertices; v++) {
                    cat_verts[(voff + v) * 3 + 0] = m->vertices[v * 3 + 2];
                    cat_verts[(voff + v) * 3 + 1] = m->vertices[v * 3 + 1];
                    cat_verts[(voff + v) * 3 + 2] = m->vertices[v * 3 + 0];
                }
                for (size_t f = 0; f < m->n_faces; f++) {
                    cat_faces[(foff + f) * 3 + 0] =
                        m->faces[f * 3 + 0] + (int32_t)voff;
                    cat_faces[(foff + f) * 3 + 1] =
                        m->faces[f * 3 + 2] + (int32_t)voff;
                    cat_faces[(foff + f) * 3 + 2] =
                        m->faces[f * 3 + 1] + (int32_t)voff;
                }
                voff += m->n_vertices;
                foff += m->n_faces;
            }
        }
        RunCtx_progress_check("weld", 0.05);

        /* 2. Merge duplicate vertices (halo overlap), drop degenerates. */
        float  *wverts = NULL;
        size_t  wnv = 0, wnf = cat_nf;
        Weld_verts(g.arena, cat_verts, cat_nv, NULL, cat_faces, cat_nf, &wnf,
                   cfg->weld_eps > 0.f ? cfg->weld_eps : 1e-4f,
                   true, &wverts, &wnv, NULL);
        if (rep)
            rep->verts_welded = cat_nv - wnv;

        /* 3. Drop faces duplicated by the overlap. */
        wnf = sf_dedup_faces(g.arena, cat_faces, wnf);
        RunCtx_progress_check("weld", 0.2);

        /* 4. BPA bridge across cube-boundary seam planes. */
        int32_t *bfaces = cat_faces;
        size_t   bnf = wnf, n_bridge = 0;
        if (SeamWeld_bridge(g.arena, wverts, wnv, cat_faces, wnf,
                            cfg->cube_size > 0.f ? cfg->cube_size : 128.f,
                            cfg->rho > 0.f ? cfg->rho : 1.5f,
                            cfg->band > 0.f ? cfg->band : 6.0f,
                            &bfaces, &bnf, &n_bridge) != 0) {
            rc = SF_ERROR;
        }
        if (rep)
            rep->bridge_faces = n_bridge;
        RunCtx_progress_check("weld", 0.4);

        ComponentMesh cm;
        memset(&cm, 0, sizeof cm);
        if (rc == SF_OK) {
            cm.verts = wverts;
            cm.faces = bfaces;
            cm.nv = wnv;
            cm.nf = bnf;
            cm.comp_id = 1;
            if (PCA_normal(cm.verts, cm.nv, cm.pca_normal, cm.centroid) != 0) {
                cm.pca_normal[0] = 1.0f;
                cm.pca_normal[1] = 0.0f;
                cm.pca_normal[2] = 0.0f;
            }
            cm.self = &cm;
        }

        /* 5+6. Close pinholes, then geometrically-interior holes. */
        if (rc == SF_OK && cfg->fill_holes) {
            size_t splits = 0, filled = 0, tris = 0, skipped = 0;
            if (PinholeFill_process(g.arena, &cm, 1, 0,
                                    &splits, &filled, &tris, &skipped) != 0)
                rc = SF_ERROR;
            else if (rep)
                rep->holes_filled += filled;
            RunCtx_progress_check("weld", 0.55);

            if (rc == SF_OK) {
                float   *hv = cm.verts;
                int32_t *hf = cm.faces;
                size_t   hnv = cm.nv, hnf = cm.nf;
                size_t   loops = 0, interior = 0, filled2 = 0;
                if (HoleFill_process_ex(g.arena, &hv, &hf, &hnv, &hnf, NULL,
                                        1, &loops, &interior, &filled2) != 0) {
                    rc = SF_ERROR;
                } else {
                    cm.verts = hv;
                    cm.faces = hf;
                    cm.nv = hnv;
                    cm.nf = hnf;
                    if (rep)
                        rep->holes_filled += filled2;
                }
            }
        }
        RunCtx_progress_check("weld", 0.7);

        /* 7. Cross-component orientation vote. */
        if (rc == SF_OK) {
            size_t flipped = 0;
            if (OrientWeld_components(g.arena, cm.verts, cm.nv,
                                      cm.faces, cm.nf, 3.0f, &flipped) != 0)
                rc = SF_ERROR;
        }

        /* 8+9. Sliver cleanup + manifold guard (winding preserved). */
        if (rc == SF_OK && cfg->cleanup) {
            WeldCleanupParams p;
            WeldCleanup_default_params(&p);
            if (WeldCleanup_process(g.arena, &cm, &p, NULL) != 0)
                rc = SF_ERROR;
            RunCtx_progress_check("weld", 0.85);
            if (rc == SF_OK) {
                if (ManifoldGuard_process(g.arena, &cm, 1, 0, NULL) != 0)
                    rc = SF_ERROR;
            }
        }

        if (rc == SF_OK) {
            RunCtx_progress_check("weld", 1.0);
            /* Reconstruction-class output: no per-vertex provenance. */
            rc = sf_mesh_out_cm(out, &cm, NULL);
        }
    }
    SF_EXCEPT_TAIL(g, rc)
    END_TRY;

    if (rc != SF_OK)
        sf_mesh_free(out);
    sf_op_end(&g);
    return rc;
}
