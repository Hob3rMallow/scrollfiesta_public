#include "pipeline_cube.h"

#include "../common/ves_platform.h"
#include "../common/qem.h"
#include "../common/mesh_trim.h"
#include "../common/dump_obj.h"
#include "../common/pipeline_constants.h"
#include "../extract/mesh_extract.h"
#include "../split/bridge_cut.h"
#include "../split/overlap_sep.h"
#include "../holefill/hole_fill.h"
#include "../remesh/orient_mesh.h"
#include "../remesh/pinhole_fill.h"
#include "../remesh/ball_pivot.h"
#include "../common/mls_project.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double now_sec(void) { return ves_clock_sec(); }
static double elapsed_since(double t0) { return now_sec() - t0; }

static int dump_stage(Arena_T arena,
                      const char *dump_dir, const char *cube_id,
                      const char *stage,
                      ComponentMesh *meshes, size_t n_meshes)
{
    if (!dump_dir || !cube_id || n_meshes == 0) return 0;
    char cube_dir[512];
    char stage_dir[1024];
    snprintf(cube_dir, sizeof(cube_dir), "%s/%s", dump_dir, cube_id);
    snprintf(stage_dir, sizeof(stage_dir), "%s/%s_%s",
             cube_dir, cube_id, stage);
    DumpObj_ensure_dir(dump_dir);
    DumpObj_ensure_dir(cube_dir);
    DumpObj_ensure_dir(stage_dir);
    DumpObj_write_meshes(arena, stage_dir, cube_id, stage,
                         meshes, n_meshes);
    return 0;
}

/* CDT/Liepa fill for the 4+ closed loops HoleFill_meshes leaves after the 3-loop
 * pass. Wraps src/holefill's HoleFill_process — this is the only TU where the
 * pipeline links Triangle/Clipper2. No pin mask: the cross-cube weld is now BPA
 * over the merged sheets, not a bit-exact seam join, so every cleanly-fillable
 * closed loop is filled wherever it sits. Tiny components skip CDT (cull). */
static int pipeline_cdt_fill(Arena_T arena, ComponentMesh *cm, void *user)
{
    (void)user;
    if (cm->nf < 100) return -1;
    size_t n_loops = 0, n_interior = 0, n_filled = 0;
    return HoleFill_process(arena, &cm->verts, &cm->faces, &cm->nv, &cm->nf,
                            NULL,
                            &n_loops, &n_interior, &n_filled);
}

/* Re-LOP + re-BPA a piece that the split stage carved out of a merged
 * multi-wrap component. After a split the verts sit at merged-LOP positions
 * (each was pulled toward the other wrap by the joint LOP) and the cut left a
 * ragged edge; projecting the now-isolated cloud (it sees only its own sheet's
 * neighbours, so it lands on its own surface) and rolling a fresh BPA gives a
 * clean single-sheet manifold. cell_origin {0,0,0}: this is a within-cube
 * re-surface of an interior sheet, so cross-cube hash alignment is moot.
 * Rewrites cm->verts/faces/nf in place and returns 0 on success; on any
 * failure returns -1 and leaves cm untouched (the split mesh is kept as-is). */
static int resurface_piece(Arena_T arena, ComponentMesh *cm, int n_threads)
{
    (void)n_threads;
    size_t nv = cm->nv;
    if (nv < 1000) return -1;            /* below the threading/benefit cutoff */
    const float zero_origin[3] = { 0.0f, 0.0f, 0.0f };
    float *pv = (float *)ARENA_ALLOC(arena, (long)(nv * 3 * sizeof(float)));
    float *pn = (float *)ARENA_ALLOC(arena, (long)(nv * 3 * sizeof(float)));
    float *sn = (float *)ARENA_ALLOC(arena, (long)(nv * 3 * sizeof(float)));
    float *sv = (MLS_PROJECT_ITERS >= 2)
        ? (float *)ARENA_ALLOC(arena, (long)(nv * 3 * sizeof(float))) : NULL;
    const float *src = cm->verts;
    float *dst = pv;
    for (int it = 0; it < MLS_PROJECT_ITERS; it++) {
        int is_last = (it == MLS_PROJECT_ITERS - 1);
        float *ndst = is_last ? pn : sn;
        MLS_project_verts(arena, src, nv, MLS_PROJECT_RADIUS_VOX,
                          zero_origin, dst, ndst);
        src = dst;
        dst = (dst == pv) ? sv : pv;
    }
    if (src != pv) memcpy(pv, src, nv * 3 * sizeof(float));
    int32_t *faces = NULL; size_t nf = 0;
    int brc = BallPivot_reconstruct(arena, pv, pn, nv, BPA_RHO_VOX, &faces, &nf);
    if (brc != 0 || nf == 0) return -1;  /* keep the split mesh on failure */
    cm->verts = pv;
    cm->faces = faces;
    cm->nf = nf;
    cm->self = cm;
    return 0;
}

int pipeline_process_cube(Arena_T arena,
                          const PipelineInput *in,
                          PipelineOutput *out)
{
    assert(arena && in && out);

    memset(out, 0, sizeof(*out));

    /* MeshExtract_run takes the cube-scoped dump dir (<dump_dir>/<cube_id>),
     * not the top-level dump dir. Materialize that path here and ensure
     * the directory exists. */
    char cube_dump_dir[1024] = {0};
    const char *mesh_dump_dir = NULL;
    if (in->dump_dir && in->cube_id) {
        snprintf(cube_dump_dir, sizeof(cube_dump_dir), "%s/%s",
                 in->dump_dir, in->cube_id);
        DumpObj_ensure_dir(in->dump_dir);
        DumpObj_ensure_dir(cube_dump_dir);
        mesh_dump_dir = cube_dump_dir;
    }

    /* ---- Step 0: MC + LOP + per-vert backface cull. ---- */
    double t0 = now_sec();
    int rc = MeshExtract_run(arena, in->tiff_path,
                             in->pred_dir,
                             in->halo_voxels,
                             (size_t)in->cube_D, (size_t)in->cube_H,
                             (size_t)in->cube_W, in->n_threads,
                             mesh_dump_dir,
                             (in->halo_voxels > 0 || in->dump_dir)
                                 ? in->cube_id : NULL,
                             in->skip_qem,
                             &out->meshes, &out->n_meshes);
    out->t_extract = elapsed_since(t0);
    if (rc != 0 || out->n_meshes == 0) {
        fprintf(stderr, "  pipeline_cube: Step 0 produced no meshes\n");
        return (rc != 0) ? rc : 1;
    }
    fprintf(stderr, "  Extract: %zu components (%.3fs)\n",
            out->n_meshes, out->t_extract);
    dump_stage(arena, in->dump_dir, in->cube_id, "step0",
               out->meshes, out->n_meshes);

    /* ---- Step 2: Bridge cut (oracle-gated; recursive Edmonds-Karp) ----
     * MC + per-vert-cull can leave two physically distinct sheets joined
     * by a thin bridge (nnUNet predicts a few stray voxels connecting two
     * wraps; MC then produces them as one connected component). The
     * Oracle's UV-grid raycast counts sheets per component; for any
     * multi-sheet component, BridgeCut runs vertex-split max-flow to find
     * the thinnest neck and severs it (with a CUT_GAP_DEPTH exclusion zone
     * so no ragged edges remain). Recursive up to BRIDGE_MAX_DEPTH levels.
     * Single-sheet components pass through unchanged (Oracle early-exit). */
    {
        double tb = now_sec();
        size_t total_in = out->n_meshes;
        size_t cap = total_in * 8 + 16;
        ComponentMesh *split_meshes = (ComponentMesh *)ARENA_ALLOC(
            arena, (long)(cap * sizeof(ComponentMesh)));
        size_t total_out = 0;
        size_t n_bridge = 0, n_ovl = 0, n_resurf = 0;
        for (size_t i = 0; i < total_in; i++) {
            ComponentMesh *cm = &out->meshes[i];
            ComponentMesh *subs = NULL;
            size_t n_subs = 0;
            int brc = BridgeCut_process(arena, cm,
                                         ORACLE_UV_GRID_SIZE,
                                         30.0,
                                         &subs, &n_subs);
            if (brc != 0 || n_subs == 0) {
                subs = cm;
                n_subs = 1;
            }
            if (n_subs > 1) n_bridge++;
            for (size_t j = 0; j < n_subs; j++) {
                /* Overlap separator (Limper lifted multicut): split the broad
                 * close-wrap overlaps that BridgeCut's thin-neck flow and the
                 * oracle's 10-vox gap both miss -- two wraps ~2-3 vox apart
                 * that BPA welded with many small triangles. Self-gates:
                 * returns the piece unchanged when it finds no projected
                 * triangle overlaps (a single sheet doesn't self-overlap along
                 * its own PCA normal). */
                ComponentMesh *osub = NULL;
                size_t n_osub = 0;
                int orc = OverlapSep_process(arena, &subs[j], 0,
                                             in->n_threads, 30.0,
                                             &osub, &n_osub);
                if (orc != 0 || n_osub == 0) { osub = &subs[j]; n_osub = 1; }
                if (n_osub > 1) n_ovl++;
                int piece_split = (n_subs > 1) || (n_osub > 1);
                for (size_t k = 0; k < n_osub; k++) {
                    if (osub[k].nf < 16) continue;   /* drop cut-zone slivers */
                    if (total_out >= cap) {
                        size_t new_cap = cap * 2;
                        ComponentMesh *grown = (ComponentMesh *)ARENA_ALLOC(
                            arena, (long)(new_cap * sizeof(ComponentMesh)));
                        memcpy(grown, split_meshes,
                               total_out * sizeof(ComponentMesh));
                        split_meshes = grown;
                        cap = new_cap;
                    }
                    split_meshes[total_out] = osub[k];
                    split_meshes[total_out].comp_id = (int)(total_out + 1);
                    split_meshes[total_out].pin_mask = NULL;
                    split_meshes[total_out].self = &split_meshes[total_out];
                    /* A split piece was surfaced from a merged point cloud and
                     * carries a ragged cut edge -- re-LOP + re-BPA it into a
                     * clean single-sheet surface (see resurface_piece). */
                    if (piece_split &&
                        resurface_piece(arena, &split_meshes[total_out],
                                        in->n_threads) == 0)
                        n_resurf++;
                    total_out++;
                }
            }
        }
        out->meshes = split_meshes;
        out->n_meshes = total_out;
        double t_b = elapsed_since(tb);
        fprintf(stderr,
            "  Split: %zu -> %zu components (%zu bridge-cut, %zu overlap, "
            "%zu re-surfaced, %.3fs)\n",
            total_in, total_out, n_bridge, n_ovl, n_resurf, t_b);

        dump_stage(arena, in->dump_dir, in->cube_id, "step2_bridge_cut",
                   out->meshes, out->n_meshes);
    }

    /* ---- Unified hole fill ----
     * One pass dispatched purely by loop size: pinch-split (bowtie verts) +
     * EXACT single-triangle (3-loop) fills (the CDT-free core), then every
     * remaining 4+ closed loop via CDT/Liepa (pipeline_cdt_fill over
     * src/holefill). No interior/exterior classification, no pin gate — a
     * closed loop is a hole to fill if cleanly fillable, else left open. The
     * cross-cube weld is a BPA pass over the merged sheets, so there is no seam
     * vert to protect: holes that graze a cube face are filled like any other. */
    {
        double tf = now_sec();
        /* Color the verts the fill adds (split copies + CDT Steiner) bright blue
         * in per-cube dumps (see dump_obj.c FILL_BLUE); cleared below for any
         * mesh that gained nothing. */
        for (size_t i = 0; i < out->n_meshes; i++) {
            out->meshes[i].nv_pre_fill = out->meshes[i].nv;
        }
        size_t pinch = 0, tri_filled = 0, tri_added = 0, skipped = 0,
               cdt_filled = 0;
        HoleFill_meshes(arena, out->meshes, out->n_meshes,
                        0 /* respect_pins: pins removed */,
                        pipeline_cdt_fill, NULL,
                        &pinch, &tri_filled, &tri_added, &skipped, &cdt_filled);
        for (size_t i = 0; i < out->n_meshes; i++) {
            ComponentMesh *cm = &out->meshes[i];
            if (cm->nv == cm->nv_pre_fill) { cm->nv_pre_fill = 0; }
        }
        double t_f = elapsed_since(tf);
        fprintf(stderr,
            "  HoleFill: %zu pinch splits, %zu 3-loop fills (+%zu tris), "
            "%zu loops skipped, %zu comps CDT-filled (%.3fs)\n",
            pinch, tri_filled, tri_added, skipped, cdt_filled, t_f);
        dump_stage(arena, in->dump_dir, in->cube_id, "holefill",
                   out->meshes, out->n_meshes);
    }

    /* ---- QEM simplification ---- */
    if (!in->skip_qem) {
        double tq = now_sec();
        float ratio = (in->qem_target_ratio > 0.0f)
                          ? in->qem_target_ratio
                          : QEM_TARGET_RATIO;
        size_t total_in_nv = 0, total_in_nf = 0;
        size_t total_out_nv = 0, total_out_nf = 0;
        size_t n_simplified = 0;
        for (size_t i = 0; i < out->n_meshes; i++) {
            ComponentMesh *cm = &out->meshes[i];
            if (cm->nf <= QEM_MIN_FACES_FOR_SIMPLIFY) continue;
            size_t target_nf = (size_t)((float)cm->nf * ratio);
            if (target_nf < 100) target_nf = 100;

            float *new_v = NULL;
            int32_t *new_f = NULL;
            size_t new_nv = 0, new_nf = 0;
            int qrc = QEM_simplify_pinned(arena, cm->verts, cm->nv,
                                          cm->faces, cm->nf,
                                          NULL,
                                          cm->pca_normal,
                                          target_nf,
                                          &new_v, &new_nv,
                                          &new_f, &new_nf,
                                          NULL);
            if (qrc == 0 && new_nf > 0) {
                total_in_nv += cm->nv;
                total_in_nf += cm->nf;
                total_out_nv += new_nv;
                total_out_nf += new_nf;
                cm->verts = new_v;
                cm->faces = new_f;
                cm->nv = new_nv;
                cm->nf = new_nf;
                n_simplified++;
            }
        }
        /* QEM edge collapses can locally mis-wind faces (its built-in
         * winding repair is incomplete), re-introducing same-direction
         * interior edges that the post-BPA orientation pass had removed.
         * Re-assert consistent winding per component. No per-vertex normals
         * survive QEM, so OrientMesh anchors each component to the majority
         * of its own faces — safe within a cube; the cross-cube seam sign was
         * already fixed at the post-BPA pass and preserved through collapse. */
        size_t total_qem_oflips = 0;
        for (size_t i = 0; i < out->n_meshes; i++) {
            ComponentMesh *cm = &out->meshes[i];
            size_t of = 0, oc = 0, orsd = 0;
            OrientMesh_consistent(arena, cm->verts, cm->nv, NULL,
                                  cm->faces, cm->nf, &of, &oc, &orsd);
            total_qem_oflips += of;
        }
        out->t_qem = elapsed_since(tq);
        fprintf(stderr,
            "  QEM: %zu/%zu comps simplified, %zu/%zu -> %zu/%zu v/f "
            "(orient %zu flips) (%.3fs)\n",
            n_simplified, out->n_meshes,
            total_in_nv, total_in_nf,
            total_out_nv, total_out_nf, total_qem_oflips, out->t_qem);
        dump_stage(arena, in->dump_dir, in->cube_id, "qem",
                   out->meshes, out->n_meshes);
    }

    /* ---- Trim to owned region (halo mode only) ---- */
    if (in->halo_voxels > 0) {
        double tt = now_sec();
        /* Inset the owned box by BPA_OWNED_TRIM_INSET so the FINAL per-cube mesh
         * (raw_snap = the weld input) never reaches a cube face — regardless of
         * which path produced it (extract, bridge-cut split, or resurface_piece,
         * which re-BPAs split pieces WITHOUT the step0 cloud inset). This is the
         * catch-all that guarantees adjacent cubes' charts don't touch: a grazing
         * wrap otherwise lands in both cubes at the shared plane and z-fight-
         * doubles. The 2*INSET gap is spanned by the seam bridge. */
        float owned_lo = (float)BPA_OWNED_TRIM_INSET;
        float owned_hi = (float)in->cube_D - (float)BPA_OWNED_TRIM_INSET;  /* cubic cubes */
        ComponentMesh *trimmed = (ComponentMesh *)ARENA_CALLOC(
            arena, (long)out->n_meshes, (long)sizeof(ComponentMesh));
        size_t n_trim_kept = 0;
        size_t total_in_nv = 0, total_in_nf = 0;
        size_t total_out_nv = 0, total_out_nf = 0;
        for (size_t i = 0; i < out->n_meshes; i++) {
            ComponentMesh *src = &out->meshes[i];
            if (src->nf == 0) continue;
            ComponentMesh *dst = &trimmed[n_trim_kept];
            float *t_v = NULL;
            int32_t *t_f = NULL;
            size_t t_nv = 0, t_nf = 0;
            int trc = Mesh_trim_to_owned_box(arena,
                src->verts, src->nv, src->faces, src->nf,
                NULL,
                owned_lo, owned_hi,
                &t_v, &t_nv, &t_f, &t_nf,
                NULL);
            if (trc != 0 || t_nf == 0) continue;
            dst->verts = t_v;
            dst->faces = t_f;
            dst->pin_mask = NULL;
            dst->vert_normals = NULL;
            dst->nv = t_nv;
            dst->nf = t_nf;
            dst->comp_id = src->comp_id;
            dst->pca_normal[0] = src->pca_normal[0];
            dst->pca_normal[1] = src->pca_normal[1];
            dst->pca_normal[2] = src->pca_normal[2];
            dst->centroid[0] = src->centroid[0];
            dst->centroid[1] = src->centroid[1];
            dst->centroid[2] = src->centroid[2];
            dst->nv_pre_fill = 0;
            dst->self = dst;
            total_in_nv += src->nv;
            total_in_nf += src->nf;
            total_out_nv += t_nv;
            total_out_nf += t_nf;
            n_trim_kept++;
        }
        out->trimmed = trimmed;
        out->n_trimmed = n_trim_kept;
        out->t_trim = elapsed_since(tt);
        fprintf(stderr,
            "  Trim: %zu/%zu comps kept, %zu/%zu -> %zu/%zu v/f (%.3fs)\n",
            n_trim_kept, out->n_meshes,
            total_in_nv, total_in_nf,
            total_out_nv, total_out_nf, out->t_trim);

        /* Dump under "raw_snap" stage name so grid_weld scans find the
         * trimmed output without modification. The raw_snap algorithm
         * itself is not run; the name is a directory convention. */
        double td = now_sec();
        dump_stage(arena, in->dump_dir, in->cube_id, "raw_snap",
                   trimmed, n_trim_kept);
        out->t_dump = elapsed_since(td);
    } else {
        /* No halo: trimmed == meshes by reference. */
        out->trimmed = out->meshes;
        out->n_trimmed = out->n_meshes;
        out->t_trim = 0.0;
    }

    return 0;
}
