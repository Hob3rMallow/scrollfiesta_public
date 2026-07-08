#include "pipeline_cube.h"

#include "../common/ves_platform.h"
#include "../common/qem.h"
#include "../common/mesh_trim.h"
#include "../common/dump_obj.h"
#include "../common/pipeline_constants.h"
#include "../extract/mesh_extract.h"
#include "../extract/mesh_resplit.h"
#include "../split/depth_peel.h"
#include "../split/dev_cut.h"
#include "../split/bridge_cut.h"
#include "../split/overlap_sep.h"
#include "../holefill/hole_fill.h"
#include "../remesh/orient_mesh.h"
#include "../remesh/pinhole_fill.h"
#include "../remesh/ball_pivot.h"
#include "../remesh/manifold_guard.h"
#include "../remesh/component_cull.h"
#include "../remesh/dev_gate.h"
#include "../common/mesh_manifold.h"
#include "../common/mls_project.h"
#include "../topology/seam_cut.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double now_sec(void) { return ves_clock_sec(); }
static double elapsed_since(double t0) { return now_sec() - t0; }

/* Read an env var as a double, falling back to `dflt` when unset/unparsable.
 * Lets the Step 1c developability-cut knobs be swept without a rebuild. */
static double env_d(const char *name, double dflt)
{
    const char *s = getenv(name);
    if (!s || !*s) return dflt;
    char *end = NULL;
    double v = strtod(s, &end);
    return (end && end != s) ? v : dflt;
}

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

/* Dump the per-component point clouds (the LOP/MLS images) of a stage as point
 * OBJs, one per component, into <dump_dir>/<cube_id>/<cube_id>_<stage>/. World-
 * shifted via cube_id (same convention as dump_stage). */
static void dump_cloud_stage(Arena_T arena,
                             const char *dump_dir, const char *cube_id,
                             const char *stage,
                             const MeshResplitCloud *clouds, size_t n)
{
    if (!dump_dir || !cube_id || !clouds || n == 0) return;
    char cube_dir[512], stage_dir[1024];
    snprintf(cube_dir, sizeof cube_dir, "%s/%s", dump_dir, cube_id);
    snprintf(stage_dir, sizeof stage_dir, "%s/%s_%s", cube_dir, cube_id, stage);
    DumpObj_ensure_dir(dump_dir);
    DumpObj_ensure_dir(cube_dir);
    DumpObj_ensure_dir(stage_dir);
    for (size_t i = 0; i < n; i++) {
        if (clouds[i].n == 0) continue;
        char path[1280];
        /* post-MLS: the LOP/MLS image that feeds BPA. */
        if (clouds[i].lop_pts) {
            snprintf(path, sizeof path, "%s/%s_%s_%03zu.obj",
                     stage_dir, cube_id, stage, i);
            DumpObj_write_points_world(arena, path, cube_id,
                                       clouds[i].lop_pts, NULL, clouds[i].n);
        }
        /* pre-MLS: the raw voxel-center cloud before LOP (sibling _pre file). */
        if (clouds[i].orig_pts) {
            snprintf(path, sizeof path, "%s/%s_%s_pre%03zu.obj",
                     stage_dir, cube_id, stage, i);
            DumpObj_write_points_world(arena, path, cube_id,
                                       clouds[i].orig_pts, NULL, clouds[i].n);
        }
    }
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
    MeshResplitCloud *clouds = NULL;
    int rc = MeshExtract_run(arena, in->tiff_path,
                             in->pred_dir,
                             in->halo_voxels,
                             (size_t)in->cube_D, (size_t)in->cube_H,
                             (size_t)in->cube_W, in->n_threads,
                             mesh_dump_dir,
                             (in->halo_voxels > 0 || in->dump_dir)
                                 ? in->cube_id : NULL,
                             in->skip_qem,
                             in->vol_in, in->p_size_in,
                             in->vol_in ? in->cube_origin_zyx : NULL,
                             &out->meshes, &out->n_meshes, &clouds);
    out->t_extract = elapsed_since(t0);
    if (rc != 0 || out->n_meshes == 0) {
        fprintf(stderr, "  pipeline_cube: Step 0 produced no meshes\n");
        return (rc != 0) ? rc : 1;
    }
    fprintf(stderr, "  Extract: %zu components (%.3fs)\n",
            out->n_meshes, out->t_extract);
    dump_cloud_stage(arena, in->dump_dir, in->cube_id, "step0_mls",
                     clouds, out->n_meshes);
    dump_stage(arena, in->dump_dir, in->cube_id, "step1_bpa",
               out->meshes, out->n_meshes);

    /* ---- Step 1b: connectivity re-split. A voxel-CC can hold several
     * disconnected surface sheets surfaced together; the merged LOP leaves
     * cross-sheet artifacts that make Step 2 over-split. Separate each into its
     * real connectivity-components and re-LOP/re-BPA each from the ORIGINAL
     * voxel-center points so Step 2 sees clean single sheets. The per-output
     * labelled point clouds (rclouds, index-aligned with out->meshes) flow into
     * Step 2 so the bridge-cut/overlap output can be re-LOP'd from originals too
     * (Step 4). ---- */
    MeshResplitCloud *rclouds = NULL;
    {
        double tr = now_sec();
        char resplit_dump_dir[1024] = {0};
        MeshResplitDump rdump = {0};
        if (in->dump_dir && in->cube_id) {
            snprintf(resplit_dump_dir, sizeof resplit_dump_dir, "%s/%s",
                     in->dump_dir, in->cube_id);
            rdump.dir = resplit_dump_dir; rdump.cube_id = in->cube_id;
            rdump.cc_stage = "step2_cc"; rdump.mls_stage = "step3_mls";
        }
        ComponentMesh *rs = NULL; size_t n_rs = 0;
        if (MeshResplit_run(arena, out->meshes, out->n_meshes, clouds,
                            in->n_threads, in->dump_dir ? &rdump : NULL,
                            &rs, &n_rs, &rclouds) == 0 && n_rs > 0) {
            if (n_rs != out->n_meshes)
                fprintf(stderr, "  Resplit: %zu -> %zu components (%.3fs)\n",
                        out->n_meshes, n_rs, elapsed_since(tr));
            out->meshes = rs;
            out->n_meshes = n_rs;
        } else {
            rclouds = NULL;   /* keep rclouds aligned with out->meshes */
        }
        dump_stage(arena, in->dump_dir, in->cube_id, "step4_bpa",
                   out->meshes, out->n_meshes);
    }

    /* ---- Step 1b-peel: stacked-wrap separator (src/split/depth_peel.c) ----
     * BPA (rho=1.2, no inter-wrap clearance) can fuse papyrus wraps ~2-3 vox apart
     * into ONE component -- a stack of near-parallel sheets the connectivity
     * re-split above cannot break (it is one component) and developability cannot
     * see (each wrap is locally flat, lambda ~ 0). The cheap separable signal is
     * DEPTH along the surface normal: only an inter-wrap bridge edge jumps across
     * the gap. DepthPeel_process cuts those edges and returns the connected layers;
     * a single sheet (however folded) has no such jump and passes through (the
     * self-gate). Each peeled layer is re-LOP'd from its OWN verts (resurface_own --
     * NOT re-labelling the parent cloud, whose 1.5-vox margin would re-merge the
     * 2-3 vox layers). Whatever DepthPeel can't cleanly bin (a strongly curved
     * stack) falls through to Step 2's overlap-sep unchanged. The wall-guard in BPA
     * prevents many fusions upstream; this catches the compacted ones it can't.
     * Env-overridable; VES_DEPTHPEEL_OFF disables the stage. ---- */
    if (!getenv("VES_DEPTHPEEL_OFF")) {
        double tdp = now_sec();
        double dp_gap = env_d("DEPTH_PEEL_MIN_GAP_VOX", DEPTH_PEEL_MIN_GAP_VOX);
        double dp_dil = env_d("DEPTH_PEEL_GAP_DEPTH",   DEPTH_PEEL_GAP_DEPTH);

        size_t n_in = out->n_meshes;
        size_t cap  = n_in * 4 + 16, cnt = 0;
        ComponentMesh    *acc  = (ComponentMesh *)ARENA_ALLOC(
            arena, (long)(cap * sizeof(ComponentMesh)));
        MeshResplitCloud *accc = (MeshResplitCloud *)ARENA_ALLOC(
            arena, (long)(cap * sizeof(MeshResplitCloud)));
        MeshResplitCloud zerocl; memset(&zerocl, 0, sizeof zerocl);
        size_t n_peeled = 0, n_layers = 0;

        #define DP_ENSURE() do { if (cnt >= cap) { size_t nc = cap * 2;          \
            ComponentMesh *na = (ComponentMesh *)ARENA_ALLOC(arena,              \
                (long)(nc * sizeof(ComponentMesh)));                            \
            memcpy(na, acc, cnt * sizeof(ComponentMesh));                        \
            for (size_t _i = 0; _i < cnt; _i++) na[_i].self = &na[_i]; acc = na; \
            MeshResplitCloud *nq = (MeshResplitCloud *)ARENA_ALLOC(arena,        \
                (long)(nc * sizeof(MeshResplitCloud)));                         \
            memcpy(nq, accc, cnt * sizeof(MeshResplitCloud)); accc = nq;         \
            cap = nc; } } while (0)
        #define DP_EMIT(MSH, CLD) do { DP_ENSURE();                             \
            acc[cnt] = (MSH); acc[cnt].self = &acc[cnt];                        \
            acc[cnt].comp_id = (int)(cnt + 1);                                  \
            accc[cnt] = (CLD); cnt++; } while (0)

        for (size_t i = 0; i < n_in; i++) {
            ComponentMesh    *cm = &out->meshes[i];
            MeshResplitCloud *cl = rclouds ? &rclouds[i] : NULL;

            ComponentMesh *pieces = NULL; size_t n_pc = 0;
            int split = 0;
            if (DepthPeel_process(arena, cm, dp_gap, dp_dil,
                                  DEPTH_PEEL_MIN_COMP_VERTS, &pieces, &n_pc) == 0
                && n_pc > 1) {
                /* re-LOP each peeled layer from its OWN verts (clean the ragged cut
                 * boundary); the layers are already separate ComponentMeshes, so
                 * hole-fill cannot re-bridge them. Emit with an empty cloud (n=0) so
                 * Step 4's re-LOP -- which would relabel a parent cloud and re-merge
                 * the close layers -- skips them, exactly as DevCut does for its own
                 * pieces. */
                for (size_t k = 0; k < n_pc; k++) {
                    MeshResplit_resurface_own(arena, &pieces[k],
                        (cl && cl->n > 0) ? cl->cell_origin : NULL,
                        MLS_RESPLIT_ITERS);
                    DP_EMIT(pieces[k], zerocl);
                }
                n_peeled++; n_layers += n_pc; split = 1;
                fprintf(stderr, "    DepthPeel: comp %zu -> %zu layers\n", i, n_pc);
            }
            if (!split)
                DP_EMIT(*cm, cl ? *cl : zerocl);
        }
        #undef DP_EMIT
        #undef DP_ENSURE

        out->meshes = acc; out->n_meshes = cnt;
        rclouds = accc;
        fprintf(stderr,
            "  DepthPeel: %zu -> %zu components (%zu stack(s) peeled into %zu "
            "layers, %.3fs)\n", n_in, cnt, n_peeled, n_layers, elapsed_since(tdp));
        dump_stage(arena, in->dump_dir, in->cube_id, "step1b_peel",
                   out->meshes, out->n_meshes);
    }

    /* ---- Step 1c: developability cut (Stein/Grinspun/Crane 2018, sec 4.4) ----
     * The FIRST splitter, ahead of bridge-cut/overlap. Per component, threshold
     * the per-vertex Crane energy lambda_v on the raw mesh and cut THROUGH the
     * non-developable seam (DevCut_process). A clean wrap is developable
     * (lambda ~ 0, like a rolled cylinder) so it has no spanning seam and passes
     * through untouched -- the "NO intra-sheet split" self-gate. When a component
     * does split, re-LOP each new piece from its OWN labelled points (re-MLS +
     * re-BPA via MeshResplit_remesh_pieces) and KEEP the split only if every
     * re-LOP'd piece is measurably more developable than the parent (the
     * lambda-energy gate). Components left un-split (or gate-rejected) fall
     * through to Step 2's bridge-cut/overlap unchanged -- which also run on the
     * accepted pieces (they self-gate to no-ops on clean sheets). Knobs are
     * env-overridable for calibration; VES_DEVCUT_OFF disables the stage. */
    if (!getenv("VES_DEVCUT_OFF")) {
        double tdc = now_sec();
        double dc_eps    = env_d("DEV_CUT_EPS",           DEV_CUT_EPS);
        double dc_gap    = env_d("DEV_CUT_GAP_DEPTH",     DEV_CUT_GAP_DEPTH);
        double dc_margin = env_d("DEV_CUT_GATE_MARGIN",   DEV_CUT_GATE_MARGIN);
        double dc_ceil   = env_d("DEV_CUT_GATE_ABS_CEIL", DEV_CUT_GATE_ABS_CEIL);

        size_t n_in = out->n_meshes;
        size_t cap  = n_in * 4 + 16, cnt = 0;
        ComponentMesh    *acc  = (ComponentMesh *)ARENA_ALLOC(
            arena, (long)(cap * sizeof(ComponentMesh)));
        MeshResplitCloud *accc = (MeshResplitCloud *)ARENA_ALLOC(
            arena, (long)(cap * sizeof(MeshResplitCloud)));
        MeshResplitCloud zerocl; memset(&zerocl, 0, sizeof zerocl);
        size_t n_cut = 0, n_gate_rej = 0;

        /* Opt-in calibration: log the per-cube developability distribution over
         * ALL input components (not just those that form a cuttable seam), so a
         * sweep shows whether a region is genuinely non-developable and how many
         * components clear the gate margin. Env DEV_CUT_CAL=1. */
        if (getenv("DEV_CUT_CAL")) {
            double mx = 0.0, sm = 0.0; size_t nge = 0;
            for (size_t i = 0; i < n_in; i++) {
                double pf = DevCut_nondev_fraction(arena, &out->meshes[i], dc_eps, NULL);
                if (pf > mx) mx = pf;
                sm += pf;
                if (pf >= dc_margin) nge++;
            }
            fprintf(stderr,
                "  DevCut-cal %s: ncomp=%zu max_nondev=%.4f mean_nondev=%.4f "
                "n_ge_margin=%zu (eps=%.3f margin=%.3f)\n",
                in->cube_id ? in->cube_id : "?", n_in, mx,
                n_in ? sm / (double)n_in : 0.0, nge, dc_eps, dc_margin);
        }

        #define DC_ENSURE() do { if (cnt >= cap) { size_t nc = cap * 2;          \
            ComponentMesh *na = (ComponentMesh *)ARENA_ALLOC(arena,              \
                (long)(nc * sizeof(ComponentMesh)));                            \
            memcpy(na, acc, cnt * sizeof(ComponentMesh));                        \
            for (size_t _i = 0; _i < cnt; _i++) na[_i].self = &na[_i]; acc = na; \
            MeshResplitCloud *nq = (MeshResplitCloud *)ARENA_ALLOC(arena,        \
                (long)(nc * sizeof(MeshResplitCloud)));                         \
            memcpy(nq, accc, cnt * sizeof(MeshResplitCloud)); accc = nq;         \
            cap = nc; } } while (0)
        #define DC_EMIT(MSH, CLD) do { DC_ENSURE();                             \
            acc[cnt] = (MSH); acc[cnt].self = &acc[cnt];                        \
            acc[cnt].comp_id = (int)(cnt + 1);                                  \
            accc[cnt] = (CLD); cnt++; } while (0)

        for (size_t i = 0; i < n_in; i++) {
            ComponentMesh    *cm = &out->meshes[i];
            MeshResplitCloud *cl = rclouds ? &rclouds[i] : NULL;

            ComponentMesh *pieces = NULL; size_t n_pc = 0;
            int accepted = 0;
            if (DevCut_process(arena, cm, dc_eps, dc_gap, DEV_CUT_MIN_COMP_VERTS,
                               &pieces, &n_pc) == 0 && n_pc > 1) {
                /* The accept test is worst_piece <= parent_nondev - margin, and a
                 * piece fraction is always >= 0, so a parent that is ALREADY more
                 * developable than the margin can never clear it. Measure the
                 * parent first (cheap -- no re-LOP) and skip the expensive
                 * speculative re-LOP+re-BPA whenever no cut could be accepted. On a
                 * clean scroll this is the common case (calibration 2026-06-05:
                 * all 15/100 4x5x5 candidates had parent_nondev < margin), and it
                 * is the lever that keeps DevCut near-free there while preserving
                 * the full gate on genuinely non-developable parents. */
                double pf = DevCut_nondev_fraction(arena, cm, dc_eps, NULL);
                if (pf < dc_margin) {
                    n_gate_rej++;
                    fprintf(stderr,
                        "    DevCut gate: comp %zu -> %zu pieces, parent_nondev="
                        "%.4f < margin %.3f -> reject (already developable; "
                        "re-LOP skipped)\n", i, n_pc, pf, dc_margin);
                } else {
                    /* re-LOP each piece from the parent's labelled ORIGINAL points
                     * (the margin-drop prevents vacuuming an adjacent close-wrap),
                     * which also yields aligned per-piece sub-clouds for Step 2/4. */
                    ComponentMesh    *rmesh = NULL; size_t n_rm = 0;
                    MeshResplitCloud *rcl   = NULL;
                    if (cl && cl->n > 0)
                        MeshResplit_remesh_pieces(arena, cl, pieces, n_pc, NULL, NULL,
                                                  &rmesh, &n_rm, &rcl);
                    if (!rmesh || n_rm < 2) { rmesh = pieces; n_rm = n_pc; rcl = NULL; }

                    /* lambda-energy gate: keep the split only if every re-LOP'd
                     * piece is more developable than the parent (and developable
                     * enough in absolute terms). */
                    double worst = 0.0;
                    for (size_t k = 0; k < n_rm; k++) {
                        double cf = DevCut_nondev_fraction(arena, &rmesh[k], dc_eps, NULL);
                        if (cf > worst) worst = cf;
                    }
                    if (worst + dc_margin <= pf && worst < dc_ceil) {
                        for (size_t k = 0; k < n_rm; k++)
                            DC_EMIT(rmesh[k], rcl ? rcl[k] : zerocl);
                        n_cut++;
                        accepted = 1;
                    } else {
                        n_gate_rej++;
                    }
                    fprintf(stderr,
                        "    DevCut gate: comp %zu -> %zu pieces, parent_nondev=%.4f "
                        "worst_piece=%.4f (margin %.3f ceil %.3f) -> %s\n",
                        i, n_rm, pf, worst, dc_margin, dc_ceil,
                        accepted ? "ACCEPT" : "reject");
                }
            }
            if (!accepted)
                DC_EMIT(*cm, cl ? *cl : zerocl);
        }
        #undef DC_EMIT
        #undef DC_ENSURE

        out->meshes = acc; out->n_meshes = cnt;
        rclouds = accc;
        fprintf(stderr,
            "  DevCut: %zu -> %zu components (%zu cut, %zu gate-rejected, "
            "%.3fs)\n", n_in, cnt, n_cut, n_gate_rej, elapsed_since(tdc));
        dump_stage(arena, in->dump_dir, in->cube_id, "step1c_devcut",
                   out->meshes, out->n_meshes);
    }

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
        size_t *range = (size_t *)ARENA_ALLOC(
            arena, (long)((total_in + 1) * sizeof(size_t)));
        size_t total_out = 0;
        size_t n_bridge = 0, n_ovl = 0, n_resurf = 0;

        /* Per-component split deadlines. The hardcoded 30s discards the lifted
         * multicut's RESULT on huge fused-wrap clumps: the deadline is only
         * checked AFTER the solve, so on an ~888k-face stacked-wrap blob the
         * multicut runs to completion (~71s, a valid 46-cluster cut) and then
         * gets thrown away by return_unchanged -- the whole stack flows on
         * UNCUT. Env-overridable so a monster cube can be re-run without a
         * rebuild; default behaviour unchanged. */
        double bridge_timeout  = env_d("BRIDGE_TIMEOUT_SEC",  30.0);
        double overlap_timeout = env_d("OVERLAP_TIMEOUT_SEC", 30.0);

        /* Phase A: topological split (oracle-gated bridge-cut + overlap) into
         * pieces; remember each input component's piece range for the step-4
         * re-LOP below. */
        for (size_t i = 0; i < total_in; i++) {
            ComponentMesh *cm = &out->meshes[i];
            range[i] = total_out;
            ComponentMesh *subs = NULL;
            size_t n_subs = 0;
            int brc = BridgeCut_process(arena, cm,
                                         ORACLE_UV_GRID_SIZE,
                                         bridge_timeout,
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
                                             in->n_threads, overlap_timeout,
                                             &osub, &n_osub);
                if (orc != 0 || n_osub == 0) { osub = &subs[j]; n_osub = 1; }
                if (n_osub > 1) n_ovl++;
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
                    total_out++;
                }
            }
        }
        range[total_in] = total_out;
        dump_stage(arena, in->dump_dir, in->cube_id, "step5_cc",
                   split_meshes, total_out);

        /* ---- Step 4: re-LOP + re-BPA each split piece from its OWN verts.
         * Each multicut/bridge-cut piece carries residual through-thickness
         * "lumps". Earlier this re-LOP'd from the parent's ORIGINAL cloud and
         * re-BPA'd -- but re-labelling the cloud onto a piece vacuums in the
         * ADJACENT close-wrap (~2-3 vox) and the re-BPA folds the now-thin sheet
         * (the step7_cc_bpa_003 self-intersection). Re-surfacing from the piece's
         * OWN verts instead (the split already assigned them correctly) cannot
         * vacuum/fold, and -- unlike topology-preserving in-place smoothing,
         * which squishes the kept faces into slivers QEM can't collapse -- it
         * re-triangulates cleanly. Only components that actually split. ---- */
        for (size_t i = 0; i < total_in; i++) {
            size_t start = range[i], n_cp = range[i + 1] - range[i];
            if (n_cp <= 1 || !rclouds || rclouds[i].n == 0) continue;
            for (size_t k = 0; k < n_cp; k++) {
                if (MeshResplit_resurface_own(arena, &split_meshes[start + k],
                                              rclouds[i].cell_origin,
                                              MLS_RESPLIT_ITERS) == 0)
                    n_resurf++;
            }
        }

        /* The bridge-cut / overlap vertex-splits can momentarily leave a bowtie
         * at a severed neck (a cut vertex left shared by both pieces). step0 is
         * already pinch-free, so split those here -- the split stage then emits a
         * vertex-manifold mesh by itself, keeping EVERY stage manifold rather
         * than relying on hole-fill's downstream pinch-split. Input is
         * edge-manifold (cuts only remove faces; re-BPA carries the Case-2
         * guard), so PinholeFill_split_pinches' precondition holds. */
        size_t split_pinch = 0;
        PinholeFill_split_pinches(arena, split_meshes, total_out, &split_pinch);
        out->meshes = split_meshes;
        out->n_meshes = total_out;
        double t_b = elapsed_since(tb);
        fprintf(stderr,
            "  Split: %zu -> %zu components (%zu bridge-cut, %zu overlap, "
            "%zu re-surfaced, %zu pinch-split, %.3fs)\n",
            total_in, total_out, n_bridge, n_ovl, n_resurf, split_pinch, t_b);

        dump_stage(arena, in->dump_dir, in->cube_id, "step7_cc_bpa",
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
        dump_stage(arena, in->dump_dir, in->cube_id, "step8_holefill",
                   out->meshes, out->n_meshes);
    }

    /* ---- Sever short handles (thin BPA self-bridges) ----
     * BPA can roll the ball into a thin self-connection, fusing a sheet into a
     * genus>0 tangle (~4/100 4x5x5 cubes, up to genus 10, all loops < 35 vox).
     * For each component, open every non-separating loop shorter than
     * SEVER_MAX_LOOP_VOX (manifold-preserving cut surgery via seam_cut); leave
     * longer loops. Runs AFTER hole fill so the opened slits are NOT re-closed,
     * and BEFORE QEM (genus is well defined on the clean filled mesh). These
     * handles are non-separating, so the Step-2 split (which severs SEPARATING
     * necks) cannot touch them -- this is the complementary cut. */
    if (!getenv("VES_SEVER_OFF")) {
        double ts = now_sec();
        size_t n_sev_comps = 0;
        long   total_sev = 0;
        for (size_t i = 0; i < out->n_meshes; i++) {
            ComponentMesh *cm = &out->meshes[i];
            if (cm->nf < 100) continue;
            float   *sv = NULL; int32_t *sf = NULL;
            size_t   snv = 0, snf = 0; long cut = 0;
            int srv_rc = SeamCut_sever_short_handles(arena, cm->verts, cm->nv,
                                                     cm->faces, cm->nf,
                                                     SEVER_MAX_LOOP_VOX,
                                                     &sv, &snv, &sf, &snf, &cut);
            if (srv_rc == 0 && cut > 0) {
                cm->verts = sv; cm->nv = snv;
                cm->faces = sf; cm->nf = snf;
                cm->vert_normals = NULL;   /* cut duplicated verts; normals stale */
                cm->nv_pre_fill = 0;
                cm->self = cm;
                n_sev_comps++;
                total_sev += cut;
            }
        }
        fprintf(stderr,
            "  Sever: %ld short handle(s) cut across %zu/%zu comp(s) (%.3fs)\n",
            total_sev, n_sev_comps, out->n_meshes, elapsed_since(ts));
        if (total_sev > 0)
            dump_stage(arena, in->dump_dir, in->cube_id, "step9_sever",
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
        dump_stage(arena, in->dump_dir, in->cube_id, "step10_qem",
                   out->meshes, out->n_meshes);
    }

    /* ---- Kibble removal: connectivity pass + surface-area filter ----
     * After hole-fill/QEM, split into connectivity-components and drop any whose
     * area is < KIBBLE_AREA_FRAC of the total meshed cube area (stray BPA
     * islands, cut-zone crumbs). Leaves only the real sheets. ---- */
    {
        double tk = now_sec();
        ComponentMesh *kept = NULL; size_t n_kept = 0, n_cull = 0;
        if (ComponentCull_by_area(arena, out->meshes, out->n_meshes,
                                  KIBBLE_AREA_FRAC, &kept, &n_kept, &n_cull) == 0
            && n_kept > 0) {
            out->meshes = kept;
            out->n_meshes = n_kept;
        }
        fprintf(stderr, "  Kibble: %zu culled, %zu kept (%.3fs)\n",
                n_cull, out->n_meshes, elapsed_since(tk));
        dump_stage(arena, in->dump_dir, in->cube_id, "step11_kibble",
                   out->meshes, out->n_meshes);
    }

    /* ---- Trim to owned region (halo mode only) ---- */
    if (in->halo_voxels > 0) {
        double tt = now_sec();
        /* Inset the owned box by BPA_OWNED_TRIM_INSET so the FINAL per-cube mesh
         * (step12_final = the weld input) never reaches a cube face — regardless
         * of which path produced it (extract, or a re-LOP'd split piece, which
         * re-BPAs WITHOUT the step0 cloud inset). This is the
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

    } else {
        /* No halo: trimmed == meshes by reference. */
        out->trimmed = out->meshes;
        out->n_trimmed = out->n_meshes;
        out->t_trim = 0.0;
    }

    /* ---- Final manifold guard (both paths): guarantee 2-manifold per-cube
     * output. The BPA Case-2 vertex guard prevents almost all pinches at the
     * source; this mops up trim-orphaned bowties and any residue -- resolve
     * >2-face edges, split pinch vertices, re-assert consistent winding. ---- */
    {
        double tg = now_sec();
        ManifoldGuardStats mg;
        ManifoldGuard_process(arena, out->trimmed, out->n_trimmed, 1 /*reorient*/, &mg);
        fprintf(stderr,
            "  Guard: %zu NM-edge(s) (-%zu faces), %zu pinch split(s), "
            "%zu reorient flip(s) (%.3fs)\n",
            mg.nm_edges_resolved, mg.faces_deleted, mg.pinch_splits,
            mg.orient_flips, elapsed_since(tg));
    }

    /* Dump the final per-cube mesh (halo mode -- this is the weld input).
     * dump_stage no-ops when dump_dir is NULL. */
    if (in->halo_voxels > 0) {
        double td = now_sec();
        dump_stage(arena, in->dump_dir, in->cube_id, "step12_final",
                   out->trimmed, out->n_trimmed);
        out->t_dump = elapsed_since(td);
    }

    /* ---- Developability gate: flag "bad sheets" in the final output ----
     * A finished sheet is bad if it has significant non-developable INTERIOR
     * vertices, or if the oracle still reports >1 sheet (a tangle the splitters
     * could not cut after one round). Detection + report only -- no cull. */
    {
        DevGateParams gp; DevGate_default_params(&gp);
        size_t n_bad = 0;
        for (size_t i = 0; i < out->n_trimmed; i++) {
            ComponentMesh *gcm = &out->trimmed[i];
            if (gcm->nf == 0) continue;
            DevGateVerdict gv;
            DevGate_classify(arena, gcm, &gp, &gv);
            if (gv.bad) {
                n_bad++;
                fprintf(stderr,
                    "  DevGate: comp %zu BAD [%s%s] int_bad=%zu/%zu (%.2f%%) "
                    "maxK=%.2f oracle=%d  (nv=%zu nf=%zu)\n",
                    i + 1,
                    (gv.reason & DEVGATE_NONDEV) ? "nondev " : "",
                    (gv.reason & DEVGATE_ORACLE) ? "oracle" : "",
                    gv.n_interior_bad, gv.n_interior,
                    100.0 * gv.frac_interior_bad, gv.max_abs_k,
                    gv.oracle_sheets, gcm->nv, gcm->nf);
            }
        }
        out->n_bad_sheets = n_bad;
        if (n_bad)
            fprintf(stderr, "  DevGate: %zu/%zu sheets flagged bad\n",
                    n_bad, out->n_trimmed);
    }

#ifdef VESUVIUS_DEBUG
    /* Invariant (CLAUDE.md §17): every emitted component is a 2-manifold. */
    for (size_t gi = 0; gi < out->n_trimmed; gi++) {
        ComponentMesh *gcm = &out->trimmed[gi];
        if (gcm->nf == 0) continue;
        MeshManifoldStats gms = MeshManifold_audit(arena, gcm->nv,
                                                   gcm->faces, gcm->nf);
        assert(MeshManifold_ok(&gms));
    }
#endif

    return 0;
}
