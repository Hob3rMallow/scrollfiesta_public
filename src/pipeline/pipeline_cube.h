#ifndef PIPELINE_CUBE_INCLUDED
#define PIPELINE_CUBE_INCLUDED

#include <stddef.h>
#include <stdint.h>
#include "../common/arena.h"
#include "../common/mesh_types.h"

/* PipelineInput / PipelineOutput / pipeline_process_cube
 *
 * Per-cube pipeline as a library function. main.c is the thin driver
 * that parses args and calls this. The grid orchestrator (grid_pipeline)
 * will eventually spawn one subprocess per cube; each subprocess calls
 * this function exactly once.
 *
 * Pipeline stages (in order):
 *   1. Marching cubes + LOP projection + per-vertex backface cull
 *      (Step 0, via MeshExtract_run).
 *   2. CVT/RVD remeshing per connected chart
 *      (skipped only when `skip_simplify` is set for dense diagnostics).
 *   3. Mesh_trim_to_owned_box per component (when halo_voxels > 0) to
 *      yield the cube's owned region for grid_weld.
 *   4. Optional OBJ dumps, named by operation in pipeline order:
 *      step0_mls, step1_bpa, step2_cc, step3_mls, step4_bpa, step5_cc,
 *      step6_cc_mls, step7_cc_bpa, step8_holefill, step9_sever, step10_cvt,
 *      step11_kibble, step12_final. grid_weld consumes step12_final.
 */

typedef struct {
    const char *tiff_path;       /* input TIFF (per-cube) */
    const char *pred_dir;        /* parent dir of TIFF; halo neighbor lookup */
    const char *cube_id;         /* "zNNNNN_yNNNNN_xNNNNN" */
    int         halo_voxels;     /* 0 = no halo */
    int         cube_D, cube_H, cube_W;
    int         n_threads;
    float       cvt_target_ratio; /* 0 = pitch-aware default when axis is known */
    float       trim_inset;      /* owned-box inset, vox; < 0 = the
                                  * BPA_OWNED_TRIM_INSET default (1.0).
                                  * 0 = charts reach the cube faces: REQUIRED
                                  * for the incremental whole-grid unwrap
                                  * (no seam bridge to span a gap; the
                                  * mutual-nearest weld merges the near-
                                  * coincident boundary rows instead). Keep
                                  * the default for the grid_weld path -- its
                                  * seam bridge NEEDS the 2*inset gap.
                                  * Halo mode with TRIM_CLIP_AT_BOX (the
                                  * default since 2026-10-20) cuts the final
                                  * mesh at the EXACT cube box whatever this
                                  * is; it then only deepens the step-0
                                  * cloud's reach past
                                  * TRIM_CLIP_CLOUD_REACH_VOX into the halo. */
    const char *dump_dir;        /* if non-NULL, writes final meshes under
                                  * <dump_dir>/<cube_id>/<cube_id>_<stage>/ */
    int         skip_simplify;   /* dense diagnostic mode; production is CVT */
    int         surface_f0_inner; /* explicit F0/fixer inward-boundary geometry */
    int         full_dumps;     /* explicit --full-dumps opt-in: export intermediate
                                 * stages too. Zero-initialized callers write only
                                 * step12_final, the downstream mesh input. */
    int         cull_oracle_tangles; /* omit final components for which the
                                      * developability oracle still sees more
                                      * than one sheet; explicit fail-closed
                                      * recovery for an unsplittable tangle */

    /* Optional in-memory input. When vol_in != NULL the TIFF/halo loaders are
     * bypassed: vol_in is a padded (p_size_in)^3 uint8 buffer whose (0,0,0) is
     * world (cube_origin - halo_voxels), exactly like HaloLoader_load.
     * cube_origin_zyx is the owned cube origin (z,y,x). Used to stream a cube
     * directly from a remote zarr instead of reading per-cube TIFFs. */
    const uint8_t *vol_in;
    int            p_size_in;
    int64_t        cube_origin_zyx[3];
} PipelineInput;

typedef struct {
    ComponentMesh *meshes;       /* full mesh (halo-included when halo > 0) */
    size_t         n_meshes;
    ComponentMesh *trimmed;      /* owned-region trimmed mesh (halo > 0);
                                  * same as `meshes` when halo == 0 */
    size_t         n_trimmed;
    double         t_extract;    /* MC + LOP + per-vert cull */
    double         t_cvt;
    double         t_trim;
    double         t_dump;
    size_t         n_bad_sheets;  /* components flagged by the developability gate */
} PipelineOutput;

int pipeline_process_cube(Arena_T arena,
                          const PipelineInput *in,
                          PipelineOutput *out);

/* Selftests of the mesher's exactness guards (KD-tree capped queries, CVT
 * boundary BVH, CVT surface projection). Returns 0 when all pass. */
int pipeline_selftest(void);

#endif
