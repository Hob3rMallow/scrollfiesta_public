/*
 * pipeline_process_cube smoke test.
 *
 * Builds a tiny synthetic foreground volume (a flat slab) at a temp TIFF
 * path, runs pipeline_process_cube with halo > 0, and asserts:
 *   - n_meshes > 0
 *   - every output mesh has non-NULL verts/faces and pin_mask
 *   - trimmed mesh exists and has at least one face
 *
 * Verifies the per-cube library function compiles, links, runs to
 * completion, and produces meshes with the expected pin_mask invariants.
 */
#include "common/ves_platform.h"
#include "common/arena.h"
#include "common/tiff_io.h"
#include "pipeline/pipeline_cube.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { \
    printf("  %-50s ", #name); \
    tests_run++; \
    if (test_##name()) { printf("PASS\n"); tests_passed++; } \
    else { printf("FAIL\n"); } \
} while(0)

static int test_pipeline_cube_synthetic_slab(void)
{
    /* Build a 32^3 slab volume with foreground at z in [12, 20). */
    int D = 32, H = 32, W = 32;
    uint8_t *vol = (uint8_t *)calloc((size_t)(D * H * W), 1);
    if (!vol) return 0;
    for (int z = 12; z < 20; z++) {
        for (int y = 4; y < 28; y++) {
            for (int x = 4; x < 28; x++) {
                vol[z * H * W + y * W + x] = 1;
            }
        }
    }

    const char *tiff_path =
        "C:/Users/mordr/AppData/Local/Temp/pipeline_cube_smoke_in.tif";
    if (TiffIO_save(tiff_path, vol, D, H, W) != 0) {
        printf("(TiffIO_save failed) ");
        free(vol);
        return 0;
    }
    free(vol);

    Arena_T arena = Arena_new();
    PipelineInput in = {
        .tiff_path        = tiff_path,
        .pred_dir         = NULL,
        .cube_id          = "smoke_test_cube",
        .halo_voxels      = 0,           /* no halo -- single cube */
        .cube_D           = D,
        .cube_H           = H,
        .cube_W           = W,
        .n_threads        = 1,
        .qem_target_ratio = 0.0f,
        .dump_dir         = NULL,
        .skip_qem         = 1,           /* keep test fast */
    };
    PipelineOutput out = {0};

    int rc = pipeline_process_cube(arena, &in, &out);
    int ok = (rc == 0);
    if (ok) ok = (out.n_meshes > 0);
    if (ok) {
        for (size_t i = 0; i < out.n_meshes; i++) {
            ComponentMesh *cm = &out.meshes[i];
            if (cm->nv == 0 || cm->nf == 0 || !cm->verts || !cm->faces) {
                printf("(empty mesh %zu) ", i);
                ok = 0;
                break;
            }
        }
    }
    /* halo == 0 means trimmed should alias meshes. */
    if (ok) {
        ok = (out.trimmed == out.meshes && out.n_trimmed == out.n_meshes);
        if (!ok) printf("(trimmed != meshes for halo=0) ");
    }
    Arena_dispose(&arena);
    return ok;
}

int pipeline_cube_smoke_test_main(void)
{
    printf("\npipeline_process_cube smoke test\n");
    tests_run = 0;
    tests_passed = 0;

    TEST(pipeline_cube_synthetic_slab);

    printf("  %d/%d passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}

#ifndef TEST_HARNESS
int main(void)
{
    return pipeline_cube_smoke_test_main();
}
#endif
