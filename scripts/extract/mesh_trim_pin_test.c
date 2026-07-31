/*
 * Mesh_trim_to_owned_box pin_mask propagation test.
 *
 * Build a mesh whose faces straddle the owned-box boundary at z=8 with
 * pinned verts in both the kept and dropped half. The trim should
 * keep only owned-side faces, remap kept verts, and emit a pin_mask
 * indexed by the new vert indices that flags exactly the pinned verts
 * that survived the trim.
 */
#include "common/ves_platform.h"
#include "common/arena.h"
#include "common/mesh_trim.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { \
    printf("  %-50s ", #name); \
    tests_run++; \
    if (test_##name()) { printf("PASS\n"); tests_passed++; } \
    else { printf("FAIL\n"); } \
} while(0)

static int find_vert_at(const float *verts, size_t nv,
                        float z, float y, float x, int32_t *out_idx)
{
    for (size_t v = 0; v < nv; v++) {
        if (fabsf(verts[v*3+0] - z) < 1e-5f &&
            fabsf(verts[v*3+1] - y) < 1e-5f &&
            fabsf(verts[v*3+2] - x) < 1e-5f) {
            *out_idx = (int32_t)v;
            return 1;
        }
    }
    return 0;
}

static int test_pin_mask_through_trim_preserves_owned(void)
{
    /* 6 verts in the (z,y,x) plane, straddling z=8 boundary.
     *   v0 z=2  PINNED, will survive (owned side)
     *   v1 z=4  free,   will survive (owned side, used by kept face)
     *   v2 z=6  free,   will survive (owned side, used by kept face)
     *   v3 z=10 PINNED, will NOT survive (halo side, not used by kept face)
     *   v4 z=12 free,   will NOT survive (halo side)
     *   v5 z=14 free,   will NOT survive (halo side) */
    float verts[] = {
        2.0f,  0.0f, 0.0f,
        4.0f,  1.0f, 0.0f,
        6.0f,  0.0f, 1.0f,
        10.0f, 0.0f, 0.0f,
        12.0f, 1.0f, 0.0f,
        14.0f, 0.0f, 1.0f,
    };
    /* Two faces, one in each region. Centroid is the avg z. */
    int32_t faces[] = {
        0, 1, 2,   /* centroid z = (2+4+6)/3 = 4   -> kept */
        3, 4, 5    /* centroid z = (10+12+14)/3=12 -> dropped */
    };
    uint8_t pin_in[6] = {1, 0, 0, 1, 0, 0};

    Arena_T arena = Arena_new();
    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;
    uint8_t *opin = NULL;

    int rc = Mesh_trim_to_owned_box(arena, verts, 6, faces, 2,
                                     pin_in,
                                     0.0f, 8.0f,
                                     &ov, &onv, &of, &onf,
                                     &opin);
    int ok = (rc == 0);
    if (ok) ok = (onf == 1);
    if (ok) ok = (onv == 3);
    if (ok) ok = (opin != NULL);

    /* The kept pinned vert (z=2) must be flagged. */
    if (ok) {
        int32_t idx = -1;
        if (!find_vert_at(ov, onv, 2.0f, 0.0f, 0.0f, &idx)) {
            printf("(kept pin missing) ");
            ok = 0;
        } else if (!opin[idx]) {
            printf("(kept pin not flagged) ");
            ok = 0;
        }
    }

    /* No surviving vert should be unexpectedly pinned (v1, v2 free
     * verts must not have inherited pins from elsewhere). */
    if (ok) {
        int32_t idx = -1;
        if (find_vert_at(ov, onv, 4.0f, 1.0f, 0.0f, &idx)) {
            if (opin[idx]) {
                printf("(free vert incorrectly pinned) ");
                ok = 0;
            }
        }
    }

    Arena_dispose(&arena);
    return ok;
}

static int test_pin_mask_null_in_opt_out(void)
{
    /* NULL pin_mask_in OR NULL out_pin_mask should be a no-op (silently
     * not emit pins). */
    float verts[] = {0,0,0, 1,0,0, 0,1,0};
    int32_t faces[] = {0,1,2};

    Arena_T arena = Arena_new();
    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;
    uint8_t *opin = NULL;

    int rc = Mesh_trim_to_owned_box(arena, verts, 3, faces, 1,
                                     NULL,           /* no pins in */
                                     -10.0f, 10.0f,
                                     &ov, &onv, &of, &onf,
                                     &opin);
    int ok = (rc == 0 && onf == 1 && onv == 3 && opin == NULL);
    Arena_dispose(&arena);
    return ok;
}

int mesh_trim_pin_test_main(void)
{
    printf("\nMesh_trim_to_owned_box pin_mask propagation tests\n");
    tests_run = 0;
    tests_passed = 0;

    TEST(pin_mask_through_trim_preserves_owned);
    TEST(pin_mask_null_in_opt_out);

    printf("  %d/%d passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}

#ifndef TEST_HARNESS
int main(void)
{
    return mesh_trim_pin_test_main();
}
#endif
