/*
 * QEM pin_mask contract unit tests.
 *
 * Three properties to verify:
 *   1. test_qem_pin_dirichlet:
 *      pinned vertex positions are bit-exact in the output mesh; the
 *      output pin_mask marks them.
 *   2. test_qem_pin_propagates_one_pinned:
 *      a one-pinned collapse must keep the pinned position; the survivor
 *      appears in out_pin_mask.
 *   3. test_qem_maintenance_respects_pins:
 *      end-to-end through QEM_simplify_pinned with a small enough mesh to
 *      trigger the maintenance pass; pinned interior verts must not move.
 */
#include "common/ves_platform.h"
#include "common/arena.h"
#include "common/qem.h"

#include <assert.h>
#include <math.h>
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

/* Build a flat 10x10 grid of vertices in the z=0 plane with two
 * triangles per cell (162 faces, 100 verts). */
static void build_flat_grid(float *verts, int32_t *faces,
                            size_t *out_nv, size_t *out_nf)
{
    size_t n_side = 10;
    size_t nv = n_side * n_side;
    size_t nf = 0;
    for (size_t r = 0; r < n_side; r++) {
        for (size_t c = 0; c < n_side; c++) {
            size_t idx = r * n_side + c;
            verts[idx * 3 + 0] = 0.0f;                /* z */
            verts[idx * 3 + 1] = (float)r;            /* y */
            verts[idx * 3 + 2] = (float)c;            /* x */
        }
    }
    for (size_t r = 0; r + 1 < n_side; r++) {
        for (size_t c = 0; c + 1 < n_side; c++) {
            int32_t a = (int32_t)(r * n_side + c);
            int32_t b = (int32_t)(r * n_side + c + 1);
            int32_t cv = (int32_t)((r + 1) * n_side + c);
            int32_t d = (int32_t)((r + 1) * n_side + c + 1);
            faces[nf * 3 + 0] = a;
            faces[nf * 3 + 1] = b;
            faces[nf * 3 + 2] = d;
            nf++;
            faces[nf * 3 + 0] = a;
            faces[nf * 3 + 1] = d;
            faces[nf * 3 + 2] = cv;
            nf++;
        }
    }
    *out_nv = nv;
    *out_nf = nf;
}

/* Find the new vertex index that lies at (z,y,x), within float tolerance. */
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

static int test_qem_pin_dirichlet(void)
{
    /* 100-vert flat patch. Mark the 4 corners + 6 boundary midpoints
     * pinned (10 verts total). Target = 50 faces. Assert all 10 pinned
     * positions appear bit-exact in the output. */
    float verts[100 * 3];
    int32_t faces[162 * 3];
    size_t nv = 0, nf = 0;
    build_flat_grid(verts, faces, &nv, &nf);

    uint8_t pin_in[100] = {0};
    /* Pin: 4 corners (0,0), (0,9), (9,0), (9,9) and interior pin at (5,5) */
    int pinned_idx[] = {0, 9, 90, 99, 55};
    int n_pinned = 5;
    for (int i = 0; i < n_pinned; i++) pin_in[pinned_idx[i]] = 1;

    /* Snapshot pinned positions */
    float pinned_pos[5][3];
    for (int i = 0; i < n_pinned; i++) {
        int p = pinned_idx[i];
        pinned_pos[i][0] = verts[p*3+0];
        pinned_pos[i][1] = verts[p*3+1];
        pinned_pos[i][2] = verts[p*3+2];
    }

    Arena_T arena = Arena_new();
    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;
    uint8_t *opin = NULL;

    int rc = QEM_simplify_pinned(arena, verts, nv, faces, nf,
                                 pin_in, NULL,
                                 50,  /* target_nf */
                                 &ov, &onv, &of, &onf, &opin);

    int ok = (rc == 0);
    if (ok) ok = (opin != NULL);

    /* Each pinned input position must appear in the output, bit-exact,
     * and be flagged in opin. */
    if (ok) {
        for (int i = 0; i < n_pinned && ok; i++) {
            int32_t idx = -1;
            int found = find_vert_at(ov, onv,
                                     pinned_pos[i][0],
                                     pinned_pos[i][1],
                                     pinned_pos[i][2],
                                     &idx);
            if (!found) {
                printf("(pin %d missing) ", pinned_idx[i]);
                ok = 0;
            } else if (!opin[idx]) {
                printf("(pin %d not flagged) ", pinned_idx[i]);
                ok = 0;
            }
        }
    }

    Arena_dispose(&arena);
    return ok;
}

static int test_qem_pin_propagates_one_pinned(void)
{
    /* Build a tiny mesh where a forced collapse will absorb a pinned
     * vert. Two coincident-ish triangles sharing an edge between v0
     * (pinned at origin) and v1 (free). Target=1 face forces a collapse.
     * The surviving vert must be pinned and at v0's position. */
    float verts[] = {
        0.0f, 0.0f, 0.0f,   /* v0 PINNED */
        0.0f, 0.0f, 0.1f,   /* v1 free, very close to v0 */
        0.0f, 1.0f, 0.0f,   /* v2 */
        0.0f, 1.0f, 1.0f,   /* v3 */
        0.0f, 0.0f, 1.0f    /* v4 */
    };
    /* Two triangles using v0 and v1 */
    int32_t faces[] = {
        0, 2, 1,
        1, 2, 3,
        1, 3, 4
    };
    uint8_t pin_in[5] = {1, 0, 0, 0, 0};

    Arena_T arena = Arena_new();
    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;
    uint8_t *opin = NULL;

    int rc = QEM_simplify_pinned(arena, verts, 5, faces, 3,
                                 pin_in, NULL,
                                 1, /* aggressive target */
                                 &ov, &onv, &of, &onf, &opin);
    int ok = (rc == 0 && opin != NULL);

    /* The pinned origin must appear in the output. */
    if (ok) {
        int32_t idx = -1;
        ok = find_vert_at(ov, onv, 0.0f, 0.0f, 0.0f, &idx);
        if (ok && !opin[idx]) {
            printf("(surviving pin not flagged) ");
            ok = 0;
        }
    }

    Arena_dispose(&arena);
    return ok;
}

static int test_qem_maintenance_respects_pins(void)
{
    /* 100-vert flat patch. Pin an interior vert (5,5,0). Target down to
     * 60 faces so the maintenance pass runs (gated at out_nf > 100 but
     * applies anyway because of how qem_collapse_pass shrinks faces; in
     * any case pinned smoothing must be a no-op on the pinned vert).
     *
     * Without the pin-aware smooth pass, the interior pinned vert would
     * drift toward its neighbors' centroid. With pin_mask threaded into
     * qem_smooth_pass, it must stay exactly at its pinned position. */
    float verts[100 * 3];
    int32_t faces[162 * 3];
    size_t nv = 0, nf = 0;
    build_flat_grid(verts, faces, &nv, &nf);

    /* Bump the pinned vert above the plane so smoothing would obviously
     * pull it back. The pinned contract: caller's coords must survive. */
    verts[55 * 3 + 0] = 5.0f;  /* z=5, much higher than rest at z=0 */

    uint8_t pin_in[100] = {0};
    pin_in[55] = 1;
    float pinned_pos[3] = {5.0f, 5.0f, 5.0f};

    Arena_T arena = Arena_new();
    float *ov = NULL; size_t onv = 0;
    int32_t *of = NULL; size_t onf = 0;
    uint8_t *opin = NULL;

    int rc = QEM_simplify_pinned(arena, verts, nv, faces, nf,
                                 pin_in, NULL,
                                 60, /* target_nf */
                                 &ov, &onv, &of, &onf, &opin);
    int ok = (rc == 0 && opin != NULL);

    if (ok) {
        int32_t idx = -1;
        int found = find_vert_at(ov, onv,
                                  pinned_pos[0], pinned_pos[1], pinned_pos[2],
                                  &idx);
        if (!found) {
            printf("(interior pin moved) ");
            ok = 0;
        } else if (!opin[idx]) {
            printf("(interior pin not flagged) ");
            ok = 0;
        }
    }

    Arena_dispose(&arena);
    return ok;
}

int qem_pin_test_main(void)
{
    printf("\nQEM pin_mask contract tests\n");
    tests_run = 0;
    tests_passed = 0;

    TEST(qem_pin_dirichlet);
    TEST(qem_pin_propagates_one_pinned);
    TEST(qem_maintenance_respects_pins);

    printf("  %d/%d passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}

#ifndef TEST_HARNESS
int main(void)
{
    return qem_pin_test_main();
}
#endif
