/*
 * oracle_unit_test.c — Unit tests for Step 1: Oracle sheet counting.
 *
 * Loads OBJ meshes from test_sheets/ corpus and verifies Oracle_count_sheets
 * returns expected sheet counts.  No external data dependencies beyond the
 * committed test_sheets/ directory.
 *
 * Usage (standalone): ./oracle_unit_test
 * Usage (harness):    compiled with -DTEST_HARNESS, called as oracle_unit_test_main()
 */
#include "common/ves_platform.h"

#include "split/oracle.h"
#include "common/arena.h"
#include "common/obj_io.h"
#include "common/pca.h"
#include "common/pipeline_constants.h"
#include "common/mesh_types.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ================================================================
 * Test infrastructure
 * ================================================================ */

static int g_pass = 0;
static int g_fail = 0;

#define TEST_SHEETS_DIR "test_sheets"

static int run_oracle_on_obj(const char *obj_path, int grid_size)
{
    Arena_T arena = Arena_new();

    float *verts = NULL;
    int32_t *faces = NULL;
    size_t nv = 0, nf = 0;

    if (ObjIO_read(arena, obj_path, &verts, &nv, &faces, &nf) != 0) {
        fprintf(stderr, "  ERROR: cannot read %s\n", obj_path);
        Arena_dispose(&arena);
        return -1;
    }

    ComponentMesh cm;
    memset(&cm, 0, sizeof(cm));
    cm.verts = verts;
    cm.faces = faces;
    cm.nv = nv;
    cm.nf = nf;
    cm.self = &cm;
    PCA_normal(verts, nv, cm.pca_normal, cm.centroid);

    int sheets = Oracle_count_sheets(arena, &cm, grid_size);

    Arena_dispose(&arena);
    return sheets;
}

static void test_oracle(const char *name, const char *obj_path,
                         int min_expected, int max_expected)
{
    int sheets = run_oracle_on_obj(obj_path, ORACLE_UV_GRID_SIZE);
    if (sheets < 0) {
        fprintf(stderr, "  %-30s  FAIL (load error)\n", name);
        g_fail++;
        return;
    }

    int ok = (sheets >= min_expected && sheets <= max_expected);
    if (ok) {
        fprintf(stderr, "  %-30s  PASS (sheets=%d, expected %d..%d)\n",
                name, sheets, min_expected, max_expected);
        g_pass++;
    } else {
        fprintf(stderr, "  %-30s  FAIL (sheets=%d, expected %d..%d)\n",
                name, sheets, min_expected, max_expected);
        g_fail++;
    }
}

/* ================================================================
 * Main
 * ================================================================ */

#ifdef TEST_HARNESS
int oracle_unit_test_main(void)
#else
int main(void)
#endif
{
    fprintf(stderr, "Oracle Unit Tests\n");
    fprintf(stderr, "========================================\n");

    /* overlap_split_1: renamed from spectral_split_1, has multiple overlapping sheets */
    test_oracle("overlap_split_1",
                TEST_SHEETS_DIR "/overlap_split_1.obj",
                2, 99);

    /* overlap_split_2: should produce at least 3 sheets */
    test_oracle("overlap_split_2",
                TEST_SHEETS_DIR "/overlap_split_2.obj",
                3, 99);

    fprintf(stderr, "========================================\n");
    fprintf(stderr, "Results: %d/%d passed\n", g_pass, g_pass + g_fail);

    return (g_fail == 0) ? 0 : 1;
}
