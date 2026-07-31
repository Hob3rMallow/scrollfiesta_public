/*
 * test_harness.h — Unified test harness declarations.
 *
 * Each test file's main() is conditionally renamed to xxx_main()
 * when compiled with -DTEST_HARNESS.  This header declares those
 * entry points so test_harness.c can call them all.
 *
 * After the May-2026 cleanup, tests are scoped to the active code path:
 *   Step 0 marching cubes + LOP + per-vert cull + QEM(pinned) + weld.
 * Tests for quarantined steps live alongside their step in attic/.
 */
#ifndef TEST_HARNESS_H
#define TEST_HARNESS_H

int test_common_main(void);
int qem_test_main(void);
int qem_pin_test_main(void);
int mesh_trim_pin_test_main(void);
int pipeline_cube_smoke_test_main(void);
int mls_iter_test_main(void);
int manifold_guard_test_main(void);
int pred_reject_test_main(void);
int cc_color_test_main(void);

#endif /* TEST_HARNESS_H */
