/*
 * cc_color_test.c -- unit test for src/common/cc_color (connected-component
 * vertex coloring, baked into the pipeline stage dumps + grid_weld welded.obj).
 * Thin wrapper over CCColor_selftest: two disjoint triangles plus a triangle
 * sharing one vertex with the first (a bowtie) -> exactly 2 components, the
 * bowtie-joined verts share a colour, the disjoint triangle differs, and empty
 * faces grey-fill without crashing.
 *
 * Build with -DTEST_HARNESS to expose cc_color_test_main() to test_harness.c.
 */
#include "common/cc_color.h"

#include <stdio.h>

#ifdef TEST_HARNESS
int cc_color_test_main(void)
#else
int main(void)
#endif
{
    int fails = CCColor_selftest();
    if (fails != 0) fprintf(stderr, "cc_color_test: FAIL (%d failing check(s))\n", fails);
    return fails;
}
