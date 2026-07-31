/*
 * test_depth_peel.c -- unit-test driver for the stacked-wrap separator
 * (src/split/depth_peel.c). Runs the splitter's own selftest (flat / curved
 * self-gate, two-sheet wall-bridge split, trivial input).
 *
 * Build: depth_peel_test.vcxproj (Windows) / `make depth_peel_test` (Linux).
 */
#include "split/depth_peel.h"

#include <stdio.h>

int main(void)
{
    int fails = 0;
    printf("======================================\n");
    printf(" depth_peel unit tests\n");
    printf("======================================\n");

    fails += DepthPeel_selftest();  /* the stacked-wrap separator itself */

    printf("======================================\n");
    printf(" depth_peel: %d failure(s)\n", fails);
    printf("======================================\n");
    return fails ? 1 : 0;
}
