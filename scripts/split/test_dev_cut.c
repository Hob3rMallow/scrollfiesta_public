/*
 * test_dev_cut.c -- unit-test driver for the developability splitter
 * (src/split/dev_cut.c). Runs the splitter's own selftest plus its two direct
 * dependencies: the shared 3x3 eigensolver and the Crane developability energy
 * (whose new per-vertex API the splitter thresholds).
 *
 * Build: dev_cut_test.vcxproj (Windows) / `make dev_cut_test` (Linux).
 */
#include "split/dev_cut.h"
#include "flatten/developability.h"
#include "common/eig3.h"

#include <stdio.h>

int main(void)
{
    int fails = 0;
    printf("======================================\n");
    printf(" dev_cut unit tests\n");
    printf("======================================\n");

    fails += Eig3_selftest();       /* dependency: symmetric 3x3 eigensolver   */
    fails += Develop_selftest();    /* dependency: Crane energy + gradient path */
    fails += DevCut_selftest();     /* the developability splitter itself       */

    printf("======================================\n");
    printf(" dev_cut: %d failure(s)\n", fails);
    printf("======================================\n");
    return fails ? 1 : 0;
}
