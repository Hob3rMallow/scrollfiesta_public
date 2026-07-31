/*
 * pred_reject_test.c -- garbage-prediction-cube detector unit tests.
 *
 * Delegates to PredReject_selftest, which builds synthetic volumes with known
 * structure (solid slab -> reject; thin sheets / flat face-on patch / empty ->
 * keep; slab + protrusion -> reject) and asserts the verdicts. Wired into the
 * all_tests harness so it runs with every Windows build.
 */
#include "common/ves_platform.h"
#include "extract/pred_reject.h"

#ifdef TEST_HARNESS
#define main pred_reject_test_main
#endif

int main(void)
{
    return PredReject_selftest();
}
