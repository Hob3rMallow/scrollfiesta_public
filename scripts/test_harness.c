/*
 * test_harness.c — Unified test runner.
 *
 * Calls each test's entry point, measures elapsed time,
 * tallies pass/fail, and prints a summary report.
 *
 * Build with: -DTEST_HARNESS
 */
#include "common/ves_platform.h"

#include <stdio.h>
#include <string.h>

#include "test_harness.h"

typedef int (*test_fn)(void);

typedef struct {
    const char *name;
    test_fn     func;
} TestEntry;

static const TestEntry tests[] = {
    { "test_common",                test_common_main                },
    { "qem_test",                   qem_test_main                   },
    { "qem_pin_test",               qem_pin_test_main               },
    { "mesh_trim_pin_test",         mesh_trim_pin_test_main         },
    { "pipeline_cube_smoke_test",   pipeline_cube_smoke_test_main   },
    { "mls_iter_test",              mls_iter_test_main              },
    { "manifold_guard_test",        manifold_guard_test_main        },
    { "pred_reject_test",           pred_reject_test_main           },
    { "cc_color_test",              cc_color_test_main              },
};

#define N_TESTS (sizeof(tests) / sizeof(tests[0]))

int main(void)
{
    fprintf(stderr, "\n=== Vesuvius-C Test Suite ===\n");

    int n_pass = 0;
    int n_fail = 0;
    const char *failed_names[N_TESTS];
    for (size_t j = 0; j < N_TESTS; j++) failed_names[j] = NULL;

    for (size_t i = 0; i < N_TESTS; i++) {
        double t0 = ves_clock_sec();
        int rc = tests[i].func();
        double elapsed = ves_clock_sec() - t0;

        if (rc == 0) {
            fprintf(stderr, "[ PASS ] %-28s (%.3fs)\n",
                    tests[i].name, elapsed);
            n_pass++;
        } else {
            fprintf(stderr, "[ FAIL ] %-28s (%.3fs)  rc=%d\n",
                    tests[i].name, elapsed, rc);
            failed_names[n_fail] = tests[i].name;
            n_fail++;
        }
    }

    fprintf(stderr, "\n=== %d/%d passed, %d failed ===\n",
            n_pass, (int)N_TESTS, n_fail);

    if (n_fail > 0) {
        fprintf(stderr, "Failed:");
        for (int i = 0; i < n_fail; i++) {
            fprintf(stderr, " %s", failed_names[i]);
        }
        fprintf(stderr, "\n");
    }

    return (n_fail == 0) ? 0 : 1;
}
