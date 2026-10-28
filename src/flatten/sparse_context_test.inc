#ifndef SPARSE_CONTEXT_TEST_IMPLEMENTATION_INCLUDED
#define SPARSE_CONTEXT_TEST_IMPLEMENTATION_INCLUDED
int Sparse_context_selftest(void)
{
    enum { N = 70, NT = 2 * N };
    int rows[NT], cols[NT], failures = 0;
    double values[NT], rhs[N], expected[N], answer[N], reference[N];
    SparseContext_T context = NULL;
    for (int trial = 0; trial < 4; trial++) {
        int count = 0;
        memset(rhs, 0, sizeof rhs);
        for (int i = 0; i < N; i++) {
            expected[i] = sin(.07 * i) + .3 * i;
            rows[count] = cols[count] = i;
            values[count++] = 4 + .01 * i + trial;
            if (i) {
                rows[count] = i;
                cols[count] = i - 1;
                values[count++] = -.7 - .1 * trial;
            }
        }
        if (trial >= 2) {
            rows[count] = N - 1;
            cols[count] = 0;
            values[count++] = -.2;
        }
        for (int k = 0; k < count; k++) {
            rhs[rows[k]] += values[k] * expected[cols[k]];
            if (rows[k] != cols[k]) {
                rhs[cols[k]] += values[k] * expected[rows[k]];
            }
        }
        if (Sparse_context_solve(&context, N, count, rows, cols, values, rhs, answer) ||
            Sparse_solve_sym(N, count, rows, cols, values, rhs, reference, SPARSE_SPD)) {
            failures++;
        }
        for (int i = 0; i < N; i++) {
            if (fabs(answer[i] - expected[i]) > 1e-10 ||
                fabs(answer[i] - reference[i]) > 1e-10) {
                failures++;
                break;
            }
        }
        for (int i = 0; i < N; i++) rhs[i] *= -3;
        if (Sparse_context_solve_factored(context, rhs, answer)) {
            failures++;
        }
        for (int i = 0; i < N; i++) {
            if (fabs(answer[i] + 3*expected[i]) > 1e-10) {
                failures++;
                break;
            }
        }
    }
    SparseContextStats stats;
    Sparse_context_stats(context, &stats);
    if (stats.calls != 4 || stats.pattern_builds != 2 ||
        stats.pattern_hits != 2 || stats.invalidations != 1 || stats.right_hand_sides != 8) {
        failures++;
    }
    double before[N];
    memcpy(before, answer, sizeof before);
    values[0] = -1;
    if (!Sparse_context_solve(&context, N, NT, rows, cols, values, rhs, answer) ||
        memcmp(before, answer, sizeof before) ||
        !Sparse_context_solve_factored(context, rhs, answer) ||
        memcmp(before, answer, sizeof before)) {
        failures++;
    }
    values[0] = 7;
    if (Sparse_context_solve(&context, N, NT, rows, cols, values, rhs, answer)) {
        failures++;
    }
    Sparse_context_free(&context);
    if (context) {
        failures++;
    }
    {
        int row[] = {0, 0}, col[] = {0, 0};
        double value[] = {2, 3}, b[] = {15}, x[] = {17};
        if (Sparse_context_solve(&context, 1, 2, row, col, value, b, x) || x[0] != 3) {
            failures++;
        }
        value[0] = NAN;
        if (!Sparse_context_solve(&context, 1, 2, row, col, value, b, x) || x[0] != 3) {
            failures++;
        }
        Sparse_context_free(&context);
    }
    fprintf(stderr, "  native owned symbolic context: changing values, pattern invalidation, failure recovery and scalar path: %s (%d failures)\n",
        failures ? "FAIL" : "ok", failures);
    return failures;
}
#endif
