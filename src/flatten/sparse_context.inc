#ifndef SPARSE_CONTEXT_IMPLEMENTATION_INCLUDED
#define SPARSE_CONTEXT_IMPLEMENTATION_INCLUDED
/* Included by sparse_solve.c so all TAUCS ordering/factor/solve ownership
 * remains in the established adapter. No shared pattern-cache eviction. */
#include "sparse_context.h"
#include "../common/ves_platform.h"

struct SparseContext {
    SpPattern pattern;
    double *pb, *px, *solution;
    double scalar;
    int n, numeric_valid;
    SparseContextStats stats;
};

void Sparse_context_stats(SparseContext_T context, SparseContextStats *stats)
{
    if (stats) {
        if (context) {
            *stats = context->stats;
        } else {
            memset(stats, 0, sizeof *stats);
        }
    }
}

void Sparse_context_free(SparseContext_T *context)
{
    if (!context || !*context) {
        return;
    }
    SparseContext_T owned = *context;
    sp_pattern_clear(&owned->pattern);
    free(owned->pb);
    free(owned->px);
    free(owned->solution);
    free(owned);
    *context = NULL;
}

static int sp_context_buffers(SparseContext_T context, int n)
{
    if (context->n == n && context->pb && context->px && context->solution) {
        return 1;
    }
    double *pb = malloc((size_t)n * sizeof *pb);
    double *px = malloc((size_t)n * sizeof *px);
    double *solution = malloc((size_t)n * sizeof *solution);
    if (!pb || !px || !solution) {
        free(pb);
        free(px);
        free(solution);
        return 0;
    }
    free(context->pb);
    free(context->px);
    free(context->solution);
    context->pb = pb;
    context->px = px;
    context->solution = solution;
    context->n = n;
    return 1;
}

static int sp_context_pattern(SparseContext_T context, taucs_ccs_matrix *matrix)
{
    SpPattern *pattern = &context->pattern;
    if (sp_pattern_matches(pattern, matrix)) {
        context->stats.pattern_hits++;
        sp_factor_free(pattern->L, 1);
        return 1;
    }
    context->stats.invalidations += pattern->L != NULL;
    sp_pattern_clear(pattern);
    int n = matrix->n, nnz = matrix->colptr[n];
    pattern->colptr = malloc((size_t)(n + 1) * sizeof *pattern->colptr);
    pattern->rowind = malloc((size_t)(nnz ? nnz : 1) * sizeof *pattern->rowind);
    if (!pattern->colptr || !pattern->rowind) {
        sp_pattern_clear(pattern);
        return 0;
    }
    double started = ves_clock_sec();
    int rc = sp_symbolic(matrix, &pattern->perm, &pattern->invperm, &pattern->L);
    context->stats.symbolic_seconds += ves_clock_sec() - started;
    if (rc) {
        sp_pattern_clear(pattern);
        return 0;
    }
    memcpy(pattern->colptr, matrix->colptr, (size_t)(n + 1) * sizeof *pattern->colptr);
    memcpy(pattern->rowind, matrix->rowind, (size_t)nnz * sizeof *pattern->rowind);
    pattern->n = n;
    pattern->nnz = nnz;
    context->stats.pattern_builds++;
    return 1;
}

int Sparse_context_factor(SparseContext_T *address, int n, int nt,
    const int *rows, const int *cols, const double *values)
{
    if (address && *address) {
        (*address)->numeric_valid = 0;
    }
    if (!address || n <= 0 || nt < 0 || !rows || !cols || !values) {
        return -1;
    }
    for (int k = 0; k < nt; k++) {
        if (rows[k] < 0 || rows[k] >= n || cols[k] < 0 ||
            cols[k] > rows[k] || !isfinite(values[k])) {
            return -1;
        }
    }
    if (!*address) {
        *address = calloc(1, sizeof **address);
    }
    SparseContext_T context = *address;
    if (!context || !sp_context_buffers(context, n)) {
        return -1;
    }
    context->stats.calls++;
    if (n == 1) {
        double diagonal = 0;
        for (int k = 0; k < nt; k++) {
            diagonal += values[k];
        }
        if (!(diagonal > 0) || !isfinite(diagonal)) {
            return -1;
        }
        if (context->pattern.L) {
            context->stats.invalidations++;
            sp_pattern_clear(&context->pattern);
        }
        context->scalar = diagonal;
        context->numeric_valid = 1;
        return 0;
    }
    double started = ves_clock_sec();
    taucs_ccs_matrix *matrix = sp_build_ccs(n, nt, rows, cols, values);
    context->stats.assembly_seconds += ves_clock_sec() - started;
    if (!matrix) {
        return -1;
    }
    int ok = sp_context_pattern(context, matrix);
    if (ok) {
        started = ves_clock_sec();
        ok = !sp_numeric(matrix, context->pattern.perm,
            context->pattern.invperm, context->pattern.L);
        context->stats.numeric_seconds += ves_clock_sec() - started;
    }
    sp_ccs_free(matrix);
    if (!ok) {
        return -1;
    }
    context->numeric_valid = 1;
    return 0;
}

int Sparse_context_solve_factored(SparseContext_T context,
    const double *rhs, double *result)
{
    if (!context || !context->numeric_valid || !rhs || !result) {
        return -1;
    }
    int n = context->n;
    for (int i = 0; i < n; i++) {
        if (!isfinite(rhs[i])) {
            return -1;
        }
    }
    context->stats.right_hand_sides++;
    if (n == 1) {
        double solution = rhs[0] / context->scalar;
        if (!isfinite(solution)) {
            return -1;
        }
        result[0] = solution;
        return 0;
    }
    double started = ves_clock_sec();
    int rc = 0;
#ifndef SPARSE_TAUCS_THREAD_SAFE
#pragma omp critical(sparse_taucs_call)
#endif
    {
        taucs_vec_permute(n, TAUCS_DOUBLE, (void *)rhs, context->pb, context->pattern.perm);
        rc = taucs_supernodal_solve_llt(context->pattern.L, context->px, context->pb);
        if (!rc) {
            taucs_vec_ipermute(n, TAUCS_DOUBLE, context->px,
                context->solution, context->pattern.perm);
        }
    }
    context->stats.solve_seconds += ves_clock_sec() - started;
    if (rc) {
        return -1;
    }
    for (int i = 0; i < n; i++) {
        if (!isfinite(context->solution[i])) {
            return -1;
        }
    }
    memcpy(result, context->solution, (size_t)n * sizeof *result);
    return 0;
}

int Sparse_context_solve(SparseContext_T *address, int n, int nt,
    const int *rows, const int *cols, const double *values,
    const double *rhs, double *result)
{
    if (!rhs || !result ||
        Sparse_context_factor(address, n, nt, rows, cols, values)) {
        return -1;
    }
    return Sparse_context_solve_factored(*address, rhs, result);
}
#endif
