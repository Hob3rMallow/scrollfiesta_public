#include "monotone_qp.h"

#include "../flatten/sparse_solve.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct {
    int32_t v[2];
    double  c[2];
    int     n;
    int32_t bound;
} EqualityRow;

static int32_t mq_find(int32_t *parent, int32_t x)
{
    int32_t r = x;
    while (parent[r] != r) r = parent[r];
    while (parent[x] != x) {
        int32_t next = parent[x];
        parent[x] = r;
        x = next;
    }
    return r;
}

/* Returns 1 when the two sets were distinct and have now been merged. */
static int mq_union(int32_t *parent, uint8_t *rank, int32_t a, int32_t b)
{
    a = mq_find(parent, a);
    b = mq_find(parent, b);
    if (a == b) return 0;
    if (rank[a] < rank[b]) {
        parent[a] = b;
    } else {
        parent[b] = a;
        if (rank[a] == rank[b]) rank[a]++;
    }
    return 1;
}

void MonotoneQpOptions_default(MonotoneQpOptions *opts)
{
    if (opts == NULL) return;
    opts->max_active_iterations = 512;
    opts->feasibility_tolerance = 1e-10;
    opts->direction_tolerance = 1e-11;
    opts->multiplier_tolerance = 1e-10;
    opts->linear_residual_tolerance = 1e-10;
    opts->verbose = 0;
}

static int validate_problem(const MonotoneQpProblem *p)
{
    if (p == NULL || p->nvar == 0 || p->nvar > (size_t)INT_MAX) return -1;
    if ((p->nrows > 0 && (p->rows == NULL || p->coeff == NULL)) ||
        (p->nbounds > 0 && p->bounds == NULL) ||
        (p->nanchors > 0 && p->anchors == NULL)) return -1;

    for (size_t r = 0; r < p->nrows; r++) {
        const MonotoneQpRow *row = &p->rows[r];
        if (row->count <= 0 || row->first > p->ncoeff ||
            (size_t)row->count > p->ncoeff - row->first ||
            !isfinite(row->target) || !isfinite(row->weight) ||
            row->weight < 0.0) return -1;
        for (int32_t j = 0; j < row->count; j++) {
            const MonotoneQpCoeff *c = &p->coeff[row->first + (size_t)j];
            if (c->var < 0 || (size_t)c->var >= p->nvar ||
                !isfinite(c->value)) return -1;
        }
    }
    for (size_t e = 0; e < p->nbounds; e++) {
        const MonotoneQpBound *b = &p->bounds[e];
        if (b->lo < 0 || b->hi < 0 ||
            (size_t)b->lo >= p->nvar || (size_t)b->hi >= p->nvar ||
            b->lo == b->hi || !isfinite(b->lower)) return -1;
    }
    for (size_t a = 0; a < p->nanchors; a++) {
        const MonotoneQpAnchor *an = &p->anchors[a];
        if (an->var < 0 || (size_t)an->var >= p->nvar ||
            !isfinite(an->value)) return -1;
        for (size_t b = 0; b < a; b++)
            if (p->anchors[b].var == an->var &&
                fabs(p->anchors[b].value - an->value) > 1e-12)
                return -1;
    }
    return 0;
}

double MonotoneQp_objective(const MonotoneQpProblem *p,
                            const double *x)
{
    if (p == NULL || x == NULL) return DBL_MAX;
    double energy = 0.0;
    for (size_t r = 0; r < p->nrows; r++) {
        const MonotoneQpRow *row = &p->rows[r];
        double residual = -row->target;
        for (int32_t j = 0; j < row->count; j++) {
            const MonotoneQpCoeff *c = &p->coeff[row->first + (size_t)j];
            residual += c->value * x[c->var];
        }
        energy += 0.5 * row->weight * residual * residual;
    }
    return energy;
}

static void objective_gradient(const MonotoneQpProblem *p,
                               const double *x, double *gradient)
{
    memset(gradient, 0, p->nvar * sizeof(double));
    for (size_t r = 0; r < p->nrows; r++) {
        const MonotoneQpRow *row = &p->rows[r];
        double residual = -row->target;
        for (int32_t j = 0; j < row->count; j++) {
            const MonotoneQpCoeff *c = &p->coeff[row->first + (size_t)j];
            residual += c->value * x[c->var];
        }
        residual *= row->weight;
        for (int32_t j = 0; j < row->count; j++) {
            const MonotoneQpCoeff *c = &p->coeff[row->first + (size_t)j];
            gradient[c->var] += c->value * residual;
        }
    }
}

static double feasibility(const MonotoneQpProblem *p, const double *x,
                          double *out_min_slack)
{
    double max_violation = 0.0;
    double min_slack = DBL_MAX;
    for (size_t e = 0; e < p->nbounds; e++) {
        const MonotoneQpBound *b = &p->bounds[e];
        double slack = x[b->hi] - x[b->lo] - b->lower;
        if (slack < min_slack) min_slack = slack;
        if (-slack > max_violation) max_violation = -slack;
    }
    if (p->nbounds == 0) min_slack = 0.0;
    if (max_violation < 0.0) max_violation = 0.0;
    if (out_min_slack != NULL) *out_min_slack = min_slack;
    return max_violation;
}

static size_t max_row_count(const MonotoneQpProblem *p)
{
    size_t n = 1;
    for (size_t r = 0; r < p->nrows; r++)
        if ((size_t)p->rows[r].count > n) n = (size_t)p->rows[r].count;
    return n;
}

/*
 * Solve the equality-constrained search-direction problem.  Active difference
 * rows imply p_hi == p_lo and are represented by union-find blocks.  A block
 * containing an anchor is eliminated.
 */
static int solve_direction_reg(Arena_T arena,
                               const MonotoneQpProblem *p,
                               const uint8_t *active,
                               const double *gradient,
                               int32_t *parent,
                               uint8_t *rank,
                               int32_t *root_to_free,
                               double *direction,
                               double *out_linear_residual,
                               int *out_nfree,
                               double diag_reg)
{
    const size_t n = p->nvar;
    for (size_t i = 0; i < n; i++) {
        parent[i] = (int32_t)i;
        rank[i] = 0;
        root_to_free[i] = -2;
    }
    for (size_t e = 0; e < p->nbounds; e++)
        if (active[e])
            mq_union(parent, rank, p->bounds[e].lo, p->bounds[e].hi);
    for (size_t i = 0; i < n; i++) parent[i] = mq_find(parent, (int32_t)i);

    for (size_t a = 0; a < p->nanchors; a++)
        root_to_free[parent[p->anchors[a].var]] = -1;

    int nfree = 0;
    for (size_t i = 0; i < n; i++) {
        int32_t r = parent[i];
        if ((size_t)r != i || root_to_free[r] == -1) continue;
        root_to_free[r] = nfree++;
    }
    for (size_t i = 0; i < n; i++)
        if (root_to_free[parent[i]] == -2)
            root_to_free[parent[i]] = nfree++;

    memset(direction, 0, n * sizeof(double));
    *out_nfree = nfree;
    *out_linear_residual = 0.0;
    if (nfree == 0) return 0;

    Arena_Mark mark = Arena_save(arena);
    double *rhs = (double *)ARENA_CALLOC(arena, (size_t)nfree, sizeof(double));
    double *solution = (double *)ARENA_CALLOC(arena, (size_t)nfree, sizeof(double));
    double *matvec = (double *)ARENA_CALLOC(arena, (size_t)nfree, sizeof(double));

    for (size_t i = 0; i < n; i++) {
        int fi = root_to_free[parent[i]];
        if (fi >= 0) rhs[fi] -= gradient[i];
    }

    size_t max_triplets = 0;
    for (size_t r = 0; r < p->nrows; r++) {
        size_t k = (size_t)p->rows[r].count;
        if (k > 0 && k > (SIZE_MAX / (k + 1))) {
            Arena_restore(arena, mark);
            return -1;
        }
        size_t add = k * (k + 1) / 2;
        if (max_triplets > SIZE_MAX - add) {
            Arena_restore(arena, mark);
            return -1;
        }
        max_triplets += add;
    }
    if (diag_reg > 0.0) max_triplets += (size_t)nfree;
    if (max_triplets == 0 || max_triplets > (size_t)INT_MAX) {
        Arena_restore(arena, mark);
        return -1;
    }

    int *rows = (int *)ARENA_ALLOC(arena, max_triplets * sizeof(int));
    int *cols = (int *)ARENA_ALLOC(arena, max_triplets * sizeof(int));
    double *vals = (double *)ARENA_ALLOC(arena, max_triplets * sizeof(double));

    size_t max_coeff = max_row_count(p);
    int32_t *agg_index = (int32_t *)ARENA_ALLOC(arena, max_coeff * sizeof(int32_t));
    double *agg_value = (double *)ARENA_ALLOC(arena, max_coeff * sizeof(double));

    size_t nt = 0;
    for (size_t r = 0; r < p->nrows; r++) {
        const MonotoneQpRow *row = &p->rows[r];
        if (row->weight == 0.0) continue;
        size_t na = 0;
        for (int32_t j = 0; j < row->count; j++) {
            const MonotoneQpCoeff *c = &p->coeff[row->first + (size_t)j];
            int fi = root_to_free[parent[c->var]];
            if (fi < 0 || c->value == 0.0) continue;
            size_t k = 0;
            while (k < na && agg_index[k] != fi) k++;
            if (k == na) {
                agg_index[na] = fi;
                agg_value[na] = c->value;
                na++;
            } else {
                agg_value[k] += c->value;
            }
        }
        for (size_t a = 0; a < na; a++) {
            if (agg_value[a] == 0.0) continue;
            for (size_t b = 0; b <= a; b++) {
                if (agg_value[b] == 0.0) continue;
                int ia = agg_index[a], ib = agg_index[b];
                rows[nt] = ia >= ib ? ia : ib;
                cols[nt] = ia >= ib ? ib : ia;
                vals[nt] = row->weight * agg_value[a] * agg_value[b];
                nt++;
            }
        }
    }
    /*
     * Optional Tikhonov diagonal.  The strip normal matrix is only barely
     * SPD in its weakest directions (an intercept whose every candidate the
     * robust loop rejected keeps one floor-weight row), and TAUCS LLT fails
     * outright on a non-positive pivot.  A diagonal at 1e-7 is invisible to
     * every real direction (weakest genuine row weights ~1e-4) and turns
     * the no-op directions definitely positive.
     */
    if (diag_reg > 0.0)
        for (int i = 0; i < nfree; i++) {
            rows[nt] = i;
            cols[nt] = i;
            vals[nt] = diag_reg;
            nt++;
        }
    if (nt == 0 || nt > (size_t)INT_MAX) {
        fprintf(stderr, "[monotone_qp] solve_direction: triplet count %zu\n",
                nt);
        Arena_restore(arena, mark);
        return -1;
    }

    if (Sparse_solve_sym(nfree, (int)nt, rows, cols, vals,
                         rhs, solution, SPARSE_SPD) != 0) {
        fprintf(stderr,
                "[monotone_qp] solve_direction: Sparse_solve_sym failed "
                "(nfree=%d nt=%zu)\n", nfree, nt);
        Arena_restore(arena, mark);
        return -1;
    }

    for (size_t t = 0; t < nt; t++) {
        int r = rows[t], c = cols[t];
        matvec[r] += vals[t] * solution[c];
        if (r != c) matvec[c] += vals[t] * solution[r];
    }
    double r2 = 0.0, b2 = 0.0;
    for (int i = 0; i < nfree; i++) {
        double residual = matvec[i] - rhs[i];
        r2 += residual * residual;
        b2 += rhs[i] * rhs[i];
    }
    *out_linear_residual = sqrt(r2) / fmax(1.0, sqrt(b2));

    for (size_t i = 0; i < n; i++) {
        int fi = root_to_free[parent[i]];
        direction[i] = fi >= 0 ? solution[fi] : 0.0;
    }
    Arena_restore(arena, mark);
    return 0;
}

static double equality_dot(const EqualityRow *a, const EqualityRow *b)
{
    double value = 0.0;
    for (int i = 0; i < a->n; i++)
        for (int j = 0; j < b->n; j++)
            if (a->v[i] == b->v[j]) value += a->c[i] * b->c[j];
    return value;
}

/*
 * At a zero search direction, recover equality multipliers from
 *
 *     gradient + M^T lambda = 0.
 *
 * For lower bounds written x_hi-x_lo-lower >= 0, valid active multipliers
 * satisfy lambda <= 0.  Anchors have unrestricted multipliers.
 */
static int multiplier_check(Arena_T arena,
                            const MonotoneQpProblem *p,
                            const uint8_t *active,
                            const double *gradient,
                            double tolerance,
                            int32_t *out_drop,
                            double *out_stationarity)
{
    size_t nactive = 0;
    for (size_t e = 0; e < p->nbounds; e++) if (active[e]) nactive++;
    size_t ne = nactive + p->nanchors;
    *out_drop = -1;
    *out_stationarity = 0.0;

    if (ne == 0) {
        double maxg = 0.0;
        for (size_t i = 0; i < p->nvar; i++)
            if (fabs(gradient[i]) > maxg) maxg = fabs(gradient[i]);
        *out_stationarity = maxg;
        return 0;
    }
    if (ne > (size_t)INT_MAX) return -1;

    Arena_Mark mark = Arena_save(arena);
    EqualityRow *eq = (EqualityRow *)ARENA_ALLOC(arena, ne * sizeof(EqualityRow));
    size_t q = 0;
    for (size_t e = 0; e < p->nbounds; e++) {
        if (!active[e]) continue;
        eq[q].v[0] = p->bounds[e].lo;
        eq[q].c[0] = -1.0;
        eq[q].v[1] = p->bounds[e].hi;
        eq[q].c[1] = 1.0;
        eq[q].n = 2;
        eq[q].bound = (int32_t)e;
        q++;
    }
    for (size_t a = 0; a < p->nanchors; a++) {
        eq[q].v[0] = p->anchors[a].var;
        eq[q].c[0] = 1.0;
        eq[q].v[1] = -1;
        eq[q].c[1] = 0.0;
        eq[q].n = 1;
        eq[q].bound = -1;
        q++;
    }

    /*
     * M M^T is nonzero only where two equalities share a variable, and every
     * equality touches at most two.  Enumerating all ne(ne+1)/2 ordered pairs
     * to discover that costs O(ne^2) time and, worse, O(ne^2) triplet slots --
     * which capped the active set near ne = 65535 long before any real problem
     * needed to stop.  Walk the variable -> equality adjacency instead, so the
     * work and the storage both scale with the actual sparsity.
     */
    size_t *vhead = (size_t *)ARENA_CALLOC(arena, p->nvar + 1, sizeof(size_t));
    size_t ninc = 0;
    for (size_t i = 0; i < ne; i++)
        for (int k = 0; k < eq[i].n; k++) {
            vhead[(size_t)eq[i].v[k] + 1]++;
            ninc++;
        }
    for (size_t v = 0; v < p->nvar; v++) vhead[v + 1] += vhead[v];
    int32_t *vadj = (int32_t *)ARENA_ALLOC(
        arena, (ninc ? ninc : 1) * sizeof(int32_t));
    size_t *vcur = (size_t *)ARENA_ALLOC(arena, p->nvar * sizeof(size_t));
    for (size_t v = 0; v < p->nvar; v++) vcur[v] = vhead[v];
    for (size_t i = 0; i < ne; i++)
        for (int k = 0; k < eq[i].n; k++)
            vadj[vcur[(size_t)eq[i].v[k]]++] = (int32_t)i;

    /* Upper bound on stored triplets: the diagonal, plus one entry per
     * co-incident pair on each shared variable. */
    size_t cap = ne;
    for (size_t v = 0; v < p->nvar; v++) {
        size_t d = vhead[v + 1] - vhead[v];
        if (d > 1) cap += d * (d - 1) / 2;
    }
    if (cap > (size_t)INT_MAX) {
        Arena_restore(arena, mark);
        return -1;
    }
    int *rows = (int *)ARENA_ALLOC(arena, (cap ? cap : 1) * sizeof(int));
    int *cols = (int *)ARENA_ALLOC(arena, (cap ? cap : 1) * sizeof(int));
    double *vals = (double *)ARENA_ALLOC(arena, (cap ? cap : 1) * sizeof(double));
    double *rhs = (double *)ARENA_ALLOC(arena, ne * sizeof(double));
    double *lambda = (double *)ARENA_CALLOC(arena, ne, sizeof(double));
    /* Marks the last row a column was written for, so a pair sharing two
     * variables still contributes exactly one triplet. */
    size_t *seen = (size_t *)ARENA_ALLOC(arena, ne * sizeof(size_t));
    for (size_t i = 0; i < ne; i++) seen[i] = ne;

    size_t nt = 0;
    for (size_t i = 0; i < ne; i++) {
        double mg = 0.0;
        for (int k = 0; k < eq[i].n; k++)
            mg += eq[i].c[k] * gradient[eq[i].v[k]];
        rhs[i] = -mg;

        double diag = equality_dot(&eq[i], &eq[i]);
        if (diag != 0.0) {
            rows[nt] = (int)i;
            cols[nt] = (int)i;
            vals[nt] = diag;
            nt++;
        }
        for (int k = 0; k < eq[i].n; k++) {
            size_t v = (size_t)eq[i].v[k];
            for (size_t e = vhead[v]; e < vhead[v + 1]; e++) {
                size_t j = (size_t)vadj[e];
                if (j >= i || seen[j] == i) continue;
                seen[j] = i;
                double dot = equality_dot(&eq[i], &eq[j]);
                if (dot == 0.0) continue;
                rows[nt] = (int)i;
                cols[nt] = (int)j;
                vals[nt] = dot;
                nt++;
            }
        }
    }
    if (nt == 0 ||
        Sparse_solve_sym((int)ne, (int)nt, rows, cols, vals,
                         rhs, lambda, SPARSE_SPD) != 0) {
        Arena_restore(arena, mark);
        return -1;
    }

    double worst = tolerance;
    for (size_t i = 0; i < nactive; i++) {
        if (lambda[i] > worst) {
            worst = lambda[i];
            *out_drop = eq[i].bound;
        }
    }

    double *stationarity =
        (double *)ARENA_ALLOC(arena, p->nvar * sizeof(double));
    memcpy(stationarity, gradient, p->nvar * sizeof(double));
    for (size_t i = 0; i < ne; i++)
        for (int k = 0; k < eq[i].n; k++)
            stationarity[eq[i].v[k]] += eq[i].c[k] * lambda[i];
    double maxs = 0.0;
    for (size_t i = 0; i < p->nvar; i++)
        if (fabs(stationarity[i]) > maxs) maxs = fabs(stationarity[i]);
    *out_stationarity = maxs;

    Arena_restore(arena, mark);
    return 0;
}

static void trace_append(MonotoneQpTrace *trace,
                         int iteration, int nactive,
                         int added, int dropped,
                         double objective, double violation,
                         double direction_norm, double step,
                         double linear_residual)
{
    if (trace == NULL || trace->entry == NULL ||
        trace->count >= trace->capacity) return;
    MonotoneQpTraceEntry *e = &trace->entry[trace->count++];
    e->iteration = iteration;
    e->n_active = nactive;
    e->added_bound = added;
    e->dropped_bound = dropped;
    e->objective = objective;
    e->max_violation = violation;
    e->direction_norm = direction_norm;
    e->step = step;
    e->reduced_linear_residual = linear_residual;
}

/*
 * Anchors-only exact quadratic minimizer: one Newton step from x with an
 * empty active set is the global optimum of a quadratic, and solve_direction
 * already eliminates anchored variables.  No feasible start, no iteration
 * cap, no rc -2/-3/-4 modes -- the caller enforces the bounds afterwards
 * (per-stroke PAVA projection in atlas_strip.c).
 */
int MonotoneQp_solve_ls(Arena_T arena,
                        const MonotoneQpProblem *p,
                        double *x,
                        double *out_linear_residual)
{
    if (arena == NULL || x == NULL) return -1;
    if (validate_problem(p) != 0) {
        fprintf(stderr, "[monotone_qp] solve_ls: problem validation failed\n");
        return -1;
    }
    Arena_Mark mark = Arena_save(arena);
    uint8_t *active =
        (uint8_t *)ARENA_CALLOC(arena, p->nbounds ? p->nbounds : 1, 1);
    double *gradient = (double *)ARENA_ALLOC(arena, p->nvar * sizeof(double));
    double *direction = (double *)ARENA_ALLOC(arena, p->nvar * sizeof(double));
    int32_t *parent = (int32_t *)ARENA_ALLOC(arena, p->nvar * sizeof(int32_t));
    uint8_t *rank = (uint8_t *)ARENA_ALLOC(arena, p->nvar * sizeof(uint8_t));
    int32_t *root_to_free =
        (int32_t *)ARENA_ALLOC(arena, p->nvar * sizeof(int32_t));

    for (size_t a = 0; a < p->nanchors; a++)
        x[p->anchors[a].var] = p->anchors[a].value;
    objective_gradient(p, x, gradient);
    double linear_residual = 0.0;
    int nfree = 0;
    int rc = solve_direction_reg(arena, p, active, gradient, parent, rank,
                                 root_to_free, direction, &linear_residual,
                                 &nfree, 1e-7);
    if (rc == 0)
        for (size_t i = 0; i < p->nvar; i++) x[i] += direction[i];
    if (out_linear_residual != NULL)
        *out_linear_residual = linear_residual;
    Arena_restore(arena, mark);
    return rc;
}

int MonotoneQp_solve(Arena_T arena,
                     const MonotoneQpProblem *p,
                     const MonotoneQpOptions *input_opts,
                     double *x,
                     MonotoneQpTrace *trace,
                     MonotoneQpStats *stats)
{
    if (arena == NULL || x == NULL || validate_problem(p) != 0) return -1;

    MonotoneQpOptions defaults;
    MonotoneQpOptions_default(&defaults);
    const MonotoneQpOptions *opts = input_opts != NULL ? input_opts : &defaults;
    if (opts->max_active_iterations <= 0 ||
        opts->feasibility_tolerance < 0.0 ||
        opts->direction_tolerance <= 0.0 ||
        opts->multiplier_tolerance < 0.0) return -1;

    MonotoneQpStats zero;
    memset(&zero, 0, sizeof zero);
    if (stats != NULL) *stats = zero;
    if (trace != NULL) trace->count = 0;

    for (size_t a = 0; a < p->nanchors; a++)
        if (fabs(x[p->anchors[a].var] - p->anchors[a].value) >
            opts->feasibility_tolerance)
            return -2;
    double min_slack = 0.0;
    if (feasibility(p, x, &min_slack) > opts->feasibility_tolerance)
        return -2;

    Arena_Mark mark = Arena_save(arena);
    uint8_t *active =
        (uint8_t *)ARENA_CALLOC(arena, p->nbounds ? p->nbounds : 1, 1);
    double *gradient = (double *)ARENA_ALLOC(arena, p->nvar * sizeof(double));
    double *direction = (double *)ARENA_ALLOC(arena, p->nvar * sizeof(double));
    int32_t *parent = (int32_t *)ARENA_ALLOC(arena, p->nvar * sizeof(int32_t));
    uint8_t *rank = (uint8_t *)ARENA_ALLOC(arena, p->nvar * sizeof(uint8_t));
    int32_t *root_to_free =
        (int32_t *)ARENA_ALLOC(arena, p->nvar * sizeof(int32_t));

    /*
     * Seed the active set with the bounds that are already tight, but never
     * with a cycle.  The multipliers are recovered from M M^T, and the
     * incidence vectors around any cycle of difference constraints sum to
     * zero, so a cycle makes that Gram singular and the SPD solve fails (-3).
     * Chain problems can never form one; a general difference graph -- stroke
     * monotonicity plus radial ordering, where one sample is bounded both
     * along its stroke and against its neighbouring wrap -- forms them
     * routinely.  A skipped bound stays satisfied, and if it ever blocks a
     * step its endpoints are already in one contracted block, where the
     * direction difference is zero and it is never chosen again.
     */
    int nactive = 0;
    for (size_t v = 0; v < p->nvar; v++) { parent[v] = (int32_t)v; rank[v] = 0; }
    for (size_t e = 0; e < p->nbounds; e++) {
        double slack = x[p->bounds[e].hi] - x[p->bounds[e].lo] -
                       p->bounds[e].lower;
        if (fabs(slack) > opts->feasibility_tolerance) continue;
        if (!mq_union(parent, rank, p->bounds[e].lo, p->bounds[e].hi))
            continue;
        active[e] = 1;
        nactive++;
    }

    MonotoneQpStats local;
    memset(&local, 0, sizeof local);
    local.objective_initial = MonotoneQp_objective(p, x);
    local.max_reduced_linear_residual = 0.0;

    int rc = -4;
    for (int iteration = 0; iteration < opts->max_active_iterations;
         iteration++) {
        objective_gradient(p, x, gradient);

        double linear_residual = 0.0;
        int nfree = 0;
        if (solve_direction_reg(arena, p, active, gradient,
                                parent, rank, root_to_free, direction,
                                &linear_residual, &nfree, 0.0) != 0) {
            rc = -3;
            break;
        }
        local.spd_solves++;
        if (linear_residual > local.max_reduced_linear_residual)
            local.max_reduced_linear_residual = linear_residual;
        if (linear_residual > opts->linear_residual_tolerance &&
            opts->verbose)
            fprintf(stderr, "[monotone_qp] warning: reduced residual %.3e\n",
                    linear_residual);

        double p2 = 0.0, x2 = 0.0;
        for (size_t i = 0; i < p->nvar; i++) {
            p2 += direction[i] * direction[i];
            x2 += x[i] * x[i];
        }
        double pnorm = sqrt(p2);
        double small = opts->direction_tolerance * fmax(1.0, sqrt(x2));
        double violation = feasibility(p, x, NULL);

        if (pnorm <= small || nfree == 0) {
            int32_t drop = -1;
            double stationarity = 0.0;
            if (multiplier_check(arena, p, active, gradient,
                                 opts->multiplier_tolerance,
                                 &drop, &stationarity) != 0) {
                rc = -3;
                break;
            }
            local.multiplier_solves++;
            if (drop >= 0) {
                active[drop] = 0;
                nactive--;
                local.dropped_total++;
                trace_append(trace, iteration, nactive, -1, drop,
                             MonotoneQp_objective(p, x), violation,
                             pnorm, 0.0, linear_residual);
                continue;
            }
            local.stationarity_residual = stationarity;
            local.iterations = iteration + 1;
            trace_append(trace, iteration, nactive, -1, -1,
                         MonotoneQp_objective(p, x), violation,
                         pnorm, 0.0, linear_residual);
            rc = 0;
            break;
        }

        double alpha = 1.0;
        int32_t blocker = -1;
        for (size_t e = 0; e < p->nbounds; e++) {
            if (active[e]) continue;
            const MonotoneQpBound *b = &p->bounds[e];
            double delta = direction[b->hi] - direction[b->lo];
            if (delta >= -opts->direction_tolerance) continue;
            double slack = x[b->hi] - x[b->lo] - b->lower;
            double candidate = slack / (-delta);
            if (candidate < alpha) {
                alpha = fmax(0.0, candidate);
                blocker = (int32_t)e;
            }
        }

        for (size_t i = 0; i < p->nvar; i++)
            x[i] += alpha * direction[i];

        if (blocker >= 0 && alpha < 1.0 -
            opts->feasibility_tolerance) {
            const MonotoneQpBound *b = &p->bounds[blocker];
            double error = b->lower - (x[b->hi] - x[b->lo]);
            int32_t rhi = parent[b->hi], rlo = parent[b->lo];
            int hi_anchored = root_to_free[rhi] < 0;
            int32_t shift_root = hi_anchored ? rlo : rhi;
            double shift = hi_anchored ? -error : error;
            for (size_t i = 0; i < p->nvar; i++)
                if (parent[i] == shift_root) x[i] += shift;
            active[blocker] = 1;
            nactive++;
            local.added_total++;
        } else {
            blocker = -1;
        }

        violation = feasibility(p, x, NULL);
        trace_append(trace, iteration, nactive, blocker, -1,
                     MonotoneQp_objective(p, x), violation,
                     pnorm, alpha, linear_residual);
        if (opts->verbose)
            fprintf(stderr,
                    "[monotone_qp] it=%d active=%d E=%.9g step=%.3g "
                    "|p|=%.3g violation=%.3g\n",
                    iteration, nactive, MonotoneQp_objective(p, x),
                    alpha, pnorm, violation);
    }

    local.objective_final = MonotoneQp_objective(p, x);
    local.active_final = nactive;
    local.max_violation = feasibility(p, x, &local.min_slack);
    if (local.iterations == 0 && rc != -4)
        local.iterations = local.spd_solves;
    if (stats != NULL) *stats = local;
    Arena_restore(arena, mark);
    return rc;
}
