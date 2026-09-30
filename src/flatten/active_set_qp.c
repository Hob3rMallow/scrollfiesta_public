/* active_set_qp.c -- Goldfarb-Idnani dual active set on one TAUCS factor. */
#include "active_set_qp.h"
#include "sparse_solve.h"
#include "convex_qp.h"
#include "../common/ves_platform.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Working set: the dense Schur factor M and its Cholesky are cap^2 doubles
 * each.  It starts at ASQP_INITIAL_ACTIVE and doubles on demand up to
 * ASQP_MAX_ACTIVE (1 GB for both at the limit): a 211k-coordinate closure
 * level legitimately binds more than 2,048 witness half-planes. */
#define ASQP_MAX_ACTIVE 8192
/* A WALL-CLOCK BUDGET, because the step guard is on the wrong axis.  It allows
 * 8*rows + 1024 steps -- 856,808 on the 21x5x5's restoration -- so a solve whose working
 * set is heading past ASQP_MAX_ACTIVE never trips it: it grinds to the SIZE cap and is
 * then refused outright, discarding everything.  Measured 2026-09-18 on that pile: five
 * solves of 10,391 to 15,013 steps, the worst 4,773 s (9.5 s of factorisation and 2,800 s
 * of triangular solves) ending REFUSED, together about six hours of a fifteen-hour repair
 * that produced nothing.  The 4x5x5's worst is 415 steps, so this never fires there.
 * Refusing on time reaches the caller's existing fallback sooner and loses nothing that
 * the size cap was not about to throw away.  0 disables it.
 *
 * CALIBRATED 2026-09-18 by replaying the saved cases, because 600
 * was too tight and refused solves that converge.  A refused level's own case, n 120,564 /
 * 144,064 rows: the active set SOLVES it in 1,281 s with 5,118 rows binding, to a maximum
 * violation of 2.1e-10 -- under the size cap, so nothing about it is pathological.  OSQP on
 * the same case takes 346 s but lands 1.9e-6 infeasible and 3,090 vox away with a LOWER
 * objective: it cuts the corner rather than converging further, so a refusal is not a free
 * substitution.  At 1,800 the exact solve is kept and the genuinely hopeless level -- the one
 * that ran 4,773 s to 15,013 steps and was refused on SIZE -- is still cut off 2.6x sooner. */
#define ASQP_MAX_SECONDS 1800.0
static double asqp_max_seconds = ASQP_MAX_SECONDS;   /* the selftest shrinks it to force the budget */
/* Per thread: concurrent solves each keep their own. */
static VES_THREAD_LOCAL double asqp_deadline = 0;
void ActiveSetQp_deadline(double deadline) { asqp_deadline = deadline > 0 && isfinite(deadline) ? deadline : 0; }
static int asqp_expired(void) { return asqp_deadline > 0 && ves_clock_sec() >= asqp_deadline; }
#define ASQP_INITIAL_ACTIVE 256
static int asqp_initial_cap = ASQP_INITIAL_ACTIVE;   /* the selftest shrinks it to force growth */
/* A constraint whose projected normal is this small relative to its own
 * H-norm is linearly dependent on the working set: dual step only. */
#define ASQP_DEPENDENT 1e-12
/* Rows violated by less than this (relative to 1+|u|) are satisfied. */
#define ASQP_FEASIBLE 1e-9
/* Ridge used only when the exact factorisation fails, relative to each
 * coordinate's own diagonal; a coordinate with no curvature takes this
 * share of the largest entry instead.  A well-posed H is factored exactly. */
#define ASQP_RIDGE 1e-8

typedef struct AsqpRow { int first, count; } AsqpRow;   /* CSR row into ids/vals */

typedef struct Asqp {
    int n, m;
    /* constraint rows, CSR */
    AsqpRow *row; int *ids; double *vals; const double *upper;
    /* factor */
    SparseFactor_T factor;
    /* iterate */
    double *x, *lambda;            /* x[n]; lambda[cap] aligned with W */
    int *W, *inW, nW, cap, max_cap; /* working set, membership flags, current and final capacity */
    double *M, *L;                 /* dense cap x cap: M = N'H^-1N and its Cholesky */
    double *yp, *z, *rhs, *s, *r;  /* n, n, n, cap, cap */
    int hn; const int *hr, *hc; const double *half, *g;   /* H (true, lower) and g for the running objective */
    double value_limit;
    int *warm, nwarm, warm_cap;    /* warm-start rows in; final working set out */
    ActiveSetQpStats st;
} Asqp;

static double asqp_row_dot(const Asqp *q, int i, const double *v)
{
    double d = 0.0;
    for (int k = q->row[i].first; k < q->row[i].first + q->row[i].count; k++) d += q->vals[k] * v[q->ids[k]];
    return d;
}

static void asqp_row_scatter(const Asqp *q, int i, double w, double *v)
{
    for (int k = q->row[i].first; k < q->row[i].first + q->row[i].count; k++) v[q->ids[k]] += w * q->vals[k];
}

static int asqp_solve(Asqp *q, const double *b, double *x)
{
    double t0 = ves_clock_sec();
    int rc = Sparse_factor_solve(q->factor, b, x);
    q->st.solve_seconds += ves_clock_sec() - t0;
    if (rc) return 0;
    for (int i = 0; i < q->n; i++) if (!isfinite(x[i])) return 0;
    return 1;
}

/* Dense lower Cholesky of the leading k x k block of M into L (row-major, lda=cap). */
static int asqp_chol_full(Asqp *q, int k)
{
    const int lda = q->cap;
    for (int i = 0; i < k; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = q->M[(size_t)i*lda + j];
            for (int t = 0; t < j; t++) sum -= q->L[(size_t)i*lda + t] * q->L[(size_t)j*lda + t];
            if (i == j) { if (!(sum > 0.0) || !isfinite(sum)) return 0; q->L[(size_t)i*lda + i] = sqrt(sum); }
            else q->L[(size_t)i*lda + j] = sum / q->L[(size_t)j*lda + j];
        }
    }
    return 1;
}

/* Append row k (M[k][0..k] already filled) to an existing k x k Cholesky. */
static int asqp_chol_append(Asqp *q, int k)
{
    const int lda = q->cap;
    for (int j = 0; j <= k; j++) {
        double sum = q->M[(size_t)k*lda + j];
        for (int t = 0; t < j; t++) sum -= q->L[(size_t)k*lda + t] * q->L[(size_t)j*lda + t];
        if (j == k) { if (!(sum > 0.0) || !isfinite(sum)) return 0; q->L[(size_t)k*lda + k] = sqrt(sum); }
        else q->L[(size_t)k*lda + j] = sum / q->L[(size_t)j*lda + j];
    }
    return 1;
}

/* Solve M r = s through L for the leading k block. */
static void asqp_chol_solve(const Asqp *q, int k, const double *s, double *r)
{
    const int lda = q->cap;
    for (int i = 0; i < k; i++) {
        double sum = s[i];
        for (int t = 0; t < i; t++) sum -= q->L[(size_t)i*lda + t] * r[t];
        r[i] = sum / q->L[(size_t)i*lda + i];
    }
    for (int i = k - 1; i >= 0; i--) {
        double sum = r[i];
        for (int t = i + 1; t < k; t++) sum -= q->L[(size_t)t*lda + i] * r[t];
        r[i] = sum / q->L[(size_t)i*lda + i];
    }
}

/* Delete slot j from the k x k Cholesky factor L of M in place: the rows
 * below j move up, which leaves one super-diagonal entry per moved row, and
 * Givens rotations applied from the right chase it off the end.  Column j
 * is never physically removed; the rotations carry its weight rightwards.
 * O((k-j)^2) against O(k^3) for a refactorisation -- at a working set of a
 * thousand rows the refactorisation was the whole cost of a drop (840 s of
 * one 938-s solve, 2026-09-15). */
static int asqp_chol_delete(Asqp *q, int j, int k)
{
    const int lda = q->cap; double *L = q->L;
    for (int i = j; i + 1 < k; i++) memmove(&L[(size_t)i*lda], &L[(size_t)(i+1)*lda], (size_t)(i + 2) * sizeof *L);
    const int kk = k - 1;
    for (int c = j; c < kk; c++) {
        double a = L[(size_t)c*lda + c], b = L[(size_t)c*lda + c + 1];
        double r = hypot(a, b);
        if (!(r > 0.0) || !isfinite(r)) return 0;
        double cs = a / r, sn = b / r;
        for (int i = c; i < kk; i++) {
            double x = L[(size_t)i*lda + c], y = L[(size_t)i*lda + c + 1];
            L[(size_t)i*lda + c] = cs * x + sn * y;
            L[(size_t)i*lda + c + 1] = cs * y - sn * x;
        }
    }
    return 1;
}

/* Remove working-set slot j: compact W, lambda and M; update the factor. */
static int asqp_drop(Asqp *q, int j)
{
    const int lda = q->cap, k = q->nW;
    q->inW[q->W[j]] = 0;
    for (int i = j; i + 1 < k; i++) { q->W[i] = q->W[i+1]; q->lambda[i] = q->lambda[i+1]; }
    for (int i = 0; i < k; i++) if (i != j) {
        int ii = i < j ? i : i - 1;
        for (int c = 0; c < k; c++) if (c != j) {
            int cc = c < j ? c : c - 1;
            q->M[(size_t)ii*lda + cc] = q->M[(size_t)i*lda + c];
        }
    }
    q->nW = k - 1; q->st.dropped++;
    return asqp_chol_delete(q, j, k) || asqp_chol_full(q, q->nW);
}

/* Double the working-set capacity, re-striding the dense blocks. */
static int asqp_grow(Asqp *q)
{
    const int old = q->cap; int cap = old * 2; if (cap > q->max_cap) cap = q->max_cap; if (cap <= old) return 0;
    double *M = calloc((size_t)cap * (size_t)cap, sizeof *M), *L = calloc((size_t)cap * (size_t)cap, sizeof *L);
    double *lambda = realloc(q->lambda, (size_t)cap * sizeof *lambda);
    if (lambda) q->lambda = lambda;
    double *s = realloc(q->s, (size_t)cap * sizeof *s); if (s) q->s = s;
    double *r = realloc(q->r, (size_t)cap * sizeof *r); if (r) q->r = r;
    int *W = realloc(q->W, (size_t)cap * sizeof *W); if (W) q->W = W;
    if (!M || !L || !lambda || !s || !r || !W) { free(M); free(L); return 0; }
    for (int i = 0; i < q->nW; i++) {
        memcpy(M + (size_t)i * cap, q->M + (size_t)i * old, (size_t)q->nW * sizeof *M);
        memcpy(L + (size_t)i * cap, q->L + (size_t)i * old, (size_t)(i + 1) * sizeof *L);
    }
    free(q->M); free(q->L); q->M = M; q->L = L; q->cap = cap;
    return 1;
}

/* 1/2 x'Hx + g'x at the current iterate. */
static double asqp_value(const Asqp *q)
{
    double xhx = 0.0, gx = 0.0;
    for (int k = 0; k < q->hn; k++) { int i = q->hr[k], j = q->hc[k]; xhx += q->half[k] * (i == j ? q->x[i] * q->x[i] : 2.0 * q->x[i] * q->x[j]); }
    for (int i = 0; i < q->n; i++) gx += q->g[i] * q->x[i];
    return 0.5 * xhx + gx;
}

/* Seed the working set with rows that were binding last time.  M = N'H^-1 N
 * costs one solve per row, lambda = M^-1 (N'x0 - u), rows with a negative
 * multiplier are dropped until the set is dual feasible, and the primal
 * iterate is then ONE more solve, x = H^-1(-g - N lambda).  On the global
 * closure levels consecutive iterations bind nearly the same 500-900 rows,
 * and a cold start paid two solves per step to rediscover them.  A seed
 * that is wrong only costs its solve: the main loop corrects it. */
static int asqp_warm_start(Asqp *q, const double *g)
{
    const int n = q->n;
    int k = 0;
    for (int t = 0; t < q->nwarm; t++) {
        int i = q->warm[t];
        if (i < 0 || i >= q->m || q->inW[i] || q->row[i].count == 0) continue;
        if (k >= q->cap && !asqp_grow(q)) break;
        q->W[k] = i; q->inW[i] = 1; k++;
    }
    if (!k) return 1;
    const int lda = q->cap;
    for (int j = 0; j < k; j++) {
        if (asqp_expired()) { q->st.reason = ASQP_REASON_TIME; return 0; }
        memset(q->rhs, 0, (size_t)n * sizeof *q->rhs);
        asqp_row_scatter(q, q->W[j], 1.0, q->rhs);
        if (!asqp_solve(q, q->rhs, q->yp)) return 0;
        for (int i = 0; i < k; i++) q->M[(size_t)i*lda + j] = asqp_row_dot(q, q->W[i], q->yp);
    }
    for (int i = 0; i < k; i++) for (int j = 0; j < i; j++) {   /* exact symmetry for the Cholesky */
        double v = 0.5 * (q->M[(size_t)i*lda + j] + q->M[(size_t)j*lda + i]);
        q->M[(size_t)i*lda + j] = v; q->M[(size_t)j*lda + i] = v;
    }
    q->nW = k;
    if (!asqp_chol_full(q, k)) {                           /* dependent seed: start cold */
        for (int j = 0; j < k; j++) q->inW[q->W[j]] = 0;
        q->nW = 0; return 1;
    }
    for (;;) {
        k = q->nW;
        for (int j = 0; j < k; j++) q->s[j] = asqp_row_dot(q, q->W[j], q->x) - q->upper[q->W[j]];
        asqp_chol_solve(q, k, q->s, q->lambda);
        int worst = -1; double wv = 0.0;
        for (int j = 0; j < k; j++) if (q->lambda[j] < wv) { wv = q->lambda[j]; worst = j; }
        if (worst < 0) break;
        if (!asqp_drop(q, worst)) { for (int j = 0; j < q->nW; j++) q->inW[q->W[j]] = 0; q->nW = 0; return 1; }
        if (!q->nW) return 1;
    }
    k = q->nW; q->st.warm = k;
    for (int i = 0; i < n; i++) q->rhs[i] = -g[i];
    for (int j = 0; j < k; j++) asqp_row_scatter(q, q->W[j], -q->lambda[j], q->rhs);
    if (!asqp_solve(q, q->rhs, q->x)) return 0;
    return 1;
}

static int asqp_run(Asqp *q, const double *g)
{
    if (asqp_expired()) { q->st.reason = ASQP_REASON_TIME; return 0; }
    const int n = q->n, m = q->m;
    /* unconstrained minimiser */
    for (int i = 0; i < n; i++) q->rhs[i] = -g[i];
    if (!asqp_solve(q, q->rhs, q->x)) { q->st.reason = ASQP_REASON_SOLVE; return 0; }
    if (q->nwarm > 0 && !asqp_warm_start(q, g)) { q->st.reason = asqp_expired() ? ASQP_REASON_TIME : ASQP_REASON_SOLVE; return 0; }
    int guard = 0, max_steps = 8 * (m > 16 ? m : 16) + 1024;
    double t_start = asqp_max_seconds > 0.0 ? ves_clock_sec() : 0.0;
    for (;;) {
        /* most violated inactive constraint */
        if (asqp_expired()) { q->st.reason = ASQP_REASON_TIME; return 0; }
        int p = -1; double worst = 0.0;
        for (int i = 0; i < m; i++) {
            if (q->inW[i]) continue;
            double v = asqp_row_dot(q, i, q->x) - q->upper[i];
            if (!isfinite(v)) { q->st.reason = ASQP_REASON_SOLVE; return 0; }
            double tol = ASQP_FEASIBLE * (1.0 + fabs(q->upper[i]));
            if (v > tol && v > worst) { worst = v; p = i; }
        }
        if (p < 0) break;                                    /* primal feasible + dual feasible = optimal */
        if (q->value_limit < INFINITY && (q->st.iterations & 15) == 0 && asqp_value(q) > q->value_limit) { q->st.reason = ASQP_REASON_ILLPOSED; return 0; }
        if (q->row[p].count == 0) { q->st.reason = ASQP_REASON_INFEASIBLE; return 0; }   /* 0 <= u_p < 0: infeasible row */
        /* y_p = H^-1 a_p */
        memset(q->rhs, 0, (size_t)n * sizeof *q->rhs);
        asqp_row_scatter(q, p, 1.0, q->rhs);
        if (!asqp_solve(q, q->rhs, q->yp)) { q->st.reason = ASQP_REASON_SOLVE; return 0; }
        double app = asqp_row_dot(q, p, q->yp);              /* a_p' H^-1 a_p */
        if (!(app > 0.0) || !isfinite(app)) { q->st.reason = ASQP_REASON_DENSE; return 0; }
        double lambda_p = 0.0;
        for (;;) {
            if (++guard > max_steps) { q->st.reason = ASQP_REASON_STEPS; return 0; }
            if (asqp_expired()) { q->st.reason = ASQP_REASON_TIME; return 0; }
            if (asqp_max_seconds > 0.0 && (guard & 63) == 0 &&
                ves_clock_sec() - t_start > asqp_max_seconds) { q->st.reason = ASQP_REASON_TIME; return 0; }
            q->st.iterations++;
            int k = q->nW;
            /* s = N' y_p ; r = M^-1 s ; z = H^-1 (a_p - N r) */
            for (int j = 0; j < k; j++) q->s[j] = asqp_row_dot(q, q->W[j], q->yp);
            if (k) asqp_chol_solve(q, k, q->s, q->r);
            double apz;
            if (k) {
                memset(q->rhs, 0, (size_t)n * sizeof *q->rhs);
                asqp_row_scatter(q, p, 1.0, q->rhs);
                for (int j = 0; j < k; j++) asqp_row_scatter(q, q->W[j], -q->r[j], q->rhs);
                if (!asqp_solve(q, q->rhs, q->z)) { q->st.reason = ASQP_REASON_SOLVE; return 0; }
                apz = asqp_row_dot(q, p, q->z);
            } else { memcpy(q->z, q->yp, (size_t)n * sizeof *q->z); apz = app; }
            /* dual step length: first multiplier to reach zero along -r */
            double t1 = INFINITY; int drop = -1;
            for (int j = 0; j < k; j++) if (q->r[j] > 0.0) {
                double t = q->lambda[j] / q->r[j];
                if (t < t1) { t1 = t; drop = j; }
            }
            int dependent = !(apz > ASQP_DEPENDENT * app);
            double vp = asqp_row_dot(q, p, q->x) - q->upper[p];
            double t2 = dependent ? INFINITY : vp / apz;
            if (dependent && drop < 0) { q->st.reason = ASQP_REASON_INFEASIBLE; return 0; }   /* no primal step and nothing to drop */
            double t = t1 < t2 ? t1 : t2;
            if (!isfinite(t) || t < 0.0) { q->st.reason = ASQP_REASON_DENSE; return 0; }
            if (!dependent) for (int i = 0; i < n; i++) q->x[i] -= t * q->z[i];
            for (int j = 0; j < k; j++) q->lambda[j] -= t * q->r[j];
            lambda_p += t;
            if (t1 < t2 || dependent) {
                if (!asqp_drop(q, drop)) { q->st.reason = ASQP_REASON_DENSE; return 0; }
                continue;                                    /* same p, smaller working set */
            }
            /* p becomes active */
            if (k >= q->cap && !asqp_grow(q)) { q->st.reason = ASQP_REASON_CAP; return 0; }
            { const int lda = q->cap;
              for (int j = 0; j < k; j++) { q->M[(size_t)k*lda + j] = q->s[j]; q->M[(size_t)j*lda + k] = q->s[j]; }
              q->M[(size_t)k*lda + k] = app; }
            if (!asqp_chol_append(q, k)) { q->st.reason = ASQP_REASON_DENSE; return 0; }
            q->W[k] = p; q->lambda[k] = lambda_p; q->inW[p] = 1; q->nW = k + 1;
            break;
        }
    }
    if (q->value_limit < INFINITY && asqp_value(q) > q->value_limit) { q->st.reason = ASQP_REASON_ILLPOSED; return 0; }
    q->st.active = q->nW;
    if (q->warm && q->warm_cap > 0) { int c = q->nW < q->warm_cap ? q->nW : q->warm_cap; memcpy(q->warm, q->W, (size_t)c * sizeof *q->warm); q->nwarm = c; }
    return 1;
}

struct ActiveSetQpFactor {
    int n, hn;
    SparseFactor_T factor;
    double factor_seconds;
    int ridge_retries;
    int *hr, *hc; double *half;   /* H kept for the running objective */
    double value_limit;
};

void ActiveSetQp_limit(ActiveSetQpFactor_T factor, double limit) { if (factor) factor->value_limit = limit; }

int ActiveSetQp_factor(int n, int hn, const int *hr, const int *hc, const double *hv,
                       ActiveSetQpFactor_T *out, ActiveSetQpStats *stats)
{
    if (stats) { memset(stats, 0, sizeof *stats); stats->reason = ASQP_REASON_INPUT; }
    if (out) *out = NULL;
    if (asqp_expired()) { if (stats) stats->reason = ASQP_REASON_TIME; return 0; }
    if (n <= 0 || hn < 0 || !hr || !hc || !hv || !out) return 0;
    for (int k = 0; k < hn; k++) if (hr[k] < 0 || hr[k] >= n || hc[k] < 0 || hc[k] > hr[k] || !isfinite(hv[k])) return 0;
    double *half = malloc((size_t)(hn ? hn : 1) * sizeof *half);
    double *diag = calloc((size_t)n, sizeof *diag);
    ActiveSetQpFactor_T f = calloc(1, sizeof *f);
    int factored = 0;
    if (!half || !diag || !f) { if (stats) stats->reason = ASQP_REASON_MEMORY; goto done; }
    f->n = n; f->hn = hn; f->value_limit = INFINITY;
    /* true Hessian from the doubled convention */
    for (int k = 0; k < hn; k++) { half[k] = 0.5 * hv[k]; if (hr[k] == hc[k]) diag[hr[k]] += half[k]; }
    f->hr = malloc((size_t)(hn ? hn : 1) * sizeof *f->hr); f->hc = malloc((size_t)(hn ? hn : 1) * sizeof *f->hc); f->half = malloc((size_t)(hn ? hn : 1) * sizeof *f->half);
    if (!f->hr || !f->hc || !f->half) { if (stats) stats->reason = ASQP_REASON_MEMORY; goto done; }
    memcpy(f->hr, hr, (size_t)hn * sizeof *f->hr); memcpy(f->hc, hc, (size_t)hn * sizeof *f->hc); memcpy(f->half, half, (size_t)hn * sizeof *f->half);
    /* Factor: exact first, ridged only when that fails. */
    {
        double t0 = ves_clock_sec(), dmax = 0.0;
        for (int i = 0; i < n; i++) if (diag[i] > dmax) dmax = diag[i];
        /* Attempt 0 factors a well-posed H EXACTLY, ridging only coordinates
         * with no curvature at all: the C/Python parity check holds the step
         * to 1e-6, and on an ill-conditioned closure step even a 1e-8
         * relative ridge moved it by 4e-6.  A merely semidefinite H, which
         * TAUCS now refuses cleanly (its multifrontal unwind's double free is
         * fixed), is retried with a ridge RELATIVE to each coordinate's own
         * curvature, as OSQP's sigma is after equilibration -- a ridge
         * measured against the largest entry biased the weakest coordinates
         * by 1e-3 of their stiffness.  A dead coordinate takes the largest
         * entry's share so its step stays finite. */
        double floor_ridge = (dmax > 0.0 ? dmax : 1.0) * ASQP_RIDGE;
        int deficient = 0;
        for (int i = 0; i < n; i++) if (!(diag[i] > 1e-12 * dmax)) deficient++;
        const double levels[3] = { 0.0, ASQP_RIDGE, 100.0 * ASQP_RIDGE };
        for (int level = 0; level < 3 && !factored; level++) {
            if (asqp_expired()) break;
            if (level > 0 || deficient) f->ridge_retries++;
            if (level == 0 && !deficient) { factored = Sparse_factor_spd(n, hn, hr, hc, half, &f->factor) == 0; continue; }
            int hn2 = hn + n; int *rr = malloc((size_t)hn2 * sizeof *rr), *cc = malloc((size_t)hn2 * sizeof *cc); double *vv = malloc((size_t)hn2 * sizeof *vv);
            if (rr && cc && vv) {
                memcpy(rr, hr, (size_t)hn * sizeof *rr); memcpy(cc, hc, (size_t)hn * sizeof *cc); memcpy(vv, half, (size_t)hn * sizeof *vv);
                for (int i = 0; i < n; i++) { rr[hn + i] = i; cc[hn + i] = i; vv[hn + i] = diag[i] > 1e-12 * dmax ? levels[level] * diag[i] : floor_ridge; }
                factored = Sparse_factor_spd(n, hn2, rr, cc, vv, &f->factor) == 0;
            }
            free(rr); free(cc); free(vv);
        }
        f->factor_seconds = ves_clock_sec() - t0;
        if (stats) { stats->factor_seconds = f->factor_seconds; stats->ridge_retries = f->ridge_retries; stats->reason = factored ? ASQP_REASON_OK : ASQP_REASON_FACTOR; }
    }
done:
    free(half); free(diag);
    if (!factored) { ActiveSetQp_free(&f); return 0; }
    *out = f;
    return 1;
}

void ActiveSetQp_free(ActiveSetQpFactor_T *factor)
{
    if (!factor || !*factor) return;
    if ((*factor)->factor) Sparse_factor_free(&(*factor)->factor);
    free((*factor)->hr); free((*factor)->hc); free((*factor)->half);
    free(*factor); *factor = NULL;
}

int ActiveSetQp_model(ActiveSetQpFactor_T f, const double *target, double *model)
{
    if (!f || !f->factor || !target || !model) return 0;
    const int n = f->n;
    double *g = malloc((size_t)n * sizeof *g), *x = malloc((size_t)n * sizeof *x);
    int ok = g && x;
    if (ok) {
        for (int i = 0; i < n; i++) { g[i] = 0.5 * target[i]; if (!isfinite(g[i])) ok = 0; }   /* rhs = -g */
        ok = ok && Sparse_factor_solve(f->factor, g, x) == 0;
        double value = 0.0;
        for (int i = 0; ok && i < n; i++) { if (!isfinite(x[i])) ok = 0; value += -0.5 * target[i] * x[i]; }
        if (ok) *model = 0.5 * value;   /* 1/2 g'x at Hx = -g */
    }
    free(g); free(x);
    return ok;
}

int ActiveSetQp_solve_factored(ActiveSetQpFactor_T f, const double *target, int m, int an,
                               const int *ar, const int *ac, const double *av, const double *upper,
                               double *result, ActiveSetQpStats *stats)
{ return ActiveSetQp_solve_warm(f, target, m, an, ar, ac, av, upper, result, stats, NULL, NULL, 0); }

int ActiveSetQp_solve_warm(ActiveSetQpFactor_T f, const double *target, int m, int an,
                           const int *ar, const int *ac, const double *av, const double *upper,
                           double *result, ActiveSetQpStats *stats, int *rows, int *count, int capacity)
{
    ActiveSetQpStats local; memset(&local, 0, sizeof local);
    if (stats) { memset(stats, 0, sizeof *stats); stats->reason = ASQP_REASON_INPUT; }
    if (!f || !f->factor || m < 0 || an < 0 || !target || !result || (m && !upper) || (an && (!ar || !ac || !av))) return 0;
    const int n = f->n;
    for (int k = 0; k < an; k++) if (ar[k] < 0 || ar[k] >= m || ac[k] < 0 || ac[k] >= n || !isfinite(av[k])) return 0;
    for (int i = 0; i < n; i++) if (!isfinite(target[i])) return 0;
    for (int i = 0; i < m; i++) if (!isfinite(upper[i])) return 0;
    local.factor_seconds = f->factor_seconds; local.ridge_retries = f->ridge_retries;
    int ok = 0, most = m < n ? m : n; if (most > ASQP_MAX_ACTIVE) most = ASQP_MAX_ACTIVE; if (most < 1) most = 1;
    int cap = most < asqp_initial_cap ? most : asqp_initial_cap; if (cap < 1) cap = 1;
    Asqp q; memset(&q, 0, sizeof q); q.n = n; q.m = m; q.cap = cap; q.max_cap = most; q.upper = upper; q.factor = f->factor;
    q.hn = f->hn; q.hr = f->hr; q.hc = f->hc; q.half = f->half; q.value_limit = f->value_limit;
    if (rows && count && capacity > 0) { q.warm = rows; q.nwarm = *count > 0 ? *count : 0; q.warm_cap = capacity; }
    double *g = malloc((size_t)n * sizeof *g);
    int *counts = calloc((size_t)(m ? m : 1), sizeof *counts);
    q.row = malloc((size_t)(m ? m : 1) * sizeof *q.row); q.ids = malloc((size_t)(an ? an : 1) * sizeof *q.ids); q.vals = malloc((size_t)(an ? an : 1) * sizeof *q.vals);
    q.x = malloc((size_t)n * sizeof *q.x); q.yp = malloc((size_t)n * sizeof *q.yp); q.z = malloc((size_t)n * sizeof *q.z); q.rhs = malloc((size_t)n * sizeof *q.rhs);
    q.lambda = calloc((size_t)cap, sizeof *q.lambda); q.s = malloc((size_t)cap * sizeof *q.s); q.r = malloc((size_t)cap * sizeof *q.r);
    q.W = malloc((size_t)cap * sizeof *q.W); q.inW = calloc((size_t)(m ? m : 1), sizeof *q.inW);
    q.M = calloc((size_t)cap * (size_t)cap, sizeof *q.M); q.L = calloc((size_t)cap * (size_t)cap, sizeof *q.L);
    if (!g || !counts || !q.row || !q.ids || !q.vals || !q.x || !q.yp || !q.z || !q.rhs ||
        !q.lambda || !q.s || !q.r || !q.W || !q.inW || !q.M || !q.L) { local.reason = ASQP_REASON_MEMORY; goto done; }
    /* gradient from the doubled convention */
    for (int i = 0; i < n; i++) g[i] = -0.5 * target[i];
    q.g = g;
    /* constraint rows to CSR (entries kept in input order within a row) */
    for (int k = 0; k < an; k++) counts[ar[k]]++;
    { int at = 0; for (int i = 0; i < m; i++) { q.row[i].first = at; q.row[i].count = 0; at += counts[i]; } }
    for (int k = 0; k < an; k++) { AsqpRow *rw = &q.row[ar[k]]; q.ids[rw->first + rw->count] = ac[k]; q.vals[rw->first + rw->count] = av[k]; rw->count++; }
    q.st = local;
    ok = asqp_run(&q, g);
    local = q.st;
    if (ok) memcpy(result, q.x, (size_t)n * sizeof *result);
    if (count && q.warm) *count = ok ? q.nwarm : 0;
done:
    free(g); free(counts); free(q.row); free(q.ids); free(q.vals);
    free(q.x); free(q.yp); free(q.z); free(q.rhs); free(q.lambda); free(q.s); free(q.r); free(q.W); free(q.inW); free(q.M); free(q.L);
    if (stats) *stats = local;
    return ok;
}

int ActiveSetQp_solve_stats(int n, int hn, const int *hr, const int *hc, const double *hv,
                            const double *target, int m, int an, const int *ar, const int *ac,
                            const double *av, const double *upper, double *result,
                            ActiveSetQpStats *stats)
{
    ActiveSetQpFactor_T f = NULL; ActiveSetQpStats fs;
    if (!ActiveSetQp_factor(n, hn, hr, hc, hv, &f, &fs)) { if (stats) *stats = fs; return 0; }
    int ok = ActiveSetQp_solve_factored(f, target, m, an, ar, ac, av, upper, result, stats);
    ActiveSetQp_free(&f);
    return ok;
}

int ActiveSetQp_solve(int n, int hn, const int *hr, const int *hc, const double *hv,
                      const double *target, int m, int an, const int *ar, const int *ac,
                      const double *av, const double *upper, double *result)
{ return ActiveSetQp_solve_stats(n, hn, hr, hc, hv, target, m, an, ar, ac, av, upper, result, NULL); }

const char *ActiveSetQp_reason(int reason)
{
    switch (reason) {
    case ASQP_REASON_OK: return "solved";
    case ASQP_REASON_INPUT: return "invalid input";
    case ASQP_REASON_MEMORY: return "out of memory";
    case ASQP_REASON_FACTOR: return "H would not factor";
    case ASQP_REASON_SOLVE: return "triangular solve failed";
    case ASQP_REASON_INFEASIBLE: return "infeasible rows";
    case ASQP_REASON_STEPS: return "step guard exhausted";
    case ASQP_REASON_CAP: return "working set outgrew its cap";
    case ASQP_REASON_DENSE: return "working-set Cholesky lost definiteness";
    case ASQP_REASON_ILLPOSED: return "objective beyond the caller's limit";
    case ASQP_REASON_TIME: return "outran its wall-clock budget";
    default: return "?";
    }
}

/* ---- selftest -------------------------------------------------------------------- */

static unsigned asqp_rng(unsigned *s) { *s ^= *s << 13; *s ^= *s >> 17; *s ^= *s << 5; return *s; }
static double asqp_unit(unsigned *s) { return (double)(asqp_rng(s) % 1000000) / 1000000.0 - 0.5; }

/* KKT residual of x for min 1/2 x'Hx + g'x, a_i'x <= u_i; multipliers by
 * nonnegative least squares are not needed: check primal feasibility, and
 * stationarity restricted to the binding set through a tiny dense solve. */
static int asqp_kkt_ok(int n, int hn, const int *hr, const int *hc, const double *hv, const double *target,
                       int m, int an, const int *ar, const int *ac, const double *av, const double *upper,
                       const double *x, double tol)
{
    double *Hx = calloc((size_t)n, sizeof *Hx), *g = malloc((size_t)n * sizeof *g), *resid = malloc((size_t)n * sizeof *resid);
    int ok = Hx && g && resid;
    if (ok) {
        double dmax = 0.0, *dg = calloc((size_t)n, sizeof *dg);
        if (!dg) { free(Hx); free(g); free(resid); return 0; }
        for (int k = 0; k < hn; k++) { double v = 0.5 * hv[k]; Hx[hr[k]] += v * x[hc[k]]; if (hr[k] != hc[k]) Hx[hc[k]] += v * x[hr[k]]; else dg[hr[k]] += v; }
        for (int i = 0; i < n; i++) if (dg[i] > dmax) dmax = dg[i];
        /* the solver minimises the RIDGED objective: stationarity holds for H + diag(ridge_i) */
        for (int i = 0; i < n; i++) {
            double ridge = dg[i] > 1e-12 * dmax ? 0.0 : (dmax > 0.0 ? dmax : 1.0) * ASQP_RIDGE;   /* attempt 0: exact on live coordinates */
            g[i] = -0.5 * target[i]; resid[i] = Hx[i] + ridge * x[i] + g[i];
        }
        free(dg);
        /* feasibility */
        for (int i = 0; i < m && ok; i++) {
            double ax = 0; for (int k = 0; k < an; k++) if (ar[k] == i) ax += av[k] * x[ac[k]];
            if (ax > upper[i] + tol * (1 + fabs(upper[i]))) ok = 0;
        }
        /* stationarity: resid must be a nonnegative combination of binding rows.
         * With at most a handful binding, solve the small least squares and
         * check the fit and the sign of every coefficient. */
        int bind[64], nb = 0;
        for (int i = 0; i < m && nb < 64; i++) {
            double ax = 0; for (int k = 0; k < an; k++) if (ar[k] == i) ax += av[k] * x[ac[k]];
            if (fabs(ax - upper[i]) <= 1e-6 * (1 + fabs(upper[i]))) bind[nb++] = i;
        }
        double A[64][64] = {{0}}, b[64] = {0}, lam[64] = {0};
        for (int i = 0; i < nb; i++) for (int j = 0; j < nb; j++) {
            double d = 0; for (int p = 0; p < an; p++) if (ar[p] == bind[i]) for (int qq = 0; qq < an; qq++) if (ar[qq] == bind[j] && ac[qq] == ac[p]) d += av[p] * av[qq];
            A[i][j] = d + (i == j ? 1e-14 : 0);
        }
        for (int i = 0; i < nb; i++) { double d = 0; for (int p = 0; p < an; p++) if (ar[p] == bind[i]) d += av[p] * (-resid[ac[p]]); b[i] = d; }
        /* Gaussian elimination */
        for (int i = 0; i < nb && ok; i++) {
            int piv = i; for (int r = i + 1; r < nb; r++) if (fabs(A[r][i]) > fabs(A[piv][i])) piv = r;
            if (fabs(A[piv][i]) < 1e-18) { ok = 0; break; }
            if (piv != i) { for (int c = 0; c < nb; c++) { double t = A[i][c]; A[i][c] = A[piv][c]; A[piv][c] = t; } double t = b[i]; b[i] = b[piv]; b[piv] = t; }
            for (int r = i + 1; r < nb; r++) { double f = A[r][i] / A[i][i]; for (int c = i; c < nb; c++) A[r][c] -= f * A[i][c]; b[r] -= f * b[i]; }
        }
        for (int i = nb - 1; i >= 0 && ok; i--) { double sum = b[i]; for (int c = i + 1; c < nb; c++) sum -= A[i][c] * lam[c]; lam[i] = sum / A[i][i]; }
        for (int i = 0; i < nb && ok; i++) if (lam[i] < -tol) ok = 0;
        if (ok) {
            double worst = 0;
            for (int v = 0; v < n; v++) {
                double r = resid[v]; for (int i = 0; i < nb; i++) for (int p = 0; p < an; p++) if (ar[p] == bind[i] && ac[p] == v) r += lam[i] * av[p];
                if (fabs(r) > worst) worst = fabs(r);
            }
            if (worst > tol * (1 + fabs(x[0]))) ok = 0;
        }
    }
    free(Hx); free(g); free(resid);
    return ok;
}

int ActiveSetQp_selftest(void)
{
    int fails = 0, seen = 0;
#define ASQP_CASE(name) do { if (fails != seen) { fprintf(stderr, "  active_set_qp: case %s failed (%d)\n", name, fails - seen); seen = fails; } } while (0)
    /* Hand cases are checked to 1e-7: a well-posed H is factored exactly, but
     * a case with a dead coordinate carries that coordinate's ASQP_RIDGE floor. */
    /* 1. one variable, unconstrained: min 1/2*4x^2 - 8x -> x = 2 (stored H = 8, target = 16) */
    { int hr[1] = {0}, hc[1] = {0}; double hv[1] = {8.0}, target[1] = {16.0}, x = 0;
      if (!ActiveSetQp_solve(1, 1, hr, hc, hv, target, 0, 0, NULL, NULL, NULL, NULL, &x) || fabs(x - 2.0) > 1e-7) fails++; }
    ASQP_CASE("1 unconstrained");
    /* 2. same with x <= 1 binding -> x = 1; and with x <= 5 slack -> x = 2 */
    { int hr[1] = {0}, hc[1] = {0}, ar[1] = {0}, ac[1] = {0}; double hv[1] = {8.0}, target[1] = {16.0}, av[1] = {1.0}, up[1] = {1.0}, x = 0;
      if (!ActiveSetQp_solve(1, 1, hr, hc, hv, target, 1, 1, ar, ac, av, up, &x) || fabs(x - 1.0) > 1e-7) fails++;
      up[0] = 5.0; if (!ActiveSetQp_solve(1, 1, hr, hc, hv, target, 1, 1, ar, ac, av, up, &x) || fabs(x - 2.0) > 1e-7) fails++; }
    ASQP_CASE("2 binding/slack bound");
    /* 3. negative upper bound forces a nonzero correction: min 1/2 x^2, -x <= -3 -> x = 3 */
    { int hr[1] = {0}, hc[1] = {0}, ar[1] = {0}, ac[1] = {0}; double hv[1] = {2.0}, target[1] = {0.0}, av[1] = {-1.0}, up[1] = {-3.0}, x = 0;
      if (!ActiveSetQp_solve(1, 1, hr, hc, hv, target, 1, 1, ar, ac, av, up, &x) || fabs(x - 3.0) > 1e-7) fails++; }
    ASQP_CASE("3 negative upper bound");
    /* 4. infeasible pair: x <= 0 and -x <= -1 -> refused */
    { int hr[1] = {0}, hc[1] = {0}, ar[2] = {0, 1}, ac[2] = {0, 0}; double hv[1] = {2.0}, target[1] = {0.0}, av[2] = {1.0, -1.0}, up[2] = {0.0, -1.0}, x = 0;
      if (ActiveSetQp_solve(1, 1, hr, hc, hv, target, 2, 2, ar, ac, av, up, &x)) fails++; }
    ASQP_CASE("4 infeasible pair");
    /* 5. two variables, a dependent duplicate row must not break the working set */
    { int hr[2] = {0, 1}, hc[2] = {0, 1}, ar[4] = {0, 0, 1, 1}, ac[4] = {0, 1, 0, 1};
      double hv[2] = {2.0, 2.0}, target[2] = {4.0, 4.0}, av[4] = {1.0, 1.0, 1.0, 1.0}, up[2] = {1.0, 1.0}, x[2] = {0, 0};
      if (!ActiveSetQp_solve(2, 2, hr, hc, hv, target, 2, 4, ar, ac, av, up, x) || fabs(x[0] - 0.5) > 1e-7 || fabs(x[1] - 0.5) > 1e-7) fails++; }
    ASQP_CASE("5 dependent duplicate row");
    /* 5b. a SINGULAR H (a coordinate with no curvature) must neither crash the
     * factorisation nor return a non-finite step: the ridge retry handles it. */
    { int hr[3] = {0, 1, 2}, hc[3] = {0, 1, 2}; double hv[3] = {4.0, 0.0, 4.0}, target[3] = {8.0, 0.0, -8.0}, x[3] = {0, 0, 0};
      ActiveSetQpStats st;
      int ok = ActiveSetQp_solve_stats(3, 3, hr, hc, hv, target, 0, 0, NULL, NULL, NULL, NULL, x, &st);
      if (!ok || !isfinite(x[0]) || !isfinite(x[1]) || !isfinite(x[2]) || fabs(x[0] - 2.0) > 1e-7 || fabs(x[2] + 2.0) > 1e-7 || st.ridge_retries < 1) fails++;
      /* and the same coordinate absent from H altogether (no explicit zero) */
      int hr2[2] = {0, 2}, hc2[2] = {0, 2}; double hv2[2] = {4.0, 4.0};
      ok = ActiveSetQp_solve_stats(3, 2, hr2, hc2, hv2, target, 0, 0, NULL, NULL, NULL, NULL, x, &st);
      if (!ok || !isfinite(x[1]) || fabs(x[0] - 2.0) > 1e-7) fails++; }
    ASQP_CASE("5b singular H");
    /* 5c. the same deficiency inside a real supernodal tree: a 24x24 grid
     * Laplacian (+I) with three coordinates stripped of every entry, so the
     * multifrontal factorisation meets zero pivots deep in its recursion. */
    {
        enum { G = 24, N = G * G };
        int *hr = malloc((size_t)(5 * N) * sizeof *hr), *hc = malloc((size_t)(5 * N) * sizeof *hc); double *hv = malloc((size_t)(5 * N) * sizeof *hv);
        double *target = malloc((size_t)N * sizeof *target), *x = malloc((size_t)N * sizeof *x);
        int hn = 0, dead[3] = {5, 300, 570};
        for (int i = 0; i < N; i++) {
            int isdead = i == dead[0] || i == dead[1] || i == dead[2];
            if (!isdead) { hr[hn] = i; hc[hn] = i; hv[hn] = 2.0 * 5.0; hn++; }
            int r = i / G, c = i % G;
            int nb[2] = { c > 0 ? i - 1 : -1, r > 0 ? i - G : -1 };
            for (int k = 0; k < 2; k++) if (nb[k] >= 0) {
                int j = nb[k], jdead = j == dead[0] || j == dead[1] || j == dead[2];
                if (isdead || jdead) continue;
                hr[hn] = i; hc[hn] = j; hv[hn] = 2.0 * -1.0; hn++;
            }
            target[i] = 2.0 * ((i * 7919) % 13 - 6);
        }
        ActiveSetQpStats st; int ok = ActiveSetQp_solve_stats(N, hn, hr, hc, hv, target, 0, 0, NULL, NULL, NULL, NULL, x, &st);
        int finite = 1; for (int i = 0; i < N; i++) if (!isfinite(x[i])) finite = 0;
        if (!ok || !finite || st.ridge_retries < 1) fails++;
        free(hr); free(hc); free(hv); free(target); free(x);
    }
    ASQP_CASE("5c supernodal dead coordinates");
    /* 5d. indefinite H with a real supernodal tree: the factorisation must
     * FAIL CLEANLY, twice in a row, and a well-posed factorisation must still
     * succeed afterwards.  TAUCS's multifrontal unwind used to free the
     * parent's frontal matrix twice (ASAN heap-use-after-free, 2026-09-15). */
    {
        enum { G = 24, N = G * G };
        int *hr = malloc((size_t)(3 * N) * sizeof *hr), *hc = malloc((size_t)(3 * N) * sizeof *hc); double *hv = malloc((size_t)(3 * N) * sizeof *hv);
        int hn = 0;
        for (int i = 0; i < N; i++) {
            hr[hn] = i; hc[hn] = i; hv[hn] = (i == 300) ? -10.0 : 5.0; hn++;
            int r = i / G, c = i % G;
            if (c > 0) { hr[hn] = i; hc[hn] = i - 1; hv[hn] = -1.0; hn++; }
            if (r > 0) { hr[hn] = i; hc[hn] = i - G; hv[hn] = -1.0; hn++; }
        }
        SparseFactor_T f = NULL; int clean = 0;
        for (int rep = 0; rep < 2; rep++) { f = NULL; if (Sparse_factor_spd(N, hn, hr, hc, hv, &f) != 0 && f == NULL) clean++; else if (f) Sparse_factor_free(&f); }
        hv[0] = 5.0; for (int k = 0; k < hn; k++) if (hr[k] == 300 && hc[k] == 300) hv[k] = 5.0;
        double *b = malloc((size_t)N * sizeof *b), *x = malloc((size_t)N * sizeof *x);
        for (int i = 0; i < N; i++) b[i] = 1.0;
        int good = Sparse_factor_spd(N, hn, hr, hc, hv, &f) == 0 && Sparse_factor_solve(f, b, x) == 0 && isfinite(x[0]);
        if (f) Sparse_factor_free(&f);
        if (clean != 2 || !good) fails++;
        free(hr); free(hc); free(hv); free(b); free(x);
    }
    ASQP_CASE("5d indefinite clean failure");
    /* 6. random SPD problems: KKT conditions, and agreement with the OSQP path.
     * The working set starts at capacity 2 here so that every growth step
     * (re-striding M and L) is exercised and checked by the KKT residual. */
    asqp_initial_cap = 2;
    unsigned seed = 0x9E3779B9u; int agree_checked = 0, kkt_bad = 0, agree_bad = 0, refused_bad = 0;
    for (int trial = 0; trial < 24; trial++) {
        int n = 6 + (int)(asqp_rng(&seed) % 25), m = (int)(asqp_rng(&seed) % 20);
        int hn = 0, an = 0; int *hr = malloc((size_t)(n*n) * sizeof *hr), *hc = malloc((size_t)(n*n) * sizeof *hc); double *hv = malloc((size_t)(n*n) * sizeof *hv);
        int *ar = malloc((size_t)(m*n + 1) * sizeof *ar), *ac = malloc((size_t)(m*n + 1) * sizeof *ac); double *av = malloc((size_t)(m*n + 1) * sizeof *av);
        double *target = malloc((size_t)n * sizeof *target), *up = malloc((size_t)(m + 1) * sizeof *up), *x = malloc((size_t)n * sizeof *x), *y = malloc((size_t)n * sizeof *y);
        /* H = B'B + n*I built from random B, stored doubled and lower */
        double *B = malloc((size_t)n*n * sizeof *B);
        for (int i = 0; i < n*n; i++) B[i] = asqp_unit(&seed);
        for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++) {
            double d = (i == j) ? (double)n : 0.0; for (int k = 0; k < n; k++) d += B[k*n + i] * B[k*n + j];
            if (fabs(d) > 1e-12 || i == j) { hr[hn] = i; hc[hn] = j; hv[hn] = 2.0 * d; hn++; }
        }
        for (int i = 0; i < n; i++) target[i] = -2.0 * asqp_unit(&seed) * 10.0;
        for (int i = 0; i < m; i++) {
            int nz = 1 + (int)(asqp_rng(&seed) % 3);
            for (int k = 0; k < nz; k++) { ar[an] = i; ac[an] = (int)(asqp_rng(&seed) % n); av[an] = asqp_unit(&seed) * 2.0; an++; }
            up[i] = asqp_unit(&seed) * 2.0 - 0.3;   /* some negative: forced corrections */
        }
        ActiveSetQpStats st; int ok = ActiveSetQp_solve_stats(n, hn, hr, hc, hv, target, m, an, ar, ac, av, up, x, &st);
        int osqp_ok = ConvexQp_solve(n, hn, hr, hc, hv, target, m, an, ar, ac, av, up, y);
        if (ok) {
            if (!asqp_kkt_ok(n, hn, hr, hc, hv, target, m, an, ar, ac, av, up, x, 1e-7)) { fails++; kkt_bad++; }
            if (osqp_ok) {
                double diff = 0, scale = 1; for (int i = 0; i < n; i++) { diff = fmax(diff, fabs(x[i] - y[i])); scale = fmax(scale, fabs(x[i])); }
                if (diff > 1e-4 * scale) { fails++; agree_bad++; }
                agree_checked++;
            }
        } else if (osqp_ok) {
            /* the active set refused what OSQP solved: only acceptable if OSQP's answer is itself infeasible */
            int feasible = 1;
            for (int i = 0; i < m; i++) { double ax = 0; for (int k = 0; k < an; k++) if (ar[k] == i) ax += av[k] * y[ac[k]]; if (ax > up[i] + 1e-5 * (1 + fabs(up[i]))) feasible = 0; }
            if (feasible) { fails++; refused_bad++; }
        }
        free(hr); free(hc); free(hv); free(ar); free(ac); free(av); free(target); free(up); free(x); free(y); free(B);
    }
    asqp_initial_cap = ASQP_INITIAL_ACTIVE;
    ASQP_CASE("6 random");
    /* 6b. THE WALL-CLOCK BUDGET.  The step guard is 8*rows + 1024, which on the 21x5x5's
     * restoration is 856,808 and never fires; the solve grinds to the working-set SIZE cap
     * instead and is refused after 4,773 s.  With the budget shrunk to nothing, a solve that
     * needs any add/drop steps must be refused ON TIME and say so, and the same solve must
     * succeed once the budget is restored -- so the budget is the only thing that changed. */
    {
        enum { N = 160 };   /* over the 64-step check interval, so the budget is really exercised */
        int hr[N], hc[N], ar[N], ac[N];
        double hv[N], target[N], av[N], up[N], x[N];
        int i = 0;
        for (i = 0; i < N; i++) {
            hr[i] = hc[i] = i; hv[i] = 2.0; target[i] = 20.0;   /* min sum (x_i - 10)^2 */
            ar[i] = i; ac[i] = i; av[i] = 1.0; up[i] = -1.0 - (double)i;   /* every row binds */
            x[i] = 0.0;
        }
        ActiveSetQpStats st;
        ActiveSetQp_deadline(ves_clock_sec()-1);
        if(ActiveSetQp_solve_stats(N,N,hr,hc,hv,target,N,N,ar,ac,av,up,x,&st) || st.reason!=ASQP_REASON_TIME)fails++;
        for(i=0;i<N;i++)if(x[i]!=0)fails++;
        ActiveSetQp_deadline(0);
        int solved = ActiveSetQp_solve_stats(N, N, hr, hc, hv, target, N, N, ar, ac, av, up, x, &st);
        int steps = st.iterations;
        if (!solved || steps < 64) { fprintf(stderr, "  active_set_qp: the budget fixture must solve in 64+ steps (solved %d in %d)\n", solved, steps); fails++; }
        asqp_max_seconds = 1e-9;
        for (i = 0; i < N; i++) x[i] = 0.0;
        int refused = !ActiveSetQp_solve_stats(N, N, hr, hc, hv, target, N, N, ar, ac, av, up, x, &st);
        asqp_max_seconds = ASQP_MAX_SECONDS;
        if (!(refused && st.reason == ASQP_REASON_TIME)) {
            fprintf(stderr, "  active_set_qp: a zero budget gave refused %d reason %s\n", refused, ActiveSetQp_reason(st.reason));
            fails++;
        }
        for (i = 0; i < N; i++) x[i] = 0.0;
        if (!ActiveSetQp_solve_stats(N, N, hr, hc, hv, target, N, N, ar, ac, av, up, x, &st)) {
            fprintf(stderr, "  active_set_qp: the restored budget must solve it again (%s)\n", ActiveSetQp_reason(st.reason));
            fails++;
        }
    }
    ASQP_CASE("6b wall-clock budget");
    /* 7. deleting working-set slots by Givens rotations must reproduce a fresh
     * factorisation of the compacted M (the SPD Cholesky factor is unique),
     * across three successive deletes carried on the same factor. */
    { unsigned seed7 = 7u; const int kd = 14; Asqp q; memset(&q, 0, sizeof q); q.cap = kd;
      q.M = calloc((size_t)kd * kd, sizeof *q.M); q.L = calloc((size_t)kd * kd, sizeof *q.L);
      double *B = malloc((size_t)kd * kd * sizeof *B), *Lref = calloc((size_t)kd * kd, sizeof *Lref);
      if (q.M && q.L && B && Lref) {
        for (int i = 0; i < kd * kd; i++) B[i] = asqp_unit(&seed7);
        for (int i = 0; i < kd; i++) for (int j = 0; j < kd; j++) { double s = i == j ? (double)kd : 0.0; for (int t = 0; t < kd; t++) s += B[i*kd + t] * B[j*kd + t]; q.M[(size_t)i*kd + j] = s; }
        const int drops[3] = { 3, 0, kd - 3 };
        int kk = kd, bad = 0;
        if (!asqp_chol_full(&q, kk)) bad++;
        for (int d = 0; d < 3 && !bad; d++) {
          int j = drops[d];
          if (!asqp_chol_delete(&q, j, kk)) { bad++; break; }
          for (int i = 0; i < kk; i++) if (i != j) { int ii = i < j ? i : i - 1; for (int c = 0; c < kk; c++) if (c != j) { int cc = c < j ? c : c - 1; q.M[(size_t)ii*kd + cc] = q.M[(size_t)i*kd + c]; } }
          kk--;
          Asqp r = q; r.L = Lref;
          if (!asqp_chol_full(&r, kk)) { bad++; break; }
          for (int i = 0; i < kk; i++) for (int c = 0; c <= i; c++)
            if (fabs(q.L[(size_t)i*kd + c] - Lref[(size_t)i*kd + c]) > 1e-9 * (1.0 + fabs(Lref[(size_t)i*kd + c]))) bad++;
        }
        if (bad) fails++;
      } else fails++;
      free(q.M); free(q.L); free(B); free(Lref); }
    ASQP_CASE("7 chol delete");
    /* 8. one factor, two row sets: must match the one-shot solves exactly */
    { int hr[3] = {0, 1, 1}, hc[3] = {0, 0, 1}; double hv[3] = {8.0, 2.0, 8.0}, target[2] = {16.0, 4.0};
      int ar[2] = {0, 0}, ac[2] = {0, 1}; double av[2] = {1.0, 1.0}, up1[1] = {1.0}, up2[1] = {0.25};
      double a[2], b[2], c[2], d[2]; ActiveSetQpFactor_T f = NULL; ActiveSetQpStats st;
      int good = ActiveSetQp_factor(2, 3, hr, hc, hv, &f, &st) && f;
      good = good && ActiveSetQp_solve_factored(f, target, 1, 2, ar, ac, av, up1, a, &st) && ActiveSetQp_solve_stats(2, 3, hr, hc, hv, target, 1, 2, ar, ac, av, up1, b, &st);
      good = good && ActiveSetQp_solve_factored(f, target, 1, 2, ar, ac, av, up2, c, &st) && ActiveSetQp_solve_stats(2, 3, hr, hc, hv, target, 1, 2, ar, ac, av, up2, d, &st);
      good = good && ActiveSetQp_solve_factored(f, target, 0, 0, NULL, NULL, NULL, NULL, a, &st);
      ActiveSetQp_free(&f);
      if (!good || f || fabs(a[0] - 2.0) > 1e-9 || fabs(a[1] - 0.0) > 1e-9 || fabs(b[0] + b[1] - 1.0) > 1e-7 || fabs(c[0] + c[1] - 0.25) > 1e-7 || fabs(c[0] - d[0]) > 1e-12 || fabs(c[1] - d[1]) > 1e-12) fails++; }
    /* 9. the unconstrained model value: min 1/2*4x^2 - 8x has value -8 at x = 2 */
    { int hr[1] = {0}, hc[1] = {0}; double hv[1] = {8.0}, target[1] = {16.0}, model = 0; ActiveSetQpFactor_T f = NULL; ActiveSetQpStats st;
      if (!ActiveSetQp_factor(1, 1, hr, hc, hv, &f, &st) || !ActiveSetQp_model(f, target, &model) || fabs(model + 8.0) > 1e-9) fails++;
      ActiveSetQp_free(&f); }
    /* 10. the running-objective limit: the optimum of case 1 has value -8; a limit
     * below it must refuse as ill-posed, a limit above it must not interfere */
    { int hr[1] = {0}, hc[1] = {0}; double hv[1] = {8.0}, target[1] = {16.0}, x = 0; ActiveSetQpFactor_T f = NULL; ActiveSetQpStats st;
      int good = ActiveSetQp_factor(1, 1, hr, hc, hv, &f, &st);
      ActiveSetQp_limit(f, -100.0);
      good = good && !ActiveSetQp_solve_factored(f, target, 0, 0, NULL, NULL, NULL, NULL, &x, &st) && st.reason == ASQP_REASON_ILLPOSED;
      ActiveSetQp_limit(f, 0.0);
      good = good && ActiveSetQp_solve_factored(f, target, 0, 0, NULL, NULL, NULL, NULL, &x, &st) && fabs(x - 2.0) < 1e-9;
      ActiveSetQp_free(&f); if (!good) fails++; }
    /* 11. warm start: min 1/2|x|^2 - 3.1'x with x0 <= 1, x1 <= 1, x2 <= 5 binds rows 0 and 1
     * (optimum (1,1,3)).  A cold solve reports that set; seeded with it the main
     * loop has nothing to do; seeded with a wrong row the answer is unchanged. */
    { int hr[3] = {0, 1, 2}, hc[3] = {0, 1, 2}; double hv[3] = {2.0, 2.0, 2.0}, target[3] = {6.0, 6.0, 6.0};
      int ar[3] = {0, 1, 2}, ac[3] = {0, 1, 2}; double av[3] = {1.0, 1.0, 1.0}, up[3] = {1.0, 1.0, 5.0};
      double x[3]; int rows[4] = {0, 0, 0, 0}, cnt = 0; ActiveSetQpFactor_T f = NULL; ActiveSetQpStats st;
      int good = ActiveSetQp_factor(3, 3, hr, hc, hv, &f, &st);
      good = good && ActiveSetQp_solve_warm(f, target, 3, 3, ar, ac, av, up, x, &st, rows, &cnt, 4) && cnt == 2
             && fabs(x[0] - 1.0) < 1e-9 && fabs(x[1] - 1.0) < 1e-9 && fabs(x[2] - 3.0) < 1e-9 && st.warm == 0;
      good = good && ActiveSetQp_solve_warm(f, target, 3, 3, ar, ac, av, up, x, &st, rows, &cnt, 4) && cnt == 2
             && fabs(x[0] - 1.0) < 1e-9 && fabs(x[1] - 1.0) < 1e-9 && fabs(x[2] - 3.0) < 1e-9 && st.warm == 2 && st.iterations == 0;
      rows[0] = 2; rows[1] = 0; cnt = 2;   /* row 2 does not bind: its multiplier is negative and it is dropped */
      good = good && ActiveSetQp_solve_warm(f, target, 3, 3, ar, ac, av, up, x, &st, rows, &cnt, 4) && cnt == 2
             && fabs(x[0] - 1.0) < 1e-9 && fabs(x[1] - 1.0) < 1e-9 && fabs(x[2] - 3.0) < 1e-9 && st.dropped >= 1;
      ActiveSetQp_free(&f); if (!good) fails++; }
    ASQP_CASE("11 warm start");
    fprintf(stderr, "  active_set_qp selftest: %s (%d failures; OSQP agreement checked on %d random problems; random: kkt %d, agreement %d, refused-but-feasible %d)\n",
            fails ? "FAIL" : "ok", fails, agree_checked, kkt_bad, agree_bad, refused_bad);
#undef ASQP_CASE
    return fails;
}
