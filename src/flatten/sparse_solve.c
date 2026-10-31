#include "sparse_solve.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "../common/ves_omp.h"

#include "taucs.h"

/* (2026-10-28) Every allocation of this file from the calling thread's private heap, as TAUCS's own blocks
 * (taucs_malloc.c under SPARSE_TAUCS_THREAD_SAFE; the C library's functions elsewhere).  On the process heap, a
 * 32-thread assembly's workers spent 57% of their samples queued on its lock in these temporaries (the CCS build, the
 * dense path, AMD's workspace: SuiteSparse_config.c).  Allocation only: no value changes.  taucs.h
 * prototypes the stubs only for TAUCS's own sources (TAUCS_CORE): without these an implicit int return truncates them. */
void *taucs_malloc_stub(size_t size);
void *taucs_calloc_stub(size_t nmemb, size_t size);
void *taucs_realloc_stub(void *ptr, size_t size);
void  taucs_free_stub(void *ptr);
#undef malloc
#undef calloc
#undef realloc
#undef free
#define malloc(n) taucs_malloc_stub(n)
#define calloc(n, s) taucs_calloc_stub((n), (s))
#define realloc(p, n) taucs_realloc_stub((p), (n))
#define free(p) taucs_free_stub(p)

/* taucs.h does not prototype these high-level entry points; declare them. */
extern taucs_ccs_matrix *taucs_ccs_create(int m, int n, int nnz, int flags);
extern void              taucs_ccs_free(taucs_ccs_matrix *m);
extern int  taucs_linsolve(taucs_ccs_matrix *A, void **F, int nrhs,
                           void *X, void *B, char *options[], void *opt_arg[]);
/* Low-level entry points of the factor-once path (taucs_private.h). */
extern void  taucs_ccs_order(taucs_ccs_matrix *A, int **perm, int **invperm, char *which);
extern taucs_ccs_matrix *taucs_ccs_permute_symmetrically(taucs_ccs_matrix *A, int *perm, int *invperm);
extern void *taucs_ccs_factor_llt_symbolic(taucs_ccs_matrix *A);
extern int   taucs_ccs_factor_llt_numeric(taucs_ccs_matrix *A, void *L);
extern int   taucs_supernodal_solve_llt(void *L, void *x, void *b);
extern void  taucs_supernodal_factor_free(void *L);
extern void  taucs_supernodal_factor_free_numeric(void *L);
extern void  taucs_vec_permute(int n, int flags, void *v, void *pv, int p[]);
extern void  taucs_vec_ipermute(int n, int flags, void *v, void *pv, int p[]);
extern void *taucs_ccs_factor_llt_mf(taucs_ccs_matrix *A);

/* 1 = the fill-reducing ordering past 64 unknowns is AMD (2026-10-22:
 * stack samples put GENMMD's mmdupd_/mmdelm_ at ~5% of the
 * running threads, on band systems solved once each); 0 = TAUCS's GENMMD.
 * A build choice (2026-10-25): 1 only where the build compiles the vendored AMD
 * (sheet_assemble.vcxproj defines it); every other tool orders with GENMMD,
 * as it did before, and links no AMD. */
#ifndef SPARSE_ORDER_AMD
#define SPARSE_ORDER_AMD 0
#endif

/* SuiteSparse AMD, vendored with OSQP (deps/osqp/amd) and compiled into the
 * tools that link OSQP.  Declared here rather than through amd.h, which pulls
 * in osqp.h: int indices and double controls, as this build configures OSQP
 * (neither OSQP_USE_LONG nor OSQP_USE_FLOAT). */
#if defined(OSQP_USE_LONG) || defined(OSQP_USE_FLOAT)
#error "sparse_solve.c declares the vendored AMD with int indices and double controls"
#endif
#if SPARSE_ORDER_AMD
extern int  amd_order(int n, const int Ap[], const int Ai[], int P[], double Control[], double Info[]);
extern void amd_defaults(double Control[]);
#endif
#define SP_AMD_CONTROL 5
#define SP_AMD_INFO 20
#define SP_AMD_STATUS 0
#define SP_AMD_LNZ 9
#define SP_AMD_OK 0
#define SP_AMD_OK_BUT_JUMBLED 1

/* Orderings come from TAUCS's allocator and return to it (taucs_malloc.c may
 * give each thread its own heap). */
extern void *taucs_malloc_stub(size_t size);
extern void  taucs_free_stub(void *ptr);

/* GENMMD and other legacy TAUCS helpers keep mutable static state (see
 * quad_ribbon_fit.c), and CLAPACK's f2c BLAS keeps every local static, so by
 * default every TAUCS call runs under one lock. Callers that hold
 * qrf_taucs_transaction nest this lock inside theirs; nothing takes that lock
 * or the pattern cache's lock while holding this one. A build that defines
 * SPARSE_TAUCS_THREAD_SAFE (sheet_assemble) links deps/clapack/TS, the BLAS and
 * LAPACK routines TAUCS reaches with their locals automatic, and TAUCS's own
 * per-thread statics (taucs_sn_llt.c, genmmd.c), and calls TAUCS unlocked:
 * parallel repair corrections on PHerc0139 21^3 spent two thirds of their thread time
 * waiting on this lock (2026-09-26). */
static taucs_ccs_matrix *sp_ccs_create(int n, int nnz)
{
    taucs_ccs_matrix *A = NULL;
#ifndef SPARSE_TAUCS_THREAD_SAFE
#pragma omp critical(sparse_taucs_call)
#endif
    A = taucs_ccs_create(n, n, nnz, TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_LOWER);
    return A;
}

static void sp_ccs_free(taucs_ccs_matrix *A)
{
#ifndef SPARSE_TAUCS_THREAD_SAFE
#pragma omp critical(sparse_taucs_call)
#endif
    taucs_ccs_free(A);
}

static void sp_factor_free(void *L, int numeric_only)
{
#ifndef SPARSE_TAUCS_THREAD_SAFE
#pragma omp critical(sparse_taucs_call)
#endif
    {
        if (numeric_only) taucs_supernodal_factor_free_numeric(L);
        else taucs_supernodal_factor_free(L);
    }
}

/* qsort comparator for (row,val) pairs packed as parallel arrays is awkward;
 * sort an index permutation of a column segment by row instead. */
typedef struct { int row; double val; } RV;
static int cmp_rv(const void *a, const void *b)
{
    int ra = ((const RV *)a)->row, rb = ((const RV *)b)->row;
    return (ra < rb) ? -1 : (ra > rb ? 1 : 0);
}

/* Build a TAUCS lower-symmetric CCS matrix from COO triplets (row >= col,
 * duplicates summed).  Returns NULL on invalid input or allocation failure.
 * (The per-column qsort build: the selftest's reference since 2026-10-22.) */
static taucs_ccs_matrix *sp_build_ccs_qsort(int n, int nt,
                                            const int *rows, const int *cols,
                                            const double *vals)
{
    /* --- Bucket triplets by column. --- */
    int *colcnt = (int *)calloc((size_t)n + 1, sizeof(int));
    if (!colcnt) return NULL;
    for (int t = 0; t < nt; t++) {
        int c = cols[t];
        if (c < 0 || c >= n || rows[t] < c || rows[t] >= n) { free(colcnt); return NULL; }
        colcnt[c]++;
    }
    int *cstart = (int *)malloc(((size_t)n + 1) * sizeof(int));
    if (!cstart) { free(colcnt); return NULL; }
    cstart[0] = 0;
    for (int c = 0; c < n; c++) cstart[c + 1] = cstart[c] + colcnt[c];

    RV *buf = (RV *)malloc((size_t)(nt > 0 ? nt : 1) * sizeof(RV));
    int *cur = (int *)malloc((size_t)n * sizeof(int));
    if (!buf || !cur) { free(colcnt); free(cstart); free(buf); free(cur); return NULL; }
    for (int c = 0; c < n; c++) cur[c] = cstart[c];
    for (int t = 0; t < nt; t++) {
        int c = cols[t];
        buf[cur[c]].row = rows[t];
        buf[cur[c]].val = vals[t];
        cur[c]++;
    }

    /* --- Per column: sort by row, merge duplicates -> CCS. --- */
    taucs_ccs_matrix *A = sp_ccs_create(n, (nt > 0 ? nt : 1));
    if (!A) { free(colcnt); free(cstart); free(buf); free(cur); return NULL; }

    int k = 0;
    for (int c = 0; c < n; c++) {
        A->colptr[c] = k;
        int s = cstart[c], e = cstart[c + 1];
        if (e > s) {
            qsort(buf + s, (size_t)(e - s), sizeof(RV), cmp_rv);
            int i = s;
            while (i < e) {
                int r = buf[i].row;
                double acc = 0.0;
                while (i < e && buf[i].row == r) { acc += buf[i].val; i++; }
                A->rowind[k] = r;
                A->values.d[k] = acc;
                k++;
            }
        }
    }
    A->colptr[n] = k;

    free(colcnt); free(cstart); free(buf); free(cur);
    return A;
}

/* Build a TAUCS lower-symmetric CCS matrix from COO triplets (row >= col,
 * duplicates summed in their input order).  Two counting sorts -- by row, then
 * a stable pass by column, so each column's entries arrive in row order -- and
 * one linear merge of equal rows: O(nt + n), no comparisons (2026-10-22: the
 * per-column qsort was 14% of the running threads, on the band
 * systems' duplicate-heavy columns).  Returns NULL on invalid input or
 * allocation failure. */
static taucs_ccs_matrix *sp_build_ccs(int n, int nt,
                                      const int *rows, const int *cols,
                                      const double *vals)
{
    int *rcnt = (int *)calloc((size_t)n + 1, sizeof(int));
    int *ccnt = (int *)calloc((size_t)n + 1, sizeof(int));
    int *by_row = (int *)malloc((size_t)(nt > 0 ? nt : 1) * sizeof(int));
    RV *buf = (RV *)malloc((size_t)(nt > 0 ? nt : 1) * sizeof(RV));
    taucs_ccs_matrix *A = NULL;
    if (!rcnt || !ccnt || !by_row || !buf) goto out;
    for (int t = 0; t < nt; t++) {
        int c = cols[t], r = rows[t];
        if (c < 0 || c >= n || r < c || r >= n) goto out;
        rcnt[r + 1]++; ccnt[c + 1]++;
    }
    for (int i = 0; i < n; i++) { rcnt[i + 1] += rcnt[i]; ccnt[i + 1] += ccnt[i]; }
    for (int t = 0; t < nt; t++) by_row[rcnt[rows[t]]++] = t;             /* rcnt[r] ends at row r's end */
    for (int q = 0; q < nt; q++) {                                          /* ccnt[c] ends at column c's end */
        int t = by_row[q], c = cols[t];
        buf[ccnt[c]].row = rows[t]; buf[ccnt[c]].val = vals[t]; ccnt[c]++;
    }
    A = sp_ccs_create(n, (nt > 0 ? nt : 1));
    if (!A) goto out;
    int k = 0, s0 = 0;
    for (int c = 0; c < n; c++) {
        int e = ccnt[c];                                                    /* column c: buf[s0, e) */
        A->colptr[c] = k;
        for (int i = s0; i < e;) {
            int r = buf[i].row; double acc = 0.0;
            while (i < e && buf[i].row == r) { acc += buf[i].val; i++; }
            A->rowind[k] = r; A->values.d[k] = acc; k++;
        }
        s0 = e;
    }
    A->colptr[n] = k;
out:
    free(rcnt); free(ccnt); free(by_row); free(buf);
    return A;
}

#if SPARSE_ORDER_AMD
/* (2026-10-22) The block graph of A for AMD: when n = b m and, for every
 * block k, its b columns reach the same row blocks outside block k (a vector
 * field's systems: the growth's band equations have 3 x 3 blocks), the m-node
 * lower-triangle graph of those blocks in (*colptr, *rowind) (malloc'd);
 * else 0 and nothing allocated.  Within-block entries are dropped (AMD ignores
 * the diagonal).  1 = built, -1 = allocation failure. */
static int sp_block_graph(const taucs_ccs_matrix *A, int b, int **colptr, int **rowind)
{
    const int n = A->n, m = n / b;
    *colptr = NULL; *rowind = NULL;
    if (b < 2 || n % b || m < 2) return 0;
    int *cp = (int *)malloc(((size_t)m + 1) * sizeof(int));
    int *ri = (int *)malloc((size_t)(A->colptr[n] > 0 ? A->colptr[n] : 1) * sizeof(int));
    if (!cp || !ri) { free(cp); free(ri); return -1; }
    int k_out = 0;
    for (int k = 0; k < m; k++) {
        cp[k] = k_out;
        int first = 1, s0 = k_out, e0 = k_out;  /* block k's first column's distinct row blocks: ri[s0, e0) */
        for (int c = b * k; c < b * k + b; c++) {
            int last = -1, cnt = 0, pos = s0;
            for (int p = A->colptr[c]; p < A->colptr[c + 1]; p++) {
                int rb = A->rowind[p] / b;
                if (rb == k || rb == last) continue;          /* rows sorted: the row blocks arrive sorted */
                last = rb;
                if (first) ri[k_out++] = rb;
                else if (pos >= e0 || ri[pos++] != rb) { free(cp); free(ri); return 0; }
                cnt++;
            }
            if (first) { e0 = k_out; first = 0; }
            else if (cnt != e0 - s0) { free(cp); free(ri); return 0; }
        }
    }
    cp[m] = k_out;
    *colptr = cp; *rowind = ri;
    return 1;
}
#endif

/* The symmetric ordering of A (lower triangle, its diagonal included) as
 * TAUCS takes it: perm[k] = the old index of the k-th pivot, invperm its
 * inverse, both from TAUCS's allocator.  n <= 64: the identity (TAUCS's f2c
 * GENMMD could fault on tiny graphs, and no ordering reduces fill materially at
 * that scale).  Else AMD (SPARSE_ORDER_AMD), falling back to GENMMD where AMD
 * fails.  *lnz (nullable): AMD's count of L's off-diagonal nonzeros, -1 when
 * AMD did not order.  0 ok. */
static int sp_order(taucs_ccs_matrix *A, int use_amd, int **perm, int **invperm, double *lnz)
{
    const int n = A->n;
    *perm = NULL; *invperm = NULL;
    if (lnz) *lnz = -1;
    if (n <= 64 || use_amd) {
        *perm = (int *)taucs_malloc_stub((size_t)n * sizeof(int)); *invperm = (int *)taucs_malloc_stub((size_t)n * sizeof(int));
        if (!*perm || !*invperm) { taucs_free_stub(*perm); taucs_free_stub(*invperm); *perm = *invperm = NULL; return -1; }
        if (n <= 64) {
            for (int i = 0; i < n; i++) { (*perm)[i] = i; (*invperm)[i] = i; }
            return 0;
        }
#if SPARSE_ORDER_AMD
        double control[SP_AMD_CONTROL], info[SP_AMD_INFO];
        amd_defaults(control);
        int st = -1;
        for (int b = 3; b >= 2 && st < 0; b--) {      /* (2026-10-22) the block graph where A has one */
            int *bcp = NULL, *bri = NULL, got = sp_block_graph(A, b, &bcp, &bri), m = n / b;
            if (got < 0) break;
            if (!got) continue;
            int *bperm = (int *)malloc((size_t)m * sizeof(int));
            if (bperm) {
                int bst = amd_order(m, bcp, bri, bperm, control, info);
                if (bst == SP_AMD_OK || bst == SP_AMD_OK_BUT_JUMBLED) {
                    for (int k = 0; k < m; k++) for (int a = 0; a < b; a++) (*perm)[b * k + a] = b * bperm[k] + a;
                    st = SP_AMD_OK;
                    if (lnz) *lnz = -2;                 /* (the block graph's count is not L's) */
                }
            }
            free(bperm); free(bcp); free(bri);
        }
        if (st < 0) st = amd_order(n, A->colptr, A->rowind, *perm, control, info);
        if (st == SP_AMD_OK || st == SP_AMD_OK_BUT_JUMBLED) {
            for (int k = 0; k < n; k++) (*invperm)[(*perm)[k]] = k;
            if (lnz && *lnz != -2) *lnz = info[SP_AMD_LNZ];
            return 0;
        }
#endif
        taucs_free_stub(*perm); taucs_free_stub(*invperm); *perm = *invperm = NULL;   /* AMD failed (or not built): GENMMD */
    }
#ifndef SPARSE_TAUCS_THREAD_SAFE
#pragma omp critical(sparse_taucs_call)
#endif
    taucs_ccs_order(A, perm, invperm, "genmmd");
    if (!*perm || !*invperm) { taucs_free_stub(*perm); taucs_free_stub(*invperm); *perm = *invperm = NULL; return -1; }
    return 0;
}

/* One SPD solve without taucs_linsolve: the ordering, ONE symmetric
 * permutation, the multifrontal supernodal factor, the triangular solves (the
 * path taucs_linsolve takes for LLT + mf, without its per-call option parsing
 * and handle).  0 ok. */
static int sp_solve_spd_once(taucs_ccs_matrix *A, int use_amd, const double *b, double *x)
{
    const int n = A->n;
    int *perm = NULL, *invperm = NULL, rc = -1;
    if (sp_order(A, use_amd, &perm, &invperm, NULL)) return -1;
    double *pb = (double *)malloc((size_t)n * sizeof(double)), *px = (double *)malloc((size_t)n * sizeof(double));
    if (pb && px) {
#ifndef SPARSE_TAUCS_THREAD_SAFE
#pragma omp critical(sparse_taucs_call)
#endif
        {
            taucs_ccs_matrix *PAPT = taucs_ccs_permute_symmetrically(A, perm, invperm);
            void *L = PAPT ? taucs_ccs_factor_llt_mf(PAPT) : NULL;
            if (L) {
                taucs_vec_permute(n, TAUCS_DOUBLE, (void *)b, pb, perm);
                if (taucs_supernodal_solve_llt(L, px, pb) == 0) {
                    taucs_vec_ipermute(n, TAUCS_DOUBLE, px, x, perm);
                    rc = 0;
                }
                taucs_supernodal_factor_free(L);
            }
            if (PAPT) taucs_ccs_free(PAPT);
        }
    }
    free(pb); free(px);
    taucs_free_stub(perm); taucs_free_stub(invperm);
    return rc;
}

/* (2026-10-22) SPD systems of at most SPARSE_DENSE_MAX unknowns by a dense LL^T in double (the growth's band systems split
 * into ~10^7 blocks a round on the PHerc0358 10^3, most of a few nodes, where TAUCS's per-call ordering, symbolic analysis
 * and allocations cost ~1 ms each).  Duplicates summed; -1 when not positive definite (TAUCS's path then decides). */
#ifndef SPARSE_DENSE_MAX
#define SPARSE_DENSE_MAX 96
#endif
int Sparse_spd_dense_limit(void)
{ return SPARSE_DENSE_MAX > 0 ? SPARSE_DENSE_MAX : 0; }

static int sp_dense_spd(int n, int nt, const int *rows, const int *cols, const double *vals, const double *b, double *x)
{
    double *L = (double *)calloc((size_t)n * (size_t)n, sizeof(double)), *y = (double *)malloc((size_t)n * sizeof(double));
    int rc = -1;
    if (!L || !y) goto out;
    for (int t = 0; t < nt; t++) {
        int r = rows[t], c = cols[t];
        if (c < 0 || r < c || r >= n) goto out;
        L[(size_t)r * (size_t)n + (size_t)c] += vals[t];          /* the lower triangle, row-major */
    }
    for (int j = 0; j < n; j++) {
        double *lj = L + (size_t)j * (size_t)n, d = lj[j];
        for (int k = 0; k < j; k++) d -= lj[k] * lj[k];
        if (!(d > 0.0)) goto out;
        lj[j] = sqrt(d);
        for (int i = j + 1; i < n; i++) {
            double *li = L + (size_t)i * (size_t)n, v = li[j];
            for (int k = 0; k < j; k++) v -= li[k] * lj[k];
            li[j] = v / lj[j];
        }
    }
    for (int i = 0; i < n; i++) {                                   /* L y = b */
        const double *li = L + (size_t)i * (size_t)n; double v = b[i];
        for (int k = 0; k < i; k++) v -= li[k] * y[k];
        y[i] = v / li[i];
    }
    for (int i = n - 1; i >= 0; i--) {                              /* L^T x = y */
        double v = y[i];
        for (int k = i + 1; k < n; k++) v -= L[(size_t)k * (size_t)n + (size_t)i] * x[k];
        x[i] = v / L[(size_t)i * (size_t)n + (size_t)i];
    }
    rc = 0;
out:
    free(L); free(y);
    return rc;
}

int Sparse_solve_sym(int n, int nt,
                     const int *rows, const int *cols, const double *vals,
                     const double *b, double *x, SparseMode mode)
{
    if (n <= 0 || nt < 0) return -1;
    if (!rows || !cols || !vals || !b || !x) return -1;

    /* Small active sets legitimately reduce to one column.  The exact scalar
     * solve is simpler and avoids sending a degenerate tree through the
     * multifrontal machinery. */
    if (n == 1) {
        double a = 0.0;
        for (int t = 0; t < nt; t++) {
            if (rows[t] != 0 || cols[t] != 0) return -1;
            a += vals[t];
        }
        if (a == 0.0 || (mode == SPARSE_SPD && a < 0.0)) return -1;
        x[0] = b[0] / a;
        return 0;
    }

    if (mode == SPARSE_SPD && n <= SPARSE_DENSE_MAX) {       /* (2026-10-22) small: dense; else (not SPD) TAUCS decides */
        if (sp_dense_spd(n, nt, rows, cols, vals, b, x) == 0) return 0;
    }
    taucs_ccs_matrix *A = sp_build_ccs(n, nt, rows, cols, vals);
    if (!A) return -1;
    if (mode == SPARSE_SPD) {                 /* (2026-10-22) AMD, one permutation, no option parsing */
        int rc_spd = sp_solve_spd_once(A, SPARSE_ORDER_AMD, b, x);
        sp_ccs_free(A);
        return rc_spd;
    }

    /* --- Solve (the indefinite path: MINRES through taucs_linsolve). --- */
    char *opt_spd[]   = { "taucs.factor.LLT=true", "taucs.factor.mf=true", NULL };
    /* TAUCS's legacy f2c GENMMD path is not reliable on very small graphs
     * (the 13-vertex welded-scroll disk reaches it with only a handful of
     * interior unknowns and can fault under /O2).  Ordering cannot materially
     * reduce fill at this scale, so use the exact identity permutation. */
    char *opt_spd_small[] = { "taucs.factor.LLT=true", "taucs.factor.mf=true",
                              "taucs.factor.ordering=identity", NULL };
    char *opt_indef[] = { "taucs.solve.minres=true",
                          "taucs.solve.maxits=10000",
                          "taucs.solve.convergetol=1e-11", NULL };
    void *F = NULL;
    char **solve_opts = (mode == SPARSE_SPD)
                      ? ((n <= 64) ? opt_spd_small : opt_spd)
                      : opt_indef;
    int rc = -1;
#ifndef SPARSE_TAUCS_THREAD_SAFE
#pragma omp critical(sparse_taucs_call)
#endif
    {
        rc = taucs_linsolve(A, &F, 1, x, (void *)b, solve_opts, NULL);
        if (F) { void *Fp = F; taucs_linsolve(NULL, &Fp, 0, NULL, NULL, NULL, NULL); }
    }
    sp_ccs_free(A);
    return (rc == TAUCS_SUCCESS) ? 0 : -1;
}

/* ---- factor-once / solve-many ------------------------------------------- */

/* Factorising in TAUCS is an ordering, a symbolic elimination with supernode
 * amalgamation, and the numeric factorisation.  On the repair's 1-2M-unknown
 * fields the first two were a quarter of the stage's stack samples
 * (2026-09-15), and they depend only on the sparsity pattern, which the
 * repair reuses across every solve of a trial and across trials.  A small
 * pattern cache keeps the permutation and the symbolic supernodal structure;
 * a hit pays the numeric factorisation only.  A slot lent to a live factor is
 * not shared: a second caller with the same pattern gets a fresh symbolic
 * factor of its own.  Cache access is guarded for OpenMP callers. */
#define SP_CACHE_SLOTS 4
typedef struct SpPattern {
    int n, nnz;
    int *colptr, *rowind;          /* the pattern this slot serves */
    int *perm, *invperm;           /* its ordering */
    void *L;                       /* symbolic supernodal factor; numeric parts freed between uses */
    int busy;                      /* lent to a live SparseFactor */
    unsigned long long stamp;      /* last use, for eviction */
} SpPattern;
static SpPattern sp_cache[SP_CACHE_SLOTS];
static unsigned long long sp_stamp = 0, sp_hits = 0, sp_misses = 0;

static void sp_pattern_clear(SpPattern *p)
{
    if (p->L) sp_factor_free(p->L, 0);
    free(p->colptr); free(p->rowind); taucs_free_stub(p->perm); taucs_free_stub(p->invperm);
    memset(p, 0, sizeof *p);
}

static int sp_pattern_matches(const SpPattern *p, const taucs_ccs_matrix *A)
{
    return p->L && p->n == A->n && p->nnz == A->colptr[A->n] &&
           memcmp(p->colptr, A->colptr, (size_t)(A->n + 1) * sizeof(int)) == 0 &&
           memcmp(p->rowind, A->rowind, (size_t)p->nnz * sizeof(int)) == 0;
}

struct SparseFactor {
    int n;
    void *L;                /* supernodal LLT: borrowed from `slot`, or owned */
    int *perm, *invperm;    /* ordering: borrowed from `slot`, or owned */
    SpPattern *slot;        /* the cache slot lent to this factor, or NULL */
    double *pb, *px;        /* permuted right-hand side and solution */
    double scalar;          /* n==1 exact path: L==NULL */
    taucs_ccs_matrix *A;    /* Sparse_factor_spd_ccs: the owned pattern refactored in place, else NULL */
};

/* Ordering + symbolic factor for a pattern nobody has cached. */
static int sp_symbolic(taucs_ccs_matrix *A, int **perm, int **invperm, void **L)
{
    *perm = NULL; *invperm = NULL; *L = NULL;
    if (sp_order(A, SPARSE_ORDER_AMD, perm, invperm, NULL)) return -1;   /* identity <= 64, else AMD (2026-10-22) */
    taucs_ccs_matrix *PAPT = NULL;
#ifndef SPARSE_TAUCS_THREAD_SAFE
#pragma omp critical(sparse_taucs_call)
#endif
    {
        PAPT = taucs_ccs_permute_symmetrically(A, *perm, *invperm);
        if (PAPT) { *L = taucs_ccs_factor_llt_symbolic(PAPT); taucs_ccs_free(PAPT); }
    }
    if (!PAPT) { taucs_free_stub(*perm); taucs_free_stub(*invperm); *perm = *invperm = NULL; return -1; }
    if (!*L) { taucs_free_stub(*perm); taucs_free_stub(*invperm); *perm = *invperm = NULL; return -1; }
    return 0;
}

/* Numeric factorisation of A's values into the symbolic structure L, which
 * was built for A's pattern under (perm, invperm). */
static int sp_numeric(taucs_ccs_matrix *A, int *perm, int *invperm, void *L)
{
    taucs_ccs_matrix *PAPT = NULL;int rc = -1;
#ifndef SPARSE_TAUCS_THREAD_SAFE
#pragma omp critical(sparse_taucs_call)
#endif
    {
        PAPT = taucs_ccs_permute_symmetrically(A, perm, invperm);
        if (PAPT) { rc = taucs_ccs_factor_llt_numeric(PAPT, L); taucs_ccs_free(PAPT); }
    }
    if (!PAPT) return -1;
    if (rc) { sp_factor_free(L, 1); return -1; }   /* idempotent: TAUCS nulls what it frees */
    return 0;
}

int Sparse_factor_spd(int n, int nt,
                      const int *rows, const int *cols, const double *vals,
                      SparseFactor_T *out)
{
    if (out == NULL) return -1;
    *out = NULL;
    if (n <= 0 || nt < 0 || !rows || !cols || !vals) return -1;

    struct SparseFactor *f = (struct SparseFactor *)calloc(1, sizeof *f);
    if (!f) return -1;
    f->n = n;

    if (n == 1) {
        double a = 0.0;
        for (int t = 0; t < nt; t++) {
            if (rows[t] != 0 || cols[t] != 0) { free(f); return -1; }
            a += vals[t];
        }
        if (a <= 0.0) { free(f); return -1; }
        f->scalar = a;
        *out = f;
        return 0;
    }

    taucs_ccs_matrix *A = sp_build_ccs(n, nt, rows, cols, vals);
    if (!A) {
        fprintf(stderr, "Sparse_factor_spd: CCS build failed (n=%d nt=%d; likely allocation)\n", n, nt);
        free(f);
        return -1;
    }
    f->pb = (double *)malloc((size_t)n * sizeof(double));
    f->px = (double *)malloc((size_t)n * sizeof(double));
    if (!f->pb || !f->px) { free(f->pb); free(f->px); sp_ccs_free(A); free(f); return -1; }

    /* a cached pattern? */
    SpPattern *slot = NULL;
#pragma omp critical(sparse_pattern_cache)
    {
        for (int i = 0; i < SP_CACHE_SLOTS; i++)
            if (!sp_cache[i].busy && sp_pattern_matches(&sp_cache[i], A)) { slot = &sp_cache[i]; break; }
        if (slot) { slot->busy = 1; slot->stamp = ++sp_stamp; sp_hits++; } else sp_misses++;
    }
    int rc = -1;
    if (slot) {
        rc = sp_numeric(A, slot->perm, slot->invperm, slot->L);
        if (rc) {
#pragma omp critical(sparse_pattern_cache)
            slot->busy = 0;
        } else { f->L = slot->L; f->perm = slot->perm; f->invperm = slot->invperm; f->slot = slot; }
    } else {
        int *perm = NULL, *invperm = NULL; void *L = NULL;
        rc = sp_symbolic(A, &perm, &invperm, &L);
        if (!rc) rc = sp_numeric(A, perm, invperm, L);
        if (rc) { if (L) sp_factor_free(L, 0); taucs_free_stub(perm); taucs_free_stub(invperm); }
        else {
            f->L = L; f->perm = perm; f->invperm = invperm;
            /* offer the pattern to the cache: an empty slot, else the oldest idle one */
            int *colptr = (int *)malloc((size_t)(n + 1) * sizeof(int));
            int *rowind = (int *)malloc((size_t)(A->colptr[n] > 0 ? A->colptr[n] : 1) * sizeof(int));
            if (colptr && rowind) {
                memcpy(colptr, A->colptr, (size_t)(n + 1) * sizeof(int));
                memcpy(rowind, A->rowind, (size_t)A->colptr[n] * sizeof(int));
#pragma omp critical(sparse_pattern_cache)
                {
                    SpPattern *victim = NULL;
                    for (int i = 0; i < SP_CACHE_SLOTS; i++) {
                        if (sp_cache[i].busy) continue;
                        if (!sp_cache[i].L) { victim = &sp_cache[i]; break; }
                        if (!victim || sp_cache[i].stamp < victim->stamp) victim = &sp_cache[i];
                    }
                    if (victim) {
                        sp_pattern_clear(victim);
                        victim->n = n; victim->nnz = A->colptr[n]; victim->colptr = colptr; victim->rowind = rowind;
                        victim->perm = perm; victim->invperm = invperm; victim->L = L;
                        victim->busy = 1; victim->stamp = ++sp_stamp;
                        f->slot = victim; colptr = rowind = NULL;
                    }
                }
            }
            free(colptr); free(rowind);   /* not installed: the factor owns L and the ordering */
        }
    }
    sp_ccs_free(A);
    if (rc) {
        fprintf(stderr, "Sparse_factor_spd: factorization failed (n=%d nt=%d)\n", n, nt);
        free(f->pb); free(f->px); free(f);
        return -1;
    }
    *out = f;
    return 0;
}

int Sparse_factor_solve_multi(SparseFactor_T f, const double *B, double *X,
                              int nrhs)
{
    if (!f || !B || !X || nrhs <= 0) return -1;
    if (f->L == NULL) {           /* scalar path */
        for (int j = 0; j < nrhs; j++) X[j] = B[j] / f->scalar;
        return 0;
    }
    const int n = f->n;
    int rc = 0;
#ifndef SPARSE_TAUCS_THREAD_SAFE
#pragma omp critical(sparse_taucs_call)
#endif
    for (int j = 0; j < nrhs && !rc; j++) {
        taucs_vec_permute(n, TAUCS_DOUBLE, (void *)(B + (size_t)j * n), f->pb, f->perm);
        if (taucs_supernodal_solve_llt(f->L, f->px, f->pb) != 0) rc = -1;
        else taucs_vec_ipermute(n, TAUCS_DOUBLE, f->px, X + (size_t)j * n, f->perm);
    }
    return rc;
}

int Sparse_factor_solve(SparseFactor_T f, const double *b, double *x)
{
    return Sparse_factor_solve_multi(f, b, x, 1);
}

void Sparse_factor_free(SparseFactor_T *pf)
{
    if (!pf || !*pf) return;
    struct SparseFactor *f = *pf;
    if (f->slot) {
        sp_factor_free(f->L, 1);   /* the slot keeps the symbolic structure */
#pragma omp critical(sparse_pattern_cache)
        f->slot->busy = 0;
    } else if (f->L) {
        sp_factor_free(f->L, 0);
        taucs_free_stub(f->perm); taucs_free_stub(f->invperm);
    }
    if (f->A) sp_ccs_free(f->A);
    free(f->pb); free(f->px);
    free(f);
    *pf = NULL;
}

/* A caller-assembled pattern, factored once symbolically and owned by the
 * factor (never lent to or borrowed from the pattern cache: a SLIM global
 * step refactors the same pattern every iteration, and several groups do so
 * at once). */
int Sparse_factor_spd_ccs(int n, const int *colptr, const int *rowind,
                          const double *values, SparseFactor_T *out)
{
    if (!out) return -1;
    *out = NULL;
    if (n <= 1 || !colptr || !rowind || !values || colptr[0] != 0) return -1;
    const int nnz = colptr[n];
    if (nnz < n) return -1;
    for (int c = 0; c < n; c++) {
        if (colptr[c + 1] <= colptr[c] || rowind[colptr[c]] != c) return -1;   /* diagonal first */
        for (int k = colptr[c] + 1; k < colptr[c + 1]; k++)
            if (rowind[k] <= rowind[k - 1] || rowind[k] >= n) return -1;
    }
    struct SparseFactor *f = (struct SparseFactor *)calloc(1, sizeof *f);
    if (!f) return -1;
    f->n = n;
    f->A = sp_ccs_create(n, nnz);
    f->pb = (double *)malloc((size_t)n * sizeof(double));
    f->px = (double *)malloc((size_t)n * sizeof(double));
    if (!f->A || !f->pb || !f->px) { Sparse_factor_free(&f); return -1; }
    memcpy(f->A->colptr, colptr, (size_t)(n + 1) * sizeof(int));
    memcpy(f->A->rowind, rowind, (size_t)nnz * sizeof(int));
    memcpy(f->A->values.d, values, (size_t)nnz * sizeof(double));
    if (sp_symbolic(f->A, &f->perm, &f->invperm, &f->L)) { f->L = NULL; Sparse_factor_free(&f); return -1; }
    if (sp_numeric(f->A, f->perm, f->invperm, f->L)) { Sparse_factor_free(&f); return -1; }
    *out = f;
    return 0;
}

int Sparse_factor_refactor_ccs(SparseFactor_T f, const double *values)
{
    if (!f || !f->A || !f->L || !values) return -1;
    memcpy(f->A->values.d, values, (size_t)f->A->colptr[f->n] * sizeof(double));
    sp_factor_free(f->L, 1);   /* numeric parts only; the symbolic structure stays */
    return sp_numeric(f->A, f->perm, f->invperm, f->L);
}

/* ---- repeated solves of one pattern (2026-10-27) ----------------------------- */

struct SparseRepeat {
    int n, nt;
    int *rows, *cols;              /* the triplet pattern, in the caller's order */
    int *slot;                     /* [nt] each triplet's position in A's values */
    taucs_ccs_matrix *A;           /* the CCS of the latest values */
    int *perm, *invperm; void *L;  /* the ordering and the supernodal factor (numeric part replaced per solve) */
    double *pb, *px;
    int numeric;                   /* L holds a numeric factor */
};

void Sparse_repeat_free(SparseRepeat_T *cache)
{
    if (!cache || !*cache) return;
    struct SparseRepeat *r = *cache;
    if (r->L) sp_factor_free(r->L, 0);
    taucs_free_stub(r->perm); taucs_free_stub(r->invperm);
    if (r->A) sp_ccs_free(r->A);
    free(r->rows); free(r->cols); free(r->slot); free(r->pb); free(r->px); free(r);
    *cache = NULL;
}

/* the pattern's CCS (sp_build_ccs: the same sums), each triplet's slot, the ordering and the symbolic factor */
static int sp_repeat_build(struct SparseRepeat *r, int n, int nt, const int *rows, const int *cols, const double *vals)
{
    r->n = n; r->nt = nt;
    r->A = sp_build_ccs(n, nt, rows, cols, vals);
    r->rows = (int *)malloc((size_t)(nt > 0 ? nt : 1) * sizeof(int)); r->cols = (int *)malloc((size_t)(nt > 0 ? nt : 1) * sizeof(int));
    r->slot = (int *)malloc((size_t)(nt > 0 ? nt : 1) * sizeof(int));
    r->pb = (double *)malloc((size_t)n * sizeof(double)); r->px = (double *)malloc((size_t)n * sizeof(double));
    if (!r->A || !r->rows || !r->cols || !r->slot || !r->pb || !r->px) return -1;
    if (nt) { memcpy(r->rows, rows, (size_t)nt * sizeof(int)); memcpy(r->cols, cols, (size_t)nt * sizeof(int)); }
    for (int t = 0; t < nt; t++) {        /* the row in its column (sorted, one entry a row) */
        int lo = r->A->colptr[cols[t]], hi = r->A->colptr[cols[t] + 1];
        while (lo < hi) { int mid = lo + (hi - lo) / 2; if (r->A->rowind[mid] < rows[t]) lo = mid + 1; else hi = mid; }
        if (lo >= r->A->colptr[cols[t] + 1] || r->A->rowind[lo] != rows[t]) return -1;
        r->slot[t] = lo;
    }
    return sp_symbolic(r->A, &r->perm, &r->invperm, &r->L);
}

int Sparse_repeat_solve(SparseRepeat_T *cache, int n, int nt, const int *rows, const int *cols, const double *vals,
                        const double *b, double *x)
{
    if (!cache) return -1;
    if (n <= SPARSE_DENSE_MAX) return Sparse_solve_sym(n, nt, rows, cols, vals, b, x, SPARSE_SPD);   /* as ever: dense, scalar */
    if (nt < 0 || !rows || !cols || !vals || !b || !x) return -1;
    struct SparseRepeat *r = *cache;
    int same = r && r->n == n && r->nt == nt && (!nt || (!memcmp(r->rows, rows, (size_t)nt * sizeof(int)) &&
                                                      !memcmp(r->cols, cols, (size_t)nt * sizeof(int))));
    if (!same) {
        Sparse_repeat_free(cache);
        for (int t = 0; t < nt; t++) if (cols[t] < 0 || cols[t] >= n || rows[t] < cols[t] || rows[t] >= n) return -1;
        r = (struct SparseRepeat *)calloc(1, sizeof *r); if (!r) return -1;
        *cache = r;
        if (sp_repeat_build(r, n, nt, rows, cols, vals)) { Sparse_repeat_free(cache); return -1; }
    } else {   /* sp_build_ccs's sums: each slot's triplets in their order, from +0 */
        double *v = r->A->values.d; int nnz = r->A->colptr[n];
        for (int k = 0; k < nnz; k++) v[k] = 0.0;
        for (int t = 0; t < nt; t++) v[r->slot[t]] += vals[t];
        if (r->numeric) { sp_factor_free(r->L, 1); r->numeric = 0; }   /* numeric parts only; the symbolic structure stays */
    }
    if (sp_numeric(r->A, r->perm, r->invperm, r->L)) { Sparse_repeat_free(cache); return -1; }
    r->numeric = 1;
    int rc = 0;
#ifndef SPARSE_TAUCS_THREAD_SAFE
#pragma omp critical(sparse_taucs_call)
#endif
    {
        taucs_vec_permute(n, TAUCS_DOUBLE, (void *)b, r->pb, r->perm);
        if (taucs_supernodal_solve_llt(r->L, r->px, r->pb) != 0) rc = -1;
        else taucs_vec_ipermute(n, TAUCS_DOUBLE, r->px, x, r->perm);
    }
    /* the numeric factor is made again by every solve: only the ordering, the symbolic structure and the CCS slots are
     * kept between solves (a growth ring holds one such cache per large block of its band) */
    sp_factor_free(r->L, 1); r->numeric = 0;
    return rc;
}

/* (2026-10-27) Sparse_repeat_solve against Sparse_solve_sym, bit for bit: random SPD patterns with duplicate triplets
 * either side of the dense cut, several value sets each, then a changed pattern (a fresh start) and a value set that is
 * not positive definite (both refuse). */
static int sp_repeat_selftest(void)
{
    unsigned long long st = 0xD1B54A32D192ED03ull;
    int bad = 0, solves = 0;
    for (int n = 40; n <= 460 && !bad; n += 60) {
        int cap = n * 16, nt = 0;
        int *r = (int *)malloc((size_t)cap * sizeof(int)), *c = (int *)malloc((size_t)cap * sizeof(int));
        double *v = (double *)malloc((size_t)cap * sizeof(double)), *b = (double *)malloc((size_t)n * sizeof(double));
        double *x1 = (double *)malloc((size_t)n * sizeof(double)), *x2 = (double *)malloc((size_t)n * sizeof(double));
        if (!r || !c || !v || !b || !x1 || !x2) bad = 1;
        for (int i = 0; i < n && !bad; i++) {
            r[nt] = i; c[nt] = i; nt++;
            r[nt] = i; c[nt] = i; nt++;                                   /* a duplicate */
            for (int k = 1; k <= 6 && i - k >= 0; k++) {
                st ^= st << 13; st ^= st >> 7; st ^= st << 17;
                if (st % 3 == 0) continue;
                int j = (int)(st % (unsigned long long)(k * 5)); if (j > i) j = i - k;   /* near and far below the diagonal */
                if (j == i) continue;
                r[nt] = i; c[nt] = j; nt++;
            }
        }
        SparseRepeat_T cache = NULL;
        for (int round = 0; round < 4 && !bad; round++) {
            for (int t = 0; t < nt; t++) {
                st ^= st << 13; st ^= st >> 7; st ^= st << 17;
                double u = (double)(st >> 11) / 9007199254740992.0;
                v[t] = r[t] == c[t] ? 3.5 + 2.0 * u : (u - 0.5) * 0.4;
            }
            for (int i = 0; i < n; i++) { st ^= st << 13; st ^= st >> 7; st ^= st << 17; b[i] = (double)(st >> 11) / 9007199254740992.0 - .5; }
            int e1 = Sparse_solve_sym(n, nt, r, c, v, b, x1, SPARSE_SPD), e2 = Sparse_repeat_solve(&cache, n, nt, r, c, v, b, x2);
            if (e1 || e2 || memcmp(x1, x2, (size_t)n * sizeof *x1)) bad = 1;
            solves++;
        }
        if (!bad && nt > n + 2) {                                         /* a changed pattern: a fresh start, still exact */
            int keep = nt - 1;
            if (Sparse_solve_sym(n, keep, r, c, v, b, x1, SPARSE_SPD) || Sparse_repeat_solve(&cache, n, keep, r, c, v, b, x2) ||
                memcmp(x1, x2, (size_t)n * sizeof *x1)) bad = 1;
            solves++;
        }
        if (!bad) {                                                       /* not positive definite: both refuse */
            for (int t = 0; t < nt; t++) if (r[t] == c[t] && r[t] == n / 2) v[t] = -50.0;
            int e1 = Sparse_solve_sym(n, nt, r, c, v, b, x1, SPARSE_SPD), e2 = Sparse_repeat_solve(&cache, n, nt, r, c, v, b, x2);
            if (!e1 || !e2) bad = 1;
        }
        Sparse_repeat_free(&cache);
        free(r); free(c); free(v); free(b); free(x1); free(x2);
    }
    fprintf(stderr, "sparse selftest: repeated solves of one pattern == the one-shot solver bit for bit (%d solves): %s\n", solves, bad ? "FAIL" : "ok");
    return bad;
}

/* Cache statistics for the selftest and the logs. */
void Sparse_factor_cache_stats(unsigned long long *hits, unsigned long long *misses)
{
    if (hits) *hits = sp_hits;
    if (misses) *misses = sp_misses;
}

/* ---- selftest ------------------------------------------------------------ */

/* Owned nonlinear-solve contexts use this adapter's exact TAUCS assembly and
 * factor ownership. Keep the public implementation in the canonical library. */
#include "sparse_context.inc"
#include "sparse_context_test.inc"

/* Dense Schur complements reach the largest degree bucket. Dirty scratch
 * beyond that bucket must neither be read nor affect the permutation. */
static int sp_ordering_selftest(void)
{
    enum { N=65 };
    const int limits[]={-1,0,1,N};
    for(int test=0;test<4;test++){
        int n=N,delta=limits[test],maxint=2147483647,nofsub=0,at=0;
        int xadj[N+1],adj[N*(N-1)],invp[N],perm[N],head[N+1],size[N],list[N],mark[N];
        for(int i=0;i<N;i++){
            xadj[i]=at+1;
            for(int j=0;j<N;j++)if(j!=i)adj[at++]=j+1;
        }
        xadj[N]=at+1;head[N]=0x01010101;
        genmmd_(&n,xadj,adj,invp,perm,&delta,head,size,list,mark,&maxint,&nofsub);
        if(head[N]!=0x01010101)return -1;
        for(int i=0;i<N;i++)if(perm[i]<1 || perm[i]>N || invp[perm[i]-1]!=i+1)return -1;
    }
    return 0;
}

#if defined(SPARSE_TAUCS_THREAD_SAFE) && defined(_OPENMP)
/* Unlocked TAUCS must give every thread its serial answer bit for bit. Grid
 * Laplacians of 24 sizes (past the n <= 64 identity ordering, so GENMMD and
 * the supernodal BLAS run), one-shot and factored, are solved serially and
 * then all at once on 16 threads, each system twice at the same time so the
 * pattern cache hands out a busy slot, three times over. */
static int sp_concurrency_selftest(void)
{
    enum { K = 24, REPEAT = 3 };
    int n[K], nt[K], *rows[K], *cols[K], bad = 0, mismatches = 0;
    double *vals[K], *b[K], *one[K], *fac[K];
    for (int k = 0; k < K; k++) {
        int g = 9+3*k; n[k] = g*g; nt[k] = 0;
        rows[k] = malloc(3*(size_t)n[k]*sizeof(int)); cols[k] = malloc(3*(size_t)n[k]*sizeof(int));
        vals[k] = malloc(3*(size_t)n[k]*sizeof(double)); b[k] = malloc((size_t)n[k]*sizeof(double));
        one[k] = malloc((size_t)n[k]*sizeof(double)); fac[k] = malloc((size_t)n[k]*sizeof(double));
        if (!rows[k] || !cols[k] || !vals[k] || !b[k] || !one[k] || !fac[k]) { bad = 1; continue; }
        for (int i = 0; i < n[k]; i++) {
            int x = i%g, y = i/g;
            rows[k][nt[k]] = i; cols[k][nt[k]] = i; vals[k][nt[k]++] = 4.01+.001*(double)((i*7+k)%13);
            if (x > 0) { rows[k][nt[k]] = i; cols[k][nt[k]] = i-1; vals[k][nt[k]++] = -1; }
            if (y > 0) { rows[k][nt[k]] = i; cols[k][nt[k]] = i-g; vals[k][nt[k]++] = -1; }
            b[k][i] = (double)((i*31+k)%17)-8;
        }
        SparseFactor_T f = NULL;
        if (Sparse_solve_sym(n[k],nt[k],rows[k],cols[k],vals[k],b[k],one[k],SPARSE_SPD) ||
            Sparse_factor_spd(n[k],nt[k],rows[k],cols[k],vals[k],&f) || Sparse_factor_solve(f,b[k],fac[k])) bad = 1;
        Sparse_factor_free(&f);
    }
    for (int r = 0; r < REPEAT && !bad; r++) {
        int t = 0;
#pragma omp parallel for schedule(dynamic,1) num_threads(16) reduction(+:mismatches)
        for (t = 0; t < 2*K; t++) {
            int k = t%K; size_t bytes = (size_t)n[k]*sizeof(double);
            double *x = malloc(bytes), *y = malloc(bytes); SparseFactor_T f = NULL;
            if (!x || !y || Sparse_solve_sym(n[k],nt[k],rows[k],cols[k],vals[k],b[k],x,SPARSE_SPD) ||
                Sparse_factor_spd(n[k],nt[k],rows[k],cols[k],vals[k],&f) || Sparse_factor_solve(f,b[k],y) ||
                memcmp(x,one[k],bytes) || memcmp(y,fac[k],bytes)) mismatches++;
            Sparse_factor_free(&f); free(x); free(y);
        }
    }
    for (int k = 0; k < K; k++) { free(rows[k]); free(cols[k]); free(vals[k]); free(b[k]); free(one[k]); free(fac[k]); }
    fprintf(stderr, "sparse selftest: unlocked TAUCS on 16 threads, %d systems x %d rounds: %s (%d mismatches)\n",
            2*K, REPEAT, bad || mismatches ? "FAIL" : "bit-identical", mismatches);
    return bad || mismatches ? -1 : 0;
}
#endif

/* (2026-10-28) this file's blocks come from the calling thread's private heap (Windows, SPARSE_TAUCS_THREAD_SAFE): two
 * threads' blocks name two different heaps and each is freed by the other thread; elsewhere they are the C library's */
static int sp_heap_selftest(void)
{
    void *blk[2] = {NULL, NULL}; int team = 0, bad = 0;
#pragma omp parallel num_threads(2)
    {
        int t = omp_get_thread_num();
        if (t == 0) team = omp_get_num_threads();
        if (t < 2) { blk[t] = malloc(4096); if (blk[t]) memset(blk[t], t+1, 4096); }
    }
    if (!blk[0] || (team > 1 && !blk[1])) bad = 1;
#if defined(_WIN32) && defined(SPARSE_TAUCS_THREAD_SAFE) && !defined(TAUCS_PROCESS_HEAP)
    if (!bad && team > 1) {
        void *h0 = *(void **)((char *)blk[0]-16), *h1 = *(void **)((char *)blk[1]-16);
        if (!h0 || !h1 || h0 == h1) bad = 1;
    }
#endif
#pragma omp parallel num_threads(2)
    {
        int t = omp_get_thread_num(), n = omp_get_num_threads();
        if (n > 1) { if (t < 2 && blk[1-t]) free(blk[1-t]); }   /* each thread frees the other's block */
        else { free(blk[0]); free(blk[1]); }
    }
    fprintf(stderr, "sparse selftest: blocks from per-thread private heaps, freed across threads (%d threads): %s\n", team,
            bad ? "FAIL" : "ok");
    return bad;
}

int Sparse_selftest(void)
{
    if (sp_heap_selftest()) return -1;   /* (2026-10-28) */
    if(sp_ordering_selftest()){
        fprintf(stderr,"sparse selftest: dense degree ordering failed\n");return -1;
    }
    if (sp_repeat_selftest()) return -1;
    /* (2026-10-22) the dense path against TAUCS: random SPD systems (a banded random part + a dominant diagonal, with
     * duplicate triplets) of 2..SPARSE_DENSE_MAX unknowns agree to 1e-10 relative; a negative-definite one is refused
     * by both */
    {
        unsigned long long st = 0x9E3779B97F4A7C15ull;
        int bad = 0; double worst = 0;
        for (int n = 2; n <= SPARSE_DENSE_MAX && !bad; n += 7) {
            int cap = n * 12, nt = 0;
            int *r = (int *)malloc((size_t)cap * sizeof(int)), *c = (int *)malloc((size_t)cap * sizeof(int));
            double *v = (double *)malloc((size_t)cap * sizeof(double)), *b = (double *)malloc((size_t)n * sizeof(double));
            double *xd = (double *)calloc((size_t)n, sizeof(double)), *xt = (double *)calloc((size_t)n, sizeof(double));
            taucs_ccs_matrix *A = NULL;
            if (!r || !c || !v || !b || !xd || !xt) bad = 1;
            for (int i = 0; i < n && !bad; i++) {
                r[nt] = i; c[nt] = i; v[nt] = 3.0 + (double)(i % 5); nt++;
                r[nt] = i; c[nt] = i; v[nt] = 1.5; nt++;                  /* a duplicate */
                for (int k = 1; k <= 4 && i - k >= 0; k++) {
                    st ^= st << 13; st ^= st >> 7; st ^= st << 17;
                    r[nt] = i; c[nt] = i - k; v[nt] = ((double)(st >> 11) / 9007199254740992.0 - 0.5) * 0.9; nt++;
                }
                st ^= st << 13; st ^= st >> 7; st ^= st << 17;
                b[i] = (double)(st >> 11) / 9007199254740992.0 - 0.5;
            }
            if (!bad && sp_dense_spd(n, nt, r, c, v, b, xd)) bad = 1;
            if (!bad && !(A = sp_build_ccs(n, nt, r, c, v))) bad = 1;
            if (!bad && sp_solve_spd_once(A, SPARSE_ORDER_AMD, b, xt)) bad = 1;
            double scale = 0, d = 0;
            for (int i = 0; i < n && !bad; i++) { d = fmax(d, fabs(xd[i] - xt[i])); scale = fmax(scale, fabs(xt[i])); }
            if (!bad && !(d <= 1e-10 * (scale > 1 ? scale : 1))) bad = 1;
            worst = fmax(worst, d);
            if (A) sp_ccs_free(A);
            free(r); free(c); free(v); free(b); free(xd); free(xt);
        }
        {   /* not positive definite: refused */
            int rr[2] = {0, 1}, cc[2] = {0, 1}; double vv[2] = {-1.0, 2.0}, bb[2] = {1, 1}, xx[2] = {0, 0};
            if (sp_dense_spd(2, 2, rr, cc, vv, bb, xx) != -1) bad = 1;
        }
        if (bad) { fprintf(stderr, "sparse selftest: dense small-system path disagrees with TAUCS\n"); return -1; }
        fprintf(stderr, "sparse selftest: dense path for n <= %d agrees with TAUCS (worst %.1e)\n", SPARSE_DENSE_MAX, worst);
    }
#if SPARSE_ORDER_AMD
    /* (2026-10-22) AMD against GENMMD: 2-D grid systems with 3 unknowns a
     * node (the growth's band systems' shape), sides 12..90; the AMD path and
     * the GENMMD path solve to the same answer, AMD's ordering is a
     * permutation, and the ordering times print */
    {
        int bad = 0; double t_amd = 0, t_mmd = 0, lnz_sum = 0; const char *why = ""; int why_side = 0;
        for (int side = 12; side <= 90 && !bad; side += 26) {
            int nn = side * side, n3 = 3 * nn, cap = n3 * 9, nt3 = 0;
            int *r3 = (int *)malloc((size_t)cap * sizeof(int)), *c3 = (int *)malloc((size_t)cap * sizeof(int));
            double *v3 = (double *)malloc((size_t)cap * sizeof(double)), *b3 = (double *)malloc((size_t)n3 * sizeof(double));
            double *xa = (double *)calloc((size_t)n3, sizeof(double)), *xg = (double *)calloc((size_t)n3, sizeof(double));
            if (!r3 || !c3 || !v3 || !b3 || !xa || !xg) bad = 1;
            for (int j = 0; j < side && !bad; j++) for (int i = 0; i < side; i++) {
                int p = j * side + i;
                for (int d = 0; d < 3; d++) {
                    int q = 3 * p + d;
                    r3[nt3] = q; c3[nt3] = q; v3[nt3] = 4.6 + 0.1 * d; nt3++;   /* SPD: grid part > 0.6, the node block >= -0.34 */
                    if (d > 0) { r3[nt3] = q; c3[nt3] = q - 1; v3[nt3] = -0.3; nt3++; }      /* within the node */
                    if (i > 0) { r3[nt3] = q; c3[nt3] = q - 3; v3[nt3] = -1.0; nt3++; }      /* the grid's edges */
                    if (j > 0) { r3[nt3] = q; c3[nt3] = q - 3 * side; v3[nt3] = -1.0; nt3++; }
                    b3[q] = sin(0.37 * q) + 0.5;
                }
            }
            taucs_ccs_matrix *A = bad ? NULL : sp_build_ccs(n3, nt3, r3, c3, v3);
            if (!A) { bad = 1; why = "build"; }
            int *pa = NULL, *ia = NULL, *pg = NULL, *ig = NULL; double lnz = -1;
            if (!bad) {
                double t0 = (double)clock() / CLOCKS_PER_SEC;
                if (sp_order(A, 1, &pa, &ia, &lnz)) { bad = 1; why = "AMD order"; }
                double t1 = (double)clock() / CLOCKS_PER_SEC;
                if (!bad && sp_order(A, 0, &pg, &ig, NULL)) { bad = 1; why = "GENMMD order"; }
                double t2 = (double)clock() / CLOCKS_PER_SEC;
                t_amd += t1 - t0; t_mmd += t2 - t1; lnz_sum += lnz > 0 ? lnz : 0;
            }
            if (!bad) {                               /* a permutation, and invperm its inverse */
                char *seen = (char *)calloc((size_t)n3, 1);
                if (!seen) bad = 1;
                for (int k = 0; k < n3 && !bad; k++) {
                    if (pa[k] < 0 || pa[k] >= n3 || seen[pa[k]] || ia[pa[k]] != k) { bad = 1; why = "not a permutation"; }
                    else seen[pa[k]] = 1;
                }
                free(seen);
            }
            if (!bad && sp_solve_spd_once(A, 1, b3, xa)) { bad = 1; why = "AMD solve"; }
            if (!bad && sp_solve_spd_once(A, 0, b3, xg)) { bad = 1; why = "GENMMD solve"; }
            double worst = 0, scale = 0;
            for (int q = 0; q < n3 && !bad; q++) { worst = fmax(worst, fabs(xa[q] - xg[q])); scale = fmax(scale, fabs(xg[q])); }
            if (!bad && !(worst <= 1e-10 * (scale > 1 ? scale : 1))) { bad = 1; why = "answers differ"; }
            if (bad && !why_side) why_side = side;
            taucs_free_stub(pa); taucs_free_stub(ia); taucs_free_stub(pg); taucs_free_stub(ig);
            if (A) sp_ccs_free(A);
            free(r3); free(c3); free(v3); free(b3); free(xa); free(xg);
        }
        if (!bad) {                                   /* the block graph: found on a 3-block system, refused without one */
            int rr[6] = {0, 1, 2, 3, 4, 5}, cc[6] = {0, 1, 2, 3, 4, 5}; double vv[6] = {2, 2, 2, 2, 2, 2};
            int r2[9], c2[9]; double v2[9]; int k2 = 0;
            for (int q = 0; q < 6; q++) { r2[k2] = q; c2[k2] = q; v2[k2] = 4; k2++; }
            r2[k2] = 3; c2[k2] = 0; v2[k2] = -1; k2++;   /* row block 1 from column 0 only: not a block pattern */
            taucs_ccs_matrix *D = sp_build_ccs(6, 6, rr, cc, vv), *E = sp_build_ccs(6, k2, r2, c2, v2);
            int *bc = NULL, *bri = NULL, gd = D ? sp_block_graph(D, 3, &bc, &bri) : -1;
            free(bc); free(bri); bc = NULL; bri = NULL;
            int ge = E ? sp_block_graph(E, 3, &bc, &bri) : -1;
            free(bc); free(bri);
            if (D) sp_ccs_free(D);
            if (E) sp_ccs_free(E);
            if (gd != 1 || ge != 0) { bad = 1; why = "block graph detection"; }
        }
        if (bad) { fprintf(stderr, "sparse selftest: AMD ordering / solve against GENMMD failed (%s, side %d)\n", why, why_side); return -1; }
        fprintf(stderr, "sparse selftest: AMD ordering agrees with GENMMD on 3 grid systems (orderings %.3f s AMD, %.3f s GENMMD; "
                        "AMD's L %.0f off-diagonal nonzeros)\n", t_amd, t_mmd, lnz_sum);
    }
#endif
    /* (2026-10-22) the counting-sort CCS build against the qsort one: random
     * lower-triangle triplets with many duplicates and integer values (sums are
     * exact in any order), so the two matrices must match exactly */
    {
        enum { TN = 57, TT = 4000 };
        int *tr = (int *)malloc(TT * sizeof(int)), *tc = (int *)malloc(TT * sizeof(int));
        double *tv = (double *)malloc(TT * sizeof(double));
        unsigned long long st = 0x2545F4914F6CDD1Dull;
        int bad = !tr || !tc || !tv;
        for (int t = 0; !bad && t < TT; t++) {
            st ^= st << 13; st ^= st >> 7; st ^= st << 17;
            int a = (int)((st >> 8) % TN), b = (int)((st >> 24) % TN);
            tr[t] = a > b ? a : b; tc[t] = a > b ? b : a; tv[t] = (double)((int)((st >> 40) % 199) - 99);
        }
        taucs_ccs_matrix *A1 = bad ? NULL : sp_build_ccs(TN, TT, tr, tc, tv);
        taucs_ccs_matrix *A2 = bad ? NULL : sp_build_ccs_qsort(TN, TT, tr, tc, tv);
        bad = !A1 || !A2;
        for (int c = 0; !bad && c <= TN; c++) bad |= A1->colptr[c] != A2->colptr[c];
        for (int k = 0; !bad && k < A1->colptr[TN]; k++) bad |= A1->rowind[k] != A2->rowind[k] || A1->values.d[k] != A2->values.d[k];
        int bad_input = 0;                                   /* an upper-triangle triplet is refused */
        if (!bad) { int r0 = 0, c0 = 1; double v0 = 1.0; taucs_ccs_matrix *A3 = sp_build_ccs(2, 1, &r0, &c0, &v0); bad_input = A3 != NULL; if (A3) sp_ccs_free(A3); }
        if (A1) sp_ccs_free(A1);
        if (A2) sp_ccs_free(A2);
        free(tr); free(tc); free(tv);
        if (bad || bad_input) { fprintf(stderr, "sparse selftest: counting-sort CCS build differs from the qsort build\n"); return -1; }
    }
    /* 6-node SPD: path-graph Laplacian + 0.5 I, two right-hand sides.  The
     * factored solves must agree with the one-shot solver to solver noise. */
    enum { N = 6 };
    int rows[2 * N], cols[2 * N];
    double vals[2 * N];
    int nt = 0;
    for (int i = 0; i < N; i++) {
        double d = 0.5 + (i > 0 ? 1.0 : 0.0) + (i + 1 < N ? 1.0 : 0.0);
        rows[nt] = i; cols[nt] = i; vals[nt] = d; nt++;
        if (i > 0) { rows[nt] = i; cols[nt] = i - 1; vals[nt] = -1.0; nt++; }
    }
    double b1[N], b2[N], x_ref[N], x_fac[N];
    for (int i = 0; i < N; i++) {
        b1[i] = 1.0 + (double)i;
        b2[i] = ((i % 2) == 0) ? 3.0 : -2.0;
    }

    SparseFactor_T f = NULL;
    if (Sparse_factor_spd(N, nt, rows, cols, vals, &f) != 0) {
        fprintf(stderr, "sparse selftest: factor failed\n");
        return -1;
    }
    const double *bs[2] = { b1, b2 };
    for (int k = 0; k < 2; k++) {
        memset(x_ref, 0, sizeof x_ref);
        memset(x_fac, 0, sizeof x_fac);
        if (Sparse_solve_sym(N, nt, rows, cols, vals, bs[k], x_ref,
                             SPARSE_SPD) != 0 ||
            Sparse_factor_solve(f, bs[k], x_fac) != 0) {
            fprintf(stderr, "sparse selftest: solve %d failed\n", k);
            Sparse_factor_free(&f);
            return -1;
        }
        for (int i = 0; i < N; i++) {
            if (fabs(x_ref[i] - x_fac[i]) > 1e-9) {
                fprintf(stderr,
                        "sparse selftest: rhs %d mismatch at %d: %g vs %g\n",
                        k, i, x_ref[i], x_fac[i]);
                Sparse_factor_free(&f);
                return -1;
            }
        }
    }
    Sparse_factor_free(&f);

    /* The same pattern again is a cache hit (numeric only): identical answers.
     * An indefinite matrix on that pattern must fail cleanly and leave the
     * slot usable, and a live factor must not share its slot. */
    {
        unsigned long long h0, m0, h1, m1;
        Sparse_factor_cache_stats(&h0, &m0);
        SparseFactor_T g = NULL, g2 = NULL;
        double x_hit[N], vbad[2 * N];
        memcpy(vbad, vals, sizeof vbad); vbad[0] = -5.0;   /* negative pivot */
        SparseFactor_T bad = NULL;
        int good = Sparse_solve_sym(N, nt, rows, cols, vals, b1, x_ref, SPARSE_SPD) == 0;   /* x_ref: the b1 reference */
        good = good && Sparse_factor_spd(N, nt, rows, cols, vals, &g) == 0 && Sparse_factor_solve(g, b1, x_hit) == 0;
        for (int i = 0; good && i < N; i++) if (fabs(x_hit[i] - x_ref[i]) > 1e-9) good = 0;
        good = good && Sparse_factor_spd(N, nt, rows, cols, vals, &g2) == 0 && g2;   /* slot busy: a second, owned factor */
        Sparse_factor_free(&g); Sparse_factor_free(&g2);
        good = good && Sparse_factor_spd(N, nt, rows, cols, vbad, &bad) != 0 && bad == NULL;
        good = good && Sparse_factor_spd(N, nt, rows, cols, vals, &g) == 0 && Sparse_factor_solve(g, b1, x_hit) == 0;
        for (int i = 0; good && i < N; i++) if (fabs(x_hit[i] - x_ref[i]) > 1e-9) good = 0;
        Sparse_factor_free(&g);
        Sparse_factor_cache_stats(&h1, &m1);
        if (!good || h1 - h0 < 2) {
            fprintf(stderr, "sparse selftest: pattern cache failed (good %d, hits +%llu, misses +%llu)\n", good, h1 - h0, m1 - m0);
            return -1;
        }
    }

    /* n==1 scalar path */
    {
        int r0 = 0, c0 = 0;
        double v0 = 4.0, b0 = 10.0, x0 = 0.0;
        SparseFactor_T s = NULL;
        if (Sparse_factor_spd(1, 1, &r0, &c0, &v0, &s) != 0 ||
            Sparse_factor_solve(s, &b0, &x0) != 0 ||
            fabs(x0 - 2.5) > 1e-12) {
            fprintf(stderr, "sparse selftest: scalar path failed\n");
            Sparse_factor_free(&s);
            return -1;
        }
        Sparse_factor_free(&s);
    }
    { unsigned long long h, m; Sparse_factor_cache_stats(&h, &m);
      fprintf(stderr, "sparse selftest: factor reuse OK (2 rhs + scalar); pattern cache hits %llu misses %llu\n", h, m); }
#if defined(SPARSE_TAUCS_THREAD_SAFE) && defined(_OPENMP)
    if (sp_concurrency_selftest()) return -1;
#endif
    return 0;
}
