#include "sparse_solve.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "../common/ves_omp.h"

#include "taucs.h"

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
 * duplicates summed).  Returns NULL on invalid input or allocation failure. */
static taucs_ccs_matrix *sp_build_ccs(int n, int nt,
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

    taucs_ccs_matrix *A = sp_build_ccs(n, nt, rows, cols, vals);
    if (!A) return -1;

    /* --- Solve. --- */
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
};

/* Ordering + symbolic factor for a pattern nobody has cached. */
static int sp_symbolic(taucs_ccs_matrix *A, int **perm, int **invperm, void **L)
{
    const int n = A->n;
    *perm = NULL; *invperm = NULL; *L = NULL;
    if (n <= 64) {
        /* TAUCS's legacy f2c GENMMD path is not reliable on very small
         * graphs (a 13-vertex disk can fault under /O2), and ordering cannot
         * materially reduce fill at this scale: the identity permutation. */
        *perm = (int *)taucs_malloc_stub((size_t)n * sizeof(int)); *invperm = (int *)taucs_malloc_stub((size_t)n * sizeof(int));
        if (!*perm || !*invperm) { taucs_free_stub(*perm); taucs_free_stub(*invperm); *perm = *invperm = NULL; return -1; }
        for (int i = 0; i < n; i++) { (*perm)[i] = i; (*invperm)[i] = i; }
    } else {
#ifndef SPARSE_TAUCS_THREAD_SAFE
#pragma omp critical(sparse_taucs_call)
#endif
        taucs_ccs_order(A, perm, invperm, "genmmd");
        if (!*perm || !*invperm) { taucs_free_stub(*perm); taucs_free_stub(*invperm); *perm = *invperm = NULL; return -1; }
    }
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
    free(f->pb); free(f->px);
    free(f);
    *pf = NULL;
}

/* Cache statistics for the selftest and the logs. */
void Sparse_factor_cache_stats(unsigned long long *hits, unsigned long long *misses)
{
    if (hits) *hits = sp_hits;
    if (misses) *misses = sp_misses;
}

/* ---- selftest ------------------------------------------------------------ */

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

int Sparse_selftest(void)
{
    if(sp_ordering_selftest()){
        fprintf(stderr,"sparse selftest: dense degree ordering failed\n");return -1;
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
