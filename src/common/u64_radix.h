#ifndef VES_U64_RADIX_INCLUDED
#define VES_U64_RADIX_INCLUDED

/*
 * Deterministic multicore radix sorting for the key-first 16-byte records used
 * by region-scale mesh edge indices.  Equal keys retain their input order
 * independently of thread count: workers own contiguous input bands and their
 * bucket offsets are assigned in worker order.
 */
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct {
    uint64_t key;
    uint64_t data;
} VesU64Record16;

static int ves_u64_radix_threads(size_t n)
{
#ifdef _OPENMP
    int nt = omp_get_max_threads();
    size_t useful = (n + 262143u) / 262144u;
    if (useful < 1) useful = 1;
    if ((size_t)nt > useful) nt = (int)useful;
    return nt > 0 ? nt : 1;
#else
    (void)n;
    return 1;
#endif
}

static int ves_u64_record16_sort(VesU64Record16 *a, size_t n)
{
    enum { NB = 65536 };
    VesU64Record16 *b;
    size_t *hist;
    int nt;
    if (n < 2) return 0;
    if (n > SIZE_MAX / sizeof(*b)) return -1;
    nt = ves_u64_radix_threads(n);
    b = (VesU64Record16 *)malloc(n * sizeof(*b));
    hist = (size_t *)calloc((size_t)nt * (size_t)NB, sizeof(*hist));
    if (b == NULL || hist == NULL) {
        free(b); free(hist);
        return -1;
    }
    {
        VesU64Record16 *src = a, *dst = b;
        int pass;
        for (pass = 0; pass < 4; pass++) {
            unsigned shift = (unsigned)pass * 16u;
            size_t sum = 0;
            memset(hist, 0,
                   (size_t)nt * (size_t)NB * sizeof(*hist));
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
            {
                int tid = 0;
#ifdef _OPENMP
                tid = omp_get_thread_num();
#endif
                size_t lo = n * (size_t)tid / (size_t)nt;
                size_t hi = n * (size_t)(tid + 1) / (size_t)nt;
                size_t *h = hist + (size_t)tid * (size_t)NB;
                size_t i;
                for (i = lo; i < hi; i++)
                    h[(unsigned)((src[i].key >> shift) &
                                 UINT64_C(0xffff))]++;
            }
            {
                size_t k;
                for (k = 0; k < (size_t)NB; k++) {
                    int tid;
                    for (tid = 0; tid < nt; tid++) {
                        size_t *p = &hist[(size_t)tid * (size_t)NB + k];
                        size_t count = *p;
                        *p = sum;
                        sum += count;
                    }
                }
            }
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
            {
                int tid = 0;
#ifdef _OPENMP
                tid = omp_get_thread_num();
#endif
                size_t lo = n * (size_t)tid / (size_t)nt;
                size_t hi = n * (size_t)(tid + 1) / (size_t)nt;
                size_t *h = hist + (size_t)tid * (size_t)NB;
                size_t i;
                for (i = lo; i < hi; i++) {
                    unsigned k = (unsigned)
                        ((src[i].key >> shift) & UINT64_C(0xffff));
                    dst[h[k]++] = src[i];
                }
            }
            {
                VesU64Record16 *tmp = src;
                src = dst;
                dst = tmp;
            }
        }
        /* Four passes return to a; retain this guard if the pass count changes. */
        if (src != a) memcpy(a, src, n * sizeof(*a));
    }
    free(b); free(hist);
    return 0;
}

#endif
