/* kdtree_test — exhaustive validation of the KD-tree used by the CVT/RVD projection
 * and the RVD candidate/fast-path queries. Every KDTree_nearest / KDTree_ball_query
 * result is checked against an O(n) brute-force oracle over many random point sets,
 * distributions, and edge cases. Returns 0 on pass. */
#include "kdtree.h"
#include "arena.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>

static uint32_t RNG = 0x1234567u;
static uint32_t xr(void) { uint32_t x=RNG; x^=x<<13; x^=x>>17; x^=x<<5; RNG=x; return x; }
static float    fr(void) { return (float)(xr()>>8) * (1.0f/16777216.0f); }         /* [0,1) */
static float    frr(float lo, float hi) { return lo + (hi-lo)*fr(); }

static float d2f(const float a[3], const float b[3]) {
    float dz=a[0]-b[0], dy=a[1]-b[1], dx=a[2]-b[2];
    return dz*dz+dy*dy+dx*dx;
}

/* brute nearest: returns index of min-distance point, *out_d2 = that squared dist */
static size_t brute_nearest(const float *pts, size_t n, const float q[3], float *out_d2) {
    float best = INFINITY; size_t bi = 0;
    for (size_t i=0;i<n;i++) { float d=d2f(&pts[i*3], q); if (d<best){best=d;bi=i;} }
    *out_d2 = best; return bi;
}

/* brute ball: writes sorted indices with d2 <= r2 into out (cap must hold all), returns count */
static int cmp_i32(const void *a, const void *b){ int32_t x=*(const int32_t*)a,y=*(const int32_t*)b; return x<y?-1:(x>y?1:0); }
static size_t brute_ball(const float *pts, size_t n, const float c[3], float r2, int32_t *out) {
    size_t k=0;
    for (size_t i=0;i<n;i++) if (d2f(&pts[i*3], c) <= r2) out[k++]=(int32_t)i;
    qsort(out, k, sizeof(int32_t), cmp_i32);
    return k;
}

static int g_fail = 0;
#define REQUIRE(cond, ...) do { if (!(cond)) { fprintf(stderr, "FAIL: "); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n"); g_fail=1; } } while (0)

/* Relative tolerance for a squared-distance comparison (float math, both sides). */
static int close_d2(float a, float b) {
    float d = fabsf(a-b);
    return d <= 1e-3f * (fabsf(a)+fabsf(b)+1.0f);
}

/* ---- one random configuration: build tree, hammer nearest + ball vs brute ---- */
static void test_config(Arena_T a, const float *pts, size_t n, float lo, float hi, int nq) {
    KDTree_T kd = KDTree_new(a, pts, n);
    int32_t *kbuf = (int32_t*)ARENA_ALLOC(a, (long)((n+1)*sizeof(int32_t)));
    int32_t *bbuf = (int32_t*)ARENA_ALLOC(a, (long)((n+1)*sizeof(int32_t)));

    for (int t=0;t<nq;t++) {
        /* query: mostly inside the box, sometimes well outside */
        float q[3];
        if (t % 4 == 0) { q[0]=frr(lo-3*(hi-lo), hi+3*(hi-lo)); q[1]=frr(lo-2, hi+2); q[2]=frr(lo-2, hi+2); }
        else            { q[0]=frr(lo, hi); q[1]=frr(lo, hi); q[2]=frr(lo, hi); }
        if (t % 7 == 1) { size_t j = xr()%n; q[0]=pts[j*3]; q[1]=pts[j*3+1]; q[2]=pts[j*3+2]; } /* on a point */

        /* --- nearest --- */
        float kd_d2 = -1.0f;
        size_t ki = KDTree_nearest(kd, q, &kd_d2);
        float bf_d2; size_t bi = brute_nearest(pts, n, q, &bf_d2);
        (void)bi;
        REQUIRE(ki < n, "nearest index %zu out of range n=%zu", ki, n);
        REQUIRE(close_d2(kd_d2, bf_d2), "nearest dist mismatch: kd=%g brute=%g (n=%zu q=%g,%g,%g)",
                kd_d2, bf_d2, n, q[0],q[1],q[2]);
        /* returned index must actually be at the returned distance (original-index convention) */
        REQUIRE(close_d2(d2f(&pts[ki*3], q), kd_d2), "nearest index<->dist inconsistency: idx d2=%g out=%g",
                d2f(&pts[ki*3], q), kd_d2);

        /* --- ball query --- (random radius; a few chosen to bracket the nearest) */
        float rr;
        int mode = t % 5;
        if      (mode==0) rr = 0.0f;                        /* only coincident points */
        else if (mode==1) rr = sqrtf(bf_d2) * 0.999f;       /* just inside -> excludes nearest */
        else if (mode==2) rr = sqrtf(bf_d2) * 1.001f + 1e-4f;/* just outside -> includes nearest */
        else if (mode==3) rr = (hi-lo) * 100.0f;            /* everything */
        else              rr = frr(0.0f, (hi-lo)*1.5f);     /* random */
        float r2 = rr*rr;
        size_t kn = KDTree_ball_query(kd, q, r2, kbuf, n+1);
        size_t bn = brute_ball(pts, n, q, r2, bbuf);
        qsort(kbuf, kn, sizeof(int32_t), cmp_i32);
        REQUIRE(kn == bn, "ball count mismatch: kd=%zu brute=%zu (n=%zu r=%g)", kn, bn, n, rr);
        if (kn == bn) {
            for (size_t i=0;i<kn;i++)
                REQUIRE(kbuf[i]==bbuf[i], "ball set mismatch at %zu: kd=%d brute=%d (n=%zu)", i, kbuf[i], bbuf[i], n);
            /* every returned point must genuinely be within the radius */
            for (size_t i=0;i<kn;i++)
                REQUIRE(d2f(&pts[(size_t)kbuf[i]*3], q) <= r2 + 1e-3f*(r2+1.0f),
                        "ball returned a point outside radius (n=%zu)", n);
        }
    }
}

int main(void) {
    Arena_T a = Arena_new();

    /* Suite 1: uniform random clouds across a range of sizes */
    size_t sizes[] = {1,2,3,4,5,8,17,64,257,1000,5000,20000};
    for (size_t si=0; si<sizeof(sizes)/sizeof(sizes[0]); si++) {
        size_t n = sizes[si];
        for (int trial=0; trial<8; trial++) {
            float *pts = (float*)ARENA_ALLOC(a, (long)(n*3*sizeof(float)));
            for (size_t i=0;i<n*3;i++) pts[i] = frr(-50.0f, 50.0f);
            int nq = n < 1000 ? 400 : 80;
            test_config(a, pts, n, -50.0f, 50.0f, nq);
        }
    }

    /* Suite 2: degenerate distributions */
    {
        /* all identical points */
        size_t n=200; float *pts=(float*)ARENA_ALLOC(a,(long)(n*3*sizeof(float)));
        for (size_t i=0;i<n;i++){pts[i*3]=3.0f;pts[i*3+1]=-1.0f;pts[i*3+2]=2.0f;}
        test_config(a, pts, n, -5.0f, 5.0f, 300);
    }
    {
        /* collinear points along X */
        size_t n=500; float *pts=(float*)ARENA_ALLOC(a,(long)(n*3*sizeof(float)));
        for (size_t i=0;i<n;i++){pts[i*3]=0.0f;pts[i*3+1]=0.0f;pts[i*3+2]=(float)i*0.1f;}
        test_config(a, pts, n, -1.0f, 51.0f, 400);
    }
    {
        /* coplanar grid (z=0) */
        int M=40; size_t n=(size_t)M*M; float *pts=(float*)ARENA_ALLOC(a,(long)(n*3*sizeof(float)));
        for (int i=0;i<M;i++)for(int j=0;j<M;j++){size_t k=(size_t)i*M+j; pts[k*3]=0.0f; pts[k*3+1]=(float)i; pts[k*3+2]=(float)j;}
        test_config(a, pts, n, -2.0f, 42.0f, 500);
    }
    {
        /* tightly-clustered points (dense, many near-ties) */
        size_t n=800; float *pts=(float*)ARENA_ALLOC(a,(long)(n*3*sizeof(float)));
        for (size_t i=0;i<n*3;i++) pts[i]=frr(0.0f, 0.01f);
        test_config(a, pts, n, 0.0f, 0.01f, 500);
    }
    {
        /* huge coordinate magnitudes (float precision stress) */
        size_t n=400; float *pts=(float*)ARENA_ALLOC(a,(long)(n*3*sizeof(float)));
        for (size_t i=0;i<n*3;i++) pts[i]=frr(-1e6f, 1e6f);
        test_config(a, pts, n, -1e6f, 1e6f, 400);
    }

    /* Suite 3: determinism — identical build+query twice must agree exactly */
    {
        size_t n=2000; float *pts=(float*)ARENA_ALLOC(a,(long)(n*3*sizeof(float)));
        for (size_t i=0;i<n*3;i++) pts[i]=frr(-10.0f,10.0f);
        KDTree_T k1=KDTree_new(a,pts,n), k2=KDTree_new(a,pts,n);
        for (int t=0;t<500;t++){
            float q[3]={frr(-12,12),frr(-12,12),frr(-12,12)};
            float d1,d2; size_t i1=KDTree_nearest(k1,q,&d1), i2=KDTree_nearest(k2,q,&d2);
            REQUIRE(i1==i2 && d1==d2, "nondeterministic nearest");
            int32_t b1[2049], b2[2049];
            float r2 = frr(0,20); r2*=r2;
            size_t c1=KDTree_ball_query(k1,q,r2,b1,2049), c2=KDTree_ball_query(k2,q,r2,b2,2049);
            qsort(b1,c1,sizeof(int32_t),cmp_i32); qsort(b2,c2,sizeof(int32_t),cmp_i32);
            REQUIRE(c1==c2, "nondeterministic ball count");
            if (c1==c2) for (size_t i=0;i<c1;i++) REQUIRE(b1[i]==b2[i], "nondeterministic ball set");
        }
    }

    Arena_dispose(&a);
    fprintf(stderr, "kdtree_test: %s\n", g_fail ? "FAIL" : "OK");
    return g_fail ? 1 : 0;
}
