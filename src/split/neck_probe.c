#include "neck_probe.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define NP_RAD2DEG 57.29577951308232f

/* Undirected half-edge carrying its owning face, for edge->face pairing. */
typedef struct { int32_t lo, hi, face; } HEdge;

static int he_cmp(const void *pa, const void *pb)
{
    const HEdge *a = (const HEdge *)pa, *b = (const HEdge *)pb;
    if (a->lo != b->lo) return (a->lo < b->lo) ? -1 : 1;
    if (a->hi != b->hi) return (a->hi < b->hi) ? -1 : 1;
    return 0;
}

/* Unit face normal in (z,y,x); zero vector for a degenerate triangle. */
static void np_face_normal(const float *verts, const int32_t *f, float n[3])
{
    const float *a = &verts[(size_t)f[0] * 3];
    const float *b = &verts[(size_t)f[1] * 3];
    const float *c = &verts[(size_t)f[2] * 3];
    float e1[3], e2[3];
    int k;
    for (k = 0; k < 3; k++) { e1[k] = b[k] - a[k]; e2[k] = c[k] - a[k]; }
    n[0] = e1[1] * e2[2] - e1[2] * e2[1];
    n[1] = e1[2] * e2[0] - e1[0] * e2[2];
    n[2] = e1[0] * e2[1] - e1[1] * e2[0];
    float L = sqrtf(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    if (L > 1e-12f) { n[0] /= L; n[1] /= L; n[2] /= L; }
    else { n[0] = 0.0f; n[1] = 0.0f; n[2] = 0.0f; }
}

void NeckProbe_default_opts(NeckProbeOpts *opts)
{
    if (!opts) return;
    opts->deg_thresh     = 60.0f;
    opts->min_neighbours = 2;
}

int NeckProbe_scan(Arena_T arena,
                   const float *verts, size_t nv,
                   const int32_t *faces, size_t nf,
                   const NeckProbeOpts *opts,
                   NeckProbeResult *out)
{
    assert(arena && out);
    memset(out, 0, sizeof(*out));
    if (!verts || !faces || nv == 0 || nf == 0) return -1;

    NeckProbeOpts o;
    NeckProbe_default_opts(&o);
    if (opts) o = *opts;
    float cos_thresh = cosf(o.deg_thresh / NP_RAD2DEG);

    /* Persistent outputs (allocated below the scratch mark so they survive). */
    uint8_t *face_flag = (uint8_t *)ARENA_CALLOC(arena, (long)nf, 1L);
    float   *face_dev  = (float *)  ARENA_ALLOC(arena, (long)(nf * sizeof(float)));

    Arena_Mark mark = Arena_save(arena);

    /* Face normals. */
    float *fn = (float *)ARENA_ALLOC(arena, (long)(nf * 3 * sizeof(float)));
    for (size_t fi = 0; fi < nf; fi++)
        np_face_normal(verts, &faces[fi * 3], &fn[fi * 3]);

    /* Half-edges, sorted -> runs of equal (lo,hi) are the faces sharing an edge. */
    size_t hcap = nf * 3;
    HEdge *he = (HEdge *)ARENA_ALLOC(arena, (long)(hcap * sizeof(HEdge)));
    size_t hk = 0;
    for (size_t fi = 0; fi < nf; fi++) {
        const int32_t *f = &faces[fi * 3];
        for (int e = 0; e < 3; e++) {
            int32_t a = f[e], b = f[(e + 1) % 3];
            if (a < 0 || b < 0 || (size_t)a >= nv || (size_t)b >= nv || a == b) continue;
            he[hk].lo   = (a < b) ? a : b;
            he[hk].hi   = (a < b) ? b : a;
            he[hk].face = (int32_t)fi;
            hk++;
        }
    }
    qsort(he, hk, sizeof(HEdge), he_cmp);

    /* For each face, accumulate the (summed) unit normals of its edge neighbours
     * and a neighbour count. The mean-of-others per shared edge is added in O(run). */
    float   *nbr_sum = (float *)  ARENA_CALLOC(arena, (long)(nf * 3), (long)sizeof(float));
    int32_t *nbr_cnt = (int32_t *)ARENA_CALLOC(arena, (long)nf, (long)sizeof(int32_t));

    size_t i = 0;
    while (i < hk) {
        size_t j = i + 1;
        while (j < hk && he[j].lo == he[i].lo && he[j].hi == he[i].hi) j++;
        size_t run = j - i;
        if (run >= 2) {
            float s[3] = { 0.0f, 0.0f, 0.0f };
            for (size_t t = i; t < j; t++) {
                size_t fc = (size_t)he[t].face;
                s[0] += fn[fc * 3 + 0];
                s[1] += fn[fc * 3 + 1];
                s[2] += fn[fc * 3 + 2];
            }
            for (size_t t = i; t < j; t++) {
                size_t fc = (size_t)he[t].face;
                nbr_sum[fc * 3 + 0] += s[0] - fn[fc * 3 + 0];
                nbr_sum[fc * 3 + 1] += s[1] - fn[fc * 3 + 1];
                nbr_sum[fc * 3 + 2] += s[2] - fn[fc * 3 + 2];
                nbr_cnt[fc] += (int32_t)(run - 1);
            }
        }
        i = j;
    }

    /* Per-face deviation from the neighbourhood mean normal (acute, [0,90]). */
    size_t n_side = 0, n_tips = 0;
    float  maxdev = 0.0f;
    for (size_t fi = 0; fi < nf; fi++) {
        face_dev[fi] = 0.0f;
        int32_t cnt = nbr_cnt[fi];
        if (cnt == 0) continue;                       /* isolated -> no judgement */
        if (cnt < o.min_neighbours) n_tips++;         /* low-support (dangling) tip */

        float m[3] = { nbr_sum[fi * 3 + 0], nbr_sum[fi * 3 + 1], nbr_sum[fi * 3 + 2] };
        float L = sqrtf(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
        if (L < 1e-12f) continue;                     /* neighbours cancel -> skip */
        m[0] /= L; m[1] /= L; m[2] /= L;

        float d  = fn[fi * 3 + 0] * m[0] + fn[fi * 3 + 1] * m[1] + fn[fi * 3 + 2] * m[2];
        float ad = (d < 0.0f) ? -d : d;               /* winding-sign independent */
        if (ad > 1.0f) ad = 1.0f;
        float dev = acosf(ad) * NP_RAD2DEG;
        face_dev[fi] = dev;
        if (dev > maxdev) maxdev = dev;
        if (ad <= cos_thresh) { face_flag[fi] = 1; n_side++; }
    }

    Arena_restore(arena, mark);                       /* frees scratch only */

    out->nf           = nf;
    out->n_sidefacing = n_side;
    out->n_tips       = n_tips;
    out->max_dev_deg  = maxdev;
    out->face_flag    = face_flag;
    out->face_dev     = face_dev;
    return 0;
}

/* ============================================================================
 * Self-test.
 * ==========================================================================*/

static int np_chk(const char *name, int cond)
{
    fprintf(stderr, "  [neckprobe %s] %s\n", name, cond ? "ok" : "FAIL");
    return cond ? 0 : 1;
}

int NeckProbe_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    NeckProbeResult r;

    /* Case 1: flat coplanar quad (two triangles, all normals +z) -> 0 side-facing. */
    {
        float V[] = { 0,0,0,  0,0,1,  0,1,0,  0,1,1 };   /* z=0 plane */
        int32_t F[] = { 0,1,2,  1,3,2 };
        int rc = NeckProbe_scan(arena, V, 4, F, 2, NULL, &r);
        fprintf(stderr, "  flat-quad: rc=%d side=%zu tips=%zu maxdev=%.1f\n",
                rc, r.n_sidefacing, r.n_tips, (double)r.max_dev_deg);
        fails += np_chk("flat-quad", rc == 0 && r.n_sidefacing == 0);
    }

    /* Case 2: flat quad + a third triangle folded 90deg up from the shared edge
     * 1-3 (a perpendicular flap). The flap shares one edge with the sheet and
     * must be flagged side-facing; the two flat triangles must not. */
    {
        /* verts 0..3 in z=0 plane; vert 4 lifted in +z over edge (1,3). */
        float V[] = { 0,0,0,  0,0,1,  0,1,0,  0,1,1,  1,0,1 };
        int32_t F[] = { 0,1,2,  1,3,2,  1,4,3 };   /* F2 = 1,4,3 is the flap */
        int rc = NeckProbe_scan(arena, V, 5, F, 3, NULL, &r);
        int flap_flagged = (rc == 0 && r.face_flag && r.face_flag[2] == 1);
        int flats_clean  = (rc == 0 && r.face_flag &&
                            r.face_flag[0] == 0 && r.face_flag[1] == 0);
        fprintf(stderr, "  flap: rc=%d side=%zu dev2=%.1f flag=[%d %d %d]\n",
                rc, r.n_sidefacing, (double)(r.face_dev ? r.face_dev[2] : -1.0f),
                r.face_flag ? r.face_flag[0] : -1,
                r.face_flag ? r.face_flag[1] : -1,
                r.face_flag ? r.face_flag[2] : -1);
        fails += np_chk("flap-flagged", flap_flagged);
        fails += np_chk("flats-clean", flats_clean);
    }

    /* Case 3: a BROAD 90deg crease -- two 2-triangle quads meeting along a shared
     * full edge. Each fold triangle has an in-plane neighbour, so the neighbourhood
     * mean is only ~45deg off and nothing should be flagged (a broad crease is not
     * a thin bridge). This is the false-positive guard. */
    {
        /* Quad A in z=0 plane: verts 0,1,2,3.  Quad B in y=0 plane sharing edge
         * (0,1): verts 0,1,4,5. */
        float V[] = {
            0,0,0,   0,0,1,   0,1,0,   0,1,1,     /* A: z=0 */
            1,0,0,   1,0,1                          /* B: lifted in +z, y=0 */
        };
        int32_t F[] = { 0,1,2, 1,3,2,    /* quad A */
                        0,4,1, 1,4,5 };  /* quad B shares edge (0,1) */
        int rc = NeckProbe_scan(arena, V, 6, F, 4, NULL, &r);
        fprintf(stderr, "  broad-crease: rc=%d side=%zu maxdev=%.1f\n",
                rc, r.n_sidefacing, (double)r.max_dev_deg);
        fails += np_chk("broad-crease-clean", rc == 0 && r.n_sidefacing == 0);
    }

    Arena_dispose(&arena);
    fprintf(stderr, "=== NeckProbe_selftest %s (%d failure%s) ===\n",
            fails ? "FAILED" : "PASSED", fails, fails == 1 ? "" : "s");
    return fails;
}
