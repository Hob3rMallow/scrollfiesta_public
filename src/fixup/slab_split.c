/* slab_split.c — geometric mid-plane carver for pred-fused wrap slabs.
 * See slab_split.h for the algorithm contract. */
#include "slab_split.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

enum { SS_MAX_THREADS = 128 };

static inline size_t ss_idx(int y, int x, int W)
{
    return (size_t)y * (size_t)W + (size_t)x;
}

void SlabSplit_params_default(SlabSplitParams *p)
{
    assert(p);
    memset(p, 0, sizeof(*p));
    p->thick_r = 2.25f;
    p->ext_r = 1.60f;
    p->grow_max = 8;
    p->rounds = 2;
    p->z_window = 3;
    p->z_support = 2;
    p->plane_cap_frac = 0.10f;
}

/* 3-4 chamfer distance to background, /3 => half-width in px (bridge_scan's
 * metric, reimplemented file-private per house style). */
static void ss_chamfer(const uint8_t *mask, int H, int W, float *out)
{
    const float BIG = 1e6f;
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            out[ss_idx(y, x, W)] = mask[ss_idx(y, x, W)] ? BIG : 0.0f;
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            float v = out[ss_idx(y, x, W)];
            if (v == 0.0f) continue;
            if (y > 0) {
                float c = out[ss_idx(y - 1, x, W)] + 3.0f;
                if (c < v) v = c;
                if (x > 0) {
                    c = out[ss_idx(y - 1, x - 1, W)] + 4.0f;
                    if (c < v) v = c;
                }
                if (x < W - 1) {
                    c = out[ss_idx(y - 1, x + 1, W)] + 4.0f;
                    if (c < v) v = c;
                }
            }
            if (x > 0) {
                float c = out[ss_idx(y, x - 1, W)] + 3.0f;
                if (c < v) v = c;
            }
            out[ss_idx(y, x, W)] = v;
        }
    }
    for (int y = H - 1; y >= 0; y--) {
        for (int x = W - 1; x >= 0; x--) {
            float v = out[ss_idx(y, x, W)];
            if (v == 0.0f) continue;
            if (y < H - 1) {
                float c = out[ss_idx(y + 1, x, W)] + 3.0f;
                if (c < v) v = c;
                if (x > 0) {
                    c = out[ss_idx(y + 1, x - 1, W)] + 4.0f;
                    if (c < v) v = c;
                }
                if (x < W - 1) {
                    c = out[ss_idx(y + 1, x + 1, W)] + 4.0f;
                    if (c < v) v = c;
                }
            }
            if (x < W - 1) {
                float c = out[ss_idx(y, x + 1, W)] + 3.0f;
                if (c < v) v = c;
            }
            out[ss_idx(y, x, W)] = v;
        }
    }
    size_t n = (size_t)H * (size_t)W;
    for (size_t i = 0; i < n; i++) out[i] /= 3.0f;
}

/* Candidate classification for one plane: 0 none, 1 extension, 2 primary.
 * A candidate is a foreground pixel at or above the tier's half-width floor
 * that is an axis local max of the half-width field (medial band). */
static void ss_candidates(const uint8_t *mask, const float *hw, int H, int W,
                          float thick_r, float ext_r, uint8_t *cand)
{
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            size_t i = ss_idx(y, x, W);
            cand[i] = 0;
            if (!mask[i]) continue;
            float v = hw[i];
            if (v < ext_r) continue;
            float up = (y > 0) ? hw[ss_idx(y - 1, x, W)] : 0.0f;
            float dn = (y < H - 1) ? hw[ss_idx(y + 1, x, W)] : 0.0f;
            float lf = (x > 0) ? hw[ss_idx(y, x - 1, W)] : 0.0f;
            float rt = (x < W - 1) ? hw[ss_idx(y, x + 1, W)] : 0.0f;
            /* strict ridge: a max across the slab normal, never a plateau
             * along the tangent (which would flood the whole interior) */
            int ridge_y = (v >= up && v >= dn) && (v > up || v > dn);
            int ridge_x = (v >= lf && v >= rt) && (v > lf || v > rt);
            if (!ridge_y && !ridge_x) continue;
            cand[i] = (uint8_t)(v >= thick_r ? 2 : 1);
        }
    }
}

/* Chebyshev-2 dilation of (cand == 2) into dil (0/1), separable OR. */
static void ss_dilate_primary(const uint8_t *cand, int H, int W,
                              uint8_t *tmp, uint8_t *dil)
{
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            uint8_t v = 0;
            for (int dy = -2; dy <= 2; dy++) {
                int yy = y + dy;
                if (yy < 0 || yy >= H) continue;
                if (cand[ss_idx(yy, x, W)] == 2) { v = 1; break; }
            }
            tmp[ss_idx(y, x, W)] = v;
        }
    }
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            uint8_t v = 0;
            for (int dx = -2; dx <= 2; dx++) {
                int xx = x + dx;
                if (xx < 0 || xx >= W) continue;
                if (tmp[ss_idx(y, xx, W)]) { v = 1; break; }
            }
            dil[ss_idx(y, x, W)] = v;
        }
    }
}

/* One full pass over the volume. Returns carved voxel count. */
static int64_t ss_pass(uint8_t *vol, uint8_t *carve_mask, int D, int H, int W,
                       const SlabSplitParams *p, SlabSplitStats *st,
                       uint8_t *cand_vol, uint8_t *dil_vol,
                       Arena_T *scratch, int n_threads, int round)
{
    size_t plane_sz = (size_t)H * (size_t)W;

    /* stage 1: per-plane candidates + primary dilation */
    {
        int z = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4)
#endif
        for (z = 0; z < D; z++) {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num() % SS_MAX_THREADS;
#endif
            Arena_Mark mark = Arena_save(scratch[tid]);
            float *hw = ARENA_ALLOC(scratch[tid], plane_sz * sizeof(float));
            uint8_t *tmp = ARENA_ALLOC(scratch[tid], plane_sz);
            const uint8_t *pl = vol + (size_t)z * plane_sz;
            uint8_t *cd = cand_vol + (size_t)z * plane_sz;
            ss_chamfer(pl, H, W, hw);
            ss_candidates(pl, hw, H, W, p->thick_r, p->ext_r, cd);
            ss_dilate_primary(cd, H, W, tmp, dil_vol + (size_t)z * plane_sz);
            Arena_restore(scratch[tid], mark);
        }
    }
    (void)n_threads;

    /* stage 2: z-support + hysteresis growth + carve, per plane */
    int64_t carved_total = 0, primary_total = 0, supported_total = 0;
    int64_t planes_carved = 0, planes_capped = 0;
    {
        int z = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4) \
    reduction(+ : carved_total, primary_total, supported_total, \
              planes_carved, planes_capped)
#endif
        for (z = 0; z < D; z++) {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num() % SS_MAX_THREADS;
#endif
            Arena_Mark mark = Arena_save(scratch[tid]);
            uint8_t *cd = cand_vol + (size_t)z * plane_sz;
            uint8_t *seed = ARENA_ALLOC(scratch[tid], plane_sz);
            int32_t *queue = ARENA_ALLOC(scratch[tid],
                                         plane_sz * sizeof(int32_t));
            uint8_t *pl = vol + (size_t)z * plane_sz;

            /* seeds: primaries with z-support */
            size_t n_seed = 0, n_primary = 0;
            memset(seed, 0, plane_sz);
            for (size_t i = 0; i < plane_sz; i++) {
                if (cd[i] != 2) continue;
                n_primary++;
                int sup = 0;
                for (int dz = -p->z_window; dz <= p->z_window; dz++) {
                    int zz = z + dz;
                    if (dz == 0 || zz < 0 || zz >= D) continue;
                    if (dil_vol[(size_t)zz * plane_sz + i]) sup++;
                }
                if (sup >= p->z_support) {
                    seed[i] = 1;
                    queue[n_seed++] = (int32_t)i;
                }
            }
            primary_total += (int64_t)n_primary;
            supported_total += (int64_t)n_seed;
            if (n_seed == 0) {
                Arena_restore(scratch[tid], mark);
                continue;
            }

            /* hysteresis: grow seeds through any candidate tier (8-conn),
             * depth-capped so growth can bridge a fork tip but can never
             * run away down a normal sheet's medial axis */
            size_t head = 0, tail = n_seed, level_end = n_seed;
            int depth = 0;
            while (head < tail && depth < p->grow_max) {
                int32_t i = queue[head++];
                int y = (int)((size_t)i / (size_t)W);
                int x = (int)((size_t)i % (size_t)W);
                for (int dy = -1; dy <= 1; dy++) {
                    for (int dx = -1; dx <= 1; dx++) {
                        if (dy == 0 && dx == 0) continue;
                        int yy = y + dy, xx = x + dx;
                        if (yy < 0 || yy >= H || xx < 0 || xx >= W) continue;
                        size_t j = ss_idx(yy, xx, W);
                        if (seed[j] || cd[j] == 0) continue;
                        seed[j] = 1;
                        queue[tail++] = (int32_t)j;
                    }
                }
                if (head == level_end) {
                    depth++;
                    level_end = tail;
                }
            }

            /* cap check against this plane's foreground */
            size_t fg = 0;
            for (size_t i = 0; i < plane_sz; i++) fg += pl[i];
            if (fg > 0 &&
                (double)tail > (double)p->plane_cap_frac * (double)fg) {
                planes_capped++;
                Arena_restore(scratch[tid], mark);
                continue;
            }

            /* carve */
            size_t carved = 0;
            for (size_t k = 0; k < tail; k++) {
                size_t i = (size_t)queue[k];
                if (!pl[i]) continue;
                pl[i] = 0;
                if (carve_mask) carve_mask[(size_t)z * plane_sz + i] = 1;
                carved++;
            }
            carved_total += (int64_t)carved;
            if (carved > 0 && round == 0) planes_carved++;
            Arena_restore(scratch[tid], mark);
        }
    }

    st->primary_px += primary_total;
    st->supported_px += supported_total;
    st->carved_px += carved_total;
    st->planes_carved += planes_carved;
    st->planes_capped += planes_capped;
    return carved_total;
}

int SlabSplit_run(Arena_T arena, uint8_t *vol, uint8_t *carve_mask,
                  int D, int H, int W,
                  const SlabSplitParams *p, SlabSplitStats *st)
{
    assert(vol && p && st);
    memset(st, 0, sizeof(*st));
    if (D <= 0 || H < 3 || W < 3) return 0;
    if (!(p->thick_r > 0.0f) || !(p->ext_r > 0.0f) ||
        p->ext_r > p->thick_r || p->grow_max < 0 || p->rounds < 1 ||
        p->z_window < 0 || p->z_support < 0 ||
        !(p->plane_cap_frac > 0.0f))
        return -1;

    size_t plane_sz = (size_t)H * (size_t)W;
    size_t vol_sz = (size_t)D * plane_sz;
    uint8_t *cand_vol = ARENA_ALLOC(arena, vol_sz);
    uint8_t *dil_vol = ARENA_ALLOC(arena, vol_sz);

    int n_threads = 1;
#ifdef _OPENMP
    n_threads = omp_get_max_threads();
    if (n_threads > SS_MAX_THREADS) n_threads = SS_MAX_THREADS;
#endif
    Arena_T scratch[SS_MAX_THREADS];
    for (int t = 0; t < n_threads; t++) scratch[t] = Arena_new();

    for (int r = 0; r < p->rounds; r++) {
        int64_t carved = ss_pass(vol, carve_mask, D, H, W, p, st,
                                 cand_vol, dil_vol, scratch, n_threads, r);
        if (carved > 0) st->rounds_run = r + 1;
        if (carved == 0) break;
    }

    for (int t = 0; t < n_threads; t++) Arena_dispose(&scratch[t]);
    return 0;
}

/* ---------------------------------------------------------------- */
/* Selftest                                                         */
/* ---------------------------------------------------------------- */

/* 2D 4-conn component count over a column range of one plane. */
static int ss_cc_count(Arena_T a, const uint8_t *pl, int H, int W,
                       int x0, int x1)
{
    Arena_Mark mark = Arena_save(a);
    int32_t *lbl = ARENA_ALLOC(a, (size_t)H * (size_t)W * sizeof(int32_t));
    int32_t *stack = ARENA_ALLOC(a, (size_t)H * (size_t)W * sizeof(int32_t));
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) lbl[ss_idx(y, x, W)] = 0;
    int n_cc = 0;
    for (int y = 0; y < H; y++) {
        for (int x = x0; x < x1; x++) {
            size_t i = ss_idx(y, x, W);
            if (!pl[i] || lbl[i]) continue;
            n_cc++;
            size_t top = 0;
            stack[top++] = (int32_t)i;
            lbl[i] = n_cc;
            while (top > 0) {
                int32_t c = stack[--top];
                int cy = (int)((size_t)c / (size_t)W);
                int cx = (int)((size_t)c % (size_t)W);
                const int ndy[4] = { -1, 1, 0, 0 };
                const int ndx[4] = { 0, 0, -1, 1 };
                for (int k = 0; k < 4; k++) {
                    int yy = cy + ndy[k], xx = cx + ndx[k];
                    if (yy < 0 || yy >= H || xx < x0 || xx >= x1) continue;
                    size_t j = ss_idx(yy, xx, W);
                    if (!pl[j] || lbl[j]) continue;
                    lbl[j] = n_cc;
                    stack[top++] = (int32_t)j;
                }
            }
        }
    }
    Arena_restore(a, mark);
    return n_cc;
}

int SlabSplit_selftest(void)
{
    enum { D = 9, H = 64, W = 96, SPAN0 = 30, SPAN1 = 70 };
    Arena_T a = Arena_new();
    size_t plane_sz = (size_t)H * (size_t)W;
    size_t vol_sz = (size_t)D * plane_sz;
    uint8_t *vol = ARENA_CALLOC(a, vol_sz, 1);
    uint8_t *mask = ARENA_CALLOC(a, vol_sz, 1);
    int fails = 0;

    /* Sheet A rows 20-22 everywhere; sheet B rows 26-28 outside the span but
     * pressed to rows 23-25 inside it (6-thick fused slab).  Control sheet
     * rows 40-42 everywhere (never thick). */
    for (int z = 0; z < D; z++) {
        uint8_t *pl = vol + (size_t)z * plane_sz;
        for (int x = 0; x < W; x++) {
            for (int y = 20; y <= 22; y++) pl[ss_idx(y, x, W)] = 1;
            int in_span = (x >= SPAN0 && x < SPAN1);
            int b0 = in_span ? 23 : 26, b1 = in_span ? 25 : 28;
            for (int y = b0; y <= b1; y++) pl[ss_idx(y, x, W)] = 1;
            for (int y = 40; y <= 42; y++) pl[ss_idx(y, x, W)] = 1;
        }
    }

    SlabSplitParams p;
    SlabSplit_params_default(&p);
    p.rounds = 1;
    p.plane_cap_frac = 0.5f;   /* the synthetic plane is mostly slab */
    /* the synthetic slab is exactly 6 thick => interior half-width 3.0 */
    SlabSplitStats st;
    if (SlabSplit_run(a, vol, mask, D, H, W, &p, &st) != 0) {
        fprintf(stderr, "[selftest] slab_split: run failed\n");
        fails++;
    }
    if (st.carved_px <= 0) {
        fprintf(stderr, "[selftest] slab_split: carved nothing\n");
        fails++;
    }
    /* the fused span must split into two sheets on every plane (checked away
     * from the span ends where the carve hands over to the genuine fork) */
    for (int z = 0; z < D && !fails; z++) {
        int cc = ss_cc_count(a, vol + (size_t)z * plane_sz, H, W,
                             SPAN0 + 4, SPAN1 - 4);
        /* rows 20-25 fused slab must now be 2 comps; control sheet is a 3rd */
        if (cc != 3) {
            fprintf(stderr,
                    "[selftest] slab_split: plane %d span cc=%d want 3\n",
                    z, cc);
            fails++;
        }
    }
    /* control sheet untouched */
    for (int z = 0; z < D; z++) {
        const uint8_t *pl = vol + (size_t)z * plane_sz;
        for (int x = 0; x < W; x++)
            for (int y = 40; y <= 42; y++)
                if (!pl[ss_idx(y, x, W)]) {
                    fprintf(stderr,
                            "[selftest] slab_split: control sheet carved at "
                            "z=%d y=%d x=%d\n", z, y, x);
                    fails++;
                    z = D; x = W; break;
                }
    }
    /* deep in the span the carve is exactly the medial rows (22-23); near
     * the span-edge forks the interface bends, so rows 21-24 are allowed
     * there plus the bounded hysteresis leak; skins (rows 20, 25) and
     * everything else must be untouched */
    {
        int margin = p.grow_max + 2;
        int bad = 0;
        for (int z = 0; z < D && !bad; z++) {
            const uint8_t *cm = mask + (size_t)z * plane_sz;
            for (int y = 0; y < H && !bad; y++) {
                for (int x = 0; x < W; x++) {
                    if (!cm[ss_idx(y, x, W)]) continue;
                    int deep = (x >= SPAN0 + margin && x < SPAN1 - margin);
                    int row_ok = deep ? (y >= 22 && y <= 23)
                                      : (y >= 21 && y <= 24);
                    int leak_ok = (x >= SPAN0 - margin && x < SPAN1 + margin);
                    if (!row_ok || !leak_ok) {
                        fprintf(stderr,
                                "[selftest] slab_split: carve out of bounds "
                                "at z=%d y=%d x=%d\n", z, y, x);
                        bad = 1;
                        break;
                    }
                }
            }
        }
        fails += bad;
    }

    /* the default cap refuses a pathological plane outright */
    {
        uint8_t *volc = ARENA_CALLOC(a, vol_sz, 1);
        memset(volc, 0, vol_sz);
        for (int z = 0; z < D; z++) {
            uint8_t *pl = volc + (size_t)z * plane_sz;
            for (int x = 0; x < W; x++) {
                for (int y = 20; y <= 22; y++) pl[ss_idx(y, x, W)] = 1;
                int in_span = (x >= SPAN0 && x < SPAN1);
                int b0 = in_span ? 23 : 26, b1 = in_span ? 25 : 28;
                for (int y = b0; y <= b1; y++) pl[ss_idx(y, x, W)] = 1;
            }
        }
        SlabSplitParams pc = p;
        pc.plane_cap_frac = 0.001f;
        SlabSplitStats stc;
        if (SlabSplit_run(a, volc, NULL, D, H, W, &pc, &stc) != 0 ||
            stc.carved_px != 0 || stc.planes_capped == 0) {
            fprintf(stderr,
                    "[selftest] slab_split: cap did not hold (carved=%lld "
                    "capped=%lld)\n", (long long)stc.carved_px,
                    (long long)stc.planes_capped);
            fails++;
        }
    }

    /* z-support: a single-plane slab must NOT be carved */
    {
        uint8_t *vol2 = ARENA_CALLOC(a, vol_sz, 1);
        for (int z = 0; z < D; z++) {
            uint8_t *pl = vol2 + (size_t)z * plane_sz;
            for (int x = 0; x < W; x++) {
                for (int y = 20; y <= 22; y++) pl[ss_idx(y, x, W)] = 1;
                int fused = (z == 4 && x >= SPAN0 && x < SPAN1);
                int b0 = fused ? 23 : 26, b1 = fused ? 25 : 28;
                for (int y = b0; y <= b1; y++) pl[ss_idx(y, x, W)] = 1;
            }
        }
        SlabSplitStats st2;
        if (SlabSplit_run(a, vol2, NULL, D, H, W, &p, &st2) != 0 ||
            st2.carved_px != 0) {
            fprintf(stderr,
                    "[selftest] slab_split: single-plane slab carved "
                    "(%lld px, want 0)\n", (long long)st2.carved_px);
            fails++;
        }
    }

    /* rounds: a 12-thick slab (4 sheets) sheds more in round 2 */
    {
        uint8_t *vol3 = ARENA_CALLOC(a, vol_sz, 1);
        for (int z = 0; z < D; z++) {
            uint8_t *pl = vol3 + (size_t)z * plane_sz;
            for (int x = 0; x < W; x++)
                for (int y = 20; y <= 31; y++) pl[ss_idx(y, x, W)] = 1;
        }
        SlabSplitParams p2 = p;
        p2.rounds = 2;
        p2.plane_cap_frac = 0.9f;   /* the synthetic plane is ALL slab */
        SlabSplitStats st3;
        if (SlabSplit_run(a, vol3, NULL, D, H, W, &p2, &st3) != 0 ||
            st3.rounds_run < 2) {
            fprintf(stderr,
                    "[selftest] slab_split: 4-sheet slab rounds_run=%d "
                    "(want >= 2)\n", st3.rounds_run);
            fails++;
        }
        int cc = ss_cc_count(a, vol3 + 4 * plane_sz, H, W, 8, W - 8);
        if (cc < 3) {
            fprintf(stderr,
                    "[selftest] slab_split: 4-sheet slab cc=%d (want >= 3)\n",
                    cc);
            fails++;
        }
    }

    Arena_dispose(&a);
    if (fails == 0) fprintf(stderr, "[selftest] slab_split: PASS\n");
    return fails ? -1 : 0;
}
