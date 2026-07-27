/*
 * pred_reject.c -- detect garbage (solid-slab) prediction cubes.
 *
 * See pred_reject.h for the rationale. The verdict is the conjunction of a 3D
 * thickness signal (erosion-survival) and a per-frame solid-rectangle signal,
 * so a dense-but-thin real tangle (fails thickness) and a flat sheet seen
 * face-on (fails the rectangle-run) are both KEPT; only a true solid slab,
 * which is thick in 3D AND a filled rectangle over a long run of frames, is
 * rejected.
 */
#include "pred_reject.h"
#include "../common/pipeline_constants.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- 3D 6-connected erosion: dst[i] = src[i] && all 6 face-neighbours set.
 * Out-of-bounds neighbours count as background, so faces erode inward.
 * Returns the number of surviving voxels. ---- */
static size_t erode6(const uint8_t *src, uint8_t *dst, int D, int H, int W)
{
    size_t HW = (size_t)H * (size_t)W;
    size_t survivors = 0;
    for (int z = 0; z < D; z++) {
        for (int y = 0; y < H; y++) {
            for (int x = 0; x < W; x++) {
                size_t i = (size_t)z * HW + (size_t)y * (size_t)W + (size_t)x;
                uint8_t keep = src[i];
                if (keep) {
                    if (z == 0 || z == D - 1 || y == 0 || y == H - 1 ||
                        x == 0 || x == W - 1) {
                        keep = 0;
                    } else {
                        keep = (uint8_t)(src[i - HW] && src[i + HW] &&
                                         src[i - (size_t)W] && src[i + (size_t)W] &&
                                         src[i - 1] && src[i + 1]);
                    }
                }
                dst[i] = keep;
                if (keep) survivors++;
            }
        }
    }
    return survivors;
}

/* ---- Largest 4-connected FG component in a da*db binary slice (slc[a*db+b]).
 * Writes its pixel area and the fraction of its bounding box it fills (extent).
 * `stack` and `visited` are caller scratch, each >= da*db elements. ---- */
static void slice_largest_rect(const uint8_t *slc, int da, int db,
                               int32_t *stack, uint8_t *visited,
                               int *out_area, double *out_bbox_fill)
{
    size_t sa = (size_t)da * (size_t)db;
    memset(visited, 0, sa);
    int best_area = 0;
    double best_fill = 0.0;

    for (int a0 = 0; a0 < da; a0++) {
        for (int b0 = 0; b0 < db; b0++) {
            size_t s0 = (size_t)a0 * (size_t)db + (size_t)b0;
            if (!slc[s0] || visited[s0]) continue;

            int top = 0;
            stack[top++] = (int32_t)s0;
            visited[s0] = 1;
            int area = 0;
            int amin = a0, amax = a0, bmin = b0, bmax = b0;

            while (top > 0) {
                int32_t s = stack[--top];
                int a = (int)(s / db);
                int b = (int)(s % db);
                area++;
                if (a < amin) amin = a;
                if (a > amax) amax = a;
                if (b < bmin) bmin = b;
                if (b > bmax) bmax = b;

                if (a > 0) {
                    size_t t = (size_t)s - (size_t)db;
                    if (slc[t] && !visited[t]) { visited[t] = 1; stack[top++] = (int32_t)t; }
                }
                if (a < da - 1) {
                    size_t t = (size_t)s + (size_t)db;
                    if (slc[t] && !visited[t]) { visited[t] = 1; stack[top++] = (int32_t)t; }
                }
                if (b > 0) {
                    size_t t = (size_t)s - 1;
                    if (slc[t] && !visited[t]) { visited[t] = 1; stack[top++] = (int32_t)t; }
                }
                if (b < db - 1) {
                    size_t t = (size_t)s + 1;
                    if (slc[t] && !visited[t]) { visited[t] = 1; stack[top++] = (int32_t)t; }
                }
            }

            if (area > best_area) {
                best_area = area;
                double bbox = (double)(amax - amin + 1) * (double)(bmax - bmin + 1);
                best_fill = (bbox > 0.0) ? (double)area / bbox : 0.0;
            }
        }
    }
    *out_area = best_area;
    *out_bbox_fill = best_fill;
}

/* ---- Scan one axis' slices for solid-rectangle frames. `axis` selects the
 * slicing: 0 -> (y,x) slices indexed by z, 1 -> (z,x) by y, 2 -> (z,y) by x.
 * Fills the fraction of rect frames and the longest consecutive run. ---- */
static void axis_rect_scan(const uint8_t *bin, int D, int H, int W, int axis,
                           uint8_t *slc, int32_t *stack, uint8_t *visited,
                           double *out_rect_frac, int *out_max_run)
{
    size_t HW = (size_t)H * (size_t)W;
    int n_slices = 0, da = 0, db = 0;
    if (axis == 0)      { n_slices = D; da = H; db = W; }
    else if (axis == 1) { n_slices = H; da = D; db = W; }
    else                { n_slices = W; da = D; db = H; }

    double slice_area = (double)da * (double)db;
    int rect_frames = 0, run = 0, max_run = 0;

    for (int s = 0; s < n_slices; s++) {
        /* Extract slice s into the contiguous da*db buffer slc[a*db + b]. */
        for (int a = 0; a < da; a++) {
            for (int b = 0; b < db; b++) {
                size_t vi;
                if (axis == 0)      vi = (size_t)s * HW + (size_t)a * (size_t)W + (size_t)b;
                else if (axis == 1) vi = (size_t)a * HW + (size_t)s * (size_t)W + (size_t)b;
                else                vi = (size_t)a * HW + (size_t)b * (size_t)W + (size_t)s;
                slc[(size_t)a * (size_t)db + (size_t)b] = bin[vi];
            }
        }

        int area = 0;
        double fill = 0.0;
        slice_largest_rect(slc, da, db, stack, visited, &area, &fill);
        double area_frac = (slice_area > 0.0) ? (double)area / slice_area : 0.0;

        int is_rect = (area_frac >= (double)GARBAGE_RECT_AREA_FRAC) &&
                      (fill >= (double)GARBAGE_RECT_FILL);
        if (is_rect) {
            rect_frames++;
            run++;
            if (run > max_run) max_run = run;
        } else {
            run = 0;
        }
    }

    *out_rect_frac = (n_slices > 0) ? (double)rect_frames / (double)n_slices : 0.0;
    *out_max_run = max_run;
}

/* ---- Largest 6-connected foreground component size (voxels). `bin` is 0/1.
 * Iterative flood fill; caller scratch: `visited` (>= n uint8) and `stack`
 * (>= n int32). n = D*H*W fits int32 for any cube up to ~1290^3, far beyond
 * the 128^3 / halo sizes used here. ---- */
static size_t largest_cc6(const uint8_t *bin, int D, int H, int W,
                          uint8_t *visited, int32_t *stack)
{
    size_t HW = (size_t)H * (size_t)W;
    size_t n  = (size_t)D * HW;
    memset(visited, 0, n);
    size_t best = 0;
    for (size_t s = 0; s < n; s++) {
        if (!bin[s] || visited[s]) continue;
        size_t sp = 0, count = 0;
        stack[sp++] = (int32_t)s;
        visited[s] = 1;
        while (sp > 0) {
            size_t i = (size_t)stack[--sp];
            count++;
            int z = (int)(i / HW);
            int rem = (int)(i % HW);
            int y = rem / W, x = rem % W;
            size_t j;
            if (z > 0)     { j = i - HW; if (bin[j] && !visited[j]) { visited[j] = 1; stack[sp++] = (int32_t)j; } }
            if (z < D - 1) { j = i + HW; if (bin[j] && !visited[j]) { visited[j] = 1; stack[sp++] = (int32_t)j; } }
            if (y > 0)     { j = i - (size_t)W; if (bin[j] && !visited[j]) { visited[j] = 1; stack[sp++] = (int32_t)j; } }
            if (y < H - 1) { j = i + (size_t)W; if (bin[j] && !visited[j]) { visited[j] = 1; stack[sp++] = (int32_t)j; } }
            if (x > 0)     { j = i - 1; if (bin[j] && !visited[j]) { visited[j] = 1; stack[sp++] = (int32_t)j; } }
            if (x < W - 1) { j = i + 1; if (bin[j] && !visited[j]) { visited[j] = 1; stack[sp++] = (int32_t)j; } }
        }
        if (count > best) best = count;
    }
    return best;
}

int PredReject_is_garbage(Arena_T arena, const uint8_t *vol,
                          int D, int H, int W, PredRejectStats *out)
{
    assert(arena && vol && D > 0 && H > 0 && W > 0);

    PredRejectStats st;
    memset(&st, 0, sizeof(st));
    st.n_vox = (size_t)D * (size_t)H * (size_t)W;
    st.reason = "keep";

    Arena_Mark mark = Arena_save(arena);

    /* Local binary copy (vol is const; any nonzero -> 1). */
    uint8_t *bin = (uint8_t *)ARENA_ALLOC(arena, (long)st.n_vox);
    size_t fg = 0;
    for (size_t i = 0; i < st.n_vox; i++) {
        uint8_t v = (uint8_t)(vol[i] ? 1 : 0);
        bin[i] = v;
        fg += (size_t)v;
    }
    st.fg_count = fg;
    st.fill_frac = (st.n_vox > 0) ? (double)fg / (double)st.n_vox : 0.0;

    /* ---- EMPTY garbage: no 6-connected component big enough to mesh. Step 0
     * keeps only components with >= MIN_CC_SIZE voxels, so if the LARGEST
     * component is below that floor, Step 0 produces nothing and meshing is a
     * certain empty FAIL. Reject pre-spawn as a sibling of the solid-slab
     * reject. This is the exact Step-0 criterion -- stricter than a total-FG
     * test, since sub-threshold scatter can sum past MIN_CC_SIZE while no single
     * component is meshable -- so it never drops a real sheet (a genuine patch
     * forms one component well above the floor). ---- */
    {
        Arena_Mark cc_mark = Arena_save(arena);
        uint8_t *cc_vis = (uint8_t *)ARENA_ALLOC(arena, (long)st.n_vox);
        int32_t *cc_stk = (int32_t *)ARENA_ALLOC(arena, (long)(st.n_vox * sizeof(int32_t)));
        st.largest_cc = largest_cc6(bin, D, H, W, cc_vis, cc_stk);
        Arena_restore(arena, cc_mark);
    }
    if (st.largest_cc < (size_t)MIN_CC_SIZE) {
        st.verdict = 1;
        st.empty = 1;
        st.reason = "reject: empty (no meshable component)";
        if (out) *out = st;
        Arena_restore(arena, mark);
        return 1;
    }

    /* ---- Signal 1: 3D erosion-survival (thickness). ---- */
    uint8_t *cur = (uint8_t *)ARENA_ALLOC(arena, (long)st.n_vox);
    uint8_t *nxt = (uint8_t *)ARENA_ALLOC(arena, (long)st.n_vox);
    memcpy(cur, bin, st.n_vox);
    size_t cur_count = fg;
    st.max_thickness = 1;            /* depth 1 = "is foreground" */
    st.interior_count = 0;
    for (int pass = 1; pass <= GARBAGE_ERODE_MAXPASS && cur_count > 0; pass++) {
        size_t nc = erode6(cur, nxt, D, H, W);
        if (nc > 0) st.max_thickness = pass + 1;
        if (pass == GARBAGE_ERODE_R) st.interior_count = nc;
        uint8_t *tmp = cur; cur = nxt; nxt = tmp;
        cur_count = nc;
    }
    st.interior_frac = (double)st.interior_count / (double)fg;

    /* ---- Signal 2: per-frame solid rectangle (strongest of the 3 axes). ---- */
    /* One scratch buffer set sized to the largest slice across all axes. */
    size_t max_slice = (size_t)H * (size_t)W;
    if ((size_t)D * (size_t)W > max_slice) max_slice = (size_t)D * (size_t)W;
    if ((size_t)D * (size_t)H > max_slice) max_slice = (size_t)D * (size_t)H;
    uint8_t *slc     = (uint8_t *)ARENA_ALLOC(arena, (long)max_slice);
    uint8_t *visited = (uint8_t *)ARENA_ALLOC(arena, (long)max_slice);
    int32_t *stack   = (int32_t *)ARENA_ALLOC(arena, (long)(max_slice * sizeof(int32_t)));

    int any_rect_pass = 0;
    st.rect_axis = 0;
    st.rect_frac = 0.0;
    st.max_rect_run = -1;
    for (int axis = 0; axis < 3; axis++) {
        double rfrac = 0.0;
        int mrun = 0;
        axis_rect_scan(bin, D, H, W, axis, slc, stack, visited, &rfrac, &mrun);
        if (rfrac >= (double)GARBAGE_RECT_FRAC && mrun >= GARBAGE_RECT_RUN) {
            any_rect_pass = 1;
        }
        /* Report the axis with the longest rectangle run. */
        if (mrun > st.max_rect_run) {
            st.max_rect_run = mrun;
            st.rect_frac = rfrac;
            st.rect_axis = axis;
        }
    }
    if (st.max_rect_run < 0) st.max_rect_run = 0;

    /* ---- Conservative AND decision. ---- */
    int thick = (st.interior_frac >= (double)GARBAGE_INTERIOR_FRAC) &&
                (st.interior_count >= (size_t)GARBAGE_INTERIOR_MIN);
    if (thick && any_rect_pass) {
        st.verdict = 1;
        st.reason = "reject: solid slab";
    } else if (!thick) {
        st.verdict = 0;
        st.reason = "keep: thin (low erosion interior)";
    } else {
        st.verdict = 0;
        st.reason = "keep: not a persistent rectangle";
    }

    if (out) *out = st;
    Arena_restore(arena, mark);
    return st.verdict;
}

/* =====================================================================
 * Self-test: synthetic volumes with KNOWN structure.
 * ===================================================================== */
static void fill_block(uint8_t *v, int D, int H, int W,
                       int z0, int z1, int y0, int y1, int x0, int x1)
{
    (void)D;
    size_t HW = (size_t)H * (size_t)W;
    for (int z = z0; z < z1; z++)
        for (int y = y0; y < y1; y++)
            for (int x = x0; x < x1; x++)
                v[(size_t)z * HW + (size_t)y * (size_t)W + (size_t)x] = 255;
}

int PredReject_selftest(void)
{
    int fails = 0;
    const int N = 64;
    size_t nv = (size_t)N * N * N;
    printf("=== pred_reject selftest ===\n");

    Arena_T arena = Arena_new();
    uint8_t *v = (uint8_t *)malloc(nv);
    if (!v) { printf("  FAIL: malloc\n"); Arena_dispose(&arena); return 3; }

    /* case 1: solid 40^3 slab -> REJECT */
    {
        memset(v, 0, nv);
        fill_block(v, N, N, N, 12, 52, 12, 52, 12, 52);
        PredRejectStats s;
        int r = PredReject_is_garbage(arena, v, N, N, N, &s);
        printf("[case1] solid 40^3 slab: verdict=%d interior_frac=%.3f run=%d (%s)\n",
               r, s.interior_frac, s.max_rect_run, s.reason);
        if (r != 1) { printf("  FAIL: expected REJECT\n"); fails++; }
    }

    /* case 2: two thin (2-vox) curved-ish sheets -> KEEP */
    {
        memset(v, 0, nv);
        /* sheet A: y in {20,21} across all z,x (a thin flat sheet) */
        fill_block(v, N, N, N, 0, N, 20, 22, 0, N);
        /* sheet B: x in {40,41} across all z,y */
        fill_block(v, N, N, N, 0, N, 0, N, 40, 42);
        PredRejectStats s;
        int r = PredReject_is_garbage(arena, v, N, N, N, &s);
        printf("[case2] two 2-vox sheets: verdict=%d interior_frac=%.3f run=%d (%s)\n",
               r, s.interior_frac, s.max_rect_run, s.reason);
        if (r != 0) { printf("  FAIL: expected KEEP\n"); fails++; }
    }

    /* case 3: flat patch seen face-on (50x50, 2 z-slices thick) -> KEEP */
    {
        memset(v, 0, nv);
        fill_block(v, N, N, N, 30, 32, 7, 57, 7, 57);
        PredRejectStats s;
        int r = PredReject_is_garbage(arena, v, N, N, N, &s);
        printf("[case3] flat face-on patch: verdict=%d interior_frac=%.3f run=%d (%s)\n",
               r, s.interior_frac, s.max_rect_run, s.reason);
        if (r != 0) { printf("  FAIL: expected KEEP\n"); fails++; }
    }

    /* case 4: solid slab + thin protrusion ("rectangle with a bit sticking up") -> REJECT */
    {
        memset(v, 0, nv);
        fill_block(v, N, N, N, 10, 50, 10, 50, 10, 45);   /* solid block */
        fill_block(v, N, N, N, 28, 31, 28, 31, 45, 60);   /* spike sticking out */
        PredRejectStats s;
        int r = PredReject_is_garbage(arena, v, N, N, N, &s);
        printf("[case4] slab + protrusion: verdict=%d interior_frac=%.3f run=%d (%s)\n",
               r, s.interior_frac, s.max_rect_run, s.reason);
        if (r != 1) { printf("  FAIL: expected REJECT\n"); fails++; }
    }

    /* case 5: empty volume -> REJECT (empty garbage) */
    {
        memset(v, 0, nv);
        PredRejectStats s;
        int r = PredReject_is_garbage(arena, v, N, N, N, &s);
        printf("[case5] empty: verdict=%d empty=%d (%s)\n", r, s.empty, s.reason);
        if (r != 1 || s.empty != 1) { printf("  FAIL: expected REJECT-empty\n"); fails++; }
    }

    /* case 6: a tiny speck (7^3 = 343 FG < MIN_CC_SIZE) -> REJECT (empty) */
    {
        memset(v, 0, nv);
        fill_block(v, N, N, N, 20, 27, 20, 27, 20, 27);
        PredRejectStats s;
        int r = PredReject_is_garbage(arena, v, N, N, N, &s);
        printf("[case6] tiny speck fg=%zu: verdict=%d empty=%d (%s)\n",
               s.fg_count, r, s.empty, s.reason);
        if (r != 1 || s.empty != 1) { printf("  FAIL: expected REJECT-empty\n"); fails++; }
    }

    /* case 7: a small 1-vox-thick sheet just above the floor (25x25 = 625 FG
     * >= MIN_CC_SIZE) -> KEEP. Guards the empty floor from over-rejecting a
     * genuine, if small, thin sheet. */
    {
        memset(v, 0, nv);
        fill_block(v, N, N, N, 30, 31, 10, 35, 10, 35);
        PredRejectStats s;
        int r = PredReject_is_garbage(arena, v, N, N, N, &s);
        printf("[case7] small sheet fg=%zu cc=%zu: verdict=%d empty=%d (%s)\n",
               s.fg_count, s.largest_cc, r, s.empty, s.reason);
        if (r != 0) { printf("  FAIL: expected KEEP\n"); fails++; }
    }

    /* case 8: scattered sub-threshold specks -- total FG (3x343 = 1029) exceeds
     * MIN_CC_SIZE but every 6-conn component (343) is below it, so Step 0 meshes
     * nothing. REJECT (empty). This is the real z04352_y04480_x02816 case and
     * guards against a total-FG-only test passing it through to a mesh FAIL. */
    {
        memset(v, 0, nv);
        fill_block(v, N, N, N,  4, 11,  4, 11,  4, 11);
        fill_block(v, N, N, N,  4, 11,  4, 11, 40, 47);
        fill_block(v, N, N, N, 40, 47, 40, 47, 40, 47);
        PredRejectStats s;
        int r = PredReject_is_garbage(arena, v, N, N, N, &s);
        printf("[case8] scatter fg=%zu cc=%zu: verdict=%d empty=%d (%s)\n",
               s.fg_count, s.largest_cc, r, s.empty, s.reason);
        if (r != 1 || s.empty != 1) { printf("  FAIL: expected REJECT-empty\n"); fails++; }
    }

    free(v);
    Arena_dispose(&arena);
    printf("=== selftest %s (%d failure%s) ===\n",
           fails ? "FAILED" : "PASSED", fails, fails == 1 ? "" : "s");
    return fails ? 3 : 0;
}
