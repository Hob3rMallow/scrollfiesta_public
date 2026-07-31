/* weld_track.c -- tracked-weld dataset IO. See weld_track.h. */
#include "weld_track.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/ves_platform.h"   /* ves_ensure_parent_dir */

int WeldTrack_write(const WeldTrack *w, const char *path)
{
    assert(w && path);
    ves_ensure_parent_dir(path);
    FILE *f = fopen(path, "wb");
    if (f == NULL)
        return -1;
    int ok = 1;
    ok &= fwrite(WELD_TRACK_MAGIC, 1, 4, f) == 4;
    uint32_t ver = WELD_TRACK_VERSION;
    ok &= fwrite(&ver, sizeof(ver), 1, f) == 1;
    double sp[3] = { w->spiral_a, w->spiral_b, w->pitch };
    ok &= fwrite(sp, sizeof(double), 3, f) == 3;
    ok &= fwrite(w->axis_point, sizeof(float), 3, f) == 3;
    uint64_t nc = (uint64_t)w->n_cubes;
    ok &= fwrite(&nc, sizeof(nc), 1, f) == 1;
    if (w->n_cubes > 0)
        ok &= fwrite(w->ids, WTRK_ID_LEN, w->n_cubes, f) == w->n_cubes;
    uint64_t nk = (uint64_t)w->n_corr;
    ok &= fwrite(&nk, sizeof(nk), 1, f) == 1;
    if (w->n_corr > 0)
        ok &= fwrite(w->corr, sizeof(WtrkCorr), w->n_corr, f) == w->n_corr;
    if (fclose(f) != 0)
        ok = 0;
    return ok ? 0 : -1;
}

int WeldTrack_read(Arena_T arena, const char *path, WeldTrack *out)
{
    assert(arena && path && out);
    memset(out, 0, sizeof(*out));
    FILE *f = fopen(path, "rb");
    if (f == NULL)
        return -1;
    int ok = 1;
    char magic[4];
    ok &= fread(magic, 1, 4, f) == 4;
    if (ok && memcmp(magic, WELD_TRACK_MAGIC, 4) != 0)
        ok = 0;
    uint32_t ver = 0;
    ok &= fread(&ver, sizeof(ver), 1, f) == 1;
    if (ok && ver != WELD_TRACK_VERSION)
        ok = 0;
    double sp[3] = { 0.0, 0.0, 0.0 };
    ok &= fread(sp, sizeof(double), 3, f) == 3;
    ok &= fread(out->axis_point, sizeof(float), 3, f) == 3;
    uint64_t nc = 0;
    ok &= fread(&nc, sizeof(nc), 1, f) == 1;
    if (!ok) {
        fclose(f);
        return -1;
    }
    out->version = ver;
    out->spiral_a = sp[0];
    out->spiral_b = sp[1];
    out->pitch = sp[2];
    out->n_cubes = (size_t)nc;
    if (out->n_cubes > 0) {
        out->ids = (char (*)[WTRK_ID_LEN])ARENA_ALLOC(arena,
                       out->n_cubes * WTRK_ID_LEN);
        ok &= fread(out->ids, WTRK_ID_LEN, out->n_cubes, f) == out->n_cubes;
    }
    uint64_t nk = 0;
    ok &= fread(&nk, sizeof(nk), 1, f) == 1;
    if (ok) {
        out->n_corr = (size_t)nk;
        if (out->n_corr > 0) {
            out->corr = (WtrkCorr *)ARENA_ALLOC(arena,
                            out->n_corr * sizeof(WtrkCorr));
            ok &= fread(out->corr, sizeof(WtrkCorr), out->n_corr, f)
                  == out->n_corr;
        }
    }
    fclose(f);
    if (!ok)
        return -1;
    out->self = out;
    return 0;
}

int32_t WeldTrack_cube_index(const WeldTrack *w, const char *id)
{
    assert(w && id);
    for (size_t i = 0; i < w->n_cubes; i++)
        if (strncmp(w->ids[i], id, WTRK_ID_LEN) == 0)
            return (int32_t)i;
    return -1;
}

int WeldTrack_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    const char *path = "output/_selftest_weld_track/wt.wtrk";

    WeldTrack w;
    memset(&w, 0, sizeof(w));
    w.version = WELD_TRACK_VERSION;
    w.spiral_a = 0.5;
    w.spiral_b = -9.5;
    w.pitch = 9.5;
    w.axis_point[0] = 0.0f;
    w.axis_point[1] = 3405.0f;
    w.axis_point[2] = 2878.0f;
    char ids[3][WTRK_ID_LEN];
    memset(ids, 0, sizeof(ids));
    snprintf(ids[0], WTRK_ID_LEN, "z04352_y03072_x02560");
    snprintf(ids[1], WTRK_ID_LEN, "z04352_y03072_x02688");
    snprintf(ids[2], WTRK_ID_LEN, "z04480_y03072_x02560");
    w.n_cubes = 3;
    w.ids = ids;
    WtrkCorr cs[2];
    memset(cs, 0, sizeof(cs));
    cs[0].cube_a = 0; cs[0].cube_b = 1; cs[0].vert_a = 10; cs[0].vert_b = 20;
    cs[0].gid_a = 0; cs[0].gid_b = 2; cs[0].dr = 0.3f; cs[0].dphi = 0.01f;
    cs[0].dist3d = 1.2f; cs[0].seam_axis = 2; cs[0].conf = 2;
    cs[1].cube_a = 0; cs[1].cube_b = 2; cs[1].vert_a = 5; cs[1].vert_b = 7;
    cs[1].gid_a = 1; cs[1].gid_b = 1; cs[1].dr = 0.1f; cs[1].dphi = -0.02f;
    cs[1].dist3d = 0.9f; cs[1].seam_axis = 0; cs[1].conf = 1;
    w.n_corr = 2;
    w.corr = cs;

    if (WeldTrack_write(&w, path) != 0) {
        fprintf(stderr, "[weld_track selftest]   FAIL: write\n");
        fails++;
    } else {
        WeldTrack r;
        if (WeldTrack_read(arena, path, &r) != 0) {
            fprintf(stderr, "[weld_track selftest]   FAIL: read\n");
            fails++;
        } else {
            int ok = r.version == WELD_TRACK_VERSION
                  && r.n_cubes == 3 && r.n_corr == 2
                  && fabs(r.spiral_b - (-9.5)) < 1e-9
                  && fabs((double)r.axis_point[1] - 3405.0) < 1e-4
                  && strncmp(r.ids[2], "z04480_y03072_x02560", WTRK_ID_LEN) == 0
                  && r.corr[0].gid_b == 2 && r.corr[0].seam_axis == 2
                  && r.corr[0].conf == 2 && r.corr[1].seam_axis == 0
                  && r.corr[1].vert_b == 7 && r.self == &r;
            if (!ok) {
                fprintf(stderr, "[weld_track selftest]   FAIL: roundtrip "
                        "mismatch\n");
                fails++;
            }
            if (WeldTrack_cube_index(&r, "z04352_y03072_x02688") != 1
                || WeldTrack_cube_index(&r, "nope") != -1) {
                fprintf(stderr, "[weld_track selftest]   FAIL: cube_index\n");
                fails++;
            }
        }
    }
    Arena_dispose(&arena);
    fprintf(stderr, "[weld_track selftest] %s (%d failure%s)\n",
            fails == 0 ? "PASSED" : "FAILED", fails, fails == 1 ? "" : "s");
    return fails;
}
