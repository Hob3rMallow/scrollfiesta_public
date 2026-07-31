/* piece_set.c -- concat the placed dir into one flat mesh. See piece_set.h. */
#include "../common/ves_platform.h"

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _MSC_VER
#include <windows.h>
#else
#include <dirent.h>
#endif

#include "../common/obj_io.h"
#include "../common/raw_sample.h"   /* vertex_normals */
#include "piece_set.h"

/* ---- id enumeration (from *_uvphi.f32) ------------------------------------ */

typedef struct { char id[48]; } PsId;

static int ps_scan_ids(Arena_T arena, const char *dir, long z_lo, long z_hi,
                       PsId **out, size_t *out_n)
{
    PsId *ids = NULL;
    size_t n = 0, cap = 0;
    const char *suffix = "_uvphi.f32";
    const size_t slen = strlen(suffix);
#ifdef _MSC_VER
    char pattern[1024];
    snprintf(pattern, sizeof(pattern), "%s/*_uvphi.f32", dir);
    WIN32_FIND_DATAA fd;
    HANDLE h = FindFirstFileA(pattern, &fd);
    if (h == INVALID_HANDLE_VALUE) return -1;
    do {
        size_t len = strlen(fd.cFileName);
        if (len <= slen || fd.cFileName[0] != 'z') continue;
        if (len - slen >= sizeof(ids[0].id)) continue;
        {   /* z-slab filter: id is "z#####_..."; keep z in [z_lo, z_hi) */
            long zo = strtol(fd.cFileName + 1, NULL, 10);
            if (zo < z_lo || zo >= z_hi) continue;
        }
        if (n >= cap) {
            size_t nc = cap == 0 ? 64 : cap * 2;
            PsId *ni = (PsId *)ARENA_ALLOC(arena, nc * sizeof(PsId));
            if (ids) memcpy(ni, ids, n * sizeof(PsId));
            ids = ni;
            cap = nc;
        }
        memcpy(ids[n].id, fd.cFileName, len - slen);
        ids[n].id[len - slen] = '\0';
        n++;
    } while (FindNextFileA(h, &fd));
    FindClose(h);
#else
    DIR *d = opendir(dir);
    if (d == NULL) return -1;
    struct dirent *ent = NULL;
    while ((ent = readdir(d)) != NULL) {
        size_t len = strlen(ent->d_name);
        if (len <= slen || ent->d_name[0] != 'z') continue;
        if (strcmp(ent->d_name + len - slen, suffix) != 0) continue;
        if (len - slen >= sizeof(ids[0].id)) continue;
        {   /* z-slab filter: id is "z#####_..."; keep z in [z_lo, z_hi) */
            long zo = strtol(ent->d_name + 1, NULL, 10);
            if (zo < z_lo || zo >= z_hi) continue;
        }
        if (n >= cap) {
            size_t nc = cap == 0 ? 64 : cap * 2;
            PsId *ni = (PsId *)ARENA_ALLOC(arena, nc * sizeof(PsId));
            if (ids) memcpy(ni, ids, n * sizeof(PsId));
            ids = ni;
            cap = nc;
        }
        memcpy(ids[n].id, ent->d_name, len - slen);
        ids[n].id[len - slen] = '\0';
        n++;
    }
    closedir(d);
#endif
    *out = ids;
    *out_n = n;
    return n > 0 ? 0 : -1;
}

static void *ps_read_blob(Arena_T arena, const char *path, size_t *out_nbytes)
{
    *out_nbytes = 0;
    FILE *f = fopen(path, "rb");
    if (f == NULL) return NULL;
    if (fseek(f, 0, SEEK_END) != 0) { fclose(f); return NULL; }
    long sz = ftell(f);
    if (sz < 0) { fclose(f); return NULL; }
    rewind(f);
    void *buf = ARENA_ALLOC(arena, (size_t)sz + 1);
    if (sz > 0 && fread(buf, 1, (size_t)sz, f) != (size_t)sz) {
        fclose(f);
        return NULL;
    }
    fclose(f);
    *out_nbytes = (size_t)sz;
    return buf;
}

int PieceSet_build(Arena_T arena, const char *placed_dir, PieceSet *out)
{
    return PieceSet_build_z(arena, placed_dir, LONG_MIN, LONG_MAX, out);
}

int PieceSet_build_z(Arena_T arena, const char *placed_dir,
                     long z_lo, long z_hi, PieceSet *out)
{
    assert(arena && placed_dir && out);
    memset(out, 0, sizeof(*out));
    out->u_min = out->v_min = 1e300;
    out->u_max = out->v_max = -1e300;

    PsId *ids = NULL;
    size_t nid = 0;
    if (ps_scan_ids(arena, placed_dir, z_lo, z_hi, &ids, &nid) != 0) return -1;

    /* pass 1: total sizes from the sidecar files (facekeep counts kept faces) */
    size_t tot_nv = 0, tot_nf = 0, n_used = 0;
    {
        Arena_Mark mark = Arena_save(arena);
        for (size_t c = 0; c < nid; c++) {
            char path[1024];
            size_t ub = 0, kb = 0;
            snprintf(path, sizeof(path), "%s/%s_uvphi.f32", placed_dir, ids[c].id);
            FILE *f = fopen(path, "rb");
            if (f == NULL) continue;
            fseek(f, 0, SEEK_END);
            long usz = ftell(f);
            fclose(f);
            if (usz <= 0 || (size_t)usz % (3 * sizeof(float)) != 0) continue;
            ub = (size_t)usz / (3 * sizeof(float));
            snprintf(path, sizeof(path), "%s/%s_facekeep.u8", placed_dir, ids[c].id);
            size_t nb = 0;
            uint8_t *keep = (uint8_t *)ps_read_blob(arena, path, &nb);
            if (keep == NULL) continue;
            for (size_t t = 0; t < nb; t++) if (keep[t]) kb++;
            tot_nv += ub;
            tot_nf += kb;
            n_used++;
        }
        Arena_restore(arena, mark);
    }
    if (n_used == 0 || tot_nv == 0 || tot_nf == 0) return -1;

    out->verts = (float *)ARENA_ALLOC(arena, tot_nv * 3 * sizeof(float));
    out->uv = (float *)ARENA_ALLOC(arena, tot_nv * 2 * sizeof(float));
    out->phi = (float *)ARENA_ALLOC(arena, tot_nv * sizeof(float));
    out->normals = (float *)ARENA_ALLOC(arena, tot_nv * 3 * sizeof(float));
    out->faces = (int32_t *)ARENA_ALLOC(arena, tot_nf * 3 * sizeof(int32_t));
    out->face_cube = (int32_t *)ARENA_ALLOC(arena, tot_nf * sizeof(int32_t));
    out->gid = (int32_t *)ARENA_ALLOC(arena, tot_nv * sizeof(int32_t));
    out->ids = (char (*)[48])ARENA_ALLOC(arena, nid * sizeof(*out->ids));
    out->cube_voff = (size_t *)ARENA_ALLOC(arena, (nid + 1) * sizeof(size_t));
    out->cube_org = (long (*)[3])ARENA_ALLOC(arena, nid * sizeof(*out->cube_org));

    /* pass 2: load + renumber */
    size_t voff = 0, foff = 0;
    for (size_t c = 0; c < nid; c++) {
        Arena_Mark mark = Arena_save(arena);
        char path[1024];
        size_t nb = 0;
        snprintf(path, sizeof(path), "%s/%s_uvphi.f32", placed_dir, ids[c].id);
        float *uvphi = (float *)ps_read_blob(arena, path, &nb);
        if (uvphi == NULL || nb % (3 * sizeof(float)) != 0 || nb == 0) {
            Arena_restore(arena, mark);
            continue;
        }
        size_t nv = nb / (3 * sizeof(float));

        snprintf(path, sizeof(path), "%s/%s_facekeep.u8", placed_dir, ids[c].id);
        size_t kb = 0;
        uint8_t *keep = (uint8_t *)ps_read_blob(arena, path, &kb);
        if (keep == NULL) { Arena_restore(arena, mark); continue; }

        float *mv = NULL;
        int32_t *mf = NULL;
        size_t mnv = 0, mnf = 0;
        snprintf(path, sizeof(path), "%s/%s_mesh.obj", placed_dir, ids[c].id);
        if (ObjIO_read(arena, path, &mv, &mnv, &mf, &mnf) != 0 ||
            mnv != nv || mnf != kb) {
            fprintf(stderr, "piece_set: WARN %s mesh/sidecar mismatch "
                    "(nv %zu vs %zu, nf %zu vs %zu) -- skipped\n",
                    ids[c].id, mnv, nv, mnf, kb);
            Arena_restore(arena, mark);
            continue;
        }
        if (voff + nv > tot_nv) { Arena_restore(arena, mark); break; }

        /* kept faces, renumbered; provenance = this cube */
        size_t nkept = 0;
        for (size_t t = 0; t < mnf; t++) {
            if (!keep[t]) continue;
            if (foff + nkept >= tot_nf) break;
            out->faces[(foff + nkept) * 3 + 0] = (int32_t)((size_t)mf[t * 3 + 0] + voff);
            out->faces[(foff + nkept) * 3 + 1] = (int32_t)((size_t)mf[t * 3 + 1] + voff);
            out->faces[(foff + nkept) * 3 + 2] = (int32_t)((size_t)mf[t * 3 + 2] + voff);
            out->face_cube[foff + nkept] = (int32_t)out->n_cubes;
            nkept++;
        }

        /* per-cube normals over the KEPT faces (bad links would skew rims) */
        {
            int32_t *kf = (int32_t *)ARENA_ALLOC(arena,
                                                 (nkept ? nkept : 1) * 3 *
                                                 sizeof(int32_t));
            size_t q = 0;
            for (size_t t = 0; t < mnf; t++) {
                if (!keep[t]) continue;
                kf[q * 3 + 0] = mf[t * 3 + 0];
                kf[q * 3 + 1] = mf[t * 3 + 1];
                kf[q * 3 + 2] = mf[t * 3 + 2];
                q++;
                if (q >= nkept) break;
            }
            float *nrm = vertex_normals(mv, nv, kf, nkept);   /* malloc'd */
            if (nrm != NULL) {
                memcpy(&out->normals[voff * 3], nrm, nv * 3 * sizeof(float));
                free(nrm);
            } else {
                memset(&out->normals[voff * 3], 0, nv * 3 * sizeof(float));
            }
        }

        memcpy(&out->verts[voff * 3], mv, nv * 3 * sizeof(float));
        for (size_t i = 0; i < nv; i++) {
            out->uv[(voff + i) * 2 + 0] = uvphi[i * 3 + 0];
            out->uv[(voff + i) * 2 + 1] = uvphi[i * 3 + 1];
            out->phi[voff + i] = uvphi[i * 3 + 2];
        }

        /* winding-group sidecar (optional; -1 fill when absent/mismatched).
         * Loaded into the persistent range BEFORE the scratch restore. */
        {
            for (size_t i = 0; i < nv; i++)
                out->gid[voff + i] = -1;
            char gp[1024];
            snprintf(gp, sizeof(gp), "%s/%s_group.i32", placed_dir, ids[c].id);
            size_t gb = 0;
            int32_t *g = (int32_t *)ps_read_blob(arena, gp, &gb);
            if (g != NULL && gb == nv * sizeof(int32_t))
                memcpy(&out->gid[voff], g, nv * sizeof(int32_t));
        }

        snprintf(out->ids[out->n_cubes], sizeof(out->ids[0]), "%s", ids[c].id);
        {
            long oz = 0, oy = 0, ox = 0;
            if (sscanf(ids[c].id, "z%05ld_y%05ld_x%05ld", &oz, &oy, &ox) != 3) {
                oz = oy = ox = -1;   /* non-canonical id (selftests) */
            }
            out->cube_org[out->n_cubes][0] = oz;
            out->cube_org[out->n_cubes][1] = oy;
            out->cube_org[out->n_cubes][2] = ox;
        }
        out->cube_voff[out->n_cubes] = voff;
        out->n_cubes++;
        out->cube_voff[out->n_cubes] = voff + nv;
        voff += nv;
        foff += nkept;
        Arena_restore(arena, mark);
    }
    out->nv = voff;
    out->nf = foff;

    /* UV extents over vertices REFERENCED by kept faces (unmapped verts carry
     * u=0 which would drag the canvas origin) */
    for (size_t t = 0; t < out->nf; t++) {
        for (int e = 0; e < 3; e++) {
            size_t i = (size_t)out->faces[t * 3 + e];
            double u = (double)out->uv[i * 2 + 0], v = (double)out->uv[i * 2 + 1];
            if (u < out->u_min) out->u_min = u;
            if (u > out->u_max) out->u_max = u;
            if (v < out->v_min) out->v_min = v;
            if (v > out->v_max) out->v_max = v;
        }
    }
    return out->nf > 0 ? 0 : -1;
}

size_t PieceSet_apply_facekeep(PieceSet *ps, const uint8_t *keep)
{
    assert(ps && keep);
    size_t w = 0;
    for (size_t f = 0; f < ps->nf; f++) {
        if (!keep[f]) continue;
        if (w != f) {
            ps->faces[w * 3 + 0] = ps->faces[f * 3 + 0];
            ps->faces[w * 3 + 1] = ps->faces[f * 3 + 1];
            ps->faces[w * 3 + 2] = ps->faces[f * 3 + 2];
            ps->face_cube[w] = ps->face_cube[f];
        }
        w++;
    }
    ps->nf = w;
    return w;
}

void PieceSet_refresh_normals(PieceSet *ps)
{
    assert(ps);
    if (ps->nv == 0)
        return;
    float *nrm = vertex_normals(ps->verts, ps->nv, ps->faces, ps->nf);
    if (nrm != NULL) {
        memcpy(ps->normals, nrm, ps->nv * 3 * sizeof(float));
        free(nrm);
    }
}

/* ============================================================================
 * Self-test.
 * ==========================================================================*/

static int ps_check(int cond, const char *what, int *fails)
{
    if (!cond) {
        fprintf(stderr, "[piece_set selftest]   FAIL: %s\n", what);
        (*fails)++;
    }
    return cond;
}

static int ps_write_blob(const char *path, const void *data, size_t nbytes)
{
    FILE *f = fopen(path, "wb");
    if (f == NULL) return -1;
    size_t w = nbytes ? fwrite(data, 1, nbytes, f) : 0;
    fclose(f);
    return w == nbytes ? 0 : -1;
}

int PieceSet_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    const char *dir = "output/_selftest_scroll_unroll";
    char path[1024];
    snprintf(path, sizeof(path), "%s/x", dir);
    ves_ensure_parent_dir(path);

    /* two fake cubes: 4 verts / 2 faces each; cube B drops face 1 */
    const char *cid[2] = { "z00000_y00000_x00000", "z00000_y00000_x00128" };
    for (int c = 0; c < 2; c++) {
        float v[4 * 3] = { 0, 0, 0, 0, 0, 10, 0, 10, 0, 0, 10, 10 };
        for (int i = 0; i < 4; i++) v[i * 3 + 2] += (float)(128 * c);
        int32_t f[2 * 3] = { 0, 1, 2, 1, 3, 2 };
        snprintf(path, sizeof(path), "%s/%s_mesh.obj", dir, cid[c]);
        ps_check(ObjIO_write(path, v, 4, f, 2) == 0, "write mesh", &fails);
        float uvphi[4 * 3];
        for (int i = 0; i < 4; i++) {
            uvphi[i * 3 + 0] = (float)(100 * c) + v[i * 3 + 2];  /* u */
            uvphi[i * 3 + 1] = v[i * 3 + 1];                     /* v */
            uvphi[i * 3 + 2] = (float)(7 * c + i);               /* phi */
        }
        snprintf(path, sizeof(path), "%s/%s_uvphi.f32", dir, cid[c]);
        ps_check(ps_write_blob(path, uvphi, sizeof(uvphi)) == 0, "write uvphi",
                 &fails);
        uint8_t keep[2] = { 1, (uint8_t)(c == 0 ? 1 : 0) };
        snprintf(path, sizeof(path), "%s/%s_facekeep.u8", dir, cid[c]);
        ps_check(ps_write_blob(path, keep, 2) == 0, "write facekeep", &fails);
        if (c == 0) {   /* gid sidecar only for cube A: B exercises -1 fill */
            int32_t g[4] = { 0, 0, 1, -1 };
            snprintf(path, sizeof(path), "%s/%s_group.i32", dir, cid[c]);
            ps_check(ps_write_blob(path, g, sizeof(g)) == 0, "write group",
                     &fails);
        } else {
            snprintf(path, sizeof(path), "%s/%s_group.i32", dir, cid[c]);
            ves_unlink(path);
        }
    }

    PieceSet ps;
    int rc = PieceSet_build(arena, dir, &ps);
    fprintf(stderr, "[piece_set selftest] build rc=%d cubes=%zu nv=%zu nf=%zu "
            "u=[%.0f,%.0f] v=[%.0f,%.0f]\n", rc, ps.n_cubes, ps.nv, ps.nf,
            ps.u_min, ps.u_max, ps.v_min, ps.v_max);
    ps_check(rc == 0, "build rc", &fails);
    ps_check(ps.n_cubes == 2 && ps.nv == 8 && ps.nf == 3,
             "counts (2 cubes, 8 verts, 3 kept faces)", &fails);
    if (ps.nf == 3) {
        ps_check(ps.face_cube[0] == 0 && ps.face_cube[1] == 0 &&
                 ps.face_cube[2] == 1, "face provenance", &fails);
        int in_range = 1;
        for (size_t t = 0; t < ps.nf * 3; t++)
            if (ps.faces[t] < 0 || (size_t)ps.faces[t] >= ps.nv) in_range = 0;
        ps_check(in_range, "renumbered faces in range", &fails);
        /* cube B's kept face references verts 4..7 */
        ps_check(ps.faces[2 * 3 + 0] >= 4, "cube B faces offset", &fails);
    }
    ps_check(ps.u_min >= -1e-3 && ps.u_max <= 238.001, "u extent", &fails);
    /* gid roundtrip: cube A carries its sidecar, cube B fills -1 */
    if (ps.nv == 8) {
        ps_check(ps.gid[0] == 0 && ps.gid[2] == 1 && ps.gid[3] == -1,
                 "cube A gid roundtrip", &fails);
        ps_check(ps.gid[4] == -1 && ps.gid[7] == -1,
                 "cube B gid -1 fill", &fails);
        ps_check(ps.phi[2] == 2.0f && ps.phi[5] == 8.0f,
                 "phi channel roundtrip", &fails);
    }
    /* cube ranges + origins */
    if (ps.n_cubes == 2) {
        ps_check(ps.cube_voff[0] == 0 && ps.cube_voff[1] == 4 &&
                 ps.cube_voff[2] == 8, "cube_voff ranges", &fails);
        ps_check(ps.cube_org[0][0] == 0 && ps.cube_org[0][2] == 0 &&
                 ps.cube_org[1][2] == 128, "cube_org parsed from id", &fails);
    }
    /* facekeep compaction keeps provenance aligned */
    if (ps.nf == 3) {
        uint8_t km[3] = { 1, 0, 1 };
        size_t nn = PieceSet_apply_facekeep(&ps, km);
        ps_check(nn == 2 && ps.nf == 2, "apply_facekeep count", &fails);
        ps_check(ps.face_cube[0] == 0 && ps.face_cube[1] == 1,
                 "apply_facekeep provenance", &fails);
        ps_check(ps.faces[1 * 3 + 0] >= 4, "apply_facekeep face indices",
                 &fails);
    }
    /* normals: every vertex referenced by a kept face has a unit normal */
    {
        int nrm_ok = 1;
        for (size_t t = 0; t < ps.nf; t++) {
            for (int e = 0; e < 3; e++) {
                size_t i = (size_t)ps.faces[t * 3 + e];
                double n2 = 0.0;
                for (int k = 0; k < 3; k++)
                    n2 += (double)ps.normals[i * 3 + k] *
                          (double)ps.normals[i * 3 + k];
                if (fabs(n2 - 1.0) > 1e-3) nrm_ok = 0;
            }
        }
        ps_check(nrm_ok, "kept-face verts carry unit normals", &fails);
    }

    Arena_dispose(&arena);
    fprintf(stderr, "[piece_set selftest] %s (%d failure%s)\n",
            fails == 0 ? "PASSED" : "FAILED", fails, fails == 1 ? "" : "s");
    return fails;
}
