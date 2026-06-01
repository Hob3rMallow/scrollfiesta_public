/*
 * halo_loader.c -- load a padded volume from a cube + its 26 neighbors.
 * See halo_loader.h for interface contract.
 */
#include "ves_platform.h"

#include "halo_loader.h"
#include "tiff_io.h"
#include "arena.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int parse_cube_id(const char *cube_id, int *vz, int *vy, int *vx)
{
    if (sscanf(cube_id, "z%d_y%d_x%d", vz, vy, vx) != 3) return -1;
    return 0;
}

static int file_exists(const char *path)
{
    FILE *fp = fopen(path, "rb");
    if (!fp) return 0;
    fclose(fp);
    return 1;
}

int HaloLoader_load(Arena_T arena,
                    const char *pred_dir,
                    const char *cube_id,
                    int cube_size,
                    int halo_voxels,
                    uint8_t **out_vol,
                    int *out_p_size,
                    int64_t out_origin[3])
{
    assert(arena);
    assert(pred_dir);
    assert(cube_id);
    assert(cube_size > 0);
    assert(halo_voxels >= 0);
    assert(out_vol && out_p_size && out_origin);

    int vz = 0, vy = 0, vx = 0;
    if (parse_cube_id(cube_id, &vz, &vy, &vx) != 0) {
        fprintf(stderr, "HaloLoader: cannot parse cube id '%s'\n", cube_id);
        return -1;
    }
    out_origin[0] = (int64_t)vz - (int64_t)halo_voxels;
    out_origin[1] = (int64_t)vy - (int64_t)halo_voxels;
    out_origin[2] = (int64_t)vx - (int64_t)halo_voxels;

    int p = cube_size + 2 * halo_voxels;
    *out_p_size = p;

    size_t pvol = (size_t)p * (size_t)p * (size_t)p;
    uint8_t *padded = (uint8_t *)ARENA_CALLOC(arena, (long)pvol, 1L);
    *out_vol = padded;

    /* Scratch arena for per-neighbor TIFF loads; reset between neighbors. */
    Arena_T scratch = Arena_new();

    int loaded = 0, missing = 0;
    int self_ok = 0;

    for (int dz = -1; dz <= 1; dz++)
    for (int dy = -1; dy <= 1; dy++)
    for (int dx = -1; dx <= 1; dx++)
    {
        int nvz = vz + dz * cube_size;
        int nvy = vy + dy * cube_size;
        int nvx = vx + dx * cube_size;

        /* For halo == 0, only the self cube has any nonempty padded region;
         * skip neighbors to avoid wasted I/O. */
        if (halo_voxels == 0 && (dz != 0 || dy != 0 || dx != 0)) continue;

        if (nvz < 0 || nvy < 0 || nvx < 0) {
            if (dz != 0 || dy != 0 || dx != 0) missing++;
            continue;
        }

        char nid[128];
        snprintf(nid, sizeof(nid), "z%05d_y%05d_x%05d", nvz, nvy, nvx);
        char npath[1024];
        snprintf(npath, sizeof(npath), "%s/%s.tif", pred_dir, nid);

        if (!file_exists(npath)) {
            if (dz == 0 && dy == 0 && dx == 0) {
                fprintf(stderr, "HaloLoader: self TIFF missing: %s\n", npath);
                Arena_dispose(&scratch);
                return -1;
            }
            missing++;
            continue;
        }

        Arena_free(scratch);
        uint8_t *nvol = NULL;
        int nD = 0, nH = 0, nW = 0;
        if (TiffIO_load(scratch, npath, &nvol, &nD, &nH, &nW) != 0) {
            if (dz == 0 && dy == 0 && dx == 0) {
                fprintf(stderr, "HaloLoader: self TIFF load failed: %s\n", npath);
                Arena_dispose(&scratch);
                return -1;
            }
            fprintf(stderr, "HaloLoader: neighbor TIFF load failed: %s\n", npath);
            missing++;
            continue;
        }
        if (nD != cube_size || nH != cube_size || nW != cube_size) {
            fprintf(stderr,
                    "HaloLoader: %s shape %dx%dx%d != expected %d^3\n",
                    npath, nD, nH, nW, cube_size);
            if (dz == 0 && dy == 0 && dx == 0) {
                Arena_dispose(&scratch);
                return -1;
            }
            missing++;
            continue;
        }

        /* Determine padded subregion and neighbor subregion for this offset.
         * Each axis is one of:
         *   d == -1 : padded [0, halo), neighbor [cube_size - halo, cube_size)
         *   d ==  0 : padded [halo, halo + cube_size), neighbor [0, cube_size)
         *   d == +1 : padded [halo + cube_size, p), neighbor [0, halo)
         */
        int pz_lo, pz_hi, nz_lo;
        if      (dz == -1) { pz_lo = 0;                       pz_hi = halo_voxels;              nz_lo = cube_size - halo_voxels; }
        else if (dz ==  0) { pz_lo = halo_voxels;             pz_hi = halo_voxels + cube_size;  nz_lo = 0; }
        else               { pz_lo = halo_voxels + cube_size; pz_hi = p;                        nz_lo = 0; }

        int py_lo, py_hi, ny_lo;
        if      (dy == -1) { py_lo = 0;                       py_hi = halo_voxels;              ny_lo = cube_size - halo_voxels; }
        else if (dy ==  0) { py_lo = halo_voxels;             py_hi = halo_voxels + cube_size;  ny_lo = 0; }
        else               { py_lo = halo_voxels + cube_size; py_hi = p;                        ny_lo = 0; }

        int px_lo, px_hi, nx_lo;
        if      (dx == -1) { px_lo = 0;                       px_hi = halo_voxels;              nx_lo = cube_size - halo_voxels; }
        else if (dx ==  0) { px_lo = halo_voxels;             px_hi = halo_voxels + cube_size;  nx_lo = 0; }
        else               { px_lo = halo_voxels + cube_size; px_hi = p;                        nx_lo = 0; }

        int dz_extent = pz_hi - pz_lo;
        int dy_extent = py_hi - py_lo;
        int dx_extent = px_hi - px_lo;
        if (dz_extent <= 0 || dy_extent <= 0 || dx_extent <= 0) {
            /* Halo of zero on at least one axis for an off-center neighbor -- skip. */
            continue;
        }

        for (int k = 0; k < dz_extent; k++) {
            int pz = pz_lo + k;
            int nz = nz_lo + k;
            for (int j = 0; j < dy_extent; j++) {
                int py = py_lo + j;
                int ny = ny_lo + j;
                uint8_t *dst = padded
                    + (size_t)pz * (size_t)p * (size_t)p
                    + (size_t)py * (size_t)p
                    + (size_t)px_lo;
                const uint8_t *src = nvol
                    + (size_t)nz * (size_t)cube_size * (size_t)cube_size
                    + (size_t)ny * (size_t)cube_size
                    + (size_t)nx_lo;
                memcpy(dst, src, (size_t)dx_extent);
            }
        }

        if (dz == 0 && dy == 0 && dx == 0) self_ok = 1;
        loaded++;
    }

    Arena_dispose(&scratch);

    if (!self_ok) {
        fprintf(stderr, "HaloLoader: self cube did not load (cube_id=%s)\n",
                cube_id);
        return -1;
    }

    fprintf(stderr,
            "HaloLoader: cube %s p_size=%d halo=%d  loaded=%d/27 missing=%d\n",
            cube_id, p, halo_voxels, loaded, missing);
    return 0;
}
