/*
 * dump_obj.c — OBJ dump helpers for intermediate pipeline output.
 *
 * Layout: dump_dir/cube_id/cube_id_stepX/cube_id_stepX_all.obj
 *                                        cube_id_stepX_001.obj
 *         dump_dir/cube_id/cube_id_all_obj/cube_id_stepX_all.obj  (copies)
 */
#include "ves_platform.h"
#include "dump_obj.h"
#include "obj_io.h"
#include "obj_colors.h"

#include <stdio.h>
#include <string.h>
#include <errno.h>

#define COMP_PATH_LEN 1024

/* ----------------------------------------------------------------
 * Utilities (file-private)
 * ---------------------------------------------------------------- */

/* Color index for a component — use comp_id so colors are stable across steps. */
static int comp_color_idx(const ComponentMesh *cm)
{
    return (cm->comp_id > 0 ? cm->comp_id - 1 : 0) % OBJ_NUM_COLORS;
}

/* Parse cube origin from cube_id format "z{vz}_y{vy}_x{vx}".
 * Returns 0 on success (origin populated in (z, y, x) order to match the
 * pipeline's vertex order); -1 on parse failure (origin left at 0,0,0). */
static int parse_cube_origin_zyx(const char *cube_id, float origin[3])
{
    origin[0] = 0.0f;
    origin[1] = 0.0f;
    origin[2] = 0.0f;
    if (!cube_id) return -1;
    int iz = 0, iy = 0, ix = 0;
    if (sscanf(cube_id, "z%d_y%d_x%d", &iz, &iy, &ix) != 3) return -1;
    origin[0] = (float)iz;
    origin[1] = (float)iy;
    origin[2] = (float)ix;
    return 0;
}

/* Add origin to every vertex in src, write result to dst. Both arrays
 * are flat [nv*3] in (z, y, x) order matching the pipeline convention.
 * dst and src may not alias. */
static void offset_verts(const float *src, float *dst, size_t nv,
                         const float origin[3])
{
    for (size_t i = 0; i < nv; i++) {
        dst[i * 3 + 0] = src[i * 3 + 0] + origin[0];
        dst[i * 3 + 1] = src[i * 3 + 1] + origin[1];
        dst[i * 3 + 2] = src[i * 3 + 2] + origin[2];
    }
}

/* Copy a file byte-for-byte.  Returns 0 on success. */
static int copy_file(const char *src, const char *dst)
{
    FILE *in = fopen(src, "rb");
    if (!in) return -1;
    FILE *out = fopen(dst, "wb");
    if (!out) { fclose(in); return -1; }
    char buf[8192];
    size_t n = 0;
    while ((n = fread(buf, 1, sizeof(buf), in)) > 0)
        fwrite(buf, 1, n, out);
    fclose(in);
    fclose(out);
    return 0;
}

/* Copy the combined all.obj into cube_dir/cube_id_all_obj/.
 * `dir` is "cube_dir/cube_id_stepX"; parent is derived by strrchr. */
static void copy_to_all_obj_dir(const char *dir, const char *cube_id,
                                 const char *all_path)
{
    /* Find parent of dir (cube_dir) — handle both / and \ separators */
    char parent[1024];
    snprintf(parent, sizeof(parent), "%s", dir);
    char *slash_fwd = strrchr(parent, '/');
    char *slash_bck = strrchr(parent, '\\');
    char *slash = slash_fwd;
    if (slash_bck && (!slash || slash_bck > slash)) slash = slash_bck;
    if (!slash) return;
    *slash = '\0';

    /* Create cube_id_all_obj/ under parent */
    char all_dir[1024];
    snprintf(all_dir, sizeof(all_dir), "%.512s/%.256s_all_obj", parent, cube_id);
    if (ves_mkdir(all_dir) != 0 && errno != EEXIST) return;

    /* Derive destination filename from all_path basename */
    const char *base = ves_path_basename(all_path);
    char dst[1024];
    snprintf(dst, sizeof(dst), "%.768s/%.200s", all_dir, base);
    copy_file(all_path, dst);
}

/* Combine individual OBJ files into one, adjusting face indices.
 * Reads each file line by line; "v" lines pass through, "f" lines
 * get vertex offset added.  An "o comp_NNN" line is inserted before
 * each component. */
static void combine_objs(const char *const *paths, const size_t *nv_arr,
                          size_t count, const char *out_path)
{
    FILE *out = fopen(out_path, "w");
    if (!out) return;

    size_t vert_offset = 0;
    for (size_t c = 0; c < count; c++) {
        fprintf(out, "o comp_%03d\n", (int)(c + 1));
        FILE *in = fopen(paths[c], "r");
        if (!in) { vert_offset += nv_arr[c]; continue; }

        char line[256];
        while (fgets(line, (int)sizeof(line), in)) {
            if (line[0] == 'v' && line[1] == ' ') {
                fputs(line, out);
            } else if (line[0] == 'f' && line[1] == ' ') {
                int a = 0, b = 0, c2 = 0;
                if (sscanf(line + 2, "%d %d %d", &a, &b, &c2) == 3)
                    fprintf(out, "f %zu %zu %zu\n",
                            (size_t)a + vert_offset,
                            (size_t)b + vert_offset,
                            (size_t)c2 + vert_offset);
            }
        }
        fclose(in);
        vert_offset += nv_arr[c];
    }
    fclose(out);
}

/* ----------------------------------------------------------------
 * Public API
 * ---------------------------------------------------------------- */

void DumpObj_extract_cube_id(const char *input_path, char *out, size_t out_sz)
{
    const char *base = ves_path_basename(input_path);
    size_t len = strlen(base);
    /* Strip .tif / .tiff extension */
    const char *dot = strrchr(base, '.');
    if (dot && dot > base) len = (size_t)(dot - base);
    if (len >= out_sz) len = out_sz - 1;
    memcpy(out, base, len);
    out[len] = '\0';
}

int DumpObj_ensure_dir(const char *path)
{
    if (ves_mkdir(path) != 0 && errno != EEXIST) {
        fprintf(stderr, "WARNING: mkdir(%s): %s\n", path, strerror(errno));
        return -1;
    }
    return 0;
}

void DumpObj_write_meshes(Arena_T arena, const char *dir,
                          const char *cube_id, const char *step_tag,
                          const ComponentMesh *meshes, size_t count)
{
    if (count == 0) return;

    Arena_Mark mark = Arena_save(arena);

    char *path_buf = (char *)ARENA_ALLOC(arena,
        (long)(count * COMP_PATH_LEN));
    size_t *nv_arr = (size_t *)ARENA_ALLOC(arena,
        (long)(count * sizeof(size_t)));
    const char **path_ptrs = (const char **)ARENA_ALLOC(arena,
        (long)(count * sizeof(const char *)));

    /* Write verts in WORLD source-voxel coords (cube-local + cube origin
     * parsed from cube_id). This makes per-cube dumps from different cubes
     * line up the same way grid_weld would stitch them, so loading
     * multiple cube_id_stepX_all.obj files in MeshLab gives the same
     * spatial layout as the final welded mesh. Without the offset, every
     * cube's dump lands at the origin and they all overlap. */
    float origin[3];
    parse_cube_origin_zyx(cube_id, origin);

    /* Bright blue for hole-fill verts (red carries an error/"bad" connotation).
     * Green is reserved by grid_weld for cross-cube weld highlights in the final
     * welded.obj — keep the three concerns visually distinct (per-cube fills =
     * blue, cross-cube weld points = green). */
    static const float FILL_BLUE[3] = {0.10f, 0.45f, 1.00f};

    for (size_t i = 0; i < count; i++) {
        char *comp_path = path_buf + (ptrdiff_t)(i * COMP_PATH_LEN);
        snprintf(comp_path, COMP_PATH_LEN,
                 "%.512s/%.256s_%.128s_%03d.obj",
                 dir, cube_id, step_tag, (int)(i + 1));

        float *shifted = (float *)ARENA_ALLOC(arena,
            (long)(meshes[i].nv * 3 * sizeof(float)));
        offset_verts(meshes[i].verts, shifted, meshes[i].nv, origin);

        if (meshes[i].nv_pre_fill > 0) {
            ObjIO_write_twotone(comp_path, shifted, meshes[i].nv,
                                meshes[i].faces, meshes[i].nf,
                                meshes[i].nv_pre_fill,
                                OBJ_COLORS[comp_color_idx(&meshes[i])],
                                FILL_BLUE);
        } else {
            ObjIO_write_colored(comp_path, shifted, meshes[i].nv,
                                meshes[i].faces, meshes[i].nf,
                                OBJ_COLORS[comp_color_idx(&meshes[i])]);
        }
        nv_arr[i] = meshes[i].nv;
        path_ptrs[i] = comp_path;
    }

    char all_path[1024];
    snprintf(all_path, sizeof(all_path), "%.512s/%.256s_%.128s_all.obj",
             dir, cube_id, step_tag);
    combine_objs(path_ptrs, nv_arr, count, all_path);
    copy_to_all_obj_dir(dir, cube_id, all_path);
    fprintf(stderr, "  [dump] %zu OBJs + combined -> %s\n", count, dir);

    Arena_restore(arena, mark);
}

void DumpObj_write_one_world(Arena_T arena, const char *path,
                             const char *cube_id,
                             const float *verts, size_t nv,
                             const int32_t *faces, size_t nf,
                             const float color[3])
{
    if (nv == 0 || !path) return;

    float origin[3];
    parse_cube_origin_zyx(cube_id, origin);

    Arena_Mark mark = Arena_save(arena);
    float *shifted = (float *)ARENA_ALLOC(arena,
        (long)(nv * 3 * sizeof(float)));
    offset_verts(verts, shifted, nv, origin);
    ObjIO_write_colored(path, shifted, nv, faces, nf, color);
    Arena_restore(arena, mark);
}

void DumpObj_write_points_world(Arena_T arena, const char *path,
                                const char *cube_id,
                                const float *verts, const float *normals,
                                size_t nv)
{
    if (nv == 0 || !path) return;

    float origin[3];
    parse_cube_origin_zyx(cube_id, origin);

    Arena_Mark mark = Arena_save(arena);
    float *shifted = (float *)ARENA_ALLOC(arena,
        (long)(nv * 3 * sizeof(float)));
    offset_verts(verts, shifted, nv, origin);

    FILE *fp = fopen(path, "w");
    if (!fp) { Arena_restore(arena, mark); return; }
    /* Verts are stored in (z, y, x). Write them as-is to match the OBJ
     * convention used everywhere else in this codebase (ObjIO_write_colored
     * does the same — MeshLab interprets the first number as X, but
     * downstream tools in this project all assume this layout). */
    for (size_t i = 0; i < nv; i++) {
        fprintf(fp, "v %.6f %.6f %.6f\n",
                (double)shifted[i * 3 + 0],
                (double)shifted[i * 3 + 1],
                (double)shifted[i * 3 + 2]);
    }
    if (normals) {
        for (size_t i = 0; i < nv; i++) {
            fprintf(fp, "vn %.6f %.6f %.6f\n",
                    (double)normals[i * 3 + 0],
                    (double)normals[i * 3 + 1],
                    (double)normals[i * 3 + 2]);
        }
    }
    fclose(fp);
    Arena_restore(arena, mark);
}

void DumpObj_write_mesh_ptrs(Arena_T arena, const char *dir,
                             const char *cube_id, const char *step_tag,
                             const ComponentMesh **meshes, size_t count)
{
    if (count == 0) return;

    Arena_Mark mark = Arena_save(arena);

    char *path_buf = (char *)ARENA_ALLOC(arena,
        (long)(count * COMP_PATH_LEN));
    size_t *nv_arr = (size_t *)ARENA_ALLOC(arena,
        (long)(count * sizeof(size_t)));
    const char **path_ptrs = (const char **)ARENA_ALLOC(arena,
        (long)(count * sizeof(const char *)));

    /* World-coord offset (see DumpObj_write_meshes for rationale). */
    float origin[3];
    parse_cube_origin_zyx(cube_id, origin);

    for (size_t i = 0; i < count; i++) {
        char *comp_path = path_buf + (ptrdiff_t)(i * COMP_PATH_LEN);
        snprintf(comp_path, COMP_PATH_LEN,
                 "%.512s/%.256s_%.128s_%03d.obj",
                 dir, cube_id, step_tag, (int)(i + 1));
        float *shifted = (float *)ARENA_ALLOC(arena,
            (long)(meshes[i]->nv * 3 * sizeof(float)));
        offset_verts(meshes[i]->verts, shifted, meshes[i]->nv, origin);
        ObjIO_write_colored(comp_path, shifted, meshes[i]->nv,
                            meshes[i]->faces, meshes[i]->nf,
                            OBJ_COLORS[comp_color_idx(meshes[i])]);
        nv_arr[i] = meshes[i]->nv;
        path_ptrs[i] = comp_path;
    }

    char all_path[1024];
    snprintf(all_path, sizeof(all_path), "%.512s/%.256s_%.128s_all.obj",
             dir, cube_id, step_tag);
    combine_objs(path_ptrs, nv_arr, count, all_path);
    copy_to_all_obj_dir(dir, cube_id, all_path);
    fprintf(stderr, "  [dump] %zu OBJs + combined -> %s\n", count, dir);

    Arena_restore(arena, mark);
}
