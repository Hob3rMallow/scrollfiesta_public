#include "obj_io.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int ObjIO_write(const char *path, const float *verts, size_t nv,
                const int32_t *faces, size_t nf)
{
    assert(path);

    FILE *f = fopen(path, "w");
    if (f == NULL) {
        return -1;
    }

    /* Write vertices */
    if (verts) {
        for (size_t i = 0; i < nv; i++) {
            fprintf(f, "v %.6f %.6f %.6f\n",
                    (double)verts[i * 3 + 0],
                    (double)verts[i * 3 + 1],
                    (double)verts[i * 3 + 2]);
        }
    }

    /* Write faces (OBJ is 1-indexed) */
    if (faces) {
        for (size_t i = 0; i < nf; i++) {
            fprintf(f, "f %d %d %d\n",
                    faces[i * 3 + 0] + 1,
                    faces[i * 3 + 1] + 1,
                    faces[i * 3 + 2] + 1);
        }
    }

    fclose(f);
    return 0;
}

int ObjIO_write_colored(const char *path, const float *verts, size_t nv,
                        const int32_t *faces, size_t nf,
                        const float color[3])
{
    assert(path);
    assert(color);

    FILE *f = fopen(path, "w");
    if (f == NULL) {
        return -1;
    }

    /* Write vertices with RGB color */
    if (verts) {
        for (size_t i = 0; i < nv; i++) {
            fprintf(f, "v %.6f %.6f %.6f %.4f %.4f %.4f\n",
                    (double)verts[i * 3 + 0],
                    (double)verts[i * 3 + 1],
                    (double)verts[i * 3 + 2],
                    (double)color[0],
                    (double)color[1],
                    (double)color[2]);
        }
    }

    /* Write faces (OBJ is 1-indexed) */
    if (faces) {
        for (size_t i = 0; i < nf; i++) {
            fprintf(f, "f %d %d %d\n",
                    faces[i * 3 + 0] + 1,
                    faces[i * 3 + 1] + 1,
                    faces[i * 3 + 2] + 1);
        }
    }

    fclose(f);
    return 0;
}

int ObjIO_write_twotone(const char *path, const float *verts, size_t nv,
                        const int32_t *faces, size_t nf,
                        size_t nv_split,
                        const float color1[3], const float color2[3])
{
    assert(path);
    assert(color1);
    assert(color2);

    FILE *f = fopen(path, "w");
    if (f == NULL) {
        return -1;
    }

    if (verts) {
        for (size_t i = 0; i < nv; i++) {
            const float *c = (i < nv_split) ? color1 : color2;
            fprintf(f, "v %.6f %.6f %.6f %.4f %.4f %.4f\n",
                    (double)verts[i * 3 + 0],
                    (double)verts[i * 3 + 1],
                    (double)verts[i * 3 + 2],
                    (double)c[0], (double)c[1], (double)c[2]);
        }
    }

    if (faces) {
        for (size_t i = 0; i < nf; i++) {
            fprintf(f, "f %d %d %d\n",
                    faces[i * 3 + 0] + 1,
                    faces[i * 3 + 1] + 1,
                    faces[i * 3 + 2] + 1);
        }
    }

    fclose(f);
    return 0;
}

int ObjIO_write_per_vertex_color(const char *path,
                                 const float *verts, size_t nv,
                                 const int32_t *faces, size_t nf,
                                 const float *colors)
{
    assert(path);
    assert(colors);

    FILE *f = fopen(path, "w");
    if (f == NULL) return -1;

    if (verts) {
        for (size_t i = 0; i < nv; i++) {
            fprintf(f, "v %.6f %.6f %.6f %.4f %.4f %.4f\n",
                    (double)verts[i * 3 + 0],
                    (double)verts[i * 3 + 1],
                    (double)verts[i * 3 + 2],
                    (double)colors[i * 3 + 0],
                    (double)colors[i * 3 + 1],
                    (double)colors[i * 3 + 2]);
        }
    }

    if (faces) {
        for (size_t i = 0; i < nf; i++) {
            fprintf(f, "f %d %d %d\n",
                    faces[i * 3 + 0] + 1,
                    faces[i * 3 + 1] + 1,
                    faces[i * 3 + 2] + 1);
        }
    }

    fclose(f);
    return 0;
}

int ObjIO_read(Arena_T arena, const char *path,
               float **out_verts, size_t *out_nv,
               int32_t **out_faces, size_t *out_nf)
{
    assert(arena);
    assert(path);
    assert(out_verts && out_nv && out_faces && out_nf);

    FILE *f = fopen(path, "r");
    if (f == NULL) {
        return -1;
    }

    /* First pass: count vertices and faces */
    size_t nv = 0, nf = 0;
    char line[1024];

    while (fgets(line, (int)sizeof(line), f) != NULL) {
        if (line[0] == 'v' && line[1] == ' ') {
            nv++;
        } else if (line[0] == 'f' && line[1] == ' ') {
            nf++;
        }
    }

    /* Allocate */
    float *verts = NULL;
    int32_t *faces = NULL;

    if (nv > 0) {
        verts = (float *)ARENA_ALLOC(arena,
                                      (long)(nv * 3 * sizeof(float)));
    }
    if (nf > 0) {
        faces = (int32_t *)ARENA_ALLOC(arena,
                                        (long)(nf * 3 * sizeof(int32_t)));
    }

    /* Second pass: read data */
    rewind(f);
    size_t vi = 0, fi = 0;

    while (fgets(line, (int)sizeof(line), f) != NULL) {
        if (line[0] == 'v' && line[1] == ' ') {
            float a = 0.0f, b = 0.0f, c = 0.0f;
            sscanf(line + 2, "%f %f %f", &a, &b, &c);
            verts[vi * 3 + 0] = a;
            verts[vi * 3 + 1] = b;
            verts[vi * 3 + 2] = c;
            vi++;
        } else if (line[0] == 'f' && line[1] == ' ') {
            /* Handle "f v1 v2 v3" and "f v1/vt1 v2/vt2 v3/vt3" */
            int i0 = 0, i1 = 0, i2 = 0;
            const char *p = line + 2;

            i0 = (int)strtol(p, (char **)&p, 10);
            /* Skip any /vt/vn */
            while (*p != '\0' && *p != ' ' && *p != '\t' && *p != '\n') { p++; }
            while (*p == ' ' || *p == '\t') { p++; }

            i1 = (int)strtol(p, (char **)&p, 10);
            while (*p != '\0' && *p != ' ' && *p != '\t' && *p != '\n') { p++; }
            while (*p == ' ' || *p == '\t') { p++; }

            i2 = (int)strtol(p, (char **)&p, 10);

            /* Convert from 1-indexed to 0-indexed */
            faces[fi * 3 + 0] = i0 - 1;
            faces[fi * 3 + 1] = i1 - 1;
            faces[fi * 3 + 2] = i2 - 1;
            fi++;
        }
    }

    fclose(f);

    *out_verts = verts;
    *out_nv = nv;
    *out_faces = faces;
    *out_nf = nf;
    return 0;
}

int ObjIO_write_uv(const char *path, const float *verts, size_t nv,
                   const int32_t *faces, size_t nf, const float *uv)
{
    assert(path);
    assert(uv || nv == 0);

    FILE *f = fopen(path, "w");
    if (f == NULL) {
        return -1;
    }

    /* Positions (native (z,y,x) order, matching the rest of the pipeline). */
    if (verts) {
        for (size_t i = 0; i < nv; i++) {
            fprintf(f, "v %.6f %.6f %.6f\n",
                    (double)verts[i * 3 + 0],
                    (double)verts[i * 3 + 1],
                    (double)verts[i * 3 + 2]);
        }
    }

    /* Texture coordinates (u, v). */
    if (uv) {
        for (size_t i = 0; i < nv; i++) {
            fprintf(f, "vt %.6f %.6f\n",
                    (double)uv[i * 2 + 0],
                    (double)uv[i * 2 + 1]);
        }
    }

    /* Faces with shared vertex/uv indices (OBJ is 1-indexed). */
    if (faces) {
        for (size_t i = 0; i < nf; i++) {
            int a = faces[i * 3 + 0] + 1;
            int b = faces[i * 3 + 1] + 1;
            int c = faces[i * 3 + 2] + 1;
            fprintf(f, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
        }
    }

    fclose(f);
    return 0;
}

int ObjIO_write_uv_masked(const char *path, const float *verts, size_t nv,
                          const int32_t *faces, size_t nf, const float *uv,
                          const uint8_t *keep)
{
    assert(path);
    assert(uv || nv == 0);

    FILE *f = fopen(path, "w");
    if (f == NULL) {
        return -1;
    }

    if (verts) {
        for (size_t i = 0; i < nv; i++) {
            fprintf(f, "v %.6f %.6f %.6f\n",
                    (double)verts[i * 3 + 0],
                    (double)verts[i * 3 + 1],
                    (double)verts[i * 3 + 2]);
        }
    }
    if (uv) {
        for (size_t i = 0; i < nv; i++) {
            fprintf(f, "vt %.6f %.6f\n",
                    (double)uv[i * 2 + 0],
                    (double)uv[i * 2 + 1]);
        }
    }
    if (faces) {
        for (size_t i = 0; i < nf; i++) {
            if (keep != NULL && keep[i] == 0) continue;
            int a = faces[i * 3 + 0] + 1;
            int b = faces[i * 3 + 1] + 1;
            int c = faces[i * 3 + 2] + 1;
            fprintf(f, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
        }
    }

    fclose(f);
    return 0;
}
