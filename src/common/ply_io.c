#include "ply_io.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int PlyIO_write_points(const char *path,
                       const float *verts, const float *normals, size_t n)
{
    assert(path);

    if (n == 0) {
        return -1;
    }
    assert(verts && normals);

    FILE *f = fopen(path, "wb");
    if (f == NULL) {
        return -1;
    }

    /* Write header */
    fprintf(f, "ply\n");
    fprintf(f, "format binary_little_endian 1.0\n");
    fprintf(f, "element vertex %zu\n", n);
    fprintf(f, "property float x\n");
    fprintf(f, "property float y\n");
    fprintf(f, "property float z\n");
    fprintf(f, "property float nx\n");
    fprintf(f, "property float ny\n");
    fprintf(f, "property float nz\n");
    fprintf(f, "end_header\n");

    /* Buffer and write binary data in one write when possible.
     * Each vertex: 6 floats = 24 bytes. */
    for (size_t i = 0; i < n; i++) {
        float data[6];
        data[0] = verts[i * 3 + 0];
        data[1] = verts[i * 3 + 1];
        data[2] = verts[i * 3 + 2];
        data[3] = normals[i * 3 + 0];
        data[4] = normals[i * 3 + 1];
        data[5] = normals[i * 3 + 2];
        if (fwrite(data, sizeof(float), 6, f) != 6) {
            fclose(f);
            return -1;
        }
    }

    fclose(f);
    return 0;
}

/* Helper: read a line from file, return length or -1 on EOF */
static int read_line(FILE *f, char *buf, int bufsize)
{
    if (fgets(buf, bufsize, f) == NULL) {
        return -1;
    }
    return (int)strlen(buf);
}

int PlyIO_read_mesh(Arena_T arena, const char *path,
                    float **out_verts, size_t *out_nv,
                    int32_t **out_faces, size_t *out_nf,
                    float **out_density)
{
    assert(arena);
    assert(path);
    assert(out_verts && out_nv && out_faces && out_nf);

    FILE *f = fopen(path, "rb");
    if (f == NULL) {
        return -1;
    }

    char line[512];
    size_t n_verts = 0;
    size_t n_faces = 0;
    int has_density = 0;
    int n_vert_props = 0;   /* count of float properties per vertex */
    int in_vertex_element = 0;

    /* Parse header */
    while (read_line(f, line, (int)sizeof(line)) >= 0) {
        if (strncmp(line, "end_header", 10) == 0) {
            break;
        }
        if (strncmp(line, "element vertex", 14) == 0) {
            sscanf(line, "element vertex %zu", &n_verts);
            in_vertex_element = 1;
            n_vert_props = 0;
        } else if (strncmp(line, "element face", 12) == 0) {
            sscanf(line, "element face %zu", &n_faces);
            in_vertex_element = 0;
        } else if (strncmp(line, "element", 7) == 0) {
            in_vertex_element = 0;
        } else if (strncmp(line, "property float", 14) == 0 && in_vertex_element) {
            n_vert_props++;
            if (strstr(line, "quality") || strstr(line, "value") ||
                strstr(line, "density")) {
                has_density = 1;
            }
        }
    }

    if (n_verts == 0) {
        fclose(f);
        *out_nv = 0;
        *out_nf = 0;
        *out_verts = NULL;
        *out_faces = NULL;
        if (out_density) { *out_density = NULL; }
        return 0;
    }

    /* Allocate output arrays */
    float *verts = (float *)ARENA_ALLOC(arena,
                                         (long)(n_verts * 3 * sizeof(float)));
    float *density = NULL;
    if (has_density && out_density) {
        density = (float *)ARENA_ALLOC(arena,
                                        (long)(n_verts * sizeof(float)));
    }

    /* Read vertex data */
    for (size_t i = 0; i < n_verts; i++) {
        float props[16]; /* up to 16 properties */
        int to_read = (n_vert_props < 16) ? n_vert_props : 16;
        if (fread(props, sizeof(float), (size_t)to_read, f) != (size_t)to_read) {
            fclose(f);
            return -1;
        }
        verts[i * 3 + 0] = props[0];
        verts[i * 3 + 1] = props[1];
        verts[i * 3 + 2] = props[2];

        /* Density is typically the 7th property (after x,y,z,nx,ny,nz) */
        if (density && n_vert_props >= 7) {
            density[i] = props[6];
        } else if (density && n_vert_props >= 4) {
            density[i] = props[3];
        }
    }

    /* Detect face list count type from header.
     * PoissonRecon uses "property list int int vertex_indices" (4-byte count).
     * Other tools may use "property list uchar int vertex_indices" (1-byte). */
    int face_count_bytes = 1; /* default: uchar */
    /* Re-read header to check face list format */
    {
        long saved_pos = ftell(f);
        rewind(f);
        char hdr_line[512];
        while (read_line(f, hdr_line, (int)sizeof(hdr_line)) >= 0) {
            if (strncmp(hdr_line, "end_header", 10) == 0) break;
            if (strncmp(hdr_line, "property list", 13) == 0) {
                if (strstr(hdr_line, "list int ") != NULL ||
                    strstr(hdr_line, "list int32 ") != NULL) {
                    face_count_bytes = 4;
                }
            }
        }
        fseek(f, saved_pos, SEEK_SET);
    }

    /* Read face data */
    int32_t *faces = NULL;
    if (n_faces > 0) {
        faces = (int32_t *)ARENA_ALLOC(arena,
                                        (long)(n_faces * 3 * sizeof(int32_t)));
        for (size_t i = 0; i < n_faces; i++) {
            /* Read face vertex count (1 or 4 bytes depending on PLY format) */
            int32_t count = 0;
            if (face_count_bytes == 4) {
                if (fread(&count, sizeof(int32_t), 1, f) != 1) {
                    fclose(f);
                    return -1;
                }
            } else {
                uint8_t count8 = 0;
                if (fread(&count8, 1, 1, f) != 1) {
                    fclose(f);
                    return -1;
                }
                count = (int32_t)count8;
            }

            if (count == 3) {
                if (fread(&faces[i * 3], sizeof(int32_t), 3, f) != 3) {
                    fclose(f);
                    return -1;
                }
            } else {
                /* Non-triangle: read and discard, fill degenerate */
                for (int32_t j = 0; j < count; j++) {
                    int32_t dummy = 0;
                    if (fread(&dummy, sizeof(int32_t), 1, f) != 1) {
                        fclose(f);
                        return -1;
                    }
                }
                faces[i * 3 + 0] = 0;
                faces[i * 3 + 1] = 0;
                faces[i * 3 + 2] = 0;
            }
        }
    }

    fclose(f);

    *out_verts = verts;
    *out_nv = n_verts;
    *out_faces = faces;
    *out_nf = n_faces;
    if (out_density) {
        *out_density = density;
    }

    return 0;
}

int PlyIO_write_winding_mesh(const char *path,
                              const float *verts, size_t nv,
                              const int32_t *faces, size_t nf,
                              const float *winding,
                              const float *confidence,
                              const uint8_t *defect)
{
    assert(path && verts && faces && winding && confidence && defect);
    if (nv == 0 || nf == 0) return -1;

    FILE *f = fopen(path, "wb");
    if (!f) return -1;

    fprintf(f, "ply\n");
    fprintf(f, "format binary_little_endian 1.0\n");
    fprintf(f, "comment vesuvius-c winding diagnostic\n");
    fprintf(f, "element vertex %zu\n", nv);
    fprintf(f, "property float x\n");
    fprintf(f, "property float y\n");
    fprintf(f, "property float z\n");
    fprintf(f, "property float winding\n");
    fprintf(f, "property float confidence\n");
    fprintf(f, "property uchar defect\n");
    fprintf(f, "element face %zu\n", nf);
    fprintf(f, "property list uchar int vertex_indices\n");
    fprintf(f, "end_header\n");

    /* Vertices: project verts are (z, y, x); emit as (x, y, z) for
     * MeshLab-correct axis orientation. */
    for (size_t i = 0; i < nv; i++) {
        float vz = verts[i * 3 + 0];
        float vy = verts[i * 3 + 1];
        float vx = verts[i * 3 + 2];
        float w  = winding[i];
        float c  = confidence[i];
        uint8_t d = defect[i];
        if (fwrite(&vx, sizeof(float), 1, f) != 1 ||
            fwrite(&vy, sizeof(float), 1, f) != 1 ||
            fwrite(&vz, sizeof(float), 1, f) != 1 ||
            fwrite(&w,  sizeof(float), 1, f) != 1 ||
            fwrite(&c,  sizeof(float), 1, f) != 1 ||
            fwrite(&d,  sizeof(uint8_t), 1, f) != 1) {
            fclose(f);
            return -1;
        }
    }

    /* Faces: PLY list of 3 int32 vertex indices, with leading uchar
     * count = 3. */
    for (size_t i = 0; i < nf; i++) {
        uint8_t count = 3;
        if (fwrite(&count, sizeof(uint8_t), 1, f) != 1 ||
            fwrite(faces + i * 3, sizeof(int32_t), 3, f) != 3) {
            fclose(f);
            return -1;
        }
    }

    fclose(f);
    return 0;
}
