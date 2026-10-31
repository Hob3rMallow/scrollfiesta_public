#include "mesh_normals.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* One face's unnormalized winding normal, added to its three corners. */
static void mn_face_add(const float *verts, const int32_t *faces, size_t face, float *normal)
{
    size_t a = (size_t)faces[face * 3];
    size_t b = (size_t)faces[face * 3 + 1];
    size_t c = (size_t)faces[face * 3 + 2];
    double e1[3], e2[3], cross[3];
    int axis = 0;
    for (axis = 0; axis < 3; axis++) {
        e1[axis] = (double)verts[b * 3 + (size_t)axis]
                 - (double)verts[a * 3 + (size_t)axis];
        e2[axis] = (double)verts[c * 3 + (size_t)axis]
                 - (double)verts[a * 3 + (size_t)axis];
    }
    cross[0] = e1[1] * e2[2] - e1[2] * e2[1];
    cross[1] = e1[2] * e2[0] - e1[0] * e2[2];
    cross[2] = e1[0] * e2[1] - e1[1] * e2[0];
    for (axis = 0; axis < 3; axis++) {
        normal[a * 3 + (size_t)axis] += (float)cross[axis];
        normal[b * 3 + (size_t)axis] += (float)cross[axis];
        normal[c * 3 + (size_t)axis] += (float)cross[axis];
    }
}

static void mn_vertex_unit(float *normal, size_t vertex)
{
    double x = (double)normal[vertex * 3];
    double y = (double)normal[vertex * 3 + 1];
    double z = (double)normal[vertex * 3 + 2];
    double length = sqrt(x * x + y * y + z * z);
    if (length > 1.0e-12) {
        normal[vertex * 3] = (float)(x / length);
        normal[vertex * 3 + 1] = (float)(y / length);
        normal[vertex * 3 + 2] = (float)(z / length);
    } else {
        normal[vertex * 3] = 0.0f;
        normal[vertex * 3 + 1] = 0.0f;
        normal[vertex * 3 + 2] = 0.0f;
    }
}

float *MeshNormals_compute(const float *verts, size_t nv,
                           const int32_t *faces, size_t nf)
{
    float *normal = (float *)calloc(nv > 0 ? nv * 3 : 1, sizeof(*normal));
    size_t face = 0, vertex = 0;
    if (normal == NULL) {
        fprintf(stderr, "mesh_normals: out of memory (%zu vertices)\n", nv);
        return NULL;
    }

    for (face = 0; face < nf; face++) {
        size_t a = (size_t)faces[face * 3];
        size_t b = (size_t)faces[face * 3 + 1];
        size_t c = (size_t)faces[face * 3 + 2];
        if (a >= nv || b >= nv || c >= nv) {
            free(normal);
            return NULL;
        }
        mn_face_add(verts, faces, face, normal);
    }
    for (vertex = 0; vertex < nv; vertex++) mn_vertex_unit(normal, vertex);
    return normal;
}

void MeshNormals_run(const float *verts, const int32_t *faces, size_t f0, size_t f1, size_t v0, size_t v1,
                     float *normal)
{
    size_t face = 0, vertex = 0;
    memset(normal + v0 * 3, 0, (v1 - v0) * 3 * sizeof *normal);
    for (face = f0; face < f1; face++) mn_face_add(verts, faces, face, normal);
    for (vertex = v0; vertex < v1; vertex++) mn_vertex_unit(normal, vertex);
}

float *MeshNormals_compute_runs(const float *verts, size_t nv, const int32_t *faces,
                                const size_t *vofs, const size_t *fofs, size_t nrun)
{
    float *normal = NULL;
    long r = 0;
    int bad = 0;
    if (nrun == 0 || vofs[0] != 0 || fofs[0] != 0 || vofs[nrun] != nv)
        return MeshNormals_compute(verts, nv, faces, nrun ? fofs[nrun] : 0);
    for (r = 0; r < (long)nrun; r++)
        if (vofs[r + 1] < vofs[r] || fofs[r + 1] < fofs[r])
            return MeshNormals_compute(verts, nv, faces, fofs[nrun]);
    normal = (float *)calloc(nv > 0 ? nv * 3 : 1, sizeof(*normal));
    if (normal == NULL) {
        fprintf(stderr, "mesh_normals: out of memory (%zu vertices)\n", nv);
        return NULL;
    }
    /* A run's faces touch only its own vertices, so every vertex receives its faces in file order, as in the
     * serial pass; a face that leaves its run sends the whole mesh down the serial pass instead. */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) reduction(| : bad)
#endif
    for (r = 0; r < (long)nrun; r++) {
        size_t face = 0, vertex = 0, v0 = vofs[r], v1 = vofs[r + 1];
        int k = 0, inside = 1;
        for (face = fofs[r]; face < fofs[r + 1] && inside; face++)
            for (k = 0; k < 3; k++) {
                size_t v = (size_t)faces[face * 3 + (size_t)k];
                if (v < v0 || v >= v1) inside = 0;
            }
        if (!inside) { bad |= 1; continue; }
        for (face = fofs[r]; face < fofs[r + 1]; face++) mn_face_add(verts, faces, face, normal);
        for (vertex = v0; vertex < v1; vertex++) mn_vertex_unit(normal, vertex);
    }
    if (bad) {
        free(normal);
        return MeshNormals_compute(verts, nv, faces, fofs[nrun]);
    }
    return normal;
}

/* The run-parallel normals against the serial pass, bit for bit, on 3000 runs of random fans (with degenerate and
 * repeated faces), then a mesh with one face across two runs (the serial fallback). */
int MeshNormals_selftest(void)
{
    enum { RUNS = 3000 };
    size_t *vofs = (size_t *)malloc((RUNS + 1) * sizeof *vofs), *fofs = (size_t *)malloc((RUNS + 1) * sizeof *fofs);
    float *verts = NULL, *a = NULL, *b = NULL;
    int32_t *faces = NULL;
    size_t nv = 0, nf = 0, r = 0, k = 0;
    uint32_t rng = 0x9e3779b9u;
    int fails = 0, test = 0;
    if (vofs == NULL || fofs == NULL) { free(vofs); free(fofs); return 1; }
    for (r = 0; r < RUNS; r++) {
        vofs[r] = nv; fofs[r] = nf;
        rng = rng * 1664525u + 1013904223u;
        nv += 3 + (rng >> 27);
        nf += 1 + (rng >> 26);
    }
    vofs[RUNS] = nv; fofs[RUNS] = nf;
    verts = (float *)malloc(nv * 3 * sizeof *verts);
    faces = (int32_t *)malloc(nf * 3 * sizeof *faces);
    if (verts == NULL || faces == NULL) fails++;
    for (k = 0; !fails && k < nv * 3; k++) {
        rng = rng * 1664525u + 1013904223u;
        verts[k] = (float)(1000.0 + (double)(rng >> 8) / 65536.0);
    }
    for (r = 0; !fails && r < RUNS; r++)
        for (k = fofs[r]; k < fofs[r + 1]; k++) {
            const size_t span = vofs[r + 1] - vofs[r];
            int c = 0;
            for (c = 0; c < 3; c++) {
                rng = rng * 1664525u + 1013904223u;
                faces[k * 3 + (size_t)c] = (int32_t)(vofs[r] + (rng >> 8) % span);
            }
        }
    for (test = 0; !fails && test < 2; test++) {
        if (test == 1) faces[(fofs[7] + 0) * 3] = (int32_t)vofs[9];     /* a face of run 7 on a vertex of run 9 */
        a = MeshNormals_compute(verts, nv, faces, nf);
        b = MeshNormals_compute_runs(verts, nv, faces, vofs, fofs, RUNS);
        if (a == NULL || b == NULL || memcmp(a, b, nv * 3 * sizeof *a) != 0) fails++;
        if (test == 0 && b != NULL) {   /* the runs one at a time, in a scrambled order, over stale values */
            for (k = 0; k < nv * 3; k++) b[k] = 7.0f;
            for (r = 0; r < RUNS; r++) {
                const size_t rr = (r * 1237) % RUNS;
                MeshNormals_run(verts, faces, fofs[rr], fofs[rr + 1], vofs[rr], vofs[rr + 1], b);
            }
            if (memcmp(a, b, nv * 3 * sizeof *a) != 0) fails++;
        }
        free(a); free(b); a = b = NULL;
    }
    fprintf(stderr, "  mesh_normals selftest: run-parallel normals vs the serial pass (%zu verts, %zu faces): %s\n",
            nv, nf, fails ? "FAIL" : "ok");
    free(vofs); free(fofs); free(verts); free(faces);
    return fails;
}
