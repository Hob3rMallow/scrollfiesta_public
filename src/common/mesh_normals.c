#include "mesh_normals.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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
        double e1[3], e2[3], cross[3];
        int axis = 0;
        if (a >= nv || b >= nv || c >= nv) {
            free(normal);
            return NULL;
        }
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
    for (vertex = 0; vertex < nv; vertex++) {
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
    return normal;
}
