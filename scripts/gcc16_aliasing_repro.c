/* Reduced testcase: gcc 16.1.0 -O3 miscompiles this correct, UB-free
 * interleaved-store loop unless -fno-strict-aliasing (or -O2) is used.
 *
 *   gcc -O3 -o a.exe gcc16_aliasing_repro.c && ./a.exe
 *   gcc -O3 -fno-strict-aliasing -o b.exe gcc16_aliasing_repro.c && ./b.exe
 *
 * The two runs print different face hashes / arrays; the -fno-strict-aliasing
 * output is the correct one (slot 26 must be 12: cell (i=4,j=0) tri0 is
 * (4,5,12) on an 8x8 grid). Verified with winlibs and MSYS2 gcc 16.1.0,
 * x86_64-w64-mingw32, July 2026. This is why the build adds
 * -fno-strict-aliasing for GNU compilers (CMakeLists.txt / src/Makefile).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef struct { float *vertices; int32_t *faces; size_t n_vertices, n_faces; int32_t *vmap; } sf_mesh;
static void make_grid_sheet(int n, sf_mesh *m)
{
    size_t nv = (size_t)n * (size_t)n;
    size_t nf = 2u * (size_t)(n - 1) * (size_t)(n - 1);
    memset(m, 0, sizeof *m);
    m->vertices = malloc(nv * 3 * sizeof(float));
    m->faces    = malloc(nf * 3 * sizeof(int32_t));
    m->n_vertices = nv;
    m->n_faces    = nf;
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) {
            size_t v = (size_t)j * (size_t)n + (size_t)i;
            m->vertices[v * 3 + 0] = (float)i;
            m->vertices[v * 3 + 1] = (float)j;
            m->vertices[v * 3 + 2] = 0.0f;
        }
    size_t f = 0;
    for (int j = 0; j + 1 < n; j++)
        for (int i = 0; i + 1 < n; i++) {
            int32_t v00 = j * n + i, v10 = j * n + i + 1;
            int32_t v01 = (j + 1) * n + i, v11 = (j + 1) * n + i + 1;
            m->faces[f * 3 + 0] = v00; m->faces[f * 3 + 1] = v10; m->faces[f * 3 + 2] = v01; f++;
            m->faces[f * 3 + 0] = v10; m->faces[f * 3 + 1] = v11; m->faces[f * 3 + 2] = v01; f++;
        }
}
int main(void) {
    sf_mesh m; make_grid_sheet(8, &m);
    unsigned long hv = 5381, hf = 5381;
    for (size_t i = 0; i < m.n_vertices * 3; i++) hv = hv * 33 + (unsigned)(m.vertices[i] * 7);
    for (size_t i = 0; i < m.n_faces * 3; i++) hf = hf * 33 + (unsigned)m.faces[i];
    printf("nv=%zu nf=%zu hv=%lu hf=%lu\n", m.n_vertices, m.n_faces, hv, hf);
    for (size_t i = 0; i < 12; i++) printf("%d ", m.faces[i]);
    printf("\n");
    return 0;
}
