/* Diagnostic: load an OBJ into a ComponentMesh and run PinholeFill_process,
 * so the PINHOLE_DEBUG instrumentation in pinhole_fill.c prints, per loop it
 * could not fan-fill, how many of its non-adjacent vertex pairs are already
 * mesh edges (blocking chords). Tells us whether the residual micro-holes are
 * fundamentally not fan-fillable (need ear-clipping) or a tracer/guard issue.
 *
 * Build (debug, with the instrumentation):
 *   cl /nologo /DPINHOLE_DEBUG /Fe:output\tools\pinhole_probe.exe /Fo:output\tools\ \
 *      scripts\pinhole_probe.c src\remesh\pinhole_fill.c src\common\arena.c \
 *      src\common\except.c
 *
 * Reads "v z y x" and "f a b c" (1-based). Plain C99. */
#include "../src/remesh/pinhole_fill.h"
#include "../src/common/mesh_types.h"
#include "../src/common/arena.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s mesh.obj\n", argv[0]); return 1; }
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror(argv[1]); return 1; }

    size_t vcap = 1 << 16, nv = 0, fcap = 1 << 16, nf = 0;
    float *V = (float *)malloc(vcap * 3 * sizeof(float));
    int32_t *F = (int32_t *)malloc(fcap * 3 * sizeof(int32_t));
    char line[256];
    while (fgets(line, sizeof line, f)) {
        if (line[0] == 'v' && line[1] == ' ') {
            float z, y, x;
            if (sscanf(line, "v %f %f %f", &z, &y, &x) == 3) {
                if (nv >= vcap) { vcap *= 2; V = (float *)realloc(V, vcap * 3 * sizeof(float)); }
                V[nv*3+0] = z; V[nv*3+1] = y; V[nv*3+2] = x; nv++;
            }
        } else if (line[0] == 'f' && line[1] == ' ') {
            int a, b, c;
            if (sscanf(line, "f %d %d %d", &a, &b, &c) == 3) {
                if (nf >= fcap) { fcap *= 2; F = (int32_t *)realloc(F, fcap * 3 * sizeof(int32_t)); }
                F[nf*3+0]=a-1; F[nf*3+1]=b-1; F[nf*3+2]=c-1; nf++;
            }
        }
    }
    fclose(f);
    fprintf(stderr, "loaded %zu v %zu f\n", nv, nf);

    Arena_T arena = Arena_new();
    /* copy into arena (PinholeFill reallocs from the arena) */
    float *av = (float *)ARENA_ALLOC(arena, (long)(nv * 3 * sizeof(float)));
    int32_t *af = (int32_t *)ARENA_ALLOC(arena, (long)(nf * 3 * sizeof(int32_t)));
    memcpy(av, V, nv * 3 * sizeof(float));
    memcpy(af, F, nf * 3 * sizeof(int32_t));

    ComponentMesh cm; memset(&cm, 0, sizeof cm);
    cm.verts = av; cm.faces = af; cm.nv = nv; cm.nf = nf;
    cm.comp_id = 1; cm.self = &cm;

    size_t sp = 0, fl = 0, ad = 0, sk = 0;
    PinholeFill_process(arena, &cm, 1, 0, &sp, &fl, &ad, &sk);
    fprintf(stderr, "pinhole: splits=%zu filled=%zu tris+=%zu skipped=%zu "
                    "-> %zu faces, %zu verts\n", sp, fl, ad, sk, cm.nf, cm.nv);

    Arena_dispose(&arena);
    free(V); free(F);
    return 0;
}
