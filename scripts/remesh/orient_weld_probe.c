/* orient_weld_probe -- run OrientWeld_components_axis standalone on an OBJ
 * and report what it decides (set ORIENT_WELD_DEBUG=1 for per-component
 * pairs/net/radial lines). Optionally writes the corrected OBJ.
 *
 *   orient_weld_probe <in.obj> [--out <fixed.obj>]
 *                     [--axis-point Z Y X] [--axis-dir Z Y X] [--radius F]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "../../src/common/arena.h"
#include "../../src/remesh/orient_weld.h"

static void *xgrow(void *p, size_t need, size_t *cap, size_t elem)
{
    size_t ncap = *cap;
    if (need <= ncap) return p;
    ncap = (ncap == 0) ? 65536 : ncap;
    while (ncap < need) ncap *= 2;
    p = realloc(p, ncap * elem);
    if (p == NULL) { fprintf(stderr, "OOM\n"); exit(1); }
    *cap = ncap;
    return p;
}

int main(int argc, char **argv)
{
    const char *in_path = NULL, *out_path = NULL;
    float axp[3] = { 0.0f, 3405.0f, 2878.0f };
    float axd[3] = { 1.0f, 0.0f, 0.0f };
    float radius = 3.0f;
    float *verts = NULL;
    int32_t *faces = NULL;
    size_t nv = 0, nf = 0, cap_v = 0, cap_f = 0;
    size_t flipped = 0, radial = 0;
    FILE *f = NULL;
    char line[1024];
    int i = 0;
    Arena_T arena = NULL;

    if (argc < 2) { fprintf(stderr, "usage: %s <in.obj> [--out fixed.obj] "
                            "[--axis-point Z Y X] [--axis-dir Z Y X] "
                            "[--radius F]\n", argv[0]); return 1; }
    in_path = argv[1];
    for (i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--out") == 0 && i + 1 < argc) out_path = argv[++i];
        else if (strcmp(argv[i], "--axis-point") == 0 && i + 3 < argc) {
            axp[0] = (float)atof(argv[++i]);
            axp[1] = (float)atof(argv[++i]);
            axp[2] = (float)atof(argv[++i]);
        } else if (strcmp(argv[i], "--axis-dir") == 0 && i + 3 < argc) {
            axd[0] = (float)atof(argv[++i]);
            axd[1] = (float)atof(argv[++i]);
            axd[2] = (float)atof(argv[++i]);
        } else if (strcmp(argv[i], "--radius") == 0 && i + 1 < argc) {
            radius = (float)atof(argv[++i]);
        }
    }
    f = fopen(in_path, "rb");
    if (f == NULL) { fprintf(stderr, "cannot open %s\n", in_path); return 1; }
    while (fgets(line, sizeof line, f) != NULL) {
        if (line[0] == 'v' && line[1] == ' ') {
            double a = 0, b = 0, c = 0;
            char *s = line + 2, *e = NULL;
            a = strtod(s, &e); s = e; b = strtod(s, &e); s = e; c = strtod(s, &e);
            verts = (float *)xgrow(verts, (nv + 1) * 3, &cap_v, sizeof(float));
            verts[nv*3+0] = (float)a; verts[nv*3+1] = (float)b;
            verts[nv*3+2] = (float)c; nv++;
        } else if (line[0] == 'f' && line[1] == ' ') {
            long idx[8]; int ni = 0;
            char *s = line + 2;
            while (ni < 8) {
                char *e = NULL; long v = 0;
                while (*s == ' ' || *s == '\t') s++;
                if (*s == '\0' || *s == '\n' || *s == '\r') break;
                v = strtol(s, &e, 10);
                if (e == s) break;
                idx[ni++] = v; s = e;
                while (*s && *s != ' ' && *s != '\t' && *s != '\n' && *s != '\r') s++;
            }
            for (i = 2; i < ni; i++) {
                faces = (int32_t *)xgrow(faces, (nf + 1) * 3, &cap_f, sizeof(int32_t));
                faces[nf*3+0] = (int32_t)(idx[0]-1);
                faces[nf*3+1] = (int32_t)(idx[i-1]-1);
                faces[nf*3+2] = (int32_t)(idx[i]-1);
                nf++;
            }
        }
    }
    fclose(f);
    fprintf(stderr, "[probe] %s: %zu verts %zu faces; axis (%.0f,%.0f,%.0f) "
            "radius %.1f\n", in_path, nv, nf,
            (double)axp[0], (double)axp[1], (double)axp[2], (double)radius);

    arena = Arena_new();
    OrientWeld_components_axis(arena, verts, nv, faces, nf, radius,
                               axp, axd, &flipped, &radial);
    fprintf(stderr, "[probe] flipped=%zu radial_decided=%zu\n", flipped, radial);

    if (out_path != NULL) {
        FILE *o = fopen(out_path, "wb");
        size_t k = 0;
        if (o == NULL) { fprintf(stderr, "cannot write %s\n", out_path); return 1; }
        for (k = 0; k < nv; k++)
            fprintf(o, "v %.4f %.4f %.4f\n", (double)verts[k*3+0],
                    (double)verts[k*3+1], (double)verts[k*3+2]);
        for (k = 0; k < nf; k++)
            fprintf(o, "f %d %d %d\n", faces[k*3+0]+1, faces[k*3+1]+1,
                    faces[k*3+2]+1);
        fclose(o);
        fprintf(stderr, "[probe] wrote %s\n", out_path);
    }
    free(verts); free(faces);
    Arena_dispose(&arena);
    return 0;
}
