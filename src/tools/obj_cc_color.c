/* obj_cc_color.c -- colour an OBJ mesh by CONNECTED COMPONENT + report a
 * fragmentation breakdown. Thin CLI over the shared core in src/common/cc_color.
 *
 * A clean scroll mesh is a small number of large connected components (the
 * wraps): one wrap should be ONE component -> ONE colour. If a decimator has
 * "exploded" the geometry -- shattered a wrap into islands -- that wrap turns
 * into confetti, and the fragmentation summary below (how many components to
 * cover 50/90/99% of faces, dust count) quantifies it. The same coloring is now
 * baked into the pipeline stage dumps and grid_weld's welded.obj.
 *
 * Output verts carry per-vertex RGB ("v x y z r g b"); faces re-emitted without
 * vt/vn. Own OBJ parser + writer (standalone); colours from src/common/cc_color.
 *
 *   obj_cc_color <in.obj> <out.obj> [--min-faces M=0] [--sat S=0.62] [--top N]
 *   obj_cc_color --selftest
 *
 *   --min-faces components with fewer faces than this are coloured dim grey
 *               (dust), so real structure stands out (default 0 = colour all)
 *   --sat       colour saturation 0..1 (default 0.62)
 *   --top       accepted for back-compat (ignored; the summary replaces the table)
 *
 * Exit: 0 ok / selftest pass, 1 IO error, 2 usage error, 3 selftest fail.
 */
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "../common/cc_color.h"

/* ---------- growing arrays ---------- */
typedef struct { float *v; size_t n, cap; } FVec;   /* verts: 3 floats each */
typedef struct { int32_t *f; size_t n, cap; } IVec; /* faces: 3 ints each   */

static int fv_push(FVec *a, float x, float y, float z)
{
    if (a->n == a->cap) { size_t nc = a->cap ? a->cap * 2 : 1u << 16;
        float *np = (float *)realloc(a->v, nc * 3 * sizeof(float)); if (!np) return -1; a->v = np; a->cap = nc; }
    a->v[a->n*3+0] = x; a->v[a->n*3+1] = y; a->v[a->n*3+2] = z; a->n++; return 0;
}
static int iv_push(IVec *a, int32_t x, int32_t y, int32_t z)
{
    if (a->n == a->cap) { size_t nc = a->cap ? a->cap * 2 : 1u << 16;
        int32_t *np = (int32_t *)realloc(a->f, nc * 3 * sizeof(int32_t)); if (!np) return -1; a->f = np; a->cap = nc; }
    a->f[a->n*3+0] = x; a->f[a->n*3+1] = y; a->f[a->n*3+2] = z; a->n++; return 0;
}

static int read_obj(const char *path, FVec *V, IVec *F)
{
    FILE *fp = fopen(path, "rb");
    if (!fp) { fprintf(stderr, "obj_cc_color: cannot open %s\n", path); return -1; }
    size_t lcap = 1u << 16; char *line = (char *)malloc(lcap);
    int rc = 0;
    if (!line) { fclose(fp); return -1; }
    while (fgets(line, (int)lcap, fp)) {
        while (!strchr(line, '\n') && !feof(fp)) {
            size_t len = strlen(line);
            char *nl = (char *)realloc(line, lcap * 2); if (!nl) { rc = -1; goto done; }
            line = nl; lcap *= 2;
            if (!fgets(line + len, (int)(lcap - len), fp)) break;
        }
        if (line[0] == 'v' && (line[1] == ' ' || line[1] == '\t')) {
            double x = 0, y = 0, z = 0;
            if (sscanf(line + 2, "%lf %lf %lf", &x, &y, &z) == 3)
                if (fv_push(V, (float)x, (float)y, (float)z) != 0) { rc = -1; goto done; }
        } else if (line[0] == 'f' && (line[1] == ' ' || line[1] == '\t')) {
            char *tok = strtok(line + 2, " \t\r\n");
            int32_t idx[3]; int n = 0;
            while (tok && n < 3) { idx[n++] = (int32_t)atol(tok); tok = strtok(NULL, " \t\r\n"); }
            if (n == 3) {
                int32_t a = idx[0], b = idx[1], c = idx[2];
                if (a < 0) a = (int32_t)V->n + a + 1;
                if (b < 0) b = (int32_t)V->n + b + 1;
                if (c < 0) c = (int32_t)V->n + c + 1;
                if (iv_push(F, a - 1, b - 1, c - 1) != 0) { rc = -1; goto done; }
            }
        }
    }
done:
    free(line); fclose(fp);
    return rc;
}

static int run(const char *in, const char *out, size_t min_faces, double sat)
{
    FVec V = {0}; IVec F = {0};
    if (read_obj(in, &V, &F) != 0) { free(V.v); free(F.f); return 1; }
    if (V.n == 0 || F.n == 0) { fprintf(stderr, "obj_cc_color: empty mesh\n"); free(V.v); free(F.f); return 1; }

    float *col = (float *)malloc(V.n * 3 * sizeof(float));
    if (!col) { fprintf(stderr, "obj_cc_color: OOM\n"); free(V.v); free(F.f); return 1; }
    CCColorOpts opts = { sat, min_faces };
    CCColorStats st;
    CCColor_compute(V.n, F.f, F.n, &opts, col, &st);

    FILE *ofp = fopen(out, "wb");
    if (!ofp) { fprintf(stderr, "obj_cc_color: cannot write %s\n", out); free(col); free(V.v); free(F.f); return 1; }
    fprintf(ofp, "# obj_cc_color: %zu components, coloured by size rank\n", st.ncomp);
    for (size_t i = 0; i < V.n; i++)
        fprintf(ofp, "v %.5g %.5g %.5g %.4f %.4f %.4f\n",
                V.v[i*3+0], V.v[i*3+1], V.v[i*3+2], col[i*3+0], col[i*3+1], col[i*3+2]);
    for (size_t f = 0; f < F.n; f++)
        fprintf(ofp, "f %d %d %d\n", F.f[f*3+0]+1, F.f[f*3+1]+1, F.f[f*3+2]+1);
    fclose(ofp);

    printf("obj_cc_color: %s -> %s\n", in, out);
    printf("  verts=%zu  faces=%zu  components=%zu\n", st.nv, st.nf, st.ncomp);
    printf("  largest = %zu faces (%.1f%%);  to cover 50%%=%zu  90%%=%zu  99%%=%zu comps\n",
           st.largest_faces, st.nf ? 100.0 * (double)st.largest_faces / (double)st.nf : 0.0,
           st.cover50, st.cover90, st.cover99);
    if (min_faces > 0)
        printf("  dust (< %zu faces): %zu comps holding %zu faces\n", min_faces, st.dust_comps, st.dust_faces);

    free(col); free(V.v); free(F.f);
    return 0;
}

int main(int argc, char **argv)
{
    if (argc == 2 && !strcmp(argv[1], "--selftest")) return CCColor_selftest() ? 3 : 0;
    if (argc < 3) {
        fprintf(stderr, "usage: %s <in.obj> <out.obj> [--min-faces M] [--sat S] [--top N]\n"
                        "       %s --selftest\n", argv[0], argv[0]);
        return 2;
    }
    size_t min_faces = 0; double sat = 0.62;
    for (int i = 3; i < argc; i++) {
        if      (!strcmp(argv[i], "--min-faces") && i + 1 < argc) min_faces = (size_t)atoll(argv[++i]);
        else if (!strcmp(argv[i], "--sat") && i + 1 < argc) sat = atof(argv[++i]);
        else if (!strcmp(argv[i], "--top") && i + 1 < argc) ++i;   /* accepted, ignored */
        else { fprintf(stderr, "obj_cc_color: unknown arg %s\n", argv[i]); return 2; }
    }
    return run(argv[1], argv[2], min_faces, sat);
}
