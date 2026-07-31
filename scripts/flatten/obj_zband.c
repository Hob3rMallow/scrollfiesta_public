/* obj_zband -- extract faces of an OBJ within a z (or any-axis) coordinate
 * band, writing a 3D OBJ and, if the input carries per-vertex vt, a flat
 * (u,v,0) OBJ. For localizing seam/slice defects (e.g. the ribbon u-smear
 * slab at a cube seam). Vertices are all emitted; only in-band faces kept.
 *
 *   obj_zband <in.obj> <zlo> <zhi> <out3d.obj> [outflat.obj] [--axis 0|1|2]
 * axis 0 = z (default, first vertex coord in our z,y,x order), 1 = y, 2 = x.
 * A face is in-band if its centroid on that axis is within [zlo,zhi].
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static float *V = NULL, *T = NULL;
static long *F = NULL;
static size_t nv = 0, nvt = 0, nf = 0, capv = 0, capt = 0, capf = 0;

static void *grow(void *p, size_t need, size_t *cap, size_t el)
{
    size_t c = *cap;
    if (need <= c) return p;
    c = c ? c : 4096;
    while (c < need) c *= 2;
    p = realloc(p, c * el);
    if (!p) { fprintf(stderr, "oom\n"); exit(1); }
    *cap = c; return p;
}

int main(int argc, char **argv)
{
    const char *in = NULL, *o3 = NULL, *of = NULL;
    double zlo = 0, zhi = 0;
    int axis = 0;
    FILE *fp = NULL, *w = NULL;
    char line[1024];
    size_t i = 0, kept = 0;
    int a = 0;

    if (argc < 5) {
        fprintf(stderr, "usage: %s <in.obj> <zlo> <zhi> <out3d.obj> "
                "[outflat.obj] [--axis 0|1|2]\n", argv[0]);
        return 1;
    }
    in = argv[1]; zlo = atof(argv[2]); zhi = atof(argv[3]); o3 = argv[4];
    for (a = 5; a < argc; a++) {
        if (strcmp(argv[a], "--axis") == 0 && a + 1 < argc) axis = atoi(argv[++a]);
        else of = argv[a];
    }
    fp = fopen(in, "rb");
    if (!fp) { fprintf(stderr, "cannot open %s\n", in); return 1; }
    while (fgets(line, sizeof line, fp)) {
        if (line[0] == 'v' && line[1] == ' ') {
            char *s = line + 2, *e = NULL;
            double x = strtod(s, &e); s = e;
            double y = strtod(s, &e); s = e;
            double z = strtod(s, &e);
            V = grow(V, (nv + 1) * 3, &capv, sizeof(float));
            V[nv*3+0] = (float)x; V[nv*3+1] = (float)y; V[nv*3+2] = (float)z;
            nv++;
        } else if (line[0] == 'v' && line[1] == 't' && line[2] == ' ') {
            char *s = line + 3, *e = NULL;
            double u = strtod(s, &e); s = e;
            double v = strtod(s, &e);
            T = grow(T, (nvt + 1) * 2, &capt, sizeof(float));
            T[nvt*2+0] = (float)u; T[nvt*2+1] = (float)v; nvt++;
        } else if (line[0] == 'f' && line[1] == ' ') {
            long idx[16]; int ni = 0, k = 0;
            char *s = line + 2;
            while (ni < 16) {
                char *e = NULL; long vi = 0;
                while (*s == ' ' || *s == '\t') s++;
                if (*s == '\0' || *s == '\n' || *s == '\r') break;
                vi = strtol(s, &e, 10);
                if (e == s) break;
                idx[ni++] = vi; s = e;
                while (*s && *s != ' ' && *s != '\t' && *s != '\n' && *s != '\r') s++;
            }
            for (k = 2; k < ni; k++) {
                F = grow(F, (nf + 1) * 3, &capf, sizeof(long));
                F[nf*3+0] = idx[0]-1; F[nf*3+1] = idx[k-1]-1; F[nf*3+2] = idx[k]-1;
                nf++;
            }
        }
    }
    fclose(fp);
    fprintf(stderr, "loaded %zu v, %zu vt, %zu f; band axis %d [%.1f,%.1f]\n",
            nv, nvt, nf, axis, zlo, zhi);

    w = fopen(o3, "wb");
    if (!w) { fprintf(stderr, "cannot write %s\n", o3); return 1; }
    for (i = 0; i < nv; i++)
        fprintf(w, "v %.4f %.4f %.4f\n", V[i*3+0], V[i*3+1], V[i*3+2]);
    for (i = 0; i < nf; i++) {
        double c = (V[(size_t)F[i*3+0]*3+axis] + V[(size_t)F[i*3+1]*3+axis]
                    + V[(size_t)F[i*3+2]*3+axis]) / 3.0;
        if (c < zlo || c > zhi) continue;
        fprintf(w, "f %ld %ld %ld\n", F[i*3+0]+1, F[i*3+1]+1, F[i*3+2]+1);
        kept++;
    }
    fclose(w);
    fprintf(stderr, "wrote %s (%zu faces in band)\n", o3, kept);

    if (of && nvt == nv) {
        w = fopen(of, "wb");
        if (w) {
            for (i = 0; i < nv; i++)
                fprintf(w, "v %.4f %.4f 0\n", T[i*2+0], T[i*2+1]);
            for (i = 0; i < nf; i++) {
                double c = (V[(size_t)F[i*3+0]*3+axis] + V[(size_t)F[i*3+1]*3+axis]
                            + V[(size_t)F[i*3+2]*3+axis]) / 3.0;
                if (c < zlo || c > zhi) continue;
                fprintf(w, "f %ld %ld %ld\n", F[i*3+0]+1, F[i*3+1]+1, F[i*3+2]+1);
            }
            fclose(w);
            fprintf(stderr, "wrote %s (flat u,v)\n", of);
        }
    }
    free(V); free(T); free(F);
    return 0;
}
