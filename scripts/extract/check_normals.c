/* check_normals.c — compare face normal consistency pre/post QEM
 *
 * Usage: check_normals <pre_dir> <post_dir> <prefix> <num_components>
 * Example: check_normals output/1013184726/1013184726_step0_pre_simplify \
 *                        output/1013184726/1013184726_step0_post_simplify \
 *                        1013184726 20
 *
 * For each component, loads OBJ, computes area-weighted mean normal,
 * counts faces whose normal opposes the mean (dot < 0).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_VERTS  2000000
#define MAX_FACES  2000000

static float  verts[MAX_VERTS * 3];
static int    faces[MAX_FACES * 3];

static int load_obj(const char *path, int *nv_out, int *nf_out)
{
    FILE *f = fopen(path, "r");
    int nv = 0, nf = 0;
    char line[512];
    if (!f) { fprintf(stderr, "Cannot open %s\n", path); return -1; }

    while (fgets(line, sizeof(line), f)) {
        if (line[0] == 'v' && line[1] == ' ') {
            if (nv >= MAX_VERTS) { fclose(f); return -1; }
            sscanf(line + 2, "%f %f %f",
                   &verts[nv*3+0], &verts[nv*3+1], &verts[nv*3+2]);
            nv++;
        } else if (line[0] == 'f' && line[1] == ' ') {
            if (nf >= MAX_FACES) { fclose(f); return -1; }
            int a, b, c;
            /* OBJ faces may have v/vt/vn format — just grab first int */
            if (sscanf(line + 2, "%d/%*d/%*d %d/%*d/%*d %d/%*d/%*d", &a, &b, &c) == 3 ||
                sscanf(line + 2, "%d//%*d %d//%*d %d//%*d", &a, &b, &c) == 3 ||
                sscanf(line + 2, "%d/%*d %d/%*d %d/%*d", &a, &b, &c) == 3 ||
                sscanf(line + 2, "%d %d %d", &a, &b, &c) == 3) {
                faces[nf*3+0] = a - 1;  /* OBJ is 1-based */
                faces[nf*3+1] = b - 1;
                faces[nf*3+2] = c - 1;
                nf++;
            }
        }
    }
    fclose(f);
    *nv_out = nv;
    *nf_out = nf;
    return 0;
}

static void cross(const float *a, const float *b, float *out)
{
    out[0] = a[1]*b[2] - a[2]*b[1];
    out[1] = a[2]*b[0] - a[0]*b[2];
    out[2] = a[0]*b[1] - a[1]*b[0];
}

static void analyze(int nv, int nf, int *inv_out, double *inv_pct_out,
                    float mean_n[3])
{
    double sum[3] = {0, 0, 0};
    int i;
    double len;

    /* Area-weighted mean normal */
    for (i = 0; i < nf; i++) {
        int a = faces[i*3+0], b = faces[i*3+1], c = faces[i*3+2];
        float e1[3], e2[3], n[3];
        e1[0] = verts[b*3+0] - verts[a*3+0];
        e1[1] = verts[b*3+1] - verts[a*3+1];
        e1[2] = verts[b*3+2] - verts[a*3+2];
        e2[0] = verts[c*3+0] - verts[a*3+0];
        e2[1] = verts[c*3+1] - verts[a*3+1];
        e2[2] = verts[c*3+2] - verts[a*3+2];
        cross(e1, e2, n);
        sum[0] += n[0]; sum[1] += n[1]; sum[2] += n[2];
    }
    len = sqrt(sum[0]*sum[0] + sum[1]*sum[1] + sum[2]*sum[2]);
    if (len > 1e-12) {
        mean_n[0] = (float)(sum[0]/len);
        mean_n[1] = (float)(sum[1]/len);
        mean_n[2] = (float)(sum[2]/len);
    } else {
        mean_n[0] = mean_n[1] = 0.0f; mean_n[2] = 1.0f;
    }

    /* Count inverted faces */
    int inv = 0;
    for (i = 0; i < nf; i++) {
        int a = faces[i*3+0], b = faces[i*3+1], c = faces[i*3+2];
        float e1[3], e2[3], n[3];
        e1[0] = verts[b*3+0] - verts[a*3+0];
        e1[1] = verts[b*3+1] - verts[a*3+1];
        e1[2] = verts[b*3+2] - verts[a*3+2];
        e2[0] = verts[c*3+0] - verts[a*3+0];
        e2[1] = verts[c*3+1] - verts[a*3+1];
        e2[2] = verts[c*3+2] - verts[a*3+2];
        cross(e1, e2, n);
        float dot = n[0]*mean_n[0] + n[1]*mean_n[1] + n[2]*mean_n[2];
        if (dot < 0.0f) inv++;
    }
    *inv_out = inv;
    *inv_pct_out = nf > 0 ? 100.0 * inv / nf : 0.0;

    (void)nv;
}

int main(int argc, char **argv)
{
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <pre_dir> <post_dir> <prefix> <num_comps>\n",
                argv[0]);
        return 1;
    }
    const char *pre_dir = argv[1];
    const char *post_dir = argv[2];
    const char *prefix = argv[3];
    int num_comps = atoi(argv[4]);
    int i;
    char path[1024];

    printf("%-5s  %8s %8s %6s  %8s %8s %6s  %s\n",
           "Comp", "Pre_NV", "Pre_NF", "Inv%",
           "Post_NV", "Post_NF", "Inv%", "MeanN dot");

    for (i = 1; i <= num_comps; i++) {
        int pre_nv = 0, pre_nf = 0, post_nv = 0, post_nf = 0;
        int pre_inv = 0, post_inv = 0;
        double pre_pct = 0, post_pct = 0;
        float pre_mn[3] = {0}, post_mn[3] = {0};
        float ndot;

        snprintf(path, sizeof(path), "%s/%s_pre_simplify_comp%03d.obj",
                 pre_dir, prefix, i);
        if (load_obj(path, &pre_nv, &pre_nf) != 0) {
            printf("%3d    MISSING PRE\n", i);
            continue;
        }
        analyze(pre_nv, pre_nf, &pre_inv, &pre_pct, pre_mn);

        snprintf(path, sizeof(path), "%s/%s_post_simplify_comp%03d.obj",
                 post_dir, prefix, i);
        if (load_obj(path, &post_nv, &post_nf) != 0) {
            printf("%3d    MISSING POST\n", i);
            continue;
        }
        analyze(post_nv, post_nf, &post_inv, &post_pct, post_mn);

        ndot = pre_mn[0]*post_mn[0] + pre_mn[1]*post_mn[1] + pre_mn[2]*post_mn[2];

        printf("%3d    %7d %7d %5.1f%%  %7d %7d %5.1f%%  %.4f%s\n",
               i, pre_nv, pre_nf, pre_pct, post_nv, post_nf, post_pct,
               ndot, ndot < 0.0f ? "  *** FLIPPED ***" : "");
    }

    return 0;
}
