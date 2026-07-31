/* Local PCA thickness analyzer for LOP point-cloud dumps.
 *
 * For each vertex, finds neighbors within radius R, computes weighted
 * PCA, and reports the smallest-eigenvalue stddev (the "thickness").
 *
 * Single-layer collapsed: thickness < 0.5 vox
 * Partial collapse:       thickness 0.5-1.5 vox
 * Failed collapse:        thickness > 1.5 vox (still two layers visible)
 *
 * Reports a histogram so we can tell whether LOP collapsed most of the
 * cloud and there's a tail of failures, or if the whole thing is uncolla-
 * psed.
 *
 * Usage: lop_thickness_analyze <obj-file> [radius_vox]
 */
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct { float z, y, x; } Vec3;

static int cmp_float(const void *a, const void *b)
{
    float fa = *(const float*)a, fb = *(const float*)b;
    if (fa < fb) return -1;
    if (fa > fb) return 1;
    return 0;
}

static void jacobi_3x3(double a[3][3], double evals[3])
{
    double m[3][3];
    memcpy(m, a, sizeof(m));
    for (int iter = 0; iter < 50; iter++) {
        double max_off = 0;
        int p = 0, q = 1;
        for (int i = 0; i < 3; i++)
            for (int j = i+1; j < 3; j++) {
                double v = fabs(m[i][j]);
                if (v > max_off) { max_off = v; p = i; q = j; }
            }
        if (max_off < 1e-14) break;
        double theta, diff = m[q][q] - m[p][p];
        if (fabs(diff) < 1e-30) theta = M_PI/4.0;
        else theta = 0.5 * atan2(2.0 * m[p][q], diff);
        double c = cos(theta), s = sin(theta);
        double app = m[p][p], aqq = m[q][q], apq = m[p][q];
        m[p][p] = c*c*app - 2*s*c*apq + s*s*aqq;
        m[q][q] = s*s*app + 2*s*c*apq + c*c*aqq;
        m[p][q] = 0; m[q][p] = 0;
        for (int r = 0; r < 3; r++) {
            if (r != p && r != q) {
                double rp = m[r][p], rq = m[r][q];
                m[r][p] = c*rp - s*rq; m[p][r] = m[r][p];
                m[r][q] = s*rp + c*rq; m[q][r] = m[r][q];
            }
        }
    }
    evals[0] = m[0][0]; evals[1] = m[1][1]; evals[2] = m[2][2];
    /* Sort ascending (insertion sort, 3 elems). */
    for (int i = 1; i < 3; i++) {
        double v = evals[i]; int j = i-1;
        while (j >= 0 && evals[j] > v) { evals[j+1] = evals[j]; j--; }
        evals[j+1] = v;
    }
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "usage: %s <obj-file> [radius_vox=6.0]\n", argv[0]);
        return 1;
    }
    float R = (argc >= 3) ? (float)atof(argv[2]) : 6.0f;

    FILE *f = fopen(argv[1], "r");
    if (!f) { perror("fopen"); return 1; }

    /* Read all verts. */
    size_t cap = 1<<16, n = 0;
    Vec3 *V = malloc(cap * sizeof(*V));
    char line[256];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] != 'v' || line[1] != ' ') continue;
        Vec3 v;
        if (sscanf(line, "v %f %f %f", &v.z, &v.y, &v.x) != 3) continue;
        if (n >= cap) { cap *= 2; V = realloc(V, cap * sizeof(*V)); }
        V[n++] = v;
    }
    fclose(f);
    fprintf(stderr, "Loaded %zu verts from %s (R=%.1f vox)\n", n, argv[1], (double)R);

    /* Build voxel grid for fast neighbor lookup. */
    float zmin=V[0].z, zmax=V[0].z, ymin=V[0].y, ymax=V[0].y, xmin=V[0].x, xmax=V[0].x;
    for (size_t i = 1; i < n; i++) {
        if (V[i].z < zmin) zmin = V[i].z; if (V[i].z > zmax) zmax = V[i].z;
        if (V[i].y < ymin) ymin = V[i].y; if (V[i].y > ymax) ymax = V[i].y;
        if (V[i].x < xmin) xmin = V[i].x; if (V[i].x > xmax) xmax = V[i].x;
    }
    int gz = (int)ceilf((zmax - zmin) / R) + 1;
    int gy = (int)ceilf((ymax - ymin) / R) + 1;
    int gx = (int)ceilf((xmax - xmin) / R) + 1;
    size_t ncells = (size_t)gz * (size_t)gy * (size_t)gx;
    /* Cell linked list: head[cell] -> idx -> next[idx] -> ... */
    int *head = malloc(ncells * sizeof(*head));
    int *next = malloc(n * sizeof(*next));
    for (size_t i = 0; i < ncells; i++) head[i] = -1;
    for (size_t i = 0; i < n; i++) {
        int cz = (int)floorf((V[i].z - zmin) / R);
        int cy = (int)floorf((V[i].y - ymin) / R);
        int cx = (int)floorf((V[i].x - xmin) / R);
        size_t c = (size_t)cz * gy * gx + (size_t)cy * gx + (size_t)cx;
        next[i] = head[c];
        head[c] = (int)i;
    }

    /* For each vertex, compute local PCA. */
    double R2 = R * R;
    int bins[20] = {0};   /* thickness stddev in 0.1-vox buckets, 0..2.0 */
    int total = 0;
    double thick_sum = 0;

    /* Sample 1 in N for speed. */
    size_t stride = (n > 5000) ? n / 5000 : 1;
    for (size_t i = 0; i < n; i += stride) {
        Vec3 vi = V[i];
        int cz = (int)floorf((vi.z - zmin) / R);
        int cy = (int)floorf((vi.y - ymin) / R);
        int cx = (int)floorf((vi.x - xmin) / R);

        double cz_s = 0, cy_s = 0, cx_s = 0, w_sum = 0;
        size_t n_used = 0;
        for (int dz = -1; dz <= 1; dz++)
            for (int dy = -1; dy <= 1; dy++)
                for (int dx = -1; dx <= 1; dx++) {
                    int qz = cz+dz, qy = cy+dy, qx = cx+dx;
                    if (qz<0||qz>=gz||qy<0||qy>=gy||qx<0||qx>=gx) continue;
                    size_t c = (size_t)qz*gy*gx + (size_t)qy*gx + qx;
                    for (int j = head[c]; j != -1; j = next[j]) {
                        double dvz = V[j].z - vi.z, dvy = V[j].y - vi.y, dvx = V[j].x - vi.x;
                        double r2 = dvz*dvz + dvy*dvy + dvx*dvx;
                        if (r2 >= R2) continue;
                        double w = 1.0;  /* uniform weight for diagnostic */
                        w_sum += w;
                        cz_s += w * V[j].z;
                        cy_s += w * V[j].y;
                        cx_s += w * V[j].x;
                        n_used++;
                    }
                }
        if (n_used < 8) continue;
        double inv_w = 1.0 / w_sum;
        double mz = cz_s * inv_w, my = cy_s * inv_w, mx = cx_s * inv_w;
        double m[3][3] = {{0}};
        for (int dz = -1; dz <= 1; dz++)
            for (int dy = -1; dy <= 1; dy++)
                for (int dx = -1; dx <= 1; dx++) {
                    int qz = cz+dz, qy = cy+dy, qx = cx+dx;
                    if (qz<0||qz>=gz||qy<0||qy>=gy||qx<0||qx>=gx) continue;
                    size_t c = (size_t)qz*gy*gx + (size_t)qy*gx + qx;
                    for (int j = head[c]; j != -1; j = next[j]) {
                        double dvz = V[j].z - vi.z, dvy = V[j].y - vi.y, dvx = V[j].x - vi.x;
                        double r2 = dvz*dvz + dvy*dvy + dvx*dvx;
                        if (r2 >= R2) continue;
                        double pz = V[j].z - mz, py = V[j].y - my, px = V[j].x - mx;
                        m[0][0] += pz*pz; m[0][1] += pz*py; m[0][2] += pz*px;
                        m[1][1] += py*py; m[1][2] += py*px;
                        m[2][2] += px*px;
                    }
                }
        m[1][0] = m[0][1]; m[2][0] = m[0][2]; m[2][1] = m[1][2];
        for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) m[a][b] *= inv_w;
        double ev[3];
        jacobi_3x3(m, ev);
        double thick = sqrt(fmax(ev[0], 0.0));
        thick_sum += thick;
        int b = (int)(thick / 0.1);
        if (b < 0) b = 0; if (b >= 20) b = 19;
        bins[b]++;
        total++;
    }

    printf("Local thickness (smallest-eigenvalue stddev) histogram over %d sampled verts:\n", total);
    printf("(neighbors within R=%.1f vox, weights uniform)\n\n", (double)R);
    printf("  bin (vox)    count   %%\n");
    for (int b = 0; b < 20; b++) {
        if (bins[b] == 0) continue;
        printf("  %.1f - %.1f    %5d  %5.1f%%\n", b*0.1, (b+1)*0.1, bins[b], 100.0*bins[b]/total);
    }
    printf("\nMean thickness: %.3f vox\n", thick_sum / total);

    free(head); free(next); free(V);
    return 0;
}
