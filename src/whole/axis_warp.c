#include "axis_warp.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { double z, y, x; } AxisWarpRow;

void AxisWarp_init(AxisWarp *warp)
{
    if (warp != NULL) memset(warp, 0, sizeof(*warp));
}

void AxisWarp_dispose(AxisWarp *warp)
{
    if (warp == NULL) return;
    free(warp->z);
    free(warp->y);
    free(warp->x);
    AxisWarp_init(warp);
}

int AxisWarp_valid(const AxisWarp *warp)
{
    return warp != NULL && warp->n >= 2 && warp->z != NULL &&
           warp->y != NULL && warp->x != NULL;
}

void AxisWarp_eval(const AxisWarp *warp, double z, double *y, double *x)
{
    size_t lo = 0;
    if (!AxisWarp_valid(warp)) {
        if (y != NULL) *y = 0.0;
        if (x != NULL) *x = 0.0;
        return;
    }
    if (z <= warp->z[0]) lo = 0;
    else if (z >= warp->z[warp->n - 1]) lo = warp->n - 2;
    else {
        size_t hi = warp->n - 1;
        while (hi - lo > 1) {
            size_t mid = lo + (hi - lo) / 2;
            if (warp->z[mid] <= z) lo = mid;
            else hi = mid;
        }
    }
    {
        double dz = warp->z[lo + 1] - warp->z[lo];
        double t = (z - warp->z[lo]) / dz;
        if (y != NULL)
            *y = warp->y[lo] + t * (warp->y[lo + 1] - warp->y[lo]);
        if (x != NULL)
            *x = warp->x[lo] + t * (warp->x[lo + 1] - warp->x[lo]);
    }
}

int AxisWarp_load_csv(AxisWarp *warp, const char *path)
{
    FILE *f;
    AxisWarpRow *row = NULL;
    size_t n = 0, cap = 0;
    char line[1024];
    AxisWarp fresh;
    if (warp == NULL || path == NULL) return -1;
    AxisWarp_init(&fresh);
    f = fopen(path, "r");
    if (f == NULL) return -1;
    while (fgets(line, sizeof(line), f) != NULL) {
        AxisWarpRow r;
        if (sscanf(line, " %lf , %lf , %lf", &r.z, &r.y, &r.x) != 3)
            continue;
        if (!isfinite(r.z) || !isfinite(r.y) || !isfinite(r.x) ||
            (n > 0 && r.z <= row[n - 1].z)) {
            free(row); fclose(f); return -1;
        }
        if (n == cap) {
            size_t next_cap = cap ? cap * 2 : 16;
            AxisWarpRow *next = (AxisWarpRow *)realloc(
                row, next_cap * sizeof(*row));
            if (next == NULL) { free(row); fclose(f); return -1; }
            row = next; cap = next_cap;
        }
        row[n++] = r;
    }
    if (fclose(f) != 0 || n < 2) { free(row); return -1; }
    fresh.z = (double *)malloc(n * sizeof(double));
    fresh.y = (double *)malloc(n * sizeof(double));
    fresh.x = (double *)malloc(n * sizeof(double));
    if (fresh.z == NULL || fresh.y == NULL || fresh.x == NULL) {
        free(row); AxisWarp_dispose(&fresh); return -1;
    }
    fresh.n = n;
    for (size_t i = 0; i < n; i++) {
        fresh.z[i] = row[i].z;
        fresh.y[i] = row[i].y;
        fresh.x[i] = row[i].x;
    }
    free(row);
    AxisWarp_eval(&fresh, 0.5 * (fresh.z[0] + fresh.z[n - 1]),
                  &fresh.reference_y, &fresh.reference_x);
    AxisWarp_dispose(warp);
    *warp = fresh;
    return 0;
}

void AxisWarp_straighten_point(const AxisWarp *warp,
                               const float in_zyx[3], float out_zyx[3])
{
    double y, x;
    if (!AxisWarp_valid(warp)) {
        if (out_zyx != in_zyx) memcpy(out_zyx, in_zyx, 3 * sizeof(float));
        return;
    }
    AxisWarp_eval(warp, (double)in_zyx[0], &y, &x);
    out_zyx[0] = in_zyx[0];
    out_zyx[1] = (float)((double)in_zyx[1] - y + warp->reference_y);
    out_zyx[2] = (float)((double)in_zyx[2] - x + warp->reference_x);
}

void AxisWarp_straighten_vertices(const AxisWarp *warp,
                                  const float *in_zyx, float *out_zyx,
                                  size_t nvertices)
{
    if (in_zyx == NULL || out_zyx == NULL) return;
    for (size_t i = 0; i < nvertices; i++)
        AxisWarp_straighten_point(warp, &in_zyx[i * 3], &out_zyx[i * 3]);
}

int AxisWarp_selftest(void)
{
    AxisWarp w;
    double z[3] = { 0.0, 10.0, 20.0 };
    double y[3] = { 2.0, 7.0, 12.0 };
    double x[3] = { -3.0, -1.0, 1.0 };
    float p[3] = { 15.0f, 11.5f, 2.0f }, q[3];
    double ey, ex;
    int failures = 0;
    AxisWarp_init(&w);
    w.z = z; w.y = y; w.x = x; w.n = 3;
    w.reference_y = 7.0; w.reference_x = -1.0;
    AxisWarp_eval(&w, 15.0, &ey, &ex);
    if (fabs(ey - 9.5) > 1e-12 || fabs(ex) > 1e-12) failures++;
    AxisWarp_straighten_point(&w, p, q);
    if (fabs((double)q[0] - 15.0) > 1e-6 ||
        fabs((double)q[1] - 9.0) > 1e-6 ||
        fabs((double)q[2] - 1.0) > 1e-6) failures++;
    AxisWarp_eval(&w, -10.0, &ey, &ex);
    if (fabs(ey - (-3.0)) > 1e-12 || fabs(ex - (-5.0)) > 1e-12)
        failures++;
    /* Stack-backed samples must not be disposed. */
    AxisWarp_init(&w);
    if (failures != 0)
        fprintf(stderr, "AxisWarp_selftest: FAIL (%d)\n", failures);
    return failures;
}
