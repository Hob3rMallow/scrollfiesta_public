/* asm_axis.c -- the scroll axis as a polyline (see asm_axis.h). */
#include "asm_axis.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define AA_DENSE_DZ 8.0       /* z step of the densified polyline, vox */
#define AA_SMOOTH_SIGMA 640.0 /* the LOESS window, vox: the derived table carries ~6-vox structure at 200-400 vox wavelengths (curvature radius ~400 vox, inside the scroll); at 640 the curvature radius is thousands of vox and the frame is rigid across the sheet */

static int aa_cmp_row(const void *x, const void *y)
{
    const double *a = x, *b = y;
    return (a[0] > b[0]) - (a[0] < b[0]);
}

/* arc length per row, and a CONTINUOUS unit tangent per row (the mean of the two adjacent
 * segments): the frame must not jump at a row, or every lattice row crossing it reads a step of
 * r x (tangent change) -- 8 vox at r 500 for one degree, 50,253 edges over the gate on the
 * familiar (dev60, 2026-09-09) */
static void aa_finish(AsmAxis *a)
{
    a->s[0] = 0.0;
    for (size_t i = 1; i < a->n; i++) {
        double dz = a->p[i*3] - a->p[(i-1)*3], dy = a->p[i*3+1] - a->p[(i-1)*3+1], dx = a->p[i*3+2] - a->p[(i-1)*3+2];
        a->s[i] = a->s[i-1] + sqrt(dz*dz + dy*dy + dx*dx);
    }
    for (size_t i = 0; i < a->n; i++) {
        double t[3] = { 0.0, 0.0, 0.0 };
        for (int side = 0; side < 2; side++) {
            size_t j0 = side == 0 ? (i > 0 ? i - 1 : i) : i, j1 = j0 + 1;
            if (j1 >= a->n) continue;
            if (side == 0 && i == 0) continue;
            double e[3] = { a->p[j1*3] - a->p[j0*3], a->p[j1*3+1] - a->p[j0*3+1], a->p[j1*3+2] - a->p[j0*3+2] };
            double L = sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
            if (L < 1e-12) continue;
            for (int k = 0; k < 3; k++) t[k] += e[k] / L;
        }
        double L = sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]);
        if (L < 1e-12) { t[0] = 1.0; t[1] = 0.0; t[2] = 0.0; L = 1.0; }
        for (int k = 0; k < 3; k++) a->t[i*3+k] = t[k] / L;
    }
}

AsmAxis *AsmAxis_load(Arena_T arena, const char *path)
{
    if (path == NULL || path[0] == 0) return NULL;
    FILE *f = fopen(path, "r");
    if (f == NULL) return NULL;
    size_t cap = 256, n = 0;
    double *rows = ARENA_ALLOC(arena, cap * 3 * sizeof(double));
    char line[512];
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '#' || line[0] == 0 || line[0] == '\n' || line[0] == '\r') continue;
        double z, y, x;
        if (sscanf(line, "%lf,%lf,%lf", &z, &y, &x) != 3) continue;
        if (n == cap) { double *nr = ARENA_ALLOC(arena, cap * 2 * 3 * sizeof(double)); memcpy(nr, rows, n * 3 * sizeof(double)); rows = nr; cap *= 2; }
        rows[n*3] = z; rows[n*3+1] = y; rows[n*3+2] = x; n++;
    }
    fclose(f);
    return AsmAxis_from_rows(arena, rows, NULL, n);
}

AsmAxis *AsmAxis_from_rows(Arena_T arena, const double *in_rows, const double *in_w, size_t n)
{
    if (n < 2) return NULL;
    /* sort the rows by z (their weights ride along), drop duplicate z */
    double *rows = ARENA_ALLOC(arena, n * 4 * sizeof(double));
    for (size_t i = 0; i < n; i++) { rows[i*4] = in_rows[i*3]; rows[i*4+1] = in_rows[i*3+1]; rows[i*4+2] = in_rows[i*3+2]; rows[i*4+3] = in_w ? in_w[i] : 1.0; }
    qsort(rows, n, 4 * sizeof(double), aa_cmp_row);
    size_t w = 0;
    for (size_t i = 0; i < n; i++) {
        if (w > 0 && fabs(rows[i*4] - rows[(w-1)*4]) < 1e-9) continue;
        if (w != i) memcpy(&rows[w*4], &rows[i*4], 4 * sizeof(double));
        w++;
    }
    n = w;
    if (n < 2) return NULL;
    double *rw = ARENA_ALLOC(arena, n * sizeof(double));
    for (size_t i = 0; i < n; i++) { rw[i] = rows[i*4+3]; }
    /* back to (z,y,x) triples for the fit below */
    double *r3 = ARENA_ALLOC(arena, n * 3 * sizeof(double));
    for (size_t i = 0; i < n; i++) { r3[i*3] = rows[i*4]; r3[i*3+1] = rows[i*4+1]; r3[i*3+2] = rows[i*4+2]; }
    rows = r3;
    /* DENSIFY by LOESS: at every AA_DENSE_DZ step a Gaussian-weighted (sigma AA_SMOOTH_SIGMA) local
     * quadratic through the rows gives y(z), x(z).  It follows the table's real wander (a cubic
     * missed the derived table by 10 vox rms) and removes the sub-vox wiggles between rows whose
     * curvature radius (~400 vox for 1 vox over 64) folds the nearest-point frame inside the scroll. */
    double zlo = rows[0], zhi = rows[(n-1)*3];
    size_t nd = (size_t)((zhi - zlo) / AA_DENSE_DZ) + 2;
    if (nd < 2) nd = 2;
    double *dense = ARENA_ALLOC(arena, nd * 3 * sizeof(double));
    for (size_t k = 0; k < nd; k++) {
        double z = zlo + (zhi - zlo) * (double)k / (double)(nd - 1);
        dense[k*3] = z;
        for (int c = 0; c < 2; c++) {
            double M[9] = { 0 }, R[3] = { 0 };
            for (size_t i = 0; i < n; i++) {
                double dz = (rows[i*3] - z) / AA_SMOOTH_SIGMA, wt = exp(-0.5 * dz * dz) * rw[i];
                double f[3] = { 1.0, dz, dz * dz };
                for (int r = 0; r < 3; r++) { for (int q = 0; q < 3; q++) M[r*3+q] += wt * f[r] * f[q]; R[r] += wt * f[r] * rows[i*3+1+c]; }
            }
            /* solve the 3x3 (a local quadratic needs three rows in reach; fall back to the mean) */
            double val = 0.0;
            if (M[0] > 1e-9) {
                double A[9]; double B[3];
                memcpy(A, M, sizeof A); memcpy(B, R, sizeof B);
                int fine = 1;
                for (int col = 0; col < 3 && fine; col++) {
                    int piv = col;
                    for (int r = col + 1; r < 3; r++) if (fabs(A[r*3+col]) > fabs(A[piv*3+col])) piv = r;
                    for (int q = 0; q < 3; q++) { double t = A[col*3+q]; A[col*3+q] = A[piv*3+q]; A[piv*3+q] = t; }
                    { double t = B[col]; B[col] = B[piv]; B[piv] = t; }
                    if (fabs(A[col*3+col]) < 1e-9 * M[0]) { fine = 0; break; }
                    for (int r = 0; r < 3; r++) {
                        if (r == col) continue;
                        double fct = A[r*3+col] / A[col*3+col];
                        for (int q = 0; q < 3; q++) A[r*3+q] -= fct * A[col*3+q];
                        B[r] -= fct * B[col];
                    }
                }
                val = fine ? B[0] / A[0] : R[0] / M[0];
            }
            dense[k*3+1+c] = val;
        }
    }
    AsmAxis *a = ARENA_ALLOC(arena, sizeof(AsmAxis));
    a->n = nd; a->p = dense; a->s = ARENA_ALLOC(arena, nd * sizeof(double)); a->t = ARENA_ALLOC(arena, nd * 3 * sizeof(double)); a->is_line = 0;
    a->n_rows = n;
    aa_finish(a);
    /* the fit's residual against the rows and the polyline's minimum curvature radius */
    a->fit_rms = 0.0; a->r_curv_min = 1e300;
    for (size_t i = 0; i < n; i++) {
        double q[3] = { rows[i*3], rows[i*3+1], rows[i*3+2] }, rr = 0.0;
        AsmAxis_project(a, q, NULL, &rr, NULL, NULL);
        a->fit_rms += rr * rr;
    }
    a->fit_rms = sqrt(a->fit_rms / (double)n);
    for (size_t i = 1; i + 1 < nd; i++) {
        double d1[3] = { dense[(i+1)*3]-dense[i*3], dense[(i+1)*3+1]-dense[i*3+1], dense[(i+1)*3+2]-dense[i*3+2] };
        double d0[3] = { dense[i*3]-dense[(i-1)*3], dense[i*3+1]-dense[(i-1)*3+1], dense[i*3+2]-dense[(i-1)*3+2] };
        double l1 = sqrt(d1[0]*d1[0]+d1[1]*d1[1]+d1[2]*d1[2]), l0 = sqrt(d0[0]*d0[0]+d0[1]*d0[1]+d0[2]*d0[2]);
        if (l0 < 1e-9 || l1 < 1e-9) continue;
        double cosang = (d0[0]*d1[0]+d0[1]*d1[1]+d0[2]*d1[2]) / (l0 * l1);
        if (cosang > 1.0) cosang = 1.0;
        if (cosang < -1.0) cosang = -1.0;
        double ang = acos(cosang);
        if (ang > 1e-9) { double rc = 0.5 * (l0 + l1) / ang; if (rc < a->r_curv_min) a->r_curv_min = rc; }
    }
    return a;
}

AsmAxis *AsmAxis_line(Arena_T arena, const double point[3], const double dir[3], double z0, double z1)
{
    AsmAxis *a = ARENA_ALLOC(arena, sizeof(AsmAxis));
    a->n = 2; a->p = ARENA_ALLOC(arena, 6 * sizeof(double)); a->s = ARENA_ALLOC(arena, 2 * sizeof(double)); a->t = ARENA_ALLOC(arena, 6 * sizeof(double)); a->is_line = 1; a->n_rows = 2;
    a->fit_rms = 0.0; a->r_curv_min = 1e300;
    double dz = dir[0];
    if (fabs(dz) < 1e-9) dz = 1e-9;
    double t0 = (z0 - point[0]) / dz, t1 = (z1 - point[0]) / dz;
    if (t1 <= t0) t1 = t0 + 1.0;
    for (int k = 0; k < 3; k++) { a->p[k] = point[k] + t0 * dir[k]; a->p[3+k] = point[k] + t1 * dir[k]; }
    aa_finish(a);
    return a;
}

void AsmAxis_project(const AsmAxis *a, const double p[3], double *s, double *r, double dir[3], double foot[3])
{
    double best_d2 = 1e300, best_s = 0.0, best_dir[3] = { 1.0, 0.0, 0.0 }, best_foot[3] = { p[0], p[1], p[2] };
    /* the rows are monotone in z: only the segments within a window of the point's z can be nearest
     * (the window covers a tilt of 45 degrees over the largest radius the scrolls have) */
    size_t lo = 0, hi = a->n - 1;
    if (a->n > 64) {
        size_t l = 0, h = a->n - 1;
        while (l < h) { size_t mid = (l + h) / 2; if (a->p[mid*3] < p[0]) l = mid + 1; else h = mid; }
        double win = 2000.0;
        size_t span = (size_t)(win / AA_DENSE_DZ) + 2;
        lo = l > span ? l - span : 0;
        hi = l + span < a->n - 1 ? l + span : a->n - 1;
    }
    for (size_t i = lo; i + 1 <= hi && i + 1 < a->n; i++) {
        const double *A = &a->p[i*3], *B = &a->p[(i+1)*3];
        double e[3] = { B[0]-A[0], B[1]-A[1], B[2]-A[2] };
        double L2 = e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
        if (L2 < 1e-18) continue;
        double t = ((p[0]-A[0])*e[0] + (p[1]-A[1])*e[1] + (p[2]-A[2])*e[2]) / L2;
        if (i > 0 && t < 0.0) t = 0.0;            /* extrapolate only beyond the ends */
        if (i + 2 < a->n && t > 1.0) t = 1.0;
        double q[3] = { A[0] + t*e[0], A[1] + t*e[1], A[2] + t*e[2] };
        double d2 = (p[0]-q[0])*(p[0]-q[0]) + (p[1]-q[1])*(p[1]-q[1]) + (p[2]-q[2])*(p[2]-q[2]);
        if (d2 < best_d2) {
            double L = sqrt(L2);
            best_d2 = d2; best_s = a->s[i] + t * L;
            /* the tangent interpolated between the rows' continuous tangents (the end rows beyond the ends) */
            double tt = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
            double dd[3] = { (1.0 - tt) * a->t[i*3] + tt * a->t[(i+1)*3], (1.0 - tt) * a->t[i*3+1] + tt * a->t[(i+1)*3+1], (1.0 - tt) * a->t[i*3+2] + tt * a->t[(i+1)*3+2] };
            double dl = sqrt(dd[0]*dd[0] + dd[1]*dd[1] + dd[2]*dd[2]);
            if (dl < 1e-12) { dd[0] = e[0] / L; dd[1] = e[1] / L; dd[2] = e[2] / L; dl = 1.0; }
            for (int k = 0; k < 3; k++) { best_dir[k] = dd[k] / dl; best_foot[k] = q[k]; }
        }
    }
    if (s) *s = best_s;
    if (r) *r = sqrt(best_d2);
    if (dir) for (int k = 0; k < 3; k++) dir[k] = best_dir[k];
    if (foot) for (int k = 0; k < 3; k++) foot[k] = best_foot[k];
}

void AsmAxis_frame(const AsmAxis *a, const double p[3], double out[3])
{
    double s, dir[3], foot[3];
    AsmAxis_project(a, p, &s, NULL, dir, foot);
    /* e1 = world +y made perpendicular to the tangent; e2 = dir x e1 (right-handed (dir, e1, e2)) */
    double e1[3] = { -dir[1]*dir[0], 1.0 - dir[1]*dir[1], -dir[1]*dir[2] };
    double l = sqrt(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]);
    if (l < 1e-9) { e1[0] = 0.0; e1[1] = 0.0; e1[2] = 1.0; l = 1.0; }   /* the tangent along +y: use +x */
    for (int k = 0; k < 3; k++) e1[k] /= l;
    double e2[3] = { dir[1]*e1[2] - dir[2]*e1[1], dir[2]*e1[0] - dir[0]*e1[2], dir[0]*e1[1] - dir[1]*e1[0] };
    double d[3] = { p[0]-foot[0], p[1]-foot[1], p[2]-foot[2] };
    out[0] = s;
    out[1] = d[0]*e1[0] + d[1]*e1[1] + d[2]*e1[2];
    out[2] = d[0]*e2[0] + d[1]*e2[1] + d[2]*e2[2];
}

void AsmAxis_mean_dir(const AsmAxis *a, double z0, double z1, double dir[3])
{
    double sum[3] = { 0.0, 0.0, 0.0 };
    size_t n = 0;
    for (int pass = 0; pass < 2 && n == 0; pass++) {
        for (size_t i = 0; i + 1 < a->n; i++) {
            double zm = 0.5 * (a->p[i*3] + a->p[(i+1)*3]);
            if (pass == 0 && (zm < z0 || zm > z1)) continue;
            double e[3] = { a->p[(i+1)*3]-a->p[i*3], a->p[(i+1)*3+1]-a->p[i*3+1], a->p[(i+1)*3+2]-a->p[i*3+2] };
            double L = sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
            if (L < 1e-9) continue;
            for (int k = 0; k < 3; k++) sum[k] += e[k] / L;
            n++;
        }
    }
    double L = sqrt(sum[0]*sum[0] + sum[1]*sum[1] + sum[2]*sum[2]);
    if (n == 0 || L < 1e-9) { dir[0] = 1.0; dir[1] = 0.0; dir[2] = 0.0; return; }
    for (int k = 0; k < 3; k++) dir[k] = sum[k] / L;
}

int AsmAxis_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    /* a straight tilted line: s is the distance along it, r the distance from it, the frame is rigid */
    {
        double pt[3] = { 100.0, 200.0, 300.0 }, dir[3] = { cos(0.5), sin(0.5), 0.0 };
        AsmAxis *a = AsmAxis_line(arena, pt, dir, 0.0, 1000.0);
        double q[3] = { pt[0] + 400.0 * dir[0] + 30.0 * (-dir[1]), pt[1] + 400.0 * dir[1] + 30.0 * dir[0], pt[2] + 40.0 };
        double s, r, d[3], f[3];
        AsmAxis_project(a, q, &s, &r, d, f);
        double s_expect = 400.0 - (0.0 - pt[0]) / dir[0] * 0.0;   /* s is measured from the line's z0 row */
        (void)s_expect;
        if (fabs(r - 50.0) > 1e-6 || fabs(d[0] - dir[0]) > 1e-9 || fabs(d[1] - dir[1]) > 1e-9) { fprintf(stderr, "  asm_axis selftest FAIL: line projection r %.4f dir (%.4f, %.4f)\n", r, d[0], d[1]); fails++; }
        double o1[3], o2[3];
        double q2[3] = { q[0] + 10.0 * dir[0], q[1] + 10.0 * dir[1], q[2] };
        AsmAxis_frame(a, q, o1); AsmAxis_frame(a, q2, o2);
        if (fabs((o2[0] - o1[0]) - 10.0) > 1e-6 || fabs(o2[1] - o1[1]) > 1e-6 || fabs(o2[2] - o1[2]) > 1e-6 || fabs(hypot(o1[1], o1[2]) - 50.0) > 1e-6) {
            fprintf(stderr, "  asm_axis selftest FAIL: line frame (%.3f %.3f %.3f) -> (%.3f %.3f %.3f)\n", o1[0], o1[1], o1[2], o2[0], o2[1], o2[2]); fails++;
        }
    }
    /* a bent polyline: the foot moves to the right segment, arc length accumulates, extrapolation past the end */
    {
        AsmAxis *a = ARENA_ALLOC(arena, sizeof(AsmAxis));
        a->n = 3; a->p = ARENA_ALLOC(arena, 9 * sizeof(double)); a->s = ARENA_ALLOC(arena, 3 * sizeof(double)); a->t = ARENA_ALLOC(arena, 9 * sizeof(double)); a->is_line = 0; a->n_rows = 3;
        a->fit_rms = 0.0; a->r_curv_min = 1e300;
        double rows[9] = { 0, 0, 0,  100, 0, 0,  200, 100, 0 };
        memcpy(a->p, rows, sizeof rows);
        aa_finish(a);
        double s, r, d[3];
        double q[3] = { 50.0, 20.0, 0.0 };
        AsmAxis_project(a, q, &s, &r, d, NULL);
        /* the row tangents are continuous (row 1 = the mean of both segments), so the tangent at the
         * middle of the first segment leans toward the bend: only its sense is checked */
        if (fabs(s - 50.0) > 1e-9 || fabs(r - 20.0) > 1e-9 || d[0] < 0.9) { fprintf(stderr, "  asm_axis selftest FAIL: polyline first segment s %.3f r %.3f dir %.3f\n", s, r, d[0]); fails++; }
        double q3[3] = { 150.0, 50.0, 0.0 };   /* on the second segment's line, at its middle */
        AsmAxis_project(a, q3, &s, &r, d, NULL);
        if (fabs(s - (100.0 + 0.5 * sqrt(2.0) * 100.0)) > 1e-6 || r > 1e-9 || d[0] < 0.7 || d[0] > 0.95) { fprintf(stderr, "  asm_axis selftest FAIL: polyline second segment s %.3f r %.3f dir %.3f\n", s, r, d[0]); fails++; }
        double q4[3] = { 300.0, 200.0, 0.0 };   /* beyond the end, on the extended line */
        AsmAxis_project(a, q4, &s, &r, d, NULL);
        if (fabs(s - (100.0 + 2.0 * sqrt(2.0) * 100.0)) > 1e-6 || r > 1e-9) { fprintf(stderr, "  asm_axis selftest FAIL: extrapolation s %.3f r %.3f\n", s, r); fails++; }   /* t = 2 on the last segment */
        double dm[3];
        AsmAxis_mean_dir(a, 0.0, 100.0, dm);
        if (fabs(dm[0] - 1.0) > 1e-9) { fprintf(stderr, "  asm_axis selftest FAIL: mean dir over the first segment (%.3f)\n", dm[0]); fails++; }
    }
    Arena_dispose(&arena);
    if (fails == 0) fprintf(stderr, "  asm_axis selftest: all passed\n");
    return fails;
}
