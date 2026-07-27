#define _USE_MATH_DEFINES
#include "eig3.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Cyclic Jacobi for a symmetric 3x3 matrix. Input `a` is copied (not modified).
 * eigenvalues ascending; eigenvectors as columns. This is pca.c::jacobi_3x3 with
 * a const input so it can be shared without exposing pca.c's normal-sign rule. */
static void jacobi_3x3(const double a[3][3], double eigenvalues[3],
                       double eigenvectors[3][3])
{
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            eigenvectors[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    double m[3][3];
    memcpy(m, a, sizeof(m));

    for (int iter = 0; iter < 100; iter++) {
        /* Largest off-diagonal magnitude. */
        double max_off = 0.0;
        int p = 0, q = 1;
        for (int i = 0; i < 3; i++) {
            for (int j = i + 1; j < 3; j++) {
                double absv = fabs(m[i][j]);
                if (absv > max_off) { max_off = absv; p = i; q = j; }
            }
        }
        if (max_off < 1e-15) break;

        double app = m[p][p], aqq = m[q][q], apq = m[p][q];
        double diff = aqq - app;
        double theta = (fabs(diff) < 1e-30) ? M_PI / 4.0
                                            : 0.5 * atan2(2.0 * apq, diff);
        double c = cos(theta), s = sin(theta);

        /* M' = G^T M G, computed into a fresh copy to avoid aliasing. */
        double nm[3][3];
        memcpy(nm, m, sizeof(nm));
        nm[p][p] = c * c * m[p][p] - 2.0 * s * c * m[p][q] + s * s * m[q][q];
        nm[q][q] = s * s * m[p][p] + 2.0 * s * c * m[p][q] + c * c * m[q][q];
        nm[p][q] = 0.0;
        nm[q][p] = 0.0;
        for (int r = 0; r < 3; r++) {
            if (r != p && r != q) {
                nm[r][p] = c * m[r][p] - s * m[r][q];
                nm[p][r] = nm[r][p];
                nm[r][q] = s * m[r][p] + c * m[r][q];
                nm[q][r] = nm[r][q];
            }
        }
        memcpy(m, nm, sizeof(m));

        for (int r = 0; r < 3; r++) {
            double vp = eigenvectors[r][p];
            double vq = eigenvectors[r][q];
            eigenvectors[r][p] = c * vp - s * vq;
            eigenvectors[r][q] = s * vp + c * vq;
        }
    }

    /* Sort ascending by eigenvalue (insertion sort over an index permutation). */
    double evals[3] = { m[0][0], m[1][1], m[2][2] };
    int order[3] = { 0, 1, 2 };
    for (int i = 1; i < 3; i++) {
        int k = i;
        while (k > 0 && evals[order[k - 1]] > evals[order[k]]) {
            int t = order[k - 1];
            order[k - 1] = order[k];
            order[k] = t;
            k--;
        }
    }
    double sorted_evecs[3][3];
    for (int i = 0; i < 3; i++) {
        eigenvalues[i] = evals[order[i]];
        for (int j = 0; j < 3; j++) {
            sorted_evecs[j][i] = eigenvectors[j][order[i]];
        }
    }
    memcpy(eigenvectors, sorted_evecs, sizeof(sorted_evecs));
}

void Eig3_sym(const double a[3][3], double evals[3], double evecs[3][3])
{
    jacobi_3x3(a, evals, evecs);
}

void Eig3_smallest(const double a[3][3], double *out_eval, double out_vec[3])
{
    double evals[3], evecs[3][3];
    jacobi_3x3(a, evals, evecs);
    if (out_eval) *out_eval = evals[0];

    double v[3] = { evecs[0][0], evecs[1][0], evecs[2][0] };
    double len = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (len > 1e-300) { v[0] /= len; v[1] /= len; v[2] /= len; }

    /* Deterministic sign: first component with magnitude > 1e-12 made positive. */
    for (int i = 0; i < 3; i++) {
        if (fabs(v[i]) > 1e-12) {
            if (v[i] < 0.0) { v[0] = -v[0]; v[1] = -v[1]; v[2] = -v[2]; }
            break;
        }
    }
    out_vec[0] = v[0];
    out_vec[1] = v[1];
    out_vec[2] = v[2];
}

/* ============================================================================
 * Self-test.
 * ==========================================================================*/

static double resid_norm(const double a[3][3], double lambda, const double v[3])
{
    double r0 = a[0][0] * v[0] + a[0][1] * v[1] + a[0][2] * v[2] - lambda * v[0];
    double r1 = a[1][0] * v[0] + a[1][1] * v[1] + a[1][2] * v[2] - lambda * v[1];
    double r2 = a[2][0] * v[0] + a[2][1] * v[1] + a[2][2] * v[2] - lambda * v[2];
    return sqrt(r0 * r0 + r1 * r1 + r2 * r2);
}

int Eig3_selftest(void)
{
    int fails = 0;

    /* (1) Diagonal: eigenvalues are the (sorted) diagonal. */
    {
        double a[3][3] = { { 5, 0, 0 }, { 0, 1, 0 }, { 0, 0, 3 } };
        double ev[3], evec[3][3];
        Eig3_sym(a, ev, evec);
        int ok = (fabs(ev[0] - 1.0) < 1e-9) && (fabs(ev[1] - 3.0) < 1e-9)
              && (fabs(ev[2] - 5.0) < 1e-9);
        if (!ok) { fprintf(stderr, "[eig3] DIAG FAIL ev=%.4f %.4f %.4f\n",
                           ev[0], ev[1], ev[2]); fails++; }
    }

    /* (2) General symmetric: residual ||A v - lambda v|| tiny for every pair,
     *     eigenvectors unit-length, eigenvalues ascending. */
    {
        double a[3][3] = { { 4, 1, 2 }, { 1, 3, 0 }, { 2, 0, 5 } };
        double ev[3], evec[3][3];
        Eig3_sym(a, ev, evec);
        int ok = (ev[0] <= ev[1] + 1e-12) && (ev[1] <= ev[2] + 1e-12);
        for (int c = 0; c < 3; c++) {
            double v[3] = { evec[0][c], evec[1][c], evec[2][c] };
            double n = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            ok = ok && (fabs(n - 1.0) < 1e-9);
            ok = ok && (resid_norm(a, ev[c], v) < 1e-9);
        }
        if (!ok) { fprintf(stderr, "[eig3] GENERAL FAIL ev=%.6f %.6f %.6f\n",
                           ev[0], ev[1], ev[2]); fails++; }
    }

    /* (3) Repeated eigenvalue + smallest-eigenpair API: no NaN, unit vector,
     *     deterministic sign across two calls. */
    {
        double a[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 5 } };
        double l1, v1[3], l2, v2[3];
        Eig3_smallest(a, &l1, v1);
        Eig3_smallest(a, &l2, v2);
        int ok = (fabs(l1 - 1.0) < 1e-9);
        double n = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
        ok = ok && (fabs(n - 1.0) < 1e-9);
        ok = ok && !(v1[0] != v1[0] || v1[1] != v1[1] || v1[2] != v1[2]);
        ok = ok && (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2]);
        if (!ok) { fprintf(stderr, "[eig3] REPEATED FAIL l=%.4f v=%.4f %.4f %.4f\n",
                           l1, v1[0], v1[1], v1[2]); fails++; }
    }

    /* (4) Zero matrix: all eigenvalues 0, no NaN, vectors unit. */
    {
        double a[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
        double ev[3], evec[3][3];
        Eig3_sym(a, ev, evec);
        int ok = (fabs(ev[0]) < 1e-12) && (fabs(ev[2]) < 1e-12);
        for (int c = 0; c < 3; c++) {
            double v[3] = { evec[0][c], evec[1][c], evec[2][c] };
            double nn = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            ok = ok && (fabs(nn - 1.0) < 1e-9);
        }
        if (!ok) { fprintf(stderr, "[eig3] ZERO FAIL\n"); fails++; }
    }

    if (fails == 0) fprintf(stderr, "[eig3 selftest] ok\n");
    return fails;
}
