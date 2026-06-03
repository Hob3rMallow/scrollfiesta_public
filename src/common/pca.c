#define _USE_MATH_DEFINES
#include "pca.h"
#include <assert.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* 3x3 symmetric Jacobi eigensolver.
 * Input:  a[3][3] symmetric matrix (double).
 * Output: eigenvalues[3] in ascending order, eigenvectors[3][3] as columns.
 * Uses cyclic Jacobi rotations. */
static void jacobi_3x3(double a[3][3], double eigenvalues[3],
                        double eigenvectors[3][3])
{
    /* Initialize eigenvectors to identity */
    int i = 0, j = 0;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            eigenvectors[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    /* Copy input since we modify it */
    double m[3][3];
    memcpy(m, a, sizeof(m));

    for (int iter = 0; iter < 100; iter++) {
        /* Find largest off-diagonal element */
        double max_off = 0.0;
        int p = 0, q = 1;
        for (i = 0; i < 3; i++) {
            for (j = i + 1; j < 3; j++) {
                double absv = fabs(m[i][j]);
                if (absv > max_off) {
                    max_off = absv;
                    p = i;
                    q = j;
                }
            }
        }

        /* Convergence check */
        if (max_off < 1e-15) {
            break;
        }

        /* Compute rotation angle */
        double app = m[p][p];
        double aqq = m[q][q];
        double apq = m[p][q];

        double theta = 0.0;
        double diff = aqq - app;
        if (fabs(diff) < 1e-30) {
            theta = M_PI / 4.0;
        } else {
            theta = 0.5 * atan2(2.0 * apq, diff);
        }

        double c = cos(theta);
        double s = sin(theta);

        /* Apply rotation to matrix: M' = G^T M G */
        double new_m[3][3];
        memcpy(new_m, m, sizeof(new_m));

        new_m[p][p] = c * c * m[p][p] - 2.0 * s * c * m[p][q] + s * s * m[q][q];
        new_m[q][q] = s * s * m[p][p] + 2.0 * s * c * m[p][q] + c * c * m[q][q];
        new_m[p][q] = 0.0;
        new_m[q][p] = 0.0;

        for (int r = 0; r < 3; r++) {
            if (r != p && r != q) {
                new_m[r][p] = c * m[r][p] - s * m[r][q];
                new_m[p][r] = new_m[r][p];
                new_m[r][q] = s * m[r][p] + c * m[r][q];
                new_m[q][r] = new_m[r][q];
            }
        }

        memcpy(m, new_m, sizeof(m));

        /* Update eigenvectors */
        for (int r = 0; r < 3; r++) {
            double vp = eigenvectors[r][p];
            double vq = eigenvectors[r][q];
            eigenvectors[r][p] = c * vp - s * vq;
            eigenvectors[r][q] = s * vp + c * vq;
        }
    }

    /* Extract eigenvalues and sort ascending */
    double evals[3] = { m[0][0], m[1][1], m[2][2] };
    int order[3] = { 0, 1, 2 };

    /* Simple insertion sort for 3 elements */
    for (i = 1; i < 3; i++) {
        int k = i;
        while (k > 0 && evals[order[k - 1]] > evals[order[k]]) {
            int t = order[k - 1];
            order[k - 1] = order[k];
            order[k] = t;
            k--;
        }
    }

    double sorted_evecs[3][3];
    for (i = 0; i < 3; i++) {
        eigenvalues[i] = evals[order[i]];
        for (j = 0; j < 3; j++) {
            sorted_evecs[j][i] = eigenvectors[j][order[i]];
        }
    }
    memcpy(eigenvectors, sorted_evecs, sizeof(sorted_evecs));
}

int PCA_normal(const float *points, size_t n,
               float out_normal[3], float out_centroid[3])
{
    assert(out_normal);

    if (n < 3) {
        out_normal[0] = 0.0f;
        out_normal[1] = 0.0f;
        out_normal[2] = 1.0f;
        if (out_centroid) {
            out_centroid[0] = 0.0f;
            out_centroid[1] = 0.0f;
            out_centroid[2] = 0.0f;
        }
        return -1;
    }

    assert(points);

    /* Compute centroid in double precision */
    double cx = 0.0, cy = 0.0, cz = 0.0;
    for (size_t i = 0; i < n; i++) {
        cx += (double)points[i * 3 + 0];
        cy += (double)points[i * 3 + 1];
        cz += (double)points[i * 3 + 2];
    }
    double inv_n = 1.0 / (double)n;
    cx *= inv_n;
    cy *= inv_n;
    cz *= inv_n;

    if (out_centroid) {
        out_centroid[0] = (float)cx;
        out_centroid[1] = (float)cy;
        out_centroid[2] = (float)cz;
    }

    /* Compute 3x3 covariance matrix in double precision */
    double cov[3][3];
    memset(cov, 0, sizeof(cov));

    for (size_t i = 0; i < n; i++) {
        double dz = (double)points[i * 3 + 0] - cx;
        double dy = (double)points[i * 3 + 1] - cy;
        double dx = (double)points[i * 3 + 2] - cz;

        cov[0][0] += dz * dz;
        cov[0][1] += dz * dy;
        cov[0][2] += dz * dx;
        cov[1][1] += dy * dy;
        cov[1][2] += dy * dx;
        cov[2][2] += dx * dx;
    }

    /* Divide by N (matching Python: cov = centered.T @ centered / N) */
    for (int i = 0; i < 3; i++) {
        for (int j = i; j < 3; j++) {
            cov[i][j] *= inv_n;
            cov[j][i] = cov[i][j];
        }
    }

    /* Jacobi eigensolver */
    double eigenvalues[3];
    double eigenvectors[3][3];
    jacobi_3x3(cov, eigenvalues, eigenvectors);

    /* Smallest eigenvalue is first (index 0) — its eigenvector is column 0 */
    double normal[3] = {
        eigenvectors[0][0],
        eigenvectors[1][0],
        eigenvectors[2][0]
    };

    /* Normalize */
    double len = sqrt(normal[0] * normal[0] +
                      normal[1] * normal[1] +
                      normal[2] * normal[2]);
    if (len > 1e-15) {
        normal[0] /= len;
        normal[1] /= len;
        normal[2] /= len;
    } else {
        normal[0] = 0.0;
        normal[1] = 0.0;
        normal[2] = 1.0;
    }

    /* Sign convention: orient toward a fixed reference vector (1, 1, 1).
     * This is robust to point-cloud differences across cubes: any two
     * eigenvectors (one per cube view) that agree in axis but differ only
     * in sign will both end up pointing toward the same hemisphere of R^3
     * after this normalization, so backface-cull keeps the same shell side
     * in both cubes. The legacy "largest-magnitude positive" rule was
     * deterministic per point cloud but not across slightly-different
     * point clouds (halo overlap causes minor numerical differences in
     * the covariance matrix that can flip which axis is the largest).
     *
     * Tie-breaks (normal nearly perpendicular to (1,1,1)) fall back to
     * "largest-magnitude positive" for stability when the (1,1,1) dot is
     * within 1e-9 of zero. */
    double ref_dot = normal[0] + normal[1] + normal[2];
    if (ref_dot < -1e-9) {
        normal[0] = -normal[0];
        normal[1] = -normal[1];
        normal[2] = -normal[2];
    } else if (ref_dot < 1e-9) {
        /* Near-degenerate vs (1,1,1) -- fall back to largest-magnitude. */
        int max_idx = 0;
        double max_abs = fabs(normal[0]);
        for (int i = 1; i < 3; i++) {
            double a = fabs(normal[i]);
            if (a > max_abs) {
                max_abs = a;
                max_idx = i;
            }
        }
        if (normal[max_idx] < 0.0) {
            normal[0] = -normal[0];
            normal[1] = -normal[1];
            normal[2] = -normal[2];
        }
    }

    out_normal[0] = (float)normal[0];
    out_normal[1] = (float)normal[1];
    out_normal[2] = (float)normal[2];

    return 0;
}

/* Shared: centroid + ascending eigenvalues + eigenvectors (as columns).
 * Component order follows the input layout. Returns 0, or -1 if n < 3. */
static int pca_eigh(const float *points, size_t n,
                    double centroid[3], double evals[3], double evecs[3][3])
{
    if (n < 3) {
        return -1;
    }
    assert(points);

    double cx = 0.0, cy = 0.0, cz = 0.0;
    for (size_t i = 0; i < n; i++) {
        cx += (double)points[i * 3 + 0];
        cy += (double)points[i * 3 + 1];
        cz += (double)points[i * 3 + 2];
    }
    double inv_n = 1.0 / (double)n;
    centroid[0] = cx * inv_n;
    centroid[1] = cy * inv_n;
    centroid[2] = cz * inv_n;

    double cov[3][3];
    memset(cov, 0, sizeof(cov));
    for (size_t i = 0; i < n; i++) {
        double d0 = (double)points[i * 3 + 0] - centroid[0];
        double d1 = (double)points[i * 3 + 1] - centroid[1];
        double d2 = (double)points[i * 3 + 2] - centroid[2];
        cov[0][0] += d0 * d0;
        cov[0][1] += d0 * d1;
        cov[0][2] += d0 * d2;
        cov[1][1] += d1 * d1;
        cov[1][2] += d1 * d2;
        cov[2][2] += d2 * d2;
    }
    for (int i = 0; i < 3; i++) {
        for (int j = i; j < 3; j++) {
            cov[i][j] *= inv_n;
            cov[j][i] = cov[i][j];
        }
    }

    jacobi_3x3(cov, evals, evecs);
    return 0;
}

/* Extract eigenvector column `col`, normalize, orient toward (1,1,1). */
static void axis_from_col(double evecs[3][3], int col, float out_axis[3])
{
    double a[3] = { evecs[0][col], evecs[1][col], evecs[2][col] };
    double len = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    if (len > 1e-15) {
        a[0] /= len; a[1] /= len; a[2] /= len;
    } else {
        a[0] = 1.0; a[1] = 0.0; a[2] = 0.0;
    }
    if (a[0] + a[1] + a[2] < 0.0) {   /* deterministic sign, matches PCA_normal */
        a[0] = -a[0]; a[1] = -a[1]; a[2] = -a[2];
    }
    out_axis[0] = (float)a[0];
    out_axis[1] = (float)a[1];
    out_axis[2] = (float)a[2];
}

static void axis_degenerate(float out_axis[3], float out_centroid[3])
{
    out_axis[0] = 1.0f; out_axis[1] = 0.0f; out_axis[2] = 0.0f;
    if (out_centroid) {
        out_centroid[0] = 0.0f; out_centroid[1] = 0.0f; out_centroid[2] = 0.0f;
    }
}

int PCA_principal_axis(const float *points, size_t n,
                       float out_axis[3], float out_centroid[3])
{
    assert(out_axis);
    double cen[3], evals[3], evecs[3][3];
    if (pca_eigh(points, n, cen, evals, evecs) != 0) {
        axis_degenerate(out_axis, out_centroid);
        return -1;
    }
    axis_from_col(evecs, 2, out_axis);   /* largest eigenvalue */
    if (out_centroid) {
        out_centroid[0] = (float)cen[0];
        out_centroid[1] = (float)cen[1];
        out_centroid[2] = (float)cen[2];
    }
    return 0;
}

int PCA_scroll_axis(const float *points, size_t n,
                    float out_axis[3], float out_centroid[3])
{
    assert(out_axis);
    double cen[3], evals[3], evecs[3][3];
    if (pca_eigh(points, n, cen, evals, evecs) != 0) {
        axis_degenerate(out_axis, out_centroid);
        return -1;
    }
    /* evals ascending. The umbilicus is the "odd one out": the eigenvalue
     * farthest from the median (evals[1]). For a roughly circular wound scroll
     * the two cross-section axes have similar spread, so the axis is the
     * outlier -- robust whether the region is taller than the disk or wider. */
    int col = (evals[2] - evals[1] >= evals[1] - evals[0]) ? 2 : 0;
    axis_from_col(evecs, col, out_axis);
    if (out_centroid) {
        out_centroid[0] = (float)cen[0];
        out_centroid[1] = (float)cen[1];
        out_centroid[2] = (float)cen[2];
    }
    return 0;
}

void PCA_project(const float *points, size_t n, const float dir[3],
                 float *out_proj, float *out_min, float *out_max)
{
    assert(points || n == 0);
    assert(dir);
    assert(out_proj || n == 0);

    if (n == 0) {
        if (out_min) { *out_min = 0.0f; }
        if (out_max) { *out_max = 0.0f; }
        return;
    }

    float mn = 1e30f;
    float mx = -1e30f;

    for (size_t i = 0; i < n; i++) {
        float dot = points[i * 3 + 0] * dir[0]
                  + points[i * 3 + 1] * dir[1]
                  + points[i * 3 + 2] * dir[2];
        out_proj[i] = dot;
        if (dot < mn) { mn = dot; }
        if (dot > mx) { mx = dot; }
    }

    if (out_min) { *out_min = mn; }
    if (out_max) { *out_max = mx; }
}

void PCA_orthonormal_basis(const float normal[3],
                           float out_u[3], float out_v[3])
{
    assert(normal);
    assert(out_u);
    assert(out_v);

    /* Hughes-Moller method: pick axis least aligned with normal */
    float ax = fabsf(normal[0]);
    float ay = fabsf(normal[1]);
    float az = fabsf(normal[2]);

    float t[3];
    if (ax <= ay && ax <= az) {
        t[0] = 1.0f; t[1] = 0.0f; t[2] = 0.0f;
    } else if (ay <= ax && ay <= az) {
        t[0] = 0.0f; t[1] = 1.0f; t[2] = 0.0f;
    } else {
        t[0] = 0.0f; t[1] = 0.0f; t[2] = 1.0f;
    }

    /* u = normalize(t - (t . normal) * normal) */
    float dot = t[0] * normal[0] + t[1] * normal[1] + t[2] * normal[2];
    out_u[0] = t[0] - dot * normal[0];
    out_u[1] = t[1] - dot * normal[1];
    out_u[2] = t[2] - dot * normal[2];

    float len = sqrtf(out_u[0] * out_u[0] +
                      out_u[1] * out_u[1] +
                      out_u[2] * out_u[2]);
    if (len > 1e-10f) {
        out_u[0] /= len;
        out_u[1] /= len;
        out_u[2] /= len;
    }

    /* v = normal x u */
    out_v[0] = normal[1] * out_u[2] - normal[2] * out_u[1];
    out_v[1] = normal[2] * out_u[0] - normal[0] * out_u[2];
    out_v[2] = normal[0] * out_u[1] - normal[1] * out_u[0];
}
