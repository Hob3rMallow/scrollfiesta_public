/*
 * seam_planes.c -- shared seam-plane detection + classification queries.
 * Bodies moved VERBATIM from seam_weld.c (2026-07); see seam_planes.h.
 */
#include "seam_planes.h"

#include <math.h>

/* Beyond the export-halo overhang, both sides of a TRUE seam carry mesh.
 * A cube's genuine outer boundary only has a ~1-vox export sliver on its
 * outside, so requiring mesh past this margin rejects those. Ties to
 * EXPORT_HALO_VOX (=1) in pipeline_constants.h. Must stay > the per-cube
 * trim inset (BPA_OWNED_TRIM_INSET=1.0) so a cube's own trimmed shell can
 * never satisfy the "mesh on both sides" test at a grid-edge plane. */
#define SEAM_MIN_SIDE 2.0

/* A seam-FACING boundary edge lies IN the seam plane; a perpendicular perimeter
 * edge crosses it. We keep only edges whose endpoints sit at nearly the same
 * plane-axis coordinate (|delta| < this). Feeding perimeter edges to the bridge
 * makes the ball roll ALONG the perimeter instead of across the gap, which
 * sprays fragments and leaves the seam open. */
#define SEAM_PARALLEL_EPS 0.75

size_t SeamPlanes_detect(const float *verts, size_t nv,
                         const uint8_t *used, double cube_size,
                         double band, SeamPlane *planes, size_t max_planes)
{
    double mn[3] = { 1e30, 1e30, 1e30 }, mx[3] = { -1e30, -1e30, -1e30 };
    for (size_t v = 0; v < nv; v++) {
        if (!used[v]) continue;
        for (int a = 0; a < 3; a++) {
            double c = verts[v*3+(size_t)a];
            if (c < mn[a]) mn[a] = c;
            if (c > mx[a]) mx[a] = c;
        }
    }
    size_t np = 0;
    for (int a = 0; a < 3 && np < max_planes; a++) {
        long k0 = (long)ceil(mn[a] / cube_size);
        long k1 = (long)floor(mx[a] / cube_size);
        for (long k = k0; k <= k1 && np < max_planes; k++) {
            double c = (double)k * cube_size;
            int has_lo = 0, has_hi = 0;
            for (size_t v = 0; v < nv; v++) {
                if (!used[v]) continue;
                double cc = verts[v*3+(size_t)a];
                if (cc > c - band && cc < c - SEAM_MIN_SIDE) has_lo = 1;
                if (cc > c + SEAM_MIN_SIDE && cc < c + band) has_hi = 1;
                if (has_lo && has_hi) break;
            }
            if (has_lo && has_hi) {
                planes[np].axis = a;
                planes[np].coord = c;
                np++;
            }
        }
    }
    return np;
}

int SeamPlanes_edge_in(const float *verts, int32_t va, int32_t vb,
                       const SeamPlane *planes, size_t np, double band)
{
    for (size_t p = 0; p < np; p++) {
        int ax = planes[p].axis;
        double ca = (double)verts[(size_t)va*3+(size_t)ax];
        double cb = (double)verts[(size_t)vb*3+(size_t)ax];
        if (fabs(ca - planes[p].coord) < band &&
            fabs(cb - planes[p].coord) < band &&
            fabs(ca - cb) < SEAM_PARALLEL_EPS)
            return 1;
    }
    return 0;
}

int SeamPlanes_face_straddle(const float *verts, int32_t a, int32_t b, int32_t c,
                             const SeamPlane *planes, size_t np, double band)
{
    for (size_t p = 0; p < np; p++) {
        int ax = planes[p].axis; double co = planes[p].coord;
        double ca = (double)verts[(size_t)a*3+(size_t)ax];
        double cb = (double)verts[(size_t)b*3+(size_t)ax];
        double cc = (double)verts[(size_t)c*3+(size_t)ax];
        double mn = ca < cb ? (ca < cc ? ca : cc) : (cb < cc ? cb : cc);
        double mx = ca > cb ? (ca > cc ? ca : cc) : (cb > cc ? cb : cc);
        if (mn < co && mx > co && (co - mn) < band && (mx - co) < band)
            return (int)p;   /* straddles this seam plane, within the band */
    }
    return -1;
}

double SeamPlanes_vert_dist(const float *verts, int32_t v,
                            const SeamPlane *planes, size_t np)
{
    double best = 1e30;
    for (size_t p = 0; p < np; p++) {
        double d = fabs((double)verts[(size_t)v*3+(size_t)planes[p].axis] - planes[p].coord);
        if (d < best) best = d;
    }
    return best;
}
