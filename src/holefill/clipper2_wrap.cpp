/*
 * clipper2_wrap.cpp — C++ wrapper calling Clipper2 Union (EvenOdd).
 *
 * Scales double coords by SCALE → int64, calls BooleanOp64,
 * scales result back to double.  All output is malloc'd.
 */
#include "clipper2_wrap.h"
#include "clipper2/clipper.h"

#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <vector>

using namespace Clipper2Lib;

static const double SCALE = 10000.0;

extern "C"
int Clipper2_union(const double *pts_xy, size_t n,
                   double **out_pts, int **out_counts,
                   int *out_n_polys)
{
    if (!pts_xy || n < 3 || !out_pts || !out_counts || !out_n_polys) {
        return -1;
    }

    *out_pts     = nullptr;
    *out_counts  = nullptr;
    *out_n_polys = 0;

    /* Build int64 path from input doubles */
    Path64 subject;
    subject.reserve(n);
    for (size_t i = 0; i < n; i++) {
        int64_t x = static_cast<int64_t>(std::round(pts_xy[i * 2 + 0] * SCALE));
        int64_t y = static_cast<int64_t>(std::round(pts_xy[i * 2 + 1] * SCALE));
        subject.emplace_back(x, y);
    }

    Paths64 subjects;
    subjects.push_back(subject);

    /* Union with EvenOdd fill — resolves self-intersections */
    Paths64 solution;
    solution = Union(subjects, FillRule::EvenOdd);

    if (solution.empty()) {
        return 0; /* valid but empty result */
    }

    /* Count total vertices and polygons */
    int n_polys = static_cast<int>(solution.size());
    size_t total_verts = 0;
    for (const auto &path : solution) {
        total_verts += path.size();
    }

    if (total_verts == 0) {
        return 0;
    }

    /* Allocate output arrays */
    double *pts = static_cast<double *>(malloc(total_verts * 2 * sizeof(double)));
    int *counts = static_cast<int *>(malloc(static_cast<size_t>(n_polys) * sizeof(int)));
    if (!pts || !counts) {
        free(pts);
        free(counts);
        return -1;
    }

    /* Fill output */
    size_t vi = 0;
    for (int p = 0; p < n_polys; p++) {
        const Path64 &path = solution[static_cast<size_t>(p)];
        counts[p] = static_cast<int>(path.size());
        for (const auto &pt : path) {
            pts[vi * 2 + 0] = static_cast<double>(pt.x) / SCALE;
            pts[vi * 2 + 1] = static_cast<double>(pt.y) / SCALE;
            vi++;
        }
    }

    *out_pts     = pts;
    *out_counts  = counts;
    *out_n_polys = n_polys;
    return 0;
}
