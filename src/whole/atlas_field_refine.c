#include "atlas_field_refine.h"

#include "../common/union_find.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <string.h>

static double afr_triangle_area(const float *a, const float *b,
                                const float *c)
{
    double x0 = (double)b[0] - (double)a[0];
    double x1 = (double)b[1] - (double)a[1];
    double x2 = (double)b[2] - (double)a[2];
    double y0 = (double)c[0] - (double)a[0];
    double y1 = (double)c[1] - (double)a[1];
    double y2 = (double)c[2] - (double)a[2];
    double cx = x1 * y2 - x2 * y1;
    double cy = x2 * y0 - x0 * y2;
    double cz = x0 * y1 - x1 * y0;
    return 0.5 * sqrt(cx * cx + cy * cy + cz * cz);
}

static double afr_row_value(const AtlasFieldConstraint *row,
                            const AtlasFieldConstraintCoeff *coeff,
                            const double *x)
{
    double value = 0.0;
    for (int32_t j = 0; j < row->count; j++) {
        const AtlasFieldConstraintCoeff *c =
            &coeff[row->first + (size_t)j];
        value += c->coefficient * x[c->variable];
    }
    return value;
}

int AtlasFieldRefine_solve(
    Arena_T arena,
    const float *vertices, size_t nvertices,
    const int32_t *triangles, size_t ntriangles,
    const int32_t *vertex_component, size_t ncomponents,
    const double *base_u,
    const AtlasFieldConstraint *constraint, size_t nconstraint,
    const AtlasFieldConstraintCoeff *constraint_coeff,
    size_t nconstraint_coeff,
    size_t nauxiliary, const double *auxiliary_initial,
    const MonotoneQpOptions *qp_options,
    double *out_u, double *out_auxiliary,
    AtlasFieldRefineStats *stats)
{
    if (arena == NULL || vertices == NULL || nvertices == 0 ||
        nvertices > (size_t)INT32_MAX || triangles == NULL ||
        ntriangles == 0 || vertex_component == NULL || ncomponents == 0 ||
        ncomponents > (size_t)INT32_MAX || base_u == NULL ||
        (nconstraint > 0 &&
         (constraint == NULL || constraint_coeff == NULL)) ||
        nauxiliary > (size_t)INT32_MAX - nvertices ||
        (nauxiliary > 0 && auxiliary_initial == NULL) ||
        out_u == NULL || stats == NULL)
        return -1;
    memset(stats, 0, sizeof *stats);
    stats->mesh_components = ncomponents;
    stats->auxiliary_variables = nauxiliary;
    stats->constraints = nconstraint;
    size_t nvar = nvertices + nauxiliary;

    double *position = (double *)ARENA_ALLOC(
        arena, nvertices * 3 * sizeof(double));
    double *component_area = (double *)ARENA_CALLOC(
        arena, ncomponents, sizeof(double));
    for (size_t i = 0; i < nvertices; i++) {
        if (!isfinite(base_u[i])) return -1;
        int32_t component = vertex_component[i];
        if (component < -1 ||
            (component >= 0 && (size_t)component >= ncomponents))
            return -1;
        for (int d = 0; d < 3; d++) {
            double value = (double)vertices[i * 3 + (size_t)d];
            if (!isfinite(value)) return -1;
            position[i * 3 + (size_t)d] = value;
        }
    }

    UnionFind quotient = UF_new(arena, (int32_t)nvar);
    for (size_t f = 0; f < ntriangles; f++) {
        int32_t a = triangles[f * 3];
        int32_t b = triangles[f * 3 + 1];
        int32_t c = triangles[f * 3 + 2];
        if (a < 0 || b < 0 || c < 0 ||
            (size_t)a >= nvertices || (size_t)b >= nvertices ||
            (size_t)c >= nvertices)
            return -1;
        int32_t component = vertex_component[a];
        if (component < 0 || vertex_component[b] != component ||
            vertex_component[c] != component)
            return -1;
        uf_union(&quotient, a, b);
        uf_union(&quotient, a, c);
        component_area[component] += afr_triangle_area(
            &vertices[(size_t)a * 3], &vertices[(size_t)b * 3],
            &vertices[(size_t)c * 3]);
    }

    for (size_t r = 0; r < nconstraint; r++) {
        const AtlasFieldConstraint *row = &constraint[r];
        if (row->count <= 0 || row->first > nconstraint_coeff ||
            (size_t)row->count > nconstraint_coeff - row->first ||
            !isfinite(row->target) || !isfinite(row->weight) ||
            row->weight < 0.0)
            return -1;
        double sum = 0.0;
        int32_t first = -1;
        for (int32_t j = 0; j < row->count; j++) {
            const AtlasFieldConstraintCoeff *c =
                &constraint_coeff[row->first + (size_t)j];
            if (c->variable < 0 || (size_t)c->variable >= nvar ||
                !isfinite(c->coefficient))
                return -1;
            sum += c->coefficient;
            if (c->coefficient == 0.0) continue;
            if (first < 0) first = c->variable;
            else if (row->weight > 0.0)
                uf_union(&quotient, first, c->variable);
        }
        if (first < 0 || fabs(sum) > 1e-8) return -1;
    }

    int32_t *best = (int32_t *)ARENA_ALLOC(
        arena, nvar * sizeof(int32_t));
    double *best_area = (double *)ARENA_ALLOC(
        arena, nvar * sizeof(double));
    for (size_t i = 0; i < nvar; i++) {
        best[i] = -1;
        best_area[i] = -DBL_MAX;
    }
    for (size_t i = 0; i < nvar; i++) {
        int32_t root = uf_find(&quotient, (int32_t)i);
        double area = -1.0;
        if (i < nvertices && vertex_component[i] >= 0)
            area = component_area[vertex_component[i]];
        if (best[root] < 0 || area > best_area[root]) {
            best[root] = (int32_t)i;
            best_area[root] = area;
        }
    }
    MonotoneQpAnchor *anchor = (MonotoneQpAnchor *)ARENA_ALLOC(
        arena, nvar * sizeof(MonotoneQpAnchor));
    size_t nanchor = 0;
    for (size_t i = 0; i < nvar; i++) {
        if (uf_find(&quotient, (int32_t)i) != (int32_t)i || best[i] < 0)
            continue;
        int32_t variable = best[i];
        anchor[nanchor].var = variable;
        anchor[nanchor].value = (size_t)variable < nvertices
            ? base_u[variable]
            : auxiliary_initial[(size_t)variable - nvertices];
        anchor[nanchor].component = (int32_t)nanchor;
        nanchor++;
    }
    stats->quotient_components = nanchor;
    stats->anchors = nanchor;

    AtlasFieldProblem problem;
    memset(&problem, 0, sizeof problem);
    problem.position = position;
    problem.u0 = base_u;
    problem.nvertex = nvertices;
    problem.nauxiliary = nauxiliary;
    problem.triangle = triangles;
    problem.ntriangle = ntriangles;
    problem.constraint = constraint;
    problem.nconstraint = nconstraint;
    problem.constraint_coeff = constraint_coeff;
    problem.nconstraint_coeff = nconstraint_coeff;
    problem.anchor = anchor;
    problem.nanchor = nanchor;
    AtlasFieldSystem system;
    if (AtlasField_build(arena, &problem, &system) != 0) return -1;
    stats->skipped_degenerate_triangles =
        system.skipped_degenerate_triangles;

    double *x = (double *)ARENA_ALLOC(arena, nvar * sizeof(double));
    memcpy(x, base_u, nvertices * sizeof(double));
    for (size_t i = 0; i < nauxiliary; i++)
        x[nvertices + i] = auxiliary_initial[i];
    double *initial = (double *)ARENA_ALLOC(arena, nvar * sizeof(double));
    memcpy(initial, x, nvar * sizeof(double));
    if (MonotoneQp_solve(arena, &system.qp, qp_options, x, NULL,
                         &stats->qp) != 0)
        return -1;

    double correction2 = 0.0;
    for (size_t i = 0; i < nvertices; i++) {
        double d = x[i] - base_u[i];
        correction2 += d * d;
        if (fabs(d) > stats->correction_max)
            stats->correction_max = fabs(d);
        out_u[i] = x[i];
    }
    stats->correction_rms = sqrt(correction2 / (double)nvertices);
    if (out_auxiliary != NULL)
        memcpy(out_auxiliary, &x[nvertices],
               nauxiliary * sizeof(double));

    double edge2 = 0.0;
    size_t nedge = 0;
    for (size_t f = 0; f < ntriangles; f++) {
        for (int e = 0; e < 3; e++) {
            int32_t a = triangles[f * 3 + (size_t)e];
            int32_t b = triangles[f * 3 + (size_t)((e + 1) % 3)];
            double dx = (double)vertices[(size_t)b * 3] -
                        (double)vertices[(size_t)a * 3];
            double dy = (double)vertices[(size_t)b * 3 + 1] -
                        (double)vertices[(size_t)a * 3 + 1];
            double dz = (double)vertices[(size_t)b * 3 + 2] -
                        (double)vertices[(size_t)a * 3 + 2];
            double length = sqrt(dx * dx + dy * dy + dz * dz);
            if (length <= 1e-12) continue;
            double delta = ((x[b] - x[a]) -
                           (base_u[b] - base_u[a])) / length;
            edge2 += delta * delta;
            nedge++;
            if (fabs(delta) > stats->edge_gradient_delta_max)
                stats->edge_gradient_delta_max = fabs(delta);
        }
    }
    stats->edge_gradient_delta_rms =
        nedge ? sqrt(edge2 / (double)nedge) : 0.0;

    double seam_before2 = 0.0, seam_after2 = 0.0;
    double order_before2 = 0.0, order_after2 = 0.0;
    for (size_t r = 0; r < nconstraint; r++) {
        const AtlasFieldConstraint *row = &constraint[r];
        double before = afr_row_value(row, constraint_coeff, initial) -
                        row->target;
        double after = afr_row_value(row, constraint_coeff, x) - row->target;
        if (row->kind == ATLAS_FIELD_ROW_SEAM) {
            seam_before2 += before * before;
            seam_after2 += after * after;
            stats->seam_rows++;
            if (fabs(after) > stats->seam_max_after)
                stats->seam_max_after = fabs(after);
        } else if (row->kind == ATLAS_FIELD_ROW_ORDER) {
            order_before2 += before * before;
            order_after2 += after * after;
            stats->order_rows++;
            if (fabs(after) > stats->order_max_after)
                stats->order_max_after = fabs(after);
        }
    }
    if (stats->seam_rows) {
        stats->seam_rms_before =
            sqrt(seam_before2 / (double)stats->seam_rows);
        stats->seam_rms_after =
            sqrt(seam_after2 / (double)stats->seam_rows);
    }
    if (stats->order_rows) {
        stats->order_rms_before =
            sqrt(order_before2 / (double)stats->order_rows);
        stats->order_rms_after =
            sqrt(order_after2 / (double)stats->order_rows);
    }
    return 0;
}
