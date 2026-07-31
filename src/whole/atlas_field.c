#include "atlas_field.h"

#include <limits.h>
#include <math.h>
#include <string.h>

static double dot3_field(const double a[3], const double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static int triangle_gradients(const double *p0, const double *p1,
                              const double *p2, double gradient[3][2],
                              double *area)
{
    double e01[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
    double e02[3] = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
    double l01 = sqrt(dot3_field(e01, e01));
    if (!isfinite(l01) || l01 <= 1e-14) return -1;
    double axis0[3] = {e01[0] / l01, e01[1] / l01, e01[2] / l01};
    double x2 = dot3_field(e02, axis0);
    double orth[3] = {e02[0] - x2 * axis0[0],
                      e02[1] - x2 * axis0[1],
                      e02[2] - x2 * axis0[2]};
    double y2 = sqrt(dot3_field(orth, orth));
    if (!isfinite(y2) || y2 <= 1e-14) return -1;

    double det = l01 * y2;
    *area = 0.5 * det;
    gradient[0][0] = -1.0 / l01;
    gradient[0][1] = (x2 - l01) / det;
    gradient[1][0] = 1.0 / l01;
    gradient[1][1] = -x2 / det;
    gradient[2][0] = 0.0;
    gradient[2][1] = 1.0 / y2;
    return 0;
}

static int validate_field(const AtlasFieldProblem *p)
{
    size_t nvar = 0;
    if (p == NULL || p->position == NULL || p->u0 == NULL ||
        p->nvertex == 0 || p->nvertex > (size_t)INT_MAX ||
        p->nauxiliary > (size_t)INT_MAX - p->nvertex ||
        p->triangle == NULL || p->ntriangle == 0 ||
        p->ntriangle > (size_t)INT_MAX ||
        (p->nobservation > 0 && p->observation == NULL) ||
        (p->nconstraint > 0 &&
         (p->constraint == NULL || p->constraint_coeff == NULL)) ||
        (p->nanchor > 0 && p->anchor == NULL)) return -1;
    nvar = p->nvertex + p->nauxiliary;

    for (size_t i = 0; i < p->nvertex; i++) {
        if (!isfinite(p->u0[i])) return -1;
        for (int k = 0; k < 3; k++)
            if (!isfinite(p->position[i * 3 + (size_t)k])) return -1;
    }
    for (size_t t = 0; t < p->ntriangle; t++) {
        if (p->triangle_weight != NULL &&
            (!isfinite(p->triangle_weight[t]) ||
             p->triangle_weight[t] < 0.0)) return -1;
        for (int k = 0; k < 3; k++) {
            int32_t v = p->triangle[t * 3 + (size_t)k];
            if (v < 0 || (size_t)v >= p->nvertex) return -1;
        }
    }
    for (size_t o = 0; o < p->nobservation; o++) {
        const AtlasFieldObservation *obs = &p->observation[o];
        if (!isfinite(obs->target) || !isfinite(obs->weight) ||
            obs->weight < 0.0) return -1;
        double sum = 0.0;
        for (int k = 0; k < 3; k++) {
            if (obs->vertex[k] < 0 ||
                (size_t)obs->vertex[k] >= p->nvertex ||
                !isfinite(obs->bary[k]) || obs->bary[k] < -1e-8 ||
                obs->bary[k] > 1.0 + 1e-8) return -1;
            sum += obs->bary[k];
        }
        if (fabs(sum - 1.0) > 1e-8) return -1;
    }
    for (size_t r = 0; r < p->nconstraint; r++) {
        const AtlasFieldConstraint *row = &p->constraint[r];
        if (row->count <= 0 || row->first > p->nconstraint_coeff ||
            (size_t)row->count > p->nconstraint_coeff - row->first ||
            !isfinite(row->target) || !isfinite(row->weight) ||
            row->weight < 0.0)
            return -1;
        for (int32_t j = 0; j < row->count; j++) {
            const AtlasFieldConstraintCoeff *c =
                &p->constraint_coeff[row->first + (size_t)j];
            if (c->variable < 0 || (size_t)c->variable >= nvar ||
                !isfinite(c->coefficient))
                return -1;
        }
    }
    for (size_t a = 0; a < p->nanchor; a++)
        if (p->anchor[a].var < 0 ||
            (size_t)p->anchor[a].var >= nvar ||
            !isfinite(p->anchor[a].value)) return -1;
    return 0;
}

int AtlasField_build(Arena_T arena, const AtlasFieldProblem *p,
                     AtlasFieldSystem *system)
{
    if (arena == NULL || system == NULL || validate_field(p) != 0) return -1;

    size_t valid_triangle = 0, skipped = 0;
    for (size_t t = 0; t < p->ntriangle; t++) {
        int32_t a = p->triangle[t * 3 + 0];
        int32_t b = p->triangle[t * 3 + 1];
        int32_t c = p->triangle[t * 3 + 2];
        double gradient[3][2], area = 0.0;
        if (triangle_gradients(&p->position[(size_t)a * 3],
                               &p->position[(size_t)b * 3],
                               &p->position[(size_t)c * 3],
                               gradient, &area) == 0 && area > 0.0)
            valid_triangle++;
        else
            skipped++;
    }
    if (valid_triangle == 0) return -1;

    if (p->nobservation > SIZE_MAX / 3) return -1;
    size_t observation_coeff = 3 * p->nobservation;
    size_t constraint_coeff = 0;
    for (size_t r = 0; r < p->nconstraint; r++) {
        size_t count = (size_t)p->constraint[r].count;
        if (constraint_coeff > SIZE_MAX - count) return -1;
        constraint_coeff += count;
    }
    if (p->nobservation > SIZE_MAX - p->nconstraint ||
        valid_triangle >
            (SIZE_MAX - p->nobservation - p->nconstraint) / 2 ||
        observation_coeff > SIZE_MAX - constraint_coeff ||
        valid_triangle >
            (SIZE_MAX - observation_coeff - constraint_coeff) / 6)
        return -1;
    size_t nrows = 2 * valid_triangle + p->nobservation + p->nconstraint;
    size_t ncoeff = 6 * valid_triangle + observation_coeff + constraint_coeff;
    if (nrows > SIZE_MAX / sizeof(MonotoneQpRow) ||
        ncoeff > SIZE_MAX / sizeof(MonotoneQpCoeff)) return -1;
    MonotoneQpRow *rows = (MonotoneQpRow *)ARENA_ALLOC(
        arena, nrows * sizeof(MonotoneQpRow));
    MonotoneQpCoeff *coeff = (MonotoneQpCoeff *)ARENA_ALLOC(
        arena, ncoeff * sizeof(MonotoneQpCoeff));

    size_t nr = 0, nk = 0;
    for (size_t t = 0; t < p->ntriangle; t++) {
        int32_t vertex[3] = {
            p->triangle[t * 3 + 0],
            p->triangle[t * 3 + 1],
            p->triangle[t * 3 + 2]
        };
        double gradient[3][2], area = 0.0;
        if (triangle_gradients(&p->position[(size_t)vertex[0] * 3],
                               &p->position[(size_t)vertex[1] * 3],
                               &p->position[(size_t)vertex[2] * 3],
                               gradient, &area) != 0 || area <= 0.0)
            continue;
        double weight = area * (p->triangle_weight != NULL
                              ? p->triangle_weight[t] : 1.0);
        for (int axis = 0; axis < 2; axis++) {
            MonotoneQpRow *row = &rows[nr++];
            row->first = nk;
            row->count = 3;
            row->target = 0.0;
            row->weight = weight;
            row->kind = ATLAS_FIELD_ROW_GRADIENT;
            row->owner = (int32_t)t;
            for (int k = 0; k < 3; k++) {
                coeff[nk].var = vertex[k];
                coeff[nk].value = gradient[k][axis];
                row->target += gradient[k][axis] * p->u0[vertex[k]];
                nk++;
            }
        }
    }
    for (size_t o = 0; o < p->nobservation; o++) {
        const AtlasFieldObservation *obs = &p->observation[o];
        MonotoneQpRow *row = &rows[nr++];
        row->first = nk;
        row->count = 3;
        row->target = obs->target;
        row->weight = obs->weight;
        row->kind = ATLAS_FIELD_ROW_OBSERVATION;
        row->owner = obs->source;
        for (int k = 0; k < 3; k++) {
            coeff[nk].var = obs->vertex[k];
            coeff[nk].value = obs->bary[k];
            nk++;
        }
    }
    for (size_t r = 0; r < p->nconstraint; r++) {
        const AtlasFieldConstraint *src = &p->constraint[r];
        MonotoneQpRow *row = &rows[nr++];
        row->first = nk;
        row->count = src->count;
        row->target = src->target;
        row->weight = src->weight;
        row->kind = src->kind;
        row->owner = src->source;
        for (int32_t j = 0; j < src->count; j++) {
            const AtlasFieldConstraintCoeff *c =
                &p->constraint_coeff[src->first + (size_t)j];
            coeff[nk].var = c->variable;
            coeff[nk].value = c->coefficient;
            nk++;
        }
    }
    if (nr != nrows || nk != ncoeff) return -1;

    memset(system, 0, sizeof *system);
    system->gradient_rows = 2 * valid_triangle;
    system->observation_rows = p->nobservation;
    system->constraint_rows = p->nconstraint;
    system->skipped_degenerate_triangles = skipped;
    system->qp.nvar = p->nvertex + p->nauxiliary;
    system->qp.rows = rows;
    system->qp.nrows = nrows;
    system->qp.coeff = coeff;
    system->qp.ncoeff = ncoeff;
    system->qp.anchors = p->anchor;
    system->qp.nanchors = p->nanchor;
    return 0;
}
