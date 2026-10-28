#include "fiber_rigid.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const double fr_period = 1.57079632679489661923;
static const double fr_degrees = 57.2957795130823208768;

typedef struct FrSum { double value, correction; } FrSum;
typedef struct FrEvent { double mean, weight; } FrEvent;

static void fr_add(FrSum *sum, double value)
{
    double y = value - sum->correction, t = sum->value + y;
    sum->correction = (t - sum->value) - y;
    sum->value = t;
}

static double fr_wrap(double value, double period)
{ return value - period * floor((value + .5 * period) / period); }

static int fr_compare(const void *a, const void *b)
{
    const FrEvent *x = a, *y = b;
    /* Descending means are ascending branch-crossing angles. */
    if (x->mean != y->mean) return x->mean > y->mean ? -1 : 1;
    return x->weight < y->weight ? -1 : x->weight > y->weight;
}

static double fr_dot(const double a[3], const double b[3])
{ return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }

static int fr_project(const float *xyz, const int32_t *face, const double *uv,
                      const SlimRefineDirection *direction, double *area,
                      double *mean, double *delta, double *score)
{
    size_t a = (size_t)face[0], b = (size_t)face[1], c = (size_t)face[2];
    double x[3], y[3], other[3], length, along, height, j[4];
    for (size_t k = 0; k < 3; k++) {
        x[k] = (double)xyz[3*b+k] - xyz[3*a+k];
        other[k] = (double)xyz[3*c+k] - xyz[3*a+k];
    }
    length = sqrt(fr_dot(x, x));
    if (!(length > 0) || !isfinite(length)) return -1;
    for (size_t k = 0; k < 3; k++) x[k] /= length;
    along = fr_dot(x, other);
    for (size_t k = 0; k < 3; k++) y[k] = other[k] - along*x[k];
    height = sqrt(fr_dot(y, y));
    if (!(height > 0) || !isfinite(height)) return -1;
    *area = .5 * length * height;
    if (!(direction->weight > 0)) return 0;
    for (size_t k = 0; k < 3; k++) y[k] /= height;
    /* Differences first: do not cancel large packed-canvas coordinates. */
    for (size_t k = 0; k < 2; k++) {
        j[2*k] = (uv[2*b+k] - uv[2*a+k]) / length;
        j[2*k+1] = (uv[2*c+k] - uv[2*a+k] - along*j[2*k]) / height;
    }
    double hx = fr_dot(direction->h, x), hy = fr_dot(direction->h, y);
    double vx = fr_dot(direction->v, x), vy = fr_dot(direction->v, y);
    double hu = j[0]*hx + j[1]*hy, hv = j[2]*hx + j[3]*hy;
    double vu = j[0]*vx + j[1]*vy, vv = j[2]*vx + j[3]*vy;
    if (!(hypot(hu, hv) > 0) || !(hypot(vu, vv) > 0) ||
            !isfinite(hu+hv+vu+vv)) return -1;
    double ah = atan2(hv, hu), av = atan2(vv, vu);
    *mean = fr_wrap(.5*(ah + av - fr_period), fr_period);
    *delta = fr_wrap(ah - av + fr_period, 2*fr_period);
    /* Independent complementary-axis score, also valid for reflections/shear. */
    double h = atan2(fabs(hv), fabs(hu)), v = atan2(fabs(vu), fabs(vv));
    *score = .5*fmin(h*h + v*v,
                    (fr_period-h)*(fr_period-h) + (fr_period-v)*(fr_period-v));
    return 0;
}

static void fr_interval(double mass, double linear, double square,
                         double lo, double hi, double *best, double *angle)
{
    double t = fmax(lo, fmin(hi, -linear/mass));
    double value = square + t*(2*linear + t*mass);
    if (value < *best) { *best = value; *angle = t; }
}

int FiberRigid_fit(Arena_T arena, const float *xyz, size_t nv,
                   const int32_t *faces, size_t nf, const double *uv,
                   const SlimRefineDirection *directions, FiberRigidStats *stats)
{
    if (stats) memset(stats, 0, sizeof *stats);
    if (!arena || !xyz || !faces || !uv || !directions || !stats || !nv || !nf ||
            nv > SIZE_MAX/(3*sizeof(float)) || nv > SIZE_MAX/(2*sizeof(double)) ||
            nf > SIZE_MAX/(3*sizeof(int32_t)) || nf > SIZE_MAX/sizeof(FrEvent)) return -1;
    for (size_t k = 0; k < 3*nv; k++) if (!isfinite(xyz[k])) return -1;
    for (size_t k = 0; k < 2*nv; k++) if (!isfinite(uv[k])) return -1;
    for (size_t k = 0; k < 3*nf; k++)
        if (faces[k] < 0 || (size_t)faces[k] >= nv) return -1;
    Arena_Mark mark = Arena_save(arena);
    FrEvent *events = ARENA_ALLOC(arena, nf*sizeof *events);
    FrSum mass = {0}, linear = {0}, square = {0}, constant = {0}, direct = {0};
    FrSum area = {0}, supported = {0}, centroid[2] = {{0}, {0}};
    size_t count = 0;
    int rc = -1;
    for (size_t f = 0; f < nf; f++) {
        const SlimRefineDirection *d = directions+f;
        double fa = 0, mean = 0, delta = 0, score = 0;
        if (!isfinite(d->weight) || d->weight < 0) goto done;
        if (d->weight > 0) for (size_t k = 0; k < 3; k++)
            if (!isfinite(d->h[k]) || !isfinite(d->v[k])) goto done;
        if (fr_project(xyz, faces+3*f, uv, d, &fa, &mean, &delta, &score)) goto done;
        fr_add(&area, fa);
        for (size_t k = 0; k < 2; k++) {
            double center = 0;
            for (size_t q = 0; q < 3; q++) center += (uv[2*(size_t)faces[3*f+q]+k] - uv[k])/3;
            fr_add(centroid+k, fa*center);
        }
        if (!(d->weight > 0)) continue;
        double weight = fa*d->weight;
        if (!(weight > 0) || !isfinite(weight)) goto done;
        events[count].mean = mean; events[count].weight = weight; count++;
        fr_add(&supported, fa); fr_add(&mass, weight);
        fr_add(&linear, weight*mean); fr_add(&square, weight*mean*mean);
        fr_add(&constant, .25*weight*delta*delta); fr_add(&direct, weight*score);
    }
    if (!(mass.value > 0) || !isfinite(mass.value+area.value+square.value+constant.value)) goto done;
    if (fabs(direct.value-square.value-constant.value) > 1e-10*mass.value) goto done;
    stats->faces = nf; stats->supported_faces = count;
    stats->area = area.value; stats->supported_area = supported.value;
    stats->evidence_mass = mass.value;
    for (size_t k = 0; k < 2; k++) stats->center[k] = uv[k] + centroid[k].value/area.value;
    stats->rms_before_degrees = sqrt(direct.value/mass.value)*fr_degrees;
    stats->signed_mean_before_degrees = linear.value/mass.value*fr_degrees;
    fprintf(stderr, "[fiber rigid] %zu/%zu supported faces; RMS %.9g degrees; sorting scalar events\n",
            count, nf, stats->rms_before_degrees);
    qsort(events, count, sizeof *events, fr_compare);
    double best = square.value, angle = 0, lo = 0;
    for (size_t e = 0; e < count; e++) {
        const FrEvent *v = events+e;
        double hi = .5*fr_period-v->mean;
        fr_interval(mass.value, linear.value, square.value, lo, hi, &best, &angle);
        fr_add(&linear, -v->weight*fr_period);
        fr_add(&square, v->weight*(fr_period*fr_period-2*fr_period*v->mean));
        lo = hi;
    }
    fr_interval(mass.value, linear.value, square.value, lo, fr_period, &best, &angle);
    stats->radians = fr_wrap(angle, fr_period);
    /* Reevaluate the retained angle instead of trusting subtractive sweep loss. */
    FrSum after = {0}, signed_after = {0};
    for (size_t e = 0; e < count; e++) {
        double mean = fr_wrap(events[e].mean+stats->radians, fr_period);
        fr_add(&after, events[e].weight*mean*mean);
        fr_add(&signed_after, events[e].weight*mean);
    }
    stats->rms_after_degrees = sqrt((after.value+constant.value)/mass.value)*fr_degrees;
    stats->signed_mean_after_degrees = signed_after.value/mass.value*fr_degrees;
    if (!isfinite(stats->radians+stats->rms_after_degrees) ||
            stats->rms_after_degrees > stats->rms_before_degrees+1e-9) goto done;
    rc = 0;
done:
    Arena_restore(arena, mark);
    if (rc) memset(stats, 0, sizeof *stats);
    return rc;
}

int FiberRigid_apply(const double *uv, size_t nv, const FiberRigidStats *stats, double *out)
{
    if (!uv || !out || !stats || !nv || nv > SIZE_MAX/(2*sizeof(double)) ||
            !isfinite(stats->radians+stats->center[0]+stats->center[1])) return -1;
    for (size_t k = 0; k < 2*nv; k++) if (!isfinite(uv[k])) return -1;
    double co = cos(stats->radians), si = sin(stats->radians);
    for (size_t v = 0; v < nv; v++) {
        double x = uv[2*v]-stats->center[0], y = uv[2*v+1]-stats->center[1];
        out[2*v] = stats->radians == 0 ? uv[2*v] : stats->center[0]+co*x-si*y;
        out[2*v+1] = stats->radians == 0 ? uv[2*v+1] : stats->center[1]+si*x+co*y;
    }
    return 0;
}

#include "fiber_rigid_test.inc"
