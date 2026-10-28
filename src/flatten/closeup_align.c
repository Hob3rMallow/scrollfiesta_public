/* Physical-field proposal, followed by unchanged SLIM and metric selection.
 * This dual-frame Poisson initializer is not a full RSP implementation. */
#include "closeup_align.h"
#include "../common/pipeline_constants.h"
#include "../common/mesh_manifold.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct CaFace {
    double g[6], initial[4], frame[4], dual[4], target[4];
    double area, evidence, projection_weight;
    int invalid_evidence;
} CaFace;

typedef struct CaEdge { int32_t a, b; double weight; } CaEdge;
typedef struct CaHalfEdge { uint64_t key; int32_t a, b; } CaHalfEdge;

typedef struct CaDomain {
    Arena_T arena;
    const float *xyz;
    const int32_t *faces;
    const SlimRefineDirection *directions;
    const uint8_t *protected_faces;
    size_t nv, nf;
    CaFace *face;
    const SlimRefineDirection *solver_directions;
    int32_t *boundary_next, *boundary_root;
    size_t boundary_loops;
} CaDomain;

static double ca_dot(const double *a, const double *b)
{ return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }

static int32_t ca_root(int32_t *p, int32_t v)
{
    while (p[v] != v) { p[v] = p[p[v]]; v = p[v]; }
    return v;
}

static int ca_arguments(Arena_T arena, const float *xyz, size_t nv,
                        const int32_t *faces, size_t nf, const double *uv)
{
    if (!arena || !xyz || !faces || !uv || nv < 3 || !nf ||
        nv > (size_t)INT_MAX/2 || nf > (size_t)INT_MAX/3 ||
        nv > SIZE_MAX/(8*sizeof(double)) ||
        nf > SIZE_MAX/(3*sizeof(CaEdge)) || nf > SIZE_MAX/sizeof(CaFace))
        return -1;
    for (size_t v = 0; v < nv; v++) {
        for (int k = 0; k < 3; k++) if (!isfinite(xyz[3*v+(size_t)k])) return -1;
        if (!isfinite(uv[2*v]) || !isfinite(uv[2*v+1])) return -1;
    }
    for (size_t f = 0; f < nf; f++) for (int k = 0; k < 3; k++)
        if (faces[3*f+(size_t)k] < 0 || (size_t)faces[3*f+(size_t)k] >= nv) return -1;
    return 0;
}

/* Differences keep translated UV charts numerically equivalent. */
static void ca_jacobian(const CaDomain *d, size_t f, const double *uv, double j[4])
{
    const int32_t *tri = d->faces+3*f;
    const double *g = d->face[f].g;
    for (size_t a = 0; a < 2; a++) {
        double u = uv[2*(size_t)tri[1]+a]-uv[2*(size_t)tri[0]+a];
        double v = uv[2*(size_t)tri[2]+a]-uv[2*(size_t)tri[0]+a];
        j[2*a] = u*g[2]+v*g[4];
        j[2*a+1] = u*g[3]+v*g[5];
    }
}

static int ca_setup(CaDomain *d, const double *uv)
{
    for (size_t f = 0; f < d->nf; f++) {
        const int32_t *tri = d->faces+3*f;
        CaFace *q = d->face+f;
        double e[3], b[3], basis[6], len, proj, height;
        for (size_t k = 0; k < 3; k++) {
            e[k] = (double)d->xyz[3*(size_t)tri[1]+k]-d->xyz[3*(size_t)tri[0]+k];
            b[k] = (double)d->xyz[3*(size_t)tri[2]+k]-d->xyz[3*(size_t)tri[0]+k];
        }
        len = sqrt(ca_dot(e,e));
        if (!(len > 0) || !isfinite(len)) return -1;
        for (size_t k = 0; k < 3; k++) basis[k] = e[k]/len;
        proj = ca_dot(b,basis);
        for (size_t k = 0; k < 3; k++) basis[3+k] = b[k]-proj*basis[k];
        height = sqrt(ca_dot(basis+3,basis+3));
        if (!(height > 0) || !isfinite(height)) return -1;
        for (size_t k = 0; k < 3; k++) basis[3+k] /= height;
        q->area = .5*len*height;
        q->g[2] = 1/len; q->g[3] = -proj/(len*height);
        q->g[4] = 0; q->g[5] = 1/height;
        q->g[0] = -q->g[2]; q->g[1] = -q->g[3]-q->g[5];
        ca_jacobian(d,f,uv,q->initial);
        memcpy(q->target,q->initial,sizeof q->target);
        if (!d->directions) continue;
        if (!isfinite(d->directions[f].weight)) { q->invalid_evidence = 1; continue; }
        if (!(d->directions[f].weight > 0)) continue;
        const SlimRefineDirection *arms = d->directions+f;
        int valid = isfinite(arms->weight);
        for (size_t a = 0; valid && a < 2; a++) {
            const double *v = a ? arms->v : arms->h;
            double norm3 = sqrt(ca_dot(v,v));
            double x = ca_dot(v,basis), y = ca_dot(v,basis+3), norm = hypot(x,y);
            if (!isfinite(norm3) || !isfinite(norm) || !(norm > 1e-10*norm3)) valid = 0;
            else { q->frame[a] = x/norm; q->frame[2+a] = y/norm; }
        }
        if (!valid || !(fabs(q->frame[0]*q->frame[3]-q->frame[1]*q->frame[2]) > 1e-8)) {
            q->invalid_evidence = 1; continue;
        }
        q->evidence = arms->weight;
        double frame[4], *j = q->initial;
        memcpy(frame,q->frame,sizeof frame);
        double hx = j[0]*frame[0]+j[1]*frame[2], hy = j[2]*frame[0]+j[3]*frame[2];
        double vx = j[0]*frame[1]+j[1]*frame[3], vy = j[2]*frame[1]+j[3]*frame[3];
        double direct = hy*hy/(hx*hx+hy*hy+DBL_MIN)+vx*vx/(vx*vx+vy*vy+DBL_MIN);
        if (direct > 1) {
            double t = frame[0]; frame[0] = frame[1]; frame[1] = t;
            t = frame[2]; frame[2] = frame[3]; frame[3] = t;
        }
        if (j[0]*frame[0]+j[1]*frame[2] < 0) { frame[0] = -frame[0]; frame[2] = -frame[2]; }
        if (j[2]*frame[1]+j[3]*frame[3] < 0) { frame[1] = -frame[1]; frame[3] = -frame[3]; }
        double det = frame[0]*frame[3]-frame[1]*frame[2];
        /* A poor dual frame may be omitted from initialization, never scoring. */
        if (!(det > .2)) continue;
        q->dual[0] = frame[3]/det; q->dual[1] = -frame[1]/det;
        q->dual[2] = -frame[2]/det; q->dual[3] = frame[0]/det;
        memcpy(q->target,q->dual,sizeof q->target);
        q->projection_weight = arms->weight;
    }
    return 0;
}

static int ca_connected(CaDomain *d)
{
    int32_t *parent = ARENA_ALLOC(d->arena,d->nv*sizeof *parent);
    uint8_t *incident = ARENA_CALLOC(d->arena,d->nv,sizeof *incident);
    for (size_t v = 0; v < d->nv; v++) parent[v] = (int32_t)v;
    for (size_t f = 0; f < d->nf; f++) for (size_t k = 0; k < 3; k++) {
        int32_t v = d->faces[3*f+k], a = ca_root(parent,d->faces[3*f]), b = ca_root(parent,v);
        incident[v] = 1;
        if (a != b) parent[a > b ? a : b] = a < b ? a : b;
    }
    int32_t root = ca_root(parent,0);
    for (size_t v = 0; v < d->nv; v++)
        if (!incident[v] || ca_root(parent,(int32_t)v) != root) return -1;
    return 0;
}

static int ca_halfedge_compare(const void *a, const void *b)
{
    uint64_t x = ((const CaHalfEdge *)a)->key, y = ((const CaHalfEdge *)b)->key;
    return x < y ? -1 : x > y;
}

static int ca_topology(CaDomain *d, CloseupAlignStats *report)
{
    MeshManifoldStats topology = MeshManifold_audit(d->arena,d->nv,d->faces,d->nf);
    if (topology.nm_edges || topology.nm_verts || topology.same_dir_edges || !topology.boundary_edges)
        return -1;
    CaHalfEdge *edges = ARENA_ALLOC(d->arena,3*d->nf*sizeof *edges);
    d->boundary_next = ARENA_ALLOC(d->arena,d->nv*sizeof *d->boundary_next);
    d->boundary_root = ARENA_ALLOC(d->arena,d->nv*sizeof *d->boundary_root);
    uint8_t *seen = ARENA_CALLOC(d->arena,d->nv,sizeof *seen);
    for (size_t v = 0; v < d->nv; v++) d->boundary_next[v] = -1;
    for (size_t f = 0; f < d->nf; f++) for (size_t k = 0; k < 3; k++) {
        int32_t a = d->faces[3*f+k], b = d->faces[3*f+(k+1)%3];
        edges[3*f+k].a = a; edges[3*f+k].b = b;
        edges[3*f+k].key = ((uint64_t)(a < b ? a : b)<<32)|(uint32_t)(a > b ? a : b);
    }
    qsort(edges,3*d->nf,sizeof *edges,ca_halfedge_compare);
    for (size_t e = 0; e < 3*d->nf;) {
        size_t end = e+1;
        while (end < 3*d->nf && edges[end].key == edges[e].key) end++;
        if (end-e == 1) {
            if (d->boundary_next[edges[e].a] >= 0) return -1;
            d->boundary_next[edges[e].a] = edges[e].b;
        }
        e = end;
    }
    size_t traversed = 0;
    for (size_t v = 0; v < d->nv; v++) if (d->boundary_next[v] >= 0 && !seen[v]) {
        int32_t at = (int32_t)v; size_t length = 0;
        do {
            if (at < 0 || (size_t)at >= d->nv || seen[at] || length >= d->nv) return -1;
            seen[at] = 1; length++; at = d->boundary_next[at];
        } while (at != (int32_t)v);
        if (length < 3) return -1;
        traversed += length;
        d->boundary_root[d->boundary_loops++] = (int32_t)v;
    }
    int64_t chi = (int64_t)d->nv-(int64_t)(topology.boundary_edges+topology.manifold_edges)+(int64_t)d->nf;
    if (traversed != topology.boundary_edges || chi != 2-(int64_t)d->boundary_loops) return -1;
    report->boundary_loops = d->boundary_loops;
    report->input_manifold_disk = d->boundary_loops == 1;
    report->input_topology_supported = 1;
    return 0;
}

/* On a consistently oriented manifold, positive triangles and nonintersecting
 * simple boundary cycles are supplemented by boundary degree one: one positive
 * outer cycle, all other cycles negative. Nested negative holes or holes outside
 * that outer cycle would imply negative degree somewhere, impossible for the
 * positive triangles. Thus the degree is at most one. See Lipman (2014),
 * Bijective Mappings of Meshes with Boundary, Theorem 1 / degree argument.
 * Existing native guards still check all boundary pairs and every triangle. */
static int ca_boundary_image(const CaDomain *d, const double *uv)
{
    size_t positive = 0;
    for (size_t loop = 0; loop < d->boundary_loops; loop++) {
        int32_t root = d->boundary_root[loop], at = root;
        const double *origin = uv+2*(size_t)root;
        double area = 0, correction = 0;
        do {
            int32_t next = d->boundary_next[at];
            double ax = uv[2*(size_t)at]-origin[0], ay = uv[2*(size_t)at+1]-origin[1];
            double bx = uv[2*(size_t)next]-origin[0], by = uv[2*(size_t)next+1]-origin[1];
            double term = ax*by-ay*bx-correction, sum = area+term;
            correction = (sum-area)-term; area = sum; at = next;
        } while (at != root);
        if (!isfinite(area) || area == 0) return -1;
        positive += area > 0;
    }
    return positive == 1 ? 0 : -1;
}

static int ca_measure(const CaDomain *d, const double *uv, int orientation, CloseupAlignMetric *out)
{
    CloseupAlignMetric m = {0};
    double angle = 0, sd = 0;
    m.faces = d->nf; m.min_det = m.min_sigma = DBL_MAX;
    m.fiber_rms_degrees = NAN;
    for (size_t f = 0; f < d->nf; f++) {
        const CaFace *q = d->face+f;
        double j[4]; ca_jacobian(d,f,uv,j);
        double det = j[0]*j[3]-j[1]*j[2];
        double hi = .5*(hypot(j[0]+j[3],j[2]-j[1])+hypot(j[0]-j[3],j[2]+j[1]));
        double lo = hi > 0 ? fabs(det)/hi : 0;
        if (!isfinite(det) || !isfinite(hi) || !isfinite(lo)) return -1;
        int bad10 = hi > 1.1 || lo < 1/1.1;
        int bad25 = hi > 1.25 || lo < 1/1.25;
        m.area += q->area;
        m.inverted_faces += !(orientation*det > 0);
        m.min_det = fmin(m.min_det,orientation*det);
        m.min_sigma = fmin(m.min_sigma,lo); m.max_sigma = fmax(m.max_sigma,hi);
        sd += q->area*(hi*hi+lo*lo+1/(hi*hi)+1/(lo*lo));
        if (bad10) m.bad10_area += q->area;
        if (bad25) m.bad25_area += q->area;
        if (!d->protected_faces || d->protected_faces[f]) {
            m.protected_area += q->area;
            if (bad10) m.protected_bad10_area += q->area;
            if (bad25) { m.protected_bad25_area += q->area; m.protected_bad25_faces++; }
        }
        m.invalid_evidence_faces += q->invalid_evidence;
        if (!(q->evidence > 0)) continue;
        const double *a = q->frame;
        double hx = j[0]*a[0]+j[1]*a[2], hy = j[2]*a[0]+j[3]*a[2];
        double vx = j[0]*a[1]+j[1]*a[3], vy = j[2]*a[1]+j[3]*a[3];
        if (!(hypot(hx,hy) > 0 && hypot(vx,vy) > 0)) return -1;
        double ah = atan2(fabs(hy),fabs(hx)), av = atan2(fabs(vx),fabs(vy));
        double bh = atan2(fabs(hx),fabs(hy)), bv = atan2(fabs(vy),fabs(vx));
        double mass = q->area*q->evidence;
        m.direction_energy += mass*fmin(hy*hy+vx*vx,hx*hx+vy*vy);
        angle += mass*.5*fmin(ah*ah+av*av,bh*bh+bv*bv);
        m.confidence_mass += mass; m.measured_area += q->area; m.measured_faces++;
    }
    if (!(m.area > 0) || !isfinite(m.area) || !isfinite(m.confidence_mass) || !isfinite(m.direction_energy)) return -1;
    m.sd_mean = sd/m.area;
    if (m.confidence_mass > 0) m.fiber_rms_degrees = sqrt(angle/m.confidence_mass)*57.2957795130823208768;
    *out = m;
    return 0;
}

int CloseupAlign_measure(Arena_T arena, const float *xyz, size_t nv,
                         const int32_t *faces, size_t nf, const double *uv,
                         const SlimRefineDirection *directions,
                         const uint8_t *protected_faces, int orientation,
                         CloseupAlignMetric *metric)
{
    if (!metric || (orientation != 1 && orientation != -1) || ca_arguments(arena,xyz,nv,faces,nf,uv))
        return CLOSEUP_ALIGN_INVALID;
    Arena_Mark mark = Arena_save(arena);
    volatile int rc = CLOSEUP_ALIGN_INVALID;
    TRY {
        CaDomain d = {arena,xyz,faces,directions,protected_faces,nv,nf,NULL};
        CloseupAlignMetric measured;
        d.face = ARENA_CALLOC(arena,nf,sizeof *d.face);
        if (!ca_setup(&d,uv) && !ca_measure(&d,uv,orientation,&measured)) { *metric = measured; rc = 0; }
    } FINALLY { Arena_restore(arena,mark); } END_TRY;
    return rc;
}

void CloseupAlign_defaults(CloseupAlignOpts *opts)
{
    if (!opts) return;
    memset(opts,0,sizeof *opts);
    SlimRefine_defaults(&opts->slim);
    opts->slim.iterations = 40;
    opts->slim.strip_pack = 0; opts->slim.guard_boundary = 1;
    opts->slim.local_steps = 1; opts->slim.local_step_min_faces = 1;
    opts->slim.local_relaxation = 2; opts->slim.local_relax_boundary = 1;
    opts->field_initialize = 1; opts->target_rms_degrees = CLOSEUP_FIBER_RMS_DEG;
}

static void ca_centroid(const CaDomain *d, const double *uv, double center[2])
{
    double total = 0;
    center[0] = center[1] = 0;
    for (size_t f = 0; f < d->nf; f++) {
        double w = d->face[f].area/3;
        for (size_t k = 0; k < 3; k++) for (size_t a = 0; a < 2; a++)
            center[a] += w*(uv[2*(size_t)d->faces[3*f+k]+a]-uv[a]);
        total += d->face[f].area;
    }
    for (size_t a = 0; a < 2; a++) center[a] = uv[a]+center[a]/total;
}

static void ca_recenter(const CaDomain *d, double *uv, const double center[2])
{
    double current[2]; ca_centroid(d,uv,current);
    for (size_t v = 0; v < d->nv; v++) for (size_t a = 0; a < 2; a++)
        uv[2*v+a] += center[a]-current[a];
}

/* Reuse the exact native feasibility guards; no new/weaker collision test.
 * Zero iterations, no guide, no turn and no pack are a read-only certificate. */
static int ca_slim(CaDomain *d, const CloseupAlignOpts *options,
                    const double *in, double *out, double weight, int measure_only,
                    SlimRefineStats *stats)
{
    Arena_Mark mark = Arena_save(d->arena);
    volatile int rc = -1;
    SlimRefineOpts opts = options->slim;
    opts.directions = measure_only ? NULL : d->solver_directions;
    opts.direction_weight = weight; opts.direction_rigid_gauge = 0;
    opts.homotopy_stages = 1; opts.final_iterations = 0;
    opts.guide_uv = NULL; opts.metric_uv = NULL; opts.guide_warm_start = 0;
    opts.scaffold_uv = NULL; opts.fixed_vertices = NULL;
    opts.strip_pack = 0; opts.guard_boundary = 1;
    if (measure_only) { opts.iterations = 0; opts.local_steps = 0; opts.local_relaxation = 0; opts.verbose = 0; }
    TRY {
        rc = SlimRefine_run_double(d->arena,d->xyz,d->nv,d->faces,d->nf,in,&opts,out,stats);
        if (!rc && (stats->components != 1 || stats->flips_after || stats->boundary_intersections_after)) rc = -1;
        if (!rc && ca_boundary_image(d,out)) rc = -1;
        /* A metric-guard rollback returns SLIM's certified input. Preserve
         * that distinction from a failed solve and keep trying lower weights. */
        if (!rc && stats->reverted) rc = 1;
    } FINALLY { Arena_restore(d->arena,mark); } END_TRY;
    return rc;
}

static int ca_edge_compare(const void *a, const void *b)
{
    const CaEdge *x = a, *y = b;
    if (x->a != y->a) return x->a < y->a ? -1 : 1;
    return x->b < y->b ? -1 : x->b > y->b;
}

static int ca_project_body(CaDomain *d, const CloseupAlignOpts *opts, double *uv)
{
    const double regularization = .05;
    size_t cap = 3*d->nf, ne = 0, at = 0;
    CaEdge *edges = ARENA_ALLOC(d->arena,cap*sizeof *edges);
    double *diagonal = ARENA_CALLOC(d->arena,3*d->nv,sizeof *diagonal);
    double *rhs = ARENA_ALLOC(d->arena,2*d->nv*sizeof *rhs);
    for (size_t f = 0; f < d->nf; f++) {
        const CaFace *q = d->face+f;
        double weight = q->area*(regularization+q->projection_weight);
        for (size_t a = 0; a < 3; a++) {
            int32_t x = d->faces[3*f+a];
            if (x) {
                double w = weight*(q->g[2*a]*q->g[2*a]+q->g[2*a+1]*q->g[2*a+1]);
                diagonal[3*(size_t)x] += w; diagonal[3*(size_t)x+2] += w;
            }
            for (size_t b = 0; b < a; b++) {
                int32_t y = d->faces[3*f+b];
                if (!x || !y) continue;
                edges[at].a = x < y ? x : y; edges[at].b = x > y ? x : y;
                edges[at++].weight = weight*(q->g[2*a]*q->g[2*b]+q->g[2*a+1]*q->g[2*b+1]);
            }
        }
    }
    diagonal[0] = diagonal[2] = 1;
    qsort(edges,at,sizeof *edges,ca_edge_compare);
    for (size_t k = 0; k < at; k++) {
        if (ne && edges[k].a == edges[ne-1].a && edges[k].b == edges[ne-1].b)
            edges[ne-1].weight += edges[k].weight;
        else edges[ne++] = edges[k];
    }
    int32_t *edge_a = ARENA_ALLOC(d->arena,(ne ? ne : 1)*sizeof *edge_a);
    int32_t *edge_b = ARENA_ALLOC(d->arena,(ne ? ne : 1)*sizeof *edge_b);
    double *blocks = ARENA_CALLOC(d->arena,3*(ne ? ne : 1),sizeof *blocks);
    for (size_t e = 0; e < ne; e++) {
        edge_a[e] = edges[e].a; edge_b[e] = edges[e].b;
        blocks[3*e] = blocks[3*e+2] = edges[e].weight;
    }
    for (int pass = 0; pass < 8; pass++) {
        memset(rhs,0,2*d->nv*sizeof *rhs);
        for (size_t f = 0; f < d->nf; f++) {
            CaFace *q = d->face+f;
            if (pass && q->projection_weight > 0) {
                double j[4]; ca_jacobian(d,f,uv,j);
                for (size_t a = 0; a < 2; a++) {
                    const double *r = q->dual+2*a;
                    double scale = (r[0]*j[2*a]+r[1]*j[2*a+1]+1)/(r[0]*r[0]+r[1]*r[1]+1);
                    scale = fmax(.5,fmin(2,scale));
                    q->target[2*a] = scale*r[0]; q->target[2*a+1] = scale*r[1];
                }
            }
            for (size_t k = 0; k < 3; k++) for (size_t a = 0; a < 2; a++) {
                size_t v = (size_t)d->faces[3*f+k];
                for (size_t b = 0; b < 2; b++) rhs[2*v+a] += q->area*q->g[2*k+b]*
                    (q->projection_weight*q->target[2*a+b]+regularization*q->initial[2*a+b]);
            }
        }
        rhs[0] = rhs[1] = 0;
        Arena_Mark mark = Arena_save(d->arena);
        volatile int rc = -1;
        TRY {
            rc = opts->slim.linear_solve(d->arena,opts->slim.linear_context,d->nv,ne,diagonal,edge_a,edge_b,blocks,rhs,uv);
        } FINALLY { Arena_restore(d->arena,mark); } END_TRY;
        if (rc) return -1;
        for (size_t k = 0; k < 2*d->nv; k++) if (!isfinite(uv[k])) return -1;
    }
    return 0;
}

static int ca_project(CaDomain *d, const CloseupAlignOpts *opts, double *uv)
{
    volatile int rc = -1, acquired = 0;
    TRY {
        if (!opts->projection_begin) rc = ca_project_body(d,opts,uv);
        else if (!opts->projection_begin(opts->slim.linear_context)) {
            acquired = 1;
            rc = ca_project_body(d,opts,uv);
        }
    } FINALLY {
        if (acquired) opts->projection_end(opts->slim.linear_context);
    } END_TRY;
    return rc;
}

static int ca_nonregression(const CloseupAlignMetric *m, const CloseupAlignMetric *base)
{
    double tolerance = 1e-10*fmax(base->area,1);
    return !m->inverted_faces && isfinite(m->sd_mean) &&
        m->sd_mean <= base->sd_mean*(1+1e-10) &&
        m->bad10_area <= base->bad10_area+tolerance && m->bad25_area <= base->bad25_area+tolerance &&
        m->protected_bad10_area <= base->protected_bad10_area+tolerance &&
        m->protected_bad25_area <= base->protected_bad25_area+tolerance &&
        m->protected_bad25_faces <= base->protected_bad25_faces;
}

static int ca_better(const CloseupAlignMetric *m, const CloseupAlignMetric *best, double target)
{
    if (!isfinite(m->fiber_rms_degrees)) return 0;
    int meet = m->fiber_rms_degrees <= target, best_meet = best->fiber_rms_degrees <= target;
    if (meet != best_meet) return meet;
    if (m->fiber_rms_degrees < best->fiber_rms_degrees-1e-7) return 1;
    return fabs(m->fiber_rms_degrees-best->fiber_rms_degrees) <= 1e-7 && m->sd_mean < best->sd_mean-1e-10;
}

static int ca_local_budget(const CloseupAlignMetric *m,const CloseupAlignMetric *base)
{
    double tolerance=1e-10*fmax(1,base->area);
    return !m->inverted_faces&&isfinite(m->sd_mean)&&
        m->bad10_area<=.05*m->area+tolerance&&
        m->protected_bad10_area<=.05*m->protected_area+tolerance&&
        m->bad25_area<=base->bad25_area+tolerance&&
        m->protected_bad25_area<=base->protected_bad25_area+tolerance&&
        m->protected_bad25_faces<=base->protected_bad25_faces;
}

/* A centered solve certificate does not survive rounded translation by
 * assumption. Select only the exact bytes that will be returned to the caller.
 * Keep trial/anchor in the centered solver frame; best is in the output frame. */
static int ca_output_candidate(CaDomain *d, const CloseupAlignOpts *opts,
                               const double *trial, const double center[2],
                               int orientation, double *best,
                               CloseupAlignStats *report, CloseupAlignTrial *record)
{
    Arena_Mark mark = Arena_save(d->arena);
    volatile int feasible = 0;
    TRY {
        double *restored = ARENA_ALLOC(d->arena,2*d->nv*sizeof *restored);
        double *checked = ARENA_ALLOC(d->arena,2*d->nv*sizeof *checked);
        for (size_t v = 0; v < d->nv; v++) {
            restored[2*v] = trial[2*v]+center[0];
            restored[2*v+1] = orientation*trial[2*v+1]+center[1];
        }
        ca_recenter(d,restored,center);
        SlimRefineStats certificate = {0};
        if (!ca_measure(d,restored,orientation,&record->metric) &&
            !ca_slim(d,opts,restored,checked,0,1,&certificate)) {
            feasible = record->feasible = 1;
            record->scale_nonregression = ca_nonregression(&record->metric,&report->before);
            int metric_ok=record->scale_nonregression||(record->seed==3&&ca_local_budget(&record->metric,&report->before));
            if (metric_ok)
                record->candidate_guard_pass = !opts->candidate_guard ||
                    opts->candidate_guard(opts->candidate_context,d->nv,restored) > 0;
            if (metric_ok && record->candidate_guard_pass &&
                (record->seed!=2 || record->metric.sd_mean*record->metric.area+
                    CLOSEUP_FIBER_DIRECTION_WEIGHT*record->metric.direction_energy <=
                    report->after.sd_mean*report->after.area+
                    CLOSEUP_FIBER_DIRECTION_WEIGHT*report->after.direction_energy+1e-8) &&
                ca_better(&record->metric,&report->after,opts->target_rms_degrees)) {
                for (size_t k = 0; k < report->trials; k++) report->trial[k].selected = 0;
                record->selected = 1;
                memcpy(best,restored,2*d->nv*sizeof *best);
                report->after = record->metric; report->selected_seed = record->seed;
                report->selected_direction_weight = record->direction_weight;
            }
        } else if (opts->slim.verbose) {
            fprintf(stderr,"[closeup align] seed=%d weight=%.9g output-frame guard refused\n",
                    record->seed,record->direction_weight);
        }
    } FINALLY { Arena_restore(d->arena,mark); } END_TRY;
    return feasible;
}

static int ca_trial(CaDomain *d, const CloseupAlignOpts *opts, const double *seed,
                     int seed_id, double weight, const double center[2], int orientation,
                     double *trial, double *best, CloseupAlignStats *report)
{
    if (report->trials >= CLOSEUP_ALIGN_MAX_TRIALS) return -1;
    CloseupAlignTrial *record = report->trial+report->trials++;
    SlimRefineStats slim = {0};
    record->seed = seed_id; record->direction_weight = weight;
    record->candidate_guard_pass = opts->candidate_guard ? 0 : 1;
    record->metric.fiber_rms_degrees = NAN;
    record->solver_status = ca_slim(d,opts,seed,trial,weight,0,&slim);
    record->iterations_run = slim.iterations_run;
    record->accepted_steps = slim.accepted_steps;
    if (record->solver_status < 0) return -1;
    double zero[2] = {0,0};
    ca_recenter(d,trial,zero);
    /* A restored-frame refusal does not corrupt the certified centered anchor
     * or erase an earlier output-feasible winner. It is not solver failure. */
    (void)ca_output_candidate(d,opts,trial,center,orientation,best,report,record);
    /* Positive status is a feasible rollback, not a new proposal
     * that may end the bounded lower-weight search. */
    return record->solver_status;
}

#include "closeup_align_rigid.inc"
#ifdef VES_CLOSEUP_LOCAL_GLOBAL
#include "closeup_align_local_global.inc"
#endif

static int ca_propose(CaDomain *d, const double *in, const CloseupAlignOpts *opts,
                       double *out, CloseupAlignStats *report, int rigid_only)
{
    d->face = ARENA_CALLOC(d->arena,d->nf,sizeof *d->face);
    if (ca_setup(d,in) || ca_connected(d)) return CLOSEUP_ALIGN_INVALID;
    if (ca_topology(d,report)) return CLOSEUP_ALIGN_TOPOLOGY;
    double j[4]; ca_jacobian(d,0,in,j);
    int orientation = j[0]*j[3]-j[1]*j[2] > 0 ? 1 : -1;
    if (ca_measure(d,in,orientation,&report->before)) return CLOSEUP_ALIGN_INVALID;
    if (report->before.inverted_faces) return CLOSEUP_ALIGN_INFEASIBLE;
    double center[2]; ca_centroid(d,in,center);
    double *initial = ARENA_ALLOC(d->arena,2*d->nv*sizeof *initial);
    double *best = ARENA_ALLOC(d->arena,2*d->nv*sizeof *best);
    double *field = ARENA_ALLOC(d->arena,2*d->nv*sizeof *field);
    double *trial = ARENA_ALLOC(d->arena,2*d->nv*sizeof *trial);
    double *anchor = ARENA_ALLOC(d->arena,2*d->nv*sizeof *anchor);
    double zero[2] = {0,0};
    for (size_t v = 0; v < d->nv; v++) {
        initial[2*v] = in[2*v]-center[0];
        initial[2*v+1] = orientation*(in[2*v+1]-center[1]);
    }
    /* Rebuild the target branches in positive-handed, centered coordinates. */
    memset(d->face,0,d->nf*sizeof *d->face);
    if (ca_setup(d,initial)) return CLOSEUP_ALIGN_INVALID;
    SlimRefineStats certificate = {0};
    if (ca_slim(d,opts,initial,trial,0,1,&certificate)) return CLOSEUP_ALIGN_INFEASIBLE;
    /* The unchanged fallback is returned in this frame, not the centered one. */
    if (ca_slim(d,opts,in,trial,0,1,&certificate)) return CLOSEUP_ALIGN_INFEASIBLE;
    report->input_feasible = 1;
    report->after = report->before; report->selected_seed = -1;
    memcpy(best,in,2*d->nv*sizeof *best);
    if (!report->before.measured_faces) {
        memcpy(out,in,2*d->nv*sizeof *out);
        report->scale_nonregression = 1;
        return CLOSEUP_ALIGN_NO_EVIDENCE;
    }
    /* Invalid/absent physical arms stay unknown. Never synthesize evidence
     * from UV, and never mutate the caller's pinned input array. */
    SlimRefineDirection *usable = ARENA_ALLOC(d->arena,d->nf*sizeof *usable);
    memcpy(usable,d->directions,d->nf*sizeof *usable);
    for (size_t f = 0; f < d->nf; f++) if (!(d->face[f].evidence > 0)) usable[f].weight = 0;
    d->solver_directions = usable;
    if (rigid_only) {
        ca_try_rigid(d,opts,initial,2,center,orientation,trial,best,report);
        goto rigid_finalize;
    }
    if (opts->field_initialize) {
        report->field_attempted = 1;
        for (size_t f = 0; f < d->nf; f++) report->field_faces += d->face[f].projection_weight > 0;
        memcpy(field,initial,2*d->nv*sizeof *field);
        report->field_status = ca_project(d,opts,field);
        if (!report->field_status) {
            ca_recenter(d,field,zero);
            report->field_status = ca_measure(d,field,1,&report->field);
            if (!report->field_status && !report->field.inverted_faces)
                report->field_status = ca_slim(d,opts,field,trial,0,1,&certificate);
            else if (!report->field_status) report->field_status = CLOSEUP_ALIGN_INFEASIBLE;
            report->field_feasible = !report->field_status;
        }
    }
    int successes = 0;
    /* Bounded, component-independent policy. Every lower-weight attempt starts
     * from the same weight-64 solution, retaining an honest old-seed control. */
    const double weights[4] = {64,16,4,1};
    for (int route = 0; route < 2; route++) {
        if (route && !report->field_feasible) continue;
        const double *seed = route ? field : initial;
        double high = 0, low = 0;
        int previous_failed = 0;
        for (size_t w = 0; w < 4; w++) {
            int rc = ca_trial(d,opts,w ? anchor : seed,route,weights[w],center,orientation,trial,best,report);
            if (!w) memcpy(anchor,rc ? seed : trial,2*d->nv*sizeof *anchor);
            if (rc >= 0) successes++;
            const CloseupAlignTrial *record = report->trial+report->trials-1;
            int valid = !rc && record->scale_nonregression && record->candidate_guard_pass;
            if (!valid) { high = weights[w]; previous_failed = 1; }
            else {
                if (previous_failed) low = weights[w];
                /* Highest scale-safe member found. Refine its bracket; later
                 * coarse weights cannot improve a strictly decreasing weight
                 * penalty's accuracy systematically, but every tested result
                 * remains selected by actual measured metrics. */
                break;
            }
        }
        if (low > 0 && high > low) for (int refine = 0; refine < 3; refine++) {
            double weight = sqrt(low*high);
            int rc = ca_trial(d,opts,anchor,route,weight,center,orientation,trial,best,report);
            if (rc >= 0) successes++;
            const CloseupAlignTrial *record = report->trial+report->trials-1;
            if (!rc && record->scale_nonregression && record->candidate_guard_pass) low = weight;
            else high = weight;
        }
    }
    if (!successes) return CLOSEUP_ALIGN_SOLVER_FAILED;
#ifdef VES_CLOSEUP_LOCAL_GLOBAL
    if(opts->local_global_iterations>0) {
        for(size_t v=0;v<d->nv;v++) {
            anchor[2*v]=best[2*v]-center[0];anchor[2*v+1]=orientation*(best[2*v+1]-center[1]);
        }
        ca_try_rigid(d,opts,anchor,2,center,orientation,trial,best,report);
        ca_try_local_global(d,opts,in,initial,center,orientation,trial,best,report);
    }
#endif
rigid_finalize:
    report->scale_nonregression = ca_nonregression(&report->after,&report->before);
    report->target_met = report->after.fiber_rms_degrees < opts->target_rms_degrees;
    if (report->selected_seed < 0) { memcpy(out,in,2*d->nv*sizeof *out); return CLOSEUP_ALIGN_UNCHANGED; }
    /* best already contains the exact output-frame bytes checked at selection.
     * Retain the final certificate; do not translate/recenter a second time. */
    if (ca_slim(d,opts,best,trial,0,1,&certificate) ||
        ca_measure(d,best,orientation,&report->after)) return CLOSEUP_ALIGN_INFEASIBLE;
    memcpy(out,best,2*d->nv*sizeof *out);
    return CLOSEUP_ALIGN_PROPOSED;
}

static int ca_overlap(const void *a, size_t an, const void *b, size_t bn)
{
    if (!a || !b) return 0;
    uintptr_t x = (uintptr_t)a, y = (uintptr_t)b;
    return x <= y ? y-x < an : x-y < bn;
}

static int ca_run(Arena_T arena, const float *xyz, size_t nv,
                         const int32_t *faces, size_t nf, const double *uv_in,
                         const SlimRefineDirection *directions,
                         const uint8_t *protected_faces,
                         const CloseupAlignOpts *opts, double *uv_out,
                         CloseupAlignStats *stats, int rigid_only)
{
    if (!opts || !uv_out || !stats || ca_arguments(arena,xyz,nv,faces,nf,uv_in) ||
        !opts->slim.linear_solve || opts->slim.iterations <= 0 ||
        (!!opts->projection_begin != !!opts->projection_end) ||
        opts->slim.guide_uv || opts->slim.metric_uv || opts->slim.scaffold_uv || opts->slim.fixed_vertices ||
        !(opts->target_rms_degrees > 0 && opts->target_rms_degrees <= 45) ||
        ca_overlap(uv_out,2*nv*sizeof *uv_out,uv_in,2*nv*sizeof *uv_in) ||
        ca_overlap(uv_out,2*nv*sizeof *uv_out,xyz,3*nv*sizeof *xyz) ||
        ca_overlap(uv_out,2*nv*sizeof *uv_out,faces,3*nf*sizeof *faces) ||
        ca_overlap(uv_out,2*nv*sizeof *uv_out,directions,nf*sizeof *directions) ||
        ca_overlap(uv_out,2*nv*sizeof *uv_out,protected_faces,nf)) return CLOSEUP_ALIGN_INVALID;
    Arena_Mark mark = Arena_save(arena);
    volatile int rc = CLOSEUP_ALIGN_INVALID;
    TRY {
        CaDomain d = {arena,xyz,faces,directions,protected_faces,nv,nf,NULL};
        CloseupAlignStats report = {0};
        double *proposal = ARENA_ALLOC(arena,2*nv*sizeof *proposal);
        report.before.fiber_rms_degrees = report.field.fiber_rms_degrees = report.after.fiber_rms_degrees = NAN;
        rc = ca_propose(&d,uv_in,opts,proposal,&report,rigid_only);
        report.status = rc;
        if (rc >= 0) memcpy(uv_out,proposal,2*nv*sizeof *uv_out);
        *stats = report;
    } FINALLY { Arena_restore(arena,mark); } END_TRY;
    return rc;
}

int CloseupAlign_propose(Arena_T arena, const float *xyz, size_t nv,
                         const int32_t *faces, size_t nf, const double *uv_in,
                         const SlimRefineDirection *directions,
                         const uint8_t *protected_faces,
                         const CloseupAlignOpts *opts, double *uv_out,
                         CloseupAlignStats *stats)
{
    return ca_run(arena,xyz,nv,faces,nf,uv_in,directions,protected_faces,opts,uv_out,stats,0);
}

int CloseupAlign_rotate(Arena_T arena, const float *xyz, size_t nv,
                         const int32_t *faces, size_t nf, const double *uv_in,
                         const SlimRefineDirection *directions,
                         const uint8_t *protected_faces,
                         const CloseupAlignOpts *opts, double *uv_out,
                         CloseupAlignStats *stats)
{
    return ca_run(arena,xyz,nv,faces,nf,uv_in,directions,protected_faces,opts,uv_out,stats,1);
}

#include "closeup_align_rigid_test.inc"
