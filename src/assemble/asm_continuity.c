#include "asm_continuity.h"
#include "../common/pipeline_constants.h"
#include "asm_boundary.h"
#include "asm_flatten.h"
#include "asm_metric.h"
#include "asm_relate.h"
#include "../common/union_find.h"
#include "../common/ves_platform.h"
#include "../common/kdtree.h"
#include "../flatten/sparse_solve.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include "../common/ves_omp.h"

/* A displacement field has two coordinates but one scalar mesh operator. */
typedef struct ActVec { double x, y; } ActVec;
typedef struct ActSpring { int32_t a, b; float weight, length; } ActSpring;
typedef struct ActMatch {
    ActVec target;
    double coefficient[4], weight;
    int32_t id[4];
    int count;
} ActMatch;

static void act_match_add(ActMatch *m, int32_t id, double coefficient)
{
    if (id < 0 || coefficient == 0) return;
    for (int i = 0; i < m->count; i++) if (m->id[i] == id) {
        m->coefficient[i] += coefficient; return;
    }
    assert(m->count < 4);
    m->id[m->count] = id; m->coefficient[m->count++] = coefficient;
}

static ActMatch act_match_pair(int32_t a, int32_t b, double x, double y)
{
    ActMatch m = {0}; m.target = (ActVec){x,y}; m.weight = 1;
    act_match_add(&m,a,1); act_match_add(&m,b,-1); return m;
}

static void act_point(const AsmChart *c, int32_t v, const float *gap, double p[2])
{
    if (c->placed_uv) {
        AsmChart_point(c,(size_t)v,p);
        if (gap) {
            double u = gap[0], w = gap[1], ct = cos(c->pose_theta), sn = sin(c->pose_theta);
            if (c->flags & ASM_CHART_MIRROR) u = -u;
            p[0] += ct*u-sn*w; p[1] += sn*u+ct*w;
        }
        return;
    }
    double u = c->uv[(size_t)v*2], w = c->uv[(size_t)v*2+1];
    if (gap) { u += gap[0]; w += gap[1]; }
    if (c->flags & ASM_CHART_MIRROR) u = -u;
    double ct = cos(c->pose_theta), sn = sin(c->pose_theta);
    p[0] = ct*u - sn*w + c->pose_x;
    p[1] = sn*u + ct*w + c->pose_y;
}

typedef struct ActFrame { double map[2][3]; } ActFrame;

/* Remove local UV stretch while retaining its two tangent directions. */
static int act_polar_frame(const ActFrame *j, ActFrame *out)
{
    double aa = 0, ab = 0, bb = 0;
    for (int q = 0; q < 3; q++) {
        aa += j->map[0][q]*j->map[0][q]; ab += j->map[0][q]*j->map[1][q]; bb += j->map[1][q]*j->map[1][q];
    }
    double dd = aa*bb-ab*ab;
    if (!(dd > 1e-12*aa*bb)) return 0;
    dd = sqrt(dd);
    double den = dd*sqrt(aa+bb+2*dd);
    if (!(den > 0) || !isfinite(den)) return 0;
    for (int q = 0; q < 3; q++) {
        out->map[0][q] = ((bb+dd)*j->map[0][q]-ab*j->map[1][q])/den;
        out->map[1][q] = ((aa+dd)*j->map[1][q]-ab*j->map[0][q])/den;
        if (!isfinite(out->map[0][q]) || !isfinite(out->map[1][q])) return 0;
    }
    return 1;
}

/* Solve the 3-D triangle Gram system, then carry its barycentric
 * differential into an orthonormal UV tangent frame. */
static int act_face_frame(const AsmChart *c, int32_t face, ActFrame *out, double *area)
{
    if (face < 0 || (size_t)face >= c->nf) return 0;
    const int32_t *f = &c->faces[(size_t)face*3];
    for (int k = 0; k < 3; k++) if (f[k] < 0 || (size_t)f[k] >= c->nv) return 0;
    double a = 0, b = 0, e = 0;
    for (int k = 0; k < 3; k++) {
        double x = c->xyz[(size_t)f[1]*3+k] - c->xyz[(size_t)f[0]*3+k];
        double y = c->xyz[(size_t)f[2]*3+k] - c->xyz[(size_t)f[0]*3+k];
        a += x*x; b += x*y; e += y*y;
    }
    double det = a*e - b*b;
    if (!(det > 1e-10*a*e)) return 0;
    ActFrame j;
    for (int k = 0; k < 2; k++) for (int q = 0; q < 3; q++) {
        double x = c->uv[(size_t)f[1]*2+k] - c->uv[(size_t)f[0]*2+k];
        double y = c->uv[(size_t)f[2]*2+k] - c->uv[(size_t)f[0]*2+k];
        double e1 = c->xyz[(size_t)f[1]*3+q] - c->xyz[(size_t)f[0]*3+q];
        double e2 = c->xyz[(size_t)f[2]*3+q] - c->xyz[(size_t)f[0]*3+q];
        j.map[k][q] = (x*(e*e1-b*e2)+y*(a*e2-b*e1))/det;
    }
    if (area) *area = 0.5*sqrt(det);
    return act_polar_frame(&j,out);
}

static int act_frame_gap(const ActFrame *frame, const double d[3], float uv[2])
{
    for (int k = 0; k < 2; k++) {
        double v = 0;
        for (int q = 0; q < 3; q++) v += frame->map[k][q]*d[q];
        uv[k] = (float)v; if (!isfinite(uv[k])) return 0;
    }
    return 1;
}

static int act_gap(const AsmChart *c, int32_t face, const double d[3], float uv[2])
{
    ActFrame frame;
    return act_face_frame(c,face,&frame,NULL) && act_frame_gap(&frame,d,uv);
}

static int act_face_before(const AsmChart *c, int32_t a, int32_t b)
{
    if (b < 0) return 1;
    int32_t x[3], y[3]; memcpy(x,c->faces+3*(size_t)a,sizeof x); memcpy(y,c->faces+3*(size_t)b,sizeof y);
    for (int i = 0; i < 2; i++) for (int j = i+1; j < 3; j++) {
        if (x[i] > x[j]) { int32_t t = x[i]; x[i] = x[j]; x[j] = t; }
        if (y[i] > y[j]) { int32_t t = y[i]; y[i] = y[j]; y[j] = t; }
    }
    for (int i = 0; i < 3; i++) if (x[i] != y[i]) return x[i] < y[i];
    return 0;
}

/* A correspondence is at a vertex, not in whichever face happens to be
 * stored first. Average valid incident tangent frames by physical area,
 * then re-orthonormalize. Storage is released after this one chart. A
 * cancelled average uses the largest valid face, with a topology-based
 * tie break, so averaging cannot erase an existing source witness. */
static ActFrame *act_vertex_frames(Arena_T arena, const AsmChart *c, uint8_t *valid)
{
    ActFrame *frame = ARENA_CALLOC(arena,c->nv,sizeof *frame);
    double *mass = ARENA_CALLOC(arena,c->nv,sizeof *mass), *largest = ARENA_CALLOC(arena,c->nv,sizeof *largest);
    int32_t *best = ARENA_ALLOC(arena,c->nv*sizeof *best);
    for (size_t v = 0; v < c->nv; v++) best[v] = -1;
    for (size_t f = 0; f < c->nf; f++) {
        ActFrame face; double area;
        if (!act_face_frame(c,(int32_t)f,&face,&area) || !(area > 0) || !isfinite(area)) continue;
        for (int k = 0; k < 3; k++) {
            int32_t v = c->faces[3*f+k];
            if (!c->boundary || !c->boundary[v]) continue;
            mass[v] += area;
            for (int u = 0; u < 2; u++) for (int q = 0; q < 3; q++) frame[v].map[u][q] += area*face.map[u][q];
            if (area > largest[v] || (area == largest[v] && act_face_before(c,(int32_t)f,best[v]))) {
                largest[v] = area; best[v] = (int32_t)f;
            }
        }
    }
    for (size_t v = 0; v < c->nv; v++) if (best[v] >= 0) {
        ActFrame average = frame[v];
        for (int u = 0; u < 2; u++) for (int q = 0; q < 3; q++) average.map[u][q] /= mass[v];
        valid[v] = (uint8_t)(act_polar_frame(&average,&frame[v]) || act_face_frame(c,best[v],&frame[v],NULL));
    }
    return frame;
}

int AsmContinuity_chart_gaps(Arena_T arena, const AsmChart *chart,
                             const int32_t *vertices, const double *displacements,
                             size_t count, float *gaps, uint8_t *valid)
{
    if (!arena || !chart || !chart->xyz || !chart->uv || !chart->faces ||
        !chart->boundary || !chart->nv || !chart->nf ||
        (count && (!vertices || !displacements || !gaps || !valid))) return -1;
    for (size_t k = 0; k < count; k++) {
        if (vertices[k] < 0 || (size_t)vertices[k] >= chart->nv) return -1;
        for (int q = 0; q < 3; q++) if (!isfinite(displacements[3*k+q])) return -1;
    }
    if (!count) return 0;
    Arena_Mark mark = Arena_save(arena);
    uint8_t *frame_valid = ARENA_CALLOC(arena,chart->nv,1);
    ActFrame *frames = act_vertex_frames(arena,chart,frame_valid);
    for (size_t k = 0; k < count; k++) {
        int32_t v = vertices[k];
        gaps[2*k] = gaps[2*k+1] = 0;
        valid[k] = (uint8_t)(frame_valid[v] && act_frame_gap(&frames[v],displacements+3*k,gaps+2*k));
    }
    Arena_restore(arena,mark);
    return 0;
}

/* Select the largest stable original incident face, as in the independent
 * Python evaluator. Source geometry and the original vertex frame determine
 * this map once; subsequent UV changes must never redefine the observation. */
static AsmTrim *act_trim_bases(Arena_T arena, const AsmChart *c, const ActFrame *frame, const uint8_t *valid)
{
    AsmTrim *trim = ARENA_CALLOC(arena,c->nv,sizeof *trim);
    double *largest = ARENA_CALLOC(arena,c->nv,sizeof *largest);
    for (size_t f = 0; f < c->nf; f++) {
        const int32_t *vertices = c->faces+3*f;
        for (int k = 0; k < 3; k++) if (vertices[k] < 0 || (size_t)vertices[k] >= c->nv) return NULL;
        if (f > INT32_MAX) return NULL;
        double e[2][3], cross[3];
        for (int k = 0; k < 2; k++) for (int q = 0; q < 3; q++)
            e[k][q] = (double)c->xyz[3*(size_t)vertices[k+1]+q]-c->xyz[3*(size_t)vertices[0]+q];
        for (int q = 0; q < 3; q++) cross[q] = e[0][(q+1)%3]*e[1][(q+2)%3]-e[0][(q+2)%3]*e[1][(q+1)%3];
        double area = hypot(hypot(cross[0],cross[1]),cross[2]);
        for (int k = 0; k < 3; k++) {
            int32_t v = vertices[k];
            if (!valid[v] || !(area > largest[v])) continue;
            double b[4] = {0};
            for (int row = 0; row < 2; row++) for (int col = 0; col < 2; col++)
                for (int q = 0; q < 3; q++) b[2*row+col] += frame[v].map[row][q]*e[col][q];
            double det = b[0]*b[3]-b[1]*b[2];
            double hi = .5*(hypot(b[0]+b[3],b[2]-b[1])+hypot(b[0]-b[3],b[2]+b[1]));
            if (!(hi > 0) || !isfinite(hi) || !(fabs(det) >= 1e-8*hi*hi)) continue;
            AsmTrim t = {0}; t.face = (int32_t)f; t.valid = 1;
            t.inverse[0] = b[3]/det; t.inverse[1] = -b[1]/det;
            t.inverse[2] = -b[2]/det; t.inverse[3] = b[0]/det;
            trim[v] = t; largest[v] = area;
        }
    }
    return trim;
}

static int act_trim_weights(const AsmChart *c, int32_t vertex, const AsmTrim *trim,
                            const float gap[2], double sign, int32_t ids[3], double weights[3])
{
    if (!trim->valid) { ids[0] = vertex; weights[0] = 1; return 1; }
    if (trim->face < 0 || (size_t)trim->face >= c->nf) return 0;
    const int32_t *face = c->faces+3*(size_t)trim->face;
    double a = sign*(trim->inverse[0]*gap[0]+trim->inverse[1]*gap[1]);
    double b = sign*(trim->inverse[2]*gap[0]+trim->inverse[3]*gap[1]);
    int found = 0;
    weights[0] = -a-b; weights[1] = a; weights[2] = b;
    for (int k = 0; k < 3; k++) {
        ids[k] = face[k];
        if (ids[k] < 0 || (size_t)ids[k] >= c->nv || !isfinite(weights[k])) return 0;
        if (ids[k] == vertex) { weights[k] += 1; found++; }
    }
    return found == 1 ? 3 : 0;
}

static int act_trim_point(const AsmChart *c, int32_t vertex, const AsmTrim *trim,
                          const float gap[2], double sign, double p[2])
{
    if (!c->uv || vertex < 0 || (size_t)vertex >= c->nv ||
        !isfinite(gap[0]) || !isfinite(gap[1])) return 0;
    if (!trim->valid) {
        /* Explicit point/offset constraints used by callers without source
         * triangles. prepare() fills every valid physical source witness. */
        float offset[2] = {(float)(sign*gap[0]),(float)(sign*gap[1])};
        act_point(c,vertex,offset,p); return 1;
    }
    int32_t ids[3]; double weights[3], local[2] = {0};
    int count = act_trim_weights(c,vertex,trim,gap,sign,ids,weights);
    if (!count) return 0;
    if (c->placed_uv) {
        p[0] = p[1] = 0;
        for (int k = 0; k < count; k++) for (int q = 0; q < 2; q++) p[q] += weights[k]*c->placed_uv[2*(size_t)ids[k]+q];
        return isfinite(p[0]) && isfinite(p[1]);
    }
    for (int k = 0; k < count; k++) for (int q = 0; q < 2; q++) local[q] += weights[k]*c->uv[2*(size_t)ids[k]+q];
    if (c->flags & ASM_CHART_MIRROR) local[0] = -local[0];
    double ct = cos(c->pose_theta), sn = sin(c->pose_theta);
    p[0] = ct*local[0]-sn*local[1]+c->pose_x; p[1] = sn*local[0]+ct*local[1]+c->pose_y;
    return isfinite(p[0]) && isfinite(p[1]);
}

static double act_dist(const AsmChart *c, int32_t a, int32_t b)
{
    double d2 = 0;
    for (int k = 0; k < 3; k++) {
        double d = c->xyz[(size_t)a*3+k] - c->xyz[(size_t)b*3+k]; d2 += d*d;
    }
    return sqrt(d2);
}

/* Sparse accepted seams are still placement observations. A missing
 * six-sample run is an unresolved certificate, not geometric evidence that
 * the clean2 layout should be pulled apart. Measure every original pair;
 * the temporary run label permits evaluation only and never certifies it. */
static int act_incumbent_metric(const AsmRun *run, const AsmRelation *r, double *rms)
{
    *rms = 1e300;
    if (r->corr_first < 0 || r->corr_count < 3 ||
        (size_t)r->corr_first+(size_t)r->corr_count > run->n_corr ||
        (r->flags & (ASM_REL_DROPPED | ASM_REL_CONTACT))) return 0;
    const AsmChart *a = &run->charts[r->a], *b = &run->charts[r->b];
    if (!AsmChart_in_layout(a) || !AsmChart_in_layout(b) ||
        ((a->flags ^ b->flags) & ASM_CHART_MIRROR) !=
        ((r->flags & ASM_REL_PARITY) ? ASM_CHART_MIRROR : 0)) return 0;
    double squared = 0, gate = ASM_CONTINUITY_GATE(r->rms);
    for (int32_t k = 0; k < r->corr_count; k++) {
        AsmCorr point = run->corr[(size_t)r->corr_first+k];
        if (point.valid != 3) return 0;
        point.run = 0;
        double error = AsmContinuity_pair_error(run,r,&point);
        if (!isfinite(error) || error > gate) return 0;
        squared += error*error;
    }
    *rms = sqrt(squared/(double)r->corr_count);
    return isfinite(*rms);
}

static int act_preserves_incumbent(const AsmRun *run, const uint8_t *changed)
{
    for (size_t i = 0; i < run->n_rels; i++) {
        const AsmRelation *r = &run->rels[i];
        if (!(r->continuity & ASM_CONT_INCUMBENT)) continue;
        if (changed && !changed[r->a] && !changed[r->b]) continue;
        double rms;
        if (!act_incumbent_metric(run,r,&rms) || rms > r->continuity_rms+1e-4) return 0;
    }
    return 1;
}

int AsmContinuity_preserves_incumbent(const AsmRun *run)
{
    return act_preserves_incumbent(run,NULL);
}

int AsmContinuity_preserves_incumbent_trial(const AsmRun *run, const AsmChart *saved, size_t n)
{
    /* A transaction holds every chart outside saved fixed. An untouched
     * seam has exactly its previously audited geometry, so only incident
     * seams need their trim frames evaluated again. Include both endpoints
     * even when just one side of the seam belongs to this transaction. */
    Arena_Mark mark = Arena_save(run->arena);
    uint8_t *changed = ARENA_CALLOC(run->arena,run->n_charts ? run->n_charts : 1,1);
    for (size_t i = 0; i < n; i++) {
        if (saved[i].id < 0 || (size_t)saved[i].id >= run->n_charts) {
            Arena_restore(run->arena,mark); return 0;
        }
        changed[saved[i].id] = 1;
    }
    int ok = act_preserves_incumbent(run,changed);
    Arena_restore(run->arena,mark); return ok;
}

#include "asm_source_reconstruction.inc"

int AsmContinuity_prepare(AsmRun *run)
{
    if (!run || !run->arena) return -1;
    /* Frozen physical frames and selected obligations belong to this run.
     * Resume/re-entry must not rebuild them from a later UV state. */
    if (run->continuity_ready) return 0;
    Arena_T arena = run->arena;
    Arena_Mark mark = Arena_save(arena);
    size_t nc = run->n_charts;
    size_t trim_measured = 0, trim_changed = 0, trim_gained = 0, trim_lost = 0;
    double phase_started = ves_clock_sec(), frame_seconds = 0, order_seconds = 0, incumbent_seconds = 0;
    double trim_error = 0, trim_max = 0;
    size_t *off = ARENA_CALLOC(arena, nc+1, sizeof(size_t));
    for (size_t r = 0; r < run->n_rels; r++) {
        AsmRelation *R = &run->rels[r]; R->continuity = 0; R->continuity_rms = 0;
        if (R->corr_first < 0 || R->corr_count <= 0) continue;
        if ((size_t)R->corr_first + (size_t)R->corr_count > run->n_corr) { Arena_restore(arena, mark); return -1; }
        off[R->a+1]++; off[R->b+1]++;
        for (int32_t k = 0; k < R->corr_count; k++) {
            AsmCorr *c = &run->corr[(size_t)R->corr_first+k]; c->valid = 0; c->run = -1;
            memset(&c->trim_a,0,sizeof c->trim_a); memset(&c->trim_b,0,sizeof c->trim_b);
        }
    }
    for (size_t i = 0; i < nc; i++) off[i+1] += off[i];
    size_t *adj = ARENA_ALLOC(arena, (off[nc] ? off[nc] : 1)*sizeof(size_t));
    size_t *fill = ARENA_CALLOC(arena, nc ? nc : 1, sizeof(size_t));
    for (size_t r = 0; r < run->n_rels; r++) {
        const AsmRelation *R = &run->rels[r];
        if (R->corr_first < 0 || R->corr_count <= 0) continue;
        adj[off[R->a]+fill[R->a]++] = r; adj[off[R->b]+fill[R->b]++] = r;
    }
    for (size_t i = 0; i < nc; i++) {
        const AsmChart *C = &run->charts[i];
        if (off[i] == off[i+1] || !AsmChart_in_layout(C)) continue;
        Arena_Mark cm = Arena_save(arena);
        uint8_t *valid = ARENA_CALLOC(arena,C->nv ? C->nv : 1,1);
        ActFrame *frame = act_vertex_frames(arena,C,valid);
        AsmTrim *trim = act_trim_bases(arena,C,frame,valid);
        if (!trim) { Arena_restore(arena,mark); return -1; }
        /* Record the old first-face measurement for diagnosis only. It
         * cannot affect the new frames, witnesses, or acceptance gates. */
        int32_t *first_face = ARENA_ALLOC(arena,C->nv*sizeof *first_face);
        for (size_t v = 0; v < C->nv; v++) first_face[v] = -1;
        for (size_t f = 0; f < C->nf; f++) for (int k = 0; k < 3; k++) {
            int32_t v = C->faces[3*f+k];
            if (v >= 0 && (size_t)v < C->nv && first_face[v] < 0) first_face[v] = (int32_t)f;
        }
        for (size_t j = off[i]; j < off[i+1]; j++) {
            const AsmRelation *R = &run->rels[adj[j]];
            const AsmChart *A = &run->charts[R->a], *B = &run->charts[R->b];
            for (int32_t k = 0; k < R->corr_count; k++) {
                AsmCorr *c = &run->corr[(size_t)R->corr_first+k];
                if (c->va < 0 || c->vb < 0 || (size_t)c->va >= A->nv || (size_t)c->vb >= B->nv) continue;
                if (!A->boundary || !B->boundary || !A->boundary[c->va] || !B->boundary[c->vb]) continue;
                double d[3];
                for (int q = 0; q < 3; q++) d[q] = B->xyz[(size_t)c->vb*3+q] - A->xyz[(size_t)c->va*3+q];
                int is_a = (int32_t)i == R->a;
                int32_t v = is_a ? c->va : c->vb;
                float legacy[2], *gap = is_a ? c->gap_a : c->gap_b;
                int old_ok = act_gap(C,first_face[v],d,legacy);
                int new_ok = valid[v] && act_frame_gap(&frame[v],d,gap);
                if (new_ok) {
                    if (!trim[v].valid) {
                        fprintf(stderr,"[assemble continuity] no stable original trim transport for chart %zu vertex %d\n",i,v);
                        Arena_restore(arena,mark); return -1;
                    }
                    if (is_a) c->trim_a = trim[v]; else c->trim_b = trim[v];
                    c->valid |= is_a ? 1u : 2u;
                }
                trim_gained += new_ok && !old_ok; trim_lost += old_ok && !new_ok;
                if (old_ok && new_ok) {
                    double error = hypot(gap[0]-legacy[0],gap[1]-legacy[1]);
                    trim_measured++; trim_changed += error > 0.25;
                    trim_error += error*error; trim_max = fmax(trim_max,error);
                }
            }
        }
        Arena_restore(arena, cm);
    }
    frame_seconds = ves_clock_sec() - phase_started; phase_started = ves_clock_sec();
    size_t nruns = 0, short_run_samples = 0, duplicate_samples = 0;
    for (size_t r = 0; r < run->n_rels; r++) {
        AsmRelation *R = &run->rels[r];
        if (R->corr_first < 0 || R->corr_count < 6) continue;
        const AsmChart *A = &run->charts[R->a], *B = &run->charts[R->b];
        if (!AsmChart_in_layout(A) || !AsmChart_in_layout(B)) continue;
        AsmBoundaryStats ordered = {0};
        if (AsmBoundary_order(arena,A,B,run->corr+(size_t)R->corr_first,(size_t)R->corr_count,&ordered)) {
            fprintf(stderr,"[assemble continuity] invalid source boundary topology or ambiguous witnesses on charts %d/%d\n",R->a,R->b);
            Arena_restore(arena,mark); return -1;
        }
        if (ordered.supported_runs) {
            R->continuity |= ASM_CONT_SOURCE;
            nruns += ordered.supported_runs; short_run_samples += ordered.short_run_samples;
        }
        duplicate_samples += ordered.duplicate_samples;
    }
    order_seconds = ves_clock_sec() - phase_started; phase_started = ves_clock_sec();
    size_t incumbent = 0;
    for (size_t i = 0; i < run->n_rels; i++) {
        AsmRelation *r = &run->rels[i]; double rms;
        if (!(r->continuity & ASM_CONT_SOURCE) && AsmRel_is_join(r) &&
            act_incumbent_metric(run,r,&rms)) {
            r->continuity |= ASM_CONT_INCUMBENT; r->continuity_rms = rms; incumbent++;
        }
    }
    fprintf(stderr,"[assemble continuity] %zu sparse clean2 seams retained as provisional placement bundles; source certificates remain unresolved\n",incumbent);
    fprintf(stderr, "[assemble continuity] %zu ordered source boundary runs; trim offsets measured from incident triangles\n", nruns);
    fprintf(stderr,"[assemble continuity] source-edge ordering retains %zu short-run observations; %zu duplicate observations cannot add votes\n",short_run_samples,duplicate_samples);
    fprintf(stderr,"[assemble continuity] trim frames versus first face: %zu common correspondence endpoints, %zu newly valid, %zu lost; %zu change > 0.25 vox, RMS %.6g, max %.6g vox\n",
            trim_measured,trim_gained,trim_lost,trim_changed,trim_measured ? sqrt(trim_error/trim_measured) : 0,trim_max);
    Arena_restore(arena, mark);
    incumbent_seconds = ves_clock_sec() - phase_started; phase_started = ves_clock_sec();
    int rc = AsmContinuity_reconstruct(run);
    fprintf(stderr, "[assemble continuity] seconds: trim frames %.1f, boundary order %.1f, incumbent seams %.1f, source reconstruction %.1f\n",
            frame_seconds, order_seconds, incumbent_seconds, ves_clock_sec() - phase_started);
    if (!rc) run->continuity_ready = 1;
    return rc;
}

int AsmContinuity_target(const AsmRun *run, const AsmRelation *r, const AsmCorr *c, int reverse, double target[2])
{
    if (c->valid != 3 || c->run < 0 || (r->flags & (ASM_REL_CONTACT | ASM_REL_DROPPED))) return 0;
    if (!run->charts || r->a < 0 || r->b < 0 || (size_t)r->a >= run->n_charts || (size_t)r->b >= run->n_charts) return 0;
    if (reverse) return act_trim_point(&run->charts[r->b],c->vb,&c->trim_b,c->gap_b,-1,target);
    return act_trim_point(&run->charts[r->a],c->va,&c->trim_a,c->gap_a,1,target);
}

double AsmContinuity_pair_error(const AsmRun *run, const AsmRelation *r, const AsmCorr *c)
{
    if (c->valid != 3 || c->run < 0 || r->a < 0 || r->b < 0 ||
        (size_t)r->a >= run->n_charts || (size_t)r->b >= run->n_charts) return 1e300;
    const AsmChart *A = &run->charts[r->a], *B = &run->charts[r->b];
    if (!A->uv || !B->uv || c->va < 0 || c->vb < 0 || (size_t)c->va >= A->nv || (size_t)c->vb >= B->nv) return 1e300;
    double a[2], b[2], ta[2], tb[2];
    act_point(A,c->va,NULL,a); act_point(B,c->vb,NULL,b);
    if (!AsmContinuity_target(run,r,c,0,tb) || !AsmContinuity_target(run,r,c,1,ta)) return 1e300;
    double ea = hypot(a[0]-ta[0],a[1]-ta[1]), eb = hypot(b[0]-tb[0],b[1]-tb[1]);
    return isfinite(ea) && isfinite(eb) ? fmax(ea,eb) : 1e300;
}

int AsmContinuity_measure(const AsmRun *run, const AsmRelation *r, double *rms, size_t *support)
{
    *rms = 1e300; *support = 0;
    if (!(r->continuity & ASM_CONT_SOURCE) || (r->flags & (ASM_REL_CONTACT | ASM_REL_DROPPED))) return 0;
    const AsmChart *A = &run->charts[r->a], *B = &run->charts[r->b];
    if (((A->flags ^ B->flags) & ASM_CHART_MIRROR) != ((r->flags & ASM_REL_PARITY) ? ASM_CHART_MIRROR : 0)) return 0;
    double e2 = 0; size_t good = 0;
    double gate = ASM_CONTINUITY_GATE(r->rms);
    for (int32_t k = 0; k < r->corr_count; k++) {
        const AsmCorr *c = &run->corr[(size_t)r->corr_first+k];
        if (c->valid != 3 || c->run < 0) continue;
        double e = AsmContinuity_pair_error(run,r,c); e2 += e*e; (*support)++;
        if (e <= gate) good++;
    }
    if (*support) *rms = sqrt(e2 / *support);
    return *support >= 6 && good*10 >= *support*9 && *rms <= gate;
}

typedef AsmContinuityHalfspace ActHalfspace;

/* Bound inverse images by bytes as well as row count. The working set's
 * triangular matrix uses at most half that many entries plus its diagonal:
 * its dimension cannot exceed either the variable count or the row count. */
static size_t act_halfspace_capacity(size_t variables)
{
    if (!variables) return 0;
    size_t capacity = (1u << 24)/variables;
    return capacity < 8192 ? capacity : 8192;
}

/* A feasible primal working set. Blocking rows join the equality solve;
 * rows with negative multipliers leave it. Normalize in the H metric and
 * reuse the existing inverse images, so only the binding Schur system is
 * factored. Failure leaves the caller's unconstrained step untouched. */
static int act_halfspace_active(Arena_T arena, int n, int nt, const int *rows, const int *cols,
                                const double *vals, const ActHalfspace *constraint, size_t m,
                                const double *images, const double *diagonal, double *x, double *lambda,
                                int limit)
{
    if (limit <= 0) return 0;
    for (size_t i = 0; i < m; i++) if (constraint[i].rhs > 0) return 0;
    Arena_Mark mark = Arena_save(arena);
    size_t cap = m < (size_t)n ? m : (size_t)n, triangles = cap*(cap+1)/2;
    size_t *active = ARENA_ALLOC(arena,cap*sizeof *active), count = 0;
    uint8_t *member = ARENA_CALLOC(arena,m,1);
    double *current = ARENA_CALLOC(arena,(size_t)n,sizeof *current);
    double *target = ARENA_ALLOC(arena,(size_t)n*sizeof *target), *direction = ARENA_ALLOC(arena,(size_t)n*sizeof *direction);
    double *scale = ARENA_CALLOC(arena,(size_t)n,sizeof *scale), *norm = ARENA_ALLOC(arena,m*sizeof *norm);
    int *sr = ARENA_ALLOC(arena,triangles*sizeof *sr), *sc = ARENA_ALLOC(arena,triangles*sizeof *sc);
    double *sv = ARENA_ALLOC(arena,triangles*sizeof *sv), *b = ARENA_ALLOC(arena,cap*sizeof *b), *mu = ARENA_ALLOC(arena,cap*sizeof *mu);
    double *correction = ARENA_ALLOC(arena,cap*sizeof *correction);
    for (int j = 0; j < nt; j++) if (rows[j] == cols[j]) scale[rows[j]] += vals[j];
    for (int j = 0; j < n; j++) scale[j] = sqrt(scale[j]);
    for (size_t i = 0; i < m; i++) norm[i] = sqrt(diagonal[i]);
    int ok = 0, iteration = 0; const char *failure = "iteration limit";
    for (; iteration < limit; iteration++) {
        memcpy(target,x,(size_t)n*sizeof *target);
        if (count) {
            size_t triplets = 0;
            for (size_t u = 0; u < count; u++) {
                size_t iu = active[u]; const ActHalfspace *a = &constraint[iu]; double ax = 0;
                for (int k = 0; k < a->count; k++) ax += a->values[k]*(current[a->ids[k]]-x[a->ids[k]]);
                b[u] = ax/norm[iu];
                for (size_t v = 0; v <= u; v++) {
                    size_t iv = active[v]; const ActHalfspace *c = &constraint[iv]; double uv = 0, vu = 0;
                    for (int k = 0; k < a->count; k++) uv += a->values[k]*images[(size_t)n*iv+a->ids[k]];
                    for (int k = 0; k < c->count; k++) vu += c->values[k]*images[(size_t)n*iu+c->ids[k]];
                    sr[triplets] = (int)u; sc[triplets] = (int)v;
                    sv[triplets++] = .5*(uv+vu)/(norm[iu]*norm[iv]);
                }
            }
            memset(mu,0,count*sizeof *mu);
            SparseFactor_T schur = NULL;
            if (Sparse_factor_spd((int)count,(int)triplets,sr,sc,sv,&schur) != 0 || Sparse_factor_solve(schur,b,mu) != 0) {
                Sparse_factor_free(&schur); failure = "dependent working set"; break;
            }
            for (size_t u = 0; u < count; u++) {
                double weight = mu[u]/norm[active[u]];
                for (int j = 0; j < n; j++) target[j] += weight*images[(size_t)n*active[u]+j];
            }
            /* Keep the direction tangent to the actual feasible iterate,
             * rather than correcting small accumulated equality slack in
             * the direction itself. Refine its residual against the same
             * factor before a nearly dependent inactive row can block it. */
            int refined = 1;
            for (int pass = 0; pass < 2; pass++) {
                for (size_t u = 0; u < count; u++) {
                    const ActHalfspace *a = &constraint[active[u]]; double ap = 0;
                    for (int k = 0; k < a->count; k++) ap += a->values[k]*(target[a->ids[k]]-current[a->ids[k]]);
                    b[u] = -ap/norm[active[u]];
                }
                if (Sparse_factor_solve(schur,b,correction) != 0) { refined = 0; break; }
                for (size_t u = 0; u < count; u++) {
                    mu[u] += correction[u]; double weight = correction[u]/norm[active[u]];
                    for (int j = 0; j < n; j++) target[j] += weight*images[(size_t)n*active[u]+j];
                }
            }
            Sparse_factor_free(&schur);
            if (!refined) { failure = "equality refinement"; break; }
        }
        double length2 = 0, current2 = 0;
        for (int j = 0; j < n; j++) {
            direction[j] = target[j]-current[j];
            double v = scale[j]*direction[j], w = scale[j]*current[j]; length2 += v*v; current2 += w*w;
        }
        if (!isfinite(length2) || !isfinite(current2)) { failure = "nonfinite direction"; break; }
        if (sqrt(length2) <= 1e-10*(1+sqrt(current2))) {
            size_t drop = SIZE_MAX; double minimum = -1e-8;
            for (size_t u = 0; u < count; u++) if (mu[u] < minimum) { minimum = mu[u]; drop = u; }
            if (drop != SIZE_MAX) {
                member[active[drop]] = 0;
                memmove(active+drop,active+drop+1,(count-drop-1)*sizeof *active); count--; continue;
            }
            ok = 1; break;
        }
        double alpha = 1, threshold = 64*DBL_EPSILON*(1+sqrt(length2)); size_t blocking = SIZE_MAX;
        for (size_t i = 0; i < m; i++) if (!member[i]) {
            const ActHalfspace *c = &constraint[i]; double ax = 0, velocity = 0;
            for (int k = 0; k < c->count; k++) { ax += c->values[k]*current[c->ids[k]]; velocity += c->values[k]*direction[c->ids[k]]; }
            velocity /= norm[i];
            if (velocity < -threshold) {
                double ratio = fmax(0,(ax-c->rhs)/norm[i])/-velocity;
                if (ratio < alpha) { alpha = ratio; blocking = i; }
            }
        }
        if (alpha == 1) memcpy(current,target,(size_t)n*sizeof *current);
        else for (int j = 0; j < n; j++) current[j] += alpha*direction[j];
        int feasible = 1;
        for (size_t i = 0; i < m; i++) {
            const ActHalfspace *c = &constraint[i]; double ax = 0;
            for (int k = 0; k < c->count; k++) ax += c->values[k]*current[c->ids[k]];
            if (!isfinite(ax) || c->rhs-ax > 1e-8*(1+fabs(c->rhs))) { feasible = 0; break; }
        }
        if (!feasible) { failure = "primal residual"; break; }
        if (blocking != SIZE_MAX) {
            if (count == cap) { failure = "working set dimension"; break; }
            active[count++] = blocking; member[blocking] = 1;
        }
    }
    if (ok) {
        memcpy(x,current,(size_t)n*sizeof *x); memset(lambda,0,m*sizeof *lambda);
        for (size_t u = 0; u < count; u++) lambda[active[u]] = fmax(0,mu[u]/norm[active[u]]);
        fprintf(stderr,"  [continuity halfspace] active set: %zu of %zu constraints, %d iterations\n",count,m,iteration+1);
    } else fprintf(stderr,"  [continuity halfspace] active set deferred: %s, %zu of %zu constraints, %d iterations\n",failure,count,m,iteration+1);
    Arena_restore(arena,mark); return ok;
}

/* Minimize .5*x'Hx-b'x subject to A*x >= d. Projected dual coordinate
 * updates reuse one SPD factor. Return 1 for a converged solution or 2 for
 * a feasible improving step obtained by projecting onto the half-spaces.
 * Neither outcome certifies geometry: the caller still runs full audits.
 * Infeasible rows cannot license a step. Inverse images use at most 128 MiB. */
static int act_halfspace_solve(Arena_T arena, int n, int nt, const int *rows, const int *cols,
                               const double *vals, const double *rhs, double *x,
                               const ActHalfspace *constraint, size_t m, int max_sweeps, int active_limit)
{
    if (n <= 0 || m > act_halfspace_capacity((size_t)n) || max_sweeps < 1 || max_sweeps > 4096 || active_limit < 0 || active_limit > 2048) return 0;
    if (!m) return Sparse_solve_sym(n,nt,rows,cols,vals,rhs,x,SPARSE_SPD) == 0;
    Arena_Mark mark = Arena_save(arena);
    double *images = ARENA_ALLOC(arena,(size_t)n*m*sizeof *images);
    double *work = ARENA_CALLOC(arena,(size_t)n,sizeof *work);
    double *diagonal = ARENA_ALLOC(arena,m*sizeof *diagonal), *lambda = ARENA_CALLOC(arena,m,sizeof *lambda);
    SparseFactor_T factor = NULL; int ok = 0;
    if (Sparse_factor_spd(n,nt,rows,cols,vals,&factor) != 0 || Sparse_factor_solve(factor,rhs,x) != 0) goto done;
    for (size_t i = 0; i < m; i++) {
        const ActHalfspace *c = &constraint[i];
        if (c->count < 1 || c->count > 6 || !isfinite(c->rhs)) goto done;
        memset(work,0,(size_t)n*sizeof *work);
        for (int k = 0; k < c->count; k++) {
            if (c->ids[k] < 0 || c->ids[k] >= n || !isfinite(c->values[k])) goto done;
            work[c->ids[k]] += c->values[k];
        }
        double *z = images+(size_t)n*i;
        if (Sparse_factor_solve(factor,work,z) != 0) goto done;
        double d = 0;
        for (int k = 0; k < c->count; k++) d += c->values[k]*z[c->ids[k]];
        if (!(d > 0) || !isfinite(d)) goto done;
        diagonal[i] = d;
    }
    ok = act_halfspace_active(arena,n,nt,rows,cols,vals,constraint,m,images,diagonal,x,lambda,active_limit);
    for (int sweep = 0; sweep < max_sweeps && !ok; sweep++) {
        for (size_t i = 0; i < m; i++) {
            const ActHalfspace *c = &constraint[i]; double ax = 0;
            for (int k = 0; k < c->count; k++) ax += c->values[k]*x[c->ids[k]];
            double next = fmax(0,lambda[i]+(c->rhs-ax)/diagonal[i]);
            if (!isfinite(next) || !isfinite(ax)) goto done;
            double delta = next-lambda[i]; lambda[i] = next;
            if (delta) for (int j = 0; j < n; j++) x[j] += delta*images[(size_t)n*i+j];
        }
        ok = 1;
        for (size_t i = 0; i < m; i++) {
            const ActHalfspace *c = &constraint[i]; double ax = 0;
            for (int k = 0; k < c->count; k++) ax += c->values[k]*x[c->ids[k]];
            double error = c->rhs-ax;
            if (!isfinite(error) || (lambda[i] > 0 ? fabs(error) : error) > 1e-8*(1+fabs(c->rhs))) { ok = 0; break; }
        }
    }
    if (!ok) {
        double primal = 0, complementarity = 0;
        for (size_t i = 0; i < m; i++) {
            const ActHalfspace *c = &constraint[i]; double ax = 0;
            for (int k = 0; k < c->count; k++) ax += c->values[k]*x[c->ids[k]];
            double error = c->rhs-ax;
            primal = fmax(primal,error);
            complementarity = fmax(complementarity,lambda[i] > 0 ? fabs(error) : fmax(0,error));
        }
        fprintf(stderr,"  [continuity halfspace] %d variables, %zu constraints: %d sweeps without convergence; primal %.9g vox, complementarity %.9g vox\n",
                n,m,max_sweeps,primal,complementarity);
        /* Stop optimizing the dual and first seek primal feasibility.
         * Cyclic H-metric projections only correct violated half-spaces;
         * they do not remove earlier multipliers to approach the optimum.
         * Zero must be a known feasible step, so these sets intersect.
         * A tiny cushion stays inside that known slack. */
        int zero_feasible = 1, projected = 0;
        for (size_t i = 0; i < m; i++) if (constraint[i].rhs > 0) zero_feasible = 0;
        if (zero_feasible) for (int sweep = 0; sweep < 512 && !projected; sweep++) {
            for (size_t i = 0; i < m; i++) {
                const ActHalfspace *c = &constraint[i]; double ax = 0;
                for (int k = 0; k < c->count; k++) ax += c->values[k]*x[c->ids[k]];
                if (!isfinite(ax)) goto done;
                double delta = fmax(0,(c->rhs+fmin(1e-7,-.25*c->rhs)-ax)/diagonal[i]);
                if (!isfinite(delta)) goto done;
                if (delta) for (int j = 0; j < n; j++) x[j] += delta*images[(size_t)n*i+j];
            }
            projected = 1;
            for (size_t i = 0; i < m; i++) {
                const ActHalfspace *c = &constraint[i]; double ax = 0;
                for (int k = 0; k < c->count; k++) ax += c->values[k]*x[c->ids[k]];
                if (!isfinite(ax) || ax < c->rhs) { projected = 0; break; }
            }
        }
        /* If projection remains incomplete, the segment toward zero still
         * gives a computable feasible endpoint. Any fallback must improve
         * the actual quadratic, then pass the full nonlinear pose audits. */
        double alpha = 1;
        for (size_t i = 0; i < m; i++) {
            const ActHalfspace *c = &constraint[i]; double ax = 0;
            if (c->rhs > 0) { alpha = 0; break; }
            for (int k = 0; k < c->count; k++) ax += c->values[k]*x[c->ids[k]];
            if (!isfinite(ax)) { alpha = 0; break; }
            if (ax < c->rhs) alpha = fmin(alpha,.999*c->rhs/ax);
        }
        if (alpha > 1e-8) {
            for (int j = 0; j < n; j++) x[j] *= alpha;
            double objective = 0, scale = 1;
            for (int j = 0; j < n; j++) { double v = rhs[j]*x[j]; objective -= v; scale += fabs(v); }
            for (int t = 0; t < nt; t++) {
                double v = (rows[t] == cols[t] ? .5 : 1)*vals[t]*x[rows[t]]*x[cols[t]];
                objective += v; scale += fabs(v);
            }
            ok = isfinite(objective) && isfinite(scale) && objective < -1e-10*scale ? 2 : 0;
            for (size_t i = 0; i < m && ok; i++) {
                const ActHalfspace *c = &constraint[i]; double ax = 0;
                for (int k = 0; k < c->count; k++) ax += c->values[k]*x[c->ids[k]];
                if (!isfinite(ax) || ax < c->rhs) ok = 0;
            }
            if (ok) fprintf(stderr,"  [continuity halfspace] feasible %s %.9g, quadratic improvement %.9g; full pose audits still required\n",
                            projected ? "projection" : "retreat",alpha,-objective);
        }
    }
    if (ok == 1) {
        /* Check stationarity independently of the accumulated inverse
         * images, including the original sparse matrix's duplicate rows. */
        double *scale = ARENA_ALLOC(arena,(size_t)n*sizeof *scale);
        for (int j = 0; j < n; j++) { work[j] = -rhs[j]; scale[j] = 1+fabs(rhs[j]); }
        for (int t = 0; t < nt; t++) {
            int a = rows[t], b = cols[t]; double v = vals[t]*x[b]; work[a] += v; scale[a] += fabs(v);
            if (a != b) { v = vals[t]*x[a]; work[b] += v; scale[b] += fabs(v); }
        }
        for (size_t i = 0; i < m; i++) for (int k = 0; k < constraint[i].count; k++) {
            int j = constraint[i].ids[k]; double v = constraint[i].values[k]*lambda[i]; work[j] -= v; scale[j] += fabs(v);
        }
        for (int j = 0; j < n; j++) if (!isfinite(x[j]) || !isfinite(work[j]) || fabs(work[j]) > 1e-7*scale[j]) ok = 0;
    }
done:
    Sparse_factor_free(&factor); Arena_restore(arena,mark); return ok;
}

int AsmContinuity_project_step(Arena_T arena, int n, int nt,
    const int *rows, const int *cols, const double *values, const double *rhs,
    double *step, const AsmContinuityHalfspace *constraints, size_t count)
{
    return act_halfspace_solve(arena,n,nt,rows,cols,values,rhs,step,
                              constraints,count,1024,256);
}

typedef struct ActPosePlane {
    AsmContinuityUvPair pair;
    double nx, ny, distance;
} ActPosePlane;

static void act_pose_plane_row(const AsmRun *run, const int32_t *index, const ActPosePlane *plane, ActHalfspace *row)
{
    const AsmChart *a = &run->charts[plane->pair.chart_a], *b = &run->charts[plane->pair.chart_b];
    double pa[2], pb[2]; act_point(a,plane->pair.vertex_a,NULL,pa); act_point(b,plane->pair.vertex_b,NULL,pb);
    double separation = plane->nx*(pb[0]-pa[0])+plane->ny*(pb[1]-pa[1]);
    row->count = 0; row->rhs = plane->distance-separation;
    int ia = index[a->id], ib = index[b->id];
    if (ia >= 0) {
        int k = row->count; row->ids[k] = 3*ia; row->values[k++] = -plane->nx;
        row->ids[k] = 3*ia+1; row->values[k++] = -plane->ny;
        row->ids[k] = 3*ia+2; row->values[k++] = plane->nx*(pa[1]-a->pose_y)-plane->ny*(pa[0]-a->pose_x); row->count = k;
    }
    if (ib >= 0) {
        int k = row->count; row->ids[k] = 3*ib; row->values[k++] = plane->nx;
        row->ids[k] = 3*ib+1; row->values[k++] = plane->ny;
        row->ids[k] = 3*ib+2; row->values[k++] = -plane->nx*(pb[1]-b->pose_y)+plane->ny*(pb[0]-b->pose_x); row->count = k;
    }
}

static int act_pose_planes_step(const AsmRun *run, const int32_t *index, const ActPosePlane *planes,
                                size_t n, const double *step, double alpha)
{
    for (size_t j = 0; j < n; j++) {
        const ActPosePlane *p = &planes[j];
        AsmChart a = run->charts[p->pair.chart_a], b = run->charts[p->pair.chart_b];
        int ia = index[a.id], ib = index[b.id];
        if (ia >= 0) { a.pose_x += alpha*step[3*ia]; a.pose_y += alpha*step[3*ia+1]; a.pose_theta += alpha*step[3*ia+2]; }
        if (ib >= 0) { b.pose_x += alpha*step[3*ib]; b.pose_y += alpha*step[3*ib+1]; b.pose_theta += alpha*step[3*ib+2]; }
        double pa[2], pb[2]; act_point(&a,p->pair.vertex_a,NULL,pa); act_point(&b,p->pair.vertex_b,NULL,pb);
        double separation = p->nx*(pb[0]-pa[0])+p->ny*(pb[1]-pa[1]);
        if (!isfinite(separation) || separation < p->distance) return 0;
    }
    return 1;
}

typedef struct ActPoseQuadratic {
    int n, nt;
    const int *rows, *cols;
    const double *vals, *rhs;
} ActPoseQuadratic;

/* Hold the proposed rotations fixed and correct only translations. With
 * fixed rotations these inequalities are exactly linear, so cyclic
 * projections can restore feasibility without discarding a tangent
 * rotation merely because its second-order curvature crosses a plane.
 * This is only a proposal: verify all nonlinear planes before exposing it,
 * and leave the input untouched on failure. Scratch is O(variables+planes). */
static int act_pose_planes_project(const AsmRun *run, const int32_t *index, const ActPosePlane *planes,
                                   size_t n_planes, int variables, double *step, double *maximum)
{
    Arena_T arena = run->arena; Arena_Mark mark = Arena_save(arena);
    ActHalfspace *rows = ARENA_CALLOC(arena,n_planes,sizeof *rows);
    double *diagonal = ARENA_CALLOC(arena,n_planes,sizeof *diagonal);
    double *delta = ARENA_CALLOC(arena,(size_t)variables,sizeof *delta);
    double *candidate = ARENA_ALLOC(arena,(size_t)variables*sizeof *candidate);
    int ok = 0;
    for (size_t i = 0; i < n_planes; i++) {
        const ActPosePlane *p = &planes[i]; ActHalfspace *row = &rows[i];
        AsmChart a = run->charts[p->pair.chart_a], b = run->charts[p->pair.chart_b];
        int ia = index[a.id], ib = index[b.id];
        if (ia >= 0) { a.pose_x += step[3*ia]; a.pose_y += step[3*ia+1]; a.pose_theta += step[3*ia+2]; }
        if (ib >= 0) { b.pose_x += step[3*ib]; b.pose_y += step[3*ib+1]; b.pose_theta += step[3*ib+2]; }
        double pa[2], pb[2]; act_point(&a,p->pair.vertex_a,NULL,pa); act_point(&b,p->pair.vertex_b,NULL,pb);
        double separation = p->nx*(pb[0]-pa[0])+p->ny*(pb[1]-pa[1]);
        if (!isfinite(separation)) goto done;
        if (ia == ib) { if (separation < p->distance) goto done; continue; }
        row->rhs = p->distance-separation+1e-8*(1+fabs(p->distance));
        if (ia >= 0) {
            row->ids[row->count] = 3*ia; row->values[row->count++] = -p->nx;
            row->ids[row->count] = 3*ia+1; row->values[row->count++] = -p->ny;
        }
        if (ib >= 0) {
            row->ids[row->count] = 3*ib; row->values[row->count++] = p->nx;
            row->ids[row->count] = 3*ib+1; row->values[row->count++] = p->ny;
        }
        for (int k = 0; k < row->count; k++) diagonal[i] += row->values[k]*row->values[k];
        if (!(diagonal[i] > 0) || !isfinite(diagonal[i]) || !isfinite(row->rhs)) goto done;
    }
    for (int sweep = 0; sweep < 128; sweep++) {
        for (size_t i = 0; i < n_planes; i++) if (rows[i].count) {
            const ActHalfspace *row = &rows[i]; double value = 0;
            for (int k = 0; k < row->count; k++) value += row->values[k]*delta[row->ids[k]];
            if (!isfinite(value)) goto done;
            if (value < row->rhs) {
                double amount = (row->rhs-value)/diagonal[i];
                for (int k = 0; k < row->count; k++) delta[row->ids[k]] += amount*row->values[k];
            }
        }
        int feasible = 1;
        for (size_t i = 0; i < n_planes; i++) {
            const ActHalfspace *row = &rows[i]; double value = 0;
            for (int k = 0; k < row->count; k++) value += row->values[k]*delta[row->ids[k]];
            /* Half the projection cushion still leaves room for the
             * independent evaluation's different floating-point order. */
            if (!isfinite(value) || value < row->rhs-5e-9*(1+fabs(planes[i].distance))) { feasible = 0; break; }
        }
        if (!feasible) continue;
        for (int j = 0; j < variables; j++) candidate[j] = step[j]+delta[j];
        if (!act_pose_planes_step(run,index,planes,n_planes,candidate,1)) goto done;
        double largest = 0;
        for (int j = 0; j < variables; j += 3) largest = fmax(largest,hypot(delta[j],delta[j+1]));
        if (!isfinite(largest)) goto done;
        *maximum = largest; memcpy(step,candidate,(size_t)variables*sizeof *step); ok = 1; break;
    }
done:
    Arena_restore(arena,mark); return ok;
}

static int act_pose_planes_direction(const AsmRun *run, const int32_t *index, const ActPosePlane *planes,
                                     size_t n_planes, const ActPoseQuadratic *q, double *step)
{
    Arena_T arena = run->arena; Arena_Mark mark = Arena_save(arena);
    double *candidate = ARENA_ALLOC(arena,(size_t)q->n*sizeof *candidate); int ok = 0;
    for (double alpha = 1; alpha >= 1.0/4096; alpha *= .5) {
        for (int j = 0; j < q->n; j++) candidate[j] = alpha*step[j];
        double correction = 0;
        if (!act_pose_planes_step(run,index,planes,n_planes,candidate,1) &&
            !act_pose_planes_project(run,index,planes,n_planes,q->n,candidate,&correction)) continue;
        /* A geometric correction must still improve this exact quadratic;
         * the caller additionally certifies the actual SOURCE residuals,
         * original pose bounds, and every original good geometry pair. */
        double objective = 0, scale = 1;
        for (int j = 0; j < q->n; j++) { double v = q->rhs[j]*candidate[j]; objective -= v; scale += fabs(v); }
        for (int t = 0; t < q->nt; t++) {
            double v = (q->rows[t] == q->cols[t] ? .5 : 1)*q->vals[t]*candidate[q->rows[t]]*candidate[q->cols[t]];
            objective += v; scale += fabs(v);
        }
        if (!isfinite(objective) || !isfinite(scale) || objective >= -1e-10*scale) continue;
        memcpy(step,candidate,(size_t)q->n*sizeof *step); ok = 1;
        if (correction || alpha < 1) fprintf(stderr,"  [continuity halfspace] nonlinear feasible step: alpha %.9g, maximum translation correction %.9g, quadratic improvement %.9g\n",
                                             alpha,correction,-objective);
        break;
    }
    Arena_restore(arena,mark); return ok;
}

static int act_close(AsmRun *run, const AsmChart *saved, size_t n, const uint8_t *selected,
                     const AsmContinuityTarget *targets, size_t n_targets, int boundary_fit, int restrict_selected,
                     const ActPosePlane *planes, size_t n_planes, const uint8_t *held)
{
    if (!n || n > 100000 || (restrict_selected && !selected)) return 0;
    Arena_T arena = run->arena; Arena_Mark mark = Arena_save(arena);
    int32_t *index = ARENA_ALLOC(arena, run->n_charts*sizeof(int32_t));
    for (size_t i = 0; i < run->n_charts; i++) index[i] = -1;
    for (size_t i = 0; i < n; i++) index[saved[i].id] = (int32_t)i;
    size_t count = 0, fixed = 0;
    for (size_t r = 0; r < run->n_rels; r++) {
        const AsmRelation *R = &run->rels[r];
        if (!(R->continuity & ASM_CONT_SOURCE) || (R->flags & (ASM_REL_CONTACT | ASM_REL_DROPPED))) continue;
        if (restrict_selected && !selected[r]) continue;
        if (selected && !AsmRel_is_join(R) && !(R->continuity & ASM_CONT_CUT) && !selected[r]) continue;
        int a = index[R->a] >= 0, b = index[R->b] >= 0;
        if (!a && !b) continue;
        if ((!a && run->charts[R->a].placement_state == ASM_PLACE_NONE && (!held || !held[R->a])) ||
            (!b && run->charts[R->b].placement_state == ASM_PLACE_NONE && (!held || !held[R->b]))) continue;
        if (a != b) fixed++;
        count += (size_t)R->corr_count;
    }
    if ((!fixed && n_targets < 6) || n_targets > 100000000 || count > (2147483647u-3*n-6*n_targets)/42) {
        fprintf(stderr, "  [continuity solve] unanchored/oversize: %zu charts, %zu fixed seams, %zu targets, %zu pairs\n", n, fixed, n_targets, count);
        Arena_restore(arena, mark); return 0;
    }
    size_t cap = 42*count+3*n+6*n_targets;
    int *rows = ARENA_ALLOC(arena, cap*sizeof(int)), *cols = ARENA_ALLOC(arena, cap*sizeof(int));
    double *vals = ARENA_ALLOC(arena, cap*sizeof(double));
    double *rhs = ARENA_ALLOC(arena, 3*n*sizeof(double)), *step = ARENA_ALLOC(arena, 3*n*sizeof(double));
    double *start = ARENA_ALLOC(arena, 3*n*sizeof(double));
    ActHalfspace *constraints = n_planes ? ARENA_ALLOC(arena,n_planes*sizeof *constraints) : NULL;
    for (size_t i = 0; i < n; i++) {
        const AsmChart *c = &run->charts[saved[i].id];
        start[3*i] = c->pose_x; start[3*i+1] = c->pose_y; start[3*i+2] = c->pose_theta;
    }
    int ok = 1;
    for (int it = 0; it < (boundary_fit ? 16 : 8) && ok; it++) {
        memset(rhs, 0, 3*n*sizeof(double)); memset(step, 0, 3*n*sizeof(double));
        size_t nt = 0;
        for (size_t r = 0; r < run->n_rels; r++) {
            const AsmRelation *R = &run->rels[r];
            if (!(R->continuity & ASM_CONT_SOURCE) || (R->flags & (ASM_REL_CONTACT | ASM_REL_DROPPED))) continue;
            if (restrict_selected && !selected[r]) continue;
            if (selected && !AsmRel_is_join(R) && !(R->continuity & ASM_CONT_CUT) && !selected[r]) continue;
            int32_t ia = index[R->a], ib = index[R->b];
            const AsmChart *A = &run->charts[R->a], *B = &run->charts[R->b];
            if (ia < 0 && ib < 0) continue;
            if ((ia < 0 && A->placement_state == ASM_PLACE_NONE && (!held || !held[R->a])) ||
                (ib < 0 && B->placement_state == ASM_PLACE_NONE && (!held || !held[R->b]))) continue;
            if (((A->flags ^ B->flags) & ASM_CHART_MIRROR) != ((R->flags & ASM_REL_PARITY) ? ASM_CHART_MIRROR : 0)) {
                fprintf(stderr, "  [continuity solve] parity: relation %d-%d, iteration %d\n", R->a, R->b, it);
                ok = 0; break;
            }
            for (int32_t k = 0; k < R->corr_count; k++) {
                const AsmCorr *c = &run->corr[(size_t)R->corr_first+k];
                if (c->valid != 3 || c->run < 0) continue;
                /* Both tangent frames define the measured trim. Optimizing
                 * only A's prediction makes the answer depend on which
                 * chart happened to be stored first, while certification
                 * checks both directions. Keep the total pair weight. */
                for (int direction = 0; direction < 2; direction++) {
                double a[2], b[2];
                act_point(A,c->va,NULL,a); act_point(B,c->vb,NULL,b);
                if (!AsmContinuity_target(run,R,c,direction,direction ? b : a)) { ok = 0; break; }
                double ex = b[0]-a[0], ey = b[1]-a[1];
                double w = 0.5 / fmax(1.0, hypot(ex,ey)/4.0);
                double hxx = w, hxy = 0, hyy = w;
                if (boundary_fit) {
                    /* The robust fit may sacrifice the boundary tail even
                     * when a rigid pose can satisfy every measured point.
                     * Minimize squared distance outside the certification
                     * balls, with a small interior least-squares tie break.
                     * Use the radial Hessian, not IRLS weights that alternate
                     * between opposite violated ends of a seam. */
                    double e = hypot(ex,ey), gate = 0.95*ASM_CONTINUITY_GATE(R->rms);
                    double outside = e > gate ? 1.0-gate/e : 0;
                    w = 0.5*(0.001+outside);
                    hxx = hyy = w; hxy = 0;
                    if (e > gate) {
                        double nx = ex/e, ny = ey/e, radial = 0.5*gate/e;
                        hxx += radial*nx*nx; hxy = radial*nx*ny; hyy += radial*ny*ny;
                    }
                }
                int ids[6], m = 0; double jx[6], jy[6];
                if (ia >= 0) {
                    ids[m] = 3*ia; jx[m] = -1; jy[m++] = 0;
                    ids[m] = 3*ia+1; jx[m] = 0; jy[m++] = -1;
                    ids[m] = 3*ia+2; jx[m] = a[1]-A->pose_y; jy[m++] = -(a[0]-A->pose_x);
                }
                if (ib >= 0) {
                    ids[m] = 3*ib; jx[m] = 1; jy[m++] = 0;
                    ids[m] = 3*ib+1; jx[m] = 0; jy[m++] = 1;
                    ids[m] = 3*ib+2; jx[m] = -(b[1]-B->pose_y); jy[m++] = b[0]-B->pose_x;
                }
                for (int u = 0; u < m; u++) {
                    rhs[ids[u]] -= w*(jx[u]*ex+jy[u]*ey);
                    for (int v = 0; v <= u; v++) {
                        rows[nt] = ids[u] > ids[v] ? ids[u] : ids[v]; cols[nt] = ids[u] < ids[v] ? ids[u] : ids[v];
                        vals[nt++] = boundary_fit ? hxx*jx[u]*jx[v]+hxy*(jx[u]*jy[v]+jy[u]*jx[v])+hyy*jy[u]*jy[v]
                                                 : w*(jx[u]*jx[v]+jy[u]*jy[v]);
                    }
                }
                }
            }
        }
        for (size_t t = 0; t < n_targets; t++) {
            const AsmContinuityTarget *target = &targets[t];
            int32_t i = index[target->chart]; if (i < 0 || !(target->weight > 0)) continue;
            const AsmChart *c = &run->charts[target->chart];
            double u = c->flags & ASM_CHART_MIRROR ? -target->u : target->u;
            double v = target->v, ct = cos(c->pose_theta), sn = sin(c->pose_theta);
            double px = ct*u-sn*v, py = sn*u+ct*v;
            double ex = target->x-px-c->pose_x, ey = target->y-py-c->pose_y;
            double w = target->weight / fmax(1.0, hypot(ex,ey)/16.0);
            double jx[3] = {1,0,-py}, jy[3] = {0,1,px};
            for (int a = 0; a < 3; a++) {
                rhs[3*i+a] += w*(jx[a]*ex+jy[a]*ey);
                for (int b = 0; b <= a; b++) { rows[nt] = 3*i+a; cols[nt] = 3*i+b; vals[nt++] = w*(jx[a]*jx[b]+jy[a]*jy[b]); }
            }
        }
        /* Only a numerical tether for underobserved modes, never a source
         * of placement support. Data residuals decide acceptance below. */
        for (size_t i = 0; i < 3*n; i++) { rows[nt] = cols[nt] = (int)i; vals[nt++] = i%3 == 2 ? 0.01 : 1e-6; }
        if (!ok) break;
        for (size_t j = 0; j < n_planes; j++) act_pose_plane_row(run,index,&planes[j],&constraints[j]);
        int solved = act_halfspace_solve(arena,(int)(3*n),(int)nt,rows,cols,vals,rhs,step,constraints,n_planes,4096,2048);
        if (!solved) {
            fprintf(stderr, "  [continuity solve] solver failed: %zu charts, %zu fixed seams, %zu targets, %zu pairs, %zu half-planes, iteration %d\n", n, fixed, n_targets, count, n_planes, it);
            ok = 0; break;
        }
        if (n_planes) {
            ActPoseQuadratic quadratic = {(int)(3*n),(int)nt,rows,cols,vals,rhs};
            if (!act_pose_planes_direction(run,index,planes,n_planes,&quadratic,step)) { ok = 0; break; }
        }
        /* A full Newton step can overshoot even when the rigid solution is
         * inside the original trust region. Shorten the shared direction
         * before applying any chart; never relax the per-step/total bounds
         * or a nonlinear contact plane to admit that shorter step. */
        double alpha=1;
        for (size_t i=0; i<3*n; i++) if (!isfinite(step[i])) { ok=0; break; }
        if (!ok) break;
        for (; alpha>=1.0/4096; alpha*=0.5) {
            int bounded=1;
            for (size_t i=0; i<n; i++) {
                const AsmChart *c=&run->charts[saved[i].id];
                double dx=alpha*step[3*i],dy=alpha*step[3*i+1],dt=alpha*step[3*i+2];
                if (fabs(dt)>0.25 || hypot(dx,dy)>256.0 ||
                    hypot(c->pose_x+dx-start[3*i],c->pose_y+dy-start[3*i+1])>256.0 ||
                    fabs(c->pose_theta+dt-start[3*i+2])>0.25) { bounded=0; break; }
            }
            if (bounded && (!n_planes || act_pose_planes_step(run,index,planes,n_planes,step,alpha))) break;
        }
        if (alpha<1.0/4096) {
            fprintf(stderr,"  [continuity solve] no bounded feasible step: iteration %d (%zu charts, %zu targets)\n",it,n,n_targets);
            ok=0; break;
        }
        if (alpha<1) fprintf(stderr,"  [continuity solve] bounded step: alpha %.9g, iteration %d (%zu charts, %zu targets)\n",alpha,it,n,n_targets);
        double largest = 0;
        for (size_t i = 0; i < n; i++) {
            AsmChart *c = &run->charts[saved[i].id];
            double dx = alpha*step[3*i], dy = alpha*step[3*i+1], dt = alpha*step[3*i+2];
            c->pose_x += dx; c->pose_y += dy; c->pose_theta += dt;
            double m = hypot(dx,dy)+128*fabs(dt); if (m > largest) largest = m;
        }
        if (largest < 1e-5) break;
    }
    Arena_restore(arena, mark); return ok;
}

int AsmContinuity_close(AsmRun *run, const AsmChart *saved, size_t n, const uint8_t *selected,
                        const AsmContinuityTarget *targets, size_t n_targets)
{
    return act_close(run,saved,n,selected,targets,n_targets,0,0,NULL,0,NULL);
}

int AsmContinuity_close_boundaries(AsmRun *run, const AsmChart *saved, size_t n, const uint8_t *selected)
{
    return act_close(run,saved,n,selected,NULL,0,1,0,NULL,0,NULL);
}

/* A supporting half-plane of the forbidden UV-distance disk contains the
 * original pair and implies the original isometry test. Pairs outside the
 * 24-voxel rim may instead stay outside it. Constraints use the unchanged
 * checkpoint, including the query chart's own layer tolerance. */
static size_t act_pose_planes_add(const AsmRun *run, const AsmChart *const *checkpoint,
                                  ActPosePlane *planes, size_t *count, size_t capacity,
                                  const AsmContinuityUvPair *pairs, size_t n_pairs)
{
    size_t before = *count;
    for (size_t j = 0; j < n_pairs && *count < capacity; j++) {
        AsmContinuityUvPair pair = pairs[j];
        if (pair.chart_a < 0 || (size_t)pair.chart_a >= run->n_charts ||
            pair.chart_b < 0 || (size_t)pair.chart_b >= run->n_charts) return SIZE_MAX;
        const AsmChart *a = checkpoint[pair.chart_a], *b = checkpoint[pair.chart_b];
        if (!a || !b || !a->uv || !b->uv || !a->xyz || !b->xyz || pair.vertex_a < 0 || pair.vertex_b < 0 ||
            (size_t)pair.vertex_a >= a->nv || (size_t)pair.vertex_b >= b->nv) return SIZE_MAX;
        if (pair.chart_a == pair.chart_b) continue; /* rigid motion cannot change this distance */
        int duplicate = 0;
        for (size_t k = 0; k < *count; k++) if (!memcmp(&pair,&planes[k].pair,sizeof pair)) { duplicate = 1; break; }
        if (duplicate) continue;
        double pa[2], pb[2]; act_point(a,pair.vertex_a,NULL,pa); act_point(b,pair.vertex_b,NULL,pb);
        double dx = pb[0]-pa[0], dy = pb[1]-pa[1], distance = hypot(dx,dy);
        const float *p = a->xyz+3*(size_t)pair.vertex_a, *q = b->xyz+3*(size_t)pair.vertex_b;
        double physical = sqrt((q[0]-p[0])*(q[0]-p[0])+(q[1]-p[1])*(q[1]-p[1])+(q[2]-p[2])*(q[2]-p[2]));
        double dl = run->chart_layer_d ? run->chart_layer_d[pair.chart_a] : 0;
        if (dl <= 0) dl = run->layer_d_global;
        double tol = fmax(.5*dl,7.0), required = fmin(physical-tol,24.0);
        if (!(distance > 0) || !isfinite(distance) || !isfinite(physical) ||
            (distance <= 24 && physical-distance > tol)) return SIZE_MAX;
        if (required <= 0) continue;
        double margin = fmin(1e-5,.25*fmax(0,distance-required));
        planes[(*count)++] = (ActPosePlane){pair,dx/distance,dy/distance,required+margin};
    }
    return *count-before;
}

/* A kept pose-graph path may put both sides of a measured, unkept seam in
 * the same component. Try that nearby SOURCE evidence before placement,
 * when rigid chart poses can still move. The current registration remains
 * the incumbent; no relation flag or source measurement is rewritten. */
static int act_register_hypotheses(AsmRun *run, const AsmChart *entry, size_t n, AsmContinuityUvAudit audit)
{
    if (n < 2 || n > 100001 || !audit) return 0;
    Arena_T arena = run->arena; Arena_Mark mark = Arena_save(arena);
    uint8_t *member = ARENA_CALLOC(arena,run->n_charts,1);
    AsmChart *saved = ARENA_ALLOC(arena,n*sizeof *saved);
    for (size_t i = 0; i < n; i++) { saved[i] = run->charts[entry[i].id]; member[entry[i].id] = 1; }
    uint8_t *selected = ARENA_CALLOC(arena,run->n_rels ? run->n_rels : 1,1);
    uint8_t *hypothesis = ARENA_CALLOC(arena,run->n_rels ? run->n_rels : 1,1);
    uint8_t *protected = ARENA_CALLOC(arena,run->n_rels ? run->n_rels : 1,1);
    uint8_t *measured = ARENA_CALLOC(arena,run->n_rels ? run->n_rels : 1,1);
    size_t nh = 0; double before = 0;
    for (size_t j = 0; j < run->n_rels; j++) {
        const AsmRelation *r = &run->rels[j];
        if ((!member[r->a] && !member[r->b]) || !AsmChart_in_layout(&run->charts[r->a]) || !AsmChart_in_layout(&run->charts[r->b])) continue;
        double rms; size_t support; int passes = AsmContinuity_measure(run,r,&rms,&support);
        protected[j] = (uint8_t)passes;
        if (!member[r->a] || !member[r->b] || support < 6 || !isfinite(rms)) continue;
        measured[j] = 1; before += rms*rms*support;
        if (rms <= 2*ASM_CONTINUITY_GATE(r->rms)) {
            selected[j] = 1;
            if (!passes) { hypothesis[j] = 1; nh++; }
        }
    }
    if (!nh || !isfinite(before)) { Arena_restore(arena,mark); return 0; }
    /* A failed seam is a local repair request. Keep the rest of the
     * component exact: a globally improved registration can otherwise
     * perturb distant placement anchors and lose material downstream.
     * Move the failed endpoints and one ring of measured SOURCE neighbours;
     * the original largest chart remains fixed even inside that ring. */
    uint8_t *seed = ARENA_CALLOC(arena,run->n_charts,1), *moving = ARENA_CALLOC(arena,run->n_charts,1);
    uint8_t *held = ARENA_CALLOC(arena,run->n_charts,1);
    for (size_t j = 0; j < run->n_rels; j++) if (hypothesis[j]) seed[run->rels[j].a] = seed[run->rels[j].b] = 1;
    memcpy(moving,seed,run->n_charts);
    for (size_t j = 0; j < run->n_rels; j++) if (selected[j] && (seed[run->rels[j].a] || seed[run->rels[j].b]))
        moving[run->rels[j].a] = moving[run->rels[j].b] = 1;
    moving[saved[n-1].id] = 0;
    AsmChart *variables = ARENA_ALLOC(arena,n*sizeof *variables); size_t n_variables = 0;
    for (size_t i = 0; i < n; i++) {
        if (moving[saved[i].id]) variables[n_variables++] = saved[i];
        else held[saved[i].id] = 1;
    }
    if (!n_variables) { Arena_restore(arena,mark); return 0; }
    fprintf(stderr,"  [registration neighbourhood] component %d: %zu moving charts, %zu held at their incumbent poses\n",saved[0].component,n_variables,n-n_variables);
    /* The final saved chart is the fixed anchor; the caller temporarily
     * keeps its ROOT state. Every retry begins at the same incumbent. */
    const AsmChart **checkpoint = ARENA_CALLOC(arena,run->n_charts,sizeof *checkpoint);
    for (size_t i = 0; i < n; i++) checkpoint[saved[i].id] = &saved[i];
    size_t plane_capacity = act_halfspace_capacity(3*n_variables);
    size_t feedback_capacity = plane_capacity < 512 ? plane_capacity : 512;
    ActPosePlane *planes = ARENA_ALLOC(arena,plane_capacity*sizeof *planes); size_t n_planes = 0;
    AsmContinuityUvPair *feedback = ARENA_ALLOC(arena,feedback_capacity*sizeof *feedback);
    AsmContinuityAuditCache cache = {0}; cache.scratch = arena;
    cache.guard_charts = saved; cache.guard_count = n;
    cache.pairs = feedback; cache.pairs_capacity = feedback_capacity;
    int ok = 0;
    size_t repaired = 0; double after = 0;
    for (int attempt = 0; attempt < 32; attempt++) {
    for (size_t i = 0; i < n; i++) run->charts[saved[i].id] = saved[i];
    int solved = act_close(run,variables,n_variables,selected,NULL,0,0,1,planes,n_planes,held);
    /* A later numerical or step-limit failure can leave an earlier useful
     * pose in place. Solver completion is not certification. Independently
     * validate the current pose, including every original bound, before
     * either accepting it or collecting further geometry constraints. */
    ok = 1;
    repaired = 0; after = 0;
    for (size_t j = 0; j < run->n_rels; j++) if (protected[j] || measured[j]) {
        double rms; size_t support; int passes = AsmContinuity_measure(run,&run->rels[j],&rms,&support);
        if (protected[j] && !passes) ok = 0;
        if (measured[j]) {
            if (support < 6 || !isfinite(rms)) ok = 0;
            else after += rms*rms*support;
        }
        if (hypothesis[j] && passes) repaired++;
    }
    ok = ok && repaired && isfinite(after) && after < before-fmax(1e-8,before*1e-6);
    for (size_t i = 0; i < n && ok; i++) {
        const AsmChart *c = &run->charts[entry[i].id];
        if (!isfinite(c->pose_x) || !isfinite(c->pose_y) || !isfinite(c->pose_theta) ||
            hypot(c->pose_x-entry[i].pose_x,c->pose_y-entry[i].pose_y) > 256 ||
            fabs(c->pose_theta-entry[i].pose_theta) > .25) ok = 0;
    }
    if (!ok) break;
    if (!solved) fprintf(stderr,"  [registration hypotheses] component %d: solver stopped; checking the current improved pose against full geometry\n",saved[0].component);
    ok = audit(run,saved,n,&cache);
    if (ok) break;
    size_t added = act_pose_planes_add(run,checkpoint,planes,&n_planes,plane_capacity,feedback,cache.pairs_count);
    if (!added || added == SIZE_MAX) break;
    fprintf(stderr,"  [registration constraints] component %d: %zu original pair half-planes, retry %d\n",saved[0].component,n_planes,attempt+1);
    }
    if (!ok) for (size_t i = 0; i < n; i++) run->charts[saved[i].id] = saved[i];
    fprintf(stderr,"  [registration hypotheses] component %d, %zu charts: %zu nearby failed seams, %zu repaired in trial, source energy %.9g -> %.9g; %s\n",
            saved[0].component,n,nh,repaired,before,after,ok ? "accepted" : "incumbent retained");
    Arena_restore(arena,mark); return ok;
}

size_t AsmContinuity_register_components(AsmRun *run, AsmContinuityUvAudit audit)
{
    if (!audit) return 0;
    Arena_T arena = run->arena; Arena_Mark mark = Arena_save(arena);
    size_t nc = run->n_charts, ng = 0;
    for (size_t i = 0; i < nc; i++) if (AsmChart_in_layout(&run->charts[i]) && run->charts[i].component >= 0 && (size_t)run->charts[i].component >= ng) ng = (size_t)run->charts[i].component+1;
    size_t *off = ARENA_CALLOC(arena, ng+1, sizeof(size_t)), *fill = ARENA_CALLOC(arena, ng ? ng : 1, sizeof(size_t));
    for (size_t i = 0; i < nc; i++) if (AsmChart_in_layout(&run->charts[i]) && run->charts[i].component >= 0) off[run->charts[i].component+1]++;
    for (size_t i = 0; i < ng; i++) off[i+1] += off[i];
    int32_t *ids = ARENA_ALLOC(arena, (nc ? nc : 1)*sizeof(int32_t));
    for (size_t i = 0; i < nc; i++) if (AsmChart_in_layout(&run->charts[i]) && run->charts[i].component >= 0) {
        int32_t g = run->charts[i].component; ids[off[g]+fill[g]++] = (int32_t)i;
    }
    uint8_t *selected = ARENA_CALLOC(arena, run->n_rels ? run->n_rels : 1, 1);
    size_t kept = 0, refused = 0;
    for (size_t g = 0; g < ng; g++) {
        size_t n = off[g+1]-off[g]; if (n < 2) continue;
        Arena_Mark gm = Arena_save(arena);
        int32_t fixed = ids[off[g]];
        for (size_t k = off[g]+1; k < off[g+1]; k++) if (run->charts[ids[k]].area3d > run->charts[fixed].area3d) fixed = ids[k];
        AsmChart *saved = ARENA_ALLOC(arena, n*sizeof(AsmChart)); size_t ns = 0;
        for (size_t k = off[g]; k < off[g+1]; k++) if (ids[k] != fixed) saved[ns++] = run->charts[ids[k]];
        saved[ns] = run->charts[fixed];
        /* The least-squares registration can sacrifice a sparse seam to
         * improve denser ones. As in UV repair, every incident source seam
         * that already passes must survive, including seams outside this
         * component and hypotheses not used by the registration solve. */
        size_t *protected = ARENA_ALLOC(arena,(run->n_rels ? run->n_rels : 1)*sizeof *protected), np = 0;
        for (size_t r = 0; r < run->n_rels; r++) {
            const AsmRelation *R = &run->rels[r];
            const AsmChart *a = &run->charts[R->a], *b = &run->charts[R->b];
            if (!AsmChart_in_layout(a) || !AsmChart_in_layout(b) ||
                (a->component != (int32_t)g && b->component != (int32_t)g)) continue;
            double rms; size_t support;
            if (AsmContinuity_measure(run,R,&rms,&support)) protected[np++] = r;
        }
        int state = run->charts[fixed].placement_state; run->charts[fixed].placement_state = ASM_PLACE_ROOT;
        AsmContinuityAuditCache audit_cache = {0}; audit_cache.scratch = arena;
        int ok = AsmContinuity_close(run,saved,ns,selected,NULL,0);
        for (size_t j = 0; j < np && ok; j++) {
            const AsmRelation *r = &run->rels[protected[j]]; double rms; size_t support;
            if (!AsmContinuity_measure(run,r,&rms,&support)) {
                fprintf(stderr,"  [registration source] component %zu breaks previously passing seam %d-%d (%.6g vox, %zu pairs); registration reverted\n",
                        g,r->a,r->b,rms,support); ok = 0;
            }
        }
        if (ok && audit(run,saved,n,&audit_cache)) kept++;
        else { for (size_t k = 0; k < n; k++) run->charts[saved[k].id] = saved[k]; refused++; }
        run->charts[fixed].placement_state = state;
        Arena_restore(arena, gm);
    }
    fprintf(stderr, "[assemble continuity] source registration: %zu components solved, %zu trials reverted; one fixed chart per component\n", kept, refused);
    Arena_restore(arena, mark); return kept;
}

size_t AsmContinuity_register_failed(AsmRun *run, const AsmChart *entry, AsmContinuityUvAudit audit)
{
    if (!entry || !audit) return 0;
    Arena_T arena = run->arena; Arena_Mark mark = Arena_save(arena);
    size_t nc = run->n_charts, ng = 0, kept = 0;
    for (size_t i = 0; i < nc; i++) if (AsmChart_in_layout(&run->charts[i]) && run->charts[i].component >= 0 && (size_t)run->charts[i].component >= ng)
        ng = (size_t)run->charts[i].component+1;
    size_t *off = ARENA_CALLOC(arena,ng+1,sizeof *off), *fill = ARENA_CALLOC(arena,ng ? ng : 1,sizeof *fill);
    for (size_t i = 0; i < nc; i++) if (AsmChart_in_layout(&run->charts[i]) && run->charts[i].component >= 0) off[run->charts[i].component+1]++;
    for (size_t i = 0; i < ng; i++) off[i+1] += off[i];
    int32_t *ids = ARENA_ALLOC(arena,(nc ? nc : 1)*sizeof *ids);
    for (size_t i = 0; i < nc; i++) if (AsmChart_in_layout(&run->charts[i]) && run->charts[i].component >= 0) {
        int32_t g = run->charts[i].component; ids[off[g]+fill[g]++] = (int32_t)i;
    }
    for (size_t g = 0; g < ng; g++) {
        size_t n = off[g+1]-off[g]; if (n < 2) continue;
        Arena_Mark gm = Arena_save(arena);
        int32_t fixed = ids[off[g]];
        for (size_t k = off[g]+1; k < off[g+1]; k++) if (run->charts[ids[k]].area3d > run->charts[fixed].area3d) fixed = ids[k];
        AsmChart *original = ARENA_ALLOC(arena,n*sizeof *original); size_t ns = 0;
        for (size_t k = off[g]; k < off[g+1]; k++) if (ids[k] != fixed) original[ns++] = entry[ids[k]];
        original[ns] = entry[fixed];
        int state = run->charts[fixed].placement_state; run->charts[fixed].placement_state = ASM_PLACE_ROOT;
        kept += act_register_hypotheses(run,original,n,audit);
        run->charts[fixed].placement_state = state;
        Arena_restore(arena,gm);
    }
    fprintf(stderr,"[assemble continuity] post-UV source registration: %zu components improved; original pose bounds and passing seams preserved\n",kept);
    Arena_restore(arena,mark); return kept;
}

/* Keep the correction local in physical distance, independent of mesh
 * resolution. The heap holds each vertex at most once. */
static void act_heap_up(int32_t *heap, int32_t *pos, const double *dist, size_t at)
{
    int32_t v = heap[at];
    while (at) {
        size_t parent = (at-1)/2;
        if (dist[heap[parent]] <= dist[v]) break;
        heap[at] = heap[parent]; pos[heap[at]] = (int32_t)at; at = parent;
    }
    heap[at] = v; pos[v] = (int32_t)at;
}

static int32_t act_heap_pop(int32_t *heap, int32_t *pos, const double *dist, size_t *n)
{
    int32_t out = heap[0], v = heap[--*n]; size_t at = 0;
    pos[out] = -2;
    if (!*n) return out;
    while (2*at+1 < *n) {
        size_t child = 2*at+1;
        if (child+1 < *n && dist[heap[child+1]] < dist[heap[child]]) child++;
        if (dist[v] <= dist[heap[child]]) break;
        heap[at] = heap[child]; pos[heap[at]] = (int32_t)at; at = child;
    }
    heap[at] = v; pos[v] = (int32_t)at; return out;
}

static void act_uv_operator(const ActSpring *spring, size_t ns, const ActMatch *match, size_t nm,
                            const ActVec *x, ActVec *y, size_t n)
{
    /* Numerical tether fixes the translation gauge of a small chart whose
     * whole surface is in the band. It contributes no continuity evidence. */
    for (size_t i = 0; i < n; i++) y[i] = (ActVec){1e-4*x[i].x, 1e-4*x[i].y};
    for (size_t i = 0; i < ns; i++) {
        int32_t a = spring[i].a, b = spring[i].b;
        ActVec d = {0,0};
        if (a >= 0) d = x[a];
        if (b >= 0) { d.x -= x[b].x; d.y -= x[b].y; }
        d.x *= spring[i].weight; d.y *= spring[i].weight;
        if (a >= 0) { y[a].x += d.x; y[a].y += d.y; }
        if (b >= 0) { y[b].x -= d.x; y[b].y -= d.y; }
    }
    for (size_t i = 0; i < nm; i++) {
        ActVec d = {0,0};
        for (int k = 0; k < match[i].count; k++) if (match[i].id[k] >= 0) {
            int32_t v = match[i].id[k]; double a = match[i].coefficient[k];
            d.x += a*x[v].x; d.y += a*x[v].y;
        }
        for (int k = 0; k < match[i].count; k++) if (match[i].id[k] >= 0) {
            int32_t v = match[i].id[k]; double a = match[i].weight*match[i].coefficient[k];
            y[v].x += a*d.x; y[v].y += a*d.y;
        }
    }
}

typedef struct ActIcEntry { int32_t col; double value; } ActIcEntry;
typedef struct ActIc { size_t *off; ActIcEntry *entry; double *diagonal; } ActIc;

static int act_ic_compare(const void *a, const void *b)
{
    int32_t x = ((const ActIcEntry *)a)->col, y = ((const ActIcEntry *)b)->col;
    return (x > y)-(x < y);
}

/* Incomplete Cholesky with no fill. Affine trim rows can introduce positive
 * off-diagonal entries, so incomplete pivots are not guaranteed positive.
 * The caller retains its positive diagonal preconditioner if this fails.
 * Storage is bounded by the actual mesh edges and four-vertex row products. */
static int act_ic_build(Arena_T arena, const ActSpring *spring, size_t ns,
                        const ActMatch *match, size_t nm, const double *diag, size_t n, ActIc *ic)
{
    ic->off = ARENA_CALLOC(arena,n+1,sizeof *ic->off);
    ic->entry = ARENA_ALLOC(arena,(ns+6*nm ? ns+6*nm : 1)*sizeof *ic->entry);
    ic->diagonal = ARENA_ALLOC(arena,n*sizeof *ic->diagonal);
    size_t *fill = ARENA_ALLOC(arena,n*sizeof *fill);
    for (int pass = 0; pass < 2; pass++) {
        for (size_t i = 0; i < ns; i++) {
            int32_t a = spring[i].a, b = spring[i].b;
            double w = spring[i].weight;
            if (a < 0 || b < 0 || a == b || w == 0) continue;
            if (a < b) { int32_t t = a; a = b; b = t; }
            if (!pass) ic->off[a+1]++;
            else ic->entry[fill[a]++] = (ActIcEntry){b,-w};
        }
        for (size_t i = 0; i < nm; i++) for (int k = 0; k < match[i].count; k++)
            for (int j = 0; j < k; j++) {
                int32_t a = match[i].id[k], b = match[i].id[j];
                double w = match[i].weight*match[i].coefficient[k]*match[i].coefficient[j];
                if (a < 0 || b < 0 || a == b || w == 0) continue;
                if (a < b) { int32_t v = a; a = b; b = v; }
                if (!pass) ic->off[a+1]++;
                else ic->entry[fill[a]++] = (ActIcEntry){b,w};
            }
        if (!pass) {
            for (size_t i = 0; i < n; i++) ic->off[i+1] += ic->off[i];
            memcpy(fill,ic->off,n*sizeof *fill);
        }
    }
    size_t used = 0;
    for (size_t i = 0; i < n; i++) {
        size_t begin = ic->off[i], end = ic->off[i+1];
        qsort(ic->entry+begin,end-begin,sizeof *ic->entry,act_ic_compare);
        ic->off[i] = used;
        for (size_t j = begin; j < end; j++) {
            ActIcEntry e = ic->entry[j];
            if (used > ic->off[i] && ic->entry[used-1].col == e.col) ic->entry[used-1].value += e.value;
            else ic->entry[used++] = e;
        }
    }
    ic->off[n] = used;
    for (size_t i = 0; i < n; i++) {
        double d = diag[i];
        for (size_t j = ic->off[i]; j < ic->off[i+1]; j++) {
            ActIcEntry *e = &ic->entry[j]; int32_t col = e->col;
            double v = e->value;
            size_t a = ic->off[i], b = ic->off[col];
            while (a < j && b < ic->off[col+1]) {
                int32_t ca = ic->entry[a].col, cb = ic->entry[b].col;
                if (ca == cb) { v -= ic->entry[a].value*ic->entry[b].value; a++; b++; }
                else if (ca < cb) a++; else b++;
            }
            v /= ic->diagonal[col]; e->value = v; d -= v*v;
        }
        if (!(d > 0) || !isfinite(d)) return 0;
        ic->diagonal[i] = sqrt(d);
    }
    return 1;
}

static double act_uv_precondition(const ActIc *ic, const double *diag,
                                  const ActVec *r, ActVec *z, size_t n)
{
    if (!ic) for (size_t i = 0; i < n; i++) z[i] = (ActVec){r[i].x/diag[i],r[i].y/diag[i]};
    else {
        for (size_t i = 0; i < n; i++) {
            ActVec v = r[i];
            for (size_t j = ic->off[i]; j < ic->off[i+1]; j++) {
                const ActIcEntry *e = &ic->entry[j];
                v.x -= e->value*z[e->col].x; v.y -= e->value*z[e->col].y;
            }
            z[i] = (ActVec){v.x/ic->diagonal[i],v.y/ic->diagonal[i]};
        }
        for (size_t i = n; i-- > 0;) {
            z[i].x /= ic->diagonal[i]; z[i].y /= ic->diagonal[i];
            for (size_t j = ic->off[i]; j < ic->off[i+1]; j++) {
                const ActIcEntry *e = &ic->entry[j];
                z[e->col].x -= e->value*z[i].x; z[e->col].y -= e->value*z[i].y;
            }
        }
    }
    double rz = 0;
    for (size_t i = 0; i < n; i++) rz += r[i].x*z[i].x+r[i].y*z[i].y;
    return rz;
}

static void act_uv_terms(const ActSpring *spring, size_t ns, const ActMatch *match, size_t nm,
                          double *diag, ActVec *rhs, size_t n)
{
    for (size_t i = 0; i < n; i++) { diag[i] = 1e-4; rhs[i] = (ActVec){0,0}; }
    for (size_t i = 0; i < ns; i++) {
        if (spring[i].a >= 0) diag[spring[i].a] += spring[i].weight;
        if (spring[i].b >= 0) diag[spring[i].b] += spring[i].weight;
    }
    for (size_t i = 0; i < nm; i++) {
        for (int k = 0; k < match[i].count; k++) if (match[i].id[k] >= 0) {
            int32_t v = match[i].id[k]; double a = match[i].coefficient[k], w = match[i].weight;
            diag[v] += w*a*a; rhs[v].x += w*a*match[i].target.x; rhs[v].y += w*a*match[i].target.y;
        }
    }
}

static int act_uv_solve(Arena_T arena, const ActSpring *spring, size_t ns,
                        const ActMatch *match, size_t nm, ActVec *x, size_t n, int *iterations)
{
    ActVec *r = ARENA_ALLOC(arena,n*sizeof *r), *p = ARENA_ALLOC(arena,n*sizeof *p);
    ActVec *ap = ARENA_ALLOC(arena,n*sizeof *ap);
    double *diag = ARENA_ALLOC(arena,n*sizeof *diag);
    act_uv_terms(spring,ns,match,nm,diag,r,n);
    double initial = 0;
    for (size_t i = 0; i < n; i++) initial += (r[i].x*r[i].x+r[i].y*r[i].y)/diag[i];
    *iterations = 0;
    if (!isfinite(initial)) return 0;
    if (initial < 1e-20) return 1;
    ActIc factor;
    const ActIc *ic = act_ic_build(arena,spring,ns,match,nm,diag,n,&factor) ? &factor : NULL;
    ActVec *z = ARENA_ALLOC(arena,n*sizeof *z);
    double rz = act_uv_precondition(ic,diag,r,p,n);
    if (!(rz > 0) || !isfinite(rz)) return 0;
    for (int it = 0; it < 500; it++) {
        act_uv_operator(spring, ns, match, nm, p, ap, n);
        double pap = 0;
        for (size_t i = 0; i < n; i++) pap += p[i].x*ap[i].x + p[i].y*ap[i].y;
        if (!(pap > 0) || !isfinite(pap)) return 0;
        double alpha = rz/pap, error = 0;
        for (size_t i = 0; i < n; i++) {
            x[i].x += alpha*p[i].x; x[i].y += alpha*p[i].y;
            r[i].x -= alpha*ap[i].x; r[i].y -= alpha*ap[i].y;
            error += (r[i].x*r[i].x+r[i].y*r[i].y)/diag[i];
        }
        *iterations = it+1;
        if (!isfinite(error)) return 0;
        if (error <= initial*1e-12) {
            /* Judge the true residual in the original diagonal norm. An
             * approximate factor or recursive residual cannot relax the
             * original stopping tolerance. Restart if roundoff drifted. */
            act_uv_operator(spring,ns,match,nm,x,r,n);
            for (size_t i = 0; i < n; i++) { r[i].x = -r[i].x; r[i].y = -r[i].y; }
            for (size_t i = 0; i < nm; i++) {
                for (int k = 0; k < match[i].count; k++) if (match[i].id[k] >= 0) {
                    int32_t v = match[i].id[k]; double w = match[i].weight*match[i].coefficient[k];
                    r[v].x += w*match[i].target.x; r[v].y += w*match[i].target.y;
                }
            }
            error = 0;
            for (size_t i = 0; i < n; i++) error += (r[i].x*r[i].x+r[i].y*r[i].y)/diag[i];
            if (!isfinite(error)) return 0;
            if (error <= initial*1e-12) return 1;
            rz = act_uv_precondition(ic,diag,r,p,n);
            if (!(rz > 0) || !isfinite(rz)) return 0;
            continue;
        }
        double next = act_uv_precondition(ic,diag,r,z,n);
        if (!(next > 0) || !isfinite(next)) return 0;
        double beta = next/rz; rz = next;
        for (size_t i = 0; i < n; i++) {
            p[i].x = z[i].x+beta*p[i].x;
            p[i].y = z[i].y+beta*p[i].y;
        }
    }
    return 0;
}

static ActVec act_uv_project(ActVec x, ActVec offset)
{
    double u = x.x+offset.x, v = x.y+offset.y, d2 = u*u+v*v;
    if (d2 > 64.0) { double scale = 8.0/sqrt(d2); x = (ActVec){scale*u-offset.x,scale*v-offset.y}; }
    return x;
}

static double act_uv_quadratic(const ActVec *x, const ActVec *ax, const ActVec *rhs, size_t n)
{
    double energy = 0;
    for (size_t i = 0; i < n; i++) energy += x[i].x*(.5*ax[i].x-rhs[i].x)+x[i].y*(.5*ax[i].y-rhs[i].y);
    return energy;
}

/* The same convex quadratic, constrained to each vertex's total movement
 * disk. Twice the diagonal majorizes this positive-edge Laplacian, so a
 * diagonal proximal step is valid and its projection is a Euclidean disk.
 * Acceleration restarts whenever energy rises. Certify the true projected
 * gradient before returning a direction; a capped iteration count is a
 * deferred optional refinement, never a converged solve. */
static int act_uv_bounded_direction(Arena_T arena, const ActSpring *spring, size_t ns,
                                    const ActMatch *match, size_t nm, const ActVec *linear,
                                    const ActVec *offset, ActVec *x, size_t n, int *iterations)
{
    Arena_Mark mark = Arena_save(arena);
    double *diag = ARENA_ALLOC(arena,n*sizeof *diag);
    ActVec *rhs = ARENA_ALLOC(arena,n*sizeof *rhs), *y = ARENA_ALLOC(arena,n*sizeof *y);
    ActVec *next = ARENA_ALLOC(arena,n*sizeof *next), *ax = ARENA_ALLOC(arena,n*sizeof *ax);
    ActVec *work = ARENA_ALLOC(arena,n*sizeof *work);
    act_uv_terms(spring,ns,match,nm,diag,rhs,n);
    double initial = 0;
    for (size_t i = 0; i < n; i++) {
        initial += (rhs[i].x*rhs[i].x+rhs[i].y*rhs[i].y)/diag[i];
        y[i] = x[i] = act_uv_project(linear[i],offset[i]);
    }
    act_uv_operator(spring,ns,match,nm,x,ax,n);
    double energy = act_uv_quadratic(x,ax,rhs,n), t = 1;
    int ok = 0; *iterations = 0;
    if (!isfinite(initial) || !isfinite(energy)) { Arena_restore(arena,mark); return 0; }
    for (int it = 0; it <= 500; it++) {
        double error = 0;
        for (size_t i = 0; i < n; i++) {
            ActVec step = {x[i].x-(ax[i].x-rhs[i].x)/(2*diag[i]),x[i].y-(ax[i].y-rhs[i].y)/(2*diag[i])};
            step = act_uv_project(step,offset[i]);
            double du = step.x-x[i].x, dv = step.y-x[i].y;
            error += 2*diag[i]*(du*du+dv*dv);
        }
        if (isfinite(error) && error <= fmax(1e-20,initial*1e-12)) { ok = 1; break; }
        if (it == 500 || !isfinite(error)) break;
        double candidate = 0;
        for (int restart = 0; restart < 2; restart++) {
            if (restart) { memcpy(y,x,n*sizeof *x); t = 1; }
            act_uv_operator(spring,ns,match,nm,y,work,n);
            for (size_t i = 0; i < n; i++) {
                ActVec step = {y[i].x-(work[i].x-rhs[i].x)/(2*diag[i]),y[i].y-(work[i].y-rhs[i].y)/(2*diag[i])};
                next[i] = act_uv_project(step,offset[i]);
            }
            act_uv_operator(spring,ns,match,nm,next,work,n);
            candidate = act_uv_quadratic(next,work,rhs,n);
            if (isfinite(candidate) && candidate <= energy+fmax(1e-12,fabs(energy)*1e-14)) break;
        }
        if (!isfinite(candidate) || candidate > energy+fmax(1e-12,fabs(energy)*1e-14)) break;
        double nt = .5*(1+sqrt(1+4*t*t)), beta = (t-1)/nt; t = nt;
        for (size_t i = 0; i < n; i++) {
            y[i] = (ActVec){next[i].x+beta*(next[i].x-x[i].x),next[i].y+beta*(next[i].y-x[i].y)};
            x[i] = next[i]; ax[i] = work[i];
        }
        energy = candidate; *iterations = it+1;
    }
    Arena_restore(arena,mark); return ok;
}

/* Singular values of UV relative to the triangle's physical metric. */
static double act_uv_quality(const AsmChart *c, const int32_t *f, const float *uv, double *signed_area)
{
    double a = 0, b = 0, e = 0;
    for (int k = 0; k < 3; k++) {
        double x = c->xyz[(size_t)f[1]*3+k]-c->xyz[(size_t)f[0]*3+k];
        double y = c->xyz[(size_t)f[2]*3+k]-c->xyz[(size_t)f[0]*3+k];
        a += x*x; b += x*y; e += y*y;
    }
    double x = uv[(size_t)f[1]*2]-uv[(size_t)f[0]*2], y = uv[(size_t)f[1]*2+1]-uv[(size_t)f[0]*2+1];
    double u = uv[(size_t)f[2]*2]-uv[(size_t)f[0]*2], v = uv[(size_t)f[2]*2+1]-uv[(size_t)f[0]*2+1];
    *signed_area = x*v-y*u;
    double det = a*e-b*b;
    if (!(det > 1e-12*a*e)) return 1e300;
    double trace = ((x*x+y*y)*e+(u*u+v*v)*a-2*(x*u+y*v)*b)/det;
    double product = *signed_area * *signed_area / det;
    double hi = 0.5*(trace+sqrt(fmax(0,trace*trace-4*product))), lo = hi > 0 ? product/hi : 0;
    return lo > 0 && isfinite(hi) ? fmax(sqrt(hi),1/sqrt(lo)) : 1e300;
}

static double act_uv_energy(const AsmRun *run, const size_t *rels, size_t nr, size_t *pairs)
{
    double energy = 0; *pairs = 0;
    for (size_t j = 0; j < nr; j++) {
        const AsmRelation *r = &run->rels[rels[j]];
        for (int32_t k = 0; k < r->corr_count; k++) {
            const AsmCorr *c = &run->corr[(size_t)r->corr_first+k];
            if (c->valid != 3 || c->run < 0) continue;
            double a[2], b[2], ta[2], tb[2];
            act_point(&run->charts[r->a], c->va, NULL, a); act_point(&run->charts[r->b], c->vb, NULL, b);
            AsmContinuity_target(run, r, c, 0, tb); AsmContinuity_target(run, r, c, 1, ta);
            for (int q = 0; q < 2; q++) energy += 0.5*((a[q]-ta[q])*(a[q]-ta[q])+(b[q]-tb[q])*(b[q]-tb[q]));
            (*pairs)++;
        }
    }
    return energy;
}

static size_t act_uv_failed_sources(const AsmRun *run, const size_t *rels, size_t nr)
{
    size_t failed = 0;
    for (size_t j = 0; j < nr; j++) {
        double rms; size_t support;
        failed += !AsmContinuity_measure(run,&run->rels[rels[j]],&rms,&support);
    }
    return failed;
}

static void act_uv_log_conflicts(const AsmRun *run, const AsmChart *guard, size_t n,
                                  const AsmContinuityAuditCache *cache, int direction, double alpha)
{
    for (size_t j = 0; j < cache->pairs_count; j++) {
        const AsmContinuityUvPair *p = &cache->pairs[j];
        const AsmChart *a = NULL, *b = NULL;
        for (size_t i = 0; i < n; i++) {
            if (guard[i].id == p->chart_a) a = &guard[i];
            if (guard[i].id == p->chart_b) b = &guard[i];
        }
        if (!a || !b || p->vertex_a < 0 || p->vertex_b < 0 ||
            (size_t)p->vertex_a >= a->nv || (size_t)p->vertex_b >= b->nv) continue;
        double u[2], v[2], tu[2], tv[2];
        act_point(a,p->vertex_a,NULL,u); act_point(b,p->vertex_b,NULL,v);
        act_point(&run->charts[p->chart_a],p->vertex_a,NULL,tu);
        act_point(&run->charts[p->chart_b],p->vertex_b,NULL,tv);
        const float *x = a->xyz+3*(size_t)p->vertex_a, *y = b->xyz+3*(size_t)p->vertex_b;
        double d3 = sqrt((x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2]));
        double dl = run->chart_layer_d ? run->chart_layer_d[p->chart_a] : 0;
        if (dl <= 0) dl = run->layer_d_global;
        fprintf(stderr,"  [continuity UV conflict] direction %d alpha %.6g: %d:%d / %d:%d, UV %.9g -> %.9g, physical %.9g, tolerance %.9g; source-closing trial rejected\n",
                direction,alpha,p->chart_a,p->vertex_a,p->chart_b,p->vertex_b,
                hypot(v[0]-u[0],v[1]-u[1]),hypot(tv[0]-tu[0],tv[1]-tu[1]),d3,fmax(.5*dl,7.0));
    }
}

static void act_uv_stats(AsmChart *c, double band)
{
    size_t stressed = 0;
    c->area_uv = 0; c->sigma_lo = 1e300; c->sigma_hi = 0; c->n_flipped = 0;
    for (size_t j = 0; j < c->nf; j++) {
        const int32_t *f = &c->faces[3*j]; double lo, hi; int8_t sign;
        AsmFlatten_face_singular(c->xyz,f,1,c->uv,&lo,&hi,&sign);
        if (lo < c->sigma_lo) c->sigma_lo = lo;
        if (hi > c->sigma_hi) c->sigma_hi = hi;
        c->n_flipped += sign < 0;
        stressed += lo < 1-band || hi > 1+band || sign <= 0;
        double a = c->uv[(size_t)f[1]*2]-c->uv[(size_t)f[0]*2], b = c->uv[(size_t)f[1]*2+1]-c->uv[(size_t)f[0]*2+1];
        double x = c->uv[(size_t)f[2]*2]-c->uv[(size_t)f[0]*2], y = c->uv[(size_t)f[2]*2+1]-c->uv[(size_t)f[0]*2+1];
        c->area_uv += 0.5*fabs(a*y-b*x);
    }
    c->stress_frac = c->nf ? (double)stressed/c->nf : 0;
    if (stressed) c->flags |= ASM_CHART_SUSPECT;
}

/* Project one vertex onto its original correction disk. A large request
 * elsewhere in the component must not scale every feasible seam correction.
 * Round inward if storing the projected point in float crosses the bound. */
static void act_uv_bound(float point[2], const float origin[2])
{
    double du = (double)point[0]-origin[0], dv = (double)point[1]-origin[1];
    double length = hypot(du,dv);
    if (!(length > 8.0) || !isfinite(length)) return;
    point[0] = (float)(origin[0]+du*(8.0/length));
    point[1] = (float)(origin[1]+dv*(8.0/length));
    while (hypot((double)point[0]-origin[0],(double)point[1]-origin[1]) > 8.0) {
        point[0] = nextafterf(point[0],origin[0]);
        point[1] = nextafterf(point[1],origin[1]);
    }
}

/* Form the stored float trial in the chart frame. Both the ordinary line
 * search and the constrained direction use this exact rounding and bound. */
static double act_uv_trial_chart(const AsmChart *c, const AsmChart *origin, const int32_t *map,
                                 const ActVec *proposal, double alpha, float *uv)
{
    memcpy(uv,c->uv,2*c->nv*sizeof *uv);
    double ct = cos(c->pose_theta), sn = sin(c->pose_theta), largest = 0;
    for (size_t v = 0; v < c->nv; v++) if (map[v] >= 0) {
        ActVec d = proposal[map[v]];
        double du = ct*d.x+sn*d.y, dv = -sn*d.x+ct*d.y;
        if (c->flags & ASM_CHART_MIRROR) du = -du;
        uv[v*2] += (float)(alpha*du); uv[v*2+1] += (float)(alpha*dv);
        act_uv_bound(uv+v*2,(origin ? origin->uv : c->uv)+v*2);
        if (!isfinite(uv[v*2]) || !isfinite(uv[v*2+1]) ||
            (origin && hypot((double)uv[v*2]-origin->uv[v*2],(double)uv[v*2+1]-origin->uv[v*2+1]) > 8.0)) return INFINITY;
        largest = fmax(largest,hypot((double)uv[v*2]-c->uv[v*2],(double)uv[v*2+1]-c->uv[v*2+1]));
    }
    return largest;
}

typedef struct ActChartMetric { AsmMetricFace *ref; AsmMetricStats before; } ActChartMetric;

/* Reweight the existing displacement energy when a seam correction spends an
 * original chart's distortion budget. Translation remains free of strain;
 * increasing a small chart's stiffness asks its neighbours to share the
 * correction. This computes another direction in the same native solver.
 * The caller still tests stored float UVs, every retained source seam,
 * original-chart metric, movement bounds and complete collision context. */
static int act_uv_metric_direction(Arena_T arena, const AsmChart *saved, const AsmChart *origin,
                                   size_t n, const size_t *off, const int32_t *map,
                                   const ActChartMetric *metric, const ActSpring *spring,
                                   const int32_t *spring_chart, size_t ns,
                                   const ActMatch *match, size_t nm, const ActVec *initial,
                                   const ActVec *offset, ActVec *out, size_t active,
                                   size_t *passes, double *maximum_stiffness)
{
    Arena_Mark mark = Arena_save(arena);
    ActSpring *weighted = ARENA_ALLOC(arena,ns*sizeof *weighted);
    double *scale = ARENA_ALLOC(arena,n*sizeof *scale);
    float *trial = ARENA_ALLOC(arena,off[n]*2*sizeof *trial);
    ActVec *bounded = ARENA_ALLOC(arena,active*sizeof *bounded);
    UnionFind components = UF_new(arena,(int32_t)active);
    uint8_t *changed = ARENA_CALLOC(arena,active,1);
    for (size_t j = 0; j < ns; j++) if (spring[j].a >= 0 && spring[j].b >= 0)
        uf_union(&components,spring[j].a,spring[j].b);
    for (size_t j = 0; j < nm; j++) {
        int32_t first = -1;
        for (int k = 0; k < match[j].count; k++) if (match[j].id[k] >= 0 && match[j].coefficient[k] != 0) {
            if (first < 0) first = match[j].id[k]; else uf_union(&components,first,match[j].id[k]);
        }
    }
    for (size_t i = 0; i < n; i++) scale[i] = 1;
    memcpy(out,initial,active*sizeof *out); *passes = 0; *maximum_stiffness = 1;
    int ok = 1, solved = 0;
    for (size_t pass = 0; pass < 8; pass++) {
        size_t damaged = 0;
        for (size_t i = 0; i < n; i++) {
            AsmMetricStats measured;
            float *uv = trial+2*off[i];
            if (!isfinite(act_uv_trial_chart(&saved[i],origin ? &origin[i] : NULL,map+off[i],out,1,uv))) { ok = 0; break; }
            AsmMetric_measure(&saved[i],metric[i].ref,uv,&measured);
            if (!AsmMetric_preserved(&metric[i].before,&measured)) {
                scale[i] *= 4; *maximum_stiffness = fmax(*maximum_stiffness,scale[i]); damaged++;
            }
        }
        if (!ok || !damaged) break;
        for (size_t j = 0; j < ns; j++) {
            weighted[j] = spring[j];
            double weight = (double)spring[j].weight*scale[spring_chart[j]];
            if (!isfinite(weight) || weight > FLT_MAX) { ok = 0; break; }
            weighted[j].weight = (float)weight;
            if (scale[spring_chart[j]] != 1) {
                if (spring[j].a >= 0) changed[uf_find(&components,spring[j].a)] = 1;
                if (spring[j].b >= 0) changed[uf_find(&components,spring[j].b)] = 1;
            }
        }
        if (!ok) break;
        Arena_Mark solve_mark = Arena_save(arena);
        int iterations = 0;
        memset(out,0,active*sizeof *out);
        ok = act_uv_solve(arena,weighted,ns,match,nm,out,active,&iterations);
        Arena_restore(arena,solve_mark);
        if (!ok) break;
        int needs_bound = 0;
        for (size_t j = 0; j < active; j++)
            needs_bound |= hypot(out[j].x+offset[j].x,out[j].y+offset[j].y) > 8;
        if (needs_bound) {
            solve_mark = Arena_save(arena);
            ok = act_uv_bounded_direction(arena,weighted,ns,match,nm,out,offset,bounded,active,&iterations);
            if (ok) memcpy(out,bounded,active*sizeof *out);
            Arena_restore(arena,solve_mark);
            if (!ok) break;
        }
        /* The quadratic and movement disks factor across disconnected
         * components. Preserve the already solved independent components
         * exactly; a stiffer distant chart cannot change their gauge through
         * a different global iterative stopping norm. */
        for (size_t j = 0; j < active; j++) if (!changed[uf_find(&components,(int32_t)j)]) out[j] = initial[j];
        *passes = pass+1; solved = 1;
    }
    Arena_restore(arena,mark); return ok && solved;
}

/* A strained face or new overlap can reject every global step while other
 * boundaries still have room to close. Hold the affected vertices at the checkpoint and
 * solve the same mesh/match quadratic for the remaining vertices. These
 * are Dirichlet constraints, not a discontinuous post-solve clamp. Each
 * solve must converge; the caller still applies all original audits and
 * retains its independently saved incumbent if this direction loses. */
/* Both constrained directions solve the same component quadratic. Only
 * their held-vertex set can differ. Retain one converged solution across
 * those directions; the cache belongs to this component's immutable base
 * springs/matches and is allocated outside either helper's scratch mark. */
typedef struct ActUvDirectionCache {
    uint8_t *fixed;
    ActVec *solution;
    size_t hits;
    int valid;
} ActUvDirectionCache;

static int act_uv_locked_direction(Arena_T arena, const AsmChart *saved, const AsmChart *origin,
                                   size_t n, const size_t *off, const int32_t *map,
                                   const ActSpring *spring, size_t ns, const ActMatch *match, size_t nm,
                                   const ActVec *linear, ActVec *out, size_t active, size_t *locked,
                                   AsmRun *run, AsmContinuityUvAudit audit, const AsmContinuityAuditCache *cache,
                                   int hold_chart, ActUvDirectionCache *solve_cache, const ActChartMetric *metric)
{
    Arena_Mark mark = Arena_save(arena);
    uint8_t *fixed = ARENA_CALLOC(arena,active,1);
    float *trial = ARENA_ALLOC(arena,off[n]*2*sizeof *trial);
    ActSpring *springs = ARENA_ALLOC(arena,ns*sizeof *springs);
    ActMatch *matches = ARENA_ALLOC(arena,nm*sizeof *matches);
    /* A first audit may allocate its baseline inside this helper's mark.
     * Keep a local cache so no such pointer escapes the scratch restore. */
    AsmContinuityAuditCache feedback = {0};
    if (cache) feedback = *cache;
    feedback.scratch = arena;
    float **previous_uv = NULL;
    if (run && audit) {
        feedback.blocked_charts = ARENA_CALLOC(arena,run->n_charts ? run->n_charts : 1,1);
        feedback.blocked_charts_count = run->n_charts;
        previous_uv = ARENA_ALLOC(arena,n*sizeof *previous_uv);
        for (size_t i = 0; i < n; i++) previous_uv[i] = run->charts[saved[i].id].uv;
    }
    memcpy(out,linear,active*sizeof *out); *locked = 0;
    int ok = 1;
    for (int pass = 0; pass < 4; pass++) {
        size_t added = 0; int shape_ok = 1;
        for (size_t i = 0; i < n && ok; i++) {
            const AsmChart *c = &saved[i]; float *uv = trial+2*off[i]; int chart_ok = 1;
            if (!isfinite(act_uv_trial_chart(c,origin ? &origin[i] : NULL,map+off[i],out,1.0,uv))) { ok = 0; break; }
            for (size_t f = 0; f < c->nf; f++) {
                const int32_t *face = c->faces+3*f;
                if (map[off[i]+face[0]] < 0 && map[off[i]+face[1]] < 0 && map[off[i]+face[2]] < 0) continue;
                double old_area, new_area;
                double old_q = act_uv_quality(c,face,c->uv,&old_area), new_q = act_uv_quality(c,face,uv,&new_area);
                if (isfinite(new_area) && (fabs(old_area) <= 1e-10 || old_area*new_area > 0) &&
                    new_q <= fmax(2.0,old_q)*(1+1e-6) &&
                    (!metric || AsmMetric_local_face_preserved(c,f,metric[i].ref+f,uv))) continue;
                shape_ok = chart_ok = 0;
                for (int k = 0; k < 3; k++) {
                    int32_t q = map[off[i]+face[k]];
                    if (q >= 0 && !fixed[q]) { fixed[q] = 1; added++; }
                }
            }
            if (metric) {
                AsmMetricStats measured;
                AsmMetric_measure(c,metric[i].ref,uv,&measured);
                if (!AsmMetric_preserved(&metric[i].before,&measured)) {
                    shape_ok = chart_ok = 0;
                    int band10 = measured.within10 < fmin(.95*measured.area,metric[i].before.within10);
                    int band25 = measured.within25 < fmin(.99*measured.area,metric[i].before.within25);
                    for (size_t f = 0; f < c->nf; f++) {
                        AsmChart triangle = *c; AsmMetricStats was, now;
                        triangle.faces = c->faces+3*f; triangle.nf = 1;
                        AsmMetric_measure(&triangle,metric[i].ref+f,c->uv,&was);
                        AsmMetric_measure(&triangle,metric[i].ref+f,uv,&now);
                        if (!(now.invalid || (band10 && now.within10 < was.within10) ||
                                              (band25 && now.within25 < was.within25))) continue;
                        /* Preserve the newly damaged original face, then
                         * resolve the same source system elsewhere. This is
                         * a held-variable constraint, never a UV clamp or
                         * removal of a source observation. */
                        for (size_t k = 0; k < 3; k++) {
                            int32_t q = map[off[i]+triangle.faces[k]];
                            if (q >= 0 && !fixed[q]) { fixed[q] = 1; added++; }
                        }
                    }
                }
            }
            if (!chart_ok && hold_chart) for (size_t v = 0; v < c->nv; v++) {
                int32_t q = map[off[i]+v];
                if (q >= 0 && !fixed[q]) { fixed[q] = 1; added++; }
            }
        }
        if (ok && shape_ok && run && audit) {
            for (size_t i = 0; i < n; i++) run->charts[saved[i].id].uv = trial+2*off[i];
            int geometry_ok = audit(run,saved,n,&feedback);
            for (size_t i = 0; i < n; i++) run->charts[saved[i].id].uv = previous_uv[i];
            if (!geometry_ok) {
                /* Holding both sides of a reported collision prevents the
                 * re-solve from moving its other side into the same conflict. */
                for (size_t i = 0; i < n; i++) if (feedback.blocked_charts[saved[i].id]) {
                    for (size_t v = 0; v < saved[i].nv; v++) {
                        int32_t q = map[off[i]+v];
                        if (q >= 0 && !fixed[q]) { fixed[q] = 1; added++; }
                    }
                }
                if (!added) ok = 0; /* no usable feedback: retain the incumbent */
            }
        } else if (ok && !shape_ok && !added) ok = 0;
        if (!ok || !added) break;
        *locked += added;
        for (size_t j = 0; j < ns; j++) {
            springs[j] = spring[j];
            if (springs[j].a >= 0 && fixed[springs[j].a]) springs[j].a = -1;
            if (springs[j].b >= 0 && fixed[springs[j].b]) springs[j].b = -1;
        }
        for (size_t j = 0; j < nm; j++) {
            matches[j] = match[j];
            for (int k = 0; k < matches[j].count; k++)
                if (matches[j].id[k] >= 0 && fixed[matches[j].id[k]]) matches[j].id[k] = -1;
        }
        if (solve_cache && solve_cache->valid && !memcmp(solve_cache->fixed,fixed,active)) {
            memcpy(out,solve_cache->solution,active*sizeof *out);
            solve_cache->hits++;
        } else {
            memset(out,0,active*sizeof *out);
            Arena_Mark sm = Arena_save(arena); int iterations = 0;
            ok = act_uv_solve(arena,springs,ns,matches,nm,out,active,&iterations);
            Arena_restore(arena,sm);
            if (ok && solve_cache) {
                memcpy(solve_cache->fixed,fixed,active);
                memcpy(solve_cache->solution,out,active*sizeof *out);
                solve_cache->valid = 1;
            }
        }
        if (!ok) break;
    }
    Arena_restore(arena,mark); return ok && *locked > 0;
}

static int act_uv_component_guarded(AsmRun *run, Arena_T arena, const int32_t *ids, size_t n,
                            const size_t *rels, size_t nr, AsmContinuityUvAudit audit, double band, FILE *fp,
                            const AsmChart *origin, const AsmChart *guard, size_t n_guard, int first_step)
{
    if (!n || !nr) return 0;
    double started = ves_clock_sec();
    AsmChart *saved = ARENA_ALLOC(arena, n*sizeof *saved);
    size_t *off = ARENA_ALLOC(arena, (n+1)*sizeof *off); off[0] = 0;
    int32_t *index = ARENA_ALLOC(arena, run->n_charts*sizeof *index);
    for (size_t i = 0; i < run->n_charts; i++) index[i] = -1;
    size_t nf = 0, nm = 0;
    for (size_t i = 0; i < n; i++) {
        saved[i] = run->charts[ids[i]]; index[ids[i]] = (int32_t)i;
        if (saved[i].nv > INT32_MAX-off[i] || saved[i].nf > INT32_MAX/3-nf) return 0;
        off[i+1] = off[i]+saved[i].nv; nf += saved[i].nf;
    }
    size_t nv = off[n];
    if (!nv || !nf) return 0;
    for (size_t i = 0; i < nr; i++) {
        int32_t count = run->rels[rels[i]].corr_count;
        if (count < 0 || (size_t)count > (INT32_MAX-nm)/2) return 0;
        nm += 2*(size_t)count;
    }
    if (!nm || nm > INT32_MAX) return 0;
    /* Scratch is released after each component. Do not allocate an unbounded
     * corpus-wide mesh system or silently accept a budget-limited solve. */
    if ((double)nv*280+(double)nf*(120+sizeof(AsmMetricFace)+3*sizeof(int32_t))+(double)nm*(sizeof(ActMatch)+6*sizeof(ActIcEntry)) > 2147483648.0) {
        fprintf(stderr, "  [continuity UV] component %d deferred: scratch bound\n", saved[0].component); return 0;
    }
    ActChartMetric *metric = ARENA_ALLOC(arena,n*sizeof *metric);
    for (size_t i = 0; i < n; i++) {
        metric[i].ref = ARENA_ALLOC(arena,saved[i].nf*sizeof *metric[i].ref);
        if (AsmMetric_prepare(&saved[i],metric[i].ref)) return 0;
        AsmMetric_measure(&saved[i],metric[i].ref,saved[i].uv,&metric[i].before);
        if (metric[i].before.invalid) return 0;
    }
    ActMatch *match = ARENA_ALLOC(arena, nm*sizeof *match); nm = 0;
    for (size_t j = 0; j < nr; j++) {
        const AsmRelation *r = &run->rels[rels[j]];
        for (int32_t k = 0; k < r->corr_count; k++) {
            const AsmCorr *c = &run->corr[(size_t)r->corr_first+k];
            if (c->valid != 3 || c->run < 0) continue;
            if ((index[r->a] < 0 && index[r->b] < 0) || c->va < 0 || c->vb < 0 ||
                (size_t)c->va >= run->charts[r->a].nv || (size_t)c->vb >= run->charts[r->b].nv) return 0;
            double a[2], b[2], ta[2], tb[2];
            act_point(&run->charts[r->a], c->va, NULL, a); act_point(&run->charts[r->b], c->vb, NULL, b);
            if (!AsmContinuity_target(run,r,c,0,tb) || !AsmContinuity_target(run,r,c,1,ta)) return 0;
            /* Each physical witness contributes two half-weight rows. Their
             * affine endpoint coefficients differ as the charts deform, so
             * averaging their targets into one vertex pair is incorrect. */
            for (int direction = 0; direction < 2; direction++) {
                ActMatch m = {0}; m.weight = .5;
                m.target = direction ? (ActVec){ta[0]-a[0],ta[1]-a[1]} : (ActVec){b[0]-tb[0],b[1]-tb[1]};
                for (int side = 0; side < 2; side++) {
                    int32_t chart = side ? r->b : r->a, vertex = side ? c->vb : c->va;
                    int32_t vertices[3] = {vertex,0,0}; double weights[3] = {1,0,0}; int count = 1;
                    if (side == direction) count = act_trim_weights(&run->charts[chart],vertex,
                        side ? &c->trim_b : &c->trim_a,side ? c->gap_b : c->gap_a,side ? -1 : 1,vertices,weights);
                    if (!count) return 0;
                    if (index[chart] >= 0) for (int k = 0; k < count; k++)
                        act_match_add(&m,(int32_t)(off[index[chart]]+vertices[k]),(side ? -1 : 1)*weights[k]);
                }
                match[nm++] = m;
            }
        }
    }
    if (!nm) return 0;
    ActSpring *spring = ARENA_ALLOC(arena, nf*3*sizeof *spring); size_t ns = 0;
    int32_t *spring_chart = ARENA_ALLOC(arena,nf*3*sizeof *spring_chart);
    size_t *adj_off = ARENA_CALLOC(arena, nv+1, sizeof *adj_off);
    for (size_t i = 0; i < n; i++) {
        const AsmChart *c = &saved[i];
        for (size_t j = 0; j < c->nf; j++) {
            const int32_t *f = &c->faces[j*3];
            for (int k = 0; k < 3; k++) if (f[k] < 0 || (size_t)f[k] >= c->nv) return 0;
            for (int k = 0; k < 3; k++) {
                int32_t a = f[k], b = f[(k+1)%3], v = f[(k+2)%3];
                double dot = 0, aa = 0, bb = 0;
                for (int q = 0; q < 3; q++) {
                    double x = c->xyz[(size_t)a*3+q]-c->xyz[(size_t)v*3+q];
                    double y = c->xyz[(size_t)b*3+q]-c->xyz[(size_t)v*3+q];
                    dot += x*y; aa += x*x; bb += y*y;
                }
                double cross = sqrt(fmax(0,aa*bb-dot*dot));
                double weight = cross > 1e-10 ? 0.5*fmax(0,dot)/cross : 0;
                double length = act_dist(c, a, b);
                if (!isfinite(weight) || !isfinite(length) || weight > 1e10) return 0;
                a += (int32_t)off[i]; b += (int32_t)off[i];
                spring_chart[ns] = (int32_t)i;
                spring[ns++] = (ActSpring){a,b,(float)weight,(float)length}; adj_off[a+1]++; adj_off[b+1]++;
            }
        }
    }
    for (size_t i = 0; i < nv; i++) adj_off[i+1] += adj_off[i];
    int32_t *adj = ARENA_ALLOC(arena, 2*ns*sizeof *adj);
    size_t *fill = ARENA_CALLOC(arena, nv, sizeof *fill);
    for (size_t i = 0; i < ns; i++) {
        adj[adj_off[spring[i].a]+fill[spring[i].a]++] = (int32_t)i;
        adj[adj_off[spring[i].b]+fill[spring[i].b]++] = (int32_t)i;
    }
    double *dist = ARENA_ALLOC(arena, nv*sizeof *dist);
    int32_t *map = ARENA_ALLOC(arena, nv*sizeof *map), *heap = ARENA_ALLOC(arena, nv*sizeof *heap);
    for (size_t i = 0; i < nv; i++) { dist[i] = 1e300; map[i] = -1; }
    size_t nh = 0;
    for (size_t i = 0; i < nm; i++) for (int k = 0; k < match[i].count; k++) {
        int32_t v = match[i].id[k];
        if (match[i].coefficient[k] == 0) continue;
        if (v < 0 || map[v] >= 0) continue;
        dist[v] = 0; heap[nh] = v; map[v] = (int32_t)nh; nh++;
    }
    while (nh) {
        int32_t v = act_heap_pop(heap, map, dist, &nh);
        for (size_t j = adj_off[v]; j < adj_off[v+1]; j++) {
            const ActSpring *s = &spring[adj[j]]; int32_t w = s->a == v ? s->b : s->a;
            double d = dist[v]+s->length;
            if (map[w] == -2 || d > 24.0 || d >= dist[w]) continue;
            dist[w] = d;
            if (map[w] < 0) { map[w] = (int32_t)nh; heap[nh++] = w; }
            act_heap_up(heap, map, dist, (size_t)map[w]);
        }
    }
    size_t active = 0;
    for (size_t i = 0; i < nv; i++) map[i] = dist[i] <= 24.0 ? (int32_t)active++ : -1;
    size_t keep = 0;
    for (size_t i = 0; i < ns; i++) {
        ActSpring s = spring[i]; s.a = map[s.a]; s.b = map[s.b];
        if ((s.a >= 0 || s.b >= 0) && s.weight > 0) { spring_chart[keep] = spring_chart[i]; spring[keep++] = s; }
    }
    ns = keep;
    for (size_t i = 0; i < nm; i++) {
        for (int k = 0; k < match[i].count; k++) if (match[i].id[k] >= 0)
            match[i].id[k] = map[match[i].id[k]];
    }
    ActVec *delta = ARENA_CALLOC(arena, active, sizeof *delta); int iterations = 0;
    double solve_started = ves_clock_sec(), prepare_sec = solve_started-started;
    Arena_Mark solve_mark = Arena_save(arena);
    int converged = act_uv_solve(arena, spring, ns, match, nm, delta, active, &iterations);
    Arena_restore(arena,solve_mark);
    if (!converged) {
        fprintf(stderr, "  [continuity UV] component %d deferred: solve did not converge (%d iterations; prepare %.3f s, solve %.3f s)\n",
                saved[0].component,iterations,prepare_sec,ves_clock_sec()-solve_started); return 0;
    }
    double max_delta = 0;
    for (size_t i = 0; i < active; i++) {
        double d = hypot(delta[i].x,delta[i].y); if (!isfinite(d)) return 0;
        if (d > max_delta) max_delta = d;
    }
    if (max_delta < 1e-4) return 0;
    ActVec *offset = ARENA_CALLOC(arena,active,sizeof *offset), *bounded = NULL;
    int needs_bound = 0, bound_iterations = 0;
    for (size_t i = 0; i < n; i++) {
        const AsmChart *c = &saved[i]; double ct = cos(c->pose_theta), sn = sin(c->pose_theta);
        for (size_t v = 0; v < c->nv; v++) if (map[off[i]+v] >= 0) {
            int32_t q = map[off[i]+v];
            if (origin) {
                double u = (double)c->uv[2*v]-origin[i].uv[2*v], w = (double)c->uv[2*v+1]-origin[i].uv[2*v+1];
                if (c->flags & ASM_CHART_MIRROR) u = -u;
                offset[q] = (ActVec){ct*u-sn*w,sn*u+ct*w};
            }
            needs_bound |= hypot(delta[q].x+offset[q].x,delta[q].y+offset[q].y) > 8.0;
        }
    }
    if (needs_bound) {
        bounded = ARENA_ALLOC(arena,active*sizeof *bounded);
        double bounded_start = ves_clock_sec();
        int converged = act_uv_bounded_direction(arena,spring,ns,match,nm,delta,offset,bounded,active,&bound_iterations);
        fprintf(stderr,"  [continuity UV] bounded direction: component %d, %zu vertices, %d iterations, %.3f s; %s\n",
                saved[0].component,active,bound_iterations,ves_clock_sec()-bounded_start,converged ? "converged" : "deferred, keeping incumbent direction");
        if (!converged) bounded = NULL;
    }
    double solved = ves_clock_sec(), solve_sec = solved-solve_started, audit_sec = 0;
    /* Protect every previously passing incident SOURCE seam, including
     * constraints not used in the local fit. */
    uint8_t *passed = ARENA_CALLOC(arena, run->n_rels ? run->n_rels : 1, 1);
    for (size_t j = 0; j < run->n_rels; j++) {
        const AsmRelation *r = &run->rels[j];
        if (index[r->a] < 0 && index[r->b] < 0) continue;
        if (!AsmChart_in_layout(&run->charts[r->a]) || !AsmChart_in_layout(&run->charts[r->b])) continue;
        double rms; size_t support; passed[j] = (uint8_t)AsmContinuity_measure(run, r, &rms, &support);
    }
    float *trial = ARENA_ALLOC(arena, nv*2*sizeof *trial);
    size_t pairs; double before = act_uv_energy(run, rels, nr, &pairs), after = before;
    double accepted_delta = 0, best_after = before; int accepted = 0, tries = 0, kept_tries = 0;
    float *best_trial = ARENA_ALLOC(arena,nv*2*sizeof *best_trial);
    ActVec *local = NULL, *metric_delta = NULL;
    ActUvDirectionCache solve_cache = {0};
    int shape_limited = 0, geometry_limited = 0, metric_limited = 0;
    AsmContinuityAuditCache audit_cache = {0}; audit_cache.scratch = arena;
    audit_cache.guard_charts = guard; audit_cache.guard_count = n_guard;
    size_t failed_before = guard ? act_uv_failed_sources(run,rels,nr) : 0;
    AsmContinuityUvPair feedback_pairs[8];
    if (failed_before) { audit_cache.pairs = feedback_pairs; audit_cache.pairs_capacity = 8; }
    for (int direction = 0; direction < 5; direction++) {
    /* A diagnostic may stop after a fully checked improving direction.
     * Later directions seek a better objective, not feasibility. */
    if (first_step && accepted) break;
    if (direction == 1 && !bounded) continue;
    if (direction == 2) {
        if (!metric_limited || (double)nv*400+(double)nf*(232+3*sizeof(int32_t))+
            (double)nm*2*(sizeof(ActMatch)+6*sizeof(ActIcEntry)) > 2147483648.0) continue;
        metric_delta = ARENA_ALLOC(arena,active*sizeof *metric_delta);
        size_t metric_passes; double stiffness, metric_started = ves_clock_sec();
        int metric_ok = act_uv_metric_direction(arena,saved,origin,n,off,map,metric,spring,spring_chart,
            ns,match,nm,bounded ? bounded : delta,offset,metric_delta,active,&metric_passes,&stiffness);
        fprintf(stderr,"  [continuity UV] metric direction: component %d, %zu solves, maximum stiffness %.6g, %.3f s; %s\n",
            saved[0].component,metric_passes,stiffness,ves_clock_sec()-metric_started,
            metric_ok ? "solved; original chart limits remain required" : "deferred, keeping incumbent direction");
        if (!metric_ok) continue;
    }
    if (direction >= 3) {
        if (!(shape_limited || geometry_limited) || (double)nv*432+(double)nf*(200+sizeof(AsmMetricFace)+3*sizeof(int32_t))+
            (double)nm*2*(sizeof(ActMatch)+6*sizeof(ActIcEntry)) > 2147483648.0) continue;
        if (!local) {
            local = ARENA_ALLOC(arena,active*sizeof *local);
            /* The extra cache is optional under the existing scratch cap.
             * A large component can still use the original uncached solve. */
            if ((double)nv*432+(double)active*(sizeof(ActVec)+1)+(double)nf*(200+sizeof(AsmMetricFace)+3*sizeof(int32_t))+
                (double)nm*2*(sizeof(ActMatch)+6*sizeof(ActIcEntry)) <= 2147483648.0) {
                solve_cache.fixed = ARENA_ALLOC(arena,active);
                solve_cache.solution = ARENA_ALLOC(arena,active*sizeof *solve_cache.solution);
            }
        }
        size_t locked = 0;
        double local_started = ves_clock_sec();
        int local_ok = act_uv_locked_direction(arena,saved,origin,n,off,map,spring,ns,match,nm,
                                               bounded ? bounded : delta,local,active,&locked,run,audit,&audit_cache,direction == 4,
                                               solve_cache.fixed ? &solve_cache : NULL,metric);
        fprintf(stderr,"  [continuity UV] %s direction: component %d, %zu vertices held, %.3f s; %s\n",
                direction == 4 ? "chart-constrained" : "constrained",saved[0].component,locked,ves_clock_sec()-local_started,
                local_ok ? "solved" : "deferred, keeping incumbent direction");
        if (!local_ok) continue;
    }
    size_t rejected[4] = {0};
    size_t rejected_metric = 0;
    int32_t bad_chart = -1; size_t bad_face = 0;
    double bad_old_quality = 0, bad_new_quality = 0, direction_before = best_after;
    int reported_conflicts = 0;
    double alpha = 1.0;
    const ActVec *proposal = direction >= 3 ? local : direction == 2 ? metric_delta : direction == 1 ? bounded : delta;
    tries = 0;
    for (; tries < 8; tries++, alpha *= 0.5) {
        int ok = 1; double attempt_delta = 0;
        for (size_t i = 0; i < n; i++) {
            const AsmChart *c = &saved[i]; float *uv = trial+2*off[i];
            double step = act_uv_trial_chart(c,origin ? &origin[i] : NULL,map+off[i],proposal,alpha,uv);
            if (!isfinite(step)) ok = 0;
            attempt_delta = fmax(attempt_delta,step);
            for (size_t f = 0; f < c->nf && ok; f++) {
                const int32_t *face = &c->faces[f*3];
                if (map[off[i]+face[0]] < 0 && map[off[i]+face[1]] < 0 && map[off[i]+face[2]] < 0) continue;
                double old_area, new_area;
                double old_q = act_uv_quality(c,face,c->uv,&old_area), new_q = act_uv_quality(c,face,uv,&new_area);
                if (!isfinite(new_area) || (fabs(old_area) > 1e-10 && old_area*new_area <= 0) ||
                    new_q > fmax(2.0,old_q)*(1+1e-6) ||
                    !AsmMetric_local_face_preserved(c,f,metric[i].ref+f,uv)) {
                    ok = 0;
                    if (bad_chart < 0) { bad_chart = c->id; bad_face = f; bad_old_quality = old_q; bad_new_quality = new_q; }
                }
            }
            if (ok) {
                AsmMetricStats measured;
                AsmMetric_measure(c,metric[i].ref,uv,&measured);
                if (!AsmMetric_preserved(&metric[i].before,&measured)) {
                    ok = 0; metric_limited = 1; rejected_metric++;
                }
            }
            run->charts[ids[i]].uv = uv;
        }
        if (!ok) rejected[0]++;
        if (ok) for (size_t j = 0; j < run->n_rels; j++) if (passed[j]) {
            double rms; size_t support;
            if (!AsmContinuity_measure(run, &run->rels[j], &rms, &support)) { ok = 0; rejected[1]++; break; }
        }
        if (ok) {
            after = act_uv_energy(run,rels,nr,&pairs);
            ok = isfinite(after) && after < best_after-fmax(1e-8,before*1e-6);
            if (!ok) rejected[2]++;
        }
        if (ok && audit) {
            double audit_started = ves_clock_sec();
            ok = audit(run,saved,n,&audit_cache); audit_sec += ves_clock_sec()-audit_started;
            if (!ok) rejected[3]++;
            if (!ok && failed_before && !reported_conflicts && act_uv_failed_sources(run,rels,nr) < failed_before) {
                act_uv_log_conflicts(run,guard,n_guard,&audit_cache,direction,alpha); reported_conflicts = 1;
            }
        }
        if (ok) {
            accepted = 1; accepted_delta = attempt_delta; best_after = after; kept_tries = tries;
            memcpy(best_trial,trial,nv*2*sizeof *trial);
            break;
        }
    }
    shape_limited |= rejected[0] > 0;
    geometry_limited |= rejected[3] > 0;
    if (bounded || !accepted || direction >= 2) fprintf(stderr,"  [continuity UV trials] component %d, %s direction: %s; rejected shape %zu source %zu energy %zu geometry %zu (chart metric %zu); first shape chart %d face %zu quality %.6g -> %.6g\n",
                                     saved[0].component,direction == 4 ? "chart-constrained" : direction == 3 ? "constrained" : direction == 2 ? "metric" : direction == 1 ? "bounded" : "clipped",best_after < direction_before ? "improved" : "unchanged",
                                     rejected[0],rejected[1],rejected[2],rejected[3],rejected_metric,bad_chart,bad_face,bad_old_quality,bad_new_quality);
    }
    if (solve_cache.hits) fprintf(stderr,"  [continuity UV] component %d: reused %zu identical constrained solves; all trial audits retained\n",
                                 saved[0].component,solve_cache.hits);
    after = best_after; if (accepted) tries = kept_tries;
    for (size_t i = 0; i < n; i++) {
        if (accepted) memcpy(saved[i].uv,best_trial+2*off[i],saved[i].nv*2*sizeof(float));
        run->charts[ids[i]] = saved[i];
        if (accepted) act_uv_stats(&run->charts[ids[i]],band);
    }
    double check_sec = ves_clock_sec()-solved-audit_sec;
    if (fp) fprintf(fp,"%d,%zu,%zu,%zu,%d,%d,%.9g,%.9g,%.9g,%.6f,%.6f,%.6f,%.6f\n",saved[0].component,n,active,nm,iterations,accepted,
                    sqrt(before/nm),sqrt((accepted ? after : before)/nm),accepted_delta,
                    prepare_sec,solve_sec,check_sec,audit_sec);
    else fprintf(stderr,"  [continuity UV] %zu charts, %zu active vertices, %zu pairs: RMS %.3f -> %.3f, maximum correction %.3f vox, line search %d; %s, %d iterations, prepare %.3f s solve %.3f s checks %.3f s audit %.3f s\n",
                 n,active,nm,sqrt(before/nm),sqrt((accepted ? after : before)/nm),accepted_delta,tries,
                 accepted ? "accepted" : "reverted",iterations,prepare_sec,solve_sec,check_sec,audit_sec);
    return accepted;
}

static int act_uv_component(AsmRun *run, Arena_T arena, const int32_t *ids, size_t n,
                            const size_t *rels, size_t nr, AsmContinuityUvAudit audit, double band, FILE *fp,
                            const AsmChart *origin)
{
    return act_uv_component_guarded(run,arena,ids,n,rels,nr,audit,band,fp,origin,NULL,0,0);
}

static int act_relax_boundary_trial(AsmRun *run, const AsmChart *saved, size_t n,
                                       const uint8_t *selected, AsmContinuityUvAudit audit, double band, int first_step)
{
    if (!n) return 0;
    Arena_T arena = Arena_new();
    uint8_t *variable = ARENA_CALLOC(arena,run->n_charts,1);
    int32_t *ids = ARENA_ALLOC(arena,n*sizeof *ids);
    for (size_t i = 0; i < n; i++) { ids[i] = saved[i].id; variable[ids[i]] = 1; }
    size_t *rels = ARENA_ALLOC(arena,(run->n_rels ? run->n_rels : 1)*sizeof *rels);
    size_t nr = 0, fixed = 0; int feasible = 1;
    for (size_t i = 0; i < run->n_rels; i++) {
        const AsmRelation *r = &run->rels[i];
        if (!(r->continuity & ASM_CONT_SOURCE) || (r->flags & (ASM_REL_CONTACT | ASM_REL_DROPPED))) continue;
        if (!AsmRel_is_join(r) && !(r->continuity & ASM_CONT_CUT) && (!selected || !selected[i])) continue;
        int a = variable[r->a], b = variable[r->b];
        if (!a && !b) continue;
        const AsmChart *A = &run->charts[r->a], *B = &run->charts[r->b];
        if (!AsmChart_in_layout(A) || !AsmChart_in_layout(B)) continue;
        if ((!a && A->placement_state == ASM_PLACE_NONE) || (!b && B->placement_state == ASM_PLACE_NONE)) continue;
        double rms; size_t support; AsmContinuity_measure(run,r,&rms,&support);
        /* Local UV correction cannot rescue another wrap. Every required
         * seam must already be within twice its ordinary residual gate. */
        if (support < 6 || !isfinite(rms) || rms > 2*ASM_CONTINUITY_GATE(r->rms)) { feasible = 0; break; }
        rels[nr++] = i; fixed += a != b;
    }
    int accepted = feasible && fixed && nr &&
        act_uv_component_guarded(run,arena,ids,n,rels,nr,audit,band,NULL,saved,NULL,0,first_step);
    Arena_dispose(&arena); return accepted;
}

int AsmContinuity_relax_boundary_trial(AsmRun *run, const AsmChart *saved, size_t n,
                                       const uint8_t *selected, AsmContinuityUvAudit audit, double band)
{ return act_relax_boundary_trial(run,saved,n,selected,audit,band,0); }

int AsmContinuity_relax_boundary_step(AsmRun *run, const AsmChart *saved, size_t n,
                                      const uint8_t *selected, AsmContinuityUvAudit audit, double band)
{ return act_relax_boundary_trial(run,saved,n,selected,audit,band,1); }

static size_t act_write_witnesses(const AsmRun *run, FILE *fp)
{
    size_t written = 0;
    if (fp) fputs("relation,chart_a,chart_b,flags,source,va,vb,valid,run,gap_au,gap_av,gap_bu,gap_bv\n",fp);
    if (!run->rels || !run->corr) return 0;
    for (size_t j = 0; j < run->n_rels; j++) {
        const AsmRelation *r = &run->rels[j];
        if (r->corr_first < 0 || r->corr_count <= 0 || (size_t)r->corr_first > run->n_corr ||
            (size_t)r->corr_count > run->n_corr-(size_t)r->corr_first) continue;
        for (int32_t k = 0; k < r->corr_count; k++) {
            const AsmCorr *c = &run->corr[(size_t)r->corr_first+k];
            if (c->valid != 3) continue;
            if (fp) fprintf(fp,"%zu,%d,%d,%u,%u,%d,%d,%u,%d,%.9g,%.9g,%.9g,%.9g\n",
                            j,r->a,r->b,r->flags,(unsigned)((r->continuity & ASM_CONT_SOURCE) != 0),c->va,c->vb,(unsigned)c->valid,c->run,
                            c->gap_a[0],c->gap_a[1],c->gap_b[0],c->gap_b[1]);
            written++;
        }
    }
    return written;
}

size_t AsmContinuity_relax_uv(AsmRun *run, AsmContinuityUvAudit audit, double band, const char *diag_dir)
{
    /* Exact original constraints for independent reconstruction/debugging.
     * Registration does not mutate these trim frames or witness indices. */
    if (diag_dir) {
        char path[2048]; snprintf(path,sizeof path,"%s/stage5_source_witnesses.csv",diag_dir);
        FILE *witnesses = fopen(path,"wb");
        if (witnesses) {
            act_write_witnesses(run,witnesses);
            fclose(witnesses);
        }
    }
    Arena_T arena = Arena_new(); size_t nc = run->n_charts, ng = 0;
    for (size_t i = 0; i < nc; i++) if (AsmChart_in_layout(&run->charts[i]) && run->charts[i].component >= 0 &&
        (size_t)run->charts[i].component >= ng) ng = (size_t)run->charts[i].component+1;
    size_t *off = ARENA_CALLOC(arena,ng+1,sizeof *off), *fill = ARENA_CALLOC(arena,ng+1,sizeof *fill);
    for (size_t i = 0; i < nc; i++) if (AsmChart_in_layout(&run->charts[i]) && run->charts[i].component >= 0) off[run->charts[i].component+1]++;
    for (size_t i = 0; i < ng; i++) off[i+1] += off[i];
    int32_t *ids = ARENA_ALLOC(arena,(nc ? nc : 1)*sizeof *ids);
    for (size_t i = 0; i < nc; i++) if (AsmChart_in_layout(&run->charts[i]) && run->charts[i].component >= 0) {
        int32_t g = run->charts[i].component; ids[off[g]+fill[g]++] = (int32_t)i;
    }
    size_t *rels = ARENA_ALLOC(arena,(run->n_rels ? run->n_rels : 1)*sizeof *rels);
    FILE *fp = NULL;
    if (diag_dir) {
        char path[2048]; snprintf(path,sizeof path,"%s/stage5_uv_repair.csv",diag_dir); fp = fopen(path,"wb");
        if (fp) fputs("component,charts,active_vertices,pairs,iterations,accepted,rms_before,rms_after,max_displacement,prepare_sec,solve_sec,checks_sec,audit_sec\n",fp);
    }
    size_t accepted = 0, attempted = 0;
    for (size_t g = 0; g < ng; g++) {
        if (off[g+1]-off[g] < 2) continue;
        size_t nr = 0;
        for (size_t j = 0; j < run->n_rels; j++) {
            const AsmRelation *r = &run->rels[j];
            const AsmChart *a = &run->charts[r->a], *b = &run->charts[r->b];
            if (!AsmChart_in_layout(a) || !AsmChart_in_layout(b) || a->component != (int32_t)g || b->component != (int32_t)g) continue;
            double rms; size_t support; AsmContinuity_measure(run,r,&rms,&support);
            /* Candidate eligibility is not certification. As in the fixed
             * boundary trial, let a nearby failing SOURCE seam reach the
             * bounded solve. Eight-voxel disks, preserved passing seams,
             * shape and geometry still decide whether any edit is kept. */
            if (support >= 6 && isfinite(rms) && rms <= 2*ASM_CONTINUITY_GATE(r->rms)) rels[nr++] = j;
        }
        if (!nr) continue;
        Arena_Mark mark = Arena_save(arena); attempted++;
        accepted += act_uv_component(run,arena,ids+off[g],off[g+1]-off[g],rels,nr,audit,band,fp,NULL);
        Arena_restore(arena,mark);
    }
    if (fp) fclose(fp);
    fprintf(stderr,"[assemble continuity] local UV repair: %zu/%zu components accepted; XYZ, faces and measured trim gaps preserved\n",accepted,attempted);
    Arena_dispose(&arena); return accepted;
}

typedef struct ActUvWork { size_t rel; double error; } ActUvWork;

static int act_uv_work_compare(const void *p, const void *q)
{
    const ActUvWork *a = p, *b = q;
    if (a->error != b->error) return a->error > b->error ? -1 : 1;
    return (a->rel > b->rel)-(a->rel < b->rel);
}

static int act_uv_is_placed(const AsmChart *c)
{
    return AsmChart_in_layout(c) && c->placed && c->placement_state != ASM_PLACE_NONE;
}

static void act_uv_bounds(const AsmChart *c, double b[4])
{
    b[0] = b[1] = 1e300; b[2] = b[3] = -1e300;
    for (size_t v = 0; v < c->nv; v++) {
        double p[2]; act_point(c,(int32_t)v,NULL,p);
        for (int k = 0; k < 2; k++) { b[k] = fmin(b[k],p[k]); b[k+2] = fmax(b[k+2],p[k]); }
    }
}

static int act_uv_bounds_near(const double a[4], const double b[4])
{
    /* A chart previously polished in this pass can move back across the
     * entire diameter of its original eight-voxel disk: 16, not 8. Every
     * changed 24-voxel rim contact is therefore inside this neighbourhood. */
    const double reach = 24.0+16.0+1e-4;
    return a[0] <= b[2]+reach && b[0] <= a[2]+reach &&
           a[1] <= b[3]+reach && b[1] <= a[3]+reach;
}

static size_t act_relax_placed_uv(AsmRun *run, AsmContinuityUvAudit audit, double band, const char *diag_dir,
                                  size_t patch_vertices)
{
    if (!audit || !run->continuity_ready || run->n_charts < 2) return 0;
    Arena_T arena = Arena_new(); size_t nc = run->n_charts;
    AsmChart *origin = ARENA_ALLOC(arena,(nc ? nc : 1)*sizeof *origin);
    memcpy(origin,run->charts,nc*sizeof *origin);
    uint8_t *placed = ARENA_CALLOC(arena,nc ? nc : 1,1), *copied = ARENA_CALLOC(arena,nc ? nc : 1,1);
    double *bounds = ARENA_ALLOC(arena,(nc ? nc : 1)*4*sizeof *bounds);
    for (size_t i = 0; i < nc; i++) {
        placed[i] = (uint8_t)act_uv_is_placed(&run->charts[i]);
        if (placed[i]) act_uv_bounds(&run->charts[i],bounds+4*i);
    }
    ActUvWork *work = ARENA_ALLOC(arena,(run->n_rels ? run->n_rels : 1)*sizeof *work); size_t nw = 0;
    size_t *off = ARENA_CALLOC(arena,nc+1,sizeof *off);
    for (size_t j = 0; j < run->n_rels; j++) {
        const AsmRelation *r = &run->rels[j];
        if (!placed[r->a] || !placed[r->b]) continue;
        double rms; size_t support; AsmContinuity_measure(run,r,&rms,&support);
        double gate = ASM_CONTINUITY_GATE(r->rms);
        if (support < 6 || !isfinite(rms) || rms > 2*gate) continue;
        work[nw++] = (ActUvWork){j,rms/gate}; off[r->a+1]++; off[r->b+1]++;
    }
    qsort(work,nw,sizeof *work,act_uv_work_compare);
    for (size_t i = 0; i < nc; i++) off[i+1] += off[i];
    size_t *adj = ARENA_ALLOC(arena,(2*nw ? 2*nw : 1)*sizeof *adj), *fill = ARENA_ALLOC(arena,(nc ? nc : 1)*sizeof *fill);
    memcpy(fill,off,nc*sizeof *fill);
    for (size_t j = 0; j < nw; j++) {
        const AsmRelation *r = &run->rels[work[j].rel]; adj[fill[r->a]++] = j; adj[fill[r->b]++] = j;
    }
    uint8_t *done = ARENA_CALLOC(arena,nw ? nw : 1,1), *member = ARENA_CALLOC(arena,nc ? nc : 1,1);
    int32_t *ids = ARENA_ALLOC(arena,(nc ? nc : 1)*sizeof *ids);
    FILE *fp = NULL;
    if (diag_dir) {
        char path[2048]; snprintf(path,sizeof path,"%s/stage5_placed_uv_repair.csv",diag_dir); fp = fopen(path,"wb");
        if (fp) fputs("component,charts,active_vertices,pairs,iterations,accepted,rms_before,rms_after,max_displacement,prepare_sec,solve_sec,checks_sec,audit_sec\n",fp);
    }
    const size_t core_vertices = patch_vertices-patch_vertices/4;
    size_t attempted = 0, accepted = 0;
    for (size_t seed = 0; seed < nw; seed++) if (!done[seed]) {
        memset(member,0,nc); size_t n = 0, nv = 0;
        const AsmRelation *start = &run->rels[work[seed].rel];
        int32_t ends[2] = {start->a,start->b};
        for (int k = 0; k < 2; k++) if (!member[ends[k]]) {
            int32_t id = ends[k]; ids[n++] = id; member[id] = 1; nv += run->charts[id].nv;
        }
        /* Limit the core to one source hop from the seed, then add a
         * one-ring halo. Distant strained faces must not scale an unrelated
         * seam's whole line search. Every unvisited edge gets its own seed
         * later, so both of its charts eventually participate together. */
        size_t seed_charts = n;
        for (size_t at = 0; at < seed_charts; at++) for (size_t j = off[ids[at]]; j < off[ids[at]+1]; j++) {
            const AsmRelation *r = &run->rels[work[adj[j]].rel]; int32_t id = r->a == ids[at] ? r->b : r->a;
            if (member[id] || copied[id] || nv > core_vertices || run->charts[id].nv > core_vertices-nv) continue;
            ids[n++] = id; member[id] = 1; nv += run->charts[id].nv;
        }
        size_t core = n;
        for (size_t at = 0; at < core; at++) for (size_t j = off[ids[at]]; j < off[ids[at]+1]; j++) {
            const AsmRelation *r = &run->rels[work[adj[j]].rel]; int32_t id = r->a == ids[at] ? r->b : r->a;
            if (member[id] || nv > patch_vertices || run->charts[id].nv > patch_vertices-nv) continue;
            ids[n++] = id; member[id] = 1; nv += run->charts[id].nv;
        }
        for (size_t i = 0; i < n; i++) if (!copied[ids[i]]) {
            int32_t id = ids[i]; origin[id].uv = ARENA_ALLOC(arena,origin[id].nv*2*sizeof(float));
            memcpy(origin[id].uv,run->charts[id].uv,origin[id].nv*2*sizeof(float)); copied[id] = 1;
        }
        Arena_Mark mark = Arena_save(arena);
        AsmChart *checkpoint = ARENA_ALLOC(arena,n*sizeof *checkpoint), *guard = ARENA_ALLOC(arena,(nc ? nc : 1)*sizeof *guard);
        size_t ng = 0;
        double box[4] = {1e300,1e300,-1e300,-1e300};
        for (size_t i = 0; i < n; i++) {
            checkpoint[i] = origin[ids[i]];
            for (int k = 0; k < 2; k++) { box[k] = fmin(box[k],bounds[4*(size_t)ids[i]+k]); box[k+2] = fmax(box[k+2],bounds[4*(size_t)ids[i]+k+2]); }
        }
        for (size_t i = 0; i < nc; i++) if (placed[i] && (member[i] || act_uv_bounds_near(box,bounds+4*i))) guard[ng++] = run->charts[i];
        size_t *rels = ARENA_ALLOC(arena,(run->n_rels ? run->n_rels : 1)*sizeof *rels), nr = 0;
        for (size_t j = 0; j < run->n_rels; j++) {
            const AsmRelation *r = &run->rels[j];
            if ((!member[r->a] && !member[r->b]) || !placed[r->a] || !placed[r->b]) continue;
            double rms; size_t support; AsmContinuity_measure(run,r,&rms,&support);
            if (support >= 6 && isfinite(rms) && rms <= 2*ASM_CONTINUITY_GATE(r->rms)) rels[nr++] = j;
        }
        attempted++;
        int kept = act_uv_component_guarded(run,arena,ids,n,rels,nr,audit,band,fp,checkpoint,guard,ng,0);
        accepted += kept;
        fprintf(stderr,"  [continuity placed UV] group %zu, %zu charts / %zu vertices, %zu guard charts, %zu source constraints: %s\n",
                attempted,n,nv,ng,nr,kept ? "accepted" : "unchanged");
        if (kept) for (size_t i = 0; i < n; i++) act_uv_bounds(&run->charts[ids[i]],bounds+4*(size_t)ids[i]);
        for (size_t j = 0; j < nw; j++) {
            const AsmRelation *r = &run->rels[work[j].rel]; if (member[r->a] && member[r->b]) done[j] = 1;
        }
        Arena_restore(arena,mark);
    }
    if (fp) fclose(fp);
    double phase_displacement = 0;
    for (size_t i = 0; i < nc; i++) if (copied[i]) for (size_t v = 0; v < origin[i].nv; v++) {
        double du = (double)run->charts[i].uv[2*v]-origin[i].uv[2*v];
        double dv = (double)run->charts[i].uv[2*v+1]-origin[i].uv[2*v+1];
        phase_displacement = fmax(phase_displacement,hypot(du,dv));
    }
    fprintf(stderr,"[assemble continuity] placed UV repair: %zu/%zu groups accepted; maximum total correction %.9g vox; original source evidence and per-chart movement disks preserved\n",accepted,attempted,phase_displacement);
    Arena_dispose(&arena); return accepted;
}

size_t AsmContinuity_relax_placed_uv(AsmRun *run, AsmContinuityUvAudit audit, double band, const char *diag_dir)
{
    return act_relax_placed_uv(run,audit,band,diag_dir,524288);
}

size_t AsmContinuity_components(AsmRun *run)
{
    Arena_T arena = run->arena; Arena_Mark mark = Arena_save(arena);
    size_t nc = run->n_charts, deferred = 0, incumbent = 0;
    UnionFind uf = UF_new(arena, (int32_t)nc);
    for (size_t i = 0; i < run->n_rels; i++) {
        AsmRelation *r = &run->rels[i];
        if (!AsmRel_is_join(r) && !(r->continuity & ASM_CONT_INCUMBENT)) continue;
        if (!AsmChart_in_layout(&run->charts[r->a]) || !AsmChart_in_layout(&run->charts[r->b])) continue;
        double rms; size_t support;
        if (AsmContinuity_measure(run, r, &rms, &support)) {
            r->continuity |= ASM_CONT_VERIFIED; r->continuity_rms = rms;
            uf_union(&uf, r->a, r->b);
        } else if (r->continuity & ASM_CONT_INCUMBENT) {
            double current;
            if (act_incumbent_metric(run,r,&current) && current <= r->continuity_rms+1e-4) {
                uf_union(&uf,r->a,r->b); incumbent++;
            } else deferred++;
            /* Bundle membership only: material identity remains unresolved. */
            r->flags |= ASM_REL_SWITCHED;
            r->continuity &= ~(uint32_t)ASM_CONT_VERIFIED;
        } else {
            /* An uncertain pose-graph edge remains a recoverable source
             * hypothesis. Its two sides must not enter the initial frame
             * as trusted anchors merely because the pose solver grouped them. */
            r->flags |= ASM_REL_SWITCHED; r->continuity &= ~(uint32_t)ASM_CONT_VERIFIED; deferred++;
        }
    }
    int32_t *ids = ARENA_ALLOC(arena, (nc ? nc : 1)*sizeof(int32_t));
    for (size_t i = 0; i < nc; i++) ids[i] = -1;
    int32_t count = 0;
    for (size_t i = 0; i < nc; i++) if (AsmChart_in_layout(&run->charts[i])) {
        int32_t root = uf_find(&uf, (int32_t)i);
        if (ids[root] < 0) ids[root] = count++;
        run->charts[i].component = ids[root];
    }
    fprintf(stderr, "[assemble continuity] %d placement bundles; %zu sparse incumbent seams retained without source certification; %zu pose-graph edges deferred for recovery\n", count, incumbent, deferred);
    Arena_restore(arena, mark); return deferred;
}

size_t AsmContinuity_confirm(AsmRun *run)
{
    size_t added = 0;
    for (size_t i = 0; i < run->n_rels; i++) {
        AsmRelation *r = &run->rels[i];
        const AsmChart *a = &run->charts[r->a], *b = &run->charts[r->b];
        if (!AsmChart_in_layout(a) || !AsmChart_in_layout(b) || a->placement_state == ASM_PLACE_NONE || b->placement_state == ASM_PLACE_NONE || a->component != b->component) continue;
        double rms; size_t n;
        if (!AsmContinuity_measure(run, r, &rms, &n)) continue;
        if (!(r->continuity & ASM_CONT_VERIFIED)) added++;
        r->continuity |= ASM_CONT_VERIFIED; r->continuity_rms = rms;
        /* Confirming a join UNDOES the gate that rejected this relation, so it
         * has to clear the bar the audit will hold the seam to (AsmField_seam:
         * uses the shared seam bound). Sharing a component is not that bar. Without this a
         * seam the isometry gate threw out at 11.1 vox came back as a join and
         * the lattice drew quads across it: on the familiar 4x5x5 six such
         * promoted pairs carried all 73 of the emitted long edges and the one
         * cross-wrap step, which is the whole gate difference against the
         * accepted sheet's ancestor (which readmitted nothing).  A relation
         * that misses the bar keeps its measurement and stays placement
         * evidence. */
        if (!AsmRel_is_join(r) && rms <= ASM_CONTINUITY_CONFIRM_RMS) {
            r->flags &= ~(uint32_t)(ASM_REL_SWITCHED | ASM_REL_WEAK | ASM_REL_SHORT);
            r->flags |= ASM_REL_READMIT;
        }
    }
    return added;
}

typedef struct ActEdge { int32_t a, b; uint32_t good, bad; } ActEdge;
typedef struct ActEdgeSlot { uint64_t key; uint32_t good, bad; } ActEdgeSlot;
static int act_cmp_edge(const void *p, const void *q)
{
    const ActEdge *a = p, *b = q;
    if (a->a != b->a) return a->a < b->a ? -1 : 1;
    return (a->b > b->b) - (a->b < b->b);
}

void AsmContinuity_audit(const AsmRun *run, int32_t primary, const int32_t *chart, const float *xyz, size_t nv,
                         const int32_t *faces, size_t nf, const char *path, AsmContinuityStats *st)
{
    memset(st, 0, sizeof *st);
    Arena_T arena = run->arena; Arena_Mark mark = Arena_save(arena);
    size_t nc = run->n_charts;
    uint8_t *emitted = ARENA_CALLOC(arena, nc ? nc : 1, 1);
    uint8_t *unsupported = ARENA_CALLOC(arena, nc ? nc : 1, 1);
    UnionFind uf = UF_new(arena, (int32_t)nc);
    /* Sparse chart pairs from actual emitted triangle edges, independent of
     * component/lineage labels. A metadata join with no face cannot pass. */
    size_t cap = 4096, unique = 0;
    ActEdgeSlot *table = calloc(cap, sizeof *table);
    if (!table) { Arena_restore(arena, mark); return; }
    for (size_t v = 0; v < nv; v++) if (chart[v] >= 0 && (size_t)chart[v] < nc) emitted[chart[v]] = 1;
    for (size_t f = 0; f < nf; f++) for (int k = 0; k < 3; k++) {
        int32_t x = faces[f*3+k], y = faces[f*3+(k+1)%3];
        if (x < 0 || y < 0 || (size_t)x >= nv || (size_t)y >= nv) continue;
        int32_t a = chart[x], b = chart[y];
        if (a == b || a < 0 || b < 0 || (size_t)a >= nc || (size_t)b >= nc) continue;
        double d2 = 0;
        for (int j = 0; j < 3; j++) { double d = xyz[(size_t)x*3+j]-xyz[(size_t)y*3+j]; d2 += d*d; }
        int good = isfinite(d2) && d2 <= 64.0;
        if (good) uf_union(&uf, a, b);
        if (a > b) { int32_t t = a; a = b; b = t; }
        if (2*(unique+1) >= cap) {
            size_t next_cap = cap*2; ActEdgeSlot *next = calloc(next_cap, sizeof *next);
            if (!next) { free(table); Arena_restore(arena, mark); return; }
            for (size_t j = 0; j < cap; j++) if (table[j].key) {
                size_t h = (size_t)((table[j].key*UINT64_C(11400714819323198485)) >> 16) & (next_cap-1);
                while (next[h].key) h = (h+1)&(next_cap-1);
                next[h] = table[j];
            }
            free(table); table = next; cap = next_cap;
        }
        uint64_t key = ((uint64_t)(a+1)<<32) | (uint32_t)(b+1);
        size_t h = (size_t)((key*UINT64_C(11400714819323198485)) >> 16) & (cap-1);
        while (table[h].key && table[h].key != key) h = (h+1)&(cap-1);
        if (!table[h].key) { table[h].key = key; unique++; }
        if (good) table[h].good = 1; else table[h].bad = 1;
    }
    ActEdge *edges = malloc((unique ? unique : 1)*sizeof *edges);
    if (!edges) { free(table); Arena_restore(arena, mark); return; }
    size_t ne = 0;
    for (size_t i = 0; i < cap; i++) if (table[i].key) {
        edges[ne++] = (ActEdge){ (int32_t)(table[i].key>>32)-1, (int32_t)(uint32_t)table[i].key-1, table[i].good, table[i].bad };
        st->wrong_placements += table[i].bad != 0;
    }
    free(table); qsort(edges, ne, sizeof *edges, act_cmp_edge);
    FILE *fp = path ? fopen(path, "wb") : NULL;
    if (fp) fputs("kind,chart_a,chart_b,support,residual,area\n", fp);
    if (fp) for (size_t i = 0; i < unique; i++) if (edges[i].bad)
        fprintf(fp, "emitted_long_edge,%d,%d,0,8,0\n", edges[i].a, edges[i].b);
    for (size_t i = 0; i < nc; i++) {
        const AsmChart *c = &run->charts[i];
        if (!emitted[i]) continue;
        if (c->placement_state == ASM_PLACE_NONE) {
            st->unsupported_charts++; st->unsupported_area += c->area3d;
            unsupported[uf_find(&uf, (int32_t)i)] = 1;
            if (fp) fprintf(fp, "unsupported,%zu,-1,0,0,%.6e\n", i, c->area3d);
        }
    }
    for (size_t i = 0; i < nc; i++) st->unsupported_islands += unsupported[i] != 0;
    for (size_t i = 0; i < run->n_rels; i++) {
        const AsmRelation *r = &run->rels[i];
        const AsmChart *a = &run->charts[r->a], *b = &run->charts[r->b];
        int ea = emitted[r->a], eb = emitted[r->b];
        /* The retained source is the clean-stage contract. Recorded relations
         * to material already excluded upstream do not describe a hole made
         * by the assembly pass. */
        if (!AsmChart_in_layout(a) || !AsmChart_in_layout(b)) continue;
        if (!ea && !eb) continue;
        if (r->flags & (ASM_REL_DROPPED | ASM_REL_CONTACT)) continue;
        ActEdge key = { r->a < r->b ? r->a : r->b, r->a < r->b ? r->b : r->a };
        const ActEdge *edge = bsearch(&key, edges, unique, sizeof *edges, act_cmp_edge);
        int linked = ea && eb && edge && edge->good && !edge->bad;
        double rms = 1e300; size_t support = 0;
        int measured = ea && eb && AsmContinuity_measure(run, r, &rms, &support);
        if (!(ea && eb)) { rms = 1e300; support = 0; }
        const char *kind = NULL;
        if (AsmRel_is_join(r)) {
            if (!linked) { st->metadata_only_joins++; kind = "metadata_only"; }
            else if ((r->continuity & ASM_CONT_SOURCE) && !measured) { st->wrong_placements++; kind = "wrong_placement"; }
            else st->verified_joins++;
        }
        if (r->continuity & ASM_CONT_SOURCE) {
            int32_t last = -1;
            for (int32_t k = 0; k < r->corr_count; k++) {
                int32_t id = run->corr[(size_t)r->corr_first+k].run;
                if (id >= 0 && id != last) { st->source_runs++; last = id; }
            }
            if (!linked || !measured) { st->unexplained_breaks++; if (!kind) kind = "source_break"; }
        } else if (!AsmRel_is_join(r)) { st->unresolved_hypotheses++; if (!kind) kind = "unresolved"; }
        if (kind && fp) fprintf(fp, "%s,%d,%d,%zu,%.6g,%.6e\n", kind, r->a, r->b, support, rms, (!ea ? a->area3d : 0) + (!eb ? b->area3d : 0));
    }
    if (fp) fclose(fp);
    st->audit_complete = primary >= 0 && nv > 0 && nf > 0;
    st->confetti_free = st->audit_complete && !st->unsupported_islands && !st->unexplained_breaks && !st->metadata_only_joins && !st->wrong_placements && !st->unresolved_hypotheses;
    free(edges); Arena_restore(arena, mark);
}

static int act_uv_reject_test(AsmRun *run, const AsmChart *saved, size_t n, AsmContinuityAuditCache *cache)
{
    (void)run; (void)saved; (void)n; (void)cache; return 0;
}

static int act_registration_accept_test(AsmRun *run, const AsmChart *saved, size_t n, AsmContinuityAuditCache *cache)
{
    (void)run; (void)saved; (void)n; (void)cache; return 1;
}

static int act_registration_source_selftest(void)
{
    int fails = 0;
    for (int mirror = 0; mirror < 2; mirror++) for (int inconsistent = 0; inconsistent < 2; inconsistent++) {
        /* Three measured seams start inside the ordinary residual gate.
         * Inconsistent flattened frames make the two dense seams pull the
         * sparse seam apart. Approval of the separate geometry callback
         * must not authorize breaking a previously passing source seam. */
        AsmRun run = {0}; run.arena = Arena_new();
        float uv[60], xyz[90];
        for (int v = 0; v < 30; v++) {
            uv[2*v] = xyz[3*v] = (float)(2*(v%5)-4);
            uv[2*v+1] = xyz[3*v+1] = (float)(2*(v/5)-5); xyz[3*v+2] = 0;
        }
        for (int i = 0; i < 3; i++) {
            AsmChart c = {0}; c.id = i; c.nv = 30; c.uv = uv; c.xyz = xyz;
            c.area3d = i ? 1000 : 1001; c.flags = mirror ? ASM_CHART_MIRROR : 0;
            AsmRun_push_chart(&run,&c);
        }
        const int sparse[6] = {0,5,10,19,24,29};
        for (int i = 0; i < 3; i++) {
            AsmRelation r = {0}; r.a = i == 2 ? 1 : 0; r.b = i == 0 ? 1 : 2;
            r.corr_first = (int32_t)run.n_corr; r.corr_count = i ? 30 : 6;
            r.n_corr = (size_t)r.corr_count; r.continuity = ASM_CONT_SOURCE;
            for (int j = 0; j < r.corr_count; j++) {
                AsmCorr c = {0}; c.va = c.vb = i ? j : sparse[j]; c.valid = 3;
                c.gap_a[0] = c.gap_b[0] = i ? (i == 2 && inconsistent ? -7.5f : 7.5f) : 0;
                AsmRun_push_corr(&run,&c);
            }
            AsmRun_push_rel(&run,&r);
        }
        AsmChart saved[3]; memcpy(saved,run.charts,sizeof saved);
        AsmRelation saved_rel[3]; memcpy(saved_rel,run.rels,sizeof saved_rel);
        AsmCorr saved_corr[66]; memcpy(saved_corr,run.corr,sizeof saved_corr);
        int ok = 1; double rms = 0; size_t support = 0;
        for (int i = 0; i < 3; i++) ok = ok && AsmContinuity_measure(&run,&run.rels[i],&rms,&support);
        size_t kept = AsmContinuity_register_components(&run,act_registration_accept_test);
        size_t broken = 0; double worst = 0;
        for (int i = 0; i < 3; i++) {
            broken += !AsmContinuity_measure(&run,&run.rels[i],&rms,&support); worst = fmax(worst,rms);
        }
        ok = ok && kept == (size_t)!inconsistent && !broken &&
            !memcmp(run.rels,saved_rel,sizeof saved_rel) && !memcmp(run.corr,saved_corr,sizeof saved_corr);
        if (inconsistent) ok = ok && !memcmp(run.charts,saved,sizeof saved);
        else ok = ok && worst < 1e-5 && !memcmp(&run.charts[0],&saved[0],sizeof saved[0]);
        if (!ok) {
            fprintf(stderr,"  continuity registration selftest FAIL: mirror %d inconsistent %d, kept %zu, broken %zu, worst %.9g\n",
                    mirror,inconsistent,kept,broken,worst); fails++;
        }
        Arena_dispose(&run.arena);
    }
    return fails;
}

static int act_registration_hypothesis_reject_test(AsmRun *run, const AsmChart *saved, size_t n, AsmContinuityAuditCache *cache)
{
    (void)run; (void)saved; (void)n; return cache->guard_charts == NULL;
}

static int act_pose_plane_feasibility_selftest(void)
{
    int fails = 0;
    for (int sign = -1; sign <= 1; sign += 2) {
        float uv[] = {0,0}; AsmChart charts[2] = {{0}};
        charts[0].id = 0; charts[1].id = 1; charts[0].uv = charts[1].uv = uv;
        charts[0].nv = charts[1].nv = 1; charts[1].pose_x = sign*10;
        AsmRun run = {0}; run.charts = charts; run.n_charts = 2;
        int32_t index[] = {-1,0};
        ActPosePlane plane = {{0,0,1,0},(double)sign,0,10};
        double step[] = {-sign*5e-10,0,0};
        /* The line search's former 1e-9 allowance accepted this step,
         * then a positive linear RHS made the next zero-step fallback
         * ineligible. This is the same mismatch as real pair76:835/277:3322. */
        int ok = !act_pose_planes_step(&run,index,&plane,1,step,1);
        step[0] = 0; ok = ok && act_pose_planes_step(&run,index,&plane,1,step,1);
        step[0] = sign*1e-9; ok = ok && act_pose_planes_step(&run,index,&plane,1,step,1);
        charts[1].pose_x = sign*(10-5e-10);
        ActHalfspace row; act_pose_plane_row(&run,index,&plane,&row);
        ok = ok && row.rhs > 0 && row.rhs < 1e-9;
        /* Identical boundary poses must give zero slack in both checks.
         * (distance-x)-y can be positive when distance-(x+y) is exactly
         * zero, even though the nonlinear line search accepts that pose. */
        charts[1].pose_x = sign*1250; charts[1].pose_y = 1.0/6;
        plane.nx = sign*.8; plane.ny = .6;
        plane.distance = plane.nx*charts[1].pose_x+plane.ny*charts[1].pose_y;
        step[0] = 0;
        act_pose_plane_row(&run,index,&plane,&row);
        ok = ok && act_pose_planes_step(&run,index,&plane,1,step,1) && row.rhs == 0;
        if (!ok) { fprintf(stderr,"  continuity plane feasibility selftest FAIL: sign %d\n",sign); fails++; }
    }
    return fails;
}

static int act_pose_plane_curvature_selftest(void)
{
    int fails = 0;
    for (int sign = -1; sign <= 1; sign += 2) {
        float anchor_uv[] = {0,0}, moving_uv[] = {10,0};
        AsmChart charts[2] = {{0}};
        charts[0].id = 0; charts[0].uv = anchor_uv; charts[0].nv = 1;
        charts[1].id = 1; charts[1].uv = moving_uv; charts[1].nv = 1;
        if (sign < 0) charts[1].flags = ASM_CHART_MIRROR;
        AsmChart original[2]; memcpy(original,charts,sizeof original);
        AsmRun run = {0}; run.arena = Arena_new(); run.charts = charts; run.n_charts = 2;
        int32_t index[] = {-1,0}; ActPosePlane plane = {{0,0,1,0},(double)sign,0,10};
        int rows[] = {0,1,2}, cols[] = {0,1,2}; double vals[] = {1,1,1}, rhs[] = {0,0,.1};
        ActPoseQuadratic q = {3,3,rows,cols,vals,rhs};
        double result[] = {1234567,0,0,.1,7654321};
        /* The tangent rotation satisfies the linear row exactly but every
         * positive line-search scale moves 10*cos(theta) inside x >= 10.
         * A translation of 10*(1-cos(theta)) restores exact feasibility
         * while retaining the useful full rotation and reducing energy. */
        ActHalfspace row; act_pose_plane_row(&run,index,&plane,&row);
        double tangent = 0;
        for (int k = 0; k < row.count; k++) tangent += row.values[k]*result[1+row.ids[k]];
        int ok = tangent == row.rhs && !act_pose_planes_step(&run,index,&plane,1,result+1,1);
        int solved = act_pose_planes_direction(&run,index,&plane,1,&q,result+1);
        double energy = .5*(result[1]*result[1]+result[2]*result[2]+result[3]*result[3])-.1*result[3];
        ok = ok && solved && act_pose_planes_step(&run,index,&plane,1,result+1,1) &&
             fabs(result[1]-sign*10*(1-cos(.1))) < 1e-6 && result[2] == 0 && result[3] == .1 && energy < 0 &&
             result[0] == 1234567 && result[4] == 7654321 && !memcmp(charts,original,sizeof original);
        /* Feasible geometry must not license an uphill quadratic step,
         * and an impossible translation interval must leave the input exact. */
        double unchanged[] = {1234567,0,0,.1,7654321};
        memcpy(result,unchanged,sizeof result); rhs[2] = 0;
        ok = ok && !act_pose_planes_direction(&run,index,&plane,1,&q,result+1) && !memcmp(result,unchanged,sizeof result);
        rhs[2] = .1;
        ActPosePlane impossible[] = {plane,{{0,0,1,0},(double)-sign,0,-9}};
        double correction = -1;
        ok = ok && !act_pose_planes_project(&run,index,impossible,2,3,result+1,&correction) &&
             !memcmp(result,unchanged,sizeof result) && correction == -1;
        /* A longer lever makes the full translation correction too costly.
         * Backtracking must use the original direction, retain the rotation
         * of the chosen scale, and reduce the actual quadratic. */
        moving_uv[0] = 100; plane.distance = 100;
        solved = act_pose_planes_direction(&run,index,&plane,1,&q,result+1);
        energy = .5*(result[1]*result[1]+result[2]*result[2]+result[3]*result[3])-.1*result[3];
        ok = ok && solved && result[3] == .025 && result[2] == 0 && energy < 0 &&
             fabs(result[1]-sign*100*(1-cos(.025))) < 2e-6 &&
             act_pose_planes_step(&run,index,&plane,1,result+1,1) &&
             result[0] == 1234567 && result[4] == 7654321 && !memcmp(charts,original,sizeof original);
        if (!ok) { fprintf(stderr,"  continuity plane curvature selftest FAIL: sign %d, solved %d\n",sign,solved); fails++; }
        Arena_dispose(&run.arena);
    }
    return fails;
}

static int act_halfspace_selftest(void)
{
    Arena_T arena = Arena_new(); int fails = 0;
    int rows[] = {0,1,1}, cols[] = {0,0,1}; double vals[] = {2,1,2}, rhs[] = {-3,-3};
    const double expected[6][2] = {{-1,-1},{0,-1.5},{0,0},{1,1},{0,-1.5},{-1,-1}};
    for (int test = 0; test < 7; test++) {
        ActHalfspace c[2] = {{{0},1,{1},0},{{1},1,{1},0}}; size_t count = 0;
        if (test == 1) count = 1;
        if (test == 2) count = 2;
        if (test == 3) { c[0] = (ActHalfspace){{0,1},2,{1,1},2}; count = 1; }
        if (test == 4) { c[1] = c[0]; count = 2; }
        if (test == 5) { c[0].rhs = c[1].rhs = -2; count = 2; }
        if (test == 6) { c[0].rhs = 1; c[1] = (ActHalfspace){{0},1,{-1},0}; count = 2; }
        double result[] = {1234567,0,0,7654321};
        int solved = act_halfspace_solve(arena,2,3,rows,cols,vals,rhs,result+1,c,count,4096,0);
        int ok = solved == (test != 6) && result[0] == 1234567 && result[3] == 7654321;
        if (test != 6) ok = ok && hypot(result[1]-expected[test][0],result[2]-expected[test][1]) < 1e-6;
        solved = act_halfspace_solve(arena,2,3,rows,cols,vals,rhs,result+1,c,count,4096,2048);
        ok = ok && solved == (test != 6) && result[0] == 1234567 && result[3] == 7654321;
        if (test != 6) ok = ok && hypot(result[1]-expected[test][0],result[2]-expected[test][1]) < 1e-6;
        if (!ok) { fprintf(stderr,"  continuity halfspace selftest FAIL: case %d, solved %d, x %.9g %.9g\n",test,solved,result[1],result[2]); fails++; }
    }
    ActHalfspace c = {{0},1,{1},-1}; double b = -4, x = 0;
    int ok = act_halfspace_solve(arena,1,1,rows,cols,vals,&b,&x,&c,1,4096,0) && fabs(x+1) < 1e-8;
    c.ids[0] = 1;
    ok = ok && !act_halfspace_solve(arena,1,1,rows,cols,vals,&b,&x,&c,1,4096,0) &&
         !act_halfspace_solve(arena,0,1,rows,cols,vals,&b,&x,&c,1,4096,0);
    if (!ok) { fprintf(stderr,"  continuity halfspace selftest FAIL: scalar or bounds\n"); fails++; }
    /* One coordinate sweep is deliberately incomplete. Zero is feasible;
     * the fallback must satisfy both original inequalities and improve the
     * actual quadratic. A positive lower bound makes zero infeasible and
     * must prevent that fallback. */
    ActHalfspace limited[] = {{{0},1,{1},-.5},{{1},1,{1},-.5}};
    double result[] = {1234567,0,0,7654321};
    int solved = act_halfspace_solve(arena,2,3,rows,cols,vals,rhs,result+1,limited,2,1,0);
    double a = result[1], d = result[2], energy = a*a+a*d+d*d+3*(a+d);
    ok = solved == 2 && a >= -.5 && d >= -.5 && energy < 0 && result[0] == 1234567 && result[3] == 7654321;
    limited[0].rhs = .5;
    ok = ok && !act_halfspace_solve(arena,2,3,rows,cols,vals,rhs,result+1,limited,2,1,0);
    if (!ok) { fprintf(stderr,"  continuity halfspace selftest FAIL: feasible retreat\n"); fails++; }
    limited[0].rhs = -.5;
    solved = act_halfspace_solve(arena,2,3,rows,cols,vals,rhs,result+1,limited,2,1,2048);
    ok = solved == 1 && hypot(result[1]+.5,result[2]+.5) < 1e-7;
    /* The vertical wall blocks first, but its multiplier later turns
     * negative. Removing it reaches the analytic projection onto the
     * tilted wall. Include an exactly redundant copy of that wall. */
    int ir[] = {0,1}, ic[] = {0,1}; double iv[] = {1,1}, ib[] = {10,1};
    ActHalfspace walls[] = {{{0},1,{-1},-1},{{0,1},2,{-1,-.05},-1.02},{{0,1},2,{-1,-.05},-1.02}};
    solved = act_halfspace_solve(arena,2,2,ir,ic,iv,ib,result+1,walls,3,1,2048);
    double t = 9.03/1.0025;
    ok = ok && solved == 1 && hypot(result[1]-(10-t),result[2]-(1-.05*t)) < 1e-7 &&
         result[1] < 1 && result[0] == 1234567 && result[3] == 7654321;
    if (!ok) { fprintf(stderr,"  continuity halfspace selftest FAIL: active set or row removal\n"); fails++; }
    /* A binding row after the old 1024-row limit must still participate.
     * The earlier rows are loose duplicates; the coupled SPD projection
     * onto x >= -.5 is analytically (-.5,-1.25). */
    const size_t many_count = 1500;
    ActHalfspace *many = ARENA_ALLOC(arena,many_count*sizeof *many);
    for (size_t i = 0; i < many_count; i++) many[i] = (ActHalfspace){{0},1,{1},-2};
    many[many_count-1].rhs = -.5;
    solved = act_halfspace_solve(arena,2,3,rows,cols,vals,rhs,result+1,many,many_count,1,2048);
    ok = solved == 1 && hypot(result[1]+.5,result[2]+1.25) < 1e-7 &&
         result[0] == 1234567 && result[3] == 7654321 &&
         act_halfspace_capacity(0) == 0 && act_halfspace_capacity(1) == 8192 &&
         act_halfspace_capacity(4096) == 4096 && act_halfspace_capacity((1u << 24)+1) == 0;
    if (!ok) { fprintf(stderr,"  continuity halfspace selftest FAIL: many rows or memory bound\n"); fails++; }
    Arena_dispose(&arena); return fails;
}

static int act_registration_hypothesis_selftest(void)
{
    int fails = 0;
    for (int mirror = 0; mirror < 2; mirror++) for (int reject = 0; reject < 2; reject++) for (int distant = 0; distant < 2; distant++) {
        const int nc = 3+distant;
        AsmRun run = {0}; run.arena = Arena_new();
        float uv[4][12] = {{0}}, xyz[4][18] = {{0}}, original_uv[4][12], original_xyz[4][18];
        int32_t faces[] = {0,1,3,0,3,2,2,3,5,2,5,4};
        for (int i = 0; i < nc; i++) {
            for (int v = 0; v < 6; v++) {
                uv[i][2*v] = xyz[i][3*v] = (float)(2*(v%2)-1);
                uv[i][2*v+1] = xyz[i][3*v+1] = (float)(2*(v/2)-2); xyz[i][3*v+2] = 0;
            }
            AsmChart chart = {0}; chart.id = i; chart.uv = uv[i]; chart.xyz = xyz[i]; chart.nv = 6; chart.nf = 4; chart.faces = faces;
            chart.area3d = i ? 1000 : 1001; chart.pose_theta = .4; chart.flags = mirror ? ASM_CHART_MIRROR : 0;
            if (i == 3) { chart.pose_x = 1000*cos(.4); chart.pose_y = 1000*sin(.4); }
            AsmRun_push_chart(&run,&chart);
        }
        for (int i = 0; i < nc; i++) {
            AsmRelation rel = {0}; rel.a = i == 1 ? 1 : 0; rel.b = i == 0 ? 1 : i == 3 ? 3 : 2;
            rel.corr_first = (int32_t)run.n_corr; rel.corr_count = rel.n_corr = 6; rel.continuity = ASM_CONT_SOURCE;
            rel.flags = i == 2 ? ASM_REL_WEAK : 0;
            for (int v = 0; v < 6; v++) {
                AsmCorr corr = {0}; corr.va = corr.vb = v; corr.valid = 3;
                corr.gap_a[0] = corr.gap_b[0] = i == 2 ? 10 : 0; AsmRun_push_corr(&run,&corr);
            }
            AsmRun_push_rel(&run,&rel);
        }
        AsmChart original[4]; memcpy(original,run.charts,nc*sizeof *original);
        AsmRelation original_rel[4]; memcpy(original_rel,run.rels,nc*sizeof *original_rel);
        AsmCorr original_corr[24]; memcpy(original_corr,run.corr,6*nc*sizeof *original_corr);
        memcpy(original_uv,uv,sizeof uv); memcpy(original_xyz,xyz,sizeof xyz);
        double rms; size_t support;
        int ok = AsmContinuity_measure(&run,&run.rels[0],&rms,&support) && AsmContinuity_measure(&run,&run.rels[1],&rms,&support) &&
                 !AsmContinuity_measure(&run,&run.rels[2],&rms,&support);
        size_t ordinary = AsmContinuity_register_components(&run,act_registration_accept_test);
        size_t kept = AsmContinuity_register_failed(&run,original,reject ? act_registration_hypothesis_reject_test : act_registration_accept_test);
        ok = ok && ordinary == (size_t)!distant && kept == (size_t)!reject &&
             AsmContinuity_measure(&run,&run.rels[0],&rms,&support) && AsmContinuity_measure(&run,&run.rels[1],&rms,&support) &&
             AsmContinuity_measure(&run,&run.rels[2],&rms,&support) == !reject &&
             !memcmp(run.rels,original_rel,nc*sizeof *original_rel) && !memcmp(run.corr,original_corr,6*nc*sizeof *original_corr) &&
             !memcmp(uv,original_uv,sizeof uv) && !memcmp(xyz,original_xyz,sizeof xyz) && !memcmp(run.charts,original,sizeof original[0]);
        if (reject) ok = ok && !memcmp(run.charts,original,nc*sizeof *original);
        /* A distant failed kept edge must not override the exact fit mask
         * and force a thousand-voxel step. It stays in the full energy audit,
         * and the unsupported chart itself must remain exactly unchanged. */
        if (distant) ok = ok && !memcmp(run.charts+3,original+3,sizeof original[3]);
        /* A second trial cannot reset the original 256-voxel bound after
         * ordinary registration has already used 253 voxels of it. */
        memcpy(run.charts,original,nc*sizeof *original);
        run.charts[0].placement_state = ASM_PLACE_ROOT;
        AsmChart checkpoint[4];
        for (int i = 1; i < nc; i++) checkpoint[i-1] = run.charts[i];
        checkpoint[nc-1] = run.charts[0];
        for (int i = 0; i < nc-1; i++) {
            checkpoint[i].pose_x -= (mirror ? -253 : 253)*cos(.4);
            checkpoint[i].pose_y -= (mirror ? -253 : 253)*sin(.4);
        }
        AsmChart current[4]; memcpy(current,run.charts,nc*sizeof *current);
        ok = ok && !act_register_hypotheses(&run,checkpoint,(size_t)nc,act_registration_accept_test) && !memcmp(run.charts,current,nc*sizeof *current);
        if (!ok) { fprintf(stderr,"  continuity hypothesis registration selftest FAIL: mirror %d reject %d distant %d, kept %zu\n",mirror,reject,distant,kept); fails++; }
        Arena_dispose(&run.arena);
    }
    return fails;
}

static int act_trim_order_selftest(void)
{
    /* Reordering the same curved surface must not change its measured
     * trim. A zero-area incident face must not mask the valid neighbours. */
    enum { NV = 27, NF = 32 };
    AsmRun run = {0}; run.arena = Arena_new();
    float xyz[2][3*NV], uv[2][2*NV]; uint8_t boundary[NV];
    int32_t original[3*NF], faces[3*(NF+NV)];
    for (int j = 0; j < 9; j++) for (int u = 0; u < 3; u++) {
        int v = 3*j+u; boundary[v] = 1;
        for (int side = 0; side < 2; side++) {
            double q = 5*u+14*side, z = 4*j;
            xyz[side][3*v] = (float)z; xyz[side][3*v+1] = (float)q;
            xyz[side][3*v+2] = (float)(0.03*q*z);
            uv[side][2*v] = (float)(5*u); uv[side][2*v+1] = (float)z;
        }
        if (j < 8 && u < 2) {
            int32_t f[] = {v,v+1,v+4,v,v+4,v+3}; memcpy(original+6*(2*j+u),f,sizeof f);
        }
    }
    for (int side = 0; side < 2; side++) {
        AsmChart c = {0}; c.id = side; c.nv = NV; c.nf = NF;
        c.xyz = xyz[side]; c.uv = uv[side]; c.faces = faces; c.boundary = boundary;
        AsmRun_push_chart(&run,&c);
    }
    AsmRelation rel = {0}; rel.a = 0; rel.b = 1; rel.rms = 0.1; rel.corr_count = rel.n_corr = 9;
    for (int j = 0; j < 9; j++) { AsmCorr c = {3*j+2,3*j}; AsmRun_push_corr(&run,&c); }
    AsmRun_push_rel(&run,&rel);
    AsmCorr reference[9]; int fails = 0;
    for (int order = 0; order < 3; order++) {
        size_t prefix = order == 2 ? NV : 0;
        for (size_t i = 0; i < prefix; i++) faces[3*i] = faces[3*i+1] = faces[3*i+2] = (int32_t)i;
        for (int f = 0; f < NF; f++) memcpy(faces+3*(prefix+f),original+3*(order ? NF-1-f : f),3*sizeof *faces);
        run.charts[0].nf = run.charts[1].nf = NF+prefix;
        int ok = AsmContinuity_prepare(&run) == 0 && (run.rels[0].continuity & ASM_CONT_SOURCE);
        double maximum = 0;
        for (int j = 0; j < 9; j++) {
            const AsmCorr *c = &run.corr[j];
            ok = ok && c->valid == 3 && c->run == 0;
            if (order) for (int k = 0; k < 2; k++) {
                maximum = fmax(maximum,fabs(c->gap_a[k]-reference[j].gap_a[k]));
                maximum = fmax(maximum,fabs(c->gap_b[k]-reference[j].gap_b[k]));
            }
        }
        if (!order) memcpy(reference,run.corr,sizeof reference);
        if (!ok || maximum > 1e-5) {
            fprintf(stderr,"  continuity trim selftest FAIL: face order %d, valid %d, maximum change %.9g vox\n",order,ok,maximum); fails++;
        }
    }
    Arena_dispose(&run.arena); return fails;
}

static int act_uv_grid_solver_selftest(void)
{
    /* Cycles make the no-fill factor approximate. Compare the actual solve
     * to a direct reference on a weighted 2-D mesh, including fixed band
     * vertices and a seam linking two distant boundary vertices. */
    enum { SIDE = 8, N = SIDE*SIDE };
    Arena_T arena = Arena_new();
    ActSpring spring[4*N]; size_t ns = 0;
    ActMatch match[] = {act_match_pair(0,-1,10,-7),act_match_pair(N-1,-1,3,8),act_match_pair(7,56,2,-4),{0}};
    match[3].target = (ActVec){3,-9}; match[3].weight = .5;
    act_match_add(&match[3],5,-1.25); act_match_add(&match[3],6,1.75);
    act_match_add(&match[3],33,.5); act_match_add(&match[3],57,-1);
    for (int y = 0; y < SIDE; y++) for (int u = 0; u < SIDE; u++) {
        int i = y*SIDE+u;
        if (u+1 < SIDE) spring[ns++] = (ActSpring){i,i+1,1,1};
        if (y+1 < SIDE) spring[ns++] = (ActSpring){i+SIDE,i,0.25f+0.5f*(i%5),1};
        if (u+1 < SIDE && y+1 < SIDE) spring[ns++] = (ActSpring){i,i+SIDE+1,0.125f,1.414f};
    }
    spring[ns++] = (ActSpring){15,-1,2,1};
    spring[ns++] = (ActSpring){-1,42,0.5f,1};
    int rows[6*N], cols[6*N], nt = 0;
    double vals[6*N], diag[N], rhs[2*N] = {0}, ref[2*N];
    for (int i = 0; i < N; i++) diag[i] = 1e-4;
    for (size_t i = 0; i < ns; i++) {
        int a = spring[i].a, b = spring[i].b;
        double w = spring[i].weight;
        if (a >= 0) diag[a] += w;
        if (b >= 0) diag[b] += w;
        if (a >= 0 && b >= 0) { rows[nt] = a > b ? a : b; cols[nt] = a < b ? a : b; vals[nt++] = -w; }
    }
    for (int i = 0; i < 4; i++) for (int k = 0; k < match[i].count; k++) {
        int a = match[i].id[k]; double ca = match[i].coefficient[k], w = match[i].weight;
        if (a < 0) continue;
        diag[a] += w*ca*ca; rhs[a] += w*ca*match[i].target.x; rhs[N+a] += w*ca*match[i].target.y;
        for (int j = 0; j < k; j++) {
            int b = match[i].id[j]; if (b < 0) continue;
            rows[nt] = a > b ? a : b; cols[nt] = a < b ? a : b;
            vals[nt++] = w*ca*match[i].coefficient[j];
        }
    }
    for (int i = 0; i < N; i++) { rows[nt] = cols[nt] = i; vals[nt++] = diag[i]; }
    SparseFactor_T factor = NULL;
    int ok = Sparse_factor_spd(N,nt,rows,cols,vals,&factor) == 0;
    if (ok) ok = Sparse_factor_solve_multi(factor,rhs,ref,2) == 0;
    Sparse_factor_free(&factor);
    ActVec x[N] = {{0}}, ax[N]; int iterations = 0;
    ok = ok && act_uv_solve(arena,spring,ns,match,4,x,N,&iterations);
    double error = 0, residual = 0, initial = 0;
    if (ok) {
        act_uv_operator(spring,ns,match,4,x,ax,N);
        for (int i = 0; i < N; i++) {
            error = fmax(error,hypot(x[i].x-ref[i],x[i].y-ref[N+i]));
            residual += ((ax[i].x-rhs[i])*(ax[i].x-rhs[i])+(ax[i].y-rhs[N+i])*(ax[i].y-rhs[N+i]))/diag[i];
            initial += (rhs[i]*rhs[i]+rhs[N+i]*rhs[N+i])/diag[i];
        }
        ok = error < 1e-4 && residual <= initial*1e-12;
    }
    if (!ok) fprintf(stderr,"  continuity UV solver selftest FAIL: mesh, %d iterations, error %.9g, residual %.9g\n",iterations,error,residual);
    Arena_dispose(&arena); return !ok;
}

static int act_uv_solver_selftest(void)
{
    /* A long, narrow band has a slow translation mode. It must solve the
     * same equations as the direct reference, even with repeated mesh edges,
     * oppositely directed matches, and both coordinates active. */
    enum { N = 2048 };
    Arena_T arena = Arena_new();
    ActSpring *spring = ARENA_ALLOC(arena,2*(N-1)*sizeof *spring);
    ActMatch match[] = {act_match_pair(0,-1,10,-7),act_match_pair(-1,N-1,-3,-8)};
    ActVec *x = ARENA_CALLOC(arena,N,sizeof *x), *ax = ARENA_ALLOC(arena,N*sizeof *ax);
    int *rows = ARENA_ALLOC(arena,2*N*sizeof *rows), *cols = ARENA_ALLOC(arena,2*N*sizeof *cols);
    double *vals = ARENA_ALLOC(arena,2*N*sizeof *vals), *rhs = ARENA_CALLOC(arena,2*N,sizeof *rhs);
    double *ref = ARENA_ALLOC(arena,2*N*sizeof *ref);
    int nt = 0;
    for (int i = 0; i < N; i++) {
        rows[nt] = cols[nt] = i; vals[nt++] = 2.0001;
        if (i) {
            spring[2*(i-1)] = (ActSpring){i-1,i,0.25f,1};
            spring[2*(i-1)+1] = (ActSpring){i,i-1,0.75f,1};
            rows[nt] = i; cols[nt] = i-1; vals[nt++] = -1;
        }
    }
    rhs[0] = 10; rhs[N-1] = 3; rhs[N] = -7; rhs[2*N-1] = 8;
    SparseFactor_T factor = NULL;
    int ok = Sparse_factor_spd(N,nt,rows,cols,vals,&factor) == 0;
    if (ok) ok = Sparse_factor_solve_multi(factor,rhs,ref,2) == 0;
    Sparse_factor_free(&factor);
    int iterations = 0;
    ok = ok && act_uv_solve(arena,spring,2*(N-1),match,2,x,N,&iterations);
    double error = 0, residual = 0;
    if (ok) {
        act_uv_operator(spring,2*(N-1),match,2,x,ax,N);
        for (int i = 0; i < N; i++) {
            error = fmax(error,hypot(x[i].x-ref[i],x[i].y-ref[N+i]));
            residual += (ax[i].x-rhs[i])*(ax[i].x-rhs[i])+(ax[i].y-rhs[N+i])*(ax[i].y-rhs[N+i]);
        }
        ok = error < 1e-4 && residual <= (100+49+9+64)*1e-12;
    }
    if (!ok) fprintf(stderr,"  continuity UV solver selftest FAIL: long band, %d iterations, error %.9g, residual %.9g\n",iterations,error,residual);
    Arena_dispose(&arena); return !ok;
}

static int act_uv_metric_preserved_test(const AsmChart *c, const float *before)
{
    Arena_T arena = Arena_new();
    AsmMetricFace *ref = ARENA_ALLOC(arena,c->nf*sizeof *ref);
    AsmMetricStats initial, after; int ok = AsmMetric_prepare(c,ref) == 0;
    if (ok) {
        AsmMetric_measure(c,ref,before,&initial); AsmMetric_measure(c,ref,c->uv,&after);
        ok = AsmMetric_preserved(&initial,&after);
        if (!ok) fprintf(stderr,"  continuity original-chart metric FAIL: chart %d, within10 %.9g -> %.9g, within25 %.9g -> %.9g, invalid %zu\n",
            c->id,initial.within10/initial.area,after.within10/after.area,
            initial.within25/initial.area,after.within25/after.area,after.invalid);
    }
    Arena_dispose(&arena); return ok;
}

static int act_uv_selftest(void)
{
    int fails = 0;
    for (int first_step = 0; first_step < 2; first_step++) for (int mirror = 0; mirror < 2; mirror++) {
        AsmRun run = {0}; run.arena = Arena_new();
        float xyz[2][243], uv[2][162], original[2][162]; int32_t faces[384]; uint8_t boundary[81];
        memset(boundary,0,sizeof boundary);
        for (int j = 0; j < 9; j++) for (int u = 0; u < 9; u++) {
            int v = j*9+u; boundary[v] = j == 0 || j == 8 || u == 0 || u == 8;
            for (int side = 0; side < 2; side++) {
                xyz[side][v*3] = (float)(j*4); xyz[side][v*3+1] = (float)(u*5+side*44); xyz[side][v*3+2] = 0;
                uv[side][v*2] = (float)(u*5); uv[side][v*2+1] = (float)(j*4);
            }
            if (j < 8 && u < 8) {
                int32_t f[6] = {v,v+1,v+10,v,v+10,v+9}; memcpy(faces+(j*8+u)*6,f,sizeof f);
            }
        }
        for (int side = 0; side < 2; side++) {
            AsmChart c = {0}; c.id = side; c.component = 0; c.nv = 81; c.nf = 128;
            c.xyz = xyz[side]; c.uv = uv[side]; c.faces = faces; c.boundary = boundary; c.area3d = 1280;
            c.pose_theta = 0.4; c.pose_x = side*(mirror ? -44 : 44)*cos(0.4); c.pose_y = side*(mirror ? -44 : 44)*sin(0.4);
            if (mirror) c.flags = ASM_CHART_MIRROR;
            AsmRun_push_chart(&run,&c);
        }
        AsmRelation r = {0}; r.a = 0; r.b = 1; r.corr_count = 9; r.n_corr = 9; r.rms = 0.1;
        for (int j = 0; j < 9; j++) { AsmCorr c = {j*9+8,j*9}; AsmRun_push_corr(&run,&c); }
        AsmRun_push_rel(&run,&r); AsmContinuity_prepare(&run);
        AsmCorr gaps[9]; memcpy(gaps,run.corr,sizeof gaps);
        /* Add a nonrigid boundary bulge after measuring the physical trim.
         * A single pose cannot remove this error along the whole seam. */
        for (int j = 0; j < 9; j++) uv[1][j*18] += (float)(4*sin(j*3.141592653589793/8));
        memcpy(original,uv,sizeof uv);
        double before, after = 1e300; size_t support;
        if (!AsmContinuity_measure(&run,&run.rels[0],&before,&support) || support != 9) fails++;
        AsmChart chart_before[2] = {run.charts[0],run.charts[1]};
        run.charts[1].component = 1;
        if (act_relax_boundary_trial(&run,&run.charts[1],1,NULL,NULL,0.05,first_step) || memcmp(uv,original,sizeof uv)) fails++;
        run.charts[0].placement_state = ASM_PLACE_ROOT;
        if (act_relax_boundary_trial(&run,&run.charts[1],1,NULL,act_uv_reject_test,0.05,first_step) || memcmp(uv,original,sizeof uv)) fails++;
        if (!act_relax_boundary_trial(&run,&run.charts[1],1,NULL,NULL,0.05,first_step) ||
            !AsmContinuity_measure(&run,&run.rels[0],&after,&support) || after >= before*0.8 ||
            memcmp(uv[0],original[0],sizeof uv[0])) {
            fprintf(stderr,"  continuity fixed-boundary selftest: residual %.9g -> %.9g, mirror %d first step %d\n",before,after,mirror,first_step); fails++;
        }
        for (int j = 0; j < 9; j++) if (memcmp(uv[1]+j*18+16,original[1]+j*18+16,2*sizeof(float))) fails++;
        memcpy(uv,original,sizeof uv); run.charts[0] = chart_before[0]; run.charts[1] = chart_before[1];
        if (AsmContinuity_relax_uv(&run,act_uv_reject_test,0.05,NULL) || memcmp(uv,original,sizeof uv)) fails++;
        if (AsmContinuity_relax_uv(&run,NULL,0.05,NULL) != 1 ||
            !AsmContinuity_measure(&run,&run.rels[0],&after,&support) || after >= before*0.5) {
            fprintf(stderr,"  continuity UV selftest: residual %.9g -> %.9g, mirror %d\n",before,after,mirror); fails++;
        }
        if (memcmp(gaps,run.corr,sizeof gaps)) fails++;
        for (int side = 0; side < 2; side++) for (int j = 0; j < 9; j++) {
            int v = j*9+(side ? 8 : 0); /* Forty voxels from the seam: fixed. */
            if (memcmp(uv[side]+v*2,original[side]+v*2,2*sizeof(float))) fails++;
            for (int u = 0; u < 9; u++) {
                v = j*9+u;
                if (xyz[side][v*3] != j*4 || xyz[side][v*3+1] != u*5+side*44 || xyz[side][v*3+2] != 0) fails++;
            }
        }
        for (int j = 0; j < 9; j++) {
            double a[2], b[2]; act_point(&run.charts[0],j*9+8,NULL,a); act_point(&run.charts[1],j*9,NULL,b);
            /* The original physical gap is four voxels. Require the
             * two-voxel physical-witness accuracy as well as a nonzero gap;
             * an arbitrary 3.5-voxel lower bound is not the trim contract. */
            if (hypot(a[0]-b[0],a[1]-b[1]) < 2 || AsmContinuity_pair_error(&run,&run.rels[0],&run.corr[j]) > 2) fails++;
        }
        for (int side = 0; side < 2; side++) if (!act_uv_metric_preserved_test(&run.charts[side],original[side])) fails++;
        Arena_dispose(&run.arena);
    }
    return fails;
}

/* A source seam just outside the ordinary residual gate is still repairable
 * inside the existing displacement disks. It must reach the trial solver;
 * all original evidence, movement and acceptance checks remain required. */
static int act_uv_failed_source_selftest(void)
{
    int fails = 0;
    for (int mirror = 0; mirror < 2; mirror++) for (int placed_only = 0; placed_only < 2; placed_only++) {
        AsmRun run = {0}; run.arena = Arena_new();
        float xyz[2][243], uv[2][162], original[2][162];
        int32_t faces[384]; uint8_t boundary[81];
        for (int j = 0; j < 9; j++) for (int u = 0; u < 9; u++) {
            int v = j*9+u; boundary[v] = j == 0 || j == 8 || u == 0 || u == 8;
            for (int side = 0; side < 2; side++) {
                xyz[side][3*v] = (float)(4*j); xyz[side][3*v+1] = (float)(5*u+44*side); xyz[side][3*v+2] = 0;
                uv[side][2*v] = (float)(5*u); uv[side][2*v+1] = (float)(4*j);
            }
            if (j < 8 && u < 8) {
                int32_t f[6] = {v,v+1,v+10,v,v+10,v+9}; memcpy(faces+(j*8+u)*6,f,sizeof f);
            }
        }
        for (int side = 0; side < 2; side++) {
            AsmChart c = {0}; c.id = side; c.component = 0; c.nv = 81; c.nf = 128;
            c.xyz = xyz[side]; c.uv = uv[side]; c.faces = faces; c.boundary = boundary; c.area3d = 1280;
            c.flags = mirror ? ASM_CHART_MIRROR : 0;
            c.pose_theta = .4;
            c.pose_x = side*(mirror ? -54 : 54)*cos(.4); c.pose_y = side*(mirror ? -54 : 54)*sin(.4);
            AsmRun_push_chart(&run,&c);
        }
        AsmRelation r = {0}; r.a = 0; r.b = 1; r.flags = ASM_REL_WEAK; r.rms = .1;
        r.corr_count = r.n_corr = 9;
        for (int j = 0; j < 9; j++) { AsmCorr pair = {j*9+8,j*9}; AsmRun_push_corr(&run,&pair); }
        AsmRun_push_rel(&run,&r); AsmContinuity_prepare(&run);
        AsmRelation relation = run.rels[0]; AsmCorr witnesses[9]; memcpy(witnesses,run.corr,sizeof witnesses);
        memcpy(original,uv,sizeof original);
        double before, after = 1e300; size_t support;
        int ok = !AsmContinuity_measure(&run,&run.rels[0],&before,&support) && support == 9 && fabs(before-10) < 1e-6;
        if (placed_only) {
            for (int side = 0; side < 2; side++) { run.charts[side].placed = 1; run.charts[side].placement_state = side ? ASM_PLACE_NONE : ASM_PLACE_ROOT; }
            ok = ok && !AsmContinuity_relax_placed_uv(&run,act_registration_accept_test,.05,NULL) && !memcmp(original,uv,sizeof original);
            run.charts[1].placed = 0; run.charts[1].placement_state = ASM_PLACE_SEAM;
            ok = ok && !AsmContinuity_relax_placed_uv(&run,act_registration_accept_test,.05,NULL) && !memcmp(original,uv,sizeof original);
            run.charts[1].placed = 1;
            ok = ok && !AsmContinuity_relax_placed_uv(&run,NULL,.05,NULL) && !memcmp(original,uv,sizeof original);
        }
        size_t refused = placed_only ? AsmContinuity_relax_placed_uv(&run,act_uv_reject_test,.05,NULL)
                                     : AsmContinuity_relax_uv(&run,act_uv_reject_test,.05,NULL);
        ok = ok && !refused && !memcmp(original,uv,sizeof original);
        int accepted = (placed_only ? AsmContinuity_relax_placed_uv(&run,act_registration_accept_test,.05,NULL)
                                    : AsmContinuity_relax_uv(&run,NULL,.05,NULL)) == 1;
        /* The outer chart edges are fixed, so a ten-voxel pose error cannot
         * be halved freely by straining the 24-voxel boundary band. Require
         * a substantial repair inside the native source gate AND the
         * original per-chart metric budgets. Full two-voxel closure of this
         * constrained example still needs a larger movable neighbourhood. */
        ok = ok && accepted && AsmContinuity_measure(&run,&run.rels[0],&after,&support) && support == 9 && after < .6*before;
        for (int side = 0; side < 2; side++) ok = act_uv_metric_preserved_test(&run.charts[side],original[side]) && ok;
        ok = ok && !memcmp(&relation,&run.rels[0],sizeof relation) && !memcmp(witnesses,run.corr,sizeof witnesses);
        for (int side = 0; side < 2; side++) for (int v = 0; v < 81; v++) {
            ok = ok && hypot((double)uv[side][2*v]-original[side][2*v],(double)uv[side][2*v+1]-original[side][2*v+1]) <= 8;
            ok = ok && xyz[side][3*v] == (v/9)*4 && xyz[side][3*v+1] == (v%9)*5+side*44 && xyz[side][3*v+2] == 0;
        }
        for (int side = 0; side < 2; side++) for (int j = 0; j < 9; j++) {
            int v = j*9+(side ? 8 : 0);
            ok = ok && !memcmp(uv[side]+2*v,original[side]+2*v,2*sizeof(float));
        }
        if (!ok) { fprintf(stderr,"  continuity UV selftest FAIL: repairable failed source, mirror %d, placed-only %d, accepted %d, RMS %.9g -> %.9g\n",mirror,placed_only,accepted,before,after); fails++; }
        Arena_dispose(&run.arena);
    }
    return fails;
}

static int act_uv_local_bound_selftest(void)
{
    /* A long chain needs corrections over the local eight-voxel bound.
     * Its clipping must not scale the feasible correction of an independent
     * seam in the same trial. Compare that seam to its separate solve. */
    AsmRun run = {0}; run.arena = Arena_new();
    float xyz[8][54], uv[8][36], original[8][36], reference[2][36];
    uint8_t boundary[18]; int32_t faces[48], ids[8]; size_t rels[6];
    memset(boundary,1,sizeof boundary);
    for (int j = 0; j < 8; j++) { int q = 2*j; int32_t f[] = {q,q+1,q+3,q,q+3,q+2}; memcpy(faces+6*j,f,sizeof f); }
    for (int i = 0; i < 8; i++) {
        int k = i < 2 ? i : i-2; double base = i < 2 ? 0 : 1000;
        for (int j = 0; j < 9; j++) for (int u = 0; u < 2; u++) {
            int v = 2*j+u;
            xyz[i][3*v] = (float)(4*j); xyz[i][3*v+1] = (float)(base+14*k+10*u); xyz[i][3*v+2] = 0;
            uv[i][2*v] = (float)(10*u); uv[i][2*v+1] = (float)(4*j);
        }
        AsmChart c = {0}; c.id = i; c.component = 0; c.nv = 18; c.nf = 16; c.area3d = 320;
        c.xyz = xyz[i]; c.uv = uv[i]; c.faces = faces; c.boundary = boundary;
        c.pose_x = base+(i < 2 ? 18 : 20)*k; ids[i] = i;
        AsmRun_push_chart(&run,&c);
        if (i != 1 && i != 7) {
            AsmRelation r = {0}; r.a = i; r.b = i+1; r.rms = .1;
            r.corr_first = (int32_t)run.n_corr; r.corr_count = 9; r.n_corr = 9;
            for (int j = 0; j < 9; j++) { AsmCorr pair = {2*j+1,2*j}; AsmRun_push_corr(&run,&pair); }
            size_t ri = AsmRun_push_rel(&run,&r); rels[ri] = ri;
        }
    }
    AsmContinuity_prepare(&run); memcpy(original,uv,sizeof uv);
    Arena_T scratch = Arena_new();
    int ok = act_uv_component(&run,scratch,ids,2,rels,1,NULL,.05,NULL,NULL);
    memcpy(reference,uv,sizeof reference); memcpy(uv,original,sizeof uv); Arena_free(scratch);
    ok = ok && act_uv_component(&run,scratch,ids,8,rels,6,NULL,.05,NULL,NULL);
    double difference = 0, displacement = 0;
    for (int i = 0; i < 8; i++) for (int v = 0; v < 18; v++) {
        displacement = fmax(displacement,hypot(uv[i][2*v]-original[i][2*v],uv[i][2*v+1]-original[i][2*v+1]));
        if (i < 2) difference = fmax(difference,hypot(uv[i][2*v]-reference[i][2*v],uv[i][2*v+1]-reference[i][2*v+1]));
    }
    size_t pairs; double residual = sqrt(act_uv_energy(&run,rels,6,&pairs)/54);
    /* Uniform translations of the six-chain charts from +8 to -8 leave
     * only 2.8 vox per chain seam, and the independent seam can close.
     * Their combined RMS is 2.5561, with no chart strain or lost source.
     * Mere pointwise clipping leaves 3.510 and misses that feasible fit. */
    ok = ok && difference < 1e-3 && displacement <= 8 && displacement > 7.9 && residual < 2.7;
    for (int i = 0; i < 8; i++) ok = act_uv_metric_preserved_test(&run.charts[i],original[i]) && ok;
    if (!ok) fprintf(stderr,"  continuity UV selftest FAIL: independent bounded seam differs %.9g vox, maximum displacement %.9g, residual %.9g\n",difference,displacement,residual);
    /* Force overlapping two-chart groups along the same six-chart chain.
     * Every later group must retain each chart's original movement disk,
     * previously passing incident source seams, and immutable witnesses. */
    memcpy(uv,original,sizeof uv);
    AsmCorr witnesses[54]; AsmRelation relations[6];
    memcpy(witnesses,run.corr,sizeof witnesses); memcpy(relations,run.rels,sizeof relations);
    for (int i = 0; i < 8; i++) { run.charts[i].placed = 1; run.charts[i].placement_state = i ? ASM_PLACE_SEAM : ASM_PLACE_ROOT; }
    int groups_ok = !act_relax_placed_uv(&run,act_uv_reject_test,.05,NULL,36) && !memcmp(uv,original,sizeof uv);
    size_t accepted = act_relax_placed_uv(&run,act_registration_accept_test,.05,NULL,36);
    displacement = 0;
    for (int i = 0; i < 8; i++) for (int v = 0; v < 18; v++)
        displacement = fmax(displacement,hypot((double)uv[i][2*v]-original[i][2*v],(double)uv[i][2*v+1]-original[i][2*v+1]));
    groups_ok = groups_ok && accepted >= 2 && displacement <= 8 && displacement > 1 &&
                !memcmp(witnesses,run.corr,sizeof witnesses) && !memcmp(relations,run.rels,sizeof relations);
    for (size_t j = 0; j < run.n_rels; j++) {
        double rms; size_t support; groups_ok = groups_ok && AsmContinuity_measure(&run,&run.rels[j],&rms,&support);
    }
    if (!groups_ok) fprintf(stderr,"  continuity UV selftest FAIL: overlapping placed groups, %zu accepted, total correction %.9g\n",accepted,displacement);
    ok = ok && groups_ok;
    Arena_dispose(&scratch); Arena_dispose(&run.arena); return !ok;
}

static int act_uv_locked_solver_selftest(void)
{
    int fails = 0;
    for (int mirror = 0; mirror < 2; mirror++) {
        Arena_T arena = Arena_new();
        float xyz[] = {0,0,0, 1,0,0, 0,1,0, 10,0,0, 11,0,0, 10,1,0};
        float uv[] = {0,0, .4f,0, 0,1, 10,0, 11,0, 10,1}, trial[12] = {0};
        int32_t faces[] = {0,1,2, 3,4,5}, map[] = {-1,0,-1,1,2,3}; size_t off[] = {0,6};
        AsmChart c = {0}; c.nv = 6; c.nf = 2; c.xyz = xyz; c.uv = uv; c.faces = faces;
        c.pose_theta = .37; if (mirror) c.flags = ASM_CHART_MIRROR;
        double ct = (mirror ? -1 : 1)*cos(c.pose_theta), sn = (mirror ? -1 : 1)*sin(c.pose_theta);
        ActSpring spring[] = {{0,-1,1,1},{0,-1,1,1},{1,2,1,1},{2,3,1,1},{3,1,1,1}};
        ActMatch match[] = {act_match_pair(0,-1,-.6*ct,-.6*sn),act_match_pair(1,-1,ct,sn),act_match_pair(2,-1,ct,sn),act_match_pair(3,-1,ct,sn)};
        ActVec linear[4] = {0}, local[4] = {0}; int iterations = 0; size_t locked = 0;
        int ok = act_uv_solve(arena,spring,5,match,4,linear,4,&iterations);
        /* The first triangle already has compression 2.5. Its free vertex
         * would increase that compression, so even a short global step is
         * blocked. The other triangle's uniform unit translation is exact
         * up to the numerical tether and must remain independently possible. */
        ok = ok && act_uv_locked_direction(arena,&c,NULL,1,off,map,spring,5,match,4,linear,local,4,&locked,NULL,NULL,NULL,0,NULL,NULL);
        double error = hypot(local[0].x,local[0].y), area, old_area;
        for (int i = 1; i < 4; i++) error = fmax(error,hypot(local[i].x-ct/1.0001,local[i].y-sn/1.0001));
        ok = ok && isfinite(act_uv_trial_chart(&c,NULL,map,local,1.0,trial)) && locked == 1 && error < 1e-7;
        for (int f = 0; f < 2; f++) {
            double old_q = act_uv_quality(&c,faces+3*f,uv,&old_area), new_q = act_uv_quality(&c,faces+3*f,trial,&area);
            ok = ok && old_area*area > 0 && new_q <= old_q*(1+1e-6);
        }
        ok = ok && !memcmp(uv,trial,6*sizeof(float));
        if (!ok) { fprintf(stderr,"  continuity constrained solver selftest FAIL: mirror %d, %zu held, error %.9g\n",mirror,locked,error); fails++; }
        Arena_dispose(&arena);
    }
    return fails;
}

/* Geometry-feedback fixture: charts0/1 must stay fixed, but an independent
 * chart2 can translate. Allocate the first baseline just like the actual
 * audit so the test also checks that helper scratch never escapes. */
static int act_uv_feedback_test(AsmRun *run, const AsmChart *saved, size_t n, AsmContinuityAuditCache *cache)
{
    if (!cache->baseline_ready) {
        cache->baseline_violations = ARENA_CALLOC(cache->scratch,1,1);
        cache->baseline_violations[0] = 0x7a; cache->baseline_queries = 1; cache->baseline_ready = 1;
    }
    if (!cache->baseline_violations || cache->baseline_queries != 1 || cache->baseline_violations[0] != 0x7a ||
        !cache->blocked_charts || cache->blocked_charts_count != run->n_charts) return 0;
    memset(cache->blocked_charts,0,run->n_charts);
    for (size_t i = 0; i < n; i++) if (saved[i].id != 2 &&
        memcmp(saved[i].uv,run->charts[saved[i].id].uv,2*saved[i].nv*sizeof(float))) {
        cache->blocked_charts[0] = cache->blocked_charts[1] = 1; return 0;
    }
    return 1;
}

static int act_uv_feedback_solver_selftest(void)
{
    int fails = 0;
    for (int mirror = 0; mirror < 2; mirror++) for (int warm = 0; warm < 2; warm++) {
        AsmRun run = {0}; run.arena = Arena_new(); Arena_T scratch = Arena_new();
        float xyz[] = {0,0,0, 1,0,0, 0,1,0}, uv[] = {0,0, 1,0, 0,1};
        int32_t faces[] = {0,1,2}, map[] = {0,1,2,3,4,5,6,7,8}, order[] = {2,0,1};
        size_t off[] = {0,3,6,9};
        for (int i = 0; i < 3; i++) {
            AsmChart c = {0}; c.id = i; c.nv = 3; c.nf = 1; c.xyz = xyz; c.uv = uv; c.faces = faces;
            c.pose_x = 10*i; c.pose_theta = .37; if (mirror) c.flags = ASM_CHART_MIRROR;
            AsmRun_push_chart(&run,&c);
        }
        AsmChart saved[3], original[3]; memcpy(original,run.charts,sizeof original);
        for (int i = 0; i < 3; i++) saved[i] = run.charts[order[i]];
        ActSpring spring[9]; ActMatch match[9];
        double ct = (mirror ? -1 : 1)*cos(.37), sn = (mirror ? -1 : 1)*sin(.37);
        for (int i = 0; i < 9; i++) {
            spring[i] = (ActSpring){i,3*(i/3)+(i+1)%3,1,1};
            match[i] = act_match_pair(i,-1,ct,sn);
        }
        AsmContinuityAuditCache cache = {0}; cache.scratch = run.arena;
        if (warm) {
            cache.baseline_ready = 1; cache.baseline_queries = 1;
            cache.baseline_violations = ARENA_CALLOC(run.arena,1,1); cache.baseline_violations[0] = 0x7a;
        }
        AsmContinuityAuditCache original_cache = cache;
        ActVec linear[9] = {0}, local[9] = {0}; int iterations = 0; size_t locked = 0;
        uint8_t held[9]; ActVec solution[9], first[9];
        ActUvDirectionCache solve_cache = {held,solution,0,0};
        int ok = act_uv_solve(scratch,spring,9,match,9,linear,9,&iterations) &&
                 act_uv_locked_direction(scratch,saved,NULL,3,off,map,spring,9,match,9,linear,local,9,&locked,
                                         &run,act_uv_feedback_test,&cache,0,&solve_cache,NULL);
        memcpy(first,local,sizeof first);
        /* Both directions hold the same two charts for geometry here. The
         * cached answer still undergoes the callback and must leave its
         * independent third chart free, with identical rollback/lifetimes. */
        ok = ok && act_uv_locked_direction(scratch,saved,NULL,3,off,map,spring,9,match,9,linear,local,9,&locked,
                                           &run,act_uv_feedback_test,&cache,1,&solve_cache,NULL) &&
             solve_cache.hits == 1 && !memcmp(first,local,sizeof first);
        double error = 0;
        for (int i = 0; i < 9; i++) error = fmax(error,hypot(local[i].x-(i < 3 ? ct/1.0001 : 0),
                                                                         local[i].y-(i < 3 ? sn/1.0001 : 0)));
        ok = ok && locked == 6 && error < 1e-7 && !memcmp(original,run.charts,sizeof original) &&
             !memcmp(&cache,&original_cache,sizeof cache) && (!warm || cache.baseline_violations[0] == 0x7a);
        if (!ok) {
            fprintf(stderr,"  continuity geometry-constrained solver selftest FAIL: mirror %d warm %d, %zu held, error %.9g\n",
                    mirror,warm,locked,error); fails++;
        }
        Arena_dispose(&scratch); Arena_dispose(&run.arena);
    }
    return fails;
}

static int act_uv_held_chart_selftest(void)
{
    int fails = 0;
    for (int mirror = 0; mirror < 2; mirror++) {
        Arena_T arena = Arena_new();
        float xyz[] = {0,0,0, 1,0,0, 0,1,0, 10,0,0, 11,0,0, 10,1,0, 20,0,0, 21,0,0, 20,1,0};
        float uv[] = {0,0, .4f,0, 0,1, 10,0, 11,0, 10,1, 20,0, 21,0, 20,1}, trial[18] = {0};
        int32_t faces[] = {0,1,2, 3,4,5}, map[] = {-1,0,-1,1,2,3,4,5,6}; size_t off[] = {0,6,9};
        AsmChart saved[2] = {0};
        for (int i = 0; i < 2; i++) {
            saved[i].id = i; saved[i].nv = i ? 3 : 6; saved[i].nf = i ? 1 : 2;
            saved[i].xyz = xyz+3*off[i]; saved[i].uv = uv+2*off[i]; saved[i].faces = faces;
            saved[i].pose_theta = .37; if (mirror) saved[i].flags = ASM_CHART_MIRROR;
        }
        double ct = (mirror ? -1 : 1)*cos(.37), sn = (mirror ? -1 : 1)*sin(.37);
        ActSpring spring[] = {{0,-1,1,1},{0,-1,1,1},{1,2,1,1},{2,3,1,1},{3,1,1,1},
                              {4,5,1,1},{5,6,1,1},{6,4,1,1}};
        ActMatch match[7];
        for (int i = 0; i < 7; i++) match[i] = act_match_pair(i,-1,(i ? 1 : -.6)*ct,(i ? 1 : -.6)*sn);
        ActVec linear[7] = {0}, local[7] = {0}; int iterations = 0; size_t locked = 0;
        uint8_t held[7]; ActVec solution[7]; ActUvDirectionCache solve_cache = {held,solution,0,0};
        int ok = act_uv_solve(arena,spring,8,match,7,linear,7,&iterations) &&
                 act_uv_locked_direction(arena,saved,NULL,2,off,map,spring,8,match,7,linear,local,7,&locked,NULL,NULL,NULL,0,&solve_cache,NULL);
        ok = ok && locked == 1 && hypot(local[1].x,local[1].y) > .9;
        /* A different held set must not reuse the previous solution: the
         * whole-chart direction now also holds the other sound triangle. */
        ok = ok && act_uv_locked_direction(arena,saved,NULL,2,off,map,spring,8,match,7,linear,local,7,&locked,NULL,NULL,NULL,1,&solve_cache,NULL) &&
             solve_cache.hits == 0;
        /* Chart0 has both a strained and a sound face: the whole-chart
         * alternative must hold BOTH. Chart1 must still be free to close. */
        double error = 0;
        for (int i = 0; i < 7; i++) error = fmax(error,hypot(local[i].x-(i >= 4 ? ct/1.0001 : 0),
                                                                         local[i].y-(i >= 4 ? sn/1.0001 : 0)));
        for (int i = 0; i < 2; i++) ok = ok && isfinite(act_uv_trial_chart(&saved[i],NULL,map+off[i],local,1,trial+2*off[i]));
        ok = ok && locked == 4 && error < 1e-7 && !memcmp(uv,trial,12*sizeof(float));
        if (!ok) { fprintf(stderr,"  continuity held-chart selftest FAIL: mirror %d, %zu held, error %.9g\n",mirror,locked,error); fails++; }
        Arena_dispose(&arena);
    }
    return fails;
}

static int act_uv_bounded_solver_selftest(void)
{
    int fails = 0;
    for (int shifted = 0; shifted < 2; shifted++) for (int rotated = 0; rotated < 2; rotated++) {
        Arena_T arena = Arena_new();
        double angle = rotated ? .37 : 0, ct = cos(angle), sn = sin(angle), shift = shifted ? 2 : 0;
        ActMatch match[5]; ActVec linear[6] = {0}, offset[6], bounded[6] = {0}; int iterations = 0;
        for (int i = 0; i < 6; i++) offset[i] = (ActVec){shift*ct,shift*sn};
        for (int i = 0; i < 5; i++) match[i] = act_match_pair(i,i+1,6*ct,6*sn);
        int ok = act_uv_solve(arena,NULL,0,match,5,linear,6,&iterations) &&
            act_uv_bounded_direction(arena,NULL,0,match,5,linear,offset,bounded,6,&iterations);
        /* Endpoint constraints are active. Solve the four free unknowns
         * independently, then check both coordinates and the disk bounds. */
        int rows[7] = {0,1,2,3,1,2,3}, cols[7] = {0,1,2,3,0,1,2};
        double values[7] = {2.0001,2.0001,2.0001,2.0001,-1,-1,-1};
        double rhs[4] = {8-shift,0,0,-8-shift}, interior[4] = {0}, error = 0;
        ok = ok && Sparse_solve_sym(4,7,rows,cols,values,rhs,interior,SPARSE_SPD) == 0;
        for (int i = 0; i < 6; i++) {
            double expected = i == 0 ? 8-shift : i == 5 ? -8-shift : interior[i-1];
            error = fmax(error,hypot(bounded[i].x-expected*ct,bounded[i].y-expected*sn));
            ok = ok && hypot(bounded[i].x+offset[i].x,bounded[i].y+offset[i].y) <= 8+1e-10;
        }
        ok = ok && error < 1e-4;
        if (!ok) { fprintf(stderr,"  continuity bounded solver selftest FAIL: shifted %d rotated %d, error %.9g, iterations %d\n",shifted,rotated,error,iterations); fails++; }
        Arena_dispose(&arena);
    }
    return fails;
}

static int act_trim_transport_selftest(void)
{
    int failures = 0;
    for (int mirror = 0; mirror < 2; mirror++) for (int rotated = 0; rotated < 2; rotated++) {
        AsmRun run = {0}; run.arena = Arena_new();
        float xyz[2][54], uv[2][36]; double original[2][36];
        uint8_t boundary[18]; int32_t faces[48]; memset(boundary,1,sizeof boundary);
        for (int j = 0; j < 8; j++) {
            int q = 2*j, f[6] = {q,q+1,q+3,q,q+3,q+2}; memcpy(faces+6*j,f,sizeof f);
        }
        for (int side = 0; side < 2; side++) {
            for (int j = 0; j < 9; j++) for (int k = 0; k < 2; k++) {
                int v = 2*j+k; xyz[side][3*v] = (float)(14*side+10*k);
                xyz[side][3*v+1] = (float)(4*j); xyz[side][3*v+2] = 0;
                uv[side][2*v] = xyz[side][3*v]; uv[side][2*v+1] = xyz[side][3*v+1];
            }
            AsmChart c = {0}; c.id = side; c.nv = 18; c.nf = 16; c.area3d = 320;
            c.xyz = xyz[side]; c.uv = uv[side]; c.faces = faces; c.boundary = boundary;
            c.flags = mirror ? ASM_CHART_MIRROR : 0; c.pose_x = 100; c.pose_y = -31;
            AsmRun_push_chart(&run,&c);
            for (int v = 0; v < 18; v++) act_point(&run.charts[side],v,NULL,original[side]+2*v);
        }
        AsmRelation r = {0}; r.a = 0; r.b = 1; r.corr_count = r.n_corr = 9; r.rms = .1;
        for (int j = 0; j < 9; j++) { AsmCorr c = {2*j+1,2*j}; AsmRun_push_corr(&run,&c); }
        AsmRun_push_rel(&run,&r);
        int ok = !AsmContinuity_prepare(&run); double before=INFINITY, after=INFINITY, world_error = 0; size_t support;
        ok = ok && AsmContinuity_measure(&run,&run.rels[0],&before,&support) && before < 1e-5 && support == 9;
        /* Change only each chart's coordinate frame. Its assembled points,
         * faces, source data and physical seam must remain unchanged. */
        for (int side = 0; side < 2; side++) {
            double angle = rotated ? (side ? -.47 : .3) : (side ? 1.1 : -.8);
            double ct = cos(angle), sn = sin(angle);
            for (int v = 0; v < 18; v++) {
                double x = uv[side][2*v], y = uv[side][2*v+1];
                uv[side][2*v] = (float)(ct*x-sn*y); uv[side][2*v+1] = (float)(sn*x+ct*y);
            }
            run.charts[side].pose_theta = mirror ? angle : -angle;
            for (int v = 0; v < 18; v++) {
                double p[2]; act_point(&run.charts[side],v,NULL,p);
                world_error = fmax(world_error,hypot(p[0]-original[side][2*v],p[1]-original[side][2*v+1]));
            }
        }
        AsmContinuity_measure(&run,&run.rels[0],&after,&support);
        ok = ok && world_error < 1e-5 && after < 1e-5 && support == 9;
        /* A common affine deformation must also carry both physical trim
         * directions. The two charts retain different local frames. */
        for (int side = 0; side < 2; side++) {
            AsmChart *c = &run.charts[side]; double ct = cos(c->pose_theta), sn = sin(c->pose_theta);
            for (int v = 0; v < 18; v++) {
                double p[2]; act_point(c,v,NULL,p);
                double x = p[0]-c->pose_x, y = p[1]-c->pose_y;
                double a = 1.03*x+.04*y, b = -.02*x+.98*y;
                double u = ct*a+sn*b;
                c->uv[2*v] = (float)(mirror ? -u : u); c->uv[2*v+1] = (float)(-sn*a+ct*b);
            }
        }
        double affine;
        AsmContinuity_measure(&run,&run.rels[0],&affine,&support);
        ok = ok && affine < 1e-5 && support == 9;
        fprintf(stderr,"  trim transport control: mirror %d rotation %d, unchanged-world error %.9g, seam %.9g -> %.9g; %s\n",
                mirror,rotated,world_error,before,after,ok ? "PASS" : "FAIL");
        fprintf(stderr,"  trim affine control: mirror %d rotation %d, seam %.9g; %s\n",mirror,rotated,affine,ok ? "PASS" : "FAIL");
        failures += !ok; Arena_dispose(&run.arena);
    }
    return failures;
}

static int act_sparse_incumbent_selftest(void)
{
    int fail = 0;
    for (int veto = 0; veto < 4; veto++) {
        Arena_T arena = Arena_new(); AsmRun run = {0}; run.arena = arena;
        float xyz[3][18], uv[12];
        int32_t faces[] = {0,1,3,0,3,2,2,3,5,2,5,4};
        uint8_t boundary[] = {1,1,1,1,1,1};
        for (int row = 0; row < 3; row++) for (int col = 0; col < 2; col++) {
            int v = 2*row+col; uv[2*v] = 20.0f*col; uv[2*v+1] = 10.0f*row;
            for (int c = 0; c < 3; c++) {
                xyz[c][3*v] = 0; xyz[c][3*v+1] = uv[2*v]+22.0f*c; xyz[c][3*v+2] = uv[2*v+1];
            }
        }
        for (int i = 0; i < 3; i++) {
            AsmChart c = {0}; c.id = i; c.xyz = xyz[i]; c.uv = uv;
            c.faces = faces; c.boundary = boundary; c.nv = 6; c.nf = 4;
            c.area3d = 400; c.placed = 1; c.pose_x = 22.0*i;
            /* Chart 2 shares a label but has no seam: it must separate. */
            AsmRun_push_chart(&run,&c);
        }
        AsmRelation r = {0}; r.a = 0; r.b = 1; r.rms = 1;
        r.corr_count = 3; r.n_corr = 3; r.seam_len = 20;
        if (veto == 1) r.flags = ASM_REL_DROPPED;
        if (veto == 2) r.flags = ASM_REL_CONTACT;
        if (veto == 3) run.charts[1].pose_x += 100;
        AsmRun_push_rel(&run,&r);
        for (int i = 0; i < 3; i++) {
            AsmCorr c = {0}; c.va = 2*i+1; c.vb = 2*i; AsmRun_push_corr(&run,&c);
        }
        int ok = AsmContinuity_prepare(&run) == 0;
        double rms; size_t support;
        ok = ok && !AsmContinuity_measure(&run,&run.rels[0],&rms,&support);
        AsmContinuity_components(&run);
        ok = ok && ((run.charts[0].component == run.charts[1].component) == (veto == 0));
        ok = ok && run.charts[2].component != run.charts[0].component &&
             !(run.rels[0].continuity & ASM_CONT_VERIFIED);
        /* Rebuilding components again must not reinterpret the unresolved
         * certificate's SWITCHED flag as permission to fragment its bundle. */
        AsmContinuity_components(&run);
        ok = ok && ((run.charts[0].component == run.charts[1].component) == (veto == 0));
        if (!veto) {
            ok = ok && !AsmRel_is_join(&run.rels[0]) &&
                 (run.rels[0].continuity & ASM_CONT_INCUMBENT) && AsmContinuity_preserves_incumbent(&run);
            AsmChart saved[2] = {run.charts[0],run.charts[1]};
            ok = ok && AsmContinuity_preserves_incumbent_trial(&run,saved,2);
            for (int i = 0; i < 2; i++) {
                double x = run.charts[i].pose_x;
                run.charts[i].pose_theta = .37;
                run.charts[i].pose_x = cos(.37)*x+800;
                run.charts[i].pose_y = sin(.37)*x-120;
            }
            ok = ok && AsmContinuity_preserves_incumbent(&run);
            ok = ok && AsmContinuity_preserves_incumbent_trial(&run,saved,2);
            run.charts[1].pose_x += 20;
            ok = ok && !AsmContinuity_preserves_incumbent(&run);
            ok = ok && !AsmContinuity_preserves_incumbent_trial(&run,saved,2) &&
                 !AsmContinuity_preserves_incumbent_trial(&run,saved,1) &&
                 !AsmContinuity_preserves_incumbent_trial(&run,saved+1,1);
            for (int i = 0; i < 3; i++) ok = ok && run.corr[i].run == -1;
        } else ok = ok && !(run.rels[0].continuity & ASM_CONT_INCUMBENT);
        if (!ok) { fprintf(stderr,"  sparse incumbent selftest FAIL: veto %d\n",veto); fail++; }
        Arena_dispose(&arena);
    }
    return fail;
}

static int act_bounded_pose_selftest(void)
{
    /* A ten-percent spread mismatch makes the first GN rotation overshoot
     * the trust bound, although the best rigid fit lies strictly inside it.
     * The material itself must stay rigid; targets are observations only. */
    int failures=0;
    for (int mirror=0; mirror<2; mirror++) {
        AsmRun run={0}; run.arena=Arena_new();
        float uv[24]; AsmContinuityTarget targets[12];
        AsmChart chart={0}; chart.id=0; chart.nv=12; chart.uv=uv;
        chart.flags=mirror ? ASM_CHART_MIRROR : 0;
        double angle=0.24,scale=1.1;
        for (int j=0; j<12; j++) {
            double theta=2*3.14159265358979323846*j/12,u=40*cos(theta),v=40*sin(theta);
            uv[2*j]=(float)(mirror ? -u : u); uv[2*j+1]=(float)v;
            memset(&targets[j],0,sizeof targets[j]); targets[j].chart=0;
            targets[j].u=uv[2*j]; targets[j].v=uv[2*j+1]; targets[j].weight=1;
            double x=mirror ? -uv[2*j] : uv[2*j],y=uv[2*j+1];
            targets[j].x=scale*(cos(angle)*x-sin(angle)*y);
            targets[j].y=scale*(sin(angle)*x+cos(angle)*y);
        }
        AsmRun_push_chart(&run,&chart);
        int solved=AsmContinuity_close(&run,&chart,1,NULL,targets,12);
        AsmChart *result=&run.charts[0];
        double error=fabs(result->pose_theta-angle),shift=hypot(result->pose_x,result->pose_y);
        int ok=solved && error<1e-5 && shift<1e-5 && fabs(result->pose_theta)<=0.25 && result->uv==uv;
        if (!ok) {
            fprintf(stderr,"  continuity bounded pose control FAIL: mirror %d solved %d angle %.9g error %.9g shift %.9g\n",mirror,solved,result->pose_theta,error,shift);
            failures++;
        }
        Arena_dispose(&run.arena);
    }
    return failures;
}

int AsmContinuity_selftest(void)
{
    int fail = 0;
#define ACT_CHECK(test) do { int failures = test(); fail += failures; if (failures) fprintf(stderr,"  %s: %d failures\n",#test,failures); } while (0)
    ACT_CHECK(AsmMetric_selftest); ACT_CHECK(AsmBoundary_selftest);
    ACT_CHECK(act_trim_transport_selftest);
    ACT_CHECK(act_pose_plane_feasibility_selftest); ACT_CHECK(act_pose_plane_curvature_selftest);
    ACT_CHECK(act_halfspace_selftest); ACT_CHECK(act_registration_hypothesis_selftest);
    ACT_CHECK(act_uv_failed_source_selftest); ACT_CHECK(act_uv_selftest);
    ACT_CHECK(act_uv_solver_selftest); ACT_CHECK(act_uv_grid_solver_selftest);
    ACT_CHECK(act_trim_order_selftest); ACT_CHECK(act_uv_local_bound_selftest);
    ACT_CHECK(act_uv_bounded_solver_selftest); ACT_CHECK(act_uv_locked_solver_selftest);
    ACT_CHECK(act_uv_feedback_solver_selftest); ACT_CHECK(act_uv_held_chart_selftest);
    ACT_CHECK(act_registration_source_selftest);
    ACT_CHECK(act_sparse_incumbent_selftest);
    ACT_CHECK(act_bounded_pose_selftest);
#undef ACT_CHECK
    {
        /* A previously used displacement disk permits16 vox of movement
         * from the current point. Its newly touched neighbour at40 must be
         * included; a chart beyond that entire rim neighbourhood is fixed. */
        double moving[4] = {0,0,0,0}, positive[4] = {40,0,40,0}, negative[4] = {-40,0,-40,0};
        double outside[4] = {40.01,0,41,0};
        if (!act_uv_bounds_near(moving,positive) || !act_uv_bounds_near(moving,negative) ||
            act_uv_bounds_near(moving,outside)) fail++;
    }
    {
        /* Empty or absent correspondence spans are valid ledger entries.
         * Invalid offsets/counts must never become diagnostic reads. */
        AsmCorr corr[1] = {{0}}; corr[0].valid = 3;
        AsmRelation rels[5] = {{0}};
        for (int i = 0; i < 5; i++) rels[i].corr_count = 1;
        rels[1].corr_first = -1; rels[2].corr_first = INT32_MAX;
        rels[3].corr_count = 2; rels[4].corr_count = -1;
        AsmRun run = {0}; run.rels = rels; run.n_rels = 5; run.corr = corr; run.n_corr = 1;
        if (act_write_witnesses(&run,NULL) != 1) fail++;
        run.corr = NULL; run.n_corr = 0;
        if (act_write_witnesses(&run,NULL)) fail++;
    }
    AsmChart c; memset(&c, 0, sizeof c);
    float xyz[] = { 0,0,0, 2,0,0, 0,2,0 };
    float uv[] = { 0,0, 0,2, -2,0 };
    int32_t f[] = { 0,1,2 }; c.xyz = xyz; c.uv = uv; c.faces = f; c.nf = 1; c.nv = 3;
    double d[] = { 3,4,5 }; float gap[2];
    if (!act_gap(&c, 0, d, gap) || fabs(gap[0]+4) > 1e-6 || fabs(gap[1]-3) > 1e-6) fail++;
    c.flags = ASM_CHART_MIRROR; c.pose_theta = 0.5; c.pose_x = 31; c.pose_y = -8;
    double p[2], q[2]; act_point(&c, 0, NULL, p); act_point(&c, 0, gap, q);
    if (fabs(hypot(q[0]-p[0],q[1]-p[1])-5) > 1e-6) fail++;
    f[2] = 1; if (act_gap(&c, 0, d, gap)) fail++;
    {
        /* Two strips separated by a physical four-voxel trim. Recover a
         * translated/rotated B against fixed A; then independently remove
         * the emitted bridge and ensure metadata cannot hide the break. */
        AsmRun run; memset(&run, 0, sizeof run); run.arena = Arena_new();
        float xyz2[2][54], uv2[2][36]; uint8_t boundary[18]; int32_t fs[48];
        memset(boundary, 1, sizeof boundary);
        for (int j = 0; j < 8; j++) {
            int q = j*2; int t[6] = {q,q+1,q+3,q,q+3,q+2}; memcpy(fs+6*j, t, sizeof t);
        }
        for (int side = 0; side < 2; side++) {
            AsmChart ch; memset(&ch, 0, sizeof ch);
            for (int j = 0; j < 9; j++) for (int u = 0; u < 2; u++) {
                int v = 2*j+u;
                xyz2[side][3*v] = (float)(j*4); xyz2[side][3*v+1] = (float)(10*u+(side ? 4 : -10)); xyz2[side][3*v+2] = 0;
                uv2[side][2*v] = (float)(10*u); uv2[side][2*v+1] = (float)(j*4);
            }
            ch.id = side; ch.component = 0; ch.placed = 1; ch.placement_state = side ? ASM_PLACE_SEAM : ASM_PLACE_ROOT;
            ch.nv = 18; ch.nf = 16; ch.xyz = xyz2[side]; ch.uv = uv2[side]; ch.faces = fs; ch.boundary = boundary; ch.area3d = 320;
            ch.pose_x = side ? 14 : 0; AsmRun_push_chart(&run, &ch);
        }
        AsmRelation r; memset(&r, 0, sizeof r); r.a = 0; r.b = 1; r.flags = ASM_REL_WEAK;
        r.corr_first = 0; r.corr_count = 9; r.n_corr = 9; r.rms = 0.1;
        for (int j = 0; j < 9; j++) { AsmCorr cr = {2*j+1,2*j}; AsmRun_push_corr(&run, &cr); }
        AsmRun_push_rel(&run, &r);
        if (AsmContinuity_prepare(&run)) fail++;
        double rms; size_t support;
        if (!AsmContinuity_measure(&run, &run.rels[0], &rms, &support) || rms > 1e-6 || support != 9) fail++;
        AsmChart saved = run.charts[1], fixed = run.charts[0];
        run.charts[1].pose_x += 12; run.charts[1].pose_y -= 7; run.charts[1].pose_theta += 0.08;
        if (!AsmContinuity_close(&run, &saved, 1, NULL, NULL, 0) || !AsmContinuity_measure(&run, &run.rels[0], &rms, &support) || rms > 0.01) fail++;
        if (memcmp(&fixed, &run.charts[0], sizeof fixed) != 0) fail++;
        /* A local pre-placement registration fixes a measured boundary
         * chart without changing its placement state. An explicit held
         * mask is required; a NONE chart must not silently become an anchor. */
        run.charts[0].placement_state = ASM_PLACE_NONE;
        AsmChart held_boundary = run.charts[0];
        run.charts[1] = saved; run.charts[1].pose_x += 12; run.charts[1].pose_y -= 7; run.charts[1].pose_theta += .08;
        AsmChart before_held = run.charts[1]; uint8_t held[] = {1,0}, selected[] = {1};
        if (act_close(&run,&saved,1,selected,NULL,0,0,1,NULL,0,NULL) ||
            memcmp(&before_held,&run.charts[1],sizeof before_held)) fail++;
        if (!act_close(&run,&saved,1,selected,NULL,0,0,1,NULL,0,held) ||
            !AsmContinuity_measure(&run,&run.rels[0],&rms,&support) || rms > .01 ||
            memcmp(&held_boundary,&run.charts[0],sizeof held_boundary)) fail++;
        run.charts[0] = fixed; run.charts[1] = saved;
        /* Different measured tangent directions need a symmetric fit.
         * The midpoint has six-voxel residual in either direction; a
         * one-sided solve moves it to twelve in the other direction. */
        for (int j = 0; j < 9; j++) run.corr[j].gap_b[0] = 16;
        run.charts[1].pose_x = 20; run.charts[1].pose_y = run.charts[1].pose_theta = 0;
        if (!AsmContinuity_close(&run, &saved, 1, NULL, NULL, 0) ||
            !AsmContinuity_measure(&run, &run.rels[0], &rms, &support) || fabs(rms-6) > 1e-5 ||
            fabs(run.charts[1].pose_x-20) > 1e-5) fail++;
        run.rels[0].a = 1; run.rels[0].b = 0;
        for (int j = 0; j < 9; j++) {
            AsmCorr *cr = &run.corr[j]; int32_t v = cr->va; cr->va = cr->vb; cr->vb = v;
            AsmTrim trim = cr->trim_a; cr->trim_a = cr->trim_b; cr->trim_b = trim;
            for (int k = 0; k < 2; k++) { float g = cr->gap_a[k]; cr->gap_a[k] = -cr->gap_b[k]; cr->gap_b[k] = -g; }
        }
        if (!AsmContinuity_close(&run, &saved, 1, NULL, NULL, 0) ||
            !AsmContinuity_measure(&run, &run.rels[0], &rms, &support) || fabs(rms-6) > 1e-5 ||
            fabs(run.charts[1].pose_x-20) > 1e-5) fail++;
        run.rels[0].a = 0; run.rels[0].b = 1;
        for (int j = 0; j < 9; j++) {
            AsmCorr *cr = &run.corr[j]; int32_t v = cr->va; cr->va = cr->vb; cr->vb = v;
            AsmTrim trim = cr->trim_a; cr->trim_a = cr->trim_b; cr->trim_b = trim;
            cr->gap_a[0] = cr->gap_b[0] = 4; cr->gap_a[1] = cr->gap_b[1] = 0;
        }
        run.charts[1] = saved;
        /* Seven central points and two displaced ends have a feasible
         * four-voxel correction. Huber fitting instead sacrifices both
         * ends: its RMS passes, but fewer than 90% of points do. */
        /* Translate each complete end edge. Moving only its witness would
         * also stretch the physical trim extension, changing this pose-only
         * robust-fit control into a different seam feasibility problem. */
        uv2[1][0] = uv2[1][32] = -11.5f;
        uv2[1][2] = uv2[1][34] = -1.5f;
        for (int mirror = 0; mirror < 2; mirror++) {
            run.charts[0] = fixed; run.charts[1] = saved;
            if (mirror) {
                run.charts[0].flags |= ASM_CHART_MIRROR;
                run.charts[1].flags |= ASM_CHART_MIRROR;
                run.charts[1].pose_x = -14;
            }
            AsmChart initial = run.charts[1], fixed_boundary = run.charts[0];
            run.charts[1].pose_x += mirror ? -4 : 4;
            if (!AsmContinuity_measure(&run,&run.rels[0],&rms,&support)) fail++;
            run.charts[1] = initial;
            if (!AsmContinuity_close(&run,&initial,1,NULL,NULL,0) ||
                AsmContinuity_measure(&run,&run.rels[0],&rms,&support) || rms >= 8) fail++;
            run.charts[1] = initial;
            if (!AsmContinuity_close_boundaries(&run,&initial,1,NULL) ||
                !AsmContinuity_measure(&run,&run.rels[0],&rms,&support) || support != 9 ||
                memcmp(&run.charts[0],&fixed_boundary,sizeof fixed_boundary)) fail++;
        }
        uv2[1][0] = uv2[1][32] = 0;
        uv2[1][2] = uv2[1][34] = 10;
        run.charts[0] = fixed; run.charts[1] = saved;
        if (AsmContinuity_confirm(&run) != 1 || !AsmRel_is_join(&run.rels[0])) fail++;
        float emitted_xyz[108]; memcpy(emitted_xyz, xyz2, sizeof emitted_xyz);
        int32_t chart[36]; for (int j = 0; j < 36; j++) chart[j] = j/18;
        int32_t bridge[] = {1,3,18, 3,20,18};
        AsmContinuityStats st;
        AsmContinuity_audit(&run, 0, chart, emitted_xyz, 36, bridge, 2, NULL, &st);
        if (!st.confetti_free || st.verified_joins != 1) fail++;
        int32_t islands[] = {0,1,3, 18,19,21};
        AsmContinuity_audit(&run, 0, chart, emitted_xyz, 36, islands, 2, NULL, &st);
        if (st.confetti_free || st.metadata_only_joins != 1 || st.unexplained_breaks != 1) fail++;
        run.charts[1].placement_state = ASM_PLACE_NONE;
        AsmContinuity_audit(&run, 0, chart, emitted_xyz, 36, islands, 2, NULL, &st);
        if (st.unsupported_islands != 1 || st.unsupported_area != 320) fail++;
        run.charts[1].placement_state = ASM_PLACE_SEAM;
        run.charts[1].pose_x += 1000;
        if (AsmContinuity_measure(&run, &run.rels[0], &rms, &support)) fail++;
        run.charts[1] = saved;
        for (int k = 0; k < 2; k++) {
            run.rels[0].flags = k ? ASM_REL_DROPPED : ASM_REL_CONTACT;
            if (AsmContinuity_measure(&run, &run.rels[0], &rms, &support) || AsmContinuity_confirm(&run)) fail++;
        }
        run.rels[0].flags = ASM_REL_PARITY;
        if (AsmContinuity_measure(&run, &run.rels[0], &rms, &support)) fail++;
        Arena_dispose(&run.arena);
    }
    fprintf(stderr, "  asm_continuity selftest: %d failure(s)\n", fail);
    return fail;
}
