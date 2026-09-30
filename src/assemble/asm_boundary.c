#include "asm_boundary.h"
#include "../common/csr.h"
#include "../common/pipeline_constants.h"
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct AbLinks { int32_t neighbour[2]; double length[2]; } AbLinks;
typedef struct AbStart { int32_t index, va, vb; int priority; } AbStart;

static double ab_distance(const AsmChart *c, int32_t a, int32_t b)
{
    double squared = 0;
    for (size_t k = 0; k < 3; k++) {
        double delta = (double)c->xyz[3*(size_t)a+k] - (double)c->xyz[3*(size_t)b+k];
        squared += delta*delta;
    }
    return sqrt(squared);
}

static int32_t ab_slot(const int32_t *off, const int32_t *col, int32_t a, int32_t b)
{
    int32_t lo = off[a], hi = off[(size_t)a+1];
    while (lo < hi) {
        int32_t mid = lo + (hi-lo)/2;
        if (col[mid] < b) lo = mid+1; else hi = mid;
    }
    return lo < off[(size_t)a+1] && col[lo] == b ? lo : -1;
}

static void ab_add_link(AbLinks *links, int32_t a, int32_t b, double length)
{
    if (a == b || length > ASM_CONT_RUN_MAX_STEP) return;
    int ka = links[a].neighbour[0] < 0 ? 0 : 1;
    int kb = links[b].neighbour[0] < 0 ? 0 : 1;
    assert(links[a].neighbour[ka] < 0 && links[b].neighbour[kb] < 0);
    links[a].neighbour[ka] = b; links[a].length[ka] = length;
    links[b].neighbour[kb] = a; links[b].length[kb] = length;
}

/* Only the boundary CSR slots participate. Walk through unsampled original
 * vertices, accumulating physical length; never bridge two separate loops. */
static int ab_chart_links(Arena_T arena, const AsmChart *c, const int32_t *sample_at, AbLinks *links, double *point_mass)
{
    Arena_Mark mark = Arena_save(arena);
    int status = -1;
    CSR_T graph = CSR_from_faces(arena,c->faces,c->nf,c->nv);
    const int32_t *off = CSR_offset(graph), *col = CSR_target(graph);
    size_t nnz = (size_t)CSR_nnz(graph);
    uint8_t *incidence = ARENA_CALLOC(arena,nnz,1);
    int32_t *boundary = ARENA_ALLOC(arena,2*c->nv*sizeof *boundary);
    uint8_t *degree = ARENA_CALLOC(arena,c->nv,1), *seen = ARENA_CALLOC(arena,c->nv,1);
    for (size_t i = 0; i < 2*c->nv; i++) boundary[i] = -1;
    for (size_t face = 0; face < c->nf; face++) {
        const int32_t *f = c->faces+3*face;
        /* Repeated-index zero faces supply neither a frame nor an edge. */
        if (f[0] == f[1] || f[0] == f[2] || f[1] == f[2]) continue;
        for (size_t k = 0; k < 3; k++) {
            int32_t a = f[k], b = f[(k+1)%3];
            int32_t ka = ab_slot(off,col,a,b), kb = ab_slot(off,col,b,a);
            if (ka < 0 || kb < 0 || incidence[ka] == 2 || incidence[kb] == 2) goto done;
            incidence[ka]++; incidence[kb]++;
        }
    }
    for (size_t v = 0; v < c->nv; v++) {
        for (int32_t k = off[v]; k < off[v+1]; k++) if (incidence[k] == 1) {
            if (degree[v] == 2) goto done;
            boundary[2*v+degree[v]] = col[k]; degree[v]++;
        }
        if (degree[v] == 1 || (sample_at[v] >= 0 && degree[v] != 2)) goto done;
        if (point_mass && sample_at[v] >= 0)
            point_mass[sample_at[v]] = .5*(ab_distance(c,(int32_t)v,boundary[2*v])+
                                           ab_distance(c,(int32_t)v,boundary[2*v+1]));
    }
    for (size_t start = 0; start < c->nv; start++) {
        if (!degree[start] || seen[start]) continue;
        int32_t current = (int32_t)start, previous = -1, first_sample = -1, last_sample = -1;
        double length = 0, first_position = 0, last_position = 0;
        size_t observations = 0;
        while (!seen[current]) {
            seen[current] = 1;
            int32_t sample = sample_at[current];
            if (sample >= 0) {
                if (first_sample < 0) { first_sample = sample; first_position = length; }
                if (last_sample >= 0) ab_add_link(links,last_sample,sample,length-last_position);
                last_sample = sample; last_position = length; observations++;
            }
            int32_t next = boundary[2*(size_t)current];
            if (next == previous) next = boundary[2*(size_t)current+1];
            double step = ab_distance(c,current,next);
            if (!(step > 0) || !isfinite(step)) goto done;
            length += step; previous = current; current = next;
        }
        if (current != (int32_t)start) goto done;
        if (observations >= 3) ab_add_link(links,last_sample,first_sample,length-last_position+first_position);
        else if (observations == 2) {
            /* Two observations do not identify which arc of a cycle they
             * represent. Keep the point obligations without inventing a run. */
            links[first_sample].neighbour[0] = links[first_sample].neighbour[1] = -1;
            links[last_sample].neighbour[0] = links[last_sample].neighbour[1] = -1;
        }
    }
    status = 0;
done:
    Arena_restore(arena,mark); return status;
}

static int ab_start_compare(const void *one, const void *two)
{
    const AbStart *a = one, *b = two;
    if (a->priority != b->priority) return a->priority < b->priority ? -1 : 1;
    if (a->va != b->va) return a->va < b->va ? -1 : 1;
    if (a->vb != b->vb) return a->vb < b->vb ? -1 : 1;
    return (a->index > b->index) - (a->index < b->index);
}

int AsmBoundary_order(Arena_T arena, const AsmChart *a, const AsmChart *b,
                      AsmCorr *samples, size_t count, AsmBoundaryStats *stats)
{
    assert(arena && a && b && stats && (samples || !count));
    memset(stats,0,sizeof *stats);
    if (!count) return 0;
    if (count > INT32_MAX || a->nv > INT32_MAX || b->nv > INT32_MAX ||
        a->nf > INT32_MAX/6 || b->nf > INT32_MAX/6) return -1;
    const AsmChart *charts[2] = {a,b};
    for (size_t side = 0; side < 2; side++) {
        const AsmChart *c = charts[side];
        if (!c->nv || !c->nf || !c->xyz || !c->faces) return -1;
        for (size_t k = 0; k < 3*c->nf; k++)
            if (c->faces[k] < 0 || (size_t)c->faces[k] >= c->nv) return -1;
    }
    Arena_Mark mark = Arena_save(arena);
    int status = -1;
    int32_t *at_a = ARENA_ALLOC(arena,a->nv*sizeof *at_a), *at_b = ARENA_ALLOC(arena,b->nv*sizeof *at_b);
    uint8_t *active = ARENA_CALLOC(arena,count,1), *seen = ARENA_CALLOC(arena,count,1);
    AbLinks *left = ARENA_CALLOC(arena,count,sizeof *left), *right = ARENA_CALLOC(arena,count,sizeof *right);
    AbLinks *joint = ARENA_CALLOC(arena,count,sizeof *joint);
    AbStart *starts = ARENA_ALLOC(arena,count*sizeof *starts);
    AsmCorr *ordered = ARENA_ALLOC(arena,count*sizeof *ordered);
    for (size_t i = 0; i < a->nv; i++) at_a[i] = -1;
    for (size_t i = 0; i < b->nv; i++) at_b[i] = -1;
    for (size_t i = 0; i < count; i++) {
        left[i].neighbour[0] = left[i].neighbour[1] = -1;
        right[i].neighbour[0] = right[i].neighbour[1] = -1;
        joint[i].neighbour[0] = joint[i].neighbour[1] = -1;
        if (samples[i].valid != 3) continue;
        int32_t va = samples[i].va, vb = samples[i].vb;
        if (va < 0 || vb < 0 || (size_t)va >= a->nv || (size_t)vb >= b->nv) goto done;
        int32_t ia = at_a[va], ib = at_b[vb];
        if (ia >= 0 || ib >= 0) {
            if (ia != ib) goto done;
            stats->duplicate_samples++; continue;
        }
        at_a[va] = at_b[vb] = (int32_t)i; active[i] = 1;
    }
    if (ab_chart_links(arena,a,at_a,left,NULL) || ab_chart_links(arena,b,at_b,right,NULL)) goto done;
    for (size_t i = 0; i < count; i++) if (active[i]) {
        for (size_t ka = 0; ka < 2; ka++) {
            int32_t j = left[i].neighbour[ka];
            if (j < 0 || (size_t)j <= i) continue;
            for (size_t kb = 0; kb < 2; kb++) if (right[i].neighbour[kb] == j &&
                fabs(left[i].length[ka]-right[i].length[kb]) <= ASM_CONT_RUN_MAX_LENGTH_ERROR) {
                ab_add_link(joint,(int32_t)i,j,.5*(left[i].length[ka]+right[i].length[kb]));
            }
        }
    }
    size_t nstart = 0, written = 0;
    for (size_t i = 0; i < count; i++) if (active[i]) {
        int degree = (joint[i].neighbour[0] >= 0) + (joint[i].neighbour[1] >= 0);
        starts[nstart++] = (AbStart){(int32_t)i,samples[i].va,samples[i].vb,degree == 1 ? 0 : degree == 2 ? 1 : 2};
    }
    qsort(starts,nstart,sizeof *starts,ab_start_compare);
    for (size_t seed = 0; seed < nstart; seed++) {
        int32_t current = starts[seed].index, previous = -1;
        if (seen[current]) continue;
        size_t begin = written;
        double length = 0;
        while (!seen[current]) {
            seen[current] = 1;
            ordered[written] = samples[current]; ordered[written].run = (int32_t)stats->runs; written++;
            int next_slot = -1;
            for (int k = 0; k < 2; k++) {
                int32_t neighbour = joint[current].neighbour[k];
                if (neighbour < 0 || neighbour == previous) continue;
                if (next_slot < 0 || samples[neighbour].va < samples[joint[current].neighbour[next_slot]].va) next_slot = k;
            }
            if (next_slot < 0) break;
            length += joint[current].length[next_slot];
            previous = current; current = joint[current].neighbour[next_slot];
        }
        if (written-begin >= ASM_CONT_RUN_MIN_SAMPLES && length >= ASM_CONT_RUN_MIN_LENGTH) {
            stats->supported_runs++; stats->support_length += length;
        } else stats->short_run_samples += written-begin;
        stats->runs++;
    }
    stats->ordered_samples = written;
    for (size_t i = 0; i < count; i++) if (!active[i]) {
        ordered[written] = samples[i]; ordered[written].run = -1; written++;
    }
    if (written != count) goto done;
    if (!stats->supported_runs) for (size_t i = 0; i < count; i++) ordered[i].run = -1;
    memcpy(samples,ordered,count*sizeof *samples); status = 0;
done:
    Arena_restore(arena,mark); return status;
}

/* The reference: one CSR build and full boundary walk per seam (selftest). */
static int ab_weights_direct(Arena_T arena, const AsmChart *a, const AsmChart *b,
                             const AsmCorr *samples, size_t count, double *weights)
{
    if (!arena || !a || !b || !samples || !weights || !count || count > INT32_MAX ||
        a->nv > INT32_MAX || b->nv > INT32_MAX || a->nf > INT32_MAX/6 || b->nf > INT32_MAX/6) return -1;
    const AsmChart *charts[2] = {a,b};
    for (int side = 0; side < 2; side++) {
        const AsmChart *c = charts[side];
        if (!c->nv || !c->nf || !c->xyz || !c->faces) return -1;
        for (size_t k = 0; k < 3*c->nf; k++) if (c->faces[k] < 0 || (size_t)c->faces[k] >= c->nv) return -1;
    }
    memset(weights,0,count*sizeof *weights);
    Arena_Mark mark = Arena_save(arena); int status = -1;
    int32_t *at[2]; AbLinks *links[2]; double *mass[2];
    for (int side = 0; side < 2; side++) {
        at[side] = ARENA_ALLOC(arena,charts[side]->nv*sizeof(int32_t));
        links[side] = ARENA_CALLOC(arena,count,sizeof(AbLinks));
        mass[side] = ARENA_CALLOC(arena,count,sizeof(double));
        for (size_t v = 0; v < charts[side]->nv; v++) at[side][v] = -1;
        for (size_t i = 0; i < count; i++) {
            links[side][i].neighbour[0] = links[side][i].neighbour[1] = -1;
            int32_t v = side ? samples[i].vb : samples[i].va;
            if (samples[i].valid != 3 || v < 0 || (size_t)v >= charts[side]->nv || at[side][v] >= 0) goto done;
            at[side][v] = (int32_t)i;
        }
        if (ab_chart_links(arena,charts[side],at[side],links[side],mass[side])) goto done;
    }
    for (size_t i = 0; i < count; i++) {
        for (int ka = 0; ka < 2; ka++) {
            int32_t j = links[0][i].neighbour[ka];
            if (j < 0 || (size_t)j <= i) continue;
            for (int kb = 0; kb < 2; kb++) if (links[1][i].neighbour[kb] == j &&
                fabs(links[0][i].length[ka]-links[1][i].length[kb]) <= ASM_CONT_RUN_MAX_LENGTH_ERROR) {
                double half = .25*(links[0][i].length[ka]+links[1][i].length[kb]);
                weights[i] += half; weights[j] += half;
            }
        }
    }
    for (size_t i = 0; i < count; i++) {
        if (weights[i] == 0) weights[i] = .5*(mass[0][i]+mass[1][i]);
        if (!(weights[i] > 0) || !isfinite(weights[i])) goto done;
    }
    status = 0;
done:
    Arena_restore(arena,mark); return status;
}

/* ---- per-chart loops, shared by every seam of the chart ------------------- */

struct AsmBoundaryLoops {
    int valid;          /* faces in range, manifold edges, every boundary loop closes with positive steps */
    int32_t *boundary;  /* [2*nv] the two boundary neighbours of a boundary vertex, else -1 */
    uint8_t *degree;    /* [nv] 0 or 2 */
    int32_t *start;     /* [nv] lowest vertex of the vertex's boundary loop, else -1 */
};

/* The seam-independent part of ab_chart_links, run once per chart: its checks
 * fail a chart for every seam alike, so they move here unchanged. */
AsmBoundaryLoops *AsmBoundary_loops(Arena_T arena, const AsmChart *c)
{
    AsmBoundaryLoops *L = ARENA_ALLOC(arena,sizeof *L); memset(L,0,sizeof *L);
    if (!c || !c->nv || !c->nf || !c->xyz || !c->faces || c->nv > INT32_MAX || c->nf > INT32_MAX/6) return L;
    for (size_t k = 0; k < 3*c->nf; k++) if (c->faces[k] < 0 || (size_t)c->faces[k] >= c->nv) return L;
    L->boundary = ARENA_ALLOC(arena,2*c->nv*sizeof *L->boundary);
    L->degree = ARENA_CALLOC(arena,c->nv,1);
    L->start = ARENA_ALLOC(arena,c->nv*sizeof *L->start);
    for (size_t i = 0; i < 2*c->nv; i++) L->boundary[i] = -1;
    for (size_t v = 0; v < c->nv; v++) L->start[v] = -1;
    Arena_Mark mark = Arena_save(arena);
    CSR_T graph = CSR_from_faces(arena,c->faces,c->nf,c->nv);
    const int32_t *off = CSR_offset(graph), *col = CSR_target(graph);
    size_t nnz = (size_t)CSR_nnz(graph);
    uint8_t *incidence = ARENA_CALLOC(arena,nnz?nnz:1,1);
    for (size_t face = 0; face < c->nf; face++) {
        const int32_t *f = c->faces+3*face;
        if (f[0] == f[1] || f[0] == f[2] || f[1] == f[2]) continue;
        for (size_t k = 0; k < 3; k++) {
            int32_t a = f[k], b = f[(k+1)%3];
            int32_t ka = ab_slot(off,col,a,b), kb = ab_slot(off,col,b,a);
            if (ka < 0 || kb < 0 || incidence[ka] == 2 || incidence[kb] == 2) goto done;
            incidence[ka]++; incidence[kb]++;
        }
    }
    for (size_t v = 0; v < c->nv; v++) {
        for (int32_t k = off[v]; k < off[v+1]; k++) if (incidence[k] == 1) {
            if (L->degree[v] == 2) goto done;
            L->boundary[2*v+L->degree[v]] = col[k]; L->degree[v]++;
        }
        if (L->degree[v] == 1) goto done;
    }
    for (size_t start = 0; start < c->nv; start++) {
        if (!L->degree[start] || L->start[start] >= 0) continue;
        int32_t current = (int32_t)start, previous = -1;
        while (L->start[current] < 0) {
            L->start[current] = (int32_t)start;
            int32_t next = L->boundary[2*(size_t)current];
            if (next == previous) next = L->boundary[2*(size_t)current+1];
            double step = ab_distance(c,current,next);
            if (!(step > 0) || !isfinite(step)) goto done;
            previous = current; current = next;
        }
        if (current != (int32_t)start) goto done;
    }
    L->valid = 1;
done:
    Arena_restore(arena,mark); return L;
}

static int ab_int32_order(const void *one, const void *two)
{
    int32_t a = *(const int32_t *)one, b = *(const int32_t *)two;
    return (a > b) - (a < b);
}

/* The seam part of ab_chart_links: the sampled vertices' checks and point
 * masses, then only the loops that hold samples, lowest start first, each
 * walked from the same vertex in the same direction as the full walk. */
static int ab_chart_links_loops(Arena_T arena, const AsmChart *c, const AsmBoundaryLoops *L, const int32_t *sample_at,
                                const int32_t *vertex, size_t count, AbLinks *links, double *point_mass)
{
    if (!L->valid) return -1;
    int32_t *starts = ARENA_ALLOC(arena,(count ? count : 1)*sizeof *starts); size_t n = 0;
    for (size_t i = 0; i < count; i++) {
        int32_t v = vertex[i];
        if (L->degree[v] != 2) return -1;
        point_mass[i] = .5*(ab_distance(c,v,L->boundary[2*(size_t)v])+ab_distance(c,v,L->boundary[2*(size_t)v+1]));
        starts[n++] = L->start[v];
    }
    qsort(starts,n,sizeof *starts,ab_int32_order);
    for (size_t i = 0; i < n; i++) {
        if (i && starts[i] == starts[i-1]) continue;
        int32_t start = starts[i], current = start, previous = -1, first_sample = -1, last_sample = -1;
        double length = 0, first_position = 0, last_position = 0;
        size_t observations = 0;
        do {
            int32_t sample = sample_at[current];
            if (sample >= 0) {
                if (first_sample < 0) { first_sample = sample; first_position = length; }
                if (last_sample >= 0) ab_add_link(links,last_sample,sample,length-last_position);
                last_sample = sample; last_position = length; observations++;
            }
            int32_t next = L->boundary[2*(size_t)current];
            if (next == previous) next = L->boundary[2*(size_t)current+1];
            length += ab_distance(c,current,next); previous = current; current = next;
        } while (current != start);
        if (observations >= 3) ab_add_link(links,last_sample,first_sample,length-last_position+first_position);
        else if (observations == 2) {
            links[first_sample].neighbour[0] = links[first_sample].neighbour[1] = -1;
            links[last_sample].neighbour[0] = links[last_sample].neighbour[1] = -1;
        }
    }
    return 0;
}

int AsmBoundary_weights_loops(Arena_T arena, const AsmChart *a, const AsmBoundaryLoops *la,
                              const AsmChart *b, const AsmBoundaryLoops *lb,
                              const AsmCorr *samples, size_t count, double *weights)
{
    if (!arena || !a || !b || !la || !lb || !samples || !weights || !count || count > INT32_MAX ||
        a->nv > INT32_MAX || b->nv > INT32_MAX || a->nf > INT32_MAX/6 || b->nf > INT32_MAX/6) return -1;
    const AsmChart *charts[2] = {a,b}; const AsmBoundaryLoops *loops[2] = {la,lb};
    for (int side = 0; side < 2; side++) {
        const AsmChart *c = charts[side];
        if (!c->nv || !c->nf || !c->xyz || !c->faces || !loops[side]->boundary) return -1;
    }
    memset(weights,0,count*sizeof *weights);
    Arena_Mark mark = Arena_save(arena); int status = -1;
    int32_t *at[2], *vertex[2]; AbLinks *links[2]; double *mass[2];
    for (int side = 0; side < 2; side++) {
        at[side] = ARENA_ALLOC(arena,charts[side]->nv*sizeof(int32_t));
        vertex[side] = ARENA_ALLOC(arena,count*sizeof(int32_t));
        links[side] = ARENA_CALLOC(arena,count,sizeof(AbLinks));
        mass[side] = ARENA_CALLOC(arena,count,sizeof(double));
        for (size_t v = 0; v < charts[side]->nv; v++) at[side][v] = -1;
        for (size_t i = 0; i < count; i++) {
            links[side][i].neighbour[0] = links[side][i].neighbour[1] = -1;
            int32_t v = side ? samples[i].vb : samples[i].va;
            if (samples[i].valid != 3 || v < 0 || (size_t)v >= charts[side]->nv || at[side][v] >= 0) goto done;
            at[side][v] = (int32_t)i; vertex[side][i] = v;
        }
        if (ab_chart_links_loops(arena,charts[side],loops[side],at[side],vertex[side],count,links[side],mass[side])) goto done;
    }
    for (size_t i = 0; i < count; i++) {
        for (int ka = 0; ka < 2; ka++) {
            int32_t j = links[0][i].neighbour[ka];
            if (j < 0 || (size_t)j <= i) continue;
            for (int kb = 0; kb < 2; kb++) if (links[1][i].neighbour[kb] == j &&
                fabs(links[0][i].length[ka]-links[1][i].length[kb]) <= ASM_CONT_RUN_MAX_LENGTH_ERROR) {
                double half = .25*(links[0][i].length[ka]+links[1][i].length[kb]);
                weights[i] += half; weights[j] += half;
            }
        }
    }
    for (size_t i = 0; i < count; i++) {
        if (weights[i] == 0) weights[i] = .5*(mass[0][i]+mass[1][i]);
        if (!(weights[i] > 0) || !isfinite(weights[i])) goto done;
    }
    status = 0;
done:
    Arena_restore(arena,mark); return status;
}

int AsmBoundary_weights(Arena_T arena, const AsmChart *a, const AsmChart *b,
                        const AsmCorr *samples, size_t count, double *weights)
{
    if (!arena || !a || !b) return -1;
    Arena_Mark mark = Arena_save(arena);
    AsmBoundaryLoops *la = AsmBoundary_loops(arena,a), *lb = AsmBoundary_loops(arena,b);
    int status = AsmBoundary_weights_loops(arena,a,la,b,lb,samples,count,weights);
    Arena_restore(arena,mark); return status;
}

int AsmBoundary_selftest(void)
{
    enum { N = 65, NV = 2*N, NF = 2*(N-1) };
    Arena_T arena = Arena_new();
    float xyz[2][3*NV]; int32_t faces[3*NF]; AsmCorr points[N];
    AsmChart charts[2] = {{0}}; AsmBoundaryStats stats = {0};
    int failures = 0;
    double expected = 0;
    for (size_t i = 0; i < N; i++) {
        double theta = 1.8*3.141592653589793*(double)i/(N-1);
        for (size_t side = 0; side < 2; side++) for (size_t edge = 0; edge < 2; edge++) {
            size_t v = edge*N+i;
            xyz[side][3*v] = (float)(64+5*cos(theta)); xyz[side][3*v+1] = (float)(64+5*sin(theta));
            xyz[side][3*v+2] = (float)(124+2.6*(double)edge+5.4*(double)side);
        }
        points[i] = (AsmCorr){0}; points[i].va = (int32_t)(N+i); points[i].vb = (int32_t)i; points[i].valid = 3;
        if (i+1 < N) {
            int32_t k = (int32_t)i;
            int32_t face[6] = {k,N+k,N+k+1,k,N+k+1,k+1}; memcpy(faces+6*i,face,sizeof face);
        }
    }
    for (size_t side = 0; side < 2; side++) {
        charts[side].xyz = xyz[side]; charts[side].faces = faces; charts[side].nv = NV; charts[side].nf = NF;
    }
    for (int32_t i = 1; i < N; i++) expected += ab_distance(&charts[0],i-1,i);
    if (AsmBoundary_order(arena,&charts[0],&charts[1],NULL,0,&stats) || stats.runs) failures++;
    /* Input order and a missing UV array cannot change source order. */
    for (size_t i = 0; i < N/2; i++) { AsmCorr p = points[i]; points[i] = points[N-1-i]; points[N-1-i] = p; }
    if (AsmBoundary_order(arena,&charts[0],&charts[1],points,N,&stats) || stats.supported_runs != 1 ||
        stats.ordered_samples != N || fabs(stats.support_length-expected) > 1e-6) failures++;
    for (size_t i = 0; i < N; i++) if (points[i].va != (int32_t)(N+i) || points[i].vb != (int32_t)i || points[i].run != 0) failures++;
    /* Removing one strip interval leaves two actual boundary components. */
    int32_t split_faces[3*(NF-2)]; size_t kept = 0;
    for (size_t f = 0; f < NF; f++) if (f/2 != 31) {
        memcpy(split_faces+3*kept,faces+3*f,3*sizeof *faces); kept++;
    }
    charts[0].faces = charts[1].faces = split_faces; charts[0].nf = charts[1].nf = kept;
    if (AsmBoundary_order(arena,&charts[0],&charts[1],points,N,&stats) || stats.supported_runs != 2 || stats.ordered_samples != N) failures++;
    if (points[31].run == points[32].run) failures++;
    /* Bad endpoint indices fail atomically, including existing run labels. */
    points[0].vb = NV; AsmCorr original[N]; memcpy(original,points,sizeof points);
    if (!AsmBoundary_order(arena,&charts[0],&charts[1],points,N,&stats) || memcmp(original,points,sizeof points)) failures++;
    /* Seam weights from per-chart loops equal the direct per-seam walk, bit
     * for bit: one and two boundary loops, all or some samples, a bad
     * endpoint, and one loops object reused by two different seams. */
    for (int variant = 0; variant < 4; variant++) {
        AsmCorr use[N]; size_t n = 0;
        charts[0].faces = charts[1].faces = variant == 0 ? faces : split_faces;
        charts[0].nf = charts[1].nf = variant == 0 ? (size_t)NF : kept;
        for (size_t i = 0; i < N; i++) if (variant != 2 || i%3 == 0) {
            use[n] = (AsmCorr){0}; use[n].va = (int32_t)(N+i); use[n].vb = (int32_t)i; use[n].valid = 3; n++;
        }
        if (variant == 3) use[0].vb = NV;
        double direct[N], cached[N], reused[N];
        int rd = ab_weights_direct(arena,&charts[0],&charts[1],use,n,direct);
        int rl = AsmBoundary_weights(arena,&charts[0],&charts[1],use,n,cached);
        if (rd != rl || (!rd && memcmp(direct,cached,n*sizeof *direct)) || (variant == 3) != (rd != 0)) failures++;
        Arena_Mark mark = Arena_save(arena);
        AsmBoundaryLoops *la = AsmBoundary_loops(arena,&charts[0]), *lb = AsmBoundary_loops(arena,&charts[1]);
        for (int again = 0; again < 2; again++) {
            size_t m = again ? n/2 : n;
            int r1 = ab_weights_direct(arena,&charts[0],&charts[1],use,m,direct);
            int r2 = AsmBoundary_weights_loops(arena,&charts[0],la,&charts[1],lb,use,m,reused);
            if (r1 != r2 || (!r1 && memcmp(direct,reused,m*sizeof *direct))) failures++;
        }
        Arena_restore(arena,mark);
    }
    Arena_dispose(&arena);
    fprintf(stderr,"  assembly source boundary order: %s (%d failures)\n",failures ? "FAIL" : "ok",failures);
    return failures;
}
