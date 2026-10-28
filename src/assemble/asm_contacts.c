#ifdef _WIN32
#include <intrin.h>
#endif
#include "asm_contacts.h"
#include "../common/pipeline_constants.h"
#include "../common/union_find.h"
#include "../common/ves_pool.h"
#include "../common/ves_platform.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double AsmContacts_material_tolerance(void) { return ASM_MATERIAL_CONTACT_TOLERANCE_VOX; }
const char *AsmContacts_material_policy(void) { return "original_boundary_overlap_width_v1"; }

#include "../common/ves_omp.h"

/* One leaf-leaf node pair, recorded in exactly the order the serial walk
 * visits them, and one contact it produced.  `lp`/`seq` restore that order
 * after a parallel pass, so the reduction adds the same numbers in the same
 * sequence as the serial walk and the answers are bit-identical. */
typedef struct AcLeafPair { int32_t a, b; } AcLeafPair;
typedef struct AcHit { int32_t lp, seq; size_t a, b; double area, raw_area, ga[6], gb[6]; } AcHit;
#define AC_UNMEASURED_SAMPLE 8   /* the first unmeasured pairs a walk keeps for its ledger */
typedef struct AcBuf { AcHit *hit; size_t n, cap; int failed; size_t unmeasured, first_a, first_b; } AcBuf;

/* WALK MEMO.  The sequence of hits a walk commits is a function of the
 * tree's partition, the coordinates, the active mask, the frame labels, the
 * walk's mode (certificate measure or soft objective) and whether gradients
 * are requested -- nothing else.  The solver re-measures unchanged
 * coordinates constantly (restoring the best iterate, the weighted second
 * pass of an objective, repeated merit checks), so a complete masked walk
 * keeps its committed sequence, one slot per kind.  A request identical in
 * every input replays that sequence through ac_commit in the same order:
 * each tally, gradient entry and callback receives exactly the values and
 * calls a fresh walk would give it.  The boxes are still refitted for the
 * request, so the index state after a replay is the fresh walk's too. */
typedef struct AcMemo {
    int valid;
    unsigned version;            /* partition generation the hits were enumerated on */
    const int32_t *frames;
    double *uv;                  /* 2*nv: the exact coordinates, bytewise */
    uint8_t *mask;               /* nf: the exact active mask, bytewise */
    AcHit *hit; size_t n, cap;   /* committed order */
    size_t leaf_pairs;
} AcMemo;
enum { AC_MEMO_KINDS = 4 };      /* (objective mode) x (gradient) */

typedef struct AcNode {
    double lo[2], hi[2];
    size_t begin, end;
    int32_t left, right;
    int active;
} AcNode;
struct AsmContacts {
    size_t nv, nf, nodes, capacity;
    const int32_t *faces;
    const int32_t *frames;
    const uint8_t *active_faces;
    size_t *order;
    AcNode *tree;
    struct AcMaterial *material;
    struct AcBandWork *band_work;
    int n_band_work;
    double tolerance;
    const double *uv;
    double *gradient;
    AsmContactVisit visit;
    void *context;
    AsmContactStats stats;
    double objective_area;
    int objective_mode;
    AsmContactWeight objective_weight;
    void *objective_context;
    int failed;
    size_t unmeasured_ab[2*AC_UNMEASURED_SAMPLE], n_unmeasured_ab;   /* the first unmeasured pairs of the walk */
    double *built;        /* the geometry this tree was built on (2*nv), or NULL */
    int     have_built;
    int     threads;      /* <= 1 keeps the exact serial walk */
    Arena_T arena;        /* the field's arena: grow-only buffers live and die with it */
    AcLeafPair *lp; size_t n_lp, cap_lp;
    AcBuf  *buf; int n_buf;
    AcHit  *merged; size_t cap_merged;
    size_t n_merged;         /* hits of the last walk, in commit order (threaded, or serial while recording) */
    int    recording;        /* the serial walk appends its committed hits to merged */
    unsigned version;        /* partition generation: every build bumps it */
    AcMemo *memo;            /* AC_MEMO_KINDS slots, allocated on first masked walk */
    int    memo_off;         /* controls only: measure every request afresh */
    int    memo_verify;      /* re-measure each replay dry and compare (ASM_CONTACT_MEMO_VERIFY) */
    int    dry;              /* record hits without committing them */
    size_t memo_replays, memo_walks, memo_mismatches;
    size_t leaves, budget;   /* leaf count of the current tree; leaf-pair budget of one walk (0 = unbounded) */
    int    budget_explicit;  /* the caller set the budget: rebuilds keep it */
    int    over_budget;      /* the current walk stopped at the budget */
    uint8_t *probe_nodes;
    size_t probe_capacity;
    const uint8_t *probe_exclude;
    struct AcCharts *charts;   /* embedded-chart pruning (AsmContacts_set_charts), or NULL */
    const double *center;      /* [2*nf] face centroids during ac_build_all, else NULL */
};

/* Bump arenas cannot realloc: grow by allocating the larger block and copying.
 * The old block is dead until the arena is disposed, which for a trial field
 * is the end of its attempt; nothing here is malloc'd, so a field that goes
 * away with its arena leaks nothing. */
static void *ac_grow(Arena_T arena, const void *old, size_t old_bytes, size_t new_bytes)
{
    void *grown = ARENA_ALLOC(arena, new_bytes);
    if (old_bytes) memcpy(grown, old, old_bytes);
    return grown;
}

static double ac_cross(double ax, double ay, double bx, double by)
{ return ax*by-ay*bx; }

/* fmin/fmax without the CRT call on the common path. MSVC calls both out of
 * line, and in the pair tests they were 30% of the busy samples of parallel
 * repair corrections on a 10^3 block (2026-09-26). Unequal ordered inputs take the
 * comparison; equal or unordered ones (signed zeros, NaN) still reach the
 * CRT function, so every result is bit-identical to fmin/fmax. */
/* Equal inputs with identical bits (a vertex shared by faces of one leaf meets
 * its own coordinate) are their own fmin/fmax, so only signed zeros and NaN
 * reach the CRT (2026-10-05: fmin+fmax were 1.6% of a large assembly's busy samples). */
static inline int ac_same(double a, double b)
{
    uint64_t x, y; memcpy(&x,&a,sizeof x); memcpy(&y,&b,sizeof y);
    return a == b && x == y;
}
static inline double ac_min(double a, double b) { return a < b ? a : b < a ? b : ac_same(a,b) ? a : fmin(a,b); }
static inline double ac_max(double a, double b) { return a > b ? a : b > a ? b : ac_same(a,b) ? a : fmax(a,b); }

#include "asm_contacts_band.inc"
#include "asm_contacts_material.inc"

static void ac_positive(const double in[6], double out[6], int order[3], const double origin[2])
{
    int reverse = ac_cross(in[2]-in[0],in[3]-in[1],in[4]-in[0],in[5]-in[1]) < 0;
    order[0] = 0; order[1] = reverse ? 2 : 1; order[2] = reverse ? 1 : 2;
    for (int k = 0; k < 3; k++) for (int d = 0; d < 2; d++) out[2*k+d] = in[2*order[k]+d]-origin[d];
}

/* Strict separating-axis predicate: tolerance is distance, not area.
 * A shared edge or vertex is not an overlapping material region. */
static int ac_overlap(const double a[6], const double b[6], double tolerance)
{
    for (int d = 0; d < 2; d++) {
        double alo = ac_min(a[d],ac_min(a[2+d],a[4+d])), ahi = ac_max(a[d],ac_max(a[2+d],a[4+d]));
        double blo = ac_min(b[d],ac_min(b[2+d],b[4+d])), bhi = ac_max(b[d],ac_max(b[2+d],b[4+d]));
        if (ac_min(ahi,bhi)-ac_max(alo,blo) <= tolerance) return 0;
    }
    for (int side = 0; side < 2; side++) {
        const double *p = side ? b : a;
        for (int k = 0; k < 3; k++) {
            int j = (k+1)%3;
            double nx = p[2*j+1]-p[2*k+1], ny = p[2*k]-p[2*j];
            double alo = DBL_MAX, ahi = -DBL_MAX, blo = DBL_MAX, bhi = -DBL_MAX;
            for (int v = 0; v < 3; v++) {
                double x = nx*a[2*v]+ny*a[2*v+1], y = nx*b[2*v]+ny*b[2*v+1];
                alo = ac_min(alo,x); ahi = ac_max(ahi,x); blo = ac_min(blo,y); bhi = ac_max(bhi,y);
            }
            double gap = ac_min(ahi,bhi)-ac_max(alo,blo);
            /* A nonpositive gap separates under any nonnegative tolerance, so
             * the edge length (a CRT hypot call) is needed only past it; a
             * degenerate edge projects to a zero gap and still separates. */
            if (tolerance >= 0 && gap <= 0) return 0;
            /* The edge length matters only near the threshold.  A negative
             * tolerance never separates a positive gap; a nonnegative one
             * cannot separate a gap beyond tolerance*(|nx|+|ny|) (>= the
             * length, with a margin for the roundoff of both products), so
             * neither needs the hypot call -- the same decision either way.
             * Touching neighbours reach here twice per objective pair. */
            if (gap > 0 && (tolerance < 0 || gap > tolerance*(fabs(nx)+fabs(ny))*(1+16*DBL_EPSILON))) continue;
            double norm = hypot(nx,ny);
            if (!(norm > 0)) return 0;
            if (gap <= tolerance*norm) return 0;
        }
    }
    return 1;
}

/* Sutherland-Hodgman audit area. The derivative below uses a separate
 * boundary integral so its own arithmetic cannot certify its objective. */
static double ac_clip_polygon(const double a[6], const double b[6],double *clipped,int *vertices)
{
    double poly[24], next[24]; memcpy(poly,a,6*sizeof(double)); int n = 3;
    for (int k = 0; k < 3 && n; k++) {
        int j = (k+1)%3, m = 0;
        double ex = b[2*j]-b[2*k], ey = b[2*j+1]-b[2*k+1];
        for (int v = 0; v < n; v++) {
            int w = (v+1)%n;
            double da = ac_cross(ex,ey,poly[2*v]-b[2*k],poly[2*v+1]-b[2*k+1]);
            double db = ac_cross(ex,ey,poly[2*w]-b[2*k],poly[2*w+1]-b[2*k+1]);
            if (da >= 0) { next[2*m] = poly[2*v]; next[2*m+1] = poly[2*v+1]; m++; }
            if ((da < 0 && db > 0) || (da > 0 && db < 0)) {
                double t = da/(da-db);
                next[2*m] = poly[2*v]+t*(poly[2*w]-poly[2*v]);
                next[2*m+1] = poly[2*v+1]+t*(poly[2*w+1]-poly[2*v+1]); m++;
            }
        }
        memcpy(poly,next,2*(size_t)m*sizeof(double)); n = m;
    }
    double area = 0;
    for (int k = 1; k+1 < n; k++)
        area += ac_cross(poly[2*k]-poly[0],poly[2*k+1]-poly[1],poly[2*(k+1)]-poly[0],poly[2*(k+1)+1]-poly[1]);
    if(clipped)memcpy(clipped,poly,2*(size_t)n*sizeof(double));
    if(vertices)*vertices=n;
    return .5*fabs(area);
}

/* Near-parallel boundary clipping divides two small determinants. Computing
 * each as a difference of rounded products can give incompatible intersection
 * parameters on opposite triangle banks, even producing a net translation
 * force. Retain the product residuals; this changes arithmetic accuracy, not
 * the contact tolerance or the treatment of coincident edges. */
static double ac_boundary_cross(double x,double y,double u,double v)
{
    double p=x*v,q=y*u,head=p-q;
    double bv=p-head,av=head+bv;
    double difference_tail=(p-av)+(bv-q);
    return head+(difference_tail+fma(x,v,-p)-fma(y,u,-q));
}

static void ac_boundary_gradient(const double a[6], const double b[6], double g[6])
{
    memset(g,0,6*sizeof(double));
    for (int k = 0; k < 3; k++) {
        int j = (k+1)%3;
        double ex = a[2*j]-a[2*k], ey = a[2*j+1]-a[2*k+1], lo = 0, hi = 1;
        int possible = 1, coincident = 0;
        for (int v = 0; v < 3; v++) {
            int w = (v+1)%3;
            double bx = b[2*w]-b[2*v], by = b[2*w+1]-b[2*v+1];
            double start = ac_boundary_cross(bx,by,a[2*k]-b[2*v],a[2*k+1]-b[2*v+1]);
            double slope = ac_boundary_cross(bx,by,ex,ey);
            if (slope > 0) lo = ac_max(lo,-start/slope);
            else if (slope < 0) hi = ac_min(hi,-start/slope);
            else { if (start < 0) possible = 0; if (start == 0) coincident = 1; }
        }
        double width = possible ? ac_max(hi-lo,0) : 0;
        if (coincident) width *= .5;
        double middle = .5*(lo+hi);
        g[2*k] += ey*width*(1-middle); g[2*k+1] -= ex*width*(1-middle);
        g[2*j] += ey*width*middle; g[2*j+1] -= ex*width*middle;
    }
}

double AsmContacts_pair(const double in_a[6], const double in_b[6], double tolerance,
                        double ga[6], double gb[6])
{
    if (ga) memset(ga,0,6*sizeof(double)); if (gb) memset(gb,0,6*sizeof(double));
    if (!isfinite(tolerance) || tolerance < 0) return NAN;
    for (int k = 0; k < 6; k++) if (!isfinite(in_a[k]) || !isfinite(in_b[k])) return NAN;
    double a[6], b[6], origin[2] = {in_a[0],in_a[1]}; int ia[3], ib[3];
    ac_positive(in_a,a,ia,origin); ac_positive(in_b,b,ib,origin);
    if (!ac_overlap(a,b,tolerance)) return 0;
    double area = ac_clip_polygon(a,b,NULL,NULL);
    if (area > 0 && (ga || gb)) {
        double g[6];
        if (ga) { ac_boundary_gradient(a,b,g); for (int k = 0; k < 3; k++) for (int d = 0; d < 2; d++) ga[2*ia[k]+d] = g[2*k+d]; }
        if (gb) { ac_boundary_gradient(b,a,g); for (int k = 0; k < 3; k++) for (int d = 0; d < 2; d++) gb[2*ib[k]+d] = g[2*k+d]; }
    }
    return area;
}

#include "asm_contacts_query.inc"

static double ac_center(const AsmContacts *q, size_t face, int axis)
{
    const int32_t *f = q->faces+3*face;
    return (q->uv[2*(size_t)f[0]+axis]+q->uv[2*(size_t)f[1]+axis]+q->uv[2*(size_t)f[2]+axis])/3;
}

/* The build's centroids (ac_build_all computes them once with ac_center, so
 * every comparison sees the values it used to recompute per call). */
static int ac_less(const AsmContacts *q, size_t a, size_t b, int axis)
{
    double x = q->center[2*a+(size_t)axis], y = q->center[2*b+(size_t)axis];
    return x < y || (x == y && a < b);
}

static void ac_select(AsmContacts *q, size_t begin, size_t end, size_t nth, int axis)
{
    while (end-begin > 1) {
        size_t pivot = q->order[begin+(end-begin)/2], i = begin, j = end-1;
        for (;;) {
            while (ac_less(q,q->order[i],pivot,axis)) i++;
            while (ac_less(q,pivot,q->order[j],axis)) j--;
            if (i >= j) break;
            size_t t = q->order[i]; q->order[i++] = q->order[j]; q->order[j--] = t;
        }
        size_t cut = i;
        if (cut == begin) cut++; /* unique face-id tie breaker guarantees progress */
        if (nth < cut) end = cut; else begin = cut;
    }
}

/* Returns 0 when the leaf reads a non-finite coordinate. */
static int ac_leaf_box(AsmContacts *q, AcNode *node)
{
    int finite = 1;
    node->lo[0] = node->lo[1] = DBL_MAX; node->hi[0] = node->hi[1] = -DBL_MAX;
    node->active = 0;
    for (size_t j = node->begin; j < node->end; j++) {
        if (!q->active_faces || q->active_faces[q->order[j]]) node->active = 1;
        const int32_t *f = q->faces+3*q->order[j];
        for (int k = 0; k < 3; k++) for (int d = 0; d < 2; d++) {
            double x = q->uv[2*(size_t)f[k]+d];
            if (!isfinite(x)) finite = 0;
            node->lo[d] = ac_min(node->lo[d],x); node->hi[d] = ac_max(node->hi[d],x);
        }
    }
    return finite;
}

/* The caller's mask decides which faces may have moved, so recompute it on
 * every leaf; it is a handful of byte reads. */
static int ac_leaf_active(const AsmContacts *q, const AcNode *node)
{
    if (!q->active_faces) return 1;
    for (size_t j = node->begin; j < node->end; j++) if (q->active_faces[q->order[j]]) return 1;
    return 0;
}

/* Nodes in the subtree over a range of s faces.  The split is a pure
 * function of the range size (mid = begin + s/2, leaves at 16), so the
 * preorder id of every node is known before it is built: the left child is
 * id+1 and the right child id+1+count(left).  Sizes at one depth differ by
 * at most one, so two slots per depth suffice and the count is O(log s). */
static size_t ac_count(size_t s)
{
    size_t total = 0, size[2] = { s, 0 }, count[2] = { 1, 0 };
    while (count[0] || count[1]) {
        size_t nsize[2] = { 0, 0 }, ncount[2] = { 0, 0 };
        for (int i = 0; i < 2; i++) {
            if (!count[i]) continue;
            total += count[i];
            if (size[i] <= 16) continue;
            size_t child[2] = { size[i]/2, size[i] - size[i]/2 };
            for (int c = 0; c < 2; c++) {
                int slot = -1;
                for (int j = 0; j < 2; j++) if (ncount[j] && nsize[j] == child[c]) slot = j;
                if (slot < 0) for (int j = 0; j < 2; j++) if (!ncount[j]) { slot = j; nsize[j] = child[c]; break; }
                if (slot < 0) return 0;
                ncount[slot] += count[i];
            }
        }
        memcpy(size, nsize, sizeof size); memcpy(count, ncount, sizeof count);
    }
    return total;
}

/* Build the subtree over [begin,end) at the preassigned preorder id. */
static int ac_build_at(AsmContacts *q, size_t begin, size_t end, int32_t id)
{
    AcNode *node = q->tree+id;
    node->begin = begin; node->end = end; node->left = node->right = -1; ac_leaf_box(q,node);
    if (end-begin > 16) {
        size_t mid = begin+(end-begin)/2;
        int axis = node->hi[1]-node->lo[1] > node->hi[0]-node->lo[0];
        ac_select(q,begin,end,mid,axis);
        node->left = id+1; node->right = id+1+(int32_t)ac_count(mid-begin);
        if ((size_t)node->right >= q->nodes) return -1;
        if (ac_build_at(q,begin,mid,node->left) < 0 || ac_build_at(q,mid,end,node->right) < 0) return -1;
    }
    return 0;
}

/* The top of the tree is built serially down to subtrees of `grain` faces,
 * which are then built in parallel: every subtree owns its id block and its
 * slice of the face order, so the result is the serial tree byte for byte
 * (the control in the selftest holds it to that). */
typedef struct AcTask { size_t begin, end; int32_t id; } AcTask;
#define AC_PARALLEL_MIN_FACES 4096
#define AC_MAX_TASKS 1024
static int ac_build_top(AsmContacts *q, size_t begin, size_t end, int32_t id, size_t grain, AcTask *tasks, int *ntasks)
{
    if (end-begin <= 16) return ac_build_at(q,begin,end,id);
    if (end-begin <= grain || *ntasks >= AC_MAX_TASKS) {
        if (*ntasks < AC_MAX_TASKS) { tasks[*ntasks].begin = begin; tasks[*ntasks].end = end; tasks[*ntasks].id = id; (*ntasks)++; return 0; }
        return ac_build_at(q,begin,end,id);
    }
    AcNode *node = q->tree+id;
    node->begin = begin; node->end = end; node->left = node->right = -1; ac_leaf_box(q,node);
    size_t mid = begin+(end-begin)/2;
    int axis = node->hi[1]-node->lo[1] > node->hi[0]-node->lo[0];
    ac_select(q,begin,end,mid,axis);
    node->left = id+1; node->right = id+1+(int32_t)ac_count(mid-begin);
    if ((size_t)node->right >= q->nodes) return -1;
    if (ac_build_top(q,begin,mid,node->left,grain,tasks,ntasks) < 0) return -1;
    return ac_build_top(q,mid,end,node->right,grain,tasks,ntasks);
}

static int ac_build_tasks(AsmContacts *q)
{
    int nt = q->threads;
    if (nt <= 1 || q->nf < AC_PARALLEL_MIN_FACES) return ac_build_at(q,0,q->nf,0);
    AcTask *tasks = malloc(AC_MAX_TASKS * sizeof *tasks); int ntasks = 0, bad = 0, t = 0;
    if (!tasks) return -1;
    size_t grain = q->nf/(size_t)(8*nt); if (grain < 1024) grain = 1024;
    if (ac_build_top(q,0,q->nf,0,grain,tasks,&ntasks) < 0) { free(tasks); return -1; }
#pragma omp parallel for schedule(dynamic,1) num_threads(nt) reduction(|:bad)
    for (t = 0; t < ntasks; t++) if (ac_build_at(q,tasks[t].begin,tasks[t].end,tasks[t].id) < 0) bad |= 1;
    free(tasks);
    return bad ? -1 : 0;
}

/* The quickselect compared two centroids per step, each recomputed from three
 * coordinates and a division (ac_select was 2.1% of a large assembly's busy samples,
 * 2026-10-06); they are computed once per build instead, with the same
 * expression, so the comparisons and the tree are unchanged. */
static int ac_build_all(AsmContacts *q)
{
    double *center = malloc((q->nf ? 2*q->nf : 1)*sizeof *center);
    if (!center) return -1;
    for (size_t f = 0; f < q->nf; f++) { center[2*f] = ac_center(q,f,0); center[2*f+1] = ac_center(q,f,1); }
    q->center = center;
    int rc = ac_build_tasks(q);
    q->center = NULL;
    free(center);
    return rc;
}

AsmContacts *AsmContacts_new(Arena_T arena, size_t nv, const int32_t *faces,
                             size_t nf, const double *uv, double tolerance)
{ return AsmContacts_new_threads(arena,nv,faces,nf,uv,tolerance,1); }

AsmContacts *AsmContacts_new_threads(Arena_T arena, size_t nv, const int32_t *faces,
                                     size_t nf, const double *uv, double tolerance, int threads)
{ return AsmContacts_new_source(arena,nv,faces,nf,uv,tolerance,threads,NULL); }

AsmContacts *AsmContacts_new_source(Arena_T arena, size_t nv, const int32_t *faces,
                                    size_t nf, const double *uv, double tolerance, int threads,
                                    const int32_t *material_vertex)
{
    if (!arena || !nv || !nf || !faces || !uv || !isfinite(tolerance) || tolerance < 0 ||
        nv > INT32_MAX || nf > SIZE_MAX/3 || nf > SIZE_MAX/sizeof(size_t) ||
        nf/4+64 > INT32_MAX || nf/4+64 > SIZE_MAX/sizeof(AcNode)) return NULL;
    for (size_t k = 0; k < 3*nf; k++) if (faces[k] < 0 || (size_t)faces[k] >= nv) return NULL;
    for (size_t k = 0; k < 2*nv; k++) if (!isfinite(uv[k])) return NULL;
    if (material_vertex) {
        for (size_t v = 0; v < nv; v++) if (material_vertex[v] < 0 || (size_t)material_vertex[v] >= nv) return NULL;
        for (size_t v = 0; v < nv; v++) if (material_vertex[material_vertex[v]] != material_vertex[v]) return NULL;
    }
    AsmContacts *q = ARENA_CALLOC(arena,1,sizeof *q);
    q->nv = nv; q->nf = nf; q->faces = faces; q->uv = uv; q->tolerance = tolerance; q->threads = threads;
    q->memo_verify = ASM_CONTACT_MEMO_VERIFY;
    q->capacity = nf/4+64;
    q->tree = ARENA_ALLOC(arena,q->capacity*sizeof *q->tree);
    q->order = ARENA_ALLOC(arena,nf*sizeof *q->order);
    q->built = ARENA_ALLOC(arena,2*nv*sizeof *q->built);
    q->arena = arena;
    q->material=acm_new(arena,nv,faces,nf,uv,material_vertex);
    if(!q->material)return NULL;
    AsmContacts_threads(q,threads);
    return AsmContacts_rebuild(q,uv) ? NULL : q;
}

/* The walk budget of a tree with this many leaves.  A single-cover sheet
 * visits ~8 overlapping leaves per leaf (the boxes are closed intervals), a
 * double cover ~30; the 21x5x5 admission of 2026-09-17 was at ~4,000 when it
 * asked the arena for 256 GiB of leaf pairs (a tree built with the candidate
 * at the uv origin and refitted after the seam fit moved it). */
static size_t ac_default_budget(size_t leaves)
{
    size_t per = (size_t)ASM_CONTACT_MAX_LEAF_PAIRS_PER_LEAF * leaves;
    return per < (size_t)ASM_CONTACT_MAX_LEAF_PAIRS ? per : (size_t)ASM_CONTACT_MAX_LEAF_PAIRS;
}

/* Rebuilding is O(n log n) over the WHOLE field and the repair asks for it
 * once per admission attempt, but an attempt that is refused restores the
 * incumbent geometry, so most of those rebuilds reproduce the tree they
 * already had.  A tree is a function of the geometry it was built on: when
 * that geometry is unchanged the existing partition IS the rebuild, and the
 * refit inside every measure re-establishes the boxes.  The comparison is one
 * linear pass against an O(n log n) median selection at every node. */
static void acc_nodes(AsmContacts *q);   /* asm_contacts_charts.inc */
int AsmContacts_rebuild(AsmContacts *q,const double *uv)
{
    if(!q || !uv)return -1;
    q->probe_exclude=NULL;
    if(q->have_built && q->nodes && !memcmp(q->built,uv,2*q->nv*sizeof *q->built))return AsmContacts_refit(q,uv);
    for(size_t k=0;k<2*q->nv;k++)if(!isfinite(uv[k]))return -1;
    q->uv=uv;q->nodes=ac_count(q->nf);
    if(!q->nodes || q->nodes>q->capacity)return -1;
    for(size_t f=0;f<q->nf;f++)q->order[f]=f;
    q->version++;   /* even a failed build leaves no partition a memo may assume */
    if(ac_build_all(q)<0)return -1;
    acc_nodes(q);
    if(!acm_refit(q->material,uv,NULL))return -1;
    if(q->built){memcpy(q->built,uv,2*q->nv*sizeof *q->built);q->have_built=1;}
    q->leaves=0;for(size_t i=0;i<q->nodes;i++)q->leaves+=q->tree[i].left<0;
    if(!q->budget_explicit)q->budget=ac_default_budget(q->leaves);
    return 0;
}

void AsmContacts_leaf_pair_budget(AsmContacts *q, size_t budget)
{ if (q) { q->budget = budget; q->budget_explicit = 1; } }
size_t AsmContacts_budget(const AsmContacts *q) { return q ? q->budget : 0; }
size_t AsmContacts_leaves(const AsmContacts *q) { return q ? q->leaves : 0; }

#include "asm_contacts_charts.inc"

static int ac_boxes(const AcNode *a, const AcNode *b)
{
    return a->lo[0] <= b->hi[0] && b->lo[0] <= a->hi[0] &&
           a->lo[1] <= b->hi[1] && b->lo[1] <= a->hi[1];
}

/* Original shared-edge neighbors with opposite third vertices stay disjoint
 * under coherent vertex motion. Their boundary derivatives cancel exactly;
 * folded neighbors on the same side still require complete measurement. */
static int ac_shared_edge_clear(const AsmContacts *q,size_t a,size_t b)
{
    const int32_t *fa=q->faces+3*a,*fb=q->faces+3*b;int shared[2],count=0,oa=-1,ob=-1;
    for(int i=0;i<3;i++){
        int found=0;for(int j=0;j<3;j++)found|=fa[i]==fb[j];
        if(found){if(count==2)return 0;shared[count++]=fa[i];}else oa=fa[i];
    }
    if(count!=2 || oa<0)return 0;
    for(int i=0;i<3;i++)if(fb[i]!=shared[0] && fb[i]!=shared[1])ob=fb[i];
    if(ob<0)return 0;
    const double *p=q->uv+2*(size_t)shared[0],*r=q->uv+2*(size_t)shared[1];
    double ax=q->uv[2*(size_t)oa]-p[0],ay=q->uv[2*(size_t)oa+1]-p[1];
    double bx=q->uv[2*(size_t)ob]-p[0],by=q->uv[2*(size_t)ob+1]-p[1];
    double x=ac_cross(r[0]-p[0],r[1]-p[1],ax,ay),y=ac_cross(r[0]-p[0],r[1]-p[1],bx,by);
    return isfinite(x) && isfinite(y) && ((x<0 && y>0) || (x>0 && y<0));
}

/* Two triangles whose boxes are apart by more than 256 ulp of the pair's
 * largest coordinate: four times the roundoff ac_material_pair allows (64 ulp
 * of the same value), so its separating-axis test would find them apart in
 * either direction and return no area, no raw area and no derivative. Most
 * pairs of two overlapping leaves are such pairs; this skips the pair setup
 * for them. A non-finite coordinate leaves the decision to the full test. */
static int ac_separated(const double x[6], const double y[6])
{
    double sum = 0, scale = 1, lo[2][2], hi[2][2];
    for (int k = 0; k < 6; k++) sum += x[k]+y[k];
    if (!isfinite(sum)) return 0;
    for (int d = 0; d < 2; d++) {
        lo[0][d] = hi[0][d] = x[d]; lo[1][d] = hi[1][d] = y[d];
        for (int k = 1; k < 3; k++) {
            double u = x[2*k+d], v = y[2*k+d];
            if (u < lo[0][d]) lo[0][d] = u; else if (u > hi[0][d]) hi[0][d] = u;
            if (v < lo[1][d]) lo[1][d] = v; else if (v > hi[1][d]) hi[1][d] = v;
        }
        for (int t = 0; t < 2; t++) {
            double m = fabs(lo[t][d]), n = fabs(hi[t][d]);
            if (m > scale) scale = m;
            if (n > scale) scale = n;
        }
    }
    double margin = 256*DBL_EPSILON*scale;
    for (int d = 0; d < 2; d++) if (lo[1][d]-hi[0][d] > margin || lo[0][d]-hi[1][d] > margin) return 1;
    return 0;
}

/* Diagnostic pair-test counters (process-wide, interlocked; AsmContacts_pair_counts).
 * Compiled out by default: interlocked adds in the pair test cost more than
 * most pairs (E18 profiling, 2026-10-05). */
#ifndef AC_PAIR_COUNTERS
#define AC_PAIR_COUNTERS 0
#endif
static volatile long long ac_count_tests,ac_count_separated,ac_count_edge_clear,ac_count_vertex,ac_count_edge,ac_count_material,ac_count_hits;
static volatile long long ac_count_leafpairs,ac_count_culls,ac_count_selfleaf;
static volatile long long ac_count_vhit_area,ac_count_vhit_raw,ac_count_vhit_derivative,ac_count_ehit,ac_count_nhit;
static void ac_count_add(volatile long long *p,long long v)
{
    if(!AC_PAIR_COUNTERS)return;
#ifdef _WIN32
    _InterlockedExchangeAdd64((volatile __int64 *)p,v);
#else
    __sync_fetch_and_add(p,v);
#endif
}
void AsmContacts_pair_counts(long long out[7])
{
    out[0]=ac_count_tests;out[1]=ac_count_separated;out[2]=ac_count_edge_clear;out[3]=ac_count_vertex;
    out[4]=ac_count_edge;out[5]=ac_count_material;out[6]=ac_count_hits;
}
void AsmContacts_hit_classes(long long out[5])
{out[0]=ac_count_vhit_area;out[1]=ac_count_vhit_raw;out[2]=ac_count_vhit_derivative;out[3]=ac_count_ehit;out[4]=ac_count_nhit;}
void AsmContacts_leaf_counts(long long out[3]){out[0]=ac_count_leafpairs;out[1]=ac_count_culls;out[2]=ac_count_selfleaf;}

/* Geometry only: no shared state is touched, so this is what a worker runs.
 * Returns 1 when the pair contacts, -1 when its area is not finite. */
static int ac_pair_area(const AsmContacts *q,size_t a,size_t b,double *area_out,double *raw_out,double *ga,double *gb,AcBandWork *work)
{
    if (q->frames && q->frames[a] != q->frames[b]) return 0;
    if (q->active_faces && !q->active_faces[a] && !q->active_faces[b]) return 0;
    if (acc_skip_faces(q,a,b)) return 0;
    ac_count_add(&ac_count_tests,1);
    int shared=0;
    if(AC_PAIR_COUNTERS){   const int32_t *fa=q->faces+3*a,*fb=q->faces+3*b;
        for(int i=0;i<3;i++)for(int j=0;j<3;j++)shared+=fa[i]==fb[j];
        if(shared==1)ac_count_add(&ac_count_vertex,1);else if(shared==2)ac_count_add(&ac_count_edge,1);
    }
    double x[6],y[6];ac_face_xy(q,q->uv,a,x);ac_face_xy(q,q->uv,b,y);
    if(ac_separated(x,y)){ac_count_add(&ac_count_separated,1);return 0;}
    if(q->objective_mode && ac_shared_edge_clear(q,a,b)){ac_count_add(&ac_count_edge_clear,1);return 0;}
    ac_count_add(&ac_count_material,1);
    double area=ac_material_pair(work,q->material,a,q->uv,x,q->material,b,q->uv,y,q->tolerance,ga,gb,q->objective_mode?raw_out:NULL);
    if (!isfinite(area)) return -1;
    int derivative=0;
    if(q->objective_mode && ga && gb)for(int k=0;k<6;k++)derivative|=ga[k]!=0 || gb[k]!=0;
    if (!(area > 0) && !(q->objective_mode && (*raw_out>0 || derivative))) return 0;
    *area_out = area;ac_count_add(&ac_count_hits,1);
    if(!AC_PAIR_COUNTERS){}
    else if(shared==1){
        if(area>0)ac_count_add(&ac_count_vhit_area,1);
        else if(q->objective_mode && *raw_out>0)ac_count_add(&ac_count_vhit_raw,1);
        else ac_count_add(&ac_count_vhit_derivative,1);
    }else if(shared>=2)ac_count_add(&ac_count_ehit,1);else ac_count_add(&ac_count_nhit,1);
    return 1;
}

/* Applies one contact to the shared tallies, in the walk's own order. */
static void ac_commit(AsmContacts *q, size_t a, size_t b, double area,double raw_area,const double *ga,const double *gb)
{
    const int32_t *fa = q->faces+3*a, *fb = q->faces+3*b;
    if(area>0){q->stats.pairs++;q->stats.area+=area;}
    double weight=q->objective_weight?q->objective_weight(q->objective_context,a,b,raw_area):1;
    if(!isfinite(weight) || weight<0){q->failed=1;return;}
    q->objective_area+=weight*raw_area;
    if (q->gradient) for (int k = 0; k < 3; k++) for (int d = 0; d < 2; d++) {
        q->gradient[2*(size_t)fa[k]+d] += weight*ga[2*k+d]; q->gradient[2*(size_t)fb[k]+d] += weight*gb[2*k+d];
    }
    if (q->visit && area>0) q->visit(q->context,a,b,area);
}

/* A pair whose material-boundary test could not finish (its gather overflowed
 * ACB_SEGMENTS): a certificate walk counts it UNMEASURED and goes on, so the
 * audit completes and says so (2026-10-23: PHerc1218 chart 63586, laid over
 * itself, killed its audit with no message); an objective walk fails as before. */
static void ac_unmeasured(AsmContacts *q, size_t count, size_t a, size_t b)
{
    q->stats.unmeasured += count;
    if (q->n_unmeasured_ab < AC_UNMEASURED_SAMPLE) {
        q->unmeasured_ab[2*q->n_unmeasured_ab] = a; q->unmeasured_ab[2*q->n_unmeasured_ab+1] = b; q->n_unmeasured_ab++;
    }
}

/* The serial walk runs the same two halves the workers do, so the two paths
 * cannot drift apart. */
static void ac_pair(AsmContacts *q, size_t a, size_t b)
{
    double area = 0.0,raw_area=0,ga[6],gb[6];
    int hit = ac_pair_area(q,a,b,&area,&raw_area,q->gradient ? ga : NULL,q->gradient ? gb : NULL,ac_work(q));
    if (hit < 0) { if (q->objective_mode) q->failed = 1; else ac_unmeasured(q,1,a,b); return; }
    if (!hit) return;
    if (q->recording) {
        if (q->n_merged == q->cap_merged) {
            size_t cap = q->cap_merged ? 2*q->cap_merged : 1024;
            q->merged = (AcHit *)ac_grow(q->arena, q->merged, q->n_merged*sizeof *q->merged, cap*sizeof *q->merged);
            q->cap_merged = cap;
        }
        AcHit *h = q->merged+q->n_merged++;
        h->lp = h->seq = 0; h->a = a; h->b = b; h->area = area; h->raw_area = raw_area;
        if (q->gradient) { memcpy(h->ga,ga,sizeof h->ga); memcpy(h->gb,gb,sizeof h->gb); }
    }
    if (!q->dry) ac_commit(q,a,b,area,raw_area,ga,gb);
}

void AsmContacts_separate_frames(AsmContacts *q, const int32_t *frames)
{ if (q) q->frames = frames; }

void AsmContacts_active_faces(AsmContacts *q, const uint8_t *moving)
{ if (q) q->active_faces = moving; }

void AsmContacts_memo_counts(const AsmContacts *q, size_t *measured, size_t *replayed)
{
    if (measured) *measured = q ? q->memo_walks : 0;
    if (replayed) *replayed = q ? q->memo_replays : 0;
}

/* One leaf-leaf visit against the budget, in the walk's own order: the serial
 * walk, the threaded enumeration and the touching walk all count here, so they
 * stop at the same pair whatever the thread count.  Returns 1 to stop. */
static int ac_leaf_pair_budget(AsmContacts *q)
{
    q->stats.leaf_pairs++;
    if (q->budget && q->stats.leaf_pairs > q->budget) { q->over_budget = 1; return 1; }
    return 0;
}

/* Leaf-pair culling. Each triangle's box is taken once per leaf pair; a pair
 * whose boxes are apart by more than 256 ulp of the leaf pair's largest
 * coordinate is a pair ac_separated rejects (its own margin is no larger),
 * so it is skipped before the pair test and the walk commits exactly the
 * same hits. A triangle with a non-finite coordinate gets an unbounded box
 * and is always tested. Leaves above AC_CULL_FACES are not culled. */
#define AC_CULL_FACES 32
typedef struct AcCull { int on; double margin; double box[2][AC_CULL_FACES][4]; } AcCull;

static void ac_cull_leaf(const AsmContacts *q, const AcNode *node, double (*box)[4], double *scale)
{
    for (size_t i = node->begin; i < node->end; i++) {
        const int32_t *f = q->faces+3*q->order[i]; double *b = box[i-node->begin], sum = 0;
        for (int d = 0; d < 2; d++) {
            double lo = q->uv[2*(size_t)f[0]+d], hi = lo; sum += lo;
            for (int k = 1; k < 3; k++) {
                double x = q->uv[2*(size_t)f[k]+d]; sum += x;
                if (x < lo) lo = x; else if (x > hi) hi = x;
            }
            b[d] = lo; b[2+d] = hi;
        }
        if (!isfinite(sum)) { b[0] = b[1] = -INFINITY; b[2] = b[3] = INFINITY; continue; }
        for (int k = 0; k < 4; k++) if (fabs(b[k]) > *scale) *scale = fabs(b[k]);
    }
}

static void ac_cull_leaves(const AsmContacts *q, const AcNode *one, const AcNode *two, int self, AcCull *c)
{
    c->on = one->end-one->begin <= AC_CULL_FACES && two->end-two->begin <= AC_CULL_FACES;
    if (!c->on) return;
    double scale = 1;
    ac_cull_leaf(q,one,c->box[0],&scale);
    if (self) memcpy(c->box[1],c->box[0],(one->end-one->begin)*sizeof c->box[0][0]);
    else ac_cull_leaf(q,two,c->box[1],&scale);
    c->margin = 256*DBL_EPSILON*scale;
}

static int ac_cull_apart(const AcCull *c, size_t i, size_t j)
{
    if (!c->on) return 0;
    const double *a = c->box[0][i], *b = c->box[1][j];
    return b[0]-a[2] > c->margin || a[0]-b[2] > c->margin || b[1]-a[3] > c->margin || a[1]-b[3] > c->margin;
}

/* One leaf pair of the serial walk, out of line so the recursive walk's
 * frame does not carry the cull boxes. */
static AC_NOINLINE void ac_walk_leaves(AsmContacts *q, const AcNode *one, const AcNode *two, int self)
{
    AcCull cull; ac_cull_leaves(q,one,two,self,&cull);
    {size_t n1=one->end-one->begin,n2=two->end-two->begin;ac_count_add(&ac_count_leafpairs,1);ac_count_add(&ac_count_culls,(long long)(self?n1*(n1-1)/2:n1*n2));if(self)ac_count_add(&ac_count_selfleaf,1);}
    for (size_t i = one->begin; i < one->end; i++)
        for (size_t j = self ? i+1 : two->begin; j < two->end; j++)
            if (!ac_cull_apart(&cull,i-one->begin,j-two->begin)) ac_pair(q,q->order[i],q->order[j]);
}

static void ac_walk(AsmContacts *q, int32_t a, int32_t b)
{
    if (q->over_budget) return;
    AcNode *one = q->tree+a, *two = q->tree+b;
    if (!one->active && !two->active) return;
    if (acc_skip_nodes(q,a,b)) return;
    if (!ac_boxes(one,two)) return;
    if (one->left < 0 && two->left < 0) {
        if (ac_leaf_pair_budget(q)) return;
        ac_walk_leaves(q,one,two,a == b);
    } else if (a == b) {
        ac_walk(q,one->left,one->left); ac_walk(q,one->left,one->right); ac_walk(q,one->right,one->right);
    } else if (two->left < 0 || (one->left >= 0 && one->end-one->begin >= two->end-two->begin)) {
        ac_walk(q,one->left,b); ac_walk(q,one->right,b);
    } else { ac_walk(q,a,two->left); ac_walk(q,a,two->right); }
}

/* The same descent as ac_walk with the same pruning, recording leaf-leaf node
 * pairs instead of testing them.  The recorded order IS the walk order. */
static int ac_collect(AsmContacts *q, int32_t a, int32_t b)
{
    if (q->over_budget) return -2;
    AcNode *one = q->tree+a, *two = q->tree+b;
    if (!one->active && !two->active) return 0;
    if (acc_skip_nodes(q,a,b)) return 0;
    if (!ac_boxes(one,two)) return 0;
    if (one->left < 0 && two->left < 0) {
        if (ac_leaf_pair_budget(q)) return -2;
        if (q->n_lp == q->cap_lp) {
            size_t cap = q->cap_lp ? 2*q->cap_lp : 4096;
            q->lp = (AcLeafPair *)ac_grow(q->arena, q->lp, q->n_lp*sizeof *q->lp, cap*sizeof *q->lp);
            q->cap_lp = cap;
        }
        q->lp[q->n_lp].a = a; q->lp[q->n_lp].b = b; q->n_lp++;
        return 0;
    }
    int rc;
    if (a == b) {
        if ((rc = ac_collect(q,one->left,one->left)) != 0) return rc;
        if ((rc = ac_collect(q,one->left,one->right)) != 0) return rc;
        return ac_collect(q,one->right,one->right);
    }
    if (two->left < 0 || (one->left >= 0 && one->end-one->begin >= two->end-two->begin)) {
        if ((rc = ac_collect(q,one->left,b)) != 0) return rc;
        return ac_collect(q,one->right,b);
    }
    if ((rc = ac_collect(q,a,two->left)) != 0) return rc;
    return ac_collect(q,a,two->right);
}

/* Touches only its own buffer, so many of these run at once.  The buffer is
 * thread-local heap for the duration of one measure (an arena cannot be
 * shared between workers) and is released as soon as the merge has read it. */
static void ac_leafpair_run(AsmContacts *q, int32_t lp, AcBuf *buf)
{
    const AcNode *one = q->tree+q->lp[lp].a, *two = q->tree+q->lp[lp].b;
    int self = q->lp[lp].a == q->lp[lp].b, seq = 0;
    AcCull cull; ac_cull_leaves(q,one,two,self,&cull);
    {size_t n1=one->end-one->begin,n2=two->end-two->begin;ac_count_add(&ac_count_leafpairs,1);ac_count_add(&ac_count_culls,(long long)(self?n1*(n1-1)/2:n1*n2));if(self)ac_count_add(&ac_count_selfleaf,1);}
    for (size_t i = one->begin; i < one->end; i++)
        for (size_t j = self ? i+1 : two->begin; j < two->end; j++) {
            if (ac_cull_apart(&cull,i-one->begin,j-two->begin)) continue;
            double area = 0.0,raw_area=0,ga[6],gb[6];
            int hit = ac_pair_area(q,q->order[i],q->order[j],&area,&raw_area,q->gradient ? ga : NULL,q->gradient ? gb : NULL,q->band_work+(buf-q->buf));
            if (hit < 0) {
                if (q->objective_mode) buf->failed = 1;
                else { if (!buf->unmeasured) { buf->first_a = q->order[i]; buf->first_b = q->order[j]; } buf->unmeasured++; }
                continue;
            }
            if (!hit) continue;
            if (buf->n == buf->cap) {
                size_t cap = buf->cap ? 2*buf->cap : 1024;
                AcHit *grown = (AcHit *)realloc(buf->hit, cap*sizeof *grown);
                if (!grown) { buf->failed = 1; return; }
                buf->hit = grown; buf->cap = cap;
            }
            AcHit *h = &buf->hit[buf->n++];
            h->lp = lp; h->seq = seq++; h->a = q->order[i]; h->b = q->order[j]; h->area = area;
            h->raw_area=raw_area;
            if (q->gradient) { memcpy(h->ga,ga,sizeof h->ga); memcpy(h->gb,gb,sizeof h->gb); }
        }
}

static int ac_cmp_hit(const void *x, const void *y)
{
    const AcHit *a = (const AcHit *)x, *b = (const AcHit *)y;
    if (a->lp < b->lp) return -1;
    if (a->lp > b->lp) return 1;
    return (a->seq > b->seq) - (a->seq < b->seq);
}

/* The per-thread hit buffers of one walk, merged into the serial walk's
 * order (leaf pair, then sequence within it) and committed in that order. */
static int ac_merge_commit(AsmContacts *q, int nt)
{
    size_t total = 0;
    for (int t = 0; t < nt; t++) {
        total += q->buf[t].n; if (q->buf[t].failed) q->failed = 1;
        if (q->buf[t].unmeasured) ac_unmeasured(q,q->buf[t].unmeasured,q->buf[t].first_a,q->buf[t].first_b);
    }
    if (q->cap_merged < total) {
        q->merged = (AcHit *)ac_grow(q->arena, q->merged, 0, total*sizeof *q->merged);
        q->cap_merged = total;
    }
    size_t at = 0;
    for (int t = 0; t < nt; t++) {
        if (q->buf[t].n) memcpy(q->merged+at, q->buf[t].hit, q->buf[t].n*sizeof *q->merged);
        at += q->buf[t].n;
        free(q->buf[t].hit); q->buf[t].hit = NULL; q->buf[t].n = q->buf[t].cap = 0;
    }
    if (!total) return 0;
    qsort(q->merged, total, sizeof *q->merged, ac_cmp_hit);
    q->n_merged = total;
    if (q->dry) return 0;
    for (size_t i = 0; i < total; i++) ac_commit(q,q->merged[i].a,q->merged[i].b,q->merged[i].area,q->merged[i].raw_area,q->merged[i].ga,q->merged[i].gb);
    return 0;
}

/* Enumerate serially, measure in parallel, then commit in the serial walk's
 * own order: the tallies receive the same values in the same sequence, so the
 * threaded answer is bit-identical to the serial one. */
static int ac_walk_threaded(AsmContacts *q)
{
    int nt = q->threads;
    q->n_lp = 0; q->n_merged = 0;
    if (ac_collect(q,0,0)) return q->over_budget ? -2 : -1;
    if (!q->n_lp) return 0;
    if (nt > (int)q->n_lp) nt = (int)q->n_lp;
    if (nt < 1) nt = 1;
    if (q->n_buf < nt) {
        q->buf = (AcBuf *)ac_grow(q->arena, q->buf, (size_t)q->n_buf*sizeof *q->buf, (size_t)nt*sizeof *q->buf);
        q->n_buf = nt;
    }
    for (int t = 0; t < nt; t++) { q->buf[t].hit = NULL; q->buf[t].n = 0; q->buf[t].cap = 0; q->buf[t].failed = 0; q->buf[t].unmeasured = 0; }
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
    {
        AcBuf *b = &q->buf[omp_get_thread_num()];
        long count = (long)q->n_lp, k;
#pragma omp for schedule(dynamic,16)
        for (k = 0; k < count; k++) ac_leafpair_run(q,(int32_t)k,b);
    }
#else
    for (size_t k = 0; k < q->n_lp; k++) ac_leafpair_run(q,(int32_t)k,&q->buf[0]);
#endif
    return ac_merge_commit(q, nt);
}

/* POOLED WALK.  A walk on one thread inside a parallel region (a repair
 * correction among the jobs of a wave) borrows the cores the other jobs leave
 * idle (ves_pool.h): the same serial enumeration and the same order-restoring
 * merge as ac_walk_threaded, so the committed sequence is the serial walk's
 * bit for bit.  On 2026-10-05 a 240-cube assembly's repair waves kept 5.6 of 16
 * workers busy while their slowest corrections ran single-threaded walks
 * (contact measurement = 55% of the busy samples). */
#define AC_POOL_MIN_LEAF_PAIRS 64
#define AC_POOL_CHUNK 16

static void ac_pool_leafpairs(void *ctx, size_t begin, size_t end, int slot)
{
    AsmContacts *q = (AsmContacts *)ctx;
    for (size_t k = begin; k < end; k++) ac_leafpair_run(q, (int32_t)k, &q->buf[slot]);
}

static int ac_walk_pooled(AsmContacts *q, int helpers)
{
    q->n_lp = 0; q->n_merged = 0;
    int rc = ac_collect(q,0,0);
    if (rc == -2 || (!rc && q->n_lp < AC_POOL_MIN_LEAF_PAIRS)) {
        /* Over its budget (the serial walk commits up to the stop) or too
         * small to share: walk serially, exactly as before. */
        q->over_budget = 0; q->stats.leaf_pairs = 0; q->n_lp = 0;
        ac_walk(q,0,0);
        return 0;
    }
    if (rc) return -1;
    int nt = 1 + helpers;
    if (q->n_buf < nt) {
        q->buf = (AcBuf *)ac_grow(q->arena, q->buf, (size_t)q->n_buf*sizeof *q->buf, (size_t)nt*sizeof *q->buf);
        q->n_buf = nt;
    }
    if (q->n_band_work < nt) {
        q->band_work = ARENA_ALLOC(q->arena, (size_t)nt*sizeof *q->band_work);
        q->n_band_work = nt;
    }
    for (int t = 0; t < nt; t++) { q->buf[t].hit = NULL; q->buf[t].n = 0; q->buf[t].cap = 0; q->buf[t].failed = 0; q->buf[t].unmeasured = 0; }
    VesPool_run(q->n_lp, AC_POOL_CHUNK, ac_pool_leafpairs, q, helpers);
    return ac_merge_commit(q, nt);
}

void AsmContacts_threads(AsmContacts *q,int threads)
{
    if(!q)return;q->threads=threads;int count=threads>1?threads:1;
    if(count>q->n_band_work){q->band_work=ARENA_ALLOC(q->arena,(size_t)count*sizeof *q->band_work);q->n_band_work=count;}
}

void AsmContacts_no_memo(AsmContacts *q) { if (q) q->memo_off = 1; }

/* A refit that costs the WHOLE field on every objective evaluation is what
 * made a one-chart admission as expensive as a whole-sheet repair: the
 * validation scan alone reads 2*nv coordinates, and every leaf re-reads its
 * faces, so the repair's per-attempt cost never fell with the size of the set
 * it was allowed to move.  Geometry outside the caller's active mask cannot
 * have moved since the last refit, so its boxes still stand and only its
 * active flag needs recomputing.  With no mask every node is active and this
 * is the original full refit. */
static int ac_refit(AsmContacts *q,const double *uv)
{
    if(!q || !uv)return 0;
    q->probe_exclude=NULL;
    if(!q->active_faces){for(size_t i=0;i<2*q->nv;i++)if(!isfinite(uv[i]))return 0;}
    q->uv=uv;
    int finite = 1;
    for (size_t i = q->nodes; i-- > 0;) {
        AcNode *node = q->tree+i;
        if (node->left < 0) {
            if (ac_leaf_active(q,node)) { if (!ac_leaf_box(q,node)) finite = 0; }
            else node->active = 0;
        } else {
            node->active = q->tree[node->left].active || q->tree[node->right].active;
            if (node->active) for (int d = 0; d < 2; d++) {
            node->lo[d] = ac_min(q->tree[node->left].lo[d],q->tree[node->right].lo[d]);
            node->hi[d] = ac_max(q->tree[node->left].hi[d],q->tree[node->right].hi[d]);
            }
        }
    }
    return acm_refit(q->material,uv,q->active_faces) && finite;
}
/* The memo slot of this request's kind, or NULL.  Only masked walks are
 * kept: those are the solver's, whose coordinates repeat; a complete walk is
 * an audit or a one-off starting measurement. */
static AcMemo *ac_memo_slot(AsmContacts *q,int gradient)
{
    if(q->memo_off || !q->active_faces)return NULL;
    if(!q->memo)q->memo=ARENA_CALLOC(q->arena,AC_MEMO_KINDS,sizeof *q->memo);
    return q->memo+(q->objective_mode?2:0)+(gradient?1:0);
}

static int ac_memo_hit(const AsmContacts *q,const AcMemo *m,const double *uv)
{
    if(!m->valid || m->version!=q->version || m->frames!=q->frames)return 0;
    if(q->budget && m->leaf_pairs>q->budget)return 0;   /* the fresh walk stops at its budget */
    return !memcmp(m->mask,q->active_faces,q->nf) && !memcmp(m->uv,uv,2*q->nv*sizeof *uv);
}

static void ac_memo_store(AsmContacts *q,AcMemo *m,const double *uv)
{
    if(!m->uv){m->uv=ARENA_ALLOC(q->arena,2*q->nv*sizeof *m->uv);m->mask=ARENA_ALLOC(q->arena,q->nf);}
    if(m->cap<q->n_merged){
        size_t cap=m->cap?m->cap:1024;while(cap<q->n_merged)cap*=2;
        m->hit=(AcHit *)ac_grow(q->arena,NULL,0,cap*sizeof *m->hit);m->cap=cap;
    }
    if(q->n_merged)memcpy(m->hit,q->merged,q->n_merged*sizeof *m->hit);
    memcpy(m->uv,uv,2*q->nv*sizeof *uv);memcpy(m->mask,q->active_faces,q->nf);
    m->n=q->n_merged;m->leaf_pairs=q->stats.leaf_pairs;m->version=q->version;m->frames=q->frames;m->valid=1;
}

/* Diagnostic (ASM_CONTACT_MEMO_VERIFY): re-measure a replay's geometry
 * without committing anything and compare the hit list with the memo, bit
 * for bit and in order.  Tallies and flags are left as the caller set them. */
static int ac_memo_verified(AsmContacts *q,const AcMemo *m,int gradient)
{
    AsmContactStats saved=q->stats;int failed=q->failed,over=q->over_budget,rc=0;
    q->stats.leaf_pairs=0;q->dry=1;q->recording=1;q->n_merged=0;q->failed=0;q->over_budget=0;
    if(q->threads>1)rc=ac_walk_threaded(q);else ac_walk(q,0,0);
    q->dry=0;q->recording=0;
    int same=!rc && !q->failed && !q->over_budget && q->n_merged==m->n && q->stats.leaf_pairs==m->leaf_pairs;
    for(size_t i=0;same && i<m->n;i++){
        const AcHit *x=q->merged+i,*y=m->hit+i;
        same=x->a==y->a && x->b==y->b && !memcmp(&x->area,&y->area,sizeof x->area) &&
             !memcmp(&x->raw_area,&y->raw_area,sizeof x->raw_area) &&
             (!gradient || (!memcmp(x->ga,y->ga,sizeof x->ga) && !memcmp(x->gb,y->gb,sizeof x->gb)));
    }
    q->stats=saved;q->failed=failed;q->over_budget=over;
    return same;
}

static int ac_measure(AsmContacts *q,const double *uv,double *gradient,
                      AsmContactVisit visit,void *context,AsmContactStats *out,double *raw_area,
                      AsmContactWeight weight,void *weight_context)
{
    if(out)memset(out,0,sizeof *out);
    if(raw_area)*raw_area=NAN;
    if(!out || !ac_refit(q,uv))return -1;
    acc_certify_all(q);
    q->gradient=gradient;q->visit=visit;q->context=context;
    q->objective_weight=weight;q->objective_context=weight_context;
    q->objective_mode=raw_area!=NULL;q->objective_area=0;
    memset(&q->stats,0,sizeof q->stats);q->failed=0;q->over_budget=0;q->n_unmeasured_ab=0;
    if(gradient)memset(gradient,0,2*q->nv*sizeof(double));
    AcMemo *memo=ac_memo_slot(q,gradient!=NULL);int fresh=1,replay=memo && ac_memo_hit(q,memo,uv);
    if(replay && q->memo_verify && !ac_memo_verified(q,memo,gradient!=NULL)){
        q->memo_mismatches++;memo->valid=0;replay=0;
        fprintf(stderr,"[contacts] memo replay differs from its geometry (%zu hits, %zu leaf pairs); measured afresh\n",memo->n,memo->leaf_pairs);
    }
    if(replay){
        /* Same leaf pairs visited, same hits committed in the same order. */
        q->stats.leaf_pairs=memo->leaf_pairs;q->memo_replays++;fresh=0;
        for(size_t i=0;i<memo->n;i++){const AcHit *h=memo->hit+i;ac_commit(q,h->a,h->b,h->area,h->raw_area,h->ga,h->gb);}
    }else{
        q->n_merged=0;q->recording=memo!=NULL;
        int helpers = q->threads > 1 ? 0 : VesPool_share();
        if (q->threads > 1 || helpers > 0) {
            int rc = q->threads > 1 ? ac_walk_threaded(q) : ac_walk_pooled(q,helpers);
            if (rc && rc != -2) {
                q->recording=0;q->objective_mode=0;q->objective_weight=NULL;q->objective_context=NULL;return -1;
            }
        }
        else ac_walk(q,0,0);
        q->recording=0;if(memo)q->memo_walks++;
    }
    q->objective_mode=0;q->objective_weight=NULL;q->objective_context=NULL;
    *out = q->stats; out->complete = !q->over_budget && !q->failed && !q->stats.unmeasured;
    if (q->over_budget) return -2;
    if (gradient) for (size_t i = 0; i < 2*q->nv; i++) if (!isfinite(gradient[i])) return -1;
    if(q->failed || !isfinite(out->area) || !isfinite(q->objective_area))return -1;
    if(q->stats.unmeasured)
        fprintf(stderr,"[contacts] %zu pair(s) UNMEASURED: the material-boundary gather overflowed its %d segments (first pair: faces %zu, %zu); the walk went on, its overlaps are a lower bound\n",
                q->stats.unmeasured,(int)ACB_SEGMENTS,q->unmeasured_ab[0],q->unmeasured_ab[1]);
    if(memo && fresh && !q->stats.unmeasured)ac_memo_store(q,memo,uv);
    if(raw_area)*raw_area=q->objective_area;return 0;
}

int AsmContacts_measure(AsmContacts *q,const double *uv,double *gradient,AsmContactVisit visit,void *context,AsmContactStats *out)
{return ac_measure(q,uv,gradient,visit,context,out,NULL,NULL,NULL);}
size_t AsmContacts_unmeasured_sample(const AsmContacts *q,size_t *ab,size_t cap)
{
    size_t n = q ? q->n_unmeasured_ab : 0; if (n > cap) n = cap;
    if (ab && n) memcpy(ab,q->unmeasured_ab,2*n*sizeof *ab);
    return n;
}
int AsmContacts_objective(AsmContacts *q,const double *uv,double *gradient,AsmContactVisit visit,void *context,AsmContactStats *out,double *raw_area)
{if(!raw_area)return -1;return ac_measure(q,uv,gradient,visit,context,out,raw_area,NULL,NULL);}
int AsmContacts_weighted_objective(AsmContacts *q,const double *uv,double *gradient,
    AsmContactVisit visit,void *context,AsmContactWeight weight,void *weight_context,
    AsmContactStats *out,double *raw_area)
{if(!raw_area)return -1;return ac_measure(q,uv,gradient,visit,context,out,raw_area,weight,weight_context);}

int AsmContacts_refit(AsmContacts *q, const double *uv)
{
    if (!q || !uv) return -1;
    const uint8_t *save = q->active_faces;
    q->active_faces = NULL;
    int fine = ac_refit(q,uv);
    q->active_faces = save;
    return fine ? 0 : -1;
}

int AsmContacts_refit_changed(AsmContacts *q,const double *uv,const uint8_t *changed)
{
    if(!q || !uv)return -1;
    const uint8_t *save=q->active_faces;q->active_faces=changed;
    int fine=ac_refit(q,uv);q->active_faces=save;
    return fine?0:-1;
}

/* One external triangle against the subtree at `node`.  `box` is its (lo, hi). */
static void ac_probe_node(AsmContacts *q, int32_t node, const double box[4], const double tri[6],
                          size_t ext, const uint8_t *exclude, AsmContactVisit visit, void *context,
                          const AsmContacts *source,const double *source_uv,size_t source_face)
{
    if (q->over_budget) return;
    if(exclude && exclude==q->probe_exclude && !q->probe_nodes[node])return;
    const AcNode *nd = q->tree+node;
    for (int d = 0; d < 2; d++) if (nd->lo[d] > box[2+d] || box[d] > nd->hi[d]) return;
    if (nd->left < 0) {
        if (ac_leaf_pair_budget(q)) return;
        for (size_t i = nd->begin; i < nd->end; i++) {
            size_t face = q->order[i];
            if (exclude && exclude[face]) continue;
            const int32_t *fb = q->faces+3*face;
            double y[6];
            for (int k = 0; k < 3; k++) for (int d = 0; d < 2; d++) y[2*k+d] = q->uv[2*(size_t)fb[k]+d];
            double area=ac_material_pair(ac_work(q),source?source->material:NULL,source_face,source_uv,tri,
                q->material,face,q->uv,y,q->tolerance,NULL,NULL,NULL);
            if (!isfinite(area)) { q->failed = 1; return; }
            if (!(area > 0)) continue;
            q->stats.pairs++; q->stats.area += area;
            if (visit) visit(context,ext,face,area);
        }
        return;
    }
    ac_probe_node(q,nd->left,box,tri,ext,exclude,visit,context,source,source_uv,source_face);
    ac_probe_node(q,nd->right,box,tri,ext,exclude,visit,context,source,source_uv,source_face);
}

static int ac_probe(AsmContacts *q, const double *uv, size_t n_ext, const double *tri,
                     const uint8_t *exclude_face, AsmContactVisit visit, void *context, AsmContactStats *out,int refit,
                     AsmContacts *source,const double *source_uv,const size_t *source_faces)
{
    if(out)memset(out,0,sizeof *out);
    if(!q || !uv || !out || (n_ext && !tri && !source) ||
       (source && (!source_uv || (source==q && source_uv!=uv))))return -1;
    const uint8_t *save = q->active_faces;
    q->active_faces = NULL;
    int fine = refit ? ac_refit(q,uv) : 1;
    if(source && !acm_refit(source->material,source_uv,NULL))fine=0;
    memset(&q->stats,0,sizeof q->stats); q->failed = 0; q->over_budget = 0; q->n_unmeasured_ab = 0;
    for (size_t e = 0; fine && e < n_ext && !q->over_budget && !q->failed; e++) {
        double xy[6];const double *t;size_t source_face=source_faces?source_faces[e]:e;
        if(source){
            if(source_face>=source->nf){q->failed=1;break;}
            ac_face_xy(source,source_uv,source_face,xy);t=xy;
        }else t=tri+6*e;
        int finite = 1;
        for (int k = 0; k < 6; k++) if (!isfinite(t[k])) finite = 0;
        if (!finite) { q->failed = 1; break; }
        double box[4] = { ac_min(t[0],ac_min(t[2],t[4])), ac_min(t[1],ac_min(t[3],t[5])),
                          ac_max(t[0],ac_max(t[2],t[4])), ac_max(t[1],ac_max(t[3],t[5])) };
        ac_probe_node(q,0,box,t,e,exclude_face,visit,context,source,source_uv,source_face);
    }
    *out = q->stats; out->complete = fine && !q->over_budget && !q->failed;
    q->active_faces = save;
    if (!fine) return -1;
    if (q->over_budget) return -2;
    return q->failed || !isfinite(out->area) ? -1 : 0;
}

int AsmContacts_probe(AsmContacts *q,const double *uv,size_t n_ext,const double *tri,
                      const uint8_t *exclude_face,AsmContactVisit visit,void *context,AsmContactStats *out)
{return ac_probe(q,uv,n_ext,tri,exclude_face,visit,context,out,1,NULL,NULL,NULL);}
int AsmContacts_probe_prepared(AsmContacts *q,size_t n_ext,const double *tri,
                               const uint8_t *exclude_face,AsmContactVisit visit,void *context,AsmContactStats *out)
{return ac_probe(q,q?q->uv:NULL,n_ext,tri,exclude_face,visit,context,out,0,NULL,NULL,NULL);}
int AsmContacts_probe_mesh(AsmContacts *q,const double *uv,AsmContacts *source,const double *source_uv,
    size_t n_ext,const size_t *source_faces,const uint8_t *exclude,AsmContactVisit visit,void *context,AsmContactStats *out)
{
    if(!source){if(out)memset(out,0,sizeof *out);return -1;}
    return ac_probe(q,uv,n_ext,NULL,exclude,visit,context,out,1,source,source_uv,source_faces);
}
int AsmContacts_probe_mesh_prepared(AsmContacts *q,AsmContacts *source,const double *source_uv,
    size_t n_ext,const size_t *source_faces,const uint8_t *exclude,AsmContactVisit visit,void *context,AsmContactStats *out)
{
    if(!source){if(out)memset(out,0,sizeof *out);return -1;}
    return ac_probe(q,q?q->uv:NULL,n_ext,NULL,exclude,visit,context,out,0,source,source_uv,source_faces);
}
void AsmContacts_prepare_exclusion(AsmContacts *q,const uint8_t *exclude)
{
    if(!q)return;q->probe_exclude=NULL;if(!exclude)return;
    if(q->nodes>q->probe_capacity){q->probe_nodes=ARENA_ALLOC(q->arena,q->nodes);q->probe_capacity=q->nodes;}
    for(size_t i=q->nodes;i-->0;){
        const AcNode *node=q->tree+i;int active=0;
        if(node->left<0){for(size_t j=node->begin;j<node->end && !active;j++)active=!exclude[q->order[j]];}
        else active=q->probe_nodes[node->left] || q->probe_nodes[node->right];
        q->probe_nodes[i]=(uint8_t)active;
    }
    q->probe_exclude=exclude;
}

#include "asm_contacts_touch.inc"

/* A permutation of complete, overlapping triangle pairs preserves every
 * geometric answer but destroys the old spatial partition. Rebuilding must
 * recover compact leaves without changing coordinates, masks or gradients. */
static int ac_repartition_control(void)
{
    enum { CELLS=256,NF=2*CELLS,NV=3*NF };
    Arena_T arena=Arena_new();double uv[2*NV],saved[2*NV],before_gradient[2*NV],after_gradient[2*NV];
    int32_t faces[3*NF],frames[NF];uint8_t active[NF];
    const double triangle[]={0,0,1,0,0,1};AsmContacts *q=NULL;int fails=0;size_t expected=0;
    for(int phase=0;phase<2;phase++){
        for(int i=0;i<NF;i++){
            int cell=i/2,where=cell;
            if(phase){where=0;for(int bit=0;bit<8;bit++)where|=((cell>>bit)&1)<<(7-bit);}
            for(int k=0;k<3;k++){
                int v=3*i+k;faces[v]=v;
                uv[2*v]=4*(where%16)+triangle[2*k]+.25*(i%2);
                uv[2*v+1]=4*(where/16)+triangle[2*k+1]+.125*(i%2);
            }
            frames[i]=cell%4+(i%10==1?100:0);active[i]=(uint8_t)(i%3==0);
        }
        if(!phase){q=AsmContacts_new(arena,NV,faces,NF,uv,1e-7);if(!q){Arena_dispose(&arena);return 1;}}
    }
    for(int cell=0;cell<CELLS;cell++)expected+=frames[2*cell]==frames[2*cell+1] && (active[2*cell] || active[2*cell+1]);
    memcpy(saved,uv,sizeof uv);AsmContacts_active_faces(q,active);AsmContacts_separate_frames(q,frames);
    AsmContactStats before={0},after={0},touch={0};double loose=0,compact=0;
    if(AsmContacts_measure(q,uv,before_gradient,NULL,NULL,&before))fails++;
    for(size_t n=0;n<q->nodes;n++)if(q->tree[n].left<0)
        loose+=(q->tree[n].hi[0]-q->tree[n].lo[0])*(q->tree[n].hi[1]-q->tree[n].lo[1]);
    if(AsmContacts_rebuild(q,uv) || AsmContacts_measure(q,uv,after_gradient,NULL,NULL,&after))fails++;
    for(size_t n=0;n<q->nodes;n++)if(q->tree[n].left<0)
        compact+=(q->tree[n].hi[0]-q->tree[n].lo[0])*(q->tree[n].hi[1]-q->tree[n].lo[1]);
    if(!(compact<.25*loose) || before.pairs!=expected || after.pairs!=expected ||
       fabs(before.area-.1953125*expected)>1e-10 || fabs(before.area-after.area)>1e-10 || memcmp(saved,uv,sizeof uv))fails++;
    for(int k=0;k<2*NV;k++)if(fabs(before_gradient[k]-after_gradient[k])>1e-10)fails++;
    for(int k=0;k<3*NF;k++)if(faces[k]!=k)fails++;
    /* Rebuilding on UNCHANGED geometry must reproduce the same partition and
     * the same answers: the existing tree already is that rebuild, and the
     * repair asks for one per admission attempt. */
    {
        size_t nodes_before=q->nodes;double again=0;AsmContactStats repeat={0};
        if(AsmContacts_rebuild(q,uv) || AsmContacts_measure(q,uv,NULL,NULL,NULL,&repeat))fails++;
        for(size_t n=0;n<q->nodes;n++)if(q->tree[n].left<0)
            again+=(q->tree[n].hi[0]-q->tree[n].lo[0])*(q->tree[n].hi[1]-q->tree[n].lo[1]);
        if(q->nodes!=nodes_before || fabs(again-compact)>1e-10 ||
           repeat.pairs!=after.pairs || fabs(repeat.area-after.area)>1e-10)fails++;
    }
    if(AsmContacts_touching(q,uv,NULL,NULL,&touch) || touch.pairs!=expected || fabs(touch.area-after.area)>1e-10)fails++;
    AsmContacts_active_faces(q,NULL);AsmContacts_separate_frames(q,NULL);
    if(AsmContacts_measure(q,uv,NULL,NULL,NULL,&after) || after.pairs!=CELLS || fabs(after.area-.1953125*CELLS)>1e-10)fails++;
    Arena_dispose(&arena);
    fprintf(stderr,"  contact repartition after chart permutation: %s (%d failures; leaf area %.9g -> %.9g)\n",fails?"FAIL":"ok",fails,loose,compact);
    return fails;
}

/* The threaded walk must be BIT-identical to the serial one, not merely close:
 * the repair accepts or rejects trials on these numbers, so a reduction that
 * depended on the thread count would make the pipeline irreproducible.  The
 * exact comparisons below are the contract, not a tolerance oversight. */
/* Diagnostic micro-benchmark (sheet_assemble --contact-bench): objective-mode
 * walks with gradients over a synthetic two-layer field of `cells` jittered
 * cells, single-threaded with the memo off.  Prints seconds per walk and a
 * checksum of the area, raw area and gradient so builds can be compared for
 * bit-identical answers as well as speed. */
int AsmContacts_bench(int cells, int reps)
{
    if (cells < 16) cells = 16;
    if (reps < 1) reps = 1;
    /* A connected jittered grid (shared vertices, as a chart is) with every
     * 61st cell folded over its neighbour: vertex- and edge-sharing
     * neighbours dominate the pair tests, a few real overlaps exercise the
     * material band at the production 4-vox width. */
    int side = 1; while (side*side < cells) side++;
    size_t nv = (size_t)(side+1)*(size_t)(side+1), nf = 2*(size_t)side*(size_t)side;
    Arena_T arena = Arena_new();
    double *uv = ARENA_ALLOC(arena,2*nv*sizeof(double)), *g = ARENA_ALLOC(arena,2*nv*sizeof(double));
    int32_t *faces = ARENA_ALLOC(arena,3*nf*sizeof(int32_t));
    uint8_t *active = ARENA_ALLOC(arena,nf);
    uint64_t seed = 0x9E3779B97F4A7C15ull;
    for (size_t y = 0; y <= (size_t)side; y++) for (size_t x = 0; x <= (size_t)side; x++) {
        seed = seed*6364136223846793005ull+1442695040888963407ull;
        double jx = (double)(seed>>40)/16777216.0-.5, jy = (double)((seed>>16)&0xFFFFFF)/16777216.0-.5;
        size_t v = y*(size_t)(side+1)+x;
        uv[2*v] = 1.5*(double)x+.5*jx; uv[2*v+1] = 1.5*(double)y+.5*jy;
        if ((v % 61) == 7) { uv[2*v] += 2.2; uv[2*v+1] += 1.1; }   /* a fold */
    }
    for (size_t y = 0; y < (size_t)side; y++) for (size_t x = 0; x < (size_t)side; x++) {
        size_t c = y*(size_t)side+x, a = y*(size_t)(side+1)+x, b = a+1, d = a+(size_t)(side+1), e = d+1;
        int32_t *f = faces+6*c;
        f[0]=(int32_t)a; f[1]=(int32_t)b; f[2]=(int32_t)e; f[3]=(int32_t)a; f[4]=(int32_t)e; f[5]=(int32_t)d;
        active[2*c] = active[2*c+1] = (uint8_t)(c%5 != 4);
    }
    AsmContacts *q = AsmContacts_new(arena,nv,faces,nf,uv,4.0);
    if (!q) { Arena_dispose(&arena); return 1; }
    q->memo_off = 1; AsmContacts_threads(q,1); AsmContacts_active_faces(q,active);
    int saved = VesPool_cores(); VesPool_configure(1);
    double area = 0, raw = 0, best = 1e300; uint64_t h = 1469598103934665603ull; size_t pairs = 0, lp = 0;
    for (int r = 0; r < reps; r++) {
        AsmContactStats st; double ra = 0;
        double t0 = ves_clock_sec();
        if (AsmContacts_objective(q,uv,g,NULL,NULL,&st,&ra)) { fprintf(stderr,"contact bench: measure failed\n"); break; }
        double dt = ves_clock_sec()-t0; if (dt < best) best = dt;
        area = st.area; raw = ra; pairs = st.pairs; lp = st.leaf_pairs;
    }
    double best_value = 1e300;
    for (int r = 0; r < reps; r++) {
        AsmContactStats st; double ra = 0;
        double t0 = ves_clock_sec();
        if (AsmContacts_objective(q,uv,NULL,NULL,NULL,&st,&ra)) break;
        double dt = ves_clock_sec()-t0; if (dt < best_value) best_value = dt;
    }
    fprintf(stderr,"contact bench: value-only objective best %.6f s/walk\n",best_value);
    for (size_t k = 0; k < 2*nv; k++) { uint64_t b; memcpy(&b,&g[k],sizeof b); h = (h^b)*1099511628211ull; }
    fprintf(stderr,"contact bench: %d cells, %zu faces, %zu leaf pairs, %zu pairs, best %.6f s/walk; area %.17g raw %.17g gradient hash %016llx\n",
            cells,nf,lp,pairs,best,area,raw,(unsigned long long)h);
    VesPool_configure(saved);
    Arena_dispose(&arena);
    return 0;
}

/* A pooled walk (one thread inside an outer job, helpers lent by ves_pool)
 * commits the serial walk's sequence: same pairs, area and gradient, bit for
 * bit, masked and complete; the memo is off so every request really walks. */
static int ac_pooled_control(void)
{
    enum { CELLS=1024, NF=2*CELLS, NV=3*NF };
    Arena_T arena=Arena_new();
    double *uv=ARENA_ALLOC(arena,2*NV*sizeof(double));
    double *serial=ARENA_ALLOC(arena,2*NV*sizeof(double)),*pooled=ARENA_ALLOC(arena,2*NV*sizeof(double));
    int32_t *faces=ARENA_ALLOC(arena,3*NF*sizeof(int32_t));
    uint8_t *active=ARENA_ALLOC(arena,NF);
    const double triangle[]={0,0,1,0,0,1};
    for(int i=0;i<NF;i++){
        int cell=i/2;
        for(int k=0;k<3;k++){
            int v=3*i+k;faces[v]=v;
            uv[2*v]=4*(cell%32)+triangle[2*k]+.25*(i%2);
            uv[2*v+1]=4*(cell/32)+triangle[2*k+1]+.125*(i%2);
        }
        active[i]=(uint8_t)(cell%3!=2);
    }
    AsmContacts *q=AsmContacts_new(arena,NV,faces,NF,uv,1e-7);
    if(!q){Arena_dispose(&arena);return 1;}
    q->memo_off=1;
    int fails=0,saved=VesPool_cores();size_t runs0=0,runs1=0,lp=0;
    VesPool_counts(&runs0,NULL);
    for(int masked=1;masked>=0;masked--){
        AsmContactStats one={0},many={0};
        AsmContacts_active_faces(q,masked?active:NULL);AsmContacts_threads(q,1);
        VesPool_configure(1);
        if(AsmContacts_measure(q,uv,serial,NULL,NULL,&one))fails++;
        VesPool_configure(8);VesPool_job_enter();
        if(AsmContacts_measure(q,uv,pooled,NULL,NULL,&many))fails++;
        VesPool_job_leave();
        lp=one.leaf_pairs;
        if(one.pairs!=many.pairs || one.area!=many.area || one.leaf_pairs!=many.leaf_pairs || !(one.pairs>0))fails++;
        for(int k=0;k<2*NV;k++)if(serial[k]!=pooled[k]){fails++;break;}
    }
    VesPool_counts(&runs1,NULL);
    if(runs1<runs0+2)fails++;              /* both measurements actually borrowed helpers */
    VesPool_configure(saved);AsmContacts_active_faces(q,NULL);
    fprintf(stderr,"  contact walk pooled vs serial: %s (%d failures; %zu leaf pairs, %zu pooled runs)\n",
            fails?"FAIL":"ok",fails,lp,runs1-runs0);
    Arena_dispose(&arena);
    return fails;
}

static int ac_threaded_control(void)
{
    enum { CELLS=512, NF=2*CELLS, NV=3*NF };
    Arena_T arena=Arena_new();
    double *uv=ARENA_ALLOC(arena,2*NV*sizeof(double));
    double *serial=ARENA_ALLOC(arena,2*NV*sizeof(double)),*threaded=ARENA_ALLOC(arena,2*NV*sizeof(double));
    int32_t *faces=ARENA_ALLOC(arena,3*NF*sizeof(int32_t));
    uint8_t *active=ARENA_ALLOC(arena,NF);
    const double triangle[]={0,0,1,0,0,1};
    for(int i=0;i<NF;i++){
        int cell=i/2;
        for(int k=0;k<3;k++){
            int v=3*i+k;faces[v]=v;
            uv[2*v]=4*(cell%32)+triangle[2*k]+.25*(i%2);
            uv[2*v+1]=4*(cell/32)+triangle[2*k+1]+.125*(i%2);
        }
        active[i]=(uint8_t)(cell%3!=2);   /* whole cells, so a masked pair really is skipped */
    }
    AsmContacts *q=AsmContacts_new(arena,NV,faces,NF,uv,1e-7);
    if(!q){Arena_dispose(&arena);return 1;}
    int fails=0;AsmContactStats masked_one={0},masked_many={0},all_one={0},all_many={0};
    for(int masked=1;masked>=0;masked--){
        AsmContactStats *one=masked?&masked_one:&all_one,*many=masked?&masked_many:&all_many;
        AsmContacts_active_faces(q,masked?active:NULL);
        AsmContacts_threads(q,1);
        if(AsmContacts_measure(q,uv,serial,NULL,NULL,one))fails++;
        AsmContacts_threads(q,8);
        if(AsmContacts_measure(q,uv,threaded,NULL,NULL,many))fails++;
        if(one->pairs!=many->pairs || one->area!=many->area)fails++;
        int gradient_differs=0;
        for(int k=0;k<2*NV;k++)if(serial[k]!=threaded[k])gradient_differs=1;
        fails+=gradient_differs;
    }
    if(!(all_one.pairs>masked_one.pairs) || !(masked_one.pairs>0))fails++;
    AsmContacts_threads(q,1);
    fprintf(stderr,"  contact walk threaded vs serial: %s (%d failures; %zu pairs masked, %zu complete)\n",
            fails?"FAIL":"ok",fails,masked_one.pairs,all_one.pairs);
    Arena_dispose(&arena);
    return fails;
}

/* A replayed walk must BE the fresh walk: identical tallies, objective area,
 * gradient bits and callback sequences, for every kind and thread count; and
 * a coordinate or mask edited in place, a new partition, or a budget the walk
 * would exceed must all measure afresh. */
typedef struct AcMemoTrace { uint64_t hash; size_t calls; } AcMemoTrace;
typedef struct AcMemoResult { int rc; AsmContactStats stats; double raw; AcMemoTrace visit, weight; double *gradient; } AcMemoResult;

static uint64_t ac_memo_mix(uint64_t h, uint64_t v)
{ return h ^ (v+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2)); }

static void ac_memo_trace(AcMemoTrace *t, size_t a, size_t b, double value)
{
    uint64_t bits; memcpy(&bits,&value,sizeof bits);
    t->hash = ac_memo_mix(ac_memo_mix(ac_memo_mix(t->hash,(uint64_t)a),(uint64_t)b),bits); t->calls++;
}

static void ac_memo_visit(void *context, size_t a, size_t b, double area)
{ ac_memo_trace(context,a,b,area); }

static double ac_memo_weight(void *context, size_t a, size_t b, double raw_area)
{ ac_memo_trace(context,a,b,raw_area); return (a+b)%3 ? 1.0 : 0.5; }

static void ac_memo_run(AsmContacts *q, const double *uv, int kind, AcMemoResult *r)
{
    int objective = kind >= 2, gradient = kind & 1;
    memset(&r->visit,0,sizeof r->visit); memset(&r->weight,0,sizeof r->weight); r->raw = 0;
    memset(r->gradient,0,2*q->nv*sizeof *r->gradient);
    r->rc = objective ?
        AsmContacts_weighted_objective(q,uv,gradient ? r->gradient : NULL,ac_memo_visit,&r->visit,ac_memo_weight,&r->weight,&r->stats,&r->raw) :
        AsmContacts_measure(q,uv,gradient ? r->gradient : NULL,ac_memo_visit,&r->visit,&r->stats);
}

static int ac_memo_same(const AcMemoResult *a, const AcMemoResult *b, size_t n)
{
    return a->rc == b->rc && a->stats.pairs == b->stats.pairs && !memcmp(&a->stats.area,&b->stats.area,sizeof a->stats.area) &&
           a->stats.leaf_pairs == b->stats.leaf_pairs && a->stats.complete == b->stats.complete &&
           !memcmp(&a->raw,&b->raw,sizeof a->raw) && a->visit.hash == b->visit.hash && a->visit.calls == b->visit.calls &&
           a->weight.hash == b->weight.hash && a->weight.calls == b->weight.calls && !memcmp(a->gradient,b->gradient,n*sizeof *a->gradient);
}

/* One request measured afresh (memo off) and through the memo; `replay`
 * says whether the memo request must have been served from its slot. */
static int ac_memo_check(AsmContacts *q, const double *uv, int kind, AcMemoResult *fresh, AcMemoResult *memo, int replay)
{
    q->memo_off = 1; ac_memo_run(q,uv,kind,fresh); q->memo_off = 0;
    size_t before = q->memo_replays; ac_memo_run(q,uv,kind,memo);
    return (q->memo_replays != before+(size_t)replay) + !ac_memo_same(fresh,memo,2*q->nv);
}

static int ac_memo_control(void)
{
    enum { CELLS=256, NF=2*CELLS, NV=3*NF };
    Arena_T arena = Arena_new(); int fails = 0; size_t replays = 0, hits = 0;
    double *base = ARENA_ALLOC(arena,2*NV*sizeof(double)), *uv = ARENA_ALLOC(arena,2*NV*sizeof(double));
    double *moved = ARENA_ALLOC(arena,2*NV*sizeof(double));
    int32_t *faces = ARENA_ALLOC(arena,3*NF*sizeof(int32_t)); uint8_t *active = ARENA_ALLOC(arena,NF);
    AcMemoResult fresh = {0}, memo = {0};
    fresh.gradient = ARENA_ALLOC(arena,2*NV*sizeof(double)); memo.gradient = ARENA_ALLOC(arena,2*NV*sizeof(double));
    const double triangle[] = {0,0,1,0,0,1};
    for (int i = 0; i < NF; i++) {
        int cell = i/2;
        for (int k = 0; k < 3; k++) {
            int v = 3*i+k; faces[v] = v;
            base[2*v] = 4*(cell%16)+triangle[2*k]+.25*(i%2)+.01*(cell%5);
            base[2*v+1] = 4*(cell/16)+triangle[2*k+1]+.125*(i%2);
        }
        active[i] = (uint8_t)(cell%3 != 2);
    }
    for (int threads = 1; threads <= 8; threads *= 8) {
        memcpy(uv,base,2*NV*sizeof *uv);
        AsmContacts *q = AsmContacts_new_threads(arena,NV,faces,NF,uv,1e-7,threads);
        if (!q) { fails++; continue; }
        AsmContacts_active_faces(q,active); q->memo_verify = 1;
        for (int kind = 0; kind < AC_MEMO_KINDS; kind++) {
            fails += ac_memo_check(q,uv,kind,&fresh,&memo,0);   /* stores */
            fails += ac_memo_check(q,uv,kind,&fresh,&memo,1);   /* replays */
            hits += fresh.stats.pairs+fresh.weight.calls;
            uv[2*(3*2)] += 1e-3;                                  /* face 2 (cell 1) is active */
            fails += ac_memo_check(q,uv,kind,&fresh,&memo,0);
            fails += ac_memo_check(q,uv,kind,&fresh,&memo,1);
            active[4] ^= 1;
            fails += ac_memo_check(q,uv,kind,&fresh,&memo,0);
            active[4] ^= 1;
            fails += ac_memo_check(q,uv,kind,&fresh,&memo,0);   /* the slot now holds the flipped mask */
        }
        /* A new partition invalidates every slot, even over the same geometry. */
        memcpy(moved,uv,2*NV*sizeof *moved); for (int v = 0; v < NV; v++) moved[2*v] += 40*(v%2);
        if (AsmContacts_rebuild(q,moved) || AsmContacts_rebuild(q,uv)) fails++;
        for (int kind = 0; kind < AC_MEMO_KINDS; kind++) fails += ac_memo_check(q,uv,kind,&fresh,&memo,0);
        /* A budget below the stored walk's leaf pairs must stop it as it
         * stops a fresh walk; an unbounded one replays again. */
        size_t leaf_pairs = memo.stats.leaf_pairs;
        if (leaf_pairs < 2) { fails++; leaf_pairs = 2; }
        AsmContacts_leaf_pair_budget(q,leaf_pairs-1);
        fails += ac_memo_check(q,uv,AC_MEMO_KINDS-1,&fresh,&memo,0);
        if (memo.rc != -2 || memo.stats.complete) fails++;
        AsmContacts_leaf_pair_budget(q,0);
        fails += ac_memo_check(q,uv,AC_MEMO_KINDS-1,&fresh,&memo,1);
        /* Verification: every replay above matched its geometry; a stored
         * hit that no longer does is caught and the walk measured afresh. */
        if (q->memo_mismatches) fails++;
        AcMemo *slot = q->memo+AC_MEMO_KINDS-1;
        if (!slot->valid || !slot->n) fails++;
        else {
            slot->hit[0].area = nextafter(slot->hit[0].area,INFINITY);
            fails += ac_memo_check(q,uv,AC_MEMO_KINDS-1,&fresh,&memo,0);
            if (q->memo_mismatches != 1) fails++;
            fails += ac_memo_check(q,uv,AC_MEMO_KINDS-1,&fresh,&memo,1);   /* stored afresh */
        }
        /* A complete (unmasked) walk is never kept. */
        AsmContacts_active_faces(q,NULL);
        fails += ac_memo_check(q,uv,0,&fresh,&memo,0);
        fails += ac_memo_check(q,uv,0,&fresh,&memo,0);
        replays += q->memo_replays;
    }
    if (!hits) fails++;
    fprintf(stderr,"  contact walk memo: %s (%d failures; %zu replays, %zu hits+weights compared)\n",
            fails ? "FAIL" : "ok",fails,replays,hits);
    Arena_dispose(&arena);
    return fails;
}

/* The parallel build must be the serial tree byte for byte: ids, ranges,
 * boxes and the face order.  Node padding is uninitialised, so compare
 * fields, not bytes. */
static int ac_parallel_build_control(void)
{
    enum { CELLS=8192, NF=2*CELLS, NV=3*NF };
    Arena_T arena=Arena_new();
    double *uv=ARENA_ALLOC(arena,2*NV*sizeof(double));int32_t *faces=ARENA_ALLOC(arena,3*NF*sizeof(int32_t));
    const double triangle[]={0,0,1,0,0,1};
    for(int i=0;i<NF;i++){int cell=i/2;for(int k=0;k<3;k++){int v=3*i+k;faces[v]=v;uv[2*v]=4*(cell%128)+triangle[2*k]+.25*(i%2);uv[2*v+1]=4*(cell/128)+triangle[2*k+1]+.125*(i%2);}}
    AsmContacts *one=AsmContacts_new_threads(arena,NV,faces,NF,uv,1e-7,1),*many=AsmContacts_new_threads(arena,NV,faces,NF,uv,1e-7,8);
    int fails=0;size_t differ=0;
    if(!one || !many || one->nodes!=many->nodes || one->nodes!=ac_count(NF))fails++;
    else{
        for(size_t n=0;n<one->nodes;n++){const AcNode *a=one->tree+n,*b=many->tree+n;
            differ+=a->begin!=b->begin || a->end!=b->end || a->left!=b->left || a->right!=b->right ||
                    a->lo[0]!=b->lo[0] || a->lo[1]!=b->lo[1] || a->hi[0]!=b->hi[0] || a->hi[1]!=b->hi[1];}
        for(size_t f=0;f<NF;f++)differ+=one->order[f]!=many->order[f];
        AsmContactStats a={0},b={0};
        if(differ || AsmContacts_measure(one,uv,NULL,NULL,NULL,&a) || AsmContacts_measure(many,uv,NULL,NULL,NULL,&b) || a.pairs!=b.pairs || a.area!=b.area)fails++;
    }
    fprintf(stderr,"  contact tree parallel vs serial build: %s (%zu nodes over %d faces, %zu differences)\n",fails?"FAIL":"ok",one?one->nodes:(size_t)0,NF,differ);
    Arena_dispose(&arena);return fails;
}

/* The walk budget: a tree whose leaves were stretched by a move it was not
 * rebuilt for must stop at the same leaf pair on the serial and the threaded
 * walk and report the count; a rebuild recovers compact leaves and an exact
 * answer under the same budget; an unbounded budget reproduces the sound
 * answer on the stretched tree (a different visit order, so the area sum is
 * compared to 1e-9 while the pair count is exact). */
static int ac_budget_control(void)
{
    enum { CELLS=256, NF=2*CELLS, NV=3*NF };
    Arena_T arena=Arena_new();
    double *uv=ARENA_ALLOC(arena,2*NV*sizeof(double));int32_t *faces=ARENA_ALLOC(arena,3*NF*sizeof(int32_t));
    const double triangle[]={0,0,1,0,0,1};
    /* two interleaved charts: the odd faces (chart B) start on top of the even ones (chart A) */
    for(int i=0;i<NF;i++){int cell=i/2;for(int k=0;k<3;k++){int v=3*i+k;faces[v]=v;uv[2*v]=4*(cell%16)+triangle[2*k]+.25*(i%2);uv[2*v+1]=4*(cell/16)+triangle[2*k+1]+.125*(i%2);}}
    AsmContacts *q=AsmContacts_new(arena,NV,faces,NF,uv,1e-7);int fails=0;
    if(!q){Arena_dispose(&arena);return 1;}
    AsmContactStats base={0};
    if(AsmContacts_measure(q,uv,NULL,NULL,NULL,&base) || !base.complete || base.pairs!=CELLS || !base.leaf_pairs)fails++;
    if(AsmContacts_leaves(q)==0 || AsmContacts_budget(q)!=(size_t)ASM_CONTACT_MAX_LEAF_PAIRS_PER_LEAF*AsmContacts_leaves(q))fails++;
    /* move chart B far away WITHOUT a rebuild: every leaf that mixed A and B faces now spans the move */
    for(int i=1;i<NF;i+=2)for(int k=0;k<3;k++)uv[2*(3*i+k)]+=1e5;
    AsmContactStats serial={0},threaded={0};
    AsmContacts_leaf_pair_budget(q,base.leaf_pairs+8);
    AsmContacts_threads(q,1);int rs=AsmContacts_measure(q,uv,NULL,NULL,NULL,&serial);
    AsmContacts_threads(q,8);int rt=AsmContacts_measure(q,uv,NULL,NULL,NULL,&threaded);
    if(rs!=-2 || rt!=-2 || serial.complete || threaded.complete || serial.leaf_pairs!=threaded.leaf_pairs || serial.leaf_pairs!=base.leaf_pairs+9)fails++;
    /* the rebuild recovers compact leaves: the walk completes under the default budget and is exact */
    AsmContactStats after={0},brute={0};
    AsmContacts_leaf_pair_budget(q,(size_t)ASM_CONTACT_MAX_LEAF_PAIRS_PER_LEAF*AsmContacts_leaves(q));
    if(AsmContacts_rebuild(q,uv) || AsmContacts_measure(q,uv,NULL,NULL,NULL,&after) || !after.complete)fails++;
    for(int i=0;i<NF;i++)for(int j=i+1;j<NF;j++){double ar=AsmContacts_pair(uv+6*i,uv+6*j,1e-7,NULL,NULL);if(ar>0){brute.pairs++;brute.area+=ar;}}
    if(after.pairs!=brute.pairs || brute.pairs!=0)fails++;
    /* an unbounded budget reproduces the sound answer on the stretched tree */
    for(int i=1;i<NF;i+=2)for(int k=0;k<3;k++)uv[2*(3*i+k)]-=1e5;
    AsmContacts_leaf_pair_budget(q,0);AsmContactStats again={0};
    AsmContacts_threads(q,1);
    if(AsmContacts_measure(q,uv,NULL,NULL,NULL,&again) || !again.complete || again.pairs!=base.pairs || fabs(again.area-base.area)>1e-9)fails++;
    /* the touching walk stops at the same budget */
    AsmContacts_leaf_pair_budget(q,base.leaf_pairs+8);AsmContactStats touch={0};
    if(AsmContacts_touching(q,uv,NULL,NULL,&touch)!=-2 || touch.complete || touch.leaf_pairs!=base.leaf_pairs+9)fails++;
    fprintf(stderr,"  contact walk budget: %s (%d failures; %zu leaf pairs sound, stopped at %zu on the stretched tree)\n",fails?"FAIL":"ok",fails,base.leaf_pairs,serial.leaf_pairs);
    Arena_dispose(&arena);return fails;
}

/* The external probe: a LOCAL field's candidate against the placed material its own field leaves
 * out.  Feeding the tree's own faces back in as external triangles must reproduce brute force
 * exactly -- with those faces excluded (the local field owns them) and without (it does not). */
static int ac_probe_control(void)
{
    enum { NF = 61, NV = 3*NF };
    Arena_T arena = Arena_new();
    double *uv = ARENA_ALLOC(arena,2*NV*sizeof(double)); int32_t *faces = ARENA_ALLOC(arena,3*NF*sizeof(int32_t));
    uint8_t *exclude = ARENA_CALLOC(arena,NF,1);
    const double base[] = {0,0, 2,0, 0,2};
    for (int i = 0; i < NF; i++) for (int k = 0; k < 3; k++) {
        faces[3*i+k] = 3*i+k;
        uv[6*i+2*k] = (i*13%17)*0.5+base[2*k]; uv[6*i+2*k+1] = (i*7%19)*0.5+base[2*k+1];
    }
    AsmContacts *q = AsmContacts_new(arena,NV,faces,NF,uv,1e-7);
    if (!q) { Arena_dispose(&arena); return 1; }
    int fails = 0;
    /* the external set: every third face, excluded from the tree side */
    size_t next = 0;
    double *tri = ARENA_ALLOC(arena,(size_t)NF*6*sizeof(double));
    int32_t *src = ARENA_ALLOC(arena,(size_t)NF*sizeof(int32_t));
    for (int i = 0; i < NF; i += 3) { memcpy(tri+6*next,uv+6*i,6*sizeof(double)); src[next++] = i; exclude[i] = 1; }
    for (int pass = 0; pass < 2; pass++) {
        AsmContactStats brute = {0}, got = {0};
        for (size_t e = 0; e < next; e++) for (int j = 0; j < NF; j++) {
            if (pass == 0 && exclude[j]) continue;
            double area = AsmContacts_pair(tri+6*e,uv+6*j,1e-7,NULL,NULL);
            if (area > 0) { brute.pairs++; brute.area += area; }
        }
        int rc = AsmContacts_probe(q,uv,next,tri,pass == 0 ? exclude : NULL,NULL,NULL,&got);
        if (rc || !got.complete || got.pairs != brute.pairs || fabs(got.area-brute.area) > 1e-9) fails++;
        AsmContactStats prepared={0};
        if(AsmContacts_probe_prepared(q,next,tri,pass==0?exclude:NULL,NULL,NULL,&prepared) ||
           !prepared.complete || prepared.pairs!=got.pairs || prepared.area!=got.area || prepared.leaf_pairs!=got.leaf_pairs)fails++;
        AsmContacts_prepare_exclusion(q,pass==0?exclude:NULL);
        if(AsmContacts_probe_prepared(q,next,tri,pass==0?exclude:NULL,NULL,NULL,&prepared) ||
           !prepared.complete || prepared.pairs!=got.pairs || prepared.area!=got.area || prepared.leaf_pairs>got.leaf_pairs)fails++;
        AsmContacts_prepare_exclusion(q,NULL);
        if (pass == 1 && got.pairs <= brute.pairs - next) fails++;   /* without the mask it also finds itself */
    }
    /* the budget stops it and says so */
    AsmContactStats bounded = {0};
    AsmContacts_leaf_pair_budget(q,1);
    if (AsmContacts_probe(q,uv,next,tri,NULL,NULL,NULL,&bounded) != -2 || bounded.complete || bounded.leaf_pairs != 2) fails++;
    if(AsmContacts_probe_prepared(q,next,tri,NULL,NULL,NULL,&bounded)!=-2 || bounded.complete || bounded.leaf_pairs!=2)fails++;
    AsmContacts_leaf_pair_budget(q,0);
    /* a measure afterwards is unaffected: the probe restored the mask and the boxes */
    uint8_t *moving = ARENA_CALLOC(arena,NF,1); moving[0] = 1;
    AsmContacts_active_faces(q,moving);
    AsmContactStats masked = {0}, complete = {0};
    if (AsmContacts_measure(q,uv,NULL,NULL,NULL,&masked)) fails++;
    AsmContactStats between = {0};
    AsmContacts_probe(q,uv,next,tri,NULL,NULL,NULL,&between);
    AsmContactStats masked2 = {0};
    if (AsmContacts_measure(q,uv,NULL,NULL,NULL,&masked2) || masked2.pairs != masked.pairs || masked2.area != masked.area) fails++;
    AsmContacts_active_faces(q,NULL);
    if (AsmContacts_measure(q,uv,NULL,NULL,NULL,&complete) || complete.pairs < masked.pairs) fails++;
    fprintf(stderr,"  contact external probe: %s (%d failures; %zu external triangles)\n",fails?"FAIL":"ok",fails,next);
    Arena_dispose(&arena); return fails;
}

static int ac_changed_refit_control(void)
{
    enum {NF=97,NV=3*NF};int fails=0;Arena_T arena=Arena_new();
    double *uv=ARENA_ALLOC(arena,2*NV*sizeof *uv);int32_t *faces=ARENA_ALLOC(arena,3*NF*sizeof *faces);
    uint8_t *changed=ARENA_CALLOC(arena,NF,1),*exclude=ARENA_CALLOC(arena,NF,1),*active=ARENA_CALLOC(arena,NF,1);
    const double base[6]={0,0,20,0,0,20};active[0]=1;
    for(int i=0;i<NF;i++)for(int k=0;k<3;k++){
        faces[3*i+k]=3*i+k;uv[6*i+2*k]=(i*13%17)*5+base[2*k];uv[6*i+2*k+1]=(i*7%19)*5+base[2*k+1];
    }
    AsmContacts *partial=AsmContacts_new_threads(arena,NV,faces,NF,uv,AsmContacts_material_tolerance(),4);
    AsmContacts *full=AsmContacts_new(arena,NV,faces,NF,uv,AsmContacts_material_tolerance());
    if(!partial || !full){Arena_dispose(&arena);return 1;}
    /* This exhaustive external cross-product can exceed the ordinary walk's
     * leaf budget. Test complete answers here, then budget refusal below. */
    AsmContacts_leaf_pair_budget(partial,0);AsmContacts_leaf_pair_budget(full,0);
    AsmContacts_active_faces(partial,active);AsmContacts_active_faces(full,active);
    for(int pass=0;pass<5;pass++){
        memset(changed,0,NF);memset(exclude,1,NF);AsmContacts_prepare_exclusion(partial,exclude);
        for(int i=0;i<NF;i++)if(pass && (pass==4 || i%7==pass-1)){
            changed[i]=1;
            for(int k=0;k<3;k++){
                uv[6*i+2*k]+=(pass==1?1000:pass==2?-70:3);
                uv[6*i+2*k+1]+=(pass==3?-27:2);
            }
        }
        if(AsmContacts_refit_changed(partial,uv,changed) || AsmContacts_refit(full,uv))fails++;
        if(partial->active_faces!=active || partial->probe_exclude)fails++;
        for(size_t i=0;i<full->nodes;i++)for(int d=0;d<2;d++)
            if(partial->tree[i].lo[d]!=full->tree[i].lo[d] || partial->tree[i].hi[d]!=full->tree[i].hi[d])fails++;
        memset(exclude,0,NF);AsmContactStats got={0},expected={0},brute={0};
        if(AsmContacts_probe_prepared(partial,NF,uv,exclude,NULL,NULL,&got) ||
           AsmContacts_probe_prepared(full,NF,uv,exclude,NULL,NULL,&expected) ||
           !got.complete || got.pairs!=expected.pairs || got.area!=expected.area || got.leaf_pairs!=expected.leaf_pairs){
            fails++;fprintf(stderr,"  changed refit pass %d partial/full: %zu / %zu pairs %.17g / %.17g area, %d / %d complete\n",pass,got.pairs,expected.pairs,got.area,expected.area,got.complete,expected.complete);
        }
        for(int i=0;i<NF;i++)for(int j=0;j<NF;j++){
            double area=AsmContacts_pair_faces(full,uv,(size_t)i,(size_t)j,NULL,NULL);
            if(area>0){brute.pairs++;brute.area+=area;}
        }
        if(got.pairs!=brute.pairs || fabs(got.area-brute.area)>1e-8){
            fails++;fprintf(stderr,"  changed refit pass %d probe: %zu / %zu pairs, %.17g / %.17g area\n",pass,got.pairs,brute.pairs,got.area,brute.area);
        }
        if(AsmContacts_measure(partial,uv,NULL,NULL,NULL,&got) || AsmContacts_measure(full,uv,NULL,NULL,NULL,&expected) ||
           got.pairs!=expected.pairs || got.area!=expected.area){
            fails++;fprintf(stderr,"  changed refit pass %d measure: %zu / %zu pairs, %.17g / %.17g area\n",pass,got.pairs,expected.pairs,got.area,expected.area);
        }
    }
    memset(changed,0,NF);changed[0]=1;double saved=uv[0];uv[0]=NAN;
    if(AsmContacts_refit_changed(partial,uv,changed)!=-1){fails++;fprintf(stderr,"  changed refit accepted nonfinite coordinate\n");}
    uv[0]=saved;if(AsmContacts_refit_changed(partial,uv,changed))fails++;
    if(AsmContacts_refit_changed(NULL,uv,changed)!=-1 || AsmContacts_refit_changed(partial,NULL,changed)!=-1 ||
       AsmContacts_refit_changed(partial,uv,NULL))fails++;
    AsmContacts_leaf_pair_budget(partial,1);AsmContacts_leaf_pair_budget(full,1);
    AsmContactStats limited={0},reference={0};
    if(AsmContacts_refit_changed(partial,uv,changed) || AsmContacts_refit(full,uv) ||
       AsmContacts_probe_prepared(partial,NF,uv,NULL,NULL,NULL,&limited)!=-2 ||
       AsmContacts_probe_prepared(full,NF,uv,NULL,NULL,NULL,&reference)!=-2 ||
       limited.complete || reference.complete || limited.leaf_pairs!=2 || reference.leaf_pairs!=2){fails++;fprintf(stderr,"  changed refit budget: %zu / %zu leaves, %d / %d complete\n",limited.leaf_pairs,reference.leaf_pairs,limited.complete,reference.complete);}
    fprintf(stderr,"  changed-face contact refit: %s (%d failures)\n",fails?"FAIL":"ok",fails);
    Arena_dispose(&arena);return fails;
}

#include "asm_contacts_material_test.inc"

/* The box pre-reject never refuses a pair the full test would measure:
 * random triangle pairs at unit and 1e6 scale, placed on both sides of the
 * margin (from 200 ulp overlapping to 600 ulp apart), in objective mode with
 * derivatives and the 4-vox material tolerance. */
static uint64_t ac_test_next(uint64_t *seed)
{ *seed = *seed*6364136223846793005ull+1442695040888963407ull; return *seed>>11; }

static int ac_separated_control(void)
{
    int failures = 0; uint64_t seed = 12345; size_t rejected = 0, trials = 20000;
    for (size_t trial = 0; trial < trials; trial++) {
        double base = trial%2 ? 1e6 : 0, x[6], y[6], unit = DBL_EPSILON*(base+8), xmax = -DBL_MAX, ymin = DBL_MAX;
        for (int k = 0; k < 6; k++) {
            x[k] = base+4*(double)ac_test_next(&seed)/9007199254740992.0;
            y[k] = base+4*(double)ac_test_next(&seed)/9007199254740992.0;
        }
        for (int k = 0; k < 3; k++) { xmax = x[2*k] > xmax ? x[2*k] : xmax; ymin = y[2*k] < ymin ? y[2*k] : ymin; }
        double delta = ((double)(ac_test_next(&seed)%801)-200)*unit;
        for (int k = 0; k < 3; k++) y[2*k] += xmax-ymin+delta;
        if (!ac_separated(x,y)) continue;
        rejected++;
        double ga[6], gb[6], raw = -1;
        double area = ac_material_pair(NULL,NULL,0,NULL,x,NULL,0,NULL,y,4.0,ga,gb,&raw);
        int zero = area == 0 && raw == 0;
        for (int k = 0; k < 6; k++) zero &= ga[k] == 0 && gb[k] == 0;
        if (!zero) failures++;
    }
    if (!rejected || rejected == trials) failures++;
    fprintf(stderr,"  contact box pre-reject: %zu of %zu pairs rejected, each measured apart by the full test: %s\n",
            rejected,trials,failures ? "FAIL" : "ok");
    return failures;
}

/* ac_min/ac_max bitwise == the comparison-then-CRT definitions they
 * replaced, on signed zeros, NaN, infinities, equal and random values. */
static int ac_minmax_control(void)
{
    const double special[] = {0.0, -0.0, 1.0, -1.0, 0.5, INFINITY, -INFINITY, NAN, -NAN, DBL_MIN, -DBL_MIN, 1e-310, DBL_MAX};
    const size_t ns = sizeof special / sizeof special[0];
    double values[64];
    uint64_t state = 0x9E3779B97F4A7C15ULL;
    int fails = 0;
    for (size_t i = 0; i < 64; i++) {
        state ^= state << 13; state ^= state >> 7; state ^= state << 17;
        values[i] = i < ns ? special[i] : (double)(int64_t)(state >> 40) / 4096.0;
    }
    for (size_t i = 0; i < 64; i++) for (size_t j = 0; j < 64; j++) {
        double a = values[i], b = values[j];
        double lo = a < b ? a : b < a ? b : fmin(a,b), hi = a > b ? a : b > a ? b : fmax(a,b);
        double got_lo = ac_min(a,b), got_hi = ac_max(a,b);
        if (memcmp(&lo,&got_lo,sizeof lo) || memcmp(&hi,&got_hi,sizeof hi)) {
            if (fails < 3) fprintf(stderr, "ac_minmax_control: (%g, %g) differs\n", a, b);
            fails++;
        }
    }
    return fails ? 1 : 0;
}

int AsmContacts_selftest(void)
{
    int fails = ac_minmax_control()+ac_separated_control()+acb_selftest()+acm_selftest()+ac_changed_refit_control()+ac_material_policy_control()+act_control()+ac_repartition_control()+ac_threaded_control()+ac_pooled_control()+ac_memo_control()+ac_parallel_build_control()+ac_budget_control()+ac_probe_control()+acc_control();
    double a[] = {0,0, 2,0, 0,2}, b[] = {.3,.2, 2.3,.2, .3,2.2}, ga[6], gb[6];
    double area = AsmContacts_pair(a,b,1e-11,ga,gb), sum = gb[0]+gb[2]+gb[4];
    for (int k = 0; k < 3; k++) b[2*k] += 1e-6;
    double plus = AsmContacts_pair(a,b,1e-11,NULL,NULL);
    for (int k = 0; k < 3; k++) b[2*k] -= 2e-6;
    double minus = AsmContacts_pair(a,b,1e-11,NULL,NULL);
    if (fabs(area-1.125) > 1e-12 || fabs((plus-minus)/2e-6-sum) > 1e-8) fails++;
    /* Independent finite differences for every shape coordinate. */
    double original_b[] = {.3,.2,2.3,.2,.3,2.2}; memcpy(b,original_b,sizeof b);
    AsmContacts_pair(a,b,1e-11,ga,gb);
    for (int side = 0; side < 2; side++) for (int k = 0; k < 6; k++) {
        double *p = side ? b : a, saved = p[k];
        p[k] = saved+1e-6; plus = AsmContacts_pair(a,b,1e-11,NULL,NULL);
        p[k] = saved-1e-6; minus = AsmContacts_pair(a,b,1e-11,NULL,NULL); p[k] = saved;
        if (fabs((plus-minus)/2e-6-(side ? gb[k] : ga[k])) > 2e-7) fails++;
    }
    double tiny[] = {.5,.5,.5+1e-8,.5,.5,.5+1e-8};
    if (!(AsmContacts_pair(a,tiny,1e-11,NULL,NULL) > 0) || AsmContacts_pair(a,tiny,1e-7,NULL,NULL) != 0) fails++;
    double touch[] = {2,0,3,0,2,1}; if (AsmContacts_pair(a,touch,0,NULL,NULL) != 0) fails++;
    double shared[] = {0,0,2,0,1,1}; if (!(AsmContacts_pair(a,shared,1e-7,NULL,NULL) > 0)) fails++;
    /* Refinement must not multiply physical overlap area or its rigid force.
     * The 4, 64 and 1024 triangle versions cover the exact same two squares. */
    for (int n = 1; n <= 16; n *= 4) {
        Arena_T ar = Arena_new(); size_t side_nv = (size_t)(n+1)*(n+1), nv = 2*side_nv, nf = (size_t)4*n*n;
        double *u = ARENA_ALLOC(ar,2*nv*sizeof(double)), *gradient = ARENA_ALLOC(ar,2*nv*sizeof(double));
        int32_t *tri = ARENA_ALLOC(ar,3*nf*sizeof(int32_t)); size_t at = 0;
        for (int side = 0; side < 2; side++) for (int y = 0; y <= n; y++) for (int x = 0; x <= n; x++) {
            size_t v = side*side_nv+(size_t)y*(n+1)+x;
            u[2*v] = (double)x/n+.2*side; u[2*v+1] = (double)y/n+.1*side;
            if (x < n && y < n) { int32_t f6[] = {(int32_t)v,(int32_t)v+1,(int32_t)v+n+2,(int32_t)v,(int32_t)v+n+2,(int32_t)v+n+1}; memcpy(tri+at,f6,sizeof f6); at += 6; }
        }
        AsmContacts *query = AsmContacts_new(ar,nv,tri,nf,u,1e-11); AsmContactStats measured;
        if (!query || AsmContacts_measure(query,u,gradient,NULL,NULL,&measured)) fails++;
        else {
            double gx = 0, gy = 0; for (size_t v = side_nv; v < nv; v++) { gx += gradient[2*v]; gy += gradient[2*v+1]; }
            if (fabs(measured.area-.72) > 1e-11 || fabs(gx+.9) > 1e-10 || fabs(gy+.8) > 1e-10) fails++;
        }
        Arena_dispose(&ar);
    }
    /* Complete BVH against independent all-pairs enumeration; then refit
     * after nonuniform moves, including shared-vertex/self-chart triangles. */
    enum { NF = 97, NV = 3*NF };
    double uv[2*NV]; int32_t f[3*NF];
    for (int i = 0; i < NF; i++) for (int k = 0; k < 3; k++) {
        f[3*i+k] = 3*i+k;
        uv[6*i+2*k] = (i*17%29)*.25+a[2*k]; uv[6*i+2*k+1] = (i*11%23)*.25+a[2*k+1];
    }
    Arena_T arena = Arena_new(); AsmContacts *q = AsmContacts_new(arena,NV,f,NF,uv,1e-7);
    if (!q) fails++;
    else for (int pass = 0; pass < 3; pass++) {
        AsmContactStats measured, brute = {0};
        if (pass==1) for (int i = 0; i < NV; i++) uv[2*i] += .013*(i%7);
        if (pass==2){
            for(int i=0;i<NV;i++)uv[2*i]+=100*(i/3%5);
            if(AsmContacts_rebuild(q,uv))fails++;
        }
        for (int i = 0; i < NF; i++) for (int j = i+1; j < NF; j++) {
            double ar = AsmContacts_pair(uv+6*i,uv+6*j,1e-7,NULL,NULL);
            if (ar > 0) { brute.pairs++; brute.area += ar; }
        }
        if (AsmContacts_measure(q,uv,NULL,NULL,NULL,&measured) || measured.pairs != brute.pairs || fabs(measured.area-brute.area) > 1e-9) fails++;
    }
    if (q) {
        uint8_t moving[NF]; AsmContactStats expected={0},affected,complete;
        for(int i=0;i<NF;i++)moving[i]=(uint8_t)(i%7==0);
        for(int i=0;i<NF;i++)for(int j=i+1;j<NF;j++)if(moving[i]||moving[j]){
            double area=AsmContacts_pair(uv+6*i,uv+6*j,1e-7,NULL,NULL);
            if(area>0){expected.pairs++;expected.area+=area;}
        }
        AsmContacts_active_faces(q,moving);
        if(AsmContacts_rebuild(q,uv))fails++;
        if(AsmContacts_measure(q,uv,NULL,NULL,NULL,&affected) || affected.pairs!=expected.pairs || fabs(affected.area-expected.area)>1e-9)fails++;
        AsmContacts_active_faces(q,NULL);
        if(AsmContacts_measure(q,uv,NULL,NULL,NULL,&complete) || complete.pairs<=affected.pairs)fails++;
    }
    Arena_dispose(&arena);
    fprintf(stderr,"  native complete triangle contacts: %s (%d failures)\n",fails ? "FAIL" : "ok",fails);
    return fails;
}
