/* asm_verdict.c -- two lattices from the assembled sheet + the ribbon_verdict spawn.
 *
 * The assembled sheet is a triangle mesh with continuous (u,v).  Two consumers
 * want a du x dv lattice of it, and they want different things:
 *
 *   1. THE ASSEMBLED RIBBON (assembled_ribbon.*) is the solid stage's input
 *      contract (`quadribbon --solidify`): ROWS ARE Z-SLICES (quad_strip and
 *      the solid stage measure per z-plane), one claimant per cell (the
 *      largest chart), the cube-seam rows closed by interpolation between
 *      seam-related charts, faces only under the wrap gate, world frame,
 *      with the sidecars the solid stage reads (support, lifted phase,
 *      provenance, lane, material, stats).
 *
 *   2. THE VERDICT LATTICE (verdict_lattice.*) is what ribbon_verdict
 *      measures, and it must expose the assembly's own errors instead of
 *      tidying them away (reviewer, 2026-09-09: the 4x21x21 carried 899,175
 *      contested cells and passed WINDING/single-cover):
 *        - rows are the sheet's OWN v, so a chart placed at the wrong v is
 *          measured there, not re-levelled to its world z;
 *        - every claimant a layer apart from the cell's primary claimant is
 *          its own vertex: a wrap drawn on top of another fails single-cover;
 *        - the material identity is the LINEAGE of measured joins (the chart
 *          graph over accepted, kept, unswitched, non-weak relations): a face
 *          joining two lineages is a raster coincidence, never a join;
 *        - faces are emitted wherever the assembly CLAIMS continuity (the
 *          same chart, or an accepted seam join) whatever the 3-D step, so a
 *          wrong join reaches TANGLE; two unrelated charts are joined only as
 *          a geometric continuation under the wrap gate;
 *        - radial measurements use the local axis frame of the configured
 *          table. A curved axis changes distances, so source-continuity
 *          certification uses the emitted points in world coordinates.
 */
#include "asm_verdict.h"

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/mesh_bin.h"
#include "../common/pipeline_constants.h"
#include "../common/union_find.h"
#include "../common/ves_platform.h"

#define AV_PI 3.14159265358979323846

void AsmVerdict_default_opts(AsmVerdictOpts *o)
{
    memset(o, 0, sizeof *o);
    o->du = 2.0; o->dv = 2.0;
    o->umb_y = 0.0; o->umb_x = 0.0;
    o->pitch = 9.5;
    o->core_radius = 0.0;
    o->verdict_exe = "build/Release/ribbon_verdict.exe";
    o->include_blobs_in_source = 0;
    o->axis = NULL;
}

typedef struct AvClaim { int64_t key; int32_t chart; float p[3]; float n[3]; double area; } AvClaim;

static uint32_t av_float_order(float value)
{
    uint32_t bits; memcpy(&bits,&value,sizeof bits);
    return bits & UINT32_C(0x80000000) ? ~bits : bits ^ UINT32_C(0x80000000);
}

static int av_cmp(const void *x, const void *y)
{
    const AvClaim *a = x, *b = y;
    if (a->key != b->key) return a->key < b->key ? -1 : 1;
    if (a->area != b->area) return a->area > b->area ? -1 : 1;   /* largest chart first */
    if (a->chart != b->chart) return (a->chart > b->chart) - (a->chart < b->chart);
    /* One chart can have several physical claimants at the same raster
     * cell. qsort is unstable: an unrelated appended strip used to change
     * the observed winner by up to 34 vox. Give distinct samples a total
     * order, including their normals, before any claimant is selected. */
    for (int k = 0; k < 6; k++) {
        uint32_t aa = av_float_order(k < 3 ? a->p[k] : a->n[k-3]);
        uint32_t bb = av_float_order(k < 3 ? b->p[k] : b->n[k-3]);
        if (aa != bb) return aa < bb ? -1 : 1;
    }
    return 0;
}

static int av_cmp_double(const void *x, const void *y)
{
    double a = *(const double *)x, b = *(const double *)y;
    return (a > b) - (a < b);
}

static int av_cmp_i64(const void *x, const void *y)
{
    int64_t a = *(const int64_t *)x, b = *(const int64_t *)y;
    return (a > b) - (a < b);
}

static int av_write_raw(const char *path, const void *v, size_t esize, size_t n)
{
    FILE *fp = fopen(path, "wb");
    if (!fp) return -1;
    size_t w = fwrite(v, esize, n, fp);
    fclose(fp);
    return w == n ? 0 : -1;
}

static int64_t av_key(long row, long col) { return ((int64_t)(row + (1 << 20)) << 24) | (int64_t)(col + (1 << 22)); }
static long av_key_row(int64_t key) { return (long)(key >> 24) - (1 << 20); }
static long av_key_col(int64_t key) { return (long)(key & 0xFFFFFF) - (1 << 22); }

static double av_dist(const float *a, const float *b)
{
    return sqrt(((double)a[0]-b[0])*((double)a[0]-b[0]) + ((double)a[1]-b[1])*((double)a[1]-b[1]) + ((double)a[2]-b[2])*((double)a[2]-b[2]));
}

static double av_normal_distance(const AvClaim *a, const AvClaim *b)
{
    double dot = 0, norm2 = 0;
    for (int k = 0; k < 3; k++) {
        double normal = a->n[k];
        dot += ((double)b->p[k]-a->p[k])*normal; norm2 += normal*normal;
    }
    /* Interpolated unit vertex normals generally have length below one.
     * Normalize the measurement, leaving the sample and its ordering intact.
     * With no normal direction, do not erase a distinct physical claimant. */
    return norm2 > 0 ? fabs(dot)/sqrt(norm2) : av_dist(a->p,b->p);
}

/* lower bound of `key` in cells sorted by key */
static size_t av_lower(const AvClaim *cells, size_t lo, size_t hi, int64_t key)
{
    while (lo < hi) { size_t mid = (lo + hi) / 2; if (cells[mid].key < key) lo = mid + 1; else hi = mid; }
    return lo;
}

/* Lifted phase: atan2 about the axis, made continuous along every row (the
 * prediction across a gap is the arc-length step over the local radius, in
 * the row's measured winding sense), then rows aligned to their predecessor
 * by the median of the overlap.  The solid stage checks the lift against
 * atan2 modulo 2 pi and uses it as the angular coordinate of its fills. */
static void av_lift_phase(const AvClaim *cells, size_t n, double ay, double ax, double du, float *phase)
{
    size_t i = 0;
    /* pass 1: per-row continuity */
    while (i < n) {
        size_t j = i;
        long row = av_key_row(cells[i].key);
        while (j < n && av_key_row(cells[j].key) == row) j++;
        /* winding sense of the row: sign of wrapped atan2 steps between close cells */
        double sense = 0.0;
        for (size_t k = i + 1; k < j; k++) {
            long dc = av_key_col(cells[k].key) - av_key_col(cells[k-1].key);
            if (dc > 3) continue;
            double a0 = atan2(cells[k-1].p[2] - ax, cells[k-1].p[1] - ay), a1 = atan2(cells[k].p[2] - ax, cells[k].p[1] - ay);
            double d = a1 - a0;
            while (d > AV_PI) d -= 2.0 * AV_PI;
            while (d < -AV_PI) d += 2.0 * AV_PI;
            sense += d;
        }
        sense = sense >= 0.0 ? 1.0 : -1.0;
        double prev = 0.0, prev_u = 0.0, prev_r = 1.0;
        int have = 0;
        for (size_t k = i; k < j; k++) {
            double raw = atan2(cells[k].p[2] - ax, cells[k].p[1] - ay);
            double r = hypot(cells[k].p[2] - ax, cells[k].p[1] - ay);
            double u = (double)av_key_col(cells[k].key) * du;
            double lifted = raw;
            if (have) {
                double predicted = prev + sense * (u - prev_u) / fmax(0.5 * (r + prev_r), 1.0);
                lifted = raw + 2.0 * AV_PI * nearbyint((predicted - raw) / (2.0 * AV_PI));
            }
            phase[k] = (float)lifted;
            prev = lifted; prev_u = u; prev_r = r; have = 1;
        }
        i = j;
    }
    /* pass 2: align each row to the previous one over their common columns */
    size_t ri = 0;
    size_t prev_i = 0, prev_j = 0;
    int have_prev = 0;
    while (ri < n) {
        size_t rj = ri;
        long row = av_key_row(cells[ri].key);
        while (rj < n && av_key_row(cells[rj].key) == row) rj++;
        if (have_prev) {
            /* median of (this - prev) over matching columns, quantized to 2 pi */
            double diffs[4096];
            size_t nd = 0, a = prev_i, b = ri;
            while (a < prev_j && b < rj && nd < 4096) {
                long ca = av_key_col(cells[a].key), cb = av_key_col(cells[b].key);
                if (ca < cb) a++;
                else if (cb < ca) b++;
                else { diffs[nd++] = (double)phase[b] - (double)phase[a]; a++; b++; }
            }
            if (nd > 0) {
                for (size_t x = 1; x < nd; x++) { double v = diffs[x]; size_t y = x; while (y > 0 && diffs[y-1] > v) { diffs[y] = diffs[y-1]; y--; } diffs[y] = v; }
                double shift = 2.0 * AV_PI * nearbyint(diffs[nd/2] / (2.0 * AV_PI));
                if (shift != 0.0) for (size_t k = ri; k < rj; k++) phase[k] = (float)((double)phase[k] - shift);
            }
        }
        prev_i = ri; prev_j = rj; have_prev = 1;
        ri = rj;
    }
}

/* ---- shared pieces ------------------------------------------------------------ */

/* Rasterize the sheet's faces into (column = u/du, row = w/dv) cells, w = the
 * sheet's v (rows_v) or the world z (z-slices): every cell whose centre lies
 * inside a face's (u, w) footprint becomes one claim with the 3-D point and
 * normal interpolated on the face (on a z-slice the point sits exactly on the
 * plane).  Claims are returned unsorted. */
static AvClaim *av_rasterize(Arena_T arena, const MeshBinData *sheet, const int32_t *sheet_chart, const float *sheet_nrm,
                             const AsmRun *run, double du, double dv, int rows_v,
                             size_t *out_n, size_t *out_cap, long *out_row_min, long *out_row_max)
{
    size_t cap = 1u << 16, n = 0;
    AvClaim *cl = ARENA_ALLOC(arena, cap * sizeof(AvClaim));
    long row_min = 1L << 30, row_max = -(1L << 30);
    for (size_t f = 0; f < sheet->nf; f++) {
        const int32_t *fv = &sheet->faces[f*3];
        double u[3], w[3];
        for (int k = 0; k < 3; k++) {
            u[k] = MeshBin_uv(sheet,(size_t)fv[k]*2);
            w[k] = rows_v ? MeshBin_uv(sheet,(size_t)fv[k]*2+1) : sheet->verts[(size_t)fv[k]*3];
        }
        double area = (u[1]-u[0])*(w[2]-w[0]) - (w[1]-w[0])*(u[2]-u[0]);
        if (fabs(area) < 1e-9) continue;
        long c0 = (long)ceil(fmin(u[0], fmin(u[1], u[2])) / du), c1 = (long)floor(fmax(u[0], fmax(u[1], u[2])) / du);
        long r0 = (long)ceil(fmin(w[0], fmin(w[1], w[2])) / dv), r1 = (long)floor(fmax(w[0], fmax(w[1], w[2])) / dv);
        double chart_area = run->charts[(size_t)sheet_chart[fv[0]]].area3d;
        for (long r = r0; r <= r1; r++) {
            double pw = (double)r * dv;
            for (long c = c0; c <= c1; c++) {
                double pu = (double)c * du;
                double l0 = ((u[1]-pu)*(w[2]-pw) - (w[1]-pw)*(u[2]-pu)) / area;
                double l1 = ((u[2]-pu)*(w[0]-pw) - (w[2]-pw)*(u[0]-pu)) / area;
                double l2 = 1.0 - l0 - l1;
                if (l0 < -1e-6 || l1 < -1e-6 || l2 < -1e-6) continue;
                if (n == cap) { AvClaim *nc = ARENA_ALLOC(arena, cap * 2 * sizeof(AvClaim)); memcpy(nc, cl, n * sizeof(AvClaim)); cl = nc; cap *= 2; }
                AvClaim *q = &cl[n++];
                q->key = av_key(r, c);
                q->chart = sheet_chart[fv[0]];
                q->area = chart_area;
                for (int d = 0; d < 3; d++)
                    q->p[d] = (float)(l0 * sheet->verts[(size_t)fv[0]*3+(size_t)d] + l1 * sheet->verts[(size_t)fv[1]*3+(size_t)d] + l2 * sheet->verts[(size_t)fv[2]*3+(size_t)d]);
                if (!rows_v) q->p[0] = (float)pw;   /* exactly on the slice plane */
                for (int d = 0; d < 3; d++)
                    q->n[d] = (float)(l0 * sheet_nrm[(size_t)fv[0]*3+(size_t)d] + l1 * sheet_nrm[(size_t)fv[1]*3+(size_t)d] + l2 * sheet_nrm[(size_t)fv[2]*3+(size_t)d]);
                if (r < row_min) row_min = r;
                if (r > row_max) row_max = r;
            }
        }
    }
    *out_n = n; *out_cap = cap; *out_row_min = row_min; *out_row_max = row_max;
    return cl;
}

static double av_sample_tolerance(const float *a, const float *b, double step)
{
    double magnitude = 0;
    for (int k = 0; k < 3; k++) magnitude = fmax(magnitude,fmax(fabs(a[k]),fabs(b[k])));
    return fmax(.001*step,4*FLT_EPSILON*magnitude);
}

/* A measured trim strip is a surface between corresponding source runs.
 * Its raster footprint may cross four empty cells at 45 degrees even when
 * the physical trim is only 4.4 vox. Construct that surface before sampling
 * it, independently of the arbitrary lattice phase. Original samples win
 * ownership; generated samples carry support zero in both consumers. */
static int av_strip_shape(const MeshBinData *sheet, const size_t v[3], double *orientation)
{
    double ax = MeshBin_uv(sheet,2*v[1])-MeshBin_uv(sheet,2*v[0]), ay = MeshBin_uv(sheet,2*v[1]+1)-MeshBin_uv(sheet,2*v[0]+1);
    double bx = MeshBin_uv(sheet,2*v[2])-MeshBin_uv(sheet,2*v[0]), by = MeshBin_uv(sheet,2*v[2]+1)-MeshBin_uv(sheet,2*v[0]+1);
    double det = ax*by-ay*bx;
    *orientation = det;
    if (!isfinite(det) || fabs(det) < 1e-9) return 0;
    double g00 = 0, g01 = 0, g11 = 0;
    for (int k = 0; k < 3; k++) {
        double p = (double)sheet->verts[3*v[1]+k]-sheet->verts[3*v[0]+k];
        double q = (double)sheet->verts[3*v[2]+k]-sheet->verts[3*v[0]+k];
        double a = (p*by-q*ay)/det, b = (q*ax-p*bx)/det;
        g00 += a*a; g01 += a*b; g11 += b*b;
    }
    double spread = hypot(g00-g11,2*g01), lo = .5*(g00+g11-spread), hi = .5*(g00+g11+spread);
    return isfinite(lo) && isfinite(hi) && lo >= .25 && hi <= 4;
}

static void av_raster_strip_triangle(Arena_T arena, const MeshBinData *sheet, const float *normals,
                                     const size_t vertex[3], const int32_t owner[3], int32_t chart_a, int32_t chart_b,
                                     double du, double dv, int rows_v, AvClaim **claims, size_t *n, size_t *cap,
                                     long *row_min, long *row_max)
{
    double u[3], w[3];
    for (int k = 0; k < 3; k++) {
        u[k] = MeshBin_uv(sheet,2*vertex[k]);
        w[k] = rows_v ? MeshBin_uv(sheet,2*vertex[k]+1) : sheet->verts[3*vertex[k]];
    }
    double area = (u[1]-u[0])*(w[2]-w[0])-(w[1]-w[0])*(u[2]-u[0]);
    if (fabs(area) < 1e-9) return;
    long c0 = (long)ceil(fmin(u[0],fmin(u[1],u[2]))/du), c1 = (long)floor(fmax(u[0],fmax(u[1],u[2]))/du);
    long r0 = (long)ceil(fmin(w[0],fmin(w[1],w[2]))/dv), r1 = (long)floor(fmax(w[0],fmax(w[1],w[2]))/dv);
    for (long r = r0; r <= r1; r++) for (long c = c0; c <= c1; c++) {
        double pu = c*du, pw = r*dv;
        double l[3]; l[0] = ((u[1]-pu)*(w[2]-pw)-(w[1]-pw)*(u[2]-pu))/area;
        l[1] = ((u[2]-pu)*(w[0]-pw)-(w[2]-pw)*(u[0]-pu))/area; l[2] = 1-l[0]-l[1];
        if (l[0] < -1e-6 || l[1] < -1e-6 || l[2] < -1e-6) continue;
        if (*n == *cap) {
            AvClaim *next = ARENA_ALLOC(arena,2*(*cap)*sizeof *next);
            memcpy(next,*claims,(*n)*sizeof *next); *claims = next; *cap *= 2;
        }
        AvClaim *q = &(*claims)[(*n)++]; q->key = av_key(r,c); q->area = -1;
        double b_weight = 0;
        for (int k = 0; k < 3; k++) if (owner[k] == chart_b) b_weight += l[k];
        q->chart = b_weight > .5 ? chart_b : chart_a;
        for (int d = 0; d < 3; d++) {
            q->p[d] = (float)(l[0]*sheet->verts[3*vertex[0]+d]+l[1]*sheet->verts[3*vertex[1]+d]+l[2]*sheet->verts[3*vertex[2]+d]);
            q->n[d] = (float)(l[0]*normals[3*vertex[0]+d]+l[1]*normals[3*vertex[1]+d]+l[2]*normals[3*vertex[2]+d]);
        }
        if (!rows_v) q->p[0] = (float)pw;
        if (r < *row_min) *row_min = r; if (r > *row_max) *row_max = r;
    }
}

static void av_source_strips(Arena_T arena, const AsmRun *run, const MeshBinData *sheet,
                              const int32_t *sheet_chart, const float *normals, double du, double dv, int rows_v,
                              AvClaim **claims, size_t *n, size_t *cap, long *row_min, long *row_max)
{
    size_t *offset = ARENA_ALLOC(arena,(run->n_charts ? run->n_charts : 1)*sizeof *offset);
    for (size_t i = 0; i < run->n_charts; i++) offset[i] = SIZE_MAX;
    for (size_t i = 0; i < sheet->nv; i++) if (sheet_chart[i] >= 0 && (size_t)sheet_chart[i] < run->n_charts && offset[sheet_chart[i]] == SIZE_MAX)
        offset[sheet_chart[i]] = i;
    size_t before = *n, strips = 0;
    for (size_t i = 0; i < run->n_rels; i++) {
        const AsmRelation *r = &run->rels[i];
        if (!AsmRel_is_join(r) || !(r->continuity & ASM_CONT_SOURCE) || offset[r->a] == SIZE_MAX || offset[r->b] == SIZE_MAX) continue;
        double rms; size_t support;
        if (!AsmContinuity_measure(run,r,&rms,&support)) continue;
        const AsmCorr *previous = NULL;
        for (int32_t k = 0; k < r->corr_count; k++) {
            const AsmCorr *pair = &run->corr[(size_t)r->corr_first+k];
            if (pair->valid != 3 || pair->run < 0 || pair->va < 0 || pair->vb < 0 ||
                (size_t)pair->va >= run->charts[r->a].nv || (size_t)pair->vb >= run->charts[r->b].nv) { previous = NULL; continue; }
            if (!(AsmContinuity_pair_error(run,r,pair) <= ASM_CONTINUITY_GATE(r->rms))) { previous = NULL; continue; }
            size_t a = offset[r->a]+(size_t)pair->va, b = offset[r->b]+(size_t)pair->vb;
            if (a >= sheet->nv || b >= sheet->nv || av_dist(sheet->verts+3*a,sheet->verts+3*b) > 8) { previous = NULL; continue; }
            if (previous && previous->run == pair->run) {
                const size_t first[3] = {offset[r->a]+(size_t)previous->va,offset[r->b]+(size_t)previous->vb,b};
                const size_t second[3] = {first[0],b,a};
                const int32_t first_owner[3] = {r->a,r->b,r->b}, second_owner[3] = {r->a,r->b,r->a};
                double d1, d2;
                if (av_strip_shape(sheet,first,&d1) && av_strip_shape(sheet,second,&d2) && d1*d2 > 0) {
                    av_raster_strip_triangle(arena,sheet,normals,first,first_owner,r->a,r->b,du,dv,rows_v,claims,n,cap,row_min,row_max);
                    av_raster_strip_triangle(arena,sheet,normals,second,second_owner,r->a,r->b,du,dv,rows_v,claims,n,cap,row_min,row_max);
                    strips++;
                }
            }
            previous = pair;
        }
    }
    if (strips) fprintf(stderr,"  [verdict] SOURCE trim surfaces: %zu strips, %zu generated claims, rows %s\n",strips,*n-before,rows_v ? "v" : "z");
}

#define av_rel_is_join AsmRel_is_join

/* sorted (min, max) chart-pair keys of the joins */
static int64_t *av_join_keys(Arena_T arena, const AsmRun *run, size_t *out_n)
{
    int64_t *relkey = ARENA_ALLOC(arena, (run->n_rels ? run->n_rels : 1) * sizeof(int64_t));
    size_t nrel = 0;
    for (size_t r = 0; r < run->n_rels; r++) {
        const AsmRelation *rl = &run->rels[r];
        if (!av_rel_is_join(rl)) continue;
        int32_t a = rl->a < rl->b ? rl->a : rl->b, b = rl->a < rl->b ? rl->b : rl->a;
        relkey[nrel++] = ((int64_t)a << 32) | (int64_t)b;
    }
    qsort(relkey, nrel, sizeof(int64_t), av_cmp_i64);
    *out_n = nrel;
    return relkey;
}

static int av_joined(const int64_t *relkey, size_t nrel, int32_t x, int32_t y)
{
    if (x == y) return 1;
    int32_t a = x < y ? x : y, b = x < y ? y : x;
    int64_t key = ((int64_t)a << 32) | (int64_t)b;
    size_t lo = 0, hi = nrel;
    while (lo < hi) { size_t mid = (lo + hi) / 2; if (relkey[mid] < key) lo = mid + 1; else hi = mid; }
    return lo < nrel && relkey[lo] == key;
}

typedef struct AvSeamClaim { int64_t key; const AsmRelation *relation; } AvSeamClaim;

static int av_cmp_seam_claim(const void *p, const void *q)
{
    int64_t a = ((const AvSeamClaim *)p)->key, b = ((const AvSeamClaim *)q)->key;
    return (a > b)-(a < b);
}

static AvSeamClaim *av_seam_claims(Arena_T arena, const AsmRun *run, size_t *n)
{
    AvSeamClaim *claims = ARENA_ALLOC(arena,(run->n_rels ? run->n_rels : 1)*sizeof *claims); *n = 0;
    for (size_t i = 0; i < run->n_rels; i++) {
        const AsmRelation *r = &run->rels[i];
        /* A raster-derived layout join has no ordered boundary witnesses.
         * It still establishes a lineage, whose contacts must pass the
         * ordinary geometric edge check. It cannot override that check or
         * a failed source relation between the same charts. */
        if (!AsmRel_is_join(r) || !(r->continuity & ASM_CONT_SOURCE)) continue;
        int32_t a = r->a < r->b ? r->a : r->b, b = r->a < r->b ? r->b : r->a;
        claims[(*n)++] = (AvSeamClaim){((int64_t)a<<32)|b,r};
    }
    qsort(claims,*n,sizeof *claims,av_cmp_seam_claim); return claims;
}

static int av_near_source_run(const AsmRun *run, const AsmRelation *r, int32_t chart, const float p[3], double radius2)
{
    int is_a = chart == r->a;
    const AsmChart *c = &run->charts[chart];
    const AsmCorr *previous = NULL;
    for (int32_t k = 0; k < r->corr_count; k++) {
        const AsmCorr *pair = &run->corr[(size_t)r->corr_first+k];
        if (pair->valid != 3 || pair->run < 0) { previous = NULL; continue; }
        int32_t v = is_a ? pair->va : pair->vb;
        if (v < 0 || (size_t)v >= c->nv) { previous = NULL; continue; }
        const float *q = c->xyz+(size_t)v*3;
        double d2 = 0;
        for (int j = 0; j < 3; j++) { double d = (double)p[j]-q[j]; d2 += d*d; }
        if (d2 <= radius2) return 1;
        if (previous && previous->run == pair->run) {
            int32_t u = is_a ? previous->va : previous->vb;
            const float *a = c->xyz+(size_t)u*3;
            double dot = 0, len2 = 0;
            for (int j = 0; j < 3; j++) { double e = (double)q[j]-a[j]; len2 += e*e; dot += ((double)p[j]-a[j])*e; }
            double t = len2 > 0 ? fmax(0,fmin(1,dot/len2)) : 0; d2 = 0;
            for (int j = 0; j < 3; j++) { double d = (double)p[j]-a[j]-t*((double)q[j]-a[j]); d2 += d*d; }
            if (d2 <= radius2) return 1;
        }
        previous = pair;
    }
    return 0;
}

/* A chart-level relation identifies a particular boundary run. It does not
 * authorize connections wherever the interiors of those charts meet in UV.
 * This tests only where the source claim applies: a badly placed run still
 * emits its long edges, which the continuity/tangle audits must reject. */
static int av_seam_claimed_at(const AsmRun *run, const AvSeamClaim *claims, size_t n,
                             const AvClaim *a, const AvClaim *b, double reach)
{
    int32_t lo = a->chart < b->chart ? a->chart : b->chart, hi = a->chart < b->chart ? b->chart : a->chart;
    int64_t key = ((int64_t)lo<<32)|hi;
    size_t l = 0, h = n;
    while (l < h) { size_t m = l+(h-l)/2; if (claims[m].key < key) l = m+1; else h = m; }
    for (; l < n && claims[l].key == key; l++) {
        const AsmRelation *r = claims[l].relation;
        if (av_near_source_run(run,r,a->chart,a->p,reach*reach) &&
            av_near_source_run(run,r,b->chart,b->p,reach*reach)) return 1;
    }
    return 0;
}

/* Lineage: the component of the chart graph over the joins -- the material
 * identity the measured seams establish.  A component placed by layer hops,
 * or a cut half, is its own lineage. */
static int32_t *av_lineages(Arena_T arena, const AsmRun *run)
{
    size_t nc = run->n_charts;
    UnionFind uf = UF_new(arena, (int32_t)nc);
    for (size_t r = 0; r < run->n_rels; r++) {
        const AsmRelation *rl = &run->rels[r];
        if (!av_rel_is_join(rl)) continue;
        if (!AsmChart_in_layout(&run->charts[(size_t)rl->a]) || !AsmChart_in_layout(&run->charts[(size_t)rl->b])) continue;
        uf_union(&uf, rl->a, rl->b);
    }
    int32_t *lin = ARENA_ALLOC(arena, (nc ? nc : 1) * sizeof(int32_t));
    int32_t *root_id = ARENA_ALLOC(arena, (nc ? nc : 1) * sizeof(int32_t));
    for (size_t i = 0; i < nc; i++) root_id[i] = -1;
    int32_t n = 0;
    for (size_t i = 0; i < nc; i++) {
        int32_t rt = uf_find(&uf, (int32_t)i);
        if (root_id[rt] < 0) root_id[rt] = n++;
        lin[i] = root_id[rt];
    }
    return lin;
}

/* p <- its axis frame (s, e1 . (p - foot), e2 . (p - foot)), in place */
static void av_to_frame(const AsmAxis *axis, float *p)
{
    double q[3] = { p[0], p[1], p[2] }, out[3];
    AsmAxis_frame(axis, q, out);
    for (int k = 0; k < 3; k++) p[k] = (float)out[k];
}

/* ---- 1. the assembled ribbon (rows = z-slices; the --solidify contract) ------- */

static size_t av_solid_quad(const AvClaim *cells, const size_t q[4], double lim, int32_t faces[6])
{
    const float *p0 = cells[q[0]].p, *p1 = cells[q[1]].p, *p2 = cells[q[2]].p, *p3 = cells[q[3]].p;
    if (!(av_dist(p0,p1) < lim && av_dist(p0,p2) < lim && av_dist(p1,p3) < lim && av_dist(p2,p3) < lim)) return 0;
    /* A long chosen diagonal does not establish a hole when the other
     * diagonal and all four sides pass the same physical step gate. Keep
     * existing triangulations; use the alternate only to recover that cell. */
    int flip = !(av_dist(p0,p3) < lim);
    if (flip && !(av_dist(p1,p2) < lim)) return 0;
    faces[0] = (int32_t)q[0]; faces[1] = (int32_t)q[1]; faces[2] = (int32_t)q[flip ? 2 : 3];
    faces[3] = (int32_t)q[flip ? 1 : 0]; faces[4] = (int32_t)q[3]; faces[5] = (int32_t)q[2];
    return 2;
}

static int av_partial_triangle(const AvClaim *cells, size_t nvert, size_t i, int left, size_t tri[3]);

static void av_append_triangle(Arena_T arena, int32_t **faces, size_t *nf, size_t *cap, const size_t tri[3])
{
    if (*nf == *cap) {
        size_t next = *cap*2;
        int32_t *grown = ARENA_ALLOC(arena,next*3*sizeof(int32_t));
        memcpy(grown,*faces,*nf*3*sizeof(int32_t)); *faces = grown; *cap = next;
    }
    for (int k = 0; k < 3; k++) (*faces)[3*(*nf)+k] = (int32_t)tri[k];
    (*nf)++;
}

static int32_t *av_solid_faces(Arena_T arena, const AvClaim *cells, size_t nvert, double lim, size_t *out_nf)
{
    size_t nf = 0, cap = nvert*2+1;
    int32_t *lf = ARENA_ALLOC(arena,cap*3*sizeof(int32_t));
    uint8_t *partial = ARENA_CALLOC(arena,nvert ? nvert : 1,1);
    int left_present = 0;
    for (size_t i = 0; i < nvert; i++) {
        int64_t key = cells[i].key;
        if (i == 0 || cells[i-1].key != key) left_present = i > 0 && cells[i-1].key == key-1;
        int64_t right = key+1, down = key+((int64_t)1<<24), diag = down+1;
        size_t ir = (size_t)-1, id = (size_t)-1, ig = (size_t)-1;
        if (i+1 < nvert && cells[i+1].key == right) ir = i+1;
        size_t lo = av_lower(cells,0,nvert,down);
        if (lo < nvert && cells[lo].key == down) { id = lo; if (lo+1 < nvert && cells[lo+1].key == diag) ig = lo+1; }
        if (ir == SIZE_MAX || id == SIZE_MAX || ig == SIZE_MAX) partial[i] |= 1;
        if (!left_present && id != SIZE_MAX) partial[i] |= 2;
        if (ir != (size_t)-1 && id != (size_t)-1 && ig != (size_t)-1) {
            const size_t q[4] = {i,ir,id,ig};
            nf += av_solid_quad(cells,q,lim,lf+3*nf);
        }
    }
    /* A boundary can cover three corners of a cell. Its triangle has the
     * same physical edge checks as a complete cell; no sample is invented. */
    for (size_t i = 0; i < nvert; i++) for (int left = 0; left < 2; left++) {
        size_t tri[3];
        if (!(partial[i] & (1u<<left))) continue;
        if (!av_partial_triangle(cells,nvert,i,left,tri)) continue;
        int ok = 1;
        for (int e = 0; e < 3; e++) if (!(av_dist(cells[tri[e]].p,cells[tri[(e+1)%3]].p) < lim)) ok = 0;
        if (ok) av_append_triangle(arena,&lf,&nf,&cap,tri);
    }
    *out_nf = nf; return lf;
}

static int av_write_solid_ribbon(Arena_T arena, const AsmRun *run, const AsmVerdictOpts *o, const char *out_dir,
                                 const MeshBinData *sheet, const int32_t *sheet_chart, const float *sheet_nrm,
                                 const int32_t *lineage, const int64_t *relkey, size_t nrel, AsmVerdictStats *st)
{
    char path[2048];
    size_t n = 0, cap = 0;
    long row_min = 0, row_max = 0;
    AvClaim *cl = av_rasterize(arena, sheet, sheet_chart, sheet_nrm, run, o->du, o->dv, 0, &n, &cap, &row_min, &row_max);
    av_source_strips(arena,run,sheet,sheet_chart,sheet_nrm,o->du,o->dv,0,&cl,&n,&cap,&row_min,&row_max);
    if (n == 0) return -1;
    qsort(cl, n, sizeof(AvClaim), av_cmp);
    /* one vertex per cell; contested = another claimant a layer apart along the normal */
    size_t nvert = 0, strip_conflicts_reported = 0;
    AvClaim *cells = ARENA_ALLOC(arena, n * sizeof(AvClaim));
    double dl = run->layer_d_global > 0.0 ? run->layer_d_global : 10.0;
    for (size_t i = 0; i < n;) {
        size_t j = i;
        while (j < n && cl[j].key == cl[i].key) j++;
        cells[nvert] = cl[i];
        for (size_t k = i + 1; k < j; k++) {
            double dn = av_normal_distance(cl+i,cl+k);
            if (dn > 0.5 * dl) {
                st->solid_contested++;
                if (cl[k].area < 0 && strip_conflicts_reported++ < 16)
                    fprintf(stderr,"  [verdict] generated strip conflict: row %ld col %ld, chart %d at %.6f %.6f %.6f versus %d at %.6f %.6f %.6f, normal separation %.9g\n",
                            av_key_row(cl[i].key),av_key_col(cl[i].key),cl[i].chart,cl[i].p[0],cl[i].p[1],cl[i].p[2],
                            cl[k].chart,cl[k].p[0],cl[k].p[1],cl[k].p[2],dn);
                break;
            }
        }
        nvert++;
        i = j;
    }
    /* SEAM ROWS.  Charts of z-adjacent cubes never touch (the per-cube trim
     * leaves ~4 vox), so with rows on z-planes the cube faces at z = k*128 are
     * EMPTY rows, and the solid stage reads every empty row as an unjoinable
     * band boundary.  The seam JOINS are the measured evidence that the two
     * charts continue across the gap, so a gap of up to 3 rows between two
     * cells of joined charts is closed by interpolation, labelled generated
     * (support 0, provenance 1), never observed. */
    size_t n_bridged = 0;
    {
        size_t bcap = 1u << 14, nb = 0;
        AvClaim *br = ARENA_ALLOC(arena, bcap * sizeof(AvClaim));
        for (size_t i = 0; i < nvert; i++) {
            long r = av_key_row(cells[i].key), c = av_key_col(cells[i].key);
            /* is (r+1, c) empty? */
            int64_t k1 = av_key(r + 1, c);
            size_t lo = av_lower(cells, i + 1, nvert, k1);
            if (lo < nvert && cells[lo].key == k1) continue;
            size_t j = (size_t)-1;
            int gap = 0;
            for (int g = 2; g <= 4 && j == (size_t)-1; g++) {
                int64_t kg = av_key(r + g, c);
                size_t l2 = av_lower(cells, lo, nvert, kg);
                if (l2 < nvert && cells[l2].key == kg) { j = l2; gap = g; }
            }
            if (j == (size_t)-1) continue;
            if (cells[i].chart == cells[j].chart) continue;   /* a hole in one chart stays a hole */
            if (!av_joined(relkey, nrel, cells[i].chart, cells[j].chart)) continue;
            /* the two cells must be one sheet: the 3-D step per row under the wrap gate */
            if (av_dist(cells[j].p, cells[i].p) > 2.5 * (double)gap * o->dv) continue;
            for (int g = 1; g < gap; g++) {
                double t = (double)g / (double)gap;
                if (nb == bcap) { AvClaim *nc = ARENA_ALLOC(arena, bcap * 2 * sizeof(AvClaim)); memcpy(nc, br, nb * sizeof(AvClaim)); br = nc; bcap *= 2; }
                AvClaim *q = &br[nb++];
                q->key = av_key(r + g, c);
                q->chart = -1 - cells[i].chart;   /* negative marks generated; the lane is the lower chart */
                q->area = 0.0;
                q->p[0] = (float)((double)(r + g) * o->dv);
                for (int d = 1; d < 3; d++) q->p[d] = (float)((1.0 - t) * cells[i].p[d] + t * cells[j].p[d]);
                for (int d = 0; d < 3; d++) q->n[d] = (float)((1.0 - t) * cells[i].n[d] + t * cells[j].n[d]);
            }
        }
        if (nb > 0) {
            AvClaim *merged = ARENA_ALLOC(arena, (nvert + nb) * sizeof(AvClaim));
            memcpy(merged, cells, nvert * sizeof(AvClaim));
            memcpy(merged + nvert, br, nb * sizeof(AvClaim));
            qsort(merged, nvert + nb, sizeof(AvClaim), av_cmp);
            /* drop duplicate keys (two bridges meeting): keep the first */
            size_t w = 0;
            for (size_t i = 0; i < nvert + nb; i++) if (w == 0 || merged[i].key != merged[w-1].key) merged[w++] = merged[i];
            cells = merged;
            n_bridged = w - nvert;
            nvert = w;
        }
    }
    st->cells_bridged = n_bridged;
    for (size_t i = 0; i < nvert; i++) if (cells[i].area < 0) st->cells_bridged++;

    /* lattice vertices: u = col * du, v = (row - row_min) * dv, z = row * dv */
    float *lv = ARENA_ALLOC(arena, nvert * 3 * sizeof(float));
    float *luv = ARENA_ALLOC(arena, nvert * 2 * sizeof(float));
    int32_t *lane = ARENA_ALLOC(arena, nvert * sizeof(int32_t));
    int32_t *mat = ARENA_ALLOC(arena, nvert * sizeof(int32_t));
    uint8_t *support = ARENA_ALLOC(arena, nvert);
    uint8_t *prov = ARENA_ALLOC(arena, nvert);
    float *phase = ARENA_ALLOC(arena, nvert * sizeof(float));
    long col_min = 1L << 30, col_max = -(1L << 30);
    for (size_t i = 0; i < nvert; i++) {
        long r = av_key_row(cells[i].key), c = av_key_col(cells[i].key);
        memcpy(&lv[i*3], cells[i].p, 3 * sizeof(float));
        luv[i*2] = (float)((double)c * o->du); luv[i*2+1] = (float)((double)(r - row_min) * o->dv);
        int32_t chart = cells[i].chart >= 0 ? cells[i].chart : -1 - cells[i].chart;
        if (cells[i].chart >= 0 && cells[i].area >= 0) { lane[i] = chart; support[i] = 255; prov[i] = 0; }
        else { lane[i] = chart; support[i] = 0; prov[i] = 1; }   /* generated seam row */
        mat[i] = lineage[(size_t)chart];
        if (c < col_min) col_min = c;
        if (c > col_max) col_max = c;
    }
    av_lift_phase(cells, nvert, o->umb_y, o->umb_x, o->du, phase);
    /* Solid connections stay under the physical wrap gate; the verdict
     * lattice below is the consumer that exposes long claimed connections. */
    size_t nf;
    int32_t *lf = av_solid_faces(arena,cells,nvert,2.99*fmax(o->du,o->dv),&nf);
    st->solid_verts = nvert; st->solid_faces = nf;
    st->solid_rows = (size_t)(row_max - row_min + 1);
    st->solid_cols = (size_t)(col_max - col_min + 1);
    snprintf(path, sizeof path, "%s/assembled_ribbon.vmesh", out_dir);
    MeshBin_write(path, lv, nvert, lf, nf, luv);
    snprintf(path, sizeof path, "%s/assembled_ribbon_lane.i32", out_dir);
    av_write_raw(path, lane, sizeof(int32_t), nvert);
    snprintf(path, sizeof path, "%s/assembled_ribbon_material_identity.i32", out_dir);
    av_write_raw(path, mat, sizeof(int32_t), nvert);
    snprintf(path, sizeof path, "%s/assembled_ribbon_support.u8", out_dir);
    av_write_raw(path, support, 1, nvert);
    snprintf(path, sizeof path, "%s/assembled_ribbon_provenance.u8", out_dir);
    av_write_raw(path, prov, 1, nvert);
    snprintf(path, sizeof path, "%s/assembled_ribbon_phase.f32", out_dir);
    av_write_raw(path, phase, sizeof(float), nvert);
    snprintf(path, sizeof path, "%s/assembled_ribbon_stats.json", out_dir);
    {
        FILE *fp = fopen(path, "wb");
        if (fp) {
            fprintf(fp, "{\n  \"id\": \"assembled\",\n  \"n_verts\": %zu,\n  \"n_faces\": %zu,\n"
                        "  \"axis_point_zyx\": [0.0, %.4f, %.4f],\n  \"axis_dir_zyx\": [1.0, 0.0, 0.0],\n"
                        "  \"opts\": { \"slice_h\": %.3f, \"grid_u\": %.3f, \"ownership\": \"assembled\" },\n"
                        "  \"slicing\": { \"n_slices\": %zu, \"z_min\": %.3f, \"n_bands\": 1 },\n"
                        "  \"lattice\": { \"col_min\": %ld, \"col_max\": %ld, \"cells_contested\": %zu, \"cells_bridged\": %zu }\n}\n",
                    nvert, nf, o->umb_y, o->umb_x, o->dv, o->du, st->solid_rows, (double)row_min * o->dv,
                    col_min, col_max, st->solid_contested, st->cells_bridged);
            fclose(fp);
        }
    }
    return 0;
}

/* ---- 2. the verdict lattice (rows = v; what ribbon_verdict measures) ---------- */

/* the vertex of cell `key` belonging to `chart`, else the cell's primary
 * claimant, else -1 (cells sorted by key, the primary first within a key) */
static size_t av_pick(const AvClaim *cells, size_t nvert, int64_t key, int32_t chart)
{
    size_t lo = av_lower(cells, 0, nvert, key);
    if (lo >= nvert || cells[lo].key != key) return (size_t)-1;
    for (size_t m = lo; m < nvert && cells[m].key == key; m++) if (cells[m].chart == chart) return m;
    return lo;
}

static size_t av_pick_near(const AvClaim *cells, size_t nvert, int64_t key, const AvClaim *reference)
{
    size_t best = av_pick(cells,nvert,key,reference->chart);
    if (best == (size_t)-1 || cells[best].chart != reference->chart) return best;
    double distance = av_dist(reference->p,cells[best].p);
    for (size_t i = best+1; i < nvert && cells[i].key == key; i++) if (cells[i].chart == reference->chart) {
        double next = av_dist(reference->p,cells[i].p);
        if (next < distance) { best = i; distance = next; }
    }
    return best;
}

/* Visit a three-corner cell exactly once, in positive raster winding.
 * Ordinarily i is its top-left corner. If that corner is absent, its
 * top-right corner owns the triangle to the left instead. Complete cells
 * retain their existing triangulation; fewer than three samples stay open. */
static int av_partial_triangle(const AvClaim *cells, size_t nvert, size_t i, int left, size_t tri[3])
{
    const int64_t stride = (int64_t)1<<24;
    int64_t key = cells[i].key;
    size_t q[4];
    if (left) {
        size_t lo = av_lower(cells,0,i,key-1);
        if (lo < nvert && cells[lo].key == key-1) return 0;
        q[0] = SIZE_MAX; q[1] = i;
        q[2] = av_pick_near(cells,nvert,key+stride-1,cells+i);
        q[3] = av_pick_near(cells,nvert,key+stride,cells+i);
    } else {
        q[0] = i;
        q[1] = av_pick_near(cells,nvert,key+1,cells+i);
        q[2] = av_pick_near(cells,nvert,key+stride,cells+i);
        q[3] = av_pick_near(cells,nvert,key+stride+1,cells+i);
    }
    int n = 0;
    const int order[4] = {0,1,3,2};
    for (int k = 0; k < 4; k++) if (q[order[k]] != SIZE_MAX) {
        if (n == 3) return 0;
        tri[n++] = q[order[k]];
    }
    return n == 3;
}

typedef struct AvLong { int64_t pair; double step; } AvLong;
typedef struct AvLatticeFaces { int32_t *faces; size_t count; AvLong *longs; size_t nlong; } AvLatticeFaces;

enum { AV_EDGE_LINEAGE = -1, AV_EDGE_REFUSED, AV_EDGE_CHART, AV_EDGE_SOURCE, AV_EDGE_GEOMETRY };

static int av_lattice_edge(const AsmRun *run, const AvClaim *a, const AvClaim *b,
                           const int32_t *lineage, const int64_t *relkey, size_t nrel,
                           const AvSeamClaim *claims, size_t nclaims, double step, double lim)
{
    if (a->chart == b->chart) return AV_EDGE_CHART;
    if (av_joined(relkey,nrel,a->chart,b->chart) && av_seam_claimed_at(run,claims,nclaims,a,b,lim)) return AV_EDGE_SOURCE;
    /* Lineages are joined by the stitcher, never by raster coincidence. */
    if (lineage[(size_t)a->chart] != lineage[(size_t)b->chart]) return AV_EDGE_LINEAGE;
    return step < lim ? AV_EDGE_GEOMETRY : AV_EDGE_REFUSED;
}

static AvLatticeFaces av_lattice_faces(Arena_T arena, const AsmRun *run,
                                       const AvClaim *cells, size_t nvert, const int32_t *lineage,
                                       const int64_t *relkey, size_t nrel,
                                       const AvSeamClaim *claims, size_t nclaims, double lim, AsmVerdictStats *st)
{
    /* faces: for every vertex, the right / down / diagonal cells' vertex of the
     * same chart, else that cell's primary.  A quad is emitted when every one
     * of its five edges is a CLAIM of continuity (one chart, or an accepted
     * seam join -- whatever the 3-D step, so a wrong join reaches TANGLE) or a
     * geometric continuation of unrelated charts under the wrap gate. */
    size_t fcap = nvert * 2 + 1, nf = 0;
    int32_t *lf = ARENA_ALLOC(arena, fcap * 3 * sizeof(int32_t));
    /* the long-join ledger: relation edges over the gate, per chart pair */
    size_t lcap = 1u << 10, nlong = 0;
    AvLong *longs = ARENA_ALLOC(arena, lcap * sizeof(AvLong));
    /* Remember boundary candidates during the existing lookups instead of
     * repeating them for every dense interior cell in the partial pass. */
    uint8_t *partial = ARENA_CALLOC(arena,nvert ? nvert : 1,1);
    int left_present = 0;
    for (size_t i = 0; i < nvert; i++) {
        int64_t key = cells[i].key;
        if (i == 0 || cells[i-1].key != key) left_present = i > 0 && cells[i-1].key == key-1;
        int32_t chart = cells[i].chart;
        size_t ir = av_pick_near(cells, nvert, key + 1, cells+i);
        size_t id = av_pick_near(cells, nvert, key + ((int64_t)1 << 24), cells+i);
        size_t ig = av_pick_near(cells, nvert, key + ((int64_t)1 << 24) + 1, cells+i);
        if (ir == SIZE_MAX || id == SIZE_MAX || ig == SIZE_MAX) partial[i] |= 1;
        if (!left_present && id != SIZE_MAX) partial[i] |= 2;
        if (ir == (size_t)-1 || id == (size_t)-1 || ig == (size_t)-1) continue;
        const size_t q[4] = { i, ir, id, ig };
        int pairs[5][2] = { {0, 1}, {0, 2}, {1, 3}, {2, 3}, {0, 3} };
        int flip = 0;
        if (chart != cells[ir].chart || chart != cells[id].chart || chart != cells[ig].chart) {
            /* A fixed raster diagonal can create a long seam edge even
             * when all four sides are sound. Choose the shorter physical
             * diagonal if it has the same continuity authorization required
             * of every emitted edge. Keep all corners and both triangles;
             * excessive perimeter steps still reach the original audits. */
            double other = av_dist(cells[ir].p,cells[id].p);
            if (other < av_dist(cells[i].p,cells[ig].p) &&
                av_lattice_edge(run,&cells[ir],&cells[id],lineage,relkey,nrel,claims,nclaims,other,lim) > 0) {
                pairs[4][0] = 1; pairs[4][1] = 2; flip = 1;
            }
        }
        AvLong quad_longs[5]; size_t n_quad_longs = 0;
        int ok = 1, by_geom = 0, cross = 0, cross_lineage = 0;
        for (int e = 0; e < 5 && ok; e++) {
            const AvClaim *a = &cells[q[pairs[e][0]]], *b = &cells[q[pairs[e][1]]];
            if (a->chart == b->chart) continue;
            cross = 1;
            double step = av_dist(a->p, b->p);
            int kind = av_lattice_edge(run,a,b,lineage,relkey,nrel,claims,nclaims,step,lim);
            if (kind == AV_EDGE_SOURCE) {
                if (step >= lim) {
                    int32_t x = a->chart < b->chart ? a->chart : b->chart, y = a->chart < b->chart ? b->chart : a->chart;
                    quad_longs[n_quad_longs++] = (AvLong){((int64_t)x << 32) | (int64_t)y,step};
                }
                continue;
            }
            if (kind == AV_EDGE_GEOMETRY) { by_geom = 1; continue; }
            ok = 0; cross_lineage = kind == AV_EDGE_LINEAGE;
        }
        if (!ok) { if (cross_lineage) st->quads_refused_lineage++; else st->quads_refused++; continue; }
        /* Only emitted edges belong in the lattice ledger. A later edge
         * may reject the whole quad after an earlier long seam passed. */
        for (size_t e = 0; e < n_quad_longs; e++) {
            if (nlong == lcap) { AvLong *nl = ARENA_ALLOC(arena, lcap * 2 * sizeof(AvLong)); memcpy(nl, longs, nlong * sizeof(AvLong)); longs = nl; lcap *= 2; }
            longs[nlong++] = quad_longs[e];
            if (quad_longs[e].step <= 8.0) st->long_join_edges_tail++;
        }
        lf[nf*3] = (int32_t)i; lf[nf*3+1] = (int32_t)ir; lf[nf*3+2] = (int32_t)(flip ? id : ig); nf++;
        lf[nf*3] = (int32_t)(flip ? ir : i); lf[nf*3+1] = (int32_t)ig; lf[nf*3+2] = (int32_t)id; nf++;
        if (by_geom) st->faces_by_geometry += 2;
        else if (cross) st->faces_by_relation += 2;
    }
    for (size_t i = 0; i < nvert; i++) for (int left = 0; left < 2; left++) {
        size_t tri[3];
        if (!(partial[i] & (1u<<left))) continue;
        if (!av_partial_triangle(cells,nvert,i,left,tri)) continue;
        AvLong triangle_longs[3]; size_t n_triangle_longs = 0;
        int ok = 1, by_geom = 0, cross = 0;
        for (int e = 0; e < 3 && ok; e++) {
            const AvClaim *a = cells+tri[e], *b = cells+tri[(e+1)%3];
            if (a->chart == b->chart) continue;
            cross = 1;
            double step = av_dist(a->p,b->p);
            int kind = av_lattice_edge(run,a,b,lineage,relkey,nrel,claims,nclaims,step,lim);
            if (kind == AV_EDGE_SOURCE) {
                if (step >= lim) {
                    int32_t x = a->chart < b->chart ? a->chart : b->chart, y = a->chart < b->chart ? b->chart : a->chart;
                    triangle_longs[n_triangle_longs++] = (AvLong){((int64_t)x<<32)|(int64_t)y,step};
                }
            } else if (kind == AV_EDGE_GEOMETRY) by_geom = 1;
            else ok = 0;
        }
        if (!ok) continue;
        av_append_triangle(arena,&lf,&nf,&fcap,tri);
        for (size_t e = 0; e < n_triangle_longs; e++) {
            if (nlong == lcap) { AvLong *nl = ARENA_ALLOC(arena,lcap*2*sizeof(AvLong)); memcpy(nl,longs,nlong*sizeof(AvLong)); longs = nl; lcap *= 2; }
            longs[nlong++] = triangle_longs[e];
            if (triangle_longs[e].step <= 8.0) st->long_join_edges_tail++;
        }
        if (by_geom) st->faces_by_geometry++;
        else if (cross) st->faces_by_relation++;
    }
    return (AvLatticeFaces){lf,nf,longs,nlong};
}

static AvClaim *av_lattice_cells(Arena_T arena, const AsmRun *run, const AsmVerdictOpts *o,
                                  const MeshBinData *sheet, const int32_t *sheet_chart, const float *sheet_nrm,
                                  const int32_t *lineage, size_t *out_n, long *out_row_min, long *out_row_max,
                                  AsmVerdictStats *st)
{
    size_t n = 0, cap = 0;
    long row_min = 0, row_max = 0;
    AvClaim *cl = av_rasterize(arena, sheet, sheet_chart, sheet_nrm, run, o->du, o->dv, 1, &n, &cap, &row_min, &row_max);
    av_source_strips(arena,run,sheet,sheet_chart,sheet_nrm,o->du,o->dv,1,&cl,&n,&cap,&row_min,&row_max);
    if (n == 0) return NULL;
    qsort(cl, n, sizeof(AvClaim), av_cmp);
    double dl = run->layer_d_global > 0.0 ? run->layer_d_global : 10.0;
    /* Separate chart overlaps use the existing normal-distance policy.
     * Within one source chart, distinct physical samples at the same UV
     * are a self-overlap even when the displacement is tangential or less
     * than half a layer. Only coincident triangle samples are duplicates. */
    AvClaim *cells = ARENA_ALLOC(arena, n * sizeof(AvClaim));
    size_t nvert = 0;
    for (size_t i = 0; i < n;) {
        size_t j = i, first = nvert, chart_first = i;
        while (j < n && cl[j].key == cl[i].key) j++;
        cells[nvert++] = cl[i];
        int contested = 0;
        for (size_t k = i + 1; k < j; k++) {
            /* Original claims of a chart have one area and are contiguous
             * in this ordering. Remember even a first sample that folded
             * into another chart, so that chart's own second sheet cannot
             * disappear merely because it was not the cell's primary. */
            if (cl[k].chart != cl[chart_first].chart || cl[k].area != cl[chart_first].area) chart_first = k;
            double tolerance = av_sample_tolerance(cl[k].p,cl[chart_first].p,fmin(o->du,o->dv));
            int self_overlap = cl[k].area >= 0 && av_dist(cl[k].p,cl[chart_first].p) > tolerance;
            double dn = av_normal_distance(cl+i,cl+k);
            if (dn <= 0.5 * dl && !self_overlap) continue;
            int seen = 0;
            for (size_t m = first; m < nvert; m++) if (cells[m].chart == cl[k].chart) {
                double limit = cells[m].area >= 0 && cl[k].area >= 0
                    ? av_sample_tolerance(cells[m].p,cl[k].p,fmin(o->du,o->dv)) : .5*dl;
                if (av_dist(cells[m].p,cl[k].p) <= limit) { seen = 1; break; }
            }
            if (seen) continue;
            cells[nvert++] = cl[k];
            st->stacked_verts++;
            contested = 1;
        }
        st->cells_contested += (size_t)contested;
        i = j;
    }

    for (size_t i = 0; i < nvert; i++) if (cells[i].area < 0) st->cells_bridged_uv++;
    /* BRIDGES: a gap of up to ASM_VERDICT_BRIDGE_CELLS empty cells between a cell and the next
     * occupied cell in u or in v, when the two cells belong to two JOINED charts and the 3-D
     * step is a sheet's (under 2.5 du per cell), is closed by generated cells carrying the
     * first cell's chart (support 0).  A gap inside one chart stays a hole. */
    if (ASM_VERDICT_BRIDGE_CELLS > 0) {
        size_t bcap = 1u << 12, nb = 0;
        AvClaim *br = ARENA_ALLOC(arena, bcap * sizeof(AvClaim));
        for (size_t i = 0; i < nvert; i++) {
            long r = av_key_row(cells[i].key), c = av_key_col(cells[i].key);
            for (int dir = 0; dir < 2; dir++) {
                int64_t next = dir == 0 ? av_key(r, c + 1) : av_key(r + 1, c);
                size_t lo = av_lower(cells, i + 1, nvert, next);
                if (lo < nvert && cells[lo].key == next) continue;   /* occupied: no gap */
                size_t j = (size_t)-1; int gap = 0;
                for (int g = 2; g <= ASM_VERDICT_BRIDGE_CELLS + 1 && j == (size_t)-1; g++) {
                    int64_t kg = dir == 0 ? av_key(r, c + g) : av_key(r + g, c);
                    size_t l2 = av_lower(cells, lo, nvert, kg);
                    if (l2 < nvert && cells[l2].key == kg) { j = l2; gap = g; }
                }
                if (j == (size_t)-1) continue;
                /* the far cell's claimant of the same chart if any, else its primary */
                size_t jj = av_pick_near(cells, nvert, cells[j].key, cells+i);
                if (jj != (size_t)-1) j = jj;
                if (cells[i].chart == cells[j].chart) continue;
                /* one LINEAGE (joined directly or through other charts), as the face rule */
                if (lineage[(size_t)cells[i].chart] != lineage[(size_t)cells[j].chart]) continue;
                if (av_dist(cells[j].p, cells[i].p) > 2.5 * (double)gap * fmax(o->du, o->dv)) continue;
                for (int g = 1; g < gap; g++) {
                    double t = (double)g / (double)gap;
                    if (nb == bcap) { AvClaim *nc = ARENA_ALLOC(arena, bcap * 2 * sizeof(AvClaim)); memcpy(nc, br, nb * sizeof(AvClaim)); br = nc; bcap *= 2; }
                    AvClaim *q = &br[nb++];
                    q->key = dir == 0 ? av_key(r, c + g) : av_key(r + g, c);
                    q->chart = cells[i].chart;
                    q->area = -2.0;   /* new row/column bridge, below SOURCE strips */
                    for (int d = 0; d < 3; d++) { q->p[d] = (float)((1.0 - t) * cells[i].p[d] + t * cells[j].p[d]); q->n[d] = (float)((1.0 - t) * cells[i].n[d] + t * cells[j].n[d]); }
                }
            }
        }
        if (nb > 0) {
            AvClaim *merged = ARENA_ALLOC(arena, (nvert + nb) * sizeof(AvClaim));
            memcpy(merged, cells, nvert * sizeof(AvClaim));
            memcpy(merged + nvert, br, nb * sizeof(AvClaim));
            qsort(merged, nvert + nb, sizeof(AvClaim), av_cmp);   /* generated (area -1) sort after any observed claimant */
            size_t w = 0;
            for (size_t i = 0; i < nvert + nb; i++) {
                /* Existing cells include retained SOURCE-strip layers.
                 * Only a newly proposed axis bridge loses to an occupied
                 * cell; never erase an already selected physical layer. */
                if (merged[i].area == -2.0 && w > 0 && merged[i].key == merged[w-1].key) continue;
                merged[w++] = merged[i];
            }
            cells = merged;
            st->cells_bridged_uv += w - nvert;
            nvert = w;
        }
    }
    *out_n = nvert; *out_row_min = row_min; *out_row_max = row_max;
    return cells;
}

static int av_write_verdict_lattice(Arena_T arena, const AsmRun *run, const AsmVerdictOpts *o, const char *out_dir,
                                    const MeshBinData *sheet, const int32_t *sheet_chart, const float *sheet_nrm,
                                    const int32_t *lineage, const int64_t *relkey, size_t nrel,
                                    const AsmAxis *axis, AsmVerdictStats *st)
{
    char path[2048];
    size_t nvert = 0; long row_min = 0, row_max = 0;
    AvClaim *cells = av_lattice_cells(arena,run,o,sheet,sheet_chart,sheet_nrm,lineage,&nvert,&row_min,&row_max,st);
    if (!cells) return -1;
    size_t nclaims; AvSeamClaim *claims = av_seam_claims(arena,run,&nclaims);
    float *lv = ARENA_ALLOC(arena, nvert * 3 * sizeof(float));
    float *luv = ARENA_ALLOC(arena, nvert * 2 * sizeof(float));
    int32_t *lane = ARENA_ALLOC(arena, nvert * sizeof(int32_t));
    int32_t *mat = ARENA_ALLOC(arena, nvert * sizeof(int32_t));
    uint8_t *support = ARENA_ALLOC(arena, nvert);
    uint8_t *lin_seen = ARENA_CALLOC(arena, run->n_charts ? run->n_charts : 1, 1);
    long col_min = 1L << 30, col_max = -(1L << 30);
    for (size_t i = 0; i < nvert; i++) {
        long r = av_key_row(cells[i].key), c = av_key_col(cells[i].key);
        memcpy(&lv[i*3], cells[i].p, 3 * sizeof(float));
        luv[i*2] = (float)((double)c * o->du); luv[i*2+1] = (float)((double)(r - row_min) * o->dv);
        lane[i] = cells[i].chart;
        support[i] = cells[i].area < 0.0 ? 0 : 255;
        mat[i] = lineage[(size_t)cells[i].chart];
        if (!lin_seen[(size_t)mat[i]]) { lin_seen[(size_t)mat[i]] = 1; st->lineages++; }
        if (c < col_min) col_min = c;
        if (c > col_max) col_max = c;
    }
    AvLatticeFaces emitted = av_lattice_faces(arena,run,cells,nvert,lineage,relkey,nrel,claims,nclaims,
                                              2.99*fmax(o->du,o->dv),st);
    int32_t *lf = emitted.faces; size_t nf = emitted.count;
    AvLong *longs = emitted.longs; size_t nlong = emitted.nlong;
    const double lim = 2.99*fmax(o->du,o->dv);
    st->lattice_verts = nvert; st->lattice_faces = nf;
    {
        char audit_path[2048];
        snprintf(audit_path, sizeof audit_path, "%s/continuity.csv", out_dir);
        int32_t primary = -1;
        for (size_t i = 0; i < run->n_charts; i++) if (run->charts[i].placement_state == ASM_PLACE_ROOT) { primary = run->charts[i].component; break; }
        AsmContinuity_audit(run, primary, lane, lv, nvert, lf, nf, audit_path, &st->continuity);
        fprintf(stderr, "[assemble continuity] emitted geometry: unsupported %zu islands / %zu charts (%.3e area), source breaks %zu, metadata-only joins %zu, wrong placements %zu, unresolved %zu; confetti_free=%d\n",
                st->continuity.unsupported_islands, st->continuity.unsupported_charts, st->continuity.unsupported_area,
                st->continuity.unexplained_breaks, st->continuity.metadata_only_joins, st->continuity.wrong_placements,
                st->continuity.unresolved_hypotheses, st->continuity.confetti_free);
    }
    if (axis != NULL) for (size_t i = 0; i < nvert; i++) av_to_frame(axis,&lv[i*3]);
    st->lattice_rows = (size_t)(row_max - row_min + 1);
    st->lattice_cols = (size_t)(col_max - col_min + 1);
    /* GAPS: every empty run between two occupied cells of a row, split by the radius of the two
     * bounding cells (the frame's (y, x) about the axis when there is one, else about the umbilicus) */
    {
        double dl = run->layer_d_global > 0.0 ? run->layer_d_global : 10.0;
        size_t cap_g = 1u << 16, nw = 0, nb2 = 0;
        double *gw = ARENA_ALLOC(arena, cap_g * sizeof(double)), *gb = ARENA_ALLOC(arena, cap_g * sizeof(double));
        st->gap_measured = axis != NULL || (o->umb_y > 0.0 && o->umb_x > 0.0);
        for (size_t i = 0; i + 1 < nvert && st->gap_measured; i++) {
            if (av_key_row(cells[i].key) != av_key_row(cells[i+1].key)) continue;
            long gapc = av_key_col(cells[i+1].key) - av_key_col(cells[i].key) - 1;
            if (gapc <= 0) continue;
            double gap = (double)gapc * o->du;
            double ri, rj;
            if (axis != NULL) { ri = hypot((double)lv[i*3+1], (double)lv[i*3+2]); rj = hypot((double)lv[(i+1)*3+1], (double)lv[(i+1)*3+2]); }
            else { ri = hypot((double)cells[i].p[1] - o->umb_y, (double)cells[i].p[2] - o->umb_x); rj = hypot((double)cells[i+1].p[1] - o->umb_y, (double)cells[i+1].p[2] - o->umb_x); }
            double d3 = av_dist(cells[i].p, cells[i+1].p);
            if (fabs(ri - rj) < 0.5 * dl && d3 <= 1.5 * gap + 8.0) {
                if (nw < cap_g) gw[nw] = gap;
                nw++; st->gap_within_vox += gap;
            } else if (fabs(ri - rj) >= 0.5 * dl) {
                if (nb2 < cap_g) gb[nb2] = gap;
                nb2++; st->gap_between_vox += gap;
            }
        }
        st->gap_within_n = nw; st->gap_between_n = nb2;
        size_t mw = nw < cap_g ? nw : cap_g, mb = nb2 < cap_g ? nb2 : cap_g;
        if (mw > 0) { qsort(gw, mw, sizeof(double), av_cmp_double); st->gap_within_p50 = gw[mw / 2]; st->gap_within_p90 = gw[(size_t)((double)(mw - 1) * 0.9)]; }
        if (mb > 0) { qsort(gb, mb, sizeof(double), av_cmp_double); st->gap_between_p50 = gb[mb / 2]; st->gap_between_p90 = gb[(size_t)((double)(mb - 1) * 0.9)]; }
    }
    /* the long-join ledger: which JOINS carry lattice edges over the wrap gate, with the
     * relation's own evidence beside them (TANGLE fails name these) */
    st->long_join_edges = nlong;
    {
        for (size_t x = 1; x < nlong; x++) { AvLong v = longs[x]; size_t y = x; while (y > 0 && (longs[y-1].pair > v.pair || (longs[y-1].pair == v.pair && longs[y-1].step < v.step))) { longs[y] = longs[y-1]; y--; } longs[y] = v; }
        snprintf(path, sizeof path, "%s/verdict_lattice_long_joins.csv", out_dir);
        FILE *fp = fopen(path, "wb");
        if (fp) fprintf(fp, "chart_a,chart_b,edges_over_gate,max_step,rel_flags,rel_rms,rel_residual,rel_robust_w,rel_seam_len,rel_normal_gap,rel_n_corr\n");
        size_t printed = 0;
        for (size_t i = 0; i < nlong;) {
            size_t j = i;
            while (j < nlong && longs[j].pair == longs[i].pair) j++;
            int32_t a = (int32_t)(longs[i].pair >> 32), b = (int32_t)(longs[i].pair & 0xFFFFFFFF);
            const AsmRelation *R = NULL;
            for (size_t r = 0; r < run->n_rels; r++) {
                const AsmRelation *T = &run->rels[r];
                if (!av_rel_is_join(T)) continue;
                if ((T->a == a && T->b == b) || (T->a == b && T->b == a)) { R = T; break; }
            }
            st->long_join_pairs++;
            if (fp) fprintf(fp, "%d,%d,%zu,%.2f,%u,%.3f,%.3f,%.3f,%.1f,%.3f,%zu\n", a, b, j - i, longs[i].step,
                            R ? R->flags : 0u, R ? R->rms : 0.0, R ? R->residual : 0.0, R ? R->robust_w : 0.0,
                            R ? R->seam_len : 0.0, R ? R->normal_gap : 0.0, R ? R->n_corr : (size_t)0);
            if (printed < 6) {
                fprintf(stderr, "  [verdict] long join %d-%d: %zu lattice edges over %.2f vox (max %.2f) | relation flags %u rms %.2f residual %.2f w %.2f seam %.0f n %zu\n",
                        a, b, j - i, lim, longs[i].step, R ? R->flags : 0u, R ? R->rms : 0.0, R ? R->residual : 0.0, R ? R->robust_w : 0.0,
                        R ? R->seam_len : 0.0, R ? R->n_corr : (size_t)0);
                printed++;
            }
            i = j;
        }
        if (fp) fclose(fp);
    }
    snprintf(path, sizeof path, "%s/verdict_lattice.vmesh", out_dir);
    MeshBin_write(path, lv, nvert, lf, nf, luv);
    snprintf(path, sizeof path, "%s/verdict_lattice_lane.i32", out_dir);
    av_write_raw(path, lane, sizeof(int32_t), nvert);
    snprintf(path, sizeof path, "%s/verdict_lattice_material_identity.i32", out_dir);
    av_write_raw(path, mat, sizeof(int32_t), nvert);
    snprintf(path, sizeof path, "%s/verdict_lattice_support.u8", out_dir);
    av_write_raw(path, support, 1, nvert);
    snprintf(path, sizeof path, "%s/verdict_lattice_stats.json", out_dir);
    {
        FILE *fp = fopen(path, "wb");
        if (fp) {
            fprintf(fp, "{\n  \"id\": \"verdict_lattice\",\n  \"rows\": \"sheet v\",\n  \"frame\": \"%s\",\n"
                        "  \"n_verts\": %zu,\n  \"n_faces\": %zu,\n  \"rows_n\": %zu,\n  \"cols_n\": %zu,\n"
                        "  \"cells_contested\": %zu,\n  \"stacked_verts\": %zu,\n  \"lineages\": %zu,\n"
                        "  \"faces_by_relation\": %zu,\n  \"faces_by_geometry\": %zu,\n  \"quads_refused\": %zu,\n  \"quads_refused_lineage\": %zu,\n  \"long_join_edges\": %zu,\n  \"long_join_pairs\": %zu,\n  \"cells_bridged\": %zu,\n"
                        "  \"axis_rows\": %zu\n}\n",
                    axis != NULL ? (axis->is_line ? "axis line" : "axis table") : "world",
                    nvert, nf, st->lattice_rows, st->lattice_cols,
                    st->cells_contested, st->stacked_verts, st->lineages,
                    st->faces_by_relation, st->faces_by_geometry, st->quads_refused, st->quads_refused_lineage, st->long_join_edges, st->long_join_pairs, st->cells_bridged_uv,
                    axis != NULL ? axis->n_rows : (size_t)0);
            fclose(fp);
        }
    }
    return 0;
}

/* ---- sources ------------------------------------------------------------------ */

/* source.vmesh (world, the --solidify contract) and verdict_source.vmesh (the
 * verdict's frame): every chart that is papyrus (blobs excluded unless asked) */
static void av_write_sources(Arena_T arena, const AsmRun *run, const AsmVerdictOpts *o, const char *out_dir,
                             const AsmAxis *axis, AsmVerdictStats *st)
{
    char path[2048];
    size_t snv = 0, snf = 0;
    for (size_t i = 0; i < run->n_charts; i++) {
        const AsmChart *c = &run->charts[i];
        if ((c->flags & ASM_CHART_BLOB) && !o->include_blobs_in_source) continue;
        snv += c->nv; snf += c->nf;
    }
    float *sv = ARENA_ALLOC(arena, (snv ? snv : 1) * 3 * sizeof(float));
    int32_t *sf = ARENA_ALLOC(arena, (snf ? snf : 1) * 3 * sizeof(int32_t));
    size_t pv = 0, pf = 0;
    for (size_t i = 0; i < run->n_charts; i++) {
        const AsmChart *c = &run->charts[i];
        if ((c->flags & ASM_CHART_BLOB) && !o->include_blobs_in_source) continue;
        memcpy(&sv[pv*3], c->xyz, c->nv * 3 * sizeof(float));
        for (size_t f = 0; f < c->nf; f++) for (int e = 0; e < 3; e++) sf[(pf+f)*3 + (size_t)e] = c->faces[f*3 + (size_t)e] + (int32_t)pv;
        pv += c->nv; pf += c->nf;
    }
    snprintf(path, sizeof path, "%s/source.vmesh", out_dir);
    MeshBin_write(path, sv, snv, sf, snf, NULL);
    if (axis != NULL) for (size_t i = 0; i < snv; i++) av_to_frame(axis, &sv[i*3]);
    snprintf(path, sizeof path, "%s/verdict_source.vmesh", out_dir);
    MeshBin_write(path, sv, snv, sf, snf, NULL);
    st->source_verts = snv; st->source_faces = snf;
}

/* ---- entry -------------------------------------------------------------------- */

static int av_copy_sheet_normals(const AsmRun *run, const MeshBinData *sheet,
                                  const int32_t *sheet_chart, float *sheet_nrm)
{
    /* Validate the whole export before the first write. Emission can fail
     * in a reused directory, leaving an older sheet and sidecar behind. A
     * final count check after memcpy is too late to protect this buffer. */
    if (sheet->nv && (!sheet_chart || !sheet_nrm)) return 0;
    size_t pv = 0;
    for (size_t i = 0; i < run->n_charts; i++) {
        const AsmChart *c = &run->charts[i];
        if (!AsmChart_registered(c)) continue;
        if (c->nv > sheet->nv-pv || (c->nv && !c->nrm) || c->id < 0 || (size_t)c->id != i) return 0;
        for (size_t v = 0; v < c->nv; v++) if (sheet_chart[pv+v] != c->id) return 0;
        pv += c->nv;
    }
    if (pv != sheet->nv) return 0;
    pv = 0;
    for (size_t i = 0; i < run->n_charts; i++) {
        const AsmChart *c = &run->charts[i];
        if (!AsmChart_registered(c)) continue;
        if (c->nv) memcpy(sheet_nrm+pv*3,c->nrm,c->nv*3*sizeof(float));
        pv += c->nv;
    }
    return 1;
}

const char *AsmVerdict_gate_word(const char *line, int *pass, int *measured)
{
    const char *w;
    *pass = 0; *measured = 1;
    if ((w = strstr(line, "PASS")) != NULL) { *pass = 1; return w + 4; }
    if ((w = strstr(line, "FAIL")) != NULL) return w + 4;
    *measured = 0;
    if ((w = strstr(line, "ABSENT")) != NULL) return w + 6;
    return NULL;
}

static int av_gate_word_selftest(void)
{
    int fails = 0, pass = -1, measured = -1; size_t n = 0;
    const char *v = AsmVerdict_gate_word("  TANGLE/cross-wrap        PASS  0 emitted edges step", &pass, &measured);
    if (!v || !pass || !measured || sscanf(v, " %zu", &n) != 1 || n != 0) fails++;
    v = AsmVerdict_gate_word("  TANGLE/cross-wrap        FAIL  187 emitted edges step", &pass, &measured);
    if (!v || pass || !measured || sscanf(v, " %zu", &n) != 1 || n != 187) fails++;
    v = AsmVerdict_gate_word("  TANGLE/cross-wrap        ABSENT  NOT MEASURED: pass --umb-y/--umb-x", &pass, &measured);
    if (!v || pass || measured || sscanf(v, " %zu", &n) == 1) fails++;
    v = AsmVerdict_gate_word("  support                        12 measured vertices", &pass, &measured);
    if (v || pass || measured) fails++;
    fprintf(stderr, "  verdict gate words PASS / FAIL / ABSENT: %s (%d failures)\n", fails ? "FAIL" : "ok", fails);
    return fails;
}

int AsmVerdict_run(AsmRun *run, const AsmVerdictOpts *o, const char *out_dir, AsmVerdictStats *st)
{
    memset(st, 0, sizeof *st);
    double t0 = ves_clock_sec();
    Arena_T arena = run->arena;
    Arena_Mark mark = Arena_save(arena);
    /* the emitted sheet carries the gauge: read it back */
    char path[2048];
    snprintf(path, sizeof path, "%s/sheet.vmesh", out_dir);
    MeshBinData sheet;
    if (MeshBin_read_precise_arena(arena,path,&sheet) != 0 || (sheet.uv == NULL && sheet.uv64 == NULL)) { Arena_restore(arena, mark); return -1; }
    snprintf(path, sizeof path, "%s/sheet_chart.i32", out_dir);
    int32_t *sheet_chart = ARENA_ALLOC(arena, sheet.nv * sizeof(int32_t));
    {
        FILE *fp = fopen(path, "rb");
        if (!fp || fread(sheet_chart, sizeof(int32_t), sheet.nv, fp) != sheet.nv) { if (fp) fclose(fp); Arena_restore(arena, mark); return -1; }
        fclose(fp);
    }
    /* per-vertex normals of the sheet from the charts */
    float *sheet_nrm = ARENA_ALLOC(arena, sheet.nv * 3 * sizeof(float));
    if (!av_copy_sheet_normals(run,&sheet,sheet_chart,sheet_nrm)) {
        fprintf(stderr,"[assemble verdict] sheet/chart sidecars or normals do not match the current charts\n");
        Arena_restore(arena, mark); return -1;
    }

    size_t nrel = 0;
    int64_t *relkey = av_join_keys(arena, run, &nrel);
    int32_t *lineage = av_lineages(arena, run);
    const AsmAxis *axis = o->axis;
    if (axis != NULL) st->axis_frame = 1;

    if (av_write_solid_ribbon(arena, run, o, out_dir, &sheet, sheet_chart, sheet_nrm, lineage, relkey, nrel, st) != 0) { Arena_restore(arena, mark); return -1; }
    if (av_write_verdict_lattice(arena, run, o, out_dir, &sheet, sheet_chart, sheet_nrm, lineage, relkey, nrel, axis, st) != 0) { Arena_restore(arena, mark); return -1; }
    av_write_sources(arena, run, o, out_dir, axis, st);

    /* spawn the verdict on the verdict lattice, in its frame */
    {
        char lattice[2048], lanep[2048], matp[2048], supp[2048], src[2048], report[2048], log[2048];
        char uy[64], ux[64], pitch[64], core[64], dus[64], dvs[64];
        snprintf(lattice, sizeof lattice, "%s/verdict_lattice.vmesh", out_dir);
        snprintf(lanep, sizeof lanep, "%s/verdict_lattice_lane.i32", out_dir);
        snprintf(matp, sizeof matp, "%s/verdict_lattice_material_identity.i32", out_dir);
        snprintf(supp, sizeof supp, "%s/verdict_lattice_support.u8", out_dir);
        snprintf(src, sizeof src, "%s/verdict_source.vmesh", out_dir);
        snprintf(report, sizeof report, "%s/verdict.json", out_dir);
        snprintf(log, sizeof log, "%s/verdict.log", out_dir);
        /* The radial frame: the lattice is already in the axis frame (origin on the axis) when a
         * table exists, else the configured umbilicus.  With NEITHER the umbilicus flags are not
         * passed at all, so the verdict reports its cross-wrap and depth-seam gates ABSENT (a
         * failure, honestly) instead of measuring radii about the world origin (2026-09-17). */
        st->frame_measured = axis != NULL || (o->umb_y > 0.0 && o->umb_x > 0.0);
        double uyv = axis != NULL ? 0.0 : o->umb_y, uxv = axis != NULL ? 0.0 : o->umb_x;
        snprintf(uy, sizeof uy, "%.3f", uyv); snprintf(ux, sizeof ux, "%.3f", uxv);
        snprintf(pitch, sizeof pitch, "%.3f", o->pitch); snprintf(core, sizeof core, "%.1f", o->core_radius);
        snprintf(dus, sizeof dus, "%.3f", o->du); snprintf(dvs, sizeof dvs, "%.3f", o->dv);
        const char *argv[32]; size_t na = 0;
        argv[na++] = o->verdict_exe; argv[na++] = lattice; argv[na++] = "--labels"; argv[na++] = lanep;
        argv[na++] = "--materials"; argv[na++] = matp; argv[na++] = "--support"; argv[na++] = supp;
        if (st->frame_measured) { argv[na++] = "--umb-y"; argv[na++] = uy; argv[na++] = "--umb-x"; argv[na++] = ux; }
        else fprintf(stderr, "[assemble verdict] no axis table and no umbilicus: the cross-wrap and depth-seam gates are ABSENT (counted as failures), not measured about the world origin\n");
        argv[na++] = "--pitch"; argv[na++] = pitch; argv[na++] = "--core-radius"; argv[na++] = core;
        argv[na++] = "--du"; argv[na++] = dus; argv[na++] = "--dv"; argv[na++] = dvs;
        argv[na++] = "--source"; argv[na++] = src; argv[na++] = "--report"; argv[na++] = report; argv[na] = NULL;
        remove(log);
        st->verdict_rc = ves_run_subprocess_logged(o->verdict_exe, argv, 0.0, log);
    }
    st->sec = ves_clock_sec() - t0;
    Arena_restore(arena, mark);
    return 0;
}

/* ---- selftest ----------------------------------------------------------------- */

static int av_partial_cell_selftest(void)
{
    int fails = 0;
    /* Each input is an entire planar triangle, exactly half of a raster
     * cell. No fourth sample exists or is needed to retain its surface. */
    for (int absent = 0; absent < 4; absent++) for (int mirror = 0; mirror < 2; mirror++) {
        Arena_T arena = Arena_new(); AsmRun run = {0}; run.arena = arena;
        float xyz[9], uv[6], normals[9] = {0}; int32_t labels[3] = {0}, face[3] = {0,1,2};
        for (int k = 0, v = 0; k < 4; k++) if (k != absent) {
            float u = (float)(100+(mirror ? -2 : 2)*(k%2)), w = (float)(20+2*(k/2));
            xyz[3*v] = w; xyz[3*v+1] = u; xyz[3*v+2] = 0;
            uv[2*v] = u; uv[2*v+1] = w; normals[3*v+2] = 1; v++;
        }
        AsmChart chart = {0}; chart.nv = 3; chart.nf = 1; chart.xyz = xyz; chart.uv = uv; chart.nrm = normals; chart.faces = face; chart.area3d = 2;
        run.charts = &chart; run.n_charts = 1;
        MeshBinData sheet = {0}; sheet.nv = 3; sheet.nf = 1; sheet.verts = xyz; sheet.uv = uv; sheet.faces = face;
        AsmVerdictOpts opts; AsmVerdict_default_opts(&opts); AsmVerdictStats st = {0}; int32_t lineage[1] = {0};
        size_t n; long rmin,rmax;
        AvClaim *cells = av_lattice_cells(arena,&run,&opts,&sheet,labels,normals,lineage,&n,&rmin,&rmax,&st);
        AvClaim before[3]; int samples_ok = n == 3;
        if (samples_ok) memcpy(before,cells,sizeof before);
        AvLatticeFaces verdict = av_lattice_faces(arena,&run,cells,n,lineage,NULL,0,NULL,0,5.98,&st);
        size_t solid_nf; int32_t *solid = av_solid_faces(arena,cells,n,5.98,&solid_nf);
        for (int consumer = 0; consumer < 2; consumer++) {
            size_t nf = consumer ? solid_nf : verdict.count;
            const int32_t *out = consumer ? solid : verdict.faces;
            int ok = samples_ok && nf == 1 && !memcmp(before,cells,sizeof before);
            if (nf == 1) {
                int32_t a = out[0], b = out[1], c = out[2];
                ok = ok && a >= 0 && b >= 0 && c >= 0 && (size_t)a < n && (size_t)b < n && (size_t)c < n;
                if (ok) {
                    double ax = (double)av_key_col(cells[b].key)-av_key_col(cells[a].key), ay = (double)av_key_row(cells[b].key)-av_key_row(cells[a].key);
                    double bx = (double)av_key_col(cells[c].key)-av_key_col(cells[a].key), by = (double)av_key_row(cells[c].key)-av_key_row(cells[a].key);
                    ok = ax*by-ay*bx == 1;
                }
            }
            if (!ok) { fprintf(stderr,"  asm_verdict selftest FAIL: partial planar cell absent %d mirror %d consumer %d loses its surface: %zu vertices, %zu faces\n",absent,mirror,consumer,n,nf); fails++; }
        }
        if (samples_ok) for (int test = 0; test < 5; test++) {
            AvClaim trial[3]; memcpy(trial,before,sizeof trial);
            AsmChart test_charts[2] = {chart,chart}; run.charts = test_charts; run.n_charts = 2;
            int32_t test_lineage[2] = {0,test == 1 ? 1 : 0};
            if (test == 1 || test == 2 || test == 4) trial[2].chart = 1;
            if (test >= 2) trial[2].p[2] += 20;
            AsmCorr pair = {0}; pair.valid = 3; pair.run = 0;
            AsmRelation relation = {0}; relation.a = 0; relation.b = 1; relation.corr_count = 1; relation.continuity = ASM_CONT_SOURCE;
            test_charts[0].nv = test_charts[1].nv = 1; test_charts[0].xyz = trial[0].p; test_charts[1].xyz = trial[2].p;
            run.corr = &pair; run.n_corr = 1;
            int64_t key = 1; AvSeamClaim claim = {1,&relation};
            size_t count = test == 0 ? 2 : 3, relations = test == 4 ? 1 : 0;
            AsmVerdictStats audit = {0};
            AvLatticeFaces checked = av_lattice_faces(arena,&run,trial,count,test_lineage,&key,relations,&claim,relations,5.98,&audit);
            size_t solid_count; av_solid_faces(arena,trial,count,5.98,&solid_count);
            size_t expected_verdict = test >= 3 ? 1 : 0, expected_solid = test == 1 ? 1 : 0;
            size_t expected_longs = test == 4 ? 2 : 0;
            int ok = checked.count == expected_verdict && solid_count == expected_solid && checked.nlong == expected_longs;
            if (!ok) { fprintf(stderr,"  asm_verdict selftest FAIL: partial cell authorization absent %d mirror %d case %d: verdict %zu solid %zu long claims %zu\n",absent,mirror,test,checked.count,solid_count,checked.nlong); fails++; }
        }
        Arena_dispose(&arena);
    }
    return fails;
}

static int av_claim_order_selftest(void)
{
    /* Adding an unrelated generated claim, or changing triangle traversal,
     * cannot choose a different physical point for an existing grid cell. */
    AvClaim reference[6] = {{0}};
    for (int k = 0; k < 6; k++) {
        reference[k].key = av_key(0,k < 4 ? 0 : k); reference[k].chart = 0;
        reference[k].area = k < 4 ? 100 : -1; reference[k].p[1] = (float)(10*k);
    }
    qsort(reference,4,sizeof *reference,av_cmp);
    int ok = 1;
    for (int extra = 0; extra < 3; extra++) for (int shift = 0; shift < 4; shift++) {
        AvClaim cells[6]; memcpy(cells,reference,sizeof cells);
        for (int k = 0; k < 4; k++) cells[k] = reference[(k+shift)%4];
        qsort(cells,(size_t)(4+extra),sizeof *cells,av_cmp);
        for (int k = 0; k < 4; k++) if (memcmp(cells+k,reference+k,sizeof *cells)) ok = 0;
    }
    if (!ok) fprintf(stderr,"  asm_verdict selftest FAIL: observed claimant changes with order or unrelated generated samples\n");
    return !ok;
}

static int av_folded_claim_selftest(void)
{
    int fails = 0;
    for (int mirror = 0; mirror < 2; mirror++) for (int direction = 0; direction < 3; direction++) {
        Arena_T arena = Arena_new(); AsmRun run = {0}; run.arena = arena; run.layer_d_global = 10;
        float xyz[24], uv[16], normals[24]; int32_t labels[8] = {0};
        int32_t faces[] = {0,1,3,0,3,2,4,5,7,4,7,6};
        for (int k = 0; k < 8; k++) {
            float u = (float)(4*(k%2)), v = (float)(4*((k%4)/2));
            xyz[3*k] = u+(direction == 1 ? 12*(k/4) : 0); xyz[3*k+1] = v;
            xyz[3*k+2] = (float)((direction == 1 ? 0 : direction == 2 ? 2 : 12)*(k/4));
            uv[2*k] = mirror ? -u : u; uv[2*k+1] = v;
            normals[3*k] = normals[3*k+1] = 0; normals[3*k+2] = 1;
        }
        AsmChart chart = {0}; chart.nv = 8; chart.nf = 4; chart.xyz = xyz; chart.uv = uv; chart.nrm = normals; chart.faces = faces; chart.area3d = 32;
        run.charts = &chart; run.n_charts = 1;
        MeshBinData sheet = {0}; sheet.nv = 8; sheet.nf = 4; sheet.verts = xyz; sheet.uv = uv; sheet.faces = faces;
        AsmVerdictOpts opts; AsmVerdict_default_opts(&opts); AsmVerdictStats st = {0}; int32_t lineage[] = {0};
        size_t n; long rmin,rmax;
        AvClaim *cells = av_lattice_cells(arena,&run,&opts,&sheet,labels,normals,lineage,&n,&rmin,&rmax,&st);
        AvLatticeFaces output = av_lattice_faces(arena,&run,cells,n,lineage,NULL,0,NULL,0,5.98,&st);
        int ok = n == 18 && st.stacked_verts == 9 && st.cells_contested == 9 && output.count == 16;
        for (size_t f = 0; f < output.count; f++) for (int e = 0; e < 3; e++) {
            int32_t a = output.faces[3*f+e], b = output.faces[3*f+(e+1)%3];
            ok = ok && cells[a].p[2] == cells[b].p[2] && av_dist(cells[a].p,cells[b].p) < 5.98;
            if (direction == 1) ok = ok && (cells[a].p[0] >= 12) == (cells[b].p[0] >= 12);
        }
        if (!ok) { fprintf(stderr,"  asm_verdict selftest FAIL: folded chart loses a layer or crosses layers, mirror %d direction %d, %zu vertices, %zu stacks, %zu faces\n",mirror,direction,n,st.stacked_verts,output.count); fails++; }
        Arena_dispose(&arena);
    }
    return fails;
}

static int av_normal_claim_selftest(void)
{
    int fails = 0;
    {
        AvClaim a = {0}, b = {0}; b.p[2] = 6;
        if (av_normal_distance(&a,&b) != 6 || av_normal_distance(&a,&a) != 0) {
            fprintf(stderr,"  asm_verdict selftest FAIL: missing normal erases physical separation\n"); fails++;
        }
    }
    /* Unit vertex normals need not stay unit after interpolation. At the
     * middle column the opposite transverse components cancel, leaving a
     * length-0.5 normal. A six-voxel normal separation is still six voxels. */
    for (int mirror = 0; mirror < 2; mirror++) for (int bent = 0; bent < 2; bent++) for (int separated = 0; separated < 2; separated++) {
        Arena_T arena = Arena_new(); AsmRun run = {0}; run.arena = arena; run.layer_d_global = 10;
        float xyz[24], uv[16], normals[24]; int32_t labels[8], faces[] = {0,1,3,0,3,2,4,5,7,4,7,6};
        for (int k = 0; k < 8; k++) {
            float u = (float)(4*(k%2)), v = (float)(4*((k%4)/2));
            xyz[3*k] = u; xyz[3*k+1] = v; xyz[3*k+2] = (float)((separated ? 6 : 4)*(k/4));
            uv[2*k] = mirror ? -u : u; uv[2*k+1] = v; labels[k] = k/4;
            normals[3*k] = bent ? (float)((k%2 ? -1 : 1)*sqrt(.75)) : 0;
            normals[3*k+1] = 0; normals[3*k+2] = bent ? .5f : 1;
        }
        AsmChart charts[2] = {{0}};
        for (int k = 0; k < 2; k++) {
            charts[k].nv = 4; charts[k].nf = 2; charts[k].xyz = xyz+12*k;
            charts[k].uv = uv+8*k; charts[k].nrm = normals+12*k; charts[k].faces = faces;
            charts[k].area3d = 16;
        }
        run.charts = charts; run.n_charts = 2;
        MeshBinData sheet = {0}; sheet.nv = 8; sheet.nf = 4; sheet.verts = xyz; sheet.uv = uv; sheet.faces = faces;
        AsmVerdictOpts opts; AsmVerdict_default_opts(&opts); AsmVerdictStats st = {0}; int32_t lineage[] = {0,1};
        size_t n; long rmin,rmax;
        av_lattice_cells(arena,&run,&opts,&sheet,labels,normals,lineage,&n,&rmin,&rmax,&st);
        size_t expected = separated ? (bent ? 3 : 9) : 0;
        if (n != 9+expected || st.stacked_verts != expected || st.cells_contested != expected) {
            fprintf(stderr,"  asm_verdict selftest FAIL: interpolated normal hides physical separation, mirror %d bent %d separated %d: %zu vertices, %zu stacks, expected %zu\n",mirror,bent,separated,n,st.stacked_verts,expected); fails++;
        }
        Arena_dispose(&arena);
    }
    return fails;
}

static int av_strip_merge_selftest(void)
{
    Arena_T arena = Arena_new(); AsmRun run = {0}; run.arena = arena; run.layer_d_global = 10;
    AsmChart charts[3] = {{0}}; AsmCorr corr[9] = {{0}};
    float xyz[120], uv[80], normals[120]; int32_t faces[102], labels[40];
    for (int side = 0; side < 3; side++) {
        AsmChart *c = charts+side; c->id = side; c->nv = side == 2 ? 4 : 18; c->nf = side == 2 ? 2 : 16;
        size_t offset = (size_t)(18*side); c->xyz = xyz+3*offset; c->uv = uv+2*offset; c->nrm = normals+3*offset;
        c->placed = 1; c->area3d = 320;
        for (int j = 0; j < (side == 2 ? 2 : 9); j++) for (int col = 0; col < 2; col++) {
            int q = (int)offset+2*j+col;
            double s = side == 2 ? 1.5+3*col : (side ? 4.4 : -10)+10*col;
            double t = side == 2 ? 14+4*j : 4*j;
            xyz[3*q] = (float)t; xyz[3*q+1] = (float)s; xyz[3*q+2] = side == 2 ? 12 : 0;
            uv[2*q] = (float)s; uv[2*q+1] = (float)t;
            normals[3*q] = normals[3*q+1] = 0; normals[3*q+2] = 1; labels[q] = side;
            if (col == 0 && j < (side == 2 ? 1 : 8)) { int32_t f[] = {q,q+1,q+3,q,q+3,q+2}; memcpy(faces+48*side+6*j,f,sizeof f); }
        }
    }
    for (int j = 0; j < 9; j++) {
        corr[j].va = 2*j+1; corr[j].vb = 2*j; corr[j].valid = 3; corr[j].run = j == 0 || j == 8 ? -1 : 0;
        corr[j].gap_a[0] = corr[j].gap_b[0] = 4.4f;
    }
    AsmRelation rel = {0}; rel.a = 0; rel.b = 1; rel.continuity = ASM_CONT_SOURCE; rel.corr_count = rel.n_corr = 9; rel.rms = .1;
    run.charts = charts; run.n_charts = 3; run.rels = &rel; run.n_rels = 1; run.corr = corr; run.n_corr = 9;
    MeshBinData sheet = {0}; sheet.nv = 40; sheet.nf = 34; sheet.verts = xyz; sheet.uv = uv; sheet.faces = faces;
    AsmVerdictOpts opts; AsmVerdict_default_opts(&opts); AsmVerdictStats st = {0}; int32_t lineage[] = {0,0,1};
    size_t n; long rmin,rmax;
    AvClaim *cells = av_lattice_cells(arena,&run,&opts,&sheet,labels,normals,lineage,&n,&rmin,&rmax,&st);
    size_t observed = 0, strip = 0, bridge = 0, generated = 0;
    for (size_t i = 0; i < n; i++) {
        generated += cells[i].area < 0; bridge += cells[i].area == -2;
        if (cells[i].key != av_key(8,1)) continue;
        observed += cells[i].chart == 2 && cells[i].p[2] == 12 && cells[i].area >= 0;
        strip += cells[i].chart == 0 && cells[i].p[2] == 0 && cells[i].area == -1;
    }
    int ok = observed == 1 && strip == 1 && bridge > 0 && generated == st.cells_bridged_uv;
    if (!ok) fprintf(stderr,"  asm_verdict selftest FAIL: axis bridge merge erases a SOURCE strip layer, observed %zu strip %zu bridges %zu generated %zu reported %zu\n",observed,strip,bridge,generated,st.cells_bridged_uv);
    Arena_dispose(&arena); return !ok;
}

static int av_strip_evidence_selftest(void)
{
    int fails = 0;
    for (int mirror = 0; mirror < 2; mirror++) for (int rows_v = 0; rows_v < 2; rows_v++) for (int test = 0; test < 12; test++) {
        Arena_T arena = Arena_new(); AsmRun run = {0}; run.arena = arena;
        AsmChart charts[2] = {{0}}; AsmCorr corr[20] = {{0}};
        float xyz[120], uv[80], local[80], normals[120]; int32_t labels[40];
        for (int side = 0; side < 2; side++) {
            AsmChart *c = charts+side; c->id = side; c->nv = 20; c->placed = 1;
            c->uv = local+40*side; c->xyz = xyz+60*side; c->nrm = normals+60*side;
            c->flags = mirror ? ASM_CHART_MIRROR : 0;
            for (int j = 0; j < 20; j++) {
                int v = 20*side+j; double s = .2+4.4*side, t = .2+4*j;
                xyz[3*v] = (float)t; xyz[3*v+1] = (float)s; xyz[3*v+2] = 0;
                local[2*v] = (float)s; local[2*v+1] = (float)t;
                uv[2*v] = (float)(mirror ? -s : s); uv[2*v+1] = (float)t;
                normals[3*v] = normals[3*v+1] = 0; normals[3*v+2] = 1; labels[v] = side;
            }
        }
        for (int j = 0; j < 20; j++) { corr[j].va = corr[j].vb = j; corr[j].valid = 3; corr[j].gap_a[0] = corr[j].gap_b[0] = 4.4f; }
        AsmRelation rel = {0}; rel.a = 0; rel.b = 1; rel.continuity = ASM_CONT_SOURCE; rel.corr_count = rel.n_corr = 20; rel.rms = .1;
        run.charts = charts; run.n_charts = 2; run.rels = &rel; run.n_rels = 1; run.corr = corr; run.n_corr = 20;
        if (test == 1) rel.continuity = 0;
        if (test == 2) rel.flags = ASM_REL_CONTACT;
        if (test == 3) rel.flags = ASM_REL_DROPPED;
        if (test == 4) rel.flags = ASM_REL_SWITCHED;
        if (test == 5) rel.flags = ASM_REL_WEAK;
        if (test == 6) charts[1].pose_x = 12;
        if (test == 7) for (int j = 20; j < 40; j++) xyz[3*j+1] += 20;
        if (test == 8) for (int j = 0; j < 20; j++) corr[j].run = -1;
        /* A relation may pass its 90-percent support gate while one local
         * witness fails. That tail cannot authorize a generated surface. */
        if (test == 9) corr[9].gap_a[0] += 12;
        if (test == 10) for (int j = 10; j < 20; j++) corr[j].run = 1;
        if (test == 11) for (int j = 20; j < 40; j++) uv[2*j] = mirror ? -.4f : .4f;
        MeshBinData sheet = {0}; sheet.verts = xyz; sheet.uv = uv; sheet.nv = 40;
        float saved_xyz[120], saved_uv[80], saved_local[80], saved_normals[120];
        AsmCorr saved_corr[20]; AsmChart saved_charts[2];
        memcpy(saved_xyz,xyz,sizeof xyz); memcpy(saved_uv,uv,sizeof uv); memcpy(saved_local,local,sizeof local);
        memcpy(saved_normals,normals,sizeof normals); memcpy(saved_corr,corr,sizeof corr); memcpy(saved_charts,charts,sizeof charts);
        size_t n = 0, cap = 64; long rmin = 100000, rmax = -100000;
        AvClaim *claims = ARENA_ALLOC(arena,cap*sizeof *claims);
        av_source_strips(arena,&run,&sheet,labels,normals,2,2,rows_v,&claims,&n,&cap,&rmin,&rmax);
        int expected = test == 0 || test == 9 || test == 10;
        int ok = (n > 0) == expected;
        double rms; size_t support;
        if (test == 9) ok = ok && AsmContinuity_measure(&run,&rel,&rms,&support) && support == 20;
        for (size_t i = 0; i < n; i++) {
            ok = ok && claims[i].area < 0 && claims[i].chart >= 0 && claims[i].chart < 2 &&
                 claims[i].p[1] > .19f && claims[i].p[1] < 4.61f && claims[i].p[2] == 0;
            if (test == 9 && claims[i].p[0] > 32.2 && claims[i].p[0] < 40.2) ok = 0;
            if (test == 10 && claims[i].p[0] > 36.2 && claims[i].p[0] < 40.2) ok = 0;
            if (!rows_v) ok = ok && claims[i].p[0] == 2*av_key_row(claims[i].key);
        }
        ok = ok && !memcmp(saved_xyz,xyz,sizeof xyz) && !memcmp(saved_uv,uv,sizeof uv) && !memcmp(saved_local,local,sizeof local) &&
             !memcmp(saved_normals,normals,sizeof normals) && !memcmp(saved_corr,corr,sizeof corr) && !memcmp(saved_charts,charts,sizeof charts);
        if (!ok) { fprintf(stderr,"  asm_verdict selftest FAIL: trim evidence case %d mirror %d rows_v %d, %zu claims\n",test,mirror,rows_v,n); fails++; }
        Arena_dispose(&arena);
    }
    return fails;
}

static int av_seam_phase_selftest(void)
{
    int fails = 0;
    const double shifts[] = {-.1,0,.1,.3,.7,1.3,1.7,1.9};
    for (int mirror = 0; mirror < 2; mirror++) for (size_t phase = 0; phase < sizeof shifts/sizeof shifts[0]; phase++) {
        Arena_T arena = Arena_new(); AsmRun run = {0}; run.arena = arena; run.layer_d_global = 12; run.continuity_ready = 1;
        AsmChart charts[2] = {{0}}; float xyz[108], uv[72], local[2][36], normals[108];
        int32_t faces[96], labels[36]; AsmCorr corr[9] = {{0}};
        double cs = sqrt(.5), sign = mirror ? -1 : 1;
        for (int side = 0; side < 2; side++) {
            AsmChart *c = &charts[side]; c->id = side; c->nv = 18; c->nf = 16; c->placed = 1; c->placement_state = ASM_PLACE_ROOT;
            c->uv = local[side]; c->xyz = xyz+54*side; c->nrm = normals+54*side; c->area3d = 320;
            c->flags = mirror ? ASM_CHART_MIRROR : 0; c->pose_x = shifts[phase];
            for (int j = 0; j < 9; j++) for (int col = 0; col < 2; col++) {
                int v = 2*j+col, q = 18*side+v;
                double s = (side ? 4.4 : -10)+10*col, t = 4*j;
                xyz[3*q] = (float)t; xyz[3*q+1] = (float)s; xyz[3*q+2] = 0;
                normals[3*q] = normals[3*q+1] = 0; normals[3*q+2] = 1;
                local[side][2*v] = (float)(cs*(s-t)); local[side][2*v+1] = (float)(cs*(s+t));
                uv[2*q] = (float)(sign*local[side][2*v]+shifts[phase]); uv[2*q+1] = local[side][2*v+1]; labels[q] = side;
                if (j < 8 && col == 0) { int32_t f[] = {q,q+1,q+3,q,q+3,q+2}; memcpy(faces+48*side+6*j,f,sizeof f); }
            }
        }
        for (int j = 0; j < 9; j++) {
            corr[j].va = 2*j+1; corr[j].vb = 2*j; corr[j].valid = 3;
            corr[j].gap_a[0] = corr[j].gap_a[1] = corr[j].gap_b[0] = corr[j].gap_b[1] = (float)(4.4*cs);
        }
        AsmRelation rel = {0}; rel.a = 0; rel.b = 1; rel.continuity = ASM_CONT_SOURCE; rel.corr_count = rel.n_corr = 9; rel.rms = .1;
        run.charts = charts; run.n_charts = 2; run.rels = &rel; run.n_rels = 1; run.corr = corr; run.n_corr = 9;
        MeshBinData sheet = {0}; sheet.verts = xyz; sheet.uv = uv; sheet.faces = faces; sheet.nv = 36; sheet.nf = 32;
        AsmVerdictOpts opts; AsmVerdict_default_opts(&opts); AsmVerdictStats st = {0};
        int32_t lineage[] = {0,0}; size_t n, nrel, nclaims; long rmin, rmax;
        int64_t *keys = av_join_keys(arena,&run,&nrel); AvSeamClaim *claims = av_seam_claims(arena,&run,&nclaims);
        AvClaim *cells = av_lattice_cells(arena,&run,&opts,&sheet,labels,normals,lineage,&n,&rmin,&rmax,&st);
        AvLatticeFaces mesh = av_lattice_faces(arena,&run,cells,n,lineage,keys,nrel,claims,nclaims,5.98,&st);
        float *points = ARENA_ALLOC(arena,3*n*sizeof(float)); int32_t *owners = ARENA_ALLOC(arena,n*sizeof(int32_t));
        for (size_t i = 0; i < n; i++) { memcpy(points+3*i,cells[i].p,3*sizeof(float)); owners[i] = cells[i].chart; }
        AsmContinuityStats audit; AsmContinuity_audit(&run,0,owners,points,n,mesh.faces,mesh.count,NULL,&audit);
        double rms; size_t support;
        int ok = AsmContinuity_measure(&run,&rel,&rms,&support) && support == 9 && rms < 1e-4 && audit.confetti_free;
        if (!ok) { fprintf(stderr,"  asm_verdict selftest FAIL: planar trim phase %.2f mirror %d, bridges %zu, source RMS %.9g, metadata %zu, wrong %zu\n",
                          shifts[phase],mirror,st.cells_bridged_uv,rms,audit.metadata_only_joins,audit.wrong_placements); fails++; }
        Arena_dispose(&arena);
    }
    return fails;
}

int AsmVerdict_selftest(void)
{
    int fails = av_normal_claim_selftest()+av_partial_cell_selftest()+av_claim_order_selftest()+av_folded_claim_selftest()+av_strip_merge_selftest()+av_seam_phase_selftest()+av_strip_evidence_selftest()+av_gate_word_selftest();
    {
        /* A reused or mismatched export must be rejected before writing
         * beyond its declared vertex count or trusting wrong chart IDs. */
        AsmRun run = {0}; AsmChart charts[3] = {{0}};
        float normals[9] = {0,0,1,0,1,0,1,0,0}, uv[6] = {0}, output[12], original[12];
        int32_t labels[3] = {0,0,1}; MeshBinData sheet = {0};
        for (int i = 0; i < 2; i++) {
            charts[i].id = i; charts[i].placed = 1; charts[i].nv = i ? 1 : 2;
            charts[i].placement_state=ASM_PLACE_SEAM;charts[i].component=17*i;
            charts[i].nrm = normals+6*i; charts[i].uv = uv+4*i;
        }
        /* Registered admission may retain a different legacy component.
         * Floating local poses remain extras, including incomplete normals. */
        charts[2].placed=1;charts[2].uv=uv;charts[2].nv=1;charts[2].id=2;
        run.charts = charts; run.n_charts = 3;
        for (int i = 0; i < 12; i++) output[i] = original[i] = (float)(100+i);
        sheet.nv = 2;
        if (av_copy_sheet_normals(&run,&sheet,labels,output) || memcmp(output,original,sizeof original)) {
            fprintf(stderr,"  asm_verdict selftest FAIL: stale sheet count changes output or its canaries\n"); fails++;
        }
        memcpy(output,original,sizeof output); sheet.nv = 3; labels[2] = 0;
        if (av_copy_sheet_normals(&run,&sheet,labels,output) || memcmp(output,original,sizeof original)) {
            fprintf(stderr,"  asm_verdict selftest FAIL: wrong sheet chart IDs accepted\n"); fails++;
        }
        labels[2] = 1; memcpy(output,original,sizeof output);
        charts[1].nrm = NULL;
        if (av_copy_sheet_normals(&run,&sheet,labels,output) || memcmp(output,original,sizeof original)) {
            fprintf(stderr,"  asm_verdict selftest FAIL: missing chart normals accepted\n"); fails++;
        }
        charts[1].nrm = normals+6;
        if (!av_copy_sheet_normals(&run,&sheet,labels,output) || memcmp(output,normals,sizeof normals) ||
            memcmp(output+9,original+9,3*sizeof(float))) {
            fprintf(stderr,"  asm_verdict selftest FAIL: valid sheet normals or canaries changed\n"); fails++;
        }
    }
    {
        /* Every perimeter step is sound in a sheared planar cell. Only
         * one diagonal exceeds the solid gate: choosing the other must
         * retain both triangles and all four corners. Real perimeter gaps
         * and cells with two excessive diagonals still remain holes. */
        for (int scenario = 0; scenario < 4; scenario++) for (int mirror = 0; mirror < 2; mirror++) {
            AvClaim cells[13] = {{0}}; const size_t q[4] = {3,5,9,12};
            double shift = scenario == 0 ? 3 : scenario == 2 ? 8 : 0;
            double width = scenario == 3 ? 5 : 2, height = scenario == 3 ? 5 : 4;
            for (int k = 0; k < 4; k++) {
                cells[q[k]].p[0] = (float)((k/2)*shift+width*((k%2)^mirror));
                cells[q[k]].p[1] = (float)((k/2)*height);
            }
            int32_t output[8] = {1234567,-1,-1,-1,-1,-1,-1,7654321};
            size_t count = av_solid_quad(cells,q,5.98,output+1);
            int ok = count == (size_t)(scenario < 2 ? 2 : 0) && output[0] == 1234567 && output[7] == 7654321;
            double area = 0; unsigned seen = 0;
            for (size_t f = 0; f < count; f++) {
                int index[3];
                for (int k = 0; k < 3; k++) {
                    int32_t vertex = output[1+3*f+k]; index[k] = -1;
                    for (int j = 0; j < 4; j++) if (vertex == (int32_t)q[j]) index[k] = j;
                    if (index[k] >= 0) seen |= 1u << index[k]; else ok = 0;
                }
                if (!ok) break;
                int a = index[0], b = index[1], c = index[2];
                area += .5*((b%2-a%2)*(c/2-a/2)-(b/2-a/2)*(c%2-a%2));
                for (int k = 0; k < 3; k++) ok = ok && av_dist(cells[q[index[k]]].p,cells[q[index[(k+1)%3]]].p) < 5.98;
            }
            if (count) ok = ok && seen == 15 && area == 1;
            else for (int k = 1; k <= 6; k++) ok = ok && output[k] == -1;
            if (!ok) {
                fprintf(stderr,"  asm_verdict selftest FAIL: solid cell diagonal scenario %d mirror %d, %zu faces, %.1f area\n",
                        scenario,mirror,count,area); fails++;
            }
        }
    }
    {
        /* A planar seam with sound perimeter edges must not acquire a bad
         * diagonal merely from the direction in which UV columns run. A
         * truly excessive seam step remains visible in either direction. */
        Arena_T arena = Arena_new(); AsmRun run = {0}; run.arena = arena;
        AsmChart charts[2] = {{0}}; AsmCorr corr[9] = {{0}}; float xyz[2][27] = {{0}};
        for (int side = 0; side < 2; side++) {
            charts[side].id = side; charts[side].nv = 9; charts[side].xyz = xyz[side];
            for (int k = 0; k < 9; k++) { xyz[side][3*k] = (float)(4*k); xyz[side][3*k+1] = (float)(5*side); }
        }
        for (int k = 0; k < 9; k++) { corr[k].va = corr[k].vb = k; corr[k].valid = 3; corr[k].run = 0; }
        AsmRelation rel = {0}; rel.a = 0; rel.b = 1; rel.continuity = ASM_CONT_SOURCE; rel.corr_count = 9;
        run.charts = charts; run.n_charts = 2; run.corr = corr; run.n_corr = 9; run.rels = &rel; run.n_rels = 1;
        size_t nclaims, nrel; AvSeamClaim *claims = av_seam_claims(arena,&run,&nclaims);
        int64_t *keys = av_join_keys(arena,&run,&nrel); int32_t lineage[2] = {0,0};
        for (int shift = 6; shift <= 10; shift += 4) for (int mirror = 0; mirror < 2; mirror++) {
            for (int k = 0; k < 9; k++) xyz[1][3*k] = (float)(4*k+shift);
            AvClaim cells[4] = {{0}};
            for (int k = 0; k < 4; k++) {
                cells[k].key = av_key(k/2,k%2); cells[k].chart = k/2;
                cells[k].p[0] = (float)((k/2)*shift+2*((k%2)^mirror)); cells[k].p[1] = (float)((k/2)*5);
            }
            AsmVerdictStats st = {0};
            AvLatticeFaces mesh = av_lattice_faces(arena,&run,cells,4,lineage,keys,nrel,claims,nclaims,5.98,&st);
            size_t long_edges = 0; double uv_area = 0; unsigned seen = 0;
            for (size_t f = 0; f < mesh.count; f++) {
                int32_t a = mesh.faces[3*f], b = mesh.faces[3*f+1], c = mesh.faces[3*f+2];
                if (a < 0 || a >= 4 || b < 0 || b >= 4 || c < 0 || c >= 4) { fails++; continue; }
                uv_area += .5*((b%2-a%2)*(c/2-a/2)-(b/2-a/2)*(c%2-a%2));
                for (int e = 0; e < 3; e++) {
                    int32_t u = mesh.faces[3*f+e], v = mesh.faces[3*f+(e+1)%3]; seen |= 1u << u;
                    if (cells[u].chart != cells[v].chart && av_dist(cells[u].p,cells[v].p) > 8) long_edges++;
                }
            }
            if (mesh.count != 2 || seen != 15 || uv_area != 1 || (shift == 6 ? long_edges != 0 : long_edges == 0)) {
                fprintf(stderr,"  asm_verdict selftest FAIL: planar seam diagonal shift %d mirror %d: %zu faces, %.1f area, %zu edges over 8 vox\n",
                        shift,mirror,mesh.count,uv_area,long_edges); fails++;
            }
        }
        Arena_dispose(&arena);
    }
    {
        /* A rejected quad has no emitted edges. Conversely, a kept source
         * seam must retain its long edges in both the mesh and the ledger. */
        Arena_T arena = Arena_new(); AsmRun run = {0}; run.arena = arena;
        AsmChart charts[3] = {{0}}; AsmCorr corr[9] = {{0}}; float xyz[2][27] = {{0}};
        for (int side = 0; side < 2; side++) {
            charts[side].id = side; charts[side].nv = 9; charts[side].xyz = xyz[side];
            for (int k = 0; k < 9; k++) xyz[side][3*k+1] = (float)(4*k);
        }
        for (int k = 0; k < 9; k++) { corr[k].va = corr[k].vb = k; corr[k].valid = 3; corr[k].run = 0; }
        AsmRelation rel = {0}; rel.a = 0; rel.b = 1; rel.continuity = ASM_CONT_SOURCE; rel.corr_count = 9;
        run.charts = charts; run.n_charts = 3; run.corr = corr; run.n_corr = 9; run.rels = &rel; run.n_rels = 1;
        size_t nclaims, nrel; AvSeamClaim *claims = av_seam_claims(arena,&run,&nclaims);
        int64_t *keys = av_join_keys(arena,&run,&nrel);
        for (int gap = 7; gap <= 10; gap += 3) for (int rejected = 0; rejected < 3; rejected++) {
            for (int k = 0; k < 9; k++) xyz[1][3*k] = (float)gap;
            AvClaim cells[4] = {{0}}; int32_t lineage[3] = {0,0,rejected == 1 ? 1 : 0};
            for (int k = 0; k < 4; k++) {
                cells[k].key = av_key(k/2,k%2); cells[k].chart = k%2;
                cells[k].p[0] = (float)((k%2)*gap); cells[k].p[1] = (float)((k/2)*2);
            }
            if (rejected) cells[3].chart = 2;
            if (rejected == 2) cells[3].p[0] = 100; /* same lineage, failed geometry */
            AsmVerdictStats st = {0};
            AvLatticeFaces mesh = av_lattice_faces(arena,&run,cells,4,lineage,keys,nrel,claims,nclaims,5.98,&st);
            int bad = rejected ? mesh.count != 0 || mesh.nlong != 0 || st.long_join_edges_tail != 0
                               : mesh.count != 2 || mesh.nlong != 3 || st.long_join_edges_tail != (gap == 7 ? 3 : 0);
            if (!rejected && mesh.count == 2) {
                const int32_t expected[6] = {0,1,3,0,3,2};
                bad |= memcmp(mesh.faces,expected,sizeof expected) != 0;
                for (size_t i = 0; i < mesh.nlong; i++) bad |= mesh.longs[i].pair != 1 || mesh.longs[i].step < gap;
            }
            bad |= st.quads_refused != (rejected == 2) || st.quads_refused_lineage != (rejected == 1);
            if (bad) {
                fprintf(stderr,"  asm_verdict selftest FAIL: quad ledger gap %d refusal %d: %zu faces, %zu long edges, %zu tail edges\n",
                        gap,rejected,mesh.count,mesh.nlong,st.long_join_edges_tail); fails++;
            }
        }
        Arena_dispose(&arena);
    }
    {
        /* A valid join spans the measured boundary, not an interior island
         * of the same chart that happens to touch it in the raster. */
        Arena_T arena = Arena_new(); AsmRun run = {0}; run.arena = arena;
        float xyz[2][27]; AsmChart charts[2] = {{0}}; AsmCorr corr[9] = {{0}};
        for (int side = 0; side < 2; side++) {
            charts[side].id = side; charts[side].nv = 9; charts[side].xyz = xyz[side];
            for (int k = 0; k < 9; k++) { xyz[side][k*3] = (float)(4*k); xyz[side][k*3+1] = (float)(4*side); xyz[side][k*3+2] = 0; }
        }
        for (int k = 0; k < 9; k++) { corr[k].va = corr[k].vb = k; corr[k].valid = 3; corr[k].run = 0; }
        AsmRelation rel = {0}; rel.a = 0; rel.b = 1; rel.continuity = ASM_CONT_SOURCE; rel.corr_count = 9;
        run.charts = charts; run.n_charts = 2; run.corr = corr; run.n_corr = 9; run.rels = &rel; run.n_rels = 1;
        size_t n; AvSeamClaim *claims = av_seam_claims(arena,&run,&n);
        AvClaim a = {0}, b = {0}; a.chart = 0; b.chart = 1;
        a.p[0] = 3; a.p[1] = 1; b.p[0] = 5; b.p[1] = 3;
        if (!av_seam_claimed_at(&run,claims,n,&a,&b,6) || !av_seam_claimed_at(&run,claims,n,&b,&a,6)) fails++;
        b.p[1] = 20;
        if (av_seam_claimed_at(&run,claims,n,&a,&b,6)) fails++;
        /* Keep long edges on the claimed run visible to the audit. This
         * must not turn into a distance-based filter for bad placements. */
        b.p[0] = 32; b.p[1] = 4;
        if (av_dist(a.p,b.p) <= 8 || !av_seam_claimed_at(&run,claims,n,&a,&b,6)) fails++;
        /* A layout relation must not make a rejected source seam, or an
         * unrelated interior contact, an unconditional raster connection. */
        AsmRelation rels[2] = {rel,{0}};
        rels[0].flags = ASM_REL_WEAK;
        rels[1].a = 0; rels[1].b = 1; rels[1].flags = ASM_REL_LAYOUT;
        run.rels = rels; run.n_rels = 2;
        claims = av_seam_claims(arena,&run,&n);
        if (n != 0 || av_seam_claimed_at(&run,claims,n,&a,&b,6)) fails++;
        /* When the source seam is accepted, its locality and the audit's
         * visibility of true long edges survive the duplicate layout join. */
        rels[0].flags = 0;
        claims = av_seam_claims(arena,&run,&n);
        if (n != 1 || !av_seam_claimed_at(&run,claims,n,&a,&b,6)) fails++;
        b.p[1] = 20;
        if (av_seam_claimed_at(&run,claims,n,&a,&b,6) || av_seam_claimed_at(&run,claims,n,&b,&a,6)) fails++;
        Arena_dispose(&arena);
    }
    /* the axis frame through a tilted straight axis: the axis lands on +z through the origin and
     * distances are preserved */
    {
        Arena_T arena = Arena_new();
        double a[3] = { cos(0.6), sin(0.6) * 0.8, sin(0.6) * 0.6 };
        double pt[3] = { 100.0, 200.0, 300.0 };
        AsmAxis *axis = AsmAxis_line(arena, pt, a, 0.0, 1000.0);
        float on[3] = { (float)(pt[0] + 250.0 * a[0]), (float)(pt[1] + 250.0 * a[1]), (float)(pt[2] + 250.0 * a[2]) };
        av_to_frame(axis, on);
        if (fabs(on[1]) > 1e-3 || fabs(on[2]) > 1e-3) { fprintf(stderr, "  asm_verdict selftest FAIL: a point on the axis is off +z in the frame (%.4f %.4f)\n", on[1], on[2]); fails++; }
        float p[3] = { 12.5f, -3.0f, 7.25f }, q[3] = { 1.0f, 2.0f, 3.0f };
        double d0 = av_dist(p, q);
        av_to_frame(axis, p); av_to_frame(axis, q);
        if (fabs(av_dist(p, q) - d0) > 1e-3) { fprintf(stderr, "  asm_verdict selftest FAIL: the frame changes a distance (%.6f vs %.6f)\n", av_dist(p, q), d0); fails++; }
        Arena_dispose(&arena);
    }
    /* joins: sorted pair keys answer membership both ways, contacts / weak / switched are not joins */
    {
        Arena_T arena = Arena_new();
        AsmRun run;
        memset(&run, 0, sizeof run);
        run.arena = arena;
        AsmRelation rels[4];
        memset(rels, 0, sizeof rels);
        rels[0].a = 3; rels[0].b = 1;
        rels[1].a = 2; rels[1].b = 5; rels[1].flags = ASM_REL_WEAK;
        rels[2].a = 4; rels[2].b = 6; rels[2].flags = ASM_REL_SWITCHED;
        rels[3].a = 7; rels[3].b = 8; rels[3].flags = ASM_REL_CONTACT;
        run.rels = rels; run.n_rels = 4;
        size_t nrel = 0;
        int64_t *keys = av_join_keys(arena, &run, &nrel);
        if (nrel != 1 || !av_joined(keys, nrel, 1, 3) || !av_joined(keys, nrel, 3, 1) || av_joined(keys, nrel, 2, 5) ||
            av_joined(keys, nrel, 4, 6) || av_joined(keys, nrel, 7, 8) || !av_joined(keys, nrel, 9, 9)) {
            fprintf(stderr, "  asm_verdict selftest FAIL: join keys (%zu)\n", nrel); fails++;
        }
        /* lineages follow the joins only */
        AsmChart charts[9];
        float uv[2] = { 0.0f, 0.0f };
        memset(charts, 0, sizeof charts);
        for (int i = 0; i < 9; i++) charts[i].uv = uv;   /* in the layout */
        run.charts = charts; run.n_charts = 9;
        int32_t *lin = av_lineages(arena, &run);
        if (lin[1] != lin[3] || lin[2] == lin[5] || lin[4] == lin[6] || lin[7] == lin[8]) { fprintf(stderr, "  asm_verdict selftest FAIL: lineages\n"); fails++; }
        Arena_dispose(&arena);
    }
    /* the pick: the same chart's vertex in a cell wins over the primary */
    {
        AvClaim cells[3];
        memset(cells, 0, sizeof cells);
        cells[0].key = av_key(0, 0); cells[0].chart = 1;
        cells[1].key = av_key(0, 1); cells[1].chart = 1;   /* primary */
        cells[2].key = av_key(0, 1); cells[2].chart = 4;   /* stacked */
        if (av_pick(cells, 3, av_key(0, 1), 4) != 2 || av_pick(cells, 3, av_key(0, 1), 9) != 1 || av_pick(cells, 3, av_key(0, 2), 1) != (size_t)-1) {
            fprintf(stderr, "  asm_verdict selftest FAIL: pick\n"); fails++;
        }
    }
    if (fails == 0) fprintf(stderr, "  asm_verdict selftest: all passed\n");
    return fails;
}
