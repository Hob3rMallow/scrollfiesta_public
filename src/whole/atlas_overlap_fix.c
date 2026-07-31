#include "atlas_overlap_fix.h"
#include "../common/ves_platform.h"
#include "atlas_collision_register.h"

#include "atlas_register.h"
#include "cube_register.h"

#include "../common/union_find.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define AOF_2PI 6.283185307179586476925286766559

void AtlasOverlapFixOptions_default(AtlasOverlapFixOptions *opts)
{
    if (opts == NULL) return;
    memset(opts, 0, sizeof(*opts));
    opts->aspect_ratio_kill = 8.0;
    opts->floater_max_faces = 0;
    opts->floater_group_max_faces = 0;
    opts->chart_qc_kill = 0.0;
    /* The XYZ seam weld bridges at rho up to 3.0, so it joins surfaces up to
     * 2*rho = 6 vox apart -- and that ceiling exists because the inter-wrap
     * clearance is 7 vox.  Matching it here accepts every weld the real weld
     * would make while staying under the distance at which a "neighbour" could
     * be the next wrap in.  Raising this past 7 would fuse wraps into one
     * group and freeze a wrong placement. */
    opts->neighbour_distance = 6.0;
    opts->neighbour_normal_dot = 0.7;
    opts->neighbour_min_shared = 3;
    opts->max_shift_wraps = 4;
    opts->overlap_cell_size = 4.0;
    opts->park_margin = 16.0;
    opts->park_unplaceable = 1;
    opts->wind_correct = 1;
    opts->wind_only = 0;
    opts->wind_residual_limit = 0.25;
    opts->wind_cap = 8;
    opts->phase_targets = 0;
    opts->tabu = 0;
    opts->tabu_span = 4;
    opts->tabu_tenure = 40;
    opts->tabu_max_iters = 2000;
    opts->tabu_stall = 300;
    opts->tabu_cell = 4.0;
    opts->tabu_ladder_override = 0;
    opts->tabu_ladder_step = 0.0;
    opts->tabu_ladder_u0 = 0.0;
    opts->tabu_active_core = 0;
    opts->lambda_ov = 1.0;
    opts->lambda_nb = 4.0;
    opts->lambda_rad = 2.0;
    opts->lambda_gauge = 0.0;
    opts->nb_unhappy_w = 0.25;
    opts->lateral_w = 8.0;
    opts->phase_all_groups = 0;
    opts->phase_lateral_targets = 0;
    opts->phase_regions = 0;
    opts->phase_flat_lock = 0;
    opts->lateral_u_gap = 200.0;
    opts->family_w = 8.0;
    opts->relayout_seam_cap = 25.0;
}

/* ========================================================================== */
/* Axis geometry                                                              */
/* ========================================================================== */

static double aof_axis_coord(const ScaffoldCalib *cal, const float *p)
{
    double d0 = (double)p[0] - (double)cal->axis_point[0];
    double d1 = (double)p[1] - (double)cal->axis_point[1];
    double d2 = (double)p[2] - (double)cal->axis_point[2];
    return d0 * (double)cal->axis_dir[0] +
           d1 * (double)cal->axis_dir[1] +
           d2 * (double)cal->axis_dir[2];
}

static double aof_radius(const ScaffoldCalib *cal, const float *p)
{
    double d0 = (double)p[0] - (double)cal->axis_point[0];
    double d1 = (double)p[1] - (double)cal->axis_point[1];
    double d2 = (double)p[2] - (double)cal->axis_point[2];
    double along = d0 * (double)cal->axis_dir[0] +
                   d1 * (double)cal->axis_dir[1] +
                   d2 * (double)cal->axis_dir[2];
    double r0 = d0 - along * (double)cal->axis_dir[0];
    double r1 = d1 - along * (double)cal->axis_dir[1];
    double r2 = d2 - along * (double)cal->axis_dir[2];
    return sqrt(r0 * r0 + r1 * r1 + r2 * r2);
}

/* ========================================================================== */
/* Sorted-bucket spatial index.  Records are (cell key, item); sorting by key  */
/* makes every cell a contiguous run found by binary search.  This avoids a    */
/* dense grid, whose u extent over a whole scroll would be hundreds of         */
/* thousands of cells wide.                                                   */
/* ========================================================================== */

typedef struct {
    uint64_t key;
    int32_t  item;
} AofRecord;

static int aof_compare_record(const void *pa, const void *pb)
{
    const AofRecord *a = (const AofRecord *)pa;
    const AofRecord *b = (const AofRecord *)pb;
    if (a->key != b->key) return a->key < b->key ? -1 : 1;
    return a->item < b->item ? -1 : (a->item > b->item ? 1 : 0);
}

/* First index whose key >= target, or n when there is none. */
static size_t aof_lower_bound(const AofRecord *record, size_t n, uint64_t target)
{
    size_t lo = 0, hi = n;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (record[mid].key < target) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

static uint64_t aof_key2(int32_t cu, int32_t cv)
{
    return ((uint64_t)(uint32_t)cu << 32) | (uint64_t)(uint32_t)cv;
}

static uint64_t aof_key3(int32_t cz, int32_t cy, int32_t cx)
{
    return ((uint64_t)(uint32_t)cz << 42) |
           ((uint64_t)((uint32_t)cy & 0x1FFFFFu) << 21) |
           (uint64_t)((uint32_t)cx & 0x1FFFFFu);
}

static int32_t aof_cell_of(double value, double origin, double cell)
{
    double c = floor((value - origin) / cell);
    if (c < -1.0e9) c = -1.0e9;
    if (c > 1.0e9) c = 1.0e9;
    return (int32_t)c;
}

/* ========================================================================== */
/* Step 1: charts = connected components of the retained mesh                  */
/* ========================================================================== */

static size_t aof_label_charts(Arena_T arena, const PieceSet *ps,
                               const uint8_t *face_keep, int32_t *vertex_chart)
{
    UnionFind uf = UF_new(arena, (int32_t)ps->nv);
    for (size_t f = 0; f < ps->nf; f++) {
        if (!face_keep[f]) continue;
        int32_t a = ps->faces[f * 3];
        int32_t b = ps->faces[f * 3 + 1];
        int32_t c = ps->faces[f * 3 + 2];
        uf_union(&uf, a, b);
        uf_union(&uf, a, c);
    }
    /* Only vertices carried by a surviving face get a chart. */
    uint8_t *used = (uint8_t *)ARENA_CALLOC(arena, ps->nv, sizeof(uint8_t));
    for (size_t f = 0; f < ps->nf; f++) {
        if (!face_keep[f]) continue;
        for (int k = 0; k < 3; k++) used[ps->faces[f * 3 + (size_t)k]] = 1;
    }
    int32_t *root_chart = (int32_t *)ARENA_ALLOC(
        arena, ps->nv * sizeof(int32_t));
    for (size_t i = 0; i < ps->nv; i++) root_chart[i] = -1;

    int32_t ncharts = 0;
    for (size_t i = 0; i < ps->nv; i++) {
        if (!used[i]) { vertex_chart[i] = -1; continue; }
        int32_t root = uf_find(&uf, (int32_t)i);
        if (root_chart[root] < 0) root_chart[root] = ncharts++;
        vertex_chart[i] = root_chart[root];
    }
    return (size_t)ncharts;
}

/* ========================================================================== */
/* Step 2: kill triangles stretched hopelessly by the parameterization         */
/* ========================================================================== */

static size_t aof_kill_bad_triangles(const PieceSet *ps, const double *u,
                                     const double *v, double max_stretch,
                                     uint8_t *face_keep)
{
    size_t killed = 0;
    for (size_t f = 0; f < ps->nf; f++) {
        if (!face_keep[f]) continue;
        int32_t vi[3] = {
            ps->faces[f * 3], ps->faces[f * 3 + 1], ps->faces[f * 3 + 2]
        };
        int valid = 1;
        for (int k = 0; k < 3; k++) {
            if (vi[k] < 0 || (size_t)vi[k] >= ps->nv ||
                !isfinite(u[vi[k]]) || !isfinite(v[vi[k]]))
                valid = 0;
        }
        if (!valid) { face_keep[f] = 0; killed++; continue; }

        int overstretched = 0;
        for (int e = 0; e < 3; e++) {
            int32_t a = vi[e], b = vi[(e + 1) % 3];
            double du = u[b] - u[a], dv = v[b] - v[a];
            double uv_length = sqrt(du * du + dv * dv);
            const float *pa = &ps->verts[(size_t)a * 3];
            const float *pb = &ps->verts[(size_t)b * 3];
            double dx = (double)pb[0] - (double)pa[0];
            double dy = (double)pb[1] - (double)pa[1];
            double dz = (double)pb[2] - (double)pa[2];
            double xyz_length = sqrt(dx * dx + dy * dy + dz * dz);
            if (!(xyz_length > 1.0e-9) || !isfinite(uv_length) ||
                uv_length > max_stretch * xyz_length) {
                overstretched = 1;
                break;
            }
        }
        if (overstretched) {
            face_keep[f] = 0;
            killed++;
        }
    }
    return killed;
}
/* Sander et al. 2001 quasi-conformal stretch of the XYZ-isometric triangle
 * mapped into (u,v).  A source sliver mapped isometrically remains qc=1; only
 * degeneration introduced by the projection makes the ratio explode. */
static double aof_face_qc(const PieceSet *ps, const double *u, const double *v,
                          int32_t a, int32_t b, int32_t c)
{
    const float *p0 = &ps->verts[(size_t)a * 3];
    const float *p1 = &ps->verts[(size_t)b * 3];
    const float *p2 = &ps->verts[(size_t)c * 3];
    double e1[3] = {p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
    double e2[3] = {p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]};
    double len = sqrt(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]);
    if (!(len > 1.0e-12)) return DBL_MAX;

    double d = (e1[0]*e2[0] + e1[1]*e2[1] + e1[2]*e2[2]) / len;
    double h2 = e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2] - d*d;
    if (!(h2 > 1.0e-24)) return DBL_MAX;
    double h = sqrt(h2);

    double du1 = u[b] - u[a], du2 = u[c] - u[a];
    double dv1 = v[b] - v[a], dv2 = v[c] - v[a];
    double j00 = du1 / len;
    double j01 = (du2 - du1*d/len) / h;
    double j10 = dv1 / len;
    double j11 = (dv2 - dv1*d/len) / h;
    double e = j00*j00 + j10*j10;
    double g = j01*j01 + j11*j11;
    double f = j00*j01 + j10*j11;
    double disc = sqrt(fmax(0.0, 0.25*(e-g)*(e-g) + f*f));
    double lmax = 0.5*(e+g) + disc;
    double lmin = 0.5*(e+g) - disc;
    if (!(lmin > 1.0e-24) || !isfinite(lmax)) return DBL_MAX;
    return sqrt(lmax / lmin);
}



/* ========================================================================== */
/* Step 3: chart geometry                                                      */
/* ========================================================================== */

typedef struct {
    double normal[3];
    double mean_radius;
    double qc_max;
    size_t nfaces;
    size_t nvertices;
} AofChartInfo;

static void aof_chart_info(const PieceSet *ps, const ScaffoldCalib *cal,
                           const int32_t *vertex_chart, size_t ncharts,
                           const double *u, const double *v,
                           const uint8_t *face_keep, AofChartInfo *info)
{
    memset(info, 0, ncharts * sizeof(*info));
    for (size_t i = 0; i < ps->nv; i++) {
        int32_t chart = vertex_chart[i];
        if (chart < 0) continue;
        info[chart].mean_radius += aof_radius(cal, &ps->verts[i * 3]);
        info[chart].nvertices++;
    }
    for (size_t chart = 0; chart < ncharts; chart++)
        if (info[chart].nvertices > 0)
            info[chart].mean_radius /= (double)info[chart].nvertices;

    for (size_t f = 0; f < ps->nf; f++) {
        if (!face_keep[f]) continue;
        int32_t a = ps->faces[f * 3];
        int32_t b = ps->faces[f * 3 + 1];
        int32_t c = ps->faces[f * 3 + 2];
        int32_t chart = vertex_chart[a];
        if (chart < 0) continue;
        const float *pa = &ps->verts[(size_t)a * 3];
        const float *pb = &ps->verts[(size_t)b * 3];
        const float *pc = &ps->verts[(size_t)c * 3];
        double ab[3] = {pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]};
        double ac[3] = {pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]};
        info[chart].normal[0] += ab[1] * ac[2] - ab[2] * ac[1];
        info[chart].normal[1] += ab[2] * ac[0] - ab[0] * ac[2];
        info[chart].normal[2] += ab[0] * ac[1] - ab[1] * ac[0];
        info[chart].nfaces++;
        double qc = aof_face_qc(ps, u, v, a, b, c);
        if (qc > info[chart].qc_max) info[chart].qc_max = qc;

    }
    for (size_t chart = 0; chart < ncharts; chart++) {
        double *n = info[chart].normal;
        double length = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        if (length > 1.0e-12) {
            n[0] /= length; n[1] /= length; n[2] /= length;
        }
    }
}

static size_t aof_chart_size_bin(size_t nfaces)
{
    static const size_t hi[9] = {1, 3, 7, 15, 31, 63, 127, 255, 511};
    for (size_t b = 0; b < 9; b++) if (nfaces <= hi[b]) return b;
    return 9;
}

static size_t aof_chart_qc_bin(double qc)
{
    double hi = 2.0;
    for (size_t b = 0; b < 8; b++, hi *= 2.0)
        if (qc <= hi) return b;
    return 8;
}

/* Cull before neighbour construction: these components must never enter the
 * winding graph, occupancy hash, or continuous fit.  The two predicates stay
 * separately switchable so each ablation has an attributable result. */
static size_t aof_kill_bad_charts(Arena_T arena, const PieceSet *ps,
                                  const int32_t *vertex_chart, size_t ncharts,
                                  const AofChartInfo *chart,
                                  const AtlasOverlapFixOptions *opts,
                                  uint8_t *face_keep,
                                  AtlasOverlapFixStats *stats)
{
    uint8_t *kill = (uint8_t *)ARENA_CALLOC(
        arena, ncharts ? ncharts : 1, sizeof(uint8_t));
    for (size_t c = 0; c < ncharts; c++) {
        stats->chart_size_hist[aof_chart_size_bin(chart[c].nfaces)]++;
        stats->chart_qc_hist[aof_chart_qc_bin(chart[c].qc_max)]++;
        int small = opts->floater_max_faces > 0 &&
                    chart[c].nfaces <= opts->floater_max_faces;
        int bad_qc = opts->chart_qc_kill > 1.0 &&
                     chart[c].qc_max > opts->chart_qc_kill;
        if (small) stats->charts_killed_small++;
        if (bad_qc) stats->charts_killed_qc++;
        kill[c] = (uint8_t)(small || bad_qc);
    }

    size_t killed = 0;
    for (size_t f = 0; f < ps->nf; f++) {
        if (!face_keep[f]) continue;
        int32_t c = vertex_chart[ps->faces[f * 3]];
        if (c >= 0 && kill[c]) {
            face_keep[f] = 0;
            killed++;
        }
    }
    return killed;
}

static size_t aof_group_size_bin(size_t nfaces)
{
    static const size_t hi[11] = {
        1, 3, 7, 15, 31, 63, 127, 255, 511, 1023, 2047
    };
    for (size_t b = 0; b < 11; b++) if (nfaces <= hi[b]) return b;
    return 11;
}

/* A tiny disconnected chart can be a valid edge piece once authentic welds
 * attach it to a real sheet.  Cull on the TOTAL welded-group population after
 * those links exist: this removes isolated toys without punching holes in the
 * major sheets.  The production histogram has a clean 456 -> 1218 face gap. */
static size_t aof_kill_small_groups(Arena_T arena, const PieceSet *ps,
                                    const int32_t *vertex_chart,
                                    const int32_t *chart_group,
                                    size_t ngroups, size_t max_faces,
                                    uint8_t *face_keep,
                                    AtlasOverlapFixStats *stats)
{
    size_t *nfaces = (size_t *)ARENA_CALLOC(
        arena, ngroups ? ngroups : 1, sizeof(size_t));
    for (size_t f = 0; f < ps->nf; f++) {
        if (!face_keep[f]) continue;
        int32_t c = vertex_chart[ps->faces[f * 3]];
        if (c < 0) continue;
        int32_t g = chart_group[c];
        if (g >= 0) nfaces[g]++;
    }

    uint8_t *kill = (uint8_t *)ARENA_CALLOC(
        arena, ngroups ? ngroups : 1, sizeof(uint8_t));
    stats->groups_before_cull = ngroups;
    for (size_t g = 0; g < ngroups; g++) {
        stats->group_size_hist[aof_group_size_bin(nfaces[g])]++;
        if (max_faces > 0 && nfaces[g] <= max_faces) {
            kill[g] = 1;
            stats->groups_killed_small++;
        }
    }

    size_t killed = 0;
    for (size_t f = 0; f < ps->nf; f++) {
        if (!face_keep[f]) continue;
        int32_t c = vertex_chart[ps->faces[f * 3]];
        if (c < 0) continue;
        int32_t g = chart_group[c];
        if (g >= 0 && kill[g]) {
            face_keep[f] = 0;
            killed++;
        }
    }
    stats->faces_killed_groups = killed;
    return killed;
}

/* ========================================================================== */
/* Step 4: authentic neighbours -- chart pairs the XYZ weld would join         */
/* ========================================================================== */

/*
 * Emit one (lo, hi) chart-pair key per close cross-chart vertex pair whose
 * LOCAL normals agree.  Called twice: once with pair == NULL to size the
 * buffer, once to fill it.  Both passes walk identical logic, so the counts
 * agree by construction.
 *
 * The normal test is deliberately per-vertex rather than per-chart.  A chart
 * covers a whole 128-vox cube, so at radius r it subtends 128/r radians of its
 * wrap -- 73 degrees at r=100.  Two fragments of the SAME wrap in adjacent
 * cubes therefore have chart-mean normals that disagree by more than any
 * sensible threshold, and gating on those means rejects exactly the welds this
 * stage exists to find, worst where the scroll is tightest.  The normals of the
 * two vertices actually in contact do not have that problem.
 *
 * |dot| rather than dot: a fragment's winding may be flipped relative to its
 * neighbour, but a weld only cares that the two sheets are parallel.
 */
/*
 * One close cross-chart vertex pair that passed the normal test.  du is the
 * parameter mismatch measured AT the contact, oriented low chart minus high
 * chart.  Two welded vertices are the same physical point, so a correct
 * parameterization would put du at zero; whatever it actually is, is the
 * registration error, and radius converts it into turns.
 */
typedef struct {
    uint64_t key;      /* (lo << 32) | hi chart pair */
    double   du;
    double   dphi;     /* phi[lo] - phi[hi] at the same fixed XYZ contact */
    double   radius;
} AofContact;

static int aof_compare_contact(const void *pa, const void *pb)
{
    const AofContact *a = (const AofContact *)pa;
    const AofContact *b = (const AofContact *)pb;
    if (a->key != b->key) return a->key < b->key ? -1 : 1;
    /* Sorted within the run so the median is just the middle element. */
    return a->du < b->du ? -1 : (a->du > b->du ? 1 : 0);
}

static int aof_compare_double(const void *pa, const void *pb)
{
    double a = *(const double *)pa, b = *(const double *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static size_t aof_scan_close_pairs(const PieceSet *ps,
                                   const ScaffoldCalib *cal,
                                   const int32_t *vertex_chart,
                                   const int32_t *vertex_cube,
                                   const double *u,
                                   const AofRecord *record, size_t nrecords,
                                   double origin[3], double cell,
                                   double max_distance, double min_dot,
                                   AofContact *contact)
{
    double limit2 = max_distance * max_distance;
    size_t count = 0;
    for (size_t i = 0; i < ps->nv; i++) {
        int32_t chart_i = vertex_chart[i];
        if (chart_i < 0) continue;
        const float *pi = &ps->verts[i * 3];
        const float *ni = &ps->normals[i * 3];
        int32_t base[3];
        for (int k = 0; k < 3; k++)
            base[k] = aof_cell_of((double)pi[k], origin[k], cell);
        for (int dz = -1; dz <= 1; dz++)
        for (int dy = -1; dy <= 1; dy++)
        for (int dx = -1; dx <= 1; dx++) {
            uint64_t key = aof_key3(base[0] + dz, base[1] + dy, base[2] + dx);
            size_t at = aof_lower_bound(record, nrecords, key);
            for (; at < nrecords && record[at].key == key; at++) {
                size_t j = (size_t)record[at].item;
                if (j <= i) continue;
                int32_t chart_j = vertex_chart[j];
                if (chart_j < 0 || chart_j == chart_i) continue;
                /* Two disconnected charts emitted by the SAME cube cannot be
                 * opposite sides of a cube seam.  Treating their close,
                 * parallel surfaces as weld evidence is precisely how nearby
                 * scroll layers became one black/peach super-component.  A
                 * same-cube gap may still receive the weaker lateral coupling
                 * later; it must never define physical weld connectivity. */
                if (vertex_cube != NULL && vertex_cube[i] >= 0 &&
                    vertex_cube[i] == vertex_cube[j])
                    continue;
                const float *pj = &ps->verts[j * 3];
                double d0 = (double)pi[0] - (double)pj[0];
                double d1 = (double)pi[1] - (double)pj[1];
                double d2 = (double)pi[2] - (double)pj[2];
                if (d0 * d0 + d1 * d1 + d2 * d2 > limit2) continue;
                const float *nj = &ps->normals[j * 3];
                double dot = (double)ni[0] * (double)nj[0] +
                             (double)ni[1] * (double)nj[1] +
                             (double)ni[2] * (double)nj[2];
                /* A zero normal means undefined, not disagreeing; PieceSet
                 * leaves those on vertices with no kept faces. */
                if (fabs(dot) < min_dot &&
                    !(ni[0] == 0.0f && ni[1] == 0.0f && ni[2] == 0.0f) &&
                    !(nj[0] == 0.0f && nj[1] == 0.0f && nj[2] == 0.0f))
                    continue;
                if (contact != NULL) {
                    int lo_is_i = chart_i < chart_j;
                    int32_t lo = lo_is_i ? chart_i : chart_j;
                    int32_t hi = lo_is_i ? chart_j : chart_i;
                    contact[count].key = ((uint64_t)(uint32_t)lo << 32) |
                                         (uint64_t)(uint32_t)hi;
                    contact[count].du = lo_is_i ? u[i] - u[j] : u[j] - u[i];
                    contact[count].dphi =
                        ps->phi != NULL && isfinite((double)ps->phi[i]) &&
                        isfinite((double)ps->phi[j])
                        ? (lo_is_i ? (double)ps->phi[i] - (double)ps->phi[j]
                                   : (double)ps->phi[j] - (double)ps->phi[i])
                        : NAN;
                    contact[count].radius =
                        0.5 * (aof_radius(cal, pi) + aof_radius(cal, pj));
                }
                count++;
            }
        }
    }
    return count;
}

/*
 * One authentic weld between two charts, carrying its own measurement of how
 * far apart the two sides are in u and how many whole turns that is.
 */
typedef struct {
    int32_t chart0;        /* the lower chart id */
    int32_t chart1;
    size_t  contacts;
    double  du_median;     /* u[chart0] - u[chart1] at the contact */
    double  radius;        /* mean contact radius */
    double  circumference; /* 2*pi*radius: one turn of u here */
    int32_t wraps;         /* round(du_median / circumference) */
    double  residual;      /* |du/circ - wraps|; 0 = a clean whole turn */
    int32_t phase_target;   /* required outward depth k[chart0]-k[chart1],
                             * from phase at the fixed contact */
    double  phase_residual; /* distance of the phase relation from its nearest
                             * integer turn */
} AofNeighbour;

static int aof_build_neighbours(Arena_T arena, const PieceSet *ps,
                                const ScaffoldCalib *cal,
                                const int32_t *vertex_chart,
                                const double *u,
                                const AtlasOverlapFixOptions *opts,
                                AofNeighbour **out_edge, size_t *out_nedges)
{
    *out_edge = NULL;
    *out_nedges = 0;
    double cell = opts->neighbour_distance;
    if (!(cell > 0.0)) return 0;

    double origin[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
    for (size_t i = 0; i < ps->nv; i++) {
        if (vertex_chart[i] < 0) continue;
        for (int k = 0; k < 3; k++)
            if ((double)ps->verts[i * 3 + (size_t)k] < origin[k])
                origin[k] = (double)ps->verts[i * 3 + (size_t)k];
    }
    if (origin[0] == DBL_MAX) return 0;

    /* PieceSet loads every cube's vertices contiguously.  Preserve that
     * provenance through the proximity scan so a within-cube near miss cannot
     * masquerade as a cross-cube seam weld. */
    int32_t *vertex_cube = NULL;
    if (ps->cube_voff != NULL && ps->n_cubes > 0) {
        vertex_cube = (int32_t *)ARENA_ALLOC(
            arena, ps->nv * sizeof(*vertex_cube));
        for (size_t i = 0; i < ps->nv; i++) vertex_cube[i] = -1;
        for (size_t q = 0; q < ps->n_cubes; q++) {
            size_t first = ps->cube_voff[q];
            size_t last = ps->cube_voff[q + 1];
            if (first > last || last > ps->nv) return -1;
            for (size_t i = first; i < last; i++)
                vertex_cube[i] = (int32_t)q;
        }
    }

    size_t nrecords = 0;
    for (size_t i = 0; i < ps->nv; i++)
        if (vertex_chart[i] >= 0) nrecords++;
    if (nrecords == 0) return 0;

    AofRecord *record = (AofRecord *)ARENA_ALLOC(
        arena, nrecords * sizeof(*record));
    size_t at = 0;
    for (size_t i = 0; i < ps->nv; i++) {
        if (vertex_chart[i] < 0) continue;
        const float *p = &ps->verts[i * 3];
        record[at].key = aof_key3(aof_cell_of((double)p[0], origin[0], cell),
                                  aof_cell_of((double)p[1], origin[1], cell),
                                  aof_cell_of((double)p[2], origin[2], cell));
        record[at].item = (int32_t)i;
        at++;
    }
    qsort(record, nrecords, sizeof(*record), aof_compare_record);

    size_t npairs = aof_scan_close_pairs(
        ps, cal, vertex_chart, vertex_cube, u, record, nrecords, origin, cell,
        opts->neighbour_distance, opts->neighbour_normal_dot, NULL);
    if (npairs == 0) return 0;

    AofContact *contact = (AofContact *)ARENA_ALLOC(
        arena, npairs * sizeof(*contact));
    size_t filled = aof_scan_close_pairs(
        ps, cal, vertex_chart, vertex_cube, u, record, nrecords, origin, cell,
        opts->neighbour_distance, opts->neighbour_normal_dot, contact);
    if (filled != npairs) return -1;
    qsort(contact, npairs, sizeof(*contact), aof_compare_contact);

    /* Run-length encode.  Each surviving contact already passed the local
     * normal test, so all that remains is to require enough of them: one stray
     * close vertex is noise, a real seam contact produces many.  The median du
     * over the run is the edge's measurement -- median, not mean, because a
     * seam that grazes a second chart contributes outliers. */
    AofNeighbour *edge = (AofNeighbour *)ARENA_ALLOC(
        arena, npairs * sizeof(*edge));
    double *phase_sample = (double *)ARENA_ALLOC(
        arena, npairs * sizeof(*phase_sample));
    size_t nedges = 0;
    for (size_t first = 0; first < npairs; ) {
        size_t last = first + 1;
        while (last < npairs && contact[last].key == contact[first].key) last++;
        size_t run = last - first;
        if (run >= opts->neighbour_min_shared) {
            double radius = 0.0;
            for (size_t i = first; i < last; i++) radius += contact[i].radius;
            radius /= (double)run;
            double du = contact[first + run / 2].du;
            double circumference = AOF_2PI * radius;
            edge[nedges].chart0 = (int32_t)(uint32_t)(contact[first].key >> 32);
            edge[nedges].chart1 =
                (int32_t)(uint32_t)(contact[first].key & 0xFFFFFFFFu);
            edge[nedges].contacts = run;
            edge[nedges].du_median = du;
            edge[nedges].radius = radius;
            edge[nedges].circumference = circumference;
            if (circumference > 0.0) {
                double turns = du / circumference;
                double nearest = floor(turns + 0.5);
                edge[nedges].wraps = (int32_t)nearest;
                edge[nedges].residual = fabs(turns - nearest);
            } else {
                edge[nedges].wraps = 0;
                edge[nedges].residual = 0.0;
            }
            size_t nphase = 0;
            for (size_t i = first; i < last; i++)
                if (isfinite(contact[i].dphi))
                    phase_sample[nphase++] = contact[i].dphi;
            if (nphase > 0) {
                qsort(phase_sample, nphase, sizeof(*phase_sample),
                      aof_compare_double);
                double dphi = phase_sample[nphase / 2];
                double turns = -(double)cal->sense * dphi / AOF_2PI;
                double nearest = floor(turns + 0.5);
                edge[nedges].phase_target = (int32_t)nearest;
                edge[nedges].phase_residual = fabs(turns - nearest);
            } else {
                edge[nedges].phase_target = 0;
                edge[nedges].phase_residual = DBL_MAX;
            }
            nedges++;
        }
        first = last;
    }
    *out_edge = edge;
    *out_nedges = nedges;
    return 0;
}

/* ========================================================================== */
/* Step 5b: integrate the per-edge wind errors into a per-chart correction     */
/*                                                                            */
/* Each trusted weld says s[chart0] - s[chart1] = wraps * circumference, where */
/* s is the amount of u to remove from that chart.  Those are difference       */
/* equations on a graph, so a spanning tree fixes every chart relative to its  */
/* root by simple accumulation -- no solve.  Each step subtracts a whole       */
/* number of turns measured AT THAT CONTACT, so the correction stays an        */
/* integer-wrap move locally even though the circumference varies across the   */
/* group.                                                                     */
/*                                                                            */
/* Non-tree edges are then a free consistency check: if going around a cycle   */
/* does not return to the same wind, the weld graph disagrees with itself      */
/* there, which localizes a genuine fusion defect rather than a gauge error.   */
/* ========================================================================== */

typedef struct {
    size_t trusted_edges;
    size_t untrusted_edges;   /* residual too large to call it a whole turn */
    size_t oversized_edges;   /* |wraps| beyond the cap: not a gauge error */
    size_t tree_edges;
    size_t cycle_edges;
    size_t cycle_consistent;
    size_t cycle_inconsistent;
    size_t charts_corrected;
    size_t gauge_islands;     /* independent trees; > groups means cuts */
    int32_t wraps_min;
    int32_t wraps_max;
} AofWindStats;

static int aof_integrate_winds(Arena_T arena, size_t ncharts,
                               const AofNeighbour *edge, size_t nedges,
                               double residual_limit, int wrap_cap,
                               double *chart_shift, AofWindStats *stats)
{
    memset(stats, 0, sizeof(*stats));
    for (size_t c = 0; c < ncharts; c++) chart_shift[c] = 0.0;
    if (ncharts == 0) return 0;

    /* Trust gate first: an edge only carries gauge evidence if its measured
     * mismatch really is close to a whole number of turns. */
    uint8_t *trusted = (uint8_t *)ARENA_CALLOC(arena, nedges ? nedges : 1,
                                               sizeof(*trusted));
    for (size_t e = 0; e < nedges; e++) {
        if (edge[e].residual > residual_limit) {
            stats->untrusted_edges++;
        } else if (wrap_cap > 0 && abs(edge[e].wraps) > wrap_cap) {
            stats->oversized_edges++;
        } else {
            trusted[e] = 1;
            stats->trusted_edges++;
            if (edge[e].wraps < stats->wraps_min)
                stats->wraps_min = edge[e].wraps;
            if (edge[e].wraps > stats->wraps_max)
                stats->wraps_max = edge[e].wraps;
        }
    }

    /* CSR adjacency over the trusted edges. */
    size_t *row = (size_t *)ARENA_CALLOC(arena, ncharts + 1, sizeof(*row));
    for (size_t e = 0; e < nedges; e++) {
        if (!trusted[e]) continue;
        row[edge[e].chart0 + 1]++;
        row[edge[e].chart1 + 1]++;
    }
    for (size_t c = 0; c < ncharts; c++) row[c + 1] += row[c];
    size_t nadj = row[ncharts];
    int32_t *adj_edge = (int32_t *)ARENA_ALLOC(
        arena, (nadj ? nadj : 1) * sizeof(*adj_edge));
    size_t *fill = (size_t *)ARENA_ALLOC(arena, ncharts * sizeof(*fill));
    for (size_t c = 0; c < ncharts; c++) fill[c] = row[c];
    for (size_t e = 0; e < nedges; e++) {
        if (!trusted[e]) continue;
        adj_edge[fill[edge[e].chart0]++] = (int32_t)e;
        adj_edge[fill[edge[e].chart1]++] = (int32_t)e;
    }

    uint8_t *visited = (uint8_t *)ARENA_CALLOC(arena, ncharts,
                                               sizeof(*visited));
    int32_t *queue = (int32_t *)ARENA_ALLOC(arena, ncharts * sizeof(*queue));
    uint8_t *in_tree = (uint8_t *)ARENA_CALLOC(arena, nedges ? nedges : 1,
                                               sizeof(*in_tree));

    for (size_t root = 0; root < ncharts; root++) {
        if (visited[root]) continue;
        stats->gauge_islands++;
        size_t head = 0, tail = 0;
        visited[root] = 1;
        chart_shift[root] = 0.0;   /* the root keeps the place it already has */
        queue[tail++] = (int32_t)root;
        while (head < tail) {
            int32_t c = queue[head++];
            for (size_t a = row[c]; a < row[c + 1]; a++) {
                size_t e = (size_t)adj_edge[a];
                int32_t other = edge[e].chart0 == c ? edge[e].chart1
                                                    : edge[e].chart0;
                double step = (double)edge[e].wraps * edge[e].circumference;
                /* s[chart0] - s[chart1] = wraps * circumference */
                double target = edge[e].chart0 == c ? chart_shift[c] - step
                                                    : chart_shift[c] + step;
                if (!visited[other]) {
                    visited[other] = 1;
                    in_tree[e] = 1;
                    stats->tree_edges++;
                    chart_shift[other] = target;
                    queue[tail++] = other;
                }
            }
        }
    }

    for (size_t e = 0; e < nedges; e++) {
        if (!trusted[e] || in_tree[e]) continue;
        stats->cycle_edges++;
        double implied = chart_shift[edge[e].chart0] -
                         chart_shift[edge[e].chart1];
        double expected = (double)edge[e].wraps * edge[e].circumference;
        double slack = edge[e].circumference > 0.0
                     ? fabs(implied - expected) / edge[e].circumference : 0.0;
        if (slack < 0.5) stats->cycle_consistent++;
        else stats->cycle_inconsistent++;
    }
    for (size_t c = 0; c < ncharts; c++)
        if (chart_shift[c] != 0.0) stats->charts_corrected++;
    return 0;
}

/* ========================================================================== */
/* Step 5: groups = connected components of the authentic-neighbour graph      */
/* ========================================================================== */

static size_t aof_label_groups(Arena_T arena, size_t ncharts,
                               const AofNeighbour *edge, size_t nedges,
                               int32_t *chart_group)
{
    UnionFind uf = UF_new(arena, (int32_t)ncharts);
    for (size_t e = 0; e < nedges; e++)
        uf_union(&uf, edge[e].chart0, edge[e].chart1);
    int32_t *root_group = (int32_t *)ARENA_ALLOC(
        arena, ncharts * sizeof(int32_t));
    for (size_t i = 0; i < ncharts; i++) root_group[i] = -1;
    int32_t ngroups = 0;
    for (size_t chart = 0; chart < ncharts; chart++) {
        int32_t root = uf_find(&uf, (int32_t)chart);
        if (root_group[root] < 0) root_group[root] = ngroups++;
        chart_group[chart] = root_group[root];
    }
    return (size_t)ngroups;
}

/* ========================================================================== */
/* Step 6: overlap in (u, v)                                                   */
/* ========================================================================== */

/* Separating-axis test on two triangles.  Sharing an edge is NOT an overlap:
 * the comparisons are strict, so adjacent charts that merely touch pass. */
static int aof_triangles_overlap(const double *ua, const double *va,
                                 const double *ub, const double *vb)
{
    for (int t = 0; t < 2; t++) {
        const double *pu = t == 0 ? ua : ub;
        const double *pv = t == 0 ? va : vb;
        for (int e = 0; e < 3; e++) {
            int i = e, j = (e + 1) % 3;
            double nx = -(pv[j] - pv[i]);
            double ny = pu[j] - pu[i];
            if (fabs(nx) < 1.0e-15 && fabs(ny) < 1.0e-15) continue;
            double min_a = DBL_MAX, max_a = -DBL_MAX;
            double min_b = DBL_MAX, max_b = -DBL_MAX;
            for (int k = 0; k < 3; k++) {
                double da = nx * ua[k] + ny * va[k];
                double db = nx * ub[k] + ny * vb[k];
                if (da < min_a) min_a = da;
                if (da > max_a) max_a = da;
                if (db < min_b) min_b = db;
                if (db > max_b) max_b = db;
            }
            if (max_a <= min_b || max_b <= min_a) return 0;
        }
    }
    return 1;
}

typedef struct {
    const PieceSet *ps;
    const double *u;
    const double *v;
    const int32_t *face_group;
    AofRecord *record;
    size_t nrecords;
    double origin_u, origin_v;
    double cell;
} AofUvIndex;

static void aof_face_box(const PieceSet *ps, const double *u, const double *v,
                         size_t f, double *box)
{
    box[0] = DBL_MAX; box[1] = -DBL_MAX;
    box[2] = DBL_MAX; box[3] = -DBL_MAX;
    for (int k = 0; k < 3; k++) {
        int32_t vertex = ps->faces[f * 3 + (size_t)k];
        double fu = u[vertex], fv = v[vertex];
        if (fu < box[0]) box[0] = fu;
        if (fu > box[1]) box[1] = fu;
        if (fv < box[2]) box[2] = fv;
        if (fv > box[3]) box[3] = fv;
    }
}

/*
 * Index every kept face whose group is not `exclude`.  Cells are sized so a
 * typical face lands in a handful of them; a face wider than the clamp is
 * indexed by its corner cells only, which can only ever cost recall on
 * pathological geometry that step 2 has already mostly removed.
 */
enum { AOF_MAX_CELL_SPAN = 64 };

static int aof_uv_index_build(Arena_T arena, const PieceSet *ps,
                              const double *u, const double *v,
                              const uint8_t *face_keep,
                              const int32_t *face_group, int32_t exclude,
                              double cell_floor, AofUvIndex *index)
{
    memset(index, 0, sizeof(*index));
    index->ps = ps;
    index->u = u;
    index->v = v;
    index->face_group = face_group;

    double origin_u = DBL_MAX, origin_v = DBL_MAX;
    double extent = 0.0;
    size_t counted = 0;
    for (size_t f = 0; f < ps->nf; f++) {
        if (!face_keep[f] || face_group[f] < 0 || face_group[f] == exclude)
            continue;
        double box[4];
        aof_face_box(ps, u, v, f, box);
        if (box[0] < origin_u) origin_u = box[0];
        if (box[2] < origin_v) origin_v = box[2];
        extent += (box[1] - box[0]) + (box[3] - box[2]);
        counted++;
    }
    if (counted == 0) return 0;

    double cell = 0.5 * extent / (double)counted;
    if (!(cell > cell_floor)) cell = cell_floor;
    if (!(cell > 0.0)) cell = 1.0;
    index->origin_u = origin_u;
    index->origin_v = origin_v;
    index->cell = cell;

    size_t nrecords = 0;
    for (int pass = 0; pass < 2; pass++) {
        size_t at = 0;
        for (size_t f = 0; f < ps->nf; f++) {
            if (!face_keep[f] || face_group[f] < 0 || face_group[f] == exclude)
                continue;
            double box[4];
            aof_face_box(ps, u, v, f, box);
            int32_t cu0 = aof_cell_of(box[0], origin_u, cell);
            int32_t cu1 = aof_cell_of(box[1], origin_u, cell);
            int32_t cv0 = aof_cell_of(box[2], origin_v, cell);
            int32_t cv1 = aof_cell_of(box[3], origin_v, cell);
            if (cu1 - cu0 > AOF_MAX_CELL_SPAN) cu1 = cu0 + AOF_MAX_CELL_SPAN;
            if (cv1 - cv0 > AOF_MAX_CELL_SPAN) cv1 = cv0 + AOF_MAX_CELL_SPAN;
            for (int32_t cu = cu0; cu <= cu1; cu++)
                for (int32_t cv = cv0; cv <= cv1; cv++) {
                    if (pass == 1) {
                        index->record[at].key = aof_key2(cu, cv);
                        index->record[at].item = (int32_t)f;
                    }
                    at++;
                }
        }
        if (pass == 0) {
            nrecords = at;
            if (nrecords == 0) return 0;
            index->record = (AofRecord *)ARENA_ALLOC(
                arena, nrecords * sizeof(AofRecord));
        } else if (at != nrecords) {
            return -1;
        }
    }
    index->nrecords = nrecords;
    qsort(index->record, nrecords, sizeof(AofRecord), aof_compare_record);
    return 0;
}

/*
 * Count indexed faces that overlap face f, whose vertices are read at
 * (u, v) plus (shift_u, shift_v).
 *
 * A pair that shares several broad-phase cells must be reported ONCE.  The
 * cell holding the lower corner of the two AABBs' intersection is the unique
 * cell both faces index and both cell walks visit, so that is the one cell
 * allowed to report the pair.  Without this every count is inflated by the
 * number of shared cells, which is silent and scales with cell size.
 *
 * only_greater is for symmetric self-scans, where the index also contains f
 * and each unordered pair would otherwise be seen from both ends.
 * stop_at_first short-circuits the placement search, which only needs to know
 * whether the count is zero.
 */
typedef void (*AofPairFn)(void *context, size_t face_a, size_t face_b);

static size_t aof_query_face(const AofUvIndex *index, size_t f,
                             double shift_u, double shift_v,
                             int only_greater, int stop_at_first,
                             AofPairFn fn, void *context)
{
    if (index->nrecords == 0) return 0;
    const PieceSet *ps = index->ps;
    double ua[3], va[3];
    double box[4] = {DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX};
    for (int k = 0; k < 3; k++) {
        int32_t vertex = ps->faces[f * 3 + (size_t)k];
        ua[k] = index->u[vertex] + shift_u;
        va[k] = index->v[vertex] + shift_v;
        if (ua[k] < box[0]) box[0] = ua[k];
        if (ua[k] > box[1]) box[1] = ua[k];
        if (va[k] < box[2]) box[2] = va[k];
        if (va[k] > box[3]) box[3] = va[k];
    }
    int32_t cu0 = aof_cell_of(box[0], index->origin_u, index->cell);
    int32_t cu1 = aof_cell_of(box[1], index->origin_u, index->cell);
    int32_t cv0 = aof_cell_of(box[2], index->origin_v, index->cell);
    int32_t cv1 = aof_cell_of(box[3], index->origin_v, index->cell);
    if (cu1 - cu0 > AOF_MAX_CELL_SPAN) cu1 = cu0 + AOF_MAX_CELL_SPAN;
    if (cv1 - cv0 > AOF_MAX_CELL_SPAN) cv1 = cv0 + AOF_MAX_CELL_SPAN;

    size_t hits = 0;
    for (int32_t cu = cu0; cu <= cu1; cu++)
    for (int32_t cv = cv0; cv <= cv1; cv++) {
        uint64_t key = aof_key2(cu, cv);
        size_t at = aof_lower_bound(index->record, index->nrecords, key);
        for (; at < index->nrecords && index->record[at].key == key; at++) {
            size_t g = (size_t)index->record[at].item;
            if (g == f || (only_greater && g <= f)) continue;
            double ub[3], vb[3];
            double gbox[4] = {DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX};
            for (int k = 0; k < 3; k++) {
                int32_t vertex = ps->faces[g * 3 + (size_t)k];
                ub[k] = index->u[vertex];
                vb[k] = index->v[vertex];
                if (ub[k] < gbox[0]) gbox[0] = ub[k];
                if (ub[k] > gbox[1]) gbox[1] = ub[k];
                if (vb[k] < gbox[2]) gbox[2] = vb[k];
                if (vb[k] > gbox[3]) gbox[3] = vb[k];
            }
            if (box[1] <= gbox[0] || gbox[1] <= box[0] ||
                box[3] <= gbox[2] || gbox[3] <= box[2])
                continue;
            /* Report only from the intersection's lower-corner cell. */
            double lo_u = box[0] > gbox[0] ? box[0] : gbox[0];
            double lo_v = box[2] > gbox[2] ? box[2] : gbox[2];
            if (aof_cell_of(lo_u, index->origin_u, index->cell) != cu ||
                aof_cell_of(lo_v, index->origin_v, index->cell) != cv)
                continue;
            if (!aof_triangles_overlap(ua, va, ub, vb)) continue;
            hits++;
            if (fn != NULL) fn(context, f, g);
            if (stop_at_first) return hits;
        }
    }
    return hits;
}

/* Attributes each reported pair to its group(s).  Pairs inside one group are
 * counted separately: a rigid group shift cannot separate them (only the
 * per-chart tabu re-gauge can). */
typedef struct {
    const int32_t *faces;
    const int32_t *face_group;
    const int32_t *vertex_chart;
    const double *face_radius;   /* optional; enables the radial histogram */
    double pitch;
    AtlasOverlapFixGroup *group; /* optional; NULL counts without attribution */
    size_t cross_pairs;
    size_t intra_pairs;
    size_t *radial_hist;         /* optional, 6 bins on |dr| / pitch */
    int    cross_only;           /* 1 = do not count intra-group pairs */
    int    after;                /* attribute cross pairs to overlap_after */
} AofCountContext;

/*
 * Bin an overlapping pair by how far apart the two faces are RADIALLY, in
 * units of the wrap pitch.  This separates the two ways an atlas can overlap:
 * |dr| ~ 0 means one sheet doubled onto itself, while |dr| ~ 1 pitch means two
 * DIFFERENT wraps landed on the same u -- the signature of a parameterization
 * that fails to advance a full circumference per turn.
 */
static void aof_bin_radial(AofCountContext *c, size_t face_a, size_t face_b)
{
    if (c->radial_hist == NULL || c->face_radius == NULL ||
        !(c->pitch > 0.0))
        return;
    double dr = fabs(c->face_radius[face_a] - c->face_radius[face_b]) /
                c->pitch;
    size_t bin = dr < 0.25 ? 0 : dr < 0.75 ? 1 : dr < 1.25 ? 2
               : dr < 1.75 ? 3 : dr < 2.5 ? 4 : 5;
    c->radial_hist[bin]++;
}

static void aof_count_pair(void *context, size_t face_a, size_t face_b)
{
    AofCountContext *c = (AofCountContext *)context;
    int32_t ga = c->face_group[face_a], gb = c->face_group[face_b];
    if (ga < 0 || gb < 0) return;
    if (ga == gb) {
        if (c->cross_only) return;
        int32_t chart_a = c->vertex_chart[c->faces[face_a * 3]];
        int32_t chart_b = c->vertex_chart[c->faces[face_b * 3]];
        if (chart_a == chart_b) return;   /* one chart folded on itself */
        c->intra_pairs++;
        if (c->group != NULL && !c->after) {
            c->group[ga].intra_overlap++;
            if (c->face_radius != NULL && c->pitch > 0.0 &&
                fabs(c->face_radius[face_a] - c->face_radius[face_b]) >=
                    0.75 * c->pitch)
                c->group[ga].intra_wrongwind++;
        }
        aof_bin_radial(c, face_a, face_b);
        return;
    }
    c->cross_pairs++;
    if (c->group == NULL) return;
    if (c->after) {
        c->group[ga].overlap_after++;
        c->group[gb].overlap_after++;
    } else {
        c->group[ga].overlap_before++;
        c->group[gb].overlap_before++;
    }
}

/* One symmetric pass over every kept face against an index of all of them. */
/* Wall-clock origin for the phase log, reset per solve. */
static double g_aof_t0 = 0.0;

/* Phase marker. Cheap and few, and the only way to tell a slow exact
 * scan from a hang: on a 1000-cube grid the solve otherwise prints one
 * line per round, and that line can be an hour away. */
static void aof_phase(const char *what)
{
    fprintf(stderr, "[overlap_fix]   %-26s t=%.1fs\n", what,
            ves_clock_sec() - g_aof_t0);
    fflush(stderr);
}

static size_t aof_scan_all(const PieceSet *ps, const uint8_t *face_keep,
                           const int32_t *face_group, const AofUvIndex *index,
                           AofCountContext *context)
{
    if (index->nrecords == 0) return 0;
    size_t total = 0;
    /* This scan is O(faces x candidates) and dominates at scale;
     * without a heartbeat it is indistinguishable from a deadlock. */
    const size_t aof_tick = 100000;
    double aof_t_scan = ves_clock_sec();
    for (size_t f = 0; f < ps->nf; f++) {
        if (f && (f % aof_tick) == 0) {
            double dt = ves_clock_sec() - aof_t_scan;
            fprintf(stderr,
                    "[overlap_fix]     scan %zu/%zu faces, %zu pairs, "
                    "%.1fs (%.0f faces/s)\n",
                    f, ps->nf, total, dt, dt > 0.0 ? (double)f / dt : 0.0);
            fflush(stderr);
        }
        if (!face_keep[f] || face_group[f] < 0) continue;
        total += aof_query_face(index, f, 0.0, 0.0, 1, 0,
                                aof_count_pair, context);
    }
    return total;
}

/* ========================================================================== */
/* Step 7: placement order -- furthest from the axis first                     */
/* ========================================================================== */

typedef struct {
    int32_t group;
    double  radius;
} AofOrderEntry;

static int aof_compare_order(const void *pa, const void *pb)
{
    const AofOrderEntry *a = (const AofOrderEntry *)pa;
    const AofOrderEntry *b = (const AofOrderEntry *)pb;
    if (a->radius > b->radius) return -1;
    if (a->radius < b->radius) return 1;
    return a->group < b->group ? -1 : (a->group > b->group ? 1 : 0);
}

/* ========================================================================== */
/* Tabu re-gauge: one integer winding depth per chart                          */
/*                                                                            */
/* The rigid group placement can fix a group that is WHOLLY mis-wound, but    */
/* not a group that disagrees with itself: a connected spiral whose u under-  */
/* advances per turn lands wrap k+1 back on wrap k, and every rigid move      */
/* preserves that.  The intra-group radial histogram says ~85% of the         */
/* residual overlap has exactly that wrong-wrap signature.  So the discrete   */
/* repair is per chart: every chart gets its own integer winding depth, welds */
/* may tear where the accumulated drift demands it, and each tear is paid     */
/* for explicitly:                                                            */
/*                                                                            */
/*   E = Lov * (occupancy cells shared between charts)          [collision]   */
/*     + sum over usable weld edges  w_e * |k[a]-k[b]-target|   [neighbour]   */
/*     + sum over charts  radE[c][k[c]]                         [radial]      */
/*                                                                            */
/* Everything is fixed-point integer (x1024) so incremental deltas cannot     */
/* drift and the search is bit-reproducible: charts scanned ascending, depths */
/* ascending, strictly-better wins, no RNG anywhere.                          */
/*                                                                            */
/* Collision is a cell proxy, not exact SAT: each chart rasterizes its kept   */
/* faces once into cell-sized (u, v) buckets; a depth move is an integer add  */
/* on the packed key.  The proxy never reaches zero (touching charts share a  */
/* boundary strip), so termination is stagnation of the best energy, and      */
/* exact SAT picks among the best assignments at the end -- the place-search  */
/* lesson: the exact test is the final verifier, never the inner loop.  At    */
/* ~10x the chart count the full move sweep goes slow; the fix then is a      */
/* per-chart candidate-vector cache invalidated by u-interval disturbance,    */
/* which is why evaluation is factored per chart.                             */
/* ========================================================================== */

enum {
    AOF_TABU_EPOCH = 50,          /* iterations between happiness refreshes */
    AOF_TABU_AUDIT_EVERY = 256,   /* incremental-vs-recount self check */
    AOF_TABU_VERIFY_MAX = 64,     /* diverse archived states SAT-verified */
    AOF_TABU_CONTACT_CAP = 32     /* weld contacts counted per edge, max */
};

#define AOF_TABU_FIX 1024.0      /* fixed-point scale of the energy */

/* Open-addressing occupancy counts.  Keys are (v cell, u cell) packed with a
 * 2^30 bias per half so a real key is never 0 and a u shift can never carry
 * into the v half; capacity is a power of two and live load stays under 1/2
 * (a rebuild reclaims slots stranded at count 0). */
typedef struct {
    uint64_t *key;
    uint32_t *count;
    size_t cap;
    size_t used;
} AofTabuHash;

typedef struct {
    int32_t face;
    double umin, umax, vmin, vmax;
} AofTabuFaceBox;

/* Memoized exact face-pair counts for one (chart, depth) placement pair. */
typedef struct {
    uint64_t *key;
    uint32_t *count;
    size_t cap;
    size_t used;
} AofTabuExactHash;

static uint64_t aof_tabu_mix(uint64_t x)
{
    x ^= x >> 30; x *= 0xBF58476D1CE4E5B9ull;
    x ^= x >> 27; x *= 0x94D049BB133111EBull;
    x ^= x >> 31;
    return x;
}

static size_t aof_tabu_slot(const AofTabuHash *h, uint64_t key)
{
    size_t mask = h->cap - 1;
    size_t s = (size_t)aof_tabu_mix(key) & mask;
    while (h->key[s] != 0 && h->key[s] != key) s = (s + 1) & mask;
    return s;
}

static uint64_t aof_tabu_cell_key(double u, double v, double cell)
{
    double cu = floor(u / cell), cv = floor(v / cell);
    double lim = (double)(1 << 29);
    if (cu < -lim) cu = -lim;
    if (cu > lim) cu = lim;
    if (cv < -lim) cv = -lim;
    if (cv > lim) cv = lim;
    return ((uint64_t)(uint32_t)((int32_t)cv + (1 << 30)) << 32) |
           (uint64_t)(uint32_t)((int32_t)cu + (1 << 30));
}

typedef struct {
    int32_t a, b;        /* chart ids, a < b */
    int32_t target;      /* required k[a] - k[b] for the weld to hold */
    int64_t w;           /* canonical weight, fixed point */
    int64_t w_unhappy;   /* steering weight when the partner is colliding */
    uint8_t family;      /* both endpoints in one happy family (boosted) */
    uint8_t phase_trusted; /* lateral target has low-residual source support */
    double  du;          /* pre-tabu u[a] - u[b]: weld contact median, or the
                          * lateral pair's chart-median gap (relayout base) */
} AofTabuEdge;

typedef struct {
    size_t ncharts, ngroups;
    int span, nk;
    int64_t Lov;              /* lambda_ov, fixed point, per exact face pair */
    double cell;
    Arena_T arena;

    /* The FIELD's wind ladder, measured from same-cube radial chart pairs.
     * ladder_step is the signed du per (2*pi*r) of one wind outward: this
     * field's u is not obliged to follow the spiral map's sign or scale
     * (measured here: u GROWS outward while the spiral u shrinks), and using
     * the spiral map unmeasured is exactly the bug that sent every re-gauged
     * chart one circumference the WRONG WAY. */
    int    ladder_measured;
    int    ladder_fixed;
    size_t ladder_pairs;
    double ladder_step;
    double ladder_u0;
    double dr_dw;             /* radius growth per wind = spiral_b * sense */
    size_t  *crow;            /* [n_cubes + 1] cube -> chart CSR (may be NULL) */
    int32_t *cmember;

    const int32_t *chart_group;
    int32_t *chart_cube;      /* [ncharts] index into ps->ids, -1 unknown */
    uint8_t *active;          /* [ncharts] movable core; halo charts are fixed */
    int8_t  *kmin;            /* [ncharts] lowest admissible depth */
    double  *u_med, *r_med;   /* [ncharts] medians */
    double  *w_geo;           /* [ncharts] radius-implied wind */
    double  *w_phi0, *w_rad;  /* [ncharts] winds in turns; w_rad recentred */
    double  *w_phase;         /* [ncharts] source phase after integer gauge */
    int32_t *phase_gauge;     /* [ncharts] source-phase graph potential */
    int32_t *phase_island;    /* [ncharts] gauge island id, -1 = no phase */
    double  *wind;            /* [ncharts] absolute continuous wind, turns */
    int8_t  *phase_bin;       /* [ncharts] median-centred missing wind */
    uint8_t *phase_trusted_group; /* [ngroups] complete, connected, consistent */
    uint8_t *phase_flat_group;/* [ngroups] trusted lift has one relative bin */
    size_t phase_gauge_islands;
    size_t phase_gauge_inconsistent;
    double  *shift;           /* [ncharts * nk] u shift per depth */
    int64_t *cell_shift;      /* [ncharts * nk] packed-key delta per depth */
    int64_t *radE;            /* [ncharts * nk] radial + gauge prior, fixed pt */

    double  *w_pool;          /* [ncharts] pooled u-wind the shifts step from */

    /* Graded rewind support: a group whose winds COLLAPSED onto each other
     * (wrong-wind intra pairs) needs wind w moved by k and wind w+1 by k+1 --
     * a per-chart search only ever finds that through uphill cascades and
     * usually scatters instead.  gbin is the member's geometric wind relative
     * to its group's innermost wind (lateral-cluster-median radius, so the
     * axis wander is voted away); the graded move assigns k = d + gbin. */
    int8_t  *gbin;            /* [ncharts] wind bin within the group */
    uint8_t *graded_group;    /* [ngroups] group is graded-move eligible */


    /* Strong cross-group phase links define a relative-gauge region.  Every
     * group retains region_offset + gbin while the region searches one shared
     * base depth.  Groups outside a region keep their ordinary move classes. */
    size_t n_region;
    size_t n_region_links;
    int32_t *region_id;       /* [ngroups], -1 outside an active region */
    int32_t *region_offset;   /* [ngroups], centred relative group base */
    size_t *region_row;       /* [n_region + 1], region -> group CSR */
    int32_t *region_group;

    /* Exact collision engine.  A chart pair is AABB/sweep-pruned, then every
     * surviving face pair is decided by triangle SAT.  Results are cached by
     * both integer depths, making repeat tabu evaluations hash lookups. */
    const PieceSet *ps;
    const double *base_u, *base_v;
    const int32_t *vertex_chart;
    size_t *frow;              /* [ncharts + 1], sorted face-box CSR */
    AofTabuFaceBox *fbox;
    double *chart_box;         /* [ncharts * 4]: umin,umax,vmin,vmax */
    double *chart_face_umax;   /* max face u width per chart */
    /* Every chart pair whose AABBs can overlap at ANY admissible pair of
     * integer depths.  Exact SAT only needs this sparse, provably complete
     * graph; iterating/allocating ncharts^2 made large grids impossible. */
    size_t *prow;              /* [ncharts + 1] potential-pair CSR */
    int32_t *pmember;          /* [2 * npotential_pairs] */
    size_t npotential_pairs;
    AofTabuExactHash exact;
    size_t exact_queries;
    size_t exact_cache_hits;
    size_t exact_cache_misses;
    size_t exact_aabb_rejects;
    size_t exact_sat_tests;

    /* Legacy sparse occupancy remains only as a proposal diagnostic while the
     * exact engine owns acceptance and the canonical collision energy. */
    uint64_t *ckey;           /* base cells, sorted + deduped per chart */
    size_t  *coff;            /* [ncharts + 1] */

    AofTabuEdge *edge;
    size_t nedges;
    size_t first_lateral;     /* edges [first_lateral..nedges) are lateral */
    size_t phase_lateral_trusted;
    size_t phase_lateral_nonzero;
    size_t  *erow;            /* [ncharts + 1] incident-edge CSR */
    int32_t *eadj;

    size_t  *grow;            /* [ngroups + 1] member CSR */
    int32_t *gmember;

    /* Neighbourhoods: connected components of the weld + lateral graph.
     * The compound move class that carries a coherent jigsaw region as one
     * unit -- no edge crosses a neighbourhood boundary, so these moves have
     * exactly zero neighbour cost, where a chart-by-chart cascade would pay
     * a lateral tear at every step and stall. */
    size_t n_nbhd;
    int32_t *nbid;            /* [ncharts] */
    size_t  *nbrow;           /* [n_nbhd + 1] member CSR */
    int32_t *nbmember;

    AofTabuHash hash;
    int8_t  *k;               /* [ncharts] current depth */
    uint8_t *happy;           /* [ncharts] collision-free right now */
    int64_t e_ov;             /* exact overlapping face pairs, unscaled */
    int64_t e_nb;             /* fixed point */
    int64_t e_rad;            /* fixed point */
    size_t rebuilds;
} AofTabu;

static int aof_tabu_compare_face_box(const void *pa, const void *pb)
{
    const AofTabuFaceBox *a = (const AofTabuFaceBox *)pa;
    const AofTabuFaceBox *b = (const AofTabuFaceBox *)pb;
    if (a->umin < b->umin) return -1;
    if (a->umin > b->umin) return 1;
    return a->face < b->face ? -1 : (a->face > b->face ? 1 : 0);
}

static size_t aof_tabu_face_lower_bound(const AofTabuFaceBox *box,
                                        size_t lo, size_t hi, double umin)
{
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (box[mid].umin < umin) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

static double aof_tabu_shift_at(const AofTabu *T, int32_t c, int k)
{
    if (k == 0) return 0.0;
    return T->shift[(size_t)c * (size_t)T->nk +
                    (size_t)(k + T->span)];
}

typedef struct {
    double umin, umax, vmin, vmax;
    int vb0, vb1;
} AofTabuPotentialBox;

typedef struct {
    double umin;
    int32_t chart;
} AofTabuPotentialOrder;

static int aof_tabu_compare_potential_order(const void *pa, const void *pb)
{
    const AofTabuPotentialOrder *a = (const AofTabuPotentialOrder *)pa;
    const AofTabuPotentialOrder *b = (const AofTabuPotentialOrder *)pb;
    if (a->umin < b->umin) return -1;
    if (a->umin > b->umin) return 1;
    return a->chart < b->chart ? -1 : (a->chart > b->chart ? 1 : 0);
}

/* Sweep expanded chart boxes in u, with a few bounded v bins keeping the
 * active set sparse.  A pair is emitted in exactly one shared v bin: the bin
 * containing max(vmin[a],vmin[b]).  `degree` selects the count pass;
 * `cursor/member` select the fill pass. */
static int aof_tabu_potential_sweep(
    Arena_T arena, const AofTabu *T, const AofTabuPotentialBox *box,
    const AofTabuPotentialOrder *order, size_t nvbins,
    size_t *degree, size_t *cursor, int32_t *member)
{
    size_t ncharts = T->ncharts;
    if (nvbins == 0 || nvbins > SIZE_MAX / ncharts) return -1;
    size_t *bin_count = (size_t *)ARENA_CALLOC(
        arena, nvbins, sizeof(*bin_count));
    int32_t *bin_member = (int32_t *)ARENA_ALLOC(
        arena, nvbins * ncharts * sizeof(*bin_member));

    for (size_t oi = 0; oi < ncharts; oi++) {
        int32_t c = order[oi].chart;
        const AofTabuPotentialBox *bc = &box[c];
        for (int vb = bc->vb0; vb <= bc->vb1; vb++) {
            size_t base = (size_t)vb * ncharts;
            size_t keep = 0;
            size_t count = bin_count[vb];
            for (size_t q = 0; q < count; q++) {
                int32_t other = bin_member[base + q];
                const AofTabuPotentialBox *bo = &box[other];
                /* Expired entries are compacted away while the list is hot. */
                if (bo->umax <= bc->umin) continue;
                bin_member[base + keep++] = other;

                if (bo->vmax <= bc->vmin || bc->vmax <= bo->vmin) continue;
                int canonical_vb = bo->vmin > bc->vmin ? bo->vb0 : bc->vb0;
                if (vb != canonical_vb) continue;

                if (degree != NULL) {
                    if (degree[c] == SIZE_MAX || degree[other] == SIZE_MAX)
                        return -1;
                    degree[c]++;
                    degree[other]++;
                } else {
                    if (cursor == NULL || member == NULL ||
                        cursor[c] >= T->prow[(size_t)c + 1] ||
                        cursor[other] >= T->prow[(size_t)other + 1])
                        return -1;
                    member[cursor[c]++] = other;
                    member[cursor[other]++] = c;
                }
            }
            bin_count[vb] = keep;
        }
        for (int vb = bc->vb0; vb <= bc->vb1; vb++) {
            size_t at = (size_t)vb * ncharts + bin_count[vb]++;
            bin_member[at] = c;
        }
    }
    return 0;
}

/* Build the exact engine's sparse universe.  Each base chart box is expanded
 * by the minimum/maximum u shift over every admissible integer depth.  Boxes
 * that do not intersect after that expansion can never collide in any state,
 * so dropping those pairs is exact, not an approximation. */
static int aof_tabu_build_potential_pairs(Arena_T arena, AofTabu *T)
{
    size_t ncharts = T->ncharts;
    AofTabuPotentialBox *box = (AofTabuPotentialBox *)ARENA_ALLOC(
        arena, ncharts * sizeof(*box));
    AofTabuPotentialOrder *order = (AofTabuPotentialOrder *)ARENA_ALLOC(
        arena, ncharts * sizeof(*order));
    double v_origin = DBL_MAX, v_limit = -DBL_MAX;

    for (size_t c = 0; c < ncharts; c++) {
        double shift_min = DBL_MAX, shift_max = -DBL_MAX;
        for (int k = T->kmin[c]; k <= T->span; k++) {
            double s = aof_tabu_shift_at(T, (int32_t)c, k);
            if (s < shift_min) shift_min = s;
            if (s > shift_max) shift_max = s;
        }
        const double *cb = &T->chart_box[c * 4];
        if (!isfinite(shift_min) || !isfinite(shift_max) ||
            !isfinite(cb[0]) || !isfinite(cb[1]) ||
            !isfinite(cb[2]) || !isfinite(cb[3]) ||
            cb[1] < cb[0] || cb[3] < cb[2])
            return -1;
        box[c].umin = cb[0] + shift_min;
        box[c].umax = cb[1] + shift_max;
        box[c].vmin = cb[2];
        box[c].vmax = cb[3];
        order[c].umin = box[c].umin;
        order[c].chart = (int32_t)c;
        if (box[c].vmin < v_origin) v_origin = box[c].vmin;
        if (box[c].vmax > v_limit) v_limit = box[c].vmax;
    }
    qsort(order, ncharts, sizeof(*order), aof_tabu_compare_potential_order);

    double v_range = v_limit - v_origin;
    size_t nvbins = v_range > 0.0 ? (size_t)ceil(v_range / 64.0) : 1;
    if (nvbins < 1) nvbins = 1;
    if (nvbins > 64) nvbins = 64;
    double v_cell = v_range > 0.0 ? v_range / (double)nvbins : 1.0;
    if (!(v_cell > 0.0) || !isfinite(v_cell)) return -1;

    for (size_t c = 0; c < ncharts; c++) {
        double vend = nextafter(box[c].vmax, -DBL_MAX);
        if (vend < box[c].vmin) vend = box[c].vmin;
        int b0 = (int)floor((box[c].vmin - v_origin) / v_cell);
        int b1 = (int)floor((vend - v_origin) / v_cell);
        if (b0 < 0) b0 = 0;
        if (b1 < b0) b1 = b0;
        if (b0 >= (int)nvbins) b0 = (int)nvbins - 1;
        if (b1 >= (int)nvbins) b1 = (int)nvbins - 1;
        box[c].vb0 = b0;
        box[c].vb1 = b1;
    }

    T->prow = (size_t *)ARENA_CALLOC(
        arena, ncharts + 1, sizeof(*T->prow));
    Arena_Mark mark = Arena_save(arena);
    if (aof_tabu_potential_sweep(arena, T, box, order, nvbins,
                                 &T->prow[1], NULL, NULL) != 0) {
        Arena_restore(arena, mark);
        return -1;
    }
    Arena_restore(arena, mark);

    for (size_t c = 0; c < ncharts; c++) {
        if (T->prow[c + 1] > SIZE_MAX - T->prow[c]) return -1;
        T->prow[c + 1] += T->prow[c];
    }
    if ((T->prow[ncharts] & 1u) != 0) return -1;
    T->npotential_pairs = T->prow[ncharts] / 2;
    T->pmember = (int32_t *)ARENA_ALLOC(
        arena, (T->prow[ncharts] ? T->prow[ncharts] : 1) *
               sizeof(*T->pmember));
    size_t *cursor = (size_t *)ARENA_ALLOC(
        arena, ncharts * sizeof(*cursor));
    memcpy(cursor, T->prow, ncharts * sizeof(*cursor));
    mark = Arena_save(arena);
    if (aof_tabu_potential_sweep(arena, T, box, order, nvbins,
                                 NULL, cursor, T->pmember) != 0) {
        Arena_restore(arena, mark);
        return -1;
    }
    Arena_restore(arena, mark);
    for (size_t c = 0; c < ncharts; c++)
        if (cursor[c] != T->prow[c + 1]) return -1;

    long double dense = ncharts > 1
        ? (long double)ncharts * (long double)(ncharts - 1) / 2.0L : 1.0L;
    fprintf(stderr,
            "[overlap_fix]     exact potential graph: %zu pairs "
            "(%.4Lf%% dense), %zu v bins\n",
            T->npotential_pairs,
            100.0L * (long double)T->npotential_pairs / dense,
            nvbins);
    fflush(stderr);
    return 0;
}

static size_t aof_tabu_exact_slot(const AofTabuExactHash *h, uint64_t key)
{
    size_t mask = h->cap - 1;
    size_t s = (size_t)aof_tabu_mix(key) & mask;
    while (h->key[s] != 0 && h->key[s] != key) s = (s + 1) & mask;
    return s;
}

static void aof_tabu_exact_grow(AofTabu *T)
{
    AofTabuExactHash *h = &T->exact;
    size_t old_cap = h->cap;
    uint64_t *old_key = h->key;
    uint32_t *old_count = h->count;
    size_t cap = old_cap ? old_cap * 2 : 64;
    h->key = (uint64_t *)ARENA_CALLOC(T->arena, cap, sizeof(uint64_t));
    h->count = (uint32_t *)ARENA_CALLOC(T->arena, cap, sizeof(uint32_t));
    h->cap = cap;
    h->used = 0;
    for (size_t i = 0; i < old_cap; i++) {
        if (old_key[i] == 0) continue;
        size_t s = aof_tabu_exact_slot(h, old_key[i]);
        h->key[s] = old_key[i];
        h->count[s] = old_count[i];
        h->used++;
    }
}

static uint64_t aof_tabu_exact_key(const AofTabu *T, int32_t ca, int ka,
                                   int32_t cb, int kb)
{
    if (ca > cb) {
        int32_t tc = ca; ca = cb; cb = tc;
        int tk = ka; ka = kb; kb = tk;
    }
    uint64_t key = (uint64_t)(uint32_t)ca;
    key = key * (uint64_t)T->ncharts + (uint64_t)(uint32_t)cb;
    key = key * (uint64_t)T->nk + (uint64_t)(ka + T->span);
    key = key * (uint64_t)T->nk + (uint64_t)(kb + T->span);
    return key + 1;
}

/* Exact overlaps between two charts at two integer placements.  The chart and
 * per-face boxes reject almost everything; SAT is called only for surviving
 * face pairs. */
static int aof_tabu_exact_aabb_separate(const AofTabu *T, int32_t ca, int ka,
                                        int32_t cb, int kb)
{
    double sa = aof_tabu_shift_at(T, ca, ka);
    double sb = aof_tabu_shift_at(T, cb, kb);
    const double *ba = &T->chart_box[(size_t)ca * 4];
    const double *bb = &T->chart_box[(size_t)cb * 4];
    return ba[1] + sa <= bb[0] + sb || bb[1] + sb <= ba[0] + sa ||
           ba[3] <= bb[2] || bb[3] <= ba[2];
}

static uint32_t aof_tabu_exact_pair_uncached(AofTabu *T, int32_t ca, int ka,
                                             int32_t cb, int kb)
{
    double sa = aof_tabu_shift_at(T, ca, ka);
    double sb = aof_tabu_shift_at(T, cb, kb);
    size_t na = T->frow[(size_t)ca + 1] - T->frow[ca];
    size_t nb = T->frow[(size_t)cb + 1] - T->frow[cb];
    if (na > nb) {
        int32_t tc = ca; ca = cb; cb = tc;
        int tk = ka; ka = kb; kb = tk;
        double ts = sa; sa = sb; sb = ts;
    }

    const PieceSet *ps = T->ps;
    uint32_t count = 0;
    size_t a0 = T->frow[ca], a1 = T->frow[(size_t)ca + 1];
    size_t b0 = T->frow[cb], b1 = T->frow[(size_t)cb + 1];
    double bwidth = T->chart_face_umax[cb];
    for (size_t ia = a0; ia < a1; ia++) {
        const AofTabuFaceBox *fa = &T->fbox[ia];
        double au0 = fa->umin + sa, au1 = fa->umax + sa;
        size_t lo = aof_tabu_face_lower_bound(
            T->fbox, b0, b1, au0 - sb - bwidth);
        size_t hi = aof_tabu_face_lower_bound(T->fbox, b0, b1, au1 - sb);
        double ua[3], va[3];
        for (int k = 0; k < 3; k++) {
            int32_t vertex = ps->faces[(size_t)fa->face * 3 + (size_t)k];
            ua[k] = T->base_u[vertex] + sa;
            va[k] = T->base_v[vertex];
        }
        for (size_t ib = lo; ib < hi; ib++) {
            const AofTabuFaceBox *fb = &T->fbox[ib];
            if (fb->umax + sb <= au0 || au1 <= fb->umin + sb ||
                fa->vmax <= fb->vmin || fb->vmax <= fa->vmin)
                continue;
            double ub[3], vb[3];
            for (int k = 0; k < 3; k++) {
                int32_t vertex = ps->faces[(size_t)fb->face * 3 + (size_t)k];
                ub[k] = T->base_u[vertex] + sb;
                vb[k] = T->base_v[vertex];
            }
            T->exact_sat_tests++;
            if (!aof_triangles_overlap(ua, va, ub, vb)) continue;
            if (count != UINT32_MAX) count++;
        }
    }
    return count;
}

static uint32_t aof_tabu_exact_pair(AofTabu *T, int32_t ca, int ka,
                                    int32_t cb, int kb)
{
    if (ca == cb) return 0;
    T->exact_queries++;
    if (ca > cb) {
        int32_t tc = ca; ca = cb; cb = tc;
        int tk = ka; ka = kb; kb = tk;
    }
    if (aof_tabu_exact_aabb_separate(T, ca, ka, cb, kb)) {
        T->exact_aabb_rejects++;
        return 0;
    }
    if (T->exact.cap == 0) aof_tabu_exact_grow(T);
    uint64_t key = aof_tabu_exact_key(T, ca, ka, cb, kb);
    size_t s = aof_tabu_exact_slot(&T->exact, key);
    if (T->exact.key[s] == key) {
        T->exact_cache_hits++;
        return T->exact.count[s];
    }
    if (T->exact.used * 2 >= T->exact.cap) {
        aof_tabu_exact_grow(T);
        s = aof_tabu_exact_slot(&T->exact, key);
    }
    T->exact_cache_misses++;
    uint32_t count = aof_tabu_exact_pair_uncached(T, ca, ka, cb, kb);
    T->exact.key[s] = key;
    T->exact.count[s] = count;
    T->exact.used++;
    return count;
}

static int32_t aof_tabu_find_edge(const AofTabu *T, int32_t ca, int32_t cb)
{
    for (size_t q = T->erow[ca]; q < T->erow[(size_t)ca + 1]; q++) {
        int32_t e = T->eadj[q];
        int32_t other = T->edge[e].a == ca ? T->edge[e].b : T->edge[e].a;
        if (other == cb) return e;
    }
    return -1;
}

static int64_t aof_tabu_exact_total(AofTabu *T)
{
    int64_t total = 0;
    for (size_t a = 0; a < T->ncharts; a++) {
        for (size_t q = T->prow[a]; q < T->prow[a + 1]; q++) {
            int32_t b = T->pmember[q];
            if ((size_t)b <= a) continue;
            total += (int64_t)aof_tabu_exact_pair(
                T, (int32_t)a, T->k[a], b, T->k[b]);
        }
    }
    return total;
}

static void aof_tabu_exact_chart_pairs(AofTabu *T, int64_t *out)
{
    memset(out, 0, T->ncharts * sizeof(*out));
    for (size_t a = 0; a < T->ncharts; a++) {
        for (size_t q = T->prow[a]; q < T->prow[a + 1]; q++) {
            int32_t b = T->pmember[q];
            if ((size_t)b <= a) continue;
            int64_t n = (int64_t)aof_tabu_exact_pair(
                T, (int32_t)a, T->k[a], b, T->k[b]);
            out[a] += n;
            out[b] += n;
        }
    }
}

static void aof_tabu_exact_refresh_happy(AofTabu *T)
{
    int64_t *pairs = (int64_t *)ARENA_ALLOC(
        T->arena, T->ncharts * sizeof(int64_t));
    aof_tabu_exact_chart_pairs(T, pairs);
    for (size_t c = 0; c < T->ncharts; c++) T->happy[c] = pairs[c] == 0;
}

/* Exact collision change for an arbitrary compound placement. */
static int64_t aof_tabu_exact_delta(AofTabu *T, const int32_t *member,
                                    const int8_t *new_k, size_t nmember,
                                    int32_t *moved)
{
    int64_t delta = 0;
    for (size_t i = 0; i < nmember; i++) moved[member[i]] = (int32_t)i + 1;
    for (size_t i = 0; i < nmember; i++) {
        int32_t a = member[i];
        for (size_t q = T->prow[a]; q < T->prow[(size_t)a + 1]; q++) {
            int32_t b = T->pmember[q];
            int32_t moved_b = moved[b];
            if (moved_b == 0) {
                delta += (int64_t)aof_tabu_exact_pair(
                             T, a, new_k[i], b, T->k[b]) -
                         (int64_t)aof_tabu_exact_pair(
                             T, a, T->k[a], b, T->k[b]);
            } else if (i < (size_t)(moved_b - 1)) {
                size_t j = (size_t)(moved_b - 1);
                delta += (int64_t)aof_tabu_exact_pair(
                             T, a, new_k[i], b, new_k[j]) -
                         (int64_t)aof_tabu_exact_pair(
                             T, a, T->k[a], b, T->k[b]);
            }
        }
    }
    for (size_t i = 0; i < nmember; i++) moved[member[i]] = 0;
    return delta;
}

static int64_t aof_tabu_energy(const AofTabu *T)
{
    return T->Lov * T->e_ov + T->e_nb + T->e_rad;
}

/* Place chart c's cells at depth k.  Returns the shared-cell pairs the
 * placement creates: the sum of prior occupants over its cells. */
static int64_t aof_tabu_insert(AofTabu *T, int32_t c, int k)
{
    size_t row = (size_t)c * (size_t)T->nk + (size_t)(k + T->span);
    int64_t shift = T->cell_shift[row];
    int64_t gained = 0;
    for (size_t at = T->coff[c]; at < T->coff[(size_t)c + 1]; at++) {
        uint64_t key = (uint64_t)((int64_t)T->ckey[at] + shift);
        size_t s = aof_tabu_slot(&T->hash, key);
        if (T->hash.key[s] == 0) { T->hash.key[s] = key; T->hash.used++; }
        gained += (int64_t)T->hash.count[s];
        T->hash.count[s]++;
    }
    return gained;
}

/* Remove chart c's cells from depth k.  Returns the shared-cell pairs the
 * removal destroys.  Keys stay behind at count 0 until the next rebuild. */
static int64_t aof_tabu_remove(AofTabu *T, int32_t c, int k)
{
    size_t row = (size_t)c * (size_t)T->nk + (size_t)(k + T->span);
    int64_t shift = T->cell_shift[row];
    int64_t lost = 0;
    for (size_t at = T->coff[c]; at < T->coff[(size_t)c + 1]; at++) {
        uint64_t key = (uint64_t)((int64_t)T->ckey[at] + shift);
        size_t s = aof_tabu_slot(&T->hash, key);
        T->hash.count[s]--;
        lost += (int64_t)T->hash.count[s];
    }
    return lost;
}

/* Occupants chart c would share at depth k.  Its own cells must be absent. */
static int64_t aof_tabu_probe(const AofTabu *T, int32_t c, int k)
{
    size_t row = (size_t)c * (size_t)T->nk + (size_t)(k + T->span);
    int64_t shift = T->cell_shift[row];
    int64_t sum = 0;
    for (size_t at = T->coff[c]; at < T->coff[(size_t)c + 1]; at++) {
        uint64_t key = (uint64_t)((int64_t)T->ckey[at] + shift);
        sum += (int64_t)T->hash.count[aof_tabu_slot(&T->hash, key)];
    }
    return sum;
}

/* Reinsert everything at the current depths: reclaims stranded slots and
 * recomputes e_ov exactly.  Doubles as the audit primitive. */
static void aof_tabu_rebuild(AofTabu *T)
{
    memset(T->hash.key, 0, T->hash.cap * sizeof(*T->hash.key));
    memset(T->hash.count, 0, T->hash.cap * sizeof(*T->hash.count));
    T->hash.used = 0;
    int64_t ov = 0;
    for (size_t c = 0; c < T->ncharts; c++)
        ov += aof_tabu_insert(T, (int32_t)c, T->k[c]);
    T->e_ov = ov;
    T->rebuilds++;
}

static int64_t aof_tabu_edge_cost(const AofTabuEdge *e, int ka, int kb)
{
    int d = ka - kb - e->target;
    return e->w * (int64_t)(d < 0 ? -d : d);
}

static int64_t aof_tabu_nb_total(const AofTabu *T)
{
    int64_t total = 0;
    for (size_t e = 0; e < T->nedges; e++)
        total += aof_tabu_edge_cost(&T->edge[e], T->k[T->edge[e].a],
                                    T->k[T->edge[e].b]);
    return total;
}

static int64_t aof_tabu_rad_total(const AofTabu *T)
{
    int64_t total = 0;
    for (size_t c = 0; c < T->ncharts; c++)
        total += T->radE[c * (size_t)T->nk + (size_t)(T->k[c] + T->span)];
    return total;
}

/* Neighbour-term change if chart c moves to depth k.  steer swaps in the
 * happiness-downweighted edge weights; that variant only ranks moves -- the
 * canonical term (steer = 0) owns best-tracking, aspiration and stopping,
 * because a state-dependent energy would make iterations incomparable. */
static int64_t aof_tabu_nb_delta(const AofTabu *T, int32_t c, int k, int steer)
{
    int64_t delta = 0;
    for (size_t a = T->erow[c]; a < T->erow[(size_t)c + 1]; a++) {
        const AofTabuEdge *e = &T->edge[T->eadj[a]];
        int32_t other = e->a == c ? e->b : e->a;
        int ka = e->a == c ? k : T->k[e->a];
        int kb = e->b == c ? k : T->k[e->b];
        int dn = ka - kb - e->target;
        int dc = T->k[e->a] - T->k[e->b] - e->target;
        int64_t w = (steer && !T->happy[other]) ? e->w_unhappy : e->w;
        delta += w * (int64_t)((dn < 0 ? -dn : dn) - (dc < 0 ? -dc : dc));
    }
    return delta;
}

static void aof_tabu_refresh_happy(AofTabu *T)
{
    for (size_t c = 0; c < T->ncharts; c++) {
        aof_tabu_remove(T, (int32_t)c, T->k[c]);
        T->happy[c] = aof_tabu_probe(T, (int32_t)c, T->k[c]) == 0;
        aof_tabu_insert(T, (int32_t)c, T->k[c]);
    }
}

/* Per-chart shared-cell count at the current assignment. */
static void aof_tabu_chart_cells(AofTabu *T, int64_t *out)
{
    for (size_t c = 0; c < T->ncharts; c++) {
        aof_tabu_remove(T, (int32_t)c, T->k[c]);
        out[c] = aof_tabu_probe(T, (int32_t)c, T->k[c]);
        aof_tabu_insert(T, (int32_t)c, T->k[c]);
    }
}

/* The incremental bookkeeping re-derived from scratch must agree exactly --
 * the arithmetic is all integer, so any difference is a bug, not drift. */
static int aof_tabu_audit(AofTabu *T)
{
    int64_t ov = T->e_ov, nb = T->e_nb, rad = T->e_rad;
    aof_tabu_rebuild(T);
    T->e_nb = aof_tabu_nb_total(T);
    T->e_rad = aof_tabu_rad_total(T);
    return (ov == T->e_ov && nb == T->e_nb && rad == T->e_rad) ? 0 : -1;
}

typedef struct {
    int32_t chart;
    float   value;
} AofTabuSample;

static int aof_tabu_compare_sample(const void *pa, const void *pb)
{
    const AofTabuSample *a = (const AofTabuSample *)pa;
    const AofTabuSample *b = (const AofTabuSample *)pb;
    if (a->chart != b->chart) return a->chart < b->chart ? -1 : 1;
    if (a->value < b->value) return -1;
    if (a->value > b->value) return 1;
    return 0;
}

static int aof_tabu_compare_double(const void *pa, const void *pb)
{
    double a = *(const double *)pa, b = *(const double *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static int aof_tabu_compare_u64(const void *pa, const void *pb)
{
    uint64_t a = *(const uint64_t *)pa, b = *(const uint64_t *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

typedef struct {
    int32_t ga, gb;
    int32_t delta;              /* required base[ga] - base[gb], ga < gb */
    int64_t weight;
} AofPhaseRegionObs;

typedef struct {
    int32_t ga, gb;
    int32_t delta;
    size_t votes;
} AofPhaseRegionLink;

static int aof_tabu_compare_region_obs(const void *pa, const void *pb)
{
    const AofPhaseRegionObs *a = (const AofPhaseRegionObs *)pa;
    const AofPhaseRegionObs *b = (const AofPhaseRegionObs *)pb;
    if (a->ga != b->ga) return a->ga < b->ga ? -1 : 1;
    if (a->gb != b->gb) return a->gb < b->gb ? -1 : 1;
    if (a->delta != b->delta) return a->delta < b->delta ? -1 : 1;
    return 0;
}

/* Build a graph over welded GROUPS, not charts.  A lateral edge a--b supplies
 *     base[group(a)] - base[group(b)]
 *       = target - gbin[a] + gbin[b].
 * Requiring three observations and 70% modal agreement rejects isolated
 * coincidences while accepting the measured production relations (23/32 for
 * black--peach, and unanimous 29/29, 18/18, 45/45 for the coloured ribbon).
 * Only cycle-consistent components become regions. */
static void aof_tabu_build_phase_regions(Arena_T arena, AofTabu *T,
                                         const AtlasOverlapFixOptions *opts)
{
    size_t ngroups = T->ngroups;
    T->region_id = (int32_t *)ARENA_ALLOC(
        arena, (ngroups ? ngroups : 1) * sizeof(int32_t));
    T->region_offset = (int32_t *)ARENA_CALLOC(
        arena, ngroups ? ngroups : 1, sizeof(int32_t));
    for (size_t g = 0; g < ngroups; g++) T->region_id[g] = -1;

    if (!opts->phase_regions || !opts->phase_lateral_targets || ngroups < 2)
        return;

    size_t cap = T->nedges - T->first_lateral;
    AofPhaseRegionObs *obs = (AofPhaseRegionObs *)ARENA_ALLOC(
        arena, (cap ? cap : 1) * sizeof(AofPhaseRegionObs));
    size_t nobs = 0;
    for (size_t e = T->first_lateral; e < T->nedges; e++) {
        const AofTabuEdge *te = &T->edge[e];
        if (!te->phase_trusted) continue;
        int32_t ga = T->chart_group[te->a];
        int32_t gb = T->chart_group[te->b];
        if (ga == gb) continue;
        int32_t delta = te->target - (int32_t)T->gbin[te->a] +
                        (int32_t)T->gbin[te->b];
        if (ga > gb) {
            int32_t swap = ga; ga = gb; gb = swap;
            delta = -delta;
        }
        obs[nobs].ga = ga;
        obs[nobs].gb = gb;
        obs[nobs].delta = delta;
        obs[nobs].weight = te->w;
        nobs++;
    }
    qsort(obs, nobs, sizeof(*obs), aof_tabu_compare_region_obs);

    AofPhaseRegionLink *link = (AofPhaseRegionLink *)ARENA_ALLOC(
        arena, (nobs ? nobs : 1) * sizeof(AofPhaseRegionLink));
    size_t nlink = 0;
    for (size_t p = 0; p < nobs; ) {
        size_t pair_end = p + 1;
        while (pair_end < nobs && obs[pair_end].ga == obs[p].ga &&
               obs[pair_end].gb == obs[p].gb)
            pair_end++;
        size_t total_votes = pair_end - p;
        size_t best_votes = 0;
        int64_t best_weight = -1;
        int32_t best_delta = 0;
        for (size_t q = p; q < pair_end; ) {
            size_t run_end = q + 1;
            int64_t run_weight = obs[q].weight;
            while (run_end < pair_end &&
                   obs[run_end].delta == obs[q].delta) {
                run_weight += obs[run_end].weight;
                run_end++;
            }
            size_t run_votes = run_end - q;
            if (run_weight > best_weight ||
                (run_weight == best_weight && run_votes > best_votes)) {
                best_weight = run_weight;
                best_votes = run_votes;
                best_delta = obs[q].delta;
            }
            q = run_end;
        }
        if (total_votes >= 3 && best_votes * 10 >= total_votes * 7) {
            link[nlink].ga = obs[p].ga;
            link[nlink].gb = obs[p].gb;
            link[nlink].delta = best_delta;
            link[nlink].votes = best_votes;
            nlink++;
        }
        p = pair_end;
    }

    uint8_t *seen = (uint8_t *)ARENA_CALLOC(
        arena, ngroups, sizeof(uint8_t));
    int32_t *potential = (int32_t *)ARENA_CALLOC(
        arena, ngroups, sizeof(int32_t));
    int32_t *queue = (int32_t *)ARENA_ALLOC(
        arena, ngroups * sizeof(int32_t));
    for (size_t root = 0; root < ngroups; root++) {
        int incident = 0;
        for (size_t l = 0; l < nlink && !incident; l++)
            incident = link[l].ga == (int32_t)root ||
                       link[l].gb == (int32_t)root;
        if (seen[root] || !incident) continue;
        size_t head = 0, tail = 0;
        int bad = 0;
        seen[root] = 1;
        potential[root] = 0;
        queue[tail++] = (int32_t)root;
        while (head < tail) {
            int32_t g = queue[head++];
            for (size_t l = 0; l < nlink; l++) {
                int32_t other, candidate;
                if (link[l].ga == g) {
                    other = link[l].gb;
                    candidate = potential[g] - link[l].delta;
                } else if (link[l].gb == g) {
                    other = link[l].ga;
                    candidate = potential[g] + link[l].delta;
                } else continue;
                if (!seen[other]) {
                    seen[other] = 1;
                    potential[other] = candidate;
                    queue[tail++] = other;
                } else if (potential[other] != candidate) {
                    bad = 1;
                }
            }
        }
        if (bad || tail < 2) continue;
        double mean = 0.0, mass = 0.0;
        for (size_t q = 0; q < tail; q++) {
            int32_t g = queue[q];
            double w = (double)(T->grow[g + 1] - T->grow[g]);
            mean += w * (double)potential[g];
            mass += w;
        }
        int32_t centre = (int32_t)floor(mean / mass + 0.5);
        int32_t rid = (int32_t)T->n_region++;
        for (size_t q = 0; q < tail; q++) {
            int32_t g = queue[q];
            T->region_id[g] = rid;
            T->region_offset[g] = potential[g] - centre;
        }
    }

    for (size_t l = 0; l < nlink; l++)
        if (T->region_id[link[l].ga] >= 0 &&
            T->region_id[link[l].ga] == T->region_id[link[l].gb])
            T->n_region_links++;

    T->region_row = (size_t *)ARENA_CALLOC(
        arena, T->n_region + 1, sizeof(size_t));
    for (size_t g = 0; g < ngroups; g++)
        if (T->region_id[g] >= 0) T->region_row[T->region_id[g] + 1]++;
    for (size_t r = 0; r < T->n_region; r++)
        T->region_row[r + 1] += T->region_row[r];
    T->region_group = (int32_t *)ARENA_ALLOC(
        arena, (T->region_row[T->n_region] ?
                T->region_row[T->n_region] : 1) * sizeof(int32_t));
    size_t *fill = (size_t *)ARENA_ALLOC(
        arena, (T->n_region ? T->n_region : 1) * sizeof(size_t));
    for (size_t r = 0; r < T->n_region; r++) fill[r] = T->region_row[r];
    for (size_t g = 0; g < ngroups; g++)
        if (T->region_id[g] >= 0)
            T->region_group[fill[T->region_id[g]]++] = (int32_t)g;
}

/* Median per chart of one per-vertex quantity: one sort, runs by chart.
 * Charts with no finite sample keep their previous value. */
static void aof_tabu_chart_median(AofTabuSample *sample, size_t n,
                                  double *out_median)
{
    qsort(sample, n, sizeof(*sample), aof_tabu_compare_sample);
    for (size_t first = 0; first < n; ) {
        size_t last = first + 1;
        while (last < n && sample[last].chart == sample[first].chart) last++;
        out_median[sample[first].chart] =
            (double)sample[first + (last - first) / 2].value;
        first = last;
    }
}

static int aof_tabu_prepare(Arena_T arena, const PieceSet *ps,
                            const ScaffoldCalib *cal,
                            const double *u, const double *v,
                            const uint8_t *face_keep,
                            const int32_t *vertex_chart, size_t ncharts,
                            const int32_t *chart_group, size_t ngroups,
                            const AofNeighbour *edge, size_t nedges,
                            const AtlasOverlapFixGroup *ginfo,
                            int wind_corrected,
                            const AtlasOverlapFixOptions *opts,
                            size_t hash_cap_override,
                            AofTabu *T, size_t *out_edges_excluded)
{
    memset(T, 0, sizeof(*T));
    *out_edges_excluded = 0;
    if (ncharts == 0 || opts->tabu_span < 1 || opts->tabu_span > 100 ||
        !(opts->tabu_cell > 0.0) ||
        (opts->tabu_active_core &&
         (!(opts->tabu_core_z_lo < opts->tabu_core_z_hi) ||
          !(opts->tabu_core_y_lo < opts->tabu_core_y_hi) ||
          !(opts->tabu_core_x_lo < opts->tabu_core_x_hi))) ||
        (opts->tabu_ladder_override &&
         (!isfinite(opts->tabu_ladder_step) ||
          fabs(opts->tabu_ladder_step) <= 1.0e-9 ||
          !isfinite(opts->tabu_ladder_u0))))
        return -1;
    T->ncharts = ncharts;
    T->ngroups = ngroups;
    T->span = opts->tabu_span;
    T->nk = 2 * opts->tabu_span + 1;
    T->arena = arena;
    T->ps = ps;
    T->base_u = u;
    T->base_v = v;
    T->vertex_chart = vertex_chart;
    T->Lov = (int64_t)floor(opts->lambda_ov * AOF_TABU_FIX + 0.5);
    T->cell = opts->tabu_cell;
    T->chart_group = chart_group;
    size_t nk = (size_t)T->nk;

    /* Group membership CSR (needed below for the radial recentring). */
    T->grow = (size_t *)ARENA_CALLOC(arena, ngroups + 1, sizeof(size_t));
    for (size_t c = 0; c < ncharts; c++) T->grow[chart_group[c] + 1]++;
    for (size_t g = 0; g < ngroups; g++) T->grow[g + 1] += T->grow[g];
    T->gmember = (int32_t *)ARENA_ALLOC(arena, ncharts * sizeof(int32_t));
    {
        size_t *fill = (size_t *)ARENA_ALLOC(arena, ngroups * sizeof(size_t));
        for (size_t g = 0; g < ngroups; g++) fill[g] = T->grow[g];
        for (size_t c = 0; c < ncharts; c++)
            T->gmember[fill[chart_group[c]]++] = (int32_t)c;
    }

    /* --- per-chart medians of u, radius, and raw phi --------------------- */
    T->u_med = (double *)ARENA_CALLOC(arena, ncharts, sizeof(double));
    T->r_med = (double *)ARENA_CALLOC(arena, ncharts, sizeof(double));
    T->w_phi0 = (double *)ARENA_CALLOC(arena, ncharts, sizeof(double));
    T->w_rad = (double *)ARENA_CALLOC(arena, ncharts, sizeof(double));
    double *phi_raw = (double *)ARENA_CALLOC(arena, ncharts, sizeof(double));
    uint8_t *phi_valid = (uint8_t *)ARENA_CALLOC(
        arena, ncharts, sizeof(uint8_t));
    {
        AofTabuSample *sample = (AofTabuSample *)ARENA_ALLOC(
            arena, ps->nv * sizeof(*sample));
        size_t n = 0;
        for (size_t i = 0; i < ps->nv; i++) {
            if (vertex_chart[i] < 0 || !isfinite(u[i])) continue;
            sample[n].chart = vertex_chart[i];
            sample[n].value = (float)u[i];
            n++;
        }
        aof_tabu_chart_median(sample, n, T->u_med);
        n = 0;
        for (size_t i = 0; i < ps->nv; i++) {
            if (vertex_chart[i] < 0) continue;
            sample[n].chart = vertex_chart[i];
            sample[n].value = (float)aof_radius(cal, &ps->verts[i * 3]);
            n++;
        }
        aof_tabu_chart_median(sample, n, T->r_med);
        n = 0;
        for (size_t i = 0; i < ps->nv; i++) {
            if (vertex_chart[i] < 0 || ps->phi == NULL ||
                !isfinite((double)ps->phi[i]))
                continue;
            sample[n].chart = vertex_chart[i];
            sample[n].value = ps->phi[i];
            phi_valid[vertex_chart[i]] = 1;
            n++;
        }
        aof_tabu_chart_median(sample, n, phi_raw);
    }

    /* --- chart -> cube provenance (the ladder pairs by cube) -------------- */
    T->chart_cube = (int32_t *)ARENA_ALLOC(arena, ncharts * sizeof(int32_t));
    for (size_t c = 0; c < ncharts; c++) T->chart_cube[c] = -1;
    if (ps->face_cube != NULL) {
        for (size_t f = 0; f < ps->nf; f++) {
            if (!face_keep[f]) continue;
            int32_t c = vertex_chart[ps->faces[f * 3]];
            if (c >= 0 && T->chart_cube[c] < 0)
                T->chart_cube[c] = ps->face_cube[f];
        }
    }
    T->active = (uint8_t *)ARENA_ALLOC(arena, ncharts * sizeof(uint8_t));
    for (size_t c = 0; c < ncharts; c++) {
        int active = 1;
        if (opts->tabu_active_core) {
            int32_t q = T->chart_cube[c];
            active = q >= 0 && (size_t)q < ps->n_cubes &&
                ps->cube_org[q][0] >= opts->tabu_core_z_lo &&
                ps->cube_org[q][0] <  opts->tabu_core_z_hi &&
                ps->cube_org[q][1] >= opts->tabu_core_y_lo &&
                ps->cube_org[q][1] <  opts->tabu_core_y_hi &&
                ps->cube_org[q][2] >= opts->tabu_core_x_lo &&
                ps->cube_org[q][2] <  opts->tabu_core_x_hi;
        }
        T->active[c] = (uint8_t)active;
    }

    /* --- winds ------------------------------------------------------------ */
    double a = cal->spiral_a, b = cal->spiral_b;
    double pi = AOF_2PI * 0.5;
    double sense = cal->sense < 0 ? -1.0 : 1.0;
    T->dr_dw = b * sense;                  /* radius gained per wind outward */
    T->w_geo = (double *)ARENA_CALLOC(arena, ncharts, sizeof(double));
    for (size_t c = 0; c < ncharts; c++)
        T->w_geo[c] = fabs(T->dr_dw) > 1.0e-9
                    ? (T->r_med[c] - a) / T->dr_dw : 0.0;

    /* Measure the FIELD's wind ladder.  Same-cube chart pairs about one
     * pitch apart in radius share an angular sector, so axis wander cancels
     * and du/(2*pi*r_mid) is the field's own per-wind u step at that spot.
     * The signed median over usable pairs gives the direction (this field's
     * u GROWS outward -- the opposite of the spiral map, which is exactly
     * the reprojection bug this measurement retires); the scale then comes
     * from the sign-agreeing majority only, because mis-wound charts -- the
     * disease this solver treats -- produce sign-flipped pairs. */
    /* cube -> charts CSR: the ladder and the lateral coupling both pair
     * charts through cube identity/adjacency. */
    size_t *crow = NULL;
    int32_t *cmember = NULL;
    if (ps->n_cubes > 0) {
        crow = (size_t *)ARENA_CALLOC(arena, ps->n_cubes + 1, sizeof(size_t));
        for (size_t c = 0; c < ncharts; c++)
            if (T->chart_cube[c] >= 0) crow[T->chart_cube[c] + 1]++;
        for (size_t q = 0; q < ps->n_cubes; q++) crow[q + 1] += crow[q];
        cmember = (int32_t *)ARENA_ALLOC(
            arena, (ncharts ? ncharts : 1) * sizeof(int32_t));
        size_t *fill = (size_t *)ARENA_ALLOC(
            arena, ps->n_cubes * sizeof(size_t));
        for (size_t q = 0; q < ps->n_cubes; q++) fill[q] = crow[q];
        for (size_t c = 0; c < ncharts; c++)
            if (T->chart_cube[c] >= 0)
                cmember[fill[T->chart_cube[c]]++] = (int32_t)c;
    }
    T->crow = crow;
    T->cmember = cmember;

    T->ladder_measured = 0;
    T->ladder_fixed = 0;
    T->ladder_step = 0.0;
    T->ladder_u0 = 0.0;
    T->ladder_pairs = 0;
    if (opts->tabu_ladder_override) {
        T->ladder_measured = 1;
        T->ladder_fixed = 1;
        T->ladder_step = opts->tabu_ladder_step;
        T->ladder_u0 = opts->tabu_ladder_u0;
    } else if (crow != NULL && fabs(T->dr_dw) > 1.0e-9) {
        double lo = 0.7 * fabs(T->dr_dw), hi = 1.4 * fabs(T->dr_dw);
        size_t ratio_cap = ncharts * 8;
        double *ratio = (double *)ARENA_ALLOC(
            arena, (ratio_cap ? ratio_cap : 1) * sizeof(double));
        size_t nratio = 0;
        for (size_t q = 0; q < ps->n_cubes; q++) {
            for (size_t i = crow[q]; i < crow[q + 1]; i++)
            for (size_t j = i + 1; j < crow[q + 1]; j++) {
                int32_t ci = cmember[i], cj = cmember[j];
                double dr = T->r_med[cj] - T->r_med[ci];
                double adr = fabs(dr);
                if (adr < lo || adr > hi) continue;
                double r_mid = 0.5 * (T->r_med[ci] + T->r_med[cj]);
                if (!(r_mid > 1.0e-9)) continue;
                double du = dr > 0.0 ? T->u_med[cj] - T->u_med[ci]
                                     : T->u_med[ci] - T->u_med[cj];
                double rat = du / (AOF_2PI * r_mid);
                if (fabs(rat) > 0.3 && fabs(rat) < 1.6 &&
                    nratio < ratio_cap)
                    ratio[nratio++] = rat;
            }
        }
        if (nratio >= 8) {
            qsort(ratio, nratio, sizeof(double), aof_tabu_compare_double);
            double sig = ratio[nratio / 2] < 0.0 ? -1.0 : 1.0;
            size_t nsame = 0;
            for (size_t i = 0; i < nratio; i++)
                if (ratio[i] * sig > 0.0) ratio[nsame++] = ratio[i];
            /* the filter preserves the sort, so the median is direct */
            T->ladder_step = ratio[nsame / 2];
            T->ladder_measured = 1;
            T->ladder_pairs = nratio;
        }
    }

    /* Wind implied by u.  Measured ladder: invert
     * u = u0 + step * arc(w),  arc(w) = 2*pi*(a*w + dr_dw*w^2/2),
     * with u0 the median offset over charts at their radius-implied wind.
     * Fallback (too few pairs -- the synthetic fixtures): invert the spiral
     * map itself at the chart's median u, root nearest the raw median phi
     * (raw phi's whole turns are per-cube garbage; it only picks branches). */
    if (T->ladder_measured) {
        double u0 = T->ladder_u0;
        if (!T->ladder_fixed) {
            double *off_u = (double *)ARENA_ALLOC(
                arena, ncharts * sizeof(double));
            for (size_t c = 0; c < ncharts; c++) {
                double w = T->w_geo[c];
                double arc = AOF_2PI * (a * w + 0.5 * T->dr_dw * w * w);
                off_u[c] = T->u_med[c] - T->ladder_step * arc;
            }
            qsort(off_u, ncharts, sizeof(double), aof_tabu_compare_double);
            u0 = off_u[ncharts / 2];
        }
        T->ladder_u0 = u0;
        for (size_t c = 0; c < ncharts; c++) {
            double w = T->w_geo[c];
            double rhs = (T->u_med[c] - u0) / (AOF_2PI * T->ladder_step);
            double half_p = 0.5 * T->dr_dw;   /* half_p*w^2 + a*w = rhs */
            if (fabs(half_p) > 1.0e-12) {
                double disc = a * a + 4.0 * half_p * rhs;
                if (disc >= 0.0) {
                    double root = sqrt(disc);
                    double p1 = (-a + root) / (2.0 * half_p);
                    double p2 = (-a - root) / (2.0 * half_p);
                    w = fabs(p1 - T->w_geo[c]) <= fabs(p2 - T->w_geo[c])
                      ? p1 : p2;
                }
            } else if (fabs(a) > 1.0e-12) {
                w = rhs / a;
            }
            T->w_phi0[c] = w;
        }
    } else {
        for (size_t c = 0; c < ncharts; c++) {
            double phi = phi_raw[c];
            if (fabs(b) > 1.0e-12) {
                double disc = a * a + b * T->u_med[c] / pi;
                if (disc >= 0.0) {
                    double root = sqrt(disc);
                    double p1 = (-a + root) * AOF_2PI / b;
                    double p2 = (-a - root) * AOF_2PI / b;
                    phi = fabs(p1 - phi_raw[c]) <= fabs(p2 - phi_raw[c])
                        ? p1 : p2;
                }
            } else if (fabs(a) > 1.0e-12) {
                phi = T->u_med[c] / a;
            }
            T->w_phi0[c] = sense * phi / AOF_2PI;
        }
    }

    /* Radius-implied wind, recentred.  The axis wanders ~85 vox, so absolute
     * radius wind is only ~57% reliable while nearby RELATIVE radius is fine:
     * subtract the median disagreement (w_rad - w_phi0) per welded group
     * (>= 4 charts; smaller groups have no robust local median and use the
     * global one), so the prior pulls on relative order, not on the wander. */
    {
        double pitch = cal->pitch;
        double *off = (double *)ARENA_ALLOC(arena, ncharts * sizeof(double));
        double *scratch = (double *)ARENA_ALLOC(arena,
                                                ncharts * sizeof(double));
        for (size_t c = 0; c < ncharts; c++) {
            double w = pitch > 0.0 ? (T->r_med[c] - a) / pitch : T->w_phi0[c];
            off[c] = w - T->w_phi0[c];
            T->w_rad[c] = w;
        }
        memcpy(scratch, off, ncharts * sizeof(double));
        qsort(scratch, ncharts, sizeof(double), aof_tabu_compare_double);
        double global_off = scratch[ncharts / 2];
        for (size_t g = 0; g < ngroups; g++) {
            size_t m0 = T->grow[g], m1 = T->grow[g + 1];
            double use = global_off;
            if (m1 - m0 >= 4) {
                for (size_t m = m0; m < m1; m++)
                    scratch[m - m0] = off[T->gmember[m]];
                qsort(scratch, m1 - m0, sizeof(double),
                      aof_tabu_compare_double);
                use = scratch[(m1 - m0) / 2];
            }
            for (size_t m = m0; m < m1; m++)
                T->w_rad[T->gmember[m]] -= use;
        }
    }

    /*
     * Source-phase lift.  phase_target is not a desired tabu depth difference:
     * it is the integer gauge needed to make independently branched per-chart
     * source angles agree at a fixed XYZ contact.  Integrate that gauge first,
     * then compare the lifted source winding with the winding already encoded
     * by the incoming u field.  Their median-centred difference is the number
     * of turns actually MISSING from the initial atlas.
     */
    T->w_phase = (double *)ARENA_ALLOC(arena, ncharts * sizeof(double));
    T->phase_gauge = (int32_t *)ARENA_CALLOC(
        arena, ncharts, sizeof(int32_t));
    T->phase_island = (int32_t *)ARENA_ALLOC(
        arena, ncharts * sizeof(int32_t));
    for (size_t c = 0; c < ncharts; c++) T->phase_island[c] = -1;
    T->phase_bin = (int8_t *)ARENA_CALLOC(arena, ncharts, sizeof(int8_t));
    T->phase_trusted_group = (uint8_t *)ARENA_CALLOC(
        arena, ngroups ? ngroups : 1, sizeof(uint8_t));
    T->phase_flat_group = (uint8_t *)ARENA_CALLOC(
        arena, ngroups ? ngroups : 1, sizeof(uint8_t));
    {
        uint8_t *seen = (uint8_t *)ARENA_CALLOC(
            arena, ncharts, sizeof(uint8_t));
        int32_t *queue = (int32_t *)ARENA_ALLOC(
            arena, ncharts * sizeof(int32_t));
        double *scratch = (double *)ARENA_ALLOC(
            arena, ncharts * sizeof(double));

        for (size_t g = 0; g < ngroups; g++) {
            size_t m0 = T->grow[g], m1 = T->grow[g + 1];
            size_t island_before = T->phase_gauge_islands;
            for (size_t m = m0; m < m1; m++) {
                int32_t root = T->gmember[m];
                if (!phi_valid[root] || seen[root]) continue;
                size_t head = 0, tail = 0;
                seen[root] = 1;
                T->phase_gauge[root] = 0;
                T->phase_island[root] = (int32_t)T->phase_gauge_islands;
                queue[tail++] = root;
                T->phase_gauge_islands++;
                while (head < tail) {
                    int32_t c = queue[head++];
                    for (size_t e = 0; e < nedges; e++) {
                        int32_t other;
                        int32_t candidate;
                        if (edge[e].phase_residual >
                                opts->wind_residual_limit ||
                            (opts->wind_cap > 0 &&
                             abs(edge[e].phase_target) > opts->wind_cap))
                            continue;
                        if (edge[e].chart0 == c) {
                            other = edge[e].chart1;
                            candidate = T->phase_gauge[c] -
                                        edge[e].phase_target;
                        } else if (edge[e].chart1 == c) {
                            other = edge[e].chart0;
                            candidate = T->phase_gauge[c] +
                                        edge[e].phase_target;
                        } else {
                            continue;
                        }
                        if (T->chart_group[other] != (int32_t)g ||
                            !phi_valid[other] || seen[other])
                            continue;
                        seen[other] = 1;
                        T->phase_gauge[other] = candidate;
                        T->phase_island[other] = T->phase_island[c];
                        queue[tail++] = other;
                    }
                }
            }

            size_t nvalid = 0;
            for (size_t m = m0; m < m1; m++) {
                int32_t c = T->gmember[m];
                if (!phi_valid[c]) {
                    T->w_phase[c] = NAN;
                    continue;
                }
                T->w_phase[c] = sense * phi_raw[c] / AOF_2PI +
                                (double)T->phase_gauge[c];
                scratch[nvalid++] = T->w_phase[c] - T->w_phi0[c];
            }
            if (nvalid > 0) {
                qsort(scratch, nvalid, sizeof(double),
                      aof_tabu_compare_double);
                double centre = scratch[nvalid / 2];
                for (size_t m = m0; m < m1; m++) {
                    int32_t c = T->gmember[m];
                    if (!phi_valid[c]) continue;
                    double d = T->w_phase[c] - T->w_phi0[c] - centre;
                    int bin = (int)floor(d + 0.5);
                    if (bin < -127) bin = -127;
                    if (bin > 127) bin = 127;
                    T->phase_bin[c] = (int8_t)bin;
                }
            }
            /* A hard coherence relation is trusted only when every member has
             * phase, the usable phase edges connect the whole group, and all
             * members agree that the incoming atlas is missing the same number
             * of winds.  Cycle consistency is checked immediately below. */
            if (m1 - m0 >= 2 && nvalid == m1 - m0 &&
                T->phase_gauge_islands == island_before + 1) {
                T->phase_trusted_group[g] = 1;
                T->phase_flat_group[g] = 1;
                for (size_t m = m0; m < m1; m++)
                    if (T->phase_bin[T->gmember[m]] != 0)
                        T->phase_flat_group[g] = 0;
            }
        }
        for (size_t e = 0; e < nedges; e++) {
            if (edge[e].phase_residual > opts->wind_residual_limit ||
                (opts->wind_cap > 0 &&
                 abs(edge[e].phase_target) > opts->wind_cap) ||
                !phi_valid[edge[e].chart0] || !phi_valid[edge[e].chart1])
                continue;
            if (T->phase_gauge[edge[e].chart0] -
                    T->phase_gauge[edge[e].chart1] != edge[e].phase_target) {
                T->phase_gauge_inconsistent++;
                T->phase_trusted_group[
                    T->chart_group[edge[e].chart0]] = 0;
                T->phase_flat_group[
                    T->chart_group[edge[e].chart0]] = 0;
            }
        }
    }

    /* --- base occupancy cells --------------------------------------------- */
    size_t *frow = (size_t *)ARENA_CALLOC(arena, ncharts + 1, sizeof(size_t));
    for (size_t f = 0; f < ps->nf; f++) {
        if (!face_keep[f]) continue;
        int32_t c = vertex_chart[ps->faces[f * 3]];
        if (c >= 0) frow[c + 1]++;
    }
    for (size_t c = 0; c < ncharts; c++) frow[c + 1] += frow[c];
    size_t nkept = frow[ncharts];
    int32_t *fidx = (int32_t *)ARENA_ALLOC(
        arena, (nkept ? nkept : 1) * sizeof(int32_t));
    {
        size_t *fill = (size_t *)ARENA_ALLOC(arena, ncharts * sizeof(size_t));
        for (size_t c = 0; c < ncharts; c++) fill[c] = frow[c];
        for (size_t f = 0; f < ps->nf; f++) {
            if (!face_keep[f]) continue;
            int32_t c = vertex_chart[ps->faces[f * 3]];
            if (c < 0) continue;
            fidx[fill[c]++] = (int32_t)f;
        }
    }

    /* Exact broad phase: face boxes sorted by u inside each chart. */
    T->frow = frow;
    T->fbox = (AofTabuFaceBox *)ARENA_ALLOC(
        arena, (nkept ? nkept : 1) * sizeof(AofTabuFaceBox));
    T->chart_box = (double *)ARENA_ALLOC(
        arena, (ncharts ? ncharts * 4 : 1) * sizeof(double));
    T->chart_face_umax = (double *)ARENA_CALLOC(
        arena, (ncharts ? ncharts : 1), sizeof(double));
    for (size_t c = 0; c < ncharts; c++) {
        double *cb = &T->chart_box[c * 4];
        cb[0] = DBL_MAX; cb[1] = -DBL_MAX;
        cb[2] = DBL_MAX; cb[3] = -DBL_MAX;
        for (size_t q = frow[c]; q < frow[c + 1]; q++) {
            AofTabuFaceBox *fb = &T->fbox[q];
            fb->face = fidx[q];
            double box[4];
            aof_face_box(ps, u, v, (size_t)fb->face, box);
            fb->umin = box[0]; fb->umax = box[1];
            fb->vmin = box[2]; fb->vmax = box[3];
            if (box[0] < cb[0]) cb[0] = box[0];
            if (box[1] > cb[1]) cb[1] = box[1];
            if (box[2] < cb[2]) cb[2] = box[2];
            if (box[3] > cb[3]) cb[3] = box[3];
            double width = box[1] - box[0];
            if (width > T->chart_face_umax[c]) T->chart_face_umax[c] = width;
        }
        qsort(&T->fbox[frow[c]], frow[c + 1] - frow[c],
              sizeof(AofTabuFaceBox), aof_tabu_compare_face_box);
    }

    T->ckey = (uint64_t *)ARENA_ALLOC(
        arena, (nkept ? nkept * 4 : 1) * sizeof(uint64_t));
    T->coff = (size_t *)ARENA_ALLOC(arena, (ncharts + 1) * sizeof(size_t));
    size_t at = 0;
    for (size_t c = 0; c < ncharts; c++) {
        T->coff[c] = at;
        size_t start = at;
        for (size_t q = frow[c]; q < frow[c + 1]; q++) {
            size_t f = (size_t)fidx[q];
            double cu = 0.0, cv = 0.0;
            for (int corner = 0; corner < 3; corner++) {
                int32_t vertex = ps->faces[f * 3 + (size_t)corner];
                T->ckey[at++] = aof_tabu_cell_key(u[vertex], v[vertex],
                                                  T->cell);
                cu += u[vertex];
                cv += v[vertex];
            }
            /* Vertices alone under-cover faces wider than a cell. */
            T->ckey[at++] = aof_tabu_cell_key(cu / 3.0, cv / 3.0, T->cell);
        }
        qsort(&T->ckey[start], at - start, sizeof(uint64_t),
              aof_tabu_compare_u64);
        size_t w = start;
        for (size_t r = start; r < at; r++)
            if (r == start || T->ckey[r] != T->ckey[w - 1])
                T->ckey[w++] = T->ckey[r];
        at = w;
    }
    T->coff[ncharts] = at;

    /* --- usable weld edges + lateral same-wind coupling ------------------- */
    int64_t Lnb = (int64_t)floor(opts->lambda_nb * AOF_TABU_FIX + 0.5);

    /* Happy families: groups already assembled correctly, whose only problem
     * is WHERE they sit.  "Correct" means no WRONG-WIND internal overlap --
     * fold pairs (near-zero radial separation) are doubled geometry that is
     * parameterized correctly, and counting them as guilt shredded exactly
     * the groups the eye says are fine.  Family bonds are boosted and exempt
     * from the unhappy-steering discount, so members cannot cheaply defect
     * one at a time -- the whole-group move, which pays no internal cost,
     * becomes the cheap escape. */
    uint8_t *fam = (uint8_t *)ARENA_CALLOC(arena, ngroups, sizeof(uint8_t));
    if (ginfo != NULL && opts->family_w > 1.0) {
        for (size_t g = 0; g < ngroups; g++)
            fam[g] = ginfo[g].intra_wrongwind == 0 && ginfo[g].ncharts >= 2;
        /* A member displaced by whole winds RELATIVE TO ITS OWN GROUP shows
         * up as a nonzero-wraps internal weld (it need not overlap anything
         * to be wrong).  Such a group is not happy: rigid treatment would
         * freeze the displacement in. */
        if (!wind_corrected)
            for (size_t e = 0; e < nedges; e++)
                if (edge[e].wraps != 0)
                    fam[chart_group[edge[e].chart0]] = 0;
    }

    /* Lateral pairs: same radius (same wind, whatever the shared gauge
     * error), adjacent cubes, side by side in u, NOT already welded -- one
     * sheet interrupted by a hole or a trimmed seam.  The weld graph cannot
     * see these, and without them the search tears jigsaw neighbours apart
     * for free (measured: 60% of such pairs torn, every one via unequal k). */
    uint64_t *lat_pair = NULL;
    size_t nlat = 0;
    if (opts->lateral_w > 0.0 && crow != NULL && ps->cube_org != NULL &&
        fabs(T->dr_dw) > 1.0e-9) {
        uint64_t *weld_key = (uint64_t *)ARENA_ALLOC(
            arena, (nedges ? nedges : 1) * sizeof(uint64_t));
        for (size_t e = 0; e < nedges; e++)
            weld_key[e] = ((uint64_t)(uint32_t)edge[e].chart0 << 32) |
                          (uint64_t)(uint32_t)edge[e].chart1;
        qsort(weld_key, nedges, sizeof(uint64_t), aof_tabu_compare_u64);
        size_t lat_cap = ncharts * 8;
        lat_pair = (uint64_t *)ARENA_ALLOC(
            arena, (lat_cap ? lat_cap : 1) * sizeof(uint64_t));
        double r_gate = 0.6 * fabs(T->dr_dw);
        for (size_t q1 = 0; q1 < ps->n_cubes; q1++)
        for (size_t q2 = q1; q2 < ps->n_cubes; q2++) {
            if (labs(ps->cube_org[q1][0] - ps->cube_org[q2][0]) > 128 ||
                labs(ps->cube_org[q1][1] - ps->cube_org[q2][1]) > 128 ||
                labs(ps->cube_org[q1][2] - ps->cube_org[q2][2]) > 128)
                continue;
            for (size_t i = crow[q1]; i < crow[q1 + 1]; i++) {
                size_t j0 = q1 == q2 ? i + 1 : crow[q2];
                for (size_t j = j0; j < crow[q2 + 1]; j++) {
                    int32_t ca = cmember[i], cb = cmember[j];
                    if (fabs(T->r_med[ca] - T->r_med[cb]) > r_gate) continue;
                    if (fabs(T->u_med[ca] - T->u_med[cb]) >
                        opts->lateral_u_gap) continue;
                    if (ca > cb) { int32_t t = ca; ca = cb; cb = t; }
                    uint64_t key = ((uint64_t)(uint32_t)ca << 32) |
                                   (uint64_t)(uint32_t)cb;
                    if (nedges > 0 &&
                        bsearch(&key, weld_key, nedges, sizeof(uint64_t),
                                aof_tabu_compare_u64) != NULL)
                        continue;   /* welded already: stronger evidence */
                    if (nlat < lat_cap) lat_pair[nlat++] = key;
                }
            }
        }
    }

    T->edge = (AofTabuEdge *)ARENA_ALLOC(
        arena, (nedges + nlat ? nedges + nlat : 1) * sizeof(AofTabuEdge));
    size_t ne = 0;
    for (size_t e = 0; e < nedges; e++) {
        /* After the wind correction the measured mismatch at every trusted
         * weld is zero turns and the target is 0.  Without it, the edge's
         * own du says what k[a] - k[b] must be for the weld to close -- in
         * units of the FIELD's per-wind step at the contact radius,
         * whichever direction and scale that runs. */
        int target = 0;
        int usable = 1;
        if (opts->phase_targets) {
            target = edge[e].phase_target;
            if (edge[e].phase_residual > opts->wind_residual_limit ||
                (opts->wind_cap > 0 && abs(target) > opts->wind_cap))
                usable = 0;
        } else if (wind_corrected) {
            if (edge[e].residual > opts->wind_residual_limit ||
                (opts->wind_cap > 0 && abs(edge[e].wraps) > opts->wind_cap))
                usable = 0;
        } else {
            double phi_e = fabs(b) > 1.0e-12
                         ? (edge[e].radius - a) * AOF_2PI / b : 0.0;
            double q = T->ladder_measured
                     ? T->ladder_step * AOF_2PI * edge[e].radius
                     : CubeReg_deltaU(a, b, phi_e, cal->sense);
            double t = fabs(q) > 1.0e-9 ? -edge[e].du_median / q : 0.0;
            target = (int)floor(t + 0.5);
            if (fabs(t - (double)target) > opts->wind_residual_limit ||
                (opts->wind_cap > 0 && abs(target) > opts->wind_cap))
                usable = 0;
        }
        if (!usable) {
            (*out_edges_excluded)++;
            continue;
        }
        size_t contacts = edge[e].contacts;
        if (contacts > AOF_TABU_CONTACT_CAP) contacts = AOF_TABU_CONTACT_CAP;
        T->edge[ne].a = edge[e].chart0;
        T->edge[ne].b = edge[e].chart1;
        T->edge[ne].target = target;
        T->edge[ne].w = Lnb * (int64_t)contacts;
        T->edge[ne].phase_trusted = 0;
        T->edge[ne].du = edge[e].du_median;
        /* Welds only exist inside one group, so the family test is direct. */
        T->edge[ne].family = fam[chart_group[edge[e].chart0]];
        if (T->edge[ne].family) {
            T->edge[ne].w = (int64_t)floor(
                (double)T->edge[ne].w * opts->family_w + 0.5);
            T->edge[ne].w_unhappy = T->edge[ne].w;
        } else {
            T->edge[ne].w_unhappy = (int64_t)floor(
                (double)T->edge[ne].w * opts->nb_unhappy_w + 0.5);
        }
        ne++;
    }
    T->first_lateral = ne;
    if (nlat > 0) {
        int64_t w_lat = (int64_t)floor(
            opts->lambda_nb * opts->lateral_w * AOF_TABU_FIX + 0.5);
        int64_t w_lat_un = (int64_t)floor(
            (double)w_lat * opts->nb_unhappy_w + 0.5);
        for (size_t l = 0; l < nlat; l++) {
            int32_t ca = (int32_t)(lat_pair[l] >> 32);
            int32_t cb = (int32_t)(uint32_t)lat_pair[l];
            T->edge[ne].a = ca;
            T->edge[ne].b = cb;
            int target = 0;
            int phase_trusted = 0;
            double da = T->w_phase[ca] - T->w_phi0[ca];
            double db = T->w_phase[cb] - T->w_phi0[cb];
            if (isfinite(da) && isfinite(db)) {
                double q = da - db;
                int qi = (int)floor(q + 0.5);
                if (fabs(q - (double)qi) <= opts->wind_residual_limit &&
                    (opts->wind_cap <= 0 || abs(qi) <= opts->wind_cap)) {
                    phase_trusted = 1;
                    T->phase_lateral_trusted++;
                    if (qi != 0) T->phase_lateral_nonzero++;
                    if (opts->phase_lateral_targets) target = qi;
                }
            }
            /* Equal radius does not imply equal atlas gauge.  When requested,
             * the source lift supplies the integer relation needed to keep a
             * continuous lateral region aligned across welded-group gauges. */
            T->edge[ne].target = target;
            T->edge[ne].phase_trusted = (uint8_t)phase_trusted;
            T->edge[ne].du = T->u_med[ca] - T->u_med[cb];
            T->edge[ne].family = chart_group[ca] == chart_group[cb] &&
                                 fam[chart_group[ca]];
            if (T->edge[ne].family) {
                T->edge[ne].w = (int64_t)floor(
                    (double)w_lat * opts->family_w + 0.5);
                T->edge[ne].w_unhappy = T->edge[ne].w;
            } else {
                T->edge[ne].w = w_lat;
                T->edge[ne].w_unhappy = w_lat_un;
            }
            ne++;
        }
    }
    T->nedges = ne;
    T->erow = (size_t *)ARENA_CALLOC(arena, ncharts + 1, sizeof(size_t));
    for (size_t e = 0; e < ne; e++) {
        T->erow[T->edge[e].a + 1]++;
        T->erow[T->edge[e].b + 1]++;
    }
    for (size_t c = 0; c < ncharts; c++) T->erow[c + 1] += T->erow[c];
    T->eadj = (int32_t *)ARENA_ALLOC(
        arena, (T->erow[ncharts] ? T->erow[ncharts] : 1) * sizeof(int32_t));
    {
        size_t *fill = (size_t *)ARENA_ALLOC(arena, ncharts * sizeof(size_t));
        for (size_t c = 0; c < ncharts; c++) fill[c] = T->erow[c];
        for (size_t e = 0; e < ne; e++) {
            T->eadj[fill[T->edge[e].a]++] = (int32_t)e;
            T->eadj[fill[T->edge[e].b]++] = (int32_t)e;
        }
    }

    /* --- neighbourhoods: components of the LATERAL graph only ------------- */
    /* One contiguous same-wind jigsaw region.  The graded repair a collapsed
     * spiral needs (wind w by +1, wind w+1 by +2, ...) tears WELDS between
     * winds but never a lateral -- laterals only bind same-wind pairs -- so
     * the unit that must ride whole is the lateral cluster, with its weld
     * boundary priced like any other move.  (Clustering on weld+lateral
     * together was tried first and rebuilt the monster component inside the
     * solver: its compound moves looked collision-free and drifted the whole
     * field.) */
    {
        UnionFind uf = UF_new(arena, (int32_t)ncharts);
        for (size_t e = T->first_lateral; e < T->nedges; e++)
            if (T->edge[e].target == 0)
                uf_union(&uf, T->edge[e].a, T->edge[e].b);
        int32_t *root_id = (int32_t *)ARENA_ALLOC(
            arena, ncharts * sizeof(int32_t));
        for (size_t c = 0; c < ncharts; c++) root_id[c] = -1;
        T->nbid = (int32_t *)ARENA_ALLOC(arena, ncharts * sizeof(int32_t));
        int32_t nn = 0;
        for (size_t c = 0; c < ncharts; c++) {
            int32_t root = uf_find(&uf, (int32_t)c);
            if (root_id[root] < 0) root_id[root] = nn++;
            T->nbid[c] = root_id[root];
        }
        T->n_nbhd = (size_t)nn;
        T->nbrow = (size_t *)ARENA_CALLOC(arena, T->n_nbhd + 1,
                                          sizeof(size_t));
        for (size_t c = 0; c < ncharts; c++) T->nbrow[T->nbid[c] + 1]++;
        for (size_t n = 0; n < T->n_nbhd; n++) T->nbrow[n + 1] += T->nbrow[n];
        T->nbmember = (int32_t *)ARENA_ALLOC(arena,
                                             ncharts * sizeof(int32_t));
        size_t *fill = (size_t *)ARENA_ALLOC(arena,
                                             T->n_nbhd * sizeof(size_t));
        for (size_t n = 0; n < T->n_nbhd; n++) fill[n] = T->nbrow[n];
        for (size_t c = 0; c < ncharts; c++)
            T->nbmember[fill[T->nbid[c]]++] = (int32_t)c;
    }

    /* --- per-depth tables ------------------------------------------------- */
    /* The step for k winds is evaluated at the chart's U-LADDER position
     * (w_phi0 -- where its u currently SITS), not at its radius: a
     * mis-gauged chart must cross the arc between its current u-wind and
     * the destination, and that differs from the radius-wind arc by the
     * pitch growth of the mis-gauge.  Pooled to the lateral cluster's
     * median wind so same-wind jigsaw neighbours at equal depth move by
     * exactly the same amount -- per-chart wobble was smearing them apart. */
    double *w_pool = (double *)ARENA_ALLOC(arena, ncharts * sizeof(double));
    {
        double *scr = (double *)ARENA_ALLOC(arena, ncharts * sizeof(double));
        for (size_t c = 0; c < ncharts; c++) w_pool[c] = T->w_phi0[c];
        for (size_t n = 0; n < T->n_nbhd; n++) {
            size_t m0 = T->nbrow[n], m1 = T->nbrow[n + 1];
            if (m1 - m0 < 2) continue;
            for (size_t m = m0; m < m1; m++)
                scr[m - m0] = T->w_phi0[T->nbmember[m]];
            qsort(scr, m1 - m0, sizeof(double), aof_tabu_compare_double);
            double wm = scr[(m1 - m0) / 2];
            for (size_t m = m0; m < m1; m++)
                w_pool[T->nbmember[m]] = wm;
        }
        /* Families and explicitly locked phase-flat groups move RIGIDLY: one
         * pooled wind for the whole group, one shift per depth.  A
         * correct-but-misplaced region keeps its incoming layout exactly and
         * is only reordered in the strip. */
        for (size_t g = 0; g < ngroups; g++) {
            if (!fam[g] &&
                !(opts->phase_flat_lock && T->phase_flat_group[g]))
                continue;
            size_t m0 = T->grow[g], m1 = T->grow[g + 1];
            if (m1 - m0 < 2) continue;
            for (size_t m = m0; m < m1; m++)
                scr[m - m0] = T->w_phi0[T->gmember[m]];
            qsort(scr, m1 - m0, sizeof(double), aof_tabu_compare_double);
            double wm = scr[(m1 - m0) / 2];
            for (size_t m = m0; m < m1; m++)
                w_pool[T->gmember[m]] = wm;
        }
    }
    T->w_pool = w_pool;

    /* --- graded-rewind bins ----------------------------------------------- */
    /* Wind bins come from the WELD CHAIN, not from absolute radius: the axis
     * wanders ~85 vox (+/-7 winds), so pointwise rounded radius mis-bins
     * members and a graded candidate would land them wrong.  Between WELDED
     * neighbours the wander is common-mode and cancels, so integrating each
     * edge's radial climb dr/dr_dw over a BFS tree gives a wander-free
     * RELATIVE wind; lateral edges vote "same wind" exactly. */
    T->gbin = (int8_t *)ARENA_CALLOC(arena, ncharts, sizeof(int8_t));
    T->graded_group = (uint8_t *)ARENA_CALLOC(arena, ngroups,
                                              sizeof(uint8_t));
    if (fabs(T->dr_dw) > 1.0e-9) {
        double *w_rel = (double *)ARENA_ALLOC(arena,
                                              ncharts * sizeof(double));
        int32_t *bq = (int32_t *)ARENA_ALLOC(arena,
                                             ncharts * sizeof(int32_t));
        uint8_t *seen = (uint8_t *)ARENA_CALLOC(arena, ncharts,
                                                sizeof(uint8_t));
        for (size_t g = 0; g < ngroups; g++) {
            size_t m0 = T->grow[g], m1 = T->grow[g + 1];
            if (m1 - m0 < 3 || fam[g]) continue;
            if (ginfo == NULL || ginfo[g].intra_wrongwind == 0) continue;
            size_t head = 0, tail = 0;
            int32_t root = T->gmember[m0];
            seen[root] = 1;
            w_rel[root] = 0.0;
            bq[tail++] = root;
            while (head < tail) {
                int32_t c = bq[head++];
                for (size_t x = T->erow[c];
                     x < T->erow[(size_t)c + 1]; x++) {
                    const AofTabuEdge *e2 = &T->edge[T->eadj[x]];
                    int32_t o = e2->a == c ? e2->b : e2->a;
                    if (T->chart_group[o] != (int32_t)g || seen[o]) continue;
                    double dw = (size_t)T->eadj[x] < T->first_lateral
                              ? (T->r_med[o] - T->r_med[c]) / T->dr_dw
                              : 0.0;
                    seen[o] = 1;
                    w_rel[o] = w_rel[c] + dw;
                    bq[tail++] = o;
                }
            }
            double lo_w = DBL_MAX, hi_w = -DBL_MAX;
            int all_seen = 1;
            for (size_t m = m0; m < m1; m++) {
                int32_t c = T->gmember[m];
                if (!seen[c]) { all_seen = 0; continue; }
                if (w_rel[c] < lo_w) lo_w = w_rel[c];
                if (w_rel[c] > hi_w) hi_w = w_rel[c];
            }
            int32_t span_bins = all_seen && hi_w > lo_w
                              ? (int32_t)floor(hi_w - lo_w + 0.5) : 0;
            if (all_seen && span_bins >= 1 && span_bins <= 100) {
                for (size_t m = m0; m < m1; m++) {
                    int32_t c = T->gmember[m];
                    int32_t b = (int32_t)floor(w_rel[c] - lo_w + 0.5);
                    if (b < 0) b = 0;
                    T->gbin[c] = (int8_t)(b > 127 ? 127 : b);
                }
                T->graded_group[g] = 1;
            }
            for (size_t m = m0; m < m1; m++) seen[T->gmember[m]] = 0;
        }
    }

    /* Controlled replacement for the radius-chain hypothesis above.  The
     * radius construction telescopes to endpoint radius, so axial centreline
     * wander is accumulated rather than cancelled (black: nine false winds;
     * peach: six).  The lifted source phase removes only each chart's integer
     * branch gauge, and subtracting w_phi0 avoids re-unwinding turns already
     * present in the incoming atlas.  Keep signed, median-centred bins: the
     * graded move's free d is the group gauge. */
    if (opts->phase_graded) {
        memset(T->gbin, 0, ncharts * sizeof(int8_t));
        memset(T->graded_group, 0, ngroups * sizeof(uint8_t));

        for (size_t g = 0; g < ngroups; g++) {
            size_t m0 = T->grow[g], m1 = T->grow[g + 1];
            if (opts->phase_all_groups) {
                if (m1 - m0 < 2 || !T->phase_trusted_group[g]) continue;
            } else {
                if (m1 - m0 < 3 || fam[g]) continue;
                if (ginfo == NULL || ginfo[g].intra_wrongwind == 0) continue;
            }
            int lo = 127, hi = -127;
            for (size_t m = m0; m < m1; m++) {
                int32_t c = T->gmember[m];
                int b = (int)T->phase_bin[c];
                if (b < lo) lo = b;
                if (b > hi) hi = b;
            }
            if (hi - lo < 1 || hi - lo > 100) continue;
            for (size_t m = m0; m < m1; m++) {
                int32_t c = T->gmember[m];
                T->gbin[c] = T->phase_bin[c];
            }
            T->graded_group[g] = 1;
        }
    }

    /* The graded hypothesis and the neighbour objective must describe the
     * same topology.  k=d+gbin implies dk=gbin[a]-gbin[b]; retaining target 0
     * calls every intended winding boundary a tear and excludes it from the
     * continuous relayout.  These targets are node-potential differences, so
     * they are cycle-consistent by construction.  Lateral edges stay at zero:
     * their endpoints are same-wind neighbours. */
    if (opts->graded_targets) {
        for (size_t e = 0; e < T->first_lateral; e++) {
            int32_t ca = T->edge[e].a, cb = T->edge[e].b;
            int32_t g = chart_group[ca];
            if (chart_group[cb] != g || !T->graded_group[g]) continue;
            T->edge[e].target = (int32_t)T->gbin[ca] -
                                (int32_t)T->gbin[cb];
        }
    }


    aof_tabu_build_phase_regions(arena, T, opts);
    T->kmin = (int8_t *)ARENA_ALLOC(arena, ncharts * sizeof(int8_t));
    T->shift = (double *)ARENA_ALLOC(arena, ncharts * nk * sizeof(double));
    T->cell_shift = (int64_t *)ARENA_ALLOC(arena,
                                           ncharts * nk * sizeof(int64_t));
    T->radE = (int64_t *)ARENA_ALLOC(arena, ncharts * nk * sizeof(int64_t));
    for (size_t c = 0; c < ncharts; c++) {
        int lo = -(int)floor(T->w_phi0[c]);
        if (lo > 0) lo = 0;                 /* never forbid the initial 0 */
        if (lo < -T->span) lo = -T->span;
        T->kmin[c] = (int8_t)lo;
        double r_u = a + T->dr_dw * w_pool[c];   /* ladder radius at w_u */
        int gauge_ref = T->graded_group[T->chart_group[c]]
                      ? (int)T->gbin[c] : 0;
        double phi_u = AOF_2PI * sense * w_pool[c];
        for (int k = -T->span; k <= T->span; k++) {
            size_t row = c * nk + (size_t)(k + T->span);
            double s = 0.0;
            if (k != 0) {
                s = T->ladder_measured
                  ? T->ladder_step * AOF_2PI * (double)k *
                    (r_u + 0.5 * (double)k * T->dr_dw)
                  : CubeReg_deltaU(a, b, phi_u, (int32_t)(k * cal->sense));
            }
            T->shift[row] = s;
            T->cell_shift[row] = (int64_t)floor(s / T->cell + 0.5);
            double dev = T->w_phi0[c] + (double)k - T->w_rad[c];
            double ax = fabs(dev);
            double h = ax <= 1.0 ? 0.5 * ax * ax : ax - 0.5;
            if (h > 2.5) h = 2.5;
            double gd = (double)(k - gauge_ref);
            double gauge = 0.5 * gd * gd;
            T->radE[row] = (int64_t)floor(
                (opts->lambda_rad * h + opts->lambda_gauge * gauge) *
                AOF_TABU_FIX + 0.5);
        }
    }

    if (aof_tabu_build_potential_pairs(arena, T) != 0) return -1;

    /* --- occupancy hash + initial state ----------------------------------- */
    size_t total_cells = T->coff[ncharts];
    size_t cap = 64;
    while (cap < 4 * total_cells) cap <<= 1;
    if (hash_cap_override > 0) {
        cap = 64;
        while (cap < hash_cap_override) cap <<= 1;
        size_t min_cap = 64;
        while (min_cap < 2 * total_cells) min_cap <<= 1;
        if (cap < min_cap) cap = min_cap;
    }
    T->hash.key = (uint64_t *)ARENA_CALLOC(arena, cap, sizeof(uint64_t));
    T->hash.count = (uint32_t *)ARENA_CALLOC(arena, cap, sizeof(uint32_t));
    T->hash.cap = cap;
    size_t exact_cap = 64;
    /* Cache only AABB-surviving potential pairs.  Reserve a useful working
     * set, then let the existing half-load grow path expand on demand. */
    size_t exact_entries = T->npotential_pairs;
    size_t soft_entries = ncharts <= SIZE_MAX / 32 ? ncharts * 32 : SIZE_MAX;
    if (exact_entries > soft_entries) exact_entries = soft_entries;
    if (exact_entries < 32) exact_entries = 32;
    if (exact_entries > SIZE_MAX / 2) return -1;
    size_t exact_target = exact_entries * 2;
    while (exact_cap < exact_target) {
        if (exact_cap > SIZE_MAX / 2) return -1;
        exact_cap <<= 1;
    }
    T->exact.key = (uint64_t *)ARENA_CALLOC(
        arena, exact_cap, sizeof(uint64_t));
    T->exact.count = (uint32_t *)ARENA_CALLOC(
        arena, exact_cap, sizeof(uint32_t));
    T->exact.cap = exact_cap;

    T->k = (int8_t *)ARENA_CALLOC(arena, ncharts, sizeof(int8_t));
    T->happy = (uint8_t *)ARENA_ALLOC(arena, ncharts * sizeof(uint8_t));
    aof_tabu_rebuild(T);
    T->rebuilds = 0;          /* the initial fill is not a rebuild */
    T->e_ov = aof_tabu_exact_total(T);
    T->e_nb = aof_tabu_nb_total(T);
    T->e_rad = aof_tabu_rad_total(T);
    aof_tabu_exact_refresh_happy(T);
    return 0;
}

/*
 * The search itself.  Every iteration evaluates every admissible single-chart
 * move and every whole-group move (a group move has zero neighbour cost by
 * construction -- no weld edge crosses a group boundary -- so it is the rigid
 * placement's proven move class expressed inside the same energy, and it lets
 * the search re-place a 324-chart group without 324 uphill single moves).
 * The best steering-ranked move is ALWAYS accepted -- going uphill is how
 * tabu escapes -- while the tabu table forbids returning to a departed value
 * for `tenure` iterations, unless doing so would beat the best energy ever
 * seen (aspiration).  Stop on stagnation of the canonical best, never on the
 * collision proxy reaching zero (it cannot).
 */
static int aof_tabu_run(Arena_T arena, AofTabu *T,
                        const AtlasOverlapFixOptions *opts,
                        AtlasOverlapFixTabuMove *move_log, size_t move_cap,
                        size_t *out_nmoves,
                        int8_t *archive_k, int64_t *archive_e,
                        int64_t *archive_ov, int64_t *archive_nb,
                        int64_t *archive_rad, size_t archive_cap,
                        size_t *out_narchive,
                        AtlasOverlapFixStats *stats)
{
    size_t ncharts = T->ncharts;
    size_t nk = (size_t)T->nk;
    int tenure = opts->tabu_tenure > 0 ? opts->tabu_tenure : 1;
    int max_iters = opts->tabu_max_iters;
    int stall = opts->tabu_stall > 0 ? opts->tabu_stall : 1;
    int steer_on = opts->nb_unhappy_w < 1.0;

    int32_t *tabu_until = (int32_t *)ARENA_CALLOC(arena, ncharts * nk,
                                                  sizeof(int32_t));
    int32_t *gtabu_until = (int32_t *)ARENA_CALLOC(arena, T->ngroups * nk,
                                                   sizeof(int32_t));
    int32_t *ntabu_until = (int32_t *)ARENA_CALLOC(
        arena, (T->n_nbhd ? T->n_nbhd : 1) * nk, sizeof(int32_t));
    int32_t *dgtabu_until = (int32_t *)ARENA_CALLOC(
        arena, (T->ngroups ? T->ngroups : 1) * nk, sizeof(int32_t));
    int32_t *rtabu_until = (int32_t *)ARENA_CALLOC(
        arena, (T->n_region ? T->n_region : 1) * nk, sizeof(int32_t));
    uint8_t *graded_lock = (uint8_t *)ARENA_CALLOC(
        arena, T->ngroups ? T->ngroups : 1, sizeof(uint8_t));
    int32_t *exact_moved = (int32_t *)ARENA_CALLOC(
        arena, ncharts ? ncharts : 1, sizeof(int32_t));
    int32_t *exact_member = (int32_t *)ARENA_ALLOC(
        arena, (ncharts ? ncharts : 1) * sizeof(int32_t));
    int8_t *exact_k = (int8_t *)ARENA_ALLOC(
        arena, (ncharts ? ncharts : 1) * sizeof(int8_t));

    /* A trusted phase potential is topology, not a move suggestion: do not
     * let singles or mixed lateral neighborhoods dismantle it while waiting
     * for the first graded assignment to be selected. */
    if (opts->lock_graded)
        for (size_t g = 0; g < T->ngroups; g++)
            graded_lock[g] = T->graded_group[g];

    int64_t best_e = aof_tabu_energy(T);
    stats->tabu_energy_initial = (double)best_e / AOF_TABU_FIX;
    size_t best_iter = 0, nmoves = 0;
    /* Archive EVERY accepted state.  Search is capped, so this is bounded
     * (2001 x 638 bytes in the production case), and no coherent early state
     * can disappear through a four-slot circular buffer again. */
    memcpy(&archive_k[0], T->k, ncharts);
    archive_e[0] = best_e;
    archive_ov[0] = T->e_ov;
    archive_nb[0] = T->e_nb;
    archive_rad[0] = T->e_rad;

    for (int iter = 1; iter <= max_iters; iter++) {
        /* Iteration heartbeat: each iteration evaluates every chart at
         * every candidate depth with an exact SAT test, so on a large
         * grid one iteration is seconds and 2000 of them is hours. Print
         * the rate so the cost is visible instead of looking like a hang. */
        if (iter == 1 || (iter % 10) == 0) {
            double dt = ves_clock_sec() - g_aof_t0;
            fprintf(stderr,
                    "[overlap_fix]     tabu iter %d/%d, best e=%.1f @%zu, "
                    "%zu moves, t=%.1fs (%.2f s/iter)\n",
                    iter, max_iters, (double)best_e / AOF_TABU_FIX,
                    best_iter, nmoves, dt, dt / (double)iter);
            fflush(stderr);
        }
        if (iter % AOF_TABU_EPOCH == 0 && steer_on)
            aof_tabu_exact_refresh_happy(T);
        if (iter % AOF_TABU_AUDIT_EVERY == 0) {
            int64_t exact = aof_tabu_exact_total(T);
            if (exact != T->e_ov || aof_tabu_nb_total(T) != T->e_nb ||
                aof_tabu_rad_total(T) != T->e_rad)
                return -1;
        }

        int64_t E = aof_tabu_energy(T);
        int have = 0;
        int64_t bs = 0, b_pairs = 0, b_ov = 0, b_nb = 0, b_rad = 0;
        int bkind = 0, bfrom = 0, bto = 0;
        int32_t bid = 0;

        for (size_t c = 0; c < ncharts; c++) {
            if (!T->active[c]) continue;
            int32_t g = T->chart_group[c];
            if ((opts->phase_regions && T->region_id[g] >= 0) ||
                (opts->lock_graded && graded_lock[g] == 1) ||
                (opts->phase_flat_lock && T->phase_flat_group[g]))
                continue;
            int cur = T->k[c];
            exact_member[0] = (int32_t)c;
            for (int k = T->kmin[c]; k <= T->span; k++) {
                if (k == cur) continue;
                exact_k[0] = (int8_t)k;
                int64_t d_pairs = aof_tabu_exact_delta(
                    T, exact_member, exact_k, 1, exact_moved);
                int64_t d_ov = T->Lov * d_pairs;
                int64_t d_nb = aof_tabu_nb_delta(T, (int32_t)c, k, 0);
                int64_t d_rad =
                    T->radE[c * nk + (size_t)(k + T->span)] -
                    T->radE[c * nk + (size_t)(cur + T->span)];
                int64_t d_can = d_ov + d_nb + d_rad;
                int banned =
                    tabu_until[c * nk + (size_t)(k + T->span)] >= iter;
                if (banned && !(E + d_can < best_e)) continue;
                int64_t d_st = steer_on
                    ? d_ov + aof_tabu_nb_delta(T, (int32_t)c, k, 1) + d_rad
                    : d_can;
                if (!have || d_st < bs) {
                    have = 1; bs = d_st;
                    bkind = 0; bid = (int32_t)c; bfrom = cur; bto = k;
                    b_pairs = d_pairs;
                    b_ov = d_ov; b_nb = d_nb; b_rad = d_rad;
                }
            }
        }

        for (size_t g = 0; g < T->ngroups; g++) {
            size_t m0 = T->grow[g], m1 = T->grow[g + 1];
            if (opts->phase_regions && T->region_id[g] >= 0) continue;
            size_t nmember = 0;
            for (size_t m = m0; m < m1; m++) {
                int32_t c = T->gmember[m];
                if (T->active[c]) exact_member[nmember++] = c;
            }
            if (nmember < 2) continue;       /* singleton == chart move */
            for (int d = -T->span; d <= T->span; d++) {
                if (d == 0) continue;
                int ok = 1;
                int64_t d_rad = 0;
                for (size_t i = 0; i < nmember && ok; i++) {
                    int32_t c = exact_member[i];
                    int knew = T->k[c] + d;
                    exact_k[i] = (int8_t)knew;
                    if (knew < T->kmin[c] || knew > T->span) ok = 0;
                    else d_rad +=
                        T->radE[(size_t)c * nk + (size_t)(knew + T->span)] -
                        T->radE[(size_t)c * nk + (size_t)(T->k[c] + T->span)];
                }
                if (!ok) continue;
                int64_t d_pairs = aof_tabu_exact_delta(
                    T, exact_member, exact_k, nmember, exact_moved);
                int64_t d_ov = T->Lov * d_pairs;
                /* Weld edges never cross a group boundary, but LATERAL
                 * edges can: charge the ones with exactly one endpoint
                 * moving.  Edges wholly inside the group are invariant. */
                int64_t d_nb = 0, d_nb_st = 0;
                for (size_t i = 0; i < nmember; i++) {
                    int32_t c = exact_member[i];
                    for (size_t x = T->erow[c];
                         x < T->erow[(size_t)c + 1]; x++) {
                        const AofTabuEdge *e2 = &T->edge[T->eadj[x]];
                        int32_t other = e2->a == c ? e2->b : e2->a;
                        if (T->chart_group[other] == (int32_t)g &&
                            T->active[other])
                            continue;
                        int ka = e2->a == c ? T->k[c] + d : T->k[e2->a];
                        int kb = e2->b == c ? T->k[c] + d : T->k[e2->b];
                        int dn = ka - kb - e2->target;
                        int dc = T->k[e2->a] - T->k[e2->b] - e2->target;
                        int64_t dd = (int64_t)((dn < 0 ? -dn : dn) -
                                               (dc < 0 ? -dc : dc));
                        d_nb += e2->w * dd;
                        d_nb_st += (T->happy[other] ? e2->w
                                                    : e2->w_unhappy) * dd;
                    }
                }
                int64_t d_can = d_ov + d_nb + d_rad;
                int banned =
                    gtabu_until[g * nk + (size_t)(d + T->span)] >= iter;
                if (banned && !(E + d_can < best_e)) continue;
                int64_t d_st = steer_on ? d_ov + d_nb_st + d_rad : d_can;
                if (!have || d_st < bs) {
                    have = 1; bs = d_st;
                    bkind = 1; bid = (int32_t)g; bfrom = 0; bto = d;
                    b_pairs = d_pairs;
                    b_ov = d_ov; b_nb = d_nb; b_rad = d_rad;
                }
            }
        }

        /* Whole-neighbourhood moves: one contiguous same-wind jigsaw region
         * travels as a unit.  Laterals never cross the boundary (they define
         * it); welds to other winds do, and are priced exactly like any
         * other move's tear.  The size cap guards against degenerate huge
         * clusters whose cross-cluster collision is too empty to steer. */
        size_t nbhd_cap = ncharts / 8 > 8 ? ncharts / 8 : 8;
        for (size_t n = 0; n < T->n_nbhd; n++) {
            size_t m0 = T->nbrow[n], m1 = T->nbrow[n + 1];
            if (m1 - m0 < 2 || m1 - m0 > nbhd_cap) continue;
            int all_active = 1;
            for (size_t m = m0; m < m1 && all_active; m++)
                all_active = T->active[T->nbmember[m]] != 0;
            if (!all_active) continue;
            if (opts->phase_regions || opts->lock_graded || opts->phase_flat_lock) {
                int locked = 0;
                for (size_t m = m0; m < m1 && !locked; m++) {
                    int32_t c = T->nbmember[m];
                    int32_t g = T->chart_group[c];
                    locked = (opts->phase_regions && T->region_id[g] >= 0) ||
                             (opts->lock_graded && graded_lock[g] == 1) ||
                             (opts->phase_flat_lock &&
                              T->phase_flat_group[g]);
                }
                if (locked) continue;
            }
            size_t nmember = m1 - m0;
            for (size_t i = 0; i < nmember; i++) exact_member[i] = T->nbmember[m0 + i];
            for (int d = -T->span; d <= T->span; d++) {
                if (d == 0) continue;
                int ok = 1;
                int64_t d_rad = 0;
                for (size_t m = m0; m < m1 && ok; m++) {
                    int32_t c = T->nbmember[m];
                    int knew = T->k[c] + d;
                    exact_k[m - m0] = (int8_t)knew;
                    if (knew < T->kmin[c] || knew > T->span) ok = 0;
                    else d_rad +=
                        T->radE[(size_t)c * nk + (size_t)(knew + T->span)] -
                        T->radE[(size_t)c * nk + (size_t)(T->k[c] + T->span)];
                }
                if (!ok) continue;
                int64_t d_pairs = aof_tabu_exact_delta(
                    T, exact_member, exact_k, nmember, exact_moved);
                int64_t d_ov = T->Lov * d_pairs;
                int64_t d_nb = 0, d_nb_st = 0;
                for (size_t m = m0; m < m1; m++) {
                    int32_t c = T->nbmember[m];
                    for (size_t x = T->erow[c];
                         x < T->erow[(size_t)c + 1]; x++) {
                        const AofTabuEdge *e2 = &T->edge[T->eadj[x]];
                        int32_t other = e2->a == c ? e2->b : e2->a;
                        if (T->nbid[other] == (int32_t)n) continue;
                        int ka = e2->a == c ? T->k[c] + d : T->k[e2->a];
                        int kb = e2->b == c ? T->k[c] + d : T->k[e2->b];
                        int dn = ka - kb - e2->target;
                        int dc = T->k[e2->a] - T->k[e2->b] - e2->target;
                        int64_t dd = (int64_t)((dn < 0 ? -dn : dn) -
                                               (dc < 0 ? -dc : dc));
                        d_nb += e2->w * dd;
                        d_nb_st += (T->happy[other] ? e2->w
                                                    : e2->w_unhappy) * dd;
                    }
                }
                int64_t d_can = d_ov + d_nb + d_rad;
                int banned =
                    ntabu_until[n * nk + (size_t)(d + T->span)] >= iter;
                if (banned && !(E + d_can < best_e)) continue;
                int64_t d_st = steer_on ? d_ov + d_nb_st + d_rad : d_can;
                if (!have || d_st < bs) {
                    have = 1; bs = d_st;
                    bkind = 2; bid = (int32_t)n; bfrom = 0; bto = d;
                    b_pairs = d_pairs;
                    b_ov = d_ov; b_nb = d_nb; b_rad = d_rad;
                }
            }
        }

        /* Graded rewinds: wind bin j of a collapsed group moves to depth
         * d + j in ONE move -- the unrolling that a per-chart search can only
         * reach through uphill cascades and usually scatters instead.  It
         * tears exactly the between-bin welds (priced below) and preserves
         * everything inside each bin. */
        for (size_t g = 0; g < T->ngroups; g++) {
            if (!T->graded_group[g]) continue;
            if (opts->phase_regions && T->region_id[g] >= 0) continue;
            if (opts->phase_flat_lock && T->phase_flat_group[g])
                continue;

            size_t m0 = T->grow[g], m1 = T->grow[g + 1];
            int all_active = 1;
            for (size_t m = m0; m < m1 && all_active; m++)
                all_active = T->active[T->gmember[m]] != 0;
            if (!all_active) continue;
            size_t nmember = m1 - m0;
            for (size_t i = 0; i < nmember; i++) exact_member[i] = T->gmember[m0 + i];
            for (int d = -T->span; d <= T->span; d++) {
                int ok = 1, any = 0;
                int64_t d_rad = 0;
                for (size_t m = m0; m < m1 && ok; m++) {
                    int32_t c = T->gmember[m];
                    int knew = d + (int)T->gbin[c];
                    exact_k[m - m0] = (int8_t)knew;
                    if (knew < T->kmin[c] || knew > T->span) { ok = 0; break; }
                    if (knew != T->k[c]) any = 1;
                    d_rad +=
                        T->radE[(size_t)c * nk + (size_t)(knew + T->span)] -
                        T->radE[(size_t)c * nk + (size_t)(T->k[c] + T->span)];
                }
                if (!ok || !any) continue;
                int64_t d_pairs = aof_tabu_exact_delta(
                    T, exact_member, exact_k, nmember, exact_moved);
                int64_t d_ov = T->Lov * d_pairs;
                int64_t d_nb = 0, d_nb_st = 0;
                for (size_t m = m0; m < m1; m++) {
                    int32_t c = T->gmember[m];
                    for (size_t x = T->erow[c];
                         x < T->erow[(size_t)c + 1]; x++) {
                        const AofTabuEdge *e2 = &T->edge[T->eadj[x]];
                        int32_t other = e2->a == c ? e2->b : e2->a;
                        int other_in = T->chart_group[other] == (int32_t)g;
                        if (other_in && other < c) continue;  /* once */
                        int na = T->chart_group[e2->a] == (int32_t)g
                               ? d + (int)T->gbin[e2->a] : T->k[e2->a];
                        int nb2 = T->chart_group[e2->b] == (int32_t)g
                               ? d + (int)T->gbin[e2->b] : T->k[e2->b];
                        int dn = na - nb2 - e2->target;
                        int dc = T->k[e2->a] - T->k[e2->b] - e2->target;
                        int64_t dd = (int64_t)((dn < 0 ? -dn : dn) -
                                               (dc < 0 ? -dc : dc));
                        d_nb += e2->w * dd;
                        d_nb_st += (T->happy[other] ? e2->w
                                                    : e2->w_unhappy) * dd;
                    }
                }
                int64_t d_can = d_ov + d_nb + d_rad;
                int banned =
                    dgtabu_until[g * nk + (size_t)(d + T->span)] >= iter;
                if (banned && !(E + d_can < best_e)) continue;
                int64_t d_st = steer_on ? d_ov + d_nb_st + d_rad : d_can;
                if (!have || d_st < bs) {
                    have = 1; bs = d_st;
                    bkind = 3; bid = (int32_t)g; bfrom = 0; bto = d;
                    b_pairs = d_pairs;
                    b_ov = d_ov; b_nb = d_nb; b_rad = d_rad;
                }
            }
        }

        /* Relative-gauge regions: strong cross-group lateral evidence fixes
         * each member group's offset, while one common base d remains free.
         * This is the move missing from independent group search: e.g. groups
         * 5/6/7 can close all 92 mutually consistent contacts without walking
         * any one group through a colliding intermediate placement. */
        for (size_t r = 0; r < T->n_region; r++) {
            size_t q0 = T->region_row[r], q1 = T->region_row[r + 1];
            size_t nmember = 0;
            int all_active = 1;
            for (size_t q = q0; q < q1; q++) {
                int32_t g = T->region_group[q];
                for (size_t m = T->grow[g]; m < T->grow[g + 1]; m++) {
                    int32_t c = T->gmember[m];
                    if (!T->active[c]) all_active = 0;
                    exact_member[nmember++] = c;
                }
            }
            if (!all_active) continue;

            /* Exact SAT evaluates the whole region directly. */

            for (int d = -T->span; d <= T->span; d++) {
                int ok = 1, any = 0;
                int64_t d_rad = 0;
                size_t ei = 0;
                for (size_t q = q0; q < q1 && ok; q++) {
                    int32_t g = T->region_group[q];
                    for (size_t m = T->grow[g]; m < T->grow[g + 1]; m++) {
                        int32_t c = T->gmember[m];
                        int knew = d + T->region_offset[g] +
                                   (int)T->gbin[c];
                        if (knew < T->kmin[c] || knew > T->span) {
                            ok = 0;
                            break;
                        }
                        exact_k[ei++] = (int8_t)knew;
                        if (knew != T->k[c]) any = 1;
                        d_rad += T->radE[(size_t)c * nk +
                                              (size_t)(knew + T->span)] -
                                 T->radE[(size_t)c * nk +
                                              (size_t)(T->k[c] + T->span)];
                    }
                }
                if (!ok || !any) continue;
                int64_t d_pairs = aof_tabu_exact_delta(
                    T, exact_member, exact_k, nmember, exact_moved);
                int64_t d_ov = T->Lov * d_pairs;
                int64_t d_nb = 0, d_nb_st = 0;
                for (size_t q = q0; q < q1; q++) {
                    int32_t g = T->region_group[q];
                    for (size_t m = T->grow[g]; m < T->grow[g + 1]; m++) {
                        int32_t c = T->gmember[m];
                        for (size_t x = T->erow[c];
                             x < T->erow[(size_t)c + 1]; x++) {
                            const AofTabuEdge *e2 = &T->edge[T->eadj[x]];
                            int32_t other = e2->a == c ? e2->b : e2->a;
                            int32_t go = T->chart_group[other];
                            int other_in = T->region_id[go] == (int32_t)r;
                            if (other_in && other < c) continue;
                            int32_t ga = T->chart_group[e2->a];
                            int32_t gb = T->chart_group[e2->b];
                            int na = T->region_id[ga] == (int32_t)r
                                   ? d + T->region_offset[ga] +
                                     (int)T->gbin[e2->a] : T->k[e2->a];
                            int nb2 = T->region_id[gb] == (int32_t)r
                                    ? d + T->region_offset[gb] +
                                      (int)T->gbin[e2->b] : T->k[e2->b];
                            int dn = na - nb2 - e2->target;
                            int dc = T->k[e2->a] - T->k[e2->b] - e2->target;
                            int64_t dd = (int64_t)((dn < 0 ? -dn : dn) -
                                                   (dc < 0 ? -dc : dc));
                            d_nb += e2->w * dd;
                            d_nb_st += (T->happy[other] ? e2->w
                                                        : e2->w_unhappy) * dd;
                        }
                    }
                }
                int64_t d_can = d_ov + d_nb + d_rad;
                int banned = rtabu_until[r * nk +
                                  (size_t)(d + T->span)] >= iter;
                if (banned && !(E + d_can < best_e)) continue;
                int64_t d_st = steer_on ? d_ov + d_nb_st + d_rad : d_can;
                if (!have || d_st < bs) {
                    have = 1; bs = d_st;
                    bkind = 4; bid = (int32_t)r; bfrom = 0; bto = d;
                    b_pairs = d_pairs;
                    b_ov = d_ov; b_nb = d_nb; b_rad = d_rad;
                }
            }
        }

        if (!have) break;   /* everything banned and nothing aspirates */

        if (bkind == 0) {
            T->k[bid] = (int8_t)bto;
            tabu_until[(size_t)bid * nk + (size_t)(bfrom + T->span)] =
                iter + tenure;
            stats->tabu_moves_chart++;
        } else if (bkind == 4) {
            size_t q0 = T->region_row[bid], q1 = T->region_row[bid + 1];
            int old_base = 0, have_old_base = 0, old_conforming = 1;
            for (size_t q = q0; q < q1; q++) {
                int32_t g = T->region_group[q];
                for (size_t m = T->grow[g]; m < T->grow[g + 1]; m++) {
                    int32_t c = T->gmember[m];
                    int base = T->k[c] - T->region_offset[g] -
                               (int)T->gbin[c];
                    if (!have_old_base) {
                        old_base = base;
                        have_old_base = 1;
                    } else if (base != old_base) {
                        old_conforming = 0;
                    }
                }
            }
            for (size_t q = q0; q < q1; q++) {
                int32_t g = T->region_group[q];
                for (size_t m = T->grow[g]; m < T->grow[g + 1]; m++) {
                    int32_t c = T->gmember[m];
                    T->k[c] = (int8_t)(bto + T->region_offset[g] +
                                        (int)T->gbin[c]);
                }
            }
            if (old_conforming && old_base >= -T->span &&
                old_base <= T->span)
                rtabu_until[(size_t)bid * nk +
                            (size_t)(old_base + T->span)] = iter + tenure;
            stats->tabu_moves_region++;
        } else {
            const int32_t *member = bkind == 2 ? T->nbmember : T->gmember;
            size_t m0 = bkind == 2 ? T->nbrow[bid] : T->grow[bid];
            size_t m1 = bkind == 2 ? T->nbrow[bid + 1] : T->grow[bid + 1];
            for (size_t m = m0; m < m1; m++) {
                int32_t c = member[m];
                if (bkind == 1 && !T->active[c]) continue;
                T->k[c] = bkind == 3 ? (int8_t)(bto + (int)T->gbin[c])
                                     : (int8_t)(T->k[c] + bto);
            }
            if (bkind == 1) {
                gtabu_until[(size_t)bid * nk + (size_t)(-bto + T->span)] =
                    iter + tenure;
                stats->tabu_moves_group++;
            } else if (bkind == 2) {
                ntabu_until[(size_t)bid * nk + (size_t)(-bto + T->span)] =
                    iter + tenure;
                stats->tabu_moves_nbhd++;
            } else {
                dgtabu_until[(size_t)bid * nk + (size_t)(bto + T->span)] =
                    iter + tenure;
                if (opts->lock_graded && opts->postgraded_refine)
                    graded_lock[bid] = 2;
                stats->tabu_moves_graded++;
            }
        }
        T->e_ov += b_pairs;
        T->e_nb += b_nb;
        T->e_rad += b_rad;

        int64_t e_after = aof_tabu_energy(T);
        if (nmoves < move_cap) {
            AtlasOverlapFixTabuMove *mv = &move_log[nmoves];
            mv->iter = iter;
            mv->kind = bkind;
            mv->id = bid;
            mv->from = bfrom;
            mv->to = bto;
            mv->d_overlap = (double)b_ov / AOF_TABU_FIX;
            mv->d_neighbour = (double)b_nb / AOF_TABU_FIX;
            mv->d_radial = (double)b_rad / AOF_TABU_FIX;
            mv->energy_after = (double)e_after / AOF_TABU_FIX;
        }
        nmoves++;
        size_t aslot = nmoves;
        if (aslot < archive_cap) {
            memcpy(&archive_k[aslot * ncharts], T->k, ncharts);
            archive_e[aslot] = e_after;
            archive_ov[aslot] = T->e_ov;
            archive_nb[aslot] = T->e_nb;
            archive_rad[aslot] = T->e_rad;
        }

        stats->tabu_iterations = (size_t)iter;

        if (e_after < best_e) {
            best_e = e_after;
            best_iter = (size_t)iter;
        }
        if ((size_t)iter - best_iter >= (size_t)stall) break;
    }

    stats->tabu_best_iter = best_iter;
    stats->tabu_hash_rebuilds = T->rebuilds;
    *out_nmoves = nmoves < move_cap ? nmoves : move_cap;
    stats->tabu_exact_queries = T->exact_queries;
    stats->tabu_exact_cache_entries = T->exact.used;
    stats->tabu_exact_cache_hits = T->exact_cache_hits;
    stats->tabu_exact_cache_misses = T->exact_cache_misses;
    stats->tabu_exact_aabb_rejects = T->exact_aabb_rejects;
    stats->tabu_exact_sat_tests = T->exact_sat_tests;
    *out_narchive = 1 + (nmoves < archive_cap - 1 ? nmoves : archive_cap - 1);
    return 0;
}

static void aof_tabu_verify_add(size_t index, const int8_t *archive_k,
                                size_t ncharts, size_t *verify,
                                size_t *nverify, int force)
{
    if (*nverify >= AOF_TABU_VERIFY_MAX) return;
    const int8_t *candidate = &archive_k[index * ncharts];
    size_t distinct_min = ncharts / 128;
    if (distinct_min < 1) distinct_min = 1;
    for (size_t v = 0; v < *nverify; v++) {
        const int8_t *prior = &archive_k[verify[v] * ncharts];
        size_t different = 0;
        for (size_t c = 0; c < ncharts; c++)
            if (candidate[c] != prior[c]) different++;
        if (different == 0 || (!force && different < distinct_min)) return;
    }
    verify[(*nverify)++] = index;
}

/* Down-select the full move archive without chronological eviction.  Pinned
 * states cover the input, canonical best, endpoint, every graded transition,
 * and the extrema of the proxy/topology tradeoff.  The remaining slots sample
 * that tradeoff and time uniformly, rejecting near-duplicate assignments. */
static size_t aof_tabu_choose_verify(const AofTabu *T,
                                     const AtlasOverlapFixTabuMove *moves,
                                     size_t nmoves, const int8_t *archive_k,
                                     const int64_t *archive_e,
                                     const int64_t *archive_ov,
                                     const int64_t *archive_nb,
                                     const int64_t *archive_rad,
                                     size_t narchive, size_t best_iter,
                                     size_t verify[AOF_TABU_VERIFY_MAX])
{
    size_t nverify = 0;
    aof_tabu_verify_add(0, archive_k, T->ncharts,
                        verify, &nverify, 1);
    if (best_iter < narchive)
        aof_tabu_verify_add(best_iter, archive_k, T->ncharts,
                            verify, &nverify, 1);
    if (narchive > 1)
        aof_tabu_verify_add(narchive - 1, archive_k, T->ncharts,
                            verify, &nverify, 1);
    for (size_t m = 0; m < nmoves && m + 1 < narchive; m++)
        if (moves[m].kind == 3)
            aof_tabu_verify_add(m + 1, archive_k, T->ncharts,
                                verify, &nverify, 1);

    size_t min_e = 0, min_ov = 0;
    for (size_t a = 1; a < narchive; a++) {
        if (archive_e[a] < archive_e[min_e]) min_e = a;
        if (archive_ov[a] < archive_ov[min_ov]) min_ov = a;
    }
    aof_tabu_verify_add(min_e, archive_k, T->ncharts,
                        verify, &nverify, 1);
    aof_tabu_verify_add(min_ov, archive_k, T->ncharts,
                        verify, &nverify, 1);

    int64_t ov_lo = archive_ov[min_ov], ov_hi = archive_ov[0];
    int64_t ov_mid = ov_lo + (ov_hi - ov_lo) / 2;
    size_t best_nb = (size_t)-1, best_rad = (size_t)-1;
    for (size_t a = 1; a < narchive; a++) {
        if (archive_ov[a] > ov_mid) continue;
        if (best_nb == (size_t)-1 || archive_nb[a] < archive_nb[best_nb])
            best_nb = a;
        if (best_rad == (size_t)-1 || archive_rad[a] < archive_rad[best_rad])
            best_rad = a;
    }
    if (best_nb != (size_t)-1)
        aof_tabu_verify_add(best_nb, archive_k, T->ncharts,
                            verify, &nverify, 1);
    if (best_rad != (size_t)-1)
        aof_tabu_verify_add(best_rad, archive_k, T->ncharts,
                            verify, &nverify, 1);

    for (size_t b = 0; b < 16 && nverify < AOF_TABU_VERIFY_MAX; b++) {
        int64_t lo = ov_lo + (ov_hi - ov_lo) * (int64_t)b / 16;
        int64_t hi = ov_lo + (ov_hi - ov_lo) * (int64_t)(b + 1) / 16;
        size_t best = (size_t)-1;
        for (size_t a = 0; a < narchive; a++) {
            if (archive_ov[a] < lo || (b < 15 && archive_ov[a] >= hi))
                continue;
            if (best == (size_t)-1 ||
                archive_nb[a] + archive_rad[a] <
                archive_nb[best] + archive_rad[best])
                best = a;
        }
        if (best != (size_t)-1)
            aof_tabu_verify_add(best, archive_k, T->ncharts,
                                verify, &nverify, 0);
    }
    for (size_t q = 1; q < AOF_TABU_VERIFY_MAX * 4 &&
                       nverify < AOF_TABU_VERIFY_MAX; q++) {
        size_t a = narchive > 1
                 ? q * (narchive - 1) / (AOF_TABU_VERIFY_MAX * 4 - 1) : 0;
        aof_tabu_verify_add(a, archive_k, T->ncharts,
                            verify, &nverify, 0);
    }
    for (size_t a = 0; a < narchive && nverify < AOF_TABU_VERIFY_MAX; a++)
        aof_tabu_verify_add(a, archive_k, T->ncharts,
                            verify, &nverify, 0);
    return nverify;
}

/* ========================================================================== */
/* Relayout: continuous re-registration with the windings held fixed          */
/*                                                                            */
/* Tabu owns the integers; this stage owns the metric.  Given the winning     */
/* depth assignment, every surviving relationship becomes a difference        */
/* equation on one continuous u shift per chart (AtlasRegister):              */
/*   - an INTACT weld (dk == target) becomes a closure spring pulling its     */
/*     contact mismatch to zero;                                              */
/*   - a SAME-DEPTH lateral pair becomes a keep-together spring;              */
/*   - a radial pair m winds apart is pulled to m measured ladder steps --    */
/*     the constraint that stops the continuous stage from re-collapsing      */
/*     the winds tabu just separated;                                         */
/*   - TORN welds and split laterals are deliberately absent: they are the    */
/*     winding decision, and springs across them would fight it.              */
/* Priors of zero keep each evidence island where tabu left it on average.    */
/* ========================================================================== */

typedef struct {
    int32_t chart;
    double u;
} AofRelayoutRank;

static int aof_relayout_rank_compare(const void *left, const void *right)
{
    const AofRelayoutRank *a = (const AofRelayoutRank *)left;
    const AofRelayoutRank *b = (const AofRelayoutRank *)right;
    if (a->u < b->u) return -1;
    if (a->u > b->u) return 1;
    if (a->chart < b->chart) return -1;
    if (a->chart > b->chart) return 1;
    return 0;
}

static void aof_relayout_edge_rms(const AtlasRegisterEdge *edge, size_t nedges,
                                  AtlasRegisterEdgeKind kind,
                                  const double *shift, double *out_rms)
{
    double square = 0.0, weight_sum = 0.0;
    for (size_t e = 0; e < nedges; e++) {
        if (edge[e].kind != kind || !(edge[e].weight > 0.0)) continue;
        double residual = shift[edge[e].component1] -
                          shift[edge[e].component0] - edge[e].target;
        double weight = edge[e].weight *
                        (kind == ATLAS_REGISTER_RADIAL ? 0.25 : 1.0);
        square += weight * residual * residual;
        weight_sum += weight;
    }
    *out_rms = weight_sum > 0.0 ? sqrt(square / weight_sum) : 0.0;
}

static double aof_relayout_width(const PieceSet *ps,
                                 const int32_t *vertex_chart,
                                 const double *u)
{
    double lo = DBL_MAX, hi = -DBL_MAX;
    for (size_t i = 0; i < ps->nv; i++) {
        if (vertex_chart[i] < 0 || !isfinite(u[i])) continue;
        if (u[i] < lo) lo = u[i];
        if (u[i] > hi) hi = u[i];
    }
    return lo <= hi ? hi - lo : 0.0;
}

static int aof_relayout(Arena_T arena, const PieceSet *ps,
                        const AofTabu *T, const int8_t *win_k,
                        const int32_t *vertex_chart,
                        const uint8_t *face_keep, const double *out_v,
                        double seam_cap, double *out_u,
                        double *out_metric_shift, double *out_chart_shift,
                        size_t *out_chart_rank,
                        AtlasOverlapFixTabuReport *report,
                        AtlasOverlapFixStats *stats)
{
    size_t ncharts = T->ncharts;
    size_t nk = (size_t)T->nk;
    stats->relayout_failed = 1;

    /* Post-tabu chart medians. */
    double *pu = (double *)ARENA_ALLOC(arena, ncharts * sizeof(double));
    for (size_t c = 0; c < ncharts; c++) {
        double s = win_k[c] != 0
                 ? T->shift[c * nk + (size_t)(win_k[c] + T->span)] : 0.0;
        pu[c] = T->u_med[c] + s;
    }

    size_t cap = T->nedges + ncharts * 8;
    AtlasRegisterEdge *redge = (AtlasRegisterEdge *)ARENA_ALLOC(
        arena, (cap ? cap : 1) * sizeof(*redge));
    size_t nre = 0;

    for (size_t e = 0; e < T->nedges; e++) {
        const AofTabuEdge *te = &T->edge[e];
        int dk = win_k[te->a] - win_k[te->b];
        double sa = win_k[te->a] != 0
            ? T->shift[(size_t)te->a * nk + (size_t)(win_k[te->a] + T->span)]
            : 0.0;
        double sb = win_k[te->b] != 0
            ? T->shift[(size_t)te->b * nk + (size_t)(win_k[te->b] + T->span)]
            : 0.0;
        double du_now = te->du + sa - sb;
        /* Same-bin contacts retain their good incoming metric offset.  A
         * nonzero target is an intentional winding cut, so hold the discrete
         * one-step jump instead of collapsing adjacent sheets back together.
         * Exact collision bounds below decide any additional admissible slide. */
        if (e < T->first_lateral) {
            if (dk != te->target) continue;      /* torn by design */
            redge[nre].component0 = te->a;
            redge[nre].component1 = te->b;
            redge[nre].target = te->target == 0 ? du_now - te->du : 0.0;
            redge[nre].weight = (double)te->w / AOF_TABU_FIX;
            redge[nre].robust_scale = 8.0;
            redge[nre].kind = ATLAS_REGISTER_WELD;
            nre++;
            stats->relayout_edges_weld++;
        } else {
            if (dk != te->target) continue;      /* split by design */
            redge[nre].component0 = te->a;
            redge[nre].component1 = te->b;
            redge[nre].target = te->target == 0 ? du_now - te->du : 0.0;
            redge[nre].weight = (double)te->w / AOF_TABU_FIX;
            redge[nre].robust_scale = 16.0;
            redge[nre].kind = ATLAS_REGISTER_LOCAL;
            nre++;
            stats->relayout_edges_lateral++;
        }
    }

    /* Radial separation at the measured ladder step. */
    if (T->ladder_measured && T->crow != NULL && ps->cube_org != NULL) {
        double p = fabs(T->dr_dw);
        for (size_t q1 = 0; q1 < ps->n_cubes; q1++)
        for (size_t q2 = q1; q2 < ps->n_cubes; q2++) {
            if (labs(ps->cube_org[q1][0] - ps->cube_org[q2][0]) > 128 ||
                labs(ps->cube_org[q1][1] - ps->cube_org[q2][1]) > 128 ||
                labs(ps->cube_org[q1][2] - ps->cube_org[q2][2]) > 128)
                continue;
            for (size_t i = T->crow[q1]; i < T->crow[q1 + 1]; i++) {
                size_t j0 = q1 == q2 ? i + 1 : T->crow[q2];
                for (size_t j = j0; j < T->crow[q2 + 1]; j++) {
                    int32_t ca = T->cmember[i], cb = T->cmember[j];
                    double dr = T->r_med[cb] - T->r_med[ca];
                    double m_real = fabs(dr) / p;
                    int m = (int)floor(m_real + 0.5);
                    if (m < 1 || m > 3) continue;
                    if (fabs(m_real - (double)m) > 0.3) continue;
                    if (nre >= cap) continue;
                    int32_t inner = dr > 0.0 ? ca : cb;
                    int32_t outer = dr > 0.0 ? cb : ca;
                    double r_mid = 0.5 * (T->r_med[ca] + T->r_med[cb]);
                    double target_du =
                        (double)m * T->ladder_step * AOF_2PI * r_mid;
                    double err = (pu[outer] - pu[inner]) - target_du;
                    /* GUARD, not demand: stabilize only pairs already near
                     * their correct separation -- near in ABSOLUTE terms, a
                     * fraction of ONE wind step, never of the m-wind span
                     * (0.35 * a 3-wind target admitted +/-2900 vox "guards"
                     * that moved charts by thousands; the unguarded version
                     * before that sheared the welds to 2000+ vox rms).
                     * Still-collapsed pairs are the next tabu round's work. */
                    double step_one =
                        fabs(T->ladder_step) * AOF_2PI * r_mid;
                    if (fabs(err) > 0.15 * step_one) continue;
                    redge[nre].component0 = inner;
                    redge[nre].component1 = outer;
                    redge[nre].target = -err;
                    redge[nre].weight = 1.0;
                    redge[nre].robust_scale = 50.0;
                    redge[nre].kind = ATLAS_REGISTER_RADIAL;
                    nre++;
                    stats->relayout_edges_radial++;
                }
            }
        }
    }

    double *prior_shift = (double *)ARENA_CALLOC(arena, ncharts,
                                                 sizeof(double));
    double *prior_weight = (double *)ARENA_ALLOC(arena,
                                                 ncharts * sizeof(double));
    size_t *chart_faces = (size_t *)ARENA_CALLOC(
        arena, ncharts, sizeof(*chart_faces));
    size_t total_chart_faces = 0;
    double edge_length_sum = 0.0;
    size_t edge_length_count = 0;
    for (size_t f = 0; f < ps->nf; f++) {
        if (face_keep != NULL && !face_keep[f]) continue;
        int32_t c = vertex_chart[ps->faces[f * 3]];
        if (c < 0 || (size_t)c >= ncharts) return -1;
        chart_faces[c]++;
        total_chart_faces++;
        for (int k = 0; k < 3; k++) {
            int32_t a = ps->faces[f * 3 + (size_t)k];
            int32_t b = ps->faces[f * 3 + (size_t)((k + 1) % 3)];
            double du = out_u[a] - out_u[b];
            double dv = out_v[a] - out_v[b];
            double length = sqrt(du * du + dv * dv);
            if (!isfinite(length)) return -1;
            edge_length_sum += length;
            edge_length_count++;
        }
    }
    double mean_faces = ncharts > 0
                      ? (double)total_chart_faces / (double)ncharts : 1.0;
    if (!(mean_faces > 0.0)) mean_faces = 1.0;
    for (size_t c = 0; c < ncharts; c++) {
        prior_weight[c] = (double)chart_faces[c] / mean_faces;
        if (prior_weight[c] < 0.25) prior_weight[c] = 0.25;
    }
    AtlasRegisterProblem problem;
    memset(&problem, 0, sizeof(problem));
    problem.ncomponents = ncharts;
    problem.edges = redge;
    problem.nedges = nre;
    problem.prior_shift = prior_shift;
    problem.prior_weight = prior_weight;
    AtlasRegisterOptions ropt;
    AtlasRegisterOptions_default(&ropt);
    double *shift = (double *)ARENA_CALLOC(arena, ncharts, sizeof(double));
    AtlasRegisterStats rstats;
    if (AtlasRegister_solve(arena, &problem, &ropt, shift, &rstats) != 0 ||
        rstats.solve_failed)
        return 0;   /* u stays as tabu placed it; the failure is reported */

    stats->relayout_metric_shift_rms = rstats.shift_rms;
    stats->relayout_metric_shift_max = rstats.shift_max;
    memcpy(out_metric_shift, shift, ncharts * sizeof(*out_metric_shift));


    /* Collision direction is part of the discrete answer.  Rank every chart
     * by its post-tabu center before looking at the metric preference; the
     * continuous solve may slide charts, but may not reinterpret their order. */
    AofRelayoutRank *order = (AofRelayoutRank *)ARENA_ALLOC(
        arena, ncharts * sizeof(*order));
    for (size_t c = 0; c < ncharts; c++) {
        order[c].chart = (int32_t)c;
        order[c].u = pu[c];
    }
    qsort(order, ncharts, sizeof(*order), aof_relayout_rank_compare);
    for (size_t r = 0; r < ncharts; r++)
        out_chart_rank[order[r].chart] = r;

    /* An intact discrete edge is topology, not free space.  Exact overlap
     * between its observations is retained for a later ownership cut instead
     * of paying for it by tearing the edge in this metric-only stage. */
    if (ncharts > SIZE_MAX / ncharts) return -1;
    size_t npairs = ncharts * ncharts;
    uint8_t *allowed_pair = (uint8_t *)ARENA_CALLOC(
        arena, npairs, sizeof(*allowed_pair));
    for (size_t e = 0; e < T->nedges; e++) {
        const AofTabuEdge *te = &T->edge[e];
        int dk = (int)win_k[te->a] - (int)win_k[te->b];
        if (dk != te->target) continue;
        allowed_pair[(size_t)te->a * ncharts + (size_t)te->b] = 1;
        allowed_pair[(size_t)te->b * ncharts + (size_t)te->a] = 1;
    }

    double mean_edge = edge_length_count > 0
                     ? edge_length_sum / (double)edge_length_count : 1.0;
    AtlasCollisionRegisterProblem cproblem;
    memset(&cproblem, 0, sizeof(cproblem));
    cproblem.faces = ps->faces;
    cproblem.nfaces = ps->nf;
    cproblem.nvertices = ps->nv;
    cproblem.face_keep = face_keep;
    cproblem.xyz = ps->verts;
    cproblem.base_u = out_u;
    cproblem.v = out_v;
    cproblem.vertex_chart = vertex_chart;
    cproblem.ncharts = ncharts;
    cproblem.allowed_chart_pair = allowed_pair;
    cproblem.chart_radius = T->r_med;
    cproblem.same_sheet_radius_tolerance = 0.6 * fabs(T->dr_dw);
    cproblem.same_sheet_xyz_tolerance = 8.0;
    cproblem.chart_rank = out_chart_rank;
    cproblem.desired_shift = shift;
    cproblem.desired_weight = prior_weight;
    cproblem.edges = redge;
    cproblem.nedges = nre;
    cproblem.topology_deviation_limit = seam_cap;
    cproblem.collision_margin = 1.0e-4 * fmax(1.0, mean_edge);
    cproblem.max_outer_iterations = 64;

    double *collision_shift = (double *)ARENA_CALLOC(
        arena, ncharts, sizeof(*collision_shift));
    AtlasCollisionRegisterStats cstats;
    if (AtlasCollisionRegister_solve(arena, &cproblem, collision_shift,
                                     &cstats) != 0)
        return -1;
    stats->relayout_exact_pairs_before = cstats.exact_pairs_before;
    stats->relayout_exact_pairs_after = cstats.exact_pairs_after;
    stats->relayout_cross_pairs_before = cstats.exact_cross_before;
    stats->relayout_cross_pairs_after = cstats.exact_cross_after;
    stats->relayout_allowed_pairs_before = cstats.exact_allowed_before;
    stats->relayout_allowed_pairs_after = cstats.exact_allowed_after;
    stats->relayout_hard_pairs_before = cstats.exact_hard_before;
    stats->relayout_hard_pairs_after = cstats.exact_hard_after;
    stats->relayout_topology_limited_before =
        cstats.exact_topology_limited_before;
    stats->relayout_topology_limited_after =
        cstats.exact_topology_limited_after;
    stats->relayout_total_bounds = cstats.total_bounds;
    stats->relayout_topology_bounds = cstats.topology_bounds;
    stats->relayout_collision_bounds_rejected =
        cstats.collision_bounds_rejected;
    stats->relayout_collision_bounds = cstats.collision_bounds;
    stats->relayout_collision_bounds_added = cstats.collision_bounds_added;
    stats->relayout_outer_iterations = cstats.outer_iterations;
    stats->relayout_qp_iterations = cstats.qp_iterations;
    stats->relayout_qp_active_bounds = cstats.qp_active_bounds;
    stats->relayout_qp_objective = cstats.qp_objective;
    stats->relayout_qp_stationarity = cstats.qp_stationarity;
    stats->relayout_qp_max_violation = cstats.qp_max_violation;
    stats->relayout_topology_infeasible = cstats.topology_infeasible;
    stats->relayout_collision_failed = cstats.collision_failed;
    if (cstats.solve_failed || cstats.collision_failed) return 0;

    if (report != NULL) {
        AtlasOverlapFixRelayoutBound *bound =
            (AtlasOverlapFixRelayoutBound *)ARENA_ALLOC(
                arena, (cstats.nbound ? cstats.nbound : 1) * sizeof(*bound));
        for (size_t i = 0; i < cstats.nbound; i++) {
            bound[i].chart_lo = cstats.bound[i].chart_lo;
            bound[i].chart_hi = cstats.bound[i].chart_hi;
            bound[i].face_lo = cstats.bound[i].face_lo;
            bound[i].face_hi = cstats.bound[i].face_hi;
            bound[i].lower = cstats.bound[i].lower;
            bound[i].final_delta = cstats.bound[i].final_delta;
            bound[i].slack = cstats.bound[i].slack;
            bound[i].updates = cstats.bound[i].updates;
            bound[i].topology = cstats.bound[i].topology;
            bound[i].collision = cstats.bound[i].collision;
            bound[i].initial = cstats.bound[i].initial;
            bound[i].active = cstats.bound[i].active;
        }
        AtlasOverlapFixRelayoutResidual *residual =
            (AtlasOverlapFixRelayoutResidual *)ARENA_ALLOC(
                arena, (cstats.nresidual ? cstats.nresidual : 1) *
                       sizeof(*residual));
        for (size_t i = 0; i < cstats.nresidual; i++) {
            residual[i].chart0 = cstats.residual[i].chart0;
            residual[i].chart1 = cstats.residual[i].chart1;
            residual[i].face0 = cstats.residual[i].face0;
            residual[i].face1 = cstats.residual[i].face1;
            residual[i].allowed = cstats.residual[i].allowed;
            residual[i].reason = cstats.residual[i].reason;
            residual[i].min_xyz = cstats.residual[i].min_xyz;
            residual[i].radius_delta = cstats.residual[i].radius_delta;
        }
        report->relayout_bound = bound;
        report->nrelayout_bound = cstats.nbound;
        report->relayout_residual = residual;
        report->nrelayout_residual = cstats.nresidual;
    }

    aof_relayout_edge_rms(redge, nre, ATLAS_REGISTER_WELD, prior_shift,
                          &stats->relayout_weld_rms_before);
    aof_relayout_edge_rms(redge, nre, ATLAS_REGISTER_WELD, collision_shift,
                          &stats->relayout_weld_rms_after);
    aof_relayout_edge_rms(redge, nre, ATLAS_REGISTER_RADIAL, prior_shift,
                          &stats->relayout_radial_rms_before);
    aof_relayout_edge_rms(redge, nre, ATLAS_REGISTER_RADIAL, collision_shift,
                          &stats->relayout_radial_rms_after);
    stats->relayout_shift_rms = cstats.shift_rms;
    stats->relayout_shift_max = cstats.shift_max;
    stats->relayout_moved_charts = cstats.moved_charts;
    stats->relayout_atlas_width_before =
        aof_relayout_width(ps, vertex_chart, out_u);
    memcpy(out_chart_shift, collision_shift,
           ncharts * sizeof(*out_chart_shift));
    stats->relayout_failed = 0;
    for (size_t i = 0; i < ps->nv; i++) {
        int32_t c = vertex_chart[i];
        if (c >= 0) out_u[i] += collision_shift[c];
    }
    stats->relayout_atlas_width_after =
        aof_relayout_width(ps, vertex_chart, out_u);
    return 0;
}

static double aof_chart_min_xyz(const PieceSet *ps, const size_t *vrow,
                                 const int32_t *vmember, int32_t ca, int32_t cb)
{
    double best2 = DBL_MAX;
    for (size_t ia = vrow[ca]; ia < vrow[(size_t)ca + 1]; ia++) {
        const float *a = &ps->verts[(size_t)vmember[ia] * 3];
        for (size_t ib = vrow[cb]; ib < vrow[(size_t)cb + 1]; ib++) {
            const float *b = &ps->verts[(size_t)vmember[ib] * 3];
            double dx = (double)a[0] - (double)b[0];
            double dy = (double)a[1] - (double)b[1];
            double dz = (double)a[2] - (double)b[2];
            double d2 = dx * dx + dy * dy + dz * dz;
            if (d2 < best2) best2 = d2;
        }
    }
    return best2 < DBL_MAX ? sqrt(best2) : DBL_MAX;
}

/* Absolute continuous wind per chart.  The source-phase lift carries the
 * wander-free structure within each gauge island; a single rounded-median
 * offset anchors the island to the wind the current field implies
 * (w_phi0 + k).  Charts without phase use that ladder wind directly.
 * Registered u is never inverted through the spiral map. */
static void aof_compute_wind(Arena_T arena, AofTabu *T, size_t ncharts,
                             const int8_t *win_k)
{
    T->wind = (double *)ARENA_ALLOC(arena, ncharts * sizeof(*T->wind));
    double *island_offset = NULL;
    if (T->phase_gauge_islands > 0) {
        island_offset = (double *)ARENA_CALLOC(
            arena, T->phase_gauge_islands, sizeof(*island_offset));
        double *gather = (double *)ARENA_ALLOC(
            arena, ncharts * sizeof(*gather));
        for (size_t i = 0; i < T->phase_gauge_islands; i++) {
            size_t n = 0;
            for (size_t c = 0; c < ncharts; c++)
                if (T->phase_island[c] == (int32_t)i)
                    gather[n++] = T->w_phi0[c] + (double)win_k[c] -
                                  T->w_phase[c];
            if (n == 0) continue;
            qsort(gather, n, sizeof(double), aof_tabu_compare_double);
            island_offset[i] = (double)lround(gather[n / 2]);
        }
    }
    for (size_t c = 0; c < ncharts; c++)
        T->wind[c] = T->phase_island[c] >= 0 && island_offset != NULL
                   ? T->w_phase[c] + island_offset[T->phase_island[c]]
                   : T->w_phi0[c] + (double)win_k[c];
}

/* Tabu placement: search per-chart depths, SAT-verify the leaderboard, apply
 * the true winner to out_u, and record the per-chart decisions. */
static int aof_tabu_place(Arena_T arena, const PieceSet *ps,
                          const ScaffoldCalib *cal,
                          const AtlasOverlapFixOptions *opts,
                          const int32_t *vertex_chart, size_t ncharts,
                          const int32_t *chart_group, size_t ngroups,
                          const AofNeighbour *edge, size_t nedges,
                          const AofChartInfo *chart,
                          const int32_t *face_group,
                          AtlasOverlapFixGroup *group,
                          double *out_u, const double *out_v,
                          const uint8_t *face_keep,
                          AtlasOverlapFixTabuReport *tabu_report,
                          AtlasOverlapFixStats *stats)
{
    AofTabu T;
    size_t excluded = 0;
    if (aof_tabu_prepare(arena, ps, cal, out_u, out_v, face_keep,
                         vertex_chart, ncharts, chart_group, ngroups,
                         edge, nedges, group, opts->wind_correct, opts, 0,
                         &T, &excluded) != 0)
        return -1;
    stats->tabu_edges_excluded = excluded;
    stats->tabu_ladder_step = T.ladder_step;
    stats->tabu_ladder_u0 = T.ladder_u0;
    stats->tabu_ladder_fixed = T.ladder_fixed;
    stats->tabu_ladder_pairs = T.ladder_pairs;
    for (size_t c = 0; c < T.ncharts; c++)
        if (T.active[c]) stats->tabu_active_charts++;
    stats->tabu_phase_gauge_islands = T.phase_gauge_islands;
    stats->tabu_phase_gauge_inconsistent = T.phase_gauge_inconsistent;
    stats->tabu_phase_lateral_trusted = T.phase_lateral_trusted;
    stats->tabu_phase_lateral_nonzero = T.phase_lateral_nonzero;
    stats->tabu_phase_regions = T.n_region;
    stats->tabu_phase_region_links = T.n_region_links;
    for (size_t g = 0; g < ngroups; g++) {
        if (T.graded_group[g]) stats->tabu_graded_groups++;
        if (T.phase_flat_group[g]) stats->tabu_phase_flat_groups++;
    }

    int64_t *cells_initial = (int64_t *)ARENA_ALLOC(
        arena, ncharts * sizeof(int64_t));
    aof_tabu_exact_chart_pairs(&T, cells_initial);

    size_t move_cap = opts->tabu_max_iters > 0
                    ? (size_t)opts->tabu_max_iters : 1;
    size_t archive_cap = move_cap + 1;
    AtlasOverlapFixTabuMove *moves = (AtlasOverlapFixTabuMove *)ARENA_ALLOC(
        arena, move_cap * sizeof(*moves));
    int8_t *archive_k = (int8_t *)ARENA_ALLOC(
        arena, archive_cap * ncharts);
    int64_t *archive_e = (int64_t *)ARENA_ALLOC(
        arena, archive_cap * sizeof(int64_t));
    int64_t *archive_ov = (int64_t *)ARENA_ALLOC(
        arena, archive_cap * sizeof(int64_t));
    int64_t *archive_nb = (int64_t *)ARENA_ALLOC(
        arena, archive_cap * sizeof(int64_t));
    int64_t *archive_rad = (int64_t *)ARENA_ALLOC(
        arena, archive_cap * sizeof(int64_t));
    size_t nmoves = 0, narchive = 0;
    if (aof_tabu_run(arena, &T, opts, moves, move_cap, &nmoves,
                     archive_k, archive_e, archive_ov, archive_nb, archive_rad,
                     archive_cap, &narchive, stats) != 0)
        return -1;

    size_t verify[AOF_TABU_VERIFY_MAX];
    size_t nverify = aof_tabu_choose_verify(
        &T, moves, nmoves, archive_k, archive_e, archive_ov, archive_nb,
        archive_rad, narchive, stats->tabu_best_iter, verify);
    stats->tabu_archive_states = narchive;

    /* Independently re-scan the archived states with the production UV index.
     * The search itself is exact now, so this is both leaderboard attribution
     * (cross versus intra) and an assertion that its separate AABB/sweep broad
     * phase never missed a SAT face pair. */
    size_t nk = (size_t)T.nk;
    double *try_u = (double *)ARENA_ALLOC(arena, ps->nv * sizeof(double));
    size_t *exact_pairs = (size_t *)ARENA_ALLOC(
        arena, nverify * sizeof(size_t));
    size_t *exact_cross = (size_t *)ARENA_ALLOC(
        arena, nverify * sizeof(size_t));
    size_t *exact_intra = (size_t *)ARENA_ALLOC(
        arena, nverify * sizeof(size_t));
    size_t *weld_jumps_25 = (size_t *)ARENA_CALLOC(
        arena, nverify, sizeof(size_t));
    size_t *weld_jumps_100 = (size_t *)ARENA_CALLOC(
        arena, nverify, sizeof(size_t));
    double *weld_mean_abs = (double *)ARENA_CALLOC(
        arena, nverify, sizeof(double));
    double *weld_rms = (double *)ARENA_CALLOC(
        arena, nverify, sizeof(double));
    double *weld_max = (double *)ARENA_CALLOC(
        arena, nverify, sizeof(double));
    for (size_t V = 0; V < nverify; V++) {
        size_t A = verify[V];
        const int8_t *kk = &archive_k[A * ncharts];
        for (size_t i = 0; i < ps->nv; i++) {
            int32_t c = vertex_chart[i];
            try_u[i] = out_u[i];
            if (c >= 0 && kk[c] != 0)
                try_u[i] +=
                    T.shift[(size_t)c * nk + (size_t)(kk[c] + T.span)];
        }
        Arena_Mark vmark = Arena_save(arena);
        AofUvIndex vidx;
        if (aof_uv_index_build(arena, ps, try_u, out_v, face_keep, face_group,
                               -1, opts->overlap_cell_size, &vidx) != 0) {
            Arena_restore(arena, vmark);
            return -1;
        }
        AofCountContext vcount;
        memset(&vcount, 0, sizeof(vcount));
        vcount.faces = ps->faces;
        vcount.face_group = face_group;
        vcount.vertex_chart = vertex_chart;
        vcount.cross_only = 0;
        vcount.after = 1;
        aof_scan_all(ps, face_keep, face_group, &vidx, &vcount);
        exact_cross[V] = vcount.cross_pairs;
        exact_intra[V] = vcount.intra_pairs;
        exact_pairs[V] = exact_cross[V] + exact_intra[V];
        if (archive_ov[A] < 0 || (uint64_t)archive_ov[A] != exact_pairs[V])
            return -1;
        double sum_abs = 0.0, sum_sq = 0.0, max_abs = 0.0;
        for (size_t e = 0; e < T.first_lateral; e++) {
            int32_t ca = T.edge[e].a, cb = T.edge[e].b;
            double sa = kk[ca] != 0
                ? T.shift[(size_t)ca * nk +
                          (size_t)(kk[ca] + T.span)] : 0.0;
            double sb = kk[cb] != 0
                ? T.shift[(size_t)cb * nk +
                          (size_t)(kk[cb] + T.span)] : 0.0;
            double adu = fabs(T.edge[e].du + sa - sb);
            sum_abs += adu;
            sum_sq += adu * adu;
            if (adu > max_abs) max_abs = adu;
            if (adu > 25.0) weld_jumps_25[V]++;
            if (adu > 100.0) weld_jumps_100[V]++;
        }
        if (T.first_lateral > 0) {
            weld_mean_abs[V] = sum_abs / (double)T.first_lateral;
            weld_rms[V] = sqrt(sum_sq / (double)T.first_lateral);
            weld_max[V] = max_abs;
        }
        Arena_restore(arena, vmark);
    }

    int any_zero = 0;
    for (size_t V = 0; V < nverify; V++)
        if (exact_pairs[V] == 0) any_zero = 1;
    size_t winner = (size_t)-1, win_pairs = 0, win_layout = 0;
    int64_t win_topology = 0;
    double win_score = 0.0;
    AtlasOverlapFixTabuCandidate *candidate =
        (AtlasOverlapFixTabuCandidate *)ARENA_CALLOC(
            arena, nverify, sizeof(*candidate));
    for (size_t V = 0; V < nverify; V++) {
        size_t A = verify[V];
        size_t pairs = exact_pairs[V];
        if (any_zero && pairs != 0) continue;
        int64_t topology = archive_nb[A] + archive_rad[A];
        size_t layout = 0;
        const int8_t *kk = &archive_k[A * ncharts];
        for (size_t c = 0; c < ncharts; c++)
            layout += (size_t)(kk[c] < 0 ? -kk[c] : kk[c]);
        double score = (double)T.Lov * (double)pairs +
                       (double)topology;
        candidate[V].state = A;
        candidate[V].exact_cross = exact_cross[V];
        candidate[V].exact_intra = exact_intra[V];
        candidate[V].layout_abs = layout;
        candidate[V].weld_jumps_25 = weld_jumps_25[V];
        candidate[V].weld_jumps_100 = weld_jumps_100[V];
        candidate[V].search_overlap = (double)archive_ov[A];
        candidate[V].neighbour = (double)archive_nb[A] / AOF_TABU_FIX;
        candidate[V].radial = (double)archive_rad[A] / AOF_TABU_FIX;
        candidate[V].energy = (double)archive_e[A] / AOF_TABU_FIX;
        candidate[V].exact_score = score / AOF_TABU_FIX;
        candidate[V].weld_mean_abs = weld_mean_abs[V];
        candidate[V].weld_rms = weld_rms[V];
        candidate[V].weld_max = weld_max[V];
        int better = winner == (size_t)-1;
        if (!better && any_zero) {
            better = topology < win_topology ||
                    (topology == win_topology && layout < win_layout);
        } else if (!better) {
            better = score < win_score ||
                    (score == win_score && pairs < win_pairs) ||
                    (score == win_score && pairs == win_pairs &&
                     topology < win_topology) ||
                    (score == win_score && pairs == win_pairs &&
                     topology == win_topology && layout < win_layout);
        }
        if (better) {
            winner = V;
            win_pairs = pairs;
            win_topology = topology;
            win_layout = layout;
            win_score = score;
        }
    }
    candidate[winner].winner = 1;
    stats->tabu_leaderboard = nverify;
    stats->tabu_winner = winner;
    stats->tabu_winner_iter = verify[winner];
    stats->tabu_exact_score = win_score / AOF_TABU_FIX;

    const int8_t *win_k = &archive_k[verify[winner] * ncharts];
    for (size_t i = 0; i < ps->nv; i++) {
        int32_t c = vertex_chart[i];
        if (c >= 0 && win_k[c] != 0)
            out_u[i] +=
                T.shift[(size_t)c * nk + (size_t)(win_k[c] + T.span)];
    }

    aof_compute_wind(arena, &T, ncharts, win_k);

    /* The pure tabu placement, kept for the relayout revert and for the
     * per-winding-layer diagnostic dumps. */
    double *u_tabu = (double *)ARENA_ALLOC(arena, ps->nv * sizeof(double));
    memcpy(u_tabu, out_u, ps->nv * sizeof(double));

    double *relayout_metric_shift = (double *)ARENA_CALLOC(
        arena, ncharts, sizeof(*relayout_metric_shift));
    double *relayout_shift = (double *)ARENA_CALLOC(
        arena, ncharts, sizeof(*relayout_shift));
    size_t *relayout_rank = (size_t *)ARENA_ALLOC(
        arena, ncharts * sizeof(*relayout_rank));
    for (size_t c = 0; c < ncharts; c++) relayout_rank[c] = c;

    /* Stage C: continuous U re-registration with winding and collision order
     * fixed.  Compatible intersections become exact bounds; incompatible
     * ones remain explicitly classified topology-limited residuals. */
    if (opts->relayout) {
        if (aof_relayout(arena, ps, &T, win_k, vertex_chart, face_keep,
                         out_v, opts->relayout_seam_cap, out_u,
                         relayout_metric_shift, relayout_shift, relayout_rank,
                         tabu_report, stats) != 0)
            return -1;
        if (!stats->relayout_failed) {
            Arena_Mark rmark = Arena_save(arena);
            AofUvIndex ridx;
            if (aof_uv_index_build(arena, ps, out_u, out_v, face_keep,
                                   face_group, -1, opts->overlap_cell_size,
                                   &ridx) != 0) {
                Arena_restore(arena, rmark);
                return -1;
            }
            AofCountContext rcount;
            memset(&rcount, 0, sizeof(rcount));
            rcount.faces = ps->faces;
            rcount.face_group = face_group;
            rcount.vertex_chart = vertex_chart;
            rcount.cross_only = 0;
            rcount.after = 1;
            aof_scan_all(ps, face_keep, face_group, &ridx, &rcount);
            size_t post = rcount.cross_pairs + rcount.intra_pairs;
            Arena_restore(arena, rmark);
            if (post > win_pairs) {
                memcpy(out_u, u_tabu, ps->nv * sizeof(double));
                stats->relayout_reverted = 1;
                memset(relayout_shift, 0,
                       ncharts * sizeof(*relayout_shift));
                if (tabu_report != NULL) {
                    tabu_report->relayout_bound = NULL;
                    tabu_report->nrelayout_bound = 0;
                    tabu_report->relayout_residual = NULL;
                    tabu_report->nrelayout_residual = 0;
                }
            }
        }
    }

    /* Re-derive the final bookkeeping at the winning assignment. */
    memcpy(T.k, win_k, ncharts);
    T.e_ov = aof_tabu_exact_total(&T);
    T.e_nb = aof_tabu_nb_total(&T);
    T.e_rad = aof_tabu_rad_total(&T);
    stats->tabu_energy_final = (double)aof_tabu_energy(&T) / AOF_TABU_FIX;
    int64_t *cells_final = (int64_t *)ARENA_ALLOC(
        arena, ncharts * sizeof(int64_t));
    aof_tabu_exact_chart_pairs(&T, cells_final);

    for (size_t c = 0; c < ncharts; c++) {
        if (win_k[c] == 0) continue;
        stats->tabu_charts_moved++;
        group[chart_group[c]].action = ATLAS_OVERLAP_FIX_SHIFTED;
    }
    stats->tabu_lateral_edges = T.nedges - T.first_lateral;
    for (size_t e = 0; e < T.nedges; e++) {
        if (T.edge[e].family) stats->tabu_family_edges++;
        if (T.k[T.edge[e].a] - T.k[T.edge[e].b] == T.edge[e].target)
            continue;
        if (T.edge[e].family) stats->tabu_family_torn++;
        if (e < T.first_lateral) stats->tabu_welds_torn++;
        else stats->tabu_lateral_torn++;
    }

    if (tabu_report != NULL) {
        AtlasOverlapFixChart *cr = (AtlasOverlapFixChart *)ARENA_ALLOC(
            arena, ncharts * sizeof(*cr));
        for (size_t c = 0; c < ncharts; c++) {
            cr[c].chart = (int32_t)c;
            cr[c].group = chart_group[c];
            cr[c].neighbourhood = T.nbid[c];
            cr[c].phase_flat = T.phase_flat_group[chart_group[c]];
            cr[c].phase_region = T.region_id[chart_group[c]];
            cr[c].region_offset = T.region_id[chart_group[c]] >= 0
                                ? T.region_offset[chart_group[c]] : 0;
            cr[c].cube = T.chart_cube[c];
            cr[c].nfaces = chart[c].nfaces;
            cr[c].nvertices = chart[c].nvertices;
            cr[c].qc_max = chart[c].qc_max;
            cr[c].radius_med = T.r_med[c];
            cr[c].u_med = T.u_med[c];
            cr[c].w_phi0 = T.w_phi0[c];
            cr[c].w_rad = T.w_rad[c];
            cr[c].w_phase = T.w_phase[c];
            cr[c].phase_gauge = T.phase_gauge[c];
            cr[c].phase_bin = T.phase_bin[c];
            cr[c].wind = T.wind[c];
            cr[c].k = win_k[c];
            cr[c].gbin = T.gbin[c];
            cr[c].shift_u = win_k[c] != 0
                ? T.shift[c * nk + (size_t)(win_k[c] + T.span)] : 0.0;
            cr[c].cells_initial = cells_initial[c];
            cr[c].relayout_rank = relayout_rank[c];
            cr[c].relayout_metric_shift_u = relayout_metric_shift[c];
            cr[c].relayout_shift_u = relayout_shift[c];
            cr[c].cells_final = cells_final[c];
        }
        AtlasOverlapFixTabuEdge *er = (AtlasOverlapFixTabuEdge *)ARENA_ALLOC(
            arena, (T.nedges ? T.nedges : 1) * sizeof(*er));
        for (size_t e = 0; e < T.nedges; e++) {
            int32_t ca = T.edge[e].a, cb = T.edge[e].b;
            double da = T.w_phase[ca] - T.w_phi0[ca];
            double db = T.w_phase[cb] - T.w_phi0[cb];
            double q = da - db;
            double qi = floor(q + 0.5);
            er[e].edge = (int32_t)e;
            er[e].chart0 = ca;
            er[e].chart1 = cb;
            er[e].target = T.edge[e].target;
            er[e].final_delta = (int32_t)win_k[ca] - (int32_t)win_k[cb];
            er[e].lateral = e >= T.first_lateral;
            er[e].family = T.edge[e].family;
            er[e].weight = (double)T.edge[e].w / AOF_TABU_FIX;
            er[e].du = T.edge[e].du;
            er[e].final_du = T.edge[e].du + cr[ca].shift_u -
                             cr[cb].shift_u;
            er[e].phase_turn = q;
            er[e].phase_residual = fabs(q - qi);
        }
        tabu_report->chart = cr;
        tabu_report->nchart = ncharts;
        tabu_report->move = moves;
        size_t *vrow = (size_t *)ARENA_CALLOC(
            arena, ncharts + 1, sizeof(*vrow));
        size_t nused = 0;
        for (size_t i = 0; i < ps->nv; i++) {
            int32_t c = vertex_chart[i];
            if (c < 0) continue;
            vrow[(size_t)c + 1]++;
            nused++;
        }
        for (size_t c = 0; c < ncharts; c++) vrow[c + 1] += vrow[c];
        int32_t *vmember = (int32_t *)ARENA_ALLOC(
            arena, (nused ? nused : 1) * sizeof(*vmember));
        size_t *vfill = (size_t *)ARENA_ALLOC(
            arena, ncharts * sizeof(*vfill));
        memcpy(vfill, vrow, ncharts * sizeof(*vfill));
        for (size_t i = 0; i < ps->nv; i++) {
            int32_t c = vertex_chart[i];
            if (c >= 0) vmember[vfill[c]++] = (int32_t)i;
        }

        size_t nop = 0;
        for (size_t ca = 0; ca < ncharts; ca++) {
            for (size_t q = T.prow[ca]; q < T.prow[ca + 1]; q++) {
                int32_t cb = T.pmember[q];
                if ((size_t)cb <= ca) continue;
                if (aof_tabu_exact_pair(&T, (int32_t)ca, win_k[ca],
                                        cb, win_k[cb]) > 0)
                    nop++;
            }
        }
        AtlasOverlapFixTabuOverlap *op =
            (AtlasOverlapFixTabuOverlap *)ARENA_ALLOC(
                arena, (nop ? nop : 1) * sizeof(*op));
        size_t oi = 0;
        for (size_t ca = 0; ca < ncharts; ca++) {
            for (size_t q = T.prow[ca]; q < T.prow[ca + 1]; q++) {
                int32_t cb = T.pmember[q];
                if ((size_t)cb <= ca) continue;
                uint32_t np = aof_tabu_exact_pair(
                    &T, (int32_t)ca, win_k[ca], cb, win_k[cb]);
                if (np == 0) continue;
                AtlasOverlapFixTabuOverlap *o = &op[oi++];
                int32_t e = aof_tabu_find_edge(&T, (int32_t)ca, cb);
                o->chart0 = (int32_t)ca;
                o->chart1 = cb;
                o->face_pairs = np;
                o->edge = e;
                o->target = e >= 0 ? T.edge[e].target : 0;
                o->final_delta = (int32_t)win_k[ca] - (int32_t)win_k[cb];
                o->lateral = e >= 0 ? (e >= (int32_t)T.first_lateral) : -1;
                o->radius_delta = fabs(T.r_med[ca] - T.r_med[cb]);
                o->min_xyz = aof_chart_min_xyz(
                    ps, vrow, vmember, (int32_t)ca, cb);
            }
        }
        if (oi != nop) return -1;
        tabu_report->nmove = nmoves;
        tabu_report->state_k = archive_k;
        tabu_report->nstate = narchive;
        tabu_report->candidate = candidate;
        tabu_report->ncandidate = nverify;
        tabu_report->edge = er;
        tabu_report->nedge = T.nedges;
        tabu_report->overlap = op;
        tabu_report->noverlap = nop;
        tabu_report->vertex_chart = vertex_chart;
        tabu_report->u_post_tabu = u_tabu;
    }
    return 0;
}

/* ========================================================================== */
/* Solve                                                                       */
/* ========================================================================== */

int AtlasOverlapFix_solve(Arena_T arena,
                          const PieceSet *ps,
                          const ScaffoldCalib *cal,
                          const float *registered_u,
                          const AtlasOverlapFixOptions *input_opts,
                          double *out_u,
                          double *out_v,
                          uint8_t *face_keep,
                          int32_t *out_vertex_group,
                          AtlasOverlapFixGroup **out_group,
                          size_t *out_ngroups,
                          AtlasOverlapFixTabuReport *tabu_report,
                          AtlasOverlapFixStats *stats)
{
    g_aof_t0 = ves_clock_sec();
    if (arena == NULL || ps == NULL || cal == NULL || registered_u == NULL ||
        out_u == NULL || out_v == NULL || face_keep == NULL || stats == NULL ||
        ps->nv == 0 || ps->nf == 0)
        return -1;

    AtlasOverlapFixOptions defaults;
    AtlasOverlapFixOptions_default(&defaults);
    const AtlasOverlapFixOptions *opts =
        input_opts != NULL ? input_opts : &defaults;
    if (!(opts->aspect_ratio_kill > 1.0) || opts->max_shift_wraps < 0)
        return -1;
    if (opts->chart_qc_kill != 0.0 && !(opts->chart_qc_kill > 1.0))
        return -1;

    memset(stats, 0, sizeof(*stats));
    stats->total_faces = ps->nf;
    if (out_group != NULL) *out_group = NULL;
    if (out_ngroups != NULL) *out_ngroups = 0;
    if (tabu_report != NULL) memset(tabu_report, 0, sizeof(*tabu_report));

    for (size_t i = 0; i < ps->nv; i++) {
        out_u[i] = (double)registered_u[i];
        out_v[i] = aof_axis_coord(cal, &ps->verts[i * 3]);
    }
    for (size_t f = 0; f < ps->nf; f++) {
        int ok = 1;
        for (int k = 0; k < 3; k++) {
            int32_t vertex = ps->faces[f * 3 + (size_t)k];
            if (vertex < 0 || (size_t)vertex >= ps->nv) ok = 0;
        }
        face_keep[f] = (uint8_t)ok;
    }

    /* 1. Kill hopeless triangles before anything measures anything. */
    stats->faces_killed_stretch = aof_kill_bad_triangles(
        ps, out_u, out_v, opts->aspect_ratio_kill, face_keep);
    stats->faces_killed = stats->faces_killed_stretch;

    /* 2. Label once to measure/cull floaters, then relabel the survivors. */
    int32_t *vertex_chart = (int32_t *)ARENA_ALLOC(
        arena, ps->nv * sizeof(int32_t));
    aof_phase("label charts");
    size_t ncharts = aof_label_charts(arena, ps, face_keep, vertex_chart);
    stats->charts_before_cull = ncharts;
    if (ncharts == 0) return 0;

    AofChartInfo *chart = (AofChartInfo *)ARENA_ALLOC(
        arena, ncharts * sizeof(*chart));
    aof_chart_info(ps, cal, vertex_chart, ncharts, out_u, out_v,
                   face_keep, chart);
    stats->faces_killed_charts = aof_kill_bad_charts(
        arena, ps, vertex_chart, ncharts, chart, opts, face_keep, stats);
    stats->faces_killed += stats->faces_killed_charts;
    if (stats->faces_killed_charts > 0) {
        ncharts = aof_label_charts(arena, ps, face_keep, vertex_chart);
        if (ncharts == 0) return 0;
        chart = (AofChartInfo *)ARENA_ALLOC(arena,
                                            ncharts * sizeof(*chart));
        aof_chart_info(ps, cal, vertex_chart, ncharts, out_u, out_v,
                       face_keep, chart);
    }
    stats->total_charts = ncharts;

    /* 3. Authentic neighbours, then groups. */
    AofNeighbour *edge = NULL;
    size_t nedges = 0;
    aof_phase("build neighbours");
    if (aof_build_neighbours(arena, ps, cal, vertex_chart, out_u, opts,
                             &edge, &nedges) != 0)
        return -1;
    stats->authentic_neighbour_edges = nedges;

    /* How close each weld's measured mismatch is to a whole number of turns.
     * This is the falsifiable part of the premise: if the errors really are
     * per-cube integer wind choices, this histogram piles up at 0. */
    for (size_t e = 0; e < nedges; e++) {
        double r = edge[e].residual;
        size_t bin = r < 0.05 ? 0 : r < 0.15 ? 1 : r < 0.25 ? 2
                   : r < 0.35 ? 3 : 4;
        stats->weld_residual_hist[bin]++;
        if (edge[e].wraps != 0) stats->weld_edges_nonzero_wrap++;
        r = edge[e].phase_residual;
        bin = r < 0.05 ? 0 : r < 0.15 ? 1 : r < 0.25 ? 2
            : r < 0.35 ? 3 : 4;
        stats->weld_phase_residual_hist[bin]++;
        if (edge[e].phase_target != 0) stats->weld_phase_edges_nonzero++;
    }

    int32_t *chart_group = (int32_t *)ARENA_ALLOC(
        arena, ncharts * sizeof(int32_t));
    aof_phase("label groups");
    size_t ngroups = aof_label_groups(arena, ncharts, edge, nedges,
                                      chart_group);

    size_t group_faces_killed = aof_kill_small_groups(
        arena, ps, vertex_chart, chart_group, ngroups,
        opts->floater_group_max_faces, face_keep, stats);
    if (group_faces_killed > 0) {
        stats->faces_killed += group_faces_killed;
        ncharts = aof_label_charts(arena, ps, face_keep, vertex_chart);
        if (ncharts == 0) return 0;
        chart = (AofChartInfo *)ARENA_ALLOC(arena,
                                            ncharts * sizeof(*chart));
        aof_chart_info(ps, cal, vertex_chart, ncharts, out_u, out_v,
                       face_keep, chart);
        edge = NULL;
        nedges = 0;
        if (aof_build_neighbours(arena, ps, cal, vertex_chart, out_u, opts,
                                 &edge, &nedges) != 0)
            return -1;
        chart_group = (int32_t *)ARENA_ALLOC(
            arena, ncharts * sizeof(int32_t));
        ngroups = aof_label_groups(arena, ncharts, edge, nedges, chart_group);

        /* The first histogram included edges belonging to the culled toys;
         * replace it with evidence from exactly the mesh entering the solve. */
        memset(stats->weld_residual_hist, 0,
               sizeof(stats->weld_residual_hist));
        memset(stats->weld_phase_residual_hist, 0,
               sizeof(stats->weld_phase_residual_hist));
        stats->weld_edges_nonzero_wrap = 0;
        stats->weld_phase_edges_nonzero = 0;
        for (size_t e = 0; e < nedges; e++) {
            double residual = edge[e].residual;
            size_t bin = residual < 0.05 ? 0 : residual < 0.15 ? 1
                       : residual < 0.25 ? 2 : residual < 0.35 ? 3 : 4;
            stats->weld_residual_hist[bin]++;
            if (edge[e].wraps != 0) stats->weld_edges_nonzero_wrap++;
            residual = edge[e].phase_residual;
            bin = residual < 0.05 ? 0 : residual < 0.15 ? 1
                : residual < 0.25 ? 2 : residual < 0.35 ? 3 : 4;
            stats->weld_phase_residual_hist[bin]++;
            if (edge[e].phase_target != 0)
                stats->weld_phase_edges_nonzero++;
        }
    }
    stats->total_charts = ncharts;
    stats->authentic_neighbour_edges = nedges;
    stats->total_groups = ngroups;

    int32_t *face_group = (int32_t *)ARENA_ALLOC(
        arena, ps->nf * sizeof(int32_t));
    for (size_t f = 0; f < ps->nf; f++) {
        face_group[f] = -1;
        if (!face_keep[f]) continue;
        int32_t c = vertex_chart[ps->faces[f * 3]];
        if (c >= 0) face_group[f] = chart_group[c];
    }
    if (out_vertex_group != NULL)
        for (size_t i = 0; i < ps->nv; i++)
            out_vertex_group[i] = vertex_chart[i] >= 0
                                ? chart_group[vertex_chart[i]] : -1;

    AtlasOverlapFixGroup *group = (AtlasOverlapFixGroup *)ARENA_CALLOC(
        arena, ngroups, sizeof(*group));
    for (size_t g = 0; g < ngroups; g++) {
        group[g].group = (int32_t)g;
        group[g].action = ATLAS_OVERLAP_FIX_SETTLED;
    }

    /* 3b. Overlap exactly as it arrived, so the wind correction below can be
     *     judged against it rather than asserted.  A global wind-only seed
     *     deliberately skips this expensive exact scan: it performs no atlas
     *     optimization and exists only to establish one coherent gauge before
     *     bounded core/halo subproblems begin. */
    if (!opts->wind_only) {
        Arena_Mark raw_mark = Arena_save(arena);
        AofUvIndex raw;
        aof_phase("uv index: raw");
        aof_phase("uv index: before");
        if (aof_uv_index_build(arena, ps, out_u, out_v, face_keep, face_group,
                               -1, opts->overlap_cell_size, &raw) != 0) {
            Arena_restore(arena, raw_mark);
            return -1;
        }
        AofCountContext raw_count;
        memset(&raw_count, 0, sizeof(raw_count));
        raw_count.faces = ps->faces;
        raw_count.face_group = face_group;
        raw_count.vertex_chart = vertex_chart;
        raw_count.group = group;
        aof_phase("scan: raw");
        aof_scan_all(ps, face_keep, face_group, &raw, &raw_count);
        stats->overlap_pairs_raw = raw_count.cross_pairs;
        stats->intra_overlap_pairs_raw = raw_count.intra_pairs;
        Arena_restore(arena, raw_mark);
        for (size_t g = 0; g < ngroups; g++) {
            group[g].overlap_raw = group[g].overlap_before;
            group[g].overlap_before = 0;
            group[g].intra_overlap = 0;
            group[g].intra_wrongwind = 0;
        }
    }

    /* 3c. Reconcile the per-cube integer wind gauges against the welds.  This
     * is the correction the incoming field never received: its cross-cube
     * winding registration reports zero accepted moves, so every cube's turn
     * index was rounded on its own and welded neighbours can sit a whole
     * circumference apart in u while touching in XYZ. */
    double *chart_shift = (double *)ARENA_ALLOC(
        arena, ncharts * sizeof(*chart_shift));
    AofWindStats wind;
    if (opts->wind_correct) {
        aof_phase("integrate winds");
        if (aof_integrate_winds(arena, ncharts, edge, nedges,
                                opts->wind_residual_limit, opts->wind_cap,
                                chart_shift, &wind) != 0)
            return -1;
        for (size_t i = 0; i < ps->nv; i++)
            if (vertex_chart[i] >= 0)
                out_u[i] -= chart_shift[vertex_chart[i]];
    } else {
        memset(&wind, 0, sizeof(wind));
        for (size_t c = 0; c < ncharts; c++) chart_shift[c] = 0.0;
    }
    stats->weld_edges_trusted = wind.trusted_edges;
    stats->weld_edges_untrusted = wind.untrusted_edges;
    stats->weld_edges_oversized = wind.oversized_edges;
    stats->gauge_tree_edges = wind.tree_edges;
    stats->gauge_cycle_edges = wind.cycle_edges;
    stats->gauge_cycle_consistent = wind.cycle_consistent;
    stats->gauge_cycle_inconsistent = wind.cycle_inconsistent;
    stats->gauge_islands = wind.gauge_islands;
    stats->charts_wind_corrected = wind.charts_corrected;
    stats->wind_wraps_min = wind.wraps_min;
    stats->wind_wraps_max = wind.wraps_max;
    if (opts->wind_only) {
        if (out_group != NULL) *out_group = group;
        if (out_ngroups != NULL) *out_ngroups = ngroups;
        return 0;
    }
    double *radius_sum = (double *)ARENA_CALLOC(arena, ngroups,
                                                sizeof(*radius_sum));
    double *radius_lo = (double *)ARENA_ALLOC(
        arena, ngroups * sizeof(*radius_lo));
    double *radius_hi = (double *)ARENA_ALLOC(
        arena, ngroups * sizeof(*radius_hi));
    for (size_t g = 0; g < ngroups; g++) {
        radius_lo[g] = DBL_MAX;
        radius_hi[g] = -DBL_MAX;
    }
    for (size_t c = 0; c < ncharts; c++) {
        int32_t g = chart_group[c];
        group[g].ncharts++;
        group[g].nfaces += chart[c].nfaces;
        group[g].nvertices += chart[c].nvertices;
        if (group[g].ncharts > stats->largest_group_charts)
            stats->largest_group_charts = group[g].ncharts;
    }
    for (size_t i = 0; i < ps->nv; i++) {
        if (vertex_chart[i] < 0) continue;
        int32_t g = chart_group[vertex_chart[i]];
        double r = aof_radius(cal, &ps->verts[i * 3]);
        radius_sum[g] += r;
        if (r < radius_lo[g]) radius_lo[g] = r;
        if (r > radius_hi[g]) radius_hi[g] = r;
    }
    for (size_t g = 0; g < ngroups; g++) {
        if (group[g].nvertices > 0)
            group[g].mean_radius = radius_sum[g] / (double)group[g].nvertices;
        group[g].radius_spread = radius_hi[g] > radius_lo[g]
                               ? radius_hi[g] - radius_lo[g] : 0.0;
        group[g].circumference = AOF_2PI * group[g].mean_radius;
    }

    /* Mean radius per face, so overlap pairs can be binned by radial
     * separation.  Killed faces keep 0; they are never queried. */
    double *face_radius = (double *)ARENA_CALLOC(arena, ps->nf,
                                                 sizeof(*face_radius));
    for (size_t f = 0; f < ps->nf; f++) {
        if (!face_keep[f]) continue;
        double r = 0.0;
        for (int k = 0; k < 3; k++)
            r += aof_radius(
                cal, &ps->verts[(size_t)ps->faces[f * 3 + (size_t)k] * 3]);
        face_radius[f] = r / 3.0;
    }

    /* 4. Baseline overlap.  Cross-group pairs are what a shift can fix; pairs
     *    inside one group move together and are reported, not chased. */
    Arena_Mark mark = Arena_save(arena);
    AofUvIndex all;
    if (aof_uv_index_build(arena, ps, out_u, out_v, face_keep, face_group,
                           -1, opts->overlap_cell_size, &all) != 0) {
        Arena_restore(arena, mark);
        return -1;
    }
    AofCountContext before_count;
    memset(&before_count, 0, sizeof(before_count));
    before_count.faces = ps->faces;
    before_count.face_group = face_group;
    before_count.vertex_chart = vertex_chart;
    before_count.group = group;
    before_count.cross_only = 0;
    before_count.face_radius = face_radius;
    before_count.pitch = cal->pitch;
    before_count.radial_hist = stats->intra_radial_hist;
    aof_phase("scan: before");
    aof_scan_all(ps, face_keep, face_group, &all, &before_count);
    stats->overlap_pairs_before = before_count.cross_pairs;
    stats->intra_overlap_pairs = before_count.intra_pairs;
    Arena_restore(arena, mark);
    for (size_t g = 0; g < ngroups; g++)
        if (group[g].overlap_before > 0) stats->groups_overlapping_before++;

    /* Park slots start below everything the atlas currently occupies. */
    double atlas_v_min = DBL_MAX;
    for (size_t i = 0; i < ps->nv; i++)
        if (vertex_chart[i] >= 0 && out_v[i] < atlas_v_min)
            atlas_v_min = out_v[i];
    if (atlas_v_min == DBL_MAX) atlas_v_min = 0.0;
    stats->park_v_base = atlas_v_min;
    double park_cursor = atlas_v_min;

    /* 5. Place: per-chart tabu re-gauge on request, else the rigid greedy
     *    group placement (furthest from the axis first) with parking. */
    if (opts->tabu) {
        aof_phase("tabu place");
        if (aof_tabu_place(arena, ps, cal, opts, vertex_chart, ncharts,
                           chart_group, ngroups, edge, nedges, chart,
                           face_group, group, out_u, out_v, face_keep,
                           tabu_report, stats) != 0)
            return -1;
    } else {
    size_t nwork = 0;
    for (size_t g = 0; g < ngroups; g++)
        if (group[g].overlap_before > 0) nwork++;
    AofOrderEntry *order = (AofOrderEntry *)ARENA_ALLOC(
        arena, (nwork ? nwork : 1) * sizeof(*order));
    size_t at_order = 0;
    for (size_t g = 0; g < ngroups; g++) {
        if (group[g].overlap_before == 0) continue;
        order[at_order].group = (int32_t)g;
        order[at_order].radius = group[g].mean_radius;
        at_order++;
    }
    qsort(order, nwork, sizeof(*order), aof_compare_order);

    for (size_t w = 0; w < nwork; w++) {
        int32_t g = order[w].group;
        double circumference = group[g].circumference;

        mark = Arena_save(arena);
        AofUvIndex rest;
        if (aof_uv_index_build(arena, ps, out_u, out_v, face_keep, face_group,
                               g, opts->overlap_cell_size, &rest) != 0) {
            Arena_restore(arena, mark);
            return -1;
        }

        /* Still overlapping after earlier groups moved? */
        size_t current = 0;
        for (size_t f = 0; f < ps->nf && current == 0; f++) {
            if (!face_keep[f] || face_group[f] != g) continue;
            current = aof_query_face(&rest, f, 0.0, 0.0, 0, 1, NULL, NULL);
        }
        if (current == 0) {
            group[g].action = ATLAS_OVERLAP_FIX_SETTLED;
            Arena_restore(arena, mark);
            continue;
        }

        /* Try -1, +1, -2, +2, ... wraps and take the first clean placement. */
        int placed = 0;
        if (circumference > 0.0) {
            for (int magnitude = 1;
                 magnitude <= opts->max_shift_wraps && !placed; magnitude++) {
                for (int sign = -1; sign <= 1 && !placed; sign += 2) {
                    double shift = (double)(sign * magnitude) * circumference;
                    size_t hits = 0;
                    for (size_t f = 0; f < ps->nf && hits == 0; f++) {
                        if (!face_keep[f] || face_group[f] != g) continue;
                        hits = aof_query_face(&rest, f, shift, 0.0, 0, 1,
                                              NULL, NULL);
                    }
                    if (hits > 0) continue;
                    for (size_t i = 0; i < ps->nv; i++)
                        if (vertex_chart[i] >= 0 &&
                            chart_group[vertex_chart[i]] == g)
                            out_u[i] += shift;
                    group[g].shift_u = shift;
                    group[g].wraps = sign * magnitude;
                    group[g].action = ATLAS_OVERLAP_FIX_SHIFTED;
                    placed = 1;
                }
            }
        }
        Arena_restore(arena, mark);
        if (placed) continue;

        if (!opts->park_unplaceable) {
            group[g].action = ATLAS_OVERLAP_FIX_STUCK;
            continue;
        }

        /* Park the whole group into a private v band under the atlas.  Each
         * band is disjoint from every other, so parked groups can never
         * overlap the atlas or each other regardless of their u. */
        double v_lo = DBL_MAX, v_hi = -DBL_MAX;
        for (size_t i = 0; i < ps->nv; i++) {
            if (vertex_chart[i] < 0 || chart_group[vertex_chart[i]] != g)
                continue;
            if (out_v[i] < v_lo) v_lo = out_v[i];
            if (out_v[i] > v_hi) v_hi = out_v[i];
        }
        if (v_lo == DBL_MAX) { v_lo = 0.0; v_hi = 0.0; }
        double shift_v = (park_cursor - opts->park_margin) - v_hi;
        for (size_t i = 0; i < ps->nv; i++)
            if (vertex_chart[i] >= 0 && chart_group[vertex_chart[i]] == g)
                out_v[i] += shift_v;
        park_cursor = park_cursor - opts->park_margin - (v_hi - v_lo);
        group[g].shift_v = shift_v;
        group[g].action = ATLAS_OVERLAP_FIX_PARKED;
    }
    }

    /* 6. Re-measure. */
    mark = Arena_save(arena);
    AofUvIndex after;
    if (aof_uv_index_build(arena, ps, out_u, out_v, face_keep, face_group,
                           -1, opts->overlap_cell_size, &after) != 0) {
        Arena_restore(arena, mark);
        return -1;
    }
    AofCountContext after_count;
    memset(&after_count, 0, sizeof(after_count));
    after_count.faces = ps->faces;
    after_count.face_group = face_group;
    after_count.vertex_chart = vertex_chart;
    after_count.group = group;
    after_count.after = 1;
    /* The greedy cannot change intra-group overlap, so it skips counting it.
     * The tabu re-gauge exists to change it, so it must be re-measured. */
    after_count.cross_only = opts->tabu ? 0 : 1;
    if (opts->tabu) {
        after_count.face_radius = face_radius;
        after_count.pitch = cal->pitch;
        after_count.radial_hist = stats->intra_radial_hist_after;
    }
    aof_phase("scan: after");
    aof_scan_all(ps, face_keep, face_group, &after, &after_count);
    stats->overlap_pairs_after = after_count.cross_pairs;
    stats->intra_overlap_pairs_after = after_count.intra_pairs;
    Arena_restore(arena, mark);

    for (size_t g = 0; g < ngroups; g++) {
        if (group[g].overlap_after > 0) stats->groups_overlapping_after++;
        switch (group[g].action) {
        case ATLAS_OVERLAP_FIX_SHIFTED: stats->shifted_groups++; break;
        case ATLAS_OVERLAP_FIX_PARKED:  stats->parked_groups++;  break;
        case ATLAS_OVERLAP_FIX_STUCK:   stats->stuck_groups++;   break;
        default:                        stats->settled_groups++; break;
        }
    }
    if (out_group != NULL) *out_group = group;
    if (out_ngroups != NULL) *out_ngroups = ngroups;
    return 0;
}

/* ========================================================================== */
/* Self-test                                                                   */
/* ========================================================================== */

static int aof_triangle_cull_selftest(void)
{
    /* Face 0 is genuinely skinny in both XYZ and UV and must survive.  Face
     * 1 is regular in XYZ but has one UV edge stretched 100x and must die. */
    float verts[18] = {
        0, 0, 0,  100, 0, 0,  0, 1, 0,
        0, 0, 0,    1, 0, 0,  0, 1, 0
    };
    int32_t faces[6] = {0, 1, 2, 3, 4, 5};
    double u[6] = {0, 100, 0, 0, 100, 0};
    double v[6] = {0, 0, 1, 0, 0, 1};
    uint8_t keep[2] = {1, 1};
    PieceSet ps;
    memset(&ps, 0, sizeof(ps));
    ps.verts = verts;
    ps.faces = faces;
    ps.nv = 6;
    ps.nf = 2;
    size_t killed = aof_kill_bad_triangles(&ps, u, v, 8.0, keep);
    double qc0 = aof_face_qc(&ps, u, v, 0, 1, 2);
    double qc1 = aof_face_qc(&ps, u, v, 3, 4, 5);
    if (fabs(qc0 - 1.0) > 1.0e-9 || qc1 < 50.0) {
        fprintf(stderr,
            "[overlap_fix selftest] FAIL: Sander QC source-sliver=%.9g "
            "projected-sliver=%.9g\n", qc0, qc1);
        return -1;
    }

    if (killed != 1 || keep[0] != 1 || keep[1] != 0) {
        fprintf(stderr,
            "[overlap_fix selftest] FAIL: UV/XYZ stretch cull killed=%zu "
            "keep=%u,%u\n", killed, (unsigned)keep[0], (unsigned)keep[1]);
        return -1;
    }
    return 0;
}

/* ---- tabu fixtures ------------------------------------------------------ */
/* All tabu tests run on the calibrated spiral a=100, b=-9.5, sense=-1 (the
 * 0139 sign convention, deliberately, to lock the sense arithmetic): wind w
 * lives at radius 100 + 9.5w, and u DECREASES outward.  Charts are small
 * triangles near wind 2 (r=119, u ~= -1376). */

static void aof_tabu_test_calib(ScaffoldCalib *cal)
{
    memset(cal, 0, sizeof(*cal));
    cal->spiral_a = 100.0;
    cal->spiral_b = -9.5;
    cal->pitch = 9.5;
    cal->sense = -1;
    cal->axis_dir[0] = 1.0f;   /* axis along z, so v == z, radius from (y,x) */
}

/* Fill the boilerplate PieceSet fields for n_charts 3-vertex charts, one per
 * cube.  verts/registered_u/uv/phi/normals/faces etc. are caller storage. */
static void aof_tabu_test_pieceset(PieceSet *ps, float *verts, float *uv,
                                   float *phi, float *normals, int32_t *faces,
                                   int32_t *face_cube, int32_t *gid,
                                   size_t *cube_voff, long (*cube_org)[3],
                                   char (*ids)[48], size_t n_charts)
{
    for (size_t c = 0; c < n_charts; c++) {
        for (int k = 0; k < 3; k++)
            faces[c * 3 + (size_t)k] = (int32_t)(c * 3 + (size_t)k);
        face_cube[c] = (int32_t)c;
        cube_voff[c] = c * 3;
        cube_org[c][0] = 0;
        cube_org[c][1] = 0;
        cube_org[c][2] = (long)c * 128;
        snprintf(ids[c], 48, "z00000_y00000_x%05zu", c * 128);
    }
    cube_voff[n_charts] = n_charts * 3;
    memset(normals, 0, n_charts * 9 * sizeof(float));
    for (size_t i = 0; i < n_charts * 3; i++) {
        gid[i] = -1;
        uv[i * 2] = 0.0f;
        uv[i * 2 + 1] = 0.0f;
    }
    memset(ps, 0, sizeof(*ps));
    ps->verts = verts;
    ps->uv = uv;
    ps->phi = phi;
    ps->normals = normals;
    ps->faces = faces;
    ps->face_cube = face_cube;
    ps->gid = gid;
    ps->nv = n_charts * 3;
    ps->nf = n_charts;
    ps->ids = ids;
    ps->cube_voff = cube_voff;
    ps->cube_org = cube_org;
    ps->n_cubes = n_charts;
}

/* A proximity weld exists only across cube provenance.  Two disconnected
 * components from one cube can be adjacent scroll layers: close, parallel,
 * and even aligned in u, but never two halves of a cube seam. */
static int aof_same_cube_weld_selftest(void)
{
    Arena_T arena = Arena_new();
    if (arena == NULL) return -1;

    float verts[] = {
        0.0f, 119.0f, 0.0f,  6.0f, 119.0f, 0.0f,  3.0f, 119.0f, 2.0f,
        0.0f, 119.0f, 4.0f,  6.0f, 119.0f, 4.0f,  3.0f, 119.0f, 6.0f
    };
    float normals[18] = {0};
    float phi[6] = {0};
    int32_t faces[] = {0, 1, 2, 3, 4, 5};
    double u[] = {0.0, 6.0, 3.0, 4.0, 10.0, 7.0};
    uint8_t face_keep[] = {1, 1};
    size_t one_cube[] = {0, 6};
    size_t two_cubes[] = {0, 3, 6};

    PieceSet ps;
    memset(&ps, 0, sizeof(ps));
    ps.verts = verts;
    ps.normals = normals;
    ps.phi = phi;
    ps.faces = faces;
    ps.nv = 6;
    ps.nf = 2;
    ps.cube_voff = one_cube;
    ps.n_cubes = 1;

    ScaffoldCalib cal;
    aof_tabu_test_calib(&cal);
    AtlasOverlapFixOptions opts;
    AtlasOverlapFixOptions_default(&opts);
    opts.neighbour_min_shared = 3;

    int32_t vertex_chart[6];
    size_t ncharts = aof_label_charts(
        arena, &ps, face_keep, vertex_chart);
    AofNeighbour *edge = NULL;
    size_t nedges = 0;
    int failed = 0;
    if (ncharts != 2 ||
        aof_build_neighbours(arena, &ps, &cal, vertex_chart, u, &opts,
                             &edge, &nedges) != 0 || nedges != 0) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: same-cube surfaces "
                        "became %zu weld edges\n", nedges);
        failed = 1;
    }

    ps.cube_voff = two_cubes;
    ps.n_cubes = 2;
    edge = NULL;
    nedges = 0;
    if (!failed &&
        (aof_build_neighbours(arena, &ps, &cal, vertex_chart, u, &opts,
                              &edge, &nedges) != 0 || nedges != 1)) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: cross-cube seam has "
                        "%zu weld edges\n", nedges);
        failed = 1;
    }
    if (!failed)
        fprintf(stderr, "[overlap_fix selftest] same-cube weld rejection ok\n");
    Arena_dispose(&arena);
    return failed ? -1 : 0;
}

/*
 * T1 -- a displaced chart comes home through the welds.
 *
 * Charts A, B, C sit side by side on wind 2 and weld in a chain A-B-C (>= 3
 * close vertex pairs per seam).  B's registered u was mis-gauged one wind
 * OUT, which no rigid group move can see, and the wind integration is OFF so
 * only the tabu search can repair it: the A-B and B-C welds measure wraps of
 * one whole turn, the neighbour term prices the tear, the radial prior
 * agrees.  Expect k = {0,-1,0}, no torn welds, and a bit-identical rerun.
 */
static int aof_tabu_return_home_selftest(void)
{
    Arena_T arena = Arena_new();
    if (arena == NULL) return -1;

    double pi = AOF_2PI * 0.5;
    double phi2 = -2.0 * AOF_2PI;                    /* wind 2: phi = -4pi */
    double u_w2 = 100.0 * phi2 + (-9.5) * phi2 * phi2 / (4.0 * pi);
    double du_out = CubeReg_deltaU(100.0, -9.5, phi2, -1);  /* one wind OUT */

    float verts[] = {
        0.0f, 119.0f, 0.0f,   6.0f, 119.0f, 0.0f,   3.0f, 119.0f, 2.0f,
        0.0f, 119.0f, 4.0f,   6.0f, 119.0f, 4.0f,   3.0f, 119.0f, 6.0f,
        0.0f, 119.0f, 8.0f,   6.0f, 119.0f, 8.0f,   3.0f, 119.0f, 10.0f
    };
    double x_off[] = {0, 0, 2,  4, 4, 6,  8, 8, 10};
    float registered_u[9], phi[9], uv[18], normals[27];
    int32_t faces[9], face_cube[3], gid[9];
    size_t cube_voff[4];
    long cube_org[3][3];
    char ids[3][48];
    for (int i = 0; i < 9; i++) {
        double base = u_w2 + x_off[i];
        if (i >= 3 && i < 6) base += du_out;         /* B mis-gauged */
        registered_u[i] = (float)base;
        phi[i] = (float)phi2;
    }
    PieceSet ps;
    aof_tabu_test_pieceset(&ps, verts, uv, phi, normals, faces, face_cube,
                           gid, cube_voff, cube_org, ids, 3);
    ScaffoldCalib cal;
    aof_tabu_test_calib(&cal);

    AtlasOverlapFixOptions opts;
    AtlasOverlapFixOptions_default(&opts);
    opts.tabu = 1;
    opts.wind_correct = 0;       /* isolate the tabu repair */
    opts.tabu_span = 2;
    opts.tabu_tenure = 5;
    opts.tabu_max_iters = 60;
    opts.tabu_stall = 25;

    int failed = 0;
    double energy_first = 0.0;
    size_t nmove_first = 0;
    int32_t k_first[3] = {0, 0, 0};
    for (int run = 0; run < 2 && !failed; run++) {
        double out_u[9], out_v[9];
        uint8_t face_keep[3];
        AtlasOverlapFixTabuReport report;
        AtlasOverlapFixStats stats;
        int rc = AtlasOverlapFix_solve(arena, &ps, &cal, registered_u, &opts,
                                       out_u, out_v, face_keep, NULL,
                                       NULL, NULL, &report, &stats);
        if (rc != 0 || report.nchart != 3) {
            fprintf(stderr, "[overlap_fix selftest] FAIL: return-home rc=%d "
                            "ncharts=%zu\n", rc, report.nchart);
            failed = 1;
            break;
        }
        if (report.chart[0].k != 0 || report.chart[1].k != -1 ||
            report.chart[2].k != 0) {
            fprintf(stderr, "[overlap_fix selftest] FAIL: return-home "
                            "k={%d,%d,%d}, expected {0,-1,0}\n",
                    report.chart[0].k, report.chart[1].k, report.chart[2].k);
            failed = 1;
            break;
        }
        if (stats.tabu_welds_torn != 0 || stats.tabu_charts_moved != 1) {
            fprintf(stderr, "[overlap_fix selftest] FAIL: return-home "
                            "torn=%zu moved=%zu\n",
                    stats.tabu_welds_torn, stats.tabu_charts_moved);
            failed = 1;
            break;
        }
        if (run == 0) {
            energy_first = stats.tabu_energy_final;
            nmove_first = report.nmove;
            for (int c = 0; c < 3; c++) k_first[c] = report.chart[c].k;
        } else if (stats.tabu_energy_final != energy_first ||
                   report.nmove != nmove_first ||
                   report.chart[0].k != k_first[0] ||
                   report.chart[1].k != k_first[1] ||
                   report.chart[2].k != k_first[2]) {
            fprintf(stderr, "[overlap_fix selftest] FAIL: return-home rerun "
                            "diverged (determinism)\n");
            failed = 1;
        }
    }
    if (!failed)
        fprintf(stderr, "[overlap_fix selftest] tabu return-home ok\n");
    Arena_dispose(&arena);
    return failed ? -1 : 0;
}

/*
 * The winding handoff: absolute wind must survive convergence.
 *
 * Round 1 repairs B (one wind out) exactly like T1, with an extra lone
 * chart D correctly gauged on wind 1.  Round 2 re-solves from round 1's
 * output the way the tool's --rounds loop does (wind_correct off, u fed
 * back); it must converge with every k == 0 -- the exact state in which
 * the old handoff serialized all-zero windings -- while the reported
 * absolute wind still says A,B,C sit one full turn out from D.
 */
static int aof_rounds_handoff_selftest(void)
{
    Arena_T arena = Arena_new();
    if (arena == NULL) return -1;

    double pi = AOF_2PI * 0.5;
    double phi2 = -2.0 * AOF_2PI;
    double phi1 = -1.0 * AOF_2PI;
    double u_w2 = 100.0 * phi2 + (-9.5) * phi2 * phi2 / (4.0 * pi);
    double u_w1 = 100.0 * phi1 + (-9.5) * phi1 * phi1 / (4.0 * pi);
    double du_out = CubeReg_deltaU(100.0, -9.5, phi2, -1);

    float verts[] = {
        0.0f, 119.0f, 0.0f,   6.0f, 119.0f, 0.0f,   3.0f, 119.0f, 2.0f,
        0.0f, 119.0f, 4.0f,   6.0f, 119.0f, 4.0f,   3.0f, 119.0f, 6.0f,
        0.0f, 119.0f, 8.0f,   6.0f, 119.0f, 8.0f,   3.0f, 119.0f, 10.0f,
        0.0f, 109.5f, 40.0f,  6.0f, 109.5f, 40.0f,  3.0f, 109.5f, 42.0f
    };
    double x_off[] = {0, 0, 2, 4, 4, 6, 8, 8, 10, 40, 40, 42};
    float registered_u[12], phi[12], uv[24], normals[36];
    int32_t faces[12], face_cube[4], gid[12];
    size_t cube_voff[5];
    long cube_org[4][3];
    char ids[4][48];
    for (int i = 0; i < 12; i++) {
        int on_wind1 = i >= 9;
        double base = on_wind1 ? u_w1 + (x_off[i] - 40.0)
                               : u_w2 + x_off[i];
        if (i >= 3 && i < 6) base += du_out;
        registered_u[i] = (float)base;
        phi[i] = (float)(on_wind1 ? phi1 : phi2);
    }
    PieceSet ps;
    aof_tabu_test_pieceset(&ps, verts, uv, phi, normals, faces, face_cube,
                           gid, cube_voff, cube_org, ids, 4);
    ScaffoldCalib cal;
    aof_tabu_test_calib(&cal);

    AtlasOverlapFixOptions opts;
    AtlasOverlapFixOptions_default(&opts);
    opts.tabu = 1;
    opts.wind_correct = 0;
    opts.tabu_span = 2;
    opts.tabu_tenure = 5;
    opts.tabu_max_iters = 60;
    opts.tabu_stall = 25;

    int failed = 0;
    double out_u[12], out_v[12];
    uint8_t face_keep[4];
    float u_iter[12];
    memcpy(u_iter, registered_u, sizeof u_iter);
    AtlasOverlapFixTabuReport report;
    AtlasOverlapFixStats stats;
    memset(&report, 0, sizeof report);
    memset(&stats, 0, sizeof stats);
    for (int round = 1; round <= 2 && !failed; round++) {
        if (AtlasOverlapFix_solve(arena, &ps, &cal, u_iter, &opts,
                                  out_u, out_v, face_keep, NULL,
                                  NULL, NULL, &report, &stats) != 0 ||
            report.nchart != 4) {
            fprintf(stderr, "[overlap_fix selftest] FAIL: handoff round %d "
                            "failed (ncharts=%zu)\n", round, report.nchart);
            failed = 1;
        }
        for (size_t i = 0; i < 12 && !failed; i++)
            u_iter[i] = (float)out_u[i];
    }
    if (!failed && stats.tabu_charts_moved != 0) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: handoff round 2 "
                        "still moved %zu charts\n", stats.tabu_charts_moved);
        failed = 1;
    }
    static const double want[4] = {2.0, 2.0, 2.0, 1.0};
    for (int c = 0; c < 4 && !failed; c++) {
        if (report.chart[c].k != 0 ||
            fabs(report.chart[c].wind - want[c]) > 0.25) {
            fprintf(stderr, "[overlap_fix selftest] FAIL: handoff chart %d "
                            "k=%d wind=%.3f want %.1f\n",
                    c, report.chart[c].k, report.chart[c].wind, want[c]);
            failed = 1;
        }
    }
    if (!failed)
        fprintf(stderr, "[overlap_fix selftest] rounds handoff ok "
                        "(wind {%.2f,%.2f,%.2f,%.2f} with k=0)\n",
                report.chart[0].wind, report.chart[1].wind,
                report.chart[2].wind, report.chart[3].wind);
    Arena_dispose(&arena);
    return failed ? -1 : 0;
}

/*
 * T2 -- a happy configuration stays put via best-restore.  Same chain as T1
 * with no displacement: every move is uphill, tabu still wanders (always-
 * accept), and the winner must be the untouched initial assignment.
 */
static int aof_tabu_stay_put_selftest(void)
{
    Arena_T arena = Arena_new();
    if (arena == NULL) return -1;

    double pi = AOF_2PI * 0.5;
    double phi2 = -2.0 * AOF_2PI;
    double u_w2 = 100.0 * phi2 + (-9.5) * phi2 * phi2 / (4.0 * pi);

    float verts[] = {
        0.0f, 119.0f, 0.0f,   6.0f, 119.0f, 0.0f,   3.0f, 119.0f, 2.0f,
        0.0f, 119.0f, 4.0f,   6.0f, 119.0f, 4.0f,   3.0f, 119.0f, 6.0f,
        0.0f, 119.0f, 8.0f,   6.0f, 119.0f, 8.0f,   3.0f, 119.0f, 10.0f
    };
    double x_off[] = {0, 0, 2,  4, 4, 6,  8, 8, 10};
    float registered_u[9], phi[9], uv[18], normals[27];
    int32_t faces[9], face_cube[3], gid[9];
    size_t cube_voff[4];
    long cube_org[3][3];
    char ids[3][48];
    for (int i = 0; i < 9; i++) {
        registered_u[i] = (float)(u_w2 + x_off[i]);
        phi[i] = (float)phi2;
    }
    PieceSet ps;
    aof_tabu_test_pieceset(&ps, verts, uv, phi, normals, faces, face_cube,
                           gid, cube_voff, cube_org, ids, 3);
    ScaffoldCalib cal;
    aof_tabu_test_calib(&cal);

    AtlasOverlapFixOptions opts;
    AtlasOverlapFixOptions_default(&opts);
    opts.tabu = 1;
    opts.wind_correct = 0;
    opts.tabu_span = 2;
    opts.tabu_tenure = 5;
    opts.tabu_max_iters = 40;
    opts.tabu_stall = 15;

    double out_u[9], out_v[9];
    uint8_t face_keep[3];
    AtlasOverlapFixTabuReport report;
    AtlasOverlapFixStats stats;
    int rc = AtlasOverlapFix_solve(arena, &ps, &cal, registered_u, &opts,
                                   out_u, out_v, face_keep, NULL,
                                   NULL, NULL, &report, &stats);
    int failed = 0;
    if (rc != 0) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: stay-put rc=%d\n", rc);
        failed = 1;
    }
    if (!failed && (report.chart[0].k != 0 || report.chart[1].k != 0 ||
                    report.chart[2].k != 0)) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: stay-put moved "
                        "k={%d,%d,%d}\n", report.chart[0].k,
                report.chart[1].k, report.chart[2].k);
        failed = 1;
    }
    if (!failed && (stats.tabu_best_iter != 0 ||
                    stats.tabu_energy_final != stats.tabu_energy_initial)) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: stay-put best_iter=%zu "
                        "E %.4f -> %.4f\n", stats.tabu_best_iter,
                stats.tabu_energy_initial, stats.tabu_energy_final);
        failed = 1;
    }
    if (!failed)
        fprintf(stderr, "[overlap_fix selftest] tabu stay-put ok\n");
    Arena_dispose(&arena);
    return failed ? -1 : 0;
}

/*
 * T3 -- the radial prior decides WHO moves.  Two identical-u charts stacked
 * one pitch apart in radius (no weld: 9.5 vox > the 6-vox weld distance),
 * plus a third settled chart anchoring the global median offset.  Clearing
 * the collision is symmetric -- inner out or outer out both work -- and only
 * the radial prior breaks the tie: the chart whose RADIUS says wind 3 while
 * its u says wind 2 is the one that must move out.
 */
static int aof_tabu_radial_outer_selftest(void)
{
    Arena_T arena = Arena_new();
    if (arena == NULL) return -1;

    double pi = AOF_2PI * 0.5;
    double phi2 = -2.0 * AOF_2PI;
    double u_w2 = 100.0 * phi2 + (-9.5) * phi2 * phi2 / (4.0 * pi);

    float verts[] = {
        /* I1: wind 2 in both radius and u */
        0.0f, 119.0f, 0.0f,   6.0f, 119.0f, 0.0f,   3.0f, 119.0f, 2.0f,
        /* O: radius one pitch OUT, u ON TOP of I1 */
        0.0f, 128.5f, 0.0f,   6.0f, 128.5f, 0.0f,   3.0f, 128.5f, 2.0f,
        /* I2: wind 2, elsewhere in u AND in v (z offset, NOT x -- an x
         * offset would inflate its RADIUS and poison the recentring
         * median, which is exactly the mistake this fixture once made) */
        40.0f, 119.0f, 0.0f,  46.0f, 119.0f, 0.0f,  43.0f, 119.0f, 2.0f
    };
    double u_off[] = {0, 0, 2,  0, 0, 2,  40, 40, 42};
    float registered_u[9], phi[9], uv[18], normals[27];
    int32_t faces[9], face_cube[3], gid[9];
    size_t cube_voff[4];
    long cube_org[3][3];
    char ids[3][48];
    for (int i = 0; i < 9; i++) {
        registered_u[i] = (float)(u_w2 + u_off[i]);
        phi[i] = (float)phi2;
    }
    PieceSet ps;
    aof_tabu_test_pieceset(&ps, verts, uv, phi, normals, faces, face_cube,
                           gid, cube_voff, cube_org, ids, 3);
    ScaffoldCalib cal;
    aof_tabu_test_calib(&cal);

    AtlasOverlapFixOptions opts;
    AtlasOverlapFixOptions_default(&opts);
    opts.tabu = 1;
    opts.tabu_span = 2;
    opts.tabu_tenure = 5;
    opts.tabu_max_iters = 60;
    opts.tabu_stall = 25;

    double out_u[9], out_v[9];
    uint8_t face_keep[3];
    AtlasOverlapFixTabuReport report;
    AtlasOverlapFixStats stats;
    int rc = AtlasOverlapFix_solve(arena, &ps, &cal, registered_u, &opts,
                                   out_u, out_v, face_keep, NULL,
                                   NULL, NULL, &report, &stats);
    int failed = 0;
    if (rc != 0 || stats.authentic_neighbour_edges != 0) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: radial-outer rc=%d "
                        "edges=%zu (expected none)\n", rc,
                rc == 0 ? stats.authentic_neighbour_edges : (size_t)0);
        failed = 1;
    }
    if (!failed && (report.chart[1].k != 1 || report.chart[0].k != 0 ||
                    report.chart[2].k != 0)) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: radial-outer "
                        "k={%d,%d,%d}, expected the OUTER chart at +1\n",
                report.chart[0].k, report.chart[1].k, report.chart[2].k);
        for (size_t c = 0; c < report.nchart; c++)
            fprintf(stderr, "    chart %zu: w_phi0=%.4f w_rad=%.4f u_med=%.2f"
                            " r_med=%.2f cells %lld->%lld\n",
                    c, report.chart[c].w_phi0, report.chart[c].w_rad,
                    report.chart[c].u_med, report.chart[c].radius_med,
                    (long long)report.chart[c].cells_initial,
                    (long long)report.chart[c].cells_final);
        for (size_t m = 0; m < report.nmove && m < 6; m++)
            fprintf(stderr, "    move %zu: iter %d %s %d: %d -> %d  dov=%.1f"
                            " dnb=%.1f drad=%.1f E=%.1f\n",
                    m, report.move[m].iter,
                    report.move[m].kind ? "group" : "chart",
                    report.move[m].id, report.move[m].from, report.move[m].to,
                    report.move[m].d_overlap, report.move[m].d_neighbour,
                    report.move[m].d_radial, report.move[m].energy_after);
        failed = 1;
    }
    if (!failed && (stats.overlap_pairs_after != 0 ||
                    stats.intra_overlap_pairs_after != 0)) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: radial-outer overlap "
                        "after %zu/%zu\n", stats.overlap_pairs_after,
                stats.intra_overlap_pairs_after);
        failed = 1;
    }
    if (!failed)
        fprintf(stderr, "[overlap_fix selftest] tabu radial-outer ok\n");
    Arena_dispose(&arena);
    return failed ? -1 : 0;
}

/* A fixed XYZ weld whose two phase gauges differ by one turn must retain that
 * source relation even though the current u values happen to meet. */
static int aof_tabu_phase_target_selftest(void)
{
    Arena_T arena = Arena_new();
    if (arena == NULL) return -1;

    double phi2 = -2.0 * AOF_2PI;
    float verts[] = {
        0.0f, 119.0f, 0.0f,  6.0f, 119.0f, 0.0f,  3.0f, 119.0f, 2.0f,
        0.0f, 119.0f, 4.0f,  6.0f, 119.0f, 4.0f,  3.0f, 119.0f, 6.0f
    };
    float registered_u[] = {0.0f, 6.0f, 3.0f, 4.0f, 10.0f, 7.0f};
    float phi[6], uv[12], normals[18];
    int32_t faces[6], face_cube[2], gid[6];
    size_t cube_voff[3];
    long cube_org[2][3];
    char ids[2][48];
    for (int i = 0; i < 6; i++)
        phi[i] = (float)(phi2 + (i >= 3 ? AOF_2PI : 0.0));

    PieceSet ps;
    aof_tabu_test_pieceset(&ps, verts, uv, phi, normals, faces, face_cube,
                           gid, cube_voff, cube_org, ids, 2);
    ScaffoldCalib cal;
    aof_tabu_test_calib(&cal);
    AtlasOverlapFixOptions opts;
    AtlasOverlapFixOptions_default(&opts);
    opts.phase_targets = 1;
    opts.tabu_span = 2;
    opts.lateral_w = 0.0;

    double u[6], v[6];
    uint8_t face_keep[2] = {1, 1};
    int32_t vertex_chart[6];
    for (int i = 0; i < 6; i++) {
        u[i] = (double)registered_u[i];
        v[i] = aof_axis_coord(&cal, &ps.verts[(size_t)i * 3]);
    }
    size_t ncharts = aof_label_charts(arena, &ps, face_keep, vertex_chart);
    AofNeighbour *edge = NULL;
    size_t nedges = 0;
    int failed = 0;
    if (ncharts != 2 ||
        aof_build_neighbours(arena, &ps, &cal, vertex_chart, u, &opts,
                             &edge, &nedges) != 0 || nedges != 1 ||
        edge[0].phase_target != -1 || edge[0].phase_residual > 1.0e-5) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: phase target charts=%zu "
                        "edges=%zu target=%d residual=%.6g\n",
                ncharts, nedges, nedges ? edge[0].phase_target : 99,
                nedges ? edge[0].phase_residual : -1.0);
        failed = 1;
    }
    if (!failed) {
        int32_t chart_group[2];
        size_t ngroups = aof_label_groups(arena, ncharts, edge, nedges,
                                          chart_group);
        AofTabu T;
        size_t excluded = 0;
        if (aof_tabu_prepare(arena, &ps, &cal, u, v, face_keep, vertex_chart,
                             ncharts, chart_group, ngroups, edge, nedges,
                             NULL, 1, &opts, 64, &T, &excluded) != 0 ||
            T.nedges != 1 || T.edge[0].target != -1 || excluded != 0) {
            fprintf(stderr, "[overlap_fix selftest] FAIL: phase target was "
                            "not preserved by tabu prepare\n");
            failed = 1;
        }
    }
    if (!failed)
        fprintf(stderr, "[overlap_fix selftest] tabu phase-target ok\n");
    Arena_dispose(&arena);
    return failed ? -1 : 0;
}

/*
 * T4 -- incremental bookkeeping == recount.  Drive the T1 problem through a
 * scripted move sequence against a deliberately tiny hash (forcing collision
 * chains and a mid-stream rebuild) and demand exact integer agreement between
 * the incrementally-maintained energy terms and a from-scratch recount after
 * every single move.
 */
static int aof_tabu_audit_selftest(void)
{
    Arena_T arena = Arena_new();
    if (arena == NULL) return -1;

    double pi = AOF_2PI * 0.5;
    double phi2 = -2.0 * AOF_2PI;
    double u_w2 = 100.0 * phi2 + (-9.5) * phi2 * phi2 / (4.0 * pi);
    double du_out = CubeReg_deltaU(100.0, -9.5, phi2, -1);

    float verts[] = {
        0.0f, 119.0f, 0.0f,   6.0f, 119.0f, 0.0f,   3.0f, 119.0f, 2.0f,
        0.0f, 119.0f, 4.0f,   6.0f, 119.0f, 4.0f,   3.0f, 119.0f, 6.0f,
        0.0f, 119.0f, 8.0f,   6.0f, 119.0f, 8.0f,   3.0f, 119.0f, 10.0f
    };
    double x_off[] = {0, 0, 2,  4, 4, 6,  8, 8, 10};
    float registered_u[9], phi[9], uv[18], normals[27];
    int32_t faces[9], face_cube[3], gid[9];
    size_t cube_voff[4];
    long cube_org[3][3];
    char ids[3][48];
    for (int i = 0; i < 9; i++) {
        double base = u_w2 + x_off[i];
        if (i >= 3 && i < 6) base += du_out;
        registered_u[i] = (float)base;
        phi[i] = (float)phi2;
    }
    PieceSet ps;
    aof_tabu_test_pieceset(&ps, verts, uv, phi, normals, faces, face_cube,
                           gid, cube_voff, cube_org, ids, 3);
    ScaffoldCalib cal;
    aof_tabu_test_calib(&cal);

    AtlasOverlapFixOptions opts;
    AtlasOverlapFixOptions_default(&opts);
    opts.tabu = 1;
    opts.tabu_span = 2;

    double u[9], v[9];
    uint8_t face_keep[3] = {1, 1, 1};
    int32_t vertex_chart[9];
    for (int i = 0; i < 9; i++) {
        u[i] = (double)registered_u[i];
        v[i] = aof_axis_coord(&cal, &ps.verts[(size_t)i * 3]);
    }
    size_t ncharts = aof_label_charts(arena, &ps, face_keep, vertex_chart);
    AofNeighbour *edge = NULL;
    size_t nedges = 0;
    int failed = 0;
    if (ncharts != 3 ||
        aof_build_neighbours(arena, &ps, &cal, vertex_chart, u, &opts,
                             &edge, &nedges) != 0 || nedges != 2) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: audit fixture charts=%zu"
                        " edges=%zu\n", ncharts, nedges);
        failed = 1;
    }
    if (!failed) {
        int32_t chart_group[3];
        size_t ngroups = aof_label_groups(arena, ncharts, edge, nedges,
                                          chart_group);
        AofTabu T;
        size_t excluded = 0;
        if (aof_tabu_prepare(arena, &ps, &cal, u, v, face_keep, vertex_chart,
                             ncharts, chart_group, ngroups, edge, nedges,
                             NULL, 0 /* wind_corrected */, &opts, 64,
                             &T, &excluded) != 0) {
            fprintf(stderr, "[overlap_fix selftest] FAIL: audit prepare\n");
            failed = 1;
        } else {
            static const int script[10][2] = {
                {1, -1}, {0, 1}, {1, 0}, {2, -1}, {0, 0},
                {2, 1}, {1, -2}, {2, 0}, {1, 0}, {0, -1}
            };
            for (int m = 0; m < 10 && !failed; m++) {
                int32_t c = script[m][0];
                int k = script[m][1];
                if (k == T.k[c]) k = k == 0 ? 1 : 0;   /* always a real move */
                int64_t d_nb = aof_tabu_nb_delta(&T, c, k, 0);
                int64_t d_rad =
                    T.radE[(size_t)c * (size_t)T.nk +
                           (size_t)(k + T.span)] -
                    T.radE[(size_t)c * (size_t)T.nk +
                           (size_t)(T.k[c] + T.span)];
                int64_t lost = aof_tabu_remove(&T, c, T.k[c]);
                T.k[c] = (int8_t)k;
                int64_t gained = aof_tabu_insert(&T, c, k);
                T.e_ov += gained - lost;
                T.e_nb += d_nb;
                T.e_rad += d_rad;
                if (m == 5) aof_tabu_rebuild(&T);   /* mid-stream reclaim */
                if (aof_tabu_audit(&T) != 0) {
                    fprintf(stderr, "[overlap_fix selftest] FAIL: audit "
                                    "mismatch after move %d (chart %d -> "
                                    "%d)\n", m, c, k);
                    failed = 1;
                }
            }
        }
    }
    if (!failed)
        fprintf(stderr, "[overlap_fix selftest] tabu audit ok\n");
    Arena_dispose(&arena);
    return failed ? -1 : 0;
}

/*
 * T5 -- the tenure ban.  One collision-free chart with lambda_rad = 0 makes
 * every move cost exactly zero, so a broken tabu list would ping-pong home
 * immediately.  The correct walk departs 0, is REFUSED the return while the
 * tenure holds (aspiration cannot fire: energy never improves), takes the
 * only other free value, then runs out of admissible moves and stops.
 */
static int aof_tabu_tenure_selftest(void)
{
    Arena_T arena = Arena_new();
    if (arena == NULL) return -1;

    double pi = AOF_2PI * 0.5;
    double phi2 = -2.0 * AOF_2PI;
    double u_w2 = 100.0 * phi2 + (-9.5) * phi2 * phi2 / (4.0 * pi);

    float verts[] = {
        0.0f, 119.0f, 0.0f,   6.0f, 119.0f, 0.0f,   3.0f, 119.0f, 2.0f
    };
    double x_off[] = {0, 0, 2};
    float registered_u[3], phi[3], uv[6], normals[9];
    int32_t faces[3], face_cube[1], gid[3];
    size_t cube_voff[2];
    long cube_org[1][3];
    char ids[1][48];
    for (int i = 0; i < 3; i++) {
        registered_u[i] = (float)(u_w2 + x_off[i]);
        phi[i] = (float)phi2;
    }
    PieceSet ps;
    aof_tabu_test_pieceset(&ps, verts, uv, phi, normals, faces, face_cube,
                           gid, cube_voff, cube_org, ids, 1);
    ScaffoldCalib cal;
    aof_tabu_test_calib(&cal);

    AtlasOverlapFixOptions opts;
    AtlasOverlapFixOptions_default(&opts);
    opts.tabu = 1;
    opts.tabu_span = 1;
    opts.tabu_tenure = 3;
    opts.tabu_max_iters = 10;
    opts.tabu_stall = 10;
    opts.lambda_rad = 0.0;      /* every move is free: pure tabu mechanics */

    double out_u[3], out_v[3];
    uint8_t face_keep[1];
    AtlasOverlapFixTabuReport report;
    AtlasOverlapFixStats stats;
    int rc = AtlasOverlapFix_solve(arena, &ps, &cal, registered_u, &opts,
                                   out_u, out_v, face_keep, NULL,
                                   NULL, NULL, &report, &stats);
    int failed = 0;
    if (rc != 0 || report.nmove != 2) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: tenure rc=%d moves=%zu "
                        "(expected exactly 2: depart, blocked-return "
                        "sidestep, then all banned)\n",
                rc, rc == 0 ? report.nmove : (size_t)0);
        failed = 1;
    }
    if (!failed && report.move[1].to == 0) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: tenure ban ignored -- "
                        "returned to the departed depth at iter %d\n",
                report.move[1].iter);
        failed = 1;
    }
    if (!failed && (stats.tabu_iterations != 2 || report.chart[0].k != 0)) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: tenure iters=%zu "
                        "final k=%d (best-restore must return 0)\n",
                stats.tabu_iterations, report.chart[0].k);
        failed = 1;
    }
    if (!failed)
        fprintf(stderr, "[overlap_fix selftest] tabu tenure ok\n");
    Arena_dispose(&arena);
    return failed ? -1 : 0;
}

/*
 * T10 -- the graded rewind unrolls a collapsed group in one move.  A/R/B/C
 * form one welded chain climbing from wind 2 to wind 3, but the registered u
 * puts ALL of them in wind 2's band: B sits exactly on A and C on R --
 * wrong-wind intra pairs, so the group is denied family status and flagged
 * graded.  Per-chart escapes each tear two welds through worse intermediate
 * states; the graded move d=0 applies bins {0,0,1,1} in one step, tearing
 * only the R-B bin boundary.  Expected: k = {0,0,+1,+1}, zero overlap after,
 * at least one graded move fired.
 */
static int aof_tabu_graded_selftest(void)
{
    Arena_T arena = Arena_new();
    if (arena == NULL) return -1;

    double pi = AOF_2PI * 0.5;
    double phi2 = -2.0 * AOF_2PI;
    double u_w2 = 100.0 * phi2 + (-9.5) * phi2 * phi2 / (4.0 * pi);

    /* Triangles are LARGE (20-vox extent) so the collision term carries
     * realistic weight against the weld tears -- with toy 2-vox faces the
     * tear cost swamps the relief and nothing can afford to unroll. */
    /* Radii chosen so the rounded-wind bins match the physical truth:
     * A, A2 (120) and R (123) round to wind 2, B and C (129) to wind 3, and
     * the weld chain closes (each hop <= 6 vox).  B collapsed onto A and C
     * onto R; the graded bins {0,0,0,1,1} move the wind-3 pair out as one,
     * tearing only the R-B bin boundary.  A2 exists so the correctly-placed
     * members are the MAJORITY: the radial recentring anchors on them and
     * the prior then picks the outward gauge over the identical-energy
     * inward one (with a 2-vs-2 group the median is genuinely ambiguous). */
    float verts[] = {
        /* A: wind 2 */
        0.0f, 120.0f, 0.0f,   20.0f, 120.0f, 0.0f,  0.0f, 120.0f, 6.0f,
        /* A2: wind 2, welded to A, correctly placed */
        0.0f, 120.0f, 6.0f,   20.0f, 120.0f, 6.0f,  0.0f, 120.0f, 12.0f,
        /* R: wind 2's outer edge, welded to A and B */
        0.0f, 123.0f, 0.0f,   20.0f, 123.0f, 0.0f,  0.0f, 123.0f, 6.0f,
        /* B: wind 3, u collapsed onto A */
        0.0f, 129.0f, 0.0f,   20.0f, 129.0f, 0.0f,  0.0f, 129.0f, 6.0f,
        /* C: wind 3, welded to B, u collapsed onto R */
        0.0f, 129.0f, 6.0f,   20.0f, 129.0f, 6.0f,  0.0f, 129.0f, 12.0f
    };
    double u_off[] = {0, 20, 10,  48, 68, 58,  24, 44, 34,
                      0, 20, 10,  24, 44, 34};
    float registered_u[15], phi[15], uv[30], normals[45];
    int32_t faces[15], face_cube[5], gid[15];
    size_t cube_voff[6];
    long cube_org[5][3];
    char ids[5][48];
    for (int i = 0; i < 15; i++) {
        registered_u[i] = (float)(u_w2 + u_off[i]);
        phi[i] = (float)phi2;
    }
    PieceSet ps;
    aof_tabu_test_pieceset(&ps, verts, uv, phi, normals, faces, face_cube,
                           gid, cube_voff, cube_org, ids, 5);
    ScaffoldCalib cal;
    aof_tabu_test_calib(&cal);

    AtlasOverlapFixOptions opts;
    AtlasOverlapFixOptions_default(&opts);
    opts.tabu = 1;
    opts.graded_targets = 1;
    opts.lock_graded = 1;
    /* Exercise enough speculative offsets to require mid-sweep occupancy
     * reclamation; the production phase-bin failure first appeared at 16. */
    opts.tabu_span = 16;
    opts.tabu_tenure = 5;
    opts.tabu_max_iters = 80;
    opts.tabu_stall = 30;
    /* A 3-vertex fixture face occupies ~4 proxy cells no matter its area
     * (cells come from vertices + centroid), so the collision relief here is
     * fixture-scale, not chart-scale.  Rebalance the weld weight to the
     * PROPORTIONS a real chart sees: collision must be able to pay for one
     * bin-boundary tear.  The radial prior stays at its default -- it is
     * what picks the outward gauge over the identical-energy inward one. */
    opts.lambda_nb = 1.0;

    double out_u[15], out_v[15];
    uint8_t face_keep[5];
    AtlasOverlapFixGroup *garr = NULL;
    size_t ngr = 0;
    AtlasOverlapFixTabuReport report;
    AtlasOverlapFixStats stats;
    int rc = AtlasOverlapFix_solve(arena, &ps, &cal, registered_u, &opts,
                                   out_u, out_v, face_keep, NULL,
                                   &garr, &ngr, &report, &stats);
    int failed = 0;
    if (rc != 0 || stats.total_groups != 1 ||
        stats.tabu_graded_groups != 1) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: graded rc=%d groups=%zu"
                        " eligible=%zu (expected one collapsed group; "
                        "intra=%zu wrongwind=%zu edges=%zu killed=%zu "
                        "charts=%zu raw=%zu/%zu before=%zu/%zu)\n", rc,
                rc == 0 ? stats.total_groups : (size_t)0,
                rc == 0 ? stats.tabu_graded_groups : (size_t)0,
                rc == 0 && ngr > 0 ? garr[0].intra_overlap : (size_t)0,
                rc == 0 && ngr > 0 ? garr[0].intra_wrongwind : (size_t)0,
                rc == 0 ? stats.authentic_neighbour_edges : (size_t)0,
                rc == 0 ? stats.faces_killed : (size_t)0,
                rc == 0 ? stats.total_charts : (size_t)0,
                rc == 0 ? stats.overlap_pairs_raw : (size_t)0,
                rc == 0 ? stats.intra_overlap_pairs_raw : (size_t)0,
                rc == 0 ? stats.overlap_pairs_before : (size_t)0,
                rc == 0 ? stats.intra_overlap_pairs : (size_t)0);
        failed = 1;
    }
    if (!failed && (report.chart[0].k != 0 || report.chart[1].k != 0 ||
                    report.chart[2].k != 0 || report.chart[3].k != 1 ||
                    report.chart[4].k != 1)) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: graded "
                        "k={%d,%d,%d,%d,%d}, expected {0,0,0,1,1}\n",
                report.chart[0].k, report.chart[1].k, report.chart[2].k,
                report.chart[3].k, report.chart[4].k);
        for (size_t m = 0; m < report.nmove && m < 8; m++)
            fprintf(stderr, "    move %zu: iter %d kind %d id %d %d->%d "
                            "dov=%.1f dnb=%.1f drad=%.1f E=%.1f\n",
                    m, report.move[m].iter, report.move[m].kind,
                    report.move[m].id, report.move[m].from,
                    report.move[m].to, report.move[m].d_overlap,
                    report.move[m].d_neighbour, report.move[m].d_radial,
                    report.move[m].energy_after);
        failed = 1;
    }
    if (!failed && (stats.overlap_pairs_after != 0 ||
                    stats.intra_overlap_pairs_after != 0)) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: graded overlap after "
                        "%zu cross / %zu intra\n", stats.overlap_pairs_after,
                stats.intra_overlap_pairs_after);
        failed = 1;
    }
    if (!failed && stats.tabu_welds_torn != 0) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: graded target left %zu "
                        "welds torn (expected zero)\n",
                stats.tabu_welds_torn);
        failed = 1;
    }
    if (!failed && stats.tabu_moves_graded < 1) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: graded solved without "
                        "a graded move (%zu) -- the move class is dead\n",
                stats.tabu_moves_graded);
        failed = 1;
    }
    if (!failed)
        fprintf(stderr, "[overlap_fix selftest] tabu graded ok\n");
    Arena_dispose(&arena);
    return failed ? -1 : 0;
}

/*
 * T9 -- the relayout runs after the tabu repair and must not damage it.
 * Same displaced chain as T1 with relayout on: B comes home through the
 * welds (its u-ladder wind is 3, so the fallback step is the exact
 * wind-3-to-2 arc), then the register pass sees only healthy springs and
 * holds the chain assembled at its physical x offsets.
 */
static int aof_tabu_relayout_selftest(void)
{
    Arena_T arena = Arena_new();
    if (arena == NULL) return -1;

    double pi = AOF_2PI * 0.5;
    double phi2 = -2.0 * AOF_2PI;
    double u_w2 = 100.0 * phi2 + (-9.5) * phi2 * phi2 / (4.0 * pi);
    double du_out = CubeReg_deltaU(100.0, -9.5, phi2, -1);

    float verts[] = {
        0.0f, 119.0f, 0.0f,   6.0f, 119.0f, 0.0f,   3.0f, 119.0f, 2.0f,
        0.0f, 119.0f, 4.0f,   6.0f, 119.0f, 4.0f,   3.0f, 119.0f, 6.0f,
        0.0f, 119.0f, 8.0f,   6.0f, 119.0f, 8.0f,   3.0f, 119.0f, 10.0f
    };
    double x_off[] = {0, 0, 2,  4, 4, 6,  8, 8, 10};
    float registered_u[9], phi[9], uv[18], normals[27];
    int32_t faces[9], face_cube[3], gid[9];
    size_t cube_voff[4];
    long cube_org[3][3];
    char ids[3][48];
    for (int i = 0; i < 9; i++) {
        double base = u_w2 + x_off[i];
        if (i >= 3 && i < 6) base += du_out;
        registered_u[i] = (float)base;
        phi[i] = (float)phi2;
    }
    PieceSet ps;
    aof_tabu_test_pieceset(&ps, verts, uv, phi, normals, faces, face_cube,
                           gid, cube_voff, cube_org, ids, 3);
    ScaffoldCalib cal;
    aof_tabu_test_calib(&cal);

    AtlasOverlapFixOptions opts;
    AtlasOverlapFixOptions_default(&opts);
    opts.tabu = 1;
    opts.wind_correct = 0;
    opts.tabu_span = 2;
    opts.tabu_tenure = 5;
    opts.tabu_max_iters = 60;
    opts.tabu_stall = 25;
    opts.relayout = 1;

    double out_u[9], out_v[9];
    uint8_t face_keep[3];
    AtlasOverlapFixTabuReport report;
    AtlasOverlapFixStats stats;
    int rc = AtlasOverlapFix_solve(arena, &ps, &cal, registered_u, &opts,
                                   out_u, out_v, face_keep, NULL,
                                   NULL, NULL, &report, &stats);
    int failed = 0;
    if (rc != 0 || report.chart[1].k != -1 || stats.relayout_failed) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: relayout rc=%d k1=%d "
                        "failed=%d\n", rc,
                rc == 0 ? report.chart[1].k : 0,
                rc == 0 ? stats.relayout_failed : -1);
        failed = 1;
    }
    if (!failed && stats.relayout_edges_weld != 2) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: relayout weld springs "
                        "%zu != 2\n", stats.relayout_edges_weld);
        failed = 1;
    }
    if (!failed && (fabs((out_u[3] - out_u[0]) - 4.0) > 3.0 ||
                    fabs((out_u[6] - out_u[0]) - 8.0) > 3.0)) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: relayout chain gaps "
                        "%.2f / %.2f, expected ~4 / ~8\n",
                out_u[3] - out_u[0], out_u[6] - out_u[0]);
        fprintf(stderr, "    relayout reverted=%d exact=%zu/%zu cross=%zu/%zu "
                        "rank={%zu,%zu,%zu} du={%.2f,%.2f,%.2f}\n",
                stats.relayout_reverted, stats.relayout_exact_pairs_before,
                stats.relayout_exact_pairs_after,
                stats.relayout_cross_pairs_before,
                stats.relayout_cross_pairs_after, report.chart[0].relayout_rank,
                report.chart[1].relayout_rank, report.chart[2].relayout_rank,
                report.chart[0].relayout_shift_u,
                report.chart[1].relayout_shift_u,
                report.chart[2].relayout_shift_u);
        failed = 1;
    }
    if (!failed)
        fprintf(stderr, "[overlap_fix selftest] tabu relayout ok\n");
    Arena_dispose(&arena);
    return failed ? -1 : 0;
}

/*
 * T7 -- lateral coupling decides who moves.  L1 and L2 sit at the SAME
 * radius in adjacent cubes, side by side in u but 14+ vox apart in XYZ (a
 * hole in the sheet) -- no weld exists, only a lateral edge.  O duplicates
 * L1's (u, v) footprint from a non-adjacent cube with the radial prior
 * switched off, so clearing the collision by moving L1 or by moving O costs
 * exactly the same in every term EXCEPT the lateral tear -- and the scan
 * order would pick L1 (lower id) if the lateral edge did not exist.
 */
static int aof_tabu_lateral_selftest(void)
{
    Arena_T arena = Arena_new();
    if (arena == NULL) return -1;

    double pi = AOF_2PI * 0.5;
    double phi2 = -2.0 * AOF_2PI;
    double u_w2 = 100.0 * phi2 + (-9.5) * phi2 * phi2 / (4.0 * pi);

    float verts[] = {
        /* L1 */
        0.0f, 119.0f, 0.0f,   6.0f, 119.0f, 0.0f,   3.0f, 119.0f, 2.0f,
        /* L2: same radius, z-offset (a hole away in XYZ, adjacent in u) */
        20.0f, 119.0f, 0.0f,  26.0f, 119.0f, 0.0f,  23.0f, 119.0f, 2.0f,
        /* O: duplicates L1's (u, v) from a NON-adjacent cube */
        0.0f, 119.0f, 40.0f,  6.0f, 119.0f, 40.0f,  3.0f, 119.0f, 42.0f
    };
    double u_off[] = {0, 0, 2,  4, 4, 6,  0, 0, 2};
    float registered_u[9], phi[9], uv[18], normals[27];
    int32_t faces[9], face_cube[3], gid[9];
    size_t cube_voff[4];
    long cube_org[3][3];
    char ids[3][48];
    for (int i = 0; i < 9; i++) {
        registered_u[i] = (float)(u_w2 + u_off[i]);
        phi[i] = (float)phi2;
    }
    PieceSet ps;
    aof_tabu_test_pieceset(&ps, verts, uv, phi, normals, faces, face_cube,
                           gid, cube_voff, cube_org, ids, 3);
    /* O's verts sit at v = 0..6 like L1's (z is v), but its cube must not be
     * grid-adjacent to L1's or the lateral scan would couple them too. */
    cube_org[2][2] = 384;
    ScaffoldCalib cal;
    aof_tabu_test_calib(&cal);

    AtlasOverlapFixOptions opts;
    AtlasOverlapFixOptions_default(&opts);
    opts.tabu = 1;
    opts.tabu_span = 2;
    opts.tabu_tenure = 5;
    opts.tabu_max_iters = 60;
    opts.tabu_stall = 25;
    opts.lambda_rad = 0.0;    /* neutral prior: only the lateral tear differs */

    double out_u[9], out_v[9];
    uint8_t face_keep[3];
    AtlasOverlapFixTabuReport report;
    AtlasOverlapFixStats stats;
    int rc = AtlasOverlapFix_solve(arena, &ps, &cal, registered_u, &opts,
                                   out_u, out_v, face_keep, NULL,
                                   NULL, NULL, &report, &stats);
    int failed = 0;
    if (rc != 0 || stats.authentic_neighbour_edges != 0 ||
        stats.tabu_lateral_edges != 1) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: lateral rc=%d welds=%zu"
                        " lateral=%zu (expected 0 welds, 1 lateral)\n", rc,
                rc == 0 ? stats.authentic_neighbour_edges : (size_t)0,
                rc == 0 ? stats.tabu_lateral_edges : (size_t)0);
        failed = 1;
    }
    if (!failed && (report.chart[2].k == 0 || report.chart[0].k != 0 ||
                    report.chart[1].k != 0)) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: lateral k={%d,%d,%d} "
                        "-- O must move, the coupled pair must not\n",
                report.chart[0].k, report.chart[1].k, report.chart[2].k);
        failed = 1;
    }
    if (!failed && (stats.tabu_lateral_torn != 0 ||
                    stats.overlap_pairs_after != 0)) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: lateral torn=%zu "
                        "overlap after=%zu\n", stats.tabu_lateral_torn,
                stats.overlap_pairs_after);
        failed = 1;
    }
    if (!failed)
        fprintf(stderr, "[overlap_fix selftest] tabu lateral ok\n");
    Arena_dispose(&arena);
    return failed ? -1 : 0;
}

/*
 * Four charts on an axis along z, so v is just the z coordinate.
 *
 *   chart 0  (u,v) = (0,0) (10,0) (0,10)      at radius ~100
 *   chart 1  (u,v) = (5,0) (15,0) (5,10)      at radius ~200  -- overlaps 0
 *   chart 2  (u,v) = (100,0) (110,0) (100,10) at radius ~150  -- clean
 *   chart 3  (u,v) = (110,0) (120,0) (110,10) at radius ~152  -- clean, and
 *                                             2 vox from chart 2 in XYZ, so
 *                                             the two are authentic neighbours
 *                                             and must land in one group.
 *
 * Expected: 4 charts, 3 groups (2 and 3 merge), one overlapping pair, and the
 * furthest-out offender (chart 1, r=200) shifts a whole wrap -- 2*pi*200 is far
 * larger than the atlas, so one wrap clears it and nothing is left overlapping.
 */
int AtlasOverlapFix_selftest(void)
{
    Arena_T arena = Arena_new();
    if (arena == NULL) return -1;
    if (AtlasCollisionRegister_selftest() != 0) {
        Arena_dispose(&arena);
        return -1;
    }
    if (aof_triangle_cull_selftest() != 0 ||
        aof_same_cube_weld_selftest() != 0 ||
        aof_tabu_return_home_selftest() != 0 ||
        aof_rounds_handoff_selftest() != 0 ||
        aof_tabu_stay_put_selftest() != 0 ||
        aof_tabu_radial_outer_selftest() != 0 ||
        aof_tabu_phase_target_selftest() != 0 ||
        aof_tabu_audit_selftest() != 0 ||
        aof_tabu_tenure_selftest() != 0 ||
        aof_tabu_lateral_selftest() != 0 ||
        aof_tabu_relayout_selftest() != 0 ||
        aof_tabu_graded_selftest() != 0) {
        Arena_dispose(&arena);
        return -1;
    }

    float verts[] = {
        /* chart 0: cube 0, y ~ 100 */
        0.0f, 100.0f, 0.0f,   0.0f, 100.0f, 10.0f,  10.0f, 100.0f, 0.0f,
        /* chart 1: cube 1, y ~ 200 */
        0.0f, 200.0f, 0.0f,   0.0f, 200.0f, 10.0f,  10.0f, 200.0f, 0.0f,
        /* chart 2: cube 2, y ~ 150 */
        0.0f, 150.0f, 0.0f,   0.0f, 150.0f, 10.0f,  10.0f, 150.0f, 0.0f,
        /* chart 3: cube 3, y ~ 152 -- within a weld of chart 2 */
        0.0f, 152.0f, 0.0f,   0.0f, 152.0f, 10.0f,  10.0f, 152.0f, 0.0f
    };
    int32_t faces[] = {0,1,2,  3,4,5,  6,7,8,  9,10,11};
    int32_t face_cube[] = {0, 1, 2, 3};
    size_t cube_voff[] = {0, 3, 6, 9, 12};
    long cube_org[][3] = {{0,0,0}, {0,128,0}, {0,256,0}, {0,384,0}};
    char ids[4][48];
    for (int i = 0; i < 4; i++)
        snprintf(ids[i], 48, "z00000_y%05d_x00000", i * 128);
    float uv[24];
    float normals[36];
    float phi[12];
    int32_t gid[12];
    memset(uv, 0, sizeof(uv));
    memset(normals, 0, sizeof(normals));
    memset(phi, 0, sizeof(phi));
    for (int i = 0; i < 12; i++) gid[i] = -1;

    float registered_u[] = {
        0.0f, 10.0f, 0.0f,
        5.0f, 15.0f, 5.0f,
        100.0f, 110.0f, 100.0f,
        110.0f, 120.0f, 110.0f
    };
    for (int i = 0; i < 12; i++) uv[i * 2] = registered_u[i];

    PieceSet ps;
    memset(&ps, 0, sizeof(ps));
    ps.verts = verts;
    ps.uv = uv;
    ps.phi = phi;
    ps.normals = normals;
    ps.faces = faces;
    ps.face_cube = face_cube;
    ps.gid = gid;
    ps.nv = 12;
    ps.nf = 4;
    ps.ids = ids;
    ps.cube_voff = cube_voff;
    ps.cube_org = cube_org;
    ps.n_cubes = 4;

    ScaffoldCalib cal;
    memset(&cal, 0, sizeof(cal));
    cal.pitch = 9.5;
    cal.spiral_a = 1.0;
    cal.spiral_b = 9.5;
    cal.sense = -1;
    cal.axis_dir[0] = 1.0f;   /* axis along z, so v == z */

    double out_u[12], out_v[12];
    uint8_t face_keep[4];
    int32_t vertex_group[12];
    AtlasOverlapFixGroup *group = NULL;
    size_t ngroups = 0;
    AtlasOverlapFixStats stats;

    AtlasOverlapFixOptions opts;
    AtlasOverlapFixOptions_default(&opts);
    opts.neighbour_min_shared = 3;

    int rc = AtlasOverlapFix_solve(arena, &ps, &cal, registered_u, &opts,
                                   out_u, out_v, face_keep, vertex_group,
                                   &group, &ngroups, NULL, &stats);

    int failed = 0;
    if (rc != 0) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: rc=%d\n", rc);
        failed = 1;
    }
    if (!failed && stats.faces_killed != 0) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: killed %zu good faces\n",
                stats.faces_killed);
        failed = 1;
    }
    if (!failed && stats.total_charts != 4) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: charts %zu != 4\n",
                stats.total_charts);
        failed = 1;
    }
    if (!failed && stats.authentic_neighbour_edges != 1) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: neighbour edges %zu != 1\n",
                stats.authentic_neighbour_edges);
        failed = 1;
    }
    if (!failed && stats.total_groups != 3) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: groups %zu != 3 "
                        "(charts 2 and 3 must weld into one)\n",
                stats.total_groups);
        failed = 1;
    }
    if (!failed && stats.overlap_pairs_before != 1) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: overlap before %zu != 1\n",
                stats.overlap_pairs_before);
        failed = 1;
    }
    if (!failed && stats.overlap_pairs_after != 0) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: overlap after %zu != 0\n",
                stats.overlap_pairs_after);
        failed = 1;
    }
    if (!failed && stats.shifted_groups != 1) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: shifted %zu != 1\n",
                stats.shifted_groups);
        failed = 1;
    }
    if (!failed) {
        /* The offender must be the one furthest from the axis. */
        int32_t moved = -1;
        for (size_t g = 0; g < ngroups; g++)
            if (group[g].action == ATLAS_OVERLAP_FIX_SHIFTED)
                moved = (int32_t)g;
        if (moved < 0 || group[moved].mean_radius < 150.0) {
            fprintf(stderr, "[overlap_fix selftest] FAIL: moved group %d at "
                            "radius %.1f, expected the ~200 one\n",
                    moved, moved >= 0 ? group[moved].mean_radius : -1.0);
            failed = 1;
        } else if (abs(group[moved].wraps) != 1) {
            fprintf(stderr, "[overlap_fix selftest] FAIL: wraps %d, expected "
                            "+/-1\n", group[moved].wraps);
            failed = 1;
        }
    }
    if (!failed && stats.parked_groups != 0) {
        fprintf(stderr, "[overlap_fix selftest] FAIL: parked %zu, expected 0 "
                        "(a wrap shift should have sufficed)\n",
                stats.parked_groups);
        failed = 1;
    }

    /* T6 -- the same fixture through the tabu strategy.  The 0/1 collision
     * must clear by re-gauging exactly one of the two (the radial prior
     * favours the r=200 outlier), the 2-3 weld must hold, nothing parks. */
    if (!failed) {
        AtlasOverlapFixOptions topts;
        AtlasOverlapFixOptions_default(&topts);
        topts.neighbour_min_shared = 3;
        topts.tabu = 1;
        topts.tabu_span = 2;
        topts.tabu_tenure = 5;
        topts.tabu_max_iters = 80;
        topts.tabu_stall = 30;
        double tab_u[12], tab_v[12];
        uint8_t tab_keep[4];
        AtlasOverlapFixTabuReport report;
        AtlasOverlapFixStats tstats;
        int trc = AtlasOverlapFix_solve(arena, &ps, &cal, registered_u,
                                        &topts, tab_u, tab_v, tab_keep,
                                        NULL, NULL, NULL, &report, &tstats);
        int moved01 = 0;
        if (trc == 0 && report.nchart == 4)
            moved01 = (report.chart[0].k != 0) + (report.chart[1].k != 0);
        if (trc != 0 || report.nchart != 4) {
            fprintf(stderr, "[overlap_fix selftest] FAIL: tabu pipeline "
                            "rc=%d ncharts=%zu\n", trc,
                    trc == 0 ? report.nchart : (size_t)0);
            failed = 1;
        } else if (tstats.overlap_pairs_after != 0 ||
                   tstats.intra_overlap_pairs_after != 0) {
            fprintf(stderr, "[overlap_fix selftest] FAIL: tabu pipeline "
                            "overlap after %zu cross / %zu intra\n",
                    tstats.overlap_pairs_after,
                    tstats.intra_overlap_pairs_after);
            failed = 1;
        } else if (moved01 != 1 || report.chart[2].k != 0 ||
                   report.chart[3].k != 0) {
            fprintf(stderr, "[overlap_fix selftest] FAIL: tabu pipeline "
                            "k={%d,%d,%d,%d}: exactly one of charts 0/1 "
                            "must move, 2/3 must not\n",
                    report.chart[0].k, report.chart[1].k,
                    report.chart[2].k, report.chart[3].k);
            failed = 1;
        } else if (tstats.tabu_welds_torn != 0 || tstats.parked_groups != 0) {
            fprintf(stderr, "[overlap_fix selftest] FAIL: tabu pipeline "
                            "torn=%zu parked=%zu\n",
                    tstats.tabu_welds_torn, tstats.parked_groups);
            failed = 1;
        } else if (tstats.tabu_family_edges != 1 ||
                   tstats.tabu_family_torn != 0) {
            fprintf(stderr, "[overlap_fix selftest] FAIL: tabu pipeline "
                            "family bonds=%zu torn=%zu (charts 2+3 are a "
                            "happy family)\n", tstats.tabu_family_edges,
                    tstats.tabu_family_torn);
            failed = 1;
        } else {
            fprintf(stderr, "[overlap_fix selftest] tabu pipeline ok "
                            "(moved chart %d to k=%d)\n",
                    report.chart[0].k != 0 ? 0 : 1,
                    report.chart[0].k != 0 ? report.chart[0].k
                                           : report.chart[1].k);
        }
    }

    if (!failed)
        fprintf(stderr,
            "[overlap_fix selftest] PASS: charts=%zu groups=%zu edges=%zu "
            "overlap=%zu->%zu settled=%zu shifted=%zu parked=%zu\n",
            stats.total_charts, stats.total_groups,
            stats.authentic_neighbour_edges, stats.overlap_pairs_before,
            stats.overlap_pairs_after, stats.settled_groups,
            stats.shifted_groups, stats.parked_groups);
    Arena_dispose(&arena);
    return failed ? -1 : 0;
}
