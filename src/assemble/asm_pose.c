/* asm_pose.c -- rigid SE(2) pose graph over charts (TAUCS back-end). */
#include "asm_pose.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/pipeline_constants.h"
#include "../common/union_find.h"
#include "../common/ves_platform.h"
#include "../flatten/sparse_solve.h"

#define AP_PI 3.14159265358979323846
/* WORLD-Z PRIOR (2026-09-09): the pose graph carried only relative constraints and one anchor
 * per component, so a long chain of seams drifts freely in rotation and v -- measured on the
 * 4x21x21 and 10x10x10 primaries as local tilts of -8..+6 deg along u and 54-98 vox (p50) of z
 * wander from the best plane, against 27 vox on the familiar.  On a vertical sheet world z is
 * a Euclidean coordinate of the surface, so every chart's own z-gradient direction is an
 * absolute measurement of its rotation and its z an absolute measurement of its position
 * along that direction.  Both enter the Gauss-Newton solve as priors, for charts whose
 * gradient magnitude says the sheet is vertical there (a tilted region gets none). */
/* A relation's information weights come from its correspondences' scatter alone (w_xy =
 * n / var), so a stiff 100-pair seam claims a 0.05-vox precision and the 5-sigma switch
 * fires at 0.25 vox: on the familiar 4x5x5 it switched 44 RIGHT joins, and with the
 * components following the switching (bug 3) their pieces were re-placed one band off
 * (46,985 stacked (u,z) cells), while on the 10x10x10 the switched relations were WRONG
 * joins (392,409 stacked cells when they were kept).  A systematic floor on every seam's
 * precision -- the measured seam floor is ~1.6 vox rms, flattening distortion included --
 * makes the switch a disagreement in vox: both weights scale by 1 / (1 + floor^2 w_xy).
 * 0 = off (A/B, 2026-09-09). */
#define AP_SEAM_FLOOR_VOX 1.0
/* Where the floor applies.  0: to the SWITCH decision only (the robust weight is computed
 * from the floored residual, the least squares keeps every seam at its measured precision);
 * 1: to the solve's weights as well -- measured on the familiar (dev47, 2026-09-09) it
 * loosened every seam (xy residual p50 0.13 -> 0.53 vox, lattice edges over 6 vox
 * 449 -> 3,636) for the same switching gain. */
#define AP_SEAM_FLOOR_IN_SOLVE 0
/* POST-PLACEMENT REFINE (AsmPose_refine, 2026-09-10): plain least squares for the first iterations
 * (a piece 100 vox off its seam has a normalized residual in the hundreds at iteration 0, so IRLS
 * would silence the very seam that must act), then Cauchy; a hold prior on the warm start keeps the
 * system SPD and leaves a chart's rotation with its z-gradient orientation unless a long seam says
 * otherwise (a 30-vox seam's w_theta ~ 18,000 wins, a 12-vox SHORT seam's ~900 barely). */
#define AP_REFINE_ITERS 12
#define AP_REFINE_LS_ITERS 3
#define AP_REFINE_HOLD_XY_VOX 50.0
#define AP_REFINE_HOLD_THETA_RAD 0.05
#define AP_COMPONENTS_AFTER_SWITCH 1   /* A/B: 0 = components over the pre-switch graph (reviewer bug 3, 2026-09-09) */
#define AP_ZPRIOR 1
#define AP_ZPRIOR_MIN_GRAD 0.85
#define AP_ZPRIOR_W_THETA 200.0     /* 1 / (0.07 rad)^2 */
#define AP_ZPRIOR_W_Z 0.04          /* 1 / (5 vox)^2, per seam relation (v2: relative) */
#define AP_ZPRIOR_ABS 0             /* 1 = the absolute centroid term of v1 (fought the seams on tall boxes) */
/* The frame direction of a component under the axial prior is its anchor chart's own z-gradient
 * direction (the anchor's rotation is 0 by convention).  When the anchor carries no usable gradient
 * (|grad| outside [AP_ZPRIOR_MIN_GRAD, 1.15]) every prior term of the WHOLE component was silent:
 * measured 2026-09-17 on the 21x5x5, whose 2,337-chart largest component (52% of the area) read
 * "axial plane: v = 0.0 + 0.0000 s" for exactly that reason.  1: such a component takes the
 * area-weighted circular mean of its charts' gradient directions rotated by their current poses
 * (recomputed per GN iteration); a component whose anchor is under the prior is unchanged, so the
 * 4x5x5 is bit-identical.  0 = the anchor only. */
#define AP_ZPRIOR_FRAME 1
/* WEAK ABSOLUTE AXIAL PRIOR (v3).  Per component and GN iteration the plane v = a + b s + c u is
 * fitted over the component's charts (V, U = the chart centroid's frame coordinates along the
 * anchor's gradient direction and its perpendicular, s = the axial coordinate), then every chart
 * is pulled toward it with weight AP_VPRIOR_W = 1/(10 vox)^2: four orders under a seam's w_xy, so
 * a seam never feels it locally; the crossover wavelength (8 pi^2 w_theta d^2 / w)^(1/4) is
 * ~1,900 vox, so the 2,000-vox undulation of the relative prior's random walk (familiar primary
 * p50 20 / p90 41 vox, 4x21x21 60 / 112, 2026-09-09) is what it removes.  The v1 absolute term
 * fought the seams because its target (the anchor's gradient direction, divided by the mean
 * gradient magnitude) was SYSTEMATIC: a residual that does not average out integrates into a
 * distributed load along a chain.  A least-squares plane has a zero-mean residual orthogonal to
 * 1, s and u.  The prior carries no rho and never enters the switch residual.  0 = off (A/B). */
#define AP_VPRIOR 1
#define AP_VPRIOR_W 0.01
#define AP_VPRIOR_MIN_CHARTS 8
#define AP_VPRIOR_MIN_SPREAD 64.0   /* vox of s spread below which b is fixed at 0 (a one-cube-tall strip) */
#define AP_VPRIOR_FIT_U 1           /* the c u term (a slightly rotated frame); 0 = A/B */

void AsmPose_default_opts(AsmPoseOpts *o)
{
    memset(o, 0, sizeof *o);
    o->gn_iters = 12;
    o->cauchy_c = 2.5;
    o->switch_off = ASM_SWITCH_OFF;
    o->axis_dir[0] = 1.0; o->axis_dir[1] = 0.0; o->axis_dir[2] = 0.0;
    o->damping = 1e-6;
}

static double ap_wrap(double a)
{
    while (a > AP_PI) a -= 2.0 * AP_PI;
    while (a < -AP_PI) a += 2.0 * AP_PI;
    return a;
}

/* ---- COO lower-triangle assembly --------------------------------------------- */

typedef struct ApCoo {
    Arena_T arena;
    int *rows, *cols; double *vals;
    size_t n, cap;
} ApCoo;

static void ap_coo_reset(ApCoo *c) { c->n = 0; }

static void ap_coo_add(ApCoo *c, int i, int j, double v)
{
    if (i < j) { int t = i; i = j; j = t; }
    if (c->n == c->cap) {
        size_t ncap = c->cap ? c->cap * 2 : 1 << 16;
        int *nr = ARENA_ALLOC(c->arena, ncap * sizeof(int));
        int *nc = ARENA_ALLOC(c->arena, ncap * sizeof(int));
        double *nv = ARENA_ALLOC(c->arena, ncap * sizeof(double));
        if (c->n) { memcpy(nr, c->rows, c->n * sizeof(int)); memcpy(nc, c->cols, c->n * sizeof(int)); memcpy(nv, c->vals, c->n * sizeof(double)); }
        c->rows = nr; c->cols = nc; c->vals = nv; c->cap = ncap;
    }
    c->rows[c->n] = i; c->cols[c->n] = j; c->vals[c->n] = v; c->n++;
}

static int ap_solve(const ApCoo *c, int n, const double *b, double *x)
{
    int rc = 0;
#ifdef _OPENMP
#pragma omp critical(qrf_taucs_transaction)
#endif
    {
        rc = Sparse_solve_sym(n, (int)c->n, c->rows, c->cols, c->vals, b, x, SPARSE_SPD);
    }
    return rc;
}

/* ---- effective relation under the mirror bits --------------------------------- */

typedef struct ApRel {
    int32_t a, b;        /* active indices */
    double theta, tx, ty;
    double point_a[2], point_b[2]; /* fitted seam centres, relative to chart centroids */
    double w_xy, w_theta;
    double w_xy_s, w_theta_s;   /* the floored weights the switch decision uses */
    double rho;          /* IRLS weight */
    double e;            /* normalized residual */
    double rxy;          /* translation residual, vox */
    size_t rel_index;
} ApRel;

static void ap_effective(const AsmRelation *r, int mirror_a, double *theta, double *tx, double *ty)
{
    if (!mirror_a) { *theta = r->theta; *tx = r->tx; *ty = r->ty; }
    else { *theta = -r->theta; *tx = -r->tx; *ty = r->ty; }
}

/* Solve translations at a material point, not at an arbitrary UV origin.
 * This also makes the damping and warm-start hold independent of that origin.
 * The original chart coordinates and source geometry are never modified. */
static void ap_chart_centres(const AsmRun *run, const int32_t *chart_of, size_t n,
                             double *cu, double *cv)
{
    for (size_t i=0; i<n; i++) {
        const AsmChart *c=&run->charts[chart_of[i]];
        double u=0,v=0;
        if (c->uv && c->nv) {
            for (size_t k=0; k<c->nv; k++) { u+=c->uv[2*k]; v+=c->uv[2*k+1]; }
            u/=(double)c->nv; v/=(double)c->nv;
        }
        cu[i]=(c->flags & ASM_CHART_MIRROR) ? -u : u; cv[i]=v;
    }
}

static void ap_seam_centres(const AsmRun *run, const AsmRelation *r, ApRel *e,
                            const double *cu, const double *cv)
{
    const AsmChart *b=&run->charts[r->b];
    double u=cu[e->b],v=cv[e->b];
    /* Match the centre whose spread supplied w_theta in ar_group_relation:
     * the mean of the two directed, trim-compensated observations. Use the
     * same all-or-nothing transport fallback as the original rigid fit. */
    if (b->uv && r->corr_first>=0 && r->corr_count>0 &&
        (size_t)r->corr_first<=run->n_corr &&
        (size_t)r->corr_count<=run->n_corr-(size_t)r->corr_first) {
        int valid=1,transported=1;
        for (int32_t k=0; k<r->corr_count; k++) {
            const AsmCorr *c=&run->corr[(size_t)r->corr_first+(size_t)k];
            if (c->vb<0 || (size_t)c->vb>=b->nv) valid=0;
            if (c->valid!=3) transported=0;
        }
        if (valid) {
            u=v=0;
            for (int32_t k=0; k<r->corr_count; k++) {
                const AsmCorr *c=&run->corr[(size_t)r->corr_first+(size_t)k];
                u+=b->uv[2*(size_t)c->vb]-(transported ? 0.5*c->gap_b[0] : 0);
                v+=b->uv[2*(size_t)c->vb+1]-(transported ? 0.5*c->gap_b[1] : 0);
            }
            u/=(double)r->corr_count; v/=(double)r->corr_count;
            if (b->flags & ASM_CHART_MIRROR) u=-u;
        }
    }
    double co=cos(e->theta),si=sin(e->theta);
    e->point_a[0]=co*u-si*v+e->tx-cu[e->a];
    e->point_a[1]=si*u+co*v+e->ty-cv[e->a];
    e->point_b[0]=u-cu[e->b]; e->point_b[1]=v-cv[e->b];
}

/* The fit supplies independent translation/rotation information about its
 * seam centre. Measuring translation at B's UV origin adds an arbitrary
 * rotation lever arm and can switch a different edge after a UV recenter. */
static void ap_seam_delta(const ApRel *e, double pa, double pb, double *dx, double *dy,
                           double J[3][6])
{
    double ca=cos(pa),sa=sin(pa),cb=cos(pb),sb=sin(pb);
    double ax=ca*e->point_a[0]-sa*e->point_a[1],ay=sa*e->point_a[0]+ca*e->point_a[1];
    double bx=cb*e->point_b[0]-sb*e->point_b[1],by=sb*e->point_b[0]+cb*e->point_b[1];
    *dx=ax-bx; *dy=ay-by;
    if (J) {
        memset(J,0,18*sizeof(double));
        J[0][0]=-1; J[0][2]=ay; J[0][3]=1; J[0][5]=-by;
        J[1][1]=-1; J[1][2]=-ax; J[1][4]=1; J[1][5]=bx;
        J[2][2]=-1; J[2][5]=1;
    }
}

static int ap_cmp_double(const void *x, const void *y)
{
    double a = *(const double *)x, b = *(const double *)y;
    return (a > b) - (a < b);
}

/* The v3 prior's plane per component: V ~ a + b s + c U over the component's charts, area
 * weighted (V, U = the chart centroid in the anchor's frame direction e and its perpendicular).
 * Returns the number of components with a fit; fits[comp*3..] hold (a, b, c), ok[comp] = 1. */
static size_t ap_vprior_fit(size_t na, const int32_t *comp, size_t ncomp, const int32_t *chart_of, const AsmRun *run,
                            const uint8_t *zg_ok, const double *Vc, const double *Uc, const double *sc,
                            double *fits, uint8_t *ok)
{
    Arena_Mark mark = Arena_save(run->arena);
    double *A = ARENA_CALLOC(run->arena, ncomp * 9, sizeof(double)), *B = ARENA_CALLOC(run->arena, ncomp * 3, sizeof(double));
    double *ws = ARENA_CALLOC(run->arena, ncomp, sizeof(double)), *ss = ARENA_CALLOC(run->arena, ncomp, sizeof(double)), *sss = ARENA_CALLOC(run->arena, ncomp, sizeof(double));
    size_t *cnt = ARENA_CALLOC(run->arena, ncomp, sizeof(size_t));
    for (size_t i = 0; i < na; i++) {
        if (!zg_ok[i]) continue;
        size_t c = (size_t)comp[i];
        double w = run->charts[chart_of[i]].area3d;
        if (w <= 0.0) continue;
        double f[3] = { 1.0, sc[i], Uc[i] };
        for (int r = 0; r < 3; r++) { for (int k = 0; k < 3; k++) A[c*9 + (size_t)(r*3+k)] += w * f[r] * f[k]; B[c*3 + (size_t)r] += w * f[r] * Vc[i]; }
        ws[c] += w; ss[c] += w * sc[i]; sss[c] += w * sc[i] * sc[i]; cnt[c]++;
    }
    size_t nfit = 0;
    for (size_t c = 0; c < ncomp; c++) {
        ok[c] = 0;
        if (cnt[c] < (size_t)AP_VPRIOR_MIN_CHARTS || ws[c] <= 0.0) continue;
        double spread = sqrt(fmax(sss[c] / ws[c] - (ss[c] / ws[c]) * (ss[c] / ws[c]), 0.0));
        int use_b = spread >= AP_VPRIOR_MIN_SPREAD, use_c = AP_VPRIOR_FIT_U;
        double M[9], R[3];
        memcpy(M, &A[c*9], sizeof M); memcpy(R, &B[c*3], sizeof R);
        /* drop a column by pinning its unknown to 0: a unit row */
        if (!use_b) { for (int k = 0; k < 3; k++) { M[1*3+k] = 0.0; M[k*3+1] = 0.0; } M[4] = 1.0; R[1] = 0.0; }
        if (!use_c) { for (int k = 0; k < 3; k++) { M[2*3+k] = 0.0; M[k*3+2] = 0.0; } M[8] = 1.0; R[2] = 0.0; }
        int fine = 1;
        for (int col = 0; col < 3 && fine; col++) {
            int piv = col;
            for (int r = col + 1; r < 3; r++) if (fabs(M[r*3+col]) > fabs(M[piv*3+col])) piv = r;
            for (int k = 0; k < 3; k++) { double t = M[col*3+k]; M[col*3+k] = M[piv*3+k]; M[piv*3+k] = t; }
            { double t = R[col]; R[col] = R[piv]; R[piv] = t; }
            if (fabs(M[col*3+col]) < 1e-12) { fine = 0; break; }
            for (int r = 0; r < 3; r++) {
                if (r == col) continue;
                double f = M[r*3+col] / M[col*3+col];
                for (int k = 0; k < 3; k++) M[r*3+k] -= f * M[col*3+k];
                R[r] -= f * R[col];
            }
        }
        if (!fine) continue;
        fits[c*3] = R[0] / M[0]; fits[c*3+1] = R[1] / M[4]; fits[c*3+2] = R[2] / M[8];
        ok[c] = 1; nfit++;
    }
    Arena_restore(run->arena, mark);
    return nfit;
}

/* |V - plane| percentiles over one component (the largest), for the ledger */
static void ap_vplane_stats(size_t na, const int32_t *comp, int32_t which, const AsmRun *run, const uint8_t *zg_ok,
                            const double *Vc, const double *Uc, const double *sc, const double *fit, double *p50, double *p90)
{
    Arena_Mark mark = Arena_save(run->arena);
    double *r = ARENA_ALLOC(run->arena, (na ? na : 1) * sizeof(double));
    size_t n = 0;
    for (size_t i = 0; i < na; i++) {
        if (!zg_ok[i] || comp[i] != which) continue;
        r[n++] = fabs(Vc[i] - (fit[0] + fit[1] * sc[i] + fit[2] * Uc[i]));
    }
    *p50 = 0.0; *p90 = 0.0;
    if (n > 0) { qsort(r, n, sizeof(double), ap_cmp_double); *p50 = r[n/2]; *p90 = r[(size_t)((double)(n-1) * 0.9)]; }
    Arena_restore(run->arena, mark);
}

/* The frame direction of every component under the axial prior (AP_ZPRIOR_FRAME): the anchor's own
 * z-gradient direction when the anchor is under the prior (cf_ok 1), else the area-weighted circular
 * mean of the component's charts' gradient directions rotated by their current poses (cf_ok 2);
 * cf_ok 0 = no chart under the prior. */
static void ap_component_frames(size_t na, size_t ncomp, const int32_t *comp, const int32_t *anchor, const int32_t *chart_of,
                                const AsmRun *run, const uint8_t *zg_ok, const double *zg_a, const double *zg_b, const double *phi,
                                double *cf_phi, uint8_t *cf_ok, double *sx, double *sy)
{
    for (size_t c = 0; c < ncomp; c++) {
        size_t A = (size_t)anchor[c];
        cf_ok[c] = zg_ok[A] ? 1 : 0; cf_phi[c] = zg_ok[A] ? atan2(zg_b[A], zg_a[A]) : 0.0; sx[c] = sy[c] = 0.0;
    }
    if (!AP_ZPRIOR_FRAME) return;
    for (size_t i = 0; i < na; i++) {
        size_t c = (size_t)comp[i];
        if (!zg_ok[i] || cf_ok[c]) continue;
        double w = run->charts[chart_of[i]].area3d, ang = phi[i] + atan2(zg_b[i], zg_a[i]);
        sx[c] += w * cos(ang); sy[c] += w * sin(ang);
    }
    for (size_t c = 0; c < ncomp; c++)
        if (!cf_ok[c] && (sx[c] != 0.0 || sy[c] != 0.0)) { cf_phi[c] = atan2(sy[c], sx[c]); cf_ok[c] = 2; }
}

int AsmPose_solve(AsmRun *run, const AsmPoseOpts *o, AsmPoseStats *st)
{
    memset(st, 0, sizeof *st);
    double t0 = ves_clock_sec();
    Arena_T arena = run->arena;
    Arena_Mark mark = Arena_save(arena);
    size_t nc = run->n_charts;

    /* active charts */
    int32_t *act = ARENA_ALLOC(arena, nc * sizeof(int32_t));   /* chart -> active index or -1 */
    size_t na = 0;
    for (size_t i = 0; i < nc; i++) {
        act[i] = AsmChart_in_layout(&run->charts[i]) ? (int32_t)na++ : -1;
        run->charts[i].component = -1;
        run->charts[i].placed = 0;
        run->charts[i].flags &= ~(uint32_t)ASM_CHART_MIRROR;
    }
    st->charts_active = na;
    if (na == 0) { Arena_restore(arena, mark); return 0; }
    int32_t *chart_of = ARENA_ALLOC(arena, na * sizeof(int32_t));
    for (size_t i = 0; i < nc; i++) if (act[i] >= 0) chart_of[act[i]] = (int32_t)i;

    /* active relations (parity resolved below) */
    size_t nr_all = run->n_rels;
    int32_t *rel_ok = ARENA_CALLOC(arena, nr_all ? nr_all : 1, sizeof(int32_t));
    for (size_t r = 0; r < nr_all; r++) {
        AsmRelation *R = &run->rels[r];
        R->flags &= ~(uint32_t)ASM_REL_SWITCHED;
        if (R->flags & (ASM_REL_DROPPED | ASM_REL_CONTACT | ASM_REL_PLACEMENT_ONLY)) continue;
        if (act[R->a] < 0 || act[R->b] < 0) continue;
        rel_ok[r] = 1;
    }

    /* adjacency over active relations (CSR) */
    int32_t *deg = ARENA_CALLOC(arena, na + 1, sizeof(int32_t));
    for (size_t r = 0; r < nr_all; r++) if (rel_ok[r]) { deg[act[run->rels[r].a] + 1]++; deg[act[run->rels[r].b] + 1]++; }
    for (size_t i = 0; i < na; i++) deg[i+1] += deg[i];
    int32_t *adj = ARENA_ALLOC(arena, ((size_t)deg[na] + 1) * sizeof(int32_t));   /* relation index */
    int32_t *fill = ARENA_CALLOC(arena, na, sizeof(int32_t));
    for (size_t r = 0; r < nr_all; r++) {
        if (!rel_ok[r]) continue;
        int32_t ia = act[run->rels[r].a], ib = act[run->rels[r].b];
        adj[deg[ia] + fill[ia]++] = (int32_t)r;
        adj[deg[ib] + fill[ib]++] = (int32_t)r;
    }

    /* parity 2-colouring by BFS; conflicts drop the closing relation */
    int8_t *mirror = ARENA_ALLOC(arena, na);
    for (size_t i = 0; i < na; i++) mirror[i] = -1;
    int32_t *queue = ARENA_ALLOC(arena, na * sizeof(int32_t));
    for (size_t s = 0; s < na; s++) {
        if (mirror[s] >= 0) continue;
        mirror[s] = 0;
        size_t qh = 0, qt = 0;
        queue[qt++] = (int32_t)s;
        while (qh < qt) {
            int32_t u = queue[qh++];
            for (int32_t k = deg[u]; k < deg[u+1]; k++) {
                int32_t r = adj[k];
                if (!rel_ok[r]) continue;
                const AsmRelation *R = &run->rels[r];
                int32_t v = (act[R->a] == u) ? act[R->b] : act[R->a];
                int8_t want = (int8_t)(mirror[u] ^ ((R->flags & ASM_REL_PARITY) ? 1 : 0));
                if (mirror[v] < 0) { mirror[v] = want; queue[qt++] = v; }
                else if (mirror[v] != want) {
                    rel_ok[r] = 0;
                    run->rels[r].flags |= ASM_REL_DROPPED;
                    st->parity_conflicts++;
                }
            }
        }
    }
    for (size_t i = 0; i < na; i++) if (mirror[i]) run->charts[chart_of[i]].flags |= ASM_CHART_MIRROR;
    double *centre_u=ARENA_ALLOC(arena,na*sizeof(double)), *centre_v=ARENA_ALLOC(arena,na*sizeof(double));
    ap_chart_centres(run,chart_of,na,centre_u,centre_v);

    /* components + anchors */
    UnionFind uf = UF_new(arena, (int32_t)na);
    for (size_t r = 0; r < nr_all; r++) if (rel_ok[r]) uf_union(&uf, act[run->rels[r].a], act[run->rels[r].b]);
    int32_t *comp = ARENA_ALLOC(arena, na * sizeof(int32_t));
    int32_t *root_comp = ARENA_ALLOC(arena, na * sizeof(int32_t));
    for (size_t i = 0; i < na; i++) root_comp[i] = -1;
    size_t ncomp = 0;
    for (size_t i = 0; i < na; i++) {
        int32_t rt = uf_find(&uf, (int32_t)i);
        if (root_comp[rt] < 0) root_comp[rt] = (int32_t)ncomp++;
        comp[i] = root_comp[rt];
    }
    int32_t *anchor = ARENA_ALLOC(arena, ncomp * sizeof(int32_t));
    double *comp_area = ARENA_CALLOC(arena, ncomp, sizeof(double));
    size_t *comp_n = ARENA_CALLOC(arena, ncomp, sizeof(size_t));
    for (size_t c = 0; c < ncomp; c++) anchor[c] = -1;
    for (size_t i = 0; i < na; i++) {
        const AsmChart *ch = &run->charts[chart_of[i]];
        size_t c = (size_t)comp[i];
        comp_area[c] += ch->area3d; comp_n[c]++;
        if (anchor[c] < 0 || ch->area3d > run->charts[chart_of[anchor[c]]].area3d) anchor[c] = (int32_t)i;
    }
    st->components = ncomp;
    for (size_t i = 0; i < na; i++) {
        int d = deg[i+1] - deg[i];
        if (d > 5) d = 5;
        st->degree_hist[d]++;
        st->degree_area[d] += run->charts[chart_of[i]].area3d;
    }
    for (size_t c = 0; c < ncomp; c++)
        if (comp_n[c] > st->largest_component_charts) { st->largest_component_charts = comp_n[c]; st->largest_component_area = comp_area[c]; }
    uint8_t *is_anchor = ARENA_CALLOC(arena, na, sizeof(uint8_t));
    for (size_t c = 0; c < ncomp; c++) is_anchor[anchor[c]] = 1;

    /* per chart: z ~ a u + b v + c over its vertices (u mirrored when the parity says so) */
    double *zg_a = ARENA_CALLOC(arena, na, sizeof(double)), *zg_b = ARENA_CALLOC(arena, na, sizeof(double));
    double *zg_z = ARENA_CALLOC(arena, na, sizeof(double)), *zg_cu = ARENA_CALLOC(arena, na, sizeof(double)), *zg_cv = ARENA_CALLOC(arena, na, sizeof(double));
    uint8_t *zg_ok = ARENA_CALLOC(arena, na, sizeof(uint8_t));
    size_t n_prior = 0;
    for (size_t i = 0; i < na; i++) {
        const AsmChart *ch = &run->charts[chart_of[i]];
        if (ch->nv < 30 || ch->uv == NULL) continue;
        size_t stride = ch->nv / 400; if (stride < 1) stride = 1;
        /* the chart's axial coordinate: the table's arc length at its centroid's foot plus the
         * local tangent's offset (exact to first order within one chart); world z without a table */
        double ax_s = 0.0, ax_dir[3] = { o->axis_dir[0], o->axis_dir[1], o->axis_dir[2] }, ax_foot[3] = { 0.0, 0.0, 0.0 };
        if (o->axis) { double cen[3] = { ch->centroid[0], ch->centroid[1], ch->centroid[2] }; AsmAxis_project(o->axis, cen, &ax_s, NULL, ax_dir, ax_foot); }
        double su = 0, sv = 0, sz = 0, n = 0;
        for (size_t v = 0; v < ch->nv; v += stride) {
            double uu = mirror[i] ? -ch->uv[v*2] : ch->uv[v*2], vv = ch->uv[v*2+1];
            double zz = ax_s + ax_dir[0]*(ch->xyz[v*3]-ax_foot[0]) + ax_dir[1]*(ch->xyz[v*3+1]-ax_foot[1]) + ax_dir[2]*(ch->xyz[v*3+2]-ax_foot[2]);
            su += uu; sv += vv; sz += zz; n += 1.0;
        }
        su /= n; sv /= n; sz /= n;
        double suu = 0, svv = 0, suv = 0, suz = 0, svz = 0;
        for (size_t v = 0; v < ch->nv; v += stride) {
            double uu = (mirror[i] ? -ch->uv[v*2] : ch->uv[v*2]) - su, vv = ch->uv[v*2+1] - sv;
            double zz = ax_s + ax_dir[0]*(ch->xyz[v*3]-ax_foot[0]) + ax_dir[1]*(ch->xyz[v*3+1]-ax_foot[1]) + ax_dir[2]*(ch->xyz[v*3+2]-ax_foot[2]) - sz;
            suu += uu*uu; svv += vv*vv; suv += uu*vv; suz += uu*zz; svz += vv*zz;
        }
        double det = suu*svv - suv*suv;
        if (fabs(det) < 1e-6) continue;
        double a = (suz*svv - svz*suv) / det, b = (svz*suu - suz*suv) / det;
        double gmag = hypot(a, b);
        if (gmag < AP_ZPRIOR_MIN_GRAD || gmag > 1.15) continue;
        zg_a[i] = a; zg_b[i] = b; zg_z[i] = sz;
        zg_cu[i] = su-centre_u[i]; zg_cv[i] = sv-centre_v[i]; zg_ok[i] = 1; n_prior++;
    }
    st->z_prior_charts = n_prior;
    if (o->verbose_prior) fprintf(stderr, "  [pose] axial prior on %zu of %zu charts (gradient magnitude >= %.2f, axis (%.3f, %.3f, %.3f) zyx)\n", n_prior, na, AP_ZPRIOR_MIN_GRAD, o->axis_dir[0], o->axis_dir[1], o->axis_dir[2]);

    /* effective relations */
    size_t nrel = 0;
    for (size_t r = 0; r < nr_all; r++) nrel += rel_ok[r];
    ApRel *rel = ARENA_ALLOC(arena, (nrel ? nrel : 1) * sizeof(ApRel));
    {
        size_t q = 0;
        for (size_t r = 0; r < nr_all; r++) {
            if (!rel_ok[r]) continue;
            const AsmRelation *R = &run->rels[r];
            ApRel *e = &rel[q++];
            e->a = act[R->a]; e->b = act[R->b];
            ap_effective(R, mirror[e->a], &e->theta, &e->tx, &e->ty);
            ap_seam_centres(run,R,e,centre_u,centre_v);
            e->w_xy = R->w_xy > 0 ? R->w_xy : 1.0;
            e->w_theta = R->w_theta > 1e-6 ? R->w_theta : 1e-6;
            e->w_xy_s = e->w_xy; e->w_theta_s = e->w_theta;
            if (AP_SEAM_FLOOR_VOX > 0.0) {
                double f = 1.0 / (1.0 + AP_SEAM_FLOOR_VOX * AP_SEAM_FLOOR_VOX * e->w_xy);
                e->w_xy_s *= f; e->w_theta_s *= f;
                if (AP_SEAM_FLOOR_IN_SOLVE) { e->w_xy *= f; e->w_theta *= f; }
            }
            e->rho = 1.0; e->e = 0.0; e->rxy = 0.0; e->rel_index = r;
        }
    }
    st->rels_active = nrel;

    double *phi = ARENA_CALLOC(arena, na, sizeof(double));
    double *px = ARENA_CALLOC(arena, na, sizeof(double));
    double *py = ARENA_CALLOC(arena, na, sizeof(double));
    ApCoo coo; memset(&coo, 0, sizeof coo); coo.arena = arena;

    if (nrel > 0) {
        /* 1. chordal rotations: z_b = e^{i theta} z_a, anchors z = 1 */
        int n2 = (int)(2 * na);
        double *b2 = ARENA_CALLOC(arena, (size_t)n2, sizeof(double));
        double *x2 = ARENA_CALLOC(arena, (size_t)n2, sizeof(double));
        ap_coo_reset(&coo);
        for (size_t q = 0; q < nrel; q++) {
            const ApRel *e = &rel[q];
            double w = e->w_theta > 1.0 ? 1.0 + log(e->w_theta) : e->w_theta; /* tame the dynamic range */
            double c = cos(e->theta), s = sin(e->theta);
            /* residual r = [pb - (c pa - s qa); qb - (s pa + c qa)]; J columns: pa, qa, pb, qb */
            double J[2][4] = { { -c, s, 1, 0 }, { -s, -c, 0, 1 } };
            int idx[4] = { 2*e->a, 2*e->a + 1, 2*e->b, 2*e->b + 1 };
            for (int r = 0; r < 4; r++) for (int cc = 0; cc <= r; cc++) {
                double v = w * (J[0][r]*J[0][cc] + J[1][r]*J[1][cc]);
                if (v != 0.0) ap_coo_add(&coo, idx[r], idx[cc], v);
            }
        }
        /* anchors: fix z = (1, 0) by identity rows, moving columns to the rhs */
        for (size_t i = 0; i < na; i++) if (is_anchor[i]) { x2[2*i] = 1.0; x2[2*i+1] = 0.0; }
        for (size_t k = 0; k < coo.n; k++) {
            int r = coo.rows[k], c = coo.cols[k];
            int ra = is_anchor[r/2], ca = is_anchor[c/2];
            if (ra && ca) { coo.vals[k] = 0.0; continue; }
            if (ra) { b2[c] -= coo.vals[k] * x2[r]; coo.vals[k] = 0.0; continue; }
            if (ca) { b2[r] -= coo.vals[k] * x2[c]; coo.vals[k] = 0.0; continue; }
        }
        for (size_t i = 0; i < na; i++) {
            if (is_anchor[i]) { b2[2*i] = 1.0; b2[2*i+1] = 0.0; ap_coo_add(&coo, 2*(int)i, 2*(int)i, 1.0); ap_coo_add(&coo, 2*(int)i+1, 2*(int)i+1, 1.0); }
            else { ap_coo_add(&coo, 2*(int)i, 2*(int)i, 1e-9); ap_coo_add(&coo, 2*(int)i+1, 2*(int)i+1, 1e-9); }
        }
        /* isolated non-anchor charts (no relations) get identity rows */
        for (size_t i = 0; i < na; i++) if (!is_anchor[i] && deg[i+1] == deg[i]) { ap_coo_add(&coo, 2*(int)i, 2*(int)i, 1.0); ap_coo_add(&coo, 2*(int)i+1, 2*(int)i+1, 1.0); }
        if (ap_solve(&coo, n2, b2, x2) != 0) { fprintf(stderr, "  [pose] chordal solve failed\n"); Arena_restore(arena, mark); return -1; }
        for (size_t i = 0; i < na; i++) {
            double re = x2[2*i], im = x2[2*i+1];
            phi[i] = (re*re + im*im > 1e-20) ? atan2(im, re) : 0.0;
        }

        /* 2. translations at fixed rotations (x and y decouple) */
        int n1 = (int)na;
        double *bx = ARENA_CALLOC(arena, na, sizeof(double));
        double *by = ARENA_CALLOC(arena, na, sizeof(double));
        ap_coo_reset(&coo);
        for (size_t q = 0; q < nrel; q++) {
            const ApRel *e = &rel[q];
            double w = e->w_xy;
            double dx,dy;
            ap_seam_delta(e,phi[e->a],phi[e->b],&dx,&dy,NULL);
            ap_coo_add(&coo, e->a, e->a, w); ap_coo_add(&coo, e->b, e->b, w); ap_coo_add(&coo, e->b, e->a, -w);
            bx[e->a] -= w * dx; bx[e->b] += w * dx;
            by[e->a] -= w * dy; by[e->b] += w * dy;
        }
        for (size_t i = 0; i < na; i++) {
            if (is_anchor[i] || deg[i+1] == deg[i]) ap_coo_add(&coo, (int)i, (int)i, 1e6);
            else ap_coo_add(&coo, (int)i, (int)i, 1e-9);
        }
        if (ap_solve(&coo, n1, bx, px) != 0 || ap_solve(&coo, n1, by, py) != 0) { fprintf(stderr, "  [pose] translation solve failed\n"); Arena_restore(arena, mark); return -1; }
        for (size_t i = 0; i < na; i++) if (is_anchor[i]) { px[i] = 0.0; py[i] = 0.0; }

        /* 3. Gauss-Newton + IRLS on (x, y, phi) */
        int n3 = (int)(3 * na);
        double *g = ARENA_CALLOC(arena, (size_t)n3, sizeof(double));
        double *dx = ARENA_CALLOC(arena, (size_t)n3, sizeof(double));
        int it;
        /* v3 prior scratch: per chart the frame coordinates of its centroid, per component the plane */
        double *vp_V = ARENA_CALLOC(arena, na, sizeof(double)), *vp_U = ARENA_CALLOC(arena, na, sizeof(double));
        double *vp_fit = ARENA_CALLOC(arena, ncomp * 3, sizeof(double));
        uint8_t *vp_ok = ARENA_CALLOC(arena, ncomp, sizeof(uint8_t));
        int32_t vp_largest = 0;
        for (size_t c = 1; c < ncomp; c++) if (comp_n[c] > comp_n[(size_t)vp_largest]) vp_largest = (int32_t)c;
        /* the frame direction of every component (AP_ZPRIOR_FRAME), refreshed per iteration */
        double *cf_phi = ARENA_CALLOC(arena, ncomp, sizeof(double)), *cf_sx = ARENA_CALLOC(arena, ncomp, sizeof(double)), *cf_sy = ARENA_CALLOC(arena, ncomp, sizeof(double));
        uint8_t *cf_ok = ARENA_CALLOC(arena, ncomp, sizeof(uint8_t));
        for (it = 0; it < o->gn_iters; it++) {
            ap_component_frames(na, ncomp, comp, anchor, chart_of, run, zg_ok, zg_a, zg_b, phi, cf_phi, cf_ok, cf_sx, cf_sy);
            if (it == 0) {
                for (size_t c = 0; c < ncomp; c++) { if (cf_ok[c] == 1) st->frames_from_anchor++; else if (cf_ok[c] == 2) st->frames_from_mean++; else st->frames_none++; }
                st->largest_frame = cf_ok[(size_t)vp_largest];
                if (o->verbose_prior) {
                    /* the three largest components: where their frame direction comes from */
                    int32_t shown[3] = { -1, -1, -1 };
                    for (int rank = 0; rank < 3; rank++) {
                        int32_t best = -1;
                        for (size_t c = 0; c < ncomp; c++) {
                            int taken = 0;
                            for (int r2 = 0; r2 < rank; r2++) if (shown[r2] == (int32_t)c) taken = 1;
                            if (!taken && (best < 0 || comp_n[c] > comp_n[(size_t)best])) best = (int32_t)c;
                        }
                        if (best < 0) break;
                        shown[rank] = best;
                        size_t under = 0, A = (size_t)anchor[best];
                        for (size_t i = 0; i < na; i++) if (comp[i] == best && zg_ok[i]) under++;
                        fprintf(stderr, "  [pose] component %d: %zu charts, anchor chart %d under the prior: %s; %zu charts under the prior; frame direction from %s\n",
                                best, comp_n[(size_t)best], chart_of[A], zg_ok[A] ? "yes" : "NO (gradient outside the band)", under,
                                cf_ok[(size_t)best] == 1 ? "the anchor" : cf_ok[(size_t)best] == 2 ? "the charts' mean" : "NOWHERE (the prior is silent here)");
                    }
                }
            }
            if (AP_VPRIOR && AP_ZPRIOR) {
                /* the chart centroids in the frame: V along the component's frame direction, U across it */
                for (size_t i = 0; i < na; i++) {
                    vp_V[i] = 0.0; vp_U[i] = 0.0;
                    if (!zg_ok[i] || !cf_ok[comp[i]]) continue;
                    double Phi = cf_phi[comp[i]], ex = cos(Phi), ey = sin(Phi);
                    double c = cos(phi[i]), sn = sin(phi[i]);
                    double Px = px[i] + c*zg_cu[i] - sn*zg_cv[i], Py = py[i] + sn*zg_cu[i] + c*zg_cv[i];
                    vp_V[i] = Px * ex + Py * ey; vp_U[i] = -Px * ey + Py * ex;
                }
                ap_vprior_fit(na, comp, ncomp, chart_of, run, zg_ok, vp_V, vp_U, zg_z, vp_fit, vp_ok);
                if (it == 0 && vp_ok[(size_t)vp_largest])
                    ap_vplane_stats(na, comp, vp_largest, run, zg_ok, vp_V, vp_U, zg_z, &vp_fit[(size_t)vp_largest*3], &st->vplane_p50_before, &st->vplane_p90_before);
            }
            /* residuals + robust weights */
            for (size_t q = 0; q < nrel; q++) {
                ApRel *e = &rel[q];
                double sx,sy; ap_seam_delta(e,phi[e->a],phi[e->b],&sx,&sy,NULL);
                double rx = px[e->b] - px[e->a] - sx;
                double ry = py[e->b] - py[e->a] - sy;
                double rt = ap_wrap(phi[e->b] - phi[e->a] - e->theta);
                double e2 = (e->w_xy_s * (rx*rx + ry*ry) + e->w_theta_s * rt*rt) / 3.0;
                e->e = sqrt(e2);
                e->rxy = sqrt(rx*rx + ry*ry);
                double u = e->e / o->cauchy_c;
                e->rho = 1.0 / (1.0 + u*u);
            }
            ap_coo_reset(&coo);
            memset(g, 0, (size_t)n3 * sizeof(double));
            for (size_t q = 0; q < nrel; q++) {
                const ApRel *e = &rel[q];
                double sx,sy,J[3][6]; ap_seam_delta(e,phi[e->a],phi[e->b],&sx,&sy,J);
                double rx = px[e->b] - px[e->a] - sx;
                double ry = py[e->b] - py[e->a] - sy;
                double rt = ap_wrap(phi[e->b] - phi[e->a] - e->theta);
                double W[3] = { e->rho * e->w_xy, e->rho * e->w_xy, e->rho * e->w_theta };
                double res[3] = { rx, ry, rt };
                int idx[6] = { 3*e->a, 3*e->a+1, 3*e->a+2, 3*e->b, 3*e->b+1, 3*e->b+2 };
                for (int r = 0; r < 6; r++) {
                    double gi = 0.0;
                    for (int k = 0; k < 3; k++) gi += J[k][r] * W[k] * res[k];
                    g[idx[r]] += gi;
                    for (int cc = 0; cc <= r; cc++) {
                        double v = 0.0;
                        for (int k = 0; k < 3; k++) v += J[k][r] * W[k] * J[k][cc];
                        if (v != 0.0) ap_coo_add(&coo, idx[r], idx[cc], v);
                    }
                }
            }
            if (AP_ZPRIOR) {
                for (size_t i = 0; i < na; i++) {
                    if (!zg_ok[i] || is_anchor[i] || deg[i+1] == deg[i]) continue;
                    size_t A = (size_t)anchor[comp[i]];
                    if (!cf_ok[comp[i]]) continue;
                    /* the frame's z direction: the anchor's own (its rotation is 0 by convention), or the component's mean */
                    double Phi = cf_phi[comp[i]], ex = cos(Phi), ey = sin(Phi);
                    double phz = atan2(zg_b[i], zg_a[i]);
                    /* rotation prior: this chart's z direction, rotated by its pose, must be the frame's */
                    double rt = ap_wrap(phi[i] + phz - Phi);
                    g[3*i+2] += AP_ZPRIOR_W_THETA * rt;
                    ap_coo_add(&coo, 3*(int)i+2, 3*(int)i+2, AP_ZPRIOR_W_THETA);
                    if (AP_ZPRIOR_ABS && zg_ok[A]) {
                        /* position prior along that direction: the centroid's z difference to the anchor */
                        double c = cos(phi[i]), sn = sin(phi[i]);
                        double Px = px[i] + c*zg_cu[i] - sn*zg_cv[i], Py = py[i] + sn*zg_cu[i] + c*zg_cv[i];
                        double PAx = zg_cu[A], PAy = zg_cv[A];
                        double gbar = 0.5 * (hypot(zg_a[i], zg_b[i]) + hypot(zg_a[A], zg_b[A]));
                        double rz = (Px - PAx) * ex + (Py - PAy) * ey - (zg_z[i] - zg_z[A]) / gbar;
                        double Jx = ex, Jy = ey, Jp = ex * (-sn*zg_cu[i] - c*zg_cv[i]) + ey * (c*zg_cu[i] - sn*zg_cv[i]);
                        double J3[3] = { Jx, Jy, Jp };
                        for (int r = 0; r < 3; r++) {
                            g[3*i + (size_t)r] += J3[r] * AP_ZPRIOR_W_Z * rz;
                            for (int cc = 0; cc <= r; cc++) { double v = J3[r] * AP_ZPRIOR_W_Z * J3[cc]; if (v != 0.0) ap_coo_add(&coo, 3*(int)i + r, 3*(int)i + cc, v); }
                        }
                    }
                }
                /* relative position term per seam relation: (P_b - P_a) . e = (z_b - z_a) / gbar */
                for (size_t q = 0; q < nrel; q++) {
                    const ApRel *e = &rel[q];
                    size_t ia = (size_t)e->a, ib = (size_t)e->b;
                    if (!zg_ok[ia] || !zg_ok[ib] || !cf_ok[comp[ia]]) continue;
                    double Phi = cf_phi[comp[ia]], ex = cos(Phi), ey = sin(Phi);
                    double ca = cos(phi[ia]), sa = sin(phi[ia]), cb = cos(phi[ib]), sb = sin(phi[ib]);
                    double Pax = px[ia] + ca*zg_cu[ia] - sa*zg_cv[ia], Pay = py[ia] + sa*zg_cu[ia] + ca*zg_cv[ia];
                    double Pbx = px[ib] + cb*zg_cu[ib] - sb*zg_cv[ib], Pby = py[ib] + sb*zg_cu[ib] + cb*zg_cv[ib];
                    double gbar = 0.5 * (hypot(zg_a[ia], zg_b[ia]) + hypot(zg_a[ib], zg_b[ib]));
                    double rz = (Pbx - Pax) * ex + (Pby - Pay) * ey - (zg_z[ib] - zg_z[ia]) / gbar;
                    double W = AP_ZPRIOR_W_Z * e->rho;
                    /* J columns: xa ya pa xb yb pb */
                    double J6[6] = { -ex, -ey, -(ex * (-sa*zg_cu[ia] - ca*zg_cv[ia]) + ey * (ca*zg_cu[ia] - sa*zg_cv[ia])),
                                      ex,  ey,  (ex * (-sb*zg_cu[ib] - cb*zg_cv[ib]) + ey * (cb*zg_cu[ib] - sb*zg_cv[ib])) };
                    int idx[6] = { 3*e->a, 3*e->a+1, 3*e->a+2, 3*e->b, 3*e->b+1, 3*e->b+2 };
                    for (int r = 0; r < 6; r++) {
                        g[idx[r]] += J6[r] * W * rz;
                        for (int cc = 0; cc <= r; cc++) { double v = J6[r] * W * J6[cc]; if (v != 0.0) ap_coo_add(&coo, idx[r], idx[cc], v); }
                    }
                }
                /* v3: the weak absolute plane term per chart, V_i - (a + b s_i + c U_i) -> 0 */
                if (AP_VPRIOR) {
                    size_t nvp = 0;
                    for (size_t i = 0; i < na; i++) {
                        if (!zg_ok[i] || is_anchor[i] || deg[i+1] == deg[i]) continue;
                        if (!cf_ok[comp[i]] || !vp_ok[(size_t)comp[i]]) continue;
                        const double *f = &vp_fit[(size_t)comp[i] * 3];
                        double Phi = cf_phi[comp[i]], ex = cos(Phi), ey = sin(Phi);
                        double c = cos(phi[i]), sn = sin(phi[i]);
                        /* the plane's residual: V - c U - (a + b s); its direction d = e - c e_perp with e_perp = (-ey, ex) */
                        double dxp = ex - f[2] * (-ey), dyp = ey - f[2] * ex;
                        double rv = vp_V[i] - f[2] * vp_U[i] - (f[0] + f[1] * zg_z[i]);
                        double Jp = dxp * (-sn*zg_cu[i] - c*zg_cv[i]) + dyp * (c*zg_cu[i] - sn*zg_cv[i]);
                        double J3[3] = { dxp, dyp, Jp };
                        for (int r = 0; r < 3; r++) {
                            g[3*i + (size_t)r] += J3[r] * AP_VPRIOR_W * rv;
                            for (int cc = 0; cc <= r; cc++) { double v = J3[r] * AP_VPRIOR_W * J3[cc]; if (v != 0.0) ap_coo_add(&coo, 3*(int)i + r, 3*(int)i + cc, v); }
                        }
                        nvp++;
                    }
                    if (it == 0) st->vprior_charts = nvp;
                }
            }
            for (size_t i = 0; i < na; i++) {
                double d = (is_anchor[i] || deg[i+1] == deg[i]) ? 1e6 : o->damping;
                for (int k = 0; k < 3; k++) ap_coo_add(&coo, 3*(int)i + k, 3*(int)i + k, d);
                if (is_anchor[i] || deg[i+1] == deg[i]) for (int k = 0; k < 3; k++) g[3*i + (size_t)k] = 0.0;
            }
            for (int k = 0; k < n3; k++) g[k] = -g[k];
            memset(dx, 0, (size_t)n3 * sizeof(double));
            if (ap_solve(&coo, n3, g, dx) != 0) { fprintf(stderr, "  [pose] GN solve failed at iter %d\n", it); break; }
            double maxstep = 0.0;
            for (size_t i = 0; i < na; i++) {
                if (is_anchor[i]) continue;
                px[i] += dx[3*i]; py[i] += dx[3*i+1]; phi[i] = ap_wrap(phi[i] + dx[3*i+2]);
                double m = fabs(dx[3*i]) + fabs(dx[3*i+1]) + 100.0 * fabs(dx[3*i+2]);
                if (m > maxstep) maxstep = m;
            }
            if (maxstep < 1e-4) { it++; break; }
        }
        st->iters = it;
        if (AP_VPRIOR && AP_ZPRIOR) {
            ap_component_frames(na, ncomp, comp, anchor, chart_of, run, zg_ok, zg_a, zg_b, phi, cf_phi, cf_ok, cf_sx, cf_sy);
            for (size_t i = 0; i < na; i++) {
                vp_V[i] = 0.0; vp_U[i] = 0.0;
                if (!zg_ok[i] || !cf_ok[comp[i]]) continue;
                double Phi = cf_phi[comp[i]], ex = cos(Phi), ey = sin(Phi);
                double c = cos(phi[i]), sn = sin(phi[i]);
                double Px = px[i] + c*zg_cu[i] - sn*zg_cv[i], Py = py[i] + sn*zg_cu[i] + c*zg_cv[i];
                vp_V[i] = Px * ex + Py * ey; vp_U[i] = -Px * ey + Py * ex;
            }
            ap_vprior_fit(na, comp, ncomp, chart_of, run, zg_ok, vp_V, vp_U, zg_z, vp_fit, vp_ok);
            if (vp_ok[(size_t)vp_largest]) {
                ap_vplane_stats(na, comp, vp_largest, run, zg_ok, vp_V, vp_U, zg_z, &vp_fit[(size_t)vp_largest*3], &st->vplane_p50_after, &st->vplane_p90_after);
                for (int k = 0; k < 3; k++) st->vprior_fit[k] = vp_fit[(size_t)vp_largest*3 + (size_t)k];
            }
            if (o->verbose_prior)
                fprintf(stderr, "  [pose] axial plane prior on %zu charts: largest component v = %.1f + %.4f s + %.4f u; |v - plane| p50 %.1f p90 %.1f -> %.1f %.1f vox\n",
                        st->vprior_charts, st->vprior_fit[0], st->vprior_fit[1], st->vprior_fit[2], st->vplane_p50_before, st->vplane_p90_before, st->vplane_p50_after, st->vplane_p90_after);
        }
        /* final residuals */
        double *elist = ARENA_ALLOC(arena, nrel * sizeof(double));
        double *xlist = ARENA_ALLOC(arena, nrel * sizeof(double));
        for (size_t q = 0; q < nrel; q++) {
            ApRel *e = &rel[q];
            double sx,sy; ap_seam_delta(e,phi[e->a],phi[e->b],&sx,&sy,NULL);
            double rx = px[e->b] - px[e->a] - sx;
            double ry = py[e->b] - py[e->a] - sy;
            double rt = ap_wrap(phi[e->b] - phi[e->a] - e->theta);
            e->e = sqrt((e->w_xy_s * (rx*rx + ry*ry) + e->w_theta_s * rt*rt) / 3.0);
            e->rxy = sqrt(rx*rx + ry*ry);
            double u = e->e / o->cauchy_c;
            e->rho = 1.0 / (1.0 + u*u);
            AsmRelation *R = &run->rels[e->rel_index];
            R->robust_w = e->rho; R->residual = e->e;
            if (e->rho < o->switch_off) { R->flags |= ASM_REL_SWITCHED; st->switched_off++; }
            elist[q] = e->e; xlist[q] = e->rxy;
        }
        /* components AFTER switching: charts held together only by switched-off relations are
         * separate components for everything downstream (reviewer, 2026-09-09) */
        if (AP_COMPONENTS_AFTER_SWITCH) {
            UnionFind uf2 = UF_new(arena, (int32_t)na);
            for (size_t q = 0; q < nrel; q++) { const ApRel *e = &rel[q]; if (e->rho >= o->switch_off) uf_union(&uf2, e->a, e->b); }
            for (size_t i = 0; i < na; i++) root_comp[i] = -1;
            size_t ncomp2 = 0;
            for (size_t i = 0; i < na; i++) { int32_t rt = uf_find(&uf2, (int32_t)i); if (root_comp[rt] < 0) root_comp[rt] = (int32_t)ncomp2++; comp[i] = root_comp[rt]; }
            st->components_after_switch = ncomp2;
        } else st->components_after_switch = ncomp;
        qsort(elist, nrel, sizeof(double), ap_cmp_double);
        qsort(xlist, nrel, sizeof(double), ap_cmp_double);
        st->resid_p50 = elist[nrel/2]; st->resid_p95 = elist[(size_t)((double)(nrel-1)*0.95)]; st->resid_max = elist[nrel-1];
        st->xy_p50 = xlist[nrel/2]; st->xy_p95 = xlist[(size_t)((double)(nrel-1)*0.95)];
    }

    /* write poses */
    for (size_t i = 0; i < na; i++) {
        AsmChart *ch = &run->charts[chart_of[i]];
        double co=cos(phi[i]),si=sin(phi[i]);
        ch->pose_x = px[i]-co*centre_u[i]+si*centre_v[i];
        ch->pose_y = py[i]-si*centre_u[i]-co*centre_v[i]; ch->pose_theta = phi[i];
        ch->component = comp[i];
        ch->placed = 1;
    }
    st->sec = ves_clock_sec() - t0;
    Arena_restore(arena, mark);
    return 0;
}

/* ---- selftest ------------------------------------------------------------------ */

int AsmPose_refine(AsmRun *run, const AsmPoseOpts *o, const int32_t *rel_idx, size_t n_rel, const uint8_t *fixed, AsmPoseRefineStats *st)
{
    memset(st, 0, sizeof *st);
    double t0 = ves_clock_sec();
    Arena_T arena = run->arena;
    Arena_Mark mark = Arena_save(arena);
    size_t nc = run->n_charts;
    int32_t *act = ARENA_ALLOC(arena, (nc ? nc : 1) * sizeof(int32_t));
    for (size_t i = 0; i < nc; i++) act[i] = -1;
    /* the listed relations whose parity agrees with the charts' mirror bits; their endpoints are the active charts */
    uint8_t *use = ARENA_CALLOC(arena, n_rel ? n_rel : 1, 1);
    size_t na = 0, nrel = 0;
    for (size_t q = 0; q < n_rel; q++) {
        const AsmRelation *R = &run->rels[(size_t)rel_idx[q]];
        const AsmChart *A = &run->charts[(size_t)R->a], *B = &run->charts[(size_t)R->b];
        int ma = (A->flags & ASM_CHART_MIRROR) != 0, mb = (B->flags & ASM_CHART_MIRROR) != 0;
        if ((ma ^ mb) != ((R->flags & ASM_REL_PARITY) != 0)) { st->rels_parity_refused++; continue; }
        use[q] = 1; nrel++;
        if (act[R->a] < 0) act[R->a] = (int32_t)na++;
        if (act[R->b] < 0) act[R->b] = (int32_t)na++;
    }
    st->charts_active = na; st->rels = nrel;
    if (na == 0 || nrel == 0) { Arena_restore(arena, mark); return 0; }
    int32_t *chart_of = ARENA_ALLOC(arena, na * sizeof(int32_t));
    for (size_t i = 0; i < nc; i++) if (act[i] >= 0) chart_of[act[i]] = (int32_t)i;
    double *centre_u=ARENA_ALLOC(arena,na*sizeof(double)), *centre_v=ARENA_ALLOC(arena,na*sizeof(double));
    ap_chart_centres(run,chart_of,na,centre_u,centre_v);
    ApRel *rel = ARENA_ALLOC(arena, nrel * sizeof(ApRel));
    {
        size_t k = 0;
        for (size_t q = 0; q < n_rel; q++) {
            if (!use[q]) continue;
            const AsmRelation *R = &run->rels[(size_t)rel_idx[q]];
            ApRel *e = &rel[k++];
            e->a = act[R->a]; e->b = act[R->b];
            ap_effective(R, (run->charts[(size_t)R->a].flags & ASM_CHART_MIRROR) != 0, &e->theta, &e->tx, &e->ty);
            ap_seam_centres(run,R,e,centre_u,centre_v);
            e->w_xy = R->w_xy > 0 ? R->w_xy : 1.0;
            e->w_theta = R->w_theta > 1e-6 ? R->w_theta : 1e-6;
            e->w_xy_s = e->w_xy; e->w_theta_s = e->w_theta;
            if (AP_SEAM_FLOOR_VOX > 0.0) {
                double f = 1.0 / (1.0 + AP_SEAM_FLOOR_VOX * AP_SEAM_FLOOR_VOX * e->w_xy);
                e->w_xy_s *= f; e->w_theta_s *= f;
                if (AP_SEAM_FLOOR_IN_SOLVE) { e->w_xy *= f; e->w_theta *= f; }
            }
            e->rho = 1.0; e->e = 0.0; e->rxy = 0.0; e->rel_index = (size_t)rel_idx[q];
        }
    }
    /* pins: the fixed charts, plus the largest chart of every connected component without one (its gauge) */
    uint8_t *pin = ARENA_CALLOC(arena, na, 1);
    {
        UnionFind uf = UF_new(arena, (int32_t)na);
        for (size_t q = 0; q < nrel; q++) uf_union(&uf, rel[q].a, rel[q].b);
        uint8_t *has_fixed = ARENA_CALLOC(arena, na, 1);
        int32_t *largest = ARENA_ALLOC(arena, na * sizeof(int32_t));
        for (size_t i = 0; i < na; i++) largest[i] = -1;
        for (size_t i = 0; i < na; i++) {
            int32_t rt = uf_find(&uf, (int32_t)i);
            if (fixed != NULL && fixed[chart_of[i]]) { pin[i] = 1; has_fixed[rt] = 1; }
            if (largest[rt] < 0 || run->charts[chart_of[i]].area3d > run->charts[chart_of[largest[rt]]].area3d) largest[rt] = (int32_t)i;
        }
        for (size_t i = 0; i < na; i++) { int32_t rt = uf_find(&uf, (int32_t)i); if (!has_fixed[rt] && largest[rt] == (int32_t)i) pin[i] = 1; }
    }
    for (size_t i = 0; i < na; i++) { if (pin[i]) st->charts_pinned++; else st->charts_free++; }
    /* warm start */
    double *phi = ARENA_ALLOC(arena, na * sizeof(double)), *px = ARENA_ALLOC(arena, na * sizeof(double)), *py = ARENA_ALLOC(arena, na * sizeof(double));
    double *phi0 = ARENA_ALLOC(arena, na * sizeof(double)), *px0 = ARENA_ALLOC(arena, na * sizeof(double)), *py0 = ARENA_ALLOC(arena, na * sizeof(double));
    for (size_t i = 0; i < na; i++) {
        const AsmChart *ch = &run->charts[chart_of[i]];
        phi[i] = phi0[i] = ch->pose_theta;
        double co=cos(phi[i]),si=sin(phi[i]);
        px[i] = px0[i] = ch->pose_x+co*centre_u[i]-si*centre_v[i];
        py[i] = py0[i] = ch->pose_y+si*centre_u[i]+co*centre_v[i];
    }
    ApCoo coo; memset(&coo, 0, sizeof coo); coo.arena = arena;
    int n3 = (int)(3 * na);
    double *g = ARENA_CALLOC(arena, (size_t)n3, sizeof(double));
    double *dx = ARENA_CALLOC(arena, (size_t)n3, sizeof(double));
    double hold_xy = 1.0 / (AP_REFINE_HOLD_XY_VOX * AP_REFINE_HOLD_XY_VOX), hold_th = 1.0 / (AP_REFINE_HOLD_THETA_RAD * AP_REFINE_HOLD_THETA_RAD);
    int it;
    for (it = 0; it < AP_REFINE_ITERS; it++) {
        for (size_t q = 0; q < nrel; q++) {
            ApRel *e = &rel[q];
            double sx,sy; ap_seam_delta(e,phi[e->a],phi[e->b],&sx,&sy,NULL);
            double rx = px[e->b] - px[e->a] - sx;
            double ry = py[e->b] - py[e->a] - sy;
            double rt = ap_wrap(phi[e->b] - phi[e->a] - e->theta);
            e->e = sqrt((e->w_xy_s * (rx*rx + ry*ry) + e->w_theta_s * rt*rt) / 3.0);
            e->rxy = sqrt(rx*rx + ry*ry);
            double u = e->e / o->cauchy_c;
            e->rho = it < AP_REFINE_LS_ITERS ? 1.0 : 1.0 / (1.0 + u*u);
        }
        ap_coo_reset(&coo);
        memset(g, 0, (size_t)n3 * sizeof(double));
        for (size_t q = 0; q < nrel; q++) {
            const ApRel *e = &rel[q];
            double sx,sy,J[3][6]; ap_seam_delta(e,phi[e->a],phi[e->b],&sx,&sy,J);
            double rx = px[e->b] - px[e->a] - sx;
            double ry = py[e->b] - py[e->a] - sy;
            double rt = ap_wrap(phi[e->b] - phi[e->a] - e->theta);
            double W[3] = { e->rho * e->w_xy, e->rho * e->w_xy, e->rho * e->w_theta };
            double res[3] = { rx, ry, rt };
            int idx[6] = { 3*e->a, 3*e->a+1, 3*e->a+2, 3*e->b, 3*e->b+1, 3*e->b+2 };
            for (int r = 0; r < 6; r++) {
                double gi = 0.0;
                for (int k = 0; k < 3; k++) gi += J[k][r] * W[k] * res[k];
                g[idx[r]] += gi;
                for (int cc = 0; cc <= r; cc++) {
                    double v = 0.0;
                    for (int k = 0; k < 3; k++) v += J[k][r] * W[k] * J[k][cc];
                    if (v != 0.0) ap_coo_add(&coo, idx[r], idx[cc], v);
                }
            }
        }
        for (size_t i = 0; i < na; i++) {
            if (pin[i]) {
                for (int k = 0; k < 3; k++) { ap_coo_add(&coo, 3*(int)i + k, 3*(int)i + k, 1e6); g[3*i + (size_t)k] = 0.0; }
                continue;
            }
            /* the hold prior on the warm start */
            g[3*i] += hold_xy * (px[i] - px0[i]); g[3*i+1] += hold_xy * (py[i] - py0[i]); g[3*i+2] += hold_th * ap_wrap(phi[i] - phi0[i]);
            ap_coo_add(&coo, 3*(int)i, 3*(int)i, hold_xy + o->damping);
            ap_coo_add(&coo, 3*(int)i+1, 3*(int)i+1, hold_xy + o->damping);
            ap_coo_add(&coo, 3*(int)i+2, 3*(int)i+2, hold_th + o->damping);
        }
        for (int k = 0; k < n3; k++) g[k] = -g[k];
        memset(dx, 0, (size_t)n3 * sizeof(double));
        if (ap_solve(&coo, n3, g, dx) != 0) { fprintf(stderr, "  [refine] GN solve failed at iter %d\n", it); break; }
        double maxstep = 0.0;
        for (size_t i = 0; i < na; i++) {
            if (pin[i]) continue;
            px[i] += dx[3*i]; py[i] += dx[3*i+1]; phi[i] = ap_wrap(phi[i] + dx[3*i+2]);
            double m = fabs(dx[3*i]) + fabs(dx[3*i+1]) + 100.0 * fabs(dx[3*i+2]);
            if (m > maxstep) maxstep = m;
        }
        if (maxstep < 1e-4 && it >= AP_REFINE_LS_ITERS) { it++; break; }
    }
    st->iters = it;
    /* final residuals, robust weights, switches */
    double *elist = ARENA_ALLOC(arena, nrel * sizeof(double)), *xlist = ARENA_ALLOC(arena, nrel * sizeof(double));
    for (size_t q = 0; q < nrel; q++) {
        ApRel *e = &rel[q];
        double sx,sy; ap_seam_delta(e,phi[e->a],phi[e->b],&sx,&sy,NULL);
        double rx = px[e->b] - px[e->a] - sx;
        double ry = py[e->b] - py[e->a] - sy;
        double rt = ap_wrap(phi[e->b] - phi[e->a] - e->theta);
        e->e = sqrt((e->w_xy_s * (rx*rx + ry*ry) + e->w_theta_s * rt*rt) / 3.0);
        e->rxy = sqrt(rx*rx + ry*ry);
        double u = e->e / o->cauchy_c;
        e->rho = 1.0 / (1.0 + u*u);
        AsmRelation *R = &run->rels[e->rel_index];
        R->robust_w = e->rho; R->residual = e->e;
        R->flags &= ~(uint32_t)ASM_REL_SWITCHED;
        if (e->rho < o->switch_off) { R->flags |= ASM_REL_SWITCHED; st->switched++; }
        elist[q] = e->e; xlist[q] = e->rxy;
    }
    qsort(elist, nrel, sizeof(double), ap_cmp_double);
    qsort(xlist, nrel, sizeof(double), ap_cmp_double);
    st->resid_p50 = elist[nrel/2]; st->resid_p95 = elist[(size_t)((double)(nrel-1)*0.95)];
    st->xy_p50 = xlist[nrel/2]; st->xy_p95 = xlist[(size_t)((double)(nrel-1)*0.95)];
    for (size_t i = 0; i < na; i++) {
        if (pin[i]) continue;
        AsmChart *ch = &run->charts[chart_of[i]];
        double co=cos(phi[i]),si=sin(phi[i]);
        ch->pose_x = px[i]-co*centre_u[i]+si*centre_v[i];
        ch->pose_y = py[i]-si*centre_u[i]-co*centre_v[i]; ch->pose_theta = phi[i];
    }
    st->sec = ves_clock_sec() - t0;
    Arena_restore(arena, mark);
    return 0;
}

static void ap_test_rel(AsmRun *run, int a, int b, const double *gx, const double *gy, const double *gphi,
                        double noise, uint32_t flags, double w)
{
    /* true relative pose: p_a = R(theta) p_b + t with theta = phi_b - phi_a, t = R(-phi_a)(x_b - x_a) */
    AsmRelation r; memset(&r, 0, sizeof r);
    r.a = a; r.b = b;
    double th = gphi[b] - gphi[a];
    double dxg = gx[b] - gx[a], dyg = gy[b] - gy[a];
    double c = cos(-gphi[a]), s = sin(-gphi[a]);
    r.theta = th + noise; r.tx = c*dxg - s*dyg + noise; r.ty = s*dxg + c*dyg - noise;
    r.w_xy = w; r.w_theta = w * 400.0; r.robust_w = 1.0; r.flags = flags; r.n_corr = 10;
    AsmRun_push_rel(run, &r);
}

/* Changing a chart's local origin must not change any material position,
 * robust weight or switch. Exercise noisy cycles in both entry points and
 * with mirrored charts; the source vertices and observations stay fixed. */
static int ap_test_origin_invariance(int refine, int parity, int prior)
{
    enum { N = 8, E = 11 };
    const double gx[N] = {0,100,200,300,0,100,200,300}, gy[N] = {0,0,0,0,100,100,100,100};
    const double phi[N] = {0,0.1,-0.2,0.3,0.05,-0.1,0.15,-0.05};
    const int edges[E][2] = {{0,1},{1,2},{2,3},{4,5},{5,6},{6,7},{0,4},{1,5},{2,6},{3,7},{3,4}};
    const float local[8] = {-20,-30,20,-30,20,30,-20,30};
    double answer[2][N*3], weights[2][E]; uint32_t flags[2][E];
    int rc = 0;
    for (int shifted = 0; shifted < 2; shifted++) {
        AsmRun run; memset(&run,0,sizeof run); run.arena = Arena_new();
        double du[N], dv[N];
        for (int i = 0; i < N; i++) {
            du[i] = shifted*(i%2 ? 1 : -1)*(4000+200*i);
            dv[i] = shifted*(1200-300*i);
            int mirrored = parity && i%2;
            AsmChart c; memset(&c,0,sizeof c); c.id=i; c.nv=prior ? 40 : 4; c.area3d=1000+i;
            c.uv = ARENA_ALLOC(run.arena,c.nv*2*sizeof(float)); c.placed=1; c.component=0;
            c.xyz = ARENA_CALLOC(run.arena,c.nv*3,sizeof(float));
            for (size_t v=0; v<c.nv; v++) {
                double u=local[2*(v%4)],w=local[2*(v%4)+1];
                c.uv[2*v]=(float)((mirrored ? -1 : 1)*(u+du[i]));
                c.uv[2*v+1]=(float)(w+dv[i]);
                c.xyz[3*v]=(float)(gy[i]+sin(phi[i])*u+cos(phi[i])*w);
            }
            if (mirrored) c.flags |= ASM_CHART_MIRROR;
            c.pose_theta=phi[i]+0.01*(i%3);
            double co=cos(c.pose_theta),si=sin(c.pose_theta);
            c.pose_x=gx[i]+(i ? 3 : 0)-co*du[i]+si*dv[i];
            c.pose_y=gy[i]-(i ? 2 : 0)-si*du[i]-co*dv[i];
            AsmRun_push_chart(&run,&c);
        }
        for (int k=0; k<E; k++) {
            int a=edges[k][0],b=edges[k][1];
            ap_test_rel(&run,a,b,gx,gy,phi,0.002*(k%3),0,4);
            AsmRelation *r=&run.rels[k];
            if (k==E-1) { r->theta=0; r->tx=0; r->ty=0; }
            double co=cos(r->theta),si=sin(r->theta);
            r->tx += du[a]-co*du[b]+si*dv[b];
            r->ty += dv[a]-si*du[b]-co*dv[b];
            if (parity && a%2) { r->theta=-r->theta; r->tx=-r->tx; }
            if (parity && (a%2 != b%2)) r->flags |= ASM_REL_PARITY;
            r->corr_first=(int32_t)run.n_corr; r->corr_count=2;
            for (int j=0; j<2; j++) {
                AsmCorr cr; memset(&cr,0,sizeof cr); cr.va=1+j; cr.vb=j ? 3 : 0;
                AsmRun_push_corr(&run,&cr);
            }
        }
        AsmPoseOpts o; memset(&o,0,sizeof o); AsmPose_default_opts(&o);
        if (refine) {
            int32_t ids[E]={0,1,2,3,4,5,6,7,8,9,10}; uint8_t fixed[N]={1,0,0,0,0,0,0,0};
            AsmPoseRefineStats st; rc |= AsmPose_refine(&run,&o,ids,E,fixed,&st);
        } else { AsmPoseStats st; rc |= AsmPose_solve(&run,&o,&st); }
        double points[N*2];
        for (int i=0; i<N; i++) {
            const AsmChart *c=&run.charts[i]; double co=cos(c->pose_theta),si=sin(c->pose_theta);
            points[2*i]=c->pose_x+co*du[i]-si*dv[i];
            points[2*i+1]=c->pose_y+si*du[i]+co*dv[i];
        }
        double co=cos(run.charts[0].pose_theta),si=sin(run.charts[0].pose_theta);
        for (int i=0; i<N; i++) {
            double x=points[2*i]-points[0],y=points[2*i+1]-points[1];
            answer[shifted][3*i]=co*x+si*y; answer[shifted][3*i+1]=-si*x+co*y;
            answer[shifted][3*i+2]=ap_wrap(run.charts[i].pose_theta-run.charts[0].pose_theta);
        }
        for (int k=0; k<E; k++) { weights[shifted][k]=run.rels[k].robust_w; flags[shifted][k]=run.rels[k].flags; }
        Arena_dispose(&run.arena);
    }
    double error=0,weight_error=0; int switches=0;
    for (int i=0; i<N*3; i++) error=fmax(error,fabs(answer[0][i]-answer[1][i]));
    for (int k=0; k<E; k++) {
        weight_error=fmax(weight_error,fabs(weights[0][k]-weights[1][k]));
        switches += flags[0][k]!=flags[1][k];
    }
    if (rc || error>1e-5 || weight_error>1e-7 || switches) {
        fprintf(stderr,"  asm_pose selftest FAIL: origin invariance refine=%d parity=%d prior=%d material %.9g weight %.9g switches %d\n",refine,parity,prior,error,weight_error,switches);
        return 1;
    }
    return 0;
}

int AsmPose_selftest(void)
{
    int fails = 0;
    {
        AsmPoseOpts options; memset(&options,0x5a,sizeof options); AsmPose_default_opts(&options);
        if (options.axis!=NULL || options.verbose_prior!=0) {
            fprintf(stderr,"  asm_pose selftest FAIL: defaults leave optional fields uninitialized\n"); fails++;
        }
    }
    AsmRun run; memset(&run, 0, sizeof run);
    run.arena = Arena_new();
    enum { N = 6 };
    double gx[N] = { 0, 100, 200, 0, 100, 200 }, gy[N] = { 0, 0, 0, 100, 100, 100 };
    double gphi[N] = { 0, 0.1, -0.2, 0.3, 0.05, -0.1 };
    float dummy_uv[2] = { 0, 0 };
    for (int i = 0; i < N; i++) {
        AsmChart c; memset(&c, 0, sizeof c);
        c.id = i; c.cube = i; c.nv = 1; c.nf = 0; c.uv = dummy_uv; c.area3d = 1000.0 + i;
        AsmRun_push_chart(&run, &c);
    }
    /* grid edges with tiny noise */
    int edges[7][2] = { {0,1},{1,2},{3,4},{4,5},{0,3},{1,4},{2,5} };
    for (int k = 0; k < 7; k++) ap_test_rel(&run, edges[k][0], edges[k][1], gx, gy, gphi, 0.001 * (k % 3), 0, 4.0);
    /* one grossly wrong relation (a sheet switch): claims 2 and 3 coincide */
    {
        AsmRelation r; memset(&r, 0, sizeof r);
        r.a = 2; r.b = 3; r.theta = 0.0; r.tx = 0.0; r.ty = 0.0; r.w_xy = 4.0; r.w_theta = 1600.0; r.robust_w = 1.0; r.n_corr = 10;
        AsmRun_push_rel(&run, &r);
    }
    AsmPoseOpts o; AsmPose_default_opts(&o);
    AsmPoseStats st;
    if (AsmPose_solve(&run, &o, &st) != 0) { fprintf(stderr, "  asm_pose selftest FAIL: solve rc\n"); fails++; }
    /* compare relative poses to the truth (global gauge is free) */
    const AsmChart *A = &run.charts[0];
    double c0 = cos(A->pose_theta), s0 = sin(A->pose_theta);
    double maxerr = 0.0;
    for (int i = 1; i < N; i++) {
        const AsmChart *B = &run.charts[i];
        /* bring B's pose into A's frame and compare with the truth in chart 0's frame */
        double dx = B->pose_x - A->pose_x, dy = B->pose_y - A->pose_y;
        double lx = c0*dx + s0*dy, ly = -s0*dx + c0*dy;   /* R(-phi_a) d */
        double tc = cos(-gphi[0]), ts = sin(-gphi[0]);
        double tx = tc*(gx[i]-gx[0]) - ts*(gy[i]-gy[0]), ty = ts*(gx[i]-gx[0]) + tc*(gy[i]-gy[0]);
        double e = hypot(lx - tx, ly - ty) + 100.0 * fabs(ap_wrap((B->pose_theta - A->pose_theta) - (gphi[i] - gphi[0])));
        if (e > maxerr) maxerr = e;
    }
    if (maxerr > 1.0) { fprintf(stderr, "  asm_pose selftest FAIL: pose error %g\n", maxerr); fails++; }
    if (!(run.rels[7].flags & ASM_REL_SWITCHED)) { fprintf(stderr, "  asm_pose selftest FAIL: wrong relation not switched off (w %g)\n", run.rels[7].robust_w); fails++; }
    for (int k = 0; k < 7; k++) if (run.rels[k].flags & ASM_REL_SWITCHED) { fprintf(stderr, "  asm_pose selftest FAIL: good relation %d switched off\n", k); fails++; }
    if (st.components != 1) { fprintf(stderr, "  asm_pose selftest FAIL: components %zu\n", st.components); fails++; }
    /* parity: a mirrored chart 1 must be recovered when its relations carry the parity bit */
    {
        AsmRun run2; memset(&run2, 0, sizeof run2); run2.arena = Arena_new();
        for (int i = 0; i < 3; i++) { AsmChart c; memset(&c, 0, sizeof c); c.id = i; c.nv = 1; c.uv = dummy_uv; c.area3d = 1000.0 - i; AsmRun_push_chart(&run2, &c); }
        double x3[3] = { 0, 100, 200 }, y3[3] = { 0, 0, 0 }, p3[3] = { 0, 0, 0 };
        /* chart 1 mirrored: relations with 1 carry PARITY; the stored transform is in unmirrored frames,
         * which for a pure translation with theta 0 is the same */
        ap_test_rel(&run2, 0, 1, x3, y3, p3, 0.0, ASM_REL_PARITY, 4.0);
        ap_test_rel(&run2, 1, 2, x3, y3, p3, 0.0, ASM_REL_PARITY, 4.0);
        AsmPoseStats s2;
        AsmPose_solve(&run2, &o, &s2);
        if (!(run2.charts[1].flags & ASM_CHART_MIRROR) || (run2.charts[0].flags & ASM_CHART_MIRROR) || (run2.charts[2].flags & ASM_CHART_MIRROR)) {
            fprintf(stderr, "  asm_pose selftest FAIL: parity colouring\n"); fails++;
        }
        if (s2.parity_conflicts != 0) { fprintf(stderr, "  asm_pose selftest FAIL: spurious parity conflict\n"); fails++; }
        Arena_dispose(&run2.arena);
    }
    Arena_dispose(&run.arena);
    /* REFINE (2026-09-10): four charts on a line, chart 0 fixed; charts 1 and 2 warm-started 30 vox off
     * their seams, chart 3 right; a wrong relation 0-3 claiming 3 sits 100 vox nearer: the seams pull 1
     * and 2 back within 0.1 vox, chart 0 stays, the wrong relation ends SWITCHED; a detached pair
     * without a fixed chart (4-5) solves with its larger chart pinned */
    {
        AsmRun run3; memset(&run3, 0, sizeof run3); run3.arena = Arena_new();
        double x6[6] = { 0, 100, 200, 300, 1000, 1100 }, y6[6] = { 0, 0, 0, 0, 0, 0 }, p6[6] = { 0, 0, 0, 0, 0, 0 };
        for (int i = 0; i < 6; i++) { AsmChart c; memset(&c, 0, sizeof c); c.id = i; c.nv = 1; c.uv = dummy_uv; c.area3d = 1000.0 + (i == 5 ? 500.0 : 0.0); c.placed = 1; c.component = 0; AsmRun_push_chart(&run3, &c); }
        ap_test_rel(&run3, 0, 1, x6, y6, p6, 0.0, 0, 100.0);
        ap_test_rel(&run3, 1, 2, x6, y6, p6, 0.0, 0, 100.0);
        ap_test_rel(&run3, 2, 3, x6, y6, p6, 0.0, 0, 100.0);
        ap_test_rel(&run3, 4, 5, x6, y6, p6, 0.0, 0, 100.0);
        /* the wrong relation carries a readmitted seam's weight (fewer correspondences): in a 4-cycle of equal
         * weights a single outlier cannot be told from the chain (every residual reads 25 vox) */
        { AsmRelation r; memset(&r, 0, sizeof r); r.a = 0; r.b = 3; r.theta = 0.0; r.tx = 200.0; r.ty = 0.0; r.w_xy = 10.0; r.w_theta = 4000.0; r.robust_w = 1.0; r.n_corr = 10; AsmRun_push_rel(&run3, &r); }
        for (int i = 0; i < 6; i++) { run3.charts[i].pose_x = x6[i]; run3.charts[i].pose_y = 0.0; run3.charts[i].pose_theta = 0.0; }
        run3.charts[1].pose_x += 30.0; run3.charts[2].pose_x += 30.0; run3.charts[5].pose_x += 20.0;
        uint8_t fixed6[6] = { 1, 0, 0, 0, 0, 0 };
        int32_t idx5[5] = { 0, 1, 2, 3, 4 };
        AsmPoseRefineStats rs;
        AsmPose_refine(&run3, &o, idx5, 5, fixed6, &rs);
        double e1 = fabs(run3.charts[1].pose_x - 100.0) + fabs(run3.charts[1].pose_y), e2 = fabs(run3.charts[2].pose_x - 200.0) + fabs(run3.charts[2].pose_y), e0 = fabs(run3.charts[0].pose_x) + fabs(run3.charts[0].pose_y);
        double e45 = fabs((run3.charts[5].pose_x - run3.charts[4].pose_x) - 100.0);
        if (e1 > 0.1 || e2 > 0.1 || e0 > 1e-9 || !(run3.rels[4].flags & ASM_REL_SWITCHED) || (run3.rels[0].flags & ASM_REL_SWITCHED) || e45 > 0.1 || rs.charts_pinned != 2 || fabs(run3.charts[5].pose_x - 1120.0) > 1e-9) {
            fprintf(stderr, "  asm_pose selftest FAIL: refine (chart 1 err %.3f, chart 2 err %.3f, chart 0 moved %.3g, wrong rel switched %d, pair err %.3f, pinned %zu, chart 5 x %.1f)\n",
                    e1, e2, e0, (run3.rels[4].flags & ASM_REL_SWITCHED) != 0, e45, rs.charts_pinned, run3.charts[5].pose_x);
            fails++;
        }
        Arena_dispose(&run3.arena);
    }
    for (int refine=0; refine<2; refine++) for (int parity=0; parity<2; parity++) for (int prior=0; prior<2; prior++)
        fails += ap_test_origin_invariance(refine,parity,prior);
    if (fails == 0) fprintf(stderr, "  asm_pose selftest: all passed\n");
    return fails;
}
