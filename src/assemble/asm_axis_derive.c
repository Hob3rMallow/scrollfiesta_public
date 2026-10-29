/* asm_axis_derive.c -- the scroll axis derived from the material (see asm_axis_derive.h). */
#include "asm_axis_derive.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/eig3.h"
#include "../common/pipeline_constants.h"
#include "../common/ves_platform.h"
#include "asm_report.h"

#define AD_PI 3.14159265358979323846

void AsmAxisDerive_default_opts(AsmAxisDeriveOpts *o)
{
    memset(o, 0, sizeof *o);
    o->slab_vox = ASM_AXIS_DERIVE_SLAB_VOX;
    o->iters = ASM_AXIS_DERIVE_ITERS;
    o->per_chart = ASM_AXIS_DERIVE_SAMPLES_PER_CHART;
    o->prior = NULL;
}

int AsmAxis_normal_direction(const float *nrm, const float *w, size_t n, double dir[3], double *cond, double *frac)
{
    if (!nrm || !w || !n || !dir) return -1;
    double M[3][3] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } }, sw = 0.0;
    for (size_t i = 0; i < n; i++) {
        double nn[3] = { nrm[i*3], nrm[i*3+1], nrm[i*3+2] }, wi = w[i];
        double L = sqrt(nn[0]*nn[0] + nn[1]*nn[1] + nn[2]*nn[2]);
        if (!(L > 1e-9) || !(wi > 0.0) || !isfinite(wi)) continue;
        for (int k = 0; k < 3; k++) nn[k] /= L;
        for (int r = 0; r < 3; r++) for (int c = 0; c < 3; c++) M[r][c] += wi * nn[r] * nn[c];
        sw += wi;
    }
    if (!(sw > 0.0)) return -1;
    double evals[3], evecs[3][3];
    Eig3_sym(M, evals, evecs);   /* ascending: evals[0] is the null direction */
    for (int k = 0; k < 3; k++) dir[k] = evecs[k][0];
    double L = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
    if (!(L > 1e-9)) return -1;
    for (int k = 0; k < 3; k++) dir[k] /= L;
    if (dir[0] < 0.0) for (int k = 0; k < 3; k++) dir[k] = -dir[k];   /* the scroll axis points up */
    if (cond) *cond = evals[1] > 0.0 ? fmax(0.0, evals[0]) / evals[1] : 1.0;
    if (frac) {
        double inside = 0.0;
        for (size_t i = 0; i < n; i++) {
            double nn[3] = { nrm[i*3], nrm[i*3+1], nrm[i*3+2] };
            double l = sqrt(nn[0]*nn[0] + nn[1]*nn[1] + nn[2]*nn[2]);
            if (!(l > 1e-9) || !(w[i] > 0.0)) continue;
            double d = fabs((nn[0]*dir[0] + nn[1]*dir[1] + nn[2]*dir[2]) / l);
            if (d <= 0.2) inside += w[i];
        }
        *frac = inside / sw;
    }
    return 0;
}

/* ---- projection onto the polyline, by local descent from the z-nearest row ---------------------
 * AsmAxis_project scans a 2,000-vox window of segments per call; a million samples per iteration
 * need the O(1) version: the polyline is smooth and monotone in z, so the nearest segment is a
 * local descent away from the row at the sample's z. */

static double ad_seg_d2(const AsmAxis *a, size_t i, const double p[3], double *t_out)
{
    const double *A = &a->p[i*3], *B = &a->p[(i+1)*3];
    double e[3] = { B[0]-A[0], B[1]-A[1], B[2]-A[2] };
    double L2 = e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
    if (L2 < 1e-18) { *t_out = 0.0; return 1e300; }
    double t = ((p[0]-A[0])*e[0] + (p[1]-A[1])*e[1] + (p[2]-A[2])*e[2]) / L2;
    if (i > 0 && t < 0.0) t = 0.0;
    if (i + 2 < a->n && t > 1.0) t = 1.0;
    double q[3] = { A[0] + t*e[0], A[1] + t*e[1], A[2] + t*e[2] };
    *t_out = t;
    return (p[0]-q[0])*(p[0]-q[0]) + (p[1]-q[1])*(p[1]-q[1]) + (p[2]-q[2])*(p[2]-q[2]);
}

static void ad_project(const AsmAxis *a, const double p[3], double *s, double dir[3], double foot[3])
{
    if (a->n <= 64) { AsmAxis_project(a, p, s, NULL, dir, foot); return; }
    size_t l = 0, h = a->n - 1;
    while (l < h) { size_t mid = (l + h) / 2; if (a->p[mid*3] < p[0]) l = mid + 1; else h = mid; }
    size_t seg = l > 0 ? l - 1 : 0;
    if (seg > a->n - 2) seg = a->n - 2;
    double t = 0.0, best = ad_seg_d2(a, seg, p, &t), tt;
    size_t bseg = seg;
    for (size_t k = seg; k > 0; k--) { double d2 = ad_seg_d2(a, k-1, p, &tt); if (d2 < best) { best = d2; bseg = k-1; t = tt; } else break; }
    for (size_t k = seg; k + 2 < a->n; k++) { double d2 = ad_seg_d2(a, k+1, p, &tt); if (d2 < best) { best = d2; bseg = k+1; t = tt; } else break; }
    const double *A = &a->p[bseg*3], *B = &a->p[(bseg+1)*3];
    double e[3] = { B[0]-A[0], B[1]-A[1], B[2]-A[2] };
    double L = sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
    double tc = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
    double dd[3] = { (1.0-tc)*a->t[bseg*3] + tc*a->t[(bseg+1)*3], (1.0-tc)*a->t[bseg*3+1] + tc*a->t[(bseg+1)*3+1], (1.0-tc)*a->t[bseg*3+2] + tc*a->t[(bseg+1)*3+2] };
    double dl = sqrt(dd[0]*dd[0] + dd[1]*dd[1] + dd[2]*dd[2]);
    if (dl < 1e-12) { dd[0] = e[0] / L; dd[1] = e[1] / L; dd[2] = e[2] / L; dl = 1.0; }
    if (s) *s = a->s[bseg] + t * L;
    if (dir) for (int k = 0; k < 3; k++) dir[k] = dd[k] / dl;
    if (foot) for (int k = 0; k < 3; k++) foot[k] = A[k] + t * e[k];
}

/* the frame of the polyline at arc length s: the point, the tangent and (e1, e2) as AsmAxis_frame */
static void ad_frame_at(const AsmAxis *a, double s, double c0[3], double t[3], double e1[3], double e2[3])
{
    size_t l = 0, h = a->n - 1;
    while (l < h) { size_t mid = (l + h) / 2; if (a->s[mid] < s) l = mid + 1; else h = mid; }
    size_t i = l > 0 ? l - 1 : 0;
    if (i > a->n - 2) i = a->n - 2;
    double L = a->s[i+1] - a->s[i];
    double f = L > 1e-12 ? (s - a->s[i]) / L : 0.0;
    for (int k = 0; k < 3; k++) {
        c0[k] = a->p[i*3+k] + f * (a->p[(i+1)*3+k] - a->p[i*3+k]);
        t[k] = (1.0-f) * a->t[i*3+k] + f * a->t[(i+1)*3+k];
    }
    double tl = sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]);
    if (tl < 1e-12) { t[0] = 1.0; t[1] = 0.0; t[2] = 0.0; tl = 1.0; }
    for (int k = 0; k < 3; k++) t[k] /= tl;
    double v1[3] = { -t[1]*t[0], 1.0 - t[1]*t[1], -t[1]*t[2] };
    double l1 = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
    if (l1 < 1e-9) { v1[0] = 0.0; v1[1] = 0.0; v1[2] = 1.0; l1 = 1.0; }
    for (int k = 0; k < 3; k++) e1[k] = v1[k] / l1;
    e2[0] = t[1]*e1[2] - t[2]*e1[1]; e2[1] = t[2]*e1[0] - t[0]*e1[2]; e2[2] = t[0]*e1[1] - t[1]*e1[0];
}

/* ---- the 2-D robust normal-line intersection ------------------------------------------------ */

typedef struct AdFit {
    double c[2];
    size_t n, kept;
    double lam_min, lam_max, cond, resid_p50, unc_major, unc_minor, aperture_deg, void_ratio, sum_w;
    int    sectors, one_sided, ok;
} AdFit;

static int ad_cmp_double(const void *x, const void *y)
{ double a = *(const double *)x, b = *(const double *)y; return (a > b) - (a < b); }

/* weighted median of v[n] with weights w[n] (both permuted by the sort of the pairs) */
static double ad_wmedian(double *pairs /* [2n]: value, weight */, size_t n)
{
    if (!n) return 0.0;
    qsort(pairs, n, 2 * sizeof(double), ad_cmp_double);
    double total = 0.0; for (size_t i = 0; i < n; i++) total += pairs[2*i+1];
    double acc = 0.0;
    for (size_t i = 0; i < n; i++) { acc += pairs[2*i+1]; if (acc >= 0.5 * total) return pairs[2*i]; }
    return pairs[2*(n-1)];
}

/* Minimise sum_i w_i rho(d_i) over the centre c, d_i = the distance from c to the line through
 * (pa_i, pb_i) along the unit in-plane normal (ma_i, mb_i), Cauchy IRLS with the scale annealed
 * geometrically from sig0 to sig1 over `rounds`; then every diagnostic of the fit.  `pairs` is
 * scratch of 2n doubles.  Returns 0 when the 2x2 system was solvable at every round. */
static int ad_fit2d(const double *pa, const double *pb, const double *ma, const double *mb, const double *w, size_t n,
                    const double c0[2], double sig0, double sig1, int rounds, AdFit *f, double *pairs)
{
    memset(f, 0, sizeof *f);
    f->c[0] = c0[0]; f->c[1] = c0[1]; f->n = n;
    if (n < 3) return -1;
    double a11 = 0.0, a12 = 0.0, a22 = 0.0, med = 1e300;
    for (int r = 0; r < rounds; r++) {
        /* the scale: the annealing schedule, never below 1.5x the previous round's median residual
         * (a section whose normal lines miss the centre by 50 vox is fitted at that scale, as the
         * tracker's 60% trim does; a fixed 12-vox scale lets a handful of lines dominate) */
        double sig = rounds > 1 ? sig0 * pow(sig1 / sig0, (double)r / (double)(rounds - 1)) : sig1;
        if (med < 1e299 && 1.5 * med > sig) sig = 1.5 * med;
        double b1 = 0.0, b2 = 0.0;
        a11 = a12 = a22 = 0.0;
        for (size_t i = 0; i < n; i++) {
            double dy = f->c[0] - pa[i], dx = f->c[1] - pb[i];
            double along = dy*ma[i] + dx*mb[i];
            double ry = dy - along*ma[i], rx = dx - along*mb[i];
            double d = sqrt(ry*ry + rx*rx), u = d / sig;
            double wi = w[i] / (1.0 + u*u);
            pairs[2*i] = d; pairs[2*i+1] = w[i];
            double p11 = 1.0 - ma[i]*ma[i], p12 = -ma[i]*mb[i], p22 = 1.0 - mb[i]*mb[i];
            a11 += wi * p11; a12 += wi * p12; a22 += wi * p22;
            b1 += wi * (p11 * pa[i] + p12 * pb[i]);
            b2 += wi * (p12 * pa[i] + p22 * pb[i]);
        }
        med = ad_wmedian(pairs, n);
        double det = a11 * a22 - a12 * a12;
        if (!(det > 1e-12)) return -1;
        f->c[0] = ( a22 * b1 - a12 * b2) / det;
        f->c[1] = (-a12 * b1 + a11 * b2) / det;
    }
    /* the diagnostics at the final centre, with the final scale */
    double tr = a11 + a22, disc = sqrt(fmax(0.0, 0.25 * (a11 - a22) * (a11 - a22) + a12 * a12));
    f->lam_max = 0.5 * tr + disc; f->lam_min = 0.5 * tr - disc;
    f->cond = f->lam_max > 0.0 ? f->lam_min / f->lam_max : 0.0;
    size_t nk = 0; double sw = 0.0, a_in = 0.0, a_ann = 0.0;
    int hit[12] = { 0 };
    for (size_t i = 0; i < n; i++) {
        double dy = f->c[0] - pa[i], dx = f->c[1] - pb[i];
        double along = dy*ma[i] + dx*mb[i];
        double ry = dy - along*ma[i], rx = dx - along*mb[i];
        double d = sqrt(ry*ry + rx*rx);
        double rr = sqrt(dy*dy + dx*dx);
        if (rr < ASM_AXIS_DERIVE_VOID_R) a_in += w[i]; else if (rr < 3.0 * ASM_AXIS_DERIVE_VOID_R) a_ann += w[i];
        pairs[2*i] = d; pairs[2*i+1] = w[i]; sw += w[i];
        if (d > 2.0 * fmax(sig1, 1.5 * med)) continue;   /* a Cauchy weight below one fifth at the final scale: the diagnostic 'kept' and the sectors */
        nk++;
        int sct = (int)((atan2(-dx, -dy) + AD_PI) * (6.0 / AD_PI));
        if (sct < 0) sct = 0;
        if (sct > 11) sct = 11;
        hit[sct] = 1;
    }
    f->kept = nk; f->sum_w = sw;
    /* the residual median over the WHOLE window: real sections are ovals whose normal lines
     * miss the centre by their own eccentricity, so a threshold on the kept count would abstain
     * on every honest slab (the 4x5x5's own table reads 15-22 vox here) */
    f->resid_p50 = ad_wmedian(pairs, n);
    for (int k = 0; k < 12; k++) f->sectors += hit[k];
    f->one_sided = f->sectors <= 6;
    {
        double d_in = a_in / (AD_PI * ASM_AXIS_DERIVE_VOID_R * ASM_AXIS_DERIVE_VOID_R);
        double d_ann = a_ann / (AD_PI * 8.0 * ASM_AXIS_DERIVE_VOID_R * ASM_AXIS_DERIVE_VOID_R);
        f->void_ratio = d_ann > 1e-12 ? d_in / d_ann : 1.0;
    }
    /* the aperture: the angular span of the kept normal LINES (angles modulo 180 degrees) */
    {
        size_t na = 0;
        for (size_t i = 0; i < n; i++) {
            double dy = f->c[0] - pa[i], dx = f->c[1] - pb[i];
            double along = dy*ma[i] + dx*mb[i];
            double ry = dy - along*ma[i], rx = dx - along*mb[i];
            if (sqrt(ry*ry + rx*rx) > 2.0 * fmax(sig1, 1.5 * med)) continue;
            double ang = atan2(mb[i], ma[i]);
            if (ang < 0.0) ang += AD_PI;
            pairs[na++] = ang;
        }
        if (na >= 2) {
            qsort(pairs, na, sizeof(double), ad_cmp_double);
            double gap = pairs[0] + AD_PI - pairs[na-1];
            for (size_t i = 1; i < na; i++) if (pairs[i] - pairs[i-1] > gap) gap = pairs[i] - pairs[i-1];
            f->aperture_deg = (AD_PI - gap) * 180.0 / AD_PI;
        }
    }
    double wbar = n ? sw / (double)n : 1.0;
    f->unc_major = f->lam_min > 1e-12 ? f->resid_p50 * sqrt(wbar / f->lam_min) : 1e300;
    f->unc_minor = f->lam_max > 1e-12 ? f->resid_p50 * sqrt(wbar / f->lam_max) : 1e300;
    f->ok = n >= (size_t)ASM_AXIS_DERIVE_MIN_SAMPLES && f->cond >= ASM_AXIS_DERIVE_MIN_COND && f->resid_p50 <= ASM_AXIS_DERIVE_MAX_RESID;
    return 0;
}

/* The windowed fit: only the samples within R of the running centre feed ad_fit2d, and the
 * window follows the estimate to a fixed point (the tracker's doctrine: the whole-slab normal
 * scatter of flattened wraps converges on their evolute cusps, and hundreds of thousands of
 * distant samples outvote any robust weight).  The first rounds use 2R so a prior a few hundred
 * vox off still reaches the material.  s* are scratch of m doubles, pairs of 2m.  Returns 0
 * when a fit converged (f holds it), -1 when the window ran dry or the system was singular. */
static int ad_fit_windowed(const double *pa, const double *pb, const double *ma, const double *mb, const double *ww, size_t m,
                           const double c0[2], double R, AdFit *f,
                           double *sa, double *sb, double *sma, double *smb, double *sw, double *pairs)
{
    double c[2] = { c0[0], c0[1] };
    int fine = 0;
    memset(f, 0, sizeof *f);
    f->c[0] = c0[0]; f->c[1] = c0[1];
    for (int it = 0; it < ASM_AXIS_DERIVE_WINDOW_ROUNDS; it++) {
        double Rw = it < 2 ? 2.0 * R : R;
        size_t k = 0;
        for (size_t i = 0; i < m; i++) {
            double dy = pa[i] - c[0], dx = pb[i] - c[1];
            if (dy*dy + dx*dx > Rw*Rw) continue;
            sa[k] = pa[i]; sb[k] = pb[i]; sma[k] = ma[i]; smb[k] = mb[i]; sw[k] = ww[i]; k++;
        }
        if (k < (size_t)ASM_AXIS_DERIVE_MIN_SAMPLES) { f->n = k; return -1; }
        if (ad_fit2d(sa, sb, sma, smb, sw, k, c, it < 2 ? 2.0 * ASM_AXIS_DERIVE_CAUCHY_VOX0 : ASM_AXIS_DERIVE_CAUCHY_VOX0, ASM_AXIS_DERIVE_CAUCHY_VOX, ASM_AXIS_DERIVE_IRLS_ROUNDS, f, pairs)) return -1;
        double move = hypot(f->c[0] - c[0], f->c[1] - c[1]);
        c[0] = f->c[0]; c[1] = f->c[1]; fine = 1;
        if (move < 1.0 && it >= 2) break;
    }
    return fine ? 0 : -1;
}

/* ---- the seed on z-slabs (no prior) ------------------------------------------------------------ */

typedef struct AdCluster { double y, x, resid, void_ratio; size_t members; } AdCluster;

/* The seed: slabs PERPENDICULAR TO `dir` (the normals' null direction, measured without a
 * centre), each fitted from a grid of starts, the void-qualified lowest-residual attractor
 * winning.  Slicing by world z instead would smear every section of a tilted scroll by
 * slab x tan(tilt) -- 42 vox at 18 degrees over a 128-vox slab, which is four wraps. */
static AsmAxis *ad_seed(Arena_T arena, Arena_T scratch, const float *p, const float *nrm, const float *w, size_t n,
                        const double dir[3], AsmAxisDeriveReport *rep)
{
    double e1[3] = { -dir[1]*dir[0], 1.0 - dir[1]*dir[1], -dir[1]*dir[2] };
    double l1 = sqrt(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]);
    if (l1 < 1e-9) { e1[0] = 0.0; e1[1] = 0.0; e1[2] = 1.0; l1 = 1.0; }
    for (int k = 0; k < 3; k++) e1[k] /= l1;
    double e2[3] = { dir[1]*e1[2] - dir[2]*e1[1], dir[2]*e1[0] - dir[0]*e1[2], dir[0]*e1[1] - dir[1]*e1[0] };
    double c0[3] = { 0.0, 0.0, 0.0 };
    for (size_t i = 0; i < n; i++) for (int k = 0; k < 3; k++) c0[k] += p[i*3+k];
    for (int k = 0; k < 3; k++) c0[k] /= (double)n;
    double zlo = 1e300, zhi = -1e300;
    for (size_t i = 0; i < n; i++) {
        double s = (p[i*3]-c0[0])*dir[0] + (p[i*3+1]-c0[1])*dir[1] + (p[i*3+2]-c0[2])*dir[2];
        if (s < zlo) zlo = s; if (s > zhi) zhi = s;
    }
    if (!(zhi > zlo)) return NULL;
    size_t nslab = (size_t)((zhi - zlo) / ASM_AXIS_DERIVE_SEED_SLAB) + 1;
    double *rows = ARENA_ALLOC(arena, nslab * 3 * sizeof(double));
    size_t nrow = 0;
    Arena_Mark mark = Arena_save(scratch);
    double *sax = ARENA_ALLOC(scratch, n * sizeof(double));   /* every sample's axial coordinate */
    for (size_t i = 0; i < n; i++)
        sax[i] = (p[i*3]-c0[0])*dir[0] + (p[i*3+1]-c0[1])*dir[1] + (p[i*3+2]-c0[2])*dir[2];
    double *pa = ARENA_ALLOC(scratch, n * sizeof(double)), *pb = ARENA_ALLOC(scratch, n * sizeof(double));
    double *ma = ARENA_ALLOC(scratch, n * sizeof(double)), *mb = ARENA_ALLOC(scratch, n * sizeof(double));
    double *ww = ARENA_ALLOC(scratch, n * sizeof(double)), *pairs = ARENA_ALLOC(scratch, 2 * n * sizeof(double));
    double *sa = ARENA_ALLOC(scratch, n * sizeof(double)), *sb = ARENA_ALLOC(scratch, n * sizeof(double));
    double *sma = ARENA_ALLOC(scratch, n * sizeof(double)), *smb = ARENA_ALLOC(scratch, n * sizeof(double)), *sw = ARENA_ALLOC(scratch, n * sizeof(double));
    AdCluster cl[64];
    for (size_t sl = 0; sl < nslab; sl++) {
        double z0 = zlo + (double)sl * ASM_AXIS_DERIVE_SEED_SLAB, z1 = z0 + ASM_AXIS_DERIVE_SEED_SLAB;
        /* the slab's in-plane samples (normals not along the axis), subsampled to the seed budget */
        size_t m = 0, in_slab = 0;
        for (size_t i = 0; i < n; i++) {
            double nt = nrm[i*3]*dir[0] + nrm[i*3+1]*dir[1] + nrm[i*3+2]*dir[2];
            if (sax[i] >= z0 && sax[i] < z1 && fabs(nt) <= ASM_AXIS_DERIVE_MAX_AXIAL) in_slab++;
        }
        size_t stride = in_slab > (size_t)ASM_AXIS_DERIVE_SEED_SAMPLES ? in_slab / (size_t)ASM_AXIS_DERIVE_SEED_SAMPLES + 1 : 1;
        size_t seen = 0;
        double ylo = 1e300, yhi = -1e300, xlo = 1e300, xhi = -1e300;
        for (size_t i = 0; i < n; i++) {
            double nt = nrm[i*3]*dir[0] + nrm[i*3+1]*dir[1] + nrm[i*3+2]*dir[2];
            if (!(sax[i] >= z0 && sax[i] < z1) || fabs(nt) > ASM_AXIS_DERIVE_MAX_AXIAL) continue;
            if (seen++ % stride) continue;
            double my = nrm[i*3]*e1[0] + nrm[i*3+1]*e1[1] + nrm[i*3+2]*e1[2];
            double mx = nrm[i*3]*e2[0] + nrm[i*3+1]*e2[1] + nrm[i*3+2]*e2[2];
            double ml = sqrt(my*my + mx*mx);
            if (ml < 0.3) continue;
            double d[3] = { p[i*3]-c0[0], p[i*3+1]-c0[1], p[i*3+2]-c0[2] };
            pa[m] = d[0]*e1[0] + d[1]*e1[1] + d[2]*e1[2];
            pb[m] = d[0]*e2[0] + d[1]*e2[1] + d[2]*e2[2];
            ma[m] = my / ml; mb[m] = mx / ml; ww[m] = w[i] * (double)stride; m++;
            if (pa[m-1] < ylo) ylo = pa[m-1]; if (pa[m-1] > yhi) yhi = pa[m-1];
            if (pb[m-1] < xlo) xlo = pb[m-1]; if (pb[m-1] > xhi) xhi = pb[m-1];
        }
        if (m < (size_t)ASM_AXIS_DERIVE_MIN_SAMPLES) continue;
        size_t nc = 0, starts = 0, converged = 0;
        AdFit reject; memset(&reject, 0, sizeof reject); reject.resid_p50 = 1e300;
        for (double y = ylo + 0.5 * ASM_AXIS_DERIVE_SEED_GRID; y < yhi; y += ASM_AXIS_DERIVE_SEED_GRID)
        for (double x = xlo + 0.5 * ASM_AXIS_DERIVE_SEED_GRID; x < xhi; x += ASM_AXIS_DERIVE_SEED_GRID) {
            double c[2] = { y, x }; AdFit f;
            starts++;
            if (ad_fit_windowed(pa, pb, ma, mb, ww, m, c, ASM_AXIS_DERIVE_LOCAL_R, &f, sa, sb, sma, smb, sw, pairs)) continue;
            converged++;
            if (!f.ok) { if (f.n >= (size_t)ASM_AXIS_DERIVE_MIN_SAMPLES && f.resid_p50 < reject.resid_p50) reject = f; continue; }
            c[0] = f.c[0]; c[1] = f.c[1];
            size_t k = 0;
            for (k = 0; k < nc; k++) if (hypot(cl[k].y - c[0], cl[k].x - c[1]) <= 30.0) break;
            if (k == nc) { if (nc >= 64) continue; cl[nc].y = c[0]; cl[nc].x = c[1]; cl[nc].resid = f.resid_p50; cl[nc].void_ratio = f.void_ratio; cl[nc].members = 0; nc++; }
            cl[k].members++;
            if (f.resid_p50 < cl[k].resid) { cl[k].resid = f.resid_p50; cl[k].void_ratio = f.void_ratio; cl[k].y = c[0]; cl[k].x = c[1]; }
        }
        int best = -1, best_any = -1;
        for (size_t k = 0; k < nc; k++) {
            if (best_any < 0 || cl[k].resid < cl[(size_t)best_any].resid) best_any = (int)k;
            if (cl[k].void_ratio >= ASM_AXIS_DERIVE_VOID_MAX) continue;
            if (best < 0 || cl[k].resid < cl[(size_t)best].resid) best = (int)k;
        }
        fprintf(stderr, "  [axis seed] s %.0f..%.0f: %zu samples, %zu of %zu starts converged, %zu attractors; %s\n", z0, z1, m, converged, starts, nc,
                best >= 0 ? "" : (best_any >= 0 ? "none in a void" : "none locked"));
        if (nc == 0 && reject.resid_p50 < 1e299)
            fprintf(stderr, "    best rejected fit (%.0f, %.0f): %zu in the window, cond %.3f, residual p50 %.1f, void %.2f\n", reject.c[0], reject.c[1], reject.n, reject.cond, reject.resid_p50, reject.void_ratio);
        for (size_t k = 0; k < nc; k++)
            fprintf(stderr, "    attractor (%.0f, %.0f): %zu starts, residual p50 %.1f, void %.2f%s\n", cl[k].y, cl[k].x, cl[k].members, cl[k].resid, cl[k].void_ratio, (int)k == best ? " <- the seed" : "");
        if (best < 0) continue;
        double sm = 0.5 * (z0 + z1), ca = cl[(size_t)best].y, cb = cl[(size_t)best].x;
        for (int k = 0; k < 3; k++) rows[nrow*3+k] = c0[k] + sm*dir[k] + ca*e1[k] + cb*e2[k];
        nrow++;
    }
    Arena_restore(scratch, mark);
    rep->seed_slabs = nrow;
    if (nrow == 0) return NULL;
    if (nrow == 1) {
        double pt[3] = { rows[0], rows[1], rows[2] };
        return AsmAxis_line(arena, pt, dir, pt[0] - (zhi - zlo), pt[0] + (zhi - zlo));
    }
    return AsmAxis_from_rows(arena, rows, NULL, nrow);
}

/* ---- the derivation ---------------------------------------------------------------------------- */

AsmAxis *AsmAxis_derive(Arena_T arena, const float *p, const float *nrm, const float *w, size_t n,
                        const AsmAxisDeriveOpts *o, AsmAxisDeriveReport *rep)
{
    double t0 = ves_clock_sec();
    memset(rep, 0, sizeof *rep);
    rep->prior_rms = -1.0;
    if (!p || !nrm || !w || n < (size_t)ASM_AXIS_DERIVE_MIN_SAMPLES || !o || o->slab_vox <= 0.0 || o->iters < 1) return NULL;
    Arena_T scratch = Arena_new();
    AsmAxis *axis = NULL;
    /* the direction the normals alone give, before any centre exists: the seed's slab normal and
     * an independent check on whatever curve comes out */
    double ndir[3] = { 1.0, 0.0, 0.0 };
    if (AsmAxis_normal_direction(nrm, w, n, ndir, &rep->normal_cond, &rep->normal_frac) == 0) {
        for (int k = 0; k < 3; k++) rep->normal_dir[k] = ndir[k];
        rep->normal_tilt_deg = acos(ndir[0] > 1.0 ? 1.0 : ndir[0]) * 180.0 / AD_PI;
    } else { rep->normal_dir[0] = 1.0; rep->normal_cond = 1.0; }
    if (o->prior) { axis = (AsmAxis *)o->prior; rep->seeded = 0; }
    else {
        axis = ad_seed(arena, scratch, p, nrm, w, n, ndir, rep);
        rep->seeded = axis ? 1 : -1;
        if (!axis) {
            /* the normals' direction through the samples' centroid: an honest last resort, reported as such */
            double c[3] = { 0.0, 0.0, 0.0 }, zlo = 1e300, zhi = -1e300;
            for (size_t i = 0; i < n; i++) { for (int k = 0; k < 3; k++) c[k] += p[i*3+k]; if (p[i*3] < zlo) zlo = p[i*3]; if (p[i*3] > zhi) zhi = p[i*3]; }
            for (int k = 0; k < 3; k++) c[k] /= (double)n;
            axis = AsmAxis_line(arena, c, ndir, zlo - 128.0, zhi + 128.0);
        }
    }
    const AsmAxis *prior = o->prior;
    rep->samples = n;
    double *s = ARENA_ALLOC(scratch, n * sizeof(double));
    int32_t *slab_of = ARENA_ALLOC(scratch, n * sizeof(int32_t));
    AsmAxisDeriveRow *rows = NULL; size_t nrows = 0;
    for (int it = 0; it < o->iters; it++) {
        Arena_Mark mark = Arena_save(scratch);
        double smin = 1e300, smax = -1e300;
        for (size_t i = 0; i < n; i++) { double q[3] = { p[i*3], p[i*3+1], p[i*3+2] }; ad_project(axis, q, &s[i], NULL, NULL); if (s[i] < smin) smin = s[i]; if (s[i] > smax) smax = s[i]; }
        size_t nslab = (size_t)((smax - smin) / o->slab_vox) + 1;
        size_t *count = ARENA_CALLOC(scratch, nslab + 1, sizeof(size_t));
        for (size_t i = 0; i < n; i++) { size_t k = (size_t)((s[i] - smin) / o->slab_vox); if (k >= nslab) k = nslab - 1; slab_of[i] = (int32_t)k; count[k+1]++; }
        for (size_t k = 0; k < nslab; k++) count[k+1] += count[k];
        size_t *order = ARENA_ALLOC(scratch, n * sizeof(size_t)), *fill = ARENA_CALLOC(scratch, nslab, sizeof(size_t));
        for (size_t i = 0; i < n; i++) order[count[(size_t)slab_of[i]] + fill[(size_t)slab_of[i]]++] = i;
        AsmAxisDeriveRow *rws = ARENA_CALLOC(arena, nslab, sizeof(AsmAxisDeriveRow));
        double *lock_rows = ARENA_ALLOC(arena, nslab * 3 * sizeof(double)), *lock_w = ARENA_ALLOC(arena, nslab * sizeof(double));
        size_t nlock = 0;
        size_t maxm = 0; for (size_t k = 0; k < nslab; k++) if (count[k+1] - count[k] > maxm) maxm = count[k+1] - count[k];
        double *pa = ARENA_ALLOC(scratch, (maxm ? maxm : 1) * sizeof(double)), *pb = ARENA_ALLOC(scratch, (maxm ? maxm : 1) * sizeof(double));
        double *ma = ARENA_ALLOC(scratch, (maxm ? maxm : 1) * sizeof(double)), *mb = ARENA_ALLOC(scratch, (maxm ? maxm : 1) * sizeof(double));
        double *ww = ARENA_ALLOC(scratch, (maxm ? maxm : 1) * sizeof(double)), *pairs = ARENA_ALLOC(scratch, 2 * (maxm ? maxm : 1) * sizeof(double));
        double *sa = ARENA_ALLOC(scratch, (maxm ? maxm : 1) * sizeof(double)), *sb = ARENA_ALLOC(scratch, (maxm ? maxm : 1) * sizeof(double));
        double *sma = ARENA_ALLOC(scratch, (maxm ? maxm : 1) * sizeof(double)), *smb = ARENA_ALLOC(scratch, (maxm ? maxm : 1) * sizeof(double)), *sw = ARENA_ALLOC(scratch, (maxm ? maxm : 1) * sizeof(double));
        for (size_t k = 0; k < nslab; k++) {
            AsmAxisDeriveRow *row = &rws[k];
            row->s = smin + ((double)k + 0.5) * o->slab_vox;
            double c0[3], t[3], e1[3], e2[3];
            ad_frame_at(axis, row->s, c0, t, e1, e2);
            row->z = c0[0]; row->y = c0[1]; row->x = c0[2];
            size_t m = 0;
            for (size_t q = count[k]; q < count[k+1]; q++) {
                size_t i = order[q];
                double nt = nrm[i*3]*t[0] + nrm[i*3+1]*t[1] + nrm[i*3+2]*t[2];
                if (fabs(nt) > ASM_AXIS_DERIVE_MAX_AXIAL) continue;
                double d[3] = { p[i*3]-c0[0], p[i*3+1]-c0[1], p[i*3+2]-c0[2] };
                double na = nrm[i*3]*e1[0] + nrm[i*3+1]*e1[1] + nrm[i*3+2]*e1[2], nb = nrm[i*3]*e2[0] + nrm[i*3+1]*e2[1] + nrm[i*3+2]*e2[2];
                double nl = sqrt(na*na + nb*nb);
                if (nl < 0.3) continue;
                pa[m] = d[0]*e1[0] + d[1]*e1[1] + d[2]*e1[2]; pb[m] = d[0]*e2[0] + d[1]*e2[1] + d[2]*e2[2];
                ma[m] = na / nl; mb[m] = nb / nl; ww[m] = w[i]; m++;
            }
            row->n = m;
            if (m < (size_t)ASM_AXIS_DERIVE_MIN_SAMPLES) continue;
            double czero[2] = { 0.0, 0.0 }; AdFit f;
            int fitted = ad_fit_windowed(pa, pb, ma, mb, ww, m, czero, ASM_AXIS_DERIVE_LOCAL_R, &f, sa, sb, sma, smb, sw, pairs) == 0;
            row->n = f.n ? f.n : m;   /* the samples inside the final window */
            if (!fitted) continue;
            row->kept = f.kept; row->cond = f.cond; row->aperture_deg = f.aperture_deg; row->resid_p50 = f.resid_p50;
            row->unc_major = f.unc_major; row->unc_minor = f.unc_minor; row->void_ratio = f.void_ratio; row->sectors = f.sectors; row->one_sided = f.one_sided;
            row->z = c0[0] + f.c[0]*e1[0] + f.c[1]*e2[0]; row->y = c0[1] + f.c[0]*e1[1] + f.c[1]*e2[1]; row->x = c0[2] + f.c[0]*e1[2] + f.c[1]*e2[2];
            if (!f.ok) continue;
            row->status = 1;
            lock_rows[nlock*3] = row->z; lock_rows[nlock*3+1] = row->y; lock_rows[nlock*3+2] = row->x;
            double u = f.unc_major > 1e-3 ? f.unc_major : 1e-3;
            lock_w[nlock] = 1.0 / (u * u);
            nlock++;
        }
        rows = rws; nrows = nslab; rep->iters = it + 1;
        rep->n_lock = nlock;
        if (nlock < 2) { Arena_restore(scratch, mark); break; }
        AsmAxis *next = AsmAxis_from_rows(arena, lock_rows, lock_w, nlock);
        if (!next) { Arena_restore(scratch, mark); break; }
        /* how far the axis moved: the locked centres against the previous polyline */
        double move = 0.0;
        for (size_t k = 0; k < nlock; k++) { double q[3] = { lock_rows[k*3], lock_rows[k*3+1], lock_rows[k*3+2] }, rr = 0.0; AsmAxis_project(axis, q, NULL, &rr, NULL, NULL); if (rr > move) move = rr; }
        rep->last_move = move;
        axis = next;
        Arena_restore(scratch, mark);
        if (move < 1.0) break;
    }
    rep->rows = rows; rep->n_rows = nrows;
    if (rep->n_lock < 2) { Arena_dispose(&scratch); rep->sec = ves_clock_sec() - t0; return NULL; }
    /* the report's summary */
    {
        double zlo = 1e300, zhi = -1e300;
        for (size_t i = 0; i < n; i++) { if (p[i*3] < zlo) zlo = p[i*3]; if (p[i*3] > zhi) zhi = p[i*3]; }
        double dir[3]; AsmAxis_mean_dir(axis, zlo, zhi, dir);
        double cz = dir[0] > 1.0 ? 1.0 : (dir[0] < -1.0 ? -1.0 : dir[0]);
        rep->tilt_deg = acos(cz) * 180.0 / AD_PI;
        double dot = fabs(dir[0]*rep->normal_dir[0] + dir[1]*rep->normal_dir[1] + dir[2]*rep->normal_dir[2]);
        rep->normal_disagree_deg = acos(dot > 1.0 ? 1.0 : dot) * 180.0 / AD_PI;
        rep->fit_rms = axis->fit_rms; rep->r_curv_min = axis->r_curv_min;
        rep->cond_min = 1e300; rep->aperture_min_deg = 1e300; rep->resid_p50_max = 0.0;
        double pr = 0.0; size_t np = 0;
        for (size_t k = 0; k < nrows; k++) {
            const AsmAxisDeriveRow *row = &rows[k];
            if (!row->status) continue;
            if (row->one_sided) rep->one_sided++;
            if (row->cond < rep->cond_min) rep->cond_min = row->cond;
            if (row->aperture_deg < rep->aperture_min_deg) rep->aperture_min_deg = row->aperture_deg;
            if (row->resid_p50 > rep->resid_p50_max) rep->resid_p50_max = row->resid_p50;
            if (prior) { double q[3] = { row->z, row->y, row->x }, rr = 0.0; AsmAxis_project(prior, q, NULL, &rr, NULL, NULL); pr += rr * rr; np++; }
        }
        if (np) rep->prior_rms = sqrt(pr / (double)np);
    }
    rep->accepted = rep->n_rows > 0 && (double)rep->n_lock >= ASM_AXIS_DERIVE_ACCEPT_LOCK * (double)rep->n_rows && rep->fit_rms <= ASM_AXIS_DERIVE_ACCEPT_RMS;
    Arena_dispose(&scratch);
    rep->sec = ves_clock_sec() - t0;
    return axis;
}

AsmAxis *AsmAxis_derive_run(Arena_T arena, const AsmRun *run, const AsmAxisDeriveOpts *o, AsmAxisDeriveReport *rep)
{
    size_t cap = 0;
    for (size_t c = 0; c < run->n_charts; c++) if (AsmChart_in_layout(&run->charts[c])) cap += run->charts[c].nf < o->per_chart ? run->charts[c].nf : o->per_chart;
    Arena_T scratch = Arena_new();
    float *p = ARENA_ALLOC(scratch, (cap ? cap : 1) * 3 * sizeof(float)), *nrm = ARENA_ALLOC(scratch, (cap ? cap : 1) * 3 * sizeof(float)), *w = ARENA_ALLOC(scratch, (cap ? cap : 1) * sizeof(float));
    size_t n = 0;
    for (size_t c = 0; c < run->n_charts; c++) {
        const AsmChart *ch = &run->charts[c];
        if (!AsmChart_in_layout(ch) || !ch->nf || !ch->xyz || !ch->faces) continue;
        size_t stride = ch->nf > o->per_chart ? ch->nf / o->per_chart + 1 : 1;
        for (size_t f = 0; f < ch->nf && n < cap; f += stride) {
            const int32_t *fv = &ch->faces[f*3];
            const float *a = &ch->xyz[(size_t)fv[0]*3], *b = &ch->xyz[(size_t)fv[1]*3], *d = &ch->xyz[(size_t)fv[2]*3];
            double e1[3] = { b[0]-a[0], b[1]-a[1], b[2]-a[2] }, e2[3] = { d[0]-a[0], d[1]-a[1], d[2]-a[2] };
            double cr[3] = { e1[1]*e2[2]-e1[2]*e2[1], e1[2]*e2[0]-e1[0]*e2[2], e1[0]*e2[1]-e1[1]*e2[0] };
            double L = sqrt(cr[0]*cr[0] + cr[1]*cr[1] + cr[2]*cr[2]);
            if (L < 1e-9) continue;
            for (int k = 0; k < 3; k++) { p[n*3+k] = (float)((a[k] + b[k] + d[k]) / 3.0); nrm[n*3+k] = (float)(cr[k] / L); }
            w[n] = (float)(0.5 * L * (double)stride);
            n++;
        }
    }
    AsmAxis *axis = AsmAxis_derive(arena, p, nrm, w, n, o, rep);
    Arena_dispose(&scratch);
    return axis;
}

/* ---- ledger and picture ------------------------------------------------------------------------- */

int AsmAxisDerive_write_csv(const AsmAxisDeriveReport *rep, const char *path)
{
    FILE *fp = fopen(path, "wb");
    if (!fp) return -1;
    fprintf(fp, "# derived axis (asm_axis_derive.c): slabs %zu locked %zu seeded %d iters %d last_move %.2f tilt %.2f fit_rms %.2f prior_rms %.2f one_sided %zu samples %zu\n",
            rep->n_rows, rep->n_lock, rep->seeded, rep->iters, rep->last_move, rep->tilt_deg, rep->fit_rms, rep->prior_rms, rep->one_sided, rep->samples);
    fprintf(fp, "# z,y,x,status,s,n,kept,cond,aperture_deg,resid_p50,unc_major,unc_minor,sectors,void_ratio,one_sided   (abstained rows are commented out: only locked rows load)\n");
    for (size_t k = 0; k < rep->n_rows; k++) {
        const AsmAxisDeriveRow *r = &rep->rows[k];
        fprintf(fp, "%s%.2f,%.2f,%.2f,%d,%.1f,%zu,%zu,%.4f,%.1f,%.2f,%.2f,%.2f,%d,%.3f,%d\n", r->status ? "" : "# ",
                r->z, r->y, r->x, r->status, r->s, r->n, r->kept, r->cond, r->aperture_deg, r->resid_p50, r->unc_major, r->unc_minor, r->sectors, r->void_ratio, r->one_sided);
    }
    return fclose(fp) ? -1 : 0;
}

int AsmAxisDerive_write_png(Arena_T arena, const AsmAxis *derived, const AsmAxis *prior, const AsmAxisDeriveReport *rep, const char *path)
{
    if (!derived || !rep || !rep->n_rows) return -1;
    Arena_Mark mark = Arena_save(arena);
    const int W = 1200, H = 600, PAD = 40;
    AsmCanvas cv = AsmCanvas_new(arena, W, H, 24, 24, 28);
    double zlo = 1e300, zhi = -1e300, lo[2] = { 1e300, 1e300 }, hi[2] = { -1e300, -1e300 };
    for (size_t k = 0; k < rep->n_rows; k++) {
        const AsmAxisDeriveRow *r = &rep->rows[k];
        if (r->z < zlo) zlo = r->z; if (r->z > zhi) zhi = r->z;
        if (!r->status) continue;
        if (r->y < lo[0]) lo[0] = r->y; if (r->y > hi[0]) hi[0] = r->y;
        if (r->x < lo[1]) lo[1] = r->x; if (r->x > hi[1]) hi[1] = r->x;
    }
    if (!(zhi > zlo) || !(hi[0] >= lo[0]) || !(hi[1] >= lo[1])) { Arena_restore(arena, mark); return -1; }
    for (int c = 0; c < 2; c++) { double m = 0.1 * (hi[c] - lo[c]) + 20.0; lo[c] -= m; hi[c] += m; }
    const uint8_t grey[3] = { 120, 120, 120 }, green[3] = { 80, 220, 120 }, red[3] = { 240, 80, 80 }, white[3] = { 230, 230, 230 };
    for (int c = 0; c < 2; c++) {
        int y0 = PAD + c * (H / 2), y1 = y0 + H / 2 - 2 * PAD;
        /* the prior and the derived curves, sampled along z */
        for (int pass = 0; pass < 2; pass++) {
            const AsmAxis *a = pass ? derived : prior;
            if (!a) continue;
            double px = -1.0, py = -1.0;
            for (int i = 0; i <= 400; i++) {
                double z = zlo + (zhi - zlo) * (double)i / 400.0, q[3] = { z, 0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]) }, foot[3];
                AsmAxis_project(a, q, NULL, NULL, NULL, foot);
                double v = foot[1 + c];
                double X = PAD + (W - 2 * PAD) * (z - zlo) / (zhi - zlo), Y = y1 - (y1 - y0) * (v - lo[c]) / (hi[c] - lo[c]);
                if (px >= 0.0) AsmCanvas_line(&cv, px, py, X, Y, pass ? green : grey);
                px = X; py = Y;
            }
        }
        for (size_t k = 0; k < rep->n_rows; k++) {
            const AsmAxisDeriveRow *r = &rep->rows[k];
            double X = PAD + (W - 2 * PAD) * (r->z - zlo) / (zhi - zlo);
            if (!r->status) { AsmCanvas_rect(&cv, X - 1, y1 + 2, X + 1, y1 + 8, red); continue; }
            double v = c ? r->x : r->y, Y = y1 - (y1 - y0) * (v - lo[c]) / (hi[c] - lo[c]);
            AsmCanvas_rect(&cv, X - 2, Y - 2, X + 2, Y + 2, r->one_sided ? red : white);
        }
        AsmCanvas_line(&cv, PAD, y1, W - PAD, y1, grey);
    }
    int rc = AsmCanvas_write(&cv, path);
    Arena_restore(arena, mark);
    return rc;
}

/* ---- selftest --------------------------------------------------------------------------------- */

static uint32_t ad_rng(uint32_t *s) { *s = *s * 1664525u + 1013904223u; return *s; }
static double ad_unit(uint32_t *s) { return (double)(ad_rng(s) % 20001) / 10000.0 - 1.0; }

/* Face samples of an Archimedean spiral (r0 + pitch * theta / 2pi, `turns` turns, theta covering
 * [0, cover_deg) of every turn, ellipticity `ellip` along e1) about the straight axis through
 * `origin` along the unit `dir`, over s in [0, length) every `ds`.  Normals are radial.  Returns n. */
static size_t ad_synth(float *p, float *nrm, float *w, size_t cap, const double origin[3], const double dir[3],
                       double r0, double pitch, int turns, double cover_deg, double ellip, double length, double ds, double jitter, uint32_t *rng)
{
    double e1[3] = { -dir[1]*dir[0], 1.0 - dir[1]*dir[1], -dir[1]*dir[2] };
    double l1 = sqrt(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]);
    for (int k = 0; k < 3; k++) e1[k] /= l1;
    double e2[3] = { dir[1]*e1[2] - dir[2]*e1[1], dir[2]*e1[0] - dir[0]*e1[2], dir[0]*e1[1] - dir[1]*e1[0] };
    size_t n = 0;
    double total = 2.0 * AD_PI * (double)turns;
    for (double s = 0.0; s < length; s += ds) {
        for (double th = 0.0; th < total && n < cap; th += 8.0 / fmax(r0 + pitch * th / (2.0 * AD_PI), 1.0)) {
            double lap = fmod(th, 2.0 * AD_PI) * 180.0 / AD_PI;
            if (lap >= cover_deg) continue;
            double r = r0 + pitch * th / (2.0 * AD_PI);
            double ca = ellip * cos(th), cb = sin(th);
            double q[3], nn[3];
            for (int k = 0; k < 3; k++) {
                q[k] = origin[k] + s * dir[k] + r * (ca * e1[k] + cb * e2[k]) + jitter * ad_unit(rng);
                nn[k] = cos(th) * e1[k] + sin(th) * e2[k];   /* the radial direction (exact for a circle) */
            }
            for (int k = 0; k < 3; k++) { p[n*3+k] = (float)q[k]; nrm[n*3+k] = (float)nn[k]; }
            w[n] = 1.0f; n++;
        }
    }
    return n;
}

/* the derived polyline against the true line over the sampled z range: the largest distance */
static double ad_truth_error(const AsmAxis *a, const double origin[3], const double dir[3], double length, double *tilt_err_deg)
{
    double worst = 0.0, sumd[3] = { 0.0, 0.0, 0.0 };
    for (int i = 0; i <= 40; i++) {
        double s = length * (double)i / 40.0, q[3] = { origin[0] + s*dir[0], origin[1] + s*dir[1], origin[2] + s*dir[2] }, rr = 0.0, d[3];
        AsmAxis_project(a, q, NULL, &rr, d, NULL);
        if (rr > worst) worst = rr;
        for (int k = 0; k < 3; k++) sumd[k] += d[k];
    }
    double L = sqrt(sumd[0]*sumd[0] + sumd[1]*sumd[1] + sumd[2]*sumd[2]);
    double cosang = L > 0.0 ? (sumd[0]*dir[0] + sumd[1]*dir[1] + sumd[2]*dir[2]) / L : 0.0;
    if (cosang > 1.0) cosang = 1.0;
    if (cosang < -1.0) cosang = -1.0;
    *tilt_err_deg = acos(fabs(cosang)) * 180.0 / AD_PI;
    return worst;
}

int AsmAxisDerive_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    uint32_t rng = 7u;
    const size_t cap = 600000;
    float *p = ARENA_ALLOC(arena, cap * 3 * sizeof(float)), *nrm = ARENA_ALLOC(arena, cap * 3 * sizeof(float)), *w = ARENA_ALLOC(arena, cap * sizeof(float));
    double tau = 25.0 * AD_PI / 180.0;
    double origin[3] = { 1000.0, 3000.0, 3000.0 }, dir[3] = { cos(tau), sin(tau), 0.0 };
    AsmAxisDeriveOpts o; AsmAxisDerive_default_opts(&o);
    /* 0. the direction from the normals alone: exact on a spiral, and a FLAT pile of parallel
     * sheets must report no separation so the caller knows the direction is meaningless */
    {
        size_t n = ad_synth(p, nrm, w, cap, origin, dir, 60.0, 12.0, 30, 360.0, 1.0, 1600.0, 32.0, 1.0, &rng);
        double got[3], cond = 0.0, frac = 0.0;
        int rc = AsmAxis_normal_direction(nrm, w, n, got, &cond, &frac);
        double dot = fabs(got[0]*dir[0] + got[1]*dir[1] + got[2]*dir[2]);
        double err = acos(dot > 1.0 ? 1.0 : dot) * 180.0 / AD_PI;
        int ok = rc == 0 && err <= 0.5 && cond <= 0.05 && frac >= 0.95;
        fprintf(stderr, "  axis direction from the normals, spiral: %s (%.3f deg from the truth, separation %.4f, %.0f%% perpendicular)\n", ok ? "ok" : "FAIL", err, cond, 100.0 * frac);
        fails += !ok;
        /* parallel flat sheets: every normal is the same, so the null space is a PLANE, not a line */
        for (size_t i = 0; i < 4096; i++) {
            p[i*3] = (float)(i % 64); p[i*3+1] = (float)(i / 64); p[i*3+2] = (float)(12 * (i % 7));
            nrm[i*3] = 0.0f; nrm[i*3+1] = 0.0f; nrm[i*3+2] = 1.0f; w[i] = 1.0f;
        }
        double flat[3], fcond = 0.0;
        int rc2 = AsmAxis_normal_direction(nrm, w, 4096, flat, &fcond, NULL);
        int ok2 = rc2 == 0 && fcond > 0.5;
        fprintf(stderr, "  axis direction from the normals, parallel sheets: %s (separation %.4f, no single direction)\n", ok2 ? "ok" : "FAIL", fcond);
        fails += !ok2;
    }
    /* 1. a tilted circular spiral, full aperture, from a prior 40 vox off and 10 degrees wrong */
    {
        size_t n = ad_synth(p, nrm, w, cap, origin, dir, 60.0, 12.0, 30, 360.0, 1.0, 1600.0, 32.0, 1.0, &rng);
        double bad_dir[3] = { cos(tau + 10.0 * AD_PI / 180.0), sin(tau + 10.0 * AD_PI / 180.0), 0.0 }, bad_origin[3] = { origin[0], origin[1] + 40.0, origin[2] - 20.0 };
        AsmAxis *prior = AsmAxis_line(arena, bad_origin, bad_dir, 900.0, 2600.0);
        o.prior = prior; AsmAxisDeriveReport rep;
        AsmAxis *a = AsmAxis_derive(arena, p, nrm, w, n, &o, &rep);
        double tilt_err = 0.0, err = a ? ad_truth_error(a, origin, dir, 1600.0, &tilt_err) : 1e300;
        int ok = a && err <= 3.0 && tilt_err <= 0.5 && rep.n_lock >= 20 && rep.one_sided == 0;
        fprintf(stderr, "  derived axis, tilted spiral from a wrong prior: %s (%zu samples, %zu/%zu slabs locked, worst %.2f vox, tilt error %.2f deg, cond min %.3f, resid max %.2f, prior rms %.1f, %.1f s)\n",
                ok ? "ok" : "FAIL", n, rep.n_lock, rep.n_rows, err, tilt_err, rep.cond_min, rep.resid_p50_max, rep.prior_rms, rep.sec);
        fails += !ok;
        /* 5. the ledger round trip: the locked rows load as a table whose projections agree */
        if (a) {
            const char *path = "asm_axis_derive_selftest.csv";
            AsmAxis *back = AsmAxisDerive_write_csv(&rep, path) == 0 ? AsmAxis_load(arena, path) : NULL;
            double e2 = 0.0;
            if (back) { double tilt2; e2 = ad_truth_error(back, origin, dir, 1600.0, &tilt2); }
            int ok2 = back && back->n_rows == rep.n_lock && e2 <= 3.5;
            fprintf(stderr, "  derived axis ledger round trip: %s (%zu rows, worst %.2f vox)\n", ok2 ? "ok" : "FAIL", back ? back->n_rows : (size_t)0, e2);
            fails += !ok2;
            remove(path);
        }
    }
    /* 2. a one-sided aperture (110 degrees of every turn) still locks, flagged one-sided, within 5 vox */
    {
        size_t n = ad_synth(p, nrm, w, cap, origin, dir, 60.0, 12.0, 30, 110.0, 1.0, 1600.0, 32.0, 1.0, &rng);
        double bad_origin[3] = { origin[0], origin[1] + 30.0, origin[2] + 30.0 };
        o.prior = AsmAxis_line(arena, bad_origin, dir, 900.0, 2600.0); AsmAxisDeriveReport rep;
        AsmAxis *a = AsmAxis_derive(arena, p, nrm, w, n, &o, &rep);
        double tilt_err = 0.0, err = a ? ad_truth_error(a, origin, dir, 1600.0, &tilt_err) : 1e300;
        int ok = a && err <= 5.0 && tilt_err <= 1.0 && rep.n_lock >= 20 && rep.one_sided == rep.n_lock;
        fprintf(stderr, "  derived axis, one-sided aperture: %s (%zu/%zu locked, %zu one-sided, worst %.2f vox, tilt error %.2f deg, aperture min %.0f deg, cond min %.3f)\n",
                ok ? "ok" : "FAIL", rep.n_lock, rep.n_rows, rep.one_sided, err, tilt_err, rep.aperture_min_deg, rep.cond_min);
        fails += !ok;
    }
    /* 3. an elliptical one-sided section: the direction holds, the along-bisector uncertainty is the larger axis */
    {
        size_t n = ad_synth(p, nrm, w, cap, origin, dir, 120.0, 12.0, 20, 110.0, 0.85, 1600.0, 32.0, 1.0, &rng);
        o.prior = AsmAxis_line(arena, origin, dir, 900.0, 2600.0); AsmAxisDeriveReport rep;
        AsmAxis *a = AsmAxis_derive(arena, p, nrm, w, n, &o, &rep);
        double tilt_err = 0.0, err = a ? ad_truth_error(a, origin, dir, 1600.0, &tilt_err) : 1e300;
        size_t wider = 0;
        for (size_t k = 0; k < rep.n_rows; k++) if (rep.rows[k].status && rep.rows[k].unc_major >= rep.rows[k].unc_minor) wider++;
        int ok = a && tilt_err <= 1.5 && rep.n_lock >= 20 && wider == rep.n_lock;
        fprintf(stderr, "  derived axis, elliptical one-sided section: %s (bias %.1f vox reported, tilt error %.2f deg, unc major/minor of the first locked row %.1f/%.1f)\n",
                ok ? "ok" : "FAIL", err, tilt_err, rep.n_lock ? rep.rows[0].unc_major : 0.0, rep.n_lock ? rep.rows[0].unc_minor : 0.0);
        fails += !ok;
    }
    /* 4. no prior: the seed found on z-slabs from a grid of starts */
    {
        size_t n = ad_synth(p, nrm, w, cap, origin, dir, 60.0, 12.0, 30, 360.0, 1.0, 1600.0, 32.0, 1.0, &rng);
        o.prior = NULL; AsmAxisDeriveReport rep;
        AsmAxis *a = AsmAxis_derive(arena, p, nrm, w, n, &o, &rep);
        double tilt_err = 0.0, err = a ? ad_truth_error(a, origin, dir, 1600.0, &tilt_err) : 1e300;
        int ok = a && rep.seeded == 1 && err <= 3.0 && tilt_err <= 0.5;
        fprintf(stderr, "  derived axis, seeded without a prior: %s (seeded %d from %zu z-slabs, worst %.2f vox, tilt error %.2f deg)\n", ok ? "ok" : "FAIL", rep.seeded, rep.seed_slabs, err, tilt_err);
        fails += !ok;
    }
    Arena_dispose(&arena);
    if (fails == 0) fprintf(stderr, "  asm_axis_derive selftest: all passed\n");
    return fails;
}
