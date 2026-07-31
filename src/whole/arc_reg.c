#include "arc_reg.h"
#include "cube_register.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef ARC_PI
#define ARC_PI 3.14159265358979323846
#endif

void ArcRegOpts_default(ArcRegOpts *o)
{
    assert(o);
    o->n_knots = 5;
    o->lambda_smooth = 4.0;
    o->lambda_ridge = 0.02;
    o->cg_iters = 200;
    o->relax_iters = 5;
    o->final_iters = 2;
    o->sigma_floor = 2.0;
    o->max_warp = 40.0;
    o->min_node_corr = 8;
    o->conf_min = 1;
    o->max_mismatch = 60.0;
    o->verbose = 0;
}

double ArcReg_eval(const ArcRegWarp *w, double phi)
{
    if (w == NULL || w->n_knots <= 0) return 0.0;
    if (w->n_knots == 1 || w->dphi <= 0.0) return w->delta[0];
    double t = (phi - w->phi0) / w->dphi;
    if (t <= 0.0) return w->delta[0];
    if (t >= (double)(w->n_knots - 1)) return w->delta[w->n_knots - 1];
    int i = (int)floor(t);
    double f = t - (double)i;
    return w->delta[i] * (1.0 - f) + w->delta[i + 1] * f;
}

/* ---- (cube,gid) -> node hash ---- */
typedef struct {
    uint64_t *key; int32_t *node; size_t cap, mask;
} NodeHash;

static void nh_init(Arena_T ar, NodeHash *h, size_t hint)
{
    size_t cap = 64; while (cap < hint * 2) cap <<= 1;
    h->cap = cap; h->mask = cap - 1;
    h->key = (uint64_t *)ARENA_ALLOC(ar, (size_t)(cap * sizeof(uint64_t)));
    h->node = (int32_t *)ARENA_ALLOC(ar, (size_t)(cap * sizeof(int32_t)));
    for (size_t i = 0; i < cap; i++) { h->key[i] = UINT64_MAX; h->node[i] = -1; }
}
static int32_t nh_get_or_add(NodeHash *h, int32_t cube, int32_t gid, int32_t *pn)
{
    uint64_t k = ((uint64_t)(uint32_t)cube << 20) | (uint64_t)((uint32_t)gid & 0xFFFFFu);
    uint64_t i = (k * 1099511628211ull) & h->mask;
    for (;;) {
        if (h->key[i] == k) return h->node[i];
        if (h->key[i] == UINT64_MAX) {
            h->key[i] = k; h->node[i] = (*pn)++;
            return h->node[i];
        }
        i = (i + 1) & h->mask;
    }
}

typedef struct {
    int32_t cube, gid;
    double  phi_min, phi_max;
    int32_t n_corr;
    int32_t n_knots, base;   /* var base index in x */
    double  phi0, dphi;
} NodeInfo;

typedef struct {
    int32_t va_lo, va_hi, vb_lo, vb_hi;
    double  wa_lo, wa_hi, wb_lo, wb_hi;
    double  r0, base_weight, like;
} Spring;

/* knot interpolation of node n at phi -> (lo var, w_lo, hi var, w_hi). */
static void knot_weights(const NodeInfo *ni, double phi,
                         int32_t *lo, double *wlo, int32_t *hi, double *whi)
{
    if (ni->n_knots <= 1 || ni->dphi <= 0.0) {
        *lo = ni->base; *wlo = 1.0; *hi = ni->base; *whi = 0.0; return;
    }
    double t = (phi - ni->phi0) / ni->dphi;
    if (t <= 0.0) { *lo = ni->base; *wlo = 1.0; *hi = ni->base; *whi = 0.0; return; }
    if (t >= (double)(ni->n_knots - 1)) {
        *lo = ni->base + ni->n_knots - 1; *wlo = 1.0; *hi = *lo; *whi = 0.0; return;
    }
    int i = (int)floor(t); double f = t - (double)i;
    *lo = ni->base + i;     *wlo = 1.0 - f;
    *hi = ni->base + i + 1; *whi = f;
}

/* Apply the current IRLS normal matrix. This is the continuous analogue of
 * StrokeStrip's global step: seam candidates are cross-section pairs, `like`
 * is their current WTS-adjacency likelihood, knot differences regularize a
 * smooth correction function, and the ridge retains Ribbon's raw chart. */
static void arc_apply_normal(const Spring *springs, size_t nsp,
                             const NodeInfo *nodes, int32_t n_nodes,
                             int total_vars, const ArcRegOpts *opts,
                             const double *x, double *y)
{
    memset(y, 0, (size_t)total_vars * sizeof(double));
    for (size_t s = 0; s < nsp; s++) {
        const Spring *sp = &springs[s];
        double Ax = x[sp->va_lo]*sp->wa_lo + x[sp->va_hi]*sp->wa_hi
                  - x[sp->vb_lo]*sp->wb_lo - x[sp->vb_hi]*sp->wb_hi;
        double cc = sp->base_weight * sp->like * Ax;
        y[sp->va_lo] += sp->wa_lo * cc;
        y[sp->va_hi] += sp->wa_hi * cc;
        y[sp->vb_lo] -= sp->wb_lo * cc;
        y[sp->vb_hi] -= sp->wb_hi * cc;
    }
    for (int32_t i = 0; i < n_nodes; i++) {
        const NodeInfo *nd = &nodes[i];
        for (int kk = 0; kk + 1 < nd->n_knots; kk++) {
            double d = x[nd->base+kk+1] - x[nd->base+kk];
            double cc = opts->lambda_smooth * d;
            y[nd->base+kk+1] += cc;
            y[nd->base+kk] -= cc;
        }
    }
    for (int i = 0; i < total_vars; i++)
        y[i] += opts->lambda_ridge * x[i];
}

static int cmp_dbl(const void *a, const void *b)
{ double x = *(const double *)a, y = *(const double *)b; return x < y ? -1 : (x > y ? 1 : 0); }

int ArcReg_solve(Arena_T arena,
                 const WeldTrack *wt, const int32_t *wt_to_local,
                 SkinVert *const *skins, const size_t *nskin,
                 const PlacedReg *regs, size_t n_cubes,
                 double spiral_a, double spiral_b,
                 const ArcRegOpts *opts_in,
                 ArcRegWarp *out_warps, size_t *out_n_warps,
                 ArcRegStats *stats)
{
    assert(arena && wt && wt_to_local && skins && nskin && regs);
    ArcRegOpts opts; if (opts_in) opts = *opts_in; else ArcRegOpts_default(&opts);
    if (opts.n_knots < 1) opts.n_knots = 1;
    if (opts.n_knots > ARC_REG_MAX_KNOTS) opts.n_knots = ARC_REG_MAX_KNOTS;
    ArcRegStats z; memset(&z, 0, sizeof z);

    Arena_Mark mark = Arena_save(arena);

    /* helper to fetch a corr endpoint's registered u + raw phi + node key. */
    NodeHash nh; nh_init(arena, &nh, wt->n_corr + 64);
    NodeInfo *nodes = (NodeInfo *)ARENA_CALLOC(arena,
                        (size_t)(wt->n_corr + 1), (size_t)sizeof(NodeInfo));
    int32_t n_nodes = 0;

    /* endpoint resolution: returns 0 ok with (reg_u, phi_raw, node); -1 skip. */
    #define RESOLVE(CUBE_WT, VERT, GID_OUT, REGU_OUT, PHI_OUT, NODE_OUT) do { \
        int32_t _cl = wt_to_local[(CUBE_WT)];                                 \
        if (_cl < 0 || (size_t)_cl >= n_cubes) { ok = 0; break; }             \
        if ((size_t)(VERT) >= nskin[_cl]) { ok = 0; break; }                  \
        SkinVert _s = skins[_cl][(VERT)];                                     \
        SkinVert _r = CubeReg_apply_vert(_s, &regs[_cl], spiral_a, spiral_b); \
        (REGU_OUT) = (double)_r.u; (PHI_OUT) = (double)_s.phi;                \
        (GID_OUT) = _s.gid;                                                   \
        int32_t _nid = nh_get_or_add(&nh, _cl, _s.gid, &n_nodes);             \
        if (_nid >= (int32_t)(wt->n_corr + 1)) { ok = 0; break; }             \
        if (nodes[_nid].n_corr == 0) {          /* first touch (CALLOC'd) */  \
            nodes[_nid].cube = _cl; nodes[_nid].gid = _s.gid;                 \
            nodes[_nid].phi_min = 1e300; nodes[_nid].phi_max = -1e300;        \
        }                                                                     \
        if ((PHI_OUT) < nodes[_nid].phi_min) nodes[_nid].phi_min = (PHI_OUT); \
        if ((PHI_OUT) > nodes[_nid].phi_max) nodes[_nid].phi_max = (PHI_OUT); \
        nodes[_nid].n_corr++;                                                 \
        (NODE_OUT) = _nid;                                                    \
    } while (0)

    /* pass 1: enumerate nodes + accumulate phi ranges. */
    Spring *springs = (Spring *)ARENA_ALLOC(arena,
                        (size_t)((wt->n_corr + 1) * sizeof(Spring)));
    /* we store (node_a,node_b,phi_a,phi_b,r0,weight) transiently, then finalize
     * knot weights in pass 2 after knots are assigned. */
    typedef struct { int32_t na, nb; double pa, pb, r0, w; } RawSpring;
    RawSpring *raw = (RawSpring *)ARENA_ALLOC(arena,
                        (size_t)((wt->n_corr + 1) * sizeof(RawSpring)));
    size_t nsp = 0, n_gated = 0;

    for (size_t c = 0; c < wt->n_corr; c++) {
        const WtrkCorr *cr = &wt->corr[c];
        if ((int)cr->conf < opts.conf_min) continue;
        int ok = 1;
        int32_t gida = 0, gidb = 0, na = -1, nb = -1;
        double ua = 0, ub = 0, pa = 0, pb = 0;
        RESOLVE(cr->cube_a, cr->vert_a, gida, ua, pa, na);
        if (!ok) continue;
        RESOLVE(cr->cube_b, cr->vert_b, gidb, ub, pb, nb);
        if (!ok) continue;
        if (na == nb) continue;                 /* self-node spring: no info */
        double r0 = ua - ub;
        if (fabs(r0) > opts.max_mismatch) { n_gated++; continue; } /* turn-off */
        raw[nsp].na = na; raw[nsp].nb = nb;
        raw[nsp].pa = pa; raw[nsp].pb = pb;
        raw[nsp].r0 = r0;
        raw[nsp].w  = (double)cr->conf;         /* conf 1 or 2 */
        nsp++;
    }
    #undef RESOLVE

    if (nsp == 0 || n_nodes == 0) {
        if (stats) *stats = z;
        Arena_restore(arena, mark);
        return 0;
    }

    /* assign knots + variable offsets. */
    int total_vars = 0;
    for (int32_t i = 0; i < n_nodes; i++) {
        NodeInfo *nd = &nodes[i];
        int nk = opts.n_knots;
        if (nd->n_corr < opts.min_node_corr) nk = 1;
        if (nd->phi_max <= nd->phi_min) nk = 1;
        nd->n_knots = nk;
        nd->base = total_vars;
        total_vars += nk;
        if (nk > 1) {
            nd->phi0 = nd->phi_min;
            nd->dphi = (nd->phi_max - nd->phi_min) / (double)(nk - 1);
        } else { nd->phi0 = nd->phi_min; nd->dphi = 0.0; }
    }

    /* pass 2: flatten springs with knot weights. */
    for (size_t s = 0; s < nsp; s++) {
        knot_weights(&nodes[raw[s].na], raw[s].pa,
                     &springs[s].va_lo, &springs[s].wa_lo,
                     &springs[s].va_hi, &springs[s].wa_hi);
        knot_weights(&nodes[raw[s].nb], raw[s].pb,
                     &springs[s].vb_lo, &springs[s].wb_lo,
                     &springs[s].vb_hi, &springs[s].wb_hi);
        springs[s].r0 = raw[s].r0;
        springs[s].base_weight = raw[s].w;
        springs[s].like = 1.0;
    }

    /* ---- StrokeStrip-style alternating local/global fit ------------------
     *
     * With likelihoods fixed, CG solves
     *   (A' W A + ls S' S + lr I)x = -A' W r0.
     * The local step then updates each candidate correspondence likelihood
     * from its current parameter residual. Broad rounds retain plausible
     * large seam drift; the final rounds tighten around the inlier noise.
     * This prevents one false cross-ply nearest-neighbour from dragging an
     * entire cube group while still fitting a smooth non-constant correction. */
    double *x  = (double *)ARENA_CALLOC(arena, (size_t)total_vars, (size_t)sizeof(double));
    double *b  = (double *)ARENA_CALLOC(arena, (size_t)total_vars, (size_t)sizeof(double));
    double *r  = (double *)ARENA_ALLOC(arena, (size_t)(total_vars * sizeof(double)));
    double *p  = (double *)ARENA_ALLOC(arena, (size_t)(total_vars * sizeof(double)));
    double *Mp = (double *)ARENA_ALLOC(arena, (size_t)(total_vars * sizeof(double)));
    double *abs_res = (double *)ARENA_ALLOC(arena, nsp * sizeof(double));
    int rounds = opts.relax_iters + opts.final_iters;
    if (rounds < 1) rounds = 1;
    for (int round = 0; round < rounds; round++) {
        memset(b, 0, (size_t)total_vars * sizeof(double));
        for (size_t s = 0; s < nsp; s++) {
            const Spring *sp = &springs[s];
            double w = sp->base_weight * sp->like;
            double c = -w * sp->r0;
            b[sp->va_lo] += sp->wa_lo * c;
            b[sp->va_hi] += sp->wa_hi * c;
            b[sp->vb_lo] -= sp->wb_lo * c;
            b[sp->vb_hi] -= sp->wb_hi * c;
        }
        arc_apply_normal(springs, nsp, nodes, n_nodes, total_vars, &opts, x, Mp);
        double rr = 0.0;
        for (int i = 0; i < total_vars; i++) {
            r[i] = b[i] - Mp[i];
            p[i] = r[i];
            rr += r[i] * r[i];
        }
        double rr0 = rr;
        for (int it = 0; it < opts.cg_iters
                           && rr > 1e-14 * (rr0 + 1e-30); it++) {
            arc_apply_normal(springs, nsp, nodes, n_nodes, total_vars,
                             &opts, p, Mp);
            double pMp = 0.0;
            for (int i = 0; i < total_vars; i++) pMp += p[i] * Mp[i];
            if (pMp <= 1e-30) break;
            double alpha = rr / pMp;
            for (int i = 0; i < total_vars; i++) {
                x[i] += alpha * p[i];
                r[i] -= alpha * Mp[i];
            }
            double rr_new = 0.0;
            for (int i = 0; i < total_vars; i++) rr_new += r[i] * r[i];
            double beta = rr_new / rr;
            for (int i = 0; i < total_vars; i++) p[i] = r[i] + beta * p[i];
            rr = rr_new;
        }

        for (int i = 0; i < total_vars; i++) {
            if (x[i] >  opts.max_warp) x[i] =  opts.max_warp;
            if (x[i] < -opts.max_warp) x[i] = -opts.max_warp;
        }
        for (size_t s = 0; s < nsp; s++) {
            const Spring *sp = &springs[s];
            double dx = x[sp->va_lo]*sp->wa_lo + x[sp->va_hi]*sp->wa_hi
                      - x[sp->vb_lo]*sp->wb_lo - x[sp->vb_hi]*sp->wb_hi;
            abs_res[s] = fabs(sp->r0 + dx);
        }
        qsort(abs_res, nsp, sizeof(double), cmp_dbl);
        double med = abs_res[nsp / 2];
        double scale = round < opts.relax_iters ? 3.0 : 1.5;
        double sigma = scale * med;
        if (sigma < opts.sigma_floor) sigma = opts.sigma_floor;
        double inv2s2 = 1.0 / (2.0 * sigma * sigma);
        for (size_t s = 0; s < nsp; s++) {
            const Spring *sp = &springs[s];
            double dx = x[sp->va_lo]*sp->wa_lo + x[sp->va_hi]*sp->wa_hi
                      - x[sp->vb_lo]*sp->wb_lo - x[sp->vb_hi]*sp->wb_hi;
            double e = sp->r0 + dx;
            springs[s].like = exp(-e * e * inv2s2);
        }
    }
    z.reweight_rounds = rounds;

    /* Magnitude bound per node: delta remains a small refinement of Ribbon's
     * already-monotone dense parameter field. */
    size_t clamps = 0; double maxabs = 0.0;
    for (int32_t i = 0; i < n_nodes; i++) {
        const NodeInfo *nd = &nodes[i];
        for (int kk = 0; kk < nd->n_knots; kk++) {
            double *v = &x[nd->base + kk];
            if (*v >  opts.max_warp) { *v =  opts.max_warp; clamps++; }
            if (*v < -opts.max_warp) { *v = -opts.max_warp; clamps++; }
            if (fabs(*v) > maxabs) maxabs = fabs(*v);
        }
    }

    /* stats: RMS + p95 of |residual| before/after. */
    double *res0 = (double *)ARENA_ALLOC(arena, (size_t)(nsp * sizeof(double)));
    double *res1 = (double *)ARENA_ALLOC(arena, (size_t)(nsp * sizeof(double)));
    double s0 = 0.0, s1 = 0.0;
    for (size_t s = 0; s < nsp; s++) {
        const Spring *sp = &springs[s];
        double dx = x[sp->va_lo]*sp->wa_lo + x[sp->va_hi]*sp->wa_hi
                  - x[sp->vb_lo]*sp->wb_lo - x[sp->vb_hi]*sp->wb_hi;
        double r_after = sp->r0 + dx;
        res0[s] = fabs(sp->r0); res1[s] = fabs(r_after);
        s0 += sp->r0 * sp->r0; s1 += r_after * r_after;
    }
    z.rms_before = sqrt(s0 / (double)nsp);
    z.rms_after  = sqrt(s1 / (double)nsp);
    qsort(res0, nsp, sizeof(double), cmp_dbl);
    qsort(res1, nsp, sizeof(double), cmp_dbl);
    z.med_before = res0[nsp / 2];
    z.med_after  = res1[nsp / 2];
    z.p95_before = res0[(size_t)((double)nsp * 0.95)];
    z.p95_after  = res1[(size_t)((double)nsp * 0.95)];
    z.n_corr_used = nsp;
    z.n_gated = n_gated;
    for (size_t s = 0; s < nsp; s++)
        if (springs[s].like < 0.1) z.n_downweighted++;
    z.n_nodes = (size_t)n_nodes;
    for (int32_t i = 0; i < n_nodes; i++) if (nodes[i].n_knots == 1) z.n_scalar_nodes++;
    z.monotone_clamps = clamps;
    z.max_abs_warp = maxabs;

    /* emit warps. */
    if (out_warps != NULL && out_n_warps != NULL) {
        size_t cap = *out_n_warps, w = 0;
        for (int32_t i = 0; i < n_nodes && w < cap; i++) {
            const NodeInfo *nd = &nodes[i];
            ArcRegWarp *ow = &out_warps[w++];
            ow->cube = nd->cube; ow->gid = nd->gid;
            ow->n_knots = nd->n_knots; ow->phi0 = nd->phi0; ow->dphi = nd->dphi;
            for (int kk = 0; kk < nd->n_knots; kk++) ow->delta[kk] = x[nd->base + kk];
        }
        *out_n_warps = w;
    }

    if (stats) *stats = z;
    Arena_restore(arena, mark);
    return 0;
}

/* ------------------------------------------------------------------ selftest */

int ArcReg_selftest(void)
{
    int fails = 0;
    Arena_T ar = Arena_new();
    double a = 0.85, b = -9.5;

    /* Two cubes, one group each; a phi-LINEAR arc-length drift f(phi)=alpha+
     * beta*phi between them that a constant du cannot remove. Build skins +
     * a WeldTrack of same-papyrus corrs, verify the warp drives the mismatch
     * toward zero. */
    const int N = 60;
    SkinVert *sk0 = (SkinVert *)ARENA_ALLOC(ar, (size_t)(N * sizeof(SkinVert)));
    SkinVert *sk1 = (SkinVert *)ARENA_ALLOC(ar, (size_t)(N * sizeof(SkinVert)));
    double alpha = 3.0, beta = 1.7;
    for (int i = 0; i < N; i++) {
        double phi = -(double)i * 0.05 - 1.0;   /* spread in phi */
        sk0[i].p[0]=sk0[i].p[1]=sk0[i].p[2]=0;
        sk0[i].u = (float)(100.0 + 0.3*(double)i);
        sk0[i].v = (float)i; sk0[i].phi = (float)phi; sk0[i].gid = 0;
        /* cube1's copy of the same papyrus, but u drifts by f(phi). */
        sk1[i] = sk0[i];
        sk1[i].u = (float)((double)sk0[i].u + alpha + beta * phi);
    }
    SkinVert *skins[2] = { sk0, sk1 };
    size_t nskin[2] = { (size_t)N, (size_t)N };
    /* identity regs (k=0, du=0): registered u == raw u here. */
    PlacedReg regs[2]; memset(regs, 0, sizeof regs);

    WeldTrack wt; memset(&wt, 0, sizeof wt);
    wt.version = 1; wt.spiral_a = a; wt.spiral_b = b; wt.pitch = 9.5;
    wt.n_cubes = 2;
    wt.ids = (char (*)[48])ARENA_ALLOC(ar, (size_t)(2 * 48));
    memcpy(wt.ids[0], "cubeA", 6); memcpy(wt.ids[1], "cubeB", 6);
    wt.n_corr = (size_t)N;
    wt.corr = (WtrkCorr *)ARENA_ALLOC(ar, (size_t)((size_t)N * sizeof(WtrkCorr)));
    for (int i = 0; i < N; i++) {
        WtrkCorr *c = &wt.corr[i];
        memset(c, 0, sizeof *c);
        c->cube_a = 0; c->cube_b = 1; c->vert_a = (uint32_t)i; c->vert_b = (uint32_t)i;
        c->gid_a = 0; c->gid_b = 0; c->seam_axis = 2; c->conf = 2;
    }
    int32_t wt2local[2] = { 0, 1 };

    ArcRegOpts o; ArcRegOpts_default(&o); o.n_knots = 5; o.lambda_ridge = 0.01;
    ArcRegWarp warps[8]; size_t nw = 8;
    ArcRegStats st;
    ArcReg_solve(ar, &wt, wt2local, skins, nskin, regs, 2, a, b, &o,
                 warps, &nw, &st);

    if (st.n_corr_used != (size_t)N) { fprintf(stderr, "[arc_reg] used %zu want %d\n", st.n_corr_used, N); fails++; }
    if (st.n_nodes != 2) { fprintf(stderr, "[arc_reg] nodes %zu want 2\n", st.n_nodes); fails++; }
    /* the linear drift is exactly representable by piecewise-linear warps ->
     * mismatch should collapse. */
    if (!(st.rms_after < 0.05 * st.rms_before)) {
        fprintf(stderr, "[arc_reg] rms %.4f -> %.4f (want >20x drop)\n",
                st.rms_before, st.rms_after); fails++;
    }
    if (st.rms_before < 1.0) { fprintf(stderr, "[arc_reg] rms_before %.3f too small\n", st.rms_before); fails++; }

    /* eval sanity: a warp evaluated at its phi0 returns delta[0]. */
    if (nw >= 1) {
        double e = ArcReg_eval(&warps[0], warps[0].phi0);
        if (fabs(e - warps[0].delta[0]) > 1e-9) { fprintf(stderr, "[arc_reg] eval endpoint\n"); fails++; }
    }

    /* scalar fallback: a node with < min_node_corr springs gets 1 knot. */
    {
        WeldTrack wt2 = wt; wt2.n_corr = 3;   /* only 3 corrs => scalar nodes */
        ArcRegStats st2; ArcRegWarp w2[8]; size_t nw2 = 8;
        ArcReg_solve(ar, &wt2, wt2local, skins, nskin, regs, 2, a, b, &o,
                     w2, &nw2, &st2);
        if (st2.n_scalar_nodes != 2) { fprintf(stderr, "[arc_reg] scalar nodes %zu want 2\n", st2.n_scalar_nodes); fails++; }
    }

    Arena_dispose(&ar);
    if (fails == 0) fprintf(stderr, "[arc_reg] selftest PASSED\n");
    else fprintf(stderr, "[arc_reg] selftest FAILED (%d)\n", fails);
    return fails;
}
