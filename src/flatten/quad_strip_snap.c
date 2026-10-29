/* ============================================================================
 * quad_strip_snap.c -- see quad_strip_snap.h.  Progressive relaxation of the
 * solid quad strip onto the real papyrus: one soft iterated energy combining a
 * CT ridge snap, an init/data anchor, slice-to-slice continuity, and a QUASI-
 * developable ARAP prior.  Movable = unresolved fill cells + fitted sites proven
 * off-ridge with a stronger reachable target; later active-set passes may softly
 * anchor fill cells already recovered onto a credible recto transition.  Reuses
 * the RAW sampler (raw_sample.h), the geometric vertex normals (mesh_normals.h),
 * and the symmetric 3x3 eigensolver (eig3.h).
 * ==========================================================================*/

#include "quad_strip_snap.h"

#include "../common/raw_sample.h"
#include "../common/mesh_normals.h"
#include "../common/csr.h"
#include "../common/eig3.h"
#include "../remesh/intersection_cleanup.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void QuadStripSnap_defaults(QuadStripSnapOpts *o) {
    memset(o, 0, sizeof *o);
    o->reach = 6.0;        /* < wrap pitch (~9.5) so a snap can't cross a wrap */
    o->step = 0.5;
    o->bright_min = 0.0;   /* auto: fitted-cell p1 + 0.30*(p99-p1) */
    o->band_frac = 0.30;
    o->min_gain = 20.0;
    o->ridge_lock = 0.75;
    o->dark_probe_frac = 0.20;
    o->guided_reach = 12.0;
    o->winding_tube = 3.0;
    o->patch_recover = 0;
    o->patch_buffer = 2;
    o->patch_min_cells = 1;
    o->patch_ray_reach = 48.0;
    o->patch_ray_margin = 1.0;
    o->patch_relax_sweeps = 64;
    o->patch_target_weight = 16.0;
    o->patch_fitted_dark_slack = 512;
    o->local_contrast = 1;
    o->snap_side = QUAD_SNAP_SIDE_NEAREST;
    o->preserve_axial = 0;
    o->prevent_intersection_growth = 1;
    o->filled_target_anchor = 0.0;
    o->fitted_source_anchor = 0.0;
    o->axis_y = o->axis_x = 0.0;
    o->quilt_smooth = 4.0;
    o->quilt_seam_smooth = 1.0;
    o->quilt_anchor = 0.01;
    o->quilt_huber = 1.5;
    o->quilt_conf_min = 0.02;
    o->quilt_sweeps = 4;
    o->w_target = 1.0;
    o->w_orig = 0.25;
    o->w_smooth = 0.3;
    o->w_dev = 0.4;
    o->rounds = 8;
    o->solve_sweeps = 2;
    o->solve_omega = 0.5;
    o->multigrid_levels = 1;
    o->multigrid_cycles = 1;
    o->multigrid_patch = 16;
    o->bake_window_low = 47.0;
    o->bake_window_high = 193.0;
    o->bake_dark_u8 = 13;
    o->remesh_motion = 0.5;
    o->remesh_motion_max = 4.0;
    o->remesh_sander = 2.0;
    o->remesh_over4 = 0.0001;
    o->require_complete_raw = 1;
    o->cell_lit_accum = NULL;
    o->cell_candidate_accum = NULL;
    o->verbose = 0;
}

/* True for either oriented target policy (recto face or band midline): the
 * profile is oriented inward and the winding machinery is armed. */
static int snap_oriented(const QuadStripSnapOpts *o) {
    return o->snap_side != QUAD_SNAP_SIDE_NEAREST;
}

/* Closest rotation to a 3x3 matrix S (polar factor) via the symmetric eigen-
 * decomposition of M = S^T S:  R = S * M^-1/2, with a determinant flip so the
 * result is a proper rotation.  (Same helper as quad_strip.c's ARAP.) */
static void closest_rotation(const double S[9], double R[9]) {
    double M[3][3];
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) {
            double s = 0.0;
            for (int k = 0; k < 3; k++) s += S[k * 3 + a] * S[k * 3 + b];
            M[a][b] = s;
        }
    double ev[3], Q[3][3];
    Eig3_sym(M, ev, Q);
    double is[3];
    for (int k = 0; k < 3; k++) is[k] = ev[k] > 1e-18 ? 1.0 / sqrt(ev[k]) : 0.0;
    double Ms[3][3];
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) {
            double s = 0.0;
            for (int k = 0; k < 3; k++) s += Q[a][k] * is[k] * Q[b][k];
            Ms[a][b] = s;
        }
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) {
            double s = 0.0;
            for (int k = 0; k < 3; k++) s += S[a * 3 + k] * Ms[k][b];
            R[a * 3 + b] = s;
        }
    double det = R[0]*(R[4]*R[8]-R[5]*R[7]) - R[1]*(R[3]*R[8]-R[5]*R[6]) + R[2]*(R[3]*R[7]-R[4]*R[6]);
    if (det < 0.0) {
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++) Ms[a][b] -= 2.0 * is[0] * Q[a][0] * Q[b][0];
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++) {
                double s = 0.0;
                for (int k = 0; k < 3; k++) s += S[a * 3 + k] * Ms[k][b];
                R[a * 3 + b] = s;
            }
    }
}

/* The solid ribbon is a structured quad lattice.  Its position/ARAP energy must
 * use the four axial lattice edges, not the six edges induced by choosing one
 * arbitrary triangle diagonal in every quad.  The latter imprints a directional
 * bias into fills and also makes red/black smoothing impossible.  Keep a CSR
 * fallback for callers that pass a compact SPAN mesh. */
typedef struct QuadAdj {
    int structured;
    int H, W;
    const int32_t *off;
    const int32_t *tgt;
} QuadAdj;

enum {
    QUILT_EDGE_FITTED_FITTED = 0,
    QUILT_EDGE_FILLED_FILLED = 1,
    QUILT_EDGE_FITTED_FILLED = 2,
    QUILT_EDGE_CLASSES = 3
};

typedef struct QuiltEdgeAccum {
    double sum_sq[QUILT_EDGE_CLASSES];
    size_t edges[QUILT_EDGE_CLASSES];
} QuiltEdgeAccum;

static __inline int quilt_edge_class(uint8_t filled_a, uint8_t filled_b) {
    if (filled_a != filled_b) return QUILT_EDGE_FITTED_FILLED;
    return filled_a ? QUILT_EDGE_FILLED_FILLED : QUILT_EDGE_FITTED_FITTED;
}

static void quilt_edge_accum(QuiltEdgeAccum *a, int edge_class, double jump) {
    a->sum_sq[edge_class] += jump * jump;
    a->edges[edge_class]++;
}

static void quilt_edge_finish(const QuiltEdgeAccum *a,
                              QuadStripSnapRegionEdgeStats *out) {
    double all_sum_sq = 0.0;
    size_t all_edges = 0;
    memset(out, 0, sizeof *out);
    for (int k = 0; k < QUILT_EDGE_CLASSES; k++) {
        all_sum_sq += a->sum_sq[k];
        all_edges += a->edges[k];
    }
    out->all_edges = all_edges;
    out->all_rms = all_edges ? sqrt(all_sum_sq / (double)all_edges) : 0.0;
    out->fitted_fitted_edges = a->edges[QUILT_EDGE_FITTED_FITTED];
    out->filled_filled_edges = a->edges[QUILT_EDGE_FILLED_FILLED];
    out->fitted_filled_edges = a->edges[QUILT_EDGE_FITTED_FILLED];
    out->fitted_fitted_rms = out->fitted_fitted_edges
        ? sqrt(a->sum_sq[QUILT_EDGE_FITTED_FITTED] /
               (double)out->fitted_fitted_edges) : 0.0;
    out->filled_filled_rms = out->filled_filled_edges
        ? sqrt(a->sum_sq[QUILT_EDGE_FILLED_FILLED] /
               (double)out->filled_filled_edges) : 0.0;
    out->fitted_filled_rms = out->fitted_filled_edges
        ? sqrt(a->sum_sq[QUILT_EDGE_FITTED_FILLED] /
               (double)out->fitted_filled_edges) : 0.0;
}

static __inline int qadj_degree(const QuadAdj *a, size_t v) {
    if (!a->structured) return (int)(a->off[v + 1] - a->off[v]);
    int r = (int)(v / (size_t)a->W), c = (int)(v % (size_t)a->W);
    return (r > 0) + (r + 1 < a->H) + (c > 0) + (c + 1 < a->W);
}

static __inline size_t qadj_neighbor(const QuadAdj *a, size_t v, int k) {
    if (!a->structured) return (size_t)a->tgt[a->off[v] + k];
    int r = (int)(v / (size_t)a->W), c = (int)(v % (size_t)a->W);
    if (r > 0) { if (k-- == 0) return v - (size_t)a->W; }
    if (r + 1 < a->H) { if (k-- == 0) return v + (size_t)a->W; }
    if (c > 0) { if (k-- == 0) return v - 1; }
    return v + 1;
}

static __inline int qadj_color(const QuadAdj *a, size_t v) {
    if (!a->structured) return 0;
    return ((int)(v / (size_t)a->W) + (int)(v % (size_t)a->W)) & 1;
}

/* Keep an in-place red/black vertex update on the same oriented side of both
 * emitted triangles in every incident structured quad as at the start of this
 * proximal pass.  The averaged-quad check is retained for broad conditioning,
 * but it cannot see an alternating skew whose two edge estimates cancel.
 * Requiring a fraction of each reference triangle's projected area prevents a
 * legal-looking quad average from hiding a collapsed or reversed half-quad. */
static int structured_update_safe(const float *verts,const float *orig,
                                  int H,int W,size_t moved,
                                  const double proposal[3],
                                  double min_area_fraction,
                                  double max_edge_scale) {
    int vr=(int)(moved/(size_t)W),vc=(int)(moved%(size_t)W);
    for(int r=vr-1;r<=vr;r++)for(int c=vc-1;c<=vc;c++){
        if(r<0||c<0||r+1>=H||c+1>=W)continue;
        size_t id[4]={(size_t)r*W+c,(size_t)r*W+c+1,
                      (size_t)(r+1)*W+c,(size_t)(r+1)*W+c+1};
        double p[4][3];
        for(int q=0;q<4;q++)for(int k=0;k<3;k++)
            p[q][k]=id[q]==moved?proposal[k]:(double)verts[id[q]*3+k];
        {
            int tri[2][3]={{0,1,3},{0,3,2}};
            for(int t=0;t<2;t++) {
                int ia=tri[t][0],ib=tri[t][1],ic=tri[t][2];
                double ce0[3],ce1[3],oe0[3],oe1[3],cn[3],on[3];
                for(int k=0;k<3;k++) {
                    ce0[k]=p[ib][k]-p[ia][k];
                    ce1[k]=p[ic][k]-p[ia][k];
                    oe0[k]=(double)orig[id[ib]*3+k]-orig[id[ia]*3+k];
                    oe1[k]=(double)orig[id[ic]*3+k]-orig[id[ia]*3+k];
                }
                cn[0]=ce0[1]*ce1[2]-ce0[2]*ce1[1];
                cn[1]=ce0[2]*ce1[0]-ce0[0]*ce1[2];
                cn[2]=ce0[0]*ce1[1]-ce0[1]*ce1[0];
                on[0]=oe0[1]*oe1[2]-oe0[2]*oe1[1];
                on[1]=oe0[2]*oe1[0]-oe0[0]*oe1[2];
                on[2]=oe0[0]*oe1[1]-oe0[1]*oe1[0];
                double ref2=on[0]*on[0]+on[1]*on[1]+on[2]*on[2];
                if(ref2>1e-18&&cn[0]*on[0]+cn[1]*on[1]+cn[2]*on[2]
                                  <=min_area_fraction*ref2)
                    return 0;
            }
        }
        double eu[3],ev[3],ou[3],ov[3];
        for(int k=0;k<3;k++){
            eu[k]=0.5*((p[1][k]-p[0][k])+(p[3][k]-p[2][k]));
            ev[k]=0.5*((p[2][k]-p[0][k])+(p[3][k]-p[1][k]));
            ou[k]=0.5*(((double)orig[id[1]*3+k]-orig[id[0]*3+k])
                      +((double)orig[id[3]*3+k]-orig[id[2]*3+k]));
            ov[k]=0.5*(((double)orig[id[2]*3+k]-orig[id[0]*3+k])
                      +((double)orig[id[3]*3+k]-orig[id[1]*3+k]));
        }
        double cn[3]={eu[1]*ev[2]-eu[2]*ev[1],
                      eu[2]*ev[0]-eu[0]*ev[2],
                      eu[0]*ev[1]-eu[1]*ev[0]};
        double on[3]={ou[1]*ov[2]-ou[2]*ov[1],
                      ou[2]*ov[0]-ou[0]*ov[2],
                      ou[0]*ov[1]-ou[1]*ov[0]};
        double ref2=on[0]*on[0]+on[1]*on[1]+on[2]*on[2];
        if(ref2>1e-18&&cn[0]*on[0]+cn[1]*on[1]+cn[2]*on[2]
                          <=min_area_fraction*ref2)
            return 0;
        if(max_edge_scale>1.0){
            double eu2=0.0,ev2=0.0,ou2=0.0,ov2=0.0;
            for(int k=0;k<3;k++){
                eu2+=eu[k]*eu[k];ev2+=ev[k]*ev[k];
                ou2+=ou[k]*ou[k];ov2+=ov[k]*ov[k];
            }
            double hi=max_edge_scale*max_edge_scale,lo=1.0/hi;
            if((ou2>1e-18&&(eu2<lo*ou2||eu2>hi*ou2))
             ||(ov2>1e-18&&(ev2<lo*ov2||ev2>hi*ov2)))return 0;
        }
    }
    return 1;
}

/* The quilted CT target is reach-bounded, but a coupled position solve can
 * otherwise accumulate a pathological tangential excursion in a weakly
 * anchored region.  Keep every proposal inside the same physical trust radius
 * as the widest supported profile search. */
static void clamp_position_from_orig(double p[3],const float *orig,size_t v,
                                     double max_displacement) {
    if(max_displacement<=0.0)return;
    double d[3]={p[0]-(double)orig[v*3],p[1]-(double)orig[v*3+1],
                 p[2]-(double)orig[v*3+2]};
    double dn=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
    if(dn<=max_displacement||dn<=1e-12)return;
    double a=max_displacement/dn;
    for(int k=0;k<3;k++)p[k]=(double)orig[v*3+(size_t)k]+a*d[k];
}

static __inline void preserve_axial_position(double p[3],const float *orig,
                                              size_t v,
                                              const QuadStripSnapOpts *o) {
    if(o->preserve_axial)p[0]=(double)orig[v*3];
}

static __inline void clamp_position_to_winding(
        double p[3],const float *orig,size_t v,const float *winding_limit,
        const QuadStripSnapOpts *o) {
    if(winding_limit&&snap_oriented(o)) {
        double oy=(double)orig[v*3+1]-o->axis_y;
        double ox=(double)orig[v*3+2]-o->axis_x;
        double py=p[1]-o->axis_y,px=p[2]-o->axis_x;
        double r0=hypot(oy,ox),r=hypot(py,px),g=(double)winding_limit[v];
        if(r>1e-9&&r0>1e-9&&g>0.0) {
            double bounded=r;
            if(bounded<r0-g)bounded=r0-g;
            if(bounded>r0+g)bounded=r0+g;
            if(bounded<1e-6)bounded=1e-6;
            if(fabs(bounded-r)>1e-12) {
                double scale=bounded/r;
                p[1]=o->axis_y+scale*py;
                p[2]=o->axis_x+scale*px;
            }
        }
    }
    preserve_axial_position(p,orig,v,o);
}

/* --------------------------------------------------------------------------
 * Structured aggregate multigrid
 *
 * The fine operator in both scalar quilting and fixed-rotation ARAP is a
 * screened weighted four-neighbour Laplacian.  A 2x2 piecewise-constant
 * prolongation gives a particularly small Galerkin hierarchy: fine edges
 * inside an aggregate cancel, crossing edges add, and all point/Dirichlet
 * terms add to the coarse diagonal.  Every coarse graph remains a structured
 * four-neighbour bipartite grid, so red/black Gauss-Seidel is valid at every
 * level.  Only coarse levels are stored; the million-vertex fine operator is
 * evaluated from its existing masks/weights in place. */

#define QMG_MAX_LEVELS 16
#define QMG_PRE_SWEEPS 2
#define QMG_POST_SWEEPS 2
#define QMG_COARSE_SWEEPS 32

typedef struct QMGLevel {
    int H, W, channels;
    size_t n, active_count;
    uint8_t *active;
    float *diag, *right, *down;
    float *x, *rhs, *residual; /* channel-major: [channel*n + vertex] */
} QMGLevel;

typedef struct QMGHierarchy {
    int nlevel, channels;       /* stored coarse levels; fine is implicit */
    QMGLevel level[QMG_MAX_LEVELS - 1];
} QMGHierarchy;

typedef struct QMGFields {
    float *channel[4];
    size_t stride;
    int channels;
} QMGFields;

typedef struct QMGRunStats {
    int levels_used, corrections, backtracks, rejects;
    int patches_accepted, patches_rejected;
    double residual_before, residual_after;
} QMGRunStats;

static int qmg_level_alloc(QMGLevel *L, int H, int W, int channels) {
    memset(L, 0, sizeof *L);
    if (H < 1 || W < 1 || channels < 1 || channels > 4 ||
        (size_t)H > SIZE_MAX / (size_t)W) return -1;
    size_t n = (size_t)H * (size_t)W;
    if (n > SIZE_MAX / (size_t)channels || n > SIZE_MAX / sizeof(float) ||
        n * (size_t)channels > SIZE_MAX / sizeof(float)) return -1;
    L->H = H; L->W = W; L->n = n; L->channels = channels;
    L->active = (uint8_t *)calloc(n, 1);
    L->diag = (float *)calloc(n, sizeof *L->diag);
    L->right = (float *)calloc(n, sizeof *L->right);
    L->down = (float *)calloc(n, sizeof *L->down);
    L->x = (float *)calloc(n * (size_t)channels, sizeof *L->x);
    L->rhs = (float *)calloc(n * (size_t)channels, sizeof *L->rhs);
    L->residual = (float *)calloc(n * (size_t)channels, sizeof *L->residual);
    if (!L->active || !L->diag || !L->right || !L->down ||
        !L->x || !L->rhs || !L->residual) {
        free(L->active); free(L->diag); free(L->right); free(L->down);
        free(L->x); free(L->rhs); free(L->residual);
        memset(L, 0, sizeof *L); return -1;
    }
    return 0;
}

static void qmg_level_free(QMGLevel *L) {
    free(L->active); free(L->diag); free(L->right); free(L->down);
    free(L->x); free(L->rhs); free(L->residual);
    memset(L, 0, sizeof *L);
}

static void qmg_hierarchy_free(QMGHierarchy *M) {
    for (int l = 0; l < M->nlevel; l++) qmg_level_free(&M->level[l]);
    memset(M, 0, sizeof *M);
}

static size_t qmg_parent_index(int parent_w, int r, int c) {
    return (size_t)(r >> 1) * (size_t)parent_w + (size_t)(c >> 1);
}

/* Add one active-active conductance to a Galerkin aggregate level. */
static int qmg_level_add_edge(QMGLevel *L, int ar, int ac, int br, int bc,
                              double weight) {
    int par = ar >> 1, pac = ac >> 1, pbr = br >> 1, pbc = bc >> 1;
    size_t a = (size_t)par * (size_t)L->W + (size_t)pac;
    size_t b = (size_t)pbr * (size_t)L->W + (size_t)pbc;
    if (a == b) {
        L->diag[a] = (float)((double)L->diag[a] - 2.0 * weight);
        return 0;
    }
    if (par == pbr && abs(pac - pbc) == 1) {
        size_t left = pac < pbc ? a : b;
        L->right[left] = (float)((double)L->right[left] + weight);
        return 0;
    }
    if (pac == pbc && abs(par - pbr) == 1) {
        size_t top = par < pbr ? a : b;
        L->down[top] = (float)((double)L->down[top] + weight);
        return 0;
    }
    return -1;
}

static int qmg_coarsen_level(const QMGLevel *fine, QMGLevel *coarse) {
    int H = (fine->H + 1) / 2, W = (fine->W + 1) / 2;
    if (qmg_level_alloc(coarse, H, W, fine->channels) != 0) return -1;
    for (int r = 0; r < fine->H; r++) for (int c = 0; c < fine->W; c++) {
        size_t v = (size_t)r * (size_t)fine->W + (size_t)c;
        if (!fine->active[v]) continue;
        size_t p = qmg_parent_index(W, r, c);
        if (!coarse->active[p]) { coarse->active[p] = 1; coarse->active_count++; }
        coarse->diag[p] = (float)((double)coarse->diag[p] + fine->diag[v]);
    }
    for (int r = 0; r < fine->H; r++) for (int c = 0; c < fine->W; c++) {
        size_t v = (size_t)r * (size_t)fine->W + (size_t)c;
        if (!fine->active[v]) continue;
        if (c + 1 < fine->W && fine->right[v] > 0.0f && fine->active[v + 1] &&
            qmg_level_add_edge(coarse, r, c, r, c + 1, fine->right[v]) != 0)
            return -1;
        if (r + 1 < fine->H && fine->down[v] > 0.0f &&
            fine->active[v + (size_t)fine->W] &&
            qmg_level_add_edge(coarse, r, c, r + 1, c, fine->down[v]) != 0)
            return -1;
    }
    return 0;
}

static int qmg_finish_hierarchy(QMGHierarchy *M, int requested_levels) {
    while (1 + M->nlevel < requested_levels && M->nlevel < QMG_MAX_LEVELS - 1) {
        QMGLevel *fine = &M->level[M->nlevel - 1];
        if (fine->active_count <= 32 || (fine->H == 1 && fine->W == 1)) break;
        QMGLevel *coarse = &M->level[M->nlevel];
        if (qmg_coarsen_level(fine, coarse) != 0) {
            qmg_level_free(coarse); return -1;
        }
        M->nlevel++;
    }
    return 0;
}

static void qmg_level_smooth(QMGLevel *L, int sweeps, double omega) {
    for (int sweep = 0; sweep < sweeps; sweep++) for (int color = 0; color < 2; color++)
        for (int r = 0; r < L->H; r++) for (int c = 0; c < L->W; c++) {
            size_t v = (size_t)r * (size_t)L->W + (size_t)c;
            if (!L->active[v] || ((r + c) & 1) != color || L->diag[v] <= 1e-12f) continue;
            for (int ch = 0; ch < L->channels; ch++) {
                size_t cv = (size_t)ch * L->n + v;
                double num = L->rhs[cv];
                if (c > 0 && L->active[v - 1])
                    num += (double)L->right[v - 1] * L->x[cv - 1];
                if (c + 1 < L->W && L->active[v + 1])
                    num += (double)L->right[v] * L->x[cv + 1];
                if (r > 0 && L->active[v - (size_t)L->W])
                    num += (double)L->down[v - (size_t)L->W] * L->x[cv - (size_t)L->W];
                if (r + 1 < L->H && L->active[v + (size_t)L->W])
                    num += (double)L->down[v] * L->x[cv + (size_t)L->W];
                double old = L->x[cv], opt = num / (double)L->diag[v];
                L->x[cv] = (float)(old + omega * (opt - old));
            }
        }
}

static void qmg_level_residual(QMGLevel *L) {
    for (int ch = 0; ch < L->channels; ch++) for (int r = 0; r < L->H; r++)
        for (int c = 0; c < L->W; c++) {
            size_t v = (size_t)r * (size_t)L->W + (size_t)c;
            size_t cv = (size_t)ch * L->n + v;
            if (!L->active[v]) { L->residual[cv] = 0.0f; continue; }
            double ax = (double)L->diag[v] * L->x[cv];
            if (c > 0 && L->active[v - 1])
                ax -= (double)L->right[v - 1] * L->x[cv - 1];
            if (c + 1 < L->W && L->active[v + 1])
                ax -= (double)L->right[v] * L->x[cv + 1];
            if (r > 0 && L->active[v - (size_t)L->W])
                ax -= (double)L->down[v - (size_t)L->W] * L->x[cv - (size_t)L->W];
            if (r + 1 < L->H && L->active[v + (size_t)L->W])
                ax -= (double)L->down[v] * L->x[cv + (size_t)L->W];
            L->residual[cv] = (float)((double)L->rhs[cv] - ax);
        }
}

static void qmg_restrict_residual(const QMGLevel *fine, QMGLevel *coarse) {
    memset(coarse->rhs, 0, coarse->n * (size_t)coarse->channels * sizeof *coarse->rhs);
    memset(coarse->x, 0, coarse->n * (size_t)coarse->channels * sizeof *coarse->x);
    for (int ch = 0; ch < fine->channels; ch++) for (int r = 0; r < fine->H; r++)
        for (int c = 0; c < fine->W; c++) {
            size_t v = (size_t)r * (size_t)fine->W + (size_t)c;
            if (!fine->active[v]) continue;
            size_t p = qmg_parent_index(coarse->W, r, c);
            coarse->rhs[(size_t)ch * coarse->n + p] +=
                fine->residual[(size_t)ch * fine->n + v];
        }
}

static void qmg_prolong_add(QMGLevel *fine, const QMGLevel *coarse) {
    for (int ch = 0; ch < fine->channels; ch++) for (int r = 0; r < fine->H; r++)
        for (int c = 0; c < fine->W; c++) {
            size_t v = (size_t)r * (size_t)fine->W + (size_t)c;
            if (!fine->active[v]) continue;
            size_t p = qmg_parent_index(coarse->W, r, c);
            fine->x[(size_t)ch * fine->n + v] +=
                coarse->x[(size_t)ch * coarse->n + p];
        }
}

static void qmg_vcycle(QMGHierarchy *M, int level_index) {
    QMGLevel *L = &M->level[level_index];
    qmg_level_smooth(L, QMG_PRE_SWEEPS, 1.0);
    if (level_index + 1 == M->nlevel) {
        qmg_level_smooth(L, QMG_COARSE_SWEEPS, 1.0);
        return;
    }
    qmg_level_residual(L);
    qmg_restrict_residual(L, &M->level[level_index + 1]);
    qmg_vcycle(M, level_index + 1);
    qmg_prolong_add(L, &M->level[level_index + 1]);
    qmg_level_smooth(L, QMG_POST_SWEEPS, 1.0);
}

static void qmg_correction_terms(const QMGLevel *L, double *out_rdot,
                                 double *out_eae) {
    double rdot = 0.0, eae = 0.0;
    for (int ch = 0; ch < L->channels; ch++) for (int r = 0; r < L->H; r++)
        for (int c = 0; c < L->W; c++) {
            size_t v = (size_t)r * (size_t)L->W + (size_t)c;
            if (!L->active[v]) continue;
            size_t cv = (size_t)ch * L->n + v;
            double x = L->x[cv];
            rdot += (double)L->rhs[cv] * x;
            eae += (double)L->diag[v] * x * x;
            if (c + 1 < L->W && L->active[v + 1])
                eae -= 2.0 * (double)L->right[v] * x * L->x[cv + 1];
            if (r + 1 < L->H && L->active[v + (size_t)L->W])
                eae -= 2.0 * (double)L->down[v] * x * L->x[cv + (size_t)L->W];
        }
    *out_rdot = rdot; *out_eae = eae;
}

static double qmg_quilt_weight(const uint8_t *filled, size_t a, size_t b,
                               double smooth, double seam_smooth) {
    return smooth * (filled[a] != filled[b] ? seam_smooth : 1.0);
}

/* Robustify the raw target term before solving the quilt.  A wrong ridge
 * family presents as one physical offset vector disagreeing with the vectors
 * around it; signed scalar depths can hide that disagreement when neighboring
 * profile directions rotate or flip.  Huber weights leave coherent candidates
 * at unit strength and let the sheet energy bridge isolated family switches. */
static double quilt_small_median(double value[4],int n) {
    for(int i=1;i<n;i++) {
        double x=value[i];int j=i-1;
        while(j>=0&&value[j]>x){value[j+1]=value[j];j--;}
        value[j+1]=x;
    }
    return n&1?value[n/2]:0.5*(value[n/2-1]+value[n/2]);
}

static void quilt_candidate_weights(const uint8_t *cand, const float *cand_t,
                                     const float *guide,
                                     const uint8_t *anchored,
                                     const QuadAdj *adj, size_t nv,
                                     double huber, float *weight,
                                     int *out_downweighted,
                                     double *out_mean) {
    int downweighted=0;double sum=0.0;size_t count=0;
    int maxdeg=0;
    for(size_t v=0;v<nv;v++) {
        int deg=qadj_degree(adj,v);if(deg>maxdeg)maxdeg=deg;
    }
    double *neighbor=maxdeg>0?(double *)malloc((size_t)3*(size_t)maxdeg*
                                                sizeof *neighbor):NULL;
    for(size_t v=0;v<nv;v++) {
        if(!cand[v]) { weight[v]=0.0f;continue; }
        double w=1.0;
        if(huber>0.0&&anchored[v]!=1&&neighbor) {
            int n=0;
            int deg=qadj_degree(adj,v);
            for(int e=0;e<deg;e++) {
                size_t j=qadj_neighbor(adj,v,e);if(!cand[j])continue;
                for(int k=0;k<3;k++)
                    neighbor[(size_t)k*(size_t)maxdeg+(size_t)n]=
                        (double)cand_t[j]*(double)guide[j*3+(size_t)k];
                n++;
            }
            if(n>=2) {
                double d2=0.0;
                for(int k=0;k<3;k++) {
                    double d=(double)cand_t[v]*(double)guide[v*3+(size_t)k]
                            -quilt_small_median(neighbor+(size_t)k*(size_t)maxdeg,n);
                    d2+=d*d;
                }
                double d=sqrt(d2);
                if(d>huber)w=fmax(0.05,huber/d);
            }
        }
        weight[v]=(float)w;sum+=w;count++;
        if(w<0.999)downweighted++;
    }
    free(neighbor);
    if(out_downweighted)*out_downweighted=downweighted;
    if(out_mean)*out_mean=count?sum/(double)count:0.0;
}

static int qmg_build_quilt(QMGHierarchy *M, const uint8_t *frozen,
                           const uint8_t *filled, const float *cand_w,
                           int H, int W, double smooth, double seam_smooth,
                           double anchor, int requested_levels) {
    memset(M, 0, sizeof *M); M->channels = 4;
    if (requested_levels <= 1) return 0;
    QMGLevel *L = &M->level[0];
    if (qmg_level_alloc(L, (H + 1) / 2, (W + 1) / 2, 4) != 0) {
        qmg_hierarchy_free(M); return -1;
    }
    M->nlevel = 1;
    for (int r = 0; r < H; r++) for (int c = 0; c < W; c++) {
        size_t v = (size_t)r * (size_t)W + (size_t)c;
        if (frozen[v]) continue;
        size_t p = qmg_parent_index(L->W, r, c);
        if (!L->active[p]) { L->active[p] = 1; L->active_count++; }
        double diag = (double)cand_w[v] + anchor;
        if (c > 0) diag += qmg_quilt_weight(filled, v, v - 1, smooth, seam_smooth);
        if (c + 1 < W) diag += qmg_quilt_weight(filled, v, v + 1, smooth, seam_smooth);
        if (r > 0) diag += qmg_quilt_weight(filled, v, v - (size_t)W, smooth, seam_smooth);
        if (r + 1 < H) diag += qmg_quilt_weight(filled, v, v + (size_t)W, smooth, seam_smooth);
        L->diag[p] = (float)((double)L->diag[p] + diag);
    }
    for (int r = 0; r < H; r++) for (int c = 0; c < W; c++) {
        size_t v = (size_t)r * (size_t)W + (size_t)c;
        if (frozen[v]) continue;
        if (c + 1 < W && !frozen[v + 1]) {
            double w = qmg_quilt_weight(filled, v, v + 1, smooth, seam_smooth);
            if (qmg_level_add_edge(L, r, c, r, c + 1, w) != 0) goto fail;
        }
        if (r + 1 < H && !frozen[v + (size_t)W]) {
            double w = qmg_quilt_weight(filled, v, v + (size_t)W, smooth, seam_smooth);
            if (qmg_level_add_edge(L, r, c, r + 1, c, w) != 0) goto fail;
        }
    }
    if (qmg_finish_hierarchy(M, requested_levels) != 0) goto fail;
    return 0;
fail:
    qmg_hierarchy_free(M); return -1;
}

static void qmg_quilt_fine_smooth(float *offset, float *conf,
                                  const uint8_t *frozen, const uint8_t *filled,
                                  const float *cand_w, const float *cand_t,
                                  const float *guide,
                                  int H, int W, double smooth,
                                  double seam_smooth, double anchor, int sweeps) {
    size_t nv = (size_t)H * (size_t)W;
    for (int sweep = 0; sweep < sweeps; sweep++) for (int color = 0; color < 2; color++) {
        ptrdiff_t vv;
        /* Within one color the 4-neighbour stencil reads only the frozen
         * opposite color, so the parallel sweep is bit-exact Gauss-Seidel. */
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (vv = 0; vv < (ptrdiff_t)nv; vv++) {
            size_t v = (size_t)vv;
            int r = (int)(v / (size_t)W), c = (int)(v % (size_t)W);
            if (frozen[v] || ((r + c) & 1) != color) continue;
            double wd = (double)cand_w[v];
            double den = wd + anchor;
            double nd[3] = {wd*(double)cand_t[v]*(double)guide[v*3+0],
                            wd*(double)cand_t[v]*(double)guide[v*3+1],
                            wd*(double)cand_t[v]*(double)guide[v*3+2]};
            double nc = wd;
#define QMG_QUILT_NBR(j_) do { \
                size_t j = (j_); \
                double ws = qmg_quilt_weight(filled, v, j, smooth, seam_smooth); \
                den += ws; \
                if (!frozen[j]) { \
                    for(int k_=0;k_<3;k_++)nd[k_] += ws*(double)offset[j*3+(size_t)k_]; \
                    nc += ws * (double)conf[j]; \
                } \
            } while (0)
            if (c > 0) QMG_QUILT_NBR(v - 1);
            if (c + 1 < W) QMG_QUILT_NBR(v + 1);
            if (r > 0) QMG_QUILT_NBR(v - (size_t)W);
            if (r + 1 < H) QMG_QUILT_NBR(v + (size_t)W);
#undef QMG_QUILT_NBR
            for(int k=0;k<3;k++)
                offset[v*3+(size_t)k] = den > 1e-12 ? (float)(nd[k] / den) : 0.0f;
            conf[v] = den > 1e-12 ? (float)(nc / den) : 0.0f;
        }
    }
}

static double qmg_quilt_restrict_residual(QMGHierarchy *M, const float *offset,
                                          const float *conf,
                                          const uint8_t *frozen,
                                          const uint8_t *filled,
                                          const float *cand_w,
                                          const float *cand_t, const float *guide,
                                          int H, int W, double smooth,
                                          double seam_smooth, double anchor) {
    QMGLevel *L = &M->level[0];
    memset(L->rhs, 0, L->n * 4u * sizeof *L->rhs);
    double r2 = 0.0, b2 = 0.0;
    for (int r = 0; r < H; r++) for (int c = 0; c < W; c++) {
        size_t v = (size_t)r * (size_t)W + (size_t)c;
        if (frozen[v]) continue;
        double wd = (double)cand_w[v], diag = wd + anchor;
        double ax[4] = {0.0, 0.0, 0.0, 0.0};
        double b[4] = {
            wd*(double)cand_t[v]*(double)guide[v*3+0],
            wd*(double)cand_t[v]*(double)guide[v*3+1],
            wd*(double)cand_t[v]*(double)guide[v*3+2], wd};
#define QMG_QUILT_RES_NBR(j_) do { \
                size_t j = (j_); \
                double ws = qmg_quilt_weight(filled, v, j, smooth, seam_smooth); \
                diag += ws; \
                if (!frozen[j]) { \
                    for(int k_=0;k_<3;k_++)ax[k_] -= ws*(double)offset[j*3+(size_t)k_]; \
                    ax[3] -= ws*(double)conf[j]; \
                } \
            } while (0)
        if (c > 0) QMG_QUILT_RES_NBR(v - 1);
        if (c + 1 < W) QMG_QUILT_RES_NBR(v + 1);
        if (r > 0) QMG_QUILT_RES_NBR(v - (size_t)W);
        if (r + 1 < H) QMG_QUILT_RES_NBR(v + (size_t)W);
#undef QMG_QUILT_RES_NBR
        size_t p = qmg_parent_index(L->W, r, c);
        for (int ch = 0; ch < 4; ch++) {
            double value = ch < 3 ? (double)offset[v*3+(size_t)ch] : (double)conf[v];
            ax[ch] += diag * value;
            double res = b[ch] - ax[ch];
            L->rhs[(size_t)ch * L->n + p] += (float)res;
            r2 += res * res; b2 += b[ch] * b[ch];
        }
    }
    return sqrt(r2) / fmax(1e-30, sqrt(b2));
}

static int qmg_solve_quilt(float *offset, float *conf,
                           const uint8_t *frozen, const uint8_t *filled,
                           const float *cand_w, const float *cand_t,
                           const float *guide,
                           int H, int W, const QuadStripSnapOpts *o,
                           QMGRunStats *stats) {
    memset(stats, 0, sizeof *stats); stats->levels_used = 1;
    qmg_quilt_fine_smooth(offset, conf, frozen, filled, cand_w, cand_t, guide, H, W,
                          o->quilt_smooth, o->quilt_seam_smooth,
                          o->quilt_anchor, o->quilt_sweeps);
    if (o->multigrid_levels <= 1 || o->multigrid_cycles <= 0) return 0;
    QMGHierarchy M;
    if (qmg_build_quilt(&M, frozen, filled, cand_w, H, W, o->quilt_smooth,
                        o->quilt_seam_smooth, o->quilt_anchor,
                        o->multigrid_levels) != 0) return -1;
    stats->levels_used = 1 + M.nlevel;
    if (!M.nlevel || !M.level[0].active_count) { qmg_hierarchy_free(&M); return 0; }
    for (int cycle = 0; cycle < o->multigrid_cycles; cycle++) {
        double rel = qmg_quilt_restrict_residual(&M, offset, conf, frozen, filled,
                                                  cand_w, cand_t, guide, H, W,
                                                  o->quilt_smooth,
                                                  o->quilt_seam_smooth,
                                                  o->quilt_anchor);
        if (cycle == 0) stats->residual_before = rel;
        memset(M.level[0].x, 0,
               M.level[0].n * 4u * sizeof *M.level[0].x);
        qmg_vcycle(&M, 0);
        double rdot = 0.0, eae = 0.0;
        qmg_correction_terms(&M.level[0], &rdot, &eae);
        if (rdot > 1e-20 && eae > 1e-20) {
            double alpha = fmin(1.0, rdot / eae);
            QMGLevel *L = &M.level[0];
            for (int r = 0; r < H; r++) for (int c = 0; c < W; c++) {
                size_t v = (size_t)r * (size_t)W + (size_t)c;
                if (frozen[v]) continue;
                size_t p = qmg_parent_index(L->W, r, c);
                double d[3];
                for(int k=0;k<3;k++)
                    d[k]=(double)offset[v*3+(size_t)k]+alpha*L->x[(size_t)k*L->n+p];
                double q = (double)conf[v] + alpha * L->x[3u*L->n + p];
                double dn=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
                double max_reach=o->guided_reach>o->reach?o->guided_reach:o->reach;
                if(dn>max_reach&&dn>1e-12)
                    for(int k=0;k<3;k++)d[k]*=max_reach/dn;
                if (q < 0.0) q = 0.0; if (q > 1.0) q = 1.0;
                for(int k=0;k<3;k++)offset[v*3+(size_t)k]=(float)d[k];
                conf[v] = (float)q;
            }
            stats->corrections++;
        } else stats->rejects++;
        qmg_quilt_fine_smooth(offset, conf, frozen, filled, cand_w, cand_t, guide, H, W,
                              o->quilt_smooth, o->quilt_seam_smooth,
                              o->quilt_anchor, 1);
    }
    stats->residual_after = qmg_quilt_restrict_residual(
        &M, offset, conf, frozen, filled, cand_w, cand_t, guide, H, W,
        o->quilt_smooth, o->quilt_seam_smooth, o->quilt_anchor);
    qmg_hierarchy_free(&M);
    return 0;
}

/* anchored state: 0 ordinary, 1 recovered fill zero-depth anchor,
 * 2 upstream fitted support with a high soft source attachment. */
static double snap_orig_anchor_multiplier(const QuadStripSnapOpts *o,
                                          uint8_t anchored) {
    if(anchored==2)return o->fitted_source_anchor;
    if(anchored==1)return o->filled_target_anchor;
    return 1.0;
}

static double snap_orig_anchor_coordinate(const QuadStripSnapOpts *o,
                                          const float *pass_orig,size_t v,
                                          int channel,uint8_t anchored) {
    const float *position=anchored==2&&o->fitted_source_positions
        ?o->fitted_source_positions:pass_orig;
    return (double)position[v*3+(size_t)channel];
}

static int qmg_build_position(QMGHierarchy *M, const uint8_t *frozen,
                              const uint8_t *anchored, const uint8_t *hast,
                              const float *conf, int H, int W,
                              const QuadStripSnapOpts *o) {
    memset(M, 0, sizeof *M); M->channels = 3;
    if (o->multigrid_levels <= 1) return 0;
    double wn = o->w_smooth + o->w_dev;
    if (wn < 0.0) return -1;
    QMGLevel *L = &M->level[0];
    if (qmg_level_alloc(L, (H + 1) / 2, (W + 1) / 2, 3) != 0) {
        qmg_hierarchy_free(M); return -1;
    }
    M->nlevel = 1;
    for (int r = 0; r < H; r++) for (int c = 0; c < W; c++) {
        size_t v = (size_t)r * (size_t)W + (size_t)c;
        if (frozen[v]) continue;
        size_t p = qmg_parent_index(L->W, r, c);
        if (!L->active[p]) { L->active[p] = 1; L->active_count++; }
        int deg = (r > 0) + (r + 1 < H) + (c > 0) + (c + 1 < W);
        double wo = o->w_orig * snap_orig_anchor_multiplier(o,anchored[v]);
        double wt = hast[v] ? o->w_target * (double)conf[v] : 0.0;
        L->diag[p] = (float)((double)L->diag[p] + wo + wt + (double)deg * wn);
    }
    if (wn > 0.0) for (int r = 0; r < H; r++) for (int c = 0; c < W; c++) {
        size_t v = (size_t)r * (size_t)W + (size_t)c;
        if (frozen[v]) continue;
        if (c + 1 < W && !frozen[v + 1] &&
            qmg_level_add_edge(L, r, c, r, c + 1, wn) != 0) goto fail;
        if (r + 1 < H && !frozen[v + (size_t)W] &&
            qmg_level_add_edge(L, r, c, r + 1, c, wn) != 0) goto fail;
    }
    if (qmg_finish_hierarchy(M, o->multigrid_levels) != 0) goto fail;
    return 0;
fail:
    qmg_hierarchy_free(M); return -1;
}

/* Right-hand side and diagonal of the fixed-rotation screened ARAP system.
 * Movable neighbours remain off-diagonal unknowns; frozen-neighbour values are
 * transferred to the right-hand side. */
static double qmg_position_rhs_diag(const QuadAdj *adj, const float *verts,
                                    const float *orig, const float *target,
                                    const float *uv, const double *R,
                                    const uint8_t *frozen,
                                    const uint8_t *anchored,
                                    const uint8_t *hast, const float *conf,
                                    size_t v, const QuadStripSnapOpts *o,
                                    double rhs[3]) {
    double wo = o->w_orig * snap_orig_anchor_multiplier(o,anchored[v]);
    double wt = hast[v] ? o->w_target * (double)conf[v] : 0.0;
    double wn = o->w_smooth + o->w_dev;
    rhs[0] = wo * snap_orig_anchor_coordinate(o,orig,v,0,anchored[v]);
    rhs[1] = wo * snap_orig_anchor_coordinate(o,orig,v,1,anchored[v]);
    rhs[2] = wo * snap_orig_anchor_coordinate(o,orig,v,2,anchored[v]);
    if (hast[v]) for (int k = 0; k < 3; k++)
        rhs[k] += wt * (double)target[v*3+k];
    double ui = (double)uv[v*2+0], vi = (double)uv[v*2+1];
    const double *Rv = &R[v*9];
    int deg = qadj_degree(adj, v);
    for (int e = 0; e < deg; e++) {
        size_t j = qadj_neighbor(adj, v, e);
        if (frozen[j]) for (int k = 0; k < 3; k++)
            rhs[k] += wn * (double)verts[j*3+k];
        if (o->w_dev != 0.0) {
            double er[2] = {ui - (double)uv[j*2+0], vi - (double)uv[j*2+1]};
            const double *Rj = &R[j*9];
            for (int a = 0; a < 3; a++) {
                double re = 0.5 * (Rv[a*3+0] + Rj[a*3+0]) * er[0]
                          + 0.5 * (Rv[a*3+1] + Rj[a*3+1]) * er[1];
                rhs[a] += o->w_dev * re;
            }
        }
    }
    return wo + wt + (double)deg * wn;
}

static void qmg_position_fine_smooth(float *verts, const float *orig,
                                     const float *target, const float *uv,
                                     const double *R, const uint8_t *frozen,
                                     const uint8_t *anchored,
                                     const uint8_t *hast, const float *conf,
                                     const float *winding_limit,
                                     const QuadAdj *adj,
                                     const QuadStripSnapOpts *o, int sweeps,
                                     int *barrier_backtracks,
                                     int *barrier_rejects) {
    size_t nv = (size_t)adj->H * (size_t)adj->W;
    double wn = o->w_smooth + o->w_dev;
    for (int sweep = 0; sweep < sweeps; sweep++) for (int color = 0; color < 2; color++)
        for (size_t v = 0; v < nv; v++) {
            if (frozen[v] || qadj_color(adj, v) != color) continue;
            double num[3], den = qmg_position_rhs_diag(
                adj, verts, orig, target, uv, R, frozen, anchored, hast, conf,
                v, o, num);
            int deg = qadj_degree(adj, v);
            for (int e = 0; e < deg; e++) {
                size_t j = qadj_neighbor(adj, v, e);
                if (!frozen[j]) for (int k = 0; k < 3; k++)
                    num[k] += wn * (double)verts[j*3+k];
            }
            if (den < 1e-9) den = 1.0;
            double cur[3] = {verts[v*3+0], verts[v*3+1], verts[v*3+2]};
            double opt[3] = {num[0]/den, num[1]/den, num[2]/den};
            double alpha = o->solve_omega, prop[3]; int bt = 0;
            for (;;) {
                for (int k = 0; k < 3; k++) prop[k] = cur[k] + alpha * (opt[k] - cur[k]);
                clamp_position_from_orig(prop,orig,v,o->guided_reach);
                clamp_position_to_winding(prop,orig,v,winding_limit,o);
                if (structured_update_safe(verts, orig, adj->H, adj->W, v,
                                           prop, 0.05, 0.0)) break;
                if (++bt >= 10) {
                    alpha = 0.0;
                    for (int k = 0; k < 3; k++) prop[k] = cur[k];
                    break;
                }
                alpha *= 0.5;
            }
            *barrier_backtracks += bt;
            if (alpha == 0.0) (*barrier_rejects)++;
            for (int k = 0; k < 3; k++) verts[v*3+k] = (float)prop[k];
        }
}

static double qmg_position_restrict_residual(
        QMGHierarchy *M, const float *verts, const float *orig,
        const float *target, const float *uv, const double *R,
        const uint8_t *frozen, const uint8_t *anchored,
        const uint8_t *hast, const float *conf, const QuadAdj *adj,
        const QuadStripSnapOpts *o) {
    QMGLevel *L = &M->level[0];
    memset(L->rhs, 0, L->n * 3u * sizeof *L->rhs);
    double r2 = 0.0, b2 = 0.0, wn = o->w_smooth + o->w_dev;
    size_t nv = (size_t)adj->H * (size_t)adj->W;
    for (size_t v = 0; v < nv; v++) {
        if (frozen[v]) continue;
        double b[3], diag = qmg_position_rhs_diag(
            adj, verts, orig, target, uv, R, frozen, anchored, hast, conf,
            v, o, b);
        double ax[3] = {diag * (double)verts[v*3+0],
                        diag * (double)verts[v*3+1],
                        diag * (double)verts[v*3+2]};
        int deg = qadj_degree(adj, v);
        for (int e = 0; e < deg; e++) {
            size_t j = qadj_neighbor(adj, v, e);
            if (!frozen[j]) for (int k = 0; k < 3; k++)
                ax[k] -= wn * (double)verts[j*3+k];
        }
        int r = (int)(v / (size_t)adj->W), c = (int)(v % (size_t)adj->W);
        size_t p = qmg_parent_index(L->W, r, c);
        for (int k = 0; k < 3; k++) {
            double res = b[k] - ax[k];
            L->rhs[(size_t)k * L->n + p] += (float)res;
            r2 += res * res; b2 += b[k] * b[k];
        }
    }
    return sqrt(r2) / fmax(1e-30, sqrt(b2));
}

static void qmg_corrected_position(const float *verts, const uint8_t *frozen,
                                   const QMGLevel *L, int fine_w, size_t v,
                                   double alpha, double p[3]) {
    int r = (int)(v / (size_t)fine_w), c = (int)(v % (size_t)fine_w);
    size_t parent = qmg_parent_index(L->W, r, c);
    for (int k = 0; k < 3; k++)
        p[k] = (double)verts[v*3+k] + (frozen[v] ? 0.0
              : alpha * (double)L->x[(size_t)k * L->n + parent]);
}

static int qmg_parent_in_patch(const QMGLevel *L, int fine_w, size_t v,
                               int pr0, int pr1, int pc0, int pc1) {
    int r = (int)(v / (size_t)fine_w), c = (int)(v % (size_t)fine_w);
    int pr = r >> 1, pc = c >> 1;
    return pr >= pr0 && pr < pr1 && pc >= pc0 && pc < pc1 &&
           pr < L->H && pc < L->W;
}

static void qmg_corrected_position_patch(
        const float *verts, const uint8_t *frozen, const QMGLevel *L,
        const float *orig,const float *winding_limit,const QuadStripSnapOpts *o,
        int fine_w,size_t v,double alpha,
        int pr0, int pr1, int pc0, int pc1, double p[3]) {
    int use = !frozen[v] && qmg_parent_in_patch(L, fine_w, v,
                                                pr0, pr1, pc0, pc1);
    int r = (int)(v / (size_t)fine_w), c = (int)(v % (size_t)fine_w);
    size_t parent = qmg_parent_index(L->W, r, c);
    for (int k = 0; k < 3; k++)
        p[k] = (double)verts[v*3+k] + (use ?
            alpha * (double)L->x[(size_t)k * L->n + parent] : 0.0);
    if(use)clamp_position_from_orig(p,orig,v,o->guided_reach);
    clamp_position_to_winding(p,orig,v,winding_limit,o);
}

/* Check only the fine quads touched by one coarse patch.  Previously a single
 * unsafe triangle anywhere on a million-vertex ribbon vetoed the entire low-
 * frequency correction.  Sequential patch commits make the final state safe:
 * each later patch is checked against all earlier committed neighbours. */
static int qmg_position_patch_safe(const float *verts, const float *orig,
                                   const uint8_t *frozen,
                                   const QMGLevel *L, int H, int W,
                                   double alpha,const float *winding_limit,
                                   const QuadStripSnapOpts *o,
                                   int pr0, int pr1, int pc0, int pc1) {
    int r0 = 2*pr0-1, r1 = 2*pr1;
    int c0 = 2*pc0-1, c1 = 2*pc1;
    if (r0 < 0) r0 = 0; if (c0 < 0) c0 = 0;
    if (r1 > H-1) r1 = H-1; if (c1 > W-1) c1 = W-1;
    for (int r = r0; r < r1; r++) for (int c = c0; c < c1; c++) {
        size_t id[4] = {(size_t)r*W+c, (size_t)r*W+c+1,
                        (size_t)(r+1)*W+c, (size_t)(r+1)*W+c+1};
        double p[4][3];
        for (int q = 0; q < 4; q++) {
            qmg_corrected_position_patch(verts,frozen,L,orig,winding_limit,o,
                                         W,id[q],alpha,
                                         pr0, pr1, pc0, pc1, p[q]);
            for (int k = 0; k < 3; k++) if (!isfinite(p[q][k])) return 0;
        }
        int tri[2][3] = {{0,1,3},{0,3,2}};
        for (int t = 0; t < 2; t++) {
            int ia=tri[t][0], ib=tri[t][1], ic=tri[t][2];
            double ce0[3],ce1[3],oe0[3],oe1[3],cn[3],on[3];
            for (int k = 0; k < 3; k++) {
                ce0[k]=p[ib][k]-p[ia][k]; ce1[k]=p[ic][k]-p[ia][k];
                oe0[k]=(double)orig[id[ib]*3+k]-orig[id[ia]*3+k];
                oe1[k]=(double)orig[id[ic]*3+k]-orig[id[ia]*3+k];
            }
            cn[0]=ce0[1]*ce1[2]-ce0[2]*ce1[1];
            cn[1]=ce0[2]*ce1[0]-ce0[0]*ce1[2];
            cn[2]=ce0[0]*ce1[1]-ce0[1]*ce1[0];
            on[0]=oe0[1]*oe1[2]-oe0[2]*oe1[1];
            on[1]=oe0[2]*oe1[0]-oe0[0]*oe1[2];
            on[2]=oe0[0]*oe1[1]-oe0[1]*oe1[0];
            double ref2=on[0]*on[0]+on[1]*on[1]+on[2]*on[2];
            if (ref2>1e-18 && cn[0]*on[0]+cn[1]*on[1]+cn[2]*on[2] <= 0.05*ref2)
                return 0;
        }
        double eu[3],ev[3],ou[3],ov[3];
        for (int k=0;k<3;k++) {
            eu[k]=0.5*((p[1][k]-p[0][k])+(p[3][k]-p[2][k]));
            ev[k]=0.5*((p[2][k]-p[0][k])+(p[3][k]-p[1][k]));
            ou[k]=0.5*(((double)orig[id[1]*3+k]-orig[id[0]*3+k])
                      +((double)orig[id[3]*3+k]-orig[id[2]*3+k]));
            ov[k]=0.5*(((double)orig[id[2]*3+k]-orig[id[0]*3+k])
                      +((double)orig[id[3]*3+k]-orig[id[1]*3+k]));
        }
        double cn[3]={eu[1]*ev[2]-eu[2]*ev[1],eu[2]*ev[0]-eu[0]*ev[2],eu[0]*ev[1]-eu[1]*ev[0]};
        double on[3]={ou[1]*ov[2]-ou[2]*ov[1],ou[2]*ov[0]-ou[0]*ov[2],ou[0]*ov[1]-ou[1]*ov[0]};
        double ref2=on[0]*on[0]+on[1]*on[1]+on[2]*on[2];
        if (ref2>1e-18 && cn[0]*on[0]+cn[1]*on[1]+cn[2]*on[2] <= 0.05*ref2)
            return 0;
    }
    return 1;
}

static int qmg_apply_position_correction(float *verts, const float *orig,
                                          const uint8_t *frozen,
                                          const QuadAdj *adj, QMGLevel *L,
                                          const float *winding_limit,
                                          const QuadStripSnapOpts *o,
                                          QMGRunStats *stats) {
    double rdot = 0.0, eae = 0.0;
    qmg_correction_terms(L, &rdot, &eae);
    if (!(rdot > 1e-20) || !(eae > 1e-20)) { stats->rejects++; return 0; }
    double base_alpha = fmin(1.0, rdot / eae);
    int patch = o->multigrid_patch;
    if (patch <= 0) patch = L->H > L->W ? L->H : L->W;
    int accepted = 0;
    for (int pr0 = 0; pr0 < L->H; pr0 += patch) {
        int pr1 = pr0 + patch; if (pr1 > L->H) pr1 = L->H;
        for (int pc0 = 0; pc0 < L->W; pc0 += patch) {
            int pc1 = pc0 + patch; if (pc1 > L->W) pc1 = L->W;
            int useful = 0;
            for (int pr=pr0; pr<pr1 && !useful; pr++) for (int pc=pc0; pc<pc1; pc++) {
                size_t p=(size_t)pr*(size_t)L->W+(size_t)pc;
                if(!L->active[p])continue;
                for(int k=0;k<3;k++) if(fabs((double)L->x[(size_t)k*L->n+p])>1e-12){
                    useful=1;break;
                }
                if(useful)break;
            }
            if(!useful)continue;
            double alpha = base_alpha;
            int bt = 0;
            while (bt < 12 && !qmg_position_patch_safe(
                       verts,orig,frozen,L,adj->H,adj->W,alpha,winding_limit,o,
                       pr0, pr1, pc0, pc1)) {
                alpha *= 0.5; bt++;
            }
            stats->backtracks += bt;
            if (bt >= 12 || alpha < 1.0/4096.0) {
                stats->patches_rejected++;
                continue;
            }
            int fr0 = 2*pr0, fr1 = 2*pr1;
            int fc0 = 2*pc0, fc1 = 2*pc1;
            if (fr1 > adj->H) fr1 = adj->H;
            if (fc1 > adj->W) fc1 = adj->W;
            for (int r = fr0; r < fr1; r++) for (int c = fc0; c < fc1; c++) {
                size_t v = (size_t)r*(size_t)adj->W+(size_t)c;
                if (frozen[v]) continue;
                size_t p = qmg_parent_index(L->W, r, c);
                double proposal[3];
                for(int k=0;k<3;k++)proposal[k]=(double)verts[v*3+k]+
                    alpha*(double)L->x[(size_t)k*L->n+p];
                clamp_position_from_orig(proposal,orig,v,o->guided_reach);
                clamp_position_to_winding(proposal,orig,v,winding_limit,o);
                for(int k=0;k<3;k++)verts[v*3+k]=(float)proposal[k];
            }
            stats->patches_accepted++;
            accepted++;
        }
    }
    if (!accepted) { stats->rejects++; return 0; }
    stats->corrections++;
    return 1;
}

static int qmg_solve_position(float *verts, const float *orig,
                              const float *target, const float *uv,
                              const double *R, const uint8_t *frozen,
                              const uint8_t *anchored, const uint8_t *hast,
                              const float *conf,const float *winding_limit,
                              const QuadAdj *adj,
                              const QuadStripSnapOpts *o,
                              int *barrier_backtracks, int *barrier_rejects,
                              QMGRunStats *stats) {
    memset(stats,0,sizeof *stats); stats->levels_used=1;
    qmg_position_fine_smooth(verts,orig,target,uv,R,frozen,anchored,hast,conf,winding_limit,
                             adj,o,o->solve_sweeps,barrier_backtracks,
                             barrier_rejects);
    if(o->multigrid_levels<=1||o->multigrid_cycles<=0)return 0;
    QMGHierarchy M;
    if(qmg_build_position(&M,frozen,anchored,hast,conf,adj->H,adj->W,o)!=0)return -1;
    stats->levels_used=1+M.nlevel;
    if(!M.nlevel||!M.level[0].active_count){qmg_hierarchy_free(&M);return 0;}
    for(int cycle=0;cycle<o->multigrid_cycles;cycle++){
        double rel=qmg_position_restrict_residual(&M,verts,orig,target,uv,R,
            frozen,anchored,hast,conf,adj,o);
        if(cycle==0)stats->residual_before=rel;
        memset(M.level[0].x,0,M.level[0].n*3u*sizeof *M.level[0].x);
        qmg_vcycle(&M,0);
        qmg_apply_position_correction(verts,orig,frozen,adj,&M.level[0],winding_limit,o,stats);
        qmg_position_fine_smooth(verts,orig,target,uv,R,frozen,anchored,hast,conf,winding_limit,
                                 adj,o,1,barrier_backtracks,barrier_rejects);
    }
    stats->residual_after=qmg_position_restrict_residual(&M,verts,orig,target,uv,R,
        frozen,anchored,hast,conf,adj,o);
    qmg_hierarchy_free(&M);return 0;
}

/* Move one vertex to the nearest bright ridge along its normal.  Samples the
 * intensity profile over [-reach, +reach] and picks the local maximum nearest
 * the current position that clears bright_min and gains >= min_gain over the
 * vertex's own sample.  Returns 1 and sets *out_t on a hit, 0 if none. */
/* Midline band finder tuning (profile-internal, not user knobs). */
#define QSS_MIDLINE_PLATEAU_TOL 0.5    /* gray levels: plateau-midpoint run */
#define QSS_MIDLINE_MIN_HALF_STEPS 1   /* min band half-width, profile steps */

typedef struct RidgeProfile {
    double current;
    double visible_current; /* bake-equivalent max within +/-2 vox */
    double target_t;
    double target_value;
    double visible_target;  /* bake-equivalent max within +/-2 of target */
    double nearest_peak_t;
    double nearest_peak_value;
    double transition_drop; /* recto: edge drop; midline: min(edge drops) */
    double band_half_width; /* midline only: half the bright-band thickness */
    int has_target;
    int on_ridge;
    int used_recto;
    int used_midline;       /* both dark boundaries found; target = center */
    int used_band_fallback; /* merged/clipped band; smoothed bounded maximum */
    int used_local_contrast;
    int used_guided;
    int used_wide;
    int used_outward;
} RidgeProfile;

/* Inspect one signed normal profile.  `on_ridge` is deliberately relative: it
 * asks whether ANY discrete local maximum is already within ridge_lock of the
 * current point, regardless of its absolute CT value.  That is the trusted-data
 * gate.  `has_target` is stricter: the peak must also clear the scroll-adaptive
 * intensity floor and improve on the current sample by min_gain. */
static RidgeProfile ridge_profile(CubeTable *ct, const double p[3], const double n[3],
                                  double reach, double step, double bright_min,
                                  double min_gain, double ridge_lock,
                                  int recto_sign, int local_contrast,
                                  int bidirectional_recto,
                                  int side_mode, double corridor) {
    RidgeProfile out;
    memset(&out, 0, sizeof out);
    out.current = -1.0;
    out.visible_current = -1.0;
    out.target_t = 0.0;
    out.target_value = -1.0;
    out.visible_target = -1.0;
    out.nearest_peak_t = 1e30;
    out.nearest_peak_value = -1.0;
    if (step <= 0.0) step = 0.5;
    if (reach <= 0.0) reach = 6.0;
    int nsteps = (int)(2.0 * reach / step) + 1;
    if (nsteps < 3) nsteps = 3;
    if (nsteps > 512) nsteps = 512;
    double prof[512];
    double profile_min = 1e30, profile_max = -1.0;
    for (int k = 0; k < nsteps; k++) {
        double t = -reach + (2.0 * reach * (double)k) / (double)(nsteps - 1);
        prof[k] = sample_trilinear(ct, p[0] + t * n[0], p[1] + t * n[1], p[2] + t * n[2]);
        if (prof[k] >= 0.0) {
            if (prof[k] < profile_min) profile_min = prof[k];
            if (prof[k] > profile_max) profile_max = prof[k];
        }
        if (fabs(t) <= 2.0 + 1e-9 && prof[k] > out.visible_current)
            out.visible_current = prof[k];
    }
    double I0 = sample_trilinear(ct, p[0], p[1], p[2]);
    out.current = I0;
    double effective_bright = bright_min, effective_gain = min_gain;
    if (local_contrast && profile_max >= profile_min &&
        profile_max - profile_min >= 12.0) {
        double span = profile_max - profile_min;
        double local_bright = profile_min + 0.55 * span;
        double local_gain = fmax(8.0, 0.20 * span);
        if (local_bright < effective_bright) effective_bright = local_bright;
        if (local_gain < effective_gain) effective_gain = local_gain;
    }
    int best_k = -1;
    double best_abs = 1e30;
    for (int k = 1; k < nsteps - 1; k++) {
        if (prof[k] < 0.0) continue;
        if (prof[k] < prof[k - 1] || prof[k] < prof[k + 1]) continue;  /* local max */
        double t = -reach + (2.0 * reach * (double)k) / (double)(nsteps - 1);
        double a = fabs(t);
        if (a < fabs(out.nearest_peak_t)) {
            out.nearest_peak_t = t;
            out.nearest_peak_value = prof[k];
        }
        if (a <= ridge_lock) out.on_ridge = 1;
        if (prof[k] < effective_bright) continue;
        if (I0 >= 0.0 && prof[k] - I0 < effective_gain) continue;
        if (a < best_abs) { best_abs = a; best_k = k; }
    }
    if (best_k >= 0) {
        out.target_t = -reach + (2.0 * reach * (double)best_k) / (double)(nsteps - 1);
        out.target_value = prof[best_k];
        out.has_target = 1;
        out.used_local_contrast = prof[best_k] < bright_min ||
            (I0 >= 0.0 && prof[best_k] - I0 < min_gain);
    }

    /* A local intensity maximum locates the papyrus core but does not identify
     * which FACE of the sheet it belongs to.  When the caller can orient the
     * profile, increasing signed distance is inward (toward the scroll axis).
     * The recto face is then the bright-to-dark transition into the inward gap.
     * Detect it over a two-sample (roughly one-voxel at the default step) chord;
     * this is the native-CT analogue of audit_ink9um_depth_contract.py's
     * strongest negative derivative.  Prefer a transition inward of the
     * current estimate, since a surface initialized near the papyrus core must
     * cross its own recto face in that direction.  A slightly overshot estimate
     * may still lock to a transition within ridge_lock behind it. */
    if (recto_sign != 0 && side_mode == QUAD_SNAP_SIDE_MIDLINE) {
        /* Band-center ("midline") targeting: locate the contiguous bright
         * papyrus band the vertex belongs to (or the nearest credible one
         * within the winding corridor) and aim for the midpoint between its
         * two dark boundaries -- the sheet's middle, NOT the recto face.
         * Both edges reuse the recto branch's two-sample-chord drop
         * detector, one run inward and one outward from a corridor-bounded
         * bright seed.  A band merged with its neighbour (no dark separator
         * within the corridor) or clipped by the corridor falls back to the
         * bounded 3-tap-smoothed maximum, and an all-dark corridor abstains
         * so the vertex relaxes under the sheet terms alone. */
        double tcap = corridor > 0.0 ? corridor : reach;
        if (tcap > reach) tcap = reach;
        int seed_k = -1;
        double seed_abs = 1e30;
        for (int k = 0; k < nsteps; k++) {
            if (prof[k] < effective_bright) continue;
            double t = -reach + (2.0 * reach * (double)k) / (double)(nsteps - 1);
            double a = fabs(t);
            if (a > tcap + 1e-9) continue;
            if (a < seed_abs - 1e-12 ||
                (fabs(a - seed_abs) <= 1e-12 &&
                 (double)recto_sign * t > 0.0)) {
                seed_abs = a;
                seed_k = k;
            }
        }
        if (seed_k >= 0) {
            int in_k = -1, out_k = -1;
            for (int k = seed_k; k - recto_sign >= 0 &&
                                 k + recto_sign >= 0 &&
                                 k - recto_sign < nsteps &&
                                 k + recto_sign < nsteps; k += recto_sign) {
                int ko = k - recto_sign, ki = k + recto_sign;
                if (prof[ko] < 0.0 || prof[ki] < 0.0) break;
                double t = -reach + (2.0 * reach * (double)k) / (double)(nsteps - 1);
                if (fabs(t) > tcap + 1e-9) break;
                if (prof[ko] >= effective_bright &&
                    prof[ko] - prof[ki] >= effective_gain) { in_k = k; break; }
            }
            for (int k = seed_k; k - recto_sign >= 0 &&
                                 k + recto_sign >= 0 &&
                                 k - recto_sign < nsteps &&
                                 k + recto_sign < nsteps; k -= recto_sign) {
                int ko = k + recto_sign, ki = k - recto_sign;
                if (prof[ko] < 0.0 || prof[ki] < 0.0) break;
                double t = -reach + (2.0 * reach * (double)k) / (double)(nsteps - 1);
                if (fabs(t) > tcap + 1e-9) break;
                if (prof[ko] >= effective_bright &&
                    prof[ko] - prof[ki] >= effective_gain) { out_k = k; break; }
            }
            double t_step = 2.0 * reach / (double)(nsteps - 1);
            if (in_k >= 0 && out_k >= 0 &&
                abs(in_k - out_k) >= 2 * QSS_MIDLINE_MIN_HALF_STEPS) {
                double t_in = -reach + t_step * (double)in_k;
                double t_out = -reach + t_step * (double)out_k;
                int center_k = (in_k + out_k) / 2;
                double drop_in = prof[in_k - recto_sign] - prof[in_k + recto_sign];
                double drop_out = prof[out_k + recto_sign] - prof[out_k - recto_sign];
                out.target_t = 0.5 * (t_in + t_out);
                out.target_value = prof[center_k];
                out.transition_drop = drop_in < drop_out ? drop_in : drop_out;
                out.band_half_width = 0.5 * fabs(t_in - t_out);
                out.has_target = 1;
                out.used_midline = 1;
                out.on_ridge = fabs(out.target_t) <= ridge_lock;
                out.used_local_contrast =
                    prof[center_k] < bright_min || out.transition_drop < min_gain;
            } else {
                /* Merged or clipped band: bounded smoothed maximum with a
                 * deterministic plateau midpoint. */
                int best_sm_k = -1;
                double best_sm = -1.0, best_sm_abs = 1e30;
                for (int k = 0; k < nsteps; k++) {
                    if (prof[k] < 0.0) continue;
                    double t = -reach + t_step * (double)k;
                    if (fabs(t) > tcap + 1e-9) continue;
                    double sm = prof[k], wt = 1.0;
                    if (k > 0 && prof[k-1] >= 0.0) { sm += prof[k-1]; wt += 1.0; }
                    if (k+1 < nsteps && prof[k+1] >= 0.0) { sm += prof[k+1]; wt += 1.0; }
                    sm /= wt;
                    if (sm < effective_bright) continue;
                    if (sm > best_sm + 1e-12 ||
                        (fabs(sm - best_sm) <= 1e-12 && fabs(t) < best_sm_abs)) {
                        best_sm = sm; best_sm_k = k; best_sm_abs = fabs(t);
                    }
                }
                if (best_sm_k >= 0) {
                    int lo_k = best_sm_k, hi_k = best_sm_k;
                    while (lo_k > 0 && prof[lo_k-1] >= 0.0 &&
                           prof[lo_k-1] >= prof[best_sm_k] - QSS_MIDLINE_PLATEAU_TOL &&
                           fabs(-reach + t_step*(double)(lo_k-1)) <= tcap + 1e-9)
                        lo_k--;
                    while (hi_k+1 < nsteps && prof[hi_k+1] >= 0.0 &&
                           prof[hi_k+1] >= prof[best_sm_k] - QSS_MIDLINE_PLATEAU_TOL &&
                           fabs(-reach + t_step*(double)(hi_k+1)) <= tcap + 1e-9)
                        hi_k++;
                    int center_k = (lo_k + hi_k) / 2;
                    out.target_t = -reach + t_step * (double)center_k;
                    out.target_value = prof[center_k];
                    out.transition_drop = 0.0;
                    out.band_half_width = 0.5 * t_step * (double)(hi_k - lo_k);
                    out.has_target = 1;
                    out.used_band_fallback = 1;
                    out.on_ridge = fabs(out.target_t) <= ridge_lock;
                    out.used_local_contrast = out.target_value < bright_min;
                }
                /* else: all-dark corridor -- abstain; generic bookkeeping and
                 * any unbounded nearest-max target from above remain, and the
                 * caller's winding veto polices that fallback. */
            }
        } else if (out.has_target && fabs(out.target_t) > tcap + 1e-9) {
            /* No bright sample inside the corridor: abstain rather than keep
             * a generic maximum the corridor would veto anyway. */
            out.has_target = 0;
        }
    }
    else if (recto_sign != 0) {
        int recto_k = -1, behind_k = -1;
        double recto_abs = 1e30, behind_abs = 1e30;
        for (int k = 1; k < nsteps - 1; k++) {
            int ko = k - recto_sign; /* sample toward the sheet / verso */
            int ki = k + recto_sign; /* sample toward the inward recto gap */
            if (prof[ko] < 0.0 || prof[ki] < 0.0) continue;
            double drop = prof[ko] - prof[ki];
            if (prof[ko] < effective_bright || drop < effective_gain) continue;
            double t = -reach + (2.0 * reach * (double)k) / (double)(nsteps - 1);
            double inward_t = (double)recto_sign * t;
            if (bidirectional_recto) {
                if (fabs(t) < recto_abs) { recto_abs = fabs(t); recto_k = k; }
            } else if (inward_t >= -ridge_lock) {
                if (fabs(t) < recto_abs) { recto_abs = fabs(t); recto_k = k; }
            } else if (fabs(t) < behind_abs) {
                behind_abs = fabs(t); behind_k = k;
            }
        }
        if (recto_k < 0) recto_k = behind_k;
        if (recto_k >= 0) {
            double t = -reach + (2.0 * reach * (double)recto_k) / (double)(nsteps - 1);
            out.target_t = t;
            out.target_value = prof[recto_k];
            {
                int ko = recto_k - recto_sign;
                int ki = recto_k + recto_sign;
                out.transition_drop = prof[ko] - prof[ki];
                out.used_local_contrast = prof[ko] < bright_min ||
                                          out.transition_drop < min_gain;
            }
            out.has_target = 1;
            out.used_recto = 1;
            out.used_outward = (double)recto_sign*t < -ridge_lock;
            out.on_ridge = fabs(t) <= ridge_lock;
        }
    }
    if (out.has_target) {
        for (int k = 0; k < nsteps; k++) {
            double t = -reach + (2.0 * reach * (double)k) / (double)(nsteps - 1);
            if (fabs(t - out.target_t) <= 2.0 + 1e-9 &&
                prof[k] > out.visible_target)
                out.visible_target = prof[k];
        }
    }
    return out;
}

/* Route bake-dark evidence without making darkness itself a geometry rule.  A
 * reported alternative must be a supported recto transition strictly inward
 * beyond the ordinary ridge lock.  A crack or ink patch with no such alternative
 * remains an honest-dark candidate for the later image-space classifier. */
static int dark_inward_recto_target(const RidgeProfile *rp, int recto_sign,
                                    double dark_max, double bright_min,
                                    double min_gain, double ridge_lock) {
    if (!rp || recto_sign == 0 || !rp->has_target || !rp->used_recto)
        return 0;
    if ((double)recto_sign * rp->target_t <= ridge_lock) return 0;
    if (rp->visible_current < 0.0 || rp->visible_current > dark_max) return 0;
    if (rp->visible_target < bright_min) return 0;
    return rp->visible_target - rp->visible_current >= min_gain;
}

/* Midline sibling of dark_inward_recto_target: a dark site reports a band
 * center as its alternative regardless of side -- the midline is symmetric,
 * so "strictly inward" is replaced by "not already centered". */
static int dark_band_center_target(const RidgeProfile *rp, int oriented_sign,
                                   double dark_max, double bright_min,
                                   double min_gain, double ridge_lock) {
    if (!rp || oriented_sign == 0 || !rp->has_target ||
        !(rp->used_midline || rp->used_band_fallback))
        return 0;
    if (fabs(rp->target_t) <= ridge_lock) return 0;
    if (rp->visible_current < 0.0 || rp->visible_current > dark_max) return 0;
    if (rp->visible_target < bright_min) return 0;
    return rp->visible_target - rp->visible_current >= min_gain;
}

/* Return +1 when +n is inward, -1 when -n is inward, or 0 when the radial
 * orientation is too ambiguous to label the sheet's faces.  Normals and points
 * are z/y/x; the scroll axis is parallel to z. */
static int inward_profile_sign(const QuadStripSnapOpts *o,
                               const double p[3], const double n[3]) {
    if (!snap_oriented(o)) return 0;
    double ry = p[1] - o->axis_y, rx = p[2] - o->axis_x;
    double rr = sqrt(ry*ry + rx*rx);
    double nn = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    if (rr < 1.0 || nn < 0.5) return 0;
    /* dot(+n, inward unit radial); require a modest orientation margin. */
    double d = -(n[1]*ry + n[2]*rx) / (nn*rr);
    if (fabs(d) < 0.15) return 0;
    return d > 0.0 ? 1 : -1;
}

static int normalize3(const double in[3], double out[3]) {
    double q = sqrt(in[0]*in[0] + in[1]*in[1] + in[2]*in[2]);
    if (q < 1e-12) return 0;
    out[0] = in[0]/q; out[1] = in[1]/q; out[2] = in[2]/q;
    return 1;
}

/* A guided correction is allowed to be wider than the historical fixed
 * three-voxel tube, but it must not change winding identity.  Measure the
 * adjacent-wrap spacing directly on the initializer: unwrap atan2 along each
 * structured row, match phase +/-2*pi, and retain 45% of the smaller local
 * radial separation.  Even if both neighbouring turns move toward one another,
 * a positive ordering gap remains.  Limits are always anchored to `orig`, not
 * to the iterated position, so repeated proximal rounds cannot ratchet across a
 * wrap. */
typedef struct WindingPoint {
    double phase, radius;
} WindingPoint;

static int winding_point_cmp(const void *aa,const void *bb) {
    const WindingPoint *a=(const WindingPoint *)aa;
    const WindingPoint *b=(const WindingPoint *)bb;
    return a->phase<b->phase?-1:(a->phase>b->phase?1:0);
}

static int winding_double_cmp(const void *aa,const void *bb) {
    double a=*(const double *)aa,b=*(const double *)bb;
    return a<b?-1:(a>b?1:0);
}

static size_t winding_lower_bound(const WindingPoint *p,size_t n,double phase) {
    size_t lo=0,hi=n;
    while(lo<hi){size_t mid=lo+(hi-lo)/2;if(p[mid].phase<phase)lo=mid+1;else hi=mid;}
    return lo;
}

static double winding_quantile(const double *v,size_t n,double q) {
    if(!n)return 0.0;
    double x=q*(double)(n-1);size_t i=(size_t)floor(x);
    double t=x-(double)i;
    return i+1<n?(1.0-t)*v[i]+t*v[i+1]:v[i];
}

static float *adaptive_winding_limits(
        const float *orig,int H,int W,double axis_y,double axis_x,
        double minimum_tube,double maximum_reach,size_t *out_supported,
        double out_percentile[3]) {
    const double TWO_PI=6.28318530717958647693;
    size_t nv=(size_t)H*(size_t)W;
    float *limit=(float *)malloc(nv*sizeof *limit);
    double *pitch=(double *)malloc(nv*sizeof *pitch);
    double *valid=(double *)malloc(nv*sizeof *valid);
    double *phase=(double *)malloc((size_t)W*sizeof *phase);
    double *radius=(double *)malloc((size_t)W*sizeof *radius);
    WindingPoint *sorted=(WindingPoint *)malloc((size_t)W*sizeof *sorted);
    if(!limit||!pitch||!valid||!phase||!radius||!sorted) {
        free(limit);free(pitch);free(valid);free(phase);free(radius);free(sorted);
        return NULL;
    }
    for(size_t v=0;v<nv;v++)pitch[v]=NAN;
    size_t nvalid=0;
    for(int r=0;r<H;r++) {
        for(int c=0;c<W;c++) {
            size_t v=(size_t)r*(size_t)W+(size_t)c;
            double y=(double)orig[v*3+1]-axis_y;
            double x=(double)orig[v*3+2]-axis_x;
            double a=atan2(x,y);
            radius[c]=hypot(y,x);
            if(c==0)phase[c]=a;
            else {
                double prev=atan2((double)orig[(v-1)*3+2]-axis_x,
                                  (double)orig[(v-1)*3+1]-axis_y);
                double d=a-prev;
                while(d>3.14159265358979323846)d-=TWO_PI;
                while(d<-3.14159265358979323846)d+=TWO_PI;
                phase[c]=phase[c-1]+d;
            }
            sorted[c].phase=phase[c];sorted[c].radius=radius[c];
        }
        qsort(sorted,(size_t)W,sizeof *sorted,winding_point_cmp);
        for(int c=0;c<W;c++) {
            double best_pitch=1e300,best_error=1e300;
            for(int side=-1;side<=1;side+=2) {
                double target=phase[c]+(double)side*TWO_PI;
                size_t j=winding_lower_bound(sorted,(size_t)W,target);
                for(int dj=-1;dj<=0;dj++) {
                    if((dj<0&&j==0)||(dj==0&&j==(size_t)W))continue;
                    size_t k=dj<0?j-1:j;
                    double error=fabs(sorted[k].phase-target);
                    double p=fabs(sorted[k].radius-radius[c]);
                    if(error<=0.15&&p>=2.0&&p<=256.0&&
                       (error<best_error-1e-12||
                        (fabs(error-best_error)<=1e-12&&p<best_pitch))) {
                        best_error=error;best_pitch=p;
                    }
                }
            }
            if(best_pitch<1e299) {
                size_t v=(size_t)r*(size_t)W+(size_t)c;
                pitch[v]=best_pitch;valid[nvalid++]=best_pitch;
            }
        }
    }
    double fallback=minimum_tube>0.0?minimum_tube/0.45:maximum_reach/0.45;
    if(nvalid) {
        qsort(valid,nvalid,sizeof *valid,winding_double_cmp);
        fallback=winding_quantile(valid,nvalid,0.50);
    }
    size_t nsupported=0;
    for(size_t v=0;v<nv;v++) {
        double p=isfinite(pitch[v])?pitch[v]:fallback;
        if(isfinite(pitch[v]))nsupported++;
        double g=0.45*p;
        if(g<minimum_tube&&p>2.0*minimum_tube/0.9)g=minimum_tube;
        if(g>maximum_reach)g=maximum_reach;
        if(g<0.25)g=0.25;
        limit[v]=(float)g;valid[v]=g;
    }
    qsort(valid,nv,sizeof *valid,winding_double_cmp);
    if(out_percentile) {
        out_percentile[0]=winding_quantile(valid,nv,0.05);
        out_percentile[1]=winding_quantile(valid,nv,0.50);
        out_percentile[2]=winding_quantile(valid,nv,0.95);
    }
    if(out_supported)*out_supported=nsupported;
    free(pitch);free(valid);free(phase);free(radius);free(sorted);
    return limit;
}

static double ridge_candidate_score(const RidgeProfile *rp) {
    if (!rp || !rp->has_target) return -1e30;
    /* Genuine two-edged band centers earn the oriented bonus alongside recto
     * transitions; merged/clipped-band fallbacks compete on brightness and
     * distance alone so real centers outrank smoothed maxima. */
    return ((rp->used_recto || rp->used_midline) ? 1000.0 : 0.0) +
           2.0 * rp->transition_drop + rp->visible_target -
           1.5 * fabs(rp->target_t);
}

static int ridge_candidate_inside_winding(
        const QuadStripSnapOpts *o,const double reference[3],
        double winding_limit,const double p[3],const double dir[3],
        const RidgeProfile *rp) {
    if(!rp||!rp->has_target)return 0;
    double candidate[3];
    for(int k=0;k<3;k++)candidate[k]=p[k]+rp->target_t*dir[k];
    if(o->preserve_axial&&fabs(candidate[0]-reference[0])>1e-4)return 0;
    double d2=0.0;for(int k=0;k<3;k++){double d=candidate[k]-reference[k];d2+=d*d;}
    if(o->guided_reach>0.0&&d2>o->guided_reach*o->guided_reach+1e-9)return 0;
    if(snap_oriented(o)&&winding_limit>0.0) {
        double r0=hypot(reference[1]-o->axis_y,reference[2]-o->axis_x);
        double r1=hypot(candidate[1]-o->axis_y,candidate[2]-o->axis_x);
        if(fabs(r1-r0)>winding_limit+1e-9)return 0;
    }
    return 1;
}

/* Search a small family of directions only after the ordinary normal profile
 * is unresolved or ambiguous.  Every direction is oriented so positive depth
 * is inward, which makes the vector quilt coherent even when mesh normals have
 * inconsistent signs.  Wider candidates must remain in a radial tube around
 * the predicted winding phase; with the adaptive corridor (45% of the locally
 * measured adjacent-wrap pitch -- ~23 vox on this dataset, so p50 corridor
 * ~10 vox) they cannot jump to the neighbouring wrap. */
static RidgeProfile guided_ridge_profile(
        CubeTable *ct, const QuadStripSnapOpts *o,
        const double p[3], const double n_in[3],
        double bright_min, double min_gain, double ridge_lock,
        double guided_dark_max, int force_guided,
        const double reference[3],double winding_limit,
        double out_dir[3], int *out_recto_sign, int *out_winding_rejected) {
    RidgeProfile best;
    double n[3];
    memset(&best, 0, sizeof best);
    best.current = best.visible_current = best.target_value =
        best.visible_target = best.nearest_peak_value = -1.0;
    best.nearest_peak_t = 1e30;
    if (!normalize3(n_in, n)) {
        out_dir[0]=out_dir[1]=out_dir[2]=0.0;
        if (out_recto_sign) *out_recto_sign = 0;
        return best;
    }
    if(o->preserve_axial) {
        n[0]=0.0;
        if(!normalize3(n,n)) {
            double radial[3]={0.0,o->axis_y-p[1],o->axis_x-p[2]};
            if(!normalize3(radial,n)) {
                out_dir[0]=out_dir[1]=out_dir[2]=0.0;
                if(out_recto_sign)*out_recto_sign=0;
                return best;
            }
        }
    }
    int rsign = inward_profile_sign(o, p, n);
    double base_dir[3] = {n[0],n[1],n[2]};
    int base_recto = 0;
    if (rsign != 0) {
        for (int k=0;k<3;k++) base_dir[k] *= (double)rsign;
        base_recto = 1;
    }
    /* Midline mode scans the whole band, which at the measured ~23-voxel
     * pitch can be thicker than the default 6-voxel reach; widen toward the
     * corridor while never exceeding the guided reach. */
    double base_reach = o->reach;
    if (o->snap_side == QUAD_SNAP_SIDE_MIDLINE && base_recto &&
        winding_limit > 0.0) {
        double widened = winding_limit + 2.0;
        if (widened > o->guided_reach) widened = o->guided_reach;
        if (widened > base_reach) base_reach = widened;
    }
    best = ridge_profile(ct,p,base_dir,base_reach,o->step,bright_min,
                         min_gain,ridge_lock,base_recto,o->local_contrast,
                         force_guided,o->snap_side,winding_limit);
    if(best.has_target&&!ridge_candidate_inside_winding(
            o,reference,winding_limit,p,base_dir,&best)) {
        best.has_target=0;
        if(out_winding_rejected)(*out_winding_rejected)++;
    }
    for (int k=0;k<3;k++) out_dir[k]=base_dir[k];
    if (out_recto_sign) *out_recto_sign=base_recto;

    /* A dim local shoulder is not sufficient evidence to stop searching.  The
     * old early return accepted any recto-labelled transition, including a
     * low-contrast resin crack, and consequently never inspected the brighter
     * same-wrap surface a few voxels away.  Raster-dark triangle interiors also
     * force this coherent alternative search even when their corner samples are
     * individually bright. */
    int profile_dark = guided_dark_max >= 0.0 && best.visible_current >= 0.0 &&
                       best.visible_current <= guided_dark_max;
    if (!snap_oriented(o) || winding_limit <= 0.0 ||
        o->guided_reach <= o->reach + 1e-9 ||
        (!force_guided && !profile_dark)) return best;

    double ry=p[1]-o->axis_y, rx=p[2]-o->axis_x;
    double radius=sqrt(ry*ry+rx*rx);
    if (radius < 1.0) return best;
    double inward[3]={0.0,-ry/radius,-rx/radius};
    double aligned[3]={n[0],n[1],n[2]};
    if (aligned[0]*inward[0]+aligned[1]*inward[1]+aligned[2]*inward[2] < 0.0)
        for(int k=0;k<3;k++)aligned[k] = -aligned[k];
    const double beta[3]={0.35,0.65,1.0};
    double best_score=ridge_candidate_score(&best);
    for(int d=0;d<3;d++) {
        double rawdir[3]={
            (1.0-beta[d])*aligned[0]+beta[d]*inward[0],
            (1.0-beta[d])*aligned[1]+beta[d]*inward[1],
            (1.0-beta[d])*aligned[2]+beta[d]*inward[2]
        };
        if(o->preserve_axial)rawdir[0]=0.0;
        double dir[3]; if(!normalize3(rawdir,dir))continue;
        RidgeProfile rp=ridge_profile(ct,p,dir,o->guided_reach,o->step,
            bright_min,min_gain,ridge_lock,1,o->local_contrast,force_guided,
            o->snap_side,winding_limit);
        int oriented_hit=rp.used_recto||rp.used_midline||rp.used_band_fallback;
        double lock_t=o->snap_side==QUAD_SNAP_SIDE_MIDLINE
            ?fabs(rp.target_t):rp.target_t;
        if(!rp.has_target||!oriented_hit||
           (!force_guided&&lock_t<=ridge_lock))continue;
        if(!ridge_candidate_inside_winding(
                o,reference,winding_limit,p,dir,&rp)) {
            if(out_winding_rejected)(*out_winding_rejected)++;
            continue;
        }
        double score=ridge_candidate_score(&rp);
        int coherent_upgrade = force_guided && rp.visible_target >= bright_min &&
            rp.transition_drop + 5.0 >= best.transition_drop &&
            rp.visible_target + 5.0 >= best.visible_target;
        if(score>best_score+1e-9 || coherent_upgrade) {
            best=rp;best.used_guided=1;
            best.used_wide=fabs(best.target_t)>o->reach+1e-9;
            for(int k=0;k<3;k++)out_dir[k]=dir[k];
            if(out_recto_sign)*out_recto_sign=1;
            best_score=score;
        }
    }
    return best;
}

static int snap_one_vertex(CubeTable *ct, const double p[3], const double n[3],
                           double reach, double step, double bright_min,
                           double min_gain, double ridge_lock, double *out_t) {
    RidgeProfile rp = ridge_profile(ct, p, n, reach, step, bright_min,
                                    min_gain, ridge_lock, 0, 0, 0,
                                    QUAD_SNAP_SIDE_NEAREST, 0.0);
    if (!rp.has_target) return 0;
    *out_t = rp.target_t;
    return 1;
}

static int hist_percentile(const size_t hist[256], size_t count, double pct) {
    if (count == 0) return 0;
    if (pct < 0.0) pct = 0.0;
    if (pct > 100.0) pct = 100.0;
    size_t rank = (size_t)floor((pct / 100.0) * (double)(count - 1));
    size_t acc = 0;
    for (int i = 0; i < 256; i++) {
        acc += hist[i];
        if (acc > rank) return i;
    }
    return 255;
}

typedef struct SnapBakeCounts {
    size_t fitted_samples, filled_samples;
    size_t fitted_dark, filled_dark, missing;
} SnapBakeCounts;

static int fixed_window_u8(double value, double lo, double hi) {
    double q = (value-lo)*255.0/(hi-lo);
    if(q<0.0)q=0.0;if(q>255.0)q=255.0;
    return (int)lround(q);
}

/* A write-free form of the production bake metric for the regular ribbon.
 * Pixel centres are sampled in the same two emitted triangles and with the
 * same +/-2 normal maximum as rawtex_bake.  The ordinary grid is 2x1 vox per
 * cell, hence exactly two samples per quad; integer custom spacings remain
 * supported. */
static SnapBakeCounts structured_bake_scan(
        CubeTable *ct, const float *verts, const float *normals,
        const uint8_t *filled, const float *uv, int H, int W,
        double lo, double hi, int dark_u8, uint8_t *dark_vertices,
        uint8_t *cell_lit, uint8_t *cell_dark) {
    SnapBakeCounts out={0};
    if(!ct||!verts||!normals||!filled||!uv||H<2||W<2||hi<=lo)return out;
    if(dark_vertices)memset(dark_vertices,0,(size_t)H*(size_t)W);
    if(cell_dark)memset(cell_dark,0,(size_t)(H-1)*(size_t)(W-1));
    double du=fabs((double)uv[2]-(double)uv[0]);
    double dv=fabs((double)uv[(size_t)W*2+1]-(double)uv[1]);
    int nu=(int)lround(du), nv=(int)lround(dv);
    if(nu<1)nu=1;if(nu>16)nu=16;if(nv<1)nv=1;if(nv>16)nv=16;
    size_t missing=0,filled_samples=0,filled_dark=0,fitted_samples=0,
           fitted_dark=0;
    /* Even quad rows then odd quad rows: same-phase rows touch disjoint
     * vertex-row pairs, so the idempotent dark_vertices stores never race and
     * CubeTable reads are documented thread-safe after the prewarm. */
    for(int phase=0;phase<2;phase++){
        int r;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,4) \
    reduction(+:missing,filled_samples,filled_dark,fitted_samples,fitted_dark)
#endif
        for(r=phase;r<H-1;r+=2)for(int c=0;c+1<W;c++){
        size_t a=(size_t)r*(size_t)W+(size_t)c;
        size_t b=a+1, cc=a+(size_t)W, d=cc+1;
        int fill_cell=(int)filled[a]+(int)filled[b]+(int)filled[cc]+(int)filled[d]>=2;
        int all_fill=filled[a]&&filled[b]&&filled[cc]&&filled[d];
        int cell_samples=0,cell_dark_samples=0;
        for(int iy=0;iy<nv;iy++)for(int ix=0;ix<nu;ix++){
            double x=((double)ix+0.5)/(double)nu;
            double y=((double)iy+0.5)/(double)nv;
            size_t id[3];double w[3];
            if(x>=y){
                id[0]=a;id[1]=b;id[2]=d;
                w[0]=1.0-x;w[1]=x-y;w[2]=y;
            }else{
                id[0]=a;id[1]=d;id[2]=cc;
                w[0]=1.0-y;w[1]=x;w[2]=y-x;
            }
            float p[3]={0,0,0},n[3]={0,0,0};
            for(int q=0;q<3;q++)for(int k=0;k<3;k++){
                p[k]+=(float)(w[q]*(double)verts[id[q]*3+k]);
                n[k]+=(float)(w[q]*(double)normals[id[q]*3+k]);
            }
            double nn=sqrt((double)n[0]*n[0]+(double)n[1]*n[1]+(double)n[2]*n[2]);
            if(nn>1e-6){for(int k=0;k<3;k++)n[k]=(float)((double)n[k]/nn);}
            double value=sample_vertex(ct,p,nn>0.5?n:NULL,2.0,5);
            if(value<0.0){missing++;continue;}
            int dark=fixed_window_u8(value,lo,hi)<=dark_u8;
            cell_samples++;cell_dark_samples+=dark;
            if(!dark&&cell_lit!=NULL)
                cell_lit[(size_t)r*(size_t)(W-1)+(size_t)c]=1;
            if(dark&&dark_vertices)
                for(int q=0;q<3;q++)dark_vertices[id[q]]=1;
            if(fill_cell){filled_samples++;if(dark)filled_dark++;}
            else{fitted_samples++;if(dark)fitted_dark++;}
        }
        /* A raster-dark continuation core is conservative: all four vertices
         * must lack source texture and at least 75% of the exact bake samples
         * must be dark.  Partly fitted seam cells belong to the later collar,
         * never to the movable core. */
        if(cell_dark&&all_fill&&cell_samples>0&&
           4*cell_dark_samples>=3*cell_samples)
            cell_dark[(size_t)r*(size_t)(W-1)+(size_t)c]=1;
        }
    }
    out.missing=missing;out.filled_samples=filled_samples;
    out.filled_dark=filled_dark;out.fitted_samples=fitted_samples;
    out.fitted_dark=fitted_dark;
    return out;
}

static SnapBakeCounts structured_bake_counts(
        CubeTable *ct, const float *verts, const float *normals,
        const uint8_t *filled, const float *uv, int H, int W,
        double lo, double hi, int dark_u8) {
    return structured_bake_scan(ct,verts,normals,filled,uv,H,W,
                                lo,hi,dark_u8,NULL,NULL,NULL);
}

typedef struct SnapFaceTargetCounts {
    size_t dark_samples, candidates, outward_candidates, no_candidate;
} SnapFaceTargetCounts;

/* Search at the locations that actually emit dark pixels.  The historical
 * snap searched only the three triangle corners, even though a triangle with
 * three bright corners can cut through dark resin between them.  For every
 * fixed-window-dark raster sample, find a supported recto target on either
 * side of the bad initializer and distribute that PHYSICAL vector to its
 * movable triangle corners. */
static SnapFaceTargetCounts structured_dark_face_targets(
        CubeTable *ct,const QuadStripSnapOpts *o,
        const float *verts,const float *orig,const float *normals,
        const uint8_t *filled,const uint8_t *frozen,const float *uv,
        const float *winding_limit,int H,int W,double ridge_min,
        float *offset_sum,float *weight_sum,int *winding_rejected,
        uint8_t *cell_candidate) {
    SnapFaceTargetCounts out={0};
    size_t nverts=(size_t)H*(size_t)W;
    memset(offset_sum,0,nverts*3*sizeof *offset_sum);
    memset(weight_sum,0,nverts*sizeof *weight_sum);
    if(!ct||!o||!verts||!orig||!normals||!filled||!frozen||!uv||
       !winding_limit||H<2||W<2)return out;
    double du=fabs((double)uv[2]-(double)uv[0]);
    double dv=fabs((double)uv[(size_t)W*2+1]-(double)uv[1]);
    int nu=(int)lround(du),nv=(int)lround(dv);
    if(nu<1)nu=1;if(nu>16)nu=16;if(nv<1)nv=1;if(nv>16)nv=16;
    double dark_raw=o->bake_window_low+
        ((double)o->bake_dark_u8/255.0)*(o->bake_window_high-o->bake_window_low);
    size_t dark_samples=0,candidates=0,outward_candidates=0,no_candidate=0;
    int wrej=0;
    /* Two-phase quad-row coloring: same-phase rows scatter into disjoint
     * vertex-row pairs, so the corner accumulators never race.  A vertex row
     * still receives its two contributing cell rows in a fixed phase order,
     * so results are run-to-run deterministic (the accumulation order differs
     * from the serial build only across phases, a last-ulp effect). */
    for(int phase=0;phase<2;phase++){
        int r;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,2) \
    reduction(+:dark_samples,candidates,outward_candidates,no_candidate,wrej)
#endif
        for(r=phase;r<H-1;r+=2)for(int c=0;c+1<W;c++) {
        size_t a=(size_t)r*(size_t)W+(size_t)c;
        size_t b=a+1,cc=a+(size_t)W,d=cc+1;
        for(int iy=0;iy<nv;iy++)for(int ix=0;ix<nu;ix++) {
            double x=((double)ix+0.5)/(double)nu;
            double y=((double)iy+0.5)/(double)nv;
            size_t id[3];double w[3];
            if(x>=y){id[0]=a;id[1]=b;id[2]=d;
                w[0]=1.0-x;w[1]=x-y;w[2]=y;}
            else{id[0]=a;id[1]=d;id[2]=cc;
                w[0]=1.0-y;w[1]=x;w[2]=y-x;}
            double p[3]={0,0,0},reference[3]={0,0,0},n[3]={0,0,0};
            double guard=0.0;
            for(int q=0;q<3;q++)for(int k=0;k<3;k++) {
                p[k]+=w[q]*(double)verts[id[q]*3+(size_t)k];
                reference[k]+=w[q]*(double)orig[id[q]*3+(size_t)k];
                n[k]+=w[q]*(double)normals[id[q]*3+(size_t)k];
            }
            for(int q=0;q<3;q++)guard+=w[q]*(double)winding_limit[id[q]];
            double nn=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
            if(nn<1e-6)continue;
            for(int k=0;k<3;k++)n[k]/=nn;
            float pf[3]={(float)p[0],(float)p[1],(float)p[2]};
            float nf[3]={(float)n[0],(float)n[1],(float)n[2]};
            double visible=sample_vertex(ct,pf,nf,2.0,5);
            if(visible<0.0||fixed_window_u8(visible,o->bake_window_low,
                                           o->bake_window_high)>o->bake_dark_u8)
                continue;
            dark_samples++;
            double dir[3];int rsign=0,wr_local=0;
            RidgeProfile rp=guided_ridge_profile(
                ct,o,p,n,ridge_min,o->min_gain,o->ridge_lock,dark_raw,1,
                reference,guard,dir,&rsign,&wr_local);
            wrej+=wr_local;
            /* Compare against the ACTUAL bake-normal sample.  A radial guided
             * profile can see the target within its own +/-2 support even while
             * the emitted triangle remains black along its geometric normal. */
            int supported=rp.has_target&&rp.visible_target>dark_raw+1e-6&&
                rp.visible_target-visible>=fmax(4.0,0.20*o->min_gain);
            if(supported&&cell_candidate!=NULL)
                cell_candidate[(size_t)r*(size_t)(W-1)+(size_t)c]=1;
            if(!supported){no_candidate++;continue;}
            int movable=0;for(int q=0;q<3;q++)movable+=!frozen[id[q]];
            if(!movable){no_candidate++;continue;}
            double delta[3]={rp.target_t*dir[0],rp.target_t*dir[1],
                             rp.target_t*dir[2]};
            if(o->preserve_axial)delta[0]=0.0;
            for(int q=0;q<3;q++)if(!frozen[id[q]]) {
                /* A face translation should be visible to every movable
                 * corner, including the low-barycentric corner. */
                double evidence=0.25+w[q];
                weight_sum[id[q]]+=(float)evidence;
                for(int k=0;k<3;k++)
                    offset_sum[id[q]*3+(size_t)k]+=(float)(evidence*delta[k]);
            }
            candidates++;
            if(rp.used_outward)outward_candidates++;
        }
        }
    }
    out.dark_samples=dark_samples;out.candidates=candidates;
    out.outward_candidates=outward_candidates;out.no_candidate=no_candidate;
    if(winding_rejected)*winding_rejected+=wrej;
    return out;
}

/* Sparse bounded topology uses the same two-triangle lattice but compact vertex
 * indexing.  One centroid sample per emitted triangle preserves the fixed-window
 * transaction gate and dark-corner seeding without pretending absent alpha
 * regions are a rectangular grid. */
static SnapBakeCounts sparse_bake_scan(
        CubeTable *ct,const float *verts,const float *normals,
        const int32_t *faces,size_t nf,const uint8_t *filled,
        size_t nv,double lo,double hi,int dark_u8,uint8_t *dark_vertices,
        uint8_t *dark_faces) {
    SnapBakeCounts out={0};
    if(!ct||!verts||!normals||!faces||!filled||hi<=lo)return out;
    if(dark_vertices)memset(dark_vertices,0,nv);
    if(dark_faces)memset(dark_faces,0,nf);
    for(size_t f=0;f<nf;f++) {
        size_t id[3]={(size_t)faces[f*3],(size_t)faces[f*3+1],
                      (size_t)faces[f*3+2]};
        float p[3]={0,0,0},n[3]={0,0,0};
        int fill_cell=(int)filled[id[0]]+(int)filled[id[1]]+
                      (int)filled[id[2]]>=2;
        for(int q=0;q<3;q++)for(int k=0;k<3;k++) {
            p[k]+=verts[id[q]*3+(size_t)k]/3.0f;
            n[k]+=normals[id[q]*3+(size_t)k]/3.0f;
        }
        double nn=sqrt((double)n[0]*n[0]+(double)n[1]*n[1]+(double)n[2]*n[2]);
        if(nn>1e-6)for(int k=0;k<3;k++)n[k]=(float)((double)n[k]/nn);
        double value=sample_vertex(ct,p,nn>0.5?n:NULL,2.0,5);
        if(value<0.0){out.missing++;continue;}
        int dark=fixed_window_u8(value,lo,hi)<=dark_u8;
        if(dark&&dark_vertices)for(int q=0;q<3;q++)dark_vertices[id[q]]=1;
        if(dark_faces&&dark&&filled[id[0]]&&filled[id[1]]&&filled[id[2]])
            dark_faces[f]=1;
        if(fill_cell){out.filled_samples++;if(dark)out.filled_dark++;}
        else{out.fitted_samples++;if(dark)out.fitted_dark++;}
    }
    return out;
}

typedef struct RemeshAudit {
    double displacement_rms;
    double tangential_p95, tangential_max;
    double sander_p95, sander_p99, sander_max;
    double quad_sander_p95, quad_sander_p99, quad_sander_max;
    size_t triangles, quads, over2, over4, flipped, degenerate;
    size_t quad_over2, quad_over4, quad_flipped, quad_degenerate;
    int advised;
} RemeshAudit;

/* Match rawtex_bake's u8 singular-value diagnostic and sander_tint's decoder.
 * The bake deliberately has a finite representable range; transaction gates
 * should see the same tail that the review PNG sees, not a different continuous
 * approximation. */
static double bake_diag_sigma(double sigma) {
    double t;
    int code;
    if (sigma < 1e-9) code = 0;
    else {
        t = 128.0 + 42.5 * (log(sigma) / log(2.0));
        if (t < 0.0) t = 0.0;
        if (t > 255.0) t = 255.0;
        code = (int)(t + 0.5);
    }
    return exp2(((double)code - 128.0) / 42.5);
}

static double triangle_sander_bake(const float *verts, const float *uv,
                                   size_t a, size_t b, size_t c,
                                   int *out_degenerate) {
    double s1=(double)uv[a*2],t1=(double)uv[a*2+1];
    double s2=(double)uv[b*2],t2=(double)uv[b*2+1];
    double s3=(double)uv[c*2],t3=(double)uv[c*2+1];
    double A2=(s2-s1)*(t3-t1)-(s3-s1)*(t2-t1);
    double smax=1e6,smin=0.0;
    int deg=fabs(A2)<1e-12;
    if (!deg) {
        double Ss[3],St[3],aa=0.0,bb=0.0,cc=0.0;
        for(int k=0;k<3;k++) {
            double q1=(double)verts[a*3+k],q2=(double)verts[b*3+k];
            double q3=(double)verts[c*3+k];
            Ss[k]=(q1*(t2-t3)+q2*(t3-t1)+q3*(t1-t2))/A2;
            St[k]=(q1*(s3-s2)+q2*(s1-s3)+q3*(s2-s1))/A2;
            aa+=Ss[k]*Ss[k];bb+=Ss[k]*St[k];cc+=St[k]*St[k];
        }
        double disc=sqrt(fmax(0.0,(aa-cc)*(aa-cc)+4.0*bb*bb));
        double lmax=0.5*(aa+cc+disc),lmin=0.5*(aa+cc-disc);
        smax=sqrt(fmax(0.0,lmax));
        smin=sqrt(fmax(0.0,lmin));
        if (lmin <= 1e-12) deg=1;
    }
    if(out_degenerate)*out_degenerate=deg;
    smax=bake_diag_sigma(smax);smin=bake_diag_sigma(smin);
    {
        double fwd=sqrt(0.5*(smax*smax+smin*smin));
        double inv=sqrt(0.5*(1.0/(smax*smax)+1.0/(smin*smin)));
        return fwd>inv?fwd:inv;
    }
}

static void sander_hist_add(size_t hist[512],double D) {
    int b=D>0.0?(int)floor(log(D)/log(2.0)*64.0):0;
    if(b<0)b=0;if(b>511)b=511;hist[b]++;
}

static void sander_hist_percentiles(const size_t hist[512],size_t count,
                                    double *p95,double *p99) {
    size_t want95=count?(size_t)ceil(0.95*(double)count):0;
    size_t want99=count?(size_t)ceil(0.99*(double)count):0;
    size_t acc=0;int b95=0,b99=0,got95=0;
    if(!count){*p95=*p99=0.0;return;}
    for(int b=0;b<512;b++) {
        acc+=hist[b];
        if(!got95&&acc>=want95){b95=b;got95=1;}
        if(acc>=want99){b99=b;break;}
    }
    *p95=exp2(((double)b95+0.5)/64.0);
    *p99=exp2(((double)b99+0.5)/64.0);
}

static double triangle_cross_dot(const float *verts,const float *orig,
                                 size_t a,size_t b,size_t c) {
    double ce0[3],ce1[3],oe0[3],oe1[3],cn[3],on[3];
    for(int k=0;k<3;k++) {
        ce0[k]=(double)verts[b*3+k]-verts[a*3+k];
        ce1[k]=(double)verts[c*3+k]-verts[a*3+k];
        oe0[k]=(double)orig[b*3+k]-orig[a*3+k];
        oe1[k]=(double)orig[c*3+k]-orig[a*3+k];
    }
    cn[0]=ce0[1]*ce1[2]-ce0[2]*ce1[1];
    cn[1]=ce0[2]*ce1[0]-ce0[0]*ce1[2];
    cn[2]=ce0[0]*ce1[1]-ce0[1]*ce1[0];
    on[0]=oe0[1]*oe1[2]-oe0[2]*oe1[1];
    on[1]=oe0[2]*oe1[0]-oe0[0]*oe1[2];
    on[2]=oe0[0]*oe1[1]-oe0[1]*oe1[0];
    return cn[0]*on[0]+cn[1]*on[1]+cn[2]*on[2];
}

/* Diagnose when the old quad sampling has ceased to be a good material
 * parameterization.  Normal motion alone is expected from CT snap; accumulated
 * TANGENTIAL motion is what makes vertices exchange their natural sampling
 * neighborhoods.  The primary Jacobian audit scores the exact two triangles
 * emitted for each structured quad and round-trips their singular values through
 * obj_bake_raw's diagnostic encoding before applying the same two-sided Sander
 * L2 metric used by the review PNG.  The former quad-average remains diagnostic.
 * This is intentionally a gate only: the next alternating stage will rebuild a
 * bounded snip and transfer its fit/recto confidence when `advised` trips. */
static RemeshAudit remesh_audit(const float *verts, const float *orig,
                                const float *base_nrm, const float *uv,
                                const int32_t *faces,size_t nf,
                                size_t nv, int H, int W,
                                double motion_trigger, double motion_max_trigger,
                                double sander_trigger, double over4_trigger) {
    RemeshAudit a;
    memset(&a, 0, sizeof a);
    size_t mh[512] = {0}, dh[512] = {0}, qdh[512] = {0};
    double move2 = 0.0, cell_scale = 1.0;
    if (nv > 1) {
        double du = 0.0, dv = 0.0;
        if (W > 1) du = fabs((double)uv[2] - (double)uv[0]);
        if (H > 1 && (size_t)W < nv)
            dv = fabs((double)uv[(size_t)W*2+1] - (double)uv[1]);
        if (du > 1e-6 && dv > 1e-6) cell_scale = du < dv ? du : dv;
        else if (du > 1e-6) cell_scale = du;
        else if (dv > 1e-6) cell_scale = dv;
        else if (faces && nf) {
            double best=1e30;
            for(size_t f=0;f<nf;f++)for(int e=0;e<3;e++) {
                size_t ia=(size_t)faces[f*3+(size_t)e];
                size_t ib=(size_t)faces[f*3+(size_t)((e+1)%3)];
                double ds=(double)uv[ia*2]-uv[ib*2];
                double dt=(double)uv[ia*2+1]-uv[ib*2+1];
                double d=sqrt(ds*ds+dt*dt);
                if(d>1e-6&&d<best)best=d;
            }
            if(best<1e29)cell_scale=best;
        }
    }
    for (size_t v = 0; v < nv; v++) {
        double d[3] = { (double)verts[v*3+0] - orig[v*3+0],
                        (double)verts[v*3+1] - orig[v*3+1],
                        (double)verts[v*3+2] - orig[v*3+2] };
        double n[3] = { base_nrm[v*3+0], base_nrm[v*3+1], base_nrm[v*3+2] };
        double dn = d[0]*n[0] + d[1]*n[1] + d[2]*n[2];
        double d2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
        double tang = sqrt(fmax(0.0, d2 - dn*dn));
        move2 += d2;
        if (tang > a.tangential_max) a.tangential_max = tang;
        int b = (int)floor(100.0 * tang / cell_scale);
        if (b < 0) b = 0; if (b > 511) b = 511;
        mh[b]++;
    }
    a.displacement_rms = nv ? sqrt(move2 / (double)nv) : 0.0;
    {
        size_t want = nv ? (size_t)ceil(0.95 * (double)nv) : 0, acc = 0;
        int b = 0;
        for (; b < 512; b++) { acc += mh[b]; if (acc >= want) break; }
        if (b > 511) b = 511;
        a.tangential_p95 = ((double)b + 0.5) * 0.01 * cell_scale;
    }

    if (H > 1 && W > 1 && (size_t)H * (size_t)W == nv) {
        for (int r = 0; r + 1 < H; r++) for (int c = 0; c + 1 < W; c++) {
            size_t i00=(size_t)r*W+c, i01=i00+1, i10=i00+W, i11=i10+1;
            double du = fabs(0.5*((double)uv[i01*2]-(double)uv[i00*2]
                                 +(double)uv[i11*2]-(double)uv[i10*2]));
            double dv = fabs(0.5*((double)uv[i10*2+1]-(double)uv[i00*2+1]
                                 +(double)uv[i11*2+1]-(double)uv[i01*2+1]));
            if (du < 1e-9) du = 1.0; if (dv < 1e-9) dv = 1.0;
            double eu[3], ev[3], ou[3], ov[3];
            for (int k = 0; k < 3; k++) {
                eu[k] = 0.5*((double)verts[i01*3+k]-verts[i00*3+k]
                            +(double)verts[i11*3+k]-verts[i10*3+k])/du;
                ev[k] = 0.5*((double)verts[i10*3+k]-verts[i00*3+k]
                            +(double)verts[i11*3+k]-verts[i01*3+k])/dv;
                ou[k] = 0.5*((double)orig[i01*3+k]-orig[i00*3+k]
                            +(double)orig[i11*3+k]-orig[i10*3+k]);
                ov[k] = 0.5*((double)orig[i10*3+k]-orig[i00*3+k]
                            +(double)orig[i11*3+k]-orig[i01*3+k]);
            }
            double g00=0.0,g01=0.0,g11=0.0;
            for (int k=0;k<3;k++){g00+=eu[k]*eu[k];g01+=eu[k]*ev[k];g11+=ev[k]*ev[k];}
            double tr=g00+g11, disc=sqrt(fmax(0.0,(g00-g11)*(g00-g11)+4.0*g01*g01));
            double lmax=0.5*(tr+disc), lmin=0.5*(tr-disc), D=1e9;
            if (lmin > 1e-12) {
                double fwd=sqrt(0.5*(lmax+lmin));
                double inv=sqrt(0.5*(1.0/lmax+1.0/lmin));
                D=fwd>inv?fwd:inv;
            } else a.quad_degenerate++;
            if (D > a.quad_sander_max) a.quad_sander_max=D;
            if (D >= 2.0) a.quad_over2++; if (D >= 4.0) a.quad_over4++;
            sander_hist_add(qdh,D);
            double cn[3]={eu[1]*ev[2]-eu[2]*ev[1],eu[2]*ev[0]-eu[0]*ev[2],eu[0]*ev[1]-eu[1]*ev[0]};
            double on[3]={ou[1]*ov[2]-ou[2]*ov[1],ou[2]*ov[0]-ou[0]*ov[2],ou[0]*ov[1]-ou[1]*ov[0]};
            if (cn[0]*on[0]+cn[1]*on[1]+cn[2]*on[2] <= 0.0) a.quad_flipped++;
            a.quads++;

            /* Score the exact (a,b,d) / (a,d,c) diagonal emitted by
             * QuadStrip_build and painted by rawtex_bake.  The averaged quad
             * above is retained only as a broad-drift diagnostic. */
            {
                size_t tri[2][3]={{i00,i01,i11},{i00,i11,i10}};
                for(int t=0;t<2;t++) {
                    int deg=0;
                    double td=triangle_sander_bake(verts,uv,tri[t][0],tri[t][1],tri[t][2],&deg);
                    if(td>a.sander_max)a.sander_max=td;
                    if(td>=2.0)a.over2++;if(td>=4.0)a.over4++;
                    if(deg)a.degenerate++;
                    if(triangle_cross_dot(verts,orig,tri[t][0],tri[t][1],tri[t][2])<=0.0)
                        a.flipped++;
                    sander_hist_add(dh,td);a.triangles++;
                }
            }
        }
        sander_hist_percentiles(dh,a.triangles,&a.sander_p95,&a.sander_p99);
        sander_hist_percentiles(qdh,a.quads,&a.quad_sander_p95,&a.quad_sander_p99);
    } else if (faces && nf) {
        for(size_t f=0;f<nf;f++) {
            size_t a0=(size_t)faces[f*3],b0=(size_t)faces[f*3+1];
            size_t c0=(size_t)faces[f*3+2];
            int deg=0;
            double td=triangle_sander_bake(verts,uv,a0,b0,c0,&deg);
            if(td>a.sander_max)a.sander_max=td;
            if(td>=2.0)a.over2++;if(td>=4.0)a.over4++;
            if(deg)a.degenerate++;
            if(triangle_cross_dot(verts,orig,a0,b0,c0)<=0.0)a.flipped++;
            sander_hist_add(dh,td);a.triangles++;
        }
        sander_hist_percentiles(dh,a.triangles,&a.sander_p95,&a.sander_p99);
    }
    a.advised = (motion_trigger > 0.0 && a.tangential_p95 >= motion_trigger*cell_scale)
             || (motion_max_trigger > 0.0 && a.tangential_max >= motion_max_trigger*cell_scale)
             || (sander_trigger > 0.0 && a.sander_p95 >= sander_trigger)
             || (over4_trigger > 0.0 && a.triangles > 0
                 && (double)a.over4/(double)a.triangles >= over4_trigger)
             || a.flipped > 0 || a.degenerate > 0;
    return a;
}

/* Build one monotone reparameterization shared by every curve along an axis.
 * Independent per-row arc maps created visible horizontal shear bands: the
 * same material feature landed in a different column on each neighboring row.
 * The ordinary metric averages edge length over the orthogonal axis.  The
 * robust alternative winsorizes at two standard deviations before averaging,
 * preventing a few crack-spanning edges from reallocating columns across an
 * otherwise coherent full-height strip. */
static double remesh_axis_edge_length(const float *src,int H,int W,int along_u,
                                      int edge_index,int line) {
    size_t stride=along_u?1u:(size_t)W;
    size_t base=along_u?(size_t)line*W:(size_t)line;
    const float *a=src+(base+(size_t)edge_index*stride)*3;
    const float *b=src+(base+(size_t)(edge_index+1)*stride)*3;
    double dz=(double)b[0]-a[0],dy=(double)b[1]-a[1],dx=(double)b[2]-a[2];
    (void)H;
    return sqrt(dz*dz+dy*dy+dx*dx);
}

static void remesh_axis_map(const float *src, int H, int W, int along_u,
                            int robust, double *map, double *edge, double *cum) {
    int n=along_u?W:H,lines=along_u?H:W;
    for(int i=0;i+1<n;i++) {
        double sum=0.0,sum2=0.0;
        for(int l=0;l<lines;l++) {
            double d=remesh_axis_edge_length(src,H,W,along_u,i,l);
            sum+=d;sum2+=d*d;
        }
        edge[i]=sum/(double)lines;
        if(robust&&lines>=8) {
            double mean=edge[i];
            double sigma=sqrt(fmax(0.0,sum2/(double)lines-mean*mean));
            double lo=fmax(0.0,mean-2.0*sigma),hi=mean+2.0*sigma,rsum=0.0;
            for(int l=0;l<lines;l++) {
                double d=remesh_axis_edge_length(src,H,W,along_u,i,l);
                if(d<lo)d=lo;if(d>hi)d=hi;rsum+=d;
            }
            edge[i]=rsum/(double)lines;
        }
    }
    cum[0]=0.0;for(int i=1;i<n;i++)cum[i]=cum[i-1]+edge[i-1];
    if(cum[n-1]<1e-9){for(int i=0;i<n;i++)map[i]=(double)i;return;}
    int k=0;
    for(int i=0;i<n;i++){
        double target=cum[n-1]*(double)i/(double)(n-1);
        while(k+1<n-1&&cum[k+1]<target)k++;
        double den=cum[k+1]-cum[k],t=den>1e-12?(target-cum[k])/den:0.0;
        map[i]=(double)k+t;
    }
}

static void remesh_curve_map(const float *src,const uint8_t *src_fill,
                             float *dst,uint8_t *dst_fill,int n,
                             size_t stride,const double *map){
    for(int i=0;i<n;i++){
        int k=(int)floor(map[i]);if(k<0)k=0;if(k>n-2)k=n-2;
        double t=map[i]-(double)k;if(i==n-1){k=n-2;t=1.0;}
        const float *a=src+(size_t)k*stride*3,*b=src+(size_t)(k+1)*stride*3;
        float *q=dst+(size_t)i*stride*3;
        for(int ch=0;ch<3;ch++)q[ch]=(float)((1.0-t)*(double)a[ch]+t*(double)b[ch]);
        int j=t<0.5?k:k+1;dst_fill[(size_t)i*stride]=src_fill[(size_t)j*stride];
    }
}

typedef struct RemeshQuality {
    double p95,p99,max,edge_log_rms;
    size_t triangles,over2,over4,flipped,degenerate;
} RemeshQuality;

static double remesh_edge_log_rms(const float *verts,const float *uv,int H,int W) {
    double sum=0.0,sum2=0.0;size_t n=0;
    for(int r=0;r<H;r++)for(int c=0;c<W;c++) {
        size_t a=(size_t)r*W+c;
        size_t nb[2]={a+1,a+(size_t)W};
        int count=(c+1<W)+(r+1<H);
        if(c+1>=W&&r+1<H)nb[0]=nb[1];
        for(int j=0;j<count;j++) {
            size_t b=nb[j];double l2=0.0,d2=0.0;
            for(int k=0;k<3;k++){double d=(double)verts[b*3+k]-verts[a*3+k];l2+=d*d;}
            for(int k=0;k<2;k++){double d=(double)uv[b*2+k]-uv[a*2+k];d2+=d*d;}
            double q=log(fmax(1e-12,sqrt(l2))/fmax(1e-12,sqrt(d2)));
            sum+=q;sum2+=q*q;n++;
        }
    }
    if(!n)return 0.0;
    {double mean=sum/(double)n;return sqrt(fmax(0.0,sum2/(double)n-mean*mean));}
}

static RemeshQuality remesh_triangle_quality(const float *verts,const float *reference,
                                             const float *uv,int H,int W) {
    RemeshQuality q;size_t hist[512]={0};memset(&q,0,sizeof q);
    for(int r=0;r+1<H;r++)for(int c=0;c+1<W;c++) {
        size_t i00=(size_t)r*W+c,i01=i00+1,i10=i00+W,i11=i10+1;
        size_t tri[2][3]={{i00,i01,i11},{i00,i11,i10}};
        for(int t=0;t<2;t++) {
            int deg=0;double d=triangle_sander_bake(verts,uv,tri[t][0],tri[t][1],tri[t][2],&deg);
            if(d>q.max)q.max=d;if(d>=2.0)q.over2++;if(d>=4.0)q.over4++;
            if(deg)q.degenerate++;
            if(reference&&triangle_cross_dot(verts,reference,tri[t][0],tri[t][1],tri[t][2])<=0.0)
                q.flipped++;
            sander_hist_add(hist,d);q.triangles++;
        }
    }
    sander_hist_percentiles(hist,q.triangles,&q.p95,&q.p99);
    q.edge_log_rms=remesh_edge_log_rms(verts,uv,H,W);
    return q;
}

static int remesh_candidate_safe(const float *cand,const float *base,int H,int W) {
    for(int r=0;r+1<H;r++)for(int c=0;c+1<W;c++) {
        size_t i00=(size_t)r*W+c,i01=i00+1,i10=i00+W,i11=i10+1;
        size_t tri[2][3]={{i00,i01,i11},{i00,i11,i10}};
        for(int t=0;t<2;t++) {
            size_t a=tri[t][0],b=tri[t][1],d=tri[t][2];
            double ce0[3],ce1[3],oe0[3],oe1[3],cn[3],on[3];
            for(int k=0;k<3;k++) {
                ce0[k]=(double)cand[b*3+k]-cand[a*3+k];
                ce1[k]=(double)cand[d*3+k]-cand[a*3+k];
                oe0[k]=(double)base[b*3+k]-base[a*3+k];
                oe1[k]=(double)base[d*3+k]-base[a*3+k];
            }
            cn[0]=ce0[1]*ce1[2]-ce0[2]*ce1[1];
            cn[1]=ce0[2]*ce1[0]-ce0[0]*ce1[2];
            cn[2]=ce0[0]*ce1[1]-ce0[1]*ce1[0];
            on[0]=oe0[1]*oe1[2]-oe0[2]*oe1[1];
            on[1]=oe0[2]*oe1[0]-oe0[0]*oe1[2];
            on[2]=oe0[0]*oe1[1]-oe0[1]*oe1[0];
            double ref2=on[0]*on[0]+on[1]*on[1]+on[2]*on[2];
            if(ref2>1e-18&&cn[0]*on[0]+cn[1]*on[1]+cn[2]*on[2]<=0.10*ref2)return 0;
        }
    }
    for(int r=0;r<H;r++)for(int c=0;c<W;c++) {
        size_t a=(size_t)r*W+c;
        if(c+1<W) {
            size_t b=a+1;double x=0.0,o=0.0;
            for(int k=0;k<3;k++){double dc=(double)cand[b*3+k]-cand[a*3+k];double db=(double)base[b*3+k]-base[a*3+k];x+=dc*dc;o+=db*db;}
            if(o>1e-18&&(x<0.25*o||x>4.0*o))return 0;
        }
        if(r+1<H) {
            size_t b=a+(size_t)W;double x=0.0,o=0.0;
            for(int k=0;k<3;k++){double dc=(double)cand[b*3+k]-cand[a*3+k];double db=(double)base[b*3+k]-base[a*3+k];x+=dc*dc;o+=db*db;}
            if(o>1e-18&&(x<0.25*o||x>4.0*o))return 0;
        }
    }
    return 1;
}

static int remesh_quality_admissible(const RemeshQuality *cand,const RemeshQuality *base) {
    const double eps=1e-12;
    if(cand->flipped||cand->degenerate)return 0;
    if(cand->over4>base->over4||cand->over2>base->over2)return 0;
    if(cand->p95>base->p95+eps||cand->p99>base->p99+eps||cand->max>base->max+eps)return 0;
    return cand->over4<base->over4||cand->over2<base->over2||
           cand->p95+eps<base->p95||cand->p99+eps<base->p99||cand->max+eps<base->max||
           cand->edge_log_rms<base->edge_log_rms*0.995;
}

static int remesh_quality_better(const RemeshQuality *a,const RemeshQuality *b) {
    const double eps=1e-12;
    if(a->over4!=b->over4)return a->over4<b->over4;
    if(a->over2!=b->over2)return a->over2<b->over2;
    if(fabs(a->p99-b->p99)>eps)return a->p99<b->p99;
    if(fabs(a->p95-b->p95)>eps)return a->p95<b->p95;
    if(fabs(a->max-b->max)>eps)return a->max<b->max;
    return a->edge_log_rms<b->edge_log_rms;
}

static void remesh_blend_map(double *map,int n,double alpha) {
    for(int i=0;i<n;i++)map[i]=(double)i+alpha*(map[i]-(double)i);
}

static void remesh_make_candidate(const float *base,const uint8_t *base_fill,
                                  float *tmp,uint8_t *tmp_fill,
                                  float *cand,uint8_t *cand_fill,
                                  int H,int W,int robust,int v_first,double alpha,
                                  double *map,double *edge,double *cum) {
    if(!v_first) {
        remesh_axis_map(base,H,W,1,robust,map,edge,cum);remesh_blend_map(map,W,alpha);
        for(int r=0;r<H;r++)remesh_curve_map(base+(size_t)r*W*3,base_fill+(size_t)r*W,
                                             tmp+(size_t)r*W*3,tmp_fill+(size_t)r*W,W,1,map);
        remesh_axis_map(tmp,H,W,0,robust,map,edge,cum);remesh_blend_map(map,H,alpha);
        for(int c=0;c<W;c++)remesh_curve_map(tmp+(size_t)c*3,tmp_fill+(size_t)c,
                                             cand+(size_t)c*3,cand_fill+(size_t)c,H,(size_t)W,map);
    } else {
        remesh_axis_map(base,H,W,0,robust,map,edge,cum);remesh_blend_map(map,H,alpha);
        for(int c=0;c<W;c++)remesh_curve_map(base+(size_t)c*3,base_fill+(size_t)c,
                                             tmp+(size_t)c*3,tmp_fill+(size_t)c,H,(size_t)W,map);
        remesh_axis_map(tmp,H,W,1,robust,map,edge,cum);remesh_blend_map(map,W,alpha);
        for(int r=0;r<H;r++)remesh_curve_map(tmp+(size_t)r*W*3,tmp_fill+(size_t)r*W,
                                             cand+(size_t)r*W*3,cand_fill+(size_t)r*W,W,1,map);
    }
}

int QuadStripSnap_remesh_structured(float *verts, uint8_t *filled,
                                    const float *uv, int H, int W, int sweeps,
                                    QuadStripRemeshStats *out_stats) {
    if(out_stats)memset(out_stats,0,sizeof *out_stats);
    if(!verts||!filled||!uv||H<2||W<2||sweeps<1)return -1;
    size_t nv=(size_t)H*(size_t)W;
    if((size_t)H>SIZE_MAX/(size_t)W||nv>SIZE_MAX/(3*sizeof(float)))return -1;
    float *base=(float *)malloc(nv*3*sizeof *base);
    float *tmp=(float *)malloc(nv*3*sizeof *tmp);
    float *cand=(float *)malloc(nv*3*sizeof *cand);
    float *best=(float *)malloc(nv*3*sizeof *best);
    uint8_t *base_fill=(uint8_t *)malloc(nv),*tf=(uint8_t *)malloc(nv);
    uint8_t *cf=(uint8_t *)malloc(nv),*best_fill=(uint8_t *)malloc(nv);
    int maxn=H>W?H:W;
    double *work=(double *)malloc((size_t)maxn*3*sizeof *work);
    if(!base||!tmp||!cand||!best||!base_fill||!tf||!cf||!best_fill||!work){
        free(base);free(tmp);free(cand);free(best);free(base_fill);free(tf);free(cf);
        free(best_fill);free(work);return -1;
    }
    double *map=work,*edge=work+maxn,*cum=work+2*maxn;
    static const double alpha[]={1.0,0.5,0.25,0.125,0.0625,0.03125};
    RemeshQuality initial={0},last={0};int have_initial=0,applied=0;
    for(int sweep=0;sweep<sweeps;sweep++) {
        memcpy(base,verts,nv*3*sizeof *base);memcpy(base_fill,filled,nv);
        RemeshQuality before=remesh_triangle_quality(base,base,uv,H,W),bestq={0};
        int have_best=0,best_robust=0,best_v_first=0;double best_alpha=0.0;
        if(!have_initial){initial=before;have_initial=1;}
        for(int robust=0;robust<2;robust++)for(int v_first=0;v_first<2;v_first++)
            for(size_t ai=0;ai<sizeof alpha/sizeof alpha[0];ai++) {
                remesh_make_candidate(base,base_fill,tmp,tf,cand,cf,H,W,robust,v_first,
                                      alpha[ai],map,edge,cum);
                if(!remesh_candidate_safe(cand,base,H,W))continue;
                RemeshQuality cq=remesh_triangle_quality(cand,base,uv,H,W);
                if(!remesh_quality_admissible(&cq,&before))continue;
                if(!have_best||remesh_quality_better(&cq,&bestq)) {
                    memcpy(best,cand,nv*3*sizeof *best);memcpy(best_fill,cf,nv);
                    bestq=cq;have_best=1;best_robust=robust;
                    best_v_first=v_first;best_alpha=alpha[ai];
                }
            }
        if(!have_best){last=before;break;}
        memcpy(verts,best,nv*3*sizeof *verts);memcpy(filled,best_fill,nv);
        last=bestq;applied++;
        if(out_stats){out_stats->used_robust_metric=best_robust;
                      out_stats->used_v_first=best_v_first;out_stats->alpha=best_alpha;}
    }
    if(!have_initial)initial=remesh_triangle_quality(verts,verts,uv,H,W);
    if(!applied)last=initial;
    if(out_stats){
        out_stats->applied_sweeps=applied;out_stats->before_p95=initial.p95;
        out_stats->before_p99=initial.p99;out_stats->after_p95=last.p95;
        out_stats->after_p99=last.p99;out_stats->before_edge_log_rms=initial.edge_log_rms;
        out_stats->after_edge_log_rms=last.edge_log_rms;out_stats->triangles=last.triangles;
        out_stats->before_over2=initial.over2;out_stats->before_over4=initial.over4;
        out_stats->after_over2=last.over2;out_stats->after_over4=last.over4;
    }
    free(base);free(tmp);free(cand);free(best);free(base_fill);free(tf);free(cf);
    free(best_fill);free(work);return 0;
}

static int snap_exact_intersections(const float *verts,size_t nv,
                                    const int32_t *faces,size_t nf,
                                    size_t *face_conflict_degree,
                                    IntersectionCleanupStats *out) {
    IntersectionCleanupParams p;
    IntersectionCleanup_default_params(&p);
    p.gap_max=1.0;p.parallel_angle_deg=20.0;
    p.max_conflicts=nf>1000000u?nf:1000000u;p.include_hinges=0;
    return IntersectionCleanup_audit(verts,nv,faces,nf,NULL,&p,
                                     face_conflict_degree,out);
}

/* Mark emitted triangles whose current normal has crossed, or nearly crossed,
 * the initializer normal.  The same 5% oriented-area margin is used by the
 * position solver's fine/coarse barriers.  Returning the face count separately
 * lets the exact collision transaction reject a locally feathered candidate
 * before the outer publication gate has to discard the entire snap pass. */
static size_t snap_orientation_violations(const float *verts,const float *orig,
                                          size_t nv,const int32_t *faces,size_t nf,
                                          float *motion_scale,
                                          size_t *newly_zeroed) {
    size_t bad=0,nnew=0;
    for(size_t f=0;f<nf;f++) {
        int32_t ia=faces[f*3],ib=faces[f*3+1],ic=faces[f*3+2];
        if(ia<0||ib<0||ic<0||(size_t)ia>=nv||(size_t)ib>=nv||(size_t)ic>=nv) {
            bad++;continue;
        }
        size_t a=(size_t)ia,b=(size_t)ib,c=(size_t)ic;
        double oe0[3],oe1[3],on[3];
        for(int k=0;k<3;k++) {
            oe0[k]=(double)orig[b*3+(size_t)k]-orig[a*3+(size_t)k];
            oe1[k]=(double)orig[c*3+(size_t)k]-orig[a*3+(size_t)k];
        }
        on[0]=oe0[1]*oe1[2]-oe0[2]*oe1[1];
        on[1]=oe0[2]*oe1[0]-oe0[0]*oe1[2];
        on[2]=oe0[0]*oe1[1]-oe0[1]*oe1[0];
        double ref2=on[0]*on[0]+on[1]*on[1]+on[2]*on[2];
        if(ref2>1e-18&&triangle_cross_dot(verts,orig,a,b,c)>0.05*ref2)continue;
        bad++;
        if(motion_scale) {
            size_t id[3]={a,b,c};
            for(int k=0;k<3;k++)if(motion_scale[id[k]]>0.0f) {
                motion_scale[id[k]]=0.0f;nnew++;
            }
        }
    }
    if(newly_zeroed)*newly_zeroed=nnew;
    return bad;
}

/* --------------------------------------------------------------------------
 * Connected dark-patch recovery
 *
 * The ordinary snap is intentionally pointwise-conservative.  Its remaining
 * failure mode is correspondingly non-pointwise: a whole harmonic patch can
 * occupy the gap beside the correct winding, so every short profile either
 * abstains or finds a locally plausible but globally inconsistent sheet.  The
 * recovery below treats an exact-bake-dark connected component as one decision:
 *
 *   1. grow a small lattice collar around the dark cells;
 *   2. ask whether all certified vertices on that collar belong to one
 *      4-connected, XYZ-continuous source region;
 *   3. use a coherent rim's strict ray majority to choose its winding side; if
 *      the rim is ambiguous, require the WHOLE patch to choose one side;
 *   4. screen the agreed ray targets through a local Laplacian displacement
 *      solve.  Certified vertices and everything outside the collar are zero-
 *      displacement Dirichlet anchors, giving the requested smooth join from
 *      recovered unknown area into known papyrus.
 *
 * This routine only builds a proposal.  QuadStripSnap_run commits it later as
 * one bake-improving, orientation-preserving, exact-intersection transaction.
 * -------------------------------------------------------------------------- */

typedef struct PatchRegion {
    int32_t label;
    size_t cells;
    int r0,r1,c0,c1;              /* inclusive dark-cell bounds */
} PatchRegion;

typedef struct PatchRecoverRunStats {
    size_t regions,regions_small;
    size_t border_consensus,ray_consensus;
    size_t border_selected,ambiguous_selected;
    size_t wobbly,no_consensus,overlap_skipped;
    size_t proposed,applied,attenuated,geometry_rejected;
    size_t vertices,vertices_retained;
    size_t dark_before,dark_after;
    size_t fitted_dark_before,missing_before;
    int transaction_accepted;
    double retained_scale;
} PatchRecoverRunStats;

typedef struct PatchRayHit {
    int hit;
    double t;
} PatchRayHit;

static int patch_region_largest_first(const void *aa,const void *bb) {
    const PatchRegion *a=(const PatchRegion *)aa,*b=(const PatchRegion *)bb;
    if(a->cells>b->cells)return -1;
    if(a->cells<b->cells)return 1;
    return a->label<b->label?-1:(a->label>b->label);
}

static int patch_known_edge_coherent(const float *verts,const float *uv,
                                     size_t a,size_t b) {
    double dz=(double)verts[a*3]-(double)verts[b*3];
    double dy=(double)verts[a*3+1]-(double)verts[b*3+1];
    double dx=(double)verts[a*3+2]-(double)verts[b*3+2];
    double du=(double)uv[a*2]-(double)uv[b*2];
    double dv=(double)uv[a*2+1]-(double)uv[b*2+1];
    double rest=hypot(du,dv);
    /* A lattice connection that jumps many source voxels is not evidence that
     * two rim arcs occupy the same physical winding.  The generous floor keeps
     * ordinary stretched-but-valid source cells connected. */
    double limit=fmax(8.0,4.0*rest+1.0);
    return dz*dz+dy*dy+dx*dx<=limit*limit;
}

static int patch_label_trusted_regions(const float *verts,const uint8_t *filled,
                                       const float *uv,int H,int W,
                                       int32_t *label,size_t *queue) {
    size_t nv=(size_t)H*(size_t)W;
    for(size_t v=0;v<nv;v++)label[v]=-1;
    int32_t next=0;
    QuadAdj a={0};a.structured=1;a.H=H;a.W=W;
    for(size_t seed=0;seed<nv;seed++) {
        if(filled[seed]||label[seed]>=0)continue;
        if(next==INT32_MAX)return -1;
        size_t qh=0,qt=0;queue[qt++]=seed;label[seed]=next;
        while(qh<qt) {
            size_t v=queue[qh++];int deg=qadj_degree(&a,v);
            for(int e=0;e<deg;e++) {
                size_t j=qadj_neighbor(&a,v,e);
                if(filled[j]||label[j]>=0||
                   !patch_known_edge_coherent(verts,uv,v,j))continue;
                label[j]=next;queue[qt++]=j;
            }
        }
        next++;
    }
    return (int)next;
}

static int patch_label_dark_regions(const uint8_t *dark,int H,int W,
                                    int32_t *label,size_t *queue,
                                    PatchRegion **out_region,size_t *out_n) {
    int CH=H-1,CW=W-1;size_t nc=(size_t)CH*(size_t)CW;
    for(size_t k=0;k<nc;k++)label[k]=-1;
    PatchRegion *region=NULL;size_t nr=0,cap=0;
    int32_t next=0;
    for(size_t seed=0;seed<nc;seed++) {
        if(!dark[seed]||label[seed]>=0)continue;
        if(next==INT32_MAX){free(region);return -1;}
        if(nr==cap) {
            size_t ncap=cap?cap*2u:64u;
            PatchRegion *grown=(PatchRegion *)realloc(region,ncap*sizeof *grown);
            if(!grown){free(region);return -1;}
            region=grown;cap=ncap;
        }
        int sr=(int)(seed/(size_t)CW),sc=(int)(seed%(size_t)CW);
        PatchRegion pr={next,0,sr,sr,sc,sc};
        size_t qh=0,qt=0;queue[qt++]=seed;label[seed]=next;
        while(qh<qt) {
            size_t k=queue[qh++];int r=(int)(k/(size_t)CW),c=(int)(k%(size_t)CW);
            pr.cells++;if(r<pr.r0)pr.r0=r;if(r>pr.r1)pr.r1=r;
            if(c<pr.c0)pr.c0=c;if(c>pr.c1)pr.c1=c;
            const int dr[4]={-1,1,0,0},dc[4]={0,0,-1,1};
            for(int e=0;e<4;e++) {
                int rr=r+dr[e],cc=c+dc[e];
                if(rr<0||rr>=CH||cc<0||cc>=CW)continue;
                size_t j=(size_t)rr*(size_t)CW+(size_t)cc;
                if(!dark[j]||label[j]>=0)continue;
                label[j]=next;queue[qt++]=j;
            }
        }
        region[nr++]=pr;next++;
    }
    qsort(region,nr,sizeof *region,patch_region_largest_first);
    *out_region=region;*out_n=nr;return 0;
}

/* First supported material band on ONE open half-ray.  Unlike ridge_profile,
 * its gain is measured against the dark patch origin, not the midpoint of a
 * shifted bidirectional profile.  Returning the band midpoint matches the
 * production midline contract and is insensitive to a noisy single maximum. */
static PatchRayHit patch_first_half_ray(CubeTable *ct,const double p[3],
                                        const double dir[3],double reach,
                                        double step,double bright_min,
                                        double min_gain,int local_contrast) {
    PatchRayHit out={0,0.0};
    if(step<=0.0)step=0.5;if(reach<=step)return out;
    int n=(int)ceil(reach/step);if(n<2)n=2;if(n>256)n=256;
    double dt=reach/(double)n,prof[257];
    double pmin=1e30,pmax=-1.0;
    for(int k=0;k<=n;k++) {
        double t=dt*(double)k;
        prof[k]=sample_trilinear(ct,p[0]+t*dir[0],p[1]+t*dir[1],p[2]+t*dir[2]);
        if(prof[k]>=0.0){if(prof[k]<pmin)pmin=prof[k];if(prof[k]>pmax)pmax=prof[k];}
    }
    double I0=prof[0];if(I0<0.0||pmax<pmin)return out;
    double threshold=bright_min,gain=min_gain;
    if(local_contrast&&pmax-pmin>=12.0) {
        double span=pmax-pmin;
        double local_bright=pmin+0.55*span;
        double local_gain=fmax(8.0,0.20*span);
        if(local_bright<threshold)threshold=local_bright;
        if(local_gain<gain)gain=local_gain;
    }
    for(int k=1;k<=n;k++) {
        if(prof[k]<threshold||prof[k]-I0<gain)continue;
        int first=k,last=k;
        while(last+1<=n&&prof[last+1]>=threshold&&prof[last+1]-I0>=gain)last++;
        /* Reject a one-sample sparkle; real papyrus at 0.5-voxel sampling has
         * support on at least two consecutive samples. */
        if(last-first+1>=2) {
            out.hit=1;out.t=0.5*dt*(double)(first+last);return out;
        }
        k=last;
    }
    return out;
}

/* Orient +dir toward the scroll axis, so a patch-wide +/- vote has one physical
 * meaning despite arbitrary triangle-normal signs.  Near-tangent normals are
 * deliberately undecidable: enough of them make the patch report wobbly/no
 * consensus instead of silently choosing a winding. */
static int patch_inward_normal(const QuadStripSnapOpts *o,const float *verts,
                               const float *normals,size_t v,double dir[3]) {
    double raw[3]={normals[v*3],normals[v*3+1],normals[v*3+2]};
    if(o->preserve_axial)raw[0]=0.0;
    if(!normalize3(raw,dir))return 0;
    double ry=(double)verts[v*3+1]-o->axis_y;
    double rx=(double)verts[v*3+2]-o->axis_x,rr=hypot(ry,rx);
    if(rr<1.0)return 0;
    double d=-(dir[1]*ry+dir[2]*rx)/rr;
    if(fabs(d)<0.15)return 0;
    if(d<0.0)for(int k=0;k<3;k++)dir[k]=-dir[k];
    return 1;
}

/* Return +1/-1 for agreement, 0 for abstention.  A coherent trusted XYZ rim is
 * already the winding decision: its rays only choose which oriented side leads
 * back to that winding, so use their strict majority and tolerate noisy or
 * tangent samples.  Only an ambiguous rim invokes the expensive fallback from
 * the recovery contract, where every unknown core vertex must cast a decisive
 * vote and one contradictory sign makes the minimal patch wobbly. */
static int patch_global_side(size_t plus,size_t minus,size_t undecided,
                             size_t core_unknown,int border_consensus,
                             int *wobbly) {
    if(wobbly)*wobbly=0;
    if(border_consensus) {
        (void)undecided;(void)core_unknown;
        if(plus==minus)return 0;
        return plus>minus?1:-1;
    }
    if(plus&&minus){if(wobbly)*wobbly=1;return 0;}
    if(!plus&&!minus)return 0;
    if(plus+minus!=core_unknown||undecided)return 0;
    return plus?1:-1;
}

static void patch_clear_box(uint8_t *a,int H,int W,int r0,int r1,int c0,int c1) {
    (void)H;
    for(int r=r0;r<=r1;r++)memset(a+(size_t)r*(size_t)W+(size_t)c0,0,
                                  (size_t)(c1-c0+1));
}

static void patch_clear_box_float(float *a,int H,int W,int channels,
                                  int r0,int r1,int c0,int c1) {
    (void)H;
    for(int r=r0;r<=r1;r++)memset(a+((size_t)r*(size_t)W+(size_t)c0)*(size_t)channels,
        0,(size_t)(c1-c0+1)*(size_t)channels*sizeof *a);
}

static int patch_recover_propose(
        CubeTable *ct,const QuadStripSnapOpts *o,const float *verts,
        const float *normals,const uint8_t *filled,const float *uv,
        int H,int W,double ridge_min,float *proposal,uint8_t *moved,
        int32_t *moved_patch,PatchRecoverRunStats *stats) {
    size_t nv=(size_t)H*(size_t)W,nc=(size_t)(H-1)*(size_t)(W-1);
    memset(stats,0,sizeof *stats);memcpy(proposal,verts,nv*3*sizeof *proposal);
    memset(moved,0,nv);
    for(size_t v=0;v<nv;v++)moved_patch[v]=-1;
    uint8_t *cell_dark=(uint8_t *)calloc(nc,1);
    int32_t *cell_label=(int32_t *)malloc(nc*sizeof *cell_label);
    int32_t *trusted_label=(int32_t *)malloc(nv*sizeof *trusted_label);
    size_t *queue=(size_t *)malloc((nv>nc?nv:nc)*sizeof *queue);
    uint8_t *core=(uint8_t *)calloc(nv,1),*zone=(uint8_t *)calloc(nv,1);
    uint8_t *work=(uint8_t *)calloc(nv,1),*claimed=(uint8_t *)calloc(nv,1);
    uint8_t *hit=(uint8_t *)calloc(nv,1);
    float *ray_t=(float *)calloc(nv*2,sizeof *ray_t);
    float *direct=(float *)calloc(nv*3,sizeof *direct);
    float *disp=(float *)calloc(nv*3,sizeof *disp);
    PatchRegion *region=NULL;size_t nr=0;
    if(!cell_dark||!cell_label||!trusted_label||!queue||!core||!zone||!work||
       !claimed||!hit||!ray_t||!direct||!disp)goto oom;
    SnapBakeCounts bake=structured_bake_scan(ct,verts,normals,filled,uv,H,W,
        o->bake_window_low,o->bake_window_high,o->bake_dark_u8,
        NULL,NULL,cell_dark);
    stats->dark_before=bake.filled_dark;
    stats->fitted_dark_before=bake.fitted_dark;
    stats->missing_before=bake.missing;
    if(patch_label_dark_regions(cell_dark,H,W,cell_label,queue,&region,&nr)!=0)goto oom;
    if(patch_label_trusted_regions(verts,filled,uv,H,W,trusted_label,queue)<0)goto oom;
    stats->regions=nr;
    QuadAdj adj={0};adj.structured=1;adj.H=H;adj.W=W;
    for(size_t ri=0;ri<nr;ri++) {
        PatchRegion *pr=&region[ri];
        if(pr->cells<(size_t)o->patch_min_cells){stats->regions_small++;continue;}
        int r0=pr->r0,r1=pr->r1+1,c0=pr->c0,c1=pr->c1+1;
        for(int r=pr->r0;r<=pr->r1;r++)for(int c=pr->c0;c<=pr->c1;c++) {
            size_t ck=(size_t)r*(size_t)(W-1)+(size_t)c;
            if(cell_label[ck]!=pr->label)continue;
            size_t a=(size_t)r*(size_t)W+(size_t)c;
            core[a]=core[a+1]=core[a+(size_t)W]=core[a+(size_t)W+1]=1;
        }
        for(int r=r0;r<=r1;r++)for(int c=c0;c<=c1;c++) {
            size_t v=(size_t)r*(size_t)W+(size_t)c;zone[v]=core[v];
        }
        for(int ring=0;ring<o->patch_buffer;ring++) {
            int nr0=r0>0?r0-1:r0,nr1=r1+1<H?r1+1:r1;
            int nc0=c0>0?c0-1:c0,nc1=c1+1<W?c1+1:c1;
            for(int r=nr0;r<=nr1;r++)for(int c=nc0;c<=nc1;c++) {
                size_t v=(size_t)r*(size_t)W+(size_t)c;int on=zone[v]!=0;
                if(!on&&r>0)on=zone[v-(size_t)W]!=0;
                if(!on&&r+1<H)on=zone[v+(size_t)W]!=0;
                if(!on&&c>0)on=zone[v-1]!=0;
                if(!on&&c+1<W)on=zone[v+1]!=0;
                work[v]=(uint8_t)on;
            }
            for(int r=nr0;r<=nr1;r++)for(int c=nc0;c<=nc1;c++) {
                size_t v=(size_t)r*(size_t)W+(size_t)c;zone[v]=work[v];work[v]=0;
            }
            r0=nr0;r1=nr1;c0=nc0;c1=nc1;
        }

        int border_label=-1,border_mixed=0;size_t border_vertices=0;
        int overlap=0;
        for(int r=r0;r<=r1;r++)for(int c=c0;c<=c1;c++) {
            size_t v=(size_t)r*(size_t)W+(size_t)c;if(!zone[v])continue;
            if(filled[v]){if(claimed[v])overlap=1;continue;}
            int32_t lab=trusted_label[v];if(lab<0)continue;
            border_vertices++;
            if(border_label<0)border_label=lab;else if(border_label!=lab)border_mixed=1;
        }
        int border_consensus=border_vertices>=2&&border_label>=0&&!border_mixed;
        if(border_consensus)stats->border_consensus++;
        if(overlap) {
            stats->overlap_skipped++;
            patch_clear_box(core,H,W,r0,r1,c0,c1);patch_clear_box(zone,H,W,r0,r1,c0,c1);
            continue;
        }

        size_t plus=0,minus=0,undecided=0,core_unknown=0;
        for(int r=r0;r<=r1;r++)for(int c=c0;c<=c1;c++) {
            size_t v=(size_t)r*(size_t)W+(size_t)c;
            if(!core[v]||!filled[v])continue;
            core_unknown++;ray_t[v*2]=ray_t[v*2+1]=-1.0f;
            double dir[3];
            if(!patch_inward_normal(o,verts,normals,v,dir)){undecided++;continue;}
            double p[3]={verts[v*3],verts[v*3+1],verts[v*3+2]};
            PatchRayHit hp=patch_first_half_ray(ct,p,dir,o->patch_ray_reach,o->step,
                ridge_min,o->min_gain,o->local_contrast);
            double neg[3]={-dir[0],-dir[1],-dir[2]};
            PatchRayHit hm=patch_first_half_ray(ct,p,neg,o->patch_ray_reach,o->step,
                ridge_min,o->min_gain,o->local_contrast);
            if(hp.hit)ray_t[v*2]=(float)hp.t;if(hm.hit)ray_t[v*2+1]=(float)hm.t;
            if(hp.hit&&!hm.hit)plus++;
            else if(!hp.hit&&hm.hit)minus++;
            else if(hp.hit&&hm.hit&&hm.t-hp.t>=o->patch_ray_margin)plus++;
            else if(hp.hit&&hm.hit&&hp.t-hm.t>=o->patch_ray_margin)minus++;
            else undecided++;
        }
        int wobbly=0;
        int side=patch_global_side(plus,minus,undecided,core_unknown,
                                   border_consensus,&wobbly);
        if(wobbly)stats->wobbly++;
        if(!side) {
            if(!wobbly)stats->no_consensus++;
            patch_clear_box(core,H,W,r0,r1,c0,c1);patch_clear_box(zone,H,W,r0,r1,c0,c1);
            continue;
        }
        stats->ray_consensus++;
        if(border_consensus)stats->border_selected++;
        else stats->ambiguous_selected++;
        size_t nhit=0;
        for(int r=r0;r<=r1;r++)for(int c=c0;c<=c1;c++) {
            size_t v=(size_t)r*(size_t)W+(size_t)c;
            if(!zone[v]||!filled[v])continue;
            disp[v*3]=disp[v*3+1]=disp[v*3+2]=0.0f;
            if(!core[v])continue;
            float t=ray_t[v*2+(side<0)];if(t<0.0f)continue;
            double dir[3];if(!patch_inward_normal(o,verts,normals,v,dir))continue;
            for(int k=0;k<3;k++)direct[v*3+(size_t)k]=(float)((double)side*dir[k]*(double)t);
            if(o->preserve_axial)direct[v*3]=0.0f;
            hit[v]=1;nhit++;
        }
        if(!nhit) {
            stats->no_consensus++;
            patch_clear_box(core,H,W,r0,r1,c0,c1);patch_clear_box(zone,H,W,r0,r1,c0,c1);
            patch_clear_box(hit,H,W,r0,r1,c0,c1);
            patch_clear_box_float(direct,H,W,3,r0,r1,c0,c1);
            continue;
        }

        double tw=o->patch_target_weight,omega=0.85;
        for(int sweep=0;sweep<o->patch_relax_sweeps;sweep++)for(int color=0;color<2;color++)
            for(int r=r0;r<=r1;r++)for(int c=c0;c<=c1;c++) {
                if(((r+c)&1)!=color)continue;
                size_t v=(size_t)r*(size_t)W+(size_t)c;
                if(!zone[v]||!filled[v])continue;
                double rhs[3]={0,0,0};int deg=qadj_degree(&adj,v);
                for(int e=0;e<deg;e++) {
                    size_t j=qadj_neighbor(&adj,v,e);
                    if(zone[j]&&filled[j])for(int k=0;k<3;k++)rhs[k]+=disp[j*3+(size_t)k];
                    /* Known or outside-zone neighbours contribute the fixed
                     * zero displacement while still counting in deg. */
                }
                double den=(double)deg+(hit[v]?tw:0.0);
                if(hit[v])for(int k=0;k<3;k++)rhs[k]+=tw*(double)direct[v*3+(size_t)k];
                if(den<=0.0)continue;
                for(int k=0;k<3;k++) {
                    if(k==0&&o->preserve_axial){disp[v*3]=0.0f;continue;}
                    size_t q=v*3+(size_t)k;
                    disp[q]=(float)((1.0-omega)*(double)disp[q]+omega*rhs[k]/den);
                }
            }
        size_t changed=0;
        for(int r=r0;r<=r1;r++)for(int c=c0;c<=c1;c++) {
            size_t v=(size_t)r*(size_t)W+(size_t)c;
            if(!zone[v]||!filled[v])continue;
            double d2=0.0;for(int k=0;k<3;k++)d2+=(double)disp[v*3+(size_t)k]*disp[v*3+(size_t)k];
            if(d2<1e-8)continue;
            for(int k=0;k<3;k++)proposal[v*3+(size_t)k]=verts[v*3+(size_t)k]+disp[v*3+(size_t)k];
            if(o->preserve_axial)proposal[v*3]=verts[v*3];
            moved[v]=claimed[v]=1;moved_patch[v]=pr->label;changed++;
        }
        if(changed){stats->proposed++;stats->vertices+=changed;}
        patch_clear_box(core,H,W,r0,r1,c0,c1);patch_clear_box(zone,H,W,r0,r1,c0,c1);
        patch_clear_box(hit,H,W,r0,r1,c0,c1);
        patch_clear_box_float(ray_t,H,W,2,r0,r1,c0,c1);
        patch_clear_box_float(direct,H,W,3,r0,r1,c0,c1);
        patch_clear_box_float(disp,H,W,3,r0,r1,c0,c1);
    }
    free(cell_dark);free(cell_label);free(trusted_label);free(queue);free(core);free(zone);
    free(work);free(claimed);free(hit);free(ray_t);free(direct);free(disp);free(region);
    return 0;
oom:
    free(cell_dark);free(cell_label);free(trusted_label);free(queue);free(core);free(zone);
    free(work);free(claimed);free(hit);free(ray_t);free(direct);free(disp);free(region);
    return -1;
}

/* Compact ragged ribbons have the same patch semantics but no row-major vertex
 * identity.  Treat exact-bake-dark emitted triangles as the patch cells, grow
 * the collar through the actual face graph, and run the same whole-component
 * ray vote plus screened graph-Laplacian displacement solve.  This deliberately
 * uses only topology which is really emitted: alpha gaps cannot become phantom
 * neighbours merely because two vertices have nearby lattice coordinates. */
static int patch_label_trusted_graph(const float *verts,const uint8_t *filled,
                                     const float *uv,const QuadAdj *adj,size_t nv,
                                     int32_t *label,size_t *queue) {
    for(size_t v=0;v<nv;v++)label[v]=-1;
    int32_t next=0;
    for(size_t seed=0;seed<nv;seed++) {
        if(filled[seed]||label[seed]>=0)continue;
        if(next==INT32_MAX)return -1;
        size_t qh=0,qt=0;queue[qt++]=seed;label[seed]=next;
        while(qh<qt) {
            size_t v=queue[qh++];int deg=qadj_degree(adj,v);
            for(int e=0;e<deg;e++) {
                size_t j=qadj_neighbor(adj,v,e);
                if(j>=nv||filled[j]||label[j]>=0||
                   !patch_known_edge_coherent(verts,uv,v,j))continue;
                label[j]=next;queue[qt++]=j;
            }
        }
        next++;
    }
    return (int)next;
}

static int patch_face_incidence(const int32_t *faces,size_t nf,size_t nv,
                                size_t **out_off,int32_t **out_face) {
    if(nf>(size_t)INT32_MAX||nf>SIZE_MAX/3u)return -1;
    size_t *degree=(size_t *)calloc(nv,sizeof *degree);
    size_t *off=(size_t *)malloc((nv+1u)*sizeof *off);
    size_t *cursor=(size_t *)malloc(nv*sizeof *cursor);
    int32_t *inc=(int32_t *)malloc((nf?nf*3u:1u)*sizeof *inc);
    if(!degree||!off||!cursor||!inc){free(degree);free(off);free(cursor);free(inc);return -1;}
    for(size_t f=0;f<nf;f++)for(int k=0;k<3;k++) {
        int32_t q=faces[f*3u+(size_t)k];
        if(q<0||(size_t)q>=nv){free(degree);free(off);free(cursor);free(inc);return -1;}
        degree[(size_t)q]++;
    }
    off[0]=0;for(size_t v=0;v<nv;v++)off[v+1]=off[v]+degree[v];
    memcpy(cursor,off,nv*sizeof *cursor);
    for(size_t f=0;f<nf;f++)for(int k=0;k<3;k++) {
        size_t v=(size_t)faces[f*3u+(size_t)k];inc[cursor[v]++]=(int32_t)f;
    }
    free(degree);free(cursor);*out_off=off;*out_face=inc;return 0;
}

static int patch_label_dark_faces(const uint8_t *dark,const int32_t *faces,
                                  size_t nf,const size_t *voff,
                                  const int32_t *vface,int32_t *label,
                                  size_t *queue,PatchRegion **out_region,
                                  size_t *out_n) {
    for(size_t f=0;f<nf;f++)label[f]=-1;
    PatchRegion *region=NULL;size_t nr=0,cap=0;int32_t next=0;
    for(size_t seed=0;seed<nf;seed++) {
        if(!dark[seed]||label[seed]>=0)continue;
        if(next==INT32_MAX){free(region);return -1;}
        if(nr==cap) {
            size_t ncap=cap?cap*2u:64u;
            PatchRegion *grown=(PatchRegion *)realloc(region,ncap*sizeof *grown);
            if(!grown){free(region);return -1;}region=grown;cap=ncap;
        }
        PatchRegion pr={next,0,0,0,0,0};
        size_t qh=0,qt=0;queue[qt++]=seed;label[seed]=next;
        while(qh<qt) {
            size_t f=queue[qh++];pr.cells++;
            for(int k=0;k<3;k++) {
                size_t v=(size_t)faces[f*3u+(size_t)k];
                for(size_t p=voff[v];p<voff[v+1];p++) {
                    size_t j=(size_t)vface[p];
                    if(!dark[j]||label[j]>=0)continue;
                    label[j]=next;queue[qt++]=j;
                }
            }
        }
        region[nr++]=pr;next++;
    }
    qsort(region,nr,sizeof *region,patch_region_largest_first);
    *out_region=region;*out_n=nr;return 0;
}

static void patch_graph_clear(const size_t *list,size_t n,uint8_t *core,
                              uint8_t *zone,uint8_t *hit,float *ray_t,
                              float *direct,float *disp,float *next_disp) {
    for(size_t q=0;q<n;q++) {
        size_t v=list[q];core[v]=zone[v]=hit[v]=0;
        ray_t[v*2]=ray_t[v*2+1]=0.0f;
        for(int k=0;k<3;k++)direct[v*3+(size_t)k]=
            disp[v*3+(size_t)k]=next_disp[v*3+(size_t)k]=0.0f;
    }
}

static int patch_recover_graph_propose(
        CubeTable *ct,const QuadStripSnapOpts *o,const float *verts,
        const float *normals,const uint8_t *filled,const float *uv,
        const int32_t *faces,size_t nf,const QuadAdj *adj,size_t nv,
        double ridge_min,float *proposal,uint8_t *moved,
        int32_t *moved_patch,PatchRecoverRunStats *stats) {
    memset(stats,0,sizeof *stats);memcpy(proposal,verts,nv*3u*sizeof *proposal);
    memset(moved,0,nv);for(size_t v=0;v<nv;v++)moved_patch[v]=-1;
    uint8_t *face_dark=(uint8_t *)calloc(nf?nf:1u,1);
    int32_t *face_label=(int32_t *)malloc((nf?nf:1u)*sizeof *face_label);
    int32_t *trusted_label=(int32_t *)malloc(nv*sizeof *trusted_label);
    size_t *face_queue=(size_t *)malloc((nf?nf:1u)*sizeof *face_queue);
    size_t *vertex_queue=(size_t *)malloc(nv*sizeof *vertex_queue);
    size_t *voff=NULL;int32_t *vface=NULL;
    uint8_t *core=(uint8_t *)calloc(nv,1),*zone=(uint8_t *)calloc(nv,1);
    uint8_t *claimed=(uint8_t *)calloc(nv,1),*hit=(uint8_t *)calloc(nv,1);
    float *ray_t=(float *)calloc(nv*2u,sizeof *ray_t);
    float *direct=(float *)calloc(nv*3u,sizeof *direct);
    float *disp=(float *)calloc(nv*3u,sizeof *disp);
    float *next_disp=(float *)calloc(nv*3u,sizeof *next_disp);
    int32_t *head=NULL,*next_face=NULL;PatchRegion *region=NULL;size_t nr=0;
    if(!face_dark||!face_label||!trusted_label||!face_queue||!vertex_queue||
       !core||!zone||!claimed||!hit||!ray_t||!direct||!disp||!next_disp)goto oom;
    SnapBakeCounts bake=sparse_bake_scan(ct,verts,normals,faces,nf,filled,nv,
        o->bake_window_low,o->bake_window_high,o->bake_dark_u8,NULL,face_dark);
    stats->dark_before=bake.filled_dark;stats->fitted_dark_before=bake.fitted_dark;
    stats->missing_before=bake.missing;
    if(patch_face_incidence(faces,nf,nv,&voff,&vface)!=0)goto oom;
    if(patch_label_dark_faces(face_dark,faces,nf,voff,vface,face_label,
                              face_queue,&region,&nr)!=0)goto oom;
    if(patch_label_trusted_graph(verts,filled,uv,adj,nv,trusted_label,
                                 vertex_queue)<0)goto oom;
    stats->regions=nr;
    head=(int32_t *)malloc((nr?nr:1u)*sizeof *head);
    next_face=(int32_t *)malloc((nf?nf:1u)*sizeof *next_face);
    if(!head||!next_face)goto oom;
    for(size_t k=0;k<nr;k++)head[k]=-1;
    for(size_t f=0;f<nf;f++) {
        int32_t lab=face_label[f];next_face[f]=-1;
        if(lab<0)continue;next_face[f]=head[(size_t)lab];head[(size_t)lab]=(int32_t)f;
    }
    for(size_t ri=0;ri<nr;ri++) {
        PatchRegion *pr=&region[ri];
        if(pr->cells<(size_t)o->patch_min_cells){stats->regions_small++;continue;}
        size_t nzone=0;
        for(int32_t fi=head[(size_t)pr->label];fi>=0;fi=next_face[(size_t)fi])
            for(int k=0;k<3;k++) {
                size_t v=(size_t)faces[(size_t)fi*3u+(size_t)k];
                core[v]=1;if(!zone[v]){zone[v]=1;vertex_queue[nzone++]=v;}
            }
        size_t frontier0=0,frontier1=nzone;
        for(int ring=0;ring<o->patch_buffer;ring++) {
            for(size_t q=frontier0;q<frontier1;q++) {
                size_t v=vertex_queue[q];int deg=qadj_degree(adj,v);
                for(int e=0;e<deg;e++) {
                    size_t j=qadj_neighbor(adj,v,e);if(j>=nv||zone[j])continue;
                    zone[j]=1;vertex_queue[nzone++]=j;
                }
            }
            frontier0=frontier1;frontier1=nzone;
        }
        int border_label=-1,border_mixed=0,overlap=0;size_t border_vertices=0;
        for(size_t q=0;q<nzone;q++) {
            size_t v=vertex_queue[q];
            if(filled[v]){if(claimed[v])overlap=1;continue;}
            int32_t lab=trusted_label[v];if(lab<0)continue;
            border_vertices++;if(border_label<0)border_label=lab;
            else if(border_label!=lab)border_mixed=1;
        }
        int border_consensus=border_vertices>=2&&border_label>=0&&!border_mixed;
        if(border_consensus)stats->border_consensus++;
        if(overlap) {
            stats->overlap_skipped++;
            patch_graph_clear(vertex_queue,nzone,core,zone,hit,ray_t,direct,disp,next_disp);
            continue;
        }
        size_t plus=0,minus=0,undecided=0,core_unknown=0;
        for(size_t q=0;q<nzone;q++) {
            size_t v=vertex_queue[q];if(!core[v]||!filled[v])continue;
            core_unknown++;ray_t[v*2]=ray_t[v*2+1]=-1.0f;
            double dir[3];
            if(!patch_inward_normal(o,verts,normals,v,dir)){undecided++;continue;}
            double p[3]={verts[v*3],verts[v*3+1],verts[v*3+2]};
            PatchRayHit hp=patch_first_half_ray(ct,p,dir,o->patch_ray_reach,o->step,
                ridge_min,o->min_gain,o->local_contrast);
            double neg[3]={-dir[0],-dir[1],-dir[2]};
            PatchRayHit hm=patch_first_half_ray(ct,p,neg,o->patch_ray_reach,o->step,
                ridge_min,o->min_gain,o->local_contrast);
            if(hp.hit)ray_t[v*2]=(float)hp.t;if(hm.hit)ray_t[v*2+1]=(float)hm.t;
            if(hp.hit&&!hm.hit)plus++;else if(!hp.hit&&hm.hit)minus++;
            else if(hp.hit&&hm.hit&&hm.t-hp.t>=o->patch_ray_margin)plus++;
            else if(hp.hit&&hm.hit&&hp.t-hm.t>=o->patch_ray_margin)minus++;
            else undecided++;
        }
        int wobbly=0,side=patch_global_side(plus,minus,undecided,core_unknown,
                                            border_consensus,&wobbly);
        if(wobbly)stats->wobbly++;
        if(!side) {
            if(!wobbly)stats->no_consensus++;
            patch_graph_clear(vertex_queue,nzone,core,zone,hit,ray_t,direct,disp,next_disp);
            continue;
        }
        stats->ray_consensus++;if(border_consensus)stats->border_selected++;
        else stats->ambiguous_selected++;
        size_t nhit=0;
        for(size_t q=0;q<nzone;q++) {
            size_t v=vertex_queue[q];if(!zone[v]||!filled[v]||!core[v])continue;
            float t=ray_t[v*2+(side<0)];if(t<0.0f)continue;
            double dir[3];if(!patch_inward_normal(o,verts,normals,v,dir))continue;
            for(int k=0;k<3;k++)direct[v*3+(size_t)k]=
                (float)((double)side*dir[k]*(double)t);
            if(o->preserve_axial)direct[v*3]=0.0f;hit[v]=1;nhit++;
        }
        if(!nhit) {
            stats->no_consensus++;
            patch_graph_clear(vertex_queue,nzone,core,zone,hit,ray_t,direct,disp,next_disp);
            continue;
        }
        double tw=o->patch_target_weight,omega=0.85;
        for(int sweep=0;sweep<o->patch_relax_sweeps;sweep++) {
            for(size_t q=0;q<nzone;q++) {
                size_t v=vertex_queue[q];
                if(!filled[v]){for(int k=0;k<3;k++)next_disp[v*3+(size_t)k]=0.0f;continue;}
                double rhs[3]={0,0,0};int deg=qadj_degree(adj,v);
                for(int e=0;e<deg;e++) {
                    size_t j=qadj_neighbor(adj,v,e);
                    if(j<nv&&zone[j]&&filled[j])for(int k=0;k<3;k++)
                        rhs[k]+=(double)disp[j*3+(size_t)k];
                }
                double den=(double)deg+(hit[v]?tw:0.0);
                if(hit[v])for(int k=0;k<3;k++)rhs[k]+=tw*(double)direct[v*3+(size_t)k];
                for(int k=0;k<3;k++) {
                    size_t z=v*3+(size_t)k;
                    next_disp[z]=(k==0&&o->preserve_axial)||den<=0.0?0.0f:
                        (float)((1.0-omega)*(double)disp[z]+omega*rhs[k]/den);
                }
            }
            for(size_t q=0;q<nzone;q++) {
                size_t v=vertex_queue[q];if(!filled[v])continue;
                for(int k=0;k<3;k++)disp[v*3+(size_t)k]=next_disp[v*3+(size_t)k];
            }
        }
        size_t changed=0;
        for(size_t q=0;q<nzone;q++) {
            size_t v=vertex_queue[q];if(!filled[v])continue;
            double d2=0.0;for(int k=0;k<3;k++)d2+=(double)disp[v*3+(size_t)k]*disp[v*3+(size_t)k];
            if(d2<1e-8)continue;
            for(int k=0;k<3;k++)proposal[v*3+(size_t)k]=
                verts[v*3+(size_t)k]+disp[v*3+(size_t)k];
            if(o->preserve_axial)proposal[v*3]=verts[v*3];
            moved[v]=claimed[v]=1;moved_patch[v]=pr->label;changed++;
        }
        if(changed){stats->proposed++;stats->vertices+=changed;}
        patch_graph_clear(vertex_queue,nzone,core,zone,hit,ray_t,direct,disp,next_disp);
    }
    free(face_dark);free(face_label);free(trusted_label);free(face_queue);free(vertex_queue);
    free(voff);free(vface);free(core);free(zone);free(claimed);free(hit);free(ray_t);
    free(direct);free(disp);free(next_disp);free(head);free(next_face);free(region);
    return 0;
oom:
    free(face_dark);free(face_label);free(trusted_label);free(face_queue);free(vertex_queue);
    free(voff);free(vface);free(core);free(zone);free(claimed);free(hit);free(ray_t);
    free(direct);free(disp);free(next_disp);free(head);free(next_face);free(region);
    return -1;
}

int QuadStripSnap_run(Arena_T arena, float *verts, size_t nv,
                      const int32_t *faces, size_t nf, const uint8_t *filled,
                      const float *uv, int grid_h, int grid_w,
                      const char *raw_dir, long chunk,
                      const QuadStripSnapOpts *opts,
                      int *out_snapped, int *out_nopap,
                      int *out_remesh_advised,
                      QuadStripSnapStats *out_stats) {
    if (out_stats) memset(out_stats, 0, sizeof *out_stats);
    if (verts == NULL || filled == NULL || uv == NULL || nv == 0 || raw_dir == NULL)
        return QUAD_STRIP_SNAP_ERROR;
    QuadStripSnapOpts o = *opts;
    if (o.rounds < 1) o.rounds = 1;
    if (o.solve_sweeps < 1) o.solve_sweeps = 1;
    if (o.solve_omega <= 0.0) o.solve_omega = 0.5;
    if (o.solve_omega > 1.0) o.solve_omega = 1.0;
    if (o.multigrid_levels < 1) o.multigrid_levels = 1;
    if (o.multigrid_levels > QMG_MAX_LEVELS) o.multigrid_levels = QMG_MAX_LEVELS;
    if (o.multigrid_cycles < 0) o.multigrid_cycles = 0;
    if (o.multigrid_patch < 0) o.multigrid_patch = 0;
    if (o.guided_reach < o.reach) o.guided_reach = o.reach;
    if (o.guided_reach > 64.0) o.guided_reach = 64.0;
    if (o.winding_tube < 0.0) o.winding_tube = 0.0;
    if (o.patch_buffer < 0) o.patch_buffer = 0;
    if (o.patch_buffer > 8) o.patch_buffer = 8;
    if (o.patch_min_cells < 1) o.patch_min_cells = 1;
    if (o.patch_ray_reach < 2.0) o.patch_ray_reach = 2.0;
    if (o.patch_ray_reach > 64.0) o.patch_ray_reach = 64.0;
    if (o.patch_ray_margin < 0.0) o.patch_ray_margin = 0.0;
    if (o.patch_relax_sweeps < 1) o.patch_relax_sweeps = 1;
    if (o.patch_relax_sweeps > 512) o.patch_relax_sweeps = 512;
    if (o.patch_target_weight <= 0.0) o.patch_target_weight = 16.0;
    if (o.patch_fitted_dark_slack < 0) o.patch_fitted_dark_slack = 0;
    if (o.bake_window_high <= o.bake_window_low)
        { o.bake_window_low = 47.0; o.bake_window_high = 193.0; }
    if (o.bake_dark_u8 < 0) o.bake_dark_u8 = 0;
    if (o.bake_dark_u8 > 255) o.bake_dark_u8 = 255;
    if (o.quilt_smooth < 0.0) o.quilt_smooth = 0.0;
    if (o.quilt_seam_smooth < 0.0) o.quilt_seam_smooth = 0.0;
    if (o.quilt_anchor < 0.0) o.quilt_anchor = 0.0;
    if (o.quilt_huber < 0.0) o.quilt_huber = 0.0;
    if (o.quilt_conf_min < 0.0) o.quilt_conf_min = 0.0;
    if (o.quilt_conf_min > 1.0) o.quilt_conf_min = 1.0;
    if (o.quilt_sweeps < 1) o.quilt_sweeps = 1;

    CubeTable ct;
    double table_reach = o.guided_reach > o.reach ? o.guided_reach : o.reach;
    if(o.patch_recover&&o.patch_ray_reach>table_reach)table_reach=o.patch_ray_reach;
    if (cubetable_init(&ct, arena, raw_dir, chunk, verts, nv, table_reach + 2.0) != 0)
        return QUAD_STRIP_SNAP_ERROR;
    int nload = cubetable_prewarm_all(&ct);
    size_t expected = cubetable_expected_chunks(&ct);
    if (out_stats) {
        out_stats->raw_chunks_expected = expected;
        out_stats->raw_chunks_loaded = (size_t)ct.n_loaded;
        out_stats->raw_chunks_missing = (size_t)ct.n_missing;
        out_stats->raw_chunks_outside_volume = ct.n_outside;
    }
    if (o.verbose)
        fprintf(stderr, "[relax] RAW coverage %d/%zu in-volume chunks "
                "(%d missing, %zu out-of-volume bbox chunks excluded)\n",
                nload, expected, ct.n_missing, ct.n_outside);
    if (o.require_complete_raw && !cubetable_is_complete(&ct)) {
        fprintf(stderr,
                "ERROR: incomplete RAW source for snap bbox: %d/%zu chunks loaded, "
                "%d missing, %zu out-of-volume bbox chunks excluded (source=%s). "
                "Refusing to relax or emit a production "
                "mesh. Use the authoritative full RAW source; "
                "--allow-incomplete-raw is diagnostics-only.\n",
                nload, expected, ct.n_missing, ct.n_outside, raw_dir);
        return QUAD_STRIP_SNAP_RAW_INCOMPLETE;
    }

    QuadAdj adj;
    memset(&adj, 0, sizeof adj);
    if (grid_h > 0 && grid_w > 0 &&
        (size_t)grid_h <= SIZE_MAX / (size_t)grid_w &&
        (size_t)grid_h * (size_t)grid_w == nv) {
        adj.structured = 1; adj.H = grid_h; adj.W = grid_w;
    } else {
        CSR_T csr = CSR_from_faces(arena, faces, nf, nv);
        adj.off = CSR_offset(csr); adj.tgt = CSR_target(csr);
    }
    if (o.verbose)
        fprintf(stderr, "[relax] adjacency=%s%s\n",
                adj.structured ? "structured-4" : "face-CSR",
                adj.structured ? ", red/black smoother" : ", Jacobi fallback");

    float   *orig   = (float *)malloc(nv * 3 * sizeof *orig);    /* init (fit / harmonic) */
    float   *target = (float *)malloc(nv * 3 * sizeof *target);  /* quilted offset, then absolute target */
    uint8_t *hast   = (uint8_t *)malloc(nv);
    uint8_t *frozen = (uint8_t *)malloc(nv);
    uint8_t *anchored = (uint8_t *)malloc(nv);                  /* recovered recto support */
    uint8_t *bake_dark_seed = (uint8_t *)malloc(nv);             /* dark emitted-triangle corners */
    uint8_t *cand   = (uint8_t *)malloc(nv);                    /* raw profile has a ridge datum */
    float   *cand_t = (float *)malloc(nv * sizeof *cand_t);     /* raw signed ridge depth */
    float   *cand_w = (float *)malloc(nv * sizeof *cand_w);     /* robust raw datum weight */
    float   *guide  = (float *)malloc(nv * 3 * sizeof *guide);  /* oriented/profile search direction */
    float   *conf   = (float *)malloc(nv * sizeof *conf);       /* diffused target confidence */
    float   *face_offset = adj.structured
        ? (float *)malloc(nv * 3 * sizeof *face_offset) : NULL;
    float   *face_weight = adj.structured
        ? (float *)malloc(nv * sizeof *face_weight) : NULL;
    float   *winding_limit = NULL;
    float   *intersection_before = o.prevent_intersection_growth
        ? (float *)malloc(nv*3*sizeof *intersection_before) : NULL;
    float   *intersection_proposed = o.prevent_intersection_growth
        ? (float *)malloc(nv*3*sizeof *intersection_proposed) : NULL;
    size_t  *intersection_degree_current = o.prevent_intersection_growth
        ? (size_t *)calloc(nf,sizeof *intersection_degree_current) : NULL;
    size_t  *intersection_degree_candidate = o.prevent_intersection_growth
        ? (size_t *)calloc(nf,sizeof *intersection_degree_candidate) : NULL;
    double  *R      = (double *)malloc(nv * 9 * sizeof *R);      /* per-cell rotation */
    float   *nxt    = adj.structured ? NULL : (float *)malloc(nv * 3 * sizeof *nxt);
    float   *scalar_nxt = adj.structured ? NULL : (float *)malloc(nv * 4 * sizeof *scalar_nxt);
    if (!orig || !target || !hast || !frozen || !anchored || !bake_dark_seed ||
        !cand || !cand_t || !cand_w || !guide || !conf ||
        !R || (adj.structured&&(!face_offset||!face_weight)) ||
        (o.prevent_intersection_growth&&(!intersection_before||!intersection_proposed||
                                         !intersection_degree_current||
                                         !intersection_degree_candidate)) ||
        (!adj.structured && (!nxt || !scalar_nxt))) {
        free(orig); free(target); free(hast); free(frozen); free(anchored); free(bake_dark_seed);
        free(cand); free(cand_t); free(cand_w);
        free(guide); free(conf); free(face_offset); free(face_weight);
        free(intersection_before);free(intersection_proposed);
        free(intersection_degree_current);free(intersection_degree_candidate);
        free(R); free(nxt); free(scalar_nxt); return -1;
    }
    memcpy(orig, verts, nv * 3 * sizeof *orig);
    memset(hast, 0, nv);
    memset(anchored, 0, nv);
    memset(bake_dark_seed, 0, nv);
    size_t winding_supported=0;double winding_pct[3]={0,0,0};
    if(adj.structured&&snap_oriented(&o)) {
        winding_limit=adaptive_winding_limits(
            orig,adj.H,adj.W,o.axis_y,o.axis_x,o.winding_tube,o.guided_reach,
            &winding_supported,winding_pct);
        if(!winding_limit) {
            free(orig);free(target);free(hast);free(frozen);free(anchored);
            free(bake_dark_seed);free(cand);free(cand_t);free(cand_w);
            free(guide);free(conf);free(face_offset);free(face_weight);
            free(intersection_before);free(intersection_proposed);
            free(intersection_degree_current);free(intersection_degree_candidate);
            free(R);free(nxt);free(scalar_nxt);return -1;
        }
        if(o.verbose)fprintf(stderr,
            "[relax] adaptive winding corridor: local +/-2pi support %zu/%zu, "
            "half-width p05/p50/p95 %.3f/%.3f/%.3f vox (anchored to initializer)\n",
            winding_supported,nv,winding_pct[0],winding_pct[1],winding_pct[2]);
    }

    /* Calibrate the ridge eligibility floor to this scroll.  Absolute 130 was
     * badly wrong on PHerc0139 (direct fitted median is below it) and therefore
     * unpinned most known-good proxy geometry.  Keep a user override, but derive
     * the default in the same window-relative spirit as surface_snap. */
    size_t hfit[256] = {0}, hfill[256] = {0};
    size_t nhfit = 0, nhfill = 0;
    for (size_t v = 0; v < nv; v++) {
        double I = sample_trilinear(&ct, verts[v*3+0], verts[v*3+1], verts[v*3+2]);
        if (I < 0.0) continue;
        int bin = (int)lround(I); if (bin < 0) bin = 0; if (bin > 255) bin = 255;
        if (filled[v]) { hfill[bin]++; nhfill++; }
        else { hfit[bin]++; nhfit++; }
    }
    int fit_p1 = hist_percentile(hfit, nhfit, 1.0);
    int fit_p25 = hist_percentile(hfit, nhfit, 25.0);
    int fit_p50 = hist_percentile(hfit, nhfit, 50.0);
    int fit_p75 = hist_percentile(hfit, nhfit, 75.0);
    int fit_p99 = hist_percentile(hfit, nhfit, 99.0);
    double ridge_min = o.bright_min;
    if (ridge_min <= 0.0) ridge_min = (double)fit_p1 + o.band_frac * (double)(fit_p99 - fit_p1);
    if (ridge_min < 0.0) ridge_min = 0.0;
    if (ridge_min > 255.0) ridge_min = 255.0;
    double dark_max = -1.0;
    if (o.dark_probe_frac > 0.0) {
        dark_max = (double)fit_p1 + o.dark_probe_frac * (double)(fit_p99 - fit_p1);
        if (dark_max < 0.0) dark_max = 0.0;
        if (dark_max > ridge_min) dark_max = ridge_min;
    }
    if (o.verbose) {
        fprintf(stderr, "[relax] direct fitted CT p1/p25/p50/p75/p99=%d/%d/%d/%d/%d "
                        "(sampled %zu), ridge_min=%.1f%s, dark-probe<=%.1f%s\n",
                fit_p1, fit_p25, fit_p50, fit_p75, fit_p99, nhfit, ridge_min,
                o.bright_min <= 0.0 ? " auto" : " explicit", dark_max,
                dark_max >= 0.0 ? " normal-max" : " disabled");
        if (nhfill > 0)
            fprintf(stderr, "[relax] direct fill CT p1/p25/p50/p75/p99=%d/%d/%d/%d/%d "
                            "(sampled %zu)\n",
                    hist_percentile(hfill, nhfill, 1.0), hist_percentile(hfill, nhfill, 25.0),
                    hist_percentile(hfill, nhfill, 50.0), hist_percentile(hfill, nhfill, 75.0),
                    hist_percentile(hfill, nhfill, 99.0), nhfill);
    }

    /* `filled == 0` is certified upstream geometry, not a brightness guess.
     * Freeze it unconditionally.  The previous implementation reinterpreted a
     * dark triangle or an off-ridge profile as permission to move a certified
     * vertex; that converted clean incoming papyrus tracks into CT crack-wall
     * cuts.  Only bounded-hole loft vertices may move here.  On later passes a
     * filled cell already on a credible recto transition becomes a stronger
     * zero-depth soft anchor while retaining sheet coupling. */
    float *init_nrm = MeshNormals_compute(verts, nv, faces, nf);
    if (init_nrm == NULL) {
        free(orig); free(target); free(hast); free(frozen); free(anchored); free(bake_dark_seed);
        free(cand); free(cand_t); free(cand_w);
        free(guide); free(conf); free(face_offset); free(face_weight); free(winding_limit);
        free(intersection_before); free(intersection_proposed);
        free(intersection_degree_current); free(intersection_degree_candidate);
        free(R); free(nxt); free(scalar_nxt); return -1;
    }
    IntersectionCleanupStats intersection_current={0};
    IntersectionCleanupStats intersection_input={0};
    int intersection_guard_audits=0,intersection_guard_backtracks=0;
    int intersection_guard_round_rollbacks=0;
    int intersection_guard_orientation_rejects=0;
    size_t intersection_guard_locally_rolled_back_vertices=0;
    if(o.prevent_intersection_growth) {
        if(snap_exact_intersections(verts,nv,faces,nf,
                                    intersection_degree_current,
                                    &intersection_current)!=0) {
            fprintf(stderr,"ERROR: exact intersection audit failed before snap; "
                           "refusing unguarded geometry motion\n");
            free(init_nrm);
            free(orig); free(target); free(hast); free(frozen); free(anchored); free(bake_dark_seed);
            free(cand); free(cand_t); free(cand_w);
            free(guide); free(conf); free(face_offset); free(face_weight); free(winding_limit);
            free(intersection_before); free(intersection_proposed);
            free(intersection_degree_current); free(intersection_degree_candidate);
            free(R); free(nxt); free(scalar_nxt); return -1;
        }
        intersection_input=intersection_current;
        intersection_guard_audits=1;
        if(o.verbose)fprintf(stderr,
            "[relax] exact round guard input: conflicts %zu "
            "(overlap %zu, stab %zu, fold %zu)\n",
            intersection_current.conflicts,intersection_current.overlap_pairs,
            intersection_current.stab_pairs,intersection_current.fold_pairs);
    }
    RemeshAudit input_ra={0};
    SnapBakeCounts input_bake={0};
    if(adj.structured) {
        input_bake=structured_bake_scan(&ct,verts,init_nrm,filled,uv,
            adj.H,adj.W,o.bake_window_low,o.bake_window_high,o.bake_dark_u8,
            bake_dark_seed,o.cell_lit_accum,NULL);
    } else {
        input_bake=sparse_bake_scan(&ct,verts,init_nrm,faces,nf,filled,nv,
            o.bake_window_low,o.bake_window_high,o.bake_dark_u8,bake_dark_seed,NULL);
    }
    if(out_stats){
        input_ra=remesh_audit(verts,verts,init_nrm,uv,faces,nf,nv,grid_h,grid_w,
            o.remesh_motion,o.remesh_motion_max,o.remesh_sander,o.remesh_over4);
    }
    size_t nmov = 0, nfill_mov = 0, nfill_recto = 0, nfill_badnormal = 0;
    size_t nfit_ridge = 0, nfit_repair = 0, nfit_dark_target = 0;
    size_t nfit_soft = 0;
    size_t nfit_dark_forced = 0;
    size_t nfit_ambiguous = 0, nfit_badnormal = 0;
    size_t nfit_dark = 0, nfill_dark = 0, nfill_dark_target = 0;
    int winding_rejected_total = 0;
    ptrdiff_t vv;
    /* Per-vertex classification: each iteration writes only frozen[v] and
     * anchored[v]; CT reads are thread-safe after the prewarm; every counter
     * reduces commutatively. */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,256) \
    reduction(+:nmov,nfill_mov,nfill_recto,nfill_badnormal,nfit_ridge,nfit_soft, \
                nfit_repair,nfit_dark_target,nfit_dark_forced,nfit_ambiguous, \
                nfit_badnormal,nfit_dark,nfill_dark,nfill_dark_target, \
                winding_rejected_total)
#endif
    for (vv = 0; vv < (ptrdiff_t)nv; vv++) {
        size_t v = (size_t)vv;
        int fitted_soft=!filled[v]&&o.fitted_source_anchor>0.0;
        if (!filled[v]&&!fitted_soft) {
            double I = sample_trilinear(&ct, verts[v*3+0], verts[v*3+1], verts[v*3+2]);
            if (dark_max >= 0.0 && I >= 0.0 && I <= dark_max) nfit_dark++;
            frozen[v] = 1;
            nfit_ridge++; /* counter now means certified/frozen fitted geometry */
            continue;
        }
        double n[3] = { init_nrm[v*3+0], init_nrm[v*3+1], init_nrm[v*3+2] };
        double nn = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        if (nn < 0.5) {
            double I = sample_trilinear(&ct, verts[v*3+0], verts[v*3+1], verts[v*3+2]);
            if (dark_max >= 0.0 && I >= 0.0 && I <= dark_max) {
                if (filled[v]) nfill_dark++; else nfit_dark++;
            }
            if (filled[v]) {
                frozen[v] = 0; nmov++; nfill_mov++; nfill_badnormal++;
            } else if (fitted_soft) {
                frozen[v] = 0; anchored[v] = 2; nmov++; nfit_soft++;
                nfit_badnormal++;
            } else if (bake_dark_seed[v]) {
                frozen[v] = 0; nmov++; nfit_repair++; nfit_dark_forced++;
                nfit_badnormal++;
            } else {
                frozen[v] = 1; nfit_badnormal++;
            }
            continue;
        }
        double p[3] = { verts[v*3+0], verts[v*3+1], verts[v*3+2] };
        double reference[3] = { orig[v*3+0], orig[v*3+1], orig[v*3+2] };
        double search_dir[3]; int rsign = 0, wr_local = 0;
        RidgeProfile rp = guided_ridge_profile(&ct,&o,p,n,ridge_min,
            o.min_gain,o.ridge_lock,dark_max,bake_dark_seed[v] != 0,
            reference,winding_limit?(double)winding_limit[v]:o.winding_tube,
            search_dir,&rsign,
            &wr_local);
        winding_rejected_total += wr_local;
        int is_dark = dark_max >= 0.0 && rp.visible_current >= 0.0 &&
                      rp.visible_current <= dark_max;
        int dark_target = o.snap_side == QUAD_SNAP_SIDE_MIDLINE
            ? dark_band_center_target(&rp, rsign, dark_max, ridge_min,
                                      o.min_gain, o.ridge_lock)
            : dark_inward_recto_target(&rp, rsign, dark_max, ridge_min,
                                       o.min_gain, o.ridge_lock);
        if (is_dark) {
            if (filled[v]) nfill_dark++; else nfit_dark++;
        }
        if(fitted_soft) {
            /* Upstream claims are strong observations, not immutable geometry.
             * Keep every one in the iterated sheet solve with a high source
             * attachment; CT evidence may still make a small coherent move. */
            frozen[v]=0;anchored[v]=2;nmov++;nfit_soft++;
            if(bake_dark_seed[v]) {
                nfit_repair++;nfit_dark_forced++;
                if(dark_target)nfit_dark_target++;
            } else if(rp.on_ridge&&!rp.used_guided)nfit_ridge++;
            else if(rp.has_target) {
                nfit_repair++;if(dark_target)nfit_dark_target++;
            } else nfit_ambiguous++;
            continue;
        }
        if (filled[v]) {
            if (dark_target) nfill_dark_target++;
            if (!bake_dark_seed[v] && (rp.used_recto || rp.used_midline) &&
                rp.on_ridge) {
                anchored[v] = 1; nfill_recto++;
            }
            frozen[v] = 0; nmov++; nfill_mov++;
        }
        else if (bake_dark_seed[v]) {
            /* Vertex-bright / triangle-dark is the signature of a ridge-family
             * seam.  Free all three triangle corners together even when their
             * individual profiles happen to peak at zero depth. */
            frozen[v] = 0; nfit_repair++; nfit_dark_forced++; nmov++;
            if (dark_target) nfit_dark_target++;
        }
        else if (rp.on_ridge && !rp.used_guided) { frozen[v] = 1; nfit_ridge++; }
        else if (rp.has_target) {
            frozen[v] = 0; nfit_repair++; nmov++;
            if (dark_target) nfit_dark_target++;
        }
        else { frozen[v] = 1; nfit_ambiguous++; }
    }
    if (o.verbose)
        fprintf(stderr, "[relax] movable %zu/%zu = bounded fill %zu + fitted repair %zu "
                        "(dark-triangle forced %zu); "
                        "soft-anchored fill: recto %zu, bad-normal %zu; "
                        "fitted source soft-anchor %zu; fitted on-ridge/frozen: "
                        "%zu, ambiguous %zu, bad-normal %zu%s; "
                        "bake-dark fitted/fill %zu/%zu, inward targets %zu/%zu\n",
                nmov, nv, nfill_mov, nfit_repair, nfit_dark_forced,
                nfill_recto, nfill_badnormal,nfit_soft,
                nfit_ridge, nfit_ambiguous, nfit_badnormal,
                o.filled_target_anchor > 0.0 ? " [active-set]" : "",
                nfit_dark, nfill_dark, nfit_dark_target, nfill_dark_target);

    int snapped = 0, nopap = 0;
    int final_ncand = 0, final_ncand_recto = 0, final_ncand_near = 0;
    int final_band_fallback = 0, final_band_abstain = 0;
    int final_raw_no_candidate = 0, final_raw_applied = 0, final_quilt_only = 0;
    int final_guided = 0, final_wide = 0, final_local = 0;
    int final_dark_seed_quilt_only = 0;
    int final_candidate_downweighted = 0;
    double final_candidate_weight_mean = 0.0;
    double final_eraw = 0.0, final_equilt = 0.0;
    size_t final_neraw = 0, final_nequilt = 0;
    QuiltEdgeAccum final_raw_regions = {{0}, {0}};
    QuiltEdgeAccum final_solved_regions = {{0}, {0}};
    QuiltEdgeAccum final_applied_regions = {{0}, {0}};
    QMGRunStats final_quilt_mg = {0}, final_position_mg = {0};
    SnapFaceTargetCounts final_face_counts={0};
    PatchRecoverRunStats patch_stats={0};
    for (int round = 0; round < o.rounds; round++) {
        if(o.prevent_intersection_growth)
            memcpy(intersection_before,verts,nv*3*sizeof *intersection_before);
        float *nrm = MeshNormals_compute(verts, nv, faces, nf);
        if (nrm == NULL) {
            free(init_nrm);
            free(orig); free(target); free(hast); free(frozen); free(anchored); free(bake_dark_seed);
            free(cand); free(cand_t); free(cand_w);
            free(guide); free(conf); free(face_offset); free(face_weight); free(winding_limit);
            free(intersection_before); free(intersection_proposed);
            free(intersection_degree_current); free(intersection_degree_candidate);
            free(R); free(nxt); free(scalar_nxt); return -1;
        }

        if(adj.structured)
            (void)structured_bake_scan(&ct,verts,nrm,filled,uv,adj.H,adj.W,
                o.bake_window_low,o.bake_window_high,o.bake_dark_u8,
                bake_dark_seed,o.cell_lit_accum,NULL);

        SnapFaceTargetCounts face_counts={0};
        if(adj.structured&&winding_limit)
            face_counts=structured_dark_face_targets(
                &ct,&o,verts,orig,nrm,filled,frozen,uv,winding_limit,
                adj.H,adj.W,ridge_min,face_offset,face_weight,
                &winding_rejected_total,o.cell_candidate_accum);

        /* (1) local closest-rotation of each cell's current edges to the flat-UV
         * rest -- for ALL cells; the frozen fitted rim supplies the real wrap
         * curvature the quasi-developable term continues.  Each cell writes
         * only its own frame from read-only neighbours. */
        {
            ptrdiff_t ii;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (ii = 0; ii < (ptrdiff_t)nv; ii++) {
                size_t i = (size_t)ii;
                double S[9] = {0,0,0,0,0,0,0,0,0};
                double ui = (double)uv[i*2+0], vi = (double)uv[i*2+1];
                int deg = qadj_degree(&adj, i);
                for (int e = 0; e < deg; e++) {
                    size_t j = qadj_neighbor(&adj, i, e);
                    double ec[3] = { verts[i*3+0]-verts[j*3+0], verts[i*3+1]-verts[j*3+1], verts[i*3+2]-verts[j*3+2] };
                    double er[3] = { ui-(double)uv[j*2+0], vi-(double)uv[j*2+1], 0.0 };
                    for (int a = 0; a < 3; a++) for (int c = 0; c < 3; c++) S[a*3+c] += ec[a]*er[c];
                }
                closest_rotation(S, &R[i*9]);
            }
        }

        /* (2) Detect a raw signed normal offset, convert it to a physical 3-D
         * offset, then QUILT that vector field before constructing targets.
         * Independent bright
         * vertices are not enough: adjacent vertices on different ridges make
         * the face interior pass through dark resin (the observed vertex-bright /
         * raster-dark failure).  A vector energy sees a target-family switch
         * even when independently oriented scalar depths happen to agree.
         * Confidence diffuses with the same structured
         * smoother, so no-ridge cells may follow coherent nearby evidence while
         * wholly unsupported regions remain governed by orig/smooth/dev only. */
        int ncand = 0, ncand_near = 0, ncand_recto = 0;
        int ncand_guided = 0, ncand_wide = 0, ncand_local = 0;
        int ndark_seed_quilt_only = 0;
        int nband_fallback = 0, nband_abstain = 0;
        /* Per-vertex target search is the round's dominant CT cost (tens of
         * trilinear reads per movable vertex); each iteration writes only its
         * own cand/conf/target/guide slots. */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,256) \
    reduction(+:ncand,ncand_near,ncand_recto,ncand_guided,ncand_wide, \
                ncand_local,ndark_seed_quilt_only,nband_fallback, \
                nband_abstain,winding_rejected_total)
#endif
        for (vv = 0; vv < (ptrdiff_t)nv; vv++) {
            size_t v = (size_t)vv;
            cand[v] = 0; cand_t[v] = 0.0f; cand_w[v]=0.0f;
            conf[v] = 0.0f; hast[v] = 0;
            for(int k=0;k<3;k++){target[v*3+(size_t)k]=0.0f;guide[v*3+(size_t)k]=0.0f;}
            if (frozen[v]) continue;
            double n[3] = { nrm[v*3+0], nrm[v*3+1], nrm[v*3+2] };
            double nn = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
            if (nn < 0.5) {
                guide[v*3+0]=guide[v*3+1]=guide[v*3+2]=0.0f;
                continue;
            }
            double p[3] = { verts[v*3+0], verts[v*3+1], verts[v*3+2] };
            double reference[3] = { orig[v*3+0], orig[v*3+1], orig[v*3+2] };
            double search_dir[3]; int rsign = 0, wr_local = 0;
            RidgeProfile rp = guided_ridge_profile(&ct,&o,p,n,ridge_min,
                o.min_gain,o.ridge_lock,dark_max,bake_dark_seed[v] != 0,
                reference,winding_limit?(double)winding_limit[v]:o.winding_tube,
                search_dir,&rsign,
                &wr_local);
            winding_rejected_total += wr_local;
            for(int k=0;k<3;k++)guide[v*3+k]=(float)search_dir[k];
            if (o.snap_side == QUAD_SNAP_SIDE_MIDLINE && !rp.has_target)
                nband_abstain++;
            if (anchored[v]==1) {
                /* A recovered recto site is evidence that zero further normal
                 * motion is correct.  Keep it in both quilting and position
                 * coupling, rather than making a hard active-set boundary. */
                cand[v] = 1; cand_t[v] = 0.0f; conf[v] = 1.0f;
                ncand++; ncand_near++; ncand_recto++;
                continue;
            }
            double t = 0.0;
            int accept = 0;
            int dark_seed = bake_dark_seed[v] != 0;
            int supported_dark_alternative = rp.has_target &&
                (rp.used_guided || fabs(rp.target_t) > o.ridge_lock) &&
                rp.visible_target >= ridge_min;
            if (rp.has_target && (!dark_seed || supported_dark_alternative)) {
                t = rp.target_t; accept = 1;
            }
            else if (!dark_seed && rp.on_ridge && rp.nearest_peak_value >= ridge_min) {
                /* Once at the ridge, min_gain intentionally vanishes.  Retain a
                 * near-zero datum so the next round cannot alternate away/back. */
                t = rp.nearest_peak_t; accept = 1; ncand_near++;
            } else if (dark_seed) ndark_seed_quilt_only++;
            if (accept) {
                cand[v] = 1; cand_t[v] = (float)t; conf[v] = 1.0f; ncand++;
                for(int k=0;k<3;k++)
                    target[v*3+(size_t)k]=(float)(t*(double)guide[v*3+(size_t)k]);
                if (rp.used_recto || rp.used_midline) ncand_recto++;
                if (rp.used_band_fallback) nband_fallback++;
                if (rp.used_guided) ncand_guided++;
                if (rp.used_wide) ncand_wide++;
                if (rp.used_local_contrast) ncand_local++;
            }
        }

        /* A dark face target describes the emitted sample that is actually
         * wrong, so it supersedes a conflicting corner-only profile.  Store it
         * in the same physical-vector representation consumed by the robust
         * quilt. */
        if(adj.structured&&winding_limit)for(size_t v=0;v<nv;v++) {
            if(frozen[v]||face_weight[v]<=0.0f)continue;
            double delta[3];double d2=0.0;
            for(int k=0;k<3;k++) {
                delta[k]=(double)face_offset[v*3+(size_t)k]/(double)face_weight[v];
                if(o.preserve_axial&&k==0)delta[k]=0.0;
                d2+=delta[k]*delta[k];
            }
            double dn=sqrt(d2);if(dn<=1e-6)continue;
            int had=cand[v]!=0;
            if(!had) {
                ncand++;
                if(bake_dark_seed[v]&&ndark_seed_quilt_only>0)
                    ndark_seed_quilt_only--;
            }
            cand[v]=1;cand_t[v]=(float)dn;conf[v]=1.0f;
            for(int k=0;k<3;k++) {
                guide[v*3+(size_t)k]=(float)(delta[k]/dn);
                target[v*3+(size_t)k]=(float)delta[k];
            }
        }

        int candidate_downweighted=0;double candidate_weight_mean=0.0;
        quilt_candidate_weights(cand,cand_t,guide,anchored,&adj,nv,o.quilt_huber,
                                 cand_w,&candidate_downweighted,
                                 &candidate_weight_mean);

        QMGRunStats quilt_mg = {0};
        if (adj.structured) {
            if (qmg_solve_quilt(target, conf, frozen, filled, cand_w, cand_t,guide,
                                adj.H, adj.W, &o, &quilt_mg) != 0) {
                quilt_mg.levels_used = 1;
                if (o.verbose) fprintf(stderr,
                    "[relax] quilt multigrid allocation failed; retaining fine-grid solution\n");
            }
        } else {
            for (int sweep = 0; sweep < o.quilt_sweeps; sweep++) {
                for (size_t v = 0; v < nv; v++) {
                    if (frozen[v]) continue;
                    double wd = (double)cand_w[v];
                    double den = wd + o.quilt_anchor;
                    double nd[3]={wd*(double)cand_t[v]*(double)guide[v*3+0],
                                  wd*(double)cand_t[v]*(double)guide[v*3+1],
                                  wd*(double)cand_t[v]*(double)guide[v*3+2]};
                    double nc = wd;
                    int deg = qadj_degree(&adj, v);
                    for (int e = 0; e < deg; e++) {
                        size_t j = qadj_neighbor(&adj, v, e);
                        double ws = o.quilt_smooth;
                        if (filled[v] != filled[j]) ws *= o.quilt_seam_smooth;
                        den += ws;
                        if (!frozen[j]) {
                            for(int k=0;k<3;k++)nd[k]+=ws*(double)target[j*3+(size_t)k];
                            nc += ws * (double)conf[j];
                        }
                    }
                    for(int k=0;k<3;k++)
                        scalar_nxt[v*4+(size_t)k]=den>1e-12?(float)(nd[k]/den):0.0f;
                    scalar_nxt[v*4+3]=den>1e-12?(float)(nc/den):0.0f;
                }
                for (size_t v = 0; v < nv; v++) if (!frozen[v]) {
                    for(int k=0;k<3;k++)target[v*3+(size_t)k]=scalar_nxt[v*4+(size_t)k];
                    conf[v]=scalar_nxt[v*4+3];
                }
            }
            quilt_mg.levels_used = 1;
        }
        final_quilt_mg = quilt_mg;

        snapped = 0;
        int raw_no_candidate = 0, raw_applied = 0, quilt_only = 0;
        double mean_abs_raw = 0.0, mean_abs_quilt = 0.0;
        for (size_t v = 0; v < nv; v++) {
            if (cand[v]) mean_abs_raw += fabs((double)cand_t[v]);
            if (!frozen[v] && !cand[v]) raw_no_candidate++;
            if (frozen[v] || conf[v] < o.quilt_conf_min) { hast[v] = 0; continue; }
            double q2=0.0;for(int k=0;k<3;k++)q2+=(double)target[v*3+(size_t)k]*target[v*3+(size_t)k];
            hast[v] = 1; snapped++; mean_abs_quilt += sqrt(q2);
            if (cand[v]) raw_applied++; else quilt_only++;
        }
        /* `nopap` historically subtracted post-diffusion confidence from the
         * movable population, so multigrid could make it nearly vanish without
         * discovering one new CT ridge.  Preserve the field for compatibility,
         * but give it the honest raw-no-candidate semantics. */
        nopap = raw_no_candidate;

        /* Explicit quilting-energy diagnostics.  Preserve the transaction
         * gate's historical raw/raw and active/active RMS, then additionally
         * partition raw, complete solved (frozen = zero Dirichlet), and actually
         * applied (inactive = zero) fields by fitted/fill edge class. */
        double eraw = 0.0, equilt = 0.0; size_t neraw = 0, nequilt = 0;
        QuiltEdgeAccum raw_regions = {{0}, {0}};
        QuiltEdgeAccum solved_regions = {{0}, {0}};
        QuiltEdgeAccum applied_regions = {{0}, {0}};
        for (size_t v = 0; v < nv; v++) {
            int deg = qadj_degree(&adj, v);
            for (int e = 0; e < deg; e++) {
                size_t j = qadj_neighbor(&adj, v, e);
                if (j <= v) continue;
                int edge_class = quilt_edge_class(filled[v], filled[j]);
                if (cand[v] && cand[j]) {
                    double d2=0.0;for(int k=0;k<3;k++){
                        double d=(double)cand_t[v]*(double)guide[v*3+(size_t)k]
                                -(double)cand_t[j]*(double)guide[j*3+(size_t)k];
                        d2+=d*d;
                    }
                    double d=sqrt(d2);eraw+=d2;neraw++;
                    quilt_edge_accum(&raw_regions, edge_class, d);
                }
                if (hast[v] && hast[j]) {
                    double d2=0.0;for(int k=0;k<3;k++){
                        double d=(double)target[v*3+(size_t)k]-target[j*3+(size_t)k];d2+=d*d;
                    }
                    equilt+=d2; nequilt++;
                }
                if (!frozen[v] || !frozen[j]) {
                    double d2=0.0;for(int k=0;k<3;k++){
                        double dv=frozen[v]?0.0:(double)target[v*3+(size_t)k];
                        double dj=frozen[j]?0.0:(double)target[j*3+(size_t)k];
                        double d=dv-dj;d2+=d*d;
                    }
                    quilt_edge_accum(&solved_regions,edge_class,sqrt(d2));
                }
                if (hast[v] || hast[j]) {
                    double d2=0.0;for(int k=0;k<3;k++){
                        double dv=hast[v]?(double)target[v*3+(size_t)k]:0.0;
                        double dj=hast[j]?(double)target[j*3+(size_t)k]:0.0;
                        double d=dv-dj;d2+=d*d;
                    }
                    quilt_edge_accum(&applied_regions,edge_class,sqrt(d2));
                }
            }
        }
        final_ncand = ncand;
        final_ncand_recto = ncand_recto;
        final_band_fallback = nband_fallback;
        final_band_abstain = nband_abstain;
        final_ncand_near = ncand_near;
        final_raw_no_candidate = raw_no_candidate;
        final_raw_applied = raw_applied;
        final_quilt_only = quilt_only;
        final_guided = ncand_guided;
        final_wide = ncand_wide;
        final_local = ncand_local;
        final_dark_seed_quilt_only = ndark_seed_quilt_only;
        final_face_counts=face_counts;
        final_eraw = eraw; final_equilt = equilt;
        final_neraw = neraw; final_nequilt = nequilt;
        final_raw_regions = raw_regions;
        final_solved_regions = solved_regions;
        final_applied_regions = applied_regions;
        final_candidate_downweighted=candidate_downweighted;
        final_candidate_weight_mean=candidate_weight_mean;
        for(size_t v=0;v<nv;v++)if(hast[v])for(int k=0;k<3;k++)
            target[v*3+(size_t)k]=(float)((double)verts[v*3+(size_t)k]
                                         +(double)target[v*3+(size_t)k]);
        free(nrm);

        /* (3) Minimize the combined soft position energy over movable cells:
         *   w_orig  (data anchor) + w_target (ridge snap, if found)
         *   + per-neighbour [ w_smooth (continuity) + w_dev (quasi-developable
         *     ARAP: pull toward v_nbr + Rbar*(rest_v - rest_nbr)) ].
         * Structured ribbons get in-place red/black Gauss-Seidel sweeps plus a
         * Galerkin aggregate correction when requested.  Coarse position
         * proposals still pass the emitted-triangle barrier before commit. */
        int barrier_backtracks = 0, barrier_rejects = 0;
        QMGRunStats position_mg = {0};
        if (adj.structured) {
            if (qmg_solve_position(verts,orig,target,uv,R,frozen,anchored,hast,
                                   conf,winding_limit,&adj,&o,&barrier_backtracks,
                                   &barrier_rejects,&position_mg) != 0) {
                position_mg.levels_used=1;
                if(o.verbose)fprintf(stderr,
                    "[relax] position multigrid allocation failed; retaining fine-grid solution\n");
            }
        } else {
            for (int sweep = 0; sweep < o.solve_sweeps; sweep++) {
                for (size_t v = 0; v < nv; v++) {
                    if (frozen[v]) continue;
                    double ui = (double)uv[v*2+0], vi = (double)uv[v*2+1];
                    double num[3] = {0,0,0}, den = 0.0;
                    double wo = o.w_orig * snap_orig_anchor_multiplier(&o,anchored[v]);
                    num[0]+=wo*snap_orig_anchor_coordinate(&o,orig,v,0,anchored[v]);
                    num[1]+=wo*snap_orig_anchor_coordinate(&o,orig,v,1,anchored[v]);
                    num[2]+=wo*snap_orig_anchor_coordinate(&o,orig,v,2,anchored[v]);
                    den += wo;
                    if (hast[v]) {
                        double wt = o.w_target * (double)conf[v];
                        num[0]+=wt*target[v*3+0]; num[1]+=wt*target[v*3+1]; num[2]+=wt*target[v*3+2];
                        den += wt;
                    }
                    double wn = o.w_smooth + o.w_dev;
                    const double *Rv = &R[v*9];
                    int deg = qadj_degree(&adj, v);
                    for (int e = 0; e < deg; e++) {
                        size_t j = qadj_neighbor(&adj, v, e);
                        num[0]+=wn*verts[j*3+0]; num[1]+=wn*verts[j*3+1]; num[2]+=wn*verts[j*3+2]; den += wn;
                        if (o.w_dev != 0.0) {
                            double er[3] = { ui-(double)uv[j*2+0], vi-(double)uv[j*2+1], 0.0 };
                            const double *Rj = &R[j*9];
                            for (int a = 0; a < 3; a++) {
                                double re = 0.5*(Rv[a*3+0]+Rj[a*3+0])*er[0] + 0.5*(Rv[a*3+1]+Rj[a*3+1])*er[1];
                                num[a] += o.w_dev * re;
                            }
                        }
                    }
                    if (den < 1e-9) den = 1.0;
                    double proposal[3]={num[0]/den,num[1]/den,num[2]/den};
                    clamp_position_from_orig(proposal,orig,v,o.guided_reach);
                    clamp_position_to_winding(proposal,orig,v,winding_limit,&o);
                    for(int k=0;k<3;k++)nxt[v*3+(size_t)k]=(float)proposal[k];
                }
                for (size_t v = 0; v < nv; v++)
                    if (!frozen[v]) { verts[v*3+0]=nxt[v*3+0]; verts[v*3+1]=nxt[v*3+1]; verts[v*3+2]=nxt[v*3+2]; }
            }
            position_mg.levels_used=1;
        }
        if(o.preserve_axial||winding_limit)for(size_t v=0;v<nv;v++) {
            double p[3]={verts[v*3],verts[v*3+1],verts[v*3+2]};
            clamp_position_to_winding(p,orig,v,winding_limit,&o);
            for(int k=0;k<3;k++)verts[v*3+(size_t)k]=(float)p[k];
        }
        if(o.prevent_intersection_growth) {
            /* The local emitted-triangle barrier prevents folds inside a UV
             * neighbourhood, but it cannot see a collision with a remote turn
             * of the scroll.  Treat the complete round as a transaction and
             * backtrack its displacement against the exact whole-mesh detector.
             * Counts (rather than pair identities) are compared because a tiny
             * legal motion inside an already-overlapping thick region can change
             * which tessellation pair represents the same physical conflict. */
            memcpy(intersection_proposed,verts,nv*3*sizeof *intersection_proposed);
            IntersectionCleanupStats accepted_stats=intersection_current;
            int accepted=0;
            size_t accepted_zeroed=0;
            double accepted_scale_mean=0.0;
            for(size_t v=0;v<nv;v++)cand_t[v]=1.0f;
            for(int attempt=0;attempt<=12;attempt++) {
                size_t zeroed=0;
                double scale_mean=1.0;
                if(attempt>0) {
                    intersection_guard_backtracks++;
                    size_t marked=0;
                    for(size_t f=0;f<nf;f++)
                        if(intersection_degree_candidate[f]>
                           intersection_degree_current[f])
                            for(int k=0;k<3;k++) {
                                int32_t vi=faces[f*3+(size_t)k];
                                if(vi>=0&&(size_t)vi<nv&&cand_t[(size_t)vi]>0.0f) {
                                    cand_t[(size_t)vi]=0.0f;marked++;
                                }
                            }
                    size_t orientation_marked=0;
                    (void)snap_orientation_violations(verts,orig,nv,faces,nf,
                                                      cand_t,&orientation_marked);
                    marked+=orientation_marked;
                    /* If pair identities churn without exposing a new-degree
                     * face, or local rollback has not converged after six
                     * tries, halve the remaining global motion as a fail-safe.
                     * Usually the degree-growth mask is sufficient and the
                     * conflict-free majority remains at scale one. */
                    if(marked==0||attempt>6)
                        for(size_t v=0;v<nv;v++)cand_t[v]*=0.5f;
                    /* Feather a zeroed face through four UV/mesh rings.  A hard
                     * zero/full boundary creates exactly the fold this guard is
                     * meant to avoid.  The relaxation is monotone: later repair
                     * attempts can only reduce a vertex's retained motion. */
                    for(int ring=0;ring<4;ring++) {
                        for(size_t v=0;v<nv;v++) {
                            float a=cand_t[v];
                            int deg=qadj_degree(&adj,v);
                            for(int e=0;e<deg;e++) {
                                size_t j=qadj_neighbor(&adj,v,e);
                                float b=cand_t[j]+0.25f;
                                if(b<a)a=b;
                            }
                            cand_w[v]=a<1.0f?a:1.0f;
                        }
                        memcpy(cand_t,cand_w,nv*sizeof *cand_t);
                    }
                    double scale_sum=0.0;
                    for(size_t v=0;v<nv;v++) {
                        if(cand_t[v]<=0.0f)zeroed++;
                        double a=(double)cand_t[v];
                        scale_sum+=a;
                        for(int k=0;k<3;k++) {
                            size_t i=v*3+(size_t)k;
                            verts[i]=(float)((double)intersection_before[i]
                                +a*((double)intersection_proposed[i]
                                    -(double)intersection_before[i]));
                        }
                    }
                    scale_mean=nv?scale_sum/(double)nv:0.0;
                    if(zeroed>intersection_guard_locally_rolled_back_vertices)
                        intersection_guard_locally_rolled_back_vertices=zeroed;
                    if(o.preserve_axial||winding_limit)for(size_t v=0;v<nv;v++) {
                        double p[3]={verts[v*3],verts[v*3+1],verts[v*3+2]};
                        clamp_position_to_winding(p,orig,v,winding_limit,&o);
                        for(int k=0;k<3;k++)verts[v*3+(size_t)k]=(float)p[k];
                    }
                }
                size_t orientation_bad=snap_orientation_violations(
                    verts,orig,nv,faces,nf,NULL,NULL);
                IntersectionCleanupStats candidate_stats={0};
                if(snap_exact_intersections(verts,nv,faces,nf,
                                            intersection_degree_candidate,
                                            &candidate_stats)!=0) {
                    memcpy(verts,intersection_before,nv*3*sizeof *verts);
                    fprintf(stderr,"ERROR: exact intersection audit failed during snap "
                                   "round %d; restored pre-round geometry\n",round+1);
                    free(init_nrm);
                    free(orig); free(target); free(hast); free(frozen); free(anchored); free(bake_dark_seed);
                    free(cand); free(cand_t); free(cand_w);
                    free(guide); free(conf); free(face_offset); free(face_weight); free(winding_limit);
                    free(intersection_before); free(intersection_proposed);
                    free(intersection_degree_current); free(intersection_degree_candidate);
                    free(R); free(nxt); free(scalar_nxt); return -1;
                }
                intersection_guard_audits++;
                if(candidate_stats.conflicts<=intersection_current.conflicts&&
                   candidate_stats.fold_pairs<=intersection_current.fold_pairs&&
                   orientation_bad==0) {
                    accepted_stats=candidate_stats;
                    accepted_zeroed=zeroed;
                    accepted_scale_mean=scale_mean;
                    accepted=1;
                    break;
                }
                if(orientation_bad)intersection_guard_orientation_rejects++;
                if(o.verbose)fprintf(stderr,
                    "[relax] round %d exact guard rejects %s attempt %d "
                    "(mean retained motion %.6f, local zero %zu): "
                    "conflicts %zu->%zu, folds %zu->%zu, "
                    "orientation violations %zu\n",
                    round+1,attempt?"repaired":"full",attempt,
                    scale_mean,zeroed,intersection_current.conflicts,
                    candidate_stats.conflicts,intersection_current.fold_pairs,
                    candidate_stats.fold_pairs,orientation_bad);
            }
            if(accepted) {
                if(o.verbose)fprintf(stderr,
                    "[relax] round %d exact guard accepts mean retained motion "
                    "%.6f (local zero %zu): "
                    "conflicts %zu->%zu (overlap %zu, stab %zu, fold %zu)\n",
                    round+1,accepted_scale_mean,accepted_zeroed,
                    intersection_current.conflicts,
                    accepted_stats.conflicts,accepted_stats.overlap_pairs,
                    accepted_stats.stab_pairs,accepted_stats.fold_pairs);
                intersection_current=accepted_stats;
                size_t *degree_swap=intersection_degree_current;
                intersection_degree_current=intersection_degree_candidate;
                intersection_degree_candidate=degree_swap;
            } else {
                memcpy(verts,intersection_before,nv*3*sizeof *verts);
                intersection_guard_round_rollbacks++;
                if(o.verbose)fprintf(stderr,
                    "[relax] round %d exact guard: no safe nonzero step; "
                    "restored pre-round geometry (%zu conflicts, %zu folds)\n",
                    round+1,intersection_current.conflicts,
                    intersection_current.fold_pairs);
            }
        }
        final_position_mg=position_mg;
        if (o.verbose)
            fprintf(stderr, "[relax] round %d/%d: raw-ridge %d (%d recto, %d near-zero), "
                            "quilt-target %d (raw-applied %d, quilt-only %d), raw-no-candidate %d; "
                            "guided/wide/local %d/%d/%d, dark-seed quilt-only %d; "
                            "dark-face samples/candidates/outward/missed %zu/%zu/%zu/%zu; "
                            "|t| raw %.2f -> quilt %.2f, "
                            "edge-rms %.3f -> %.3f, barrier backtracks/rejects %d/%d\n",
                    round+1, o.rounds, ncand, ncand_recto, ncand_near,
                    snapped,raw_applied,quilt_only,raw_no_candidate,
                    ncand_guided,ncand_wide,ncand_local,ndark_seed_quilt_only,
                    face_counts.dark_samples,face_counts.candidates,
                    face_counts.outward_candidates,face_counts.no_candidate,
                    ncand ? mean_abs_raw/(double)ncand : 0.0,
                    snapped ? mean_abs_quilt/(double)snapped : 0.0,
                    neraw ? sqrt(eraw/(double)neraw) : 0.0,
                    nequilt ? sqrt(equilt/(double)nequilt) : 0.0,
                    barrier_backtracks, barrier_rejects);
        if(o.verbose&&(quilt_mg.levels_used>1||position_mg.levels_used>1))
            fprintf(stderr,
                "[relax] multigrid: quilt L%d corr %d residual %.3e->%.3e; "
                "position L%d corr %d residual %.3e->%.3e coarse backtracks/rejects %d/%d, "
                "patches accepted/rejected %d/%d\n",
                quilt_mg.levels_used,quilt_mg.corrections,
                quilt_mg.residual_before,quilt_mg.residual_after,
                position_mg.levels_used,position_mg.corrections,
                position_mg.residual_before,position_mg.residual_after,
                position_mg.backtracks,position_mg.rejects,
                position_mg.patches_accepted,position_mg.patches_rejected);
    }
    /* The patch stage is deliberately after all short-range quilt rounds.  It
     * therefore sees only their residual raster-dark components, and its wider
     * rays cannot influence the ordinary pointwise solve.  Commit all accepted
     * patches together so global agreement is tested against exact geometry as
     * one coherent surface motion.  Exact-conflict repair attenuates and then,
     * only if necessary, rejects WHOLE patch labels rather than isolated
     * vertices: every retained vertex therefore keeps the selected direction
     * and the already-solved smooth collar. */
    if(o.patch_recover&&snap_oriented(&o)) {
        float *patch_before=intersection_before;
        float *patch_proposed=intersection_proposed;
        int own_patch_arrays=0,patch_failed=0,have_patch_before=0;
        if(!patch_before||!patch_proposed) {
            patch_before=(float *)malloc(nv*3*sizeof *patch_before);
            patch_proposed=(float *)malloc(nv*3*sizeof *patch_proposed);
            own_patch_arrays=1;
        }
        uint8_t *patch_moved=(uint8_t *)malloc(nv);
        int32_t *patch_id=(int32_t *)malloc(nv*sizeof *patch_id);
        uint8_t *patch_alive=NULL;
        uint8_t *patch_reduce=NULL;
        float *patch_scale=NULL;
        float *patch_nrm=MeshNormals_compute(verts,nv,faces,nf);
        if(!patch_before||!patch_proposed||!patch_moved||!patch_id||!patch_nrm)patch_failed=1;
        if(!patch_failed) {
            memcpy(patch_before,verts,nv*3*sizeof *patch_before);
            have_patch_before=1;
            int patch_rc=adj.structured
                ?patch_recover_propose(&ct,&o,patch_before,patch_nrm,filled,uv,
                    adj.H,adj.W,ridge_min,patch_proposed,patch_moved,patch_id,
                    &patch_stats)
                :patch_recover_graph_propose(&ct,&o,patch_before,patch_nrm,filled,
                    uv,faces,nf,&adj,nv,ridge_min,patch_proposed,patch_moved,
                    patch_id,&patch_stats);
            if(patch_rc!=0)
                patch_failed=1;
        }
        free(patch_nrm);
        if(!patch_failed&&patch_stats.proposed) {
            patch_alive=(uint8_t *)calloc(patch_stats.regions,1);
            patch_reduce=(uint8_t *)calloc(patch_stats.regions,1);
            patch_scale=(float *)calloc(patch_stats.regions,sizeof *patch_scale);
            if(!patch_alive||!patch_reduce||!patch_scale)patch_failed=1;
            for(size_t v=0;!patch_failed&&v<nv;v++)
                if(patch_moved[v]&&patch_id[v]>=0&&
                   (size_t)patch_id[v]<patch_stats.regions) {
                    patch_alive[(size_t)patch_id[v]]=1;
                    patch_scale[(size_t)patch_id[v]]=1.0f;
                }
        }
        if(!patch_failed&&patch_stats.proposed) {
            /* Judge the patch against the beginning of the whole proximal pass,
             * not merely against the already-relaxed pre-patch surface.  Two
             * individually sub-90-degree rotations can otherwise compose into
             * an absolute inversion that this local transaction misses and the
             * outer component transaction must roll back wholesale. */
            size_t orientation_base=snap_orientation_violations(
                patch_before,orig,nv,faces,nf,NULL,NULL);
            int accepted=0;
            IntersectionCleanupStats accepted_ix=intersection_current;
            SnapBakeCounts accepted_bake={0};
            for(size_t v=0;v<nv;v++) {
                int32_t id=patch_id[v];
                cand_t[v]=(patch_moved[v]&&id>=0&&
                    (size_t)id<patch_stats.regions)?patch_scale[(size_t)id]:0.0f;
            }
            size_t max_attempts=patch_stats.proposed+8u;
            if(max_attempts<8u)max_attempts=8u;if(max_attempts>64u)max_attempts=64u;
            for(size_t attempt=0;attempt<max_attempts;attempt++) {
                if(attempt>0) {
                    size_t reduced=0,killed=0;
                    memset(patch_reduce,0,patch_stats.regions);
                    if(o.prevent_intersection_growth)for(size_t f=0;f<nf;f++)
                        if(intersection_degree_candidate[f]>
                           intersection_degree_current[f])
                            for(int k=0;k<3;k++) {
                                int32_t vi=faces[f*3+(size_t)k];
                                if(vi<0||(size_t)vi>=nv||!patch_moved[(size_t)vi])continue;
                                int32_t id=patch_id[(size_t)vi];
                                if(id>=0&&(size_t)id<patch_stats.regions&&patch_alive[(size_t)id])
                                    patch_reduce[(size_t)id]=1;
                            }
                    memcpy(cand_w,cand_t,nv*sizeof *cand_w);
                    size_t orientation_marked=0;
                    (void)snap_orientation_violations(verts,orig,nv,faces,nf,
                                                       cand_w,&orientation_marked);
                    if(orientation_marked)for(size_t v=0;v<nv;v++)
                        if(patch_moved[v]&&cand_t[v]>0.0f&&cand_w[v]<=0.0f) {
                            int32_t id=patch_id[v];
                            if(id>=0&&(size_t)id<patch_stats.regions&&patch_alive[(size_t)id])
                                patch_reduce[(size_t)id]=1;
                        }
                    /* Preserve the patch-wide direction and smooth collar while
                     * attenuating a geometrically guilty patch.  Four failed
                     * scales (1, 1/2, 1/4, 1/8) are enough evidence to discard
                     * it.  If exact pair churn or the bake gate exposes no
                     * unique culprit, attenuate every survivor together. */
                    int have_mark=0;
                    for(size_t k=0;k<patch_stats.regions;k++)
                        if(patch_alive[k]&&patch_reduce[k]){have_mark=1;break;}
                    if(!have_mark)for(size_t k=0;k<patch_stats.regions;k++)
                        if(patch_alive[k])patch_reduce[k]=1;
                    for(size_t k=0;k<patch_stats.regions;k++)if(patch_alive[k]&&patch_reduce[k]) {
                        patch_scale[k]*=0.5f;reduced++;
                        if(patch_scale[k]<0.1249f) {
                            patch_scale[k]=0.0f;patch_alive[k]=0;killed++;
                        }
                    }
                    for(size_t v=0;v<nv;v++) {
                        int32_t id=patch_id[v];
                        cand_t[v]=(patch_moved[v]&&id>=0&&
                            (size_t)id<patch_stats.regions&&patch_alive[(size_t)id])
                            ?patch_scale[(size_t)id]:0.0f;
                    }
                    (void)reduced;(void)killed;
                }
                for(size_t v=0;v<nv;v++)for(int k=0;k<3;k++) {
                    size_t q=v*3+(size_t)k;
                    verts[q]=patch_moved[v]?(float)((double)patch_before[q]+(double)cand_t[v]*
                        ((double)patch_proposed[q]-(double)patch_before[q])):patch_before[q];
                }
                float *candidate_nrm=MeshNormals_compute(verts,nv,faces,nf);
                if(!candidate_nrm){patch_failed=1;break;}
                SnapBakeCounts candidate_bake=adj.structured
                    ?structured_bake_counts(&ct,verts,candidate_nrm,filled,uv,
                        adj.H,adj.W,o.bake_window_low,o.bake_window_high,
                        o.bake_dark_u8)
                    :sparse_bake_scan(&ct,verts,candidate_nrm,faces,nf,filled,nv,
                        o.bake_window_low,o.bake_window_high,o.bake_dark_u8,
                        NULL,NULL);
                free(candidate_nrm);
                size_t orientation_bad=snap_orientation_violations(
                    verts,orig,nv,faces,nf,NULL,NULL);
                /* Fitted threshold samples can toggle solely from normal
                 * interpolation next to a moved unknown face.  The explicit
                 * small allowance prevents that boundary raster noise from
                 * globally attenuating otherwise geometry-safe patches. */
                int bake_better=candidate_bake.filled_dark<patch_stats.dark_before&&
                    candidate_bake.fitted_dark<=patch_stats.fitted_dark_before+
                        (size_t)o.patch_fitted_dark_slack&&
                    candidate_bake.missing<=patch_stats.missing_before;
                int exact_ok=1;IntersectionCleanupStats candidate_ix={0};
                if(o.prevent_intersection_growth) {
                    if(snap_exact_intersections(verts,nv,faces,nf,
                        intersection_degree_candidate,&candidate_ix)!=0) {
                        patch_failed=1;break;
                    }
                    intersection_guard_audits++;
                    exact_ok=candidate_ix.conflicts<=intersection_current.conflicts&&
                             candidate_ix.fold_pairs<=intersection_current.fold_pairs;
                }
                if(bake_better&&orientation_bad<=orientation_base&&exact_ok) {
                    accepted=1;accepted_ix=candidate_ix;accepted_bake=candidate_bake;
                    break;
                }
                intersection_guard_backtracks++;
                if(orientation_bad>orientation_base)intersection_guard_orientation_rejects++;
                size_t active_patches=0,active_vertices=0;double scale_sum=0.0;
                for(size_t k=0;k<patch_stats.regions;k++)active_patches+=patch_alive[k]!=0;
                for(size_t v=0;v<nv;v++)if(patch_moved[v]&&cand_t[v]>0.0f) {
                    active_vertices++;scale_sum+=(double)cand_t[v];
                }
                size_t rolled=patch_stats.vertices-active_vertices;
                if(rolled>intersection_guard_locally_rolled_back_vertices)
                    intersection_guard_locally_rolled_back_vertices=rolled;
                if(o.verbose)fprintf(stderr,
                    "[relax] dark-patch transaction rejects attempt %zu: "
                    "retained patches/vertices %zu/%zu (mean scale %.3f); "
                    "filled dark %zu->%zu, fitted dark %zu->%zu, missing %zu->%zu; "
                    "conflicts %zu->%zu, folds %zu->%zu, orientation %zu->%zu\n",
                    attempt,active_patches,active_vertices,
                    active_vertices?scale_sum/(double)active_vertices:0.0,
                    patch_stats.dark_before,candidate_bake.filled_dark,
                    patch_stats.fitted_dark_before,candidate_bake.fitted_dark,
                    patch_stats.missing_before,candidate_bake.missing,
                    intersection_current.conflicts,candidate_ix.conflicts,
                    intersection_current.fold_pairs,candidate_ix.fold_pairs,
                    orientation_base,orientation_bad);
                if(!active_vertices)break;
            }
            if(accepted) {
                patch_stats.transaction_accepted=1;
                double scale_sum=0.0;
                for(size_t k=0;k<patch_stats.regions;k++)if(patch_alive[k]) {
                    patch_stats.applied++;
                    patch_stats.attenuated+=patch_scale[k]<0.999f;
                }
                for(size_t v=0;v<nv;v++)if(patch_moved[v]&&cand_t[v]>0.0f) {
                    patch_stats.vertices_retained++;scale_sum+=(double)cand_t[v];
                }
                patch_stats.geometry_rejected=patch_stats.proposed-patch_stats.applied;
                patch_stats.retained_scale=patch_stats.vertices?
                    scale_sum/(double)patch_stats.vertices:0.0;
                patch_stats.dark_after=accepted_bake.filled_dark;
                if(o.prevent_intersection_growth) {
                    intersection_current=accepted_ix;
                    size_t *swap=intersection_degree_current;
                    intersection_degree_current=intersection_degree_candidate;
                    intersection_degree_candidate=swap;
                }
            } else {
                memcpy(verts,patch_before,nv*3*sizeof *verts);
                patch_stats.geometry_rejected=patch_stats.proposed;
                patch_stats.dark_after=patch_stats.dark_before;
                if(!patch_failed)intersection_guard_round_rollbacks++;
            }
        } else if(!patch_failed)patch_stats.dark_after=patch_stats.dark_before;
        if(o.verbose&&!patch_failed)fprintf(stderr,
            "[relax] dark-patch recovery: regions %zu (%zu small), "
            "border/ray consensus %zu/%zu (selected border/ambiguous %zu/%zu), "
            "wobbly/no-consensus/overlap %zu/%zu/%zu; "
            "proposed/applied/attenuated/geometry-rejected %zu/%zu/%zu/%zu patches "
            "(%zu proposed, %zu retained vertices), filled-dark %zu->%zu, "
            "transaction %s scale %.3f\n",
            patch_stats.regions,patch_stats.regions_small,
            patch_stats.border_consensus,patch_stats.ray_consensus,
            patch_stats.border_selected,patch_stats.ambiguous_selected,
            patch_stats.wobbly,patch_stats.no_consensus,patch_stats.overlap_skipped,
            patch_stats.proposed,patch_stats.applied,patch_stats.attenuated,
            patch_stats.geometry_rejected,
            patch_stats.vertices,patch_stats.vertices_retained,
            patch_stats.dark_before,patch_stats.dark_after,
            patch_stats.transaction_accepted?"accepted":"no-op/rolled-back",
            patch_stats.retained_scale);
        if(patch_failed) {
            fprintf(stderr,"ERROR: connected dark-patch recovery failed; restored pre-patch geometry\n");
            if(have_patch_before)memcpy(verts,patch_before,nv*3*sizeof *verts);
            free(patch_moved);free(patch_id);free(patch_alive);free(patch_reduce);free(patch_scale);
            if(own_patch_arrays){free(patch_before);free(patch_proposed);}
            free(init_nrm);
            free(orig); free(target); free(hast); free(frozen); free(anchored); free(bake_dark_seed);
            free(cand); free(cand_t); free(cand_w);
            free(guide); free(conf); free(face_offset); free(face_weight); free(winding_limit);
            free(intersection_before); free(intersection_proposed);
            free(intersection_degree_current); free(intersection_degree_candidate);
            free(R); free(nxt); free(scalar_nxt); return -1;
        }
        free(patch_moved);free(patch_id);free(patch_alive);free(patch_reduce);free(patch_scale);
        if(own_patch_arrays){free(patch_before);free(patch_proposed);}
    } else if(o.patch_recover&&o.verbose) {
        fprintf(stderr,"[relax] dark-patch recovery skipped: requires cylindrical/oriented ribbon\n");
    }
    {
        SnapBakeCounts output_bake={0};
        float *final_nrm=out_stats?MeshNormals_compute(verts,nv,faces,nf):NULL;
        if(final_nrm&&adj.structured)output_bake=structured_bake_counts(&ct,verts,
            final_nrm,filled,uv,adj.H,adj.W,o.bake_window_low,o.bake_window_high,
            o.bake_dark_u8);
        else if(final_nrm)output_bake=sparse_bake_scan(&ct,verts,final_nrm,faces,nf,
            filled,nv,o.bake_window_low,o.bake_window_high,o.bake_dark_u8,NULL,NULL);
        free(final_nrm);
        RemeshAudit ra = remesh_audit(verts, orig, init_nrm, uv, faces, nf, nv,
                                      grid_h, grid_w,
                                      o.remesh_motion, o.remesh_motion_max,
                                      o.remesh_sander, o.remesh_over4);
        if (out_remesh_advised) *out_remesh_advised = ra.advised;
        if (out_stats) {
            out_stats->displacement_rms = ra.displacement_rms;
            out_stats->tangential_p95 = ra.tangential_p95;
            out_stats->tangential_max = ra.tangential_max;
            out_stats->sander_p95 = ra.sander_p95;
            out_stats->sander_p99 = ra.sander_p99;
            out_stats->sander_max = ra.sander_max;
            out_stats->sander_samples = ra.triangles;
            out_stats->quads = ra.quads;
            out_stats->over2_fraction = ra.triangles ? (double)ra.over2/(double)ra.triangles : 0.0;
            out_stats->over4_fraction = ra.triangles ? (double)ra.over4/(double)ra.triangles : 0.0;
            out_stats->quad_sander_p95 = ra.quad_sander_p95;
            out_stats->quad_sander_p99 = ra.quad_sander_p99;
            out_stats->quad_sander_max = ra.quad_sander_max;
            out_stats->quad_over2_fraction = ra.quads ? (double)ra.quad_over2/(double)ra.quads : 0.0;
            out_stats->quad_over4_fraction = ra.quads ? (double)ra.quad_over4/(double)ra.quads : 0.0;
            out_stats->flipped = ra.flipped;
            out_stats->degenerate = ra.degenerate;
            out_stats->dark_probe_threshold = dark_max;
            out_stats->fitted_bake_dark = nfit_dark;
            out_stats->fitted_dark_recto_target = nfit_dark_target;
            out_stats->filled_bake_dark = nfill_dark;
            out_stats->filled_dark_recto_target = nfill_dark_target;
            out_stats->input_bake_fitted_samples=input_bake.fitted_samples;
            out_stats->input_bake_filled_samples=input_bake.filled_samples;
            out_stats->input_bake_fitted_dark=input_bake.fitted_dark;
            out_stats->input_bake_filled_dark=input_bake.filled_dark;
            out_stats->input_bake_missing=input_bake.missing;
            out_stats->output_bake_fitted_samples=output_bake.fitted_samples;
            out_stats->output_bake_filled_samples=output_bake.filled_samples;
            out_stats->output_bake_fitted_dark=output_bake.fitted_dark;
            out_stats->output_bake_filled_dark=output_bake.filled_dark;
            out_stats->output_bake_missing=output_bake.missing;
            out_stats->input_sander_p95=input_ra.sander_p95;
            out_stats->input_sander_p99=input_ra.sander_p99;
            out_stats->input_over4_fraction=input_ra.triangles?
                (double)input_ra.over4/(double)input_ra.triangles:0.0;
            out_stats->quilt_mg_levels = final_quilt_mg.levels_used;
            out_stats->quilt_mg_corrections = final_quilt_mg.corrections;
            out_stats->quilt_mg_residual_before = final_quilt_mg.residual_before;
            out_stats->quilt_mg_residual_after = final_quilt_mg.residual_after;
            out_stats->position_mg_levels = final_position_mg.levels_used;
            out_stats->position_mg_corrections = final_position_mg.corrections;
            out_stats->position_mg_backtracks = final_position_mg.backtracks;
            out_stats->position_mg_rejects = final_position_mg.rejects;
            out_stats->position_mg_patches_accepted=final_position_mg.patches_accepted;
            out_stats->position_mg_patches_rejected=final_position_mg.patches_rejected;
            out_stats->position_mg_residual_before = final_position_mg.residual_before;
            out_stats->position_mg_residual_after = final_position_mg.residual_after;
        }
        if (o.verbose)
        fprintf(stderr,
                "[relax] remesh audit: displacement-rms %.3f vox, tangential p95/max "
                "%.3f/%.3f vox; triangle Sander p95/p99/max %.3f/%.3f/%.3f, "
                ">=2x %.3f%%, >=4x %.3f%% (%zu samples), flipped %zu, degenerate %zu; "
                "quad-average p95/p99 %.3f/%.3f -> %s\n",
                ra.displacement_rms, ra.tangential_p95, ra.tangential_max,
                ra.sander_p95, ra.sander_p99, ra.sander_max,
                ra.triangles ? 100.0*(double)ra.over2/(double)ra.triangles : 0.0,
                ra.triangles ? 100.0*(double)ra.over4/(double)ra.triangles : 0.0,
                ra.triangles, ra.flipped, ra.degenerate,
                ra.quad_sander_p95, ra.quad_sander_p99,
                ra.advised ? "REMESH ADVISED" : "keep current lattice");
    }

    /* Re-sample the mesh that actually leaves the pass.  The calibration
     * histogram above describes the pass input and is therefore unsuitable for
     * deciding whether an iterative candidate earned its geometric cost. */
    if (out_stats) {
        size_t hf[256] = {0}, hh[256] = {0}, nfitted = 0, nfilled = 0;
        for (size_t v = 0; v < nv; v++) {
            double I = sample_trilinear(&ct, verts[v*3+0], verts[v*3+1], verts[v*3+2]);
            if (I < 0.0) continue;
            int bin = (int)lround(I); if (bin < 0) bin = 0; if (bin > 255) bin = 255;
            if (filled[v]) { hh[bin]++; nfilled++; }
            else { hf[bin]++; nfitted++; }
        }
        out_stats->fitted_ct_samples = nfitted;
        out_stats->filled_ct_samples = nfilled;
        out_stats->fitted_ct_p1 = hist_percentile(hf, nfitted, 1.0);
        out_stats->fitted_ct_p25 = hist_percentile(hf, nfitted, 25.0);
        out_stats->fitted_ct_p50 = hist_percentile(hf, nfitted, 50.0);
        out_stats->fitted_ct_p75 = hist_percentile(hf, nfitted, 75.0);
        out_stats->fitted_ct_p99 = hist_percentile(hf, nfitted, 99.0);
        out_stats->filled_ct_p1 = hist_percentile(hh, nfilled, 1.0);
        out_stats->filled_ct_p25 = hist_percentile(hh, nfilled, 25.0);
        out_stats->filled_ct_p50 = hist_percentile(hh, nfilled, 50.0);
        out_stats->filled_ct_p75 = hist_percentile(hh, nfilled, 75.0);
        out_stats->filled_ct_p99 = hist_percentile(hh, nfilled, 99.0);
        out_stats->snapped = snapped;
        out_stats->nopap = nopap;
        out_stats->raw_ridge = final_ncand;
        out_stats->recto_ridge = final_ncand_recto;
        out_stats->band_fallback = final_band_fallback;
        out_stats->band_abstain = final_band_abstain;
        out_stats->near_zero_ridge = final_ncand_near;
        out_stats->raw_no_candidate = final_raw_no_candidate;
        out_stats->raw_candidate_applied = final_raw_applied;
        out_stats->quilt_only_applied = final_quilt_only;
        out_stats->guided_ridge = final_guided;
        out_stats->wide_ridge = final_wide;
        out_stats->local_contrast_ridge = final_local;
        out_stats->winding_rejected = winding_rejected_total;
        out_stats->dark_seed_forced_movable = (int)nfit_dark_forced;
        out_stats->dark_seed_quilt_only = final_dark_seed_quilt_only;
        out_stats->dark_face_samples=final_face_counts.dark_samples;
        out_stats->dark_face_candidates=final_face_counts.candidates;
        out_stats->dark_face_outward_candidates=final_face_counts.outward_candidates;
        out_stats->dark_face_no_candidate=final_face_counts.no_candidate;
        out_stats->adaptive_winding_vertices=winding_supported;
            out_stats->patch_regions=patch_stats.regions;
            out_stats->patch_regions_small=patch_stats.regions_small;
            out_stats->patch_border_consensus=patch_stats.border_consensus;
            out_stats->patch_ray_consensus=patch_stats.ray_consensus;
            out_stats->patch_border_selected=patch_stats.border_selected;
            out_stats->patch_ambiguous_selected=patch_stats.ambiguous_selected;
            out_stats->patch_wobbly=patch_stats.wobbly;
            out_stats->patch_no_consensus=patch_stats.no_consensus;
            out_stats->patch_overlap_skipped=patch_stats.overlap_skipped;
            out_stats->patch_proposed=patch_stats.proposed;
            out_stats->patch_applied=patch_stats.applied;
            out_stats->patch_attenuated=patch_stats.attenuated;
            out_stats->patch_geometry_rejected=patch_stats.geometry_rejected;
            out_stats->patch_vertices=patch_stats.vertices;
            out_stats->patch_vertices_retained=patch_stats.vertices_retained;
            out_stats->patch_dark_before=patch_stats.dark_before;
            out_stats->patch_dark_after=patch_stats.dark_after;
            out_stats->patch_transaction_accepted=patch_stats.transaction_accepted;
            out_stats->patch_retained_scale=patch_stats.retained_scale;
            out_stats->adaptive_winding_p05=winding_pct[0];
            out_stats->adaptive_winding_p50=winding_pct[1];
            out_stats->adaptive_winding_p95=winding_pct[2];
            out_stats->intersection_guard_audits=intersection_guard_audits;
            out_stats->intersection_guard_backtracks=intersection_guard_backtracks;
            out_stats->intersection_guard_round_rollbacks=intersection_guard_round_rollbacks;
            out_stats->intersection_guard_orientation_rejects=
                intersection_guard_orientation_rejects;
            out_stats->intersection_guard_locally_rolled_back_vertices=
                intersection_guard_locally_rolled_back_vertices;
            if(o.prevent_intersection_growth) {
                out_stats->intersection_audit_complete=1;
                out_stats->input_intersections=intersection_input.conflicts;
                out_stats->input_overlap=intersection_input.overlap_pairs;
                out_stats->input_stab=intersection_input.stab_pairs;
                out_stats->input_fold=intersection_input.fold_pairs;
                out_stats->output_intersections=intersection_current.conflicts;
                out_stats->output_overlap=intersection_current.overlap_pairs;
                out_stats->output_stab=intersection_current.stab_pairs;
                out_stats->output_fold=intersection_current.fold_pairs;
            }
            for(size_t v=0;v<nv;v++) {
            double drift=fabs((double)verts[v*3]-(double)orig[v*3]);
            if(drift>out_stats->axial_drift_max)out_stats->axial_drift_max=drift;
        }
        out_stats->quilt_downweighted = final_candidate_downweighted;
        out_stats->quilt_candidate_weight_mean = final_candidate_weight_mean;
        out_stats->quilt_raw_edges = final_neraw;
        out_stats->quilt_edges = final_nequilt;
        out_stats->quilt_raw_edge_rms = final_neraw ? sqrt(final_eraw/(double)final_neraw) : 0.0;
        out_stats->quilt_edge_rms = final_nequilt ? sqrt(final_equilt/(double)final_nequilt) : 0.0;
        out_stats->quilt_seam_smooth = o.quilt_seam_smooth;
        quilt_edge_finish(&final_raw_regions, &out_stats->quilt_raw_regions);
        quilt_edge_finish(&final_solved_regions, &out_stats->quilt_solved_regions);
        quilt_edge_finish(&final_applied_regions, &out_stats->quilt_applied_regions);
        if (o.verbose) {
            fprintf(stderr, "[relax] final direct fitted CT p1/p25/p50/p75/p99=%d/%d/%d/%d/%d (%zu)\n",
                    out_stats->fitted_ct_p1, out_stats->fitted_ct_p25,
                    out_stats->fitted_ct_p50, out_stats->fitted_ct_p75,
                    out_stats->fitted_ct_p99, nfitted);
            fprintf(stderr, "[relax] final direct fill CT p1/p25/p50/p75/p99=%d/%d/%d/%d/%d (%zu)\n",
                    out_stats->filled_ct_p1, out_stats->filled_ct_p25,
                    out_stats->filled_ct_p50, out_stats->filled_ct_p75,
                    out_stats->filled_ct_p99, nfilled);
            fprintf(stderr,
                    "[relax] quilt solved RMS all/fit/fill/seam %.3f/%.3f/%.3f/%.3f "
                    "(%zu/%zu/%zu/%zu edges); applied %.3f/%.3f/%.3f/%.3f "
                    "(%zu/%zu/%zu/%zu edges)\n",
                    out_stats->quilt_solved_regions.all_rms,
                    out_stats->quilt_solved_regions.fitted_fitted_rms,
                    out_stats->quilt_solved_regions.filled_filled_rms,
                    out_stats->quilt_solved_regions.fitted_filled_rms,
                    out_stats->quilt_solved_regions.all_edges,
                    out_stats->quilt_solved_regions.fitted_fitted_edges,
                    out_stats->quilt_solved_regions.filled_filled_edges,
                    out_stats->quilt_solved_regions.fitted_filled_edges,
                    out_stats->quilt_applied_regions.all_rms,
                    out_stats->quilt_applied_regions.fitted_fitted_rms,
                    out_stats->quilt_applied_regions.filled_filled_rms,
                    out_stats->quilt_applied_regions.fitted_filled_rms,
                    out_stats->quilt_applied_regions.all_edges,
                    out_stats->quilt_applied_regions.fitted_fitted_edges,
                    out_stats->quilt_applied_regions.filled_filled_edges,
                    out_stats->quilt_applied_regions.fitted_filled_edges);
        }
    }
    free(init_nrm);
    free(orig); free(target); free(hast); free(frozen); free(anchored); free(bake_dark_seed);
    free(cand); free(cand_t); free(cand_w);
    free(guide); free(conf); free(face_offset); free(face_weight); free(winding_limit);
    free(intersection_before); free(intersection_proposed);
    free(intersection_degree_current); free(intersection_degree_candidate);
    free(R); free(nxt); free(scalar_nxt);
    if (out_snapped) *out_snapped = snapped;
    if (out_nopap) *out_nopap = nopap;
    (void)nload;
    return 0;
}

/* ============================================================================
 * Self-test -- inject a synthetic 1-cube RAW table with a bright plane and check
 * an off-plane filled vertex snaps onto it, while a vertex in uniform dark does
 * not.  The CubeTable is public precisely so tests can inject slots.
 * ==========================================================================*/
#define CK(c, m) do { if (!(c)) { fprintf(stderr, "  FAIL: %s\n", (m)); fails++; } \
                      else fprintf(stderr, "  ok: %s\n", (m)); } while (0)

int QuadStripSnap_selftest(void) {
    int fails = 0;
    Arena_T arena = Arena_new();
    fprintf(stderr, "[selftest] quad_strip_snap\n");

    const long C = 16;
    CubeTable ct;
    memset(&ct, 0, sizeof ct);
    ct.arena = arena; ct.chunk = C;
    ct.cz0 = ct.cy0 = ct.cx0 = 0; ct.nz = ct.ny = ct.nx = 1;
    ct.slot = (uint8_t **)ARENA_CALLOC(arena, 1, sizeof(uint8_t *));
    uint8_t *vol = ARENA_ALLOC(arena, (size_t)C * C * C);
    for (long z = 0; z < C; z++) for (long y = 0; y < C; y++) for (long x = 0; x < C; x++) {
        int val = 40;                 /* dark resin */
        if (x == 8) val = 210;        /* bright papyrus ridge at x=8 */
        else if (x == 7 || x == 9) val = 120;
        vol[((size_t)z * C + (size_t)y) * C + (size_t)x] = (uint8_t)val;
    }
    ct.slot[0] = vol; ct.n_loaded = 1;

    /* t0: fitted claims remain finite-weight observations.  A high configured
     * multiplier is deliberately distinct from both ordinary vertices and the
     * later-pass recovered-fill anchor; zero retains legacy hard-freeze mode. */
    {
        QuadStripSnapOpts ao; QuadStripSnap_defaults(&ao);
        float persistent[6]={10.0f,11.0f,12.0f,20.0f,21.0f,22.0f};
        float pass[6]={30.0f,31.0f,32.0f,40.0f,41.0f,42.0f};
        ao.filled_target_anchor=4.0;ao.fitted_source_anchor=64.0;
        ao.fitted_source_positions=persistent;
        CK(fabs(snap_orig_anchor_multiplier(&ao,0)-1.0)<1e-12&&
           fabs(snap_orig_anchor_multiplier(&ao,1)-4.0)<1e-12&&
           fabs(snap_orig_anchor_multiplier(&ao,2)-64.0)<1e-12&&
           snap_orig_anchor_coordinate(&ao,pass,1,2,2)==22.0&&
           snap_orig_anchor_coordinate(&ao,pass,1,2,0)==42.0,
           "t0 ordinary, recovered-fill, and fitted-source anchors have distinct soft weights");
    }

    /* t1: a filled vertex 3 vox off the ridge, normal along +x, snaps to x~8. */
    {
        double p[3] = { 8.0, 8.0, 5.0 }, n[3] = { 0.0, 0.0, 1.0 }, t = 0.0;
        int hit = snap_one_vertex(&ct, p, n, 6.0, 0.5, 130.0, 20.0, 0.75, &t);
        CK(hit == 1, "t1 off-ridge vertex finds a ridge");
        CK(fabs((p[2] + t) - 8.0) < 0.6, "t1 snaps onto the bright plane (x~8)");
    }
    /* t2: nearest-ridge preference -- start between two planes, pick the closer. */
    {
        /* add a second ridge at x=2 */
        for (long z = 0; z < C; z++) for (long y = 0; y < C; y++)
            vol[((size_t)z * C + (size_t)y) * C + 2] = 210;
        double p[3] = { 8.0, 8.0, 6.0 }, n[3] = { 0.0, 0.0, 1.0 }, t = 0.0;  /* x=6: nearer x=8 */
        int hit = snap_one_vertex(&ct, p, n, 6.0, 0.5, 130.0, 20.0, 0.75, &t);
        CK(hit == 1 && fabs((p[2] + t) - 8.0) < 0.6, "t2 picks the nearer ridge (x=8 not x=2)");
    }
    /* t3: uniform-dark region -> no ridge, no move. */
    {
        for (long z = 0; z < C; z++) for (long y = 0; y < C; y++) for (long x = 0; x < C; x++)
            vol[((size_t)z * C + (size_t)y) * C + (size_t)x] = 60;
        double p[3] = { 8.0, 8.0, 8.0 }, n[3] = { 0.0, 0.0, 1.0 }, t = 0.0;
        int hit = snap_one_vertex(&ct, p, n, 6.0, 0.5, 130.0, 20.0, 0.75, &t);
        CK(hit == 0, "t3 uniform dark -> no snap (honest-black candidate)");
    }
    /* t4: fitted-data trust is ridge-relative, not an absolute brightness gate. */
    {
        for (long z = 0; z < C; z++) for (long y = 0; y < C; y++) for (long x = 0; x < C; x++)
            vol[((size_t)z * C + (size_t)y) * C + (size_t)x] = (uint8_t)(x == 8 ? 95 : 60);
        double p[3] = { 8.0, 8.0, 8.0 }, n[3] = { 0.0, 0.0, 1.0 };
        RidgeProfile rp = ridge_profile(&ct, p, n, 6.0, 0.5, 130.0, 20.0, 0.75, 0, 0, 0, QUAD_SNAP_SIDE_NEAREST, 0.0);
        CK(rp.on_ridge == 1, "t4 dim fitted ridge remains trusted below bright_min");
        CK(rp.has_target == 0, "t4 dim ridge is not misreported as a snap target");
    }
    /* t5: an inward-oriented profile selects the recto bright->dark shoulder,
     * rather than the arbitrary nearest maximum of a thick papyrus band. */
    {
        for (long z = 0; z < C; z++) for (long y = 0; y < C; y++) for (long x = 0; x < C; x++)
            vol[((size_t)z * C + (size_t)y) * C + (size_t)x] =
                (uint8_t)(x >= 6 && x <= 8 ? 210 : 40);
        double p[3] = { 8.0, 8.0, 7.0 }, n[3] = { 0.0, 0.0, 1.0 };
        RidgeProfile rp = ridge_profile(&ct, p, n, 4.0, 0.5, 130.0, 20.0, 0.75, 1, 0, 0, QUAD_SNAP_SIDE_RECTO, 0.0);
        CK(rp.has_target == 1 && rp.used_recto == 1,
           "t5 inward profile identifies a recto transition");
        CK(p[2] + rp.target_t >= 7.9,
           "t5 recto target is on the inward shoulder of the bright band");
    }
    /* t5b: distinguish a dim local ridge from a substantially brighter recto
     * sheet farther inward.  The diagnostic feeds explicit bake/refit masks;
     * it does not unpin geometry merely because an intensity is dark. */
    {
        for (long z = 0; z < C; z++) for (long y = 0; y < C; y++) for (long x = 0; x < C; x++) {
            int val = 40;
            if (x == 4) val = 80;
            if (x >= 7 && x <= 8) val = 210;
            vol[((size_t)z * C + (size_t)y) * C + (size_t)x] = (uint8_t)val;
        }
        double p[3] = { 8.0, 8.0, 4.0 }, n[3] = { 0.0, 0.0, 1.0 };
        RidgeProfile rp = ridge_profile(&ct, p, n, 6.0, 0.5, 130.0, 20.0, 0.75, 1, 0, 0, QUAD_SNAP_SIDE_RECTO, 0.0);
        CK(fabs(rp.nearest_peak_t) <= 0.5 && rp.has_target == 1,
           "t5b dim local ridge coexists with a brighter recto target");
        CK(dark_inward_recto_target(&rp, 1, 85.0, 130.0, 20.0, 0.75) == 1,
           "t5b bake-dark probe reports the farther inward target");
        CK(dark_inward_recto_target(&rp, -1, 85.0, 130.0, 20.0, 0.75) == 0,
           "t5b probe never reports the same target outward");
    }
    /* t5c: burnt material can preserve a strong local transition even when its
     * absolute CT values sit below the scroll-wide bright floor. */
    {
        for (long z=0;z<C;z++) for (long y=0;y<C;y++) for (long x=0;x<C;x++) {
            int val=50;
            if(x==8)val=100;else if(x==7||x==9)val=78;
            vol[((size_t)z*C+(size_t)y)*C+(size_t)x]=(uint8_t)val;
        }
        double p[3]={8.0,8.0,5.0},n[3]={0.0,0.0,1.0};
        RidgeProfile global=ridge_profile(&ct,p,n,6.0,0.5,130.0,20.0,0.75,0,0,0,QUAD_SNAP_SIDE_NEAREST,0.0);
        RidgeProfile local=ridge_profile(&ct,p,n,6.0,0.5,130.0,20.0,0.75,0,1,0,QUAD_SNAP_SIDE_NEAREST,0.0);
        CK(!global.has_target&&local.has_target&&local.used_local_contrast,
           "t5c robust local contrast recovers a dim supported ridge");
    }
    /* t5d: a bad/tangential supplied normal may miss the sheet entirely.  The
     * wider search blends toward radial inward, but the initializer-anchored
     * winding corridor vetoes a target beyond the same-wrap half-width. */
    {
        for (long z=0;z<C;z++) for (long y=0;y<C;y++) for (long x=0;x<C;x++) {
            int val=40;
            if(x>=6&&x<=7)val=210;
            vol[((size_t)z*C+(size_t)y)*C+(size_t)x]=(uint8_t)val;
        }
        QuadStripSnapOpts go;QuadStripSnap_defaults(&go);
        go.snap_side=QUAD_SNAP_SIDE_RECTO;go.axis_y=8.0;go.axis_x=8.0;
        go.reach=1.0;go.guided_reach=6.0;go.winding_tube=3.1;
        double p[3]={8.0,8.0,4.0},n[3]={1.0,0.0,0.0},dir[3];
        int sign=0,rejected=0;
        RidgeProfile guided=guided_ridge_profile(&ct,&go,p,n,130.0,20.0,0.75,
                                                  80.0,0,p,3.1,dir,&sign,&rejected);
        CK(guided.has_target&&guided.used_guided&&sign==1&&dir[2]>0.5,
           "t5d guided search recovers a dark site with a tangential normal");
        go.winding_tube=0.5;rejected=0;
        RidgeProfile vetoed=guided_ridge_profile(&ct,&go,p,n,130.0,20.0,0.75,
                                                  80.0,0,p,0.5,dir,&sign,&rejected);
        CK(!vetoed.has_target&&rejected>0,
           "t5d winding corridor rejects a neighbouring-phase candidate");
    }
    /* t5f: when a bad initializer lies in the inward gap, its own recto is
     * OUTWARD.  The former rp.target_t>ridge_lock gate threw this supported
     * candidate away and preferred darkness (or a different inward wrap). */
    {
        for (long z=0;z<C;z++) for (long y=0;y<C;y++) for (long x=0;x<C;x++) {
            int val=(x>=3&&x<=4)?210:40;
            vol[((size_t)z*C+(size_t)y)*C+(size_t)x]=(uint8_t)val;
        }
        QuadStripSnapOpts go;QuadStripSnap_defaults(&go);
        go.snap_side=QUAD_SNAP_SIDE_RECTO;go.preserve_axial=1;
        go.axis_y=8.0;go.axis_x=8.0;go.reach=1.0;go.guided_reach=6.0;
        double p[3]={8.0,8.0,6.0},n[3]={1.0,0.0,0.0},dir[3];
        int sign=0,rejected=0;
        RidgeProfile rp=guided_ridge_profile(&ct,&go,p,n,130.0,20.0,0.75,
            80.0,1,p,3.0,dir,&sign,&rejected);
        CK(rp.has_target&&rp.used_guided&&rp.used_outward&&rp.target_t<0.0,
           "t5f bake-dark search retains a supported outward recto");
        CK(fabs(dir[0])<1e-9&&p[2]+rp.target_t*dir[2]<5.0,
           "t5f cylindrical search changes radius/phase but not axial z");
    }
    /* t5g: local winding limits expand beyond three voxels when the measured
     * pitch permits it, yet remain strictly below half the adjacent-turn gap. */
    {
        enum { H=3,W=150,N=H*W };
        float init[N*3];
        for(int r=0;r<H;r++)for(int c=0;c<W;c++) {
            size_t v=(size_t)r*W+(size_t)c;
            double a=0.10*(double)c;
            double rad=30.0+10.0*a/6.28318530717958647693;
            init[v*3]=(float)r;
            init[v*3+1]=(float)(rad*cos(a));
            init[v*3+2]=(float)(rad*sin(a));
        }
        size_t supported=0;double pct[3]={0};
        float *guard=adaptive_winding_limits(init,H,W,0.0,0.0,3.0,12.0,
                                              &supported,pct);
        CK(guard&&supported>(size_t)N/2,
           "t5g adaptive winding guard finds adjacent-turn support");
        CK(pct[1]>4.0&&pct[1]<5.0,
           "t5g adaptive corridor is 45 percent of local pitch");
        free(guard);
    }
    /* t5h: a face can bake dark even though corner-only profiles are useless.
     * Search the two real raster samples and distribute their physical target
     * vectors to the quad corners. */
    {
        for(long z=0;z<C;z++)for(long y=0;y<C;y++)for(long x=0;x<C;x++) {
            int val=(x>=6&&x<=7)?210:40;
            vol[((size_t)z*C+(size_t)y)*C+(size_t)x]=(uint8_t)val;
        }
        float v[12]={8,7,4, 8,9,4, 9,7,4, 9,9,4};
        float n[12]={1,0,0, 1,0,0, 1,0,0, 1,0,0};
        float uv[8]={0,0, 2,0, 0,1, 2,1};
        uint8_t fill[4]={1,1,1,1},frozen[4]={0,0,0,0};
        float guard[4]={4,4,4,4},sum[12],weight[4];
        QuadStripSnapOpts go;QuadStripSnap_defaults(&go);
        go.snap_side=QUAD_SNAP_SIDE_RECTO;go.preserve_axial=1;go.axis_y=8;go.axis_x=8;
        go.reach=1;go.guided_reach=6;go.bright_min=130;
        int rejected=0;
        SnapFaceTargetCounts fc=structured_dark_face_targets(
            &ct,&go,v,v,n,fill,frozen,uv,guard,2,2,130.0,sum,weight,&rejected,
            NULL);
        CK(fc.dark_samples==2&&fc.candidates==2,
           "t5h dark emitted-triangle samples find face targets");
        CK(weight[0]>0.0f&&fabs((double)sum[0])<1e-9&&sum[2]>0.0f,
           "t5h face target is distributed at fixed axial z");
    }
    /* t12a: midline mode targets the CENTER of a bright band, not its recto
     * shoulder.  Band x in [6,8]; vertex at x=6 must move to x~7. */
    {
        for (long z=0;z<C;z++) for (long y=0;y<C;y++) for (long x=0;x<C;x++)
            vol[((size_t)z*C+(size_t)y)*C+(size_t)x] =
                (uint8_t)(x>=6&&x<=8?210:40);
        double p[3]={8.0,8.0,6.0},n[3]={0.0,0.0,1.0};
        RidgeProfile rp=ridge_profile(&ct,p,n,6.0,0.5,130.0,20.0,0.75,1,0,0,
                                      QUAD_SNAP_SIDE_MIDLINE,5.0);
        CK(rp.has_target==1&&rp.used_midline==1&&rp.transition_drop>0.0,
           "t12a midline finds the two-edged band center");
        CK(fabs(p[2]+rp.target_t-7.0)<0.6,
           "t12a target is the band midpoint (x~7), not the recto shoulder");
        CK(rp.band_half_width>0.5,
           "t12a band half width is measured");
    }
    /* t12b: a band merged with its neighbour (no dark separator within the
     * corridor) falls back to the bounded smoothed plateau midpoint. */
    {
        for (long z=0;z<C;z++) for (long y=0;y<C;y++) for (long x=0;x<C;x++)
            vol[((size_t)z*C+(size_t)y)*C+(size_t)x] =
                (uint8_t)(x>=3&&x<=11?210:40);
        double p[3]={8.0,8.0,7.0},n[3]={0.0,0.0,1.0};
        RidgeProfile rp=ridge_profile(&ct,p,n,6.0,0.5,130.0,20.0,0.75,1,0,0,
                                      QUAD_SNAP_SIDE_MIDLINE,3.0);
        CK(rp.has_target==1&&rp.used_band_fallback==1&&fabs(rp.target_t)<0.6,
           "t12b merged band falls back to the corridor-bounded plateau midpoint");
    }
    /* t12c: a band clipped by the corridor keeps a bounded fallback target. */
    {
        for (long z=0;z<C;z++) for (long y=0;y<C;y++) for (long x=0;x<C;x++)
            vol[((size_t)z*C+(size_t)y)*C+(size_t)x] =
                (uint8_t)(x>=6&&x<=8?210:40);
        double p[3]={8.0,8.0,6.0},n[3]={0.0,0.0,1.0};
        RidgeProfile rp=ridge_profile(&ct,p,n,6.0,0.5,130.0,20.0,0.75,1,0,0,
                                      QUAD_SNAP_SIDE_MIDLINE,1.5);
        CK(rp.has_target==1&&rp.used_band_fallback==1&&
           fabs(rp.target_t)<=1.5+1e-6,
           "t12c clipped band keeps a corridor-bounded fallback target");
    }
    /* t12d: an all-dark corridor abstains entirely. */
    {
        memset(vol,40,(size_t)C*C*C);
        double p[3]={8.0,8.0,8.0},n[3]={0.0,0.0,1.0};
        RidgeProfile rp=ridge_profile(&ct,p,n,6.0,0.5,130.0,20.0,0.75,1,0,0,
                                      QUAD_SNAP_SIDE_MIDLINE,4.0);
        CK(rp.has_target==0,
           "t12d all-dark corridor abstains (honest-dark candidate)");
    }
    /* t12e: an all-bright corridor centers in place (already on the midline). */
    {
        memset(vol,210,(size_t)C*C*C);
        double p[3]={8.0,8.0,8.0},n[3]={0.0,0.0,1.0};
        RidgeProfile rp=ridge_profile(&ct,p,n,6.0,0.5,130.0,20.0,0.75,1,0,0,
                                      QUAD_SNAP_SIDE_MIDLINE,4.0);
        CK(rp.has_target==1&&rp.used_band_fallback==1&&rp.on_ridge==1,
           "t12e all-bright plateau centers in place (on_ridge)");
    }
    /* t12f: corridor bounding -- a mid-gap vertex with both neighbouring bands
     * outside the corridor abstains; widening the corridor admits a bounded
     * target. */
    {
        for (long z=0;z<C;z++) for (long y=0;y<C;y++) for (long x=0;x<C;x++) {
            int val=40;
            if((x>=2&&x<=3)||(x>=9&&x<=10))val=210;
            vol[((size_t)z*C+(size_t)y)*C+(size_t)x]=(uint8_t)val;
        }
        double p[3]={8.0,8.0,6.0},n[3]={0.0,0.0,1.0};
        RidgeProfile tight=ridge_profile(&ct,p,n,6.0,0.5,130.0,20.0,0.75,1,0,0,
                                         QUAD_SNAP_SIDE_MIDLINE,2.0);
        RidgeProfile wide=ridge_profile(&ct,p,n,6.0,0.5,130.0,20.0,0.75,1,0,0,
                                        QUAD_SNAP_SIDE_MIDLINE,4.0);
        CK(tight.has_target==0,
           "t12f corridor rejects bands beyond the winding half-width");
        CK(wide.has_target==1&&fabs(wide.target_t)<=4.0+1e-6,
           "t12f widened corridor admits a bounded band target");
    }
    /* t12g: the dark-site probe reports a distant band center and never an
     * already-centered vertex. */
    {
        for (long z=0;z<C;z++) for (long y=0;y<C;y++) for (long x=0;x<C;x++)
            vol[((size_t)z*C+(size_t)y)*C+(size_t)x] =
                (uint8_t)(x>=6&&x<=8?210:40);
        double pd[3]={8.0,8.0,3.0},n[3]={0.0,0.0,1.0};
        RidgeProfile far_rp=ridge_profile(&ct,pd,n,6.0,0.5,130.0,20.0,0.75,1,0,0,
                                          QUAD_SNAP_SIDE_MIDLINE,6.0);
        CK(dark_band_center_target(&far_rp,1,85.0,130.0,20.0,0.75)==1,
           "t12g dark probe reports the distant band center");
        double pc[3]={8.0,8.0,7.0};
        RidgeProfile ctr=ridge_profile(&ct,pc,n,6.0,0.5,130.0,20.0,0.75,1,0,0,
                                       QUAD_SNAP_SIDE_MIDLINE,5.0);
        CK(dark_band_center_target(&ctr,1,85.0,130.0,20.0,0.75)==0,
           "t12g probe never reports an already-centered vertex");
    }
    /* t5e: the in-core fixed-window counter samples the actual emitted
     * triangle interiors and preserves fitted/fill identity. */
    {
        memset(vol,40,(size_t)C*C*C);
        float v[12]={8,8,8, 8,8,9, 8,9,8, 8,9,9};
        float n[12]={1,0,0, 1,0,0, 1,0,0, 1,0,0};
        float uv[8]={0,0, 2,0, 0,1, 2,1};
        uint8_t fill[4]={1,1,0,0};
        uint8_t dark_seed[4]={0};
        uint8_t dark_cell[1]={0};
        SnapBakeCounts bc=structured_bake_scan(&ct,v,n,fill,uv,2,2,47,193,13,
                                                dark_seed,NULL,dark_cell);
        CK(bc.filled_samples==2&&bc.filled_dark==2&&bc.fitted_samples==0&&
           bc.missing==0,
           "t5e fixed-window bake audit counts dark filled triangle samples");
        CK(dark_seed[0]&&dark_seed[1]&&dark_seed[2]&&dark_seed[3],
           "t5e dark triangle interiors seed their fitted corners for rescan");
        CK(!dark_cell[0],
           "t5e partly fitted seam is never a movable dark-patch core");
        memset(fill,1,sizeof fill);dark_cell[0]=0;
        (void)structured_bake_scan(&ct,v,n,fill,uv,2,2,47,193,13,
                                   NULL,NULL,dark_cell);
        CK(dark_cell[0],
           "t5e all-fill 75%-dark cell becomes a conservative patch core");
        memset(fill,0,sizeof fill);
        bc=structured_bake_counts(&ct,v,n,fill,uv,2,2,47,193,13);
        CK(bc.fitted_samples==2&&bc.fitted_dark==2&&bc.filled_samples==0,
           "t5e fixed-window bake audit separates fitted samples");
    }
    /* t5h: a coherent trusted rim already identifies the winding, so its ray
     * samples select a side by strict majority and the Laplacian absorbs local
     * abstentions/noise.  Only an ambiguous rim requires unanimous global
     * agreement; opposite signs there mean a wobbly minimal patch. */
    {
        int wobble=0;
        CK(patch_global_side(7,0,3,10,1,&wobble)==1&&!wobble,
           "t5h one trusted rim permits same-side rays plus Laplace-filled abstentions");
        CK(patch_global_side(7,0,3,10,0,&wobble)==0&&!wobble,
           "t5h ambiguous rim requires every unknown vertex to vote");
        CK(patch_global_side(7,1,2,10,1,&wobble)==1&&!wobble,
           "t5h coherent rim tolerates minority ray noise");
        CK(patch_global_side(5,5,0,10,1,&wobble)==0&&!wobble,
           "t5h coherent rim still abstains on a tied side vote");
        CK(patch_global_side(7,1,2,10,0,&wobble)==0&&wobble,
           "t5h ambiguous opposite votes classify the minimal patch as wobbly");
        CK(patch_global_side(10,0,0,10,0,&wobble)==1&&!wobble,
           "t5h ambiguous rim accepts unanimous whole-patch agreement");
    }
    /* t6c: a raw target on the wrong ridge family is a physical-vector
     * outlier even if a scalar depth representation could hide the direction
     * change.  Robust quilting must lower only that datum's strength. */
    {
        enum { H=3,W=3,N=9 };
        QuadAdj qa={0};qa.structured=1;qa.H=H;qa.W=W;
        uint8_t cand0[N],anchor0[N];float ct0[N],g0[N*3],w0[N];
        memset(cand0,1,sizeof cand0);memset(anchor0,0,sizeof anchor0);
        memset(g0,0,sizeof g0);
        for(int i=0;i<N;i++){ct0[i]=1.0f;g0[i*3]=1.0f;}
        ct0[4]=8.0f;
        int down=0;double mean=0.0;
        quilt_candidate_weights(cand0,ct0,g0,anchor0,&qa,N,1.5,w0,&down,&mean);
        CK(down==1&&w0[4]<0.30f&&w0[0]>0.99f,
           "t6c robust vector quilt downweights an isolated ridge-family switch");
        CK(mean>0.90&&mean<1.0,
           "t6c coherent target population retains nearly full data weight");
    }
    /* t6d: compact triangle meshes can have more than four neighbours.  The
     * robust-vector scratch space must follow the CSR degree rather than the
     * four-neighbour assumption of the structured lattice. */
    {
        enum { N=7 };
        const int32_t off[N+1]={0,6,7,8,9,10,11,12};
        const int32_t tgt[12]={1,2,3,4,5,6, 0,0,0,0,0,0};
        QuadAdj qa={0};qa.off=off;qa.tgt=tgt;
        uint8_t cand0[N],anchor0[N];float ct0[N],g0[N*3],w0[N];
        memset(cand0,1,sizeof cand0);memset(anchor0,0,sizeof anchor0);
        memset(g0,0,sizeof g0);
        for(int i=0;i<N;i++){ct0[i]=1.0f;g0[i*3]=1.0f;}
        ct0[0]=8.0f;
        int down=0;double mean=0.0;
        quilt_candidate_weights(cand0,ct0,g0,anchor0,&qa,N,1.5,w0,&down,&mean);
        CK(qadj_degree(&qa,0)==6,
           "t6d compact CSR exposes a degree-six triangle-mesh vertex");
        CK(down==1&&w0[0]<0.30f&&w0[1]>0.99f&&isfinite(mean),
           "t6d robust vector quilt handles CSR degree greater than four");
    }
    /* t6: the structured energy has exactly four axial neighbours and every
     * edge crosses red/black parity (triangle diagonals must never leak in). */
    {
        QuadAdj qa; memset(&qa, 0, sizeof qa); qa.structured = 1; qa.H = 3; qa.W = 4;
        CK(qadj_degree(&qa, 0) == 2, "t6 structured corner degree = 2");
        CK(qadj_degree(&qa, 5) == 4, "t6 structured interior degree = 4");
        int axial = 1, bipartite = 1;
        for (size_t v = 0; v < 12; v++) {
            int deg = qadj_degree(&qa, v);
            for (int e = 0; e < deg; e++) {
                size_t j = qadj_neighbor(&qa, v, e);
                size_t d = v > j ? v - j : j - v;
                if (d != 1 && d != 4) axial = 0;
                if (qadj_color(&qa, v) == qadj_color(&qa, j)) bipartite = 0;
            }
        }
        CK(axial, "t6 structured adjacency has no triangle diagonals");
        CK(bipartite, "t6 every structured edge crosses red/black parity");
    }
    /* t6b: quilting diagnostics are a disjoint, count-weighted partition. */
    {
        QuiltEdgeAccum ea = {{0}, {0}};
        QuadStripSnapRegionEdgeStats es;
        quilt_edge_accum(&ea, quilt_edge_class(0, 0), 3.0);
        quilt_edge_accum(&ea, quilt_edge_class(1, 1), 4.0);
        quilt_edge_accum(&ea, quilt_edge_class(0, 1), 12.0);
        quilt_edge_finish(&ea, &es);
        CK(es.all_edges == 3 && es.fitted_fitted_edges == 1 &&
           es.filled_filled_edges == 1 && es.fitted_filled_edges == 1,
           "t6b quilt edge classes are disjoint and exhaustive");
        CK(fabs(es.all_rms - sqrt(169.0/3.0)) < 1e-12 &&
           fabs(es.fitted_fitted_rms - 3.0) < 1e-12 &&
           fabs(es.filled_filled_rms - 4.0) < 1e-12 &&
           fabs(es.fitted_filled_rms - 12.0) < 1e-12,
           "t6b quilt RMS union is count weighted");
    }
    /* t7: accumulated tangential motion, not a harmless normal offset, trips
     * the event-driven remesh gate on an otherwise isometric grid. */
    {
        float ov[12*3], vv[12*3], uv[12*2], bn[12*3];
        for (int r=0;r<3;r++) for (int c=0;c<4;c++) {
            int i=r*4+c;
            ov[i*3+0]=vv[i*3+0]=(float)c;
            ov[i*3+1]=vv[i*3+1]=(float)r;
            ov[i*3+2]=vv[i*3+2]=0.0f;
            uv[i*2+0]=(float)c; uv[i*2+1]=(float)r;
            bn[i*3+0]=0.0f; bn[i*3+1]=0.0f; bn[i*3+2]=1.0f;
        }
        RemeshAudit ra0=remesh_audit(vv,ov,bn,uv,NULL,0,12,3,4,0.75,4.0,2.0,0.0001);
        CK(!ra0.advised && fabs(ra0.sander_p95-1.0)<0.02,
           "t7 fresh isometric lattice does not request remesh");
        for(int i=0;i<12;i++) vv[i*3+0]+=1.0f;
        RemeshAudit ra1=remesh_audit(vv,ov,bn,uv,NULL,0,12,3,4,0.75,4.0,2.0,0.0001);
        CK(ra1.advised && ra1.tangential_p95>0.9,
           "t7 one-cell tangential drift requests remesh");
    }
    /* t7b: a sparse severe tail must not hide below a healthy p95.  Pure
     * normal motion keeps the motion gates quiet, so this specifically proves
     * the >=4x population gate added for localized damage. */
    {
        enum { TH=20, TW=20, TN=TH*TW };
        float ov[TN*3], vv[TN*3], uv[TN*2], bn[TN*3];
        for (int r=0;r<TH;r++) for (int c=0;c<TW;c++) {
            int i=r*TW+c;
            ov[i*3]=vv[i*3]=(float)c;
            ov[i*3+1]=vv[i*3+1]=(float)r;
            ov[i*3+2]=vv[i*3+2]=0.0f;
            uv[i*2]=(float)c;uv[i*2+1]=(float)r;
            bn[i*3]=0.0f;bn[i*3+1]=0.0f;bn[i*3+2]=1.0f;
        }
        vv[((TH/2)*TW+TW/2)*3+2]=20.0f;
        RemeshAudit ra=remesh_audit(vv,ov,bn,uv,NULL,0,TN,TH,TW,
                                    100.0,100.0,100.0,0.005);
        CK(ra.advised && ra.over4>0 && ra.sander_p95<2.0,
           "t7 sparse >=4x distortion tail requests remesh below p95 gate");
    }
    /* t7c: a checkerboard/skew mode cancels exactly in the old averaged-quad
     * Jacobian even though one emitted triangle is severely compressed.  The
     * primary bake-faithful statistic must see it and fire the tail gate. */
    {
        float ov[4*3],vv[4*3],uv[4*2],bn[4*3];
        for(int r=0;r<2;r++)for(int c=0;c<2;c++) {
            int i=r*2+c;
            ov[i*3]=vv[i*3]=(float)c;
            ov[i*3+1]=vv[i*3+1]=(float)r;
            ov[i*3+2]=vv[i*3+2]=0.0f;
            vv[i*3]+=((r+c)&1)?-0.45f:0.45f;
            uv[i*2]=(float)c;uv[i*2+1]=(float)r;
            bn[i*3]=0.0f;bn[i*3+1]=0.0f;bn[i*3+2]=1.0f;
        }
        RemeshAudit ra=remesh_audit(vv,ov,bn,uv,NULL,0,4,2,2,
                                    100.0,100.0,100.0,0.25);
        CK(ra.quads==1&&ra.triangles==2&&ra.quad_sander_p95<1.02,
           "t7c averaged quad hides the alternating skew");
        CK(ra.advised&&ra.over4>0&&ra.sander_p95>4.0,
           "t7c bake-triangle metric exposes severe half-quad distortion");
    }
    /* t8: the objective-selected structured correspondence step equalizes a
     * badly sampled surface curve while preserving its geometric boundary and
     * transferring fill identity.  Its exact emitted-triangle tail must be
     * Pareto-nonworse, not merely its mean edge spacing. */
    {
        float v[8*3],uv[8*2]; uint8_t fl[8]; QuadStripRemeshStats rs;
        const float xx[4]={0.0f,0.1f,2.9f,3.0f};
        for(int r=0;r<2;r++)for(int c=0;c<4;c++){
            int i=r*4+c;v[i*3]=xx[c];v[i*3+1]=(float)r;v[i*3+2]=0.0f;
            uv[i*2]=(float)c;uv[i*2+1]=(float)r;
            fl[i]=(uint8_t)(c>=2);
        }
        int rc=QuadStripSnap_remesh_structured(v,fl,uv,2,4,1,&rs);
        double e0=v[3]-v[0],e1=v[6]-v[3],e2=v[9]-v[6];
        CK(rc==0 && rs.applied_sweeps==1 &&
           fabs(v[0])<1e-6 && fabs(v[9]-3.0)<1e-6,
           "t8 remesh preserves row endpoints");
        CK(e0>0.1 && e1<2.8 && e2>0.1,
           "t8 bounded remesh improves physical row sampling");
        CK(rs.after_over2<=rs.before_over2 && rs.after_over4<=rs.before_over4 &&
           rs.after_p95<=rs.before_p95+1e-12 &&
           rs.after_p99<=rs.before_p99+1e-12,
           "t8 accepted correspondence is triangle-tail Pareto-nonworse");
        CK(fl[0]==0 && fl[3]==1,
           "t8 remesh transfers fitted/fill identity to new samples");
    }
    /* t8b: an already uniform lattice is an explicit no-op.  This guards
     * against gratuitous resampling when no correspondence objective improves. */
    {
        float v[12*3],before[12*3],uv[12*2];uint8_t fl[12];
        QuadStripRemeshStats rs;
        for(int r=0;r<3;r++)for(int c=0;c<4;c++){
            int i=r*4+c;v[i*3]=(float)c;v[i*3+1]=(float)r;v[i*3+2]=0.0f;
            uv[i*2]=(float)c;uv[i*2+1]=(float)r;fl[i]=(uint8_t)((r+c)&1);
        }
        memcpy(before,v,sizeof v);
        int rc=QuadStripSnap_remesh_structured(v,fl,uv,3,4,2,&rs);
        CK(rc==0&&rs.applied_sweeps==0&&memcmp(v,before,sizeof v)==0,
           "t8b uniform lattice stays byte-identical when search has no gain");
    }
    /* t9: the structured line-search barrier accepts an ordinary update and
     * rejects a proposal that would reverse the quad Jacobian. */
    {
        float v[12]={0,0,0, 1,0,0, 0,1,0, 1,1,0};
        float o[12];memcpy(o,v,sizeof v);
        double good[3]={0.1,0.1,0.0},bad[3]={2.0,2.0,0.0};
        CK(structured_update_safe(v,o,2,2,0,good,0.05,0.0),
           "t9 orientation barrier accepts a safe update");
        CK(!structured_update_safe(v,o,2,2,0,bad,0.05,0.0),
           "t9 orientation barrier rejects a flipped quad");
        double far[3]={30.0,40.0,0.0};
        clamp_position_from_orig(far,o,0,12.0);
        CK(fabs(sqrt(far[0]*far[0]+far[1]*far[1])-12.0)<1e-9,
           "t9 position proposal is capped at guided physical reach");
    }

    /* t10: a low-frequency vector quilt mode is the failure case that a few
     * fine red/black sweeps cannot remove.  With equal fine-sweep work (three
     * baseline versus two pre + one post), the aggregate V-cycle must reduce
     * the true screened-Laplacian residual substantially more. */
    {
        enum { GH=65, GW=65, GN=GH*GW };
        float *db=(float *)calloc((size_t)GN*3u,sizeof(float));
        float *cb=(float *)calloc(GN,sizeof(float));
        float *dm=(float *)calloc((size_t)GN*3u,sizeof(float));
        float *cm=(float *)calloc(GN,sizeof(float));
        float *qt=(float *)calloc(GN,sizeof(float));
        float *qg=(float *)calloc((size_t)GN*3u,sizeof(float));
        float *qw=(float *)calloc(GN,sizeof(float));
        uint8_t *fr=(uint8_t *)calloc(GN,1),*fi=(uint8_t *)calloc(GN,1);
        uint8_t *ca=(uint8_t *)calloc(GN,1);
        int alloc_ok=db&&cb&&dm&&cm&&qt&&qg&&qw&&fr&&fi&&ca;
        CK(alloc_ok,"t10 vector multigrid fixture allocates");
        if(alloc_ok){
            for(int r=0;r<GH;r++)for(int c=0;c<GW;c++){
                int i=r*GW+c;fr[i]=(uint8_t)(r==0||c==0||r+1==GH||c+1==GW);
                if(!fr[i]&&r>=GH/2-2&&r<=GH/2+2&&c>=12&&c<GW-12){
                    ca[i]=1;qw[i]=1.0f;qg[(size_t)i*3]=1.0f;
                    qt[i]=(float)(0.75+0.25*sin(0.08*(double)c));
                    db[(size_t)i*3]=dm[(size_t)i*3]=qt[i];cb[i]=cm[i]=1.0f;
                }
            }
            QuadStripSnapOpts qb,qm;QuadStripSnap_defaults(&qb);qm=qb;
            qb.quilt_anchor=0.001;qb.quilt_sweeps=3;qb.multigrid_levels=1;
            qm.quilt_anchor=0.001;qm.quilt_sweeps=2;qm.multigrid_levels=8;
            qm.multigrid_cycles=1;
            QMGRunStats sb,sm;int rb=qmg_solve_quilt(db,cb,fr,fi,qw,qt,qg,GH,GW,&qb,&sb);
            QMGHierarchy check;int hc=qmg_build_quilt(&check,fr,fi,qw,GH,GW,
                qb.quilt_smooth,qb.quilt_seam_smooth,qb.quilt_anchor,8);
            double base_rel=hc==0?qmg_quilt_restrict_residual(&check,db,cb,fr,fi,
                qw,qt,qg,GH,GW,qb.quilt_smooth,qb.quilt_seam_smooth,qb.quilt_anchor):HUGE_VAL;
            if(hc==0)qmg_hierarchy_free(&check);
            int rm=qmg_solve_quilt(dm,cm,fr,fi,qw,qt,qg,GH,GW,&qm,&sm);
            CK(rb==0&&rm==0&&hc==0&&sm.levels_used>=4&&sm.corrections==1,
               "t10 vector hierarchy executes a real coarse correction");
            CK(isfinite(base_rel)&&sm.residual_after<0.70*base_rel,
               "t10 vector V-cycle beats equal fine-grid smoothing on low modes");
            int bounded=1;for(int i=0;i<GN;i++){
                double d2=0.0;for(int k=0;k<3;k++)d2+=(double)dm[(size_t)i*3+k]*dm[(size_t)i*3+k];
                if(cm[i]<-1e-6f||cm[i]>1.000001f||sqrt(d2)>qm.guided_reach+1e-6)bounded=0;
            }
            CK(bounded,"t10 vector correction preserves offset/confidence bounds");

            memset(dm,0,(size_t)GN*3u*sizeof *dm);memset(cm,0,GN*sizeof *cm);
            memset(ca,0,GN);memset(qw,0,GN*sizeof *qw);memset(qt,0,GN*sizeof *qt);
            QMGRunStats sz;int rz=qmg_solve_quilt(dm,cm,fr,fi,qw,qt,qg,GH,GW,&qm,&sz);
            int zero=1;for(int i=0;i<GN;i++)if(dm[(size_t)i*3]!=0.0f||
                dm[(size_t)i*3+1]!=0.0f||dm[(size_t)i*3+2]!=0.0f||cm[i]!=0.0f)zero=0;
            CK(rz==0&&zero&&sz.corrections==0,
               "t10 zero-residual vector hierarchy is an exact no-op");
        }
        free(db);free(cb);free(dm);free(cm);free(qt);free(qg);free(qw);
        free(fr);free(fi);free(ca);
    }

    /* t11: the same hierarchy corrects the fixed-rotation 3-D position system.
     * A smooth z displacement is deliberately under-converged by fine sweeps;
     * the V-cycle must lower residual while the global emitted-triangle barrier
     * retains orientation. */
    {
        enum { PH=33, PW=33, PN=PH*PW };
        float *orig=(float *)malloc(PN*3*sizeof(float));
        float *vb=(float *)malloc(PN*3*sizeof(float));
        float *vm=(float *)malloc(PN*3*sizeof(float));
        float *tar=(float *)malloc(PN*3*sizeof(float));
        float *puv=(float *)malloc(PN*2*sizeof(float));
        float *pcf=(float *)calloc(PN,sizeof(float));
        double *pr=(double *)calloc(PN*9,sizeof(double));
        uint8_t *pfr=(uint8_t *)calloc(PN,1),*pan=(uint8_t *)calloc(PN,1);
        uint8_t *pha=(uint8_t *)calloc(PN,1);
        int alloc_ok=orig&&vb&&vm&&tar&&puv&&pcf&&pr&&pfr&&pan&&pha;
        CK(alloc_ok,"t11 position multigrid fixture allocates");
        if(alloc_ok){
            for(int r=0;r<PH;r++)for(int c=0;c<PW;c++){
                int i=r*PW+c;
                orig[i*3]=(float)c;orig[i*3+1]=(float)r;orig[i*3+2]=0.0f;
                memcpy(&vb[i*3],&orig[i*3],3*sizeof(float));
                memcpy(&vm[i*3],&orig[i*3],3*sizeof(float));
                memcpy(&tar[i*3],&orig[i*3],3*sizeof(float));
                puv[i*2]=(float)c;puv[i*2+1]=(float)r;
                pr[i*9]=pr[i*9+4]=pr[i*9+8]=1.0;
                pfr[i]=(uint8_t)(r==0||c==0||r+1==PH||c+1==PW);
                if(!pfr[i]&&r>=PH/2-1&&r<=PH/2+1&&c>=PW/2-3&&c<=PW/2+3){
                    pha[i]=1;pcf[i]=1.0f;tar[i*3+2]=1.5f;
                }
            }
            QuadAdj pa;memset(&pa,0,sizeof pa);pa.structured=1;pa.H=PH;pa.W=PW;
            QuadStripSnapOpts pb,pm;QuadStripSnap_defaults(&pb);pm=pb;
            pb.w_orig=0.001;pb.w_target=1.0;pb.w_smooth=0.5;pb.w_dev=0.0;
            pb.solve_sweeps=3;pb.solve_omega=0.8;pb.multigrid_levels=1;
            pm=pb;pm.solve_sweeps=2;pm.multigrid_levels=7;pm.multigrid_cycles=1;
            pm.multigrid_patch=4;
            int btb=0,brb=0,btm=0,brm=0;QMGRunStats sb,sm;
            int rb=qmg_solve_position(vb,orig,tar,puv,pr,pfr,pan,pha,pcf,NULL,&pa,&pb,
                                      &btb,&brb,&sb);
            QMGHierarchy check;int hc=qmg_build_position(&check,pfr,pan,pha,pcf,
                                                          PH,PW,&pm);
            double base_rel=hc==0?qmg_position_restrict_residual(&check,vb,orig,
                tar,puv,pr,pfr,pan,pha,pcf,&pa,&pm):HUGE_VAL;
            if(hc==0)qmg_hierarchy_free(&check);
            int rm=qmg_solve_position(vm,orig,tar,puv,pr,pfr,pan,pha,pcf,NULL,&pa,&pm,
                                      &btm,&brm,&sm);
            RemeshQuality pq=remesh_triangle_quality(vm,orig,puv,PH,PW);
            CK(rb==0&&rm==0&&hc==0&&sm.levels_used>=4&&sm.corrections==1&&
               sm.patches_accepted>1,
               "t11 position hierarchy executes localized barriered corrections");
            CK(isfinite(base_rel)&&sm.residual_after<0.75*base_rel,
               "t11 position V-cycle beats equal fine-grid smoothing on low modes");
            CK(pq.flipped==0&&pq.degenerate==0,
               "t11 position coarse correction preserves emitted triangles");
        }
        free(orig);free(vb);free(vm);free(tar);free(puv);free(pcf);free(pr);
        free(pfr);free(pan);free(pha);
    }

    /* t13: end-to-end connected-patch proposal.  A 3x3 unknown island sits in
     * dark resin at x=8, with its one connected trusted collar already on the
     * x=4..5 sheet.  Bidirectional rays also see a farther x=13..14 sheet; the
     * whole patch must select the nearer inward side, snap toward it, and decay
     * displacement smoothly to the fixed collar without axial drift. */
    {
        enum { H=7,W=7,N=H*W };
        memset(vol,40,(size_t)C*C*C);
        for(long z=0;z<C;z++)for(long y=0;y<C;y++)for(long x=0;x<C;x++)
            if((x>=4&&x<=5)||(x>=13&&x<=14))
                vol[((size_t)z*C+(size_t)y)*C+(size_t)x]=210;
        float pv[N*3],pn[N*3],puv[N*2],proposal[N*3];
        uint8_t pf[N],moved[N];int32_t moved_patch[N];
        for(int r=0;r<H;r++)for(int c=0;c<W;c++) {
            int v=r*W+c;int unknown=r>=2&&r<=4&&c>=2&&c<=4;
            pv[v*3]=(float)(4+r);pv[v*3+1]=(float)(4+c);
            pv[v*3+2]=unknown?8.0f:5.0f;
            pn[v*3]=pn[v*3+1]=0.0f;pn[v*3+2]=1.0f;
            puv[v*2]=(float)(2*c);puv[v*2+1]=(float)r;
            pf[v]=(uint8_t)unknown;
        }
        QuadStripSnapOpts po;QuadStripSnap_defaults(&po);
        po.snap_side=QUAD_SNAP_SIDE_MIDLINE;po.preserve_axial=1;
        po.axis_y=8.0;po.axis_x=0.0;po.patch_buffer=2;
        po.patch_min_cells=1;po.patch_ray_reach=8.0;
        po.patch_ray_margin=1.0;po.patch_relax_sweeps=80;
        PatchRecoverRunStats ps={0};
        int rc=patch_recover_propose(&ct,&po,pv,pn,pf,puv,H,W,130.0,
                                     proposal,moved,moved_patch,&ps);
        int known_fixed=1,axial_fixed=1;
        for(int v=0;v<N;v++) {
            if(!pf[v]&&moved[v])known_fixed=0;
            if(fabs((double)proposal[v*3]-(double)pv[v*3])>1e-7)axial_fixed=0;
        }
        int mid=3*W+3,corner=2*W+2;
        double mid_move=(double)pv[mid*3+2]-(double)proposal[mid*3+2];
        double corner_move=(double)pv[corner*3+2]-(double)proposal[corner*3+2];
        CK(rc==0&&ps.regions==1&&ps.border_consensus==1&&
           ps.ray_consensus==1&&ps.proposed==1,
           "t13 connected dark patch earns one trusted-rim/global-ray proposal");
        CK(moved[mid]&&mid_move>2.0&&proposal[mid*3+2]>4.0f,
           "t13 agreed patch snaps substantially toward the nearer winding");
        CK(known_fixed&&axial_fixed&&mid_move>corner_move,
           "t13 local Laplacian blends snapped unknowns into a fixed known collar");

        enum { NF=2*(H-1)*(W-1) };
        int32_t gf[NF*3];size_t fn=0;
        for(int r=0;r+1<H;r++)for(int c=0;c+1<W;c++) {
            int32_t a=(int32_t)(r*W+c),b=a+1,cc=a+W,d=cc+1;
            gf[fn*3]=a;gf[fn*3+1]=b;gf[fn*3+2]=d;fn++;
            gf[fn*3]=a;gf[fn*3+1]=d;gf[fn*3+2]=cc;fn++;
        }
        CSR_T gcsr=CSR_from_faces(arena,gf,fn,N);
        QuadAdj ga={0};ga.off=CSR_offset(gcsr);ga.tgt=CSR_target(gcsr);
        PatchRecoverRunStats gps={0};
        rc=patch_recover_graph_propose(&ct,&po,pv,pn,pf,puv,gf,fn,&ga,N,
            130.0,proposal,moved,moved_patch,&gps);
        mid_move=(double)pv[mid*3+2]-(double)proposal[mid*3+2];
        known_fixed=1;axial_fixed=1;
        for(int v=0;v<N;v++) {
            if(!pf[v]&&moved[v])known_fixed=0;
            if(fabs((double)proposal[v*3]-(double)pv[v*3])>1e-7)axial_fixed=0;
        }
        CK(rc==0&&gps.regions==1&&gps.border_consensus==1&&
           gps.ray_consensus==1&&gps.proposed==1&&moved[mid]&&mid_move>2.0,
           "t13b compact face-graph patch makes one global winding decision");
        CK(known_fixed&&axial_fixed,
           "t13b compact Laplacian keeps trusted and axial coordinates fixed");
    }

    Arena_dispose(&arena);
    fprintf(stderr, "[selftest] %s (%d failures)\n", fails == 0 ? "ALL PASS" : "FAILURES", fails);
    return fails;
}
