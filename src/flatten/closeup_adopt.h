/* ============================================================================
 * closeup_adopt.h -- the reduced placement of an ink close-up: pool charts
 * (other pieces of the assembly, kibble the minimum-piece filter dropped)
 * that continue the fused sheet are joined to it by their 3-D contact.
 *
 * A pool chart touches the fused surface (closeup_surface.h) where its
 * vertices lie within on_vox of it along the normal, with their foot on data
 * (or on a fill node within foot_fill_px of data: the trim gap between cube
 * charts is ~1.25 vox) and a normal within normal_deg.  Those contacts are
 * correspondences between the chart's own flattening (its layout UV, any
 * rigid frame; reflected first when it shows the recto mirror-imaged) and the
 * frame UV of their feet; a 2-D rigid fit (Procrustes) places the chart.  It
 * is ADOPTED when the fit leaves at most rms_px, the contacts span at least
 * min_extent_px, at least agree_frac of its vertices that land on the surface
 * lie within conflict_vox of it (another wrap runs 12-14 vox away: a chart of
 * it placed through a touch point lands on the surface far off), and at least
 * new_frac land where the surface has no data (it adds papyrus). Conflict
 * validation uses the same supported band as contacts: data and fill within
 * foot_fill_px of data. More distant extrapolation is not an observation and
 * cannot veto a measured continuation.  Where it adds papyrus it must also
 * continue the surface around it: each vertex landing on no supported node
 * lies within continuity_slope x its frame distance + conflict_vox of the
 * nearest supported node within continuity_px (and, among the charts accepted
 * in one round, of the nearest earlier claim).  The chart keeps its own flattening,
 * so it must also be metrically sound: faces whose flattening stretches or
 * squeezes beyond metric_ratio may carry at most metric_frac of its 3-D area
 * (2026-10-07: squeezed pool ribbons rendered as a smeared bridge).  And it must
 * not stack onto the surface: at most layer_frac of its vertices may have
 * surface data within layer_vox across their normal that lies far away along
 * the surface (2026-10-15: two routes laid charts on neighbouring sheets) --
 * nor onto a chart accepted earlier in the same round (that one waits).  Charts
 * accepted in one round claim their new nodes in order of support; a later one
 * whose claim disagrees in 3-D waits for the next round.
 * ==========================================================================*/
#ifndef CLOSEUP_ADOPT_H
#define CLOSEUP_ADOPT_H

#include <stddef.h>
#include <stdint.h>

#include "closeup_surface.h"

/* The pool: a UV mesh whose vertices and faces are grouped by chart. */
typedef struct CloseupPool {
    size_t         nv, nf, nch;
    const float   *verts;     /* [nv*3] world (z, y, x) */
    const float   *nrm;       /* [nv*3] unit vertex normals (any sign) */
    const float   *inward;    /* [nv*3] unit direction toward the axis */
    const double  *uv;        /* [nv*2] the chart's own flattening */
    const int32_t *faces;     /* [nf*3] */
    const size_t  *vofs;      /* [nch+1] chart k owns vertices vofs[k] .. vofs[k+1]-1 */
    const size_t  *fofs;      /* [nch+1] ... and faces fofs[k] .. fofs[k+1]-1 */
} CloseupPool;

typedef struct CloseupAdoptOpts {
    double on_vox, foot_fill_px, min_extent_px, rms_px, conflict_vox, agree_frac, normal_deg, new_frac;
    int    min_corr;
    double metric_ratio, metric_frac;   /* the metric gate (CLOSEUP_ADOPT_METRIC_*) */
    double continuity_px, continuity_slope;   /* the continuity gate (CLOSEUP_ADOPT_CONTINUITY_*); px 0: off */
    double layer_vox, layer_frac, layer_margin_px;   /* the layer gate (CLOSEUP_ADOPT_LAYER_*); vox 0: off */
} CloseupAdoptOpts;

void CloseupAdoptOpts_default(CloseupAdoptOpts *o);

/* frame uv = R(theta) (reflect ? -u : u, v) + (tx, ty), R = [c -s; s c] */
typedef struct CloseupPlacement {
    size_t chart;            /* pool chart index */
    int    reflect;
    double c, s, tx, ty;
    size_t corr, landed, agree, fresh;   /* contacts; vertices on the surface; ...within conflict_vox; on no data */
    double rms, extent;
} CloseupPlacement;

typedef struct CloseupAdoptStats {
    size_t candidates, touched, placed;
    size_t refused_extent, refused_rms, refused_conflict, refused_nothing_new, deferred;
    size_t refused_metric, refused_continuity;
    size_t refused_layer;
    size_t deferred_layer;   /* of deferred: stacked onto a chart accepted earlier in the same round */
} CloseupAdoptStats;

/* One round against the current surface: every pool chart with done[k] == 0 is tried; *out (malloc'd) lists
 * the adopted charts in acceptance order.  0 on success, -1 on no memory. */
int CloseupAdopt_round(const CloseupGrid *g, const CloseupPool *pool, const uint8_t *done,
                       const CloseupAdoptOpts *o, CloseupPlacement **out, size_t *nout, CloseupAdoptStats *st);

/* Read-only diagnostics, delivered serially after placement.  Reason: 0 no sufficient contact,
 * 1 extent, 2 RMS, 3 conflict, 4 nothing new, 5 placed, 6 deferred by another chart, 7 metric (its own
 * flattening stretches beyond metric_ratio over more than metric_frac of its area: metric_bad_share), 8 continuity
 * (where it adds papyrus, more than 1 - agree_frac of its vertices lie farther in 3-D from the nearest supported
 * node within continuity_px than continuity_slope x their frame distance + conflict_vox: cont_bad of cont_checked),
 * 9 layered (more than layer_frac of its sampled vertices have surface data within layer_vox across their normal, far
 * along the surface: layer_bad of layer_checked).
 * Split conflict counts are populated only when the contact and RMS gates passed.
 * The observer must not retain the borrowed record or modify the surface/pool; nonzero aborts the call. */
typedef struct CloseupAdoptAudit {
    CloseupPlacement p;
    size_t screen_hits, data_landed, data_agree, fill_landed, fill_agree;
    size_t ignored_fill_landed, ignored_fill_agree;
    double data_max_error, fill_max_error;
    double ignored_fill_max_error;
    int reason;
    double metric_bad_share;    /* 3-D area share of faces beyond metric_ratio (every chart that was evaluated) */
    size_t cont_checked, cont_bad;   /* the continuity gate's vertices (see reason 8) */
    size_t layer_checked, layer_bad;  /* the layer gate's vertices (see reason 9) */
} CloseupAdoptAudit;

typedef int (*CloseupAdoptObserver)(const CloseupAdoptAudit *record, void *context);

int CloseupAdopt_round_audited(const CloseupGrid *g, const CloseupPool *pool, const uint8_t *done,
                              const CloseupAdoptOpts *o, CloseupPlacement **out, size_t *nout,
                              CloseupAdoptStats *st, CloseupAdoptObserver observer, void *context);

/* The placement applied to a chart-UV point. */
void CloseupPlacement_apply(const CloseupPlacement *p, double u, double v, double *fu, double *fv);

int CloseupAdopt_selftest(void);

#endif /* CLOSEUP_ADOPT_H */
