#ifndef ASM_RELATE_INCLUDED
#define ASM_RELATE_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "asm_types.h"

/* ============================================================================
 * asm_relate.h -- relations between charts.
 *
 * SEAM JOINS: boundary vertices of charts in adjacent cubes that are mutual
 * nearest neighbours across the trim gap (<= ASM_SEAM_MAX_VOX), with
 * agreeing normals (absolute value: per-cube orientation is independent, the
 * sign becomes the relation's PARITY), a tangential displacement (the
 * TANGENCY gate: a pair that meets along the normal is two wraps touching),
 * and enough matched length.  Each accepted chart pair gets one rigid 2-D
 * transform fitted on source-transported correspondences in both trim
 * directions. Pairwise gates compare physical midpoints, with one vote
 * per original observation. The fit is checked for isometry and
 * unimodality (a sheet switch inside a seam is bimodal), with information
 * weights from the residual.
 *
 * LAYER NEIGHBOURS: rays along +/- the vertex normal from a sample of each
 * chart's vertices; the first other-chart surface hit gives the local layer
 * distance and the layer-neighbour pair (negative evidence, and the ruler
 * used wherever a pitch constant would otherwise appear).
 * ==========================================================================*/

#define ASM_RELATE_HIST 16

typedef struct AsmRelateOpts {
    double cube_size;         /* pile cube edge, vox (128) */
    double seam_max_vox;      /* mutual-nearest pair gate */
    double seam_band_vox;     /* boundary vertices this close to the neighbour's box are candidates */
    double seam_min_len;      /* matched seam length gate, vox */
    double normal_dot_min;    /* |n_a . n_b| >= this */
    double tangency_max_vox;  /* per-pair median |d . n| <= this */
    double iso_tol;           /* isometry check tolerance (fraction) */
    double seam_noise_vox;    /* position noise budget of one correspondence (ASM_SEAM_NOISE_VOX) */
    double layer_probe_vox;   /* ray reach along the normal */
    int    layer_samples;     /* vertices sampled per chart for the probes */
    const char *diag_csv;     /* optional per-chart-pair gate dump */
} AsmRelateOpts;

typedef struct AsmRelateStats {
    size_t cube_pairs;
    size_t candidates;        /* mutual-nearest pairs passing distance + normal */
    size_t chart_pairs;       /* distinct chart pairs with >= 3 candidates */
    size_t rej_len, rej_tangency, rej_iso, rej_unimodal, rej_fit;
    size_t wrap_on_iso;       /* isometry-rejected pairs the unimodality gate ALSO rejects (ASM_WRAP_MEASURE_WEAK) */
    size_t wrap_anchors;      /* ASM_REL_CROSSWRAP relations still handed to placement as anchors (ASM_WRAP_VETO_ANCHOR 0) */
    size_t weak_kept, short_kept;            /* gate-rejected relations kept as placement evidence (WEAK / SHORT) */
    size_t bridge_demoted;                   /* load-bearing seams demoted by the bridge premium */
    size_t bridge_corroborated;              /* marginal bridges spared because another seam joins the same two blocks (ASM_SEAM_BRIDGE_CORROBORATED) */
    size_t bridge_leaves;                    /* marginal bridges spared because one block is below ASM_SEAM_BRIDGE_LEAF_AREA */
    size_t contact_ledgered;                 /* tangency rejections kept as ASM_REL_CONTACT ledger rows (ASM_RELATE_LEDGER_CONTACT) */
    size_t short_rej_corr, short_rej_rms;    /* short seams not kept: too few correspondences / rigid rms over the limit */
    size_t accepted;
    size_t corr_total;
    size_t parity_flipped;    /* accepted relations with opposite orientation */
    double rms_p50, rms_p95;
    double iso_disc_sum;       /* sum of per-pair median relative seam discrepancy */
    size_t iso_disc_n;
    size_t layer_hits;
    size_t layer_hist[ASM_RELATE_HIST];  /* hits per 2-vox bin of probe distance */
    size_t gap_hist[ASM_RELATE_HIST];    /* mutual-nearest seam pairs per 1-vox bin of 3-D distance (all, ungated) */
    size_t gap_pairs;
    size_t layer_pairs;
    size_t charts_with_layer;
    double layer_d_global;
    double seam_sec, layer_sec;
} AsmRelateStats;

void AsmRelate_default_opts(AsmRelateOpts *o);

/* Fills run->rels / run->corr / run->layers / run->chart_layer_d.  Charts
 * must be loaded (stage 1).  Parallel over cube pairs with `threads`. */
int AsmRelate_run(AsmRun *run, const AsmRelateOpts *opts, int threads,
                  AsmRelateStats *stats);

/* 2-D rigid fit b -> a on paired points (qa, qb are [n*2]), weighted by w
 * (NULL = 1).  Returns theta, t and the rms; 0 on success. */
int AsmRelate_rigid_fit(const double *qa, const double *qb, const double *w, size_t n,
                        double *theta, double *tx, double *ty, double *rms);

int AsmRelate_selftest(void);

#endif
