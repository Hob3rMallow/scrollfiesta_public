#ifndef ASM_AXIS_DERIVE_INCLUDED
#define ASM_AXIS_DERIVE_INCLUDED

#include <stddef.h>

#include "../common/arena.h"
#include "asm_axis.h"
#include "asm_types.h"

/* ============================================================================
 * asm_axis_derive.h -- the scroll axis DERIVED from the material (2026-09-17).
 *
 * The wrap normals of a rolled sheet all point through its centre.  Given
 * the cleaned charts (face centroids, face normals, face areas), the axis is
 * found slab by slab: samples are projected onto the current polyline,
 * binned by arc length s, expressed in the LOCAL tangent frame of their
 * slab (so a tilted scroll is not smeared by z-slabs), and the point closest
 * to the in-plane normal lines is fitted robustly (Cauchy IRLS, the scale
 * annealed from ASM_AXIS_DERIVE_CAUCHY_VOX0 down to ASM_AXIS_DERIVE_CAUCHY_VOX).
 * The accepted slab centres are smoothed through the same LOESS as a loaded
 * table (AsmAxis_from_rows) and the whole is iterated from a prior -- the
 * configured table, the umbilicus line, or a seed found on z-slabs from a
 * grid of starts (the void-qualified, lowest-residual attractor wins:
 * evolute cusps of flattened wraps have MORE support, never rank by it).
 *
 * Why not the standalone tracker: its lock rule needs the centre inside the
 * material (>= 1,500 samples within 250 vox, >= 9 of 12 sectors, slope
 * <= 0.5) and z-perpendicular slabs, which abstains by construction on a
 * column beside the umbilicus (the 21x5x5) and smears a tilted scroll's
 * sections (pherc343: 1 lock of 10 slabs).  Here every slab reports its
 * conditioning and aperture instead of abstaining on a sector count, and a
 * one-sided aperture is a documented state (the direction is robust, the
 * radius along the bisector is not), not a failure.
 * ==========================================================================*/

typedef struct AsmAxisDeriveOpts {
    double        slab_vox;    /* arc-length slab thickness (ASM_AXIS_DERIVE_SLAB_VOX) */
    int           iters;       /* projection / fit / smooth iterations (ASM_AXIS_DERIVE_ITERS) */
    size_t        per_chart;   /* face samples per chart (ASM_AXIS_DERIVE_SAMPLES_PER_CHART) */
    const AsmAxis *prior;      /* nullable: the table or line to start from; NULL = seed */
} AsmAxisDeriveOpts;

typedef struct AsmAxisDeriveRow {
    double s, z, y, x;         /* slab centre along the CURRENT axis (s) and the fitted centre (world) */
    size_t n;                  /* samples that fed the fit (in-plane normals, |n.t| <= the axial limit) */
    size_t kept;               /* ... with a final Cauchy weight above one half */
    double cond;               /* lambda_min / lambda_max of the weighted normal scatter */
    double aperture_deg;       /* angular span of the kept normal LINES (mod 180) */
    double resid_p50;          /* weighted median distance of the centre from the kept normal lines, vox */
    double unc_major, unc_minor;   /* resid_p50 x sqrt(sum w / lambda): the uncertainty ellipse's axes, vox */
    double void_ratio;         /* area density inside 50 vox of the centre over the annulus to 150 vox */
    int    sectors;            /* 30-degree position sectors around the centre holding kept samples (0..12) */
    int    one_sided;          /* the kept samples surround less than half of the centre */
    int    status;             /* 0 abstained, 1 locked */
} AsmAxisDeriveRow;

typedef struct AsmAxisDeriveReport {
    AsmAxisDeriveRow *rows; size_t n_rows, n_lock;
    int    seeded;             /* 1 = the prior was a seed found on z-slabs, 0 = the caller's prior, -1 = the world-z fallback */
    size_t seed_slabs;         /* z-slabs that yielded a void-qualified seed */
    int    iters;              /* iterations run */
    double last_move;          /* the largest slab-centre move of the last iteration, vox */
    double tilt_deg;           /* mean direction off +z */
    double fit_rms;            /* the smoothed curve against its rows, vox */
    double r_curv_min;
    double prior_rms;          /* rms distance of the derived curve from the prior over the rows (-1 without one) */
    size_t one_sided;
    double cond_min, aperture_min_deg, resid_p50_max;
    size_t samples;
    /* the independent direction measurement (AsmAxis_normal_direction) and the derived curve against it */
    double normal_dir[3], normal_cond, normal_frac, normal_tilt_deg, normal_disagree_deg;
    int    accepted;           /* the derived curve is usable: at least ASM_AXIS_DERIVE_ACCEPT_LOCK of the slabs locked and the
                                * rows sit within ASM_AXIS_DERIVE_ACCEPT_RMS of the smoothed curve (else the caller's frame stands) */
    double sec;
} AsmAxisDeriveReport;

void AsmAxisDerive_default_opts(AsmAxisDeriveOpts *o);

/* THE AXIS DIRECTION FROM THE NORMALS ALONE.  A rolled sheet is everywhere
 * parallel to its axis, so every surface normal is perpendicular to it: the
 * direction minimising sum_i w_i (n_i . t)^2 is the smallest eigenvector of
 * sum_i w_i n_i n_i^T.  It needs no centre, no circularity and no prior -- the
 * measurement that survives the folds and waves a normal-LINE fit chokes on.
 * *cond = lambda_min / lambda_mid (how well the null direction is separated
 * from the sheet's own spread; a flat pile of parallel sheets reads ~1 and the
 * direction means nothing), *frac = the weighted fraction of samples with
 * |n . t| <= 0.2.  The sign is made positive in z.  Returns 0 on success. */
int AsmAxis_normal_direction(const float *nrm, const float *w, size_t n, double dir[3], double *cond, double *frac);

/* Samples: p (z,y,x) [n*3], unit normals [n*3] (sign irrelevant), weights [n]
 * (face areas).  Returns the derived polyline, or NULL when no slab locked. */
AsmAxis *AsmAxis_derive(Arena_T arena, const float *p, const float *nrm, const float *w, size_t n,
                        const AsmAxisDeriveOpts *o, AsmAxisDeriveReport *rep);

/* Gathers the samples from the run's in-layout charts (face centroids, face
 * normals, face areas, at most o->per_chart faces per chart) and derives. */
AsmAxis *AsmAxis_derive_run(Arena_T arena, const AsmRun *run, const AsmAxisDeriveOpts *o, AsmAxisDeriveReport *rep);

/* The ledger: '#' provenance header, then z,y,x,status,s,n,kept,cond,aperture_deg,
 * resid_p50,unc_major,unc_minor,sectors,void_ratio per slab -- loadable through
 * AsmAxis_load (its first three columns).  Returns 0 on success. */
int AsmAxisDerive_write_csv(const AsmAxisDeriveReport *rep, const char *path);

/* A viewable bitmap: the derived y(z) and x(z) with the prior's, and the rows' residuals. */
int AsmAxisDerive_write_png(Arena_T arena, const AsmAxis *derived, const AsmAxis *prior, const AsmAxisDeriveReport *rep, const char *path);

int AsmAxisDerive_selftest(void);

#endif
