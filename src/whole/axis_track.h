#ifndef VESUVIUS_AXIS_TRACK_H
#define VESUVIUS_AXIS_TRACK_H

#include <stddef.h>
#include <stdint.h>

/* ============================================================================
 * axis_track.h -- derive the scroll's umbilicus CURVE from its wrap meshes.
 *
 * The umbilicus of PHerc0139 is not a line: measured 2026-09-03 it drifts
 * ~625 vox over 2,600 z (slope up to 0.24).  Every theta/r stage of the
 * flatten lane (certificate, fit, solid fills, verdict, core routing, tile
 * selector) needs the LOCAL centre c(z).  This module estimates c(z) per
 * z-slab from two independent observations of the wraps, so the curve is
 * derived, never copied from a manifest:
 *
 *   PRIMARY (orientation): the point closest to the in-plane wrap-normal
 *   lines, 60%-trimmed IRLS, restricted to samples within
 *   AXIS_TRACE_LOCAL_R of a running prior.  Wrap normals of a concentric
 *   spiral all point through its centre; the local window keeps the
 *   whole-slab centre-of-mass bias of a flattened scroll out.
 *
 *   SECONDARY (position): polar first-harmonic recentring.  Rasterize the
 *   slab in (theta, r) about the current centre, extract wrap ridges that
 *   span >= AXIS_TRACE_RIDGE_MIN_ARC of theta, fit
 *       r(theta) = R + k*theta + a*cos(theta) + b*sin(theta)
 *   per ridge and move the centre by the mean (a, b) until it stops moving.
 *   hypot(a, b) is the eccentricity of the wraps about the centre (d1);
 *   the signed mean slope k is the scroll's winding pitch per radian and
 *   its SIGN is the physical winding sense, measured without any lift.
 *
 * MEASURED 2026-09-02 on the PHerc0139 4x5x5 pile: the wraps are nowhere
 * near circular (V-folds and waves of +-50..100 vox in the polar raster),
 * so the harmonic model cannot chain a ridge through a fold and locks on
 * few slabs.  The PRIMARY estimator is therefore the position authority
 * (it is the estimator the eyes validated on this scroll, residual p50
 * ~22 vox); the secondary is an optional consistency check reported when
 * it locks.  A slab where the primary abstains is left to interpolation,
 * never held.
 * ==========================================================================*/

typedef enum {
    AXIS_ROW_ABSTAIN = 0,   /* no estimate: interpolated or dropped */
    AXIS_ROW_LOCK    = 1,   /* primary normal-line position locked */
    AXIS_ROW_INTERP  = 2,   /* linear interpolation across a short gap */
    AXIS_ROW_EXTRAP  = 3    /* explicit extrapolation past the last lock */
} AxisRowStatus;

/* One z-slab's sample reservoir.  Filled by the streaming driver one cube
 * at a time (reservoir sampling keeps the memory bounded), consumed by the
 * estimators.  Face samples carry the face centroid (y, x) and the unit
 * in-plane normal (ny, nx); vertex samples carry (y, x). */
typedef struct {
    double z0, z1;
    double *fy, *fx, *fny, *fnx, *fa;   /* fa = triangle area, vox^2 */
    size_t  nf, nf_cap, nf_seen;
    double *vy, *vx;
    size_t  nv, nv_cap, nv_seen;
    uint32_t rng;
} AxisSlab;

/* One table row. */
typedef struct {
    double z;                /* slab centre */
    double y, x;             /* the curve: fused estimate */
    double y1, x1;           /* primary (normal-line) centre */
    double y2, x2;           /* secondary (polar harmonic) centre */
    int    status;           /* AxisRowStatus */
    int    primary_ok, secondary_ok;
    size_t n_in;             /* primary samples inside the local window */
    double cond;             /* lambda_min / lambda_max of the in-window normal scatter */
    int    sectors;          /* 30-degree sectors around the centre holding kept samples (0..12) */
    double void_ratio;       /* wrap-area density inside r < AXIS_TRACE_VOID_R over the
                              * density in the annulus to 3x that radius: ~0 at the
                              * umbilicus void, ~1 at an evolute cusp inside the stack */
    double resid_p50;        /* normal-line residual median, vox */
    double d1;               /* secondary eccentricity at the final centre, vox */
    size_t ridges;           /* ridges spanning >= the minimum arc */
    double pitch;            /* signed mean ridge slope, vox per turn (sign = sense) */
    double r_wall;           /* innermost ridge radius, vox (0 = none) */
    double agree;            /* |primary - secondary|, vox */
    double jump;             /* |primary - prior| the driver passed, vox */
} AxisRow;

int  AxisSlab_init(AxisSlab *s, double z0, double z1,
                   size_t face_cap, size_t vert_cap, uint32_t seed);
void AxisSlab_dispose(AxisSlab *s);
void AxisSlab_add_face(AxisSlab *s, const float a_zyx[3],
                       const float b_zyx[3], const float c_zyx[3]);
void AxisSlab_add_vertex(AxisSlab *s, const float p_zyx[3]);

/* Primary estimator about a prior centre.  Fills row->y1/x1, n_in, cond,
 * resid_p50, primary_ok.  Returns 0 when the lock rule holds, -1 to abstain
 * (the row's diagnostics are still filled). */
int AxisTrack_primary(const AxisSlab *s, double prior_y, double prior_x,
                      AxisRow *row);

/* Secondary estimator from a start centre.  Fills row->y2/x2, d1, ridges,
 * pitch, r_wall, secondary_ok.  Returns 0 on lock, -1 to abstain. */
int AxisTrack_secondary(const AxisSlab *s, double start_y, double start_x,
                        AxisRow *row);

/* Derive a seed on one slab from a grid of starts (spacing
 * AXIS_TRACE_SEED_GRID over the slab's vertex bounding box): every start is
 * iterated through the primary estimator (the local window follows the
 * estimate) until it stops moving, converged centres are clustered (within
 * 30 vox) and the cluster with the most in-window support and lowest
 * residual wins.  Reports up to max_clusters clusters (y, x, members,
 * resid_p50 in cl_d1, support in cl_ridges) into the caller's arrays.
 * MEASURED 2026-09-02 on the full 21x grid: the flattened wraps' evolute
 * cusps are normal-line attractors as isotropic as the umbilicus (cond
 * 0.42 vs 0.21) with more support; what separates the umbilicus is the
 * VOID it sits in (wrap-area density inside 50 vox ~0 vs ~1 in the stack).
 * Candidates need void_ratio < AXIS_TRACE_VOID_MAX; the lowest residual
 * among them wins.
 * Returns 0 with the winner in *seed_y/*seed_x, -1 when no start locked. */
int AxisTrack_seed(const AxisSlab *s, double *seed_y, double *seed_x,
                   double *cl_y, double *cl_x, size_t *cl_members,
                   double *cl_d1, size_t *cl_ridges, size_t max_clusters,
                   size_t *n_clusters);

/* Same, with the start grid restricted to the square of half-width win_r
 * about (win_y, win_x) (a prior curve's centre); win_r <= 0 = the whole
 * slab bounding box. */
int AxisTrack_seed_in(const AxisSlab *s, double win_y, double win_x, double win_r,
                      double *seed_y, double *seed_x,
                      double *cl_y, double *cl_x, size_t *cl_members,
                      double *cl_d1, size_t *cl_ridges, size_t max_clusters,
                      size_t *n_clusters);

/* Post-process a finished row array (in z order): retain primary positions
 * (secondary disagreement remains diagnostic), slope limit (abstains, never clamps),
 * interpolate gaps of <= AXIS_TABLE_MAX_GAP_ROWS, mark rows beyond the last
 * lock EXTRAP.  Returns the number of locked rows. */
size_t AxisTrack_finish(AxisRow *rows, size_t n);

/* Write the provenance CSV (header lines start with '#'; data rows are
 * z,y,x,status,n,cond,resid_p50,d1,y1,x1,y2,x2,pitch,r_wall so the first
 * three columns load through AxisWarp_load_csv). */
int AxisTrack_write_csv(const char *path, const AxisRow *rows, size_t n,
                        const char *header_lines);

/* Winding sense from the locked rows' signed pitch: +1 / -1 / 0 (no
 * evidence).  *agree = fraction of locked rows agreeing with the majority. */
int AxisTrack_sense(const AxisRow *rows, size_t n, double *agree);

/* 0 quiet, 1 per-call summaries, 2 per-iteration traces (stderr). */
extern int AxisTrack_verbose;

/* Debug view: the area-weighted (theta, r) raster of one slab about
 * (cy, cx), theta across (2 deg per column from -pi), r down (1 vox per
 * row, up to r_max), detected peaks in red, peaks that joined a ridge of
 * the minimum arc in green.  Wraps must be horizontal when the centre is
 * right.  Returns 0 on success. */
int AxisTrack_polar_debug(const AxisSlab *s, double cy, double cx,
                          int r_max, const char *png_path);

int AxisTrack_selftest(void);

#endif
