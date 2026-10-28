#ifndef FIBER_FRAME_INCLUDED
#define FIBER_FRAME_INCLUDED

#include <stddef.h>
#include <stdint.h>
#include "../common/arena.h"
#include "fiber_map.h"

/* ============================================================================
 * fiber_frame.h -- fibre-aligned re-flattening of a tifxyz canvas from a 3-D fibre field (2026-10-15; the user, on
 * PHerc0211 line 11's 19,959 px close-up page, a 2,089 px arc by the fibres: "Do all the things").
 *
 * The fibre field is a model's per-voxel fibre presence p and sign-free unit direction d (villa fiber3d SD2 as an
 * uncompressed uint8 zarr: <root>/presence, <root>/nx, <root>/ny; nx = round(127 d_x + 128), d folded to d_z >= 0,
 * the official field's encoding).  On the canvas:
 *
 *   field    every FIBER_FLATTEN_SAMPLE_PX px of the window, at -HALF_DEPTH..+HALF_DEPTH vox along the render
 *            normal -(dP/du x dP/dv) (vc_render --flip-normals), each fibre direction is expressed in an orthonormal
 *            physical tangent frame as an angle psi; z4 = sum p exp(4 i psi) on a grid of
 *            FIBER_FLATTEN_CELL_PX px, Gaussian FIBER_FLATTEN_SIGMA_PX -> the cross field.  arg(z4) / 4 selects
 *            the arm nearest the physical rows' axis (H, round the scroll); it is then pulled back to canvas
 *            coordinates for theta. The other physical arm V is perpendicular to H.
 *   solve    the new coordinates (u', v') at the grid's nodes: du'=H.dP, dv'=V.dP, using the unit physical arms
 *            and the full source metric, including shear, by finite differences on node edges, weighted least
 *            squares (TAUCS Cholesky).  The gauge fixes only the translation (FIBER_FLATTEN_EPS): the page may turn
 *            and bend as far as the fibres say.  Narrow cracks and enclosed holes (closing FIBER_FLATTEN_GAP_CELLS +
 *            hole fill) are carried as nodes tied to the OLD layout with FIBER_FLATTEN_GAP_WEIGHT, so pieces
 *            joined by thin bridges cannot float apart.
 *   audit    H from the rows (old: the field; new: its image under the map), the fibre climb (a fibre's height as
 *            it crosses the page: the column-median H angle integrated, old by old column, new by u'), the map's
 *            rigid part and what is left (page, FIBER_FLATTEN_TILE_PX patches), area ratio, folded nodes.
 *   inverse  each new pixel's old canvas position (quasi-Newton on the bilinear forward map, from the best affine;
 *            past the solved nodes the map carries the nearest node's displacement).
 * ==========================================================================*/

typedef struct FiberVolume *FiberVolume_T;

/* Open <root>/presence, <root>/nx, <root>/ny (each an uncompressed uint8 zarr v2 array at <name>/0).  0 on success. */
int  FiberVolume_open(const char *root, FiberVolume_T *out);
void FiberVolume_close(FiberVolume_T *fv);
/* Presence in [0, 1] at the voxel nearest (x, y, z) (scan frame, voxels) and its unit direction d (x, y, z, d_z >= 0);
 * 0 (and d = 0) where the field is absent.  Not thread-safe (a chunk cache per array). */
double FiberVolume_sample(FiberVolume_T fv, double x, double y, double z, double d[3]);

typedef struct {
    /* inputs: the canvas (H x W nodes, x/y/z per node in the scan frame) and the window [x0, x1) x [y0, y1) */
    int W, H;
    int x0, y0, x1, y1;
    /* the grid: gh x gw nodes, node (i, j) at canvas row gy0 + i * cell, column gx0 + j * cell */
    int gh, gw, cell, gy0, gx0;
    float   *theta;    /* [gh*gw] H direction in canvas coordinates; arm selected in the physical frame */
    float   *weight;   /* [gh*gw] smoothed evidence (sum of presence) */
    float   *coherence; /* [gh*gw] magnitude of the mean fourfold phase, [0,1]; diagnostic only */
    uint8_t *real;     /* [gh*gw] 1 = the canvas is valid at the node and its tangents are finite */
    uint8_t *closed;   /* [gh*gw] real + gap nodes (the solve's unknowns) */
    double  *pu, *pv;  /* [gh*gw*3] the canvas tangents at the nodes (x, y, z per canvas px) */
    double  *U, *V;    /* [gh*gw] the solved new coordinates (NaN outside closed) */
    size_t samples, samples_with_fibres;
} FiberFrame;

typedef struct {
    double p50, p90, p99;
} FiberQuant;

typedef struct {
    /* the solve */
    size_t nodes_real, nodes_gap, edges;
    double rms_residual_u, rms_residual_v;   /* px per node edge, real edges */
    double L_h_p50, L_v_p50;                 /* vox per canvas px along the arms */
    double evidence_share;                   /* real nodes with evidence */
    /* H from the rows */
    double H_old_signed_mean, H_new_signed_mean;
    FiberQuant H_old_abs, H_new_abs;         /* deg */
    /* the fibre climb */
    double climb_old_range, climb_new_range;         /* px */
    double climb_old_bow, climb_new_bow;             /* px, after the best straight line */
    double climb_old_end, climb_new_end;             /* px, end to end */
    double climb_old_width, climb_new_width;         /* px */
    double col_angle_old_p10, col_angle_old_p90, col_angle_new_p10, col_angle_new_p90;
    /* the map */
    double rotation_deg;
    FiberQuant nonrigid_px;  double nonrigid_max;
    FiberQuant tile_nonrigid_px, tile_rotation_deg;
    size_t tiles, flipped;
    double area_p1, area_p10, area_p50, area_p90, area_p99;
    /* the affine fit used by the inverse: (u', v') = A (col, row) + b */
    double A[4], b[2];
} FiberFrameStats;

/* Alignment of measured physical fibre arms under an externally retained
 * bilinear map. Quantiles pair before/after at the same supported evidence
 * nodes. Each node takes the worst angle over all incident map cells, so a
 * refined map cannot hide opposing one-sided slopes in a coarse difference.
 * This samples the field grid; it is not a whole-mesh geometric certificate. */
typedef struct FiberFrameMapStats {
    size_t real_nodes, evidence_nodes, supported_nodes;
    size_t missing_nodes, invalid_nodes, folded_nodes, cell_sides;
    FiberQuant H_old_abs, H_new_abs, V_old_abs, V_new_abs;
    double determinant_min, determinant_max;
} FiberFrameMapStats;

typedef struct FiberFrameBranchExample {
    int col_a, row_a, col_b, row_b, quarter;
    double selected_deg, residual_deg, normal_deg, coherence, maximum_world_step;
} FiberFrameBranchExample;

/* Read-only comparison across supported canvas-grid paths. Tangent directions
 * are transported by the shortest normal rotation. This is NOT a certificate
 * of original-mesh adjacency or same-leaf identity. Antipodal normals are
 * explicitly unmeasurable. No branch is changed and no acceptance gate uses
 * these counts. Examples retain at most eight switches with greatest coherence. */
typedef struct FiberFrameBranchStats {
    size_t evidence_nodes, neighbor_pairs, unsupported_paths, invalid_frames;
    size_t antipodal_normals, compared_pairs, nearest_quarter[4], switches;
    FiberQuant selected_deg, residual_deg, switch_normal_deg, switch_coherence;
    FiberQuant switch_maximum_world_step;
    size_t examples;
    FiberFrameBranchExample example[8];
} FiberFrameBranchStats;

/* The cross field on the window.  X: [H*W*3] float (x, y, z per node; invalid nodes have valid[i] == 0).
 * Arena-allocates the grid arrays.  0 on success, -1 on bad input. */
int FiberFrame_field(Arena_T arena, const float *X, const uint8_t *valid, int H, int W,
                     int x0, int y0, int x1, int y1, FiberVolume_T fv, FiberFrame *ff);

/* Solve the new coordinates (TAUCS).  0 on success. */
int FiberFrame_solve(Arena_T arena, FiberFrame *ff, FiberFrameStats *st);

/* The audit (fills the rest of st).  Writes the two climb profiles to profile_csv when it is not NULL. */
int FiberFrame_audit(Arena_T arena, const FiberFrame *ff, FiberFrameStats *st, const char *profile_csv);

/* Read-only audit of the actual forward map, with exact bilinear derivatives.
 * map coordinates use the original canvas frame; output translations do not
 * affect the result. A zero return means the audit ran, not that the map passed.
 * Missing/invalid/folded counts must be inspected, including when no node can
 * be measured (then the angle quantiles are NaN). ff need not have been solved. */
int FiberFrame_audit_map(Arena_T arena, const FiberFrame *ff,
                         const FiberMapResult *map, FiberFrameMapStats *stats);

int FiberFrame_audit_branches(Arena_T arena, const FiberFrame *ff,
                              const float *X, const uint8_t *valid, int H, int W,
                              FiberFrameBranchStats *stats);

/* Save measured physical fibre arms for the constrained alignment solve.
 * Only real nodes meeting the same evidence rule as the audit are written. */
int FiberFrame_write_directions(Arena_T arena, const FiberFrame *ff, const char *path);

/* The new page: Hn x Wn pixels, pixel (r, c) at (u', v') = (u_min + c, v_min + r).  old_rc[(r*Wn + c)*2 + {0,1}] =
 * the old canvas (row, col), NaN where the pixel has no valid preimage.  threads <= 0 -> 1. */
int FiberFrame_inverse(Arena_T arena, const FiberFrame *ff, const FiberFrameStats *st, int threads,
                       int *Hn, int *Wn, double *u_min, double *v_min, float **old_rc);

/* Bilinear resample of the canvas X at the old positions: Xn [Hn*Wn*3] (-1 where invalid: a corner invalid, the
 * position outside the window, or no preimage), vn [Hn*Wn] (1 valid), then vn eroded by one pixel (vc_render blanks
 * pixels next to invalid nodes).  Returns the valid count in *n_valid. */
int FiberFrame_resample(Arena_T arena, const float *X, const uint8_t *valid, int H, int W, int x0, int y0, int x1, int y1,
                        const float *old_rc, int Hn, int Wn, float **Xn, uint8_t **vn, size_t *n_valid);

/* Unit tests: helpers (closing, hole fill, smoothing, circular median) and a planted bowed field on a tilted plane
 * that must come out straight.  0 = pass. */
int FiberFrame_selftest(void);

#endif
