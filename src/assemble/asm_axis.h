#ifndef ASM_AXIS_INCLUDED
#define ASM_AXIS_INCLUDED

#include <stddef.h>

#include "../common/arena.h"

/* ============================================================================
 * asm_axis.h -- the scroll axis as a polyline, evaluated PER POINT.
 *
 * A derived axis table (configs/axis.csv: z,y,x per z-slab) is a
 * polyline; a single straight line (the configured umbilicus along +z) is
 * the two-row degenerate case.  Everything the assembler asks of the axis
 * comes from projecting a point onto that polyline:
 *
 *   s     arc length along the axis to the foot of the point (the AXIAL
 *         coordinate: the sheet is perpendicular to the local axis, so its v
 *         is s up to a constant);
 *   r     distance from the axis line (the radius: one turn is 2 pi r);
 *   dir   the unit tangent at the foot (the local axis direction);
 *   frame (s, e1 . (p - foot), e2 . (p - foot)) with (dir, e1, e2) a
 *         right-handed orthonormal frame, e1 the world +y direction made
 *         perpendicular to the tangent: the measurement frame in which the
 *         axis is straight along +z through the origin.
 *
 * A single point per box (the umbilicus, or a line fitted over the box) was
 * measured wrong by up to 190 vox at the ends of the 21x5x5 and drifts 236
 * vox in y across the 10x10x10 (changelog 2026-09-08/09); every consumer --
 * the pose prior, placement, the verdict frame -- evaluates the table at its
 * own point instead.
 * ==========================================================================*/

typedef struct AsmAxis {
    size_t  n;         /* polyline vertices (the table densified to 8-vox z steps), at least 2 */
    size_t  n_rows;    /* rows of the table it came from */
    double  fit_rms;   /* the cubic's residual against the rows, vox */
    double  r_curv_min;/* the polyline's smallest curvature radius, vox: the frame is rigid only for r << this */
    double *p;         /* [n*3] (z,y,x) */
    double *s;         /* [n] cumulative arc length from row 0 */
    double *t;         /* [n*3] continuous unit tangent per row (mean of the adjacent segments) */
    int     is_line;   /* the straight fallback (point + direction) */
} AsmAxis;

/* Reads a z,y,x table ('#' lines skipped, rows sorted by z, duplicates
 * dropped).  NULL when the file is absent or holds fewer than two rows. */
AsmAxis *AsmAxis_load(Arena_T arena, const char *path);

/* The polyline through z,y,x rows ([n*3], any order; duplicate z dropped):
 * densified at 8-vox z steps by the Gaussian-weighted local quadratic
 * (LOESS) of AsmAxis_load, with optional per-row weights ([n], NULL = 1).
 * NULL with fewer than two distinct rows.  AsmAxis_load is this on the
 * file's rows, so a derived axis and a loaded table go through one path. */
AsmAxis *AsmAxis_from_rows(Arena_T arena, const double *rows, const double *row_w, size_t n);

/* The straight axis through `point` along the unit `dir` (z,y,x), spanning
 * z0..z1 (the extrapolation is linear beyond either end anyway). */
AsmAxis *AsmAxis_line(Arena_T arena, const double point[3], const double dir[3], double z0, double z1);

/* Nearest point of the polyline to p (linear extrapolation beyond the ends).
 * Any output may be NULL. */
void AsmAxis_project(const AsmAxis *a, const double p[3], double *s, double *r, double dir[3], double foot[3]);

/* The measurement frame of p: out = (s, e1 . (p - foot), e2 . (p - foot)). */
void AsmAxis_frame(const AsmAxis *a, const double p[3], double out[3]);

/* Mean unit tangent over the rows within z0..z1 (all rows when none fall in). */
void AsmAxis_mean_dir(const AsmAxis *a, double z0, double z1, double dir[3]);

int AsmAxis_selftest(void);

#endif
