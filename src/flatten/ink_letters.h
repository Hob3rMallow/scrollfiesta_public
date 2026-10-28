/* ============================================================================
 * ink_letters.h -- letter boxes, text rows and the row rectification warp of
 * an ink map (ink_closeup's letters stage).
 *
 * Letters: the ink probability map (uint8, p = v / 255), Gaussian-smoothed,
 * hysteresis-thresholded (seeds p >= t_hi grown through p >= t_lo, inside the
 * mask), 8-connected; components closer than close_px merge into one letter
 * (a letter's strokes); letters outside [min_area, max_area] px^2 are dropped.
 * Each letter carries a p-weighted PCA oriented box (OBB: 2-98% extents on its
 * axes) and its axis-aligned box.
 *
 * Rows: Greek book hands write letters of one height along straight lines, so
 * letters whose centres lie within row_tol_h letter heights of each other
 * vertically and row_link_h heights horizontally chain into a row; a row of at
 * least row_min letters takes a robust quadratic y = f(x) of its centres.
 *
 * Rectification: a displacement field D on a coarse grid (warp_px) maps the
 * map's (x, y) to (x, y) + D so every row's letter centres land on one height
 * (a free level per row), while D stays as conformal as possible (the
 * Cauchy-Riemann residual of grad D -- letters rotate with their row instead
 * of shearing), smooth (second differences) and small (a weak anchor).  A
 * linear least-squares problem, solved on TAUCS.  Applied to the UV of the
 * close-up surface, it re-renders the text in straight rows from 3-D.
 * ==========================================================================*/
#ifndef INK_LETTERS_H
#define INK_LETTERS_H

#include <stddef.h>
#include <stdint.h>

typedef struct InkLettersOpts {
    double smooth_px;        /* Gaussian sigma before thresholding */
    double t_hi, t_lo;       /* hysteresis on p */
    int    close_px;         /* components within this distance merge */
    double min_area, max_area;
    double row_link_h, row_tol_h;   /* in median letter heights */
    int    row_min;          /* letters a row needs */
    double warp_px;          /* rectification grid step */
    double w_row, w_conf, w_smooth, w_anchor;   /* rectification weights */
} InkLettersOpts;

void InkLettersOpts_default(InkLettersOpts *o);

typedef struct InkLetter {
    double cx, cy;           /* p-weighted centre (px, pixel centres at index + 0.5) */
    double angle;            /* OBB major axis angle (rad, x toward y) */
    double ha, hb;           /* OBB half extents along the major / minor axis */
    double bx0, by0, bx1, by1;   /* axis-aligned box [bx0, bx1) x [by0, by1) */
    double area, mean_p, max_p;
    int    row;              /* row index, -1 none */
} InkLetter;

typedef struct InkRow {
    double a, b, c;          /* y = a + b x + c x^2 */
    double x0, x1;           /* extent of its letter centres */
    double height;           /* median letter box height */
    double rms;              /* residual of its centres, px */
    size_t n;                /* letters */
} InkRow;

/* Detect letters.  mask (nullable) limits the search.  *letters is malloc'd.  0 on success. */
int InkLetters_detect(const uint8_t *ink, const uint8_t *mask, size_t W, size_t H, const InkLettersOpts *o,
                      InkLetter **letters, size_t *n);

/* Chain letters into rows (sets letters[i].row); *rows is malloc'd; *h_med the median letter height. */
int InkLetters_rows(InkLetter *letters, size_t n, const InkLettersOpts *o, InkRow **rows, size_t *nrows,
                    double *h_med);

typedef struct InkWarp {
    size_t nx, ny;           /* nodes */
    double step;             /* px */
    double *dx, *dy;         /* [nx*ny] displacement at node (i, j) = (i step, j step) */
    double *level;           /* [nrows] the straightened height of each row */
    size_t nrows;
    double rms_before, rms_after;   /* letter centres about their row's level/curve, px */
    double scale_min, scale_max;    /* extreme local scale of the warp over the letters' cells */
    double rot_max_deg;             /* largest local rotation over the letters' cells */
} InkWarp;

/* Solve the rectification warp over a W x H map (see the file comment).  0 on success. */
int  InkLetters_warp(const InkLetter *letters, size_t n, const InkRow *rows, size_t nrows, size_t W, size_t H,
                     const InkLettersOpts *o, InkWarp *w);
void InkWarp_free(InkWarp *w);
/* Bilinear displacement at (x, y) (clamped to the grid). */
void InkWarp_eval(const InkWarp *w, double x, double y, double *dx, double *dy);

/* ---- drawing on an RGB raster (overlays, galleries) ------------------------------------------------- */
void InkDraw_line(uint8_t *rgb, size_t W, size_t H, double x0, double y0, double x1, double y1, int thick,
                  uint8_t r, uint8_t g, uint8_t b);
void InkDraw_obb(uint8_t *rgb, size_t W, size_t H, const InkLetter *L, int thick, uint8_t r, uint8_t g, uint8_t b);
/* A 5x7 bitmap text (digits, lower-case letters, space, '.', '-', '=', ':', '/', '%') scaled by `scale`,
 * top-left at (x, y); characters outside the set draw as blanks. */
void InkDraw_text(uint8_t *rgb, size_t W, size_t H, int x, int y, int scale, const char *s,
                  uint8_t r, uint8_t g, uint8_t b);

int InkLetters_selftest(void);

#endif /* INK_LETTERS_H */
