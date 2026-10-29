#ifndef ASM_REPORT_INCLUDED
#define ASM_REPORT_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "asm_types.h"
#include "asm_conflict.h"

/* ============================================================================
 * asm_report.h -- PNG previews and ledgers for every assembly stage.
 *
 * Every stage drops a viewable bitmap (house rule).  The rasterizer is a
 * tiny software triangle filler on an RGB canvas; colours are per cube (so
 * seams read as colour changes) with a red tint proportional to a chart's
 * flattening stress and magenta for flipped faces.
 * ==========================================================================*/

typedef struct AsmCanvas {
    uint8_t *rgb;    /* [w*h*3] */
    int      w, h;
} AsmCanvas;

AsmCanvas AsmCanvas_new(Arena_T arena, int w, int h, uint8_t r, uint8_t g, uint8_t b);
void      AsmCanvas_tri(AsmCanvas *c, double x0, double y0, double x1, double y1,
                        double x2, double y2, const uint8_t rgb[3]);
void      AsmCanvas_line(AsmCanvas *c, double x0, double y0, double x1, double y1,
                         const uint8_t rgb[3]);
void      AsmCanvas_rect(AsmCanvas *c, double x0, double y0, double x1, double y1,
                         const uint8_t rgb[3]);
int       AsmCanvas_write(const AsmCanvas *c, const char *path);

/* Deterministic colour for a cube / chart index. */
void AsmReport_colour(uint32_t key, uint8_t rgb[3]);

/* Flattened charts packed into a grid of cells (largest first, at most
 * max_cells charts), each scaled uniformly to its cell.  Excluded charts are
 * not drawn.  Returns 0 on success. */
int AsmReport_charts_atlas_png(Arena_T arena, const char *path,
                               const AsmChart *charts, size_t n,
                               int cell_px, size_t max_cells);

/* Global layout: components stacked top to bottom (largest first, at most
 * max_components of them), every placed chart drawn through its pose in its
 * component's frame.  Contradiction cells (nullable) are painted red on top.
 * The canvas is scaled so the longer side fits max_px (never upsampled).
 * Returns 0 on success. */
int AsmReport_layout_png(Arena_T arena, const char *path,
                         const AsmChart *charts, size_t n,
                         const AsmConflictCell *cells, size_t n_cells, double cell_size,
                         int max_px, size_t max_components);

/* Pose ledger of every chart in the layout (stage,chart,component,u,v,theta,
 * mirror,cu,cv,area,cube,cz,cy,cx): the frame origin, the chart's layout
 * centroid, its cube and 3-D centroid.  Returns 0 on success. */
int AsmReport_poses_csv(const AsmRun *run, const char *path, const char *stage);

#endif
