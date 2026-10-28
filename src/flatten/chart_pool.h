/* ============================================================================
 * chart_pool.h -- an assembler chart store as an ink_closeup adoption pool.
 *
 * The assembler's stage 1 (sheet_assemble --stop-after clean, region_job's
 * `cleaner`) leaves one chart store per cube, <store>/<cube id>.asc
 * (src/assemble/asm_store.c, format ASMCHT02): every connected component of
 * the cube's mesh with its own intrinsic flattening.  An ink_closeup pool is a
 * layout directory: stage6_repaired.vmesh (vertices, faces, per-vertex UV)
 * with stage6_repaired_chart.i32 (chart id per vertex), each chart's vertices
 * and faces contiguous, and its UV any rigid frame of an isometric map
 * (closeup_adopt places every chart by a 2-D rigid fit of its contacts).
 *
 * ChartPool_from_store writes such a pool from a chart store: the charts the
 * assembler keeps for layout (no BLOB, EXTRAS, NONMANIFOLD, FLAT_FAILED or
 * TINY flag, with a flattening), their own UV widened to float64, pool chart
 * ids 0..n-1 in cube order, and the <vmesh>.chartidx sidecar, so a close-up
 * pages in only the charts near its surface.  The meshes of new papyrus
 * (region jobs outside an assembled box) become adoptable without an
 * assembly: no relation, pose or placement is computed or needed.
 *
 * Coordinates as stored: vertices (z, y, x) world voxels (float32).
 * ==========================================================================*/
#ifndef CHART_POOL_H
#define CHART_POOL_H

#include <stddef.h>

typedef struct ChartPoolStats {
    size_t cubes;          /* <cube id>.asc files in the box */
    size_t cubes_failed;   /* unreadable / malformed files (skipped whole) */
    size_t charts;         /* charts written */
    size_t charts_skipped; /* charts the assembler keeps out of its layout */
    size_t verts, faces;   /* written */
    double seconds;
} ChartPoolStats;

/* The pool of the chart store store_dir (every <z#####_y#####_x#####>.asc whose cube origin lies in box =
 * {z0, z1, y0, y1, x0, x1}, half-open; NULL = every cube) into out_dir (created): stage6_repaired.vmesh,
 * stage6_repaired_chart.i32, the .chartidx sidecar and chart_pool_report.json.  0 on success, -1 on error
 * (nothing published), 1 when no chart qualifies (nothing written). */
int ChartPool_from_store(const char *store_dir, const char *out_dir, const long *box, ChartPoolStats *st);

int ChartPool_selftest(void);

#endif /* CHART_POOL_H */
