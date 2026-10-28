/* ============================================================================
 * closeup_f0bridge.h -- a close-up surface carried across a mesh gap along the
 * F0 prediction.
 *
 * Where the sheets of a crushed scroll converge (the corners of the lens around
 * an eye) the per-cube mesher cuts the touching sheets apart with a gap
 * exclusion zone, so the charts of one sheet end 10-30 vox before the charts
 * that continue it.  No chart touches across such a gap, and the fitted fill
 * extrapolates straight while the sheet bends: closeup_adopt (3-D contact)
 * cannot cross it.  The F0 prediction does trace the sheet through the gap.
 *
 * One RING of the bridge: every FILL node 4-adjacent to the data (or to a node
 * snapped earlier in the same pass), inside the target box and no deeper than
 * depth_px beyond the CHART data (data observed by a mesh chart, not by an
 * earlier bridge), reads the F0 profile along its normal through its predicted
 * position (the fill, moved by the mean correction of its snapped neighbours);
 * the foreground run whose centre is nearest becomes the snapped point when
 *   - its centre lies within shift_vox of the prediction,
 *   - it is no thicker than max_run_vox (thicker = touching sheets, ambiguous),
 *   - the segment from the neighbour it continues to the point stays in the
 *     foreground (the same sheet, not a jump across a gap to the next one).
 *   - it does not lie on the surface already there (data, or a node snapped
 *     earlier in the pass) within dup_vox while far from that node along the
 *     surface (more than 3 x their distance + dup_margin_px apart in the frame):
 *     two fronts of the bridge meeting round an obstacle, or the bridge coming
 *     back over the seed, would otherwise lay the same papyrus twice.  Layers
 *     of a real fold lie at least a papyrus thickness apart and pass.
 *   - it does not stack onto the surface already there: no data (or earlier
 *     snap) within layer_vox of it along its normal, more across than beside
 *     it, far from that node along the surface (the same frame test) -- a fold
 *     back over the surface or a step onto the next wrap.  Round a bend tighter
 *     than the frame margin the two arms are near along the surface and pass.
 * A pass runs `rings` rings, then the snapped nodes with the data around them
 * become observation faces (grid quads) with chart id f0_chart, for the
 * caller's refit.  Nothing is snapped where the prediction is missing, merged
 * or broken: the bridge stops there, as the sheet's evidence does.
 *
 * F0 volume: 128^3 uint8 cubes, <dir>/z#####_y#####_x#####.tif (the region
 * jobs' cubes_PRED, cube id = voxel origin), nonzero = papyrus; several
 * directories are searched in order.  Coordinates (z, y, x) voxels.
 * ==========================================================================*/
#ifndef CLOSEUP_F0BRIDGE_H
#define CLOSEUP_F0BRIDGE_H

#include <stddef.h>
#include <stdint.h>

#include "closeup_surface.h"

typedef struct CloseupF0 *CloseupF0;

/* A cube cache over dirs (copied), at most max_cubes decoded cubes (2 MiB each) kept.  NULL on no memory. */
CloseupF0 CloseupF0_open(const char *const *dirs, size_t ndirs, size_t max_cubes);
void      CloseupF0_close(CloseupF0 *f);
/* Load (serially) every cube meeting the box lo..hi (voxels, inclusive) not yet cached; returns how many cubes of
 * the box exist in some directory, -1 on error. */
long      CloseupF0_preload(CloseupF0 f, const double lo[3], const double hi[3]);
/* The voxel nearest (z, y, x): 1 papyrus, 0 not, -1 no cube there (or not preloaded).  Read-only: safe from many
 * threads once the cubes are preloaded. */
int       CloseupF0_at(const CloseupF0 f, double z, double y, double x);

typedef struct CloseupF0BridgeOpts {
    double window_vox;    /* the profile reaches this far either side of the prediction along the normal */
    double step_vox;      /* profile sample step */
    double shift_vox;     /* the chosen run's centre at most this far from the prediction */
    double max_run_vox;   /* a run thicker than this (touching sheets) is never snapped */
    double depth_px;      /* no snap deeper than this beyond the chart data */
    double dup_vox;       /* the duplicate-layer guard's reach (<= 0: off) */
    double dup_margin_px; /* ... and its frame margin */
    double ambig_vox;     /* a second run less than this farther than the nearest: ambiguous, no snap (<= 0: off) */
    double layer_vox;     /* the layer guard's reach along the normal (<= 0: off) */
    int    rings;         /* rings per pass (one refit per pass) */
    double box[4];        /* target box in frame UV: u0, v0, u1, v1 */
} CloseupF0BridgeOpts;

void CloseupF0BridgeOpts_default(CloseupF0BridgeOpts *o);

typedef struct CloseupF0BridgeStats {
    size_t candidates;    /* ring nodes tried */
    size_t snapped;       /* snapped nodes */
    size_t no_cube;       /* no F0 cube under the profile */
    size_t no_run;        /* no foreground in the profile */
    size_t shifted;       /* nearest run beyond shift_vox */
    size_t thick;         /* nearest run thicker than max_run_vox */
    size_t broken;        /* the segment to the neighbour leaves the foreground */
    size_t deep;          /* ring nodes beyond depth_px (not tried) */
    size_t dup;           /* the snapped point would lie on the surface already there, far away along it */
    size_t ambiguous;     /* a second run nearly as near as the nearest (converging sheets) */
    size_t layered;       /* the snapped point would stack onto the surface along its normal, far away along it */
    size_t faces;         /* observation faces appended */
    size_t observations;  /* observations appended */
} CloseupF0BridgeStats;

/* One pass (o->rings rings) on the fitted grid g: snapped nodes -> observation faces appended to obs with chart id
 * f0_chart.  Data nodes whose g->chart is f0_chart count as bridge data (depth grows through them), every other
 * DATA node as chart data (depth 0).  The F0 cubes near the boundary are preloaded here.  0 on success (st filled;
 * st->snapped may be 0), -1 on error. */
int CloseupF0Bridge_pass(const CloseupGrid *g, CloseupF0 f0, const CloseupF0BridgeOpts *o, int32_t f0_chart,
                         CloseupObs *obs, CloseupF0BridgeStats *st);

int CloseupF0Bridge_selftest(void);

#endif /* CLOSEUP_F0BRIDGE_H */
