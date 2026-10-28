#ifndef FIBER_MAP_INCLUDED
#define FIBER_MAP_INCLUDED

#include <stddef.h>
#include <stdint.h>
#include "../common/arena.h"

/* Measured unit PHYSICAL fibre arms, pulled back to canvas coordinates.
 * Applying a canvas->physical-UV Jacobian to h/v gives their output vectors.
 * Weight is physical source-area quadrature times relative evidence. */
typedef struct FiberMapDirection {
    double column, row, h[2], v[2], weight;
} FiberMapDirection;

typedef struct FiberMapDirections {
    size_t count;
    const FiberMapDirection *sample;
} FiberMapDirections;

/* Native fibre-directions-v1 CSV; positive finite mass, independent positive
 * arm basis, exact original canvas coordinates. Nonzero on missing/bad data. */
int FiberMap_read_directions(Arena_T arena, const char *path, FiberMapDirections *out);

/* A continuous forward map for applying a fibre frame to an existing mesh.
 * Missing grid nodes are harmonic extensions of the measured displacement.
 * A failed cell or source triangle releases nearby map constraints, then the
 * same field is solved again. If that exhausts the samples, a constrained
 * least-squares fit projects every original face's physical Jacobian and every
 * bilinear cell's orientation path while solving the same shared field.
 * Original material vertices/faces are untouched.
 * Success certifies every bilinear cell and every source triangle throughout
 * the straight path from the baseline frame, plus each original triangle's
 * physical metric at the endpoint. It is not a material-contact or fibre-quality
 * certificate. A baseline that worsens source metric is replaced by a rigid
 * physical frame; the measured constraints still define the requested map.
 * Candidate basis refinement preserves the padded domain and measured nodes;
 * additional unknowns have no invented fibre targets. */
typedef struct FiberMapGrid {
    int rows, cols;
    double row0, col0, cell;
    const float *u, *v;
    double affine[6];
    double source_scale[2];      /* canvas pixels per original physical UV unit */
    const FiberMapDirections *directions; /* optional direct alignment objective */
    /* Optional caller work budget/checkpoint. Zero preserves existing policy.
     * A checkpoint sees only a full original-face/cell certified incumbent;
     * return 1 when the caller's independently measured target is met, 0 to
     * continue, or -1 on failure. A work limit never means convergence. */
    int max_iterations;
    double max_seconds;
    void *checkpoint_context;
    int (*checkpoint)(void *,const double *,size_t,int);
} FiberMapGrid;

typedef struct FiberMapResult {
    int rows, cols;
    double row0, col0, cell;
    double *uv;                  /* [rows*cols*2], all finite on success */
} FiberMapResult;

typedef struct FiberMapStats {
    size_t source_nodes, extended_nodes, relaxed_nodes;
    size_t baseline_pinned_nodes;
    size_t unsafe_cells_before, unsafe_faces_before;
    size_t metric_faces_before, metric_faces_after, source_metric_bad, affine_metric_bad;
    int rigid_baseline;
    int passes;
    double node_correction_rms, node_correction_max;
    double seconds;
    int metric_iterations;
    double metric_fit_before, metric_fit_after;
    size_t direction_samples;
    int direction_fit;
    double direction_angle_before, direction_angle_after; /* weighted radians^2, both arms */
} FiberMapStats;

/* xyz is the unchanged [nv*3] physical source geometry. canvas_uv and output_uv
 * are [nv*2], in the input/output map's pixel frames (output is physical units).
 * The arena owns result.uv and scratch. Nonzero leaves no certified result. */
int FiberMap_fit(Arena_T arena, const FiberMapGrid *grid,
                  const float *xyz, const double *canvas_uv, size_t nv,
                 const int32_t *faces, size_t nf, double *output_uv,
                 FiberMapResult *result, FiberMapStats *stats);

int FiberMap_selftest(void);

#endif
