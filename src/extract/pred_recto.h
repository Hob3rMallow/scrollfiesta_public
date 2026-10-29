#ifndef PRED_RECTO_INCLUDED
#define PRED_RECTO_INCLUDED

#include "../common/arena.h"
#include "../common/mesh_types.h"
#include "../remesh/ball_pivot.h"

/* The inward boundary of the supplied F0/fixer mask, not a CT intensity
 * ridge. Coordinates follow the mesher: voxel i is centred at i + 0.5.
 * The volume and axis samples are borrowed for one cube's lifetime. */
typedef struct PredRectoField {
    const uint8_t *vol;
    int d, h, w, halo;
    BpaReconGate axis;
} PredRectoField;

typedef struct PredRectoStats {
    size_t tested, supported, outside, ambiguous, missing, broad, reverted;
    double sum_abs, sum_signed, max_abs;
} PredRectoStats;

/* Independent, compact-support queries give matching targets in a shared
 * halo. Unsupported rays do not borrow another wrap or clamp a distant edge.
 * Returns 1 with signed inward displacement, or 0 with a refusal counted. */
int PredRecto_target(const PredRectoField *field, const float pos[3],
                     const float normal[3], float inward[3], float *distance,
                     PredRectoStats *stats);
void PredRecto_cloud(const PredRectoField *field, float *verts,
                     const float *normals, size_t nv, PredRectoStats *stats);
/* Re-anchor after remeshing. Accept/revert each whole chart atomically if
 * any face would reverse or collapse; never partly clip its displacements. */
void PredRecto_meshes(Arena_T arena, const PredRectoField *field,
                      ComponentMesh *meshes, size_t count, const char *stage);
void PredRecto_report(const char *stage, const PredRectoStats *stats);
int PredRecto_selftest(void);

#endif
