#ifndef CUBE_SCHEDULE_INCLUDED
#define CUBE_SCHEDULE_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"

/* ============================================================================
 * cube_schedule.h -- discover the per-cube dump grid and order it for
 * cross-seam registration: an outgoing spiral from a clean mid-shell seed.
 *
 * Registration (cube_register.h) corrects each cube's integer winding against
 * its ALREADY-PLACED neighbors, so the order must (a) start where the radius
 * anchor is reliable (mid-shell, not the chaotic near-umbilicus core), and
 * (b) only ever visit a cube adjacent to placed ones. Best-first flood:
 * priority = |w_phase - w_phase(seed)| ascending, where w_phase =
 * r/pitch - theta/2pi is the continuous winding phase of the cube center --
 * an annulus-by-annulus spiral sweep outward AND inward from the seed band,
 * reaching the core last (it then registers against many placed neighbors).
 * ==========================================================================*/

typedef struct {
    char     id[48];       /* z#####_y#####_x##### */
    int64_t  oz, oy, ox;   /* origin (source-space voxel) */
    double   r, theta;     /* cube-center polar coords about the axis */
    double   w_phase;      /* r/pitch - theta/(2*pi) at the cube center */
    int32_t  nbr[6];       /* 6-neighbor node index (-z,+z,-y,+y,-x,+x); -1 = absent */
} CubeNode;

/* Scan the z*_y*_x* dirs under dump_dir, keeping cubes that have
 * <id>/<id>_<leaf_stage>/<id>_<leaf_stage>_all.obj, fill CubeNode geometry
 * about axis_point (z,y,x order, axis = +z), link 6-neighbors (origins one
 * chunk apart), and emit the registration order. seed_id (may be NULL) forces
 * the seed cube; the default seed is the present cube with the MEDIAN center
 * radius. Disconnected grid components are re-seeded and tagged in out_comp
 * (0 = main flood). pitch <= 0 requests a pitch-free provisional radial order;
 * callers can calibrate a pitch from that seed and call link_and_order again.
 * Outputs arena-allocated. Returns 0; -1 if no cubes. */
int CubeSched_build(Arena_T arena, const char *dump_dir, const char *leaf_stage,
                    int64_t chunk, const float axis_point[3], double pitch,
                    const char *seed_id,
                    CubeNode **out_nodes, size_t *out_n,
                    int32_t **out_order, int32_t **out_comp);

/* Link + order over an already-filled node array (id/oz/oy/ox set; r/theta/
 * w_phase computed here too). Exposed so the selftest and audit can run the
 * pure logic without a dump dir on disk. seed_idx < 0 = median-radius seed.
 * pitch <= 0 uses |r-r(seed)| priority without assuming a wrap spacing. */
int CubeSched_link_and_order(Arena_T arena, CubeNode *nodes, size_t n,
                             int64_t chunk, const float axis_point[3],
                             double pitch, int32_t seed_idx,
                             int32_t **out_order, int32_t **out_comp);

/* In-process unit tests (grid link, flood order invariants, hole, island).
 * Returns 0 if all pass, else the number of failures. Logs to stderr. */
int CubeSched_selftest(void);

#endif
