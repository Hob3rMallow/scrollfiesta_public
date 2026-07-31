#ifndef WELD_TRACK_INCLUDED
#define WELD_TRACK_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"

/* ============================================================================
 * weld_track.h -- the tracked-weld dataset: cross-cube boundary CORRESPONDENCES
 * (Pass 1), persisted as one versioned binary sidecar keyed to the placed set
 * by (cube-id string, local skin-vertex index). Pass-2 bridge/fill GEOMETRY is
 * a later addition to the same file.
 *
 * "weld everything at full resolution but track the welds in their own dataset,
 * patch them in as needed": a correspondence asserts "these two boundary verts
 * are the SAME physical papyrus at the same turn." It is a CONSTRAINT (winding
 * + arc-length continuity) and the connectivity backbone for reassembly -- it
 * is NEVER a destructive merge of the geometry the solvers see.
 *
 * seam_axis records the cube-adjacency direction: 0 = z (ALONG the scroll axis,
 * winding-trivial -- these are the ministrip-stacking seams), 1 = y, 2 = x (the
 * RADIAL plane, where winding lives). Winding consumes the y/x correspondences
 * per z-layer; z correspondences drive ministrip reassembly (weld + hole-fill +
 * arc-length continuity).
 * ==========================================================================*/

#define WELD_TRACK_MAGIC   "WTRK"
#define WELD_TRACK_VERSION 1u
#define WTRK_ID_LEN        48

/* one cross-cube boundary correspondence (40 bytes, no compiler padding). */
typedef struct {
    uint32_t cube_a, cube_b;   /* index into WeldTrack.ids (u32: full scrolls
                                  exceed 65k cubes) */
    uint32_t vert_a, vert_b;   /* LOCAL skin-vertex index within that cube */
    int32_t  gid_a, gid_b;     /* per-cube winding-group id (NOT cross-comparable) */
    float    dr;               /* |r_a - r_b| about the axis, vox */
    float    dphi;             /* wrap_pi(phi_b - phi_a), rad (near 0 = same wrap) */
    float    dist3d;           /* 3D separation, vox */
    uint8_t  seam_axis;        /* 0 = z, 1 = y, 2 = x */
    uint8_t  conf;             /* 2 = mutual-nearest + margin, 1 = mutual only */
    uint8_t  _pad[2];
} WtrkCorr;

typedef struct WeldTrack {
    uint32_t version;
    double   spiral_a, spiral_b, pitch;   /* DeltaU params (from placed_index) */
    float    axis_point[3];               /* umbilicus (z,y,x) */
    size_t   n_cubes;
    char   (*ids)[WTRK_ID_LEN];           /* [n_cubes] cube-id strings */
    size_t   n_corr;
    WtrkCorr *corr;                       /* [n_corr] */
    void    *self;                        /* validation sentinel: self == this */
} WeldTrack;

/* Binary sidecar IO (little-endian; the codebase assumes LE). Read allocates
 * ids/corr from arena. Return 0; -1 on IO / format / version error. */
int WeldTrack_write(const WeldTrack *w, const char *path);
int WeldTrack_read(Arena_T arena, const char *path, WeldTrack *out);

/* cube-id string -> table index, or -1 if absent. */
int32_t WeldTrack_cube_index(const WeldTrack *w, const char *id);

/* In-process unit tests (roundtrip write/read, cube lookup). Writes scratch
 * under output/_selftest_weld_track/. Returns 0 if all pass, else failures. */
int WeldTrack_selftest(void);

#endif
