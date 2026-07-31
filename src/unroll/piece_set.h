#ifndef PIECE_SET_INCLUDED
#define PIECE_SET_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"

/* ============================================================================
 * piece_set.h -- concatenate a scroll_whole placed dir (per-cube registered
 * meshes + UV sidecars) into one flat renumbered mesh for the global strip
 * raster. Faces the front half dropped (bad fusion links, unmapped verts) are
 * excluded here; vertices are all kept so sidecar indices stay aligned.
 * Per-vertex normals (for the normal-max RAW sampling) are computed per cube
 * over its kept faces.
 * ==========================================================================*/

typedef struct {
    float   *verts;      /* [nv*3] world (z,y,x) */
    float   *uv;         /* [nv*2] registered global (u,v), vox */
    float   *phi;        /* [nv] registered winding angle (third channel of
                            _uvphi.f32) -- SeamOwn's rehome needs it for the
                            exact DeltaU relocation */
    float   *normals;    /* [nv*3] unit (zero where undefined) */
    int32_t *faces;      /* [nf*3] kept faces, renumbered into the union */
    int32_t *face_cube;  /* [nf] index into ids[] (provenance) */
    int32_t *gid;        /* [nv] per-cube winding-group id from _group.i32;
                            -1 = unmapped or sidecar absent. NOT comparable
                            across cubes -- same-cube grouping only. */
    size_t   nv, nf;
    char   (*ids)[48];   /* [n_cubes] cube ids in load order */
    size_t  *cube_voff;  /* [n_cubes+1] vertex ranges (cubes load contiguous) */
    long   (*cube_org)[3];  /* [n_cubes] (oz,oy,ox) voxel origin from the id */
    size_t   n_cubes;
    double   u_min, u_max, v_min, v_max;   /* over verts of kept faces */
} PieceSet;

/* Load every <id>_uvphi.f32 + <id>_mesh.obj + <id>_facekeep.u8 (+ optional
 * <id>_group.i32) under placed_dir (the scroll_whole output). All outputs
 * arena-allocated. Returns 0; -1 if the dir has no complete placed cubes. */
int PieceSet_build(Arena_T arena, const char *placed_dir, PieceSet *out);

/* Same, but only cubes whose z-origin (parsed from the id, "z#####_...") lies
 * in [z_lo, z_hi) are loaded -- the z-slab / strip selector. Cubes share the
 * one global solve (same placed dir, same calibration), so strips exported
 * this way co-register on the u lattice and stitch in v (world z) trivially.
 * z_lo = LONG_MIN, z_hi = LONG_MAX loads everything (== PieceSet_build). */
int PieceSet_build_z(Arena_T arena, const char *placed_dir,
                     long z_lo, long z_hi, PieceSet *out);

/* Compact faces in place to keep[f] != 0, keeping face_cube aligned (the
 * step2 SeamOwn interop: downstream stages then see kept-only faces with no
 * mask plumbing). Vertices are untouched so sidecar indices stay aligned.
 * Returns the new nf (also written back to ps->nf). */
size_t PieceSet_apply_facekeep(PieceSet *ps, const uint8_t *keep);

/* Recompute per-vertex normals over the CURRENT faces (verts never shared
 * across cubes, so this stays per-cube). Call after apply_facekeep or after
 * a stage mutates verts. */
void PieceSet_refresh_normals(PieceSet *ps);

/* In-process unit tests (tiny synthetic placed dir under
 * output/_selftest_scroll_unroll/). Returns 0 if all pass. */
int PieceSet_selftest(void);

#endif
