#ifndef PINHOLE_FILL_INCLUDED
#define PINHOLE_FILL_INCLUDED

#include <stddef.h>
#include "../common/arena.h"
#include "../common/mesh_types.h"

/*
 * Dedicated micro-hole closer for Ball-Pivoting output. SEPARATE from
 * src/holefill (the CDT path), which is tuned for large holes and degenerates
 * at 3-vertex loops (it inserts Steiner points / slivers). The BPA surface
 * instead has many tiny gaps the roll skipped: ~one missing triangle per gap,
 * plus a handful of bowtie ("pinch") vertices where two sheet corners touch.
 *
 * Run this AFTER orientation (so boundary half-edges are consistently wound)
 * and BEFORE HoleFill_process (so the CDT path only ever sees genuine large
 * holes). Two phases per component:
 *
 *  1. Split non-manifold vertices. A vertex whose incident faces fall into
 *     >=2 edge-connected fans is a bowtie; each extra fan gets its own copy
 *     of the vertex (position / normal / pin status duplicated) and that
 *     fan's faces are repointed to it. This makes the mesh vertex-manifold
 *     and, crucially, makes the open boundary a set of simple cycles so the
 *     loop tracer is unambiguous.
 *
 *  2. Fill EXACT single-triangle (3-loop) holes. Boundary half-edges (edges in
 *     exactly one face) are chained by direction into closed loops; a 3-loop is
 *     closed with its one triangle, wound OPPOSITE the boundary half-edges so
 *     the patch is consistent with the surrounding faces — no PCA, no CDT, no
 *     Steiner points, no smoothing. A larger loop (4+) is a real hole left to
 *     HoleFill_process (CDT). A 3-loop is filled only if (i) it fits inside a
 *     PINHOLE_MAX_DIAM_VOX ball (no-merger: never bridge the wide gap between
 *     two scroll wraps), (ii) its face-component owns at least one OTHER
 *     boundary loop (a tiny island's sole rim is never capped into a bubble),
 *     and (iii) its triangle is not a degenerate collinear sliver (min altitude
 *     >= HOLEFILL_MIN_ALT_VOX) — the "chewed" boundary near a cube face is a row
 *     of collinear verts whose fill would be zero-area. Loops failing any gate,
 *     or larger than 3 verts, are left for HoleFill_process.
 *
 * If `respect_pins` is non-zero and a component has a pin_mask, any loop that
 * touches a pinned (halo seam) vertex is left unfilled — the cross-cube seam
 * weld owns that geometry.
 *
 * Reallocates each component's verts / vert_normals / pin_mask (phase 1 may
 * add verts) and faces (phase 2 adds faces) from `arena`. Counts are optional.
 *
 * Returns 0 on success.
 */
int PinholeFill_process(Arena_T arena,
                        ComponentMesh *meshes, size_t n_meshes,
                        int respect_pins,
                        size_t *out_pinch_splits,
                        size_t *out_loops_filled,
                        size_t *out_tris_added,
                        size_t *out_loops_skipped);

/* Phase-1 only: split pinch (bowtie) vertices in each component, no hole filling.
 * The vertex-manifold half of PinholeFill, exposed for the post-trim manifold
 * guard (a trimmed mesh has no pins and needs no 3-loop fills, only de-pinching).
 * Requires EDGE-manifold input (resolve >2-face edges first). Returns 0 on
 * success; out_splits (optional) = total vertices duplicated. */
int PinholeFill_split_pinches(Arena_T arena,
                              ComponentMesh *meshes, size_t n_meshes,
                              size_t *out_splits);

/* Per-mesh CDT/Liepa fill for the 4+ loops PinholeFill leaves (returns 0 on
 * success). The pipeline implements it over src/holefill (the only TU that
 * links Triangle/Clipper2); grid_weld passes NULL, so the weld binary stays
 * CDT-free and 4+ loops are left open for the seam weld. */
typedef int (*HoleFill_cdt_fn)(Arena_T arena, ComponentMesh *cm, void *user);

/* Unified hole-fill entry: the pinch-split + EXACT 3-loop fills of
 * PinholeFill_process (the CDT-free core), then every remaining 4+ closed loop
 * is handed to cdt_fill (NULL => left open). One call closes all per-cube holes,
 * dispatched purely by loop size — no interior/exterior classification. */
int HoleFill_meshes(Arena_T arena,
                    ComponentMesh *meshes, size_t n_meshes,
                    int respect_pins,
                    HoleFill_cdt_fn cdt_fill, void *cdt_user,
                    size_t *out_pinch_splits,
                    size_t *out_loops_filled,
                    size_t *out_tris_added,
                    size_t *out_loops_skipped,
                    size_t *out_cdt_filled);

#define PINHOLE_MAX_LOOP 12   /* largest boundary loop the loop tracer follows */

/* A loop is only filled if its vertices fit inside a ball of this diameter
 * (voxels). This is the structural no-merger guard: a genuine BPA pinhole is a
 * few missing triangles spanning a couple of voxels, whereas the gap between
 * two scroll wraps is far wider. Refusing wider loops means PinholeFill can
 * never bridge two wraps, regardless of how few vertices the loop has. */
#define PINHOLE_MAX_DIAM_VOX 4.5f

/* A single-triangle (3-loop) fill is rejected if its min altitude (2*area /
 * longest edge, voxels) is below this — i.e. the triangle is a near-collinear
 * zero-area sliver. The "chewed" ragged boundary along a cube face is a row of
 * ~1-vox-spaced collinear verts; filling such a 3-loop yields a degenerate
 * triangle overlapping its own boundary edge (observed altitude ~0.002 vox)
 * that corrupts the weld. Healthy fills are ~0.4 vox; leave the notch open. */
#define HOLEFILL_MIN_ALT_VOX 0.05f

/* PinholeFill phase 0 (close_bowtie_gaps): resolve a non-manifold "bowtie"
 * vertex IN PLACE -- keep the single vertex and fan-triangulate the angular
 * gaps between its fans -- instead of duplicating it (split_pinch_verts), which
 * on a near-coincident cross-cube seam spawns a doubled, z-fighting fill. Only
 * applied when the fans form one well-defined surface: every incident face's
 * normal agrees with the area-weighted vertex normal to within
 * BOWTIE_COHERENCE_COS, AND each inter-fan gap is narrower than
 * BOWTIE_GAP_MAX_VOX. The gap cap sits well below the inter-wrap clearance
 * (~CUT_GAP_DEPTH=7) so two wraps touching at a point are NEVER welded -- they
 * fail the gap test and fall through to the split. */
#define BOWTIE_COHERENCE_COS 0.2f
#define BOWTIE_GAP_MAX_VOX   2.0f

#endif
