/*
 * hole_fill.h — Interior hole detection and filling.
 *
 * Finds boundary loops, classifies interior vs. exterior holes,
 * fills interior holes with CDT + cotangent Laplacian CG smoothing.
 */
#ifndef HOLE_FILL_INCLUDED
#define HOLE_FILL_INCLUDED

#include "../common/arena.h"
#include <stdint.h>
#include <stddef.h>

/*
 * HoleFill_process — detect and fill interior holes in a triangle mesh.
 *
 * Parameters:
 *   arena         Arena for all allocations
 *   verts         Pointer to vertex array [nv*3], may be reallocated
 *   faces         Pointer to face array [nf*3], may be reallocated
 *   nv            Pointer to vertex count, updated on output
 *   nf            Pointer to face count, updated on output
 *   cube_shape    Volume dimensions [D, H, W] for interior/exterior classification
 *   out_n_loops   (optional) Total boundary-loop count
 *   out_n_interior (optional) Interior-loop count (fill candidates)
 *   out_n_filled  (optional) Successfully-filled count
 *
 * Returns 0 on success, -1 on failure.
 */
int HoleFill_process(Arena_T arena,
                     float **verts, int32_t **faces,
                     size_t *nv, size_t *nf,
                     const int cube_shape[3],
                     size_t *out_n_loops,
                     size_t *out_n_interior,
                     size_t *out_n_filled);

/*
 * HoleFill_process_ex — same as HoleFill_process, with an interior_only switch.
 *
 *   interior_only != 0 : fill ONLY geometrically-interior holes (loops the
 *                        surface surrounds), determined by a signed-area /
 *                        winding test relative to the averaged incident-face
 *                        normal. Outer perimeters and still-open boundary bays
 *                        are left untouched. Correct for multi-component meshes
 *                        (each component's own perimeter is recognised), so it
 *                        is the mode the post-weld joined mesh wants.
 *   interior_only == 0 : fill every cleanly-fillable closed 4+ loop
 *                        (HoleFill_process wraps this).
 *
 * cube_shape is unused (the interior test is geometric); pass NULL.
 */
int HoleFill_process_ex(Arena_T arena,
                        float **verts, int32_t **faces,
                        size_t *nv, size_t *nf,
                        const int cube_shape[3],
                        int interior_only,
                        size_t *out_n_loops,
                        size_t *out_n_interior,
                        size_t *out_n_filled);

#endif /* HOLE_FILL_INCLUDED */
