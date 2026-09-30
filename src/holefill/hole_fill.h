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

/* Legacy displacement allowance for lineage readers of older hole-fill
 * artifacts. Current chart filling retains source coordinates and vertex IDs;
 * an uncertifiable cycle stays open. */
#define HOLEFILL_CHART_PINCH_MAX_STEP_VOX 0.25f

/* Chart fills must not create a closed source fan that already makes the
 * assembler's 0.75..1.25 stretch band impossible. This is a necessary angle
 * bound, independent of triangulation density, coordinate axes, and UV seeds.
 * An incompatible opening stays open; original material is never discarded
 * merely to make a generated cap pass. */
#define HOLEFILL_CHART_MAX_CONDITION (1.25 / 0.75)

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
 * HoleFill_process_ex — same as HoleFill_process, with a classification mode.
 *
 *   interior_only == 1 : fill ONLY geometrically-interior holes (loops the
 *                        surface surrounds), determined by a signed-area /
 *                        winding test relative to the averaged incident-face
 *                        normal. Outer perimeters and still-open boundary bays
 *                        are left untouched. Correct for multi-component meshes
 *                        (each component's own perimeter is recognised), so it
 *                        is the mode the post-weld joined mesh wants.
 *   interior_only == 2 : chart mode.  Treat each connected component as one
 *                        source chart, preserve that chart's longest boundary
 *                        loop as its perimeter, and fill every other cleanly-
 *                        fillable loop. New interior fans must also pass the
 *                        necessary bound HOLEFILL_CHART_MAX_CONDITION. This is
 *                        the post-weld mode when chart
 *                        provenance is authoritative and avoids misclassifying
 *                        steep 3-D punctures by projected winding.
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
