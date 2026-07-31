#ifndef ATLAS_PLACE_SEARCH_INCLUDED
#define ATLAS_PLACE_SEARCH_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"

/* ============================================================================
 * atlas_place_search.h -- put each welded piece on a whole wind, by search.
 *
 * The move itself is settled: a single RIGID translation in u per piece, v
 * untouched.  Rigid is not a compromise, it is the correct move -- measured,
 * a rigid shift reproduces the incoming cube-adjacency numbers to the digit,
 * while the exact phi-dependent shear degrades them 5x because its middle term
 * b*k*phi reads phi, and phi's integer part is a per-cube rounding that was
 * never repaired (622 winding-group nodes, 258 edges, 373 components, 235
 * frustrated, zero corrections applied on the 4x5x5 fixture).
 *
 * What is NOT settled is which wind each piece belongs on.  Both signals that
 * could supply the integer are compromised by the same ~85-vox axis wander:
 * phi is unregistered across cube seams, and radius agrees with it only 56.8%
 * of the time.  So this module does not ask either of them.  It SEARCHES the
 * integer and scores every candidate against the thing we actually want --
 * does this piece land on top of anything? -- so a placement is confirmed by
 * measurement instead of inherited from a rounding.
 *
 * Scoring is a sparse occupancy raster of the atlas in (u, v).  Cells are
 * cheap enough to re-score thousands of candidates, and because v never moves
 * and u moves by a constant, a candidate is just an integer offset along one
 * axis of the key -- no re-rasterization per candidate.  The exact tri-tri
 * audit stays where it belongs, as the verifier of the final answer rather
 * than the inner loop of a search.
 *
 * The overlap constraint is GLOBAL on purpose.  Moving one piece a wind out
 * can clear its own neighbours and land it on a piece from a different support
 * component entirely, so every candidate is scored against the whole atlas,
 * with everything not being moved held fixed as background.
 * ==========================================================================*/

/*
 * One radial-adjacency assertion, in the only form the search needs: these two
 * pieces must end up `target` apart in u, and in the incoming field they are
 * `du_base` apart.  piece_a/piece_b of -1 mean an end held fixed as background.
 *
 * This is what stops the search from merely PACKING.  Collision-freedom alone
 * underdetermines the answer -- any permutation that tiles without overlap
 * scores perfectly -- and measured on support 3 the collision-only objective
 * picked an order the radial evidence likes markedly less than phi's (3,031
 * pairs correctly separated against 5,001).  The radial pairs are 86%
 * intra-cube, so unlike phi's integer they are not poisoned by the unrepaired
 * cross-cube registration, which makes them the right second term.
 */
typedef struct {
    int32_t piece_a;   /* inner */
    int32_t piece_b;   /* outer */
    double  du_base;   /* u_outer - u_inner in the incoming field */
    double  target;    /* signed u the pair asserts for (outer - inner) */
} AtlasPlacePair;

typedef struct {
    double cell;         /* raster cell size in u and v, vox (4.0) */
    double lambda_move;  /* cells of collision one wind of travel is worth;
                          * keeps the search from spreading pieces further
                          * than the overlap actually requires (8.0) */
    double lambda_pair;  /* cells of collision one unsatisfied radial pair is
                          * worth (1.0) */
    int    rounds;       /* coordinate-descent sweeps (6) */
} AtlasPlaceSearchOptions;

void AtlasPlaceSearchOptions_default(AtlasPlaceSearchOptions *opts);

typedef struct {
    size_t pieces;
    size_t candidates_evaluated;
    size_t static_cells;
    size_t moving_cells;
    size_t collisions_before;
    size_t collisions_after;
    size_t collisions_static_after;  /* against the rest of the atlas */
    size_t moved_pieces;
    int    rounds_run;
    int    converged;                /* a full sweep changed nothing */

    size_t pairs;
    size_t pairs_intra;              /* both ends in one piece: unreachable */
    size_t pairs_satisfied_before;
    size_t pairs_satisfied_after;
} AtlasPlaceSearchStats;

/*
 * u, v            [nv] the incoming atlas
 * vertex_piece    [nv] index of the piece a vertex belongs to, or -1 to hold
 *                 it fixed as background
 * npieces, ncand  the candidate table is [npieces * ncand]
 * shift           [npieces * ncand] u translation of candidate j for piece p
 * cost            [npieces * ncand] travel cost of that candidate, in winds
 *                 (|k|); scaled by lambda_move and added to the collisions
 * out_choice      [npieces] chosen candidate index
 */
int AtlasPlaceSearch_solve(Arena_T arena,
                           const double *u, const double *v, size_t nv,
                           const int32_t *vertex_piece,
                           size_t npieces, size_t ncand,
                           const double *shift, const double *cost,
                           const AtlasPlacePair *pairs, size_t npairs,
                           const AtlasPlaceSearchOptions *opts,
                           int32_t *out_choice,
                           AtlasPlaceSearchStats *stats);

int AtlasPlaceSearch_selftest(void);

#endif
