#ifndef QUADRIBBON_UNTANGLE_INCLUDED
#define QUADRIBBON_UNTANGLE_INCLUDED

#include <stddef.h>
#include <stdint.h>

/* Collision untangler for a parameterized quadribbon, extracted from
 * solid_quad_ribbon.c's structured-lattice shell (turn-order repair +
 * elastic-shell contact solve + settle-toward-rest) and generalized to a
 * plain triangle mesh with UV: the lattice quad-cell becomes an occupied
 * UV bucket (~one column step wide, one row pitch tall), the grid
 * 4-neighbourhood becomes the mesh vertex adjacency, and red/black
 * Gauss-Seidel becomes two-buffer Jacobi with doubled sweeps.  Everything
 * that gates acceptance is unchanged: the exact BVH audit
 * (IntersectionCleanup_audit_visit_parallel), the baseline-relative
 * orientation preflight, the trust-region alpha ladder with feathered local
 * rollbacks, phase-violation progress under bounded pair churn, the
 * persistent contact ledger with return-path adoption, and the
 * lexicographically-best partial publication.
 *
 * Motion is RADIAL ONLY about the scroll axis (axis parallel to +Z in zyx
 * coordinates, the house frame): one scalar displacement per vertex from the
 * rest geometry.  Topology and UV are never modified -- re-run the metric
 * parameterization afterwards to restore honest arc length on the moved
 * geometry. */

typedef struct QuadribbonUntangleOpts {
    int require_collision_free;   /* 1: roll back unless fully untangled */
    int collision_patience;       /* accepted rounds without a new best */
    int collision_rounds;         /* active-set round cap */
    double collision_collar;      /* speculative contact collar, vox */
    double displacement_bound;    /* |radial displacement| cap, vox */
    int settle_rounds;            /* settle-toward-rest rounds (0 = off) */
    double settle_beta;           /* per-round pull toward the rest coil */
    double minimum_u_separation;  /* long-range conflict threshold; <=0 =>
                                   * max(100, 64 * median in-row u step) */
    double wrap_pitch_hint;       /* |pitch| prior (vox) used only when the
                                   * measured inference is ambiguous; sign is
                                   * always measured.  <=0 => none */
} QuadribbonUntangleOpts;

typedef struct QuadribbonUntangleStats {
    int attempted, accepted_rounds, complete, retained_partial;
    size_t input_conflicts, input_long_conflicts;
    size_t output_conflicts, output_long_conflicts;
    size_t contacts_built, hard_contacts;
    double pitch, clearance, movement_rms, movement_max;
    int settle_rounds_run, settle_accepted;
    double settle_rms_before, settle_rms_after;
    int turn_order_rounds;
    size_t turn_order_input_long, turn_order_output_long;
} QuadribbonUntangleStats;

void QuadribbonUntangle_defaults(QuadribbonUntangleOpts *o);

/* Untangles verts in place (movable == NULL means every vertex may move).
 * Returns 0 on success (possibly a retained partial state; see stats),
 * 1 when require_collision_free rolled the transaction back, -1 on error. */
int QuadribbonUntangle_run(float *verts, size_t nv,
                           const int32_t *faces, size_t nf,
                           const float *uv, const uint8_t *movable,
                           double axis_y, double axis_x,
                           const QuadribbonUntangleOpts *opts,
                           QuadribbonUntangleStats *out);

int QuadribbonUntangle_selftest(void);

#endif
