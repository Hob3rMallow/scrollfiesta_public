#ifndef PATCH_REPAIR_INCLUDED
#define PATCH_REPAIR_INCLUDED

#include <stddef.h>
#include <stdint.h>
#include "../common/arena.h"

/*
 * PatchRepair — local tangent-plane Delaunay re-triangulation of BPA tear
 * patches ("lightning-bolt" slits).
 *
 * Diagnosis this addresses (changelog 2026-07-09): BPA's empty-ball test at
 * rho=1.2 on ~1.0-vox in-plane spacing dead-ends in patches (gap=0, blockers
 * COPLANAR, |height| ~ 0.001 vox), the front self-collides, and the mesh is
 * left with small jagged interior boundary loops, often self-pinching. The
 * points are fine; the triangulation is torn. Global rho changes were tried
 * and REVERTED (they re-scale the whole edge economy the weld is calibrated
 * to, and connected wrong sheets). This module repairs LOCALLY instead:
 *
 *   per connected sheet:                          (sheet = face-connectivity
 *     find interior tear clusters                  component; repairs are an
 *       (boundary edges NOT on the perimeter)      individual-sheet operation)
 *     grow a K-ring face patch around each; merge overlapping patches
 *     gate: small, planar, nowhere near the sheet perimeter, single simple ring
 *     delete the patch faces; re-triangulate the ring interior in the patch's
 *       tangent plane: ear-clip the ring + insert the patch's interior verts +
 *       Lawson flips to (constrained) Delaunay
 *     orient new faces to the deleted faces' mean normal AND verify every ring
 *       edge is used exactly once, in reverse of its kept outside face
 *     commit only if every check passes; otherwise the patch is left untouched
 *
 * Topology invariants (by construction, not by tolerance):
 *   - Vertices are never moved, added, or removed.
 *   - A patch is a connected face-neighbourhood of ONE component; the
 *     triangulation sees only that patch's vertices. No proximity queries.
 *     New edges connect only vertices of the same patch => the repair can
 *     NEVER connect two components (sheets), and never bridges within a sheet
 *     beyond the patch ring.
 *   - Patches whose ring is not one simple cycle, that touch the sheet
 *     perimeter, that fail the planarity gate, or whose triangulation fails
 *     any validity check are skipped (mesh unchanged there).
 *
 * Output faces are allocated on `arena`; per the Weld_verts convention the
 * scratch is NOT restored (caller's arena lifetime subsumes the call).
 * Deterministic for identical input (index-ordered iteration throughout).
 */

typedef struct PatchRepairStats {
    size_t n_components;      /* face-connectivity components seen            */
    size_t n_tear_clusters;   /* non-perimeter boundary clusters found        */
    size_t n_patches;         /* patches after overlap merging                */
    size_t n_repaired;        /* patches committed                            */
    size_t n_skip_gate;       /* skipped: size/perimeter/planarity gates      */
    size_t n_skip_ring;       /* skipped: ring not a single simple cycle      */
    size_t n_skip_project;    /* skipped: 2D projection degeneracies          */
    size_t n_skip_tri;        /* skipped: triangulation/orientation invalid   */
    size_t faces_deleted;     /* faces removed by committed repairs           */
    size_t faces_added;       /* faces added by committed repairs             */
} PatchRepairStats;

/* Re-triangulate tear patches of (verts, faces). On success *out_faces /
 * *out_nf hold the repaired face list (arena-allocated; == input copy when
 * nothing was repairable). Returns 0 on success (including "nothing to do"),
 * non-zero only on contract violations. `stats` may be NULL. */
int PatchRepair_repair(Arena_T arena,
                       const float *verts, size_t nv,
                       const int32_t *faces, size_t nf,
                       int32_t **out_faces, size_t *out_nf,
                       PatchRepairStats *stats);

#endif
