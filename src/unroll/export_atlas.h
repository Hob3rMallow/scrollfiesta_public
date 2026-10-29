#ifndef EXPORT_ATLAS_INCLUDED
#define EXPORT_ATLAS_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"
#include "piece_set.h"
#include "scaffold.h"

/* ============================================================================
 * export_atlas.h -- L3 EXPORT ATLAS. Split the registered piece set into
 * (wrap k x z-slab) pieces by scaffold labels and export each as its own VC3D
 * tifxyz segment on a shared du/dv lattice, plus an atlas.json manifest.
 *
 * Each piece is a contiguous band of the registered strip: one wrap (or a few)
 * is a contiguous u-band, one z-slab is a contiguous v-band, so pieces tile
 * the strip disjointly BY PHYSICS -- no UV packing. Every segment records its
 * origin_uv, and because canvas origins are floor(u/du+0.5) all pieces share
 * one integer lattice, so adjacent segments' pixels co-register.
 *
 * v1 exports the REGISTERED UVs directly (the decoupling: export UV is
 * whatever is in the piece set, independent of repair history). The
 * ribbon_relax polish (opts.relax) is the A-3 quality pass.
 * ==========================================================================*/

typedef struct {
    int sign;
    int curved_axis;
    double offset;
    double coherence;
    size_t groups, ambiguous_groups;
    size_t cycle_conflicts, isolated_vertices;
} AtlasPolarStats;

typedef struct {
    int    piece_wraps;    /* wraps (k) per piece (default 1) */
    double slab_v;         /* z-slab height in vox; >= v-extent => one slab (default 4096) */
    double du, dv;         /* segment grid step, vox/px (default 1) */
    int    write_winding;  /* emit winding.tif per segment (default 1 in atlas mode) */
    /* face gates (tifxyz defaults; real/source max_edge3d is disabled by
     * default because valid edge scale follows the upstream remesher) */
    double max_edge3d, stretch_ratio, stretch_floor, conflict_dist;
    /* A-3 polish: per-piece symmetric-Dirichlet relax of the piece's UVs,
     * pinning the rim to the shared frame. v1 default 0 (export raw). */
    int    relax;
    int    relax_sweeps;
    int    phase_u;        /* export u=F(phi), retaining registered v */
    int    polar_u;        /* lifted polar angle, registered phi supplies turn */
    int    input_polar_u;  /* input PieceSet already uses the lifted chart */
    AtlasPolarStats input_polar_stats; /* provenance for that prebuilt chart */
    int    verbose;
} AtlasOpts;

void AtlasOpts_default(AtlasOpts *o);

/* Build a smooth covering-space coordinate: local polar angle supplies the
 * differential phase, while registered phi selects its integer 2pi lift.
 * Returned uv is arena-owned; v is copied from the input piece set. */
int Atlas_lift_polar_u(Arena_T arena, const PieceSet *ps,
                       const ScaffoldCalib *c, float **out_uv,
                       AtlasPolarStats *stats);

typedef struct {
    size_t n_pieces;          /* pieces with >=1 kept face */
    size_t n_written;         /* segments successfully written */
    size_t n_empty;           /* (k,slab) cells with no faces (skipped) */
    size_t n_over_cap;        /* pieces that exceeded the TIFF band cap */
    size_t n_quarantine_faces;/* total exclusions below (legacy aggregate) */
    size_t n_phase_quarantine_faces; /* faces spanning >=2 wraps */
    size_t n_u_outlier_faces; /* faces outside the robust per-piece u window */
    size_t total_valid_px, total_multi_px, total_conflict_px;
    size_t total_conflict_same_cube_px, total_conflict_cross_cube_px;
    size_t total_conflict_unknown_cube_px;
    size_t total_conflict_d_le4, total_conflict_d_le8;
    size_t total_conflict_d_le16, total_conflict_d_le32;
    size_t total_conflict_d_gt32;
    double total_conflict_d_sum, max_conflict_d;
    size_t total_conflict_turn_known, total_conflict_turn_le025;
    size_t total_conflict_turn_le05, total_conflict_turn_le1;
    size_t total_conflict_turn_gt1;
    size_t total_conflict_face_shared_edge;
    size_t total_conflict_face_shared_vertex, total_conflict_face_disjoint;
    size_t total_conflict_centroid_d_le8, total_conflict_centroid_d_gt8;
    size_t total_conflict_centroid_du_le2, total_conflict_centroid_du_le4;
    size_t total_conflict_centroid_du_le8, total_conflict_centroid_du_gt8;
    size_t total_conflict_centroid_dv_le2, total_conflict_centroid_dv_le4;
    size_t total_conflict_centroid_dv_le8, total_conflict_centroid_dv_gt8;
    double conflict_frac;     /* total_conflict_px / total_valid_px */
    int    polar_sign, polar_curved_axis;
    double polar_offset, polar_coherence;
    size_t polar_groups, polar_ambiguous_groups;
    size_t polar_cycle_conflicts, polar_isolated_vertices;
    double seconds;
} AtlasStats;

/* Enumerate pieces from ps's scaffold (calib->sense gives k = floor(sense*phi/
 * 2pi); v == world z), export each to seg_root/<prefix>_w<k>_z<v>/, and write
 * seg_root/atlas.json. opts NULL -> defaults. Returns 0, -1 on bad input. */
int ExportAtlas_run(Arena_T arena, const PieceSet *ps, const ScaffoldCalib *c,
                    const char *seg_root, const char *prefix,
                    const AtlasOpts *opts, AtlasStats *out);

int ExportAtlas_selftest(void);

#endif
