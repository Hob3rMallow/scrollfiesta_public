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
    int    verbose;
} AtlasOpts;

void AtlasOpts_default(AtlasOpts *o);

typedef struct {
    size_t n_pieces;          /* pieces with >=1 kept face */
    size_t n_written;         /* segments successfully written */
    size_t n_empty;           /* (k,slab) cells with no faces (skipped) */
    size_t n_over_cap;        /* pieces that exceeded the TIFF band cap */
    size_t n_quarantine_faces;/* faces spanning >=2 wraps (contested; excluded) */
    size_t total_valid_px, total_multi_px, total_conflict_px;
    double conflict_frac;     /* total_conflict_px / total_valid_px */
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
