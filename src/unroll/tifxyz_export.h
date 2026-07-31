#ifndef TIFXYZ_EXPORT_INCLUDED
#define TIFXYZ_EXPORT_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"
#include "piece_set.h"

/* ============================================================================
 * tifxyz_export.h -- write the unrolled strip as a VC3D tifxyz segment.
 *
 * A tifxyz segment is a directory of three same-size float32 TIFFs
 * (x.tif / y.tif / z.tif) where pixel (row, col) holds one surface sample's
 * 3D position in ABSOLUTE full-resolution voxel coordinates of the scroll
 * volume, plus meta.json (VC3D's loader requires "scale" = grid px per voxel
 * = [1/du, 1/dv]) and mask.tif (uint8; VC3D keeps a pixel only when the mask
 * value is >= 255). Invalid samples are (-1,-1,-1) in all three bands; the
 * VC3D loader additionally drops any sample with z <= 0.
 *
 * Grid convention: cols = u (arc length around the scroll), rows = v (world
 * z along the scroll axis) -- the same canvas mapping as scroll_raster, so
 * tifxyz pixel (row, col) lines up 1:1 with the stage rawtex/diagclass TIFs
 * when du/dv match. Axis order swaps on write: our verts are native (z,y,x),
 * their bands are x.tif <- vert[2], y.tif <- vert[1], z.tif <- vert[0].
 *
 * Coverage is first-cover-wins in ascending face order (never a positional
 * average: a multi-cover pixel whose covers straddle two wraps must not
 * blend into mid-air). Later covers only count multi/conflict stats.
 * ==========================================================================*/

typedef struct {
    double du, dv;            /* vox per px (grid step, default 1.0);
                                 meta "scale" = [1/du, 1/dv] */
    double stretch_ratio;     /* face gates, same semantics + defaults as
                                 RasterOpts: uv edge > ratio*3d ... */
    double stretch_floor;     /*   ... and > floor px  => skip face (uv) */
    double max_edge3d;        /* REAL/source face 3D-edge cap; default 0
                                 (disabled because mesh resolution varies) */
    double conflict_dist;     /* a later cover farther than this (vox) from
                                 the stored position marks the px conflicted
                                 (default 2.0) */
    int    flip_u;            /* mirror columns (u runs right-to-left) */
    int    flip_v;            /* mirror rows */
    int    write_winding;     /* emit winding.tif = phi / 2pi (default 0) */
    int    write_provenance;  /* emit provenance.tif = per-px cover state
                                 0=empty 1=single 2=overlap 3=contested
                                 (winding collapse: two covers > conflict_dist
                                 apart in 3D). Non-destructive: mask stays 255
                                 on contested px. Default 1. */
    const uint8_t *face_keep; /* [nf] 0 = face dropped (NULL = keep all) */
    /* Forced canvas window (vox). When have_window != 0 the canvas spans
     * [win_u0,win_u1] x [win_v0,win_v1] instead of the whole piece set's
     * uv range -- the atlas uses this so every piece's segment is cropped to
     * its own wrap x slab, and adjacent segments land on the SAME integer
     * du/dv lattice (so their pixels co-register). */
    int    have_window;
    double win_u0, win_u1, win_v0, win_v1;
} TifxyzOpts;

void TifxyzOpts_default(TifxyzOpts *o);

typedef struct {
    size_t W, H;              /* grid px */
    double u0, v0;            /* canvas origin (vox, scroll_raster mapping) */
    size_t valid;             /* covered px */
    size_t multi;             /* px with >1 cover */
    size_t conflicts;         /* px with a cover > conflict_dist from the
                                 first (parameterization defect) */
    size_t skip_uv, skip_3d, skip_own;   /* gated faces */
    double bbox_lo[3], bbox_hi[3];       /* world (x,y,z) over valid px */
    double area_vx2;          /* valid 2x2 quads * du*dv */
} TifxyzStats;

/* Write <seg_dir>/{x,y,z}.tif + mask.tif (+ winding.tif + provenance.tif) +
 * meta.json with meta uuid = `uuid` (VC3D convention: uuid == dir basename).
 * Creates seg_dir. Returns 0 on success; -1 on error (bad input, canvas
 * unreasonable or over the classic-TIFF 4 GB band limit -- raise du/dv --
 * zero covered pixels, or I/O failure). Nothing is written on -1 except
 * possibly partial band files from a failed late write. */
int TifxyzExport_run(Arena_T arena, const PieceSet *ps, const char *seg_dir,
                     const char *uuid, const TifxyzOpts *opts,
                     TifxyzStats *out);

/* In-process unit tests (synthetic piece sets under
 * output/_selftest_scroll_unroll/). Returns 0 if all pass. */
int TifxyzExport_selftest(void);

#endif
