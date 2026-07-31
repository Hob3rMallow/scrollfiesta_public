#ifndef SCROLL_RASTER_INCLUDED
#define SCROLL_RASTER_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"
#include "../common/raw_sample.h"
#include "piece_set.h"

/* ============================================================================
 * scroll_raster.h -- whole-grid strip raster, column-tiled in u.
 *
 * Same per-pixel semantics as Rawtex_write_tif (rawtex_bake.c): every covered
 * pixel center is located in its UV triangle, the 3D point + normal
 * barycentrically interpolated, and the RAW volume freshly sampled there
 * (normal-max multi-tap); UV-stretch and long-3D-edge face gates. The absolute
 * 3D-edge gate is disabled for real source faces by default because their
 * valid edge scale follows the upstream remesher; synthetic fill keeps the
 * conservative 6-voxel membrane gate. Differences
 * that make it whole-grid capable:
 *   - the (u,v) frame is ABSOLUTE (can be negative/huge): the canvas origin
 *     is the piece-set u/v minimum, not 0;
 *   - the uint8 canvas is allocated at full size (W*H, ~0.5 GB at 4x21x21
 *     scale) but the double/uint32 ACCUMULATORS live per u-band
 *     (band_cols x H), so peak memory stays bounded at any strip width --
 *     this removes rawtex_bake's W*H > 2^26 cap;
 *   - one lazy CubeTable over cubes_RAW serves all bands.
 * Gate: band_cols=0 (single-shot) must be byte-identical to any banding.
 * ==========================================================================*/

typedef struct {
    double du, dv;           /* vox per pixel (default 1.0) */
    size_t band_cols;        /* accumulator band width in px; 0 = single band */
    double range;            /* normal-max sampling reach, vox (default 2.0) */
    int    nsteps;           /* samples along the normal (default 5) */
    double stretch_ratio;    /* face gate: uv edge > ratio*3d (default 4) */
    double stretch_floor;    /*   ... and > floor vox (default 25) */
    double max_edge3d;       /* REAL face gate: 3D edge > this => skip;
                                default 0 (disabled; source resolution varies) */
    double synth_max_edge3d; /* SYNTH face membrane gate (default 6) */
    double pct_lo, pct_hi;   /* contrast percentiles (default 1/99) */
    int    write_diag;       /* emit <prefix>_diagclass.tif (default 1) */
    int    write_xyzmap;     /* emit <prefix>_xyzmap.png (+ _strip): every
                                covered px colored by its world (x,y,z) as RGB
                                (R<-x, G<-y, B<-z), normalized over the REAL-face
                                vert bbox so all stages share one frame. Locally
                                smooth where the fit is smooth; kinks/jumps mark
                                where step6 fill blows up. Default 1. */
    int    strip_w;          /* stacked-strip target width (<=0 default) */
    CubeTable *ct;           /* shared RAW table (NULL = build own from
                                raw_dir). NEVER Arena_restore across it. */
    const uint8_t *face_keep;/* [nf] ownership mask (NULL = keep all); 0 =
                                face skipped, its sole-cover px -> class 6 */
    double lo_fixed, hi_fixed;/* contrast window override (hi>lo skips the
                                sampling pass -- keeps stage TIFs comparable) */
    size_t synth_f0;         /* faces >= synth_f0 are synthetic fill (step6
                                param_extend): they paint ONLY px no real face
                                claimed (not even as a gate diagnostic), as
                                class 7, so real coverage and its stats stay
                                byte-identical. (size_t)-1 = no synth faces. */
} RasterOpts;

void RasterOpts_default(RasterOpts *o);

typedef struct {
    size_t W, H;             /* canvas px */
    double u0, v0;           /* canvas origin (vox) */
    double lo, hi;           /* contrast window used */
    double fill;             /* covered / total px */
    size_t filled, multi;    /* covered px; px with >1 cover (overlaps) */
    size_t filled_synth;     /* class-7 px painted by synth faces (excluded
                                from filled/fill so stage rows stay
                                comparable) */
    size_t synth_nosample;   /* class-7 px where synth geometry exists but
                                the RAW sample failed (outside loaded cubes)
                                -- gray stays 0, tint still marks the ribbon */
    size_t skip_uv, skip_3d; /* gated REAL faces (synth gates don't count) */
    size_t skip_own;         /* faces dropped by opts->face_keep */
    size_t n_bands;
    int    raw_loaded, raw_missing;   /* CubeTable traffic */
    double xyz_lo[3], xyz_hi[3];      /* per-coord (z,y,x) normalization bbox
                                         actually used for the xyzmap, over
                                         verts of REAL faces (f < synth_f0) --
                                         stable across stages for comparison */
} RasterStats;

/* Rasterize the piece set against the RAW cube dir and write
 *   <prefix>_rawtex.tif (+ .png when narrow enough)
 *   <prefix>_rawtex_strip.tif/.png   (stacked review image)
 *   <prefix>_diagclass.tif  (0 bg|1 ok|2 skip-uv|3 skip-3d|4 multi|6 skip-own|
 *                            7 synth-fill; 5 stays reserved for rawtex_bake's
 *                            "dark")
 * Returns 0 on success, -1 on error. */
int ScrollRaster_run(Arena_T arena, const PieceSet *ps,
                     const char *raw_dir, long raw_chunk,
                     const char *prefix, const RasterOpts *opts,
                     RasterStats *out);

/* In-process unit tests (synthetic injected CubeTable; banding byte-identity;
 * face gates). Returns 0 if all pass. */
int ScrollRaster_selftest(void);

#endif
