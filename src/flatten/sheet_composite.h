#ifndef SHEET_COMPOSITE_INCLUDED
#define SHEET_COMPOSITE_INCLUDED

#include <stddef.h>
#include <stdint.h>

/* First-cover composite of quadribbon layer-sheet bakes into one big-sheet
 * view (C port of python/scripts/composite_layer_sheets.py, which leaves the
 * deliverable chain).  Sheet 0 is the primary surface; sheets 1..2 hold
 * charts that physically stack on it.  Painting is per-pixel first-cover in
 * sheet order with DIRECT face coverage from any layer beating an earlier
 * layer's void-filled crack pixel; a bounded composite-level crack fill then
 * closes 1-4 px lane seams that per-sheet fill cannot see across.  This is a
 * VIEW over the unchanged layer bakes -- no coordinate, face, or assignment
 * is modified -- and the provenance image records which layer (and whether
 * void-fill) supplied every pixel. */

enum { SHEET_COMPOSITE_MAX_LAYERS = 3 };

typedef struct SheetCompositeLayer {
    const uint8_t *tex;   /* [w*h] grayscale bake */
    const uint8_t *cov;   /* [w*h] coverage: 255 direct, 128 void-fill, 0 bg */
    size_t w, h;
    long off_u, off_v;    /* integer placement in the common frame, >= 0 */
} SheetCompositeLayer;

typedef struct SheetCompositeStats {
    size_t width, height;
    size_t painted_px;
    size_t composite_crack_fill_px;
    size_t empty_px;
    size_t direct_px[SHEET_COMPOSITE_MAX_LAYERS];
    size_t void_fill_px[SHEET_COMPOSITE_MAX_LAYERS];
    size_t occluded_px[SHEET_COMPOSITE_MAX_LAYERS];
} SheetCompositeStats;

/* Composites the layers and writes <out_tex_png> as a texture-bearing RGB
 * provenance view (layer 0 neutral, layer 1 orange, layer 2 blue) and
 * <out_prov_png> as the flat 24-bit provenance key; either path may be NULL.
 * If out_tex_data is non-NULL it receives the same tinted RGB canvas
 * (malloc'd, width*height*3 bytes; caller frees) so callers can build stacked
 * comparison views without re-reading PNGs.  Returns 0 on success. */
int SheetComposite_run(const SheetCompositeLayer *layers, size_t n_layers,
                       const char *out_tex_png, const char *out_prov_png,
                       uint8_t **out_tex_data, SheetCompositeStats *out);

/* Write one texture-bearing, full-common-canvas PNG per provenance layer.
 * Pixels outside that layer's own coverage are black.  Paths for existing
 * layers must be non-NULL; the array has SHEET_COMPOSITE_MAX_LAYERS slots. */
int SheetComposite_write_layer_views(
    const SheetCompositeLayer *layers, size_t n_layers,
    const char *const out_layer_png[SHEET_COMPOSITE_MAX_LAYERS]);

int SheetComposite_selftest(void);

#endif
