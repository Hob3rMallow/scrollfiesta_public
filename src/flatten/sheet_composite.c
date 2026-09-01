#include "sheet_composite.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/ves_png.h"

/* provenance palette (RGB): per layer, [direct, void-filled]; must match the
 * retired python compositor so historical provenance PNGs stay comparable */
static const uint8_t SC_PROV[SHEET_COMPOSITE_MAX_LAYERS][2][3] = {
    { { 235, 235, 235 }, { 150, 140, 140 } },
    { { 235, 140,  60 }, { 150,  85,  40 } },
    { {  80, 180, 220 }, {  50, 110, 140 } },
};
static const uint8_t SC_BG[3]    = { 12, 12, 18 };
static const uint8_t SC_FILLC[3] = { 210, 200, 120 };

/* Luminance-preserving chroma blend.  The auxiliary palette is normalized by
 * its integer BT.709 luminance before 3/5 of its chroma is mixed into the gray
 * texture.  Fibers therefore remain readable at approximately the same
 * brightness instead of becoming a flat provenance swatch. */
enum { SC_TINT_NUM = 3, SC_TINT_DEN = 5 };

static void sc_tint_pixel(uint8_t gray, size_t layer, uint8_t rgb[3])
{
    if (layer == 0 || layer >= SHEET_COMPOSITE_MAX_LAYERS) {
        rgb[0] = rgb[1] = rgb[2] = gray;
        return;
    }
    const uint8_t *base = SC_PROV[layer][0];
    int luma = (54 * (int)base[0] + 183 * (int)base[1] +
                19 * (int)base[2] + 128) / 256;
    int denominator = SC_TINT_DEN * luma;
    for (int c = 0; c < 3; c++) {
        int scale = (SC_TINT_DEN - SC_TINT_NUM) * luma +
                    SC_TINT_NUM * (int)base[c];
        int value = ((int)gray * scale + denominator / 2) / denominator;
        rgb[c] = (uint8_t)(value > 255 ? 255 : value);
    }
}

static int sc_canvas_size(const SheetCompositeLayer *layers, size_t n_layers,
                          size_t *out_w, size_t *out_h)
{
    size_t W = 0, H = 0;
    if (layers == NULL || n_layers == 0 ||
        n_layers > SHEET_COMPOSITE_MAX_LAYERS)
        return -1;
    for (size_t l = 0; l < n_layers; l++) {
        const SheetCompositeLayer *L = &layers[l];
        if (L->tex == NULL || L->cov == NULL || L->w == 0 || L->h == 0 ||
            L->off_u < 0 || L->off_v < 0)
            return -1;
        if ((size_t)L->off_u > SIZE_MAX - L->w ||
            (size_t)L->off_v > SIZE_MAX - L->h)
            return -1;
        if ((size_t)L->off_u + L->w > W) W = (size_t)L->off_u + L->w;
        if ((size_t)L->off_v + L->h > H) H = (size_t)L->off_v + L->h;
    }
    if (W == 0 || H == 0 || W > SIZE_MAX / H || W * H > SIZE_MAX / 3)
        return -1;
    *out_w = W;
    *out_h = H;
    return 0;
}

static void *sc_calloc(size_t n, size_t sz)
{
    void *p = calloc(n, sz);
    if (p == NULL) fprintf(stderr, "sheet_composite: out of memory\n");
    return p;
}

int SheetComposite_run(const SheetCompositeLayer *layers, size_t n_layers,
                       const char *out_tex_png, const char *out_prov_png,
                       uint8_t **out_tex_data, SheetCompositeStats *out)
{
    size_t W = 0, H = 0;
    if (out_tex_data != NULL) *out_tex_data = NULL;
    if (out == NULL || sc_canvas_size(layers, n_layers, &W, &H) != 0)
        return -1;
    memset(out, 0, sizeof *out);

    uint8_t *tex = (uint8_t *)sc_calloc(W * H, 1);
    uint8_t *prov = (uint8_t *)sc_calloc(W * H, 3);
    uint8_t *painted = (uint8_t *)sc_calloc(W * H, 1);
    uint8_t *owner = (uint8_t *)sc_calloc(W * H, 1);
    if (tex == NULL || prov == NULL || painted == NULL || owner == NULL) {
        free(tex); free(prov); free(painted); free(owner);
        return -1;
    }
    for (size_t p = 0; p < W * H; p++) {
        prov[p * 3 + 0] = SC_BG[0];
        prov[p * 3 + 1] = SC_BG[1];
        prov[p * 3 + 2] = SC_BG[2];
        owner[p] = UINT8_MAX;
    }

    /* Two passes: DIRECT face coverage from any layer beats a void-filled
     * crack pixel from an earlier layer; layer order breaks remaining ties. */
    for (int want_direct = 1; want_direct >= 0; want_direct--) {
        for (size_t l = 0; l < n_layers; l++) {
            const SheetCompositeLayer *L = &layers[l];
            size_t taken = 0;
            for (size_t y = 0; y < L->h; y++) {
                size_t drow = ((size_t)L->off_v + y) * W + (size_t)L->off_u;
                size_t srow = y * L->w;
                for (size_t x = 0; x < L->w; x++) {
                    uint8_t c = L->cov[srow + x];
                    int covered = want_direct ? (c == 255) : (c > 0);
                    if (!covered || painted[drow + x]) continue;
                    tex[drow + x] = L->tex[srow + x];
                    memcpy(&prov[(drow + x) * 3],
                           SC_PROV[l][want_direct ? 0 : 1], 3);
                    owner[drow + x] = (uint8_t)l;
                    painted[drow + x] = 1;
                    taken++;
                }
            }
            if (want_direct) {
                out->direct_px[l] = taken;
            } else {
                out->void_fill_px[l] = taken;
                /* occluded = covered pixels an earlier layer already owns */
                size_t occl = 0;
                for (size_t y = 0; y < L->h; y++) {
                    size_t drow = ((size_t)L->off_v + y) * W
                                + (size_t)L->off_u;
                    size_t srow = y * L->w;
                    for (size_t x = 0; x < L->w; x++)
                        if (L->cov[srow + x] > 0 && painted[drow + x]) occl++;
                }
                occl -= out->direct_px[l] + out->void_fill_px[l];
                out->occluded_px[l] = occl;
            }
        }
    }

    /* Composite-level bounded crack fill (the baker's own rule, applied where
     * BOTH sides of a lane seam are visible -- per-sheet fill cannot see
     * across sheets): an unpainted pixel with >= 5/8 painted neighbours, or
     * bracketed by painted pixels within 2 above AND below (or left AND
     * right), takes their mean.  Two passes close 1-4 px seams; real holes
     * stay honest.  Tagged separately in the provenance. */
    uint8_t *take = (uint8_t *)sc_calloc(W * H, 1);
    uint8_t *tval = (uint8_t *)sc_calloc(W * H, 1);
    if (take == NULL || tval == NULL) {
        free(tex); free(prov); free(painted); free(owner);
        free(take); free(tval);
        return -1;
    }
    for (int pass = 0; pass < 2; pass++) {
        size_t nfill = 0;
        long lw = (long)W, lh = (long)H;
        for (long y = 0; y < lh; y++) {
            for (long x = 0; x < lw; x++) {
                size_t p = (size_t)y * W + (size_t)x;
                if (painted[p]) { take[p] = 0; continue; }
                unsigned int acc = 0;
                int nsupp = 0;
                for (int dy = -1; dy <= 1; dy++)
                    for (int dx = -1; dx <= 1; dx++) {
                        if (dy == 0 && dx == 0) continue;
                        long yy = y + dy, xx = x + dx;
                        if (yy < 0 || yy >= lh || xx < 0 || xx >= lw)
                            continue;
                        size_t q = (size_t)yy * W + (size_t)xx;
                        if (painted[q]) { nsupp++; acc += tex[q]; }
                    }
                if (nsupp == 0) { take[p] = 0; continue; }
                int bracket = 0;
                {
                    int up = 0, dn = 0, lt = 0, rt = 0;
                    for (int d = 1; d <= 2 && !bracket; d++) {
                        if (y - d >= 0 &&
                            painted[(size_t)(y - d) * W + (size_t)x]) up = 1;
                        if (y + d < lh &&
                            painted[(size_t)(y + d) * W + (size_t)x]) dn = 1;
                        if (x - d >= 0 &&
                            painted[(size_t)y * W + (size_t)(x - d)]) lt = 1;
                        if (x + d < lw &&
                            painted[(size_t)y * W + (size_t)(x + d)]) rt = 1;
                        if ((up && dn) || (lt && rt)) bracket = 1;
                    }
                }
                if (nsupp >= 5 || bracket) {
                    take[p] = 1;
                    tval[p] = (uint8_t)(acc / (unsigned int)nsupp);
                    nfill++;
                } else {
                    take[p] = 0;
                }
            }
        }
        if (nfill == 0) break;
        for (size_t p = 0; p < W * H; p++) {
            if (!take[p]) continue;
            tex[p] = tval[p];
            memcpy(&prov[p * 3], SC_FILLC, 3);
            painted[p] = 1;
        }
        out->composite_crack_fill_px += nfill;
    }
    free(take);
    free(tval);

    size_t np = 0;
    for (size_t p = 0; p < W * H; p++) np += painted[p];
    out->width = W;
    out->height = H;
    out->painted_px = np;
    out->empty_px = W * H - np;

    int rc = 0;
    uint8_t *view = NULL;
    if (out_tex_png != NULL || out_tex_data != NULL) {
        view = (uint8_t *)sc_calloc(W * H, 3);
        if (view == NULL) {
            rc = -1;
        } else {
            for (size_t p = 0; p < W * H; p++)
                sc_tint_pixel(tex[p], owner[p], &view[p * 3]);
            if (out_tex_png != NULL &&
                VesPng_write_rgb(out_tex_png, view, (int)W, (int)H) != 0)
                rc = -1;
        }
    }
    if (out_prov_png != NULL &&
        VesPng_write_rgb(out_prov_png, prov, (int)W, (int)H) != 0)
        rc = -1;
    if (out_tex_data != NULL && rc == 0) *out_tex_data = view;
    else free(view);
    free(tex);
    free(prov);
    free(painted);
    free(owner);
    return rc;
}

int SheetComposite_write_layer_views(
    const SheetCompositeLayer *layers, size_t n_layers,
    const char *const out_layer_png[SHEET_COMPOSITE_MAX_LAYERS])
{
    size_t W = 0, H = 0;
    if (out_layer_png == NULL ||
        sc_canvas_size(layers, n_layers, &W, &H) != 0)
        return -1;
    for (size_t l = 0; l < n_layers; l++) {
        const SheetCompositeLayer *L = &layers[l];
        if (out_layer_png[l] == NULL) return -1;
        uint8_t *view = (uint8_t *)sc_calloc(W * H, 3);
        if (view == NULL) return -1;
        for (size_t y = 0; y < L->h; y++) {
            size_t drow = ((size_t)L->off_v + y) * W + (size_t)L->off_u;
            size_t srow = y * L->w;
            for (size_t x = 0; x < L->w; x++) {
                if (L->cov[srow + x] == 0) continue;
                sc_tint_pixel(L->tex[srow + x], l,
                              &view[(drow + x) * 3]);
            }
        }
        int rc = VesPng_write_rgb(out_layer_png[l], view, (int)W, (int)H);
        free(view);
        if (rc != 0) return -1;
    }
    return 0;
}

int SheetComposite_selftest(void)
{
    /* Layer 0: 8x4 with a 2-px vertical unpainted crack at x=3..4 and a
     * void-fill pixel; layer 1: 4x4 at off_u=2 covering the crack columns
     * with direct pixels only on x=2 (global 4).  Expect: direct-first
     * beats layer-0 void where layer 1 is direct; crack closes by the
     * bounded fill. */
    enum { W0 = 8, H0 = 4, W1 = 4, H1 = 4 };
    uint8_t tex0[W0 * H0], cov0[W0 * H0];
    uint8_t tex1[W1 * H1], cov1[W1 * H1];
    int fails = 0;
    for (int y = 0; y < H0; y++)
        for (int x = 0; x < W0; x++) {
            int p = y * W0 + x;
            int crack = (x == 3 || x == 4);
            tex0[p] = 100;
            cov0[p] = crack ? 0 : 255;
        }
    cov0[1 * W0 + 5] = 128;   /* layer-0 void-fill pixel */
    tex0[1 * W0 + 5] = 90;
    for (int y = 0; y < H1; y++)
        for (int x = 0; x < W1; x++) {
            int p = y * W1 + x;
            tex1[p] = 200;
            cov1[p] = (x == 2) ? 255 : 0;   /* global x = 4 */
        }
    /* give layer 1 a void pixel at local x=3,y=1 (global x=5,y=1): layer 0
     * void competes with layer 1 void there; layer order wins for layer 0 */
    cov1[1 * W1 + 3] = 128;

    SheetCompositeLayer L[2];
    memset(L, 0, sizeof L);
    L[0].tex = tex0; L[0].cov = cov0; L[0].w = W0; L[0].h = H0;
    L[1].tex = tex1; L[1].cov = cov1; L[1].w = W1; L[1].h = H1;
    L[1].off_u = 2; L[1].off_v = 0;

    SheetCompositeStats S;
    {
        uint8_t neutral[3], orange[3], blue[3];
        sc_tint_pixel(100, 0, neutral);
        sc_tint_pixel(100, 1, orange);
        sc_tint_pixel(100, 2, blue);
        if (neutral[0] != 100 || neutral[1] != 100 || neutral[2] != 100 ||
            !(orange[0] > orange[1] && orange[1] > orange[2]) ||
            !(blue[2] > blue[1] && blue[1] > blue[0])) {
            fprintf(stderr, "  FAIL: provenance texture tint palette\n");
            fails++;
        }
    }
    int rc = SheetComposite_run(L, 2, NULL, NULL, NULL, &S);
    if (rc != 0) { fprintf(stderr, "  FAIL: composite rc\n"); fails++; }
    if (S.width != 8 || S.height != 4) {
        fprintf(stderr, "  FAIL: composite canvas %zux%zu\n",
                S.width, S.height);
        fails++;
    }
    /* layer 0 direct = 8x4 minus 2x4 crack minus 1 void = 23 */
    if (S.direct_px[0] != 23) {
        fprintf(stderr, "  FAIL: layer0 direct %zu != 23\n", S.direct_px[0]);
        fails++;
    }
    /* layer 1 direct = its x==2 column (global 4, in the crack) = 4 */
    if (S.direct_px[1] != 4) {
        fprintf(stderr, "  FAIL: layer1 direct %zu != 4\n", S.direct_px[1]);
        fails++;
    }
    /* layer 0's void pixel survives (layer order beats layer 1 void) */
    if (S.void_fill_px[0] != 1 || S.void_fill_px[1] != 0) {
        fprintf(stderr, "  FAIL: void px %zu/%zu != 1/0\n",
                S.void_fill_px[0], S.void_fill_px[1]);
        fails++;
    }
    /* the remaining crack column (global x=3) is 1 px wide with painted
     * pixels left and right: the bounded fill closes all 4 */
    if (S.composite_crack_fill_px != 4) {
        fprintf(stderr, "  FAIL: crack fill %zu != 4\n",
                S.composite_crack_fill_px);
        fails++;
    }
    if (S.painted_px != 32) {
        fprintf(stderr, "  FAIL: painted %zu != 32\n", S.painted_px);
        fails++;
    }
    fprintf(stderr, "[sheet_composite selftest] %s (%d failure(s))\n",
            fails == 0 ? "PASS" : "FAIL", fails);
    return fails == 0 ? 0 : -1;
}
