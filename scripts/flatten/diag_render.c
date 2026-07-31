/* diag_render -- colorize an obj_bake_raw diagnostic (or the grayscale
 * texture) map, stack it into N vertical strips, and write a PNG DIRECTLY
 * (vendored stb_image_write). Done in C because GDI+ silently truncates the
 * very wide (24k+ px) source TIFFs on read -- libtiff loads them whole.
 *
 *   diag_render <in.tif> <out.png> <mode> [--strips N] [--gap G] [--rows y0 y1]
 *   mode: class | signed | verr | gray
 * class : 0 bg->navy 1 ok->gray 2 skip-uv->red 3 skip-3d->orange
 *         4 multi->steelblue 5 dark->magenta
 * signed: 0->navy; <128 blue (stretch); >128 red (compress); 128 gray
 * verr  : 0->navy; <1vox green; 1-3vox yellow; >3vox red  (value = |dv|*32)
 * gray  : 0->navy (holes), else grayscale passthrough (the baked texture)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../src/common/arena.h"
#include "../../src/common/tiff_io.h"
#include "../../src/common/ves_png.h"

static void colormap(const char *mode, int v, unsigned char *rgb)
{
    if (strcmp(mode, "gray") == 0) {
        if (v == 0) { rgb[0]=12; rgb[1]=14; rgb[2]=40; }   /* hole -> navy */
        else { rgb[0]=(unsigned char)v; rgb[1]=(unsigned char)v; rgb[2]=(unsigned char)v; }
    } else if (strcmp(mode, "snap") == 0) {
        /* surface_snap classes: 0 bg 1 good 2 fixable 3 crack 4 anchorless 5 mixed */
        switch (v) {
        case 0: rgb[0]=12;  rgb[1]=14;  rgb[2]=40;  break;  /* bg navy */
        case 1: rgb[0]=110; rgb[1]=110; rgb[2]=110; break;  /* good gray */
        case 2: rgb[0]=33;  rgb[1]=204; rgb[2]=64;  break;  /* fixable green */
        case 3: rgb[0]=224; rgb[1]=31;  rgb[2]=31;  break;  /* crack red */
        case 4: rgb[0]=245; rgb[1]=140; rgb[2]=25;  break;  /* anchorless orange */
        case 5: rgb[0]=230; rgb[1]=217; rgb[2]=38;  break;  /* mixed yellow */
        default: rgb[0]=0; rgb[1]=0; rgb[2]=0; break;
        }
    } else if (strcmp(mode, "class") == 0) {
        switch (v) {
        case 0: rgb[0]=12;  rgb[1]=14;  rgb[2]=40;  break;  /* bg navy */
        case 1: rgb[0]=150; rgb[1]=150; rgb[2]=150; break;  /* ok gray */
        case 2: rgb[0]=220; rgb[1]=30;  rgb[2]=30;  break;  /* skip-uv red */
        case 3: rgb[0]=240; rgb[1]=140; rgb[2]=20;  break;  /* skip-3d orange */
        case 4: rgb[0]=60;  rgb[1]=90;  rgb[2]=150; break;  /* multi steelblue */
        case 5: rgb[0]=230; rgb[1]=20;  rgb[2]=210; break;  /* dark magenta */
        case 6: rgb[0]=30;  rgb[1]=200; rgb[2]=190; break;  /* skip-own teal */
        case 7: rgb[0]=110; rgb[1]=230; rgb[2]=90;  break;  /* synth green */
        default: rgb[0]=0; rgb[1]=0; rgb[2]=0; break;
        }
    } else if (strcmp(mode, "signed") == 0) {
        if (v == 0) { rgb[0]=12; rgb[1]=14; rgb[2]=40; }
        else if (v < 128) {
            int t = 128 - v;
            rgb[0] = 40;
            rgb[1] = (unsigned char)(60 + t > 255 ? 255 : 60 + t);
            rgb[2] = (unsigned char)(80 + t*2 > 255 ? 255 : 80 + t*2);
        } else {
            int t = v - 128;
            rgb[0] = (unsigned char)(80 + t*2 > 255 ? 255 : 80 + t*2);
            rgb[1] = (unsigned char)(120 - t < 0 ? 0 : 120 - t);
            rgb[2] = 40;
        }
    } else {   /* verr */
        if (v == 0) { rgb[0]=12; rgb[1]=14; rgb[2]=40; }
        else if (v < 32) { rgb[0]=40;  rgb[1]=200; rgb[2]=60; }
        else if (v < 96) { rgb[0]=230; rgb[1]=210; rgb[2]=30; }
        else             { rgb[0]=230; rgb[1]=30;  rgb[2]=30; }
    }
}

/* snaptex compositing: tint a surface_snap class `c` over the grayscale RAW
 * texture value `g`, preserving the fiber luminance (multiply hue*g/255).
 *   c 1 good -> plain gray; 2 fixable green; 3 crack red; 4 anchorless orange;
 *   5 mixed yellow; texture hole (g==0) -> solid class hue if dark else navy.
 * `tmin` floors the tinted classes so the dominant channel is >= tmin even over
 * dark texture (fixable areas ARE dark, so without this the green is near-black);
 * above the floor, g still modulates so the fibers stay visible. */
static void tint_snap(int c, int g, int tmin, unsigned char *rgb)
{
    int tinted = (c >= 2 && c <= 5);
    if (g == 0) {                       /* texture hole: show class solid, else navy */
        if (tinted) colormap("snap", c, rgb);
        else { rgb[0]=12; rgb[1]=14; rgb[2]=40; }
        return;
    }
    if (tinted) {
        unsigned char hue[3];
        int hmax = 0, gmin = 0, ge = g;
        colormap("snap", c, hue);       /* solid class hue */
        hmax = hue[0];
        if (hue[1] > hmax) hmax = hue[1];
        if (hue[2] > hmax) hmax = hue[2];
        /* smallest g whose dominant channel reaches tmin (ceil div) */
        gmin = (hmax > 0) ? (tmin * 255 + hmax - 1) / hmax : 255;
        if (ge < gmin) ge = gmin;
        if (ge > 255) ge = 255;
        rgb[0] = (unsigned char)(hue[0] * ge / 255);
        rgb[1] = (unsigned char)(hue[1] * ge / 255);
        rgb[2] = (unsigned char)(hue[2] * ge / 255);
    } else {                            /* good / uncovered under texture -> plain gray */
        if (g < 40) g = 40;             /* floor so faint fibers still show */
        rgb[0] = rgb[1] = rgb[2] = (unsigned char)g;
    }
}

int main(int argc, char **argv)
{
    Arena_T arena = NULL;
    uint8_t *vol = NULL, *ovl = NULL;
    int D = 0, H = 0, W = 0, OD = 0, OH = 0, OW = 0;
    const char *mode = NULL, *overlay = NULL;
    int nstrips = 8, gap = 6, tint_min = 64;
    int sw = 0, oh = 0, s = 0, y = 0, x = 0;
    int y0 = 0, y1 = 0, bandH = 0, a = 0;
    unsigned char *out = NULL;

    if (argc < 4) {
        fprintf(stderr, "usage: %s <in.tif> <out.png> <class|snap|snaptex|signed|verr|gray>"
                " [--strips N] [--gap G] [--rows y0 y1] [--overlay class.tif] [--tint-min N]\n"
                "  --overlay: in gray mode, tint MULTI-cover (overlap) pixels "
                "red over the grayscale texture\n"
                "  snaptex: <in.tif> is the RAW texture, --overlay is a surface_snap\n"
                "           class map; tints each class hue over the grayscale fibers\n", argv[0]);
        return 1;
    }
    mode = argv[3];
    for (a = 4; a < argc; a++) {
        if (strcmp(argv[a], "--rows") == 0 && a + 2 < argc) {
            y0 = atoi(argv[a+1]); y1 = atoi(argv[a+2]); a += 2;
        } else if (strcmp(argv[a], "--strips") == 0 && a + 1 < argc) {
            nstrips = atoi(argv[++a]);
        } else if (strcmp(argv[a], "--gap") == 0 && a + 1 < argc) {
            gap = atoi(argv[++a]);
        } else if (strcmp(argv[a], "--overlay") == 0 && a + 1 < argc) {
            overlay = argv[++a];
        } else if (strcmp(argv[a], "--tint-min") == 0 && a + 1 < argc) {
            tint_min = atoi(argv[++a]);   /* snaptex: floor tinted classes to stay visible */
        }
    }
    if (nstrips < 1) nstrips = 1;

    arena = Arena_new();
    if (TiffIO_load(arena, argv[1], &vol, &D, &H, &W) != 0) {
        fprintf(stderr, "cannot load %s\n", argv[1]); return 1;
    }
    if (overlay != NULL) {
        if (TiffIO_load(arena, overlay, &ovl, &OD, &OH, &OW) != 0
            || OW != W || OH != H) {
            fprintf(stderr, "overlay %s missing or size mismatch (%dx%d vs %dx%d)\n",
                    overlay ? overlay : "", OW, OH, W, H);
            return 1;
        }
    }
    if (strcmp(mode, "snaptex") == 0 && ovl == NULL) {
        fprintf(stderr, "snaptex requires --overlay <snapclass.tif>\n");
        return 1;
    }
    /* optional row band: render only y in [y0,y1) as a single full-width row
     * (nstrips forced to 1) so a thin v-band is legible across all of u. */
    if (y1 > y0) {
        if (y0 < 0) y0 = 0;
        if (y1 > H) y1 = H;
        bandH = y1 - y0;   /* stacked into the usual strips, just this v-range */
    } else {
        y0 = 0; bandH = H;
    }
    sw = (W + nstrips - 1) / nstrips;               /* strip width */
    oh = nstrips * (bandH + gap) - gap;              /* stacked height */
    fprintf(stderr, "loaded %dx%d -> %d strips of %dx%d, out %dx%d\n",
            W, H, nstrips, sw, H, sw, oh);
    out = (unsigned char *)calloc((size_t)sw * oh * 3, 1);
    if (out == NULL) { fprintf(stderr, "oom\n"); return 1; }

    for (s = 0; s < nstrips; s++) {
        int x0 = s * sw, oy0 = s * (bandH + gap);
        for (y = 0; y < bandH; y++) {
            for (x = 0; x < sw; x++) {
                int sx = x0 + x;
                size_t oidx = ((size_t)(oy0 + y) * sw + x) * 3;
                size_t sidx;
                if (sx >= W) continue;               /* ragged last strip */
                sidx = (size_t)(y0 + y) * W + sx;
                /* snaptex: base = RAW texture (vol), tint by snap class (ovl) */
                if (strcmp(mode, "snaptex") == 0) {
                    tint_snap((int)ovl[sidx], (int)vol[sidx], tint_min, &out[oidx]);
                    continue;
                }
                colormap(mode, vol[sidx], &out[oidx]);
                /* overlap overlay: MULTI-cover pixels (class 4) keep the
                 * papyrus in the red channel and are damped in green/blue, so
                 * overlaps read as red-tinted texture -- you see BOTH the
                 * papyrus and exactly where two wraps share a (u,v). */
                if (ovl != NULL && ovl[sidx] == 4) {
                    int g = vol[sidx];
                    if (g < 40) g = 40;              /* floor so faint overlaps show */
                    out[oidx+0] = (unsigned char)g;              /* R = texture */
                    out[oidx+1] = (unsigned char)(g * 3 / 10);   /* G damped */
                    out[oidx+2] = (unsigned char)(g * 3 / 10);   /* B damped */
                }
            }
        }
    }

    if (VesPng_write_rgb(argv[2], out, sw, oh) != 0) {
        fprintf(stderr, "cannot write %s\n", argv[2]); free(out); return 1;
    }
    free(out);
    Arena_dispose(&arena);
    fprintf(stderr, "wrote %s (%dx%d RGB PNG)\n", argv[2], sw, oh);
    return 0;
}
