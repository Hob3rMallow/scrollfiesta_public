/* tif_probe -- dump per-column-band value histogram of a single-page uint8
 * TIFF (obj_bake_raw diag maps). Authoritative check that GDI+/PowerShell
 * viewers are not lying about content.
 *   tif_probe <file.tif> [nbands]        class-map mode (values 0..5)
 *   tif_probe <file.tif> [nbands] --cont continuous mode (0..255 buckets) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../src/common/arena.h"
#include "../../src/common/tiff_io.h"

int main(int argc, char **argv)
{
    Arena_T arena = NULL;
    uint8_t *vol = NULL;
    int D = 0, H = 0, W = 0;
    int nbands = 8;
    int b = 0, band_w = 0;
    long total[6] = { 0, 0, 0, 0, 0, 0 };
    long overflow = 0;
    int cont = 0, i = 0;

    if (argc < 2) { fprintf(stderr, "usage: %s <file.tif> [nbands] [--cont]\n", argv[0]); return 1; }
    for (i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--cont") == 0) cont = 1;
        else nbands = atoi(argv[i]);
    }
    if (nbands < 1) nbands = 1;
    arena = Arena_new();
    if (TiffIO_load(arena, argv[1], &vol, &D, &H, &W) != 0) {
        fprintf(stderr, "cannot load %s\n", argv[1]); return 1;
    }
    fprintf(stderr, "loaded %s: D=%d H=%d W=%d\n", argv[1], D, H, W);
    if (cont) {
        /* continuous 0..255: report nonzero pixels split into value ranges.
         * For verr (*32): <32=<1vox 32..96=1..3vox 96..160=3..5vox >160=>5vox.
         * For stretch/squash (128=iso): report the same byte buckets. */
        long nz = 0, b_lo = 0, b_mid = 0, b_hi = 0, b_vhi = 0;
        long below128 = 0, at128 = 0, above128 = 0, far = 0;
        int y = 0, x = 0;
        for (y = 0; y < H; y++) for (x = 0; x < W; x++) {
            int v = vol[(size_t)y * W + x];
            if (v == 0) continue;
            nz++;
            if (v < 32) b_lo++; else if (v < 96) b_mid++;
            else if (v < 160) b_hi++; else b_vhi++;
            if (v < 124) below128++; else if (v <= 132) at128++; else above128++;
            if (v < 86 || v > 170) far++;   /* >~1 octave from iso */
        }
        fprintf(stderr, "nonzero=%ld of %ld (%.1f%%)\n", nz, (long)W*H,
                100.0*(double)nz/(double)(W*H));
        fprintf(stderr, "  verr-view:  <1vox=%ld (%.1f%%) 1-3vox=%ld (%.1f%%) "
                "3-5vox=%ld (%.1f%%) >5vox=%ld (%.1f%%)\n",
                b_lo, 100.0*b_lo/(double)(nz?nz:1), b_mid, 100.0*b_mid/(double)(nz?nz:1),
                b_hi, 100.0*b_hi/(double)(nz?nz:1), b_vhi, 100.0*b_vhi/(double)(nz?nz:1));
        fprintf(stderr, "  sigma-view: <iso=%ld (%.1f%%) ~iso=%ld (%.1f%%) "
                ">iso=%ld (%.1f%%)  |>1 octave|=%ld (%.1f%%)\n",
                below128, 100.0*below128/(double)(nz?nz:1), at128, 100.0*at128/(double)(nz?nz:1),
                above128, 100.0*above128/(double)(nz?nz:1), far, 100.0*far/(double)(nz?nz:1));
        Arena_dispose(&arena);
        return 0;
    }
    band_w = (W + nbands - 1) / nbands;
    fprintf(stderr, "per %d-col band: [class 0..5 counts] (values>=6 flagged)\n",
            band_w);
    for (b = 0; b < nbands; b++) {
        long cnt[6] = { 0, 0, 0, 0, 0, 0 }, over = 0;
        int x0 = b * band_w, x1 = (b + 1) * band_w, y = 0, x = 0;
        long nz = 0;
        double vmax = 0.0, vsum = 0.0;
        if (x1 > W) x1 = W;
        for (y = 0; y < H; y++) {
            for (x = x0; x < x1; x++) {
                int v = vol[(size_t)y * W + x];
                if (v < 6) cnt[v]++; else { over++; overflow++; }
                if (v > 0) { nz++; vsum += v; if (v > vmax) vmax = v; }
                if (v < 6) total[v]++;
            }
        }
        fprintf(stderr, "  u[%6d,%6d]: 0=%ld 1=%ld 2=%ld 3=%ld 4=%ld 5=%ld  "
                ">=6:%ld  nz=%ld maxv=%.0f meannz=%.1f\n",
                (int)(x0), (int)(x1), cnt[0], cnt[1], cnt[2], cnt[3], cnt[4],
                cnt[5], over, nz, vmax, nz ? vsum / (double)nz : 0.0);
    }
    fprintf(stderr, "TOTAL: 0=%ld 1=%ld 2=%ld 3=%ld 4=%ld 5=%ld  >=6:%ld\n",
            total[0], total[1], total[2], total[3], total[4], total[5], overflow);
    Arena_dispose(&arena);
    return 0;
}
