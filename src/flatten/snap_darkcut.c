/* ============================================================================
 * snap_darkcut.c -- see snap_darkcut.h. Binary DARK/GOOD graph cut via GCO.
 * ==========================================================================*/

#include "snap_darkcut.h"
#include "../common/gco_wrap.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#define SD_DARK 0
#define SD_GOOD 1
/* Keep ALL GCO costs O(1e4): from experience GCO v3.0 misbehaves with large
 * per-site costs (an INF-scale anchor over ~1M verts cascades the whole mesh to
 * one label). SCALE is chosen so the max data cost SCALE*(hi-dk) stays under the
 * anchor; only the data:smoothness RATIO drives the cut, so scaling down is free. */
#define SD_SCALE   64.0         /* raw-unit -> integer cost scale */
#define SD_ANCHOR  10000        /* forced-GOOD data cost (strongest cost; bounded) */
#define SD_WMAX    1073741824.0 /* 2^30 neighbour-weight clamp (GCO int range) */

static int32_t clampi(double x)
{
    if (x < 0.0) x = 0.0;
    if (x > (double)SD_ANCHOR) x = (double)SD_ANCHOR;   /* data never exceeds the anchor */
    return (int32_t)x;
}

int SnapDark_segment(const double *cv, const uint8_t *has, size_t nv,
                     CSR_T adj, const SnapDarkOpts *opts, uint8_t *out_dark)
{
    double dk = opts->dk;
    double lambda = opts->lambda;
    double sigma = (opts->sigma > 1e-6) ? opts->sigma : 1.0;
    int32_t *data = NULL, smooth[4] = { 0, 1, 1, 0 };  /* Potts */
    int *labels = NULL;
    const int32_t *off = NULL, *tgt = NULL;
    GCO_Handle gc = NULL;
    size_t i = 0;

    if (nv == 0 || nv > (size_t)INT_MAX) return -1;

    /* data[site*2 + label]: pivot at dk so lambda=0 == (cv < dk).
     *   cv <  dk : DARK free, GOOD costs SCALE*(dk-cv)  -> DARK
     *   cv >= dk : GOOD free, DARK costs SCALE*(cv-dk)+1 -> GOOD (ties -> GOOD) */
    data = (int32_t *)malloc(nv * 2 * sizeof(int32_t));
    if (data == NULL) return -1;
    for (i = 0; i < nv; i++) {
        if (!has[i] || (opts->band_hi > 0.0 && cv[i] > opts->band_hi)) {
            data[i*2 + SD_DARK] = SD_ANCHOR; /* no sample, or clearly bright -> force GOOD */
            data[i*2 + SD_GOOD] = 0;
            continue;
        }
        {
            double d = cv[i] - dk;
            if (d < 0.0) {
                data[i*2 + SD_DARK] = 0;
                data[i*2 + SD_GOOD] = clampi(SD_SCALE * (-d));
            } else {
                data[i*2 + SD_DARK] = clampi(SD_SCALE * d) + 1;
                data[i*2 + SD_GOOD] = 0;
            }
        }
    }

    gc = GCO_create((int)nv, 2);
    if (gc == NULL) { free(data); return -1; }
    GCO_set_data_cost(gc, data);
    GCO_set_smooth_cost(gc, smooth);

    /* contrast-sensitive n-links on the 1-ring (u<v, each edge once). Skip edges
     * touching a no-sample vertex: it is forced GOOD and must not drag neighbours. */
    if (lambda > 0.0) {
        off = CSR_offset(adj);
        tgt = CSR_target(adj);
        for (i = 0; i < nv; i++) {
            int32_t e = 0;
            for (e = off[i]; e < off[i+1]; e++) {
                int32_t v = tgt[e];
                double dcv = 0.0, w = 0.0, ws = 0.0;
                if ((size_t)v <= i) continue;
                if (!has[i] || !has[v]) continue;
                dcv = cv[i] - cv[(size_t)v];
                w = lambda * exp(-(dcv * dcv) / (2.0 * sigma * sigma));
                ws = SD_SCALE * w;
                if (ws < 1.0) ws = 1.0;
                if (ws > SD_WMAX) ws = SD_WMAX;
                GCO_set_neighbor(gc, (int)i, (int)v, (int)ws);
            }
        }
    }

    GCO_expansion(gc, -1);

    labels = (int *)malloc(nv * sizeof(int));
    if (labels == NULL) { GCO_destroy(gc); free(data); return -1; }
    GCO_get_labels(gc, labels, (int)nv);
    for (i = 0; i < nv; i++) out_dark[i] = (labels[i] == SD_DARK) ? (uint8_t)1 : (uint8_t)0;

    free(labels);
    GCO_destroy(gc);
    free(data);
    return 0;
}
