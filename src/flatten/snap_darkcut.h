#ifndef SNAP_DARKCUT_INCLUDED
#define SNAP_DARKCUT_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/csr.h"

/* ============================================================================
 * snap_darkcut.h -- binary Boykov-Jolly graph cut (GCO alpha-expansion, exact
 * for 2 labels) segmenting mesh vertices into DARK vs GOOD, replacing a hard
 * per-vertex intensity threshold. The data term pivots at `dk` so that
 * lambda==0 reduces EXACTLY to `cv < dk`; the contrast-sensitive Potts
 * smoothness on the mesh 1-ring grows dark cores into adjacent DIM area up to a
 * sharp intensity edge (but stops at bright) -- capturing the true extent of
 * dark regions, not just the darkest cores. Calls GCO through its extern-C
 * wrapper; this stays C.
 * ==========================================================================*/

typedef struct SnapDarkOpts {
    double dk;       /* dark pivot (raw units) = data-term crossover == old threshold */
    double lambda;   /* smoothness weight; 0 => pure threshold (parity/AB) */
    double sigma;    /* contrast scale (raw units) in exp(-(dcv)^2 / 2 sigma^2) */
    double band_hi;  /* bright anchor: cv > band_hi => forced GOOD (data[DARK]=INF),
                      * so the cut cannot bleed into clearly-bright papyrus; <=0 off */
} SnapDarkOpts;

/* Fill out_dark[nv] with 1=DARK, 0=GOOD. cv[nv] = per-vertex intensity, has[nv]
 * = sampled flag (0 => forced GOOD). adj = mesh 1-ring (CSR_from_faces). Returns
 * 0 on success, -1 on failure (caller should fall back to the hard threshold). */
int SnapDark_segment(const double *cv, const uint8_t *has, size_t nv,
                     CSR_T adj, const SnapDarkOpts *opts, uint8_t *out_dark);

#endif
