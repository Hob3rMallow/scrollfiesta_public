#ifndef ASM_FIBRE_FRAME_INCLUDED
#define ASM_FIBRE_FRAME_INCLUDED

#include <stddef.h>
#include "asm_types.h"
#include "../common/zarr_u8.h"

/* asm_fibre_frame.h -- the papyrus fibre cross of assembled charts, measured on the CT along each chart, and the turn
 * that lays a piece's fibres along the page axes (2026-10-21, the user: "align all blocks going in ... such that the
 * principle fiber axis directions align with the X and Y axes (and such that the main lengthwise fiber direction gets
 * the X axis) ... this should be the default behaviour going forwards").
 *
 * Measured: PHerc0358's fibres cross at about +-45 deg to the winding ribbon's rows while z runs down its columns, so
 * a close-up's fibre alignment had to turn every patch by up to 45 deg and folded the page; other scrolls' fibres sit
 * within ~5 deg of the rows.
 *
 * A chart is rasterized on a 1 vox grid of its current UV frame; each covered node samples the CT (nearest voxel,
 * averaged over -1, 0, +1 vox along the interpolated normal), the image is smoothed twice with a 3x3 binomial, and
 * every node with a complete neighbourhood adds its oriented second-derivative energy at spacings 2 and 4 vox:
 * the second derivative across direction theta is D = P + Q cos 2theta + R sin 2theta (P = (Suu + Svv)/2,
 * Q = (Suu - Svv)/2, R = Suv) and D^2 = w + Re(m2 e^{-i2theta}) + Re(m4 e^{-i4theta}) with
 *     w += P^2 + (Q^2 + R^2)/2,  m2 += 2P (Q + iR),  m4 += (Q + iR)^2 / 2.
 * Two superposed fibre families add to D linearly (cos^2 and sin^2 of the angle to their normals), so the cross
 * reinforces the 4-fold term: arg(m4) / 4 is the cross and the 2-fold term says which family is stronger.  (The first
 * selftest measured gradient moments instead: the vector sum of two families' gradients cancels their 4-fold term,
 * coherence 0.02 on a clean crosshatch.)  The moments add over charts that share a frame and map exactly through a
 * reflection and a rotation (AsmFibre_transform), so the region pack sums them per piece in its presentation frame. */

typedef struct AsmFibreMoments { double w, m2[2], m4[2]; size_t px; } AsmFibreMoments;

typedef struct AsmFibreStats { size_t charts, px, skipped; size_t chunk_misses, chunk_missing; double seconds; } AsmFibreStats;

/* Every registered chart of run into out[run->n_charts] (zero for the others), from the uncompressed uint8 zarr v2
 * at raw (<raw>/0/.zarray).  0 on success, -1 when the volume cannot be opened or memory runs out. */
int AsmFibre_measure(const AsmRun *run, const char *raw, AsmFibreMoments *out, AsmFibreStats *st);

/* The same CT cross measurement for one caller-owned physical tangent patch.
 * Sampling coordinates describe the measurement plane, not a flattening result.
 * The caller keeps a volume/cache open across patches. No geometry is changed.
 * face_samples[nf], when present, counts valid derivative centers supported by
 * each triangle. Distinct world points projecting onto the same sample make
 * that sample unknown, never a first-face-wins measurement. */
void AsmFibre_measure_patch(ZarrU8 *volume, Arena_T arena, const AsmChart *chart,
                            AsmFibreMoments *out, size_t *face_samples);

/* The moments in the frame p -> R(cosine, sine) M p, where M reflects u -> -u when mirror is set (anp_frame's order). */
void AsmFibre_transform(AsmFibreMoments *m, int mirror, double cosine, double sine);

/* One piece's decision, from its moments in its presentation frame:
 *   near   the fibre family nearest +u, at cross_deg in (-45, 45]; far = near + 90 deg
 *   decided_by 1 (lengthwise): |cross_deg| <= 45 - ASM_FIBRE_TIE_DEG, the near family takes +u (the winding direction
 *              the pack already laid along +u, where a rolled papyrus's lengthwise fibres run);
 *   decided_by 2 (main): a near-diagonal cross, the family with more fibre energy takes +u;
 *   decided_by 0: no measurable cross (coherence |m4| / w under ASM_FIBRE_MIN_COHERENCE or too few nodes), no turn.
 * turn_deg: the rotation (counter-clockwise from +u toward +v) that lays the chosen family along +u. */
typedef struct AsmFibreTurn {
    double coherence, cross_deg, energy_near, energy_far, turn_deg;
    int decided_by;
} AsmFibreTurn;
void AsmFibre_turn(const AsmFibreMoments *m, AsmFibreTurn *out);

int AsmFibre_selftest(void);

#endif
