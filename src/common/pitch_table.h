#ifndef PITCH_TABLE_INCLUDED
#define PITCH_TABLE_INCLUDED

/* pitch_table -- measured radial wrap-pitch table for winding gates.
 *
 * Every winding gate/detector evaluates dw = dr/pitch - dtheta/(2*pi).  The
 * historical scalar pitch (9.5 vox) is wrong by >2x against measurement
 * (p50 ~ 20 vox mid-sheet on PHerc0139, 28-33 at r 60-120, ~14 outside
 * r 410), which makes real cross-wrap hops under-score and tangential
 * same-wrap chords over-score.  This module supplies pitch(r) as a
 * piecewise-linear table measured from the prediction
 * (scripts/radial_pitch_probe.py --pitch-table-out) and the exact winding
 * increment between two radii as the analytic integral of dr/pitch(r).
 *
 * CSV format ('#' comment lines and an optional "r,pitch..." header):
 *     r,pitch[,count]
 * with r ascending.  Outside the table the end pitch values extend as
 * constants.
 *
 * The FIT side keeps its scalar wrap spacing (shared world-u scale);
 * this table is for GATES and DETECTORS only. */

typedef struct PitchTable *PitchTable;

/* Load a table; returns 0 and sets *out, nonzero on IO/parse failure
 * (*out = NULL).  Free with PitchTable_free. */
int PitchTable_load(const char *path, PitchTable *out);
void PitchTable_free(PitchTable table);

/* Piecewise-linear pitch at radius r (clamped to the table ends).
 * table == NULL is invalid here; gates use PitchTable_dturns which
 * falls back to the scalar. */
double PitchTable_pitch_at(PitchTable table, double r);

/* Winding turns spanned between two radii: integral of dr/pitch(r) from
 * r_a to r_b (signed).  table == NULL yields exactly
 * (r_b - r_a) / scalar_fallback -- bit-compatible with the historical
 * expression, so an unarmed gate is unchanged. */
double PitchTable_dturns(PitchTable table, double r_a, double r_b,
                         double scalar_fallback);

/* Process-global table from the VES_WRAP_PITCH_TABLE environment variable.
 * Loaded once on first call (NULL when unset or unloadable; a load failure
 * is reported to stderr once).  Intended for gate arming paths. */
PitchTable PitchTable_from_env(void);

/* Number of knots (0 for NULL). */
int PitchTable_knots(PitchTable table);

/* Deterministic synthetic tests; returns number of failures. */
int PitchTable_selftest(void);

#endif
