#ifndef ATLAS_RIBBON_TEXTURE_INCLUDED
#define ATLAS_RIBBON_TEXTURE_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "../unroll/piece_set.h"
#include "atlas_ribbon_fit.h"
#include "atlas_solution.h"

typedef struct {
    int search_columns;       /* +/- source observation columns (8) */
    size_t min_samples;       /* paired RAW samples required (24) */
    double min_correlation;   /* Pearson NCC acceptance floor (.20) */
    double min_margin;        /* best-vs-non-neighbor peak (.02) */
    double normal_range;      /* normal-max RAW sampling reach (2) */
    int normal_steps;         /* normal-max taps (5) */
} AtlasRibbonTextureOptions;

typedef struct {
    int32_t row0;
    int32_t row1;
    int32_t chart;
    int32_t group;
    int32_t winding;
    int32_t lag_columns;
    uint32_t samples;
    int accepted;
    double correlation;
    double zero_correlation;
    double margin;
    double desired_shift_delta;
    double current_shift_delta;
    double residual;
} AtlasRibbonTexturePair;

typedef struct {
    AtlasRibbonTexturePair *pair;
    size_t npairs;
    size_t accepted_pairs;
    size_t sampled_layer_points;
    size_t missing_layer_points;
    int cubes_loaded;
    int cubes_missing;
    double accepted_correlation_mean;
    double accepted_margin_mean;
    double current_delta_rms;
    double desired_delta_rms;
    double residual_rms;
    double residual_p95;
    double residual_max;
} AtlasRibbonTextureResult;

void AtlasRibbonTextureOptions_default(AtlasRibbonTextureOptions *opts);

/* Measure RAW texture phase between consecutive rows of each retained chart.
 * No U field is changed. desired_shift_delta is the row1-row0 chart-shift
 * difference that aligns the measured source-space signals. */
int AtlasRibbonTexture_measure(
    Arena_T arena,
    const PieceSet *ps,
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *observations,
    const AtlasRibbonFitResult *fit,
    const char *raw_dir,
    const AtlasRibbonTextureOptions *opts,
    AtlasRibbonTextureResult *out);

/* Recompute current row-delta and residual statistics after the U(v) field
 * changes. RAW-derived lags, correlations, margins, and acceptance stay fixed. */
int AtlasRibbonTexture_refresh(
    Arena_T scratch,
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *observations,
    const AtlasRibbonFitResult *fit,
    AtlasRibbonTextureResult *phase);

/* Deterministic sign/phase fixture for the one-dimensional matcher. */
int AtlasRibbonTexture_selftest(void);

#endif
