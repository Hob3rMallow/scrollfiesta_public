#ifndef WINDING_MRF_INCLUDED
#define WINDING_MRF_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"

/* Integer winding-label MRF.
 *
 * Each site receives an integer label k.  A field observation contributes a
 * robust unary around center, while an edge contributes
 *
 *     weight * |(k_b - k_a) - target|
 *
 * This lets imperfect generalized-winding measurements vote without being
 * rounded into hard truth.  Multiple, even contradictory, observations for
 * the same site pair are retained as separate factors. */
typedef struct {
    double center;          /* generalized-winding estimate */
    double sigma;           /* observation scale; <= 0 disables unary */
    double weight;          /* non-negative unary weight */
    int32_t initial_label;  /* INT32_MIN chooses round(center) */
    int fixed;              /* non-zero fixes the site to initial_label */
} WindingMRFSite;

typedef struct {
    int32_t a, b;
    int32_t target;         /* preferred label[b] - label[a] */
    double weight;          /* non-negative observation weight */
} WindingMRFEdge;

/* Additional sparse unary terms.  A penalty is paid only when `site` takes
 * `label`.  This is deliberately separate from the robust field unary: an
 * outer conflict-resolution pass can rule out a previously occupied winding
 * without pretending that it measured a new Gaussian center. */
typedef struct {
    int32_t site;
    int32_t label;
    double weight;
} WindingMRFUnaryPenalty;

typedef struct {
    int32_t label_min, label_max; /* inclusive common label window */
    int max_iterations;           /* alpha-beta swap passes */
    double cost_scale;            /* natural energy -> graph-cut integer */
    double unary_huber_delta;     /* robust z-score transition */
    double confidence_temperature;
    double abstain_threshold;
} WindingMRFOptions;

typedef struct {
    size_t sites, input_edges, active_edges, unique_pairs, labels;
    size_t changed_labels, abstained_sites;
    int32_t solution_min, solution_max;
    int64_t quantized_energy_before, quantized_energy_after;
    double energy_before, energy_after;
    double mean_confidence, min_confidence;
} WindingMRFStats;

void WindingMRF_default_options(WindingMRFOptions *options);

/* Outputs are arena-owned arrays of nsites entries.  Confidence is a local
 * conditional posterior sharpness in [0,1], with all neighboring MAP labels
 * held fixed; it is deliberately not presented as an exact global marginal. */
int WindingMRF_solve(
    Arena_T arena,
    const WindingMRFSite *sites, size_t nsites,
    const WindingMRFEdge *edges, size_t nedges,
    const WindingMRFOptions *options,
    int32_t **out_labels,
    float **out_confidence,
    WindingMRFStats *stats);

/* Penalty-aware form used by iterative conflict exclusion.  Duplicate
 * (site,label) entries add.  The ordinary entry point above is exactly this
 * solve with no sparse penalties. */
int WindingMRF_solve_with_penalties(
    Arena_T arena,
    const WindingMRFSite *sites, size_t nsites,
    const WindingMRFEdge *edges, size_t nedges,
    const WindingMRFUnaryPenalty *penalties, size_t npenalties,
    const WindingMRFOptions *options,
    int32_t **out_labels,
    float **out_confidence,
    WindingMRFStats *stats);

/* Synthetic noisy-loop, duplicate-factor, abstention, and the real five-site
 * PHerc0139 Artifact C seam-alias regression. */
int WindingMRF_selftest(void);

#endif
