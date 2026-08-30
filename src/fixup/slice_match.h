#ifndef SLICE_MATCH_INCLUDED
#define SLICE_MATCH_INCLUDED

#include <stdint.h>
#include "../common/arena.h"
#include "slice_trace.h"

/* slice_match — per-plane endpoint matching for the 2D prediction gap fixup.
 *
 * Builds candidate endpoint pairs, applies the safety gates, scores each pair
 * with ONE symmetric score, and accepts pairs greedily in score order while
 * both endpoints are free and the join's Bezier crosses no already-accepted
 * join. For a single symmetric preference this greedy IS the unique stable
 * matching — chosen over the 2026-05 prototype's two-sided Gale-Shapley,
 * whose bipartite top/bottom partition forced wrong pairings (see the
 * recovered MATCHING_ANALYSIS.md post-mortem).
 *
 * Score = W_DIST * gauss(d) + W_DIR * facing + W_ADJ * s_adj
 *       + W_TRACK * s_track,
 * where s_adj is LOCAL adjacent-plane connectivity evidence: a flood fill
 * restricted to the join's corridor in planes z +/- 1 (weight 1) and z +/- 2
 * (weight 0.5) must walk from a's fragment to b's fragment. Global per-plane
 * CC labels are useless for this — a healthy spiral cross-section is ONE
 * connected curve — so the flood is corridor-local by construction.
 *
 * Hard gates, cheap to expensive:
 *   excluded endpoints; same skeleton CC or same mask CC (would close a
 *   loop); reach cap; tangent facing + anti-parallel opposition; radial gate
 *   |r_a - r_b| <= radial_dr_max about the umbilicus (armed when umbilicus
 *   given — a cross-wrap join moves ~one pitch radially); far-tier evidence
 *   (d > reach_safe REQUIRES s_adj >= FIXUP_ADJ_MIN_EVIDENCE); score floor;
 *   Bezier arc-ratio; anti-merger third-CC margin; cross-sheet corridor in
 *   adjacent planes; no crossing of accepted joins.
 *
 * Deterministic: candidates sorted by (score desc, a asc, b asc). */

typedef struct SliceMatchParams {
    float reach_safe;     /* d <= this: standard gates only */
    float reach_max;      /* hard cap */
    float min_score;
    int   radial_armed;   /* 1 = umbilicus known, radial gate active */
    float umb_py, umb_px; /* umbilicus in PLANE coords (world - volume origin) */
    float radial_dr_max;
    int   provisional;    /* 1 = round 1 with a track round coming: skip the
                           * far-tier evidence gate (round-1 joins are never
                           * painted; they only seed ConnectionTracks). In the
                           * FINAL round the far tier passes on corridor
                           * evidence OR moderate track support — a PERSISTENT
                           * gap is broken in z±1 too, so corridor evidence
                           * only ever validates flicker gaps; persistence is
                           * exactly what the track term measures. */
} SliceMatchParams;

enum {
    SLICEMATCH_REJ_SAME_CC  = 1,  /* same skeleton or mask component */
    SLICEMATCH_REJ_TANGENT  = 2,  /* facing / opposition gate */
    SLICEMATCH_REJ_RADIAL   = 3,
    SLICEMATCH_REJ_EVIDENCE = 4,  /* far tier without adjacent-plane evidence */
    SLICEMATCH_REJ_SCORE    = 5,
    SLICEMATCH_REJ_ARC      = 6,
    SLICEMATCH_REJ_MERGER   = 7,
    SLICEMATCH_REJ_XSHEET   = 8,
    SLICEMATCH_REJ_CROSSING = 9,  /* would cross an accepted join */
    SLICEMATCH_REJ_OCCUPIED = 10, /* an endpoint already took a better join */
    SLICEMATCH_REJ_N_REASONS = 11
};

typedef struct SliceMatchReject {
    int32_t a, b;
    uint8_t reason;               /* SLICEMATCH_REJ_* */
} SliceMatchReject;

typedef struct SliceMatchJoin {
    int32_t a, b;                 /* endpoint indices, a < b */
    float   score, dist, s_adj, s_track;
    uint8_t far_tier;             /* 1 = dist > reach_safe */
} SliceMatchJoin;

typedef struct SliceMatchStats {
    int32_t n_pairs_in_reach;
    int32_t n_candidates;         /* survived every pair gate */
    int32_t n_accepted;
    int32_t rej_count[SLICEMATCH_REJ_N_REASONS];
} SliceMatchStats;

/* Match one plane z of the D*H*W binary volume `vol` (0/1, pristine — no
 * paint). mask_labels is plane z's CC labels; all_labels is the whole
 * volume's per-plane CC labels (for the cross-sheet gate). s_track is an
 * optional [n_eps*n_eps] matrix (round 2) or NULL (round 1).
 * Outputs (arena): *out_joins accepted joins; *out_rejects gated pairs (for
 * the viewers; capped, only pairs that were within reach). match_out is
 * caller-allocated [n_eps], -1 = unmatched. stats may be NULL.
 * Returns 0 on success; n_eps < 2 succeeds with zero joins. */
int SliceMatch_run(Arena_T arena,
                   const SliceTraceEndpoint *eps, int32_t n_eps,
                   const uint8_t *vol, const int32_t *all_labels,
                   int D, int H, int W, int z,
                   const int32_t *mask_labels,
                   const float *s_track,
                   const SliceMatchParams *params,
                   int32_t *match_out,
                   SliceMatchJoin **out_joins, int32_t *out_n_joins,
                   SliceMatchReject **out_rejects, int32_t *out_n_rejects,
                   SliceMatchStats *stats);

int SliceMatch_selftest(void);

#endif
