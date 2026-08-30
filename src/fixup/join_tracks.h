#ifndef JOIN_TRACKS_INCLUDED
#define JOIN_TRACKS_INCLUDED

#include <stdint.h>
#include "../common/arena.h"
#include "slice_trace.h"

/* join_tracks — cross-plane consistency for the 2D prediction gap fixup
 * (the TRACK_MATCHING.md design recovered from 744d0ba^, arena-ported).
 *
 * EndpointTracks associate "the same endpoint" across consecutive planes
 * (spatial distance < FIXUP_TRACK_SPATIAL, tangent agreement, staleness
 * bound FIXUP_TRACK_GAP_TOL). ConnectionTracks aggregate, per pair of
 * endpoint tracks, how often round-1 matching actually paired them
 * (n_confirmed) out of the planes where both were present (n_possible);
 * confidence = confirmed / possible.
 *
 * Round 2 re-runs the per-plane matcher with the s_track term:
 *   +FIXUP_TRACK_BOOST_STRONG  conf >= CONF_STRONG and confirmed >= MIN_SUPPORT
 *   +FIXUP_TRACK_BOOST_MOD     conf >= CONF_MOD
 *    FIXUP_TRACK_ANTI          either endpoint's BEST connection is to a
 *                              different track (contradicts the pattern)
 *    0                         no track evidence
 *
 * The final support filter keeps only joins whose ConnectionTrack was
 * confirmed in >= FIXUP_MIN_SUPPORT planes (mirrors skeleton_stack.py's
 * junction_tube_min_planes = 3): a join that exists in one plane only is
 * skeleton noise, a real sheet gap persists in z.
 *
 * Build is a serial z-scan (deterministic); all queries are read-only and
 * thread-safe afterwards. */

typedef struct FixupPlaneView {
    const SliceTraceEndpoint *eps;   /* [n_eps] */
    const int32_t            *match; /* [n_eps], partner index or -1 */
    int32_t                   n_eps;
} FixupPlaneView;

typedef struct JoinTracks *JoinTracks_T;

/* Build the registry from per-plane trace + match results. */
JoinTracks_T JoinTracks_build(Arena_T arena, const FixupPlaneView *planes,
                              int32_t n_planes);

/* s_track matrix [n_eps*n_eps] for plane z (arena-allocated, symmetric). */
float *JoinTracks_s_track_matrix(Arena_T arena, JoinTracks_T jt, int32_t z,
                                 const FixupPlaneView *plane);

/* Number of planes in which the (a,b) pair's ConnectionTrack was confirmed
 * (0 if either endpoint has no track or the pair was never matched). */
int32_t JoinTracks_pair_support(JoinTracks_T jt, int32_t z,
                                int32_t a, int32_t b);

/* ConnectionTrack id for the pair, or -1. Stable across queries; usable as
 * a join-track identifier in manifests. */
int32_t JoinTracks_pair_conn_id(JoinTracks_T jt, int32_t z,
                                int32_t a, int32_t b);

int32_t JoinTracks_n_endpoint_tracks(JoinTracks_T jt);
int32_t JoinTracks_n_connection_tracks(JoinTracks_T jt);

int JoinTracks_selftest(void);

#endif
