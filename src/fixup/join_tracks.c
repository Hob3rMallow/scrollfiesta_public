/* join_tracks.c — cross-plane ConnectionTrack consistency.
 * See join_tracks.h. Ported from 744d0ba^:.../gap_fill_2d_v1.c (the
 * EndpointTrack/ConnectionTrack machinery), arena-allocated, member lists
 * dropped (only the aggregates are needed for scoring). */

#include "join_tracks.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "../common/pipeline_constants.h"

typedef struct {
    float   mean_y, mean_x;          /* running means */
    float   mean_tan_dy, mean_tan_dx;
    int32_t last_z;
    int32_t n_members;
    int32_t best_conn;               /* strongest ConnectionTrack, or -1 */
    float   best_conn_conf;
} jt_track;

typedef struct {
    int32_t track_a, track_b;        /* canonical a < b */
    int32_t n_confirmed;
    int32_t n_possible;
    float   confidence;
} jt_conn;

struct JoinTracks {
    Arena_T   arena;
    int32_t   n_planes;
    int32_t **plane_track_ids;       /* [n_planes][n_eps(z)], -1 = none */
    jt_track *tracks;
    int32_t   n_tracks;
    jt_conn  *conns;
    int32_t   n_conns;
    int32_t  *conn_map;              /* open-addressing hash, -1 = empty */
    int32_t   conn_map_cap;          /* power of two */
};

static uint32_t jt_hash(int32_t a, int32_t b)
{
    return (uint32_t)a * 2654435761u + (uint32_t)b * 40503u;
}

static int32_t jt_conn_lookup(const struct JoinTracks *jt,
                              int32_t a, int32_t b)
{
    if (jt->conn_map_cap == 0) return -1;
    uint32_t mask = (uint32_t)(jt->conn_map_cap - 1);
    uint32_t h = jt_hash(a, b) & mask;
    for (int32_t probe = 0; probe < jt->conn_map_cap; probe++) {
        uint32_t slot = (h + (uint32_t)probe) & mask;
        int32_t idx = jt->conn_map[slot];
        if (idx < 0) return -1;
        if (jt->conns[idx].track_a == a && jt->conns[idx].track_b == b)
            return idx;
    }
    return -1;
}

static void jt_conn_insert(struct JoinTracks *jt, int32_t a, int32_t b,
                           int32_t conn_idx)
{
    uint32_t mask = (uint32_t)(jt->conn_map_cap - 1);
    uint32_t h = jt_hash(a, b) & mask;
    for (int32_t probe = 0; probe < jt->conn_map_cap; probe++) {
        uint32_t slot = (h + (uint32_t)probe) & mask;
        if (jt->conn_map[slot] < 0) {
            jt->conn_map[slot] = conn_idx;
            return;
        }
    }
    assert(0 && "jt conn hash full");
}

/* Grow-by-copy in the arena (counts are small; waste is bounded). */
static jt_track *jt_grow_tracks(Arena_T arena, jt_track *old,
                                int32_t n, int32_t *cap)
{
    int32_t new_cap = (*cap == 0) ? 1024 : *cap * 2;
    jt_track *g = ARENA_ALLOC(arena, (size_t)new_cap * sizeof(*g));
    if (n > 0) memcpy(g, old, (size_t)n * sizeof(*g));
    *cap = new_cap;
    return g;
}

static jt_conn *jt_grow_conns(Arena_T arena, jt_conn *old,
                              int32_t n, int32_t *cap)
{
    int32_t new_cap = (*cap == 0) ? 2048 : *cap * 2;
    jt_conn *g = ARENA_ALLOC(arena, (size_t)new_cap * sizeof(*g));
    if (n > 0) memcpy(g, old, (size_t)n * sizeof(*g));
    *cap = new_cap;
    return g;
}

JoinTracks_T JoinTracks_build(Arena_T arena, const FixupPlaneView *planes,
                              int32_t n_planes)
{
    assert(planes || n_planes == 0);
    struct JoinTracks *jt = NULL;
    ARENA_NEW(arena, jt);
    memset(jt, 0, sizeof(*jt));
    jt->arena = arena;
    jt->n_planes = n_planes;
    jt->plane_track_ids =
        ARENA_ALLOC(arena, ((size_t)n_planes > 0 ? (size_t)n_planes : 1) *
                    sizeof(int32_t *));

    int32_t track_cap = 0;
    jt->tracks = jt_grow_tracks(arena, NULL, 0, &track_cap);

    /* --- Pass A: endpoint-track association (serial z-scan). ---
     * The per-plane id tables AND the growable track array are RESULTS that
     * later passes and the queries read — allocate the tables up front and
     * use NO scratch mark here (the tracks array grows by arena copy inside
     * the loop, so any mark spanning the loop would free the live copy: the
     * slice_trace result-inside-scratch-mark bug, third sighting — at 10x
     * scale the recycled bytes were Pass C floats read back as track ids). */
    for (int32_t z = 0; z < n_planes; z++) {
        const FixupPlaneView *pv = &planes[z];
        jt->plane_track_ids[z] =
            ARENA_ALLOC(arena, ((size_t)pv->n_eps > 0 ? (size_t)pv->n_eps
                                                      : 1) *
                        sizeof(int32_t));
    }
    int32_t live_cap = 4096;
    int32_t *live = ARENA_ALLOC(arena, (size_t)live_cap * sizeof(int32_t));

    for (int32_t z = 0; z < n_planes; z++) {
        const FixupPlaneView *pv = &planes[z];
        int32_t *ids = jt->plane_track_ids[z];

        /* live tracks: recently seen */
        int32_t n_live = 0;
        for (int32_t t = 0; t < jt->n_tracks; t++) {
            if (z - jt->tracks[t].last_z > FIXUP_TRACK_GAP_TOL) continue;
            if (n_live >= live_cap) {
                int32_t new_cap = live_cap * 2;
                int32_t *g = ARENA_ALLOC(arena,
                                         (size_t)new_cap * sizeof(int32_t));
                memcpy(g, live, (size_t)n_live * sizeof(int32_t));
                live = g;
                live_cap = new_cap;
            }
            live[n_live++] = t;
        }

        for (int32_t i = 0; i < pv->n_eps; i++) {
            ids[i] = -1;
            if (pv->eps[i].excluded) continue;

            float ey = (float)pv->eps[i].y;
            float ex = (float)pv->eps[i].x;
            float ety = pv->eps[i].tan_dy;
            float etx = pv->eps[i].tan_dx;

            int32_t best_track = -1;
            float best_d2 = FIXUP_TRACK_SPATIAL * FIXUP_TRACK_SPATIAL;
            for (int32_t li = 0; li < n_live; li++) {
                jt_track *trk = &jt->tracks[live[li]];
                if (trk->last_z == z) continue;   /* one endpoint per plane */
                float dy = ey - trk->mean_y;
                float dx = ex - trk->mean_x;
                float d2 = dy * dy + dx * dx;
                if (d2 >= best_d2) continue;
                if (ety * trk->mean_tan_dy + etx * trk->mean_tan_dx < 0.0f)
                    continue;
                best_d2 = d2;
                best_track = live[li];
            }

            if (best_track >= 0) {
                jt_track *trk = &jt->tracks[best_track];
                trk->n_members++;
                float alpha = 1.0f / (float)trk->n_members;
                trk->mean_y = trk->mean_y * (1.0f - alpha) + ey * alpha;
                trk->mean_x = trk->mean_x * (1.0f - alpha) + ex * alpha;
                trk->mean_tan_dy =
                    trk->mean_tan_dy * (1.0f - alpha) + ety * alpha;
                trk->mean_tan_dx =
                    trk->mean_tan_dx * (1.0f - alpha) + etx * alpha;
                trk->last_z = z;
                ids[i] = best_track;
            } else {
                if (jt->n_tracks >= track_cap)
                    jt->tracks = jt_grow_tracks(arena, jt->tracks,
                                                jt->n_tracks, &track_cap);
                jt_track *trk = &jt->tracks[jt->n_tracks];
                trk->mean_y = ey;
                trk->mean_x = ex;
                trk->mean_tan_dy = ety;
                trk->mean_tan_dx = etx;
                trk->last_z = z;
                trk->n_members = 1;
                trk->best_conn = -1;
                trk->best_conn_conf = 0.0f;
                ids[i] = jt->n_tracks;
                jt->n_tracks++;
            }
        }
    }
    /* (live[] stays in the arena — a few hundred KB at worst; a restore
     * here would free the grown tracks array and every id table) */

    /* --- Pass B: ConnectionTracks from matches. --- */
    int32_t conn_cap = 0;
    jt->conns = jt_grow_conns(arena, NULL, 0, &conn_cap);
    for (int32_t z = 0; z < n_planes; z++) {
        const FixupPlaneView *pv = &planes[z];
        const int32_t *ids = jt->plane_track_ids[z];
        for (int32_t i = 0; i < pv->n_eps; i++) {
            int32_t j = pv->match ? pv->match[i] : -1;
            if (j <= i) continue;
            int32_t ta = ids[i], tb = ids[j];
            if (ta < 0 || tb < 0) continue;
            if (ta > tb) { int32_t tmp = ta; ta = tb; tb = tmp; }
            if (jt->conn_map_cap > 0 && jt_conn_lookup(jt, ta, tb) >= 0)
                continue;
            if (jt->n_conns >= conn_cap ||
                jt->conn_map_cap < 4 * (jt->n_conns + 1)) {
                /* grow storage and rebuild the hash */
                if (jt->n_conns >= conn_cap)
                    jt->conns = jt_grow_conns(arena, jt->conns, jt->n_conns,
                                              &conn_cap);
                int32_t want = 4 * conn_cap;
                int32_t cap2 = 1;
                while (cap2 < want) cap2 *= 2;
                jt->conn_map_cap = cap2;
                jt->conn_map = ARENA_ALLOC(arena,
                                           (size_t)cap2 * sizeof(int32_t));
                memset(jt->conn_map, 0xFF, (size_t)cap2 * sizeof(int32_t));
                for (int32_t c = 0; c < jt->n_conns; c++)
                    jt_conn_insert(jt, jt->conns[c].track_a,
                                   jt->conns[c].track_b, c);
            }
            int32_t ci = jt->n_conns++;
            jt->conns[ci].track_a = ta;
            jt->conns[ci].track_b = tb;
            jt->conns[ci].n_confirmed = 0;
            jt->conns[ci].n_possible = 0;
            jt->conns[ci].confidence = 0.0f;
            jt_conn_insert(jt, ta, tb, ci);
        }
    }

    /* Count confirmations. */
    for (int32_t z = 0; z < n_planes; z++) {
        const FixupPlaneView *pv = &planes[z];
        const int32_t *ids = jt->plane_track_ids[z];
        for (int32_t i = 0; i < pv->n_eps; i++) {
            int32_t j = pv->match ? pv->match[i] : -1;
            if (j <= i) continue;
            int32_t ta = ids[i], tb = ids[j];
            if (ta < 0 || tb < 0) continue;
            if (ta > tb) { int32_t tmp = ta; ta = tb; tb = tmp; }
            int32_t ci = jt_conn_lookup(jt, ta, tb);
            if (ci >= 0) jt->conns[ci].n_confirmed++;
        }
    }

    /* --- Pass C: n_possible via per-plane presence stamps. --- */
    if (jt->n_conns > 0) {
        Arena_Mark stamp_mark = Arena_save(arena);
        int32_t *stamp = ARENA_ALLOC(arena,
                                     ((size_t)jt->n_tracks > 0 ?
                                      (size_t)jt->n_tracks : 1) *
                                     sizeof(int32_t));
        for (int32_t t = 0; t < jt->n_tracks; t++) stamp[t] = -1;
        for (int32_t z = 0; z < n_planes; z++) {
            const FixupPlaneView *pv = &planes[z];
            const int32_t *ids = jt->plane_track_ids[z];
            for (int32_t i = 0; i < pv->n_eps; i++)
                if (ids[i] >= 0) stamp[ids[i]] = z;
            for (int32_t c = 0; c < jt->n_conns; c++)
                if (stamp[jt->conns[c].track_a] == z &&
                    stamp[jt->conns[c].track_b] == z)
                    jt->conns[c].n_possible++;
        }
        Arena_restore(arena, stamp_mark);
    }

    /* Confidence + best-connection links. */
    for (int32_t c = 0; c < jt->n_conns; c++) {
        jt_conn *ct = &jt->conns[c];
        ct->confidence = (ct->n_possible > 0)
            ? (float)ct->n_confirmed / (float)ct->n_possible : 0.0f;
        jt_track *ta = &jt->tracks[ct->track_a];
        jt_track *tb = &jt->tracks[ct->track_b];
        if (ct->confidence > ta->best_conn_conf) {
            ta->best_conn_conf = ct->confidence;
            ta->best_conn = c;
        }
        if (ct->confidence > tb->best_conn_conf) {
            tb->best_conn_conf = ct->confidence;
            tb->best_conn = c;
        }
    }

    return jt;
}

static int32_t jt_pair_conn(JoinTracks_T jt, int32_t z, int32_t a, int32_t b)
{
    assert(jt);
    if (z < 0 || z >= jt->n_planes) return -1;
    const int32_t *ids = jt->plane_track_ids[z];
    int32_t ta = ids[a], tb = ids[b];
    assert(ta < jt->n_tracks && tb < jt->n_tracks);
    if (ta < 0 || tb < 0 || ta >= jt->n_tracks || tb >= jt->n_tracks)
        return -1;
    if (ta > tb) { int32_t tmp = ta; ta = tb; tb = tmp; }
    return jt_conn_lookup(jt, ta, tb);
}

float *JoinTracks_s_track_matrix(Arena_T arena, JoinTracks_T jt, int32_t z,
                                 const FixupPlaneView *plane)
{
    assert(jt && plane);
    int32_t n = plane->n_eps;
    float *m = ARENA_CALLOC(arena, ((size_t)n > 0 ? (size_t)n * (size_t)n : 1),
                            sizeof(float));
    if (z < 0 || z >= jt->n_planes) return m;
    const int32_t *ids = jt->plane_track_ids[z];

    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = i + 1; j < n; j++) {
            int32_t ta = ids[i], tb = ids[j];
            /* defensive: an out-of-range id means a corrupted table — treat
             * as untracked rather than dereference garbage */
            assert(ta < jt->n_tracks && tb < jt->n_tracks);
            if (ta >= jt->n_tracks) ta = -1;
            if (tb >= jt->n_tracks) tb = -1;
            float s = 0.0f;
            if (ta >= 0 && tb >= 0) {
                int32_t qa = ta < tb ? ta : tb;
                int32_t qb = ta < tb ? tb : ta;
                int32_t ci = jt_conn_lookup(jt, qa, qb);
                if (ci >= 0) {
                    const jt_conn *ct = &jt->conns[ci];
                    if (ct->confidence >= FIXUP_TRACK_CONF_STRONG &&
                        ct->n_confirmed >= FIXUP_MIN_SUPPORT)
                        s = FIXUP_TRACK_BOOST_STRONG;
                    else if (ct->confidence >= FIXUP_TRACK_CONF_MOD)
                        s = FIXUP_TRACK_BOOST_MOD;
                }
                if (s <= 0.0f) {
                    /* contradiction: either track's best connection is a
                     * strong link to a DIFFERENT track */
                    int contradicts = 0;
                    const jt_track *trk_a = &jt->tracks[ta];
                    const jt_track *trk_b = &jt->tracks[tb];
                    if (trk_a->best_conn >= 0 &&
                        trk_a->best_conn_conf >= FIXUP_TRACK_CONF_STRONG) {
                        const jt_conn *bc = &jt->conns[trk_a->best_conn];
                        int32_t other = (bc->track_a == ta) ? bc->track_b
                                                            : bc->track_a;
                        if (other != tb) contradicts = 1;
                    }
                    if (trk_b->best_conn >= 0 &&
                        trk_b->best_conn_conf >= FIXUP_TRACK_CONF_STRONG) {
                        const jt_conn *bc = &jt->conns[trk_b->best_conn];
                        int32_t other = (bc->track_a == tb) ? bc->track_b
                                                            : bc->track_a;
                        if (other != ta) contradicts = 1;
                    }
                    if (contradicts) s = FIXUP_TRACK_ANTI;
                }
            }
            m[(size_t)i * (size_t)n + (size_t)j] = s;
            m[(size_t)j * (size_t)n + (size_t)i] = s;
        }
    }
    return m;
}

int32_t JoinTracks_pair_support(JoinTracks_T jt, int32_t z,
                                int32_t a, int32_t b)
{
    int32_t ci = jt_pair_conn(jt, z, a, b);
    return (ci >= 0) ? jt->conns[ci].n_confirmed : 0;
}

int32_t JoinTracks_pair_conn_id(JoinTracks_T jt, int32_t z,
                                int32_t a, int32_t b)
{
    return jt_pair_conn(jt, z, a, b);
}

int32_t JoinTracks_n_endpoint_tracks(JoinTracks_T jt)
{
    assert(jt);
    return jt->n_tracks;
}

int32_t JoinTracks_n_connection_tracks(JoinTracks_T jt)
{
    assert(jt);
    return jt->n_conns;
}

/* ---------------------------------------------------------------- */
/* Selftest                                                         */
/* ---------------------------------------------------------------- */

int JoinTracks_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();

    /* 8 planes. Endpoints 0 and 1 sit at fixed positions with facing
     * tangents, matched on planes 0..5 (6 planes). Endpoint 2 is a drifter
     * matched with endpoint 3 only on plane 3. */
    enum { NP = 8 };
    SliceTraceEndpoint *eps = ARENA_CALLOC(arena, (size_t)NP * 4,
                                           sizeof(SliceTraceEndpoint));
    int32_t *match = ARENA_ALLOC(arena, (size_t)NP * 4 * sizeof(int32_t));
    FixupPlaneView planes[NP];

    for (int32_t z = 0; z < NP; z++) {
        SliceTraceEndpoint *pe = &eps[z * 4];
        int32_t *pm = &match[z * 4];
        for (int k = 0; k < 4; k++) pm[k] = -1;
        /* pair A: (30,40)->(30,50), tangents facing */
        pe[0].y = 30; pe[0].x = 40; pe[0].tan_dy = 0.0f; pe[0].tan_dx = 1.0f;
        pe[1].y = 30; pe[1].x = 50; pe[1].tan_dy = 0.0f; pe[1].tan_dx = -1.0f;
        /* pair B: far corner */
        pe[2].y = 80; pe[2].x = 20; pe[2].tan_dy = 0.0f; pe[2].tan_dx = 1.0f;
        pe[3].y = 80; pe[3].x = 28; pe[3].tan_dy = 0.0f; pe[3].tan_dx = -1.0f;
        if (z < 6) { pm[0] = 1; pm[1] = 0; }
        if (z == 3) { pm[2] = 3; pm[3] = 2; }
        planes[z].eps = pe;
        planes[z].match = pm;
        planes[z].n_eps = 4;
    }

    JoinTracks_T jt = JoinTracks_build(arena, planes, NP);

    /* Track association should give 4 stable tracks. */
    {
        int ok = (JoinTracks_n_endpoint_tracks(jt) == 4 &&
                  JoinTracks_n_connection_tracks(jt) == 2);
        fprintf(stderr,
                "[selftest] join_tracks build (tracks=%d conns=%d) -> %s\n",
                (int)JoinTracks_n_endpoint_tracks(jt),
                (int)JoinTracks_n_connection_tracks(jt), ok ? "ok" : "FAIL");
        if (!ok) fails++;
    }

    /* Persistent pair: support 6, s_track strong. Outlier pair: support 1. */
    {
        int32_t supA = JoinTracks_pair_support(jt, 2, 0, 1);
        int32_t supB = JoinTracks_pair_support(jt, 3, 2, 3);
        float *m = JoinTracks_s_track_matrix(arena, jt, 2, &planes[2]);
        float sA = m[0 * 4 + 1];
        int ok = (supA == 6 && supB == 1 &&
                  fabsf(sA - FIXUP_TRACK_BOOST_STRONG) < 1e-6f);
        fprintf(stderr,
                "[selftest] join_tracks support (A=%d B=%d sA=%.2f) -> %s\n",
                (int)supA, (int)supB, (double)sA, ok ? "ok" : "FAIL");
        if (!ok) fails++;
    }

    /* Contradiction: pairing endpoint 0 with endpoint 3 (whose track's best
     * connection is elsewhere... 0's best is with 1) -> ANTI. */
    {
        float *m = JoinTracks_s_track_matrix(arena, jt, 3, &planes[3]);
        float s = m[0 * 4 + 3];
        int ok = (fabsf(s - FIXUP_TRACK_ANTI) < 1e-6f);
        fprintf(stderr,
                "[selftest] join_tracks contradiction (s=%.2f) -> %s\n",
                (double)s, ok ? "ok" : "FAIL");
        if (!ok) fails++;
    }

    /* Empty registry queries behave. */
    {
        JoinTracks_T empty = JoinTracks_build(arena, NULL, 0);
        int ok = (JoinTracks_n_endpoint_tracks(empty) == 0 &&
                  JoinTracks_pair_support(empty, 0, 0, 1) == 0);
        fprintf(stderr, "[selftest] join_tracks empty -> %s\n",
                ok ? "ok" : "FAIL");
        if (!ok) fails++;
    }

    Arena_dispose(&arena);
    fprintf(stderr, "=== join_tracks selftest %s (%d failure%s) ===\n",
            fails ? "FAILED" : "passed", fails, fails == 1 ? "" : "s");
    return fails ? 3 : 0;
}
